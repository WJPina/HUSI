library(reticulate)
library(dplyr)
library(magrittr)
library(Seurat)
library(mclust)
library(tibble)

setwd("/mnt/data1/wangj/AgingScore/GSE163530_COVID-19/GSE171668_scnRNA/")
set.seed(233)

mm_l2 = readRDS("/home/wangjing/wangj/AgingScore/Data/Bulk_TrainModel/mm_l2.rds")
sc <- import("scanpy")
###### preprocess data
### raw data
anndata = sc$read_h5ad("Covid_Lung_Endo_1010.h5ad")
table(anndata$obs$Cluster)

anndata_use = anndata
Counts <- anndata_use$raw$to_adata()$copy()$T$to_df()
Meta <- anndata_use$obs

Endo.m <- CreateSeuratObject(counts=Counts,
                             meta.data=Meta[,c('Cluster','SubCluster','Viral+','donor','disease','hUSI')],
                             min.cells = 0,min.features = 0)
Endo.m <- NormalizeData(Endo.m)
hist(colSums(Endo.m@assays$RNA@data))

table(Endo.m$SubCluster)

### load clinical data
ClinicMeta = read.csv('/home/wangjing/wangj/AgingScore/GSE163530_COVID-19/GSE162911_GeoMx/ClinicMeta.csv',row.names = 1)
ClinicMeta <- ClinicMeta[complete.cases(ClinicMeta$Days_to_death),]
median(ClinicMeta$Days_to_death)
ClinicMeta$Progress = ifelse(ClinicMeta$Days_to_death < 15,'severe','moderate')
table(ClinicMeta$Progress)

Endo.m$Donor = factor(gsub('_[12345]','',Endo.m$donor))
Endo.m$Age = factor(ClinicMeta[Endo.m$Donor,'Age'],levels = c("30-35","40-45","50-55","55-60","60-65","65-70","75-80","80-85",">89"),ordered = T)
Endo.m$Progress = factor(ClinicMeta[Endo.m$Donor,'Progress'],levels = c('moderate','severe'),ordered = T)
Endo.m$IMV = ClinicMeta[Endo.m$Donor,'IMV_days']

celltype = Endo.m$SubCluster
celltype = celltype[!grepl("doublet|mix",tolower(celltype))] %>% droplevels()
table(celltype)

Endo.m <- Endo.m[,names(celltype)]
Endo.m$celltype <- celltype

Endo.m = ScaleData(Endo.m) %>% FindVariableFeatures() %>% RunPCA()
ElbowPlot(Endo.m)

Endo.m = RunUMAP(Endo.m,dims = 1:20)
DimPlot(Endo.m,group.by = 'celltype')

FeaturePlot(Endo.m,features = c('hUSI','nCount_RNA'))

cor(Endo.m$hUSI,colSums(Endo.m@assays$RNA@data))


### Capillary Aerocytes
library(monocle)
library(clusterProfiler)
VlnPlot(Endo.m,features = c("CA4","EDNRB","FIBIN","TBX2","CDKN2B","RPRML","CHST1","APLN"),group.by = 'celltype')


subEndo = subset(Endo.m,celltype == 'Capillary Aerocytes')
subEndo = FindVariableFeatures(subEndo) 

expr_matrix <- Counts[,colnames(subEndo)] %>% as.matrix() %>% as('sparseMatrix')
p_data <- subEndo@meta.data
f_data <- data.frame(gene_short_name = rownames(expr_matrix),row.names = rownames(expr_matrix))
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds) %>% estimateDispersions()

ordergene <- VariableFeatures(subEndo)
cds <- setOrderingFilter(cds, ordergene)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree',norm_method = 'log')
cds <- orderCells(cds)

plot_cell_trajectory(cds,color_by="Pseudotime", size=2,show_backbone=TRUE) 

### marker expression
subEndo$State = factor(ifelse(cds[,colnames(subEndo)]$State == 1,'Root',
                                  ifelse(cds[,colnames(subEndo)]$State == 2,'Normal','Senescent')),
                           levels = c('Root','Normal','Senescent'))

markers = c('VWF','NANOG','CD44','CD9','CDK6','CDKN1A','CDKN2A','CDKN2B')
# my_genes <- cds[row.names(subset(fData(cds), gene_short_name %in% markers)),]
# plot_genes_violin(my_genes,ncol=2, min_expr=0.1)
mat = AverageExpression(subEndo,features = markers,group.by = 'age_state',assays = 'RNA',slot = 'data')$RNA
pheatmap::pheatmap(mat,scale = 'row',cluster_cols = F)

### branch degs
BEAM_res <- BEAM(cds, branch_point = 1, cores = 2)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

genes = BEAM_res[BEAM_res$qval<1e-20,'gene_short_name']

heatmap = plot_genes_branched_heatmap(cds[genes,],
                                      branch_point = 1,
                                      num_clusters = 2,
                                      cores = 1,
                                      use_gene_short_name = T,
                                      show_rownames = F,
                                      return_heatmap = T,
                                      branch_labels = c('Normal','Senescence'))
genes_cluster = split(heatmap[["ph"]][["tree_row"]][["labels"]],heatmap[["annotation_row"]][["Cluster"]])
names(genes_cluster) = paste('Cluster',names(genes_cluster),sep='_')

####### validate state marker
marker_set = genes_cluster
names(marker_set) = c('Normal','Senescent')
lapply(marker_set,length)

### bulk degs
bulk = read.csv('/mnt/data1/wangj/AgingScore/GSE206677_humanECs/GSE206677_fpkm.txt',sep = '\t',header = T)
bulk = bulk[bulk$gene_symbol != '',]
bulk = bulk[-1]
exp_bulk = bulk %>%
              mutate(sum = rowSums(.[-1])) %>%
              filter(sum != 0) %>%
              group_by(gene_symbol) %>% 
              mutate(max = max(sum)) %>% 
              filter(sum == max) %>% 
              dplyr::select(-c(sum,max)) %>%
              column_to_rownames("gene_symbol") %>%
              data.matrix

geneID <- intersect(rownames(exp_bulk),rownames(subEndo@assays$RNA))
exp_bulk <- exp_bulk[geneID,c(7:12)]

deg.tab <- NULL
pheno <- rep(c(0,1),each = 3)
for(i in 1:nrow(exp_bulk)){
  lm.model <- lm(exp_bulk[i,]~pheno)
  line <- c(summary(lm.model)$coefficients[2,4],summary(lm.model)$coefficients[2,1])
  deg.tab <- rbind(deg.tab,line)
}
rownames(deg.tab) <- rownames(exp_bulk)
colnames(deg.tab) <- c("pval","coefficients")
deg.tab <- as.data.frame(deg.tab)
deg.tab$p.adjust <- p.adjust(deg.tab$pval,method="BH")
deg.tab <- deg.tab[deg.tab[,3]<0.05,]

### enrich single cell state marker in bulk degs
enrich_geneList=list()
enrich_reList=list()
N = nrow(exp_bulk)
for(de in c("up","down")){
  if(de == "up"){
    bk = rownames(filter(deg.tab,coefficients>0)) 
  }
  else {
    bk = rownames(filter(deg.tab,coefficients<0))
  }
  M = length(bk)
  for (state in names(marker_set)) {
    genes = marker_set[[state]]
    n = length(genes)
    k = length(intersect(bk,genes))
    enrich_reList[[paste(state,de,sep = '_')]] = phyper(k-1,M, N-M, n, lower.tail=FALSE)
    enrich_geneList[[paste(state,de,sep = '_')]][['bulk']] = bk
    enrich_geneList[[paste(state,de,sep = '_')]][['single-cell']] = genes
  }
}
enrich_reList


### enrichment of state marker
subEndo <- ScaleData(subEndo)
Idents(subEndo) <- 'age_state'
State_markers_list <- FindAllMarkers(subEndo,only.pos = T,assay = 'RNA')
State_markers = filter(State_markers,p_val_adj < 0.05)
State_markers = State_markers[order(State_markers$cluster,State_markers$avg_log2FC,decreasing = T),]
State_markers_list = split(State_markers$gene,State_markers$cluster)
lapply(State_markers_list, length)

### enrich branch degs go pathway 
pathways = read.gmt("/mnt/data1/wangj/GeneSets/c5.go.v2023.1.Hs.symbols.gmt")
pathways = pathways[grep('GOBP',pathways$term),]
genes =  genes_cluster$Cluster_1
cluster_enrich = enricher(gene = genes, 
                          universe = rownames(cds),
                          TERM2GENE = pathways,
                          pvalueCutoff=0.05)
cluster_enrich@result = filter(cluster_enrich@result,qvalue < 0.05)

df_plot = cluster_enrich@result
df_plot$class = substr(df_plot$ID,1,4)
df_plot = filter(df_plot,qvalue < 0.05)
df_plot = df_plot[order(df_plot$Count,decreasing = T),]
df_plot$label = tolower(df_plot$ID) %>% gsub('^go[bp|cc|mf]*_','',.) %>% gsub('^kegg_','',.) %>% gsub("_",' ',.)
df_plot$label


### load WTA ROI
library(stringr)
library(stringi)
setwd("/mnt/data1/wangj/AgingScore/GSE163530_COVID-19/GSE162911_GeoMx/")
wta.counts <- read.csv("Broad-COVID_WTA_Q3Norm_TargetCountMatrix.txt", row.names=1,header = T,sep = '\t')
wta.assay <- CreateAssayObject(wta.counts) %>% NormalizeData()

df_segments <- read.csv("Broad-COVID_WTA_SegmentProperties.txt", row.names=1,header = T,sep = '\t')
rownames(df_segments) <- str_replace_all(rownames(df_segments), '-', '.')

df_tissue <- read.table("annotation_file_wta.txt", sep = "\t",row.names=1,header = T)
rownames(df_tissue) <- str_replace_all(rownames(df_tissue), '-', '.')

meta <- cbind(df_segments,df_tissue)
meta$donor <- as.character(unlist(lapply(strsplit(meta$scan.name,"-"),"[",1)))

meta$condition <- substr(meta$donor,1,1)
meta$donor <- stri_replace_all_regex(meta$donor,c("C01","C2","C3","S01","S02","S03","S09","S10","S11","S16","S18","S28"),
                                     c("D22","D23","D24","D18","D19","D20","D21","D8","D9","D10","D11","D12"),vectorize = F)

meta$Progress = factor(ClinicMeta[meta$donor,'Progress'],levels = c('moderate','severe'))
meta = meta[complete.cases(meta$donor),]

Bulk <- as.matrix(wta.assay@data[,rownames(meta)])
results = corto::ssgsea(Bulk, State_markers_list)
colnames(results) = colnames(Bulk)




