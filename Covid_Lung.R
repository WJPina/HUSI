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
anndata = sc$read_h5ad("Covid_Lung_Endo.h5ad")
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
VlnPlot(Endo.m,features = c("APLN","EDNRB", "TBX2" ,"FOXP2"),group.by = 'celltype')

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

plot_cell_trajectory(cds,color_by="State", size=2,show_backbone=TRUE) 

### branch degs
BEAM_res <- BEAM(cds, branch_point = 1, cores = 10)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

genes = BEAM_res[BEAM_res$qval<1e-10,'gene_short_name']

heatmap = plot_genes_branched_heatmap(cds[genes,],
                                      branch_point = 1,
                                      num_clusters = 2,
                                      cores = 1,
                                      use_gene_short_name = T,
                                      show_rownames = F,
                                      return_heatmap = T)
genes_cluster = split(heatmap[["ph"]][["tree_row"]][["labels"]],heatmap[["annotation_row"]][["Cluster"]])
names(genes_cluster) = paste('Cluster',names(genes_cluster),sep='_')

### enrich branch degs go pathway 
pathways = read.gmt("/mnt/data1/wangj/GeneSets/KEGG.gmt")
# pathways = pathways[grep('GOBP',pathways$term),]
cluster_enrich = enricher(gene = genes_cluster$Cluster_2, 
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

####### validate state marker
marker_set = genes_cluster
names(marker_set) = c('Senescent','Stem-like')
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


### enrichment of branch start cell group
Pseudotime = cds$Pseudotime
names(Pseudotime) <- colnames(cds)
branches = split(Pseudotime,cds$State)
State_1_cells = names(sort(branches[[1]],decreasing = F))[1:100]
State_2_cells = names(sort(branches[[2]],decreasing = T))[1:100]
State_3_cells = names(sort(branches[[3]],decreasing = T))[1:100]

cell_use = c(State_1_cells,State_2_cells,State_3_cells)
label = rep(c('Normal','Stem-like','Senescent'),each = 100)
obj = subEndo[,cell_use]
obj$State = factor(label)
Idents(obj) <- 'State'
State_markers = FindAllMarkers(obj,assay = 'RNA',logfc.threshold = 1,only.pos = T) 
State_markers = filter(State_markers,p_val_adj < 0.05)
State_markers = State_markers[order(State_markers$cluster,State_markers$avg_log2FC,decreasing = T),]
State_markers_list = split(State_markers$gene,State_markers$cluster)
lapply(State_markers_list, length)

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

### marker expression
subEndo$age_state = factor(ifelse(cds[,colnames(subEndo)]$State == 1,'Normal',
                                  ifelse(cds[,colnames(subEndo)]$State == 2,'Stem-like','Senescent')),
                           levels = c('Normal','Stem-like','Senescent'))

markers = c('NANOG','CD44','CD9','CDK6','CDKN1A','CDKN2A','CDKN2B')
df_plot = AverageExpression(subEndo,features = markers,group.by = 'age_state',slot = 'data',assays = 'RNA')$RNA %>% 
          as.matrix()
df_plot = df_plot[rowSums(df_plot)!=0,]
pheatmap::pheatmap(df_plot,scale = 'row',cluster_cols = F)


### validat in mouse stem Cap data
setwd('/home/wangjing/wangj/AgingScore/GSE211335_mouseCap')
files = list.files('.')
files = files[grep('h5',files)]
sceList <-lapply(files,function(file){CreateSeuratObject(counts = Read10X_h5(file),
                                                         project = substr(file,12,17),
                                                         min.cells = 0,min.features = 0)})
meta = read.csv('GSE211335_Endothelial_seurat_metadata.csv',row.names = 1)
rownames(meta) = paste(rownames(meta),'-1',sep = '')
sce <- merge(x = sceList[[1]],y = sceList[c(2,3)],add.cell.ids = c('PoolA','PoolB','PoolC'))
sce_endo = sce[,rownames(meta)]
sce_endo$cluster = factor(paste('Cluster',meta$seurat_clusters,sep = ''))
Idents(sce_endo) <- 'cluster'
sce_endo = NormalizeData(sce_endo)
stem_markers = FindMarkers(sce_endo,ident.1 = 'Cluster7',assay = 'RNA',only.pos = T)

stem_markers = filter(stem_markers,p_val_adj<0.05)
stem_markers = stem_markers[order(stem_markers$avg_log2FC,decreasing = T),]

library(biomaRt)
mouse <- useMart('ensembl',dataset = "mmusculus_gene_ensembl")
human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl")
m2h.g <- getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",
                values = rownames(stem_markers)[1:50],mart = mouse,
                attributesL = c("hgnc_symbol"),
                martL = human,uniqueRows = T)







