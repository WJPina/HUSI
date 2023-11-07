


### create lost RNA-seq
setwd("~/wangj/AgingScore/Data/Bulk_RNA-seq/")
RNAList = list()
RNAList[['GSE60340']] <- fread("GSE60340_induced_sene_TPM.csv") %>% column_to_rownames("Gene Name") %>% as.matrix
mat_GSE60340 <- do.call(cbind,list.files("GSE130306/", "GSM", full.names = T) %>% 
                              lapply(function(x) {fread(x) %>% filter( !duplicated(gene_name) ) %>% column_to_rownames("gene_name")})) %>% as.matrix
colnames(mat_GSE60340) <- list.files("GSE130306/", "GSM", full.names = T) %>% gsub(".*RNAseq_", "", .) %>% gsub(".txt.gz", "", .)
RNAList[['GSE60340_OIS']] <- mat_GSE60340[,grep('OIS',colnames(mat_GSE60340))]
RNAList[['GSE60340_RS']] <- mat_GSE60340[,grep('RS',colnames(mat_GSE60340))]

for(i in 1:length(RNAList)){
  mat = RNAList[[i]]
  rates = list() 
  for (rate in lost_rate) {
    if(rate == 0){tmat = mat}
    else{
      corrdinates = expand.grid(1:nrow(mat), 1:ncol(mat))
      lost = corrdinates[sample(1:nrow(corrdinates),round(nrow(corrdinates)*rate),replace = F),]
      tmat = mat
      for(j in 1:nrow(lost)){
        tmat[lost[,1][j],lost[,2][j]] <- 0
      }
    }
    rates[[paste('rate_',rate,sep = '')]] = tmat %>% apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )})
  }
  Lost_ScoreList[[names(RNAList)[i]]] = rates
}


library(msigdbr)
set.seed(233)
msigdb_all = msigdbr()
pathways_sene = msigdb_all %>% 
                  filter( (gs_cat %in% "H" ) ) %>% 
                  mutate( gs_name=gsub("HALLMARK_", "", gs_name) ) %>% 
                  plyr::dlply(.variables = "gs_name", .fun = function(x) x$gene_symbol )
res_fgsea_sene = fgseaMultilevel(pathways = pathways_sene, stats = sort(mm_l2$w,decreasing=T), nPermSimple = 10000)

df = filter(res_fgsea_sene,pval < 0.05) %>% 
  .[order(.$NES,decreasing=T),]
df[,c('pathway','NES')]


# ### Enge2017 GSE81547 human pancreas organisms age
# Enge2017List = list.files("Enge2017", pattern = "*.gz",full.names = T)
# meta = read.csv("Enge2017/metadata.txt",row.names = 1)
# table(meta$DONOR_AGE)
# meta = filter(meta, DONOR_AGE %in% c(1,54))

# samples = lapply(strsplit(basename(Enge2017List),"_"),'[',1) %>% unlist 
# files = Enge2017List[samples %in% meta$Sample.Name]

# scdat = data.frame(gene = '')
# for(file in files){
#     scol = fread(file) 
#     colnames(scol) = c('gene',strsplit(basename(file),"_")[[1]][1])
#     scdat = right_join(scdat,scol,by = 'gene')           
# }
# scdat = column_to_rownames(scdat,"gene") %>% data.matrix
# dim(scdat)
# Enge2017 = CreateSeuratObject(scdat) %>% 
#                 # NormalizeData(normalization.method = 'RC',scale.factor = 1e6)
#                 NormalizeData()
# Enge2017$Condition = ifelse(meta$DONOR_AGE %in% c(54), 'Senescence', 'Growing')
# hUSI = GetAssayData(Enge2017) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}

# dat_laf = calc_scores(lafSet,Enge2017,'laf')
# dat_marker = calc_scores(SenMarkers,Enge2017,'marker')
# dat_ssgsea = calc_scores(EnrichSet,Enge2017,'ssgsea')

# Enge2017List = list(laf=dat_laf,marker=dat_marker,ssgsea=dat_ssgsea)

# Results$Enge2017 = Enge2017List
# Results$Enge2017$hUSI = hUSI

### melanoma heatmap
mat <- exp[unlist(gene_set),rownames(meta)]
mat=mat[!apply(mat, 1, sd)==0,]
mat=Matrix::t(scale(Matrix::t(mat),center=TRUE))
mat=mat[is.na(row.names(mat)) == FALSE,]
mat[is.nan(mat)] = 0
mat[mat>2] = 2
mat[mat< -2] = -2


library(ComplexHeatmap)
library(circlize)
top_anno <- HeatmapAnnotation(df = data.frame(Condition = meta$condition),
                              col = list(Condition = c('young'='#007E99','senescent' = '#FF745A')),
                              show_legend = F)
left_anno = rowAnnotation(df = data.frame('State'= rep(c("Cycling","Moderate senescent","Senescent"),times=lapply(gene_set,length))),
                          show_legend = T,
                          col = list(State = bar))
col <- colorRamp2(c(-2,0,2), c("blue","white", "red"), space = "LAB")
png('/home/wangjing/codebase/HUSI/Figures/Melanoma_microarray.png',width = 1000,height = 800,res= 200)
Heatmap(mat,
        show_column_names = F,
        show_row_names = F,
        row_title = NULL,
        col = col,
        cluster_rows = T,
        cluster_row_slices = FALSE,
        cluster_columns = F,
        top_annotation = top_anno,
        left_annotation = left_anno,
        row_names_gp = gpar(fontsize = 12),
        column_split = c(rep('young',4),rep('senescent',4)),
        row_split = rep(c("Cycling","Moderate senescent","Senescent"),times=lapply(gene_set,length)))
dev.off()


png('/home/wangjing/wangj/codebase/HUSI/Figures/model/compare_marker.png',width = 4000,height = 2500,res = 400)
df_plot %>%
  mutate(dataset = factor(lapply(strsplit(rownames(df_plot),'\\.'),'[',1) %>% unlist)) %>%
  mutate(condition = recode(dataset,'Teo2019'='OIS','Tang2019'='RS & IRIS','Aarts2017'='OIS','Georgilis2018'='OIS')) %>%
  reshape2::melt(value.name = "AUC",variable.name = 'Method',ids = 'dataset') %>%
  group_by(dataset) %>%
  mutate(method=forcats::fct_reorder(Method, AUC, mean, .desc = T) ) %>% 
  ggplot(aes(method, AUC)) + 
  geom_boxplot(aes(color=method), outlier.shape = NA, width=0.5) + 
  geom_jitter(aes(color=method), width = 0.2,size = 0.1) + 
  labs(x='marker', y="AUC") +
  theme_classic()+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(size = 16))+
  facet_wrap(~dataset,ncol = 2,scales = "free_y")
# scale_color_manual(values = c('#e63946','#023e8a','#457b9d','#f77f00','#fcbf49','#eae2b7'))
dev.off()


### age state fraction in TCGA SKCM patients
library(ComplexHeatmap)
library(circlize)
marker_plot = unlist(lapply(marker_set, function(x) {rownames(x[order(x$avg_log2FC,decreasing = T),][0:200,])})) %>% unique
genes = intersect(rownames(Bulk),marker_plot)
mat <- as.matrix(Bulk[genes,])
mat <- apply(mat,1,scale)
mat <- t(mat)
rownames(mat) <- genes
colnames(mat) <- colnames(Bulk)

Frac_scale = data.frame(Frac)[colnames(mat),]
### scale Frac into 0-1 by column
Frac_scale <- apply(Frac_scale,2,function(x) {(x-min(x)) / (max(x)-min(x))})

top_anno <- HeatmapAnnotation(df = Frac_scale,
                              col = list(Cycling = colorRamp2(c(0, 1), c("white", "#5cb85c")),
                                         Transition = colorRamp2(c(0, 1), c( "white", "#428bca")),
                                         Senescent= colorRamp2(c(0, 1), c("white", "#d9534f"))),
                              show_legend = F)
col <- colorRamp2(c(-2,0,2), c("blue","white", "red"), space = "LAB")
png('/home/wangjing/wangj/codebase/HUSI/Figures/TCGA_SKCM_heatmap.png',width = 1500,height = 1000,res= 200)
Heatmap(mat,
        show_column_names = F,
        show_row_names = F,
        row_title = NULL,
        col = col,
        cluster_rows = T,
        cluster_row_slices = FALSE,
        cluster_columns = T,
        column_km=3,
        row_km = 3,
        top_annotation = top_anno,
        heatmap_legend_param = list(title = "Expression",
                                         at = c(-2,-1,0,1,2),
                                         labels = c("-2", "-1", "0","1", "2"),
                                         labels_gp = gpar(fontsize = 12),
                                         title_gp = gpar(fontsize = 12, fontface = "bold"),
                                         legend_width = unit(30, "mm")))
dev.off()


### load spatial ROI data
setwd("/mnt/data1/wangj/AgingScore/GSE163530_COVID-19/GSE162911_GeoMx/")
############################## Load WTA data ###################################
wta.counts <- read.csv("Broad-COVID_WTA_Q3Norm_TargetCountMatrix.txt", row.names=1,header = T,sep = '\t')
df_segments <- read.csv("Broad-COVID_WTA_SegmentProperties.txt", row.names=1,header = T,sep = '\t')
rownames(df_segments) <- str_replace_all(rownames(df_segments), '-', '.')
 
df_tissue <- read.table("annotation_file_wta.txt", sep = "\t",row.names=1,header = T)
rownames(df_tissue) <- str_replace_all(rownames(df_tissue), '-', '.')

meta <- cbind(df_segments,df_tissue) %>% rename("tissue" = "Primary_Morph")
meta$donor <- as.character(unlist(lapply(strsplit(meta$scan.name,"-"),"[",1)))

### Split into Patient (S) and Control (C) groups
meta$group <- substr(meta$donor,1,1)
meta$donor <- stri_replace_all_regex(meta$donor,c("C01","C2","C3","S01","S02","S03","S09","S10","S11","S16","S18","S28"),
                                     c("D22","D23","D24","D18","D19","D20","D21","D8","D9","D10","D11","D12"),vectorize = F)
### SARS-Cov-2 signature score
Vss <- read.csv('SARS-CoV-2 signature score.csv',header = T)
Vss$barcodekey <- gsub('-','\\.',Vss$barcodekey)
wta <- wta.counts[,Vss$barcodekey]
meta <- meta[Vss$barcodekey,]
meta$Virus_score <- Vss$Virus


###### decovolution in Covid-19 bulk data
### load bulk data
Bulk <- read.table("/mnt/data1/wangj/AgingScore/GSE150316_COVID-19_bulk/GSE150316_RawCounts_Final.txt")
Bulk <- Bulk[,grepl('lung',colnames(Bulk))]

clinicalMatrix <- read.csv("/mnt/data1/wangj/AgingScore/GSE150316_COVID-19_bulk/clinical.csv",header=TRUE,stringsAsFactors = FALSE,row.names = 1)
clinicalMatrix <- clinicalMatrix[complete.cases(clinicalMatrix[, 'DTD']), ] 
rownames(clinicalMatrix) <- gsub('Case ','case',rownames(clinicalMatrix))
cases = unlist(lapply(strsplit(colnames(Bulk),'\\.'),'[',1))

Bulk <- Bulk[,cases %in% rownames(clinicalMatrix)]
dim(Bulk)


### SARS-Cov-2 signature score
Vss <- read.csv('SARS-CoV-2 signature score.csv',header = T)
Vss$barcodekey <- gsub('-','\\.',Vss$barcodekey)
wta <- wta.counts[,Vss$barcodekey]
meta <- meta[Vss$barcodekey,]
meta$Virus_score <- Vss$Virus


### rename cell type
# "B+Plasma","T+NK","Fibroblasts","Epithelial"

# celltype = lung.obj$SubCluster
# levels(celltype)[which(levels(celltype)%in% c("AT1","AT2","KRT8+ PATS/ADI/DATPs","Secretory"))] <- 'Epithelial'
# levels(celltype)[which(levels(celltype)%in% c("Proliferative fibroblast","Fibroblast","Myofibroblast"))] <- 'Fibroblasts'
# levels(celltype)[which(levels(celltype)%in% c("NK cells","NK/NKT","B cells","Plasma cells PRDM1/BLIMP hi","Plasma cells PRDM1/BLIMP int","Plasmablasts"))] <- 'B cells'
# levels(celltype)[which(levels(celltype)%in% c("CD8+ T cells","CD4+ T cells metabolically active","CD4+ Treg","Doublets CD4 T/ CD8 T/ NK cells"))] <- 'T+NK cells'
# levels(celltype)[which(levels(celltype)%in% c("Capillary Aerocytes","Lymphatic EC","Capillary 1","Capillary","Aerocytes","Capillary 2","Artery EC","Vein EC"))] <- 'Endothelial'
# celltype = celltype[!grepl("doublet",tolower(celltype))] %>% droplevels()
# table(celltype)



### expression of SASP
sasp = read.csv("/mnt/data1/wangj/AgingScore/AgingScorePro/SASP.csv")[[1]]
sasp = sasp[sasp %in% rownames(Endo.m)]

cor_list = list()
for(gene in sasp){
  cor_list[[gene]] = cor(Endo.m@assays$RNA@data[gene,],Endo.m$hSI)
}

sort(unlist(cor_list),decreasing = T)



# Endo.m <- FindVariableFeatures(object = Endo.m, selection.method = 'vst')
# Endo.m <- ScaleData(Endo.m,features=rownames(Endo.m))
# Endo.m <- RunPCA(Endo.m,npcs = 30,verbose = F,features=rownames(Endo.m))
# ElbowPlot(Endo.m,ndims=30,reduction="pca")
# 
# Endo.m <- RunTSNE(Endo.m, reduction = "pca", dims = 1:30)
# Endo.m <- RunUMAP(Endo.m,reduction = "pca",dims = 1:30)
# DimPlot(Endo.m,group.by = 'celltype')


EDNRB, TBX2, FOXP2, CLEC4E, SPON2, PRKG, CHRM2, S100A4, EDA, HPGD, SOSTDC1 and do not express Weibel–Palade body-associated genes.
endothelial-specific Weibel-Palade bodies (vWF, SELP, EDN)


APLN,APLNR, EPCR, CD93 
genes = c('THSD4')


### enrichment pathway 
enrichplot::cnetplot(cluster_enrich,circular=F,
                     colorEdge = TRUE,showCategory = 5,
                     cex_gene = 0.5,cex_category = 0.8,
                     cex_label_gene = 0.6,node_label = 'gene') + 
  theme(legend.position = 'none')



### single-cell stem-like deg
dir = '/mnt/data1/wangj/AgingScore/GSE211335_mouseCap/'
files = list.files(dir)
files = files[grep('h5',files)]
mouse.capList = lapply(files, function(file){CreateSeuratObject(counts = Read10X_h5(paste(dir,file,sep = '')),min.cells = 0,min.features = 0,project = substr(file,12,17))})
meta = read.csv('/mnt/data1/wangj/AgingScore/GSE211335_mouseCap/GSE211335_Endothelial_seurat_metadata.csv',row.names = 1)
rownames(meta) = paste(rownames(meta),'1',sep = '-')
mouse.cap = merge(x = mouse.capList[[1]],y = c(mouse.capList[[2]],mouse.capList[[3]]),add.cell.ids=c('PoolA','PoolB','PoolC'))
mouse.cap = mouse.cap[,rownames(meta)]
mouse.cap$cluster = factor(meta$seurat_clusters)
mouse.cap = NormalizeData(mouse.cap)
Idents(mouse.cap) <- 'cluster'
marker_7 = FindMarkers(mouse.cap,ident.1 = 7,only.pos = T,assay = 'RNA')
marker_7_use = filter(marker_7,p_val_adj < 0.05&avg_log2FC>0.5)


### enrich single cell state marker in bulk degs
enrich_geneList=list()
enrich_reList=list()
N = nrow(exp_bulk)
for(de in c("Senescence",'Stemness')){
  if(de == "Senescence"){
    bk = rownames(filter(deg.tab,coefficients>0)) 
  }
  else {
    bk = homologene(rownames(marker_7_use),inTax = 10090, outTax = 9606)[['9606']]
  }
  M = length(bk)
  for (state in names(marker_set)) {
    genes = marker_set[[state]]
    n = length(genes)
    k = length(intersect(bk,genes))
    enrich_reList[[paste(state,de,sep = '_')]] = phyper(k-1,M, N-M, n, lower.tail=FALSE)
    enrich_geneList[[paste(state,de,sep = '_')]][['bk']] = bk
    enrich_geneList[[paste(state,de,sep = '_')]][['fate']] = genes
  }
}


pwc <- df_plot %>% 
        group_by(Progress) %>%
        wilcox_test(ssGSEA ~ State, paired = F)%>%
        add_significance("p") %>%
        add_xy_position(x = "Progress")
pwc

# features = rowMeans(Endo.m@assays$RNA@data)
# features = names(features)[features > quantile(features,0.75)]
# AgeScore  = Endo.m@assays$RNA@data[features,] %>% {apply( ., 2, function(z) {cor(z, mm_l2$w[ features ], method="sp", use="complete.obs" )})}
# Endo.m$hUSI <- AgeScore[colnames(Endo.m)]



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
sce_endo$condition = factor(meta$Sample)
sce_endo = subset(sce_endo,cluster == 'Cluster2')

Idents(sce_endo) <- 'condition'
sce_endo = NormalizeData(sce_endo)
stem_markers = FindAllMarkers(sce_endo,
                              assay = 'RNA',
                              only.pos = T,logfc.threshold = 0.5)

stem_markers = filter(stem_markers,p_val_adj<0.05)
stem_markers = stem_markers[order(stem_markers$cluster,stem_markers$avg_log2FC,decreasing = T),]

# mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
convert_mouse_to_human <- function(gene_list){
  output = c()
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  
  return (output)
}
markers = convert_mouse_to_human("DPT")

df_plot = AverageExpression(subEndo,features = markers,group.by = 'age_state',slot = 'data',assays = 'RNA')$RNA %>% as.matrix()
df_plot = df_plot[rowSums(df_plot)!=0,]
pheatmap::pheatmap(df_plot,scale = 'row',cluster_cols = F)

VlnPlot(subEndo,features = markers,group.by = 'age_state')


# Pseudotime = cds$Pseudotime
# names(Pseudotime) <- colnames(cds)
# branches = split(Pseudotime,cds$State)
# State_1_cells = names(sort(branches[[1]],decreasing = F))[1:100]
# State_2_cells = names(sort(branches[[2]],decreasing = T))[1:100]
# State_3_cells = names(sort(branches[[3]],decreasing = T))[1:100]
# 
# cell_use = c(State_1_cells,State_2_cells,State_3_cells)
# label = rep(c('Normal','Anti_senescence','Senescent'),each = 100)
# obj = subEndo[,cell_use]
# obj$State = factor(label)
# Idents(obj) <- 'State'
# State_markers = FindAllMarkers(obj,assay = 'RNA',logfc.threshold = 1,only.pos = T) 
# State_markers = filter(State_markers,p_val_adj < 0.05)
# State_markers = State_markers[order(State_markers$cluster,State_markers$avg_log2FC,decreasing = T),]
# State_markers_list = split(State_markers$gene,State_markers$cluster)
# lapply(State_markers_list, length)





### new platform
p1.2 <- ArrayList [["GSE16058"]][[1]] %>% 
  pData %>% 
  mutate(condition=gsub("growth status: ", "", characteristics_ch1.3) ) %>% 
  mutate(passage = `passage:ch1`) %>%
  inner_join(
    data.frame(sene_score=ScoreList[["s_16058"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("Growing","Senescent")),] %>% 
  mutate(condition = factor(condition,levels = c("Growing","Senescent"),ordered = T)) %>% 
  ggplot(aes(condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = condition,shape =passage),width = 0.1) +
  geom_signif(comparisons = list(c("Growing", "Senescent")),test = "t.test",test.args = c("less"), map_signif_level = T) + 
  mytheme()+
  guides(color = 'none')+
  ylab('hUSI')+
  xlab("RS of HMEC")



### Capillary Aerocytes
library(monocle)
library(clusterProfiler)
VlnPlot(covid.m,features = c("CA4","EDNRB","FIBIN","TBX2","CDKN2B","RPRML","CHST1","APLN"),group.by = 'celltype')


subEndo = subset(covid.m,celltype == 'Capillary Aerocytes')
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


# ### Using PHATE to predict the Transitional process of tumor aging
# Epi.data <- t((GetAssayData(EpiExp.m)))
# Epi.data <- as.data.frame(Epi.data)
# Epi.phate <- phate(Epi.data)
# branch <- EpiExp.m$age_class
# 
# EpiExp.m[['phate']] <- CreateDimReducObject(Epi.phate$embedding,key="phate_")
# DimPlot(EpiExp.m, reduction = 'phate', group.by = 'State',label=F,pt.size=2)

use_python("/home/tools/anaconda3/envs/sc/bin/python3", required = T)



### ranking the genes based on each gene expression correlation with the aging score
corAging <- function(x,agingScore){cor <- cor(x,agingScore);cor}
cor.genes <- apply(GetAssayData(subEndo),1,corAging,agingScore=subEndo$hUSI)
cor.genes[is.na(cor.genes)] <- 0
features <- names(cor.genes)[order(abs(cor.genes),decreasing=TRUE)][1:1500]
features 

subEndo.matrix <- as.matrix(GetAssayData(subEndo))[features,]
subEndo[['AgingExp']] <- CreateAssayObject(subEndo.matrix)
DefaultAssay(subEndo) <- "AgingExp"
subEndo$batch <- rep("batch",ncol(subEndo))

### Running the ICAnet
library(ICAnet)
source('~/wangj/codebase/HUSI/getPPI_String.R')
Ica.epi <- ICAcomputing(subEndo,ICA.type="JADE",RMT=TRUE,two.stage=FALSE)
rownames(Ica.epi$ica.pooling) <- features
# Ica.filter <- CrossBatchGrouping(Ica.epi$ica.pooling)
PPI <- readRDS('PPI_subEndo.rds')
# PPI <- getPPI_String(subEndo,species=9606)
subEndo <- RunICAnet(subEndo,Ica.epi$ica.pooling,PPI.net = PPI,scale=FALSE,
                     ModuleSignificance = FALSE,cores = 1,aucMaxRank=500)

subEndo <- RunPCA(subEndo, features = rownames(subEndo), npcs = 30, 
                  reduction.name = "IcaNetPCA", reduction.key = "IcaNetPCA_", verbose = F)
ElbowPlot(subEndo,reduction = "IcaNetPCA")

subEndo <- RunUMAP(subEndo, reduction = "IcaNetPCA", dims = 1:30, 
                   reduction.name = "IcaNetUAMP",  reduction.key = "IcaNetUAMP_")
subEndo <- FindNeighbors(subEndo, reduction="IcaNetPCA",dims = 1:30,graph.name = "IcaNet_snn") %>%
  FindClusters(algorithm=2,graph.name = "IcaNet_snn",resolution = 0.3)
DimPlot(subEndo, reduction = 'IcaNetUAMP', pt.size = 1) +
  FeaturePlot(subEndo,features = 'hUSI', reduction = 'IcaNetUAMP')

VlnPlot(subEndo,features = 'hUSI')



### trajectory plot
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_new/CA_fate.png',width =800,height = 600,res = 300)
plot_cell_trajectory(cds,color_by="State", cell_size=0.5,show_backbone=TRUE) + 
  scale_color_manual(values = c("#2a9d8f","#e9c46a",'#264653'))+
  theme(legend.position = 'none')
# plot_cell_trajectory(cds,color_by="hUSI", cell_size=0.5,show_backbone=TRUE)+
#   scale_colour_gradientn(colors = colorRampPalette(c( 'blue',"white",'red'))(62))
dev.off()


### branch gene heatmap
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_new/CA_fate_genes.png',width = 1300,height = 1500,res = 300)
plot_genes_branched_heatmap(cds[genes,],
                            branch_point = 1,
                            num_clusters = 2,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = F,
                            hmcols = colorRampPalette(rev(brewer.pal(9, "Spectral")))(62),
                            branch_colors = c("#2a9d8f","#e9c46a",'#264653'),
                            branch_labels = c('Anti-senescence','Senescence'))
dev.off()

bar = c("#2a9d8f","#e9c46a",'#264653')
names(bar) = c('Root','Cell fate 1','Cell fate 2')

### middle group state percentage
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_new/CA_percent.png',width = 800,height = 1000,res = 300)
cds@phenoData@data %>% 
  mutate(Fate = ifelse(State == '1','Root',ifelse(State == '2','Cell fate 1','Cell fate 2'))) %>%
  ggplot()+
  geom_bar(aes(x = Progress,fill = Fate),stat = 'count',position = 'fill')+
  theme_classic()+
  theme(axis.title = element_text(face = 'bold'),
        text = element_text(size = 16),
        legend.position = 'none')+
  # theme(text = element_text(size = 16),
  #       axis.title = element_blank(),
  #       legend.position = c(0.15,1.25),legend.direction = 'horizontal',
  #       plot.margin=unit(c(2,1,1,1),'lines'))+
  scale_fill_manual(values = bar)
dev.off()

### hUSI fate
# cds$Fate = ifelse(cds$State == '1','Root',ifelse(cds$State == '2','Cell fate 1','Cell fate 2'))
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_new/HUSI_fate.png',width = 1000,height = 1200,res = 300)
ggviolin(subEndo@meta.data,
         x="State", y="hUSI", fill = "State", 
         palette =  c("#2a9d8f",'#264653',"#e9c46a"),
         add = "boxplot", 
         add.params = list(fill="white"),
)+
  theme(axis.title = element_text(face = 'bold'),
        legend.position = 'none')+
  stat_compare_means(comparisons = list(c('Normal','Root'),c('Senescent','Noraml'),c('Senescent','Root')),
                     method = 't.test',
                     method.args = list(alternative = "greater"))
dev.off()


### pathway network
library(enrichplot)
cols = mypalette$palette3[1:length(ids)]
# id_1 = df_plot$ID
# id_2 = df_plot$ID
# ids = unique(c(id_1,id_2))
names(cols) = ids
cols

foldChange = BEAM_res[genes,]$qval
names(foldChange) = genes
p = cnetplot(cluster_enrich,
             circular=F,
             colorEdge = F,
             showCategory = 20,
             cex_gene = 0.8,
             cex_category = 0.8,
             foldChange = -log10(foldChange),
             cex_label_gene = 0.6,
             node_label = 'gene') + 
  theme(legend.position = 'none')

p[["layers"]][[2]][["mapping"]][["colour_new"]][[2]][[2]] <- cols[cluster_enrich@result$ID]
p[["data"]][["color"]][which(!(p[["data"]][["name"]] %in% names(foldChange)) == T)] <- cols[cluster_enrich@result$ID]

png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_new/Ant-senescence.png',width = 1000,height = 1000,res = 300)
print(p)
dev.off()

### plot legend
# labels_1 = df_plot$label
# labels_2 = df_plot$label
# labels = unique(c(labels_1,labels_2))
# names(cols) = labels

comm = c("focal adhesion","ecm-receptor interaction",
         "human papillomavirus infection","vascular smooth muscle contraction")
sene = c("regulation of actin cytoskeleton","glycosaminoglycan binding proteins",
         "arrhythmogenic right ventricular cardiomyopathy","hypertrophic cardiomyopathy",
         "dilated cardiomyopathy","platelet activation","apelin signaling pathway")
anti = c("pi3k-akt signaling pathway","cd molecules",
         "rap1 signaling pathway","cell adhesion molecules",
         "amoebiasis","tgf-beta signaling pathway","hippo signaling pathway",
         "axon guidance","african trypanosomiasis" ,"adherens junction",
         "small cell lung cancer", "malaria")

png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_new/legend.png',width = 2000,height = 2000,res = 300)
bar = c(cols[anti],cols[sene],cols[comm])
palette(bar)
plot(y = 1:length(bar),x=rep(0.2,length(bar)), col = 1:length(bar), pch = 19,xlim = c(0,2),cex = 2) 
labels <- names(bar)
text(y = 1:length(bar),x=rep(0.2,length(bar)), labels = labels, pos = 4)
dev.off()

palette("default")




### gene overlap
colorList = list()
colorList[['Anti_senescence_up']] <- c("#FF745A","#e9c46a")
colorList[['Anti_senescence_down']] <- c("#007E99","#e9c46a")
colorList[['Senescent_up']] <- c("#FF745A",'#264653')
colorList[['Senescent_down']] <- c("#007E99",'#264653')

pList=list()
library(ggvenn)
for(set in names(enrich_reList)){
  pList[[set]] <- enrich_geneList[[set]] %>% 
    ggvenn(show_percentage = T,
           show_elements = F,
           label_sep = ",",
           digits = 1,
           set_name_size = 0,
           stroke_color = "white",
           fill_color = colorList[[set]],
           text_size = 4,
           fill_alpha = 0.8)+
    labs(title = paste("Padj value:",as.character(signif(enrich_reList[[set]],2))))+
    theme(plot.title = element_text(hjust = 0.5,vjust = 0,size = 15))+
    scale_y_continuous(limits = c(-1, 1))
  
}
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_new/CA_gene_enrich.png',width = 1200,height = 800,res = 200)
ggarrange(pList[['Anti_senescence_up']],
          pList$Senescent_up,
          pList[['Anti_senescence_down']],
          pList$Senescent_down,
          ncol = 2,nrow = 2,legend = "none")
dev.off()


### ROI bulk ssGSEA
library(rstatix)
names(bar) = c('Normal','Anti_senescence','Senescent')

cbind(meta,t(results[,rownames(meta)]))%>%
  dplyr::select(c('Normal','Senescent','Anti_senescence','Primary_Morph','Progress')) %>%
  reshape2::melt(variable.name = 'State',value.name = 'ssGSEA') %>%
  ggboxplot(x = 'Progress',y = 'ssGSEA',fill = 'State')+
  scale_fill_manual(values = bar)+
  theme_classic()


### Anti_senescence or senescence marker
### stem
VlnPlot(subEndo,features = c('CD9','CD44','CDKN1A','CDKN2A'),group.by = 'age_state',pt.size = 1.5,col = bar)

for(sl in names(EnrichSet)){
  score = GSVA::gsva(Counts_gtex, list(signature=EnrichSet[[sl]]), method="ssgsea", ssgsea.norm = TRUE, verbose = TRUE)[1,]
  scoreList_GTEx[[sl]] = score
}

cormat_gtex <- do.call(data.frame, scoreList_GTEx) %>% cor(method = 'spearman')

for(sl in names(EnrichSet)){
  score = GSVA::gsva(Counts_tcga, list(signature=EnrichSet[[sl]]), method="ssgsea", ssgsea.norm = TRUE, verbose = TRUE)[1,]
  scoreList_TCGA[[sl]] = score
}
cormat_tcga <- do.call(data.frame, scoreList_TCGA) %>% cor(method = 'spearman')
mat = cormat_tcga
png('/home/wangjing/wangj/codebase/HUSI/Figures/model/compare_cor_TCGA.png',width = 2000,height = 1800,res = 350)
col <- colorRamp2(c(-1,0,1), c("#1d3557","white", "#e63946"), space = "LAB",transparency = 0)
Heatmap(mat,
        show_column_names = T,
        show_row_names = T,
        row_title = NULL,
        col = col,
        cluster_rows = T,
        cluster_columns = T,
        row_names_gp = gpar(fontsize = 14),
        heatmap_legend_param = list(title = "Spearman\ncofficient"),
        rect_gp = gpar(col = "white", lwd = 1),
        # column_title = "Scores of 7123 samples in GTEx")
        column_title = "Scores of 10495 samples in TCGA")
dev.off()
