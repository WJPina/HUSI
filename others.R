### Georgilis2018
Georgilis2018 = CreateSeuratObject(
    fread("Georgilis2018/valid_TPM_dataset.tsv") %>% column_to_rownames("Gene Name") %>% data.matrix, 
    meta.data = fread("Georgilis2018/filereport_read_run_PRJNA395386_tsv.txt", header = T) %>% column_to_rownames("sample_title")
)
Georgilis2018 = subset(Georgilis2018,read_count>1e6)
Georgilis2018 = NormalizeData(Georgilis2018, scale.factor = 1e4)
hSI = GetAssayData(Georgilis2018) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}

dat = calc_scores(ScoreList,Georgilis2018)
dat_Geo = dat %>% 
        rownames_to_column('sample') %>%
        mutate(condition=case_when(
        grepl("Allstars|Water|Hs|Extras[5-8]", sample) ~ "Senescence", 
        grepl("Noninduced|Extras[1-4]", sample) ~ "Growing", 
        TRUE ~ "other")) %>% 
        mutate(hSI=hSI) %>%
        filter( !condition %in% "other" ) %>% 
        mutate( condition=as.factor(condition))
dat_Geo = dat_Geo[,complete.cases(t(dat_Geo))]
auc_Geo <- calc_auc(dat_Geo,'marker',ScoreList)

auc %>%
    melt(value.name = "Accuracy") %>%
    mutate( feature=forcats::fct_reorder(Var1, Accuracy, mean, .desc = T) ) %>%
    Rmisc::summarySE(measurevar = "Accuracy", groupvars = "feature") %>%
    {
        ggplot(data = ., aes(feature, Accuracy, fill=feature)) +
        geom_bar(stat = "identity", position = position_dodge(), width = 0.4) +
        geom_errorbar(aes(ymin = Accuracy - sd, ymax = {Accuracy + sd} %>% ifelse(. >1, 1, .)),
        width = 0.2, position = position_dodge(0.9)) +
        labs(x=NULL, y="AUC") +
        theme(legend.position = "none")+
        theme_classic()
    }

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