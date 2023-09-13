library(reticulate)
library(dplyr)
library(magrittr)
library(Seurat)
library(mclust)

setwd("/mnt/data1/wangj/AgingScore/GSE163530_COVID-19/GSE171668_scnRNA/")

mm_l2 = readRDS("/home/wangjing/wangj/AgingScore/Data/Bulk_TrainModel/mm_l2.rds")
sc <- import("scanpy")
###### preprocess data
### raw data
anndata = sc$read_h5ad("lung.h5ad")
### subset data
anndata_use = anndata[(anndata$obs$method  == 'nuclei') & 
                        (anndata$obs$doublet == F) & 
                        (anndata$obs$Cluster %in% c("Endothelial","B+Plasma","T+NK","Fibroblasts","Epithelial"))]
anndata_use
Counts <- anndata_use$raw$to_adata()$copy()$T$to_df()
Meta <- anndata_use$obs

lung.obj <- CreateSeuratObject(counts=Counts ,meta.data=Meta[,c('Cluster','SubCluster','Viral+','donor','disease')])
lung.obj <- NormalizeData(lung.obj)
lung.obj <- FindVariableFeatures(object = lung.obj, selection.method = 'vst')
lung.obj <- ScaleData(lung.obj,features=rownames(lung.obj))
lung.obj <- RunPCA(lung.obj,npcs = 20,verbose = F,features=rownames(lung.obj)) 
ElbowPlot(lung.obj,ndims=20,reduction="pca")

lung.obj <- RunTSNE(lung.obj, reduction = "pca", dims = 1:20)
DimPlot(lung.obj, reduction = 'tsne', group.by = 'SubCluster',label=1)

### rename cell type
celltype = lung.obj$SubCluster
celltype = celltype[!grepl("mix|doublet",tolower(celltype))] %>% droplevels()
levels(celltype)[which(levels(celltype)%in% c("AT1","AT2","KRT8+ PATS/ADI/DATPs","Secretory"))] <- 'Epithelial'
levels(celltype)[which(levels(celltype)%in% c("Proliferative fibroblast","Fibroblast","Myofibroblast"))] <- 'Fibroblasts'
levels(celltype)[which(levels(celltype)%in% c("NK cells","NK/NKT"))] <- 'NK cells'
levels(celltype)[which(levels(celltype)%in% c("B cells","Plasma cells PRDM1/BLIMP hi","Plasma cells PRDM1/BLIMP int","Plasmablasts"))] <- 'B cells'
levels(celltype)[which(levels(celltype)%in% c("CD8+ T cells","CD4+ T cells metabolically active","CD4+ Treg"))] <- 'T cells'
table(celltype)

lung.obj <- lung.obj[,names(celltype)]
lung.obj$celltype <- celltype
DimPlot(lung.obj,group.by = 'Cluster')
save(lung.obj,file = 'mianCells_9.11.RData')
### explore endothelial cells age state
Endo.m = subset(lung.obj,Cluster == 'Endothelial')
AgeScore  = Endo.m@assays$RNA@data[] %>% {apply( ., 2, function(z) {cor(z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )})}
Endo.m$hSI <- AgeScore[colnames(Endo.m)]
gaussian = Endo.m@meta.data$hSI %>% {log2((1+ .)/(1- .))} %>% Mclust(G=3)
Endo.m$age_class <- gaussian$classification
DimPlot(Endo.m, reduction = 'tsne', group.by = 'age_class',label=1)

corAging <- function(x,agingScore){cor <- cor(x,agingScore);cor}
cor.genes <- apply(GetAssayData(Endo.m),1,corAging,agingScore=Endo.m$hSI)
cor.genes[is.na(cor.genes)] <- 0
features <- rownames(Endo.m)[order(abs(cor.genes),decreasing=TRUE)[1:1500]]
features

Endo.matrix <- as.matrix(GetAssayData(Endo.m))[features,]
Endo.m[['AgingExp']] <- CreateAssayObject(Endo.matrix)
DefaultAssay(Endo.m) <- "AgingExp"

Endo.m <- ScaleData(Endo.m,features=rownames(Endo.m))
Endo.m <- RunPCA(Endo.m,npcs = 50,verbose = F,features=rownames(Endo.m)) 
ElbowPlot(Endo.m,ndims=50,reduction="pca")

Endo.m$age_state = factor(ifelse(Endo.m$age_class == 1,"Normal",ifelse(Endo.m$age_class == 2,"Transition","Senescent")),levels = c("Normal","Transition","Senescent"))
DimPlot(Endo.m, reduction = 'pca', group.by = 'age_state',label=F,pt.size = 1.5)

### expression of SASP
DefaultAssay(Endo.m) <- 'RNA'
Endo.m <- ScaleData(Endo.m,features = rownames(Endo.m))
sasp = read.csv("/mnt/data1/wangj/AgingScore/AgingScorePro/SASP.csv")[[1]]
sasp = sasp[sasp %in% rownames(Endo.m)]

markers = FindMarkers(Endo.m,
                      ident.1 = 'Senescent',
                      ident.2 = 'Normal',
                        assay = 'RNA',
                        features = rownames(Endo.m),
                        only.pos = T,
                        group.by = 'age_state',
                        logfc.threshold = 0.1,
                        min.diff.pct = 0)
sasp_plot = rownames(markers)[rownames(markers) %in% sasp]

### age state marker enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
### remove ribosomal protein
ribosomal = read.table("/mnt/data1/wangj/MouseBrain/Ribosome.txt",stringsAsFactors=FALSE)
res = markers[!rownames(markers) %in% ribosomal$V1,]
# GO
pathways = read.gmt("/mnt/data1/wangj/GeneSets/c5.go.v2023.1.Hs.symbols.gmt")
go <- enricher(gene = rownames(res), universe = rownames(Endo.m@assays$RNA),TERM2GENE = pathways,pvalueCutoff=0.05)

# KEGG
pathways = read.gmt("/mnt/data1/wangj/GeneSets/KEGG.gmt")
kegg <- enricher(gene = rownames(res), universe = rownames(Endo.m@assays$RNA),TERM2GENE = pathways)

##### validate in microarray
marker_set = list()
for(i in unique(Endo.m$age_state)){
  marker_set[[i]] <- FindMarkers(Endo.m,ident.1 = i,group.by="age_state",only.pos = TRUE,assay = "RNA",logfc.threshold = 0.1,min.diff.pct = 0)
  marker_set[[i]] <- filter(marker_set[[i]],p_val_adj<0.05)
}
lapply(marker_set,dim)
### bulk degs
library(annaffy)
load('/home/wangjing/wangj/AgingScore/Data/Bulk_Microarray/Valid.RData')

GPL = c("GPL570" = "hgu133plus2.db","GPL3921" = "hthgu133a.db","GPL11532" = "hugene11sttranscriptcluster.db")
gene_id = aafSymbol( rownames(ArrayList$GSE77239[[1]]), GPL[ArrayList$GSE77239[[1]]@annotation]) %>% as.character
exp = ArrayList$GSE77239[[1]] %>% exprs %>% set_rownames( gene_id )
exp = exp[rownames(exp)!="character(0)",]
head(exp)

meta <- ArrayList [["GSE77239"]][[1]] %>% 
        pData %>% 
        mutate(condition= `cells:ch1`) %>% 
        mutate(condition = factor( ifelse(grepl('young',condition),'young','Senescent'),levels = c("young","Senescent"),ordered = T)) 
meta[,'condition']

geneID <- intersect(rownames(exp),rownames(Endo.m@assays$RNA))
exp_bulk <- exp[geneID,rownames(meta)]

deg.tab <- NULL
pheno <- ifelse(meta[,'condition'] == 'young',0,1)
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
    for (state in levels(Endo.m$age_state)) {
      genes = rownames(marker_set[[state]])
      n = length(genes)
      k = length(intersect(bk,genes))
      enrich_reList[[paste(state,de,sep = '_')]] = phyper(k-1,M, N-M, n, lower.tail=FALSE)
      enrich_geneList[[paste(state,de,sep = '_')]][['bulk']] = bk
      enrich_geneList[[paste(state,de,sep = '_')]][['single-cell']] = genes
  }
}

### correlation between senlung.objnct percentage and DTD/age
ClinicMeta = read.csv('/home/wangjing/wangj/AgingScore/GSE163530_COVID-19/GSE162911_GeoMx/ClinicMeta.csv',row.names = 1)
### calculate Senescent percentage of each donor
result <- Endo.m@meta.data %>%
            mutate(Donor = gsub('_[123]','',donor))
df = aggregate(result$age_state, by = list(result$Donor), FUN = table)
df = cbind(df['Group.1'],df['x']/rowSums(df['x']))
colnames(df) = c('Donor','Normal','Transition','Senescent')
df = column_to_rownames(df,'Donor')
df$DTD = ClinicMeta[rownames(df),]$Days_to_death
df$Age = ClinicMeta[rownames(df),]$Age
df$Age = gsub('>89','90-100',df$Age)
df$Age_min = strsplit(df$Age,'-') %>% sapply(function(x) {as.numeric(x[1])})
df$Age_max = strsplit(df$Age,'-') %>% sapply(function(x) {as.numeric(x[2])})
df$Age_mean = (df$Age_min+df$Age_max)/2
df$group = ifelse(df$DTD>16.5,ifelse(df$Age_mean < 70,'middle_moderate','old_moderate'),ifelse(df$Age_mean < 70,'middle_severe','old_severe'))
table(df$group)

df_plot = cbind(reshape2::melt(df[,c("Normal","Transition","Senescent")]),Group = rep(df$group,3))
df_plot = aggregate(df_plot$value, by = list(df_plot$Group,df_plot$variable), FUN = median)
colnames(df_plot) = c('Group','age_state','value')
df_plot$Group = factor(df_plot$Group,levels = c('middle_moderate','middle_severe','old_moderate','old_severe'))
df_plot

Endo.m$celltype_state <- paste(as.character(Endo.m$age_state),as.character(Endo.m$celltype),sep = ":")
Endo.m$group <- df[gsub('_[123]','',Endo.m$donor),]$group


###### deconvolution on spatial
### calculate major cell type profile
lung.obj$subtype <- as.character(lung.obj$Cluster)
lung.obj$subtype[names(Endo.m$age_state)] <- as.character(Endo.m$age_state)
table(lung.obj$subtype)
profile = AverageExpression(lung.obj,assays = 'RNA',group.by = 'subtype',slot = 'data',)$RNA %>% as.matrix()
head(profile)

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
### deconvolution
wta <- wta[rowSums(wta)>0,]

features = intersect(rownames(wta),rownames(profile))
Bulk <- log2(wta[features,]+1)
profile <- profile[features,]

library(TED)
ted <- run.Ted(ref.dat= t(profile),
                X=t(Bulk),
                cell.type.labels= colnames(profile),
                input.type="GEP",
                n.cores=10)
Frac <- ted$res$final.gibbs.theta
colMeans(Frac)

Hmisc::rcorr(Frac,type='spearman')



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













###### cell chat analysis
library(CellChat)
library(ggalluvial)

lung.obj$subtype = lung.obj$celltype
lung.obj$subtype[colnames(EpiExp.m)] = as.character(EpiExp.m$age_state)
lung.obj$subtype = gsub('\\.','',lung.obj$subtype) %>% factor
table(lung.obj$subtype)

cellchat <- createCellChat(object = lung.obj,group.by = 'subtype')
cellchat@DB <- CellChatDB.human 
table(cellchat@idents)

future::plan("multicore", workers = 10) 
options(future.globals.maxSize = 8000 * 1024^2)

cellchat <- subsetData(cellchat,rownames(lung.obj))
cellchat <- identifyOverExpressedGenes(cellchat) 

cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = T, population.size=TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

