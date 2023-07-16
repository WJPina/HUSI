library(dplyr)
library(annaffy)
library(magrittr)
library(plyr)
library(mclust)
library(ICAnet)
library(phateR)
library(reticulate)
library(Seurat)
use_python("/home/tools/anaconda3/envs/scpy-env/bin/python3", required = T)


mm_l2 = readRDS("/home/wangjing/wangj/AgingScore/Data/Bulk_TrainModel/mm_l2.rds")
setwd('/home/wangjing/wangj/AgingScore/Data/scRNA_melanome')
############################### application on melanoma ########################
### load melanoma data
melanoma <- read.table("GSE72056_melanoma_single_cell_revised_v2.txt",sep="\t",header=TRUE)
### load meta data
sample <- read.csv("sample.csv",row.names=1,stringsAsFactors=FALSE)
metadata <- t(melanoma[c(2,3),-1])
metadata <- metadata[metadata[,1]!=0,]
metadata.sub <- metadata[metadata[,2]==0,]
metadata.sub <- metadata.sub[metadata.sub[,1]==2,]
metadata.sub2 <- metadata[metadata[,2]!=0,]
metadata.sub2 <- metadata.sub2[metadata.sub2[,1]!=2,]
metadata <- rbind(metadata.sub,metadata.sub2)
celltype <- c(rep("Tumor",nrow(metadata.sub)),metadata.sub2[,2])
celltype[celltype=="1"] <- "T cell"
celltype[celltype=="2"] <- "B cell"
celltype[celltype=="3"] <- "Macro. cell"
celltype[celltype=="4"] <- "Endo. cell"
celltype[celltype=="5"] <- "CAF cell"
celltype[celltype=="6"] <- "NK cell"

names(celltype) <- c(rownames(metadata.sub),rownames(metadata.sub2))
melanoma <- melanoma[!duplicated(melanoma[,1]),]
label <- melanoma[,1]
melanoma <- melanoma[-c(1:3),c(rownames(metadata.sub),rownames(metadata.sub2))]
rownames(melanoma) <- label[-c(1:3)]
library(Seurat)
melanoma_obj <- CreateSeuratObject(melanoma)

### calculate aging score
AgeScore  = melanoma_obj@assays$RNA@data[] %>% {apply( ., 2, function(z) {cor(z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )})}

melanoma_obj$hSI <- AgeScore[colnames(melanoma_obj)]
melanoma_obj$celltype <- celltype

melanoma_obj <- FindVariableFeatures(object = melanoma_obj, selection.method = 'vst', nfeatures =1500)
melanoma_obj <- ScaleData(melanoma_obj,features=rownames(melanoma_obj))
melanoma_obj <- RunPCA(melanoma_obj,npcs = 30,verbose = F,features=rownames(melanoma_obj)) 
# ElbowPlot(melanoma_obj,ndims=30,reduction="pca")
# melanoma_obj <- RunTSNE(melanoma_obj, reduction = "pca", dims = 1:20) 
DimPlot(melanoma_obj, reduction = 'pca', group.by = 'celltype',label=1)

### extract tumor cells
EpiExp.m <- melanoma_obj[,melanoma_obj$celltype=="Tumor"]

### ranking the genes based on each gene expression correlation with the aging score
corAging <- function(x,agingScore){cor <- cor(x,agingScore);cor}
cor.genes <- apply(GetAssayData(EpiExp.m),1,corAging,agingScore=EpiExp.m$hSI)
cor.genes[is.na(cor.genes)] <- 0
features <- rownames(EpiExp.m)[order(abs(cor.genes),decreasing=TRUE)[1:1500]]
features 

### divided tumor cells by gaussion 
gaussian = EpiExp.m@meta.data$hSI %>% {log2((1+ .)/(1- .))} %>% Mclust(G = 3)
EpiExp.m$age_class <- gaussian$classification
DimPlot(EpiExp.m, reduction = 'pca', group.by = 'age_class',label=1)

###Running the ICAnet
source('/home/wangjing/wangj/AgingScore/AgingScorePro/Data1_Scripts/getPPI_String.R')

EpiExp.matrix <- as.matrix(GetAssayData(EpiExp.m,))[features,]
EpiExp.m[['AgingExp']] <- CreateAssayObject(EpiExp.matrix)
DefaultAssay(EpiExp.m) <- "AgingExp"

EpiExp.m$batch <- rep("batch",ncol(EpiExp.m))
Ica.epi <- ICAcomputing(EpiExp.m,ICA.type="JADE",RMT=TRUE,two.stage=FALSE)
Ica.filter <- CrossBatchGrouping(Ica.epi$ica.pooling)
dev.off()
# PPI <- readRDS("PPI_feature.rds")
PPI <- getPPI_String(EpiExp.m,species=9606)
EpiExp.m <- RunICAnet(EpiExp.m,Ica.epi$ica.pooling,PPI.net = PPI,scale=FALSE,ModuleSignificance = FALSE,cores = 1,aucMaxRank=900)

### Using PHATE to predict the transition process of tumor aging
Epi.data <- t((GetAssayData(EpiExp.m)))
Epi.data <- as.data.frame(Epi.data)
Epi.phate <- phate(Epi.data)
branch <- EpiExp.m$age_class
palette(rainbow(3))
plot(Epi.phate, col = branch,pch=16)

EpiExp.m[['phate']] <- CreateDimReducObject(Epi.phate$embedding,key="phate_")
EpiExp.m$age_state = factor(ifelse(EpiExp.m$age_class == 1,"Cycling",ifelse(EpiExp.m$age_class == 2,"Moderate_senescent","Senescent")),levels = c("Cycling","Moderate_senescent","Senescent"))
save(EpiExp.m,file = paste("tumor_",Sys.Date(),'.RData', sep = ""))
### aging marker genes
DefaultAssay(EpiExp.m)='RNA'
SenMarkers = c("CDKN1A", "SERPINE1")

##### validate in microarray
### age state degs
marker_set = list()
for(i in unique(EpiExp.m$age_state)){
  marker_set[[i]] <- FindMarkers(EpiExp.m,ident.1 = i,group.by="age_state",only.pos = TRUE,assay = "RNA")
  marker_set[[i]] <- filter(marker_set[[i]],p_val_adj<0.05)
}
lapply(marker_set,dim)
### bulk degs
load('/home/wangjing/wangj/AgingScore/Data/Bulk_Microarray/Valid.RData')

GPL = c("GPL570" = "hgu133plus2.db","GPL3921" = "hthgu133a.db","GPL11532" = "hugene11sttranscriptcluster.db")
gene_id = aafSymbol( rownames(ArrayList$GSE83922[[1]]), GPL[ArrayList$GSE83922[[1]]@annotation]) %>% as.character
exp = ArrayList$GSE83922[[1]] %>% exprs %>% set_rownames( gene_id )
exp = exp[rownames(exp)!="character(0)",]
head(exp)

meta <- ArrayList [["GSE83922"]][[1]] %>% 
        pData %>% 
        mutate(condition= `cell phenotype:ch1`) %>% 
        .[which(.$condition %in% c("young","senescent")),] %>% 
        mutate(condition = factor(condition,levels = c("young","senescent"),ordered = T)) 
meta[,'cell phenotype:ch1']

geneID <- intersect(rownames(exp),rownames(EpiExp.m))
exp_bulk <- exp[geneID,rownames(meta)]

deg.tab <- NULL
pheno <- c(rep(0,4),rep(1,4))
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
    for (state in levels(EpiExp.m$age_state)) {
      genes = rownames(marker_set[[state]])
      n = length(genes)
      k = length(intersect(bk,genes))
      enrich_reList[[paste(state,de,sep = '_')]] = phyper(k-1,M, N-M, n, lower.tail=FALSE)
      enrich_geneList[[paste(state,de,sep = '_')]][['bulk']] = bk
      enrich_geneList[[paste(state,de,sep = '_')]][['single-cell']] = genes
  }
}

###### load TCGA SKCM bulk data
clinicalMatrix <- read.csv("SKCM_clinicalMatrix.csv",header=TRUE,stringsAsFactors = FALSE,row.names=1)
tcga_melanoma <- read.table("melanoma_tcga",header=TRUE,row.names=1,stringsAsFactors=FALSE)

sample_tcga <- rownames(clinicalMatrix)
sample_tcga <- gsub("-",".",sample_tcga)
rownames(clinicalMatrix) <- sample_tcga
sample_tcga <- intersect(sample_tcga,colnames(tcga_melanoma))
clinicalMatrix <- clinicalMatrix[sample_tcga,]

### filter samples by tumor index
tumor_index <- unlist(lapply(strsplit(rownames(clinicalMatrix),split="\\."),'[',4)) %>% as.numeric
clinicalMatrix <- clinicalMatrix[tumor_index<11,]
tcga_melanoma <- tcga_melanoma[,rownames(clinicalMatrix)]

# markerID <- sapply(marker_set,function(x){rownames(x)}) %>% unlist %>% as.character %>% unique
markerID = unlist(lapply(marker_set, function(x) {rownames(x[order(x$avg_log2FC,decreasing = T),][1:200,])})) %>% unique
tcga_melanoma <- tcga_melanoma[intersect(rownames(tcga_melanoma),markerID),]
tcga_melanoma <- tcga_melanoma[rowSums(tcga_melanoma)>0,]



### deconvolute TCGA SKCM sample by ENIGMA
source("/home/wangjing/wangj/ENIGMA/ENIGMASpatialPro/scripts/ENIGMA.R")
source("/home/wangjing/wangj/ENIGMA/ENIGMASpatialPro/scripts/ENIGMA_revise.R")
Bulk = tcga_melanoma
profile = AverageExpression(EpiExp.m, features = markerID , slot = 'data',assays = 'RNA',group.by = 'age_state')$RNA %>% data.frame
Frac <- get_proportion(Bulk, profile)
Frac <- Frac$theta
colMeans(Frac)

