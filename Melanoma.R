library(dplyr)
library(annaffy)
library(magrittr)
library(plyr)
library(mclust)
library(ICAnet)
library(phateR)
library(reticulate)
library(Seurat)

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
DimPlot(melanoma_obj, reduction = 'tsne', group.by = 'celltype',label=1)

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
DimPlot(EpiExp.m, reduction = 'tsne', group.by = 'age_class',label=1)

###Running the ICAnet
source('/home/wangjing/wangj/AgingScore/AgingScorePro/Data1_Scripts/getPPI_String.R')

EpiExp.matrix <- as.matrix(GetAssayData(EpiExp.m,))[features,]
EpiExp.m[['AgingExp']] <- CreateAssayObject(EpiExp.matrix)
DefaultAssay(EpiExp.m) <- "AgingExp"

EpiExp.m$batch <- rep("batch",ncol(EpiExp.m))
Ica.epi <- ICAcomputing(EpiExp.m,ICA.type="JADE",RMT=TRUE,two.stage=FALSE)
Ica.filter <- CrossBatchGrouping(Ica.epi$ica.pooling)
# PPI <- readRDS("PPI_feature.rds")
PPI <- getPPI_String(EpiExp.m,species=9606)

EpiExp.m <- RunICAnet(EpiExp.m,Ica.epi$ica.pooling,PPI.net = PPI,scale=FALSE,ModuleSignificance = FALSE,cores = 1,aucMaxRank=900)
dev.off()

### Using PHATE to predict the transition process of tumor aging
# use_python("/home/tools/anaconda3/envs/scpy-env/bin/python3", required = T)

Epi.data <- t((GetAssayData(EpiExp.m)))
Epi.data <- as.data.frame(Epi.data)
Epi.phate <- phate(Epi.data)
branch <- EpiExp.m$age_class
palette(rainbow(3))
plot(Epi.phate, col = branch,pch=16)

EpiExp.m[['phate']] <- CreateDimReducObject(Epi.phate$embedding,key="phate_")

EpiExp.m$age_state = factor(ifelse(EpiExp.m$age_class == 1,"Cycling",ifelse(EpiExp.m$age_class == 2,"Moderate senescent","Senescent")),levels = c("Cycling","Moderate senescent","Senescent"))

### aging marker genes
SenMarkers = c("CDKN1A", "SERPINE1")

### valida in microarray
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

marker_set = list()
for(i in unique(EpiExp.m$age_state)){
  marker_set[[i]] <- FindMarkers(EpiExp.m,ident.1 = i,group.by="age_state",only.pos = TRUE,assay = "RNA")
  marker_set[[i]] <- filter(marker_set[[i]],p_val_adj<0.05)
}
lapply(marker_set,dim)

gene_set <- lapply(marker_set,function(x){x <- rownames(x)[rownames(x)%in%rownames(exp)]})

mat <- exp[unlist(gene_set),rownames(meta)]
mat=mat[!apply(mat, 1, sd)==0,]
mat=Matrix::t(scale(Matrix::t(mat),center=TRUE))
mat=mat[is.na(row.names(mat)) == FALSE,]
mat[is.nan(mat)] = 0
mat[mat>2] = 2
mat[mat< -2] = -2