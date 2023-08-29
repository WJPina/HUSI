library(dplyr)
library(annaffy)
library(magrittr)
library(plyr)
library(mclust)
library(ICAnet)
library(phateR)
library(reticulate)
library(Seurat)
library(tibble)
library(data.table)

setwd('/home/wangjing/wangj/AgingScore/Data/scRNA_melanome')

use_python("/home/tools/anaconda3/envs/scpy-env/bin/python3", required = T)
mm_l2 = readRDS("/home/wangjing/wangj/AgingScore/Data/Bulk_TrainModel/mm_l2.rds")
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
melanoma_obj$celltype <- celltype

### calculate aging score
AgeScore  = melanoma_obj@assays$RNA@data[] %>% {apply( ., 2, function(z) {cor(z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )})}
melanoma_obj$hSI <- AgeScore[colnames(melanoma_obj)]

melanoma_obj <- FindVariableFeatures(object = melanoma_obj, selection.method = 'vst', nfeatures =1500)
melanoma_obj <- ScaleData(melanoma_obj,features=rownames(melanoma_obj))
melanoma_obj <- RunPCA(melanoma_obj,npcs = 30,verbose = F,features=rownames(melanoma_obj)) 
ElbowPlot(melanoma_obj,ndims=30,reduction="pca")
melanoma_obj <- RunTSNE(melanoma_obj, reduction = "pca", dims = 1:20) 
DimPlot(melanoma_obj, reduction = 'pca', group.by = 'celltype',label=1)

### extract tumor cells
EpiExp.m <- melanoma_obj[,melanoma_obj$celltype=="Tumor"]
### divided tumor cells by gaussion 
gaussian = EpiExp.m@meta.data$hSI %>% {log2((1+ .)/(1- .))} %>% Mclust()
EpiExp.m$age_class <- gaussian$classification
DimPlot(EpiExp.m, reduction = 'pca', group.by = 'age_class',label=1)
### ranking the genes based on each gene expression correlation with the aging score
corAging <- function(x,agingScore){cor <- cor(x,agingScore);cor}
cor.genes <- apply(GetAssayData(EpiExp.m),1,corAging,agingScore=EpiExp.m$hSI)
cor.genes[is.na(cor.genes)] <- 0
features <- rownames(EpiExp.m)[order(abs(cor.genes),decreasing=TRUE)[1:1500]]
features 
### Running the ICAnet
source('/home/wangjing/wangj/AgingScore/AgingScorePro/Data1_Scripts/getPPI_String.R')

EpiExp.matrix <- as.matrix(GetAssayData(EpiExp.m))[features,]
EpiExp.m[['AgingExp']] <- CreateAssayObject(EpiExp.matrix)
DefaultAssay(EpiExp.m) <- "AgingExp"

EpiExp.m$batch <- rep("batch",ncol(EpiExp.m))
Ica.epi <- ICAcomputing(EpiExp.m,ICA.type="JADE",RMT=TRUE,two.stage=FALSE)
Ica.filter <- CrossBatchGrouping(Ica.epi$ica.pooling)
dev.off()
PPI <- readRDS("PPI_melanoma.RDS")
# PPI <- getPPI_String(EpiExp.m,species=9606)
EpiExp.m <- RunICAnet(EpiExp.m,Ica.epi$ica.pooling,PPI.net = PPI,scale=FALSE,ModuleSignificance = FALSE,cores = 1,aucMaxRank=900)

### Using PHATE to predict the transition process of tumor aging
Epi.data <- t((GetAssayData(EpiExp.m)))
Epi.data <- as.data.frame(Epi.data)
Epi.phate <- phate(Epi.data)
branch <- EpiExp.m$age_class
palette(rainbow(3))
plot(Epi.phate, col = branch,pch=16)

EpiExp.m[['phate']] <- CreateDimReducObject(Epi.phate$embedding,key="phate_")
EpiExp.m$age_state = factor(ifelse(EpiExp.m$age_class == 1,"Cycling",ifelse(EpiExp.m$age_class == 2,"Transition","Senescent")),levels = c("Cycling","Transition","Senescent"))
save(EpiExp.m,file = paste("tumor_",Sys.Date(),'.RData', sep = ""))

### age state degs
marker_set = list()
for(i in unique(EpiExp.m$age_state)){
  marker_set[[i]] <- FindMarkers(EpiExp.m,ident.1 = i,group.by="age_state",only.pos = TRUE,assay = "RNA")
  marker_set[[i]] <- filter(marker_set[[i]],p_val_adj<0.05)
}
lapply(marker_set,dim)
###### enrichment of each state marker 
library(clusterProfiler)
library(org.Hs.eg.db)
### remove ribosomal protein
ribosomal = read.table("/mnt/data1/wangj/MouseBrain/Ribosome.txt",stringsAsFactors=FALSE)
marker_set_mr = lapply(marker_set,function(x) {x <- x[!rownames(x) %in% ribosomal$V1,] %>% rownames_to_column('gene')})
res = rbindlist(marker_set_mr,idcol = 'state')
ids=bitr(res$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
sce.markers=merge(res,ids,by.x='gene',by.y='SYMBOL')
gcSample=split(sce.markers$ENTREZID, sce.markers$state)
bg=bitr(rownames(EpiExp.m@assays$RNA),'SYMBOL','ENTREZID','org.Hs.eg.db')
# GO
go <- compareCluster(gcSample, fun= "enrichGO",ont = "BP",OrgDb= "org.Hs.eg.db", pvalueCutoff=0.05,pAdjustMethod = "fdr",universe = bg$ENTREZID)

# KEGG
kegg <- compareCluster(gcSample, fun= "enrichKEGG",organism="hsa", pvalueCutoff=0.05,pAdjustMethod = "fdr",universe = bg$ENTREZID)

##### validate in microarray
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

geneID <- intersect(rownames(exp),rownames(EpiExp.m@assays$RNA))
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

##### analysis on TCGA SKCM
### load TCGA SKCM bulk data
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

### only keep the degs of each single cell state
# features = unlist(lapply(marker_set, function(x) {rownames(x[order(x$avg_log2FC,decreasing = T),])})) %>% unique %>% intersect(rownames(tcga_melanoma),.)
features = intersect(rownames(tcga_melanoma),rownames(EpiExp.m@assays$RNA))
Bulk <- tcga_melanoma[features,]
Bulk <- Bulk[rowSums(Bulk)>0,]
dim(Bulk)
profile = AverageExpression(EpiExp.m, features = rownames(Bulk), slot = 'count',assays = 'RNA',group.by = 'age_state')$RNA %>% as.matrix()
head(profile)
### deconvolute TCGA SKCM sample by ENIGMA
source("/home/wangjing/wangj/ENIGMA/ENIGMASpatialPro/scripts/ENIGMA.R")
source("/home/wangjing/wangj/ENIGMA/ENIGMASpatialPro/scripts/ENIGMA_revise.R")
Frac <- get_proportion(Bulk, profile)
Frac <- Frac$theta
colMeans(Frac)

### using CIBERSORT to deconvolute TCGA SKCM sample
# source("/home/wangjing/wangj/codebase/HUSI/CIBERSORT.R")
# out = tcga_melanoma
# out = rbind(ID = colnames(out),out)
# out[1:5,1:5]
# write.table(out,file="tcga_melanoma.txt",sep="\t",quote=F,col.names=F)  

# SKCM_CIBER = CIBERSORT("LM22.txt", "tcga_melanoma.txt", perm=100, QN=TRUE)
# write.csv(SKCM_CIBER,"tcga_melanoma_CIBERSORT.csv")

SKCM_CIBER = read.csv('tcga_melanoma_CIBERSORT.csv',row.names = 1)
head(SKCM_CIBER)
### calculate the correlation between each single cell state and each immune cell type
corList = list()
for(i in 1:ncol(Frac)){
  cors = c()
  for(j in 1:(ncol(SKCM_CIBER)-3)){
    cors = c(cors,cor(Frac[rownames(SKCM_CIBER),i],SKCM_CIBER[,j],method='spearman'))
  }
  corList[[colnames(Frac)[i]]] = cors
}

cor_mat = as.data.frame(corList)
rownames(cor_mat) = colnames(SKCM_CIBER)[1:(ncol(SKCM_CIBER)-3)]
cor_mat = t(cor_mat)

### survival analysis
library(survival)
library(survminer)

survival_data <- clinicalMatrix[, c("days_to_last_followup", "vital_status")]
### remove samples with NA in survival data
survival_data <- survival_data[!survival_data$days_to_last_followup%in%c("","[Discrepancy]"),]
### relpace LIVING with 0 and DECEASED with 1
survival_data$vital_status <- as.numeric(survival_data$vital_status == "DECEASED")
survival_data$days_to_last_followup <- as.numeric(survival_data$days_to_last_followup)

data = data.frame(Frac[rownames(survival_data),],OS = survival_data$vital_status,OS.time = survival_data$days_to_last_followup)

###### cell cat analysis
library(CellChat)
library(ggalluvial)

melanoma_obj$subtype = melanoma_obj$celltype
melanoma_obj$subtype[colnames(EpiExp.m)] = as.character(EpiExp.m$age_state)
melanoma_obj$subtype = gsub('\\.','',melanoma_obj$subtype) %>% factor
table(melanoma_obj$subtype)

cellchat <- createCellChat(object = melanoma_obj,group.by = 'subtype')
cellchat@DB <- CellChatDB.human 
table(cellchat@idents)

future::plan("multicore", workers = 10) 
options(future.globals.maxSize = 8000 * 1024^2)

cellchat <- subsetData(cellchat,rownames(melanoma_obj))
cellchat <- identifyOverExpressedGenes(cellchat) 

cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = T, population.size=TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

### pathway compare
sig = cellchat@LR$LRsig
pathways = names(cellchat@netP$centr)
df = data.frame(
  pathway = pathways,
  sene_out = unlist(lapply(cellchat@netP$centr,function(x){x[['outdeg']]['Senescent']})) ,
  sene_income = unlist(lapply(cellchat@netP$centr,function(x){x[['indeg']]['Senescent']})) ,
  cycle_out = unlist(lapply(cellchat@netP$centr,function(x){x[['outdeg']]['Cycling']})) ,
  cycle_income = unlist(lapply(cellchat@netP$centr,function(x){x[['indeg']]['Cycling']})) ,
  row.names = pathways)

df = filter(df,rowSums(df[,-1])!=0)
pathways = df[df$sene_income!=0&df$cycle_income == 0,]$pathway
pathways = c("CSPG4","CD6","BMP","CCL","TGFb" )
cellchat@idents <- factor(cellchat@idents,levels = c('Endo cell','B cell','T cell','NK cell','Macro cell','CAF cell',"Cycling","Transition","Senescent"))



netVisual_bubble(cellchat, sources.use = 'Senescent', targets.use = 'T cell', remove.isolate = FALSE)

nPatterns = 6
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
pattern <- cellchat@netP$pattern$outgoing$pattern$signaling
[cellchat@netP$pattern$outgoing$pattern$signaling$Contribution>0.5,]
pattern_use = pattern[pattern$Pattern %in% c("Pattern 1","Pattern 6"),]
netAnalysis_river(cellchat, pattern = "outgoing",cutoff = 0.5,sources.use = c('Cycling','Transition','Senescent'),signaling = as.character(pattern_use$Signaling))

### survival analysis of receptor-ligand pairs
data = data.frame(t(tcga_melanoma[c('BMPR1B','BMPR2','TGFBR1','TGFBR2'),rownames(survival_data)]),OS = survival_data$vital_status,OS.time = survival_data$days_to_last_followup)

