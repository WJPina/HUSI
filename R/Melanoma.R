##### This scripts is used to analyze the single cell melanoma data 
################################# library ######################################
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
library(destiny)

source('R/functions.R')
source('R/functions_new.R')

setwd('hUSI/')
############################# load melanoma data ###############################
melanoma <- read.table("GSE72056_melanoma_single_cell_revised_v2.txt",sep="\t",header=TRUE) ### download from GEO
melanoma[1:5,1:5]

sample = as.numeric(melanoma[1,-1])
sample = paste('Mel',sample,sep='')
malignant = as.numeric(melanoma[2,-1])
malignant <- ifelse(malignant==0,'unresolved',ifelse(malignant==1,'non-malignant','malignant'))
celltype = as.numeric(melanoma[3,-1])
celltype <- stringi::stri_replace_all_regex(celltype ,c(0:6),c('Tumor','T','B','Macro','Endo','CAF','NK'),vectorize=F)

melanoma_meta = data.frame(cell = colnames(melanoma)[-1],sample,malignant,celltype)
melanoma_meta$treatment = ifelse(melanoma_meta$sample %in% c('Mel53', 'Mel59', 'Mel65', 'Mel67', 'Mel71', 'Mel79', 'Mel80', 'Mel81', 'Mel82', 'Mel84', 'Mel89'),'untreated','treated')

### filter unresolved cells
melanoma_meta = melanoma_meta[melanoma_meta$malignant!='unresolved',]
### filter tumor cells with non-malignant state or normal cells with malignant state
melanoma_meta = melanoma_meta[(melanoma_meta$malignant=='malignant'&melanoma_meta$celltype =='Tumor')|((melanoma_meta$malignant=='non-malignant'&melanoma_meta$celltype !='Tumor')),]
rownames(melanoma_meta) = melanoma_meta$cell

### filter cells and genes
melanoma_filter <- melanoma[-c(1:3),rownames(melanoma_meta)]
genes = melanoma[-c(1:3),1]
melanoma_filter <- melanoma_filter[!duplicated(genes),]
rownames(melanoma_filter) <- genes[!duplicated(genes)]
dim(melanoma_filter)
melanoma_filter[1:5,1:5]

### build object
melanoma_obj <- CreateSeuratObject(melanoma_filter,meta.data = melanoma_meta)
melanoma_obj <- FindVariableFeatures(object = melanoma_obj, selection.method = 'vst', nfeatures =1500)
melanoma_obj <- ScaleData(melanoma_obj,features=rownames(melanoma_obj))
melanoma_obj <- RunPCA(melanoma_obj,npcs = 30,verbose = F,features=rownames(melanoma_obj)) 
ElbowPlot(melanoma_obj,ndims=30,reduction="pca")
melanoma_obj <- RunTSNE(melanoma_obj, reduction = "pca", dims = 1:20) 
DimPlot(melanoma_obj,group.by = 'celltype',label = T)
########################### calculate aging score ##############################
load('Data/SenOCLR_l2=1_drop.rdata')

hUSI  = scoreOCLR(GetAssayData(melanoma_obj),SenOCLR,'spearman')[[1]] %>% minmax()
melanoma_obj$hUSI = hUSI

################################ only tumor cells ##############################
EpiExp.m <- melanoma_obj[,melanoma_obj$celltype=="Tumor"]

### divided tumor cells by gaussion
set.seed(123)

data = EpiExp.m$hUSI
model = data %>% Mclust()
meanV = model$parameters$mean

SenClass = ifelse(data<meanV[1],1,ifelse(data<meanV[2],2,3))

EpiExp.m$SenClass <- factor(paste0("C",SenClass),levels = c('C1','C2','C3'))
table(EpiExp.m$SenClass)
################################### ICAnet #####################################
### correlation between gene expression and hUSI
corAging <- function(x,agingScore){cor <- cor(x,agingScore);cor}
cor.genes <- apply(GetAssayData(EpiExp.m),1,corAging,agingScore=EpiExp.m$hUSI)
cor.genes[is.na(cor.genes)] <- 0
features <- names(cor.genes)[order(abs(cor.genes),decreasing=TRUE)][1:2000]
features 

### Running the ICAnet
EpiExp.matrix <- as.matrix(GetAssayData(EpiExp.m))[features,]
EpiExp.m[['SenExp']] <- CreateAssayObject(EpiExp.matrix)
DefaultAssay(EpiExp.m) <- "SenExp"

EpiExp.m$batch <- rep("batch",ncol(EpiExp.m))
Ica.epi <- ICAcomputing(EpiExp.m,ICA.type="JADE",RMT=TRUE,two.stage=FALSE)
rownames(Ica.epi$ica.pooling) <- features 
# Ica.filter <- CrossBatchGrouping(Ica.epi$ica.pooling)
# dev.off()

PPI <- getPPI_String(EpiExp.m,species=9606)

EpiExp.m <- RunICAnet(EpiExp.m,Ica.epi$ica.pooling,PPI.net = PPI,scale=FALSE,ModuleSignificance = FALSE,cores = 1,aucMaxRank=500)
# dev.off()

### Using PHATE to predict the transition process of tumor aging
library(phateR)

EpiExp.m <- RunPCA(EpiExp.m,features = rownames(EpiExp.m),reduction.key = "pcaNet_",reduction.name = "pcaNet",npcs = 5) 
Epi.data <- EpiExp.m@reductions$pcaNet@cell.embeddings
Epi.data <- as.data.frame(Epi.data)
Epi.phate <- phate(Epi.data)
EpiExp.m[['phate']] <- CreateDimReducObject(Epi.phate$embedding,key="phate_")

DimPlot(EpiExp.m, reduction = 'phate', group.by = 'SenClass',label=T)

EpiExp.m$State = factor(ifelse(EpiExp.m$SenClass == 'C1',"Cycling",ifelse(EpiExp.m$SenClass == 'C2',"Transitional","Senescent")),levels = c("Cycling","Transitional","Senescent"))

# save(EpiExp.m,file='EpiExp.m.Rdata')
load("EpiExp.m.Rdata")
########################### Senescence state degs ##############################
### state marker genes
marker_set = list()
marker_set[['Cycling']] <- FindMarkers(EpiExp.m,ident.1 = 'Cycling',ident.2 = 'Senescent',group.by="State",only.pos = TRUE,assay = "RNA")
marker_set[['Transitional']] <- FindMarkers(EpiExp.m,ident.1 = 'Transitional',ident.2 = 'Cycling',group.by="State",only.pos = TRUE,assay = "RNA")
marker_set[['Senescent']] <- FindMarkers(EpiExp.m,ident.1 = 'Senescent',ident.2 = 'Cycling',group.by="State",only.pos = TRUE,assay = "RNA")
marker_set = lapply(marker_set,function(x){x<-filter(x,p_val_adj<0.05)})
lapply(marker_set,dim)

### bulk  microarray degs
GPL = c("GPL570" = "hgu133plus2.db","GPL3921" = "hthgu133a.db","GPL11532" = "hugene11sttranscriptcluster.db")
Array <- GEOquery::getGEO('GSE83922', destdir="./Data/",AnnotGPL = F,getGPL = F)[[1]]
gene_id = aafSymbol(rownames(Array), GPL[Array@annotation]) %>% as.character
exp = Array %>% exprs %>% set_rownames( gene_id )
exp = exp[rownames(exp)!="character(0)",]
head(exp)

meta <- Array %>% 
        pData %>% 
        mutate(condition= `cell phenotype:ch1`) %>% 
        .[which(.$condition %in% c("young","senescent")),] %>% 
        mutate(condition = factor(condition,levels = c("young","senescent"),ordered = T)) 
meta[,'cell phenotype:ch1']

geneID <- intersect(rownames(exp),rownames(EpiExp.m@assays$RNA))
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
  for (state in levels(EpiExp.m$State)) {
    genes = rownames(marker_set[[state]])
    n = length(genes)
    k = length(intersect(bk,genes))
    enrich_reList[[paste(state,de,sep = '_')]] = phyper(k-1,M, N-M, n, lower.tail=FALSE)
    enrich_geneList[[paste(state,de,sep = '_')]][['bulk']] = bk
    enrich_geneList[[paste(state,de,sep = '_')]][['single-cell']] = genes
  }
}
####################### enrichment of each state marker ########################
library(clusterProfiler)
library(org.Hs.eg.db)
### remove ribosomal protein
ribosomal = read.table("Data/Ribosome.txt",stringsAsFactors=FALSE)
marker_set_mr = lapply(marker_set,function(x) {x <- x[!rownames(x) %in% ribosomal$V1,] %>% rownames_to_column('gene')})
res = rbindlist(marker_set_mr,idcol = 'State')
res = filter(res,State%in%c("Cycling","Transitional","Senescent") & avg_log2FC>0.5)
res$State = factor(res$State,levels = c("Senescent","Transitional","Cycling"),ordered = T)

### GOBP
ids=bitr(res$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
sce.markers=merge(res,ids,by.x='gene',by.y='SYMBOL')
gcSample=split(sce.markers$ENTREZID, sce.markers$State)
bg=bitr(rownames(EpiExp.m@assays$RNA),'SYMBOL','ENTREZID','org.Hs.eg.db')
go <- compareCluster(gcSample, fun= "enrichGO",ont = "BP",OrgDb= "org.Hs.eg.db",
                     pvalueCutoff=0.05,pAdjustMethod = "fdr",universe = bg$ENTREZID)
############################### Analysis on TCGA SKCM ##########################
### load TCGA SKCM bulk data RPKM
clinicalMatrix <- read.csv("Data/SKCM_clinicalMatrix.csv",header=TRUE,stringsAsFactors = FALSE,row.names=1)
tcga_melanoma <- read.table("Data/melanoma_tcga",header=TRUE,row.names=1,stringsAsFactors=FALSE)

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
res = rbindlist(marker_set_mr,idcol = 'State')
res = filter(res,State%in%c("Cycling","Transitional","Senescent") & avg_log2FC>0)

features = intersect(rownames(tcga_melanoma),rownames(EpiExp.m@assays$RNA))
Bulk <- tcga_melanoma[features,]
Bulk <- Bulk[rowSums(Bulk)>0,]
dim(Bulk)

profile = AverageExpression(EpiExp.m, features = rownames(Bulk), slot = 'counts',assays = 'RNA',group.by = 'State')$RNA %>% as.matrix()
head(profile)

profile = profile[rownames(profile)[rownames(profile) %in% res$gene],]
dim(profile)

### deconvolute TCGA SKCM sample 
library(EpiDISH)
Frac_epi <- epidish(beta.m = Bulk, ref.m = profile, method = "RPC") 
Frac <- Frac_epi$estF
colMeans(Frac)

### using CIBERSORT to deconvolute TCGA SKCM sample
source("R/CIBERSORT.R")
out = tcga_melanoma
out = rbind(ID = colnames(out),out)
out[1:5,1:5]

SKCM_CIBER = CIBERSORT("Data/LM22.txt", "Data/tcga_melanoma.txt", perm=100, QN=TRUE)
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
cor_mat

### survival analysis
library(survival)
library(survminer)

survival_data <- clinicalMatrix[, c("days_to_last_followup", "days_to_death", "vital_status")]
survival_data <- survival_data[!survival_data$days_to_last_followup%in%c("","[Discrepancy]"),]
survival_data$vital_status <- as.numeric(survival_data$vital_status == "DECEASED")
survival_data[is.na(survival_data)] <- 0
survival_data$os.time = as.numeric(survival_data$days_to_last_followup)
survival_data$os.time[survival_data$vital_status==1] <- as.numeric(survival_data$days_to_death[survival_data$vital_status==1]) 

data = data.frame(Frac[rownames(survival_data),],OS = survival_data$vital_status,OS.time = survival_data$os.time)

############################# cell chat analysis ############################### 
library(CellChat)
library(ggalluvial)

subtype = melanoma_obj$celltype
subtype[colnames(EpiExp.m)] <- as.character(EpiExp.m$State)
table(subtype)

melanoma_obj$subtype = factor(subtype,levels = c('Endo','B','T','NK','Macro','CAF',"Cycling","Transitional","Senescent"))
melanoma_obj$subtype[colnames(EpiExp.m)] = as.character(EpiExp.m$State)
melanoma_obj$subtype = gsub('\\.','',melanoma_obj$subtype) %>% factor
table(melanoma_obj$subtype)

cellchat <- createCellChat(object = melanoma_obj,group.by = 'subtype')
cellchat@DB <- CellChatDB.human 
table(cellchat@idents)
cellchat@idents <- factor(cellchat@idents,levels = c('Endo','B','T','NK','Macro','CAF',"Cycling","Transitional","Senescent"))

future::plan("multicore", workers = 20) 
options(future.globals.maxSize = 8000 * 1024^2)

cellchat <- subsetData(cellchat,rownames(melanoma_obj))
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat)
# cellchat <- projectData(cellchat, PPI.human)
cellchat <- smoothData(cellchat, adj =PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = T, population.size=TRUE)#!!!
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

### unique receive pathway in senescent
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
pathways

### survival analysis of receptor-ligand pairs
data = data.frame(t(tcga_melanoma[c('TGFBR2','ACVR1','ALCAM','CCR10'),rownames(survival_data)]),OS = survival_data$vital_status,OS.time = survival_data$os.time)

### find output patterns
nPatterns = 6
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
pattern <- cellchat@netP$pattern$outgoing$pattern$signaling[cellchat@netP$pattern$outgoing$pattern$signaling$Contribution>0.5,]
pattern_use = pattern[pattern$Pattern %in% c("Pattern 1","Pattern 6"),]
pattern_use$Signaling
