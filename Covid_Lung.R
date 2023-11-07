library(reticulate)
library(dplyr)
library(magrittr)
library(Seurat)
library(mclust)
library(tibble)

setwd("/mnt/data1/wangj/AgingScore/GSE163530_COVID-19/GSE171668_scnRNA/")

mm_l2 = readRDS("/home/wangjing/wangj/AgingScore/Data/Bulk_TrainModel/mm_l2.rds")
use_python('/home/tools/anaconda3/envs/velocyto-env/bin/python')
sc <- import("scanpy")
###### preprocess data
### raw data
anndata = sc$read_h5ad("Covid_Lung_Endo.h5ad")
table(anndata$obs$SubCluster)

Counts <- anndata$T$to_df()
Meta <- anndata$obs
covid.m <- CreateSeuratObject(counts=Counts,
                              meta.data=Meta[,c('Cluster','SubCluster','Viral+','donor','disease','hUSI')])
table(covid.m$SubCluster)


### load clinical data
ClinicMeta = read.csv('/home/wangjing/wangj/AgingScore/GSE163530_COVID-19/GSE162911_GeoMx/ClinicMeta.csv',row.names = 1)
ClinicMeta <- ClinicMeta[complete.cases(ClinicMeta$Days_to_death),]
median(ClinicMeta$Days_to_death)
ClinicMeta$Progress = ifelse(ClinicMeta$Days_to_death < 15,'severe','moderate')
table(ClinicMeta$Progress)

covid.m$Donor = factor(gsub('_[12345]','',covid.m$donor))
covid.m$Age = factor(ClinicMeta[covid.m$Donor,'Age'],levels = c("30-35","40-45","50-55","55-60","60-65","65-70","75-80","80-85",">89"),ordered = T)
covid.m$Progress = factor(ClinicMeta[covid.m$Donor,'Progress'],levels = c('moderate','severe'),ordered = T)
covid.m$DTD = ClinicMeta[covid.m$Donor,'Days_to_death']

### umap of all endo cells
covid.m = ScaleData(covid.m) %>% FindVariableFeatures() %>% RunPCA()
ElbowPlot(covid.m)

covid.m = RunUMAP(covid.m,dims = 1:20)
DimPlot(covid.m,group.by = 'SubCluster')


### artery cells
subEndo = covid.m[,covid.m$SubCluster == 'Artery EC']
### ranking the genes based on each gene expression correlation with the aging score
corAging <- function(x,agingScore){cor <- cor(x,agingScore);cor}
cor.genes <- apply(GetAssayData(subEndo),1,corAging,agingScore=subEndo$hUSI)
cor.genes[is.na(cor.genes)] <- 0
features <- names(cor.genes)[order(abs(cor.genes),decreasing=TRUE)][1:1500]

subEndo.matrix <- as.matrix(GetAssayData(subEndo))[features,]
subEndo[['AgingExp']] <- CreateAssayObject(subEndo.matrix)
DefaultAssay(subEndo) <- "AgingExp"
subEndo$batch <- rep("batch",ncol(subEndo))

# ### Running the ICAnet
library(ICAnet)
source('~/wangj/codebase/HUSI/getPPI_String.R')

# PPI <- readRDS('PPI_subEndo.rds')
PPI <- getPPI_String(subEndo,species=9606)

Ica.epi <- ICAcomputing(subEndo,ICA.type="JADE",RMT=TRUE,two.stage=FALSE)
rownames(Ica.epi$ica.pooling) <- features
png('/home/wangjing/wangj/codebase/HUSI/Figures/ICA.png')
subEndo <- RunICAnet(subEndo,Ica.epi$ica.pooling,PPI.net = PPI,scale=FALSE,
                     ModuleSignificance = FALSE,cores = 1,aucMaxRank=300)
dev.off()
### diffusion map
library(destiny)
exp = subEndo@assays$IcaNet@data %>% as.matrix()
pd <- new('AnnotatedDataFrame', data = as.data.frame(subEndo@meta.data))
fData <- data.frame(gene_short_name = row.names(exp), row.names = row.names(exp))
fd <- new('AnnotatedDataFrame', data = fData)
myExpressionSet <- ExpressionSet(assayData=as.matrix(exp),
                                 phenoData=pd,
                                 annotation="EC")
set.seed(233)
dm <- DiffusionMap(myExpressionSet,n_pcs = 5)
dev.off()
plot(eigenvalues(dm), ylim = 0:1, pch = 20,xlab = 'Diffusion component (DC)', ylab = 'Eigenvalue')
dpt <- DPT(dm)
plot(dpt)

### diffusion embedding
subEndo[['diffusion']] <- CreateDimReducObject(data.frame(DC1=dm$DC1,DC2=dm$DC2,row.names = colnames(exp)) %>% as.matrix(),key="DC_")

root_cell = paste('DPT',grep(names(sort(dm$DC1+dm$DC2,decreasing = F)[1]),names(dm$DC1)),sep = '')
subEndo$DPT <- dpt[[root_cell]]

### perform tradeseq
library(tradeSeq)

counts = subEndo@assays$IcaNet@counts %>% as.matrix()
cellWeights = data.frame(lineage = rep(1,ncol(counts)),row.names = colnames(counts))
pseudotime <- matrix(subEndo$DPT, ncol = ncol(cellWeights),nrow = ncol(counts), byrow = FALSE)
ts <- fitGAM(counts = counts ,
             pseudotime = pseudotime,
             cellWeights = cellWeights,
             genes = rownames(counts),
             nknots = 6)

### select trajectory associated genes
assoRes <- associationTest(ts)
assoRes <- filter(assoRes,pvalue < 0.05)
assoRes <- assoRes[order(assoRes$meanLogFC,decreasing = T),]

### Icanet
ICAdown <- c("ICAnet-10-3-6", "ICAnet-10-1-6", "ICAnet-27-1-3", "ICAnet-16-3-3", "ICAnet-1-1-10", "ICAnet-21-4-4", "ICAnet-21-3-3", "ICAnet-3-28-8", "ICAnet-3-3-9", "ICAnet-3-1-8")

ICAdown_genes <- unique(unlist(subEndo@misc$IcaNet_geneSets[ICAdown]))
ICAdown_genes

# ICAup <- c('ICAnet-19-2-4','ICAnet-3-2-3','ICAnet-23-2-8','ICAnet-4-1-5')
ICAup <- c('ICAnet-22-1-4', 'ICAnet-2-6-7', 'ICAnet-17-7-6', 'ICAnet-17-3-5', 'ICAnet-15-1-5', 'ICAnet-15-13-4', 'ICAnet-17-6-5', 'ICAnet-4-3-4', 'ICAnet-1-1-11', 'ICAnet-4-4-7', 'ICAnet-2-5-5', 'ICAnet-15-4-6', 'ICAnet15-2-5', 'ICAnet-15-10-5', 'ICAnet-3-1-4', 'ICAnet14-2-3')
ICAup_genes <- unique(unlist(subEndo@misc$IcaNet_geneSets[ICAup]))
ICAup_genes









### other covid-19 data
meta <- read.csv('/home/wangjing/wangj/AgingScore/GSE171524_NormalLung/lung_metaData.txt',sep = '\t')
meta = meta[-1,]
rownames(meta) <- meta$NAME
cells_use = meta[meta$cell_type_fine == 'Arterial endothelial cells','NAME']
counts = read.csv('/home/wangjing/wangj/AgingScore/GSE171524_NormalLung/GSE171524_processed_data.csv',sep = ',')
counts_use = counts[,gsub('-','.',cells_use)]
rownames(counts_use) <- counts$X

df_plot = data.frame(t(counts_use[ICAup_genes,])) %>% rownames_to_column('cell') %>% 
          reshape2::melt(variable.name = 'gene',value.name = 'expression')
df_plot$cell <- gsub('\\.','-',df_plot$cell)
df_plot$gene <- gsub('\\.','-',df_plot$gene)

df_plot$condition = meta[df_plot$cell,'disease__ontology_label']
df_plot$DTD = meta[df_plot$cell,'interval_death_symptoms_onset_days']
df_plot$Donor = meta[df_plot$cell,'biosample_id']

ggplot(df_plot)+
  geom_boxplot(aes(x = gene,y = expression,fill = condition))

filter(df_plot,condition == 'COVID-19' & !is.na(DTD)) %>%
  mutate(DTD = as.numeric(DTD)) %>%
  mutate(Progress = ifelse(DTD<median(DTD),'severe','moderate')) %>%
  ggplot()+
  geom_boxplot(aes(x = gene,y = expression,fill = Progress))



