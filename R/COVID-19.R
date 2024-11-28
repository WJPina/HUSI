### This script is used to analyze the COVID-19 data

library(Seurat)
library(mclust)

setwd('hUSI/')

### Load the expression data
files = list.files('GSE171524_NormalLung/Counts/') ### download the data from GEO
countsList = lapply(files, function(x) read.csv(paste0('GSE171524_NormalLung/Counts/',x), header = T, row.names = 1))
genes = Reduce(intersect, lapply(countsList, rownames))
countsList = lapply(countsList, function(x) x[genes,])
counts = do.call(cbind, countsList)
dim(counts)
counts[1:5,1:5]
### Load meta data 
meta = read.csv('GSE171524_NormalLung/lung_metaData.txt', header = T,row.names=1,sep = '\t')
dim(meta)
meta[1:5,1:5]
meta = meta[-1,]
rownames(meta) = gsub('-','\\.',rownames(meta))
table(rownames(meta) %in% colnames(counts))
colnames(counts)[!colnames(counts) %in% rownames(meta)]

### create object
covid <- CreateSeuratObject(counts = counts[,rownames(meta)], meta.data = meta, project = 'COVID-19')
covid <- NormalizeData(covid)

covid <- FindVariableFeatures(covid, selection.method = "vst", nfeatures = 2000)
covid <- ScaleData(covid)
covid <- RunPCA(covid, npcs = 30, verbose = F)
covid <- RunUMAP(covid, reduction = "pca", dims = 1:30)

### load umap data
umap = read.csv('Data/lung_clusterfile.txt', header = T, row.names = 1, sep = '\t')
dim(umap)
umap[1:5,]
umap = umap[-1,]
rownames(umap) = gsub('-','\\.',rownames(umap))
umap = umap[colnames(covid),]
colnames(umap) = c('UMAP_1','UMAP_2')
umap$UMAP_1 = as.numeric(as.character(umap[,1]))
umap$UMAP_2 = as.numeric(as.character(umap[,2]))
umap = as.matrix(umap) 
covid[['umap']] <- CreateDimReducObject(umap,key="UMAP_")
covid$cell_type_main = factor(covid$cell_type_main,levels = c("APC-like","Endothelial cells","Epithelial cells","Fibroblasts","Mast cells","Myeloid","Neuronal cells","B cells","T cells" ),ordered = T)

DimPlot(covid, group.by = 'cell_type_main', label = T) + NoLegend()

meta0=covid@meta.data
covid@meta.data=covid@meta.data[,!(colnames(covid@meta.data) %in% c('hUSI','SenClass_fine','cell_type_SenClass_fine','cell_type_SenClass','SenClass_merge'))]

### calculate aging score
source('R/functions_new.R')
load('Data/SenOCLR_l2=1_drop.rdata')

hUSI  = scoreOCLR(GetAssayData(covid),SenOCLR,'spearman')[[1]] %>% minmax()

covid$hUSI = hUSI

myorder=rev(names(table(covid$cell_type_main)))
covid$cell_type_main = factor(covid$cell_type_main,levels = myorder,ordered = T)
df_plot = covid@meta.data[order(covid@meta.data$cell_type_main),]
covid$cell_type_fine <- factor(covid$cell_type_fine,levels = unique(df_plot$cell_type_fine),ordered = T)

### divided cells by gaussion
df = covid@meta.data
SenClass = c()
for(c in unique(df$cell_type_main)){
  print(c)
  data = filter(df,cell_type_main==c)
  husi = data$hUSI
  set.seed(123)
  model = husi %>% Mclust()
  meanV = model$parameters$mean
  classification = cut(husi, breaks = c(-Inf, meanV, Inf), include.lowest = TRUE) %>% as.integer()
  print(table(classification))
  names(classification) <- rownames(data)
  SenClass <- c(SenClass,classification)
}

covid$SenClass_main = paste('C',SenClass[colnames(covid)],sep = '')


### divided cells by gaussion in cell type fine
df = covid@meta.data
SenClass = c()
SenClass_merge = c()
for(c in unique(df$cell_type_fine)){
  print(c)
  data = filter(df,cell_type_fine==c)
  husi = data$hUSI

  classification = class_hUSI(husi)
  names(classification) <- rownames(data)
  
  SenClass <- c(SenClass,classification)
  print(table(classification))
}

covid$SenClass_fine = ifelse(SenClass[colnames(covid)]=='0','Normal cells','Senescent cells')

################################# DEGs #########################################
DEGs = {}
celltype_use =levels(covid$cell_type_fine)
for(c in celltype_use){
  print(c)
  covid_sub = subset(covid,cell_type_fine==c)
  print(table(covid_sub$SenClass_fine))
  Idents(covid_sub) = 'SenClass_fine'
  
  DEG = FindMarkers(covid_sub,ident.1 = 'Senescent cells',ident.2='Normal cells',logfc.threshold = 0.5,only.pos = T,min.diff.pct = 0)
  DEG = filter(DEG,p_val_adj<0.05)
  if(nrow(DEG)>0){
    DEG$cluster = c
    DEG$gene = rownames(DEG)
    DEGs[[c]] <- DEG
  }
}
sapply(DEGs, dim)


ribosomal = read.table("Data/Ribosome.txt",stringsAsFactors=FALSE)

pathways = read.gmt('Data/KEGG_hsa.gmt')
pathwaysList = list()
for(c in celltype_use){
  print(c)
  genes = DEGs[[c]]$gene
  genes = genes[!genes%in%ribosomal$V1]
  tryCatch({
    enrich = enricher(gene = genes,pvalueCutoff = 0.05,
                      universe = rownames(covid@assays$RNA),TERM2GENE = pathways)
    re = enrich@result
    re$celltype = c
    pathwaysList[[c]] <- re
  },error = function(e){print('No enrichment!!!')})
}

############################## proportion changes ##############################
df = covid@meta.data %$% table(group,cell_type_fine,SenClass_fine) %>% as.array() %>% data.frame()
df = filter(df,Freq !=0)

head(df)

df_summary <- df %>%
  group_by(group,cell_type_fine) %>%         
  mutate(Total = sum(Freq), SenFrac = Freq / Total) %>%  
  ungroup() 

df_diff <- df_summary %>%
  group_by(cell_type_fine,SenClass_fine) %>%
  mutate(Frac_diff = SenFrac[2]-SenFrac[1]) %>%
  ungroup()

df_diff <- df_diff[complete.cases(df_diff),]
df_diff= filter(df_diff,Frac_diff!=0)
head(df_diff[order(df_diff$Frac_diff,decreasing = T),],20)

########## correlation between days to death and senescent cells ###############
info = covid@meta.data[,c('interval_death_symptoms_onset_days','age','donor_id')]
info = info[!duplicated(info),]
info = info[complete.cases(info),]
info$DTD = as.numeric(info$interval_death_symptoms_onset_days)
info$Age = as.numeric(info$age)
rownames(info) = info$donor_id
info

df = covid@meta.data %$% table(cell_type_fine,SenClass_fine,donor_id) %>% as.array() %>% data.frame()
head(df)
df = filter(df,Freq !=0)
df = filter(df,donor_id %in% info$donor_id)
df$DTD = info$DTD[match(df$donor_id,info$donor_id)]
df$Age = info$Age[match(df$donor_id,info$donor_id)]
head(df)

df_summary <- df %>%
  group_by(donor_id,cell_type_fine) %>%         
  mutate(Total = sum(Freq), SenFrac = Freq / Total) %>%  
  ungroup()     

df = covid@meta.data
df = filter(df,group=='COVID-19'&interval_death_symptoms_onset_days!=''&age!='')
df$DTD = as.numeric(df$interval_death_symptoms_onset_days)

table(df$DTD)
DTD_s = unique(df$DTD) %>% quantile(0.25)
DTD_m = unique(df$DTD) %>% quantile(0.75)
df$Prognosis = factor(ifelse(df$DTD<=DTD_s ,'Severe',ifelse(df$DTD>=DTD_m,'Moderate','Other')))
df = filter(df,DTD!='Other')
df$Prognosis = factor(df$Prognosis ,levels = c('Severe','Moderate'),ordered = T)
table(df$Prognosis)
table(df$Prognosis,df$donor_id)

df = df %$% table(Prognosis,cell_type_fine,SenClass_fine) %>% as.array() %>% data.frame()
df = filter(df,Freq !=0)
head(df)

df_summary <- df %>%
  group_by(Prognosis,cell_type_fine) %>%         
  mutate(Total = sum(Freq), SenFrac = Freq / Total) %>%  
  ungroup() 

df_diff <- df_summary %>%
  group_by(cell_type_fine,SenClass_fine) %>%
  mutate(Frac_diff = SenFrac[2]-SenFrac[1]) %>%
  ungroup()

df_diff <- df_diff[complete.cases(df_diff),]
df_diff= filter(df_diff,Frac_diff!=0)
head(df_diff[order(df_diff$Frac_diff,decreasing = T),],20)



