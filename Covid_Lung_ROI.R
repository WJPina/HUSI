### deconvolute spatial ROI data
library(stringr)
library(stringi)
library(SpatialDecon)
library(Seurat)
library(dplyr)

load('/mnt/data1/wangj/AgingScore/GSE163530_COVID-19/GSE171668_scnRNA/Lung.RData')

profile <- create_profile_matrix(mtx = as.matrix(obj@assays$RNA@counts),            
                                 cellAnnots = obj@meta.data %>% mutate(CellID = rownames(.)),  
                                 cellTypeCol = "State",  
                                 cellNameCol = "CellID",           
                                 matrixName = "custom_state_profile", 
                                 outDir = NULL,                    
                                 normalize = FALSE,                
                                 minCellNum = 5,                   
                                 minGenes = 10,                    
                                 scalingFactor = 5,               
                                 discardCellTypes = TRUE)          


### Load ROI WTA data
setwd("/mnt/data1/wangj/AgingScore/GSE163530_COVID-19/GSE162911_GeoMx/")
wta.counts <- read.csv("Broad-COVID_WTA_Q3Norm_TargetCountMatrix.txt", row.names=1,header = T,sep = '\t')
wta.assay <- CreateAssayObject(wta.counts) %>% NormalizeData()
Bulk <- as.matrix(wta.assay@data)

df_segments <- read.csv("Broad-COVID_WTA_SegmentProperties.txt", row.names=1,header = T,sep = '\t')
rownames(df_segments) <- str_replace_all(rownames(df_segments), '-', '.')

df_tissue <- read.table("annotation_file_wta.txt", sep = "\t",row.names=1,header = T)
rownames(df_tissue) <- str_replace_all(rownames(df_tissue), '-', '.')

meta <- cbind(df_segments,df_tissue)
meta$donor <- as.character(unlist(lapply(strsplit(meta$scan.name,"-"),"[",1)))

meta$condition <- substr(meta$donor,1,1)
meta$donor <- stri_replace_all_regex(meta$donor,c("C01","C2","C3","S01","S02","S03","S09","S10","S11","S16","S18","S28"),
                                     c("D22","D23","D24","D18","D19","D20","D21","D8","D9","D10","D11","D12"),vectorize = F)

meta$DTD = ClinicMeta[meta$donor,'Days_to_death']
meta$Age = ClinicMeta[meta$donor,'Age']

meta$Age = gsub('>89','90-100',meta$Age)
meta$Age_min = strsplit(meta$Age,'-') %>% sapply(function(x) {as.numeric(x[1])})
meta$Age_max = strsplit(meta$Age,'-') %>% sapply(function(x) {as.numeric(x[2])})
meta$Age_mean = (meta$Age_min+meta$Age_max)/2
meta$group = ifelse(meta$DTD>16.5,ifelse(meta$Age_mean < 67.5,'middle_moderate','old_moderate'),
                    ifelse(meta$Age_mean < 67.5,'middle_severe','old_severe'))
table(meta$group)

meta.use = filter(meta,condition == 'S')
meta.use = filter(meta.use,Primary_Morph %in% c('Lymph node','Artery','Vein'))
table(meta.use$Primary_Morph)

Bulk.use = Bulk[,rownames(meta.use)]
bg = derive_GeoMx_background(norm = Bulk.use,
                             probepool = rep(1, nrow(Bulk.use)),
                             negnames = c("Neg Probe","SARS-CoV-2 Neg"))
res = spatialdecon(norm = Bulk.use,
                   bg = bg,
                   X = profile,
                   align_genes = TRUE)




