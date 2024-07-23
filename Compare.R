mm_l2 = readRDS("~/wangj/AgingScore/Data/Bulk_TrainModel/mm_l2.rds")
setwd('/home/wangjing/wangj/AgingScore/Comparison')
######################## Comparision in three single-cell dataset  ###################################
library(Seurat)
library(data.table)
library(tibble)
library(sparseMatrixStats)
library(org.Hs.eg.db)
library(dplyr)
library(caret)
library(magrittr)
library(reshape2)
source('~/wangj/codebase/HUSI/functions.R')

set.seed(233)

calc_laf <- function(mat, genes) {
    genes = genes %>% filter(Gene %in% rownames(mat))
    x = mat[genes$Gene,]
    x = t( t(x)/colMedians(mat) ) > 1
    res = apply(x, 2, function(y) {
        purrr::invoke_map_lgl(
            .f = as.list( ifelse(genes$Direction=="Increase", "isTRUE", "isFALSE") ), 
            .x = as.list(y)
        )
    })
    colSums(res)/nrow(res)
}

calc_scores <- function(scorelist=NULL,object=NULL,type=''){
    data = object@meta.data
    if(type == 'method'){
        for(sl in names(scorelist)){
            if(sl == 'laf'){
                for (s in names(scorelist[[sl]])) {
                    score = calc_laf(GetAssayData(object) %>% as.matrix, scorelist[[sl]][[s]])
                    data[[s]] = score
                }
            }
            if(sl == 'lassoCS'){
                mat = object[names(lassoCS),] %>% GetAssayData %>% as.matrix 
                mat = apply(mat,2,function(x) {x * lassoCS[rownames(mat)]})
                tryCatch({score = colSums(mat) %>% data.frame 
                                    rownames(score) = colnames(object)}, 
                            error = function(e) {
                                score = mat %>% data.frame})
                colnames(score) = 'lassoCS'   
                data = cbind(data,score)
            }
            if(sl == 'CSS'){
                mat = object[names(CSS),] %>% GetAssayData %>% as.matrix 
                mat = apply(mat,2,function(x) {x * CSS[rownames(mat)]}) 
                tryCatch({score = colSums(mat) %>% data.frame 
                                    rownames(score) = colnames(object)}, 
                            error = function(e) {
                                score = mat %>% data.frame})
                colnames(score) = 'CSS'
                data = cbind(data,score)
            }
        }
    }
    if(type == 'ssgsea'){
        for(sl in names(scorelist)){
            score = GSVA::gsva(GetAssayData(object) %>% as.matrix, list(signature=scorelist[[sl]]), method="ssgsea", ssgsea.norm = TRUE, verbose = TRUE)[1,]
            data[[sl]] = score
        }
    }
    if(type == 'marker'){
        score = object[scorelist,] %>% GetAssayData %>% as.matrix %>% t %>% data.frame 
        rownames(score) = colnames(object)
        data = cbind(data,score)
    }
    return(data)
} 


### score list
SenMarkers <<- c("GLB1", "TP53", "CDKN1A", "CDKN2A", "LMNB1", "IL1A", "RB1", "CDK1", "CDK4","CDK6", "MKI67", "CDKN2B",'SERPINE1')

EnrichSet<<-cogena::gmt2list("gene_50signatures_merge.gmt")

EnrichSet=EnrichSet[43:49]
names(EnrichSet)
rep_sene_genes = readxl::read_excel("SigRS.xls", sheet = 5, skip = 2)$Gene_names
EnrichSet$SigRS = rep_sene_genes
names(EnrichSet) 
names(EnrichSet) <- c("SenMayo","CellAge", "GenAge", "ASIG", "SASP","AgingAtlas", "SenUp","SigRS")

laf_DAS = fread("LaffertyWhyte2010_DAS.csv")
laf_mSS = fread("LaffertyWhyte2010_mSS.csv")
lafSet = list(DAS =laf_DAS,mSS=laf_mSS,`DAS+mSS`=rbind(laf_DAS, laf_mSS))
lassoCS = c(-0.158,-0.153,0.17,0.347,-0.215,0.288,-0.196,-0.245,-0.103,0.127)
names(lassoCS) = c("ITGA8","SEMA3G","DPYSL3","IFITM1","ZNF521", "SOCS3","PCSK6", "DUSP1", "SLC44A4","IL20RB")
CSS =  c(0.2921,0.1810,0.0524,0.1285)
names(CSS) = c('C1QTNF6','SQOR','LYPD3','FAM83B')
Methods <<- list('laf' = lafSet,'lassoCS'=lassoCS,'CSS'=CSS)

### Teo2019  GSE115301 IMR90 Smart-seq2 OIS,3',counts
Teo2019 = CreateSeuratObject(
    fread("Teo2019/GSE115301_Growing_Sen_10x_count.txt.gz") %>% column_to_rownames("V1") %>% data.matrix, 
    meta.data = fread("Teo2019/GSE115301_Growing_Sen_10x_metadata.txt.gz", header = T) %>% column_to_rownames("V1")
) %>% NormalizeData()


### cell states in Teo dataset
Teo2019 <- FindVariableFeatures(Teo2019) %>% ScaleData() %>% RunPCA() %>% RunTSNE(dims=1:15)
Teo2019 <-  FindNeighbors(Teo2019,dims = 1:15,reduction = "pca") %>% FindClusters(resolution = 0.1)
DimPlot(Teo2019,group.by = 'Condition2')
Teo2019$group = ifelse(as.character(Teo2019$seurat_clusters) == 0,'cluster0:Growing',
                       ifelse(as.character(Teo2019$seurat_clusters)== 1,'cluster1:OIS','cluster2:Secondary senescence'))


###  auc
Teo2019$Condition = as.factor(ifelse(Teo2019$Condition2 == 'RIS','Senescence','Growing'))

hUSI = GetAssayData(Teo2019) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )} %>% minmax


dat_method = calc_scores(Methods,Teo2019,'method')
dat_marker = calc_scores(SenMarkers,Teo2019,'marker')
dat_ssgsea = calc_scores(EnrichSet,Teo2019,'ssgsea')

dat_Teo2019 = cbind(dat_method,dat_marker,dat_ssgsea)
dat_Teo2019$hUSI = hUSI

### Tang2019 GSE119807 HCA2 fibroblast cell RS,IRIS,3',counts
Tang2019List = list.files("Tang2019", full.names = T) %>% 
  lapply(function(x) {
    scdat = fread(x) %>% 
      column_to_rownames("GENE") %>% 
      data.matrix }) %>% set_names( list.files("Tang2019") %>% gsub(".*_", "", .) %>% gsub("\\..*", "", .) )

genes = Reduce(intersect,lapply(Tang2019List, function(x)rownames(x)))
Tang2019 = do.call(cbind,lapply(Tang2019List,function(x) x <- x[genes,]))
Tang2019 = CreateSeuratObject(Tang2019) %>% NormalizeData()
Tang2019$Condition = rep(names(Tang2019List),each=400)
Tang2019$Condition = as.factor(ifelse(Tang2019$Condition %in% c("senescence","LowPD50Gy"),'Senescence','Growing'))


hUSI = GetAssayData(Tang2019) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )} %>% minmax

dat_method = calc_scores(Methods,Tang2019,'method')
dat_marker = calc_scores(SenMarkers,Tang2019,'marker')
dat_ssgsea = calc_scores(EnrichSet,Tang2019,'ssgsea')

dat_Tang2019 = cbind(dat_method,dat_marker,dat_ssgsea)
dat_Tang2019$hUSI = hUSI


### Aarts2017 GSE94980 IMR90 OSKM-expressing reprogramming-induced senescence,3',counts
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
Aarts2017List = list.files("Aarts2017", pattern = "*.gz",full.names = T)
scdat = data.frame(gene = '')
for(file in Aarts2017List){
    scol = fread(file) 
    colnames(scol) = c('gene',paste(strsplit(basename(file),"_")[[1]][c(2:5)],collapse='_'))
    scdat = right_join(scdat,scol,by = 'gene')           
}
ids=bitr(scdat$gene,'ENSEMBL','SYMBOL','org.Hs.eg.db')
exp=merge(scdat,ids,by.x='gene',by.y='ENSEMBL')
exp = exp[-1]
exp = exp[!duplicated(exp$SYMBOL),]
rownames(exp) <- exp$SYMBOL
exp = exp[-ncol(exp)]
dim(exp)

Aarts2017 = CreateSeuratObject(exp) %>% NormalizeData()
Aarts2017$Condition = ifelse(grepl('OSKM',colnames(Aarts2017)), 'Senescence', 'Growing')
hUSI = GetAssayData(Aarts2017) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )} %>% minmax

dat_method = calc_scores(Methods,Aarts2017,'method')
dat_marker = calc_scores(SenMarkers,Aarts2017,'marker')
dat_ssgsea = calc_scores(EnrichSet,Aarts2017,'ssgsea')

dat_Aarts2017 = cbind(dat_method,dat_marker,dat_ssgsea)
dat_Aarts2017$hUSI = hUSI

### calculate auc
dat = list(dat_Teo2019,dat_Tang2019,dat_Aarts2017)
names(dat) = c('Teo2019','Tang2019','Aarts2017')

### only senescence markers

AUClist=list()
for(dataset in names(dat)){
  data = dat[[dataset]]
  fold = createMultiFolds(data$Condition, k = 10, times = 3)
  auc = sapply(fold, function(sampling){
    df = data[sampling,]
    auc_1 = c('hUSI',"DAS","mSS",'DAS+mSS','lassoCS','CSS',names(EnrichSet)) %>%
      sapply(function(idx) {pROC::roc(df$Condition, df[,idx], levels=c("Growing", "Senescence"),direction = '<') %>% pROC::auc()})
    auc_2 = c(SenMarkers[SenMarkers %in% colnames(df)]) %>% 
      sapply(function(idx) {pROC::roc(df$Condition, df[,idx], levels=c("Growing", "Senescence"))%>% pROC::auc()})
    auc = c(auc_1,auc_2)
    return(auc)
  })
  AUClist[[dataset]] = auc
}

#################################### comparision in GTEx and TCGA #############################
##### GTEx bulk data
### load cs score
TCSER_gtex = read.csv('Senescence_quantification_GTEX.csv')
dim(TCSER_gtex)
TCSER_gtex = TCSER_gtex[!duplicated(TCSER_gtex$Sample),]
rownames(TCSER_gtex) = TCSER_gtex$Sample
library(CePa)
GTEx_exp <- read.gct("bulk-gex_v8_rna-seq_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
dim(GTEx_exp)
GTEx_exp[1:5,1:5]
table(rownames(TCSER_gtex) %in% gsub('\\.','-',colnames(GTEx_exp)))
colnames(GTEx_exp) <- gsub('\\.','-',colnames(GTEx_exp))
Counts_gtex = GTEx_exp[,intersect(rownames(TCSER_gtex),colnames(GTEx_exp))]
Counts_gtex = esembl2symbol(counts = Counts_gtex,version = '26')
dim(Counts_gtex)
Counts_gtex[1:5,1:5]

Counts_gtex = log2(Counts_gtex + .00001) %>% scale()
hUSI_gtex <- Counts_gtex %>% apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )
### gene sets
scoreList_GTEx = list()
scoreList_GTEx[['TCSER']] = TCSER_gtex[colnames(Counts_gtex),]$CS.score
scoreList_GTEx[['hUSI']] = hUSI_gtex

### gene sets
for(sl in names(EnrichSet)){
  score = GSVA::gsva(Counts_gtex, list(signature=EnrichSet[[sl]]), method="ssgsea", ssgsea.norm = TRUE, verbose = TRUE)[1,]
  scoreList_GTEx[[sl]] = score
}

###### TCGA bulk data
### load cs score
TCSER_tcga = read.csv('Senescence_quantification_TCGA.csv')
dim(TCSER_tcga)
TCSER_tcga = TCSER_tcga[!duplicated(TCSER_tcga$sample),]
rownames(TCSER_tcga) = TCSER_tcga$sample

TCGA_exp = fread('tcga_RSEM_gene_tpm.gz') 
TCGA_exp = data.frame(TCGA_exp)
TCGA_exp = column_to_rownames(TCGA_exp,'sample')
dim(TCGA_exp)
TCGA_exp[1:5,1:5]

rownames(TCSER_tcga) = gsub('-','.',TCSER_tcga$sample)
Counts_tcga = TCGA_exp[,rownames(TCSER_tcga)] %>% as.matrix
Counts_tcga = esembl2symbol(version = '23',counts = Counts_tcga)
dim(Counts_tcga)
Counts_tcga[1:5,1:5]

hUSI_tcga <- Counts_tcga %>% apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )

for(sl in names(EnrichSet)){
  score = GSVA::gsva(Counts_tcga, list(signature=EnrichSet[[sl]]), method="ssgsea", ssgsea.norm = TRUE, verbose = TRUE)[1,]
  scoreList_TCGA[[sl]] = score
}

cormat_tcga <- do.call(data.frame, scoreList_TCGA) %>% cor(method = 'spearman')

### validate hUSI in GTEx and TCGA by meta data
### Age
GTEx_sample <- read.csv("/mnt/data1/wangj/MyProject/Data/GTEx/Bulk_SnRNA/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",sep = "\t",row.names = 1,header = T)
GTEx_sample<- GTEx_sample[,c("SMTS","SMTSD")]
GTEx_sample$Sample <- rownames(GTEx_sample)
GTEx_sample$Patient <- lapply(strsplit(rownames(GTEx_sample),'-'),function(x) {paste(x[1],x[2],sep = '-')}) %>% unlist
length(unique(GTEx_sample$Patient))

GTEx_Pheno <- read.table("/mnt/data1/wangj/MyProject/Data/GTEx/Bulk_SnRNA/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",sep = "\t",header = T)
length(unique(GTEx_Pheno$SUBJID))
GTEx_meta <- merge(GTEx_Pheno,GTEx_sample,by.x = "SUBJID",by.y = "Patient")
length(unique(GTEx_meta$SUBJID))

table(names(scoreList_GTEx$hUSI) %in% GTEx_meta$Sample)
rownames(GTEx_meta) = GTEx_meta$Sample
