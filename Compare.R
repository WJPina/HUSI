mm_l2 = readRDS("/home/wangjing/wangj/AgingScore/Data/Bulk_TrainModel/mm_l2.rds")
######################## Comparision in 4 single-cell dataset  ###################################
setwd("/home/wangjing/wangj/AgingScore/Comparison")

library(Seurat)
library(data.table)
library(tibble)
library(sparseMatrixStats)
library(org.Hs.eg.db)
library(dplyr)
library(caret)
library(magrittr)
library(reshape2)

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
    if(type == 'laf'){
        for(sl in names(scorelist)){
                score = calc_laf(GetAssayData(object) %>% as.matrix, scorelist[[sl]])
                data[[sl]] = score
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

calc_auc <- function(data=data.frame(),type = '',scorelist=NULL){
    fold = createMultiFolds(data$condition, k = 10, times = 3)
    auc = sapply(fold, function(sampling) {
        df = data[sampling,]
        if(type=="laf"){
            c("DAS","mSS",'DAS_mSS','hUSI') %>% 
            sapply(function(idx) {pROC::roc(df$condition, df[,idx], levels=c("Growing", "Senescence")) %>% pROC::auc()})
        }
        else if (type=='marker'){
            c('hUSI',scorelist[scorelist %in% colnames(data)]) %>% 
            sapply(function(idx) {pROC::roc(df$condition, df[,idx], levels=c("Growing", "Senescence")) %>% pROC::auc()})
        }
        else{
            c('hUSI',names(scorelist)) %>% 
            sapply(function(idx) {pROC::roc(df$condition, df[,idx], levels=c("Growing", "Senescence")) %>% pROC::auc()})
        }
    })
    return(auc)
}

getGenes <- function(type) {
    if(type=='laf'){
        scorelist=lafSet
    }
    else if(type=='ssgsea'){
        scorelist=EnrichSet
    }
    else if(type=='marker'){
        scorelist=SenMarkers
    }
    else{
        scorelist=NULL
    }
    return(scorelist)
}
### score list
laf_DAS = fread("LaffertyWhyte2010_DAS.csv")
laf_mSS = fread("LaffertyWhyte2010_mSS.csv")
lafSet <<- list(DAS = laf_DAS,mSS=laf_mSS,DAS_mSS=rbind(laf_DAS, laf_mSS))

SenMarkers <<- c("GLB1", "TP53", "CDKN1A", "CDKN2A", "LMNB1", "IL1A", "RB1", "CDK1", "CDK4", "MKI67", "CDKN2B")

EnrichSet<<-cogena::gmt2list("/mnt/data2/zhouxl/Pan_Cancer/Data/Signatures/gene_50signatures_merge.gmt")
EnrichSet=EnrichSet[43:49]
names(EnrichSet)
rep_sene_genes = readxl::read_excel("SigRS.xls", sheet = 5, skip = 2)$Gene_names
EnrichSet$Sig.RS = rep_sene_genes
names(EnrichSet) 
names(EnrichSet) <- c("SenMayo","CellAge", "GenAge", "ASIG", "SASP","AgingAtlas", "SenUp","SigRS")

Results = list()

### Teo2019  GSE115301 IMR90 Smart-seq2 OIS 
Teo2019 = CreateSeuratObject(
    fread("Teo2019/GSE115301_Growing_Sen_10x_count.txt.gz") %>% column_to_rownames("V1") %>% data.matrix, 
    meta.data = fread("Teo2019/GSE115301_Growing_Sen_10x_metadata.txt.gz", header = T) %>% column_to_rownames("V1")
)
### log
Teo2019 = NormalizeData(Teo2019)
### CPM
# Teo2019 = NormalizeData(Teo2019,normalization.method = 'RC',scale.factor = 1e6)

hUSI = GetAssayData(Teo2019) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}
dat_laf = calc_scores(lafSet,Teo2019,'laf')
dat_marker = calc_scores(SenMarkers,Teo2019,'marker')
dat_ssgsea = calc_scores(EnrichSet,Teo2019,'ssgsea')

Teo2019List = list(laf=dat_laf,marker=dat_marker,ssgsea=dat_ssgsea)

Results$Teo2019 = Teo2019List
Results$Teo2019$hUSI = hUSI

### Tang2019 GSE119807 HCA2 fibroblast cell RS,IRIS
Tang2019List = list.files("Tang2019", full.names = T) %>% 
    lapply(function(x) {
        scdat = fread(x) %>% 
            column_to_rownames("GENE") %>% 
            data.matrix %>% 
            CreateSeuratObject() %>% 
            # NormalizeData(normalization.method = 'RC',scale.factor = 1e6)
            NormalizeData()

        hUSI = GetAssayData(scdat) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}
        data <- list(laf=data.frame(calc_scores(lafSet,scdat,'laf'),'hUSI'=hUSI),
                        marker=data.frame(calc_scores(SenMarkers,scdat,'marker'),'hUSI'=hUSI),
                        ssgsea=data.frame(calc_scores(EnrichSet,scdat,'ssgsea'),'hUSI'=hUSI))
        return(data)
    }) %>% 
    set_names( list.files("Tang2019") %>% gsub(".*_", "", .) %>% gsub("\\..*", "", .) )

Results$Tang2019 = Tang2019List

### Aarts2017 GSE94980 IMR90 OSKM-expressing reprogramming-induced senescence
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

Aarts2017 = CreateSeuratObject(exp) %>% 
                # NormalizeData(normalization.method = 'RC',scale.factor = 1e6)
                NormalizeData()
Aarts2017$Condition = ifelse(grepl('OSKM',colnames(Aarts2017)), 'Senescence', 'Growing')
hUSI = GetAssayData(Aarts2017) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}

dat_laf = calc_scores(lafSet,Aarts2017,'laf')
dat_marker = calc_scores(SenMarkers,Aarts2017,'marker')
dat_ssgsea = calc_scores(EnrichSet,Aarts2017,'ssgsea')

Aarts2017List = list(laf=dat_laf,marker=dat_marker,ssgsea=dat_ssgsea)

Results$Aarts2017 = Aarts2017List
Results$Aarts2017$hUSI = hUSI

### Enge2017 GSE81547 human pancreas organisms age
Enge2017List = list.files("Enge2017", pattern = "*.gz",full.names = T)
meta = read.csv("Enge2017/metadata.txt",row.names = 1)
table(meta$DONOR_AGE)
meta = filter(meta, DONOR_AGE %in% c(1,54))

samples = lapply(strsplit(basename(Enge2017List),"_"),'[',1) %>% unlist 
files = Enge2017List[samples %in% meta$Sample.Name]

scdat = data.frame(gene = '')
for(file in files){
    scol = fread(file) 
    colnames(scol) = c('gene',strsplit(basename(file),"_")[[1]][1])
    scdat = right_join(scdat,scol,by = 'gene')           
}
scdat = column_to_rownames(scdat,"gene") %>% data.matrix
dim(scdat)
Enge2017 = CreateSeuratObject(scdat) %>% 
                # NormalizeData(normalization.method = 'RC',scale.factor = 1e6)
                NormalizeData()
Enge2017$Condition = ifelse(meta$DONOR_AGE %in% c(54), 'Senescence', 'Growing')
hUSI = GetAssayData(Enge2017) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}

dat_laf = calc_scores(lafSet,Enge2017,'laf')
dat_marker = calc_scores(SenMarkers,Enge2017,'marker')
dat_ssgsea = calc_scores(EnrichSet,Enge2017,'ssgsea')

Enge2017List = list(laf=dat_laf,marker=dat_marker,ssgsea=dat_ssgsea)

Results$Enge2017 = Enge2017List
Results$Enge2017$hUSI = hUSI


AUClist = list()
for(dataset in names(Results)){
    for(type in c('laf','marker','ssgsea')){
        scorelist = getGenes(type)
        if(dataset == 'Teo2019'){
            dat = Results[[dataset]][[type]] %>% 
                    mutate(condition=gsub("[0-9]$", "", Condition1) ) %>% 
                    mutate(hUSI=Results[[dataset]][['hUSI']])
            dat = dat[,complete.cases(t(dat))]

            auc <- calc_auc(dat,type,scorelist)
        }
        else if (dataset == 'Tang2019') {
            dat = data.table::rbindlist(lapply(Results[[dataset]],function(x){x[[type]]}),use.names=TRUE, fill=TRUE, idcol="sample") %>% data.frame
            dat = dat[,complete.cases(t(dat))]
            dat$condition = as.factor(ifelse(dat$sample %in% c("senescence","LowPD50Gy"),'Senescence','Growing'))

            auc <- calc_auc(dat,type,scorelist)
        }
        else if (dataset %in% c('Aarts2017','Enge2017')) {
            dat = Results[[dataset]][[type]] %>% 
                    mutate(condition=Condition) %>% 
                    mutate(hUSI=Results[[dataset]][['hUSI']])
            dat = dat[,complete.cases(t(dat))]

            auc <- calc_auc(dat,type,scorelist)
        }
        AUClist[[dataset]][[type]] = auc
    }
}







