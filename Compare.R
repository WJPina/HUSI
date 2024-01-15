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

calc_auc <- function(data=data.frame(),type = '',scorelist=NULL){
    fold = createMultiFolds(data$condition, k = 10, times = 3)
    auc = sapply(fold, function(sampling) {
        df = data[sampling,]
        if(type == c("method")){
            c("DAS","mSS",'DAS+mSS','lassoCS','CSS','hUSI') %>% 
            sapply(function(idx) {pROC::roc(df$condition, df[,idx], levels=c("Growing", "Senescence"),direction = '<') %>% pROC::auc()})
        }
        else if (type=='marker'){
            c('hUSI',scorelist[scorelist %in% colnames(data)]) %>% 
            sapply(function(idx) {pROC::roc(df$condition, df[,idx], levels=c("Growing", "Senescence")) %>% pROC::auc()})
        }
        else{
            c('hUSI',names(scorelist)) %>% 
            sapply(function(idx) {pROC::roc(df$condition, df[,idx], levels=c("Growing", "Senescence"),direction = '<') %>% pROC::auc()})
        }
    })
    return(auc)
}

getGenes <- function(type) {
    if(type=='method'){
        scorelist=Methods
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
SenMarkers <<- c("GLB1", "TP53", "CDKN1A", "CDKN2A", "LMNB1", "IL1A", "RB1", "CDK1", "CDK4","CDK6", "MKI67", "CDKN2B")

EnrichSet<<-cogena::gmt2list("/mnt/data2/zhouxl/Pan_Cancer/Data/Signatures/gene_50signatures_merge.gmt")
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

Results = list()

### Teo2019  GSE115301 IMR90 Smart-seq2 OIS,3',counts
Teo2019 = CreateSeuratObject(
    fread("Teo2019/GSE115301_Growing_Sen_10x_count.txt.gz") %>% column_to_rownames("V1") %>% data.matrix, 
    meta.data = fread("Teo2019/GSE115301_Growing_Sen_10x_metadata.txt.gz", header = T) %>% column_to_rownames("V1")
) %>% NormalizeData()

hUSI = GetAssayData(Teo2019) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}
dat_method = calc_scores(Methods,Teo2019,'method')
dat_marker = calc_scores(SenMarkers,Teo2019,'marker')
dat_ssgsea = calc_scores(EnrichSet,Teo2019,'ssgsea')

Teo2019List = list(method=dat_method,marker=dat_marker,ssgsea=dat_ssgsea)
Results$Teo2019 = Teo2019List
Results$Teo2019$hUSI = hUSI

### Tang2019 GSE119807 HCA2 fibroblast cell RS,IRIS,3',counts
Tang2019List = list.files("Tang2019", full.names = T) %>% 
    lapply(function(x) {
        scdat = fread(x) %>% 
            column_to_rownames("GENE") %>% 
            data.matrix %>% 
            CreateSeuratObject() %>% 
            NormalizeData()

        hUSI = GetAssayData(scdat) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}
        data <- list(method=data.frame(calc_scores(Methods,scdat,'method'),'hUSI'=hUSI),
                        marker=data.frame(calc_scores(SenMarkers,scdat,'marker'),'hUSI'=hUSI),
                        ssgsea=data.frame(calc_scores(EnrichSet,scdat,'ssgsea'),'hUSI'=hUSI))
        return(data)
    }) %>% 
    set_names( list.files("Tang2019") %>% gsub(".*_", "", .) %>% gsub("\\..*", "", .) )

Results$Tang2019 = Tang2019List

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
hUSI = GetAssayData(Aarts2017) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}

dat_method = calc_scores(Methods,Aarts2017,'method')
dat_marker = calc_scores(SenMarkers,Aarts2017,'marker')
dat_ssgsea = calc_scores(EnrichSet,Aarts2017,'ssgsea')

Aarts2017List = list(method=dat_method,marker=dat_marker,ssgsea=dat_ssgsea)

Results$Aarts2017 = Aarts2017List
Results$Aarts2017$hUSI = hUSI

### Georgilis2018 GSE101766 IMR90 OIS,full-length,TPM
Georgilis2018 = CreateSeuratObject(
    counts = fread("Georgilis2018/valid_TPM_dataset.tsv") %>% column_to_rownames("Gene Name") %>% data.matrix, 
    meta.data = fread("Georgilis2018/filereport_read_run_PRJNA395386_tsv.txt", header = T) %>% column_to_rownames("sample_title")) %>% 
    NormalizeData()
Georgilis2018 = subset(Georgilis2018,read_count>1e6)
Georgilis2018$Condition = case_when(grepl("Allstars|Water|Hs|Extras[5-8]", colnames(Georgilis2018)) ~ "Senescence", 
                                    grepl("Noninduced|Extras[1-4]", colnames(Georgilis2018)) ~ "Growing", TRUE ~ "other")
Georgilis2018 = subset(Georgilis2018,Condition!='other')
Georgilis2018$Condition = factor(Georgilis2018$Condition)

hUSI = GetAssayData(Georgilis2018) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}

dat_method = calc_scores(Methods,Georgilis2018,'method')
dat_marker = calc_scores(SenMarkers,Georgilis2018,'marker')
dat_ssgsea = calc_scores(EnrichSet,Georgilis2018,'ssgsea')

Georgilis2018List = list(method=dat_method,marker=dat_marker,ssgsea=dat_ssgsea)

Results$Georgilis2018 = Georgilis2018List
Results$Georgilis2018$hUSI = hUSI


### AUC
AUClist = list()
for(dataset in c('Teo2019',"Tang2019",'Aarts2017','Georgilis2018')){
    for(type in c('method','marker','ssgsea')){
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
            colnames(dat) = gsub("DAS.mSS", "DAS+mSS", colnames(dat))
            dat$condition = as.factor(ifelse(dat$sample %in% c("senescence","LowPD50Gy"),'Senescence','Growing'))

            auc <- calc_auc(dat,type,scorelist)
        }
        else if (dataset %in% c('Aarts2017','Georgilis2018')) {
            dat = Results[[dataset]][[type]] %>% 
                    mutate(condition=Condition) %>% 
                    mutate(hUSI=Results[[dataset]][['hUSI']])
            dat = dat[,complete.cases(t(dat))]

            auc <- calc_auc(dat,type,scorelist)
        }
        AUClist[[dataset]][[type]] = auc
    }
}

### write auc results
### method
df_plot = do.call(rbind,lapply(AUClist,function(x) {data.frame(t(x[['method']]))}))
### marker
comm_markers = c("hUSI","GLB1","TP53","CDKN1A","CDKN2A","LMNB1","RB1","CDK1","CDK4","CDK6","MKI67","CDKN2B")
df_plot = do.call(rbind,lapply(AUClist,function(x) {data.frame(t(x[['marker']][comm_markers,]))}))
### ssGSEA
df_plot = do.call(rbind,lapply(AUClist,function(x) {data.frame(t(x[['ssgsea']]))}))

#################################### comparision in GTEx and TCGA #############################
setwd('/home/wangjing/wangj/AgingScore/Data/Bulk_TrainModel/CS_score_of_single_cell_datasets/')
source("/mnt/data1/wangj/codebase/HUSI/transID.R")
##### GTEx bulk data
### load cs score
TCSER_gtex = read.csv('Senescence_quantification_GTEX.csv')
dim(TCSER_gtex)
TCSER_gtex = TCSER_gtex[!duplicated(TCSER_gtex$Sample),]
rownames(TCSER_gtex) = TCSER_gtex$Sample
library(CePa)
GTEx_exp <- read.gct("/home/wangjing/wangj/AgingScore/Data/Bulk_TrainModel/CS_score_of_single_cell_datasets/bulk-gex_v8_rna-seq_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
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

### marker genes
for(sl in SenMarkers){
  if(sl %in% rownames(Counts_gtex)){
    score = Counts_gtex[sl,] 
    names(score) = colnames(Counts_gtex)
    scoreList_GTEx[[sl]] = score
  }
}

### Methods
for(sl in names(Methods)){
    if(sl == 'laf'){
        for (s in names(Methods[[sl]])) {
            score = calc_laf(Counts_gtex %>% as.matrix, Methods[[sl]][[s]])
            names(score) = colnames(Counts_gtex)
            scoreList_GTEx[[s]] = score
            }
        }
        if(sl == 'lassoCS'){
            genes = intersect(rownames(Counts_gtex),names(lassoCS))
            mat = Counts_gtex[genes,] 
            mat = apply(mat,2,function(x) {x * lassoCS[rownames(mat)]})
            tryCatch({score = colSums(mat) %>% data.frame 
            rownames(score) = colnames(Counts_gtex)}, 
            error = function(e) {
            score = mat %>% data.frame})
            rownames(score) = colnames(Counts_gtex)
            scoreList_GTEx[[sl]] = score   
        }
        if(sl == 'CSS'){
            genes = intersect(rownames(Counts_gtex),names(CSS))
            mat = Counts_gtex[genes,] 
            mat = apply(mat,2,function(x) {x * CSS[rownames(mat)]}) 
            tryCatch({score = colSums(mat) %>% data.frame 
            rownames(score) = colnames(Counts_gtex)}, 
            error = function(e) {
            score = mat %>% data.frame})
            rownames(score) = colnames(Counts_gtex)
            scoreList_GTEx[[sl]] = score
    }

}

cormat_gtex <- do.call(data.frame, scoreList_GTEx) %>% cor(method = 'spearman')
###### TCGA bulk data
### load cs score
TCSER_tcga = read.csv('Senescence_quantification_TCGA.csv')
dim(TCSER_tcga)
TCSER_tcga = TCSER_tcga[!duplicated(TCSER_tcga$sample),]
rownames(TCSER_tcga) = TCSER_tcga$sample

TCGA_exp = fread('/mnt/data1/wangj/MyProject/Data/TCGA_pan_tpm/tcga_RSEM_gene_tpm.gz') 
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
