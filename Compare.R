mm_l2 = readRDS("/home/wangjing/wangj/AgingScore/BulkData/Bulk_TrainModel/mm_l2.rds")
############################## Comparision ######################################
setwd("/home/wangjing/wangj/AgingScore/Comparison")

library(Seurat)
library(data.table)
library(tibble)
library(sparseMatrixStats)
library(org.Hs.eg.db)
library(dplyr)
library(caret)
library(magrittr)

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

calc_scores <- function(scorelist=NULL,object=NULL,type='laf'){
    data = object@meta.data
    if(type == 'laf'){
        for(sl in names(scorelist)){
                score = calc_laf(GetAssayData(object) %>% expm1 %>% as.matrix, scorelist[[sl]])
                data[[sl]] = score
        }
    }
    if(type == 'ssgsea'){
        for(sl in names(scorelist)){
            score = GSVA::gsva(GetAssayData(object) %>% expm1 %>% as.matrix, list(signature=scorelist[[sl]]), method="ssgsea", ssgsea.norm = TRUE, verbose = TRUE)[1,]
            data[[sl]] = score
        }
    }
    if(type == 'marker'){
        score = object[scorelist,] %>% GetAssayData %>% expm1 %>% as.matrix %>% t %>% data.frame 
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
            c("DAS","mSS",'DAS_mSS','hSI') %>% 
            sapply(function(idx) {pROC::roc(df$condition, df[,idx], levels=c("Growing", "Senescence")) %>% pROC::auc()})
        }
        else if (type=='marker'){
            c('hSI',scorelist[scorelist %in% colnames(data)]) %>% 
            sapply(function(idx) {pROC::roc(df$condition, df[,idx], levels=c("Growing", "Senescence")) %>% pROC::auc()})
        }
        else{
            c('hSI',names(scorelist)) %>% 
            sapply(function(idx) {pROC::roc(df$condition, df[,idx], levels=c("Growing", "Senescence")) %>% pROC::auc()})
        }
    })
    return(auc)
}
### score list
laf_DAS = fread("LaffertyWhyte2010_DAS.csv")
laf_mSS = fread("LaffertyWhyte2010_mSS.csv")
lafSet = list(DAS = laf_DAS,mSS=laf_mSS,DAS_mSS=rbind(laf_DAS, laf_mSS))

SenMarkers = c("GLB1", "TP53", "CDKN1A", "CDKN2A", "LMNB1", "IL1A", "RB1", "CDK1", "CDK4", "MKI67", "CDKN2B")

EnrichSet=cogena::gmt2list("/mnt/data2/zhouxl/Pan_Cancer/Data/Signatures/gene_50signatures_merge.gmt")
EnrichSet=EnrichSet[43:50]
rep_sene_genes = readxl::read_excel("SigRS.xls", sheet = 5, skip = 2)$Gene_names
EnrichSet$Sig.RS = rep_sene_genes

### Teo2019
Teo2019 = CreateSeuratObject(
    fread("Teo2019/GSE115301_Growing_Sen_10x_count.txt.gz") %>% column_to_rownames("V1") %>% data.matrix, 
    meta.data = fread("Teo2019/GSE115301_Growing_Sen_10x_metadata.txt.gz", header = T) %>% column_to_rownames("V1")
)
Teo2019 = NormalizeData(Teo2019, scale.factor = 1e4)
hSI = GetAssayData(Teo2019) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}

dat_laf = calc_scores(lafSet,Teo2019,'laf')
dat_marker = calc_scores(SenMarkers,Teo2019,'marker')
dat_ssgsea = calc_scores(EnrichSet,Teo2019,'ssgsea')

Teo2019List = list(laf=dat_laf,marker=dat_marker,ssgase=dat_ssgsea)

dat_Teo = dat_ssgsea %>% 
            mutate( condition=gsub("[0-9]$", "", Condition1) ) %>% 
            mutate(hSI=hSI)
dat_Teo = dat_Teo[,complete.cases(t(dat_Teo))]

# auc <- calc_auc(dat_Teo,'laf',lafSet)
# auc <- calc_auc(dat_Teo,'marker',SenMarkers)
auc <- calc_auc(dat_Teo,'ssgsea',EnrichSet)

###Tang2019
Tang2019List = list.files("Tang2019", full.names = T) %>% 
    lapply(function(x) {
        scdat = fread(x) %>% 
            column_to_rownames("GENE") %>% 
            data.matrix %>% 
            CreateSeuratObject() %>% 
            NormalizeData(scale.factor=1e4)

        hSI = GetAssayData(scdat) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}
        data <- list(laf=data.frame(calc_scores(lafSet,scdat,'laf'),'hSI'=hSI),
                        marker=data.frame(calc_scores(SenMarkers,scdat,'marker'),'hSI'=hSI),
                        ssgsea=data.frame(calc_scores(EnrichSet,scdat,'ssgsea'),'hSI'=hSI))
        return(data)
    }) %>% 
    set_names( list.files("Tang2019") %>% gsub(".*_", "", .) %>% gsub("\\..*", "", .) )

type = 'ssgsea'
dat_Tang = data.table::rbindlist(lapply(Tang2019List,function(x){x[[type]]}),use.names=TRUE, fill=TRUE, idcol="sample") %>% data.frame
dat_Tang = dat_Tang[,complete.cases(t(dat_Tang))]
dat_Tang$condition = as.factor(ifelse(dat_Tang$sample %in% c("senescence","LowPD50Gy"),'Senescence','Growing'))

# auc <- calc_auc(dat_Tang,'laf',lafSet)
# auc <- calc_auc(dat_Tang,'marker',SenMarkers)
auc <- calc_auc(dat_Tang,'ssgsea',EnrichSet)




