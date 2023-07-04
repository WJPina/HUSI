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

calc_scores <- function(scorelist=list(),object=NULL){
    data = object@meta.data
    for(sl in names(scorelist)){
        if(sl %in% c("DAS","mSS","DAS_mSS")){
            score = calc_laf(GetAssayData(object) %>% expm1 %>% as.matrix, scorelist[[sl]])
            data[[sl]] = score
        }
        if(sl == 'Sig.RS'){
            score = GSVA::gsva(GetAssayData(object) %>% expm1 %>% as.matrix, 
                            list(signature1311=scorelist[[sl]], 
                                signature1311_transform=AnnotationDbi::select(org.Hs.eg.db, scorelist[[sl]], columns = "SYMBOL", keytype = "ENSEMBL")$SYMBOL), 
                            method="ssgsea", ssgsea.norm = TRUE, verbose = TRUE)[1,]
            data[[sl]] = score
        }
        if(sl == 'SenMarkers'){
            score = object[scorelist[[sl]],] %>% GetAssayData %>% expm1 %>% as.matrix %>% t %>% data.frame 
            rownames(score) = colnames(object)
            data = cbind(data,score)
        }
    }
    return(data)
} 

calc_auc <- function(data=data.frame(),type = '',scorelist){
    fold = createMultiFolds(data$condition, k = 10, times = 3)
    auc = sapply(fold, function(sampling) {
        df = data[sampling,]
        if(type=="method"){
            c("DAS","mSS",'Sig.RS','DAS_mSS','hSI') %>% 
            sapply(function(idx) {pROC::roc(df$condition, df[,idx], levels=c("Growing", "Senescence")) %>% pROC::auc()})
        }
        else {
            c('hSI',scorelist[['SenMarkers']][scorelist[['SenMarkers']] %in% colnames(data)]) %>% 
            sapply(function(idx) {pROC::roc(df$condition, df[,idx], levels=c("Growing", "Senescence")) %>% pROC::auc()})
        }
    })
    return(auc)
}
### score list
laf_DAS = fread("LaffertyWhyte2010_DAS.csv")
laf_mSS = fread("LaffertyWhyte2010_mSS.csv")
rep_sene_genes = readxl::read_excel("SigRS.xls", sheet = 5, skip = 2)$Ensembl_ID
SenMarkers = c("GLB1", "TP53", "CDKN1A", "CDKN2A", "LMNB1", "IL1A", "RB1", "CDK1", "CDK4", "MKI67", "CDKN2B")

ScoreList = list(DAS = laf_DAS,mSS=laf_mSS,DAS_mSS=rbind(laf_DAS, laf_mSS),Sig.RS=rep_sene_genes,SenMarkers=SenMarkers)
### Teo2019
Teo2019 = CreateSeuratObject(
    fread("Teo2019/GSE115301_Growing_Sen_10x_count.txt.gz") %>% column_to_rownames("V1") %>% data.matrix, 
    meta.data = fread("Teo2019/GSE115301_Growing_Sen_10x_metadata.txt.gz", header = T) %>% column_to_rownames("V1")
)
Teo2019 = NormalizeData(Teo2019, scale.factor = 1e4)
hSI = GetAssayData(Teo2019) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}

dat = calc_scores(ScoreList,Teo2019)
dat_Teo = dat %>% 
        mutate( condition=gsub("[0-9]$", "", Condition1) ) %>% 
        mutate(hSI=hSI)
dat_Teo = dat_Teo[,complete.cases(t(dat_Teo))]
auc_Teo <- calc_auc(dat_Teo,'method',ScoreList)

### Georgilis2018
Georgilis2018 = CreateSeuratObject(
    fread("Georgilis2018/valid_TPM_dataset.tsv") %>% column_to_rownames("Gene Name") %>% data.matrix, 
    meta.data = fread("Georgilis2018/filereport_read_run_PRJNA395386_tsv.txt", header = T) %>% column_to_rownames("sample_title")
)
Georgilis2018 = subset(Georgilis2018,read_count>1e6)
Georgilis2018 = NormalizeData(Georgilis2018, scale.factor = 1e4)
hSI = GetAssayData(Georgilis2018) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}

dat = calc_scores(ScoreList,Georgilis2018)
dat_Geo = dat %>% 
        rownames_to_column('sample') %>%
        mutate(condition=case_when(
        grepl("Allstars|Water|Hs|Extras[5-8]", sample) ~ "Senescence", 
        grepl("Noninduced|Extras[1-4]", sample) ~ "Growing", 
        TRUE ~ "other")) %>% 
        mutate(hSI=hSI) %>%
        filter( !condition %in% "other" ) %>% 
        mutate( condition=as.factor(condition))
dat_Geo = dat_Geo[,complete.cases(t(dat_Geo))]
auc_Geo <- calc_auc(dat_Geo,'marker',ScoreList)

###Tang2019
Tang2019List = list.files("Tang2019", full.names = T) %>% 
    lapply(function(x) {
        scdat = fread(x) %>% 
            column_to_rownames("GENE") %>% 
            data.matrix %>% 
            CreateSeuratObject() %>% 
            NormalizeData(scale.factor=1e4)

        hSI = GetAssayData(scdat) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}
        data <- calc_scores(ScoreList,scdat)
        data = data.frame(data,'hSI'=hSI)
        return(data)
    }) %>% 
    set_names( list.files("Tang2019") %>% gsub(".*_", "", .) %>% gsub("\\..*", "", .) )

dat_Tang = data.table::rbindlist(Tang2019List,use.names=TRUE, fill=TRUE, idcol="sample") %>% data.frame
dat_Tang = dat_Tang[,complete.cases(t(dat_Tang))]
dat_Tang$condition = as.factor(ifelse(dat_Tang$sample %in% c("senescence","LowPD50Gy"),'Senescence','Growing'))
auc_Tang <- calc_auc(dat_Tang,'method',ScoreList)




