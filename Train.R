library(Seurat)
library(dplyr)
library(tibble)
library(data.table)
################################ Pre-preocess raw trian data ###############################
setwd("~/wangj/AgingScore/Data/Bulk_TrainModel")
raw_TPM = cbind(
    fread("dermal_for_training_TPM.csv") %>% column_to_rownames("Gene Name") %>% data.matrix, 
    fread("senescence_raw_TPM_for_training.tsv") %>% column_to_rownames("Gene Name") %>% as.matrix
)
gtf = rtracklayer::import("gencode.v31.annotation.gtf")
# Remove non protein coding genes
idx = rownames(raw_TPM) %in% as.data.table(gtf)[type %in% "gene" & gene_type %in% "protein_coding"]$gene_name
exprData = raw_TPM[idx,]
sprintf('RNA-Seq was filtered from %d features to %d protein coding genes', nrow(raw_TPM), nrow(exprData))
# Filter out genes with close to zero expression
idx = apply(exprData, 1, function(v) sum(v<=3) )/ncol(exprData) > 0.99
sprintf('There were %d genes with 99%% of samples having TPM < 3', sum(idx))
sprintf(' Of these, %d genes the mean expression was %0.2g and median %0.2g', sum(idx), mean(exprData[idx,]), median(exprData[idx,]))
exprData = exprData[!idx,]
dim(exprData)
# Log-transform data and standardize
X = log2(exprData + .00001) %>% scale()
############################## Pre-preocess raw bulk meata data ##############################
meta_list = list()
## samples over 10 days are labeled senescent
meta_list$CurrBio17 = fread("E-MTAB-5403.sdrf.txt", select = c("Source Name", "Factor Value[cell type]", "Factor Value[irradiate]")) %>% unique %>% 
    mutate(condition = if_else( grepl("day [12]0", `Factor Value[irradiate]`), "senescent", "other" ) ) %>% 
    dplyr::rename(sample_title = `Source Name`) %>% mutate( study_accession = "E-MTAB-5403" ) %>% dplyr::select(study_accession, sample_title, condition)

meta_list$Alspach_GSE56293 = fread("PRJNA243087.txt", select = c("study_accession", "sample_title")) %>% 
    mutate( condition = if_else( grepl("late", sample_title), "senescent", "other" ) )

meta_list$Herranz_GSE61130 = fread("PRJNA260300.txt", select = c("study_accession", "sample_title")) %>% 
    mutate( sample_title = gsub(" ", "_", sample_title) ) %>% 
    mutate( condition = if_else( grepl("EV_-_", sample_title), "other", "senescent" ) )

meta_list$Marthandan_GSE63577 = plyr::ldply(c("PRJNA259739.txt", "PRJNA268292.txt"), .fun = fread, select = c("study_accession", "sample_title") ) %>% 
    filter(sample_title %in% colnames(X)) %>% 
    mutate( condition = rep( c(rep("other", 3), rep("senescent", 3)), 5) )

meta_list$Marthandan_GSE64553 = fread("PRJNA271291.txt", select = c("study_accession", "sample_title")) %>% 
    filter(sample_title %in% colnames(X)) %>% 
    mutate( condition = if_else( grepl("22", sample_title), "other", "senescent" ) )

meta_list$Rai_GSE53356 = fread("PRJNA231833.txt", select = c("study_accession", "sample_title")) %>% 
    mutate( condition = rep( c(rep("other", 2), rep("senescent", 2)), 1) )

meta_list$Crowe_GSE58910 = fread("PRJNA253866.txt", select = c("study_accession", "sample_title")) %>% 
    mutate( condition = rep( c(rep("other", 2), rep("senescent", 2)), 1) )

meta_list$Casella_GSE130727 = fread("PRJNA541183.txt", select = c("study_accession", "sample_title")) %>% 
    mutate( sample_title = gsub(" ", "_", sample_title) %>% gsub("\\(", "_", .) %>% gsub("\\)", "", .) ) %>% 
    mutate( condition = if_else(
        grepl("\\+", sample_title) | grepl("PDL52", sample_title) | grepl("HRAS", sample_title) | grepl("Replicative_Exhaustion", sample_title), 
        "senescent", "other"
    ) )
# aggregate to a dataframe
metadata_P = do.call(rbind, meta_list) %>% rownames_to_column("study_name") %>% mutate( study_name = gsub("\\..*", "", study_name) )
metadata_P
# aggregate dermal meta
meta_GB_2018 = data.frame(
    sample_title=fread("dermal_allAges_TPM.csv") %>% colnames %>% .[-1]
) %>% 
    tidyr::separate(sample_title, into=c("id", "age_raw", "gender_raw", "race"), remove=F, sep="_") %>% 
    mutate(disease=if_else(grepl("mos", age_raw) | grepl("HGPS", gender_raw), "HGPS", "normal") %>% factor(levels = c("normal", "HGPS"))) %>% 
    mutate(gender=if_else( grepl("^M", gender_raw, ignore.case = T) , "M", "F") %>% factor(levels = c("M", "F"))) %>% 
    mutate(age=gsub("[Yy].*", "", age_raw) %>% as.numeric) %>% 
    
    filter( disease == "normal" ) %>% 
    mutate( condition=case_when(
        age>80 ~ "senescent", 
        age<25 ~ "other"
    ) ) %>% 
    tidyr::drop_na(condition) %>% 
    mutate( study_accession="PRJNA454681" ) %>% 
    dplyr::select( study_accession, sample_title, condition )
meta_GB_2018$study_name <- rownames(meta_GB_2018)
meta_GB_2018
# aggregate all
metadata <- rbind(metadata_P,meta_GB_2018)
metadata$study_name <- rownames(metadata)
metadata
save(X,metadata,file = 'ModelTrainData.RData')
############################# Train model ######################################
load("ModelTrainData.RData")
### mean center and split to train and background dataset
X_centre = X - (apply(X, 1, mean))
idx_senescent = colnames(X_centre) %in% filter(metadata, condition %in% "senescent")$sample_title
X_tr = X_centre[,idx_senescent]
X_bk = X_centre[,!idx_senescent]
# training model
mm_l2 = gelnet( t(X_tr), NULL, 0, 1 )
# mm_l1 = gelnet( t(X_tr), NULL, 0.1, 0 )
saveRDS(mm_l2,file="mm_l2.rds")
# saveRDS(mm_l1,file="mm_l1.rds")

### Leave-one-out cross validation
auc <- c()
for(i in 1:ncol(X_tr)){
  ## Train a model on non-left-out data
  X1 <- X_tr[,-i]
  m1 <- gelnet( t(X1), NULL, 0, 1 )
  ## Score the left-out sample against the background
  s_bk <- apply( X_bk, 2, function(z) {cor( m1$w, z, method="sp" )} )
  s1 <- cor( m1$w, X_tr[,i], method="sp" )
  ## AUC = P( left-out sample is scored above the background )
  auc[i] <- sum( s1 > s_bk ) / length(s_bk)  
  cat( "Current AUC: ", auc[i], "\n" )
  cat( "Average AUC: ", mean(auc), "\n" )
}




