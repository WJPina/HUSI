library(Seurat)
setwd("~/wangj/AgingScore/BulkData/")
# # load model
mm_l2 = readRDS("/home/wangjing/wangj/AgingScore/BulkData/Bulk_TrainModel/mm_l2.rds")
library(dplyr)
library(tibble)
library(data.table)
################################ Pre-preocess raw trian data ###############################
setwd("~/wangj/AgingScore/Data1/Bulk_TrainModel")
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
# mean center and split to train and background dataset
X_centre = X - (apply(X, 1, mean))
idx_senescent = colnames(X_centre) %in% filter(metadata, condition %in% "senescent")$sample_title
X_tr = X_centre[,idx_senescent]
X_bk = X_centre[,!idx_senescent]
# training model
mm_l2 = gelnet( t(X_tr), NULL, 0, 1 )
mm_l1 = gelnet( t(X_tr), NULL, 0.1, 0 )
saveRDS(mm_l2,file="mm_l2.rds")
saveRDS(mm_l1,file="mm_l1.rds")
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
################### Valid model in  RNA-seq and microarray #######################
library(annaffy)
library(magrittr)
library(ggplot2)
library(reshape2)
library(GEOquery)
library(tibble)
library(data.table)
library(dplyr)

### RNA-seq
setwd("~/wangj/AgingScore/Data1/Bulk_RNA-seq/")
# GSE60340
s_IS = fread("GSE60340_induced_sene_TPM.csv") %>% 
    column_to_rownames("Gene Name") %>% 
    data.matrix %>% 
    {
        apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )
    }
# GSE130306
library(magrittr)
s_RS = list.files("GSE130306/", "GSM", full.names = T) %>% 
    sapply(function(x) {
        fread(x) %>% 
        filter( !duplicated(gene_name) ) %>% 
        column_to_rownames("gene_name") %>% 
        data.matrix %>% 
        {
        apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )
        }
    }) %>% 
    set_names( list.files("GSE130306/", "GSM", full.names = T) %>% gsub(".*RNAseq_", "", .) %>% gsub(".txt.gz", "", .) )


### Microarray
setwd("~/wangj/AgingScore/BulkData/Bulk_Microarray/")
# load raw data
fs = paste(c('GSE19864','GSE16058','GSE83922','GSE11954','GSE100014','GSE77239'),"eSet.Rdata",sep = "_")
ArrayList <- sapply(fs, function(x) mget(load(x)), simplify = TRUE) 
names(ArrayList) <- unlist(lapply(strsplit(names(ArrayList),"_"),"[",1))
# calculate hSI
GPL = c("GPL570" = "hgu133plus2.db","GPL3921" = "hthgu133a.db","GPL11532" = "hugene11sttranscriptcluster.db")
ScoreList = list()
for(i in 1:length(ArrayList)){
  gene_id = aafSymbol( rownames(ArrayList[[i]][[1]]), GPL[ArrayList[[i]][[1]]@annotation]) %>% as.character
  ScoreList[[i]] <- ArrayList[[i]][[1]] %>% 
                    exprs %>% 
                    set_rownames( gene_id ) %>% 
                    melt(.) %>% 
                    filter(!Var1 %in% "character(0)") %>% 
                    as.data.table %>% 
                    dcast.data.table( Var1~Var2, fun.aggregate = mean ) %>% 
                    column_to_rownames("Var1") %>% 
                    data.matrix %>% 
                    { apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )})
                    }
}
names(ScoreList) <- paste("s",gsub("GSE","",names(ArrayList)),sep = "_")
sapply(ScoreList,length)
################################# Batch Effect #################################
data.matrix <- read.table("/home/wangjing/wangj/AgingScore/Data1/Bulk_BatchEffect/batch_IMR90_4OHT.tsv",sep = "\t",header = T,row.names = 1)
s_batch <- data.matrix %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}
########################### model weight GSEA ###################################
library(fgsea)
library(clusterProfiler)

library(msigdbr)
set.seed(233)
msigdb_all = msigdbr()
pathways_sene = msigdb_all %>% 
                  filter( (gs_cat %in% "H" ) ) %>% 
                  mutate( gs_name=gsub("HALLMARK_", "", gs_name) ) %>% 
                  plyr::dlply(.variables = "gs_name", .fun = function(x) x$gene_symbol )
res_fgsea_sene = fgseaMultilevel(pathways = pathways_sene, stats = sort(mm_l2$w,decreasing=T), nPermSimple = 10000)

hallmarker = read.gmt("/mnt/data3/wangj2/GeneSets/h.all.v2023.1.Hs.symbols.gmt")
fgsea <- GSEA(sort(mm_l2$w,decreasing=T), 
                exponent = 1, 
                minGSSize = 10,
                maxGSSize = 500, 
                pvalueCutoff = 0.05, 
                pAdjustMethod = "BH", 
                TERM2GENE=hallmarker,
                seed = 233,
                by = "fgsea",
                eps=0)

df = fgsea@result
df = filter(df,pvalue < 0.05)
df = df[order(df$NES,decreasing=T),]
df[,c('ID','NES')]
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

############################### application on melanoma ########################
### load melanoma data
setwd('/home/wangjing/wangj/AgingScore/Data1/scRNA_melanome')
melanoma <- read.table("GSE72056_melanoma_single_cell_revised_v2.txt",sep="\t",header=TRUE)

### load meta data
sample <- read.csv("sample.csv",row.names=1,stringsAsFactors=FALSE)
metadata <- t(melanoma[c(2,3),-1])
metadata <- metadata[metadata[,1]!=0,]
metadata.sub <- metadata[metadata[,2]==0,]
metadata.sub <- metadata.sub[metadata.sub[,1]==2,]
metadata.sub2 <- metadata[metadata[,2]!=0,]
metadata.sub2 <- metadata.sub2[metadata.sub2[,1]!=2,]
metadata <- rbind(metadata.sub,metadata.sub2)
celltype <- c(rep("Tumor",nrow(metadata.sub)),metadata.sub2[,2])
celltype[celltype=="1"] <- "T cell"
celltype[celltype=="2"] <- "B cell"
celltype[celltype=="3"] <- "Macro. cell"
celltype[celltype=="4"] <- "Endo. cell"
celltype[celltype=="5"] <- "CAF cell"
celltype[celltype=="6"] <- "NK cell"

names(celltype) <- c(rownames(metadata.sub),rownames(metadata.sub2))
melanoma <- melanoma[!duplicated(melanoma[,1]),]
label <- melanoma[,1]
melanoma <- melanoma[-c(1:3),c(rownames(metadata.sub),rownames(metadata.sub2))]
rownames(melanoma) <- label[-c(1:3)]
library(Seurat)
melanoma_obj <- CreateSeuratObject(melanoma)

### calculate aging score
mm_l2 = readRDS('/home/wangjing/wangj/AgingScore/AgingScorePro/l2_model_add_dermal.rds')
library(dplyr)
AgeScore  = melanoma_obj@assays$RNA@data[] %>% {apply( ., 2, function(z) {cor(z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )})}

melanoma_obj$hSI <- AgeScore[colnames(melanoma_obj)]
melanoma_obj$celltype <- celltype

melanoma_obj <- FindVariableFeatures(object = melanoma_obj, selection.method = 'vst', nfeatures =1500)
melanoma_obj <- ScaleData(melanoma_obj,features=rownames(melanoma_obj))
melanoma_obj <- RunPCA(melanoma_obj,npcs = 30,verbose = F,features=rownames(melanoma_obj)) 
# ElbowPlot(melanoma_obj,ndims=30,reduction="pca")
melanoma_obj <- RunTSNE(melanoma_obj, reduction = "pca", dims = 1:20) 
DimPlot(melanoma_obj, reduction = 'tsne', group.by = 'celltype',label=1)

### extract tumor cells
EpiExp.m <- melanoma_obj[,melanoma_obj$celltype=="Tumor"]

### ranking the genes based on each gene expression correlation with the aging score
corAging <- function(x,agingScore){cor <- cor(x,agingScore);cor}
cor.genes <- apply(GetAssayData(EpiExp.m),1,corAging,agingScore=EpiExp.m$hSI)
cor.genes[is.na(cor.genes)] <- 0
features <- rownames(EpiExp.m)[order(abs(cor.genes),decreasing=TRUE)[1:1500]]
features 

### divided tumor cells by gaussion 
library(plyr)
library(mclust)
gaussian = EpiExp.m@meta.data$hSI %>% {log2((1+ .)/(1- .))} %>% Mclust(G = 3)
EpiExp.m$age_class <- gaussian$classification
DimPlot(EpiExp.m, reduction = 'tsne', group.by = 'age_class',label=1)

###Running the ICAnet
library(ICAnet)
source('/home/wangjing/wangj/AgingScore/AgingScorePro/Data1_Scripts/getPPI_String.R')

EpiExp.matrix <- as.matrix(GetAssayData(EpiExp.m,))[features,]
EpiExp.m[['AgingExp']] <- CreateAssayObject(EpiExp.matrix)
DefaultAssay(EpiExp.m) <- "AgingExp"

EpiExp.m$batch <- rep("batch",ncol(EpiExp.m))
Ica.epi <- ICAcomputing(EpiExp.m,ICA.type="JADE",RMT=TRUE,two.stage=FALSE)
Ica.filter <- CrossBatchGrouping(Ica.epi$ica.pooling)
# PPI <- readRDS("PPI_feature.rds")
PPI <- getPPI_String(EpiExp.m,species=9606)

EpiExp.m <- RunICAnet(EpiExp.m,Ica.epi$ica.pooling,PPI.net = PPI,scale=FALSE,ModuleSignificance = FALSE,cores = 1,aucMaxRank=1000)
dev.off()

### Using PHATE to predict the transition process of tumor aging
library(phateR)
library(reticulate)
use_python("/home/tools/anaconda3/envs/sc/bin/python3", required = T)

Epi.data <- t((GetAssayData(EpiExp.m)))
Epi.data <- as.data.frame(Epi.data)
Epi.phate <- phate(Epi.data)
branch <- EpiExp.m$age_class
EpiExp.m[['phate']] <- CreateDimReducObject(Epi.phate$embedding,key="phate_")

### markers in microarray
library(dplyr)
library(annaffy)
library(magrittr)
library(GEOquery)
load('/home/wangjing/wangj/AgingScore/AgingScorePro/Data1_Scripts/ModelValidData.RData')
GPL = c("GPL570" = "hgu133plus2.db","GPL3921" = "hthgu133a.db","GPL11532" = "hugene11sttranscriptcluster.db")
gene_id = aafSymbol( rownames(ArrayList$GSE83922[[1]]), GPL[ArrayList$GSE83922[[1]]@annotation]) %>% as.character
exp = ArrayList$GSE83922[[1]] %>% exprs %>% set_rownames( gene_id )
exp = exp[rownames(exp)!="character(0)",]
head(exp)

meta <- ArrayList [["GSE83922"]][[1]] %>% 
        pData %>% 
        mutate(condition= `cell phenotype:ch1`) %>% 
        .[which(.$condition %in% c("young","senescent")),] %>% 
        mutate(condition = factor(condition,levels = c("young","senescent"),ordered = T)) 

marker_set = list()
for(i in unique(EpiExp.m$age_class)){
  marker_set[[i]] <- FindMarkers(EpiExp.m,ident.1 = i,group.by="age_class",only.pos = TRUE,assay = "RNA")
  marker_set[[i]] <- filter(marker_set[[i]],p_val_adj<0.05)
}
lapply(marker_set,dim)

gene_set <- lapply(marker_set,function(x){x <- rownames(x)[rownames(x)%in%rownames(exp)]})

mat <- exp[unlist(gene_set),rownames(meta)]
mat=mat[!apply(mat, 1, sd)==0,]
mat=Matrix::t(scale(Matrix::t(mat),center=TRUE))
mat=mat[is.na(row.names(mat)) == FALSE,]
mat[is.nan(mat)] = 0
mat[mat>2] = 2
mat[mat< -2] = -2




