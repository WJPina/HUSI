##### This scripts is for model training and evaluation
#################################### library ################################### 
library(data.table)
library(dplyr)
library(gelnet)
library(openxlsx)
library(magrittr)
library(org.Hs.eg.db)
library(GEOquery)
library(limma)
library(affy)
library(DESeq2)
library(dplyr)

source('R/functions_new.R')
source('R/functions.R')

setwd('HUSI/')
################### preprocess collected bulk RNA-seq data #####################
########## read collected bulk samples
raw_TPM = read.csv('Data/GEP_TPM_hUSI.csv') ### uploded to figueshare
dim(raw_TPM)
raw_TPM[1:5,1:5]

gene_id = raw_TPM$gene_id
gene_id = unlist(lapply(strsplit(gene_id,'\\|'),'[',2))
raw_TPM=raw_TPM[,-1]
raw_TPM=raw_TPM[!duplicated(gene_id),]
gene_id=gene_id[!duplicated(gene_id)]
rownames(raw_TPM) = gene_id

raw_meta = read.xlsx('Data/trainSet_meta843.csv')
rownames(raw_meta) <- raw_meta$Data_id
dim(raw_meta)
raw_meta[1:5,1:5]

########## filter samples
meta = raw_meta[colnames(raw_TPM),]
table(meta$Stage)
meta$evidence_n=rowSums(as.matrix(meta[,c(14:24)]))

meta$label = ifelse(meta$Stage == 'Sene'& meta$evidence_n>=2,'senescent','non-senescent')
table(meta$label)
dat = table(meta$Cell_type_id,meta$label) %>% as.data.frame.array()
dat = dat[apply(dat!= 0 , 1 , all),]

meta = filter(meta,Cell_type_id %in% rownames(dat))
table(meta$label)
meta[1:5,1:5]

exprData = raw_TPM[,rownames(meta)] %>% as.matrix()
dim(exprData)
exprData[1:5,1:5]

###### Remove non protein coding genes
gtf = rtracklayer::import("gencode.v31.annotation.gtf") ### download from gencode
idx = rownames(exprData) %in% as.data.table(gtf)[type %in% "gene" & gene_type %in% "protein_coding"]$gene_name
exprData = exprData[idx,]

###### remove ribosomal and mitchodrial genes
mtGenes = grep("^MT-", rownames(exprData))
riboGenes = grep("^RP[SL]", rownames(exprData))
exprData = exprData[-c(mtGenes, riboGenes),]
sprintf('genes was filtered from %d to %d', nrow(raw_TPM), nrow(exprData))

###### Filter out genes with close to zero expression
idx = apply(exprData, 1, function(v) sum(v<=3) )/ncol(exprData) > 0.99
sprintf('There were %d genes with 99%% of samples having TPM < 3', sum(idx))
sprintf(' Of these, %d genes the mean expression was %0.2g and median %0.2g', 
        sum(idx), mean(exprData[idx,]), median(exprData[idx,]))

exprData = exprData[!idx,]
dim(exprData)
save(exprData,meta,file = 'Data/raw_full.rdata')

###################### calculate cell type marker ##############################
load('Data/raw_full.rdata')
### load raw counts data
counts = read.csv('Data/GEP_Counts_hUSI.csv') ### uploded to figueshare
gene_id = counts$gene_id
gene_id = unlist(lapply(strsplit(gene_id,'\\|'),'[',2))
gene_idx = match(rownames(exprData),gene_id)
counts = counts[,-1]
counts = counts[gene_idx,colnames(exprData)]
rownames(counts) = gene_id[gene_idx]
dim(counts)
counts[1:5,1:5]

### DEGs for non-senescent cells in each cell type
DEGs_celltypes = list()

m = filter(meta,label=='non-senescent')
Expmat=counts[,rownames(m)]
for(c in unique(meta$Cell_type)){
  print(c)
  m$condition = ifelse(m$Cell_type==c,c,'other')
  print(table(m$condition))
  dds <- DESeqDataSetFromMatrix(countData = Expmat, colData = m, design= ~condition)
  dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 20, parallel = FALSE)
  res <- results(dds1, contrast = c('condition', c,'other'))
  res = data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  res = filter(res,padj<0.05&log2FoldChange>0)
  print(nrow(res))
  DEGs_celltypes[[c]] = res
}

############### filter makers calculated in non-senescent samples ##############
load('Data/raw_full.rdata')

idx_senescent = colnames(exprData) %in% filter(meta, label == 'senescent')$Data_id
table(idx_senescent)

X_tr0 = exprData[,idx_senescent]
sen_meta=meta[idx_senescent,]

genenums=unique(c(seq(1,10,1),seq(10,100,5),120,150,200,300,0))
ari_res=NULL
for (gene_num in genenums) {
  MarkerSets <- lapply(DEGs_celltypes, function(x) {
    x<-filter(x,log2FoldChange>1,padj<0.05);
    x<-x[order(x$log2FoldChange,decreasing = T),]; 
    return(na.omit(rownames(x)[0:gene_num]))})
  
  gene=unique(unlist(MarkerSets))
  newexp_h=X_tr0
  newexp_h=newexp_h[which(!(rownames(newexp_h) %in% gene)),]
  
  sce=CreateSeuratObject(log1p(newexp_h),meta.data = sen_meta[,1:13],min.cells = 3,min.features = 200)#log1pTPM
  sce <- FindVariableFeatures(sce, selection.method ="vst", nfeatures =2000)
  sce <- ScaleData(sce,features = rownames(sce))
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- RunUMAP(sce,dims = 1:20)
  sce <- FindNeighbors(sce,dims = 1:20)
  for (res in c(0.6,0.8,1)) {
    sce=FindClusters(sce,resolution = res,algorithm = 2);table(Idents(sce))
  }
  
  for (res in c(0.6,0.8,1)) {
    resd=paste0('RNA_snn_res.',res)
    ari_new1 <- adjustedRandIndex(sce@meta.data$Cell_type, sce@meta.data[,resd])
    dat1=data.frame(gene_num=gene_num,ARI=ari_new1,resolution=res)
    ari_res=rbind(ari_res,dat1)
  }
  
}

###### Top85
gene_num=85
MarkerSets <- lapply(DEGs_celltypes, function(x) {
  x<-filter(x,log2FoldChange>1,padj<0.05);
  x<-x[order(x$log2FoldChange,decreasing = T),]; 
  return(na.omit(rownames(x)[0:gene_num]))})

gene=unique(unlist(MarkerSets))

expadj=X_tr0
expadj=expadj[which(!(rownames(expadj) %in% gene)),]

###### all samples remove cell type marker genes
exprData = exprData[rownames(expadj),]
exprData[,colnames(expadj)] <- expadj
dim(exprData)

save(exprData,meta,file='Data/raw_drop.rdata')

#################################### OOD test ##################################
load('Data/raw_drop.rdata')
##########  generate train and test data
celltypes = unique(meta$Cell_type)
celltypes
### shuffle cell types for 10 times
celltypes_shuffled <- list()
set.seed(2024)

celltypes_shuffled=replicate(3, sample(celltypes,size=length(celltypes)))
chunks <- apply(celltypes_shuffled,2,function(x) split(x, cut(seq_along(x), 10 , labels= FALSE)))

celltype_train = list()
for(t in 1:3){
  for (k in 1:10){
    s = paste0("t_",t,':k_',k)
    print(s)
    celltype_train[[s]] = celltypes[!celltypes %in% chunks[[t]][[k]]]
  }
}

train_test_data = list()
for(i in 1:length(celltype_train)){
  idx = which(meta$Cell_type %in% celltype_train[[i]])
  train_set = exprData[,idx] %>% normalize_exp()
  test_set = exprData[,-idx] %>% normalize_exp()
  train_test_data[[names(celltype_train)[i]]] = list(train_set,test_set)
}

########## training model by different l2 parameter
SenOCLR_L2 = list()
for(l2 in c(0.001,0.01,0.1,1,5,10)){
  for (i in 1:30){
    train_set = train_test_data[[i]][[1]]
    train_set_sen = train_set[,colnames(train_set) %in% rownames(meta[meta$label == 'senescent',])]
    SenOCLR = gelnet( t(train_set_sen), NULL, 0, l2)
    n = paste(paste0('l2_',as.character(l2),":",names(train_test_data)[[i]]))
    print(n)
    SenOCLR_L2[[n]] <- SenOCLR
  }
}

########## score test set by each model
hUSIs_list = list()
for (n in names(SenOCLR_L2)){
  l2 = strsplit(n,':')[[1]][1]
  t = strsplit(n,':')[[1]][2]
  k = strsplit(n,':')[[1]][3]
  print(sprintf("%s; %s; %s",l2,t,k))
  test_set = train_test_data[[paste(t,k,sep = ':')]][[2]]
  hUSIs_list_m <- list()
  for (m in c('spearman','pearson','dot','logistic','gsva','mean')){
    hUSIs = scoreOCLR(profile = test_set,model = SenOCLR_L2[[n]],m = m)
    hUSIs_list_m[[m]] <- hUSIs
  }
  hUSIs_list[[n]] <- hUSIs_list_m
}
###################### define l2 and quantification by AUC #####################
###### calculate AUC by cell type for each model
auc_models = {}
for(q in c('spearman','pearson','logistic','dot')){
  hUSIs = lapply(hUSIs_list, function(x){x[[q]][[1]]})
  auc_l2 = {}
  for(l2 in c(0.001,0.01,0.1,1,5,10)){
    n = names(hUSIs)[grepl(paste0('l2_',as.character(l2),":"),names(hUSIs))]
    print(n)
    auc_celltypes = {}
    for(nn in n){
      pred = hUSIs[[nn]]
      test_label = meta[names(pred),'label']
      samples = names(pred)
      celltypes = meta[samples,'Cell_type']
      for(c in unique(celltypes)){
        idx = which(celltypes == c)
        pred_c = pred[idx]
        test_label_c = test_label[idx]
        auc = pROC::roc(test_label_c, pred_c, levels=c('non-senescent','senescent'),direction = '<') %>% pROC::auc()
        if(!c %in% names(auc_celltypes)){
          auc_celltypes[[c]] = c(auc)}
        else{
          auc_celltypes[[c]] = c(auc_celltypes[[c]],auc)
        }
      }
    }
    auc_l2[[paste0('l2_',as.character(l2))]] <- auc_celltypes
  }
  auc_models[[q]] = auc_l2
}
### comfirm the AUC numbers for each celltype (3 is expected)
lapply(auc_models ,function(x) lapply(x, function(y) lapply(y, length))) %>% unlist() %>% unique()

################################ train OCLR model ##############################
###### l2=1,Q=spearman, final model
load('Data/raw_drop.rdata')
X = normalize_exp(exprData) 
### split data set
idx_senescent = colnames(X) %in% filter(meta, label == 'senescent')$Data_id
table(idx_senescent)

X_tr = X[,idx_senescent]
X_bk = X[,!idx_senescent]

SenOCLR = gelnet(t(X_tr), NULL, 0, 1)

save(SenOCLR,file = 'Data/SenOCLR_l2=1_drop.rdata') ### the final model

############### Compare all quantification methods in final model ##############
hUSIs_list = hUSIs_list[grepl(paste0('l2_',as.character(l2)),names(hUSIs_list))]
hUSIs_list = lapply(hUSIs_list, function(x){y=x$mean;z=x$gsva;
                               x$mean<-NULL;x$gsva<-NULL;
                               x<-lapply(x, function(h) h<-h[[1]])
                               n=names(x);
                               x<-c(x,y,z);
                               names(x)=c(n,paste('mean:',names(y),sep=''),paste('gsva:',names(z),sep=''));
                               return(x)})

auc_models = {}
for(q in names(hUSIs_list[[1]])){
  hUSIs = lapply(hUSIs_list, function(x){x[[q]]})
  auc_l2 = {}
  for(l2 in c(1)){
    n = names(hUSIs)[grepl(paste0('l2_',as.character(l2),":"),names(hUSIs))]
    print(n)
    auc_celltypes = {}
    for(nn in n){
      pred = hUSIs[[nn]]
      test_label = meta[names(pred),'label']
      samples = names(pred)
      celltypes = meta[samples,'Cell_type']
      for(c in unique(celltypes)){
        idx = which(celltypes == c)
        pred_c = pred[idx]
        test_label_c = test_label[idx]
        auc = pROC::roc(test_label_c, pred_c, levels=c('non-senescent','senescent'),direction = '<') %>% pROC::auc()
        # print(paste(n,c,auc))
        if(!c %in% names(auc_celltypes)){
          auc_celltypes[[c]] = c(auc)}
        else{
          auc_celltypes[[c]] = c(auc_celltypes[[c]],auc)
        }
      }
    }
    auc_l2[[paste0('l2_',as.character(l2))]] <- auc_celltypes
  }
  auc_models[[q]] = auc_l2
}

######## performance of hUSI in batch effect, sparsity, and out liner ###########
###### batch effect 
load('Data/SenOCLR_l2=1_drop.rdata')
load('Data/IMR90_OIS_RNAcountsList.rdata')
load('Data/gfd.rdata')

RNAtpmList <- lapply(RNAcountsList[c("PRJNA395363","PRJNA395378","PRJNA449912","PRJNA293984")], function(x){x<-Counts2TPM(x,gfe_70)})
RNAtpmList[['PRJNA503415']] <- Counts2TPM(RNAcountsList[['PRJNA503415']],gfe_72)

genes = Reduce(intersect,lapply(RNAtpmList,function(x)rownames(x)))
names(RNAtpmList) <- NULL
BEexp = do.call(cbind,lapply(RNAtpmList,function(x)x=x[genes,]))
BEexp =  normalize_exp(BEexp)
dim(BEexp)

BEhUSIs = lapply(c('spearman'),function(m){
  data.frame(hUSI=scoreOCLR(BEexp,SenOCLR,m)[[1]] %>% minmax(),
             label = RNAmeta$condition[match(colnames(BEexp),RNAmeta$sample_title)])})

names(BEhUSIs) <- c('spearman')

BEAUC = lapply(BEhUSIs, function(BEhUSI) {
  pROC::roc(BEhUSI$label, BEhUSI$hUSI, levels=c("senescent", "other"),direction = '>') %>% 
    pROC::auc()}) %>% unlist()

BEAUC
####### sparsity
depths=c(0.8,0.6,0.4,0.2)

RNAnormList_zero <- list()
for(i in 1:length(RNAcountsList)){
  id=names(RNAcountsList)[[i]]
  print(id)
  mat <- RNAcountsList[[i]]
  for (depth in depths) {
    tmat = mat
    print(depth)
    tmat = floor(tmat[,-ncol(tmat)] * depth)
    tmat$Symbol = mat$Symbol
    if (id %in% c("PRJNA395363","PRJNA395378","PRJNA449912","PRJNA293984")){tmat <- Counts2TPM(tmat,gfe_70) %>% normalize_exp()}
    if (id %in% c("PRJNA503415")){tmat <- Counts2TPM(tmat,gfe_72) %>% normalize_exp()}
    RNAnormList_zero[[paste(names(RNAcountsList)[[i]],':','depth_',depth,sep = '')]] <- tmat 
    
  }
}

BEhUSIs_sparse = lapply(c('spearman'),function(m){
  lapply(RNAnormList_zero, function(x){
    data.frame(hUSI=scoreOCLR(x ,SenOCLR,m)[[1]] %>% minmax(),
               label = RNAmeta$condition[match(colnames(x),RNAmeta$sample_title)])})
})

names(BEhUSIs_sparse) = c('spearman')

BEAUC_sparse = lapply(BEhUSIs_sparse, function(BEhUSI) {
  unlist(lapply(BEhUSI, function(x){
    pROC::roc(x$label, x$hUSI, levels=c("senescent", "other"),direction = '>') %>% 
      pROC::auc()})) }) %>% unlist()
BEAUC_sparse

###### out liner
lost_rate = c(0.01,0.05,0.1,0.2)

RNAnormList_outliner <- list()
for(i in 1:length(RNAcountsList)){
  print(names(RNAcountsList)[[i]])
  mat <- RNAcountsList[[i]]
  maxv = max(mat[,-ncol(mat)])
  for (rate in lost_rate) {
    print(rate)
    corrdinates = expand.grid(1:nrow(mat), 1:ncol(mat))
    set.seed(123)
    lost = corrdinates[sample(1:nrow(corrdinates),round(nrow(corrdinates)*rate),replace = F),]
    tmat = mat
    for(j in 1:nrow(lost)){
      tmat[lost[,1][j],lost[,2][j]] <- maxv*100
    }
    if (id %in% c("PRJNA395363","PRJNA395378","PRJNA449912","PRJNA293984")){tmat <- Counts2TPM(tmat,gfe_70) %>% normalize_exp()}
    if (id %in% c("PRJNA503415")){tmat <- Counts2TPM(tmat,gfe_72) %>% normalize_exp()}
    RNAnormList_outliner[[paste(names(RNAcountsList)[[i]],':','rate_',rate,sep = '')]] <- tmat 
  }
}

BEhUSIs_outliner = lapply(c('spearman'),function(m){
  lapply(RNAnormList_outliner, function(x){
    data.frame(hUSI=scoreOCLR(x ,SenOCLR,m)[[1]] %>% minmax(),
               label = RNAmeta$condition[match(colnames(x),RNAmeta$sample_title)])})
})

names(BEhUSIs_outliner) = c('spearman')

BEAUC_outliner = lapply(BEhUSIs_outliner, function(BEhUSI) {
  unlist(lapply(BEhUSI, function(x){
    pROC::roc(x$label, x$hUSI, levels=c("senescent", "other"),direction = '>') %>% 
      pROC::auc()}))}) %>% unlist()
BEAUC_outliner

#################### calculate AUC_th, FN, FP,F1 by SSE #########################
####### classify hUSI by SSE 
hUSIs_list = hUSIs_list[grep('l2_1:',names(hUSIs_list))]
hUSIs_list = lapply(hUSIs_list, function(x){x<-x$spearman[[1]]})

hUSIs_class_list = list()
for(n in names(hUSIs_list)){
  print(n)
  hUSIs = hUSIs_list[[n]]
  hUSIs_data = split(hUSIs,meta[names(hUSIs),'Cell_type'])

  hUSIs_class = lapply(hUSIs_data, function(x){class_hUSI(x)})
  names(hUSIs_class) = NULL
  hUSIs_class = unlist(hUSIs_class)

  hUSIs_class_list[[n]] <- hUSIs_class
}
###### AUC_th
AUC_list = {}
for(n in names(hUSIs_list)){
  print(n)
  hUSIs_f = hUSIs_list[[n]]
  hUSIs_class_f = hUSIs_class_list[[n]]
  hUSI_data = data.frame('hUSI'=hUSIs_f,'hUSI_class'=hUSIs_class_f[names(hUSIs_f)])
  
  hUSI_data_list = split(hUSI_data,meta[names(hUSIs_f),'Cell_type'])
  auc_list = lapply(hUSI_data_list, function(x){apply(x, 2, function(y){
    pROC::roc(meta[rownames(x),]$label, y, levels=c("senescent", "non-senescent"),direction = '>') %>% pROC::auc()
  })})
  for (c in names(auc_list)){
    if(!(c %in% names(AUC_list))){
      AUC_list[[c]] <- auc_list[[c]]
    }
    else{
      AUC_list[[c]] <- rbind(AUC_list[[c]],auc_list[[c]])
    }
  }
}

AUC = lapply(AUC_list, function(x) x <- x[,'hUSI']) %>% unlist %>% mean;AUC
AUC_th = lapply(AUC_list, function(x) x <- x[,'hUSI_class']) %>% unlist %>% mean;AUC_th

###### FN, FP for binary results
metrics_list = {}
for(n in names(hUSIs_class_list)){
  print(n)
  hUSIs_class = hUSIs_class_list[[n]]
  hUSI_data_list = split(hUSIs_class,meta[names(hUSIs_class),'Cell_type'])
  m_list <- 
    lapply(hUSI_data_list, function(x){
    pre_label = x %>% factor(levels = c(0,1));
    true_label = ifelse(meta[names(pre_label),]$label == 'senescent', 1, 0) %>% factor(levels = c(0,1));
    conf_matrix <- table(Actual = true_label, Predicted = pre_label);
    FN <- conf_matrix[2, 1];
    FP <- conf_matrix[1, 2];
    TP <- conf_matrix[2, 2];
    
    precision <- TP / (TP + FP);
    recall <- TP / (TP + FN);
    F1_score <- 2 * (precision * recall) / (precision + recall);
    if(is.na(F1_score)){F1_score = 0};
    
    return(c(FN/sum(true_label==1),FP/sum(true_label==0),F1_score))
  })
  for (c in names(m_list)){
    if(!(c %in% names(metrics_list))){
      metrics_list[[c]] <-  m_list[[c]]
    }
    else{
      metrics_list[[c]] <- rbind(metrics_list[[c]], m_list[[c]])
    }
  }
}

metrics <- do.call(rbind,metrics_list) %>% data.frame()
colnames(metrics) <- c('false_negative','false_positive','F1_score')
head(metrics)
mean(metrics$false_negative);mean(metrics$false_positive);mean(metrics$F1_score)

######################## compare with two-class model  ########################## 
load('Data/raw_drop.rdata')
meta$label_b = as.factor(ifelse(meta$label == 'non-senescent',0,1))
table(meta$label_b)
### test 1 OOD 
train_test_data = train_test_data
### test 2 different cell type numbers (k=5,10,20,30)
train_test_data = list()
celltypes = unique(meta$Cell_type)
train_test_data = list()
for(k in c(5,10,20,30)){
  set.seed(233)
  c=replicate(30, sample(celltypes,size=k))
  for(i in 1:ncol(c)){
    idx = which(meta$Cell_type %in% c[,i])
    train_set = exprData[,idx] %>% normalize_exp()
    test_set = exprData[,-idx] %>% normalize_exp()
    train_test_data[[paste("k",as.character(k),"_t",as.character(i),sep = '')]] = list(train_set,test_set)
  }
}

### supported vector machine
library(e1071)
library(caret)
set.seed(233)
SVM <- lapply(train_test_data, function(x){
  train_set = x[[1]]
  test_set = x[[2]]
  rownames(train_set) = gsub('-','_',rownames(train_set))
  rownames(test_set) = gsub('-','_',rownames(test_set))
  train_label = meta[colnames(train_set),'label_b']
  test_label = meta[colnames(test_set),'label_b']
  train_set = t(train_set) %>% as.data.frame()
  test_set = t(test_set) %>% as.data.frame()
  train_set$label = train_label
  ### fit model
  svm_model = svm(label~.,data=train_set,probability = TRUE)
  pred = predict(svm_model,test_set,probability = TRUE)
  pred = attr(x = pred, which = "probabilities")[,1]
  res = list(pred, test_label)
  return(res)})

### random forest
library(randomForest)
set.seed(233)
RF <- lapply(train_test_data, function(x){
  train_set = x[[1]]
  test_set = x[[2]]
  rownames(train_set) = gsub('-','_',rownames(train_set))
  rownames(test_set) = gsub('-','_',rownames(test_set))
  train_label = meta[colnames(train_set),'label_b']
  test_label = meta[colnames(test_set),'label_b']
  train_set = t(train_set) %>% as.data.frame()
  test_set = t(test_set) %>% as.data.frame()
  train_set$label = train_label
  ### fit model
  rf_model = randomForest(label~.,data=train_set)
  pred = predict(rf_model,test_set,type='prob')[,2]
  res = list(pred, test_label)
  return(res)})

### elastic-net
library(glmnet)
set.seed(233)
EN <- lapply(train_test_data, function(x){
  train_set = x[[1]]
  test_set = x[[2]]
  rownames(train_set) = gsub('-','_',rownames(train_set))
  rownames(test_set) = gsub('-','_',rownames(test_set))
  train_label = meta[colnames(train_set),'label_b'] %>% as.vector() %>% as.numeric() %>% as.matrix()
  test_label = meta[colnames(test_set),'label_b']
  train_set = t(train_set) %>% as.matrix()
  test_set = t(test_set) %>% as.matrix()
  ### cross-validation
  cv_model = cv.glmnet(train_set, train_label)
  best_lambda <- cv_model$lambda.min
  print(paste("Best Lambda:", best_lambda))
  ### fit model
  en_model <- glmnet(train_set, train_label, lambda = best_lambda)
  pred =  predict(en_model, s = best_lambda, newx = as.matrix(test_set))
  res = list(pred, test_label)
  return(res)})

### OCLR
library(gelnet)
set.seed(233)
OCLR <- lapply(train_test_data[c(1:3)], function(x){
  train_set = x[[1]]
  test_set = x[[2]]
  rownames(train_set) = gsub('-','_',rownames(train_set))
  rownames(test_set) = gsub('-','_',rownames(test_set))
  train_label = meta[colnames(train_set),'label_b']
  test_label = meta[colnames(test_set),'label_b']
  train_set = as.matrix(train_set[,train_label=="1"]) %>% t
  test_set = as.matrix(test_set) 
  ### fit model
  oclr = gelnet(train_set, NULL, 0, 1)
  pred = scoreOCLR(profile = test_set,model = oclr,m = 'spearman')[[1]]
  res = list(pred, test_label)
  return(res)})

### calculate AUC in OOD test (Test1)
models = list(SVM,RF,EN,OCLR)
names(models) = c('SVM','RF','EN','OCLR')

auc_models = {}
for(m in c('SVM','RF','EN','OCLR')){
  model = models[[m]]
  auc_celltypes = {}
  for(n in names(model)){
    x=model[[n]]
    pred = x[[1]]
    test_label = x[[2]]
    samples = train_test_data[[n]][[2]] %>% colnames()
    celltypes = meta[samples,'Cell_type']
    for(c in unique(celltypes)){
      idx = which(celltypes == c)
      pred_c = pred[idx]
      test_label_c = test_label[idx]
      auc = pROC::roc(test_label_c, pred_c, levels=c('0','1'),direction = '<') %>% pROC::auc()
      if(!c %in% names(auc_celltypes)){
        auc_celltypes[[c]] = c(auc)}
      else{
        auc_celltypes[[c]] = c(auc_celltypes[[c]],auc)
      }
    }
  }
  auc_models[[m]] = auc_celltypes
}

lapply(auc_models,function(x) {unlist(lapply(x,function(y)mean(y))) %>% mean})

### calculate AUC in different cell type numbers (Test2)
models = list(SVM,RF,EN,OCLR)
names(models) = c('SVM','RF','EN','OCLR')

auc_models = {}
for(m in c('SVM','RF','EN','OCLR')){
  model = models[[m]]
  auc_celltypes = {}
  for(n in names(model)){
    x=model[[n]]
    pred = x[[1]]
    test_label = x[[2]]
    samples = train_test_data[[n]][[2]] %>% colnames()
    celltypes = meta[samples,'Cell_type']
    for(c in unique(celltypes)){
      idx = which(celltypes == c)
      pred_c = pred[idx]
      test_label_c = test_label[idx]
      auc = pROC::roc(test_label_c, pred_c, levels=c('0','1'),direction = '<') %>% pROC::auc()
      ids = paste(strsplit(n,'_')[[1]][1],c,sep = '_')
      if(!ids %in% names(auc_celltypes)){
        auc_celltypes[[ids]] = c(auc)}
      else{
        auc_celltypes[[ids]] = c(auc_celltypes[[ids]],auc)
      }
    }
  }
  auc_models[[m]] = auc_celltypes
}

lapply(auc_models,function(x) {unlist(lapply(x[grepl('k5',names(x))],function(y)mean(y))) %>% mean})

################# enrichment model weights on cellmarker database ###############
load('Data/SenOCLR_l2=1_drop.rdata')
### CELLMARKER
cellmarker = read.table("Data/Human_cell_markers.txt",header = T,sep = '\t')
cellmarker = cellmarker[,c('cellName','geneSymbol')]
cellmarker = cellmarker[!is.na(cellmarker$geneSymbol),]
cellmarker = cellmarker[!grepl('^NA$|^NA,',cellmarker$geneSymbol),]

SenSetlist <- loadSetdata()
senemarkers = unique(unlist(SenSetlist))

cellmarker_gmt = data.frame()
for(i in 1:nrow(cellmarker)){
  c = cellmarker$cellName[i]
  print(c)
  genes = strsplit(cellmarker$geneSymbol[i],',')[[1]]
  for(g in genes){
    if(g !='NA'){
        gg=trimws(g, which = c("both"))
        gg=gsub('\\[','',gg)
        gg=gsub('\\]','',gg)
        if(! gg %in% senemarkers){cellmarker_gmt=rbind(cellmarker_gmt,c(c,gg))}
    }
  }
}
colnames(cellmarker_gmt) <- c('term','gene')
length(unique(cellmarker_gmt$term))

SenSetlist_sub <- SenSetlist[c('SenMayo','SenUp')]
senmarker_gmt = do.call(rbind,lapply(names(SenSetlist_sub), function(x){y=data.frame(term=x,gene = SenSetlist_sub[[x]])}))

all_gmt=rbind(cellmarker_gmt,senmarker_gmt)

set.seed(223)
fgsea <- GSEA(sort(SenOCLR$w,decreasing=T), 
              exponent = 1, 
              minGSSize = 0,
              maxGSSize = 20000, 
              pvalueCutoff = 1, 
              pAdjustMethod = "BH", 
              TERM2GENE=all_gmt,
              seed = 223,
              by = "fgsea",
              eps=0)

df_plot = fgsea@result
df_plot = df_plot[order(df_plot$NES,decreasing=T),];dim(df_plot)
df_plot[df_plot$p.adjust<0.05,c('NES','pvalue')]

############################ model weight GSEA ##################################
load('Data/SenOCLR_l2=1_drop.rdata')

library(fgsea)
library(clusterProfiler)

### databases
Hallmarker = read.gmt("Data/h.all.v2023.1.Hs.symbols.gmt")
KEGG = read.gmt('Data/KEGG_hsa.gmt')
GO=read.gmt('Data/c5.go.v2023.1.Hs.symbols.gmt')
Reactome = read.gmt("Data/c2.cp.v2023.2.Hs.symbols.gmt")
Reactome = Reactome[grep('^REACTOME', Reactome$term),]
databases=list('Hallmarker'=Hallmarker,'KEGG'=KEGG,'GO'=GO,'Reactome'=Reactome)

fgseaList=list()
for(d in names(databases)){
  print(d)
  set.seed(233)
  pathways = databases[[d]]
  pathways$term = sub("^[^_]*_", "", pathways$term)
  fgsea <- GSEA(sort(SenOCLR$w,decreasing=T), 
                exponent = 1, 
                minGSSize = 5,
                maxGSSize = 1000, 
                pvalueCutoff = 1, 
                pAdjustMethod = "BH", 
                TERM2GENE=pathways,
                seed = 233,
                by = "fgsea",
                eps=0)
  fgseaList[[d]] <- fgsea
}

re = data.frame()
for(d in names(databases)){
  print(d)
  df = fgseaList[[d]]@result
  df = filter(df, pvalue < 0.05)
  df = df[order(df$NES,decreasing=T),]
  df$Database = d
  print(nrow(df))
  re = rbind(re,df)
}




