mm_l2 = readRDS("/home/wangjing/wangj/AgingScore/Data/Bulk_TrainModel/mm_l2.rds")
################### Valid model in  RNA-seq and microarray #######################
library(annaffy)
library(magrittr)
library(ggplot2)
library(reshape2)
library(GEOquery)
library(tibble)
library(data.table)
library(dplyr)

################################# Batch Effect #################################
data.matrix <- read.table("/home/wangjing/wangj/AgingScore/Data1/Bulk_BatchEffect/batch_IMR90_4OHT.tsv",sep = "\t",header = T,row.names = 1)
s_batch <- data.matrix %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}
################################# Validation #################################
### RNA-seq
setwd("~/wangj/AgingScore/Data/Bulk_RNA-seq/")
### GSE60340
s_IS = fread("GSE60340_induced_sene_TPM.csv") %>% 
  column_to_rownames("Gene Name") %>% 
  data.matrix %>% 
  {
    apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )
  }
### GSE130306
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
setwd("~/wangj/AgingScore/Data/Bulk_Microarray/")
fs = paste(c('GSE19864','GSE16058','GSE83922','GSE11954','GSE100014','GSE77239'),"eSet.Rdata",sep = "_")
ArrayList <- sapply(fs, function(x) mget(load(x)), simplify = TRUE) 
names(ArrayList) <- unlist(lapply(strsplit(names(ArrayList),"_"),"[",1))
GPL = c("GPL570" = "hgu133plus2.db","GPL3921" = "hthgu133a.db","GPL11532" = "hugene11sttranscriptcluster.db")

### calculate huSI
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
################################# Stability #################################
lost_rate = seq(0,0.9,0.1)
set.seed(233)
Lost_ScoreList = list()

### create lost array
for(i in 1:length(ArrayList)){
  gene_id = aafSymbol( rownames(ArrayList[[i]][[1]]), GPL[ArrayList[[i]][[1]]@annotation]) %>% as.character
  mat <- ArrayList[[i]][[1]] %>% 
          exprs %>% 
          set_rownames( gene_id ) %>% 
          melt(.) %>% 
          filter(!Var1 %in% "character(0)") %>% 
          as.data.table %>% 
          dcast.data.table( Var1~Var2, fun.aggregate = mean ) %>% 
          column_to_rownames("Var1") %>% 
          data.matrix 
  rates = list() 
  for (rate in lost_rate) {
    if(rate == 0){tmat = mat}
    else{
      corrdinates = expand.grid(1:nrow(mat), 1:ncol(mat))
      lost = corrdinates[sample(1:nrow(corrdinates),round(nrow(corrdinates)*rate),replace = F),]
      tmat = mat
      for(j in 1:nrow(lost)){
        tmat[lost[,1][j],lost[,2][j]] <- 0
      }
    }
    rates[[paste('rate_',rate,sep = '')]] = tmat %>% apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )})
  }
  Lost_ScoreList[[names(ArrayList)[i]]] = rates
}

### create lost RNA-seq
setwd("~/wangj/AgingScore/Data/Bulk_RNA-seq/")
RNAList = list()
RNAList[['GSE60340']] <- fread("GSE60340_induced_sene_TPM.csv") %>% column_to_rownames("Gene Name") %>% as.matrix
mat_GSE60340 <- do.call(cbind,list.files("GSE130306/", "GSM", full.names = T) %>% 
                              lapply(function(x) {fread(x) %>% filter( !duplicated(gene_name) ) %>% column_to_rownames("gene_name")})) %>% as.matrix
colnames(mat_GSE60340) <- list.files("GSE130306/", "GSM", full.names = T) %>% gsub(".*RNAseq_", "", .) %>% gsub(".txt.gz", "", .)
RNAList[['GSE60340_OIS']] <- mat_GSE60340[,grep('OIS',colnames(mat_GSE60340))]
RNAList[['GSE60340_RS']] <- mat_GSE60340[,grep('RS',colnames(mat_GSE60340))]

for(i in 1:length(RNAList)){
  mat = RNAList[[i]]
  rates = list() 
  for (rate in lost_rate) {
    if(rate == 0){tmat = mat}
    else{
      corrdinates = expand.grid(1:nrow(mat), 1:ncol(mat))
      lost = corrdinates[sample(1:nrow(corrdinates),round(nrow(corrdinates)*rate),replace = F),]
      tmat = mat
      for(j in 1:nrow(lost)){
        tmat[lost[,1][j],lost[,2][j]] <- 0
      }
    }
    rates[[paste('rate_',rate,sep = '')]] = tmat %>% apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )})
  }
  Lost_ScoreList[[names(RNAList)[i]]] = rates
}

save(Lost_ScoreList,file= 'Valid_lost_1013.RData')
########################### model weight GSEA ###################################
library(fgsea)
library(clusterProfiler)

hallmarker = read.gmt("/mnt/data3/wangj2/GeneSets/h.all.v2023.1.Hs.symbols.gmt")
fgsea <- GSEA(sort(mm_l2$w,decreasing=T), 
                exponent = 1, 
                minGSSize = 10,
                maxGSSize = 500, 
                pvalueCutoff = 0.1, 
                pAdjustMethod = "BH", 
                TERM2GENE=hallmarker,
                seed = 233,
                by = "fgsea",
                eps=0)
fgsea@result$ID = gsub('HALLMARK_','',fgsea@result$ID)
rownames(fgsea@result) = fgsea@result$ID
fgsea@result$Description = gsub('HALLMARK_','',fgsea@result$Description)
names(fgsea@geneSets) = gsub('HALLMARK_','',names(fgsea@geneSets))

df = fgsea@result
df = filter(df,pvalue < 0.05)
df = df[order(df$NES,decreasing=T),]
df[,c('ID','NES')]

#################################### valid in CS #############################
setwd('/home/wangjing/wangj/AgingScore/Data/Bulk_TrainModel/CS_score_of_single_cell_datasets/')
source("/mnt/data1/wangj/codebase/HUSI/transID.R")
##### GTEx bulk data
### load cs score
TCSER_gtex = read.csv('Senescence_quantification_GTEX.csv')
dim(TCSER_gtex)
TCSER_gtex = TCSER_gtex[!duplicated(TCSER_gtex$Sample),]
rownames(TCSER_gtex) = TCSER_gtex$Sample
library(CePa)
GTEx_exp <- read.gct("/home/wangjing/wangj/MyProject/Data/GTEx/Bulk_SnRNA/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
colnames(GTEx_exp) <- gsub('\\.','-',colnames(GTEx_exp))

table(rownames(TCSER_gtex) %in% gsub('\\.','-',colnames(GTEx_exp)))

Counts = GTEx_exp[,intersect(rownames(TCSER_gtex),colnames(GTEx_exp))]
dim(Counts)
Counts[1:5,1:5]
Counts = esembl2symbol(version = '26',counts = Counts)

hUSI_gtex <- Counts %>% apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )


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
Counts = TCGA_exp[,rownames(TCSER_tcga)]
dim(Counts)
length(samples)

Counts = esembl2symbol(version = '23',counts = Counts)
hUSI_tcga <- Counts%>% apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )




