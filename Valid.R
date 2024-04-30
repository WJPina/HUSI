mm_l2 = readRDS("mm_l2.rds")
################## Valid model in  RNA-seq and microarray #######################
library(annaffy)
library(magrittr)
library(ggplot2)
library(reshape2)
library(GEOquery)
library(tibble)
library(data.table)
library(dplyr)
source('~/wangj/codebase/HUSI/functions.R')
########################### model weight GSEA ###################################
library(fgsea)
library(clusterProfiler)

hallmarker = read.gmt("/mnt/data1/wangj/GeneSets/h.all.v2023.1.Hs.symbols.gmt")
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
fgsea@result$ID = gsub('HALLMARK_','',fgsea@result$ID)
rownames(fgsea@result) = fgsea@result$ID
fgsea@result$Description = gsub('HALLMARK_','',fgsea@result$Description)
names(fgsea@geneSets) = gsub('HALLMARK_','',names(fgsea@geneSets))


KEGG = read.gmt('/mnt/data1/wangj/GeneSets/KEGG_hsa.gmt')
fgsea <- GSEA(sort(mm_l2$w,decreasing=T), 
              exponent = 1, 
              minGSSize = 10,
              maxGSSize = 500, 
              pvalueCutoff = 0.05, 
              pAdjustMethod = "BH", 
              TERM2GENE=KEGG,
              seed = 233,
              by = "fgsea",
              eps=0)
fgsea@result$ID = gsub('KEGG_','',fgsea@result$ID)
rownames(fgsea@result) = fgsea@result$ID
fgsea@result$Description = gsub('KEGG_','',fgsea@result$Description)
names(fgsea@geneSets) = gsub('KEGG_','',names(fgsea@geneSets))


GO = read.gmt('/mnt/data1/wangj/GeneSets/c5.go.v2023.1.Hs.symbols.gmt')
fgsea <- GSEA(sort(mm_l2$w,decreasing=T), 
              exponent = 1, 
              minGSSize = 10,
              maxGSSize = 500, 
              pvalueCutoff = 0.05, 
              pAdjustMethod = "BH", 
              TERM2GENE=GO,
              seed = 233,
              by = "fgsea",
              eps=0)
fgsea@result$ID = gsub('GO_','',fgsea@result$ID)
rownames(fgsea@result) = fgsea@result$ID
fgsea@result$Description = gsub('GO_','',fgsea@result$Description)
names(fgsea@geneSets) = gsub('GO_','',names(fgsea@geneSets))

library(reactome.db)
# gtf = rtracklayer::import("/mnt/data1/wangj/AgingScore/Data/Bulk_TrainModel/gencode.v31.annotation.gtf")
# ENST2ENTR = fread("/mnt/data1/wangj/GeneSets/gencode.v31lift37.metadata.EntrezGene.gz") %>% data.table()
# colnames(ENST2ENTR) <- c("transcript_id",'ENTREZID')
# ENST2SYM = as.data.table(gtf)
# ENST2SYM = ENST2SYM[!(ENST2SYM$transcript_id == 'NA')]
# ENST2SYM = ENST2SYM[,c('transcript_id','gene_name')]
# ids = left_join(ENST2SYM,ENST2ENTR,by='transcript_id')
# ids = ids[,c('ENTREZID','gene_name')]
# ids = ids[!duplicated(ids$gene_name)]
# ids = data.frame(ids)
# rownames(ids) = ids$gene_name 
# 
# symbols = names(mm_l2$w[!is.na(ids[names(mm_l2$w),'ENTREZID'])])
# entrezids = ids[symbols,'ENTREZID']

ids = bitr(names(mm_l2$w),'SYMBOL','ENTREZID','org.Hs.eg.db')
ids = ids[!duplicated(ids$SYMBOL),]
rownames(ids) = ids$SYMBOL

entrezids = ids[names(mm_l2$w),'ENTREZID']

Reactome = reactomePathways(entrezids %>% as.character())
Reactome = do.call(rbind,Reactome) %>% melt
Reactome = Reactome[c(1,3)]
colnames(Reactome) = c('term','gene')

ranks = mm_l2$w[symbols]
names(ranks) = entrezids

fgsea <- GSEA(sort(ranks,decreasing=T), 
              exponent = 1, 
              minGSSize = 2,
              maxGSSize = 500, 
              pvalueCutoff = 0.05, 
              pAdjustMethod = "BH", 
              TERM2GENE=Reactome,
              seed = 233,
              by = "fgsea",
              eps=0)

df = fgsea@result
df = df[order(df$NES,decreasing=T),]
df[,c('ID','NES','pvalue')]


################################# Batch Effect #################################
data.matrix <- read.table("/mnt/data1/wangj/AgingScore/Data/Bulk_BatchEffect/batch_IMR90_4OHT.tsv",sep = "\t",header = T,row.names = 1)
s_batch <- data.matrix %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )} %>% minmax
################################# Validation #################################
###### Microarray
fs = paste(c('GSE19864','GSE16058','GSE83922','GSE11954','GSE100014','GSE77239'),"eSet.Rdata",sep = "_")
ArrayList <- sapply(fs, function(x) mget(load(x)), simplify = TRUE) 
names(ArrayList) <- unlist(lapply(strsplit(names(ArrayList),"_"),"[",1))
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

###### RNA-seq
### GSE60340
s_IS = fread("/mnt/data1/wangj/AgingScore/Data/Bulk_RNA-seq/GSE60340_induced_sene_TPM.csv") %>% 
  column_to_rownames("Gene Name") %>% 
  data.matrix %>% 
  {
    apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} ) %>% minmax
  }
### GSE130306
library(magrittr)
s = list.files("/mnt/data1/wangj/AgingScore/Data/Bulk_RNA-seq/GSE130306/", "GSM", full.names = T) %>% 
  sapply(function(x) {
    fread(x) %>% 
      filter( !duplicated(gene_name) ) %>% 
      column_to_rownames("gene_name") %>% 
      data.matrix %>% 
      {
        apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )
      }
  }) %>% 
  set_names( list.files("/mnt/data1/wangj/AgingScore/Data/Bulk_RNA-seq/GSE130306/", "GSM", full.names = T) %>% 
               gsub(".*RNAseq_", "", .) %>% gsub(".txt.gz", "", .) )

s_OIS = s[grepl('OISD[0|2|4|6|10]',names(s))] %>% minmax
s_RS = s[grepl('RS',names(s))] %>% minmax
  
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
