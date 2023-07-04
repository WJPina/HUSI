mm_l2 = readRDS("/home/wangjing/wangj/AgingScore/BulkData/Bulk_TrainModel/mm_l2.rds")
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