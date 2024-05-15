mm_l2 = readRDS("~/wangj/AgingScore/Data/Bulk_TrainModel/mm_l2.rds")
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
### HALMARKER
hallmarker = read.gmt("/mnt/data1/wangj/GeneSets/h.all.v2023.1.Hs.symbols.gmt")
fgsea <- GSEA(sort(mm_l2$w,decreasing=T), 
              exponent = 1, 
              minGSSize = 10,
              maxGSSize = 500, 
              pvalueCutoff = 1, 
              pAdjustMethod = "BH", 
              TERM2GENE=hallmarker,
              seed = 233,
              by = "fgsea",
              eps=0)

df_hallmarker = fgsea@result
df_hallmarker = filter(df_hallmarker, pvalue < 0.05)
df_hallmarker = df_hallmarker[order(df_hallmarker$NES,decreasing=T),]
df_hallmarker$Database = 'HALLMARKER'
df_hallmarker[,c('ID','NES','pvalue')]

### KEGG
KEGG = read.gmt('/mnt/data1/wangj/GeneSets/KEGG_hsa.gmt')
fgsea <- GSEA(sort(mm_l2$w,decreasing=T), 
              exponent = 1, 
              minGSSize = 10,
              maxGSSize = 1000, 
              pvalueCutoff = 1, 
              pAdjustMethod = "BH", 
              TERM2GENE=KEGG,
              seed = 233,
              by = "fgsea",
              eps=0)

df_kegg = fgsea@result
df_kegg = filter(df_kegg, pvalue < 0.05)
df_kegg = df_kegg[order(df_kegg$NES,decreasing=T),]
df_kegg$Database = 'KEGG'
df_kegg[,c('ID','NES','pvalue')]

### Reactome
Reactome = read.gmt("/mnt/data1/wangj/GeneSets/c2.cp.v2023.2.Hs.symbols.gmt")
Reactome = Reactome[grep('^REACTOME', Reactome$term),]
head(Reactome)

fgsea <- GSEA(sort(mm_l2$w,decreasing=T), 
              exponent = 1, 
              minGSSize = 2,
              maxGSSize = 2000, 
              pvalueCutoff = 1, 
              pAdjustMethod = "BH", 
              TERM2GENE=Reactome,
              seed = 233,
              by = "fgsea",
              eps=0)
              
df_reactome = fgsea@result
df_reactome = filter(df_reactome, pvalue < 0.05)
df_reactome = df_reactome[order(df_reactome$NES,decreasing=T),]
df_reactome$Database = 'Reactome'
df_reactome[,c('ID','NES','pvalue')]

### senescence pathways
sene_pathways = read.gmt('/mnt/data1/wangj/GeneSets/senescence_genesets.v2023.2.Hs.gmt')

fgsea_sene <- GSEA(sort(mm_l2$w,decreasing=T), 
                   exponent = 1, 
                   minGSSize = 2,
                   maxGSSize = 1000, 
                   pvalueCutoff = 1, 
                   pAdjustMethod = "BH", 
                   TERM2GENE=sene_pathways_2,
                   seed = 233,
                   by = "fgsea",
                   eps=0)

df = fgsea_sene@result
df = df[order(df$NES,decreasing=T),]
df[,c('NES','pvalue')]

### permutation
t = 1000
Ratio=c()
set.seed(233)
for (i in 1:t){
  weights = mm_l2$w
  names(weights) = sample(names(weights))
  fgsea <- GSEA(sort(mm_l2$w,decreasing=T), 
                     exponent = 1, 
                     minGSSize = 2,
                     maxGSSize = 1000, 
                     pvalueCutoff = 1, 
                     pAdjustMethod = "BH", 
                     TERM2GENE=sene_pathways,
                     seed = 233,
                     by = "fgsea",
                     eps=0)
  df = fgsea@result
  all = nrow(df)
  print(all)
  df = filter(df,pvalue < 0.05)
  sig = nrow(df)
  print(sig)
  r= sig/all
  print(r)
  Ratio = c(Ratio,r)
}
p_value <- sum(Ratio >= 0.6595745) / length(Ratio)

###################### model weight preference ################################
### calculate cell type markers
load('/mnt/data1/wangj/AgingScore/Data/Bulk_TrainModel/ModelTrainData.RData')
idx_senescent = colnames(X) %in% filter(metadata, condition %in% "senescent")$sample_title
X_centre = X - (apply(X, 1, mean))
X_tr = X_centre[,idx_senescent]
X_bk = X_centre[,!idx_senescent]

metadata = read.csv('/home/wangjing/wangj/codebase/HUSI/Figures/revison/train_set_metadata.csv')

celltype_markers <- list()
for(state in c('non','sene')){
  if(state == 'non'){exp=X_bk} else{exp=X_tr}
  pheno = data.frame(sample=colnames(exp),celltype=metadata$celltype[match(colnames(exp),metadata$sample_title)])
  ##remove low expression
  exp = exp[rowMeans(exp)>quantile(rowMeans(exp),0.25),]
  
  cm = list()
  for(c in unique(pheno$celltype)){
    print(c)
    deg.tab <- NULL
    p <- ifelse(pheno$celltype == c,1,0)
    s <- pheno$sample[pheno$celltype == c] 
    for(i in 1:nrow(exp)){
      lm.model <- lm(exp[i,]~p)
      line <- c(summary(lm.model)$coefficients[2,4],summary(lm.model)$coefficients[2,1])
      deg.tab <- rbind(deg.tab,line)
    }
    rownames(deg.tab) <- rownames(exp)
    colnames(deg.tab) <- c("pval","coefficients")
    deg.tab <- as.data.frame(deg.tab)
    deg.tab$p.adjust <- p.adjust(deg.tab$pval,method="BH")
    deg.tab <- deg.tab[deg.tab[,3]<0.05,]
    cm[[c]] <- deg.tab
  }
  
  celltype_markers[[state]] <- cm
}

celltype_up <- lapply(celltype_markers, function(x){
  lapply(x, function(y) {y <- filter(y,coefficients>0);y<-y[order(y$coefficients,decreasing = T),];y<-rownames(y)[1:200];y})})

MarkerSets <- lapply(unique(pheno$celltype), function(x) intersect(celltype_up$non[[x]],celltype_up$sene[[x]]))
names(MarkerSets) <- unique(pheno$celltype)

EnrichSet<-cogena::gmt2list("/mnt/data1/wangj/AgingScore/Comparison/gene_50signatures_merge.gmt")
SenSet = EnrichSet[c(43,49)]
names(SenSet) <- c("SenMayo","SenUp")

all_sets = c(SenSet,MarkerSets)
df = do.call(cbind, lapply(lapply(all_sets, unlist), `length<-`, max(lengths(all_sets)))) %>% data.frame()
df[is.na(df)] <- ''
write.csv(df,file='/home/wangjing/wangj/codebase/HUSI/Figures/revison/valid_GSEA_celltype.csv')

all_sets_weight = lapply(all_sets, function(x){x<-mm_l2$w[x]})

Marker_set = do.call(rbind,lapply(names(all_sets), function(x){y=data.frame(term=x,gene = all_sets[[x]])}))
fgsea <- GSEA(sort(mm_l2$w,decreasing=T), 
              exponent = 1, 
              minGSSize = 2,
              maxGSSize = 2000, 
              pvalueCutoff = 1, 
              pAdjustMethod = "BH", 
              TERM2GENE=Marker_set,
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
################################# RNA-seq ####################################
### GSE60340
s = fread("/mnt/data1/wangj/AgingScore/Data/Bulk_RNA-seq/GSE60340_induced_sene_TPM.csv") %>% 
  column_to_rownames("Gene Name") %>% 
  data.matrix 

s =  s[,grep("Immortal|Adria|H2O2|5-aza",colnames(s))]
colnames(s)
s_CIS =  s %>% apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} ) %>% minmax

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
  
############################# Microarray ##################################
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
