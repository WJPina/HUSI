### compare to SENCAN
setwd('/mnt/data1/wangj/codebase/HUSI/Figures/revison/')

library(glmnet)
library(edgeR)
library(dplyr)
library(tibble)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(data.table)

load("/mnt/data1/wangj/AgingScore/Comparison/Methods/SENCAN/SENCAN_classifier.rda")
SENCAN <- function(df) {
  ids=bitr(rownames(df),'SYMBOL','ENSEMBL','org.Hs.eg.db')
  ids=ids[!duplicated(ids$SYMBOL),]
  ids=ids[!duplicated(ids$ENSEMBL),]
  df =df[ids$SYMBOL,]
  rownames(df)=ids$ENSEMBL
  dm <- data.matrix(df)
  # Map gene symbols of the input file to the classifier structure
  map <- match(dge$genes$ensembl_id, rownames(dm))
  
  # Impute missing values if necessary
  dmMapped <- dm[map,]
  dmImputed <- dmMapped
  for (i in 1:ncol(dmImputed)) {
    sample <- dmImputed[,i]
    if (any(is.na(sample))) {
      cors <- rep(NA, ncol(dge$counts))
      names(cors) <- rownames(dge$samples)
      for (j in 1:ncol(dge$counts)) {
        have_val <- which(!is.na(sample))
        cors[j] <- pcaPP::cor.fk(sample[have_val], dge$counts[have_val,j])
      }
      nearest <- tail(order(cors), n=2)
      need_val <- which(is.na(sample))
      dmImputed[need_val,i] <- apply(dge$counts[need_val,nearest], 1, mean)
    }
  }
  
  # Normalize the input sample towards the classifier reference
  dge_pred <- DGEList(counts = dmImputed, genes=data.frame(ensembl_id=dge$genes$ensembl_id))
  dge_pred <- calcNormFactors_WSPref(dge_pred, edger_ref_column)
  dmPredCPM <- cpm(dge_pred, log=T)
  dmPredCPM[dmPredCPM < min(dmCPM)] <- min(dmCPM)
  pred_standardized_expression <- (dmPredCPM - gene_mean) / gene_sd
  
  # Calculate classifier scores
  res <- predict(full_classifier, t(pred_standardized_expression[use_gene_ix,]), s="lambda.1se")[,1]
  scores <- 1/(1+exp(-res))
  
  return(scores)
}

### RNA-seq profiles IMR90
IMR90_OIS <- read.table("/mnt/data1/wangj/AgingScore/Data/Bulk_BatchEffect/batch_IMR90_4OHT.tsv",sep = "\t",header = T,row.names = 1)
SENCAN_scores <- SENCAN(IMR90_OIS)

### plot
library(ggsci)
png('/home/wangjing/wangj/codebase/HUSI/Figures/revison/SENCAN_BE.png',width = 1500,height = 1500,res = 300)
read.table("/mnt/data1/wangj/AgingScore/Data/Bulk_BatchEffect/batch_IMR90_4OHT_add_condition_tsv.txt",sep = "\t",header = T) %>%
  dplyr::select(c("title","study_accession","condition")) %>%
  mutate(hUSI = SENCAN_scores[title]) %>%
  column_to_rownames("title") %>%
  mutate(Condition = factor(ifelse(condition == 'other','Other','Senescent'),levels = c("Other","Senescent"),ordered = T)) %>%
  ggplot(aes(Condition, hUSI)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = study_accession),width = 0.1) +
  scale_colour_jama()+
  geom_line(aes(group = study_accession), linetype="dashed", col="skyblue") + 
  geom_signif(comparisons = list(c("Other","Senescent")),test = "t.test",test.args = c("less"),map_signif_level = T) + 
  ylab('hUSI')+
  ggtitle("OIS in IMR90")+
  theme_classic() %+replace% theme(text = element_text(size = 16),axis.text.x = element_text(size=16))
dev.off()


### RNA-seq profiles
MDAH041_CIS <- fread("/mnt/data1/wangj/AgingScore/Data/Bulk_RNA-seq/GSE60340_induced_sene_TPM.csv") %>% 
                column_to_rownames("Gene Name") %>% 
                data.matrix

WI38 = list.files("/mnt/data1/wangj/AgingScore/Data/Bulk_RNA-seq/GSE130306/", "GSM", full.names = T) %>% 
          lapply(function(x){fread(x) %>% 
          filter( !duplicated(gene_name) ) %>% 
          column_to_rownames("gene_name") %>% 
          data.matrix})
WI38_OIS_RS = sapply(WI38, function(x){x[,1]})
colnames(WI38_OIS_RS) = list.files("/mnt/data1/wangj/AgingScore/Data/Bulk_RNA-seq/GSE130306/", "GSM", full.names = T) %>% gsub(".*RNAseq_", "", .) %>% gsub(".txt.gz", "", .)

### preprocess the data
SENCAN_scores_1 = SENCAN(MDAH041_CIS)
SENCAN_scores_2 = SENCAN(WI38_OIS_RS)

### preprocess the data
MDAH041_CIS_norm <- SENCAN_norm(MDAH041_CIS)
WI38_OIS_RS_norm <- lapply(WI38, function(x){SENCAN_norm(x)})
WI38_OIS_RS_norm = sapply(WI38_OIS_RS_norm, function(x){x[,1]})
colnames(WI38_OIS_RS_norm) = colnames(WI38_OIS_RS)


### plot
mytheme <- function () { 
  theme_classic() %+replace% 
    theme(text = element_text(size = 16),legend.position = "none",
          axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))}

p1.1 <- data.frame(SENCAN_score  = SENCAN_scores_1) %>% 
  mutate(condition= gsub("_.*", "", names(SENCAN_scores_1))) %>%
  .[which(.$condition %in% c("Immortal","Adria","H2O2","5-aza")),] %>%
  mutate(Treatment = factor(condition,levels = c("Immortal","Adria","H2O2","5-aza"),ordered = T)) %>% 
  ggplot(aes(Treatment, SENCAN_score )) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Treatment),width = 0.4,size=2)+
  geom_signif(comparisons = list(c("Immortal","H2O2"),c("Immortal","5-aza")),
              test = "t.test",
              step_increase=0.1,
              map_signif_level = T,
              test.args = c("less")) + 
  mytheme()+
  ylab('SENCAN')+
  ggtitle("CIS in MDAH041")

p1.2<- data.frame(SENCAN_score = SENCAN_scores_2[grepl('OISD[0|2|4|6|10]',names(SENCAN_scores_2))]) %>% 
  mutate(condition= gsub("_.*", "", rownames(.))) %>%
  mutate(Treatment = factor(condition,levels = c("OISD0","OISD2","OISD4","OISD6",'OISD10'),ordered = T)) %>% 
  ggplot(aes(Treatment, SENCAN_score )) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Treatment),width = 0.4,size=2)+
  geom_signif(comparisons = list(c("OISD0","OISD2"),c("OISD2","OISD4"),c("OISD4","OISD6")),
              test = "t.test",
              step_increase=0.1,
              map_signif_level = T,
              test.args = c("less")) + 
  mytheme()+
  ylab('SENCAN')+
  ggtitle("OIS in WI-38")+
  ylim(0,1)

p1.3<- data.frame(SENCAN_score = SENCAN_scores_2[grepl('RS',names(SENCAN_scores_2))]) %>% 
  mutate(condition= gsub('_.*','',gsub("RS_", "", rownames(.)))) %>%
  mutate(Condition = factor(condition,levels = c("Proliferative","Senescent"),ordered = T)) %>% 
  ggplot(aes(Condition, SENCAN_score )) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Condition),width = 0.4,size=2)+
  geom_signif(comparisons = list(c("Proliferative","Senescent")),
              test = "t.test",
              step_increase=0.1,
              map_signif_level = T,
              test.args = c("less")) + 
  mytheme()+
  ylab('SENCAN')+
  ggtitle("RS in WI-38")+
  ylim(0,1)

fig <- ggarrange(p1.1,p1.2,p1.3,ncol = 3,nrow = 1,common.legend = F,widths = c(4,4,3))
png('/home/wangjing/wangj/codebase/HUSI/Figures/revison/SENCAN_RNA-seq.png',width = 3000,height = 1500,res = 300)
fig
dev.off()



