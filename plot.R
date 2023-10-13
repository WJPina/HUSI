library(ggplot2)
library(dplyr)
library(tibble)
library(data.table)
library(annaffy)
library(ggsignif)
library(stringi)
library(ggpubr)
library(ggsci)
library(ggrepel)

mytheme <- function () { 
    theme_classic() %+replace% 
    theme(text = element_text(size = 16),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "right",
            axis.text.x = element_blank())
    }

### plot train Rep.CS heatmap
load('ModelTrainData.RData')
meta = read.csv('metadata_train.csv') 
meta = tibble::column_to_rownames(meta,'sample_title')
table(colnames(X) %in% rownames(meta))

genes = read.csv('Kasit_2019_31560156_AgingGenes.csv')
genes = genes[genes$Gene.symbol %in% rownames(X),]
table(genes$signature)

groups = meta$treatment
names(groups) = rownames(meta)
my_order <- sort(table(groups),decreasing = F)
groups <- groups[order(factor(groups,levels = names(my_order)))]
groups <- factor(groups,levels = names(my_order))

library(ComplexHeatmap)
library(circlize)
mypalette = read.csv('~/scripts/colors.csv')
bar = c('#ffd7ba',"#ffcdb2", "#ffb4a2", "#e5989b", "#b5838d", "#6d6875")
names(bar) = levels(groups)
top_anno <- HeatmapAnnotation(df = data.frame(Condition = groups),show_legend = T,col = list(Condition = bar))
left_anno = rowAnnotation(df = data.frame(Signature= rep(c("up","down"),times=c(519,729))),show_legend = T,
                          col = list(Signature = c('up' = '#e63946','down' = '#023e8a')))
col <- colorRamp2(c(-1.5,0,1.5), c("#023e8a","white", "#e63946"), space = "LAB")

png('/home/wangjing/wangj/codebase/HUSI/Figures/model/CS_sinature.png',width = 2000,height = 1500,res= 300)
mat = log2(exprData[genes$Gene.symbol,names(groups)] + .00001) %>% scale()
Heatmap(mat,
        show_column_names = F,
        show_row_names = F,
        row_title = NULL,
        col = col,
        cluster_rows = T,
        cluster_row_slices = FALSE,
        cluster_columns = F,
        top_annotation = top_anno,
        left_annotation = left_anno,
        row_names_gp = gpar(fontsize = 12),
        row_split = factor(rep(c("up","down"),times=c(519,729)),levels = c('up','down'),ordered = T),
        column_split = factor(rep(c('Senescent','other'),times = c(126,122)),levels=c('Senescent','other'),ordered = T))
dev.off()

### model leave-one-out auc
load("/home/wangjing/wangj/AgingScore/Data/Bulk_TrainModel/model_auc.RData")
auc = sort(auc)

png('/home/wangjing/wangj/codebase/HUSI/Figures/model/cross-validation.png',width = 1000,height = 800,res = 300)
ggplot(aes(x = 1:length(auc),y = auc),data=data.frame(auc=auc))+
    geom_line()+
    theme_classic()+
    theme(text=element_text(size=12),axis.title.x = element_blank())+
    ggtitle('Leave-one-out cross-validation')+
    ylab("Correctly Ranks Rate")
dev.off()

### bacth effect
load('/home/wangjing/wangj/AgingScore/AgingScorePro/Data1_Scripts/ModelValidData_Batch.RData')

png('/home/wangjing/wangj/codebase/HUSI/Figures/model/valid_BE.png',width = 1500,height = 1500,res = 300)
read.table("/home/wangjing/wangj/AgingScore/Data/Bulk_BatchEffect/batch_IMR90_4OHT_add_condition_tsv.txt",sep = "\t",header = T) %>%
  dplyr::select(c("title","study_accession","condition")) %>%
  mutate(hUSI = s_batch[title]) %>%
  column_to_rownames("title") %>%
  mutate(condition = factor(condition,levels = c("other","senescent"),ordered = T)) %>%
  ggplot(aes(condition, hUSI)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = study_accession),width = 0.1) +
  scale_colour_jama()+
  geom_line(aes(group = study_accession), linetype="dashed", col="skyblue") + 
  geom_signif(comparisons = list(c("other","senescent")),test = "t.test",test.args = c("less"),map_signif_level = T) + 
  ylab('hUSI')+
  xlab('Condition')+
  ggtitle("OIS of IMR90")+
  theme_classic() %+replace% theme(text = element_text(size = 16),axis.text.x = element_text(size=16))
dev.off()

### valid data in micro-array
names(ArrayList)
p2.1 <- ArrayList [["GSE19864"]][[1]] %>% 
  pData %>% 
  dplyr::select(title, geo_accession) %>% 
  mutate(condition=case_when(
    grepl("Growing", title, ignore.case = T) ~ "Growing", 
    grepl("Senescent", title, ignore.case = T) ~ "Senescent", 
    grepl("Confluent", title, ignore.case = T) ~ "confluent", 
    grepl("0.1% FBS", title, ignore.case = T) ~ "low serum"
  )) %>% 
  .[which(.$condition %in% c("Growing","Senescent")),] %>%
  mutate(condition = factor(condition,levels = c("Growing","Senescent"),ordered = T)) %>% 
  inner_join(
    data.frame(sene_score=ScoreList[["s_19864"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  ggplot(aes(condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = condition),width = 0.1)+
  geom_signif(comparisons = list(c("Growing", "Senescent")),test = "t.test",test.args = c("less"), map_signif_level = T) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle("OIS in IMR90")

p2.2 <- ArrayList [["GSE16058"]][[1]] %>% 
  pData %>% 
  mutate(condition=gsub("growth status: ", "", characteristics_ch1.3) ) %>% 
  mutate(passage = `passage:ch1`) %>%
  inner_join(
    data.frame(sene_score=ScoreList[["s_16058"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("Growing","Senescent")),] %>% 
  mutate(condition = factor(condition,levels = c("Growing","Senescent"),ordered = T)) %>% 
  ggplot(aes(condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = condition,shape =passage),width = 0.1) +
  geom_signif(comparisons = list(c("Growing", "Senescent")),test = "t.test",test.args = c("less"), map_signif_level = T) + 
  mytheme()+
  guides(color = 'none')+
  ylab('hUSI')+
  ggtitle('RS in HMEC')

p2.3 <- ArrayList [["GSE83922"]][[1]] %>% 
  pData %>% 
  mutate(condition= `cell phenotype:ch1`) %>% 
  inner_join(
    data.frame(sene_score=ScoreList[["s_83922"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("young","senescent")),] %>% 
  mutate(condition = factor(condition,levels = c("young","senescent"),ordered = T)) %>%
  ggplot(aes(condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = condition),width = 0.1) +
  geom_signif(comparisons = list(c("young", "senescent")),test = "t.test",test.args = c("less"), map_signif_level = T) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle('RS in Melanocyte')

p2.4 <- ArrayList [["GSE11954"]][[1]] %>% 
  pData %>% 
  mutate(condition=case_when(
    grepl("growing", description, ignore.case = T) ~ "Growing", 
    grepl("senescent", description, ignore.case = T) ~ "Senescent")) %>% 
  inner_join(
    data.frame(sene_score=ScoreList[["s_11954"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("Growing","Senescent")),] %>% 
  mutate(condition = factor(condition,levels = c("Growing","Senescent"),ordered = T)) %>% 
  ggplot(aes(condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = condition),width = 0.1) +
  geom_signif(comparisons = list(c("Growing","Senescent")),test = "t.test", map_signif_level = T) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle('CIS in HSC')

p2.5 <- ArrayList [["GSE100014"]][[1]] %>% 
  pData %>% 
  mutate(condition=case_when(
    grepl("Proliferating", title, ignore.case = T) ~ "Proliferating", 
    grepl("senescence", title, ignore.case = T) ~ "Senescent")) %>% 
  inner_join(
    data.frame(sene_score=ScoreList[["s_100014"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("Proliferating","Senescent")),] %>% 
  mutate(condition = factor(condition,levels = c("Proliferating","Senescent"),ordered = T)) %>% 
  ggplot(aes(condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = condition),width = 0.1) +
  geom_signif(comparisons = list(c("Proliferating","Senescent")),test = "t.test",test.args = c("less"), map_signif_level = T) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle('CIS in HBEC')

p2.6 <- ArrayList [["GSE77239"]][[1]] %>% 
  pData %>% 
  mutate(condition=case_when(
    grepl("young", `cells:ch1`, ignore.case = T) ~ "young", 
    grepl("old", `cells:ch1`, ignore.case = T) ~ "old")) %>% 
  mutate(condition = factor(condition,levels = c("young","old"),ordered = T)) %>% 
  rename("treatment:ch1" = "treatment") %>% 
  inner_join(
    data.frame(sene_score=ScoreList[["s_77239"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  ggplot(aes(condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = condition,shape = treatment),width = 0.1) +
  geom_line(aes(group = treatment), linetype="dashed", col="skyblue") + 
  geom_signif(comparisons = list(c("young","old")),test = "t.test",test.args = c("less"), map_signif_level = T) + 
  mytheme()+
  guides(color = 'none')+
  ylab('hUSI')+
  ggtitle('RS in HCAEC')

fig <- ggarrange(p2.1,p2.2,p2.3,p2.4,p2.5,p2.6,ncol = 3,nrow = 2,common.legend = F)+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
png('/home/wangjing/wangj/codebase/HUSI/Figures/model/valid_array.png',width = 4500,height = 2400,res = 300)
fig
dev.off()

### valid data in RNA-seq
load("/home/wangjing/wangj/AgingScore/AgingScorePro/Data1_Scripts/ModelValidData.RData")

p1.1 <- data.frame(hUSI = s_IS) %>% 
  mutate(condition= gsub("_.*", "", names(s_IS))) %>%
  .[which(.$condition %in% c("Immortal","Adria","H2O2","5-aza")),] %>%
  mutate(condition = factor(condition,levels = c("Immortal","Adria","H2O2","5-aza"),ordered = T)) %>% 
  ggplot(aes(condition, hUSI)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = condition),width = 0.1)+
  geom_signif(comparisons = list(c("Immortal","Adria"),c("Immortal","H2O2"),c("Immortal","5-aza")),
              test = "t.test",
              step_increase=0.1,
              map_signif_level = T,
              test.args = c("less")) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle("CIS in MDAH041")

p1.2<- data.frame(hUSI = s_RS[grepl('OISD[0|2|4|6|10]',names(s_RS))]) %>% 
  mutate(condition= gsub("_.*", "", rownames(.))) %>%
  mutate(condition = factor(condition,levels = c("OISD0","OISD2","OISD4","OISD6",'OISD10'),ordered = T)) %>% 
  ggplot(aes(condition, hUSI)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = condition),width = 0.1)+
  geom_signif(comparisons = list(c("OISD0","OISD2"),c("OISD2","OISD4"),c("OISD4","OISD6")),
              test = "t.test",
              step_increase=0.1,
              map_signif_level = T,
              test.args = c("less")) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle("OIS in WI-38")

p1.3<- data.frame(hUSI = s_RS[grepl('RS',names(s_RS))]) %>% 
  mutate(condition= gsub('_.*','',gsub("RS_", "", rownames(.)))) %>%
  mutate(condition = factor(condition,levels = c("Proliferative","Senescent"),ordered = T)) %>% 
  ggplot(aes(condition, hUSI)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = condition),width = 0.1)+
  geom_signif(comparisons = list(c("Proliferative","Senescent")),
              test = "t.test",
              step_increase=0.1,
              map_signif_level = T,
              test.args = c("less")) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle("RS in WI-38")

fig <- ggarrange(p1.1,p1.2,p1.3,ncol = 3,nrow = 1,common.legend = F)
png('/home/wangjing/wangj/codebase/HUSI/Figures/model/valid_RNA-seq.png',width = 4500,height = 1200,res = 300)
fig
dev.off()

### GSEA
library(enrichplot)
library(gggsea)

# load('/home/wangjing/wangj/AgingScore/AgingScorePro/Data1_Scripts/Model_GSEA.RData')
cols = c(rev(colorRampPalette(c("transparent", 'red'))(10))[1:7],colorRampPalette(c("transparent", 'blue'))(10)[4:10])
names(cols) = df$ID
png('/home/wangjing/wangj/codebase/HUSI/Figures/model/valid_GSEA_DOWN.png',width = 1500,height = 1500,res = 300)
source('~/wangj/codebase/HUSI/mygseaplot2.R')
mygseaplot2(fgsea,
          geneSetID = df$ID[8:14],
          # geneSetID = df$ID[1:7],
          title = "Negatively enriched hallmarker gene sets",
          # title = "Positively enriched hallmark gene sets",
          color= cols[8:14],
          # color = cols[1:7],
          base_size = 12,
          rel_heights = c(1, 0.2, 0.4),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line"
)
dev.off()

### validate in GTEx
png('/home/wangjing/wangj/codebase/HUSI/Figures/model/valid_GTEx_TCSER.png',width = 1200,height = 900,res = 300)
TCSER_gtex[names(hUSI_gtex),] %>% mutate(hUSI_gtex = hUSI_gtex) %>%
  ggscatter(x = "hUSI_gtex", y = "CS.score",
            color = '#8d99ae',
            add = "reg.line", 
            conf.int = TRUE,
            size = 1,
            add.params = list(color = "#1d3557"),
            ggtheme = theme_classic())+ 
  stat_cor(method = "spearman",label.x = -.08, label.y = 1,color='black')+
  ggtitle('Correlation in GTEx')+
  theme_classic()+
  theme(text = element_text(size = 12))
dev.off()

### validate in TCGA
png('/home/wangjing/wangj/codebase/HUSI/Figures/model/valid_TCGA_TCSER.png',width = 1200,height = 900,res = 300)
TCSER_tcga[names(hUSI_tcga),] %>% mutate(hUSI_tcga = hUSI_tcga) %>%
  ggscatter(x = "hUSI_tcga", y = "score",
            color = '#8d99ae',
            add = "reg.line", 
            conf.int = TRUE,
            size = 1,
            add.params = list(color = "#1d3557"),
            ggtheme = theme_classic())+ 
  stat_cor(method = "spearman",label.x = -.09, label.y = 1.2,color='black')+
  ggtitle('Correlation in TCGA')+
  theme_classic()+
  ylab('CS.score')+
  theme(text = element_text(size = 12))
dev.off()


### plot lost rate stability
mytheme <- function () { 
  theme_classic() %+replace% 
    theme(text = element_text(size = 16),
          legend.position = "right")
}

p2.1 <- ArrayList [["GSE19864"]][[1]] %>% 
  pData %>% 
  dplyr::select(title, geo_accession) %>% 
  mutate(condition=case_when(
    grepl("Growing", title, ignore.case = T) ~ "Growing", 
    grepl("Senescent", title, ignore.case = T) ~ "Senescent", 
    grepl("Confluent", title, ignore.case = T) ~ "confluent", 
    grepl("0.1% FBS", title, ignore.case = T) ~ "low serum"
  )) %>% 
  .[which(.$condition %in% c("Growing","Senescent")),] %>%
  mutate(condition = factor(condition,levels = c("Growing","Senescent"),ordered = T)) %>% 
  inner_join(
    data.frame(hUSI=Lost_ScoreList[["GSE19864"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>%
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  ggplot(aes(x = rate,y = hUSI,color = condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle("GSE19864")
  

p2.2 <- ArrayList [["GSE16058"]][[1]] %>% 
  pData %>% 
  mutate(condition=gsub("growth status: ", "", characteristics_ch1.3) ) %>% 
  mutate(passage = `passage:ch1`) %>%
  inner_join(
    data.frame(hUSI=Lost_ScoreList[["GSE16058"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("Growing","Senescent")),] %>% 
  mutate(condition = factor(condition,levels = c("Growing","Senescent"),ordered = T)) %>%
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  ggplot(aes(x = rate,y = hUSI,color = condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle("GSE16058")

p2.3 <- ArrayList [["GSE83922"]][[1]] %>% 
  pData %>% 
  mutate(condition= `cell phenotype:ch1`) %>% 
  inner_join(
    data.frame(hUSI=Lost_ScoreList[["GSE83922"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("young","senescent")),] %>% 
  mutate(condition = factor(ifelse(condition == "young","Growing","Senescent"),levels = c("Growing","Senescent"),ordered = T)) %>%
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  ggplot(aes(x = rate,y = hUSI,color = condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle("GSE83922")


p2.4 <- ArrayList [["GSE11954"]][[1]] %>% 
  pData %>% 
  mutate(condition=case_when(
    grepl("growing", description, ignore.case = T) ~ "Growing", 
    grepl("senescent", description, ignore.case = T) ~ "Senescent")) %>% 
  inner_join(
    data.frame(hUSI=Lost_ScoreList[["GSE11954"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("Growing","Senescent")),] %>% 
  mutate(condition = factor(condition,levels = c("Growing","Senescent"),ordered = T)) %>% 
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  ggplot(aes(x = rate,y = hUSI,color = condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle("GSE11954")


p2.5 <- ArrayList [["GSE100014"]][[1]] %>% 
  pData %>% 
  mutate(condition=case_when(
    grepl("Proliferating", title, ignore.case = T) ~ "Proliferating", 
    grepl("senescence", title, ignore.case = T) ~ "Senescent")) %>% 
  inner_join(
    data.frame(hUSI=Lost_ScoreList[["GSE100014"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("Proliferating","Senescent")),] %>% 
  mutate(condition = factor(ifelse(condition == 'Proliferating',"Growing","Senescent"),levels = c("Growing","Senescent"),ordered = T)) %>%
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  ggplot(aes(x = rate,y = hUSI,color = condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle("GSE100014")


p2.6 <- ArrayList [["GSE77239"]][[1]] %>% 
  pData %>% 
  mutate(condition=case_when(
    grepl("young", `cells:ch1`, ignore.case = T) ~ "young", 
    grepl("old", `cells:ch1`, ignore.case = T) ~ "old")) %>% 
  mutate(condition = factor(ifelse(condition == 'young',"Growing","Senescent"),levels = c("Growing","Senescent"),ordered = T)) %>% 
  rename("treatment:ch1" = "treatment") %>% 
  inner_join(
    data.frame(hUSI=Lost_ScoreList[["GSE77239"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>%
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  ggplot(aes(x = rate,y = hUSI,color = condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle("GSE77239")


fig2 <- ggarrange(p2.1,p2.2,p2.3,p2.4,p2.5,p2.6,ncol = 2,nrow = 3,common.legend = T)
png('/home/wangjing/wangj/codebase/HUSI/Figures/model/valid_lost.png',width = 3000,height = 3000,res = 300)
fig2
dev.off()


###### comparision auc
### laf
df_plot = do.call(rbind,lapply(AUClist,function(x) {data.frame(t(x[['laf']]))}))

### marker
comm_markers = c("hUSI","GLB1","TP53","CDKN1A","CDKN2A","LMNB1","RB1","CDK1","CDK4","MKI67","CDKN2B")
df_plot = do.call(rbind,lapply(AUClist,function(x) {data.frame(t(x[['marker']][comm_markers,]))}))

### ssGSEA
df_plot = do.call(rbind,lapply(AUClist,function(x) {data.frame(t(x[['ssgsea']]))}))

png('/home/wangjing/wangj/codebase/HUSI/Figures/model/compare_ssgsea.png',width = 1500,height = 1500,res = 300)
df_plot %>%
    mutate(dataset = lapply(strsplit(rownames(df_plot),'\\.'),'[',1) %>% unlist)%>%
    reshape2::melt(value.name = "AUC",variable.name = 'Method',ids = 'dataset') %>%
    group_by(dataset) %>%
    mutate(method=forcats::fct_reorder(Method, AUC, mean, .desc = T) ) %>% 
    ggplot(aes(method, AUC)) + 
    geom_boxplot(aes(color=method), outlier.shape = NA, width=0.5) + 
    geom_jitter(aes(color=method), width = 0.2,size = 0.1) + 
    labs(x=NULL, y="AUC") +
    theme_classic()+
    theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1))+
    facet_wrap(~dataset,ncol = 2,scales = "free")
dev.off()


df_plot_list = lapply(AUClist, function(x) {lapply(x, function(y) {rank(-rowMeans(y))})}) 

df_plot = lapply(df_plot_list, function(x) {x[['marker']]})
df_plot = do.call(cbind, lapply(lapply(df_plot, unlist), `length<-`, max(lengths(df_plot)))) %>% data.frame
df_plot[is.na(df_plot)] = 0
df_plot = reshape2::melt(as.matrix(df_plot),value.name = "AUC_rank")
colnames(df_plot) = c('method','dataset','AUC_rank')
aov.mean<-aggregate(df_plot$AUC_rank,by=list(df_plot$method),FUN=mean)
aov.sd<-aggregate(df_plot$AUC_rank,by=list(df_plot$method),FUN=sd)
aov<-data.frame(aov.mean,sd=aov.sd$x)

png('/home/wangjing/wangj/codebase/HUSI/Figures/model/compare_marker_rank.png',width = 1500,height = 1500,res = 300)
ggplot(data=aov,aes(x=reorder(Group.1,x), y=x,fill = reorder(Group.1,x)))+
    geom_bar(stat="identity",position="dodge")+
    geom_errorbar(aes(ymax=x+sd,ymin=ifelse(x-sd <0,0,x-sd)),position=position_dodge(0.9),width=0.15)+
    theme_classic()+
    xlab('Marker')+
    # xlab('Gene Set')+
    ylab('AUC rank in 4 single-cell datasets')+
    scale_fill_manual(values = colorRampPalette(c("#F44336","#0D47A1"))(nrow(aov)))+
    theme(legend.position = 'none',axis.text.x = element_text(angle = 45,vjust = 1, hjust=1),text=element_text(size=16))
dev.off()

### overlap gene set plot
genes = intersect(names(mm_l2$w),unique(unlist(EnrichSet)))
mat = sapply(EnrichSet, function(set){ifelse(genes %in% set,1,0)})
mat = cbind(mat,hUSI=((mm_l2$w-min(mm_l2$w))/(max(mm_l2$w)-min(mm_l2$w)))[genes])

library(ComplexHeatmap)
library(circlize)
color <- colorRamp2(c(0, 1), c("#1d3557", "#e63946"))
png('/home/wangjing/wangj/codebase/HUSI/Figures/model/compare_geneSet_overlap.png',width = 3000,height = 3000,res = 300)
circos.clear()
circos.par(gap.after = c(40))
circos.heatmap(mat[1:20,], 
               col = color,
               dend.side = "inside",
               rownames.side = "none",
               track.height = 0.6,
               cell.border = 'white',
               cell.lwd = 0.1)
grid.draw(Legend(title = "hUSI_scaled",title_gp = gpar(fontsize = 16, fontface = "bold"),title_position = "topcenter", col_fun = color,size=unit(10,'mm')))
dev.off()
############################### menanome #######################################
load('~/wangj/AgingScore/Data/scRNA_melanome/melanoma.RData')

bar = c("#2a9d8f","#e9c46a","#e76f51")
names(bar) = c("Cycling","Transitional","Senescent")
### melanome trajectory
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/Melanoma_diffusion.png',width = 1200,height = 1000,res = 300)
# DimPlot(EpiExp.m, reduction = 'phate', group.by = 'age_state',label=F,pt.size=2,cols = bar)+ggtitle("Melanoma tumor cells")
ggplot(data = data.frame(DC1 = dm$DC1,DC2 = dm$DC2, State = dm$age_state)) +
  geom_point(aes(DC1, DC2,colour = State),size = 2)+
  scale_color_manual(values = bar)+
  theme_classic()+
  theme(legend.position = c(0.3,0.7),text = element_text(size =16))+
  ggtitle('Melanoma tumor cells')
dev.off()

### aging markers
DefaultAssay(EpiExp.m)='RNA'
SenMarkers = c("CDKN1A", "SERPINE1")
EpiExp.m[['diffusion']] <- CreateDimReducObject(data.frame(DC1=dm$DC1,DC2=dm$DC2,row.names = colnames(exp)) %>% as.matrix(),key="DC_")
DefaultAssay(EpiExp.m) <- 'RNA'

p1 = FeaturePlot(EpiExp.m,features=SenMarkers[1],reduction='diffusion',order=T,pt.size = 2) + xlim(-0.05,0.05) + ylim(-0.03,0.15)
p2 = FeaturePlot(EpiExp.m,features=SenMarkers[2],reduction='diffusion',order=T,pt.size = 2) + xlim(-0.05,0.05) + ylim(-0.03,0.15)
p3 = VlnPlot(EpiExp.m,features=SenMarkers[1],group.by='age_state',assay = 'RNA',cols = bar,pt.size = 0) + theme(axis.title.x = element_blank())
p4 = VlnPlot(EpiExp.m,features=SenMarkers[2],group.by='age_state',assay = 'RNA',cols = bar,pt.size = 0) + theme(axis.title.x = element_blank())

png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/Melanoma_markers.png',width = 2500,height = 2000,res = 300)
ggarrange(p1,p3,p2,p4,ncol = 2,nrow = 2,legend = "none") 
dev.off()

### microarray gene overlap
colorList = list()
colorList[['Cycling_up']] <- c("#e63946","#2a9d8f")
colorList[['Cycling_down']] <- c("#007E99","#2a9d8f")
colorList[['Transitional_up']] <- c("#e63946","#e9c46a")
colorList[['Transitional_down']] <- c("#007E99","#e9c46a")
colorList[['Senescent_up']] <- c("#e63946","#e76f51")
colorList[['Senescent_down']] <- c("#007E99","#e76f51")

pList=list()
library(ggvenn)
for(set in names(enrich_reList)){
    pList[[set]] <- enrich_geneList[[set]] %>% 
                    ggvenn(show_percentage = T,
                            show_elements = F,
                            label_sep = ",",
                            digits = 1,
                            set_name_size = 0,
                            stroke_color = "white",
                            fill_color = colorList[[set]],
                            text_size = 4,fill_alpha = 0.8)+
                    labs(title = paste("Padj value:",as.character(signif(enrich_reList[[set]],2))))+
                    theme(plot.title = element_text(hjust = 0.5,vjust = 0,size = 15))+
                    scale_y_continuous(limits = c(-1, 1))
    
}
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/Melanoma_enrich.png',width = 1500,height = 600,res = 160)
ggarrange(pList$Cycling_up,
            pList$Transitional_up,
            pList$Senescent_up,
            pList$Cycling_down,
            pList$Transitional_down,
            pList$Senescent_down,
            ncol = 3,nrow = 2,legend = "none")
dev.off()

### marker enrichment plot
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/Melanoma_state_enrich.png',width = 2000,height = 1500,res = 300)
dotplot(go,showCategory = 5) + theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5))+
  ggtitle("GO:BP Enrichment of State Marker Genes")+theme(plot.title = element_text(hjust = 1))
dev.off()

### draw heatmap of SKCM CIBERSORT
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/TCGA_SKCM_state_immue_cor.png',width = 1800,height = 1200,res= 300)
bk <- c(seq(-0.3,0.3,by=0.01))
pheatmap::pheatmap(cor_mat,show_colnames = T,show_rownames = T,
                   cluster_rows = F,cluster_cols = T,fontsize = 12,
                   border_color = "white",color = colorRampPalette(c("blue", "#f5f4f4", "red"))(length(bk)),
                   legend_breaks=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3),breaks=bk)
dev.off()

### survival plot
plotsurv <- function(myfit){
    p <- ggsurvplot(
    myfit,
    risk.table = T,
    pval = TRUE,
    conf.int = TRUE,
    xlim = c(0,4000),
    break.time.by = 1000,
    risk.table.y.text.col = T,
    risk.table.y.text = FALSE)
    return(p)
}

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'Cycling')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ Cycling,data = dat)
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/TCGA_SKCM_survival_Cycling.png',width = 1200,height = 1600,res = 300)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'Transitional')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ Transitional,data = dat)
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/TCGA_SKCM_survival_Transitional.png',width = 1200,height = 1600,res = 300)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'Senescent')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ Senescent,data = dat)
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/TCGA_SKCM_survival_Senescent.png',width = 1200,height = 1600,res = 300)
plotsurv(fit)
dev.off()

### all melanome cell plot
mypalette <- read.csv("~/scripts/colors.csv",header = T)
bar = mypalette$palette1[1:length(unique(melanoma_obj$subtype))]
names(bar) <- unique(melanoma_obj$subtype)

png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/Melanoma_all_cell.png',width = 2000,height = 1500,res = 300)
DimPlot(melanoma_obj, group.by = 'subtype',reduction = "tsne",label = T,
        cols = bar,repel = TRUE,label.box = T,label.color = "white",pt.size = 1.5,label.size = 6)+
  ggtitle("Melanoma cell types")
dev.off()

### cell chat strength count
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/melanome_cellchat_scatter.png',width = 1200,height = 1000,res = 300)
source('/home/wangjing/wangj/codebase/HUSI/netAnalysis_signalingRole_scatter.R')
netAnalysis_signalingRole_scatter_log(cellchat,color.use = bar,dot.alpha = 0)
dev.off()

### cell chat strength cord
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/melanome_cellchat_network.png',width = 1600,height = 1200,res = 300)
mat <- log10(cellchat@net$weight)
par(mfrow = c(2,3), xpd=TRUE,mar = c(0.2, 0.2, 0.2,0.2))
for (i in c(3,9,7)) {
  mat_plot <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat_plot[i, ] <- mat[i, ]
  netVisual_circle(mat_plot, vertex.weight = groupSize, weight.scale = T,arrow.size = 1,color.use = bar,edge.width.max = 5,title.name = rownames(mat)[i])
}
for (i in c(3,9,7)) {
  mat_plot <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat_plot[, i] <- mat[, i]
  netVisual_circle(mat_plot, vertex.weight = groupSize, weight.scale = T,arrow.size = 1,color.use = bar,edge.width.max = 5,title.name = rownames(mat)[i])
}
dev.off()

### different pattern
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/melanome_cellchat_pattern.png',width = 2000,height = 1500,res = 300)
netAnalysis_river(cellchat, pattern = "outgoing",cutoff = 0.5,sources.use = c('Cycling','Transitional','Senescent'),signaling = as.character(pattern_use$Signaling))
dev.off()

### pathway strength bubble
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/melanome_cellchat_pathway.png',width = 1500,height = 1200,res = 300)
source('/home/wangjing/wangj/codebase/HUSI/netVisual_bubble_my.R')
netVisual_bubble_my(cellchat,signaling = pathways,targets.use = c('Cycling','Transitional','Senescent'), remove.isolate = F,sources.use = c('T cell','NK cell','Macro cell','CAF cell'),thresh=0.01)
dev.off()

### pathway strength cord
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/melanome_cellchat_chord.png',width = 2000,height = 1000,res = 300)
par(mfrow = c(1,2), xpd=TRUE,mar = c(0.2, 0.2, 0.2,0.2))
for(p in c('BMP','TGFb')){
    netVisual_aggregate(cellchat, signaling = p,layout = "chord",remove.isolate = F,color.use = bar,sources.use = c('T cell','NK cell','Macro cell','CAF cell'),targets.use = c('Cycling','Transitional','Senescent'))
}
dev.off()

### pathway gene expression
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/melanome_cellchat_expression.png',width = 900,height = 1000,res = 300)
plotGeneExpression(cellchat, features=c('BMPR1B','BMPR2','TGFBR1','TGFBR2'),idents = c('Cycling','Transitional','Senescent'),color.use = bar)
dev.off()

### receptor survival plot
plotsurv <- function(myfit){
  p <- ggsurvplot(
    myfit,
    risk.table = F,
    pval = TRUE,
    conf.int = TRUE,
    xlim = c(0,4000),
    break.time.by = 1000,
    risk.table.y.text.col = T,
    risk.table.y.text = FALSE)
  return(p)
}

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'BMPR2')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ BMPR2,data = dat)
png('/home/wangjing/wangj/codebase/HUSI/Figures/TCGA_SKCM_survival_BMPR2.png',width = 1200,height = 1000,res = 300)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'BMPR1B')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ BMPR1B,data = dat)
png('/home/wangjing/wangj/codebase/HUSI/Figures/TCGA_SKCM_survival_BMPR1B.png',width = 1200,height = 1000,res = 300)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'TGFBR1')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ TGFBR1,data = dat)
png('/home/wangjing/wangj/codebase/HUSI/Figures/TCGA_SKCM_survival_TGFBR1.png',width = 1200,height = 1000,res = 300)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'TGFBR2')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ TGFBR2,data = dat)
png('/home/wangjing/wangj/codebase/HUSI/Figures/TCGA_SKCM_survival_TGFBR2.png',width = 1200,height = 1000,res = 300)
plotsurv(fit)
dev.off()