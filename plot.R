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
library(reshape2)
library(rstatix)

##################### senescent marker expression in train set
genes = c('CDKN1A','CDKN2B')
df_plot = data.frame(t(cbind(X_bk[genes,],X_tr[genes,])),Condition = c(rep('other',122),rep('senescent',126))) %>%
            melt(value.name = 'expression',variable.name = 'marker',id.var = 'Condition') 
pwc <- df_plot %>% group_by(marker) %>% 
        rstatix::pairwise_t_test(expression ~ Condition, paired = F,p.adjust.method = "bonferroni")
pwc
pwc <- pwc %>% rstatix::add_xy_position(x = "marker")
png('/home/wangjing/wangj/codebase/HUSI/Figures/model/valid_markers.png',width = 1200,height = 1500,res = 300)
df_plot %>%
  ggplot(aes(x = marker, y = expression)) + 
  geom_boxplot(aes(color = Condition),outlier.shape = NA) + 
  geom_point(aes(color = Condition),position = position_jitterdodge(jitter.width = 0.2),size = 1) +
  stat_pvalue_manual(pwc)+
  theme_classic()+
  theme(text = element_text(size = 16),axis.text.x = element_text(size=16),legend.position = 'bottom',axis.title.x = element_blank())+
  scale_color_manual(values = c('senescent' = '#e63946','other'='#023e8a'))+
  ylab('Normalized Expression')
dev.off()

############################ GSEA
library(enrichplot)
library(gggsea)
source('~/wangj/codebase/HUSI/mygseaplot2.R')
cols = c(rev(colorRampPalette(c("transparent", 'red'))(6))[1:5],colorRampPalette(c("transparent", 'blue'))(5)[c(2:5)])
names(cols) = df$ID
png('/home/wangjing/wangj/codebase/HUSI/Figures/model/valid_GSEA_DOWN.png',width = 1600,height = 1500,res = 300)
mygseaplot2(fgsea,
            geneSetID = df$ID[6:9],
            # geneSetID = df$ID[1:5],
            title = "Negatively enriched gene sets",
            # title = "Positively enriched gene sets",
            color= cols[6:9],
            # color = cols[1:5],
            base_size = 12,
            rel_heights = c(1, 0.2, 0.4),
            subplots = 1:3,
            pvalue_table = FALSE,
            ES_geom = "line"
)
dev.off()

############################### plot train Rep.CS heatmap
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
top_anno <- HeatmapAnnotation(df = data.frame(Condition = groups),show_legend = T,col = list(Condition = bar),
                              show_annotation_name = F)
left_anno = rowAnnotation(df = data.frame(Signature= rep(c("up","down"),times=c(519,729))),show_legend = T,
                          col = list(Signature = c('up' = '#e63946','down' = '#023e8a')),show_annotation_name = F)
col <- colorRamp2(c(-1.5,0,1.5), c("#023e8a","white", "#e63946"), space = "LAB")

png('/home/wangjing/wangj/codebase/HUSI/Figures/model/valid_CS_sinature.png',width = 2200,height = 1500,res= 300)
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
        heatmap_legend_param = list(title = "Scaled\nexpression"),
        row_split = factor(rep(c("up","down"),times=c(519,729)),levels = c('up','down'),ordered = T),
        column_split = factor(rep(c('Senescent','other'),times = c(126,122)),levels=c('Senescent','other'),ordered = T))
dev.off()

########################### model leave-one-out auc
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

################################## bacth effect
png('/home/wangjing/wangj/codebase/HUSI/Figures/model/valid_BE.png',width = 1500,height = 1500,res = 300)
read.table("/home/wangjing/wangj/AgingScore/Data/Bulk_BatchEffect/batch_IMR90_4OHT_add_condition_tsv.txt",sep = "\t",header = T) %>%
  dplyr::select(c("title","study_accession","condition")) %>%
  mutate(hUSI = s_batch[title]) %>%
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


################################## valid data in RNA-seq
mytheme <- function () { 
  theme_classic() %+replace% 
    theme(text = element_text(size = 16),legend.position = "none",
          axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))}

p1.1 <- data.frame(hUSI = s_IS) %>% 
  mutate(condition= gsub("_.*", "", names(s_IS))) %>%
  .[which(.$condition %in% c("Immortal","Adria","H2O2","5-aza")),] %>%
  mutate(Treatment = factor(condition,levels = c("Immortal","Adria","H2O2","5-aza"),ordered = T)) %>% 
  ggplot(aes(Treatment, hUSI)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Treatment),width = 0.1)+
  geom_signif(comparisons = list(c("Immortal","H2O2"),c("Immortal","5-aza")),
              test = "t.test",
              step_increase=0.1,
              map_signif_level = T,
              test.args = c("less")) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle("CIS in MDAH041")

p1.2<- data.frame(hUSI = s_RS[grepl('OISD[0|2|4|6|10]',names(s_RS))]) %>% 
  mutate(condition= gsub("_.*", "", rownames(.))) %>%
  mutate(Treatment = factor(condition,levels = c("OISD0","OISD2","OISD4","OISD6",'OISD10'),ordered = T)) %>% 
  ggplot(aes(Treatment, hUSI)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Treatment),width = 0.1)+
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
  mutate(Condition = factor(condition,levels = c("Proliferative","Senescent"),ordered = T)) %>% 
  ggplot(aes(Condition, hUSI)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Condition),width = 0.1)+
  geom_signif(comparisons = list(c("Proliferative","Senescent")),
              test = "t.test",
              step_increase=0.1,
              map_signif_level = T,
              test.args = c("less")) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle("RS in WI-38")

fig <- ggarrange(p1.1,p1.2,p1.3,ncol = 3,nrow = 1,common.legend = F,widths = c(4,4,3))
png('/home/wangjing/wangj/codebase/HUSI/Figures/model/valid_RNA-seq.png',width = 3000,height = 1500,res = 300)
fig
dev.off()
############################# valid data in micro-array
mytheme <- function () { 
  theme_classic() %+replace% 
    theme(text = element_text(size = 16),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "right",
          axis.text.x = element_blank())
}

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
  mutate(Condition = factor(condition,levels = c("Growing","Senescent"),ordered = T)) %>% 
  inner_join(
    data.frame(sene_score=ScoreList[["s_19864"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  ggplot(aes(Condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Condition),width = 0.1)+
  geom_signif(comparisons = list(c("Growing", "Senescent")),test = "t.test",test.args = c("less"), map_signif_level = T) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle("OIS in IMR90")+
  scale_colour_manual(values = c('#023e8a','#e63946'))

p2.2 <- ArrayList [["GSE16058"]][[1]] %>% 
  pData %>% 
  mutate(condition=gsub("growth status: ", "", characteristics_ch1.3) ) %>% 
  mutate(passage = `passage:ch1`) %>%
  inner_join(
    data.frame(sene_score=ScoreList[["s_16058"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("Growing","Senescent")),] %>% 
  mutate(Condition = factor(condition,levels = c("Growing","Senescent"),ordered = T)) %>% 
  mutate(Passage = passage) %>% 
  ggplot(aes(Condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Condition,shape = Passage),width = 0.1) +
  geom_signif(comparisons = list(c("Growing", "Senescent")),test = "t.test",test.args = c("less"), map_signif_level = T) + 
  mytheme()+
  guides(color = 'none')+
  ylab('hUSI')+
  ggtitle('RS in HMEC')+
  scale_colour_manual(values = c('#023e8a','#e63946'))

p2.3 <- ArrayList [["GSE83922"]][[1]] %>% 
  pData %>% 
  mutate(condition= `cell phenotype:ch1`) %>% 
  inner_join(
    data.frame(sene_score=ScoreList[["s_83922"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("young","senescent")),] %>% 
  mutate(Condition = factor(ifelse(condition == 'young','Young','Senescent'),levels = c("Young","Senescent"),ordered = T)) %>%
  ggplot(aes(Condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Condition),width = 0.1) +
  geom_signif(comparisons = list(c("Young", "Senescent")),test = "t.test",test.args = c("less"), map_signif_level = T) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle('RS in Melanocyte')+
  scale_colour_manual(values = c('#023e8a','#e63946'))

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
  mutate(Condition = factor(condition,levels = c("Growing","Senescent"),ordered = T)) %>% 
  ggplot(aes(Condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Condition),width = 0.1) +
  geom_signif(comparisons = list(c("Growing","Senescent")),test = "t.test", map_signif_level = T) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle('CIS in HSC')+
  scale_colour_manual(values = c('#023e8a','#e63946'))

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
  mutate(Condition = factor(condition,levels = c("Proliferating","Senescent"),ordered = T)) %>% 
  ggplot(aes(Condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Condition),width = 0.1) +
  geom_signif(comparisons = list(c("Proliferating","Senescent")),test = "t.test",test.args = c("less"), map_signif_level = T) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle('CIS in HBEC')+
  scale_colour_manual(values = c('#023e8a','#e63946'))

p2.6 <- ArrayList [["GSE77239"]][[1]] %>% 
  pData %>% 
  mutate(condition=case_when(
    grepl("young", `cells:ch1`, ignore.case = T) ~ "young", 
    grepl("old", `cells:ch1`, ignore.case = T) ~ "old")) %>% 
  mutate(Condition = factor(condition,levels = c("young","old"),ordered = T)) %>% 
  rename("Treatment"="treatment:ch1") %>% 
  inner_join(
    data.frame(sene_score=ScoreList[["s_77239"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  ggplot(aes(Condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Condition,shape = Treatment),width = 0.1) +
  geom_line(aes(group = Treatment), linetype="dashed", col="skyblue") + 
  geom_signif(comparisons = list(c("young","old")),test = "t.test",test.args = c("less"), map_signif_level = T) + 
  mytheme()+
  guides(color = 'none')+
  ylab('hUSI')+
  ggtitle('RS in HCAEC')+
  scale_colour_manual(values = c('#023e8a','#e63946'))

fig <- ggarrange(p2.1,p2.2,p2.3,p2.4,p2.5,p2.6,ncol = 3,nrow = 2,common.legend = F)+ 
        theme(plot.margin = unit(c(0,0,0,0), "cm"))
png('/home/wangjing/wangj/codebase/HUSI/Figures/model/valid_array.png',width = 4500,height = 2400,res = 300)
fig
dev.off()

################################# plot lost rate stability
load('/home/wangjing/wangj/AgingScore/Data/Bulk_TrainModel/Valid_lost_1013.RData')
mytheme <- function () { 
  theme_classic() %+replace% 
    theme(text = element_text(size = 16))
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
  mutate(condition = factor(ifelse(condition == 'Growing','Other','Senescent'),levels = c("Other","Senescent"),ordered = T)) %>% 
  inner_join(
    data.frame(hUSI=Lost_ScoreList[["GSE19864"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>%
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  mutate(Condition = condition) %>%
  ggplot(aes(x = rate,y = hUSI,color = Condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle("OIS in IMR90")+
  scale_colour_manual(values = c('#023e8a','#e63946'))+
  xlab('lost rate')
  

p2.2 <- ArrayList [["GSE16058"]][[1]] %>% 
  pData %>% 
  mutate(condition=gsub("growth status: ", "", characteristics_ch1.3) ) %>% 
  mutate(passage = `passage:ch1`) %>%
  inner_join(
    data.frame(hUSI=Lost_ScoreList[["GSE16058"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("Growing","Senescent")),] %>% 
  mutate(condition = factor(ifelse(condition == 'Growing','Other','Senescent'),levels = c("Other","Senescent"),ordered = T)) %>%
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  mutate(Condition = condition) %>%
  ggplot(aes(x = rate,y = hUSI,color = Condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle('RS in HMEC')+
  scale_colour_manual(values = c('#023e8a','#e63946'))+
  xlab('lost rate')

p2.3 <- ArrayList [["GSE83922"]][[1]] %>% 
  pData %>% 
  mutate(condition= `cell phenotype:ch1`) %>% 
  inner_join(
    data.frame(hUSI=Lost_ScoreList[["GSE83922"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("young","senescent")),] %>% 
  mutate(condition = factor(ifelse(condition == "young","Other","Senescent"),levels = c("Other","Senescent"),ordered = T)) %>%
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  mutate(Condition = condition) %>%
  ggplot(aes(x = rate,y = hUSI,color = Condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle('RS in Melanocyte')+
  scale_colour_manual(values = c('#023e8a','#e63946'))+
  xlab('lost rate')


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
  mutate(condition = factor(ifelse(condition == 'Growing','Other','Senescent'),levels = c("Other","Senescent"),ordered = T)) %>% 
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  mutate(Condition = condition) %>%
  ggplot(aes(x = rate,y = hUSI,color = Condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle('CIS in HSC')+
  scale_colour_manual(values = c('#023e8a','#e63946'))+
  xlab('lost rate')


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
  mutate(condition = factor(ifelse(condition == 'Proliferating',"Other","Senescent"),levels = c("Other","Senescent"),ordered = T)) %>%
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  mutate(Condition = condition) %>%
  ggplot(aes(x = rate,y = hUSI,color = Condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle('CIS in HBEC')+
  scale_colour_manual(values = c('#023e8a','#e63946'))+
  xlab('lost rate')


p2.6 <- ArrayList [["GSE77239"]][[1]] %>% 
  pData %>% 
  mutate(condition=case_when(
    grepl("young", `cells:ch1`, ignore.case = T) ~ "young", 
    grepl("old", `cells:ch1`, ignore.case = T) ~ "old")) %>% 
  mutate(condition = factor(ifelse(condition == 'young',"Other","Senescent"),levels = c("Other","Senescent"),ordered = T)) %>% 
  rename("treatment"="treatment:ch1") %>% 
  inner_join(
    data.frame(hUSI=Lost_ScoreList[["GSE77239"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>%
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  mutate(Condition = condition) %>%
  ggplot(aes(x = rate,y = hUSI,color = Condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle('RS in HCAEC')+
  scale_colour_manual(values = c('#023e8a','#e63946'))+
  xlab('lost rate')

fig <- ggarrange(p2.1,p2.2,p2.3,p2.4,p2.5,p2.6,ncol = 3,nrow = 2,common.legend = T,legend = 'bottom')
png('/home/wangjing/wangj/codebase/HUSI/Figures/model/valid_lost.png',width = 4500,height = 2400,res = 300)
fig
dev.off()

############################################# comparison auc
load('/home/wangjing/wangj/AgingScore/Comparison/compare.RData')
### value
df_plot_list = lapply(Results[c("Aarts2017","Georgilis2018")], function(x){data.frame(hUSI = x$hUSI,Condition = x$method$Condition)})
df_plot_list[['Teo2019']] = data.frame(hUSI = Results$Teo2019$hUSI,Condition = gsub("[0-9]$", "", Results$Teo2019$method$Condition1))
df_plot_list[["Tang2019"]] <- data.frame(hUSI=unlist(lapply(Results$Tang2019,function(x){x$method$hUSI})),Condition = rep(c('Growing','Senescence'),each=400,time = 2))
df_plot_list <- lapply(names(df_plot_list), function(x){df_plot_list[[x]]$Dataset = rep(x,nrow(df_plot_list[[x]]));return(df_plot_list[[x]])})
df_plot = data.table::rbindlist(df_plot_list) %>% data.frame()
pwc <- df_plot %>% group_by(Dataset) %>% 
        rstatix::pairwise_t_test(hUSI ~ Condition, paired = F,p.adjust.method = "bonferroni")
pwc
pwc <- pwc %>% rstatix::add_xy_position(x = "Dataset")

png('/home/wangjing/wangj/codebase/HUSI/Figures/model/compare_hUSI.png',width = 1400,height = 1500,res = 300)
df_plot %>%
  mutate(Condition = factor(Condition,levels = c('Growing','Senescence'),ordered = T)) %>%
  ggplot(aes(x = Dataset,y=hUSI,color = Condition)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_point(aes(color = Condition),position = position_jitterdodge(jitter.width = 0.3),size = 0.5) +
    scale_color_manual(values = c('#023e8a','#e63946'))+
    stat_pvalue_manual(pwc)+
    theme_classic()+
    theme(legend.position = 'top',text = element_text(size = 16),axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
dev.off()



### rank 
df_plot_list = lapply(AUClist, function(x) {lapply(x, function(y) {rank(rowMeans(y))})}) 

data = data.frame()
for(type in c('marker','method','ssgsea')) {
  df_plot = lapply(df_plot_list, function(x) {x[[type]]})
  df_plot = do.call(cbind, lapply(lapply(df_plot, unlist), `length<-`, max(lengths(df_plot)))) %>% data.frame
  df_plot[is.na(df_plot)] = 0
  df_plot = reshape2::melt(as.matrix(df_plot),value.name = "AUC_rank")
  colnames(df_plot) = c('method','dataset','AUC_rank')
  aov.mean<-aggregate(df_plot$AUC_rank,by=list(df_plot$method),FUN=mean)
  aov.sd<-aggregate(df_plot$AUC_rank,by=list(df_plot$method),FUN=sd)
  aov<-data.frame(aov.mean,sd=aov.sd$x,type = rep(type,nrow(aov.mean)))
  data <- rbind(data,aov)
}

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}

data$group <- ifelse(data$type == 'marker','Expression Level',ifelse(data$type == 'method','Senecence Score','ssGSEA Score'))
data$color <- ifelse(data$Group.1 == 'hUSI','hUSI',ifelse(data$type == 'marker','Expression Level',ifelse(data$type == 'method','Senecence Score','ssGSEA Score')))

png('/home/wangjing/wangj/codebase/HUSI/Figures/model/compare_rank.png',width = 3500,height = 1500,res = 300)
ggplot(data=data,aes(x=reorder_within(Group.1,x,group), y=x,fill = color))+
    scale_x_reordered()+
    geom_bar(stat="identity",position="dodge")+
    geom_errorbar(aes(ymax=x+sd,ymin=ifelse(x-sd <0,0,x-sd)),position=position_dodge(0.9),width=0.15)+
    theme_classic()+
    ylab('AUC rank')+
    xlab('Benchmark')+
    scale_fill_manual(values = c('#94d2bd','#e63946','#fed9b7','#00afb9'))+
    ggforce::facet_row(vars(group), scales = 'free', space = 'free')+
    theme(legend.position = 'none',axis.text.x = element_text(angle = 45,vjust = 1, hjust=1),
          text=element_text(size=16))
dev.off()

################################# overlap gene set plot
genes = intersect(names(mm_l2$w),unique(unlist(EnrichSet)))
mat = sapply(EnrichSet, function(set){ifelse(genes %in% set,1,0)})
mat = cbind(hUSI=mm_l2$w[genes],mat)

library(ComplexHeatmap)
library(circlize)
color <- colorRamp2(c(-0.03,-0.02,-0.01,0,0.01,0.02,0.03), c('#012a4a','#014f86','#468faf','#e8e8e4','#e9c46a','#e76f51','#e63946'))
png('/home/wangjing/wangj/codebase/HUSI/Figures/model/compare_geneSet_overlap.png',width = 3000,height = 3000,res = 300)
circos.clear()
circos.par(gap.after = c(40))
circos.heatmap(mat, 
               col = color,
               dend.side = "inside",
               rownames.side = "none",
               track.height = 0.6,
               cell.border = 'white',
               cell.lwd = 0.1)
grid.draw(Legend(title = "hUSI",title_gp = gpar(fontsize = 16, fontface = "bold"),
                 title_position = "topcenter", col_fun = color,
                 legend_height = unit(4, "cm"),grid_width = unit(5, "mm")))
dev.off()

#################################### GTEx/TCGA
### gtex

### tcga
TCSER = read.csv('Senescence_quantification_TCGA.csv')
TCSER <- column_to_rownames(TCSER_tcga,'sample')
rownames(TCSER) <- gsub('-','.',rownames(TCSER))
df_plot <- scoreList_TCGA[c('TCSER','hUSI')] %>% data.frame() %>% 
            mutate(CancerType = TCSER[names(scoreList_TCGA$hUSI),'cancertype'])


png('/home/wangjing/wangj/codebase/HUSI/Figures/model/TCGA_TCSER.png',width = 1000,height = 900,res = 300)
  df_plot %>%
  ggscatter(x = "hUSI", y = "TCSER",
            color = 'CancerType',
            # color = '#f07167',
            add = "reg.line", 
            conf.int = TRUE,
            size = 1,
            add.params = list(color = "#6c757d"),
            ggtheme = theme_classic())+ 
  # stat_cor(method = "spearman",label.x = -.08, label.y = 1,color='black')+
  stat_cor(method = "spearman",label.x = -.09, label.y = 1.15,color='black')+
  ggtitle('7123 samples in GTEx')+
  ggtitle("10495 samples in TCGA")+
  theme_classic()+
  theme(text = element_text(size = 12))
dev.off()


############################### menanome #######################################
bar = c("#2a9d8f","#e9c46a","#e76f51")
names(bar) = c("Cycling","Transitional","Senescent")

### histogram
png("/home/wangjing/wangj/codebase/HUSI/Figures/melanome/Melanoma_hist.png",width = 1200,height = 1200,res = 300)
ggplot(EpiExp.m@meta.data)+
  geom_histogram(aes(x=hUSI),bins = 50,fill = 'white',color = '#adb5bd',linewidth = 0.4)+
  geom_density(aes(x=hUSI,fill = State),alpha = 0.6,color = NA)+
  scale_fill_manual(values = bar)+
  theme_classic()+
  theme(text = element_text(size = 14),legend.position = c(0.25,0.7))+
  ylab('Count')
dev.off()

### state hUSI plot
pwc = EpiExp.m@meta.data %>% 
        rstatix::pairwise_t_test(hUSI ~ State, paired = F,p.adjust.method = "bonferroni") %>%
        rstatix::add_xy_position(x = "State")
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/Melanoma_hUSI.png',width = 1200,height = 1200,res = 300)
EpiExp.m@meta.data %>%
  ggplot(aes(x = State,y = hUSI,color = State))+
  geom_boxplot()+
  scale_color_manual(values = bar)+
  # labs(caption = get_pwc_label(pwc))+
  stat_pvalue_manual(pwc)+
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.position = 'none')
dev.off()

### cor gene volcano
res <- data.frame(expression = rowMeans(EpiExp.m@assays$RNA@data[names(cor.genes),]),correlation = cor.genes,row.names = names(cor.genes))
res$threshold = ifelse(rownames(res) %in% features, ifelse(res$correlation > 0,"Positive","Negative"),"Background")
label_data <- subset(res, abs(res$correlation) > 0.5 & res$expression  > 4)
png("/home/wangjing/wangj/codebase/HUSI/Figures/melanome/Melanoma_corgenes.png",width = 1000,height = 1200,res = 300)
ggplot(data = res, aes(y=expression,x=correlation,color=threshold)) +
  geom_point(alpha=0.5, size=2)+
  geom_text_repel(data = label_data,label = rownames(label_data),color = 'black',max.overlaps = 10,size = 3,show.legend = F)+
  scale_color_manual(values = c("Positive" = "#FF745A","Negative"="#007E99","Background" = "#ABABAB"))+
  ylab('Mean normalized expression')+
  xlab('Pearson cofficient ')+
  theme_classic()+
  theme(text = element_text(size = 16),legend.position = 'none')
dev.off()


### icanet heatmap
library(ComplexHeatmap)
library(circlize)

root_cell = paste('DPT',grep(names(sort(dm$DC1+dm$DC2,decreasing = T)[1]),names(dm$DC1)),sep = '')
EpiExp.m$DPT <- dpt[[root_cell]]
groups = sort(EpiExp.m$DPT,decreasing = F)
top_anno <- HeatmapAnnotation(df = data.frame(DPT = groups,State = EpiExp.m$State[names(groups)]),show_legend = T,col = list(DPT = colorRamp2(c(0,4,8), c("#f2cc8f", "#f4f1de",'#81b29a')),State = bar))
col <- colorRamp2(c(-1,0,1), c("#023e8a","white", "#e63946"), space = "LAB")
DefaultAssay(EpiExp.m) <- 'IcaNet'
EpiExp.m <- ScaleData(EpiExp.m)
mat = EpiExp.m@assays$IcaNet@scale.data %>% as.matrix()
mat = mat[,names(groups)]

png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/Melanoma_icanet.png',width = 2500,height = 1800,res= 300)
Heatmap(mat,
        show_column_names = F,
        show_row_names = T,
        row_title = NULL,
        col = col,
        cluster_rows = T,
        cluster_row_slices = FALSE,
        cluster_columns = F,
        top_annotation = top_anno,
        row_names_gp = gpar(fontsize = 12),
        heatmap_legend_param = list(title = "Scaled\nexpression"))
dev.off()

### melanome trajectory
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/Melanoma_diffusion.png',width = 1200,height = 1000)
ggplot(data = data.frame(DC1 = dm$DC1,DC2 = dm$DC2, State = dm$State)) +
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
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/TCGA_SKCM_survival_Transitional.png',width = 1400,height = 1600,res = 300)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'Senescent')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ Senescent,data = dat)
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/TCGA_SKCM_survival_Senescent.png',width = 1200,height = 1600,res = 300)
plotsurv(fit)
dev.off()

### all melanome cell plot
bar = c("#2a9d8f","#e9c46a","#e76f51","#AAC5E2","#FFC592","#A4D99D","#F4A39E","#CCB7DC","#CBAEA8")
names(bar) <- c("Cycling","Transitional","Senescent","B cell","CAF cell","Endo cell","Macro cell","NK cell","T cell")

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
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/melanome_cellchat_pattern.png',width = 1900,height = 1500,res = 300)
netAnalysis_river(cellchat, pattern = "outgoing",cutoff = 0.5,
                  sources.use = c('Cycling','Transitional','Senescent'),
                  signaling = as.character(pattern_use$Signaling))
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
    netVisual_aggregate(cellchat, signaling = p,layout = "chord",remove.isolate = F,
                        color.use = bar,sources.use = c('T cell','NK cell','Macro cell','CAF cell'),
                        targets.use = c('Cycling','Transitional','Senescent'))
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