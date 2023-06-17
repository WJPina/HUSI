library(ggplot2)
library(dplyr)
library(tibble)
library(data.table)
library(annaffy)
library(ggsignif)
library(stringi)
library(ggpubr)

mytheme <- function () { 
    theme_classic() %+replace% 
    theme(text = element_text(size = 16),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "right",
            axis.text.x = element_blank())
    }

# Fig 3.1
setwd("~/wangj/AgingScore/AgingScorePro/Data1_Scripts")
load("ModelValidData.RData")
# Fig 3.1a
p3.1a <- data.frame(sene_score = s_IS) %>% 
    mutate(condition= gsub("_.*", "", names(s_IS))) %>%
    .[which(.$condition %in% c("Immortal","Adria","H2O2","5-aza")),] %>%
    mutate(condition = factor(condition,levels = c("Immortal","Adria","H2O2","5-aza"),ordered = T)) %>% 
    ggplot(aes(condition, sene_score)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(aes(color = condition),width = 0.1)+
    geom_signif(comparisons = list(c("Immortal","Adria"),c("Immortal","H2O2"),c("Immortal","5-aza")),
                test = "t.test",
                step_increase=0.1,
                map_signif_level = T,
                test.args = c("less")) + 
    mythem()+
    ylab('hUSI')
# Fig 3.1b
p3.1b <- data.frame(sene_score = s_RS) %>%
    cbind(condition = unlist(sapply(strsplit(names(s_RS),"_"),"[",1))) %>%
    .[which(!.$condition %in% c("RS","OISD5")),] %>%
    mutate(condition = factor(condition,levels = c("OISD0", "OISD2","OISD4","OISD6","OISD10"),ordered = T)) %>%
    ggplot(aes(condition, sene_score)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(aes(color = condition),width = 0.1)+
    geom_signif(comparisons = list(c("OISD0", "OISD2"),c("OISD0","OISD4"),c("OISD0","OISD6"),c("OISD0","OISD10")),
                test = "t.test",          
                step_increase=0.1,
                map_signif_level = T,
                test.args = c("less")) + 
    mythem()+
    ylab('hUSI')
# Fig 3.1c
p3.1c <- data.frame(sene_score = s_RS) %>%
    rownames_to_column("smp") %>% 
    mutate(condition=gsub("_Rep[0-9]", "", smp)) %>%
    .[which(.$condition %in% c("RS_Proliferative","RS_Senescent")),] %>%
    mutate(condition=gsub("RS_", "", condition)) %>%
    ggplot(aes(condition, sene_score)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(aes(color = condition),width = 0.1)+
    geom_signif(comparisons = list(c("Proliferative","Senescent")),
                test = "t.test",
                step_increase=0.1,
                map_signif_level = T,
                test.args = c("less")) + 
    mythem()+           
    ylab('hUSI')

fig3.1 <- ggarrange(p3.1a,p3.1b,p3.1c,ncol = 3,nrow = 1,common.legend = F)+ theme(plot.margin = unit(c(1,1,1,1), "cm"))

png('Figures/Valid_1.png',width = 3500,height = 1200,res = 300)
fig3.1
dev.off()
# Fig 3.2
names(ArrayList)
# Fig 3.2a GSE19864
p3.2a <- ArrayList [["GSE19864"]][[1]] %>% 
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
    ylab('hUSI')
# Fig 3.2b GSE16058
p3.2b <- ArrayList [["GSE16058"]][[1]] %>% 
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
    ylab('hUSI')
# Fig 3.2c GSE83922
p3.2c <- ArrayList [["GSE83922"]][[1]] %>% 
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
    ylab('hUSI')
# Fig 3.2d GSE11954
p3.2d <- ArrayList [["GSE11954"]][[1]] %>% 
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
    ylab('hUSI')
# Fig 3.2e GSE100014
p3.2e <- ArrayList [["GSE100014"]][[1]] %>% 
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
    ylab('hUSI')
# Fig 3.2f GSE77239
p3.2f <- ArrayList [["GSE77239"]][[1]] %>% 
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
    ylab('hUSI')

fig3.2 <- ggarrange(p3.2a,p3.2b,p3.2c,p3.2d,p3.2e,p3.2f,ncol = 3,nrow = 2,common.legend = F)+ theme(plot.margin = unit(c(1,1,1,1), "cm"))

png('Figures/Valid_2.png',width = 3500,height = 2500,res = 300)
fig3.2
dev.off()
# Fig 3.3
load('ModelValidData_Batch.RData')
png('Figures/Valid_3.png',width = 1500,height = 1000,res = 300)
read.table("/home/wangjing/wangj/AgingScore/Data1/Bulk_BatchEffect/batch_IMR90_4OHT_add_condition_tsv.txt",sep = "\t",header = T) %>%
    dplyr::select(c("title","study_accession","condition")) %>%
    mutate(sene_score = s_batch[title]) %>%
    column_to_rownames("title") %>%
    mutate(condition = factor(condition,levels = c("other","senescent"),ordered = T)) %>%
    ggplot(aes(condition, sene_score)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(aes(color = study_accession),width = 0.1) +
    geom_line(aes(group = study_accession), linetype="dashed", col="skyblue") + 
    theme_classic()+
    geom_signif(comparisons = list(c("other","senescent")),test = "t.test",test.args = c("less"),map_signif_level = T) + 
    mythem()+
    ylab('hUSI')
dev.off()

### Fig 3.6
load('Model_GSEA.RData')
png('Figures/Valid_4.png',width = 1500,height = 1000,res = 200)
plotGseaTable(pathways_sene[res_fgsea_sene[padj<0.05][order(NES,decreasing = T),pathway]], mm_l2$w, res_fgsea_sene, render = F ) %>% patchwork::wrap_elements()
dev.off()

### Fig3.17
setwd('/home/wangjing/wangj/AgingScore/Data1/scRNA_melanome')
load('tumor.RData')

png('/home/wangjing/wangj/AgingScore/AgingScorePro/Data1_Scripts/Figures/melanoma_1.png',width = 1200,height = 800,res= 200)
DimPlot(EpiExp.m, reduction = 'phate', group.by = 'age_class',label=1,pt.size=2)
dev.off()

png('/home/wangjing/wangj/AgingScore/AgingScorePro/Data1_Scripts/Figures/melanoma_2.png',width = 800,height = 400,res= 200)
VlnPlot(EpiExp.m,features = c('CDKN1A','SERPINE1'), group.by = 'age_class',pt.size=0,assay = 'RNA')
dev.off()
### Fig3.18

library(ComplexHeatmap)
library(circlize)
top_anno <- HeatmapAnnotation(df = data.frame(Condition = meta$condition),
                              col = list(Condition = c('young'='#007E99','senescent' = '#FF745A')),
                              show_legend = F)

left_anno = rowAnnotation(df = data.frame('State'= rep(c('Cycling Cell','Median Senescent Cell','Senescent Cell'),
                                                        times=lapply(gene_set,length))),
                          show_legend = T)

col <- colorRamp2(c(-2,0,2), c("blue","white", "red"), space = "LAB")

png('/home/wangjing/wangj/AgingScore/AgingScorePro/Data1_Scripts/Figures/melanoma_3.png',width = 1000,height = 800,res= 200)
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
        column_split = c(rep('young',4),rep('senescent',4)),
        row_split = rep(c('Cycling Cell','Median Senescent Cell','Senescent Cell'),times=lapply(gene_set,length)))
dev.off()