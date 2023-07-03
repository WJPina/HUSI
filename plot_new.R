library(ggplot2)
library(dplyr)
library(tibble)
library(data.table)
library(annaffy)
library(ggsignif)
library(stringi)
library(ggpubr)
library(ggsci)

mytheme <- function () { 
    theme_classic() %+replace% 
    theme(text = element_text(size = 16),
            axis.ticks.x = element_blank(),
            legend.position = "right",
            axis.text.x = element_blank())
    }

### model leave-one-out auc
load("/home/wangjing/wangj/AgingScore/BulkData/Bulk_TrainModel/model_auc.RData")

png('/home/wangjing/codebase/HUSI/Figures/auc.png',width = 1000,height = 800,res = 200)
ggplot(aes(x = auc),data=data.frame(auc=auc))+
    geom_density()+
    theme_classic()+
    theme(text=element_text(size=16))+
    xlab("Leave-one-out cross-validation auc")+
    ylab("Density")
dev.off()

### valid data
## new celltype
load("/home/wangjing/wangj/AgingScore/AgingScorePro/Data1_Scripts/ModelValidData.RData")
load("/home/wangjing/wangj/AgingScore/BulkData/Bulk_Microarray/Valid.RData")

p1.1 <- data.frame(sene_score = s_IS) %>% 
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
    mytheme()+
    ylab('hUSI')+
    xlab("DIS of MDAH041")
### new platform
p1.2 <- ArrayList [["GSE16058"]][[1]] %>% 
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
    xlab("RS of HMEC")

fig1 <- ggarrange(p1.1,p1.2,ncol = 2,nrow = 1,common.legend = F)
png('/home/wangjing/codebase/HUSI/Figures/valid_new.png',width = 2000,height = 800,res = 200)
fig1
dev.off()

### bacth effect
load('/home/wangjing/wangj/AgingScore/AgingScorePro/Data1_Scripts/ModelValidData_Batch.RData')

png('/home/wangjing/codebase/HUSI/Figures/valid_BE.png',width = 1500,height = 1200,res = 300)
read.table("/home/wangjing/wangj/AgingScore/BulkData/Bulk_BatchEffect/batch_IMR90_4OHT_add_condition_tsv.txt",sep = "\t",header = T) %>%
    dplyr::select(c("title","study_accession","condition")) %>%
    mutate(sene_score = s_batch[title]) %>%
    column_to_rownames("title") %>%
    mutate(condition = factor(condition,levels = c("other","senescent"),ordered = T)) %>%
    ggplot(aes(condition, sene_score)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(aes(color = study_accession),width = 0.1) +
    scale_colour_jama()+
    geom_line(aes(group = study_accession), linetype="dashed", col="skyblue") + 
    theme_classic()+
    geom_signif(comparisons = list(c("other","senescent")),test = "t.test",test.args = c("less"),map_signif_level = T) + 
    ylab('hUSI')+
    ggtitle("OIS of IMR90")+
    theme_classic() %+replace% theme(text = element_text(size = 16),axis.title.x=element_blank())
dev.off()

### GSEA
library(enrichplot)
library(gggsea)

load('/home/wangjing/wangj/AgingScore/AgingScorePro/Data1_Scripts/Model_GSEA.RData')
cols = c("#B71C1C", "#C62828", "#D32F2F", "#E53935", "#F44336", "#EF5350", "#E57373","#0D47A1", "#1565C0", "#1976D2", "#1E88E5")
names(cols) = df$ID
png('/home/wangjing/codebase/HUSI/Figures/valid_GSEA_DOWN.png',width = 2500,height = 1500,res = 200)
gseaplot2(fgsea,geneSetID = df$ID[8:11],
          title = "Negative enriched hallmarker gene sets",
          color= cols[8:11],
          base_size = 14,
          rel_heights = c(1, 0.2, 0.4),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line"
)
dev.off()



