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


### comparision auc
png('/home/wangjing/codebase/HUSI/Figures/compare_Tang2019_ssgsea.png',width = 1800,height = 1000,res = 200)
auc %>%
    melt(value.name = "Accuracy") %>%
    mutate(feature=as.character(Var1)) %>% 
    mutate(feature=forcats::fct_reorder(feature, Accuracy, mean, .desc = T) ) %>% 
    ggplot(aes(feature, Accuracy)) + 
    geom_boxplot(aes(color=feature), outlier.shape = NA, width=0.5) + 
    geom_jitter(aes(color=feature), width = 0.2) + 
    labs(x=NULL, y="AUC") +
    theme(legend.position = "none")+
    theme_classic()
dev.off()

##### menanome
bar = c("#5cb85c","#428bca","#d9534f")
names(bar) = c("Cycling","Moderate_senescent","Senescent")
### Phate trajectory
png('/home/wangjing/codebase/HUSI/Figures/Melanoma_phate.png',width = 1800,height = 800,res = 200)
DimPlot(EpiExp.m, reduction = 'phate', group.by = 'age_state',label=F,pt.size=2,cols = bar)+ggtitle("Melanoma tumor cells")
dev.off()

### aging markers
p1 = FeaturePlot(EpiExp.m,features=SenMarkers[1],reduction='phate',order=T,pt.size = 2) + xlim(-0.04,0.05) + ylim(-0.03,0.04)
p2 = FeaturePlot(EpiExp.m,features=SenMarkers[2],reduction='phate',order=T,pt.size = 2) + xlim(-0.04,0.05) + ylim(-0.03,0.04)
p3 = VlnPlot(EpiExp.m,features=SenMarkers[1],group.by='age_state',assay = 'RNA',cols = bar,pt.size = 0) + theme(axis.title.x = element_blank())
p4 = VlnPlot(EpiExp.m,features=SenMarkers[2],group.by='age_state',assay = 'RNA',cols = bar,pt.size = 0) + theme(axis.title.x = element_blank())
png('/home/wangjing/codebase/HUSI/Figures/Melanoma_markers.png',width = 1800,height = 1200,res = 200)
ggarrange(p1,p3,p2,p4,ncol = 2,nrow = 2,legend = "none") 
dev.off()

### microarray gene overlap
colorList = list()
colorList[['Cycling_up']] <- c("#FF745A","#5cb85c")
colorList[['Cycling_down']] <- c("#007E99","#5cb85c")
colorList[['Moderate_senescent_up']] <- c("#FF745A","#428bca")
colorList[['Moderate_senescent_down']] <- c("#007E99","#428bca")
colorList[['Senescent_up']] <- c("#FF745A","#d9534f")
colorList[['Senescent_down']] <- c("#007E99","#d9534f")

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
                            text_size = 4)+
                    labs(title = paste("Padj value:",as.character(signif(enrich_reList[[set]],2))))+
                    theme(plot.title = element_text(hjust = 0.5,vjust = 0,size = 12))+
                    scale_y_continuous(limits = c(-1, 1))
    
}
png('/home/wangjing/codebase/HUSI/Figures/Melanoma_enrich.png',width = 1500,height = 1500,res = 200)
ggarrange(pList$Cycling_up,pList$Cycling_down,
            pList[['Moderate_senescent_up']],pList[['Moderate_senescent_down']],
            pList$Senescent_up,pList$Senescent_down,
            ncol = 2,nrow = 3,legend = "none")
dev.off()