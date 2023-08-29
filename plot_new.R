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
            legend.position = "right",
            axis.text.x = element_blank())
    }


### model leave-one-out auc
load("/home/wangjing/wangj/AgingScore/BulkData/Bulk_TrainModel/model_auc.RData")

png('/home/wangjing/wangj/codebase/HUSI/Figures/auc.png',width = 1000,height = 800,res = 200)
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
png('/home/wangjing/wangj/codebase/HUSI/Figures/valid_new.png',width = 2000,height = 800,res = 200)
fig1
dev.off()

### bacth effect
load('/home/wangjing/wangj/AgingScore/AgingScorePro/Data1_Scripts/ModelValidData_Batch.RData')

png('/home/wangjing/wangj/codebase/HUSI/Figures/valid_BE.png',width = 1500,height = 1200,res = 300)
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
png('/home/wangjing/wangj/codebase/HUSI/Figures/valid_GSEA_DOWN.png',width = 2500,height = 1500,res = 200)
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
png('/home/wangjing/wangj/codebase/HUSI/Figures/compare_Tang2019_ssgsea.png',width = 1800,height = 1000,res = 200)
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
names(bar) = c("Cycling","Transition","Senescent")
### Phate trajectory
png('/home/wangjing/wangj/codebase/HUSI/Figures/Melanoma_phate.png',width = 1800,height = 800,res = 200)
DimPlot(EpiExp.m, reduction = 'phate', group.by = 'age_state',label=F,pt.size=2,cols = bar)+ggtitle("Melanoma tumor cells")
dev.off()

### aging markers
DefaultAssay(EpiExp.m)='RNA'
SenMarkers = c("CDKN1A", "SERPINE1")
p1 = FeaturePlot(EpiExp.m,features=SenMarkers[1],reduction='phate',order=T,pt.size = 2) + xlim(-0.04,0.05) + ylim(-0.03,0.04)
p2 = FeaturePlot(EpiExp.m,features=SenMarkers[2],reduction='phate',order=T,pt.size = 2) + xlim(-0.04,0.05) + ylim(-0.03,0.04)
p3 = VlnPlot(EpiExp.m,features=SenMarkers[1],group.by='age_state',assay = 'RNA',cols = bar,pt.size = 0) + theme(axis.title.x = element_blank())
p4 = VlnPlot(EpiExp.m,features=SenMarkers[2],group.by='age_state',assay = 'RNA',cols = bar,pt.size = 0) + theme(axis.title.x = element_blank())
png('/home/wangjing/wangj/codebase/HUSI/Figures/Melanoma_markers.png',width = 1800,height = 1200,res = 200)
ggarrange(p1,p3,p2,p4,ncol = 2,nrow = 2,legend = "none") 
dev.off()

### marker enrichment plot
png('/home/wangjing/wangj/codebase/HUSI/Figures/Melanoma_state_enrich.png',width = 1500,height = 1000,res = 200)
dotplot(go,showCategory = 5) + theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5))+
  ggtitle("GO:BP enrichment of state marker genes")+theme(plot.title = element_text(hjust = 1))
dev.off()

### microarray gene overlap
colorList = list()
colorList[['Cycling_up']] <- c("#FF745A","#5cb85c")
colorList[['Cycling_down']] <- c("#007E99","#5cb85c")
colorList[['Transition_up']] <- c("#FF745A","#428bca")
colorList[['Transition_down']] <- c("#007E99","#428bca")
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
                    theme(plot.title = element_text(hjust = 0.5,vjust = 0,size = 15))+
                    scale_y_continuous(limits = c(-1, 1))
    
}
png('/home/wangjing/wangj/codebase/HUSI/Figures/Melanoma_enrich.png',width = 1500,height = 600,res = 160)
ggarrange(pList$Cycling_up,
            pList$Transition_up,
            pList$Senescent_up,
            pList$Cycling_down,
            pList$Transition_down,
            pList$Senescent_down,
            ncol = 3,nrow = 2,legend = "none")
dev.off()

### age state fraction in TCGA SKCM patients
library(ComplexHeatmap)
library(circlize)
marker_plot = unlist(lapply(marker_set, function(x) {rownames(x[order(x$avg_log2FC,decreasing = T),][0:200,])})) %>% unique
genes = intersect(rownames(Bulk),marker_plot)
mat <- as.matrix(Bulk[genes,])
mat <- apply(mat,1,scale)
mat <- t(mat)
rownames(mat) <- genes
colnames(mat) <- colnames(Bulk)

Frac_scale = data.frame(Frac)[colnames(mat),]
### scale Frac into 0-1 by column
Frac_scale <- apply(Frac_scale,2,function(x) {(x-min(x)) / (max(x)-min(x))})

top_anno <- HeatmapAnnotation(df = Frac_scale,
                              col = list(Cycling = colorRamp2(c(0, 1), c("white", "#5cb85c")),
                                         Transition = colorRamp2(c(0, 1), c( "white", "#428bca")),
                                         Senescent= colorRamp2(c(0, 1), c("white", "#d9534f"))),
                              show_legend = F)
col <- colorRamp2(c(-2,0,2), c("blue","white", "red"), space = "LAB")
png('/home/wangjing/wangj/codebase/HUSI/Figures/TCGA_SKCM_heatmap.png',width = 1500,height = 1000,res= 200)
Heatmap(mat,
        show_column_names = F,
        show_row_names = F,
        row_title = NULL,
        col = col,
        cluster_rows = T,
        cluster_row_slices = FALSE,
        cluster_columns = T,
        column_km=3,
        row_km = 3,
        top_annotation = top_anno,
        heatmap_legend_param = list(title = "Expression",
                                         at = c(-2,-1,0,1,2),
                                         labels = c("-2", "-1", "0","1", "2"),
                                         labels_gp = gpar(fontsize = 12),
                                         title_gp = gpar(fontsize = 12, fontface = "bold"),
                                         legend_width = unit(30, "mm")))
dev.off()

### draw heatmap of SKCM CIBERSORT
png('/home/wangjing/wangj/codebase/HUSI/Figures/TCGA_SKCM_state_immue_cor.png',width = 2000,height = 1000,res= 250)
bk <- c(seq(-0.2,0.2,by=0.01))
pheatmap::pheatmap(cor_mat,show_colnames = T,show_rownames = T,cluster_rows = F,cluster_cols = T,fontsize = 12,border_color = "white",color = colorRampPalette(c("blue", "#f5f4f4", "red"))(length(bk)),legend_breaks=seq(-0.2,0.2,by=0.1),breaks=bk)
dev.off()

### survival plot
plotsurv <- function(myfit){
    p <- ggsurvplot(
    myfit,
    risk.table = F,
    pval = TRUE,
    conf.int = FALSE,
    xlim = c(0,4000),
    break.time.by = 1000,
    risk.table.y.text.col = T,
    risk.table.y.text = FALSE)
    return(p)
}

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'Cycling')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ Cycling,data = dat)
png('/home/wangjing/wangj/codebase/HUSI/Figures/TCGA_SKCM_survival_Cycling.png',width = 1200,height = 1200,res = 300)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'Transition')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ Transition,data = dat)
p2 <- plotsurv(fit)
png('/home/wangjing/wangj/codebase/HUSI/Figures/TCGA_SKCM_survival_Transition.png',width = 1200,height = 1200,res = 300)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'Senescent')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ Senescent,data = dat)
png('/home/wangjing/wangj/codebase/HUSI/Figures/TCGA_SKCM_survival_Senescent.png',width = 1200,height = 1200,res = 300)
plotsurv(fit)
dev.off()

### all melanome cell plot
mypalette <- read.csv("~/scripts/colors.csv",header = T)
bar = mypalette$palette1[1:length(unique(melanoma_obj$subtype))]
names(bar) <- unique(melanoma_obj$subtype)

png('/home/wangjing/wangj/codebase/HUSI/Figures/Melanoma_all_cell.png',width = 1500,height = 1200,res = 200)
DimPlot(melanoma_obj, group.by = 'subtype',reduction = "tsne",label = T,cols = bar,repel = TRUE,label.box = T,label.color = "white",pt.size = 1.5,label.size = 6)+
  theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank(),legend.position = "none")+ggtitle("Melanoma cell types")
dev.off()

### cellchat plot
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome_cellchat_scatter.png',width = 1000,height = 1200,res = 200)
source('/home/wangjing/wangj/codebase/HUSI/netAnalysis_signalingRole_scatter.R')
netAnalysis_signalingRole_scatter_log(cellchat,color.use = bar,dot.alpha = 0)
dev.off()

png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome_cellchat_network.png',width = 1600,height = 1000,res = 200)
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

source('/home/wangjing/wangj/codebase/HUSI/netVisual_bubble_my.R')
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome_cellchat_pathway.png',width = 1500,height = 1200,res = 300)
netVisual_bubble_my(cellchat,signaling = pathways,targets.use = c('Cycling','Transition','Senescent'), remove.isolate = F,sources.use = c('T cell','NK cell','Macro cell','CAF cell'),thresh=0.01)
dev.off()

png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome_cellchat_chord.png',width = 2000,height = 1000,res = 300)
par(mfrow = c(1,2), xpd=TRUE,mar = c(0.2, 0.2, 0.2,0.2))
for(p in c('BMP','TGFb')){
    netVisual_aggregate(cellchat, signaling = p,layout = "chord",remove.isolate = F,color.use = bar,sources.use = c('T cell','NK cell','Macro cell','CAF cell'),targets.use = c('Cycling','Transition','Senescent'))
}
dev.off()


png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome_cellchat_expression.png',width = 900,height = 1000,res = 300)
plotGeneExpression(cellchat, features=c('BMPR1B','BMPR2','TGFBR1','TGFBR2'),idents = c('Cycling','Transition','Senescent'),color.use = bar)
dev.off()

### receptor survival plot
cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'BMPR2')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ BMPR2,data = dat)
png('/home/wangjing/wangj/codebase/HUSI/Figures/TCGA_SKCM_survival_BMPR2.png',width = 1200,height = 1200,res = 300)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'TGFBR1')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ TGFBR1,data = dat)
png('/home/wangjing/wangj/codebase/HUSI/Figures/TCGA_SKCM_survival_TGFBR1.png',width = 1200,height = 1200,res = 300)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'TGFBR2')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ TGFBR2,data = dat)
png('/home/wangjing/wangj/codebase/HUSI/Figures/TCGA_SKCM_survival_TGFBR2.png',width = 1200,height = 1200,res = 300)
plotsurv(fit)
dev.off()