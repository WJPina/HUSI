### This script is for Fig5 and FigS6 
################################ library #######################################
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
library(ComplexHeatmap)
library(circlize)
library(cowplot)

setwd("hUSI/")
############################### all cell types #################################
bar = c("#eb0101","#AAC5E2","#FFC592","#A4D99D","#F4A39E","#CCB7DC","#CBAEA8")
names(bar) <- c("Tumor","B","CAF","Endo","Macro","NK","T")

### cell types
pdf('Results/Melanoma/Melanoma_all_cell.pdf',width = 4.4,height = 4)
DimPlot(melanoma_obj, group.by = 'celltype',reduction = "tsne",label = T,
        cols = bar,repel = TRUE,label.box = T,label.color = "white",pt.size = 1.5,label.size = 6)+
  ggtitle("Melanoma cell types")+NoLegend()
dev.off()

### hUSI
pdf('Results/Melanoma/Melanoma_all_cell_hUSI.pdf',width = 4.8,height = 4)
FeaturePlot(melanoma_obj, 'hUSI',reduction = "tsne",order = T,
            repel = TRUE,label.color = "white",pt.size = 1,label.size = 6)+ 
  scale_color_gradient2(low ='#3AB370' ,mid = "#EAE7CC",high = "#FD1593",midpoint = 0.5)
dev.off()
################################# tumor hUSI ###################################
### highlight tumor cells
bar = c("#eb0101",rep("grey80",6))
names(bar) <- c("Tumor","B","CAF","Endo","Macro","NK","T")

pdf('Results/Melanoma/Melanoma_tumor.pdf',width = 4.4,height = 4)
DimPlot(melanoma_obj, group.by = 'celltype',reduction = "tsne",order=T,
        cols = bar,repel = TRUE,pt.size = 1.5,label.size = 6)+
  ggtitle("Melanoma tumor cells")+NoLegend()
dev.off()


### senescence Class
bar = c("#2a9d8f","#e9c46a","#e76f51")
names(bar) = c("C1","C2","C3")

### histogram
pdf("Results/Melanoma/Melanoma_tumor_hist.pdf",width = 4.4,height = 4)
ggplot(EpiExp.m@meta.data)+
  # geom_histogram(aes(x=hUSI),bins = 50,fill = 'white',color = '#adb5bd',linewidth = 0.4)+
  geom_density(aes(x=hUSI,fill = SenClass,color = SenClass),alpha = 0.6)+
  scale_color_manual(values = bar)+
  scale_fill_manual(values = bar)+
  theme_classic()+
  theme(text = element_text(size = 14),legend.position = c(0.25,0.7))+
  ylab('Count')
dev.off()

### SenClass hUSI level
pwc = EpiExp.m@meta.data %>%
  rstatix::pairwise_t_test(hUSI ~ SenClass, paired = F,p.adjust.method = "bonferroni") %>%
  rstatix::add_xy_position(x = "SenClass")
pwc

pdf('Results/Melanoma/Melanoma_tumor_hUSI.pdf',width = 4,height = 4)
EpiExp.m@meta.data %>%
  ggplot(aes(x = SenClass,y = hUSI,color = SenClass))+
  geom_boxplot()+
  scale_color_manual(values = bar)+
  labs(caption = get_pwc_label(pwc))+
  stat_pvalue_manual(pwc)+
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.position = 'none')
dev.off()

########################## tumor Senescence state ##############################
### icanet heatmap
library(ComplexHeatmap)
library(circlize)

bar = c("#2a9d8f","#e9c46a","#e76f51")
names(bar) = c('Cycling', 'Transitional', 'Senescent')

groups = sort(EpiExp.m$hUSI,decreasing = F)
top_anno <- HeatmapAnnotation(df = data.frame(hUSI = groups,State = EpiExp.m$State[names(groups)]),show_legend = T,
                              col = list(hUSI = colorRamp2(c(0,0.5,1), c('#81b29a',"#f4f1de","#f2cc8f")),State = bar))
col <- colorRamp2(c(-1,0,1), c("#023e8a","white", "#e63946"), space = "LAB")
DefaultAssay(EpiExp.m) <- 'IcaNet'
EpiExp.m <- ScaleData(EpiExp.m)
mat = EpiExp.m@assays$IcaNet@scale.data %>% as.matrix()
mat = mat[,names(groups)]

pdf('Results/Melanoma/Melanoma_icanet.pdf',width = 7,height =6)
Heatmap(mat,
        show_column_names = F,
        show_row_names = F,
        row_title = NULL,
        col = col,
        cluster_rows = T,
        cluster_row_slices = FALSE,
        cluster_columns = F,
        top_annotation = top_anno,
        row_names_gp = gpar(fontsize = 12),
        heatmap_legend_param = list(title = "Scaled\nexpression"))
dev.off()

### trajectory
pdf('Results/Melanoma/Melanoma_phate.pdf',width =5,height =4)
DimPlot(EpiExp.m, reduction = 'phate', group.by = 'State',label=T,label.box=T,cols=bar)+
  theme_classic()+
  theme(legend.position = 'none',text = element_text(size =14),plot.margin = margin(1,6,1,5,unit = 'mm'))+
  ggtitle('Melanoma tumor cells')
dev.off()

### aging markers
DefaultAssay(EpiExp.m)='RNA'
SenMarkers = c("CDKN1A", "SERPINE1")

p1 = FeaturePlot(EpiExp.m,features=SenMarkers[1],reduction='phate',order=F,pt.size = 2) + xlim(-0.04,0.06) + ylim(-0.04,0.06)
p2 = FeaturePlot(EpiExp.m,features=SenMarkers[2],reduction='phate',order=F,pt.size = 2) + xlim(-0.04,0.06) + ylim(-0.04,0.06)
p3 = VlnPlot(EpiExp.m,features=SenMarkers[1],group.by='State',assay = 'RNA',cols = bar,pt.size = 0) + theme(axis.title.x = element_blank()) + NoLegend()
p4 = VlnPlot(EpiExp.m,features=SenMarkers[2],group.by='State',assay = 'RNA',cols = bar,pt.size = 0) + theme(axis.title.x = element_blank())+ NoLegend()

pdf('Results/Melanoma/Melanoma_markers.pdf',width = 8.4,height = 7)
ggarrange(p1,p3,p2,p4,ncol = 2,nrow = 2) 
dev.off()

################################## immune therapy ########################
pdf('Results/Melanoma/Melanoma_tumor_hUSI_treatment.pdf',width = 3,height = 4)
EpiExp.m$treatment = factor(EpiExp.m$treatment,levels = c('untreated','treated'))
VlnPlot(EpiExp.m,features = 'hUSI',group.by = 'treatment',cols = c('untreated'='#999999','treated'='#009999'),pt.size = 0)+
  geom_boxplot(width=0.2,outliers = F)+
  theme(legend.position = 'none',axis.text.x = element_text(angle = 0,hjust = 0.5),plot.margin = margin(0,0,0,0))+
  stat_compare_means(aes(label = ..p.signif..),label.x = 1.5, label.y = 1.03)+
  ylim(c(0,1.05))+
  xlab('Immune trerapy')
dev.off()

pdf('Results/Melanoma/Melanoma_tumor_treatment_percent.pdf',width = 3.4,height = 4)
table(EpiExp.m$treatment,EpiExp.m$State) %>% as.array() %>% data.frame() %>% 
  set_names(c('Treatment','SenClass','Percent')) %>%
  ggplot(aes(x=Treatment,y=Percent,fill=SenClass))+
  geom_bar(stat = 'identity',position = 'fill')+
  scale_fill_manual(values = bar)+
  theme_bw()+
  xlab('Immune trerapy')
dev.off()

##################### Characterizing Senescence state ##########################
### validate in microarray markers
colorList = list()
colorList[['Cycling_up']] <- c("#993D00","#2a9d8f")
colorList[['Cycling_down']] <- c("#007E99","#2a9d8f")
colorList[['Transitional_up']] <- c("#993D00","#e9c46a")
colorList[['Transitional_down']] <- c("#007E99","#e9c46a")
colorList[['Senescent_up']] <- c("#993D00","#e76f51")
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
           text_size = 3,fill_alpha = 0.8)+
    labs(title = paste("Padj value:",as.character(signif(enrich_reList[[set]],2))))+
    theme(plot.title = element_text(hjust = 0.5,vjust = 0,size = 12))+
    scale_y_continuous(limits = c(-1, 1))
  
}

pdf('Results/Melanoma/Melanoma_maker_overlap.pdf',width = 8,height =4)
ggarrange(pList$Senescent_up,
          pList$Transitional_up,
          pList$Cycling_up,
          pList$Senescent_down,
          pList$Transitional_down,
          pList$Cycling_down,
          ncol = 3,nrow = 2,legend = "none")
dev.off()

### enrichment for state markers
pdf('Results/Melanoma/Melanoma_state_enrich.pdf',width = 7,height = 5)
dotplot(go,showCategory = 5) +
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5))+
  ggtitle("GO:BP Enrichment of State Marker Genes")+
  theme(plot.title = element_text(hjust = 1))+
  scale_size(range = c(3,5))+
  scale_color_gradient2(high = '#FF7856',mid = '#F72735',low = '#A50021')
dev.off()

######################### deconvolution in TCGA-SCKM ###########################
pdf('Results/Melanoma/TCGA_SKCM_state_immue_cor.pdf',width =5.8,height = 4)
bk <- c(seq(-0.3,0.3,by=0.01))
# cor_mat <- cormat_gtex
pheatmap::pheatmap(cor_mat,show_colnames = T,show_rownames = T,
                   cluster_rows = F,cluster_cols = T,fontsize = 12,
                   border_color = "white",color = colorRampPalette(c("blue", "#f5f4f4", "red"))(length(bk)),
                   legend_breaks=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3),breaks=bk)
dev.off()

### survival plot
plotsurv <- function(myfit){
  p <- ggsurvplot(
    myfit,
    risk.table = F,
    pval = TRUE,
    conf.int = T,
    xlim = c(0,4000),
    break.time.by = 1000,
    risk.table.y.text.col = T,
    risk.table.y.text = FALSE)
  return(p)
}

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'Cycling')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ Cycling,data = dat)
pdf('Results/Melanoma/TCGA_SKCM_survival_Cycling.pdf',width = 4.2,height = 4)
plotsurv(fit)+xlab('Time (days)')
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'Transitional')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ Transitional,data = dat)
pdf('Results/Melanoma/TCGA_SKCM_survival_Transitional.pdf',width = 4.2,height = 4)
plotsurv(fit)+xlab('Time (days)')
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'Senescent')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ Senescent,data = dat)
pdf('Results/Melanoma/TCGA_SKCM_survival_Senescent2.pdf',width = 4.2,height = 4)
plotsurv(fit)+xlab('Time (days)')
dev.off()

######################### cellchat analysis  ###########################
bar = c("#2a9d8f","#e9c46a","#e76f51","#AAC5E2","#FFC592","#A4D99D","#F4A39E","#CCB7DC","#CBAEA8")
names(bar) <- c("Cycling","Transitional","Senescent","B","CAF","Endo","Macro","NK","T")

### cell chat strength count
pdf('Results/Melanoma/melanome_cellchat_scatter.pdf',width = 5,height =3.8)
netAnalysis_signalingRole_scatter_log(cellchat,color.use = bar,dot.alpha = 0)
dev.off()

### cell chat strength network
groupSize <- table(cellchat@idents) %>% as.numeric()
pdf('Results/Melanoma/melanome_cellchat_network.pdf',width = 6,height = 3)
par(mfrow = c(1,2), xpd=TRUE,mar = c(0, 1, 1, 0.5))

p1=netVisual_circle(cellchat@net$weight, 
                 targets.use = c("Cycling","Transitional","Senescent"),
                 color.use = bar[colnames(cellchat@net$count)],
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge = F, 
                 title.name = "Incoming weights/strength")
p2=netVisual_circle(cellchat@net$weight, 
                 sources.use = c("Cycling","Transitional","Senescent"),
                 color.use = bar[colnames(cellchat@net$count)],
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge = F, 
                 title.name = "Outgoing weights/strength")
p1+p2
dev.off()


### pathway strength bubble
pathways=pathways[!(pathways %in% c("CSPG4","CD46"))]
pdf('Results/Melanoma/melanome_cellchat_Incoming2.pdf',width = 5.1,height = 3)
netVisual_bubble_my(cellchat,signaling = pathways,
                    targets.use = c('Cycling','Transitional','Senescent'), 
                    remove.isolate = F,
                    sources.use = c('T','NK','Macro','CAF','Endo','B'),thresh=0.01)
dev.off()

### pathway strength cord
pdf('Results/Melanoma/melanome_cellchat_chord_TGFb.pdf',width = 3,height = 3)
par(mfrow = c(1,1), xpd=TRUE,mar = c(0, 0, 0.2,0))
netVisual_aggregate(cellchat, signaling = 'TGFb',layout = "chord",remove.isolate = F,
                    color.use = bar,sources.use = c('T','NK','Macro','CAF','Endo','B'),
                    targets.use = c('Cycling','Transitional','Senescent'))
dev.off()

pdf('Results/Melanoma/melanome_cellchat_chord_CD6.pdf',width = 3,height = 3)
par(mfrow = c(1,1), xpd=TRUE,mar = c(0, 0, 0.2,0))
netVisual_aggregate(cellchat, signaling = 'CD6',layout = "chord",remove.isolate = F,
                    color.use = bar,sources.use = c('T','NK','Macro','CAF','Endo','B'),
                    targets.use = c('Cycling','Transitional','Senescent'))
dev.off()

pdf('Results/Melanoma/melanome_cellchat_chord_CCL.pdf',width = 3,height = 3)
par(mfrow = c(1,1), xpd=TRUE,mar = c(0, 0, 0.2,0))
netVisual_aggregate(cellchat, signaling = 'CCL',layout = "chord",remove.isolate = F,
                    color.use = bar,sources.use = c('T','NK','Macro','CAF','Endo','B'),
                    targets.use = c('Cycling','Transitional','Senescent'))
dev.off()

pdf('Results/Melanoma/melanome_cellchat_chord_CSPG4.pdf',width = 3,height = 3)
par(mfrow = c(1,1), xpd=TRUE,mar = c(0, 0, 0.2,0))
netVisual_aggregate(cellchat, signaling = 'CSPG4',layout = "chord",remove.isolate = F,
                    color.use = bar,sources.use = c('T','NK','Macro','CAF','Endo','B'),
                    targets.use = c('Cycling','Transitional','Senescent'))
dev.off()

pdf('Results/Melanoma/melanome_cellchat_chord_CD46.pdf',width = 3,height = 3)
par(mfrow = c(1,1), xpd=TRUE,mar = c(0, 0, 0.2,0))
netVisual_aggregate(cellchat, signaling = 'CD46',layout = "chord",remove.isolate = F,
                    color.use = bar,sources.use = c('T','NK','Macro','CAF','Endo','B'),
                    targets.use = c('Cycling','Transitional','Senescent'))
dev.off()


### pathway gene expression
pdf('Results/Melanoma/melanome_cellchat_expression_fourgenes.pdf',width = 2.5,height = 3)
plotGeneExpression(cellchat, features=c('TGFBR2','ACVR1','ALCAM','CCR10'),idents = c('Cycling','Transitional','Senescent'),color.use = bar)
dev.off()

### receptor survival plot
plotsurv <- function(myfit){
  p <- ggsurvplot(
    myfit,
    risk.table = F,
    pval = TRUE,
    conf.int = T,
    xlim = c(0,4000),
    break.time.by = 1000,
    risk.table.y.text.col = T,
    risk.table.y.text = FALSE)+xlab('Time (days)')
  return(p)
}

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'TGFBR2')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ TGFBR2,data = dat)
pdf('Results/Melanoma/TCGA_SKCM_survival_TGFBR2.pdf',width = 4.2,height = 4)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'ACVR1')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ ACVR1,data = dat)
pdf('Results/Melanoma/TCGA_SKCM_survival_ACVR1.pdf',width = 4.2,height = 4)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'ALCAM')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ ALCAM,data = dat)
pdf('Results/Melanoma/TCGA_SKCM_survival_ALCAM.pdf',width = 4.2,height = 4)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'CCR10')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ CCR10,data = dat)
pdf('Results/Melanoma/TCGA_SKCM_survival_CCR10.pdf',width = 4.2,height = 4)
plotsurv(fit)
dev.off()

##################################################
genes=list(Marker.geneset=c('TGFBR2','ACVR1','ALCAM','CCR10'))
exp = tcga_melanoma[,rownames(survival_data)] %>% as.matrix()
csore=GSVA::gsva(exp,genes,method="zscore",kcdf="Gaussian")
data = data.frame(t(tcga_melanoma[c('TGFBR2','ACVR1','ALCAM','CCR10'),rownames(survival_data)]),
                  OS = survival_data$vital_status,OS.time = survival_data$os.time,
                  Score=csore[1,] %>% as.numeric())

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'Score')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ Score,data = dat)
pdf('Results/Melanoma/TCGA_SKCM_survival_zscore_sore.pdf',width = 4.2,height = 4)
plotsurv(fit)
dev.off()


### different pattern
pdf('Results/Melanoma/melanome_cellchat_pattern_out.pdf',width = 6,height = 4)
netAnalysis_river(cellchat, pattern = "outgoing",cutoff = 0.5,
                  color.use = bar[c('Cycling','Transitional','Senescent')],
                  color.use.pattern = c('Pattern 1'='#e63946','Pattern 6'='#457b9d'),
                  sources.use = c('Cycling','Transitional','Senescent'),
                  signaling = as.character(pattern_use$Signaling))
dev.off()













































