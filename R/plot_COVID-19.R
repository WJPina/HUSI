### This script is for Fig4 and FigS5

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
library(colorspace)

mypalette = read.csv('Data/colors.csv',header = T)
setwd("hUSI/")

### cell types
pdf('Results/Covid-19/covid_all_cell.pdf',width = 5.5,height = 5)
DimPlot(covid, group.by = 'cell_type_main',reduction = "umap",label = T,raster=FALSE,
        cols = mypalette$palette3,repel = TRUE,label.box = T,label.color = "white",pt.size = 0.1,label.size = 6)+
  ggtitle("COVID-19 cell types")+NoLegend()
dev.off()

### hUSI
pdf('Results/Covid-19/covid_all_cell_hUSI.pdf',width = 6,height = 5)
FeaturePlot(covid, 'hUSI',reduction = "umap",order = T,raster=FALSE,
            repel = TRUE,label.color = "white",pt.size = 0.1)+ 
  scale_color_gradient2(low ='#3AB370' ,mid = "#EAE7CC",high = "#FD1593",midpoint = 0.5)
dev.off()

### hUSI levels
df_plot = covid@meta.data[,c('cell_type_fine','hUSI','group','age','cell_type_main')]
pwc <- df_plot %>% group_by(cell_type_fine) %>% 
  rstatix::pairwise_t_test(hUSI ~ group, paired = F,p.adjust.method = "bonferroni")
pwc
pwc <- pwc %>% rstatix::add_xy_position(x = "cell_type_fine")


pdf('Results/Covid-19/covid_all_cell_hUSI_levels.pdf',width = 11,height = 5)
df_plot %>%
  ggplot(aes(x = cell_type_fine, y = hUSI)) + 
  geom_violin(aes(color = group),outlier.shape = NA) +
  geom_boxplot(aes(color = group),width=0.3,position = position_dodge(0.9), outlier.shape = NA, alpha = 0.6) +
  stat_pvalue_manual(pwc)+
  theme_classic()+
  theme(text = element_text(size = 12),axis.text.x = element_text(size=12,angle = 45,hjust = 1),
        axis.title.x = element_blank(),plot.margin = margin(0,0,0,22))+
  scale_color_manual(values = c('Control' = '#959596','COVID-19'='#e63946'))+
  ylab('hUSI')
dev.off()

### UMAP in SenClass
pdf('Results/Covid-19/covid_all_SenClass_fine.pdf',width = 6.5,height = 4.5)
DimPlot(covid, group.by = 'SenClass_fine',reduction = "umap",label = F,raster=FALSE,
        cols = c('#91D1C2','#E64B35'),pt.size = 0.1)+
  ggtitle("COVID-19 cells senecence class")
dev.off()

### SASP expression heatmap
SASP =c("CDKN1A","IL6","CCL2","CXCL10","SERPINE1","THBS1","TIMP1",
        'UBC','RPS27A','JUN','UBB','STAT3','UBE2E1',"RPS6KA2","RPS6KA3",
        "DHFR","MKI67","LMNB1")

covid = ScaleData(covid,features = SASP,split.by = 'cell_type_fine')

Mean_sasp = AverageExpression(covid,assays = 'RNA',slot = 'scale.data',features = SASP,
                              group.by = c('cell_type_fine','SenClass_fine'))$RNA

Mean_sasp[1:5,]

myorder = rowMeans(Mean_sasp[,grepl('Senescent',colnames(Mean_sasp))]) %>% sort(decreasing = T) %>% names()
mat = Mean_sasp[myorder,]

mat = mat[rowSums(mat)!=0,]
celltype = gsub('_.*','',colnames(mat))

colnames(mat) = gsub('.*_','',colnames(mat))
left_anno = rowAnnotation(df = data.frame(genes= factor(rep(c("up","donw"),times=c(15,3)),levels = c("up","donw"),ordered = T)),show_legend = T,
                          col = list(genes = c('up' = '#B26F2C','donw' = '#6B990F')),show_annotation_name = F)
col <- colorRamp2(c(-0.05,0,0.05), c("#023e8a","#FFEBEE","#e63946"), space = "LAB")


pdf('Results/Covid-19/SASP_SenClass_fine2.pdf',width = 16,height = 5.5)
Heatmap(mat,
        show_column_names = F,
        show_row_names = T,
        row_title = NULL,
        col = col,
        cluster_rows = F,
        cluster_row_slices = FALSE,
        cluster_columns = F,
        # left_annotation = left_anno,
        row_names_gp = gpar(fontsize = 12),
        border = 1,
        border_gp = gpar(col = "black"),
        rect_gp=gpar(col="white",lwd=1),
        heatmap_legend_param = list(title = "scaled\nexpression"),
        column_split = factor(celltype,levels = levels(covid$cell_type_fine)),
        gap = unit(2, "mm"),
        column_gap = unit(3, "mm"),
        column_title_rot = 45,
        row_split = factor(rep(c("up","donw"),times=c(15,3)),levels = c("up","donw"),ordered = T))
dev.off()


#### enrichment
df_plot = do.call(rbind,pathwaysList)
df_plot = filter(df_plot,p.adjust<0.05&Count>=3)
# df_plot = df_plot %>% group_by(celltype) %>% top_n(.,n = 10,wt = -p.adjust)
df_plot$celltype = factor(df_plot$celltype,ordered = T)
df_plot = df_plot[order(df_plot$celltype),]
df_plot$pathway = gsub('KEGG_','',df_plot$Description) %>% gsub('GOBP_','',.)

pathway_show = table(df_plot$pathway) %>% sort(decreasing = T)
pathway_show = names(pathway_show[pathway_show>=9])
pathway_show

pathway_show = c("Viral myocarditis","Cellular senescence",'Antigen processing and presentation',
                 "ECM-receptor interaction","IL-17 signaling pathway","p53 signaling pathway")
df_plot = filter(df_plot,pathway %in% pathway_show )
df_plot$pathway = factor(df_plot$pathway,levels = pathway_show,ordered = T)

head(df_plot)

pdf('Results/Covid-19/covid_SenClass_enrich_KEGG.pdf',width = 5.5,height = 5.6)
ggplot(df_plot,aes(y=celltype,x=pathway,color=p.adjust,size=Count))+
  geom_point()+
  scale_color_gradient2(high = '#FF7856',mid = '#F72735',low = '#A50021')+
  scale_size(range = c(5,8))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("KEGG Enrichment for senescent cells")+
  ylab('Cell type fine')+
  xlab('Pathway')
dev.off()


df=covid@meta.data
res=df %>% group_by(cell_type_fine,SenClass_fine) %>% dplyr::summarise(num=n()) %>% 
  merge(.,df %>% group_by(cell_type_fine) %>% dplyr::summarise(all=n()),by='cell_type_fine') %>% mutate(por=num/all)

cf=res %>% filter(SenClass_fine=="Senescent cells") %>% arrange(por)
pdf("Results/Covid-19/Class_percentage_col_plot.pdf",width = 12,height = 6)
ggplot(res,aes(x=factor(cell_type_fine,levels = cf$cell_type_fine),y=por,fill=SenClass_fine))+
  geom_bar(stat = "identity",position = "fill")+theme_bw()+
  scale_fill_manual(values = c('Senescent cells' = '#E64B35' ,'Normal cells'='#91D1C2'))+
  labs(x="",y="Percentages",fill="Celltype")+
  theme(axis.text = element_text(size=13),
        axis.text.x = element_text(angle=45,hjust=1),
        axis.title = element_text(size=15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size=15))+
  geom_hline(yintercept = .5,linetype="dashed",color="grey60",size=1)
dev.off()

pdf("Results/Covid-19/Class_percentage_bar_plot.pdf",width = 14,height = 6)
ggplot(res,aes(x=factor(cell_type_fine,levels = cf$cell_type_fine),y=por,fill=SenClass_fine))+
  geom_bar(stat = "identity",position=position_dodge())+theme_bw()+
  scale_fill_manual(values = c('Senescent cells' = '#E64B35' ,'Normal cells'='#91D1C2'))+
  scale_color_manual(values = c('white'))+
  labs(x="",y="Percentages",fill="Celltype")+
  theme(axis.text = element_text(size=14),
        axis.text.x = element_text(angle=45,hjust=1),
        axis.title = element_text(size=16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size=16))
  # geom_hline(yintercept = .5,linetype="dashed",color="grey60",size=1)
dev.off()

#### frac diff COVID-19 - Control
df_diff$Changes = factor(ifelse(df_diff$Frac_diff>0,'Increased','Decreased'),levels = c('Increased','Decreased'))

pdf('Results/Covid-19/covid_SenClass_fra_fine_top.pdf',width = 10.5,height = 5.5)
filter(df_diff,group=='COVID-19' & SenClass_fine %in% c('Senescent cells')) %>%
ggplot(aes(x=reorder(cell_type_fine,Frac_diff),y=Frac_diff,color=Changes))+
  geom_point(aes(size=SenFrac))+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 12,angle = 45,hjust = 1),
        axis.title  = element_text(size = 14))+
  xlab('Cell type fine')+
  ylab('Fraction')+
  ggtitle('Fraction difference of senescent cells (COVID - Control)')+
  scale_color_manual(values = c('Decreased' = '#959596','Increased'='#e63946'))+
  geom_hline(yintercept = 0,linetype='dashed', color='grey30')+
  theme(strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),plot.margin = margin(0,0,0,40))+
  scale_size(name="Fraction in\nCOVID-19") 
dev.off()

#### frac diff Severe - Moderate
df_diff$Changes = factor(ifelse(df_diff$Frac_diff>0,'Increased','Decreased'),levels = c('Increased','Decreased'))
pdf('Results/Covid-19/covid_SenClass_fra_fine_DTD.pdf',width = 10.5,height = 5.5)
filter(df_diff,Prognosis=='Severe' & SenClass_fine %in% c('Senescent cells')) %>%
  ggplot(aes(x=reorder(cell_type_fine,Frac_diff),y=Frac_diff,color=Changes))+
  geom_point(aes(size=SenFrac))+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12,angle = 45,hjust = 1),
        axis.title  = element_text(size = 14))+
  xlab('Cell type fine')+
  ylab('Fraction')+
  ggtitle('Fraction difference of senescent cells (Severe - Moderate)')+
  scale_color_manual(values = c('Decreased' = '#959596','Increased'='#e63946'))+
  geom_hline(yintercept = 0,linetype='dashed', color='grey30')+
  theme(strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),plot.margin = margin(0,0,0,40))+
  scale_size(name="Fraction in\nSevere") 
dev.off()

