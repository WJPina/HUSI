library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(ggpubr)
library(RColorBrewer)
library(ggrepel)

mypalette <- read.csv("~/scripts/colors.csv",header = T)
cols = mypalette$palette7
names(cols) = levels(covid.m$SubCluster)

### Endo umap
png("/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_1010/Endo.png",width = 1800,height = 1500,res = 300)
DimPlot(covid.m, reduction = "umap",label = T,pt.size = 0.2,repel = F,label.box = T,group.by = 'SubCluster')+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 20),
        legend.position = "none")+
  ggtitle('Endothelial cell types')
dev.off()

### patient data 
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_1010/donor.png',width = 1800,height = 1500,res = 400)
covid.m@meta.data[,c('Donor','DTD')] %>%
  distinct() %>%
  mutate(counts = as.numeric(table(covid.m$Donor))) %>%
  mutate(Progress = ifelse(DTD<15,'severe','moderate')) %>%
  ggplot()+
  geom_point(aes(x = reorder(Donor,-DTD,),y = DTD,size = counts,color=Progress))+
  theme_bw()+
  theme(text = element_text(size = 16),axis.text.x = element_text(angle = -90,vjust = 0.5,hjust = 0))+
  scale_color_manual(values = c('#457b9d',"#bc4749"))+
  xlab('Days to Death')
dev.off()

### Progress group
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_1010/progress.png',width = 1800,height = 1500,res = 300)
ggplot(covid.m@meta.data)+
  geom_bar(aes(x=SubCluster, fill = Progress),position = 'fill')+
    theme_bw()+
    theme(legend.position = 'top',text = element_text(size = 16),axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1))+
    scale_fill_manual(values = c('#457b9d',"#bc4749"))+
  ylab('Percent')
dev.off()

### Progress hUSI
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_1010//HUSI_progress.png',width = 1600,height = 1500,res = 300)
ggviolin(covid.m@meta.data,
         x="Progress", y="hUSI", fill = "Progress", 
         palette = c('#457b9d',"#bc4749"),
         add = "boxplot", 
         add.params = list(fill="white"),
          )+
  theme(legend.position = 'top',axis.title = element_text(face = 'bold'),
        axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = list(c('severe','moderate')),
                     method = 't.test',
                     method.args = list(alternative = "greater"))+
  facet_wrap(~SubCluster)
dev.off()

### hUSI hist
png("/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_1010/hist.png",width = 1000,height = 1200,res = 300)
ggplot(subEndo@meta.data)+
  geom_histogram(aes(x=hUSI,y = ..density..),bins = 50,fill = 'white',color = '#adb5bd',linewidth = 0.4)+
  stat_function(fun = function(x) dnorm(x, mean = mean(subEndo$hUSI), 
                                           sd = sd(subEndo$hUSI)),color ='#d62828', linewidth = 1) +
  theme_classic()+
  theme(text = element_text(size = 14),legend.position = c(0.25,0.7),plot.margin = margin(10,10,10,10))+
  ylab('Count')
dev.off()

### cor gene volcano
res <- data.frame(expression = rowMeans(subEndo@assays$RNA@data[names(cor.genes),]),correlation = cor.genes,row.names = names(cor.genes))
res$threshold = ifelse(rownames(res) %in% features, ifelse(res$correlation > 0,"Positive","Negative"),"Background")
label_data <- subset(res, abs(res$correlation) > 0.1 & res$expression  > 1)
png("/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_1010/corgenes.png",width = 1000,height = 1200,res = 300)
ggplot(data = res, aes(y=expression,x=correlation,color=threshold)) +
  geom_point(alpha=0.5, size=2)+
  geom_text_repel(data = label_data,label = rownames(label_data),color = 'black',max.overlaps = 10,size = 3,show.legend = F)+
  scale_color_manual(values = c("Positive" = "#FF745A","Negative"="#007E99","Background" = "#ABABAB"))+
  ylab('Mean normalized expression')+
  xlab('Pearson cofficient ')+
  theme_classic()+
  theme(text = element_text(size = 16),legend.position = 'none')
dev.off()

### diffusion map
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_1010/DPT_diffusion.png',width = 1500,height = 1200,res = 300)
FeaturePlot(subEndo,features = 'DPT',reduction = 'diffusion',order = T,pt.size = 1.5) + 
  xlim(-0.04,0.05) + ylim(-0.07,0.05)+
  scale_color_gradient(low = 'white',high = 'red')+
  ggtitle('Artery EC')+
  theme(legend.position = 'bottom')
dev.off()

png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_1010/Marker_diffusion.png',width = 2000,height = 2400,res = 300)
DefaultAssay(subEndo) <- 'RNA'
(FeaturePlot(subEndo,features = 'IL6',reduction = 'diffusion',order = T,pt.size = 1.5) + 
  xlim(-0.04,0.05) + ylim(-0.07,0.05)+
  scale_color_gradient(low = 'white',high = 'red'))/
(FeaturePlot(subEndo,features = 'TPM4',reduction = 'diffusion',order = T,pt.size = 1.5) + 
  xlim(-0.04,0.05) + ylim(-0.07,0.05)+
  scale_color_gradient(low = 'white',high = 'red'))
dev.off()

### cor DPT hUSI
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_1010/hUSI_DPT.png',width = 1100,height = 1200,res = 300)
subEndo@meta.data %>%
  ggscatter(x = "hUSI", y = "DPT",
                color = '#8d99ae',
                add = "reg.line", 
                conf.int = TRUE,
                size = 1,
                add.params = list(color = "#1d3557"),
                 ggtheme = theme_classic())+ 
  stat_cor(method = "spearman",label.x = 0.01, label.y = 1,color='black')+
  theme_classic()+
  theme(text = element_text(size = 12),plot.margin = margin(10,10,10,10))
dev.off()

### ICAnet heatmap
DefaultAssay(subEndo) <- 'IcaNet'
subEndo <- ScaleData(subEndo)

genes = rownames(assoRes)
mat = subEndo@assays$IcaNet@scale.data[genes,names(sort(subEndo$DPT,decreasing = F))] %>% as.matrix()
top_anno <- HeatmapAnnotation(df = data.frame(DPT = subEndo$DPT[colnames(mat)]),show_legend = T,
                              col = list(DPT = colorRamp2(c(0,2,4), c("#f2cc8f", "#f4f1de",'#81b29a'))))
left_anno = rowAnnotation(df = data.frame(waldStat = assoRes[rownames(mat),'waldStat']),show_legend = T,
                          col = list(waldStat = colorRamp2(c(30,90,130), c("#f0f3bd","#00a896",'#05668d'))),
                          gp = gpar(fontsize = 12))
col <- colorRamp2(c(-1,0,1), c("#023e8a","white", "#e63946"), space = "LAB")
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_1010/ICAnet_heatmap.png',width = 2000,height = 1500,res = 300)
Heatmap(mat,
        show_column_names = F,
        show_row_names = T,
        col = col,
        top_annotation = top_anno,
        left_annotation = left_anno,
        row_title = NULL,
        cluster_rows = T,
        cluster_columns = F,
        heatmap_legend_param = list(title = "Scaled\nexpression"))
dev.off()












