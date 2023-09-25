library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(ggpubr)
library(RColorBrewer)

mypalette <- read.csv("~/scripts/colors.csv",header = T)
cols = mypalette$palette7
names(cols) = levels(Endo.m$celltype)

### Endo umap
png("/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_new/Endo.png",width = 1800,height = 1500,res = 300)
DimPlot(Endo.m, reduction = "umap",label = T,pt.size = 0.2,repel = F,label.box = T,group.by = 'celltype')+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 20),
        legend.position = "none",
        panel.border = element_rect(fill=NA,color="#afadad"))+
  ggtitle('Endothelial cell types')
dev.off()


### Progress group
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_new/HUSI_progress.png',width = 1500,height = 1500,res = 250)
ggviolin(Endo.m@meta.data,
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
  facet_wrap(~celltype)
dev.off()

### trajectory plot
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_new/CA_fate.png',width =800,height = 600,res = 300)
plot_cell_trajectory(cds,color_by="State", cell_size=0.5,show_backbone=TRUE) + 
  scale_color_manual(values = c("#2a9d8f","#e9c46a",'#264653'))+
  theme(legend.position = 'none')
# plot_cell_trajectory(cds,color_by="hUSI", cell_size=0.5,show_backbone=TRUE)+
#   scale_colour_gradientn(colors = colorRampPalette(c( 'blue',"white",'red'))(62))
dev.off()



### branch gene heatmap
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_new/CA_fate_genes.png',width = 1300,height = 1500,res = 300)
plot_genes_branched_heatmap(cds[genes,],
                            branch_point = 1,
                            num_clusters = 2,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = F,
                            hmcols = colorRampPalette(rev(brewer.pal(9, "Spectral")))(62),
                            branch_colors = c("#2a9d8f","#e9c46a",'#264653'),
                            branch_labels = c("Cell fate 1", "Cell fate 2"))
dev.off()

bar = c("#2a9d8f","#e9c46a",'#264653')
names(bar) = c('Root','Cell fate 1','Cell fate 2')

### middle group state percentage
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_new/CA_percent.png',width = 800,height = 1000,res = 300)
cds@phenoData@data %>% 
  mutate(Fate = ifelse(State == '1','Root',ifelse(State == '2','Cell fate 1','Cell fate 2'))) %>%
  ggplot()+
  geom_bar(aes(x = Progress,fill = Fate),stat = 'count',position = 'fill')+
  theme_classic()+
  theme(axis.title = element_text(face = 'bold'),
        text = element_text(size = 16),
        legend.position = 'none')+
  # theme(text = element_text(size = 16),
  #       axis.title = element_blank(),
  #       legend.position = c(0.15,1.25),legend.direction = 'horizontal',
  #       plot.margin=unit(c(2,1,1,1),'lines'))+
  scale_fill_manual(values = bar)
dev.off()

### hUSI fate
# cds$Fate = ifelse(cds$State == '1','Root',ifelse(cds$State == '2','Cell fate 1','Cell fate 2'))
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_new/HUSI_fate.png',width = 1000,height = 1200,res = 300)
ggviolin(obj@meta.data,
         x="State", y="hUSI", fill = "State", 
         palette =  c("#2a9d8f",'#264653',"#e9c46a"),
         add = "boxplot", 
         add.params = list(fill="white"),
)+
  theme(axis.title = element_text(face = 'bold'),
        legend.position = 'none')+
  stat_compare_means(comparisons = list(c('Senescent','Stem-like'),c('Stem-like','Normal'),c('Senescent','Normal')),
                     method = 't.test',
                     method.args = list(alternative = "greater"))
dev.off()


enrichplot::cnetplot(cluster_enrich,circular=F,
                     colorEdge = TRUE,showCategory = 20,
                     cex_gene = 0.5,cex_category = 0.8,
                     cex_label_gene = 0.6,node_label = 'gene') + 
  theme(legend.position = 'none')




### gene overlap
colorList = list()
colorList[['Stem-like_up']] <- c("#FF745A","#e9c46a")
colorList[['Stem-like_down']] <- c("#007E99","#e9c46a")
colorList[['Senescent_up']] <- c("#FF745A",'#264653')
colorList[['Senescent_down']] <- c("#007E99",'#264653')

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
           text_size = 4,
           fill_alpha = 0.8)+
    labs(title = paste("Padj value:",as.character(signif(enrich_reList[[set]],2))))+
    theme(plot.title = element_text(hjust = 0.5,vjust = 0,size = 15))+
    scale_y_continuous(limits = c(-1, 1))
  
}
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_new/CA_gene_enrich.png',width = 1200,height = 800,res = 200)
ggarrange(pList[['Stem-like_up']],
          pList$Senescent_up,
          pList[['Stem-like_down']],
          pList$Senescent_down,
          ncol = 2,nrow = 2,legend = "none")
dev.off()


### ROI bulk ssGSEA
library(rstatix)
names(bar) = c('Normal','Stem-like','Senescent')

cbind(meta,t(results[,rownames(meta)]))%>%
  dplyr::select(c('Normal','Senescent','Stem-like','Primary_Morph','Progress')) %>%
  reshape2::melt(variable.name = 'State',value.name = 'ssGSEA') %>%
  ggboxplot(x = 'Progress',y = 'ssGSEA',fill = 'State')+
  scale_fill_manual(values = bar)+
  theme_classic()


### stem-like or senescence marker
### stem
VlnPlot(subEndo,features = c('CD9','CD44','CDKN1A','CDKN2A'),group.by = 'age_state',pt.size = 1.5,col = bar)





