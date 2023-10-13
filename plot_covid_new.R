library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(ggpubr)
library(RColorBrewer)

mypalette <- read.csv("~/scripts/colors.csv",header = T)
cols = mypalette$palette7
names(cols) = levels(covid.m$celltype)

### Endo umap
png("/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_1010//Endo.png",width = 1800,height = 1500,res = 300)
DimPlot(covid.m, reduction = "umap",label = T,pt.size = 0.2,repel = F,label.box = T,group.by = 'celltype')+
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
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_1010//HUSI_progress.png',width = 1500,height = 1500,res = 300)
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
  facet_wrap(~Cluster)
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
                            branch_labels = c('Anti-senescence','Senescence'))
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
ggviolin(subEndo@meta.data,
         x="State", y="hUSI", fill = "State", 
         palette =  c("#2a9d8f",'#264653',"#e9c46a"),
         add = "boxplot", 
         add.params = list(fill="white"),
)+
  theme(axis.title = element_text(face = 'bold'),
        legend.position = 'none')+
  stat_compare_means(comparisons = list(c('Normal','Root'),c('Senescent','Noraml'),c('Senescent','Root')),
                     method = 't.test',
                     method.args = list(alternative = "greater"))
dev.off()


### pathway network
library(enrichplot)
cols = mypalette$palette3[1:length(ids)]
# id_1 = df_plot$ID
# id_2 = df_plot$ID
# ids = unique(c(id_1,id_2))
names(cols) = ids
cols

foldChange = BEAM_res[genes,]$qval
names(foldChange) = genes
p = cnetplot(cluster_enrich,
         circular=F,
         colorEdge = F,
         showCategory = 20,
         cex_gene = 0.8,
         cex_category = 0.8,
         foldChange = -log10(foldChange),
         cex_label_gene = 0.6,
         node_label = 'gene') + 
  theme(legend.position = 'none')

p[["layers"]][[2]][["mapping"]][["colour_new"]][[2]][[2]] <- cols[cluster_enrich@result$ID]
p[["data"]][["color"]][which(!(p[["data"]][["name"]] %in% names(foldChange)) == T)] <- cols[cluster_enrich@result$ID]

png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_new/Ant-senescence.png',width = 1000,height = 1000,res = 300)
print(p)
dev.off()

### plot legend
# labels_1 = df_plot$label
# labels_2 = df_plot$label
# labels = unique(c(labels_1,labels_2))
# names(cols) = labels

comm = c("focal adhesion","ecm-receptor interaction",
         "human papillomavirus infection","vascular smooth muscle contraction")
sene = c("regulation of actin cytoskeleton","glycosaminoglycan binding proteins",
         "arrhythmogenic right ventricular cardiomyopathy","hypertrophic cardiomyopathy",
         "dilated cardiomyopathy","platelet activation","apelin signaling pathway")
anti = c("pi3k-akt signaling pathway","cd molecules",
         "rap1 signaling pathway","cell adhesion molecules",
         "amoebiasis","tgf-beta signaling pathway","hippo signaling pathway",
         "axon guidance","african trypanosomiasis" ,"adherens junction",
         "small cell lung cancer", "malaria")

png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19_new/legend.png',width = 2000,height = 2000,res = 300)
bar = c(cols[anti],cols[sene],cols[comm])
palette(bar)
plot(y = 1:length(bar),x=rep(0.2,length(bar)), col = 1:length(bar), pch = 19,xlim = c(0,2),cex = 2) 
labels <- names(bar)
text(y = 1:length(bar),x=rep(0.2,length(bar)), labels = labels, pos = 4)
dev.off()

palette("default")




### gene overlap
colorList = list()
colorList[['Anti_senescence_up']] <- c("#FF745A","#e9c46a")
colorList[['Anti_senescence_down']] <- c("#007E99","#e9c46a")
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
ggarrange(pList[['Anti_senescence_up']],
          pList$Senescent_up,
          pList[['Anti_senescence_down']],
          pList$Senescent_down,
          ncol = 2,nrow = 2,legend = "none")
dev.off()


### ROI bulk ssGSEA
library(rstatix)
names(bar) = c('Normal','Anti_senescence','Senescent')

cbind(meta,t(results[,rownames(meta)]))%>%
  dplyr::select(c('Normal','Senescent','Anti_senescence','Primary_Morph','Progress')) %>%
  reshape2::melt(variable.name = 'State',value.name = 'ssGSEA') %>%
  ggboxplot(x = 'Progress',y = 'ssGSEA',fill = 'State')+
  scale_fill_manual(values = bar)+
  theme_classic()


### Anti_senescence or senescence marker
### stem
VlnPlot(subEndo,features = c('CD9','CD44','CDKN1A','CDKN2A'),group.by = 'age_state',pt.size = 1.5,col = bar)





