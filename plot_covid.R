library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(ggpubr)

### age state
bar = c("#5cb85c","#428bca","#d9534f")
names(bar) = c("Normal","Transition","Senescent")
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19/Endo_pca.png',width = 1800,height = 800,res = 200)
DimPlot(Endo.m, reduction = 'pca', group.by = 'age_state',label=F,pt.size=2,cols = bar)+ggtitle('COVID-19 Lung Endothelial cell')
dev.off()
### sasp map
pList = list()
for (gene in sasp_plot) {
   pList[[gene]] <- FeaturePlot(Endo.m,features = gene,reduction = 'pca',pt.size=2) 
}
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19/Endo_marker.png',width = 2500,height = 1500,res = 200)
ggarrange(pList[[1]],pList[[2]],pList[[3]],pList[[4]],ncol = 2,nrow = 2,legend = "none") 
dev.off()

pList = list()
for (gene in sasp_plot) {
    pList[[gene]] = VlnPlot(Endo.m,features=gene,group.by='age_state',assay = 'RNA',cols = bar,pt.size = 0) + theme(axis.title.x = element_blank())
}

png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19/Endo_markers_vln.png',width = 1800,height = 1200,res = 200)
ggarrange(pList[[1]],pList[[2]],pList[[3]],pList[[4]],ncol = 2,nrow = 2,legend = "none") 
dev.off()


### marker enrichment plot
EnrichBarPlot <- function(data) {
    ggplot(data) +
    geom_bar(aes(y = reorder(ID, -log10(p.adjust)), x = -log10(p.adjust), fill = Count), stat = "identity") +
    scale_fill_distiller(palette='Reds',direction =1)+
    theme(axis.title.y = element_blank(),legend.position = c(-1.2, 0.85),legend.key.size = unit(0.8, "lines"))
}

df = as.data.frame(go@result)
df = df[order(df$p.adjust,decreasing = F),][1:10,]
df$ID = gsub("GOCC_","",df$ID) %>% gsub("GOBP_","",.) %>% gsub("GOMF_","",.) %>% tolower

png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19/Endo_state_GOenrich.png',width = 1500,height = 1000,res = 300)
EnrichBarPlot(df) + ggtitle('GO')
dev.off()

### microarray gene overlap
colorList = list()
colorList[['Normal_up']] <- c("#FF745A","#5cb85c")
colorList[['Normal_down']] <- c("#007E99","#5cb85c")
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
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19/Endo_enrich.png',width = 1500,height = 600,res = 160)
ggarrange(pList$Normal_up,
            pList$Transition_up,
            pList$Senescent_up,
            pList$Normal_down,
            pList$Transition_down,
            pList$Senescent_down,
            ncol = 3,nrow = 2,legend = "none")
dev.off()

### age state percentage in each group
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19/Endo_state_percent.png',width = 1000,height = 800,res = 250)
ggplot(data = df_plot) +
    geom_bar(aes(x = Group, y = value, fill = age_state), stat = "identity",position = "fill")+
    scale_fill_manual(values = bar) +
    theme_bw()+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = -15, hjust = 0),text = element_text(size = 12))
dev.off()

### age celltyp state percentage in each group
mat = Endo.m@meta.data %$% table(celltype_state,group) %>% as.data.frame.matrix %>% as.matrix
mat = mat/rowSums(mat)
mat = t(scale(mat))
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19/Endo_state_celltype_percent.png',width = 2000,height = 1000,res = 250)
col <- colorRamp2(c(-2,0,2), c("blue","white", "red"), space = "LAB")
ComplexHeatmap::Heatmap(mat,col = col,cluster_rows = F,cluster_columns = T,show_row_names = T,show_column_names = T,heatmap_legend_param = list(title = "Scaled percent"))
dev.off()


























