library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(ggpubr)
library(RColorBrewer)

bar = c("#5cb85c","#428bca","#d9534f")
names(bar) = c("Normal","Transition","Senescent")

### age state
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
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19/Endo_state_percent.png',width = 1200,height = 800,res = 250)
ggplot(data = df_plot) +
    geom_boxplot(aes(x = Group, y = percent, fill = age_state))+
    scale_fill_manual(values = bar) +
    theme_bw()+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = -15, hjust = 0),text = element_text(size = 12))
dev.off()


Frac[rownames(meta),] %>%
    data.frame() %>%
    mutate(Group = meta$group) %>%
    reshape2::melt(id.vars = 'Group',variable.name = 'age_state',value.name = 'percent')%>%
    filter(age_state == 'Senescent') %>%
    ggplot() +
    geom_boxplot(aes(x = Group, y = percent))+
    theme_bw()+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = -15, hjust = 0),text = element_text(size = 12))




### age celltyp state percentage in each group
mat = Endo.m@meta.data %$% table(celltype_state,group) %>% as.data.frame.matrix %>% as.matrix
mat = mat/rowSums(mat)
mat = t(scale(mat))
col <- colorRamp2(c(-2,0,2), c("blue","white", "red"), space = "LAB")

png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19/Endo_state_celltype_percent.png',width = 2000,height = 1000,res = 250)
ComplexHeatmap::Heatmap(mat,col = col,cluster_rows = F,cluster_columns = T,show_row_names = T,show_column_names = T,heatmap_legend_param = list(title = "Scaled percent"))
dev.off()

### lymphatic EC fate
### trajectory plot
cds$LymphaticEC = cds$age_state
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19/Endo_lymfate_trajectory.png',width = 1500,height = 1000,res = 300)
plot_cell_trajectory(cds,color_by="LymphaticEC", size=3,show_backbone=TRUE) + 
    scale_color_manual(values = bar)
dev.off()

### branch gene heatmap
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19/Endo_lymfate_genes.png',width = 1500,height = 1800,res = 300)
plot_genes_branched_heatmap(cds[genes_cluster$Cluster_1,],
                            branch_point = 1,
                            num_clusters = 1,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T,
                            hmcols = colorRampPalette(rev(brewer.pal(9, "Spectral")))(62),
                            branch_colors = as.character(bar),
                            branch_labels = c("Cell fate 1", "Cell fate 2"))
dev.off()

### erich go pathway 
df_plot = cluster_enrich@result
df_plot$class = substr(df_plot$ID,1,4)
df_plot = filter(df_plot,qvalue < 0.05)
df_plot = df_plot[order(df_plot$class,df_plot$qvalue,decreasing = F),] %>% group_by(class) %>% do(head(., n = 10)) %>% ungroup()
df_plot
tolower(df_plot$ID)

### middle group state percentage
png('/home/wangjing/wangj/codebase/HUSI/Figures/covid-19/Endo_lymfate_percent.png',width = 1200,height = 500,res = 300)
cds@phenoData@data %>% 
    # filter(group %in% c('middle_moderate','middle_severe')) %>%
    mutate(Fate = ifelse(State == '1','Root',ifelse(State == '2','Cell fate 1','Cell fate 2'))) %>%
    ggplot()+
    geom_bar(aes(y = group,fill = Fate),stat = 'count',position = 'fill')+
    theme_classic()+
    theme(text = element_text(size = 16),
          axis.title = element_blank(),
          legend.position = c(0.15,1.25),legend.direction = 'horizontal',
          plot.margin=unit(c(2,1,1,1),'lines'))+
    scale_fill_manual(values = brewer.pal(3, "Reds"))
dev.off()


















