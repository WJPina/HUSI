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


##################### train set counts
### histogram
png("HUSI/Figures/melanome/hist.png",width = 1200,height = 1200,res = 300)
data.frame(Counts = colSums(exprData[,metadata$sample_title]),Dataset = metadata$study_accession) %>%
ggplot()+
  geom_histogram(aes(x=Counts),bins = 50,color = '#adb5bd',linewidth = 0.4)+
  # scale_fill_manual()+
  theme_classic()+
  theme(text = element_text(size = 14),legend.position = c(0.25,0.7))+
  ylab('Count')
dev.off()

##################### senescent marker expression in train set
# genes = c('CDKN1A','CDKN2B')
genes = c('CDKN1A','CDKN2B','SERPINE1' )
df_plot = data.frame(t(cbind(X_bk[genes,],X_tr[genes,])),Condition = c(rep('other',122),rep('senescent',126))) %>%
            melt(value.name = 'expression',variable.name = 'marker',id.var = 'Condition') 
pwc <- df_plot %>% group_by(marker) %>% 
        rstatix::pairwise_t_test(expression ~ Condition, paired = F,p.adjust.method = "bonferroni")
pwc
pwc <- pwc %>% rstatix::add_xy_position(x = "marker")
# png('HUSI/Figures/model/valid_markers.png',width = 1200,height = 1500,res = 300)
pdf('HUSI/Figures/model/valid_markers.pdf',width = 4,height = 5)
df_plot %>%
  ggplot(aes(x = marker, y = expression)) + 
  geom_boxplot(aes(color = Condition),outlier.shape = NA) + 
  geom_point(aes(color = Condition),position = position_jitterdodge(jitter.width = 0.2),size = 1) +
  stat_pvalue_manual(pwc)+
  theme_classic()+
  theme(text = element_text(size = 16),axis.text.x = element_text(size=16),legend.position = 'bottom',axis.title.x = element_blank())+
  scale_color_manual(values = c('senescent' = '#e63946','other'='#023e8a'))+
  ylab('Normalized Expression')
dev.off()

##################### gene expression and weights
res <- data.frame(expression = rowMeans(X[names(mm_l2$w),]),weight = mm_l2$w,gene = names(mm_l2$w),row.names = names(mm_l2$w))
res$threshold = ifelse(res$weight > 0,"Positive","Negative")
label_data <- subset(res, abs(res$weight) > 0.02)
pdf("HUSI/Figures/model/model_genes.pdf",width = 8,height = 6)
ggplot(data = res, aes(y=expression,x=weight,color=threshold)) +
  geom_point(alpha=0.5, size=2)+
  geom_text_repel(data = label_data,label = rownames(label_data),color = 'black',max.overlaps = 10,size = 3,show.legend = F)+
  scale_color_manual(values = c("Positive" = "#FF745A","Negative"="#007E99"))+
  ylab('Mean normalized expression')+
  xlab('Weight in model')+
  theme_classic()+
  theme(text = element_text(size = 16),legend.position = 'none')
dev.off()


genes = c('CDKN1A','CDKN2B','SERPINE1' )
top_pos_genes = sort(mm_l2$w,decreasing = T)[1:10] %>% names()

mat = X_tr[c(genes,top_pos_genes),]

left_anno = rowAnnotation(df = data.frame(genes= rep(c("Canonical markers","Top weighted genes"),times=c(3,10))),show_legend = T,
                          col = list(genes = c('Canonical markers' = '#023e8a','Top weighted genes' = '#e63946')),show_annotation_name = F)
col <- colorRamp2(c(-1.5,0,1.5), c("#023e8a","white", "#e63946"), space = "LAB")

png('/home/wangjing/wangj/codebase/HUSI/Figures/revison',)
Heatmap(mat,
        show_column_names = T,
        show_row_names = T,
        row_title = NULL,
        col = col,
        cluster_rows = T,
        cluster_row_slices = FALSE,
        cluster_columns = F,
        left_annotation = left_anno,
        row_names_gp = gpar(fontsize = 12),
        heatmap_legend_param = list(title = "Scaled\nexpression"),
        row_split = factor(rep(c("Canonical markers","Top weighted genes"),times=c(3,10))))

############################ GSEA
library(enrichplot)
library(gggsea)
source('/home/wangjing/wangj/codebase/HUSI/functions.R')

cols = c(rev(colorRampPalette(c("transparent", 'red'))(6))[1:5],colorRampPalette(c("transparent", 'blue'))(5)[c(2:5)])
names(cols) = df$ID
# png('HUSI/Figures/model/valid_GSEA_DOWN.png',width = 1600,height = 1500,res = 300)
pdf('HUSI/Figures/model/valid_GSEA_UP.pdf',width = 6,height = 5)
mygseaplot2(fgsea,
            # geneSetID = df$ID[6:9],
            geneSetID = df$ID[1:5],
            # title = "Negatively enriched gene sets",
            title = "Positively enriched gene sets",
            # color= cols[6:9],
            color = cols[1:5],
            base_size = 12,
            rel_heights = c(1, 0.2, 0.4),
            subplots = 1:3,
            pvalue_table = FALSE,
            ES_geom = "line"
)
dev.off()

### plot dot for enrichment results
df = rbind(df_hallmarker,df_kegg,df_reactome)
dim(df)
write.csv(df,'/home/wangjing/wangj/codebase/HUSI/fgsea.csv',row.names = F)

### only keep top 20 positive and negative gene sets of Reactome
kegg_sene = c('KEGG_Cytokines and growth factors', 'KEGG_Chemokine signaling pathway', 'KEGG_Lysosome', 'KEGG_IL-17 signaling pathway', 'KEGG_Cytokine-cytokine receptor interaction', 'KEGG_Cell adhesion molecules', 'KEGG_Natural killer cell mediated cytotoxicity', 'KEGG_Wnt signaling pathway', 
              'KEGG_DNA repair and recombination proteins', 'KEGG_Signaling pathways regulating pluripotency of stem cells', 'KEGG_Cell cycle', 'KEGG_DNA replication', 'KEGG_DNA replication proteins', 'KEGG_ATP-dependent chromatin remodeling', 'KEGG_Fanconi anemia pathway', 'KEGG_Homologous recombination')
reactome_sene <-c ('REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES', 'REACTOME_INTERFERON_ALPHA_BETA_SIGNALING', 'REACTOME_INTERLEUKIN_10_SIGNALING', 'REACTOME_ONCOGENIC_MAPK_SIGNALING', 'REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM', 'REACTOME_INNATE_IMMUNE_SYSTEM', 'REACTOME_MAPK_FAMILY_SIGNALING_CASCADES', 
                   'REACTOME_DNA_REPAIR', 'REACTOME_CELL_CYCLE_MITOTIC', 'REACTOME_DNA_REPLICATION', 'REACTOME_CELL_CYCLE_CHECKPOINTS', 'REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR', 'REACTOME_CHROMOSOME_MAINTENANCE', 'REACTOME_CELL_CYCLE', 'REACTOME_G2_M_CHECKPOINTS', 'REACTOME_ACTIVATION_OF_THE_PRE_REPLICATIVE_COMPLEX', 'REACTOME_DNA_STRAND_ELONGATION')

df = rbind(df_hallmarker,
           df_kegg[df_kegg$Description %in% kegg_sene,],
           df_reactome[df_reactome$Description %in% reactome_sene,])
dim(df)

p1 <- ggplot(df[df$NES>0,],aes(x = NES,y=reorder(toupper(Description),NES),size = -log10(pvalue)))+
  geom_point(color = "#FF3F3F")+
  theme_bw()+
  theme(text = element_text(size = 16))+
  xlab('Normalized Enrichment Score (NES)')+
  ylab('Pathway')+
  ggtitle('Positively enriched gene sets')+
  ggforce::facet_col(~Database, scales = 'free', space = 'free')

p2 <- ggplot(df[df$NES<0,],aes(x = NES,y=reorder(toupper(Description),-NES),size = -log10(pvalue)))+
  geom_point(color = "#3F3FFF")+
  theme_bw()+
  theme(text = element_text(size = 16))+
  xlab('Normalized Enrichment Score (NES)')+
  ylab('Pathway')+
  ggtitle('Negatively enriched gene sets')+
  ggforce::facet_col(~Database, scales = 'free', space = 'free')

png('/home/wangjing/wangj/codebase/HUSI/Figures/revison/valid_GSEA_dot.png',width = 3500,height = 4200,res = 300)
plot_grid(p1,p2,ncol = 1,nrow = 2,rel_widths = c(1,1.5))
dev.off()


### senescence pathways
pdf('/home/wangjing/wangj/codebase/HUSI/Figures/revison/valid_GSEA_senePathways_2.pdf',width = 15)
pheatmap(data.frame(NES=df$NES,row.names = rownames(df)),show_rownames = T,
         display_numbers = T,number_color = "black",cellwidth = 25,legend = F)
dev.off()


png('/home/wangjing/wangj/codebase/HUSI/Figures/revison/valid_GSEA_SenoMyo.png',width = 5000,height = 400,res = 300)
weights = sort(mm_l2$w[filter(sene_pathways,term == 'SAUL_SEN_MAYO')$gene],decreasing = T)
pheatmap(data.frame('weightx100' = weights[1:50]*100) %>% t,cluster_cols = F,
         display_numbers = T,number_format = '%.2f',number_color = "black",
         legend = F,color = colorRampPalette(c("white", "red"))(50))
dev.off()

### permutaion test
p <- 
data.frame(Ratio=Ratio, times=1:length(Ratio)) %>%
ggplot(aes(x=Ratio)) + 
  geom_histogram(bins = 20,fill='grey',color='black')+
  geom_vline(xintercept =0.6595745,colour="red")+
  theme_classic()+ggtitle('Ratios of permutations (times=1000)')

png('/home/wangjing/wangj/codebase/HUSI/Figures/revison/valid_GSEA_permutation.png',width = 2000,height = 1500,res = 300)
p+ ggbreak::scale_x_break(c(0.4,0.6),space = 0.3)
dev.off()

######################### weights preference ###################################
### pheatmap markers
groups = metadata$celltype
names(groups) = metadata$sample_title
my_order <- sort(table(groups),decreasing = F)
groups <- groups[order(factor(groups,levels = names(my_order)))]
groups <- factor(groups,levels = names(my_order))
top_anno <- HeatmapAnnotation(df = metadata[,c('condition','celltype')],show_legend = T,show_annotation_name = F)
col <- colorRamp2(c(-1.5,0,1.5), c("#023e8a","white", "#e63946"), space = "LAB")

genes = unlist(MarkerSets$HFF_fibroblast)
genes = genes[!is.na(genes)]
genes = genes[genes %in% rownames(X)]
mat = X_centre[genes,names(groups)] %>% t
mat = aggregate(mat,by = list(celltype = groups),mean) %>% data.frame() %>% column_to_rownames('celltype') %>% t

Heatmap(mat,
        show_column_names = T,
        show_row_names = T,
        row_title = NULL,
        col = col,
        cluster_rows = F,
        cluster_row_slices = FALSE,
        cluster_columns = F,
        row_names_gp = gpar(fontsize = 12),
        heatmap_legend_param = list(title = "Scaled\nexpression"))


### enrichment
df$set_type = ifelse(df$Description %in% c('SenUp','SenMayo'),'Senescence','Cell type')
df$Enrichment = ifelse(df$pvalue<=0.05,'Significant','non-significant')
write.csv(df,file='/home/wangjing/wangj/codebase/HUSI/Figures/revison/valid_GSEA_celltype_1.csv')

png('/home/wangjing/wangj/codebase/HUSI/Figures/revison/valid_GSEA_celltype_1.png',width = 2000,height = 1200,res = 300)
ggplot(df,aes(x = NES,y=reorder(Description,NES),size = -log10(pvalue),color = Enrichment,shape=set_type))+
  geom_point()+
  scale_color_manual(values = c('#1d3557','#e63946'))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  xlab('Normalized Enrichment Score (NES)')+
  ylab('Pathway')+
  geom_vline(xintercept = 0)
dev.off()

### distribution
df_plot = do.call(cbind, lapply(lapply(all_sets_weight, unlist), `length<-`, max(lengths(all_sets_weight)))) %>% data.frame
df_plot = melt(df_plot)
colnames(df_plot) <- c('Set','Weight')
df_plot = df_plot[!is.na(df_plot$Weight),]
write.csv(df_plot,file='/home/wangjing/wangj/codebase/HUSI/Figures/revison/valid_GSEA_celltype_2.csv')

df_plot$type = factor(ifelse(df_plot$Set %in% names(SenSet),'SenSet','CelltypeSet'),
                      levels = c('SenSet','CelltypeSet'),ordered = T)
sdf = aggregate(df_plot$Weight,by=list(df_plot$Set),mean) %>% data.frame()
colnames(sdf) <- c('Set','mean')
sdf$sd = aggregate(df_plot$Weight,by=list(df_plot$Set),sd)[,'x']
sdf = sdf[order(sdf$mean,decreasing = T),]
df_plot$Set <- factor(df_plot$Set,levels = sdf$Set,ordered = T)

png('/home/wangjing/wangj/codebase/HUSI/Figures/revison/valid_GSEA_celltype_2.png',width = 2000,height = 1500,res = 300)
ggplot(df_plot,aes(x=Set,y=Weight,color=type))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = 'top',text = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1))+
  scale_color_manual(values = c('#e63946','#1d3557'))
dev.off()

### validation on microarray
df_IMR90 <- 
  ArrayList [["GSE19864"]][[1]] %>% 
  pData %>% 
  dplyr::select(title, geo_accession) %>% 
  mutate(condition=case_when(
    grepl("Growing", title, ignore.case = T) ~ "Growing", 
    grepl("Senescent", title, ignore.case = T) ~ "Senescent", 
    grepl("Confluent", title, ignore.case = T) ~ "confluent", 
    grepl("0.1% FBS", title, ignore.case = T) ~ "low serum"
  )) %>% 
  .[which(.$condition %in% c("Growing","Senescent")),] %>%
  mutate(Condition = factor(condition,levels = c("Growing","Senescent"),ordered = T)) %>% 
  inner_join(
    data.frame(sene_score=ScoreList[["s_19864"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession") 
df_IMR90$celltype = 'IMR90_fibroblast'

df_Melanocyte <- 
ArrayList [["GSE83922"]][[1]] %>% 
  pData %>% 
  mutate(condition= `cell phenotype:ch1`) %>% 
  inner_join(
    data.frame(sene_score=ScoreList[["s_83922"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("young","senescent")),] %>% 
  mutate(Condition = factor(ifelse(condition == 'young','Young','Senescent'),levels = c("Young","Senescent"),ordered = T)) 
df_Melanocyte$celltype<- 'Melanocyte'


df_plot <- rbind(df_IMR90[,c('Condition','sene_score','celltype')],df_Melanocyte[,c('Condition','sene_score','celltype')])
colnames(df_plot) <- c('Condition','hUSI','celltype')
df_plot$sample <- paste(df_plot$Condition,df_plot$celltype,sep = '_')
df_plot$hUSI <- minmax(df_plot$hUSI)

png('/home/wangjing/wangj/codebase/HUSI/Figures/revison/valid_GSEA_celltype_3.png',width = 2000,height = 2500,res = 300)
ggplot(df_plot,aes(x=sample,y=hUSI,color=celltype))+
  geom_boxplot(outlier.color = NA)+ 
  theme_bw()+
  theme(legend.position = 'top',text = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1))+
  scale_color_manual(values = c('#e63946','#1d3557'))
dev.off()



########################### model leave-one-out auc
auc = sort(auc)

png('HUSI/Figures/model/cross-validation.png',width = 1000,height = 800,res = 300)
# pdf('HUSI/Figures/model/cross-validation.pdf',width = 5,height = 4)
ggplot(aes(x = 1:length(auc),y = auc),data=data.frame(auc=auc))+
    geom_line()+
    geom_hline(yintercept = mean(auc),linetype = 'dashed',color = 'red')+
    geom_text(aes(x = 15,y = mean(auc),label = 'mean: 0.9'),vjust = -1)+
    theme_classic()+
    theme(text=element_text(size=12),axis.title.x = element_blank())+
    ggtitle('LOOCV')+
    ylab("AUC")
dev.off()

################################## bacth effect
# pdf('HUSI/Figures/model/valid_BE.pdf',width = 5,height = 5)
p.be <- 
read.table("/mnt/data1/wangj/AgingScore/Data/Bulk_BatchEffect/batch_IMR90_4OHT_add_condition_tsv.txt",sep = "\t",header = T) %>%
  dplyr::select(c("title","study_accession","condition")) %>%
  mutate(hUSI = s_batch[title]) %>%
  column_to_rownames("title") %>%
  mutate(Condition = factor(ifelse(condition == 'other','Other','Senescent'),levels = c("Other","Senescent"),ordered = T)) %>%
  ggplot(aes(Condition, hUSI)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = study_accession),width = 0.2) +
  scale_colour_jama()+
  geom_line(aes(group = study_accession), linetype="dashed", col="skyblue") + 
  geom_signif(comparisons = list(c("Other","Senescent")),test = "t.test",test.args = c("less"),map_signif_level = T) + 
  ylab('hUSI')+
  ggtitle("OIS in IMR90")+
  mytheme() +ylim(0,1.3) + theme(legend.position = 'right')

################################## valid data in RNA-seq
mytheme <- function () { 
  theme_classic() %+replace% 
    theme(text = element_text(size = 16),legend.position = "none",
          axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))}

p1.1 <- data.frame(hUSI = s_CIS) %>% 
  mutate(condition= gsub("_.*", "", names(s_CIS))) %>%
  mutate(Treatment = factor(condition,levels = c("Immortal","Adria","H2O2","5-aza"),ordered = T)) %>% 
  ggplot(aes(Treatment, hUSI)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Treatment),width = 0.2,height = 0,size=2)+
  geom_signif(comparisons = list(c("Immortal","H2O2"),c("Immortal","5-aza")),
              test = "t.test",
              step_increase=0.1,
              map_signif_level = T,
              test.args = c("less")) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle("CIS in MDAH041")+ ylim(0,1.3)

p1.2 <- data.frame(hUSI = s_OIS) %>% 
  mutate(condition= gsub("_.*", "", rownames(.))) %>%
  mutate(Treatment = factor(condition,levels = c("OISD0","OISD2","OISD4","OISD6",'OISD10'),ordered = T)) %>% 
  ggplot(aes(Treatment, hUSI)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Treatment),width = 0.3,height = 0,size=2)+
  geom_signif(comparisons = list(c("OISD0","OISD2"),c("OISD2","OISD4"),c("OISD4","OISD6")),
              test = "t.test",
              step_increase=0.1,
              map_signif_level = T,
              test.args = c("less")) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle("OIS in WI-38")+ ylim(0,1.3)

p1.3 <- data.frame(hUSI = s_RS) %>% 
  mutate(condition= gsub('_.*','',gsub("RS_", "", rownames(.)))) %>%
  mutate(Condition = factor(condition,levels = c("Proliferative","Senescent"),ordered = T)) %>% 
  ggplot(aes(Condition, hUSI)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Condition),width = 0.2,height = 0,size=2)+
  geom_signif(comparisons = list(c("Proliferative","Senescent")),
              test = "t.test",
              y_position=1.1,
              map_signif_level = T,
              test.args = c("less")) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle("RS in WI-38")+ ylim(0,1.3)

fig <- plot_grid(p.be,p1.1,p1.2,p1.3,ncol = 4,nrow = 1, align = "h",rel_widths = c(6,4,5,3))
png('/home/wangjing/wangj/codebase/HUSI/Figures/revison/valid_RNA-seq.png',width = 4200,height = 1800,res = 300)
# pdf('HUSI/Figures/model/valid_RNA-seq.pdf',width = 10,height = 5)
fig
dev.off()
############################# valid data in micro-array
mytheme <- function () { 
  theme_classic() %+replace% 
    theme(text = element_text(size = 16),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "right",
          axis.text.x = element_blank())
}

names(ArrayList)
p2.1 <- ArrayList [["GSE19864"]][[1]] %>% 
  pData %>% 
  dplyr::select(title, geo_accession) %>% 
  mutate(condition=case_when(
    grepl("Growing", title, ignore.case = T) ~ "Growing", 
    grepl("Senescent", title, ignore.case = T) ~ "Senescent", 
    grepl("Confluent", title, ignore.case = T) ~ "confluent", 
    grepl("0.1% FBS", title, ignore.case = T) ~ "low serum"
  )) %>% 
  .[which(.$condition %in% c("Growing","Senescent")),] %>%
  mutate(Condition = factor(condition,levels = c("Growing","Senescent"),ordered = T)) %>% 
  inner_join(
    data.frame(sene_score=ScoreList[["s_19864"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  ggplot(aes(Condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Condition),width = 0.1,size=4)+
  geom_signif(comparisons = list(c("Growing", "Senescent")),test = "t.test",test.args = c("less"), map_signif_level = T) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle("OIS in IMR90")+
  scale_colour_manual(values = c('#023e8a','#e63946'))

p2.2 <- ArrayList [["GSE16058"]][[1]] %>% 
  pData %>% 
  mutate(condition=gsub("growth status: ", "", characteristics_ch1.3) ) %>% 
  mutate(passage = `passage:ch1`) %>%
  inner_join(
    data.frame(sene_score=ScoreList[["s_16058"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("Growing","Senescent")),] %>% 
  mutate(Condition = factor(condition,levels = c("Growing","Senescent"),ordered = T)) %>% 
  mutate(Passage = passage) %>% 
  ggplot(aes(Condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Condition,shape = Passage),width = 0.1,size=4) +
  geom_signif(comparisons = list(c("Growing", "Senescent")),test = "t.test",test.args = c("less"), map_signif_level = T) + 
  mytheme()+
  guides(color = 'none')+
  ylab('hUSI')+
  ggtitle('RS in HMEC')+
  scale_colour_manual(values = c('#023e8a','#e63946'))

p2.3 <- ArrayList [["GSE83922"]][[1]] %>% 
  pData %>% 
  mutate(condition= `cell phenotype:ch1`) %>% 
  inner_join(
    data.frame(sene_score=ScoreList[["s_83922"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("young","senescent")),] %>% 
  mutate(Condition = factor(ifelse(condition == 'young','Young','Senescent'),levels = c("Young","Senescent"),ordered = T)) %>%
  ggplot(aes(Condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Condition),width = 0.1,size=4) +
  geom_signif(comparisons = list(c("Young", "Senescent")),test = "t.test",test.args = c("less"), map_signif_level = T) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle('RS in Melanocyte')+
  scale_colour_manual(values = c('#023e8a','#e63946'))

p2.4 <- ArrayList [["GSE11954"]][[1]] %>% 
  pData %>% 
  mutate(condition=case_when(
    grepl("growing", description, ignore.case = T) ~ "Growing", 
    grepl("senescent", description, ignore.case = T) ~ "Senescent")) %>% 
  inner_join(
    data.frame(sene_score=ScoreList[["s_11954"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("Growing","Senescent")),] %>% 
  mutate(Condition = factor(condition,levels = c("Growing","Senescent"),ordered = T)) %>% 
  ggplot(aes(Condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Condition),width = 0.1,size=4) +
  geom_signif(comparisons = list(c("Growing","Senescent")),test = "t.test", map_signif_level = T) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle('CIS in HSC')+
  scale_colour_manual(values = c('#023e8a','#e63946'))

p2.5 <- ArrayList [["GSE100014"]][[1]] %>% 
  pData %>% 
  mutate(condition=case_when(
    grepl("Proliferating", title, ignore.case = T) ~ "Proliferating", 
    grepl("senescence", title, ignore.case = T) ~ "Senescent")) %>% 
  inner_join(
    data.frame(sene_score=ScoreList[["s_100014"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("Proliferating","Senescent")),] %>% 
  mutate(Condition = factor(condition,levels = c("Proliferating","Senescent"),ordered = T)) %>% 
  ggplot(aes(Condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Condition),width = 0.1,size=4) +
  geom_signif(comparisons = list(c("Proliferating","Senescent")),test = "t.test",test.args = c("less"), map_signif_level = T) + 
  mytheme()+
  ylab('hUSI')+
  ggtitle('CIS in HBEC')+
  scale_colour_manual(values = c('#023e8a','#e63946'))

p2.6 <- ArrayList [["GSE77239"]][[1]] %>% 
  pData %>% 
  mutate(condition=case_when(
    grepl("young", `cells:ch1`, ignore.case = T) ~ "young", 
    grepl("old", `cells:ch1`, ignore.case = T) ~ "old")) %>% 
  mutate(Condition = factor(condition,levels = c("young","old"),ordered = T)) %>% 
  rename("Treatment"="treatment:ch1") %>% 
  inner_join(
    data.frame(sene_score=ScoreList[["s_77239"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  ggplot(aes(Condition, sene_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Condition,shape = Treatment),width = 0.1,size=4) +
  geom_line(aes(group = Treatment), linetype="dashed", col="skyblue") + 
  geom_signif(comparisons = list(c("young","old")),test = "t.test",test.args = c("less"), map_signif_level = T) + 
  mytheme()+
  guides(color = 'none')+
  ylab('hUSI')+
  ggtitle('RS in HCAEC')+
  scale_colour_manual(values = c('#023e8a','#e63946'))

fig <- ggarrange(p2.1,p2.2,p2.3,p2.4,p2.5,p2.6,ncol = 3,nrow = 2,common.legend = F)+ 
        theme(plot.margin = unit(c(0,0,0,0), "cm"))
# png('HUSI/Figures/model/valid_array.png',width = 4500,height = 2400,res = 300)
pdf('HUSI/Figures/model/valid_array.pdf',width = 15,height = 8)
fig
dev.off()

################################# plot lost rate stability
mytheme <- function () { 
  theme_classic() %+replace% 
    theme(text = element_text(size = 16))
}

p2.1 <- ArrayList [["GSE19864"]][[1]] %>% 
  pData %>% 
  dplyr::select(title, geo_accession) %>% 
  mutate(condition=case_when(
    grepl("Growing", title, ignore.case = T) ~ "Growing", 
    grepl("Senescent", title, ignore.case = T) ~ "Senescent", 
    grepl("Confluent", title, ignore.case = T) ~ "confluent", 
    grepl("0.1% FBS", title, ignore.case = T) ~ "low serum"
  )) %>% 
  .[which(.$condition %in% c("Growing","Senescent")),] %>%
  mutate(condition = factor(ifelse(condition == 'Growing','Other','Senescent'),levels = c("Other","Senescent"),ordered = T)) %>% 
  inner_join(
    data.frame(hUSI=Lost_ScoreList[["GSE19864"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>%
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  mutate(Condition = condition) %>%
  ggplot(aes(x = rate,y = hUSI,color = Condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle("OIS in IMR90")+
  scale_colour_manual(values = c('#023e8a','#e63946'))+
  xlab('Zeroing-out Rate')


p2.2 <- ArrayList [["GSE16058"]][[1]] %>% 
  pData %>% 
  mutate(condition=gsub("growth status: ", "", characteristics_ch1.3) ) %>% 
  mutate(passage = `passage:ch1`) %>%
  inner_join(
    data.frame(hUSI=Lost_ScoreList[["GSE16058"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("Growing","Senescent")),] %>% 
  mutate(condition = factor(ifelse(condition == 'Growing','Other','Senescent'),levels = c("Other","Senescent"),ordered = T)) %>%
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  mutate(Condition = condition) %>%
  ggplot(aes(x = rate,y = hUSI,color = Condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle('RS in HMEC')+
  scale_colour_manual(values = c('#023e8a','#e63946'))+
  xlab('Zeroing-out Rate')

p2.3 <- ArrayList [["GSE83922"]][[1]] %>% 
  pData %>% 
  mutate(condition= `cell phenotype:ch1`) %>% 
  inner_join(
    data.frame(hUSI=Lost_ScoreList[["GSE83922"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("young","senescent")),] %>% 
  mutate(condition = factor(ifelse(condition == "young","Other","Senescent"),levels = c("Other","Senescent"),ordered = T)) %>%
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  mutate(Condition = condition) %>%
  ggplot(aes(x = rate,y = hUSI,color = Condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle('RS in Melanocyte')+
  scale_colour_manual(values = c('#023e8a','#e63946'))+
  xlab('Zeroing-out Rate')


p2.4 <- ArrayList [["GSE11954"]][[1]] %>% 
  pData %>% 
  mutate(condition=case_when(
    grepl("growing", description, ignore.case = T) ~ "Growing", 
    grepl("senescent", description, ignore.case = T) ~ "Senescent")) %>% 
  inner_join(
    data.frame(hUSI=Lost_ScoreList[["GSE11954"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("Growing","Senescent")),] %>% 
  mutate(condition = factor(ifelse(condition == 'Growing','Other','Senescent'),levels = c("Other","Senescent"),ordered = T)) %>% 
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  mutate(Condition = condition) %>%
  ggplot(aes(x = rate,y = hUSI,color = Condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle('CIS in HSC')+
  scale_colour_manual(values = c('#023e8a','#e63946'))+
  xlab('Zeroing-out Rate')


p2.5 <- ArrayList [["GSE100014"]][[1]] %>% 
  pData %>% 
  mutate(condition=case_when(
    grepl("Proliferating", title, ignore.case = T) ~ "Proliferating", 
    grepl("senescence", title, ignore.case = T) ~ "Senescent")) %>% 
  inner_join(
    data.frame(hUSI=Lost_ScoreList[["GSE100014"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>% 
  .[which(.$condition %in% c("Proliferating","Senescent")),] %>% 
  mutate(condition = factor(ifelse(condition == 'Proliferating',"Other","Senescent"),levels = c("Other","Senescent"),ordered = T)) %>%
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  mutate(Condition = condition) %>%
  ggplot(aes(x = rate,y = hUSI,color = Condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle('CIS in HBEC')+
  scale_colour_manual(values = c('#023e8a','#e63946'))+
  xlab('Zeroing-out Rate')


p2.6 <- ArrayList [["GSE77239"]][[1]] %>% 
  pData %>% 
  mutate(condition=case_when(
    grepl("young", `cells:ch1`, ignore.case = T) ~ "young", 
    grepl("old", `cells:ch1`, ignore.case = T) ~ "old")) %>% 
  mutate(condition = factor(ifelse(condition == 'young',"Other","Senescent"),levels = c("Other","Senescent"),ordered = T)) %>% 
  rename("treatment:ch1"="treatment") %>% 
  inner_join(
    data.frame(hUSI=Lost_ScoreList[["GSE77239"]]) %>% rownames_to_column("geo_accession"), 
    by="geo_accession"
  ) %>%
  .[c('condition',paste('hUSI',names(Lost_ScoreList[[1]]),sep = '.'))] %>% 
  reshape2::melt(id.vars = 'condition',variable.name ='rate',value.name = 'hUSI') %>%
  mutate(rate = as.numeric(gsub('hUSI\\.rate_','',rate))) %>%
  mutate(Condition = condition) %>%
  ggplot(aes(x = rate,y = hUSI,color = Condition)) + 
  geom_smooth()+
  mytheme()+
  ggtitle('RS in HCAEC')+
  scale_colour_manual(values = c('#023e8a','#e63946'))+
  xlab('Zeroing-out Rate')

fig <- ggarrange(p2.1,p2.2,p2.3,p2.4,p2.5,p2.6,ncol = 3,nrow = 2,common.legend = T,legend = 'bottom')
# png('HUSI/Figures/model/valid_lost.png',width = 4500,height = 2400,res = 300)
pdf('HUSI/Figures/model/valid_lost.pdf',width = 15,height = 8)
fig
dev.off()

########################### hSUI validate ######################################
### value
df_plot_list = lapply(Results[c("Aarts2017","Georgilis2018")], function(x){data.frame(hUSI = x$hUSI,Condition = x$method$Condition)})
df_plot_list[['Teo2019']] = data.frame(hUSI = Results$Teo2019$hUSI,Condition = gsub("[0-9]$", "", Results$Teo2019$method$Condition1))
df_plot_list[["Tang2019"]] <- data.frame(hUSI=unlist(lapply(Results$Tang2019,function(x){x$method$hUSI})),Condition = rep(c('Growing','Senescence'),each=400,time = 2))
df_plot_list <- lapply(names(df_plot_list), function(x){df_plot_list[[x]]$Dataset = rep(x,nrow(df_plot_list[[x]]));return(df_plot_list[[x]])})
df_plot = data.table::rbindlist(df_plot_list) %>% data.frame()
pwc <- df_plot %>% group_by(Dataset) %>% 
        rstatix::pairwise_t_test(hUSI ~ Condition, paired = F,p.adjust.method = "bonferroni")
pwc
pwc <- pwc %>% rstatix::add_xy_position(x = "Dataset")

# png('HUSI/Figures/model/compare_hUSI.png',width = 1400,height = 1500,res = 300)
pdf('HUSI/Figures/model/compare_hUSI.pdf',width = 7,height = 7)
df_plot %>%
  mutate(Condition = factor(Condition,levels = c('Growing','Senescence'),ordered = T)) %>%
  ggplot(aes(x = Dataset,y=hUSI,color = Condition)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_point(aes(color = Condition),position = position_jitterdodge(jitter.width = 0.3),size = 0.5) +
    scale_color_manual(values = c('#023e8a','#e63946'))+
    stat_pvalue_manual(pwc)+
    theme_classic()+
    theme(legend.position = 'top',text = element_text(size = 16),axis.text.x = element_text(angle = 30,vjust = 1, hjust=1))
dev.off()

################################################################################
######################## comparison score ######################################
################################################################################
### Teo dataset
exp = Teo2019@assays$RNA@data[c('CDKN1A','CDKN2B','IL1B','IL6','IL8'),] %>% data.frame() %>% t %>% 
  data.frame() %>% mutate(cell = colnames(Teo2019),hUSI=hUSI[colnames(Teo2019)]) %>% melt() %>% setnames(c('cell','gene','expression'))

cbind(exp,cluster=Teo2019$group[exp$cell]) %>%
  ggplot(aes(y=expression,x = cluster,fill = cluster)) +
  geom_boxplot(outlier.shape = NA) + 
  theme_classic()+
  theme(legend.position = 'top',text = element_text(size = 16),axis.text.x = element_text(angle = 30,vjust = 1, hjust=1))+
  facet_wrap(~gene,scales = 'free_y',ncol = 2)
############################### mean auc values ################################
### plot data
df_plot_list = lapply(AUClist[c("Teo2019","Tang2019","Aarts2017")], function(x) {melt(x) %>% set_names(c('method','unit','AUC'))})

plot_list = list()
for (dataset in names(df_plot_list)) {
  df_plot = df_plot_list[[dataset]]
  df_plot$method_type <- ifelse(df_plot$method %in% c("DAS","mSS",'DAS+mSS','lassoCS','CSS'),'Senecence Score',
                          ifelse(df_plot$method %in% names(EnrichSet),'Gene set',
                                 ifelse(df_plot$method %in% SenMarkers,'Marker','hUSI')))
  df_plot$method_type <- factor(df_plot$method_type,levels = c('hUSI','Marker','Gene set','Senecence Score'),ordered = T)
  df_plot$dataset = dataset
  
  plot_list[[dataset]] <-
    ggplot(df_plot,aes(x=reorder(method,AUC), y=AUC,color = method_type))+
      geom_jitter(size=1)+
      geom_boxplot(color='black',outlier.color = NA,fill=NA,linewidth=0.2)+
      theme_classic()+
      ylab('AUC values')+
      xlab('Quantification method')+
      scale_color_manual(values = c('#e63946','#94d2bd','#00afb9','#fed9b7'))+
      theme(panel.grid.major.y = element_line(),
            axis.text.x = element_text(angle = 45,vjust = 1, hjust=1),
            text=element_text(size=16))+
      ggtitle(dataset)+ylim(0,1)
}
  
  
fig <- ggarrange(plotlist = plot_list,ncol = 3,nrow = 1,common.legend = T)

png('~/wangj/codebase/HUSI/Figures/revison/compare_auc.png',width = 6000,height = 2500,res = 300)
fig
dev.off()

############################### mean auc rank ##################################
### only keep common genes
method_used <- c('hUSI','GLB1','TP53','CDKN1A','CDKN2A','LMNB1','RB1','CDK1','CDK4','CDK6','MKI67','CDKN2B',
                 "DAS","mSS",'DAS+mSS','lassoCS','CSS',
                 names(EnrichSet))
df_plot_list = lapply(AUClist[c("Teo2019","Tang2019","Aarts2017")], function(x) {x = x[method_used,]}) 
df_plot_list = lapply(df_plot_list, function(x) {rank(rowMeans(x))})


df_plot = do.call(cbind, df_plot_list)
df_plot = reshape2::melt(as.matrix(df_plot),value.name = "AUC_rank")
colnames(df_plot) = c('method','dataset','AUC_rank')
aov.mean<-aggregate(df_plot$AUC_rank,by=list(df_plot$method),FUN=mean)
aov.sd<-aggregate(df_plot$AUC_rank,by=list(df_plot$method),FUN=sd)
aov<-data.frame(aov.mean,sd=aov.sd$x)
df_plot_mean <- aov
colnames(df_plot_mean) <- c('method','mean','sd')
df_plot_mean$method_type <- ifelse(df_plot_mean$method %in% c("DAS","mSS",'DAS+mSS','lassoCS','CSS'),'Senecence Score',
                     ifelse(df_plot_mean$method %in% names(EnrichSet),'Gene set',
                            ifelse(df_plot_mean$method %in% SenMarkers,'Marker','hUSI')))
df_plot_mean$group <- factor(df_plot_mean$method_type,levels = c('hUSI','Marker','Gene set','Senecence Score'),ordered = T)

png('~/wangj/codebase/HUSI/Figures/revison/compare_rank.png',width = 4000,height = 2000,res = 300)
ggplot()+
    geom_bar(data=df_plot_mean,aes(x=reorder(method,mean), y=mean,fill = group),stat="identity",position="dodge")+
    geom_jitter(data=df_plot,aes(x = method,y=AUC_rank,color=dataset),width = 0.3,size=2,height = 0) +
    geom_errorbar(data=df_plot_mean,aes(x=reorder(method,mean),y=mean,ymax=mean+sd,ymin=ifelse(mean-sd <0,0,mean-sd)),
                  position=position_dodge(0.9),width=0.4)+
    theme_classic()+
    ylab('AUC rank')+
    xlab('Quantification method')+
    scale_fill_manual(values = c('#e63946','#94d2bd','#00afb9','#fed9b7'))+
    scale_color_manual(values = c('#277da1','#f9c74f','#7209b7'))+
    theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1),
          text=element_text(size=16))
dev.off()

#################################### GTEx/TCGA #################################
### gtex
TCSER = read.csv('CS_score_of_single_cell_datasets/Senescence_quantification_GTEX.csv')
TCSER = TCSER[!duplicated(TCSER$Sample),]
rownames(TCSER) = TCSER$Sample
df_plot <- scoreList_GTEx[c('TCSER','hUSI')] %>% data.frame() 
df_plot$Group <- TCSER[names(scoreList_GTEx$hUSI),'primary_site']

### tcga
TCSER = read.csv('CS_score_of_single_cell_datasets/Senescence_quantification_TCGA.csv')
TCSER <- column_to_rownames(TCSER,'sample')
rownames(TCSER) <- gsub('-','.',rownames(TCSER))
df_plot <- scoreList_TCGA[c('TCSER','hUSI')] %>% data.frame() 
df_plot$Group <- TCSER[names(scoreList_TCGA$hUSI),'cancertype']

png('/home/wangjing/wangj/codebase/HUSI/Figures/revison/TCGA_TCSER.png',width =4000,height = 3000,res = 300)
# pdf('HUSI/Figures/model/GTEx_TCSER.pdf',width = 4.5,height = 4)
  df_plot %>%
  ggscatter(x = "hUSI", y = "TCSER",
            # color = '#2ec4b6',
            color = '#f07167',
            add = "reg.line", 
            conf.int = TRUE,
            size = 1,
            add.params = list(color = "#6c757d"),
            ggtheme = theme_classic())+ 
  # stat_cor(method = "spearman",label.x = -.08, label.y = 1,color='black')+
  # stat_cor(method = "spearman",label.x = -.09, label.y = 1.15,color='black')+
  # ggtitle('7123 samples in GTEx')+
  ggtitle("10495 samples in TCGA")+
  theme_classic()+
  theme(text = element_text(size = 12))+
    facet_wrap(~Group,scales = 'free')
dev.off()

### split cancer type or tissue
coeffient = unlist(lapply(split(df_plot[,c('TCSER','hUSI')],df_plot$Group),FUN = function(x) cor(x[,1],x[,2],method = 'sp'))) %>% sort()
hUSI_mat =  unlist(lapply(split(df_plot[,'hUSI'],df_plot$Group),FUN = function(x) mean(x)))[names(cor_mat)]
TCSER_mat =  unlist(lapply(split(df_plot[,'TCSER'],df_plot$Group),FUN = function(x) mean(x)))[names(cor_mat)]

# png('HUSI/Figures/model/GTEx_Tissue.png',width = 3800,height = 1000,res = 300)
# png('HUSI/Figures/model/TCGA_CancerType.png',width = 3900,height = 800,res = 300)
# pdf('HUSI/Figures/model/GTEx_Tissue.pdf',width = 15,height = 3.5)
pdf('HUSI/Figures/model/TCGA_CancerType.pdf',width = 15,height = 2.8)
Heatmap(t(coeffient),cluster_rows = F,cluster_columns = F,
        cell_fun = function(j, i, x, y, w, h, col) {grid.text(round(t(coeffient),2)[i, j], x, y) },
        show_heatmap_legend = F,
        border = T,
        rect_gp = gpar(col = "black", lwd = 1),
        row_title = 'Spearman\nCoeffient',
        row_title_rot = 0,
        top_annotation = HeatmapAnnotation(hUSI = anno_boxplot(split(df_plot[,'hUSI'],df_plot$Group), width = unit(5, "cm"),outline = F, gp = gpar(fill = "#d62828")),
        TCSER = anno_boxplot(split(df_plot[,'TCSER'],df_plot$Group), width = unit(5, "cm"),outline = F, gp = gpar(fill = "#1d3557"))))
dev.off()

### GTEx by Age
df = data.frame(TCSER = scoreList_GTEx$TCSER,Age = GTEx_meta[names(scoreList_GTEx$hUSI),'AGE'],Tissue = GTEx_meta[names(scoreList_GTEx$hUSI),'SMTSD'])
# pdf('HUSI/Figures/model/GTEx_Age.pdf',width = 4.5,height = 4)
png('/home/wangjing/wangj/codebase/HUSI/Figures/revison/GTEx_Age_TCSER.png',width = 1200,height = 1200,res=300)
ggplot(df,aes(x = Age,y = TCSER,fill = Age)) + 
  geom_boxplot() + 
  scale_fill_brewer(palette = "Reds")+
  geom_signif(comparisons = list(c("20-29","30-39"),c("30-39","40-49"),c("40-49","50-59"),c("50-59","60-69"),c("60-69","70-79")),
            test = "t.test",
            step_increase=0.1,
            map_signif_level = T,
            test.args = c("less"))+
  theme_classic()+
  theme(text = element_text(size = 16),legend.position = "none")
dev.off()


################################################################################
################################# melanoma #####################################
################################################################################

### all melanome cell plot
bar = c("#2a9d8f","#e9c46a","#e76f51","#AAC5E2","#FFC592","#A4D99D","#F4A39E","#CCB7DC","#CBAEA8")
names(bar) <- c("Cycling","Transitional","Senescent","B cell","CAF cell","Endo cell","Macro cell","NK cell","T cell")

# ### integrate type
# melanoma_obj$subtype = factor(ifelse(as.character(melanoma_obj$subtype) %in% c("Cycling","Transitional","Senescent"),
#                                      'Tumor cell',as.character(melanoma_obj$subtype)),
#                               levels = c('Tumor cell','B cell','CAF cell','Endo cell','Macro cell','NK cell','T cell'))

# bar = c("#eb0101","#AAC5E2","#FFC592","#A4D99D","#F4A39E","#CCB7DC","#CBAEA8")
# # bar = c("#eb0101",rep('grey',6))
# names(bar) <- c("Tumor cell","B cell","CAF cell","Endo cell","Macro cell","NK cell","T cell")

png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/Melanoma_all_cell_tumor.png',width = 1800,height = 1500,res = 300)
DimPlot(melanoma_obj, group.by = 'subtype',reduction = "tsne",label = T,
        cols = bar,repel = TRUE,label.box = T,label.color = "white",pt.size = 1.5,label.size = 6)+
  ggtitle("Melanoma cell types")+NoLegend()
dev.off()


bar = c("#2a9d8f","#e9c46a","#e76f51")
names(bar) = c("Cycling","Transitional","Senescent")

### histogram
png("/home/wangjing/wangj/codebase/HUSI/Figures/revison/Melanoma_hist.png",width = 1200,height = 1200,res = 300)
# ggplot(EpiExp.m@meta.data)+
filter(df,celltype == 'Tumor') %>% mutate(subpopulation = as.character(subtype)) %>%
ggplot()+
  # geom_histogram(aes(x=hUSI),bins = 50,fill = 'white',color = '#adb5bd',linewidth = 0.4)+
  geom_density(aes(x=hUSI,fill = subpopulation,color = subpopulation),alpha = 0.6)+
  scale_color_manual(values = bar)+
  scale_fill_manual(values = bar)+
  theme_classic()+
  theme(text = element_text(size = 14),legend.position = c(0.25,0.7))+
  ylab('Count')
dev.off()

### state hUSI plot
pwc = EpiExp.m@meta.data %>% 
        rstatix::pairwise_t_test(hUSI ~ State, paired = F,p.adjust.method = "bonferroni") %>%
        rstatix::add_xy_position(x = "State")
png('HUSI/Figures/melanome/Melanoma_hUSI.png',width = 1200,height = 1200,res = 300)
EpiExp.m@meta.data %>%
  ggplot(aes(x = State,y = hUSI,color = State))+
  geom_boxplot()+
  scale_color_manual(values = bar)+
  # labs(caption = get_pwc_label(pwc))+
  stat_pvalue_manual(pwc)+
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.position = 'none')
dev.off()

### proportion of states
png('/home/wangjing/wangj/codebase/HUSI/Figures/revison/melanoma_state_percentage.png',res=300,width = 800,height = 1000)
# melanoma_obj$subtype <- factor(melanoma_obj$subtype,levels = table(melanoma_obj$subtype) %>% sort(decreasing = F) %>% names(),ordered = T)
# ggplot(melanoma_obj@meta.data,aes(x = orig.ident,fill = subtype))+
filter(df,celltype == 'Tumor') %>% mutate(subpopulation = as.character(subtype)) %>%
ggplot(aes(x = celltype,fill = subpopulation))+
  geom_bar(position = 'fill')+
  scale_fill_manual(values = bar)+
  theme_bw()+
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  # xlab('Tumor cells (All)')+
  xlab('Tumor cells (Melanocytes)')+
  ylab('Percentage')
dev.off()

### cor gene volcano
res <- data.frame(expression = rowMeans(EpiExp.m@assays$RNA@data[names(cor.genes),]),correlation = cor.genes,row.names = names(cor.genes))
res$threshold = ifelse(rownames(res) %in% features, ifelse(res$correlation > 0,"Positive","Negative"),"Background")
label_data <- subset(res, abs(res$correlation) > 0.5 & res$expression  > 4)
png("HUSI/Figures/melanome/Melanoma_corgenes.png",width = 1000,height = 1200,res = 300)
ggplot(data = res, aes(y=expression,x=correlation,color=threshold)) +
  geom_point(alpha=0.5, size=2)+
  geom_text_repel(data = label_data,label = rownames(label_data),color = 'black',max.overlaps = 10,size = 3,show.legend = F)+
  scale_color_manual(values = c("Positive" = "#FF745A","Negative"="#007E99","Background" = "#ABABAB"))+
  ylab('Mean normalized expression')+
  xlab('Pearson cofficient ')+
  theme_classic()+
  theme(text = element_text(size = 16),legend.position = 'none')
dev.off()


### icanet heatmap
library(ComplexHeatmap)
library(circlize)

root_cell = paste('DPT',grep(names(sort(dm$DC1+dm$DC2,decreasing = T)[1]),names(dm$DC1)),sep = '')
EpiExp.m$DPT <- dpt[[root_cell]]
groups = sort(EpiExp.m$DPT,decreasing = F)
top_anno <- HeatmapAnnotation(df = data.frame(DPT = groups,State = EpiExp.m$State[names(groups)]),show_legend = T,col = list(DPT = colorRamp2(c(0,4,8), c("#f2cc8f", "#f4f1de",'#81b29a')),State = bar))
col <- colorRamp2(c(-1,0,1), c("#023e8a","white", "#e63946"), space = "LAB")
DefaultAssay(EpiExp.m) <- 'IcaNet'
EpiExp.m <- ScaleData(EpiExp.m)
mat = EpiExp.m@assays$IcaNet@scale.data %>% as.matrix()
mat = mat[,names(groups)]

png('HUSI/Figures/melanome/Melanoma_icanet.png',width = 2500,height = 1800,res= 300)
Heatmap(mat,
        show_column_names = F,
        show_row_names = T,
        row_title = NULL,
        col = col,
        cluster_rows = T,
        cluster_row_slices = FALSE,
        cluster_columns = F,
        top_annotation = top_anno,
        row_names_gp = gpar(fontsize = 12),
        heatmap_legend_param = list(title = "Scaled\nexpression"))
dev.off()

### melanome trajectory
png('HUSI/Figures/melanome/Melanoma_diffusion.png',width = 1400,height = 1000,res = 300)
ggplot(data = data.frame(DC1 = dm$DC1,DC2 = dm$DC2, State = dm$State)) +
  geom_point(aes(DC1, DC2,colour = State),size = 2)+
  scale_color_manual(values = bar)+
  theme_classic()+
  theme(legend.position = c(0.3,0.7),text = element_text(size =16))+
  ggtitle('Melanoma tumor cells')
dev.off()

### aging markers
DefaultAssay(EpiExp.m)='RNA'
SenMarkers = c("CDKN1A", "SERPINE1")
EpiExp.m[['diffusion']] <- CreateDimReducObject(data.frame(DC1=dm$DC1,DC2=dm$DC2,row.names = colnames(exp)) %>% as.matrix(),key="DC_")
DefaultAssay(EpiExp.m) <- 'RNA'

p1 = FeaturePlot(EpiExp.m,features=SenMarkers[1],reduction='diffusion',order=T,pt.size = 2) + xlim(-0.05,0.05) + ylim(-0.03,0.15)
p2 = FeaturePlot(EpiExp.m,features=SenMarkers[2],reduction='diffusion',order=T,pt.size = 2) + xlim(-0.05,0.05) + ylim(-0.03,0.15)
p3 = VlnPlot(EpiExp.m,features=SenMarkers[1],group.by='State',assay = 'RNA',cols = bar,pt.size = 0) + theme(axis.title.x = element_blank()) + NoLegend()
p4 = VlnPlot(EpiExp.m,features=SenMarkers[2],group.by='State',assay = 'RNA',cols = bar,pt.size = 0) + theme(axis.title.x = element_blank())+ NoLegend()

png('HUSI/Figures/melanome/Melanoma_markers.png',width = 2500,height = 2000,res = 300)
ggarrange(p1,p3,p2,p4,ncol = 2,nrow = 2) 
dev.off()

### microarray gene overlap
colorList = list()
colorList[['Cycling_up']] <- c("#e63946","#2a9d8f")
colorList[['Cycling_down']] <- c("#007E99","#2a9d8f")
colorList[['Transitional_up']] <- c("#e63946","#e9c46a")
colorList[['Transitional_down']] <- c("#007E99","#e9c46a")
colorList[['Senescent_up']] <- c("#e63946","#e76f51")
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
                            text_size = 4,fill_alpha = 0.8)+
                    labs(title = paste("Padj value:",as.character(signif(enrich_reList[[set]],2))))+
                    theme(plot.title = element_text(hjust = 0.5,vjust = 0,size = 15))+
                    scale_y_continuous(limits = c(-1, 1))
    
}
png('HUSI/Figures/melanome/Melanoma_enrich.png',width = 1500,height = 600,res = 160)
ggarrange(pList$Cycling_up,
            pList$Transitional_up,
            pList$Senescent_up,
            pList$Cycling_down,
            pList$Transitional_down,
            pList$Senescent_down,
            ncol = 3,nrow = 2,legend = "none")
dev.off()

### marker enrichment plot
png('HUSI/Figures/melanome/Melanoma_state_enrich.png',width = 2000,height = 1500,res = 300)
dotplot(go,showCategory = 5) + theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5))+
  ggtitle("GO:BP Enrichment of State Marker Genes")+theme(plot.title = element_text(hjust = 1))
dev.off()

### draw heatmap of SKCM CIBERSORT
png('HUSI/Figures/melanome/TCGA_SKCM_state_immue_cor.png',width = 1800,height = 1200,res= 300)
bk <- c(seq(-0.3,0.3,by=0.01))
cor_mat <- cormat_gtex
pheatmap::pheatmap(cor_mat,show_colnames = T,show_rownames = T,
                   cluster_rows = F,cluster_cols = T,fontsize = 12,
                   border_color = "white",color = colorRampPalette(c("blue", "#f5f4f4", "red"))(length(bk)),
                   legend_breaks=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3),breaks=bk)
dev.off()

### survival plot
plotsurv <- function(myfit){
    p <- ggsurvplot(
    myfit,
    risk.table = T,
    pval = TRUE,
    conf.int = F,
    xlim = c(0,4000),
    break.time.by = 1000,
    risk.table.y.text.col = T,
    risk.table.y.text = FALSE)
    return(p)
}

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'Cycling')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ Cycling,data = dat)
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/TCGA_SKCM_survival_Cycling.png',width = 1200,height = 1600,res = 300)
plotsurv(fit)+xlab('Time (days)')
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'Transitional')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ Transitional,data = dat)
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/TCGA_SKCM_survival_Transitional.png',width = 1200,height = 1600,res = 300)
plotsurv(fit)+xlab('Time (days)')
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'Senescent')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ Senescent,data = dat)
png('/home/wangjing/wangj/codebase/HUSI/Figures/melanome/TCGA_SKCM_survival_Senescent.png',width = 1200,height = 1600,res = 300)
plotsurv(fit)+xlab('Time (days)')
dev.off()

### cell chat strength count
png('HUSI/Figures/melanome/melanome_cellchat_scatter.png',width = 1200,height = 1000,res = 300)
netAnalysis_signalingRole_scatter_log(cellchat,color.use = bar,dot.alpha = 0)
dev.off()

### cell chat strength network
groupSize <- table(cellchat@idents) %>% as.numeric()
png('HUSI/Figures/melanome/melanome_cellchat_network.png',width = 1600,height = 800,res = 300)
par(mfrow = c(1,2), xpd=TRUE,mar = c(0, 1, 1, 0.5))
netVisual_circle(cellchat@net$weight, 
                 sources.use = c("Cycling","Transitional","Senescent"),
                 color.use =  bar[colnames(cellchat@net$count)],
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge = F, 
                 title.name = "Outgoing weights/strength")+
netVisual_circle(cellchat@net$weight, 
                   targets.use = c("Cycling","Transitional","Senescent"),
                   color.use =  bar[colnames(cellchat@net$count)],
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   label.edge = F, 
                   title.name = "Incoming weights/strength")
dev.off()

### different pattern
png('HUSI/Figures/melanome/melanome_cellchat_pattern_2.png',width = 2000,height = 1000,res = 300)
netAnalysis_river(cellchat, pattern = "outgoing",cutoff = 0.5,
                  color.use = bar[c('Cycling','Transitional','Senescent')],
                  color.use.pattern = c('Pattern 1'='#e63946','Pattern 6'='#457b9d'),
                  sources.use = c('Cycling','Transitional','Senescent'),
                  signaling = as.character(pattern_use$Signaling))
dev.off()


netVisual_bubble_my(cellchat,signaling = 'AGRN',
                   sources.use = c('Cycling','Transitional','Senescent'), 
                   remove.isolate = F,
                   targets.use = c('T cell','NK cell','Macro cell','CAF cell','Endo cell'),thresh=0.01)


### pathway strength bubble
# png('HUSI/Figures/melanome/melanome_cellchat_pathway.png',width = 1500,height = 1200,res = 300)
pdf('HUSI/Figures/melanome/melanome_cellchat_pathway.pdf',width = 6,height = 5)
netVisual_bubble_my(cellchat,signaling = pathways,
                    targets.use = c('Cycling','Transitional','Senescent'), 
                    remove.isolate = F,
                    sources.use = c('T cell','NK cell','Macro cell','CAF cell','Endo cell'),thresh=0.01)
dev.off()

### pathway strength cord
png('HUSI/Figures/melanome/melanome_cellchat_chord.png',width = 2000,height = 1000,res = 300)
par(mfrow = c(1,2), xpd=TRUE,mar = c(0.2, 0.2, 0.2,0.2))
for(p in c('BMP','TGFb')){
    netVisual_aggregate(cellchat, signaling = p,layout = "chord",remove.isolate = F,
                        color.use = bar,sources.use = c('T cell','NK cell','Macro cell','CAF cell'),
                        targets.use = c('Cycling','Transitional','Senescent'))
}
dev.off()

pdf('HUSI/Figures/melanome/melanome_cellchat_chord_2.pdf',width = 6,height = 6)
par(mfrow = c(2,2), xpd=TRUE,mar = c(0.2, 0.2, 0.2,0.2))
for(p in c("CCL","CSPG4","LEP","CD6")){
    netVisual_aggregate(cellchat, signaling = p,layout = "chord",remove.isolate = F,
                        color.use = bar,sources.use = c('T cell','NK cell','Macro cell','CAF cell'),
                        targets.use = c('Cycling','Transitional','Senescent'))
}
dev.off()


### pathway gene expression
png('HUSI/Figures/melanome/melanome_cellchat_expression.png',width = 900,height = 1000,res = 300)
plotGeneExpression(cellchat, features=c('BMPR1B','BMPR2','TGFBR1','TGFBR2'),idents = c('Cycling','Transitional','Senescent'),color.use = bar)
dev.off()

png('HUSI/Figures/melanome/melanome_cellchat_expression_2.png',width = 1000,height = 1000,res = 400)
plotGeneExpression(cellchat, features=c('LEPR','ITGA2','ITGB1'),idents = c('Cycling','Transitional','Senescent'),color.use = bar)
dev.off()

png('HUSI/Figures/melanome/melanome_cellchat_expression_CD6.png',width = 1000,height = 1000,res = 400)
plotGeneExpression(cellchat, features=c('ALCAM'),idents = c('Cycling','Transitional','Senescent'),color.use = bar)
dev.off()

png('HUSI/Figures/melanome/melanome_cellchat_expression_CCL.png',width = 1000,height = 1000,res = 400)
plotGeneExpression(cellchat, features=c('CCR10'),idents = c('Cycling','Transitional','Senescent'),color.use = bar)
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



cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'BMPR2')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ BMPR2,data = dat)
png('HUSI/Figures/melanome/TCGA_SKCM_survival_BMPR2.png',width = 1400,height = 1500,res = 300)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'BMPR1B')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ BMPR1B,data = dat)
png('HUSI/Figures/melanome/TCGA_SKCM_survival_BMPR1B.png',width = 1400,height = 1500,res = 300)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'TGFBR1')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ TGFBR1,data = dat)
png('HUSI/Figures/melanome/TCGA_SKCM_survival_TGFBR1.png',width = 1400,height = 1500,res = 300)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'TGFBR2')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ TGFBR2,data = dat)
png('HUSI/Figures/melanome/TCGA_SKCM_survival_TGFBR2.png',width = 1400,height = 1500,res = 300)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'CCR10')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ CCR10,data = dat)
png('HUSI/Figures/melanome/TCGA_SKCM_survival_CCR10.png',width = 1400,height = 1500,res = 300)
plotsurv(fit)
dev.off()

cut <- surv_cutpoint(data,time = "OS.time",event = "OS",variables = 'ALCAM')
dat <- surv_categorize(cut)
fit <- survfit(Surv(OS.time, OS) ~ ALCAM,data = dat)
png('HUSI/Figures/melanome/TCGA_SKCM_survival_ALCAM.png',width = 1400,height = 1500,res = 300)
plotsurv(fit)
dev.off()



### re-cut group
# c('BMPR2','BMPR1B','TGFBR1','TGFBR2')
gene = 'BMPR2'

library(cutoff)
cut_2 <- logrank(data = data,
                 time = "OS.time",  
                 y = "OS",         
                 x = gene,  
                 cut.numb = 1,      
                 n.per = 0.1,       
                 y.per = 0.4,       
                 p.cut = 0.1,       
                 round = 5) 

cut_2[order(cut_2$pvalue, decreasing = F), ]

# data$Group = ifelse(data[gene] < cut_2[order(cut_2$pvalue, decreasing = F), ][1, 1],"Low",
#                     ifelse(data[gene] <= cut_2[order(cut_2$pvalue, decreasing = F), ][1, 2], 'Middle','high'))

data$Group = ifelse(data[gene] < cut_2[order(cut_2$pvalue, decreasing = F), ][1, 1],"Low",'high')

surv_obj <- Surv(time = data$OS.time, event = data$OS)
surv_fit <- survfit(surv_obj ~ data$Group)

ggsurvplot(surv_fit, data = data, surv.median.line = "hv",
           pval = TRUE,                   
           xlab = "Time (days)", ylab = "Survival Probability", 
           legend.title = "",            
           break.x.by = 1000, conf.int = TRUE,           
           color = "strata")             


p <- ggsurvplot(
  myfit,
  risk.table = F,
  pval = TRUE,
  conf.int = TRUE,
  xlim = c(0,4000),
  break.time.by = 1000,
  risk.table.y.text.col = T,
  risk.table.y.text = FALSE)











