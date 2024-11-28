##### This scripts is for Figure 1,S1
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
library(patchwork)
library(gg.gap)

setwd('hUSI/')
################################ preprocess train set #############################
load('Data/raw_full.rdata')
X = exprData
### senescent marker expression in train set
genes = c('CDKN1A','CDKN2B','SERPINE1' )
df_plot = data.frame(t(X[genes,]),label = meta$label) %>%
  melt(value.name = 'expression',variable.name = 'marker',id.var = 'label') 
df_plot$expression=log(df_plot$expression+1)

pwc <- df_plot %>% group_by(marker) %>% 
  rstatix::pairwise_t_test(expression ~ label, paired = F,p.adjust.method = "bonferroni")
pwc
pwc <- pwc %>% rstatix::add_xy_position(x = "marker")


pdf('Results/Model/valid_markers.pdf',width = 3,height = 4)
df_plot %>%
  ggplot(aes(x = marker, y = expression)) + 
  geom_boxplot(aes(color = label),outlier.shape = NA) + 
  stat_pvalue_manual(pwc)+
  theme_classic()+
  theme(text = element_text(size = 14),axis.text.x = element_text(size=14,angle = 45,hjust=1),legend.position = 'top',axis.title.x = element_blank())+
  scale_color_manual(values = c('senescent' = '#e63946','non-senescent'='#959596'))+
  ylab('Normalized Expression')
dev.off()

############################### merge cell type ################################
####### merge cell type
df_plot = table(meta$Cell_type,meta$Cell_type_id) %>% as.array() %>% data.frame()
colnames(df_plot) <- c('Cell_type','Cell_type_id','Percent')

pdf('Results/Model/bulk_barplot_celltype.pdf',width =9,height = 6)
ggplot(df_plot,aes(x=reorder(Cell_type,-Percent),y=Percent,fill=Cell_type_id))+
  geom_bar(stat='identity',position = 'fill')+
  scale_fill_manual(values = bar)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  xlab('Merge cell type')
dev.off()

####### plot DEGs
DEGs_celltypes = DEGs_celltypes$`non-senescent`
gene_num=85
MarkerSets <- lapply(DEGs_celltypes, function(x) {
  x<-filter(x,log2FoldChange>1,padj<0.05);
  x<-x[order(x$log2FoldChange,decreasing = T),]; 
  return(na.omit(rownames(x)[0:gene_num]))})

gene=unique(unlist(MarkerSets))
mat = sapply(DEGs_celltypes, function(c){c[gene,'log2FoldChange']})
rownames(mat) = gene

pdf('Results/Model/merge_celltype_DEG2.pdf',width =10,height = 6)
col <- colorRamp2(c(-4,0,4), c("#023e8a","white", "#e63946"), space = "LAB")
ComplexHeatmap::Heatmap(t(mat),cluster_rows = F,cluster_columns = F,show_column_names = F,
                        col=col,heatmap_legend_param = list(title = "LFC"),)
dev.off()

####### plot remove marker
res0=0.6
pdf('Results/Model/ARI_celltype_markers.pdf',width =4,height = 3)
ggplot(ari_res %>% filter(resolution==res0),aes(x=gene_num,y=ARI))+
  geom_point()+theme_bw()+
  geom_smooth(method = 'loess')+
  geom_hline(yintercept = (ari_res %>% filter(gene_num==0,resolution==res0))$ARI,color='red',linetype='dashed')+
  geom_vline(xintercept = 85,color='red',linetype='dashed')+
  labs(title=paste0("Resolution=",res0))+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12))+theme_classic()
dev.off()

########################## define l2 and quantification ########################
df_plot = lapply(auc_models,function(x) {lapply(x,function(y) mean(unlist(y))) })
df_plot = rbindlist(df_plot) %>% data.frame()
rownames(df_plot) <- names(auc_models)
df_plot$quntification = rownames(df_plot)
df_plot = melt(df_plot,variable.name = c('l2'),value.name = 'AUC')
df_plot$l2 = gsub('l2_','',df_plot$l2) 
df_plot$l2 = factor(df_plot$l2 ,levels = c("0.001","0.01","0.1","1","5","10"))
head(df_plot)

pdf('/Results/Model/AUC_l2.pdf',width = 5,height =4)#6,4.5
p <- ggplot(df_plot,aes(fill=quntification,y=AUC,x=l2 ))+
  geom_bar(stat = 'identity',position = 'dodge')+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))+
  scale_fill_manual(values = c("#91D1C2","#FFBF80","#B2E5FF","#EE0000FF"))+
  ylab('mean AUC')

gg.gap(plot = p,
         segments = c(0.1, 0.9),
         tick_width = 0.05,
         rel_heights = c(0.1, 0, 0.1),
         ylim = c(0,0.98),margin = margin(0,2.2,0,0.5,unit = 'cm'))

add.legend(plot = p,
           margin = c(top=100,right=100,bottom=100,left=320,unit = 'cm'))

dev.off()

pdf('Results/Model/AUC_l2_2.pdf',width = 5,height = 4)
p <- ggplot(df_plot,aes(x=quntification,y=AUC,fill=l2 ))+
  geom_bar(stat = 'identity',position = 'dodge')+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))+
  scale_fill_manual(values =  colorRampPalette(brewer.pal(9, "Reds"))(7)[-1])+
  ylab('mean AUC')

gg.gap::gg.gap(plot = p,
       segments = c(0.1, 0.9),
       tick_width = 0.05,
       rel_heights = c(0.1, 0, 0.1),
       ylim = c(0,0.98),margin = margin(0,2.2,0,0.5,unit = 'cm'))

gg.gap::add.legend(plot = p,
           margin = c(top=100,right=100,bottom=100,left=320,unit = 'cm'))

dev.off()

################### Compare all quantification methods in final model ############################
df_plot = lapply(auc_models,function(x) {lapply(x,function(y) mean(unlist(y))) })
df_plot = rbindlist(df_plot) %>% data.frame()
rownames(df_plot) <- names(auc_models)
df_plot$quntification = rownames(df_plot)
df_plot = melt(df_plot,variable.name = c('l2'),value.name = 'AUC')
df_plot$l2 = gsub('l2_','',df_plot$l2) 

df_plot = df_plot[order(df_plot$AUC,decreasing = T),]
df_plot$rank = 1:nrow(df_plot)
df_plot$quntification = gsub('spearman','spearman\n(hUSI)',df_plot$quntification)

pdf("Results/Model/AUC_k-folds_OOD_drop_l2=1.pdf",width = 5,height = 3.5)
ggplot(df_plot[1:50,], aes(x=-rank, y=AUC)) +
geom_point(aes(color = rank <= 1),size=2) +  
geom_text_repel(aes(label=ifelse(rank <= 10, quntification, "")), 
                vjust=-0.5, hjust=0.8,max.overlaps = 20) +  
scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey80")) + 
theme_bw() %+replace% 
theme(text = element_text(size=14),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 20,r = 20,b = 10,l = 10))+
ylab('mean AUC') +
xlab('Quntification')
dev.off()

################################## simulation ##################################
###### batch effect
df_plot=do.call(rbind,BEhUSIs)
df_plot$quantification = unlist(lapply(strsplit(rownames(df_plot),'\\.'),'[',1)) %>% 
  factor(levels = c('dot','logistic','pearson','spearman'))
df_plot$study = RNAmeta$study_accession[match(colnames(BEexp),RNAmeta$sample_title)]
labels = as_labeller(function(x) paste('AUC: ',as.character(round(BEAUC[x],3)),sep=''))

pdf('Results/Model/simulation_BE2.pdf',width = 4.5,height = 4)
set.seed(223)
ggplot(df_plot %>% filter(quantification=='spearman')) +
  geom_jitter(aes(y=hUSI,x=quantification,color=label,shape=study))+
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)) +
  scale_color_manual(values = c('senescent' = '#e63946','other'='#959596'))+
  facet_grid(~quantification,scales = 'free_x',labeller = labels)+
  theme_bw()+
  theme(strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        text = element_text(size=14))
dev.off()

###### sparsity
df_plot = data.frame(AUC=BEAUC_sparse,
                     quantifcation=unlist(lapply(strsplit(names(BEAUC_sparse),'\\.'),'[',1)),
                     depth=unlist(lapply(strsplit(names(BEAUC_sparse),'_'),'[',2)))
df_plot = aggregate(df_plot$AUC,by=list(df_plot$quantifcation,df_plot$depth),mean)
names(df_plot) <- c('quantifcation','depth','mean_auc')
df_plot$quantifcation = factor(df_plot$quantifcation,levels = c('spearman'))
df_plot$depth = factor(df_plot$depth,levels = c('0.8','0.6','0.4','0.2'))

pdf('Results/Model/simulation_sparsity.pdf',width = 2,height = 3)
ggplot(df_plot %>% filter(quantifcation=="spearman")) +
  geom_bar(aes(y=mean_auc,x=depth,fill=depth),stat = 'identity',position = 'dodge')+
  scale_fill_npg(alpha = 0.6)+
  ylab('mean AUC')+theme_classic()+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size=12),
        axis.title  = element_text(size=14),legend.position = "none")
dev.off()

###### out liner
df_plot = data.frame(AUC=BEAUC_outliner,
                     quantifcation=unlist(lapply(strsplit(names(BEAUC_outliner),'\\.'),'[',1)),
                     percent =unlist(lapply(strsplit(names(BEAUC_outliner),'_'),'[',2)))

df_plot = aggregate(df_plot$AUC,by=list(df_plot$quantifcation,df_plot$percent),mean)
names(df_plot) <- c('quantifcation','percent','mean_auc')
df_plot$quantifcation = factor(df_plot$quantifcation,levels = c('spearman'))
df_plot$percent = factor(df_plot$percent ,levels = c("0.01","0.05","0.1","0.2"))

pdf('Results/Model/simulation_outliner.pdf',width = 2,height = 3)
ggplot(df_plot %>% filter(quantifcation=="spearman")) +
  geom_bar(aes(y=mean_auc,x=percent,fill=percent),stat = 'identity',position = 'dodge')+
  # scale_fill_manual(values = c("#91D1C2","#FFBF80","#B2E5FF","#EE0000FF"))+
  scale_fill_npg(alpha = 0.6)+
  ylab('mean AUC')+theme_classic()+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title  = element_text(size=14),legend.position = "none")
dev.off()

########################### AUC_th,FN,FP,F1 score ########################
### FN,FP,F1 score
pdf('Results/Model/Metrics.pdf',width = 2.5,height = 3)
data.frame(mean(metrics$false_negative),mean(metrics$false_positive),mean(metrics$F1_score),AUC_th) %>% 
  set_names(c('FN','FP','F1','AUC_binary')) %>% melt %>% mutate(variable=factor(variable,levels = c('FN','FP','F1','AUC_binary'))) %>%
ggplot(aes(x=variable,y=value,fill=variable)) +
  geom_bar(stat = 'identity')+
  geom_text(aes(label=sprintf("%.2f", value)), vjust=1.5, color="black") +
  theme_bw()+
  theme(legend.position = 'none',axis.title.x = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),text = element_text(size=12))+
  scale_fill_manual(values = c('#51C3CC','#72D9FF','#FFAD72','#5fc48c'))#+labs(y="Score")
dev.off()

############################## two-class model AUC #############################
### test1
df_plot = lapply(auc_models,function(x) {lapply(x,function(y)mean(y))})
df_plot = lapply(df_plot,function(x) {data.frame(celltype=names(x),auc=unlist(x))})
df_plot = lapply(names(df_plot),function(x) {df_plot[[x]]$model = x;return(df_plot[[x]])})
df_plot = do.call(rbind,df_plot)
df_plot$model = factor(df_plot$model,levels = c('RF','EN','SVM','OCLR'),ordered = T)
head(df_plot)

pdf('Results/Model/OOD_two-class_celltype_test1.pdf',width =10,height =4.5)
ggplot(df_plot %>% mutate(celltype=gsub("t_","",celltype)),aes(x=reorder(celltype,auc),y=auc,fill=model))+
  geom_bar(stat = 'identity',position = 'dodge')+
  theme_bw()+
  theme(axis.text = element_text(size=12),
    axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_fill_manual(values = c("#7EC3E5","#C3E57E","#E5B17E","#EE0000FF"))+
  xlab('Cell type')+
  ylab("mean AUC")
dev.off()


df_plot = lapply(auc_models,function(x) unlist(x) %>% data.frame() %>% set_names('auc'))
df_plot = lapply(names(df_plot),function(x) {df_plot[[x]]$model = x;return(df_plot[[x]])})
df_plot = do.call(rbind,df_plot)
df_plot = aggregate(auc~model,data = df_plot,mean)
df_plot$model = factor(df_plot$model,levels = c('RF','EN','SVM','OCLR'),ordered = T)
head(df_plot)

pdf('Results/Model/OOD_two-class_mean_test1.pdf',width = 3,height = 4)
p<-ggplot(df_plot,aes(x=reorder(model,auc),y=auc,fill=model))+
  geom_bar(stat = 'identity',position = 'dodge')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  scale_fill_manual(values = c("#7EC3E5","#C3E57E","#E5B17E","#EE0000FF"))+
  xlab('model')+
  ylab("mean AUC")+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14))
gg.gap(plot = p,
       segments = c(0.1, 0.9),
       tick_width = 0.05,
       rel_heights = c(0.1, 0, 0.1),
       ylim = c(0,0.98),margin = margin(0,0.5,0,0.5,unit = 'cm'))
dev.off()

### test2
df_plot = lapply(auc_models,function(x) {lapply(x,function(y)mean(y))})
df_plot = lapply(df_plot,function(x) {data.frame(celltype=names(x),auc=unlist(x))})
df_plot = lapply(names(df_plot),function(x) {df_plot[[x]]$model = x;return(df_plot[[x]])})
df_plot = do.call(rbind,df_plot)
df_plot$model = factor(df_plot$model,levels = c('RF','EN','SVM','OCLR'),ordered = T)
df_plot$Train_Celltypes = unlist(lapply(strsplit(df_plot$celltype,'_'),'[',1))
df_plot = aggregate(auc~model+Train_Celltypes,data = df_plot,mean)
df_plot$Train_Celltypes=gsub("k","",df_plot$Train_Celltypes) %>% paste0(.,"ct")
df_plot$Train_Celltypes = factor(df_plot$Train_Celltypes,levels = c('30ct','20ct','10ct','5ct'))
head(df_plot)

pdf('Results/Model/two-class_mean_test2.pdf',width = 5,height = 4)
p<-ggplot(df_plot,aes(x=Train_Celltypes,y=auc,fill=model))+
  geom_bar(stat = 'identity',position = 'dodge')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0,hjust = 1))+
  scale_fill_manual(values = c("#7EC3E5","#C3E57E","#E5B17E","#EE0000FF"))+
  xlab('model')+
  ylab("mean AUC")+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14))
gg.gap(plot = p,
       segments = c(0.1, 0.85),
       tick_width = 0.05,
       rel_heights = c(0.1, 0, 0.1),
       ylim = c(0,0.98),margin = margin(0,1.5,0,0.5,unit = 'cm'))
add.legend(plot = p,
           margin = c(top=100,right=100,bottom=100,left=390,unit = 'cm'))
dev.off()

############################## enrichment cellmarker ###########################
df_plot$set_type = ifelse(df_plot$Description %in% c('SenUp','SenMayo'),'Senescence','Cell type')
df_plot$Enrichment = ifelse(df_plot$p.adjust<=0.05,'Significant','non-significant')

df_plot2=df_plot %>% filter(Description!="SenMayo")## has 115 pathways
pdf('Results/Model/CellMarkerEnrich2.pdf',width =4.5,height =3)
ggplot(df_plot2,aes(x = NES,y=-log10(p.adjust),color = Enrichment,shape=set_type))+
  geom_point(size=2,alpha=0.8)+
  # scale_color_distiller(palette = 'Reds',direction = 1)+
  scale_color_manual(values = c('grey60','#e63946'))+
  scale_shape_manual(values = c(19, 17))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  xlab('NES')+
  ylab('-log10(padj)')+
  geom_vline(xintercept = 0,linetype="dashed",color="grey80")
dev.off()

################################### GSEA #######################################
library(enrichplot)
library(gggsea)

d='Hallmarker'
df = filter(re,Database==d)
df[,c('NES','p.adjust')]

### hallmarker pathway
ids = c('TNFA_SIGNALING_VIA_NFKB',
        'INFLAMMATORY_RESPONSE',
        'P53_PATHWAY',
        'KRAS_SIGNALING_UP',
        'INTERFERON_GAMMA_RESPONSE',
        'INTERFERON_ALPHA_RESPONSE')
pdf('Results/Model/enrichment_up2.pdf',width =7,height =5)
enrich_plot(gsea=fgseaList[[d]],data=df,pathway_id=ids,direction = 'up',save = F,database=d)
dev.off()

ids = c('E2F_TARGETS', 
        'G2M_CHECKPOINT', 
        'MYC_TARGETS_V1', 
        'MITOTIC_SPINDLE')
pdf('Results/Model/enrichment_down2.pdf',width =7,height =5)
enrich_plot(gsea=fgseaList[[d]],data=df,pathway_id=ids,direction = 'down',save = F,database=d)
dev.off()


### other databases
top_re = re %>% group_by(Database) %>% top_n(5,NES)
top_re = filter(top_re,NES>0)
top_re$Description = paste(top_re$Database,":",top_re$Description,sep = '')

pdf('Results/Model/DatabaseEnrich_up.pdf',width = 7,height = 6)
ggplot(top_re,aes(y = NES,x=reorder(Description,NES)))+
  geom_bar(stat='identity',fill="#E57E7E")+
  theme_classic()+
  theme(legend.position = 'none',
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_text(color="black",size=10),
        axis.line.x = element_line(color='black'),
        axis.line = element_line(colour = "black",linewidth = rel(1)),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  coord_flip()+
  geom_text(data = top_re,aes(label = reorder(Description, NES),y= 0.01),size = 3,hjust = 0)+
  xlab('Description')+
  ylab('NES')
dev.off()

top_re = re %>% group_by(Database) %>% top_n(5,-NES)
top_re = filter(top_re,NES<0)
top_re$Description = paste(top_re$Database,":",top_re$Description,sep = '')

pdf('Results/Model/DatabaseEnrich_down.pdf',width = 7,height = 6)
ggplot(top_re,aes(y = -NES,x=reorder(Description,-NES)))+
  geom_bar(stat='identity',fill="#7EC3E5")+
  theme_classic()+
  theme(legend.position = 'none',
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_text(color="black",size=10),
        axis.line.x = element_line(color='black'),
        axis.line = element_line(colour = "black",linewidth = rel(1)),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  coord_flip()+
  geom_text(data = top_re,aes(label = reorder(Description, NES),y= 0.01),size = 3,hjust = 0)+
  xlab('Description')+
  ylab('-NES')
dev.off()



### model weight of classical cycling markers
load('SenOCLR_l2=1_drop.rdata')
SenSetlist <<- loadSetdata()
SenSet <- SenSetlist$SenMayo
SenSet

SenSet = SenSet[SenSet %in% names(SenOCLR$w)]

weights = sort(SenOCLR$w,decreasing = T) 

genes = names(weights)[which(names(weights) %in% SenSet)]
rank = which(names(weights) %in% SenSet)
genes = genes[rank/length(weights) < 0.1]
genes
weights[which(names(weights) %in% c('SERPINE1','CDKN1A','CDKN2A','MKI67',"LMNB1"))]

pdf('Results/Model/weight_density2.pdf',width=6.8,height=4)
ggplot(data.frame(weight = weights))+
  geom_density(aes(x=weight))+
  geom_vline(xintercept = c(weights[floor(length(weights)*0.1)],weights[floor(length(weights)*0.9)]),color='#2b2a2a',linetype="dashed")+
  geom_vline(xintercept = c(weights[which(names(weights) %in% c('CDKN1A'))]),color='#ac0f0f')+
  geom_vline(xintercept = c(weights[which(names(weights) %in% c('CDKN2A'))]),color='red',size=0.1)+
  geom_vline(xintercept = c(weights[which(names(weights) %in% c('SERPINE1'))]),color='black',size=0.1)+
  geom_vline(xintercept = c(weights[which(names(weights) %in% c('MKI67',"LMNB1"))]),color='#202269')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab('Weight')+
  ylab('Density')
dev.off()





