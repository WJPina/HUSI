### This script is for Fig2 and FigS3 in the paper

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
library(magrittr)

################################# plot auc microarray #####################################
### only keep common methods 
method_used <- Reduce(intersect,lapply(AUClist, names))
df_plot_list = lapply(AUClist, function(x) {x = x[method_used]}) 

df_plot = do.call(cbind, df_plot_list)
df_plot = melt(df_plot)
colnames(df_plot) <- c('method','dataset','AUC')

df_plot$method_type <- ifelse(df_plot$method %in% Methods$TraditionScore,'TraditionScore',
                               ifelse(df_plot$method %in% Methods$SenSet,'SenSet',
                                      ifelse(df_plot$method %in% Methods$SenMarker,'SenMarker',
                                             ifelse(df_plot$method %in% Methods$MachineScore,'MachineScore','hUSI'))))

df_plot$method_type <- factor(df_plot$method_type,levels = c('hUSI','SenMarker','SenSet','TraditionScore','MachineScore'),ordered = T)

df_plot_mean = aggregate(df_plot$AUC,by=list(df_plot$method),mean)
df_plot_mean$sd = aggregate(df_plot$AUC,by=list(df_plot$method),sd)$x
colnames(df_plot_mean) = c('method','mean','sd')
df_plot_mean$method_type = df_plot$method_type[match(df_plot_mean$method,df_plot$method)]


pdf("Results/Validation/compare_auc_array2.pdf",width =8,height = 4)
ggplot()+
  geom_bar(data=df_plot_mean,aes(x=reorder(method,mean), y=mean,fill = method_type),stat="identity",position="dodge")+
  geom_errorbar(data=df_plot_mean,aes(x=reorder(method,mean),y=mean,ymax=mean+sd,ymin=ifelse(mean-sd <0,0,mean-sd)),
                position=position_dodge(0.9),width=0.4)+
  geom_jitter(data=df_plot,aes(x = method,y=AUC,color=dataset),width = 0.3,size=1.1,height = 0) +
  theme_classic()+
  ylab('AUC')+
  xlab('Senescence scoring method')+
  # scale_fill_manual(values = c('#e63946','#94d2bd','#00afb9','#fed9b7','#CFDCEF'))+
  scale_fill_manual(values = c('#D2352C','#71ABB6','#BAAFD1','#FFD699','#CFDCEF'))+
  # scale_color_manual(values = c('#277da1','#f9c74f','#7209b7','#68001D','#003560','#5C5C5C'))+
  scale_color_manual(values = c('#277da1','#fc8002','#8481BA','pink','#ADDB88','#F14454'))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1),
        text=element_text(size=14))+scale_y_continuous(breaks = seq(0, 1, by = 0.25)) 
dev.off()


################### visualization single cell data sets ######################## add significance test !!!
d="Teo2019"
obj = scDataSets[[d]]
obj <- FindVariableFeatures(obj) %>% ScaleData() %>% RunPCA() %>% RunTSNE(dims=1:15)
obj$hUSI = datlist[[d]][colnames(obj),]$hUSI
metadata = fread("Data/GSE115301_Growing_Sen_10x_metadata.txt.gz", header = T) %>% column_to_rownames("V1")
obj$group = metadata[colnames(obj),'Condition2']
obj$group = ifelse(obj$group == 'RIS','Primary Senecent',ifelse(obj$group == 'GFP','Sencondary Senecent','Growing'))
obj$group = factor(obj$group,levels = c('Growing','Sencondary Senecent','Primary Senecent'))

pdf(paste('Results/Validation/benchmark_',d,'_tsne.pdf',sep=''),width = 13,height = 3.5)
p1 <- DimPlot(obj,group.by = 'group',reduction = 'tsne',cols = c('#91D1C2','#FFD699','#E64B35'),label = T,label.box = T,repel = T)+ggtitle(d)+NoLegend()
p2 <- FeaturePlot(obj,features = 'hUSI',)+scale_color_gradient2(low ='#3AB370' ,mid = "#EAE7CC",high = "#FD1593",midpoint = 0.5)
p3 <- VlnPlot(obj,features = 'hUSI',group.by = 'group',cols = c('#91D1C2','#FFD699','#E64B35'),pt.size = 0)+
  geom_boxplot(width=0.2)+theme(axis.title = element_blank(),axis.text.x = element_blank())+
    stat_compare_means(aes(label = ..p.signif..),label.x = 1.5, label.y = 1.03)
p1 + p2 + p3
dev.off()

df_plot = obj@meta.data[,c('hUSI','group')]
pwc <- df_plot %>% 
  rstatix::pairwise_t_test(hUSI ~ group, paired = F,p.adjust.method = "bonferroni")
pwc

d="Teo2019O2"
obj = scDataSets[[d]]
obj <- FindVariableFeatures(obj) %>% ScaleData() %>% RunPCA() %>% RunTSNE(dims=1:15)
obj$hUSI = datlist[[d]][colnames(obj),]$hUSI

obj$group = factor(ifelse(obj$Condition=='senescent','O2-induced Senescence','proliferation'),levels = c('proliferation','O2-induced Senescence'))

pdf(paste('Results/Validation/benchmark_',d,'_tsne.pdf',sep=''),width = 13,height = 3.5)
p1= DimPlot(obj,group.by = 'group',reduction = 'tsne',cols = c('#91D1C2','#E64B35'),label = T,label.box = T,repel = T)+ggtitle(d)+NoLegend()
p2= FeaturePlot(obj,features = 'hUSI',)+scale_color_gradient2(low ='#3AB370' ,mid = "#EAE7CC",high = "#FD1593",midpoint = 0.5)
p3= VlnPlot(obj,features = 'hUSI',group.by = 'group',cols = c('#91D1C2','#E64B35'),pt.size = 0)+
     geom_boxplot(width=0.2)+theme(axis.title = element_blank(),axis.text.x = element_blank())
p1+p2+p3
dev.off()


d="Tang2019"
obj = scDataSets[[d]]
obj <- FindVariableFeatures(obj) %>% ScaleData() %>% RunPCA() %>% RunTSNE(dims=1:15)
obj$hUSI = datlist[[d]][colnames(obj),]$hUSI
obj$group = factor(obj$group,levels = c('LowPDCtrl','HighPDCtrl','senescence','LowPD50Gy'))

pdf(paste('Results/Validation/benchmark_',d,'_tsne.pdf',sep=''),width = 13,height = 3.5)
p1= DimPlot(obj,group.by = 'group',reduction = 'tsne',cols = c('#00CCCC','#91D1C2','#FFD699','#E64B35'),label = T,label.box = T,repel = T)+ggtitle(d)+NoLegend()
p2= FeaturePlot(obj,features = 'hUSI')+scale_color_gradient2(low ='#3AB370' ,mid = "#EAE7CC",high = "#FD1593",midpoint = 0.5)
p3= VlnPlot(obj,features = 'hUSI',group.by = 'group',cols = c('#00CCCC','#91D1C2','#FFD699','#E64B35'),pt.size = 0)+
     geom_boxplot(width=0.2)+theme(axis.title = element_blank(),axis.text.x = element_blank())
p1+p2+p3
dev.off()


d="Aarts2017"
obj = scDataSets[[d]]
obj <- FindVariableFeatures(obj) %>% ScaleData() %>% RunPCA() %>% RunTSNE(dims=1:15)
obj$hUSI = datlist[[d]][colnames(obj),]$hUSI
obj$group = factor(ifelse(obj$Condition=='senescent','OSKM-induced Senescence','proliferation'),levels = c('proliferation','OSKM-induced Senescence'))

pdf(paste('Results/Validation/benchmark_',d,'_tsne.pdf',sep=''),width = 13,height = 3.5)
p1= DimPlot(obj,group.by = 'group',reduction = 'tsne',cols = c('#91D1C2','#E64B35'),label = T,label.box = T,repel = T)+ggtitle(d)+NoLegend()
p2= FeaturePlot(obj,features = 'hUSI')+scale_color_gradient2(low ='#3AB370' ,mid = "#EAE7CC",high = "#FD1593",midpoint = 0.5)
p3= VlnPlot(obj,features = 'hUSI',group.by = 'group',cols = c('#91D1C2','#E64B35'),pt.size = 0)+
     geom_boxplot(width=0.2)+theme(axis.title = element_blank(),axis.text.x = element_blank())
p1+p2+p3
dev.off()


d="Zirkel2018"
obj = scDataSets[[d]]
obj <- FindVariableFeatures(obj) %>% ScaleData() %>% RunPCA() %>% RunTSNE(dims=1:15)
obj$hUSI = datlist[[d]][colnames(obj),]$hUSI
obj$group = factor(ifelse(obj$Condition=='senescent','Replicative Senescence','Growing'),levels = c('Growing','Replicative Senescence'))

pdf(paste('Results/Validation/benchmark_',d,'_tsne.pdf',sep=''),width = 13,height = 3.5)
p1= DimPlot(obj,group.by = 'group',reduction = 'tsne',cols = c('#91D1C2','#E64B35'),label = T,label.box = T,repel = T)+ggtitle(d)+NoLegend()
p2= FeaturePlot(obj,features = 'hUSI')+scale_color_gradient2(low ='#3AB370' ,mid = "#EAE7CC",high = "#FD1593",midpoint = 0.5)
p3= VlnPlot(obj,features = 'hUSI',group.by = 'group',cols = c('#91D1C2','#E64B35'),pt.size = 0)+
     geom_boxplot(width=0.2)+theme(axis.title = element_blank(),axis.text.x = element_blank())
p1+p2+p3
dev.off()


d="Noah2023"
obj = scDataSets[[d]]
DefaultAssay(obj) <- "integrated"
obj <- ScaleData(obj, verbose = FALSE) %>%
            RunPCA(npcs = 30, verbose = FALSE) %>%
            RunTSNE(dims = 1:30)

obj$hUSI = datlist[[d]][colnames(obj),]$hUSI
obj$group = factor(ifelse(obj$orig.ident=='CTRL_2','Growing',
                          ifelse(obj$orig.ident%in%c('ETO_1','ETO_2'),'CIS',
                                 ifelse(obj$orig.ident%in%c('IR_1','IR_2'),'IRIS','RS'))),levels = c('Growing','RS','CIS','IRIS'))

pdf(paste('Results/Validation/benchmark_',d,'_tsne.pdf',sep=''),width = 13,height = 3.5)
p1=DimPlot(obj,group.by = 'group',reduction = 'tsne',cols = c('#00CCCC','#FFD699','#F76D5E','#E64B35'),label = T,label.box = T,repel = T)+ggtitle(d)+NoLegend()
p2=FeaturePlot(obj,features = 'hUSI')+scale_color_gradient2(low ='#3AB370' ,mid = "#EAE7CC",high = "#FD1593",midpoint = 0.5)
p3=VlnPlot(obj,features = 'hUSI',group.by = 'group',cols = c('#00CCCC','#FFD699','#F76D5E','#E64B35'),pt.size = 0)+
     geom_boxplot(width=0.2)+theme(axis.title = element_blank(),axis.text.x = element_blank())
p1+p2+p3
dev.off()

################################# plot auc single cell ##################################### 
### only keep common methods 
method_used <- Reduce(intersect,lapply(AUClist, names))
df_plot_list = lapply(AUClist, function(x) {x = x[method_used]}) 

df_plot = do.call(cbind, df_plot_list)
df_plot = melt(df_plot)
colnames(df_plot) <- c('method','dataset','AUC')

df_plot$method_type <- ifelse(df_plot$method %in% Methods$TraditionScore,'TraditionScore',
                               ifelse(df_plot$method %in% Methods$SenSet,'SenSet',
                                      ifelse(df_plot$method %in% Methods$SenMarker,'SenMarker',
                                             ifelse(df_plot$method %in% Methods$MachineScore,'MachineScore','hUSI'))))

df_plot$method_type <- factor(df_plot$method_type,levels = c('hUSI','SenMarker','SenSet','TraditionScore','MachineScore'),ordered = T)

df_plot_mean = aggregate(df_plot$AUC,by=list(df_plot$method),mean)
df_plot_mean$sd = aggregate(df_plot$AUC,by=list(df_plot$method),sd)$x
colnames(df_plot_mean) = c('method','mean','sd')
df_plot_mean$method_type = df_plot$method_type[match(df_plot_mean$method,df_plot$method)]

pdf("Results/Validation/compare_auc_sc2.pdf",width =8.5,height = 4)
ggplot()+
  geom_bar(data=df_plot_mean,aes(x=reorder(method,mean), y=mean,fill = method_type),stat="identity",position="dodge")+
  geom_errorbar(data=df_plot_mean,aes(x=reorder(method,mean),y=mean,ymax=mean+sd,ymin=ifelse(mean-sd <0,0,mean-sd)),
                position=position_dodge(0.9),width=0.4)+
  geom_jitter(data=df_plot,aes(x = method,y=AUC,color=dataset),width = 0.3,size=1.1,height = 0) +
  theme_classic()+
  ylab('AUC')+
  xlab('Senescence scoring method')+
  # scale_fill_manual(values = c('#e63946','#94d2bd','#00afb9','#fed9b7','#CFDCEF'))+
  scale_fill_manual(values = c('#D2352C','#71ABB6','#FFD699','#BAAFD1','#CFDCEF'))+
  # scale_color_manual(values = c('#277da1','#f9c74f','#7209b7','#68001D','#003560','#5C5C5C'))+
  scale_color_manual(values = c('#277da1','#fc8002','#8481BA','pink','#ADDB88','#F14454'))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1),
        text=element_text(size=14))
dev.off()

################################ in vivo samples ################################ 
## plot
fm="GTEx"
table(dat$label)

res=do.call(rbind,apply(dat[,c(-1,-ncol(dat)+1)], 2, function(x){
  dd=cor.test(x,as.numeric(as.factor(dat$label)),method = "pearson")
  dd=data.frame(cor=dd$estimate,pval=dd$p.value)
  return(dd)
})) %>% mutate(method=rownames(.))

library(ggrepel)
pdf(paste0("Results/Validation/",fm,"_cor_with_age_volcano.pdf"),width=5,height = 3.5)
ggplot(res, aes(x = cor, y = -log10(pval))) +
  geom_point(aes(size = -log10(pval), color = cor)) +
  geom_text_repel(data = subset(res, cor > 0.05), 
                  aes(label = method), size = 4, show.legend = FALSE) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),
        # axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  scale_color_gradientn(colours = c('skyblue3', "grey80", 'tomato3')) +
  labs(x = 'Correlation between score and age', y = '-log10(P value)')
dev.off()

## boxplot
dat = cbind(dat,GTEx_meta[,5,drop=F])

pdf(paste0("Results/Validation/",fm,"_cor_with_age.pdf"),width=5,height = 4.5)
ggplot(dat,aes(x=age,y=hUSI,fill=age))+
  geom_boxplot()+theme_classic()+
  geom_signif(comparisons = list(c("20-29","30-39"),c("30-39","40-49"),c("40-49","50-59"),c("50-59","60-69"),c("60-69","70-79")),
              step_increase=0.04, map_signif_level = T, test.args = c("less"))+
  theme(panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 14))+
  labs(x="Age groups")+scale_fill_brewer(palette = "Reds")+NoLegend()+
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) 
dev.off()

## each tissue
table(rownames(dat)==rownames(GTEx_meta))

### hUSI levels
df_plot = dat %>% mutate(group=ifelse(age %in% c('20-29', '30-39', '40-49'),"20-49","50-79"))
pwc <- df_plot %>% dplyr::group_by(SMTS) %>% 
  rstatix::pairwise_t_test(hUSI ~ group, paired = F,p.adjust.method = "bonferroni")
pwc
pwc <- pwc %>% rstatix::add_xy_position(x = "SMTS")

pdf(paste0("Results/Validation/",fm,"_tissue_hUSI.pdf"),width=9.5,height = 4)
df_plot %>%
  ggplot(aes(x = SMTS, y = hUSI)) + 
  geom_boxplot(aes(color = group),width=0.6,position = position_dodge(0.9), outlier.shape = NA, alpha = 0.6) +
  stat_pvalue_manual(pwc)+
  theme_classic()+
  theme(text = element_text(size = 12),axis.text.x = element_text(size=12,angle = 45,hjust = 1),
        axis.title.x = element_blank(),plot.margin = margin(0,0,0,20))+
  scale_color_manual(values = c('20-49' = '#959596','50-79'='#e63946'))+
  ylab('hUSI')+scale_y_continuous(breaks = seq(0,1,by=0.25))
dev.off()

### TCGA
library(ggunchained)
p=ggplot(sdata,aes(x = type,y = hUSI,fill = group)) +
  geom_split_violin(alpha = .6, trim = F,color = NA,width = 1.5) +
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.4),size=0.7) +
  stat_summary(fun = 'mean',
               geom = "pointrange",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x), width = .8,
               size = 0.2,position = position_dodge(0.4)) +
  theme_bw(base_size = 14.5) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        legend.position = 'top',panel.grid = element_blank()) +
  scale_fill_manual(values = c("#8CBDA7","#D2352C"))+labs(x="",y="hUSI")+
  ggpubr::stat_compare_means(aes(group=group),
                             symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                              symbols = c("***", "**", "*", "ns")),label = "p.signif",
                             label.y = 1,size = 4)+ylim(c(0,1.05))

pdf("Results/Validation/TCGA_NormalTumor_sen_score_paired.pdf",width =5.5,height = 4.3)
p
dev.off()

####################### Bulk skin
d=fm="GSE113957" #skin

pdf(paste('Results/Validation/benchmark_',d,'_tsne.pdf',sep=''),width = 13,height = 3.5)
p1= DimPlot(obj,group.by = 'label',reduction = 'tsne',cols = c('#91D1C2','#E64B35'),label = T,label.box = T,repel = T)+ggtitle(d)+NoLegend()
p2= FeaturePlot(obj,features = 'hUSI',cols = c('#00CCCC','#FF6E00'))+scale_color_gradient2(low ='#3AB370' ,mid = "#EAE7CC",high = "#FD1593",midpoint = 0.5)
p3= VlnPlot(obj,features = 'hUSI',group.by = 'label',cols =  c('#3AB370','#FD1593'),pt.size = 0)+
     geom_boxplot(width=0.2)+theme(axis.title = element_blank(),axis.text.x = element_blank())
p1+p2+p3
dev.off()

df_plot = obj@meta.data[,c('hUSI','label')]
pwc <- df_plot %>% 
  rstatix::pairwise_t_test(hUSI ~ label, paired = F,p.adjust.method = "bonferroni")
pwc

sum=obj@meta.data
subtwo=sum[,c("hUSI","age","Stage")] %>% filter(Stage=="normal")
# cor.test(subtwo$SCEID,subtwo$age)
cor.test(subtwo$hUSI,subtwo$age,method = "spearman")
pdf(paste0("Results/Validation/",fm,"_cor_with_age_fpkm.pdf"),width=4,height = 3.5)
ggplot(subtwo,aes(x=age,y=hUSI))+
  geom_point(size=2,color="tomato3")+theme_classic()+
  geom_smooth(method = 'lm',color="#71ABB6")+
  theme(panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),
        # axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  annotate("text",x=24,y=1,label="Cor=0.42, P=6.48e-07")
dev.off()

################################ in vivo cells ################################ 
### Pancreas
pdf(paste('Results/Validation/benchmark_',d,'_tsne.pdf',sep=''),width = 13,height = 3.5)
p1=DimPlot(obj,group.by = 'group',reduction = 'tsne',cols = c('#91D1C2','#FFD699','#E64B35'),label = T,label.box = T,repel = T)+ggtitle(d)+NoLegend()
p2=FeaturePlot(obj,features = 'hUSI',reduction = 'tsne')+
  scale_color_gradient2(low ='#3AB370' ,mid = "#EAE7CC",high = "#FD1593",midpoint = 0.5)
p3=VlnPlot(obj,features = 'hUSI',group.by = 'group',cols =  c('#91D1C2','#FFD699','#E64B35'),pt.size = 0)+
     geom_boxplot(width=0.2)+theme(axis.title = element_blank(),axis.text.x = element_blank())
p1+p2+p3
dev.off()

df_plot = obj@meta.data[,c('hUSI','group')]
pwc <- df_plot %>% 
  rstatix::pairwise_t_test(hUSI ~ group, paired = F,p.adjust.method = "bonferroni")
pwc

############################# in SA-beta-Gal dataset ############################ 
pdf(paste('Results/Validation/benchmark_',d,'_tsne.pdf',sep=''),width = 13,height = 3.5)
p1=DimPlot(obj,group.by = 'group',reduction = 'tsne',cols = c('#00CCCC','#91D1C2','#FFD699','#E2AE79','pink',"#F3766D",'#E64B35','#F14454'),label = T,label.box = T,repel = T)+ggtitle(d)+NoLegend()
p2=FeaturePlot(obj,features = 'hUSI',cols = c('#00CCCC','#FF6E00'))+
  scale_color_gradient2(low ='#3AB370' ,mid = "#EAE7CC",high = "#FD1593",midpoint = 0.5)
p3=VlnPlot(obj,features = 'hUSI',group.by = 'group',cols =  c('#00CCCC','#91D1C2','#FFD699','#E2AE79','pink',"#F3766D",'#E64B35','#F14454'),pt.size = 0)+
     geom_boxplot(width=0.2)+theme(axis.title = element_blank(),axis.text.x = element_blank())
p1+p2+p3
dev.off()

df_plot = obj@meta.data[,c('hUSI','group')]
pwc <- df_plot %>% 
  rstatix::pairwise_t_test(hUSI ~ group, paired = F,p.adjust.method = "bonferroni")
pwc

obj$hUSI = class_hUSI(obj$hUSI) ### write class function in function_new.R

sum=obj@meta.data

sum2 <- sum[,c("id","hUSI")] %>% as.data.frame() %>% 
  group_by(id) %>%
  dplyr::summarise(across(everything(), ~ sum(.) / n(), .names = "por_{col}")) %>% as.data.frame()

bt=fread("Data/meta_betaGal.csv") %>% as.data.frame()
sum2=merge(sum2,bt[,c("id","percent")],by="id")
cor.test(sum2$por_hUSI,sum2$percent,method = "spearman")

sum2$id=sapply(strsplit(sum2$id,"_"),function(x){paste0(x[[1]],"_",x[[2]])})
pdf(paste0("Results/Validation/",fm,"_cor_with_age.pdf"),width=5,height = 3.5)
ggplot(sum2,aes(x=percent,y=por_hUSI,color=id))+
  geom_point(size=2.5)+theme_bw()+
  geom_smooth(method = "lm",se=F,color="#71ABB6")+
  theme(panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))+
  scale_color_manual(values = c('#00CCCC','#91D1C2','#FFD699','#E2AE79','pink',"#F3766D",'#E64B35','#F14454'))+
  # scale_color_npg(alpha = 0.9)+
  annotate("text",x=0.3,y=1,label="Cor=0.87, P=1.246e-04")+
  labs(x = 'SA-beta-Gal positive ratio', y = 'hUSI senescent ratio',color="Samples")
dev.off()
