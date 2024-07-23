mm_l2 <<- readRDS("~/wangj/AgingScore/Data/Bulk_TrainModel/mm_l2.rds")
######################## Comparision in three single-cell dataset  ###################################
library(Seurat)
library(data.table)
library(tibble)
library(sparseMatrixStats)
library(org.Hs.eg.db)
library(dplyr)
library(caret)
library(magrittr)
library(reshape2)
source('~/wangj/codebase/HUSI/functions.R')
set.seed(233)

cal_hUSI <- function(object){
  score <- GetAssayData(object) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )} %>% minmax
  return(score)
}

cal_marker <- function(object,markers){
  markers = markers[markers%in%rownames(object)]
  scores = object[markers,] %>% GetAssayData %>% as.matrix %>% t %>% data.frame %>% as.list()
  scores = lapply(scores, function(x){names(x) = colnames(object);return(x)})
  return(scores)
} 

cal_hUSI_mean <- function(object,k=500){
  exp=GetAssayData(object)
  scores <- list()
  all_genes = sort(mm_l2$w,decreasing = T) %>% names
  all_genes = all_genes[all_genes %in% rownames(object)]
  top_genes = all_genes[1:k] 
  for(i in 1:length(top_genes)){
    id = paste('top_',i,sep='')
    genes <- top_genes[1:i]
    if(i == 1){scores[[id]] <- exp[genes,]}
    else(scores[[id]] <- colMeans(exp[genes,]))
  }
  return(scores)
}

load('~/wangj/AgingScore/Data/Bulk_TrainModel/ModelTrainExp.RData')
cal_hUSI_exp <- function(object,train_exp){
  weights = rowMeans(train_exp)
  score <- GetAssayData(object) %>% {apply( ., 2, function(z) {cor( z, weights[ rownames(.) ], method="sp", use="complete.obs" )} )} %>% minmax
  return(score)
}
################################## load benchmark data 
### Teo2019  GSE115301 IMR90 Smart-seq2 OIS,3',counts
Teo2019 = CreateSeuratObject(
  fread("Teo2019/GSE115301_Growing_Sen_10x_count.txt.gz") %>% column_to_rownames("V1") %>% data.matrix, 
  meta.data = fread("Teo2019/GSE115301_Growing_Sen_10x_metadata.txt.gz", header = T) %>% column_to_rownames("V1")
) %>% NormalizeData()

Teo2019$Condition = as.factor(ifelse(Teo2019$Condition2 == 'Growing','Growing','Senescence'))

### Tang2019 GSE119807 HCA2 fibroblast cell RS,IRIS,3',counts
Tang2019List = list.files("Tang2019", full.names = T) %>% 
  lapply(function(x) {
    scdat = fread(x) %>% 
      column_to_rownames("GENE") %>% 
      data.matrix }) %>% set_names( list.files("Tang2019") %>% gsub(".*_", "", .) %>% gsub("\\..*", "", .) )

genes = Reduce(intersect,lapply(Tang2019List, function(x)rownames(x)))
Tang2019 = do.call(cbind,lapply(Tang2019List,function(x) x <- x[genes,]))
Tang2019 = CreateSeuratObject(Tang2019) %>% NormalizeData()
Tang2019$Condition = rep(names(Tang2019List),each=400)
Tang2019$Condition = as.factor(ifelse(Tang2019$Condition %in% c("senescence","LowPD50Gy"),'Senescence','Growing'))

### Aarts2017 GSE94980 IMR90 OSKM-expressing reprogramming-induced senescence,3',counts
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
Aarts2017List = list.files("Aarts2017", pattern = "*.gz",full.names = T)
scdat = data.frame(gene = '')
for(file in Aarts2017List){
  scol = fread(file) 
  colnames(scol) = c('gene',paste(strsplit(basename(file),"_")[[1]][c(2:5)],collapse='_'))
  scdat = right_join(scdat,scol,by = 'gene')           
}
ids=bitr(scdat$gene,'ENSEMBL','SYMBOL','org.Hs.eg.db')
exp=merge(scdat,ids,by.x='gene',by.y='ENSEMBL')
exp = exp[-1]
exp = exp[!duplicated(exp$SYMBOL),]
rownames(exp) <- exp$SYMBOL
exp = exp[-ncol(exp)]
dim(exp)

Aarts2017 = CreateSeuratObject(exp) %>% NormalizeData()
Aarts2017$Condition = ifelse(grepl('OSKM',colnames(Aarts2017)), 'Senescence', 'Growing')

########################################### calculate score
### score list
SenMarkers <<- c("CDKN1A", "CDKN2A", "CDKN2B",'SERPINE1')

objectList = list(Teo2019,Tang2019,Aarts2017)
names(objectList) <-  c('Teo2019','Tang2019','Aarts2017')

dat = list()
for(dataset in names(objectList)){
  obj = objectList[[dataset]]
  hUSI = cal_hUSI(obj)
  sene_marker = cal_marker(obj,SenMarkers)
  # hUSI_mean = cal_hUSI_mean(obj)
  # hUSI_exp = cal_hUSI_exp(obj,X_tr)
  dat[[dataset]] = c(list('hUSI'=hUSI),sene_marker)
}

### auc
AUClist=list()
for(dataset in names(dat)){
  data = dat[[dataset]] %>% data.frame()
  data$Condition = objectList[[dataset]][,rownames(data)]$Condition
  fold = createMultiFolds(data$Condition, k = 10, times = 3)
  auc = sapply(fold, function(sampling){
    df = data[sampling,]
    auc = names(dat[[dataset]]) %>%
      sapply(function(idx) {pROC::roc(df$Condition, df[,idx], levels=c("Growing", "Senescence"),direction = '<') %>% pROC::auc()})
    return(auc)
  })
  AUClist[[dataset]] = auc
}

### plot senesence markers
### mean auc values
df_plot_list = lapply(AUClist[c("Teo2019","Tang2019","Aarts2017")], function(x) {melt(x[c('hUSI',SenMarkers),]) %>% set_names(c('method','unit','AUC'))})

plot_list = list()
for (dataset in names(df_plot_list)) {
  df_plot = df_plot_list[[dataset]]
  df_plot$method_type <- ifelse(df_plot$method %in% SenMarkers,'Marker','hUSI')
  # df_plot$method_type <- factor(df_plot$method_type,levels = c('hUSI','Marker'),ordered = T)
  # df_plot$method_type <- factor(df_plot$method,levels = c('hUSI','hUSI_exp'),ordered = T)
  df_plot$dataset = dataset
  
  plot_list[[dataset]] <-
    ggplot(df_plot,aes(x=reorder(method,AUC), y=AUC,color = method_type))+
    geom_jitter(size=0.5)+
    geom_boxplot(color='black',outlier.color = NA,fill=NA,linewidth=0.1)+
    theme_classic()+
    ylab('AUC values')+
    # xlab('Quantification method')+
    xlab('')+
    scale_color_manual(values = c('#e63946','#94d2bd'))+
    theme(panel.grid.major.y = element_line(),
          axis.text.x = element_text(angle = 45,vjust = 1, hjust=1),
          text=element_text(size=16))+
    ggtitle(dataset)+ylim(0,1)
}


fig <- ggarrange(plotlist = plot_list,ncol = 3,nrow = 1,common.legend = T)

png('~/wangj/codebase/HUSI/Figures/revison/compare_auc_exp.png',width = 2000,height = 1200,res = 300)
fig
dev.off()

#### mean auc rank 
### only keep common genes
method_used <- c('hUSI',SenMarkers)
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
df_plot_mean$method_type <- ifelse(df_plot_mean$method %in% SenMarkers,'Marker','hUSI')
df_plot_mean$group <- factor(df_plot_mean$method_type,levels = c('hUSI','Marker'),ordered = T)

png('~/wangj/codebase/HUSI/Figures/revison/compare_rank_sene.png',width = 1200,height = 1500,res = 300)
ggplot()+
  geom_bar(data=df_plot_mean,aes(x=reorder(method,mean), y=mean,fill = group),stat="identity",position="dodge")+
  geom_jitter(data=df_plot,aes(x = method,y=AUC_rank,color=dataset),width = 0.3,size=2,height = 0) +
  geom_errorbar(data=df_plot_mean,aes(x=reorder(method,mean),y=mean,ymax=mean+sd,ymin=ifelse(mean-sd <0,0,mean-sd)),
                position=position_dodge(0.9),width=0.4)+
  theme_classic()+
  ylab('AUC rank')+
  # xlab('Quantification method')+
  xlab('')+
  scale_fill_manual(values = c('#e63946','#94d2bd','#00afb9','#fed9b7'))+
  scale_color_manual(values = c('#277da1','#f9c74f','#7209b7'))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1),
        text=element_text(size=16))
dev.off()

### plot top mean and expression
plot_list = list()
for(dataset in names(AUClist)){
  df_plot = AUClist[[dataset]][c('hUSI','hUSI_exp',names(hUSI_mean)),]
  df_plot_mean = sort(rowMeans(df_plot))
  df_plot_mean = data.frame(AUC_mean=df_plot_mean,top_mean=names(df_plot_mean)) 
  labeldata = rbind(filter(df_plot_mean,top_mean %in% c('hUSI','hUSI_exp')),
                    df_plot_mean[names(hUSI_mean),][df_plot_mean[names(hUSI_mean),]$AUC_mean == max(df_plot_mean[names(hUSI_mean),]$AUC_mean),][1,])
  labeldata = labeldata[!duplicated(labeldata$top_mean),]
  plot_list[[dataset]] <- 
  ggplot(df_plot_mean,aes(x=reorder(top_mean,AUC_mean),y=AUC_mean)) + geom_point(size=0.2,color='grey60')+
    theme_classic()+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
    geom_point(data = labeldata,size=1,aes(x=top_mean,y=AUC_mean),colour="red")+
    geom_text_repel(data = labeldata,label = rownames(labeldata),max.overlaps = 10,colour = "black",size = 5,show.legend = F)+
    xlab('Quantidcation method')+ylab('mean AUC value')+
    ylim(0,1)+ggtitle(dataset)
}

fig <- ggarrange(plotlist = plot_list,ncol = 1,nrow = 3,common.legend = T)

png('~/wangj/codebase/HUSI/Figures/revison/compare_auc_weights_expression.png',width = 4000,height = 3000,res = 300)
fig
dev.off()


### quantification methods
cal_logistic <- function(object){
  exp <- GetAssayData(object)
  genes <- intersect(rownames(exp),names(mm_l2$w))
  w = mm_l2$w[genes]
  exp = exp[genes,]
  score <- apply(exp,2,function(x){
    likelihood <- exp(w %*% x)/(1+exp(w %*% x));
    return(likelihood)
  })
  return(score)
}

objectList = list(Teo2019,Tang2019,Aarts2017)
names(objectList) <-  c('Teo2019','Tang2019','Aarts2017')

dat = list()
for(dataset in names(objectList)){
  obj = objectList[[dataset]]
  hUSI = cal_hUSI(obj)
  hUSI_logistic = cal_logistic(obj)
  dat[[dataset]] = list(cell=colnames(obj),'hUSI'=hUSI,'hUSI_logistic'=hUSI_logistic)
}


plot_list = list()
for (dataset in names(dat)) {
  df_plot = data.frame(hUSI=dat[[dataset]]$hUSI,hUSI_logistic=dat[[dataset]]$hUSI_logistic)
  plot_list[[dataset]] <-
    df_plot %>%
    ggscatter(x = "hUSI", y = "hUSI_logistic",
              color = 'grey30',
              add = "reg.line", 
              conf.int = TRUE,
              size = 1,
              add.params = list(color = "#6c757d"),
              ggtheme = theme_classic())+
    stat_cor(method = "spearman",color='black')+
    ggtitle(dataset)
}


fig <- ggarrange(plotlist = plot_list,ncol = 3,nrow = 1,common.legend = T)

png('~/wangj/codebase/HUSI/Figures/revison/compare_quntification_1.png',width = 3300,height = 1200,res = 300)
fig
dev.off()

openxlsx::write.xlsx(dat,"~/wangj/codebase/HUSI/Figures/revison/compare_quntification_1.xlsx")


### show weights distribution
EnrichSet<<-cogena::gmt2list("gene_50signatures_merge.gmt")
SeneGenes <- EnrichSet$SenMayo

SenMarkers <<- c("GLB1", "CDKN1A", "CDKN2A",  "IL1A", "CDKN2B",'SERPINE1')
SeneGenes <- unique(c(SeneGenes,SenMarkers))

weights = mm_l2$w[SeneGenes]
weights = weights[!is.na(weights)]
match(names(weights),names(mm_l2$w))

df_plot = data.frame(weight = mm_l2$w,set = ifelse(names(mm_l2$w)%in%SeneGenes,'SenMayo','Other'))
df_plot$gene <- rownames(df_plot)
top = filter(df_plot,set=='SenMayo')
top = top[order(top$weight,decreasing = T),]


  
  png('~/wangj/codebase/HUSI/Figures/revison/compare_quntification_3.png',width = 3300,height = 1000,res = 300)
  ggplot(df_plot)+
    geom_density(aes(x=weight),alpha=0.5)+
    scale_fill_manual(values = c('#457b9d','#e63946'))+
    geom_vline(xintercept =df_plot$weight[which(df_plot$gene %in%c('CDKN1A','SERPINE1','MKI67','LMNB1'))])+
    geom_vline(xintercept=quantile(df_plot$weight,c(0.1,0.9)),color='red')+
    theme_bw()
  dev.off()

### weights distribution
load('../Data/Bulk_TrainModel/ModelTrainData.RData')



