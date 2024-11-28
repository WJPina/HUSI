##### This scripts is used to validate the performance of hUSI in vitro samples and in vitro cells
#################################### library ################################### 
library(data.table)
library(dplyr)
library(gelnet)
library(openxlsx)
library(annaffy)
library(magrittr)
library(ggplot2)
library(reshape2)
library(GEOquery)
library(tibble)
library(data.table)
library(dplyr)

setwd('hUSI/')

source('R/functions_new.R')
source('R/functions.R')

load('Data/SenOCLR_l2=1_drop.rdata')

################################### Benchmark #################################
set.seed(233)
### load methods functions
laflist<<-loadlafdata()
SenSetlist <<- loadSetdata()
load("Data/SENCAN_classifier.rda")
CSsccore <- read.csv('Data/CS score_PMID31560156.csv',header = T)
CSsccorelist <<- list('Underexpressed'=filter(CSsccore,expression=='Underexpressed')$Gene.symbol,
                      'Overexpressed'=filter(CSsccore,expression=='Overexpressed')$Gene.symbol)
coefslist <<- loadpredictorFix()

################################## in vitro samples ############################
####### microarray data sets
Arraymeta = read.csv('Data/ArrayMeta.csv',header = T,row.names = 1)

## gc-RMA (Robust Multi-Array Averaging with GC content correction),Probel level linear (PLM)
Arrayid = c('GSE19864','GSE16058','GSE83922','GSE11954','GSE100014','GSE77239')
GPL = c("GPL570" = "hgu133plus2.db","GPL3921" = "hthgu133a.db","GPL11532" = "hugene11sttranscriptcluster.db")
ArrayList <- list()
for(id in Arrayid){
  Array <- getGEO(id, destdir=".",AnnotGPL = F,getGPL = F)
  ArrayList[[id]] <- Array
}

ArrayexpList <- list()
for(id in names(ArrayList)){
  print(id)
  Array <- ArrayList[[id]]
  gene_id = aafSymbol(rownames(Array[[1]]), GPL[Array[[1]]@annotation]) %>% as.character
  mat<- Array[[1]] %>% 
    exprs %>% 
    set_rownames( gene_id ) %>% 
    melt(.) %>% 
    filter(!Var1 %in% "character(0)") %>% 
    as.data.table %>% 
    dcast.data.table( Var1~Var2, fun.aggregate = mean ) %>% 
    column_to_rownames("Var1")
  mat <- mat[,colnames(mat) %in% rownames(Arraymeta)]
  print(dim(mat))
  print(dim(Arraymeta[colnames(mat),]))
  ArrayexpList[[id]] <- CreateSeuratObject(mat,meta.data = Arraymeta[colnames(mat),])
}
names(ArrayexpList) <- c('PRJNA124029','PRJNA115473','PRJNA327413','PRJNA105815','PRJNA390490','PRJNA309835')
save(Arraymeta,ArrayexpList,file='Arraydata.rdata')

### methods used to compare
Methods = list('SenMarker' = c("GLB1", "TP53", "CDKN1A", "CDKN2A", "LMNB1", "IL1A", "RB1", "CDK1", "CDK4","CDK6", "MKI67", "CDKN2B",'SERPINE1'),
               'SenSet'=c("SenMayo","CellAge", "GenAge", "ASIG", "SASP","AgingAtlas", "SenUp","SigRS"),
               'MachineScore'=c('CSS','lassoCS','ECcores36','HCCcores19'),
               'TraditionScore'=c('DAS','mSS','DAS_mSS','CS_score'))

### calculate scores
datlist = {}
for(d in names(ArrayexpList)){
  print(d)
  obj = ArrayexpList[[d]]
  dat = calc_scores(Methods,obj,Repeat=10)
  hUSI = scoreOCLR(profile = GetAssayData(obj) %>% as.matrix,model = SenOCLR,m = 'spearman')
  dat$hUSI = hUSI[[1]] %>% minmax()
  dat$label=gsub('other','non-senescent',obj$condition)
  rownames(dat) = dat$cellname
  datlist[[d]] <- dat
}

### calculate auc
AUClist=list()
for(d in names(datlist)){
  dat = datlist[[d]]
  m = c('hUSI',as.character(unlist(Methods[c('SenSet','MachineScore','TraditionScore')])))
  auc_1 =  m%>%
    sapply(function(idx) {pROC::roc(dat$label, dat[,idx], levels=c('non-senescent','senescent'),direction = '<') %>% pROC::auc()})
 
  m = Methods[['SenMarker']] [Methods[['SenMarker']] %in% colnames(dat)] 
  auc_2 = m %>% 
    sapply(function(idx) {pROC::roc(dat$label, dat[,idx], levels=c('non-senescent','senescent'))%>% pROC::auc()})
  auc = c(auc_1,auc_2)
  
  AUClist[[d]] = auc
}
lapply(AUClist,sort)

############################### in vivo samples #################################
### GTEx
library(CePa)
## calculate scores
fm="GTEx"
GTEx_exp <- read.gct("bulk-gex_v8_rna-seq_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz") ### download from GTEx
Counts_gtex = ID2symbol(counts = GTEx_exp,database = 'genecode',version = '26',dir = 'Data/',GRCh = '38')
rownames(Counts_gtex) <- Counts_gtex$Symbol
Counts_gtex$Symbol <- NULL

### meta
GTEx_sample <- read.csv("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",sep = "\t",row.names = 1,header = T)
GTEx_sample<- GTEx_sample[,c("SMTS","SMTSD")]
GTEx_sample$Sample <- rownames(GTEx_sample)
GTEx_sample$Patient <- lapply(strsplit(rownames(GTEx_sample),'-'),function(x) {paste(x[1],x[2],sep = '-')}) %>% unlist
length(unique(GTEx_sample$Patient))

GTEx_Pheno <- read.table("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",sep = "\t",header = T)
length(unique(GTEx_Pheno$SUBJID))
GTEx_meta <- merge(GTEx_Pheno,GTEx_sample,by.x = "SUBJID",by.y = "Patient")
length(unique(GTEx_meta$SUBJID))
rownames(GTEx_meta) = gsub('-','\\.',GTEx_meta$Sample)

### subset data
samples = intersect(rownames(GTEx_meta),colnames(Counts_gtex))
Counts_gtex = Counts_gtex[,samples]
dim(Counts_gtex)
GTEx_meta=GTEx_meta[samples,]

### load methods functions
set.seed(233)
laflist<<-loadlafdata()
SenSetlist <<- loadSetdata()
load("Data/SENCAN_classifier.rda")
CSsccore <- read.csv('CS score_PMID31560156.csv',header = T)
CSsccorelist <<- list('Underexpressed'=filter(CSsccore,expression=='Underexpressed')$Gene.symbol,
                      'Overexpressed'=filter(CSsccore,expression=='Overexpressed')$Gene.symbol)
coefslist <<- loadpredictorFix()

### methods used to compare
Methods = list(
  'SenMarker' = c("GLB1", "TP53", "CDKN1A", "CDKN2A", "LMNB1", "IL1A", "RB1", "CDK1", "CDK4","CDK6", "MKI67", "CDKN2B",'SERPINE1'),
  'SenSet'=c("SenMayo","CellAge", "GenAge", "ASIG", "SASP","AgingAtlas", "SenUp","SigRS"),
  'MachineScore'=c('CSS','lassoCS','ECcores36','HCCcores19'),
  'TraditionScore'=c('DAS','mSS','DAS_mSS','CS_score'))

obj = CreateSeuratObject(counts = Counts_gtex,meta.data = GTEx_meta[colnames(Counts_gtex),],project = 'GTEx')
dat = calc_scores(Methods,obj)
hUSI = scoreOCLR(profile = Counts_gtex %>% as.matrix,model = SenOCLR,m = 'spearman')
dat$hUSI = hUSI[[1]] %>% minmax()
dat$label=obj$AGE
rownames(dat) = dat$cellname
dat$hUSI=hUSI[[1]] %>% minmax()


## count for SENCAN and SENCID
count=fread("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
count=count %>% as.data.frame() %>% select(-Name)# %>% column_to_rownames("Description")
count=count[!duplicated(count$Description),]
rownames(count)=count$Description
count=count[,-1]

### meta
GTEx_sample <- read.csv("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",sep = "\t",row.names = 1,header = T)
GTEx_sample<- GTEx_sample[,c("SMTS","SMTSD")]
GTEx_sample$Sample <- rownames(GTEx_sample)
GTEx_sample$Patient <- lapply(strsplit(rownames(GTEx_sample),'-'),function(x) {paste(x[1],x[2],sep = '-')}) %>% unlist
length(unique(GTEx_sample$Patient))

GTEx_Pheno <- read.table("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",sep = "\t",header = T)
length(unique(GTEx_Pheno$SUBJID))
GTEx_meta <- merge(GTEx_Pheno,GTEx_sample,by.x = "SUBJID",by.y = "Patient")
length(unique(GTEx_meta$SUBJID))
rownames(GTEx_meta) = gsub('-','\\.',GTEx_meta$Sample)

### subset data
colnames(count)=gsub('-','\\.',colnames(count))
samples = intersect(rownames(GTEx_meta),colnames(count))
count = count[,samples]
dim(count)
GTEx_meta=GTEx_meta[samples,]

obj = CreateSeuratObject(counts = count,project = 'GTEx')
Methods = list('MachineScore'=c('SENCAN','SENCID'))
dat_2 = calc_scores(Methods,obj)
dat$SENCAN=dat_2$SENCAN
dat$SenCID = dat_2$SENCID

### TCGA
### calculate scores
TCGA_exp = fread('tcga_RSEM_gene_tpm.gz') ### download from Xena 
TCGA_exp = data.frame(TCGA_exp)
TCGA_exp = column_to_rownames(TCGA_exp,'sample')

Counts_tcga = ID2symbol(counts = TCGA_exp,database = 'genecode',version = '23')
rownames(Counts_tcga) <- Counts_tcga$Symbol
Counts_tcga$Symbol <- NULL
Counts_tcga[1:5,1:5]

hUSI = scoreOCLR(profile = Counts_tcga %>% as.matrix,model = SenOCLR,m = 'spearman')
dd=data.frame(hUSI=hUSI[[1]] %>% minmax())
rownames(dd)=names(hUSI[[1]])

mrna=Counts_tcga
dmeta=data.frame(sample=colnames(Counts_tcga)) %>% mutate(bcr_patient_barcode=sapply(strsplit(sample,"[.]"), function(x){paste0(x[[1]],".",x[[2]],".",x[[3]])}),
                                                          type=sapply(strsplit(sample,"[.]"), function(x){x[[4]]}));head(dmeta)

meta=readxl::read_xlsx("Data/TCGA-CDR-SupplementalTableS1.xlsx") %>% 
  as.data.frame(meta)
meta$bcr_patient_barcode=gsub("-",".",meta$bcr_patient_barcode)
rownames(meta)=meta$bcr_patient_barcode

## Tumor vs Normal CS distribution
# only keep samples which have gene expression and clinical information
samples=data.frame(sample=colnames(mrna),group=substr(colnames(mrna),14,15)) %>%
  mutate(bcr_patient_barcode=substr(colnames(mrna),1,12)) %>%
  merge(.,meta[,2:3],by="bcr_patient_barcode") %>%
  mutate(group=ifelse(group>=10,"Normal","Tumor"))#10496

mrna_cs<- as.data.frame(hUSI)
colnames(mrna_cs)=c("sample","hUSI")

library(ggunchained)
## paried
# Tumor
Tumor_group <- samples[which(samples$group == "Tumor"),] 
# Normal
Normal_group <- samples[which(samples$group == "Normal"),] 
# paired
Normal_group <- Normal_group[Normal_group$bcr_patient_barcode %in% Tumor_group$bcr_patient_barcode,]

paired_expr_data <- mrna_cs[substring(mrna_cs$sample,1,12) %in% Normal_group$bcr_patient_barcode,] 
paired_phe <- samples[match(paired_expr_data$sample,samples$sample),]
paired_phe <- droplevels(paired_phe)
table(paired_phe$type)

cancer=table(paired_phe$type)[which(table(paired_phe$type)>=20)] %>% names()
sdata=merge(paired_expr_data,samples,by="sample") %>% filter(type %in% cancer)
sdata$group=factor(sdata$group,levels = c("Tumor","Normal"))

####################### Bulk skin
d=fm="GSE113957" #skin
exp=fread("GSE113957_fpkm.txt.gz") %>% as.data.frame() ### download from GEO
gene=exp[,c("Annotation/Divergence","Transcript ID")]
colnames(gene)[1]="geneid"
gene=gene %>% mutate(gene=sapply(strsplit(geneid,"[|]"),function(x){x[[1]]}),type=sapply(strsplit(geneid,"[|]"),function(x){x[[length(x)]]}))
rownames(gene)=gene$`Transcript ID`
gene=gene[!duplicated(gene$gene),]

exp=exp[,-2:-8]
exp=exp %>% column_to_rownames("Transcript ID")
exp=exp[gene$`Transcript ID`,]
rownames(exp)=gene$gene

fpkmToTpm <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpm <- apply(exp,2,fpkmToTpm)

meta= fread("GSE113957_meta.csv",header = T) %>% as.data.frame() %>% 
  mutate(V1=gsub("-",".",V1))
rownames(meta) = meta$sample_title.y
meta=meta[colnames(exp),]
table(meta$sample_title.y==colnames(exp))

obj = CreateSeuratObject(counts = log1p(tpm),meta.data = meta)
obj <- FindVariableFeatures(obj) %>% ScaleData() %>% RunPCA() %>% RunTSNE(dims=1:20)

hUSI = scoreOCLR(profile = tpm %>% as.matrix,model = SenOCLR,m = 'spearman')
obj$hUSI=minmax(hUSI[[1]])
obj$label=ifelse(obj$age<51,"1-50","51-96")

################################ in vitro cells ################################ 
### methods used to compare
Methods = list('SenMarker' = c("GLB1", "TP53", "CDKN1A", "CDKN2A", "LMNB1", "IL1A", "RB1", "CDK1", "CDK4","CDK6", "MKI67", "CDKN2B",'SERPINE1'),
               'SenSet'=c("SenMayo","CellAge", "GenAge", "ASIG", "SASP","AgingAtlas", "SenUp","SigRS"),
               'MachineScore'=c('CSS','lassoCS','ECcores36','HCCcores19','SENCAN','SENCID'),
               'TraditionScore'=c('DAS','mSS','DAS_mSS','CS_score'))

### load single cell data set
scDataSets = loadSCdata()
sapply(scDataSets,dim)
### calculate scores
datlist = {}

for(d in names(scDataSets)){
  print(d)
  obj = scDataSets[[d]]
  dat = calc_scores(Methods,obj)
  hUSI = scoreOCLR(profile = GetAssayData(obj) %>% as.matrix,model = SenOCLR,m = 'spearman')
  dat$hUSI = hUSI[[1]] %>% minmax()
  dat$label=obj$Condition
  rownames(dat) = dat$cellname
  datlist[[d]] <- dat
}

### calculate auc
AUClist=list()
for(d in names(datlist)){
  dat = datlist[[d]]
  m = c('hUSI',as.character(unlist(Methods[c('SenSet','MachineScore','TraditionScore')])))
  auc_1 =  m%>%
    sapply(function(idx) {pROC::roc(dat$label, dat[,idx], levels=c('non-senescent','senescent'),direction = '<') %>% pROC::auc()})
 
  m = Methods[['SenMarker']] [Methods[['SenMarker']] %in% colnames(dat)] 
  auc_2 = m %>% 
    sapply(function(idx) {pROC::roc(dat$label, dat[,idx], levels=c('non-senescent','senescent'))%>% pROC::auc()})
  auc = c(auc_1,auc_2)
  
  AUClist[[d]] = auc
}
lapply(AUClist,sort)

################################ in vivo cells ################################ 
### Pancreas
library(tibble)
library(edgeR)
d=fm="GSE81547_pancreas" 

exp = fread("GSE81547_pancreas_count.csv") ### download from GEO
exp = as.data.frame(exp)
exp = column_to_rownames(exp,'V1')
dim(exp)
exp[1:5,1:5] 

meta= fread("GSE81547_pancreas_meta.csv",header = T) %>% as.data.frame() %>% 
  mutate(V1=gsub("-",".",V1))
rownames(meta) = meta$V1
table(meta$V1==colnames(exp))
obj=CreateSeuratObject(exp,meta.data = meta)
obj <- NormalizeData(obj)%>%  FindVariableFeatures() %>% ScaleData() %>% RunPCA()
obj=RunTSNE(obj,dims=1:10)
obj=RunUMAP(obj,dims = 1:10)

hUSI = scoreOCLR(profile = GetAssayData(obj) %>% as.matrix,model = SenOCLR,m = 'spearman')
obj$hUSI = minmax(hUSI[[1]])
obj$group=ifelse(obj$donor_age<7,"1-6",ifelse(obj$donor_age<23,"21-22","38-54"))

############################# in SA-beta-Gal dataset ############################ 
fm = 'GSE175533_SAbetaG'
exp = fread('GSE175533_SAbetaG_norm.csv') ### download from GEO
exp = as.data.frame(exp)
exp = column_to_rownames(exp,'V1')
exp=t(exp)
colnames(exp)=gsub("-",".",colnames(exp)) 
colSums(expm1(exp[,1:4]))

meta= fread("GSE175533_SAbetaG_meta.csv",header = T) %>% as.data.frame() %>% 
  mutate(V1=gsub("-",".",V1))
rownames(meta) = meta$V1
table(meta$V1==colnames(exp))

obj=CreateSeuratObject(exp,meta.data = meta)#logTP10k
obj <- FindVariableFeatures(obj) %>% ScaleData() %>% RunPCA() %>% RunTSNE(dims=1:20)
hUSI = scoreOCLR(profile = GetAssayData(obj) %>% as.matrix,model = SenOCLR,m = 'spearman')
obj$hUSI = minmax(hUSI[[1]])
obj$group = factor(obj$PDL,levels = c('htert_2','htert_7','PDL_25','PDL_29','PDL_33','PDL_37','PDL_46','PDL_50' ))

####################### new generated RNA-seq data #############################
## load TPM
samples=fread("Data/GeneratedRNA-seqData/OIS_meta_file.csv") %>% as.data.frame()

data <- fread(paste0("GeneratedRNA-seqData/4OHT_1/geneExpr.stringtie.tab")) %>% ### download from GEO GSE282274
  select(`Gene Name`, TPM) %>%  
  distinct(`Gene Name`, .keep_all = TRUE) %>%  
  as.data.frame() %>% column_to_rownames('Gene Name') 
colnames(data)='4OHT_1'

for (i in samples$run_accession[-1]) {
  tmp <- fread(paste0("GeneratedRNA-seqData/", i, "/geneExpr.stringtie.tab")) %>%
    select(`Gene Name`, TPM) %>%  
    distinct(`Gene Name`, .keep_all = TRUE) %>%  
    as.data.frame() %>% column_to_rownames('Gene Name') 
  colnames(tmp)=i
  tmp=tmp[intersect(rownames(tmp),rownames(data)),,drop=F]
  data=data[intersect(rownames(tmp),rownames(data)),,drop=F]
  table(rownames(data)==rownames(tmp))
  data=cbind(data,tmp)
  print(i)
}

hUSI  = scoreOCLR(data,SenOCLR,'spearman')[[1]] %>% minmax()

meta=data.frame(sample=samples$run_accession,group=c(rep("4OHT",3),rep("DMSO",3)),hUSI=hUSI[samples$run_accession])
meta$group=factor(meta$group,levels = c("DMSO","4OHT"))

ggplot(meta,aes(x=group,y=hUSI,fill=group))+
  geom_boxplot()+theme_bw()+
  geom_jitter(color="tomato3",width = 0.31)+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12))+scale_fill_aaas()+NoLegend()+
  ggsignif::geom_signif(comparisons = list(c("DMSO","4OHT")),test = "t.test")+
  labs(x="Group")
tpmRNA=data
