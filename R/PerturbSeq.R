### This scrpts is for novel senescence regulator detection in Perturb-seq dataset (Fig3 and FigS4)

library(ggvenn)
library(openxlsx)
library(ggrepel)

### Calculate hUSI
# scPerturb_RPE1_hUSI.ipynb

### Other analysis and visualization
dd='ReplogleWeissman2022_rpe1'
res_husi=fread(paste0('/mnt/data2/zhouxl/hUSI/Results_data/RPE1_hUSI_',dd,'_fc_SenOCLR_l2_1_drop.csv'),sep='\t')[,-1] %>% as.data.frame()
meta=fread(paste0('/mnt/data2/zhouxl/hUSI/Results_data/RPE1_hUSI_',dd,'_meta_SenOCLR_l2_1_drop.csv'),sep='\t')[,-1] %>% as.data.frame()

score=meta %>% group_by(perturbation) %>% dplyr::summarise(husi=mean(hUSI))
res_husi=merge(res_husi,score,by.x="group",by.y='perturbation') %>% arrange(desc(deltaMean))

res_husi2=res_husi %>% filter(fc>1.5,padj<0.05) %>% arrange(desc(deltaMean));dim(res_husi2);head(res_husi2,10)
openxlsx::write.xlsx(res_husi2,"./RPE1_hUSI_RPE1_fc_SenOCLR_l2_1_drop_signif_fc1.5_1057genes.csv")## Supplementary data 9

## TOP10
tmp=res_husi2$group[1:10]
submeta=meta %>% filter(perturbation %in% c(tmp,"control"))
submeta$perturbation=factor(submeta$perturbation,levels = c("control",tmp))

mycol=c(pal_nejm(alpha = 0.8)(8),'#E5D2DD', '#53A85F', '#F1BB72')
pdf("/mnt/data2/zhouxl/hUSI/Figures_new/Validation/scPerturb_seq_RPE1_boxplot2.pdf",width = 4.5,height = 3.5)
ggplot(submeta,aes(x=perturbation,y=hUSI,fill=perturbation))+
  geom_boxplot()+theme_classic()+
  # ggsignif::geom_signif(comparisons = list(c("control",tmp[1]),c("control",tmp[2]),c("control",tmp[3]),c("control",tmp[4]),
  #                                          c("control",tmp[5]),c("control",tmp[6]),c("control",tmp[7]),c("control",tmp[8]),
  #                                          c("control",tmp[9]),c("control",tmp[10])),map_signif_level = T,step_increase = 0.07)+
  theme(panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))+NoLegend()+ylim(c(0,1))+scale_fill_manual(values = mycol)+#+scale_fill_nejm(alpha = 0.8)
  labs(x="Perturbation")
dev.off()


## enrichment
KEGG = read.gmt('/mnt/data1/wangj/GeneSets/KEGG_hsa.gmt')
databases=list('KEGG'=KEGG)
res=res_husi$deltaMean;names(res)=res_husi$group

fgseaList=list()
for(d in names(databases)){
  print(d)
  set.seed(233)
  pathways = databases[[d]]
  pathways$term = sub("^[^_]*_", "", pathways$term)
  fgsea <- GSEA(sort(res,decreasing=T), 
                exponent = 1, 
                minGSSize = 3,
                maxGSSize = 1000, 
                pvalueCutoff = 1, 
                pAdjustMethod = "BH", 
                TERM2GENE=pathways,
                seed = 233,
                by = "fgsea",
                eps=0)
  fgseaList[[d]] <- fgsea
}

re = data.frame()
for(d in names(databases)){
  print(d)
  df = fgseaList[[d]]@result
  df = filter(df, pvalue < 0.05)
  df = df[order(df$NES,decreasing=T),]
  df$Database = d
  print(nrow(df))
  re = rbind(re,df)
}

d='KEGG'
df = filter(re,Database==d)
df[,c('NES','p.adjust')]

### hallmarker pathway
ids = c('DNA replication',
        'DNA replication proteins',
        'Cell cycle','p53 signaling pathway')

pdf('/mnt/data2/zhouxl/hUSI/Figures_new/Validation/scPerturb_enrichment.pdf',width =6,height =4.5)
enrich_plot(gsea=fgseaList[[d]],data=df,pathway_id=ids,direction = 'up',save = F,database=d)
dev.off()

################## Nilab RNA-seq
## load TPM
samples=fread("/mnt/data3/zhouxl/Qiushi/RNA_seq/meta_file.csv") %>% as.data.frame()

data <- fread(paste0("/mnt/data3/zhoumengyu/data/N2308074_80-1220641544_RNA_2023-06-05/aging/01.alignment/step2/4OHT_1/geneExpr.stringtie.tab")) %>%
  select(`Gene Name`, TPM) %>%  
  distinct(`Gene Name`, .keep_all = TRUE) %>%  
  as.data.frame() %>% column_to_rownames('Gene Name') 
colnames(data)='4OHT_1'

for (i in samples$run_accession[-1]) {
  tmp <- fread(paste0("/mnt/data3/zhoumengyu/data/N2308074_80-1220641544_RNA_2023-06-05/aging/01.alignment/step2/", i, "/geneExpr.stringtie.tab")) %>%
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

source('/home/wangjing/wangj/AgingScore/HUSI_new/codes/functions_new.R')
load('/home/wangjing/wangj/AgingScore/HUSI_new/data_final/SenOCLR_l2=1_drop.rdata')

hUSI  = scoreOCLR(data,SenOCLR,'spearman')[[1]] %>% minmax()

meta=data.frame(sample=samples$run_accession,group=c(rep("4OHT",3),rep("DMSO",3)),hUSI=hUSI[samples$run_accession])
meta$group=factor(meta$group,levels = c("DMSO","4OHT"))

pdf("/mnt/data2/zhouxl/hUSI/Figures_new/Validation/NiLab_RNA_seq_hUSI.pdf",width = 3,height = 4)
ggplot(meta,aes(x=group,y=hUSI,color=group))+
  geom_boxplot()+theme_classic()+scale_color_manual(values = c('#91D1C2','#E64B35'))+
  geom_jitter(color="#8481BA",width = 0.31)+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12))+scale_fill_aaas()+NoLegend()+
  ggsignif::geom_signif(comparisons = list(c("DMSO","4OHT")),test = "t.test",color="black")+
  labs(x="Group")+scale_y_continuous(breaks = seq(0,1.08,0.25))

dev.off()

tpmRNA=data
save(tpmRNA,meta,file="/mnt/data2/zhouxl/hUSI/data/Qiushi_RNA_seq/RNA_seq_TPM.RData")## can download from GSE282274








