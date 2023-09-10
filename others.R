### Georgilis2018
Georgilis2018 = CreateSeuratObject(
    fread("Georgilis2018/valid_TPM_dataset.tsv") %>% column_to_rownames("Gene Name") %>% data.matrix, 
    meta.data = fread("Georgilis2018/filereport_read_run_PRJNA395386_tsv.txt", header = T) %>% column_to_rownames("sample_title")
)
Georgilis2018 = subset(Georgilis2018,read_count>1e6)
Georgilis2018 = NormalizeData(Georgilis2018, scale.factor = 1e4)
hSI = GetAssayData(Georgilis2018) %>% {apply( ., 2, function(z) {cor( z, mm_l2$w[ rownames(.) ], method="sp", use="complete.obs" )} )}

dat = calc_scores(ScoreList,Georgilis2018)
dat_Geo = dat %>% 
        rownames_to_column('sample') %>%
        mutate(condition=case_when(
        grepl("Allstars|Water|Hs|Extras[5-8]", sample) ~ "Senescence", 
        grepl("Noninduced|Extras[1-4]", sample) ~ "Growing", 
        TRUE ~ "other")) %>% 
        mutate(hSI=hSI) %>%
        filter( !condition %in% "other" ) %>% 
        mutate( condition=as.factor(condition))
dat_Geo = dat_Geo[,complete.cases(t(dat_Geo))]
auc_Geo <- calc_auc(dat_Geo,'marker',ScoreList)

auc %>%
    melt(value.name = "Accuracy") %>%
    mutate( feature=forcats::fct_reorder(Var1, Accuracy, mean, .desc = T) ) %>%
    Rmisc::summarySE(measurevar = "Accuracy", groupvars = "feature") %>%
    {
        ggplot(data = ., aes(feature, Accuracy, fill=feature)) +
        geom_bar(stat = "identity", position = position_dodge(), width = 0.4) +
        geom_errorbar(aes(ymin = Accuracy - sd, ymax = {Accuracy + sd} %>% ifelse(. >1, 1, .)),
        width = 0.2, position = position_dodge(0.9)) +
        labs(x=NULL, y="AUC") +
        theme(legend.position = "none")+
        theme_classic()
    }

### melanoma heatmap
mat <- exp[unlist(gene_set),rownames(meta)]
mat=mat[!apply(mat, 1, sd)==0,]
mat=Matrix::t(scale(Matrix::t(mat),center=TRUE))
mat=mat[is.na(row.names(mat)) == FALSE,]
mat[is.nan(mat)] = 0
mat[mat>2] = 2
mat[mat< -2] = -2


library(ComplexHeatmap)
library(circlize)
top_anno <- HeatmapAnnotation(df = data.frame(Condition = meta$condition),
                              col = list(Condition = c('young'='#007E99','senescent' = '#FF745A')),
                              show_legend = F)
left_anno = rowAnnotation(df = data.frame('State'= rep(c("Cycling","Moderate senescent","Senescent"),times=lapply(gene_set,length))),
                          show_legend = T,
                          col = list(State = bar))
col <- colorRamp2(c(-2,0,2), c("blue","white", "red"), space = "LAB")
png('/home/wangjing/codebase/HUSI/Figures/Melanoma_microarray.png',width = 1000,height = 800,res= 200)
Heatmap(mat,
        show_column_names = F,
        show_row_names = F,
        row_title = NULL,
        col = col,
        cluster_rows = T,
        cluster_row_slices = FALSE,
        cluster_columns = F,
        top_annotation = top_anno,
        left_annotation = left_anno,
        row_names_gp = gpar(fontsize = 12),
        column_split = c(rep('young',4),rep('senescent',4)),
        row_split = rep(c("Cycling","Moderate senescent","Senescent"),times=lapply(gene_set,length)))
dev.off()


### age state fraction in TCGA SKCM patients
library(ComplexHeatmap)
library(circlize)
marker_plot = unlist(lapply(marker_set, function(x) {rownames(x[order(x$avg_log2FC,decreasing = T),][0:200,])})) %>% unique
genes = intersect(rownames(Bulk),marker_plot)
mat <- as.matrix(Bulk[genes,])
mat <- apply(mat,1,scale)
mat <- t(mat)
rownames(mat) <- genes
colnames(mat) <- colnames(Bulk)

Frac_scale = data.frame(Frac)[colnames(mat),]
### scale Frac into 0-1 by column
Frac_scale <- apply(Frac_scale,2,function(x) {(x-min(x)) / (max(x)-min(x))})

top_anno <- HeatmapAnnotation(df = Frac_scale,
                              col = list(Cycling = colorRamp2(c(0, 1), c("white", "#5cb85c")),
                                         Transition = colorRamp2(c(0, 1), c( "white", "#428bca")),
                                         Senescent= colorRamp2(c(0, 1), c("white", "#d9534f"))),
                              show_legend = F)
col <- colorRamp2(c(-2,0,2), c("blue","white", "red"), space = "LAB")
png('/home/wangjing/wangj/codebase/HUSI/Figures/TCGA_SKCM_heatmap.png',width = 1500,height = 1000,res= 200)
Heatmap(mat,
        show_column_names = F,
        show_row_names = F,
        row_title = NULL,
        col = col,
        cluster_rows = T,
        cluster_row_slices = FALSE,
        cluster_columns = T,
        column_km=3,
        row_km = 3,
        top_annotation = top_anno,
        heatmap_legend_param = list(title = "Expression",
                                         at = c(-2,-1,0,1,2),
                                         labels = c("-2", "-1", "0","1", "2"),
                                         labels_gp = gpar(fontsize = 12),
                                         title_gp = gpar(fontsize = 12, fontface = "bold"),
                                         legend_width = unit(30, "mm")))
dev.off()