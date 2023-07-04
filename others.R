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