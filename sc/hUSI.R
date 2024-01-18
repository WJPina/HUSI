#####################################################
# This script is used to calculate the hUSI score and 
# fit the hUSI score by gassian mixtrue model in R.
#####################################################

library(dplyr)
library(mclust)

set.seed(233)

cal_hUSI <- function(exp){
    mm_l2 <- read.csv("/home/wangjing/wangj/codebase/HUSI/mm_l2.csv", header = TRUE, sep = ",",row.names=1)
    genes = intersect(rownames(exp), rownames(mm_l2))
    w = mm_l2[genes,'w']
    exp = exp[genes,]
    score <- apply(exp, 2, function(x) {cor(x, w,method="sp", use="complete.obs" )})
    return(score)
}

fit_hUSI <- function(score,g=NULL){
    if(is.null(g)){
        gaussian = score %>% {log2((1+ .)/(1- .))} %>% Mclust()
    }
    else {
       gaussian = score %>% {log2((1+ .)/(1- .))} %>% Mclust(G=g)
    }
    age_class <- gaussian$classification
    return(age_class)
}