library(dplyr)

mm_l2 <- read.csv("HUSI/package/mm_l2.csv", header = TRUE, sep = ",")

cal_hUSI <- function(mat){
    score <- apply( ., 2, function(z) {cor(z, mm_l2$w[ rownames(.) ],method="sp", use="complete.obs" )})
    return(score)
}