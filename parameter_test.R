load('~/wangj/AgingScore/Data/Bulk_TrainModel/ModelTrainData.RData')
library(gelnet)
library(dplyr)
# Log-transform data and standardize
X = log2(exprData + .00001) %>%  scale()
# ### mean center and split to train and background dataset
X_centre = X - (apply(X, 1, mean))

idx_senescent = colnames(X) %in% filter(metadata, condition %in% "senescent")$sample_title
X_tr = X_centre[,idx_senescent]
X_bk = X_centre[,!idx_senescent]

### parameter test Leave-one-out cross validation
### Leave-one-out cross validation
auclist <- list()
for(l2 in c(0,0.5,1,1.5,2)){
  auc <- c()
  for(i in 1:ncol(X_tr)){
    ## Train a model on non-left-out data
    X1 <- X_tr[,-i]
    m1 <- gelnet( t(X1), NULL, 0, l2)
    ## Score the left-out sample against the background
    s_bk <- apply( X_bk, 2, function(z) {cor( m1$w, z, method="sp" )} )
    s1 <- cor( m1$w, X_tr[,i], method="sp" )
    ## CRP = P( left-out sample is scored above the background )
    auc[i] <- sum( s1 > s_bk ) / length(s_bk)  
    cat( "Current auc: ", auc[i], "\n" )
    cat( "Average auc: ", mean(auc), "\n" )
  }
  auclist[[as.character(l2)]] <- auc
}

save(auclist,file='parameter_test_auclist.RData')