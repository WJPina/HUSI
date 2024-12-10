load('sc/SenOCLR_features.rdata')

### calculate hUSI scores
minmax <- function(x){
  x_scale = (x-min(x))/(max(x)-min(x))
  return(x_scale)
}
cal_hUSI <- function(profile){
    w = SenOCLR$w
    genes = intersect(names(w),rownames(profile))
    w = w[genes]
    profile = profile[genes,] %>% as.matrix()
    scores = apply(profile, 2, function(s) {cor(w , s, method="spearman" )}) %>% minmax()
    names(scores) = colnames(profile)
    return(scores)
}
### SSE
SSE_hUSI <- function(scores){
    if(is.null(names(scores))){
        print("hUSI scores must be a named vector!!!")
    }
    else {
        th = NA
        minSSE = NA
        for (i in scores){
            low = scores[scores <= i]
            high = scores[scores > i]
            sse = sum((low-mean(low))^2) + sum((high-mean(high))^2)
            if(is.na(th) & is.na(minSSE)){th = i;minSSE = sse}
            else if (sse < minSSE) {th = i;minSSE = sse}
            else {th = th;minSSE = minSSE}
        }
        senclass = ifelse(scores <= th,0,1)
        names(senclass) = names(scores)
        return(senclass)
    }
}
### GMM
suppressMessages(library(mclust))
suppressMessages(library(dplyr))

GMM_hUSI <- function(scores){
    set.seed(123)
    if(is.null(names(scores))){
        print("hUSI scores must be a named vector!!!")
    }
    else {
        model = scores %>% Mclust()
        meanV = model$parameters$mean
        scores = sort(scores)
        senclass = c()
        meanV = sort(meanV)
       for(i in 1:length(meanV)){
            class_labels = ifelse(scores < meanV[i],i,i+1)
            names(class_labels) = names(scores)
            if(i == length(meanV)){
                senclass = c(senclass,class_labels)
            }
            else {
                senclass = c(senclass,class_labels[class_labels == i])
            }
            scores = scores[class_labels != i]
        }
        return(senclass)
    }

}
