### define score method
scoreOCLR <- function(profile,model,m='spearman'){
  w = model$w
  genes = intersect(names(w),rownames(profile))
  w = w[genes]
  profile = profile[genes,] %>% as.matrix
  if(m == 'spearman') {husis = list(apply(profile, 2, function(s) {cor(w , s, method="spearman" )}))}
  else if (m == 'pearson'){husis = list(apply(profile, 2, function(s) {cor(w, s, method="pearson" )}))}
  else if (m == 'dot'){husis = list(apply(profile,2,function(x){(w %*% x)}))}
  # else if (m == 'logistic'){husis = list(apply(profile,2,function(x){exp(w %*% x)/(1+exp(w %*% x))}))}
  else if (m == 'logistic'){husis = list(apply(profile,2,function(x){plogis(sum(w*x))}))}
  else if (m == 'gsva'){
    n_values = c(100, 200, 500, 1000)
    husis = pbapply::pblapply(n_values, function(n) {
      top_w = sort(w, decreasing = TRUE)[1:n]
      score = GSVA::gsva(profile, list(signature=names(top_w)), method="ssgsea", ssgsea.norm = TRUE, verbose = TRUE)[1,]
      return(score)
    }, cl = 1)

    names(husis) <- as.character(n_values)
  }
  else if (m == 'mean'){
    pb <- txtProgressBar(min = 0, max = 500, style = 3)
    husis = {}
    for (n in 1:500) {
      setTxtProgressBar(pb, n)

      top_w = sort(w, decreasing = TRUE)[1:n]
      if (n == 1) {score = profile[names(top_w), ]}
      else { score = colMeans(profile[names(top_w), ])}
      husis[[as.character(n)]] <- score
    }
    close(pb)
  }
  return(husis)
}

### define normalization method
normalize_exp <- function(x){
  x_log  = log2(x + .00001) %>%  scale()
  x_centre = x_log  - (apply(x_log, 1, mean))
  return(x_centre)
}
### define class by SSE method
class_hUSI <- function(scores){
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
  print(th)
  c = ifelse(scores <= th,0,1)
  return(c)
}

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}
scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}


zerout_mat <- function(mat,rate=0.5){
  set.seed(233)
  corrdinates = expand.grid(1:nrow(mat), 1:ncol(mat))
  lost = corrdinates[sample(1:nrow(corrdinates),round(nrow(corrdinates)*rate),replace = F),]
  tmat = mat
  for(j in 1:nrow(lost)){
    tmat[lost[,1][j],lost[,2][j]] <- 0
  }
  return(tmat)
}



enrich_plot <- function(gsea,data,pathway_id=NULL,direction='up',leg_x = 0.5,leg_y = 0.5,save=T,database='Hallmarker',dir='./'){
  
  if(!is.null(pathway_id)){
    pathway_n = length(pathway_id)
    if(direction=='up'){
      cols = colorRampPalette(c("#D72E3E","#FEB1A3"))(pathway_n)[1:pathway_n]
      title = "Positively enriched gene sets"}
    else if(direction=='down'){
      cols = colorRampPalette(c("#B0C0D2", "#2E5A87"))(pathway_n)[c(1:pathway_n)]
      title = "Negatively enriched gene sets"}
    else{return(message("please input up or down !!!"))}
    
    names(cols) = pathway_id
    p <-  mygseaplot2(gsea,
                      geneSetID = pathway_id,
                      title = title,
                      color = cols,
                      base_size = 14,
                      rel_heights = c(1, 0.2, 0.4),
                      subplots = 1:3,
                      pvalue_table = FALSE,
                      leg_x = leg_x,
                      leg_y = leg_y,
                      ES_geom = "line")
    if(save){
      png(paste(dir,'valid_GSEA_',direction,'_',database,'.png',sep=''),width = 2100,height = 2000,res = 300)
      print(p)
      dev.off()
    }
    else{p}
  }
  else{message("please input right pathway names and numbers !!!")}
}


library(DelayedMatrixStats)
calc_laf <- function(mat, genes) {
  genes = genes %>% filter(Gene %in% rownames(mat))
  x = mat[genes$Gene,]
  x = t( t(x)/colMedians(mat) ) > 1
  res = apply(x, 2, function(y) {
    purrr::invoke_map_lgl(
      .f = as.list( ifelse(genes$Direction=="Increase", "isTRUE", "isFALSE") ), 
      .x = as.list(y)
    )
  })
  colSums(res)/nrow(res)
}

minmax <- function(x){
  x_scale = (x-min(x))/(max(x)-min(x))
  return(x_scale)
}

library(glmnet)
library(edgeR)
library(reticulate)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(Seurat)
library(tibble)


loadSCdata <- function(datasets=''){
tryCatch({
  objectList = list()
  if(datasets==''){datasets=c('Teo2019','Teo2019O2','Tang2019','Aarts2017','Zirkel2018','Noah2023')}
  for(d in datasets){
    if(d == 'Teo2019'){
      ### Teo2019  GSE115301 IMR90 10x OIS,3',counts, downloaded from GEO
      # Teo2019 = CreateSeuratObject(
      #   fread("Teo2019_count.csv") %>% column_to_rownames("V1") %>% data.matrix, 
      #   meta.data = fread("Teo2019_meta.csv", header = T) %>% column_to_rownames("V1")
      # ) %>% NormalizeData()
      # 
      # Teo2019$Condition = as.factor(ifelse(Teo2019$label == 'senescence','senescent','non-senescent'))
      # save(Teo2019,file='Data/Teo2019.rdata')
      
      load('Data/Teo2019.rdata')

      objectList[[d]] <- Teo2019
    }
    if(d == 'Teo2019O2'){
      ### Teo2019  GSE115301 IMR90 10x OSIS,3',counts, downloaded from GEO
      # Teo2019O2 = CreateSeuratObject(
      #   fread("Teo2019O2_count.csv") %>% column_to_rownames("V1") %>% data.matrix, 
      #   meta.data = fread("Teo2019O2_meta.csv", header = T) %>% column_to_rownames("V1")
      # ) %>% NormalizeData()
      # 
      # Teo2019O2$Condition = as.factor(ifelse(Teo2019O2$label == 'senescence','senescent','non-senescent'))
      # save(Teo2019O2,file='Data/Teo2019O2.rdata')
      
      load('Data/Teo2019O2.rdata')
      objectList[[d]] <- Teo2019O2
    }
    if(d == 'Tang2019'){
      ### Tang2019 GSE119807 HCA2 fibroblast cell RS,IRIS,3',counts, downloaded from GEO
      # Tang2019 = CreateSeuratObject(
      #   fread("Tang2019_count.csv") %>% column_to_rownames("V1") %>% data.matrix, 
      #   meta.data = fread("Tang2019_meta.csv", header = T) %>% column_to_rownames("V1")
      # ) %>% NormalizeData()
      # 
      # Tang2019$Condition = as.factor(ifelse(Tang2019$label == 'senescence','senescent','non-senescent'))
      
      # Tang2019List = list.files("Tang2019", full.names = T) %>%
      #   lapply(function(x) {
      #     scdat = fread(x) %>%
      #       column_to_rownames("GENE") %>%
      #       data.matrix }) %>% set_names( list.files("Tang2019") %>% gsub(".*_", "", .) %>% gsub("\\..*", "", .) )
      # 
      # genes = Reduce(intersect,lapply(Tang2019List, function(x)rownames(x)))
      # Tang2019 = do.call(cbind,lapply(Tang2019List,function(x) x <- x[genes,]))
      # Tang2019 = CreateSeuratObject(Tang2019) %>% NormalizeData()
      # Tang2019$group = rep(names(Tang2019List),each=400)
      # Tang2019$Condition = as.factor(ifelse(Tang2019$group %in% c("senescence","LowPD50Gy"),'senescent','non-senescent'))
      # save(Tang2019,file='Data/Tang2019.rdata')
      
      load('Data/Tang2019.rdata')
      objectList[[d]] <- Tang2019
      
    }
    if(d == 'Aarts2017'){
      ### Aarts2017 GSE94980 IMR90 OSKM-expressing reprogramming-induced senescence,3',counts, downloaded from GEO
      # Aarts2017 = CreateSeuratObject(
      #   fread("GSE94980_count.csv") %>% column_to_rownames("V1") %>% data.matrix, 
      #   meta.data = fread("GSE94980_meta.csv", header = T) %>% column_to_rownames("V1")
      # ) %>% NormalizeData()
      # 
      # Aarts2017$Condition = as.factor(ifelse(Aarts2017$label == 'senescence','senescent','non-senescent'))
      # save(Aarts2017,file='Data/Aarts2017.rdata')

      load('Data/Aarts2017.rdata')
      objectList[[d]] <- Aarts2017
    }
    if(d == 'Zirkel2018'){
      ### Zirkel2018 GSE102090 HUVEC DROP-seq RS counts, downloaded from GEO
      # Zirkel2018 = CreateSeuratObject(
      #   fread("GSE102090_count.csv") %>% column_to_rownames("V1") %>% data.matrix, 
      #   meta.data = fread("GSE102090_meta.csv", header = T) %>% column_to_rownames("V1")
      # ) %>% NormalizeData()
      # 
      # Zirkel2018$Condition = as.factor(ifelse(Zirkel2018$label == 'senescence','senescent','non-senescent'))
      # save(Zirkel2018,file='Data/Zirkel2018.rdata')
      
      load('Data/Zirkel2018.rdata')
      objectList[[d]] <- Zirkel2018
    }
    if(d == 'Noah2023'){
      ### Noah2023 GSE226225 WI38 10X RS,IR,CIS counts, downloaded from GEO
      # Noah2023list = lapply(list.dirs('GSE226225')[-1],function(x){CreateSeuratObject(Read10X(data.dir = x),project = gsub('GSE226225/GSM\\d{7}_','',x))})
      # names(Noah2023list) <- lapply(Noah2023list, function(x)x@project.name) %>% unlist()
      # Noah2023list <- lapply(Noah2023list, function(x){
      #   x$percent.mt <- PercentageFeatureSet(x, pattern = "^MT-")
      #   x = subset(x,percent.mt<12 & nCount_RNA>1000 & nCount_RNA<120000 & nFeature_RNA >300)
      #   x <- NormalizeData(x) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
      # })
      # features <- SelectIntegrationFeatures(object.list = Noah2023list)
      # Anchors <- FindIntegrationAnchors(object.list =  Noah2023list,anchor.features = features,dims = 1:30)
      # Noah2023 <- IntegrateData(anchorset = Anchors, dims = 1:30)
      # DefaultAssay(Noah2023) <- "RNA"
      # Noah2023$Condition = as.factor(ifelse(Noah2023$orig.ident == 'CTRL_2','non-senescent','senescent'))
      # save(Noah2023,file='Data/Noah2023.rdata')
      
      load('Data/Noah2023.rdata')
      objectList[[d]] <- Noah2023
    }
  }
  return(objectList)},
  error = function(e) {print("Please prepared scRNA-seq data following the functions_new.R scripts!!!")})
}


loadSetdata <- function(){
  EnrichSet<-cogena::gmt2list("Data/gene_50signatures_merge.gmt")
  EnrichSet=EnrichSet[43:49]
  names(EnrichSet)
  rep_sene_genes = readxl::read_excel("Data/SigRS.xls", sheet = 5, skip = 2)$Gene_names
  EnrichSet$SigRS = rep_sene_genes
  names(EnrichSet) 
  names(EnrichSet) <- c("SenMayo","CellAge", "GenAge", "ASIG", "SASP","AgingAtlas", "SenUp","SigRS")
  return(EnrichSet)
}

loadlafdata <- function(){
  laf_DAS = fread("Data/LaffertyWhyte2010_DAS.csv")
  laf_mSS = fread("Data/LaffertyWhyte2010_mSS.csv")
  laf_DAS_mSS =rbind(laf_DAS, laf_mSS)
  laf = list(DAS = laf_DAS,mSS=laf_mSS,DAS_mSS=laf_DAS_mSS)
  return(laf)
}


loadpredictorFix <- function(){
  CSS =  c(0.2921,0.1810,0.0524,0.1285)
  names(CSS) = c('C1QTNF6','SQOR','LYPD3','FAM83B')
  lassoCS = c(-0.158,-0.153,0.17,0.347,-0.215,0.288,-0.196,-0.245,-0.103,0.127)
  names(lassoCS) = c("ITGA8","SEMA3G","DPYSL3","IFITM1","ZNF521", "SOCS3","PCSK6", "DUSP1", "SLC44A4","IL20RB")
  ECcores36 = c(0.756, 0.739, 0.629, 0.472, 0.439, 0.429, 0.371, 0.332, 0.322, 0.3, 0.281, 0.263, 0.251, 0.235, 0.189, 0.186, 0.176, 0.14, 0.134, 0.128, 0.128, 0.102, 0.095, 0.067, 0.046, 0.042, 0.04, 0.013, -0.022, -0.078, -0.162, -0.186, -0.192, -0.215, -0.464, -0.593)
  names(ECcores36) = c('IGFBP5', 'IFI27', 'PLAT', 'MX1', 'IFIT1', 'SEMA3A', 'TAGLN', 'LAMC2', 'C19orf33', 'OAS2', 'OAS1', 'LHX6', 'CD44', 'ABI3BP', 'KCNMA1', 'P3H2', 'IFITM1', 'DIRAS3', 'PLA2G4C', 'IL32', 'ADIRF', 'ISG15', 'ADGRF5', 'SAMD9L', 'LYVE1', 'IFI6', 'SERPINE2', 'DDX58', 'NUF2', 'CD200', 'BMX', 'TOP2A', 'FST', 'PEG10', 'COLEC12', 'PHGDH')
  HCCcores19 = c(-0.24208,-0.32203,-0.15081,0.66330,0.35072,-0.46527,0.64921,-0.07426,0.37272,0.16174,0.63189,-0.61190,-0.33728,0.37964,0.45727,0.06315,0.20073,-0.16337,-0.48263)
  names(HCCcores19) = c("CDCA5","CENPF","CENPW","CDCA8","SPC25","CDKN3","CENPA","BUB1","DLGAP5","IGSF3","HMMR","TOP2A","RAD54L","TTK","GINS1","PTTG1","ETV4","GINS2","PKMYT1" )
  coefs = list('CSS'=CSS,'lassoCS'=lassoCS,'ECcores36'=ECcores36,'HCCcores19'=HCCcores19)
  return(coefs)
}


calc_scores <- function(methodslist,object=NULL,type=NULL){
  data = data.frame(cellname=colnames(object))
  if(!is.null(type)){methodslist = methodslist[type]}
  else{methodslist = methodslist}
  for(m in names(methodslist)){
    print(sprintf("Start calculating senescence scores by %s !!!",m))
    if(m == 'SenMarker'){
      scorelist=intersect(methodslist[[m]],rownames(object))
      print(sprintf("%d senescence marker was found in this dataset",length(scorelist)))
      print(scorelist)
      score = object[scorelist,] %>% GetAssayData %>% as.matrix %>% t %>% data.frame 
      rownames(score) = colnames(object)
      data = cbind(data,score)
    }
    if(m == 'SenSet'){
      scorelist=methodslist[[m]]
      for(sl in scorelist){
        print(sl)
        score = GSVA::gsva(GetAssayData(object) %>% as.matrix, list(signature=SenSetlist[[sl]]), method="ssgsea", ssgsea.norm = TRUE, verbose = TRUE)[1,]
        data[[sl]] = score
      }
    }
    if(m == 'MachineScore'){
      scorelist=methodslist[[m]]
      for(sl in scorelist){
        print(sl)
        if(sl %in% c('CSS','lassoCS','ECcores36','HCCcores19')){
          if(length(names(coefslist[[sl]])[names(coefslist[[sl]]) %in%rownames(object)])==0){print(sprintf('%s: no gene was found!!!',sl))}
          else{
            mat = object[names(coefslist[[sl]]),] %>% GetAssayData %>% data.matrix 
            mat = apply(mat,2,function(x) {x * coefslist[[sl]][rownames(mat)]}) 
            tryCatch({score = colSums(mat)}, 
            error = function(e) { score = mat})
            data[[sl]] = score
            }
          }
        if(sl == 'SENCAN'){
          score = cal_SENCAN(object)
          data[[sl]] = score
        }
        if(sl == 'SENCID'){
          score = cal_SENCID(object)
          data[[sl]] = score
        }
      }
    }
    
    if(m == 'TraditionScore'){
      scorelist=methodslist[[m]]
      for(sl in scorelist){
        print(sl)
        if(sl %in% c('DAS','mSS','DAS_mSS')){
          score = calc_laf(GetAssayData(object) %>% as.matrix, laflist[[sl]])
          data[[sl]] = score
        }
        if(sl == 'CS_score'){
          scorecs = GSVA::gsva(GetAssayData(object) %>% as.matrix, CSsccorelist[c('Underexpressed','Overexpressed')], method="ssgsea", ssgsea.norm = TRUE,  verbose = TRUE)
          data[[sl]] = scorecs[2,]-scorecs[1,]
        }
      }
    }
  }
  return(data)
}

cal_SENCAN <- function(object){
  dm <- object@assays$RNA@counts %>% data.matrix
  genes = rownames(object)
  gene_ref = dge$genes$gene_name
  count_ref = dge$counts
  # Map gene symbols of the input file to the classifier structure
  map <- match(gene_ref, genes)
  # Impute missing values if necessary
  dmMapped <- dm[map,]
  dmImputed <- dmMapped
  for (i in 1:ncol(dmImputed)) {
    sample <- dmImputed[,i]
    if (any(is.na(sample))) {
      cors <- rep(NA, ncol(count_ref))
      names(cors) <- rownames(dge$samples)
      for (j in 1:ncol(count_ref)) {
        have_val <- which(!is.na(sample))
        ### calculated Kendall's tau rank correlation coefficient between input sample and dge sample among shared genes
        cors[j] <- pcaPP::cor.fk(sample[have_val], count_ref[have_val,j])
      }
      ### for missing genes in input sample, give mean gene values from samples in dge with top 2 coefficients
      nearest <- tail(order(cors), n=2)
      need_val <- which(is.na(sample))
      dmImputed[need_val,i] <- apply(count_ref[need_val,nearest], 1, mean)
    }
  }
  # Normalize the input sample towards the classifier reference
  dge_pred <- edgeR::DGEList(counts = dmImputed, genes=data.frame(ensembl_id=dge$genes$ensembl_id))
  dge_pred <- calcNormFactors_WSPref(dge_pred, edger_ref_column)
  dmPredCPM <- edgeR::cpm(dge_pred, log=T)
  dmPredCPM[dmPredCPM < min(dmCPM)] <- min(dmCPM)
  pred_standardized_expression <- (dmPredCPM - gene_mean) / gene_sd
  
  # Calculate classifier scores
  res <- predict(full_classifier, t(pred_standardized_expression[use_gene_ix,]), s="lambda.1se")[,1]
  sencan_score <- 1/(1+exp(-res))
  return(sencan_score)
}

cal_SENCID <- function(object,denoise=FALSE){
  tryCatch({
      SenCIDapi <- import("SenCID.api")
  sc <- import("scanpy")
  np <- import("numpy")
  counts = object@assays$RNA@counts %>% as.matrix 
  adata = sc$AnnData(t(counts))
  adata$var_names = rownames(counts)
  adata$obs_names = colnames(counts)

  res = SenCIDapi$SenCID(adata = adata, sidnums = as.integer(c(1,2,3,4,5,6)), denoising = denoise, binarize = TRUE, threads = as.integer(1), savetmp = TRUE)
  recSID = res[[2]]
  SID = gsub('rec_','',which.max(colMeans(recSID[,-1])) %>% names())
  print(sprintf('%s was used in this dataset',SID))
  pred_dict = res[[1]][[SID]]
  s = pred_dict$SID_Score

  return(s)
  },
  error = function(e){print("Please install SenCID in a appointed python environment!!!")})
}


calc_times <- function(object){
  times = c()
  
  startTime <- Sys.time()
  calc_scores(list('SenSet'=c("SenUp")),object)
  endTime <- Sys.time()
  t = as.numeric(endTime-startTime)
  times <- c(times,t)
  
  startTime <- Sys.time()
  calc_scores(list('MachineScore'=c('ECcores36')),object)
  endTime <- Sys.time()
  t = as.numeric(endTime-startTime)
  times <- c(times,t)
  
  startTime <- Sys.time()
  calc_scores(list('MachineScore'=c('SENCID')),object)
  endTime <- Sys.time()
  t = as.numeric(endTime-startTime)
  times <- c(times,t)
  
  startTime <- Sys.time()
  calc_scores(list('TraditionScore'=c('CS_score')),object)
  endTime <- Sys.time()
  t = as.numeric(endTime-startTime)
  times <- c(times,t)
  
  startTime <- Sys.time()
  scoreOCLR(profile = GetAssayData(object) %>% as.matrix,model = SenOCLR,m = 'logistic')
  endTime <- Sys.time()
  t = as.numeric(endTime-startTime)
  times <- c(times,t)
  
  names(times) <- c("SenUp",'ECcores36','SENCID','CS_score','hUSI')
  
  return(times)
}





