
setwd("/mnt/data1/wangj/AgingScore/GSE163530_COVID-19/GSE171668_scnRNA/")
###### load data from multior aganism 
Bulk = read.csv("GSE171668_rsem.genes.counts.matrix.txt", header = T, row.names = 1,sep = '\t')
meta = read.csv("GSE171668_bulk_metadata.csv", header = T, row.names = 1)
