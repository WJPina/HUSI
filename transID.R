library(R.utils)
esembl2symbol <- function(version='23',dir = '/mnt/data3/wangj2/GeneSets/',counts = NULL) {
  refpath = paste(dir,'esembl2symbol.human.v',version,'.csv',sep = '')
  if(!file.exists(refpath)){
    filename = paste('gencode.v',version,'.annotation.gtf',sep = '')
    filepath = paste(dir,filename,sep = '')
    if(!file.exists(filepath)){
      if(file.exists(paste(filepath,'.gz',sep = ''))){
        gunzip(paste(filepath,'.gz',sep = ''))
      }
      else{
        link = paste('https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_',version,'/',filename,'.gz',sep = '')
        download.file(link, destfile =paste(filepath,'.gz',sep = ''), method = "auto")
        gunzip(paste(filepath,'.gz',sep = ''))
      }
    }
    GTF =  read.csv(filepath,sep = '\t',header = F,stringsAsFactors = F,skip = 5)
    annotation = GTF[grepl('gene',GTF$V3),'V9']
    esembl2symbol = data.frame(symbol = sub(".*gene_name\\s(\\S+);.*", "\\1", annotation))
    rownames(esembl2symbol) = sub(".*gene_id\\s(\\S+);.*", "\\1", annotation)
    write.csv(esembl2symbol,file = refpath)
  }
  esembl2symbol = read.csv(refpath,header = T,row.names = 1)
  symbol = esembl2symbol[rownames(counts),'symbol']
  counts = counts[!duplicated(symbol),]
  rownames(counts) = symbol[!duplicated(symbol)]
  return(counts)
}