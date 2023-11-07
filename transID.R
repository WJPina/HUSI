library(R.utils)
library(data.table)
esembl2symbol <- function(counts = NULL,database = 'genecode',version='70',dir = '/mnt/data3/wangj2/GeneSets/',GRCh ='37') {
  if(database == 'genecode'){
    refpath = paste(dir,'genecode.esembl2symbol.human.v',version,'.csv',sep = '')
    filename = paste('gencode.v',version,'.annotation.gtf.gz',sep = '')
    link = paste('https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_',version,'/',filename,'.gz',sep = '')
  }
  else if (database == 'esemble') {
    refpath = paste(dir,'esemble.esembl2symbol.human.v',version,'.csv',sep = '')
    filename = paste('Homo_sapiens.GRCh',GRCh,'.',version,'.gtf.gz',sep = '')
    link = paste('https://ftp.ensembl.org/pub/release-',version,'/gtf/homo_sapiens/Homo_sapiens.GRCh',GRCh,'.',version,'.gtf.gz',sep = '')
  }
  if(!file.exists(refpath)){
    filepath = paste(dir,filename,sep = '')
    if(!file.exists(filepath) & !file.exists(gsub('.gz','',filepath))){
        download.file(link, destfile =filepath, method = "auto")
    }
    GTF =  fread(filepath)
    if (database == 'genecode'){
      annotation = GTF[grepl('gene',GTF$V3),'V9'][[1]]    
      esembl2symbol = data.frame(symbol = sub(".*gene_name\\s(\\S+);.*", "\\1", annotation))
      rownames(esembl2symbol) = sub(".*gene_id\\s(\\S+);.*", "\\1", annotation)
    }
    else if (database == 'esemble') {
      annotation = GTF$V9
      genes = sub('.*gene_name\\s"(\\S+)";.*', "\\1", annotation)
      esembl2symbol = data.frame(symbol=genes[!duplicated(genes)])
      rownames(esembl2symbol) = sub('.*gene_id\\s"(\\S+)";.*', "\\1", annotation[!duplicated(genes)])
    }
    write.csv(esembl2symbol,file = refpath)
  }
  esembl2symbol = read.csv(refpath,header = T,row.names = 1)
  counts = counts[rownames(counts) %in% rownames(esembl2symbol),]
  symbol = esembl2symbol[rownames(counts),'symbol']
  counts = counts[!duplicated(symbol),]
  rownames(counts) = symbol[!duplicated(symbol)]
  return(counts)
}
