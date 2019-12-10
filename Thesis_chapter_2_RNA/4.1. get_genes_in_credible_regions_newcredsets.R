# Annotating GWAS genes present in credible regions

# 21/02/2017
#################### getting genes from coordinates


library(biomaRt)
library(dplyr)
currentDate <- Sys.Date() # to save date in name of output files


# T2D
#FG in old version

setwd("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis")

credset_names = read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/credset/credset_names.txt",
                           header = F)

credset = list()
region = list()

for (n in credset_names$V1) {
  credset[[n]] = read.table(
    file = paste(
      "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/credset/",
      n,
      sep = ""
    ),
    header = T
  ) # 380 99% credible sets T2D 2017 HRC (Mahajan et al 2018)

  if(nrow(credset[[n]])==1){
    region[[n]]$location = paste(unique(credset[[n]]$Chr),min(credset[[n]]$Pos),min(credset[[n]]$Pos)+1,sep=":")  }
  else{
    region[[n]]$location = paste(unique(credset[[n]]$Chr),min(credset[[n]]$Pos),max(credset[[n]]$Pos),sep=":")
    
  }
}

region = unlist(region)


filterlist <- list(region)  

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") 
# grch37 used for annotation of RNA-seq background

results <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name', "start_position", "end_position"),
                 filters = c("chromosomal_region"),values = filterlist, mart = ensembl)

results=results[which(results$gene_biotype=="lincRNA" | results$gene_biotype=="protein_coding" ),]

nrow(results)

# ############ merge genes and GWAS loci they come from


# 
#   chr= as.numeric(gsub("(.*?):.*", "\\1",region))
#   minus_chr= gsub("^[^:]*:","",region)             # remove everything up to 1st colon
#   start_position =  as.numeric(gsub("(.*?):.*", "\\1",minus_chr))    # remove everything after colon
#   rm(minus_chr)
#   end_position =  as.numeric(gsub(".*:","",region))    # remove everything before last colon
#   regions = cbind(chr,start_position,end_position)
# merge with genomicRanges


# add different lengths

kbs <- list()

distance = c("0","50","100","200","500") # list of distances in kb

for(d in distance){
  
  
  chr=as.numeric(gsub("\\:[^:]*","",region))
  minus_d =  as.numeric(gsub("(^.+:)(\\d+)(:.+$)", "\\2",region))-as.numeric(d)*1000   #vector of -distance
  plus_d =  as.numeric(gsub(".*:","",region))+as.numeric(d)*1000          #vector of + distance
  
  
  kbs[[d]]$location = paste(chr,minus_d,plus_d,sep=":") #paste everything
  
  # call ensembl
  filterlist <- list(kbs[[d]]$location)  
  
  ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") 
  # grch37 used for annotation of RNA-seq background
  
  results <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name', "start_position", "end_position"),
                   filters = c("chromosomal_region"),values = filterlist, mart = ensembl)
  
  
  results=results[which(results$gene_biotype=="lincRNA" | results$gene_biotype=="protein_coding" ),]
  
  
  write.table(results,quote=F,row.names=F,file=paste(currentDate,"T2D_annotated_genes_in_credible_regions_plusminus_",d,"_kb.txt",sep=""),sep="\t")
  rm(minus_d,plus_d,chr,filterlist,results)
  }

## just DIMAS and WOOD beta cell physiological regions

# T2D

phys  = c("ADCY5","ARAP1","CDKAL1","CDKN2A","DGKB","GCK","HHEX","IGF2BP2","KCNQ1","MTNR1B","PROX1","SLC30A8","TCF7L2","THADA")
setwd("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis")

credset_names = read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/credset/credset_names.txt",
                           header = F)

credset = list()
region = list()

for (n in credset_names$V1) {
  credset[[n]] = read.table(
    file = paste(
      "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/credset/",
      n,
      sep = ""
    ),
    header = T
  ) # 380 99% credible sets T2D 2017 HRC (Mahajan et al 2018)
  if (grepl(paste(phys, collapse="|"), n)) {
    if (nrow(credset[[n]]) == 1) {
      region[[n]]$location = paste(unique(credset[[n]]$Chr),
                                   min(credset[[n]]$Pos),
                                   min(credset[[n]]$Pos) + 1,
                                   sep = ":")
    }
    else{
      region[[n]]$location = paste(unique(credset[[n]]$Chr),
                                   min(credset[[n]]$Pos),
                                   max(credset[[n]]$Pos),
                                   sep = ":")
      
    }
  }
}

region = unlist(region)


filterlist <- list(region)  

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") 
# grch37 used for annotation of RNA-seq background

results <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name', "start_position", "end_position"),
                 filters = c("chromosomal_region"),values = filterlist, mart = ensembl)

results=results[which(results$gene_biotype=="lincRNA" | results$gene_biotype=="protein_coding" ),]

nrow(results)

# ############ merge genes and GWAS loci they come from


# 
#   chr= as.numeric(gsub("(.*?):.*", "\\1",region))
#   minus_chr= gsub("^[^:]*:","",region)             # remove everything up to 1st colon
#   start_position =  as.numeric(gsub("(.*?):.*", "\\1",minus_chr))    # remove everything after colon
#   rm(minus_chr)
#   end_position =  as.numeric(gsub(".*:","",region))    # remove everything before last colon
#   regions = cbind(chr,start_position,end_position)
# merge with genomicRanges


# add different lengths

kbs <- list()

distance = c("0","50","100","200","500") # list of distances in kb

for(d in distance){
  
  
  chr=as.numeric(gsub("\\:[^:]*","",region))
  minus_d =  as.numeric(gsub("(^.+:)(\\d+)(:.+$)", "\\2",region))-as.numeric(d)*1000   #vector of -distance
  plus_d =  as.numeric(gsub(".*:","",region))+as.numeric(d)*1000          #vector of + distance
  
  
  kbs[[d]]$location = paste(chr,minus_d,plus_d,sep=":") #paste everything
  
  # call ensembl
  filterlist <- list(kbs[[d]]$location)  
  
  ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") 
  # grch37 used for annotation of RNA-seq background
  
  results <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name', "start_position", "end_position"),
                   filters = c("chromosomal_region"),values = filterlist, mart = ensembl)
  
  
  results=results[which(results$gene_biotype=="lincRNA" | results$gene_biotype=="protein_coding" ),]
  
  
  write.table(results,quote=F,row.names=F,file=paste(currentDate,"_physlociDimasWood_T2D_annotated_genes_in_credible_regions_plusminus_",d,"_kb.txt",sep=""),sep="\t")
  rm(minus_d,plus_d,chr,filterlist,results)
}
