# Annotating GWAS genes present in credible regions

# 21/02/2017
#################### getting genes from coordinates


library(biomaRt)
library(dplyr)
currentDate <- Sys.Date() # to save date in name of output files


# T2D
regions= read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set/filenames_credible_regions_Feb17.txt",header = F)

#FG
regions= read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_FG/FG_regions.txt",header = F)

regions$V1=gsub("\\chr*","",regions$V1) #remove "chr" in every element in location

region_names=sub(".*_","",regions$V1) # remove everything before last "_"
region_names=sub("\\..*","",region_names) # remove everything after dot
regions$V1=gsub("\\_[^_]*$","",regions$V1) # remove everything after last "_"
regions$V1=gsub("_", ":", regions$V1)  # replace _ for :
regions$V1=gsub("-", ":", regions$V1)  # replace - for :

regions=cbind(regions,region_names)
colnames(regions)[1]=c("location")
rm(region_names)

# chr.region=do.call(paste, c(regions, sep=":"))  # convert each row (three columns) to element in list for biomaRt call
# filterlist <- list(chr.region) 
filterlist <- list(regions$location)  

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") 
# grch37 used for annotation of RNA-seq background

results <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name', "start_position", "end_position"),
                 filters = c("chromosomal_region"),values = filterlist, mart = ensembl)


results=results[which(results$gene_biotype=="lincRNA" | results$gene_biotype=="protein_coding" ),]

# ############ merge genes and GWAS loci they come from


  regions$chr= as.numeric(gsub("(.*?):.*", "\\1",regions$location))
  minus_chr= gsub("^[^:]*:","",regions$location)             # remove everything up to 1st colon
  regions$start_position =  as.numeric(gsub("(.*?):.*", "\\1",minus_chr))    # remove everything after colon
  rm(minus_chr)
  regions$end_position =  as.numeric(gsub(".*:","",regions$location))    # remove everything before last colon
############3
  ###########
  # need to sort by order of first and then second digit for IRanges to work properly!!
  # check previous GWAS loci annotation
  


  # creating IRanges objects for easy comparison
 #  kbs_ranges=split(IRanges(regions$start_position,
 #                           regions$end_position),
 #                   regions$chr) # range for kbs
 # 
 #  GWAS_ranges=split(IRanges(results$start_position,
 #                            results$end_position,
 #                            names = results$external_gene_name),
 #                    results$chromosome_name) # range for GWAS
 # merged=data.frame()
 #   for(c in 1:20){
 #    overlap=findOverlaps(GWAS_ranges[[c]],kbs_ranges[[c]],minoverlap = 1L) # for every chromosome, relates indices in GWAS_ranges with indices in kbs_ranges lists
 #    check_query = GWAS_ranges[[c]][overlap@from]@NAMES # genes in location names from overlap indices
 #    check_subject = regions[match(kbs_ranges[[c]][overlap@to]@start,regions$start_position),]#gets positions in kbs table from start position in subject (got from overlap matches)
 #    check_subject$genes_in_loci=check_query #merges genes in loci (query table) with GWAS_loci table
 #    merged <- rbind(merged,check_subject)
 #  }
#   merged[[d]]$ensembl_gene_id <- GWAS_list[[d]][match(merged[[d]]$genes_in_loci,GWAS_list[[d]]$external_gene_name),1] ###append gene id of every gene
#   
#   write.table(merged[[d]],quote=F,row.names=F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/",
#                                                             currentDate,"_genes_and_lincRNA",d,"kbs_around_GWAS_loci_glycemic_and_T2D.tsv",sep=""),sep="\t")
#   
# }

# write.table(results,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set/",currentDate,"annotated_genes_in_credible_regions.txt",sep=""),sep="\t",row.names = F,col.names=T,quote = F)

# add different lengths

kbs <- list()

distance = c("50","100","200","500") # list of distances in kb

for(d in distance){
  
  
  chr=as.numeric(gsub("\\:[^:]*","",regions[,1]))
  minus_d =  as.numeric(gsub("(^.+:)(\\d+)(:.+$)", "\\2",regions[,1]))-as.numeric(d)*1000   #vector of -distance
  plus_d =  as.numeric(gsub(".*:","",regions[,1]))+as.numeric(d)*1000          #vector of + distance
  
  
  kbs[[d]]$location = paste(chr,minus_d,plus_d,sep=":") #paste everything
  
  # call ensembl
  filterlist <- list(kbs[[d]]$location)  
  
  ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") 
  # grch37 used for annotation of RNA-seq background
  
  results <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name', "start_position", "end_position"),
                   filters = c("chromosomal_region"),values = filterlist, mart = ensembl)
  
  
  results=results[which(results$gene_biotype=="lincRNA" | results$gene_biotype=="protein_coding" ),]
  
  
  write.table(results,quote=F,row.names=F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_FG/",
                                                      currentDate,"FG_annotated_genes_in_credible_regions_plusminus_",d,"_kb.txt",sep=""),sep="\t")
  rm(minus_d,plus_d,chr,filterlist,results)
  }

