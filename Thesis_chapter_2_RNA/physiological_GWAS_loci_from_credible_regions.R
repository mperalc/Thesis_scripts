# from physiological loci classifications (Dimas 2014 and Wood 2017, provided by Jason) 
# make lists of genes at variable windows from credible regions (DIAGRAM)

######################### get genes in credible sets and divide according to physiological loci classification #########

library(biomaRt)
library(dplyr)
currentDate <- Sys.Date() # to save date in name of output files

extendLocations=function(df){
  
  df$chr= as.numeric(gsub("(.*?):.*", "\\1",df$location))
  minus_chr= gsub("^[^:]*:","",df$location)             # remove everything up to 1st colon
  df$start_position = as.numeric(gsub("(.*?):.*", "\\1",minus_chr))    # remove everything after colon
  rm(minus_chr)
  df$end_position = as.numeric(gsub(".*:","",df$location))    # remove everything before last colon
  
  
  # add different lengths
  
  kbs <- list()
  
  distance = c("50","100","200","500") # list of distances in kb
  
  for(d in distance){
    
    
    chr=as.numeric(gsub("\\:[^:]*","",df[,1]))
    minus_d =  as.numeric(gsub("(^.+:)(\\d+)(:.+$)", "\\2",df[,1]))-as.numeric(d)*1000   #vector of -distance
    plus_d =  as.numeric(gsub(".*:","",df[,1]))+as.numeric(d)*1000          #vector of + distance
    
    
    kbs[[d]]$location = paste(chr,minus_d,plus_d,sep=":") #paste everything
    
    # call ensembl
    results = callBiomart(kbs[[d]])
    write.table(results,quote=F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                                           currentDate,
                                           "_T2D_phys_loci_",
                                           d,
                                           "_kb_from_credset_",
                                           l,
                                           ".txt",sep=""),
                col.names = T,row.names = F,sep="\t")
    
  }
  
  
}
callBiomart=function(df){
  filterlist <- list(df$location)  
  
  ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") 
  # grch37 used for annotation of RNA-seq background
  
  results <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name', "start_position", "end_position"),
                   filters = c("chromosomal_region"),values = filterlist, mart = ensembl)
  
  
  results=results[which(results$gene_biotype=="lincRNA" | results$gene_biotype=="protein_coding" ),]
  return(results)
}


regions= read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_T2D/filenames_credible_regions_Feb17.txt",header = F)

physLoci= read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/Physiological-Loci.txt",
                     header = T)

lociTypes=c("HG","IR","PI","UC","BC","Wood") # all types from Dimas

# making consistent some loci names
physLoci$Locus=as.character(physLoci$Locus)

physLoci[which(physLoci$Locus=="CDKN2A_B"),1] <- "CDKN2A-B"
physLoci[which(physLoci$Locus=="CDC123_CAMK1D"),1] <- "CDC123-CAMK1D"
physLoci[which(physLoci$Locus=="HHEX_IDE"),1] <- "HHEX-IDE"
physLoci[grep("KCNQ1",physLoci$Locus),1] <- "KCNQ1"
# keep in mind there are two KCNQ1
physLoci[which(physLoci$Locus=="TSPAN8_LGR5"),1] <- "TSPAN8-LGR5"


# trim loci names in GWAS regions
regions$V1=gsub("\\chr*","",regions$V1) #remove "chr" in every element in location

region_names=sub(".*_","",regions$V1) # remove everything before last "_"
region_names=sub("\\..*","",region_names) # remove everything after dot
regions$V1=gsub("\\_[^_]*$","",regions$V1) # remove everything after last "_"
regions$V1=gsub("_", ":", regions$V1)  # replace _ for :
regions$V1=gsub("-", ":", regions$V1)  # replace - for :

regions=cbind(regions,region_names)
colnames(regions)[1]=c("location")
rm(region_names)
regions$region_names=as.character(regions$region_names)


lociList= list()

for(l in lociTypes){
  print(l)
  if(l=="UC"){
    physLociList=physLoci[which(physLoci$Dimas2014=="UC"),"Locus"] # UC 
    # those in Wood are considered BC loci
    df=regions[grep(paste(physLociList,collapse="|"),regions$region_names),] # take credible regions
    
    lociList[[l]] = callBiomart(df)
    write.table(lociList[[l]],file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                                         currentDate,
                                         "_T2D_phys_loci_0_kb_from_credset_",
                                         l,
                                         ".txt",sep=""),
                col.names = T,row.names = F,quote=F,sep="\t")
    extendLocations(df)

    
  }
  if(l=="BC"){
    physLociList=physLoci[which(physLoci$Dimas2014=="BC"),
                          "Locus"] # BC in Dimas o
    # those in Wood are considered BC loci
    df=regions[grep(paste(physLociList,collapse="|"),regions$region_names),] # take credible regions
    
    lociList[[l]] = callBiomart(df)
    write.table(lociList[[l]],file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                                         currentDate,
                                         "_T2D_phys_loci_0_kb_from_credset_",
                                         l,
                                         ".txt",sep=""),
                col.names = T,row.names = F,quote=F,sep="\t")
    extendLocations(df)
  }
  if(l=="IR"){
    physLociList=physLoci[which(physLoci$Dimas2014=="IR"),
                          "Locus"] # IR in Dimas 
    df=regions[grep(paste(physLociList,collapse="|"),regions$region_names),] # take credible regions
    
    lociList[[l]] = callBiomart(df)
    write.table(lociList[[l]],file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                                         currentDate,
                                         "_T2D_phys_loci_0_kb_from_credset_",
                                         l,
                                         ".txt",sep=""),
                col.names = T,row.names = F,quote=F,sep="\t")
    extendLocations(df)
  }
  if(l=="HG"){
    physLociList=physLoci[which(physLoci$Dimas2014=="HG"),
                          "Locus"] # HG in Dimas 
    df=regions[grep(paste(physLociList,collapse="|"),regions$region_names),] # take credible regions
    
    lociList[[l]] = callBiomart(df)
    write.table(lociList[[l]],file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                                         currentDate,
                                         "_T2D_phys_loci_0_kb_from_credset_",
                                         l,
                                         ".txt",sep=""),
                col.names = T,row.names = F,quote=F,sep="\t")
    extendLocations(df)
  }
  if(l=="PI"){
    
    physLociList=physLoci[which(physLoci$Dimas2014=="PI"),
                          "Locus"] #PI in Dimas ("IP")
  
    df=regions[grep(paste(physLociList,collapse="|"),regions$region_names),] # take credible regions
    
    lociList[[l]] = callBiomart(df)
    write.table(lociList[[l]],file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                                         currentDate,
                                         "_T2D_phys_loci_0_kb_from_credset_",
                                         l,
                                         ".txt",sep=""),
                col.names = T,row.names = F,quote=F,sep="\t")
    extendLocations(df)
  }
  if(l=="Wood"){
    
    physLociList=physLoci[which(!is.na(physLoci$Wood2017.GenomeWide)),
                          "Locus"] #Detected genome wide in Wood, for any measure of BC function
    # gene names can be shared by other Dimas classifications
    
    df=regions[grep(paste(physLociList,collapse="|"),regions$region_names),] # take credible regions
    
    lociList[[l]] = callBiomart(df)
    write.table(lociList[[l]],file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                                         currentDate,
                                         "_T2D_phys_loci_0_kb_from_credset_",
                                         l,
                                         ".txt",sep=""),
                col.names = T,row.names = F,quote=F,sep="\t")
    extendLocations(df)
  }
}


######################### hypergeometric test #######################

## Hypergeometric test for datasets: are they enriched in T2D/FG genes?

# Takes in:
# A list of SNPs in LD with credible sets for DIAGRAM (T2D)
# I then need to take the SNPs coordinates of the extremes of each regions, and take out all genes inside (ENSEMBL)
# Datasets to test for enrichment in the above
# The background dataset to create random sets of genes. It would be the list of genes from the initial RNA-seq list (before filtering)

# Then, I test how my datasets are enriched in T2D genes compared to random sets of genes of the same size. 

library(Homo.sapiens)
library(dplyr)
library(biomaRt)
library(ggplot2)


phyperRandom <- function(myGeneList, myGeneSet, genome){
  myRandomGS <- sample( genome,size=length(myGeneSet) )  # take random sample of genes from whole genome, same size as set tested
  myX <- length(which(myGeneList %in% myRandomGS))  # number of T2D in random set
  myM <- length(which(myGeneList %in% genome)) #  T2D genes in my pool
  myN <- length(genome) - length(which(myGeneList %in% genome))  # genes in genome not T2D 
  myK <- length(myGeneSet) # size of set tested
  return(list(prob=1-phyper(q=myX-1,
                            m=myM, n=myN, k=myK), n.overlap=myX))
}


currentDate <- Sys.Date() # to save date in name of output files

stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN6", "EN7")

trait= "T2D"    # T2D or fasting glucose (FG)?
### files

#### background "genome" data
# just data used for the Differential expression analysis

load("/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/session_objects/dge_cc_no_monogenic_0kb.xz")  # 15220 genes and lincRNA

# gene test data


top="623"  # 500, 623 or all

# files with DE results (across-stages)
# # all genes, get list of ensembl gene ids

sig_stages=list()
for(s in stage){
  
  sig_stages[[s]] <- read.csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/across-stages/2017-07-03_results_without_monogenic/2017-07-03_sig_",s,"_diff_expression_maxvals_across-stages_results_logFC1.csv",sep=""),header = T) # read in DEA data for each stage
  colnames(sig_stages[[s]])[1]="GeneID"
  # sig_stages[[s]]=sig_stages[[s]][ order(sig_stages[[s]]$adj.P.Val,decreasing=F),] # order by adjusted p.value 
  sig_stages[[s]]=sig_stages[[s]][ order(sig_stages[[s]]$logFC,decreasing=T),] # order by logFC 
  
  if(top==500){# top 500 genes
  sig_stages[[s]]=sig_stages[[s]][c(1:500),]
  }
  if(top==623){# top 623 genes (min significant DE genes, in PGT stage)
  sig_stages[[s]]=sig_stages[[s]][c(1:623),]
  }
}


lociTypes=c("HG", "IR", "UC", "BC","PI","Wood")
distance = c("0","50","100","200","500") # list of distances in kb


# ####individually, for each phenotypical loci Type ##########

test <- list()


for(d in distance){
  for(l in lociTypes){
    test[[d]]=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                                    "2017-06-06",
                                    "_T2D_phys_loci_",d,"_kb_from_credset_",
                                    l,  # change when needed
                                    ".txt",sep=""),
                         header=T)
    colnames(test[[d]])[1]=c("GeneID")
    test[[d]]$GeneID=as.character(test[[d]]$GeneID)
    
    df_list <- list()
    
    overlap_list <- list()  # to write down number of T2D genes in diff expr results
    
    prob_results=list()
    sums_probs=list()
    b_test <- list()
    overlap <- list()
    
    for(z in stage){
      
      overlap[["total_in_background"]] = length(which(test[[d]]$GeneID %in% dge_cc$genes$ensembl_gene_id))
      
      overlap[[z]] = length(which(test[[d]]$GeneID %in% sig_stages[[z]]$GeneID))
      
      if(overlap[[z]]<3){  # Don't want to calculate probability if the overlap is extremely small. That significance would be inflated.
        prob_results[[z]]=1
      }
      else{
      # p-value for enrichment in T2D genes using my differentially expressed genes for each stage
      prob_results[[z]]=1-phyper(q=length(which(test[[d]]$GeneID %in% sig_stages[[z]]$GeneID))-1,
                                 m=length(which(test[[d]]$GeneID %in% dge_cc$genes$ensembl_gene_id)),
                                 n=length(dge_cc$genes$ensembl_gene_id) - length(which(test[[d]]$GeneID %in% dge_cc$genes$ensembl_gene_id)),
                                 k=length(sig_stages[[z]]$GeneID))
      }
      
      gene_list=test[[d]][which(test[[d]]$GeneID %in% sig_stages[[z]]$GeneID),]
      
      # same, for random set of genes of same size. random gene set distribution of p-values
      random=list()
      random_overlaps=numeric()
      probab_random=numeric()
      
      for(i in 1:10000){
        random[[i]] <- as.data.frame(phyperRandom( test[[d]]$GeneID, sig_stages[[z]]$GeneID, dge_cc$genes$ensembl_gene_id))  # get probs from random gene sets, for each diff expr set tested
        probab_random[i] <- random[[i]]$prob
        random_overlaps[i] <- random[[i]]$n.overlap
        
      }
      
      
      hist(random_overlaps)
      
      # empirical p-value 
      
      sums_probs[[z]]$pval=((sum(probab_random<=prob_results[[z]])+1)/(10000+1)) # fraction of probabilities from random draws that are more extreme (smaller) than the one from my tested set
      # +1 in case sum(probab_random<=prob_results[[z]]) ==0
      # that would be the permuted? p-value
      
      # calculate the 95% confidence interval from the binomial distribution
      b_test[[z]]=binom.test(sum(probab_random<=prob_results[[z]]),10000)  # If those fractions are my "successes", then the binomial test can give me a confidence interval
      sums_probs[[z]]$conf.int.low=b_test[[z]]$conf.int[1]
      sums_probs[[z]]$conf.int.up=b_test[[z]]$conf.int[2]
      
      
      print(z)
      
    }
    overlap_list[[d]] <- do.call("rbind", overlap)
    
    df_list[[d]] <- do.call("rbind", sums_probs)
    
    # ######### save tables
    # ordered by p-value
    
    
    # write.table(df_list[[d]],quote=F,row.names=T,col.names = T,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
    #                                                                      "top",top,"_DEA_timecourse_pvals_hyperg_",
    #                                                                      currentDate,
    #                                                                       "_genes_enrichment_10k_permutations_",
    #                                                                       d,
    #                                                                       "_kb_from_credible_regions_T2D_",
    #                                                                       l,
    #                                                                       ".txt",sep=""),sep="\t")
    # write.table(overlap_list[[d]],quote=F,row.names=T,col.names = F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
    #                                                                            "top",top,"_DEA_timecourse_overlaps_hyperg_",
    #                                                                            currentDate,
    #                                                                            "_genes_enrichment_10k_permutations_",
    #                                                                            d,
    #                                                                            "_kb_from_credible_regions_T2D_",
    #                                                                            l,
    #                                                                            ".txt",sep=""),sep="\t")
    # ordered by logFC
    
    
    write.table(df_list[[d]],quote=F,row.names=T,col.names = T,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                                                                          "top",top,"_DEA_timecourse_pvals_hyperg_",
                                                                          currentDate,
                                                                          "_genes_enrichment_10k_permutations_",
                                                                          d,
                                                                          "_kb_from_credible_regions_T2D_",
                                                                          l,
                                                                          "_orderedlogFC.txt",sep=""),sep="\t")
    write.table(overlap_list[[d]],quote=F,row.names=T,col.names = F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                                                                               "top",top,"_DEA_timecourse_overlaps_hyperg_",
                                                                               currentDate,
                                                                               "_genes_enrichment_10k_permutations_",
                                                                               d,
                                                                               "_kb_from_credible_regions_T2D_",
                                                                               l,
                                                                               "_orderedlogFC.txt",sep=""),sep="\t")
  }
}




############## joint beta cell function #############################

phyperRandom <- function(myGeneList, myGeneSet, genome){
  myRandomGS <- sample( genome,size=length(myGeneSet) )  # take random sample of genes from whole genome, same size as set tested
  myX <- length(which(myGeneList %in% myRandomGS))  # number of T2D in random set
  myM <- length(which(myGeneList %in% genome)) #  T2D genes in my pool
  myN <- length(genome) - length(which(myGeneList %in% genome))  # genes in genome not T2D 
  myK <- length(myGeneSet) # size of set tested
  return(list(prob=1-phyper(q=myX-1,
                            m=myM, n=myN, k=myK), n.overlap=myX))
}


currentDate <- Sys.Date() # to save date in name of output files

stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN6", "EN7")

trait= "T2D"    # T2D or fasting glucose (FG)?
### files

#### background "genome" data
# just data used for the Differential expression analysis

load("/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/session_objects/dge_cc_no_monogenic_0kb.xz")  # 15220 genes and lincRNA
dge_cc=x
# gene test data


top="all"  # 500, 623 or all

# files with DE results (across-stages)
# # all genes, get list of ensembl gene ids

sig_stages=list()
for(s in stage){
  
  sig_stages[[s]] <- read.csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/across-stages/2017-07-03_results_without_monogenic_0kb/2017-07-03_sig_",s,"_diff_expression_maxvals_across-stages_results_logFC1.csv",sep=""),header = T) # read in DEA data for each stage
  colnames(sig_stages[[s]])[1]="GeneID"
  # sig_stages[[s]]=sig_stages[[s]][ order(sig_stages[[s]]$adj.P.Val,decreasing=F),] # order by adjusted p.value 
  sig_stages[[s]]=sig_stages[[s]][ order(sig_stages[[s]]$logFC,decreasing=T),] # order by logFC 
  
  if(top==500){# top 500 genes
    sig_stages[[s]]=sig_stages[[s]][c(1:500),]
  }
  if(top==623){# top 623 genes (min significant DE genes, in PGT stage)
    sig_stages[[s]]=sig_stages[[s]][c(1:623),]
  }
}


lociTypes=c("HG", "IR", "UC", "BC","PI","Wood")
distance = c("0","50","100","200","500") # list of distances in kb


# ####individually, for each phenotypical loci Type ##########

test <- list()

# BC + HG + PI + Wood (unique, because Wood has elements of the other classifiers)
for(d in distance){
  BC=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                           "2017-06-06",
                           "_T2D_phys_loci_",d,"_kb_from_credset_",
                           "BC",  # change when needed
                           ".txt",sep=""),  header=T)
  HG=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                           "2017-06-06",
                           "_T2D_phys_loci_",d,"_kb_from_credset_",
                           "HG",  # change when needed
                           ".txt",sep=""),  header=T)
  PI=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                           "2017-06-06",
                           "_T2D_phys_loci_",d,"_kb_from_credset_",
                           "PI",  # change when needed
                           ".txt",sep=""),  header=T)
  Wood=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                             "2017-06-06",
                             "_T2D_phys_loci_",d,"_kb_from_credset_",
                             "Wood",  # change when needed
                             ".txt",sep=""),  header=T)
  
  joint=rbind(BC,HG,PI,Wood)
  joint=unique(joint) # eliminate duplicate rows created by the way the Wood df was created
  
  colnames(joint)[1]=c("GeneID")
  joint$GeneID=as.character(joint$GeneID)
  
  df_list <- list()
  
  overlap_list <- list()  # to write down number of T2D genes in diff expr results
  
  prob_results=list()
  sums_probs=list()
  b_test <- list()
  overlap <- list()
  
  for(z in stage){
    
    # p-value for enrichment in T2D genes using my differentially expressed genes for each stage
    prob_results[[z]]=1-phyper(q=length(which(joint$GeneID %in% sig_stages[[z]]$GeneID))-1,
                               m=length(which(joint$GeneID %in% dge_cc$genes$ensembl_gene_id)),
                               n=length(dge_cc$genes$ensembl_gene_id) - length(which(joint$GeneID %in% dge_cc$genes$ensembl_gene_id)),
                               k=length(sig_stages[[z]]$GeneID))
    gene_list=joint[which(joint$GeneID %in% sig_stages[[z]]$GeneID),]
    
    overlap[["total_in_background"]] = length(which(joint$GeneID %in% dge_cc$genes$ensembl_gene_id))
    
    overlap[[z]] = length(which(joint$GeneID %in% sig_stages[[z]]$GeneID))
    
    # same, for random set of genes of same size. random gene set distribution of p-values
    random=list()
    random_overlaps=numeric()
    probab_random=numeric()
    
    for(i in 1:10000){
      random[[i]] <- as.data.frame(phyperRandom( joint$GeneID, sig_stages[[z]]$GeneID, dge_cc$genes$ensembl_gene_id))  # get probs from random gene sets, for each diff expr set tested
      probab_random[i] <- random[[i]]$prob
      random_overlaps[i] <- random[[i]]$n.overlap
      
    }
    
    
    hist(random_overlaps)
    
    # empirical p-value 
    
    sums_probs[[z]]$pval=((sum(probab_random<=prob_results[[z]])+1)/(10000+1)) # fraction of probabilities from random draws that are more extreme (smaller) than the one from my tested set
    # +1 in case sum(probab_random<=prob_results[[z]]) ==0
    # that would be the permuted? p-value
    
    # calculate the 95% confidence interval from the binomial distribution
    b_test[[z]]=binom.test(sum(probab_random<=prob_results[[z]]),10000)  # If those fractions are my "successes", then the binomial test can give me a confidence interval
    sums_probs[[z]]$conf.int.low=b_test[[z]]$conf.int[1]
    sums_probs[[z]]$conf.int.up=b_test[[z]]$conf.int[2]
    
    
    print(z)
    
  }
  overlap_list[[d]] <- do.call("rbind", overlap)
  
  df_list[[d]] <- do.call("rbind", sums_probs)
  
  # ######### save tables
  # ordered by p-value
  
  # 
  # write.table(df_list[[d]],quote=F,row.names=T,col.names = T,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
  #                                                                       "top",top,"_DEA_timecourse_pvals_hyperg_",
  #                                                                       currentDate,
  #                                                                       "_genes_enrichment_10k_permutations_",
  #                                                                       d,
  #                                                                       "_kb_from_credible_regions_T2D_",
  #                                                                       "BC_Wood_HG_PI",
  #                                                                       ".txt",sep=""),sep="\t")
  # write.table(overlap_list[[d]],quote=F,row.names=T,col.names = F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
  #                                                                            "top",top,"_DEA_timecourse_overlaps_hyperg_",
  #                                                                            currentDate,
  #                                                                            "_genes_enrichment_10k_permutations_",
  #                                                                            d,
  #                                                                            "_kb_from_credible_regions_T2D_",
  #                                                                            "BC_Wood_HG_PI",
  #                                                                            ".txt",sep=""),sep="\t")
  # 
  # ordered by logFc
  
  write.table(df_list[[d]],quote=F,row.names=T,col.names = T,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                                                                        "top",top,"_DEA_timecourse_pvals_hyperg_",
                                                                        currentDate,
                                                                        "_genes_enrichment_10k_permutations_",
                                                                        d,
                                                                        "_kb_from_credible_regions_T2D_",
                                                                        "BC_Wood_HG_PI",
                                                                        "orderedlogFC.txt",sep=""),sep="\t")
  write.table(overlap_list[[d]],quote=F,row.names=T,col.names = F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                                                                             "top",top,"_DEA_timecourse_overlaps_hyperg_",
                                                                             currentDate,
                                                                             "_genes_enrichment_10k_permutations_",
                                                                             d,
                                                                             "_kb_from_credible_regions_T2D_",
                                                                             "BC_Wood_HG_PI",
                                                                             "orderedlogFC.txt",sep=""),sep="\t")
  
}



##############plot #############################
# read in to save time if I just want to plot
top="623"
classifier="IR"

df_list <- list()
for (d in distance){
  df_list[[d]]=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                                     "top",top,"/",
                                     "each_classifier_separately/",
                          "top",top,"_DEA_timecourse_pvals_hyperg_",
                          "2017-06-06",
                          "_genes_enrichment_10k_permutations_",
                          d,
                          "_kb_from_credible_regions_T2D_",
                          classifier,
                          ".txt",sep=""),header=T)
  
  df_list[[d]]$stage=rownames(df_list[[d]])
  df_list[[d]]$distance=rep(d,8)
  
}

df= do.call("rbind", df_list)
colnames(df)=c("pval","lowCI","highCI","stage","distance")

df$stage=ordered(df$stage,levels = stage)
df$distance=ordered(df$distance,levels = distance)

# plot p-values and confidence intervals
# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right
diaPalette <- c("#CADAE8", "#7883BA", "#755A91", "#CC85B1", "#C15858", "#F4B8B0", "#96665A", "#6DA567")  # Diabetologia palette


# # x axis distance
# # -log10

df_log=df
df_log$pval=-log10(df$pval)
df_log$lowCI=-log10(df$lowCI)
df_log$highCI=-log10(df$highCI)

colnames(df_log)[4]="Stages"
levels(df_log$Stages) <- c(levels(df_log$Stages), c("EN","BLC","GT","PF"))  # change name of stages, first increase levels of factor
df_log[which(df_log$Stages=="EN6"),"Stages"] <- "EN"
df_log[which(df_log$Stages=="EN7"),"Stages"] <- "BLC"
df_log[which(df_log$Stages=="PGT"),"Stages"] <- "GT"
df_log[which(df_log$Stages=="PFG"),"Stages"] <- "PF"

df_log <- droplevels(df_log)  # drop unused levels (old names)

df_log$Stages <- ordered(df_log$Stages, levels = c("iPSC", "DE", "GT", "PF", "PE", "EP","EN", "BLC") ) # order levels, for colours


p=ggplot(df_log, aes(x=distance, y=pval, colour=Stages, group=Stages)) +
  geom_errorbar(aes(ymin=lowCI, ymax=highCI), colour="black", width=.5, size=.9, position=pd) +
  geom_line(position=pd,size=1.5) +
  scale_colour_manual(values=diaPalette) +
  geom_point(position=pd, size=4, shape=21, fill="white") +
  xlab("Distance (kb)") +
  ylab("Permuted p-value (-log10)") +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05),linetype="dashed",size=0.8,col="red") +
  theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           panel.border = element_rect(fill = NA, colour = "black"),
           legend.key = element_blank(), # legend.position = c(0.5,0.5),
           axis.title.y = element_text(face="bold", angle=90, size=16, vjust=0.2),
           axis.title.x = element_text(face="bold", size=16, vjust=0),
           axis.text.x = element_text(face="bold", colour = "black", angle=90, size=16, vjust=0.2, hjust =1 ),
           axis.text.y = element_text(face="bold", colour = "black",size=16),
           axis.ticks = element_line(colour = "black"),
           axis.line = element_line(colour = "black"),
           legend.text = element_text(face="bold", colour = "black",size=14),
           legend.title = element_text(face="bold", colour = "black",size=16))

ggsave(p,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Physiological_loci/",
                  "top",top,"/",
                  "each_classifier_separately/",
                  "top",top,"_DEA_timecourse_pvals_hyperg_",
                  "2017-06-06",
                  "_genes_enrichment_10k_permutations_",
                  "all_kb_from_credible_regions_T2D_",
                  classifier,
                  ".png",sep=""),
       width=10,height=8,units="in",dpi=300)

