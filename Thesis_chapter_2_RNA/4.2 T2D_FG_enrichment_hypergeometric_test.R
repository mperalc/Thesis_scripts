## Hypergeometric test for datasets: are they enriched in T2D/FG genes?

# for each gene Diff expr in each stage
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

currentDate <- Sys.Date() # to save date in name of output files

stage = c("iPSC", "DE", "PGT", "PFG", "PE", "EP", "EN6", "EN7")

trait = "T2D"    # T2D or fasting glucose (FG)?
type = "all"     #or top623
order = "p-value"   # or logFC
in_folder_DE = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/across-stages"  # folder that contains diff expression data 
# not removing monogenic diabetes genes

### files

#### background "genome" data
# just data used for the Differential expression analysis

 load("/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/session_objects/dge_cc.xz")  # 15221 genes and lincRNA


distance = c("0","50","100","200","500") # list of distances in kb

test <- list()


for (d in distance) {
  if (trait == "FG") {
    test[[d]] = read.table(
      file = paste(
        "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
        trait,
        "/2017-05-03FG_annotated_genes_in_credible_regions_plusminus_",
        d,
        "_kb.txt",
        sep = ""
      ),
      header = T
    )
  }
  else{
    test[[d]] = read.table(
      file = paste(
        "/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/Genes_in_credible_regions_T2D_HRC_2018",
        "/2019-10-06T2D_annotated_genes_in_credible_regions_plusminus_",
        d,
        "_kb.txt",
        sep = ""
      ),
      header = T
    )
    
  }
  colnames(test[[d]])[1] = c("GeneID")
  test[[d]]$GeneID = as.character(test[[d]]$GeneID)
}


# files with DE results (across-stages)

# # all genes, get list of ensembl gene ids

sig_stages=list()
for(s in stage){
 if(order=="p-value"){
   sig_stages[[s]] <- read.csv(file=paste(in_folder_DE,"/2017-07-03_sig_",s,"_diff_expression_maxvals_across-stages_results_logFC1.csv",sep=""),header = T) # read in DEA data for each stage
   colnames(sig_stages[[s]])[1]="GeneID"
   sig_stages[[s]]=sig_stages[[s]][ order(sig_stages[[s]]$adj.P.Val,decreasing=F),] # order by adjusted p.value 
 }
  else{
    sig_stages[[s]] <- read.csv(file=paste(in_folder_DE,"/2017-07-03_sig_",s,"_diff_expression_maxvals_across-stages_results_logFC1.csv",sep=""),header = T) # read in DEA data for each stage
    colnames(sig_stages[[s]])[1]="GeneID"
    sig_stages[[s]]=sig_stages[[s]][ order(sig_stages[[s]]$logFC,decreasing=T),] # order by logFC 
  }
  if(type=="top623"){
    
    # top 623 genes (min significant DE genes, in PGT stage)
    sig_stages[[s]]=sig_stages[[s]][c(1:623),]
  }

  }

##################
# hypergeometric test with permutations

# Total balls in the urn: background dataset (all RNA-seq genes)
# x --> number of white balls drawn without replacement from an urn which contains both black and white balls (T2D genes in my tested dataset)
# m --> White balls in the urn: characteristic I'm testing for enrichment (T2D genes)
# n --> Black balls in the urn: total-white balls (non-T2D genes)
# k --> Number of balls drawn from the urn: depends on the size of the dataset tested 

# I want the p-value of getting (number of T2D genes in my DE dataset) or more white balls in a sample 
# of size (number of genes in my DE dataset) from an urn with (number of T2D genes) white balls and 
# (total genes RNAseq-T2D genes in that set) black balls.

# 1-phyper(q, m, n, k)

# q in this case is x-1 (diff exp dataset -1), the probability of x or more. It is a quantile defined as the smallest value x such that F(x) ≥ p, 
# where F is the distribution function.
# 


# comparison of probabilities

phyperRandom <- function(myGeneList, myGeneSet, genome){
  myRandomGS <- sample( genome,size=length(myGeneSet) )  # take random sample of genes from whole genome, same size as set tested
  myX <- length(which(myGeneList %in% myRandomGS))  # number of T2D in random set
  myM <- length(which(myGeneList %in% genome)) #  T2D genes in my pool
  myN <- length(genome) - length(which(myGeneList %in% genome))  # genes in genome not T2D 
  myK <- length(myGeneSet) # size of set tested
  return(list(prob=1-phyper(q=myX-1,
                            m=myM, n=myN, k=myK), n.overlap=myX))
}

df_list <- list()
distance = c("0","50","100","200","500") # list of distances in kb
overlap_list <- list()  # to write down number of T2D genes in diff expr results

for(d in distance){
  
  prob_results=list()
  sums_probs=list()
  b_test <- list()
  overlap <- list()
  
  for(z in stage){
    
    # p-value for enrichment in T2D genes using my differentially expressed genes for each stage
    prob_results[[z]]=1-phyper(q=length(which(test[[d]]$GeneID %in% sig_stages[[z]]$GeneID))-1,
                               m=length(which(test[[d]]$GeneID %in% dge_cc$genes$ensembl_gene_id)),
                               n=length(dge_cc$genes$ensembl_gene_id) - length(which(test[[d]]$GeneID %in% dge_cc$genes$ensembl_gene_id)),
                               k=length(sig_stages[[z]]$GeneID))
    gene_list=test[[d]][which(test[[d]]$GeneID %in% sig_stages[[z]]$GeneID),]
   
    overlap[["total_in_background"]] = length(which(test[[d]]$GeneID %in% dge_cc$genes$ensembl_gene_id))
    
    overlap[[z]] = length(which(test[[d]]$GeneID %in% sig_stages[[z]]$GeneID))
    
    # same, for random set of genes of same size. random gene set distribution of p-values
    random=list()
    random_overlaps=numeric()
    probab_random=numeric()
    
    for(i in 1:10000){
      random[[i]] <- as.data.frame(phyperRandom( test[[d]]$GeneID, sig_stages[[z]]$GeneID, dge_cc$genes$ensembl_gene_id))  # get probs from random gene sets, for each diff expr set tested
      probab_random[i] <- random[[i]]$prob
      random_overlaps[i] <- random[[i]]$n.overlap
      
    }
    
    
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
  
  # if all genes
  if (type=="all"){
    write.table(df_list[[d]],quote=F,row.names=T,col.names = T,file=paste("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/Genes_in_credible_regions_T2D_HRC_2018/",
                                                                          currentDate,
                                                                          "_all_DEA_timecourse_pvals_hyperg_",
                                                                          trait,
                                                                          "_genes_enrichment_10k_permutations_",
                                                                          d,"_",
                                                                          "with_monogenic",
                                                                          "_kb_from_credible_regions.txt",sep=""),sep="\t")
    write.table(overlap_list[[d]],quote=F,row.names=T,col.names = F,file=paste("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/Genes_in_credible_regions_T2D_HRC_2018/",
                                                                          currentDate,
                                                                          "_all_DEA_timecourse_overlaps_hyperg_",
                                                                          trait,
                                                                          "_genes_enrichment_",
                                                                          d,"_",
                                                                          "with_monogenic",
                                                                          "_kb_from_credible_regions.txt",sep=""),sep="\t")
  }else{
    # if top 623 genes
    #ordered by logFC
    if(order=="logFC"){
      write.table(df_list[[d]],quote=F,row.names=T,col.names = T,file=paste("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/Genes_in_credible_regions_T2D_HRC_2018/",
                                                                            currentDate,
                                                                            "_top623_orderedByLogFC_DEA_timecourse_pvals_hyperg_",
                                                                            trait,
                                                                            "_genes_enrichment_10k_permutations_",
                                                                            d,"_",
                                                                            "with_monogenic",
                                                                            "_kb_from_credible_regions.txt",sep=""),sep="\t")
      write.table(overlap_list[[d]],quote=F,row.names=T,col.names = F,file=paste("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/Genes_in_credible_regions_T2D_HRC_2018/",
                                                                                 currentDate,
                                                                                 "_top623_orderedByLogFC_DEA_timecourse_overlaps_hyperg_",
                                                                                 trait,
                                                                                 "_genes_enrichment_",
                                                                                 d,"_",
                                                                                 "with_monogenic",
                                                                                 "_kb_from_credible_regions.txt",sep=""),sep="\t")
    }else{
      # ordered by p-value
      write.table(df_list[[d]],quote=F,row.names=T,col.names = T,file=paste("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/Genes_in_credible_regions_T2D_HRC_2018/",
                                                                            currentDate,
                                                                            "_top623_DEA_timecourse_pvals_hyperg_",
                                                                            trait,
                                                                            "_genes_enrichment_10k_permutations_",
                                                                            d,"_",
                                                                            "with_monogenic",
                                                                            "_kb_from_credible_regions.txt",sep=""),sep="\t")
      write.table(overlap_list[[d]],quote=F,row.names=T,col.names = F,file=paste("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/Genes_in_credible_regions_T2D_HRC_2018/",
                                                                                 currentDate,
                                                                                 "_top623_DEA_timecourse_overlaps_hyperg_",
                                                                                 trait,
                                                                                 "_genes_enrichment_",
                                                                                 d,"_",
                                                                                 "with_monogenic",
                                                                                 "_kb_from_credible_regions.txt",sep=""),sep="\t")
    }
    
    
    
    
  }
}
