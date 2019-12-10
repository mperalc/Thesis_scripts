
# the files created by this script get larger with each iteration (doesn't overwrite files) delete them before re-running this script!!!!!!

# change for saving gene names or entrez ids

library(biomaRt)

##### variable part

setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/across-stages")


stages = c("iPSC", "DE", "PFG", "PGT", "PE", "EP", "EN6", "EN7")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))  # map ensembl id to entrez id
all = T  # use all genes or only top 500?

for (s in stages) {
    
    big_list <- read.csv(file = paste("2017-07-03_sig_", s, "_diff_expression_maxvals_across-stages_results_logFC1.csv", sep = ""), header = T)  # read in DEA data for each stage
    
    # GMT = c(s,'NA',as.character(big_list[[s]]$gene_name)) # select ordering by peak and timecourse p-value?
    # write(GMT,file='/Users/Marta/Documents/WTCHG/DPhil/Data/input_for_GSEA/DEA_peak_timecourse_for_GSEA.gmt',ncolumns =length(GMT),sep='\t',append=T) #writes vector of gene list
    # to this file line by line select = big_list[which(big_list$sig_peak.timecourse=='both' | big_list$sig_peak.timecourse=='peak' ),] select=select[
    # order(select$adjust.pval_peak,decreasing=F),] # order by adjusted p.value of 'peak' (stage-specific) method # ensebl id to entrez id
    # ensembl_genes<-as.character(select$ensembl_gene_id) # get my list of ensembl gene ids entrez_ids <- getBM( filters= 'ensembl_gene_id', attributes=
    # c('ensembl_gene_id','entrezgene'), values= ensembl_genes, mart= mart) entrez_ids <- entrez_ids[order(match(entrez_ids$ensembl_gene_id,select$ensembl_gene_id)),] # order as
    # 'select' list select$entrez_id <- entrez_ids[unique(match(entrez_ids$ensembl_gene_id,select$ensembl_gene_id)),2] # select unique matching rows
    # select=select[c(1,ncol(select),2:(ncol(select)-1))] # reorder columns select <- select[which(!select$entrez_id=='NA'),] # remove Nas GMT_peak_ids =
    # c(s,'NA',as.character(select$entrez_id)) GMT_peak = c(s,'NA',as.character(select$gene_name)) # will have more items than the id one, due to lost mapping of ensembl-entrez ids
    # GMT_peak = GMT_peak[c(1:502)] # select 500 genes + stage identifier + description ('NA') GMT_peak = GMT_peak[which(!is.na(GMT_peak))] # remove NAs GMT_peak_ids =
    # GMT_peak_ids[c(1:502)] # same for entrez ids GMT_peak_ids = GMT_peak_ids[which(!is.na(GMT_peak_ids))] # remove NAs # save with gene names: #
    # write(GMT_peak,file='/Users/Marta/Documents/WTCHG/DPhil/Data/input_for_GSEA/DEA_peak_for_GSEA.gmt',ncolumns =length(GMT_peak),sep='\t',append=T) # save with entrez ids:
    # write(GMT_peak_ids,file='/Users/Marta/Documents/WTCHG/DPhil/Data/input_for_GSEA/DEA_peak_for_GSEA_entrez_ids.gmt',ncolumns =length(GMT_peak),sep='\t',append=T)
    
    big_list = big_list[order(big_list$adj.P.Val, decreasing = F), ]  # order by adjusted p.value of 'timecourse' (across-stages) method
    
    ############### The following needs to be corrected (rank_pvals_for_GSEA_DEA works well in this part) ensebl id to entrez id ensembl_genes<-as.character(big_list$ensembl_gene_id) # get my
    ############### list of ensembl gene ids entrez_ids <- getBM( filters= 'ensembl_gene_id', attributes= c('ensembl_gene_id','entrezgene'), values= ensembl_genes, mart= mart) entrez_ids <-
    ############### entrez_ids[order(match(entrez_ids$ensembl_gene_id,big_list$ensembl_gene_id)),] # order as big list big_list$entrez_id <-
    ############### entrez_ids[unique(match(entrez_ids$ensembl_gene_id,big_list$ensembl_gene_id)),2] # select unique matching rows big_list=big_list[c(1,ncol(big_list),2:(ncol(big_list)-1))] #
    ############### reorder columns big_list <- big_list[which(!big_list$entrez_id=='NA'),] # remove Nas
    
    
    
    
    
    # GMT_timecourse_ids = c(s,'NA',as.character(big_list$entrez_id))
    GMT_timecourse = c(s, "NA", as.character(big_list$external_gene_name))
    # GMT_timecourse_ensembl=c(s,'NA',as.character(ensembl_genes))
    
    # GMT_timecourse_ids = GMT_timecourse_ids[which(!is.na(GMT_timecourse_ids))] # remove NAs
    GMT_timecourse = GMT_timecourse[which(!is.na(GMT_timecourse))]  # remove NAs
    
    # 
    if (all == F) {
        GMT_timecourse = GMT_timecourse[c(1:502)]  # select 500 genes + stage identifier + description ('NA') 
        # GMT_timecourse_ids = GMT_timecourse_ids[c(1:502)] # same for entrez ids
    }
    
    # GMT_timecourse_ensembl = GMT_timecourse_ensembl[c(1:502)] # same for ensembl ids GMT_timecourse_ensembl = GMT_timecourse_ensembl[which(!is.na(GMT_timecourse_ensembl))] #
    # remove NAs
    
    
    if (all == F) {
        # save with gene names:
        write(GMT_timecourse, file = "/Users/Marta/Documents/WTCHG/DPhil/Data/input_for_GSEA/DEA_timecourse_for_GSEA_top500.gmt", ncolumns = length(GMT_timecourse), sep = "\t", 
            append = T)
        # save with entrez ids: write(GMT_timecourse_ids,file='/Users/Marta/Documents/WTCHG/DPhil/Data/input_for_GSEA/DEA_timecourse_for_GSEA_entrez_ids_top500.gmt',ncolumns
        # =length(GMT_timecourse_ids),sep='\t',append=T)
    } else {
        # save with gene names:
        write(GMT_timecourse, file = "/Users/Marta/Documents/WTCHG/DPhil/Data/input_for_GSEA/DEA_timecourse_for_GSEA.gmt", ncolumns = length(GMT_timecourse), sep = "\t", append = T)
        # save with entrez ids: write(GMT_timecourse_ids,file='/Users/Marta/Documents/WTCHG/DPhil/Data/input_for_GSEA/DEA_timecourse_for_GSEA_entrez_ids.gmt',ncolumns
        # =length(GMT_timecourse_ids),sep='\t',append=T)
    }
    # save with ensembl ids: write(GMT_timecourse_ensembl,file='/Users/Marta/Documents/WTCHG/DPhil/Data/input_for_GSEA/DEA_timecourse_for_GSEA_ensembl_ids.gmt',ncolumns
    # =length(GMT_timecourse_ensembl),sep='\t',append=T)
    
}




