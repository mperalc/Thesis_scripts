
# the files created by this script get larger with each iteration (doesn't overwrite files) delete them before re-running this script!!!!!!

# change for saving gene names or entrez ids

library(biomaRt)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))  # map ensembl id to entrez id



big_list <- read.table(file = "/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/Magenta/AllMAGENTAGeneAssocScores_Diamante_2018_T2D_thesis_110_40_Oct06_19.allgenescores", header = T)  # read in magenta
big_list = big_list[c(1,3,2)]
colnames(big_list) = c("gene_name", "p-value", "entrezgene")

big_list = big_list[order(big_list$`p-value`, decreasing = F), ]  # rank, lower values on top

# whatever genes have infinite value (read in as 0, check), assign to them the next lowest p-value
big_list[which(big_list$`p-value` == "0"), 2] = rep(unique(big_list$`p-value`)[2], times = length(big_list[which(big_list$`p-value` == "0"), 2]))


big_list <- big_list[which(!big_list$`p-value` == "NaN"), ]  # remove NaN p-vals
big_list <- big_list[which(big_list$`p-value` < 0.05), ]  # remove p-vals that are not significant

big_list$`p-value` = -log10(big_list$`p-value`)  # convert p-values to -log10

big_list = big_list[order(big_list$`p-value`, decreasing = T), ]  # rank, higher values on top

# genes <- as.character(input$external_gene_name)  # get my list of names
# entrez_ids <- getBM(filters = "external_gene_name", attributes = c("external_gene_name", "entrezgene"), values = genes, mart = mart)
# entrez_ids <- entrez_ids[order(match(entrez_ids$external_gene_name, input$external_gene_name)), ]  # order as input
# output = merge(input, entrez_ids, by = "external_gene_name")

GMT_T2D_ids = c("T2D", "NA", as.character(big_list$entrezgene))  # 366 genes p-val <0.01
# 1125 genes with p-val<0.05

# save with entrez ids:
write(unique(GMT_T2D_ids), file = "/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/Magenta/magenta_0.05_T2D_HRC.gmt", ncolumns = length(GMT_T2D_ids), sep = "\t", append = T)

## save with gene names
GMT_T2D_name = c("T2D", "NA", as.character(big_list$gene_name))  # 366 genes p-val <0.01
write(unique(GMT_T2D_name), file = "/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/Magenta/magenta_0.05_T2D_HRC_genename.gmt", ncolumns = length(GMT_T2D_ids), sep = "\t", append = T)







