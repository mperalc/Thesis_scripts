# transform list of T2D genes from credible sets to gmt

# the files created by this script get larger with each iteration (doesn't overwrite files) delete them before re-running this script!!!!!!

# change for saving gene names or entrez ids

trait = "T2D"
library(biomaRt)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))  # map ensembl id to entrez id


distance = c("0", "50", "100", "200", "500")  # list of distances in kb
for (d in distance) {
    if (trait == "T2D") {
        input = read.table(file = paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_", trait, "/2017-02-22", trait, "annotated_genes_in_credible_regions_plusminus_", 
            d, "_kb.txt", sep = ""), header = T)
    } else {
        input = read.table(file = paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_", trait, "/2017-05-03", trait, "_annotated_genes_in_credible_regions_plusminus_", 
            d, "_kb.txt", sep = ""), header = T)
    }
    
    input$ensembl_gene_id = as.character(input$ensembl_gene_id)
    input$external_gene_name = as.character(input$external_gene_name)
    
    
    genes <- as.character(input$external_gene_name)  # get my list of names
    entrez_ids <- getBM(filters = "external_gene_name", attributes = c("external_gene_name", "entrezgene"), values = genes, mart = mart)
    entrez_ids <- entrez_ids[order(match(entrez_ids$external_gene_name, input$external_gene_name)), ]  # order as input
    output = merge(input, entrez_ids, by = "external_gene_name")
    
    output = output[which(!is.na(output$entrezgene)), ]  # remove NAs
    
    GMT_T2D_ids = c(d, "NA", as.character(output$entrezgene))
    
    # save with entrez ids:
    write(unique(GMT_T2D_ids), file = paste("/Users/Marta/Documents/WTCHG/DPhil/Data/input_for_GSEA/", trait, "_credsets.gmt", sep = ""), ncolumns = length(GMT_T2D_ids), sep = "\t", 
        append = T)
}




