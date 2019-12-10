# Rank pvals output file (from DEA) and create .txt file for GSEA

### variable

out_folder = "/Users/Marta/Documents/WTCHG/DPhil/Data/input_for_GSEA/"
select = "logFC"  # pval or logFC
############ 

library(biomaRt)

stages = c("iPSC", "DE", "PFG", "PGT", "PE", "EP", "EN6", "EN7")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))  # map ensembl id to entrez id


for (s in stages) {
    input <- read.csv(file = paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/across-stages/2017-07-03_sig_", s, "_diff_expression_maxvals_across-stages_results_logFC1.csv", 
        sep = ""), header = T)
    input$adj.P.Val = -log10(input$adj.P.Val)
    
    if (select == "pval") {
        input = input[order(input$adj.P.Val, decreasing = T), ]  # higher p-values on top
        
        write.table(input[, c("external_gene_name", "adj.P.Val")], file = paste(out_folder, "gene_names_", s, "_DEA_for_GSEA.rnk", sep = ""), quote = F, col.names = F, row.names = F, 
            sep = "\t")
    } else {
        input = input[order(input$logFC, decreasing = T), ]  # higher p-values on top
        
        write.table(input[, c("external_gene_name", "logFC")], file = paste(out_folder, "gene_names_", s, "_DEA_for_GSEA_ordered_by_logFC.rnk", sep = ""), quote = F, col.names = F, 
            row.names = F, sep = "\t")
    }
    genes <- as.character(input$external_gene_name)  # get my list of names
    entrez_ids <- getBM(filters = "external_gene_name", attributes = c("external_gene_name", "entrezgene"), values = genes, mart = mart)
    entrez_ids <- entrez_ids[order(match(entrez_ids$external_gene_name, input$external_gene_name)), ]  # order as input
    output = merge(input, entrez_ids, by = "external_gene_name")
    
    output <- output[which(!is.na(output$entrezgene)), ]  # remove Na 
    
    if (select == "pval") {
        output = output[, c("entrezgene", "adj.P.Val")]
        output = output[order(output$adj.P.Val, decreasing = T), ]  # higher p-values on top
        
        write.table(output, file = paste(out_folder, "gene_ids_", s, "_DEA_for_GSEA.rnk", sep = ""), quote = F, col.names = F, row.names = F, sep = "\t")
        
    } else {
        
        output = output[, c("entrezgene", "logFC")]
        output = output[order(output$logFC, decreasing = T), ]  # higher p-values on top
        
        write.table(output, file = paste(out_folder, "gene_ids_", s, "_DEA_for_GSEA_ordered_by_logFC.rnk", sep = ""), quote = F, col.names = F, row.names = F, sep = "\t")
        
        
    }
}
