# modify output of MAGENTA to make suitable GSEA input always remember to delete the .gmt files after each run!  variable

setwd("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/Magenta")
MAGENTA = read.table("AllMAGENTAGeneAssocScores_Diamante_2018_T2D_thesis_110_40_Oct06_19.allgenescores", header = T)
MAGENTA = MAGENTA[c("Gene_Symbol", "Entrez_ID", "Gene_p.value")]

output_filename = "AllMAGENTAGeneAssocScores_Diamante_2018_T2D_thesis_110_40_Oct06_19"

keep_nosig = F  # keep pvalues that are not significant?
pval = 0.01  # p-value threshold
######### 

# Give p-values of 0 the next lowest p-value
MAGENTA[which(MAGENTA$Gene_p.value == "0"), "Gene_p.value"] = rep(unique(MAGENTA$Gene_p.value)[2], times = length(MAGENTA[which(MAGENTA$Gene_p.value == "0"), "Gene_p.value"]))

MAGENTA <- MAGENTA[which(!MAGENTA$Gene_p.value == "NaN"), ]  # remove NaN p-vals


if (keep_nosig == F) {
    MAGENTA <- MAGENTA[which(MAGENTA$Gene_p.value < pval), ]  # remove p-vals that are not significant
}


MAGENTA$Gene_p.value = -log10(MAGENTA$Gene_p.value)  # convert p-values to -log10

MAGENTA = MAGENTA[order(MAGENTA$Gene_p.value, decreasing = T), ]  # rank, higher values on top
MAGENTA$Gene_Symbol = as.factor(MAGENTA$Gene_Symbol)
######### save rnk file save with entrez ids:

if (keep_nosig == T) {
    write.table(MAGENTA[c("Entrez_ID", "Gene_p.value")], file = paste("entrez_ids_", output_filename, ".rnk", sep = ""), sep = "\t", row.names = F, col.names = F)
} else {
    write.table(MAGENTA[c("Entrez_ID", "Gene_p.value")], file = paste("entrez_ids_", pval, "_", output_filename, ".rnk", sep = ""), sep = "\t", row.names = F, col.names = F)
}

# save with gene names

if (keep_nosig == T) {
    write.table(MAGENTA[c("Gene_Symbol", "Gene_p.value")], file = paste("gene_name_", output_filename, ".rnk", sep = ""), sep = "\t", row.names = F, col.names = F, quote = F)
} else {
    write.table(MAGENTA[c("Gene_Symbol", "Gene_p.value")], file = paste("gene_name_", pval, "_", output_filename, ".rnk", sep = ""), sep = "\t", row.names = F, col.names = F, quote = F)
}
# save with entrez ids:

######### save gmt file save with entrez ids:
GMT_T2D_ids = c("T2D", "NA", as.character(MAGENTA$Entrez_ID))

if (keep_nosig == T) {
    write(unique(GMT_T2D_ids), file = paste("entrez_ids_", output_filename, ".gmt", sep = ""), ncolumns = length(GMT_T2D_ids), sep = "\t", append = T)
    
} else {
    write(unique(GMT_T2D_ids), file = paste("entrez_ids_", pval, "_", output_filename, ".gmt", sep = ""), ncolumns = length(GMT_T2D_ids), sep = "\t", append = T)
    
}

# save with gene names

GMT_T2D_name = c("T2D", "NA", as.character(MAGENTA$Gene_Symbol))


if (keep_nosig == T) {
    write(unique(GMT_T2D_name), file = paste("gene_name_", output_filename, ".gmt", sep = ""), ncolumns = length(GMT_T2D_name), sep = "\t", append = T)
    
} else {
    write(unique(GMT_T2D_name), file = paste("gene_name_", pval, "_", output_filename, ".gmt", sep = ""), ncolumns = length(GMT_T2D_name), sep = "\t", append = T)
}
