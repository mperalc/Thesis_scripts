# 
setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC")

library(readr)
library(summarytools)
# Differentiation stages

stage = c("iPSC", "DE", "GT", "PF", "PE", "EP", "EN")
# Now loading gene stats
gene_stats = list()
for (s in stage) {
  gene_stats[[s]] = read_delim(file = paste0("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/",s,"/GenePredictionStats.txt"), delim = "\t")
  gene_stats[[s]]$CellType = rep(s,nrow(gene_stats[[s]]))
}

gene_stats = do.call("rbind",gene_stats)

# Number of distal enhancers per gene 
descr(gene_stats$nDistalEnhancersPredicted, style = "rmarkdown")

# Only for expressed genes
descr(gene_stats[gene_stats$gene_is_expressed_proxy == TRUE,"nDistalEnhancersPredicted"], style = "rmarkdown")

# # General summary
# dfSummary(subset_for_factor_analysis, plain.ascii = FALSE, style = "grid", 
#           graph.magnif = 0.75, valid.col = FALSE, tmp.img.dir = ".")


# Distal enhancers per gene and stage

descr(gene_stats[gene_stats$CellType== "iPSC" & gene_stats$gene_is_expressed_proxy == TRUE,"nDistalEnhancersPredicted"], style = "rmarkdown")
descr(gene_stats[gene_stats$CellType== "DE" & gene_stats$gene_is_expressed_proxy == TRUE,"nDistalEnhancersPredicted"], style = "rmarkdown")
descr(gene_stats[gene_stats$CellType== "GT" & gene_stats$gene_is_expressed_proxy == TRUE,"nDistalEnhancersPredicted"], style = "rmarkdown")
descr(gene_stats[gene_stats$CellType== "PF" & gene_stats$gene_is_expressed_proxy == TRUE,"nDistalEnhancersPredicted"], style = "rmarkdown")
descr(gene_stats[gene_stats$CellType== "PE" & gene_stats$gene_is_expressed_proxy == TRUE,"nDistalEnhancersPredicted"], style = "rmarkdown")
descr(gene_stats[gene_stats$CellType== "EP" & gene_stats$gene_is_expressed_proxy == TRUE,"nDistalEnhancersPredicted"], style = "rmarkdown")
descr(gene_stats[gene_stats$CellType== "EN" & gene_stats$gene_is_expressed_proxy == TRUE,"nDistalEnhancersPredicted"], style = "rmarkdown")
descr(gene_stats[gene_stats$CellType== "BLC" & gene_stats$gene_is_expressed_proxy == TRUE,"nDistalEnhancersPredicted"], style = "rmarkdown")
