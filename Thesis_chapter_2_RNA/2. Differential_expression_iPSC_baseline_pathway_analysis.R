### pathway analysis differentially expressed genes

setwd("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis")

library(GenomicRanges)
library(ChIPseeker)
library(UpSetR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)


## differential expression results
stage = c( "iPSC"  ,   "DE"   ,    "PGT"  ,    "PFG" ,     "PE" ,"EP",      "EN6", "EN7")
diff_expr = list()
for(s in stage){
  diff_expr[[s]] = read.csv(paste0("differential_expression_results/",
                                   "2018-05-09_sig_",s,"_diff_expression_maxvals_logFC1.csv"))
  diff_expr[[s]] = diff_expr[[s]]$ensembl_gene_id
}

load(file="dge_counts.xz")

ensembl38 = useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  host = "aug2017.archive.ensembl.org",
  path = "/biomart/martservice" ,
  dataset = "hsapiens_gene_ensembl"
)
all_ensembl_to_entrez = getBM(
  attributes = c('ensembl_gene_id', 'entrezgene'),
  filters = 'ensembl_gene_id',
  values = dge_counts$genes$ensembl_gene_id,
  mart = ensembl38
)  # all genes from dge object
all_ensembl_to_entrez = na.omit(all_ensembl_to_entrez$entrezgene)

diff_expr_entrez = list()
for(s in stage) {
  diff_expr_entrez[[s]] = getBM(
    attributes = c('ensembl_gene_id', 'entrezgene'),
    filters = 'ensembl_gene_id',
    values = diff_expr[[s]],
    mart = ensembl38
  )
  diff_expr_entrez[[s]] = as.character(na.omit(diff_expr_entrez[[s]]$entrezgene))
}

gene_universe_annotated = enrichPathway(all_ensembl_to_entrez, readable = T)

names(diff_expr_entrez) = c("iPSC", "DE" ,  "GT",  "PF" , "PE",   "EP" ,  "EN" , "BLC" )
compEnrich_pathway <- compareCluster(geneCluster   = diff_expr_entrez,
                                     fun           = "enrichPathway",
                                     pvalueCutoff  = 0.05,
                                     pAdjustMethod = "BH",
                                     universe = as.character(all_ensembl_to_entrez),
                                     readable = T)


compEnrich_pathway@compareClusterResult$geneID = gsub("/",",",compEnrich_pathway@compareClusterResult$geneID)


write.table(compEnrich_pathway@compareClusterResult,
            file = "pathway_analysis_diff_expr.txt",
            col.names = T, row.names = F, quote = F,sep = "\t")

png(
   'pathway_analysis_diff_expr_top10.png',units = "in",res = 400,type = "cairo",
  width = 12,
  height = 12
)

dotplot(compEnrich_pathway, showCategory = 10, title = "Reactome Pathway Enrichment Analysis")

dev.off()

