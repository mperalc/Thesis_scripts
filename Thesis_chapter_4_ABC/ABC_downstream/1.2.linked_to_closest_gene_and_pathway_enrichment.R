# Comparing linked genes to closest gene

library(GenomicRanges)
library(annotatr)
library(annotate)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(ReactomePA)
library(clusterProfiler)
library(DOSE)
library('org.Hs.eg.db')

stages = c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

## ABC output directly from run
subset_for_factor_analysis = list()
subset_for_factor_analysis.gr = list()
pred_per_stage = list()
for (s in stages) {
  pred_per_stage[[s]] = read_delim(file = paste0("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/",s,"/EnhancerPredictions.txt"), delim = "\t")
  subset_for_factor_analysis[[s]] = pred_per_stage[[s]][,c("chr", 
                                                           "start", 
                                                           "end", 
                                                           "cellType", 
                                                           "class",
                                                           "isPromoterElement",
                                                           "enhancerSymbol",
                                                           "name",
                                                           "normalized_h3K27ac",
                                                           "normalized_atac",
                                                           "distance",
                                                           "isSelfPromoter",
                                                           "TargetGene",
                                                           "TargetGeneExpression",
                                                           "TargetGeneTSS",
                                                           "TargetGeneExpression",
                                                           "TargetGenePromoterActivityQuantile",
                                                           "hic.distance.adj",
                                                           "ABC.Score")]
  subset_for_factor_analysis[[s]] = subset_for_factor_analysis[[s]][!is.na(subset_for_factor_analysis[[s]]$TargetGeneExpression),]
  subset_for_factor_analysis[[s]] = subset_for_factor_analysis[[s]][subset_for_factor_analysis[[s]]$chr %in% paste0("chr",1:22),]
  subset_for_factor_analysis.gr[[s]] = GRanges(subset_for_factor_analysis[[s]])
  
}


subset_for_factor_analysis_all = do.call("rbind",subset_for_factor_analysis )
nrow(subset_for_factor_analysis_all) # includes enhancers that map to multiple genes

## Make list to genomic ranges
subset_for_factor_analysis.gr = GRangesList(subset_for_factor_analysis.gr)
subset_for_factor_analysis_all.gr = GRanges(subset_for_factor_analysis_all)

# Those that are acting on closest gene


# I have checked and the UCSC and knowgene give slightly different annotations for some genes, mapping to promoter when they arent
refseq <- makeTxDbFromUCSC(genome="hg19", tablename="refGene")
subset_for_factor_analysis_all_closest_gene = annotatePeak(subset_for_factor_analysis_all.gr,
                                                           tssRegion = c(-2000,500),TxDb = refseq, verbose = T,
                                                           overlap = "TSS",  # VERY IMPORTANT
                                                           annoDb = "org.Hs.eg.db")
subset_for_factor_analysis_all_closest_gene = as.data.frame(subset_for_factor_analysis_all_closest_gene)

subset_for_factor_analysis_all_closest_gene = subset_for_factor_analysis_all_closest_gene[c(1:21,28,30,32)] # subset columns of interest
colnames(subset_for_factor_analysis_all_closest_gene)[c(22:24)] = c("closest_geneId","closest_distanceToTSS","closest_SYMBOL") # rename

subset_for_factor_analysis_all_closest_gene_match = subset_for_factor_analysis_all_closest_gene[subset_for_factor_analysis_all_closest_gene$TargetGene == subset_for_factor_analysis_all_closest_gene$closest_SYMBOL,]

not_closest_gene_match = subset_for_factor_analysis_all_closest_gene[subset_for_factor_analysis_all_closest_gene$TargetGene != subset_for_factor_analysis_all_closest_gene$closest_SYMBOL,]
nrow(not_closest_gene_match) / nrow(subset_for_factor_analysis_all_closest_gene)## percentage of total enhancers not acting on closest gene

nrow(unique(not_closest_gene_match[c("seqnames","start","end","cellType")])) # when taking into account only unique enhancers
# there are actually roughly the same amount of enhancers acting on the closest gene
nrow(unique(subset_for_factor_analysis_all_closest_gene_match[c("seqnames","start","end","cellType")]))
# but enhancers usually act on multiple genes


nrow(unique(subset_for_factor_analysis_all_closest_gene_match[c("seqnames","start","end","cellType")])) / nrow(unique(subset_for_factor_analysis_all_closest_gene[c("seqnames","start","end","cellType")]))## percentage of total enhancers acting on closest gene
# 
match = unique(subset_for_factor_analysis_all_closest_gene_match[c("seqnames","start","end","cellType")])
all = unique(subset_for_factor_analysis_all_closest_gene[c("seqnames","start","end","cellType")])
# # Per stage
nrow(match[match$cellType =="iPSC",])/
  nrow(all[all$cellType == "iPSC",])
nrow(match[match$cellType =="DE",])/
  nrow(all[all$cellType == "DE",])
nrow(match[match$cellType =="GT",])/
  nrow(all[all$cellType == "GT",])
nrow(match[match$cellType =="PF",])/
  nrow(all[all$cellType == "PF",])
nrow(match[match$cellType =="PE",])/
  nrow(all[all$cellType == "PE",])
nrow(match[match$cellType =="EP",])/
  nrow(all[all$cellType == "EP",])
nrow(match[match$cellType =="EN",])/
  nrow(all[all$cellType == "EN",])
nrow(match[match$cellType =="BLC",])/
  nrow(all[all$cellType == "BLC",])
mean(c( 0.6888156,0.6772246,0.6900119,0.7219877, 0.725073,0.6894132,  0.72088,0.6829296))

# Pathway analysis #   ########################



genes_target_entrez  = list()
genes_closest_entrez  = list()
for(s in stages){
genes_target = unique(subset_for_factor_analysis_all_closest_gene[subset_for_factor_analysis_all_closest_gene$cellType==s,"TargetGene"])
genes_closest = unique(subset_for_factor_analysis_all_closest_gene[subset_for_factor_analysis_all_closest_gene$cellType==s,"closest_SYMBOL"])

genes_target_entrez[[s]] = mapIds(org.Hs.eg.db, genes_target, 'ENTREZID', 'SYMBOL')
genes_closest_entrez[[s]] = mapIds(org.Hs.eg.db, genes_closest, 'ENTREZID', 'SYMBOL')
}

compEnrich_target <- compareCluster(geneCluster   = genes_target_entrez,
                             fun           = "enrichPathway",
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH",
                             readable = T)
compEnrich_closest <- compareCluster(geneCluster   = genes_closest_entrez,
                                    fun           = "enrichPathway",
                                    pvalueCutoff  = 0.05,
                                    pAdjustMethod = "BH",
                                    readable = T)

png(
  paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/", '/enrichment_reactome_top10_ABC_closest.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 12,
  height = 12
)

dotplot(compEnrich_closest, showCategory = 10, title = "Reactome Pathway Enrichment Analysis: closest gene")


dev.off()



png(
  paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/", '/enrichment_reactome_top10_ABC_target.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 12,
  height = 12
)

dotplot(compEnrich_target, showCategory = 10, title = "Reactome Pathway Enrichment Analysis: predicted gene")

dev.off()

write.csv(compEnrich_target,"/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/ABC_gene_target_pathway.csv",
            row.names = F,quote = F)
write.csv(compEnrich_closest,"/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/ABC_gene_closest_pathway.csv",
            row.names = F,quote = F)
# closest and target
write.csv(subset_for_factor_analysis_all_closest_gene,"/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/ABC_gene_closest_and_target.csv",
          row.names = F,quote = F)

# pathway analysis of those enhancers not shared across all stages (more reflective of biology of each stage?)

not_shared = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/for_fGWAS/fGWAS_ABC_output_enhancers_not_in_all_stages.bed")
colnames(not_shared) = c("chr","start","end","stage")

subset = list()
for(s in stages){
  not_shared.gr = GRanges(not_shared[not_shared$stage==s,])
  subset_for_factor_analysis_all_closest_gene.gr = GRanges(subset_for_factor_analysis_all_closest_gene[subset_for_factor_analysis_all_closest_gene$cellType==s,])
  subset[[s]] = subsetByOverlaps(subset_for_factor_analysis_all_closest_gene.gr, not_shared.gr)
  subset[[s]] =as.data.frame(subset[[s]])
}
subset= do.call("rbind",subset)
genes_target_entrez  = list()
genes_closest_entrez  = list()
for(s in stages){
  genes_target = unique(subset[subset$cellType==s,"TargetGene"])
  genes_closest = unique(subset[subset$cellType==s,"closest_SYMBOL"])
  
  genes_target_entrez[[s]] = mapIds(org.Hs.eg.db, genes_target, 'ENTREZID', 'SYMBOL')
  genes_closest_entrez[[s]] = mapIds(org.Hs.eg.db, genes_closest, 'ENTREZID', 'SYMBOL')
}

compEnrich_target <- compareCluster(geneCluster   = genes_target_entrez,
                                    fun           = "enrichPathway",
                                    pvalueCutoff  = 0.05,
                                    pAdjustMethod = "BH",
                                    readable = T)
compEnrich_closest <- compareCluster(geneCluster   = genes_closest_entrez,
                                     fun           = "enrichPathway",
                                     pvalueCutoff  = 0.05,
                                     pAdjustMethod = "BH",
                                     readable = T)

png(
  paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/", '/enrichment_reactome_top10_ABC_closest_not_shared_all_stages.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 12,
  height = 12
)

dotplot(compEnrich_closest, showCategory = 10, title = "Pathway Enrichment: ABC active enhancer mapped to closest gene")


dev.off()



png(
  paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/", '/enrichment_reactome_top10_ABC_target_not_shared_all_stages.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 12,
  height = 12
)

dotplot(compEnrich_target, showCategory = 10, title = "Pathway Enrichment: ABC active enhancer mapped to predicted gene")

dev.off()

write.csv(compEnrich_target,"/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/ABC_gene_target_pathway_not_shared_all_stages.csv",
          row.names = F,quote = F)
write.csv(compEnrich_closest,"/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/ABC_gene_closest_pathway_not_shared_all_stages.csv",
          row.names = F,quote = F)
