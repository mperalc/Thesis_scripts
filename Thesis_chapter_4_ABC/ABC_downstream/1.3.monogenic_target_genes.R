# How many monogenic diabetes genes in the output of ABC

stages = c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

## ABC output directly from run
subset_for_factor_analysis = list()
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

}


subset_for_factor_analysis_all = do.call("rbind",subset_for_factor_analysis )
nrow(subset_for_factor_analysis_all) # includes enhancers that map to multiple genes

monogenic = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Monogenic_diabetes/Monogenic_updated_2019.txt")

monogenic$V1 %in% subset_for_factor_analysis_all$TargetGene
sum(monogenic$V1 %in% subset_for_factor_analysis_all$TargetGene) / length(monogenic$V1)
#29/34
monogenic[monogenic$V1 %in% subset_for_factor_analysis_all$TargetGene,]

# analysing closest target

closest_and_target = read.csv("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/ABC_gene_closest_and_target.csv")

sum(monogenic$V1 %in% closest_and_target$closest_SYMBOL) / length(monogenic$V1)
#22/34
monogenic[monogenic$V1 %in%  closest_and_target$closest_SYMBOL,]
