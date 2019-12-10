# Calculating stats on number of enhancers per gene
# 
setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC")

library(readr)
library(summarytools)
# Differentiation stages

stage = c("iPSC", "DE", "GT", "PF", "PE", "EP", "EN")
# Load ABC model enhancer prediction results

pred_per_stage = list()
subset_for_factor_analysis = list()
for (s in stage) {
  pred_per_stage[[s]] = read_delim(file = paste0("Predictions/",s,"/EnhancerPredictions.txt"), delim = "\t")
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
                                                           "TargetGeneTSS",
                                                           "TargetGeneExpression",
                                                           "TargetGenePromoterActivityQuantile",
                                                           "hic.distance.adj",
                                                           "ABC.Score")]
  }

# Change header of ATAC and H3K27ac reps

subset_for_factor_analysis = do.call("rbind",subset_for_factor_analysis )

subset_for_factor_analysis$cellType = factor(subset_for_factor_analysis$cellType, 
                                             levels = stage, ordered = T)
# Enhancer sizes
subset_for_factor_analysis$enhancer_size = subset_for_factor_analysis$end - subset_for_factor_analysis$start
descr(subset_for_factor_analysis$enhancer_size, style = "rmarkdown")


# Number of enhancers that pass ABC score threshold
freq(subset_for_factor_analysis$cellType, style = "rmarkdown", cumul = FALSE, report.nas = F)

# Percentage of genic and intergenic enhancers
freq(subset_for_factor_analysis$class, style = "rmarkdown", cumul = FALSE, report.nas = F)

# Self promoter (T or F)
freq(subset_for_factor_analysis$isSelfPromoter, style = "rmarkdown", cumul = FALSE, report.nas = F)

# # Percentage of genic and intergenic enhancers
# print(ctable(subset_for_factor_analysis$cellType,subset_for_factor_analysis$class, prop= "r"),
#       method = "render")

# Distance to TSS
descr(subset_for_factor_analysis$distance, style = "rmarkdown")
# 
# # General summary
# dfSummary(subset_for_factor_analysis, plain.ascii = FALSE, style = "grid", 
#           graph.magnif = 0.75, valid.col = FALSE, tmp.img.dir = ".")

