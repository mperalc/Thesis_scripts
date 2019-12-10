# ABC output to homer input per stage
# no header
# format: chr start end feature
library(GenomicRanges)
library(annotatr)
library(annotate)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

stages = c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

binary = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/for_homer/binary_predictions_ABC.txt",
                    header = T) # this function needs columns "Name","Chr","Start","End","Peaks_present", in that order


nrow(binary)
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

}


subset_for_factor_analysis_all = do.call("rbind",subset_for_factor_analysis )
nrow(subset_for_factor_analysis_all) # includes enhancers that map to multiple genes

# Save all enhancers per stage FROM original ABC peaks (not merged)
subset_for_factor_analysis_all$strand = rep(".",nrow(subset_for_factor_analysis_all))
for_homer_all = subset_for_factor_analysis_all[,c("name","chr","start","end","strand","cellType")]

for_homer_all = unique(for_homer_all) # because some map to multiple genes

nrow(for_homer_all) ## Number of  enhancers, for all stages, going into homer (all enhancers, including shared across stages)


# A bit more enhancers per stage than in the binary(merged), because some small enhancers might be merged if there is an overlapping
# enhancers from a different stage in between. Example:
########iPSC####                 #######iPSC######
               ########DE ########
for (s in stages) {
write.table(for_homer_all[for_homer_all$cellType == s,c("name","chr","start","end","strand")],
            paste0("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/for_homer/",s,"_homer_ABC_output_all_enhancers.bed"),
            sep="\t",row.names = F,col.names = F,quote = F)

}

# subset those that are not shared in all stages?

