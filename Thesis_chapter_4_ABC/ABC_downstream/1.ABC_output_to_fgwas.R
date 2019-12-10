# ABC output to fGWAS input per stage
# no header
# format: chr start end feature
library(GenomicRanges)
library(annotatr)
library(annotate)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

stages = c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

binary = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/for_fGWAS/binary_predictions_ABC.txt",
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
  subset_for_factor_analysis.gr[[s]] = GRanges(subset_for_factor_analysis[[s]])
  
}


subset_for_factor_analysis_all = do.call("rbind",subset_for_factor_analysis )
nrow(subset_for_factor_analysis_all) # includes enhancers that map to multiple genes

## Make list to genomic ranges
subset_for_factor_analysis.gr = GRangesList(subset_for_factor_analysis.gr)
subset_for_factor_analysis_all.gr = GRanges(subset_for_factor_analysis_all)

# Save all enhancers per stage FROM original ABC peaks (not merged)
for_fgwas_all = subset_for_factor_analysis_all[,c("chr","start","end","cellType")]
for_fgwas_all = unique(for_fgwas_all) # because some map to multiple genes
for_fgwas_all.gr = GRanges(for_fgwas_all)

nrow(for_fgwas_all) ## Number of  enhancers, for all stages, going into fgwas (all enhancers, including shared across stages)
table(for_fgwas_all$cellType)
# A bit more enhancers per stage than in the binary(merged), because some small enhancers might be merged if there is an overlapping
# enhancers from a different stage in between. Example:
########iPSC####                 #######iPSC######
               ########DE ########

write.table(for_fgwas_all,
            "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/for_fGWAS/fGWAS_ABC_output_all_enhancers.bed",
            sep="\t",row.names = F,col.names = F,quote = F)



# Removing those that are shared across all stages

keep = binary[rowSums(binary[5:ncol(binary)])<8,"Name"]
not_all = binary[binary$Name %in% keep,]
not_all.gr = GRanges(not_all)  # GRanges of binary peaks not shared at all stages

# Coparing originbal peaks and binary peaks not shared across all stages
# and keeping those in first that overlap the second

not_all_stages.gr = subsetByOverlaps(for_fgwas_all.gr,not_all.gr)
# save for fgwas
not_all_stages = as.data.frame(not_all_stages.gr)
write.table(not_all_stages[c("seqnames","start","end","cellType")],
            "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/for_fGWAS/fGWAS_ABC_output_enhancers_not_in_all_stages.bed",
            sep="\t",row.names = F,col.names = F,quote = F)


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
# nrow(subset_for_factor_analysis_all_closest_gene_match) / nrow(subset_for_factor_analysis_all_closest_gene)## percentage of total enhancers acting on closest gene
# 
# # Per stage
# nrow(subset_for_factor_analysis_all_closest_gene_match[subset_for_factor_analysis_all_closest_gene_match$cellType =="iPSC",])/
#   nrow(subset_for_factor_analysis_all_closest_gene[subset_for_factor_analysis_all_closest_gene$cellType == "iPSC",])
# nrow(subset_for_factor_analysis_all_closest_gene_match[subset_for_factor_analysis_all_closest_gene_match$cellType =="DE",])/
#   nrow(subset_for_factor_analysis_all_closest_gene[subset_for_factor_analysis_all_closest_gene$cellType == "DE",])
# nrow(subset_for_factor_analysis_all_closest_gene_match[subset_for_factor_analysis_all_closest_gene_match$cellType =="GT",])/
#   nrow(subset_for_factor_analysis_all_closest_gene[subset_for_factor_analysis_all_closest_gene$cellType == "GT",])
# nrow(subset_for_factor_analysis_all_closest_gene_match[subset_for_factor_analysis_all_closest_gene_match$cellType =="PF",])/
#   nrow(subset_for_factor_analysis_all_closest_gene[subset_for_factor_analysis_all_closest_gene$cellType == "PF",])
# nrow(subset_for_factor_analysis_all_closest_gene_match[subset_for_factor_analysis_all_closest_gene_match$cellType =="PE",])/
#   nrow(subset_for_factor_analysis_all_closest_gene[subset_for_factor_analysis_all_closest_gene$cellType == "PE",])
# nrow(subset_for_factor_analysis_all_closest_gene_match[subset_for_factor_analysis_all_closest_gene_match$cellType =="EP",])/
#   nrow(subset_for_factor_analysis_all_closest_gene[subset_for_factor_analysis_all_closest_gene$cellType == "EP",])
# nrow(subset_for_factor_analysis_all_closest_gene_match[subset_for_factor_analysis_all_closest_gene_match$cellType =="EN",])/
#   nrow(subset_for_factor_analysis_all_closest_gene[subset_for_factor_analysis_all_closest_gene$cellType == "EN",])
# nrow(subset_for_factor_analysis_all_closest_gene_match[subset_for_factor_analysis_all_closest_gene_match$cellType =="BLC",])/
#   nrow(subset_for_factor_analysis_all_closest_gene[subset_for_factor_analysis_all_closest_gene$cellType == "BLC",])
# mean(c(0.3409788,0.3344966,0.3451742,0.3922113, 0.3709272,0.2608624, 0.3374178,0.2811362))

write.table(unique(subset_for_factor_analysis_all_closest_gene_match[c("seqnames","start","end","cellType")]),
            "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/for_fGWAS/fGWAS_ABC_output_enhancers_match_closest_gene.bed",
            sep="\t",row.names = F,col.names = F,quote = F)

not_closest_gene_match = subset_for_factor_analysis_all_closest_gene[subset_for_factor_analysis_all_closest_gene$TargetGene != subset_for_factor_analysis_all_closest_gene$closest_SYMBOL,]
nrow(not_closest_gene_match) / nrow(subset_for_factor_analysis_all_closest_gene)## percentage of total enhancers not acting on closest gene

write.table(unique(not_closest_gene_match[c("seqnames","start","end","cellType")]),
            "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/for_fGWAS/fGWAS_ABC_output_enhancers_match_not_closest_gene.bed",
            sep="\t",row.names = F,col.names = F,quote = F)

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

## number of enhancers across all stages that match / don't match closest gene?
# 
binary.gr = GRanges(binary)
match = GRanges(na.omit(subset_for_factor_analysis_all_closest_gene_match[c(1:3)]))
nomatch = GRanges(na.omit(not_closest_gene_match[c(1:3)]))

match = as.data.frame(subsetByOverlaps(binary.gr, match))
nomatch = as.data.frame(subsetByOverlaps(binary.gr, nomatch))
nrow(match)
nrow(nomatch)
