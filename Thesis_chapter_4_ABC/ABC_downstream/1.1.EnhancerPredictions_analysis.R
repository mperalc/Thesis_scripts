# Calculating stats on number of enhancers per gene
# 
setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC")

library(readr)
library(summarytools)
library(UpSetR)
library(GenomicRanges)
library(ChIPseeker)
library(UpSetR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(ReactomePA)
library(clusterProfiler)
library(DOSE)
# Differentiation stages

stage = c("iPSC", "DE", "GT", "PF", "PE", "EP", "EN")
# Load ABC model enhancer prediction results

## For all prediction sas they come out of ABC
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
                                                           "TargetGeneExpression",
                                                           "TargetGeneTSS",
                                                           "TargetGeneExpression",
                                                           "TargetGenePromoterActivityQuantile",
                                                           "hic.distance.adj",
                                                           "ABC.Score")]
  }


subset_for_factor_analysis = do.call("rbind",subset_for_factor_analysis )

subset_for_factor_analysis = subset_for_factor_analysis[!is.na(subset_for_factor_analysis$TargetGeneExpression),]
subset_for_factor_analysis = subset_for_factor_analysis[subset_for_factor_analysis$chr %in% paste0("chr",1:22),]
nrow(subset_for_factor_analysis)

subset_for_factor_analysis$cellType = factor(subset_for_factor_analysis$cellType, 
                                             levels = stage, ordered = T)
# Enhancer sizes
subset_for_factor_analysis$enhancer_size = subset_for_factor_analysis$end - subset_for_factor_analysis$start
descr(subset_for_factor_analysis$enhancer_size, style = "rmarkdown")


# Number of enhancers that pass ABC score threshold
freq(subset_for_factor_analysis$cellType, style = "rmarkdown", cumul = FALSE, report.nas = F)

# Percentage of genic and intergenic enhancers
freq(subset_for_factor_analysis$class, style = "simple", cumul = FALSE, report.nas = F)

# Self promoter (T or F)
freq(subset_for_factor_analysis$isSelfPromoter, style = "rmarkdown", cumul = FALSE, report.nas = F)

# # Percentage of genic and intergenic enhancers
# print(ctable(subset_for_factor_analysis$cellType,subset_for_factor_analysis$class, prop= "r"),
#       method = "render")

# Enhancer names per stage. To see how many genes are linked per enhancer (each enhancer appears linked to one gene just once). 
for(s in stage){
frequencies = table(subset_for_factor_analysis[subset_for_factor_analysis$cellType == s,"name"])

frequencies = as.data.frame(frequencies)
print(descr(frequencies$Freq))
}
mean(c(2.17,2.71,1.99,1.87,2.03,2.06,2.05))
# Distance to TSS
descr(subset_for_factor_analysis$distance, style = "simple")
# 
# # General summary
# dfSummary(subset_for_factor_analysis, plain.ascii = FALSE, style = "grid", 
#           graph.magnif = 0.75, valid.col = FALSE, tmp.img.dir = ".")

# for sex chr, non-expressed target filtered predictions

binary = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/for_fGWAS/binary_predictions_ABC.txt",
                    header = T) # this function needs columns "Name","Chr","Start","End","Peaks_present", in that order

png(
  paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions", '/upset_consensus_set_top15_enhancers.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 12,
  height = 8
)

upset(as.data.frame(binary[,c(5:12)]), order.by = "freq",nsets = 8,point.size = 2, line.size = 1, 
      mainbar.y.label = "Number of active enhancers", sets.x.label = "Active enhancers per stage", 
      text.scale = c(2, 2, 2, 2, 2, 2),
      sets = rev(stages), keep.order = TRUE,
      nintersects = 15)
# Only plotting top 15 intersects

dev.off()

# all active enhancers per stage
sum(binary$iPSC)
sum(binary$DE)
sum(binary$GT) 
sum(binary$PF)
sum(binary$PE)
sum(binary$EP)
sum(binary$EN)
sum(binary$BLC)
# mean of active enhancers per stage
mean(c(sum(binary$iPSC),sum(binary$DE),
       sum(binary$GT) ,sum(binary$PF),
       sum(binary$PE),sum(binary$EP),
       sum(binary$EN),sum(binary$BLC)  ))

# Unique active enhancers per stage

4097/sum(binary$iPSC) # iPSC
3123/sum(binary$DE) # DE
1382/sum(binary$GT)  # GT
1309/sum(binary$PF)  # PF
1074/sum(binary$PE)  # PE
1198/sum(binary$EP)  # EP
1202/sum(binary$EN)  # EN
1198/sum(binary$BLC)  # BLC

# mean of unique active enhancers per stage
mean(c(4097/sum(binary$iPSC),3123/sum(binary$DE),
       1382/sum(binary$GT) ,1309/sum(binary$PF),
       1074/sum(binary$PE),1198/sum(binary$EP),
       1202/sum(binary$EN),1198/sum(binary$BLC)  ))


# Now get annotation stats for closest gene AND for linked gene (distance stats)
# And the usual for the genomic feature

# I have checked and the UCSC and knowgene give slightly different annotations for some genes, mapping to promoter when they arent
refseq <- makeTxDbFromUCSC(genome="hg19", tablename="refGene")

per_stage = list()


for (s in stages) {
  subset_c = binary[,c("Name","Chr","Start","End",s)]
  per_stage[[s]] = subset_c[subset_c[s]==1,]
  per_stage[[s]] = GRanges(per_stage[[s]])
}
per_stage_GRL = GRangesList(per_stage)

peakAnnoList <- lapply(per_stage_GRL, annotatePeak, TxDb=txdb,
                       tssRegion=c(-2000, 500), verbose=FALSE)
peakAnnoList_refseq = lapply(per_stage_GRL, annotatePeak, TxDb=refseq,
                             tssRegion=c(-2000, 500), verbose=FALSE)
# stats = list()
# for(s in stages){
#   stats[[s]] = peakAnnoList[[s]]@annoStat
#   stats[[s]]$stage = rep(s,nrow(stats[[s]]))
# }
# stats = do.call("rbind",stats)
# descr(stats[stats$Feature == "Promoter (<=1kb)","Frequency"])
# descr(stats[stats$Feature == "Promoter (1-2kb)","Frequency"])
# descr(stats[stats$Feature == "Distal Intergenic","Frequency"])
# descr(stats[stats$Feature == "1st Intron","Frequency"])
# descr(stats[stats$Feature == "Other Intron","Frequency"])


# 
# png(
#   paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions", '/GenomicAnnotationBar_predictions.png', sep = ''),units = "in",res = 400,type = "cairo",
#   width = 6,
#   height = 6
# )
# 
# plotAnnoBar(peakAnnoList, title = "Distribution of active enhancers relative to genomic features")
# plotAnnoBar(peakAnnoList_refseq, title = "Distribution of active enhancers relative to genomic features")
# 
# # There are some peaks marked as promoters, but linked to a closer promoter peak
# dev.off()

png(
  paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions", '/TSSDistancenBar_closest_gene_predictions.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 6,
  height = 6
)


plotDistToTSS(peakAnnoList_refseq,
              title="Distribution of ATAC-seq peaks relative to TSS")
dev.off()

p=as.data.frame(p$data)
descr(p[p$Feature=="0-1kb" & p$sign==1,"freq"])
descr(p[p$Feature=="1-3kb" & p$sign==1,"freq"])
descr(p[p$Feature=="5-10kb" & p$sign==1,"freq"])
descr(p[p$Feature=="10-100kb" & p$sign==1,"freq"])
descr(p[p$Feature==">100kb" & p$sign==1,"freq"])

