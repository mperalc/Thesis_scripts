### H3K27ac genomic annotation consensus
# and upsetplot
library(readr)
library(UpSetR)
library(ChIPseeker)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <-TxDb.Hsapiens.UCSC.hg19.knownGene
library(GenomicAlignments)

basedir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks/merged/"

stages <- c("iPSC","DE","GT","PF","PE","EP","EN","BLC")


consensus = read_tsv(paste0(basedir,"H3K27ac_binary_optimal.bed"),col_names = T)

# total consensus peaks per stage
colSums(consensus[c(5:12)])
mean(colSums(consensus[c(5:12)]))

# Unique per stage
png(
  paste(basedir, '/upset_consensus_set.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 12,
  height = 8
)

upset(as.data.frame(consensus[,c(5:12)]), order.by = "freq",nsets = 8,point.size = 2, line.size = 1, 
      mainbar.y.label = "Number of peaks", sets.x.label = "Peaks per stage", 
      text.scale = c(2, 2, 2, 2, 2, 2),
      sets = rev(stages), keep.order = TRUE,
  nintersects = 15)
# Only plotting top 15 intersects

dev.off()

# Unique peaks per stage

22100/colSums(consensus[c(5:12)])[1] # iPSC
37406/colSums(consensus[c(5:12)])[2] # DE
40921/colSums(consensus[c(5:12)])[3] # GT
23813/colSums(consensus[c(5:12)])[4] # PF
22456/colSums(consensus[c(5:12)])[5] # PE
20552/colSums(consensus[c(5:12)])[6] # EP
21485/colSums(consensus[c(5:12)])[7] # EN
21348/colSums(consensus[c(5:12)])[8] # BLC

# mean of unique peaks per stage
mean(c(22100/colSums(consensus[c(5:12)])[1], # iPSC
       37406/colSums(consensus[c(5:12)])[2], # DE
       40921/colSums(consensus[c(5:12)])[3], # GT
       23813/colSums(consensus[c(5:12)])[4], # PF
       22456/colSums(consensus[c(5:12)])[5], # PE
       20552/colSums(consensus[c(5:12)])[6], # EP
       21485/colSums(consensus[c(5:12)])[7], # EN
       21348/colSums(consensus[c(5:12)])[8] )) # BLC

# Peak genomic annotations  ########################################################
## these are not the ones that will be included in the thesis. 
# They are only for comparison with the filtered consensus set


peakAnnoConsensus <- annotatePeak(GRanges(as.data.frame(consensus)), tssRegion=c(-2000, 500),
                                  TxDb=txdb, annoDb="org.Hs.eg.db")

# Overall distribution

plotDistToTSS(peakAnnoConsensus,
              title="Distribution of ChIP-seq peaks relative to TSS")
plotAnnoPie(peakAnnoConsensus)

# Distribution per stage

# Grab peaks per stage

consensus_per_stage = list()

consensus = as.data.frame(consensus)

for (s in stages) {
  subset_c = consensus[,c("Name","Chr","Start","End",s)]
  consensus_per_stage[[s]] = subset_c[subset_c[s]==1,]
  consensus_per_stage[[s]] = GRanges(consensus_per_stage[[s]])
}
consensus_per_stage_GRL = GRangesList(consensus_per_stage)

peakAnnoList <- lapply(consensus_per_stage_GRL, annotatePeak, TxDb=txdb,
                       tssRegion=c(-2000, 500), verbose=FALSE)

plotAnnoBar(peakAnnoList, title = "Distribution of ChIP-seq peaks relative to genomic features")


plotDistToTSS(peakAnnoList,
              title="Distribution of ChIP-seq peaks relative to TSS")

## these are not the ones that will be included in the thesis. 
# They are only for comparison with the filtered consensus set