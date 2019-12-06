# Calculate TSS enrichment in IDR peaks

library(readr)
library(ATACseqQC)
library(Rsamtools)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicAlignments)

basedir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/bams_merged_from_reps/"

stages <- c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

bamfile=BamFile(paste0(basedir,"iPSC.bam"))
seqinfo(bamfile)
# How many reads to read
yieldSize(bamfile) <- 10000000
# yieldSize(bamfile) <- NA

# Test reading bam file
# open(bamfile)
# 
# seqlev <- "chr13" ## subsample data for quick run
# which <- as(seqinfo(Hsapiens)[seqlev], "GRanges")
# 
# bamTop <- scanBam(bamfile, which=which,what = scanBamWhat() )
# 
# close(bamfile)

## end of test

# read in with GAlignments
gal <- readGAlignments(bamfile, param=ScanBamParam(what="flag"))

# Promoter/Transcript body (PT) score
# coverage of promoter divided by the coverage of its transcript body. 
#PT score will show if the signal is enriched in promoters.

txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
pt <- PTscore(gal, txs)
plot(pt$log2meanCoverage, pt$PT_score, 
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript")

# No idea what values should be

## Transcription Start Site (TSS) Enrichment Score
# ratio between aggregate distribution of reads centered on TSSs and that flanking 
# the corresponding TSSs. TSS score = the depth of TSS (1000 bp each side) / the depth 
# of end flanks (100bp each end). TSS enrichment score is calculated according to the 
# definition at https://www.encodeproject.org/data-standards/terms/#enrichment. 

tsse <- TSSEscore(gal, txs)
summary(tsse$TSS.enrichment.score)



peaks_df = do.call('rbind', peaks)

# Merge  peaks within stage to remove overlapping with different summits
peaks_subset = peaks_df[c(1:3,11:12)]
colnames(peaks_subset) = c("chr","start","end","sampid","names")
reduced_peaks = list()
for (s in stages) {
  reduced_peaks[[s]] = makeGRangesFromDataFrame(peaks_subset[peaks_subset$sampid ==s,])
  reduced_peaks[[s]] = reduce(reduced_peaks[[s]])
  reduced_peaks[[s]] = as.data.frame(reduced_peaks[[s]])
  reduced_peaks[[s]]$sampid = rep(s,nrow(reduced_peaks[[s]]))
  reduced_peaks[[s]]$names= paste(reduced_peaks[[s]]$seqnames,
                                  reduced_peaks[[s]]$start,
                                  reduced_peaks[[s]]$end, sep="_")
}


reduced_peaks = do.call('rbind', reduced_peaks)

# Number of non-overlapping peaks WITHIN SAME STAGE
#
freq(reduced_peaks$sampid)


# After merging among stages with bedtools and gathering counts with featureCounts, I get the combined (consensus) set of peaks

basedir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/"

consensus = read_tsv(paste0(basedir,"binary_conservative_normalQuality_narrowpeaks.bed"),col_names = T)

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
      sets = rev(stages), keep.order = TRUE)
    #  nintersects = 15)
# Only plotting top 15 intersects

dev.off()

# Unique peaks per stage

14096/54613 # iPSC
12209/53027 # DE
3968/49789 # GT
7252/64205 # PF
5540/57195 # PE
885/30406 # EP
8261/55803 # EN
1986/38941 # BLC

# mean of unique peaks per stage
mean(c(0.258107,0.2302412,0.07969632,0.1129507,0.0968616,0.0291061,0.1480386,0.05100023 ))

# Peak genomic annotations
peakAnnoConsensus <- annotatePeak(GRanges(as.data.frame(consensus)), tssRegion=c(-2000, 500),
                                  TxDb=txdb, annoDb="org.Hs.eg.db")

# Overall distribution
plotAnnoPie(peakAnnoConsensus)

plotDistToTSS(peakAnnoConsensus,
              title="Distribution of ATAC-seq peaks relative to TSS")

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
                       tssRegion=c(-500, 500), verbose=FALSE)

png(
  paste(basedir, '/GenomicAnnotationBar_consensus_set.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 6,
  height = 6
)

plotAnnoBar(peakAnnoList, title = "Distribution of ATAC-seq peaks relative to genomic features")

dev.off()

png(
  paste(basedir, '/TSSDistancenBar_consensus_set.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 6,
  height = 6
)

plotDistToTSS(peakAnnoList,
              title="Distribution of ATAC-seq peaks relative to TSS")
dev.off()
