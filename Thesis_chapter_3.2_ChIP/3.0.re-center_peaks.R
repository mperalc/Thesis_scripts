### Re-center peaks around summits
# stage-specific and all peaks
# for HOMER

library(readr)
library(GenomicRanges)
basedir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks/optimal/"
## re-centering around optimal peaks

currentDate <- Sys.Date() # to save date in name of output files


filelist <- list.files(basedir)
narrowPeaks = list()

for (f in filelist) {
  s = substring(f, 1, nchar(f)-27)
  narrowPeaks[[s]] = read_tsv(paste0(basedir,f),col_names = F)
  narrowPeaks[[s]]$sampid = rep(s,nrow(narrowPeaks[[s]]))
  narrowPeaks[[s]]$names= paste(narrowPeaks[[s]]$X1,narrowPeaks[[s]]$X2,narrowPeaks[[s]]$X3, sep="_")
  narrowPeaks[[s]] = narrowPeaks[[s]][c(1:3,10:12)]
  colnames(narrowPeaks[[s]]) = c("Chr","Start","End","Summit","sampid","names")
  # re-center +-50bp around summits
  narrowPeaks[[s]]$start = (narrowPeaks[[s]]$Start+narrowPeaks[[s]]$Summit)-50
  narrowPeaks[[s]]$end = (narrowPeaks[[s]]$Start+narrowPeaks[[s]]$Summit)+50
  narrowPeaks[[s]]$newnames = paste(narrowPeaks[[s]]$Chr,narrowPeaks[[s]]$start,narrowPeaks[[s]]$end, sep="_")
  narrowPeaks[[s]] = narrowPeaks[[s]][c("Chr","start","end","sampid","names","newnames")]
}


basedir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks/trimmed/"
outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ChIP-seq/HOMER/input/"

consensus = read_tsv("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks/merged/H3K27ac_binary_optimal.bed",col_names = T)
consensus = as.data.frame(consensus)
# those that pass the 1CPM filter and do not have sex chr
counts = read.table(paste0(basedir,"CPMs_H3K27ac_1CPM_trim_optimal.txt"))
consensus = consensus[consensus$Name %in% rownames(counts),]

stages <- c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

for(s in stages){
  subset = consensus[consensus[s]==1 & rowSums(consensus[5:12]) ==1,c("Name","Chr","Start","End")]
  
  ## get optimal peaks per stage that overlap the stage-specific consensus set peaks
  overlaps = findOverlaps(GRanges(narrowPeaks[[s]]),GRanges(subset))
  subset = subsetByOverlaps(GRanges(narrowPeaks[[s]]), GRanges(subset))
  subset = as.data.frame(subset)
  subset = subset[c("newnames","seqnames","start","end")]
  subset$dot =  rep(".",nrow(subset))
  write.table(
    subset,
    paste0(outFolder,s,"_stage_specific_for_HOMER_centered_summits_50bp.bed"),
    sep = "\t",
    row.names = F,
    col.names = F,
    quote = F
  )
}

## also all peaks per stage because 50% are unique anyways
## filter by CPM
for(s in stages){
  subset = consensus[consensus[s]==1,c("Name","Chr","Start","End")]
  ## get optimal peaks per stage that overlap the stage-specific consensus set peaks
  subset = subsetByOverlaps(GRanges(narrowPeaks[[s]]), GRanges(subset))
  subset = as.data.frame(subset)
  subset = subset[c("newnames","seqnames","start","end")]
  subset$dot =  rep(".",nrow(subset))
  write.table(
    subset,
    paste0(outFolder,s,"_all_peaks_for_HOMER_centered_summits_50bp.bed"),
    sep = "\t",
    row.names = F,
    col.names = F,
    quote = F
  )
}
