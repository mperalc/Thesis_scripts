## select stage-specific peaks for HOMER
library(readr)

stages <- c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

# After merging among stages with bedtools and gathering counts with featureCounts, I get the combined (consensus) set of peaks

basedir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks/trimmed/"
outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ChIP-seq/HOMER/input/"
consensus = read_tsv("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks/merged/H3K27ac_binary_optimal.bed",col_names = T)
consensus = as.data.frame(consensus)
# those that pass the 1CPM filter and do not have sex chr
counts = read.table(paste0(basedir,"CPMs_H3K27ac_1CPM_trim_optimal.txt"))
consensus = consensus[consensus$Name %in% rownames(counts),]

# HOMER input format
# chr13-113654392-113656424       chr13   113654392       113656424       .
# no header
for(s in stages){
  subset = consensus[consensus[s]==1 & rowSums(consensus[5:12]) ==1,c("Name","Chr","Start","End")]
  subset$dot =  rep(".",nrow(subset))
  write.table(
    subset,
    paste0(outFolder,s,"_stage_specific_for_HOMER.bed"),
    sep = "\t",
    row.names = F,
    col.names = F,
    quote = F
  )
}

## also all peaks per stage because 50% are unique anyways
for(s in stages){
  subset = consensus[consensus[s]==1,c("Name","Chr","Start","End")]
  subset$dot =  rep(".",nrow(subset))
  write.table(
    subset,
    paste0(outFolder,s,"_all_peaks_for_HOMER.bed"),
    sep = "\t",
    row.names = F,
    col.names = F,
    quote = F
  )
}
