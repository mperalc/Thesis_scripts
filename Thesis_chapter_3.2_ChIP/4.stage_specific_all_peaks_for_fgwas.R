# stage specific and all peaks for fgwas
library(readr)
# chr1    714177  714621  feature
# no header
stages <- c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

basedir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks/trimmed/"
outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ChIP-seq/fgwas/input/"

consensus = read_tsv("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks/merged/H3K27ac_binary_optimal.bed",col_names = T)
# those that pass the 1CPM filter and do not have sex chr
consensus = as.data.frame(consensus)

counts = read.table(paste0(basedir,"CPMs_H3K27ac_1CPM_trim_optimal.txt"))
consensus = consensus[consensus$Name %in% rownames(counts),]

subset = list()
for(s in stages){
  subset[[s]] = consensus[consensus[s]==1 & rowSums(consensus[5:12]) ==1,c("Chr","Start","End")]
  subset[[s]]$stage =  rep(s,nrow(subset[[s]]))
  
}
subset = do.call("rbind",subset)
write.table(
  subset,
  paste0(outFolder,"stage_specific_for_fgwas.bed"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)

## also all peaks per stage because 50% are unique anyways
subset = list()

for(s in stages){
  subset[[s]] = consensus[consensus[s]==1,c("Chr","Start","End")]
  subset[[s]]$stage =  rep(s,nrow(subset[[s]]))
  
}
subset = do.call("rbind",subset)
write.table(
  subset,
  paste0(outFolder,"all_peaks_for_fgwas.bed"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)