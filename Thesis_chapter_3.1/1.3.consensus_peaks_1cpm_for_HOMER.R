### save stage-specific peaks and peaks per stage (after removing peaks shared across all stages) for HOMER

# format (no header)
# id chr start end .
# one file per stage

library(readr)

stages <- c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

# After merging among stages with bedtools and gathering counts with featureCounts, I get the combined (consensus) set of peaks

basedir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/IDR_peaks/merged/"

consensus = read_tsv(paste0(basedir,"ATAC_binary_conservative_normalQuality_narrowpeaks.bed"),col_names = T)


# 1 CPM-filtered set of peaks

cpm = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/counts/CPMs_atac-seq_1CPM_trim_conservativeCounts_normalquality.txt")
# total consensus peaks per stage
colSums(consensus[c(5:12)])
mean(colSums(consensus[c(5:12)]))


consensus = as.data.frame(consensus)
rownames(consensus) = consensus$Name
cpm_filt_consensus = consensus[rownames(cpm),]

## remove peaks that are present in all stages

cpm_filt_consensus = cpm_filt_consensus[rowSums(cpm_filt_consensus[c(5:ncol(cpm_filt_consensus))])<8,]
nrow(cpm_filt_consensus)
# A total of 112,428 peaks

# save per stage
subset = list()

for(s in stages) {
  subset[[s]] = cpm_filt_consensus[cpm_filt_consensus[s] == 1,c("Name","Chr","Start","End")]
  subset[[s]]$dot = rep(".",nrow(subset[[s]]))
  write.table(subset[[s]],file = paste0("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/HOMER_results_per_stage/input/not_shared_all/",s,"_HOMER_ATAC.bed"),
              col.names = F, row.names = F, quote = F,sep = "\t")
  
}

# save only stage-specific peaks

subset = list()

for(s in stages) {
  subset[[s]] = cpm_filt_consensus[cpm_filt_consensus[s] == 1 & rowSums(cpm_filt_consensus[c(5:ncol(cpm_filt_consensus))])<2,c("Name","Chr","Start","End")]
  subset[[s]]$dot = rep(".",nrow(subset[[s]]))
  write.table(subset[[s]],file = paste0("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/HOMER_results_per_stage/input/stage-specific/",s,"_HOMER_ATAC.bed"),
              col.names = F, row.names = F, quote = F, sep = "\t")
  
}
