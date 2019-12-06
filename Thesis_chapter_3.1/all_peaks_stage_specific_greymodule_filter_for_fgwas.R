## subset ATAC peaks per stage for fgwas
library(readr)
stages <- c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

# After merging among stages with bedtools and gathering counts with featureCounts, I get the combined (consensus) set of peaks

basedir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/IDR_peaks/merged/"
outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/fGWAS/new_ATAC_peaks_2018/"
consensus = read_tsv(paste0(basedir,"ATAC_binary_conservative_normalQuality_narrowpeaks.bed"),col_names = T)


# 1 CPM-filtered set of peaks

cpm = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/counts/CPMs_atac-seq_1CPM_trim_conservativeCounts_normalquality.txt")
# total consensus peaks per stage

consensus = as.data.frame(consensus)
rownames(consensus) = consensus$Name
cpm_filt_consensus = consensus[rownames(cpm),]


## subset per stage
all = list()
for(s in stages){
  all[[s]] = cpm_filt_consensus[cpm_filt_consensus[s]==1,c("Chr","Start","End")]
  all[[s]]$stage =  rep(s,nrow(all[[s]]))
  
}
all = do.call("rbind",all)
write.table(
  all,
  paste0(outFolder,"all_ATAC_peaks_for_fgwas.bed"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)

## also stage-specific peaks
specific = list()

for(s in stages){
  specific[[s]] =   consensus[cpm_filt_consensus[s]==1 & rowSums(cpm_filt_consensus[5:12]) ==1,c("Chr","Start","End")]
  specific[[s]]$stage =  rep(s,nrow(specific[[s]]))
  
}
specific = do.call("rbind",specific)
write.table(
  specific,
  paste0(outFolder,"stage_specific_ATAC_peaks_for_fgwas.bed"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)

## remove those peaks with low variability accross stages
# for example those in the grey WGCNA module
WGCNA = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/WGCNA/ATAC_WGCNA_P12_S120_deepSplit2_merge0.20_signedHybrid_rmOutliersT_june2019_1CPM/WGCNA.gene2module.30M.txt")
WGCNA = WGCNA[WGCNA$V2=="grey",]
cpm_filt_WGCNA = cpm_filt_consensus[!(cpm_filt_consensus$Name %in% WGCNA$V1),]

## subset per stage
all = list()
for(s in stages){
  all[[s]] = cpm_filt_WGCNA[cpm_filt_WGCNA[s]==1,c("Chr","Start","End")]
  all[[s]]$stage =  rep(s,nrow(all[[s]]))
  
}
all = do.call("rbind",all)
write.table(
  all,
  paste0(outFolder,"all_ATAC_peaks_for_fgwas_grey_module_filtered.bed"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)

## also stage-specific peaks
specific = list()

for(s in stages){
  specific[[s]] =   consensus[cpm_filt_WGCNA[s]==1 & rowSums(cpm_filt_WGCNA[5:12]) ==1,c("Chr","Start","End")]
  specific[[s]]$stage =  rep(s,nrow(specific[[s]]))
  
}
specific = do.call("rbind",specific)
write.table(
  specific,
  paste0(outFolder,"stage_specific_ATAC_peaks_for_fgwas_grey_module_filtered.bed"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)
