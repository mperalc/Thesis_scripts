# important script

library(GenomicRanges)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(gridGraphics)
library(dplyr)
source("/Users/Marta/Documents/WTCHG/R scripts/atac-seq/corWGCNAmodules.R")
#source("/Users/Marta/Documents/WTCHG/R scripts/atac-seq/corAndPvalue.R") # what happened???

options(stringsAsFactors = FALSE)

basedir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/IDR_peaks/merged/"
outdir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/"

consensus = read_tsv(paste0(basedir,"ATAC_binary_conservative_normalQuality_narrowpeaks.bed"),col_names = T)

# 1 CPM-filtered set of peaks

cpm = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/counts/CPMs_atac-seq_1CPM_trim_conservativeCounts_normalquality.txt")

consensus = as.data.frame(consensus)
rownames(consensus) = consensus$Name
peaks = consensus[rownames(cpm),]


credset_names = read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/credset/credset_names.txt",
                           header = F)

credset = list()

for (n in credset_names$V1) {
  credset[[n]] = read.table(
    file = paste(
      "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/credset/",
      n,
      sep = ""
    ),
    header = T
  ) # 380 99% credible sets T2D 2017 HRC (Mahajan et al 2018)
  
}


# allocate each peak in correct credible set
# some will be empty, others will have peak
# integrate peak info and credset SNPs with findOverlaps

peaks_table.gr  = makeGRangesFromDataFrame(peaks,keep.extra.columns = T)

peaks_per_credset = list()
# n="credible_set_Eur_EYA2_20_45317678.txt"
for (n in names(credset)) {
  credset_subset = credset[[n]]
  credset_subset = credset_subset[, c("Chr", "Pos")]
  colnames(credset_subset) = c("chr", "start")
  credset_subset$end = credset_subset$start + 1
  credset_subset$chr = paste("chr", credset_subset$chr, sep = "")
  credset_subset$chr = as.character(credset_subset$chr)
  #there are duplicates within same credible set loci??
  credset_subset = unique(credset_subset)
  rownames(credset_subset) = paste("snp",
                                   credset_subset$chr,
                                   credset_subset$start,
                                   credset_subset$end,
                                   sep = "_")
  credset.gr = makeGRangesFromDataFrame(credset_subset)
  
  # Get all peaks in module overlapping credible set SNPs :
  overlapPeak = subsetByOverlaps(peaks_table.gr, credset.gr)
  hits <-
    findOverlaps(peaks_table.gr, credset.gr) # to add SNP metadata to overlapPeak table
  rsid <-
    CharacterList(split(names(credset.gr)[subjectHits(hits)], queryHits(hits)))
  rsid = paste(rsid, collapse = ",")
  
  mcols(overlapPeak) <- DataFrame(mcols(overlapPeak), rsid)
  
  
  overlapPeak = as.data.frame(overlapPeak)
  overlapPeak = overlapPeak[, c(1:3, 6)]
  colnames(overlapPeak) = c("Chr", "Start", "End", "SNPpos")
  if (nrow(overlapPeak) > 0) {
    overlapPeak$PeakID = paste("peak", overlapPeak$Chr, overlapPeak$Start, sep = "_")
  }
  peaks_per_credset[[n]] = overlapPeak
  
}
rm(credset_subset, overlapPeak)

# Get number of overlapping peaks. 
lapply(peaks_per_credset, nrow)

# how many credsets don't have any overlapping peaks?
# number of credible sets w/o any peaks over their SNPs
peaks_per_credset %>%
  lapply(nrow) %>%
  (function(x)
    x == 0) %>%
  sum

nrow(credset_names) # number credible sets
# 104 of 380 credible sets don't have any overlapping peaks (check this number every time)
(104 / 380) * 100 # 27%

# 73% do have
peaks_per_credset %>%
  lapply(nrow) %>%
  unlist %>%
  sum
# 2376 peaks in total in all credible sets.


#### get PPA for all SNPs under peaks in every credible set
# list of credible set
# SNP overlap peak
# snp, ppa, peak
# if no peak overlapping - "none" under peak col
# cbind and save with row names
SNPs_in_peaks = list()

for (n in names(credset)) {
  
  credset_subset = credset[[n]]
  credset_subset = credset_subset[, c("Chr", "Pos")]
  colnames(credset_subset) = c("chr", "start")
  credset_subset$end = credset_subset$start + 1
  credset_subset$chr = paste("chr", credset_subset$chr, sep = "")
  credset_subset$chr = as.character(credset_subset$chr)
  credset_subset = unique(credset_subset)
  rownames(credset_subset) = paste("snp",
                                   credset_subset$chr,
                                   credset_subset$start,
                                   credset_subset$end,
                                   sep = "_")
  credset.gr = makeGRangesFromDataFrame(credset_subset)
  
  # Get credible set SNPs overlapping peaks in module :
  overlapSNP = subsetByOverlaps(credset.gr, peaks_table.gr)
  hits <-
    findOverlaps(credset.gr, peaks_table.gr) # to add SNP metadata to overlapPeak table
  peakid =  unlist(split(peaks_table.gr$Name[subjectHits(hits)], queryHits(hits)))
  
  
  mcols(overlapSNP) <- DataFrame(mcols(overlapSNP), peakid)
  
  overlapSNP = as.data.frame(overlapSNP)
  if (nrow(overlapSNP) == 0) {
    credset[[n]]$IndexSNP = paste(credset[[n]]$Chr, credset[[n]]$Pos, sep =
                                    "_")
    overlapSNP = credset[[n]]
    overlapSNP$PeakID = rep("none", nrow(overlapSNP))
    SNPs_in_peaks[[n]] = overlapSNP
  }
  else {
    overlapSNP = overlapSNP[, c(1:2, 6)]
    colnames(overlapSNP) = c("Chr", "Pos", "PeakID")
    
    overlapSNP$IndexSNP = paste(gsub("\\chr*", "", overlapSNP$Chr), overlapSNP$Pos, sep =
                                  "_")
    credset[[n]]$IndexSNP = paste(credset[[n]]$Chr, credset[[n]]$Pos, sep =
                                    "_")
    
    SNPs_in_peaks[[n]] = merge(credset[[n]], overlapSNP[,c(3,4)], by = "IndexSNP", all = T)
    SNPs_in_peaks[[n]][is.na(SNPs_in_peaks[[n]]$PeakID),"PeakID"] = "none"
  }
}
SNPs_in_peaks = do.call("rbind",SNPs_in_peaks)
SNPs_in_peaks$credset = rownames(SNPs_in_peaks)
# the numbered credset names in the SNPs_in_peaks table do NOT denote order of SNPs by PPA!

# Add data on peak per stage
colnames(peaks)[1] = "PeakID"

merged_table <- merge(SNPs_in_peaks, peaks, by.x = "PeakID", by.y = "PeakID",  all.x = T)
merged_table = merged_table[c(6,2,5,1,7:9,10:17)]

write.table(merged_table,
            paste(outdir,"credset_SNPs_in_ATAC_peaks_1CPM_PPA.txt", sep = ""),
            sep = "\t", quote = F, row.names = F, col.names = T)
## write table that has only peaks overlapping SNPs
peaks_overlap_SNP = merged_table[!is.na(merged_table$Chr.y),]
peaks_overlap_SNP$PeakID = paste("peak",peaks_overlap_SNP$Chr.y,peaks_overlap_SNP$Start,sep = "_")

## percentage of peaks overlapping SNPs
(length(unique(peaks_overlap_SNP$PeakID ))/nrow(peaks))*100
# 1.5% of peaks overlap SNPs, 1975

## snps overlapping peaks
(length(unique(peaks_overlap_SNP$IndexSNP))/length(unique(merged_table$IndexSNP)))*100
#3.6% of SNPs
# 3890
## snps with PPA over or equal 10%
sum(peaks_overlap_SNP$PPAg>=0.1)
## snps with PPA over or equal 50%
sum(peaks_overlap_SNP$PPAg>=0.5)

write.table(peaks_overlap_SNP,
            paste(outdir,"credset_SNPs_in_ATAC_peaks_1CPM_PPA_NoNA.txt", sep = ""),
            sep = "\t", quote = F, row.names = F, col.names = T)
