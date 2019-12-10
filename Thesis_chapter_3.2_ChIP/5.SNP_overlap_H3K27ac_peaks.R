# saving table with information of H3K27ac peaks that overlap SNPs of the credible set

library(GenomicRanges)
library(dplyr)
library(summarytools)
options(stringsAsFactors = FALSE)
# folder to save results
Chip_dir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ChIP-seq/SNPs_in_H3K27ac_peaks/"


### check loci classification from Anubha - there are some SNPs shared between loci

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

# H3K27ac peaks 
peaks = read.table(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks/merged/H3K27ac_binary_optimal.bed",
  header = T
  )

nrow(peaks)

# trim down with count file
trim = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks/trimmed/CPMs_H3K27ac_1CPM_trim_optimal.txt")
peaks = peaks[peaks$Name %in% rownames(trim),]
rm(trim)

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
# 5246 peaks in total in all credible sets.

temp = do.call("rbind", peaks_per_credset)
length(unique(temp$PeakID)) # 4399. So there are duplicates
sum(duplicated(temp$PeakID))
dups = temp[duplicated(temp$PeakID), "PeakID"] # list of duplicated ids

temp = temp[which(temp$PeakID %in% dups),]


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
            paste(Chip_dir,"credset_SNPs_in_H3K27ac_peaks_1CPM_PPA.txt", sep = ""),
            sep = "\t", quote = F, row.names = F, col.names = T)
## write table that has only peaks overlapping SNPs
peaks_overlap_SNP = merged_table[!is.na(merged_table$Chr.y),]
peaks_overlap_SNP$PeakID = paste("peak",peaks_overlap_SNP$Chr.y,peaks_overlap_SNP$Start,sep = "_")

## percentage of peaks overlapping SNPs
(length(unique(peaks_overlap_SNP$PeakID ))/nrow(peaks))*100
# 1.6% of peaks overlap SNPs, 4399

## snps overlapping peaks
(length(unique(peaks_overlap_SNP$IndexSNP))/length(unique(merged_table$IndexSNP)))*100
#8% of SNPs
# 8841
## snps with PPA over or equal 10%
sum(peaks_overlap_SNP$PPAg>=0.1)
## snps with PPA over or equal 50%
sum(peaks_overlap_SNP$PPAg>=0.5)

## merge with HOMER TF binding info
overlaps = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/motif_overlap_HRC_credset_2018/credible_set_HRC_SNP_TFBS_p95.txt",
                      sep = "\t",header = T)

peaks_overlap_SNP$IndexSNP = gsub("_",":",peaks_overlap_SNP$IndexSNP)
colnames(overlaps)[1] = "IndexSNP"
peaks_overlap_SNP = merge(peaks_overlap_SNP,overlaps,all.x = T, by = "IndexSNP")

## merge with rsID info
rsid = read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/HRC_credset.snp_ann.txt",
                  header = F)
colnames(rsid)=  c("rsID","IndexSNP","signal")
peaks_overlap_SNP = merge(peaks_overlap_SNP,rsid, all.x = T)

## locus annotation
peaks_overlap_SNP$signal = gsub("_.*","",peaks_overlap_SNP$signal)



write.table(peaks_overlap_SNP,
            paste(Chip_dir,"credset_SNPs_in_H3K27ac_peaks_1CPM_PPA_NoNA.txt", sep = ""),
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(na.omit(peaks_overlap_SNP$rsID),
            paste(Chip_dir,"credset_SNPs_in_H3K27ac_peaks_1CPM_PPA_NoNA_rsid.txt", sep = ""),
            sep = "\t", quote = F, row.names = F, col.names = F)
