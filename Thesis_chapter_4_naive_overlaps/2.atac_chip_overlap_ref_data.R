# Overlap H3K27ac and ATAC-seq and get

# Stats of overlaps at enhancers and promoter locations for ATAC-seq and H3K27ac and compare to publicly available data


library(GenomicRanges)


folder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/ENCODE/processed/"
subfolders = c("A549","Adult_adrenal_gland","Adult_breast_epithelium","Adult_pancreas",
               "Adult_stomach","Adult_thyroid")
# Initialize lists
ATAC_list = list()
H3K27ac_list = list()

#ATAC
for (s in subfolders) {
  files = list.files(path = paste(folder, s, "/ATAC/bed_peaks", sep = ""))
  print(paste(s, "with files",files))
  for (f in files) {
    ATAC_list[[f]] = read.table(paste(folder, s, "/ATAC/bed_peaks/",f, sep = ""),
                                header = F,skip = 1)
    # skip header in gzipped, and track name in unzipped
    
  }
}

#ChIP
for (s in subfolders) {
  files = list.files(path = paste(folder, s, "/H3K27ac", sep = ""))
  for (f in files) {
    H3K27ac_list[[f]] = read.table(paste(folder, s, "/H3K27ac/",f, sep = ""),
                                   header = F,skip = 1)
    # skip header in gzipped, and track name in unzipped
    
  }
}

# Filter 1st column (Chr) so that it only contains chr 1=22 
# Make sure both ATAC and ChIP have similar filters for this comparison
chrs = paste('chr',rep(c(1:22),1),sep = "")

for (n in 1:10) {
  ATAC_list[[n]]$V1 = as.character(ATAC_list[[n]]$V1)
  H3K27ac_list[[n]]$V1 = as.character(H3K27ac_list[[n]]$V1)
  
  ATAC_list[[n]] =  ATAC_list[[n]][ATAC_list[[n]]$V1 %in% chrs,]
  H3K27ac_list[[n]] = H3K27ac_list[[n]][H3K27ac_list[[n]]$V1 %in% chrs,]
}



# Make granges objects
for (n in 1:10) {
  ATAC_list[[n]] = makeGRangesFromDataFrame(
    df = ATAC_list[[n]],
    seqnames.field = "V1",
    start.field = "V2",
    end.field = "V3",
    keep.extra.columns = T
  )
  H3K27ac_list[[n]] = makeGRangesFromDataFrame(
    df = H3K27ac_list[[n]],
    seqnames.field = "V1",
    start.field = "V2",
    end.field = "V3",
    keep.extra.columns = T
  )
}

lt = list()
chip_in_atac= list()
atac_in_chip = list()
df <- data.frame(total_peaks_atac=numeric(1),
                 unique_overlaps_atac_in_chip=character(1), 
                 atac_in_chip=numeric(1),
                 percentage_atac_in_chip = numeric(1),
                 total_peaks_chip=numeric(1),
                 unique_overlaps_chip_in_atac=character(1), 
                 percentage_chip_in_atac = numeric(1),
                 stringsAsFactors=FALSE) 

for (n in 1:10){
  peaks.gr = ATAC_list[[n]]
  chip.gr = H3K27ac_list[[n]]
  
  
  # ATAC in ChIP
  atac_in_chip[[n]] = subsetByOverlaps(peaks.gr, chip.gr)
  overlap = countOverlaps(peaks.gr, chip.gr)
  df$total_peaks_atac = length(overlap) # how many atac peaks in sample
  df$unique_overlaps_atac_in_chip = paste(unique(overlap),collapse = ",") # there are some atac peaks that overlap more than one chip peak
  df$atac_in_chip = length(overlap) - sum(overlap == 0) # atac peaks in chip peaks
  df$percentage_atac_in_chip = (df$atac_in_chip/df$total_peaks_atac)*100
  # includes atac peaks with more than one match in chip peaks (e.g. large atac peak overlapping small chip peaks)
  
  # chip in atac
  chip_in_atac[[s]] = subsetByOverlaps(chip.gr,peaks.gr)
  
  overlap = countOverlaps(chip.gr,peaks.gr)
  df$total_peaks_chip=length(overlap) # how many atac peaks in stage s 
  df$unique_overlaps_chip_in_atac=paste(unique(overlap),collapse=",") # there are some atac peaks that overlap more than one chip peak
  df$chip_in_atac = length(overlap) - sum(overlap == 0) # atac peaks in chip peaks
  df$percentage_chip_in_atac = (df$chip_in_atac/  df$total_peaks_chip)*100
  lt[[n]] = df
}
lt_df = do.call("rbind",lt)
lt_df
lt_df$Name = c("A549","A549","A549","Adult_adrenal_gland","Adult_breast_epithelium",
               "Adult_pancreas","Adult_pancreas","Adult_pancreas",
               "Adult_stomach","Adult_thyroid")

## Islets


folder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/"
subfolders = c("Islets_Agata")
# Initialize lists
ATAC_list = list()
H3K27ac_list = list()


#ATAC
for (s in subfolders) {
  files = list.files(path = paste(folder, s, "/ATAC", sep = ""))
  print(paste(s, "with files",files))
  for (f in files) {
    ATAC_list[[f]] = read.table(paste(folder, s, "/ATAC/",f, sep = ""),
                                header = F)
    # skip header in gzipped, and track name in unzipped
    
  }
}

#ChIP
for (s in subfolders) {
  files = list.files(path = paste(folder, s, "/H3K27ac", sep = ""))
  for (f in files) {
    H3K27ac_list[[f]] = read.table(paste(folder, s, "/H3K27ac/",f, sep = ""),
                                   header = F)
    # skip header in gzipped, and track name in unzipped
  }
}

# Filter 1st column (Chr) so that it only contains chr 1-22
# Make sure both ATAC and ChIP have similar filters for this comparison
chrs = paste('chr',rep(c(1:22),1),sep = "")

n=1
ATAC_list[[n]]$V1 = as.character(ATAC_list[[n]]$V1)
H3K27ac_list[[n]]$V1 = as.character(H3K27ac_list[[n]]$V1)

ATAC_list[[n]] =  ATAC_list[[n]][ATAC_list[[n]]$V1 %in% chrs,]
H3K27ac_list[[n]] = H3K27ac_list[[n]][H3K27ac_list[[n]]$V1 %in% chrs,]

# Make granges objects

ATAC_list[[n]] = makeGRangesFromDataFrame(
  df = ATAC_list[[n]],
  seqnames.field = "V1",
  start.field = "V2",
  end.field = "V3",
  keep.extra.columns = T
)
H3K27ac_list[[n]] = makeGRangesFromDataFrame(
  df = H3K27ac_list[[n]],
  seqnames.field = "V1",
  start.field = "V2",
  end.field = "V3",
  keep.extra.columns = T
)


lt = list()
chip_in_atac= list()
atac_in_chip = list()
df <- data.frame(total_peaks_atac=numeric(1),
                 unique_overlaps_atac_in_chip=character(1), 
                 atac_in_chip=numeric(1),
                 percentage_atac_in_chip = numeric(1),
                 total_peaks_chip=numeric(1),
                 unique_overlaps_chip_in_atac=character(1), 
                 percentage_chip_in_atac = numeric(1),
                 stringsAsFactors=FALSE) 

peaks.gr = ATAC_list[[n]]
chip.gr = H3K27ac_list[[n]]


# ATAC in ChIP
atac_in_chip[[n]] = subsetByOverlaps(peaks.gr, chip.gr)
overlap = countOverlaps(peaks.gr, chip.gr)
df$total_peaks_atac = length(overlap) # how many atac peaks in sample
df$unique_overlaps_atac_in_chip = paste(unique(overlap),collapse = ",") # there are some atac peaks that overlap more than one chip peak
df$atac_in_chip = length(overlap) - sum(overlap == 0) # atac peaks in chip peaks
df$percentage_atac_in_chip = (df$atac_in_chip/df$total_peaks_atac)*100
# includes atac peaks with more than one match in chip peaks (e.g. large atac peak overlapping small chip peaks)

# chip in atac
chip_in_atac[[s]] = subsetByOverlaps(chip.gr,peaks.gr)

overlap = countOverlaps(chip.gr,peaks.gr)
df$total_peaks_chip=length(overlap) # how many atac peaks in stage s 
df$unique_overlaps_chip_in_atac=paste(unique(overlap),collapse=",") # there are some atac peaks that overlap more than one chip peak
df$chip_in_atac = length(overlap) - sum(overlap == 0) # atac peaks in chip peaks
df$percentage_chip_in_atac = (df$chip_in_atac/  df$total_peaks_chip)*100
lt[[n]] = df

lt_df2 = do.call("rbind",lt)
lt_df2$Name = "islets"
rbind(lt_df,lt_df2)

## Add to dataframe the comparison of IDR ATAC and H3K27ac of my data per stage   # #######################

library(GenomicRanges)
library(UpSetR)
library(data.table)

#ATAC
atac = read.table( "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/IDR_peaks/merged/ATAC_binary_conservative_normalQuality_narrowpeaks.bed",
                   header = T
)
# filter so that it only contains autosomes
atac = atac[atac$Chr %in% chrs,]
peaks.gr = makeGRangesFromDataFrame(atac,keep.extra.columns = T)

#ChIP
datChip = read.table(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks/merged/H3K27ac_binary_optimal.bed",header = T
)
# reorder
datChip = datChip[,c("Name","Chr","Start","End","iPSC","DE","GT","PF","PE","EP","EN","BLC")]
datChip = datChip[datChip$Chr %in% chrs,]
chip.gr = makeGRangesFromDataFrame(
  df = datChip,
  keep.extra.columns = T
)
lt = list()
chip_in_atac= list()
atac_in_chip = list()
df <- data.frame(total_peaks_atac=numeric(1),
                 unique_overlaps_atac_in_chip=character(1), 
                 atac_in_chip=numeric(1),
                 percentage_atac_in_chip = numeric(1),
                 total_peaks_chip=numeric(1),
                 unique_overlaps_chip_in_atac=character(1), 
                 percentage_chip_in_atac = numeric(1),
                 stringsAsFactors=FALSE) 
stages=c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

for (s in stages){
  atac_subset = atac[atac[s] == 1,]
  peaks.gr = makeGRangesFromDataFrame(atac_subset,keep.extra.columns = T)
  
  datChip_subset = datChip[datChip[s] == 1,]
  
  chip.gr = makeGRangesFromDataFrame(df = datChip_subset,
                                     keep.extra.columns = T)
  
  # ATAC in ChIP
  atac_in_chip[[s]] = subsetByOverlaps(peaks.gr, chip.gr)
  
  overlap = countOverlaps(peaks.gr, chip.gr)
  df$total_peaks_atac=length(overlap) # how many atac peaks in stage s 
  df$unique_overlaps_atac_in_chip=paste(unique(overlap),collapse=",") # there are some atac peaks that overlap more than one chip peak
  df$atac_in_chip = length(overlap) - sum(overlap == 0) # atac peaks in chip peaks
  df$percentage_atac_in_chip = (df$atac_in_chip/df$total_peaks_atac)*100
  # includes atac peaks with more than one match in chip peaks (e.g. large atac peak overlapping small chip peaks)
  
  # chip in atac
  chip_in_atac[[s]] = subsetByOverlaps(chip.gr,peaks.gr)
  
  overlap = countOverlaps(chip.gr,peaks.gr)
  df$total_peaks_chip=length(overlap) # how many atac peaks in stage s 
  df$unique_overlaps_chip_in_atac=paste(unique(overlap),collapse=",") # there are some atac peaks that overlap more than one chip peak
  df$chip_in_atac = length(overlap) - sum(overlap == 0) # atac peaks in chip peaks
  df$percentage_chip_in_atac = (df$chip_in_atac/  df$total_peaks_chip)*100
  lt[[s]] = df
}
lt_df3 = do.call("rbind",lt)
lt_df3$Name = stages
rbind(lt_df,lt_df2,lt_df3)

write.table(rbind(lt_df,lt_df2,lt_df3),"/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/ATAC_H3K27ac_overlaps.txt",
            row.names = F,
            col.names = T,
            quote = F,
            sep = "\t")


