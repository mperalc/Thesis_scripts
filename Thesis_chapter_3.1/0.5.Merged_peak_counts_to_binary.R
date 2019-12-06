# from merged peaks to fGWAS input
source("scripts/peaks_merged_to_binary_df.R")
stages = c("iPSC", "DE", "GT", "PF", "PE", "EP", "EN", "BLC")


#################### from peaks directly####################################
####################### read in the file with all peaks called (for all samples), and stages that have them#############
peaks = read.table(file = snakemake@input[["merged"]],
                   header = F,
                   check.names = F) # file with info of start and end of peaks, and samples that have it
colnames(peaks) = c("Chr", "Start", "End", "Peaks_present")  # add header
peaks$Name = paste("peak", peaks$Chr, peaks$Start, sep = "_") # add name column needed for the binary function
peaks = peaks[c(5, 1:4)] # reorder

binary = peaks_merged_to_binary_df(peaks) # this function needs columns "Name","Chr","Start","End","Peaks_present", in that order

# sort columns
binary = binary[,c("Name","Chr","Start","End","iPSC","DE","GT","PF","PE","EP","EN","BLC")]
write.table(
    binary,
    file = snakemake@output[["binary_ATAC"]],
    sep = "\t",
    row.names = F,
    col.names = T,
    quote = F
  )
