# ABC output: significant, linked to expressed genes, in autosomal chromosomes, combined across all stages


stage = c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

merged = read.table(paste0("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/for_fGWAS/merged/merged_predictions_ABC.bed"),
           header = F, sep = "\t")

# make binary
peaks_merged_to_binary_df = function(peaks){
  if(ncol(peaks)>5){
    stop("Data frame with wrong number of columns")
  }
  # to get sample names for columns
  samples=paste(peaks$Peaks_present,collapse=",") 
  samples=unique(unlist(strsplit(samples, ",")))
  
  df=data.frame(matrix(ncol = length(samples)+4,nrow = nrow(peaks))) # make matrix to fill with n columns for samples 
  # + 4 more (Name of feature, chr, start, end)
  colnames(df)=c(c("Name","Chr","Start","End"),samples) # filling in new df
  df$Name=peaks$Name
  df$Chr=peaks$Chr
  df$Start=peaks$Start
  df$End=peaks$End
  df[,-1:-4] = 0  # all count columns to 0 (absence of peak)
  
  for (c in colnames(df[,-1:-4])){
    
    rows=peaks[grep(c,peaks$Peaks_present),"Name"] # peaks present in sample "c"
    df[which(df$Name %in% rows),which(colnames(df)==c)] = 1 # put 1 in column for "c" if peak is present in "c"
  }
  
  
  return(df)
}

colnames(merged) = c("Chr", "Start", "End", "Peaks_present")  # add header
merged = unique(merged)
merged$Name = paste("peak", merged$Chr, merged$Start, sep = "_") # add name column needed for the binary function

merged = merged[c("Name","Chr", "Start", "End","Peaks_present")]
binary = peaks_merged_to_binary_df(merged) # this function needs columns "Name","Chr","Start","End","Peaks_present", in that order
binary = binary[,c("Name","Chr","Start","End",stage)]
write.table(binary,"/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/for_fGWAS/binary_predictions_ABC.txt",sep = "\t",quote = F,
            row.names = F)
