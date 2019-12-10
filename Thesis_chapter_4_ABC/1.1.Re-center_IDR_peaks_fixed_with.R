# Re-size IDR peak files
# Output is saf file for featureCounts

input_folder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/IDR_peaks/"
output_folder = input_folder
resize = 250 # bp around summit to resize
###
stages = c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

peaks = list()
for(s in stages) {
  peaks[[s]] = read.table(paste0(input_folder,s,".conservative_peak.narrowPeak"),
                     header = F,skip = 1)
  peaks[[s]] = peaks[[s]][c(1,2,3,10)]
  
  peaks[[s]]$summit = peaks[[s]]$V2 + peaks[[s]]$V10
  peaks[[s]]$start = peaks[[s]]$summit - resize
  peaks[[s]]$end = peaks[[s]]$summit + resize
  peaks[[s]]$ID = paste(peaks[[s]]$V1,peaks[[s]]$start,sep = "_")
  peaks[[s]]$strand = rep(".",nrow(peaks[[s]]))
  peaks[[s]] = peaks[[s]][c("ID","V1","start","end","strand")]
  colnames(peaks[[s]]) = c("GeneID", "Chr", "Start", "End", "Strand")
  write.table(peaks[[s]],
              paste0(output_folder,s,".conservative_peak.saf"),
              col.names = T,
              row.names = F,
              sep = "\t",
              quote = F)
}

