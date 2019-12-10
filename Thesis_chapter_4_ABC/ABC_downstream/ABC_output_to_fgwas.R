# ABC output to fGWAS input per stage
# no header
# format chr start end feature

stages = c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

binary = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/for_fGWAS/binary_predictions_ABC.txt",
                    header = T) # this function needs columns "Name","Chr","Start","End","Peaks_present", in that order


nrow(binary)

# save all enhancers per stage
peak_per_stage=list()

for(s in stages){ # get name of stage per peak
  peak_per_stage[[s]]=binary[which(binary[s]==1),c("Chr","Start","End")]
  peak_per_stage[[s]]$Stage=rep(s,nrow(peak_per_stage[[s]]))
}
fgwas_annotation=do.call("rbind",peak_per_stage) # bind all peaks

fgwas_annotation=fgwas_annotation[order(fgwas_annotation[,1],fgwas_annotation[,2]),]
fgwas_annotation=fgwas_annotation[!(fgwas_annotation$Chr=="chrY"|fgwas_annotation$Chr=="chrX"),]  # remove sex chromosomes
write.table(fgwas_annotation,
            "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/for_fGWAS/fGWAS_ABC_output_all_enhancers.bed",
            sep="\t",row.names = F,col.names = F,quote = F)



# Removing those that are shared across all stages

keep = binary[rowSums(binary[5:ncol(binary)])<8,"Name"]
not_all = binary[binary$Name %in% keep,]

peak_per_stage=list()

for(s in stages){ # get name of stage per peak
  peak_per_stage[[s]]=not_all[which(not_all[s]==1),c("Chr","Start","End")]
  peak_per_stage[[s]]$Stage=rep(s,nrow(peak_per_stage[[s]]))
}
fgwas_annotation=do.call("rbind",peak_per_stage) # bind all peaks
fgwas_annotation=fgwas_annotation[order(fgwas_annotation[,1],fgwas_annotation[,2]),]
fgwas_annotation=fgwas_annotation[!(fgwas_annotation$Chr=="chrY"|fgwas_annotation$Chr=="chrX"),]  # remove sex chromosomes
write.table(fgwas_annotation,
            "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/for_fGWAS/fGWAS_ABC_output_enhancers_not_in_all_stages.bed",
            sep="\t",row.names = F,col.names = F,quote = F)


# Those that are acting on closest gene?