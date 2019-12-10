# make files for HOMER
# one file per stage and category
# id  chr  start  end  .
# no header

library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/"

## LMR

## DMR

stages = c("DE","GT","PF","PE","EP","EN","BLC")

annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)


comparisons = c(paste0(stages," - iPSC_DMRs_hypermethylated"),"islet - BLC_DMRs_hypermethylated")
comparisons = c(comparisons,paste0(stages," - iPSC_DMRs_hypomethylated"),"islet - BLC_DMRs_hypomethylated")

DMR = list()
for(c in comparisons){
  DMR[[c]] = read.csv(paste(outFolder,"DMR/",c,".csv",sep=""))
  DMR[[c]]$label = rep(c,nrow(DMR[[c]]))
  DMR[[c]]$label = gsub(" ", "", DMR[[c]]$label, fixed = TRUE)
  DMR[[c]]$id = paste(DMR[[c]]$seqnames, DMR[[c]]$start,DMR[[c]]$end,sep = "-")
  DMR[[c]]$dot = rep(".",nrow(DMR[[c]]))
  DMR[[c]] = DMR[[c]][c("id","seqnames","start","end","dot")]
  write.table(
    DMR[[c]],
    paste0(outFolder,c,"_for_HOMER.bed"),
    sep = "\t",
    row.names = F,
    col.names = F,
    quote = F
  )
}

stages = c("iPSC","DE","GT","PF","PE","EP","EN","BLC")
LMR = read.table(paste0(outFolder,"LMR_HMR_0.5_0.5/LMR_for_fgwas.bed"))
for(s in stages){
  subsetdf = LMR[LMR$V4==s,]
  subsetdf$id = paste(subsetdf$V1, subsetdf$V2,subsetdf$V3,sep = "-")
  subsetdf$dot = rep(".",nrow(subsetdf))
  subsetdf = subsetdf[c("id","V1","V2","V3","dot")]
  write.table(
    subsetdf,
    paste0(outFolder,s,"_LMR_for_HOMER.bed"),
    sep = "\t",
    row.names = F,
    col.names = F,
    quote = F
  )
}
