### DMP, DMR, LMR, HMR to fgwas input ############
library(GenomicRanges)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(ggplot2)
library(reshape)
`%notin%` <- function(x,y) !(x %in% y)

outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/"

beta_lower_threshold = format(round(0.5, 2), nsmall = 1)
beta_higher_threshold = format(round(0.5, 2), nsmall = 1)

###############################################################################################
# get the latest EPIC annotation data
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

stages = c("DE","GT","PF","PE","EP","EN","BLC")


### DMP #################
comparisons = c(paste0(stages," - iPSC_DMPs_hypermethylated"),"islet - BLC_DMPs_hypermethylated")
comparisons = c(comparisons,paste0(stages," - iPSC_DMPs_hypomethylated"),"islet - BLC_DMPs_hypomethylated")

DMP = list()
for(c in comparisons){
DMP[[c]] = read.csv(paste(outFolder,"DMP/",c,".csv",sep=""))
rownames(DMP[[c]]) = DMP[[c]]$CpG
DMP[[c]] = merge(DMP[[c]],annEPIC, by = 0)
DMP[[c]]$label = rep(c,nrow(DMP[[c]]))
DMP[[c]]$label = gsub(" ", "", DMP[[c]]$label, fixed = TRUE)
}

DMP = do.call("rbind", DMP)
DMP = DMP[c("chr","pos","label")]
DMP$end = DMP$pos + 1
DMP = DMP[c(1,2,4,3)]
# sort
DMP = makeGRangesFromDataFrame(DMP, keep.extra.columns = T,start.field = "pos")
DMP = GenomicRanges::sort(DMP)
subset = as.data.frame(DMP)
write.table(
  subset[c(1:3,6)],
  paste0(outFolder,"DMP_for_fgwas.bed"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)

rm(subset,DMP)
## DMR ##############
####################

comparisons = c(paste0(stages," - iPSC_DMRs_hypermethylated"),"islet - BLC_DMRs_hypermethylated")
comparisons = c(comparisons,paste0(stages," - iPSC_DMRs_hypomethylated"),"islet - BLC_DMRs_hypomethylated")

DMR = list()
for(c in comparisons){
  DMR[[c]] = read.csv(paste(outFolder,"DMR/",c,".csv",sep=""))
  DMR[[c]]$label = rep(c,nrow(DMR[[c]]))
  DMR[[c]]$label = gsub(" ", "", DMR[[c]]$label, fixed = TRUE)
  
}

DMR = do.call("rbind", DMR)
DMR = DMR[c("seqnames","start","end","label")]
# sort
DMR = makeGRangesFromDataFrame(DMR, keep.extra.columns = T)
DMR = GenomicRanges::sort(DMR)
subset = as.data.frame(DMR)
write.table(
  subset[c(1:3,6)],
  paste0(outFolder,"DMR_for_fgwas.bed"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)
rm(DMR,subset)

### DMR divided by genomic feature

hyper = read.csv(paste(outFolder,"DMR/DMRs_hypermethylated_annotated_allstages.csv",sep=""))
hypo = read.csv(paste(outFolder,"DMR/DMRs_hypomethylated_annotated_allstages.csv",sep=""))

## plot abundance barplot
df = as.data.frame(cbind(table(hyper$contrast),table(hypo$contrast)))
df$contrast = names(table(hyper$contrast))
colnames(df) = c("hypermethylated","hypomethylated","contrast")
melted = melt(df,id.vars = "contrast")

# negative values for plotting
melted[melted$variable == "hypomethylated","value"] = -melted[melted$variable == "hypomethylated","value"]
colnames(melted) = c("contrast","DMR type", "Number of DMRs")
melted$contrast = factor(melted$contrast, levels = rev(unique(hypo$contrast)), ordered = T)


p = ggplot(melted, aes( x = contrast, y = `Number of DMRs`)) +  
  geom_bar(stat = "identity", aes(fill = `DMR type`)) + 
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(filename =  paste0(outFolder,"DMRs_barplot.png"),
      plot = p,device = "png",width = 6, height = 5, units = "in", dpi = 400)

####

# change useless long intron exon downstream names
hyper$annotation = gsub(".*Promoter .*", "Promoter", hyper$annotation)
hypo$annotation = gsub(".*Promoter .*", "Promoter", hypo$annotation )

hypo_not_prom = hypo[hypo$annotation %notin% "Promoter",]
hyper_not_prom = hyper[hyper$annotation %notin% "Promoter",]

hypo_prom = hypo[hypo$annotation %in% "Promoter",]
hyper_prom = hyper[hyper$annotation %in% "Promoter",]
  
hypo_intergenic = hypo[hypo$annotation == "Distal Intergenic" ,]
hyper_intergenic = hyper[hyper$annotation =="Distal Intergenic" ,]

write.table(
  hypo_not_prom[c(1:3,24)],
  paste0(outFolder,"DMR_for_fgwas_hypo_notpromoter.bed"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)

write.table(
  hyper_not_prom[c(1:3,24)],
  paste0(outFolder,"DMR_for_fgwas_hyper_notpromoter.bed"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)
write.table(
  hypo_prom[c(1:3,24)],
  paste0(outFolder,"DMR_for_fgwas_hypo_promoter.bed"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)

write.table(
  hyper_prom[c(1:3,24)],
  paste0(outFolder,"DMR_for_fgwas_hyper_promoter.bed"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)
write.table(
  hyper_intergenic[c(1:3,24)],
  paste0(outFolder,"DMR_for_fgwas_hyper_intergenic.bed"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)

write.table(
  hypo_intergenic[c(1:3,24)],
  paste0(outFolder,"DMR_for_fgwas_hypo_intergenic.bed"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)

rm(hyper,hypo,hyper_intergenic,hyper_not_prom,hypo_intergenic,hypo_not_prom)
### LMR, HMR ##################
LMR = list()
HMR = list()
stages2 = c("iPSC",stages,"islet")

for(c in stages2){
  LMR[[c]] = read.table(paste(outFolder,"LMR_HMR_0.5_0.5",
                              "/",c,"_",beta_lower_threshold,"_meanBeta_LMR.seg.bed",sep=""),skip = 1)
  LMR[[c]]$label = rep(c,nrow(LMR[[c]]))
  HMR[[c]] = read.table(paste(outFolder,"LMR_HMR_0.5_0.5",
                              "/",c,"_",beta_higher_threshold,"_meanBeta_HMR.seg.bed",sep=""),skip = 1)
  HMR[[c]]$label = rep(c,nrow(HMR[[c]]))
}

LMR = do.call("rbind", LMR)
HMR = do.call("rbind", HMR)

LMR = LMR[c(1:3,10)]
HMR = HMR[c(1:3,10)]

# sort
LMR = makeGRangesFromDataFrame(LMR, keep.extra.columns = T,seqnames.field = "V1",start.field = "V2", end.field = "V3")
HMR = makeGRangesFromDataFrame(HMR, keep.extra.columns = T,seqnames.field = "V1",start.field = "V2", end.field = "V3")

LMR = GenomicRanges::sort(LMR)
HMR = GenomicRanges::sort(HMR)

subset1 = as.data.frame(LMR)
subset2 = as.data.frame(HMR)

## plot abundance barplot
df = as.data.frame(cbind(table(subset2$label),table(subset1$label)))
df$stage = rownames(df)
colnames(df) = c("HMR","LMR","stage")
melted = melt(df,id.vars = "stage")

colnames(melted) = c("Stage","Type", "Number of regions")
melted$Stage = factor(melted$Stage, levels = c("iPSC",stages,"islet"), ordered = T)


p = ggplot(melted, aes( x = Stage, y = `Number of regions`)) +  
  geom_bar(stat = "identity", aes(fill = Type),position=position_dodge()) + 
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(filename =  paste0(outFolder,"LMRs_barplot.png"),
       plot = p,device = "png",width = 6, height = 5, units = "in", dpi = 400)

write.table(
  subset1[c(1:3,6)],
  paste0(outFolder,"LMR_for_fgwas.bed"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)
write.table(
  subset2[c(1:3,6)],
  paste0(outFolder,"HMR_for_fgwas.bed"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)
rm(LMR,subset1)

### LMR subsets ##################

features = c("CpGislands","CpGshores","CpGshelves","CpGopensea","promoter")
for(f in features){
  LMR = read.table(paste(outFolder,"LMR_feature_subsets/",beta_lower_threshold,"_LMR_",f,".bed",sep=""))
  
  # sort
  LMR = makeGRangesFromDataFrame(LMR, keep.extra.columns = T,seqnames.field = "V1",start.field = "V2", end.field = "V3")
  
  LMR = GenomicRanges::sort(LMR)
  
  LMR = as.data.frame(LMR)

  write.table(
    LMR[c(1:3,6)],
    paste0(outFolder,"LMR_feature_subsets/LMR_",f,"for_fgwas.bed"),
    sep = "\t",
    row.names = F,
    col.names = F,
    quote = F
  )

}


### annotate promoters from my whole dataset

inFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/"

# load matrix of normalized beta values
beta= read.csv(paste(inFolder,"quantile_normalised_beta_detP_0.01_nocrossreact.csv",sep=""), header=T,row.names = 1,check.names=F)


# get the  EPIC annotation data
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

beta= read.csv(paste(inFolder,"quantile_normalised_beta_detP_0.01_nocrossreact.csv",sep=""), header=T,row.names = 1,check.names=F)

beta_ann <- annEPIC[match(rownames(beta),annEPIC$Name),
                    c(1:4)]
beta_ann$end = beta_ann$pos + 1
beta_ann = makeGRangesFromDataFrame(beta_ann,start.field = "pos")

# Overall distribution

# CpGs genomic annotations

beta_ann <- annotatePeak(beta_ann, tssRegion=c(-2000, 500),
                         TxDb=txdb, annoDb="org.Hs.eg.db")


# change useless long intron exon downstream names
beta_ann@anno$annotation = gsub(".*Intron .*", "Intron", beta_ann@anno$annotation )
beta_ann@anno$annotation = gsub(".*Exon .*", "Exon", beta_ann@anno$annotation )
beta_ann@anno$annotation = gsub(".*Downstream .*", "Downstream", beta_ann@anno$annotation )
beta_ann@anno$annotation = gsub(".*Promoter .*", "Promoter", beta_ann@anno$annotation )

beta_ann = as.data.frame(beta_ann@anno)
beta_ann = beta_ann[c(1:3,6)]

write.table(
  beta_ann,
  paste0(outFolder,"allCpGs_genomic_annotation_for_fgwas.bed"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)
# None of the provided annotations were viable (all individual results had undefined confidence intervals)