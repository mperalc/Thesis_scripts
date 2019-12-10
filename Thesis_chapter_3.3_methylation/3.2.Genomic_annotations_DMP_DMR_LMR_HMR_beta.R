library(ChIPseeker)
library(ggplot2)
library(grid)
library(gridExtra)
library(GenomicRanges)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(minfiData)
library(minfi)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(reshape2)
library(regioneR)

currentDate <- Sys.Date() # to save date in name of output files

inFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/"
outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/"

# annotation
annEPIC <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)


# load matrix of normalized beta values
beta= read.csv(paste(inFolder,"quantile_normalised_beta_detP_0.01_nocrossreact.csv",sep=""), header=T,row.names = 1,check.names=F)

beta_ann <- annEPIC[match(rownames(beta),annEPIC$Name),
                      c(1:4)]
beta_ann$end = beta_ann$pos + 1
beta_ann = makeGRangesFromDataFrame(beta_ann,start.field = "pos")

# Overall distribution

# CpGs genomic annotations

beta_ann <- annotatePeak(beta_ann, tssRegion=c(-2000, 500),
                                  TxDb=txdb, annoDb="org.Hs.eg.db")

# png(
#   paste(outFolder, "genomic_annotations_ChIPseeker", currentDate, ".png", sep = ""),
#   width = 7,
#   height = 8,
#   units = "in",
#   res = 400,
#   type = "cairo"
# )
# plotAnnoPie(beta_ann)
# 
# dev.off()

# Distribution of beta values per island feature

beta_ann <- annEPIC[match(rownames(beta),annEPIC$Name),
                    c(4,19)]
beta_ann = merge(beta, beta_ann,by = 0)
beta_ann = beta_ann[c(2:33,35)]


# melt
melted = melt(as.data.frame(beta_ann),id.vars = "Relation_to_Island")

melted$Relation_to_Island = gsub(".*Shore.*", "Shore", melted$Relation_to_Island )
melted$Relation_to_Island = gsub(".*Shelf.*", "Shelf", melted$Relation_to_Island )

p<-ggplot(melted, aes(x=Relation_to_Island, y=value)) +
  geom_violin(trim=FALSE) + 
  ggtitle("Beta values per CpG island feature") +
  xlab("CpG annotations") +
  ylab("Beta value") +
  theme_bw()

ggsave(filename = paste(outFolder, "CpGisland_annotations_beta_distribution_overall", currentDate, ".png", sep = ""),plot = p,
       width = 5,height = 5,units = "in",dpi = 400)

### Distribution of beta values per genomic feature
beta_ann <- annEPIC[match(rownames(beta),annEPIC$Name),
                    c(1:4)]
beta_ann$end = beta_ann$pos + 1
beta_ann = makeGRangesFromDataFrame(beta_ann,start.field = "pos")
beta_ann <- annotatePeak(beta_ann, tssRegion=c(-2000, 500),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

# change useless long intron exon downstream names
beta_ann@anno$annotation = gsub(".*Intron .*", "Intron", beta_ann@anno$annotation )
beta_ann@anno$annotation = gsub(".*Exon .*", "Exon", beta_ann@anno$annotation )
beta_ann@anno$annotation = gsub(".*Downstream .*", "Downstream", beta_ann@anno$annotation )
beta_ann@anno$annotation = gsub(".*Promoter .*", "Promoter", beta_ann@anno$annotation )

beta_ann = as.data.frame(beta_ann@anno)
beta_ann = merge(beta, beta_ann[c(1,2,6)],by = 0)
beta_ann = beta_ann[c(2:33,36)]
melted = melt(beta_ann,id.vars = "annotation")
melted$stage = rep( c("iPSC","iPSC","iPSC","DE","DE","DE","GT","GT","GT","PF","PF","PF","PE","PE","PE","EP","EN","EN","EN","BLC","BLC","islet","islet","islet","islet","islet","islet","islet","islet","islet","islet","islet"),each=nrow(beta_ann))


png(
  paste(outFolder, "genomic_annotations_beta_distribution_overall", currentDate, ".png", sep = ""),
  width = 6.5,
  height = 5,
  units = "in",
  res = 400,
  type = "cairo"
)
p<-ggplot(melted, aes(x=annotation, y=value)) +
  geom_violin(trim=FALSE) + 
  ggtitle("Beta values per genomic annotation") +
  xlab("Genomic annotations") +
  ylab("Beta value") +
  theme_bw()
plot(p)

dev.off()

# Distribution per feature per stage

# 
# png(
#   paste(outFolder, "genomic_annotations_beta_distribution_per_stage", currentDate, ".png", sep = ""),
#   width = 8,
#   height = 8,
#   units = "in",
#   res = 400,
#   type = "cairo"
# )
# p<-ggplot(melted, aes(x=annotation, y=value, fill=stage)) +
#   geom_violin() + 
#   ggtitle("Beta values per genomic annotation") +
#   xlab("Genomic annotations") +
#   ylab("Beta value") +
#   theme_bw()
# plot(p)
# 
# dev.off()


# Genomic annotations of DMP

# Genomic annotations of DMR
# read in data
outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMR/"


results_df_list_hyper  = read.csv(file=paste(outFolder,"DMRs_hypermethylated_annotated_allstages.csv",sep=""))
results_df_list_hypo =  read.csv(file=paste(outFolder,"DMRs_hypomethylated_annotated_allstages.csv",sep=""))
test = rbind(results_df_list_hyper,results_df_list_hypo)
summarytools::descr(test$width) # getting size of DMRs
summarytools::descr(test$meanbetafc) ## changes in beta values
summarytools::descr(results_df_list_hyper$meanbetafc) ## changes in beta values
summarytools::descr(results_df_list_hypo$meanbetafc) ## changes in beta values


results_df_list_hyper = makeGRangesFromDataFrame(results_df_list_hyper)
results_df_list_hypo = makeGRangesFromDataFrame(results_df_list_hypo)

## significance of overlap with promoter annotations
PR <- promoters(txdb, upstream=2000, downstream=500)
PR = PR[seqnames(PR) %in% paste0("chr",c(1:22)),]

promoters_hyper <- overlapPermTest(results_df_list_hyper, PR, ntimes=100, genome="hg19", count.once=TRUE, verbose = T)
promoters_hypo<- overlapPermTest(results_df_list_hypo, PR, ntimes=100, genome="hg19", count.once=TRUE, verbose = T)
promoters_hyper
promoters_hypo

##### LMR, HMR
outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/LMR_HMR_0.5_0.5/"
LMR = read.table(paste0(outFolder,"LMR_for_fgwas.bed"),sep = "\t",header = F)
HMR = read.table(paste0(outFolder,"HMR_for_fgwas.bed"),sep = "\t",header  = F)
LMR = unique(LMR[c(1:3)])
HMR = unique(HMR[c(1:3)])

LMR$size = LMR$V3 - LMR$V2
HMR$size = HMR$V3 - HMR$V2
summarytools::descr(LMR$size)
summarytools::descr(HMR$size)

colnames(LMR) = c("chr","start","end","size")

LMR.gr = makeGRangesFromDataFrame(LMR)
LMR_ann <- annotatePeak(LMR.gr, tssRegion=c(-2000, 500),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

png(
  paste(outFolder, "genomic_annotations_LMR", currentDate, ".png", sep = ""),
  width = 7,
  height = 8,
  units = "in",
  res = 400,
  type = "cairo"
)
plotAnnoPie(LMR_ann)

dev.off()
## pathway analysis
outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/LMR_HMR_0.5_0.5/"
LMR = read.table(paste0(outFolder,"LMR_for_fgwas.bed"),sep = "\t",header = F)
# HMR does not make sense because regions are too big
colnames(LMR) = c("chr","start","end","stage")

stages <- c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

## all per stage
LMR_list = list()
for(s in stages){
  LMR_list [[s]] = LMR[LMR$stage == s,]
  LMR_list [[s]] = GRanges(LMR_list[[s]])
}

per_stage_GRL = GRangesList(LMR_list)
# Get list of Entrez genes to compare all stages at same time
peakAnnoList <- lapply(per_stage_GRL, annotatePeak, TxDb=txdb,
                       tssRegion=c(-2000, 500), verbose=T)
genes_specific = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)




all_specific_pathway = compareCluster(geneCluster   = genes_specific,
                                      fun           = "enrichPathway",
                                      pvalueCutoff  = 0.05,
                                      pAdjustMethod = "BH",
                                      readable = TRUE)

png(
  paste0(outFolder, "enrichment_reactome_top20_LMR.png"),units = "in",res = 400,type = "cairo",
  width = 9,
  height = 5
)
dotplot(all_specific_pathway,showCategory = 20, title = "Reactome Pathway Enrichment Analysis") 
dev.off()
