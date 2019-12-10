# Segment methylation to get highly and lowly methylated regions

library(methylKit)
library(GenomicRanges)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(annotatr)
library(ggplot2)


### initial parameters
inFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/"
outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/"
beta_lower_threshold = 0.50
beta_higher_threshold = 0.50

###
# beta values
beta = read.csv(paste(inFolder,"quantile_normalised_beta_detP_0.01_nocrossreact.csv",sep=""), header=T,row.names = 1,check.names=F)

# get the latest EPIC annotation data
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

CpGs_anno = merge(x = beta, y = annEPIC, by = 0)
rownames(CpGs_anno) = CpGs_anno$Row.names

##### Get mean beta values at each position by stage ############
CpGs_anno_mean = list()

# iPSC
CpGs_anno_mean[["iPSC"]] = CpGs_anno[c(2:4,34,35)]
CpGs_anno_mean[["iPSC"]]$beta = rowMeans(as.matrix(CpGs_anno_mean[["iPSC"]][c(1:3)]))
CpGs_anno_mean[["iPSC"]] = makeGRangesFromDataFrame(CpGs_anno_mean[["iPSC"]][c("chr","pos","beta")],
                                                        start.field="pos",end.field = "pos",keep.extra.columns = TRUE)

# DE
CpGs_anno_mean[["DE"]] = CpGs_anno[c(5:7,34,35)]
CpGs_anno_mean[["DE"]]$beta = rowMeans(as.matrix(CpGs_anno_mean[["DE"]][c(1:3)]))
CpGs_anno_mean[["DE"]] = makeGRangesFromDataFrame(CpGs_anno_mean[["DE"]][c("chr","pos","beta")],
                                                    start.field="pos",end.field = "pos",keep.extra.columns = TRUE)

# GT
CpGs_anno_mean[["GT"]] = CpGs_anno[c(8:10,34,35)]
CpGs_anno_mean[["GT"]]$beta = rowMeans(as.matrix(CpGs_anno_mean[["GT"]][c(1:3)]))
CpGs_anno_mean[["GT"]] = makeGRangesFromDataFrame(CpGs_anno_mean[["GT"]][c("chr","pos","beta")],
                                                  start.field="pos",end.field = "pos",keep.extra.columns = TRUE)

# PF
CpGs_anno_mean[["PF"]] = CpGs_anno[c(11:13,34,35)]
CpGs_anno_mean[["PF"]]$beta = rowMeans(as.matrix(CpGs_anno_mean[["PF"]][c(1:3)]))
CpGs_anno_mean[["PF"]] = makeGRangesFromDataFrame(CpGs_anno_mean[["PF"]][c("chr","pos","beta")],
                                                  start.field="pos",end.field = "pos",keep.extra.columns = TRUE)

# PE
CpGs_anno_mean[["PE"]] = CpGs_anno[c(14:16,34,35)]
CpGs_anno_mean[["PE"]]$beta = rowMeans(as.matrix(CpGs_anno_mean[["PE"]][c(1:3)]))
CpGs_anno_mean[["PE"]] = makeGRangesFromDataFrame(CpGs_anno_mean[["PE"]][c("chr","pos","beta")],
                                                  start.field="pos",end.field = "pos",keep.extra.columns = TRUE)
# EP
CpGs_anno_mean[["EP"]] = CpGs_anno[c(17,34,35)]
CpGs_anno_mean[["EP"]]$beta = as.matrix(CpGs_anno_mean[["EP"]][c(1)])
CpGs_anno_mean[["EP"]] = makeGRangesFromDataFrame(CpGs_anno_mean[["EP"]][c("chr","pos","beta")],
                                                  start.field="pos",end.field = "pos",keep.extra.columns = TRUE)

# EN
CpGs_anno_mean[["EN"]] = CpGs_anno[c(18:20,34,35)]
CpGs_anno_mean[["EN"]]$beta = rowMeans(as.matrix(CpGs_anno_mean[["EN"]][c(1:3)]))
CpGs_anno_mean[["EN"]] = makeGRangesFromDataFrame(CpGs_anno_mean[["EN"]][c("chr","pos","beta")],
                                                  start.field="pos",end.field = "pos",keep.extra.columns = TRUE)

# BLC
CpGs_anno_mean[["BLC"]] = CpGs_anno[c(21:22,34,35)]
CpGs_anno_mean[["BLC"]]$beta = rowMeans(as.matrix(CpGs_anno_mean[["BLC"]][c(1:2)]))
CpGs_anno_mean[["BLC"]] = makeGRangesFromDataFrame(CpGs_anno_mean[["BLC"]][c("chr","pos","beta")],
                                                  start.field="pos",end.field = "pos",keep.extra.columns = TRUE)
# islet
CpGs_anno_mean[["islet"]] = CpGs_anno[c(23:35)]
CpGs_anno_mean[["islet"]]$beta = rowMeans(as.matrix(CpGs_anno_mean[["islet"]][c(1:11)]))
CpGs_anno_mean[["islet"]] = makeGRangesFromDataFrame(CpGs_anno_mean[["islet"]][c("chr","pos","beta")],
                                                   start.field="pos",end.field = "pos",keep.extra.columns = TRUE)

# segmentation analysis will usually reveal high or low methylated regions, where low methylated regions could be interesting for 
# gene regulation. The algorithm first finds segments that have CpGs with similar methylation levels, then those segments are 
# classified to segment groups based on their mean methylation levels. 
# This enables us to group segments with similar methylation levels to the same class.

# Matthias defines LMR as median methylation < 0.50 with < 30 CpGs from the pipeline methylseek
# UMR are those with more than 30
# Fernandez-Jimenez et al 2017 defines it as those not in CpG islands, that correspond to enhancers
# While UMR (inCpG islands) would correspond to promoters
# subset these and run fgwas


# CpG island annotations
annots = c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic')
# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg19', annotations = annots)

LMR_annotated= list()
LMR_annotated_all= list()

for (n in names(CpGs_anno_mean)) {
res=methSeg(CpGs_anno_mean[[n]],diagnostic.plot=TRUE,maxInt=100,minSeg=2,G=1:4)
# forcing only 4 groups does not change segments or means. It works only for classification purposes, as it learns
# the groups from the distibution of mean beta values per position per stage

# Intersect the regions we read in with the annotations
LMR_annotated[[n]] = annotate_regions(
  regions = res[which(res$seg.mean<=beta_lower_threshold),],
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)
LMR_annotated[[n]] = as.data.frame(LMR_annotated[[n]])
LMR_annotated[[n]]$stage = rep(n,nrow(LMR_annotated[[n]]))
# minimum segment length = 2 (as MAtthias)

# All to save annotations
LMR_annotated_all[[n]] = annotate_regions(
  regions = res,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)
LMR_annotated_all[[n]] = as.data.frame(LMR_annotated_all[[n]])
LMR_annotated_all[[n]]$stage = rep(n,nrow(LMR_annotated_all[[n]]))
# 
# methSeg2bed(res[which(res$seg.mean<=beta_lower_threshold),],filename=paste0(outFolder,n,"_",beta_lower_threshold,"_meanBeta_LMR.seg.bed"))
# methSeg2bed(res[which(res$seg.mean>=beta_higher_threshold),],filename=paste0(outFolder,n,"_",beta_higher_threshold,"_meanBeta_HMR.seg.bed"))

}

LMR_annotated = do.call("rbind",LMR_annotated)
# save all segments for plotting heatman of mean beta
LMR_annotated_all = do.call("rbind",LMR_annotated_all)
LMR_annotated_all = LMR_annotated_all[c(1:3,7,8,11,20:22)]
LMR_annotated_all = unique(LMR_annotated_all)

write.table(LMR_annotated_all,paste0(outFolder,"regions_annotated_all.bed"),
            row.names = F,quote = F,col.names = F)



annots_order = c('hg19_genes_promoters','hg19_genes_intergenic',
                 "hg19_cpg_islands","hg19_cpg_shores", "hg19_cpg_shelves","hg19_cpg_inter" )
LMR_plot_coannotations = plot_coannotations(
  annotated_regions = LMR_annotated,
  axes_label = 'Annotations',
  annotation_order = annots_order,
  plot_title = 'LMR: Regions in Pairs of Annotations')
print(LMR_plot_coannotations)
ggsave(paste0(outFolder,beta_lower_threshold,"_","LMR_coannotations.png"),plot = LMR_plot_coannotations,device = "png", 
       width = 6, height = 6, dpi = 400, units = "in")

stage = c("iPSC","DE","GT","PF","PE","EP","EN","BLC","islet")

for (s in stage){
LMR_plot_num_annotated= plot_numerical_coannotations(
  annotated_regions =LMR_annotated[LMR_annotated$stage == s,],
  x = 'seg.mean',
  annot1 = 'hg19_cpg_islands',
  annot2 = 'hg19_genes_promoters',
  bin_width = 0.01,
  plot_title = 'LMR in CpG islands and promoters',
  x_label = 'Mean Methylation')

ggsave(paste0(outFolder,s,"_",beta_lower_threshold,"_LMR_coannotations_cpg_promoters.png"),plot = LMR_plot_num_annotated,device = "png", 
       width = 6, height = 6, dpi = 400, units = "in")

LMR_plot_num_annotated= plot_numerical_coannotations(
  annotated_regions =LMR_annotated[LMR_annotated$stage == s,],
  x = 'seg.mean',
  annot1 = 'hg19_cpg_inter',
  annot2 = 'hg19_genes_intergenic',
  bin_width = 0.01,
  plot_title = 'LMR in CpG open sea and intergenic regions',
  x_label = 'Mean Methylation')
ggsave(paste0(outFolder,s,"_",beta_lower_threshold,"_LMR_coannotations_opensea_intergenic.png"),plot = LMR_plot_num_annotated,device = "png", 
       width = 6, height = 6, dpi = 400, units = "in")

}
# subset by genomic feaure fpr all stages
LMR_CpGislands = LMR_annotated[ LMR_annotated$annot.type =='hg19_cpg_islands',c(1:3,22)]
LMR_CpGislands = unique(LMR_CpGislands)

LMR_CpGshores= LMR_annotated[ LMR_annotated$annot.type =='hg19_cpg_shores',c(1:3,22)]
LMR_CpGshores = unique(LMR_CpGshores)

LMR_CpGshelves= LMR_annotated[ LMR_annotated$annot.type =='hg19_cpg_shelves',c(1:3,22)]
LMR_CpGshelves = unique(LMR_CpGshelves)

LMR_CpGopensea= LMR_annotated[ LMR_annotated$annot.type =='hg19_cpg_inter',c(1:3,22)]
LMR_CpGopensea = unique(LMR_CpGopensea)

LMR_promoter= LMR_annotated[ LMR_annotated$annot.type =='hg19_genes_promoters',c(1:3,22)]
LMR_promoter = unique(LMR_promoter)

LMR_not_promoter = findOverlaps(makeGRangesFromDataFrame(LMR_annotated,keep.extra.columns = T),makeGRangesFromDataFrame(LMR_promoter,keep.extra.columns = T))
exclude = unique(queryHits(LMR_not_promoter))
exclude = setdiff(rownames(LMR_annotated),rownames(LMR_annotated[exclude,]))
LMR_not_promoter = LMR_annotated[exclude,]
LMR_not_promoter= LMR_not_promoter[ ,c(1:3,22)]
LMR_not_promoter = unique(LMR_not_promoter)
# and CpGopensea NOT in promoters (exclude those)

LMR_CpGopensea_not_prom = findOverlaps(makeGRangesFromDataFrame(LMR_CpGopensea,keep.extra.columns = T),makeGRangesFromDataFrame(LMR_promoter,keep.extra.columns = T))
exclude = unique(queryHits(LMR_CpGopensea_not_prom))
exclude = setdiff(rownames(LMR_CpGopensea),rownames(LMR_CpGopensea[exclude,]))
LMR_CpGopensea_not_prom = LMR_CpGopensea[exclude,]

write.table(LMR_CpGislands,paste0(outFolder,beta_lower_threshold,"_","LMR_CpGislands.bed"),row.names = F,quote = F,col.names = F)
write.table(LMR_CpGshores,paste0(outFolder,beta_lower_threshold,"_","LMR_CpGshores.bed"),row.names = F,quote = F,col.names = F)
write.table(LMR_CpGshelves,paste0(outFolder,beta_lower_threshold,"_","LMR_CpGshelves.bed"),row.names = F,quote = F,col.names = F)
write.table(LMR_CpGopensea,paste0(outFolder,beta_lower_threshold,"_","LMR_CpGopensea.bed"),row.names = F,quote = F,col.names = F)
write.table(LMR_promoter,paste0(outFolder,beta_lower_threshold,"_","LMR_promoter.bed"),row.names = F,quote = F,col.names = F)
write.table(LMR_not_promoter,paste0(outFolder,beta_lower_threshold,"_","LMR_not_promoter.bed"),row.names = F,quote = F,col.names = F)

##
write.table(LMR_CpGopensea_not_prom,paste0(outFolder,beta_lower_threshold,"_","LMR_CPGopensea_NOT_promoter.bed"),row.names = F,quote = F,col.names = F)

## fgwas error: None of the provided annotations were viable (all individual results had undefined confidence intervals). Terminating analysis.
## more lenient: in GpGopensea and in intergenic regions
# might have some overlap with promoter in opposite strand

LMR_intergenic= LMR_annotated[ LMR_annotated$annot.type =='hg19_genes_intergenic',c(1:3,22)]
LMR_intergenic = unique(LMR_intergenic)


LMR_CpGopensea_intergenic = findOverlaps(makeGRangesFromDataFrame(LMR_CpGopensea,keep.extra.columns = T),makeGRangesFromDataFrame(LMR_intergenic,keep.extra.columns = T))
include = unique(queryHits(LMR_CpGopensea_intergenic))
LMR_CpGopensea_intergenic = LMR_CpGopensea[include,]

write.table(LMR_CpGopensea_intergenic,paste0(outFolder,beta_lower_threshold,"_","LMR_CpGopensea_intergenic.bed"),row.names = F,quote = F,col.names = F)


