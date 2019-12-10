# Genomic annotations of 1CPM filtered consensus peaks
# Also upsetplot of shared and unique peaks
# Also functional enrichment (GO, pathways, etc)

library(readr)
library(summarytools)
library(ggplot2)
library(grid)
library(gridExtra)
library(GenomicRanges)
library(ChIPseeker)
library(UpSetR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(ReactomePA)
library(clusterProfiler)
library(DOSE)

stages <- c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

# After merging among stages with bedtools and gathering counts with featureCounts, I get the combined (consensus) set of peaks

basedir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/IDR_peaks/merged/"

consensus = read_tsv(paste0(basedir,"ATAC_binary_conservative_normalQuality_narrowpeaks.bed"),col_names = T)


# 1 CPM-filtered set of peaks

cpm = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/counts/CPMs_atac-seq_1CPM_trim_conservativeCounts_normalquality.txt")
# total consensus peaks per stage
colSums(consensus[c(5:12)])
mean(colSums(consensus[c(5:12)]))

# Unique per stage
png(
  paste(basedir, '/upset_consensus_set_top15.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 12,
  height = 8
)

upset(as.data.frame(consensus[,c(5:12)]), order.by = "freq",nsets = 8,point.size = 2, line.size = 1, 
      mainbar.y.label = "Number of peaks", sets.x.label = "Peaks per stage", 
      text.scale = c(2, 2, 2, 2, 2, 2),
      sets = rev(stages), keep.order = TRUE,
      nintersects = 15)
# Only plotting top 15 intersects

dev.off()

# Unique peaks per stage

14098/54585 # iPSC
12209/53001 # DE
3965/49759 # GT
7247/64162 # PF
5536/57151 # PE
885/30390 # EP
8269/55757 # EN
1966/38826 # BLC

# mean of unique peaks per stage
mean(c(14098/54585,12209/53001,3965/49759,7247/64162,5536/57151,885/30390,8269/55757,1966/38826 ))


##### Now get distribution but for 1CPM filtered set
# Overall distribution
# Peak genomic annotations
consensus = as.data.frame(consensus)
rownames(consensus) = consensus$Name
cpm_filt_consensus = consensus[rownames(cpm),]

peakAnnoConsensus <- annotatePeak(GRanges(as.data.frame(cpm_filt_consensus)), tssRegion=c(-2000, 500),
                                  TxDb=txdb, annoDb="org.Hs.eg.db")

plotAnnoPie(peakAnnoConsensus)

plotDistToTSS(peakAnnoConsensus,
              title="Distribution of ATAC-seq peaks relative to TSS")

# dist 0-1 kb
sum(abs(peakAnnoConsensus@anno@elementMetadata$distanceToTSS)>=0 & abs(peakAnnoConsensus@anno@elementMetadata$distanceToTSS)<1000) / nrow(cpm_filt_consensus)

# dist 1-3 kb
sum(abs(peakAnnoConsensus@anno@elementMetadata$distanceToTSS)>=1000 & abs(peakAnnoConsensus@anno@elementMetadata$distanceToTSS)<3000) / nrow(cpm_filt_consensus)

# dist 3-5 kb
sum(abs(peakAnnoConsensus@anno@elementMetadata$distanceToTSS)>=3000 & abs(peakAnnoConsensus@anno@elementMetadata$distanceToTSS)<5000) / nrow(cpm_filt_consensus)

# dist 5-10 kb
sum(abs(peakAnnoConsensus@anno@elementMetadata$distanceToTSS)>=5000 & abs(peakAnnoConsensus@anno@elementMetadata$distanceToTSS)<10000) / nrow(cpm_filt_consensus)

# dist 10-100kb
sum(abs(peakAnnoConsensus@anno@elementMetadata$distanceToTSS)>=10000 & abs(peakAnnoConsensus@anno@elementMetadata$distanceToTSS)<100000) / nrow(cpm_filt_consensus)

# dist 10-100kb
sum(abs(peakAnnoConsensus@anno@elementMetadata$distanceToTSS)>=100000)/ nrow(cpm_filt_consensus)



# Distribution per stage

# Grab peaks per stage

per_stage = list()

cpm_filt_consensus = as.data.frame(cpm_filt_consensus)

for (s in stages) {
  subset_c = cpm_filt_consensus[,c("Name","Chr","Start","End",s)]
  per_stage[[s]] = subset_c[subset_c[s]==1,]
  per_stage[[s]] = GRanges(per_stage[[s]])
}
per_stage_GRL = GRangesList(per_stage)

peakAnnoList <- lapply(per_stage_GRL, annotatePeak, TxDb=txdb,
                       tssRegion=c(-2000, 500), verbose=FALSE)
stats = list()
for(s in stages){
  stats[[s]] = peakAnnoList[[s]]@annoStat
  stats[[s]]$stage = rep(s,nrow(stats[[s]]))
}
stats = do.call("rbind",stats)
descr(stats[stats$Feature == "Promoter (<=1kb)","Frequency"])
descr(stats[stats$Feature == "Promoter (1-2kb)","Frequency"])
descr(stats[stats$Feature == "Distal Intergenic","Frequency"])
descr(stats[stats$Feature == "1st Intron","Frequency"])
descr(stats[stats$Feature == "Other Intron","Frequency"])



png(
  paste(basedir, '/GenomicAnnotationBar_cpm_filt_consensus_set.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 6,
  height = 6
)

plotAnnoBar(peakAnnoList, title = "Distribution of ATAC-seq peaks relative to genomic features")

dev.off()

png(
  paste(basedir, '/TSSDistancenBar_cpm_filt_consensus_set.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 6,
  height = 6
)

p = plotDistToTSS(peakAnnoList,
                  title="Distribution of ATAC-seq peaks relative to TSS")
plot(p)
dev.off()

p=as.data.frame(p$data)
descr(p[p$Feature=="0-1kb" & p$sign==1,"freq"])
descr(p[p$Feature=="1-3kb" & p$sign==1,"freq"])
descr(p[p$Feature=="5-10kb" & p$sign==1,"freq"])
descr(p[p$Feature=="10-100kb" & p$sign==1,"freq"])
descr(p[p$Feature==">100kb" & p$sign==1,"freq"])


# Functional enrichment

# selecting stage-specific peaks
stage_specific = list()
for(s in stages){
  stage_specific[[s]] = cpm_filt_consensus[which(cpm_filt_consensus[s] == 1 & rowSums(cpm_filt_consensus[c(setdiff(stages,s))])==0),]
  stage_specific[[s]] = GRanges(stage_specific[[s]])
  
}
stage_specific_GRL = GRangesList(stage_specific)


BLCpeak_specific = seq2gene(stage_specific_GRL$BLC, tssRegion = c(-2000, 500), flankDistance = 3000, TxDb=txdb,)
BLCpeak_specific_pathway = enrichPathway(BLCpeak_specific)
dotplot(BLCpeak_specific_pathway)

emapplot(BLCpeak_specific_pathway)
cnetplot(BLCpeak_specific_pathway) # change to gene name list before calling enrichPathway for this to look nicer

# changing entrez to gene id for enrichment plot
BLCpeak_specific_pathway_readable <- setReadable(BLCpeak_specific_pathway, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(BLCpeak_specific_pathway_readable,categorySize="pvalue")
heatplot(BLCpeak_specific_pathway_readable)

# Get list of Entrez genes to compare all stages at same time
genes_specific = lapply(stage_specific_GRL, seq2gene,TxDb=txdb,
                        tssRegion=c(-2000, 500),flankDistance = 3000)

# GO analysis
all_specific_GO = compareCluster(genes_specific, fun="enrichGO", OrgDb='org.Hs.eg.db',pvalueCutoff = 0.01)
dotplot(all_specific_GO, showCategory=50)

# Pathway analysis
all_specific_pathway = compareCluster(genes_specific, fun="enrichPathway")

dotplot(all_specific_pathway) # They all look very similar

# Disease analysis
all_specific_disease <- compareCluster(genes_specific,fun = "enrichDO", pvalueCutoff  = 0.01,minGSSize     = 5)
dotplot(all_specific_disease, showCategory=50)

# What about all genes?
# 
genes_unspecific = lapply(
  per_stage_GRL,
  seq2gene,
  TxDb = txdb,
  tssRegion = c(-2000, 500),
  flankDistance = 3000
)
all_unspecific_GO = compareCluster(
  genes_unspecific,
  fun = "enrichGO",
  OrgDb = 'org.Hs.eg.db',
  pvalueCutoff = 0.01,
  minGSSize     = 5
)
dotplot(all_unspecific_GO, showCategory = 50)

all_unspecific_pathway = compareCluster(
  genes_unspecific,
  fun = "enrichPathway",
  pvalueCutoff = 0.01,
  minGSSize  = 5
)
dotplot(all_unspecific_pathway, showCategory = 50)


# Disease analysis
all_unspecific_disease <-
  compareCluster(
    genes_unspecific,
    fun = "enrichDO",
    pvalueCutoff  = 0.01,
    minGSSize     = 5
  )
dotplot(all_unspecific_disease, showCategory = 50)

# Nothing worth looking at here I'm afraid 


## without seq2 gene (that performs many to many peak and gene mapping)

genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compEnrich <- compareCluster(geneCluster   = genes,
                             fun           = "enrichPathway",
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH")


png(
  paste(basedir, '/enrichment_reactome_top20_1CPM_filt_all_peaks.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 12,
  height = 12
)

dotplot(compEnrich, showCategory = 20, title = "Reactome Pathway Enrichment Analysis")

dev.off()

## stage-specific peaks
peakAnnoList_specifc =  lapply(stage_specific_GRL, annotatePeak, TxDb=txdb,
                             tssRegion=c(-2000, 500), verbose=FALSE)
genes = lapply(peakAnnoList_specifc, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compEnrich <- compareCluster(geneCluster   = genes,
                             fun           = "enrichPathway",
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH")


png(
  paste(basedir, '/enrichment_reactome_top20_1CPM_filt_stage_specific_peaks.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 12,
  height = 12
)

dotplot(compEnrich, showCategory = 20, title = "Reactome Pathway Enrichment Analysis")

dev.off()

library(annotate)
## annotate SOME entrez ids to gene names 
# EN insulin processing
getSYMBOL(c("56605","3800","25924","55763","1363","64924","5122","55770","60412","169026","3799","3798","10640","5873","5126"), 'org.Hs.eg')
##EN regulation of beta cell development
getSYMBOL(c("2494","4760","8850","3280","3516","4825","55534","2255","222546","5078","3174","50674","5080","2308","3175","6928","208","3642","3170"),'org.Hs.eg')
