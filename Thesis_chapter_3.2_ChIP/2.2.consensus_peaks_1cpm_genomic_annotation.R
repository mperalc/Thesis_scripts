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

basedir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks/trimmed/"

consensus = read_tsv("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks/merged/H3K27ac_binary_optimal.bed",col_names = T)


# 1 CPM-filtered set of peaks

cpm = read.table(paste0(basedir,"CPMs_H3K27ac_1CPM_trim_optimal.txt"))


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
              title="Distribution of H3K27ac ChIP-seq peaks relative to TSS")

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

png(
  paste(basedir, '/GenomicAnnotationBar_1CPM_filt_consensus_set.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 6,
  height = 6
)

plotAnnoBar(peakAnnoList, title = "Distribution of H3K27ac ChIP-seq peaks relative to genomic features")

dev.off()

png(
  paste(basedir, '/TSSDistancenBar_1CPM_filt_consensus_set.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 6,
  height = 6
)

plotDistToTSS(peakAnnoList,
              title="Distribution of  H3K27ac ChIP-seq  peaks relative to TSS")
dev.off()


## stage-specific peaks

stage_specific = list()
for(s in stages){
  stage_specific[[s]] = cpm_filt_consensus[cpm_filt_consensus[s]==1 & rowSums(cpm_filt_consensus[5:12]) ==1,c("Chr","Start","End")]
  stage_specific[[s]]$stage =  rep(s,nrow(stage_specific[[s]]))
  stage_specific[[s]] = GRanges(stage_specific[[s]])
  
  
}

stage_specific_GRL = GRangesList(stage_specific)

stage_peakAnnoList <- lapply(stage_specific_GRL, annotatePeak, TxDb=txdb,
                       tssRegion=c(-2000, 500), verbose=FALSE)


png(
  paste(basedir, '/GenomicAnnotationBar_1CPM_filt_stage_specific.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 6,
  height = 6
)

plotAnnoBar(stage_peakAnnoList, title = "H3K27ac ChIP-seq peaks relative to genomic features: stage-specific")

dev.off()

png(
  paste(basedir, '/TSSDistancenBar_1CPM_filt_stage_specific.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 6,
  height = 6
)

plotDistToTSS(stage_peakAnnoList,
              title="H3K27ac ChIP-seq  peaks relative to TSS: stage-specific")
dev.off()




### calculate significance of genomic annotations 
library(regioneR)
PR <- promoters(txdb, upstream=2000, downstream=500)
PR = PR[seqnames(PR) %in% paste0("chr",c(1:22)),]

s="PF"
promoters <- overlapPermTest(per_stage_GRL[[s]], PR, ntimes=100, genome="hg19", count.once=TRUE, verbose = T)
promoters
plot(promoters, plotType = "Tailed")

png(paste0(basedir,"permtest.png"))
plot(promoters, plotType = "Tailed")
dev.off()

# Functional enrichment


# Get list of Entrez genes to compare all stages at same time
genes_all = lapply(per_stage_GRL, seq2gene,TxDb=txdb,
                        tssRegion=c(-2000, 500),flankDistance = 100000)

# GO analysis
all_GO = compareCluster(genes_all, fun="enrichGO", OrgDb='org.Hs.eg.db',pvalueCutoff = 0.01)
dotplot(all_GO, showCategory=10)

# Pathway analysis
all_pathway = compareCluster(genes_all, fun="enrichPathway")

dotplot(all_pathway,
        showCategory=30) # They all look very similar

# Disease analysis
all_disease <- compareCluster(genes_all,fun = "enrichDO", pvalueCutoff  = 0.01,minGSSize     = 5)
dotplot(all_disease, showCategory=10)

# What about stage-specific?
# 
genes_specific = lapply(
  stage_specific_GRL,
  seq2gene,
  TxDb = txdb,
  tssRegion = c(-2000, 500),
  flankDistance = 100000
)
specific_GO = compareCluster(
  genes_specific,
  fun = "enrichGO",
  OrgDb = 'org.Hs.eg.db',
  pvalueCutoff = 0.01,
  minGSSize     = 5
)
dotplot(specific_GO, showCategory = 30)

specific_pathway = compareCluster(
  genes_specific,
  fun = "enrichPathway",
  pvalueCutoff = 0.01,
  minGSSize  = 5
)

dotplot(specific_pathway, 
        showCategory = 50)

# Disease analysis
specific_disease <-
  compareCluster(
    genes_specific,
    fun = "enrichDO",
    pvalueCutoff  = 0.01,
    minGSSize     = 5
  )
dotplot(specific_disease, showCategory = 30)

# Nothing worth looking at here I'm afraid 

## without seq2 gene (that performs many to many peak and gene mapping)

genes = lapply(stage_peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compEnrich <- compareCluster(geneCluster   = genes,
                           fun           = "enrichPathway",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")


png(
  paste(basedir, '/enrichment_reactome_top20_1CPM_filt_stage_specific.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 12,
  height = 12
)

dotplot(compEnrich, showCategory = 20, title = "Reactome Pathway Enrichment Analysis")

dev.off()
