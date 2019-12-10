# Overlap H3K27ac and ATAC-seq and get
# Enhancers : overlaps not at -2000, +500 bp from TSS
# Promoters: Those within those boundaries


library(GenomicRanges)
library(UpSetR)
library(data.table)
library(annotatr)
library(annotate)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(ReactomePA)
library(clusterProfiler)
library(DOSE)
library('org.Hs.eg.db')

#ATAC
atac = read.table(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/IDR_peaks/merged/ATAC_binary_conservative_normalQuality_narrowpeaks.bed",
  header = T
)

chrs = paste0("chr", rep(1:22))

# filter so that it only contains autosomes
atac = atac[atac$Chr %in% chrs, ]
peaks.gr = makeGRangesFromDataFrame(atac, keep.extra.columns = T)

#ChIP
datChip = read.table(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks/merged/H3K27ac_binary_optimal.bed",
  header = T
)
# reorder
datChip = datChip[, c("Name",
                      "Chr",
                      "Start",
                      "End",
                      "iPSC",
                      "DE",
                      "GT",
                      "PF",
                      "PE",
                      "EP",
                      "EN",
                      "BLC")]

datChip = datChip[datChip$Chr %in% chrs, ]
chip.gr = makeGRangesFromDataFrame(df = datChip,
                                   keep.extra.columns = T)

stages = c("iPSC", "DE", "GT", "PF", "PE", "EP", "EN", "BLC")

enhancers_atac_overlap = list()
enhancers_atac_overlap_fgwas = list()
enhancers_atac_overlap.gr = list()
enhancers_atac_overlap_fgwas_not_in_all_stages = list()
enhancers_atac_overlap_fgwas_not_in_all_stages.gr = list()
# subsetting ATAC-seq per stage
for (s in stages) {
  atac_subset = atac[atac[s] == 1,]
  peaks.gr = makeGRangesFromDataFrame(atac_subset, keep.extra.columns = T)
  
  datChip_subset = datChip[datChip[s] == 1,]
  
  chip.gr = makeGRangesFromDataFrame(df = datChip_subset,
                                     keep.extra.columns = T)
  
  # ATAC in ChIP
  enhancers_atac_overlap[[s]] = subsetByOverlaps(peaks.gr, chip.gr)
  enhancers_atac_overlap[[s]] = as.data.frame(enhancers_atac_overlap[[s]])
  enhancers_atac_overlap_fgwas[[s]] = enhancers_atac_overlap[[s]]
  enhancers_atac_overlap.gr [[s]] =  GRanges(enhancers_atac_overlap[[s]])
  enhancers_atac_overlap_fgwas[[s]]$stage = rep(s,nrow(enhancers_atac_overlap[[s]]))
  enhancers_atac_overlap_fgwas_not_in_all_stages[[s]] = enhancers_atac_overlap_fgwas[[s]][rowSums(enhancers_atac_overlap_fgwas[[s]][c(7:14)])<8,]
  enhancers_atac_overlap_fgwas_not_in_all_stages.gr[[s]] = GRanges(enhancers_atac_overlap_fgwas_not_in_all_stages[[s]])
}
enhancers_atac_overlap.gr = GRangesList(enhancers_atac_overlap.gr)

# enhancers_atac_overlap = do.call("rbind",enhancers_atac_overlap)
# enhancers_atac_overlap = unique(enhancers_atac_overlap )
# nrow(enhancers_atac_overlap) # 27909 enhancers by this simple overlap method
# 
# png(
#   paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps", '/upset_consensus_set_top15_enhancers_promoters.png', sep = ''),units = "in",res = 400,type = "cairo",
#   width = 13,
#   height = 8
# )
# 
# upset(as.data.frame(enhancers_atac_overlap[,c(7:14)]), order.by = "freq",nsets = 8,point.size = 2, line.size = 1, 
#       mainbar.y.label = "Number of active enhancers", sets.x.label = "Active enhancers per stage", 
#       text.scale = c(2, 2, 2, 2, 2, 2),
#       sets = rev(stages), keep.order = TRUE,
#       nintersects = 15)
# # Only plotting top 15 intersects
# 
# dev.off()

# to merge and for fgwas
for(s in stages){
  write.table(enhancers_atac_overlap_fgwas[[s]][enhancers_atac_overlap_fgwas[[s]]$stage==s,c(1:3)],
              file = paste0("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/merged/",s,"_enhancers_and_promoter.bed"),
              sep = "\t", row.names = F, col.names = F, quote = F)
}


# For fgwas 
enhancers_atac_overlap_fgwas = do.call("rbind",enhancers_atac_overlap_fgwas)
enhancers_atac_overlap_fgwas_not_in_all_stages = do.call("rbind",enhancers_atac_overlap_fgwas_not_in_all_stages)

### write resulting tables
# write.table(enhancers_atac_overlap[c(1:3,6:14)], "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/enhancers_promoters_atac_overlap_H3K27ac_binary.txt",
#             sep = "\t", row.names = F, quote = F)

write.table(enhancers_atac_overlap_fgwas[c(1:3,15)], "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/enhancers_promoters_atac_overlap_H3K27ac_fgwas.bed",
            sep = "\t", row.names = F, col.names = F, quote = F)
write.table(enhancers_atac_overlap_fgwas[c(1:3,15)], "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/enhancers_promoters_atac_overlap_H3K27ac_fgwas.bed",
            sep = "\t", row.names = F, col.names = F, quote = F)
write.table(enhancers_atac_overlap_fgwas_not_in_all_stages[c(1:3,15)], "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/enhancers_promoters_atac_overlap_H3K27ac_fgwas_not_all_stages.bed",
            sep = "\t", row.names = F, col.names = F, quote = F)

# Save the subset of those not shared across all
# id chr start end .
# no header, in separate files per stage
homer = enhancers_atac_overlap_fgwas_not_in_all_stages[c(1:3,6:15)]
homer$dot = rep(".",nrow(homer))
for(s in stages){
  
  write.table(homer[homer$stage == s, c("Name","seqnames","start","end","dot")], 
              paste0("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/for_homer/",s,
                     "_enhancers_atac_overlap_H3K27ac_homer_not_all_stages.bed"),
              sep = "\t", row.names = F, col.names = F, quote = F)
}

# Gene pathway analysis for all, and not shared across all

# Pathway analysis #   ########################


refseq <- makeTxDbFromUCSC(genome="hg19", tablename="refGene")

enhancers_atac_overlap_fgwas_not_in_all_stages.gr = GRangesList(enhancers_atac_overlap_fgwas_not_in_all_stages.gr)

overlap_not_shared_annotated = lapply(enhancers_atac_overlap_fgwas_not_in_all_stages.gr, annotatePeak, TxDb=refseq,
                                      tssRegion=c(-2000, 500), verbose=T,annoDb = "org.Hs.eg.db")
overlap_annotated = lapply(enhancers_atac_overlap.gr, annotatePeak, TxDb=refseq,
                           tssRegion=c(-2000, 500), verbose=T,annoDb = "org.Hs.eg.db")

# plotAnnoBar(overlap_not_shared_annotated, title = "Genomic features of overlaps (not shared in all)")
# 
# plotDistToTSS(overlap_not_shared_annotated,
#               title="TSS distribution of overlaps (not shared in all)")

genes = lapply(overlap_not_shared_annotated, function(i) as.data.frame(i)$geneId)
gene_names = lapply(genes,function(x) mapIds(org.Hs.eg.db, x,  'SYMBOL','ENTREZID'))

compEnrich <- compareCluster(geneCluster   = genes,
                             fun           = "enrichPathway",
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH")


png(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/pathway_enrichment/enrichment_reactome_top10_not_shared_in_all_stages_overlaps.png",
  units = "in",res = 400,type = "cairo",
  width = 12,
  height = 12
)

dotplot(compEnrich, showCategory = 10, title = "Pathway enrichment: overlaps not shared in all stages")

dev.off()

write.csv(compEnrich,"/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/pathway_enrichment/enrichment_reactome_not_shared_in_all_stages_overlaps_pathway.csv",
          row.names = F,quote = F)

## all


genes = lapply(overlap_annotated, function(i) as.data.frame(i)$geneId)
gene_names = lapply(genes,function(x) mapIds(org.Hs.eg.db, x,  'SYMBOL','ENTREZID'))

compEnrich <- compareCluster(geneCluster   = genes,
                             fun           = "enrichPathway",
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH")


png(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/pathway_enrichment/enrichment_reactome_top10_all_stages_overlaps.png",
  units = "in",res = 400,type = "cairo",
  width = 12,
  height = 12
)

dotplot(compEnrich, showCategory = 10, title = "Pathway enrichment: overlaps shared in all stages")

dev.off()

write.csv(compEnrich,"/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/pathway_enrichment/enrichment_reactome_all_stages_overlaps_pathway.csv",
          row.names = F,quote = F)


## save those that are only in enhancers!! to be comparable with ABC #########################################
##############################################################################################################



overlap_annotated
overlap_not_shared_annotated


png(
  paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps", '/upset_consensus_set_top15_enhancers.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 13,
  height = 8
)

upset(as.data.frame(enhancers_atac_overlap[,c(7:14)]), order.by = "freq",nsets = 8,point.size = 2, line.size = 1, 
      mainbar.y.label = "Number of active enhancers", sets.x.label = "Active enhancers per stage", 
      text.scale = c(2, 2, 2, 2, 2, 2),
      sets = rev(stages), keep.order = TRUE,
      nintersects = 15)
# Only plotting top 15 intersects

dev.off()


# For fgwas
enhancers_atac_overlap_fgwas = do.call("rbind",enhancers_atac_overlap_fgwas)
enhancers_atac_overlap_fgwas_not_in_all_stages = do.call("rbind",enhancers_atac_overlap_fgwas_not_in_all_stages)

### write resulting tables
write.table(enhancers_atac_overlap[c(1:3,6:14)], "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/enhancers_atac_overlap_H3K27ac_binary.txt",
            sep = "\t", row.names = F, quote = F)

write.table(enhancers_atac_overlap_fgwas[c(1:3,15)], "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/enhancers_atac_overlap_H3K27ac_fgwas.bed",
            sep = "\t", row.names = F, col.names = F, quote = F)
write.table(enhancers_atac_overlap_fgwas_not_in_all_stages[c(1:3,15)], "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/enhancers_atac_overlap_H3K27ac_fgwas_not_all_stages.bed",
            sep = "\t", row.names = F, col.names = F, quote = F)

# Save the subset of those not shared across all for HOMER
# id chr start end .
# no header, in separate files per stage
homer = enhancers_atac_overlap_fgwas_not_in_all_stages[c(1:3,6:15)]
homer$dot = rep(".",nrow(homer))
for(s in stages){
  
  write.table(homer[homer$stage == s, c("Name","seqnames","start","end","dot")], 
              paste0("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/for_homer/",s,
                     "_enhancers_atac_overlap_H3K27ac_homer_not_all_stages.bed"),
              sep = "\t", row.names = F, col.names = F, quote = F)
}

# Gene pathway analysis for all, and not shared across all

# Pathway analysis #   ########################


refseq <- makeTxDbFromUCSC(genome="hg19", tablename="refGene")

enhancers_atac_overlap_fgwas_not_in_all_stages.gr = GRangesList(enhancers_atac_overlap_fgwas_not_in_all_stages.gr)

overlap_not_shared_annotated = lapply(enhancers_atac_overlap_fgwas_not_in_all_stages.gr, annotatePeak, TxDb=refseq,
                                      tssRegion=c(-2000, 500), verbose=T,annoDb = "org.Hs.eg.db")
overlap_annotated = lapply(enhancers_atac_overlap.gr, annotatePeak, TxDb=refseq,
                           tssRegion=c(-2000, 500), verbose=T,annoDb = "org.Hs.eg.db")

# plotAnnoBar(overlap_not_shared_annotated, title = "Genomic features of overlaps (not shared in all)")
# 
# plotDistToTSS(overlap_not_shared_annotated,
#               title="TSS distribution of overlaps (not shared in all)")

genes = lapply(overlap_not_shared_annotated, function(i) as.data.frame(i)$geneId)
gene_names = lapply(genes,function(x) mapIds(org.Hs.eg.db, x,  'SYMBOL','ENTREZID'))

compEnrich <- compareCluster(geneCluster   = genes,
                             fun           = "enrichPathway",
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH")


png(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/pathway_enrichment/enrichment_reactome_top10_not_shared_in_all_stages_overlaps.png",
  units = "in",res = 400,type = "cairo",
  width = 12,
  height = 12
)

dotplot(compEnrich, showCategory = 10, title = "Pathway enrichment: overlaps not shared in all stages")

dev.off()

write.csv(compEnrich,"/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/pathway_enrichment/enrichment_reactome_not_shared_in_all_stages_overlaps_pathway.csv",
          row.names = F,quote = F)

## all


genes = lapply(overlap_annotated, function(i) as.data.frame(i)$geneId)
gene_names = lapply(genes,function(x) mapIds(org.Hs.eg.db, x,  'SYMBOL','ENTREZID'))

compEnrich <- compareCluster(geneCluster   = genes,
                             fun           = "enrichPathway",
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH")


png(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/pathway_enrichment/enrichment_reactome_top10_all_stages_overlaps.png",
  units = "in",res = 400,type = "cairo",
  width = 12,
  height = 12
)

dotplot(compEnrich, showCategory = 10, title = "Pathway enrichment: overlaps shared in all stages")

dev.off()

write.csv(compEnrich,"/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/pathway_enrichment/enrichment_reactome_all_stages_overlaps_pathway.csv",
          row.names = F,quote = F)


