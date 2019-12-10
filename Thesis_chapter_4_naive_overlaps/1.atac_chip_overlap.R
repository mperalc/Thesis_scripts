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
enhancers_atac_overlap.gr = list()

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
}
enhancers_atac_overlap.gr = GRangesList(enhancers_atac_overlap.gr)







############################# from merged peaks to binary ###############################

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

stages = c("iPSC", "DE", "GT", "PF", "PE", "EP", "EN", "BLC")


peaks = read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/merged/merged_enhancers_and_promoters.bed",
                   header = F,
                   check.names = F) # file with info of start and end of peaks, and samples that have it
colnames(peaks) = c("Chr", "Start", "End", "Peaks_present")  # add header
peaks$Name = paste("peak", peaks$Chr, peaks$Start, sep = "_") # add name column needed for the binary function
peaks = peaks[c(5, 1:4)] # reorder

binary = peaks_merged_to_binary_df(peaks) # this function needs columns "Name","Chr","Start","End","Peaks_present", in that order
nrow(binary)

# sort columns
binary = binary[,c("Name","Chr","Start","End","iPSC","DE","GT","PF","PE","EP","EN","BLC")]
write.table(
  binary,
  file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/merged/binary_enhancers_and_promoters_ATAC_H3K27ac_overlaps.bed",
  sep = "\t",
  row.names = F,
  col.names = T,
  quote = F
)

png(
  paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps", '/upset_consensus_set_top15_enhancers_promoters.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 13,
  height = 8
)

upset(as.data.frame(binary[,c(5:12)]), order.by = "freq",nsets = 8,point.size = 2, line.size = 1,
      mainbar.y.label = "ATAC-seq overlap H3K27ac", sets.x.label = "ATAC-seq overlap H3K27ac",
      text.scale = c(2, 2, 2, 2, 2, 2),
      sets = rev(stages), keep.order = TRUE,
      nintersects = 15)
# Only plotting top 15 intersects

dev.off()

### getting only enhancers to make it comparable with ABC
binary.gr = GRanges(binary)
refseq <- makeTxDbFromUCSC(genome="hg19", tablename="refGene")

binary_annotated = annotatePeak(binary.gr,TxDb=refseq,
                                tssRegion=c(-2000, 500), verbose=T,annoDb = "org.Hs.eg.db")
annotations_not_promoter = setdiff(1:length(binary_annotated@anno$annotation),grep("Promoter",binary_annotated@anno$annotation))
binary_not_promoter = binary[annotations_not_promoter,]

## upsetplot enhancers (not promoter)
png(
  paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps", '/upset_consensus_set_top15_enhancers.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 13,
  height = 8
)

upset(as.data.frame(binary_not_promoter[,c(5:12)]), order.by = "freq",nsets = 8,point.size = 2, line.size = 1,
      mainbar.y.label = "ATAC-seq overlap H3K27ac", sets.x.label = "ATAC-seq overlap H3K27ac",
      text.scale = c(2, 2, 2, 2, 2, 2),
      sets = rev(stages), keep.order = TRUE,
      nintersects = 15)
# Only plotting top 15 intersects

dev.off()

sum(rowSums(binary_not_promoter[,c(5:12)])==8) ## no need to remove those shared across all: very few (85)

write.table(
  binary_not_promoter,
  file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/merged/binary_enhancers_ATAC_H3K27ac_overlaps.bed",
  sep = "\t",
  row.names = F,
  col.names = T,
  quote = F
)

## pathway enrichment
binary_not_promoter.gr = GRanges(binary_not_promoter)
binary_not_promoter_annotated = annotatePeak(binary_not_promoter.gr,TxDb=refseq,
                                tssRegion=c(-2000, 500), verbose=T,annoDb = "org.Hs.eg.db")

genes = list()
for(s in stages){
  subset = as.data.frame(binary_not_promoter_annotated@anno)
  genes[[s]] = subset[subset[s]==1,"geneId"]
  
}
gene_names = lapply(genes,function(x) mapIds(org.Hs.eg.db, x,  'SYMBOL','ENTREZID'))

compEnrich <- compareCluster(geneCluster   = genes,
                             fun           = "enrichPathway",
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH",readable=T)


png(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/pathway_enrichment/enrichment_reactome_top10_enhancers_overlaps.png",
  units = "in",res = 400,type = "cairo",
  width = 12,
  height = 12
)

dotplot(compEnrich, showCategory = 10, title = "Pathway enrichment: ATAC-seq and H3K27ac overlaps at enhancers")

dev.off()

write.csv(compEnrich,"/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/pathway_enrichment/enrichment_reactome_top10_enhancers_overlaps.csv",
          row.names = F,quote = F)

# analysing overlap with monogenic diabetes genes
monogenic = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Monogenic_diabetes/Monogenic_updated_2019.txt")
binary_not_promoter_annotated = as.data.frame(binary_not_promoter_annotated@anno)

sum(monogenic$V1 %in% binary_not_promoter_annotated$SYMBOL) / length(monogenic$V1)
#29/34
monogenic[monogenic$V1 %in% binary_not_promoter_annotated$SYMBOL,]


### HOMER enrichment
# format: chr start end feature

for(s in stages){
  subset = binary_not_promoter_annotated[binary_not_promoter_annotated[s]==1,c("seqnames","start","end","Name")]
  write.table(subset,paste0("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/for_homer/",s,"_enhancers_atac_overlap_H3K27ac_homer.bed"),
            row.names = F,quote = F, col.names = F, sep = "\t")
}

### fgws T2D enrichment
# format: chr start end feature
subset = list()
for(s in stages){
  subset[[s]] = binary_not_promoter_annotated[binary_not_promoter_annotated[s]==1,c("seqnames","start","end")]
  subset[[s]]$stage = rep(s, nrow(subset[[s]]))

}
subset = do.call("rbind",subset)
write.table(subset,paste0("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ATAC_H3K27ac_overlaps/for_fgwas/fGWAS_enhancers_atac_overlap_H3K27ac.bed"),
            row.names = F,quote = F, col.names = F, sep = "\t")

