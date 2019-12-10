# overlap of ref data with promoter and other features

library(regioneR)

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

# Overlap H3K27ac and ATAC-seq and get

# Stats of overlaps at enhancers and promoter locations for ATAC-seq and H3K27ac and compare to publicly available data



folder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/ENCODE/processed/"
subfolders = c("A549","Adult_adrenal_gland","Adult_breast_epithelium","Adult_pancreas",
               "Adult_stomach","Adult_thyroid")
# Initialize lists
ATAC_list = list()
H3K27ac_list = list()

#ATAC
for (s in subfolders) {
  files = list.files(path = paste(folder, s, "/ATAC/bed_peaks", sep = ""))
  print(paste(s, "with files",files))
  for (f in files) {
    ATAC_list[[f]] = read.table(paste(folder, s, "/ATAC/bed_peaks/",f, sep = ""),
                                header = F,skip = 1)
    # skip header in gzipped, and track name in unzipped
    
  }
}

#ChIP
for (s in subfolders) {
  files = list.files(path = paste(folder, s, "/H3K27ac", sep = ""))
  for (f in files) {
    H3K27ac_list[[f]] = read.table(paste(folder, s, "/H3K27ac/",f, sep = ""),
                                   header = F,skip = 1)
    # skip header in gzipped, and track name in unzipped
    
  }
}

# Filter 1st column (Chr) so that it only contains chr 1=22 
# Make sure both ATAC and ChIP have similar filters for this comparison
chrs = paste('chr',rep(c(1:22),1),sep = "")

for (n in 1:10) {
  ATAC_list[[n]]$V1 = as.character(ATAC_list[[n]]$V1)
  H3K27ac_list[[n]]$V1 = as.character(H3K27ac_list[[n]]$V1)
  
  ATAC_list[[n]] =  ATAC_list[[n]][ATAC_list[[n]]$V1 %in% chrs,]
  H3K27ac_list[[n]] = H3K27ac_list[[n]][H3K27ac_list[[n]]$V1 %in% chrs,]
}


pt_chip = list()
# Make granges objects
for (n in 1:10) {
  ATAC_list[[n]] = makeGRangesFromDataFrame(
    df = ATAC_list[[n]],
    seqnames.field = "V1",
    start.field = "V2",
    end.field = "V3",
    keep.extra.columns = T
  )
  H3K27ac_list[[n]] = makeGRangesFromDataFrame(
    df = H3K27ac_list[[n]],
    seqnames.field = "V1",
    start.field = "V2",
    end.field = "V3",
    keep.extra.columns = T
  )
 # pt_chip[[n]] <- overlapPermTest(H3K27ac_list[[n]],  ATAC_list[[n]], ntimes=100, genome="hg19", count.once=TRUE, verbose = T)

}
pt_chip

peakAnnoList_atac <- lapply(ATAC_list, annotatePeak, TxDb=txdb,
                            tssRegion=c(-2000, 500), verbose=T)
peakAnnoList_chip <- lapply(H3K27ac_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-2000, 500), verbose=T)
names(peakAnnoList_chip) =c("A549-rep1","A549-rep2","A549-rep3","Adult_adrenal_gland","Adult_breast_epithelium","Adult_pancreas-rep1",
                            "Adult_pancreas-rep2","Adult_pancreas-rep3",
                            "Adult_stomach","Adult_thyroid")
names(peakAnnoList_atac) =  c("A549-rep1","A549-rep2","A549-rep3","Adult_adrenal_gland","Adult_breast_epithelium","Adult_pancreas-rep1",
                              "Adult_pancreas-rep2","Adult_pancreas-rep3",
                              "Adult_stomach","Adult_thyroid")

png(paste0(folder,"ATAC_genomic_annotations_publicData.png"),res = 400,height = 5,width = 6,units = "in")
plotAnnoBar(peakAnnoList_atac, title = "Distribution of ATAC-seq peaks relative to genomic features")
dev.off()

png(paste0(folder,"chip_genomic_annotations_publicData.png"),res = 400,height = 5,width = 6,units = "in")
plotAnnoBar(peakAnnoList_chip, title = "Distribution of ChIP-seq peaks relative to genomic features")
dev.off()


## my data and these are not directly comparable at the point of IDR, because these have way more (just 1 dataset, so no high confidence peaks). 