## total reads, read number and fraction of reads in peaks

# total mapped counts vs peak number

library(readr)
library(summarytools)
library(ggplot2)
library(grid)
library(gridExtra)
library(GenomicRanges)
library(ggrepel)
library(ChIPseeker)


basedir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks/"

currentDate <- Sys.Date() # to save date in name of output files


filelist <- list.files(basedir)
filelist = setdiff(filelist,c("merged","trimmed")) ## remove folders
narrowPeaks = list()

for (f in filelist) {
  s = substring(f, 1, nchar(f)-64)
  narrowPeaks[[s]] = read_tsv(paste0(basedir,f),col_names = F)
  narrowPeaks[[s]]$sampid = rep(s,nrow(narrowPeaks[[s]]))
  narrowPeaks[[s]]$names= paste(narrowPeaks[[s]]$X1,narrowPeaks[[s]]$X2,narrowPeaks[[s]]$X3, sep="_")
}

narrowPeaks = do.call('rbind', narrowPeaks)

# Load read mapping QC matrix
FRIP = read.table(file='/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/ENCODE_QC_aggregate_files/FRIP.tsv', 
                             header=TRUE, sep=' ')
full_stat_frame = read.table(file='/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/ENCODE_QC_aggregate_files/stat_frame.tsv', 
                             header=TRUE, sep=' ')

# Number of peaks per stage out of the Encode pipeline
freq(narrowPeaks$sampid)
# mean
mean(table(narrowPeaks$sampid))
sd(table(narrowPeaks$sampid))

# example number of duplicated peaks in one sample
sum(duplicated(narrowPeaks[narrowPeaks$sampid=="BLC-SBAd2.1","names"]))
# o
# Seems like there are no duplicates like in the kundaje atac pipeline


# Coverage plot: example
# peaks_cp = makeGRangesFromDataFrame(narrowPeaks[narrowPeaks$sampid=="BLC-SBAd2.1",],keep.extra.columns = TRUE,
#                                     seqnames.field = "X1",start.field = "X2",end.field = "X3")
# 
# png(
#   paste('/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/QC/covplot_EP_SB21.png', sep = ''),units = "in",res = 400,type = "cairo",
#   width = 12,
#   height = 11
# )
# covplot(peaks_cp,weightCol = "X5",title = paste("ATAC-seq peaks for BLC SBad2.1 over chromosomes"))
# 
# dev.off()
# 
# 
# # Coverage plot: iPSC example
# peaks_cp = makeGRangesFromDataFrame(narrowPeaks[narrowPeaks$sampid=="iPSC_SBneo1.1",],keep.extra.columns = TRUE,
#                                     seqnames.field = "X1",start.field = "X2",end.field = "X3")
# png(
#   paste('/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/QC/covplot_iPSC_Neo1.png', sep = ''),units = "in",res = 400,type = "cairo",
#   width = 12,
#   height = 11
# )
# covplot(peaks_cp,weightCol = "X5",title = paste("ATAC-seq peaks for iPSC SBneo1.1 over chromosomes"))
# 
# dev.off()

# view(dfSummary(reduced_peaks, plain.ascii = FALSE, style = "grid", 
#                graph.magnif = 0.75, valid.col = FALSE, 
#                tmp.img.dir = "../QC"))
# At most 24 reps of same peak name, indicative of peaks shared accross stages and not repeated within same stage and sample


reduced_peaks = as.data.frame(table(narrowPeaks$sampid))
colnames(reduced_peaks) = c("sampid","Number_of_peaks")
# merge info on number of peaks and reads
full_stat_frame = full_stat_frame[c("Type","rep","stage","filtered_mapped_dedup_total_read_pairs")]
full_stat_frame$sampid = paste(full_stat_frame$stage,full_stat_frame$rep,sep="-")
full_stat_frame = full_stat_frame[full_stat_frame$Type=="sample",]
full_stat_frame$pcnt_reads_in_peaks = FRIP$pcnt_reads_in_peaks
merged = merge(full_stat_frame[,c("sampid","filtered_mapped_dedup_total_read_pairs","pcnt_reads_in_peaks")],reduced_peaks)

#plot read pairs vs peaks

# area and color is pcnt of reads in peaks

p1 <-
  ggplot(merged,
         aes(filtered_mapped_dedup_total_read_pairs,Number_of_peaks)) +
  geom_point(aes(size = pcnt_reads_in_peaks, color = pcnt_reads_in_peaks)) +
  scale_color_continuous(
    high = "#5ad2f4",
    low = "#fe6847",
    limits = c(5, 30),
    breaks = seq(5, 30, by = 5)
  ) +
scale_size_continuous(limits = c(5, 30), breaks = seq(5, 30, by = 5)) +
guides(color = guide_legend("% of reads in peaks"),
       size = guide_legend("% of reads in peaks")) +
geom_text_repel(
  aes(label = sampid),
  size = 3,
  col = 'black',
  point.padding = 0.5
) +
labs(y = 'Total peaks',
     x = 'Total read pairs passing all QC') +

theme_bw() +
theme(
  panel.grid = element_blank(),
  axis.title.x = element_text(face = "bold"),
  axis.title.y = element_text(face = "bold"),
  legend.title = element_text(face = "bold"),
  legend.position="bottom"
)

plot(p1)
ggsave(
  paste(
    "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks_vs_reads",
    currentDate,
    ".jpg",
    sep = ""
  )
  ,
  p1,
  width = 6,
  height = 6,
  units = "in",
  dpi = 400
)
