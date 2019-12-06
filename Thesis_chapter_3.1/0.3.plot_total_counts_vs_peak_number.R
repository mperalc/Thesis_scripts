# total mapped counts vs peak number

library(readr)
library(summarytools)
library(ggplot2)
library(grid)
library(gridExtra)
library(GenomicRanges)
library(ggrepel)
library(ChIPseeker)


basedir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/narrowPeaks/"

currentDate <- Sys.Date() # to save date in name of output files


filelist <- list.files(basedir)
stages <- substring(filelist, 1, nchar(filelist)-38)

narrowPeaks = list()

for (s in stages) {
  narrowPeaks[[s]] = read_tsv(paste0(basedir,s,".tn5.pval0.01.300K.bfilt.narrowPeak.gz"),col_names = F)
  narrowPeaks[[s]]$sampid = rep(s,nrow(narrowPeaks[[s]]))
  narrowPeaks[[s]]$names= paste(narrowPeaks[[s]]$X1,narrowPeaks[[s]]$X2,narrowPeaks[[s]]$X3, sep="_")
}

narrowPeaks = do.call('rbind', narrowPeaks)

# Load read mapping QC matrix
full_stat_frame = read.table(file='/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/QC/qc_summary.tsv', 
                             header=TRUE, sep='\t')


# Number of peaks out of Kundaje pipeline
freq(narrowPeaks$sampid)

# example number of duplicated peaks in one sample
sum(duplicated(narrowPeaks[narrowPeaks$sampid=="BLC_SBad2.1","names"]))
# This is for example same range with different summits
narrowPeaks[narrowPeaks$names=="chr11_57508674_57509096",]


# Coverage plot: EP example
peaks_cp = makeGRangesFromDataFrame(narrowPeaks[narrowPeaks$sampid=="BLC_SBad2.1",],keep.extra.columns = TRUE,
                                    seqnames.field = "X1",start.field = "X2",end.field = "X3")

png(
  paste('/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/QC/covplot_EP_SB21.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 12,
  height = 11
)
covplot(peaks_cp,weightCol = "X5",title = paste("ATAC-seq peaks for BLC SBad2.1 over chromosomes"))

dev.off()


# Coverage plot: iPSC example
peaks_cp = makeGRangesFromDataFrame(narrowPeaks[narrowPeaks$sampid=="iPSC_SBneo1.1",],keep.extra.columns = TRUE,
                                    seqnames.field = "X1",start.field = "X2",end.field = "X3")
png(
  paste('/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/QC/covplot_iPSC_Neo1.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 12,
  height = 11
)
covplot(peaks_cp,weightCol = "X5",title = paste("ATAC-seq peaks for iPSC SBneo1.1 over chromosomes"))

dev.off()


# Merge peaks within sample and stage
narrowPeaks_subset = narrowPeaks[c(1:3,5,11:12)]
colnames(narrowPeaks_subset) = c("chr","start","end","score","sampid","names")
reduced_peaks = list()
reduced_peaks_df = list()

for (s in stages) {
  reduced_peaks[[s]] = makeGRangesFromDataFrame(narrowPeaks_subset[narrowPeaks_subset$sampid == s,], keep.extra.columns = T)
  reduced_peaks[[s]] = reduce(reduced_peaks[[s]])
  reduced_peaks_df[[s]] = as.data.frame(reduced_peaks[[s]])
  reduced_peaks_df[[s]]$sampid = rep(s,nrow(reduced_peaks_df[[s]]))
  reduced_peaks_df[[s]]$names= paste(reduced_peaks_df[[s]]$seqnames,
                                  reduced_peaks_df[[s]]$start,
                                  reduced_peaks_df[[s]]$end, sep="_")
}



# rbind for plotting
reduced_peaks = do.call('rbind', reduced_peaks_df)

# Number of non-overlapping peaks
freq(reduced_peaks$sampid)

mean(table(reduced_peaks$sampid))
sd(table(reduced_peaks$sampid))

 # view(dfSummary(reduced_peaks, plain.ascii = FALSE, style = "grid", 
 #                graph.magnif = 0.75, valid.col = FALSE, 
 #                tmp.img.dir = "../QC"))
# At most 24 reps of same peak name, indicative of peaks shared accross stages and not repeated within same stage and sample


reduced_peaks = as.data.frame(table(reduced_peaks$sampid))
colnames(reduced_peaks) = c("sampid","Number_of_peaks")
# merge info on number of peaks and reads
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
    limits = c(10, 35),
    breaks = seq(10, 35, by = 5)
  ) +
  scale_size_continuous(limits = c(10, 35), breaks = seq(10, 35, by = 5)) +
  guides(color = guide_legend("% of reads in peaks"),
         size = guide_legend("% of reads in peaks")) +
  geom_text_repel(
    aes(label = sampid),
    size = 2.5,
    col = 'black',
    point.padding = 0.5
  ) +
  labs(y = 'Total peaks',
       x = 'Total read pairs passing all QC') +
  
  geom_vline(xintercept = 25000000, color="grey",linetype="dashed") +
  
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
    "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/QC/peaks_vs_reads",
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
