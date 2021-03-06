# log-transformed histogram of read frequency (normalized) vs fragment size
# as in fig2 Buenrostro et al. 2013

library(ggplot2)
library(scales)

currentDate <- Sys.Date() # to save date in name of output files


basedir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/fragment_length_counts_filtered_bams/"

filelist <- list.files(basedir)
samples <- substring(filelist, 1, nchar(filelist)-30)

fragments = list()
for (s in samples) {
  fragments[[s]] = read_delim(paste0(basedir,s,".bam_fragment_length_count.txt"),col_names = F,delim = " ")
  colnames(fragments[[s]]) = c("count","fragment")
  fragments[[s]]$sampid = rep(s,nrow(fragments[[s]]))
  
  # normalize: count / total counts
  fragments[[s]]["count"] = fragments[[s]]["count"] / sum(fragments[[s]]["count"])
  
}

fragments = do.call("rbind", fragments)

# samples as ordered factors for plot
samples = paste(rep(c("iPSC","DE","GT","PF","PE","EP","EN","BLC"),each=3),rep(c("SBad2.1","SBad3.1","SBneo1.1"),8),sep = "_" )

toPlotlog = fragments
toPlotlog$sampid = factor(toPlotlog$sampid,levels = samples,ordered = T)
pal2 = c(rep("#002147",13) ,
         rep("#002147",2), "red","#002147","orange", # C
         rep("#002147",3),"pink",rep("#002147",2)) # D  
p <-
  ggplot(data = toPlotlog, aes(x = fragment, y = count, group = sampid, color = sampid)) +
  #ggtitle(unique(long$GeneName)) +
  # xlab("Differentiation stages") +
  # ylab("Normalized read density x ") +
  labs(x = "Fragment length (bp)", y = expression("Normalized read density")) +
  expand_limits(y = 0) +
  scale_x_continuous(limits = c(0, 1000)) +
  scale_color_manual(values = pal2) +
  
  
  scale_y_log10(breaks = trans_breaks("log10",  function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(scaled = T, sides = "l") +
  # geom_hline(
  #   yintercept = 0,
  #   linetype = "dashed",
  #   size = 1,
  #   col = "#DCDCDC"
  # ) +
  # geom_density() +
  geom_line(
    #   #aes(linetype = sample, col = sample), 
    size = 1) +
  #scale_colour_manual(values = diaPalette) +  # diabetologia pallete
  #scale_colour_manual(values="#000000") +  # for black and white, otherwise map lines and point colours to samples
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(size = 2),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold.italic"),
    legend.text = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 13, face = "bold")
  )

plot(p) 

ggsave(
  paste(
    "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/QC/fragment_length_vs_read_density",
    currentDate,
    ".jpg",
    sep = ""
  )
  ,
  p,
  width = 10,
  height = 5,
  units = "in",
  dpi = 400
)
