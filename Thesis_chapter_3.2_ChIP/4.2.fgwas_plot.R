## fgwas plots
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)

stages <- c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ChIP-seq/fgwas/fgwas_output/"
fgwas.out.dir = outFolder



"%&%" <- function(a,b) paste0(a,b)

# Annotation plot



an_plot <- function(mydf,
                    mytitl = "",
                    mylow = -10,
                    myhigh = 10,
                    interval = 1) {
  #order by log2FE
  
  mydf <-
    within(mydf, annotation <-
             factor(mydf$annotation, levels = rev(unique(mydf$annotation))))
  plt <- ggplot(data = mydf, aes(x = annotation, y = estimate, fill = peaks, group = peaks)) +
    geom_hline(yintercept = 0, linetype = 2) +
    ylab("Log2FE") + xlab("Annotation") +
    geom_errorbar(aes(ymin = CI_lo, ymax = CI_hi,),stat = "identity", width = 0.3, position = position_dodge(.4)) +
    geom_point(shape = 21,
               size = 2.5,
               col = "black",
               position = position_dodge(.4))  +
    scale_fill_manual( values = c("#744253", "#f9b5ac")) +  
    theme_bw()  +  theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold")
    ) +
    scale_y_continuous(breaks = seq(mylow, myhigh, by = interval)) +
    coord_flip(ylim = c(mylow, myhigh)) +
    ggtitle(mytitl)
  return(plt)
}


all_df = read.table(paste0(outFolder,"all_peaks/-individual-annotation-results.txt"),header = T)
all_df$annotation = sub('\\_.*', '', all_df$parameter)
all_df$peaks = rep("all",nrow(all_df))
all_df = all_df[order(all_df$estimate,decreasing = T),]

specific_df = read.table(paste0(outFolder,"stage_specific/-individual-annotation-results.txt"), header = T)
specific_df$annotation = sub('\\_.*', '', specific_df$parameter)
specific_df$peaks = rep("stage-specific",nrow(specific_df))
specific_df = specific_df[order(specific_df$estimate,decreasing = T),]

df = rbind(all_df,specific_df)

plt1 <-
  an_plot(
    df,
    mytitl = "fGWAS enrichment of H3K27ac ChIP-seq peaks in T2D SNPs",
    mylow = -30,
    myhigh = 6,
    interval = 5
  )
plt1

ggsave(
  filename = paste0(fgwas.out.dir,"all_specific_H3K27ac_fgwas_enrich.separate.png"),
  plot = plt1,
  width = 6,
  height = 6,
  dpi = 320
)
