## fgwas plots
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)

inFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/fgwas_output/"
fgwas.out.dir = inFolder


"%&%" <- function(a,b) paste0(a,b)

# Annotation plot



an_plot_low_high <- function(mydf,
                    mytitl = "",
                    mylow = -10,
                    myhigh = 10,
                    interval = 1) {
  #order by log2FE
  
  mydf <-
    within(mydf, annotation <-
             factor(mydf$annotation, levels = rev(unique(mydf$annotation))))
  plt <- ggplot(data = mydf, aes(x = annotation, y = estimate, fill = methylation, group = methylation)) +
    geom_hline(yintercept = 0, linetype = 2) +
    ylab("Log2FE") + xlab("Annotation") +
    geom_errorbar(aes(ymin = CI_lo, ymax = CI_hi,),stat = "identity", width = 0.3, position = position_dodge(.4)) +
    geom_point(shape = 21,
               size = 2.5,
               col = "black",
               position = position_dodge(.4))  +
     scale_fill_manual( labels = c("hypermethylated", "hypomethylated"),values = c( "firebrick1","steelblue1")) +  # hypermethylated in red
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


DMP = read.table(paste0(inFolder,"DMP/DMP-individual-annotation-results.txt"),header = T)
DMP$annotation = sub('\\_.*', '', DMP$parameter)
DMP$methylation = gsub('.*_', '', DMP$parameter)
DMP = DMP[order(DMP$estimate,decreasing = T),]

DMR = read.table(paste0(inFolder,"DMR/-individual-annotation-results.txt"), header = T)
DMR$annotation = sub('\\_.*', '', DMR$parameter)
DMR$methylation = gsub('.*_', '', DMR$parameter)
DMR = DMR[order(DMR$estimate,decreasing = T),]


HMR = read.table(paste0(inFolder,"HMR/-individual-annotation-results.txt"), header = T)
HMR$annotation = sub('\\_.*', '', HMR$parameter)
HMR$methylation = "hypermethylated"
HMR = HMR[order(HMR$estimate,decreasing = T),]


LMR = read.table(paste0(inFolder,"LMR/-individual-annotation-results.txt"), header = T)
LMR$annotation = sub('\\_.*', '', LMR$parameter)
LMR$methylation = "hypomethylated"
LMR = LMR[order(LMR$estimate,decreasing = T),]



plt1 <-
  an_plot_low_high(
    DMP,
    mytitl = "fGWAS enrichment of DMPs in T2D SNPs",
    mylow = -30,
    myhigh = 6,
    interval = 5
  )
plt1


plt2 <-
  an_plot_low_high(
    DMR,
    mytitl = "fGWAS enrichment of DMRs in T2D SNPs",
    mylow = -30,
    myhigh = 5,
    interval = 5
  )
plt2

plt3 <-
  an_plot_low_high(
    rbind(LMR,HMR),
    mytitl = "fGWAS enrichment of LMRs and HMRs in T2D SNPs",
    mylow = -2,
    myhigh = 3,
    interval = 1
  )
plt3

ggsave(
  filename = paste0(fgwas.out.dir,"DMP_fgwas_enrich.separate.png"),
  plot = plt1,
  width = 6,
  height = 6,
  dpi = 320
)
ggsave(
  filename = paste0(fgwas.out.dir,"DMR_fgwas_enrich.separate.png"),
  plot = plt2,
  width = 6,
  height = 6,
  dpi = 320
)
ggsave(
  filename = paste0(fgwas.out.dir,"LMR_HMR_fgwas_enrich.separate.png"),
  plot = plt3,
  width = 6,
  height = 6,
  dpi = 320
)

## subsets of annotations

DMR_hypo_intergenic = read.table(paste0(inFolder,"DMR_hypo_intergenic/-individual-annotation-results.txt"), header = T)
DMR_hypo_intergenic$annotation = sub('\\_.*', '', DMR_hypo_intergenic$parameter)
DMR_hypo_intergenic$methylation = gsub('.*_', '', "hypomethylated")
DMR_hypo_intergenic = DMR_hypo_intergenic[order(DMR_hypo_intergenic$estimate,decreasing = T),]


DMR_hypo_promoter = read.table(paste0(inFolder,"DMR_hypo_promoter/-individual-annotation-results.txt"), header = T)
DMR_hypo_promoter$annotation = sub('\\_.*', '', DMR_hypo_promoter$parameter)
DMR_hypo_promoter$methylation = gsub('.*_', '', "hypomethylated")
DMR_hypo_promoter = DMR_hypo_promoter[order(DMR_hypo_promoter$estimate,decreasing = T),]

DMR_hypo_not_promoter = read.table(paste0(inFolder,"DMR_hypo_not_promoter/-individual-annotation-results.txt"), header = T)
DMR_hypo_not_promoter$annotation = sub('\\_.*', '', DMR_hypo_not_promoter$parameter)
DMR_hypo_not_promoter$methylation = gsub('.*_', '', "hypomethylated")
DMR_hypo_not_promoter = DMR_hypo_not_promoter[order(DMR_hypo_not_promoter$estimate,decreasing = T),]

DMR_hyper_intergenic = read.table(paste0(inFolder,"DMR_hyper_intergenic/-individual-annotation-results.txt"), header = T)
DMR_hyper_intergenic$annotation = sub('\\_.*', '', DMR_hyper_intergenic$parameter)
DMR_hyper_intergenic$methylation = gsub('.*_', '', "hypermethylated")
DMR_hyper_intergenic = DMR_hyper_intergenic[order(DMR_hyper_intergenic$estimate,decreasing = T),]


DMR_hyper_promoter = read.table(paste0(inFolder,"DMR_hyper_promoter/-individual-annotation-results.txt"), header = T)
DMR_hyper_promoter$annotation = sub('\\_.*', '', DMR_hyper_promoter$parameter)
DMR_hyper_promoter$methylation = gsub('.*_', '', "hypermethylated")
DMR_hyper_promoter = DMR_hyper_promoter[order(DMR_hyper_promoter$estimate,decreasing = T),]

DMR_hyper_not_promoter = read.table(paste0(inFolder,"DMR_hyper_not_promoter/-individual-annotation-results.txt"), header = T)
DMR_hyper_not_promoter$annotation = sub('\\_.*', '', DMR_hyper_not_promoter$parameter)
DMR_hyper_not_promoter$methylation = gsub('.*_', '', "hypermethylated")
DMR_hyper_not_promoter = DMR_hyper_not_promoter[order(DMR_hyper_not_promoter$estimate,decreasing = T),]

plt4 <-
  an_plot_low_high(
    rbind(DMR_hypo_intergenic,DMR_hyper_intergenic),
    mytitl = "fGWAS enrichment of intergenic DMRs in T2D SNPs",
    mylow = -30,
    myhigh = 5,
    interval = 5
  )
plt4
# nothing significant

plt5 <-
  an_plot_low_high(
    rbind(DMR_hypo_not_promoter,DMR_hyper_not_promoter),
    mytitl = "fGWAS enrichment of non-promoter DMRs in T2D SNPs",
    mylow = -30,
    myhigh = 5,
    interval = 5
  )
plt5

plt6 <-
  an_plot_low_high(
    rbind(DMR_hypo_promoter,DMR_hyper_promoter),
    mytitl = "fGWAS enrichment of promoter DMRs in T2D SNPs",
    mylow = -30,
    myhigh = 5,
    interval = 5
  )
plt6

ggsave(
  filename = paste0(fgwas.out.dir,"DMR_intergenic_fgwas_enrich.separate.png"),
  plot = plt4,
  width = 6,
  height = 6,
  dpi = 320
)
ggsave(
  filename = paste0(fgwas.out.dir,"DMR_not_promoter_fgwas_enrich.separate.png"),
  plot = plt5,
  width = 6,
  height = 6,
  dpi = 320
)

ggsave(
  filename = paste0(fgwas.out.dir,"DMR_promoter_fgwas_enrich.separate.png"),
  plot = plt6,
  width = 6,
  height = 6,
  dpi = 320
)




#### subset by other features
an_plot_features<- function(mydf,
                             mytitl = "",
                             mylow = -10,
                             myhigh = 10,
                             interval = 1) {
  #order by log2FE
  
  mydf <-
    within(mydf, annotation <-
             factor(mydf$annotation, levels = rev(unique(mydf$annotation))))
  plt <- ggplot(data = mydf, aes(x = annotation, y = estimate, fill = CpG, group = CpG)) +
    geom_hline(yintercept = 0, linetype = 2) +
    ylab("Log2FE") + xlab("Annotation") +
    geom_errorbar(aes(ymin = CI_lo, ymax = CI_hi,),stat = "identity", width = 0.3, position = position_dodge(.4)) +
    geom_point(shape = 21,
               size = 2.5,
               col = "black",
               position = position_dodge(.4))  +
    scale_fill_manual( labels = c("CpG islands", "CpG shores","CpG shelves","CpG open sea"),
                       values = c( "red","orange","green","blue")) +  
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


LMR_CpGislands = read.table(paste0(inFolder,"LMR_CpGislands/-individual-annotation-results.txt"), header = T)
LMR_CpGislands$annotation = sub('\\_.*', '', LMR_CpGislands$parameter)
LMR_CpGislands$CpG = gsub('.*_', '', "CpG islands")
LMR_CpGislands = LMR_CpGislands[order(LMR_CpGislands$estimate,decreasing = T),]


LMR_CpGshores= read.table(paste0(inFolder,"LMR_CpGshores/-individual-annotation-results.txt"), header = T)
LMR_CpGshores$annotation = sub('\\_.*', '', LMR_CpGshores$parameter)
LMR_CpGshores$CpG = gsub('.*_', '', "CpG shores")
LMR_CpGshores= LMR_CpGshores[order(LMR_CpGshores$estimate,decreasing = T),]

LMR_CpGshelves= read.table(paste0(inFolder,"LMR_CpGshelves/-individual-annotation-results.txt"), header = T)
LMR_CpGshelves$annotation = sub('\\_.*', '', LMR_CpGshelves$parameter)
LMR_CpGshelves$CpG = gsub('.*_', '', "CpG shelves")
LMR_CpGshelves = LMR_CpGshelves[order(LMR_CpGshelves$estimate,decreasing = T),]

LMR_CpGopensea = read.table(paste0(inFolder,"LMR_CpGopensea/-individual-annotation-results.txt"), header = T)
LMR_CpGopensea$annotation = sub('\\_.*', '', LMR_CpGopensea$parameter)
LMR_CpGopensea$CpG = gsub('.*_', '', "CpG open sea")
LMR_CpGopensea = LMR_CpGopensea[order(LMR_CpGopensea$estimate,decreasing = T),]

plt7 <-
  an_plot_features(
    rbind(LMR_CpGislands,LMR_CpGshores,LMR_CpGshelves,LMR_CpGopensea),
    mytitl = "fGWAS enrichment of LMRs at CpG features",
    mylow = -1,
    myhigh = 3,
    interval = 1
  )
plt7

ggsave(
  filename = paste0(fgwas.out.dir,"LMR_CpGfeatures_fgwas_enrich.separate.png"),
  plot = plt7,
  width = 6,
  height = 6,
  dpi = 320
)

### LMR promoter or not promoter


LMR_promoter = read.table(paste0(inFolder,"LMR_promoter/-individual-annotation-results.txt"), header = T)
LMR_promoter$annotation = sub('\\_.*', '', LMR_promoter$parameter)
LMR_promoter$methylation = gsub('.*_', '', "hypomethylated")
LMR_promoter = LMR_promoter[order(LMR_promoter$estimate,decreasing = T),]

LMR_not_promoter = read.table(paste0(inFolder,"LMR_not_promoter/-individual-annotation-results.txt"), header = T)
LMR_not_promoter$annotation = sub('\\_.*', '', LMR_not_promoter$parameter)
LMR_not_promoter$methylation = gsub('.*_', '', "hypomethylated")
LMR_not_promoter = LMR_not_promoter[order(LMR_not_promoter$estimate,decreasing = T),]



an_plot <- function(mydf,
                             mytitl = "",
                             mylow = -10,
                             myhigh = 10,
                             interval = 1) {
  #order by log2FE
  
  mydf <-
    within(mydf, annotation <-
             factor(mydf$annotation, levels = rev(unique(mydf$annotation))))
  plt <- ggplot(data = mydf, aes(x = annotation, y = estimate, fill = methylation, group = methylation)) +
    geom_hline(yintercept = 0, linetype = 2) +
    ylab("Log2FE") + xlab("Annotation") +
    geom_errorbar(aes(ymin = CI_lo, ymax = CI_hi,),stat = "identity", width = 0.3, position = position_dodge(.4)) +
    geom_point(shape = 21,
               size = 2.5,
               col = "black",
               position = position_dodge(.4))  +
    scale_fill_manual( labels = c("hypomethylated"),values = c( "steelblue1")) +  # hypermethylated in red
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

plt8 <-
  an_plot(
    LMR_promoter,
    mytitl = "fGWAS enrichment of promoter LMRs in T2D SNPs",
    mylow = -1,
    myhigh = 4,
    interval = 1
  )
plt8

plt9 <-
  an_plot(
    LMR_not_promoter,
    mytitl = "fGWAS enrichment of non-promoter LMRs in T2D SNPs",
    mylow = -30,
    myhigh = 5,
    interval = 5
  )
plt9

ggsave(
  filename = paste0(fgwas.out.dir,"LMR_promoter_fgwas_enrich.separate.png"),
  plot = plt8,
  width = 6,
  height = 6,
  dpi = 320
)

ggsave(
  filename = paste0(fgwas.out.dir,"LMR_not_promoter_fgwas_enrich.separate.png"),
  plot = plt9,
  width = 6,
  height = 6,
  dpi = 320
)

