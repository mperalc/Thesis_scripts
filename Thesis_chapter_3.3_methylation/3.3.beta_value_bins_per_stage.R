# Bins of beta values per stage
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(dplyr)

currentDate <- Sys.Date() # to save date in name of output files

inFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/"
outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/"



# load matrix of normalized beta values
beta= read.csv(paste(inFolder,"quantile_normalised_beta_detP_0.01_nocrossreact.csv",sep=""), header=T,row.names = 1,check.names=F)

melted = melt(beta)
melted$stage = rep( c("iPSC","iPSC","iPSC","DE","DE","DE","GT","GT","GT","PF","PF","PF","PE","PE","PE","EP","EN","EN","EN","BLC","BLC","islet","islet","islet","islet","islet","islet","islet","islet","islet","islet","islet"),each=nrow(beta))
melted$stage = factor(melted$stage,levels = unique(melted$stage),ordered = T)
melted[melted$value<=0.25,"category"] = "<0.25"
melted[melted$value>0.25 & melted$value<=0.50,"category"] = "0.25-0.50"
melted[melted$value>0.50 & melted$value<=0.75,"category"] = "0.50-0.75"
melted[melted$value>0.75,"category"] = ">0.75"
melted$category = factor(melted$category,levels = c("<0.25","0.25-0.50","0.50-0.75",">0.75"),ordered = T)

#frequencies
df <- melted[c(3,4)] %>%
  group_by_all() %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)*100)

p = ggplot(df, aes(x = stage, y = n, fill = category)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format())+
  ggtitle("Beta value bins per stage") +
  xlab("Stage") +
  ylab("Percentage of total CpGs") +
  guides(fill=guide_legend(title="Beta value")) +
  theme_bw()

p

png(
  paste(outFolder, "beta_vals_bin_per_stage", currentDate, ".png", sep = ""),
  width = 6,
  height = 6,
  units = "in",
  res = 400,
  type = "cairo"
)

plot(p)

dev.off()