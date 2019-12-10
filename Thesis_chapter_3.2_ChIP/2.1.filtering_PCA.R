### Filtering counts of consensus H327ac peaks
# Normalization, pearson cor and PCA

library(edgeR)
library(ggplot2)
library(reshape2)
library(plyr)

# from call to this script with snakemake
counts_file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks/merged/H3K27ac_counts_merged_all_samples_peaks_optimal.txt"
outdir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/peaks/trimmed/"


palette <- c("#CADAE8", "#7883BA", "#755A91", "#CC85B1", "#C15858", "#F4B8B0", "#96665A", "#6DA567")  # Diabetologia palette

#### function
# sort df by two strings 
# the order of the strings is important
# 1st by "first", 2nd by "second"


# mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

trim_by_stage = function(dge, CPM=CPM) {
  data = cpm(dge, normalized.lib.sizes = TRUE)
  data[data < CPM] = 0
  data[data > CPM] = 1 # convert to table where 0=0 and 1 equals any value above n
  
  c = c() # initiate vector of numbers, to save rows to take out of matrix
  
  for (r in 1:nrow(data)) {
    for (i in seq(1, 9, 3)) {
      # loop through stages iPSC to GT 3 by 3 (to jump over samples of same stage)
      if (sum(data[r, c(i, i + 1, i + 2)]) > 2) {
        break
      }
      else{
        for (i in seq(10, 11, 2)) {
          # pf stage

          if (sum(data[r, c(i, i + 1)]) > 1) {
            break
          }
          else {
            for (i in seq(12, 20, 3)) {

              if (sum(data[r, c(i, i + 1, i + 2)]) > 2) {
                break
              }
              else{
                for (i in seq(21, 22, 2)) {

                  # BLC stage
                  if (sum(data[r, c(i, i + 1)]) > 1) {
                    break
                  } 
                }
              }
            }
          }
          
        }
        
        
      }
    }

    if(i==21) {
      if (sum(data[r, c(i, i + 1)]) < 2)
        {
          c = c(c, r)  # if it gets to the end and there are no valid values, append to vector
        }
    }
  }
  
  rm(data)
  c= unique(c)
  return(c)
  
}
##

counts <-
  read.table(
    file = counts_file,
    header = T,
    check.names = F,
    skip = 1
  ) # file with feature counts & info of start and end of peaks (all peaks present in at least one sample as detected by MACS2)

sample = factor(c("SBAd2.1", "SBAd3.1", "SBNeo1.1"),
                levels = c("SBAd2.1", "SBAd3.1", "SBNeo1.1"))   # data comes from three donors

stage = factor(
  c("iPSC", "DE", "GT", "PF", "PE", "EP", "EN", "BLC"),
  levels = c("iPSC", "DE", "GT", "PF", "PE", "EP", "EN", "BLC")
) # 8 stages

rownames(counts) = counts[, 1]  # peaks are row names
genes = counts[,c(1:6)] # save info for later
counts = counts[, c(7:ncol(counts))] # getting just the counts

colnames(counts) = c(paste(rep(stage[1:3], each = 3), rep(sample, 3), sep = "-"),
                     "PF-SBAd2.1"  ,  "PF-SBAd3.1",
                     paste(rep(stage[5:7], each = 3), rep(sample, 3), sep = "-"),
                     "BLC-SBAd2.1",   "BLC-SBAd3.1" )


## Step 1: exploratory normalization

message(paste("################################## \n 
              The initial number of peaks is: \n", nrow(counts), "\n ##################################")) 

message("Plotting distribution of raw counts")

mdata = melt(counts)
#mdata_aggregated = count(mdata, c("variable", "value"))
mdata_aggregated = count(mdata, c( "value"))
colnames(mdata_aggregated) = c("raw_counts_per_peak", "number")

p1 = ggplot(mdata_aggregated, aes(x = raw_counts_per_peak, y = number)) +
  geom_line(size=2) +  scale_x_continuous(limits = c(0, 50), breaks=seq(0, 50, 5))+
  theme_minimal() +
  theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           panel.border = element_rect(fill = NA, colour = "black"),
           legend.key = element_blank(), 
           axis.title.y = element_text(face="bold", size=16),
           axis.title.x = element_blank(),
           axis.text.x = element_text(face="bold", colour = "black", size=16),
           axis.text.y = element_text(face="bold", colour = "black",size=16),
           axis.ticks = element_line(colour = "black"),
           axis.line = element_line(colour = "black"),
           legend.text = element_text(face="bold", colour = "black",size=14),
           legend.title = element_text(face="bold", colour = "black",size=16))

plot(p1)

# ggsave(
#   paste0(outdir,"H3K27ac_raw_distribution_small.png"),
#   p1,
#   width = 6,
#   height = 5,
#   units = "in",
#   dpi = 300
# )

p2 = ggplot(mdata_aggregated, aes(x = raw_counts_per_peak, y = number)) +
  geom_line(size=2) +  scale_x_continuous(limits = c(50, 1600), breaks=seq(50, 1600, 250))+
  theme_minimal() +
  theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           panel.border = element_rect(fill = NA, colour = "black"),
           legend.key = element_blank(), 
           axis.title.y = element_text(face="bold", size=16),
           axis.title.x = element_blank(),
           axis.text.x = element_text(face="bold", colour = "black", size=16),
           axis.text.y = element_text(face="bold", colour = "black",size=16),
           axis.ticks = element_line(colour = "black"),
           axis.line = element_line(colour = "black"),
           legend.text = element_text(face="bold", colour = "black",size=14),
           legend.title = element_text(face="bold", colour = "black",size=16))

plot(p2)
# ggsave(
#   paste0(outdir,"H3K27ac_raw_distribution_large.png"),
#   p2,
#   width = 6,
#   height = 5,
#   units = "in",
#   dpi = 300
# )


# Total counts vs total peaks

#The initial normalization of chip-seq data doesn't have any filtering. 
#Its purpose is to explore the behaviour of the counts and set a cutoff value for posterior filtering. 
#I'll use edgeR, using TMM normalization, as the distribution of counts per peak (which will be treated as counts per gene, 
#just like in RNA-seq data) is negative binomial (variance much larger than the mean). edgeR also manages better the presence of zero counts than limma+voom. 
#Now I will create the DGE object:


group <- factor(c(rep("iPSC",3),rep("DE",3),rep("GT",3),rep( "PF",2),rep( "PE",3),rep( "EP",3),rep( "EN",3),rep( "BLC",2)),levels = stage,ordered = T)   #group by stages
dge <- DGEList(counts = counts, group = group, genes = genes) #create dge object


#Calculate the normalization factors using TMM:

dge <- calcNormFactors(dge) # Calculating normalization factors, TMM normalization by default

CPMs = cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)  #calculate counts per million using the normalized library sizes

CPMs = as.data.frame(CPMs)

message("################################## \n Getting distribution of chip-seq CPMs \n ################################ ")

each_sample_counts_mean = colMeans(CPMs) # calculate mean CPMs per sample

each_peak_counts_mean = rowMeans(CPMs)
all_sample_counts_mean = mean(each_peak_counts_mean)    # calculate mean CPMs for all samples

appended_columns = unlist(CPMs, use.names = F)
all_sample_counts_median = median(appended_columns) # calculate median CPMs for all samples

all_sample_counts_variance = var(appended_columns) # calculate variance for CPMs for all samples
all_sample_counts_median = median(appended_columns)

mode = getmode(appended_columns)


# plotting frequency of CPMs


mdata <- melt(CPMs)

# Add stage info for plot
mdata$Stage = gsub("\\-.*","",mdata$variable)

# 
# # median per sample to plot, because calling function internally does median removing extreme values 
# median_per_sample = ddply(mdata, "variable",summarise, value = median(value))
# median_per_sample$Stage = factor(gsub("\\-.*","",median_per_sample$variable),
#                                  levels = stage,
#                                  ordered = T)
# 
# d1 = ggplot(mdata, aes(x = variable, y = value, colour = Stage)) +
#   geom_violin(trim = FALSE) +
#   geom_point(data=median_per_sample) + 
#   ylab("CPMs") + 
#   scale_color_manual(values = palette) +
#   ggtitle("Violin plot of CPMs per sample, with median (no y limits)") +
#   theme_minimal() +
#   theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#            panel.background = element_blank(),
#            panel.border = element_rect(fill = NA, colour = "black"),
#            legend.key = element_blank(), 
#            axis.title.y = element_text(face="bold", angle=90, size=16, vjust=0.2),
#            axis.title.x = element_blank(),
#            axis.text.x = element_text(face="bold", colour = "black", angle=90, size=16, vjust=0.2, hjust =1 ),
#            axis.text.y = element_text(face="bold", colour = "black",size=16),
#            axis.ticks = element_line(colour = "black"),
#            axis.line = element_line(colour = "black"),
#            legend.text = element_text(face="bold", colour = "black",size=14),
#            legend.title = element_text(face="bold", colour = "black",size=16))
# plot(d1)
# ggsave(file=large_violin_plot,
#        width=10,height=8,units="in",dpi=300)
# 
# d2 = ggplot(mdata, aes(x = variable, y = value, colour = Stage)) +
#   geom_violin(trim = FALSE, scale = "count") +
#   geom_point(data=median_per_sample) + 
#   scale_color_manual(values = palette) +
#   ggtitle("Violin plot of CPMs per sample, with median") +
#   scale_y_continuous(limits = c(0, 10), breaks = c(1:10)) + 
#   ylab("CPMs") + 
#   theme_minimal() +
#   theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#            panel.background = element_blank(),
#            panel.border = element_rect(fill = NA, colour = "black"),
#            legend.key = element_blank(), 
#            axis.title.y = element_text(face="bold", angle=90, size=16, vjust=0.2),
#            axis.title.x = element_blank(),
#            plot.title =  element_text(face="bold", size=16, vjust=0.2),
#            axis.text.x = element_text( colour = "black", angle=90, size=16, vjust=0.2, hjust =1 ),
#            axis.text.y = element_text(colour = "black",size=16),
#            axis.ticks = element_line(colour = "black"),
#            axis.line = element_line(colour = "black"),
#            legend.text = element_text( colour = "black",size=14),
#            legend.title = element_text(face="bold", colour = "black",size=16))
# plot(d2)
# ggsave(file=small_violin_plot,
#        width=10,height=8,units="in",dpi=300)

message(paste("############## \n The mean counts for all samples is: \n " , 
              round(all_sample_counts_mean,digits=2)))
message(paste("############## \n The variance of counts for all samples is: \n ", 
              round(all_sample_counts_variance,digits=2)))


message("############## \n Trimming by CPM and keeping only autosomal chromosomes \n ############## \n")

message("############## \n 1 CPM trimming \n ############## \n")
CPM = 1
original_dge = dge
original_CPMs = CPMs
# 2nd trim:
# go through the table 3 by 3
# if there are more than 2 values in three columns above 0, go to next row.
# if not, go to next 3
# if gets to end of loop and not found the values, get that row. After loop, subset matrix excluding those rows

c = trim_by_stage(dge,CPM=CPM)  # calling function above

ditch = dge$counts[c, ] # table to throw away

keep = !(rownames(dge$counts) %in% rownames(ditch))

dge <-
  dge[keep, , keep.lib.sizes = TRUE]

#keep only autosomal peaks
chr <- as.factor(paste("chr", c(1:22), sep=""))
dge <- dge[which(dge$genes$Chr %in% chr),, keep.lib.sizes = TRUE]

dim(dge)

saveRDS(dge, file = paste0(outdir,"dge_H3K27ac_1CPM_trim_optimal.xz") , compress = "xz",version = 2)   # saving the dge object

# save CPM table so I don't have to re-calculate CPMs on trimmed dge object (and therefore give lower CPMs than there really are!!

CPMs = CPMs[which(rownames(CPMs) %in% dge$genes$Geneid),]
CPMs = cbind(dge$genes[,2:4],CPMs)
write.table(CPMs, file = paste0(outdir,"CPMs_H3K27ac_1CPM_trim_optimal.txt"))

message(paste("############## \n Saved dge object after filtering for 1CPM has \n ", nrow(dge$counts), "peaks. ############## \n"))

# message(" ############## \n 10CPM trimming \n ############## \n")
# CPM = 10
# dge = original_dge
# CPMs = original_CPMs
# 
# 
# c = trim_by_stage(dge,CPM=CPM)  # calling function above
# 
# ditch = dge$counts[c, ] # table to throw away
# 
# keep = !(rownames(dge$counts) %in% rownames(ditch))
# 
# dge <-
#   dge[keep, , keep.lib.sizes = TRUE]
# 
# #keep only autosomal chromosomes
# chr <- as.factor(paste("chr", c(1:22), sep=""))
# dge<- dge[which(dge$genes$Chr %in% chr), keep.lib.sizes = TRUE]
# 
# dim(dge)
# 
# saveRDS(dge, file = paste0(outdir,"dge_H3K27ac_10CPM_trim_optimal.xz"), compress = "xz",version = 2)
# 
# CPMs = CPMs[which(rownames(CPMs) %in% dge$genes$Geneid),]
# CPMs = cbind(dge$genes[,2:4],CPMs)
# write.table(CPMs, file = paste0(outdir,"CPMs_H3K27ac_10CPM_trim_optimal.txt"))

message(paste("############## \n Saved dge object after filtering for 10CPMs has \n ", nrow(dge$counts), "peaks. ############## \n"))

message(paste("############## \n Plotting Pearson correlation matrices \ n ############## \n"))

#create the design matrix
stages= c(rep(as.character(stage[1:3]), each = 3),
          "PF"  ,  "PF",
          rep(as.character(stage[5:7]), each = 3),
          "BLC",   "BLC" )

samples=c(rep(as.character(sample[1:3]), 3),
          "SBAd2.1" , "SBAd3.1" ,
          rep(as.character(sample[1:3]), 3),
          "SBAd2.1",  "SBAd3.1" )
design <- model.matrix( ~ stages + samples)
design

# This converts the counts to log-counts per million with associated precision weights. After this, he RNA-seq data can be analyzed as if it was microarray data.

#voom= voom(original_dge, design = design, plot = TRUE)
voom = voom(dge, design = design, plot = TRUE)
#voom_10CPM = voom(dge_10CPM, design = design, plot = TRUE)

# # pearson correlation matrix


pearson_cor = function(x){
  
  
  cormat <- round(cor(x,method = "pearson"),2)
  # reorder_cormat <- function(cormat) {
  #   # Use correlation between variables as distance
  #   dd <- as.dist((1 - cormat) / 2)
  #   hc <- hclust(dd)
  #   cormat <- cormat[hc$order, hc$order]
  # }
  #cormat = reorder_cormat(cormat)
  
  melted_cormat <- melt(cormat)
  melted_cormat$Var1 = factor(melted_cormat$Var1,levels = unique(melted_cormat$Var1), ordered = T)
  melted_cormat$Var2 = factor(melted_cormat$Var2,levels = unique(melted_cormat$Var1), ordered = T)
  
  ggheatmap <- ggplot(melted_cormat, aes(Var1, Var2, fill = value)) +
    geom_tile(aes(fill = value)) + 
    scale_fill_gradientn(
      limit = c(0.0, 1.0),
      values = c(round(0.00,digits = 0), round(0.25,digits = 2),round(0.50,digits = 1),round(0.75,digits = 2),round(1.00,digits = 0)),
      colors = c("#faffff", "#002147"),
      space = "Lab",
      name = "Pearson\nCorrelation",
      position = "right",
      expand = expand_scale(add = .6)
    ) +
    scale_x_discrete(breaks = levels(melted_cormat$Var1)) +
    scale_y_discrete(breaks = levels(melted_cormat$Var1)) +
    theme_minimal() + # minimal theme
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(
        size = 11,
        hjust = 1,
        angle = 90
      ),
      legend.text = element_text(colour = "black", size = 11, margin(1)),
      legend.title = element_text(colour = "black", size = 12, face = "bold"),
      axis.text.y =  element_text( size = 11),
      legend.direction = "horizontal",
      legend.position = "top") +
    coord_fixed()
  
  
  plot(ggheatmap)
  write.table(cormat,paste0(outdir,"Pearson_cor_matrix.txt"))
}

png(
  paste0(outdir,"Pearson_cor_chip-seq_1CPM_trim.png"),
  type = "cairo",
  width = 9,
  height = 6,
  units = "in",
  res = 200,
  pointsize = 13
)
pearson_cor(voom$E)
dev.off()

pearson_cor(voom$E[,c(1:3)])

message(paste("############## \n Plotting PCA: filtered data\ n ############## \n"))

plot_pca = function(x, s = samples, st = stages) {
  pca1 <- prcomp(t(x), retx = TRUE, scale. = F)
  
  plot(pca1, type = "l") #variance vs first 10 components
  summary(pca1)          #importance of each component (important line is "proportion of variance")
  
  percentVar <- (pca1$sdev) ^ 2 / sum(pca1$sdev ^ 2)
  percentVar <- round(100 * percentVar)
  pcs <- as.data.frame(pca1$x)
  pcs <- cbind(pcs, sample = s, stage = st)
  pcs$stage <- ordered(pcs$stage, levels = stage)
  
  diaPalette <-
    c(
      "#000000",
      "#CADAE8",
      "#7883BA",
      "#755A91",
      "#CC85B1",
      "#F4B8B0",
      "#96665A",
      "#96165A"
    )  # Diabetologia pallete
  p <- ggplot(pcs, aes(PC1, PC2, colour = stage, shape = s)) +
    geom_point(size = 4) + xlab(paste0("PC1:" , percentVar[1], "% variance")) +
    scale_color_manual(values = diaPalette, name  ="Stage") +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_shape_discrete(name  ="Sample")
  p <-
    p + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black"),
      legend.key = element_blank(),
      # legend.position = c(0.5,0.5),
      axis.title.y = element_text(
        angle = 90,
        size = 18,
        vjust = 0.2,
        face = "bold"
      ),
      axis.title.x = element_text(
        size = 18,
        vjust = 0,
        face = "bold"
      ),
      axis.text.x = element_text(
        
        colour = "black",
        angle = 90,
        size = 16,
        vjust = 0.2,
        hjust = 1
      ),
      
      axis.text.y = element_text(colour = "black", size = 16),
      legend.text = element_text(colour = "black", size = 16),
      legend.title = element_text(colour = "black", size = 16, face = "bold"),
      axis.ticks = element_line(colour = "black"),
      axis.line = element_line(colour = "black")
    )
  return(p)
}

p = plot_pca(voom$E)
p
ggsave(
  paste0(outdir,"PCA_1CPM.png"),
  p,
  width = 8,
  height = 6,
  units = "in",
  dpi = 400
)

p2 = plot_pca(voom$E[,c(1:14,16:22)],s = samples[c(1:14,16:22)] ,st = stages[c(1:14,16:22)])
p2

ggsave(
  paste0(outdir,"PCA_1CPM_noEP2.1.png"),
  p2,
  width = 8,
  height = 6,
  units = "in",
  dpi = 400
)

