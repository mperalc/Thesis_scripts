
library(edgeR)
library(ggplot2)
library(reshape2)
library(plyr)

# from call to this script with snakemake
#counts_file = snakemake@input[["counts_ATAC_all"]]
counts_file="data/counts/ATAC_counts_merged_all_samples_peaks_conservative_normalquality.txt"

large_violin_plot = snakemake@output[["large_violin_plot"]] 
small_violin_plot = snakemake@output[["small_violin_plot"]] 
ATAC_raw_count_distribution_large = snakemake@output[["ATAC_raw_count_distribution_large"]]
ATAC_raw_count_distribution_small = snakemake@output[["ATAC_raw_count_distribution_small"]]
counts_ATAC_1CPM = snakemake@output[["counts_ATAC_1CPM"]]
counts_ATAC_10CPM = snakemake@output[["counts_ATAC_10CPM"]]
CPMs_ATAC_1CPM = snakemake@output[["CPMs_ATAC_1CPM"]]
CPMs_ATAC_10CPM = snakemake@output[["CPMs_ATAC_10CPM"]]
PCA_1CPM = snakemake@output[["PCA_1CPM"]]
# 

palette <- c("#CADAE8", "#7883BA", "#755A91", "#CC85B1", "#C15858", "#F4B8B0", "#96665A", "#6DA567")  # Diabetologia palette

#### function
# sort df by two strings 
# the order of the strings is important
# 1st by "first", 2nd by "second"


sort_by_2_strings=function(df,first,second){
  order=c() # to store column order
  for (f in first){
    for(s in second){
      x=intersect(grep(paste("\\b",f,"\\b",sep=""), names(df)),
                  grep(paste("\\b",s,"\\b",sep=""), names(df))) # \\b are boundary anchors that allow to match whole words contained within a pattern
      # searching column numbers that match the two patterns
      order=c(order,x)
    }
  }
  df=df[,order]
  return(df)
}


# mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

trim_by_stage = function(dge, CPM=1) {
  data = cpm(dge, normalized.lib.sizes = TRUE)
  data[data < CPM] = 0
  data[data > CPM] = 1 # convert to table where 0=0 and 1 equals any value above n
  
  c = c() # initiate vector of numbers, to save rows to take out of matrix
  
  for (r in 1:nrow(data)) {
    for (i in seq(1, 24, 3)) {
      # loop through stages 3 by 3 (to jump over samples of same stage)
      
      if (sum(data[r, c(i, i + 1, i + 2)]) > 2) {
        break
      }
      if (i == 22) {
        c = c(c, r)  # if it gets to the end and there are no valid values, append to vector
      }
      
      
    }
    
  }
  
  rm(data)
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

sample = factor(c("SBad2.1", "SBad3.1", "SBneo1.1"),
levels = c("SBad2.1", "SBad3.1", "SBneo1.1"))   # data comes from three donors

stage = factor(
c("iPSC", "DE", "GT", "PF", "PE", "EP", "EN", "BLC"),
levels = c("iPSC", "DE", "GT", "PF", "PE", "EP", "EN", "BLC")
) # 8 stages

rownames(counts) = counts[, 1]  # peaks are row names
genes = counts[,c(1:6)] # save info for later
counts = counts[, c(7:ncol(counts))] # getting just the counts

colnames(counts) = paste(rep(stage, each = 3), rep(sample, 8), sep = "-")

# sort by stage and sample
counts = sort_by_2_strings(counts, stage, sample)

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

ggsave(
ATAC_raw_count_distribution_small,
p1,
width = 6,
height = 5,
units = "in",
dpi = 300
)
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
ggsave(
  ATAC_raw_count_distribution_large,
  p2,
  width = 6,
  height = 5,
  units = "in",
  dpi = 300
)


# Total counts vs total peaks

#The initial normalization of atac-seq data doesn't have any filtering. 
#Its purpose is to explore the behaviour of the counts and set a cutoff value for posterior filtering. 
#I'll use edgeR, using TMM normalization, as the distribution of counts per peak (which will be treated as counts per gene, 
#just like in RNA-seq data) is negative binomial (variance much larger than the mean). edgeR also manages better the presence of zero counts than limma+voom. 
#Now I will create the DGE object:


group <- rep(stage, each = 3)    #group by stages
dge <- DGEList(counts = counts, group = group, genes = genes) #create dge object
  

#Calculate the normalization factors using TMM:

dge <- calcNormFactors(dge) # Calculating normalization factors, TMM normalization by default

CPMs = cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)  #calculate counts per million using the normalized library sizes

CPMs = as.data.frame(CPMs)

message("################################## \n Getting distribution of ATAC-seq CPMs \n ################################ ")

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


d1 = ggplot(mdata, aes(x = variable, y = value, colour = Stage)) +
  #geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1) 
  geom_violin(trim = FALSE)
# ggplot(mdata_aggregated,aes(x=normalized_CPMs_per_peak, y=number, colour = samples)) +
#   geom_point() +  scale_x_continuous(limits=c(0,1000))
plot(d1)

# median per sample to plot, because calling function internally does median removing extreme values 
median_per_sample = ddply(mdata, "variable",summarise, value = median(value))
median_per_sample$Stage = factor(gsub("\\-.*","",median_per_sample$variable),
                                 levels = stage,
                                 ordered = T)

d1 = ggplot(mdata, aes(x = variable, y = value, colour = Stage)) +
  geom_violin(trim = FALSE) +
  geom_point(data=median_per_sample) + 
  ylab("CPMs") + 
  scale_color_manual(values = palette) +
  ggtitle("Violin plot of CPMs per sample, with median (no y limits)") +
  theme_minimal() +
  theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           panel.border = element_rect(fill = NA, colour = "black"),
           legend.key = element_blank(), 
           axis.title.y = element_text(face="bold", angle=90, size=16, vjust=0.2),
           axis.title.x = element_blank(),
           axis.text.x = element_text(face="bold", colour = "black", angle=90, size=16, vjust=0.2, hjust =1 ),
           axis.text.y = element_text(face="bold", colour = "black",size=16),
           axis.ticks = element_line(colour = "black"),
           axis.line = element_line(colour = "black"),
           legend.text = element_text(face="bold", colour = "black",size=14),
           legend.title = element_text(face="bold", colour = "black",size=16))
plot(d1)
ggsave(file=large_violin_plot,
       width=10,height=8,units="in",dpi=300)

d2 = ggplot(mdata, aes(x = variable, y = value, colour = Stage)) +
  geom_violin(trim = FALSE, scale = "count") +
  geom_point(data=median_per_sample) + 
  scale_color_manual(values = palette) +
  ggtitle("Violin plot of CPMs per sample, with median") +
  scale_y_continuous(limits = c(0, 10), breaks = c(1:10)) + 
  ylab("CPMs") + 
  theme_minimal() +
theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         panel.border = element_rect(fill = NA, colour = "black"),
         legend.key = element_blank(), 
         axis.title.y = element_text(face="bold", angle=90, size=16, vjust=0.2),
         axis.title.x = element_blank(),
         plot.title =  element_text(face="bold", size=16, vjust=0.2),
         axis.text.x = element_text( colour = "black", angle=90, size=16, vjust=0.2, hjust =1 ),
         axis.text.y = element_text(colour = "black",size=16),
         axis.ticks = element_line(colour = "black"),
         axis.line = element_line(colour = "black"),
         legend.text = element_text( colour = "black",size=14),
         legend.title = element_text(face="bold", colour = "black",size=16))
plot(d2)
ggsave(file=small_violin_plot,
       width=10,height=8,units="in",dpi=300)

message(paste("############## \n The mean counts for all samples is: \n " , 
              round(all_sample_counts_mean,digits=2)))
message(paste("############## \n The variance of counts for all samples is: \n ", 
              round(all_sample_counts_variance,digits=2)))


## Step 2: exploratory filtering

#Plot the distribution of CPM trimming and percentage of peaks that remain from the original 
#amount (from CPM 1 to 7 or so).


old_peaks = dim(dge)[1]# save how many peaks there are initially

remainder_after_trim_by_stage_and_cpm = function(dge,
                                                 CPM = 1,
                                                 old_peaks = dim(dge)[1]) {
 
    # go through the table 3 by 3
  # if there are more than 2 values in three columns above 0, go to next row.
  # if not, go to next 3
  # if gets to end of loop and not found the values, get that row. After loop, subset matrix excluding those rows
  c = trim_by_stage(dge,CPM=CPM)  # calling function above
  
  fail = 100 * (length(c) / nrow(dge$counts))
  peaks_fail = length(c)
  ditch = dge$counts[c, ] # table to throw away
  
  keep = !(rownames(dge$counts) %in% rownames(ditch))
  
  dge <-
    dge[keep, , keep.lib.sizes = TRUE]
  
  #keep only autosomal peaks
  chr <- as.factor(paste("chr", c(1:22), sep=""))
  dge <- dge[which(dge$genes$Chr %in% chr),, keep.lib.sizes = TRUE]
  
  message(paste("The counts object now has ", nrow(dge$counts), " peaks", sep = ""))
  
  rm(c, keep, ditch)
  
  peaks_remain = ((nrow(dge) / old_peaks) * 100)
  
  return(peaks_remain)
}

# call function for all CPMs
peaks_remain = c()
for (i in 1:10) {
  peaks_remain[i] = remainder_after_trim_by_stage_and_cpm(dge, CPM = i, old_peaks =
                                                            dim(dge)[1])
}


df = data.frame(CPM = c(1:10), peaks_remain = peaks_remain)


p <- ggplot(df, aes(y = peaks_remain, x = CPM)) +
  geom_line() + xlab("CPM cutoff") + ylab("% of peaks that remain") +
  scale_x_continuous(limits = c(1, 10), breaks = c(1:10)) + ylim(c(0, 100))
p

ggsave(file="data/plots/peaks_that_remain_after_CPM_cutoff.png",
       width=10,height=8,units="in",dpi=300)


#Save dge object after trimming:
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

save(dge, file = counts_ATAC_1CPM , compress = "xz")   # saving the dge object

# save CPM table so I don't have to re-calculate CPMs on trimmed dge object (and therefore give lower CPMs than there really are!!

CPMs = CPMs[which(rownames(CPMs) %in% dge$genes$Geneid),]
CPMs = cbind(dge$genes[,2:4],CPMs)
write.table(CPMs, file = CPMs_ATAC_1CPM)

message(paste("############## \n Saved dge object after filtering for 1CPM has \n ", nrow(dge_1CPM$counts), "peaks. ############## \n"))

message(" ############## \n 10CPM trimming \n ############## \n")
CPM = 10
dge = original_dge
CPMs = original_CPMs


c = trim_by_stage(dge,CPM=CPM)  # calling function above

ditch = dge$counts[c, ] # table to throw away

keep = !(rownames(dge$counts) %in% rownames(ditch))

dge <-
  dge[keep, , keep.lib.sizes = TRUE]

#keep only autosomal chromosomes
chr <- as.factor(paste("chr", c(1:22), sep=""))
dge<- dge[which(dge$genes$Chr %in% chr), keep.lib.sizes = TRUE]

dim(dge)

save(dge, file = counts_ATAC_10CPM, compress = "xz")

CPMs = CPMs[which(rownames(CPMs) %in% dge$genes$Geneid),]
CPMs = cbind(dge$genes[,2:4],CPMs)
write.table(CPMs, file = CPMs_ATAC_10CPM)

message(paste("############## \n Saved dge object after filtering for 10CPMs has \n ", nrow(dge_10CPM$counts), "peaks. ############## \n"))

message(paste("############## \n Plotting Pearson correlation matrices \ n ############## \n"))

#create the design matrix
stages=rep(stage,each=3)
samples=rep(sample,8)
design <- model.matrix( ~ stages + samples)
design

# This converts the counts to log-counts per million with associated precision weights. After this, he RNA-seq data can be analyzed as if it was microarray data.

#voom= voom(original_dge, design = design, plot = TRUE)
voom_1CPM = voom(dge_1CPM, design = design, plot = TRUE)
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
  write.table(cormat,"data/plots/Pearson_cor_matrix.txt")
}

png(
"data/plots/Pearson_cor_atac-seq_1CPM_trim.png",
type = "cairo",
width = 9,
height = 6,
units = "in",
res = 200,
pointsize = 13
)
pearson_cor(voom_1CPM$E)
dev.off()
# 
# png(
# "data/plots/Pearson_cor_atac-seq_10CPM_trim.png",
# type = "cairo",
# width = 7,
# height = 5,
# units = "in",
# res = 200,
# pointsize = 13
# )
# pearson_cor(voom_10CPM$E)
# 
# dev.off()

message(paste("############## \n Plotting PCA: 1CPM-filtered data\ n ############## \n"))

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

p = plot_pca(voom_1CPM$E)
p

ggsave(
  PCA_1CPM,
  p,
  width = 8,
  height = 6,
  units = "in",
  dpi = 350
)

#message(paste("############## \n Plotting PCA: 10CPM-filtered data\ n ############## \n"))
# 
# p = plot_pca(voom_10CPM$E)
# p

# Plotted these transfromations from unfiltered dge and is almost identical to the 1CPM filtered
