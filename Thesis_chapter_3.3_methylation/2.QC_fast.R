# Methylation QC pipeline (quick)

library(RColorBrewer)
library(minfi)
library(ENmix)
library(GEOquery)
library(wateRmelon) # for BMIQ normalization
library(preprocessCore)  # different normalization functions
library(IlluminaHumanMethylationEPICmanifest)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(gridExtra)
library(ChAMP)

currentDate <- Sys.Date() # to save date in name of output files
outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/"

`%notin%` <- function(x,y) !(x %in% y)

# Read IDAT files
baseDir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/P160281_MethylationEPIC_NicolaBeer/03.archive/P160281_MethylationEPIC_NicolaBeer.all.tar"
sheet = read.metharray.sheet(base="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/P160281_MethylationEPIC_NicolaBeer/03.archive/3.Results", recursive = F)
# add stage column

sheet$stage= c("iPSC","iPSC","iPSC","DE","DE","DE","GT","GT","GT","PF","PF","PF","PE","PE","PE","EP","EN","EN","EN","BLC","BLC","islet","islet","islet","islet","islet","islet","islet","islet","islet","islet","islet")

# add sample column
sheet$sample=c("SBNeo1.1","SBAd2.1","SBAd3.1","SBNeo1.1","SBAd2.1","SBAd3.1","SBNeo1.1","SBAd2.1","SBAd3.1","SBNeo1.1","SBAd2.1","SBAd3.1","SBNeo1.1","SBAd2.1","SBAd3.1","SBAd3.1","SBNeo1.1","SBAd2.1","SBAd3.1","SBNeo1.1","SBAd3.1","ISL","ISL","ISL","ISL","ISL","ISL","R","R","R","R","R")

sheet$Sample_Name=gsub( "^.*?_", "", sheet$Sample_Name) # take out stage part of experiment id (sample name)

sheet$simple_id = paste(
  c(
    "iPSC", "iPSC",  "iPSC",
    "DE",   "DE",   "DE",
    "GT",   "GT",   "GT",
    "PF",   "PF",   "PF",
    "PE",   "PE",   "PE",
    "EP",   "EN",   "EN",   "EN",
    "BLC",  "BLC",
    "124I",    "141 (14-19)",    "177",    "179",    "182",    "184",    "124R",    "136",    "137",    "138",    "139"  ),
  sheet$sample,
  sep = "-"
)
rownames(sheet) <- sheet$Sample_Name # add that as rownames
sheet$Basename <- file.path(baseDir, sub("_.*", "", sheet$Sample_Name),sheet$Sample_Name)


#read cross-reactive probes
cross_reactive_EPIC <- read.csv("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/cross-reactive_probes_Pidsley_2016/EPIC_cross_reactive_probes.csv",sep=";")
rgSet <- read.metharray.exp(targets = sheet) # Reads an entire methylation array experiment 

rgSet

  # assayData: 1052641 features, 32 samples 
  # array: IlluminaHumanMethylationEPIC
  # annotation: ilm10b2.hg19

head(sampleNames(rgSet)) # we see the samples are named following a standard IDAT naming convention with a 10 digit number


rgSet <- rgSet[,sheet$Sample_Name] # reordering before merging
#pData(rgSet) <- pD  # merging into pheno data of rgSet
pData(rgSet) <- as(sheet, "DataFrame") # merging into pheno data of rgSet
rgSet

# give more descriptive names to our samples:

sampleNames(rgSet) <- sheet$simple_id


############ probes QC  ####################

# transform to methylset

mraw <- preprocessRaw(rgSet) #  A MethylSet object contains only the methylated and unmethylated signals. 


##########modifying multifreqpoly###########

bincount <- function(x,breaks){
  x <- x[!is.na(x)]
  bc <- table(.bincode(x, breaks, TRUE, TRUE))
  temp=data.frame(id=c(1: (length(breaks)-1)))
  bc=data.frame(id=as.numeric(names(bc)),counts=as.numeric(bc))
  resu=merge(temp,bc,by="id",all.x=TRUE)
  resu$counts[is.na(resu$counts)]=0
  resu$counts[order(resu$id)]
}

# multifreqpoly <- function(mat, nbreaks=100, col=1:ncol(mat), xlab="", 
#                           ylab="Frequency", legend = list(x = "topright", fill=col,inset=c(-0.2,0),cex=0.6,
#                                                           legend = if(is.null(colnames(mat))) paste(1:ncol(mat)) else 
#                                                             colnames(mat)),...)
# {
#   if(!is.matrix(mat)) stop("Warning: input data is not a numeric matrix\n")
#   if(is.null(col)) col="black"
#   col=rep(col,ceiling(ncol(mat)/length(col)))
#   if(nbreaks > nrow(mat)) nbreaks=min(15,round(nrow(mat)/2))
#   
#   breaks <- seq(min(mat,na.rm=TRUE), max(mat,na.rm=TRUE), 
#                 diff(range(mat,na.rm=TRUE))/nbreaks)
#   mids <- 0.5 * (breaks[-1] + breaks[-length(breaks)])
#   counts <- sapply(data.frame(mat),bincount,breaks=breaks)
#   plot(range(mids),c(0,max(counts)),type="n",xlab=xlab,ylab=ylab,...)
#   for(i in 1:ncol(counts)){lines(mids,counts[,i],col=col[i],...)}
#   if(is.list(legend)) do.call(graphics::legend, legend)
# }

multifreqpoly <- function(mat, nbreaks=100, col=1:ncol(mat), xlab="Beta", renaming_sheet=sheet,
                          ylab="Frequency", main = "multifreqpoly", legend = TRUE,...)
{
  if(!is.matrix(mat)) stop("Warning: input data is not a numeric matrix\n")
  if(is.null(col)) col="black"
  # col=rep(col,ceiling(ncol(mat)/length(col)))
  if(nbreaks > nrow(mat)) nbreaks=min(15,round(nrow(mat)/2))
  
  breaks <- seq(min(mat,na.rm=TRUE), max(mat,na.rm=TRUE),
                diff(range(mat,na.rm=TRUE))/nbreaks)
  mids <- 0.5 * (breaks[-1] + breaks[-length(breaks)])
  counts <- sapply(data.frame(mat),bincount,breaks=breaks)
  df = as.data.frame(cbind(counts,breaks[2:length(breaks)]))
  # plot(range(mids),c(0,max(counts)),type="n",xlab=xlab,ylab=ylab,...)
  colnames(df)=c(sheet$simple_id,"Breaks")
  melted = melt(df,id.vars=rev(names(df))[1])
  colnames(melted)=c("Breaks","Samples","Frequency")
  # melted$Stages = sub("\\..*", "", melted$Samples)
  subdf = sheet[c("simple_id","stage")]
  colnames(subdf) = c("Samples","Stages")
  melted$Stages = subdf[match(melted$Samples, subdf$Samples),"Stages"]
  melted$Stages = factor(melted$Stages,levels = unique(melted$Stages),ordered = T)
  # for(i in 1:ncol(counts)){lines(mids,counts[,i],col=col[i],...)}
  # if(is.list(legend)) do.call(graphics::legend, legend)
  ggplot(melted,aes(x=Breaks,y=Frequency, group = Samples, col = Stages)) +
    geom_line(size=1) +
    xlab(xlab) +
    ggtitle(main) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black"),
      plot.title = element_text(
        size = 20,
        face = "bold"
      ),
      axis.title.y = element_text(
        angle = 90,
        size = 18,
        face = "bold"
      ),
      axis.title.x = element_text(
        size = 18,
        vjust = 0,
        face = "bold"
      ),
      axis.text.x = element_text(
        
        colour = "black",
        size = 16


      ),
      
      axis.text.y = element_text(colour = "black", size = 16),
      legend.text = element_text(colour = "black", size = 16),
      legend.title = element_text(colour = "black", size = 16, face = "bold"),
      axis.ticks = element_line(colour = "black"),
      axis.line = element_line(colour = "black")
    )
  
}
#####################

#total intensity plot is userful for data quality inspection
#and identification of outlier samples

#  A RatioSet object is a class designed to store Beta values and/or M values instead of the methylated and unmethylated signals. 

RSet <- ratioConvert(mraw, what = "both", keepCN = TRUE)
RSet
  # Preprocessing
  # Method: Raw (no normalization or bg correction)
  # minfi version: 1.18.6
  # Manifest version: 0.3.0

# The functions getBeta, getM and getCN return respectively the Beta value matrix, M value matrix and the Copy Number matrix.

 beta <- getBeta(RSet)
 beta.noQN = beta # for later
 p1 = multifreqpoly(beta,main="All probes",xlab="Beta value")
 
 anno = getProbeType(RSet)   # get the EPIC Annotation data: probe types I or II
 
 beta1 = beta[anno == "I", ]
 beta2 = beta[anno == "II", ]
 p2 = multifreqpoly(beta1,main="Infinium I probes",xlab="Beta value")
 p3 = multifreqpoly(beta2,main="Infinium II probes",xlab="Beta value")
 
 
 
 png("multifreq_plots_beta_values_raw.png",height=6,width=18,units = "in",res = 600,type = "cairo")
 
 print(grid.arrange(p1, p2, p3, nrow = 1))
 dev.off()

# # The function mapToGenome applied to a RatioSet object will add genomic coordinates to each probe together with some additional annotation 
# # information. The output object is a GenomicRatioSet 
# 
 GRset <- mapToGenome(RSet,mergeManifest=TRUE)
 GRset

#Poor performing probes are generally filtered out prior to differential methylation analysis. 
# As the signal from these probes is unreliable, by removing them we perform fewer statistical 
# tests and thus incur a reduced multiple testing penalty. We filter out probes that have failed 
# in one or more samples based on detection p-value.

# detP <- minfi::detectionP(rgSet) # calculate detection p-vals
# 
# # # % of probes in total and for each sample that have a detection p-value above 0.01:
# failed <- detP>=0.01
# colMeans(failed)*100 # % of failed positions per sample
# sum(rowMeans(failed)>=0.01) # How many positions failed in >1% of samples?
# # 
# # # remove poor quality samples
# 
# keep <- colMeans(detP) <0.01
# rgSet = rgSet[,keep]
# # remove from detection p-val table:
# detP <- detP[,keep]
# dim(detP) # no samples lost 08/2019
# 
# # # ensure probes are in the same order in the mSetSq and detP objects 
# detP <- detP[match(featureNames(GRset),rownames(detP)),]
# # remove any probes that have failed in one or more samples
# keep <- rowSums(detP < 0.01) == ncol(GRset)
# table(keep)
# (table(keep)[1]/table(keep)[2])*100 # % of failed probes in one or more samples
# 
# # GRsetFlt <- GRset[keep,]
# # GRsetFlt # lost 4324 probes (0.6%)

# remove probes from sex chromosomes
anno=getAnnotation(GRset)

# if your data includes males and females, remove probes on the sex chromosomes 
discard <- intersect(featureNames(GRset), anno$Name[anno$chr %in% c("chrX","chrY")])
length(discard) # 19627
(length(discard)/length(featureNames(GRset)))*100 # to lose 16927 probes (2.3%)

# removal of probes where common SNPs may affect the CpG. You can either remove all probes 
# affected by SNPs (default), or only those with minor allele frequencies greater than a 
# specified value.
# Retained because I may be interested in some of those if allele is associated with T2D

# remove probes with SNPs at CpG site 
# test= getSnpInfo(GRsetFlt, snpAnno = NULL)
# 
# GRsetFlt <- dropLociWithSnps(GRsetFlt,snpAnno = NULL) # is this using the annotation from the object?
# 
# GRsetFlt # lost 28453 probes (3.4%)

# Drop cross-hybridizing probes (aprox 5%).

probeNames <- featureNames(GRset)
discard2 <- intersect(probeNames,as.character(cross_reactive_EPIC[,1])) # keep only probes that are not cross reactive
length(discard2) # 43177
(length(discard2)/length(featureNames(GRset)))*100   #5%

discard=unique(c(discard,discard2)) # 61400
(length(discard)/length(featureNames(GRset)))*100  # 7.10% of the original probes will be lost

############################################ normalization

# QN function (Matthias')

##get detection p-values:

dp <- minfi::detectionP(rgSet, type = "m+u")

## Type II probes
TypeII.Name <- getProbeInfo(rgSet, type = "II")$Name
TypeII.Green <- getGreen(rgSet)[getProbeInfo(rgSet, type = "II")$AddressA,]
TypeII.Red <- getRed(rgSet)[getProbeInfo(rgSet, type = "II")$AddressA,]
rownames(TypeII.Red) <- TypeII.Name
colnames(TypeII.Red) <- sampleNames(rgSet)
rownames(TypeII.Green) <- TypeII.Name
colnames(TypeII.Green) <- sampleNames(rgSet)

## Type I probes, split into green and red channels
TypeI.Green.Name <- getProbeInfo(rgSet, type = "I-Green")$Name
TypeI.Green.M <- getGreen(rgSet)[getProbeInfo(rgSet, type = "I-Green")$AddressB,] # address for methylated signal=B
rownames(TypeI.Green.M) <- TypeI.Green.Name
colnames(TypeI.Green.M) <- sampleNames(rgSet)
TypeI.Green.U <- getGreen(rgSet)[getProbeInfo(rgSet, type = "I-Green")$AddressA,] # address for unmethylated signal=A. 
rownames(TypeI.Green.U) <- TypeI.Green.Name
colnames(TypeI.Green.U) <- sampleNames(rgSet)

TypeI.Red.Name <- getProbeInfo(rgSet, type = "I-Red")$Name
TypeI.Red.M <- getRed(rgSet)[getProbeInfo(rgSet, type = "I-Red")$AddressB,]
rownames(TypeI.Red.M) <- TypeI.Red.Name
colnames(TypeI.Red.M) <- sampleNames(rgSet)
TypeI.Red.U <- getRed(rgSet)[getProbeInfo(rgSet, type = "I-Red")$AddressA,]
rownames(TypeI.Red.U) <- TypeI.Red.Name
colnames(TypeI.Red.U) <- sampleNames(rgSet)

##remove high missingness probes
d = ifelse(dp<0.01,1,NA) # if p-value <0.01, make it 1. Else, make it NA.
cr = data.frame(rowSums(is.na(d))/length(d[1,])) # sums NAs in each row, then divide by nº of samples.
  # In the end, it excludes probes with p-values above 0.01 in a min of 1 sample. 
exclude.badcalls = rownames(cbind(cr,rownames(cr))[cbind(cr,rownames(cr))[,1]>0.01,])
# 4548 probes
# test=names(as.matrix(cr)[which(as.matrix(cr)[,1]>0.02),])  # does the same as exclude.badcalls

exclude.sites = unique(c(exclude.badcalls,discard))  #  4548 failed probes + previous ones = 65243
# notice we're not discarding them yet

##exclude.sites = unique(rbind(as.matrix(exclude.chrX), as.matrix(exclude.chrY),as.matrix(exclude.cas),as.matrix(exclude.snps),as.matrix(crossmap),as.matrix(exclude.badcalls), as.matrix(exclude.mhc)))

##remove high missingness samples 

mind = data.frame(colSums(is.na(d))/length(d[,1])) # sums NAs in each column, and divides by nº of probes
remove.mind = rownames(cbind(mind,rownames(mind))[cbind(mind,rownames(mind))[,1]>0.01,])  # remove samples with over 1% failed probes
# test2=names(as.matrix(mind)[which(as.matrix(mind)[,1]>0.02),])  # does the same as remove.mind
# no failed samples

samples = rownames(cbind(mind,rownames(mind))[cbind(mind,rownames(mind))[,1]<0.01,])  # samples with less than 1% failed probes.

TypeII.Green =subset(TypeII.Green, select=samples)
TypeII.Red = subset(TypeII.Red, select=samples)
TypeI.Green.M = subset(TypeI.Green.M, select=samples)
TypeI.Green.U = subset(TypeI.Green.U, select=samples)
TypeI.Red.M = subset(TypeI.Red.M, select=samples)
TypeI.Red.U = subset(TypeI.Red.U, select=samples)
  # Everything stays the same.

#--------------------------------------------------------------------------------------------------------------------------------

#############QN (quantile normalization)

#exclude sites (filtering out failed probes)
TypeII.Green = TypeII.Green[rownames(TypeII.Green) %notin% as.matrix(exclude.sites),]
TypeII.Red = TypeII.Red[rownames(TypeII.Red) %notin% as.matrix(exclude.sites),]
TypeI.Green.M = TypeI.Green.M[rownames(TypeI.Green.M) %notin% as.matrix(exclude.sites),]
TypeI.Green.U = TypeI.Green.U[rownames(TypeI.Green.U) %notin% as.matrix(exclude.sites),]
TypeI.Red.M = TypeI.Red.M[rownames(TypeI.Red.M) %notin% as.matrix(exclude.sites),]
TypeI.Red.U = TypeI.Red.U[rownames(TypeI.Red.U) %notin% as.matrix(exclude.sites),]

# normalizing using quantiles

TypeII.Red.norm=normalize.quantiles(TypeII.Red)
TypeII.Green.norm=normalize.quantiles(TypeII.Green)
TypeI.Green.M.norm=normalize.quantiles(TypeI.Green.M)
TypeI.Green.U.norm=normalize.quantiles(TypeI.Green.U)
TypeI.Red.M.norm=normalize.quantiles(TypeI.Red.M)
TypeI.Red.U.norm=normalize.quantiles(TypeI.Red.U)
#rm(TypeII.Red,TypeII.Green,TypeI.Green.M,TypeI.Green.U,TypeI.Red.M,TypeI.Red.U)

#calculate betas
TypeII.betas = TypeII.Green.norm/(TypeII.Red.norm+TypeII.Green.norm+100)
rownames(TypeII.betas)=rownames(TypeII.Green)
colnames(TypeII.betas)=colnames(TypeII.Green)
TypeI.Green.betas = TypeI.Green.M.norm/(TypeI.Green.M.norm+TypeI.Green.U.norm+100)
rownames(TypeI.Green.betas)=rownames(TypeI.Green.M)
colnames(TypeI.Green.betas)=colnames(TypeI.Green.M)
TypeI.Red.betas = TypeI.Red.M.norm/(TypeI.Red.M.norm+TypeI.Red.U.norm+100)
rownames(TypeI.Red.betas)=rownames(TypeI.Red.M)
colnames(TypeI.Red.betas)=colnames(TypeI.Red.M)

beta=rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas)


p1 = multifreqpoly(beta,main="All probes",xlab="Beta value",legend=F)
p2 = multifreqpoly(rbind(TypeI.Red.betas,TypeI.Green.betas),main="Infinium I",xlab="Beta value",legend=F)
p3 = multifreqpoly(TypeII.betas,main="Infinium II",xlab="Beta value",legend=F)



png("multifreqplots_filtered_and_QN_without_CR_probes.png",height=6,width=18,units = "in",res = 600,type = "cairo")

print(grid.arrange(p1, p2, p3, nrow = 1))
dev.off()

# most changes in type II probe distribution (in abs value). Also in type I.
# before and after filtering and normalization:

p1 = multifreqpoly(beta.noQN,main="Before filtering and QN",xlab="Beta value",legend=F)
p2 = multifreqpoly(beta,main="After",xlab="Beta value",legend=F)

jpeg("multifreq_plots_beta_values_Matthias_pipeline_removing_CRprobes_comparison_un-normalized.jpg",height=9,width=9,units = "in",res = 600,type = "cairo")

print(grid.arrange(p1, p2, nrow = 2))
dev.off()


beta <- as.matrix(beta)
rm(TypeII.Red,TypeII.Green,TypeI.Green.M,TypeI.Green.U,TypeI.Red.M,TypeI.Red.U)

##########end of QNorm 

dim(beta)

# save for diff methylation analysis
 write.csv(beta, "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/quantile_normalised_beta_detP_0.01_nocrossreact.csv", col.names=T,row.names=T, quote=F)

 # write sample info
 write.csv(sheet, "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/samples_info.csv", col.names=T,row.names=F, quote=F)
 
 
 # For PCA and other analyses
m=log2(beta/(1-beta))    #M=log2(Beta/(1-Beta))

# Compare with beta
p1 = multifreqpoly(beta,main="M values after filtering and normalization",xlab="Beta value",legend=F)
p2 = multifreqpoly(m,main="M values after filtering and normalization",xlab="M value",legend=F)


jpeg("multifreq_plots_M_and_beta_values_comparison.jpg",height=9,width=9,units = "in",res = 600,type = "cairo")

print(grid.arrange(p1, p2, nrow = 2))
dev.off()


Nas=which(is.na(m))  # which rows have NA? 
write.csv(m, "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/quantile_normalised_M_detP_0.01_nocrossreact.csv", col.names=T,row.names=T, quote=F)


# PCA after filtering and normalization

plot_pca=function(x,sample,stage){
  pca1<-prcomp(t(x), retx=TRUE,scale. = F)   # transpose so that samples are in rows, and probes in columns
  
  # plot(pca1, type = "l") #variance vs first 10 components
  summary(pca1)          #importance of each component (important line is "proportion of variance")
  
  percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar <- round(100 * percentVar)
  pcs <- as.data.frame(pca1$x)
  pcs <- cbind(pcs,sample=sample,stage=stage)
  pcs$sample=factor(pcs$sample,levels=unique(sample))
  pcs$stage=factor(pcs$stage,levels=unique(stage))
  # pcs$stage <- ordered(pcs$stage, levels = stage)
  
  p <- ggplot(pcs, aes(PC1, PC2, colour=stage, shape=sample)) + 
    geom_point(size = 3) + xlab (paste0( "PC1:" ,percentVar[ 1 ],"% variance")) + 
    ylab (paste0( "PC2: ",percentVar[ 2 ],"% variance" ))
  p <- p + theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                  panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, colour = "black"), 
                  legend.key = element_blank(),# legend.position = c(0.5,0.5),
                  axis.title.y = element_text(face="bold", angle=90, size=12, vjust=0.2),
                  axis.title.x = element_text(face="bold", size=12, vjust=0),
                  axis.text.x = element_text(face="bold", colour = "black", angle=90, size=12, vjust=0.2, hjust =1 ),
                  axis.text.y = element_text(face="bold", colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  axis.line = element_line(colour = "black"),
                  #legend.position = "bottom",
                  legend.text = element_text(size=7))
  
  return(p)
}


p1=plot_pca(m,sample=sheet$sample,stage=sheet$stage)
ggsave(
  paste(
    "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/PCA_after_filtering_and_QN",
    currentDate,
    ".jpg",
    sep = ""
  ),
  p1,
  width = 6,
  height = 5,
  units = "in",
  dpi = 400
)

# plot without adult islets (R and ISL)

p_nadults=plot_pca(m[,1:21],sample=sheet$sample[1:21],stage= sheet$stage[1:21] )
ggsave(
  paste(
    "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/PCA_after_filtering_and_QN_no_adult_islets",
    currentDate,
    ".jpg",
    sep = ""
  ),
  p_nadults,
  width = 6,
  height = 5,
  units = "in",
  dpi = 400
)

# Removing batch effects from SB3.1 outliers
batch <- sheet$sample[1:21]

batch_corrected = removeBatchEffect(m[,1:21],batch)

p_corrected=plot_pca(batch_corrected,sample=sheet$sample[1:21],stage= sheet$stage[1:21] )
ggsave(
  paste(
    "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/PCA_after_filtering_and_QN_corrected_noAdults",
    currentDate,
    ".jpg",
    sep = ""
  ),
  p_corrected,
  width = 6,
  height = 5,
  units = "in",
  dpi = 400
)
# If plotted with batch correction of islets: weird location of islets


# sample distance cluster


plot_sdc=function(x){
  sampleDists <- dist(t(x))   
  #This function computes and returns the distance matrix computed by using the specified distance measure to compute the distances between the rows of a data matrix.
  # by default, "euclidean"
  sampleDistMatrix <- as.matrix(sampleDists)
 
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  #png("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/distance_matrix_voom_peaks.png", type="cairo",
  #    width=7,height=5,units="in",res=200,pointsize = 13)
  
  sdc= pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
  
  # dev.off()
  
  return(sdc)
}

sdc=plot_sdc(batch_corrected) 

# png("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/sdc_after_filtering_QN.png",height=9,width=9,units="in",res=300)
print(sdc)
# dev.off()

#heatmap with ggplot
create_heatmap_cor=function(x){
  cormat <- round(cor(x),2)
  head(cormat)
  
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  upper_tri <- get_upper_tri(cormat)
  upper_tri
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient(low = "white", high = "red",  
                         limit = c(0.8,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 11, hjust = 1),
          axis.title = element_blank())+
    coord_fixed()
  # Print the heatmap
  print(ggheatmap)
  
  # # add correlation in each cell
  # ggheatmap +
  #   geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  #   theme(
  #     axis.title.x = element_blank(),
  #     axis.title.y = element_blank(),
  #     panel.grid.major = element_blank(),
  #     panel.border = element_blank(),
  #     panel.background = element_blank(),
  #     axis.ticks = element_blank(),
  #     legend.justification = c(1, 0),
  #     legend.position = c(0.6, 0.7),
  #     legend.direction = "horizontal")+
  #   guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
  #                                title.position = "top", title.hjust = 0.5))

}

# heatmap is made with pearson correlation from m values, and doesn't provide information (as samples are too highly correlated). 

heatcor = create_heatmap_cor(m)
ggsave(
  paste(
    "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/Pearson_cor_Mvalues",
    currentDate,
    ".jpg",
    sep = ""
  ),
  heatcor,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300
)

## only islets
create_heatmap_cor2=function(x){
  cormat <- round(cor(x),2)
  head(cormat)
  
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  upper_tri <- get_upper_tri(cormat)
  upper_tri
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient(low = "white", high = "red",  
                        limit = c(0.8,1), space = "Lab", 
                        name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 11, hjust = 1),
          axis.title = element_blank())+
    coord_fixed()

  
  # add correlation in each cell
  ggheatmap = ggheatmap +
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 4)
  # Print the heatmap
  print(ggheatmap)
}

heatcor2 = create_heatmap_cor2(m[,c(22:ncol(m))])
ggsave(
  paste(
    "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/Pearson_cor_Mvalues_isletsonly",
    currentDate,
    ".jpg",
    sep = ""
  ),
  heatcor2,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300
)


### only diff stages
heatcor3= create_heatmap_cor(m[,c(1:21)])
ggsave(
  paste(
    "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/Pearson_cor_Mvalues_diffsonly",
    currentDate,
    ".jpg",
    sep = ""
  ),
  heatcor3,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300
)

#summary information of probes:
# do the same with piecharts and %

probes=rownames(beta)  # get probe names 

data("probe.features.epic") #load probe features from EPIC data

  cgi.info <- table(probe.features[probes,"cgi"])
  chromsome.info <- table(probe.features[probes,"CHR"])
  feature.info <- table(probe.features[probes,"feature"])
  type.info <- table(probe.features[probes,"Type"])
  
  ################# make plots prettier
  mycols <- c("#f796a0","#102542", "#feb95f", "#4ce0d2")
  
  labels = as.data.frame(cgi.info)
  labels = percent(labels$Freq/sum(labels$Freq))
  
  
  png(
    paste(outFolder, "probes_in_genome_CpGisland", currentDate, ".png", sep = ""),
    width = 6,
    height = 7,
    units = "in",
    res = 400,
    type = "cairo"
  )
  
  pie(cgi.info,
      labels = labels,
      col = mycols)
  legend(
    -.80,
    -.80,
    names(cgi.info),
    cex = 0.60,
    fill = mycols,
    title = "CpG annotations"
  )
  dev.off()
  

  
  labels = as.data.frame(feature.info)
  labels = percent(labels$Freq/sum(labels$Freq))
  features = c(
    "First exon",
    "3'UTR ",
    "5'UTR",
    "Gene body",
    "Exon boundary", # ExonBnd, within 20 basesof an exon boundary, i.e., the start or end of an exon;
    "Intergenic",
    "TSS1500", # TSS1500,200-1,500 bases upstream of the TSS.
    "TSS200"   # 0-200 bases upstream of the transcriptional start site (TSS); 
  )
  
  png(
    paste(outFolder, "genomic_annotations", currentDate, ".png", sep = ""),
    width = 6,
    height = 8,
    units = "in",
    res = 400,
    type = "cairo"
  )
  
  pie(feature.info,
      labels = labels,
      col = c("#9ee493","#baa898","#f0a868","#0353a4","#1be7ff","#00120b","#9ca3db","#eb9fef")) 
  legend(
    -.80,
    -.80,
    features,
    cex = 0.60,
    fill =  c("#9ee493","#baa898","#f0a868","#0353a4","#1be7ff","#00120b","#9ca3db","#eb9fef"),
    title = "Genomic annotations"
  )
  dev.off()
  
  
  
  labels = as.data.frame(type.info)
  labels = percent(labels$Freq/sum(labels$Freq))
  
  png(
    paste(outFolder, "probe_types", currentDate, ".png", sep = ""),
    width = 6,
    height = 6,
    units = "in",
    res = 400,
    type = "cairo"
  )
  pie(type.info,
      labels = labels,
      col = c("#102542", "#f796a0"))
  legend(
    -.80,
    -.80,
    c("I", "II"),
    cex = 0.60,
    fill = c("#102542", "#f796a0"),
    title = "Infinium probes type"
  )
  
  dev.off()
  