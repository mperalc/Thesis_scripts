# differential methylation analysis based on limma (DMP) and DMRcate (DMR), following methylationArrayAnalysis vignette

# load packages required for analysis
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
data(dmrcatedata)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(ChIPseeker)
library(ReactomePA)
library(clusterProfiler)
library(DOSE)
library(summarytools)


### initial parameters
inFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/"
outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/"
coreNum = 8 # Number of cores to use

###


# load matrix of normalized beta values
beta= read.csv(paste(inFolder,"quantile_normalised_beta_detP_0.01_nocrossreact.csv",sep=""), header=T,row.names = 1,check.names=F)
 # load matrix of normalized M values
m = read.csv(paste(inFolder,"quantile_normalised_M_detP_0.01_nocrossreact.csv",sep=""), header=T,row.names = 1,check.names=F)

# load matrix with sample information
pD=read.csv(paste(inFolder,"samples_info.csv",sep=""), header=T)

# get the latest EPIC annotation data
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)


stages=c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

islets=colnames(beta)[22:32]  # islets names

### differentially methylated probes (DMP)  ############

# create the design matrix, which contains the comparisons I want to make


# this is the factor of interest
stage <- pD$stage
# this is the individual effect that we need to account for
sample <- as.character(pD$sample)
## change R for ISL because otherwise  eBayes fails
sample[sample=="R"] = "ISL"
sample = factor(sample,levels = unique(sample))
pD_modified = pD[c("stage","sample")]
pD_modified[pD_modified$sample=="R","sample"] = "ISL"
pD_modified$sample = factor(pD_modified$sample,levels = unique(pD_modified$sample))

# use the above to create a design matrix
design <- model.matrix(~0+stage+sample, data=pD_modified)
colnames(design) <- c(levels(stage),levels(sample)[-1])

# fit the linear model 
fit <- lmFit(m, design)
# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(DE-iPSC,
                            GT-iPSC,
                            PF-iPSC,
                            PE-iPSC,
                            EP-iPSC,
                            EN-iPSC,
                            BLC-iPSC,
                            islet-BLC,
                            levels=design)
contMatrix


# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))


## results for every contrast
DMPs = list()
for(n in 1:8) {
  DMPs[[n]] =  topTable(
      fit2,
      number = Inf,
      coef = n,
      sort.by = "p"
    )
  DMPs[[n]] = DMPs[[n]][abs(DMPs[[n]]$logFC)>1 & DMPs[[n]]$adj.P.Val<0.01,]
 
  DMPs[[n]]$CpG = rownames(DMPs[[n]]) 
  DMPs[[n]]$contrast = rep(n,nrow(DMPs[[n]])) 
  
}

names(DMPs) = colnames(contMatrix)

for (name in names(DMPs)){
  # distinguish between increase or decrease in methylation
# write.csv(DMPs[[name]][DMPs[[name]]$logFC>0,], file=paste(outFolder,name,"_DMPs_hypermethylated.csv",sep=""), row.names=FALSE)
# write.csv(DMPs[[name]][DMPs[[name]]$logFC<0,], file=paste(outFolder,name,"_DMPs_hypomethylated.csv",sep=""), row.names=FALSE)
}

## sumamrize
DMP_summary = do.call("rbind",DMPs)
summarytools::freq(DMP_summary$contrast)
summarytools::freq(DMP_summary[DMP_summary$logFC>0,"contrast"])
summarytools::freq(DMP_summary[DMP_summary$logFC<0,"contrast"])

### Gene ontology
for(n in names(DMPs)){
sigCpGs_hyper = DMPs[[n]][DMPs[[n]]$logFC>0,"CpG"]
sigCpGs_hypo = DMPs[[n]][DMPs[[n]]$logFC<0,"CpG"]

# Get all the CpG sites used in the analysis to form the background
all <- rownames(m)

# The genes that have more CpGs associated with them will have a higher probability of being identified as differentially methylated 
# compared to genes with fewer CpGs. We can look at this bias in the data by specifying plot=TRUE

gst_hyper <- gometh(sig.cpg=sigCpGs_hyper, all.cpg=all, plot.bias=TRUE)
gst_hypo <- gometh(sig.cpg=sigCpGs_hypo, all.cpg=all, plot.bias=TRUE)

write.csv(gst_hyper, file=paste(outFolder,n,"_DMPs_GO_hypermethylated.csv",sep=""), row.names=TRUE)
write.csv(gst_hypo, file=paste(outFolder,n,"_DMPs_GO_hypomethylated.csv",sep=""), row.names=TRUE)

}
####

### annotate features
 annEPICSub <- annEPIC[match(rownames(m),annEPIC$Name),
                       c(1:4,12:19,24:ncol(annEPIC))]
# rm(annEPIC) ## free up memory
# annEPICSub = annEPICSub[!(is.na(rownames(annEPICSub))),] ## remove missing values
# annEPICSub = annEPICSub[c(1:5,11:13,21)]
# DMPs_anno[[n]] = merge(DMPs,
#                   as.data.frame(annEPICSub),
#                   by = 0,
#                   all.x = TRUE)




##### DMRs ####

# matrix of M-values is annotated with the relevant information about the probes such as their genomic position, 
# gene annotation, etc. By default, this is done using the ilmn12.hg19 annotation, but this can be substituted for 
# any argument compatible with the interface provided by the minfi package. 
# The limma pipeline is then used for differential methylation analysis to calculate moderated t-statistics.
# (DMRcate is bassed on limma)

results.ranges = list()
results_df_list = list()
for(c in colnames(contMatrix)) {
  

myAnnotation <- cpg.annotate(object = as.matrix(m), datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = c, arraytype = "EPIC",fdr = 0.01  )
# this step calculates stats for all CpGs (800840)
# so I will use these as background for enrichment

# Once we have the relevant statistics for the individual CpGs, we can then use the dmrcate function to combine them to 
# identify differentially methylated regions. 
# The main output table DMRs$results contains all of the regions found, along with their genomic annotations and p-values.

DMRs <- dmrcate(myAnnotation, lambda=1000, 
                C=2)

# As for the probe-wise analysis, it is advisable to visualise the results to ensure that they make sense.
# The regions can easily be viewed using the DMR.plot function provided in the DMRcate package 

results.ranges[[c]] <- extractRanges(DMRs, genome = "hg19")
results_df_list[[c]] = as.data.frame(results.ranges[[c]])
results_df_list[[c]]$contrast = rep(c, nrow(results_df_list[[c]]))
# write.csv(results_df_list[[c]][results_df_list[[c]]$meanbetafc>0,], file=paste(outFolder,c,"_DMRs_hypermethylated.csv",sep=""), row.names=FALSE)
# write.csv(results_df_list[[c]][results_df_list[[c]]$meanbetafc<0,], file=paste(outFolder,c,"_DMRs_hypomethylated.csv",sep=""), row.names=FALSE)

}


## summarize results_df_list[[c]] in terms of hyper and hypomethylated DMRs per stage and make table
results_df_list = do.call("rbind",results_df_list)
results_df_list_hypo = results_df_list[results_df_list$meanbetafc<0,]
results_df_list_hyper = results_df_list[results_df_list$meanbetafc>0,]

freq(results_df_list$contrast,cumul = FALSE) 
freq(results_df_list_hypo$contrast,cumul = FALSE)  ### DMR fc<0 per contrast
freq(results_df_list_hyper$contrast, cumul = FALSE)  ### DMR fc>0 per contrast

summarytools::descr(results_df_list$minfdr)

### reactome and GO analysis


results_df_list_annotated = annotatePeak(makeGRangesFromDataFrame(results_df_list), 
                                              tssRegion=c(-2000, 500),TxDb=txdb, annoDb="org.Hs.eg.db")
results_df_list_hypo_annotated = annotatePeak(makeGRangesFromDataFrame(results_df_list_hypo), 
                                         tssRegion=c(-2000, 500),TxDb=txdb, annoDb="org.Hs.eg.db")
results_df_list_hyper_annotated = annotatePeak(makeGRangesFromDataFrame(results_df_list_hyper), 
                                         tssRegion=c(-2000, 500),TxDb=txdb, annoDb="org.Hs.eg.db")
# most DMRs  are located at promoters

png(
  paste(outFolder, "annotation_pie_DMRs_hyper.png", sep = ""),
  width = 7,
  height = 7,
  units = "in",
  res = 400,
  type = "cairo"
)
plotAnnoPie(results_df_list_hyper_annotated)
dev.off()

## save annotated 
results_df_list_hypo_annotated = as.data.frame(results_df_list_hypo_annotated)
results_df_list_hyper_annotated = as.data.frame(results_df_list_hyper_annotated)
results_df_list_hypo_annotated = cbind(results_df_list_hypo_annotated, results_df_list_hypo[c(6:12)])
results_df_list_hyper_annotated = cbind(results_df_list_hyper_annotated, results_df_list_hyper[c(6:12)])

write.csv(results_df_list_hyper_annotated, file=paste(outFolder,"DMRs_hypermethylated_annotated_allstages.csv",sep=""), row.names=FALSE)
write.csv(results_df_list_hypo_annotated, file=paste(outFolder,"DMRs_hypomethylated_annotated_allstages.csv",sep=""), row.names=FALSE)


## make universe for pathway enrichment background
universe = data.frame(cbind(myAnnotation$CHR,myAnnotation$pos,(myAnnotation$pos)+1))
colnames(universe) = c("chr","start","end")
universe$chr = paste0("chr",universe$chr)
universe = annotatePeak(GRanges(universe),tssRegion=c(-2000, 500),TxDb=txdb, annoDb="org.Hs.eg.db") # just to put GRanges in right format


universe_annotated = seq2gene(
  universe@anno,
  tssRegion = c(-2000, 500),
  flankDistance = 5000,
  TxDb = txdb
)


for(n in names(results.ranges)) {
  
  results.ranges[[n]] = GRanges(results_df_list_hypo[results_df_list_hyper$contrast == n,])
  # test hypomethylated separately: no sign result
  # test hypermethylated: 
  gene <-
    seq2gene(
      results.ranges[[n]],
      tssRegion = c(-2000, 500),
      flankDistance = 5000,
      TxDb = txdb
    )
  pathway <- enrichPathway(gene,organism = "human",minGSSize = 5,readable = TRUE,
                           universe = universe_annotated)
  
  
  png(
    paste(outFolder, n, "_pathway_enrichplot.png", sep = ""),
    width = 6,
    height = 7,
    units = "in",
    res = 400,
    type = "cairo"
  )
  plot(dotplot(pathway))
  
  dev.off()

# enrich pathway without cqll to seq2gene gives less significant terms

##GO

png(
  paste(outFolder, n, "_GO_enrichplot.png", sep = ""),
  width = 6,
  height = 7,
  units = "in",
  res = 400,
  type = "cairo"
)
ego <- enrichGO(gene          = gene,
                universe      = universe_annotated,
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
# nothing found
plot(dotplot(ego))

dev.off()
}
