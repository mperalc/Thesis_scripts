## plot DMRs

### DMP, DMR, LMR, HMR to fgwas input ############
library(GenomicRanges)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(Gviz)
inFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/"
outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/"

percent = 150 # % extra space to plot
# indicate which genome is being used
gen <- "hg19"
# the index of the DMR that we will plot. Check index in DMR annotated csv files
dmrIndex <- 32286

is_hyper = FALSE # if false, plot hypomethylated DMR
###############################################################################################


# load matrix with sample information
pD=read.csv(paste(inFolder,"samples_info.csv",sep=""), header=T)

# load matrix of normalized beta values
beta= read.csv(paste(inFolder,"quantile_normalised_beta_detP_0.01_nocrossreact.csv",sep=""), header=T,row.names = 1,check.names=F)

# Epic annotation data
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)


hyper = read.csv(paste(outFolder,"DMR/DMRs_hypermethylated_annotated_allstages.csv",sep=""))
hypo = read.csv(paste(outFolder,"DMR/DMRs_hypomethylated_annotated_allstages.csv",sep=""))




# ### Plots of specific DMRs
# set up the grouping variables and colours
pal <-  c("#0b0014", "#7883BA", "#755A91", "#CC85B1", "#C15858", "#F4B8B0", "#96665A", "#6DA567","#d81159") 

if(is_hyper == TRUE) {
# extract chromosome number and location from DMR results 
chrom <- as.character(hyper$seqnames[dmrIndex])
start <- hyper$start[dmrIndex]
end <- hyper$end[dmrIndex]
# add % extra space to plot
minbase <- start - ((percent/100)*(end-start))
maxbase <- end + ((percent/100)*(end-start))

iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name=chrom)
gTrack <- GenomeAxisTrack(col="black", cex=1, 
                          name="", 
                          fontcolor="black")
rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq", 
                    from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                    rstarts="exonStarts", rends="exonEnds", gene="name", 
                    symbol="name2", transcript="name", strand="strand", 
                    fill="darkblue",stacking="squish", name="RefSeq", 
                    showId=TRUE, geneSymbol=TRUE)

## order epic annotation data
annEPIC <- annEPIC[order(annEPIC$chr,annEPIC$pos),]

beta <- beta[match(annEPIC$Name,rownames(beta)),]
beta = na.omit(beta)
annEPIC <- annEPIC[annEPIC$Name %in% rownames(beta),]

hyper = GRanges(hyper)
# create genomic ranges object from methylation data
cpgData <- GRanges(seqnames=Rle(annEPIC$chr),
                   ranges=IRanges(start=annEPIC$pos, end=annEPIC$pos),
                   strand=Rle(rep("*",nrow(annEPIC))),
                   betas=beta)
# extract data on CpGs in DMR
cpgData <- subsetByOverlaps(cpgData,hyper[dmrIndex])

# order for plot
pD$stage = factor(pD$stage,levels = unique(pD$stage), ordered = T)
# methylation data track
methTrack <- DataTrack(range=cpgData, groups=pD$stage,genome = gen,
                       chromosome=chrom, ylim=c(-0.05,1.05), col=pal,
                       type=c("a","p"), name="DNA Methylation\n(beta value)",
                       background.panel="white", legend=TRUE, cex.title=0.8,
                       cex.axis=0.8, cex.legend=0.9)



# DMR position data track
dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                            chromosome=chrom,fill="darkred")

tracks <- list(iTrack, gTrack,rTrack,  dmrTrack,methTrack)
sizes <- c(1,1,1,1,5) # set up the relative sizes of the tracks

hyper = as.data.frame(hyper)
png(
  paste(outFolder, hyper[dmrIndex,"SYMBOL"], "_DMR_plot.png", sep = ""),
  width = 6,
  height = 7,
  units = "in",
  res = 400,
  type = "cairo"
)
p = plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes=sizes, length(tracks))
dev.off()

} else {

  # extract chromosome number and location from DMR results
  chrom <- as.character(hypo$seqnames[dmrIndex])
  start <- hypo$start[dmrIndex]
  end <- hypo$end[dmrIndex]
  # add % extra space to plot
  minbase <- start - ((percent/100)*(end-start))
  maxbase <- end + ((percent/100)*(end-start))

  iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name=chrom)
  gTrack <- GenomeAxisTrack(col="black", cex=1,
                            name="",
                            fontcolor="black")
  rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq",
                      from=minbase, to=maxbase, trackType="GeneRegionTrack",
                      rstarts="exonStarts", rends="exonEnds", gene="name",
                      symbol="name2", transcript="name", strand="strand",
                      fill="darkblue",stacking="pack", name="RefSeq",
                      showId=TRUE, geneSymbol=TRUE)

  ## order epic annotation data
  annEPIC <- annEPIC[order(annEPIC$chr,annEPIC$pos),]

  beta <- beta[match(annEPIC$Name,rownames(beta)),]
  beta = na.omit(beta)
  annEPIC <- annEPIC[annEPIC$Name %in% rownames(beta),]

  hypo = GRanges(hypo)
  # create genomic ranges object from methylation data
  cpgData <- GRanges(seqnames=Rle(annEPIC$chr),
                     ranges=IRanges(start=annEPIC$pos, end=annEPIC$pos),
                     strand=Rle(rep("*",nrow(annEPIC))),
                     betas=beta)
  # extract data on CpGs in DMR
  cpgData <- subsetByOverlaps(cpgData,hypo[dmrIndex])

  # order for plot
  pD$stage = factor(pD$stage,levels = unique(pD$stage), ordered = T)
  # methylation data track
  methTrack <- DataTrack(range=cpgData, groups=pD$stage,genome = gen,
                         chromosome=chrom, ylim=c(-0.05,1.05), col=pal,
                         type=c("a","p"), name="DNA Methylation\n(beta value)",
                         background.panel="white", legend=TRUE, cex.title=0.8,
                         cex.axis=0.8, cex.legend=0.9)



  # DMR position data track
  dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR",
                              chromosome=chrom,fill="darkred")

  tracks <- list(iTrack, gTrack,rTrack,  dmrTrack,methTrack)
  sizes <- c(1,1,1,1,5) # set up the relative sizes of the tracks

  hypo = as.data.frame(hypo)
  
  png(
    paste(outFolder, hypo[dmrIndex,"SYMBOL"], "_DMR_plot.png", sep = ""),
    width = 6,
    height = 7,
    units = "in",
    res = 400,
    type = "cairo"
  )
  p = plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE,
             add35=TRUE, grid=TRUE, lty.grid=3, sizes=sizes, length(tracks))

  dev.off()


}