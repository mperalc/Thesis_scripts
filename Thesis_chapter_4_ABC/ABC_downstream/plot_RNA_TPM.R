#
#  PLOT LONGITUDINAL DAta
#

#   Data from TPM normalised counts
#   List of genes of interest
# IN: file with longitudinal data, list of elements to plot
# OUT: plot with longitudinal data

plot_RNA_TPM = function(genes, name = "test"){
  
# Load necessary libraries
require(readr)  # to read big tables fast
require(ggplot2) # to plot
require(reshape2)  # to modify dataframes for ggplot2

###############################################################
#           VARIABLE
###############################################################

directory = "data/counts/"
filename = "31.01.2017.Differentiation_v2.gene.tpm.tsv"

# Data from input columns
donor = c("Ad2.1", "Ad3.1", "Neo1.1")  # original samples, here called donor

stage = c("iPSC", "DE", "PGT", "PFG", "PE", "EP", "ENstage6", "ENstage7") # 8 differentiation stages

# alternative names to plot:
stage_2 = c("iPSC", "DE", "GT", "PF", "PE", "EP", "EN", "BLC")  #shortening EN names EN7= BLC (beta-like cells)

# rna-seq data
tpm = as.data.frame(read_tsv(
  paste(
    directory,
    filename,
    sep = ""
  )
))   # read tpm file for plotting longitudinal tpm data

output_type = "pdf"  

# colour palette
diaPalette <-
  c("#C15858", "#6DA567", "#7883BA")  # Diabetologia palette

########################################################


# QC and re-shaping of data for plotting

QC_and_reshaping = function() {
  #genes=sort(genes) # sort in alphabetical order
  plot_long = tpm[match(genes, tpm$GeneName), ]  # extracts from tpm data frame the rows that contain the genes of interest
  plot_long = na.omit(plot_long)              # remove NAs, in case there was not a match for a gene
  
  diff = setdiff(genes, tpm$GeneName)           # which ones in the genes' list are not in the table (probably have a different name)?
  
  #order the columns by stage
  
  nc = plot_long[, grepl("iPSC" , names(plot_long))]  #takes columns whose name matches x
  
  #Do the same for the other stages
  
  for (s in stage[-1])  {
    i = plot_long[, grepl(s , names(plot_long))]
    nc = cbind(nc, i)
  }
  
  plot_long = cbind(plot_long[c(1:2)], nc)
  rm(i, nc)
  
  gene_number = nrow(plot_long)   # how many genes to plot
  
  # melt data for ggplot2
  
  long = melt(plot_long, measure.vars = c(3:ncol(plot_long)))
  head(long)
  
  # rename stages and samples
  
  samples <- c(rep(donor, 8))
  
  long$variable = rep(stage_2, each = 3 * gene_number)                    # sample size times number of genes
  
  colnames(long)[which(names(long) == "variable")] <- "stage"
  long$Sample = rep(samples, each = gene_number)
  long$stage <- factor(long$stage, levels = stage_2)
  long$Sample = as.factor(long$Sample)
  long$GeneName = factor(long$GeneName, levels = genes)
  return(long)
  
} # this function takes parameters from the global environment
# Working directory:

    long = QC_and_reshaping()  # reshape for plotting
    head(long)
    
    #create plots first
    plot_list = list()
    
    for (i in unique(long$GeneName))
    {
      long2 = long[long$GeneName == i, ] # subset each gene
      p <-
        ggplot(data = long2, aes(x = stage, y = value, group = Sample)) +
        ggtitle(unique(long2$GeneName)) +
        xlab ("Differentiation stages") +
        ylab ("Expression [TPM]") +
        expand_limits(y = 0) +
        geom_hline(
          yintercept = 0,
          linetype = "dashed",
          size = 1,
          col = "#DCDCDC"
        ) +
        geom_line(aes(linetype = Sample, col = Sample), size = 1) +
        scale_colour_manual(values = diaPalette) +  # diabetologia pallete
        geom_point(size = 3, aes(shape = Sample, col = Sample)) +
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
      
      plot_list[[i]] = p
      print(i)
    }
  
    
    # I have to open pdf connection before looping so that it saves one gene on each page
    somePDFPath = paste("data/plots/",
                        name,
                        ".pdf",
                        sep = "")
    pdf(file = somePDFPath)
    for (i in unique(long$GeneName)) {
      print(plot_list[[i]])
    }
    dev.off()
    
  
    }

  

