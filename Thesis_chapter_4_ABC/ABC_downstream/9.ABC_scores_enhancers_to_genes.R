### Heatmap of ABC scores of selected genes

library(reshape2)
library(tidyverse)
library(scales)
library(GenomicRanges)
library(ComplexHeatmap)
library(tidyr)
library(ggpubr)

# from targets of enhancers overlapping T2D genes
predicted_targets = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/functional_PPA/PPA_genetic_functional_ABC_fully_annotated.txt",
                               header = T)
topPPA = 0.5
# select only top PPAf
predicted_targets = predicted_targets[predicted_targets$PPAf>topPPA,]

scores_all = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/all_expressed.txt",
                        header = T)

### for selected enhancers

 predicted_targets = predicted_targets[, c("Chr" ,
                                           "Start",
                                           "End",
                                           "ID",
                                           "IndexSNP",
                                           "credset",
                                           "target_Gene",
                                           "PPAf")]

 
 
 predicted_targets.gr = GRanges(predicted_targets)
 scores_all.gr = GRanges(scores_all)
 
 overlaps = findOverlaps(scores_all.gr, predicted_targets.gr)
 
 predicted_targets =  predicted_targets[overlaps@to,]
 subset = scores_all[overlaps@from, ]
 
 subset = as.data.frame(subset)

 
 nrow(predicted_targets) == nrow(subset) 
 subset = cbind(subset[c("name",
                         "cellType",
                         "TargetGene",
                         "TargetGeneExpression",
                         "ABC.Score")], 
                predicted_targets[c("ID", "IndexSNP", "credset","PPAf")])
 
 subset$cellType = factor(subset$cellType,levels = c("iPSC","DE","GT","PF","PE","EP","EN","BLC"), ordered = T)
 subset$group = paste(paste(subset$credset,subset$ID , sep = " - " ),subset$TargetGene,sep = " --> ")
 subset$labels = paste(subset$credset,subset$ID , sep = " - " )
 
 subset = subset[order(subset$PPAf,decreasing = T),]
 subset$group = factor(subset$group,levels = rev(unique(subset$group)),ordered = T)
 
p1= ggplot(data = subset,aes(x=cellType, y = group)) +
   geom_tile(aes(fill = ABC.Score),colour = "white") +
   scale_fill_gradient(low = "mistyrose1", high = "firebrick4") +
   theme_bw() +
  xlab("Stage")+
  ylab("Credible set variants") 
  
 p1
 ## same for gene expression in TPM - mean-corrected
 
 subset2 = subset
 subset2$TargetGeneExpression = log2( subset2$TargetGeneExpression )
 #scaled for each gene by mean
 means <- subset2 %>%
   group_by(TargetGene) %>%
   summarise(TargetGeneExpression = mean(TargetGeneExpression))

 means = as.data.frame(means)
 subset2$log2_scaled =  subset2$TargetGeneExpression - means[match(subset2$TargetGene, means$TargetGene),"TargetGeneExpression"]
 subset2$TargetGene = factor(subset2$TargetGene,levels = rev(unique(subset2$TargetGene)),ordered = T)
 
 
 
 p2 = ggplot(data = subset2,aes(x=cellType,y=TargetGene)) +
   geom_tile(aes(fill = log2_scaled),colour = "white") +
   scale_fill_gradientn(colours=c("lightblue2","white","coral2")) 
p2

 ggarrange(p1, p2, labels = c("A", "B"),
           common.legend = F)
 
 # plot gene expression throughout all stages
 directory = "/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/"
 filename = "31.01.2017.Differentiation_v2.gene.tpm.tsv"
 tpm = as.data.frame(read_tsv(
   paste(
     directory,
     filename,
     sep = ""
   )
 ))   # read tpm file for plotting longitudinal tpm data
 
 # Data from input columns
 donor = c("Ad2.1", "Ad3.1", "Neo1.1")  # original samples, here called donor
 
 stage = c("iPSC", "DE", "PGT", "PFG", "PE", "EP", "ENstage6", "ENstage7") # 8 differentiation stages
 
 # alternative names to plot:
 stage_2 = c("iPSC", "DE", "GT", "PF", "PE", "EP", "EN", "BLC")  #shortening EN names EN7= BLC (beta-like cells)
 genes = subset$TargetGene
 
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
   
   gene_number = nrow(unique(plot_long)  ) # how many genes to plot
   
   # melt data for ggplot2
   
   long = melt(unique(plot_long), measure.vars = c(3:ncol(plot_long)))
   head(long)
   
   # rename stages and samples
   
   samples <- c(rep(donor, 8))
   
   long$variable = rep(stage_2, each = 3 * gene_number)                    # sample size times number of genes
   
   colnames(long)[which(names(long) == "variable")] <- "stage"
   long$Sample = rep(samples, each = gene_number)
   long$stage <- factor(long$stage, levels = stage_2)
   long$Sample = as.factor(long$Sample)
   long$GeneName = factor(long$GeneName, levels = unique(genes))
   
   means <- long %>%
     group_by(GeneName,stage) %>%
     summarise(value = mean(value))
   
 means = as.data.frame(means)
 
 means$log2TPM= log2( means$value )
 #scaled for each gene by mean
 means_of_means <-  means %>%
   group_by(GeneName) %>%
   summarise(log2TP = mean(log2TPM))
 
 means_of_means = as.data.frame(means_of_means)
 means$log2_scaled =  means$log2TP - means[match(means$GeneName, means_of_means$GeneName),"log2TPM"]
 means$GeneName = factor( means$GeneName,levels = rev(unique(subset2$TargetGene)),ordered = T)
 means$log2TPM[means$value == 0] = "NA"
 means$log2TPM = as.numeric(means$log2TPM)
 means = na.omit(means)
 
 p3 = ggplot(data = means,aes(x=stage,y=GeneName)) +
   geom_tile(aes(fill = log2_scaled),colour = "white") +
   scale_fill_gradientn(colours=c("lightblue2","white","coral2")) +
   theme_bw()
 p3
 
 ggarrange(p1, p3, labels = c("A", "B"),
           common.legend = F)
 
 # not scaled
 p3 = ggplot(data = means,aes(x=stage,y=GeneName)) +
   geom_tile(aes(fill = log2TPM),colour = "white") +
   scale_fill_gradientn(colours=c("lightblue2","white","coral2"),breaks = c(min(means$log2TPM),1,max(means$log2TPM)),
                        labels  = c(round(min(means$log2TPM),0),0,round(max(means$log2TPM),0)),
                                    limits =c(min(means$log2TPM),max(means$log2TPM)) ) +
   xlab("Stage")+
   ylab("Target Genes")
 p3
 
 p = ggarrange(p1, p3, labels = c("A", "B"),
           common.legend = F,widths = c(1.4,1))
 
 png(paste0("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/functional_PPA/PPAf_top_",topPPA,".png"),
      height = 10,width = 10,units = "in",res = 400, type = "cairo")
 p
 dev.off()
 # manual list
 # genes = c("POU5F1","NANOG","GATA4","GATA6","FOXA2","SOX9")
 
 # 
 # # subset scores from selected genes for all enhancers
 # scores_all_subset = scores_all[scores_all$TargetGene %in% genes,c("name","cellType","TargetGene","ABC.Score")]
 # scores_all_subset$group = paste(scores_all_subset$name,scores_all_subset$TargetGene,sep = "--")
 # scores_all_subset$cellType = factor(scores_all_subset$cellType,levels = c("iPSC","DE","GT","PF","PE","EP","EN","BLC"), ordered = T)
 # 
 # # plot lineplot
 # # ggplot(data = scores_all_subset,aes(x=cellType,y=ABC.Score,col=TargetGene, group = group)) +
 # #   geom_line() +
 # #   geom_point()
 # 
 # 
 # ## majority of enhancers are unique
 # # better to calculate average
 # means <- scores_all_subset %>% 
 #   group_by(TargetGene,cellType) %>% 
 #   summarise(ABC.Score = mean(ABC.Score))
 # 
 # # plot heatmap
 # ggplot(data = means,aes(x=cellType,y=TargetGene)) +
 #   geom_tile(aes(fill = ABC.Score),colour = "white") + 
 #   scale_fill_gradient(low = "white", high = "firebrick1")
 # 
 # # scaled for each gene by mean
 # means_of_means <- means %>% 
 #   group_by(TargetGene) %>% 
 #   summarise(ABC.Score = mean(ABC.Score))
 # 
 # means = as.data.frame(means)
 # means_of_means = as.data.frame(means_of_means)
 # 
 # means$ABC.Score_normalized = means$ABC.Score - means_of_means[match(means$TargetGene, means_of_means$TargetGene),"ABC.Score"]
 # 
 # ggplot(data = means,aes(x=cellType,y=TargetGene)) +
 #   geom_tile(aes(fill = ABC.Score_normalized),colour = "white") + 
 #   scale_fill_gradientn(colours=c("steelblue","white","firebrick1"),
 #                        values  = rescale(c(min(means$ABC.Score_normalized), 0, max(means$ABC.Score_normalized))))
 