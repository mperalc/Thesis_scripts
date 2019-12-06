library(WGCNA)
library(reshape2)
library(dplyr)
library(ggplot2)

outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/WGCNA/ATAC_WGCNA_P12_S120_deepSplit2_merge0.30_signedHybrid_rmOutliersT_aug2019_1CPM"

# fgwas_enriched_modules = snakemake@input[["fgwas_enriched_modules"]] 

# load data

# fgwas_enriched_modules = read.table(paste(outFolder, fgwas_enriched_modules, sep = "/"), header = T)
# fgwas_enriched_modules = gsub("_ln", "", fgwas_enriched_modules$parameter)
# fgwas_enriched_modules = setdiff(fgwas_enriched_modules, "grey")

vstMat = read.table(paste(outFolder, "/vstMat.txt", sep = ""))
degData = t(vstMat)

load(paste(outFolder, "/net.xz", sep = ""))
moduleColors = labels2colors(net$colors)
datME = moduleEigengenes(degData, moduleColors)$eigengenes

#plot

stages = factor(c("iPSC","DE","GT","PF","PE","EP","EN","BLC"))

#tidy module names
colnames(datME) = gsub("ME", "", colnames(datME))
#datME = datME[,fgwas_enriched_modules]
datME = cbind(rownames(degData),datME)
colnames(datME)[1] = "stageXsample" 

#add stage info for plot
datME$stage = gsub("\\..*","",datME$stageXsample)
#melted <- melt(datME,id.vars = "stage", measure.vars = fgwas_enriched_modules)
melted <- melt(datME,id.vars = "stage", measure.vars = unique(labels2colors(net$colors)))
colnames(melted) = c("stage","module","value")

# calculate stats  
melted_sum <- melted %>%
  group_by(stage, module) %>%
  summarize(mean = mean(value),
            min   = min(value),
            max = max(value)) %>%
  ungroup()

#order stages for plot
melted_sum$stage = factor(melted_sum$stage,
                          levels = c("iPSC","DE","GT","PF","PE","EP","EN","BLC"),
                          ordered = T)  

png(paste(outFolder, "/Modules.ribbon_plots_all.png", sep = ""), width = 12,height = 9,units = "in",type = "cairo",
    res = 400,pointsize = 12)

p = ggplot(melted_sum, aes(x = stage, y = mean, group = module)) +
   geom_ribbon(aes(ymin = min, ymax = max),
                   alpha = .25, colour = NA) +
  geom_line(lwd = 1.5) +
  theme_bw() + 
  ylab("module eigengene") +
  ggtitle("Eigengenes for WGCNA ATAC-seq modules") +
    
    facet_wrap( ~ module, scales = "free", ncol = 4)

  
  
  p
  
  dev.off()
  
# for enriched modules, all together
  png(paste(outFolder, "/Modules.ribbon_plots_enriched.png", sep = ""), width = 5,height = 4,units = "in",type = "cairo",
      res = 400,pointsize = 12)
  
   p = ggplot(melted_sum[melted_sum$module %in% c("brown","green","yellow","blue"),], aes(x = stage, y = mean, 
                                                                                          group = module, col=module, fill=module)) +
    geom_ribbon(aes(ymin = min, ymax = max),
                alpha = .25, colour = NA) +
    scale_color_manual(values = c("blue","brown","yellow","green"), aesthetics = c("color","fill"))+
    geom_line(lwd = 1.5) +
    theme_bw() + 
    ylab("module eigengene") +
    ggtitle("Eigengenes for T2D-enriched WGCNA ATAC-seq modules") 
  
  p
  dev.off()
  # re-coded for paper
 # melted_sum$module = as.character(melted_sum$module)
 # melted_sum[melted_sum$module == "black","module"] = "10"
 # melted_sum[melted_sum$module == "turquoise","module"] = "1"
 # melted_sum[melted_sum$module == "green","module"] = "3"
 # melted_sum[melted_sum$module == "red","module"] = "2"
 # melted_sum$module = factor(melted_sum$module, levels = c("1","2","3","10"), ordered = T)
  
 # pdf(paste(outFolder, "/Modules.ribbon_plots_all_significant_modules_fgwas_recoded.pdf", sep = ""))
  
  #p = ggplot(melted_sum, aes(x = stage, y = mean, group = module, colour = module, fill = module)) +
   # geom_ribbon(aes(ymin = min, ymax = max),
    #            alpha = .25, colour = NA) +
    #geom_line(lwd = 1.5) +
    #scale_color_manual(values = setNames(c("#393e41","#508991","#317b22","#92140c"),
    #                                     unique(melted_sum$module))) +
    #scale_fill_manual(values = setNames(c("#393e41","#508991","#317b22","#92140c"),
    #                                    unique(melted_sum$module))) +
    #theme_bw() + 
    #theme(
    #  panel.grid.minor = element_blank(),
    #  panel.grid.major = element_blank(),
    #  panel.border = element_rect(size = 2),
    #  axis.text = element_text(size = 12, face = "bold"),
    #  axis.title = element_text(size = 14, face = "bold"),
    #  plot.title = element_text(size = 15, face = "bold"),
    #  legend.text = element_text(size = 11, face = "bold"),
    #  legend.title = element_text(size = 13, face = "bold")) +
    #ylab("module eigengene") +
    #ggtitle("Eigengenes from fGWAS-enriched WGCNA ATAC-seq modules")
  
  
  
#  print(p)
  
 # dev.off()
