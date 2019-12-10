# Functional PPAs from ABC enhancers
# Annotations to calculate functional PPAs come from stages enriched in T2D (best joint model, EN and PF)
# 
library(ggplot2)
library(ggrepel)

genetic_credset_PPA_folder = "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/credset_T2D_2018"
functional_PPA_folder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/functional_PPA"
variant_list = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/re-weighted_functional_data_WGCNA_fgwas_2019/variant_list.txt")
variant_list$id = paste(variant_list$V1,variant_list$V2,sep="_")

genetic_PPA = list()
functional_PPA = list()
for (file in list.files(path = genetic_credset_PPA_folder)) {
  genetic_PPA[[file]] = read.table(paste(genetic_credset_PPA_folder,file,sep = "/"), header = T)
  genetic_PPA[[file]]$credset = rep(file,nrow(genetic_PPA[[file]])) 
}
genetic_PPA = do.call("rbind",genetic_PPA)
genetic_PPA$id = paste(genetic_PPA$Chr,genetic_PPA$Pos,sep = "_")
files = list.files(path = functional_PPA_folder)
files = files[grep(".bfs.gz",files,fixed = T)]
for (file in files) {
  functional_PPA[[file]] = read.table(paste(functional_PPA_folder,file,sep = "/"), header = T)
  functional_PPA[[file]]$id=gsub("\\:", "_", functional_PPA[[file]]$id)
  functional_PPA[[file]]$id=gsub("chr", "", functional_PPA[[file]]$id)
  functional_PPA[[file]]$file = file
}

test_functional_PPA = do.call("rbind",functional_PPA)

PPA_genetic_functional = merge(genetic_PPA[c("IndexSNP","PPAg","id","credset")],
                               test_functional_PPA[c("file","id","PPA","EN","PF")])

PPA_genetic_functional = PPA_genetic_functional[PPA_genetic_functional$PPA != "NaN",]


### test
# take each SNP from the file they should belong to
keep_all = c()
for (i in variant_list$id){
  file = variant_list[variant_list$id == i, "V3"]
  file = gsub("fgwas.input","",file)
  keep = rownames(PPA_genetic_functional[PPA_genetic_functional$IndexSNP == i & grepl(file,PPA_genetic_functional$file), ])
  keep_all = c(keep_all,keep)
}
PPA_genetic_functional = PPA_genetic_functional[keep_all,]
####




length(PPA_genetic_functional$id)
length(unique(PPA_genetic_functional$id))

# There are various genetic credible sets that share the same SNPs with different PPAs
# Functional PPAs for those don't vary
# Sort by higher PPAg
PPA_genetic_functional = PPA_genetic_functional[order(PPA_genetic_functional$PPAg,decreasing = T),]
# Some variants are present in two files
# if they do not overlap any parameter, they should not be upweighted (as they are not the top variant for that conditional GWAS)


to_ditch = c()
i = 1
for (v in unique(PPA_genetic_functional$id)) {
  to_test = PPA_genetic_functional[PPA_genetic_functional$id == v, ] ## per id
  to_test = to_test[to_test$EN == 0 & to_test$PF == 0 , ]
  if (nrow(to_test) > 1) {
    to_test = to_test[which(to_test$PPA == max(to_test$PPA)), ] ## rowname of max underserving PPA is dropped
    to_ditch = c(to_ditch, rownames(to_test))
  }
  if (floor((i / length(unique(PPA_genetic_functional$id)))) %in% seq(0, 100, by =
                                                              5)) {
    message(paste0("progress: ", (i / length(
      unique(PPA_genetic_functional$id)
    )) * 100, "%"))
  }
  i = i + 1
}

to_keep = base::setdiff(as.character(rownames(PPA_genetic_functional)),to_ditch)
PPA_genetic_functional = PPA_genetic_functional[to_keep,]
PPA_genetic_functional = PPA_genetic_functional[!duplicated(PPA_genetic_functional$id),] 



PPA_genetic_functional = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/functional_PPA/PPA_genetic_functional_ABC.txt",
                                    header = T)


# No duplicated ids
length(PPA_genetic_functional$id) == length(unique(PPA_genetic_functional$id))

# plot

PPA_genetic_functional$stage = rep("not enriched", nrow(PPA_genetic_functional))
PPA_genetic_functional$stage[PPA_genetic_functional$EN == 1] = "EN"
PPA_genetic_functional$stage[PPA_genetic_functional$PF == 1] = "PF"
PPA_genetic_functional$stage[PPA_genetic_functional$PF == 1 & PPA_genetic_functional$EN == 1 ] = "both"
PPA_genetic_functional$stage = factor(PPA_genetic_functional$stage,levels=c("PF","EN","both","not enriched"), ordered = T)

# Variable transparency:
# Add bivariate density for each point
#  PPA_genetic_functional$density_new <- fields::interp.surface(
#  MASS::kde2d(PPA_genetic_functional$PPAg, PPA_genetic_functional$PPA, n=1000), 
#    PPA_genetic_functional[,c("PPAg", "PPA")])
# # # Adjusting so that 1/density is between 0 and 1
#  PPA_genetic_functional[PPA_genetic_functional$density == 0,"density"] = 4.762731e-282
#  PPA_genetic_functional[PPA_genetic_functional$density < 1,"density"] = 1

# Add transparency only to "not enriched" stage points
PPA_genetic_functional$difference = PPA_genetic_functional$PPA -PPA_genetic_functional$PPAg
PPA_genetic_functional$alpha = rep(0.3, nrow(PPA_genetic_functional))
PPA_genetic_functional[!(PPA_genetic_functional$stage=="not enriched"), "alpha"] = 1

# Add labels for top points from enriched stages
PPA_genetic_functional$labels = rep("NA", nrow(PPA_genetic_functional))
PPA_genetic_functional[PPA_genetic_functional$PPA>0.70 & !(PPA_genetic_functional$stage=="not enriched") ,
                       "labels"] = PPA_genetic_functional[PPA_genetic_functional$PPA>0.70 & !(PPA_genetic_functional$stage=="not enriched"),
                                                          "credset"]

toplot = PPA_genetic_functional[c("PPA","PPAg","stage","labels","alpha")]
# Plot
plot = ggplot(toplot, aes(x = PPAg, y = PPA,
                                          #alpha = 1/density,
                                          #label = labels,
                                          color = stage)) +
                                          
  # scale_x_continuous(expand=c(0,0)) +
  # scale_y_continuous(expand=c(0,0)) +
  
   geom_point(size=5, alpha = PPA_genetic_functional$alpha) +
   scale_color_manual(values = setNames(c("#d8d6d6","#317b22","#23d1e9","#0c0c0c"),
                                        unique(PPA_genetic_functional$stage))) +
   scale_fill_manual(values = setNames(c("#d8d6d6","#317b22","#23d1e9","#0c0c0c"),
                                       unique(PPA_genetic_functional$stage))) +
  xlab("genetic PPA") +
  ylab("functional PPA") + 
  ggtitle("Re-calculated PPAs from ABC enhancers per enriched stage") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  # geom_text_repel() +
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"))

 plot

ggsave("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/functional_PPA/functional_PPA_ABC.png", plot = plot,
        width = 7, height = 7, units = "in",dpi = 400)

write.table(PPA_genetic_functional[c(1:(ncol(PPA_genetic_functional)-2))],"/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/functional_PPA/PPA_genetic_functional_ABC.txt",
            quote = F, col.names = T, row.names = F, sep = "\t")
# Those that go over 70% PPA after running fGWAS.
PPA_genetic_functional[PPA_genetic_functional$PPA>0.70 & PPA_genetic_functional$PPAg<0.70,]


# Specific locus
PPA_genetic_functional[grep(pattern = "PROX1", x = PPA_genetic_functional$credset),]
