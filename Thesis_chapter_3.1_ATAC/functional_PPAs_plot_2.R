# Functional PPAs from ATAC-seq peaks (islet diff model)
# Annotations to calculate functional PPAs come from WGCNA modules enriched in T2D (full model)
# 1CPM significant: brown,green , yellow , blue
library(ggplot2)
library(ggrepel)

genetic_credset_PPA_folder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/functional_PPA_fGWAS_WGCNA/credset_T2D_2018"
functional_PPA_folder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/functional_PPA_fGWAS_WGCNA/1_CPM/"
variant_list = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/functional_PPA_fGWAS_WGCNA/variant_list.txt")
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
files = base::setdiff(files,"all_modules_selecting_3")
files = files[grep(".bfs.gz",files,fixed = T)]
for (file in files) {
  functional_PPA[[file]] = read.table(paste(functional_PPA_folder,file,sep = "/"), header = T)
  functional_PPA[[file]]$id=gsub("\\:", "_", functional_PPA[[file]]$id)
  functional_PPA[[file]]$id=gsub("chr", "", functional_PPA[[file]]$id)
  functional_PPA[[file]]$file = file
}

test_functional_PPA = do.call("rbind",functional_PPA)

PPA_genetic_functional = merge(genetic_PPA[c("IndexSNP","PPAg","id","credset")],
                               test_functional_PPA[c("file","id","PPA","brown","green","yellow","blue")])

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
# Some variants are present in two stages
# if they do not overlap any parameter, they should not be upweighted (as they are not the top variant for that conditional GWAS)

to_ditch = c()
for( v in unique(PPA_genetic_functional$id)) {
to_test = PPA_genetic_functional[PPA_genetic_functional$id == v,] ## per id
to_test = to_test[to_test$brown == 0 & to_test$green == 0 & to_test$yellow == 0 & to_test$blue == 0 ,]
if(nrow(to_test)>1){
to_test = to_test[which(to_test$PPA == max(to_test$PPA)),] ## rowname of max underserving PPA is dropped
to_ditch = c(to_ditch,rownames(to_test))
}
}

to_keep = base::setdiff(as.character(rownames(PPA_genetic_functional)),to_ditch)
PPA_genetic_functional = PPA_genetic_functional[to_keep,]
PPA_genetic_functional = PPA_genetic_functional[!duplicated(PPA_genetic_functional$id),] 

# No duplicated ids
length(PPA_genetic_functional$id) == length(unique(PPA_genetic_functional$id))

# plot

PPA_genetic_functional$module = rep("not enriched", nrow(PPA_genetic_functional))
PPA_genetic_functional$module[PPA_genetic_functional$blue == 1] = "blue"
PPA_genetic_functional$module[PPA_genetic_functional$yellow == 1] = "yellow"
PPA_genetic_functional$module[PPA_genetic_functional$green == 1] = "green"
PPA_genetic_functional$module[PPA_genetic_functional$brown == 1] = "brown"
PPA_genetic_functional$module = factor(PPA_genetic_functional$module,levels=c("brown","green","yellow","blue","not enriched"), ordered = T)


# Variable transparency:
# Add bivariate density for each point
#  PPA_genetic_functional$density_new <- fields::interp.surface(
#  MASS::kde2d(PPA_genetic_functional$PPAg, PPA_genetic_functional$PPA, n=1000), 
#    PPA_genetic_functional[,c("PPAg", "PPA")])
# # # Adjusting so that 1/density is between 0 and 1
#  PPA_genetic_functional[PPA_genetic_functional$density == 0,"density"] = 4.762731e-282
#  PPA_genetic_functional[PPA_genetic_functional$density < 1,"density"] = 1

# Add transparency only to "not enriched" module points
PPA_genetic_functional$difference = PPA_genetic_functional$PPA -PPA_genetic_functional$PPAg
PPA_genetic_functional$alpha = rep(0.3, nrow(PPA_genetic_functional))
PPA_genetic_functional[!(PPA_genetic_functional$module=="not enriched"), "alpha"] = 1

# Add labels for top points from enriched modules
PPA_genetic_functional$labels = rep("NA", nrow(PPA_genetic_functional))
PPA_genetic_functional[PPA_genetic_functional$PPA>0.70 & !(PPA_genetic_functional$module=="not enriched") ,
                       "labels"] = PPA_genetic_functional[PPA_genetic_functional$PPA>0.70 & !(PPA_genetic_functional$module=="not enriched"),
                                                          "credset"]

toplot = PPA_genetic_functional[c("PPA","PPAg","module","labels","alpha")]
# Plot
plot = ggplot(toplot, aes(x = PPAg, y = PPA,
                                          #alpha = 1/density,
                                          #label = labels,
                                          color = module)) +
                                          
  # scale_x_continuous(expand=c(0,0)) +
  # scale_y_continuous(expand=c(0,0)) +
  
   geom_point(size=5, alpha = PPA_genetic_functional$alpha) +
   scale_color_manual(values = setNames(c("#d8d6d6","blue","brown","yellow","green"),
                                        unique(PPA_genetic_functional$module))) +
   scale_fill_manual(values = setNames(c("#d8d6d6","blue","brown","yellow","green"),
                                       unique(PPA_genetic_functional$module))) +
  xlab("genetic PPA") +
  ylab("functional PPA") + 
  ggtitle("Re-calculated PPAs from ATAC-seq module enrichment information") + 
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

ggsave("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/functional_PPA_fGWAS_WGCNA/1_CPM/functional_PPA_ATAC_WGCNA_1CPM.png", plot = plot,
        width = 7, height = 7, units = "in",dpi = 400)

# Those that go over 50% PPA after running fGWAS.
PPA_genetic_functional[PPA_genetic_functional$PPA>0.50 & PPA_genetic_functional$PPAg<0.50,]

# Specific locus
PPA_genetic_functional[grep(pattern = "PROX1", x = PPA_genetic_functional$credset),]

### annotate with SNP ids

ids = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/HRC_credset.snp_ann.txt")
colnames(ids) = c("ID","SNP_pos","locus")
colnames(PPA_genetic_functional)[1] = "SNP_pos"
PPA_genetic_functional$SNP_pos = gsub("_",":",PPA_genetic_functional$SNP_pos)
PPA_genetic_functional$IndexSNP = gsub("_",":",PPA_genetic_functional$IndexSNP)

PPA_genetic_functional = merge(PPA_genetic_functional,ids, by = "SNP_pos", all.x = T)
PPA_genetic_functional$locus = sapply(strsplit(PPA_genetic_functional$credset, "_"), function(x) x[4])
PPA_genetic_functional = PPA_genetic_functional[c("locus","ID","SNP_pos"  ,  "IndexSNP" ,  "PPAg" ,  "PPA" , "module"   ,"difference"       )]
write.table(PPA_genetic_functional[c(1:ncol(PPA_genetic_functional))],"/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/functional_PPA_fGWAS_WGCNA/1_CPM/PPA_genetic_functional_ATAC_WGCNA_1CPM.txt",
            quote = F, col.names = T, row.names = F, sep = "\t")
