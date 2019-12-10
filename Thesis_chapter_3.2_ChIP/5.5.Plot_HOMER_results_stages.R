# display homer enrichment results in heatmap

library(ggplot2)
library(reshape2)
library(ggdendro)
library(grid)
library(stringr)

HOMER_result_path = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ChIP-seq/HOMER/output/full_size_peaks/stage-specific/"
score = 0.70


novo_motif_matches = read.table(file = paste0(HOMER_result_path,"/novo_motif_matches.txt"), sep = "\t")
novo_motif_table = read.table(file = paste0(HOMER_result_path,"/novo_motif_table.txt"))

message(paste0("Selecting motifs with p-val <= 10e-10 and score > = ", score))

novo_motif_table = novo_motif_table[novo_motif_table$p_value <= 10^-10,]

novo_motif_table$motif_module = paste(novo_motif_table$motif, novo_motif_table$module, sep = "_")
novo_motif_matches$motif_module = paste(novo_motif_matches$motif, novo_motif_matches$module, sep = "_")

novo_motif_matches = novo_motif_matches[novo_motif_matches$motif_module %in% novo_motif_table$motif_module,]
novo_motif_matches = novo_motif_matches[novo_motif_matches$score >= score,]
novo_motif_matches$p_value = novo_motif_table[match(novo_motif_matches$motif_module,novo_motif_table$motif_module),
                                                    c("p_value")]
novo_motif_matches$log_pvalue = novo_motif_table[match(novo_motif_matches$motif_module,novo_motif_table$motif_module),
                                              c("log_pvalue")]
message("Selecting minimal p-val in repeated motifs per module")

stages = list()
for(s in unique(novo_motif_matches$module)){
stages[[s]] = novo_motif_matches[novo_motif_matches$module == s, ]
stages[[s]]  = stages[[s]] [!duplicated(stages[[s]] $TF),]

}


message("meting data for ggplot")

melted = do.call("rbind",stages)
melted = melt(melted, id.vars = c("module","TF"), measure.vars = c("log_pvalue"),value.name = "p")
melted = melted[,c("module","TF","p")]
# transform to -log pval

melted$p = -(melted$p)



#specific order
melted$module <- factor(melted$module , levels = c("iPSC","DE","GT","PF","PE","EP","EN","BLC"),ordered = T)


# reorder TFs
melted.scaled = melted
melted.scaled <- reshape(melted.scaled, idvar = "TF", timevar = "module", direction = "wide")

order <- as.character(melted.scaled$TF)
melted$TF = factor(x = melted$TF,
                    levels = order, 
                    ordered = TRUE)

# all TFs, group by TF p-val and module
pdf(paste(HOMER_result_path, "/HOMER_enrichment_stages.pdf", sep = ""), height = 12)

p = ggplot( melted, aes(module, TF) ) +
  geom_tile(aes(fill = p)) +
  scale_fill_gradient(high = "black", low = "lightblue", 
                      breaks = c(min(melted$p), max(melted$p))) +
  
  theme_bw() + 
  theme(panel.grid = element_blank())

p
dev.off()


# group by TF family when there is more than 1 family member (those with number in last position)

melted$TF_family = str_replace(melted$TF, "[0-9]{2}$","")
melted$TF_family = str_replace(melted$TF_family, "[0-9]{1}$","")
melted$TF_family = paste(melted$TF_family,"family")

melted$TF = as.character(melted$TF)
melted$TF_family = as.character(melted$TF_family)


# if unique per module after regex, copy original
select_duplicated_TFs = function(melted, modules =unique(melted$module)){
  duplicated = c()
  for (m in modules) {
    module = melted[melted$module == m,]
    module = module[order(module$TF_family),]
    duplicated = append(duplicated, as.character(unique(module[duplicated(module$TF_family),"TF_family"])))
  }
 return(duplicated)
}

duplicated = unique(select_duplicated_TFs(melted))

melted[!(melted$TF_family %in% duplicated),"TF_family"] = melted[!(melted$TF_family %in% duplicated),"TF"]

# reorder TFs
melted.scaled = melted
melted.scaled <- reshape(melted.scaled, idvar = "TF_family", timevar = "module", direction = "wide")

order <- as.character(melted.scaled$TF_family)
melted$TF_family = factor(x = melted$TF_family,
                   levels = order, 
                   ordered = TRUE)
png(paste(HOMER_result_path, "/HOMER_enrichment_sign_stages_family.png", sep = ""), height = 14, width = 6,units = "in",res = 400)


p = ggplot( melted, aes(module, TF_family) ) +
  geom_tile(aes(fill = p)) +
  scale_fill_gradient(high = "black", low = "lightblue", 
                      breaks = c(min(melted$p), max(melted$p))) +
  ylab("TF family") +
  xlab("WGCNA module") +
  labs(fill = "-log10(p-value)")+
  
  theme_bw()

p
dev.off()

# only best match

best_match = novo_motif_table[which(novo_motif_table$p_value <= 10^-10 & novo_motif_table$best_match_score >= score),c("log_pvalue","TF","module")]

melted = melt(best_match,id.vars = c("module","TF"), measure.vars = "log_pvalue", value.name = "p")

# melted$module <- factor(melted$module , levels = c("2","3","10","1"))

# reorder TFs
melted.scaled = melted
melted.scaled <- reshape(melted.scaled, idvar = "TF", timevar = "module", direction = "wide")

order <- as.character(melted.scaled$TF)
melted$TF = factor(x = melted$TF,
                          levels = order, 
                          ordered = TRUE)
melted$p = -(melted$p)

pdf(paste(HOMER_result_path, "/HOMER_enrichment_sign_stages_best_match.pdf", sep = ""), height = 12)

p = ggplot( melted, aes(module, TF) ) +
  geom_tile(aes(fill = p)) +
  scale_fill_gradient(high = "black", low = "lightblue", 
                      breaks = c(min(melted$p), max(melted$p))) +
  
  theme_bw() + 
  theme(panel.grid = element_blank())


p
dev.off()


## how many monogenic diabetes genes
monogenic = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Monogenic_diabetes/Monogenic_no_lipodystrophy.txt")
monogenic = monogenic$V1
## how many T2D loci
T2D = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/HRC_locus_annotation.txt", header = T)
T2D = T2D$Nearest_gene
T2D = unique(T2D)


### taking all TFs over filter
HOMER_result_path = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ChIP-seq/HOMER/output/full_size_peaks/stage-specific/"

score = 0.70


novo_motif_matches = read.table(file = paste0(HOMER_result_path,"/novo_motif_matches.txt"), sep = "\t")
novo_motif_table = read.table(file = paste0(HOMER_result_path,"/novo_motif_table.txt"))

message(paste0("Selecting motifs with p-val <= 10e-10 and score > = ", score))

novo_motif_table = novo_motif_table[novo_motif_table$p_value <= 10^-10,]

novo_motif_table$motif_module = paste(novo_motif_table$motif, novo_motif_table$module, sep = "_")
novo_motif_matches$motif_module = paste(novo_motif_matches$motif, novo_motif_matches$module, sep = "_")

novo_motif_matches = novo_motif_matches[novo_motif_matches$motif_module %in% novo_motif_table$motif_module,]
novo_motif_matches = novo_motif_matches[novo_motif_matches$score >= score,]
novo_motif_matches$p_value = novo_motif_table[match(novo_motif_matches$motif_module,novo_motif_table$motif_module),
                                              c("p_value")]
novo_motif_matches$log_pvalue = novo_motif_table[match(novo_motif_matches$motif_module,novo_motif_table$motif_module),
                                                 c("log_pvalue")]
message("Selecting minimal p-val in repeated motifs per module")

stages = list()
for(s in unique(novo_motif_matches$module)){
  stages[[s]] = novo_motif_matches[novo_motif_matches$module == s, ]
  stages[[s]]  = stages[[s]] [!duplicated(stages[[s]] $TF),]
  
}


message("meting data for ggplot")

melted = do.call("rbind",stages)
melted = melt(melted, id.vars = c("module","TF"), measure.vars = c("log_pvalue"),value.name = "p")
melted = melted[,c("module","TF","p")]
# transform to -log pval

melted$p = -(melted$p)

TF = unique(melted$TF)

## in monogenic diabetes genes
TF[TF %in% monogenic]
length(TF[TF %in% monogenic])

## in T2D closest gene to index variant
TF[TF %in% T2D]
length(TF[TF %in% T2D])

