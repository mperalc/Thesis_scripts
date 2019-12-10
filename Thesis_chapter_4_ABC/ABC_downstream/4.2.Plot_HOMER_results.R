# display homer enrichment results in heatmap

library(ggplot2)
library(reshape2)
library(ggdendro)
library(grid)
library(stringr)
# source("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/ABC_downstream/plot_RNA_TPM.R")

HOMER_result_path = "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/HOMER/output"

score = 0.70


novo_motif_matches = read.table(file = paste0(HOMER_result_path,"/novo_motif_matches.txt"), sep = "\t")
novo_motif_table = read.table(file = paste0(HOMER_result_path,"/novo_motif_table.txt"))

message(paste0("Selecting motifs with p-val <= 10e-10 and score > = ", score))

novo_motif_table = novo_motif_table[novo_motif_table$p_value <= 10^-10,]

novo_motif_table$motif_stage = paste(novo_motif_table$motif, novo_motif_table$stage, sep = "_")
novo_motif_matches$motif_stage = paste(novo_motif_matches$motif, novo_motif_matches$stage, sep = "_")

novo_motif_matches = novo_motif_matches[novo_motif_matches$motif_stage %in% novo_motif_table$motif_stage,]
novo_motif_matches = novo_motif_matches[novo_motif_matches$score >= score,]
novo_motif_matches$p_value = novo_motif_table[match(novo_motif_matches$motif_stage,novo_motif_table$motif_stage),
                                                    "p_value"]
message("Selecting minimal p-val in repeated motifs per stage")

stage_list = list()
for(s in unique(novo_motif_matches$stage)){
  stage_list[[s]] = novo_motif_matches[novo_motif_matches$stage == s, ]
  stage_list[[s]] = stage_list[[s]][!duplicated(stage_list[[s]]$TF),]
}

message("meting data for ggplot")

melted = do.call("rbind",stage_list)
melted = melt(melted, id.vars = c("stage","TF"), measure.vars = "p_value", value.name = "p")
melted = melted[,c("stage","TF","p")]
# transform to -log pval
#melted$p = -log(melted$p)

melted$stage = as.character(melted$stage)


#specific order
melted$stage <- factor(melted$stage , levels = c("iPSC","DE","GT","PF","PE","EP","EN","BLC"), ordered = T)


# reorder TFs
melted.scaled = melted
melted.scaled <- reshape(melted.scaled, idvar = "TF", timevar = "stage", direction = "wide")
melted.scaled = melted.scaled[
  order( melted.scaled[,"p.iPSC"], melted.scaled[,"p.DE"] , melted.scaled[,"p.GT"], melted.scaled[,"p.PF"],
         melted.scaled[,"p.PE"], melted.scaled[,"p.EP"], melted.scaled[,"p.EN"], melted.scaled[,"p.BLC"]),
  ]

order <- as.character(melted.scaled$TF)
melted$TF = factor(x = melted$TF,
                    levels = order, 
                    ordered = TRUE)

# all TFs, group by TF p-val and stage
# pdf(paste(HOMER_result_path, "/HOMER_enrichment_ABC.pdf", sep = ""), height = 15)
png(paste(HOMER_result_path, "/HOMER_enrichment_ABC.png", sep = ""), height = 20, width = 5, units = "in",
    res = 400,type = "cairo")
p = ggplot( melted, aes(stage, TF) ) +
  geom_tile(aes(fill = p)) +
  scale_fill_gradient(trans = "log", 
                      high = "lightblue", low = "black") +
  
  theme_bw() 
print(p)

dev.off()

# Plot TPMs of these genes
# plot_RNA_TPM(unique(melted$TF),"HOMER_enrichment_ABC_TF_TPMs")

# group by TF family when there is more than 1 family member (those with number in last position)

melted$TF_family = str_replace(melted$TF, "[0-9]{2}$","")
melted$TF_family = str_replace(melted$TF_family, "[0-9]{1}$","")

## for some other specific families
melted$TF_family =str_replace(melted$TF_family, "ESRR[a-zA-Z]","ESRR")
melted$TF_family =str_replace(melted$TF_family, "NFY[a-zA-Z]","NFY")
melted$TF_family =str_replace(melted$TF_family, "JUN[a-zA-Z]","JUN")
melted$TF_family =str_replace(melted$TF_family, "ZBTB[a-zA-Z]","ZBTB")
melted$TF_family =str_replace(melted$TF_family, "FOS[a-zA-Z]","FOS")
melted$TF_family =str_replace(melted$TF_family, "HOX[a-zA-Z]","HOX")

melted$TF_family = paste(melted$TF_family,"family")

melted$TF = as.character(melted$TF)
melted$TF_family = as.character(melted$TF_family)


# if unique per stage after regex, copy original
select_duplicated_TFs = function(melted, stage = unique(novo_motif_matches$stage)){
  duplicated = c()
  for (m in stage) {
    stage = melted[melted$stage == m,]
    stage = stage[order(stage$TF_family),]
    duplicated = append(duplicated, as.character(unique(stage[duplicated(stage$TF_family),"TF_family"])))
  }
 return(duplicated)
}

duplicated = unique(select_duplicated_TFs(melted))

melted[!(melted$TF_family %in% duplicated),"TF_family"] = melted[!(melted$TF_family %in% duplicated),"TF"]

# reorder TFs
melted.scaled = melted
melted.scaled <- reshape(melted.scaled, idvar = "TF_family", timevar = "stage", direction = "wide")
melted.scaled = melted.scaled[
  order( melted.scaled[,"p.iPSC"], melted.scaled[,"p.DE"] , melted.scaled[,"p.GT"], melted.scaled[,"p.PF"],
         melted.scaled[,"p.PE"], melted.scaled[,"p.EP"], melted.scaled[,"p.EN"], melted.scaled[,"p.BLC"]),
  ]

order <- as.character(melted.scaled$TF_family)
melted$TF_family = factor(x = melted$TF_family,
                   levels = order, 
                   ordered = TRUE)
# pdf(paste(HOMER_result_path, "/HOMER_enrichment_sign_ABC_family.pdf", sep = ""), height = 15)
png(paste(HOMER_result_path, "/HOMER_enrichment_sign_ABC_family.png", sep = ""), height = 15, width = 5, units = "in",
    res = 400,type = "cairo")

p = ggplot( melted, aes(stage, TF_family) ) +
  geom_tile(aes(fill = p)) +
  scale_fill_gradient(trans = "log", 
                      high = "lightblue", low = "black") +
  
  theme_bw()
print(p)

dev.off()

# only best match

best_match = novo_motif_table[which(novo_motif_table$p_value <= 10^-10 & novo_motif_table$best_match_score >= score),c("p_value","TF","stage")]

melted = melt(best_match,id.vars = c("stage","TF"), measure.vars = "p_value", value.name = "p")
melted$stage = as.character(melted$stage)

melted$stage <- factor(melted$stage , levels = c("iPSC","DE","GT","PF","PE","EP","EN","BLC"), ordered = T)

# reorder TFs
melted.scaled = melted
melted.scaled <- reshape(melted.scaled, idvar = "TF", timevar = "stage", direction = "wide")
melted.scaled = melted.scaled[
  order( melted.scaled[,"p.iPSC"], melted.scaled[,"p.DE"] , melted.scaled[,"p.GT"], melted.scaled[,"p.PF"],
         melted.scaled[,"p.PE"], melted.scaled[,"p.EP"], melted.scaled[,"p.EN"], melted.scaled[,"p.BLC"]),
  ]

 
order <- as.character(melted.scaled$TF)
melted$TF = factor(x = melted$TF,
                          levels = order, 
                          ordered = TRUE)

# pdf(paste(HOMER_result_path, "/HOMER_enrichment_sign_ABC_best_match.pdf", sep = ""), height = 15)
 png(paste(HOMER_result_path, "/HOMER_enrichment_sign_ABC_best_match.png", sep = ""), height = 13, width = 5, units = "in",
      res = 400,type = "cairo")


p = ggplot( melted, aes(stage, TF) ) +
  geom_tile(aes(fill = p)) +
  scale_fill_gradient(trans = "log", 
                      high = "lightblue", low = "black") +
  
  theme_bw()  


print(p)

dev.off()


