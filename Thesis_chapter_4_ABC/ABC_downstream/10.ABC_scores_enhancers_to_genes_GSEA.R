### Analyze enhancer overlaps with cred sets with GSEA
## Ranked per stage by
## PPAf
## ABC Score
## gene expression

library(fgsea)
library(liger)
library(annotate)
library('org.Hs.eg.db')
library(GenomicRanges)
library(ReactomePA)
library(clusterProfiler)

# from targets of enhancers overlapping T2D genes
predicted_targets = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/functional_PPA/PPA_genetic_functional_ABC_fully_annotated.txt",
                               header = T)
# topPPA = 0.5
# select only top PPAf
# predicted_targets = predicted_targets[predicted_targets$PPAf>topPPA,]

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


# get lists per stage
stage = c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

# PPA_ordered_genelist = list()
# for(s in stage){
#   PPA_ordered_genelist[[s]] = subset[subset$cellType == s,c("TargetGene","PPAf")]
# }

# # Reactome pathways for the universe of gene names
# universe = unique(scores_all$TargetGene)
# universe = as.character(universe)
# universe = mapIds(org.Hs.eg.db, universe, 'ENTREZID', 'SYMBOL')
# pathways <- reactomePathways(universe)

# res = list()
# ppa = list()
# for(s in stage){
#   id = mapIds(org.Hs.eg.db,PPA_ordered_genelist[[s]]$TargetGene,'ENTREZID', 'SYMBOL')
#   names =  mapIds(org.Hs.eg.db,PPA_ordered_genelist[[s]]$TargetGene,'ENTREZID', 'SYMBOL')
# 
#   names(id) = id
#   ppa[[s]] = PPA_ordered_genelist[[s]]$PPAf
#   names(ppa[[s]]) = id
#   res[[s]] = iterative.bulk.gsea(ppa,
#                       set.list = pathways)
#   res[[s]]=res[[s]][order(res[[s]]$q.val),]
# }
# 
#   id = mapIds(org.Hs.eg.db,as.character(subset$TargetGene),'ENTREZID', 'SYMBOL')
#   names =  mapIds(org.Hs.eg.db,as.character(subset$TargetGene),'ENTREZID', 'SYMBOL')
# 
#   names(id) = id
#   ppa = subset$PPAf
#   names(ppa) = id
# res = iterative.bulk.gsea(ppa, set.list = pathways, n.rand =  c(1e2,1e3))
# gsea(ppa,pathways[["Insulin processing"]],n.rand = 1000)

### does not make any sense because PPAs are too low to give meaningful analysis
# do hypergeometric test


PPA_ordered_genelist = list()

res = list()
exp = list()
for(s in stage){
  PPA_ordered_genelist[[s]] = subset[subset$cellType == s,c("TargetGene","TargetGeneExpression")]
  id = mapIds(org.Hs.eg.db,as.character(PPA_ordered_genelist[[s]]$TargetGene),'ENTREZID', 'SYMBOL')
  names =  mapIds(org.Hs.eg.db,as.character(PPA_ordered_genelist[[s]]$TargetGene),'ENTREZID', 'SYMBOL')

  names(id) = id
  exp[[s]] = PPA_ordered_genelist[[s]]$TargetGeneExpression
  names(exp[[s]]) = id
  exp[[s]] = exp[[s]][!duplicated(exp[[s]])]
  
  # universe of target gene per stage
  # interpretation: more enriched ingenes related to significant terms than what expected by general enhancer landscape?
  universe = unique(scores_all$TargetGene[scores_all$cellType==s])
  universe = as.character(universe)
  universe = mapIds(org.Hs.eg.db, universe, 'ENTREZID', 'SYMBOL')
  universe = unique(universe)
  pathways <- reactomePathways(universe)
  
  res[[s]] = iterative.bulk.gsea(exp[[s]],
                      set.list = pathways)
  res[[s]]=res[[s]][order(res[[s]]$q.val),]
  res[[s]]$pathways = rownames(res[[s]])
}
gsea(exp[[s]],pathways[["Insulin processing"]],n.rand = 1000)
gsea(exp[["EN"]],pathways[["Regulation of gene expression in beta cells"]],n.rand = 1000)
gsea(exp[["EN"]],pathways[["Insulin processing"]],n.rand = 1000)
plotEnrichment(pathways[["Insulin processing"]],exp[[s]])



fgseaRes <- fgsea(pathways, exp[["BLC"]], minSize=1, maxSize=500, nperm=1000)
head(fgseaRes)


### nothing to see here, the gene targets are too small for this kind of analysis
## try normal hypergeometric
gene_names = lapply(exp,names)
compEnrich_target <- compareCluster(geneCluster   = gene_names,
                                    fun           = "enrichPathway",
                                    pvalueCutoff  = 0.05,
                                    pAdjustMethod = "BH",
                                    readable = T)
dotplot(compEnrich_target, showCategory = 10, title = "Pathway Enrichment: T2D SNP mapped to predicted gene")
compEnrich_target@compareClusterResult
## regulation of beta cell development in PE

## if testing all genes:
gene_names = unlist(gene_names)
gene_names = unique(gene_names)
result = enrichPathway(gene = gene_names,
                       pvalueCutoff  = 0.05,
                       pAdjustMethod = "BH",
                       readable = T )
dotplot(result, showCategory = 10, title = "Pathway Enrichment: T2D SNP mapped to predicted gene")
result
# nothing interesting