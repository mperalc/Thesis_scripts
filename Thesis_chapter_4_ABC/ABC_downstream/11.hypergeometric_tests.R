## hypergeometric tests
### enrichment of T2D genes per stage in TF targets

library(tftargets)
library(annotate)
library('org.Hs.eg.db')
library(ggplot2)
library(scales)
library(purrr)

outdir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/hypergeometric/"
# from targets of enhancers overlapping T2D genes
predicted_targets = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/functional_PPA/PPA_genetic_functional_ABC_fully_annotated.txt",
                               header = T)
# topPPA = 0.5
# select only top PPAf
# predicted_targets = predicted_targets[predicted_targets$PPAf>topPPA,]

scores_all = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/all_expressed.txt",
                        header = T)
############################################## fix
# monogenic diabetes genes
# monogenic = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Monogenic_diabetes/Monogenic_updated_2019.txt")
# other_genes = c("NANOG","SOX17","NKX6-1","NKX2-2","PROX1","MAFA","MAFB","ONECUT1","CUX6","TCF7L2")
# to_calculate = unique(c(as.character(monogenic$V1),other_genes))
# TF = to_calculate
##################################

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

## trusted databases
ENCODE= tftargets::ENCODE
TTRUST = tftargets::TRRUST
tf = list()
for (n in unique(c(names(ENCODE),names(TTRUST)))){
  if( is.null(ENCODE[[n]]) == F){
  tf[[n]] = unname(mapIds(org.Hs.eg.db,as.character(ENCODE[[n]]),'SYMBOL','ENTREZID') )
  }
  tf[[n]] = append(tf[[n]] ,TTRUST[[n]])
}

##################
# hypergeometric test with permutations
# FOR every TF
# Total balls in the urn: background dataset (all genes)
# x --> number of white balls drawn without replacement from an urn which contains both black and white balls (TF targets in my tested dataset)
# m --> White balls in the urn: characteristic I'm testing for enrichment (TF targets)
# n --> Black balls in the urn: total-white balls (non-TF targets)
# k --> Number of balls drawn from the urn: depends on the size of the dataset tested 

# I want the p-value of getting (number of TF target genes in my gene  dataset) or more white balls in a sample 
# of size (number of genes in my dataset) from an urn with (number of TF targets) white balls and 
# (total genes TF target genes in that set) black balls.

# 1-phyper(q, m, n, k)

# q in this case is x-1 (diff exp dataset -1), the probability of x or more. It is a quantile defined as the smallest value x such that F(x) ≥ p, 
# where F is the distribution function.
# 


# comparison of probabilities

phyperRandom <- function(myGeneList, myGeneSet, genome){
  myRandomGS <- sample( genome,size=length(myGeneSet) )  # take random sample of genes from whole genome, same size as set tested
  myX <- length(which(myGeneList %in% myRandomGS))  # number of T2D in random set
  myM <- length(which(myGeneList %in% genome)) #  T2D genes in my pool
  myN <- length(genome) - length(which(myGeneList %in% genome))  # genes in genome not T2D 
  myK <- length(myGeneSet) # size of set tested
  return(list(prob=1-phyper(q=myX-1,
                            m=myM, n=myN, k=myK), n.overlap=myX))
}

df_list <- list()
TF = names(tf) # list of TFs to test
overlap_list <- list()  # to write down number of TF targets per stage

all_genes = as.character(unique(scores_all$TargetGene))
for(t in TF){
  message(paste0("working on TF ",t))
  prob_results=list()
  sums_probs=list()
  b_test <- list()
  overlap <- list()
  
  for(s in stage){
    t2d_targets = as.character(unique(subset[subset$cellType == s,"TargetGene"]))
    # p-value for enrichment in TF targets using my gene targets of ABC-predicted enhancers for each stage
    prob_results[[s]]=1-phyper(q=length(which(unique(tf[[t]]) %in% t2d_targets ))-1,
                               m=length(which(unique(tf[[t]]) %in% all_genes)),
                               n=length(all_genes) - length(which(unique(tf[[t]]) %in% all_genes)),
                               k=length(t2d_targets))
    gene_list=tf[[t]][which(unique(tf[[t]]) %in% t2d_targets)]
    
    
    overlap[["total_in_background"]] = length(which(unique(tf[[t]]) %in% all_genes))
    
    overlap[[s]] = length(which(unique(tf[[t]]) %in% t2d_targets))
    
    # same, for random set of genes of same size. random gene set distribution of p-values
    random=list()
    random_overlaps=numeric()
    probab_random=numeric()
    
    for(i in 1:10000){
      random[[i]] <- as.data.frame(phyperRandom( unique(tf[[t]]), t2d_targets, all_genes))  # get probs from random gene sets, for each diff expr set tested
      probab_random[i] <- random[[i]]$prob
      random_overlaps[i] <- random[[i]]$n.overlap
      
    }
    
    
    # hist(random_overlaps)
   
    # empirical p-value 
    
    sums_probs[[s]]$pval=((sum(probab_random<=prob_results[[s]])+1)/(10000+1)) # fraction of probabilities from random draws that are more extreme (smaller) than the one from my tested set
    # +1 in case sum(probab_random<=prob_results[[s]]) ==0
    # that would be the permuted? p-value
    
    # calculate the 95% confidence interval from the binomial distribution
    b_test[[s]]=binom.test(sum(probab_random<=prob_results[[s]]),10000)  # If those fractions are my "successes", then the binomial test can give me a confidence interval
    sums_probs[[s]]$conf.int.low=b_test[[s]]$conf.int[1]
    sums_probs[[s]]$conf.int.up=b_test[[s]]$conf.int[2]
    
    
    print(s)
    
  }
  overlap_list[[t]] <- do.call("rbind", overlap)
  
  df_list[[t]] <- do.call("rbind", sums_probs)
  
  # ######### save tables
  # ordered by p-value
  


   
}
#LYL1
overlap_list2 <- map_df(overlap_list, ~as.data.frame(.x), .id="names")
overlap_list2$stages = rep(c("total_back",stage),nrow(overlap_list2)/9)

df_list2 <- map_df(df_list, ~as.data.frame(.x), .id="names")
df_list2$stages = rep(stage,nrow(df_list2)/8)
df_list2$pval = unlist(df_list2$pval)
df_list2$conf.int.low= unlist(df_list2$conf.int.low)
df_list2$conf.int.up = unlist(df_list2$conf.int.up)

df_list2 = as.data.frame(df_list2)

write.table(df_list2,quote=F,row.names=F,col.names = T,file=paste0(outdir,"/ABC_predicted_targets_from_enhancers_at_T2D_SNPs_TF_enrichments.txt"),sep="\t")
write.table(overlap_list2,quote=F,row.names=F,col.names = T,file=paste0(outdir,"/ABC_predicted_targets_from_enhancers_at_T2D_SNPs_TF_overlaps.txt"),sep="\t")



# plot TFs that show significant enrichment at any stage
df_list2 = df_list2[order(df_list2$conf.int.up,decreasing=F),]
TF_to_keep = unique(df_list2[df_list2$conf.int.up<0.05,"names"])

test = df_list2[df_list2$names %in% TF_to_keep,]

# to_plot = df_list2[df_list2$names %in% TF_to_keep,]
to_plot = rbind(to_plot,test)
to_plot = unique(to_plot)


to_plot$stages = factor(to_plot$stages,levels = stage,ordered = T)
to_plot$pval = -log10(to_plot$pval)

# to_plot =to_plot[
#   order( to_plot[,to_plot$stages == "iPSC"], to_plot[,to_plot$stages == "DE"] , to_plot[,"GT"], to_plot[,"PF"],
#          to_plot[,"PE"], to_plot[,"EP"], to_plot[,"EN"], to_plot[,"BLC"]),
#   ]

p = ggplot(data = to_plot,aes(x=stages,y=names)) + 
  geom_tile(aes(fill=pval)) +
  scale_fill_gradientn(colours=c("pink1","coral3"),
                       values  = rescale(c(min(to_plot$pval), 1.30,1.300001, max(to_plot$pval))),
                       limits = c(1.300001,max(to_plot$pval)),
                       name = "Permuted p-value") +

  
  xlab("Stage")+
  ylab("Enriched transcription factors")

png(paste0(outdir,"/ABC_predicted_targets_from_enhancers_at_T2D_SNPs_TF_enrichments_0.05_threshold.png"),
    height = 10,width = 6,units = "in",res = 400, type = "cairo")
p
dev.off()
  
  