# important script

library(GenomicRanges)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(gridGraphics)
library(dplyr)


options(stringsAsFactors = FALSE)
outdir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/"
consensus = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/for_fGWAS/binary_predictions_ABC.txt",
                    header = T) # this function needs columns "Name","Chr","Start","End","Peaks_present", in that order



credset_names = read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/credset/credset_names.txt",
                           header = F)

credset = list()

for (n in credset_names$V1) {
  credset[[n]] = read.table(
    file = paste(
      "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/credset/",
      n,
      sep = ""
    ),
    header = T
  ) # 380 99% credible sets T2D 2017 HRC (Mahajan et al 2018)
  
}


# allocate each peak in correct credible set
# some will be empty, others will have peak
# integrate peak info and credset SNPs with findOverlaps

peaks_table.gr  = makeGRangesFromDataFrame(consensus,keep.extra.columns = T)

peaks_per_credset = list()
# n="credible_set_Eur_EYA2_20_45317678.txt"
for (n in names(credset)) {
  credset_subset = credset[[n]]
  credset_subset = credset_subset[, c("Chr", "Pos")]
  colnames(credset_subset) = c("chr", "start")
  credset_subset$end = credset_subset$start + 1
  credset_subset$chr = paste("chr", credset_subset$chr, sep = "")
  credset_subset$chr = as.character(credset_subset$chr)
  #there are duplicates within same credible set loci??
  credset_subset = unique(credset_subset)
  rownames(credset_subset) = paste("snp",
                                   credset_subset$chr,
                                   credset_subset$start,
                                   credset_subset$end,
                                   sep = "_")
  credset.gr = makeGRangesFromDataFrame(credset_subset)
  
  # Get all peaks in module overlapping credible set SNPs :
  overlapPeak = subsetByOverlaps(peaks_table.gr, credset.gr)
  hits <-
    findOverlaps(peaks_table.gr, credset.gr) # to add SNP metadata to overlapPeak table
  rsid <-
    CharacterList(split(names(credset.gr)[subjectHits(hits)], queryHits(hits)))
  rsid = paste(rsid, collapse = ",")
  
  mcols(overlapPeak) <- DataFrame(mcols(overlapPeak), rsid)
  
  
  overlapPeak = as.data.frame(overlapPeak)
  overlapPeak = overlapPeak[, c(1:3, 6)]
  colnames(overlapPeak) = c("Chr", "Start", "End", "SNPpos")
  if (nrow(overlapPeak) > 0) {
    overlapPeak$PeakID = paste("peak", overlapPeak$Chr, overlapPeak$Start, sep = "_")
  }
  peaks_per_credset[[n]] = overlapPeak
  
}
rm(credset_subset, overlapPeak)

# Get number of overlapping peaks. 
lapply(peaks_per_credset, nrow)

# how many credsets don't have any overlapping peaks?
# number of credible sets w/o any peaks over their SNPs
peaks_per_credset %>%
  lapply(nrow) %>%
  (function(x)
    x == 0) %>%
  sum

nrow(credset_names) # number credible sets
# 104 of 380 credible sets don't have any overlapping peaks (check this number every time)
(204 / 380) * 100 # 54%

((380-204) / 380) * 100 # 46% do have (176)
peaks_per_credset %>%
  lapply(nrow) %>%
  unlist %>%
  sum
# 1076 active enhancers in total in all credible sets.


#### get PPA for all SNPs under peaks in every credible set
# list of credible set
# SNP overlap peak
# snp, ppa, peak
# if no peak overlapping - "none" under peak col
# cbind and save with row names
SNPs_in_peaks = list()

for (n in names(credset)) {
  
  credset_subset = credset[[n]]
  credset_subset = credset_subset[, c("Chr", "Pos")]
  colnames(credset_subset) = c("chr", "start")
  credset_subset$end = credset_subset$start + 1
  credset_subset$chr = paste("chr", credset_subset$chr, sep = "")
  credset_subset$chr = as.character(credset_subset$chr)
  credset_subset = unique(credset_subset)
  rownames(credset_subset) = paste("snp",
                                   credset_subset$chr,
                                   credset_subset$start,
                                   credset_subset$end,
                                   sep = "_")
  credset.gr = makeGRangesFromDataFrame(credset_subset)
  
  # Get credible set SNPs overlapping peaks in module :
  overlapSNP = subsetByOverlaps(credset.gr, peaks_table.gr)
  hits <-
    findOverlaps(credset.gr, peaks_table.gr) # to add SNP metadata to overlapPeak table
  peakid =  unlist(split(peaks_table.gr$Name[subjectHits(hits)], queryHits(hits)))
  
  
  mcols(overlapSNP) <- DataFrame(mcols(overlapSNP), peakid)
  
  overlapSNP = as.data.frame(overlapSNP)
  if (nrow(overlapSNP) == 0) {
    credset[[n]]$IndexSNP = paste(credset[[n]]$Chr, credset[[n]]$Pos, sep =
                                    "_")
    overlapSNP = credset[[n]]
    overlapSNP$PeakID = rep("none", nrow(overlapSNP))
    SNPs_in_peaks[[n]] = overlapSNP
  }
  else {
    overlapSNP = overlapSNP[, c(1:2, 6)]
    colnames(overlapSNP) = c("Chr", "Pos", "PeakID")
    
    overlapSNP$IndexSNP = paste(gsub("\\chr*", "", overlapSNP$Chr), overlapSNP$Pos, sep =
                                  "_")
    credset[[n]]$IndexSNP = paste(credset[[n]]$Chr, credset[[n]]$Pos, sep =
                                    "_")
    
    SNPs_in_peaks[[n]] = merge(credset[[n]], overlapSNP[,c(3,4)], by = "IndexSNP", all = T)
    SNPs_in_peaks[[n]][is.na(SNPs_in_peaks[[n]]$PeakID),"PeakID"] = "none"
  }
}
SNPs_in_peaks = do.call("rbind",SNPs_in_peaks)
SNPs_in_peaks$credset = rownames(SNPs_in_peaks)
# the numbered credset names in the SNPs_in_peaks table do NOT denote order of SNPs by PPA!

# Add data on peak per stage
colnames(consensus)[1] = "PeakID"

merged_table <- merge(SNPs_in_peaks, consensus, by.x = "PeakID", by.y = "PeakID",  all.x = T)
merged_table = merged_table[c(6,2,3:5,1,7:9,10:17)]
merged_table$SNP_pos = paste(merged_table$Chr.x,merged_table$Pos,sep = ":")
merged_table = merged_table[c(1,2,18,3:17)]

write.table(merged_table,
            paste(outdir,"credset_SNPs_in_ABC_enhancers_PPA.txt", sep = ""),
            sep = "\t", quote = F, row.names = F, col.names = T)


## write table that has only peaks overlapping SNPs
peaks_overlap_SNP = merged_table[!is.na(merged_table$Chr.y),]
peaks_overlap_SNP$PeakID = paste("peak",peaks_overlap_SNP$Chr.y,peaks_overlap_SNP$Start,sep = "_")

## percentage of peaks overlapping SNPs
(length(unique(peaks_overlap_SNP$PeakID ))/nrow(consensus))*100
# 2.80% of enhancers overlap SNPs, 852 out of 30,341

## snps overlapping peaks
(length(unique(peaks_overlap_SNP$IndexSNP))/length(unique(merged_table$IndexSNP)))*100
#1.8% of SNPs
# 1924
## snps with PPA over or equal 10%
sum(peaks_overlap_SNP$PPAg>=0.1)
## snps with PPA over or equal 50%
sum(peaks_overlap_SNP$PPAg>=0.5)

write.table(peaks_overlap_SNP,
            paste(outdir,"credset_SNPs_in_ABC_enhancers_PPA_NoNA.txt", sep = ""),
            sep = "\t", quote = F, row.names = F, col.names = T)
## merge with SNP ids

ids = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/HRC_credset.snp_ann.txt")
colnames(ids) = c("ID","SNP_pos","locus")
# peaks_overlap_SNP$IndexSNP = gsub("_",":",peaks_overlap_SNP$IndexSNP)
peaks_overlap_SNP = merge(peaks_overlap_SNP,ids,by="SNP_pos")
peaks_overlap_SNP = peaks_overlap_SNP [c(
  "credset",
  "ID",
  "SNP_pos",
  "PPAg",
  "PeakID",
  "Chr.y",  "Start",  "End",
  "iPSC",  "DE"  ,  "GT" ,  "PF"  ,  "PE" ,  "EP" ,  "EN" ,  "BLC"
)]
colnames(peaks_overlap_SNP)[colnames(peaks_overlap_SNP)=="Chr.y"] = "Chr"

# save SNPids for TF binding prediction
write.table(unique(peaks_overlap_SNP$ID), 
            paste(outdir,"credset_SNPs_in_ABC_enhancers_rsID.txt", sep = ""),
            sep = "\t", quote = F, row.names = F, col.names = F)
## add info of target gene of overlapping enhancer. Are they acting on closest gene?yes/no 
targets = read.csv(paste0(outdir,"ABC_gene_closest_and_target.csv"))

#each enhancer maps to several genes, so I must concatenate all targets per enhancer
target_gene = aggregate(TargetGene ~., targets[c(1:3,10,15)], toString)
closest_gene = aggregate(closest_SYMBOL ~., targets[c(1:3,10,24)], toString)
targets = merge(target_gene,closest_gene[c(4,5)],by="name",all.x = T)

# same region but slightly different per stage: need to combine those targets and closest for the SNP table


targets.gr = GRanges(targets)
peaks_overlap_SNP.gr = GRanges(peaks_overlap_SNP)
overlaps = findOverlaps(peaks_overlap_SNP.gr,targets.gr)
overlaps = as.data.frame(overlaps)

#dummyt variables
peaks_overlap_SNP$target_Gene= rep("NA",nrow(peaks_overlap_SNP))
peaks_overlap_SNP$closest_Gene= rep("NA",nrow(peaks_overlap_SNP))

## fill in
for(q in  unique(overlaps$queryHits)){
  t = targets[overlaps[overlaps$queryHits==q,"subjectHits"],"TargetGene"]
  c = targets[overlaps[overlaps$queryHits==q,"subjectHits"],"closest_SYMBOL"]
  # unlist, remove whitespaces
  t = unlist(strsplit(t, ","))
  t = unique(gsub(" ", "", t))
  c = unlist(strsplit(c, ","))
  c = unique(gsub(" ", "", c))
  
  peaks_overlap_SNP$target_Gene[q]=paste(t,collapse = ",")
  peaks_overlap_SNP$closest_Gene[q]=paste(c,collapse = ",")
  
}

#Are they acting on the gene giving locus name? yes/no

### simplify credible set name

peaks_overlap_SNP$credset = sapply(strsplit(peaks_overlap_SNP$credset, "_"), function(x) x[4])

# match of target to credible set locus
peaks_overlap_SNP$match_locus = rep("FALSE",nrow(peaks_overlap_SNP))
for(n in 1:nrow(peaks_overlap_SNP)){
peaks_overlap_SNP$match_locus[n] = grepl(peaks_overlap_SNP$credset[n],peaks_overlap_SNP$target_Gene[n])
}
peaks_overlap_SNP$match_locus = as.logical(peaks_overlap_SNP$match_locus)
sum(peaks_overlap_SNP$match_locus) # 874 SNPs were it matches
(sum(peaks_overlap_SNP$match_locus) / nrow(peaks_overlap_SNP)) * 100
# 38%
# percentage that match to closest gene
match_closest = mapply(grepl,peaks_overlap_SNP$closest_Gene,peaks_overlap_SNP$target_Gene)
sum(match_closest)/nrow(peaks_overlap_SNP) # 70%

# acting on multiple genes:
how_many = lapply(strsplit(peaks_overlap_SNP$target_Gene,","),length)
sum(do.call("rbind",how_many) > 1) #  1164 act on more than one, 1164 /2277 = 51%
sum(do.call("rbind",how_many) == 1) #  1113 act on 1 exactly

# annotate binding TFs

tf = read.table(paste0(outdir,"SNP2TFBS/SNP2TFBS_matches_ABC_enhancers.txt"),header = F)
tf$V7 = gsub("TF=","",tf$V7)
colnames(tf) = c("chr","start","end","ref","alt","Nmatch","TF","ScoreDiff","ID")  
tf = tf[c(4:5,7:9)]
peaks_overlap_SNP = merge(peaks_overlap_SNP,tf,by="ID",all.x = T)


## merge with rsID info
rsid = read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/HRC_credset.snp_ann.txt",
                  header = F)
colnames(rsid)=  c("rsID","IndexSNP","signal")
peaks_overlap_SNP = merge(peaks_overlap_SNP,rsid, all.x = T)

write.table(peaks_overlap_SNP,
            paste(outdir,"credset_SNPs_in_ABC_enhancers_PPA_NoNA_all_annotations.txt", sep = ""),
            sep = "\t", quote = F, row.names = F, col.names = T)



############################################################################
## fix to get index snp
# merge info from this table with the one from functional ppa