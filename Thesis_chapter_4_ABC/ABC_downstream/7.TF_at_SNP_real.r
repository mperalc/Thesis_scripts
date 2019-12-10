### creating file for HOMER to find bound motifs at SNPs in ABC-predicted enhancers

# findMotifsGenome.pl ERalpha.peaks hg18 MotifOutputDirectory/ -find motif1.motif -size 30 > outputfile.txt

predicted_targets = read.table("../data/enhancers_overlap_SNP/PPA_genetic_functional_ABC_fully_annotated.txt",
                               header = T)

# 
# SNP = as.character(predicted_targets$SNPposition)
# 
# SNP = strsplit(SNP,":")
# SNP = do.call("rbind",SNP)
# SNP = as.data.frame(SNP)
# colnames(SNP) = c("chr","start")
# SNP$chr = as.character(SNP$chr)
# SNP$start = as.character(SNP$start)
# 
# SNP$chr = paste0("chr",SNP$chr)
# SNP$start = as.numeric(SNP$start)
# 
# SNP$end = SNP$start + 1
# 
# SNP$id = predicted_targets$ID
# SNP$strand = rep(".",nrow(SNP))
# SNP = SNP[c("id","chr","start","end","strand")]
# SNP = unique(SNP)
# write.table(SNP,"../data/enhancers_overlap_SNP/SNP_overlap_enhancer_pos.bed", col.names = F, row.names = F,sep = "\t", quote = F)


## run homer and load again
TF = read.table("../ABC_output/HOMER/output/SNP_annotated/SNP_at_enhancer_annotated.txt", header=T,sep = "\t")

TF = TF[abs(TF$Offset) < 7,]
TF = TF[TF$MotifScore > 6.90,]

## load keys of TF names
keys = read.table("../ABC_output/HOMER/PWMs/Genename_Motif.txt", header=F,sep = "\t")
colnames(keys) = c("name","Motif.Name","Inference")
TF = merge(TF, keys, by="Motif.Name",all.x = T)
TF$name = as.character(TF$name)
## concatenate TF per SNP
TF = aggregate(name ~., TF[c(2,7)], toString)
write.table(TF,"../data/enhancers_overlap_SNP/TF_SNP_overlap_enhancer_filtered.txt", col.names = T, row.names = F,sep = "\t", quote = F)

predicted_targets$HOMER_TF = TF[predicted_targets$ID,"name"]

write.table(predicted_targets,"../data/enhancers_overlap_SNP/PPA_genetic_functional_ABC_fully_annotated_HOMER_TFs.txt", col.names = T, row.names = F,sep = "\t", quote = F)
