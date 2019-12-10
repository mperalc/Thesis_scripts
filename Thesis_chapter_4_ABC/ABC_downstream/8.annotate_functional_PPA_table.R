## annotate functional ppa table with tf, gene match etc

ppa = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/functional_PPA/PPA_genetic_functional_ABC.txt",
            header = T, sep = "\t")

outdir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/"

peaks_overlap_snp = read.table(paste(outdir,"credset_SNPs_in_ABC_enhancers_PPA_NoNA_all_annotations.txt", sep = ""),
                               header = T, sep = "\t")
ppa$IndexSNP = gsub("_",":",ppa$IndexSNP)
ppa$id = gsub("_",":",ppa$id)

colnames(ppa)[1] = "SNPposition"

colnames(peaks_overlap_snp)[2] = "SNPposition"
ppa = ppa[c( "SNPposition", "IndexSNP",      "PPA",     
              "difference" )]
colnames(ppa) = c( "SNPposition", "IndexSNP",      "PPAf",      
                       "difference" )
merged = merge(peaks_overlap_snp,ppa, by = "SNPposition")


write.table(merged,"/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/functional_PPA/PPA_genetic_functional_ABC_fully_annotated.txt",
                             col.names = T, row.names = F,
            quote = F,sep = "\t")
