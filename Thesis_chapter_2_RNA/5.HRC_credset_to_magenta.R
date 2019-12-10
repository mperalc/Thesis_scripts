### HRC credsets to Magenta input


library(biomaRt)
library(dplyr)
currentDate <- Sys.Date() # to save date in name of output files


# T2D

setwd("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis")

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

region = do.call("rbind",credset)


write.table(region[c("Chr","Pos","PPAg")],file = "HRC_PPAg_for_Magenta.txt",sep = "\t",col.names = F, row.names = F, quote = F)



### whole DIAMANTE GWAS
gwas = read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/HRC_GWAS/EU.2017JulyRel.Summary.bed",header = F)
gwas = gwas[c(1,2,7)]
gwas$V1 = gsub(pattern = "chr",replacement = "",gwas$V1)
write.table(gwas,file = "Diamante_gwas_for_Magenta.txt",sep = "\t",col.names = F, row.names = F, quote = F)
