# simplify ABC_output
input = list()
for (s in stage) {
  input[[s]] = read.table(paste0("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/",s,"/EnhancerPredictions.txt"),
                          header = T, sep = "\t")
  input[[s]] = input[[s]][c( 1:3,(ncol(input[[s]])-27):ncol(input[[s]]))]
  input[[s]] = input[[s]][which(!input[[s]]$TargetGeneExpression=="NA"),]
  input[[s]]$strand = rep(".",nrow(input[[s]]))
  input[[s]] = input[[s]][c( "chr" ,  "start" ,    "end", "strand","cellType")]
  ## remove those in sex chromosomes
  input[[s]] = input[[s]][input[[s]]$chr %in% paste0("chr",1:22),]
  
  write.table( input[[s]],paste0("/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/to_merge/",s,".txt"),
               col.names = F, row.names = F, quote = F, sep = "\t")

}


