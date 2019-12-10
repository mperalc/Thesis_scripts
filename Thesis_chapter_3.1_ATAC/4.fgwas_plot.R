
"%&%" <- function(a,b) paste0(a,b)
library("dplyr")
library("data.table")
library("ggplot2")
library("gridExtra")

stages <- c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

outFolder = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/fGWAS/new_ATAC_peaks_2018/ATAC_WGCNA_P12_S120_deepSplit2_merge0.30_signedHybrid_rmOutliersT_aug2019_1CPM/"
fgwas.out.dir = outFolder

best.joint.model.file.prefix <- paste0(fgwas.out.dir, "best-joint-model")

my.title <- "Enrichment in T2D signals: ATAC-seq WGCNA modules"
output.prefix <- "WGCNA_ATAC_1CPM"

remove.vec <- c("intron","utr_5","utr_3","exon","transcript","promoter")


# Processing function
annot_process <- function(fname){
  df <- fread(fname)
  df$CI_lo <- gsub("<","",df$CI_lo); df$CI_lo <- gsub("fail",NA,df$CI_lo); df$CI_lo <- as.numeric(df$CI_lo)
  df$CI_hi <- gsub(">","",df$CI_hi); df$CI_hi <- gsub("fail",NA,df$CI_hi); df$CI_hi <- as.numeric(df$CI_hi)
  return(df)
}


# Build data frames for plotting  



update_names <- function(par.df, key.df){
  par.df$parameter <- gsub("_ln","",par.df$parameter)
  for (i in 1:dim(par.df)[1]){
    annot <- par.df$parameter[i]
    if (grepl("distance_tss",annot)==TRUE){
      annot <- "distance_tss"
    }
    plotname <- key.df$PlotName[grepl(annot,key.df$Name)]
    par.df$parameter[i] <- plotname
    print(c(annot,plotname))
  }
  return(par.df)
}

build_param_df <- function(name.vec=NULL){
  files <- list.files(fgwas.out.dir)
  param.files <- files[grepl(".params",files,fixed=TRUE)]
  param.files <- param.files[!grepl("+",param.files,fixed=TRUE)]
  param.files <- param.files[!grepl("drop",param.files,fixed=TRUE)]  
  param.files <- param.files[!grepl("best-joint-model",param.files,fixed=TRUE)]  
  param.files <- param.files[!grepl("fgwas_run_loci-partition",param.files,fixed=TRUE)]  
  out.df <- c()
  for (f in param.files){
    stack.df <- annot_process(fgwas.out.dir%&%f)
    out.df <- rbind(out.df,stack.df)
  }
  out.df <- filter(out.df,parameter!="pi_region")
  out.df <- arrange(out.df,desc(estimate))
  out.df$parameter = gsub("_ln", "", out.df$parameter)
  return(out.df)
}

build_param_best_df <- function(pre){
  
  param.df <- annot_process(pre %&% ".params")
  param.df <- arrange(param.df,desc(estimate))
  param.df <- filter(param.df,parameter!="pi")
  param.df <- filter(param.df,parameter!="pi_region")
  param.df$parameter = gsub("_ln", "", param.df$parameter)
  
  return(param.df)
}


# Annotation plot

an_plot <- function(mydf,mytitl="",mylow=-10,myhigh=10,interval=1){
  mydf <- within(mydf,parameter<-factor(mydf$parameter,levels=rev(mydf$parameter)))
  mydf <- filter(mydf,parameter!="pi_region")
  plt <- ggplot(data=mydf,aes(x=parameter,y=estimate)) +
    geom_hline(yintercept=0,linetype=2) +
    ylab("lnFE") + xlab("Annotation") +
    geom_errorbar(aes(ymin=CI_lo,ymax=CI_hi),width=0.1) +
    geom_point(shape=21,size=1.5,col="black",fill = "steelblue1" )  +
    theme_bw()  +  theme(legend.position = "none",
                         panel.grid.minor=element_blank(),
                         panel.grid.major=element_blank()) +
    scale_y_continuous(breaks=seq(mylow,myhigh,by=interval)) +
    coord_flip(ylim=c(mylow,myhigh)) +
    ggtitle(mytitl)
  return(plt)
}

# Generate plot 1 : Individual annotations that individually meet marginal significance (confidence intervals not overlapping zero)

param.sing.df <- build_param_df()
plt1 <- an_plot(param.sing.df,mytitl=my.title,mylow=-1,myhigh=4,interval=1)
plt1
write.table(x=param.sing.df,file=fgwas.out.dir %&% output.prefix %&% ".fgwas_enrich.separate.txt",
            sep="\t",quote=F,col.names=TRUE,row.names=FALSE)
ggsave(filename=fgwas.out.dir %&% output.prefix %&% ".fgwas_enrich.separate.png",plot=plt1,width=5.5,height=5)


# Generate plot 2 : Joint analysis of annotations in the 'Best' model


param.df <- build_param_best_df(pre=best.joint.model.file.prefix)
plt2 <- an_plot(param.df,mytitl=my.title%&%": Best Joint Model",mylow=-1,myhigh=4,interval=1)
plt2
write.table(x=param.df,file=fgwas.out.dir %&% output.prefix %&% ".fgwas_enrich.joint.txt",
            sep="\t",quote=F,col.names=TRUE,row.names=FALSE)
ggsave(filename=fgwas.out.dir %&% output.prefix %&% ".fgwas_enrich.joint.png",plot=plt2,width=3,height=2)

