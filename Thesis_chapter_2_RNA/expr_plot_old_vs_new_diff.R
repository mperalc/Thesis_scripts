# plot expression of same gene in Old vs New data
OptimisedDate <- Sys.Date() # to save date in name of output files

library(reshape2)
library(limma)
library(edgeR)
library(ggplot2)

load(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/session_objects/dge_old_and_new_filtered_75bp.xz",verbose=TRUE)  #loading the dge object for conservative counts

donors= c("Ad2.1","Ad3.1","Neo1.1")  # data comes from three donors
stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","ENstage6", "ENstage7") # 8 stages
old_stages=c("iPSC","DE","PGT","PFG","PE","ENstage6") # 6 stages
old_donors= c("Ad2.1","Ad3.4")


# plot with counts per million

cpm=cbind(filtered_combined_commongenes$genes,cpm(filtered_combined_commongenes$counts))



###################### plot
###################################################################################################################################

list=c("NEUROG3","PDX1","INS","MAFA","ABCC8","G6PC2","GCG","SST")


################ 1st part as QC script

plot_long=cpm[match(list,cpm$external_gene_name),]  # extracts from cpm data frame the rows that contain the genes of interest
plot_long=na.omit(plot_long)              # remove NAs, in case there was not a match for a gene

diff=setdiff(list,cpm$external_gene_name)           # which ones in the list are not in the table (probably have a different name)?

#order the columns by stage

nc=plot_long[ , grepl( "iPSC" , names( plot_long ) ) ]  #takes columns whose name matches x

#Do the same for the other stages

for (s in stage[-1])  {

  i= plot_long[ , grepl( s , names( plot_long ) ) ]
  nc=cbind(nc,i)
}

plot_long=cbind(plot_long[c(1:2)],nc)
rm(i,nc)


# melt data for ggplot2

long = melt(plot_long, measure.vars = c(3:ncol(plot_long)))
head(long)
long=long[,c(2:4)]
colnames(long)=c("external_gene_name","variable","value")

# rename stages and samples
stage_2= c("iPSC", "DE", "GT", "PF", "PE", "EP","EN", "BLC")  #shortening EN names EN7= BLC (beta-like cells)
long$Stages <- c(rep(stage_2[1:5],each=5*length(list)),rep(c("EP"),each=3*length(list)),rep(c("EN"),each=5*length(list)),rep(c("BLC"),each=3*length(list)))
long$Stages <- ordered(long$Stages,levels=stage_2)

long$Sample<- c(rep(rep(c("Ad2.1","Ad3.1","Neo1.1","Ad2.1","Ad3.4"),each=1*length(list)),5),
  rep(c("Ad2.1","Ad3.1","Neo1.1"),each=1*length(list)),
  rep(c("Ad2.1","Ad3.1","Neo1.1","Ad2.1","Ad3.4"),each=1*length(list)),
  rep(c("Ad2.1","Ad3.1","Neo1.1"),each=1*length(list)))

long$Experiment= c(rep(rep(c("Optimised","Optimised","Optimised","Previous","Previous"),each=1*length(list),5)),
                   rep(c("Optimised","Optimised","Optimised"),each=1*length(list)),
                   rep(c("Optimised","Optimised","Optimised","Previous","Previous"),each=1*length(list)),
                   rep(c("Optimised","Optimised","Optimised"),each=1*length(list)))

head(long)
long$external_gene_name=factor(long$external_gene_name,levels=list)

# plot
# diaPalette <- c("#C15858", "#7883BA")  # Diabetologia palette red and blue
 diaPalette <- c("#7883BA","#CADAE1")  # Diabetologia palette darkblue and lightblue
 #diaPalette <- c("#96165A","#755A91")  # Diabetologia palette red and purple (as in PCA)
 
#bw <- c("#000000", "#A9A9A9")  # black and white for diabetologia

# png(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2/new_vs_old_cpm/",list,".png",sep=""), type="cairo",
#     width=8,height=5,units="in",res=300,pointsize = 12)
# 
# p <- ggplot(data = long, aes(x = stage, y = value, group=interaction(experiment,sample),col=experiment)) +
#   ggtitle(unique(long$external_gene_name)) +
#   xlab ("Differentiation stages") +
#   ylab ("Expression [CPM]") +
#   expand_limits(y = 0) +
#   geom_hline(yintercept=0,linetype="dashed",size=1,col="#DCDCDC") +
#   geom_line(aes(linetype = sample,col=experiment), size = 1) +
#   geom_point(size=3,aes(col=experiment,shape=experiment)) +
#   scale_colour_manual(values=diaPalette) +
#   theme_bw()+
#   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
#         panel.border=element_rect(size=2),axis.text=element_text(size=12,face="bold"),
#         axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=16,face="bold"),
#         legend.text = element_text(size=11,face="bold"),legend.title = element_text(size=13,face="bold")) 
#  
# 
# print(p)
# dev.off()




# windowsFonts(Helvetica=windowsFont("Helvetica"))  # seems like Helvetica is the default font anyways for this OS

png("/Users/Marta/Documents/WTCHG/DPhil/Plots/new_vs_old_expr_paper_mat_markers_2019.png", type="cairo", 
     width=9,height=10,units="in",res=1000,pointsize = 13)


p <- ggplot(data = long, aes(x = Stages, y = value, group=interaction(Sample,Experiment),col=Experiment)) +
  #ggtitle("Counts per gene, stage & sample") +
  # xlab ("Differentiation stages") +
  ylab ("Expression (CPM)") +
  expand_limits(y = 0) +
  geom_hline(yintercept=0,linetype="dashed",size=1,col="#DCDCDC") +
  geom_line(aes(col=Experiment), size = 1) +
  geom_point(size=3,aes(col=Experiment,shape=Sample)) +
  scale_shape_manual(values=c(15,17,18,16)) +
  scale_fill_manual(values=diaPalette,guide="none") + # guide="none" takes out legend for this alpha
  scale_color_manual(values=diaPalette) +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        panel.border=element_rect(size=2),
        axis.text.x=element_text(size=16,face="bold",family = "Helvetica",angle=45,vjust=0.55),axis.text.y=element_text(size=16,face="bold",family = "Helvetica"),
        axis.title=element_text(size=16,face="bold",family = "Helvetica"),plot.title = element_text(size=16,face="bold",family = "Helvetica"),
        axis.title.x = element_blank(),  axis.title.y = element_text(size=18,face="bold",family = "Helvetica"),
        legend.text = element_text(size=16,face="bold",family = "Helvetica"),legend.title = element_text(size=16,face="bold",family = "Helvetica"),
        strip.background = element_blank(),strip.text= element_text(size=16,face="bold.italic",family = "Helvetica") ) +  # strip controls subfigure titles
  
  facet_wrap(~external_gene_name,scales="free",ncol=2)

plot(p)

dev.off()

#}



