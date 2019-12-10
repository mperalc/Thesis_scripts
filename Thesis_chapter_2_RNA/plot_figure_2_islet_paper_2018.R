# Diabetologia Figure 2 version with dots instead of bars

# plot barplot with means and CI for Diabetologia
# Islet diff paper 2017

##
setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list")
datos=read.table("T2D_FG_for_Diabetologia.txt",stringsAsFactors = F,header = T)  # read data

head(datos)
colnames(datos)=c("Sample","Trait","GWAS loci", "Permuted p-value","Low CI","High CI", "Distance")    
datos$`Permuted p-value`=as.numeric(datos$`Permuted p-value`)
datos$`Low CI`=as.numeric(datos$`Low CI`)
datos$`High CI`=as.numeric(datos$`High CI`)
datos$Distance=factor(datos$Distance,levels = unique(datos$Distance),ordered = T)
datos$Trait=factor(datos$Trait)
datos$`GWAS loci`=as.factor(datos$`GWAS loci`)
####


library(ggplot2)
library(reshape2)  # to modify dataframes for ggplot2
library(ggsignif)  # significance lines
library(grid)  # tick marks on the inside
library(RColorBrewer)
library(gridExtra)
library(gtable)
library(stringr)  # to wrap text of labels

df_log=datos
df_log$`Permuted p-value`=-log10(datos$`Permuted p-value`)
df_log$`Low CI`=-log10(datos$`Low CI`)
df_log$`High CI`=-log10(datos$`High CI`)
df_log <- data.frame(df_log,bin=paste(df_log[,2],df_log[,3],sep="."))

df_log$bin = factor(df_log$bin,levels=c("FG.All","T2D.Beta-cell_function","T2D.All"))
#rename so I don't have to change later
levels(df_log$bin) <- c("Fasting glucose","T2D (beta-cell)","T2D (all)")
df_log$Distance=paste(df_log$Distance,"kb")
df_log$Distance=factor(df_log$Distance,levels=c("0 kb", "50 kb", "100 kb","200 kb", "500 kb"))
                          
#pal1=c("#96165a","#b2278b","#d09bef","#d794e8","#c4b1ed")   # in the same scale than the BLC stage from the pca plot
pal2=c("#96165a","#bf3b5c","#d1615e","#d8997b","#efddb1")   # in the same scale than the BLC stage from the pca plot

tiff("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/hypergeometric_test_Diabetologia_logscale.tiff", type="cairo",
    width=8,height=5,units="in",res=1200,pointsize = 12,compression = "lzw")

p1 <- ggplot(data = df_log, 
            aes(x = bin, 
                y = `Permuted.p.value`, 
                fill=Distance)) +  
  # geom_errorbar(aes(ymin=`Low.CI`, ymax=`High.CI`),
  #              # color="black",
  #               position =position_dodge(0.7),
  #              width=0.2,size=0.5) +
  geom_point(stat="identity",aes(fill=Distance),
             position=position_dodge(width = 0.7),
             size=6,
             pch=21) +  # this + fill gives coloured fill and black outline+ error bars
  
  # scale_fill_brewer("Distance (kb)",type = "seq", palette = "YlOrBr", direction = -1)  +
 # scale_color_manual(values="black") +
  scale_fill_manual(values = pal2) +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05),linetype="dashed",size=0.8,col="grey") +
  scale_y_continuous(limits=c(0,5),
                     expand=c(0,0),
                     expression(-log[10]~"permuted ("~italic(p)~"-value)"))  + # set breaks and take out space under x axis, add sec axis
  xlab("") +
  #ylab("-log10 (p-value)") +
  theme( # no background grids
    # panel.border=element_rect(size=1),  # border as rectangle
    panel.border=element_blank(), # no border
    axis.line.x = element_line(size = 0.30, linetype = "solid", colour = "black"),# axis lines
    axis.line.y = element_line(size = 0.30, linetype = "solid", colour = "black"),
    panel.grid = element_blank(),  # no grid lines
    axis.text.x = element_text(size=16 ,colour="black"),
    axis.text.y =element_text(size=16 ,colour ="black"),  # bold ticks on y, leave margin for inward ticks
    axis.title.x=element_text(size=16 ,colour="black"), 
    axis.title.y=element_text(size=16 ,colour="black"), # no title on x or y
    legend.title = element_text(size=16 ,colour="black"),
    legend.text = element_text(size=15 ,colour="black"))
plot(p1)
dev.off()

#### all stages hypergeometric test 0kb
datos=read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/T2D_FG_for_Diabetologia_allstages_0kb.txt",
                 stringsAsFactors = F,header = T)


head(datos)
colnames(datos)=c("Trait","GWAS loci", "Permuted p-value","Low CI","High CI", "Stage")    
datos$`Permuted p-value`=as.numeric(datos$`Permuted p-value`)
datos$`Low CI`=as.numeric(datos$`Low CI`)
datos$`High CI`=as.numeric(datos$`High CI`)
datos$Stage=factor(datos$Stage,levels = unique(datos$Stage),ordered = T)
datos$Trait=factor(datos$Trait)
datos$`GWAS loci`=as.factor(datos$`GWAS loci`)



df_log=datos
df_log$`Permuted p-value`=-log10(datos$`Permuted p-value`)
df_log$`Low CI`=-log10(datos$`Low CI`)
df_log$`High CI`=-log10(datos$`High CI`)
df_log <- data.frame(df_log,bin=paste(df_log[,1],df_log[,2],sep="."))

df_log$bin = factor(df_log$bin,levels=c("FG.All","T2D.Beta-cell_function","T2D.All"))
#rename so I don't have to change later
levels(df_log$bin) <- c("Fasting glucose","T2D (beta-cell)","T2D (all)")


diaPalette <- c("#000000","#CADAE8", "#7883BA", "#755A91", "#CC85B1",   "#F4B8B0",
                "#96665A", "#96165A")  # Diabetologia palette


# 
# png("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/hypergeometric_test_Diabetologia_logscale_allstages_0kb.png", type="cairo",
#     width=8,height=4,units="in",res=500,pointsize = 12)
tiff("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/hypergeometric_test_Diabetologia_logscale_allstages_0kb.tiff", type="cairo",
     width=8,height=5,units="in",res=1200,pointsize = 12,compression = "lzw")


p2 <- ggplot(data = df_log, 
            aes(x = bin, 
                y = `Permuted.p.value`, 
                fill=Stage)) +
  
 # geom_bar(stat="identity", color="black", position=position_dodge()) +
  # geom_errorbar(aes(ymin=`Low.CI`, ymax=`High.CI`),
  #               # color="black",
  #               position =position_dodge(0.7),
  #               width=0.2,size=0.5) +
  geom_point(stat="identity",aes(fill=Stage),
             position=position_dodge(width = 0.7),
             size=6,
             pch=21) +  # this + fill gives coloured fill and black outline+ error bars
  #scale_fill_brewer("Stage",type = "seq", palette = "YlOrBr", direction = 1)  +
  scale_fill_manual(values=diaPalette) +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05),linetype="dashed",size=0.8,col="grey") +
  scale_y_continuous(limits=c(0,5),
                     expand=c(0,0),
                     expression(-log[10]~"permuted ("~italic(p)~"-value)"))  + # set breaks and take out space under x axis, add sec axis
  xlab("") +
  #ylab("-log10 (p-value)") +
  theme( # no background grids
    #panel.border=element_rect(size=1),  # border as rectangle
    panel.border=element_blank(), # no border
    axis.line.x = element_line(size = 0.30, linetype = "solid", colour = "black"),# axis lines
    axis.line.y = element_line(size = 0.30, linetype = "solid", colour = "black"),
    panel.grid = element_blank(),  # no grid lines
    axis.text.x = element_text(size=16 ,colour="black"),
    axis.text.y =element_text(size=16 ,colour ="black"),  # bold ticks on y, leave margin for inward ticks
    axis.title.x=element_text(size=16 ,colour="black"), 
    axis.title.y=element_text(size=16 ,colour="black"), # no title on x or y
    legend.title = element_text(size=16 ,colour="black"),
    legend.text = element_text(size=15 ,colour="black"))
plot(p2)
dev.off()


## Now the plot of the GSEA results

setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/input_for_GSEA/")
datos=read.table("GSEA_for_islet_paper.txt",stringsAsFactors = F,header = T)  # read data


head(datos)
colnames(datos)=c("Type","Stage","q value") 
datos$`q value`=as.numeric(datos$`q value`)
datos$Stage=factor(datos$Stage,levels = unique(datos$Stage),ordered = T)
datos$`q value`=as.numeric(datos$`q value`)
datos[datos$Type=="GWAS.ranked","Type"]<- "Ranked T2D GWAS list"
datos[datos$Type=="Diff.exp.ranked","Type"]<- "Ranked differentially expressed genes"
datos$Type=factor(datos$Type,levels=unique(datos$Type),ordered=T)

df_log=datos
df_log$`q value`=-log10(datos$`q value`)



# png("/Users/Marta/Documents/WTCHG/DPhil/Data/input_for_GSEA/GSEA_q.values.png", type="cairo",
#     width=11,height=5,units="in",res=500,pointsize = 12)
tiff("/Users/Marta/Documents/WTCHG/DPhil/Data/input_for_GSEA/GSEA_q.values.tiff", type="cairo",
     width=11,height=5,units="in",res=1200,pointsize = 12,compression = "lzw")

p3 <- ggplot(data = df_log, 
            aes(x = Type, 
                y = `q value`, 
                fill=Stage)) +
  
 # geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_point(stat="identity",aes(fill=Stage),
             position=position_dodge(width = 0.7),
             size=6,
             pch=21)+
  #  scale_fill_brewer("Distance (kb)",type = "seq", palette = "YlOrBr", direction = 1)  +
  scale_fill_manual(name="Stage",values=diaPalette) +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05),linetype="dashed",size=0.8,col="grey") +
  scale_y_continuous(limits=c(0,3),
                     expand=c(0,0),
                     expression(-log[10]~"("~italic(q)~"-value)"))  + # set breaks and take out space under x axis, add sec axis
  xlab("") +
  #ylab("-log10 (p-value)") +
  theme( # no background grids
    #panel.border=element_rect(size=1),  # border as rectangle
    panel.border=element_blank(), # no border
    axis.line.x = element_line(size = 0.30, linetype = "solid", colour = "black"),# axis lines
    axis.line.y = element_line(size = 0.30, linetype = "solid", colour = "black"),
    panel.grid = element_blank(),  # no grid lines
    axis.text.x = element_text(size=16 ,colour="black"),
    axis.text.y =element_text(size=16 ,colour ="black"),  # bold ticks on y, leave margin for inward ticks
    axis.title.x=element_text(size=16 ,colour="black"), 
    axis.title.y=element_text(size=16 ,colour="black"), # no title on x or y
    legend.title = element_text(size=16 ,colour="black"),
    legend.text = element_text(size=15 ,colour="black"))
plot(p3)
dev.off()

# 
# multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#   library(grid)
#   
#   # Make a list from the ... arguments and plotlist
#   plots <- c(list(...), plotlist)
#   
#   numPlots = length(plots)
#   
#   # If layout is NULL, then use 'cols' to determine layout
#   if (is.null(layout)) {
#     # Make the panel
#     # ncol: Number of columns of plots
#     # nrow: Number of rows needed, calculated from # of cols
#     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
#                      ncol = cols, nrow = ceiling(numPlots/cols))
#   }
#   
#   if (numPlots==1) {
#     print(plots[[1]])
#     
#   } else {
#     # Set up the page
#     grid.newpage()
#     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
#     
#     # Make each plot, in the correct location
#     for (i in 1:numPlots) {
#       # Get the i,j matrix positions of the regions that contain this subplot
#       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
#       
#       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
#                                       layout.pos.col = matchidx$col))
#     }
#   }
# }

# a version with better alignment that needs gridExtra and gtable
AlignPlots <- function(...) {
  LegendWidth <- function(x) x$grobs[[8]]$grobs[[1]]$widths[[4]]
  
  plots.grobs <- lapply(list(...), ggplotGrob)
  
  max.widths <- do.call(unit.pmax, lapply(plots.grobs, "[[", "widths"))
  plots.grobs.eq.widths <- lapply(plots.grobs, function(x) {
    x$widths <- max.widths
    x
  })
  
  legends.widths <- lapply(plots.grobs, LegendWidth)
  max.legends.width <- do.call(max, legends.widths)
  plots.grobs.eq.widths.aligned <- lapply(plots.grobs.eq.widths, function(x) {
    if (is.gtable(x$grobs[[8]])) {
      x$grobs[[8]] <- gtable_add_cols(x$grobs[[8]],
                                      unit(abs(diff(c(LegendWidth(x),
                                                      max.legends.width))),
                                           "mm"))
    }
    x
  })
  
  plots.grobs.eq.widths.aligned
}


plots <- AlignPlots(p3+ggtitle("a")+scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) ,
                    p2+ggtitle("b") )
fig=do.call(grid.arrange, plots)

#fig=multiplot( p3+ggtitle("a)"), p2+ggtitle("b)")+theme(legend.position = "none"),p1+ggtitle("c)"), cols=1)
#fig=multiplot( p3+ggtitle("a)"), p2+ggtitle("b)"),p1+ggtitle("c)"), cols=1)

tiff("/Users/Marta/Documents/WTCHG/DPhil/Islet\ diff\ paper\ Marta/fig2.tiff", type="cairo",
     width=9,height=6,units="in",res=600,pointsize = 12,compression = "none")
plot(fig)
dev.off()

tiff("/Users/Marta/Documents/WTCHG/DPhil/Islet\ diff\ paper\ Marta/esm_fig_5.tiff", type="cairo",
     width=7,height=3.5,units="in",res=600,pointsize = 12,compression = "none")
plot(p1)
dev.off()






# Diabetologia Figure 2 version with bars instead of plots

# plot barplot with means and CI for Diabetologia
# Islet diff paper 2017

##
setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list")
datos=read.table("T2D_FG_for_Diabetologia.txt",stringsAsFactors = F,header = T)  # read data

head(datos)
colnames(datos)=c("Sample","Trait","GWAS loci", "Permuted p-value","Low CI","High CI", "Distance")    
datos$`Permuted p-value`=as.numeric(datos$`Permuted p-value`)
datos$`Low CI`=as.numeric(datos$`Low CI`)
datos$`High CI`=as.numeric(datos$`High CI`)
datos$Distance=factor(datos$Distance,levels = unique(datos$Distance),ordered = T)
datos$Trait=factor(datos$Trait)
datos$`GWAS loci`=as.factor(datos$`GWAS loci`)
####




df_log=datos
df_log$`Permuted p-value`=-log10(datos$`Permuted p-value`)
df_log$`Low CI`=-log10(datos$`Low CI`)
df_log$`High CI`=-log10(datos$`High CI`)
df_log <- data.frame(df_log,bin=paste(df_log[,2],df_log[,3],sep="."))

df_log$bin = factor(df_log$bin,levels=c("FG.All","T2D.Beta-cell_function","T2D.All"))
#rename so I don't have to change later
levels(df_log$bin) <- c("Fasting glucose","T2D (beta-cell)","T2D (all)")
df_log$Distance=paste(df_log$Distance,"kb")
df_log$Distance=factor(df_log$Distance,levels=c("0 kb", "50 kb", "100 kb","200 kb", "500 kb"))

#pal1=c("#96165a","#b2278b","#d09bef","#d794e8","#c4b1ed")   # in the same scale than the BLC stage from the pca plot
pal2=c("#96165a","#bf3b5c","#d1615e","#d8997b","#efddb1")   # in the same scale than the BLC stage from the pca plot

# tiff("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/hypergeometric_test_Diabetologia_logscale.tiff", type="cairo",
#      width=8,height=5,units="in",res=1200,pointsize = 12,compression = "lzw")

p4 <- ggplot(data = df_log, 
             aes(x = bin, 
                 y = `Permuted.p.value`, 
                 fill=Distance)) +  
  geom_bar(stat="identity", color="black", position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin=`Low.CI`, ymax=`High.CI`),
                # color="black",
                position =position_dodge(0.9),
                width=0.2,size=0.5) +
 
  # scale_fill_brewer("Distance (kb)",type = "seq", palette = "YlOrBr", direction = -1)  +
  # scale_color_manual(values="black") +
  scale_fill_manual(values = pal2) +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05),linetype="dashed",size=0.8,col="grey") +
  scale_y_continuous(limits=c(0,5),
                     expand=c(0,0),
                     expression(-log[10]~"permuted ("~italic(p)~"-value)"))  + # set breaks and take out space under x axis, add sec axis
  xlab("") +
  #ylab("-log10 (p-value)") +
  theme( # no background grids
    # panel.border=element_rect(size=1),  # border as rectangle
    panel.border=element_blank(), # no border
    axis.line.x = element_line(size = 0.30, linetype = "solid", colour = "black"),# axis lines
    axis.line.y = element_line(size = 0.30, linetype = "solid", colour = "black"),
    panel.grid = element_blank(),  # no grid lines
    axis.text.x = element_text(size=16 ,colour="black"),
    axis.text.y =element_text(size=16 ,colour ="black"),  # bold ticks on y, leave margin for inward ticks
    axis.title.x=element_text(size=16 ,colour="black"), 
    axis.title.y=element_text(size=16 ,colour="black"), # no title on x or y
    legend.title = element_text(size=16 ,colour="black"),
    legend.text = element_text(size=15 ,colour="black"))
plot(p4)
# dev.off()

#### all stages hypergeometric test 0kb
datos=read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/T2D_FG_for_Diabetologia_allstages_0kb.txt",
                 stringsAsFactors = F,header = T)


head(datos)
colnames(datos)=c("Trait","GWAS loci", "Permuted p-value","Low CI","High CI", "Stage")    
datos$`Permuted p-value`=as.numeric(datos$`Permuted p-value`)
datos$`Low CI`=as.numeric(datos$`Low CI`)
datos$`High CI`=as.numeric(datos$`High CI`)
datos$Stage=factor(datos$Stage,levels = unique(datos$Stage),ordered = T)
datos$Trait=factor(datos$Trait)
datos$`GWAS loci`=as.factor(datos$`GWAS loci`)



df_log=datos
df_log$`Permuted p-value`=-log10(datos$`Permuted p-value`)
df_log$`Low CI`=-log10(datos$`Low CI`)
df_log$`High CI`=-log10(datos$`High CI`)
df_log <- data.frame(df_log,bin=paste(df_log[,1],df_log[,2],sep="."))

df_log$bin = factor(df_log$bin,levels=c("FG.All","T2D.Beta-cell_function","T2D.All"))
#rename so I don't have to change later
levels(df_log$bin) <- c("Fasting glucose","T2D (beta-cell)","T2D (all)")


diaPalette <- c("#000000","#CADAE8", "#7883BA", "#755A91", "#CC85B1",   "#F4B8B0",
                "#96665A", "#96165A")  # Diabetologia palette


# 
# png("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/hypergeometric_test_Diabetologia_logscale_allstages_0kb.png", type="cairo",
#     width=8,height=4,units="in",res=500,pointsize = 12)
# tiff("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/hypergeometric_test_Diabetologia_logscale_allstages_0kb.tiff", type="cairo",
#      width=8,height=5,units="in",res=1200,pointsize = 12,compression = "lzw")
# 

p5 <- ggplot(data = df_log, 
             aes(x = bin, 
                 y = `Permuted.p.value`, 
                 fill=Stage)) +
  
   geom_bar(stat="identity", color="black", position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin=`Low.CI`, ymax=`High.CI`),
                # color="black",
                position =position_dodge(0.9),
                width=0.2,size=0.5) +
 
  #scale_fill_brewer("Stage",type = "seq", palette = "YlOrBr", direction = 1)  +
  scale_fill_manual(values=diaPalette) +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05),linetype="dashed",size=0.8,col="grey") +
  scale_y_continuous(limits=c(0,5),
                     expand=c(0,0),
                     expression(-log[10]~"permuted ("~italic(p)~"-value)"))  + # set breaks and take out space under x axis, add sec axis
  xlab("") +
  #ylab("-log10 (p-value)") +
  theme( # no background grids
    #panel.border=element_rect(size=1),  # border as rectangle
    panel.border=element_blank(), # no border
    axis.line.x = element_line(size = 0.30, linetype = "solid", colour = "black"),# axis lines
    axis.line.y = element_line(size = 0.30, linetype = "solid", colour = "black"),
    panel.grid = element_blank(),  # no grid lines
    axis.text.x = element_text(size=16 ,colour="black"),
    axis.text.y =element_text(size=16 ,colour ="black"),  # bold ticks on y, leave margin for inward ticks
    axis.title.x=element_text(size=16 ,colour="black"), 
    axis.title.y=element_text(size=16 ,colour="black"), # no title on x or y
    legend.title = element_text(size=16 ,colour="black"),
    legend.text = element_text(size=15 ,colour="black"))
 plot(p5)
# dev.off()


## Now the plot of the GSEA results

setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/input_for_GSEA/")
datos=read.table("GSEA_for_islet_paper.txt",stringsAsFactors = F,header = T)  # read data


head(datos)
colnames(datos)=c("Type","Stage","q value") 
datos$`q value`=as.numeric(datos$`q value`)
datos$Stage=factor(datos$Stage,levels = unique(datos$Stage),ordered = T)
datos$`q value`=as.numeric(datos$`q value`)
datos[datos$Type=="GWAS.ranked","Type"]<- "Ranked GWAS list"
datos[datos$Type=="Diff.exp.ranked","Type"]<- "Ranked differentially expressed genes"
datos$Type=factor(datos$Type,levels=unique(datos$Type),ordered=T)

df_log=datos
df_log$`q value`=-log10(datos$`q value`)



# png("/Users/Marta/Documents/WTCHG/DPhil/Data/input_for_GSEA/GSEA_q.values.png", type="cairo",
#     width=11,height=5,units="in",res=500,pointsize = 12)
# tiff("/Users/Marta/Documents/WTCHG/DPhil/Data/input_for_GSEA/GSEA_q.values.tiff", type="cairo",
#      width=11,height=5,units="in",res=1200,pointsize = 12,compression = "lzw")

p6 <- ggplot(data = df_log, 
             aes(x = Type, 
                 y = `q value`, 
                 fill=Stage)) +
  
   geom_bar(stat="identity", color="black", position=position_dodge(0.9)) +
 
  #  scale_fill_brewer("Distance (kb)",type = "seq", palette = "YlOrBr", direction = 1)  +
  scale_fill_manual(name="Stage",values=diaPalette) +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05),linetype="dashed",size=0.8,col="grey") +
  scale_y_continuous(limits=c(0,3),
                     expand=c(0,0),
                     expression(-log[10]~"("~italic(q)~"-value)"))  + # set breaks and take out space under x axis, add sec axis
  xlab("") +
  #ylab("-log10 (p-value)") +
  theme( # no background grids
    #panel.border=element_rect(size=1),  # border as rectangle
    panel.border=element_blank(), # no border
    axis.line.x = element_line(size = 0.30, linetype = "solid", colour = "black"),# axis lines
    axis.line.y = element_line(size = 0.30, linetype = "solid", colour = "black"),
    panel.grid = element_blank(),  # no grid lines
    axis.text.x = element_text(size=16 ,colour="black"),
    axis.text.y =element_text(size=16 ,colour ="black"),  # bold ticks on y, leave margin for inward ticks
    axis.title.x=element_text(size=16 ,colour="black"), 
    axis.title.y=element_text(size=16 ,colour="black"), # no title on x or y
    legend.title = element_text(size=16 ,colour="black"),
    legend.text = element_text(size=15 ,colour="black"))
plot(p6)
# dev.off()

# a version with better alignment that needs gridExtra and gtable
AlignPlots <- function(...) {
  LegendWidth <- function(x) x$grobs[[8]]$grobs[[1]]$widths[[4]]
  
  plots.grobs <- lapply(list(...), ggplotGrob)
  
  max.widths <- do.call(unit.pmax, lapply(plots.grobs, "[[", "widths"))
  plots.grobs.eq.widths <- lapply(plots.grobs, function(x) {
    x$widths <- max.widths
    x
  })
  
  legends.widths <- lapply(plots.grobs, LegendWidth)
  max.legends.width <- do.call(max, legends.widths)
  plots.grobs.eq.widths.aligned <- lapply(plots.grobs.eq.widths, function(x) {
    if (is.gtable(x$grobs[[8]])) {
      x$grobs[[8]] <- gtable_add_cols(x$grobs[[8]],
                                      unit(abs(diff(c(LegendWidth(x),
                                                      max.legends.width))),
                                           "mm"))
    }
    x
  })
  
  plots.grobs.eq.widths.aligned
}


plots <- AlignPlots(p6+ggtitle("a") ,
                    p5+ggtitle("b") ,
                    p4+ggtitle("c"))
fig2=do.call(grid.arrange, plots)

#fig=multiplot( p3+ggtitle("a)"), p2+ggtitle("b)")+theme(legend.position = "none"),p1+ggtitle("c)"), cols=1)
#fig=multiplot( p3+ggtitle("a)"), p2+ggtitle("b)"),p1+ggtitle("c)"), cols=1)

tiff("/Users/Marta/Documents/WTCHG/DPhil/Islet\ diff\ paper\ Marta/fig2_barplots.tiff", type="cairo",
     width=9.5,height=9,units="in",res=600,pointsize = 12,compression = "none")
plot(fig2)
dev.off()

