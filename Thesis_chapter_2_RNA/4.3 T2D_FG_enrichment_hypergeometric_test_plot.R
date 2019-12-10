# Diabetologia Figure 2 version with dots instead of bars

# plot for thesis

library(ggplot2)
library(reshape2)  # to modify dataframes for ggplot2
library(ggsignif)  # significance lines
library(grid)  # tick marks on the inside
library(RColorBrewer)
library(gridExtra)
library(gtable)
library(stringr)  # to wrap text of labels
##
setwd("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/Genes_in_credible_regions_T2D_HRC_2018")
datos=read.table("hypergeom_0kb_for_thesis.txt",stringsAsFactors = F,header = T,sep="\t")  # read data

head(datos)
colnames(datos)=c("Stage","Trait","GWAS loci", "Permuted p-value","Low CI","High CI", "Distance")    
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
df_log$Stage=factor(df_log$Stage,levels = unique(df_log$Stage),ordered = T)

diaPalette <- c("#000000","#CADAE8", "#7883BA", "#755A91", "#CC85B1",   "#F4B8B0",
                "#96665A", "#96165A")  # Diabetologia palette


p1 <- ggplot(data = df_log, 
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
plot(p1)


## Now the plot of the GSEA results

setwd("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/GSEA/HRC_T2D_data_2018")
datos=read.table("GSEA_for_thesis.txt",stringsAsFactors = F,header = T)  # read data


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



p2 <- ggplot(data = df_log, 
             aes(x = Type, 
                 y = `q value`, 
                 fill=Stage)) +
  
  geom_point(stat="identity",
             position=position_dodge(width = 0.7),
             size=6,
             pch=21) +  # this + fill gives coloured fill and black outline+ error bars
  
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
plot(p2)
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


plots <- AlignPlots(p2+ggtitle("a") ,
                    p1+ggtitle("b"))
fig2=do.call(grid.arrange, plots)

#fig=multiplot( p3+ggtitle("a)"), p2+ggtitle("b)")+theme(legend.position = "none"),p1+ggtitle("c)"), cols=1)
#fig=multiplot( p3+ggtitle("a)"), p2+ggtitle("b)"),p1+ggtitle("c)"), cols=1)

png("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/fig2.8_dots.png", type="cairo",
     width=9.5,height=9,units="in",res=600,pointsize = 12)
plot(fig2)
dev.off()

## fig 2.9 A fasting glucose distance hypergeometric
setwd("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/hypergeom_FG")
FG=read.table("FG_for_thesis_distance_BLC.txt",stringsAsFactors = F,header = T,sep="\t")  # read data
datos = FG
colnames(datos)=c("Stage","Trait","GWAS loci", "Permuted p-value","Low CI","High CI", "Distance")    
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
df_log$Stage=factor(df_log$Stage,levels = unique(df_log$Stage),ordered = T)

pal1 <- c("#96165A","#e85d75","#e16f7c","#e1aa7d","#f2d0a4")  #  palette


p3 <- ggplot(data = df_log, 
             aes(x = bin, 
                 y = `Permuted.p.value`, 
                 fill=Distance)) +
  
  # geom_bar(stat="identity", color="black", position=position_dodge()) +
  # geom_errorbar(aes(ymin=`Low.CI`, ymax=`High.CI`),
  #               # color="black",
  #               position =position_dodge(0.7),
  #               width=0.2,size=0.5) +
  geom_point(stat="identity",aes(fill=Distance),
             position=position_dodge(width = 0.7),
             size=6,
             pch=21) +  # this + fill gives coloured fill and black outline+ error bars
  #scale_fill_brewer("Stage",type = "seq", palette = "YlOrBr", direction = 1)  +
  scale_fill_manual(values=pal1) +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05),linetype="dashed",size=0.8,col="grey") +
  scale_y_continuous(limits=c(0,5),
                     expand=c(0,0),
                     expression(-log[10]~"permuted ("~italic(p)~"-value)"))  + # set breaks and take out space under x axis, add sec axis
  xlab("Fasting glucose: BLC stage") +
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

png("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/fig2.9_FG.png", type="cairo",
    width=4,height=4,units="in",res=600,pointsize = 12)
plot(p3)
dev.off()


### fig 2.9 B T2D distance hypergeometric beta cell function loci and all
setwd("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/Genes_in_credible_regions_T2D_HRC_2018")
datos=read.table("T2D_for_thesis_distance_BLC_beta_loci.txt",stringsAsFactors = F,header = T,sep="\t")  # read data
colnames(datos)=c("Test","Trait","GWAS loci", "Permuted p-value","Low CI","High CI", "Distance","Stage")    
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
df_log$Stage=factor(df_log$Stage,levels = unique(df_log$Stage),ordered = T)

pal1 <- c("#96165A","#e85d75","#e16f7c","#e1aa7d","#f2d0a4")  #  palette


p4 <- ggplot(data = df_log, 
             aes(x = GWAS.loci, 
                 y = `Permuted.p.value`, 
                 fill=Distance)) +
  
  # geom_bar(stat="identity", color="black", position=position_dodge()) +
  # geom_errorbar(aes(ymin=`Low.CI`, ymax=`High.CI`),
  #               # color="black",
  #               position =position_dodge(0.7),
  #               width=0.2,size=0.5) +
  geom_point(stat="identity",aes(fill=Distance),
             position=position_dodge(width = 0.7),
             size=6,
             pch=21) +  # this + fill gives coloured fill and black outline+ error bars
  #scale_fill_brewer("Stage",type = "seq", palette = "YlOrBr", direction = 1)  +
  scale_fill_manual(values=pal1) +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05),linetype="dashed",size=0.8,col="grey") +
  scale_y_continuous(limits=c(0,5),
                     expand=c(0,0),
                     expression(-log[10]~"permuted ("~italic(p)~"-value)"))  + # set breaks and take out space under x axis, add sec axis
  xlab("T2D (beta-cell): BLC stage") +
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
plot(p4)

png("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/fig2.9_T2D_beta_cell_function.png", type="cairo",
    width=4,height=4,units="in",res=600,pointsize = 12)
plot(p4)
dev.off()

### All T2D genes
setwd("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/Genes_in_credible_regions_T2D_HRC_2018")
datos=read.table("T2D_for_thesis_distance_EN.txt",stringsAsFactors = F,header = T,sep="\t")  # read data
colnames(datos)=c("Test","Trait","GWAS loci", "Permuted p-value","Low CI","High CI", "Distance","Stage")    
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
df_log$Stage=factor(df_log$Stage,levels = unique(df_log$Stage),ordered = T)

pal1 <- c("#96165A","#e85d75","#e16f7c","#e1aa7d","#f2d0a4")  #  palette


p5 <- ggplot(data = df_log, 
             aes(x = GWAS.loci, 
                 y = `Permuted.p.value`, 
                 fill=Distance)) +
  
  # geom_bar(stat="identity", color="black", position=position_dodge()) +
  # geom_errorbar(aes(ymin=`Low.CI`, ymax=`High.CI`),
  #               # color="black",
  #               position =position_dodge(0.7),
  #               width=0.2,size=0.5) +
  geom_point(stat="identity",aes(fill=Distance),
             position=position_dodge(width = 0.7),
             size=6,
             pch=21) +  # this + fill gives coloured fill and black outline+ error bars
  #scale_fill_brewer("Stage",type = "seq", palette = "YlOrBr", direction = 1)  +
  scale_fill_manual(values=pal1) +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05),linetype="dashed",size=0.8,col="grey") +
  scale_y_continuous(limits=c(0,5),
                     expand=c(0,0),
                     expression(-log[10]~"permuted ("~italic(p)~"-value)"))  + # set breaks and take out space under x axis, add sec axis
  xlab("T2D (all): EN stage") +
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

png("/Users/Marta/Documents/WTCHG/R scripts/Diff_v2/Thesis/fig2.9_T2D_all.png", type="cairo",
    width=4,height=4,units="in",res=600,pointsize = 12)
plot(p5)
dev.off()