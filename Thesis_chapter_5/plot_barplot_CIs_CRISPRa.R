# plot barplot with means and CI



library(ggplot2)
library(stringr)


##
setwd("/Users/Marta/Documents/WTCHG/DPhil/Lab/qPCR/qPCR_results/qPCR_PROX1_CRISPRa_iPSC_SB_12052019/")
datos=read.table("forR.txt",stringsAsFactors = F, header = T, sep = "\t")  # read data
datos = datos[,c(1,3,4,5)]
## dummy bar to remove later so that aesthetics are mainained
datos = rbind((c("empty SAM 1",1,0,0)),datos)


datos
colnames(datos)=c("Sample","Relative expression", "plus SD","minus SD")    
datos$`Relative expression`=as.numeric(datos$`Relative expression`)
datos$`minus SD`=as.numeric(datos$`minus SD`)
datos$`plus SD`=as.numeric(datos$`plus SD`)
datos$Sample=factor(datos$Sample,levels = datos$Sample)
datos

datos$x = c("empty SAM","empty SAM",
            "Activation Prom SNP", "Activation Prom SNP",
            "Activation Prom 2",  "Activation Prom 2",
            "Activation Prom 3", "Activation Prom 3",
            "Activation Enh 1",   "Activation Enh 1",
            "Activation Enh 2", "Activation Enh 2",
            "Activation Enh 5","Activation Enh 5")
datos$rep = as.factor(c(1,2,1,2,1,2,1,2,1,2,1,2,1,2))
####


png("qPCR_CRISPRa.png", type="cairo",
    width=15,height=4.5,units="in",res=300,pointsize = 12)

p <- ggplot(data = datos, aes(x = x, y = `Relative expression`,group=Sample,col=rep,fill=rep)) +
  geom_bar(stat="identity", 
           width=0.6, 
           show.legend=T, 
           size = 1.0,
           position = position_dodge(0.8)) +
  scale_fill_manual(values = c( "1" = "orange", "2" = "blue"))+
  scale_color_manual(values = c( "1" = "orange", "2" = "blue"))+
  
  geom_errorbar(aes(ymin=`minus SD`, 
                    ymax=`plus SD`), 
                width=0.1,size=0.5 ,
                position = position_dodge(0.8),
                col="black") +
  xlab("Sample") + 
  scale_y_continuous(
    #breaks = seq(0,2, by = 0.2),
    expand = c(0,0 ), 
    limits=c(0,75),
    sec.axis=sec_axis(~.,breaks=derive(),labels=NULL)) + # set breaks and take out space under x axis, add sec axis
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + ## wrap labels
  
  theme_bw()+
  
  theme(
    strip.background = element_rect(fill="gray90", colour=FALSE),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour="gray90"),
    panel.margin = unit(1, "lines"),
    plot.title=element_text(vjust=1, size = 14, face = "bold"),
    legend.position="none",
    axis.title.x=element_text(vjust=-0.5, size = 15, face = "bold"),
    axis.title.y=element_text(vjust=1, size = 15, face = "bold"),
    axis.text.x = element_text(size=13,colour="black"),
    axis.text.y = element_text(size=13,colour="black"))+  
  geom_segment(mapping=aes(x=0,xend=nrow(datos)+1,y=1,yend=1),linetype="dashed",
               col="#444444",size=0.6,show.legend=F) # line relative expression control

plot(p)
dev.off()


## Are prom act 3 reps significantly different to empty sam reps?
# by non-parametric test
stats = read.table("stats_deltaCt.txt",stringsAsFactors = F, header = T, sep = "\t")  # read data
stats

wilcox.test(x = stats[stats$s %in% c("empty_SAM 1","empty_SAM 2"),"deltaCT"],
            y = stats[stats$s %in% c("PROX1 Act Prom 3.1","PROX1 Act Prom 3.2"),"deltaCT"])
t.test(x = stats[stats$s %in% c("empty_SAM 1","empty_SAM 2"),"deltaCT"],
       y = stats[stats$s %in% c("PROX1 Act Prom 3.1","PROX1 Act Prom 3.2"),"deltaCT"])

t.test(x = stats[stats$s %in% c("empty_SAM 1","empty_SAM 2"),"deltaCT"],
       y = stats[stats$s %in% c("PROX1 Act Prom 2.1","PROX1 Act Prom 2.2"),"deltaCT"])
