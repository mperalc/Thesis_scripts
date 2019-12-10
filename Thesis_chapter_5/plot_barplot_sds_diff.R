# plot barplot with means and SDs

library(ggplot2)
library(stringr)
library(gridExtra)


##
setwd("/Users/Marta/Documents/WTCHG/DPhil/Lab/qPCR/qPCR_results/2019-06-24_EXP90_KOLFC2_KOLFC2-KRAB_DIFF/")
datos=read.table("KOLFC2_DIFF_R.txt",stringsAsFactors = F, header = T, sep = "\t",skipNul = TRUE)  # read data


head(datos)
colnames(datos)=c("Stage","Gene","2^-DeltaCT","2^-DeltaCT+sd", "2^-DeltaCT-sd")    

datos$Stage=factor(datos$Stage,levels = unique(datos$Stage), ordered = T)
datos$Gene=factor(datos$Gene,levels = c("PROX1","PDX1","NKX6.1","NKX2.2","GCG","INS"), ordered = T)
datos$Gene
datos


p1 <- ggplot(data = datos, aes(x = Stage, y=`2^-DeltaCT`, fill=Gene)) +
  geom_bar(stat = "identity",
           width=0.9, 
           show.legend=T, 
           col="black", 
           size = .2,
           position=position_dodge()) +
  geom_errorbar(aes(
                    ymin=`2^-DeltaCT-sd`, 
                    ymax=`2^-DeltaCT+sd`,
                    group = Gene), 
                width=0.5,size=0.3,
                position = position_dodge(.9)) +
   scale_fill_brewer(palette = "Set1") +
    scale_y_continuous(expand = c(0,0 ), 
                      limits=c(0,1.8)) +
  theme_bw()+
 
  theme(
    strip.background = element_rect(fill = "gray90", colour = FALSE),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "gray90"),
    panel.margin = unit(1, "lines"),
    axis.text.x=element_text(size = 14,
                             #angle = 90, hjust = 1,vjust=0.5,
                             face = "bold"),
    axis.text.y=element_text(size = 14,
                             face = "bold"),
    axis.title =element_text(size = 16,
                             face = "bold"),
    legend.position = "none"
  ) 

plot(p1)

p2 = ggplot(data = datos, aes(x = Stage, y=`2^-DeltaCT`, fill=Gene)) +
  geom_bar(stat = "identity",
           width=0.9, 
           show.legend=T, 
           col="black", 
           size = .2,
           position=position_dodge()) +
  geom_errorbar(aes(
    ymin=`2^-DeltaCT-sd`, 
    ymax=`2^-DeltaCT+sd`,
    group = Gene), 
    width=0.5,size=0.3,
    position = position_dodge(.9)) +
  scale_fill_brewer(palette = "Set1") +
  
  scale_y_continuous(expand = c(0,0 ), 
                     limits=c(0,50)) +
  theme_bw()+
  
  theme(
    strip.background = element_rect(fill = "gray90", colour = FALSE),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "gray90"),
    panel.margin = unit(1, "lines"),
    axis.text.x=element_text(size = 14,
                             #angle = 90, hjust = 1,vjust=0.5,
                             face = "bold"),
    axis.text.y=element_text(size = 14,
                             face = "bold"),
    axis.title =element_text(size = 16,
                             face = "bold"),
    legend.position = "none"
  )  
  
 plot(p2)
 
 p3 = ggplot(data = datos, aes(x = Stage, y=`2^-DeltaCT`, fill=Gene)) +
   geom_bar(stat = "identity",
            width=0.9, 
            show.legend=T, 
            col="black", 
            size = .2,
            position=position_dodge()) +
   geom_errorbar(aes(
     ymin=`2^-DeltaCT-sd`, 
     ymax=`2^-DeltaCT+sd`,
     group = Gene), 
     width=0.5,size=0.3,
     position = position_dodge(.9)) +
   scale_fill_brewer(palette = "Set1") +
   
   scale_y_continuous(expand = c(0,0),limits = c(0,4000))+
   theme_bw()+
   
   theme(
     strip.background = element_rect(fill = "gray90", colour = FALSE),
     panel.grid.minor = element_blank(),
     panel.grid.major = element_blank(),
     panel.border = element_rect(colour = "gray90"),
     panel.margin = unit(1, "lines"),
     axis.text.x=element_text(size = 14,
                              #angle = 90, hjust = 1,vjust=0.5,
                              face = "bold"),
     axis.text.y=element_text(size = 14,
                              face = "bold"),
     axis.title =element_text(size = 16,
                              face = "bold"),
     legend.text = element_text(size = 14,
                                face = "bold"),
     
     legend.position="bottom"
   )  
 
 plot(p3) 
 
 #get common legend
 g_legend<-function(a.gplot){
   tmp <- ggplot_gtable(ggplot_build(a.gplot))
   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
   legend <- tmp$grobs[[leg]]
   return(legend)}
 mylegend<-g_legend(p3)
 
 # Arrange in grid
 
 png("KOLFC2_qPCR.png", type="cairo",
     width=8.5,height=7,units="in",res=300,pointsize = 12)
 
 
 grid.arrange( arrangeGrob(p1, p2,p3 + theme(legend.position = "none"),nrow=3), 
               nrow=2,mylegend,heights=c(10, 1))
dev.off()



