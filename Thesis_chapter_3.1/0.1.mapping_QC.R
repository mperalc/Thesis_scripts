# Mapping QC
# 
# From mapping results Kundaje pipeline

library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(rjson)
library(reshape2)
library(stringr)

# Read QC file from pipeline
basedir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/Re-processed_Aug2018/QC/individual_qc/"



qclist <- list.files(basedir)
qclist <- qclist[substring(qclist, nchar(qclist), nchar(qclist))=='n']

names(qclist) <- qclist

json_summaries <- lapply(qclist, FUN=function(x) fromJSON(file= paste0(basedir,qclist[x])))

stat_getter <- function(i) {
  df <-
    data.frame(
      raw_bam_sequenced = json_summaries[[i]]$flagstat_qc[[1]]$total/2,
      read_mapping_rate = json_summaries[[i]]$flagstat_qc[[1]]$mapped_pct,
      raw_bam_total_paired_reads = json_summaries[[i]]$flagstat_qc[[1]]$paired_properly /
        2,
      pcnt_of_total_reads_paired = ((json_summaries[[i]]$flagstat_qc[[1]]$paired_properly /
                                      2) / (json_summaries[[i]]$flagstat_qc[[1]]$total/2))*100,
      filtered_bam_total_paired_reads = json_summaries[[i]]$dup_qc[[1]]$paired_reads,
      dupes_pcnt_in_filtered_bam_paired_reads = (json_summaries[[i]]$dup_qc[[1]]$dupes_pct)*100,
      filtered_mapped_dedup_total_read_pairs = json_summaries[[i]]$nodup_flagstat_qc[[1]]$read1,
      pcnt_total_raw_read_pairs_pass_all_qc = (json_summaries[[i]]$nodup_flagstat_qc[[1]]$read1 /
        (json_summaries[[i]]$flagstat_qc[[1]]$total / 2))*100,
      pcnt_of_paired_reads_pass_all_qc = (json_summaries[[i]]$nodup_flagstat_qc[[1]]$read1 /
        (json_summaries[[i]]$flagstat_qc[[1]]$paired_properly / 2))*100,
      pcnt_reads_in_peaks = (json_summaries[[i]]$frip_macs2_qc$rep1[[1]])*100,
      NRF = json_summaries[[i]]$pbc_qc[[1]]$NRF,
      PBC1 = json_summaries[[i]]$pbc_qc[[1]]$PBC1,
      PBC2 = json_summaries[[i]]$pbc_qc[[1]]$PBC2
    )
  rownames(df) <- i
  return(df)
}



full_stat_frame <- do.call('rbind', lapply(qclist, FUN=stat_getter))
full_stat_frame <- cbind(data.frame(sampid=substring(qclist, 1, nchar(qclist)-5), stringsAsFactors=FALSE), full_stat_frame)
write.table(full_stat_frame, row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE, file=paste(basedir, '/qc_summary.tsv', sep=''))



# Plot


png(
  paste(basedir, '/qc_pass_renamed_2.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 6,
  height = 6,
  pointsize = 12
)
ggplot(
  full_stat_frame,  aes(x = filtered_mapped_dedup_total_read_pairs, y = raw_bam_sequenced, label =    sampid)
) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  annotation_custom(
    grob = grobTree(textGrob(
      paste(
        "Pearson's r: ",
        round(
          cor(
            full_stat_frame$filtered_mapped_dedup_total_read_pairs,
            full_stat_frame$raw_bam_sequenced
          ),
          2)),
        x = 0.60,
        y = 0.20,
        hjust = 0,
        gp = gpar(
          col = "black",
          fontsize = 11,
          fontface = "bold"
        )
      )
    )) +
    geom_text_repel(aes(label = sampid), size = 2.5, col = 'black', point.padding = 0.5) +
      labs(
        x = 'Total read pairs passing all QC',
        
        y = 'Total raw sequenced read pairs'
      ) +
  geom_vline(xintercept = 25000000, color="grey",linetype="dashed") +
  
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))
dev.off()

# Raw mean pairs sequenced
mean(full_stat_frame$raw_bam_sequenced)
sd(full_stat_frame$raw_bam_sequenced)

# Raw mean pairs paired
mean(full_stat_frame$raw_bam_total_paired_reads)
sd(full_stat_frame$raw_bam_total_paired_reads)

# Percentage of raw total means paired
mean(full_stat_frame$pcnt_of_total_reads_paired)
sd(full_stat_frame$pcnt_of_total_reads_paired)
# Mean mapping rate
# 
mean(full_stat_frame$read_mapping_rate)

# Proportio of duplicates
mean(full_stat_frame$dupes_prop)
sd(full_stat_frame$dupes_prop)

# Filtered, deduped
mean(full_stat_frame$filtered_mapped_dedup_total_reads)
sd(full_stat_frame$filtered_mapped_dedup_total_reads)

# Proportion of properly paired reads that passed all QC
mean(full_stat_frame$prop_of_total_reads_pass_all_qc)
sd(full_stat_frame$prop_of_total_reads_pass_all_qc)

# Percentage of all raw reads that passed all QC
mean(full_stat_frame$pcnt_total_raw_read_pairs_pass_all_qc)
sd(full_stat_frame$pcnt_total_raw_read_pairs_pass_all_qc)

#Proportion of reads in peaks
mean(full_stat_frame$prop_reads_in_peaks)
sd(full_stat_frame$prop_reads_in_peaks)

# Reads that remain
for_boxplot = full_stat_frame[,c("read_mapping_rate","pcnt_of_total_reads_paired","pcnt_total_raw_read_pairs_pass_all_qc")]
colnames(for_boxplot) = c("Mapped reads", "Mapped, paired reads", "Mapped, paired, deduped reads")
melted = melt(for_boxplot)

png(
  paste(basedir, '/qc_reads_remain.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 4,
  height = 4
)
ggplot(melted, aes(x=variable, y=value, color = variable)) + 
  geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.2))  + 
  scale_color_manual(values = c("#177e89","#d8973c","#eb5160")) +
  labs(
    y = '% of total sequenced reads'
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 14)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(face="bold"),
        axis.title.y = element_text(face = "bold"),
        legend.position="none")

dev.off()

# Non-redundant fraction (NRF) scatterplot

highlight <- subset(full_stat_frame,NRF<0.9)

png(
  paste(basedir, '/library_complexity.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 6,
  height = 7
)
ggplot(full_stat_frame, aes(y=PBC2, x=PBC1)) + 
  geom_point(aes(size=NRF),alpha=1/3) +
  geom_point(data=highlight, aes(size=NRF), colour="orange") +
  geom_vline(xintercept = 0.9, color="grey",linetype="dashed")+
  labs(
    y = 'PBC2',
    x = "PBC1"
  ) +
  geom_text_repel(aes(label = sampid), size = 2.5, col = 'black', point.padding = 0.5) +
  
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(face="bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.position = "none")

dev.off()
