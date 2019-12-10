# Mapping QC H3K27ac
# 
# From mapping results encode pipeline

library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(rjson)
library(reshape2)
library(stringr)
library(tidyverse)

# Read QC file from pipeline
basedir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/ChIP-seq/ENCODE_pipeline/ENCODE_QC_aggregate_files/"


qclist <- list.files(basedir)
qclist <- qclist[substring(qclist, nchar(qclist), nchar(qclist))=='n']

names(qclist) <- qclist

json_summaries <- lapply(qclist, FUN=function(x) fromJSON(file= paste0(basedir,qclist[x])))

stat_getter <- function(i) {
  if(length(json_summaries[[i]]$flagstat_qc)>2){
  df <-
    data.frame(
      raw_bam_sequenced_rep1 = json_summaries[[i]]$flagstat_qc[[1]]$total/2,
      raw_bam_sequenced_rep2 = json_summaries[[i]]$flagstat_qc[[2]]$total/2,
      raw_bam_sequenced_rep3 = json_summaries[[i]]$flagstat_qc[[3]]$total/2,
      
      read_mapping_rate_rep1 = json_summaries[[i]]$flagstat_qc[[1]]$mapped_pct,
      read_mapping_rate_rep2 = json_summaries[[i]]$flagstat_qc[[2]]$mapped_pct,
      read_mapping_rate_rep3 = json_summaries[[i]]$flagstat_qc[[3]]$mapped_pct,
      
      raw_bam_total_paired_reads_rep1 = json_summaries[[i]]$flagstat_qc[[1]]$paired_properly /
        2,
      raw_bam_total_paired_reads_rep2 = json_summaries[[i]]$flagstat_qc[[2]]$paired_properly /
        2,
      raw_bam_total_paired_reads_rep3 = json_summaries[[i]]$flagstat_qc[[3]]$paired_properly /
        2,
      
      pcnt_of_total_reads_paired_rep1 = ((json_summaries[[i]]$flagstat_qc[[1]]$paired_properly /
                                            2) / (json_summaries[[i]]$flagstat_qc[[1]]$total/2))*100,
      pcnt_of_total_reads_paired_rep2 = ((json_summaries[[i]]$flagstat_qc[[2]]$paired_properly /
                                            2) / (json_summaries[[i]]$flagstat_qc[[2]]$total/2))*100,
      pcnt_of_total_reads_paired_rep3 = ((json_summaries[[i]]$flagstat_qc[[3]]$paired_properly /
                                            2) / (json_summaries[[i]]$flagstat_qc[[3]]$total/2))*100,
      
      filtered_bam_total_paired_reads_rep1 = json_summaries[[i]]$dup_qc[[1]]$paired_reads,
      filtered_bam_total_paired_reads_rep2 = json_summaries[[i]]$dup_qc[[2]]$paired_reads,
      filtered_bam_total_paired_reads_rep3 = json_summaries[[i]]$dup_qc[[3]]$paired_reads,
      
      dupes_pcnt_in_filtered_bam_paired_reads_rep1 = (json_summaries[[i]]$dup_qc[[1]]$dupes_pct)*100,
      dupes_pcnt_in_filtered_bam_paired_reads_rep2 = (json_summaries[[i]]$dup_qc[[2]]$dupes_pct)*100,
      dupes_pcnt_in_filtered_bam_paired_reads_rep3 = (json_summaries[[i]]$dup_qc[[3]]$dupes_pct)*100,
      
      filtered_mapped_dedup_total_read_pairs_rep1 = json_summaries[[i]]$nodup_flagstat_qc[[1]]$read1,
      filtered_mapped_dedup_total_read_pairs_rep2 = json_summaries[[i]]$nodup_flagstat_qc[[2]]$read1,
      filtered_mapped_dedup_total_read_pairs_rep3 = json_summaries[[i]]$nodup_flagstat_qc[[3]]$read1,
      
      pcnt_total_raw_read_pairs_pass_all_qc_rep1 = (json_summaries[[i]]$nodup_flagstat_qc[[1]]$read1 /
                                                 (json_summaries[[i]]$flagstat_qc[[1]]$total / 2))*100,
      pcnt_total_raw_read_pairs_pass_all_qc_rep2 = (json_summaries[[i]]$nodup_flagstat_qc[[2]]$read1 /
                                                      (json_summaries[[i]]$flagstat_qc[[2]]$total / 2))*100,
      pcnt_total_raw_read_pairs_pass_all_qc_rep3 = (json_summaries[[i]]$nodup_flagstat_qc[[3]]$read1 /
                                                      (json_summaries[[i]]$flagstat_qc[[3]]$total / 2))*100,
      
      pcnt_of_paired_reads_pass_all_qc_rep1 = (json_summaries[[i]]$nodup_flagstat_qc[[1]]$read1 /
                                            (json_summaries[[i]]$flagstat_qc[[1]]$paired_properly / 2))*100,
      pcnt_of_paired_reads_pass_all_qc_rep2 = (json_summaries[[i]]$nodup_flagstat_qc[[2]]$read1 /
                                                 (json_summaries[[i]]$flagstat_qc[[2]]$paired_properly / 2))*100,
      pcnt_of_paired_reads_pass_all_qc_rep3 = (json_summaries[[i]]$nodup_flagstat_qc[[3]]$read1 /
                                                 (json_summaries[[i]]$flagstat_qc[[3]]$paired_properly / 2))*100,
      
      pcnt_reads_in_peaks_rep1 = (json_summaries[[i]]$frip_macs2_qc$rep1[[1]])*100,
      pcnt_reads_in_peaks_rep2 = (json_summaries[[i]]$frip_macs2_qc$rep2[[1]])*100,
      pcnt_reads_in_peaks_rep3 = (json_summaries[[i]]$frip_macs2_qc$rep3[[1]])*100,
      
      NRF_rep1 = json_summaries[[i]]$pbc_qc[[1]]$NRF,
      NRF_rep2 = json_summaries[[i]]$pbc_qc[[2]]$NRF,
      NRF_rep3 = json_summaries[[i]]$pbc_qc[[3]]$NRF,
      
      PBC1_rep1 = json_summaries[[i]]$pbc_qc[[1]]$PBC1,
      PBC1_rep2 = json_summaries[[i]]$pbc_qc[[2]]$PBC1,
      PBC1_rep3 = json_summaries[[i]]$pbc_qc[[3]]$PBC1,
      
      PBC2_rep1 = json_summaries[[i]]$pbc_qc[[1]]$PBC2,
      PBC2_rep2 = json_summaries[[i]]$pbc_qc[[2]]$PBC2,
      PBC2_rep3 = json_summaries[[i]]$pbc_qc[[3]]$PBC2,
      
      ## controls  ###
      raw_bam_sequenced_ctrl_rep1= json_summaries[[i]]$ctl_flagstat_qc[[1]]$total/2,
      raw_bam_sequenced_ctrl_rep2 = json_summaries[[i]]$ctl_flagstat_qc[[2]]$total/2,
      raw_bam_sequenced_ctrl_rep3 = json_summaries[[i]]$ctl_flagstat_qc[[3]]$total/2,
      
      read_mapping_rate_ctrl_rep1 = json_summaries[[i]]$ctl_flagstat_qc[[1]]$mapped_pct,
      read_mapping_rate_ctrl_rep2= json_summaries[[i]]$ctl_flagstat_qc[[2]]$mapped_pct,
      read_mapping_rate_ctrl_rep3 = json_summaries[[i]]$ctl_flagstat_qc[[3]]$mapped_pct,
      
      raw_bam_total_paired_reads_ctrl_rep1 = json_summaries[[i]]$ctl_flagstat_qc[[1]]$paired_properly /
        2,
      raw_bam_total_paired_reads_ctrl_rep2 = json_summaries[[i]]$ctl_flagstat_qc[[2]]$paired_properly /
        2,
      raw_bam_total_paired_reads_ctrl_rep3 = json_summaries[[i]]$ctl_flagstat_qc[[3]]$paired_properly /
        2,
      
      pcnt_of_total_reads_paired_ctrl_rep1 = ((json_summaries[[i]]$ctl_flagstat_qc[[1]]$paired_properly /
                                                 2) / (json_summaries[[i]]$ctl_flagstat_qc[[1]]$total/2))*100,
      pcnt_of_total_reads_paired_ctrl_rep2 = ((json_summaries[[i]]$ctl_flagstat_qc[[2]]$paired_properly /
                                                 2) / (json_summaries[[i]]$ctl_flagstat_qc[[2]]$total/2))*100,
      pcnt_of_total_reads_paired_ctrl_rep3 = ((json_summaries[[i]]$ctl_flagstat_qc[[3]]$paired_properly /
                                                 2) / (json_summaries[[i]]$ctl_flagstat_qc[[3]]$total/2))*100,
      
      filtered_bam_total_paired_reads_ctrl_rep1 = json_summaries[[i]]$ctl_dup_qc[[1]]$paired_reads,
      filtered_bam_total_paired_reads_ctrl_rep2 = json_summaries[[i]]$ctl_dup_qc[[2]]$paired_reads,
      filtered_bam_total_paired_reads_ctrl_rep3 = json_summaries[[i]]$ctl_dup_qc[[3]]$paired_reads,
      
      dupes_pcnt_in_filtered_bam_paired_reads_ctrl_rep1 = (json_summaries[[i]]$ctl_dup_qc[[1]]$dupes_pct)*100,
      dupes_pcnt_in_filtered_bam_paired_reads_ctrl_rep2 = (json_summaries[[i]]$ctl_dup_qc[[2]]$dupes_pct)*100,
      dupes_pcnt_in_filtered_bam_paired_reads_ctrl_rep3 = (json_summaries[[i]]$ctl_dup_qc[[3]]$dupes_pct)*100,
      
      filtered_mapped_dedup_total_read_pairs_ctrl_rep1 = json_summaries[[i]]$ctl_nodup_flagstat_qc[[1]]$read1,
      filtered_mapped_dedup_total_read_pairs_ctrl_rep2 = json_summaries[[i]]$ctl_nodup_flagstat_qc[[2]]$read1,
      filtered_mapped_dedup_total_read_pairs_ctrl_rep3 = json_summaries[[i]]$ctl_nodup_flagstat_qc[[3]]$read1,
      
      pcnt_total_raw_read_pairs_pass_all_qc_ctrl_rep1 = (json_summaries[[i]]$ctl_nodup_flagstat_qc[[1]]$read1 /
                                                           (json_summaries[[i]]$ctl_flagstat_qc[[1]]$total / 2))*100,
      pcnt_total_raw_read_pairs_pass_all_qc_ctrl_rep2 = (json_summaries[[i]]$ctl_nodup_flagstat_qc[[2]]$read1 /
                                                           (json_summaries[[i]]$ctl_flagstat_qc[[2]]$total / 2))*100,
      pcnt_total_raw_read_pairs_pass_all_qc_ctrl_rep3 = (json_summaries[[i]]$ctl_nodup_flagstat_qc[[3]]$read1 /
                                                           (json_summaries[[i]]$ctl_flagstat_qc[[3]]$total / 2))*100,
      
      pcnt_of_paired_reads_pass_all_qc_ctrl_rep1 = (json_summaries[[i]]$ctl_nodup_flagstat_qc[[1]]$read1 /
                                                      (json_summaries[[i]]$ctl_flagstat_qc[[1]]$paired_properly / 2))*100,
      pcnt_of_paired_reads_pass_all_qc_ctrl_rep2 = (json_summaries[[i]]$ctl_nodup_flagstat_qc[[2]]$read1 /
                                                      (json_summaries[[i]]$ctl_flagstat_qc[[2]]$paired_properly / 2))*100,
      pcnt_of_paired_reads_pass_all_qc_ctrl_rep3 = (json_summaries[[i]]$ctl_nodup_flagstat_qc[[3]]$read1 /
                                                      (json_summaries[[i]]$ctl_flagstat_qc[[3]]$paired_properly / 2))*100,
      
      NRF_ctrl_rep1 = json_summaries[[i]]$ctl_pbc_qc[[1]]$NRF,
      NRF_ctrl_rep2 = json_summaries[[i]]$ctl_pbc_qc[[2]]$NRF,
      NRF_ctrl_rep3= json_summaries[[i]]$ctl_pbc_qc[[3]]$NRF,
      
      PBC1_ctrl_rep1= json_summaries[[i]]$ctl_pbc_qc[[1]]$PBC1,
      PBC1_ctrl_rep2 = json_summaries[[i]]$ctl_pbc_qc[[2]]$PBC1,
      PBC1_ctrl_rep3 = json_summaries[[i]]$ctl_pbc_qc[[3]]$PBC1,
      
      PBC2_ctrl_rep1 = json_summaries[[i]]$ctl_pbc_qc[[1]]$PBC2,
      PBC2_ctrl_rep2 = json_summaries[[i]]$ctl_pbc_qc[[2]]$PBC2,
      PBC2_ctrl_rep3 = json_summaries[[i]]$ctl_pbc_qc[[3]]$PBC2
      
    )
  }
  else{
    df <-
      data.frame(
        raw_bam_sequenced_rep1 = json_summaries[[i]]$flagstat_qc[[1]]$total/2,
        raw_bam_sequenced_rep2 = json_summaries[[i]]$flagstat_qc[[2]]$total/2,
        
        read_mapping_rate_rep1 = json_summaries[[i]]$flagstat_qc[[1]]$mapped_pct,
        read_mapping_rate_rep2 = json_summaries[[i]]$flagstat_qc[[2]]$mapped_pct,
        
        raw_bam_total_paired_reads_rep1 = json_summaries[[i]]$flagstat_qc[[1]]$paired_properly /  2,
        raw_bam_total_paired_reads_rep2 = json_summaries[[i]]$flagstat_qc[[2]]$paired_properly /  2,
        
        pcnt_of_total_reads_paired_rep1 = ((json_summaries[[i]]$flagstat_qc[[1]]$paired_properly /
                                         2) / (json_summaries[[i]]$flagstat_qc[[1]]$total/2))*100,
        pcnt_of_total_reads_paired_rep2 = ((json_summaries[[i]]$flagstat_qc[[2]]$paired_properly /
                                              2) / (json_summaries[[i]]$flagstat_qc[[2]]$total/2))*100,
        
        filtered_bam_total_paired_reads_rep1 = json_summaries[[i]]$dup_qc[[1]]$paired_reads,
        filtered_bam_total_paired_reads_rep2 = json_summaries[[i]]$dup_qc[[2]]$paired_reads,
        
        dupes_pcnt_in_filtered_bam_paired_reads_rep1 = (json_summaries[[i]]$dup_qc[[1]]$dupes_pct)*100,
        dupes_pcnt_in_filtered_bam_paired_reads_rep2 = (json_summaries[[i]]$dup_qc[[2]]$dupes_pct)*100,
        
        filtered_mapped_dedup_total_read_pairs_rep1 = json_summaries[[i]]$nodup_flagstat_qc[[1]]$read1,
        filtered_mapped_dedup_total_read_pairs_rep2 = json_summaries[[i]]$nodup_flagstat_qc[[2]]$read1,
        
        pcnt_total_raw_read_pairs_pass_all_qc_rep1 = (json_summaries[[i]]$nodup_flagstat_qc[[1]]$read1 /
                                                        (json_summaries[[i]]$flagstat_qc[[1]]$total / 2))*100,
        pcnt_total_raw_read_pairs_pass_all_qc_rep2 = (json_summaries[[i]]$nodup_flagstat_qc[[2]]$read1 /
                                                        (json_summaries[[i]]$flagstat_qc[[2]]$total / 2))*100,
        
        pcnt_of_paired_reads_pass_all_qc_rep1 = (json_summaries[[i]]$nodup_flagstat_qc[[1]]$read1 /
                                                   (json_summaries[[i]]$flagstat_qc[[1]]$paired_properly / 2))*100,
        pcnt_of_paired_reads_pass_all_qc_rep2 = (json_summaries[[i]]$nodup_flagstat_qc[[2]]$read1 /
                                                   (json_summaries[[i]]$flagstat_qc[[2]]$paired_properly / 2))*100,
        
        pcnt_reads_in_peaks_rep1 = (json_summaries[[i]]$frip_macs2_qc$rep1[[1]])*100,
        pcnt_reads_in_peaks_rep2 = (json_summaries[[i]]$frip_macs2_qc$rep2[[1]])*100,
        
        NRF_rep1 = json_summaries[[i]]$pbc_qc[[1]]$NRF,
        NRF_rep2 = json_summaries[[i]]$pbc_qc[[2]]$NRF,
        
        PBC1_rep1 = json_summaries[[i]]$pbc_qc[[1]]$PBC1,
        PBC1_rep2 = json_summaries[[i]]$pbc_qc[[2]]$PBC1,
        
        PBC2_rep1 = json_summaries[[i]]$pbc_qc[[1]]$PBC2,
        PBC2_rep2 = json_summaries[[i]]$pbc_qc[[2]]$PBC2,
        
        ## controls  ###
        raw_bam_sequenced_ctrl_rep1= json_summaries[[i]]$ctl_flagstat_qc[[1]]$total/2,
        raw_bam_sequenced_ctrl_rep2 = json_summaries[[i]]$ctl_flagstat_qc[[2]]$total/2,
        raw_bam_sequenced_ctrl_rep3 = json_summaries[[i]]$ctl_flagstat_qc[[3]]$total/2,
        
        read_mapping_rate_ctrl_rep1 = json_summaries[[i]]$ctl_flagstat_qc[[1]]$mapped_pct,
        read_mapping_rate_ctrl_rep2= json_summaries[[i]]$ctl_flagstat_qc[[2]]$mapped_pct,
        read_mapping_rate_ctrl_rep3 = json_summaries[[i]]$ctl_flagstat_qc[[3]]$mapped_pct,
        
        raw_bam_total_paired_reads_ctrl_rep1 = json_summaries[[i]]$ctl_flagstat_qc[[1]]$paired_properly /
          2,
        raw_bam_total_paired_reads_ctrl_rep2 = json_summaries[[i]]$ctl_flagstat_qc[[2]]$paired_properly /
          2,
        raw_bam_total_paired_reads_ctrl_rep3 = json_summaries[[i]]$ctl_flagstat_qc[[3]]$paired_properly /
          2,
        
        pcnt_of_total_reads_paired_ctrl_rep1 = ((json_summaries[[i]]$ctl_flagstat_qc[[1]]$paired_properly /
                                                   2) / (json_summaries[[i]]$ctl_flagstat_qc[[1]]$total/2))*100,
        pcnt_of_total_reads_paired_ctrl_rep2 = ((json_summaries[[i]]$ctl_flagstat_qc[[2]]$paired_properly /
                                                   2) / (json_summaries[[i]]$ctl_flagstat_qc[[2]]$total/2))*100,
        pcnt_of_total_reads_paired_ctrl_rep3 = ((json_summaries[[i]]$ctl_flagstat_qc[[3]]$paired_properly /
                                                   2) / (json_summaries[[i]]$ctl_flagstat_qc[[3]]$total/2))*100,
        
        filtered_bam_total_paired_reads_ctrl_rep1 = json_summaries[[i]]$ctl_dup_qc[[1]]$paired_reads,
        filtered_bam_total_paired_reads_ctrl_rep2 = json_summaries[[i]]$ctl_dup_qc[[2]]$paired_reads,
        filtered_bam_total_paired_reads_ctrl_rep3 = json_summaries[[i]]$ctl_dup_qc[[3]]$paired_reads,
        
        dupes_pcnt_in_filtered_bam_paired_reads_ctrl_rep1 = (json_summaries[[i]]$ctl_dup_qc[[1]]$dupes_pct)*100,
        dupes_pcnt_in_filtered_bam_paired_reads_ctrl_rep2 = (json_summaries[[i]]$ctl_dup_qc[[2]]$dupes_pct)*100,
        dupes_pcnt_in_filtered_bam_paired_reads_ctrl_rep3 = (json_summaries[[i]]$ctl_dup_qc[[3]]$dupes_pct)*100,
        
        filtered_mapped_dedup_total_read_pairs_ctrl_rep1 = json_summaries[[i]]$ctl_nodup_flagstat_qc[[1]]$read1,
        filtered_mapped_dedup_total_read_pairs_ctrl_rep2 = json_summaries[[i]]$ctl_nodup_flagstat_qc[[2]]$read1,
        filtered_mapped_dedup_total_read_pairs_ctrl_rep3 = json_summaries[[i]]$ctl_nodup_flagstat_qc[[3]]$read1,
        
        pcnt_total_raw_read_pairs_pass_all_qc_ctrl_rep1 = (json_summaries[[i]]$ctl_nodup_flagstat_qc[[1]]$read1 /
                                                             (json_summaries[[i]]$ctl_flagstat_qc[[1]]$total / 2))*100,
        pcnt_total_raw_read_pairs_pass_all_qc_ctrl_rep2 = (json_summaries[[i]]$ctl_nodup_flagstat_qc[[2]]$read1 /
                                                             (json_summaries[[i]]$ctl_flagstat_qc[[2]]$total / 2))*100,
        pcnt_total_raw_read_pairs_pass_all_qc_ctrl_rep3 = (json_summaries[[i]]$ctl_nodup_flagstat_qc[[3]]$read1 /
                                                             (json_summaries[[i]]$ctl_flagstat_qc[[3]]$total / 2))*100,
        
        pcnt_of_paired_reads_pass_all_qc_ctrl_rep1 = (json_summaries[[i]]$ctl_nodup_flagstat_qc[[1]]$read1 /
                                                        (json_summaries[[i]]$ctl_flagstat_qc[[1]]$paired_properly / 2))*100,
        pcnt_of_paired_reads_pass_all_qc_ctrl_rep2 = (json_summaries[[i]]$ctl_nodup_flagstat_qc[[2]]$read1 /
                                                        (json_summaries[[i]]$ctl_flagstat_qc[[2]]$paired_properly / 2))*100,
        pcnt_of_paired_reads_pass_all_qc_ctrl_rep3 = (json_summaries[[i]]$ctl_nodup_flagstat_qc[[3]]$read1 /
                                                        (json_summaries[[i]]$ctl_flagstat_qc[[3]]$paired_properly / 2))*100,
        
        NRF_ctrl_rep1 = json_summaries[[i]]$ctl_pbc_qc[[1]]$NRF,
        NRF_ctrl_rep2 = json_summaries[[i]]$ctl_pbc_qc[[2]]$NRF,
        NRF_ctrl_rep3= json_summaries[[i]]$ctl_pbc_qc[[3]]$NRF,
        
        PBC1_ctrl_rep1= json_summaries[[i]]$ctl_pbc_qc[[1]]$PBC1,
        PBC1_ctrl_rep2 = json_summaries[[i]]$ctl_pbc_qc[[2]]$PBC1,
        PBC1_ctrl_rep3 = json_summaries[[i]]$ctl_pbc_qc[[3]]$PBC1,
        
        PBC2_ctrl_rep1 = json_summaries[[i]]$ctl_pbc_qc[[1]]$PBC2,
        PBC2_ctrl_rep2 = json_summaries[[i]]$ctl_pbc_qc[[2]]$PBC2,
        PBC2_ctrl_rep3 = json_summaries[[i]]$ctl_pbc_qc[[3]]$PBC2
      )
  }
  rownames(df) <- i
  return(df)
}



full_stat_frame <- lapply(qclist, FUN=stat_getter)
stat_frame = as.data.frame(names(unlist(full_stat_frame)) %>% str_detect('raw_bam_sequenced') %>% keep(unlist(full_stat_frame),.))
colnames(stat_frame) = "raw_bam_sequenced"
stat_frame$rep = sub('.*\\_', '', rownames(stat_frame))
stat_frame$stage = sub('\\..*', '', rownames(stat_frame))
stat_frame[stat_frame$rep == "rep1","rep"] = "SBAd2.1"
stat_frame[stat_frame$rep == "rep2","rep"] = "SBAd3.1"
stat_frame[stat_frame$rep == "rep3","rep"] = "SBNeo1.1"
stat_frame$Type = ifelse(grepl("ctrl", rownames(stat_frame)),"input","sample")
stat_frame$sampid = paste(stat_frame$Type,stat_frame$stage,stat_frame$rep,sep = "_")
stat_frame$filtered_mapped_dedup_total_read_pairs = as.numeric(names(unlist(full_stat_frame)) %>% 
                                                                 str_detect('filtered_mapped_dedup_total_read_pairs') %>% 
                                                                 keep(unlist(full_stat_frame),.))
stat_frame$raw_bam_total_paired_reads = as.numeric(names(unlist(full_stat_frame)) %>% 
                                                                 str_detect('raw_bam_total_paired_reads') %>% 
                                                                 keep(unlist(full_stat_frame),.))
  
stat_frame$pcnt_of_total_reads_paired  = as.numeric(names(unlist(full_stat_frame)) %>% 
                                                  str_detect('pcnt_of_total_reads_paired') %>% 
                                                  keep(unlist(full_stat_frame),.))

stat_frame$read_mapping_rate = as.numeric(names(unlist(full_stat_frame)) %>% 
                                           str_detect('read_mapping_rate') %>% 
                                           keep(unlist(full_stat_frame),.))

stat_frame$dupes_pcnt = as.numeric(names(unlist(full_stat_frame)) %>% 
                                     str_detect('dupes_pcnt') %>% 
                                     keep(unlist(full_stat_frame),.))

stat_frame$filtered_mapped_dedup_total_reads = as.numeric(names(unlist(full_stat_frame)) %>% 
                                                            str_detect('filtered_mapped_dedup_total_read_pairs') %>% 
                                                            keep(unlist(full_stat_frame),.))

stat_frame$pcnt_of_paired_reads_pass_all_qc = as.numeric(names(unlist(full_stat_frame)) %>% 
                                                          str_detect('pcnt_of_paired_reads_pass_all_qc') %>% 
                                                          keep(unlist(full_stat_frame),.))

stat_frame$pcnt_total_raw_read_pairs_pass_all_qc = as.numeric(names(unlist(full_stat_frame)) %>% 
                                                                str_detect('pcnt_total_raw_read_pairs_pass_all_qc') %>% 
                                                                keep(unlist(full_stat_frame),.))

stat_frame$NRF = as.numeric(names(unlist(full_stat_frame)) %>% 
                          str_detect('NRF') %>% 
                          keep(unlist(full_stat_frame),.))

stat_frame$PBC1= as.numeric(names(unlist(full_stat_frame)) %>% 
                              str_detect('PBC1') %>% 
                              keep(unlist(full_stat_frame),.))

stat_frame$PBC2= as.numeric(names(unlist(full_stat_frame)) %>% 
                              str_detect('PBC2') %>% 
                              keep(unlist(full_stat_frame),.))


pcnt_reads_in_peaks = names(unlist(full_stat_frame)) %>% 
                                              str_detect('pcnt_reads_in_peaks') %>% 
                                              keep(unlist(full_stat_frame),.)

write.table(as.data.frame(pcnt_reads_in_peaks),file = paste0(basedir,"FRIP.tsv"),row.names = T,quote = F)
write.table(stat_frame,file = paste0(basedir,"stat_frame.tsv"),row.names = T,quote = F)


png(
  paste(basedir, '/qc_pass_renamed_2.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 10,
  height = 10,
  pointsize = 12
)
ggplot(
  stat_frame,  aes(x = filtered_mapped_dedup_total_read_pairs, y = raw_bam_sequenced, label = sampid),
) +
  geom_point(stat = "identity",aes(color = Type)) +
  geom_smooth(method = "lm", se = F) +
  annotation_custom(
    grob = grobTree(textGrob(
      paste(
        "Pearson's r: ",
        round(
          cor(
            stat_frame$filtered_mapped_dedup_total_read_pairs,
            stat_frame$raw_bam_sequenced
          ),
          2)),
      x = 0.60,
      y = 0.80,
      hjust = 0,
      gp = gpar(
        col = "black",
        fontsize = 11,
        fontface = "bold"
      )
    )
    )) +
  geom_text_repel(aes(label = sampid), size = 3, col = 'black', point.padding = 0.5) +
  labs(
    x = 'Total read pairs passing all QC',
    
    y = 'Total raw sequenced read pairs'
  ) +
  geom_vline(xintercept = 20000000, color="grey",linetype="dashed") +
  
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14),
        axis.text = element_text(size = 14))
dev.off()



# Raw mean pairs sequenced
mean(stat_frame$raw_bam_sequenced)
sd(stat_frame$raw_bam_sequenced)

# Raw mean pairs paired
mean(stat_frame$raw_bam_total_paired_reads)
sd(stat_frame$raw_bam_total_paired_reads)

# Percentage of raw total reads paired
mean(stat_frame$pcnt_of_total_reads_paired)
sd(stat_frame$pcnt_of_total_reads_paired)
# Mean mapping rate
# 
mean(stat_frame$read_mapping_rate)

# Proportio of duplicates  ###
# sample
mean(stat_frame$dupes_pcnt[stat_frame$Type == "sample"])
sd(stat_frame$dupes_pcnt[stat_frame$Type == "sample"])
# input
mean(stat_frame$dupes_pcnt[stat_frame$Type == "input"])
sd(stat_frame$dupes_pcnt[stat_frame$Type == "input"])

# Filtered, deduped ###
mean(stat_frame$filtered_mapped_dedup_total_reads[stat_frame$Type == "sample"])
sd(stat_frame$filtered_mapped_dedup_total_reads[stat_frame$Type == "sample"])
mean(stat_frame$filtered_mapped_dedup_total_reads[stat_frame$Type == "input"])
sd(stat_frame$filtered_mapped_dedup_total_reads[stat_frame$Type == "input"])

# Proportion of properly paired reads that passed all QC
mean(stat_frame$prop_of_total_reads_pass_all_qc)
sd(stat_frame$prop_of_total_reads_pass_all_qc)

# Percentage of all raw reads that passed all QC
mean(stat_frame$pcnt_total_raw_read_pairs_pass_all_qc[stat_frame$Type == "sample"])
sd(stat_frame$pcnt_total_raw_read_pairs_pass_all_qc[stat_frame$Type == "sample"])
mean(stat_frame$pcnt_total_raw_read_pairs_pass_all_qc[stat_frame$Type == "input"])
sd(stat_frame$pcnt_total_raw_read_pairs_pass_all_qc[stat_frame$Type == "input"])

#Proportion of reads in peaks
mean(pcnt_reads_in_peaks)
sd(pcnt_reads_in_peaks)

# Reads that remain
for_boxplot = stat_frame[,c("read_mapping_rate","pcnt_of_total_reads_paired","pcnt_total_raw_read_pairs_pass_all_qc","Type")]
colnames(for_boxplot) = c("Mapped reads", "Mapped, paired reads", "Mapped, paired, deduped reads","Type")
melted = melt(for_boxplot)

png(
  paste(basedir, '/qc_reads_remain.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 5.5,
  height = 4.5
)
ggplot(melted, aes(x=variable, y=value,  fill = Type)) + 
  geom_boxplot() + 
  # geom_jitter(shape=16, position=position_jitter(0.2), aes(col = Type))  + 
  # scale_color_manual(values = c("#177e89","#eb5160")) +
  labs(
    y = '% of total sequenced reads'
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 14)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(face="bold"),
        axis.title.y = element_text(face = "bold"),
        legend.position="right")

dev.off()

# Non-redundant fraction (NRF) scatterplot

highlight <- subset(stat_frame,NRF<0.7)
## all samples have NRF < 0.9. 0.7 indicates concerning levels of library complexity

png(
  paste(basedir, '/library_complexity.png', sep = ''),units = "in",res = 400,type = "cairo",
  width = 8,
  height = 8
)
ggplot(stat_frame, aes(y=PBC2, x=PBC1)) + 
  geom_point(aes(size=NRF),alpha=1/3, color = "grey") +
  geom_point(data=highlight, aes(size=NRF),alpha=1/3, colour="darkorange") +
  geom_vline(xintercept = 0.9, color="grey",linetype="dashed")+
  labs(
    y = 'PBC2',
    x = "PBC1"
  ) +
  geom_text_repel(aes(label = sampid), size = 3, col = 'black', point.padding = 0.5) +
  
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(face="bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14),
        axis.title.x = element_text(face = "bold",size = 14),
        legend.position = "none")

dev.off()


