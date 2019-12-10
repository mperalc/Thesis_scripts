# Gather all de novo and related motifs from HOMER output and plot them (per WGCNA module)

library(homerkit)
library(reshape2)
library(tidyr)
library(stringr)

key_path = "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/ABC_downstream/Genename_Motif.txt"

HOMER_result_path = "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/HOMER/output"

output_path = "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/analysis/HOMER/output"

## function
read_denovo_html_results = function(path, homer_dir = TRUE) {
  if (homer_dir == TRUE) {
    path = paste0(path, "/homerResults")
  }
  if (!file.exists(path)) {
    warning(paste("No files found"))
    return(NULL)
  }
  
  # Match correct html files
  filenames = list.files(path, pattern = "*.info.html")
  split <- strsplit(filenames, "motif")  
  split <- as.numeric(sapply(split, function(x) x <- sub(".info.html", "", x[2])))
  filenames <- filenames[order(split)]
  
  df_list = list()
  
  for (f in filenames) {
    print(f)
    ## Read in  html file
    html = readLines(paste(path,f,sep = "/"))
    
    
    # get number of motif from file
    mypattern = "motif([^<]*).info.html"
    n = gsub(mypattern,
             '\\1',
             grep(mypattern, f, value = TRUE))
    
    # Create dataframe
    df = data.frame(matrix(ncol = 13, nrow = 1))
    colnames(df) =   c(
      'motif_name',
      'consensus',
      'p_value',
      'log_p_value',
      "info_content",
      'tgt_num',
      'tgt_pct',
      'bgd_num',
      'bgd_pct',
      'tgt_pos',
      'bgd_pos',
      'strand_bias',
      'multiplicity'
    )
    # Main header
    mypattern = '<H2>Information for ([^<]*)</H2>'
    df$motif_name = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = paste('.*-([^<]*) \\(Motif ',
                      n,
                      '\\)</H2>',
                      sep = "")
    df$consensus = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    # Other characteristics: p-value, log-pvalue etc as in main function
    mypattern = '<TR><TD>p-value:</TD><TD>([^<]*)</TD></TR>'
    df$p_value = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = '<TR><TD>log p-value:</TD><TD>([^<]*)</TD></TR>'
    df$log_p_value = gsub(mypattern, '\\1', grep(mypattern, html, value =
                                                   TRUE))
    
    mypattern = '<TR><TD>Information Content per bp:</TD><TD>([^<]*)</TD></TR>'
    df$info_content = gsub(mypattern, '\\1', grep(mypattern, html, value =
                                                    TRUE))
    
    mypattern = '<TR><TD>Number of Target Sequences with motif</TD><TD>([^<]*)</TD></TR>'
    df$tgt_num = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = '<TR><TD>Percentage of Target Sequences with motif</TD><TD>([^<]*)</TD></TR>'
    df$tgt_pct = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = '<TR><TD>Number of Background Sequences with motif</TD><TD>([^<]*)</TD></TR>'
    df$bgd_num = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = '<TR><TD>Percentage of Background Sequences with motif</TD><TD>([^<]*)</TD></TR>'
    df$bgd_pct = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = '<TR><TD>Average Position of motif in Targets</TD><TD>([^<]*)</TD></TR>'
    df$tgt_pos = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = '<TR><TD>Average Position of motif in Background</TD><TD>([^<]*)</TD></TR>'
    df$bgd_pos = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = '<TR><TD>Strand Bias \\(log2 ratio \\+ to \\- strand density\\)</TD><TD>([^<]*)</TD></TR>'
    df$strand_bias = gsub(mypattern, '\\1', grep(mypattern, html, value =
                                                   TRUE))
    
    mypattern = '<TR><TD>Multiplicity \\(# of sites on avg that occur together\\)</TD><TD>([^<]*)</TD></TR>'
    df$multiplicity = gsub(mypattern, '\\1', grep(mypattern, html, value =
                                                    TRUE))
    
    df_list[[n]][["Motif_information"]] = df
    
    ########### new information
    
    mypattern = '<H4>([^<]*)</H4>'
    length_df = length(gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE)))
    
    df = data.frame(matrix(ncol = 9, nrow = length_df))
    
    colnames(df) =   c(
      'motif_name',
      'ID',
      'database',
      "rank",
      'score',
      'offset',
      'orientation',
      "original_alignment",
      "matched_alignment"
    )
    
    # Known motif matches
    mypattern = '<H4>([^<]*)</H4>'
    df$motif_name = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    
    
    
    df <- df %>%
      tidyr::separate_(
        'motif_name',
        c('motif_name', 'ID', 'database'),
        '/',
        extra = "drop",
        fill = "right"
      )
    
    
    ## Detect if parentheses are present in motif_name
    ## to break apart into motif_name vs. motif_family
    cond <- stringr::str_detect(df$motif_name, '\\(') %>%
      sum(., na.rm = TRUE) > 0
    if (cond == TRUE) {
      df <- df %>%
        tidyr::separate_(
          'motif_name',
          c('motif_name', 'motif_family'),
          '\\(',
          extra = "drop",
          fill = "right"
        )
      df$motif_family <-
        stringr::str_replace(df$motif_family, '\\)', '')
    }
    
    # Ranks
    mypattern = '<TR><TD>Match Rank:</TD><TD>([^<]*)</TD></TR>'
    df$rank = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    # Scores
    mypattern = '<TR><TD>Score:</TD><TD>([^<]*)</TD</TR>'
    df$score = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    # Offset
    mypattern = '<TR><TD>Offset:</TD><TD>([^<]*)</TD</TR>'
    df$offset = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    # Orientation
    mypattern = '<TR><TD>Orientation:</TD><TD>([^<]*)</TD></TR>'
    df$orientation = gsub(mypattern, '\\1', grep(mypattern, html, value =
                                                   TRUE))
    
    # Alignment
    # original
    
    mypattern = paste('.*-([^<]*) \\(Motif ',
                      n,
                      '\\)</H2>',
                      sep = "")
    df$original_alignment = gsub(mypattern, '\\1', grep(mypattern, html, value =
                                                          TRUE))
    
    # matched
    mypattern = '.+>([^<]+)</FONT></TD></TR></TABLE>'
    df$matched_alignment = gsub(mypattern, '\\1', grep(mypattern, html, value =
                                                         TRUE))
    
    df_list[[n]][["Matches_to_known_motifs"]] = df
    
    
  }
  return(df_list)
}

# list all module directories
dirs = list.dirs(HOMER_result_path, full.names = F, recursive = F)


results = list()
novo_motif_matches = list()
novo_motif_table = list()
known_motif_table = list()
for (d in dirs) {
  results[[d]] <- read_homer_output(folder = paste(HOMER_result_path, d, sep = "/"))
  #novo_motif_matches[[d]] = results[[d]]$novel_motif_tfs
  novo_motif_matches[[d]] = read_denovo_html_results(paste(HOMER_result_path, d, sep = "/"))
  novo_motif_matches[[d]] <- lapply(1:length(novo_motif_matches[[d]]), function (i) {novo_motif_matches[[d]][[i]]$Matches_to_known_motifs}) 
  novo_motif_matches[[d]] = do.call("rbind",novo_motif_matches[[d]])
  novo_motif_matches[[d]]$motif = rep(paste("motif",
                                            1:(nrow(novo_motif_matches[[d]])/10),
                                            sep="_"),
                                      each=10)
  novo_motif_table[[d]] = results[[d]]$novel_motif_table
  novo_motif_table[[d]]$motif = paste("motif",
                                            1:(nrow(novo_motif_table[[d]])),
                                            sep="_")
  known_motif_table[[d]] = results[[d]]$known_motif_table

  # Resnaming columns of known motif table for rbind
  colnames(known_motif_table[[d]]) = c("motif_name","consensus", "p_value", "log_p_value",
                                       "q_value_benjamini", "number_of_target_sequences_with_motif",
                                       "percent_of_target_sequences_with_motif",
                                       "number_of_background_sequences_with_motif",
                                       "percent_of_background_sequences_with_motif")
}

novo_motif_matches = do.call("rbind",novo_motif_matches)
novo_motif_table = do.call("rbind",novo_motif_table)
known_motif_table = do.call("rbind",known_motif_table)

novo_motif_matches$stage = gsub("\\..*","",rownames(novo_motif_matches))
novo_motif_table$stage = gsub("\\..*","",rownames(novo_motif_table))
known_motif_table$stage = gsub("\\..*","",rownames(known_motif_table))


# Add TF info

key = read.table(key_path)
colnames(key) = c("TF","motif","inference")
key$TF = as.character(key$TF)
key$motif = as.character(key$motif)


novo_motif_matches$TF = key[match( novo_motif_matches$motif_name, key$motif),"TF"]
novo_motif_table$TF = key[match( novo_motif_table$best_match, key$motif),"TF"]
known_motif_table$TF = key[match( known_motif_table$motif_name, key$motif),"TF"]

novo_motif_matches$inference = key[match( novo_motif_matches$motif_name, key$motif),"inference"]
novo_motif_table$inference = key[match( novo_motif_table$best_match, key$motif),"inference"]
known_motif_table$inference = key[match( known_motif_table$motif_name, key$motif),"inference"]

# fix very low p-values: higher than  exp(-700.000) go to 0 and can't be plotted later
novo_motif_table[novo_motif_table$p_value == 0.00,"p_value"] = exp(-700.000)
#
write.table(novo_motif_matches, file = paste0(output_path,"/novo_motif_matches.txt"),
            quote = F, row.names = T, sep = "\t")
write.table(novo_motif_table, file = paste0(output_path,"/novo_motif_table.txt"),
            quote = F, row.names = T, sep = "\t")
write.table(known_motif_table, file = paste0(output_path,"/known_motif_table.txt"),
            quote = F, row.names = T, sep = "\t")
