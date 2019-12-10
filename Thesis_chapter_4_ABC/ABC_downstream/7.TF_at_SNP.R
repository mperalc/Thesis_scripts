## find TF binding motif at SNPs

library(atSNP)

.structure <- function(pval_mat) {
  if(is.matrix(pval_mat)) {
    if(nrow(pval_mat) > 1) {
      id <- apply(pval_mat[, c(2, 4)], 1, which.min)
      return(cbind(pval_mat[, c(1, 3)][cbind(seq_along(id), id)],
                   pval_mat[, c(2, 4)][cbind(seq_along(id), id)]))
    }
  }
  
  pval_mat <- as.vector(pval_mat)
  id <- which.min(pval_mat[c(2, 4)])
  return(cbind(pval_mat[c(1, 3)][id],
               pval_mat[c(2, 4)][id]))
}

CheckSameLength <- function(x) {
  if(length(x) == 1) {
    return(TRUE)
  }
  return(var(unlist(sapply(x, length))) == 0)
}

myStrSplit <- function(x, split) {
  ret <- list(seq_along(x))
  for(i in seq_along(x)) {
    ret[[i]] <- x[i]
    for(sp in split) {
      ret[[i]] <- unlist(strsplit(ret[[i]], sp))
      ret[[i]] <- ret[[i]][nchar(ret[[i]]) > 0]
      if(length(ret[[i]]) == 0)
        break
    }
  }
  return(ret)
}

#' @import doParallel
startParallel <- function(ncores) {
  if(.Platform$OS.type == "unix") {
    registerDoParallel(ncores)
  } else {
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    Return("cl")
  }
}

#' @import doParallel
endParallel <- function() {
  if(.Platform$OS.type != "unix") {
    stopCluster(cl)
  }
}


LoadSNPData_updated <- function(filename = NULL, genome.lib = "BSgenome.Hsapiens.UCSC.hg19",
                        snp.lib = "SNPlocs.Hsapiens.dbSNP.20120608",
                        snpids = NULL, half.window.size = 30, default.par = FALSE,
                        mutation = FALSE, ...) {
  message("Starting analysis...")
  useFile <- FALSE
  rsid.rm <- rsid.missing <- rsid.duplicate <- rsid.na <- NULL
  if(!is.null(filename)) {
    if(file.exists(filename)) {
      useFile <- TRUE
    }
  }
  if(useFile) {
    if(!is.null(snpids)) {
      message("Warning: load SNP information from 'filename' only. The argument 'snpids' is overridden.")
    }
    tbl <- read.table(filename, header = TRUE, stringsAsFactors = FALSE, ...)
    tbl<-tbl[c("snpid", "chr", "snp", "a1", "a2")]
    ## check if the input file has the required information
    if(sum(!c("snp", "chr", "a1", "a2", "snpid") %in% names(tbl)) > 0) {
      stop("Error: 'filename' must be a table containing 'snp' and 'chr' columns.")
    }
    snpid.index <- seq(nrow(tbl))
  } else {
    message("Checking SNPids..")
    if(is.null(snpids)) {
      stop("Error: either 'snpids' should be a vector, or 'filename' should be the file name that contains the SNP information.")
    }
    snpid.index <- seq(length(snpids))
    ## load the corresponding snp library
    library(package = snp.lib, character.only = TRUE)
    rsid.missing.all <- NULL
    while(TRUE) {
      snps <- get(snp.lib)
      snp.loc <- tryCatch({snpsById(snps, snpids)}, error = function(e) return(e$message))
      ## remove rsids not included in the database
      if(class(snp.loc) == "character") {
        rsid.missing <- myStrSplit(snp.loc, split = c(": ", "\n"))[[1]][-1]
        rsid.missing <- myStrSplit(rsid.missing, split = c(",", " "))[[1]]
        if(length(rsid.missing) > 1) {
          if(as.integer(rsid.missing[length(rsid.missing)]) <= as.integer(rsid.missing[length(rsid.missing) - 1])) {
            rsid.missing <- rsid.missing[-length(rsid.missing)]
          }
        }
        rsid.missing <- paste("rs", rsid.missing, sep = "")
        rsid.missing.all <- c(rsid.missing.all, rsid.missing)
        snpids <- snpids[!snpids %in% rsid.missing]
        snp.loc <- tryCatch({snpsById(snps, snpids)}, error = function(e) return(e$message))
      } else {
        break
      }
    }
    # if(!is.null(rsid.missing.all)) {
    #   message("Warning: the following rsids are not included in the database and discarded: ")
    #   message(paste(rsid.missing.all, collapse = ", "))
    #   rsid.missing <- rsid.missing.all
    # }
    
    snp.alleles <- snpsById(snps, snpids)
    snp.alleles <- IUPAC_CODE_MAP[snp.alleles@elementMetadata$alleles_as_ambig]
    gr = snpsById(snps, snpids)
    gr = as(gr, "GRanges")
    snp.strands <- as.character(GenomicRanges::as.data.frame(gr)$strand)
    message("Checking nuumber of alleles")
    if(sum(nchar(snp.alleles) > 2) > 0) {
      message("Warning: the following SNPs have more than 2 alleles. All pairs of nucleotides are considered as pairs of the SNP and the reference allele:")
      rsid.duplicate <- snpids[nchar(snp.alleles) > 2]
      message(paste(rsid.duplicate, collapse = ", "))
    }
    ## retain only SNPs with >= 2 alleles
    message("Retain only SNPs with >= 2 alleles...")
    tbl <- NULL
    for(nalleles in 2:4) {
      ids <- which(sapply(snp.alleles, nchar) == nalleles)
      if(length(ids) == 0) {
        next
      }
      snp.loc.n = snpsById(snps, snpids)
      snp.loc.n = as(snp.loc.n, "GRanges")
      chr = snp.loc.n@seqnames
      snp.loc.n = snp.loc.n@ranges@start
      snp.alleles.n <- snp.alleles[ids]
      snp.ids.n <- snpids[ids]
      snp.alleles.n <- strsplit(snp.alleles.n, "")
      snp.strands.n <- snp.strands[ids]
      ## get all pairs of alleles
      for(i_allele1 in seq(nalleles - 1)) {
        for(i_allele2 in (i_allele1 + 1):nalleles) {
          a1 <- sapply(snp.alleles.n, function(x) x[i_allele1])
          a2 <- sapply(snp.alleles.n, function(x) x[i_allele2])
          
          ## revert the alleles on the reverse strand
          id.rev <- which(snp.strands.n != "+")
          if(length(id.rev) > 0) {
            rev.codes <- c("A", "C", "G", "T")
            names(rev.codes) <- rev(rev.codes)
            a1[id.rev] <- rev.codes[a1[id.rev]]
            a2[id.rev] <- rev.codes[a2[id.rev]]
          }
          tbl <- rbind(tbl,
                       data.frame(
                         snp = snp.loc.n,
                         chr = paste0("chr", chr),
                         a1 = as.character(a1),
                         a2 = as.character(a2),
                         snpid = as.character(snp.ids.n),
                         index = snpid.index[ids])
          )
        }
      }
    }
    
    tbl <- tbl[order(tbl$index), ]
    if(!is.null(filename)) {
      write.table(tbl[, c('snp', 'chr', 'a1', 'a2', 'snpid')], file = filename, row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
    snpid.index <- tbl$index
  }
  
  ## load the corresponding genome version
  library(package = genome.lib, character.only = TRUE)
  species<-get(strsplit(genome.lib, "[.]")[[1]][2])
  seqvec <- BSgenome::getSeq(species,
                   as.character(tbl$chr),
                   start=tbl$snp - half.window.size,
                   end=tbl$snp + half.window.size,
                   as.character=TRUE)
  codes <- seq(4)
  names(codes) <- c("A", "C", "G", "T")
  sequences <- sapply(seqvec, function(x) codes[strsplit(x, "")[[1]]])
  sequences <- matrix(sequences, ncol = length(seqvec))
  snpid.output <- as.character(tbl$snpid)
  rownames(sequences) <- NULL
  a1 <- codes[as.character(tbl$a1)]
  a2 <- codes[as.character(tbl$a2)]
  names(a1) <- names(a2) <- NULL
  keep.id <- which(apply(sequences, 2, function(x) sum(is.na(x))) == 0)
  if(length(keep.id) < nrow(tbl)) {
    message("Warning: the following rows are discarded because the reference genome sequences contain non ACGT characters:")
    rsid.na <- tbl[-keep.id, ]$snpid
    print(tbl[-keep.id, ])
  }
  ## remove sequences containing non ACGT characters
  sequences <- sequences[, keep.id, drop = FALSE]
  a1 <- a1[keep.id]
  a2 <- a2[keep.id]
  snpid.output <- snpid.output[keep.id]
  snpid.index <- snpid.index[keep.id]
  ## whether use the default parameters
  if(!default.par) {
    transition <- .Call("transition_matrix", sequences, package = "atSNP")
    prior <- apply(transition, 1, sum)
    prior <- prior / sum(prior)
    transition <- transition / apply(transition, 1, sum)
    names(prior) <- colnames(transition) <- rownames(transition) <- c("A", "C", "G", "T")
  } else {
    data(default_par)
  }
  if(!mutation) {
    ## real SNP data
    a1.ref.base.id <- which(a1 == sequences[half.window.size + 1, ])
    a2.ref.base.id <- which(a2 == sequences[half.window.size + 1, ])
    ## store SNPs that have the same base in SNP and REF alleles only once
    a2.ref.base.id <- a2.ref.base.id[!a2.ref.base.id %in% a1.ref.base.id]
    discard.id <- setdiff(seq_along(a1), c(a1.ref.base.id, a2.ref.base.id))
    if(length(discard.id) > 0) {
      message("Warning: the following sequences are discarded because the reference nucleotide matches to neither a1 nor a2:")
      rsid.rm <- as.character(tbl[keep.id[discard.id], ]$snpid)
      message("snpid\tchr\tsnp\ta1\ta2")
      message(paste(apply(tbl[keep.id[discard.id], c("snpid", "chr", "snp", "a1", "a2")], 1, function(x) paste(x, collapse = "\t")), collapse = "\n"))
    }
  } else {
    ## single nucleotide mutation data
    a1.ref.base.id <- seq_along(a1)
    a2.ref.base.id <- numeric(0)
  }
  sequences <- sequences[, c(a1.ref.base.id, a2.ref.base.id), drop = FALSE]
  snpid.output <- snpid.output[c(a1.ref.base.id, a2.ref.base.id)]
  ref.base <- c(a1[a1.ref.base.id], a2[a2.ref.base.id])
  snp.base <- c(a2[a1.ref.base.id], a1[a2.ref.base.id])
  snpid.index <- snpid.index[c(a1.ref.base.id, a2.ref.base.id)]
  ## Keep the order of SNPs as in the input file
  if(useFile) {
    output.index = seq(ncol(sequences))
  } else {
    output.index = order(snpid.index)
  }
  return(list(
    sequence_matrix= matrix(sequences[, output.index], nrow=2*half.window.size+1),
    ref_base = ref.base[output.index],
    snp_base = snp.base[output.index],
    snpids = snpid.output[output.index],
    transition = transition,
    prior = prior,
    rsid.na = rsid.na,
    rsid.rm = rsid.rm,
    rsid.duplicate = rsid.duplicate,
    rsid.missing = rsid.missing
  ))
}

# load motif library
pwms = LoadMotifLibrary(filename = "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/PWMs/PWMs_filtered.txt",
                        tag = ">",skipcols = 0, skiprows = 1, transpose = F, field = 1,
                        sep = c(">"," ","\t"), pseudocount = 0)

# load SNP data
SNP = read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Integration/ABC/Predictions/credset_SNPs_in_ABC_enhancers_rsID.txt")
SNP_info = LoadSNPData_updated(snpids = as.character(SNP$V1[1:50]),
                       genome.lib = "BSgenome.Hsapiens.UCSC.hg19",
                       snp.lib = "SNPlocs.Hsapiens.dbSNP144.GRCh37",
                       half.window.size = 30)

atsnp.scores = ComputeMotifScore(pwms,SNP_info,ncores = 2)
