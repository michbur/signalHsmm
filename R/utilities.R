#' signal.hsmm - prediction of signal peptides
#'
#' signal.hsmm predicts presence of signal peptides using the hidden
#' semi-Markov models.
#' 
#' ...
#' 
#' @importFrom seqinr read.fasta
#' @docType package
#' @name signal.hsmm



#READ UNIPROT DATA -----------------------------
#helper function to get seqs from  .txt files
preliminary_seqs <- function(all_lines, signal) {
  prot_ids <- grep("\\<ID   ", all_lines)
  seqs_start <- grep("\\<SQ   ", all_lines) + 1
  seqs_end <-  grep("^//", all_lines) - 1
  #test if protein has information regarding signal
  
  prot_sig <- rep(NA, length(prot_ids))
  
  if (signal) {
    signals <- grep("FT   SIGNAL", all_lines)
    all_ids <- sort(c(prot_ids, signals), method = "quick")
    prot_sig[which(all_ids %in% signals) - 1L:length(signals)] <- signals
  }
  
  rbind(prot_ids, seqs_start, seqs_end, prot_sig)
}

#removes proteins with probable or potential signal peptides
#removes proteins without cleavage site for signalase
#return indices of proteins which surely have signal peptide
remove_unsure <- function(all_lines, all_seqs) {
  signals <- all_seqs[4, ]
  #remove with unknown signal peptide end
  (1L:ncol(all_seqs))[-unique(unlist(lapply(c(">", "<1", "?", "Or "), function(pattern) 
    grep(pattern, all_lines[signals], fixed = TRUE))))]
}

#proteins not encoded in nucleus
remove_nonnuclear <- function(all_lines, all_seqs) {
  nucl <- grep("OG   ", all_lines)
  all_ids <- sort(c(all_seqs[1, ], nucl), method = "quick")
  setdiff(1L:ncol(all_seqs), which(all_ids %in% nucl) - 1L:length(nucl))
}

#removes noncleavable seqs
remove_notcleaved <- function(all_lines, all_seqs) {
  not_cleaved <- grep("Not cleaved", all_lines)
  all_ids <- sort(c(all_seqs[1, ], not_cleaved), method = "quick")
  setdiff(1L:ncol(all_seqs), which(all_ids %in% not_cleaved) - 1L:length(not_cleaved))
}


#remove seqs with atypical and/or not identified aas
get_atyp <- function(list_prots) {
  which(colSums(sapply(list_prots, function(prot) 
    sapply(c("B", "U", "X", "Z"), function(atyp_aa)
      atyp_aa %in% prot))) > 0)
}

#' Read data from UniProt database
#'
#' Read data saved in UniProt original flat text format.
#'
#' @param \code{\link{connection}} to UniProt data.
#' @param euk logical value if data has eukaryotic origin.
#' @keywords manip
#' @return a list of sequences. Each element has class \code{\link[seqinr]{SeqFastaAA}}.
#' Slot \code{sig} contains the range of signal peptide.
#' @export


read_uniprot <- function(connection, euk) {
  
  all_lines <- readLines(connection)
  
  all_seqs <- preliminary_seqs(all_lines, signal = TRUE) 
  
  #remove unsure signal peptides
  only_sure <- remove_unsure(all_lines, all_seqs)
  #remove notcleaved signal peptides
  cleaved <- remove_notcleaved(all_lines, all_seqs)
  sure_cleaved <- intersect(only_sure, cleaved) 
  
  if (euk) {
    #remove proteins directed to nucleus
    only_nonnuclear <- remove_nonnuclear(all_lines, all_seqs)
    sure_cleaved <- intersect(only_nonnuclear, sure_cleaved)
  } 
  
  sure_seqs <- all_seqs[, sure_cleaved]
  
  
  list_prots <- lapply(1L:ncol(sure_seqs), function(i) {
    start_seq <- sure_seqs[2,i]
    end_seq <- sure_seqs[3,i]
    
    ith_seq <- strsplit(gsub(" ", "", paste0(all_lines[start_seq:end_seq], collapse = "")), "")[[1]]
    
    class(ith_seq) <- "SeqFastaAA"
    aa_name <- strsplit(all_lines[sure_seqs[1,i]], "   ")[[1]][2]
    attr(ith_seq, "name") <- aa_name
    attr(ith_seq, "Annot") <- paste0(">", aa_name)
    attr(ith_seq, "class") <- "SeqFastaAA"
    
    #to do - think about something smarter than suppressWarnings
    sig <- suppressWarnings(as.numeric(strsplit(strsplit(all_lines[sure_seqs[4,i]], "SIGNAL       ")[[1]][2], " ")[[1]]))
    sig <- as.numeric(na.omit(sig))
    attr(ith_seq, "sig") <- sig
    #line is preserved just to have an additional source of information
    #attr(ith_seq, "line") <- all_lines[sure_seqs[4,i]]
    
    ith_seq
  })
  
  names(list_prots) <- sapply(list_prots, function(i) attr(i, "name"))
  atypical <- get_atyp(list_prots)
  if (length(atypical) > 0){
    list_prots[-atypical]
  } else {
    list_prots
  }
}


# HEURISTIC ALGORITHM

#' Localize n-, h- and c-region in signal peptide
#'
#' A heuristic algorithm able to find borders between distinct regions
#' constituting signal peptides 
#'
#' @param protein a vector of amino acids or object of class 
#' \code{\link[seqinr]{SeqFastaAA}}
#' @param signal range of signal peptide. If \code{NULL}, the attribute \code{sig}
#' of \code{protein} will be used.
#' @return a vector of length 4 containing positions of:
#' \enumerate{
#'   \item start of n-region,
#'   \item start of h-region,
#'   \item start of c-region,
#'   \item cleavage site.
#' }
#' @references Henrik Nielsen, Anders Krogh (1998). Prediction of signal peptides
#' and signal anchors by a hidden Markov model. Proc. Sixth Int. Conf. on 
#' Intelligent Systems for Molecular Biology.
#' @export


#function to find n, h and c regions in signal peptide
find_nhc <- function(protein, signal = NULL) {
  if (is.null(signal)) 
    signal <- attr(protein, "sig")
  
  sig <- protein[signal[1]:signal[2]]
  start_c <- length(sig) - 2
  
  #noh number of hydrophobic residues
  noh <- 0
  while(noh < 2 && start_c > 1) {
    start_c <- start_c - 1
    noh <- noh + ifelse(sig[start_c] %in% c("A", "I", "L", "M", "F", "W", "V"), 1, -1)
    noh <- ifelse(noh < 0, 0, noh)
  }
  start_c <- start_c + 2
  
  start_h <- start_c - 6
  #if statement to prevent negative values
  if (start_h > 1) {
    #nonh number of nonhydrophobic residues
    nonh <- 0
    #noc number of charged
    noc <- 0
    while(nonh < 3 && noc == 0 && start_h > 1) {
      start_h <- start_h - 1
      nonh <- nonh + ifelse(sig[start_h] %in% c("A", "I", "L", "M", "F", "W", "V"), -1, 1)
      nonh <- ifelse(nonh < 0, 0, nonh)
      noc <- ifelse(sig[start_h] %in% c("R", "H", "K", "D", "E"), 1, 0)
    }
  } else {
    start_h <- 1
  }
  
  #prestart_c = start_c - 1 - prevents start_h > start_c
  prestart_c <- start_c - 1
  noh <- 0
  while(noh == 0 && start_h < prestart_c) {
    start_h <- start_h + 1
    noh <- noh + ifelse(sig[start_h] %in% c("A", "I", "L", "M", "F", "W", "V"), 1, 0)
  }
  #c(start_n = signal[1], start_h = start_h, start_c = start_c, cs = signal[2])
  c(start_n = 1, start_h = start_h, start_c = start_c, cs = signal[2])
}


# SIGNAL-HSMM ------------------------------------
degenerate <- function(seq, aa_group) {
  for (i in 1L:length(aa_group)) {
    seq[seq %in% aa_group[[i]]] <- i
  }
  seq
}



calc_t <- function(list_prots, aa_list) {
  nhc <- t(vapply(list_prots, find_nhc, rep(0, 4)))
  
  n_region <- NULL
  h_region <- NULL
  c_region <- NULL
  rest <- NULL
  
  for(i in 1L:length(list_prots)){
    region_starts <- nhc[i, ]
    n_region <- c(n_region, list_prots[[i]][1:(region_starts[2] - 1)])
    h_region <- c(h_region, list_prots[[i]][region_starts[2]:(region_starts[3] - 1)])
    c_region <- c(c_region, list_prots[[i]][region_starts[3]:(region_starts[4] - 1)])
    rest <- c(rest, list_prots[[i]][region_starts[4]:length(list_prots[[i]])])
  }
  
  t1 <- rep(0, length(aa_list))
  temp <- table(degenerate(n_region, aa_list))
  t1[as.numeric(names(temp))] <- temp
  names(t1) <- 1:length(aa_list)
  
  t2 <- rep(0, length(aa_list))
  temp <- table(degenerate(h_region, aa_list))
  t2[as.numeric(names(temp))] <- temp
  names(t2) <- 1:length(aa_list)
  
  t3 <- rep(0, length(aa_list))
  temp <- table(degenerate(c_region, aa_list))
  t3[as.numeric(names(temp))] <- temp
  names(t3) <- 1:length(aa_list)
  
  t4 <- rep(0, length(aa_list))
  temp <- table(degenerate(rest, aa_list))
  t4[as.numeric(names(temp))] <- temp
  names(t4) <- 1:length(aa_list)
  
  len_c <- nhc[, "cs"] - nhc[, "start_c"]
  len_h <- nhc[, "start_c"] - nhc[, "start_h"]
  len_n <- nhc[, "start_h"] - nhc[, "start_n"]
  lengths <- matrix(c(len_n, len_h, len_c), ncol = 3)
  colnames(lengths) <- c("n", "h", "c")
  
  list(mean_cs = mean(nhc[, 4]), sd_cs = sd(nhc[, 4]), t1 = t1, t2 = t2, t3 = t3, t4 = t4, 
       lengths = lengths)
}

measure_region <- function(region, max.length = 32) {
  lengths <- table(region)
  res <- rep(0, max.length)
  lengths <- lengths[as.numeric(names(lengths))>0] #removing lengths smaller than 1
  
  start_l <- min(as.numeric(names(lengths)))
  end_l <- max(as.numeric(names(lengths)))
  if(prod(start_l:end_l %in% as.numeric(names(lengths)))){
    max_length <- length(lengths) #if all lengths are present in training set
  } else{
    max_length <- 1
    sl <- sum(lengths)
    while(sum(lengths[1:max_length])/sl <= 0.51) {
      max_length <- which.min(start_l:end_l %in% as.numeric(names(lengths))) - 1
      start_l <- start_l + 1  #to assure that max_length wouldn't be too small
      max_length <- ifelse(max_length == 0, length(lengths), max_length)
    }
  }
  max_length <- min(max_length, max.length)
  
  prop_lengths <- lengths[1:max_length]/sum(lengths[1:max_length])
  res[as.numeric(names(prop_lengths))] <- prop_lengths
  res
}


signal_hsmm_train <- function(train_data, test_data, aa_group, max.length = 32) {
  train_data <- lapply(train_data, toupper)
  ts <- calc_t(train_data, aa_group)
  
  t1 <- ts[["t1"]]
  t2 <- ts[["t2"]]
  t3 <- ts[["t3"]]
  t4 <- ts[["t4"]]
  
  overall <- t4 #table(degenerate(unlist(analized_sequences), aa5))
  overall.probs <- overall/sum(overall)          
  overall.probs.log = log(overall.probs) #for viterbi
  
  lengths <- ts[["lengths"]]
  params <- apply(lengths, 2, measure_region, max.length = max.length)
  params <- cbind(params, rep(1/max.length, max.length))
  
  #setting params for hmm -------
  additional_margin = 10
  pipar <- c(1,0,0,0)
  tpmpar <- matrix(c(0, 1, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1,
                     0, 0, 0, 0), 4, byrow = TRUE)
  od <- matrix(c((t1/sum(t1))[1:4],
                 (t2/sum(t2))[1:4],
                 (t3/sum(t3))[1:4],
                 (t4/sum(t4))[1:4]), 4, byrow = TRUE)
  
  decisions <- signal_hsmm(test_data, aa_group, pipar = pipar, tpmpar = tpmpar, od = od, 
                           overall.probs.log = overall.probs.log, params = params)
  prob.sig <- exp(decisions[,1] - decisions[,2])
  
  #change output to normal decision
    data.frame(signal.peptide = round(1-1/(1+prob.sig), 0) == 1, 
          sig.start = 1,
          sig.end = decisions[, 3])
#   data.frame(signal.peptide = prob.sig, 
#              sig.start = 1,
#              sig.end = decisions[, 3])
}

signal_hsmm <- function(list_prot, aa_group, pipar, tpmpar, 
                        od, overall.probs.log, params, max.length = 50) { 
  t(vapply(list_prot, function(prot) {
    probka <- as.numeric(degenerate(toupper(prot), aa_group)[1:max.length])
    probka <- na.omit(probka)
    viterbi.res <- duration.viterbi(probka, pipar, tpmpar, od, params)
    viterbi_path <- viterbi.res[["path"]]
    c.site <- ifelse(any(viterbi_path==3), max(which(viterbi_path==3)), length(probka))
    prob.signal <- viterbi.res[["viterbi"]][c.site, viterbi_path[c.site]]
    prob.non <- Reduce(function(x, y) x + overall.probs.log[y], probka[1:c.site], 0)
    c(prob.signal = prob.signal, prob.non = prob.non, c.site = c.site)   
  }, c(0, 0, 0)))
}
