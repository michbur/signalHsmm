#' Read data from UniProt database
#'
#' Read data saved in UniProt original flat text format.
#'
#' @param connection a \code{\link{connection}} to UniProt data in text format.
#' @param euk logical value if data has eukaryotic origin.
#' @keywords manip
#' @return a list of sequences. Each element has class \code{\link[seqinr]{SeqFastaAA}}.
#' Slot \code{sig} contains the range of signal peptide.
#' @export
#' @keywords manip

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

