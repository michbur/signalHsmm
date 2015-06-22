#' Read data from UniProt database
#'
#' Read data saved in UniProt original flat text format.
#'
#' @param connection a \code{\link{connection}} to UniProt data in text format.
#' @param what \code{NULL} or a single \code{character} determining which information 
#' should be extracted. Currently can have following values: \code{"signal"}, 
#' \code{"transit"} or \code{NULL}.
#' @param euk logical value if data has an eukaryotic origin.
#' @keywords manip
#' @return a list of sequences. Each element has a class \code{\link[seqinr]{SeqFastaAA}}.
#' Attribute \code{sig} contains the range (start and end) of signal peptide. Attributes 
#' \code{OS} and \code{OC} represents respectively OS and OC fields in the protein
#' description. Sequence with more than one cleavage site or atypical aminoacids 
#' are removed without any notice.
#' @export
#' @keywords manip

read_uniprot <- function(connection, what = "signal", euk) {
  
  all_lines <- readLines(connection)
  
  all_seqs <- preliminary_seqs(all_lines, what = what) 
  
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
    
    if(!is.null(what)) {
      sig <- suppressWarnings(as.numeric(strsplit(strsplit(all_lines[sure_seqs[4,i]], 
                                                           paste0(toupper(what), "       "))[[1]][2], " ")[[1]]))
      sig <- as.numeric(na.omit(sig))
      attr(ith_seq, "sig") <- sig
    }
    
    os <- all_lines[sure_seqs["os_start", i]:sure_seqs["os_end", i]]
    attr(ith_seq, "OS") <- paste0(vapply(strsplit(os, "OS   "), function(single_os) single_os[2], "a"), collapse = " ")
    
    oc <- all_lines[sure_seqs["oc_start", i]:sure_seqs["oc_end", i]]
    attr(ith_seq, "OC") <- paste0(vapply(strsplit(oc, "OC   "), function(single_oc) single_oc[2], "a"), collapse = " ")
    
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
preliminary_seqs <- function(all_lines, what) {
  prot_ids <- grep("\\<ID   ", all_lines)
  seqs_start <- grep("\\<SQ   ", all_lines) + 1
  seqs_end <-  grep("^//", all_lines) - 1
  #test if protein has information regarding signal
  
  prot_sig <- rep(NA, length(prot_ids))
  
  signals <- grep(paste0("FT   ", toupper(what)), all_lines)
  
  all_ids <- sort(c(prot_ids, signals), method = "quick")
  prot_sig[which(all_ids %in% signals) - 1L:length(signals)] <- signals
  
  
  os_list <- get_block(grep("\\<OS   ", all_lines))
  os_start <- vapply(os_list, function(i) i[1], 0)
  os_end <- vapply(os_list, function(i) i[length(i)], 0)
  

  oc_list <- get_block(grep("\\<OC   ", all_lines))
  oc_start <- vapply(oc_list, function(i) i[1], 0)
  oc_end <- vapply(oc_list, function(i) i[length(i)], 0)
  
  res <- rbind(prot_ids, 
               seqs_start, seqs_end, 
               prot_sig,
               os_start, os_end,
               oc_start, oc_end)
  
  rownames(res) <- c("prot_ids", "seqs_start", "seqs_end", "prot_sig", 
                     "os_start", "os_end",
                     "oc_start", "oc_end")
  res
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

#gets list of starts and ends of specific blocks in UniProt data
get_block <- function(block_indices) {
  block_list <- list()
  tmp_block <- c()
  
  for(i in block_indices) {
    if(is.null(tmp_block[length(tmp_block)])) {
      tmp_block <- i
    } else {
      if(tmp_block[length(tmp_block)] + 1 == i) {
        tmp_block <- c(tmp_block, i)
      } else {
        block_list <- c(block_list, list(tmp_block))
        tmp_block <- i
      } 
    }
  }
  
  c(block_list, list(tmp_block))
}
