#' signal.hsmm - prediction of signal peptides
#'
#' signal.hsmm predicts presence of signal peptides using the hidden
#' semi-Markov models.
#' 
#' @description Implementing hidden semi-Markov model and novel approach to sequence
#' analysis, signal.hsmm is new, highly accurate signal peptide predictor in eukaryotic
#' proteins. 
#' @importFrom seqinr read.fasta
#' @docType package
#' @name signal.hsmm
#' @examples
#' few_predictions <- run_signal.hsmm(benchmark_dat[1:3])
#' #summary all predictions
#' summary(few_predictions)
#' #summary one prediction
#' summary(few_predictions[[1]])
#' #plot one prediction
#' plot(few_predictions[[1]])
NULL



degenerate <- function(seq, aa_group) {
  for (i in 1L:length(aa_group)) {
    seq[seq %in% aa_group[[i]]] <- i
  }
  seq
}


