#' Degenerate protein sequence
#'
#' 'Degenerates' protein sequence by aggregating aminoacids to bigger groups.
#' 
#' @param seq \code{character} vector representing single aminoacid sequence.
#' @param aa_group list of aminoacid groups to which sequence should be aggregated.
#' @keywords manip
#' @return a \code{character} vector.
#' @export
#' @keywords manip
#' @examples
#' sample_seq <- sample(seqinr:::a()[-1], 30, replace = TRUE)
#' table(sample_seq)
#' 
#' #compared with aggregated sequence
#' deg_seq <- degenerate(sample_seq, aaaggregation)
#' table(deg_seq)



degenerate <- function(seq, aa_group) {
  for (i in 1L:length(aa_group)) {
    seq[seq %in% aa_group[[i]]] <- i
  }
  seq
}


