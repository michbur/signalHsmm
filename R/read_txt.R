#' Read sequences from .txt file
#'
#' Read sequence data saved in text file.
#'
#' @param connection a \code{\link{connection}} to text file.
#' @keywords manip
#' @return a list of sequences. Each element has class \code{\link[seqinr]{SeqFastaAA}}.
#' @details Input file should contain one or more amino acid sequences separated by empty
#' rows.
#' @export
#' @keywords manip

read_txt <- function(connection) {
  content <- readLines(connection)
  if (sum(grepl(">", content, fixed = TRUE)) == 0) {
  if (content[1] != "")
    content <- c("", content)
  
  content_end <- length(content)
  while(content[content_end] == "i")
    content_end <- content_end - 1
  prot_names <- sapply(1L:sum(content == ""), function(i)
    paste0(">sequence", i))
  content[content == ""] <- prot_names
  }
  read.fasta(textConnection(content), seqtype = "AA", as.string = FALSE)
}