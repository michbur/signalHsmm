#' Protein test
#'
#' Checks if an object is a protein.
#'
#' @param object \code{character} vector where each elemenents represent one amino acid.
#' @return \code{TRUE} or \code{FALSE}.
#' @export

is_protein <- function(object) {
  all(toupper(object) %in% a()[-1])
}