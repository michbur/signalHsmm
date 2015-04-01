#put aaaggregation here just to just R CHECK notes
aaaggregation <- list(c("k", "r", "h"), 
                      c("v", "i", "l", "m", "f", "w", "c", "u"),
                      c("s", "t", "n", "q"), 
                      c("d", "e", "a", "p", "y", "g"))
names(aaaggregation) = c("1", "2", "3", "4")

#' @name benchmark_dat
#' @title Benchmark data set for signal.hsmm
#' @description Lists eukaryotic proteins added to UniProt database release 
#' 2014_07 between 2011 and 2014 (140 proteins with signal peptide and 280 randomly 
#' sampled proteins without signal peptide). 
#' All proteins were used in the benchmark test to compare the performance
#' of signal.hsmm and other signal peptide predictors.
#' @docType data
#' @usage benchmark_dat
#' @format a list of \code{\link[seqinr]{SeqFastaAA}} objects. 
#' Slot \code{sig} contains the range of signal peptide (if any).
#' @source \href{http://www.uniprot.org/}{UniProt}
#' @examples summary(benchmark_dat)
#' @keywords datasets
NULL

#' @name aaaggregation
#' @title Scheme for amino acid aggregation
#' @description Amino acids are grouped together in larger sets based on their 
#' physicochemical properties important in  the recognition of signal peptide.
#' @docType data
#' @usage aaaggregation
#' @format a list of length four. Each element contains a \code{character} vector 
#' of amino acid names (one-letter abbreviations).
#' @keywords datasets
NULL