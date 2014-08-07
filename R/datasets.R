#' @name benchmark_dat
#' @title Benchmark data set for signal.hsmm
#' @description This data set lists eukaryotic proteins added to UniProt database release 
#' 2014_07 between 2011 and 2014 (140 proteins with signal peptide and 280 randomly 
#' sampled proteins without signal peptide). 
#' All proteins were used in benchmark test comparing the performance
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
#' @format a list of length four containing one-letter name of amino acid
#' @keywords datasets
NULL