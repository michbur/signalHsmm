#' Train signal.hsmm_model object
#' 
#' @param train_data training data.
#' @param aa_group method of aggregating amino acids.
#' @param max_length maximum length of signal peptide.
#' @export
#' @return object of class \code{signal.hsmm_model}

train_hsmm <- function(train_data, aa_group, max.length = 32) {
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
  
  list(aa_group = aa_group, pipar = pipar, tpmpar = tpmpar, od = od, 
       overall.probs.log = overall.probs.log, params = params)
  #change output to normal decision
  cbind(prob.sig = exp(decisions[,1] - decisions[,2]), 
        sig.end = decisions[, 3])
  #debug version of output
  #   list(prob.sig = exp(decisions[,1] - decisions[,2]), 
  #        sig.end = decisions[, 3], pipar = pipar, tpmpar = tpmpar, od = od, 
  #        overall.probs.log = overall.probs.log, params = params, lengths = lengths)
}