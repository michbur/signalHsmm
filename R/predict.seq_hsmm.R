#' Predict signal.hsmm_model object
#' 
#' @param object \code{signal.hsmm_model} object.
#' @param newdata unknown sequence of class \code{character} or \code{character}.
#' Alternatively, a \code{list} of sequences in mentioned formats.
#' @param ... further arguments passed to or from other methods.
#' @export


predict.seq_hsmm <- function(object, newdata, ...){
  if (class(newdata) == "numeric" || class(newdata) == "factor" || 
        class(newdata) == "data.frame" || class(newdata) == "matrix")
    stop("Input data must have class 'SeqFastaAA', 'character' or 'list'.")
  
  if(class(newdata) == "SeqFastaAA" || 
       class(newdata) == "character") {
    #single input
    decisions <- signal.hsmm_decision(newdata, aa_group = object[["aa_group"]], 
                                      pipar = object[["pipar"]], 
                                      tpmpar = object[["tpmpar"]], 
                                      od = object[["od"]], 
                                      overall_probs_log = object[["overall_probs_log"]], 
                                      params = object[["params"]])
    decisions <- list(decisions)
    names(decisions) <- attr(newdata, "name")
  } else {
    #list input
    decisions <- lapply(newdata, function(prot)
      signal.hsmm_decision(prot, aa_group = aaaggregation, 
                           pipar = signal.hsmm_model[["pipar"]], 
                           tpmpar = signal.hsmm_model[["tpmpar"]], 
                           od = signal.hsmm_model[["od"]], 
                           overall_probs_log = signal.hsmm_model[["overall_probs_log"]], 
                           params = signal.hsmm_model[["params"]]))
  }
  class(decisions) <- "hsmm_pred_list"
  decisions
}