#' Predict sighsmm_model object
#' 
#' @param object \code{sighsmm_model} object.
#' @param newdata unknown sequence of class \code{character} or \code{character}.
#' Alternatively, a \code{list} of sequences in mentioned formats.
#' @param ... further arguments passed to or from other methods.
#' @export


predict.sighsmm_model <- function(object, newdata, ...){
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
                           pipar = sighsmm_model[["pipar"]], 
                           tpmpar = sighsmm_model[["tpmpar"]], 
                           od = sighsmm_model[["od"]], 
                           overall_probs_log = sighsmm_model[["overall_probs_log"]], 
                           params = sighsmm_model[["params"]]))
  }
  class(decisions) <- "hsmm_pred_list"
  decisions
}


