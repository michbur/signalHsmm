# SIGNAL-HSMM ------------------------------------

#' Predict presence of signal peptide in protein
#'
#' Using the hidden semi-Markov model predict presence of signal peptide in 
#' eukaryotic proteins.
#'
#' @param test_data single protein sequence (\code{character} vector) or list of 
#' sequences. It may be an object of class \code{\link[seqinr]{SeqFastaAA}}.
#' @return An object of class \code{hsmm_pred_list}.
#' @details Function \code{signal.hsmm} returns respectively probability of presence of 
#' signal peptide, start of signal peptide and the probable cleavage site localization.
#' If input consists from more than one sequence, result is a data.frame where each column
#' contains above values for different proteins.
#' @note Currently start of signal peptide is naively set as 1 amino acid.
#' @export
#' @seealso \code{\link{hsmm_pred_list}} \code{\link{hsmm_pred}} 
#' @keywords classif
#' @examples
#' #run signal.hsmm on one sequence
#' x1 <- run_signal.hsmm(benchmark_dat[[1]])
#' 
#' #run signal.hsmm on one sequence, but input is a character vector
#' x2 <- run_signal.hsmm(c("m", "a", "g", "k", "e", "v", "i", "f", 
#' "i", "m", "a", "l", "f", "i", "a", "v", "e", "s", "s", "p", "i", 
#' "f", "s", "f", "d", "d", "l", "v", "c", "p", "s", "v", "t", "s", 
#' "l", "r", "v", "n", "v", "e", "k", "n", "e", "c", "s", "t", "k", 
#' "k", "d", "c", "g", "r", "n", "l", "c", "c", "e", "n", "q", "n", 
#' "k", "i", "n", "v", "c", "v", "g", "g", "i", "m", "p", "l", "p", 
#' "k", "p", "n", "l", "d", "v", "n", "n", "i", "g", "g", "a", "v", 
#' "s", "e", "s", "v", "k", "q", "k", "r", "e", "t", "a", "e", "s", 
#' "l"))
#' 
#' #run signal.hsmm on list of sequences
#' x3 <- run_signal.hsmm(benchmark_dat[1:3])
#' #see summary of results
#' summary(x3)
#' #print results as data frame
#' pred2df(x3)
#' #summary one result
#' summary(x3[[1]])
#' plot(x3[[1]])

run_signal.hsmm <- function(test_data) {
  signal.hsmm_model <- structure(list(pipar = c(1, 0, 0, 0), 
                                      tpmpar = structure(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0), 
                                                         .Dim = c(4L, 4L)), 
                                      od = structure(c(0.203181196494605, 0.00235316265060241, 
                                                       0.069742553064986, 0.122654208889416, 0.360670685375614, 0.72484469126506, 
                                                       0.284975325524704, 0.288121595016022, 0.190922445671445, 0.0822430346385542, 
                                                       0.299423271300315, 0.230261863691714, 0.245225672458335, 0.190559111445783, 
                                                       0.345858850109995, 0.358962332402848), 
                                                     .Dim = c(4L, 4L)), 
                                      overall_probs_log = structure(c(-2.09838619260539, -1.24437268303286, -1.4685380799115, 
                                                                      -1.02453781966768), 
                                                                    .Names = c("1", "2", "3", "4")), 
                                      params = structure(c(0.102263430597218, 0.180256340332697, 
                                                           0.166621216253068, 0.0910826288519226, 0.103081538041996, 0.0820834469593673, 
                                                           0.0526315789473684, 0.0597218434687756, 0.0349059176438506, 0.024815925824925, 
                                                           0.0231797109353695, 0.0179983637851104, 0.0106353967821107, 0.00899918189255522, 
                                                           0.0092718843741478, 0.00654485955822198, 0.00545404963185165, 
                                                           0.00299972729751841, 0.00299972729751841, 0.00190891737114808, 
                                                           0.000818107444777748, 0.00136351240796291, 0.00136351240796291, 
                                                           0.00245432233433324, 0.000818107444777748, 0.00190891737114808, 
                                                           0.0016362148895555, 0.000818107444777748, 0.000545404963185165, 
                                                           0.000272702481592582, 0.000545404963185165, 0, 0.000544069640914037, 
                                                           0.000544069640914037, 0.00190424374319913, 0.00108813928182807, 
                                                           0.0043525571273123, 0.0160500544069641, 0.0356365614798694, 0.0710010881392818, 
                                                           0.0900435255712731, 0.126768226332971, 0.173830250272035, 0.142546245919478, 
                                                           0.12132752992383, 0.0807943416757345, 0.0489662676822633, 0.0375408052230686, 
                                                           0.0176822633297062, 0.0119695321001088, 0.00843307943416757, 
                                                           0.00326441784548422, 0.00217627856365615, 0.00108813928182807, 
                                                           0.00163220892274211, 0.000272034820457019, 0.000544069640914037, 
                                                           0, 0, 0, 0, 0, 0, 0, 0, 0.139858233369684, 0.136314067611778, 
                                                           0.320338058887677, 0.112050163576881, 0.155125408942203, 0.0602508178844057, 
                                                           0.0463467829880044, 0.014721919302072, 0.00872410032715376, 0.00327153762268266, 
                                                           0.002453653217012, 0.00054525627044711, 0, 0, 0, 0, 0, 0, 0, 
                                                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.03125, 0.03125, 0.03125, 
                                                           0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 
                                                           0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 
                                                           0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 
                                                           0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 
                                                           0.03125), 
                                                         .Dim = c(32L, 4L), 
                                                         .Dimnames = list(NULL, c("n", "h", "c", "")))), 
                                 .Names = c("pipar", "tpmpar", "od", "overall_probs_log", 
                                            "params"))
  if (class(test_data) == "numeric" || class(test_data) == "factor" || 
        class(test_data) == "data.frame" || class(test_data) == "matrix")
    stop("Input data must have class 'SeqFastaAA', 'character' or 'list'.")
  
  
  if(class(test_data) == "SeqFastaAA" || 
       class(test_data) == "character") {
    #single input
    decisions <- signal.hsmm_decision(test_data, aa_group = aaaggregation, 
                                      pipar = signal.hsmm_model[["pipar"]], 
                                      tpmpar = signal.hsmm_model[["tpmpar"]], 
                                      od = signal.hsmm_model[["od"]], 
                                      overall_probs_log = signal.hsmm_model[["overall_probs_log"]], 
                                      params = signal.hsmm_model[["params"]])
    decisions <- list(decisions)
    names(decisions) <- attr(test_data, "name")
  } else {
    #list input
    decisions <- lapply(test_data, function(prot)
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

signal.hsmm_decision <- function(prot, aa_group, pipar, tpmpar, 
                                 od, overall_probs_log, params) {
  if (length(prot) == 1) {
    prot <- strsplit(prot, "")[[1]]
    if ("name" %in% names(attributes(prot)))
      attr(prot, "name") <- "undefined_name"
    if (length(prot) == 1)
      stop("Input sequence is too short.")
  }
  deg_sample <- as.numeric(degenerate(toupper(prot)[1L:50], aa_group))
  #remove atypical amino acids
  deg_sample <- na.omit(deg_sample)
  viterbi_res <- duration_viterbi(deg_sample, pipar, tpmpar, od, params)
  viterbi_path <- viterbi_res[["path"]]
  c_site <- ifelse(any(viterbi_path == 3), 
                   max(which(viterbi_path == 3)), 
                   length(deg_sample))
  #get probabilities of signal peptide model
  prob.signal <- viterbi_res[["viterbi"]][c_site, viterbi_path[c_site]]
  #get probabilities of no signal peptide model
  prob.non <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site], 0)
  prob.total <- exp(prob.signal - prob.non)
  res <- list(sp_probability = unname(1 - 1/(1 + prob.total)), 
              sp_start = 1,
              sp_end = c_site,
              struc = viterbi_path,
              prot = toupper(prot[1L:70]),
              name = attr(prot, "name"),
              str_approx = 0)
  class(res) <- "hsmm_pred"
  
  #structure approximation - if atypical (normally negative signal peptide)
  while(!all(1L:4 %in% res[["struc"]])) {
    res[["struc"]] <- c(res[["struc"]], which.min(1L:4 %in% res[["struc"]]))
    res[["str_approx"]] <- res[["str_approx"]] + 1
  }
  
  res
}


#' GUI for signal.hsmm
#'
#' A graphical user interface for predicting presence of signal peptides/
#' @return null.
#' @export
#' @seealso \code{\link{run_signal.hsmm}}
#' @examples
#' \donttest{
#' gui_signal.hsmm()
#'}

gui_signal.hsmm <- function() {
  runApp(system.file("signal_gui", package = "signal.hsmm"))
}