# SIGNAL-HSMM ------------------------------------

#' Predict presence of signal peptide in protein
#'
#' Using the hidden semi-Markov model predict presence of signal peptide in 
#' eukaryotic proteins.
#'
#' @param test_data single protein sequence (\code{character} vector) or list of 
#' sequences. It may be an object of class \code{\link[seqinr]{SeqFastaAA}}.
#' @return An object of class \code{hsmm_pred_list}.
#' @details Function \code{signalHsmm} returns respectively probability of presence of 
#' signal peptide, start of signal peptide and the probable cleavage site localization.
#' If input consists of more than one sequence, result is a data.frame where each column
#' contains above values for different proteins.
#' @note Currently start of signal peptide is naively set as 1 amino acid.
#' @useDynLib signalHsmm
#' @export
#' @seealso \code{\link{hsmm_pred_list}} \code{\link{hsmm_pred}} 
#' @keywords classif
#' @examples
#' #run signalHsmm on one sequence
#' x1 <- run_signalHsmm(benchmark_dat[[1]])
#' 
#' #run signalHsmm on one sequence, but input is a character vector
#' x2 <- run_signalHsmm(c("m", "a", "g", "k", "e", "v", "i", "f", 
#' "i", "m", "a", "l", "f", "i", "a", "v", "e", "s", "s", "p", "i", 
#' "f", "s", "f", "d", "d", "l", "v", "c", "p", "s", "v", "t", "s", 
#' "l", "r", "v", "n", "v", "e", "k", "n", "e", "c", "s", "t", "k", 
#' "k", "d", "c", "g", "r", "n", "l", "c", "c", "e", "n", "q", "n", 
#' "k", "i", "n", "v", "c", "v", "g", "g", "i", "m", "p", "l", "p", 
#' "k", "p", "n", "l", "d", "v", "n", "n", "i", "g", "g", "a", "v", 
#' "s", "e", "s", "v", "k", "q", "k", "r", "e", "t", "a", "e", "s", 
#' "l"))
#' 
#' #run signalHsmm on list of sequences
#' x3 <- run_signalHsmm(benchmark_dat[1:3])
#' #see summary of results
#' summary(x3)
#' #print results as data frame
#' pred2df(x3)
#' #summary one result
#' summary(x3[[1]])
#' plot(x3[[1]])

run_signalHsmm <- function(test_data) {
  predict.sighsmm_model(signalHsmm_main_model, test_data)
}

signalHsmm_decision <- function(prot, aa_group, pipar, tpmpar, 
                                od, overall_probs_log, params) {
  if (length(prot) == 1) {
    prot <- strsplit(prot, "")[[1]]
    if ("name" %in% names(attributes(prot)))
      attr(prot, "name") <- "undefined_name"
    if (length(prot) == 1)
      stop("Input sequence is too short.")
  }
  if(!is_protein(prot))
    stop("Atypical aminoacids detected, analysis cannot be performed.")

  deg_sample <- as.numeric(degenerate(tolower(prot)[1L:ifelse(length(prot) > 50, 50, length(prot))], aa_group))
  #remove atypical amino acids
  deg_sample <- na.omit(deg_sample)
  viterbi_res <- duration_viterbi(deg_sample - 1, pipar, tpmpar, od, params)
  viterbi_path <- viterbi_res[["path"]] + 1
  c_site <- ifelse(any(viterbi_path == 4), 
                   max(which(viterbi_path == 3)) + 1, 
                   length(deg_sample))
  #get probabilities of signal peptide model
  prob.signal <- viterbi_res[["viterbi"]][c_site, viterbi_path[c_site]]
  #get probabilities of no signal peptide model
  prob.non <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site], 0)
  prob.total <- exp(prob.signal - prob.non)
  res <- list(sp_probability = unname(1 - 1/(1 + prob.total)), #rescale(unname(1 - 1/(1 + prob.total))),
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


#' GUI for signalHsmm
#'
#' A graphical user interface for predicting presence of signal peptides.
#' @return null.
#' @export
#' @seealso \code{\link{run_signalHsmm}}
#' @note
#' Any ad-blocking software may be cause of malfunctions.

gui_signalHsmm <- function() {
  runApp(system.file("signal_gui", package = "signalHsmm"))
}

signalHsmm_main_model <- structure(list(aa_group = structure(list(`1` = c("r", "n", "d", 
                                                                          "q", "e", "h", "k"), `2` = c("g", "p", "s", "t", "y"), `3` = c("i", 
                                                                                                                                         "l", "m", "f", "w", "v"), `4` = c("a", "c")), .Names = c("1", 
                                                                                                                                                                                                  "2", "3", "4")), pipar = c(1, 0, 0, 0), 
                                        tpmpar = structure(c(0, 
                                                             0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0), .Dim = c(4L, 4L
                                                             )), 
                                        od = structure(c(0.282643251835592, 0.0130196824965156, 0.2081471295978, 
                                                         0.32818076798191, 0.286125198698055, 0.119080803616956, 0.411825369542798, 
                                                         0.311313005085296, 0.339187041102112, 0.704048679335078, 0.234788587143348, 
                                                         0.259545441021044, 0.0920445083642419, 0.16385083455145, 0.145238913716054, 
                                                         0.100960785911751), .Dim = c(4L, 4L)), 
                                        overall_probs_log = structure(c(-1.11419070051382, 
                                                                        -1.16695642571174, -1.34882348179846, -2.29302309583215), .Names = c("1", 
                                                                                                                                             "2", "3", "4")), 
                                        params = structure(c(0.0935251798561151, 0.187450039968026, 
                                                             0.211031175059952, 0.0823341326938449, 0.0975219824140687, 0.0819344524380496, 
                                                             0.0503597122302158, 0.0431654676258993, 0.0327737809752198, 0.0239808153477218, 
                                                             0.0215827338129496, 0.017585931254996, 0.0115907274180655, 0.00959232613908873, 
                                                             0.0103916866506795, 0.00559552358113509, 0.00479616306954436, 
                                                             0.00359712230215827, 0.00279776179056755, 0.00239808153477218, 
                                                             0.000799360511590727, 0.000799360511590727, 0.00199840127897682, 
                                                             0.00239808153477218, 0, 0, 0, 0, 0, 0, 0, 0, 0.00118811881188119, 
                                                             0.00118811881188119, 0.00277227722772277, 0.00158415841584158, 
                                                             0.00316831683168317, 0.0158415841584158, 0.0384158415841584, 
                                                             0.0601980198019802, 0.0994059405940594, 0.123960396039604, 0.15960396039604, 
                                                             0.131881188118812, 0.12990099009901, 0.0831683168316832, 0.0499009900990099, 
                                                             0.041980198019802, 0.0237623762376238, 0.0118811881188119, 0.0102970297029703, 
                                                             0.00475247524752475, 0.00158415841584158, 0.000792079207920792, 
                                                             0.00158415841584158, 0.000396039603960396, 0.000792079207920792, 
                                                             0, 0, 0, 0, 0, 0, 0, 0, 0.143879173290938, 0.156200317965024, 
                                                             0.284578696343402, 0.104928457869634, 0.171701112877583, 0.0600158982511924, 
                                                             0.0445151033386327, 0.0174880763116057, 0.010731319554849, 0.00278219395866455, 
                                                             0.00278219395866455, 0.000397456279809221, 0, 0, 0, 0, 0, 0, 
                                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.03125, 0.03125, 0.03125, 
                                                             0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 
                                                             0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 
                                                             0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 
                                                             0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 
                                                             0.03125), .Dim = c(32L, 4L), 
                                                           .Dimnames = list(NULL, c("n", "h", "c", "")))), 
                                   .Names = c("aa_group", "pipar", "tpmpar", "od", "overall_probs_log", "params"), class = "sighsmm_model")
