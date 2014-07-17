# SIGNAL-HSMM ------------------------------------

#' Predict presence of signal peptide in protein
#'
#' Using the hidden semi-Markov model predict presence of signal peptide in 
#' eukaryotic proteins.
#'
#' @param single protein sequence or list of sequences. May be an object of class 
#' \code{\link[seqinr]{SeqFastaAA}}
#' @return A vector or data frame (see datails).
#' @details Function \code{signal.hsmm} returns respectively probability of presence of 
#' signal peptide, start of signal peptide and the probable cleavage site localization.
#' If input consists from more than one sequence, result is a data.frame where each column
#' contains above values for different proteins.
#' @note Currently start of signal peptide is naively set as 1 amino acid.
#' @export

run.signal.hsmm <- function(test_data) {
  signal.hsmm_model <- structure(list(pipar = c(1, 0, 0, 0), 
                                      tpmpar = structure(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0), .Dim = c(4L, 4L)), 
                                      od = structure(c(0.203181196494605, 0.00235316265060241, 
                                                       0.069742553064986, 0.122654208889416, 0.360670685375614, 0.72484469126506, 
                                                       0.284975325524704, 0.288121595016022, 0.190922445671445, 0.0822430346385542, 
                                                       0.299423271300315, 0.230261863691714, 0.245225672458335, 0.190559111445783, 
                                                       0.345858850109995, 0.358962332402848), .Dim = c(4L, 4L)), 
                                      overall_probs_log = structure(c(-2.09838619260539, 
                                                                      -1.24437268303286, -1.4685380799115, -1.02453781966768), .Names = c("1", 
                                                                                                                                          "2", "3", "4")), 
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
                                                           0.03125), .Dim = c(32L, 4L), .Dimnames = list(NULL, c("n", "h", 
                                                                                                                 "c", "")))), 
                                 .Names = c("pipar", "tpmpar", "od", "overall_probs_log", 
                                            "params"))
  
  if(class(test_data) == "SeqFastaAA" || 
       class(test_data) == "character") {
    decisions <- signal.hsmm_decision(test_data, aa_group = aaaggregation, 
                                      pipar = signal.hsmm_model[["pipar"]], 
                                      tpmpar = signal.hsmm_model[["tpmpar"]], 
                                      od = signal.hsmm_model[["od"]], 
                                      overall_probs_log = signal.hsmm_model[["overall_probs_log"]], 
                                      params = signal.hsmm_model[["params"]])
  } else {
    decisions <- lapply(test_data, function(prot)
      signal.hsmm_decision(prot, aa_group = aaaggregation, 
                           pipar = signal.hsmm_model[["pipar"]], 
                           tpmpar = signal.hsmm_model[["tpmpar"]], 
                           od = signal.hsmm_model[["od"]], 
                           overall_probs_log = signal.hsmm_model[["overall_probs_log"]], 
                           params = signal.hsmm_model[["params"]]))
  }
  
  decisions
}

signal.hsmm_decision <- function(prot, aa_group, pipar, tpmpar, 
                                 od, overall_probs_log, params) {
  prot <- toupper(prot)[1L:50]
  deg_sample <- as.numeric(degenerate(prot, aa_group))
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
  prob.non <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1:c_site], 0) 
  res <- list(signal.peptide = 1 - 1/(1 + prob.signal), 
              sig.start = 1,
              sig.end = c_site,
              struc = viterbi_path,
              prot = prot)
  class(res) <- "hsmm_pred"
  res
}


#' Plot signal.hsmm predictions
#'
#' Plot objects of class \code{\link[signal.hsmm]{hsmm_pred}}.
#'
#' @param object of class \code{\link[signal.hsmm]{hsmm_pred}}.
#' @param if \code{TRUE}, legend is added to the plot.
#' @return Nothing.
#' @export

plot.hsmm_pred <- function(x, add_legend = TRUE) {
  plot(c(1, 50), c(1, 5), type="n", axes=F, ylab = "", xlab = "Amino acid index")
  axis(1, 1L:50, labels = FALSE)
  axis(1, 1L:25*2 - 1, labels = 1L:25*2 - 1)
  text(1L:50, 1, x[["prot"]])
  #get borders of regions, add 0.5 to have countinous regions
  struc <- cumsum(rle(x[["struc"]])[["lengths"]]) + 0.5
  lines(c(1, struc[1] + 0.5), c(1.5, 1.5), col = "blue", lwd = 5)
  lines(c(struc[1], struc[2]), c(1.5, 1.5), col = "red", lwd = 5)
  lines(c(struc[2], struc[3]), c(1.5, 1.5), col = "green", lwd = 5)
  lines(c(struc[3], 50), c(2.5, 2.5), col = "orange", lwd = 5)
  lines(c(struc[3], struc[3]), c(1.5, 2.5), lty = "dashed", lwd = 2)
  if (add_legend)
    legend("left", col = c("blue", "red", "green", "orange", "black"),
           lwd = c(5, 5, 5, 5, 2), lty = c(rep("solid", 4), "dashed"),
           legend = c("n-region", "h-region", "c-region", "mature protein", 
                      "cleavage site"), bty = "n")
}

#' Summarize signal.hsmm predictions
#'
#' Summarizes objects of class \code{\link[signal.hsmm]{hsmm_pred}}.
#'
#' @param object of class \code{\link[signal.hsmm]{hsmm_pred}}.
#' @return Nothing.
#' @export

summary.hsmm_pred <- function(object, ...) {
  struc <- rle(x[["struc"]])[["lengths"]]
  cstruc <- cumsum(struc)
  cat(paste0("Probability of signal peptide presence: ", 
             object[["signal.peptide"]], "\n"))
  cat(paste0("Start of signal peptide: ", object[["sig.start"]], "\n"))
  cat(paste0("End of signal peptide: ", object[["sig.end"]], "\n"))
  cat(paste0("n-region (length ", struc[1], "):\n"))
  cat(paste0(c("         ", object[["prot"]][1L:cstruc[1]]), collapse = ""))
  cat(paste0("\nh-region (length ", struc[2], "):\n"))
  cat(paste0(c("         ", object[["prot"]][(cstruc[1] + 1):cstruc[2]]), collapse = ""))
  cat(paste0("\nc-region (length ", struc[3], "):\n"))
  cat(paste0(c("         ", object[["prot"]][(cstruc[2] + 1):cstruc[3]]), collapse = ""))
}

