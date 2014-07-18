#' hsmm_pred class
#'
#' A single prediction of \code{signal.hsmm}.
#'
#' @details Always a named list of five elements
#' \enumerate{
#' \item \code{sp_probability} is a probability of signal peptide presence.
#' \item \code{sp_start} is a start of potential signal peptide (naively 1 aminoacid).
#' \item \code{sp_end} is a position of last amino acid of signal peptide.
#' \item \code{struc} is numeric vector representing predicted structure of input 
#' protein.
#' \item \code{prot} is character vector containing input sequence of amino acids.
#' }
#' @seealso \code{\link{summary.hsmm_pred}} \code{\link{plot.hsmm_pred}}
#' @name hsmm_pred
NULL

#' Plot single signal.hsmm prediction
#'
#' Plot objects of class \code{\link{hsmm_pred}}.
#'
#' @param object of class \code{\link{hsmm_pred}}.
#' @param if \code{TRUE}, legend is added to the plot.
#' @S3method plot hsmm_pred
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
    legend("left", col = c("blue", "red", "green", "orange", "black", "white"),
           lwd = c(5, 5, 5, 5, 2, 1), lty = c(rep("solid", 4), "dashed", "blank"),
           legend = c("n-region", "h-region", "c-region", "mature protein", 
                      "cleavage site", 
                      paste0("Signal peptide probability: ", 
                             round(x[["sp_probability"]], 2))), bty = "n")
}

#' Summarize single signal.hsmm prediction
#'
#' Summarizes objects of class \code{\link{hsmm_pred}}.
#'
#' @S3method summary hsmm_pred
#' @param object of class \code{\link{hsmm_pred}}.
#' @return Nothing.
#' @export

summary.hsmm_pred <- function(object, ...) {
  struc <- rle(x[["struc"]])[["lengths"]]
  cstruc <- cumsum(struc)
  cat(paste0("Probability of signal peptide presence: ", 
             object[["sp_probability"]], "\n"))
  cat(paste0("Start of signal peptide: ", object[["sp_start"]], "\n"))
  cat(paste0("End of signal peptide: ", object[["sp_end"]], "\n"))
  cat(paste0("n-region (length ", struc[1], "):\n"))
  cat(paste0(c("         ", object[["prot"]][1L:cstruc[1]]), collapse = ""))
  cat(paste0("\nh-region (length ", struc[2], "):\n"))
  cat(paste0(c("         ", object[["prot"]][(cstruc[1] + 1):cstruc[2]]), collapse = ""))
  cat(paste0("\nc-region (length ", struc[3], "):\n"))
  cat(paste0(c("         ", object[["prot"]][(cstruc[2] + 1):cstruc[3]]), collapse = ""))
}
