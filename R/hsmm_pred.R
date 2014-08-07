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
#' @param x object of class \code{\link{hsmm_pred}}.
#' @param add_legend logical, if \code{TRUE}, legend is added to the plot.
#' @param ... ignored.
#' @return Nothing.
#' @export
#' @keywords hplot methods

plot.hsmm_pred <- function(x, add_legend = TRUE, ...) {
  plot(c(1, 50), c(1, 5), type="n", axes=F, ylab = "", xlab = "Amino acid index",
       main = x[["name"]])
  axis(1, 1L:50, labels = FALSE)
  axis(1, 1L:25*2 - 1, labels = 1L:25*2 - 1)
  
  #get borders of regions, add 0.5 to have countinous regions for purpose of easy plotting
  cstruc <- cumsum(rle(x[["struc"]])[["lengths"]])
  cstruc05 <- c(1, cstruc + 0.5)
  cstruc <- c(0, cstruc)
  sig_colours <- c("blue", "red", "green", "orange")
  
  #old boring black text
  #text(1L:50, 1, x[["prot"]])
  # 
  
  for(i in 1L:4) {
    #print amino acids in color!
    text((cstruc[i] + 1):cstruc[i + 1], 1, 
         x[["prot"]][(cstruc[i] + 1):cstruc[i + 1]], 
         col = sig_colours[i])
    lines(c(cstruc05[i], cstruc05[i + 1]), c(1.5, 1.5) + ifelse(i == 4, 0, 1), 
          col = sig_colours[i], lwd = 5)
  }
  
  lines(c(cstruc05[4], cstruc05[4]), c(1.5, 2.5), lty = "dashed", lwd = 2)
  if (add_legend)
    legend("topright", 
           col = rev(c(sig_colours, "black", "white", "white")),
           lwd = rev(c(5, 5, 5, 5, 2, 1, 1)), 
           lty = rev(c(rep("solid", 4), "dashed", "blank", "blank")),
           legend = rev(c("n-region", 
                          "h-region", 
                          "c-region", 
                          "mature protein", 
                          "cleavage site", 
                          paste0("Signal peptide probability: ", 
                                 signif(x[["sp_probability"]], digits = 2)),
                          ifelse(x[["str_approx"]] > 0, 
                                 "Signal peptide structure interpolated",
                                 " "))), 
           bty = "n")
}

#' Summarize single signal.hsmm prediction
#'
#' Summarizes objects of class \code{\link{hsmm_pred}}.
#'
#' @param object of class \code{\link{hsmm_pred}}.
#' @param ... ignored
#' @return Nothing.
#' @export
#' @keywords print methods

summary.hsmm_pred <- function(object, ...) {
  struc <- rle(object[["struc"]])[["lengths"]]
  cstruc <- cumsum(struc)
  cat(paste0(object[["name"]], "\n"))
  cat(paste0("Probability of signal peptide presence: ", 
             format(object[["sp_probability"]], digits = 4), "\n"))
  cat(paste0("Start of signal peptide: ", object[["sp_start"]], "\n"))
  cat(paste0("End of signal peptide: ", object[["sp_end"]], "\n"))
  cat(paste0("n-region (length ", struc[1], "):\n"))
  cat(paste0(c("         ", object[["prot"]][1L:cstruc[1]]), collapse = ""))
  cat(paste0("\nh-region (length ", struc[2], "):\n"))
  cat(paste0(c("         ", object[["prot"]][(cstruc[1] + 1):cstruc[2]]), collapse = ""))
  cat(paste0("\nc-region (length ", struc[3], "):\n"))
  cat(paste0(c("         ", object[["prot"]][(cstruc[2] + 1):cstruc[3]], "\n"), collapse = ""))
  if(object[["str_approx"]] > 0)
    cat("Signal peptide structure interpolated.\n")
}
