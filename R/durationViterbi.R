#' Compute most probable path with extended Viterbi algorithm
#' 
#' Viterbi algorithm for Hidden Markov Model with duration
#' Valid only for special case on data
#' 
#' @param aa_sample \code{character} vector representing single aminoacid sequence.
#' @param pipar Probabilities of initial state in Markov Model
#' @param tpmpar Matrix with transition probabilities between states
#' @param od Matrix of response probabilities. Eg. od[1,2] is a probability of signal 2 in state 1
#' @param params Matrix of probability distribution for duration. 
#'              Eg. params[10,2] is probability of duration of time 10 in state 2
#' @export
#' @return A list of length four:
#' \itemize{
#' \item{a vector of most probable path (\code{path},}
#' \item{values of probability in all intermediate points (\code{viterbi}),}
#' \item{matrix that for every signal and state gives the previous state in viterbi path (\code{psi}),}
#' \item{matrix that for every signal and state gives the duration in that state on viterbi path(\code{duration})}
#' }
#' @note Currently has very restricted application to specific input

duration_viterbi <- function(aa_sample, pipar, tpmpar, od, params){
  max.duration <- dim(params)[1]
  viterbi <- matrix(nrow=length(aa_sample), ncol = 4)
  psi <- matrix(nrow=length(aa_sample), ncol = 4)
  dura <- matrix(nrow=length(aa_sample), ncol = 4)
  for(j in 1L:4) {
    viterbi[1,j] <- log(pipar[j]) + log(od[j,aa_sample[1]])
    psi[1,j] <- 1
    dura[1,j] <- 1
  }  
  for(i in 2:length(aa_sample)){
    for(j in 1L:4){
      max <- -Inf
      max.i <- 0
      max.dur <- 0 
      for(k in 1L:4){
        for(d in 1L:min(max.duration, i)){
          if(i - d == 0){
            if(j == 1){
              previous <- 1
              transition <- 1
            }
            else{
              previous <- -Inf
              transition <- 1
            }
          }
          else{
            previous <- viterbi[i - d, k]
            transition <- log(tpmpar[k, j])
          }
          duration <- log(params[d, j])
          responses <- sum(log(od[j, aa_sample[(i - d + 1):i]]))
          if(previous + transition + duration + responses > max){
            max <- previous + transition + duration + responses
            max.i <- k
            max.dur <- d
          }
        }
        
      }
      viterbi[i, j] <- max
      psi[i, j] <- max.i
      dura[i, j] <- max.dur
    }  
  }
  path <- NULL
  path[length(aa_sample)] <- which.max(viterbi[length(aa_sample),])
  i <- length(aa_sample) - 1
  last <- length(aa_sample)
  while(i > 1){
    if(last - i < dura[last, path[last]]){
      path[i] <- path[last]
    }
    else{
      path[i] <- psi[last, path[last]]
      last <- i
    }
    i <- i -1
  }
  path[1] <- 1
  
  list(path = path,
       viterbi = viterbi,
       psi = psi,
       duration = dura)
}
