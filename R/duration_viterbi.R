#' Compute most probable path with extended Viterbi algorithm
#' 
#' Viterbi algorithm for Hidden Markov Model with duration
#' Valid only for special case on data
#' 
#' @param aa_sample \code{character} vector representing single aminoacid sequence.
#' @param pipar Probabilities of initial state in Markov Model.
#' @param tpmpar Matrix with transition probabilities between states.
#' @param od Matrix of response probabilities. Eg. od[1,2] is a 
#' probability of signal 2 in state 1.
#' @param params Matrix of probability distribution for duration.
#' Eg. params[10,2] is probability of duration of time 10 in state 2.
#' @export
#' @return A list of length four:
#' \itemize{
#' \item{path}{ a vector of most probable path}
#' \item{viterbi}{ values of probability in all intermediate points,}
#' \item{psi}{ matrix that gives for every signal and state the previous state in 
#' viterbi path,}
#' \item{duration}{ matrix that gives for every signal and state gives the duration 
#' in that state on viterbi path.}
#' }
#' @note Currently has very restricted application to specific input All computations are on logarithms of probabilities

duration_viterbi <- function(aa_sample, pipar, tpmpar, od, params){
  max.duration <- dim(params)[1]
  nstates <- length(pipar)
  viterbi <- matrix(nrow=length(aa_sample), ncol = nstates) #probabiliy values for viterbi path that ends in specific state and signal
  psi <- matrix(nrow=length(aa_sample), ncol = nstates) #the previous state for viterbi path that ends in specific state and signal
  dura <- matrix(nrow=length(aa_sample), ncol = nstates) #the current duration for viterbi path that ends in specific state and signal
  #first signal is treated seperately
  for(j in 1L:nstates) {
    viterbi[1,j] <- log(pipar[j]) + log(od[j,aa_sample[1]])
    psi[1,j] <- 1
    dura[1,j] <- 1
  }  
  for(i in 2:length(aa_sample)){
    #For each state we will compute the probability of the viterbi path that ends in that state
    for(j in 1L:nstates){
      max <- -Inf
      max.i <- 0 #previous state that is maximising probability
      max.dur <- 0  #duration that is maximising probability
      for(k in 1L:nstates){
        for(d in 1L:min(max.duration, i)){
          # for every possible previous state, and for every possible duration in current state
          if(i - d == 0){ #if duration is as long as number of signal considered
            if(j == 1){ #only first state is accepted
              previous <- 1
              transition <- 1
            }
            else{ #other states are discarded with log-probability -Inf
              previous <- -Inf
              transition <- 1
            }
          }
          else{ #there is some previous state
            previous <- viterbi[i - d, k] #previous probability on this viterbi path
            transition <- log(tpmpar[k, j]) #probability of transition to current state from prevois state
          }
          duration <- log(params[d, j]) #probability of duration that lasts time d
          responses <- sum(log(od[j, aa_sample[(i - d + 1):i]])) #probability of generating d signal in current state
          if(previous + transition + duration + responses > max){ 
            #if that path is better than the best yet found store it
            max <- previous + transition + duration + responses
            max.i <- k
            max.dur <- d
          }
        }
        
      }
      #assign information about the best path that ends on signal i and state j
      viterbi[i, j] <- max 
      psi[i, j] <- max.i
      dura[i, j] <- max.dur
    }  
  }
  #now we extract information about the path. We look for the most probable path
  path <- NULL
  #the last state
  path[length(aa_sample)] <- which.max(viterbi[length(aa_sample),])
  i <- length(aa_sample) - 1
  last <- length(aa_sample)
  while(i > 1){
    if(last - i < dura[last, path[last]]){ #we are still in the same state
      path[i] <- path[last]
    }
    else{ #time to change the state for the previous, stored in matrix psi
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
