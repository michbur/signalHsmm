#' Adds k-mer hidden state to signalHsmm model
#' 
#' Changes parameters for Hidden Semi-Markov Model to
#' add k-mer
#' 
#' @param kMer \code{character} vector representing k-mer aminoacid sequence.
#' @param pipar Probabilities of initial state in Markov Model.
#' @param tpmpar Matrix with transition probabilities between states.
#' @param od Matrix of response probabilities. Eg. od[1,2] is a 
#' probability of signal 2 in state 1.
#' @param params Matrix of probability distribution for duration.
#' Eg. params[10,2] is probability of duration of time 10 in state 2.
#' @param pState number denoting hidden state right before k-mer.
#' @param nState number denoting hidden state right after k-mer.
#' @param pTrans Probability of change from pState to k-mer hidden state
#' @export
#' @return A list of length four:
#' \itemize{
#' \item{pipar}{ a vector of new probabilities of initial state in Markov Model,}
#' \item{tpmpar}{ a matrix with new transition probabilities between states,}
#' \item{od}{ matrix of new response probabilities,}
#' \item{params}{ matrix of new probability distributions for duration.}
#' }
#' @note Currently add only k-mers without distance.
add_k_mer_state <- function(kMer, pipar, tpmpar, od, params, pState, nState, pTrans){
  pipar2 <- c(pipar, rep(0, length(kMer)))
  nStates <- length(pipar)
  
  tpmpar2 <- cbind(tpmpar, matrix(0, nrow=nStates, ncol=length(kMer)))
  tpmpar2[pState,] <- (1-pTrans)*tpmpar2[pState,]
  tpmpar2[pState,nStates+1] <- pTrans
  
  od2 <- od
  params2 <- params
  for(i in 1:(length(kMer)-1)){
    sig <- unlist(kMer[i])
    #transition matrix
    transition <- c(rep(0, nStates+i), 1, rep(0, length(kMer)-i-1)) 
    tpmpar2 <- rbind(tpmpar2, transition)
    #response probability
    response <- rep(0, ncol(od))
    response[sig] <- 1/length(sig)
    od2 <- rbind(od2, response)
    #duration probability
    params2 <- cbind(params2, c(1, rep(0,nrow(params)-1)))
  }
  
  #last state

  sig <- unlist(tail(kMer,1))
  transition <- rep(0, ncol(tpmpar2))
  transition[nState] <- 1
  tpmpar2 <- rbind(tpmpar2, transition)
  #observation probs
  response <- rep(0, ncol(od))
  response[sig] <- 1/length(sig)
  od2 <- rbind(od2, response)
  #duration probability
  params2 <- cbind(params2, c(1, rep(0,nrow(params)-1)))
  
  return(list(pipar=pipar2,
              tpmpar=tpmpar2,
              od=od2,
              params=params2))
}
