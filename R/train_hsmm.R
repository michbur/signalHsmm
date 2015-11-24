#' Train sighsmm_model object
#' 
#' @param train_data training data.
#' @param aa_group method of aggregating amino acids.
#' @param max_length maximum length of signal peptide.
#' @export
#' @return object of class \code{sighsmm_model}.

train_hsmm <- function(train_data, aa_group, max_length = 32) {
  ngroups <- length(aa_group)
  train_data <- lapply(train_data, toupper)
  ts <- calc_t(train_data, aa_group)
  
  t1 <- ts[["t1"]]
  t2 <- ts[["t2"]]
  t3 <- ts[["t3"]]
  t4 <- ts[["t4"]]
  
  overall <- t4 #table(degenerate(unlist(analized_sequences), aa5))
  overall.probs <- overall/sum(overall)          
  overall_probs_log = log(overall.probs) #for viterbi
  
  lengths <- ts[["lengths"]]
  params <- apply(lengths, 2, measure_region, max_length = max_length)
  params <- cbind(params, rep(1/max_length, max_length))
  
  #setting params for hmm -------
  additional_margin = 10
  pipar <- c(1,0,0,0)
  tpmpar <- matrix(c(0, 1, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1,
                     0, 0, 0, 0), 4, byrow = TRUE)
  od <- matrix(c((t1/sum(t1))[1L:ngroups],
                 (t2/sum(t2))[1L:ngroups],
                 (t3/sum(t3))[1L:ngroups],
                 (t4/sum(t4))[1L:ngroups]), 4, byrow = TRUE)
  
  res <- list(aa_group = aa_group, pipar = pipar, tpmpar = tpmpar, od = od, 
       overall_probs_log = overall_probs_log, params = params)
  class(res) <- "sighsmm_model"
  res
}

calc_t <- function(list_prots, aa_list) {
  nhc <- t(vapply(list_prots, find_nhc, rep(0, 4)))
  
  n_region <- NULL
  h_region <- NULL
  c_region <- NULL
  rest <- NULL
  
  for(i in 1L:length(list_prots)){
    region_starts <- nhc[i, ]
    n_region <- c(n_region, list_prots[[i]][1:(region_starts[2] - 1)])
    h_region <- c(h_region, list_prots[[i]][region_starts[2]:(region_starts[3] - 1)])
    c_region <- c(c_region, list_prots[[i]][region_starts[3]:(region_starts[4] - 1)])
    rest <- c(rest, list_prots[[i]][region_starts[4]:length(list_prots[[i]])])
  }
  
  t1 <- rep(0, length(aa_list))
  temp <- table(degenerate(tolower(n_region), aa_list))
  t1[as.numeric(names(temp))] <- temp
  names(t1) <- 1:length(aa_list)
  
  t2 <- rep(0, length(aa_list))
  temp <- table(degenerate(tolower(n_region), aa_list))
  t2[as.numeric(names(temp))] <- temp
  names(t2) <- 1:length(aa_list)
  
  t3 <- rep(0, length(aa_list))
  temp <- table(degenerate(tolower(n_region), aa_list))
  t3[as.numeric(names(temp))] <- temp
  names(t3) <- 1:length(aa_list)
  
  t4 <- rep(0, length(aa_list))
  temp <- table(degenerate(tolower(rest), aa_list))
  t4[as.numeric(names(temp))] <- temp
  names(t4) <- 1:length(aa_list)
  
  len_c <- nhc[, "cs"] - nhc[, "start_c"]
  len_h <- nhc[, "start_c"] - nhc[, "start_h"]
  len_n <- nhc[, "start_h"] - nhc[, "start_n"]
  lengths <- matrix(c(len_n, len_h, len_c), ncol = 3)
  colnames(lengths) <- c("n", "h", "c")
  
  list(mean_cs = mean(nhc[, 4]), sd_cs = sd(nhc[, 4]), t1 = t1, t2 = t2, t3 = t3, t4 = t4, 
       lengths = lengths)
}

measure_region <- function(region, max_length = 32) {
  lengths <- table(region)
  res <- rep(0, max_length)
  lengths <- lengths[as.numeric(names(lengths))>0] #removing lengths smaller than 1
  
  start_l <- min(as.numeric(names(lengths)))
  end_l <- max(as.numeric(names(lengths)))
  if(prod(start_l:end_l %in% as.numeric(names(lengths)))){
    max_length <- length(lengths) #if all lengths are present in training set
  } else{
    max_length <- 1
    sl <- sum(lengths)
    while(sum(lengths[1:max_length])/sl <= 0.51) {
      max_length <- which.min(start_l:end_l %in% as.numeric(names(lengths))) - 1
      start_l <- start_l + 1  #to assure that max_length wouldn't be too small
      max_length <- ifelse(max_length == 0, length(lengths), max_length)
    }
  }
  max_length <- min(max_length, max_length)
  
  prop_lengths <- lengths[1:max_length]/sum(lengths[1:max_length])
  res[as.numeric(names(prop_lengths))] <- prop_lengths
  res
}