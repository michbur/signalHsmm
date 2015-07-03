#use cut-off to calculate new probability
rescale <- function(x)
  ifelse(x > 0.05, (x - 0.05)/0.95 * 0.5 + 0.5, x*10)
