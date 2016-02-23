#use cut-off to calculate new probability
rescale <- function(x)
  ifelse(x > 0.025, (x - 0.025)/0.975 * 0.5 + 0.5, x*20)
