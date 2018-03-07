standard_CI <- function(d, n.ints, alpha){
  # This function calculates the vector of values
  # of b and s functions evaluated at the knots such
  # that y, d and n_ints specify the standard 1-alpha
  # confidence interval for theta.
  #
  # The main use of this function is to provide
  # a starting value for the optimization problem.
  #
  # Inputs:
  # d: the b and s functions are optimized in the interval (0, d]
  # n.ints: number of intervals in (0, d]
  # c.alpha = quantile of the standard normal distribution
  #
  # Output: the vector out.
  #
  # Written by P.Kabaila in June 2008.
  # Rewritten in R by R Mainzer, March 2017

  c.alpha <- stats::qnorm(1 - alpha/2)

  out <- c(rep(0, n.ints - 1), rep(c.alpha, n.ints))

}
