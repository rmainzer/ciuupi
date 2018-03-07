integrand_sel <- function(x, gam, y, d, n.ints, alpha, s.spl){
  # This function evaluates
  # (s(|x|) - c_alpha) * phi(x-gamma)
  # for a vector x.
  #
  # Inputs:
  # x: vector at which sel integrand is to be evaluated
  # gam: parameter
  # y: contains knots values of the b and s functions
  # d: the b and s functions are optimized in the interval (0, d]
  # n.ints: number of intervals in (0, d]
  # c.alpha = quantile of the standard normal distribution
  #
  # Output:
  # A vector of values of the function with the same dimension
  # as x.
  #
  # Written by P.Kabaila in June 2008.
  # Rewritten in R by R Mainzer, March 2017

  c.alpha <- stats::qnorm(1 - alpha/2)
  tmp1 <- s.spl(x) - c.alpha
  tmp2 <- stats::dnorm(x - gam, 0, 1)
  res <- tmp1 * tmp2

}
