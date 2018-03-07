constraints_slsqp_gausslegendre <- function(gams, rho, y, d, n.ints, alpha, n.nodes, natural){
  # This function computes inequality constraints.
  #
  # Inputs:
  # gams: set of gammas at which the coverage is
  #   required to be greater than or equal to 1 - alpha (vector)
  # rho: parameter (correlation)
  # y: contains knots values of the b and s functions
  # d: the b and s functions are optimized in the interval (0, d]
  # n.ints: number of intervals in (0, d]
  # c.alpha = quantile of the standard normal distribution
  #
  # Output:
  # A vector of inequality constraints
  #
  # Written by P.Kabaila in June 2008
  # Rewritten in R by R Mainzer, March 2017

  len.gams <- length(gams)
  covs <- rep(0, len.gams)

  c.alpha <- stats::qnorm(1 - alpha/2)

  b.spl <- spline_b(y, d, n.ints, c.alpha, natural)
  s.spl <- spline_s(y, d, n.ints, c.alpha, natural)

  for(i in 1:len.gams){
    covs[i] <- compute_cov_legendre(gams[i], rho, y, d, n.ints, alpha, n.nodes, b.spl, s.spl)
  }

  out <- covs - (1 - alpha)

}
