compute_ratio_minus1 <- function(lambda, rho, alpha, gams,
                            d, n.ints, n.nodes, natural){
  # This function computes the ratio (expected gain / maximum
  # potential loss) - 1.  Another program can then us this
  # program to find the value of lambda which makes
  # (this ratio - 1) = 0.
  #
  # The coverage probability integrand is approximated by
  # gauss-legendre quadrature between the knots.
  #
  # Inputs:
  # gams = set of gammas at which the coverage is
  #   required to be greater than or equal to 1 - alpha
  # rho = a parameter of the model
  # lambda: used in specifying the objective function
  # d: The function b is not the same as for the standard
  #    1-alpha CI in the interval [-d,d]
  # n.ints = number of intervals in [0,d]
  # alpha: the desired minimum coverage probability is
  #        1 - alpha
  # n.nodes = the number of nodes/weights used for gauss
  #   legendre quadrature
  #
  # Output:
  # A vector that specifies the new confidence interval.
  # Plots of the b and s functions, and of the coverage
  # probability and scaled expected length.
  #
  # Written by P.Kabaila in June 2008.
  # Rewritten in R by R Mainzer, March 2017

  c.alpha <- stats::qnorm(1 - alpha/2)

  new.par <- optimize_knots(lambda, rho, alpha, gams,
                             d, n.ints, n.nodes, natural)

  s.spl <- spline_s(new.par, d, n.ints, c.alpha, natural)

  # Compute the required ratio
  sel.max <- stats::optimize(compute_sel, c(0, d), maximum = TRUE,
                      y = new.par, d = d, n.ints = n.ints,
                      n.nodes = n.nodes, alpha = alpha, s.spl = s.spl)$objective
  sel.min <- compute_sel(gam = 0, new.par, d, n.ints, n.nodes, alpha, s.spl = s.spl)

  expected.gain <- 1 - sel.min^2
  max.potential.loss <- sel.max^2 - 1

  # Output the required ratio minus 1
  out <- expected.gain / max.potential.loss - 1

}
