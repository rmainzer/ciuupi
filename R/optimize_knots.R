optimize_knots <- function(lambda, rho, alpha, gams,
                            d, n.ints, n.nodes, natural){
  # This function will find the value of the vector
  # (b(1),b(2),...,b(5),s(1),...,s(5) that specifies
  # the CIUUPI, given the optimized value of lambda.
  # This vector is found by numerical constrained
  # optimization.
  #
  # The coverage probability integrand is approximated by
  # Gauss-Legendre quadrature between the knots.
  #
  # Inputs:
  # gams = set of gammas at which the coverage is
  #   required to be greater than or equal to 1 - alpha
  # rho = a parameter of the model
  # lambda: used in specifying the objective function
  # d: The function b is not the same as for the standard
  #    1-alpha CI in the interval [-d,d]
  # n.ints = number of intervals in [0,d]
  # constr = 1 if cubic spline b has zero first derivative
  #          at -d and d.
  #        = 0 if the cubic spline b has no such constraint
  # alpha: the desired minimum coverage probability is
  #        1 - alpha
  # n.nodes = the number of nodes/weights used for gauss
  #   legendre quadrature
  #
  # Output:
  # A vector with the values at the knots of the b and s
  # functions
  #
  # Written by P.Kabaila in June 2008.
  # Rewritten in R by R Mainzer, March 2017

  # Specify the quantile of the standard normal distribution
  c.alpha <- stats::qnorm(1 - alpha/2)

  # Find a starting value for the optimization function
  start <- standard_CI(d, n.ints, alpha)

  # Specify lower and upper bounds on the vector of values
  # of the b and s functions evaluated at the knots
  low <- c(rep(-100, n.ints - 1), rep(0.5, n.ints))
  up <- c(rep(100, n.ints - 1), rep(200, n.ints))

  # Make the objective function a function of one argument, y
  obj_fun <- functional::Curry(objective, lambda = lambda, d = d,
                               n.ints = n.ints, n.nodes = n.nodes, alpha = alpha,
                               natural = natural)

  # Make the constraint function a function of one argument, y
  cons_fun <- functional::Curry(constraints_slsqp_gausslegendre, gams = gams, rho = rho,
                                d = d, n.ints = n.ints,
                                alpha = alpha, n.nodes = n.nodes, natural = natural)

  # Find the values of the knots using the optimization function
  res <- nloptr::slsqp(start, obj_fun, hin = cons_fun, lower = low,
               upper = up, nl.info = FALSE)
  new.par <- res$par

  # Output the vector with knot values which specifies the new
  # confidence interval
  out <- new.par

}

