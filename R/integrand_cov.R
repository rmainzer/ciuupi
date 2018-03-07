integrand_cov <- function(x, gam, rho, y, d, n.ints, alpha, b.spl, s.spl){
  # This function evaluates
  # (k(x, gam, rho) - k_dag(x, gam, rho)) * phi(x - gam)
  # for a vector x.
  #
  # Inputs:
  # x: vector at which cp integrand is to be evaluated
  # gam: parameter
  # rho: parameter (correlation)
  # y: contains knots values of the b and s functions
  # d: the b and s functions are optimized in the interval (0, d]
  # n.ints: number of intervals in (0, d]
  # c.alpha = quantile of the standard normal distribution
  #
  # Output:
  # A vector of values of the function with the same
  # dimension as x.
  #
  # Written by P.Kabaila in June 2008
  # Rewritten in R by R Mainzer in March 2017

  c.alpha <- stats::qnorm(1 - alpha/2)

  mu1 <- rho * (x - gam)
  var <- 1 - rho^2
  k.dag1 <- Psi(-c.alpha, c.alpha, mu1, var)

  term.a1 <- b.spl(x)
  term.b1 <- s.spl(x)

  lh <- term.a1 - term.b1
  uh <- term.a1 + term.b1

  k1 <- Psi(lh, uh, mu1, var)
  term1 <- stats::dnorm(x - gam, 0, 1)

  mu2 <- rho * (-x - gam)
  k.dag2 <- Psi(-c.alpha, c.alpha, mu2, var)

  term.a2 <- b.spl(-x)
  term.b2 <- s.spl(-x)

  lh2 <- term.a2 - term.b2
  uh2 <- term.a2 + term.b2

  k2 <- Psi(lh2, uh2, mu2, var)
  term2 <- stats::dnorm(x + gam, 0, 1)

  res <- (k1 - k.dag1) * term1 + (k2 - k.dag2) * term2

}
