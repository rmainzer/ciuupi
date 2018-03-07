compute_cov_legendre <- function(gam, rho, y, d, n.ints, alpha, n.nodes, b.spl, s.spl){
  # This function computes the coverage probability of J(b, s).
  # The integral from (0, d) is broken down to integrals over
  # knots.  Each integral is computed using gauss legendre quadrature.
  # The number of nodes and weights for the approximation of
  # each integral can be changed.
  #
  # Input:
  # gam: parameter
  # rho: correlation
  # y: contains information about the knots of the b and s functions
  # d: b and s functions are constant after d
  # n.ints: number of intervals between 0 and d
  # c.alpha: quantile of the standard normal distribution
  # n.nodes: number of nodes/weights for each integral
  #
  # Written by R Mainzer, May 2017

  # Specify where the knots are locatated
  knots <- seq(0, d, by = d/n.ints)

  # Set up a vector to store the results
  int <- rep(0, length(knots))

  # Find the nodes and weights of the legendre quadrature
  quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")
  nodes <- quad.info$nodes
  weights <- quad.info$weights

  for(i in 1:d){
    # Specify bounds of the integral
    a <- knots[i]
    b <- knots[i+1]

    # Find the approximate integral
    adj.nodes <- ((b - a) / 2) * nodes + (a + b) / 2
    q <- integrand_cov(adj.nodes, gam, rho, y, d, n.ints, alpha, b.spl, s.spl)
    int[i] <- ((b - a) / 2) * sum(weights * q)
  }

  cp <- (1 - alpha) + sum(int)

}
