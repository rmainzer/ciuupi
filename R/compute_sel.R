compute_sel <- function(gam, y, d, n.ints, n.nodes, alpha, s.spl){
  # This function computes the value of the scaled expected
  # length for given functions b and s.
  # In other words, this function computes
  #
  # 1 + (1/c_alpha) int_{-d}^d (s(|x|) - c_alpha) phi(x-gamma) dx
  #
  # Inputs:
  # gam: parameter
  # y: contains knots values of the b and s functions
  # d: the b and s functions are optimized in the interval (0, d]
  # n.ints: number of intervals in (0, d]
  # c.alpha = quantile of the standard normal distribution
  #
  # Output:
  # The scaled expected length for given functions b and s.
  #
  # Written by P.Kabaila in June 2008.
  # Rewritten in R by R Mainzer, March 2017

  c.alpha <- stats::qnorm(1 - alpha/2)

  # Specify where the knots are locatated
  knots <- seq(-d, d, by = d/n.ints)

  # Set up a vector to store the results
  int <- rep(0, length(knots))

  # Find the nodes and weights of the legendre quadrature
  quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")
  nodes <- quad.info$nodes
  weights <- quad.info$weights

  for(i in 1:(length(knots) - 1)){
    # Specify bounds of the integral
    a <- knots[i]
    b <- knots[i+1]

    # Find the approximate integral
    adj.nodes <- ((b - a) / 2) * nodes + (a + b) / 2
    q <- integrand_sel(adj.nodes, gam, y, d, n.ints, alpha, s.spl)
    int[i] <- ((b - a) / 2) * sum(weights * q)
  }

  out <- 1 + (sum(int) / c.alpha)

}
