objective <- function(y, lambda, d, n.ints, n.nodes, alpha, natural){
  # This function computes the value of the objective function
  # for given functions b and s.
  # In other words, this function computes
  #
  # int_0^d (s(x) - c_alpha)(lambda + phi(x)) dx
  #
  # Inputs:
  # y: contains knots values of the b and s functions
  # lambda: used in specifying the objective function
  # d: the b and s functions are optimized in the interval (0, d]
  # n.ints: number of intervals in (0, d]
  # c.alpha = quantile of the standard normal distribution
  #
  # Output:
  # The objective function.
  #
  # Written by P.Kabaila in June 2008.
  # Rewritten in R by R Mainzer, March 2017

  c.alpha <- stats::qnorm(1 - alpha/2)

  s.spl <- spline_s(y, d, n.ints, c.alpha, natural)

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
    q <- integrand_obj(adj.nodes, y, lambda, d, n.ints, alpha, s.spl)
    int[i] <- ((b - a) / 2) * sum(weights * q)
  }

  out <- sum(int)

}
