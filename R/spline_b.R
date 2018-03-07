spline_b <- function(y, d, n.ints, c.alpha, natural){
  # Return the value of the b function at x
  #
  # Input:
  # x: value at which to evaluate the b function
  # y: contains knot values of the b and s functions
  # d: b and s functions are constant after d
  # n.ints: number of intervals in (0, d]
  # c.alpha: quantile of the standard normal distribution
  #
  # Written by R Mainzer, March 2017

  b.knots <- seq(0, d, by = d/n.ints)
  y.rev <- rev(y[1:n.ints - 1])
  b.vals <- c(0, y[1:n.ints - 1], 0)

  b.knots.all <- seq(-d, d, by = d/n.ints)
  b.vals.all <- c(0, -y.rev, b.vals)

  # If natural = 1 use natural cubic spline, otherwise use clamped cubic
  # spline
  if(natural == 1){
    b.spl <- stats::splinefun(b.knots.all, b.vals.all, method = "natural")
  } else {
    b.spl.pp <- pracma::cubicspline(b.knots.all, b.vals.all, endp2nd = TRUE)
    b.spl <- function(x) pracma::ppval(b.spl.pp, x)
  }

  out <- b.spl

}
