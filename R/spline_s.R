spline_s <- function(y, d, n.ints, c.alpha, natural){
  # Return the value of the s function at x
  #
  # Input:
  # x: value at which to evaluate the s function
  # y: contains knot values of the b and s functions
  # d: b and s functions are constant after d
  # n.ints: number of intervals in (0, d]
  # c.alpha: quantile of the standard normal distribution
  #
  # Written by R Mainzer, March 2017

  s.knots <- seq(0, d, d/n.ints)
  s.vals <- c(y[n.ints:(2 * n.ints - 1)], c.alpha)

  s.knots.all <- seq(-d, d, d/n.ints)
  s.vals.all <- c(rev(s.vals), s.vals[2:(d+1)])

  if(natural == 1){
    s.spl <- stats::splinefun(s.knots.all, s.vals.all, method = "natural")
  } else {
    s.spl.pp <- pracma::cubicspline(s.knots.all, s.vals.all, endp2nd = TRUE)
    s.spl <- function(x) pracma::ppval(s.spl.pp, x)
  }

  out <- s.spl

}
