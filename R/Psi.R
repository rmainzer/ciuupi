Psi <- function(x, y, mu, variance){
  # This function calculates
  #
  # Psi(x,y,mu,variance) = P(x le Z le y)
  #                      = P(Z le y) - P(Z le x)
  #
  # where Z ~ N(mu,variance).
  #
  # Inputs:
  # x = number
  # y = number that is greater than or equal
  #      to x
  # mu = mean of normal distribution
  # variance = variance of normal distribution
  #
  # Output:
  # A single numerical value
  #
  # Written by P.Kabaila in May 2008.
  # Rewritten in R by R Mainzer in March 2017

  sigma <- sqrt(variance)
  term1 <- stats::pnorm(y, mean = mu, sd = sigma)
  term2 <- stats::pnorm(x, mean = mu, sd = sigma)
  out <- term1 - term2

}
