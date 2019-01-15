#' Compute the coverage probability of the CIUUPI
#'
#' Evaluate the coverage probability of the confidence interval that
#' utilizes uncertain prior information (CIUUPI) at \code{gam}.
#'
#' @param gam A value of gamma or vector of gamma values at which
#' the coverage probability function is evaluated
#' @param bsvec The vector (b(1),...,b(5),s(0),...,s(5)) that specifies the CIUUPI
#' @param a A vector used to specify the parameter of interest
#' @param c A vector used to specify the parameter about which
#' we have uncertain prior information
#' @param x The n by p design matrix
#' @param rho A known correlation
#' @param alpha The nominal coverage probability is 1 - \code{alpha}
#' @param natural Equal to 1 (default) if the b and s functions are obtained by
#' natural cubic spline interpolation or 0 if obtained by clamped cubic spline
#' interpolation
#'
#' @return The value(s) of the coverage probability of the CIUUPI at \code{gam}.
#'
#' @details
#' Suppose that \deqn{y = X \beta + \epsilon} where \eqn{y} is a random \eqn{n}-vector of
#' responses, \eqn{X} is a known \eqn{n} by \eqn{p} matrix with linearly
#' independent columns, \eqn{\beta} is an unknown parameter \eqn{p}-vector and
#' \eqn{\epsilon} is the random error with components that are iid normally distributed
#' with zero mean and known variance.
#' The parameter of interest is \eqn{\theta = } \code{a}' \eqn{\beta}. The uncertain
#' prior information is that \eqn{\tau = } \code{c}' \eqn{\beta} -
#' \code{t} = 0, where \code{a}
#' and \code{c} are specified linearly independent vectors and \code{t} is a specified number.
#'  \code{rho} is the known
#' correlation between the least squares estimators of \eqn{\theta} and \eqn{\tau}.
#' The user must specify either \code{a}, \code{c} and \code{x} or
#' \code{rho}.  If \code{a}, \code{c} and \code{x} are specified then
#' \code{rho} is computed.
#'
#' The CIUUPI is specified by the vector (b(1),...,b(5),s(0),...,s(5)), \code{alpha} and \code{natural}
#'
#' @seealso
#' \code{\link{ciuupi}}, \code{\link{bsciuupi}}
#'
#' @examples
#' alpha <- 0.05
#'
#' # Find the vector (b(1),b(2),...,b(5),s(0),s(1),...,s(5)) that specifies the
#' # CIUUPI: (this may take a few minutes to run)
#' \donttest{
#' bsvec <- bsciuupi(alpha, rho = 0.4)
#' }
#'
#' # The result (to 7 decimal places) is
#' bsvec <- c(0.129443483, 0.218926703, 0.125880945, 0.024672734, -0.001427343,
#'            1.792489585, 1.893870240, 2.081786492, 2.080407355,  1.986667246,
#'            1.958594824)
#'
#' # Graph the coverage probability function
#' gam <- seq(0, 10, by = 0.1)
#' cp <- cpciuupi(gam, bsvec, alpha, rho = 0.4)
#' plot(gam, cp, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
#' main = "Coverage Probability", col = "blue",
#' xlab = expression(paste("|", gamma, "|")), ylim = c(0.94999, 0.95001))
#' abline(h = 1-alpha, lty = 2)
#'
#' @export

cpciuupi <- function(gam, bsvec, alpha, natural = 1, rho = NULL, a, c, x){
  # Use this program to compute the coverage probability of the
  # new confidence interval.
  #
  # Written by R Mainzer, September 2017

  # Find rho
  if(is.null(rho)){

    # Do the QR decomposition of the X matrix and find X transpose X inverse
    qrstr <- qr(x)
    R <- qr.R(qrstr)
    XTXinv <- chol2inv(R)

    # Compute rho
    rho <- (t(a) %*% XTXinv %*% c) / sqrt( t(a) %*% XTXinv %*% a %*% t(c) %*% XTXinv %*% c)
    rho <- as.numeric(rho)

  }

  # The following inputs are needed here
  n.nodes <- 5
  d <- 6
  n.ints <- 6
  c.alpha <- stats::qnorm(1 - alpha/2)

  # Find the b and s functions
  b.spl <- spline_b(bsvec, d, n.ints, c.alpha, natural)
  s.spl <- spline_s(bsvec, d, n.ints, c.alpha, natural)

  # Compute the coverage probability
  res <- rep(0, length(gam))
  for(i in 1:length(gam)){
    res[i] <- compute_cov_legendre(gam[i], rho, bsvec, d, n.ints, alpha, n.nodes, b.spl, s.spl)
  }

  out <- res

}
