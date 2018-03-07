#' Compute the scaled expected length of the CIUUPI
#'
#' Evaluate the scaled expected length of the confidence interval that
#' utilizes uncertain prior information (CIUUPI) at \code{gam}.
#'
#' @param gam A value of gamma or vector of gamma values at which
#' the scaled expected length function is evaluated
#' @param bsvec The vector (b(1),...,b(5),s(0),...,s(5)) that specifies the CIUUPI
#' @param a A vector used to specify the parameter of interest
#' @param c A vector used to specify the parameter about which
#' we have uncertain prior information
#' @param x The n by p design matrix
#' @param rho A known correlation
#' @param alpha The minimum coverage probability is 1 - \code{alpha}
#' @param natural Equal to 1 (default) if the functions b and s are obtained by natural cubic
#' spline interpolation or 0 if obtained by clamped cubic spline interpolation
#'
#' @return The value(s) of the scaled expected length at \code{gam}.
#'
#' @details
#'
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
#' The scaled expected length is defined as the expected length of the
#' CIUUPI divided
#' by the expected length of the standard confidence interval with the same minimum coverage probability.
#'
#' @seealso
#' \code{\link{ciuupi}}, \code{\link{bsciuupi}}
#'
#' @examples
#' \dontrun{
#' # Find the optimized knots of the b and s functions
#' alpha <- 0.05
#' bsvec <- bsciuupi(alpha, rho = 0.4)
#'
#' # Graph the scaled expected length function
#' gam <- seq(0, 8, by = 0.1)
#' sel <- selciuupi(gam, bsvec, alpha, rho = 0.4)
#' plot(gam, sel, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
#' main = "Scaled Expected Length", col = "blue",
#' xlab = expression(paste("|", gamma, "|")))
#' abline(h = 1, lty = 2)
#' }
#'
#' @export

selciuupi <- function(gam, bsvec, alpha, natural = 1, rho = NULL, a, c, x){
  # Use this program to compute the scaled expected length of the
  # new confidence interval.
  #
  # Written by R Mainzer, September 2017

  # Compute rho
  if(is.null(rho)){

    # Do the QR decomposition of the X matrix and find X transpose X inverse
    qrstr <- qr(x)
    R <- qr.R(qrstr)
    XTXinv <- solve(t(R) %*% R)

    # Compute rho
    rho <- (t(a) %*% XTXinv %*% c) / sqrt( t(a) %*% XTXinv %*% a %*% t(c) %*% XTXinv %*% c)
    rho <- as.numeric(rho)

  }

  # The following inputs are needed here
  d <- 6
  n.ints <- 6
  n.nodes <- 5

  # Specify the function s
  c.alpha <- stats::qnorm(1 - alpha/2)
  s.spl <- spline_s(bsvec, d, n.ints, c.alpha, natural)

  # Compute the scaled expected length
  res <- rep(0, length(gam))
  for(i in 1:length(gam)){
    res[i] <- compute_sel(gam[i], bsvec, d, n.ints, n.nodes, alpha, s.spl)
  }

  out <- res

}
