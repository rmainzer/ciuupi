#' Compute the standard confidence interval
#'
#' @param a A vector used to specify the parameter of interest
#' @param x A known n by p matrix
#' @param y A known n-vector of responses
#' @param alpha 1 - \code{alpha} is the nominal coverage probability of the
#' confidence interval
#' @param sig Standard deviation of the random error.  If a value is not
#' specified, \code{sig} will be estimated from the data.
#'
#' @return The standard confidence interval
#'
#' @details
#'Suppose that \deqn{Y = X \beta + \epsilon} is a random \eqn{n}-vector of
#'responses, \eqn{X} is a known \eqn{n} by \eqn{p} matrix with linearly
#'independent columns, \eqn{\beta} is an unknown parameter \eqn{p}-vector and
#' \eqn{\epsilon} has a multivariate normal distribution with mean vector 0 and
#' variance \code{sig}^2 times the \eqn{n} by \eqn{n} identity matrix.
#' Then \code{cistandard} will compute the standard confidence interval for
#' \code{a}' \eqn{\beta}.
#'
#' In the example below we use the data set described in Table 7.5
#' of Box et al. (1963).  A description of the parameter of interest  is given
#' in Dicsussion 5.8, p.3426 of Kabaila and Giri (2009).
#'
#' @examples
#' y <- c(87.2, 88.4, 86.7, 89.2)
#' x1 <- c(-1, 1, -1, 1)
#' x2 <- c(-1, -1, 1, 1)
#' x <- cbind(rep(1, 4), x1, x2, x1*x2)
#' a <- c(0, 2, 0, -2)
#'
#' # Calculate the standard 95% confidence interval when sigma = 0.8
#' res <- cistandard(a, x, y, 0.05, sig = 0.8)
#' res
#'
#' @references
#' Box, G.E.P., Connor, L.R., Cousins, W.R., Davies, O.L., Hinsworth, F.R., Sillitto, G.P. (1963)
#' The Design and Analysis
#' of Industrial Experiments, 2nd edition, reprinted. Oliver and Boyd, London.
#'
#' Kabaila, P. and Giri, K. (2009) Confidence intervals in regression utilizing
#' prior information.  Journal of Statistical Planning and Inference, 139,
#' 3419 - 3429.
#'
#' @export

cistandard <- function(a, x, y, alpha, sig = NULL){

  # Do the QR decomposition of the X matrix and find X transpose X inverse
  qrstr <- qr(x)
  R <- qr.R(qrstr)
  XTXinv <- chol2inv(R)

  # Find beta hat, theta hat and tau hat
  beta.hat <- XTXinv %*% t(x) %*% y
  theta.hat <- as.numeric(t(a) %*% beta.hat)

  # If sigma is not sepcified, find an estimate of sigma
  if(is.null(sig)){
    cat("NOTE: Error variance not supplied by user.
        Error variance will be estimated from data", "\n")

    # If sig is missing, estimate it from the data
    n <- dim(x)[1]
    p <- dim(x)[2]
    sigsq <- (t(y - x %*% beta.hat) %*% (y - x %*% beta.hat)) / (n - p)
    sig <- as.numeric(sqrt(sigsq))

  }

  # Find variance of theta hat on sigma squared
  v.theta <- as.numeric(t(a) %*% XTXinv %*% a)

  # Find the standard confidence interval
  standard.ci <- theta.hat + c(-1, 1) * sig * sqrt(v.theta) * stats::qnorm(1 - alpha/2)

  data.frame(lower = standard.ci[1], upper = standard.ci[2])

}
