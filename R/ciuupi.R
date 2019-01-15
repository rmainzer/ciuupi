#' Compute the confidence interval that utilizes the uncertain prior information
#'
#' @param a A vector used to specify the parameter of interest
#' @param c A vector used to specify the parameter about which we have uncertain
#' prior information
#' @param t A number used to specify the parameter about which we have uncertain
#' prior information
#' @param x The n by p design matrix
#' @param bsvec The vector (b(1),...,b(5),s(0),...,s(5)) that specifies the CIUUPI
#' @param y The n-vector of observed responses
#' @param alpha 1 - \code{alpha} is the minimum coverage probability of the
#' confidence interval
#' @param natural Equal to 1 (default) if b and s functions are obtained by
#' natural cubic spline interpolation or 0 if obtained by clamped cubic spline
#' interpolation
#' @param sig Standard deviation of the random error.  If a value is not
#' specified then \code{sig} is estimated from the data.
#'
#' @return The confidence interval that utilizes uncertain prior information
#'
#' @details
#' Suppose that \deqn{y = X \beta + \epsilon} where \eqn{y} is a random \eqn{n}-vector of
#' responses, \eqn{X} is a known \eqn{n} by \eqn{p} matrix with linearly
#' independent columns, \eqn{\beta} is an unknown parameter \eqn{p}-vector and
#' \eqn{\epsilon} with components that are iid normally distributed
#' with zero mean and known variance.
#' Then \code{ciuupi} will compute a confidence interval for
#' \eqn{\theta=}\code{a}' \eqn{\beta} that utilizes the uncertain prior information that
#' \code{c}' \eqn{\beta} - \code{t} = 0, where \code{a} and \code{c} are specified
#' linearly independent vectors and \code{t} is a specified number.
#'
#' In the example below we use the data set described in Table 7.5
#' of Box et al. (1963).  A description of the parameter of interest and the
#' parameter about which we have uncertain prior information is given in
#' Dicsussion 5.8, p.3426 of Kabaila and Giri (2009).
#'
#' @examples
#' # Specify alpha, a, c, x
#' alpha <- 0.05
#' a <- c(0, 2, 0, -2)
#' c <- c(0, 0, 0, 1)
#' x1 <- c(-1, 1, -1, 1)
#' x2 <- c(-1, -1, 1, 1)
#' x <- cbind(rep(1, 4), x1, x2, x1*x2)
#'
#' # Find the vector (b(1),b(2),...,b(5),s(0),s(1),...,s(5)) that specifies the
#' # CIUUPI: (this may take a few minutes to run)
#' \donttest{
#' bsvec <- bsciuupi(alpha, a = a, c = c, x = x)
#' }
#'
#' # The result (to 7 decimal places) is
#' bsvec <- c(-0.03639701, -0.18051953, -0.25111411, -0.15830362, -0.04479113,
#'            1.71997203, 1.79147968, 2.03881195, 2.19926399, 2.11845381,
#'            2.00482563)
#'
#' # Specify t and y
#' t <- 0
#' y <- c(87.2, 88.4, 86.7, 89.2)
#'
#' # Find the CIUUPI
#' res <- ciuupi(alpha, a, c, x, bsvec, t, y, natural = 1, sig = 0.8)
#' res
#'
#' @references
#'
#' Box, G.E.P., Connor, L.R., Cousins, W.R., Davies, O.L., Hinsworth, F.R., Sillitto, G.P. (1963)
#' The Design and Analysis
#' of Industrial Experiments, 2nd edition, reprinted. Oliver and Boyd, London.
#'
#' Kabaila, P. and Giri, K. (2009) Confidence intervals in regression utilizing
#' prior information.  Journal of Statistical Planning and Inference, 139,
#' 3419 - 3429.
#'
#' @export

ciuupi <- function(alpha, a, c, x, bsvec, t, y, natural = 1, sig = NULL){

  # Do the QR decomposition of the X matrix and find X transpose X inverse
  qrstr <- qr(x)
  R <- qr.R(qrstr)
  XTXinv <- chol2inv(R)

  # Find beta hat, theta hat and tau hat
  beta.hat <- XTXinv %*% t(x) %*% y
  theta.hat <- as.numeric(t(a) %*% beta.hat)
  tau.hat <- as.numeric(t(c) %*% beta.hat - t)

  # If sigma is not specified, find an estimate of sigma
  if(is.null(sig)){
    cat("NOTE: Error variance not supplied by user.
        Error variance will be estimated from data.", "\n")

    if(dim(x)[1] - dim(x)[2] < 30){
      cat("WARNING: n - p < 30", "\n")
    }

    # If sig is missing, estimate it from the data
    n <- dim(x)[1]
    p <- dim(x)[2]
    sigsq <- (t(y - x %*% beta.hat) %*% (y - x %*% beta.hat)) / (n - p)
    sig <- as.numeric(sqrt(sigsq))

  }

  # Find gamma hat
  v.tau <- as.numeric(t(c) %*% XTXinv %*% c)
  gam.hat <- tau.hat / (sig * sqrt(v.tau))

  # Find variance of theta hat on sigma squared
  v.theta <- as.numeric(t(a) %*% XTXinv %*% a)

  # Compute the confidence interval
  bsfuns <- bsspline(gam.hat, bsvec, alpha, natural)

  new.ci <- theta.hat - sig * sqrt(v.theta) * bsfuns[, 2] +
             c(-1, 1) * sig * sqrt(v.theta) * bsfuns[, 3]

  data.frame(lower = new.ci[1], upper = new.ci[2],
                    row.names = c("ciuupi"))

}
