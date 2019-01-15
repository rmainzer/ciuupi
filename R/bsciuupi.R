#' Compute the vector (b(1),...,b(5),s(0),...,s(5)) that specifies the CIUUPI
#'
#' Compute the vector (b(1),...,b(5),s(0),...,s(5)) that specifies
#' the confidence interval that utilizes uncertain prior
#' information (CIUUPI).
#'
#' @param a A vector used to specify the parameter of interest
#' @param c A vector used to specify the parameter about which
#' we have uncertain prior information
#' @param x The n by p design matrix
#' @param rho A known correlation
#' @param alpha The minimum coverage probability is 1 - alpha
#' @param natural Equal to 1 (default) if the functions b and s are found by
#' natural cubic spline interpolation or 0 if these functions are found by
#' clamped cubic spline interpolation in the interval [-6,6]
#'
#' @return The vector
#' \eqn{(b(1), b(2), \dots, b(5), s(0), s(1), \dots, s(5))} that specifies the CIUUPI.
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
#' The confidence interval for \eqn{\theta}, with minimum coverage probability
#'  \code{1 - alpha}, that utilizes the uncertain prior information that
#'  \eqn{\tau = } 0 belongs to a class of confidence intervals indexed
#'  by the functions b and s.
#' The function b is an odd continuous function and the function s is an even
#' continuous function. In addition, b(x)=0 and s(x) is equal to the
#' \eqn{1 - \alpha/2}
#' quantile of the standard normal distribution for all |x| greater than
#' or equal to 6. The values of these functions in the interval \eqn{[-6,6]}
#' are specified by  \eqn{b(1), b(2), \dots, b(5)} and
#' \eqn{s(0), s(1), \dots, s(5)} as follows. By assumption, \eqn{b(0)=0}
#' and \eqn{b(-i)=-b(i)}
#' and \eqn{s(-i)=s(i)} for \eqn{i=1,...,6}.
#' The values of \eqn{b(x)} and \eqn{s(x)} for any \eqn{x} in the interval \eqn{[-6,6]}
#' are found using cube spline interpolation for the given values of \eqn{b(i)}
#' and \eqn{s(i)} for \eqn{i=-6,-5,...,0,1,...,5,6}.
#'
#' The vector  \eqn{(b(1), b(2), \dots, b(5), s(0), s(1), \dots, s(5))}
#' is found by numerical constrained optimization
#' so that the confidence interval has minimum
#' coverage probability \code{1 - alpha} and utilizes the uncertain prior information
#' through its desirable expected length properties.
#' The optimization is performed using the \code{slsqp} function
#' in the \code{nloptr} package.
#'
#' @seealso
#' \code{\link{ciuupi}}
#'
#' @examples
#' # Compute the vector (b(1),...,b(5),s(0),...,s(5)) that specifies the CIUUPI,
#' # for given alpha and rho: (may take a few minutes to run)
#' \donttest{
#' bsvec <- bsciuupi(0.05, rho = 0.4)
#' }
#'
#' # The result (to 7 decimal places) is
#' bsvec <- c(0.129443483, 0.218926703, 0.125880945, 0.024672734, -0.001427343,
#'            1.792489585, 1.893870240, 2.081786492, 2.080407355,  1.986667246,
#'            1.958594824)
#' bsvec
#'
#' # Compute the vector (b(1),...,b(5),s(0),...,s(5)) that specifies the CIUUPI,
#' # for given alpha, a, c and x
#' x1 <- c(-1, 1, -1, 1)
#' x2 <- c(-1, -1, 1, 1)
#' x <- cbind(rep(1, 4), x1, x2, x1*x2)
#' a <- c(0, 2, 0, -2)
#' c <- c(0, 0, 0, 1)
#'
#' # The following may take a few minutes to run:
#' \donttest{
#' bsvec2 <- bsciuupi(0.05, a = a, c = c, x = x)
#' }
#'
#' # The result (to 7 decimal places) is
#' bsvec2 <- c(-0.03639701, -0.18051953, -0.25111411, -0.15830362, -0.04479113,
#'            1.71997203, 1.79147968, 2.03881195, 2.19926399, 2.11845381,
#'            2.00482563)
#' bsvec2
#'
#' @export

bsciuupi <- function(alpha, natural = 1, rho = NULL, a, c, x){

  # Specify the values of the inputs to other functions
  gams <- seq(0, 8, by = 0.05) # Constrain coverage probability at these values
  n.iter <- 5
  d <- 6
  n.ints <- 6
  n.nodes <- 5

  if(is.null(rho)){

    # Do the QR decomposition of the matrix X and find X transpose X inverse
    qrstr <- qr(x)
    R <- qr.R(qrstr)
    XTXinv <- chol2inv(R)

    # Compute rho
    rho <- (t(a) %*% XTXinv %*% c) /
      sqrt( t(a) %*% XTXinv %*% a %*% t(c) %*% XTXinv %*% c)
    rho <- as.numeric(rho)

  }

  cat("Computing the vector (b(1),...,b(5),s(0),...,s(5)) that specifies the CIUUPI... ")

  if(rho == 0){

    new.par <- standard_CI(d, n.ints, alpha)

  } else {

    lambda <- compute_lambda(rho, alpha, n.iter, d, n.ints, n.nodes, gams,
                            natural)

    new.par <- optimize_knots(lambda, rho, alpha, gams,
                              d, n.ints, n.nodes, natural)

  }

  cat("DONE", "\n")

  out <- new.par

}
