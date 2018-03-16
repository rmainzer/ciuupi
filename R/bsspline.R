#' Evaluate the functions b and s at x
#'
#' Evaluate the functions b and s, as specified by (b(1),b(2),...,b(5),s(0),s(1),...,s(5)),
#' \code{alpha} and \code{natural}, at \code{x}.
#'
#' @param x A value or vector of values at which the functions b and s
#' are to be evaluated
#' @param bsvec The vector (b(1),b(2),...,b(5),s(0),s(1),...,s(5))
#' @param alpha The minimum coverage probability is 1 - \code{alpha}
#' @param natural Equal to 1 (default) for natural cubic spline interpolation
#' or 0 for clamped cubic spline interpolation
#'
#' @return A data frame containing \code{x} and the corresponding values of the
#' functions b and s.
#'
#' @details The function b is an odd continuous function and the function s is an even
#' continuous function. In addition, b(x)=0 and s(x) is equal to the
#' \eqn{1 - \alpha/2}
#' quantile of the standard normal distribution for all |x| greater than
#' or equal to 6. The values of these functions in the interval \eqn{[-6,6]}
#' are specified by the vector \eqn{(b(1), b(2), \dots, b(5),s(0), s(1), \dots, s(5))}
#' as follows. By assumption, \eqn{b(0)=0}
#' and \eqn{b(-i)=-b(i)}
#' and \eqn{s(-i)=s(i)} for \eqn{i=1,...,6}.
#' The values of \eqn{b(x)} and \eqn{s(x)} for any \eqn{x} in the interval \eqn{[-6,6]}
#' are found using cube spline interpolation for the given values of \eqn{b(i)}
#' and \eqn{s(i)} for \eqn{i=-6,-5,...,0,1,...,5,6}.
#'
#'
#' The vector (b(1),b(2),...,b(5),s(0),s(1),...,s(5)) that specifies the confidence interval
#' that utilizes uncertain prior information (CIUUPI) is obtained using \code{\link{bsciuupi}}.
#'
#' @seealso
#' \code{\link{bsciuupi}}, \code{\link{ciuupi}}
#'
#' @examples
#' alpha <- 0.05
#'
#' # Find the vector (b(1),b(2),...,b(5),s(0),s(1),...,s(5)) that specifies the
#' # CIUUPI using the code:
#' # bsvec <- bsciuupi(alpha, rho = 0.4)
#' # This may take a few minutes to run
#' # The end result (to 7 decimal places) is
#' bsvec <- c(0.129443483, 0.218926703, 0.125880945, 0.024672734, -0.001427343,
#'            1.792489585, 1.893870240, 2.081786492, 2.080407355,  1.986667246,
#'            1.958594824)
#'
#' # Graph the functions b and s
#' x <- seq(0, 8, by = 0.1)
#' xseq <- seq(0, 6, by = 1)
#' bvec <- c(0, bsvec[1:5], 0)
#' quantile <- qnorm(1-(alpha)/2, 0, 1)
#' svec <- c(bsvec[6:11], quantile)
#' splineval <- bsspline(x, bsvec, alpha)
#'
#' plot(x, splineval[, 2], type = "l", main = "b function",
#' ylab = " ", las = 1, lwd = 2, xaxs = "i", col = "blue")
#' points(xseq, bvec, pch = 19, col = "blue")
#' plot(x, splineval[, 3], type = "l", main = "s function",
#' ylab = " ", las = 1, lwd = 2, xaxs = "i",  col = "blue")
#' points(xseq, svec, pch = 19, col = "blue")
#'
#' @export

bsspline <- function(x, bsvec, alpha, natural = 1){

  # Set input
  d <- 6
  n.ints <- 6
  c.alpha <- stats::qnorm(1 - alpha/2)

  # Find b and s functions
  sspl <- spline_s(bsvec, d, n.ints, c.alpha, natural)
  bspl <- spline_b(bsvec, d, n.ints, c.alpha, natural)

  x1 <- x[which(x <= -d)]
  x2 <- x[which(x > -d & x < d)]
  x3 <- x[which(x >= d)]

  bspl.res <- c(rep(0, length(x1)), bspl(x2), rep(0, length(x3)))
  sspl.res <- c(rep(c.alpha, length(x1)), sspl(x2), rep(c.alpha, length(x3)))

  out <- data.frame(x = x, b = bspl.res, s = sspl.res)

}
