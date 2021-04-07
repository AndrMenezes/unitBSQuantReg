#' @importFrom extraDistr dfatigue pfatigue qfatigue rfatigue
#' @importFrom stats runif qnorm
#' @name dpqr-uBS
#' @aliases dpqr-uBS duBS puBS quBS ruBS
#'
#' @title The unit-Birnbaum-Saunders distribution
#'
#' @description Density function, distribution function, quantile function, random taumber generation function
#' for the unit-Birnbaum-Saunders distribution re-parameterized in terms of the \eqn{\tau}th quantile, \eqn{\tau \in (0,1)}.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu location parameter indicating the \eqn{\tau}-quantile \eqn{\tau \in (0, 1)}.
#' @param alpha nonnegative shape parameter.
#' @param tau the parameter to specify which quantile is to used.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#'
#' @return \code{duBS} gives the density, \code{puBS} gives the distribution function,
#' \code{quBS} gives the quantile function and \code{ruBS} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @author
#' Josmar Mazucheli
#'
#' André Felipe B. Menezes
#'
#' @references
#' Mazucheli, J., Menezes, A. F. B and Dey, S. (2018). The unit-Birnbaum-Saunders distribution with applications. \emph{Chilean Journal of Statistics}, \bold{9}(1), 47--57.
#'
#' @examples
#' set.seed(6969)
#' x <- ruBS(n = 1000, mu = 0.5, alpha = 1.5, tau = 0.5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by =  0.01)
#' hist(x, prob = TRUE, main = 'unit-Weibull distribution')
#' lines(S, duBS(x = S, mu = 0.5, alpha = 1.5, tau = 0.5), col = 2)
#' plot(ecdf(x))
#' lines(S, puBS(q = S, mu = 0.5, alpha = 1.5, tau = 0.5), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(quBS(p = S, mu = 0.5, alpha = 1.5, tau = 0.5), col = 2)
#'
#' @rdname dpqr-uBS
#' @export
#
duBS <- function (x, mu, alpha, tau = 0.5, log = FALSE)
{
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, alpha > 0, tau > 0, tau < 1)
  rqnorm <- qnorm(1 - tau);
  rqalpha <- (rqnorm * alpha) ^ 2;
  beta <- -log(mu) * (0.2e1 + rqalpha - rqnorm * alpha * sqrt(rqalpha + 0.4e1)) / 0.2e1;
  if(log)
  {
    dfatigue(x = -log(x), alpha, beta, log = TRUE) - log(x)
  }
  else
  {
    dfatigue(x = -log(x), alpha, beta, log = FALSE) / x
  }
}

#' @rdname dpqr-uBS
#' @export
#'
puBS <- function (q, mu, alpha, tau = 0.5, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(q > 0, q < 1, mu > 0, mu < 1, alpha > 0, tau > 0, tau < 1)
  rqnorm <- qnorm(1 - tau);
  rqalpha <- (rqnorm * alpha) ^ 2;
  beta <- -log(mu) * (0.2e1 + rqalpha - rqnorm * alpha * sqrt(rqalpha + 0.4e1)) / 0.2e1;
  if(lower.tail)
  {
    cdf <- pfatigue(-log(q), alpha, beta, lower.tail = FALSE, log.p = FALSE)
  }
  else
  {
    cdf <- pfatigue(-log(q), alpha, beta, lower.tail = TRUE, log.p = FALSE)
  }
  if(log.p) return(log(cdf)) else return(cdf)
}

#' @rdname dpqr-uBS
#' @export
#'
quBS <- function(p, mu, alpha, tau = 0.5, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(p > 0, p < 1, mu > 0, mu < 1, alpha > 0, tau > 0, tau < 1)
  rqnorm <- qnorm(1 - tau);
  rqalpha <- (rqnorm * alpha) ^ 2;
  beta <- -log(mu) * (0.2e1 + rqalpha - rqnorm * alpha * sqrt(rqalpha + 0.4e1)) / 0.2e1;
  if(lower.tail)
  {
    qtf <- exp(-qfatigue(1 - p, alpha, beta, lower.tail = TRUE, log.p = FALSE))
  }
  else
  {
    qtf <- exp(-qfatigue(1 - p, alpha, beta, lower.tail = FALSE, log.p = FALSE))
  }
  if(log.p) return(log(qtf)) else return(qtf)
}

#' @rdname dpqr-uBS
#' @export
#'
ruBS <- function(n, mu, alpha, tau = 0.5)
{
  stopifnot(n > 0, mu > 0, mu < 1, alpha > 0, tau > 0, tau < 1)
  u <- runif(n)
  quBS(p = u, mu = mu, alpha = alpha, tau = tau)
}
