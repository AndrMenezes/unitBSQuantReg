#' @title Half-Normal Probability Plots with Simulated Envelopes
#'
#' @description Generates (half-)normal probability plots for a \code{\link{unitBSQuantReg}} object.
#'
#' @param object fitted model object of class \code{\link{unitBSQuantReg}}.
#' @param nsim number of simulations used to compute envelope. Default is 99 replications.
#' @param resid.type type of residuals to be used: The currently options are \code{cox-snell} or \code{quantile}.
#' @param halfnormal logical. If \code{TRUE}, a half-normal probability plot is produced. If \code{FALSE}, a normal probability plot is produced.
#' @param plot Should the (half-)normal probability plot be plotted? Default is \code{TRUE}.
#' @param level confidence level of the simulated envelope. Default is 95\%.
#' @param ... currently not used.
#'
#' @references
#' Atkinson, A. C., (1981). Two graphical displays for outlying and influential observations in regression. \emph{Biometrika} \bold{68}(1), 13--20.
#'
#' @author
#' Josmar Mazucheli
#'
#' André Felipe B. Menezes
#'
#'
#' @examples
#' set.seed(123)
#' x <- runif(100, -1, 1)
#' eta <- 1 + 2 * x
#' mu <- exp(eta) / (1 + exp(eta))
#' data <- data.frame(y = ruBS(n = 100, mu, alpha = 1.0, tau = 0.5), x)
#' fit <- unitBSQuantReg(y ~ x, data = data, tau = 0.5)
#' hnp(fit, resid.type = "cox-snell")
#' x1 <- runif(100, -1, 1); x2 <- runif(100, -1, 1)
#' eta <- 1 + 2 * x1 + 3 * x2
#' mu <- exp(eta) / (1 + exp(eta))
#' data <- data.frame(y = ruBS(n = 100, mu, alpha = 1.0, tau = 0.5), x1, x2)
#' fit <- unitBSQuantReg(y ~ ., data = data, tau = 0.5)
#' hnp(fit, resid.type = "quantile")
#'
#' @importFrom stats qnorm qexp residuals quantile median
#' @importFrom graphics plot lines
#'
#' @name hnp
NULL

#' @rdname hnp
#' @export
hnp <- function(object, ...)
{
  UseMethod("hnp", object)
}

#' @rdname hnp
#' @export
hnp.unitBSQuantReg <- function(object, nsim = 99, halfnormal = TRUE, plot = TRUE, level = 0.95, resid.type = c("cox-snell", "quantile"), ...)
{
  mu   <- object$fitted.values$mu
  alpha <- (object$fitted.values$alpha)
  tau  <- object$tau
  init <- object$coefficients
  control <- object$control
  data     <- object$df
  formula  <- object$formula
  link     <- object$link$name
  link.alpha <- object$link.alpha$name
  n        <- nrow(data)
  id       <- as.character(formula)[2]

  ysim    <- matrix(data = ruBS(n = n * nsim, mu = mu, alpha = alpha, tau = tau), nrow = n, ncol = nsim)

  res_sim <- sapply(1:nsim, function(j)
  {
    df      <- data
    df[,id] <- ysim[,j]
    obj     <- tryCatch(expr = unitBSQuantReg(
      formula = formula,
      tau = tau,
      link = link,
      link.alpha = link.alpha,
      data = df,
      start = init,
      control = control
      ),
      error = function(e) NULL)
    if(!is.null(obj))
    {
      residuals(obj, type = resid.type)
    }
    else {
      rep(NaN, n)
    }
  })

  is.na(res_sim) <- sapply(res_sim, is.infinite)

  if (halfnormal)
  {
    res_obs <- sort(abs(residuals(object, type = resid.type)))
    res_sim <- apply(res_sim, 2, function(x) sort(abs(x), na.last = TRUE))
    if (resid.type == "quantile")
    {
      res_teo <- qnorm((1:n + n - 1/8) / (2 * n + 0.5))
    }
    if (resid.type == "cox-snell")
    {
      res_teo <- qexp((1:n + n - 1/8) / (2 * n + 0.5))
    }
  }
  else
  {
    res_obs <- sort(residuals(object, type = resid.type))
    res_sim <- apply(res_sim, 2, function(x) sort(x, na.last = TRUE))
    if (resid.type == "quantile")
    {
      res_teo <- qnorm((1:n - 3 / 8) / (n + 1 / 4))
    }
    if (resid.type == "cox-snell") {
      res_teo <- qexp((1:n - 3 / 8) / (n + 1 / 4))
    }
  }

  alpha   <- (1 - level)/2
  res_lwr <- apply(res_sim, 1, quantile, probs = alpha, na.rm = T)
  res_upr <- apply(res_sim, 1, quantile, probs = 1 - alpha, na.rm = T)
  res_mid <- apply(res_sim, 1, median, na.rm = T)

  if (plot)
  {
    Ry <- c(min(res_lwr), max(res_upr))
    Rx <- range(res_teo)

    plot(
      x = res_teo,
      y = res_obs,
      xlab = 'Theoretical quantiles',
      ylab = 'Residuals',
      xlim = Rx,
      ylim = Ry,
      bty = 'o',
      pch = 3
    )
    lines(x = res_teo, y = res_lwr)
    lines(x = res_teo, y = res_upr)
    lines(x = res_teo, y = res_mid, lty = 2)
  }

  list(
    res_obs = res_obs,
    res_teo = res_teo,
    lower   = res_lwr,
    median  = res_mid,
    upper   = res_upr
  )
}
