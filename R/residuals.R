#' @title Residuals Method for \code{unitBSQuantReg} Objects
#'
#' @description Extract Cox-Snell and Quantile residuals from a unit-Birnbaum-Saunders object.
#'
#' @author
#' Josmar Mazucheli
#' André Felipe B. Menezes
#'
#' @param object fitted model object of class \code{\link{unitBSQuantReg}}.
#' @param type character indicating type of residuals.
#' @param ... currently not used.
#'
#' @references
#' Cox, D. R. and Snell, E. J. (1968). A general definition of residuals. \emph{Journal of the Royal Statistical Society. Series B (Methodological)} \bold{2} (30), 248--275.
#' Dunn, P. K. and Smyth, G. K. (1996). Randomized Quantile Residuals. \emph{Journal of Computational and Graphical Statistics} \bold{4} (3) 236--244.
#'
#' @rdname residuals.unitBSQuantReg
#' @export
residuals.unitBSQuantReg <- function(object, type = c("cox-snell", "quantile"), ...)
{
  mu  <- object$fitted.values$mu
  alpha <- exp(object$fitted.values$alpha)
  tau <- object$tau
  y   <- object$data$y

  res <- switch(type,"cox-snell" = {-log(puBS(y, mu, alpha, tau, lower.tail = FALSE))}, "quantile" = {qnorm(puBS(y, mu, alpha, tau))})
  res
}
