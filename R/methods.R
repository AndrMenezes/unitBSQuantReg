#' @name methods-unitBSQuantReg
#' @title Methods for \code{unitBSQuantReg} Objects
#'
#' @description Methods for extracting information from fitted unit-Birnbaum-Saunders quantile regression
#' objects of class \code{\link{unitBSQuantReg}}.
#'
#' @param object,x fitted model object of class \code{\link{unitBSQuantReg}}.
#' @param digits  minimal number of _significant_ digits
#' @param correlation logical; if \code{TRUE}, the correlation matrix of
#'   the estimated parameters is returned and printed. Default is \code{FALSE}.
#' @param parm a specification of which parameters are to be given confidence intervals,
#' either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param type character indicating type of fitted values to return.
#' @param ... additional argument(s) for methods. Currently not used.
#'
#'
#' @author
#' Josmar Mazucheli
#'
#' André Felipe B. Menezes
#'
#' @importFrom stats pnorm cov2cor coef vcov printCoefmat
NULL

#' @rdname methods-unitBSQuantReg
#' @export

print.unitBSQuantReg <- function(x, digits = 4, ...)
{
  cat("\n unit-Birnbaum-Saunders quantile regression model", sep = "")

  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Mu coefficients (quantile model with ", x$link$name, " link and tau = ", x$tau,  "): \n", sep = "")

  print.default(FF(x$mu_coefficients, Digits = digits), print.gap = 2, quote = FALSE)

  cat("\n")

  if (!is.null(x$formula.alpha))
  {
    cat("Alpha coefficients (shape model with ", x$link.alpha$name, " link):", "\n", sep = "")

    print.default(FF(x$alpha_coefficients, Digits = digits), print.gap = 2, quote = FALSE)

    cat("\n")
  }
  else
  {
    cat("Model with constant shape parameter:", "\n", sep = "")

    print.default(format(x$alpha_coefficients, digits = digits), print.gap = 2, quote = FALSE)

    cat("\n")
  }
  invisible(x)
}

# Summary -----------------------------------------------------------------
#' @rdname methods-unitBSQuantReg
#' @export
summary.unitBSQuantReg <- function(object, correlation = FALSE, ...)
{
  estimates <- object$coefficients
  stderror  <- sqrt(diag(object$vcov))
  zvalue    <- estimates/stderror
  pvalue    <- 2 * pnorm(-abs(zvalue))
  table     <- cbind("Estimate"    = estimates,
                     "Std. Error"  = stderror,
                     "Z value"     = zvalue,
                     "Pr(>|z|)"    = pvalue)
  if (correlation)
  {
    correlation <- cov2cor(object$vcov)
  }

  out <- list(coeftable   = table,
              loglik      = object$loglik,
              correlation = correlation,
              call        = object$call,
              tau         = object$tau,
              link.mu     = object$link$name,
              link.alpha  = object$link.alpha$name,
              data        = object$data)
  class(out) <- "summary.unitBSQuantReg"
  out
}

# Print output summary ----------------------------------------------------
#' @export
print.summary.unitBSQuantReg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  p <- ncol(x$data$X)
  q <- ifelse(is.null(x$data$Z), 1, ncol(x$data$Z))

  cat("\n Wald-tests for the unit-Birnbaum-Saunders quantile regression model", "\n" ,sep = "")

  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Mu coefficients: (quantile model with ", x$link.mu, " link and tau = ", x$tau,  "): \n", sep = "")
  printCoefmat(x$coeftable[1:p, , drop = FALSE], digits = digits, has.Pvalue = TRUE)
  cat("\n")

  if(q > 1)
  {
    cat("alpha coefficients (shape model with ", x$link.alpha, " link):", "\n", sep = "")
    printCoefmat(x$coeftable[1:p + q, , drop = FALSE], digits = digits, has.Pvalue = TRUE)
    cat("\n")
  }

  else
  {
    id <- 1:p + q
    cat("Model with constant shape:", "\n", sep = "")
    printCoefmat(x$coeftable[id[length(id)], , drop = FALSE], digits = digits, has.Pvalue = TRUE)
    cat("\n")
  }

  if (is.matrix(x$correlation))
  {
    cat("Correlation of coefficients:", "\n", sep = "")
    corr <- x$correlation
    corr <- format(round(corr, 2L), nsmall = 2L, digits = digits)
    corr[!lower.tri(corr)] <- ""
    print(corr[-1, -ncol(corr), drop = FALSE], quote = FALSE)
    cat("\n")
  }
  invisible(x)
}

# coef function -----------------------------------------------------------
#' @rdname methods-unitBSQuantReg
#' @export
coef.unitBSQuantReg <- function(object, ...)
{
  if (!missing(...))
  {
    warning("Extra arguments discarded")
  }
  object$coefficients
}

# vcov function -----------------------------------------------------------
#' @rdname methods-unitBSQuantReg
#' @export
vcov.unitBSQuantReg <- function(object, ...)
{
  if (!missing(...))
  {
    warning("Extra arguments discarded")
  }
  object$vcov
}


# logLik function ---------------------------------------------------------
#' @rdname methods-unitBSQuantReg
#' @export
logLik.unitBSQuantReg <- function(object, ...)
{
  if (!missing(...))
  {
    warning("extra arguments discarded")
  }
  ll <- object$loglik
  attr(ll, "df")   <- object$npar
  attr(ll, "nobs") <- object$nobs
  class(ll) <- "logLik"
  ll
}

# confint function --------------------------------------------------------
#' @rdname methods-unitBSQuantReg
#' @export
confint.unitBSQuantReg <- function(object, parm, level = 0.95, ...)
{
  cf <- coef(object)
  ses <- sqrt(diag(vcov(object)))
  pnames <- names(ses)
  if (missing(parm)) {
    parm <- pnames
  }
  else if (is.numeric(parm)) {
    parm <- pnames[parm]
  }
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  pct <- format.perc(a, 3)
  ci <- array(NA_real_, dim = c(length(parm), 2L),
              dimnames = list(parm, pct))
  ci[] <- cf[parm] + ses[parm] %o% fac
  ci
}

#' @rdname methods-unitBSQuantReg
#' @export
fitted.unitBSQuantReg <- function(object, type = c("mu", "alpha", "all"),  ...)
{
  if (!missing(...)) {
    warning("Extra arguments discarded")
  }

  if (length(type) != 1) type <- "mu"

  switch (type,
          "all" = object$fitted.values,
          "mu"  = object$fitted.values$mu,
          "alpha" = object$fitted.values$alpha
  )
}


