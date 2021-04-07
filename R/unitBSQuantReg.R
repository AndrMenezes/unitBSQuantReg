#' @name unitBSQuantReg.control
#' @aliases unitBSQuantReg.control
#'
#' @title Control Parameters for the unit-Birnbaum-Saunders quantile regression
#'
#' @description Parameters that control fitting of the unit-Birnbaum-Saunders regression models using \code{unitBSQuantReg}.
#'
#' @author Josmar Mazucheli
#' @author André Felipe B. Menezes
#'
#' @param method characters string specifying the method argument passed to \code{optim}.
#' @param maxit integer specifying the \code{maxit} argument (maximal number of iterations) passed to \code{optim}.
#' @param reltol relative convergence tolerance passed to \code{optim}.
#' @param ... arguments passed to \code{optim}.
#'
#' @rdname unitBSQuantReg.control
#' @export
unitBSQuantReg.control <- function(method = "BFGS", maxit = 5000, reltol = 1e-8, ...) {
  out <- list(
    method = method,
    maxit = maxit,
    reltol = reltol
  )
  out <- c(out, list(...))
  if (!is.null(out$fnscale)) warning("fnscale must not be modified")
  out$fnscale <- 1
  out
}

#' @title The unit-Birnbaum-Saunders quantile regression
#'
#' @description Fit the unit-Birnbaum-Saunders quantile regression.
#'
#' @author
#' Josmar Mazucheli
#'
#' André Felipe B. Menezes
#'
#' @param formula symbolic description of the quantile model.
#' @param tau the quantile(s) to be estimated, number between 0 and 1.
#' @param data data.frame contain the variables in the model.
#' @param link character specification of the link function in the quantile model. Default is \code{logit}.
#' @param link.alpha character specification of the link function in the shape model. Default is \code{log}.
#' @param start an optional vector with starting values for all parameters
#' @param control a list of control arguments specified via \code{unitBSQuantReg.control}.
#' @param y numeric vector of response variable.
#' @param X numeric matrix. Regressor matrix for the quantile model.
#' @param Z numeric matrix. Regressor matrix for the shape model. Default is constant shape model.
#'
#' @details The estimated regression coefficients are in logit scale (default)
#' and estimated shape parameter in log-scale.
#'
#'
#' @examples
#' set.seed(123)
#' data <- data.frame(y = ruBS(n = 100, mu = 0.5, alpha = 1.0, tau = 0.5))
#' fit <- unitBSQuantReg(y ~ 1, data = data, tau = 0.5)
#' exp(fit$coefficients[1]) / (1 + exp(fit$coefficients[1])); exp(fit$coefficients[2])
#'
#' set.seed(123)
#' x <- runif(100, -1, 1)
#' eta <- 1 + 2 * x
#' mu <- exp(eta) / (1 + exp(eta))
#' data <- data.frame(y = ruBS(n = 100, mu, alpha = 2.0, tau = 0.5), x)
#' fit <- unitBSQuantReg(y ~ x, data = data, tau = 0.5)
#'
#' set.seed(123)
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' eta <- 1 + 2 * x1 + 3 * x2
#' mu <- exp(eta) / (1 + exp(eta))
#' data <- data.frame(y = ruBS(n = 100, mu, alpha = 1.0, tau = 0.5), x1, x2)
#' fit <- unitBSQuantReg(y ~ ., data = data, tau = 0.5)
#'
#' @importFrom stats optim make.link model.frame model.matrix model.response
#' @importFrom quantreg rq.fit
#' @importFrom Formula as.Formula
#' @name unitBSQuantReg
NULL

#' @rdname unitBSQuantReg
#' @export
unitBSQuantReg <- function(formula, tau, data, link, link.alpha = NULL, start = NULL,
                           control = unitBSQuantReg.control())
{
  if (missing(link)) link <- "logit"

  # building model matrices
  ff <- Formula::as.Formula(formula)
  mf <- model.frame(ff, data = data)
  y  <- model.response(mf, "numeric")
  X  <- model.matrix(ff, data = data, rhs = 1)
  p  <- ncol(X)

  if (length(ff)[2] == 1)
  {
    Z <- NULL
    q <- 1
  }
  else
  {
    Z  <- model.matrix(ff, data = data, rhs = 2)
    q  <- ncol(Z)
  }

  ## check response variable
  if(!(min(y) > 0 & max(y) < 1)) {
    stop("invalid dependent variable, all observations must be in (0, 1)")
  }

  fit <- unitBSQuantReg.fit(
    y = y,
    X = X,
    Z = Z,
    tau = tau,
    link = link,
    link.alpha = link.alpha,
    start = start,
    control = control
  )

  mu_coefficients     <- fit$par[1:p]
  alpha_coefficients  <- fit$par[(1:q) + p]
  vcov                <- solve(fit$hessian)

  # fitted values
  linkobj.mu    <- make.link(link)
  fitted_mu     <- linkobj.mu$linkinv(X %*% mu_coefficients)
  fitted_alpha  <- alpha_coefficients
  linkobj.alpha <- NULL
  if (!is.null(Z))
  {
    linkobj.alpha <- make.link(link.alpha)
    fitted_alpha  <- linkobj.alpha$linkinv(Z %*% alpha_coefficients)
  }

  # Output
  out <- list(
    call             = match.call(),
    formula          = ff,
    control          = control,
    link             = linkobj.mu,
    link.alpha       = linkobj.alpha,
    tau              = tau,
    loglik           = -fit$value,
    vcov             = vcov,
    coefficients     = fit$par,
    mu_coefficients  = mu_coefficients,
    alpha_coefficients = alpha_coefficients,
    fitted.values    = list(mu = as.numeric(fitted_mu),
                            alpha = as.numeric(fitted_alpha)),
    nobs             = length(y),
    npar             = length(fit$par),
    df.residual      = length(y) - length(fit$par),
    data             = list(X = X, Z = Z, y = y),
    df               = data
  )
  class(out) <- "unitBSQuantReg"
  out
}

#' @rdname unitBSQuantReg
#' @export
unitBSQuantReg.fit <- function(y, X, Z = NULL, tau, link, link.alpha, start,
                               control = unitBSQuantReg.control()) {
  n <- length(y)
  p <- ncol(X)

  linkobj <- make.link(link)
  method  <- control$method
  control$method <- NULL

  # Initial guess
  if (is.null(start)) {
    # For beta
    ystar        <- linkobj$linkfun(y)
    reg_ini      <- suppressWarnings(quantreg::rq.fit(X, ystar, tau = tau))
    start        <- reg_ini$coefficients
    names(start) <- colnames(X)

    # For alpha/gamma
    if (is.null(Z))
    {
      q     <- 1
      start <- c(start, "alpha" = 0.5)
      linkobj.alpha <- NULL
    }
    else
    {
      q     <- ncol(Z)
      gamma <- rep(0, q)
      start <- c(start, gamma)

      names(start)[1:q + p] <- colnames(Z)
      linkobj.alpha   <- make.link(link.alpha)
    }
  }

  # Maximization
  opt <- optim(
    par     = start,
    fn      = llunitBSQuantReg,
    method  = method,
    hessian = TRUE,
    control = control,
    X       = X,
    Z       = Z,
    y       = y,
    tau     = tau,
    linkinv = linkobj$linkinv,
    linkinv_alpha = linkobj.alpha$linkinv
  )
  ## check if the optim converged
  if (opt$convergence > 0)
  {
    opt$converged <- FALSE
    warning("optimization failed to converge")
  } else
  {
    opt$converged <- TRUE
  }
  opt
}
