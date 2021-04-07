# Compute the negative of the log-like
log.like <- function(y, mu, alpha, tau)
{
  rqnorm <- qnorm(0.1e1 - tau);
  rqalpha <- (rqnorm * alpha) ^ 2;
  beta <- -log(mu) * (0.2e1 + rqalpha - rqnorm * alpha * sqrt(rqalpha + 0.4e1)) / 0.2e1;
  t1 <- log(y);
  t3 <- log(alpha);
  t5 <- log(beta);
  t7 <- 0.1e1 / t1;
  t9 <- 0.1e1 * beta * t7;
  t10 <- sqrt(-t9);
  t13 <- log(t10 - t10 * t9);
  t14 <- alpha * alpha;
  t15 <- 0.1e1 / t14;
  ll <- -0.1612085714e1 - 0.1e1 * t1 - 0.1e1 * t3 - 0.1e1 * t5 + t13 + 0.1000000000e1 * t15 + 0.5000000000e0 * t15 * t1 / beta + 0.5000000000e0 * t15 * beta * t7;
  -sum(ll)
}

llunitBSQuantReg <- function(par, tau, linkinv, linkinv_alpha, X, Z, y)
{
  # location parameter (mu)
  p <- ncol(X);
  seqp <- seq.int(length.out = p);
  beta <- par[seqp]
  mu <- linkinv(X %*% beta)

  # shape parameter (alpha)
  if (is.null(Z))
  {
    alpha <- par[-seqp]
  }
  else
  {
    q <- ncol(Z)
    gamma <- par[seqp + q]
    alpha <- linkinv_alpha(Z %*% gamma)
  }
  suppressWarnings(log.like(y, mu, alpha, tau))
  #suppressWarnings(-sum(duBS(y, mu, alpha, tau, log=TRUE)))
}

# Format output ------------------------------------------------------------
format.perc <- function(probs, digits)
{
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")
}

FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}
#
#
# # Parametrization ---------------------------------------------------------
#
# muinv <- function(mu, alpha, tau)
# {
#   rqnorm = qnorm(0.1e1 - tau, 0.0, 1.0, TRUE, FALSE)
#   t1 = log(mu)
#   t3 = rqnorm * rqnorm
#   t5 = alpha * alpha
#   t11 = sqrt(t3 * t5 + 0.4e1)
#   return(-t1 - 0.5000000000e0 * t1 * t3 * t5 + 0.5000000000e0 * t1 * rqnorm * alpha * t11)
# }
#
#
