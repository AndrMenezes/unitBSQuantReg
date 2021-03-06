% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{methods-unitBSQuantReg}
\alias{methods-unitBSQuantReg}
\alias{print.unitBSQuantReg}
\alias{summary.unitBSQuantReg}
\alias{coef.unitBSQuantReg}
\alias{vcov.unitBSQuantReg}
\alias{logLik.unitBSQuantReg}
\alias{confint.unitBSQuantReg}
\alias{fitted.unitBSQuantReg}
\title{Methods for \code{unitBSQuantReg} Objects}
\usage{
\method{print}{unitBSQuantReg}(x, digits = 4, ...)

\method{summary}{unitBSQuantReg}(object, correlation = FALSE, ...)

\method{coef}{unitBSQuantReg}(object, ...)

\method{vcov}{unitBSQuantReg}(object, ...)

\method{logLik}{unitBSQuantReg}(object, ...)

\method{confint}{unitBSQuantReg}(object, parm, level = 0.95, ...)

\method{fitted}{unitBSQuantReg}(object, type = c("mu", "alpha", "all"), ...)
}
\arguments{
\item{digits}{minimal number of _significant_ digits}

\item{...}{additional argument(s) for methods. Currently not used.}

\item{object, x}{fitted model object of class \code{\link{unitBSQuantReg}}.}

\item{correlation}{logical; if \code{TRUE}, the correlation matrix of
the estimated parameters is returned and printed. Default is \code{FALSE}.}

\item{parm}{a specification of which parameters are to be given confidence intervals,
either a vector of numbers or a vector of names. If missing, all parameters are considered.}

\item{level}{the confidence level required.}

\item{type}{character indicating type of fitted values to return.}
}
\description{
Methods for extracting information from fitted unit-Birnbaum-Saunders quantile regression
objects of class \code{\link{unitBSQuantReg}}.
}
\author{
Josmar Mazucheli

André Felipe B. Menezes
}
