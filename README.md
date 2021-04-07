
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `unitBSQuantReg`: unit-Birnbaum-Saunders quantile regression

<!-- badges: start -->

<!-- badges: end -->

The goal of unitBSQuantReg is to provide tools for fitting
unit-Birnbaum-Saunders quantile regression proposed by Mazucheli et
al. (2021)

## Installation

You can install the development version of `unitBSQuantReg` from
[GitHub](https://github.com/AndrMenezes/unitBSQuantReg) with:

``` r
# install.packages("remotes")
remotes::install_github("AndrMenezes/unitBSQuantReg")
```

## Example

The packages follows the structure of `glm` object. The main function is
`unitBSQuantReg`.

``` r
library(unitBSQuantReg)
u <- "http://www.leg.ufpr.br/lib/exe/fetch.php/publications:papercompanions:qbmult_dataset.csv"
bodyfat <- read.table(u, sep = ",", header = TRUE)

bodyfat$BMI <- bodyfat$BMI / 100
bodyfat$SEX <- factor(bodyfat$SEX)
bodyfat$IPAQ <- factor(bodyfat$IPAQ)
bodyfat$ARMS <- bodyfat$ARMS / 100
taus <- c(0.10, 0.25, 0.50, 0.75, 0.90)
fits <- lapply(taus, function(t) {
 unitBSQuantReg(ARMS ~ BMI + AGE + SEX + IPAQ, data = bodyfat, tau = t) 
})
names(fits) <- paste0("tau = ", taus)

sapply(fits, coef)
#>                tau = 0.1   tau = 0.25    tau = 0.5   tau = 0.75    tau = 0.9
#> (Intercept) -3.393706436 -3.105268427 -2.800573308 -2.522304596 -2.288775811
#> BMI          9.525014787  9.045273907  8.581546439  8.161440329  7.820028174
#> AGE          0.004775823  0.004597843  0.004253661  0.004053604  0.003905743
#> SEX2        -1.000554224 -0.946797349 -0.895571787 -0.850003490 -0.813349771
#> IPAQ1       -0.127040271 -0.122241236 -0.116000193 -0.111150056 -0.108209997
#> IPAQ2       -0.271372687 -0.256254591 -0.244194061 -0.232195603 -0.222719661
#> alpha        0.159770328  0.159726267  0.159501076  0.159302666  0.159330034
```

The currently methods implemented are

``` r
methods(class = "unitBSQuantReg")
#> [1] coef      confint   fitted    hnp       logLik    print     residuals
#> [8] summary   vcov     
#> see '?methods' for accessing help and source code
```

  - Mazucheli, J., Leiva, V., Alves, B., Menezes, A. F. B. (2021) A New
    Quantile Regression for Modeling Bounded Data Under a Unit
    Birnbaum-Saunders Distribution with Applications in Medicine and
    Politics. [*Symmetry*](https://www.mdpi.com/journal/symmetry), to
    appear.
