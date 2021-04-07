library(unitBSQuantReg)
#######################################################
data(BrazilElection2014, package = "baquantreg")
ind     <- which(BrazilElection2014$UF_name == "SE")
data.se <- BrazilElection2014[ind, c("percVotes", "HDI")]

taus    <- c(0.10, 0.25, 0.50, 0.75, 0.90)
fits    <- lapply(taus, function(TAU) unitBSQuantReg(percVotes ~ HDI, data = data.se, tau = TAU))
lapply(fits,  summary)
#######################################################
bodyfat <- read.table("http://www.leg.ufpr.br/lib/exe/fetch.php/publications:papercompanions:qbmult_dataset.csv", sep = ",", header = T)

bodyfat$z1 <- bodyfat$AGE
bodyfat$z2 <- bodyfat$BMI / 100
bodyfat$z3 <- with(bodyfat, ifelse(SEX == 1, 0, 1))
bodyfat$z4 <- with(bodyfat, ifelse(IPAQ == 1, 1,0))
bodyfat$z5 <- with(bodyfat, ifelse(IPAQ == 2, 1,0))

bodyfat$arms <- bodyfat$ARMS / 100;
taus    <- c(0.10, 0.25, 0.50, 0.75, 0.90)
fits    <- lapply(taus, function(TAU) unitBSQuantReg(arms ~ z1 + z2 + z3 + z4 + z5, data = bodyfat, tau = TAU))
lapply(fits,  summary)







