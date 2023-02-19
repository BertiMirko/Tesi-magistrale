# NS per dati Poisson con StARS
# viene apportata una modifica a XMRF per poter effettuare la parallelizzazione su Windows

library(BiocParallel)

source("learn2count-master/R/datasimPois.R")
source("structure_learning-master/scripts/result_scores.R")

# carico la modifica di LPGM.network,
# dove ho sostituito mclapply con bplapply
source("LPGM.network2.R")
library(XMRF)

# per fare usare la funzione modificata in XMRF
# fonte: https://stackoverflow.com/questions/24331690/modify-package-function
environment(LPGM.network2) <- asNamespace('XMRF')
assignInNamespace("LPGM.network",LPGM.network2,ns="XMRF")

pois.stars_res <- function(n, p, SAdj, nsim, lambda, lambda_err, betaval, nCpus){
  W <- matrix(NA,nsim,11)
  colnames(W) <- c("TP","FP","FN","PPV","Se","F1","time1","time2","time3","NA1","NA2")
  set.seed(1234)
  for(z in 1:nsim){
    cat(z,"\n")
    simDat <- pois.simdata(n=n, p=p, B=SAdj, lambda=lambda, lambda.c=lambda_err)
    simDat <- t(simDat)
    pois.time = system.time(lpgm.fit <- XMRF(simDat, method="LPGM",stability="star", N=20, beta=betaval, 
                                             lmin=0.01, nlams=20, parallel=TRUE, nCpus=nCpus,th=0.001))
    stars_res <- result(SAdj,lpgm.fit$network[[lpgm.fit$opt.index]])
    W[z,] <- c(stars_res,pois.time)
  }
  return(W)
}

################################################################################
## p=10

# carichiamo i grafi
# a) sample scalefree graph
load("structure_learning-master/data/scalefreesample10.RData")
SAdj1 <- SAdj

# b) sample random graph 
load("structure_learning-master/data/randomsample10-03.RData")
SAdj2 <- SAdj

# c) hub graph
load("structure_learning-master/data/hubsample10.RData")
SAdj3 <- SAdj
rm(SAdj)


p <- 10
n <- 100
lam <- 5
lam_noise <- 0.5
nsim <- 50
nCpus <- 6
betav <- 0.5

system.time(prova <- pois.stars_res(n, p, SAdj1, nsim, lam, lam_noise, betav, nCpus))

system.time(prova2 <- pois.stars_res(n, p, SAdj2, nsim, lam, lam_noise, betav, nCpus))

system.time(prova3 <- pois.stars_res(n, p, SAdj3, nsim, lam, lam_noise, betav, nCpus))


###############################################################################
## p=100

# carichiamo i grafi
# a) sample scalefree graph
load("structure_learning-master/data/scalefreesample100.RData")
SAdj1 <- SAdj

# b) sample random graph 
load("structure_learning-master/data/randomsample100-003.RData")
SAdj2 <- SAdj

# c) hub graph
load("structure_learning-master/data/hubsample100.RData")
SAdj3 <- SAdj
rm(SAdj)


p <- 100
n <- 200
lam <- 5
lam_noise <- 0.5
nsim <- 50
nCpus <- 6
betav <- 0.05

system.time(prova <- pois.stars_res(n, p, SAdj1, nsim, lam, lam_noise, betav, nCpus))

system.time(prova2 <- pois.stars_res(n, p, SAdj2, nsim, lam, lam_noise, betav, nCpus))

system.time(prova3 <- pois.stars_res(n, p, SAdj3, nsim, lam, lam_noise, betav, nCpus))