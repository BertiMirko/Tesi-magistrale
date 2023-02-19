# NODEWISE REGRESSION per DATI POISSON trasformati con log(x+1)

library(igraph)
library(glmnet)
library(lars)
library(doSNOW)
library(doParallel)
library(doMPI)
source("learn2count-master/R/datasimPois.R")
source("structure_learning-master/scripts/result_scores.R")

pois.nr.t <- function(X,alpha,extend){
  
  p <- ncol(X)
  n <- nrow(X)
  adj <- matrix(1,p,p)
  diag(adj) <- 0
  X.tr <- log(X+1)
  X_s <- scale(X.tr)
  
  V <- foreach(i=1:p, .combine="cbind",.packages = c("glmnet")) %dopar%{
    # calcolo il lambda ottimale
    sigma_a_hat <- as.numeric(sqrt((n^-1)*(X_s[,i]%*%X_s[,i])))
    lambda_opt <- (2*sigma_a_hat/sqrt(n))*(pnorm(alpha/(2*p*n^2),lower.tail=F)^(-1))
    m.lasso=glmnet(x=X_s[,-i],y=X_s[,i],alpha=1,lambda=lambda_opt,intercept=FALSE,family="gaussian")
    indxs <- which(m.lasso$beta==0)  
    indxs[indxs>=i] <- indxs[indxs>=i]+1
    adj[indxs,i] <- 0
    return(adj[,i])
  }
  if (extend == TRUE){
    adj <- V + t(V)
    adj[which(adj != 0)] <-1
  }else{
    adj <- V * t(V)
  }
  return(adj)
}


nod_regr.t <- function(n, p, SAdj, nsim, alpha,mu,mu.nois){
  W <- matrix(NA,nsim,11)
  colnames(W) <- c("TP","FP","FN","PPV","Se","F1","time1","time2","time3","NA1","NA2")
  set.seed(123)
  for (i in 1:nsim){
    cat(i,"\n")
    ##simulo dati dal modello Poisson
    X <- pois.simdata(n, p,SAdj,mu,mu.nois)
    
    ##stima col modello Poisson 
    pois.time = system.time(adj.pois <- try(pois.nr.t(X,alpha,extend = TRUE),silent = TRUE))
    if (class (adj.pois)[1]!="try-error"){
      Nr.pois = result(SAdj,adj.pois)
    }else{
      adj.pois = matrix(NA,p,p)
      Nr.pois = result(SAdj,adj.pois)
    }
    W[i,] <- c(Nr.pois,pois.time)
  }
  return(W)
}


# per calcolo in parallelo
cl <- makeCluster(6)  # 6 = numero di cores - 2
registerDoParallel(cl)

## p=10

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

mu <- 0.5
mu_noise <- 0.5
nsim <- 50
n <- 100
p <- 10
alpha <- 0.05

system.time(prova <- nod_regr.t(n,p,SAdj1,nsim,alpha=alpha,mu,mu_noise))

system.time(prova2 <- nod_regr.t(n,p,SAdj2,nsim,alpha=alpha,mu,mu_noise))

system.time(prova3 <- nod_regr.t(n,p,SAdj3,nsim,alpha=alpha,mu,mu_noise))


#############################
## p=100

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


# 1) mu=0.5
mu <- 0.5
mu_noise <- 0.5
nsim <- 50
n <- 200
p <- 100
alpha <- 0.05

system.time(prova <- nod_regr.t(n,p,SAdj1,nsim,alpha=alpha,mu,mu_noise))

system.time(prova2 <- nod_regr.t(n,p,SAdj2,nsim,alpha=alpha,mu,mu_noise))

system.time(prova3 <- nod_regr.t(n,p,SAdj3,nsim,alpha=alpha,mu,mu_noise))

stopCluster(cl)