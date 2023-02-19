# NEIGHBORHOOD SELECTION con Lasso - dati Gaussiani

library(mvtnorm)
library(igraph)
library(matrixcalc)
library(glmnet)
library(lars)
library(doSNOW)
library(doParallel)
library(doMPI)
library(simPATHy)
source("structure_learning-master/scripts/result_scores.R")

gauss.nr <- function(X,alpha,extend){
  
  p <- ncol(X)
  n <- nrow(X)
  adj <- matrix(1,p,p)
  diag(adj) <- 0
  X_s <- scale(X)
  
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

nod_regr_gauss2 <- function(n, SAdj, Sigma, mu, nsim, alpha){
  W <- matrix(NA,nsim,11)
  colnames(W) <- c("TP","FP","FN","PPV","Se","F1","time1","time2","time3","NA1","NA2")
  set.seed(123)
  for (i in 1:nsim){
    cat(i,"\n")
    ##simulo i dati 
    X <- rmvnorm(n,mean=mu,sigma=Sigma)
    
    ##stima 
    gauss.time = system.time(adj.gauss <- try(gauss.nr(X,alpha,extend = TRUE),silent = TRUE))
    if (class (adj.gauss)[1]!="try-error"){
      Nr.gauss = result(SAdj,adj.gauss)
    }else{
      adj.gauss = matrix(NA,p,p)
      Nr.gauss = result(SAdj,adj.gauss)
    }
    W[i,] <- c(Nr.gauss,gauss.time)
  }
  return(W)
}


################################################################################
## p=10

# a) sample scalefree graph
load("structure_learning-master/data/scalefreesample10.RData")
SAdj1 <- SAdj
g1 <- graph_from_adjacency_matrix(SAdj1,mode="undirected")
g1 <- as_graphnel(g1)

# b) sample random graph 
load("structure_learning-master/data/randomsample10-03.RData")
SAdj2 <- SAdj
g2 <- graph_from_adjacency_matrix(SAdj2,mode="undirected")
g2 <- as_graphnel(g2)

# c) hub graph
load("structure_learning-master/data/hubsample10.RData")
SAdj3 <- SAdj
rm(SAdj)
g3 <- graph_from_adjacency_matrix(SAdj3,mode="undirected")
g3 <- as_graphnel(g3)

# Generazione di Omega dalla matrice di adiacenza
# e poi di Sigma invertendo Omega

S1 <- as.matrix(SAdj1)
h <- 1
set.seed(123)
for(i in 1:nrow(S1)){
  for(j in h:ncol(S1)){
    if(S1[i,j]!=0) S1[i,j] = S1[j,i] = 0.245  
  }
  h <- h+1
}
diag(S1) <- 1
Sigma1 <- solve(S1)
colnames(Sigma1)=rownames(Sigma1)=1:10
Sigma1<-fitSgraph(g1,Sigma1)
Sigma1 <- round(Sigma1,5)
is.symmetric.matrix(Sigma1)
is.positive.definite(Sigma1)

S2 <- as.matrix(SAdj2)
h <- 1
set.seed(123)
for(i in 1:nrow(S2)){
  for(j in h:ncol(S2)){
    if(S2[i,j]!=0) S2[i,j] = S2[j,i] = 0.245 
  }
  h <- h+1
}
diag(S2) <- 1
Sigma2 <- solve(S2)
colnames(Sigma2)=rownames(Sigma2)=1:10
Sigma2<-fitSgraph(g2,Sigma2)
Sigma2 <- round(Sigma2,5)
is.symmetric.matrix(Sigma2)
is.positive.definite(Sigma2)

S3 <- as.matrix(SAdj3)
h <- 1
set.seed(123)
for(i in 1:nrow(S3)){
  for(j in h:ncol(S3)){
    if(S3[i,j]!=0) S3[i,j] = S3[j,i] = 0.245 
  }
  h <- h+1
}
diag(S3) <- 1
Sigma3 <- solve(S3)
colnames(Sigma3)=rownames(Sigma3)=1:10
Sigma3<-fitSgraph(g3,Sigma3)
Sigma3 <- round(Sigma3,5)
is.symmetric.matrix(Sigma3)
is.positive.definite(Sigma3)


# per calcolo in parallelo
cl <- makeCluster(6)  # 6 = numero di cores - 2
registerDoParallel(cl)


n <- 100
nsim <- 50
alpha <- 0.05
mu <- rep(0,10)
names(mu) <- 1:10

system.time(prova <- nod_regr_gauss2(n=n,SAdj=SAdj1,Sigma=Sigma1,mu=mu,nsim=nsim,alpha=alpha))

system.time(prova2 <- nod_regr_gauss2(n=n,SAdj=SAdj2,Sigma=Sigma2,mu=mu,nsim=nsim,alpha=alpha))

system.time(prova3 <- nod_regr_gauss2(n=n,SAdj=SAdj3,Sigma=Sigma3,mu=mu,nsim=nsim,alpha=alpha))


################################################################################
## p=100

# a) sample scalefree graph
load("structure_learning-master/data/scalefreesample100.RData")
SAdj1 <- SAdj
g1 <- graph_from_adjacency_matrix(SAdj1,mode="undirected")
g1 <- as_graphnel(g1)

# b) sample random graph 
load("structure_learning-master/data/randomsample100-003.RData")
SAdj2 <- SAdj
g2 <- graph_from_adjacency_matrix(SAdj2,mode="undirected")
g2 <- as_graphnel(g2)

# c) hub graph
load("structure_learning-master/data/hubsample100.RData")
SAdj3 <- SAdj
rm(SAdj)
g3 <- graph_from_adjacency_matrix(SAdj3,mode="undirected")
g3 <- as_graphnel(g3)


S1 <- as.matrix(SAdj1)
h <- 1
set.seed(123)
for(i in 1:nrow(S1)){
  for(j in h:ncol(S1)){
    if(S1[i,j]!=0) S1[i,j] = S1[j,i] = 0.245 
  }
  h <- h+1
}
diag(S1) <- 1
Sigma1 <- solve(S1)
colnames(Sigma1)=rownames(Sigma1)=1:100
Sigma1<-fitSgraph(g1,Sigma1)
Sigma1 <- round(Sigma1,5)
is.symmetric.matrix(Sigma1)
is.positive.definite(Sigma1)

S2 <- as.matrix(SAdj2)
h <- 1
set.seed(123)
for(i in 1:nrow(S2)){
  for(j in h:ncol(S2)){
    if(S2[i,j]!=0) S2[i,j] = S2[j,i] = 0.245 
  }
  h <- h+1
}
diag(S2) <- 1
Sigma2 <- solve(S2)
colnames(Sigma2)=rownames(Sigma2)=1:100
Sigma2<-fitSgraph(g2,Sigma2)
Sigma2 <- round(Sigma2,5)
is.symmetric.matrix(Sigma2)
is.positive.definite(Sigma2)

S3 <- as.matrix(SAdj3)
h <- 1
set.seed(123)
for(i in 1:nrow(S3)){
  for(j in h:ncol(S3)){
    if(S3[i,j]!=0) S3[i,j] = S3[j,i] = 0.22    # 0.22 per garantire che sia definita positiva
  }
  h <- h+1
}
diag(S3) <- 1
Sigma3 <- solve(S3)
colnames(Sigma3)=rownames(Sigma3)=1:100
Sigma3<-fitSgraph(g3,Sigma3)
Sigma3 <- round(Sigma3,5)
is.symmetric.matrix(Sigma3)
is.positive.definite(Sigma3)


n <- 200
nsim <- 50
alpha <- 0.05
mu <- rep(0,100)
names(mu) <- 1:100

system.time(prova <- nod_regr_gauss2(n=n,SAdj=SAdj1,Sigma=Sigma1,mu=mu,nsim=nsim,alpha=alpha))

system.time(prova2 <- nod_regr_gauss2(n=n,SAdj=SAdj2,Sigma=Sigma2,mu=mu,nsim=nsim,alpha=alpha))

system.time(prova3 <- nod_regr_gauss2(n=n,SAdj=SAdj3,Sigma=Sigma3,mu=mu,nsim=nsim,alpha=alpha))

stopCluster(cl)