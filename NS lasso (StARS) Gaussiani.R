# NEIGHBORHOOD SELECTION per dati Gaussiani con StARS

library(mvtnorm)
library(igraph)
library(matrixcalc)
library(doSNOW)
library(doParallel)
library(doMPI)
library(simPATHy)
library(huge)
source("structure_learning-master/scripts/result_scores.R")

NSgauss_stars <- function(X,extend="and"){
  
  if(extend!="and" & extend!="or"){
    cat("Il valore inserito per extend non è valido","\n")
    break
  }
  n <- nrow(X)
  p <- ncol(X)
  
  # griglia di valori per il parametro di regolazione (rho)
  products <- rep(NA, p)
  rho_max <- 0
  for(j in 1:p){
    products <- apply(X, 2, function(x) abs(x%*%X[,j]))
    products[j] <- 0           
    rho_max <- max(rho_max, max(products))
  }
  
  rho_min <- 10^-4
  rho_grid <- seq(from=log(rho_max), to=log(rho_min), length.out=100) 
  rho_grid <- exp(rho_grid)
  
  huge_gauss <- huge(X, lambda=rho_grid, method="mb", sym=extend,verbose=FALSE)
  hugesel_gauss <- huge.select(huge_gauss,criterion="stars",verbose=FALSE)
  
  return(hugesel_gauss$refit)
}

result_NSgauss <- function(n, p, SAdj, nsim, mu, Sigma, extend="and"){
  
  set.seed(123)
  
  W <- foreach(k=1:nsim, .combine="rbind",.packages=c("mvtnorm","huge"))%dopar%{
    
    result = function(trueG ,estimatedG){
      TP <-  sum(trueG *estimatedG)
      FP <-  sum((estimatedG-trueG)==1)
      FN <- sum((estimatedG-trueG)==-1)
      PPV  <- TP/(TP+FP)
      Se <- TP/(TP+FN)
      F1 <- 2*(PPV*Se)/(PPV+Se)
      return (c( TP,FP,FN,PPV,Se,F1))
    }
    
    NSgauss_stars <- function(X,extend="and"){
      
      if(extend!="and" & extend!="or"){
        cat("Il valore inserito per extend non è valido","\n")
        break
      }
      n <- nrow(X)
      p <- ncol(X)
      
      products <- rep(NA, p)
      rho_max <- 0
      for(j in 1:p){
        products <- apply(X, 2, function(x) abs(x%*%X[,j]))
        products[j] <- 0          
        rho_max <- max(rho_max, max(products))
      }
      
      rho_min <- 10^-4
      rho_grid <- seq(from=log(rho_max), to=log(rho_min), length.out=100)   
      rho_grid <- exp(rho_grid)
      
      huge_gauss <- huge(X, lambda=rho_grid, method="mb", sym=extend,verbose=FALSE)
      hugesel_gauss <- huge.select(huge_gauss,criterion="stars",verbose=FALSE)
      
      return(hugesel_gauss$refit)
    }
    
    ## simulazione dati Gaussiani
    X <- rmvnorm(n, mean=mu, sigma=Sigma)
    
    ## stima
    gauss.time = system.time(adj.gauss <- try(NSgauss_stars(X,extend=extend),silent = TRUE))
    if (class (adj.gauss)[1]!="try-error"){
      Pc.gauss = result(SAdj,adj.gauss)
    }else{
      adj.gauss = matrix(NA,p,p)
      Pc.gauss = result(SAdj,adj.gauss)
    }
    return(c(Pc.gauss, gauss.time)) 
  }
  
  colnames(W) <- c("TP","FP","FN","PPV","Se","F1","time1","time2","time3","NA1","NA2")
  return(W)
}

#################################################################################
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

n <- 200
p <- 10
mu <- rep(0, p)
nsim=50

system.time(prova <- result_NSgauss(n, p, SAdj1, nsim, mu, Sigma1, extend="and"))

system.time(prova2 <- result_NSgauss(n, p, SAdj2, nsim, mu, Sigma2, extend="and"))

system.time(prova3 <- result_NSgauss(n, p, SAdj3, nsim, mu, Sigma3, extend="and"))


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
p <- 100
mu <- rep(0, p)
nsim=50

system.time(prova <- result_NSgauss(n, p, SAdj1, nsim, mu, Sigma1, extend="and"))

system.time(prova2 <- result_NSgauss(n, p, SAdj2, nsim, mu, Sigma2, extend="and"))

system.time(prova3 <- result_NSgauss(n, p, SAdj3, nsim, mu, Sigma3, extend="and"))

stopCluster(cl)