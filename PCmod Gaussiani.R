## PC modificato per dati Gaussiani

library(mvtnorm)
library(igraph)
library(matrixcalc)
library(doSNOW)
library(doParallel)
library(doMPI)
library(simPATHy)
source("structure_learning-master/scripts/result_scores.R")

gauss.wald2 <- function(X,maxcard,alpha,extend){
  
  p <- ncol(X)
  n <- nrow(X)
  adj <- matrix(1,p,p)
  diag(adj) <- 0
  ncard <- 0
  
  while(ncard <= maxcard){
    #calcolo i vicini per tutti i geni
    neighbor <- foreach(k = 1:p, .combine = "cbind") %do%{
      if(length(which(adj[,k]==1))==0){
        return(0)
      }else{
        return(which(adj[,k]==1))
      }
    }
    colnames(neighbor) <- 1:p 
    
    # ricavo il numero di vicini per ciascun gene/nodo
    numneigh <- apply(neighbor, 2, function(x) length(unique(x)))
    numneigh[which(neighbor[1,]==0)] <- 0
    genes <- which(numneigh >= (ncard+1))
    
    if(length(genes)!=0){
      V <- foreach(i = 1:p, .combine="rbind",.packages = c("glmGamPoi")) %dopar% {
        gs <- unique(neighbor[,i]) # considero solo i geni che hanno il gene i come vicino (quindi i vicini di i)
        gs <- intersect(gs,genes) # prendo solo quelli che hanno altri ncard vicini oltre a i
        if(length(gs)>1){
          util <- which(rowSums(adj[,gs]) >= 1)
        }else{
          util <- which(adj[,gs]==1)
        }
        if(length(setdiff(util,i)) >= ncard){
          if(length(setdiff(util,i))==1){
            poscomb <- c(setdiff(util,i))
            poscomb <- as.matrix(c(i,poscomb))
            
          }else { 
            poscomb <- combn(setdiff(util,i),ncard) 
            poscomb <- rbind(rep(i,ncol(poscomb)),poscomb)
          }
          
          excl <- 0  # per togliere i geni-risposta una volta che non sono piu' adiacenti al gene i
          
          for(j in 1:ncol(poscomb)){
            
            indxs <- which(apply(neighbor,2,function(x) all(poscomb[,j]%in%x)))
            indxs <- setdiff(indxs,excl)
            if(length(indxs)!=0){
              if(length(indxs)>1){
                mod <- lm(X[,indxs]~ scale(X[,poscomb[,j]]))
                for(r in 1:length(indxs)){
                    if(coef(summary(mod))[[r]][2,4] > alpha){   
                           adj[i,indxs[r]] <- 0
                           excl <- c(excl,r)
                        }
                      } 
              }else{
                mod <- lm(X[,indxs]~ scale(X[,poscomb[,j]]))
                 for(r in 1:length(indxs)){
                     if(coef(summary(mod))[2,4] > alpha){   
                            adj[i,indxs[r]] <- 0
                            excl <- c(excl,r)
                          }
                        }
              }
            }
          }
        } 
        return(adj[i,])
      }
    }
    if (extend == TRUE){
      adj <- V + t(V)
      adj[which(adj != 0)] <-1
    }else
      
      adj <- V * t(V)
    
    ncard <- ncard + 1
  } 
  
  return(adj)
}


result_PC2Gauss <- function(n, SAdj, Sigma, mu, nsim, maxcard, alpha){
  W <- matrix(NA,nsim,11)
  colnames(W) <- c("TP","FP","FN","PPV","Se","F1","time1","time2","time3","NA1","NA2")
  set.seed(123)
  for (i in 1:nsim){
    cat(i,"\n")
    ##simulo i dati 
    X <- rmvnorm(n,mean=mu,sigma=Sigma)
    
    ##stima 
    gauss.time = system.time(adj.gauss <- try(gauss.wald2(X,maxcard,alpha,extend = TRUE),silent = TRUE))
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



#######################################################################################
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
b <- 0.15
nsim <- 50
alpha <- 2*(pnorm(n^b,lower.tail=F))
m <- 8
mu <- rep(0,10)
names(mu) <- 1:10

system.time(prova <- result_PC2Gauss(n=n,SAdj=SAdj1,Sigma=Sigma1,mu=mu,nsim=nsim,maxcard=m,alpha=alpha))

system.time(prova2 <- result_PC2Gauss(n=n,SAdj=SAdj2,Sigma=Sigma2,mu=mu,nsim=nsim,maxcard=m,alpha=alpha))

system.time(prova3 <- result_PC2Gauss(n=n,SAdj=SAdj3,Sigma=Sigma3,mu=mu,nsim=nsim,maxcard=m,alpha=alpha))


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
b <- 0.225
alpha <- 2*(pnorm(n^b,lower.tail=F))
mu <- rep(0,100)
names(mu) <- 1:100
m <- 3

system.time(prova <- result_PC2Gauss(n=n,SAdj=SAdj1,Sigma=Sigma1,mu=mu,nsim=nsim,maxcard=m,alpha=alpha))

system.time(prova2 <- result_PC2Gauss(n=n,SAdj=SAdj2,Sigma=Sigma2,mu=mu,nsim=nsim,maxcard=m,alpha=alpha))

system.time(prova3 <- result_PC2Gauss(n=n,SAdj=SAdj3,Sigma=Sigma3,mu=mu,nsim=nsim,maxcard=m,alpha=alpha))

stopCluster(cl)
