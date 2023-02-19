# ALGORITMO PC modificato per dati Poisson

pois.wald <- function(X,maxcard,alpha,extend){
  
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
    
    numneigh <- apply(neighbor, 2, function(x) length(unique(x)))
    numneigh[which(neighbor[1,]==0)] <- 0
    genes <- which(numneigh >= (ncard+1))
    
    if(length(genes)!=0){
      V <- foreach(i = 1:p, .combine="rbind",.packages = c("glmGamPoi")) %dopar% {
        gs <- unique(neighbor[,i]) 
        gs <- intersect(gs,genes) 
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
          
          excl <- 0 # utile a togliere i geni-risposta una volta che non sono piu' adiacenti al gene i
          
          for(j in 1:ncol(poscomb)){
            indxs <- which(apply(neighbor,2,function(x) all(poscomb[,j]%in%x)))
            indxs <- setdiff(indxs,excl)
            
            if(length(indxs)!=0){
              
              mod <- glm_gp(t(X[,indxs]),design=model.matrix(~ scale(X[,poscomb[,j]])),overdispersion=FALSE,overdispersion_shrinkage = FALSE,do_cox_reid_adjustment = FALSE)
              # calcolo dei p-values (https://github.com/const-ae/glmGamPoi/issues/12)
              pred <- predict(mod, se.fit = TRUE, newdata = diag(nrow = ncol(mod$Beta)))
              pval <- t(2*(1-pnorm(abs(pred$fit/pred$se.fit))))
              colnames(pval) <- indxs
              for(r in indxs){
                if(pval[2,colnames(pval)==r] > alpha){
                  adj[i,r] <- 0
                  excl <- c(excl,r)
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

result_algorithms7 <- function(n, SAdj, nsim, maxcard, alpha, mu, mu.nois){
  W <- matrix(NA,nsim,11)
  colnames(W) <- c("TP","FP","FN","PPV","Se","F1","time1","time2","time3","NA1","NA2")
  set.seed(123)
  for (z in 1:nsim){
    cat(z,"\n")
    ## simulo dati Poisson
    X <- pois.simdata(n, p, SAdj, mu, mu.nois)
    
    ## stima con Poisson
    pois.time = system.time(adj.pois <- try(pois.wald7(X,maxcard,alpha,extend = TRUE),silent = TRUE))
    if (class (adj.pois)[1]!="try-error"){
      Pc.pois = result(SAdj,adj.pois)
    }else{
      adj.pois = matrix(NA,p,p)
      Pc.pois = result(SAdj,adj.pois)
    }
    
    
    W[z,] <- c(Pc.pois,pois.time)
  }
  return(W)
}


library(doSNOW)
library(doParallel)
library(doMPI)
source("learn2count-master/R/datasimPois.R")
source("structure_learning-master/scripts/result_scores.R")

cl <- makeCluster(6)  # 6 = numero di cores - 2 
registerDoParallel(cl)

######################################################################
# p=10

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


n <- 100
p <- 10
mu <- 0.5
mu_noise <- 0.5
b <- 0.15
alpha_opt <- 2*(pnorm(n^b,lower.tail=F))
nsim <- 50
m <- 8

system.time(prova <- result_algorithms7(n, SAdj1, nsim, m, alpha_opt, mu, mu_noise))

system.time(prova2 <- result_algorithms7(n, SAdj2, nsim, m, alpha_opt, mu, mu_noise))

system.time(prova3 <- result_algorithms7(n, SAdj3, nsim, m, alpha_opt, mu, mu_noise))


####################################################################
# p=100

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


n <- 200
p <- 100
mu <- 0.5
mu_noise <- 0.5
b <- 0.225
b <- 0.2
alpha_opt <- 2*(pnorm(n^b,lower.tail=F))
nsim <- 50
m <- 3

system.time(prova <- result_algorithms7(n, SAdj1, nsim, m, alpha_opt, mu, mu_noise))

system.time(prova2 <- result_algorithms7(n, SAdj2, nsim, m, alpha_opt, mu, mu_noise))

system.time(prova3 <- result_algorithms7(n, SAdj3, nsim, m, alpha_opt, mu, mu_noise))

stopCluster(cl)