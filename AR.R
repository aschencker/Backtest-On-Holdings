library(PortfolioAnalytics)
library(PerformanceAnalytics)
library(foreach)
library(parallel)
library(doParallel)
library(NlcOptim)
library(zoo)

setwd("/Users/Ari's_Folder/Desktop/Backtest/SSIF Backtest/SSIF Excel")
f.sheet <- read.csv("Factors2.csv")
max.sheet <- read.csv("Max.Sheet.csv")
min.sheet <- read.csv("Min.Sheet.csv")
p.sheet <- read.csv("Portfolio.Sheet.csv")
r.sheet <- read.csv("Returns.Sheet2.csv")
s.sheet <- read.csv("SSIF.Sector.Sheet.csv")
w.sheet <- read.csv("SSIF.Weight.Sheet.csv")
b.sheet <- f.sheet[,2] + f.sheet[,10]


row.names(max.sheet) <- (max.sheet[,1])
row.names(min.sheet) <- (min.sheet[,1])
row.names(p.sheet) <- (p.sheet[,1])
row.names(r.sheet) <- (r.sheet[,1])
row.names(s.sheet) <- (s.sheet[,1])
row.names(w.sheet) <- (w.sheet[,1])

w.sheet <- w.sheet[,-1]
p.sheet <- p.sheet[,-1]
s.sheet <- s.sheet[,-1]
max.sheet <- max.sheet[,-1]
min.sheet <- min.sheet[,-1]
r.sheet <- r.sheet[,-1]

numCores <- detectCores()
registerDoParallel(cores = numCores-1)

s.sheet_row <- function(i){
  if(i < 25){
    return(1)
  } 
  if(i < 59){
    return(2)
  } 
  if(i < 80){
    return(3)
  } 
  if(i < 105){
    return(4)
  }else{
    return(5)
  }
}

AR_QP <- function(Months,Min_Max){
  if(Min_Max == "Max"){
    Min_Max <- -1
  }
  if(Min_Max == "Min"){
    Min_Max <- 1
  }
  foreach(i = 1:145, .combine = "c") %dopar% {
    Holdings <- colnames(p.sheet[i,p.sheet[i,]==1])
    Sector <- s.sheet[s.sheet_row(i),which(p.sheet[i,]==1)]
    Location <- match(as.matrix(rownames(p.sheet[i,])),rownames(r.sheet))
    Returns <- rollapply(r.sheet[(Location-(Months*21-1)):Location,Holdings]+1,width = 21,FUN = prod,by =21,by.column = TRUE)-1
    Scalars <- c()
    Factors <- rollapply(f.sheet[(Location-(Months*21-1)):Location,-1]+1,width = 21,FUN = prod,by = 21,by.column = TRUE)-1
    W <- as.numeric(w.sheet[i,])
  
    MKT <- Factors[,1]
    SMB <- Factors[,2]
    HML <- Factors[,3]
    RMW <- Factors[,4]
    CMA <- Factors[,5]
    MOM <- Factors[,6]
    ST <- Factors[,7]
    LT <- Factors[,8]
    RF <- Factors[,9]

    COND <- which(Sector == "COND")
    CONS <- which(Sector == "CONS")
    ENRS <- which(Sector == "ENRS")
    FINL <- which(Sector == "FINL")
    HLTH <- which(Sector == "HLTH")
    INDU <- which(Sector == "INDU")
    INFT <- which(Sector == "INFT")
    MATR <- which(Sector == "MATR")
    RLST <- which(Sector == "RLST")
    TELS <- which(Sector == "TELS")
    UTIL <- which(Sector == "UTIL")
    
    x <- rep(1/length(Holdings),length(Holdings))
    lb <- as.numeric(min.sheet[i,which(p.sheet[i,]==1)])
    ub <- as.numeric(max.sheet[i,which(p.sheet[i,]==1)])
    
    Coefs <- c()
    Errors <- c()
    for(j in 1:ncol(Returns)){
      Coefs <- rbind(Coefs,lm((Returns[,j]-RF)~MKT)[["coefficients"]])
      Errors <- cbind(Errors,lm((Returns[,j]-RF)~MKT)[["residuals"]])
    }
    Alpha <- Coefs[,1] * Min_Max
    Omega <- cov(Errors) * 59
    
    confun <- function(x){
      f <- NULL
      f <- rbind(f,sum(as.numeric(x)[COND])-W[1]) 
      f <- rbind(f,sum(as.numeric(x)[CONS])-W[2])
      f <- rbind(f,sum(as.numeric(x)[ENRS])-W[3])
      f <- rbind(f,sum(as.numeric(x)[FINL])-W[4])
      f <- rbind(f,sum(as.numeric(x)[HLTH])-W[5])
      f <- rbind(f,sum(as.numeric(x)[INDU])-W[6])
      f <- rbind(f,sum(as.numeric(x)[INFT])-W[7])
      f <- rbind(f,sum(as.numeric(x)[MATR])-W[8])
      if(length(RLST)>0){
        f <- rbind(f,sum(as.numeric(x)[RLST])-W[9])
      }
      f <- rbind(f,sum(as.numeric(x)[TELS])-W[10])
      f <- rbind(f,sum(as.numeric(x)[UTIL])-W[11])
      return(list(ceq=f,c=NULL))
    }
    objfun <- function(x){
      return(((Alpha %*% x)/sqrt(t(x) %*% Omega %*% x)) * 100)
    }
    
    Weights <- solnl(x,objfun,confun, lb = lb, ub = ub)[["par"]]
    
    Portfolio_Returns <- c()
    Location2 <- match(as.matrix(rownames(p.sheet[(i+1),])),rownames(r.sheet))
    for(k in 1:nrow(as.matrix(r.sheet[(Location+1):Location2,Holdings]))){
      Portfolio_Returns[k] <- as.matrix(rowSums(r.sheet[(Location+k),Holdings]*Weights))
      Weights <- ((r.sheet[(Location+k),Holdings]+1)*Weights)/sum(((r.sheet[(Location+k),Holdings]+1)*Weights))
    }
    return(Portfolio_Returns)
  }
}

AR <- AR_QP(60,"Max")


