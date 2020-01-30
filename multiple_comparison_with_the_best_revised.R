########################################################## 
# Multiple comparison with the best 
# Modified from : Ertefaie et al (2015)
# Last modified: April 13, 2017
########################################################## 
rm(list = ls())
library(tidyverse)
library(survival)
library(DTR)
library(MASS)
library(parallel)
RNGkind("L'Ecuyer-CMRG")
setwd("/home/ycchao/jobs/Qui")

source('helper.R') # Misc. helper functions
source('simSMART.R') # simulate SMART with survival outcomes
source('calH.R') # Estimate 2 non-parametric baseline hazards for response and death
source('loglik.R') # Estimate log-likelihood

########################################################## Set parameters ####
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID") %>% as.numeric()

args <- commandArgs(trailingOnly = T)

seed <- as.numeric(args[1])
n <- as.numeric(args[2])   # sample size
Tcheck <- as.numeric(args[3])
adt <- eval(parse(text=args[4]))
delta <- as.numeric(args[5])
n_iter <- as.numeric(args[6])
bootstrapN <- as.numeric(args[7])
tseq <- eval(parse(text=args[8])) # time at which we want to compare S(t) between regimens
# delta=1.3 #results in paper show delta ranging between 0 and 2

tolerance=1e-5; maxit=50; pi.x <- 0.5; pi.z1.R <- .4;
pi.z1.NR <- .6; pi.z2.R=.6; pi.z2.NR <- .7; Tcheck=5; adt=c(3,20); 
decimal=1; bshape1=3/4; bscale1=4; bshape2=3; bscale2=4
thetaZ=c('X'); gammaZ=c('X'); etaZ=c('X','Zr'); muZ=NULL
np.t=length(thetaZ); np.g=length(gammaZ); np.e=length(etaZ); np.m=length(muZ); n.p=np.t+np.g+np.e+np.m 
TimeRDu = round(seq(0, adt[2], by=10^(-decimal)),decimal)

tbetas=c(delta, -delta, -delta/4, -delta/2) 

set.seed(seed)
s <- .Random.seed
i <- 1
while(i < task_id){
  s <- nextRNGStream(s)
  i <- i + 1
}
.Random.seed <- s

result_table <- matrix(NA,nrow=length(tseq)*n_iter,ncol=5)
mean_table <- matrix(NA,nrow=length(tseq)*n_iter,ncol=5)
# V_all <- matrix(NA,nrow=4*n_iter,ncol=4*length(tseq))
for(iter in 1:n_iter){
  
  savefit = list(NULL)
  simudata_raw <- simSMART(n, tbetas, pi.x, pi.z1.R, pi.z1.NR, pi.z2.R, pi.z2.NR, 
                       Tcheck, bshape1, bscale1, bshape2, bscale2, adt, decimal)
  
  ########################################################## Get estimate of S(t) ####
  Tind <- order(simudata_raw$T1,simudata_raw$T2)
  sorted <- simudata_raw[Tind,]
  len.et=length(TimeRDu)
  T1match= match(sorted$T1,TimeRDu)
  Tnr = which(TimeRDu==Tcheck); if(length(Tnr)==0){Tnr=len.et}
  
  thetaZs <- as.matrix(subset(sorted,select=thetaZ)) # covariate set for theta
  gammaZs <- as.matrix(subset(sorted,select=gammaZ)) # covariate set for gamma
  etaZs <- as.matrix(subset(sorted,select=etaZ)) # covariate set for gamma
  muZs <- as.matrix(subset(sorted,select=muZ)) # covariate set for gamma
  
  ## Change to counting process
  dN1r=dN1d=dN1r2d=dN1nr2d=array(0,len.et)
  for (t in 1:len.et){
    dN1r[t]=sum((sorted$T1==TimeRDu[t] & sorted$R==1))
    dN1d[t]=sum((sorted$T1==TimeRDu[t] & sorted$D==1))  
    dN1r2d[t]=sum((sorted$T2==TimeRDu[t] & sorted$R==1 & sorted$C2==1))  
    dN1nr2d[t]=sum((sorted$T2==TimeRDu[t] & sorted$NR==1 & sorted$C2==1))  
  }
  dNdeath=dN1d+dN1r2d+dN1nr2d
  
  risk.1 <- risk.1r2d <- risk.1nr2d <-array(0,dim=c(n,len.et))
  event.1r <- event.1d <- event.1r2d<- event.1nr2d <-array(0,dim=c(n,len.et))
  for (i in 1:nrow(sorted)){
    subjecti = sorted[i,]
    risk.1[i,] <- (subjecti$T1>=TimeRDu)
    risk.1r2d[i,] <- (subjecti$T1<TimeRDu)*(subjecti$R==1)*(subjecti$T2>=TimeRDu) 
    risk.1nr2d[i,] <- (subjecti$T1<TimeRDu)*(subjecti$NR==1)*(subjecti$T2>=TimeRDu)
    event.1r[i,] <- (subjecti$T1==TimeRDu)*(subjecti$R==1)
    event.1d[i,] <- (subjecti$T1==TimeRDu)*(subjecti$D==1)
    event.1r2d[i,] <- (subjecti$T2==TimeRDu)*(subjecti$R==1)*(subjecti$C2==1)
    event.1nr2d[i,] <- (subjecti$T2==TimeRDu)*(subjecti$NR==1)*(subjecti$C2==1)
  }   
  
  ## RUN MLE (may take a couple minutes for n=500)####
  test <- optim(par=rep(0.5,n.p),fn=loglik,gr=NULL, method = "L-BFGS-B", 
                lower=rep(-10,n.p), upper=rep(10,n.p), hessian=T, control=list(maxit=1000,trace=1))
  para.b <- test$p
  hessian.mat=test$hessian
  p.se <- sqrt(abs(diag(solve(-hessian.mat))))
  
  ## SURVIVAL PREDICTIONS ####
  H.fit=calH(para.b)
  H.fit.1=H.fit$H1 
  H.fit.2=H.fit$H2 
  Zt.V=0; Zt.X=para.b[1]
  Zg.V=0; Zg.X=para.b[2]
  Ze.1=0; Ze.X=para.b[3]; Ze.Zr=para.b[4]; Ze.X1Zr1=0; Ze.X1Zr0=0; Ze.X1Zr1V1=0
  Zm.1=Zm.X=Zm.Znr=Zm.XZnr=0
  
  fit.A1B1=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V), gamma=exp(Zg.V), eta=exp(Ze.1), mu=exp(Zm.1))
  fit.A1B2=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V), gamma=exp(Zg.V), eta=exp(Ze.1+Ze.Zr), mu=exp(Zm.1))
  fit.A2B1=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V+Zt.X), gamma=exp(Zg.V+Zg.X), eta=exp(Ze.1+Ze.X+Ze.X1Zr0), mu=exp(Zm.1+Zm.X))
  fit.A2B2=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V+Zt.X), gamma=exp(Zg.V+Zg.X), eta=exp(Ze.1+Ze.X+Ze.Zr+Ze.X1Zr1+Ze.X1Zr1V1), mu=exp(Zm.1+Zm.X))
  savemean = cbind(TimeRDu, fit.A1B1, fit.A1B2, fit.A2B1, fit.A2B2)
  
  
  ########################################################## Get bootstrapped SE estimate of S(t) ####
  for (bootiter in 1:bootstrapN){
    if( bootiter %% 5 == 0 ) cat(paste('boostrap',bootiter),'\n')
    sampled_id <- sample(1:n,n,replace=T)
    simudata <- simudata_raw[sampled_id,]   # bootstrapped data set
    Tind <- order(simudata$T1,simudata$T2)
    sorted <- simudata[Tind,]
    len.et=length(TimeRDu)
    T1match= match(sorted$T1,TimeRDu)
    Tnr = which(TimeRDu==Tcheck); if(length(Tnr)==0){Tnr=len.et}
    
    thetaZs <- as.matrix(subset(sorted,select=thetaZ)) # covariate set for theta
    gammaZs <- as.matrix(subset(sorted,select=gammaZ)) # covariate set for gamma
    etaZs <- as.matrix(subset(sorted,select=etaZ)) # covariate set for gamma
    muZs <- as.matrix(subset(sorted,select=muZ)) # covariate set for gamma
    
    ## Change to counting process
    dN1r=dN1d=dN1r2d=dN1nr2d=array(0,len.et)
    for (t in 1:len.et){
      dN1r[t]=sum((sorted$T1==TimeRDu[t] & sorted$R==1))
      dN1d[t]=sum((sorted$T1==TimeRDu[t] & sorted$D==1))  
      dN1r2d[t]=sum((sorted$T2==TimeRDu[t] & sorted$R==1 & sorted$C2==1))  
      dN1nr2d[t]=sum((sorted$T2==TimeRDu[t] & sorted$NR==1 & sorted$C2==1))  
    }
    dNdeath=dN1d+dN1r2d+dN1nr2d
    
    risk.1 <- risk.1r2d <- risk.1nr2d <-array(0,dim=c(n,len.et))
    event.1r <- event.1d <- event.1r2d<- event.1nr2d <-array(0,dim=c(n,len.et))
    for (i in 1:nrow(sorted)){
      subjecti = sorted[i,]
      risk.1[i,] <- (subjecti$T1>=TimeRDu)
      risk.1r2d[i,] <- (subjecti$T1<TimeRDu)*(subjecti$R==1)*(subjecti$T2>=TimeRDu) 
      risk.1nr2d[i,] <- (subjecti$T1<TimeRDu)*(subjecti$NR==1)*(subjecti$T2>=TimeRDu)
      event.1r[i,] <- (subjecti$T1==TimeRDu)*(subjecti$R==1)
      event.1d[i,] <- (subjecti$T1==TimeRDu)*(subjecti$D==1)
      event.1r2d[i,] <- (subjecti$T2==TimeRDu)*(subjecti$R==1)*(subjecti$C2==1)
      event.1nr2d[i,] <- (subjecti$T2==TimeRDu)*(subjecti$NR==1)*(subjecti$C2==1)
    }   
    
    ## RUN MLE (may take a couple minutes for n=500)####
    test <- optim(par=rep(0.5,n.p),fn=loglik,gr=NULL, method = "L-BFGS-B", 
                  lower=rep(-10,n.p), upper=rep(10,n.p), hessian=T, control=list(maxit=1000,trace=1))
    para.b <- test$p
    hessian.mat=test$hessian
    p.se <- sqrt(abs(diag(solve(-hessian.mat))))
    
    ## SURVIVAL PREDICTIONS ####
    H.fit=calH(para.b)
    H.fit.1=H.fit$H1 
    H.fit.2=H.fit$H2 
    Zt.V=0; Zt.X=para.b[1]
    Zg.V=0; Zg.X=para.b[2]
    Ze.1=0; Ze.X=para.b[3]; Ze.Zr=para.b[4]; Ze.X1Zr1=0; Ze.X1Zr0=0; Ze.X1Zr1V1=0
    Zm.1=Zm.X=Zm.Znr=Zm.XZnr=0
    
    fit.A1B1=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V), gamma=exp(Zg.V), eta=exp(Ze.1), mu=exp(Zm.1))
    fit.A1B2=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V), gamma=exp(Zg.V), eta=exp(Ze.1+Ze.Zr), mu=exp(Zm.1))
    fit.A2B1=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V+Zt.X), gamma=exp(Zg.V+Zg.X), eta=exp(Ze.1+Ze.X+Ze.X1Zr0), mu=exp(Zm.1+Zm.X))
    fit.A2B2=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V+Zt.X), gamma=exp(Zg.V+Zg.X), eta=exp(Ze.1+Ze.X+Ze.Zr+Ze.X1Zr1+Ze.X1Zr1V1), mu=exp(Zm.1+Zm.X))
    savefit[[bootiter]] = cbind(TimeRDu, fit.A1B1, fit.A1B2, fit.A2B1, fit.A2B2)
  }
  
  ########################################################## Multiple comparison with the best ####
  getC<-function(V,alpha, nrep){
    x<-mvrnorm(nrep,rep(0,4),V)
    for(i in 1:4){
      temp<-matrix(0,nrep,4)
      for(j in 1:4){
        if(j!=i){
          temp[,j]<-(x[,j]-x[,i])/sqrt(V[i,i]+V[j,j]-2*V[i,j])
        }
      }	
      temp<-temp[,-i]
      temp<-apply(temp,1,sort)
      temp<-apply(temp,2,max)
      c[i]<-quantile(temp,1-alpha)
    }
    c
  } 
  
  seqt = which(TimeRDu %in% tseq)
  for (t in 1:length(tseq)){
    St=do.call("rbind",lapply(savefit, function(x) x[seqt[t],]))[,-1]
    V=cov(St)
    # theta = St[1,]
    theta = savemean[seqt[t],-1]
    c<-rep(0,4) 
    for(i in 1:50){c<-c+getC(V,alpha=0.05, nrep=500)} 
    c<-c/50
    S.pos<-rep(1,4)
    for(i in 1:4){
      for(j in 1:4){
        S.pos[i]= S.pos[i]*(theta[i]>=theta[j]-c[i]*sqrt(V[i,i]+V[j,j]-2*V[i,j]))
      }
    }
    result_table[length(tseq) * (iter - 1) + t,] <- c(tseq[t],S.pos)
    mean_table[length(tseq) * (iter - 1) + t,] <- savemean[seqt[t],]
    # V_all[ ((iter - 1) * 4 + 1):(iter * 4), ((t - 1) * 4 + 1):(t * 4) ] <- V
    if(t==1) cat(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~','\n')
    cat(' At t = ',tseq[t],', best regimen(s) are', c('A1B1','A1B2','A2B1','A2B2')[S.pos==1],'\n')
    if(t==length(tseq)) cat(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~','\n')
  }
  
}

save_path <- Sys.getenv("RESULT")
"%+%" <- function(x,y) paste(x,y,sep="")
output_file_tag <- args[9]
output_file_name <- output_file_tag %+% sprintf("_%d_%d",n,task_id) %+% ".txt"
output_mean_name <- output_file_tag %+% sprintf("_mean_%d_%d",n,task_id) %+% ".txt"
write.table(result_table,file=file.path(save_path,output_file_name),sep=" ",col.names = F,row.names = F)
write.table(mean_table,file=file.path(save_path,output_mean_name),sep=" ",col.names = F,row.names = F)


# n <- 1000
# B <- 1000
# x <- rnorm(n*B,0,10) %>% matrix(nrow=B,ncol=n)
# sd <- apply(x,1,mean) %>% sd()
# 
# y <- matrix(NA,nrow=B,ncol=n)
# x1 <- rnorm(n,0,10)
# for(i in 1:B){
#   id <- sample(1:n,n,replace=T)
#   y[i,] <- x1[id]
# }
# sd_y <- apply(y,1,mean) %>% sd()
# 
# sd
# sd_y



# 
# beta <- c()
# for(i in 1:1000){
#   dt <- tibble(x1 = rnorm(1000,10,1),
#                    y1 = 2 * x1 + rnorm(1000,0,2))
#   model <- lm(y1~x1,dt)
#   beta <- c(beta,model$coefficients[2])
# }
# sd(beta)
# mean(beta)
# 
# beta2 <- c()
# dt <- tibble(x1 = rnorm(1000,10,1),
#              y1 = 2 * x1 + rnorm(1000,0,2))
# for(i in 1:1000){
#   sample_id <- sample(1:1000,1000,replace=T)
#   dt2 <- dt[sample_id,]
#   model2 <- lm(y1~x1,dt2)
#   beta2 <- c(beta2,model2$coefficients[2])
# }
# sd(beta2)
# mean(beta2)

