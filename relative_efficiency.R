#############################################################################
# Joint Modeling and Multiple Comparisons with the Best of Data from a SMART
# with Survival Outcomes
# Data: Simulated data from function simSMART()
# Output: 
#   - Table 2 in the manuscript
# Last modified: January 29, 2020
#############################################################################
rm(list = ls())
library(tidyverse)
library(survival)
library(DTR)
library(MASS)
source('helper.R') # Misc. helper functions
source('simSMART.R') # simulate SMART with survival outcomes
source('calH.R') # Estimate 2 non-parametric baseline hazards for response and death
source('loglik.R') # Estimate log-likelihood

########################################################## Set parameters ####
bootstrapN = 5 # number of bootstrap replication; set to 1000 in paper
seed=1234; n=500; tolerance=1e-5; maxit=50; pi.x <- 0.5; pi.z1.R <- .5; pi.z1.NR <- .5; pi.z2.R=.5; pi.z2.NR <- .5
Tcheck=5; adt=c(3,20); decimal=1; bshape1=3/4; bscale1=4; bshape2=3; bscale2=4
tbetas <- c(0.6, -1.3, -0.7, -1, 1.8)#eta;
thetaZ=c('X'); gammaZ=c('X'); etaZ=c('X','Zr','X1Zr1'); muZ=NULL
np.t=length(thetaZ); np.g=length(gammaZ); np.e=length(etaZ); np.m=length(muZ); n.p=np.t+np.g+np.e+np.m 
TimeRDu = round(seq(0, adt[2], by=10^(-decimal)),decimal)
n_iter = 1000


variance_table <- matrix(NA,nrow=length(timeseq)*n_iter,ncol=13)
for (iter in 1:n_iter){
  savefit = list(NULL)
  simudata_raw <- simSMART(n, tbetas, pi.x, pi.z1.R, pi.z1.NR, pi.z2.R, pi.z2.NR, 
                           Tcheck, bshape1, bscale1, bshape2, bscale2, adt, decimal)
  
  ## Variance of Guo & Tsiatis (2005) Weighted-risk set estimator:
  simudata_raw$TR = 0; simudata_raw$TR[simudata_raw$R==1]=simudata_raw$T1[simudata_raw$R==1]
  simudata_raw$Z=simudata_raw$Zr; simudata_raw$Z[is.na(simudata_raw$Z)]=0
  simudata_raw$U=simudata_raw$T2
  simudata_raw$delta=simudata_raw$C2
  simudata2=subset(simudata_raw, select=c(X,Z,TR,R,U,delta))
  est01 = WRSEestimate(data=simudata2); est01
  guo_tsiatis = cbind(est01$time, est01$SE11, est01$SE12, est01$SE21, est01$SE22)
  colnames(guo_tsiatis) = c('time', 'SE11', 'SE12', 'SE21', 'SE22')
  ## Variance of Lunceford et al (2002)
  est01LDT = LDTestimate(simudata2)
  lunceford = cbind(est01LDT$time, est01LDT$SE11, est01LDT$SE12, est01LDT$SE21, est01LDT$SE22)
  colnames(lunceford) = c('time', 'SE11', 'SE12', 'SE21', 'SE22')
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
    
    fit.A1B1=sapply(TimeRDu, FUN=s.km2.pred, theta=exp(Zt.V), gamma=exp(Zg.V), eta=exp(Ze.1), mu=exp(Zm.1))
    fit.A1B2=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V), gamma=exp(Zg.V), eta=exp(Ze.1+Ze.Zr), mu=exp(Zm.1))
    fit.A2B1=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V+Zt.X), gamma=exp(Zg.V+Zg.X), eta=exp(Ze.1+Ze.X+Ze.X1Zr0), mu=exp(Zm.1+Zm.X))
    fit.A2B2=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V+Zt.X), gamma=exp(Zg.V+Zg.X), eta=exp(Ze.1+Ze.X+Ze.Zr+Ze.X1Zr1+Ze.X1Zr1V1), mu=exp(Zm.1+Zm.X))
    savefit[[bootiter]] = cbind(TimeRDu, fit.A1B1, fit.A1B2, fit.A2B1, fit.A2B2)
  }
  
  TimeRDu =  savefit[[1]][,1]; len.et=length(TimeRDu)
  seqt = which(TimeRDu %in% timeseq)
  sigma = matrix(NA,nrow=length(seqt), ncol=4)
  for (t in 1:length(seqt)){
    St=do.call("rbind",lapply(savefit, function(x) x[seqt[t],]))
    StV0 = St[,2:5]
    sigma[t,] = apply(StV0,2,sd)
  }
  joint.sigma=cbind(timeseq,sigma)
  colnames(joint.sigma)[1] <- "time"
  
  variance_table[{length(timeseq) * (iter - 1) + 1}:{length(timeseq) * iter},] <- joint.sigma %>% 
    merge(lunceford,by = "time",all.x=T) %>% 
    merge(guo_tsiatis,by = "time",all.x=T) %>% as.matrix()
}

result <- variance_table %>% as.tibble
result_table <- result %>% 
  group_by(V1) %>% 
  summarize(joint_A = mean(V2,na.rm=T),
            joint_B = mean(V3,na.rm=T),
            joint_C = mean(V4,na.rm=T),
            joint_D = mean(V5,na.rm=T),
            lunce_A = mean(V6,na.rm=T),
            lunce_B = mean(V7,na.rm=T),
            lunce_C = mean(V8,na.rm=T),
            lunce_D = mean(V9,na.rm=T),
            guo_A = mean(V10,na.rm=T),
            guo_B = mean(V11,na.rm=T),
            guo_C = mean(V12,na.rm=T),
            guo_D = mean(V13,na.rm=T))



