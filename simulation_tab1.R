#############################################################################
# Joint Modeling and Multiple Comparisons with the Best of Data from a SMART
# with Survival Outcomes
# Data: Simulated data from function simSMART() or simSMART_two_censoring()
# Estimators: 
#   - Beta coefficients from NPMLE using optim()
#   - SE from hessian matrix
# Output: 
#   - Table 1 in the manuscript (independent censoring)
#   - Table S3 in the manuscript (dependent censoring)
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
seed=1234; n=200; tolerance=1e-5; maxit=50; pi.x <- 0.5; pi.z1.R <- .5; pi.z1.NR <- .5; pi.z2.R=.5; pi.z2.NR <- .5
Tcheck=5; adt=c(3,20); decimal=2; bshape1=3/4; bscale1=4; bshape2=3; bscale2=4
tbetas <- c(-1.5, 0.6, -0.5, -1.3, -1.2, -0.7, -1, 1.8, -2.4, -1, -0.5)
thetaZ=c('V','X'); gammaZ=c('V','X'); etaZ=c('1','X','Zr','X1Zr1','X1Zr1V1'); muZ=c('1','X')
np.t=length(thetaZ); np.g=length(gammaZ); np.e=length(etaZ); np.m=length(muZ); n.p=np.t+np.g+np.e+np.m 
TimeRDu = round(seq(0, adt[2], by=10^(-decimal)),decimal)
n_iter = 2 # number of iterations, set to 1000 in the manuscript
resp_cen_prob=0.15

result_table <- matrix(NA,nrow=length(tbetas),ncol=3*n_iter)
for (iter in 1:n_iter){
  ########################################################## Simulate data ####
  #### Generation of SMART with independent censoring
  simudata <- simSMART(n, tbetas, pi.x, pi.z1.R, pi.z1.NR, pi.z2.R, pi.z2.NR, 
                       Tcheck, bshape1, bscale1, bshape2, bscale2, adt, decimal)
  #### Generation of SMART with dependent censoring
  # simudata <- simSMART_two_censoring(n, tbetas, pi.x, pi.z1.R, pi.z1.NR, pi.z2.R, pi.z2.NR,
  #                                    Tcheck, bshape1, bscale1, bshape2, bscale2, adt, resp_cen_prob, 
  #                                    decimal)
  Tind <- order(simudata$T1,simudata$T2)
  sorted <- simudata[Tind,]
  TimeRDu = sort(unique(c(simudata$T1,simudata$T2)))
  len.et=length(TimeRDu)
  T1match= match(sorted$T1,TimeRDu)
  Tnr = which(TimeRDu==Tcheck); if(length(Tnr)==0){Tnr=len.et}
  
  thetaZs <- as.matrix(subset(sorted,select=thetaZ)) # covariate set for theta
  gammaZs <- as.matrix(subset(sorted,select=gammaZ)) # covariate set for gamma
  etaZs <- as.matrix(subset(sorted,select=etaZ)) # covariate set for eta
  muZs <- as.matrix(subset(sorted,select=muZ)) # covariate set for mu
  
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
  
  ########################################################## RUN MLE (may take a couple minutes for n=500)####
  test <- optim(par=rep(0.5,n.p),fn=loglik,gr=NULL, method = "L-BFGS-B", 
                lower=rep(-10,n.p), upper=rep(10,n.p), hessian=T, control=list(maxit=1000,trace=1))
  hessian.mat=test$hessian
  p.se <- sqrt(abs(diag(solve(-hessian.mat))))
  para.b <- test$p
  par.var.b <- (p.se)^2 
  para.true.b <- rep(0,length(para.b)) 
  p.value.b <- 2*pnorm(abs(para.b - para.true.b)/sqrt(abs(par.var.b)),lower.tail=F) 
  result_table[,(3 * iter -2):(3 * iter)] <- matrix(round(c(para.b,p.se,p.value.b),3),byrow = F)
  cat(iter)
}

result <- result_table %>% as.tibble
beta_estimate <- result %>% select(num_range("V",seq(1,ncol(.),3))) %>% apply(.,1,mean)
sd_beta <- result %>% select(num_range("V",seq(1,ncol(.),3))) %>% apply(.,1,sd)
mean_se <- result %>% select(num_range("V",seq(2,ncol(.),3))) %>% apply(.,1,mean)
beta_each <- result %>% select(num_range("V",seq(1,ncol(.),3)))
se_each <- result %>% select(num_range("V",seq(2,ncol(.),3)))
lower_bound <- beta_each - qnorm(0.975) * se_each
upper_bound <- beta_each + qnorm(0.975) * se_each
coverage <- (lower_bound <= tbetas) * (upper_bound >= tbetas)
CP <- apply(coverage,1,mean)

table1 <- data.frame(bias = beta_estimate-tbetas,
                           sd_beta = sd_beta,
                           mean_se = mean_se,
                           CP = CP) %>% round(.,3)

