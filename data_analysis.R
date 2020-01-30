library(gdata)
data_surv <- read.xls(file.path(yourpath,"calgbdata.xls"))
source('helper.R') # Misc. helper functions
source('calH.R') # Estimate 2 non-parametric baseline hazards for response and death
source('loglik.R') # Estimate log-likelihood

Tcheck <- 160      # If not respond until Tcheck, the patient is considered to be a non-responder.
sex <- data_surv$sex - 1 # 0=M,1=F
trt1 <- data_surv$trt1 - 1   # First stage treatment
R <- ifelse(data_surv$resp==1,1,0)  # Indicator of response
NR <- ifelse(data_surv$resp==0 & data_surv$survtime>=Tcheck,1,0) # Indicator of non-response and those who still survives at Tcheck
ED <- ifelse(data_surv$resp==0 & data_surv$survtime<Tcheck,1,0)  # early death = death before Tcheck
T1 <- round(ifelse(R==1,data_surv$resp_time/30,ifelse(NR==1,Tcheck/30,data_surv$survtime/30)),2)  
Zr <- ifelse(R==1 & data_surv$consent==1,data_surv$trt2,NA)  # Responders to first stage trt who consent to continue the trial will receive second stage trt
T2 <- round(ifelse(R==1 & data_surv$consent==0,T1,data_surv$survtime/30),2)
Death <- ifelse(R==1 & data_surv$consent==0,0,data_surv$death)  # For those who do not consent to continue, we simply count them as censored patients.
trt1trt2 <- trt1 * Zr
trt1trt2sex <- trt1 * Zr * sex
calgb_dat <- data.frame(X1=1,V=sex,X=trt1,X1Zr1=trt1trt2,X1Zr1V1=trt1trt2sex,R=R,NR=NR,D=ED,T1=T1,Zr=Zr,T2=T2,C2=Death)
colnames(calgb_dat)[1] <- "1"
thetaZ=c('V','X'); gammaZ=c('V','X'); etaZ=c('1','X','Zr','X1Zr1','X1Zr1V1'); muZ=NULL
np.t=length(thetaZ); np.g=length(gammaZ); np.e=length(etaZ); np.m=length(muZ); n.p=np.t+np.g+np.e+np.m;tolerance<-1e-5;maxit<-50

#### analysis ####
Tind <- order(calgb_dat$T1,calgb_dat$T2)
sorted <- calgb_dat[Tind,]
TimeRDu = sort(unique(c(calgb_dat$T1,calgb_dat$T2)))
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

n <- nrow(calgb_dat)
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
result <- round(data.frame(para.b,sqrt(abs(par.var.b)),p.value.b),3)
names(result) <- c('est','se','pvalue')
result

H.fit=calH(para.b)
H.fit.1=H.fit$H1 
H.fit.2=H.fit$H2 
Zt.V=para.b[1]; Zt.X=para.b[2]
Zg.V=para.b[3]; Zg.X=para.b[4]
Ze.1=para.b[5]; Ze.X=para.b[6]; Ze.Zr=para.b[7]; Ze.X1Zr1=para.b[8]; Ze.X1Zr1V1=para.b[9]
Zm.1=Zm.X=Zm.Zr=0



########################### Time-to-death: Joint model and WRSE (Guo&Tsiatis,2005)  ####
## Joint model 
fit.V0A1B1=sapply(TimeRDu, FUN=s.km2.pred, theta=exp(0), gamma=exp(0), eta=exp(Ze.1), mu=exp(Zm.1))
fit.V0A1B2=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(0), gamma=exp(0), eta=exp(Ze.1+Ze.Zr), mu=exp(Zm.1))
fit.V0A2B1=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(0+Zt.X), gamma=exp(0+Zg.X), eta=exp(Ze.1+Ze.X), mu=exp(Zm.1+Zm.X))
fit.V0A2B2=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(0+Zt.X), gamma=exp(0+Zg.X), eta=exp(Ze.1+Ze.X+Ze.Zr+Ze.X1Zr1), mu=exp(Zm.1+Zm.X))
fit.V1A1B1=sapply(TimeRDu, FUN=s.km2.pred, theta=exp(Zt.V), gamma=exp(Zg.V), eta=exp(Ze.1), mu=exp(Zm.1))
fit.V1A1B2=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V), gamma=exp(Zg.V), eta=exp(Ze.1+Ze.Zr), mu=exp(Zm.1))
fit.V1A2B1=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V+Zt.X), gamma=exp(Zg.V+Zg.X), eta=exp(Ze.1+Ze.X), mu=exp(Zm.1+Zm.X))
fit.V1A2B2=sapply(TimeRDu,FUN=s.km2.pred, theta=exp(Zt.V+Zt.X), gamma=exp(Zg.V+Zg.X), eta=exp(Ze.1+Ze.X+Ze.Zr+Ze.X1Zr1+Ze.X1Zr1V1), mu=exp(Zm.1+Zm.X))
savefit = cbind(TimeRDu, fit.V0A1B1, fit.V0A1B2, fit.V0A2B1, fit.V0A2B2,fit.V1A1B1, fit.V1A1B2, fit.V1A2B1, fit.V1A2B2)

## Survival estimator (Guo&Tsiatis) from package 'DTR':
library(DTR)
calgb_dat$TR = 0; calgb_dat$TR[calgb_dat$R==1]=calgb_dat$T1[calgb_dat$R==1]
calgb_dat$Z=calgb_dat$Zr; calgb_dat$Z[is.na(calgb_dat$Z)]=0
calgb_dat$U=calgb_dat$T2
calgb_dat$delta=calgb_dat$C2
calgb_dat2=subset(calgb_dat, select=c(V,X,Z,TR,R,U,delta))
est01 = WRSEestimate(data=calgb_dat2);


calgb_long <- savefit %>% 
  as.tibble() %>% 
  dplyr::rename(time = TimeRDu) %>% 
  gather(key=DTR,value=survival,fit.V0A1B1:fit.V1A2B2) %>% 
  mutate(V = ifelse(str_detect(DTR,"V1"),1,0),
         DTR = ifelse(str_detect(DTR,"A1B1"),"A1B1C1",
                      ifelse(str_detect(DTR,"A1B2"),"A1B2C1",
                             ifelse(str_detect(DTR,"A2B1"),"A2B1C1","A2B2C1"))),
         model = "joint")

gender_labels <- c("male", "female")
names(gender_labels) <- c(0, 1)

ggplot(calgb_long) + 
  geom_line(aes(x = time, y = survival, group = DTR, color = DTR)) + 
  scale_color_manual(breaks = c("A1B1C1", "A1B2C1", "A2B1C1", "A2B2C1"),
                     values = c("grey80", "grey60", "grey40", "grey20")) + 
  facet_grid(cols = vars(V), labeller = labeller(V = gender_labels)) + 
  xlim(0, 30) + 
  ggtitle("Predicted survival rates for males versus females using joint model") + 
  labs(x = "Month", y = "Fraction Surviving") +
  theme_bw()



calgb_wrse_long <- cbind(est01$time,est01$SURV11,est01$SURV12,est01$SURV21,est01$SURV22) %>% 
  as.tibble() %>% 
  dplyr::rename(time = V1, A1B1 = V2, A1B2 = V3, A2B1 = V4, A2B2 = V5) %>% 
  gather(key=DTR,value=survival,A1B1:A2B2) %>% 
  mutate(A = ifelse(str_detect(DTR,"A2"),2,1),
         B = ifelse(str_detect(DTR,"B2"),2,1), 
         model = "Guo")

ggplot(calgb_wrse_long) + 
  geom_line(aes(x = time, y = survival, group = DTR, color = DTR)) + 
  scale_color_manual(breaks = c("A1B1", "A1B2", "A2B1", "A2B2"),
                     values = c("grey80", "grey60", "grey40", "grey20"),
                     labels = c("A1B1C1", "A1B2C1", "A2B1C1", "A2B2C1")) + 
  xlim(0, 30) + 
  ggtitle("Predicted survival rates for all using WRSE") + 
  labs(x = "Month", y = "Fraction Surviving") + 
  theme_bw()

  
  
  