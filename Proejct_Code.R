## table 2
y = c(46.79,46.54,95.82,95.57,201.48,201.28,471.19,469.27,
      602.63,598.54,696.43,691.17,773.07,744.45,835.45,805.88)
x = c(3.17,3.48,3.56,3.86,7.14,7.39,101.27,103.65,
      281.47,286.56,637.41,643.96,1126.94,1162.55,1725.30,1761.88)
x.log = log(x)
y.log=log(y)
Ads.data <- list(y,x)
Ads.data

## figure 1
plot(x,y, xlab = "PVA,mg/l",pch=20, ylab = "Amount Adsorb")
lines(x,y,lty = 1)
M2<-lm(y ~ x.log)
summary(M2)
a = -64.60
b = 117.64
y.M2 = a+b*x.log
points(x,y.M2, col = "blue", lty = 2,pch=20)
lines(x,y.M2, col = "blue", lty = 2)

# langmuir linearization
y.l=x/y
M1<-lm(y.l~x)
summary(M1)
b.l=1/(summary(M1)$coe[2])
a.l=(summary(M1)$coe[2])/(summary(M1)$coe[1])
y.M1=(a.l*b.l*x)/(1+a.l*x)
points(x,y.M1,col="red3",pch=20)
lines(x,y.M1,col="red",lty=3)
legend("topleft",legend=c("y","Langmuir","y=a+b*log(x)")
       ,col = c("black","red","blue"),lty=c(1,2,3),cex = 1)

##################################################################################
library(stats4)
#table 4
# MLE Model1
MLE1<-function(theta){ #theta=c(alpha_star, beta_star, sigma^2)
  mle<- -sum(-log(theta[3])-(y.log-theta[1]-theta[2]-x.log+log(1+exp(theta[1])*x))^2/(2*theta[3]^2))
  }

est=optim(MLE1,p=c(-4,7,1),hessian = TRUE)
fisher_info=est$hessian
se=sqrt(diag(solve(fisher_info)))
upper=est$par+1.96*se
upper=exp(upper) #0.0443002 854.8947926   1.3597738
lower=est$par-1.96*se
lower=exp(lower) #0.02475629 625.14629018   1.16101469

# MLE Model 2 
MLE2<-function(theta){
    mle<- -sum(-log(theta[3])-(y-theta[1]-theta[2]*log(x))^2/(2*theta[3]^2))
}
est=optim(theta<-c(-60,100,20),MLE2,hessian = TRUE)
fisher_info=est$hessian
se=sqrt(diag(solve(fisher_info)))
upper=est$par+1.96*se # -45.54711 121.37915  25.12492
lower=est$par-1.96*se # -83.63541 113.90057  12.19460

##############################################
x=c(3.17,3.48,3.56,3.86,7.14, 7.39, 101.27, 103.65, 281.47, 286.56, 
    637.41, 643.96, 1126.94, 1162.55, 1725.3, 1761.88)
y=c(46.79, 46.54, 95.82, 95.57, 201.48, 201.28, 471.19, 469.27, 602.63, 
    598.54, 696.43, 691.17, 773.07, 744.45, 835.45, 805.88)

## Table 4
## Random Walk Metropolis Hastings for Model 1
library(mgcv)
library(invgamma)
library(MCMCpack)
n=length(y)
y.log=log(y)
x.log=log(x)
N=10000 #iteration time
Theta=NULL
theta=c(2,2) #initialize
V=matrix(c(0.1,0,0,0.1),2,2) #covariance matrix
Sigma2=NULL
sigma2=2
set.seed(1)

accept=0
rate=0
post<-function(alpha,beta,sigma2){ # posterior function
  posterior<- exp(-sum((y.log-alpha-beta-x.log+log(1+exp(alpha)*x))^2)/(2*sigma2))
  return(posterior)
}

for (i in 1:N){
  
  #update alpha beta
  theta.new=rmvn(1,theta,V)
  alpha = post(theta.new[1],theta.new[2],sigma2)/post(theta[1],theta[2],sigma2)
  alpha = min(1,alpha)
  u=runif(1)
  if (u < alpha) {
    theta = theta.new
    accept=accept+1}
  Theta=rbind(Theta,theta)
  
  ## update sigma
  a=n/2
  b=sum((y.log-theta[1]-theta[2]-x.log+log(1+exp(theta[1])*x))^2)/2
  sigma2=rinvgamma(1,shape=a,scale=1/b)
  Sigma2=rbind(Sigma2,sigma2)
}
rate=accept/N 
rate

## Acf
par(mfrow=c(1,3))
acf(Theta[-seq(1,2000),1])
acf(Theta[-seq(1,2000),2])
acf(Sigma2[-seq(1,2000)])

Theta=cbind(exp(Theta[-seq(1,2000),1]),exp(Theta[-seq(1,2000),2]),Sigma2[-seq(1,2000)])
N=dim(Theta)[1]

## trace plot
plot(c(1:N),Theta[1:N,1],type = "l",main = "alpha",xlab = "iteration times")
plot(c(1:N),Theta[1:N,2],type = "l",main = "beta",xlab = "iteration times")
plot(c(1:length(Sigma2)),sqrt(Sigma2),type = "l",main = "sigma",xlab = "iteration times")

## marginal posterior distributions
par(mfrow=c(1,3))
plot(density(Theta[1:N,1]),main = "marginal posterior of alpha")
plot(density(Theta[1:N,2]),main = "marginal posterior of beta")
plot(density(sqrt(Theta[1:N,3])),main = "marginal posterior of sigma",xlim=c(0,3))

## 95% percentile intervals  
alpha95=quantile(Theta[,1],c(0.025,0.975)) #0.0238,0.0472
beta95=quantile(Theta[,2],c(0.025,0.975)) #612,876
sigma95=quantile(sqrt(Theta[,3]),c(0.025,0.975)) #0.179 0.386

## posterior mean
mean(Theta[,1]) #0.0337
mean(Theta[,2]) #734.4698
mean(sqrt(Sigma2)) #0.2586

######################################################################
y=c(46.79, 46.54, 95.82, 95.57, 201.48, 201.28, 471.19, 469.27, 602.63, 
    598.54, 696.43, 691.17, 773.07, 744.45, 835.45, 805.88)
x=c(3.17,3.48,3.56,3.86,7.14, 7.39, 101.27, 103.65, 281.47, 286.56, 
    637.41, 643.96, 1126.94, 1162.55, 1725.3, 1761.88)
n=length(y)
logx=log(x)
logy=log(y)
S=10000

##Gibbs sampler for Model 2
PHI=matrix(nrow=S,ncol=3)
PHI[1,]=phi=c(0,0,var(y))
for (s in 2:S) {
  mu1=sum(y-phi[2]*logx)/n
  t1=phi[3]/n
  phi[1]=rnorm(1,mu1,sqrt(t1))
  
  mu2=sum((y-phi[1])*logx)/sum(logx^2)
  t2=phi[3]/sum(logx^2)
  phi[2]=rnorm(1,mu2,sqrt(t2))
  
  t3=1/2*sum((y-phi[1]-phi[2]*logx)^2)
  phi[3]=rinvgamma(1,n/2,t3)
  
  PHI[s,]=phi
}
alpha=PHI[,1]
beta=PHI[,2]
sigma2=PHI[,3]

mean(alpha)
mean(beta)
mean(sqrt(sigma2))
quantile(alpha,c(0.025,0.975))
quantile(beta,c(0.025,0.975))
quantile(sqrt(sigma2),c(0.025,0.975))

#Convergence
par(mfrow=c(1,3))
acf(alpha)
acf(beta)
acf(sqrt(sigma2))

effectiveSize(alpha)
effectiveSize(beta)
effectiveSize(sqrt(sigma2))

######################################################################
#predictive values for MODEL 2
set.seed(1111)
density_m2=function(m,n){
  dnorm(y[n],alpha[m]+beta[m]*log(x[n]),sqrt(sigma2[m]))
}

weights=rep(0,S)
pred_y_m2=NULL
for (i in 1:n) {
  for (j in 1:S) {
    weights[j]=density_m2(j,i)
  }
  w=weights/sum(weights)
  phi_star=w%*%PHI
  pred_y_m2=rbind(pred_y_m2,rnorm(5000,mean=phi_star[1]+phi_star[2]*logx[i],sd=sqrt(phi_star[3])))
}

PI_M2=NULL
for (i in 1:n) {
  PI_M2=rbind(PI_M2,quantile(pred_y_m2[i,],c(0.025,0.25,0.75,0.975)))
}
PI_M2=cbind(y,PI_M2)
Ind1=ifelse(PI_M2[,1]>PI_M2[,3]&PI_M2[,1]<PI_M2[,4],1,0)
Ind2=ifelse(PI_M2[,1]>PI_M2[,2]&PI_M2[,1]<PI_M2[,5],1,0)
PI_M2=cbind(PI_M2,Ind1,Ind2)

library(xtable)
xtable(PI_M2)

#predictive values for MODEL 1
set.seed(111)
trace_alpha=Theta[,1]
trace_beta=Theta[,2]
trace_sigma_2=Theta[,3]

density_m1=function(m,n){
  dnorm(logy[n],log(trace_alpha[m]*trace_beta[m]*x[n])-log(1+trace_alpha[m]*x[n]),sqrt(trace_sigma_2[m]))
}

weights=rep(0,N)
pred_y_m1=NULL
for (i in 1:n) {
  for (j in 1:N) {
    weights[j]=density_m1(j,i)
  }
  w=weights/sum(weights)
  phi_star=w%*%Theta
  pred_y_m1=rbind(pred_y_m1,rnorm(5000,log(phi_star[1]*phi_star[2]*x[i])-log(1+phi_star[1]*x[i]),sd=sqrt(phi_star[3])))
}

PI_M1=NULL
for (i in 1:n) {
  PI_M1=rbind(PI_M1,quantile(pred_y_m1[i,],c(0.025,0.25,0.75,0.975)))
}
PI_M1=cbind(y,exp(PI_M1))
Ind1=ifelse(PI_M1[,1]>PI_M1[,3]&PI_M1[,1]<PI_M1[,4],1,0)
Ind2=ifelse(PI_M1[,1]>PI_M1[,2]&PI_M1[,1]<PI_M1[,5],1,0)
PI_M1=cbind(PI_M1,Ind1,Ind2)

library(xtable)
xtable(PI_M1)

######################################################################
#Deviance vs IQR plot
deviance_m1=abs(y-apply(pred_y_m1,1,median))
deviance_m2=abs(y-apply(pred_y_m2,1,median))
plot(PI_M1[,4]-PI_M1[,3],deviance_m1,type="p",pch=1,
     xlab="IQR",ylab="Deviance",
     ylim = c(0,850))
lines(PI_M2[,4]-PI_M2[,3],deviance_m2,type="p",pch=2)
legend("topleft",legend = c("Model 1","Model 2"),pch=c(1,2))

#dr Table
density=function(m,n){dnorm(y[n],alpha[1:m]+beta[1:m]*log(x[n]),sqrt(sigma2[1:m]))}
dr2=NULL
for (i in 1:n) {
  dr2[i]=S/(sum(1/density(S,i)))
}
dr2

PHI2=matrix(nrow=S,ncol=3)
alpha0=rnorm(1,0.0339,0.0059)
beta0=rnorm(1,732.1175,67.3516)
sigma0=var(logy)
PHI2[1,]=phi2=c(alpha0,beta0,sigma0)
for (s in 2:S) {
  phi2[1]=rnorm(1,0.0339,0.0059)
  
  phi2[2]=rnorm(1,732.1175,67.3516)
  
  t4=1/2*sum((y-log(phi2[1])-log(phi2[2])-logx+log(1+phi2[1]*x))^2)
  phi2[3]=rinvgamma(1,n/2,t4)
  
  PHI2[s,]=phi2
}
alpha2=PHI2[,1]
beta2=PHI2[,2]
sigma22=PHI2[,3]

mean(alpha2)
mean(beta2)
mean(sigma22)

density2=function(m,n){dnorm(y[n],log(alpha2[2:m])+log(beta2[2:m])+log(x[n])-log(1+alpha2[2:m]*x[n]),sqrt(sigma22[2:m]))}
dr1=NULL
for (i in 1:n) {
  dr1[i]=(S-1)/(sum(1/density2((S-1),i)))
}
dr1

lr=log(dr1,10)-log(dr2,10)

#CPO plot
obs=1:n
plot(obs,dr2,pch=2,ylim=c(0,0.02),xlab="Observation Number",ylab="CPO")
points(obs,dr1)
legend("topleft",c("cpo1","cpo2"),pch=c(1,2)) 

#lr plot
plot(obs,lr,ylim=c(-2,0.5),xlab="Observation Number",ylab="lr")
abline(h=0)
