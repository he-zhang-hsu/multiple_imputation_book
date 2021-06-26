# Example 8.7

library(coda)
library(lattice)
library(Matrix)
library(mice)
library(mitools)
library(MASS)
library(survival)
###########################################################################################################
# Data Generating Function
###########################################################################################################
miaftA=function(N,p,alpha,beta1,beta2,gamma0,a,b,c,d,e)
{
## Complete data generation!  
Z = rbinom(N,1,p)
Z = sort(Z)                # To make the following generation easier
N0 = length(Z[Z==0])
N1 = length(Z[Z==1])

# Specify the parameters in the Normal distribution.
mu1 = c(1,1)
sigma1 = matrix(c(1,0.5,0.5,1),nrow=2)   # A general way to calculate the square root of a matrix
eigen1 = eigen(sigma1)
V1 = eigen1$vectors
sigma1_root =  V1%*%diag(sqrt(eigen1$values)) %*% t(V1)
mu0 = c(0,0)
sigma0 = matrix(c(1,0.5,0.5,1),ncol=2)
eigen0 = eigen(sigma0)
V0 = eigen0$vectors
sigma0_root =  V0%*%diag(sqrt(eigen0$values)) %*% t(V0)

# Generate data of X.
x01 = rnorm(N0)            # Start to generate (X1,X2) for Z=0
x02 = rnorm(N0)
x0 = cbind(x01,x02) %*% t(sigma0_root) + mu0
x11 = rnorm(N1)            # Start to generate (X1,X2) for Z=1
x12 = rnorm(N1)
x1 = cbind(x11,x12) %*% t(sigma1_root) + mu1
X = rbind(x0,x1)

# Generate error items.
epsilon = rnorm(N)

# Generate the response.
Ti = exp(alpha + beta1*X[,1]+beta2*X[,2]+gamma0*Z+epsilon)

# Decide the censoring.
# Ti is the complete survival data. Yi is censored data with indication delta.
Ci = runif(N,min=0,max=a)
Yi = Ti
delta = rep(1,N)
for(i in 1:N)
{
if(Ti[i]>Ci[i])
  {
    Yi[i] = Ci[i]
    delta[i] = 0
  }
}
logYi = log(Yi)    # logYi are the transformed observed time.

# use following to control the censoring case.
Yi2 = log(Ti)
Yi2[which(Ti>Ci)] = NA
Yi2.cen = rep(-100,N)
Yi2.cen[which(Ti>Ci)] = log(Ci)[which(Ti>Ci)]

#####Generate missing data under MAR for X1 and Z.
# Decide the missing data for X1 and Z
x1missrate = exp(b+c*X[,2]+d*logYi+e*delta)/(1+exp(b+c*X[,2]+d*logYi+e*delta))
U1 = runif(N)
sele.indX1 <- rep(1, N) 
sele.indX1[ x1missrate>U1 ] <- 0

zmissrate = c()
for(i in 1:N)
{
if (is.na(X[i,1]))	zmissrate[i] = exp(b+c*X[i,2]+d*logYi[i]+e*delta[i])/(1+exp(b+c*X[i,2]+d*logYi[i]+e*delta[i]))
if (!is.na(X[i,1]))	zmissrate[i] = exp(b+c*X[i,1]+c*X[i,2]+d*logYi[i]+e*delta[i])/(1+exp(b+c*X[i,1]+c*X[i,2]+d*logYi[i]+e*delta[i]))
}
U2 = runif(N)
sele.indZ <- rep(1, N) 
sele.indZ[ zmissrate>U2 ] <- 0

X1=X[,1]
X2=X[,2]
data <- data.frame( logYi,delta,Yi2,Yi2.cen,X1,X2,Z, sele.indX1,sele.indZ ) 
return( data )
}

###########################################################################################################
# Main Function for Imputation, Analysis and Pooling 
###########################################################################################################
mainfmiaft <-  function( logYi,delta,Yi2,Yi2.cen,X1,X2,Z, sele.indX1,sele.indZ, distribution )
  {
 ###############################################
 # The full-cohort data analysis
 ###############################################
DATAF = data.frame(logYi,delta,X1,X2,Z)
FullAnal = survreg(Surv(exp(logYi),delta)~X1+X2+Z,dist = distribution,data=DATAF, score=TRUE)
FullResult = data.frame(summary(FullAnal)$table[,1:2],Fulllow=summary(FullAnal)$table[,1]-1.96*summary(FullAnal)$table[,2],Fullupp=summary(FullAnal)$table[,1]+1.96*summary(FullAnal)$table[,2])

if (dim(FullResult)[1]==0 | sum(is.na(FullAnal$score))>0)
{RFULLb1=t(c(NA,NA,NA,NA))
 RFULLb2=t(c(NA,NA,NA,NA))
 RFULLg=t(c(NA,NA,NA,NA))} else 
{RFULLb1 = FullResult[2,]
 RFULLb2 = FullResult[3,]
 RFULLg = FullResult[4,]}
write.table( RFULLb1, "RFULLb1.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( RFULLb2, "RFULLb2.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( RFULLg, "RFULLg.txt", append=T, row.names=FALSE, col.names=FALSE )
 

###############################################
 # The complete data analysis
 ###############################################
X1[sele.indX1==0] = NA
Z[sele.indZ==0] = NA
DATAC = data.frame(logYi,delta,X1,X2,Z)
CompAnal = survreg(Surv(exp(logYi),delta)~X1+X2+Z,dist = distribution,data=DATAC, score=TRUE)
CompResult = data.frame(summary(CompAnal)$table[,1:2],Complow= summary(CompAnal)$table[,1]-1.96*summary(CompAnal)$table[,2],Compupp=summary(CompAnal)$table[,1]+1.96*summary(CompAnal)$table[,2])

if (dim(CompResult)[1]==0 | sum(is.na(CompAnal$score))>0)
{RCOMPb1=t(c(NA,NA,NA,NA))
 RCOMPb2=t(c(NA,NA,NA,NA))
 RCOMPg=t(c(NA,NA,NA,NA))} else 
{RCOMPb1 = CompResult[2,]
 RCOMPb2 = CompResult[3,]
 RCOMPg = CompResult[4,]}
write.table( RCOMPb1, "RCOMPb1.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( RCOMPb2, "RCOMPb2.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( RCOMPg, "RCOMPg.txt", append=T, row.names=FALSE, col.names=FALSE )
 ###############################################
 # RMICE ---PMM
 ###############################################
m = 20
vars=list()
Z_factor=as.factor(Z);
Data1 = data.frame(logYi,delta,X1,X2,Z_factor)
imp <- mice(Data1,m= m,maxit = 20, meth=c("","","pmm","","logreg"),seed=100);
mit=imputationList(lapply(1:m,complete,data=imp))
models=with(mit,survreg(Surv(exp(logYi),delta)~X1+X2+as.numeric(Z_factor),dist = distribution, score=TRUE))
div=vector()
for(i in 1:m)
{ div[i]=sum(is.na(models[[i]]$score)) }
if ( sum(div)==0 ){
betas<-MIextract(models,fun=coef)
for(i in 1:m)
{vars[[i]]<-MIextract(models,fun=vcov)[[i]][1:4,1:4]}
s=summary(MIcombine(betas,vars))
beta.mi=s$results
se.mi=s$se
lower.mi=s$"(lower"
upper.mi=s$"upper)"} else {
beta.mi=c(NA,NA,NA,NA)
se.mi=c(NA,NA,NA,NA)
lower.mi=c(NA,NA,NA,NA)
upper.mi=c(NA,NA,NA,NA)}

RMICEb1 = c( beta.mi[2],se.mi[2],lower.mi[2],upper.mi[2] )
RMICEb2 = c( beta.mi[3],se.mi[3],lower.mi[3],upper.mi[3] )
RMICEg = c( beta.mi[4],se.mi[4],lower.mi[4],upper.mi[4] )
write.table( t(RMICEb1), "RMICEpmmb1.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( t(RMICEb2), "RMICEpmmb2.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( t(RMICEg), "RMICEpmmg.txt", append=T, row.names=FALSE, col.names=FALSE )
 ###############################################
 # RMICE ---PMM -SEP
 ###############################################
#Data1 = data.frame(logYi,delta,X1,X2,Z)
# separate data1 into two subsets based on delta;
Data1_delta_1=Data1[delta==1,];
Data1_delta_0=Data1[delta==0,];
#m = 20
# imp <- mice(Data1,m= m,maxit = 20, meth=c("","","pmm","","logreg","","",""));
# predictormat is the matrix specifying which variables in the dataset are to be used as predictors in MICE
# If the [I,j]a€?th element is 1, it means the j-th variable is used as a predictor for the ith variable;
# if the [I,j]a€?th element is 0, it means that j-th variable is not used as a predictor for the ith variable.
# because delta is constant in each of the dataset, it shouldna€?t be used as the predictor;
predictormat=1-diag(1,ncol(Data1));
# exclude delta as the predictor in MICE 
predictormat[,2]=0;
if (sum(table(Data1_delta_1[,5])==0)==0 & sum(table(Data1_delta_0[,5])==0)==0){
# separate imputation;
Imp_delta_1 <- mice(Data1_delta_1,m= m,maxit = 20, meth=c("","","pmm","","logreg"), predictorMatrix=predictormat, seed=100);
Imp_delta_0 <- mice(Data1_delta_0,m= m,maxit = 20, meth=c("","","pmm","","logreg"), predictorMatrix=predictormat, seed=100);
imp=rbind(Imp_delta_1, Imp_delta_0)
mit=imputationList(lapply(1:m,complete,data=imp))
models=with(mit,survreg(Surv(exp(logYi),delta)~X1+X2+as.numeric(Z_factor),dist = distribution, score=TRUE))
div=vector()
for(i in 1:m)
{ div[i]=sum(is.na(models[[i]]$score)) }
if ( sum(div)==0 ){
betas<-MIextract(models,fun=coef)
vars=list()
for(i in 1:m)
{vars[[i]]<-MIextract(models,fun=vcov)[[i]][1:4,1:4]}
s=summary(MIcombine(betas,vars))
beta.mi=s$results
se.mi=s$se
lower.mi=s$"(lower"
upper.mi=s$"upper)"}} else {
beta.mi=c(NA,NA,NA,NA)
se.mi=c(NA,NA,NA,NA)
lower.mi=c(NA,NA,NA,NA)
upper.mi=c(NA,NA,NA,NA)}

RMICEb1 = c( beta.mi[2],se.mi[2],lower.mi[2],upper.mi[2] )
RMICEb2 = c( beta.mi[3],se.mi[3],lower.mi[3],upper.mi[3] )
RMICEg = c( beta.mi[4],se.mi[4],lower.mi[4],upper.mi[4] )
write.table( t(RMICEb1), "RMICEIpmmb1.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( t(RMICEb2), "RMICEIpmmb2.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( t(RMICEg), "RMICEIpmmg.txt", append=T, row.names=FALSE, col.names=FALSE )
 ###############################################
 # RMICE ---NORM
 ###############################################
imp <- mice(Data1,m= m,maxit =20, meth=c("","","norm","","logreg"),seed=100);
mit=imputationList(lapply(1:m,complete,data=imp))
models=with(mit,survreg(Surv(exp(logYi),delta)~X1+X2+as.numeric(Z_factor),dist = distribution, score=TRUE))
div=vector()
for(i in 1:m)
{ div[i]=sum(is.na(models[[i]]$score)) }
if ( sum(div)==0 ){
betas<-MIextract(models,fun=coef)
vars=list()
for(i in 1:m)
{vars[[i]]<-MIextract(models,fun=vcov)[[i]][1:4,1:4]}
s=summary(MIcombine(betas,vars))
beta.mi=s$results
se.mi=s$se
lower.mi=s$"(lower"
upper.mi=s$"upper)"} else {
beta.mi=c(NA,NA,NA,NA)
se.mi=c(NA,NA,NA,NA)
lower.mi=c(NA,NA,NA,NA)
upper.mi=c(NA,NA,NA,NA)}

RMICEb1 = c( beta.mi[2],se.mi[2],lower.mi[2],upper.mi[2] )
RMICEb2 = c( beta.mi[3],se.mi[3],lower.mi[3],upper.mi[3] )
RMICEg = c( beta.mi[4],se.mi[4],lower.mi[4],upper.mi[4] )
write.table( t(RMICEb1), "RMICEnormb1.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( t(RMICEb2), "RMICEnormb2.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( t(RMICEg), "RMICEnormg.txt", append=T, row.names=FALSE, col.names=FALSE )
 ###############################################
 # RMICE ---NORM -SEP
 ###############################################
if (sum(table(Data1_delta_1[,5])==0)==0 & sum(table(Data1_delta_0[,5])==0)==0){
# separate imputation;
Imp_delta_1 <- mice(Data1_delta_1,m= m,maxit = 20, meth=c("","","norm","","logreg"), predictorMatrix=predictormat, seed=100);
Imp_delta_0 <- mice(Data1_delta_0,m= m,maxit = 20, meth=c("","","norm","","logreg"), predictorMatrix=predictormat, seed=100);
# combine the imputations;
imp=rbind(Imp_delta_1, Imp_delta_0)
mit=imputationList(lapply(1:m,complete,data=imp))
models=with(mit,survreg(Surv(exp(logYi),delta)~X1+X2+as.numeric(Z_factor),dist = distribution, score=TRUE))
div=vector()
for(i in 1:m)
{ div[i]=sum(is.na(models[[i]]$score)) }
if ( sum(div)==0 ){
betas<-MIextract(models,fun=coef)
vars=list()
for(i in 1:m)
{vars[[i]]<-MIextract(models,fun=vcov)[[i]][1:4,1:4]}
s=summary(MIcombine(betas,vars))
beta.mi=s$results
se.mi=s$se
lower.mi=s$"(lower"
upper.mi=s$"upper)"}} else {
beta.mi=c(NA,NA,NA,NA)
se.mi=c(NA,NA,NA,NA)
lower.mi=c(NA,NA,NA,NA)
upper.mi=c(NA,NA,NA,NA)}

RMICEb1 = c( beta.mi[2],se.mi[2],lower.mi[2],upper.mi[2] )
RMICEb2 = c( beta.mi[3],se.mi[3],lower.mi[3],upper.mi[3] )
RMICEg = c( beta.mi[4],se.mi[4],lower.mi[4],upper.mi[4] )
write.table( t(RMICEb1), "RMICEInormb1.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( t(RMICEb2), "RMICEInormb2.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( t(RMICEg), "RMICEInormg.txt", append=T, row.names=FALSE, col.names=FALSE )
 ###############################################
 # RMICE ---PMM -H(t)
 ###############################################
Yi=exp(logYi)
Data1h = data.frame(Yi,delta,X1,X2,Z_factor)
Data1h$ch=nelsonaalen(Data1h, Yi, delta)
pMatrix=matrix(c(0,0,0,0,0,0, 0,0,0,0,0,0, 0,1,0,1,1,1, 0,0,0,0,0,0, 0,1,1,1,0,1, 0,0,0,0,0,0),6,6,byrow=T)
imp <- mice(Data1h,m= m,maxit = 20, meth=c("","","pmm","","logreg",""),predictorMatrix=pMatrix,seed=100)
mit=imputationList(lapply(1:m,complete,data=imp))
models=with(mit,survreg(Surv(Yi,delta)~X1+X2+as.numeric(Z_factor),dist = distribution, score=TRUE))
div=vector()
for(i in 1:m)
{ div[i]=sum(is.na(models[[i]]$score)) }
if ( sum(div)==0 ){
betas<-MIextract(models,fun=coef)
vars=list()
for(i in 1:m)
{vars[[i]]<-MIextract(models,fun=vcov)[[i]][1:4,1:4]}
s=summary(MIcombine(betas,vars))
beta.mi=s$results
se.mi=s$se
lower.mi=s$"(lower"
upper.mi=s$"upper)"} else {
beta.mi=c(NA,NA,NA,NA)
se.mi=c(NA,NA,NA,NA)
lower.mi=c(NA,NA,NA,NA)
upper.mi=c(NA,NA,NA,NA)}

RMICEb1 = c( beta.mi[2],se.mi[2],lower.mi[2],upper.mi[2] )
RMICEb2 = c( beta.mi[3],se.mi[3],lower.mi[3],upper.mi[3] )
RMICEg = c( beta.mi[4],se.mi[4],lower.mi[4],upper.mi[4] )
write.table( t(RMICEb1), "RMICEHpmmb1.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( t(RMICEb2), "RMICEHpmmb2.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( t(RMICEg), "RMICEHpmmg.txt", append=T, row.names=FALSE, col.names=FALSE )
 ###############################################
 # RMICE ---PMM -H(t) -SEP
 ###############################################
Data1_delta_1=Data1h[delta==1,];
Data1_delta_0=Data1h[delta==0,];
pMatrixI=matrix(c(0,0,0,0,0,0, 0,0,0,0,0,0, 0,1,0,1,1,1, 0,0,0,0,0,0, 0,1,1,1,0,1, 0,0,0,0,0,0),6,6,byrow=T)
pMatrixI[,2]=0;
if (sum(table(Data1_delta_1[,5])==0)==0 & sum(table(Data1_delta_0[,5])==0)==0){
Imp_delta_1 <- mice(Data1_delta_1,m= m,maxit = 20, meth=c("","","pmm","","logreg",""), predictorMatrix=pMatrixI, seed=100);
Imp_delta_0 <- mice(Data1_delta_0,m= m,maxit = 20, meth=c("","","pmm","","logreg",""), predictorMatrix=pMatrixI, seed=100);
imp=rbind(Imp_delta_1, Imp_delta_0)
mit=imputationList(lapply(1:m,complete,data=imp))
models=with(mit,survreg(Surv(Yi,delta)~X1+X2+as.numeric(Z_factor),dist = distribution, score=TRUE))
div=vector()
for(i in 1:m)
{ div[i]=sum(is.na(models[[i]]$score)) }
if ( sum(div)==0 ){
betas<-MIextract(models,fun=coef)
vars=list()
for(i in 1:m)
{vars[[i]]<-MIextract(models,fun=vcov)[[i]][1:4,1:4]}
s=summary(MIcombine(betas,vars))
beta.mi=s$results
se.mi=s$se
lower.mi=s$"(lower"
upper.mi=s$"upper)"}} else {
beta.mi=c(NA,NA,NA,NA)
se.mi=c(NA,NA,NA,NA)
lower.mi=c(NA,NA,NA,NA)
upper.mi=c(NA,NA,NA,NA)}

RMICEb1 = c( beta.mi[2],se.mi[2],lower.mi[2],upper.mi[2] )
RMICEb2 = c( beta.mi[3],se.mi[3],lower.mi[3],upper.mi[3] )
RMICEg = c( beta.mi[4],se.mi[4],lower.mi[4],upper.mi[4] )
write.table( t(RMICEb1), "RMICEHIpmmb1.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( t(RMICEb2), "RMICEHIpmmb2.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( t(RMICEg), "RMICEHIpmmg.txt", append=T, row.names=FALSE, col.names=FALSE )
 ###############################################
 # RMICE ---NORM -H(t)
 ###############################################
imp <- mice(Data1h,m= m,maxit =20, meth=c("","","norm","","logreg",""),predictorMatrix=pMatrix,seed=100)
mit=imputationList(lapply(1:m,complete,data=imp))
models=with(mit,survreg(Surv(Yi,delta)~X1+X2+as.numeric(Z_factor),dist = distribution, score=TRUE))
div=vector()
for(i in 1:m)
{ div[i]=sum(is.na(models[[i]]$score)) }
if ( sum(div)==0 ){
betas<-MIextract(models,fun=coef)
vars=list()
for(i in 1:m)
{vars[[i]]<-MIextract(models,fun=vcov)[[i]][1:4,1:4]}
s=summary(MIcombine(betas,vars))
beta.mi=s$results
se.mi=s$se
lower.mi=s$"(lower"
upper.mi=s$"upper)"} else {
beta.mi=c(NA,NA,NA,NA)
se.mi=c(NA,NA,NA,NA)
lower.mi=c(NA,NA,NA,NA)
upper.mi=c(NA,NA,NA,NA)}

RMICEb1 = c( beta.mi[2],se.mi[2],lower.mi[2],upper.mi[2] )
RMICEb2 = c( beta.mi[3],se.mi[3],lower.mi[3],upper.mi[3] )
RMICEg = c( beta.mi[4],se.mi[4],lower.mi[4],upper.mi[4] )
write.table( t(RMICEb1), "RMICEHnormb1.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( t(RMICEb2), "RMICEHnormb2.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( t(RMICEg), "RMICEHnormg.txt", append=T, row.names=FALSE, col.names=FALSE )
 ###############################################
 # RMICE ---NORM -H(t) -SEP
 ###############################################
if (sum(table(Data1_delta_1[,5])==0)==0 & sum(table(Data1_delta_0[,5])==0)==0){
Imp_delta_1 <- mice(Data1_delta_1,m= m,maxit = 20, meth=c("","","norm","","logreg",""), predictorMatrix=pMatrixI, seed=100);
Imp_delta_0 <- mice(Data1_delta_0,m= m,maxit = 20, meth=c("","","norm","","logreg",""), predictorMatrix=pMatrixI, seed=100);
imp=rbind(Imp_delta_1, Imp_delta_0)
mit=imputationList(lapply(1:m,complete,data=imp))
models=with(mit,survreg(Surv(Yi,delta)~X1+X2+as.numeric(Z_factor),dist = distribution, score=TRUE))
div=vector()
for(i in 1:m)
{ div[i]=sum(is.na(models[[i]]$score)) }
if ( sum(div)==0 ){
betas<-MIextract(models,fun=coef)
vars=list()
for(i in 1:m)
{vars[[i]]<-MIextract(models,fun=vcov)[[i]][1:4,1:4]}
s=summary(MIcombine(betas,vars))
beta.mi=s$results
se.mi=s$se
lower.mi=s$"(lower"
upper.mi=s$"upper)"}} else {
beta.mi=c(NA,NA,NA,NA)
se.mi=c(NA,NA,NA,NA)
lower.mi=c(NA,NA,NA,NA)
upper.mi=c(NA,NA,NA,NA)}

RMICEb1 = c( beta.mi[2],se.mi[2],lower.mi[2],upper.mi[2] )
RMICEb2 = c( beta.mi[3],se.mi[3],lower.mi[3],upper.mi[3] )
RMICEg = c( beta.mi[4],se.mi[4],lower.mi[4],upper.mi[4] )
write.table( t(RMICEb1), "RMICEHInormb1.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( t(RMICEb2), "RMICEHInormb2.txt", append=T, row.names=FALSE, col.names=FALSE )
write.table( t(RMICEg), "RMICEHInormg.txt", append=T, row.names=FALSE, col.names=FALSE )
  
}


###########################################################################################################
# Generate data 
###########################################################################################################
N = 500
p = 0.5
alpha=0
beta1=1.5
beta2=3
gamma0=1.5
a=150
b=0.8
c=1.3
d=-1.3
e=-3
set.seed(100)
TIMES = 1000

logYi.0 <- matrix( 0, nrow = N, ncol = TIMES )
Yi2.0 <- matrix( 0, nrow = N, ncol = TIMES )
delta.0 <- matrix( 0, nrow = N, ncol = TIMES )
Yi2.cen.0 <- matrix( 0, nrow = N, ncol = TIMES )
X1.0 <- matrix( 0, nrow = N, ncol = TIMES )
X2.0 <- matrix( 0, nrow = N, ncol = TIMES )
Z.0 <- matrix( 0, nrow = N, ncol = TIMES )
sele.indX1.0 <- matrix( 0, nrow = N, ncol = TIMES )
sele.indZ.0 <- matrix( 0, nrow = N, ncol = TIMES )
sele.ind.0 <- matrix( 0, nrow = N, ncol = TIMES )
sele.percX1.0 <- vector() 
sele.percZ.0 <- vector() 
sele.perc.0 <- vector()

for (k in 1:TIMES)
{
data <- miaftA(N,p,alpha,beta1,beta2,gamma0,a,b,c,d,e)
logYi.0[,k] <- data$logYi
delta.0[,k] <- data$delta
Yi2.0[,k] <- data$Yi2
Yi2.cen.0[,k] <- data$Yi2.cen
X1.0[,k] <- data$X1
X2.0[,k] <- data$X2
Z.0[,k] <- data$Z
sele.indX1.0[,k] <- data$sele.indX1
sele.indZ.0[,k] <- data$sele.indZ
sele.ind.0[,k] <- ifelse(data$sele.indX1==1 & data$sele.indZ==1, 1, 0)

sele.percX1.0[k] <- mean(sele.indX1.0[,k])
sele.percZ.0[k] <- mean(sele.indZ.0[,k])
sele.perc.0[k] <- mean(sele.ind.0[,k])
sele.indX1 <- sele.indX1.0[,k]
}

# selection percentages 
sele.percX1 <- mean(sele.percX1.0)
sele.percZ <- mean(sele.percZ.0)
sele.perc <- mean(sele.perc.0)
# case percentage
caseperc <- mean( apply( delta.0, 2, mean ) )

probab=c(round(caseperc, 2),round(sele.percX1, 2),round(sele.percZ, 2),round(sele.perc, 2))
write.csv(probab,"prob.csv")


###########################################################################################################
# Feed generated data to the main function to produce estimates for beta1, beta2, gamma0
###########################################################################################################
options("digits" = 4)

for( l in 1:TIMES )
  {
# obtain data from previously generated data vectors
    logYi <-  logYi.0[,l] 
    delta <- delta.0[,l]
    Yi2 <- Yi2.0[,l]
    Yi2.cen <- Yi2.cen.0[,l] 
    X1 <- X1.0[,l]
    X2 <- X2.0[,l]
    Z <- Z.0[,l]
    sele.indX1 <- sele.indX1.0[,l]
    sele.indZ <- sele.indZ.0[,l]
mainfmiaft( logYi,delta,Yi2,Yi2.cen,X1,X2,Z, sele.indX1,sele.indZ, "lognormal" )
  }


###########################################################################################################
# Retrieve the saved estimates, and summarize as reported in the paper
###########################################################################################################
RFULLb1 = read.table("RFULLb1.txt")
RFULLb2 = read.table("RFULLb2.txt")
RFULLg = read.table("RFULLg.txt")
if (length(which(is.na(RFULLb1[,1])))>0)
{RFULLb1 = RFULLb1[-which(is.na(RFULLb1[,1])),]
RFULLb2 = RFULLb2[-which(is.na(RFULLb2[,1])),]
RFULLg = RFULLg[-which(is.na(RFULLg[,1])),]}

RCOMPb1 = read.table("RCOMPb1.txt")
RCOMPb2 = read.table("RCOMPb2.txt")
RCOMPg = read.table("RCOMPg.txt")
if (length(which(is.na(RCOMPb1[,1])))>0)
{RCOMPb1 = RCOMPb1[-which(is.na(RCOMPb1[,1])),]
RCOMPb2 = RCOMPb2[-which(is.na(RCOMPb2[,1])),]
RCOMPg = RCOMPg[-which(is.na(RCOMPg[,1])),]}

RMICEpmmb1 = read.table("RMICEpmmb1.txt")
RMICEpmmb2 = read.table("RMICEpmmb2.txt")
RMICEpmmg = read.table("RMICEpmmg.txt")
if (length(which(is.na(RMICEpmmb1[,1])))>0)
{RMICEpmmb1 = RMICEpmmb1[-which(is.na(RMICEpmmb1[,1])),]
RMICEpmmb2 = RMICEpmmb2[-which(is.na(RMICEpmmb2[,1])),]
RMICEpmmg = RMICEpmmg[-which(is.na(RMICEpmmg[,1])),]}

RMICEIpmmb1 = read.table("RMICEIpmmb1.txt")
RMICEIpmmb2 = read.table("RMICEIpmmb2.txt")
RMICEIpmmg = read.table("RMICEIpmmg.txt")
if (length(which(is.na(RMICEIpmmb1[,1])))>0)
{RMICEIpmmb1 = RMICEIpmmb1[-which(is.na(RMICEIpmmb1[,1])),]
RMICEIpmmb2 = RMICEIpmmb2[-which(is.na(RMICEIpmmb2[,1])),]
RMICEIpmmg = RMICEIpmmg[-which(is.na(RMICEIpmmg[,1])),]}

RMICEnormb1 = read.table("RMICEnormb1.txt")
RMICEnormb2 = read.table("RMICEnormb2.txt")
RMICEnormg = read.table("RMICEnormg.txt")
if (length(which(is.na(RMICEnormb1[,1])))>0)
{RMICEnormb1 = RMICEnormb1[-which(is.na(RMICEnormb1[,1])),]
RMICEnormb2 = RMICEnormb2[-which(is.na(RMICEnormb2[,1])),]
RMICEnormg = RMICEnormg[-which(is.na(RMICEnormg[,1])),]}

RMICEInormb1 = read.table("RMICEInormb1.txt")
RMICEInormb2 = read.table("RMICEInormb2.txt")
RMICEInormg = read.table("RMICEInormg.txt")
if (length(which(is.na(RMICEInormb1[,1])))>0)
{RMICEInormb1 = RMICEInormb1[-which(is.na(RMICEInormb1[,1])),]
RMICEInormb2 = RMICEInormb2[-which(is.na(RMICEInormb2[,1])),]
RMICEInormg = RMICEInormg[-which(is.na(RMICEInormg[,1])),]}

RMICEHpmmb1 = read.table("RMICEHpmmb1.txt")
RMICEHpmmb2 = read.table("RMICEHpmmb2.txt")
RMICEHpmmg = read.table("RMICEHpmmg.txt")
if (length(which(is.na(RMICEHpmmb1[,1])))>0)
{RMICEHpmmb1 = RMICEHpmmb1[-which(is.na(RMICEHpmmb1[,1])),]
RMICEHpmmb2 = RMICEHpmmb2[-which(is.na(RMICEHpmmb2[,1])),]
RMICEHpmmg = RMICEHpmmg[-which(is.na(RMICEHpmmg[,1])),]}

RMICEHIpmmb1 = read.table("RMICEHIpmmb1.txt")
RMICEHIpmmb2 = read.table("RMICEHIpmmb2.txt")
RMICEHIpmmg = read.table("RMICEHIpmmg.txt")
if (length(which(is.na(RMICEHIpmmb1[,1])))>0)
{RMICEHIpmmb1 = RMICEHIpmmb1[-which(is.na(RMICEHIpmmb1[,1])),]
RMICEHIpmmb2 = RMICEHIpmmb2[-which(is.na(RMICEHIpmmb2[,1])),]
RMICEHIpmmg = RMICEHIpmmg[-which(is.na(RMICEHIpmmg[,1])),]}

RMICEHnormb1 = read.table("RMICEHnormb1.txt")
RMICEHnormb2 = read.table("RMICEHnormb2.txt")
RMICEHnormg = read.table("RMICEHnormg.txt")
if (length(which(is.na(RMICEHnormb1[,1])))>0)
{RMICEHnormb1 = RMICEHnormb1[-which(is.na(RMICEHnormb1[,1])),]
RMICEHnormb2 = RMICEHnormb2[-which(is.na(RMICEHnormb2[,1])),]
RMICEHnormg = RMICEHnormg[-which(is.na(RMICEHnormg[,1])),]}

RMICEHInormb1 = read.table("RMICEHInormb1.txt")
RMICEHInormb2 = read.table("RMICEHInormb2.txt")
RMICEHInormg = read.table("RMICEHInormg.txt")
if (length(which(is.na(RMICEHInormb1[,1])))>0)
{RMICEHInormb1 = RMICEHInormb1[-which(is.na(RMICEHInormb1[,1])),]
RMICEHInormb2 = RMICEHInormb2[-which(is.na(RMICEHInormb2[,1])),]
RMICEHInormg = RMICEHInormg[-which(is.na(RMICEHInormg[,1])),]}

case = matrix(-1,ncol = 6, nrow = 30)
colnames(case) = c("TSE","SSE","Percentage(%)","Coverage(%)","Width of 95% C.I.","converged")

case[1,] = c(mean(RFULLb1[,2]),sd(RFULLb1[,1]),(mean(RFULLb1[,1])-beta1)*100/beta1,
             length(which(beta1>RFULLb1[,3] & beta1<RFULLb1[,4]))/10,mean(RFULLb1[,4]-RFULLb1[,3]), dim(RFULLb1)[1])
case[2,] = c(mean(RCOMPb1[,2]),sd(RCOMPb1[,1]),(mean(RCOMPb1[,1])-beta1)*100/beta1,
              length(which(beta1>RCOMPb1[,3] & beta1<RCOMPb1[,4]))/10,mean(RCOMPb1[,4]-RCOMPb1[,3]), dim(RCOMPb1)[1])
case[3,] = c(mean(RMICEnormb1[,2]),sd(RMICEnormb1[,1]),(mean(RMICEnormb1[,1])-beta1)*100/beta1,
              length(which(beta1>RMICEnormb1[,3] & beta1<RMICEnormb1[,4]))/10,mean(RMICEnormb1[,4]-RMICEnormb1[,3]), dim(RMICEnormb1)[1])
case[4,] = c(mean(RMICEInormb1[,2]),sd(RMICEInormb1[,1]),(mean(RMICEInormb1[,1])-beta1)*100/beta1,
              length(which(beta1>RMICEInormb1[,3] & beta1<RMICEInormb1[,4]))/10,mean(RMICEInormb1[,4]-RMICEInormb1[,3]), dim(RMICEInormb1)[1])
case[5,] = c(mean(RMICEpmmb1[,2]),sd(RMICEpmmb1[,1]),(mean(RMICEpmmb1[,1])-beta1)*100/beta1,
              length(which(beta1>RMICEpmmb1[,3] & beta1<RMICEpmmb1[,4]))/10,mean(RMICEpmmb1[,4]-RMICEpmmb1[,3]), dim(RMICEpmmb1)[1])
case[6,] = c(mean(RMICEIpmmb1[,2]),sd(RMICEIpmmb1[,1]),(mean(RMICEIpmmb1[,1])-beta1)*100/beta1,
              length(which(beta1>RMICEIpmmb1[,3] & beta1<RMICEIpmmb1[,4]))/10,mean(RMICEIpmmb1[,4]-RMICEIpmmb1[,3]), dim(RMICEIpmmb1)[1])
case[7,] = c(mean(RMICEHnormb1[,2]),sd(RMICEHnormb1[,1]),(mean(RMICEHnormb1[,1])-beta1)*100/beta1,
              length(which(beta1>RMICEHnormb1[,3] & beta1<RMICEHnormb1[,4]))/10,mean(RMICEHnormb1[,4]-RMICEHnormb1[,3]), dim(RMICEHnormb1)[1])
case[8,] = c(mean(RMICEHInormb1[,2]),sd(RMICEHInormb1[,1]),(mean(RMICEHInormb1[,1])-beta1)*100/beta1,
              length(which(beta1>RMICEHInormb1[,3] & beta1<RMICEHInormb1[,4]))/10,mean(RMICEHInormb1[,4]-RMICEHInormb1[,3]), dim(RMICEHInormb1)[1])
case[9,] = c(mean(RMICEHpmmb1[,2]),sd(RMICEHpmmb1[,1]),(mean(RMICEHpmmb1[,1])-beta1)*100/beta1,
              length(which(beta1>RMICEHpmmb1[,3] & beta1<RMICEHpmmb1[,4]))/10,mean(RMICEHpmmb1[,4]-RMICEHpmmb1[,3]), dim(RMICEHpmmb1)[1])
case[10,] = c(mean(RMICEHIpmmb1[,2]),sd(RMICEHIpmmb1[,1]),(mean(RMICEHIpmmb1[,1])-beta1)*100/beta1,
              length(which(beta1>RMICEHIpmmb1[,3] & beta1<RMICEHIpmmb1[,4]))/10,mean(RMICEHIpmmb1[,4]-RMICEHIpmmb1[,3]), dim(RMICEHIpmmb1)[1])

case[11,] = c(mean(RFULLb2[,2]),sd(RFULLb2[,1]),(mean(RFULLb2[,1])-beta2)*100/beta2,
              length(which(beta2>RFULLb2[,3] & beta2<RFULLb2[,4]))/10,mean(RFULLb2[,4]-RFULLb2[,3]), dim(RFULLb2)[1])
case[12,] = c(mean(RCOMPb2[,2]),sd(RCOMPb2[,1]),(mean(RCOMPb2[,1])-beta2)*100/beta2,
              length(which(beta2>RCOMPb2[,3] & beta2<RCOMPb2[,4]))/10,mean(RCOMPb2[,4]-RCOMPb2[,3]), dim(RCOMPb2)[1])
case[13,] = c(mean(RMICEnormb2[,2]),sd(RMICEnormb2[,1]),(mean(RMICEnormb2[,1])-beta2)*100/beta2,
              length(which(beta2>RMICEnormb2[,3] & beta2<RMICEnormb2[,4]))/10,mean(RMICEnormb2[,4]-RMICEnormb2[,3]), dim(RMICEnormb2)[1])
case[14,] = c(mean(RMICEInormb2[,2]),sd(RMICEInormb2[,1]),(mean(RMICEInormb2[,1])-beta2)*100/beta2,
              length(which(beta2>RMICEInormb2[,3] & beta2<RMICEInormb2[,4]))/10,mean(RMICEInormb2[,4]-RMICEInormb2[,3]), dim(RMICEInormb2)[1])
case[15,] = c(mean(RMICEpmmb2[,2]),sd(RMICEpmmb2[,1]),(mean(RMICEpmmb2[,1])-beta2)*100/beta2,
              length(which(beta2>RMICEpmmb2[,3] & beta2<RMICEpmmb2[,4]))/10,mean(RMICEpmmb2[,4]-RMICEpmmb2[,3]), dim(RMICEpmmb2)[1])
case[16,] = c(mean(RMICEIpmmb2[,2]),sd(RMICEIpmmb2[,1]),(mean(RMICEIpmmb2[,1])-beta2)*100/beta2,
              length(which(beta2>RMICEIpmmb2[,3] & beta2<RMICEIpmmb2[,4]))/10,mean(RMICEIpmmb2[,4]-RMICEIpmmb2[,3]), dim(RMICEIpmmb2)[1])
case[17,] = c(mean(RMICEHnormb2[,2]),sd(RMICEHnormb2[,1]),(mean(RMICEHnormb2[,1])-beta2)*100/beta2,
              length(which(beta2>RMICEHnormb2[,3] & beta2<RMICEHnormb2[,4]))/10,mean(RMICEHnormb2[,4]-RMICEHnormb2[,3]), dim(RMICEHnormb2)[1])
case[18,] = c(mean(RMICEHInormb2[,2]),sd(RMICEHInormb2[,1]),(mean(RMICEHInormb2[,1])-beta2)*100/beta2,
              length(which(beta2>RMICEHInormb2[,3] & beta2<RMICEHInormb2[,4]))/10,mean(RMICEHInormb2[,4]-RMICEHInormb2[,3]), dim(RMICEHInormb2)[1])
case[19,] = c(mean(RMICEHpmmb2[,2]),sd(RMICEHpmmb2[,1]),(mean(RMICEHpmmb2[,1])-beta2)*100/beta2,
              length(which(beta2>RMICEHpmmb2[,3] & beta2<RMICEHpmmb2[,4]))/10,mean(RMICEHpmmb2[,4]-RMICEHpmmb2[,3]), dim(RMICEHpmmb2)[1])
case[20,] = c(mean(RMICEHIpmmb2[,2]),sd(RMICEHIpmmb2[,1]),(mean(RMICEHIpmmb2[,1])-beta2)*100/beta2,
              length(which(beta2>RMICEHIpmmb2[,3] & beta2<RMICEHIpmmb2[,4]))/10,mean(RMICEHIpmmb2[,4]-RMICEHIpmmb2[,3]), dim(RMICEHIpmmb2)[1])

case[21,] = c(mean(RFULLg[,2]),sd(RFULLg[,1]),(mean(RFULLg[,1])-gamma0)*100/gamma0,
              length(which(gamma0>RFULLg[,3] & gamma0<RFULLg[,4]))/10,mean(RFULLg[,4]-RFULLg[,3]), dim(RFULLg)[1])
case[22,] = c(mean(RCOMPg[,2]),sd(RCOMPg[,1]),(mean(RCOMPg[,1])-gamma0)*100/gamma0,
              length(which(gamma0>RCOMPg[,3] & gamma0<RCOMPg[,4]))/10,mean(RCOMPg[,4]-RCOMPg[,3]), dim(RCOMPg)[1])
case[23,] = c(mean(RMICEnormg[,2]),sd(RMICEnormg[,1]),(mean(RMICEnormg[,1])-gamma0)*100/gamma0,
              length(which(gamma0>RMICEnormg[,3] & gamma0<RMICEnormg[,4]))/10,mean(RMICEnormg[,4]-RMICEnormg[,3]), dim(RMICEnormg)[1])
case[24,] = c(mean(RMICEInormg[,2]),sd(RMICEInormg[,1]),(mean(RMICEInormg[,1])-gamma0)*100/gamma0,
              length(which(gamma0>RMICEInormg[,3] & gamma0<RMICEInormg[,4]))/10,mean(RMICEInormg[,4]-RMICEInormg[,3]), dim(RMICEInormg)[1])
case[25,] = c(mean(RMICEpmmg[,2]),sd(RMICEpmmg[,1]),(mean(RMICEpmmg[,1])-gamma0)*100/gamma0,
              length(which(gamma0>RMICEpmmg[,3] & gamma0<RMICEpmmg[,4]))/10,mean(RMICEpmmg[,4]-RMICEpmmg[,3]), dim(RMICEpmmg)[1])
case[26,] = c(mean(RMICEIpmmg[,2]),sd(RMICEIpmmg[,1]),(mean(RMICEIpmmg[,1])-gamma0)*100/gamma0,
              length(which(gamma0>RMICEIpmmg[,3] & gamma0<RMICEIpmmg[,4]))/10,mean(RMICEIpmmg[,4]-RMICEIpmmg[,3]), dim(RMICEIpmmg)[1])
case[27,] = c(mean(RMICEHnormg[,2]),sd(RMICEHnormg[,1]),(mean(RMICEHnormg[,1])-gamma0)*100/gamma0,
              length(which(gamma0>RMICEHnormg[,3] & gamma0<RMICEHnormg[,4]))/10,mean(RMICEHnormg[,4]-RMICEHnormg[,3]), dim(RMICEHnormg)[1])
case[28,] = c(mean(RMICEHInormg[,2]),sd(RMICEHInormg[,1]),(mean(RMICEHInormg[,1])-gamma0)*100/gamma0,
              length(which(gamma0>RMICEHInormg[,3] & gamma0<RMICEHInormg[,4]))/10,mean(RMICEHInormg[,4]-RMICEHInormg[,3]), dim(RMICEHInormg)[1])
case[29,] = c(mean(RMICEHpmmg[,2]),sd(RMICEHpmmg[,1]),(mean(RMICEHpmmg[,1])-gamma0)*100/gamma0,
              length(which(gamma0>RMICEHpmmg[,3] & gamma0<RMICEHpmmg[,4]))/10,mean(RMICEHpmmg[,4]-RMICEHpmmg[,3]), dim(RMICEHpmmg)[1])
case[30,] = c(mean(RMICEHIpmmg[,2]),sd(RMICEHIpmmg[,1]),(mean(RMICEHIpmmg[,1])-gamma0)*100/gamma0,
              length(which(gamma0>RMICEHIpmmg[,3] & gamma0<RMICEHIpmmg[,4]))/10,mean(RMICEHIpmmg[,4]-RMICEHIpmmg[,3]), dim(RMICEHIpmmg)[1])

write.csv(case,"summary_miaft.csv")

# results are written into the file "summary_miaft.csv"
# rows 1-6, 11-16, 21-26 correspond to the results in Table 8.9
# need to remove all the files written out to the home directory after the simulation








