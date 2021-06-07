# Example 12.1

rm(list=ls());
# library(car);
# library(mvtnorm);
library(mice);
library(norm);
# library(HI);
options(digits=4);


# set up the random seed;
set.seed(197789);

# complete-data sample size;

rowobs=1000;

# mean and variance of x;
mu_x=1;
var_x=1;

# regression parameters; 
beta0=-2;
beta1=1;

# error variance;
var_error=1;

# only one simulation;
cycle_no=1;

# number of multiple imputations;
mi_no=50;

# begin the multiple imputation;

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# normal error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the y variable;
y=beta0+beta1*x+error;

# set up the missing data as MCAR on x;
# miss_indi=runif(n=rowobs)<0.40;

# set up the missing data as MAR on x;
alpha0=-1.4;
alpha1=1;
miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x)/(1+exp(alpha0+alpha1*x)));


y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
x_obs=x[miss_indi==0];
x_mis=x[miss_indi==1];


# now impute the missing y to get estimates of the marginal;

y_completed=y_completed_IM=y_miss;

for (i in 1:mi_no)
{
# linear model imputation;

y_imputed=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);

y_completed[miss_seq]=y_imputed;


}

cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

# plot the distribution of the data;
# using the last imputation;

# MCAR ;
# Fig. 12.1 top row;
# hist(y_obs, freq=FALSE, main="MCAR");
# hist(y_imputed, freq=FALSE, xlab="y_imp", main="MCAR");

# Table 12.1
# mean(y_obs);
# var(y_obs);
# mean(y_imputed);
# var(y_imputed);

# t.test(y_obs, y_imputed);
# ks.test(y_obs, y_imputed);

# Fig. 12.2 top row;
# plot(x_obs, y_obs, xlim=c(-2,5), main="MCAR");
# abline(lm(y_obs~x_obs));
# lm(y_obs ~ x_obs);

# plot(x_mis, y_imputed, xlim=c(-2,5), pch=17, col="red", ylab="y_imp", main="MCAR");
# abline(lm(y_imputed ~ x_mis), col="red");
# lm(y_imputed ~ x_mis);




# MAR;
# Fig. 12.1 bottom row;
hist(y_obs, freq=FALSE, main="MAR");
hist(y_imputed, freq=FALSE, xlab="y_imp", main="MAR");

# Table 12.1
mean(y_obs);
var(y_obs);
mean(y_imputed);
var(y_imputed);

# t.test(y_obs, y_imputed);
ks.test(y_obs, y_imputed);

# Fig. 12.2 bottom row;
plot(x_obs, y_obs, xlim=c(-2,5), main="MAR");
abline(lm(y_obs~x_obs));
lm(y_obs ~ x_obs);

plot(x_mis, y_imputed, xlim=c(-2,5), pch=17, col="red", ylab="y_imp", main="MAR");
abline(lm(y_imputed ~ x_mis), col="red");
lm(y_imputed ~ x_mis);





