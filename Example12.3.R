# Example 12.3

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
beta2=1;

# error variance;
var_error=1;

cycle=1;

set.seed(cycle);

x1=mu_x+sqrt(var_x)*rnorm(rowobs);

x2=mu_x+sqrt(var_x)*rnorm(rowobs);

# normal error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the y variable;
y=beta0+beta1*x1+beta2*x2+error;

# set up the missing data as MAR on x1;
alpha0=-1.4;
alpha1=1;
miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x1)/(1+exp(alpha0+alpha1*x1)));


y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];

x1_obs=x1[miss_indi==0];
x1_mis=x1[miss_indi==1];

x2_obs=x2[miss_indi==0];
x2_mis=x2[miss_indi==1];


# only one imputation;

mi_no=1;

y_completed_x1=y_completed_x2=y_completed_x1_x2=y_miss;


for (i in 1:mi_no)
{
# imputation using only x1;
y_imputed_x1=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x1);


y_completed_x1[miss_seq]=y_imputed_x1;

# imputation using only x2;
y_imputed_x2=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x2);

y_completed_x2[miss_seq]=y_imputed_x2;

# imputation using both x1 and x2;
y_imputed_x1_x2=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x1,x2));

y_completed_x1_x2[miss_seq]=y_imputed_x1_x2;


}


# the end of imputation;

# plots;
# Fig. 12.3;

# First row;
plot(x1_obs, y_obs);
abline(lm(y_obs~x1_obs));
lm(y_obs ~ x1_obs);

plot(x2_obs, y_obs);
abline(lm(y_obs~x2_obs));
lm(y_obs ~ x2_obs);

# Second row;

plot(x1_mis, y_imputed_x1, pch=17, ylab="y_imp: MI-X1", col="red");
abline(lm(y_imputed_x1 ~ x1_mis), col="red");
lm(y_imputed_x1 ~ x1_mis);

plot(x2_mis, y_imputed_x1, pch=17, ylab="y_imp: MI-X1", col="red");
abline(lm(y_imputed_x1 ~ x2_mis), col="red");
lm(y_imputed_x1 ~ x2_mis);

# Third row;
plot(x1_mis, y_imputed_x2, pch=17, ylab="y_imp: MI-X2", col="red");
abline(lm(y_imputed_x2 ~ x1_mis), col="red");
lm(y_imputed_x2 ~ x1_mis);

plot(x2_mis, y_imputed_x2, pch=17, ylab="y_imp: MI-X2", col="red");
abline(lm(y_imputed_x2 ~ x2_mis), col="red");
lm(y_imputed_x2 ~ x2_mis);

# Forth row;
plot(x1_mis, y_imputed_x1_x2, pch=17, ylab="y_imp: MI-X1X2", col="red");
abline(lm(y_imputed_x1_x2 ~ x1_mis), col="red");
lm(y_imputed_x1_x2 ~ x1_mis);

plot(x2_mis, y_imputed_x1_x2, pch=17, ylab="y_imp: MI-X1X2", col="red");
abline(lm(y_imputed_x1_x2 ~ x2_mis), col="red");
lm(y_imputed_x1_x2 ~ x2_mis);


