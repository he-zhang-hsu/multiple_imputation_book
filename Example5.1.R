# Example 5.1

rm(list=ls());
library(car);
library(mvtnorm);
library(mice);
library(norm);
library(HI);
options(digits=4);

# function performance() is used to evaluate the simulation estimates;
# It will output bias, relative bias, MSE, standard deviation, 
# mean of standard errors, coverage rate, and lengths of 95% confidence intervals
# the gold standard is the mean of estimates from before-deletion method;


performance=function(mean_vec, var_vec, true)
{
bias=mean(mean_vec)-true;
rbias=bias/true;
mse=mean((mean_vec-true)^2);
# lower and upper 95% CI;
low95_vec=mean_vec-1.96*sqrt(var_vec);
up95_vec=mean_vec+1.96*sqrt(var_vec);
mean_length=mean(up95_vec-low95_vec);
mean_cov=mean((low95_vec < true)*(up95_vec > true));
sd=sqrt(var(mean_vec));
se=mean(sqrt(var_vec));
estimand=as.data.frame(cbind(bias,rbias,mse,mean_length,mean_cov,sd,se))
}

# function mi_performance() is used to evaluate the simulation estimates for 
# multiple imputation methods;
# It will output bias, relative bias, MSE, standard deviation, 
# mean of standard errors, coverage rate, and lengths of 95% confidence intervals
# the gold standard is the mean of estimates from before-deletion method;

mi_performance=function(mean_mat, var_mat, true, n)
{
cycle_no=nrow(mean_mat);
mi_mean=rep(NA, cycle_no);
mi_var=rep(NA, cycle_no);
mi_df=rep(NA, cycle_no);
mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(mean_mat[i,], var_mat[i,], n=n);
mi_mean[i]=summary$qbar;
mi_var[i]=summary$t;
mi_df[i]=summary$df;
mi_f[i]=summary$f;
}

bias=mean(mi_mean)-true;
rbias=bias/true;
mse=mean((mi_mean-true)^2);
# coverage;
low95_vec=mi_mean-qt(.975, mi_df)*sqrt(mi_var);
up95_vec=mi_mean+qt(.975, mi_df)*sqrt(mi_var);
mean_length=mean(up95_vec-low95_vec);
mean_cov=mean((low95_vec < true)*(up95_vec > true));
se=mean(sqrt(mi_var));
sd=sqrt(var(mi_mean));
frac=mean(mi_f);
estimand=as.data.frame(cbind(bias,rbias,mse,mean_length,mean_cov,sd,se,frac))
}


# set up the random seed;
set.seed(197789);

# complete-data sample size in the simulation;

rowobs=1000;

# mean and variance of X;
mu_x=10;
var_x=4;

# regression parameters; 
beta0=1;
beta1=1;

# transformation parameter;
lambda=1/3;


# error variance;
# r_square is the variance explained by the model;
r_square=0.5;
var_error=(1-r_square)/r_square * beta1^2*var_x;


# number of simulations;
cycle_no=1000;

# number of multiple imputations;
mi_no=30;

# vectors and matrices holding the parameter estimates;
# marginal mean;
# Before-deletion analysis
BD_mean_vec=BD_var_vec=rep(NA, cycle_no);

# Complete-case analysis;
CC_mean_vec=CC_var_vec=rep(NA, cycle_no);

# norml model imputation
norm_mean_mat=norm_var_mat=matrix(NA, cycle_no, mi_no);

# imputation based on transformation;
tran_mean_mat=tran_var_mat=matrix(NA, cycle_no, mi_no);

# less than 5%-tile
# Before-deletion analysis;
BD_prop_bottom_vec=rep(NA, cycle_no);
BD_prop_bottom_var_vec=rep(NA, cycle_no);

# Complete-case analysis;
CC_prop_bottom_vec=rep(NA, cycle_no);
CC_prop_bottom_var_vec=rep(NA, cycle_no);


# normal imputation;
norm_prop_bottom_mat=matrix(NA, cycle_no, mi_no);
norm_prop_bottom_var_mat=matrix(NA, cycle_no, mi_no);


# imputation based on transformation;
tran_prop_bottom_mat=matrix(NA, cycle_no, mi_no);
tran_prop_bottom_var_mat=matrix(NA, cycle_no, mi_no);



# Vector holding estimated lambdas;
lambda_observed_vec=var_lambda_vec=rep(NA, cycle_no);

# Vector holding the number of observed cases;
obsno_vec=rep(NA, cycle_no);

# generate the complete data;
# to estimate the populatoin quantiles;
y_complete_mat=matrix(NA, rowobs, cycle_no);

for (cycle in 1:cycle_no)
{
 set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# normal error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the y variable;
# take the absolute value to ensure y1 is positive;

y1=abs(beta0+beta1*x+error);

# transform y1 using box-cox transformation;
y=(1+lambda*y1)^(1/lambda);
y_complete_mat[,cycle]=y;

}

population_quantile=quantile(y_complete_mat, c(0.05, 0.25, 0.5, 0.75, 0.95));

# begin the simulation;

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# normal error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the y variable;
# take the absolute value to ensure y1 is positive;

y1=abs(beta0+beta1*x+error);

# transform y1 using box-cox transformation;
y=(1+lambda*y1)^(1/lambda);

# set up the missing data as MAR on x;
# alpha0=-17; # 30% missing;
# alpha0=-16; # 50% missing;
# alpha0=-40; # 50% missing;
alpha0=36 # 30% missing;
alpha1=-4;
miss_indi=rbinom(n=rowobs,size=1,p=1/(1+exp(-(alpha0+alpha1*x))));
# MCAR
# miss_indi=runif(n=rowobs)<0.30;
mean(miss_indi);

miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];

y_miss=y;
y_miss[miss_indi==1]=NA;

x_miss=x;
x_miss[miss_indi==1]=NA;
x_obs=x[miss_indi==0];

miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(x_miss[!is.na(x_miss)]);
mis_no=rowobs-obs_no;

y_obs=y[miss_indi==0];

# observed data is cbind(x, y_miss);

# Before-deletion analysis;
# sample mean;
y_mean=mean(y);
y_var=var(y)/rowobs;
BD_mean_vec[cycle]=y_mean;
BD_var_vec[cycle]=y_var;


# proportions less than the bottom 5%-tile of the population;
BD_prop_bottom_vec[cycle]=mean(y<population_quantile[1]);
BD_prop_bottom_var_vec[cycle]=BD_prop_bottom_vec[cycle]*(1-BD_prop_bottom_vec[cycle])/rowobs;

# different missing data methods;

# Complete-case analysis;
# mean estimands;
CC_mean=mean(y_miss, na.rm="T");
CC_mean_vec[cycle]=CC_mean;
CC_var_vec[cycle]=var(y_miss, na.rm="T")/length(y_miss[!is.na(y_miss)]);

# proportions less than the 5%-tile of the population;
CC_prop_bottom_vec[cycle]=mean(y_miss<population_quantile[1], na.rm=T);
CC_prop_bottom_var_vec[cycle]=CC_prop_bottom_vec[cycle]*(1-CC_prop_bottom_vec[cycle])/obs_no;


# obsno_vec contains the number of observed cases from each simualtion;
obsno_vec[cycle]=obs_no;

# estimate transformation parameter conditional on x;
# observed data analysis;
# first obtain the transformation parameter using observed data;
lambda_observed=powerTransform(y_miss~x);

# point estimate of lambda;
lambda_observed_vec[cycle]=lambda_observed$lambda;

# variance estimate of lambda;
var_lambda=1/lambda_observed$hessian;
var_lambda_vec[cycle]=var_lambda;

# then obtain the transformed y;
tran_y_observed=bcPower(y_miss, lambda_observed$lambda);

# now impute the missing y's to get estimates of the marginal;
# normal imputation;
y_norm_completed=y_miss;

 for (i in 1:mi_no)
{
set.seed(i);
# simple linear imputation on the original scale;
y_norm_imputed=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);

y_norm_completed[miss_seq]=y_norm_imputed;

# set the negative imputed values as zero.
y_norm_completed[y_norm_completed <0]=0;

# marginal means;
norm_mean=mean(y_norm_completed);
norm_var=var(y_norm_completed)/rowobs;

norm_mean_mat[cycle,i]=norm_mean;
norm_var_mat[cycle,i]=norm_var;

# proportions;
norm_prop_bottom_mat[cycle,i]=mean(y_norm_completed<population_quantile[1]);
norm_prop_bottom_var_mat[cycle,i]=norm_prop_bottom_mat[cycle,i]*(1-norm_prop_bottom_mat[cycle,i])/rowobs;

}


# now impute the missing y's using transformation;
y_tran_completed=y_miss;

for (i in 1:mi_no)
{
set.seed(i);
# draw lambda using normal approximation;
# ensure the draw lambad to be positive
lambda_tran=abs(lambda_observed$lambda+sqrt(var_lambda)*rnorm(1));

# imputation on the transformed scale;
tran_y_tran=bcPower(y_miss, lambda_tran);

tran_y_tran_imputed=mice.impute.norm(tran_y_tran, ry=as.logical(1-miss_indi), seed=i, x=x);

# transform back
y_tran_completed[miss_seq]=abs((1+lambda_tran*as.vector(tran_y_tran_imputed)))^(1/lambda_tran);

# marginal means;

tran_mean=mean(y_tran_completed);
tran_var=var(y_tran_completed)/rowobs;

tran_mean_mat[cycle,i]=tran_mean;
tran_var_mat[cycle,i]=tran_var;

# proportions;
tran_prop_bottom_mat[cycle,i]=mean(y_tran_completed<population_quantile[1]);
tran_prop_bottom_var_mat[cycle,i]=tran_prop_bottom_mat[cycle,i]*(1-tran_prop_bottom_mat[cycle,i])/rowobs;

}

cat("the cycle is", cycle, "\n");

}

# the end of multiple imputation;

# plot the distribution of the data from the last simulation;
# Fig. 5.1 histograms;
hist(y, xlab="y", main="Before-deletion");
hist(y_miss, xlab="y", main="Observed");
hist(y_norm_completed, xlab="y", main="Normal Imputation");
hist(y_tran_completed, xlab="y", main="Transformation Imputation");

# Fig. 5.2 scatter plot.
test=lm(y_obs~x_obs);
test2=lowess(y_obs ~x_obs);
j <- order(x_obs)
plot(x_obs,y_obs, xlab="x", lty=1, ylab="y");
# plot(test2$x, test2$y);
lines(test2$x,test2$y,col="red",lwd=3)
abline(test, col="green", lty=2, lwd=3);

# check the performance of the estimates;
# population quantity
# mean estimand;

# mean(lambda_observed_vec);
# hist(lambda_observed_vec);
# mean(var_lambda_vec);

#############################################################################
# Simulation performance
# Mean estimates
#
# For Table 5.1
true_mean=mean(BD_mean_vec);
true_mean;

# before-deletion;
BD=performance(BD_mean_vec, BD_var_vec, true_mean);
BD;

# complete-case analysis;
CC=performance(CC_mean_vec, CC_var_vec, true_mean);
CC;

# Multiple imputation analysis;
# normal imputation 

norm_MI=mi_performance(norm_mean_mat, norm_var_mat, true_mean, rowobs);
norm_MI;

# Imputation based on transformation
tran_MI=mi_performance(tran_mean_mat, tran_var_mat, true_mean, rowobs);
tran_MI;


###########################################################
# the proportion less than 5%-tile
# For Table 5.1

true_mean=mean(BD_prop_bottom_vec);
true_mean;


# before-deletion;
BD=performance(BD_prop_bottom_vec, BD_prop_bottom_var_vec, true_mean);
BD;

# complete-case analysis;
CC=performance(CC_prop_bottom_vec, CC_prop_bottom_var_vec, true_mean);
CC;

# Multiple imputation analysis;
# normal imputation 
norm_MI=mi_performance(norm_prop_bottom_mat, norm_prop_bottom_var_mat, true_mean, rowobs);
norm_MI;

# Imputation based on transformation
tran_MI=mi_performance(tran_prop_bottom_mat, tran_prop_bottom_var_mat, true_mean, rowobs);
tran_MI;
