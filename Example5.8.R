# Example 5.8

rm(list=ls());
# library(car);
# library(mvtnorm);
library(mice);
library(norm);
# library(HI);
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

# complete-data sample size;
rowobs=1000;

# mean and variance of x;
mu_x=1;
var_x=1;

# auxiliary variable z;
# mean and variance of z;
mu_z=1;
var_z=3;

# regression parameters; 
beta0=-2;
beta1=1;

# error variance;
var_error=1;

# the number of simulations;
cycle_no=1000;

# the number of multiple imputations;
mi_no=50;

# Vectors and matrices holding the parameter estimates;

BD_mean_vec=BD_mean_var_vec=BD_slope_vec=BD_slope_var_vec=rep(NA, cycle_no);
CC_mean_vec=CC_mean_var_vec=CC_slope_vec=CC_slope_var_vec=rep(NA, cycle_no);

# Exclusive imputation;
MI_mean_mat=MI_mean_var_mat=MI_slope_mat=MI_slope_var_mat=matrix(NA, cycle_no, mi_no);

# Inclusive imputation;
IM_mean_mat=IM_mean_var_mat=IM_slope_mat=IM_slope_var_mat=matrix(NA, cycle_no, mi_no);



# begin the simulation;

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# normal error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the outcome variable;
y=beta0+beta1*x+error;


# generaite auxiliary variable z;
# Scenario I
# not related to y nor related to the missingness mechanism;
# z=mu_z+sqrt(var_z)*rnorm(rowobs);
# set up the missing data as MCAR;
# miss_indi=runif(n=rowobs)<0.40;


# Scenario II
# related to Y but not related to the missingness mechanism;
# error_z=sqrt(var_error)*rnorm(rowobs);
# z=2+y+error_z;
# set up the missing data as MCAR;
# miss_indi=runif(n=rowobs)<0.40;



# Scenario III
# not related to Y but related to missingness mechanism;
# z=mu_z+sqrt(var_z)*rnorm(rowobs);
# alpha0=-1.7;
# alpha1=1
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*z)/(1+exp(alpha0+alpha1*z)));

# Scenario IV
#  related to Y and also related to missingness mechainsim;
error_z=sqrt(var_error)*rnorm(rowobs);
z=2+y+error_z;
 alpha0=-1.7;
 alpha1=1
 miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*z)/(1+exp(alpha0+alpha1*z)));

obs_indi=1-miss_indi;

y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
x_obs=x[miss_indi==0];
x_mis=x[miss_indi==1];

# compelete data analysis;
# marginal mean;
BD_mean_vec[cycle]=mean(y);
BD_mean_var_vec[cycle]=var(y)/rowobs;


# regression analysis;
BD=summary(lm(y~x));

BD_slope_vec[cycle]=BD$coeff[2,1];
BD_slope_var_vec[cycle]=BD$coeff[2,2]^2;


# different missing data methods;

# complete-case analysis;

# marginal mean;
CC_mean_vec[cycle]=mean(y_obs);
CC_mean_var_vec[cycle]=var(y_obs)/obs_no;


# regression analysis;
CC=summary(lm(y_obs~x_obs));

CC_slope_vec[cycle]=CC$coeff[2,1];
CC_slope_var_vec[cycle]=CC$coeff[2,2]^2;


# multiple imputation analysis;

# now impute the missing y's to get estimates of the marginal;
y_completed_MI=y_completed_IM=y_miss;

# the matrix holds the multiply imputed data;

y_completed_MI_mat=y_completed_IM_mat=matrix(NA, nrow=rowobs, ncol=mi_no);

for (i in 1:mi_no)
{
# Exclusive imputation;

set.seed(i);

y_imputed_MI=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);


y_completed_MI[miss_seq]=y_imputed_MI;

y_completed_MI_mat[,i]=y_completed_MI;

# completed-data analysis;

# marginal mean;
MI_mean_mat[cycle,i]=mean(y_completed_MI);
MI_mean_var_mat[cycle,i]=var(y_completed_MI)/rowobs;

# regression analysis;

IMP=summary(lm(y_completed_MI~x));

MI_slope_mat[cycle,i]=IMP$coeff[2,1];
MI_slope_var_mat[cycle,i]=IMP$coeff[2,2]^2;


# exclusive imputation;

y_imputed_IM=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x,z));

y_completed_IM[miss_seq]=y_imputed_IM;

y_completed_IM_mat[,i]=y_completed_IM;

# completed-data analysis;

# marginal mean;
IM_mean_mat[cycle,i]=mean(y_completed_IM);
IM_mean_var_mat[cycle,i]=var(y_completed_IM)/rowobs;

# regression analysis;

IMP_IM=summary(lm(y_completed_IM~x));

IM_slope_mat[cycle,i]=IMP_IM$coeff[2,1];
IM_slope_var_mat[cycle,i]=IMP_IM$coeff[2,2]^2;


}


cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;


#############################################################################
# Simulation results
# Table 5.6

# marginal mean of Y;
true_mean=mean(BD_mean_vec);
true_mean;

# before-deletion;
BD=performance(BD_mean_vec, BD_mean_var_vec, true_mean);
BD;

# complete-case analysis;
CC=performance(CC_mean_vec, CC_mean_var_vec, true_mean);
CC;

# Multiple imputation analysis;
# Exclusive imputation;
MI=mi_performance(MI_mean_mat, MI_mean_var_mat, true_mean, rowobs);
MI;

# Inclusive imputation;
IM=mi_performance(IM_mean_mat, IM_mean_var_mat, true_mean, rowobs);
IM;


##########################################################################
# slope estimand;
# 

true_slope=mean(BD_slope_vec);
true_slope;

# before-deletion;
BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;

# complete-case analysis;
CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;

# Exclusive multiple imputation;
MI=mi_performance(MI_slope_mat, MI_slope_var_mat, true_slope, rowobs);
MI;

# Inclusive multiple imputation;
IM=mi_performance(IM_slope_mat, IM_slope_var_mat, true_slope, rowobs);
IM;


