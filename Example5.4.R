# Example 5.4

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
mu_x=5;
var_x=1;

# number of simulations;
cycle_no=1000;

# number of multiple imputations;
mi_no=30;


# linear normal model;
# regression parameters;
# Scenario I 
# beta0=10;
# Scenario II
beta0=0 
beta1=1;
beta2=1;

# error variance;
var_error=1;

# generate the population;
y_complete_mat=matrix(NA, rowobs, cycle_no);

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x1=mu_x+sqrt(var_x)*rnorm(rowobs);
x2=mu_x+sqrt(var_x)*rnorm(rowobs);


# normal error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the y variable via linear normal model;
# Scenario I
# y=beta0+beta1*x1+beta2*x2+error;

# generate the y variable via transformation;
# Scenario II
y=(beta0+beta1*x1+beta2*x2+error)^4;
y_complete_mat[,cycle]=y;

}

# get the percentiles;
# quantile(y, c(0.05, 0.25, 0.50, 0.75, 0.95));
y_cutoff=quantile(y_complete_mat, 0.05);

# Vectors and matrices holding the parameter estimates;
# Before-deletion analysis;
BD_mean_vec=BD_mean_var_vec=BD_prop_vec=BD_prop_var_vec=rep(NA, cycle_no);

# Complete-case analysis;
CC_mean_vec=CC_mean_var_vec=CC_prop_vec=CC_prop_var_vec=rep(NA, cycle_no);

# Multiple imputation;
# normal model imputation;
MI_mean_mat=MI_mean_var_mat=MI_prop_mat=MI_prop_var_mat=matrix(NA, cycle_no, mi_no);

# predictive mean matching (PMM) imputation
IM_MI_mean_mat=IM_MI_mean_var_mat=IM_MI_prop_mat=IM_MI_prop_var_mat=matrix(NA, cycle_no, mi_no);


# begin the simulation;

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x1=mu_x+sqrt(var_x)*rnorm(rowobs);
x2=mu_x+sqrt(var_x)*rnorm(rowobs);


# normal error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the y variable via linear normal model;
# Scenario I
# y=beta0+beta1*x1+beta2*x2+error;

# generate the y variable via transformation;
# Scenario II
y=(beta0+beta1*x1+beta2*x2+error)^4;

# set up the missing data as MCAR on x;
miss_indi=runif(n=rowobs)<0.30;

# set up the missing data as MAR on x;
# alpha0=-1.4;
# alpha1=1;
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x)/(1+exp(alpha0+alpha1*x)));


y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
# x_obs=x[miss_indi==0];
# x_mis=x[miss_indi==1];

# compelete data analysis;
# marginal mean;
BD_mean_vec[cycle]=mean(y);
BD_mean_var_vec[cycle]=var(y)/rowobs;


# proportion analysis;
BD_prop_vec[cycle]=mean(y<y_cutoff);
BD_prop_var_vec[cycle]=var(y<y_cutoff)/rowobs;

# different missing data methods;

# complete-case analysis;

# marginal mean;
CC_mean_vec[cycle]=mean(y_obs);
CC_mean_var_vec[cycle]=var(y_obs)/obs_no;


# proportion analysis;
CC_prop_vec[cycle]=mean(y_obs<y_cutoff);
CC_prop_var_vec[cycle]=var(y_obs<y_cutoff)/obs_no;


# now impute the missing y's to get estimates of the marginal;
y_completed=y_completed_IM=y_miss;

for (i in 1:mi_no)
{
set.seed(i);
# normal model imputation;

y_imputed=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x1,x2));


y_completed[miss_seq]=y_imputed;

# completed-data analysis;

# marginal mean;
MI_mean_mat[cycle,i]=mean(y_completed);
MI_mean_var_mat[cycle,i]=var(y_completed)/rowobs;

# proportion analysis;

MI_prop_mat[cycle,i]=mean(y_completed<y_cutoff);
MI_prop_var_mat[cycle,i]=var(y_completed<y_cutoff)/rowobs;

# PMM imputation;
y_imputed_IM=mice.impute.pmm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x1,x2));
y_completed_IM[miss_seq]=y_imputed_IM;

# completed-data analysis;

# marginal mean;
IM_MI_mean_mat[cycle,i]=mean(y_completed_IM);
IM_MI_mean_var_mat[cycle,i]=var(y_completed_IM)/rowobs;

# proportion analysis;

IM_MI_prop_mat[cycle,i]=mean(y_completed_IM < y_cutoff);
IM_MI_prop_var_mat[cycle,i]=var(y_completed_IM < y_cutoff)/rowobs;



}

cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

# plot the distribution of the data;

# check the performance of the estimates;
# population quantity

#############################################################################
# Simulation results
# Table 5.3
# marginal mean estimand of y;

true_mean=mean(BD_mean_vec);
true_mean;

# before-deletion;
BD=performance(BD_mean_vec, BD_mean_var_vec, true_mean);
BD;

# complete-case analysis;
CC=performance(CC_mean_vec, CC_mean_var_vec, true_mean);
CC;

# Multiple imputation analysis;
# normal imputation;
MI=mi_performance(MI_mean_mat, MI_mean_var_mat, true_mean, rowobs);
MI;

# PMM imputation;
IM_MI=mi_performance(IM_MI_mean_mat, IM_MI_mean_var_mat, true_mean, rowobs);
IM_MI;


##########################################################################
# proportional estimand;

# complete-data analysis;

true_prop=mean(BD_prop_vec);
true_prop;

# before-deletion;
BD=performance(BD_prop_vec, BD_prop_var_vec, true_prop);
BD;

# complete-case analysis;
CC=performance(CC_prop_vec, CC_prop_var_vec, true_prop);
CC;

# Multiple imputation analysis;
# normal imputation 
MI=mi_performance(MI_prop_mat, MI_prop_var_mat, true_prop, rowobs);
MI;

# PMM imputation
IM_MI=mi_performance(IM_MI_prop_mat, IM_MI_prop_var_mat, true_prop, rowobs);
IM_MI;

################################################
# make some plots;
# using data from the last simulation
# Fig. 5.7
hist(y, freq=T, main="Before-deletion", xlab="y");
hist(y_miss, freq=T, main="Observed", xlab="y");
hist(y_completed, main="Normal Imputation", xlab="y");
hist(y_completed_IM, main="PMM Imputation", xlab="y");


# Fig. 5.8

y_plot_obs=y;
y_plot_obs[miss_indi==1]=NA;

x=x1+x2;

# for normal imputation;
y_plot_imputed=y_completed;
y_plot_imputed[miss_indi==0]=NA;
y_plot_impute=cbind(y_plot_obs, y_plot_imputed);

matplot(x[1:250], y_plot_impute[1:250,], pch=c(1,17), xlab="x1+x2", ylab="y", main="Normal Imputation");
abline(h=y_cutoff, col="green");

# for PMM imputation;

y_plot_imputed=y_completed_IM;
y_plot_imputed[miss_indi==0]=NA;
y_plot_impute=cbind(y_plot_obs, y_plot_imputed);

matplot(x[1:250], y_plot_impute[1:250,], pch=c(1,17), xlab="x1+x2", ylab="y", main="PMM Imputation");
abline(h=y_cutoff, col="green");



