# Example 14.3

rm(list=ls());
# library(car);
# library(mvtnorm);
library(mice);
library(norm);
library(boot);
library(Hmisc);
# library(HI);
options(digits=4);


# rsq function to obtain R-Squared from the data 
rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
} 

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
estimand=as.data.frame(cbind(bias,rbias,mse,mean_length,mean_cov,sd,se))
}


# set up the random seed;
set.seed(197789);

# complete-data sample size;
rowobs=100;

# mean and variance of x;
mu_x=1;
var_x=1;

# error variance;
var_error=1;

# number of simulations;
cycle_no=1000;

# number of multiple imputations;
mi_no=50;

# number of bootstrap samples;
boot_no=200;

# matrices holding the parameter estimates;
# R-square from the linear regression;
BD_mean_vec=BD_var_vec=CC_mean_vec=CC_var_vec=rep(NA, cycle_no);

# delta means delta method;
# trans means transformation method;
# bootcombine means bootstrap method;

delta_mean_mat=delta_var_mat=trans_mean_mat=trans_var_mat=bootcombine_mean_mat=bootcombine_var_mat=matrix(NA, cycle_no, mi_no);

boot_CI_matrix=matrix(NA, cycle_no,2);

# regression parameters; 
beta0=-2;
beta1_vec=c(1/sqrt(19), 1/3, 0.5, 1, 2, 3);
# beta1=1; # R-square would be 0.5;
# beta1=0.5; # R-square would be 0.20;
# beta1=2; # R-square would be 0.8;
# beta1=3; # R-square would be 0.9;
# beta1=1/3; # R-square would be 0.1;
# beta1=1/sqrt(19); # R-square would be 0.05;


# begin the simulation;
# Results from Table 14.4 will be printed out in the cycles

for (l in 1:6)

{ 

beta1=beta1_vec[l];

cat("the l is", l, "\n");
cat("the beta1 is", beta1, "\n");

for (cycle in 1:cycle_no)

{

set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# normal error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the y variable;
y=beta0+beta1*x+error;

# set up the missing data as MCAR;

miss_indi=runif(n=rowobs)<0.40;
y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
delta_no=length(y_miss[!is.na(y_miss)]);
y_delta=y[miss_indi==0];

# deltaerved data is cbind(x, y_miss);

# compelete data analysis;
# regression analysis;
reg_com=summary(lm(y~x));
y_r2=reg_com$r.squared;

# to obtain the statistics at transformed scale; harel (2009);
Q_com=1/2*log((1+sqrt(y_r2))/(1-sqrt(y_r2)));
Q_com_var=1/(rowobs-3);
BD_mean_vec[cycle]=Q_com;
BD_var_vec[cycle]=Q_com_var;

# different missing data methods;

# complete-case analysis;
# mean estimands;
reg_cc=summary(lm(y_miss~x));
y_r2_cc=reg_cc$r.squared;
Q_cc=1/2*log((1+sqrt(y_r2_cc))/(1-sqrt(y_r2_cc)));
Q_cc_var=1/(delta_no-3);

CC_mean_vec[cycle]=Q_cc;
CC_var_vec[cycle]=Q_cc_var;

# now impute the missing y's to get estimates of the marginal;
y_completed=y_miss;

# boot_results_matrix to hold the bootstrap estimates;
boot_results_matrix=matrix(NA, mi_no, boot_no);

for (i in 1:mi_no)
{

# bivariate model imputation; 

y_imputed=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);

y_completed[miss_seq]=y_imputed;

# marginal means;
reg_impute=summary(lm(y_completed~x));
y_r2_impute=reg_impute$r.squared;
Q_impute=1/2*log((1+sqrt(y_r2_impute))/(1-sqrt(y_r2_impute)));
Q_impute_var=1/(rowobs-3);

# bootstrap;
boot_data=data.frame(y_completed,x);
boot_results = boot(data=boot_data, statistic=rsq, R=boot_no, formula=y_completed~x);


# delta method;
delta_mean_mat[cycle,i]=y_r2_impute;
delta_var_mat[cycle,i]=4*y_r2_impute*(1-y_r2_impute)^2/(rowobs-3);

# transformation method;
trans_mean_mat[cycle,i]=Q_impute;
trans_var_mat[cycle,i]=Q_impute_var;

# bootstrap Rubin combine;
bootcombine_mean_mat[cycle,i]=delta_mean_mat[cycle,i];
bootcombine_var_mat[cycle,i]=var(boot_results$t);

# hold bootstrap estimates in the matrix;
boot_results_matrix[i,]=boot_results$t;

}

# obtain the lower and upper percentiles of the bootstrap estimates;
boot_CI_matrix[cycle,]=quantile(boot_results_matrix, c(.025, .975));
# cat("the cycle is", cycle, "\n");

}

BD_mean_vec_tran=((exp(2*BD_mean_vec)-1)/(exp(2*BD_mean_vec)+1))^2;
true=mean(BD_mean_vec_tran);
true;

cat("the true is", true, "\n");

# multiple imputation estimates;
# marginal transformation imputation;

delta_mi_mean=rep(NA, cycle_no);
delta_mi_var=rep(NA, cycle_no);
delta_mi_df=rep(NA, cycle_no);
delta_mi_f=rep(NA, cycle_no);

trans_mi_mean=rep(NA, cycle_no);
trans_mi_var=rep(NA, cycle_no);
trans_mi_df=rep(NA, cycle_no);
trans_mi_f=rep(NA, cycle_no);

bootcombine_mi_mean=rep(NA, cycle_no);
bootcombine_mi_var=rep(NA, cycle_no);
bootcombine_mi_df=rep(NA, cycle_no);
bootcombine_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
delta_summary=pool.scalar(delta_mean_mat[i,], delta_var_mat[i,], n=rowobs);
delta_mi_mean[i]=delta_summary$qbar;
delta_mi_var[i]=delta_summary$t;
delta_mi_df[i]=delta_summary$df;
delta_mi_f[i]=delta_summary$f;

trans_summary=pool.scalar(trans_mean_mat[i,], trans_var_mat[i,], n=rowobs);
trans_mi_mean[i]=trans_summary$qbar;
trans_mi_var[i]=trans_summary$t;
trans_mi_df[i]=trans_summary$df;
trans_mi_f[i]=trans_summary$f;

bootcombine_summary=pool.scalar(bootcombine_mean_mat[i,], bootcombine_var_mat[i,], n=rowobs);
bootcombine_mi_mean[i]=bootcombine_summary$qbar;
bootcombine_mi_var[i]=bootcombine_summary$t;
bootcombine_mi_df[i]=bootcombine_summary$df;
bootcombine_mi_f[i]=bootcombine_summary$f;



}

delta_mi_true=mean(delta_mi_mean);
delta_mi_bias=delta_mi_true-true;

delta_mi_bias;

cat("the delta_mi_bias is", delta_mi_bias, "\n");

# coverage;
delta_mean_low95_vec=delta_mi_mean-qt(.975, delta_mi_df)*sqrt(delta_mi_var);
delta_mean_up95_vec=delta_mi_mean+qt(.975, delta_mi_df)*sqrt(delta_mi_var);

mean_delta_length=mean(delta_mean_up95_vec-delta_mean_low95_vec);
mean_delta_length;

cat("the mean_delta_length is", mean_delta_length, "\n");

delta_mi_coverage=(delta_mean_low95_vec < true)*(delta_mean_up95_vec > true);
mean(delta_mi_coverage);


cat("the delta_mi_coveragee is", mean(delta_mi_coverage), "\n");


# logit transformation;
trans_back_mean=((exp(2*trans_mi_mean)-1)/(exp(2*trans_mi_mean)+1))^2;
trans_mi_true=mean(trans_back_mean);
trans_mi_bias=trans_mi_true-true;

trans_mi_bias;

cat("the trans_mi_bias is", trans_mi_bias, "\n");

# coverage;
trans_mean_low95_vec=trans_mi_mean-qt(.975, trans_mi_df)*sqrt(trans_mi_var);
trans_mean_up95_vec=trans_mi_mean+qt(.975, trans_mi_df)*sqrt(trans_mi_var);

trans_back_low95_vec=((exp(2*trans_mean_low95_vec)-1)/(exp(2*trans_mean_low95_vec)+1))^2;
trans_back_up95_vec=((exp(2*trans_mean_up95_vec)-1)/(exp(2*trans_mean_up95_vec)+1))^2;

mean_trans_length=mean(trans_back_up95_vec-trans_back_low95_vec);
mean_trans_length;


cat("the mean_trans_length is", mean_trans_length, "\n");

trans_mi_coverage=(trans_back_low95_vec < true)*(trans_back_up95_vec > true);
mean(trans_mi_coverage);


cat("the trans_mi_coverage is", mean(trans_mi_coverage), "\n");


# bootstrap first and then combine;

bootcombine_mi_true=mean(bootcombine_mi_mean);
bootcombine_mi_bias=bootcombine_mi_true-true;

bootcombine_mi_bias;

cat("the bootcombine_mi_bias is", bootcombine_mi_bias, "\n");

# coverage;
bootcombine_mean_low95_vec=bootcombine_mi_mean-qt(.975, bootcombine_mi_df)*sqrt(bootcombine_mi_var);
bootcombine_mean_up95_vec=bootcombine_mi_mean+qt(.975, bootcombine_mi_df)*sqrt(bootcombine_mi_var);

mean_bootcombine_length=mean(bootcombine_mean_up95_vec-bootcombine_mean_low95_vec);
mean_bootcombine_length;

cat("the mean_bootcombine_length is", mean_bootcombine_length, "\n");

bootcombine_mi_coverage=(bootcombine_mean_low95_vec < true)*(bootcombine_mean_up95_vec > true);
mean(bootcombine_mi_coverage);


cat("the bootcombine_mi_coverage is", mean(bootcombine_mi_coverage), "\n");

}

# the end of simulation;




