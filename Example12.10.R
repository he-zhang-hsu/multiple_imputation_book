# Example 12.10

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
estimand=as.data.frame(cbind(bias,rbias,mse,mean_length,mean_cov,sd,se))
}

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
beta3=1;

# error variance;
var_error=1;

# number of simulations;
cycle_no=1000;

# number of multiple imputations;
mi_no=50;

# Vectors and matrices holding the parameter estimates;
BD_mean_vec=BD_mean_var_vec=CC_mean_vec=CC_mean_var_vec=p_mis_vec=rep(NA, cycle_no);
MI_mean_x1_mat=MI_mean_x1_var_mat=MI_mean_x1plusx2_mat=MI_mean_x1plusx2_var_mat=MI_mean_x1x2_mat=MI_mean_x1x2_var_mat=matrix(NA, cycle_no, mi_no);
MI_rsquare_x1_mat=MI_rsquare_x1plusx2_mat=MI_rsquare_x1x2_mat=matrix(NA, cycle_no, mi_no);



# begin the simulation;

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x1=mu_x+sqrt(var_x)*rnorm(rowobs);
x2=mu_x+sqrt(var_x)*rnorm(rowobs);


# normal error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the y variable;
y=beta0+beta1*x1+beta2*x2+beta3*x1*x2+error;

# set up the missing data as MCAR on x;
miss_indi=runif(n=rowobs)<0.40;

# set up the missing data as MAR on x;
# alpha0=-1.4;
# alpha1=0.3;
# alpha2=0.7;
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x1+alpha2*x2)/(1+exp(alpha0+alpha1*x1+alpha2*x2)));


y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
# x_obs=x[miss_indi==0];
# x_mis=x[miss_indi==1];
p_mis_vec[cycle]=mis_no/rowobs;

# compelete data analysis;
# marginal mean;
BD_mean_vec[cycle]=mean(y);
BD_mean_var_vec[cycle]=var(y)/rowobs;


# different missing data methods;

# complete-case analysis;

# marginal mean;
CC_mean_vec[cycle]=mean(y_obs);
CC_mean_var_vec[cycle]=var(y_obs)/obs_no;

# now impute the missing y's;
y_completed_x1=y_completed_x1plusx2=y_completed_x1x2=y_miss;

for (i in 1:mi_no)
{
# imputation model only includes x1;
y_imputed_x1=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x1);

y_completed_x1[miss_seq]=y_imputed_x1;

# completed-data analysis;

# marginal mean;
MI_mean_x1_mat[cycle,i]=mean(y_completed_x1);
MI_mean_x1_var_mat[cycle,i]=var(y_completed_x1)/rowobs;

# r-square;
MI_rsquare_x1_mat[cycle,i]=summary(lm(y_completed_x1~x1))$r.squared;

# imputation model includes both x1 and x2;
y_imputed_x1plusx2=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x1,x2));

y_completed_x1plusx2[miss_seq]=y_imputed_x1plusx2;

# completed-data analysis;

# marginal mean;
MI_mean_x1plusx2_mat[cycle,i]=mean(y_completed_x1plusx2);
MI_mean_x1plusx2_var_mat[cycle,i]=var(y_completed_x1plusx2)/rowobs;

# r-square;
MI_rsquare_x1plusx2_mat[cycle,i]=summary(lm(y_completed_x1plusx2~x1+x2))$r.squared;

# imputation model includes both x1 and x2 and their interactions;
y_imputed_x1x2=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x1,x2, x1*x2));

y_completed_x1x2[miss_seq]=y_imputed_x1x2;

# completed-data analysis;

# marginal mean;
MI_mean_x1x2_mat[cycle,i]=mean(y_completed_x1x2);
MI_mean_x1x2_var_mat[cycle,i]=var(y_completed_x1x2)/rowobs;

# r-square;
MI_rsquare_x1x2_mat[cycle,i]=summary(lm(y_completed_x1x2~x1*x2))$r.squared;

}

cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

#############################################################################
# marginal mean estimand;
# Table 12.5
true_mean=mean(BD_mean_vec);
true_mean;

BD=performance(BD_mean_vec, BD_mean_var_vec, true_mean);
BD;

CC=performance(CC_mean_vec, CC_mean_var_vec, true_mean);
CC;

# multiple imputation estimates;
# For the model with only x1;

mean_mi_x1_mean=rep(NA, cycle_no);
mean_mi_x1_var=rep(NA, cycle_no);
mean_mi_x1_df=rep(NA, cycle_no);
mean_mi_x1_f=rep(NA, cycle_no);


for (i in 1:cycle_no)
{
mean_summary=pool.scalar(MI_mean_x1_mat[i,], MI_mean_x1_var_mat[i,], n=rowobs);
mean_mi_x1_mean[i]=mean_summary$qbar;
mean_mi_x1_var[i]=mean_summary$t;
mean_mi_x1_df[i]=mean_summary$df;
mean_mi_x1_f[i]=mean_summary$fmi;


}


# For marginal means;

mean_mi_x1_bias=mean(mean_mi_x1_mean)-true_mean;
mean_mi_x1_bias;

mean_mi_x1_bias/true_mean;

mean_mi_x1_mse=mean((mean_mi_x1_mean-true_mean)^2);
mean_mi_x1_mse;

# coverage;
MI_x1_mean_low95_vec=mean_mi_x1_mean-qt(.975, mean_mi_x1_df)*sqrt(mean_mi_x1_var);
MI_x1_mean_up95_vec=mean_mi_x1_mean+qt(.975, mean_mi_x1_df)*sqrt(mean_mi_x1_var);

MI_x1_mean_length=mean(MI_x1_mean_up95_vec-MI_x1_mean_low95_vec);

MI_x1_mean_length;

MI_x1_mean_coverage=(MI_x1_mean_low95_vec < true_mean)*(MI_x1_mean_up95_vec > true_mean);
mean(MI_x1_mean_coverage);

# variance estimates;
# SE
mean(sqrt(mean_mi_x1_var));
sqrt(var(mean_mi_x1_mean));

mean_mi_x1_rsquare=rowMeans(MI_rsquare_x1_mat);

# R-square from completed data;
mean(mean_mi_x1_rsquare);

mean_mi_x1_f_approx=p_mis_vec*(1-mean_mi_x1_rsquare)/(1-p_mis_vec*mean_mi_x1_rsquare);

# FMI calculated by Rubin's formula
mean(mean_mi_x1_f);

# FMI related to r-square;
mean(mean_mi_x1_f_approx);

# plot(mean_mi_x1_f, mean_mi_x1_f_approx);


# multiple imputation estimates;
# For the model with x1 and x2;
mean_mi_x1plusx2_mean=rep(NA, cycle_no);
mean_mi_x1plusx2_var=rep(NA, cycle_no);
mean_mi_x1plusx2_df=rep(NA, cycle_no);
mean_mi_x1plusx2_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
mean_summary=pool.scalar(MI_mean_x1plusx2_mat[i,], MI_mean_x1plusx2_var_mat[i,], n=rowobs);
mean_mi_x1plusx2_mean[i]=mean_summary$qbar;
mean_mi_x1plusx2_var[i]=mean_summary$t;
mean_mi_x1plusx2_df[i]=mean_summary$df;
mean_mi_x1plusx2_f[i]=mean_summary$fmi;


}


# For marginal means;

mean_mi_x1plusx2_bias=mean(mean_mi_x1plusx2_mean)-true_mean;
mean_mi_x1plusx2_bias;

mean_mi_x1plusx2_bias/true_mean;

mean_mi_x1plusx2_mse=mean((mean_mi_x1plusx2_mean-true_mean)^2);
mean_mi_x1plusx2_mse;

# coverage;
MI_x1plusx2_mean_low95_vec=mean_mi_x1plusx2_mean-qt(.975, mean_mi_x1plusx2_df)*sqrt(mean_mi_x1plusx2_var);
MI_x1plusx2_mean_up95_vec=mean_mi_x1plusx2_mean+qt(.975, mean_mi_x1plusx2_df)*sqrt(mean_mi_x1plusx2_var);

MI_x1plusx2_mean_length=mean(MI_x1plusx2_mean_up95_vec-MI_x1plusx2_mean_low95_vec);

MI_x1plusx2_mean_length;

MI_x1plusx2_mean_coverage=(MI_x1plusx2_mean_low95_vec < true_mean)*(MI_x1plusx2_mean_up95_vec > true_mean);
mean(MI_x1plusx2_mean_coverage);

# variance estimates;
# SE;
mean(sqrt(mean_mi_x1plusx2_var));
sqrt(var(mean_mi_x1plusx2_mean));

# R-square from completed data
mean_mi_x1plusx2_rsquare=rowMeans(MI_rsquare_x1plusx2_mat);

mean(mean_mi_x1plusx2_rsquare);

mean_mi_x1plusx2_f_approx=p_mis_vec*(1-mean_mi_x1plusx2_rsquare)/(1-p_mis_vec*mean_mi_x1plusx2_rsquare);

# FMI calculated by Rubin's formula;
mean(mean_mi_x1plusx2_f);

# FMI related to R-square
mean(mean_mi_x1plusx2_f_approx);

# plot(mean_mi_x1plusx2_f, mean_mi_x1plusx2_f_approx);

# multiple imputation estimates;
# with x1, x2, and x1*x2;

mean_mi_x1x2_mean=rep(NA, cycle_no);
mean_mi_x1x2_var=rep(NA, cycle_no);
mean_mi_x1x2_df=rep(NA, cycle_no);
mean_mi_x1x2_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
mean_summary=pool.scalar(MI_mean_x1x2_mat[i,], MI_mean_x1x2_var_mat[i,], n=rowobs);
mean_mi_x1x2_mean[i]=mean_summary$qbar;
mean_mi_x1x2_var[i]=mean_summary$t;
mean_mi_x1x2_df[i]=mean_summary$df;
mean_mi_x1x2_f[i]=mean_summary$fmi;


}


# For marginal means;

mean_mi_x1x2_bias=mean(mean_mi_x1x2_mean)-true_mean;
mean_mi_x1x2_bias;

mean_mi_x1x2_bias/true_mean;

mean_mi_x1x2_mse=mean((mean_mi_x1x2_mean-true_mean)^2);
mean_mi_x1x2_mse;

# coverage;
MI_x1x2_mean_low95_vec=mean_mi_x1x2_mean-qt(.975, mean_mi_x1x2_df)*sqrt(mean_mi_x1x2_var);
MI_x1x2_mean_up95_vec=mean_mi_x1x2_mean+qt(.975, mean_mi_x1x2_df)*sqrt(mean_mi_x1x2_var);

MI_x1x2_mean_length=mean(MI_x1x2_mean_up95_vec-MI_x1x2_mean_low95_vec);

MI_x1x2_mean_length;

MI_x1x2_mean_coverage=(MI_x1x2_mean_low95_vec < true_mean)*(MI_x1x2_mean_up95_vec > true_mean);
mean(MI_x1x2_mean_coverage);

# variance estimates;
# SE
mean(sqrt(mean_mi_x1x2_var));

sqrt(var(mean_mi_x1x2_mean));

# R-square from completed data
mean_mi_x1x2_rsquare=rowMeans(MI_rsquare_x1x2_mat);

mean(mean_mi_x1x2_rsquare);


mean_mi_x1x2_f_approx=p_mis_vec*(1-mean_mi_x1x2_rsquare)/(1-p_mis_vec*mean_mi_x1x2_rsquare);

# FMI calculated by Rubin's formula;
mean(mean_mi_x1x2_f);

# FMI related to R-square
mean(mean_mi_x1x2_f_approx);



