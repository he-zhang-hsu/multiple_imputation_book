# Example 11.1

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

# error variance;
var_error=1;

# measurement error model;
# Scenario I
# Correspond to Table 11.2
c0=0;
c1=1;
c2=1;

# Scenario II
# Correspond to Table 11.3
# c0=-0.2;
# c1=1.2;
# c2=1;

# number of simulations;
cycle_no=1000;

# number of imputations;
mi_no=80;

# Vectors and matrices holding the parameter estimates;
# xyslope means regressing x on y;

# BD means using true y;
BD_mean_vec=BD_mean_var_vec=BD_slope_vec=BD_slope_var_vec=BD_xyslope_vec=BD_xyslope_var_vec=rep(NA, cycle_no);

# CC means using y from validation data;
# The validation method;
CC_mean_vec=CC_mean_var_vec=CC_slope_vec=CC_slope_var_vec=CC_xyslope_vec=CC_xyslope_var_vec=rep(NA, cycle_no);

# NA (naive method) means using mismeasued y;
NA_mean_vec=NA_mean_var_vec=NA_slope_vec=NA_slope_var_vec=NA_xyslope_vec=NA_xyslope_var_vec=rep(NA, cycle_no);

# MI means using multiply imputed y;
MI_mean_mat=MI_mean_var_mat=MI_slope_mat=MI_slope_var_mat=MI_xyslope_mat=MI_xyslope_var_mat=matrix(NA, cycle_no, mi_no);

# begin the simulation;

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the true y variable;
y_true=beta0+beta1*x+error;

# y is observed after applying some error;
y=c0+c1*y_true+c2*rnorm(rowobs);


# set up the missing data as MCAR on x;
miss_indi=runif(n=rowobs)<0.80;

# set up the missing data as MAR on x;
# alpha0=-1.4;
# alpha1=1;
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x)/(1+exp(alpha0+alpha1*x)));


obs_indi=1-miss_indi;

y_miss=y_true;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_true_obs=y_miss[miss_indi==0];
x_obs=x[miss_indi==0];
x_mis=x[miss_indi==1];

# compelete data (using true y) analysis;
# marginal mean;
BD_mean_vec[cycle]=mean(y_true);
BD_mean_var_vec[cycle]=var(y_true)/rowobs;


# regression analysis;
BD=summary(lm(y_true~x));

BD_slope_vec[cycle]=BD$coeff[2,1];
BD_slope_var_vec[cycle]=BD$coeff[2,2]^2;

# reverse regression analysis;
BD_xy=summary(lm(x~y_true));

BD_xyslope_vec[cycle]=BD_xy$coeff[2,1];
BD_xyslope_var_vec[cycle]=BD_xy$coeff[2,2]^2;


# different missing data methods;

# complete-case analysis;

# marginal mean;
CC_mean_vec[cycle]=mean(y_true_obs);
CC_mean_var_vec[cycle]=var(y_true_obs)/obs_no;


# regression analysis;
CC=summary(lm(y_true_obs~x_obs));

CC_slope_vec[cycle]=CC$coeff[2,1];
CC_slope_var_vec[cycle]=CC$coeff[2,2]^2;

# reverse regression analysis;
CC_xy=summary(lm(x_obs~y_true_obs));

CC_xyslope_vec[cycle]=CC_xy$coeff[2,1];
CC_xyslope_var_vec[cycle]=CC_xy$coeff[2,2]^2;

# using mismeasured y method;

# marginal mean;
NA_mean_vec[cycle]=mean(y);
NA_mean_var_vec[cycle]=var(y)/rowobs;

# regression analysis;

IMP=summary(lm(y~x));

NA_slope_vec[cycle]=IMP$coeff[2,1];
NA_slope_var_vec[cycle]=IMP$coeff[2,2]^2;

# reverse regression analysis;

IMP_xy=summary(lm(x~y));

NA_xyslope_vec[cycle]=IMP_xy$coeff[2,1];
NA_xyslope_var_vec[cycle]=IMP_xy$coeff[2,2]^2;


# multiple imputation analysis;

# now impute the missing y_true's to get estimates of the marginal;
y_completed_MI=y_miss;

# the matrix holds the multiply imputed data;
# x_new is the observed covariates;
# including both x and mismeasured y;

x_new=cbind(x, y);
y_completed_MI_mat=matrix(NA, nrow=rowobs, ncol=mi_no);

for (i in 1:mi_no)
{

set.seed(i);

y_imputed_MI=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x_new);

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

# reverse regression analysis;

IMP_xy=summary(lm(x~y_completed_MI));

MI_xyslope_mat[cycle,i]=IMP_xy$coeff[2,1];
MI_xyslope_var_mat[cycle,i]=IMP_xy$coeff[2,2]^2;


}


cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;


#############################################################################
# Simulation results
# marginal mean;

true_mean=mean(BD_mean_vec);
true_mean;


BD=performance(BD_mean_vec, BD_mean_var_vec, true_mean);
BD;
CC=performance(CC_mean_vec, CC_mean_var_vec, true_mean);
CC;
NA_results=performance(NA_mean_vec, NA_mean_var_vec, true_mean);
NA_results;

MI=mi_performance(MI_mean_mat, MI_mean_var_mat, true_mean, rowobs);
MI;

##########################################################################
# slope estimand;
# regressing y on x;
# 

true_slope=mean(BD_slope_vec);
true_slope;

BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;
CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;

NA_results=performance(NA_slope_vec, NA_slope_var_vec, true_slope);
NA_results;

MI=mi_performance(MI_slope_mat, MI_slope_var_mat, true_slope, rowobs);
MI;


#######################################################################
# reverse slope;
# regressin x on y;

true_xyslope=mean(BD_xyslope_vec);
true_xyslope;

BD=performance(BD_xyslope_vec, BD_xyslope_var_vec, true_xyslope);
BD;
CC=performance(CC_xyslope_vec, CC_xyslope_var_vec, true_xyslope);
CC;

NA_results=performance(NA_xyslope_vec, NA_xyslope_var_vec, true_xyslope);
NA_results;

MI=mi_performance(MI_xyslope_mat, MI_xyslope_var_mat, true_xyslope, rowobs);
MI;




################################################################################# 

# plot;
# Fig. 11.1
# Use the last simulation from Scenario I;

# Top left;
plot(y_true[1:250], x[1:250], xlab="true y", xlim=c(-5.5, 4), ylab="x");
BD_fit=lm(x[1:250]~y_true[1:250]);
abline(BD_fit);

# Top right;
plot(y[1:250], x[1:250], ylab="x", xlim=c(-5.5, 4), pch=16, xlab="mismeasured y");
NA_fit=lm(x[1:250]~y[1:250]);
abline(NA_fit);

# The bottom
# plot the observed and imputed points;
y_plot_obs=y_miss;

y_plot_imputed=y_completed_MI;
y_plot_imputed[miss_indi==0]=NA;
y_plot_impute=cbind(y_plot_obs, y_plot_imputed);
x_1=x;
x_1[miss_indi==1]=NA;
x_2=x;
x_2[miss_indi==0]=NA;
x_plot_impute=cbind(x_1, x_2);


matplot(y_completed_MI[1:250], x_plot_impute[1:250,], pch=c(1,17), xlim=c(-5.5, 4), ylab="x", xlab="completed y");
impute_fit=lm(x[1:250]~y_completed_MI[1:250]);
abline(impute_fit);


