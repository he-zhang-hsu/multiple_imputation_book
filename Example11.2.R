# Example 11.2

rm(list=ls());
# library(car);
library(mice);
library(norm);
# adaptive metroplis rejection sampling;
library(HI);
# include the library MASS for generating from multivariate normal distribution;
library(MASS);
# include the library msm and bayesm for using the function of drawing truncated normal
# distribution;
library(msm);
library(bayesm);
library(mvtnorm);
library(sn);
library(arm);
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
var_x=1/4;

# regression parameters; 
beta0=-2;
beta1=2; 

# number of simulations;
cycle_no=1000;

# number of imputations;
mi_no=80;

# sensitivity and specificity;
sens=0.9;
spec=0.8;

# Vectors and matrices holding the parameter estimates;
BD_mean_vec=BD_mean_var_vec=rep(NA, cycle_no);

# CC is the valdiation method;
CC_mean_vec=CC_mean_var_vec=rep(NA, cycle_no);

# NA is the naive method;
NA_mean_vec=NA_mean_var_vec=rep(NA, cycle_no);

# MI is the multiple imputation method;
MI_mean_mat=MI_mean_var_mat=matrix(NA, cycle_no, mi_no);

BD_slope_vec=BD_slope_var_vec=rep(NA, cycle_no);
CC_slope_vec=CC_slope_var_vec=rep(NA, cycle_no);
NA_slope_vec=NA_slope_var_vec=rep(NA, cycle_no);
MI_slope_mat=MI_slope_var_mat=matrix(NA, cycle_no, mi_no);

# begin the simulation;

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# generate the logistic regression outcome y;
y_true=rbinom(n=rowobs, size=1, prob=exp(beta0+beta1*x)/(1+exp(beta0+beta1*x)));

# assign the misclassification to y_true;
# y is the observed binary data;
y=rep(NA, rowobs);
alpha0=log(1/spec-1);
alpha1=log(1/(1-sens)-1)-log(1/spec-1);
y=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*y_true)/(1+exp(alpha0+alpha1*y_true)));

# set up the missing data as MCAR on x;

miss_indi=runif(n=rowobs)<0.80;
y_miss=y_true;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y_true[miss_indi==0];
x_obs=x[miss_indi==0];
x_miss=x[miss_indi==1];

# observed data is cbind(x, y_miss);

# compelete data analysis;
# sample mean;
y_mean=mean(y_true);
y_var=var(y_true)/rowobs;
BD_mean_vec[cycle]=y_mean;
BD_mean_var_vec[cycle]=y_var;

# logistic regression coefficient for slope;
BD_logistic=summary(glm(y_true~x, family=binomial(link="logit")));
BD_slope_vec[cycle]=BD_logistic$coeff[2,1];
BD_slope_var_vec[cycle]=BD_logistic$coeff[2,2]^2;


# different missing data methods;

# complete-case analysis;
# use the validation sample
# mean estimands;
cc_mean=mean(y_miss, na.rm="T");
CC_mean_vec[cycle]=cc_mean;
CC_mean_var_vec[cycle]=var(y_miss, na.rm="T")/obs_no;

# logistic regression coefficient for slope;
CC_logistic=summary(glm(y_obs~x_obs, family=binomial(link="logit")));
CC_slope_vec[cycle]=CC_logistic$coeff[2,1];
CC_slope_var_vec[cycle]=CC_logistic$coeff[2,2]^2;

# use the misreported version;
NA_mean_vec[cycle]=mean(y);
NA_mean_var_vec[cycle]=var(y)/rowobs;

# logistic regression coefficient for slope;
NA_logistic=summary(glm(y~x, family=binomial(link="logit")));
NA_slope_vec[cycle]=NA_logistic$coeff[2,1];
NA_slope_var_vec[cycle]=NA_logistic$coeff[2,2]^2;

# now impute the missing y's to get estimates of the marginal;
y_completed_MI=y_miss;

# The imputation predictor uses
# both x and misreported y;
x_new=cbind(x, y);

for (i in 1:mi_no)
{

set.seed(i);

# logistic regression model imputation; 
y_imputed_MI=mice.impute.logreg(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x_new);

y_completed_MI[miss_seq]=y_imputed_MI;

# marginal means;
MI_mean=mean(y_completed_MI);
MI_var=var(y_completed_MI)/rowobs;

MI_mean_mat[cycle,i]=MI_mean;
MI_mean_var_mat[cycle,i]=MI_var;

# logistic regression coefficient for slope;
MI_logistic=summary(glm(y_completed_MI~x, family=binomial(link="logit")));
MI_slope_mat[cycle,i]=MI_logistic$coeff[2,1];
MI_slope_var_mat[cycle,i]=MI_logistic$coeff[2,2]^2;

}

cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

#############################################################################
# Simulation results
# Table 11.5
# marginal mean estimand;

true=mean(BD_mean_vec);
true;

true;


BD=performance(BD_mean_vec, BD_mean_var_vec, true);
BD;

CC=performance(CC_mean_vec, CC_mean_var_vec, true);
CC;

NA_results=performance(NA_mean_vec, NA_mean_var_vec, true);
NA_results;

MI=mi_performance(MI_mean_mat, MI_mean_var_mat, true, rowobs);
MI;





##########################################################################
# slope estimand;

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


