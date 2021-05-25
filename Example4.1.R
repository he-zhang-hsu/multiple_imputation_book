# Example 4.1

rm(list=ls());
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


# generate data;
# mean and variance of the covariate X (from a normal distribution);
mu_x=1;
var_x=1;

# regression parameters;
# both beta0, beta1 and beta2 are fixed at 
beta0=-2;
beta1=1;

# error variance;
var_error=1;


cycle_no=1000;
mi_no=50;

# matrices holding the parameter estimates;
# complete-data inferences;
# xyslope means regressing x on y;
# Before-deletion;
BD_mean_vec=BD_mean_var_vec=BD_slope_vec=BD_slope_var_vec=BD_xyslope_vec=BD_xyslope_var_vec=rep(NA, cycle_no);

# complete-case analysis;
CC_mean_vec=CC_mean_var_vec=CC_slope_vec=CC_slope_var_vec=CC_xyslope_vec=CC_xyslope_var_vec=rep(NA, cycle_no);

# regression prediction;
RP_mean_vec=RP_mean_var_vec=RP_slope_vec=RP_slope_var_vec=RP_xyslope_vec=RP_xyslope_var_vec=rep(NA, cycle_no);

# missing data indicator method;
INDI_xyslope_vec=INDI_xyslope_var_vec=rep(NA, cycle_no);

# multiple imputation using normal linear models;
MI_mean_mat=MI_mean_var_mat=MI_slope_mat=MI_slope_var_mat=MI_xyslope_mat=MI_xyslope_var_mat=matrix(NA, cycle_no, mi_no);

for (cycle in 1:cycle_no)

{

set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# normal error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the variable Y;
y=beta0+beta1*x+error;

# set up the missing data as MCAR on x;
# Results for Table 4.1
# miss_indi=runif(n=rowobs)<0.40;

# set up the missing data as MAR on x;
# Results for Table 4.2
alpha0=-1.4;
alpha1=1;
miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x)/(1+exp(alpha0+alpha1*x)));


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

# reverse regression analysis;
BD_xy=summary(lm(x~y));

BD_xyslope_vec[cycle]=BD_xy$coeff[2,1];
BD_xyslope_var_vec[cycle]=BD_xy$coeff[2,2]^2;


# different missing data methods;

# complete-case analysis;

# marginal mean;
CC_mean_vec[cycle]=mean(y_obs);
CC_mean_var_vec[cycle]=var(y_obs)/obs_no;


# regression analysis;
CC=summary(lm(y_obs~x_obs));

CC_slope_vec[cycle]=CC$coeff[2,1];
CC_slope_var_vec[cycle]=CC$coeff[2,2]^2;

# reverse regression analysis;
CC_xy=summary(lm(x_obs~y_obs));

CC_xyslope_vec[cycle]=CC_xy$coeff[2,1];
CC_xyslope_var_vec[cycle]=CC_xy$coeff[2,2]^2;

# missing data indicator method;
y_miss_0=y_miss;
y_miss_0[miss_indi==1]=0;

INDI_xy=summary(lm(x ~ obs_indi+y_miss_0));
INDI_xyslope_vec[cycle]=INDI_xy$coeff[3,1];
INDI_xyslope_var_vec[cycle]=INDI_xy$coeff[3,2]^2;



# now impute the missing y's;
y_completed=y_miss;

# regression prediction method;
y_imputed=mice.impute.norm.predict(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);


y_completed[miss_seq]=y_imputed;

# completed-data analysis;

# marginal mean;
RP_mean_vec[cycle]=mean(y_completed);
RP_mean_var_vec[cycle]=var(y_completed)/rowobs;

# regression analysis;

IMP=summary(lm(y_completed~x));

RP_slope_vec[cycle]=IMP$coeff[2,1];
RP_slope_var_vec[cycle]=IMP$coeff[2,2]^2;

# reverse regression analysis;

IMP_xy=summary(lm(x~y_completed));

RP_xyslope_vec[cycle]=IMP_xy$coeff[2,1];
RP_xyslope_var_vec[cycle]=IMP_xy$coeff[2,2]^2;


# multiple imputation analysis based on normal linear models;

# now impute the missing y's;
y_completed_MI=y_miss;

# the matrix holds the multiply imputed data;

y_completed_MI_mat=matrix(NA, nrow=rowobs, ncol=mi_no);

for (i in 1:mi_no)
{

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

# reverse regression analysis;

IMP_xy=summary(lm(x~y_completed_MI));

MI_xyslope_mat[cycle,i]=IMP_xy$coeff[2,1];
MI_xyslope_var_mat[cycle,i]=IMP_xy$coeff[2,2]^2;


}


cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

#############################################################################
# Simulation estimates;

# marginal mean of Y;
true_mean=mean(BD_mean_vec);
true_mean;

# before-deletion;
BD=performance(BD_mean_vec, BD_mean_var_vec, true_mean);
BD;

# complete-case analysis;
CC=performance(CC_mean_vec, CC_mean_var_vec, true_mean);
CC;

# RP imputation;
RP=performance(RP_mean_vec, RP_mean_var_vec, true_mean);
RP;

# Multiple imputation analysis;

MI=mi_performance(MI_mean_mat, MI_mean_var_mat, true_mean, rowobs);
MI;



##########################################################################
# slope estimand;
# Regressing Y on X;

true_slope=mean(BD_slope_vec);
true_slope;

# before-deletion;
BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;

# complete-case analysis;
CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;

# RP imputation;
RP=performance(RP_slope_vec, RP_slope_var_vec, true_slope);
RP;

# Multiple imputation;
MI=mi_performance(MI_slope_mat, MI_slope_var_mat, true_slope, rowobs);
MI;

#######################################################################
# slope estimand for regressing X on Y;


true_xyslope=mean(BD_xyslope_vec);
true_xyslope;


# before-deletion; 
BD=performance(BD_xyslope_vec, BD_xyslope_var_vec, true_xyslope);
BD;

# complete-case analysis;

CC=performance(CC_xyslope_vec, CC_xyslope_var_vec, true_xyslope);
CC;

# RP imputation;
RP=performance(RP_xyslope_vec, RP_xyslope_var_vec, true_xyslope);
RP;

# missingness indicator method;
INDI=performance(INDI_xyslope_vec, INDI_xyslope_var_vec, true_xyslope);
INDI;

# multiple imputation estimates;

MI=mi_performance(MI_xyslope_mat, MI_xyslope_var_mat, true_xyslope, rowobs);
MI;



################################################################################# 
#
# plot;
# Use the last simulation replicate;

y_plot_obs=y;
y_plot_miss=y;
y_plot_obs[miss_indi==1]=NA;
y_plot_miss[miss_indi==0]=NA;
y_plot=cbind(y_plot_obs, y_plot_miss);

# plot the observed and missing points;
# choose 250 data points;
# Scatter plots in Fig. 4.2

matplot(x[1:250], y_plot[1:250,], pch=c(1,16), xlab="x", ylab="y", main="Scatter plot before deletion");
# matplot(y_plot[1:250,], x[1:250], pch=c(1,16), xlab="y", ylab="x");

# plot multiply imputed datasets;
# the 1st imputation;
y_plot_MI_1=y_completed_MI_mat[,1];
y_plot_MI_1[miss_indi==0]=NA;
y_plot_MI_1=cbind(y_plot_obs, y_plot_MI_1);
matplot(x[1:250], y_plot_MI_1[1:250,], pch=c(1,17), xlab="x", ylab="y",
   main="Scatter plot of the 1st imputation");

# the 2nd imputation;
y_plot_MI_2=y_completed_MI_mat[,2];
y_plot_MI_2[miss_indi==0]=NA;
y_plot_MI_2=cbind(y_plot_obs, y_plot_MI_2);
matplot(x[1:250], y_plot_MI_2[1:250,], pch=c(1,17), xlab="x", ylab="y",
  main="Scatter plot of the 2nd imputation");

# the 3rd imputation;
y_plot_MI_3=y_completed_MI_mat[,3];
y_plot_MI_3[miss_indi==0]=NA;
y_plot_MI_3=cbind(y_plot_obs, y_plot_MI_3);
matplot(x[1:250], y_plot_MI_3[1:250,], pch=c(1,17), xlab="x", ylab="y",
  main="Scatter plot of the 3rd imputation");

# plot histograms;
# Fig. 4.1

hist(y_plot_miss, xlab="y_miss", main="Histogram of deleted values");
hist(y_plot_MI_1[,2], xlab="y_impute", main="Histogram of the 1st imputation");
hist(y_plot_MI_2[,2], xlab="y_impute", main="Histogram of the 2nd imputation");
hist(y_plot_MI_3[,2], xlab="y_impute", main="Histogram of the 3rd imputation");




