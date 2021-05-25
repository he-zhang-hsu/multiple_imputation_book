# Example 4.9

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

# complete-data sample size in the simulation;
rowobs=1000;

# generate data;
# mean and variance of z;
mu_z=0;
var_z=1;

# regression parameters for X given Z;
alpha0=1;
alpha1=1;

# regression parameters for Y given X and Z; 
beta0=-2;
beta1=1;
beta2=1;

# error variance;
# regression error for Y on X and Z;
var_y_xz_error=1;

# regression error for X on Z;
var_x_z_error=1;

cycle_no=1000;
mi_no=50;

# Vectors and matrices holding the parameter estimates;
# complete-data inferences;
# Before-deletion analysis;
BD_mean_vec=BD_mean_var_vec=BD_slope_vec=BD_slope_var_vec=BD_xyslope_vec=BD_xyslope_var_vec=rep(NA, cycle_no);

# Complete-case analysis;
CC_mean_vec=CC_mean_var_vec=CC_slope_vec=CC_slope_var_vec=CC_xyslope_vec=CC_xyslope_var_vec=rep(NA, cycle_no);

# Correct multiple imputation including the outcome;
MI_mean_mat=MI_mean_var_mat=MI_slope_mat=MI_slope_var_mat=MI_xyslope_mat=MI_xyslope_var_mat=matrix(NA, cycle_no, mi_no);

# Incorrect multiple imputation excluding the outcome;
IM_MI_mean_mat=IM_MI_mean_var_mat=IM_MI_slope_mat=IM_MI_slope_var_mat=IM_MI_xyslope_mat=IM_MI_xyslope_var_mat=matrix(NA, cycle_no, mi_no);

# begin the simulation;
for (cycle in 1:cycle_no)

{
set.seed(cycle);

# the distribution of covariate z;
z=mu_z+sqrt(var_z)*rnorm(rowobs);

# the distribution of covariate x;
x=alpha0+alpha1*z+sqrt(var_x_z_error)*rnorm(rowobs);

# generate the outcome variable;
y=beta0+beta1*x+beta2*z+sqrt(var_y_xz_error)*rnorm(rowobs);

# set up the missing data as MCAR on x;
miss_indi=runif(n=rowobs)<0.30;

# set up the missing data as MAR on x;
# alpha0=-1.4;
# alpha1=1;
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x)/(1+exp(alpha0+alpha1*x)));


x_mis=x;
x_mis[miss_indi==1]=NA;
x_obs=x[miss_indi==0];

miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(x_mis[!is.na(x_mis)]);
mis_no=rowobs-obs_no;

y_obs=y[miss_indi==0];
y_mis=y[miss_indi==1];
z_obs=z[miss_indi==0];
z_mis=z[miss_indi==1];

# compelete data analysis;
# marginal mean of X;
BD_mean_vec[cycle]=mean(x);
BD_mean_var_vec[cycle]=var(x)/rowobs;


# regression analysis for Y given X and Z;
BD=summary(lm(y~x+z));

# the slope on x;
BD_slope_vec[cycle]=BD$coeff[2,1];
BD_slope_var_vec[cycle]=BD$coeff[2,2]^2;

# the slope on z;
BD_xyslope_vec[cycle]=BD$coeff[3,1];
BD_xyslope_var_vec[cycle]=BD$coeff[3,2]^2;


# different missing data methods;

# complete-case analysis;

# marginal mean of X;
CC_mean_vec[cycle]=mean(x_obs);
CC_mean_var_vec[cycle]=var(x_obs)/obs_no;


# regression analysis for Y given X and Z;
CC=summary(lm(y_obs~x_obs+z_obs));

# the slope coefficient on x;
CC_slope_vec[cycle]=CC$coeff[2,1];
CC_slope_var_vec[cycle]=CC$coeff[2,2]^2;

# the slope coefficient on z;
CC_xyslope_vec[cycle]=CC$coeff[3,1];
CC_xyslope_var_vec[cycle]=CC$coeff[3,2]^2;

# now impute the missing Y's to get estimates of the marginal;
x_completed=x_completed_IM=x_mis;

for (i in 1:mi_no)
{
set.seed(i);

# correct model imputation by including y and z;
x_imputed=mice.impute.norm(x_mis, ry=as.logical(1-miss_indi), seed=i, x=cbind(y,z));

x_completed[miss_seq]=x_imputed;

# completed-data analysis;

# marginal mean of X;
MI_mean_mat[cycle,i]=mean(x_completed);
MI_mean_var_mat[cycle,i]=var(x_completed)/rowobs;

# regression analysis for Y given X and Z;

IMP=summary(lm(y~x_completed+z));

# the slope coefficient on x;
MI_slope_mat[cycle,i]=IMP$coeff[2,1];
MI_slope_var_mat[cycle,i]=IMP$coeff[2,2]^2;

# the slope coefficeint on z;
MI_xyslope_mat[cycle,i]=IMP$coeff[3,1];
MI_xyslope_var_mat[cycle,i]=IMP$coeff[3,2]^2;

# incorrect imputation ignoring y;
x_imputed_IM=mice.impute.norm(x_mis, ry=as.logical(1-miss_indi), seed=i, x=as.matrix(z));

x_completed_IM[miss_seq]=x_imputed_IM;

# completed-data analysis;

# marginal mean of X;
IM_MI_mean_mat[cycle,i]=mean(x_completed_IM);
IM_MI_mean_var_mat[cycle,i]=var(x_completed_IM)/rowobs;

# regression analysis of Y given X and Z;

IM_IMP=summary(lm(y~x_completed_IM+z));

# the slope coefficient for x;
IM_MI_slope_mat[cycle,i]=IM_IMP$coeff[2,1];
IM_MI_slope_var_mat[cycle,i]=IM_IMP$coeff[2,2]^2;

# the slope coefficient for z;

IM_MI_xyslope_mat[cycle,i]=IM_IMP$coeff[3,1];
IM_MI_xyslope_var_mat[cycle,i]=IM_IMP$coeff[3,2]^2;



}

cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

# plot the distribution of the data;

# check the performance of the estimates;
# population quantity

#############################################################################
# Simulation performance
# marginal mean for x;
# For Table 4.14
true_mean=mean(BD_mean_vec);
true_mean;

# before-deletion;
BD=performance(BD_mean_vec, BD_mean_var_vec, true_mean);
BD;

# complete-case analysis;
CC=performance(CC_mean_vec, CC_mean_var_vec, true_mean);
CC;

# Multiple imputation analysis;
# Correct multiple imputation using the outcome variable

MI=mi_performance(MI_mean_mat, MI_mean_var_mat, true_mean, rowobs);
MI;

# Incorrect multiple imputation ignoring the outcome variable
IM_MI=mi_performance(IM_MI_mean_mat, IM_MI_mean_var_mat, true_mean, rowobs);
IM_MI;


##########################################################################
# slope estimand for beta1;
# Table 4.14
# complete-data analysis;

true_slope=mean(BD_slope_vec);
true_slope;

# before-deletion;
BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;

# complete-case analysis;
CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;

# Multiple imputation analysis;
# Correct multiple imputation using the outcome variable

MI=mi_performance(MI_slope_mat, MI_slope_var_mat, true_slope, rowobs);
MI;

# Incorrect multiple imputation ignoring the outcome variable
IM_MI=mi_performance(IM_MI_slope_mat, IM_MI_slope_var_mat, true_slope, rowobs);
IM_MI;

##########################################################################
# slope estimand for beta2;
# Table 4.14
# complete-data analysis;

true_slope=mean(BD_xyslope_vec);
true_slope;

# before-deletion;
BD=performance(BD_xyslope_vec, BD_xyslope_var_vec, true_slope);
BD;

# complete-case analysis;
CC=performance(CC_xyslope_vec, CC_xyslope_var_vec, true_slope);
CC;

# Multiple imputation analysis;
# Correct multiple imputation using the outcome variable

MI=mi_performance(MI_xyslope_mat, MI_xyslope_var_mat, true_slope, rowobs);
MI;

# Incorrect multiple imputation ignoring the outcome variable
IM_MI=mi_performance(IM_MI_xyslope_mat, IM_MI_xyslope_var_mat, true_slope, rowobs);
IM_MI;


