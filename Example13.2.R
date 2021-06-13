# Example 13.2

rm(list=ls());
# library(car);
# library(mvtnorm);
library(mice);
library(norm);
# library(HI);
options(digits=4);

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
rowobs=1000;

# mean and variance of x;
mu_x=1;
var_x=1;

# auxiliary variable z;
mu_z=1;
var_z=3;

# regression parameters; 
beta0=-2;
beta1=1;

# error variance;
var_error=1;

# number of simulations;

cycle_no=1000;

# number of imputations;
mi_no=50;

# matrices holding the parameter estimates;
# xyslope means regressing x on y;
BD_mean_vec=BD_mean_var_vec=BD_slope_vec=BD_slope_var_vec=BD_xyslope_vec=BD_xyslope_var_vec=rep(NA, cycle_no);
CC_mean_vec=CC_mean_var_vec=CC_slope_vec=CC_slope_var_vec=CC_xyslope_vec=CC_xyslope_var_vec=rep(NA, cycle_no);

# exclusive imputation;
MI_mean_mat=MI_mean_var_mat=MI_slope_mat=MI_slope_var_mat=MI_xyslope_mat=MI_xyslope_var_mat=matrix(NA, cycle_no, mi_no);

# inclusive imputation;
EX_mean_mat=EX_mean_var_mat=EX_slope_mat=EX_slope_var_mat=EX_xyslope_mat=EX_xyslope_var_mat=matrix(NA, cycle_no, mi_no);

# record the logistic regression coefficient of missingness indicator on x,y,z;
logit_coef_mat=matrix(NA, cycle_no, 3);

# begin the simulation;

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# normal error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the outcome variable;
y=beta0+beta1*x+error;


# auxiliary variable related to y but not related to the missingness mechanism;
error_z=sqrt(var_error)*rnorm(rowobs);

z=2+y+error_z;

# set up the missing data as MNAR on Y only;
# alpha0=-1.4;
# alpha1=-1;
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*y)/(1+exp(alpha0+alpha1*y)));

# set up the missing data as MNAR depends on both Y and X only;
alpha0=-1.4;
alpha1=-0.75;
alpha2=0.25;
miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*y+alpha2*x)/(1+exp(alpha0+alpha1*y+alpha2*x)));


propensity_logistic=summary(glm(miss_indi~x+y+z, family=binomial(link="logit")));
logit_coef_mat[cycle,]=propensity_logistic$coeff[2:4,1];

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

# multiple imputation analysis;

# now impute the missing y's to get estimates of the marginal;
y_completed_MI=y_completed_EX=y_miss;

# the matrix holds the multiply imputed data;

y_completed_MI_mat=y_completed_EX_mat=matrix(NA, nrow=rowobs, ncol=mi_no);

for (i in 1:mi_no)
{
# exclusive imputation;

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

# Inclusive imputation;

y_imputed_EX=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x,z));

y_completed_EX[miss_seq]=y_imputed_EX;

y_completed_EX_mat[,i]=y_completed_EX;

# completed-data analysis;

# marginal mean;
EX_mean_mat[cycle,i]=mean(y_completed_EX);
EX_mean_var_mat[cycle,i]=var(y_completed_EX)/rowobs;

# regression analysis;

IMP=summary(lm(y_completed_EX~x));

EX_slope_mat[cycle,i]=IMP$coeff[2,1];
EX_slope_var_mat[cycle,i]=IMP$coeff[2,2]^2;

# reverse regression analysis;

IMP_xy=summary(lm(x~y_completed_EX));

EX_xyslope_mat[cycle,i]=IMP_xy$coeff[2,1];
EX_xyslope_var_mat[cycle,i]=IMP_xy$coeff[2,2]^2;



}


cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

# plot the distribution of the data;

# check the performance of the estimates;
# population quantity

#############################################################################
# Simulation results;
# Table 13.1
# marginal mean estimand;

true_mean=mean(BD_mean_vec);
true_mean;


BD=performance(BD_mean_vec, BD_mean_var_vec, true_mean);
BD;

CC=performance(CC_mean_vec, CC_mean_var_vec, true_mean);
CC;

# Exclusive imputation;
MI=mi_performance(MI_mean_mat, MI_mean_var_mat, true_mean, rowobs);
MI;

# Inclusive imputation
EX=mi_performance(EX_mean_mat, EX_mean_var_mat, true_mean, rowobs);
EX;


##########################################################################
# slope estimand;

# complete-data analysis;

true_slope=mean(BD_slope_vec);
true_slope;

BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;

CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;

# Exclusive imputation;
MI=mi_performance(MI_slope_mat, MI_slope_var_mat, true_slope, rowobs);
MI;

# Inclusive imputation;
EX=mi_performance(EX_slope_mat, EX_slope_var_mat, true_slope, rowobs);
EX;



#######################################################################
# reverse slope;

# complete-data analysis;

true_xyslope=mean(BD_xyslope_vec);
true_xyslope;

BD=performance(BD_xyslope_vec, BD_xyslope_var_vec, true_xyslope);
BD;
CC=performance(CC_xyslope_vec, CC_xyslope_var_vec, true_xyslope);
CC;

# Exclusive imputation;
MI=mi_performance(MI_xyslope_mat, MI_xyslope_var_mat, true_xyslope, rowobs);
MI;

# Inclusive imputation
EX=mi_performance(EX_xyslope_mat, EX_xyslope_var_mat, true_xyslope, rowobs);
EX;
