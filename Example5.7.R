# Example 5.7

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
# library(MCMCpack);
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

# complete-data sample size
rowobs=1000;

# number of simulations;
cycle_no=1000;

# number of multiple imputations;
mi_no=30;

# number of MCMC iterations;
draw_no=1000;

# Vector and matrices holding the parameter estimates;
# For mean estimates;
BD_mean_vec=BD_mean_var_vec=CC_mean_vec=CC_mean_var_vec=rep(NA, cycle_no);

# Multiple imputation;

# ROUND: linear imputation using only the main effect of x and round;
# PMM_linear: PMM imputation at the linear scale;
# PMM_logit: PMM imputation at the logit scale;
ROUND_mean_mat=ROUND_mean_var_mat=matrix(NA, cycle_no, mi_no);
PMM_linear_mean_mat=PMM_linear_mean_var_mat=PMM_logit_mean_mat=PMM_logit_mean_var_mat=matrix(NA, cycle_no, mi_no);

# For linear regression slope;
BD_slope_vec=BD_slope_var_vec=rep(NA, cycle_no);
CC_slope_vec=CC_slope_var_vec=rep(NA, cycle_no);

# Multiple imputation;
ROUND_slope_mat=ROUND_slope_var_mat=matrix(NA, cycle_no, mi_no);
PMM_linear_slope_mat=PMM_linear_slope_var_mat=PMM_logit_slope_mat=PMM_logit_slope_var_mat=matrix(NA, cycle_no, mi_no);


# begin the simulation;

for (cycle in 1:cycle_no)

{
set.seed(cycle);


# the following data generation is based on White et al. (2010) simulation;
# Also see Example 4.6.R
y=rbinom(n=rowobs, size=1, prob=0.1);
zeros_y=rowobs-sum(y); # number of ys are 0;
x=rep(0,rowobs);
x[y==0]=rbinom(n=zeros_y, size=1, prob=0.8);
u=rnorm(rowobs);
z=y+x+0.5*u+2*rnorm(rowobs);


# set up the missing data as MCAR on x;


miss_indi=runif(n=rowobs)<0.30;
y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
x_obs=x[miss_indi==0];
x_miss=x[miss_indi==1];


u_obs=u[miss_indi==0];
u_miss=u[miss_indi==1];
z_obs=z[miss_indi==0];
z_miss=z[miss_indi==1];

# observed data is cbind(x, y_miss);

# compelete data analysis;
# sample mean;
BD_mean_vec[cycle]=mean(y);
BD_mean_var_vec[cycle]=var(y)/rowobs;


# linear regression coefficient of z on y, x, and u;
BD_linear=summary(lm(z~0+y+x+u));
BD_slope_vec[cycle]=BD_linear$coef[1,1];
BD_slope_var_vec[cycle]=BD_linear$coef[1,2]^2;


# different missing data methods;

# complete-case analysis;
# mean estimands;
CC_mean=mean(y_miss, na.rm="T");
CC_mean_vec[cycle]=CC_mean;
CC_mean_var_vec[cycle]=var(y_miss, na.rm="T")/obs_no;

CC_linear=summary(lm(z_obs~0+y_obs+x_obs+u_obs));
CC_slope_vec[cycle]=CC_linear$coef[1,1];
CC_slope_var_vec[cycle]=CC_linear$coef[1,2]^2;

# now impute the missing y's to get estimates of the marginal;
y_completed_ROUND=y_completed_PMM_linear=y_completed_PMM_logit=y_miss;

for (i in 1:mi_no)
{
set.seed(i);

# posterior draws based on the MLE approximation;

mlefit=glm(y_obs~cbind(x_obs,u_obs, z_obs), family=binomial(link="logit"));
regcoef=mlefit$coef;  
covmatrix=vcov(mlefit);
regdraw=mvrnorm(n=1, mu=regcoef, Sigma=covmatrix);

beta0_draw=regdraw[1];
beta1_draw=regdraw[2];
# top code the draws;
if (beta1_draw > 20) beta1_draw=20;
beta2_draw=regdraw[3];
beta3_draw=regdraw[4];

coef_s1x_draw=c(beta0_draw, beta1_draw, beta2_draw, beta3_draw);


# PMM imputation based on a logistic model;
y_imputed_PMM_logit=rep(NA, mis_no);
obs_pred_mean=cbind(cbind(1, x_obs, u_obs, z_obs)%*%regcoef,y_obs);
miss_pred_draw=cbind(1, x_miss, u_miss, z_miss)%*%coef_s1x_draw;

# choose number of donors as 5;
donor_no=5;
for(k in 1:mis_no)
{
donor_matrix=obs_pred_mean;
donor_matrix[,1]=abs(miss_pred_draw[k]-obs_pred_mean[,1]);
donor_matrix_sorted=donor_matrix[order(donor_matrix[,1]),];
y_imputed_PMM_logit[k]=sample(donor_matrix_sorted[1:donor_no,2],1);
}


# Rounding imputation;
# based on a linear normal model;
y_imputed_ROUND=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x,u,z));
y_imputed_ROUND[y_imputed_ROUND >=0.5]=1;
y_imputed_ROUND[y_imputed_ROUND < 0.5]=0;


# PMM imputation based on a linear normal model;

y_imputed_PMM_linear=mice.impute.pmm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x,u,z));


y_completed_ROUND[miss_seq]=y_imputed_ROUND;
y_completed_PMM_logit[miss_seq]=y_imputed_PMM_logit;
y_completed_PMM_linear[miss_seq]=y_imputed_PMM_linear;

# Assign estimates;
# Rounding
# Mean estimates;
ROUND_mean_mat[cycle,i]=mean(y_completed_ROUND);
ROUND_mean_var_mat[cycle,i]=var(y_completed_ROUND)/rowobs;

# Slope;
ROUND_linear=summary(lm(z~0+y_completed_ROUND+x+u));
ROUND_slope_mat[cycle,i]=ROUND_linear$coef[1,1];
ROUND_slope_var_mat[cycle,i]=ROUND_linear$coef[1,2]^2;

# PMM based on a logistic model;
# Mean estimates;
PMM_logit_mean_mat[cycle,i]=mean(y_completed_PMM_logit);
PMM_logit_mean_var_mat[cycle,i]=var(y_completed_PMM_logit)/rowobs;

# Slope;
PMM_logit_linear=summary(lm(z~0+y_completed_PMM_logit+x+u));
PMM_logit_slope_mat[cycle,i]=PMM_logit_linear$coef[1,1];
PMM_logit_slope_var_mat[cycle,i]=PMM_logit_linear$coef[1,2]^2;

# PMM based on a linear model;
# Mean estimates;
PMM_linear_mean_mat[cycle,i]=mean(y_completed_PMM_linear);
PMM_linear_mean_var_mat[cycle,i]=var(y_completed_PMM_linear)/rowobs;

# Slope
PMM_linear_linear=summary(lm(z~0+y_completed_PMM_linear+x+u));
PMM_linear_slope_mat[cycle,i]=PMM_linear_linear$coef[1,1];
PMM_linear_slope_var_mat[cycle,i]=PMM_linear_linear$coef[1,2]^2;

}

cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

#############################################################################
# Simulation results
# Table 5.5
# The mean of y

true_mean=mean(BD_mean_vec);
true_mean;

# before-deletion;
BD=performance(BD_mean_vec, BD_mean_var_vec, true_mean);
BD;

# complete-case analysis;
CC=performance(CC_mean_vec, CC_mean_var_vec, true_mean);
CC;

# Multiple imputation analysis;

# Imputation based on rounding
ROUND_MI=mi_performance(ROUND_mean_mat, ROUND_mean_var_mat, true_mean, rowobs);
ROUND_MI;

# PMM Imputation based on a linear model;
PMM_linear_MI=mi_performance(PMM_linear_mean_mat, PMM_linear_mean_var_mat, true_mean, rowobs);
PMM_linear_MI;


# PMM Imputation based on a logistic model;
PMM_logit_MI=mi_performance(PMM_logit_mean_mat, PMM_logit_mean_var_mat, true_mean, rowobs);
PMM_logit_MI;



##########################################################################
# slope estimand;

true_slope=mean(BD_slope_vec);
true_slope;

# before-deletion;
BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;

# complete-case analysis;
CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;

# Multiple imputation analysis;

# Round imputation;
ROUND_MI=mi_performance(ROUND_slope_mat, ROUND_slope_var_mat, true_slope, rowobs);
ROUND_MI;

# PMM imputation based on a linear normal model;
PMM_linear_MI=mi_performance(PMM_linear_slope_mat, PMM_linear_slope_var_mat, true_slope, rowobs);
PMM_linear_MI;


# PMM imputation based on a logit model;
PMM_logit_MI=mi_performance(PMM_logit_slope_mat, PMM_logit_slope_var_mat, true_slope, rowobs);
PMM_logit_MI;

