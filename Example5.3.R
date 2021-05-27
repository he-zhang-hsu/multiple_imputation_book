# Example 5.3

rm(list=ls());
# library(car);
# library(mvtnorm);
library(mice);
library(norm);
library(MASS);
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

# function parameter_draw() draws regression parameter
# from a linear normal model;
parameter_draw=function(data)
{
 # obtain all the computational components;
# data=cbind(x_obs, y_obs);
 p=ncol(data);
 n_total=nrow(data);
 x=cbind(rep(1,n_total), data[,1:(p-1)]);
 y=data[,p];
 # perform least-square calculation;
 beta_hat=solve(t(x)%*%x)%*%t(x)%*%y;
 s_square=1/(n_total-p)*crossprod(y-x%*%beta_hat);
 inv_v_beta=solve(t(x)%*%x);
 
 # draw sigma_square;
 # set.seed(19760421);
 sigma_square=1/(rgamma(1, shape=n_total/2-p/2, scale=2/((n_total-p)*s_square)));
 
 # draw beta;
 beta=mvrnorm(n=1, mu=beta_hat, Sigma=sigma_square*inv_v_beta);

 parameter=c(beta_hat, beta,sigma_square);
 parameter;
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

# error variance;
var_error=1;


cycle_no=1000;
mi_no=50;

# Vectors and matrices holding the parameter estimates;
# Before-deletion analysis;
BD_mean_vec=BD_mean_var_vec=BD_slope_vec=BD_slope_var_vec=BD_xyslope_vec=BD_xyslope_var_vec=rep(NA, cycle_no);

# Complete-case analysis;
CC_mean_vec=CC_mean_var_vec=CC_slope_vec=CC_slope_var_vec=CC_xyslope_vec=CC_xyslope_var_vec=rep(NA, cycle_no);

# multiple imputation analysis;
# MI for normal linear model imputation ignoring the x^2;
MI_mean_mat=MI_mean_var_mat=MI_slope_mat=MI_slope_var_mat=MI_xyslope_mat=MI_xyslope_var_mat=matrix(NA, cycle_no, mi_no);


# IM_MI for the kernel regression imputation
IM_MI_mean_mat=IM_MI_mean_var_mat=IM_MI_slope_mat=IM_MI_slope_var_mat=IM_MI_xyslope_mat=IM_MI_xyslope_var_mat=matrix(NA, cycle_no, mi_no);

# begin the simulation;

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# normal error distribution;
error=sqrt(var_error)*rnorm(rowobs);

# generate the y variable;
y=beta0+beta1*x+beta2*x^2+error;

# set up the missing data as MCAR on x;
miss_indi=runif(n=rowobs)<0.30;

# set up the missing data as MAR on x;
# alpha0=-2;
# alpha1=1;
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x)/(1+exp(alpha0+alpha1*x)));


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
# the predictor includes x and x^2
BD=summary(lm(y~x+I(x*x)));

# linear term;
BD_slope_vec[cycle]=BD$coeff[2,1];
BD_slope_var_vec[cycle]=BD$coeff[2,2]^2;

# quardratic term;

BD_xyslope_vec[cycle]=BD$coeff[3,1];
BD_xyslope_var_vec[cycle]=BD$coeff[3,2]^2;


# different missing data methods;

# complete-case analysis;

# marginal mean;
CC_mean_vec[cycle]=mean(y_obs);
CC_mean_var_vec[cycle]=var(y_obs)/obs_no;


# regression analysis;
CC=summary(lm(y_obs~x_obs+I(x_obs*x_obs)));

# linear term;
CC_slope_vec[cycle]=CC$coeff[2,1];
CC_slope_var_vec[cycle]=CC$coeff[2,2]^2;

# quardratic term;

CC_xyslope_vec[cycle]=CC$coeff[3,1];
CC_xyslope_var_vec[cycle]=CC$coeff[3,2]^2;

# now impute the missing y's;
y_completed=y_completed_IM=y_miss;

# h is the bandwidth of the kernel imputation;
# termed as lambda in the book
# h=0.5;
h=0.1;
# h=0.2;
# h=0.05;
# h=0.01;

for (i in 1:mi_no)
{
set.seed(i);

# linear normal model imputation;
# ignoring x^2;
y_imputed=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);


y_completed[miss_seq]=y_imputed;

# completed-data analysis;

# mean of y;
MI_mean_mat[cycle,i]=mean(y_completed);
MI_mean_var_mat[cycle,i]=var(y_completed)/rowobs;

# regression analysis;

IMP=summary(lm(y_completed~x+I(x*x)));

# linear term; 
MI_slope_mat[cycle,i]=IMP$coeff[2,1];
MI_slope_var_mat[cycle,i]=IMP$coeff[2,2]^2;

# quadratic term;
MI_xyslope_mat[cycle,i]=IMP$coeff[3,1];
MI_xyslope_var_mat[cycle,i]=IMP$coeff[3,2]^2;


# kernel resampling method;

# bootstrap the observed data;
boot_index=sample(seq(1,obs_no,1), size=obs_no, replace=T);
x_obs_boot=x_obs[boot_index];
y_obs_boot=y_obs[boot_index];
# create the kernel matrix;
W_matrix=matrix(NA, mis_no, obs_no);
y_imputed_IM=rep(NA, mis_no);
for (l in 1:mis_no)
{
 K=exp(-1/2*((x_mis[l]-x_obs_boot)/h)^2)/h;
 if (sum(K)==0){W_matrix[l,]=1/obs_no;}
 else{ 
 W_matrix[l,]=K/sum(K);}
 y_imputed_IM[l]=sample(y_obs_boot, size = 1, replace = T, prob = W_matrix[l,]);
}

y_completed_IM[miss_seq]=y_imputed_IM;

# completed-data analysis;

# mean of y;
IM_MI_mean_mat[cycle,i]=mean(y_completed_IM);
IM_MI_mean_var_mat[cycle,i]=var(y_completed_IM)/rowobs;

# regression analysis;

IM_IMP=summary(lm(y_completed_IM~x+I(x*x)));

# linear term;
IM_MI_slope_mat[cycle,i]=IM_IMP$coeff[2,1];
IM_MI_slope_var_mat[cycle,i]=IM_IMP$coeff[2,2]^2;

# quardratic term;
IM_MI_xyslope_mat[cycle,i]=IM_IMP$coeff[3,1];
IM_MI_xyslope_var_mat[cycle,i]=IM_IMP$coeff[3,2]^2;



}

cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

#############################################################################
# Simulation results;
# Table 5.2
# For the mean estimand;

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

# kernel imputation;
IM_MI=mi_performance(IM_MI_mean_mat, IM_MI_mean_var_mat, true_mean, rowobs);
IM_MI;


##########################################################################
# linear slope estimand;

true_slope=mean(BD_slope_vec);
true_slope;

# before-deletion;
BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;

# complete-case analysis;
CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;

# Multiple imputation;
# normal imputation;
MI=mi_performance(MI_slope_mat, MI_slope_var_mat, true_slope, rowobs);
MI;

# kernel imputation;
IM_MI=mi_performance(IM_MI_slope_mat, IM_MI_slope_var_mat, true_slope, rowobs);
IM_MI;


#######################################################################
# quardratic slope;

true_xyslope=mean(BD_xyslope_vec);
true_xyslope;

# before-deletion; 
BD=performance(BD_xyslope_vec, BD_xyslope_var_vec, true_xyslope);
BD;

# complete-case analysis;

CC=performance(CC_xyslope_vec, CC_xyslope_var_vec, true_xyslope);
CC;

# multiple imputation estimates;
# normal imputation
MI=mi_performance(MI_xyslope_mat, MI_xyslope_var_mat, true_xyslope, rowobs);
MI;

# kernel imputation;
IM_MI=mi_performance(IM_MI_xyslope_mat, IM_MI_xyslope_var_mat, true_xyslope, rowobs);
IM_MI;



########################################################
# plots
# Fig. 5.5
# for h=0.1 and last simulation replicate

y_plot_obs=y;
y_plot_obs[miss_indi==1]=NA;

# for normal imputation;
y_plot_imputed=y_completed;
y_plot_imputed[miss_indi==0]=NA;
y_plot_impute=cbind(y_plot_obs, y_plot_imputed);

matplot(x[1:250], y_plot_impute[1:250,], pch=c(1,17), xlab="x", ylab="y", main="Normal Imputation");

# for kernel imputation;

y_plot_imputed=y_completed_IM;
y_plot_imputed[miss_indi==0]=NA;
y_plot_impute=cbind(y_plot_obs, y_plot_imputed);

matplot(x[1:250], y_plot_impute[1:250,], pch=c(1,17), xlab="x", ylab="y", main="Kernel Imputation");

# calculate the W matrix;

W_matrix_mean=matrix(NA, mis_no, obs_no);
# y_imputed_IM=rep(NA, mis_no);
for (l in 1:mis_no)
{
 K=exp(-1/2*((x_mis[l]-x_obs)/h)^2)/h;
 if (sum(K)==0){W_matrix_mean[l,]=1/obs_no;}
 else{ 
 W_matrix_mean[l,]=K/sum(K);}
# y_imputed_IM[l]=sample(y_obs_boot, size = 1, replace = T, prob = W_matrix[l,]);
}

# x is -.2059;
x_mis[1];

# Fig. 5.4
plot(x_obs, W_matrix_mean[1,], ylab="Sampling Probability");
abline(v=x_mis[1], col="red");

plot(y_obs, W_matrix_mean[1,], ylab="Sampling Probability");

prob=W_matrix_mean[1,];

prob_data=cbind(y_obs, prob);
new_prob_data=prob_data[order(prob),];


##############################
# plot the PMM imputation;
# Fig. 5.6  

# calculate the distance matrix;
D_matrix_mean=matrix(NA, mis_no, obs_no);

parameter=parameter_draw(cbind(x_obs, y_obs));
beta_hat=parameter[1:2];
beta_draw=parameter[3:4];
y_impute_pmm=cbind(1,x_mis)%*%beta_draw;
y_obs_pred=cbind(1,x_obs)%*%beta_hat;

for (i in 1:mis_no)
{
 for (j in 1:obs_no)
 {
  D_matrix_mean[i,j]=abs(y_impute_pmm[i]-y_obs_pred[j]);
 }
}


plot(x_obs, D_matrix_mean[1,], ylab="distance");
abline(v=x_mis[1], col="red");


plot(y_obs, D_matrix_mean[1,], ylab="Distance Function");

distance=D_matrix_mean[1,];

plot(y_obs[distance <.2], distance[distance<.2]);

# mean(distance[distance<.2]);
# distance_data=cbind(y_obs, distance);
# new_distance_data=distance_data[order(distance),];

# abline(v=x_mis[1], col="red");



