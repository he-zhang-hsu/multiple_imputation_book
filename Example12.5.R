# Example 12.5
rm(list=ls());
library(MASS);

# write a function of doing linear regression imputation;
# This can also be done by calling R MICE;

reg_linear_impute=function(missing_data)
{
 # obtain all the computational components;
 p=ncol(missing_data);
 n_total=nrow(missing_data);
 x=cbind(rep(1,n_total), missing_data[,1:(p-1)]);
 y=missing_data[,p];
 y_obs=y[!is.na(y)];

 obs_no=length(y_obs);
 x_obs=x[!is.na(y),];
 y_mis=y[is.na(y)];
 x_mis=x[is.na(y),];
 
 mis_no=length(y_mis);
 
 # perform least-square calculation;
 beta_hat_obs=solve(t(x_obs)%*%x_obs)%*%t(x_obs)%*%y_obs;
 s_square_obs=1/(obs_no-p)*crossprod(y_obs-x_obs%*%beta_hat_obs);
 inv_v_beta_obs=solve(t(x_obs)%*%x_obs);
 
 # draw sigma_square;
 # set.seed(197789);
 sigma_square=1/(rgamma(1, shape=obs_no/2-p/2, scale=2/((obs_no-p)*s_square_obs)));
 
 # draw beta;
 beta=mvrnorm(n=1, mu=beta_hat_obs, Sigma=sigma_square*inv_v_beta_obs);

 # impute data;
 y_impute=x_mis%*%beta+sqrt(sigma_square)*rnorm(mis_no);

 impute_data=missing_data;
 impute_data[is.na(y),p]=y_impute;

 # return impute_data;
 impute_data;
}


# function reg_quad_impute impute data from normal quadratic regression model;

reg_quad_impute=function(missing_data)
{
 # obtain all the computational components;
 # the covariate are 1, x, and x^2;
 p=ncol(missing_data)+1;
 n_total=nrow(missing_data);
 x=cbind(rep(1,n_total), missing_data[,1:(p-2)], missing_data[,1:(p-2)]^2);
 y=missing_data[,p-1];
 y_obs=y[!is.na(y)];

 obs_no=length(y_obs);
 x_obs=x[!is.na(y),];
 y_mis=y[is.na(y)];
 x_mis=x[is.na(y),];
 
 mis_no=length(y_mis);
 
 # perform least-square calculation;
 beta_hat_obs=solve(t(x_obs)%*%x_obs)%*%t(x_obs)%*%y_obs;
 s_square_obs=1/(obs_no-p)*crossprod(y_obs-x_obs%*%beta_hat_obs);
 inv_v_beta_obs=solve(t(x_obs)%*%x_obs);
 
 # draw sigma_square;
 # set.seed(197789);
 sigma_square=1/(rgamma(1, shape=obs_no/2-p/2, scale=2/((obs_no-p)*s_square_obs)));
 
 # draw beta;
 beta=mvrnorm(n=1, mu=beta_hat_obs, Sigma=sigma_square*inv_v_beta_obs);

 # impute data;
 y_impute=x_mis%*%beta+sqrt(sigma_square)*rnorm(mis_no);

 impute_data=missing_data;
 impute_data[is.na(y),p-1]=y_impute;

 # return impute_data;
 impute_data;
}

# function emp_quantile gives the empirical quantiles of x in density y;
emp_quantile=function(x, y)
{
 n=length(y);
 y_sort=sort(y);
 if (x < y_sort[1])
 {
  p_x=0;
 }
 if (x == y_sort[1])
 {
  p_x=1/n;
 }
 if (x >= y_sort[n])
 {
  p_x=1;
 }
 if (x > y_sort[1] & x < y_sort[n])
 {
 for (i in 1:(n-1))
 {
 if (x > y_sort[i] & x < y_sort[i+1])
 {
  p_x=i/n+(x-y_sort[i])/(y_sort[i+1]-y_sort[i])*1/n;
 }
 }
 }
 p_x;
}



# Create incomplete data set;
# Assume a sinlge X is fully observed;
# and a single Y is fully observed;
# Duplicate X but leave the Y as missing;
# The model of observed Y on X is the key assumption;
# Assume the regression model is a nonlinear second order;

# create complete data before deletion;
# n_ori is the number of complete-data sample size;
n_ori=1000;
# n_ori=50;

# simu_no is the number of simulations;
# simu_no is 600 for n_ori=50; simu_no is 300 for n_ori=1000;
simu_no=300;
# simu_no=600;

# ppp_simu is the matrix holding estimated P-value from simulation-based approach;

# check 17 statistics;
ppp_simu=matrix(0, nrow=simu_no, ncol=17);
ppp_mi=array(0, c(simu_no, 3, 17));

mi_estimates=matrix(0, nrow=simu_no, ncol=17);
rep_estimates=matrix(0, nrow=simu_no, ncol=17);
 
# missrate_vec holds the missingness rates estimated from each replicate;
missrate_vec=rep(0, simu_no);

for (simu in 1:simu_no)
{

# x is a normal(0,1);
set.seed(simu);

# x could be uniform distribution; set is as (-3, 3);
x=6*runif(n_ori)-3;

# y|x is a second order polynomial regression;
# regression coefficient
beta_0=1;
beta_1=1;
beta_2=0.5;
epsilon=rnorm(n_ori);
y=beta_0+beta_1*x+beta_2*x^2+epsilon;
ori_data=cbind(x,y);




# y=beta_0+beta_1*x^2+beta_2*sin(x)+epsilon;

# create missing data on Y, assuming MCAR and the missingess proportion is r;
# r=0.8;
# ori_data_missing=ori_data;
# ori_data_missing[runif(n_ori)<r,2]=NA;

# create missing data on Y, assuming MAR: the logit of missingness prob dependes on x;
# alphar=-1.30 corresponds to missingness proportion is 0.2;
# alphar=-0.40 corresponds to missingness proportion is 0.4;
# alphar=0.40 corresponds to missingness proportion is 0.6;
# alphar=1.30 corrssponds to missingness proportion is 0.8;

alphar=-1.30; 
betar=0.5;

miss=rep(1,n_ori);

while(n_ori-sum(miss) <4)
{
miss=rbinom(n=n_ori,size=1,p=1/(1+exp(-(alphar+betar*x))));
missrate_vec[simu]=mean(miss);

ori_data_missing=ori_data;
ori_data_missing[miss==1,2]=NA;
}

# par(mfrow=c(1,1));
# plot(x=ori_data_missing[,1], y=ori_data_missing[,2]);

# composite the concatenated data
concat_data=rbind(ori_data_missing, cbind(x, rep(NA, n_ori)));

# compare the distribution of observed Y and imputed Y;
# obtain all the computational components;
# n_total=nrow(missing_data);
# x=missing_data[,1];
# y=missing_data[,2];
# y_obs=y[!is.na(y)];
# x_obs=x[!is.na(y)];
# y_mis=y[is.na(y)];
# x_mis=x[is.na(y)];
# sampleindivec=seq(1, n_total);
# missindivec=sampleindivec[is.na(y)];


# fit the observed data with a model;
# in this case, the quardratic model;

# quard_yobs_xobs=lm(y_obs ~ cbind(x_obs,x_obs^2)); 
# summary(quard_yobs_xobs);

# linear_yobs_xobs=lm(y_obs ~ x_obs);
# summary(linear_yobs_xobs);




# simulate the replicates of posterior predictions;
# for each replicate, fit the quardratic model and store the regression parameters;

# rep_no is the number of replicates used in the simulation-based approach;
rep_no=5000;

# stat_yimpute_mat holds statistics from fitting the completed data;
stat_yimpute_mat=matrix(0, nrow=rep_no, ncol=17);

# stat_yrep_mat holds statistics from fitting the replicated data;
stat_yrep_mat=matrix(0, nrow=rep_no, ncol=17);

for (i in 1:rep_no)
{
# using the linear imputation model
# for data with quadratic model;
 impute_data=reg_linear_impute(concat_data);
 y_impute=impute_data[1:n_ori,2];
 y_rep=impute_data[(n_ori+1):(2*n_ori),2];

# quadratic regression of y on x;
quad_yimpute_x=lm(y_impute ~ cbind(x, x^2));
quad_yrep_x=lm(y_rep ~ cbind(x, x^2));

# quadratic regression of x on y;
quad_x_yimpute=lm(x ~ cbind(y_impute, y_impute^2));
quad_x_yrep=lm(x ~ cbind(y_rep, y_rep^2));

# linear regression of y on x;
linear_yimpute_x=lm(y_impute ~ x);
linear_yrep_x=lm(y_rep ~ x);

# linear regression of x on y;
linear_x_yimpute=lm(x ~ y_impute);
linear_x_yrep=lm(x ~ y_rep);

stat_yimpute_mat[i,]=c(mean(y_impute), var(y_impute), quantile(y_impute, probs=c(0.05, 0.25, 0.50, 0.75, 0.95)),
linear_yimpute_x$coef, linear_x_yimpute$coef, quad_yimpute_x$coef, 
quad_x_yimpute$coef);

stat_yrep_mat[i,]=c(mean(y_rep), var(y_rep), quantile(y_rep, probs=c(0.05, 0.25, 0.50, 0.75, 0.95)),
linear_yrep_x$coef, linear_x_yrep$coef, quad_yrep_x$coef, 
quad_x_yrep$coef);


}


# plot the histogram of coef_yimpute_mat and coef_yrep_mat;
# par(mfrow=c(2,3));
# hist(coef_yimpute_mat[,1], main="the histogram of imputed intercept estimate");
# hist(coef_yimpute_mat[,2], main="the histogram of imputed slope estiamte");
# hist(coef_yimpute_mat[,3], main="the histogram of imputed quadratic order estimate");
# hist(coef_yimpute_mat[,4], main="the histogram of log variance estimate");
# hist(coef_yrep_mat[,1], main="the histogram of replicated intercept estimate");
# hist(coef_yrep_mat[,2], main="the histogram of replicated slope estiamte");
# hist(coef_yrep_mat[,3], main="the histogram of replicated quadratic order estimate");
# hist(coef_yrep_mat[,4], main="the histogram of replicated log variance estiamte");



# the difference between two sets of parameter estimates;

stat_diff_mat=stat_yrep_mat-stat_yimpute_mat;
 
# plot the histogram plots of coef_diff_mat;

# par(mfrow=c(2,2));
# hist(coef_diff_mat[,1], main="the histogram of intercept estimate difference");
# hist(coef_diff_mat[,2], main="the histogram of slope estiamte difference");
# hist(coef_diff_mat[,3], main="the histogram of quadratic order estimate difference");
# hist(coef_diff_mat[,4], main="the histogram of log variance estimate difference");

# par(mfrow=c(2,2));
# qqnorm(coef_diff_mat[,1], main="the qq-plot of intercept estimate difference");
# qqline(coef_diff_mat[,1]);
# qqnorm(coef_diff_mat[,2], main="the qq-plot of slope estiamte difference");
# qqline(coef_diff_mat[,2]);
# qqnorm(coef_diff_mat[,3], main="the qq-plot of quadratic order estimate difference");
# qqline(coef_diff_mat[,3]);
# hist(coef_diff_mat[,4], main="the histogram of log variance estimate difference");


# hist(coef_diff_mat[,4]/coef_yimpute_mat[,4], main="the histogram of intercept sd estimate relative difference");
# hist(coef_diff_mat[,5]/coef_yimpute_mat[,5], main="the histogram of slope order sd estimate relative difference");
# hist(coef_diff_mat[,6]/coef_yimpute_mat[,6], main="the histogram of quadratic order sd estimate relative difference");


# quantile(coef_diff_mat[,1], c(0.05, 0.95));
# quantile(coef_diff_mat[,2], c(0.05, 0.95));
# quantile(coef_diff_mat[,3], c(0.05, 0.95));
# quantile(coef_diff_mat[,4], c(0.05, 0.95));
# quantile(coef_diff_mat[,4]/coef_yimpute_mat[,4], c(0.05, 0.95));
# quantile(coef_diff_mat[,5]/coef_yimpute_mat[,5], c(0.05, 0.95));
# quantile(coef_diff_mat[,6]/coef_yimpute_mat[,6], c(0.05, 0.95));

# Calcuate the posterior predictive p-values;

for (l in 1:17)
{ 
ppp_simu[simu,l]=emp_quantile(x=0, y=stat_diff_mat[,l]);
}

# The average of the estimates from imputations and replicates;

mi_estimates[simu,]=colMeans(stat_yimpute_mat);
rep_estimates[simu,]=colMeans(stat_yrep_mat);

cat("the simulation is", simu, "\n");


}

# the average missingness rate is
mean(missrate_vec);

# Table 12.3
# the average P-values from the simulation-based approach is
colMeans(ppp_simu);

# Estimates from completed-data and their replicates;
colMeans(mi_estimates);
colMeans(rep_estimates);

