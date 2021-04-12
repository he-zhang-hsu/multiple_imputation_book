################## data generate ##############
# Program: boxcox_impute_expolore_simulation; ############
# generate data, and try out different missing data methods;
# run some simulations to assess the biasedness of various methods;

library(car);
library(mvtnorm);
library(mice);
library(norm);
library(HI);
options(digits=4);
rm(list=ls());


bootstrap=function(n)
{
 rowobs=n;  
 donor_vector=ceiling(runif(rowobs)*rowobs);
 donor_vector;
}


bay_bootstrap=function(n)
{

# First produce n-1 uniform random numbers and sort them in increasing order 
# the minimum number is zero and the max is 1 ***/
  rowobs=n;
  x=c(0,sort(runif(rowobs-1)),1);

# then for any uniform number y[row], return index j so that y[j-1] < z <= y[j] #
y=runif(rowobs);

donor_vector=rep(0,rowobs);  

for (row in 1:rowobs)
 {
  begin_index=1;
  end_index=rowobs+1;

while ((end_index-begin_index)>1)
{
  middle_index=floor((begin_index+end_index)/2);
  if (y[row] <=x[middle_index]) end_index=middle_index
  else begin_index=middle_index
}

donor_vector[row]=end_index-1;

}
donor_vector;
}

parameter_draw=function(data)
{
 # obtain all the computational components;
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

 parameter=c(beta,sigma_square);
 parameter;
}

boxcox.like=function(mu,z,lambda,sigma_square,n)
{
 # x is the covariate;
 # beta0 and beta1 are coefficients;
 # sigma_square is the residual variance;
 # n is the sample size;
partial.like=-n/2*log(2*pi)-n/2*log(sigma_square)-1/(2*sigma_square)*crossprod((z^(lambda)-1)/lambda-mu)+(lambda-1)*(1-2/n)*sum(log(z));
# partial.like=-lambda;
 partial.like;
}


boxcox.marginal.like=function(z,lambda,n,covariate)
{
  # n is the sample size;
  beta_lambda_hat=solve(t(covariate)%*%covariate)%*%t(covariate)%*%(z^(lambda)-1)/lambda;
  cross_prod=crossprod((z^(lambda)-1)/lambda-covariate%*%beta_lambda_hat);
partial.like=-(n-2)/2*(log(cross_prod)-2/n*(lambda-1)*sum(log(z)));
# partial.like=-lambda;
 partial.like;
}


# function region.lambda gives the region of lambda
region.lambda=function(mu,z,lambda,sigma_square,n)
{
 (lambda > -3)*(lambda < 3)
}

region.marginal.lambda=function(z,lambda,n,covariate)
{
 (lambda > -3)*(lambda < 3)
}

# function boxcox_cov gives the variance-covariance estimates of the box-cox models;
boxcox_cov=function(x_mat, beta_com, y2, lambda, sigma2, n)
{ 

# constructing the fisher information matrix;
p=length(beta_com);
epsilon=(y2^lambda-1)/lambda-x_mat%*%beta_com;
e_dlambda=1/lambda^2*((log(y2)*lambda-1)*y2^lambda+1);
e_dlambda2=y2^lambda*log(y2)^2/lambda-2/lambda*e_dlambda;


cov_betabeta=t(x_mat)%*%x_mat;
cov_betalambda=-t(x_mat)%*%e_dlambda;
# cov_betasigma2=rep(0,p);
cov_betasigma2=t(x_mat)%*%epsilon/sigma2;
cov_lambdalambda=crossprod(epsilon, e_dlambda2)+crossprod(e_dlambda);
# cov_lambdasigma2=-sum(log(y2));
cov_lambdasigma2=-crossprod(epsilon,e_dlambda)/sigma2;
cov_sigma2sigma2=crossprod(epsilon)/sigma2^2-n/(2*sigma2);
final_cov=matrix(0, nrow=p+2, ncol=p+2);
final_cov[1:p, 1:p]=cov_betabeta;
final_cov[1:p,p+1]=cov_betalambda;
final_cov[1:p,p+2]=cov_betasigma2;
final_cov[p+1,p+1]=cov_lambdalambda;
final_cov[p+1,p+2]=cov_lambdasigma2;
final_cov[p+2,p+2]=cov_sigma2sigma2;
final_cov[p+1,1:p]=final_cov[1:p,p+1];
final_cov[p+2,p+1]=final_cov[p+1,p+2];
para_cov=ginv(final_cov/sigma2);
para_cov;
}


#boxcox.like=function(lambda, y)
#{
# partial.like=-n/2*log(sigma_square)-1/(2*sigma_square)*crossprod((y^(lambda)-1)/lambda-x%*%beta)+(lambda-1)*sum(log(y));
#partial.like=-lambda;
# partial.like;
#}

# function region.lambda gives the region of lambda
#region.lambda=function(lambda, y)
#{
# (lambda > -3)*(lambda < 3)
#}

# lambda_init=1;
# partial.like=-n_obs/2*log(2*pi)-n_obs/2*log(sigma_square_draw)-1/(2*sigma_square_draw)*crossprod((y2_obs^(lambda_init)-1)/lambda_init-cbind(rep(1,n_obs),x_obs)%*%beta_draw)+(lambda_init-1)*sum(log(y2_obs));

# partial.like=-n_obs/2*log(2*pi)-n_obs/2*log(sigma_square_draw)-1/(2*sigma_square_draw)*sum(((y2_obs^(lambda_init)-1)/lambda_init-cbind(rep(1,n_obs),x_obs)%*%beta_draw)^2)+(lambda_init-1)*sum(log(y2_obs));


rowobs=1000;

# set up the random seed;
set.seed(197789);

# generate data;
# distribution of x, the covariate;
# could be normal, or could be other distributions;
# fixed at mean 10, and variance 1;
mu_x=10;
# mu_x=8;
var_x=4;

# regression parameters;
# both beta0 and beta1 are fixed at 
beta0=1;
beta1=1;
# beta1=2;
# beta1=1;
# beta1=1;
# lambda;

lambda=1/3;


# error variance;
# r_square is the variance explained by the model;
r_square=0.5;
var_error=(1-r_square)/r_square * beta1^2*var_x;


cycle_no=1000;
mi_no=30;

# matrices holding the parameter estimates;
# complete-data inferences;
y2_mean_vec=y2_var_vec=cc_mean_vec=cc_mean_var_vec=rep(NA, cycle_no);

# y2_int_mean_vec=y2_int_var_vec=y2_fisher_int_var_vec=cc_int_mean_vec=cc_int_var_vec=rep(NA, cycle_no);

# y2_jacob_int_mean_vec=y2_jacob_int_var_vec=rep(NA, cycle_no);


# y2_slope_mean_vec=y2_slope_var_vec=y2_fisher_slope_var_vec=cc_slope_mean_vec=cc_slope_var_vec=rep(NA, cycle_no);

# y2_jacob_slope_mean_vec=y2_jacob_slope_var_vec=rep(NA, cycle_no);

# y2_cc_fisher_int_var_vec=y2_cc_fisher_slope_var_vec=rep(NA, cycle_no);

marginal_mean_mat=marginal_var_mat=obs_mean_mat=obs_var_mat=guess_mean_mat=guess_var_mat=draw_mean_mat=draw_var_mat=lambdatrue_mean_mat=lambdatrue_var_mat=boot_mean_mat=boot_var_mat=bayes_mean_mat=bayes_var_mat=matrix(NA, cycle_no, mi_no);

# marginal_com_int_mean_mat=marginal_com_int_var_mat=obs_com_int_mean_mat=obs_com_int_var_mat=guess_com_int_mean_mat=guess_com_int_var_mat=draw_com_int_mean_mat=draw_com_int_var_mat=lambdatrue_com_int_mean_mat=lambdatrue_com_int_var_mat=boot_com_int_mean_mat=boot_com_int_var_mat=bayes_com_int_mean_mat=matrix(NA, cycle_no, mi_no);

# marginal_com_fisher_int_var_mat=marginal_com_fisher_slope_var_mat=matrix(NA, cycle_no, mi_no);

# obs_com_fisher_int_var_mat=obs_com_fisher_slope_var_mat=matrix(NA, cycle_no, mi_no);

# guess_com_fisher_int_var_mat=guess_com_fisher_slope_var_mat=matrix(NA, cycle_no, mi_no);

# marginal_com_slope_mean_mat=marginal_com_slope_var_mat=obs_com_slope_mean_mat=obs_com_slope_var_mat=guess_com_slope_mean_mat=guess_com_slope_var_mat=draw_com_slope_mean_mat=draw_com_slope_var_mat=lambdatrue_com_slope_mean_mat=lambdatrue_com_slope_var_mat=boot_com_slope_mean_mat=boot_com_slope_var_mat=bayes_com_slope_mean_mat=matrix(NA, cycle_no, mi_no);

# draw_com_fisher_int_var_mat=draw_com_fisher_slope_var_mat=matrix(NA, cycle_no, mi_no);

# lambdatrue_com_fisher_int_var_mat=lambdatrue_com_fisher_slope_var_mat=matrix(NA, cycle_no, mi_no);

# boot_com_fisher_int_var_mat=boot_com_fisher_slope_var_mat=matrix(NA, cycle_no, mi_no);

# bayes_com_int_var_mat=bayes_com_slope_var_mat=bayes_com_fisher_int_var_mat=bayes_com_fisher_slope_var_mat=matrix(NA, cycle_no, mi_no);

# y2_quantile_mat=cc_quantile_mat=marginal_quantile_mat=obs_quantile_mat=guess_quantile_mat=draw_quantile_mat=matrix(NA, cycle_no, mi_no);
# corr_y2_x_vec=cc_corr_y2_x_vec=marginal_corr_y2_x_vec=obs_corr_y2_x_vec=guess_corr_y2_x_vec=draw_corr_y2_x_vec=rep(NA, cycle_no);
# corr_y2_x_rank_vec=cc_corr_y2_x_rank_vec=marginal_corr_y2_x_rank_vec=obs_corr_y2_x_rank_vec=guess_corr_y2_x_rank_vec=draw_corr_y2_x_rank_vec=rep(NA, cycle_no);

lambda_com_vec=lambda_marginal_vec=lambda_marginal_com_vec=lambda_observed_vec=lambda_obs_com_vec=lambda_guess_vec=lambda_draw_vec=lambda_draw_com_vec=var_lambda_vec=rep(NA, cycle_no);
# beta0_com_vec=beta0_cc_vec=beta0_marginal_vec=beta0_observed_vec=beta0_obs_com_vec=beta0_guess_vec=beta0_guess_com_vec=beta0_draw_vec=beta0_draw_com_vec=rep(NA, cycle_no);
# beta1_com_vec=beta1_cc_vec=beta1_marginal_vec=beta1_observed_vec=beta1_obs_com_vec=beta1_guess_vec=beta1_guess_com_vec=beta1_draw_vec=beta1_draw_com_vec=rep(NA, cycle_no);

# lambda_obs_com_mat=lambda_guess_com_mat=lambda_draw_com_mat=lambda_true_com_mat=lambda_boot_com_mat=lambda_bayes_com_mat=matrix(NA, cycle_no, mi_no);


# y2_correct_int_mean_vec=y2_correct_slope_mean_vec=y2_correct_int_var_vec=y2_correct_slope_var_vec=rep(NA, cycle_no);
# cc_correct_int_mean_vec=cc_correct_slope_mean_vec=cc_correct_int_var_vec=cc_correct_slope_var_vec=rep(NA, cycle_no);
# obs_com_correct_int_mean_mat=obs_com_correct_slope_mean_mat=obs_com_correct_int_var_mat=obs_com_correct_slope_var_mat=matrix(NA, cycle_no, mi_no);
# guess_com_correct_int_mean_mat=guess_com_correct_slope_mean_mat=guess_com_correct_int_var_mat=guess_com_correct_slope_var_mat=matrix(NA, cycle_no, mi_no);

# draw_com_correct_int_mean_mat=draw_com_correct_slope_mean_mat=draw_com_correct_int_var_mat=draw_com_correct_slope_var_mat=matrix(NA, cycle_no, mi_no);
# lambdatrue_com_correct_int_mean_mat=lambdatrue_com_correct_slope_mean_mat=lambdatrue_com_correct_int_var_mat=lambdatrue_com_correct_slope_var_mat=matrix(NA, cycle_no, mi_no);
# boot_com_correct_int_mean_mat=boot_com_correct_slope_mean_mat=boot_com_correct_int_var_mat=boot_com_correct_slope_var_mat=matrix(NA, cycle_no, mi_no);
# bayes_com_correct_int_mean_mat=bayes_com_correct_slope_mean_mat=bayes_com_correct_int_var_mat=bayes_com_correct_slope_var_mat=matrix(NA, cycle_no, mi_no);


obsno_vec=rep(NA, cycle_no);
# obs_corr_mat=guess_corr_mat=draw_corr_mat=lambdatrue_corr_mat=boot_corr_mat=bayes_corr_mat=matrix(NA, cycle_no, mi_no);

# greater than 95%-tile
y2_prop_up_vec=y2_cc_prop_up_vec=rep(NA, cycle_no);
obs_prop_up_mat=marginal_prop_up_mat=guess_prop_up_mat=draw_prop_up_mat=lambdatrue_prop_up_mat=boot_prop_up_mat=bayes_prop_up_mat=matrix(NA, cycle_no, mi_no);

# less than 5%-tile
y2_prop_bottom_vec=y2_cc_prop_bottom_vec=rep(NA, cycle_no);
obs_prop_bottom_mat=marginal_prop_bottom_mat=guess_prop_bottom_mat=draw_prop_bottom_mat=lambdatrue_prop_bottom_mat=boot_prop_bottom_mat=bayes_prop_bottom_mat=matrix(NA, cycle_no, mi_no);



# generate the complete data;
y2_complete_mat=matrix(NA, rowobs, cycle_no);

for (cycle in 1:cycle_no)
{
 set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# error distribution
# could be normal, or could be other distributions;
error=sqrt(var_error)*rnorm(rowobs);

# generate the latent variable;
# take the absolute value to ensure y1 is positive;

y1=abs(beta0+beta1*x+error);

# transform y1 using box-cox transformation;
y2=(1+lambda*y1)^(1/lambda);
y2_complete_mat[,cycle]=y2;

}

population_quantile=quantile(y2_complete_mat, c(0.05, 0.25, 0.5, 0.75, 0.95));

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# error distribution
# could be normal, or could be other distributions;
error=sqrt(var_error)*rnorm(rowobs);

# generate the latent variable;
# take the absolute value to ensure y1 is positive;

y1=abs(beta0+beta1*x+error);

# transform y1 using box-cox transformation;
y2=(1+lambda*y1)^(1/lambda);

# set up the missing data as MAR on x;
# alpha0=-17; # 30% missing;
# alpha0=-16; # 50% missing;
# alpha0=-40; # 50% missing;
 alpha0=36 # 30% missing;
 alpha1=-4;
miss_indi=rbinom(n=rowobs,size=1,p=1/(1+exp(-(alpha0+alpha1*x))));
# MCAR
# miss_indi=runif(n=rowobs)<0.30;
mean(miss_indi);

miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];

# the difference between observed and incomplete data;
# mean(y2[miss_indi==0]);
# mean(y2[miss_indi==1]);

# hist(y2[miss_indi==1]);
# hist(y2[miss_indi==0]);


# miss_indi=runif(n=rowobs)<0.20;
y2_miss=y2;
y2_miss[miss_indi==1]=NA;

# observed data is cbind(x, y2_miss);

# compelete data analysis;
# sample mean;
y2_mean=mean(y2);
y2_var=var(y2)/rowobs;
y2_mean_vec[cycle]=y2_mean;
y2_var_vec[cycle]=y2_var;


# sample quantiles;
# y2_quantile=quantile(y2, c(0.05, 0.25, 0.5, 0.75, 0.95));
# y2_quantile_mat[cycle,]=y2_quantile;

# proportions;
y2_prop_up_vec[cycle]=mean(y2>population_quantile[5]);
y2_prop_bottom_vec[cycle]=mean(y2<population_quantile[1]);


# pearson correlation;
# corr_y2_x=cor(x, y2, method="pearson");
# corr_y2_x_vec[cycle]=corr_y2_x;
# corr_y2_x_rank=cor(x, y2, method="spearman");
# corr_y2_x_rank_vec[cycle]=corr_y2_x_rank;

# regression analysis on the transformed scale
lambda_com=powerTransform(y2~x);
lambda_com_vec[cycle]=lambda_com$lambda;

# tran_y2 is the variable at the transformed scale;
# tran_y2=bcPower(y2, lambda_com$lambda);

# reg_com=lm(tran_y2 ~ x);

# y2_int_mean_vec[cycle]=reg_com$coef[1];
# y2_slope_mean_vec[cycle]=reg_com$coef[2];

# obtain naive variance estimates at the transformed scale;

# y2_int_var_vec[cycle]=vcov(reg_com)[1];
# y2_slope_var_vec[cycle]=vcov(reg_com)[4];


# obtain variance estimates at the transformed scale but account for the variance of lambda;
# sigma2=as.numeric(crossprod(residuals(reg_com))/rowobs);
# para_cov=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_com$coef, y2=y2, lambda=lambda_com$lambda, sigma2=sigma2, n=rowobs);

# y2_fisher_int_var_vec[cycle]=para_cov[1,1];
# y2_fisher_slope_var_vec[cycle]=para_cov[2,2];


# do the regression on the scaled;
# jacob=(exp(mean(log(y2))))^(lambda_com$lambda-1);
# tran_y2_jacob=tran_y2/jacob;

# do the regression at the transformed scale but have the correct guess of lambda;
# lambdacorrect=lambda;
# lambdacorrect=round(lambda_com$lambda);
# y2_correct=bcPower(y2, lambdacorrect);
# reg_com_correct=lm(y2_correct ~ x);
# y2_correct_int_mean_vec[cycle]=reg_com_correct$coef[1];
# y2_correct_slope_mean_vec[cycle]=reg_com_correct$coef[2];

# y2_correct_int_var_vec[cycle]=vcov(reg_com_correct)[1];
# y2_correct_slope_var_vec[cycle]=vcov(reg_com_correct)[4];

# different missing data methods;

# complete-case analysis;
# mean estimands;
cc_mean=mean(y2_miss, na.rm="T");
cc_mean_vec[cycle]=cc_mean;
cc_mean_var_vec[cycle]=var(y2_miss, na.rm="T")/length(y2_miss[!is.na(y2_miss)]);



# cc_quantile=quantile(y2_miss, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm="TRUE");
# cc_quantile_mat[cycle,]=cc_quantile;

# proportions;
y2_cc_prop_up_vec[cycle]=mean(y2_miss>population_quantile[5], na.rm=T);
y2_cc_prop_bottom_vec[cycle]=mean(y2_miss<population_quantile[1], na.rm=T);



# cc_corr_y2_x=cor(x, y2_miss, method="pearson", use="complete.obs");
# cc_corr_y2_x_vec[cycle]=cc_corr_y2_x;

# obsno_vec contains the number of observed cases from each simualtion;
obsno_vec[cycle]=length(y2_miss[!is.na(y2_miss)]);

# cc_corr_y2_x_rank=cor(x, y2_miss, method="spearman", use="complete.obs");
# cc_corr_y2_x_rank_vec[cycle]=cc_corr_y2_x_rank;

# marginal transformation;
lambda_marginal=powerTransform(y2_miss);
lambda_marginal_vec[cycle]=lambda_marginal$lambda;

# transformation conditional on x;
lambda_observed=powerTransform(y2_miss~x);


tran_y2_marginal=bcPower(y2_miss, lambda_observed$lambda);
# tran_y2_marginal=bcPower(y2_miss, lambda_com$lambda);

tran_y2_marginal_true=bcPower(y2_miss, lambda_marginal$lambda);

# reg_cc=lm(tran_y2_marginal ~ x);



# cc_int_mean_vec[cycle]=reg_cc$coef[1];
# cc_slope_mean_vec[cycle]=reg_cc$coef[2];

# obtain naive variance estimates at the transformed scale;

# cc_int_var_vec[cycle]=vcov(reg_cc)[1];
# cc_slope_var_vec[cycle]=vcov(reg_cc)[4];

# obtain variance estimates at the transformed scale but account for the variance of lambda;
# sigma2_cc=as.numeric(crossprod(residuals(reg_cc))/length(y2_miss[!is.na(y2_miss)]));
# para_cov_cc=boxcox_cov(x_mat=cbind(1,x[!is.na(y2_miss)]), beta_com=reg_cc$coef, y2=y2_miss[!is.na(y2_miss)], lambda=lambda_observed$lambda, sigma2=sigma2_cc, n=length(y2_miss[!is.na(y2_miss)]));
# para_cov_cc=boxcox_cov(x_mat=cbind(1,x[!is.na(y2_miss)]), beta_com=reg_cc$coef, y2=y2_miss[!is.na(y2_miss)], lambda=lambda_com$lambda, sigma2=sigma2_cc, n=length(y2_miss[!is.na(y2_miss)]));


# y2_cc_fisher_int_var_vec[cycle]=para_cov_cc[1,1];
# y2_cc_fisher_slope_var_vec[cycle]=para_cov_cc[2,2];

# lambdacorrect=lambda;
# lambdacorrect=round(lambda_observed$lambda);
# y2_miss_correct=bcPower(y2_miss, lambdacorrect);
# reg_cc_correct=lm(y2_miss_correct ~ x);
# cc_correct_int_mean_vec[cycle]=reg_cc_correct$coef[1];
# cc_correct_slope_mean_vec[cycle]=reg_cc_correct$coef[2];

# cc_correct_int_var_vec[cycle]=vcov(reg_cc_correct)[1];
# cc_correct_slope_var_vec[cycle]=vcov(reg_cc_correct)[4];


# imputation after marginal transformation;
y2_marginal_completed=y2_miss;

for (i in 1:mi_no)
{
set.seed(i);
tran_y2_marginal_imputed=mice.impute.norm(tran_y2_marginal_true, ry=as.logical(1-miss_indi), seed=i, x=x);
y2_marginal_completed[miss_seq]=abs((1+lambda_marginal$lambda*tran_y2_marginal_imputed))^(1/lambda_marginal$lambda);

# analysis using completed data;
marginal_mean=mean(y2_marginal_completed);
marginal_var=var(y2_marginal_completed)/rowobs;

marginal_mean_mat[cycle,i]=marginal_mean;
marginal_var_mat[cycle,i]=marginal_var;

marginal_prop_up_mat[cycle,i]=mean(y2_marginal_completed>population_quantile[5]);
marginal_prop_bottom_mat[cycle,i]=mean(y2_marginal_completed<population_quantile[1]);


# marginal transformation after imputation;

# lambda_marginal_com=powerTransform(y2_marginal_completed~x);

# tran_y2_marginal_com=bcPower(y2_marginal_completed, lambda_marginal_com$lambda);

# reg_marginal=lm(tran_y2_marginal_com ~ x);

# marginal_com_int_mean_mat[cycle,i]=reg_marginal$coef[1];
# marginal_com_slope_mean_mat[cycle,i]=reg_marginal$coef[2];

# obtain naive variance estimates at the transformed scale;

# marginal_com_int_var_mat[cycle,i]=vcov(reg_marginal)[1];
# marginal_com_slope_var_mat[cycle,i]=vcov(reg_marginal)[4];

# obtain variance estimates at the transformed scale but account for the variance of lambda;
# sigma2_marginal=as.numeric(crossprod(residuals(reg_marginal))/rowobs);
# para_cov_marginal=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_marginal$coef, y2=y2_marginal_completed, lambda=lambda_marginal_com$lambda, sigma2=sigma2_marginal, n=rowobs);

# marginal_com_fisher_int_var_mat[cycle,i]=para_cov_marginal[1,1];
# marginal_com_fisher_slope_var_mat[cycle,i]=para_cov_marginal[2,2];


}


# marginal_quantile=quantile(y2_marginal_completed, c(0.05, 0.25, 0.5, 0.75, 0.95));
# marginal_quantile_mat[cycle,]=marginal_quantile;

# marginal_corr_y2_x=cor(x, y2_marginal_completed, method="pearson");
# marginal_corr_y2_x_vec[cycle]=marginal_corr_y2_x;

# marginal_corr_y2_x_rank=cor(x, y2_marginal_completed, method="spearman");
# marginal_corr_y2_x_rank_vec[cycle]=marginal_corr_y2_x_rank;


# observed data analysis;
# first obtain the transformation parameter using observed data;
lambda_observed=powerTransform(y2_miss~x);
lambda_observed_vec[cycle]=lambda_observed$lambda;

var_lambda=1/lambda_observed$hessian;
var_lambda_vec[cycle]=var_lambda;

# then obtain the transformed y2;
tran_y2_observed=bcPower(y2_miss, lambda_observed$lambda);

# beta0_observed=lm(tran_y2_observed ~ x)$coef[1];
# beta1_observed=lm(tran_y2_observed ~ x)$coef[2];

# beta0_observed_vec[cycle]=beta0_observed;
# beta1_observed_vec[cycle]=beta1_observed;

# now impute the missing y2's to get estimates of the marginal;
y2_obs_completed=y2_miss;

for (i in 1:mi_no)
{
set.seed(i);
tran_y2_obs_imputed=mice.impute.norm(tran_y2_observed, ry=as.logical(1-miss_indi), seed=i, x=x);
y2_obs_completed[miss_seq]=abs((1+lambda_observed$lambda*tran_y2_obs_imputed))^(1/lambda_observed$lambda);

# marginal means;
obs_mean=mean(y2_obs_completed);
obs_var=var(y2_obs_completed)/rowobs;

obs_mean_mat[cycle,i]=obs_mean;
obs_var_mat[cycle,i]=obs_var;

# proportions;

obs_prop_up_mat[cycle,i]=mean(y2_obs_completed>population_quantile[5]);
obs_prop_bottom_mat[cycle,i]=mean(y2_obs_completed<population_quantile[1]);



# marginal correlation;
# obs_corr_mat[cycle,i]=cor(x, y2_completed, method="pearson");

# lambda_obs_com=powerTransform(y2_completed~x);
# lambda_obs_com_mat[cycle,i]=lambda_obs_com$lambda;
# tran_y2_obs_com=bcPower(y2_completed, lambda_obs_com$lambda);
# tran_y2_obs_com=bcPower(y2_completed, lambda_com$lambda);
# tran_y2_obs_com=bcPower(y2_completed, lambda_observed$lambda);

# reg_obs_com=lm(tran_y2_obs_com ~ x);

# obs_com_int_mean_mat[cycle,i]=reg_obs_com$coef[1];
# obs_com_slope_mean_mat[cycle,i]=reg_obs_com$coef[2];

# obtain naive variance estimates at the transformed scale;

# obs_com_int_var_mat[cycle,i]=vcov(reg_obs_com)[1];
# obs_com_slope_var_mat[cycle,i]=vcov(reg_obs_com)[4];


# obtain variance estimates at the transformed scale but account for the variance of lambda;
# sigma2_obs=as.numeric(crossprod(residuals(reg_obs_com))/rowobs);
# para_cov_obs=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_obs_com$coef, y2=y2_completed, lambda=lambda_obs_com$lambda, sigma2=sigma2_obs, n=rowobs);
# para_cov_obs=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_obs_com$coef, y2=y2_completed, lambda=lambda_com$lambda, sigma2=sigma2_obs, n=rowobs);
# para_cov_obs=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_obs_com$coef, y2=y2_completed, lambda=lambda_observed$lambda, sigma2=sigma2_obs, n=rowobs);



# obs_com_fisher_int_var_mat[cycle,i]=para_cov_obs[1,1];
# obs_com_fisher_slope_var_mat[cycle,i]=para_cov_obs[2,2];


# lambdacorrect=lambda;
# lambdacorrect=round(lambda_obs_com$lambda);
# lambdacorrect=round(lambda_observed$lambda);
# y2_completed_correct=bcPower(y2_completed, lambdacorrect);
# reg_obs_com_correct=lm(y2_completed_correct ~ x);
# obs_com_correct_int_mean_mat[cycle,i]=reg_obs_com_correct$coef[1];
# obs_com_correct_slope_mean_mat[cycle,i]=reg_obs_com_correct$coef[2];

# obs_com_correct_int_var_mat[cycle,i]=vcov(reg_obs_com_correct)[1];
# obs_com_correct_slope_var_mat[cycle,i]=vcov(reg_obs_com_correct)[4];


}

# obs_quantile=quantile(y2_completed, c(0.05, 0.25, 0.5, 0.75, 0.95));
# obs_quantile_mat[cycle,]=obs_quantile;

# obs_corr_y2_x=cor(x, y2_completed, method="pearson");
# obs_corr_y2_x_vec[cycle]=obs_corr_y2_x;

# obs_corr_y2_x_rank=cor(x, y2_completed, method="spearman");
# obs_corr_y2_x_rank_vec[cycle]=obs_corr_y2_x_rank;

# marginal transformation after imputation;
# lambda_obs_com=powerTransform(y2_completed);
# 
# beta0_obs_com=lm(tran_y2_obs_com ~ x)$coef[1];
# beta1_obs_com=lm(tran_y2_obs_com ~ x)$coef[2];

# beta0_obs_com_vec[cycle]=beta0_obs_com;
# beta1_obs_com_vec[cycle]=beta1_obs_com;


# do the proper imputation;
# lambda_draw_vec[cycle]=lambda_draw;

# then obtain the transformed y2;

# beta0_draw=lm(tran_y2_draw ~ x)$coef[1];
# beta1_draw=lm(tran_y2_draw ~ x)$coef[2];

# beta0_draw_vec[cycle]=beta0_draw;
# beta1_draw_vec[cycle]=beta1_draw;

# make a guess of lambda from the incomplete case analysis, still improper imputation

# lambda_guess=round(lambda_observed$lambda);
# lambda_guess_vec[cycle]=lambda_guess;


# then obtain the transformed y2;
# tran_y2_guess=bcPower(y2_miss, lambda_guess);


# now impute the missing y2's to get estimates of the marginal;
y2_guess_completed=y2_miss;

 for (i in 1:mi_no)
 {
set.seed(i);
# tran_y2_guess_imputed=mice.impute.norm(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=x);
# y2_guess_completed[miss_seq]=(1+lambda_guess*tran_y2_guess_imputed)^(1/lambda_guess);
# if (lambda_guess==0)
# y2_guess_completed[miss_seq]=exp(tran_y2_guess_imputed)

# simple linear imputation on the original scale;
tran_y2_guess_imputed=mice.impute.norm(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=x);

# tran_y2_guess_imputed=mice.impute.pmm(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=x);

y2_guess_completed[miss_seq]=tran_y2_guess_imputed;

y2_guess_completed[y2_guess_completed <0]=0;

# marginal means;
guess_mean=mean(y2_guess_completed);
guess_var=var(y2_guess_completed)/rowobs;

guess_mean_mat[cycle,i]=guess_mean;
guess_var_mat[cycle,i]=guess_var;

# proportions;
guess_prop_up_mat[cycle,i]=mean(y2_guess_completed>population_quantile[5]);
guess_prop_bottom_mat[cycle,i]=mean(y2_guess_completed<population_quantile[1]);



# marginal correlation;
# guess_corr_mat[cycle,i]=cor(x, y2_guess_completed, method="pearson");

# lambda_guess_com=powerTransform(y2_guess_completed~x);
# lambda_guess_com_mat[cycle,i]=lambda_guess_com$lambda;
# tran_y2_guess_com=bcPower(y2_guess_completed, lambda_guess_com$lambda);
# tran_y2_guess_com=bcPower(y2_guess_completed, lambda_com$lambda);
# tran_y2_guess_com=bcPower(y2_guess_completed, lambda_observed$lambda);

# reg_guess_com=lm(tran_y2_guess_com ~ x);

# guess_com_int_mean_mat[cycle,i]=reg_guess_com$coef[1];
# guess_com_slope_mean_mat[cycle,i]=reg_guess_com$coef[2];

# obtain naive variance estimates at the transformed scale;

# guess_com_int_var_mat[cycle,i]=vcov(reg_guess_com)[1];
# guess_com_slope_var_mat[cycle,i]=vcov(reg_guess_com)[4];


# obtain variance estimates at the transformed scale but account for the variance of lambda;
# sigma2_guess=as.numeric(crossprod(residuals(reg_guess_com))/rowobs);
# para_cov_guess=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_guess_com$coef, y2=y2_guess_completed, lambda=lambda_guess_com$lambda, sigma2=sigma2_guess, n=rowobs);
# para_cov_guess=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_guess_com$coef, y2=y2_guess_completed, lambda=lambda_obs_com$lambda, sigma2=sigma2_guess, n=rowobs);
# para_cov_guess=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_guess_com$coef, y2=y2_guess_completed, lambda=lambda_observed$lambda, sigma2=sigma2_guess, n=rowobs);



# guess_com_fisher_int_var_mat[cycle,i]=para_cov_guess[1,1];
# guess_com_fisher_slope_var_mat[cycle,i]=para_cov_guess[2,2];

# lambdacorrect=lambda;
# lambdacorrect=round(lambda_guess_com$lambda);
# lambdacorrect=round(lambda_observed$lambda);
# y2_guess_completed_correct=bcPower(y2_guess_completed, lambdacorrect);
# reg_guess_com_correct=lm(y2_guess_completed_correct ~ x);
# guess_com_correct_int_mean_mat[cycle,i]=reg_guess_com_correct$coef[1];
# guess_com_correct_slope_mean_mat[cycle,i]=reg_guess_com_correct$coef[2];

# guess_com_correct_int_var_mat[cycle,i]=vcov(reg_guess_com_correct)[1];
# guess_com_correct_slope_var_mat[cycle,i]=vcov(reg_guess_com_correct)[4];


 }


# now impute the missing y2's to get estimates of the marginal;
y2_draw_completed=y2_miss;

for (i in 1:mi_no)
{
set.seed(i);
lambda_draw=abs(lambda_observed$lambda+sqrt(var_lambda)*rnorm(1));
tran_y2_draw=bcPower(y2_miss, lambda_draw);

tran_y2_draw_imputed=mice.impute.norm(tran_y2_draw, ry=as.logical(1-miss_indi), seed=i, x=x);
y2_draw_completed[miss_seq]=abs((1+lambda_draw*as.vector(tran_y2_draw_imputed)))^(1/lambda_draw);

# marginal means;

draw_mean=mean(y2_draw_completed);
draw_var=var(y2_draw_completed)/rowobs;

draw_mean_mat[cycle,i]=draw_mean;
draw_var_mat[cycle,i]=draw_var;

# proportions;
draw_prop_up_mat[cycle,i]=mean(y2_draw_completed>population_quantile[5]);
draw_prop_bottom_mat[cycle,i]=mean(y2_draw_completed<population_quantile[1]);



# marginal correlation;
# draw_corr_mat[cycle,i]=cor(x, y2_draw_completed, method="pearson");

# lambda_draw_com=powerTransform(y2_draw_completed~x);
# lambda_draw_com_mat[cycle,i]=lambda_draw_com$lambda;

# tran_y2_draw_com=bcPower(y2_draw_completed, lambda_draw_com$lambda);
# tran_y2_draw_com=bcPower(y2_draw_completed, lambda_com$lambda);
# tran_y2_draw_com=bcPower(y2_draw_completed, lambda_observed$lambda);
# tran_y2_draw_com=bcPower(y2_draw_completed, lambda_draw);




# reg_draw_com=lm(tran_y2_draw_com ~ x);

# draw_com_int_mean_mat[cycle,i]=reg_draw_com$coef[1];
# draw_com_slope_mean_mat[cycle,i]=reg_draw_com$coef[2];

# obtain naive variance estimates at the transformed scale;


# draw_com_int_var_mat[cycle,i]=vcov(reg_draw_com)[1];
# draw_com_slope_var_mat[cycle,i]=vcov(reg_draw_com)[4];

# obtain variance estimates at the transformed scale but account for the variance of lambda;
# sigma2_draw=as.numeric(crossprod(residuals(reg_draw_com))/rowobs);
# para_cov_draw=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_draw_com$coef, y2=y2_draw_completed, lambda=lambda_draw_com$lambda, sigma2=sigma2_draw, n=rowobs);
# para_cov_draw=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_draw_com$coef, y2=y2_draw_completed, lambda=lambda_com$lambda, sigma2=sigma2_draw, n=rowobs);
# para_cov_draw=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_draw_com$coef, y2=y2_draw_completed, lambda=lambda_observed$lambda, sigma2=sigma2_draw, n=rowobs);
# para_cov_draw=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_draw_com$coef, y2=y2_draw_completed, lambda=lambda_draw, sigma2=sigma2_draw, n=rowobs);




# draw_com_fisher_int_var_mat[cycle,i]=para_cov_draw[1,1];
# draw_com_fisher_slope_var_mat[cycle,i]=para_cov_draw[2,2];

# lambdacorrect=lambda;
# lambdacorrect=round(lambda_draw_com$lambda);
# lambdacorrect=round(lambda_observed$lambda);


# y2_draw_completed_correct=bcPower(y2_draw_completed, lambdacorrect);
# reg_draw_com_correct=lm(y2_draw_completed_correct ~ x);
# draw_com_correct_int_mean_mat[cycle,i]=reg_draw_com_correct$coef[1];
# draw_com_correct_slope_mean_mat[cycle,i]=reg_draw_com_correct$coef[2];

# draw_com_correct_int_var_mat[cycle,i]=vcov(reg_draw_com_correct)[1];
# draw_com_correct_slope_var_mat[cycle,i]=vcov(reg_draw_com_correct)[4];

}

# draw_quantile=quantile(y2_draw_completed, c(0.05, 0.25, 0.5, 0.75, 0.95));
# draw_quantile_mat[cycle,]=draw_quantile;

# draw_corr_y2_x=cor(x, y2_draw_completed, method="pearson");
# draw_corr_y2_x_vec[cycle]=draw_corr_y2_x;

# draw_corr_y2_x_rank=cor(x, y2_draw_completed, method="spearman");
# draw_corr_y2_x_rank_vec[cycle]=draw_corr_y2_x_rank;

# marginal transformation after imputation;
# lambda_draw_com=powerTransform(y2_draw_completed);
# 
# beta0_draw_com=lm(tran_y2_draw_com ~ x)$coef[1];
# beta1_draw_com=lm(tran_y2_draw_com ~ x)$coef[2];

# beta0_draw_com_vec[cycle]=beta0_draw_com;
# beta1_draw_com_vec[cycle]=beta1_draw_com;

# imputation assuming lambda is known and have the right guess;
y2_lambdatrue_completed=y2_miss;

for (i in 1:mi_no)
{
set.seed(i);
lambda_true=lambda;
tran_y2_lambdatrue=bcPower(y2_miss, lambda_true);

tran_y2_lambdatrue_imputed=mice.impute.norm(tran_y2_lambdatrue, ry=as.logical(1-miss_indi), seed=i, x=x);
y2_lambdatrue_completed[miss_seq]=abs((1+lambda_true*as.vector(tran_y2_lambdatrue_imputed)))^(1/lambda_true);

# marginal means;

lambdatrue_mean=mean(y2_lambdatrue_completed);
lambdatrue_var=var(y2_lambdatrue_completed)/rowobs;

lambdatrue_mean_mat[cycle,i]=lambdatrue_mean;
lambdatrue_var_mat[cycle,i]=lambdatrue_var;

# proportions;

lambdatrue_prop_up_mat[cycle,i]=mean(y2_lambdatrue_completed>population_quantile[5]);
lambdatrue_prop_bottom_mat[cycle,i]=mean(y2_lambdatrue_completed<population_quantile[1]);



# marginal correlation;
# lambdatrue_corr_mat[cycle,i]=cor(x, y2_lambdatrue_completed, method="pearson");


# lambda_true_com=powerTransform(y2_lambdatrue_completed~x);
# lambda_true_com_mat[cycle,i]=lambda_true_com$lambda;


# tran_y2_lambdatrue_com=bcPower(y2_lambdatrue_completed, lambda_true_com$lambda);
# tran_y2_lambdatrue_com=bcPower(y2_lambdatrue_completed, lambda_com$lambda);
# tran_y2_lambdatrue_com=bcPower(y2_lambdatrue_completed, lambda_observed$lambda);

# reg_lambdatrue_com=lm(tran_y2_lambdatrue_com ~ x);

# lambdatrue_com_int_mean_mat[cycle,i]=reg_lambdatrue_com$coef[1];
# lambdatrue_com_slope_mean_mat[cycle,i]=reg_lambdatrue_com$coef[2];

# obtain naive variance estimates at the transformed scale;


# lambdatrue_com_int_var_mat[cycle,i]=vcov(reg_lambdatrue_com)[1];
# lambdatrue_com_slope_var_mat[cycle,i]=vcov(reg_lambdatrue_com)[4];

# obtain variance estimates at the transformed scale but account for the variance of lambda;

# sigma2_lambdatrue=as.numeric(crossprod(residuals(reg_lambdatrue_com))/rowobs);
# para_cov_lambdatrue=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_lambdatrue_com$coef, y2=y2_lambdatrue_completed, lambda=lambda_true_com$lambda, sigma2=sigma2_lambdatrue, n=rowobs);
# para_cov_lambdatrue=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_lambdatrue_com$coef, y2=y2_lambdatrue_completed, lambda=lambda_com$lambda, sigma2=sigma2_lambdatrue, n=rowobs);
# para_cov_lambdatrue=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_lambdatrue_com$coef, y2=y2_lambdatrue_completed, lambda=lambda_observed$lambda, sigma2=sigma2_lambdatrue, n=rowobs);



# lambdatrue_com_fisher_int_var_mat[cycle,i]=para_cov_lambdatrue[1,1];
# lambdatrue_com_fisher_slope_var_mat[cycle,i]=para_cov_lambdatrue[2,2];

# lambdacorrect=lambda;
# lambdacorrect=round(lambda_true_com$lambda);
# lambdacorrect=round(lambda_observed$lambda);


# y2_lambdatrue_completed_correct=bcPower(y2_lambdatrue_completed, lambdacorrect);
# reg_lambdatrue_com_correct=lm(y2_lambdatrue_completed_correct ~ x);
# lambdatrue_com_correct_int_mean_mat[cycle,i]=reg_lambdatrue_com_correct$coef[1];
# lambdatrue_com_correct_slope_mean_mat[cycle,i]=reg_lambdatrue_com_correct$coef[2];

# lambdatrue_com_correct_int_var_mat[cycle,i]=vcov(reg_lambdatrue_com_correct)[1];
# lambdatrue_com_correct_slope_var_mat[cycle,i]=vcov(reg_lambdatrue_com_correct)[4];


}

# imputation using estimated lambda from bootstap sample of observed cases;
y2_obs=y2_miss[miss_indi==0];
x_obs=x[miss_indi==0];
obs_no=length(y2_obs);
y2_boot_completed=y2_miss;

for (i in 1:mi_no)
{
# boot_index=bay_bootstrap(obs_no);
set.seed(i);
boot_index=bootstrap(obs_no);
lambda_boot=abs(powerTransform(y2_obs[boot_index]~x_obs[boot_index])$lambda);

tran_y2_boot=bcPower(y2_miss, lambda_boot);

tran_y2_boot_imputed=mice.impute.norm(tran_y2_boot, ry=as.logical(1-miss_indi), seed=i, x=x);
y2_boot_completed[miss_seq]=abs((1+lambda_boot*as.vector(tran_y2_boot_imputed)))^(1/lambda_boot);

# marginal means;

boot_mean=mean(y2_boot_completed);
boot_var=var(y2_boot_completed)/rowobs;

boot_mean_mat[cycle,i]=boot_mean;
boot_var_mat[cycle,i]=boot_var;

# proportions;
boot_prop_up_mat[cycle,i]=mean(y2_boot_completed>population_quantile[5]);
boot_prop_bottom_mat[cycle,i]=mean(y2_boot_completed<population_quantile[1]);



# marginal correlation;
# boot_corr_mat[cycle,i]=cor(x, y2_boot_completed, method="pearson");


# lambda_boot_com=powerTransform(y2_boot_completed~x);
# lambda_boot_com_mat[cycle,i]=lambda_boot_com$lambda;


# tran_y2_boot_com=bcPower(y2_boot_completed, lambda_boot_com$lambda);
# tran_y2_boot_com=bcPower(y2_boot_completed, lambda_com$lambda);
# tran_y2_boot_com=bcPower(y2_boot_completed, lambda_observed$lambda);
# tran_y2_boot_com=bcPower(y2_boot_completed, lambda_boot$lambda);

# reg_boot_com=lm(tran_y2_boot_com ~ x);



# boot_com_int_mean_mat[cycle,i]=reg_boot_com$coef[1];
# boot_com_slope_mean_mat[cycle,i]=reg_boot_com$coef[2];

# obtain naive variance estimates at the transformed scale;


# boot_com_int_var_mat[cycle,i]=vcov(reg_boot_com)[1];
# boot_com_slope_var_mat[cycle,i]=vcov(reg_boot_com)[4];

# obtain variance estimates at the transformed scale but account for the variance of lambda;

# sigma2_boot=as.numeric(crossprod(residuals(reg_boot_com))/rowobs);


# para_cov_boot=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_boot_com$coef, y2=y2_boot_completed, lambda=lambda_boot_com$lambda, sigma2=sigma2_boot, n=rowobs);
# para_cov_boot=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_boot_com$coef, y2=y2_boot_completed, lambda=lambda_com$lambda, sigma2=sigma2_boot, n=rowobs);
# para_cov_boot=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_boot_com$coef, y2=y2_boot_completed, lambda=lambda_observed$lambda, sigma2=sigma2_boot, n=rowobs);
# para_cov_boot=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_boot_com$coef, y2=y2_boot_completed, lambda=lambda_boot$lambda, sigma2=sigma2_boot, n=rowobs);




# boot_com_fisher_int_var_mat[cycle,i]=para_cov_boot[1,1];
# boot_com_fisher_slope_var_mat[cycle,i]=para_cov_boot[2,2];

# lambdacorrect=lambda;
# lambdacorrect=round(lambda_boot_com$lambda);
# lambdacorrect=round(lambda_observed$lambda);


# y2_boot_completed_correct=bcPower(y2_boot_completed, lambdacorrect);
# reg_boot_com_correct=lm(y2_boot_completed_correct ~ x);
# boot_com_correct_int_mean_mat[cycle,i]=reg_boot_com_correct$coef[1];
# boot_com_correct_slope_mean_mat[cycle,i]=reg_boot_com_correct$coef[2];

# boot_com_correct_int_var_mat[cycle,i]=vcov(reg_boot_com_correct)[1];
# boot_com_correct_slope_var_mat[cycle,i]=vcov(reg_boot_com_correct)[4];


}

# multiple imputation by fully bayesian;
# draw regression parameters, sigma_square, and lambda from observed cases;
y2_obs=y2_miss[miss_indi==0];
x_obs=x[miss_indi==0];

# write out the file;
# write.table(x=cbind(x_obs, y2_obs), file="C:\\Users\\WDQ7\\test.dat", row.names=FALSE, col.names=FALSE, na="NA");


# get initial values of lambda;
# a=powerTransform(y2_obs~x_obs);

# lambda_init=a$lambda;
# var_lambda_init=1/a$hessian;


# tran_y2_obs_init=bcPower(y2_obs, lambda_init);
# reg_com=lm(tran_y2_obs_init ~ x_obs);
# vcov(reg_com)[1];
# vcov(reg_com)[4];



# jacob=(exp(mean(log(y2_obs))))^(lambda_init-1);
# tran_y2_obs_jacob=tran_y2_obs_init/jacob;

# beta_sigmasquare=parameter_draw(cbind(x_obs, tran_y2_obs_jacob));

# beta_draw=beta_sigmasquare[1:2];
# sigma_square_draw=beta_sigmasquare[3];


# beta_draw_mat=matrix(NA, 10000, 2);
# sigma_square_draw_vec=rep(NA, 10000);
# lambda_draw_vec=rep(NA, 10000);
# accept_vec=rep(0, 10000);



# draw regression parameters and variance conditioning on lambda;
# normalize the transformed outcome;

# for (i in 1:10000)
# {

# beta_sigmasquare=parameter_draw(cbind(x_obs, tran_y2_obs_init));

# beta_draw=beta_sigmasquare[1:2];
# sigma_square_draw=beta_sigmasquare[3];


# lambda_init=arms(y.start=1, myldens=boxcox.marginal.like, indFunc=region.marginal.lambda, n.sample=1, covariate=cbind(rep(1,n_obs),x_obs), z=y2_obs, n=n_obs);

# lambda_init=a$lambda+sqrt(var_lambda_init)*rnorm(1);



# using MH algorithm to draw lambda;
# coef_s1x_draw is the parameter draw from last iteration;
# draw coef_s1x_draw from a proposal distribution;
# tune_s1x=0.25;
# tune_lambda=2.4^2;
# lambda_new=lambda_init+sqrt(tune_lambda*var_lambda_init)*rnorm(1);

# evaluate the likelihood ratio;
# ln_lambda_last=boxcox.like(mu=cbind(rep(1,n_obs),x_obs)%*%beta_draw,z=y2_obs,lambda=lambda_init,sigma_square=sigma_square_draw,n=n_obs);
# ln_lambda_new=boxcox.like(mu=cbind(rep(1,n_obs),x_obs)%*%beta_draw,z=y2_obs,lambda=lambda_new,sigma_square=sigma_square_draw,n=n_obs);


# ratio_lambda=exp(ln_lambda_new-ln_lambda_last);

# MH_random_lambda=runif(1);
# if (MH_random_lambda < ratio_lambda) {lambda_init=lambda_new; accept_vec[i]=1};




# tran_y2_obs_init=bcPower(y2_obs, lambda_init);

# beta_draw_mat[i,]=beta_draw;
# sigma_square_draw_vec[i]=sigma_square_draw;
# lambda_draw_vec[i]=lambda_init;

# }

# var_beta0=var(beta_draw_mat[,1]);
# var_beta1=var(beta_draw_mat[,2]);
# var_sigma_square=var(sigma_square_draw_vec);


# now impute the missing y2's to get estimates of the marginal;
y2_bayes_completed=y2_miss;

for (i in 1:mi_no)
{
set.seed(i);
lambda_bayes=abs(arms(y.start=lambda_observed$lambda, myldens=boxcox.marginal.like, indFunc=region.marginal.lambda, n.sample=1, covariate=cbind(rep(1,obs_no),x_obs), z=y2_obs, n=obs_no));

tran_y2_bayes=bcPower(y2_miss, lambda_bayes);

tran_y2_bayes_imputed=mice.impute.norm(tran_y2_bayes, ry=as.logical(1-miss_indi), seed=i, x=x);
y2_bayes_completed[miss_seq]=abs((1+lambda_bayes*as.vector(tran_y2_bayes_imputed)))^(1/lambda_bayes);

# marginal means;

bayes_mean=mean(y2_bayes_completed);
bayes_var=var(y2_bayes_completed)/rowobs;

bayes_mean_mat[cycle,i]=bayes_mean;
bayes_var_mat[cycle,i]=bayes_var;

# proportions;
bayes_prop_up_mat[cycle,i]=mean(y2_bayes_completed>population_quantile[5]);
bayes_prop_bottom_mat[cycle,i]=mean(y2_bayes_completed<population_quantile[1]);



# marginal corrleation;

# bayes_corr_mat[cycle,i]=cor(x, y2_bayes_completed, method="pearson");


# lambda_bayes_com=powerTransform(y2_bayes_completed~x);
# lambda_bayes_com_mat[cycle,i]=lambda_bayes_com$lambda;

# tran_y2_bayes_com=bcPower(y2_bayes_completed, lambda_bayes_com$lambda);
# tran_y2_bayes_com=bcPower(y2_bayes_completed, lambda_com$lambda);
# tran_y2_bayes_com=bcPower(y2_bayes_completed, lambda_observed$lambda);
# tran_y2_bayes_com=bcPower(y2_bayes_completed, lambda_bayes);

# reg_bayes_com=lm(tran_y2_bayes_com ~ x);

# bayes_com_int_mean_mat[cycle,i]=reg_bayes_com$coef[1];
# bayes_com_slope_mean_mat[cycle,i]=reg_bayes_com$coef[2];

# obtain naive variance estimates at the transformed scale;


# bayes_com_int_var_mat[cycle,i]=vcov(reg_bayes_com)[1];
# bayes_com_slope_var_mat[cycle,i]=vcov(reg_bayes_com)[4];

# obtain variance estimates at the transformed scale but account for the variance of lambda;

# sigma2_bayes=as.numeric(crossprod(residuals(reg_bayes_com))/rowobs);
# para_cov_bayes=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_bayes_com$coef, y2=y2_bayes_completed, lambda=lambda_bayes_com$lambda, sigma2=sigma2_bayes, n=rowobs);
# para_cov_bayes=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_bayes_com$coef, y2=y2_bayes_completed, lambda=lambda_com$lambda, sigma2=sigma2_bayes, n=rowobs);
# para_cov_bayes=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_bayes_com$coef, y2=y2_bayes_completed, lambda=lambda_observed$lambda, sigma2=sigma2_bayes, n=rowobs);
# para_cov_bayes=boxcox_cov(x_mat=cbind(1,x), beta_com=reg_bayes_com$coef, y2=y2_bayes_completed, lambda=lambda_bayes, sigma2=sigma2_bayes, n=rowobs);




# bayes_com_fisher_int_var_mat[cycle,i]=para_cov_bayes[1,1];
# bayes_com_fisher_slope_var_mat[cycle,i]=para_cov_bayes[2,2];

# lambdacorrect=lambda;
# lambdacorrect=round(lambda_bayes_com$lambda);
# lambdacorrect=round(lambda_observed$lambda);


# y2_bayes_completed_correct=bcPower(y2_bayes_completed, lambdacorrect);
# reg_bayes_com_correct=lm(y2_bayes_completed_correct ~ x);
# bayes_com_correct_int_mean_mat[cycle,i]=reg_bayes_com_correct$coef[1];
# bayes_com_correct_slope_mean_mat[cycle,i]=reg_bayes_com_correct$coef[2];

# bayes_com_correct_int_var_mat[cycle,i]=vcov(reg_bayes_com_correct)[1];
# bayes_com_correct_slope_var_mat[cycle,i]=vcov(reg_bayes_com_correct)[4];


}

cat("the cycle is", cycle, "\n");



}

# the end of multiple imputation;

# plot the distribution of the data;
# par(mfrow=c(2,2))
# hist(x);
# hist(y2);
# hist(y2_miss);
# plot(x,y2);
# plot(x,y2_miss);

# hist(y1);
# hist(y2);
# save as Fig. 5.1 histograms;
hist(y2, xlab="y", main="Before-deletion");
hist(y2_miss, xlab="y", main="Observed");
hist(y2_guess_completed, xlab="y", main="Normal Imputation");
hist(y2_draw_completed, xlab="y", main="Transformation Imputation");

# same as Fig. 5.1 scatter plot.
test=lm(y2_obs~x_obs);
test2=lowess(y2_obs ~x_obs);
j <- order(x_obs)
plot(x_obs,y2_obs, xlab="x", lty=1, ylab="y");
# plot(test2$x, test2$y);
lines(test2$x,test2$y,col="red",lwd=3)
abline(test, col="green", lty=2, lwd=3);

# check the performance of the estimates;
# population quantity
# mean estimand;

mean(lambda_marginal_vec);
mean(lambda_observed_vec);
mean(var_lambda_vec);

#############################################################################
## Mean estimate;

true=mean(y2_mean_vec);
true;
true_mse=mean((y2_mean_vec-true)^2);
true_mse;

# lower and upper 95% CI;
y2_mean_low95_vec=y2_mean_vec-1.96*sqrt(y2_var_vec);
y2_mean_up95_vec=y2_mean_vec+1.96*sqrt(y2_var_vec);

mean_length=mean(y2_mean_up95_vec-y2_mean_low95_vec);
mean_length;

coverage=(y2_mean_low95_vec < true)*(y2_mean_up95_vec > true);

mean(coverage);

# variance estimates;
mean(sqrt(y2_var_vec));
sqrt(var(y2_mean_vec));

# complete case analysis;
cc_true=mean(cc_mean_vec);
cc_bias=cc_true-true;
cc_bias;

cc_bias/true;

cc_mse=mean((cc_mean_vec-true)^2);
cc_mse;

# lower and upper 95% CI;
cc_mean_low95_vec=cc_mean_vec-1.96*sqrt(cc_mean_var_vec);
cc_mean_up95_vec=cc_mean_vec+1.96*sqrt(cc_mean_var_vec);

mean_cc_length=mean(cc_mean_up95_vec-cc_mean_low95_vec);
mean_cc_length;

cc_coverage=(cc_mean_low95_vec < true)*(cc_mean_up95_vec > true);

mean(cc_coverage);

# variance estimates;
mean(sqrt(cc_mean_var_vec));
sqrt(var(cc_mean_vec));


# multiple imputation estimates;
# marginal transformation imputation;

marginal_mi_mean=rep(NA, cycle_no);
marginal_mi_var=rep(NA, cycle_no);
marginal_mi_df=rep(NA, cycle_no);
marginal_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(marginal_mean_mat[i,], marginal_var_mat[i,], n=rowobs);
marginal_mi_mean[i]=summary$qbar;
marginal_mi_var[i]=summary$t;
marginal_mi_df[i]=summary$df;
marginal_mi_f[i]=summary$f;
}

marginal_mi_true=mean(marginal_mi_mean);
marginal_mi_bias=marginal_mi_true-true;
marginal_mi_bias;
marginal_mi_bias/true;

marginal_mi_mse=mean((marginal_mi_mean-true)^2);
marginal_mi_mse;

# coverage;
marginal_mean_low95_vec=marginal_mi_mean-qt(.975, marginal_mi_df)*sqrt(marginal_mi_var);
marginal_mean_up95_vec=marginal_mi_mean+qt(.975, marginal_mi_df)*sqrt(marginal_mi_var);

mean_marginal_length=mean(marginal_mean_up95_vec-marginal_mean_low95_vec);
mean_marginal_length;


marginal_mi_coverage=(marginal_mean_low95_vec < true)*(marginal_mean_up95_vec > true);
mean(marginal_mi_coverage);

# variance estimates;
mean(sqrt(marginal_mi_var));
sqrt(var(marginal_mi_mean));


# observed data transformation imputation;

obs_mi_mean=rep(NA, cycle_no);
obs_mi_var=rep(NA, cycle_no);
obs_mi_df=rep(NA, cycle_no);
obs_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(obs_mean_mat[i,], obs_var_mat[i,], n=rowobs);
obs_mi_mean[i]=summary$qbar;
obs_mi_var[i]=summary$t;
obs_mi_df[i]=summary$df;
obs_mi_f[i]=summary$f;
}

obs_mi_true=mean(obs_mi_mean);
obs_mi_bias=obs_mi_true-true;
obs_mi_bias;
obs_mi_bias/true;

obs_mi_mse=mean((obs_mi_mean-true)^2);
obs_mi_mse;

# coverage;
obs_mean_low95_vec=obs_mi_mean-qt(.975, obs_mi_df)*sqrt(obs_mi_var);
obs_mean_up95_vec=obs_mi_mean+qt(.975, obs_mi_df)*sqrt(obs_mi_var);

mean_obs_length=mean(obs_mean_up95_vec-obs_mean_low95_vec);
mean_obs_length;


obs_mi_coverage=(obs_mean_low95_vec < true)*(obs_mean_up95_vec > true);
mean(obs_mi_coverage);

# variance estimates;
mean(sqrt(obs_mi_var));
sqrt(var(obs_mi_mean));

# linear normal imputation;

guess_mi_mean=rep(NA, cycle_no);
guess_mi_var=rep(NA, cycle_no);
guess_mi_df=rep(NA, cycle_no);
guess_mi_f=rep(NA, cycle_no);

 for (i in 1:cycle_no)
 {
 summary=pool.scalar(guess_mean_mat[i,], guess_var_mat[i,], n=rowobs);
 guess_mi_mean[i]=summary$qbar;
 guess_mi_var[i]=summary$t;
 guess_mi_df[i]=summary$df;
 guess_mi_f[i]=summary$f;
 }

 guess_mi_true=mean(guess_mi_mean);
 guess_mi_bias=guess_mi_true-true;
 guess_mi_bias;
 guess_mi_bias/true;

 guess_mi_mse=mean((guess_mi_mean-true)^2);
 guess_mi_mse;

# coverage;
 guess_mean_low95_vec=guess_mi_mean-qt(.975, guess_mi_df)*sqrt(guess_mi_var);
 guess_mean_up95_vec=guess_mi_mean+qt(.975, guess_mi_df)*sqrt(guess_mi_var);

 mean_guess_length=mean(guess_mean_up95_vec-guess_mean_low95_vec);
 mean_guess_length;


 guess_mi_coverage=(guess_mean_low95_vec < true)*(guess_mean_up95_vec > true);
 mean(guess_mi_coverage);

# variance estimates;
 mean(sqrt(guess_mi_var));
 sqrt(var(guess_mi_mean));



# imputation using draws;

draw_mi_mean=rep(NA, cycle_no);
draw_mi_var=rep(NA, cycle_no);
draw_mi_df=rep(NA, cycle_no);
draw_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(draw_mean_mat[i,], draw_var_mat[i,], n=rowobs);
draw_mi_mean[i]=summary$qbar;
draw_mi_var[i]=summary$t;
draw_mi_df[i]=summary$df;
draw_mi_f[i]=summary$f;
}

draw_mi_true=mean(draw_mi_mean);
draw_mi_bias=draw_mi_true-true;
draw_mi_bias;
draw_mi_bias/true;

draw_mi_mse=mean((draw_mi_mean-true)^2);
draw_mi_mse;

# coverage;
draw_mean_low95_vec=draw_mi_mean-qt(.975, draw_mi_df)*sqrt(draw_mi_var);
draw_mean_up95_vec=draw_mi_mean+qt(.975, draw_mi_df)*sqrt(draw_mi_var);

mean_draw_length=mean(draw_mean_up95_vec-draw_mean_low95_vec);
mean_draw_length;


draw_mi_coverage=(draw_mean_low95_vec < true)*(draw_mean_up95_vec > true);
mean(draw_mi_coverage);

# variance estimates;
mean(sqrt(draw_mi_var));
sqrt(var(draw_mi_mean));

# imputation using true lambda;

lambdatrue_mi_mean=rep(NA, cycle_no);
lambdatrue_mi_var=rep(NA, cycle_no);
lambdatrue_mi_df=rep(NA, cycle_no);
lambdatrue_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(lambdatrue_mean_mat[i,], lambdatrue_var_mat[i,], n=rowobs);
lambdatrue_mi_mean[i]=summary$qbar;
lambdatrue_mi_var[i]=summary$t;
lambdatrue_mi_df[i]=summary$df;
lambdatrue_mi_f[i]=summary$f;
}

lambdatrue_mi_true=mean(lambdatrue_mi_mean);
lambdatrue_mi_bias=lambdatrue_mi_true-true;
lambdatrue_mi_bias;
lambdatrue_mi_bias/true;

lambdatrue_mi_mse=mean((lambdatrue_mi_mean-true)^2);
lambdatrue_mi_mse;

# coverage;
lambdatrue_mean_low95_vec=lambdatrue_mi_mean-qt(.975, lambdatrue_mi_df)*sqrt(lambdatrue_mi_var);
lambdatrue_mean_up95_vec=lambdatrue_mi_mean+qt(.975, lambdatrue_mi_df)*sqrt(lambdatrue_mi_var);

mean_lambdatrue_length=mean(lambdatrue_mean_up95_vec-lambdatrue_mean_low95_vec);
mean_lambdatrue_length;


lambdatrue_mi_coverage=(lambdatrue_mean_low95_vec < true)*(lambdatrue_mean_up95_vec > true);
mean(lambdatrue_mi_coverage);

# variance estimates;
mean(sqrt(lambdatrue_mi_var));
sqrt(var(lambdatrue_mi_mean));

# imputation using bootstrap;

boot_mi_mean=rep(NA, cycle_no);
boot_mi_var=rep(NA, cycle_no);
boot_mi_df=rep(NA, cycle_no);
boot_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(boot_mean_mat[i,], boot_var_mat[i,], n=rowobs);
boot_mi_mean[i]=summary$qbar;
boot_mi_var[i]=summary$t;
boot_mi_df[i]=summary$df;
boot_mi_f[i]=summary$f;
}

boot_mi_true=mean(boot_mi_mean);
boot_mi_bias=boot_mi_true-true;
boot_mi_bias;
boot_mi_bias/true;

boot_mi_mse=mean((boot_mi_mean-true)^2);
boot_mi_mse;

# coverage;
boot_mean_low95_vec=boot_mi_mean-qt(.975, boot_mi_df)*sqrt(boot_mi_var);
boot_mean_up95_vec=boot_mi_mean+qt(.975, boot_mi_df)*sqrt(boot_mi_var);

mean_boot_length=mean(boot_mean_up95_vec-boot_mean_low95_vec);
mean_boot_length;


boot_mi_coverage=(boot_mean_low95_vec < true)*(boot_mean_up95_vec > true);
mean(boot_mi_coverage);

# variance estimates;
mean(sqrt(boot_mi_var));
sqrt(var(boot_mi_mean));

# imputation using bayesian;

bayes_mi_mean=rep(NA, cycle_no);
bayes_mi_var=rep(NA, cycle_no);
bayes_mi_df=rep(NA, cycle_no);
bayes_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(bayes_mean_mat[i,], bayes_var_mat[i,], n=rowobs);
bayes_mi_mean[i]=summary$qbar;
bayes_mi_var[i]=summary$t;
bayes_mi_df[i]=summary$df;
bayes_mi_f[i]=summary$f;
}

bayes_mi_true=mean(bayes_mi_mean);
bayes_mi_bias=bayes_mi_true-true;
bayes_mi_bias;
bayes_mi_bias/true;

bayes_mi_mse=mean((bayes_mi_mean-true)^2);
bayes_mi_mse;

# coverage;
bayes_mean_low95_vec=bayes_mi_mean-qt(.975, bayes_mi_df)*sqrt(bayes_mi_var);
bayes_mean_up95_vec=bayes_mi_mean+qt(.975, bayes_mi_df)*sqrt(bayes_mi_var);

mean_bayes_length=mean(bayes_mean_up95_vec-bayes_mean_low95_vec);
mean_bayes_length;


bayes_mi_coverage=(bayes_mean_low95_vec < true)*(bayes_mean_up95_vec > true);
mean(bayes_mi_coverage);

# variance estimates;
mean(sqrt(bayes_mi_var));
sqrt(var(bayes_mi_mean));

###########################################################
# top 5\% tile
# proportion estimates;
true=mean(y2_prop_up_vec);
true;
true_mse=mean((y2_prop_up_vec-true)^2);
true_mse;
y2_prop_var_vec=y2_prop_up_vec*(1-y2_prop_up_vec)/(rowobs-1);

# lower and upper 95% CI;
y2_prop_low95_vec=y2_prop_up_vec-1.96*sqrt(y2_prop_var_vec);
y2_prop_up95_vec=y2_prop_up_vec+1.96*sqrt(y2_prop_var_vec);

mean_length=mean(y2_prop_up95_vec-y2_prop_low95_vec);
mean_length;

coverage=(y2_prop_low95_vec < true)*(y2_prop_up95_vec > true);

mean(coverage);

# variance estimates;
mean(sqrt(y2_prop_var_vec));
sqrt(var(y2_prop_up_vec));

# complete case analysis;
cc_true=mean(y2_cc_prop_up_vec);
cc_bias=cc_true-true;
cc_bias/true;


cc_mse=mean((y2_cc_prop_up_vec-true)^2);
cc_mse;
cc_var_vec=y2_cc_prop_up_vec*(1-y2_cc_prop_up_vec)/(obsno_vec-1);


# lower and upper 95% CI;
cc_mean_low95_vec=y2_cc_prop_up_vec-1.96*sqrt(cc_var_vec);
cc_mean_up95_vec=y2_cc_prop_up_vec+1.96*sqrt(cc_var_vec);

mean_cc_length=mean(cc_mean_up95_vec-cc_mean_low95_vec);
mean_cc_length;

cc_coverage=(cc_mean_low95_vec < true)*(cc_mean_up95_vec > true);

mean(cc_coverage);

# variance estimates;
mean(sqrt(cc_var_vec));
sqrt(var(y2_cc_prop_up_vec));



# marginal transformation imputation;

# observed data transformation imputation;

marginal_mi_mean=rep(NA, cycle_no);
marginal_mi_var=rep(NA, cycle_no);
marginal_mi_df=rep(NA, cycle_no);
marginal_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(marginal_prop_up_mat[i,], marginal_prop_up_mat[i,]*(1-marginal_prop_up_mat[i,])/(rowobs-1), n=rowobs);
marginal_mi_mean[i]=summary$qbar;
marginal_mi_var[i]=summary$t;
marginal_mi_df[i]=summary$df;
marginal_mi_f[i]=summary$f;
}

marginal_mi_true=mean(marginal_mi_mean);
marginal_mi_bias=marginal_mi_true-true;

marginal_mi_bias/true;

marginal_mi_mse=mean((marginal_mi_mean-true)^2);
marginal_mi_mse;

# coverage;
marginal_mean_low95_vec=marginal_mi_mean-qt(.975, marginal_mi_df)*sqrt(marginal_mi_var);
marginal_mean_up95_vec=marginal_mi_mean+qt(.975, marginal_mi_df)*sqrt(marginal_mi_var);

mean_marginal_length=mean(marginal_mean_up95_vec-marginal_mean_low95_vec);
mean_marginal_length

marginal_mi_coverage=(marginal_mean_low95_vec < true)*(marginal_mean_up95_vec > true);
mean(marginal_mi_coverage);

# variance estimates;
mean(sqrt(marginal_mi_var));
sqrt(var(marginal_mi_mean));



# observed data transformation imputation;

obs_mi_mean=rep(NA, cycle_no);
obs_mi_var=rep(NA, cycle_no);
obs_mi_df=rep(NA, cycle_no);
obs_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(obs_prop_up_mat[i,], obs_prop_up_mat[i,]*(1-obs_prop_up_mat[i,])/(rowobs-1), n=rowobs);
obs_mi_mean[i]=summary$qbar;
obs_mi_var[i]=summary$t;
obs_mi_df[i]=summary$df;
obs_mi_f[i]=summary$f;
}

obs_mi_true=mean(obs_mi_mean);
obs_mi_bias=obs_mi_true-true;

obs_mi_bias/true;

obs_mi_mse=mean((obs_mi_mean-true)^2);
obs_mi_mse;

# coverage;
obs_mean_low95_vec=obs_mi_mean-qt(.975, obs_mi_df)*sqrt(obs_mi_var);
obs_mean_up95_vec=obs_mi_mean+qt(.975, obs_mi_df)*sqrt(obs_mi_var);

mean_obs_length=mean(obs_mean_up95_vec-obs_mean_low95_vec);
mean_obs_length

obs_mi_coverage=(obs_mean_low95_vec < true)*(obs_mean_up95_vec > true);
mean(obs_mi_coverage);

# variance estimates;
mean(sqrt(obs_mi_var));
sqrt(var(obs_mi_mean));

# linear normal imputation;

guess_mi_mean=rep(NA, cycle_no);
guess_mi_var=rep(NA, cycle_no);
guess_mi_df=rep(NA, cycle_no);
guess_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(guess_prop_up_mat[i,], guess_prop_up_mat[i,]*(1-guess_prop_up_mat[i,])/(rowobs-1), n=rowobs);
guess_mi_mean[i]=summary$qbar;
guess_mi_var[i]=summary$t;
guess_mi_df[i]=summary$df;
guess_mi_f[i]=summary$f;
}

guess_mi_true=mean(guess_mi_mean);
guess_mi_bias=guess_mi_true-true;

guess_mi_bias/true;


guess_mi_mse=mean((guess_mi_mean-true)^2);
guess_mi_mse;

# coverage;
guess_mean_low95_vec=guess_mi_mean-qt(.975, guess_mi_df)*sqrt(guess_mi_var);
guess_mean_up95_vec=guess_mi_mean+qt(.975, guess_mi_df)*sqrt(guess_mi_var);

mean_guess_length=mean(guess_mean_up95_vec-guess_mean_low95_vec);
mean_guess_length

guess_mi_coverage=(guess_mean_low95_vec < true)*(guess_mean_up95_vec > true);
mean(guess_mi_coverage);

# variance estimates;
mean(sqrt(guess_mi_var));
sqrt(var(guess_mi_mean));



# imputation using draws;

draw_mi_mean=rep(NA, cycle_no);
draw_mi_var=rep(NA, cycle_no);
draw_mi_df=rep(NA, cycle_no);
draw_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(draw_prop_up_mat[i,], draw_prop_up_mat[i,]*(1-draw_prop_up_mat[i,])/(rowobs-1), n=rowobs);
draw_mi_mean[i]=summary$qbar;
draw_mi_var[i]=summary$t;
draw_mi_df[i]=summary$df;
draw_mi_f[i]=summary$f;
}

draw_mi_true=mean(draw_mi_mean);
draw_mi_bias=draw_mi_true-true;

draw_mi_bias/true;


draw_mi_mse=mean((draw_mi_mean-true)^2);
draw_mi_mse;

# coverage;
draw_mean_low95_vec=draw_mi_mean-qt(.975, draw_mi_df)*sqrt(draw_mi_var);
draw_mean_up95_vec=draw_mi_mean+qt(.975, draw_mi_df)*sqrt(draw_mi_var);

mean_draw_length=mean(draw_mean_up95_vec-draw_mean_low95_vec);
mean_draw_length

draw_mi_coverage=(draw_mean_low95_vec < true)*(draw_mean_up95_vec > true);
mean(draw_mi_coverage);

# variance estimates;
mean(sqrt(draw_mi_var));
sqrt(var(draw_mi_mean));

# imputation using true lambda;

lambdatrue_mi_mean=rep(NA, cycle_no);
lambdatrue_mi_var=rep(NA, cycle_no);
lambdatrue_mi_df=rep(NA, cycle_no);
lambdatrue_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(lambdatrue_prop_up_mat[i,], lambdatrue_prop_up_mat[i,]*(1-lambdatrue_prop_up_mat[i,])/(rowobs-1), n=rowobs);
lambdatrue_mi_mean[i]=summary$qbar;
lambdatrue_mi_var[i]=summary$t;
lambdatrue_mi_df[i]=summary$df;
lambdatrue_mi_f[i]=summary$f;
}

lambdatrue_mi_true=mean(lambdatrue_mi_mean);
lambdatrue_mi_bias=lambdatrue_mi_true-true;
lambdatrue_mi_bias/true;

lambdatrue_mi_mse=mean((lambdatrue_mi_mean-true)^2);
lambdatrue_mi_mse;

# coverage;
lambdatrue_mean_low95_vec=lambdatrue_mi_mean-qt(.975, lambdatrue_mi_df)*sqrt(lambdatrue_mi_var);
lambdatrue_mean_up95_vec=lambdatrue_mi_mean+qt(.975, lambdatrue_mi_df)*sqrt(lambdatrue_mi_var);

mean_lambdatrue_length=mean(lambdatrue_mean_up95_vec-lambdatrue_mean_low95_vec);
mean_lambdatrue_length

lambdatrue_mi_coverage=(lambdatrue_mean_low95_vec < true)*(lambdatrue_mean_up95_vec > true);
mean(lambdatrue_mi_coverage);

# variance estimates;
mean(sqrt(lambdatrue_mi_var));
sqrt(var(lambdatrue_mi_mean));

# imputation using bootstrap;

boot_mi_mean=rep(NA, cycle_no);
boot_mi_var=rep(NA, cycle_no);
boot_mi_df=rep(NA, cycle_no);
boot_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(boot_prop_up_mat[i,], boot_prop_up_mat[i,]*(1-boot_prop_up_mat[i,])/(rowobs-1), n=rowobs);
boot_mi_mean[i]=summary$qbar;
boot_mi_var[i]=summary$t;
boot_mi_df[i]=summary$df;
boot_mi_f[i]=summary$f;
}

boot_mi_true=mean(boot_mi_mean);
boot_mi_bias=boot_mi_true-true;
boot_mi_bias/true;

boot_mi_mse=mean((boot_mi_mean-true)^2);
boot_mi_mse;

# coverage;
boot_mean_low95_vec=boot_mi_mean-qt(.975, boot_mi_df)*sqrt(boot_mi_var);
boot_mean_up95_vec=boot_mi_mean+qt(.975, boot_mi_df)*sqrt(boot_mi_var);

mean_boot_length=mean(boot_mean_up95_vec-boot_mean_low95_vec);
mean_boot_length

boot_mi_coverage=(boot_mean_low95_vec < true)*(boot_mean_up95_vec > true);
mean(boot_mi_coverage);

# variance estimates;
mean(sqrt(boot_mi_var));
sqrt(var(boot_mi_mean));

# imputation using a fully bayesian;

bayes_mi_mean=rep(NA, cycle_no);
bayes_mi_var=rep(NA, cycle_no);
bayes_mi_df=rep(NA, cycle_no);
bayes_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(bayes_prop_up_mat[i,], bayes_prop_up_mat[i,]*(1-bayes_prop_up_mat[i,])/(rowobs-1), n=rowobs);
bayes_mi_mean[i]=summary$qbar;
bayes_mi_var[i]=summary$t;
bayes_mi_df[i]=summary$df;
bayes_mi_f[i]=summary$f;
}

bayes_mi_true=mean(bayes_mi_mean);
bayes_mi_bias=bayes_mi_true-true;
bayes_mi_bias/true;

bayes_mi_mse=mean((bayes_mi_mean-true)^2);
bayes_mi_mse;

# coverage;
bayes_mean_low95_vec=bayes_mi_mean-qt(.975, bayes_mi_df)*sqrt(bayes_mi_var);
bayes_mean_up95_vec=bayes_mi_mean+qt(.975, bayes_mi_df)*sqrt(bayes_mi_var);

mean_bayes_length=mean(bayes_mean_up95_vec-bayes_mean_low95_vec);
mean_bayes_length

bayes_mi_coverage=(bayes_mean_low95_vec < true)*(bayes_mean_up95_vec > true);
mean(bayes_mi_coverage);

# variance estimates;
mean(sqrt(bayes_mi_var));
sqrt(var(bayes_mi_mean));

#####
## bottom 5%-tile

# proportion estimates;
true=mean(y2_prop_bottom_vec);
true;
true_mse=mean((y2_prop_bottom_vec-true)^2);
true_mse;
y2_prop_var_vec=y2_prop_bottom_vec*(1-y2_prop_bottom_vec)/(rowobs-1);

# lower and upper 95% CI;
y2_prop_low95_vec=y2_prop_bottom_vec-1.96*sqrt(y2_prop_var_vec);
y2_prop_up95_vec=y2_prop_bottom_vec+1.96*sqrt(y2_prop_var_vec);

mean_length=mean(y2_prop_up95_vec-y2_prop_low95_vec);
mean_length

coverage=(y2_prop_low95_vec < true)*(y2_prop_up95_vec > true);

mean(coverage);

# variance estimates;
mean(sqrt(y2_prop_var_vec));
sqrt(var(y2_prop_bottom_vec));

# complete case analysis;
cc_true=mean(y2_cc_prop_bottom_vec);
cc_bias=cc_true-true;
cc_bias/true;


cc_mse=mean((y2_cc_prop_bottom_vec-true)^2);
cc_mse;
cc_var_vec=y2_cc_prop_bottom_vec*(1-y2_cc_prop_bottom_vec)/(obsno_vec-1);


# lower and upper 95% CI;
cc_mean_low95_vec=y2_cc_prop_bottom_vec-1.96*sqrt(cc_var_vec);
cc_mean_up95_vec=y2_cc_prop_bottom_vec+1.96*sqrt(cc_var_vec);

mean_cc_length=mean(cc_mean_up95_vec-cc_mean_low95_vec);
mean_cc_length

cc_coverage=(cc_mean_low95_vec < true)*(cc_mean_up95_vec > true);

mean(cc_coverage);

# variance estimates;
mean(sqrt(cc_var_vec));
sqrt(var(y2_cc_prop_bottom_vec));



# marginal transformation imputation;

# observed data transformation imputation;

marginal_mi_mean=rep(NA, cycle_no);
marginal_mi_var=rep(NA, cycle_no);
marginal_mi_df=rep(NA, cycle_no);
marginal_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(marginal_prop_bottom_mat[i,], marginal_prop_bottom_mat[i,]*(1-marginal_prop_bottom_mat[i,])/(rowobs-1), n=rowobs);
marginal_mi_mean[i]=summary$qbar;
marginal_mi_var[i]=summary$t;
marginal_mi_df[i]=summary$df;
marginal_mi_f[i]=summary$f;
}

marginal_mi_true=mean(marginal_mi_mean);
marginal_mi_bias=marginal_mi_true-true;

marginal_mi_bias/true;

marginal_mi_mse=mean((marginal_mi_mean-true)^2);
marginal_mi_mse;

# coverage;
marginal_mean_low95_vec=marginal_mi_mean-qt(.975, marginal_mi_df)*sqrt(marginal_mi_var);
marginal_mean_up95_vec=marginal_mi_mean+qt(.975, marginal_mi_df)*sqrt(marginal_mi_var);

mean_marginal_length=mean(marginal_mean_up95_vec-marginal_mean_low95_vec);
mean_marginal_length

marginal_mi_coverage=(marginal_mean_low95_vec < true)*(marginal_mean_up95_vec > true);
mean(marginal_mi_coverage);

# variance estimates;
mean(sqrt(marginal_mi_var));
sqrt(var(marginal_mi_mean));



# observed data transformation imputation;

obs_mi_mean=rep(NA, cycle_no);
obs_mi_var=rep(NA, cycle_no);
obs_mi_df=rep(NA, cycle_no);
obs_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(obs_prop_bottom_mat[i,], obs_prop_bottom_mat[i,]*(1-obs_prop_bottom_mat[i,])/(rowobs-1), n=rowobs);
obs_mi_mean[i]=summary$qbar;
obs_mi_var[i]=summary$t;
obs_mi_df[i]=summary$df;
obs_mi_f[i]=summary$f;
}

obs_mi_true=mean(obs_mi_mean);
obs_mi_bias=obs_mi_true-true;

obs_mi_bias/true;

obs_mi_mse=mean((obs_mi_mean-true)^2);
obs_mi_mse;

# coverage;
obs_mean_low95_vec=obs_mi_mean-qt(.975, obs_mi_df)*sqrt(obs_mi_var);
obs_mean_up95_vec=obs_mi_mean+qt(.975, obs_mi_df)*sqrt(obs_mi_var);

mean_obs_length=mean(obs_mean_up95_vec-obs_mean_low95_vec);
mean_obs_length

obs_mi_coverage=(obs_mean_low95_vec < true)*(obs_mean_up95_vec > true);
mean(obs_mi_coverage);

# variance estimates;
mean(sqrt(obs_mi_var));
sqrt(var(obs_mi_mean));

# linear normal imputation;

guess_mi_mean=rep(NA, cycle_no);
guess_mi_var=rep(NA, cycle_no);
guess_mi_df=rep(NA, cycle_no);
guess_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(guess_prop_bottom_mat[i,], guess_prop_bottom_mat[i,]*(1-guess_prop_bottom_mat[i,])/(rowobs-1), n=rowobs);
guess_mi_mean[i]=summary$qbar;
guess_mi_var[i]=summary$t;
guess_mi_df[i]=summary$df;
guess_mi_f[i]=summary$f;
}

guess_mi_true=mean(guess_mi_mean);
guess_mi_bias=guess_mi_true-true;

guess_mi_bias/true;


guess_mi_mse=mean((guess_mi_mean-true)^2);
guess_mi_mse;

# coverage;
guess_mean_low95_vec=guess_mi_mean-qt(.975, guess_mi_df)*sqrt(guess_mi_var);
guess_mean_up95_vec=guess_mi_mean+qt(.975, guess_mi_df)*sqrt(guess_mi_var);

mean_guess_length=mean(guess_mean_up95_vec-guess_mean_low95_vec);
mean_guess_length

guess_mi_coverage=(guess_mean_low95_vec < true)*(guess_mean_up95_vec > true);
mean(guess_mi_coverage);

# variance estimates;
mean(sqrt(guess_mi_var));
sqrt(var(guess_mi_mean));



# imputation using draws;

draw_mi_mean=rep(NA, cycle_no);
draw_mi_var=rep(NA, cycle_no);
draw_mi_df=rep(NA, cycle_no);
draw_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(draw_prop_bottom_mat[i,], draw_prop_bottom_mat[i,]*(1-draw_prop_bottom_mat[i,])/(rowobs-1), n=rowobs);
draw_mi_mean[i]=summary$qbar;
draw_mi_var[i]=summary$t;
draw_mi_df[i]=summary$df;
draw_mi_f[i]=summary$f;
}

draw_mi_true=mean(draw_mi_mean);
draw_mi_bias=draw_mi_true-true;

draw_mi_bias/true;


draw_mi_mse=mean((draw_mi_mean-true)^2);
draw_mi_mse;

# coverage;
draw_mean_low95_vec=draw_mi_mean-qt(.975, draw_mi_df)*sqrt(draw_mi_var);
draw_mean_up95_vec=draw_mi_mean+qt(.975, draw_mi_df)*sqrt(draw_mi_var);

mean_draw_length=mean(draw_mean_up95_vec-draw_mean_low95_vec);
mean_draw_length

draw_mi_coverage=(draw_mean_low95_vec < true)*(draw_mean_up95_vec > true);
mean(draw_mi_coverage);

# variance estimates;
mean(sqrt(draw_mi_var));
sqrt(var(draw_mi_mean));

# imputation using true lambda;

lambdatrue_mi_mean=rep(NA, cycle_no);
lambdatrue_mi_var=rep(NA, cycle_no);
lambdatrue_mi_df=rep(NA, cycle_no);
lambdatrue_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(lambdatrue_prop_bottom_mat[i,], lambdatrue_prop_bottom_mat[i,]*(1-lambdatrue_prop_bottom_mat[i,])/(rowobs-1), n=rowobs);
lambdatrue_mi_mean[i]=summary$qbar;
lambdatrue_mi_var[i]=summary$t;
lambdatrue_mi_df[i]=summary$df;
lambdatrue_mi_f[i]=summary$f;
}

lambdatrue_mi_true=mean(lambdatrue_mi_mean);
lambdatrue_mi_bias=lambdatrue_mi_true-true;
lambdatrue_mi_bias/true;

lambdatrue_mi_mse=mean((lambdatrue_mi_mean-true)^2);
lambdatrue_mi_mse;

# coverage;
lambdatrue_mean_low95_vec=lambdatrue_mi_mean-qt(.975, lambdatrue_mi_df)*sqrt(lambdatrue_mi_var);
lambdatrue_mean_up95_vec=lambdatrue_mi_mean+qt(.975, lambdatrue_mi_df)*sqrt(lambdatrue_mi_var);

mean_lambdatrue_length=mean(lambdatrue_mean_up95_vec-lambdatrue_mean_low95_vec);
mean_lambdatrue_length

lambdatrue_mi_coverage=(lambdatrue_mean_low95_vec < true)*(lambdatrue_mean_up95_vec > true);
mean(lambdatrue_mi_coverage);

# variance estimates;
mean(sqrt(lambdatrue_mi_var));
sqrt(var(lambdatrue_mi_mean));

# imputation using bootstrap;

boot_mi_mean=rep(NA, cycle_no);
boot_mi_var=rep(NA, cycle_no);
boot_mi_df=rep(NA, cycle_no);
boot_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(boot_prop_bottom_mat[i,], boot_prop_bottom_mat[i,]*(1-boot_prop_bottom_mat[i,])/(rowobs-1), n=rowobs);
boot_mi_mean[i]=summary$qbar;
boot_mi_var[i]=summary$t;
boot_mi_df[i]=summary$df;
boot_mi_f[i]=summary$f;
}

boot_mi_true=mean(boot_mi_mean);
boot_mi_bias=boot_mi_true-true;
boot_mi_bias/true;

boot_mi_mse=mean((boot_mi_mean-true)^2);
boot_mi_mse;

# coverage;
boot_mean_low95_vec=boot_mi_mean-qt(.975, boot_mi_df)*sqrt(boot_mi_var);
boot_mean_up95_vec=boot_mi_mean+qt(.975, boot_mi_df)*sqrt(boot_mi_var);

mean_boot_length=mean(boot_mean_up95_vec-boot_mean_low95_vec);
mean_boot_length

boot_mi_coverage=(boot_mean_low95_vec < true)*(boot_mean_up95_vec > true);
mean(boot_mi_coverage);

# variance estimates;
mean(sqrt(boot_mi_var));
sqrt(var(boot_mi_mean));

# imputation using a fully bayesian;

bayes_mi_mean=rep(NA, cycle_no);
bayes_mi_var=rep(NA, cycle_no);
bayes_mi_df=rep(NA, cycle_no);
bayes_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(bayes_prop_bottom_mat[i,], bayes_prop_bottom_mat[i,]*(1-bayes_prop_bottom_mat[i,])/(rowobs-1), n=rowobs);
bayes_mi_mean[i]=summary$qbar;
bayes_mi_var[i]=summary$t;
bayes_mi_df[i]=summary$df;
bayes_mi_f[i]=summary$f;
}

bayes_mi_true=mean(bayes_mi_mean);
bayes_mi_bias=bayes_mi_true-true;
bayes_mi_bias/true;

bayes_mi_mse=mean((bayes_mi_mean-true)^2);
bayes_mi_mse;

# coverage;
bayes_mean_low95_vec=bayes_mi_mean-qt(.975, bayes_mi_df)*sqrt(bayes_mi_var);
bayes_mean_up95_vec=bayes_mi_mean+qt(.975, bayes_mi_df)*sqrt(bayes_mi_var);

mean_bayes_length=mean(bayes_mean_up95_vec-bayes_mean_low95_vec);
mean_bayes_length

bayes_mi_coverage=(bayes_mean_low95_vec < true)*(bayes_mean_up95_vec > true);
mean(bayes_mi_coverage);

# variance estimates;
mean(sqrt(bayes_mi_var));
sqrt(var(bayes_mi_mean));





################## do not look at them for Example 6.1 #######################
# regression estimandtrue;
# using fisher information based matrix
# assume lambda is unknown
# slope;
# complete data;
true=mean(y2_slope_mean_vec);
true_mse=mean((y2_slope_mean_vec-true)^2);

true_mse

# lower and upper 95% CI;
y2_slope_mean_low95_vec=y2_slope_mean_vec-1.96*sqrt(y2_fisher_slope_var_vec);
y2_slope_mean_up95_vec=y2_slope_mean_vec+1.96*sqrt(y2_fisher_slope_var_vec);

# y2_slope_mean_low95_vec=y2_slope_mean_vec-1.96*sqrt(y2_slope_var_vec);
# y2_slope_mean_up95_vec=y2_slope_mean_vec+1.96*sqrt(y2_slope_var_vec);


mean_length=mean(y2_slope_mean_up95_vec-y2_slope_mean_low95_vec);

mean_length

coverage=(y2_slope_mean_low95_vec < true)*(y2_slope_mean_up95_vec > true);

mean(coverage);

# variance estimates;
mean(y2_fisher_slope_var_vec);
mean(y2_slope_var_vec);
var(y2_slope_mean_vec);

par(mfrow=c(1,1));
hist(y2_slope_mean_vec);


# cc analysis;
cc_true=mean(cc_slope_mean_vec);
cc_bias=cc_true-true;
cc_mse=mean((cc_slope_mean_vec-true)^2);


# lower and upper 95% CI;
cc_mean_low95_vec=cc_slope_mean_vec-1.96*sqrt(y2_cc_fisher_slope_var_vec);
cc_mean_up95_vec=cc_slope_mean_vec+1.96*sqrt(y2_cc_fisher_slope_var_vec);

mean_cc_length=mean(cc_mean_up95_vec-cc_mean_low95_vec);

cc_coverage=(cc_mean_low95_vec < true)*(cc_mean_up95_vec > true);

mean(cc_coverage);

# variance estimates;
mean(y2_cc_fisher_slope_var_vec);
var(cc_slope_mean_vec);



# multiple imputation estimates;
# marginal transformation imputation;

# marginal_mi_mean=rep(NA, cycle_no);
# marginal_mi_var=rep(NA, cycle_no);
# marginal_mi_df=rep(NA, cycle_no);
# marginal_mi_f=rep(NA, cycle_no);

# for (i in 1:cycle_no)
# {
# summary=pool.scalar(marginal_mean_mat[i,], marginal_var_mat[i,]);
# marginal_mi_mean[i]=summary$qbar;
# marginal_mi_var[i]=summary$t;
# marginal_mi_df[i]=summary$df;
# marginal_mi_f[i]=summary$f;
# }

# marginal_mi_true=mean(marginal_mi_mean);
# marginal_mi_bias=marginal_mi_true-true;

# marginal_mi_mse=mean((marginal_mi_mean-true)^2);

# coverage;
# marginal_mean_low95_vec=marginal_mi_mean-qt(.975, marginal_mi_df)*sqrt(marginal_mi_var);
# marginal_mean_up95_vec=marginal_mi_mean+qt(.975, marginal_mi_df)*sqrt(marginal_mi_var);

# mean_marginal_length=mean(marginal_mean_up95_vec-marginal_mean_low95_vec);


# marginal_mi_coverage=(marginal_mean_low95_vec < true)*(marginal_mean_up95_vec > true);
# mean(marginal_mi_coverage);

# variance estimates;
# mean(marginal_mi_var);
# var(marginal_mi_mean);

# observed data transformation imputation;

obs_mi_mean=rep(NA, cycle_no);
obs_mi_var=rep(NA, cycle_no);
obs_mi_df=rep(NA, cycle_no);
obs_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(obs_com_slope_mean_mat[i,], obs_com_fisher_slope_var_mat[i,]);
obs_mi_mean[i]=summary$qbar;
obs_mi_var[i]=summary$t;
obs_mi_df[i]=summary$df;
obs_mi_f[i]=summary$f;
}

obs_mi_true=mean(obs_mi_mean);
obs_mi_bias=obs_mi_true-true;

obs_mi_mse=mean((obs_mi_mean-true)^2);

# coverage;
obs_mean_low95_vec=obs_mi_mean-qt(.975, obs_mi_df)*sqrt(obs_mi_var);
obs_mean_up95_vec=obs_mi_mean+qt(.975, obs_mi_df)*sqrt(obs_mi_var);

mean_obs_length=mean(obs_mean_up95_vec-obs_mean_low95_vec);


obs_mi_coverage=(obs_mean_low95_vec < true)*(obs_mean_up95_vec > true);
mean(obs_mi_coverage);

# variance estimates;
mean(obs_mi_var);
var(obs_mi_mean);

# Imputation using guess of lambda;

guess_mi_mean=rep(NA, cycle_no);
guess_mi_var=rep(NA, cycle_no);
guess_mi_df=rep(NA, cycle_no);
guess_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(guess_com_slope_mean_mat[i,], guess_com_fisher_slope_var_mat[i,]);
guess_mi_mean[i]=summary$qbar;
guess_mi_var[i]=summary$t;
guess_mi_df[i]=summary$df;
guess_mi_f[i]=summary$f;
}

guess_mi_true=mean(guess_mi_mean);
guess_mi_bias=guess_mi_true-true;

guess_mi_mse=mean((guess_mi_mean-true)^2);

guess_mi_mse;

# coverage;
guess_mean_low95_vec=guess_mi_mean-qt(.975, guess_mi_df)*sqrt(guess_mi_var);
guess_mean_up95_vec=guess_mi_mean+qt(.975, guess_mi_df)*sqrt(guess_mi_var);

mean_guess_length=mean(guess_mean_up95_vec-guess_mean_low95_vec);


guess_mi_coverage=(guess_mean_low95_vec < true)*(guess_mean_up95_vec > true);
mean(guess_mi_coverage);

# variance estimates;
mean(guess_mi_var);
var(guess_mi_mean);



# imputation using draws;

draw_mi_mean=rep(NA, cycle_no);
draw_mi_var=rep(NA, cycle_no);
draw_mi_df=rep(NA, cycle_no);
draw_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(draw_com_slope_mean_mat[i,], draw_com_fisher_slope_var_mat[i,]);
# summary=pool.scalar(draw_com_slope_mean_mat[i,], draw_com_slope_var_mat[i,]);


draw_mi_mean[i]=summary$qbar;
draw_mi_var[i]=summary$t;
draw_mi_df[i]=summary$df;
draw_mi_f[i]=summary$f;
}

draw_mi_true=mean(draw_mi_mean);
draw_mi_bias=draw_mi_true-true;

draw_mi_mse=mean((draw_mi_mean-true)^2);

# coverage;
draw_mean_low95_vec=draw_mi_mean-qt(.975, draw_mi_df)*sqrt(draw_mi_var);
draw_mean_up95_vec=draw_mi_mean+qt(.975, draw_mi_df)*sqrt(draw_mi_var);

mean_draw_length=mean(draw_mean_up95_vec-draw_mean_low95_vec);


draw_mi_coverage=(draw_mean_low95_vec < true)*(draw_mean_up95_vec > true);
mean(draw_mi_coverage);

# variance estimates;
mean(draw_mi_var);
var(draw_mi_mean);

# imputation using true lambda;

lambdatrue_mi_mean=rep(NA, cycle_no);
lambdatrue_mi_var=rep(NA, cycle_no);
lambdatrue_mi_df=rep(NA, cycle_no);
lambdatrue_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(lambdatrue_com_slope_mean_mat[i,], lambdatrue_com_fisher_slope_var_mat[i,]);
lambdatrue_mi_mean[i]=summary$qbar;
lambdatrue_mi_var[i]=summary$t;
lambdatrue_mi_df[i]=summary$df;
lambdatrue_mi_f[i]=summary$f;
}

lambdatrue_mi_true=mean(lambdatrue_mi_mean);
lambdatrue_mi_bias=lambdatrue_mi_true-true;

lambdatrue_mi_mse=mean((lambdatrue_mi_mean-true)^2);

# coverage;
lambdatrue_mean_low95_vec=lambdatrue_mi_mean-qt(.975, lambdatrue_mi_df)*sqrt(lambdatrue_mi_var);
lambdatrue_mean_up95_vec=lambdatrue_mi_mean+qt(.975, lambdatrue_mi_df)*sqrt(lambdatrue_mi_var);

mean_lambdatrue_length=mean(lambdatrue_mean_up95_vec-lambdatrue_mean_low95_vec);


lambdatrue_mi_coverage=(lambdatrue_mean_low95_vec < true)*(lambdatrue_mean_up95_vec > true);
mean(lambdatrue_mi_coverage);

# variance estimates;
mean(lambdatrue_mi_var);
var(lambdatrue_mi_mean);

# imputation using bootstrap;

boot_mi_mean=rep(NA, cycle_no);
boot_mi_var=rep(NA, cycle_no);
boot_mi_df=rep(NA, cycle_no);
boot_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(boot_com_slope_mean_mat[i,], boot_com_fisher_slope_var_mat[i,]);
boot_mi_mean[i]=summary$qbar;
boot_mi_var[i]=summary$t;
boot_mi_df[i]=summary$df;
boot_mi_f[i]=summary$f;
}

boot_mi_true=mean(boot_mi_mean);
boot_mi_bias=boot_mi_true-true;

boot_mi_mse=mean((boot_mi_mean-true)^2);

# coverage;
boot_mean_low95_vec=boot_mi_mean-qt(.975, boot_mi_df)*sqrt(boot_mi_var);
boot_mean_up95_vec=boot_mi_mean+qt(.975, boot_mi_df)*sqrt(boot_mi_var);

mean_boot_length=mean(boot_mean_up95_vec-boot_mean_low95_vec);


boot_mi_coverage=(boot_mean_low95_vec < true)*(boot_mean_up95_vec > true);
mean(boot_mi_coverage);

# variance estimates;
mean(boot_mi_var);
var(boot_mi_mean);

# imputation using a fully bayesian;

bayes_mi_mean=rep(NA, cycle_no);
bayes_mi_var=rep(NA, cycle_no);
bayes_mi_df=rep(NA, cycle_no);
bayes_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(bayes_com_slope_mean_mat[i,], bayes_com_fisher_slope_var_mat[i,]);
bayes_mi_mean[i]=summary$qbar;
bayes_mi_var[i]=summary$t;
bayes_mi_df[i]=summary$df;
bayes_mi_f[i]=summary$f;
}

bayes_mi_true=mean(bayes_mi_mean);
bayes_mi_bias=bayes_mi_true-true;

bayes_mi_mse=mean((bayes_mi_mean-true)^2);

# coverage;
bayes_mean_low95_vec=bayes_mi_mean-qt(.975, bayes_mi_df)*sqrt(bayes_mi_var);
bayes_mean_up95_vec=bayes_mi_mean+qt(.975, bayes_mi_df)*sqrt(bayes_mi_var);

mean_bayes_length=mean(bayes_mean_up95_vec-bayes_mean_low95_vec);


bayes_mi_coverage=(bayes_mean_low95_vec < true)*(bayes_mean_up95_vec > true);
mean(bayes_mi_coverage);

# variance estimates;
mean(bayes_mi_var);
var(bayes_mi_mean);



# regression estimand;
# assume lambda is known and fixed at the true lambda
# slope;
# complete data;
true=mean(y2_correct_slope_mean_vec);
true_mse=mean((y2_correct_slope_mean_vec-true)^2);

# lower and upper 95% CI;
y2_correct_slope_mean_low95_vec=y2_correct_slope_mean_vec-1.96*sqrt(y2_correct_slope_var_vec);
y2_correct_slope_mean_up95_vec=y2_correct_slope_mean_vec+1.96*sqrt(y2_correct_slope_var_vec);

mean_length=mean(y2_correct_slope_mean_up95_vec-y2_correct_slope_mean_low95_vec);

coverage=(y2_correct_slope_mean_low95_vec < true)*(y2_correct_slope_mean_up95_vec > true);

mean(coverage);

# variance estimates;
mean(y2_correct_slope_var_vec);
var(y2_correct_slope_mean_vec);

# cc analysis;
cc_true=mean(cc_correct_slope_mean_vec);
cc_bias=cc_true-true;
cc_mse=mean((cc_correct_slope_mean_vec-true)^2);


# lower and upper 95% CI;
cc_mean_low95_vec=cc_correct_slope_mean_vec-1.96*sqrt(cc_correct_slope_var_vec);
cc_mean_up95_vec=cc_correct_slope_mean_vec+1.96*sqrt(cc_correct_slope_var_vec);

mean_cc_length=mean(cc_mean_up95_vec-cc_mean_low95_vec);

cc_coverage=(cc_mean_low95_vec < true)*(cc_mean_up95_vec > true);

mean(cc_coverage);

# variance estimates;
mean(cc_correct_slope_var_vec);
var(cc_correct_slope_mean_vec);



# multiple imputation estimates;
# marginal transformation imputation;

# marginal_mi_mean=rep(NA, cycle_no);
# marginal_mi_var=rep(NA, cycle_no);
# marginal_mi_df=rep(NA, cycle_no);
# marginal_mi_f=rep(NA, cycle_no);

# for (i in 1:cycle_no)
# {
# summary=pool.scalar(marginal_mean_mat[i,], marginal_var_mat[i,]);
# marginal_mi_mean[i]=summary$qbar;
# marginal_mi_var[i]=summary$t;
# marginal_mi_df[i]=summary$df;
# marginal_mi_f[i]=summary$f;
# }

# marginal_mi_true=mean(marginal_mi_mean);
# marginal_mi_bias=marginal_mi_true-true;

# marginal_mi_mse=mean((marginal_mi_mean-true)^2);

# coverage;
# marginal_mean_low95_vec=marginal_mi_mean-qt(.975, marginal_mi_df)*sqrt(marginal_mi_var);
# marginal_mean_up95_vec=marginal_mi_mean+qt(.975, marginal_mi_df)*sqrt(marginal_mi_var);

# mean_marginal_length=mean(marginal_mean_up95_vec-marginal_mean_low95_vec);


# marginal_mi_coverage=(marginal_mean_low95_vec < true)*(marginal_mean_up95_vec > true);
# mean(marginal_mi_coverage);

# variance estimates;
# mean(marginal_mi_var);
# var(marginal_mi_mean);


# observed data transformation imputation;

obs_mi_mean=rep(NA, cycle_no);
obs_mi_var=rep(NA, cycle_no);
obs_mi_df=rep(NA, cycle_no);
obs_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(obs_com_correct_slope_mean_mat[i,], obs_com_correct_slope_var_mat[i,]);
obs_mi_mean[i]=summary$qbar;
obs_mi_var[i]=summary$t;
obs_mi_df[i]=summary$df;
obs_mi_f[i]=summary$f;
}

obs_mi_true=mean(obs_mi_mean);
obs_mi_bias=obs_mi_true-true;

obs_mi_mse=mean((obs_mi_mean-true)^2);

# coverage;
obs_mean_low95_vec=obs_mi_mean-qt(.975, obs_mi_df)*sqrt(obs_mi_var);
obs_mean_up95_vec=obs_mi_mean+qt(.975, obs_mi_df)*sqrt(obs_mi_var);

mean_obs_length=mean(obs_mean_up95_vec-obs_mean_low95_vec);


obs_mi_coverage=(obs_mean_low95_vec < true)*(obs_mean_up95_vec > true);
mean(obs_mi_coverage);

# variance estimates;
mean(obs_mi_var);
var(obs_mi_mean);

# Imutation using guess of lambda;

guess_mi_mean=rep(NA, cycle_no);
guess_mi_var=rep(NA, cycle_no);
guess_mi_df=rep(NA, cycle_no);
guess_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(guess_com_correct_slope_mean_mat[i,], guess_com_correct_slope_var_mat[i,]);
guess_mi_mean[i]=summary$qbar;
guess_mi_var[i]=summary$t;
guess_mi_df[i]=summary$df;
guess_mi_f[i]=summary$f;
}

guess_mi_true=mean(guess_mi_mean);
guess_mi_bias=guess_mi_true-true;

guess_mi_mse=mean((guess_mi_mean-true)^2);

# coverage;
guess_mean_low95_vec=guess_mi_mean-qt(.975, guess_mi_df)*sqrt(guess_mi_var);
guess_mean_up95_vec=guess_mi_mean+qt(.975, guess_mi_df)*sqrt(guess_mi_var);

mean_guess_length=mean(guess_mean_up95_vec-guess_mean_low95_vec);


guess_mi_coverage=(guess_mean_low95_vec < true)*(guess_mean_up95_vec > true);
mean(guess_mi_coverage);

# variance estimates;
mean(guess_mi_var);
var(guess_mi_mean);


# imputation using draws;

draw_mi_mean=rep(NA, cycle_no);
draw_mi_var=rep(NA, cycle_no);
draw_mi_df=rep(NA, cycle_no);
draw_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(draw_com_correct_slope_mean_mat[i,], draw_com_correct_slope_var_mat[i,]);
draw_mi_mean[i]=summary$qbar;
draw_mi_var[i]=summary$t;
draw_mi_df[i]=summary$df;
draw_mi_f[i]=summary$f;
}

draw_mi_true=mean(draw_mi_mean);
draw_mi_bias=draw_mi_true-true;

draw_mi_mse=mean((draw_mi_mean-true)^2);

# coverage;
draw_mean_low95_vec=draw_mi_mean-qt(.975, draw_mi_df)*sqrt(draw_mi_var);
draw_mean_up95_vec=draw_mi_mean+qt(.975, draw_mi_df)*sqrt(draw_mi_var);

mean_draw_length=mean(draw_mean_up95_vec-draw_mean_low95_vec);


draw_mi_coverage=(draw_mean_low95_vec < true)*(draw_mean_up95_vec > true);
mean(draw_mi_coverage);

# variance estimates;
mean(draw_mi_var);
var(draw_mi_mean);

# imputation using true lambda;

lambdatrue_mi_mean=rep(NA, cycle_no);
lambdatrue_mi_var=rep(NA, cycle_no);
lambdatrue_mi_df=rep(NA, cycle_no);
lambdatrue_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(lambdatrue_com_correct_slope_mean_mat[i,], lambdatrue_com_correct_slope_var_mat[i,]);
lambdatrue_mi_mean[i]=summary$qbar;
lambdatrue_mi_var[i]=summary$t;
lambdatrue_mi_df[i]=summary$df;
lambdatrue_mi_f[i]=summary$f;
}

lambdatrue_mi_true=mean(lambdatrue_mi_mean);
lambdatrue_mi_bias=lambdatrue_mi_true-true;

lambdatrue_mi_mse=mean((lambdatrue_mi_mean-true)^2);

# coverage;
lambdatrue_mean_low95_vec=lambdatrue_mi_mean-qt(.975, lambdatrue_mi_df)*sqrt(lambdatrue_mi_var);
lambdatrue_mean_up95_vec=lambdatrue_mi_mean+qt(.975, lambdatrue_mi_df)*sqrt(lambdatrue_mi_var);

mean_lambdatrue_length=mean(lambdatrue_mean_up95_vec-lambdatrue_mean_low95_vec);


lambdatrue_mi_coverage=(lambdatrue_mean_low95_vec < true)*(lambdatrue_mean_up95_vec > true);
mean(lambdatrue_mi_coverage);

# variance estimates;
mean(lambdatrue_mi_var);
var(lambdatrue_mi_mean);

# imputation using bootstrap;

boot_mi_mean=rep(NA, cycle_no);
boot_mi_var=rep(NA, cycle_no);
boot_mi_df=rep(NA, cycle_no);
boot_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(boot_com_correct_slope_mean_mat[i,], boot_com_correct_slope_var_mat[i,]);
boot_mi_mean[i]=summary$qbar;
boot_mi_var[i]=summary$t;
boot_mi_df[i]=summary$df;
boot_mi_f[i]=summary$f;
}

boot_mi_true=mean(boot_mi_mean);
boot_mi_bias=boot_mi_true-true;

boot_mi_mse=mean((boot_mi_mean-true)^2);

# coverage;
boot_mean_low95_vec=boot_mi_mean-qt(.975, boot_mi_df)*sqrt(boot_mi_var);
boot_mean_up95_vec=boot_mi_mean+qt(.975, boot_mi_df)*sqrt(boot_mi_var);

mean_boot_length=mean(boot_mean_up95_vec-boot_mean_low95_vec);


boot_mi_coverage=(boot_mean_low95_vec < true)*(boot_mean_up95_vec > true);
mean(boot_mi_coverage);

# variance estimates;
mean(boot_mi_var);
var(boot_mi_mean);

# imputation using a fully bayesian;

bayes_mi_mean=rep(NA, cycle_no);
bayes_mi_var=rep(NA, cycle_no);
bayes_mi_df=rep(NA, cycle_no);
bayes_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(bayes_com_correct_slope_mean_mat[i,], bayes_com_correct_slope_var_mat[i,]);
bayes_mi_mean[i]=summary$qbar;
bayes_mi_var[i]=summary$t;
bayes_mi_df[i]=summary$df;
bayes_mi_f[i]=summary$f;
}

bayes_mi_true=mean(bayes_mi_mean);
bayes_mi_bias=bayes_mi_true-true;

bayes_mi_mse=mean((bayes_mi_mean-true)^2);

# coverage;
bayes_mean_low95_vec=bayes_mi_mean-qt(.975, bayes_mi_df)*sqrt(bayes_mi_var);
bayes_mean_up95_vec=bayes_mi_mean+qt(.975, bayes_mi_df)*sqrt(bayes_mi_var);

mean_bayes_length=mean(bayes_mean_up95_vec-bayes_mean_low95_vec);


bayes_mi_coverage=(bayes_mean_low95_vec < true)*(bayes_mean_up95_vec > true);
mean(bayes_mi_coverage);

# variance estimates;
mean(bayes_mi_var);
var(bayes_mi_mean);

# correlation estimates;
# complete data
corr_y2_x_vec_tran=1/2*log((1+corr_y2_x_vec)/(1-corr_y2_x_vec));

true=mean(corr_y2_x_vec_tran);
true_mse=mean((corr_y2_x_vec_tran-true)^2);

# lower and upper 95% CI;
com_corr_low95_vec_tran=corr_y2_x_vec_tran-1.96*sqrt(1/(rowobs-3));
com_corr_up95_vec_tran=corr_y2_x_vec_tran+1.96*sqrt(1/(rowobs-3));

mean_length=mean(com_corr_up95_vec_tran-com_corr_low95_vec_tran);

coverage=(com_corr_low95_vec_tran < true)*(com_corr_up95_vec_tran > true);

mean(coverage);

# CC method;
cc_corr_y2_x_vec_tran=1/2*log((1+cc_corr_y2_x_vec)/(1-cc_corr_y2_x_vec));

cc_bias=mean(cc_corr_y2_x_vec_tran)-true;
cc_mse=mean((cc_corr_y2_x_vec_tran-true)^2);

# lower and upper 95% CI;
cc_corr_low95_vec_tran=cc_corr_y2_x_vec_tran-1.96*sqrt(1/(obsno_vec-3));
cc_corr_up95_vec_tran=cc_corr_y2_x_vec_tran+1.96*sqrt(1/(obsno_vec-3));

mean_cc_length=mean(cc_corr_up95_vec_tran-cc_corr_low95_vec_tran);

cc_coverage=(cc_corr_low95_vec_tran < true)*(cc_corr_up95_vec_tran > true);

mean(cc_coverage);

# observed data transformation imputation;

obs_mi_mean=rep(NA, cycle_no);
obs_mi_var=rep(NA, cycle_no);
obs_mi_df=rep(NA, cycle_no);
obs_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(1/2*log((1+obs_corr_mat[i,])/(1-obs_corr_mat[i,])), 1/(rowobs-3));
obs_mi_mean[i]=summary$qbar;
obs_mi_var[i]=summary$t;
obs_mi_df[i]=summary$df;
obs_mi_f[i]=summary$f;
}

obs_mi_true=mean(obs_mi_mean);
obs_mi_bias=obs_mi_true-true;

obs_mi_mse=mean((obs_mi_mean-true)^2);

# coverage;
obs_mean_low95_vec=obs_mi_mean-qt(.975, obs_mi_df)*sqrt(obs_mi_var);
obs_mean_up95_vec=obs_mi_mean+qt(.975, obs_mi_df)*sqrt(obs_mi_var);

mean_obs_length=mean(obs_mean_up95_vec-obs_mean_low95_vec);


obs_mi_coverage=(obs_mean_low95_vec < true)*(obs_mean_up95_vec > true);
mean(obs_mi_coverage);

mean(obs_mi_var);
var(obs_mi_mean);

# draw transformation imputation;

draw_mi_mean=rep(NA, cycle_no);
draw_mi_var=rep(NA, cycle_no);
draw_mi_df=rep(NA, cycle_no);
draw_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(1/2*log((1+draw_corr_mat[i,])/(1-draw_corr_mat[i,])), 1/(rowobs-3));
draw_mi_mean[i]=summary$qbar;
draw_mi_var[i]=summary$t;
draw_mi_df[i]=summary$df;
draw_mi_f[i]=summary$f;
}

draw_mi_true=mean(draw_mi_mean);
draw_mi_bias=draw_mi_true-true;

draw_mi_mse=mean((draw_mi_mean-true)^2);

# coverage;
draw_mean_low95_vec=draw_mi_mean-qt(.975, draw_mi_df)*sqrt(draw_mi_var);
draw_mean_up95_vec=draw_mi_mean+qt(.975, draw_mi_df)*sqrt(draw_mi_var);

mean_draw_length=mean(draw_mean_up95_vec-draw_mean_low95_vec);


draw_mi_coverage=(draw_mean_low95_vec < true)*(draw_mean_up95_vec > true);
mean(draw_mi_coverage);

mean(draw_mi_var);
var(draw_mi_mean);

# imputation using true lambda;

lambdatrue_mi_mean=rep(NA, cycle_no);
lambdatrue_mi_var=rep(NA, cycle_no);
lambdatrue_mi_df=rep(NA, cycle_no);
lambdatrue_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(1/2*log((1+lambdatrue_corr_mat[i,])/(1-lambdatrue_corr_mat[i,])), 1/(rowobs-3));
lambdatrue_mi_mean[i]=summary$qbar;
lambdatrue_mi_var[i]=summary$t;
lambdatrue_mi_df[i]=summary$df;
lambdatrue_mi_f[i]=summary$f;
}

lambdatrue_mi_true=mean(lambdatrue_mi_mean);
lambdatrue_mi_bias=lambdatrue_mi_true-true;

lambdatrue_mi_mse=mean((lambdatrue_mi_mean-true)^2);

# coverage;
lambdatrue_mean_low95_vec=lambdatrue_mi_mean-qt(.975, lambdatrue_mi_df)*sqrt(lambdatrue_mi_var);
lambdatrue_mean_up95_vec=lambdatrue_mi_mean+qt(.975, lambdatrue_mi_df)*sqrt(lambdatrue_mi_var);

mean_lambdatrue_length=mean(lambdatrue_mean_up95_vec-lambdatrue_mean_low95_vec);


lambdatrue_mi_coverage=(lambdatrue_mean_low95_vec < true)*(lambdatrue_mean_up95_vec > true);
mean(lambdatrue_mi_coverage);

# variance estimates;
mean(lambdatrue_mi_var);
var(lambdatrue_mi_mean);

# imputation using bootstrap;

boot_mi_mean=rep(NA, cycle_no);
boot_mi_var=rep(NA, cycle_no);
boot_mi_df=rep(NA, cycle_no);
boot_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(1/2*log((1+boot_corr_mat[i,])/(1-boot_corr_mat[i,])), 1/(rowobs-3));
boot_mi_mean[i]=summary$qbar;
boot_mi_var[i]=summary$t;
boot_mi_df[i]=summary$df;
boot_mi_f[i]=summary$f;
}

boot_mi_true=mean(boot_mi_mean);
boot_mi_bias=boot_mi_true-true;

boot_mi_mse=mean((boot_mi_mean-true)^2);

# coverage;
boot_mean_low95_vec=boot_mi_mean-qt(.975, boot_mi_df)*sqrt(boot_mi_var);
boot_mean_up95_vec=boot_mi_mean+qt(.975, boot_mi_df)*sqrt(boot_mi_var);

mean_boot_length=mean(boot_mean_up95_vec-boot_mean_low95_vec);


boot_mi_coverage=(boot_mean_low95_vec < true)*(boot_mean_up95_vec > true);
mean(boot_mi_coverage);

# variance estimates;
mean(boot_mi_var);
var(boot_mi_mean);

# imputation using a fully bayesian;

bayes_mi_mean=rep(NA, cycle_no);
bayes_mi_var=rep(NA, cycle_no);
bayes_mi_df=rep(NA, cycle_no);
bayes_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(1/2*log((1+bayes_corr_mat[i,])/(1-bayes_corr_mat[i,])), 1/(rowobs-3));
bayes_mi_mean[i]=summary$qbar;
bayes_mi_var[i]=summary$t;
bayes_mi_df[i]=summary$df;
bayes_mi_f[i]=summary$f;
}

bayes_mi_true=mean(bayes_mi_mean);
bayes_mi_bias=bayes_mi_true-true;

bayes_mi_mse=mean((bayes_mi_mean-true)^2);

# coverage;
bayes_mean_low95_vec=bayes_mi_mean-qt(.975, bayes_mi_df)*sqrt(bayes_mi_var);
bayes_mean_up95_vec=bayes_mi_mean+qt(.975, bayes_mi_df)*sqrt(bayes_mi_var);

mean_bayes_length=mean(bayes_mean_up95_vec-bayes_mean_low95_vec);


bayes_mi_coverage=(bayes_mean_low95_vec < true)*(bayes_mean_up95_vec > true);
mean(bayes_mi_coverage);

# variance estimates;
mean(bayes_mi_var);
var(bayes_mi_mean);

