# Example 4.3 

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


# function bootstrap() provided the sample index of a nonparameteric sample;
bootstrap=function(n)
{
 rowobs=n;  
# donor_vector=ceiling(runif(rowobs)*rowobs);
donor_vector=sample(seq(1,n,1), replace=TRUE); 
donor_vector;
}



# set up the random seed;
set.seed(197789);


# complete-data sample size;
rowobs=1000;


# generate data;
# mean and variance of covariate X;
mu_x=1;
var_x=1/4;

# regression parameters;
# both beta0 and beta1 are fixed at 
beta0=-2;
# beta1=2; # beta1=1 creates around 30% of 1's;
# beta1=1; # beta1=2 creates around 50% of 1's;
# beta1=0.5;
# beta1=-1; # beta1=-1 creates around 7% of 1's;
# beta1=-2; # beta1=-2 creates around 6.7% of 1's;
beta1=-3; # beta1=-3 creates around 7% of 1's;

# If beta1 is less than -3, say -4 then more warning messages would appear

# number of simulations;
cycle_no=1000;

# number of multiple imputations;
mi_no=20;

# number of MCMC iterations for the Bayesian imputation;
draw_no=1000;

# vectors or matrices holding the parameter estimates;
# For marginal mean estimates;
# Before-deletion;
BD_mean_vec=BD_var_vec=rep(NA, cycle_no);
# Complete-case analysis;
CC_mean_vec=CC_var_vec=rep(NA, cycle_no);

# Bayes logistic regression imputation;
bayes_mean_mat=bayes_var_mat=matrix(NA, cycle_no, mi_no);

# MLE logistic regression imputation;
mle_mean_mat=mle_var_mat=matrix(NA, cycle_no, mi_no);

# Bootstrap logistic regression imputation;
boot_mean_mat=boot_var_mat=matrix(NA, cycle_no, mi_no);

# For the logistic regression slope coefficients;
# Before-deletion;
BD_slope_vec=BD_slope_var_vec=rep(NA, cycle_no);
# Complete-case analysis;
CC_slope_vec=CC_slope_var_vec=rep(NA, cycle_no);

# Bayes logistic regression imputation;
bayes_slope_mat=bayes_slope_var_mat=matrix(NA, cycle_no, mi_no);

# MLE logistic regression imputation;
mle_slope_mat=mle_slope_var_mat=matrix(NA, cycle_no, mi_no);

# Bootstrap logistic regression imputation;
boot_slope_mat=boot_slope_var_mat=matrix(NA, cycle_no, mi_no);


# begin the simulation of multiple imputation;

for (cycle in 1:cycle_no)

{
set.seed(cycle);

# Covariate X is a normal distribution;
x=mu_x+sqrt(var_x)*rnorm(rowobs);

# generate the logistic regression outcome y;
y=rbinom(n=rowobs, size=1, prob=exp(beta0+beta1*x)/(1+exp(beta0+beta1*x)));

# set up the missing data as MCAR on x;

miss_indi=runif(n=rowobs)<0.20;
y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
x_obs=x[miss_indi==0];
x_miss=x[miss_indi==1];

# observed data is cbind(x, y_miss);

# before-deletion data analysis;
# sample mean;
BD_mean=mean(y);
BD_var=var(y)/rowobs;
BD_mean_vec[cycle]=BD_mean;
BD_var_vec[cycle]=BD_var;

# logistic regression coefficient for slope;
BD_logistic=summary(glm(y~x, family=binomial(link="logit")));
BD_slope_vec[cycle]=BD_logistic$coeff[2,1];
BD_slope_var_vec[cycle]=BD_logistic$coeff[2,2]^2;


# different missing data methods;

# complete-case analysis;
# mean estimands;
CC_mean=mean(y_miss, na.rm="T");
CC_mean_vec[cycle]=CC_mean;
CC_var_vec[cycle]=var(y_miss, na.rm="T")/obs_no;

# logistic regression coefficient for slope;
# logistic regression coefficient for slope;
CC_logistic=summary(glm(y_obs~x_obs, family=binomial(link="logit")));
CC_slope_vec[cycle]=CC_logistic$coeff[2,1];
CC_slope_var_vec[cycle]=CC_logistic$coeff[2,2]^2;


# now impute the missing y's to get estimates of the marginal;
y_completed_mle=y_completed_boot=y_completed_bayes=y_miss;

for (i in 1:mi_no)
{

set.seed(i);

# bivariate model imputation; 
# MLE logistic regression imputation;
y_imputed_mle=mice.impute.logreg(y_miss, ry=as.logical(1-miss_indi), seed=i, x=as.matrix(x));

# Bootstarp imputation based;
# y_imputed_boot=mice.impute.logreg.boot(y_miss, ry=as.logical(1-miss_indi), seed=i, x=as.matrix(x));

boot_index=bootstrap(obs_no);
logistic_boot=summary(glm(y_obs[boot_index]~x_obs[boot_index], family=binomial(link="logit")));
beta0_boot=logistic_boot$coef[1];
beta1_boot=logistic_boot$coef[2];
y_imputed_boot=rbinom(n=mis_no, size=1, prob=exp(beta0_boot+beta1_boot*x_miss)/(1+exp(beta0_boot+beta1_boot*x_miss)));


# imputation model based on exact MCMC algorithm;
# diffuse prior;
# prior is proportional to 1
coef_s1x_posterior=bayesglm(y_obs~x_obs, family=binomial(link="logit"), prior.scale=Inf, prior.df=Inf);  

# cauchy prior with scale 2.5;
# coef_s1x_posterior=bayesglm(y_obs~x_obs, family=binomial(link="logit"), prior.scale=2.5, prior.df=1);  

coef_s1x_draw_collection=sim(coef_s1x_posterior, n.sims=draw_no)@coef;

coef_s1x_draw=coef_s1x_draw_collection[draw_no,];
beta0_bayes=coef_s1x_draw[1];
beta1_bayes=coef_s1x_draw[2];
y_imputed_bayes=rbinom(n=mis_no, size=1, prob=exp(beta0_bayes+beta1_bayes*x_miss)/(1+exp(beta0_bayes+beta1_bayes*x_miss)));


y_completed_mle[miss_seq]=y_imputed_mle;
y_completed_boot[miss_seq]=y_imputed_boot;
y_completed_bayes[miss_seq]=y_imputed_bayes;

# For MLE logistic regression imputation;
# marginal means;
mle_mean=mean(y_completed_mle);
mle_var=var(y_completed_mle)/rowobs;

mle_mean_mat[cycle,i]=mle_mean;
mle_var_mat[cycle,i]=mle_var;

# logistic regression coefficient for slope;
mle_logistic=summary(glm(y_completed_mle~x, family=binomial(link="logit")));
mle_slope_mat[cycle,i]=mle_logistic$coeff[2,1];
mle_slope_var_mat[cycle,i]=mle_logistic$coeff[2,2]^2;

# For bootstrap logistic regression imputation;
boot_mean=mean(y_completed_boot);
boot_mean_var=var(y_completed_boot)/rowobs;

boot_mean_mat[cycle,i]=boot_mean;
boot_var_mat[cycle,i]=boot_mean_var;

# logistic regression coefficient for slope;
boot_logistic=summary(glm(y_completed_boot~x, family=binomial(link="logit")));
boot_slope_mat[cycle,i]=boot_logistic$coeff[2,1];
boot_slope_var_mat[cycle,i]=boot_logistic$coeff[2,2]^2;


# For Bayes logistic regression imputation;
bayes_mean=mean(y_completed_bayes);
bayes_mean_var=var(y_completed_bayes)/rowobs;

bayes_mean_mat[cycle,i]=bayes_mean;
bayes_var_mat[cycle,i]=bayes_mean_var;

# logistic regression coefficient for slope;
bayes_logistic=summary(glm(y_completed_bayes~x, family=binomial(link="logit")));
bayes_slope_mat[cycle,i]=bayes_logistic$coeff[2,1];
bayes_slope_var_mat[cycle,i]=bayes_logistic$coeff[2,2]^2;



}

cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

#############################################################################
# Simulation results
# marginal mean estimand;

true_mean=mean(BD_mean_vec);
true_mean;

# before-deletion;
BD=performance(BD_mean_vec, BD_var_vec, true_mean);
BD;

# complete-case analysis;
CC=performance(CC_mean_vec, CC_var_vec, true_mean);
CC;

# Multiple imputation analysis;
# MLE logistic regression imputation;

mle_MI=mi_performance(mle_mean_mat, mle_var_mat, true_mean, rowobs);
mle_MI;

# Bayesian logistic regression imputation;
bayes_MI=mi_performance(bayes_mean_mat, bayes_var_mat, true_mean, rowobs);
bayes_MI;

# bootstrap logistic regression imputation;
boot_MI=mi_performance(boot_mean_mat, boot_var_mat, true_mean, rowobs);
boot_MI;


##########################################################################
# logistic regression slope estimand;
# Table 4.5

true_slope=mean(BD_slope_vec);
true_slope;

# before-deletion;
BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;

# complete-case analysis;
CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;

# Multiple imputation analysis;
# MLE logistic regression imputation;

mle_MI=mi_performance(mle_slope_mat, mle_slope_var_mat, true_slope, rowobs);
mle_MI;

# Bayesian logistic regression imputation;
bayes_MI=mi_performance(bayes_slope_mat, bayes_slope_var_mat, true_slope, rowobs);
bayes_MI;

# bootstrap logistic regression imputation;
boot_MI=mi_performance(boot_slope_mat, boot_slope_var_mat, true_slope, rowobs);
boot_MI;

