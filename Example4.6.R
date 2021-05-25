# Example 4.6

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

# function mice.impute.lda.proper() imputes the binary data
# using discriminant analysis model and is from R MICE
mice.impute.lda.proper <- function(y, ry, x, wy = NULL, ...) { 
if (is.null(wy)) wy <- !ry 
fy <- as.factor(y) 
nc <- length(levels(fy)) 
    
   #   SvB June 2009 - take bootstrap sample of training data 
      idx <- sample((1:length(y))[ry], size=sum(ry), replace=TRUE) 
      x[ry,] <- x[idx,] 
      y[ry] <- y[idx] 
   #   end bootstrap 
fy <- as.factor(y) 
nc <- length(levels(fy)) 
    
   fit <- lda(x, fy, subset = ry) 
   post <- predict(fit, x[wy, , drop = FALSE])$posterior 
   un <- rep(runif(sum(wy)), each = nc) 
   idx <- 1 + apply(un > apply(post, 1, cumsum), 2, sum) 
   return(levels(fy)[idx]) 
 } 


# complete-data sample size of the simulation;
rowobs=1000;

# set up the random seed;
set.seed(197789);

# number of simulations;
cycle_no=1000;

# number of multiple imputations;
mi_no=30;

# number of MCMC iterations for Bayesian logistic regression imputation;
draw_no=1000;

# matrices holding the parameter estimates;
# complete-data inferences;
# For marginal means;
# Before-deletion analysis;
BD_mean_vec=BD_var_vec=rep(NA, cycle_no);

# Complete-case analysis;
CC_mean_vec=CC_var_vec=rep(NA, cycle_no);

# For multiple imputation
# MLE logistic regression imputation without adjustment for
# data separation
draw_mean_mat=draw_var_mat=matrix(NA, cycle_no, mi_no);

# MLE logistic regression imputation with adjustment for 
# data separation
mle_mean_mat=mle_var_mat=matrix(NA, cycle_no, mi_no);

# bootstrap imputation;
boot_mean_mat=boot_var_mat=matrix(NA, cycle_no, mi_no);

# linear discriminant analysis imputation;
lda_mean_mat=lda_var_mat=matrix(NA, cycle_no, mi_no);

# Bayes logistic regression imputation;
bayes_mean_mat=bayes_var_mat=matrix(NA, cycle_no, mi_no);

# For linear regression coefficient;
# Before-deletion analysis
BD_slope_vec=BD_slope_var_vec=rep(NA, cycle_no);

# Complete-case analysis;
CC_slope_vec=CC_slope_var_vec=rep(NA, cycle_no);

# For multiple imputation;
# MLE logistic regression imputation without adjustment for
# data separation
draw_slope_mat=draw_slope_var_mat=matrix(NA, cycle_no, mi_no);

# MLE logistic regression imputation with adjustment for 
# data separation
mle_slope_mat=mle_slope_var_mat=matrix(NA, cycle_no, mi_no);

# bootstrap imputation;
boot_slope_mat=boot_slope_var_mat=matrix(NA, cycle_no, mi_no);

# linear discriminant analysis imputation;
lda_slope_mat=lda_slope_var_mat=matrix(NA, cycle_no, mi_no);

# Bayes logistic regression imputation;
bayes_slope_mat=bayes_slope_var_mat=matrix(NA, cycle_no, mi_no);


# begin the simulation;

for (cycle in 1:cycle_no)
{

set.seed(cycle);


# the following data generation is based on White et al. (2010) simulation;
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

# Before-deletion data analysis;
# sample mean;
BD_mean=mean(y);
BD_var=var(y)/rowobs;
BD_mean_vec[cycle]=BD_mean;
BD_var_vec[cycle]=BD_var;


# linear regression coefficient of z on y, x, and u;
BD_linear=summary(lm(z~0+y+x+u));
BD_slope_vec[cycle]=BD_linear$coef[1,1];
BD_slope_var_vec[cycle]=BD_linear$coef[1,2]^2;


# different missing data methods;

# complete-case analysis;
# mean estimands;
CC_mean=mean(y_miss, na.rm="T");
CC_mean_vec[cycle]=CC_mean;
CC_var_vec[cycle]=var(y_miss, na.rm="T")/obs_no;

CC_linear=summary(lm(z_obs~0+y_obs+x_obs+u_obs));
CC_slope_vec[cycle]=CC_linear$coef[1,1];
CC_slope_var_vec[cycle]=CC_linear$coef[1,2]^2;

# now impute the missing y's to get estimates of the marginal;
y_completed_draw=y_completed_mle=y_completed_boot=y_completed_lda=y_completed_bayes=y_miss;

for (i in 1:mi_no)
{
set.seed(i);

# MLE logistic regression imputation without adjustment;

mlefit=glm(y_obs~cbind(x_obs,u_obs, z_obs), family=binomial(link="logit"));
regcoef=mlefit$coef;  
covmatrix=vcov(mlefit);
regdraw=mvrnorm(n=1, mu=regcoef, Sigma=covmatrix);

beta0_draw=regdraw[1];
beta1_draw=regdraw[2];
# top code the very large draws;
if (beta1_draw > 20) beta1_draw=20;
beta2_draw=regdraw[3];
beta3_draw=regdraw[4];


y_imputed_draw=rbinom(n=mis_no, size=1, prob=exp(beta0_draw+beta1_draw*x_miss+beta2_draw*u_miss+beta3_draw*z_miss)/(1+exp(beta0_draw+beta1_draw*x_miss+beta2_draw*u_miss+beta3_draw*z_miss)));


# pseudo observation imputation;
# or MLE logistic regression imputation with adjustment for data separation;

y_imputed_mle=mice.impute.logreg(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x,u,z));

# bootstrap imputation;

boot_index=bootstrap(obs_no);
logistic_boot=summary(glm(y_obs[boot_index]~cbind(x_obs[boot_index], u_obs[boot_index], z_obs[boot_index]), family=binomial(link="logit")));
beta0_boot=logistic_boot$coef[1];
beta1_boot=logistic_boot$coef[2];
beta2_boot=logistic_boot$coef[3];
beta3_boot=logistic_boot$coef[4];

y_imputed_boot=rbinom(n=mis_no, size=1, prob=exp(beta0_boot+beta1_boot*x_miss+beta2_boot*u_miss+beta3_boot*z_miss)/(1+exp(beta0_boot+beta1_boot*x_miss+beta2_boot*u_miss+beta3_boot*z_miss)));

# discriminant imputation based;

y_imputed_lda=mice.impute.lda.proper(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x,u,z));

# imputation model based on exact MCMC algorithm;
# diffuse prior;
# coef_s1x_posterior=bayesglm(y_obs~x_obs, family=binomial(link="logit"), prior.scale=Inf, prior.df=Inf);  

# cauchy prior with scale 2.5;
# coef_s1x_posterior=bayesglm(y_obs~x_obs, family=binomial(link="logit"), prior.scale=2.5, prior.df=1);  

coef_s1x_posterior=bayesglm(y_obs~cbind(x_obs,u_obs, z_obs), family=binomial(link="logit"), prior.scale=2.5, prior.df=1);  

coef_s1x_draw_collection=sim(coef_s1x_posterior, n.sims=draw_no)@coef;

coef_s1x_draw=coef_s1x_draw_collection[draw_no,];
beta0_bayes=coef_s1x_draw[1];
beta1_bayes=coef_s1x_draw[2];
beta2_bayes=coef_s1x_draw[3];
beta3_bayes=coef_s1x_draw[4];


y_imputed_bayes=rbinom(n=mis_no, size=1, prob=exp(beta0_bayes+beta1_bayes*x_miss+beta2_bayes*u_miss+beta3_bayes*z_miss)/(1+exp(beta0_bayes+beta1_bayes*x_miss+beta2_bayes*u_miss+beta3_bayes*z_miss)));


y_completed_draw[miss_seq]=y_imputed_draw;
y_completed_mle[miss_seq]=y_imputed_mle;
y_completed_boot[miss_seq]=y_imputed_boot;
y_completed_lda[miss_seq]=as.numeric(y_imputed_lda);
y_completed_bayes[miss_seq]=y_imputed_bayes;

# Estimates;
# MLE logistic regression imptuation without any adjustment;

draw_mean=mean(y_completed_draw);
draw_var=var(y_completed_draw)/rowobs;

draw_mean_mat[cycle,i]=draw_mean;
draw_var_mat[cycle,i]=draw_var;

draw_linear=summary(lm(z~0+y_completed_draw+x+u));
draw_slope_mat[cycle,i]=draw_linear$coef[1,1];
draw_slope_var_mat[cycle,i]=draw_linear$coef[1,2]^2;

# pseudo observation imputation;
mle_mean=mean(y_completed_mle);
mle_var=var(y_completed_mle)/rowobs;

mle_mean_mat[cycle,i]=mle_mean;
mle_var_mat[cycle,i]=mle_var;

mle_linear=summary(lm(z~0+y_completed_mle+x+u));
mle_slope_mat[cycle,i]=mle_linear$coef[1,1];
mle_slope_var_mat[cycle,i]=mle_linear$coef[1,2]^2;

# bootstrap imputation;

boot_mean=mean(y_completed_boot);
boot_var=var(y_completed_boot)/rowobs;

boot_mean_mat[cycle,i]=boot_mean;
boot_var_mat[cycle,i]=boot_var;

boot_linear=summary(lm(z~0+y_completed_boot+x+u));
boot_slope_mat[cycle,i]=boot_linear$coef[1,1];
boot_slope_var_mat[cycle,i]=boot_linear$coef[1,2]^2;

# discriminant analysis imputation;

lda_mean=mean(y_completed_lda);
lda_mean_var=var(y_completed_lda)/rowobs;

lda_mean_mat[cycle,i]=lda_mean;
lda_var_mat[cycle,i]=lda_mean_var;

lda_linear=summary(lm(z~0+y_completed_lda+x+u));
lda_slope_mat[cycle,i]=lda_linear$coef[1,1];
lda_slope_var_mat[cycle,i]=lda_linear$coef[1,2]^2;

# Bayesian imputation;

bayes_mean=mean(y_completed_bayes);
bayes_mean_var=var(y_completed_bayes)/rowobs;

bayes_mean_mat[cycle,i]=bayes_mean;
bayes_var_mat[cycle,i]=bayes_mean_var;

bayes_linear=summary(lm(z~0+y_completed_bayes+x+u));
bayes_slope_mat[cycle,i]=bayes_linear$coef[1,1];
bayes_slope_var_mat[cycle,i]=bayes_linear$coef[1,2]^2;

}

cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

#############################################################################
# Simulation results;
# Mean estimand;
# For Table 4.9

true_mean=mean(BD_mean_vec);
true_mean;

# before-deletion;
BD=performance(BD_mean_vec, BD_var_vec, true_mean);
BD;

# complete-case analysis;
CC=performance(CC_mean_vec, CC_var_vec, true_mean);
CC;

# Multiple imputation analysis;
# MLE logistic regression imputation without any adjustment;

draw_MI=mi_performance(draw_mean_mat, draw_var_mat, true_mean, rowobs);
draw_MI;

# Pseudo observation imputation;
mle_MI=mi_performance(mle_mean_mat, mle_var_mat, true_mean, rowobs);
mle_MI;

# discriminant analysis imputation;
lda_MI=mi_performance(lda_mean_mat, lda_var_mat, true_mean, rowobs);
lda_MI;

# bootstrap imputation;
boot_MI=mi_performance(boot_mean_mat, boot_var_mat, true_mean, rowobs);
boot_MI;

# Bayes imputation;
bayes_MI=mi_performance(bayes_mean_mat, bayes_var_mat, true_mean, rowobs);
bayes_MI;




##########################################################################
# slope estimand;
# For Table 4.9
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
# MLE logistic regression imputation without any adjustment;
draw_MI=mi_performance(draw_slope_mat, draw_slope_var_mat, true_slope, rowobs);
draw_MI;

# Pseudo observation imputation;
mle_MI=mi_performance(mle_slope_mat, mle_slope_var_mat, true_slope, rowobs);
mle_MI;


# discriminant analysis imputation;
lda_MI=mi_performance(lda_slope_mat, lda_slope_var_mat, true_slope, rowobs);
lda_MI;

# bootstrap imputation;
boot_MI=mi_performance(boot_slope_mat, boot_slope_var_mat, true_slope, rowobs);
boot_MI;

# Bayes imputation;
bayes_MI=mi_performance(bayes_slope_mat, bayes_slope_var_mat, true_slope, rowobs);
bayes_MI;

# make plots usin the estimates draws from the last simulation;
# this is from the Bayes logistic regression imputation method;

hist(coef_s1x_draw_collection[,2], main="Hisogram of posterior draws", xlab="gamma1");

