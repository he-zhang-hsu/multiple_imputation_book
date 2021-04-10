################## data generate ##############
# Program: boxcox_impute_expolore_simulation; ############
# generate data, and try out different missing data methods;
# run some simulations to assess the biasedness of various methods;

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
rm(list=ls());



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

bootstrap=function(n)
{
 rowobs=n;  
 donor_vector=ceiling(runif(rowobs)*rowobs);
 donor_vector;
}







rowobs=1000;

# set up the random seed;
set.seed(197789);

# generate data;
# distribution of x, the covariate;
# could be normal, or could be other distributions;
# fixed at mean 10, and variance 1;
mu_x=1;
var_x=1/4;

# regression parameters;
# both beta0 and beta1 are fixed at 
beta0=-2;
# beta1=1; # beta1=1 creates around 30% of 1's;
# beta1=2; # beta1=2 creates around 50% of 1's;
beta1=0.5;
#beta1=0; # creates around 12% of 1's;
# beta1=-1; # beta1=-1 creates around 7% of 1's;
# beta1=-2; # beta1=-2 creates around 6.7% of 1's;
# beta1=-3; # beta1=-3 creates around 7% of 1's;
# beta1=-4; # beta1=-4 creates around 1.57% of 1's;
# beta1=-5; # beta1=-5 creates 1.25% of 1's;
# beta1=-6;



cycle_no=1000;
mi_no=20;
draw_no=1000;

# matrices holding the parameter estimates;
# complete-data inferences;
y2_mean_vec=y2_var_vec=cc_mean_vec=cc_var_vec=y2_prop_vec=y2_prop_var_vec=y2_cc_prop_vec=y2_cc_prop_var_vec=rep(NA, cycle_no);

bayes_mean_mat=bayes_var_mat=mle_mean_mat=mle_var_mat=boot_mean_mat=boot_var_mat=guess_mean_mat=guess_var_mat=draw_mean_mat=draw_var_mat=lambdatrue_mean_mat=lambdatrue_var_mat=boot_mean_mat=boot_var_mat=bayes_mean_mat=bayes_var_mat=matrix(NA, cycle_no, mi_no);

BD_slope_vec=BD_slope_var_vec=rep(NA, cycle_no);
CC_slope_vec=CC_slope_var_vec=rep(NA, cycle_no);
mle_slope_mat=mle_slope_var_mat=boot_slope_mat=boot_slope_var_mat=bayes_slope_mat=bayes_slope_var_mat=matrix(NA, cycle_no, mi_no);




for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);
# x=runif(rowobs, min=0, max=2);

# generate the logistic regression outcome y2;
y2=rbinom(n=rowobs, size=1, prob=exp(beta0+beta1*x)/(1+exp(beta0+beta1*x)));



# set up the missing data as MCAR on x;


miss_indi=runif(n=rowobs)<0.20;
y2_miss=y2;
y2_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y2_miss[!is.na(y2_miss)]);
mis_no=rowobs-obs_no;
y2_obs=y2[miss_indi==0];
x_obs=x[miss_indi==0];
x_miss=x[miss_indi==1];

# observed data is cbind(x, y2_miss);

# compelete data analysis;
# sample mean;
y2_mean=mean(y2);
y2_var=var(y2)/rowobs;
y2_mean_vec[cycle]=y2_mean;
y2_var_vec[cycle]=y2_var;

# logistic regression coefficient for slope;
BD_logistic=summary(glm(y2~x, family=binomial(link="logit")));
BD_slope_vec[cycle]=BD_logistic$coeff[2,1];
BD_slope_var_vec[cycle]=BD_logistic$coeff[2,2]^2;


# different missing data methods;

# complete-case analysis;
# mean estimands;
cc_mean=mean(y2_miss, na.rm="T");
cc_mean_vec[cycle]=cc_mean;
cc_var_vec[cycle]=var(y2_miss, na.rm="T")/obs_no;

# logistic regression coefficient for slope;
# logistic regression coefficient for slope;
CC_logistic=summary(glm(y2_obs~x_obs, family=binomial(link="logit")));
CC_slope_vec[cycle]=CC_logistic$coeff[2,1];
CC_slope_var_vec[cycle]=CC_logistic$coeff[2,2]^2;




# now impute the missing y2's to get estimates of the marginal;
y2_completed_mle=y2_completed_boot=y2_completed_bayes=y2_miss;

for (i in 1:mi_no)
{

set.seed(i);

# bivariate model imputation; 
# Rubin's imputation based on MLE;
y2_imputed_mle=mice.impute.logreg(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=as.matrix(x));

# Bootstarp imputation based;
# y2_imputed_boot=mice.impute.logreg.boot(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=as.matrix(x));

boot_index=bootstrap(obs_no);
logistic_boot=summary(glm(y2_obs[boot_index]~x_obs[boot_index], family=binomial(link="logit")));
beta0_boot=logistic_boot$coef[1];
beta1_boot=logistic_boot$coef[2];
y2_imputed_boot=rbinom(n=mis_no, size=1, prob=exp(beta0_boot+beta1_boot*x_miss)/(1+exp(beta0_boot+beta1_boot*x_miss)));


# imputation model based on exact MCMC algorithm;
# diffuse prior;
coef_s1x_posterior=bayesglm(y2_obs~x_obs, family=binomial(link="logit"), prior.scale=Inf, prior.df=Inf);  

# cauchy prior with scale 2.5;
# coef_s1x_posterior=bayesglm(y2_obs~x_obs, family=binomial(link="logit"), prior.scale=2.5, prior.df=1);  

coef_s1x_draw_collection=sim(coef_s1x_posterior, n.sims=draw_no)@coef;

coef_s1x_draw=coef_s1x_draw_collection[draw_no,];
beta0_bayes=coef_s1x_draw[1];
beta1_bayes=coef_s1x_draw[2];
y2_imputed_bayes=rbinom(n=mis_no, size=1, prob=exp(beta0_bayes+beta1_bayes*x_miss)/(1+exp(beta0_bayes+beta1_bayes*x_miss)));



# univariate model imputation
# y2_imputed=mice.impute.norm(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(rep(1,rowobs),x));


y2_completed_mle[miss_seq]=y2_imputed_mle;
y2_completed_boot[miss_seq]=y2_imputed_boot;
y2_completed_bayes[miss_seq]=y2_imputed_bayes;

# marginal means;
mle_mean=mean(y2_completed_mle);
mle_var=var(y2_completed_mle)/rowobs;

mle_mean_mat[cycle,i]=mle_mean;
mle_var_mat[cycle,i]=mle_var;

# logistic regression coefficient for slope;
mle_logistic=summary(glm(y2_completed_mle~x, family=binomial(link="logit")));
mle_slope_mat[cycle,i]=mle_logistic$coeff[2,1];
mle_slope_var_mat[cycle,i]=mle_logistic$coeff[2,2]^2;



boot_mean=mean(y2_completed_boot);
boot_mean_var=var(y2_completed_boot)/rowobs;

boot_mean_mat[cycle,i]=boot_mean;
boot_var_mat[cycle,i]=boot_mean_var;



# logistic regression coefficient for slope;
boot_logistic=summary(glm(y2_completed_boot~x, family=binomial(link="logit")));
boot_slope_mat[cycle,i]=boot_logistic$coeff[2,1];
boot_slope_var_mat[cycle,i]=boot_logistic$coeff[2,2]^2;


bayes_mean=mean(y2_completed_bayes);
bayes_mean_var=var(y2_completed_bayes)/rowobs;

bayes_mean_mat[cycle,i]=bayes_mean;
bayes_var_mat[cycle,i]=bayes_mean_var;

# logistic regression coefficient for slope;
bayes_logistic=summary(glm(y2_completed_bayes~x, family=binomial(link="logit")));
bayes_slope_mat[cycle,i]=bayes_logistic$coeff[2,1];
bayes_slope_var_mat[cycle,i]=bayes_logistic$coeff[2,2]^2;



}

cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

# plot the distribution of the data;

# check the performance of the estimates;
# population quantity
# mean estimand;
# true=pnorm((beta0+beta1*mu_x-cut_off)/(beta1^2*var_x+var_error));

#############################################################################
# marginal mean estimand;

true=mean(y2_mean_vec);
true;

mean_bias=mean(y2_mean_vec-true);
mean_bias;

mean_bias/true;

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
cc_mean_low95_vec=cc_mean_vec-1.96*sqrt(cc_var_vec);
cc_mean_up95_vec=cc_mean_vec+1.96*sqrt(cc_var_vec);

mean_cc_length=mean(cc_mean_up95_vec-cc_mean_low95_vec);
mean_cc_length;

cc_coverage=(cc_mean_low95_vec < true)*(cc_mean_up95_vec > true);

mean(cc_coverage);

# variance estimates;
mean(sqrt(cc_var_vec));
sqrt(var(cc_mean_vec));


# multiple imputation estimates;
# marginal transformation imputation;

# observed data transformation imputation;

mle_mi_mean=rep(NA, cycle_no);
mle_mi_var=rep(NA, cycle_no);
mle_mi_df=rep(NA, cycle_no);
mle_mi_f=rep(NA, cycle_no);


boot_mi_mean=rep(NA, cycle_no);
boot_mi_var=rep(NA, cycle_no);
boot_mi_df=rep(NA, cycle_no);
boot_mi_f=rep(NA, cycle_no);


bayes_mi_mean=rep(NA, cycle_no);
bayes_mi_var=rep(NA, cycle_no);
bayes_mi_df=rep(NA, cycle_no);
bayes_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
mle_summary=pool.scalar(mle_mean_mat[i,], mle_var_mat[i,], n=rowobs);
mle_mi_mean[i]=mle_summary$qbar;
mle_mi_var[i]=mle_summary$t;
mle_mi_df[i]=mle_summary$df;
mle_mi_f[i]=mle_summary$f;

boot_summary=pool.scalar(boot_mean_mat[i,], boot_var_mat[i,], n=rowobs);
boot_mi_mean[i]=boot_summary$qbar;
boot_mi_var[i]=boot_summary$t;
boot_mi_df[i]=boot_summary$df;
boot_mi_f[i]=boot_summary$f;

bayes_summary=pool.scalar(bayes_mean_mat[i,], bayes_var_mat[i,], n=rowobs);
bayes_mi_mean[i]=bayes_summary$qbar;
bayes_mi_var[i]=bayes_summary$t;
bayes_mi_df[i]=bayes_summary$df;
bayes_mi_f[i]=bayes_summary$f;




}

# mle imputation;
mle_mi_true=mean(mle_mi_mean);
mle_mi_bias=mle_mi_true-true;

mle_mi_bias;

mle_mi_bias/true;

mle_mi_mse=mean((mle_mi_mean-true)^2);
mle_mi_mse;

# coverage;
mle_mean_low95_vec=mle_mi_mean-qt(.975, mle_mi_df)*sqrt(mle_mi_var);
mle_mean_up95_vec=mle_mi_mean+qt(.975, mle_mi_df)*sqrt(mle_mi_var);

mean_mle_length=mean(mle_mean_up95_vec-mle_mean_low95_vec);
mean_mle_length;

mle_mi_coverage=(mle_mean_low95_vec < true)*(mle_mean_up95_vec > true);
mean(mle_mi_coverage);

# variance estimates;
mean(sqrt(mle_mi_var));
sqrt(var(mle_mi_mean));


# bootstrap imputation;
boot_mi_mean=boot_mi_mean[!is.na(boot_mi_f)];
boot_mi_var=boot_mi_var[!is.na(boot_mi_f)];
boot_mi_df=boot_mi_df[!is.na(boot_mi_f)];
boot_mi_f=boot_mi_f[!is.na(boot_mi_f)];

boot_mi_true=mean(boot_mi_mean);
boot_mi_bias=boot_mi_true-true;

boot_mi_bias;

boot_mi_bias/true;

boot_mi_mse=mean((boot_mi_mean-true)^2);
boot_mi_mse;

# coverage;
boot_mean_low95_vec=boot_mi_mean-qt(.975, boot_mi_df)*sqrt(boot_mi_var);
boot_mean_up95_vec=boot_mi_mean+qt(.975, boot_mi_df)*sqrt(boot_mi_var);

# trans_mean_low95_vec=trans_mean_low95_vec[!is.na(trans_mean_low95_vec)];
# trans_mean_up95_vec=trans_mean_up95_vec[!is.na(trans_mean_up95_vec)];


mean_boot_length=mean(boot_mean_up95_vec-boot_mean_low95_vec);
mean_boot_length;


boot_mi_coverage=(boot_mean_low95_vec < true)*(boot_mean_up95_vec > true);
mean(boot_mi_coverage);

# variance estimates;
mean(sqrt(boot_mi_var));
sqrt(var(boot_mi_mean));

# bayes imputation;
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



##########################################################################
# slope estimand;

# complete-data analysis;

true_slope=mean(BD_slope_vec);
true_slope;


# performance;
BD_slope_mse=mean((BD_slope_vec-true_slope)^2);
BD_slope_mse;

# lower and upper 95% CI;
BD_slope_low95_vec=BD_slope_vec-1.96*sqrt(BD_slope_var_vec);
BD_slope_up95_vec=BD_slope_vec+1.96*sqrt(BD_slope_var_vec);

BD_slope_length=mean(BD_slope_up95_vec-BD_slope_low95_vec);
BD_slope_length;

BD_slope_cov=(BD_slope_low95_vec < true_slope)*(BD_slope_up95_vec > true_slope);
mean(BD_slope_cov);

# variance estimates;
mean(sqrt(BD_slope_var_vec));
sqrt(var(BD_slope_vec));


# complete-case analysis;
CC_slope_bias=mean(CC_slope_vec)-true_slope;
CC_slope_bias;

CC_slope_bias/true_slope;

CC_slope_mse=mean((CC_slope_vec-true_slope)^2);
CC_slope_mse;

# lower and upper 95% CI;
CC_slope_low95_vec=CC_slope_vec-1.96*sqrt(CC_slope_var_vec);
CC_slope_up95_vec=CC_slope_vec+1.96*sqrt(CC_slope_var_vec);

CC_slope_length=mean(CC_slope_up95_vec-CC_slope_low95_vec);
CC_slope_length;

CC_slope_cov=(CC_slope_low95_vec < true_slope)*(CC_slope_up95_vec > true_slope);

mean(CC_slope_cov);

# variance estimates;
mean(sqrt(CC_slope_var_vec));
sqrt(var(CC_slope_vec));


# multiple imputation estimates;

mle_mi_mean=rep(NA, cycle_no);
mle_mi_var=rep(NA, cycle_no);
mle_mi_df=rep(NA, cycle_no);
mle_mi_f=rep(NA, cycle_no);

boot_mi_mean=rep(NA, cycle_no);
boot_mi_var=rep(NA, cycle_no);
boot_mi_df=rep(NA, cycle_no);
boot_mi_f=rep(NA, cycle_no);

bayes_mi_mean=rep(NA, cycle_no);
bayes_mi_var=rep(NA, cycle_no);
bayes_mi_df=rep(NA, cycle_no);
bayes_mi_f=rep(NA, cycle_no);



for (i in 1:cycle_no)
{
mle_summary=pool.scalar(mle_slope_mat[i,], mle_slope_var_mat[i,], n=rowobs);
mle_mi_mean[i]=mle_summary$qbar;
mle_mi_var[i]=mle_summary$t;
mle_mi_df[i]=mle_summary$df;
mle_mi_f[i]=mle_summary$f;

boot_summary=pool.scalar(boot_slope_mat[i,], boot_slope_var_mat[i,], n=rowobs);
boot_mi_mean[i]=boot_summary$qbar;
boot_mi_var[i]=boot_summary$t;
boot_mi_df[i]=boot_summary$df;
boot_mi_f[i]=boot_summary$f;

bayes_summary=pool.scalar(bayes_slope_mat[i,], bayes_slope_var_mat[i,], n=rowobs);
bayes_mi_mean[i]=bayes_summary$qbar;
bayes_mi_var[i]=bayes_summary$t;
bayes_mi_df[i]=bayes_summary$df;
bayes_mi_f[i]=bayes_summary$f;



}



# For slope;

# proper imputation;

mle_mi_bias=mean(mle_mi_mean)-true_slope;
mle_mi_bias;

mle_mi_bias/true_slope;

mle_mi_mse=mean((mle_mi_mean-true_slope)^2);
mle_mi_mse;

# coverage;
mle_slope_low95_vec=mle_mi_mean-qt(.975, mle_mi_df)*sqrt(mle_mi_var);
mle_slope_up95_vec=mle_mi_mean+qt(.975, mle_mi_df)*sqrt(mle_mi_var);

mle_slope_length=mean(mle_slope_up95_vec-mle_slope_low95_vec);
mle_slope_length;


mle_slope_coverage=(mle_slope_low95_vec < true_slope)*(mle_slope_up95_vec > true_slope);

mean(mle_slope_coverage);

# variance estimates;
mean(sqrt(mle_mi_var));
sqrt(var(mle_mi_mean));

boot_mi_bias=mean(boot_mi_mean)-true_slope;
boot_mi_bias;

boot_mi_bias/true_slope;

boot_mi_mse=mean((boot_mi_mean-true_slope)^2);
boot_mi_mse;

# coverage;
boot_slope_low95_vec=boot_mi_mean-qt(.975, boot_mi_df)*sqrt(boot_mi_var);
boot_slope_up95_vec=boot_mi_mean+qt(.975, boot_mi_df)*sqrt(boot_mi_var);

boot_slope_length=mean(boot_slope_up95_vec-boot_slope_low95_vec);
boot_slope_length;

boot_slope_coverage=(boot_slope_low95_vec < true_slope)*(boot_slope_up95_vec > true_slope);

mean(boot_slope_coverage);

# variance estimates;
mean(sqrt(boot_mi_var));
sqrt(var(boot_mi_mean));

bayes_mi_bias=mean(bayes_mi_mean)-true_slope;
bayes_mi_bias;

bayes_mi_bias/true_slope;

bayes_mi_mse=mean((bayes_mi_mean-true_slope)^2);
bayes_mi_mse;

# coverage;
bayes_slope_low95_vec=bayes_mi_mean-qt(.975, bayes_mi_df)*sqrt(bayes_mi_var);
bayes_slope_up95_vec=bayes_mi_mean+qt(.975, bayes_mi_df)*sqrt(bayes_mi_var);

bayes_slope_length=mean(bayes_slope_up95_vec-bayes_slope_low95_vec);
bayes_slope_length;


bayes_slope_coverage=(bayes_slope_low95_vec < true_slope)*(bayes_slope_up95_vec > true_slope);

mean(bayes_slope_coverage);

# variance estimates;
mean(sqrt(bayes_mi_var));
sqrt(var(bayes_mi_mean));

# histogram of the posterior draws of the slope;
hist(coef_s1x_draw_collection[,2]);

# the mean of logistic model
# z=rnorm(1000000, 1, sd=0.5);
# beta1=-5;
# beta1=-4
# beta1=-3;
# beta1=-2;
# beta1=-1;
# beta1=0;
# beta1=0.5;
# beta1=1;
#beta1=2;
#y2=rbinom(n=1000000, size=1, prob=exp(beta0+beta1*z)/(1+exp(beta0+beta1*z)));
# mean(y2);




# test=MCMClogit(y2_obs~x_obs, burnin=1000, mcmc=2000, tune=0.40, beta.start=c(beta0, beta1));
# x_posterior=bayesglm(hospice_ac[true==1,2]~hospice_ac[true==1,3:n_col], family=binomial(link="logit"), prior.scale=2.5, prior.df=1);  
    # Same as M3, explicitly specifying Cauchy prior with scale 2.5 coef_s1