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






rowobs=1000;

# set up the random seed;
set.seed(197789);



cycle_no=1000;
mi_no=30;
draw_no=1000;

# matrices holding the parameter estimates;
# complete-data inferences;
y2_mean_vec=y2_var_vec=cc_mean_vec=cc_var_vec=rep(NA, cycle_no);

draw_mean_mat=draw_var_mat=mle_mean_mat=mle_var_mat=boot_mean_mat=boot_var_mat=lda_mean_mat=lda_var_mat=bayes_mean_mat=bayes_var_mat=matrix(NA, cycle_no, mi_no);

BD_slope_vec=BD_slope_var_vec=rep(NA, cycle_no);
CC_slope_vec=CC_slope_var_vec=rep(NA, cycle_no);
draw_slope_mat=draw_slope_var_mat=mle_slope_mat=mle_slope_var_mat=boot_slope_mat=boot_slope_var_mat=lda_slope_mat=lda_slope_var_mat=bayes_slope_mat=bayes_slope_var_mat=matrix(NA, cycle_no, mi_no);




for (cycle in 1:cycle_no)

{
set.seed(cycle);

# x=mu_x+sqrt(var_x)*rnorm(rowobs);
# x=runif(rowobs, min=0, max=2);

# generate the logistic regression outcome y2;
# y2=rbinom(n=rowobs, size=1, prob=exp(beta0+beta1*x)/(1+exp(beta0+beta1*x)));


# the following data generation is based on White et al. (2010) simulation;
y2=rbinom(n=rowobs, size=1, prob=0.1);
zeros_y2=rowobs-sum(y2); # number of y2s are 0;
x=rep(0,rowobs);
x[y2==0]=rbinom(n=zeros_y2, size=1, prob=0.8);
u=rnorm(rowobs);
z=y2+x+0.5*u+2*rnorm(rowobs);


# set up the missing data as MCAR on x;


miss_indi=runif(n=rowobs)<0.30;
y2_miss=y2;
y2_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y2_miss[!is.na(y2_miss)]);
mis_no=rowobs-obs_no;
y2_obs=y2[miss_indi==0];
x_obs=x[miss_indi==0];
x_miss=x[miss_indi==1];


u_obs=u[miss_indi==0];
u_miss=u[miss_indi==1];
z_obs=z[miss_indi==0];
z_miss=z[miss_indi==1];

# observed data is cbind(x, y2_miss);

# compelete data analysis;
# sample mean;
y2_mean=mean(y2);
y2_var=var(y2)/rowobs;
y2_mean_vec[cycle]=y2_mean;
y2_var_vec[cycle]=y2_var;

# logistic regression coefficient for slope;
# BD_logistic=summary(glm(y2~x, family=binomial(link="logit")));
# BD_slope_vec[cycle]=BD_logistic$coeff[2,1];
# BD_slope_var_vec[cycle]=BD_logistic$coeff[2,2]^2;

# linear regression coefficient of z on y2, x, and u;
BD_logistic=summary(lm(z~0+y2+x+u));
BD_slope_vec[cycle]=BD_logistic$coef[1,1];
BD_slope_var_vec[cycle]=BD_logistic$coef[1,2]^2;


# different missing data methods;

# complete-case analysis;
# mean estimands;
cc_mean=mean(y2_miss, na.rm="T");
cc_mean_vec[cycle]=cc_mean;
cc_var_vec[cycle]=var(y2_miss, na.rm="T")/obs_no;

# logistic regression coefficient for slope;
# logistic regression coefficient for slope;
# CC_logistic=summary(glm(y2_obs~x_obs, family=binomial(link="logit")));
# CC_slope_vec[cycle]=CC_logistic$coeff[2,1];
# CC_slope_var_vec[cycle]=CC_logistic$coeff[2,2]^2;

CC_logistic=summary(lm(z_obs~0+y2_obs+x_obs+u_obs));
CC_slope_vec[cycle]=CC_logistic$coef[1,1];
CC_slope_var_vec[cycle]=CC_logistic$coef[1,2]^2;




# now impute the missing y2's to get estimates of the marginal;
y2_completed_draw=y2_completed_mle=y2_completed_boot=y2_completed_lda=y2_completed_bayes=y2_miss;

for (i in 1:mi_no)
{
set.seed(i);

# bivariate model imputation; 
# Rubin's imputation based on MLE;
# y2_imputed_mle=mice.impute.logreg(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=as.matrix(x));

# original normal approximation imputation;

mlefit=glm(y2_obs~cbind(x_obs,u_obs, z_obs), family=binomial(link="logit"));
regcoef=mlefit$coef;  
covmatrix=vcov(mlefit);
regdraw=mvrnorm(n=1, mu=regcoef, Sigma=covmatrix);

beta0_draw=regdraw[1];
beta1_draw=regdraw[2];
if (beta1_draw > 20) beta1_draw=20;
beta2_draw=regdraw[3];
beta3_draw=regdraw[4];

# test=exp(beta0_draw+beta1_draw*x_miss+beta2_draw*u_miss+beta3_draw*z_miss)/(1+exp(beta0_draw+beta1_draw*x_miss+beta2_draw*u_miss+beta3_draw*z_miss));



y2_imputed_draw=rbinom(n=mis_no, size=1, prob=exp(beta0_draw+beta1_draw*x_miss+beta2_draw*u_miss+beta3_draw*z_miss)/(1+exp(beta0_draw+beta1_draw*x_miss+beta2_draw*u_miss+beta3_draw*z_miss)));


# pseudo observation imputation;

y2_imputed_mle=mice.impute.logreg(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x,u,z));

# bootstrap imputation;

boot_index=bootstrap(obs_no);
logistic_boot=summary(glm(y2_obs[boot_index]~cbind(x_obs[boot_index], u_obs[boot_index], z_obs[boot_index]), family=binomial(link="logit")));
beta0_boot=logistic_boot$coef[1];
beta1_boot=logistic_boot$coef[2];
beta2_boot=logistic_boot$coef[3];
beta3_boot=logistic_boot$coef[4];

y2_imputed_boot=rbinom(n=mis_no, size=1, prob=exp(beta0_boot+beta1_boot*x_miss+beta2_boot*u_miss+beta3_boot*z_miss)/(1+exp(beta0_boot+beta1_boot*x_miss+beta2_boot*u_miss+beta3_boot*z_miss)));



# discriminant imputation based;

y2_imputed_lda=mice.impute.lda.proper(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x,u,z));



# imputation model based on exact MCMC algorithm;
# diffuse prior;
# coef_s1x_posterior=bayesglm(y2_obs~x_obs, family=binomial(link="logit"), prior.scale=Inf, prior.df=Inf);  

# cauchy prior with scale 2.5;
# coef_s1x_posterior=bayesglm(y2_obs~x_obs, family=binomial(link="logit"), prior.scale=2.5, prior.df=1);  

coef_s1x_posterior=bayesglm(y2_obs~cbind(x_obs,u_obs, z_obs), family=binomial(link="logit"), prior.scale=2.5, prior.df=1);  

coef_s1x_draw_collection=sim(coef_s1x_posterior, n.sims=draw_no)@coef;

coef_s1x_draw=coef_s1x_draw_collection[draw_no,];
beta0_bayes=coef_s1x_draw[1];
beta1_bayes=coef_s1x_draw[2];
beta2_bayes=coef_s1x_draw[3];
beta3_bayes=coef_s1x_draw[4];


y2_imputed_bayes=rbinom(n=mis_no, size=1, prob=exp(beta0_bayes+beta1_bayes*x_miss+beta2_bayes*u_miss+beta3_bayes*z_miss)/(1+exp(beta0_bayes+beta1_bayes*x_miss+beta2_bayes*u_miss+beta3_bayes*z_miss)));


y2_completed_draw[miss_seq]=y2_imputed_draw;
y2_completed_mle[miss_seq]=y2_imputed_mle;
y2_completed_boot[miss_seq]=y2_imputed_boot;
y2_completed_lda[miss_seq]=as.numeric(y2_imputed_lda);
y2_completed_bayes[miss_seq]=y2_imputed_bayes;

# normal draw imputation;

draw_mean=mean(y2_completed_draw);
draw_var=var(y2_completed_draw)/rowobs;

draw_mean_mat[cycle,i]=draw_mean;
draw_var_mat[cycle,i]=draw_var;

draw_logistic=summary(lm(z~0+y2_completed_draw+x+u));
draw_slope_mat[cycle,i]=draw_logistic$coef[1,1];
draw_slope_var_mat[cycle,i]=draw_logistic$coef[1,2]^2;




# pseudo observation imputation;
mle_mean=mean(y2_completed_mle);
mle_var=var(y2_completed_mle)/rowobs;

mle_mean_mat[cycle,i]=mle_mean;
mle_var_mat[cycle,i]=mle_var;

mle_logistic=summary(lm(z~0+y2_completed_mle+x+u));
mle_slope_mat[cycle,i]=mle_logistic$coef[1,1];
mle_slope_var_mat[cycle,i]=mle_logistic$coef[1,2]^2;

# bootstrap imputation;

boot_mean=mean(y2_completed_boot);
boot_var=var(y2_completed_boot)/rowobs;

boot_mean_mat[cycle,i]=boot_mean;
boot_var_mat[cycle,i]=boot_var;

boot_logistic=summary(lm(z~0+y2_completed_boot+x+u));
boot_slope_mat[cycle,i]=boot_logistic$coef[1,1];
boot_slope_var_mat[cycle,i]=boot_logistic$coef[1,2]^2;



# discriminant analysis imputation;

lda_mean=mean(y2_completed_lda);
lda_mean_var=var(y2_completed_lda)/rowobs;

lda_mean_mat[cycle,i]=lda_mean;
lda_var_mat[cycle,i]=lda_mean_var;

lda_logistic=summary(lm(z~0+y2_completed_lda+x+u));
lda_slope_mat[cycle,i]=lda_logistic$coef[1,1];
lda_slope_var_mat[cycle,i]=lda_logistic$coef[1,2]^2;


# Bayesian imputation;

bayes_mean=mean(y2_completed_bayes);
bayes_mean_var=var(y2_completed_bayes)/rowobs;

bayes_mean_mat[cycle,i]=bayes_mean;
bayes_var_mat[cycle,i]=bayes_mean_var;

bayes_logistic=summary(lm(z~0+y2_completed_bayes+x+u));
bayes_slope_mat[cycle,i]=bayes_logistic$coef[1,1];
bayes_slope_var_mat[cycle,i]=bayes_logistic$coef[1,2]^2;



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

draw_mi_mean=rep(NA, cycle_no);
draw_mi_var=rep(NA, cycle_no);
draw_mi_df=rep(NA, cycle_no);
draw_mi_f=rep(NA, cycle_no);


mle_mi_mean=rep(NA, cycle_no);
mle_mi_var=rep(NA, cycle_no);
mle_mi_df=rep(NA, cycle_no);
mle_mi_f=rep(NA, cycle_no);


boot_mi_mean=rep(NA, cycle_no);
boot_mi_var=rep(NA, cycle_no);
boot_mi_df=rep(NA, cycle_no);
boot_mi_f=rep(NA, cycle_no);

lda_mi_mean=rep(NA, cycle_no);
lda_mi_var=rep(NA, cycle_no);
lda_mi_df=rep(NA, cycle_no);
lda_mi_f=rep(NA, cycle_no);

bayes_mi_mean=rep(NA, cycle_no);
bayes_mi_var=rep(NA, cycle_no);
bayes_mi_df=rep(NA, cycle_no);
bayes_mi_f=rep(NA, cycle_no);



for (i in 1:cycle_no)
{
draw_summary=pool.scalar(draw_mean_mat[i,], draw_var_mat[i,], n=rowobs);
draw_mi_mean[i]=draw_summary$qbar;
draw_mi_var[i]=draw_summary$t;
draw_mi_df[i]=draw_summary$df;
draw_mi_f[i]=draw_summary$f;

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

lda_summary=pool.scalar(lda_mean_mat[i,], lda_var_mat[i,], n=rowobs);
lda_mi_mean[i]=lda_summary$qbar;
lda_mi_var[i]=lda_summary$t;
lda_mi_df[i]=lda_summary$df;
lda_mi_f[i]=lda_summary$f;

bayes_summary=pool.scalar(bayes_mean_mat[i,], bayes_var_mat[i,], n=rowobs);
bayes_mi_mean[i]=bayes_summary$qbar;
bayes_mi_var[i]=bayes_summary$t;
bayes_mi_df[i]=bayes_summary$df;
bayes_mi_f[i]=bayes_summary$f;




}

# original normal imputation;
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


# pseudo observation imputation;
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


# lda imputation;
lda_mi_true=mean(lda_mi_mean);
lda_mi_bias=lda_mi_true-true;

lda_mi_bias;

lda_mi_bias/true;

lda_mi_mse=mean((lda_mi_mean-true)^2);
lda_mi_mse;

# coverage;
lda_mean_low95_vec=lda_mi_mean-qt(.975, lda_mi_df)*sqrt(lda_mi_var);
lda_mean_up95_vec=lda_mi_mean+qt(.975, lda_mi_df)*sqrt(lda_mi_var);

mean_lda_length=mean(lda_mean_up95_vec-lda_mean_low95_vec);
mean_lda_length;

lda_mi_coverage=(lda_mean_low95_vec < true)*(lda_mean_up95_vec > true);
mean(lda_mi_coverage);

# variance estimates;
mean(sqrt(lda_mi_var));
sqrt(var(lda_mi_mean));



# Bayes imputation;

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

draw_mi_mean=rep(NA, cycle_no);
draw_mi_var=rep(NA, cycle_no);
draw_mi_df=rep(NA, cycle_no);
draw_mi_f=rep(NA, cycle_no);



mle_mi_mean=rep(NA, cycle_no);
mle_mi_var=rep(NA, cycle_no);
mle_mi_df=rep(NA, cycle_no);
mle_mi_f=rep(NA, cycle_no);

boot_mi_mean=rep(NA, cycle_no);
boot_mi_var=rep(NA, cycle_no);
boot_mi_df=rep(NA, cycle_no);
boot_mi_f=rep(NA, cycle_no);

lda_mi_mean=rep(NA, cycle_no);
lda_mi_var=rep(NA, cycle_no);
lda_mi_df=rep(NA, cycle_no);
lda_mi_f=rep(NA, cycle_no);



bayes_mi_mean=rep(NA, cycle_no);
bayes_mi_var=rep(NA, cycle_no);
bayes_mi_df=rep(NA, cycle_no);
bayes_mi_f=rep(NA, cycle_no);



for (i in 1:cycle_no)
{

draw_summary=pool.scalar(draw_slope_mat[i,], draw_slope_var_mat[i,], n=rowobs);
draw_mi_mean[i]=draw_summary$qbar;
draw_mi_var[i]=draw_summary$t;
draw_mi_df[i]=draw_summary$df;
draw_mi_f[i]=draw_summary$f;



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

lda_summary=pool.scalar(lda_slope_mat[i,], lda_slope_var_mat[i,], n=rowobs);
lda_mi_mean[i]=lda_summary$qbar;
lda_mi_var[i]=lda_summary$t;
lda_mi_df[i]=lda_summary$df;
lda_mi_f[i]=lda_summary$f;



bayes_summary=pool.scalar(bayes_slope_mat[i,], bayes_slope_var_mat[i,], n=rowobs);
bayes_mi_mean[i]=bayes_summary$qbar;
bayes_mi_var[i]=bayes_summary$t;
bayes_mi_df[i]=bayes_summary$df;
bayes_mi_f[i]=bayes_summary$f;



}



# For slope;
# normal approximation;

draw_mi_bias=mean(draw_mi_mean)-true_slope;
draw_mi_bias;

draw_mi_bias/true_slope;

draw_mi_mse=mean((draw_mi_mean-true_slope)^2);
draw_mi_mse;

# coverage;
draw_slope_low95_vec=draw_mi_mean-qt(.975, draw_mi_df)*sqrt(draw_mi_var);
draw_slope_up95_vec=draw_mi_mean+qt(.975, draw_mi_df)*sqrt(draw_mi_var);

draw_slope_length=mean(draw_slope_up95_vec-draw_slope_low95_vec);
draw_slope_length;


draw_slope_coverage=(draw_slope_low95_vec < true_slope)*(draw_slope_up95_vec > true_slope);

mean(draw_slope_coverage);

# variance estimates;
mean(sqrt(draw_mi_var));
sqrt(var(draw_mi_mean));




# pseudo observation imputation;

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

# bootstrap imputation;

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


# lda imputation;

lda_mi_bias=mean(lda_mi_mean)-true_slope;
lda_mi_bias;

lda_mi_bias/true_slope;

lda_mi_mse=mean((lda_mi_mean-true_slope)^2);
lda_mi_mse;

# coverage;
lda_slope_low95_vec=lda_mi_mean-qt(.975, lda_mi_df)*sqrt(lda_mi_var);
lda_slope_up95_vec=lda_mi_mean+qt(.975, lda_mi_df)*sqrt(lda_mi_var);

lda_slope_length=mean(lda_slope_up95_vec-lda_slope_low95_vec);
lda_slope_length;

lda_slope_coverage=(lda_slope_low95_vec < true_slope)*(lda_slope_up95_vec > true_slope);

mean(lda_slope_coverage);

# variance estimates;
mean(sqrt(lda_mi_var));
sqrt(var(lda_mi_mean));




# Bayesian imputation;

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


hist(coef_s1x_draw_collection[,2], main="Hisogram of posterior draws", xlab="gamma1");



# test=MCMClogit(y2_obs~x_obs, burnin=1000, mcmc=2000, tune=0.40, beta.start=c(beta0, beta1));
# x_posterior=bayesglm(hospice_ac[true==1,2]~hospice_ac[true==1,3:n_col], family=binomial(link="logit"), prior.scale=2.5, prior.df=1);  
    # Same as M3, explicitly specifying Cauchy prior with scale 2.5 coef_s1