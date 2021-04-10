################## data generate ##############
# Program: boxcox_impute_expolore_simulation; ############
# generate data, and try out different missing data methods;
# run some simulations to assess the biasedness of various methods;


# run some simulations to compare logistic regression imputation and discriminant analysis models;

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

# error=rnorm(rowobs);

# summary(glm(y2~x, family=binomial(link="logit")));
# mean(x[y2==0]);
# mean(x[y2==1]);


# set up the random seed;
set.seed(197789);

# generate data;
# distribution of x, the covariate;
# could be normal, or could be other distributions;
# fixed at mean 10, and variance 1;
# mu_x=1;
# var_x=1;

# regression parameters;
# both beta0 and beta1 are fixed at 
# beta0=-2;
# beta1=1; # beta1=1 creates around 30% of 1's;
# beta1=-1; # beta1=-1 creates around 7% of 1's;
# beta1=-2; # beta1=-2 creates around 6% of 1's;
# beta1=-3; # beta1=-3 creates around 7% of 1's;
# beta1=-4; # beta1=-4 creates around 1.57% of 1's;
# beta1=-5; # beta1=-5 creates 1.25% of 1's;
# beta1=-6;



cycle_no=1000;
mi_no=50;

p=0.5;
# p=0.1;



# c=1;
# c=2;
# c=3;
c=4;

# matrices holding the parameter estimates;
# complete-data inferences;
y2_mean_vec=y2_var_vec=cc_mean_vec=cc_var_vec=rep(NA, cycle_no);

mle_mean_mat=mle_var_mat=lda_mean_mat=lda_var_mat=matrix(NA, cycle_no, mi_no);

BD_slope_vec=BD_slope_var_vec=rep(NA, cycle_no);
CC_slope_vec=CC_slope_var_vec=rep(NA, cycle_no);
mle_slope_mat=mle_slope_var_mat=lda_slope_mat=lda_slope_var_mat=matrix(NA, cycle_no, mi_no);


for (cycle in 1:cycle_no)

{
set.seed(cycle);
# generate the model using the discriminant analysis model;

y2=rbinom(n=rowobs, size=1, prob=p);

error=rnorm(rowobs);

# generate the discriminant x;
# the mean of x is 1-0=1;
x=1+c*y2+error;


# generate the model using the logistic model;

# set up the missing data as MCAR on x;
# corresponding to the 1st discriminant analysis;
# mu_x=1.5;
# var_x=1.25;
# x=mu_x+sqrt(var_x)*rnorm(rowobs);
# regression parameters;
# both beta0 and beta1 are fixed at 
# beta0=-3/2;
# beta1=1;

# corresponding to the 2nd discriminant analysis;
# mu_x=2;
# var_x=2;
# x=mu_x+sqrt(var_x)*rnorm(rowobs);
# regression parameters;
# both beta0 and beta1 are fixed at 
# beta0=-4;
# beta1=2;

# corresponding to the 3rd discriminant analysis;
# mu_x=2.5;
# var_x=13/4;
# x=mu_x+sqrt(var_x)*rnorm(rowobs);
# regression parameters;
# both beta0 and beta1 are fixed at 
# beta0=-7.5;
# beta1=3;

# corresponding to the 4th discriminant analysis;
# mu_x=3;
# var_x=5;
# x=mu_x+sqrt(var_x)*rnorm(rowobs);
# regression parameters;
# both beta0 and beta1 are fixed at 
# beta0=-13;
# beta1=4;



# generate the logistic regression outcome y2;
# y2=rbinom(n=rowobs, size=1, prob=exp(beta0+beta1*x)/(1+exp(beta0+beta1*x)));



# 40% missing;
# MCAR;
miss_indi=runif(n=rowobs)<0.40;

# MAR;
# set up the missing data as MAR on x;
#alpha0=-1.4;
# alpha0=-3;
# alpha1=1;
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x)/(1+exp(alpha0+alpha1*x)));


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

# x_y2_0=x[y2==0];
# n_y2_0=length(x_y2_0);
# x_y2_1=x[y2==1];
# n_y2_1=length(x_y2_1);

# mean of X given y2=0;
# x_y2_0_mean=mean(x_y2_0);
# x_y2_0_var=var(x_y2_0)/n_y2_0;
# x_y2_0_mean_vec[cycle]=x_y2_0_mean;
# x_y2_0_var_vec[cycle]=x_y2_0_var;

# mean of X given y2=1;
# x_y2_1_mean=mean(x_y2_1);
# x_y2_1_var=var(x_y2_1)/n_y2_1;
# x_y2_1_mean_vec[cycle]=x_y2_1_mean;
# x_y2_1_var_vec[cycle]=x_y2_1_var;




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

# x_obs_y2_0=x_obs[y2_obs==0];
# n_obs_y2_0=length(x_obs_y2_0);
# x_obs_y2_1=x_obs[y2_obs==1];
# n_obs_y2_1=length(x_obs_y2_1);

# mean of X given y2=0;
# x_obs_y2_0_mean=mean(x_obs_y2_0);
# x_obs_y2_0_var=var(x_obs_y2_0)/n_obs_y2_0;
# x_obs_y2_0_mean_vec[cycle]=x_obs_y2_0_mean;
# x_obs_y2_0_var_vec[cycle]=x_obs_y2_0_var;

# mean of X given y2=1;
# x_obs_y2_1_mean=mean(x_obs_y2_1);
# x_obs_y2_1_var=var(x_obs_y2_1)/n_obs_y2_1;
# x_obs_y2_1_mean_vec[cycle]=x_obs_y2_1_mean;
# x_obs_y2_1_var_vec[cycle]=x_obs_y2_1_var;




# now impute the missing y2's to get estimates of the marginal;
y2_completed_mle=y2_completed_lda=y2_miss;

for (i in 1:mi_no)
{
set.seed(i);
# bivariate model imputation; 
# Rubin's imputation based on MLE;
y2_imputed_mle=mice.impute.logreg(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=as.matrix(x));

# discriminant imputation based;
# y2_imputed_boot=mice.impute.lda(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=as.matrix(x));

y2_imputed_lda=mice.impute.lda.proper(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=as.matrix(x));


# test_data=as.data.frame(cbind(as.factor(y2_miss), x));
# test_data$V1=as.factor(test_data$V1);

# u=mice(data=test_data, m=1, method=c('lda'), seed=i);

# the imputed values are from a list;

# y2_imputed_boot=as.numeric(u$imp$V1[[1]])-1;


y2_completed_mle[miss_seq]=y2_imputed_mle;
y2_completed_lda[miss_seq]=as.numeric(y2_imputed_lda);

# marginal means;
mle_mean=mean(y2_completed_mle);
mle_var=var(y2_completed_mle)/rowobs;

mle_mean_mat[cycle,i]=mle_mean;
mle_var_mat[cycle,i]=mle_var;

# logistic regression coefficient for slope;
mle_logistic=summary(glm(y2_completed_mle~x, family=binomial(link="logit")));
mle_slope_mat[cycle,i]=mle_logistic$coeff[2,1];
mle_slope_var_mat[cycle,i]=mle_logistic$coeff[2,2]^2;

# means of X given Y;

# x_y2_mle_0=x[y2_completed_mle==0];
# n_y2_mle_0=length(x_y2_mle_0);
# x_y2_mle_1=x[y2_completed_mle==1];
# n_y2_mle_1=length(x_y2_mle_1);

# mean of X given y2=0;
# x_y2_mle_0_mean=mean(x_y2_mle_0);
# x_y2_mle_0_var=var(x_y2_mle_0)/n_y2_mle_0;
# x_y2_mle_0_mean_mat[cycle,i]=x_y2_mle_0_mean;
# x_y2_mle_0_var_mat[cycle,i]=x_y2_mle_0_var;

# mean of X given y2=1;
# x_y2_mle_1_mean=mean(x_y2_mle_1);
# x_y2_mle_1_var=var(x_y2_mle_1)/n_y2_mle_1;
# x_y2_mle_1_mean_mat[cycle,i]=x_y2_mle_1_mean;
# x_y2_mle_1_var_mat[cycle,i]=x_y2_mle_1_var;



lda_mean=mean(y2_completed_lda);
lda_mean_var=var(y2_completed_lda)/rowobs;

lda_mean_mat[cycle,i]=lda_mean;
lda_var_mat[cycle,i]=lda_mean_var;



# logistic regression coefficient for slope;
lda_logistic=summary(glm(y2_completed_lda~x, family=binomial(link="logit")));
lda_slope_mat[cycle,i]=lda_logistic$coeff[2,1];
lda_slope_var_mat[cycle,i]=lda_logistic$coeff[2,2]^2;

# means of X given Y;

# x_y2_boot_0=x[y2_completed_boot==0];
# n_y2_boot_0=length(x_y2_boot_0);
# x_y2_boot_1=x[y2_completed_boot==1];
# n_y2_boot_1=length(x_y2_boot_1);

# mean of X given y2=0;
# x_y2_boot_0_mean=mean(x_y2_boot_0);
# x_y2_boot_0_var=var(x_y2_boot_0)/n_y2_boot_0;
# x_y2_boot_0_mean_mat[cycle,i]=x_y2_boot_0_mean;
# x_y2_boot_0_var_mat[cycle,i]=x_y2_boot_0_var;

# mean of X given y2=1;
# x_y2_boot_1_mean=mean(x_y2_boot_1);
# x_y2_boot_1_var=var(x_y2_boot_1)/n_y2_boot_1;
# x_y2_boot_1_mean_mat[cycle,i]=x_y2_boot_1_mean;
# x_y2_boot_1_var_mat[cycle,i]=x_y2_boot_1_var;



}

cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

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


lda_mi_mean=rep(NA, cycle_no);
lda_mi_var=rep(NA, cycle_no);
lda_mi_df=rep(NA, cycle_no);
lda_mi_f=rep(NA, cycle_no);


for (i in 1:cycle_no)
{
mle_summary=pool.scalar(mle_mean_mat[i,], mle_var_mat[i,], n=rowobs);
mle_mi_mean[i]=mle_summary$qbar;
mle_mi_var[i]=mle_summary$t;
mle_mi_df[i]=mle_summary$df;
mle_mi_f[i]=mle_summary$f;

lda_summary=pool.scalar(lda_mean_mat[i,], lda_var_mat[i,], n=rowobs);
lda_mi_mean[i]=lda_summary$qbar;
lda_mi_var[i]=lda_summary$t;
lda_mi_df[i]=lda_summary$df;
lda_mi_f[i]=lda_summary$f;


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


# discriminant imputation;
# lda_mi_mean=lda_mi_mean[!is.na(trans_mi_f)];
# lda_mi_var=t_mi_var[!is.na(trans_mi_f)];
# lda_mi_df=trans_mi_df[!is.na(trans_mi_f)];
# trans_mi_f=trans_mi_f[!is.na(trans_mi_f)];

lda_mi_true=mean(lda_mi_mean);
lda_mi_bias=lda_mi_true-true;

lda_mi_bias;

lda_mi_bias/true;

lda_mi_mse=mean((lda_mi_mean-true)^2);
lda_mi_mse;

# coverage;
lda_mean_low95_vec=lda_mi_mean-qt(.975, lda_mi_df)*sqrt(lda_mi_var);
lda_mean_up95_vec=lda_mi_mean+qt(.975, lda_mi_df)*sqrt(lda_mi_var);

# trans_mean_low95_vec=trans_mean_low95_vec[!is.na(trans_mean_low95_vec)];
# trans_mean_up95_vec=trans_mean_up95_vec[!is.na(trans_mean_up95_vec)];


mean_lda_length=mean(lda_mean_up95_vec-lda_mean_low95_vec);
mean_lda_length;


lda_mi_coverage=(lda_mean_low95_vec < true)*(lda_mean_up95_vec > true);
mean(lda_mi_coverage);

# variance estimates;
mean(sqrt(lda_mi_var));
sqrt(var(lda_mi_mean));


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

lda_mi_mean=rep(NA, cycle_no);
lda_mi_var=rep(NA, cycle_no);
lda_mi_df=rep(NA, cycle_no);
lda_mi_f=rep(NA, cycle_no);



for (i in 1:cycle_no)
{
mle_summary=pool.scalar(mle_slope_mat[i,], mle_slope_var_mat[i,], n=rowobs);
mle_mi_mean[i]=mle_summary$qbar;
mle_mi_var[i]=mle_summary$t;
mle_mi_df[i]=mle_summary$df;
mle_mi_f[i]=mle_summary$f;

lda_summary=pool.scalar(lda_slope_mat[i,], lda_slope_var_mat[i,], n=rowobs);
lda_mi_mean[i]=lda_summary$qbar;
lda_mi_var[i]=lda_summary$t;
lda_mi_df[i]=lda_summary$df;
lda_mi_f[i]=lda_summary$f;



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


###############################################################################
# conditional mean of X given Y;
# mean of x given y2 =0;

# plot the distribution of the data;

# par(mfrow=c(3,2));
# hist(x_obs[y2_obs==1]);
# hist(x_obs[y2_obs==0]);
# hist(x[y2_imputed_mle==1]);
# hist(x[y2_imputed_mle==0]);
# hist(x[y2_imputed_boot==1]);
# hist(x[y2_imputed_boot==0]);


# summary(x_obs);
# summary(x);
#summary(x_obs[y2_obs==1]);

#summary(x_obs[y2_obs==0]);

# breaks=seq(-5,9,0.25);

# par(mfrow=c(2,2));
# hist(x[y2==1], xlim=c(-4,11), col=rgb(1,0,0,0.5))
# hist(x[y2==0], col=rgb(0,0,1,0.5), add=T)

# hist(x, breaks=seq(-4,11,0.5));


# hist(x_obs[y2_obs==1], xlim=c(-4,11), breaks=breaks, xlab="x", main="Observed Data", col=rgb(1,0,0,0.5)) # red;
# hist(x_obs[y2_obs==0], col=rgb(0,0,1,0.5), breaks=breaks, add=T) # blue;

# hist(x[miss_indi==1 & y2==1], xlim=c(-4,11), xlab="x", main="Unobseved Data", breaks=breaks, col=rgb(1,0,0,0.5)) # red;
# hist(x[miss_indi==1 & y2==0], col=rgb(0,0,1,0.5), breaks=breaks, add=T) # blue;

# hist(x[y2_imputed_mle==1], xlim=c(-4,11), xlab="x", main="Logistic Imputation", breaks=breaks, col=rgb(1,0,0,0.5))
# hist(x[y2_imputed_mle==0], breaks=breaks, col=rgb(0,0,1,0.5), add=T)

# hist(x[y2_imputed_boot==1], xlim=c(-4,11), xlab="x", main="Discriminant Imputation", breaks=breaks, col=rgb(1,0,0,0.5))
# hist(x[y2_imputed_boot==0], breaks=breaks, col=rgb(0,0,1,0.5), add=T)



# par(mfrow=c(1,1));
# hist(x);
# hist(x[y2==1], xlim=c(-4,11), col=rgb(1,0,0,0.5))
# hist(x[y2==0], col=rgb(0,0,1,0.5), add=T)



