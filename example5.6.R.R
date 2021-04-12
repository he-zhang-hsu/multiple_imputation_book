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







rowobs=10000;

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
beta1=1; # beta1=1 creates around 30% of 1's;
# beta2=1/2;
beta2=1/2;
# beta1=-1; # beta1=-1 creates around 7% of 1's;
# beta1=-2; # beta1=-2 creates around 6% of 1's;
# beta1=-3; # beta1=-3 creates around 7% of 1's;
# beta1=-4; # beta1=-4 creates around 1.57% of 1's;
# beta1=-5; # beta1=-5 creates 1.25% of 1's;
# beta1=-6;



cycle_no=1000;
mi_no=30;
draw_no=1000;

# matrices holding the parameter estimates;
# complete-data inferences;
y2_mean_vec=y2_var_vec=cc_mean_vec=cc_var_vec=y2_prop_vec=y2_prop_var_vec=y2_cc_prop_vec=y2_cc_prop_var_vec=rep(NA, cycle_no);
y2_cc_prop_raghu_vec=y2_cc_prop_raghu_var_vec=y2_prop_raghu_vec=y2_prop_raghu_var_vec=rep(NA, cycle_no);

obs_prop_mat=obs_prop_var_mat=matrix(NA, cycle_no, mi_no);
raghu_prop_mat=raghu_prop_var_mat=matrix(NA, cycle_no, mi_no);

marginal_mean_mat=marginal_var_mat=obs_mean_mat=obs_var_mat=trans_mean_mat=trans_var_mat=guess_mean_mat=guess_var_mat=draw_mean_mat=draw_var_mat=lambdatrue_mean_mat=lambdatrue_var_mat=boot_mean_mat=boot_var_mat=bayes_mean_mat=bayes_var_mat=pmm_mean_mat=pmm_var_mat=matrix(NA, cycle_no, mi_no);

BD_slope_vec=BD_slope_var_vec=rep(NA, cycle_no);
CC_slope_vec=CC_slope_var_vec=rep(NA, cycle_no);
mle_slope_mat=mle_slope_var_mat=boot_slope_mat=boot_slope_var_mat=bayes_slope_mat=bayes_slope_var_mat=pmm_slope_mat=pmm_slope_var_mat=matrix(NA, cycle_no, mi_no);


for (cycle in cycle_no:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);
# x=runif(rowobs, min=0, max=2);

# generate the logistic regression outcome y2;
# Model I
# y2=rbinom(n=rowobs, size=1, prob=exp(beta0+beta1*x)/(1+exp(beta0+beta1*x)));

# Model II
y2=rbinom(n=rowobs, size=1, prob=exp(beta0+beta1*x+beta2*x^2)/(1+exp(beta0+beta1*x+beta2*x^2)));




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

# observed data is cbind(x, y2_miss);

# compelete data analysis;
# sample mean;
y2_mean=mean(y2);
y2_var=var(y2)/rowobs;
y2_mean_vec[cycle]=y2_mean;
y2_var_vec[cycle]=y2_var;

# logistic regression coefficient for slope;
# under model 1
# BD_logistic=summary(glm(y2~x, family=binomial(link="logit")));
# under model 2
BD_logistic=summary(glm(y2~x+I(x*x), family=binomial(link="logit")));

# BD_slope_vec[cycle]=BD_logistic$coeff[2,1];
# BD_slope_var_vec[cycle]=BD_logistic$coeff[2,2]^2;

BD_slope_vec[cycle]=BD_logistic$coeff[3,1];
BD_slope_var_vec[cycle]=BD_logistic$coeff[3,2]^2;



# different missing data methods;

# complete-case analysis;
# mean estimands;
cc_mean=mean(y2_miss, na.rm="T");
cc_mean_vec[cycle]=cc_mean;
cc_var_vec[cycle]=var(y2_miss, na.rm="T")/obs_no;

# logistic regression coefficient for slope;
# logistic regression coefficient for slope;
# under model 1
# CC_logistic=summary(glm(y2_obs~x_obs, family=binomial(link="logit")));
# under model 2
CC_logistic=summary(glm(y2_obs~x_obs+I(x_obs*x_obs), family=binomial(link="logit")));

# CC_slope_vec[cycle]=CC_logistic$coeff[2,1];
# CC_slope_var_vec[cycle]=CC_logistic$coeff[2,2]^2;

CC_slope_vec[cycle]=CC_logistic$coeff[3,1];
CC_slope_var_vec[cycle]=CC_logistic$coeff[3,2]^2;




# now impute the missing y2's to get estimates of the marginal;
y2_completed_mle=y2_completed_boot=y2_completed_bayes=y2_completed_pmm=y2_miss;

for (i in 1:mi_no)
{
set.seed(i);

# bivariate model imputation; 
# Rubin's imputation based on MLE;
y2_imputed_mle=mice.impute.logreg(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=as.matrix(x));

# linear normal imputation and then rounding;
y2_imputed_boot=mice.impute.norm(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=as.matrix(x));
y2_imputed_boot[y2_imputed_boot >=0.5]=1;
y2_imputed_boot[y2_imputed_boot < 0.5]=0;

# pmm imputation based on a linear normal model;
y2_imputed_bayes=mice.impute.pmm(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=as.matrix(x));


# pmm imputation based on a logit model;
y2_imputed_pmm=rep(NA, mis_no);

# under model 2;
CC_logistic=summary(glm(y2_obs~x_obs, family=binomial(link="logit")));
obs_pred_mean=cbind(cbind(1, x_obs)%*%CC_logistic$coefficients[,1],y2_obs);

# draw logistic regression coefficients from a posterior;
coef_s1x_posterior=bayesglm(y2_obs~x_obs, family=binomial(link="logit"), prior.scale=Inf, prior.df=Inf);  
coef_s1x_draw_collection=sim(coef_s1x_posterior, n.sims=draw_no)@coef;

coef_s1x_draw=coef_s1x_draw_collection[draw_no,];

miss_pred_draw=cbind(1, x_miss)%*%coef_s1x_draw;

# choose number of donors as 5;
donor_no=5;
for(k in 1:mis_no)
{
donor_matrix=obs_pred_mean;
donor_matrix[,1]=abs(miss_pred_draw[k]-obs_pred_mean[,1]);
donor_matrix_sorted=donor_matrix[order(donor_matrix[,1]),];
y2_imputed_pmm[k]=sample(donor_matrix_sorted[1:donor_no,2],1);
}

# univariate model imputation
# y2_imputed=mice.impute.norm(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(rep(1,rowobs),x));

# mle: logistic regression imputation;

y2_completed_mle[miss_seq]=y2_imputed_mle;

# boot: rounding;

y2_completed_boot[miss_seq]=y2_imputed_boot;

# bayes pmm imputation based on a linera normal model;

y2_completed_bayes[miss_seq]=y2_imputed_bayes;

# pmm: pmm imputation based on the glm model;

y2_completed_pmm[miss_seq]=y2_imputed_pmm;

# logistic regression imputation;
obs_mean=mean(y2_completed_mle);
obs_var=var(y2_completed_mle)/rowobs;

obs_mean_mat[cycle,i]=obs_mean;
obs_var_mat[cycle,i]=obs_var;

# logistic regression coefficient for slope;
# under model 1
# mle_logistic=summary(glm(y2_completed_mle~x, family=binomial(link="logit")));

# under model 2
mle_logistic=summary(glm(y2_completed_mle~x+I(x*x), family=binomial(link="logit")));

# mle_slope_mat[cycle,i]=mle_logistic$coeff[2,1];
# mle_slope_var_mat[cycle,i]=mle_logistic$coeff[2,2]^2;

mle_slope_mat[cycle,i]=mle_logistic$coeff[3,1];
mle_slope_var_mat[cycle,i]=mle_logistic$coeff[3,2]^2;



# rounding imputation;

trans_mean=mean(y2_completed_boot);
trans_mean_var=var(y2_completed_boot)/rowobs;

trans_mean_mat[cycle,i]=trans_mean;
trans_var_mat[cycle,i]=trans_mean_var;



# logistic regression coefficient for slope;
# under model 1
# boot_logistic=summary(glm(y2_completed_boot~x, family=binomial(link="logit")));
# under model 2
boot_logistic=summary(glm(y2_completed_boot~x+I(x*x), family=binomial(link="logit")));

# boot_slope_mat[cycle,i]=boot_logistic$coeff[2,1];
# boot_slope_var_mat[cycle,i]=boot_logistic$coeff[2,2]^2;

boot_slope_mat[cycle,i]=boot_logistic$coeff[3,1];
boot_slope_var_mat[cycle,i]=boot_logistic$coeff[3,2]^2;


# pmm on linear model imputation;

marginal_mean=mean(y2_completed_bayes);
marginal_mean_var=var(y2_completed_bayes)/rowobs;

marginal_mean_mat[cycle,i]=marginal_mean;
marginal_var_mat[cycle,i]=marginal_mean_var;

# logistic regression coefficient for slope;
# under model 1
# bayes_logistic=summary(glm(y2_completed_bayes~x, family=binomial(link="logit")));
# under model 2
bayes_logistic=summary(glm(y2_completed_bayes~x+I(x*x), family=binomial(link="logit")));

# bayes_slope_mat[cycle,i]=bayes_logistic$coeff[2,1];
# bayes_slope_var_mat[cycle,i]=bayes_logistic$coeff[2,2]^2;

bayes_slope_mat[cycle,i]=bayes_logistic$coeff[3,1];
bayes_slope_var_mat[cycle,i]=bayes_logistic$coeff[3,2]^2;



# pmm on logistic regression model imputation;

pmm_mean=mean(y2_completed_pmm);
pmm_mean_var=var(y2_completed_pmm)/rowobs;

pmm_mean_mat[cycle,i]=pmm_mean;
pmm_var_mat[cycle,i]=pmm_mean_var;

# logistic regression coefficient for slope;
# under model 1
# pmm_logistic=summary(glm(y2_completed_pmm~x, family=binomial(link="logit")));
# under model 2
pmm_logistic=summary(glm(y2_completed_pmm~x+I(x*x), family=binomial(link="logit")));

# pmm_slope_mat[cycle,i]=pmm_logistic$coeff[2,1];
# pmm_slope_var_mat[cycle,i]=pmm_logistic$coeff[2,2]^2;

pmm_slope_mat[cycle,i]=pmm_logistic$coeff[3,1];
pmm_slope_var_mat[cycle,i]=pmm_logistic$coeff[3,2]^2;



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

# logistic imputation;

obs_mi_mean=rep(NA, cycle_no);
obs_mi_var=rep(NA, cycle_no);
obs_mi_df=rep(NA, cycle_no);
obs_mi_f=rep(NA, cycle_no);

# rounding;

trans_mi_mean=rep(NA, cycle_no);
trans_mi_var=rep(NA, cycle_no);
trans_mi_df=rep(NA, cycle_no);
trans_mi_f=rep(NA, cycle_no);

# pmm on linear model; 

marginal_mi_mean=rep(NA, cycle_no);
marginal_mi_var=rep(NA, cycle_no);
marginal_mi_df=rep(NA, cycle_no);
marginal_mi_f=rep(NA, cycle_no);

# pmm on logistic model;

pmm_mi_mean=rep(NA, cycle_no);
pmm_mi_var=rep(NA, cycle_no);
pmm_mi_df=rep(NA, cycle_no);
pmm_mi_f=rep(NA, cycle_no);



for (i in 1:cycle_no)
{
obs_summary=pool.scalar(obs_mean_mat[i,], obs_var_mat[i,], n=rowobs);
obs_mi_mean[i]=obs_summary$qbar;
obs_mi_var[i]=obs_summary$t;
obs_mi_df[i]=obs_summary$df;
obs_mi_f[i]=obs_summary$f;

trans_summary=pool.scalar(trans_mean_mat[i,], trans_var_mat[i,], n=rowobs);
trans_mi_mean[i]=trans_summary$qbar;
trans_mi_var[i]=trans_summary$t;
trans_mi_df[i]=trans_summary$df;
trans_mi_f[i]=trans_summary$f;

marginal_summary=pool.scalar(marginal_mean_mat[i,], marginal_var_mat[i,], n=rowobs);
marginal_mi_mean[i]=marginal_summary$qbar;
marginal_mi_var[i]=marginal_summary$t;
marginal_mi_df[i]=marginal_summary$df;
marginal_mi_f[i]=marginal_summary$f;

pmm_summary=pool.scalar(pmm_mean_mat[i,], pmm_var_mat[i,], n=rowobs);
pmm_mi_mean[i]=pmm_summary$qbar;
pmm_mi_var[i]=pmm_summary$t;
pmm_mi_df[i]=pmm_summary$df;
pmm_mi_f[i]=pmm_summary$f;




}

# logistic imputation;
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


# rounding imputation;

trans_mi_true=mean(trans_mi_mean);
trans_mi_bias=trans_mi_true-true;

trans_mi_bias;

trans_mi_bias/true;

trans_mi_mse=mean((trans_mi_mean-true)^2);
trans_mi_mse;

# coverage;
trans_mean_low95_vec=trans_mi_mean-qt(.975, trans_mi_df)*sqrt(trans_mi_var);
trans_mean_up95_vec=trans_mi_mean+qt(.975, trans_mi_df)*sqrt(trans_mi_var);


mean_trans_length=mean(trans_mean_up95_vec-trans_mean_low95_vec);
mean_trans_length;


trans_mi_coverage=(trans_mean_low95_vec < true)*(trans_mean_up95_vec > true);
mean(trans_mi_coverage);

# variance estimates;
mean(sqrt(trans_mi_var));
sqrt(var(trans_mi_mean));

# pmm on linear model;

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

# pmm on logistic model;

pmm_mi_true=mean(pmm_mi_mean);
pmm_mi_bias=pmm_mi_true-true;

pmm_mi_bias;

pmm_mi_bias/true;

pmm_mi_mse=mean((pmm_mi_mean-true)^2);
pmm_mi_mse;

# coverage;
pmm_mean_low95_vec=pmm_mi_mean-qt(.975, pmm_mi_df)*sqrt(pmm_mi_var);
pmm_mean_up95_vec=pmm_mi_mean+qt(.975, pmm_mi_df)*sqrt(pmm_mi_var);


mean_pmm_length=mean(pmm_mean_up95_vec-pmm_mean_low95_vec);
mean_pmm_length;


pmm_mi_coverage=(pmm_mean_low95_vec < true)*(pmm_mean_up95_vec > true);
mean(pmm_mi_coverage);

# variance estimates;
mean(sqrt(pmm_mi_var));
sqrt(var(pmm_mi_mean));

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
# logistic imputation;

mle_mi_mean=rep(NA, cycle_no);
mle_mi_var=rep(NA, cycle_no);
mle_mi_df=rep(NA, cycle_no);
mle_mi_f=rep(NA, cycle_no);

# rounding;

boot_mi_mean=rep(NA, cycle_no);
boot_mi_var=rep(NA, cycle_no);
boot_mi_df=rep(NA, cycle_no);
boot_mi_f=rep(NA, cycle_no);

# pmm on linear model;

bayes_mi_mean=rep(NA, cycle_no);
bayes_mi_var=rep(NA, cycle_no);
bayes_mi_df=rep(NA, cycle_no);
bayes_mi_f=rep(NA, cycle_no);

# pmm on logistic model;

pmm_mi_mean=rep(NA, cycle_no);
pmm_mi_var=rep(NA, cycle_no);
pmm_mi_df=rep(NA, cycle_no);
pmm_mi_f=rep(NA, cycle_no);



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

pmm_summary=pool.scalar(pmm_slope_mat[i,], pmm_slope_var_mat[i,], n=rowobs);
pmm_mi_mean[i]=pmm_summary$qbar;
pmm_mi_var[i]=pmm_summary$t;
pmm_mi_df[i]=pmm_summary$df;
pmm_mi_f[i]=pmm_summary$f;



}



# For slope;

# proper imputation;
# logistic imputation;

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

# rounding imputation;

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

# pmm on linear model;

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

# pmm on logistic model;

pmm_mi_bias=mean(pmm_mi_mean)-true_slope;
pmm_mi_bias;

pmm_mi_bias/true_slope;

pmm_mi_mse=mean((pmm_mi_mean-true_slope)^2);
pmm_mi_mse;

# coverage;
pmm_slope_low95_vec=pmm_mi_mean-qt(.975, pmm_mi_df)*sqrt(pmm_mi_var);
pmm_slope_up95_vec=pmm_mi_mean+qt(.975, pmm_mi_df)*sqrt(pmm_mi_var);

pmm_slope_length=mean(pmm_slope_up95_vec-pmm_slope_low95_vec);
pmm_slope_length;


pmm_slope_coverage=(pmm_slope_low95_vec < true_slope)*(pmm_slope_up95_vec > true_slope);

mean(pmm_slope_coverage);

# variance estimates;
mean(sqrt(pmm_mi_var));
sqrt(var(pmm_mi_mean));

# observed data;
# logistic imputation;
plot(ecdf(x_obs[y2_obs==1]), col="red");
lines(ecdf(x_obs[y2_obs==0]), col="green");

lines(ecdf(x_miss[y2_imputed_mle==1]), col="blue");
lines(ecdf(x_miss[y2_imputed_mle==0]), col="yellow");

# rounding imputation;

plot(ecdf(x_obs[y2_obs==1]), col="red");
lines(ecdf(x_obs[y2_obs==0]), col="green");

lines(ecdf(x_miss[y2_imputed_boot==1]), col="blue");
lines(ecdf(x_miss[y2_imputed_boot==0]), col="yellow");

# PMM based on the linear normal model;

plot(ecdf(x_obs[y2_obs==1]), col="red");
lines(ecdf(x_obs[y2_obs==0]), col="green");

lines(ecdf(x_miss[y2_imputed_bayes==1]), col="blue");
lines(ecdf(x_miss[y2_imputed_bayes==0]), col="yellow");


# PMM based on the logistic model;

plot(ecdf(x_obs[y2_obs==1]), col="red");
lines(ecdf(x_obs[y2_obs==0]), col="green");

lines(ecdf(x_miss[y2_imputed_pmm==1]), col="blue");
lines(ecdf(x_miss[y2_imputed_pmm==0]), col="yellow");



# test=MCMClogit(y2_obs~x_obs, burnin=1000, mcmc=2000, tune=0.40, beta.start=c(beta0, beta1));
# x_posterior=bayesglm(hospice_ac[true==1,2]~hospice_ac[true==1,3:n_col], family=binomial(link="logit"), prior.scale=2.5, prior.df=1);  
    # Same as M3, explicitly specifying Cauchy prior with scale 2.5 coef_s1