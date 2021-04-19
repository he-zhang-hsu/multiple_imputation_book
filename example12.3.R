################## data generate ##############
# Program: boxcox_impute_expolore_simulation; ############
# generate data, and try out different missing data methods;
# run some simulations to assess the biasedness of various methods;

# library(car);
# library(mvtnorm);
library(mice);
library(norm);
# library(HI);
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


improper_imputation=function(data)
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
 # beta=mvrnorm(n=1, mu=beta_hat, Sigma=sigma_square*inv_v_beta);

 parameter=c(beta_hat,sigma_square);
 parameter;
}





rowobs=1000;

# set up the random seed;
set.seed(197789);

# generate data;
# distribution of x, the covariate;
# could be normal, or could be other distributions;
# fixed at mean 10, and variance 1;
mu_x=1;
var_x=1;

# regression parameters;
# both beta0, beta1 and beta2 are fixed at 
beta0=-2;
# beta1=0;
beta1=1;
beta2=1;

# error variance;
var_error=1;

cycle=1;

set.seed(cycle);

x1=mu_x+sqrt(var_x)*rnorm(rowobs);

x2=mu_x+sqrt(var_x)*rnorm(rowobs);


# error distribution
# could be normal, or could be other distributions;
error=sqrt(var_error)*rnorm(rowobs);

# generate the latent variable;
y=beta0+beta1*x1+beta2*x2+error;

# set up the missing data as MCAR on x;
# miss_indi=runif(n=rowobs)<0.40;

# set up the missing data as MAR on x1;
alpha0=-1.4;
alpha1=1;
miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x1)/(1+exp(alpha0+alpha1*x1)));


y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];

x1_obs=x1[miss_indi==0];
x1_mis=x1[miss_indi==1];

x2_obs=x2[miss_indi==0];
x2_mis=x2[miss_indi==1];



mi_no=1;

y_completed_x1=y_completed_x2=y_completed_x1_x2=y_miss;


for (i in 1:mi_no)
{
# imputation using only x1;
# y_imputed_x1=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x1);


# y_completed_x1[miss_seq]=y_imputed_x1;

# imputation using only x2;
# y_imputed_x2=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x2);


# y_completed_x2[miss_seq]=y_imputed_x2;

# imputation using both x1 and x2;
y_imputed_x1_x2=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x1,x2));

y_completed_x1_x2[miss_seq]=y_imputed_x1_x2;


}


# the end of multiple imputation;

# plot the distribution of the data;
# par(mfrow=c(1,1));
# hist(y_obs, freq=FALSE);
# hist(y_imputed, freq=FALSE);

plot(x1_obs, y_obs);
abline(lm(y_obs~x1_obs));
lm(y_obs ~ x1_obs);

plot(x2_obs, y_obs);
abline(lm(y_obs~x2_obs));
lm(y_obs ~ x2_obs);



plot(x1_mis, y_imputed_x1, pch=17, ylab="y_imp: MI-X1", col="red");
abline(lm(y_imputed_x1 ~ x1_mis), col="red");
lm(y_imputed_x1 ~ x1_mis);

plot(x2_mis, y_imputed_x1, pch=17, ylab="y_imp: MI-X1", col="red");
abline(lm(y_imputed_x1 ~ x2_mis), col="red");
lm(y_imputed_x1 ~ x2_mis);

plot(x1_mis, y_imputed_x2, pch=17, ylab="y_imp: MI-X2", col="red");
abline(lm(y_imputed_x2 ~ x1_mis), col="red");
lm(y_imputed_x2 ~ x1_mis);

plot(x2_mis, y_imputed_x2, pch=17, ylab="y_imp: MI-X2", col="red");
abline(lm(y_imputed_x2 ~ x2_mis), col="red");
lm(y_imputed_x2 ~ x2_mis);

plot(x1_mis, y_imputed_x1_x2, pch=17, ylab="y_imp: MI-X1X2", col="red");
abline(lm(y_imputed_x1_x2 ~ x1_mis), col="red");
lm(y_imputed_x1_x2 ~ x1_mis);

plot(x2_mis, y_imputed_x1_x2, pch=17, ylab="y_imp: MI-X1X2", col="red");
abline(lm(y_imputed_x1_x2 ~ x2_mis), col="red");
lm(y_imputed_x1_x2 ~ x2_mis);


mean(y_obs);
var(y_obs);
mean(y_imputed);
var(y_imputed);

t.test(y_obs, y_imputed);
ks.test(y_obs, y_imputed);


# check the performance of the estimates;
# population quantity

#############################################################################
# marginal mean estimand;
true_mean=mean(BD_mean_vec);
true_mean;

# performance
BD_mean_mse=mean((BD_mean_vec-true_mean)^2);
BD_mean_mse;

# lower and upper 95% CI;
BD_mean_low95_vec=BD_mean_vec-1.96*sqrt(BD_mean_var_vec);
BD_mean_up95_vec=BD_mean_vec+1.96*sqrt(BD_mean_var_vec);

BD_mean_length=mean(BD_mean_up95_vec-BD_mean_low95_vec);
BD_mean_length;

BD_mean_cov=(BD_mean_low95_vec < true_mean)*(BD_mean_up95_vec > true_mean);

mean(BD_mean_cov);

# variance estimates;
mean(BD_mean_var_vec);
var(BD_mean_vec);


# complete case analysis;

CC_mean_bias=mean(CC_mean_vec)-true_mean;
CC_mean_bias;

CC_mean_mse=mean((CC_mean_vec-true_mean)^2);
CC_mean_mse;

# lower and upper 95% CI;
CC_mean_low95_vec=CC_mean_vec-1.96*sqrt(CC_mean_var_vec);
CC_mean_up95_vec=CC_mean_vec+1.96*sqrt(CC_mean_var_vec);

CC_mean_length=mean(CC_mean_up95_vec-CC_mean_low95_vec);
CC_mean_length;

CC_mean_cov=(CC_mean_low95_vec < true_mean)*(CC_mean_up95_vec > true_mean);

mean(CC_mean_cov);

# variance estimates;
mean(CC_mean_var_vec);
var(CC_mean_vec);


# multiple imputation estimates;
mean_mi_mean=rep(NA, cycle_no);
mean_mi_var=rep(NA, cycle_no);
mean_mi_df=rep(NA, cycle_no);
mean_mi_f=rep(NA, cycle_no);

mean_IM_mi_mean=rep(NA, cycle_no);
mean_IM_mi_var=rep(NA, cycle_no);
mean_IM_mi_df=rep(NA, cycle_no);
mean_IM_mi_f=rep(NA, cycle_no);



for (i in 1:cycle_no)
{
mean_summary=pool.scalar(MI_mean_mat[i,], MI_mean_var_mat[i,]);
mean_mi_mean[i]=mean_summary$qbar;
mean_mi_var[i]=mean_summary$t;
mean_mi_df[i]=mean_summary$df;
mean_mi_f[i]=mean_summary$f;

IM_mean_summary=pool.scalar(IM_MI_mean_mat[i,], IM_MI_mean_var_mat[i,]);
mean_IM_mi_mean[i]=IM_mean_summary$qbar;
mean_IM_mi_var[i]=IM_mean_summary$t;
mean_IM_mi_df[i]=IM_mean_summary$df;
mean_IM_mi_f[i]=IM_mean_summary$f;


}


# For marginal means;
# proper imputation;

mean_mi_bias=mean(mean_mi_mean)-true_mean;
mean_mi_bias;

mean_mi_mse=mean((mean_mi_mean-true_mean)^2);
mean_mi_mse;

# coverage;
MI_mean_low95_vec=mean_mi_mean-qt(.975, mean_mi_df)*sqrt(mean_mi_var);
MI_mean_up95_vec=mean_mi_mean+qt(.975, mean_mi_df)*sqrt(mean_mi_var);

MI_mean_length=mean(MI_mean_up95_vec-MI_mean_low95_vec);

MI_mean_length;

MI_mean_coverage=(MI_mean_low95_vec < true_mean)*(MI_mean_up95_vec > true_mean);
mean(MI_mean_coverage);

# variance estimates;
mean(mean_mi_var);
var(mean_mi_mean);

# improper imputation;

mean_IM_mi_bias=mean(mean_IM_mi_mean)-true_mean;
mean_IM_mi_bias;

mean_IM_mi_mse=mean((mean_IM_mi_mean-true_mean)^2);
mean_IM_mi_mse;

# coverage;
IM_MI_mean_low95_vec=mean_IM_mi_mean-qt(.975, mean_IM_mi_df)*sqrt(mean_IM_mi_var);
IM_MI_mean_up95_vec=mean_IM_mi_mean+qt(.975, mean_IM_mi_df)*sqrt(mean_IM_mi_var);

IM_MI_mean_length=mean(IM_MI_mean_up95_vec-IM_MI_mean_low95_vec);

IM_MI_mean_length;

IM_MI_mean_coverage=(IM_MI_mean_low95_vec < true_mean)*(IM_MI_mean_up95_vec > true_mean);
mean(IM_MI_mean_coverage);

# variance estimates;
mean(mean_IM_mi_var);
var(mean_IM_mi_mean);


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
mean(BD_slope_var_vec);
var(BD_slope_vec);


# complete-case analysis;
CC_slope_bias=mean(CC_slope_vec)-true_slope;
CC_slope_bias;

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
mean(CC_slope_var_vec);
var(CC_slope_vec);


# multiple imputation estimates;

slope_mi_mean=rep(NA, cycle_no);
slope_mi_var=rep(NA, cycle_no);
slope_mi_df=rep(NA, cycle_no);
slope_mi_f=rep(NA, cycle_no);

slope_IM_mi_mean=rep(NA, cycle_no);
slope_IM_mi_var=rep(NA, cycle_no);
slope_IM_mi_df=rep(NA, cycle_no);
slope_IM_mi_f=rep(NA, cycle_no);



for (i in 1:cycle_no)
{
slope_summary=pool.scalar(MI_slope_mat[i,], MI_slope_var_mat[i,]);
slope_mi_mean[i]=slope_summary$qbar;
slope_mi_var[i]=slope_summary$t;
slope_mi_df[i]=slope_summary$df;
slope_mi_f[i]=slope_summary$f;

IM_slope_summary=pool.scalar(IM_MI_slope_mat[i,], IM_MI_slope_var_mat[i,]);
slope_IM_mi_mean[i]=IM_slope_summary$qbar;
slope_IM_mi_var[i]=IM_slope_summary$t;
slope_IM_mi_df[i]=IM_slope_summary$df;
slope_IM_mi_f[i]=IM_slope_summary$f;



}



# For slope;

# proper imputation;

slope_mi_bias=mean(slope_mi_mean)-true_slope;
slope_mi_bias;

slope_mi_mse=mean((slope_mi_mean-true_slope)^2);
slope_mi_mse;

# coverage;
MI_slope_low95_vec=slope_mi_mean-qt(.975, slope_mi_df)*sqrt(slope_mi_var);
MI_slope_up95_vec=slope_mi_mean+qt(.975, slope_mi_df)*sqrt(slope_mi_var);

MI_slope_length=mean(MI_slope_up95_vec-MI_slope_low95_vec);
MI_slope_length;


MI_slope_coverage=(MI_slope_low95_vec < true_slope)*(MI_slope_up95_vec > true_slope);

mean(MI_slope_coverage);

# variance estimates;
mean(slope_mi_var);
var(slope_mi_mean);

# improper imputation;

slope_IM_mi_bias=mean(slope_IM_mi_mean)-true_slope;
slope_IM_mi_bias;

slope_IM_mi_mse=mean((slope_IM_mi_mean-true_slope)^2);
slope_IM_mi_mse;

# coverage;
IM_MI_slope_low95_vec=slope_IM_mi_mean-qt(.975, slope_IM_mi_df)*sqrt(slope_IM_mi_var);
IM_MI_slope_up95_vec=slope_IM_mi_mean+qt(.975, slope_IM_mi_df)*sqrt(slope_IM_mi_var);

IM_MI_slope_length=mean(IM_MI_slope_up95_vec-IM_MI_slope_low95_vec);
IM_MI_slope_length;


IM_MI_slope_coverage=(IM_MI_slope_low95_vec < true_slope)*(IM_MI_slope_up95_vec > true_slope);

mean(IM_MI_slope_coverage);

# variance estimates;
mean(slope_IM_mi_var);
var(slope_IM_mi_mean);

#######################################################################
# reverse slope;

# complete-data analysis;

true_xyslope=mean(BD_xyslope_vec);
true_xyslope;


# performance;
BD_xyslope_mse=mean((BD_xyslope_vec-true_xyslope)^2);
BD_xyslope_mse;

# lower and upper 95% CI;
BD_xyslope_low95_vec=BD_xyslope_vec-1.96*sqrt(BD_xyslope_var_vec);
BD_xyslope_up95_vec=BD_xyslope_vec+1.96*sqrt(BD_xyslope_var_vec);

BD_xyslope_length=mean(BD_xyslope_up95_vec-BD_xyslope_low95_vec);
BD_xyslope_length;

BD_xyslope_cov=(BD_xyslope_low95_vec < true_xyslope)*(BD_xyslope_up95_vec > true_xyslope);
mean(BD_xyslope_cov);

# variance estimates;
mean(BD_xyslope_var_vec);
var(BD_xyslope_vec);


# complete-case analysis;
CC_xyslope_bias=mean(CC_xyslope_vec)-true_xyslope;
CC_xyslope_bias;

CC_xyslope_mse=mean((CC_xyslope_vec-true_xyslope)^2);
CC_xyslope_mse;

# lower and upper 95% CI;
CC_xyslope_low95_vec=CC_xyslope_vec-1.96*sqrt(CC_xyslope_var_vec);
CC_xyslope_up95_vec=CC_xyslope_vec+1.96*sqrt(CC_xyslope_var_vec);

CC_xyslope_length=mean(CC_xyslope_up95_vec-CC_xyslope_low95_vec);
CC_xyslope_length;

CC_xyslope_cov=(CC_xyslope_low95_vec < true_xyslope)*(CC_xyslope_up95_vec > true_xyslope);

mean(CC_xyslope_cov);

# variance estimates;
mean(CC_xyslope_var_vec);
var(CC_xyslope_vec);


# multiple imputation estimates;

xyslope_mi_mean=rep(NA, cycle_no);
xyslope_mi_var=rep(NA, cycle_no);
xyslope_mi_df=rep(NA, cycle_no);
xyslope_mi_f=rep(NA, cycle_no);

xyslope_IM_mi_mean=rep(NA, cycle_no);
xyslope_IM_mi_var=rep(NA, cycle_no);
xyslope_IM_mi_df=rep(NA, cycle_no);
xyslope_IM_mi_f=rep(NA, cycle_no);



for (i in 1:cycle_no)
{
xyslope_summary=pool.scalar(MI_xyslope_mat[i,], MI_xyslope_var_mat[i,]);
xyslope_mi_mean[i]=xyslope_summary$qbar;
xyslope_mi_var[i]=xyslope_summary$t;
xyslope_mi_df[i]=xyslope_summary$df;
xyslope_mi_f[i]=xyslope_summary$f;

IM_xyslope_summary=pool.scalar(IM_MI_xyslope_mat[i,], IM_MI_xyslope_var_mat[i,]);
xyslope_IM_mi_mean[i]=IM_xyslope_summary$qbar;
xyslope_IM_mi_var[i]=IM_xyslope_summary$t;
xyslope_IM_mi_df[i]=IM_xyslope_summary$df;
xyslope_IM_mi_f[i]=IM_xyslope_summary$f;



}



# proper imputation;

xyslope_mi_bias=mean(xyslope_mi_mean)-true_xyslope;
xyslope_mi_bias;

xyslope_mi_mse=mean((xyslope_mi_mean-true_xyslope)^2);
xyslope_mi_mse;

# coverage;
MI_xyslope_low95_vec=xyslope_mi_mean-qt(.975, xyslope_mi_df)*sqrt(xyslope_mi_var);
MI_xyslope_up95_vec=xyslope_mi_mean+qt(.975, xyslope_mi_df)*sqrt(xyslope_mi_var);

MI_xyslope_length=mean(MI_xyslope_up95_vec-MI_xyslope_low95_vec);
MI_xyslope_length;


MI_xyslope_coverage=(MI_xyslope_low95_vec < true_xyslope)*(MI_xyslope_up95_vec > true_xyslope);

mean(MI_xyslope_coverage);

# variance estimates;
mean(xyslope_mi_var);
var(xyslope_mi_mean);

# improper imputation;

xyslope_IM_mi_bias=mean(xyslope_IM_mi_mean)-true_xyslope;
xyslope_IM_mi_bias;

xyslope_IM_mi_mse=mean((xyslope_IM_mi_mean-true_xyslope)^2);
xyslope_IM_mi_mse;

# coverage;
IM_MI_xyslope_low95_vec=xyslope_IM_mi_mean-qt(.975, xyslope_IM_mi_df)*sqrt(xyslope_IM_mi_var);
IM_MI_xyslope_up95_vec=xyslope_IM_mi_mean+qt(.975, xyslope_IM_mi_df)*sqrt(xyslope_IM_mi_var);

IM_MI_xyslope_length=mean(IM_MI_xyslope_up95_vec-IM_MI_xyslope_low95_vec);
IM_MI_xyslope_length;


IM_MI_xyslope_coverage=(IM_MI_xyslope_low95_vec < true_xyslope)*(IM_MI_xyslope_up95_vec > true_xyslope);

mean(IM_MI_xyslope_coverage);

# variance estimates;
mean(xyslope_IM_mi_var);
var(xyslope_IM_mi_mean);


