################## data generate ##############
# Program: boxcox_impute_expolore_simulation; ############
# generate data, and try out different missing data methods;
# run some simulations to assess the biasedness of various methods;

# library(car);
# library(mvtnorm);
library(mice);
library(norm);
library(boot);
library(Hmisc);
# library(HI);
options(digits=4);
rm(list=ls());

date();


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

# function to obtain R-Squared from the data 
rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
} 

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
estimand=as.data.frame(cbind(bias,rbias,mse,mean_length,mean_cov,sd,se))
}






rowobs=100;

# set up the random seed;
set.seed(197789);

# generate data;
# distribution of x, the covariate;
# could be normal, or could be other distributions;
# fixed at mean 10, and variance 1;
mu_x=1;
var_x=1;


# error variance;
var_error=1;


cycle_no=1000;
mi_no=50;
boot_no=200;

# matrices holding the parameter estimates;
# complete-data inferences;
y2_mean_vec=y2_var_vec=cc_mean_vec=cc_var_vec=y2_prop_vec=y2_prop_var_vec=y2_cc_prop_vec=y2_cc_prop_var_vec=rep(NA, cycle_no);
y2_cc_prop_raghu_vec=y2_cc_prop_raghu_var_vec=y2_prop_raghu_vec=y2_prop_raghu_var_vec=rep(NA, cycle_no);

obs_prop_mat=obs_prop_var_mat=matrix(NA, cycle_no, mi_no);
raghu_prop_mat=raghu_prop_var_mat=matrix(NA, cycle_no, mi_no);

marginal_mean_mat=marginal_var_mat=obs_mean_mat=obs_var_mat=trans_mean_mat=trans_var_mat=bootcombine_mean_mat=bootcombine_var_mat=draw_mean_mat=draw_var_mat=lambdatrue_mean_mat=lambdatrue_var_mat=boot_mean_mat=boot_var_mat=bayes_mean_mat=bayes_var_mat=matrix(NA, cycle_no, mi_no);

boot_CI_matrix=matrix(NA, cycle_no,2);



# regression parameters;
# both beta0 and beta1 are fixed at 
beta0=-2;
beta1_vec=c(1/sqrt(19), 1/3, 0.5, 1, 2, 3);
# beta1=1; # R-square would be 0.5;
# beta1=0.5; # R-square would be 0.20;
# beta1=2; # R-square would be 0.8;
# beta1=3; # R-square would be 0.9;
# beta1=1/3; # R-square would be 0.1;
# beta1=1/sqrt(19); # R-square would be 0.05;


for (l in 1:6)

{ 

beta1=beta1_vec[l];

cat("the l is", l, "\n");
cat("the beta1 is", beta1, "\n");

for (cycle in 1:cycle_no)

{

set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# error distribution
# could be normal, or could be other distributions;
error=sqrt(var_error)*rnorm(rowobs);

# generate the latent variable;
y2=beta0+beta1*x+error;

# set up the missing data as MCAR on x;

miss_indi=runif(n=rowobs)<0.40;
y2_miss=y2;
y2_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y2_miss[!is.na(y2_miss)]);
y2_obs=y2[miss_indi==0];

# observed data is cbind(x, y2_miss);

# compelete data analysis;
# regression analysis;
reg_com=summary(lm(y2~x));
y2_r2=reg_com$r.squared;

# to obtain the statistics at transformed scale; harel (2009);
Q_com=1/2*log((1+sqrt(y2_r2))/(1-sqrt(y2_r2)));
Q_com_var=1/(rowobs-3);
y2_mean_vec[cycle]=Q_com;
y2_var_vec[cycle]=Q_com_var;


# different missing data methods;

# complete-case analysis;
# mean estimands;
reg_cc=summary(lm(y2_miss~x));
y2_r2_cc=reg_cc$r.squared;
Q_cc=1/2*log((1+sqrt(y2_r2_cc))/(1-sqrt(y2_r2_cc)));
Q_cc_var=1/(obs_no-3);

cc_mean_vec[cycle]=Q_cc;
cc_var_vec[cycle]=Q_cc_var;


# now impute the missing y2's to get estimates of the marginal;
y2_completed=y2_miss;

# boot_results_matrix to hold the bootstrap estimates;
boot_results_matrix=matrix(NA, mi_no, boot_no);

for (i in 1:mi_no)
{

# bivariate model imputation; 

y2_imputed=mice.impute.norm(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=x);

# univariate model imputation
# y2_imputed=mice.impute.norm(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=rep(1,rowobs));


y2_completed[miss_seq]=y2_imputed;

# marginal means;
reg_impute=summary(lm(y2_completed~x));
y2_r2_impute=reg_impute$r.squared;
Q_impute=1/2*log((1+sqrt(y2_r2_impute))/(1-sqrt(y2_r2_impute)));
Q_impute_var=1/(rowobs-3);

# bootstrap;
boot_data=data.frame(y2_completed,x);
boot_results = boot(data=boot_data, statistic=rsq, R=boot_no, formula=y2_completed~x);


# delta method;
obs_mean_mat[cycle,i]=y2_r2_impute;
obs_var_mat[cycle,i]=4*y2_r2_impute*(1-y2_r2_impute)^2/(rowobs-3);

# transformation method;
trans_mean_mat[cycle,i]=Q_impute;
trans_var_mat[cycle,i]=Q_impute_var;

# bootstrap Rubin combine;
bootcombine_mean_mat[cycle,i]=obs_mean_mat[cycle,i];
bootcombine_var_mat[cycle,i]=var(boot_results$t);

# hold bootstrap estimates in the matrix;
boot_results_matrix[i,]=boot_results$t;

}

# obtain the lower and upper percentiles of the bootstrap estimates;
boot_CI_matrix[cycle,]=quantile(boot_results_matrix, c(.025, .975));
# cat("the cycle is", cycle, "\n");

}

# the end of multiple imputation;

# plot the distribution of the data;

# check the performance of the estimates;
# population quantity
# mean estimand;
# true=pnorm((beta0+beta1*mu_x-cut_off)/(beta1^2*var_x+var_error));
y2_mean_vec_tran=((exp(2*y2_mean_vec)-1)/(exp(2*y2_mean_vec)+1))^2;
true=mean(y2_mean_vec_tran);
true;

cat("the true is", true, "\n");

# mean_bias=mean(y2_mean_vec_tran-true);
# mean_bias;

# true_mse=mean((y2_mean_vec_tran-true)^2);
# true_mse;

# lower and upper 95% CI;
# y2_mean_low95_vec=y2_mean_vec-1.96*sqrt(y2_var_vec);
# y2_mean_up95_vec=y2_mean_vec+1.96*sqrt(y2_var_vec);

# y2_mean_low95_vec_tran=((exp(2*y2_mean_low95_vec)-1)/(exp(2*y2_mean_low95_vec)+1))^2
# y2_mean_up95_vec_tran=((exp(2*y2_mean_up95_vec)-1)/(exp(2*y2_mean_up95_vec)+1))^2;


# mean_length=mean(y2_mean_up95_vec_tran-y2_mean_low95_vec_tran);
# mean_length;

# coverage=(y2_mean_low95_vec_tran < true)*(y2_mean_up95_vec_tran > true);

# mean(coverage);



# variance estimates;
# mean(y2_var_vec);
# var(y2_mean_vec);


# complete case analysis;
# cc_mean_vec_tran=((exp(2*cc_mean_vec)-1)/(exp(2*cc_mean_vec)+1))^2;
# cc_true=mean(cc_mean_vec_tran);
# cc_bias=cc_true-true;
# cc_bias;

# cc_mse=mean((cc_mean_vec_tran-true)^2);
# cc_mse;

# lower and upper 95% CI;
# cc_mean_low95_vec=cc_mean_vec-1.96*sqrt(cc_var_vec);
# cc_mean_up95_vec=cc_mean_vec+1.96*sqrt(cc_var_vec);

# cc_mean_low95_vec_tran=((exp(2*cc_mean_low95_vec)-1)/(exp(2*cc_mean_low95_vec)+1))^2
# cc_mean_up95_vec_tran=((exp(2*cc_mean_up95_vec)-1)/(exp(2*cc_mean_up95_vec)+1))^2;

# mean_cc_length=mean(cc_mean_up95_vec_tran-cc_mean_low95_vec_tran);
# mean_cc_length;

# cc_coverage=(cc_mean_low95_vec_tran < true)*(cc_mean_up95_vec_tran > true);

# mean(cc_coverage);

# variance estimates;
# mean(cc_var_vec);
# var(cc_mean_vec);


# multiple imputation estimates;
# marginal transformation imputation;

# observed data transformation imputation;

obs_mi_mean=rep(NA, cycle_no);
obs_mi_var=rep(NA, cycle_no);
obs_mi_df=rep(NA, cycle_no);
obs_mi_f=rep(NA, cycle_no);


trans_mi_mean=rep(NA, cycle_no);
trans_mi_var=rep(NA, cycle_no);
trans_mi_df=rep(NA, cycle_no);
trans_mi_f=rep(NA, cycle_no);

bootcombine_mi_mean=rep(NA, cycle_no);
bootcombine_mi_var=rep(NA, cycle_no);
bootcombine_mi_df=rep(NA, cycle_no);
bootcombine_mi_f=rep(NA, cycle_no);




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

bootcombine_summary=pool.scalar(bootcombine_mean_mat[i,], bootcombine_var_mat[i,], n=rowobs);
bootcombine_mi_mean[i]=bootcombine_summary$qbar;
bootcombine_mi_var[i]=bootcombine_summary$t;
bootcombine_mi_df[i]=bootcombine_summary$df;
bootcombine_mi_f[i]=bootcombine_summary$f;



}

obs_mi_true=mean(obs_mi_mean);
obs_mi_bias=obs_mi_true-true;

obs_mi_bias;

cat("the obs_mi_bias is", obs_mi_bias, "\n");

# obs_mi_mse=mean((obs_mi_mean-true)^2);
# obs_mi_mse;

# coverage;
obs_mean_low95_vec=obs_mi_mean-qt(.975, obs_mi_df)*sqrt(obs_mi_var);
obs_mean_up95_vec=obs_mi_mean+qt(.975, obs_mi_df)*sqrt(obs_mi_var);

# cbind(obs_mean_low95_vec, obs_mean_up95_vec);

mean_obs_length=mean(obs_mean_up95_vec-obs_mean_low95_vec);
mean_obs_length;

cat("the mean_obs_length is", mean_obs_length, "\n");



obs_mi_coverage=(obs_mean_low95_vec < true)*(obs_mean_up95_vec > true);
mean(obs_mi_coverage);


cat("the obs_mi_coveragee is", mean(obs_mi_coverage), "\n");


# variance estimates;
# mean(obs_mi_var);
# var(obs_mi_mean);


# logit transformation;
trans_back_mean=((exp(2*trans_mi_mean)-1)/(exp(2*trans_mi_mean)+1))^2;
trans_mi_true=mean(trans_back_mean);
trans_mi_bias=trans_mi_true-true;

trans_mi_bias;

cat("the trans_mi_bias is", trans_mi_bias, "\n");



# trans_mi_mse=mean((trans_back_mean-true)^2);
# trans_mi_mse;

# coverage;
trans_mean_low95_vec=trans_mi_mean-qt(.975, trans_mi_df)*sqrt(trans_mi_var);
trans_mean_up95_vec=trans_mi_mean+qt(.975, trans_mi_df)*sqrt(trans_mi_var);

trans_back_low95_vec=((exp(2*trans_mean_low95_vec)-1)/(exp(2*trans_mean_low95_vec)+1))^2;
trans_back_up95_vec=((exp(2*trans_mean_up95_vec)-1)/(exp(2*trans_mean_up95_vec)+1))^2;

# cbind(trans_back_low95_vec, trans_back_up95_vec);

mean_trans_length=mean(trans_back_up95_vec-trans_back_low95_vec);
mean_trans_length;


cat("the mean_trans_length is", mean_trans_length, "\n");



trans_mi_coverage=(trans_back_low95_vec < true)*(trans_back_up95_vec > true);
mean(trans_mi_coverage);


cat("the trans_mi_coverage is", mean(trans_mi_coverage), "\n");



# variance estimates;
# mean(trans_mi_var);
# var(trans_mi_mean);


# bootstrap first and then combine;

bootcombine_mi_true=mean(bootcombine_mi_mean);
bootcombine_mi_bias=bootcombine_mi_true-true;

bootcombine_mi_bias;

cat("the bootcombine_mi_bias is", bootcombine_mi_bias, "\n");



# bootcombine_mi_mse=mean((bootcombine_mi_mean-true)^2);
# bootcombine_mi_mse;

# coverage;
bootcombine_mean_low95_vec=bootcombine_mi_mean-qt(.975, bootcombine_mi_df)*sqrt(bootcombine_mi_var);
bootcombine_mean_up95_vec=bootcombine_mi_mean+qt(.975, bootcombine_mi_df)*sqrt(bootcombine_mi_var);

# cbind(obs_mean_low95_vec, obs_mean_up95_vec);

mean_bootcombine_length=mean(bootcombine_mean_up95_vec-bootcombine_mean_low95_vec);
mean_bootcombine_length;


cat("the mean_bootcombine_length is", mean_bootcombine_length, "\n");



bootcombine_mi_coverage=(bootcombine_mean_low95_vec < true)*(bootcombine_mean_up95_vec > true);
mean(bootcombine_mi_coverage);


cat("the bootcombine_mi_coverage is", mean(bootcombine_mi_coverage), "\n");


# variance estimates;
# mean(bootcombine_mi_var);
# var(bootcombine_mi_mean);

# bootstrap confidence interval method;
# bootci_mean_low95_vec=boot_CI_matrix[,1];
# bootci_mean_up95_vec=boot_CI_matrix[,2];

# mean_bootci_length=mean(bootci_mean_up95_vec-bootci_mean_low95_vec);
# mean_bootci_length;

# bootci_mi_coverage=(bootci_mean_low95_vec < true)*(bootci_mean_up95_vec > true);
# mean(bootci_mi_coverage);

}


# plots;




par(mfrow=c(2,2));
hist(obs_mean_mat, main="Original Scale R-square");
qqnorm(obs_mean_mat, main="Original Scale R-square");
qqline(obs_mean_mat, main="Original Scale R-square");
hist(trans_mean_mat, main="Transformed Scale R-square");
qqnorm(trans_mean_mat, main="Transformed Scale R-square");
qqline(trans_mean_mat, main="Transformed Scale R-square");

date();




