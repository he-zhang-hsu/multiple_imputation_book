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
library(MASS);
# library(HI);
options(digits=4);
rm(list=ls());

date();


parameter_draw=function(data,m)
{
 # obtain all the computational components;
 # m is the number of Gibbs draws;
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
 sigma_square=1/(rgamma(m, shape=n_total/2-p/2, scale=2/((n_total-p)*s_square)));
 
 # draw beta;
 beta=matrix(NA,m,2);
 for(i in 1:m)
 {
 beta[i,]=mvrnorm(n=1, mu=beta_hat, Sigma=sigma_square[i]*inv_v_beta);
 }
 parameter=cbind(beta,sigma_square);
 parameter;
}

# function to obtain R-Squared from the data 
rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
} 

# function to obtain coefficients and residual sigma square from the data
reg_bootstrap=function(formula, data, indices) {
  d = data[indices,]; # allows boot to select sample 
  fit = lm(formula, data=d);
  beta_estimate=fit$coef;
  anova_fit=anova(fit);
  sigma2_estimate=anova_fit["Residuals", "Mean Sq"];
 # return(c(beta_estimate, sigma2_estimate))
  return(beta_estimate[2])
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

# regression parameters;
# both beta0 and beta1 are fixed at 
beta0=-2;
beta1=1; # R-square would be 0.5;
# beta1=0.5; # R-square would be 0.20;
# beta1=2; # R-square would be 0.8;
# beta1=3; # R-square would be 0.9;
# beta1=1/3; # R-square would be 0.1;
# beta1=1/sqrt(19); # R-square would be 0.05;


# error variance;
var_error=1;

miss_rate=0.40;

cycle_no=1000;
boot_no=200;
gibbs_no=200;

mi_no_vec=c(5,10,20,50,100);

for (l in 1:length(mi_no_vec))

{

mi_no=mi_no_vec[l]
cat("the l is", l, "\n");
cat("the mi_no is", mi_no, "\n");



# mi_no=5;
# mi_no=10;
# mi_no=20;
# mi_no=50;
# boot_no=1000;
# gibbs_no=1000;


# matrices holding the parameter estimates;
# complete-data inferences;
y2_mean_vec=y2_var_vec=cc_mean_vec=cc_var_vec=gibbs_mean_vec=gibbs_var_vec=boot_mean_vec=boot_var_vec=rep(NA, cycle_no);
y2_cc_prop_raghu_vec=y2_cc_prop_raghu_var_vec=y2_prop_raghu_vec=y2_prop_raghu_var_vec=rep(NA, cycle_no);

obs_prop_mat=obs_prop_var_mat=matrix(NA, cycle_no, mi_no);
raghu_prop_mat=raghu_prop_var_mat=matrix(NA, cycle_no, mi_no);

marginal_mean_mat=marginal_var_mat=obs_mean_mat=obs_var_mat=trans_mean_mat=trans_var_mat=bootcombine_mean_mat=bootcombine_var_mat=draw_mean_mat=draw_var_mat=lambdatrue_mean_mat=lambdatrue_var_mat=boot_mean_mat=boot_var_mat=bayes_mean_mat=bayes_var_mat=matrix(NA, cycle_no, mi_no);

gibbs_CI_matrix=matrix(NA, cycle_no,2);
boot_CI_matrix=matrix(NA, cycle_no,2);




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


miss_indi=runif(n=rowobs)< miss_rate;
y2_miss=y2;
y2_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y2_miss[!is.na(y2_miss)]);
y2_obs=y2[miss_indi==0];

# observed data is cbind(x, y2_miss);

# compelete data analysis;
# regression analysis;
reg_com=summary(lm(y2~x));

# the slope coefficients;
x_coeff=reg_com$coefficients[2,1];
x_coeff_var=reg_com$coefficients[2,2]^2;

y2_mean_vec[cycle]=x_coeff;
y2_var_vec[cycle]=x_coeff_var;


# different missing data methods;

# complete-case analysis;
# mean estimands;
reg_cc=summary(lm(y2_miss~x));
x_coeff_cc=reg_cc$coefficients[2,1];
x_coeff_cc_var=reg_cc$coefficients[2,2]^2;

cc_mean_vec[cycle]=x_coeff_cc;
cc_var_vec[cycle]=x_coeff_cc_var;


# now impute the missing y2's to get estimates of the marginal;
y2_completed=y2_miss;

# gibbs_results_matrix to hold the bayes estimates;
gibbs_results_matrix=matrix(NA, mi_no, gibbs_no);

# boot_results_matrix to hold the bootstrap estimates;
boot_results_matrix=matrix(NA, mi_no, boot_no);

for (i in 1:mi_no)
{

# bivariate model imputation; 

y2_imputed=mice.impute.norm(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=x);

# univariate model imputation
# y2_imputed=mice.impute.norm(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=rep(1,rowobs));


y2_completed[miss_seq]=y2_imputed;

# Rubin's imputation approach;
reg_impute=summary(lm(y2_completed~x));
x_coeff_impute=reg_impute$coefficients[2,1];
x_coeff_impute_var=reg_impute$coefficients[2,2]^2;

obs_mean_mat[cycle,i]=x_coeff_impute;
obs_var_mat[cycle,i]=x_coeff_impute_var;

# Bayesian estimate;
gibbs_data=cbind(x, y2_completed);
gibbs_draws=parameter_draw(gibbs_data, gibbs_no);
gibbs_results_matrix[i,]=gibbs_draws[,2];


# bootstrap;
boot_data=data.frame(y2_completed,x);
boot_results = boot(data=boot_data, statistic=reg_bootstrap, R=boot_no, formula=y2_completed~x);
boot_results_matrix[i,]=boot_results$t;

}


gibbs_mean_vec[cycle]=mean(gibbs_results_matrix);
gibbs_var_vec[cycle]=var(as.vector(gibbs_results_matrix));
gibbs_CI_matrix[cycle,]=quantile(gibbs_results_matrix, c(.025, .975));


boot_mean_vec[cycle]=mean(boot_results_matrix);
boot_var_vec[cycle]=var(as.vector(boot_results_matrix));
boot_CI_matrix[cycle,]=quantile(boot_results_matrix, c(.025, .975));

# obtain the lower and upper percentiles of the Bayes and bootstrap estimates;
# cat("the cycle is", cycle, "\n");

}

# the end of multiple imputation;

# plot the distribution of the data;

# check the performance of the estimates;
# population quantity
# mean estimand;
true=mean(y2_mean_vec);
true;

cat("the true is", true, "\n");

# mean_bias=mean(y2_mean_vec-true);
# mean_bias;

# true_mse=mean((y2_mean_vec-true)^2);
# true_mse;

# lower and upper 95% CI;
# y2_mean_low95_vec=y2_mean_vec-1.96*sqrt(y2_var_vec);
# y2_mean_up95_vec=y2_mean_vec+1.96*sqrt(y2_var_vec);

# mean_length=mean(y2_mean_up95_vec-y2_mean_low95_vec);
# mean_length;

# coverage=(y2_mean_low95_vec < true)*(y2_mean_up95_vec > true);

# mean(coverage);



# variance estimates;
# mean(y2_var_vec);
# var(y2_mean_vec);


# complete case analysis;
# cc_true=mean(cc_mean_vec);
# cc_bias=cc_true-true;
# cc_bias;

# cc_mse=mean((cc_mean_vec-true)^2);
# cc_mse;

# lower and upper 95% CI;
# cc_mean_low95_vec=cc_mean_vec-1.96*sqrt(cc_var_vec);
#  cc_mean_up95_vec=cc_mean_vec+1.96*sqrt(cc_var_vec);

# mean_cc_length=mean(cc_mean_up95_vec-cc_mean_low95_vec);
# mean_cc_length;

# cc_coverage=(cc_mean_low95_vec < true)*(cc_mean_up95_vec > true);

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


for (i in 1:cycle_no)
{
obs_summary=pool.scalar(obs_mean_mat[i,], obs_var_mat[i,], n=rowobs);
obs_mi_mean[i]=obs_summary$qbar;
obs_mi_var[i]=obs_summary$t;
obs_mi_df[i]=obs_summary$df;
obs_mi_f[i]=obs_summary$f;


}

obs_mi_true=mean(obs_mi_mean);
obs_mi_bias=obs_mi_true-true;

obs_mi_bias;

cat("the obs_mi_bias is", obs_mi_bias, "\n");

obs_mi_mse=mean((obs_mi_mean-true)^2);
obs_mi_mse;

# coverage;
obs_mean_low95_vec=obs_mi_mean-qt(.975, obs_mi_df)*sqrt(obs_mi_var);
obs_mean_up95_vec=obs_mi_mean+qt(.975, obs_mi_df)*sqrt(obs_mi_var);

# cbind(obs_mean_low95_vec, obs_mean_up95_vec);

mean_obs_length=mean(obs_mean_up95_vec-obs_mean_low95_vec);
mean_obs_length;

cat("the mean_obs_length is", mean_obs_length, "\n");

obs_mi_coverage=(obs_mean_low95_vec < true)*(obs_mean_up95_vec > true);
mean(obs_mi_coverage);

cat("the obs_mi_coverage is", mean(obs_mi_coverage), "\n");

# variance estimates;
mean(sqrt(obs_mi_var));
sqrt(var(obs_mi_mean));

cat("the SE is", mean(sqrt(obs_mi_var)), "\n");
cat("the SD is", sqrt(var(obs_mi_mean)), "\n");

# Bayes confidence interval method;
gibbs_true=mean(gibbs_mean_vec);
gibbs_bias=gibbs_true-true;
gibbs_bias;

cat("the gibbs_bias is", gibbs_bias, "\n");

gibbsci_mean_low95_vec=gibbs_CI_matrix[,1];
gibbsci_mean_up95_vec=gibbs_CI_matrix[,2];

mean_gibbsci_length=mean(gibbsci_mean_up95_vec-gibbsci_mean_low95_vec);
mean_gibbsci_length;

cat("the mean_gibbsci_length is", mean_gibbsci_length, "\n");

gibbsci_mi_coverage=(gibbsci_mean_low95_vec < true)*(gibbsci_mean_up95_vec > true);
mean(gibbsci_mi_coverage);

cat("the gibbsci_mi_coverage", mean(gibbsci_mi_coverage), "\n");

# variance estimates;
mean(sqrt(gibbs_var_vec));
sqrt(var(gibbs_mean_vec));

cat("the SE is", mean(sqrt(gibbs_var_vec)), "\n");
cat("the SD is", sqrt(var(gibbs_mean_vec)), "\n");



# bootstrap confidence interval method;
boot_true=mean(boot_mean_vec);
boot_bias=boot_true-true;
boot_bias;

cat("the boot_bias is", boot_bias, "\n");


bootci_mean_low95_vec=boot_CI_matrix[,1];
bootci_mean_up95_vec=boot_CI_matrix[,2];

mean_bootci_length=mean(bootci_mean_up95_vec-bootci_mean_low95_vec);
mean_bootci_length;

cat("the mean_bootci_lenght is", mean_bootci_length, "\n");

bootci_mi_coverage=(bootci_mean_low95_vec < true)*(bootci_mean_up95_vec > true);
mean(bootci_mi_coverage);

cat("the bootci_mi_coverage is", mean(bootci_mi_coverage), "\n");

# variance estimates;
mean(sqrt(boot_var_vec));
sqrt(var(boot_mean_vec));

cat("the SE is", mean(sqrt(boot_var_vec)), "\n");
cat("the SD is", sqrt(var(boot_mean_vec)), "\n");

}

# plots;

# par(mfrow=c(2,2));
# hist(obs_mean_mat, main="Original Scale R-square");
# qqnorm(obs_mean_mat, main="Original Scale R-square");
# qqline(obs_mean_mat, main="Original Scale R-square");
# hist(trans_mean_mat, main="Transformed Scale R-square");
# qqnorm(trans_mean_mat, main="Transformed Scale R-square");
# qqline(trans_mean_mat, main="Transformed Scale R-square");

date();




