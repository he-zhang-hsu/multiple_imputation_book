# Example 13.10

rm(list=ls());
# library(car);
# library(mvtnorm);
library(mice);
library(norm);
# library(HI);
options(digits=4);

# generate complete data;
# complete-data sample size;
rowobs=1000;

# set up the random seed;
set.seed(197789);

# generate the response indicator;
R=rbinom(n=rowobs, size=1, prob=0.5);
obs_no=sum(R);
mis_no=rowobs-obs_no;

# pattern mixture model parameters;
mu_1=1;
var_1=1;

mu_0=1.5;
var_0=1.5;

# complete data;

y_obs=rnorm(obs_no, mu_1, sd=sqrt(var_1));
y_mis=rnorm(mis_no, mu_0, sd=sqrt(var_0)); 

y_com=rep(NA, rowobs);

y_com[R==1]=y_obs;
y_com[R==0]=y_mis;

CC_mean=mean(y_obs);
CC_var=var(y_obs)/obs_no;

CC_mean_low95=CC_mean-1.96*sqrt(CC_var);
CC_mean_up95=CC_mean+1.96*sqrt(CC_var);

CC_mean
CC_mean_low95;
CC_mean_up95;



BD_mean=mean(y_com);
BD_var=var(y_com)/rowobs;

BD_mean
BD_mean_low95=BD_mean-1.96*sqrt(BD_var);
BD_mean_up95=BD_mean+1.96*sqrt(BD_var);

BD_mean_low95;
BD_mean_up95;


# multiple imputation sensitivity analysis;
# missing data vector is y_delete;

y_delete=y_com;
y_delete[R==0]=NA;

# Start multiple imputation;
M=50;
# The vectors containing sensitivity parameters;
c_mu_vec=c(-0.5, -0.25, 0, 0.25, 0.5);
c_sigma_vec=c(0.25, 0.5, 1, 2, 4);

# Results from Table 13.5

for (l in 1:length(c_mu_vec))
{
for (s in 1:length(c_sigma_vec))
{

c_mu=c_mu_vec[l];
c_sigma=c_sigma_vec[s];

y_completed=y_delete;

MI_mean_vec=MI_mean_var_vec=rep(NA, M);

for (i in 1:M)
{
 set.seed(i);
# draw parameters from the model;
# sufficient statistics;
  mu_hat=mean(y_obs);
  s_square=1/(obs_no-1)*crossprod(y_obs-mu_hat);

  sigma1_square=1/(rgamma(1, shape=obs_no/2-1/2, scale=2/((obs_no-1)*s_square)));
  mu1=rnorm(n=1, mean=mu_hat, sd=sqrt(sigma1_square/obs_no));
  
  # infer about the nonresponse parameter;
  sigma0_square=c_sigma^2*sigma1_square;
  mu0=mu1*(1+c_mu);

  # draw imputations;
  y_impute=rnorm(n=mis_no, mean=mu0, sd=sqrt(sigma0_square));
  y_completed[R==0]=y_impute;

  # estimate the mean;
  MI_mean_vec[i]=mean(y_completed);
  MI_mean_var_vec[i]=var(y_completed)/rowobs;


}
 
# MI_mean_vec;

# combine the multiple imputation estimates;
mean_summary=pool.scalar(MI_mean_vec, MI_mean_var_vec, n=rowobs);
mean_mi_mean=mean_summary$qbar;
mean_mi_var=mean_summary$t;
mean_mi_df=mean_summary$df;
mean_mi_f=mean_summary$f;

MI_mean_low95=mean_mi_mean-qt(.975, mean_mi_df)*sqrt(mean_mi_var);
MI_mean_up95=mean_mi_mean+qt(.975, mean_mi_df)*sqrt(mean_mi_var);

cat("the c_mu is", c_mu, "\n");
cat("the c_sigma is", c_sigma, "\n");



mean_mi_mean;
MI_mean_low95;
MI_mean_up95;

cat("the mean_mi_mean is", mean_mi_mean, "\n");
cat("the lower bound is", MI_mean_low95, "\n");
cat("the upper bound is", MI_mean_up95, "\n");




}

}

