# Example 3.8
# The peformance of improper multiple imputations;
rm(list=ls());

library(mice);
library(norm);
# library(HI);
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



set.seed(197789);


rowobs=1000;
cycle_no=1000;


# the number of imputations;

mi_no=50;


mu_y=1;
var_y=1;

# Before-deletion analysis
BD_mean_vec=BD_mean_var_vec=rep(NA, cycle_no);

# Completec-case analysis
CC_mean_vec=CC_mean_var_vec=rep(NA, cycle_no);

# MI is the proper multiple imputation
MI_mean_mat=MI_mean_var_mat=matrix(NA, cycle_no, mi_no);

# IMP are the improper multiple imputations;
# fixing mu
IMP_MI_mu_mean_mat=IMP_MI_mu_mean_var_mat=matrix(NA, cycle_no, mi_no);

# fixing sigma square
IMP_MI_sigma_mean_mat=IMP_MI_sigma_mean_var_mat=matrix(NA, cycle_no, mi_no);

# fixing both mu and sigma square
IMP_MI_mu_sigma_mean_mat=IMP_MI_mu_sigma_mean_var_mat=matrix(NA, cycle_no, mi_no);


for (cycle in 1:cycle_no)
{

set.seed(cycle);

# generate complete data;

y=rnorm(rowobs,mean=mu_y, sd=sqrt(var_y));

# set up the missing data as MCAR ;
miss_indi=runif(n=rowobs)<0.50;


y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];

# Before-deletion analysis;
# marginal mean;
BD_mean_vec[cycle]=mean(y);
BD_mean_var_vec[cycle]=var(y)/rowobs;



# Complete-case analysis;
# marginal mean;
CC_mean_vec[cycle]=mean(y_obs);
CC_mean_var_vec[cycle]=var(y_obs)/obs_no;



# generate imputed data;

# obtain observed-data sufficient statistics;
mu_hat=mean(y_obs);
s_square=1/(obs_no-1)*crossprod(y_obs-mu_hat);

# multiple imputation;
y_completed=y_completed_imp_mu=y_completed_imp_sigma=y_completed_imp_mu_sigma=y_miss;

for (i in 1:mi_no)
{

set.seed(i);

# proper imputation;
# draw sigma_square;
sigma_square=1/(rgamma(1, shape=obs_no/2-1/2, scale=2/((obs_no-1)*s_square)));
# draw mu
mu=rnorm(1, mu_hat, sd=sqrt(sigma_square/obs_no));

y_imputed=rnorm(mis_no, mu, sd=sqrt(sigma_square));


y_completed[miss_seq]=y_imputed;


# completed-data analysis;

# marginal mean;
MI_mean_mat[cycle,i]=mean(y_completed);
MI_mean_var_mat[cycle,i]=var(y_completed)/rowobs;

# improper imputation, fixing mu at mu_hat;
y_imputed_imp_mu=rnorm(mis_no, mu_hat, sd=sqrt(sigma_square));

y_completed_imp_mu[miss_seq]=y_imputed_imp_mu;


# completed-data analysis;

# marginal mean;
IMP_MI_mu_mean_mat[cycle,i]=mean(y_completed_imp_mu);
IMP_MI_mu_mean_var_mat[cycle,i]=var(y_completed_imp_mu)/rowobs;

# improper imputation, fixing sigma at sigma_hat;
y_imputed_imp_sigma=rnorm(mis_no, mu, sd=sqrt(s_square));

y_completed_imp_sigma[miss_seq]=y_imputed_imp_sigma;

# completed-data analysis;

# marginal mean;
IMP_MI_sigma_mean_mat[cycle,i]=mean(y_completed_imp_sigma);
IMP_MI_sigma_mean_var_mat[cycle,i]=var(y_completed_imp_sigma)/rowobs;

# improper imputation, fixing mu and sigma both;

y_imputed_imp_mu_sigma=rnorm(mis_no, mu_hat, sd=sqrt(s_square));

y_completed_imp_mu_sigma[miss_seq]=y_imputed_imp_mu_sigma;


# completed-data analysis;

# marginal mean;
IMP_MI_mu_sigma_mean_mat[cycle,i]=mean(y_completed_imp_mu_sigma);
IMP_MI_mu_sigma_mean_var_mat[cycle,i]=var(y_completed_imp_mu_sigma)/rowobs;



}

cat("the cycle is", cycle, "\n");



}

#############################################################################
# Simulation results for the marginal mean estimand of Y
# Table 3.7

true_mean=mean(BD_mean_vec);
true_mean;

# Before-deletion analysis;

BD=performance(BD_mean_vec, BD_mean_var_vec, true_mean);
BD;

# Complete-case analysis;

CC=performance(CC_mean_vec, CC_mean_var_vec, true_mean);
CC;

# Proper multiple imputation analysis;

MI=mi_performance(MI_mean_mat, MI_mean_var_mat, true_mean, rowobs);
MI;


# Improper multiple imputation analysis; 
# fixing mu
IMP_MI_mu=mi_performance(IMP_MI_mu_mean_mat, IMP_MI_mu_mean_var_mat, true_mean, rowobs);
IMP_MI_mu;

# fixing sigma square
IMP_MI_sigma=mi_performance(IMP_MI_sigma_mean_mat, IMP_MI_sigma_mean_var_mat, true_mean, rowobs);
IMP_MI_sigma;


# fixing both mu and sigma square
IMP_MI_mu_sigma=mi_performance(IMP_MI_mu_sigma_mean_mat, IMP_MI_mu_sigma_mean_var_mat, true_mean, rowobs);
IMP_MI_mu_sigma;
