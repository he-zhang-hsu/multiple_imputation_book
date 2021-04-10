# example 3.9
# the improper multiple imputations;

library(mice);
library(norm);
# library(HI);
options(digits=4);
rm(list=ls());


set.seed(197789);


rowobs=1000;
cycle_no=1000;

# mi_no=2;

# mi_no=5;

# mi_no=20;

mi_no=50;

# mi_no=200;

# M=1000;

# M=10000;

mu_y=1;
var_y=1;

BD_mean_vec=BD_mean_var_vec=rep(NA, cycle_no);
CC_mean_vec=CC_mean_var_vec=rep(NA, cycle_no);
MI_mean_mat=MI_mean_var_mat=matrix(NA, cycle_no, mi_no);

# IMP are the improper multiple imputations;
IMP_MI_mu_mean_mat=IMP_MI_mu_mean_var_mat=matrix(NA, cycle_no, mi_no);
IMP_MI_sigma_mean_mat=IMP_MI_sigma_mean_var_mat=matrix(NA, cycle_no, mi_no);
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

# compelete data analysis;
# marginal mean;
BD_mean_vec[cycle]=mean(y);
BD_mean_var_vec[cycle]=var(y)/rowobs;



# complete case analysis;
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
mean(sqrt(BD_mean_var_vec));
sqrt(var(BD_mean_vec));


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
mean(sqrt(CC_mean_var_vec));
sqrt(var(CC_mean_vec));

# multiple imputation estimates;

mean_mi_mean=rep(NA, cycle_no);
mean_mi_var=rep(NA, cycle_no);
mean_mi_df=rep(NA, cycle_no);
mean_mi_f=rep(NA, cycle_no);



for (i in 1:cycle_no)
{
# mean_summary=pool.scalar(MI_mean_mat[i,], MI_mean_var_mat[i,], n=rowobs);
# for improper imputations;
# mean_summary=pool.scalar(IMP_MI_mu_mean_mat[i,], IMP_MI_mu_mean_var_mat[i,], n=rowobs);
# mean_summary=pool.scalar(IMP_MI_sigma_mean_mat[i,], IMP_MI_sigma_mean_var_mat[i,], n=rowobs);
mean_summary=pool.scalar(IMP_MI_mu_sigma_mean_mat[i,], IMP_MI_mu_sigma_mean_var_mat[i,], n=rowobs);



mean_mi_mean[i]=mean_summary$qbar;
mean_mi_var[i]=mean_summary$t;
mean_mi_df[i]=mean_summary$df;
mean_mi_f[i]=mean_summary$f;


}


# For marginal means;

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
mean(sqrt(mean_mi_var));
sqrt(var(mean_mi_mean));




# analyze a single imputation;
SI_mean_vec=MI_mean_mat[,1];
SI_mean_var_vec=MI_mean_var_mat[,1];


SI_mean_bias=mean(SI_mean_vec)-true_mean;
SI_mean_bias;

SI_mean_mse=mean((SI_mean_vec-true_mean)^2);
SI_mean_mse;

# lower and upper 95% CI;
SI_mean_low95_vec=SI_mean_vec-1.96*sqrt(SI_mean_var_vec);
SI_mean_up95_vec=SI_mean_vec+1.96*sqrt(SI_mean_var_vec);

SI_mean_length=mean(SI_mean_up95_vec-SI_mean_low95_vec);
SI_mean_length;

SI_mean_cov=(SI_mean_low95_vec < true_mean)*(SI_mean_up95_vec > true_mean);

mean(SI_mean_cov);

# variance estimates;
mean(sqrt(SI_mean_var_vec));
sqrt(var(SI_mean_vec));
