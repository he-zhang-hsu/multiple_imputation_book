# Example 14.4

rm(list=ls());
# library(car);
# library(mvtnorm);
library(mice);
library(norm);
library(survey);
library(sampling);
# library(HI);
options(digits=4);

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
frac=mean(mi_f);
estimand=as.data.frame(cbind(bias,rbias,mse,mean_length,mean_cov,sd,se,frac))
}

# function mi_fmi() calculate the fraction of missing information;

mi_fmi=function(mean_mat, var_mat, true, n)
{
cycle_no=nrow(mean_mat);
mi_mean=rep(NA, cycle_no);
mi_var=rep(NA, cycle_no);
mi_df=rep(NA, cycle_no);
mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
summary=pool.scalar(mean_mat[i,], var_mat[i,], n=n);
mi_f[i]=summary$f;
}
mi_f;

}

# function mi_deff() calculate the design effect using final MI estimates;
# for proportions;

mi_deff=function(mean_mat, var_mat, true, n)
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

deff=mi_var/(mi_mean*(1-mi_mean)/n);
deff;
}

 
# set up the finite population;
# N is the final population size;
N=100000;

# rowobs is the sample size;     
rowobs=1000;

# set up the random seed;
set.seed(197789);

# mean and variance of z;
mu_z=1;
var_z=1;

# regression parameters;
beta0=-2;
beta1=1;

# c is the cut-off value for estimating the proportions;
c=-1;

# error variance;
var_error=1;

z=mu_z+sqrt(var_z)*rnorm(N);

# error distribution
error=sqrt(var_error)*rnorm(N);

# generate the outcome variable;
y=beta0+beta1*z+error;

# generate selection probability;
alpha0=-1.4;
alpha1=1;

# selection probability proportional to z; 
sample_prob=exp(alpha0+alpha1*z)/(1+exp(alpha0+alpha1*z));

# normalize the selection probability;
true_sample_prob=rowobs*sample_prob/sum(sample_prob);

# sampling weights;
weight=1/true_sample_prob;



# number of simulations;

cycle_no=1000;

# number of imputations;

mi_no=30;

# Vectors and matrices holding estimates for proportions and design effects;
 
BD_mean_vec=BD_mean_var_vec=CC_mean_vec=CC_mean_var_vec=UW_mean_vec=UW_mean_var_vec=rep(NA, cycle_no);
obsno_vec=rep(NA, cycle_no);
BD_UW_mean_vec=BD_UW_mean_var_vec=rep(NA, cycle_no);
BD_deff_vec=CC_deff_vec=rep(NA, cycle_no);
MI_mean_mat=MI_mean_var_mat=MI_deff_mat=matrix(NA, cycle_no, mi_no);
MIW_mean_mat=MIW_mean_var_mat=MIW_deff_mat=matrix(NA, cycle_no, mi_no);
MIWE_mean_mat=MIWE_mean_var_mat=MIWE_deff_mat=matrix(NA, cycle_no, mi_no);
MI_Z_mean_mat=MI_Z_mean_var_mat=MI_Z_deff_mat=matrix(NA, cycle_no, mi_no);
MI_W_mean_mat=MI_W_mean_var_mat=MI_W_deff_mat=matrix(NA, cycle_no, mi_no);
MI_P_mean_mat=MI_P_mean_var_mat=MI_P_deff_mat=matrix(NA, cycle_no, mi_no);
MI_ZW_mean_mat=MI_ZW_mean_var_mat=MI_ZW_deff_mat=matrix(NA, cycle_no, mi_no);
MI_LOGW_mean_mat=MI_LOGW_mean_var_mat=MI_LOGW_deff_mat=matrix(NA, cycle_no, mi_no);



# begin the simulation;

for (cycle in 1:cycle_no)

{
set.seed(cycle);

# sample y and z;
indicator=sample(seq(1,N,1), replace=FALSE, size=rowobs, prob=true_sample_prob);
sum(weight[indicator])/N;

z_sample=z[indicator];
 
y_sample=y[indicator];

# focus on the proportions;
ylec_sample=(y_sample<=c);

wt_sample=weight[indicator];

# compelete data analysis;
bd.data=as.data.frame(cbind(y_sample, ylec_sample, wt_sample));
bd=svydesign(id=~1, weights=wt_sample, data=bd.data);
# q.bd=svymean(~y_sample, bd, deff=TRUE);
q.bd=svymean(~ylec_sample, bd, deff=TRUE);
BD_mean_vec[cycle]=coef(q.bd);
BD_mean_var_vec[cycle]=SE(q.bd)^2;
BD_deff_vec[cycle]=deff(q.bd);


# unweighted data analysis;
# treat them as a simple random sample;
BD_UW_mean_vec[cycle]=mean(ylec_sample);
BD_UW_mean_var_vec[cycle]=var(ylec_sample)/rowobs;

# assign missing values;
# MCAR;
miss_indi=runif(n=rowobs)<0.30;

y_miss=y_sample;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
obsno_vec[cycle]=obs_no;
mis_no=rowobs-obs_no;
y_obs=y_sample[miss_indi==0];
ylec_obs=(y_obs<=c);
wt_sample_obs=wt_sample[miss_indi==0]*rowobs/obs_no;

# complete-case analysis;
cc.data=as.data.frame(cbind(y_obs, ylec_obs, wt_sample_obs));
cc=svydesign(id=~1, weights=wt_sample_obs, data=cc.data);
# q.cc=svymean(~y_obs, cc, deff=TRUE);
q.cc=svymean(~ylec_obs, cc, deff=TRUE);
obs_no_eff=obs_no/deff(q.cc);

CC_mean_vec[cycle]=coef(q.cc);
CC_mean_var_vec[cycle]=SE(q.cc)^2;
CC_deff_vec[cycle]=deff(q.cc);


# unweighted analysis;
UW_mean_vec[cycle]=mean(ylec_obs);
UW_mean_var_vec[cycle]=var(ylec_obs)/obs_no;


# normal imputation;
# obtain observed-data sufficient statistics;
mu_hat=mean(y_obs);
s_square=1/(obs_no-1)*crossprod(y_obs-mu_hat);

# obtain weighted observed-data sufficient statistics;
mu_hat_w=q.cc
s_square_w=obs_no/(obs_no-1)*sum(wt_sample_obs*(y_obs-mu_hat_w)^2)/(sum(wt_sample_obs)-1);

for (i in 1:mi_no)
{
# normal model imputation;
set.seed(i);
# draw sigma_square;
sigma_square=1/(rgamma(1, shape=obs_no/2-1/2, scale=2/((obs_no-1)*s_square)));
mu=rnorm(1, mu_hat, sd=sqrt(sigma_square/obs_no));
y_imputed=rnorm(mis_no, mu, sd=sqrt(sigma_square));
y_completed=y_miss;
y_completed[miss_indi==1]=y_imputed;

ylec_completed=(y_completed <= c);

# marginal mean;
mi.data=as.data.frame(cbind(y_completed, ylec_completed, wt_sample));
mi=svydesign(id=~1, weights=wt_sample, data=mi.data);
# q.mi=svymean(~y_completed, mi, deff=TRUE);
q.mi=svymean(~ylec_completed, mi, deff=TRUE);
MI_mean_mat[cycle,i]=coef(q.mi);
MI_mean_var_mat[cycle,i]=SE(q.mi)^2;
MI_deff_mat[cycle,i]=deff(q.mi);


# normal model imputation but using weighted estimate;
set.seed(i);
# draw sigma_square;
sigma_square_w=1/(rgamma(1, shape=obs_no/2-1/2, scale=2/((obs_no-1)*s_square_w)));
mu_w=rnorm(1, mu_hat_w, sd=sqrt(sigma_square_w/obs_no));
y_imputed_w=rnorm(mis_no, mu_w, sd=sqrt(sigma_square_w));
y_completed_w=y_miss;
y_completed_w[miss_indi==1]=y_imputed_w;

ylec_completed_w=(y_completed_w <= c);
# marginal mean;
mi_w.data=as.data.frame(cbind(y_completed_w, ylec_completed_w, wt_sample));
mi_w=svydesign(id=~1, weights=wt_sample, data=mi_w.data);
# q.mi_w=svymean(~y_completed_w, mi_w, deff=TRUE);
q.mi_w=svymean(~ylec_completed_w, mi_w, deff=TRUE);
MIW_mean_mat[cycle,i]=coef(q.mi_w);
MIW_mean_var_mat[cycle,i]=SE(q.mi_w)^2;
MIW_deff_mat[cycle,i]=deff(q.mi_w);


# normal model imputation but using weighted estimate and using effective sample size;
set.seed(i);
# draw sigma_square;
sigma_square_we=1/(rgamma(1, shape=obs_no_eff/2-1/2, scale=2/((obs_no_eff-1)*s_square_w)));
mu_we=rnorm(1, mu_hat_w, sd=sqrt(sigma_square_we/obs_no_eff));
y_imputed_we=rnorm(mis_no, mu_we, sd=sqrt(sigma_square_we));
y_completed_we=y_miss;
y_completed_we[miss_indi==1]=y_imputed_we;

ylec_completed_we=(y_completed_we <=c);
# marginal mean;
mi_we.data=as.data.frame(cbind(y_completed_we, ylec_completed_we, wt_sample));
mi_we=svydesign(id=~1, weights=wt_sample, data=mi_we.data);
# q.mi_we=svymean(~y_completed_we, mi_we, deff=TRUE);
q.mi_we=svymean(~ylec_completed_we, mi_we, deff=TRUE);
MIWE_mean_mat[cycle,i]=coef(q.mi_we);
MIWE_mean_var_mat[cycle,i]=SE(q.mi_we)^2;
MIWE_deff_mat[cycle,i]=deff(q.mi_we);


# linear model imputation using design;
y_imputed_z=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=z_sample);
y_completed_z=y_miss;
y_completed_z[miss_indi==1]=y_imputed_z;

ylec_completed_z=(y_completed_z <=c);

# marginal mean;
mi_z.data=as.data.frame(cbind(y_completed_z, ylec_completed_z, wt_sample));
mi_z=svydesign(id=~1, weights=wt_sample, data=mi_z.data);
q.mi_z=svymean(~ylec_completed_z, mi_z, deff=TRUE);
MI_Z_mean_mat[cycle,i]=coef(q.mi_z);
MI_Z_mean_var_mat[cycle,i]=SE(q.mi_z)^2;
MI_Z_deff_mat[cycle,i]=deff(q.mi_z);


# linear model imputation using weights as the predictor;
y_imputed_wt=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=wt_sample);
y_completed_wt=y_miss;
y_completed_wt[miss_indi==1]=y_imputed_wt;

ylec_completed_wt=(y_completed_wt <=c);

# marginal mean;
mi_wt.data=as.data.frame(cbind(y_completed_wt, ylec_completed_wt, wt_sample));
mi_wt=svydesign(id=~1, weights=wt_sample, data=mi_wt.data);
# q.mi_wt=svymean(~y_completed_wt, mi_wt, deff=TRUE);
q.mi_wt=svymean(~ylec_completed_wt, mi_wt, deff=TRUE);
MI_W_mean_mat[cycle,i]=coef(q.mi_wt);
MI_W_mean_var_mat[cycle,i]=SE(q.mi_wt)^2;
MI_W_deff_mat[cycle,i]=deff(q.mi_wt);


# linear model imputation using log-weight as the predictor;

y_imputed_logw=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=log(wt_sample));
y_completed_logw=y_miss;
y_completed_logw[miss_indi==1]=y_imputed_logw;

ylec_completed_logw=(y_completed_logw <=c);

# marginal mean;
mi_logw.data=as.data.frame(cbind(y_completed_logw, ylec_completed_logw, wt_sample));
mi_logw=svydesign(id=~1, weights=wt_sample, data=mi_logw.data);
# q.mi_logw=svymean(~y_completed_logw, mi_logw, deff=TRUE);
q.mi_logw=svymean(~ylec_completed_logw, mi_logw, deff=TRUE);
MI_LOGW_mean_mat[cycle,i]=coef(q.mi_logw);
MI_LOGW_mean_var_mat[cycle,i]=SE(q.mi_logw)^2;
MI_LOGW_deff_mat[cycle,i]=deff(q.mi_logw);

# linear model imputation using selection probability as the predictor;
y_imputed_p=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=1/wt_sample);
y_completed_p=y_miss;
y_completed_p[miss_indi==1]=y_imputed_p;

ylec_completed_p=(y_completed_p <=c);

# marginal mean;
mi_p.data=as.data.frame(cbind(y_completed_p, ylec_completed_p, wt_sample));
mi_p=svydesign(id=~1, weights=wt_sample, data=mi_p.data);
# q.mi_p=svymean(~y_completed_p, mi_p, deff=TRUE);
q.mi_p=svymean(~ylec_completed_p, mi_p, deff=TRUE);
MI_P_mean_mat[cycle,i]=coef(q.mi_p);
MI_P_mean_var_mat[cycle,i]=SE(q.mi_p)^2;
MI_P_deff_mat[cycle,i]=deff(q.mi_p);


# linear model imputation using both design and weights as the predictor;
y_imputed_zw=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(z_sample, wt_sample));
y_completed_zw=y_miss;
y_completed_zw[miss_indi==1]=y_imputed_zw;

ylec_completed_zw=(y_completed_zw <=c);

# marginal mean;
mi_zw.data=as.data.frame(cbind(y_completed_zw, ylec_completed_zw, wt_sample));
mi_zw=svydesign(id=~1, weights=wt_sample, data=mi_zw.data);
q.mi_zw=svymean(~ylec_completed_zw, mi_zw, deff=TRUE);
MI_ZW_mean_mat[cycle,i]=coef(q.mi_zw);
MI_ZW_mean_var_mat[cycle,i]=SE(q.mi_zw)^2;
MI_ZW_deff_mat[cycle,i]=deff(q.mi_zw);


# y_completed_MI[miss_seq]=y_imputed_MI;

# y_completed_MI_mat[,i]=y_completed_MI;

# completed-data analysis;



}




cat("the cycle is", cycle, "\n");


}

# End of multiple imputation;

######################################################################
# Simulation results for all methods;
# Table 14.5 only includes the MI-Z imputation method

true_mean=mean(y<=c);


true_mean;

BD=performance(BD_mean_vec, BD_mean_var_vec, true_mean);
BD;
CC=performance(CC_mean_vec, CC_mean_var_vec, true_mean);
CC;
UW=performance(UW_mean_vec, UW_mean_var_vec, true_mean);
UW;
MI=mi_performance(MI_mean_mat, MI_mean_var_mat, true_mean, rowobs);
MI;
MIW=mi_performance(MIW_mean_mat, MIW_mean_var_mat, true_mean, rowobs);
MIW;
MIWE=mi_performance(MIWE_mean_mat, MIWE_mean_var_mat, true_mean, rowobs);
MIWE;
MI_Z=mi_performance(MI_Z_mean_mat, MI_Z_mean_var_mat, true_mean, rowobs);
MI_Z;
MI_W=mi_performance(MI_W_mean_mat, MI_W_mean_var_mat, true_mean, rowobs);
MI_W;
MI_LOGW=mi_performance(MI_LOGW_mean_mat, MI_LOGW_mean_var_mat, true_mean, rowobs);
MI_LOGW;
MI_P=mi_performance(MI_P_mean_mat, MI_P_mean_var_mat, true_mean, rowobs);
MI_P;
MI_ZW=mi_performance(MI_ZW_mean_mat, MI_ZW_mean_var_mat, true_mean, rowobs);
MI_ZW;

# design effect;
# effective sample size;

mean(BD_deff_vec);
mean(rowobs/BD_deff_vec);
mean(CC_deff_vec);
mean(obsno_vec/CC_deff_vec);

# within-imputation design effect;
# and effective sample size;

mean(MI_deff_mat);
mean(rowobs/MI_deff_mat);
mean(MIW_deff_mat);
mean(rowobs/MIW_deff_mat);
mean(MIWE_deff_mat);
mean(rowobs/MIWE_deff_mat);
mean(MI_Z_deff_mat);
mean(rowobs/MI_Z_deff_mat);
mean(MI_W_deff_mat);
mean(rowobs/MI_W_deff_mat);
mean(MI_LOGW_deff_mat);
mean(rowobs/MI_LOGW_deff_mat);
mean(MI_P_deff_mat);
mean(rowobs/MI_P_deff_mat);
mean(MI_ZW_deff_mat);
mean(rowobs/MI_ZW_deff_mat);


# MI design effect;
# and effective sample size;
mean(mi_deff(MI_mean_mat, MI_mean_var_mat, true_mean, rowobs));
mean(rowobs/mi_deff(MI_mean_mat, MI_mean_var_mat, true_mean, rowobs));

mean(mi_deff(MIW_mean_mat, MIW_mean_var_mat, true_mean, rowobs));
mean(rowobs/mi_deff(MIW_mean_mat, MIW_mean_var_mat, true_mean, rowobs));

mean(mi_deff(MIWE_mean_mat, MIWE_mean_var_mat, true_mean, rowobs));
mean(rowobs/mi_deff(MIWE_mean_mat, MIWE_mean_var_mat, true_mean, rowobs));

mean(mi_deff(MI_Z_mean_mat, MI_Z_mean_var_mat, true_mean, rowobs));
mean(rowobs/mi_deff(MI_Z_mean_mat, MI_Z_mean_var_mat, true_mean, rowobs));

mean(mi_deff(MI_W_mean_mat, MI_W_mean_var_mat, true_mean, rowobs));
mean(rowobs/mi_deff(MI_W_mean_mat, MI_W_mean_var_mat, true_mean, rowobs));

mean(mi_deff(MI_LOGW_mean_mat, MI_LOGW_mean_var_mat, true_mean, rowobs));
mean(rowobs/mi_deff(MI_LOGW_mean_mat, MI_LOGW_mean_var_mat, true_mean, rowobs));

mean(mi_deff(MI_P_mean_mat, MI_P_mean_var_mat, true_mean, rowobs));
mean(rowobs/mi_deff(MI_P_mean_mat, MI_P_mean_var_mat, true_mean, rowobs));


mean(mi_deff(MI_ZW_mean_mat, MI_ZW_mean_var_mat, true_mean, rowobs));
mean(rowobs/mi_deff(MI_ZW_mean_mat, MI_ZW_mean_var_mat, true_mean, rowobs));


# within-imputation design effect adjusted by fraction of missing information;
# and effective sample size;

mean(1/(1-mi_fmi(MI_mean_mat, MI_mean_var_mat, true_mean, rowobs))
*rowMeans(MI_deff_mat));

mean(rowobs/(1/(1-mi_fmi(MI_mean_mat, MI_mean_var_mat, true_mean, rowobs))
*rowMeans(MI_deff_mat)));


mean(1/(1-mi_fmi(MIW_mean_mat, MIW_mean_var_mat, true_mean, rowobs))
*rowMeans(MIW_deff_mat));

mean(rowobs/(1/(1-mi_fmi(MIW_mean_mat, MIW_mean_var_mat, true_mean, rowobs))
*rowMeans(MIW_deff_mat)));


mean(1/(1-mi_fmi(MIWE_mean_mat, MIWE_mean_var_mat, true_mean, rowobs))
*rowMeans(MIWE_deff_mat));

mean(rowobs/(1/(1-mi_fmi(MIWE_mean_mat, MIWE_mean_var_mat, true_mean, rowobs))
*rowMeans(MIWE_deff_mat)));



mean(1/(1-mi_fmi(MI_Z_mean_mat, MI_Z_mean_var_mat, true_mean, rowobs))
*rowMeans(MI_Z_deff_mat));

mean(rowobs/(1/(1-mi_fmi(MI_Z_mean_mat, MI_Z_mean_var_mat, true_mean, rowobs))
*rowMeans(MI_Z_deff_mat)));


mean(1/(1-mi_fmi(MI_W_mean_mat, MI_W_mean_var_mat, true_mean, rowobs))
*rowMeans(MI_W_deff_mat));

mean(rowobs/(1/(1-mi_fmi(MI_W_mean_mat, MI_W_mean_var_mat, true_mean, rowobs))
*rowMeans(MI_W_deff_mat)));

mean(1/(1-mi_fmi(MI_LOGW_mean_mat, MI_LOGW_mean_var_mat, true_mean, rowobs))
*rowMeans(MI_LOGW_deff_mat));

mean(rowobs/(1/(1-mi_fmi(MI_LOGW_mean_mat, MI_LOGW_mean_var_mat, true_mean, rowobs))
*rowMeans(MI_LOGW_deff_mat)));


mean(1/(1-mi_fmi(MI_P_mean_mat, MI_P_mean_var_mat, true_mean, rowobs))
*rowMeans(MI_P_deff_mat));

mean(rowobs/(1/(1-mi_fmi(MI_P_mean_mat, MI_P_mean_var_mat, true_mean, rowobs))
*rowMeans(MI_P_deff_mat)));



mean(1/(1-mi_fmi(MI_ZW_mean_mat, MI_ZW_mean_var_mat, true_mean, rowobs))
*rowMeans(MI_ZW_deff_mat));

mean(rowobs/(1/(1-mi_fmi(MI_ZW_mean_mat, MI_ZW_mean_var_mat, true_mean, rowobs))
*rowMeans(MI_ZW_deff_mat)));


