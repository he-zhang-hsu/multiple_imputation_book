################## data generate ##############
# Program: boxcox_impute_expolore_simulation; ############
# generate data, and try out different missing data methods;
# run some simulations to assess the biasedness of various methods;

# library(car);
# library(mvtnorm);
library(mice);
library(norm);
library(survey);
library(sampling);

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

weight_estimate=function(y, weight, n)
{
y_wt_mean=weighted.mean(y, weight);

# variance of the population;
y_wt_var=n/(n-1)*sum(weight*(y-y_wt_mean)^2)/(sum(weight)-1);

# variance estimate of the weighted mean;
var_y_wt_mean=sum(weight^2)*y_wt_var/sum(weight)^2;

weight_estimate=c(y_wt_mean, var_y_wt_mean);
weight_estimate;

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
frac=mean(mi_f);
estimand=as.data.frame(cbind(bias,rbias,mse,mean_length,mean_cov,sd,se,frac))
}


 
# set up the finite population;

N=100000;
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

# error variance;
var_error=1;

x=mu_x+sqrt(var_x)*rnorm(N);

# error distribution
# could be normal, or could be other distributions;
error=sqrt(var_error)*rnorm(N);

# generate the outcome variable;
y=beta0+beta1*x+error;

# generate selection probability;
alpha0=-1.4;
alpha1=1;

sample_prob=exp(alpha0+alpha1*x)/(1+exp(alpha0+alpha1*x));

true_sample_prob=rowobs*sample_prob/sum(sample_prob);
weight=1/true_sample_prob;


# indicator=sample(seq(1,N,1), size=rowobs, prob=true_sample_prob);
# sum(weight[indicator])/N;

# x_sample=x[indicator];
 
# y_sample=y[indicator];
# wt_sample=weight[indicator];
# wt_total=sum(wt_sample);

# test.data=as.data.frame(cbind(y_sample, wt_sample));
# test=svydesign(id=~1, weights=wt_sample, data=test.data);
# svymean(~y_sample, test);

# weighted mean;
# y_wt_mean=weighted.mean(y_sample, wt_sample);
# y_wt_total=sum(wt_sample);
# var_y_wt_mean=1/((y_wt_total^2)*rowobs*(rowobs-1))*sum((y_sample*wt_sample-y_wt_mean*wt_sample)^2);

# variance of the population;
# y_wt_var=rowobs/(rowobs-1)*sum(wt_sample*(y_sample-y_wt_mean)^2)/(sum(wt_sample)-1);

# e=y_sample-y_wt_mean;

# Hajekestimator(y_sample, 1/wt_sample, type="mean");

# variance estimate of the weighted mean;
# var_y_wt_mean=sum(wt_sample^2)*y_wt_var/sum(wt_sample)^2;

# deff=rowobs*sum(wt_sample^2)/sum(wt_sample)^2;

# plots;
# plot(x_sample, y_sample);
# plot(wt_sample, y_sample);
# plot(1/wt_sample, y_sample);

cycle_no=1000;
mi_no=30;

BD_mean_vec=BD_mean_var_vec=CC_mean_vec=CC_mean_var_vec=UW_mean_vec=UW_mean_var_vec=rep(NA, cycle_no);
BD_UW_mean_vec=BD_UW_mean_var_vec=rep(NA, cycle_no);
BD_deff_vec=CC_deff_vec=rep(NA, cycle_no);
MI_mean_mat=MI_mean_var_mat=MI_deff_mat=matrix(NA, cycle_no, mi_no);
MIW_mean_mat=MIW_mean_var_mat=MIW_deff_mat=matrix(NA, cycle_no, mi_no);
MIWE_mean_mat=MIWE_mean_var_mat=MIWE_deff_mat=matrix(NA, cycle_no, mi_no);
MI_X_mean_mat=MI_X_mean_var_mat=MI_X_deff_mat=matrix(NA, cycle_no, mi_no);
MI_W_mean_mat=MI_W_mean_var_mat=MI_W_deff_mat=matrix(NA, cycle_no, mi_no);
MI_P_mean_mat=MI_P_mean_var_mat=MI_P_deff_mat=matrix(NA, cycle_no, mi_no);
MI_XW_mean_mat=MI_XW_mean_var_mat=MI_XW_deff_mat=matrix(NA, cycle_no, mi_no);
MI_LOGW_mean_mat=MI_LOGW_mean_var_mat=MI_LOGW_deff_mat=matrix(NA, cycle_no, mi_no);





for (cycle in 1:cycle_no)

{
set.seed(cycle);

# sample y and x;
indicator=sample(seq(1,N,1), replace=TRUE, size=rowobs, prob=true_sample_prob);
sum(weight[indicator])/N;

x_sample=x[indicator];
 
y_sample=y[indicator];
wt_sample=weight[indicator];

# compelete data analysis;
bd.data=as.data.frame(cbind(y_sample, wt_sample));
bd=svydesign(id=~1, weights=wt_sample, data=bd.data);
q.bd=svymean(~y_sample, bd, deff=TRUE);
BD_mean_vec[cycle]=coef(q.bd);
BD_mean_var_vec[cycle]=SE(q.bd)^2;
BD_deff_vec[cycle]=deff(q.bd);


# unweighted data analysis;
# treat them as a simple random sample;
BD_UW_mean_vec[cycle]=mean(y_sample);
BD_UW_mean_var_vec[cycle]=var(y_sample)/rowobs;



# assign missing values;
# MCAR;
miss_indi=runif(n=rowobs)<0.30;

y_miss=y_sample;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y_sample[miss_indi==0];
wt_sample_obs=wt_sample[miss_indi==0]*rowobs/obs_no;

# complete-case analysis;
cc.data=as.data.frame(cbind(y_obs, wt_sample_obs));
cc=svydesign(id=~1, weights=wt_sample_obs, data=cc.data);
q.cc=svymean(~y_obs, cc, deff=TRUE);
obs_no_eff=obs_no/deff(q.cc);

CC_mean_vec[cycle]=coef(q.cc);
CC_mean_var_vec[cycle]=SE(q.cc)^2;
CC_deff_vec[cycle]=deff(q.cc);


# unweighted analysis;
UW_mean_vec[cycle]=mean(y_obs);
UW_mean_var_vec[cycle]=var(y_obs)/obs_no;


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

# marginal mean;
mi.data=as.data.frame(cbind(y_completed, wt_sample));
mi=svydesign(id=~1, weights=wt_sample, data=mi.data);
q.mi=svymean(~y_completed, mi, deff=TRUE);
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

# marginal mean;
mi_w.data=as.data.frame(cbind(y_completed_w, wt_sample));
mi_w=svydesign(id=~1, weights=wt_sample, data=mi_w.data);
q.mi_w=svymean(~y_completed_w, mi_w, deff=TRUE);
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

# marginal mean;
mi_we.data=as.data.frame(cbind(y_completed_we, wt_sample));
mi_we=svydesign(id=~1, weights=wt_sample, data=mi_we.data);
q.mi_we=svymean(~y_completed_we, mi_we, deff=TRUE);
MIWE_mean_mat[cycle,i]=coef(q.mi_we);
MIWE_mean_var_mat[cycle,i]=SE(q.mi_we)^2;
MIWE_deff_mat[cycle,i]=deff(q.mi_we);




# linear model imputation using design;
y_imputed_x=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x_sample);
y_completed_x=y_miss;
y_completed_x[miss_indi==1]=y_imputed_x;

# marginal mean;
mi_x.data=as.data.frame(cbind(y_completed_x, wt_sample));
mi_x=svydesign(id=~1, weights=wt_sample, data=mi_x.data);
q.mi_x=svymean(~y_completed_x, mi_x, deff=TRUE);
MI_X_mean_mat[cycle,i]=coef(q.mi_x);
MI_X_mean_var_mat[cycle,i]=SE(q.mi_x)^2;
MI_X_deff_mat[cycle,i]=deff(q.mi_x);


# linear model imputation using weights as the predictor;
y_imputed_wt=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=wt_sample);
y_completed_wt=y_miss;
y_completed_wt[miss_indi==1]=y_imputed_wt;

# marginal mean;
mi_wt.data=as.data.frame(cbind(y_completed_wt, wt_sample));
mi_wt=svydesign(id=~1, weights=wt_sample, data=mi_wt.data);
q.mi_wt=svymean(~y_completed_wt, mi_wt, deff=TRUE);
MI_W_mean_mat[cycle,i]=coef(q.mi_wt);
MI_W_mean_var_mat[cycle,i]=SE(q.mi_wt)^2;
MI_W_deff_mat[cycle,i]=deff(q.mi_wt);


# linear model imputation using log-weight as the predictor;

y_imputed_logw=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=log(wt_sample));
y_completed_logw=y_miss;
y_completed_logw[miss_indi==1]=y_imputed_logw;

# marginal mean;
mi_logw.data=as.data.frame(cbind(y_completed_logw, wt_sample));
mi_logw=svydesign(id=~1, weights=wt_sample, data=mi_logw.data);
q.mi_logw=svymean(~y_completed_logw, mi_logw, deff=TRUE);
MI_LOGW_mean_mat[cycle,i]=coef(q.mi_logw);
MI_LOGW_mean_var_mat[cycle,i]=SE(q.mi_logw)^2;
MI_LOGW_deff_mat[cycle,i]=deff(q.mi_logw);




# linear model imputation using selection probability as the predictor;
y_imputed_p=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=1/wt_sample);
y_completed_p=y_miss;
y_completed_p[miss_indi==1]=y_imputed_p;

# marginal mean;
mi_p.data=as.data.frame(cbind(y_completed_p, wt_sample));
mi_p=svydesign(id=~1, weights=wt_sample, data=mi_p.data);
q.mi_p=svymean(~y_completed_p, mi_p, deff=TRUE);
MI_P_mean_mat[cycle,i]=coef(q.mi_p);
MI_P_mean_var_mat[cycle,i]=SE(q.mi_p)^2;
MI_P_deff_mat[cycle,i]=deff(q.mi_p);


# linear model imputation using both design and weights as the predictor;
y_imputed_xw=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x_sample, wt_sample));
y_completed_xw=y_miss;
y_completed_xw[miss_indi==1]=y_imputed_xw;

# marginal mean;
mi_xw.data=as.data.frame(cbind(y_completed_xw, wt_sample));
mi_xw=svydesign(id=~1, weights=wt_sample, data=mi_xw.data);
q.mi_xw=svymean(~y_completed_xw, mi_xw, deff=TRUE);
MI_XW_mean_mat[cycle,i]=coef(q.mi_xw);
MI_XW_mean_var_mat[cycle,i]=SE(q.mi_xw)^2;
MI_XW_deff_mat[cycle,i]=deff(q.mi_xw);





# y_completed_MI[miss_seq]=y_imputed_MI;

# y_completed_MI_mat[,i]=y_completed_MI;

# completed-data analysis;



}




cat("the cycle is", cycle, "\n");


}

# marginal mean estimand;
true_mean=mean(y);
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
MI_X=mi_performance(MI_X_mean_mat, MI_X_mean_var_mat, true_mean, rowobs);
MI_X;
MI_W=mi_performance(MI_W_mean_mat, MI_W_mean_var_mat, true_mean, rowobs);
MI_W;
MI_LOGW=mi_performance(MI_LOGW_mean_mat, MI_LOGW_mean_var_mat, true_mean, rowobs);
MI_LOGW;
MI_P=mi_performance(MI_P_mean_mat, MI_P_mean_var_mat, true_mean, rowobs);
MI_P;
MI_XW=mi_performance(MI_XW_mean_mat, MI_XW_mean_var_mat, true_mean, rowobs);
MI_XW;

# design effect;

mean(BD_deff_vec);
mean(CC_deff_vec);
mean(MI_deff_mat);
mean(MIW_deff_mat);
mean(MIWE_deff_mat);
mean(MI_X_deff_mat);
mean(MI_W_deff_mat);
mean(MI_LOGW_deff_mat);
mean(MI_P_deff_mat);
mean(MI_XW_deff_mat);



################################################################################# 

# plot;

y_plot_obs=y_sample;
y_plot_miss=y_sample;
y_plot_obs[miss_indi==1]=NA;
y_plot_miss[miss_indi==0]=NA;
cbind(y_plot_obs, y_plot_miss);

# plot(x,y);

y_plot=cbind(y_plot_obs, y_plot_miss);

# plot the observed and missing points;

matplot(x_sample[1:250], y_plot[1:250,], pch=c(1,16), xlab="z", ylab="y");

matplot(wt_sample[1:250], y_plot[1:250,], pch=c(1,16), xlab="weight", ylab="y");

matplot(1/wt_sample[1:250], y_plot[1:250,], pch=c(1,16), xlab="selection prob", ylab="y");



# matplot(y_plot[1:250,], x[1:250], pch=c(1,16), xlab="y", ylab="x");

# unweighed imputation;
# plot the observed and imputed points;

y_plot_imputed=y_completed;
y_plot_imputed[miss_indi==0]=NA;
y_plot_impute=cbind(y_plot_obs, y_plot_imputed);

matplot(x_sample[1:250], y_plot_impute[1:250,], pch=c(1,17), xlab="z", ylab="y");

matplot(wt_sample[1:250], y_plot_impute[1:250,], pch=c(1,17), xlab="weight", ylab="y");

matplot(1/wt_sample[1:250], y_plot_impute[1:250,], pch=c(1,17), xlab="selection prob", ylab="y");

# weighted imputation;

# y_plot_imputed=y_completed_we;
# y_plot_imputed[miss_indi==0]=NA;
# y_plot_impute=cbind(y_plot_obs, y_plot_imputed);

# matplot(x_sample[1:250], y_plot_impute[1:250,], pch=c(1,15), xlab="z", ylab="y");

# matplot(wt_sample[1:250], y_plot_impute[1:250,], pch=c(1,15), xlab="weight", ylab="y");

# matplot(1/wt_sample[1:250], y_plot_impute[1:250,], pch=c(1,15), xlab="selection prob", ylab="y");

# imputation using x;

y_plot_imputed=y_completed_x;
y_plot_imputed[miss_indi==0]=NA;
y_plot_impute=cbind(y_plot_obs, y_plot_imputed);


matplot(x_sample[1:250], y_plot_impute[1:250,], pch=c(1,18), xlab="z", ylab="y");

matplot(wt_sample[1:250], y_plot_impute[1:250,], pch=c(1,18), xlab="weight", ylab="y");

matplot(1/wt_sample[1:250], y_plot_impute[1:250,], pch=c(1,18), xlab="selection prob", ylab="y");



# imputation using probability;

y_plot_imputed=y_completed_p;
y_plot_imputed[miss_indi==0]=NA;
y_plot_impute=cbind(y_plot_obs, y_plot_imputed);

matplot(x_sample[1:250], y_plot_impute[1:250,], pch=c(1,15), col=c("black", "red"), bg="red", xlab="z", ylab="y");

matplot(wt_sample[1:250], y_plot_impute[1:250,], pch=c(1,15), col=c("black", "red"), bg="red", xlab="weight", ylab="y");

matplot(1/wt_sample[1:250], y_plot_impute[1:250,], pch=c(1,15), col=c("black", "red"), bg="red", xlab="selection prob", ylab="y");




# imputation using weight;

y_plot_imputed=y_completed_wt;
y_plot_imputed[miss_indi==0]=NA;
y_plot_impute=cbind(y_plot_obs, y_plot_imputed);

matplot(x_sample[1:250], y_plot_impute[1:250,], pch=c(1,25), col=c("black", "red"), bg="red", xlab="z", ylab="y");

matplot(wt_sample[1:250], y_plot_impute[1:250,], pch=c(1,25), col=c("black", "red"), bg="red", xlab="weight", ylab="y");

matplot(1/wt_sample[1:250], y_plot_impute[1:250,], pch=c(1,25), col=c("black", "red"), bg="red", xlab="selection prob", ylab="y");



