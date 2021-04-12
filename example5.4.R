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





# set up the random seed;
set.seed(197789);





# generate data;
# distribution of x, the covariate;
# could be normal, or could be other distributions;
# fixed at mean 10, and variance 1;
rowobs=1000;
mu_x=5;
var_x=1;
cycle_no=1000;
mi_no=30;


# linear normal model;
# both beta0, beta1 and beta2 are fixed at 
beta0=10;
# beta1=0;
beta1=1;
beta2=1;
# error variance;
var_error=1;


# generate complete data and get their percentiles;
# x1=mu_x+sqrt(var_x)*rnorm(rowobs*cycle_no);
# x2=mu_x+sqrt(var_x)*rnorm(rowobs*cycle_no);

# error distribution
# could be normal, or could be other distributions;
# error=sqrt(var_error)*rnorm(rowobs*cycle_no);

# generate the latent variable;
# y=beta0+beta1*x1+beta2*x2+error;

# transformation; from schenker et al. (1996)
# both beta0, beta1 and beta2 are fixed at 
# beta0=0;
# beta1=1;
# beta2=1;
# error variance;
# var_error=1;


# generate complete data and get their percentiles;
# x1=mu_x+sqrt(var_x)*rnorm(rowobs*cycle_no);
# x2=mu_x+sqrt(var_x)*rnorm(rowobs*cycle_no);

# error distribution
# could be normal, or could be other distributions;
# error=sqrt(var_error)*rnorm(rowobs*cycle_no);

# generate the latent variable via transformation;
# y=(beta0+beta1*x1+beta2*x2+error)^4;


# generate the population;
y2_complete_mat=matrix(NA, rowobs, cycle_no);



for (cycle in 1:cycle_no)

{
set.seed(cycle);

x1=mu_x+sqrt(var_x)*rnorm(rowobs);
x2=mu_x+sqrt(var_x)*rnorm(rowobs);


# error distribution
# could be normal, or could be other distributions;
error=sqrt(var_error)*rnorm(rowobs);

# generate the latent variable via linear normal model;
y=beta0+beta1*x1+beta2*x2+error;

# generate the latent variable via transformation;
# y=(beta0+beta1*x1+beta2*x2+error)^4;
y2_complete_mat[,cycle]=y;

}





# get the percentiles;
# quantile(y, c(0.05, 0.25, 0.50, 0.75, 0.95));
y_cutoff=quantile(y2_complete_mat, 0.05);


# matrices holding the parameter estimates;
# complete-data inferences;
# xyslope means regressing x on y;
BD_mean_vec=BD_mean_var_vec=BD_slope_vec=BD_slope_var_vec=BD_xyslope_vec=BD_xyslope_var_vec=rep(NA, cycle_no);
CC_mean_vec=CC_mean_var_vec=CC_slope_vec=CC_slope_var_vec=CC_xyslope_vec=CC_xyslope_var_vec=rep(NA, cycle_no);
MI_mean_mat=MI_mean_var_mat=MI_slope_mat=MI_slope_var_mat=MI_xyslope_mat=MI_xyslope_var_mat=matrix(NA, cycle_no, mi_no);
# IM_MI indicates PMM MI
IM_MI_mean_mat=IM_MI_mean_var_mat=IM_MI_slope_mat=IM_MI_slope_var_mat=IM_MI_xyslope_mat=IM_MI_xyslope_var_mat=matrix(NA, cycle_no, mi_no);



for (cycle in cycle_no:cycle_no)

{
set.seed(cycle);

x1=mu_x+sqrt(var_x)*rnorm(rowobs);
x2=mu_x+sqrt(var_x)*rnorm(rowobs);


# error distribution
# could be normal, or could be other distributions;
error=sqrt(var_error)*rnorm(rowobs);

# generate the latent variable via linear normal model;
y=beta0+beta1*x1+beta2*x2+error;

# generate the latent variable via transformation;
# y=(beta0+beta1*x1+beta2*x2+error)^4;



# set up the missing data as MCAR on x;
miss_indi=runif(n=rowobs)<0.30;

# set up the missing data as MAR on x;
# alpha0=-1.4;
# alpha1=1;
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x)/(1+exp(alpha0+alpha1*x)));


y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
# x_obs=x[miss_indi==0];
# x_mis=x[miss_indi==1];

# compelete data analysis;
# marginal mean;
BD_mean_vec[cycle]=mean(y);
BD_mean_var_vec[cycle]=var(y)/rowobs;


# proportion analysis;
BD_slope_vec[cycle]=mean(y<y_cutoff);
BD_slope_var_vec[cycle]=var(y<y_cutoff)/rowobs;

# different missing data methods;

# complete-case analysis;

# marginal mean;
CC_mean_vec[cycle]=mean(y_obs);
CC_mean_var_vec[cycle]=var(y_obs)/obs_no;


# proportion analysis;
CC_slope_vec[cycle]=mean(y_obs<y_cutoff);
CC_slope_var_vec[cycle]=var(y_obs<y_cutoff)/obs_no;


# now impute the missing y2's to get estimates of the marginal;
y_completed=y_completed_IM=y_miss;

for (i in 1:mi_no)
{
set.seed(i);
# normal model imputation;

y_imputed=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x1,x2));


y_completed[miss_seq]=y_imputed;

# completed-data analysis;

# marginal mean;
MI_mean_mat[cycle,i]=mean(y_completed);
MI_mean_var_mat[cycle,i]=var(y_completed)/rowobs;

# proportion analysis;

MI_slope_mat[cycle,i]=mean(y_completed<y_cutoff);
MI_slope_var_mat[cycle,i]=var(y_completed<y_cutoff)/rowobs;

# PMM imputation;
y_imputed_IM=mice.impute.pmm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x1,x2));
y_completed_IM[miss_seq]=y_imputed_IM;

# completed-data analysis;

# marginal mean;
IM_MI_mean_mat[cycle,i]=mean(y_completed_IM);
IM_MI_mean_var_mat[cycle,i]=var(y_completed_IM)/rowobs;

# proportion analysis;

IM_MI_slope_mat[cycle,i]=mean(y_completed_IM < y_cutoff);
IM_MI_slope_var_mat[cycle,i]=var(y_completed_IM < y_cutoff)/rowobs;



}

cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

# plot the distribution of the data;

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
mean(sqrt(BD_mean_var_vec));
sqrt(var(BD_mean_vec));


# complete case analysis;

CC_mean_bias=mean(CC_mean_vec)-true_mean;
CC_mean_bias;
CC_mean_bias/true_mean;

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

mean_IM_mi_mean=rep(NA, cycle_no);
mean_IM_mi_var=rep(NA, cycle_no);
mean_IM_mi_df=rep(NA, cycle_no);
mean_IM_mi_f=rep(NA, cycle_no);



for (i in 1:cycle_no)
{
mean_summary=pool.scalar(MI_mean_mat[i,], MI_mean_var_mat[i,], n=rowobs);
mean_mi_mean[i]=mean_summary$qbar;
mean_mi_var[i]=mean_summary$t;
mean_mi_df[i]=mean_summary$df;
mean_mi_f[i]=mean_summary$f;

IM_mean_summary=pool.scalar(IM_MI_mean_mat[i,], IM_MI_mean_var_mat[i,], n=rowobs);
mean_IM_mi_mean[i]=IM_mean_summary$qbar;
mean_IM_mi_var[i]=IM_mean_summary$t;
mean_IM_mi_df[i]=IM_mean_summary$df;
mean_IM_mi_f[i]=IM_mean_summary$f;


}


# For marginal means;
# normal imputation;

mean_mi_bias=mean(mean_mi_mean)-true_mean;
mean_mi_bias;
mean_mi_bias/true_mean;

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

# PMM imputation;

mean_IM_mi_bias=mean(mean_IM_mi_mean)-true_mean;
mean_IM_mi_bias;
mean_IM_mi_bias/true_mean;

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
mean(sqrt(mean_IM_mi_var));
sqrt(var(mean_IM_mi_mean));


##########################################################################
# proportional estimand;

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
slope_summary=pool.scalar(MI_slope_mat[i,], MI_slope_var_mat[i,], n=rowobs);
slope_mi_mean[i]=slope_summary$qbar;
slope_mi_var[i]=slope_summary$t;
slope_mi_df[i]=slope_summary$df;
slope_mi_f[i]=slope_summary$f;

IM_slope_summary=pool.scalar(IM_MI_slope_mat[i,], IM_MI_slope_var_mat[i,], n=rowobs);
slope_IM_mi_mean[i]=IM_slope_summary$qbar;
slope_IM_mi_var[i]=IM_slope_summary$t;
slope_IM_mi_df[i]=IM_slope_summary$df;
slope_IM_mi_f[i]=IM_slope_summary$f;



}



# For proprotions;

# normal imputation;

slope_mi_bias=mean(slope_mi_mean)-true_slope;
slope_mi_bias;
slope_mi_bias/true_slope;

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
mean(sqrt(slope_mi_var));
sqrt(var(slope_mi_mean));

# PMM imputation;

slope_IM_mi_bias=mean(slope_IM_mi_mean)-true_slope;
slope_IM_mi_bias;
slope_IM_mi_bias/true_slope;

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
mean(sqrt(slope_IM_mi_var));
sqrt(var(slope_IM_mi_mean));

# make some plots;
# only for the transformation model;
# for example 5.4;
par(mfrow=c(1,1));
hist(y, freq=T, main="Before-deletion", xlab="y");
hist(y_miss, freq=T, main="Observed", xlab="y");
hist(y_completed, main="Normal Imputation", xlab="y");
hist(y_completed_IM, main="PMM Imputation", xlab="y");



y_plot_obs=y;
y_plot_obs[miss_indi==1]=NA;

x=x1+x2;

# for normal imputation;
y_plot_imputed=y_completed;
y_plot_imputed[miss_indi==0]=NA;
y_plot_impute=cbind(y_plot_obs, y_plot_imputed);

matplot(x, y_plot_impute, pch=c(1,17), xlab="x1+x2", ylab="y", main="Normal Imputation");
abline(h=y_cutoff, col="green");

# for PMM imputation;

y_plot_imputed=y_completed_IM;
y_plot_imputed[miss_indi==0]=NA;
y_plot_impute=cbind(y_plot_obs, y_plot_imputed);

matplot(x, y_plot_impute, pch=c(1,17), xlab="x1+x2", ylab="y", main="PMM Imputation");
abline(h=y_cutoff, col="green");




# x1_obs=x1[miss_indi==0];
# x1_miss=x1[miss_indi==1];

# par(mfrow=c(1,1));
# lw1 <- loess(y_obs ~ x1_obs);
# plot(y_obs ~ x1_obs, main="Scatter plot of observed Y vs. X1");
# j <- order(x1_obs)
# lines(x1_obs[j],lw1$fitted[j],col="red",lwd=3)


# lw1 <- loess(y_imputed ~ x1_miss);
# plot(y_imputed ~ x1_miss, main="Scatter plot of NM imputed Y vs. X1");
# j <- order(x1_miss)
# lines(x1_miss[j],lw1$fitted[j],col="red",lwd=3)

# lw1 <- loess(y_imputed_IM ~ x1_miss);
# plot(y_imputed_IM ~ x1_miss, main="Scatter plot of PMM imputed Y vs. X1");
# j <- order(x1_miss)
# lines(x1_miss[j],lw1$fitted[j],col="red",lwd=3)



# par(mfrow=c(2,2));
# plot(x1_obs, y_obs);
# plot(x1_miss, y_imputed);
# plot(x1_miss, y_imputed_IM);
