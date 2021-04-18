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

# measurement error model;
c0=0;
c1=1;
c2=1;

# c0=-0.2;
# c1=1.2;
# c2=1;

cycle_no=1000;
mi_no=80;

# matrices holding the parameter estimates;
# complete-data inferences;
# xyslope means regressing x on y;
# BD means using true y;
BD_mean_vec=BD_mean_var_vec=BD_slope_vec=BD_slope_var_vec=BD_xyslope_vec=BD_xyslope_var_vec=rep(NA, cycle_no);
# CC means using y from validation data;
CC_mean_vec=CC_mean_var_vec=CC_slope_vec=CC_slope_var_vec=CC_xyslope_vec=CC_xyslope_var_vec=rep(NA, cycle_no);
# SI means using mismeasued y;
SI_mean_vec=SI_mean_var_vec=SI_slope_vec=SI_slope_var_vec=SI_xyslope_vec=SI_xyslope_var_vec=rep(NA, cycle_no);
# MI means using mltiply imputed y;
MI_mean_mat=MI_mean_var_mat=MI_slope_mat=MI_slope_var_mat=MI_xyslope_mat=MI_xyslope_var_mat=matrix(NA, cycle_no, mi_no);

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# error distribution
# could be normal, or could be other distributions;
error=sqrt(var_error)*rnorm(rowobs);

# generate the latent variable;
y_true=beta0+beta1*x+error;

# y is observed after applying some error;
y=c0+c1*y_true+c2*rnorm(rowobs);

# 
# set up the missing data as MCAR on x;
miss_indi=runif(n=rowobs)<0.80;

# set up the missing data as MAR on x;
# alpha0=-1.4;
# alpha1=1;
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x)/(1+exp(alpha0+alpha1*x)));


obs_indi=1-miss_indi;

y_miss=y_true;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_true_obs=y_miss[miss_indi==0];
x_obs=x[miss_indi==0];
x_mis=x[miss_indi==1];

# compelete data (using true y) analysis;
# marginal mean;
BD_mean_vec[cycle]=mean(y_true);
BD_mean_var_vec[cycle]=var(y_true)/rowobs;


# regression analysis;
BD=summary(lm(y_true~x));

BD_slope_vec[cycle]=BD$coeff[2,1];
BD_slope_var_vec[cycle]=BD$coeff[2,2]^2;

# reverse regression analysis;
BD_xy=summary(lm(x~y_true));

BD_xyslope_vec[cycle]=BD_xy$coeff[2,1];
BD_xyslope_var_vec[cycle]=BD_xy$coeff[2,2]^2;


# different missing data methods;

# complete-case analysis;

# marginal mean;
CC_mean_vec[cycle]=mean(y_true_obs);
CC_mean_var_vec[cycle]=var(y_true_obs)/obs_no;


# regression analysis;
CC=summary(lm(y_true_obs~x_obs));

CC_slope_vec[cycle]=CC$coeff[2,1];
CC_slope_var_vec[cycle]=CC$coeff[2,2]^2;

# reverse regression analysis;
CC_xy=summary(lm(x_obs~y_true_obs));

CC_xyslope_vec[cycle]=CC_xy$coeff[2,1];
CC_xyslope_var_vec[cycle]=CC_xy$coeff[2,2]^2;

# using mismeasured y method;

# marginal mean;
SI_mean_vec[cycle]=mean(y);
SI_mean_var_vec[cycle]=var(y)/rowobs;

# regression analysis;

IMP=summary(lm(y~x));

SI_slope_vec[cycle]=IMP$coeff[2,1];
SI_slope_var_vec[cycle]=IMP$coeff[2,2]^2;

# reverse regression analysis;

IMP_xy=summary(lm(x~y));

SI_xyslope_vec[cycle]=IMP_xy$coeff[2,1];
SI_xyslope_var_vec[cycle]=IMP_xy$coeff[2,2]^2;


# multiple imputation analysis;

# now impute the missing y_true's to get estimates of the marginal;
y_completed_MI=y_miss;

# the matrix holds the multiply imputed data;
# x_new is the observed covariates;
x_new=cbind(x, y);
y_completed_MI_mat=matrix(NA, nrow=rowobs, ncol=mi_no);

for (i in 1:mi_no)
{
# proper model imputation;
set.seed(i);

y_imputed_MI=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x_new);


y_completed_MI[miss_seq]=y_imputed_MI;

y_completed_MI_mat[,i]=y_completed_MI;

# completed-data analysis;

# marginal mean;
MI_mean_mat[cycle,i]=mean(y_completed_MI);
MI_mean_var_mat[cycle,i]=var(y_completed_MI)/rowobs;

# regression analysis;

IMP=summary(lm(y_completed_MI~x));

MI_slope_mat[cycle,i]=IMP$coeff[2,1];
MI_slope_var_mat[cycle,i]=IMP$coeff[2,2]^2;

# reverse regression analysis;

IMP_xy=summary(lm(x~y_completed_MI));

MI_xyslope_mat[cycle,i]=IMP_xy$coeff[2,1];
MI_xyslope_var_mat[cycle,i]=IMP_xy$coeff[2,2]^2;


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


BD=performance(BD_mean_vec, BD_mean_var_vec, true_mean);
BD;
CC=performance(CC_mean_vec, CC_mean_var_vec, true_mean);
CC;
SI=performance(SI_mean_vec, SI_mean_var_vec, true_mean);
SI;

MI=mi_performance(MI_mean_mat, MI_mean_var_mat, true_mean, rowobs);
MI;



##########################################################################
# slope estimand;

# complete-data analysis;

true_slope=mean(BD_slope_vec);
true_slope;

BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;
CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;
SI=performance(SI_slope_vec, SI_slope_var_vec, true_slope);
SI;

MI=mi_performance(MI_slope_mat, MI_slope_var_mat, true_slope, rowobs);
MI;


#######################################################################
# reverse slope;

# complete-data analysis;

true_xyslope=mean(BD_xyslope_vec);
true_xyslope;

BD=performance(BD_xyslope_vec, BD_xyslope_var_vec, true_xyslope);
BD;
CC=performance(CC_xyslope_vec, CC_xyslope_var_vec, true_xyslope);
CC;
SI=performance(SI_xyslope_vec, SI_xyslope_var_vec, true_xyslope);
SI;

MI=mi_performance(MI_xyslope_mat, MI_xyslope_var_mat, true_xyslope, rowobs);
MI;




################################################################################# 

# plot;

plot(y_true[1:250], x[1:250], xlab="true y", xlim=c(-5.5, 4), ylab="x");
BD_fit=lm(x[1:250]~y_true[1:250]);
abline(BD_fit);

plot(y[1:250], x[1:250], ylab="x", xlim=c(-5.5, 4), pch=16, xlab="mismeasured y");
SI_fit=lm(x[1:250]~y[1:250]);
abline(SI_fit);



# plot the observed and imputed points;
y_plot_obs=y_miss;

y_plot_imputed=y_completed_MI;
y_plot_imputed[miss_indi==0]=NA;
y_plot_impute=cbind(y_plot_obs, y_plot_imputed);
x_1=x;
x_1[miss_indi==1]=NA;
x_2=x;
x_2[miss_indi==0]=NA;
x_plot_impute=cbind(x_1, x_2);


matplot(y_completed_MI[1:250], x_plot_impute[1:250,], pch=c(1,17), xlim=c(-5.5, 4), ylab="x", xlab="completed y");
impute_fit=lm(x[1:250]~y_completed_MI[1:250]);
abline(impute_fit);


cbind(y_miss[1:250], y[1:250], x[1:250]);



# matplot(y_plot_impute[1:250,], x[1:250], pch=c(1,17), xlab="y", ylab="x");

# plot multiply imputed datasets;
# the 1st imputation;
y_plot_MI_1=y_completed_MI_mat[,1];
y_plot_MI_1[miss_indi==0]=NA;
y_plot_MI_1=cbind(y_plot_obs, y_plot_MI_1);
matplot(x[1:250], y_plot_MI_1[1:250,], pch=c(1,17), xlab="x", ylab="y");
