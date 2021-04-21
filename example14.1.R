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

# auxiliary variable z;
mu_z=1;
var_z=3;

# regression parameters;
# both beta0, beta1 and beta2 are fixed at 
beta0=-2;
# beta1=0;
beta1=1;

# error variance;
var_error=1;


cycle_no=1000;
mi_no=50;

# matrices holding the parameter estimates;
# complete-data inferences;
# xyslope means regressing x on y;
BD_mean_vec=BD_mean_var_vec=BD_slope_vec=BD_slope_var_vec=BD_xyslope_vec=BD_xyslope_var_vec=rep(NA, cycle_no);
CC_mean_vec=CC_mean_var_vec=CC_slope_vec=CC_slope_var_vec=CC_xyslope_vec=CC_xyslope_var_vec=rep(NA, cycle_no);
SI_mean_vec=SI_mean_var_vec=SI_slope_vec=SI_slope_var_vec=SI_xyslope_vec=SI_xyslope_var_vec=rep(NA, cycle_no);
INDI_xyslope_vec=INDI_xyslope_var_vec=rep(NA, cycle_no);

# inclusive imputation;
MI_mean_mat=MI_mean_var_mat=MI_slope_mat=MI_slope_var_mat=MI_xyslope_mat=MI_xyslope_var_mat=matrix(NA, cycle_no, mi_no);
# exclusive imputation;
EX_mean_mat=EX_mean_var_mat=EX_slope_mat=EX_slope_var_mat=EX_xyslope_mat=EX_xyslope_var_mat=matrix(NA, cycle_no, mi_no);



for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# error distribution
# could be normal, or could be other distributions;
error=sqrt(var_error)*rnorm(rowobs);

# generate the outcome variable;
y=beta0+beta1*x+error;


# generaite auxiliary variable z;
# not related to y nor related to the missingness mechanism;
# z=mu_z+sqrt(var_z)*rnorm(rowobs);

# related to Y but not related to the missingness mechanism;
error_z=sqrt(var_error)*rnorm(rowobs);

#z=2+y+error_z;

z=-1+y+y^2+error_z;

# Z not related to Y but related to missingness mechanism;

# z=mu_z+sqrt(var_z)*rnorm(rowobs);

# alpha0=-1.7;
# alpha1=1
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*z)/(1+exp(alpha0+alpha1*z)));

# Z related to Y and also related to missingness mechainsim;
# error_z=sqrt(var_error)*rnorm(rowobs);
# z=2+y+error_z;
# alpha0=-1.7;
# alpha1=1
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*z)/(1+exp(alpha0+alpha1*z)));



# set up the missing data as MCAR on x;
miss_indi=runif(n=rowobs)<0.40;

# set up the missing data as MAR on x;
#alpha0=-1.4;
#alpha1=1;
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x)/(1+exp(alpha0+alpha1*x)));

obs_indi=1-miss_indi;

y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
x_obs=x[miss_indi==0];
x_mis=x[miss_indi==1];
z_obs=z[miss_indi==0];
z_mis=z[miss_indi==1];

# compelete data analysis;
# marginal mean;
BD_mean_vec[cycle]=mean(y);
BD_mean_var_vec[cycle]=var(y)/rowobs;


# regression analysis;
# coefficient on z;
# BD=summary(lm(y~z+x));

BD=summary(lm(y~x));


BD_slope_vec[cycle]=BD$coeff[2,1];
BD_slope_var_vec[cycle]=BD$coeff[2,2]^2;


# different missing data methods;

# complete-case analysis;

# marginal mean;
CC_mean_vec[cycle]=mean(y_obs);
CC_mean_var_vec[cycle]=var(y_obs)/obs_no;


# regression analysis;
# CC=summary(lm(y_obs~z_obs+x_obs));

CC=summary(lm(y_obs~x_obs));


CC_slope_vec[cycle]=CC$coeff[2,1];
CC_slope_var_vec[cycle]=CC$coeff[2,2]^2;



# multiple imputation analysis;

# now impute the missing y's to get estimates of the marginal;
y_completed_MI=y_completed_EX=y_miss;

# the matrix holds the multiply imputed data;

y_completed_MI_mat=y_completed_EX_mat=matrix(NA, nrow=rowobs, ncol=mi_no);

for (i in 1:mi_no)
{
# inclusive imputation;

set.seed(i);

y_imputed_MI=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x,z));

y_completed_MI[miss_seq]=y_imputed_MI;

y_completed_MI_mat[,i]=y_completed_MI;

# completed-data analysis;

# marginal mean;
MI_mean_mat[cycle,i]=mean(y_completed_MI);
MI_mean_var_mat[cycle,i]=var(y_completed_MI)/rowobs;

# regression analysis;

# IMP=summary(lm(y_completed_MI~z+x));

IMP=summary(lm(y_completed_MI~x));

MI_slope_mat[cycle,i]=IMP$coeff[2,1];
MI_slope_var_mat[cycle,i]=IMP$coeff[2,2]^2;

# exclusive imputation;


y_imputed_EX=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);


y_completed_EX[miss_seq]=y_imputed_EX;

y_completed_EX_mat[,i]=y_completed_EX;

# completed-data analysis;

# marginal mean;
EX_mean_mat[cycle,i]=mean(y_completed_EX);
EX_mean_var_mat[cycle,i]=var(y_completed_EX)/rowobs;

# regression analysis;

# IMP=summary(lm(y_completed_EX~z+x));

IMP=summary(lm(y_completed_EX~x));



EX_slope_mat[cycle,i]=IMP$coeff[2,1];
EX_slope_var_mat[cycle,i]=IMP$coeff[2,2]^2;


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

MI=mi_performance(MI_mean_mat, MI_mean_var_mat, true_mean, rowobs);
MI;

EX=mi_performance(EX_mean_mat, EX_mean_var_mat, true_mean, rowobs);
EX;




##########################################################################
# slope estimand;

# complete-data analysis;

true_slope=mean(BD_slope_vec);
true_slope;

BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;
CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;


MI=mi_performance(MI_slope_mat, MI_slope_var_mat, true_slope, rowobs);
MI;

EX=mi_performance(EX_slope_mat, EX_slope_var_mat, true_slope, rowobs);
EX;


################################################################################# 

# plot;

y_plot_obs=y;
y_plot_miss=y;
y_plot_obs[miss_indi==1]=NA;
y_plot_miss[miss_indi==0]=NA;
cbind(y_plot_obs, y_plot_miss);

# plot(x,y);

y_plot=cbind(y_plot_obs, y_plot_miss);

# plot the observed and missing points;

matplot(x[1:250], y_plot[1:250,], pch=c(1,16), xlab="x", ylab="y");
# matplot(y_plot[1:250,], x[1:250], pch=c(1,16), xlab="y", ylab="x");

# plot the observed and imputed points;

y_plot_imputed=y_completed;
y_plot_imputed[miss_indi==0]=NA;
y_plot_impute=cbind(y_plot_obs, y_plot_imputed);

matplot(x[1:250], y_plot_impute[1:250,], pch=c(1,17), xlab="x", ylab="y");
# matplot(y_plot_impute[1:250,], x[1:250], pch=c(1,17), xlab="y", ylab="x");

# plot multiply imputed datasets;
# the 1st imputation;
y_plot_MI_1=y_completed_MI_mat[,1];
y_plot_MI_1[miss_indi==0]=NA;
y_plot_MI_1=cbind(y_plot_obs, y_plot_MI_1);
matplot(x[1:250], y_plot_MI_1[1:250,], pch=c(1,17), xlab="x", ylab="y");

# the 2nd imputation;
y_plot_MI_2=y_completed_MI_mat[,2];
y_plot_MI_2[miss_indi==0]=NA;
y_plot_MI_2=cbind(y_plot_obs, y_plot_MI_2);
matplot(x[1:250], y_plot_MI_2[1:250,], pch=c(1,17), xlab="x", ylab="y");

# the 3rd imputation;
y_plot_MI_3=y_completed_MI_mat[,3];
y_plot_MI_3[miss_indi==0]=NA;
y_plot_MI_3=cbind(y_plot_obs, y_plot_MI_3);
matplot(x[1:250], y_plot_MI_3[1:250,], pch=c(1,17), xlab="x", ylab="y");


