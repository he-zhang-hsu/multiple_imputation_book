# Example 14.2

rm(list=ls());
# library(car);
# library(mvtnorm);
library(mice);
library(norm);
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
estimand=as.data.frame(cbind(bias,rbias,mse,mean_length,mean_cov,sd,se))
}

# set up the random seed;
set.seed(197789);

# complete-data sample size;
rowobs=1000;

# mean and variance of x;
mu_x=3;
var_x=1;

# regression parameters; 
beta0=-2;
beta1=1;


# error variance;
var_error=1;

# number of simulations;
cycle_no=1000;

# number of imputations;
mi_no=50;

# cut-off value for the proportions;
cut_off=1;

# cut_off=-2;

# cut_off=3;

# matrices holding the parameter estimates;
# complete-data inferences;
BD_mean_vec=BD_var_vec=CC_mean_vec=CC_var_vec=rep(NA, cycle_no);

# proportions;
BD_prop_vec=BD_prop_var_vec=BD_prop_model_vec=BD_prop_model_var_vec=rep(NA, cycle_no);
CC_prop_vec=CC_prop_var_vec=CC_prop_model_vec=CC_prop_model_var_vec=rep(NA, cycle_no);

# obs means moment estimator;
obs_prop_mat_uni=obs_prop_var_mat_uni=obs_prop_mat_biv=obs_prop_var_mat_biv=matrix(NA, cycle_no, mi_no);

# model means model based estimator;
model_prop_mat_uni=model_prop_var_mat_uni=model_prop_mat_biv=model_prop_var_mat_biv=matrix(NA, cycle_no, mi_no);

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# normal error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the y variable;
y=beta0+beta1*x+error;


# set up the missing data as MCAR on x;


miss_indi=runif(n=rowobs)<0.40;
y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
y_obs=y[miss_indi==0];
mis_no=rowobs-obs_no;

# observed data is cbind(x, y_miss);

# compelete data analysis;
# sample mean;
y_mean=mean(y);
y_var=var(y)/rowobs;
BD_mean_vec[cycle]=y_mean;
BD_var_vec[cycle]=y_var;


# proportions: methods of moment estimator;
y_prop=mean(y<cut_off);
y_prop_var=y_prop*(1-y_prop)/rowobs;
BD_prop_vec[cycle]=y_prop;
BD_prop_var_vec[cycle]=y_prop_var;



# proportions: model-based estimator;
y_prop_model=pnorm((cut_off-y_mean)/sqrt(var(y)))
y_prop_model_var=(dnorm((cut_off-y_mean)/sqrt(var(y))))^2*(1/rowobs+(cut_off-y_mean)^2/(2*var(y)^2*(rowobs-1)));
BD_prop_model_vec[cycle]=y_prop_model;
BD_prop_model_var_vec[cycle]=y_prop_model_var;




# different missing data methods;

# complete-case analysis;
# mean estimands;
CC_mean=mean(y_miss, na.rm="T");
CC_mean_vec[cycle]=CC_mean;
CC_var_vec[cycle]=var(y_miss, na.rm="T")/obs_no;

# proportions; method-of-moment estimator;

CC_prop=mean(y_miss<cut_off, na.rm=T);
CC_prop_var=CC_prop*(1-CC_prop)/obs_no;
CC_prop_vec[cycle]=CC_prop;
CC_prop_var_vec[cycle]=CC_prop_var;

# proportions: model-based estimator;

CC_prop_model=pnorm((cut_off-mean(y_obs))/sqrt(var(y_obs)))
CC_prop_model_var=(dnorm((cut_off-mean(y_obs))/sqrt(var(y_obs))))^2*(1/obs_no+(cut_off-mean(y_obs))^2/(2*var(y_obs)^2*(obs_no-1)));
CC_prop_model_vec[cycle]=CC_prop_model;
CC_prop_model_var_vec[cycle]=CC_prop_model_var;


# now impute the missing y's to get estimates of the marginal;
y_completed_uni=y_completed_biv=y_miss;

mu_hat=mean(y_obs);
s_square=1/(obs_no-1)*crossprod(y_obs-mu_hat);



for (i in 1:mi_no)
{

# univariate model imputation
# draw sigma_square;
sigma_square=1/(rgamma(1, shape=obs_no/2-1/2, scale=2/((obs_no-1)*s_square)));
mu=rnorm(1, mu_hat, sd=sqrt(sigma_square/obs_no));

y_imputed_uni=rnorm(mis_no, mu, sd=sqrt(sigma_square));
y_completed_uni[miss_seq]=y_imputed_uni;

# proportions: method of moment estimator;
obs_prop_uni=mean(y_completed_uni<cut_off);
obs_prop_var_uni=obs_prop_uni*(1-obs_prop_uni)/rowobs;
obs_prop_mat_uni[cycle,i]=obs_prop_uni;
obs_prop_var_mat_uni[cycle,i]=obs_prop_var_uni;

# marginal means;
obs_mean_uni=mean(y_completed_uni);
# obs_var_uni=var(y_completed_uni)/rowobs;

# obs_mean_mat[cycle,i]=obs_mean;
# obs_var_mat[cycle,i]=obs_var;


# proportions: model estimator;
model_prop_uni=pnorm((cut_off-obs_mean_uni)/sqrt(var(y_completed_uni)));
model_prop_var_uni=(dnorm((cut_off-obs_mean_uni)/sqrt(var(y_completed_uni))))^2*(1/rowobs+(cut_off-obs_mean_uni)^2/(2*var(y_completed_uni)^2*(rowobs-1)));
model_prop_mat_uni[cycle,i]=model_prop_uni;
model_prop_var_mat_uni[cycle,i]=model_prop_var_uni;

# bivariate model imputation; 
y_imputed_biv=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);
y_completed_biv[miss_seq]=y_imputed_biv;

# proportions: method of moment estimator;
obs_prop_biv=mean(y_completed_biv<cut_off);
obs_prop_var_biv=obs_prop_biv*(1-obs_prop_biv)/rowobs;
obs_prop_mat_biv[cycle,i]=obs_prop_biv;
obs_prop_var_mat_biv[cycle,i]=obs_prop_var_biv;

# marginal means;
obs_mean_biv=mean(y_completed_biv);
# obs_var_uni=var(y_completed_uni)/rowobs;

# obs_mean_mat[cycle,i]=obs_mean;
# obs_var_mat[cycle,i]=obs_var;


# proportions: model estimator;
model_prop_biv=pnorm((cut_off-obs_mean_biv)/sqrt(var(y_completed_biv)));
model_prop_var_biv=(dnorm((cut_off-obs_mean_biv)/sqrt(var(y_completed_biv))))^2*(1/rowobs+(cut_off-obs_mean_biv)^2/(2*var(y_completed_biv)^2*(rowobs-1)));
model_prop_mat_biv[cycle,i]=model_prop_biv;
model_prop_var_mat_biv[cycle,i]=model_prop_var_biv;




}

cat("the cycle is", cycle, "\n");



}

# the end of multiple imputation;

###############################################################
# Simulation results;
# For proportion;
# Table 14.2

true_prop=mean(BD_prop_vec);

BD_moment=performance(BD_prop_vec, BD_prop_var_vec, true_prop);
BD_moment;

BD_efficient=performance(BD_prop_model_vec, BD_prop_model_var_vec, true_prop);
BD_efficient;



CC_moment=performance(CC_prop_vec, CC_prop_var_vec, true_prop);
CC_moment;


CC_efficient=performance(CC_prop_model_vec, CC_prop_model_var_vec, true_prop);
CC_efficient;


MI_uni_moment=mi_performance(obs_prop_mat_uni, obs_prop_var_mat_uni, true_prop, rowobs);
MI_uni_moment;


MI_uni_efficient=mi_performance(model_prop_mat_uni, model_prop_var_mat_uni, true_prop, rowobs);
MI_uni_efficient;

MI_biv_moment=mi_performance(obs_prop_mat_biv, obs_prop_var_mat_biv, true_prop, rowobs);
MI_biv_moment;


MI_biv_efficient=mi_performance(model_prop_mat_biv, model_prop_var_mat_biv, true_prop, rowobs);
MI_biv_efficient;



