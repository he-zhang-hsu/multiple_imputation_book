# Example 5.6
rm(list=ls());


# library(car);
library(mice);
library(norm);
# adaptive metroplis rejection sampling;
library(HI);
# include the library MASS for generating from multivariate normal distribution;
library(MASS);
# include the library msm and bayesm for using the function of drawing truncated normal
# distribution;
library(msm);
library(bayesm);
library(mvtnorm);
library(sn);
library(arm);
# library(MCMCpack);
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


# set up the random seed;
set.seed(197789);

# complete-data sample size;
rowobs=1000;

# mean and variance of x;
mu_x=1;
var_x=1/4;

# regression parameters;
# both beta0 and beta1 are fixed at 
beta0=-2;
beta1=1; 

# Scenario I
# beta2=0;
# Scenario II
beta2=1/2;

# number of simulations;
cycle_no=1000;

# number of imputations;
mi_no=30;

# draw_no is the MCMC iterations;
draw_no=1000;

# Vectors and matrices holding the parameter estimates;
# For the marginal mean
BD_mean_vec=BD_mean_var_vec=CC_mean_vec=CC_mean_var_vec=rep(NA, cycle_no);

# For multiple imputation
# LOGIT: logistic regression imputation using only the main effect of x;
# ROUND: linear imputation using only the main effect of x and round;
# PMM_linear: PMM imputation at the linear scale;
# PMM_logit: PMM imputation at the logit scale;
LOGIT_mean_mat=LOGIT_mean_var_mat=ROUND_mean_mat=ROUND_mean_var_mat=matrix(NA, cycle_no, mi_no);

PMM_linear_mean_mat=PMM_linear_mean_var_mat=PMM_logit_mean_mat=PMM_logit_mean_var_mat=matrix(NA, cycle_no, mi_no);

# For the logistic regression coefficient
BD_slope_vec=BD_slope_var_vec=rep(NA, cycle_no);
CC_slope_vec=CC_slope_var_vec=rep(NA, cycle_no);

LOGIT_slope_mat=LOGIT_slope_var_mat=ROUND_slope_mat=ROUND_slope_var_mat=matrix(NA, cycle_no, mi_no);
PMM_linear_slope_mat=PMM_linear_slope_var_mat=PMM_logit_slope_mat=PMM_logit_slope_var_mat=matrix(NA, cycle_no, mi_no);


# begin the simulation;

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# generate the logistic regression outcome y;
# Model I
# y=rbinom(n=rowobs, size=1, prob=exp(beta0+beta1*x)/(1+exp(beta0+beta1*x)));

# Model II
y=rbinom(n=rowobs, size=1, prob=exp(beta0+beta1*x+beta2*x^2)/(1+exp(beta0+beta1*x+beta2*x^2)));

# set up the missing data as MCAR on x;

miss_indi=runif(n=rowobs)<0.30;
y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
x_obs=x[miss_indi==0];
x_miss=x[miss_indi==1];

# compelete data analysis;
# sample mean;
BD_mean_vec[cycle]=mean(y);
BD_mean_var_vec[cycle]=var(y)/rowobs;

# logistic regression coefficient for slope;
# under model 1
# BD_logistic=summary(glm(y~x, family=binomial(link="logit")));
# linear slope coefficient
# BD_slope_vec[cycle]=BD_logistic$coeff[2,1];
# BD_slope_var_vec[cycle]=BD_logistic$coeff[2,2]^2;


# under model 2
BD_logistic=summary(glm(y~x+I(x*x), family=binomial(link="logit")));
# linear term coefficient
BD_slope_vec[cycle]=BD_logistic$coeff[2,1];
BD_slope_var_vec[cycle]=BD_logistic$coeff[2,2]^2;



# different missing data methods;

# complete-case analysis;
# mean estimands;
CC_mean=mean(y_miss, na.rm="T");
CC_mean_vec[cycle]=CC_mean;
CC_mean_var_vec[cycle]=var(y_miss, na.rm="T")/obs_no;

# logistic regression coefficient for slope;
# logistic regression coefficient for slope;
# under model 1
# CC_logistic=summary(glm(y_obs~x_obs, family=binomial(link="logit")));
# linear slope coefficient;
# CC_slope_vec[cycle]=CC_logistic$coeff[2,1];
# CC_slope_var_vec[cycle]=CC_logistic$coeff[2,2]^2;


# under model 2
CC_logistic=summary(glm(y_obs~x_obs+I(x_obs*x_obs), family=binomial(link="logit")));
# linear term coefficient
CC_slope_vec[cycle]=CC_logistic$coeff[2,1];
CC_slope_var_vec[cycle]=CC_logistic$coeff[2,2]^2;


# now impute the missing y's to get estimates of the marginal;
y_completed_LOGIT=y_completed_ROUND=y_completed_PMM_linear=y_completed_PMM_logit=y_miss;

for (i in 1:mi_no)
{

set.seed(i);

# bivariate model imputation; 
# Logistic regression imputation only includeing x;
y_imputed_LOGIT=mice.impute.logreg(y_miss, ry=as.logical(1-miss_indi), seed=i, x=as.matrix(x));

# linear normal imputation and then rounding;
y_imputed_ROUND=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=as.matrix(x));
y_imputed_ROUND[y_imputed_ROUND >=0.5]=1;
y_imputed_ROUND[y_imputed_ROUND < 0.5]=0;

# pmm imputation based on a linear normal model;
y_imputed_PMM_linear=mice.impute.pmm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=as.matrix(x));

# pmm imputation based on a logit model;
y_imputed_PMM_logit=rep(NA, mis_no);


CC_logistic=summary(glm(y_obs~x_obs, family=binomial(link="logit")));
obs_pred_mean=cbind(cbind(1, x_obs)%*%CC_logistic$coefficients[,1],y_obs);

# draw logistic regression coefficients from a posterior;
coef_s1x_posterior=bayesglm(y_obs~x_obs, family=binomial(link="logit"), prior.scale=Inf, prior.df=Inf);  
coef_s1x_draw_collection=sim(coef_s1x_posterior, n.sims=draw_no)@coef;

coef_s1x_draw=coef_s1x_draw_collection[draw_no,];

miss_pred_draw=cbind(1, x_miss)%*%coef_s1x_draw;

# choose number of donors as 5;
donor_no=5;
for(k in 1:mis_no)
{
donor_matrix=obs_pred_mean;
# distance function;
donor_matrix[,1]=abs(miss_pred_draw[k]-obs_pred_mean[,1]);
donor_matrix_sorted=donor_matrix[order(donor_matrix[,1]),];
y_imputed_PMM_logit[k]=sample(donor_matrix_sorted[1:donor_no,2],1);
}

y_completed_LOGIT[miss_seq]=y_imputed_LOGIT;

y_completed_ROUND[miss_seq]=y_imputed_ROUND;

y_completed_PMM_linear[miss_seq]=y_imputed_PMM_linear;

y_completed_PMM_logit[miss_seq]=y_imputed_PMM_logit;

# Store the multiple imputation estimates;
# For means;
# logistic regression imputation;
LOGIT_mean_mat[cycle,i]=mean(y_completed_LOGIT);
LOGIT_mean_var_mat[cycle,i]=var(y_completed_LOGIT)/rowobs;

# logistic regression coefficient for slope;
# under model 1
# logit_logistic=summary(glm(y_completed_LOGIT~x, family=binomial(link="logit")));


# LOGIT_slope_mat[cycle,i]=logit_logistic$coeff[2,1];
# LOGIT_slope_var_mat[cycle,i]=logit_logistic$coeff[2,2]^2;


# under model 2
logit_logistic=summary(glm(y_completed_LOGIT~x+I(x*x), family=binomial(link="logit")));

LOGIT_slope_mat[cycle,i]=logit_logistic$coeff[2,1];
LOGIT_slope_var_mat[cycle,i]=logit_logistic$coeff[2,2]^2;



# rounding imputation;

ROUND_mean_mat[cycle,i]=mean(y_completed_ROUND);
ROUND_mean_var_mat[cycle,i]=var(y_completed_ROUND)/rowobs;

# logistic regression coefficient for slope;
# under model 1
# round_logistic=summary(glm(y_completed_ROUND~x, family=binomial(link="logit")));
# ROUND_slope_mat[cycle,i]=round_logistic$coeff[2,1];
# ROUND_slope_var_mat[cycle,i]=round_logistic$coeff[2,2]^2;


# under model 2
round_logistic=summary(glm(y_completed_ROUND~x+I(x*x), family=binomial(link="logit")));

ROUND_slope_mat[cycle,i]=round_logistic$coeff[2,1];
ROUND_slope_var_mat[cycle,i]=round_logistic$coeff[2,2]^2;

# PMM on linear model imputation;

PMM_linear_mean_mat[cycle,i]=mean(y_completed_PMM_linear);
PMM_linear_mean_var_mat[cycle,i]=var(y_completed_PMM_linear)/rowobs;

# logistic regression coefficient for slope;
# under model 1
# PMM_linear_logistic=summary(glm(y_completed_PMM_linear~x, family=binomial(link="logit")));
# PMM_linear_slope_mat[cycle,i]=PMM_linear_logistic$coeff[2,1];
# PMM_linear_slope_var_mat[cycle,i]=PMM_linear_logistic$coeff[2,2]^2;



# under model 2
PMM_linear_logistic=summary(glm(y_completed_PMM_linear~x+I(x*x), family=binomial(link="logit")));

PMM_linear_slope_mat[cycle,i]=PMM_linear_logistic$coeff[2,1];
PMM_linear_slope_var_mat[cycle,i]=PMM_linear_logistic$coeff[2,2]^2;

# pmm on logistic regression model imputation;

PMM_logit_mean_mat[cycle,i]=mean(y_completed_PMM_logit);
PMM_logit_mean_var_mat[cycle,i]=var(y_completed_PMM_logit)/rowobs;

# logistic regression coefficient for slope;
# under model 1
# PMM_logit_logistic=summary(glm(y_completed_PMM_logit~x, family=binomial(link="logit")));
# PMM_logit_slope_mat[cycle,i]=PMM_logit_logistic$coeff[2,1];
# PMM_logit_slope_var_mat[cycle,i]=PMM_logit_logistic$coeff[2,2]^2;


# under model 2
PMM_logit_logistic=summary(glm(y_completed_PMM_logit~x+I(x*x), family=binomial(link="logit")));
PMM_logit_slope_mat[cycle,i]=PMM_logit_logistic$coeff[2,1];
PMM_logit_slope_var_mat[cycle,i]=PMM_logit_logistic$coeff[2,2]^2;



}

cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

#############################################################################
# Simulation results
# Table 5.4
# For the marginal mean estimand

true_mean=mean(BD_mean_vec);
true_mean;

# before-deletion;
BD=performance(BD_mean_vec, BD_mean_var_vec, true_mean);
BD;

# complete-case analysis;
CC=performance(CC_mean_vec, CC_mean_var_vec, true_mean);
CC;

# Multiple imputation analysis;
# Logit imputation 

LOGIT_MI=mi_performance(LOGIT_mean_mat, LOGIT_mean_var_mat, true_mean, rowobs);
LOGIT_MI;

# Imputation based on rounding
ROUND_MI=mi_performance(ROUND_mean_mat, ROUND_mean_var_mat, true_mean, rowobs);
ROUND_MI;

# PMM Imputation based on a linear model;
PMM_linear_MI=mi_performance(PMM_linear_mean_mat, PMM_linear_mean_var_mat, true_mean, rowobs);
PMM_linear_MI;


# PMM Imputation based on a logistic model;
PMM_logit_MI=mi_performance(PMM_logit_mean_mat, PMM_logit_mean_var_mat, true_mean, rowobs);
PMM_logit_MI;


##########################################################################
# slope estimand;

# complete-data analysis;

true_slope=mean(BD_slope_vec);
true_slope;

# before-deletion;
BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;

# complete-case analysis;
CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;

# Multiple imputation analysis;
# Logit imputation;

LOGIT_MI=mi_performance(LOGIT_slope_mat, LOGIT_slope_var_mat, true_slope, rowobs);
LOGIT_MI;

# Round imputation;
ROUND_MI=mi_performance(ROUND_slope_mat, ROUND_slope_var_mat, true_slope, rowobs);
ROUND_MI;

# PMM imputation based on a linear normal model;
PMM_linear_MI=mi_performance(PMM_linear_slope_mat, PMM_linear_slope_var_mat, true_slope, rowobs);
PMM_linear_MI;


# PMM imputation based on a logit model;
PMM_logit_MI=mi_performance(PMM_logit_slope_mat, PMM_logit_slope_var_mat, true_slope, rowobs);
PMM_logit_MI;


#################################################################
# Some exploratory plots
# observed data;
# logistic imputation;
plot(ecdf(x_obs[y_obs==1]), col="red");
lines(ecdf(x_obs[y_obs==0]), col="green");

lines(ecdf(x_miss[y_imputed_LOGIT==1]), col="blue");
lines(ecdf(x_miss[y_imputed_LOGIT==0]), col="yellow");

# rounding imputation;

plot(ecdf(x_obs[y_obs==1]), col="red");
lines(ecdf(x_obs[y_obs==0]), col="green");

lines(ecdf(x_miss[y_imputed_ROUND==1]), col="blue");
lines(ecdf(x_miss[y_imputed_ROUND==0]), col="yellow");

# PMM based on the linear normal model;

plot(ecdf(x_obs[y_obs==1]), col="red");
lines(ecdf(x_obs[y_obs==0]), col="green");

lines(ecdf(x_miss[y_imputed_PMM_linear==1]), col="blue");
lines(ecdf(x_miss[y_imputed_PMM_linear==0]), col="yellow");


# PMM based on the logistic model;

plot(ecdf(x_obs[y_obs==1]), col="red");
lines(ecdf(x_obs[y_obs==0]), col="green");

lines(ecdf(x_miss[y_imputed_PMM_logit==1]), col="blue");
lines(ecdf(x_miss[y_imputed_PMM_logit==0]), col="yellow");


