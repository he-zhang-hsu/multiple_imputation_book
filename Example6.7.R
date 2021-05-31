# Example 6.7

rm(list=ls());
library(R2WinBUGS);
library(MASS);
library(mice);
library(norm);
library(mix);
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


# model.file is the WinBUGS syntax;
model.file <- system.file(package="R2WinBUGS", "model", "logitmodel_missing_bivariate_forR.txt")

# Let's take a look on the WinBUGS syntax
file.show(model.file)

# generate the logistic outcome;
beta0=-1;
beta1=0.5;
beta2=0.5;

# generate bivariate normal x variables;
mu=rep(0,2);
Sigma=matrix(c(1,0.5,0.5,1), nrow=2,ncol=2);

# Complete-data sample size;

rowobs=1000;

# number of simulations;
cycle_no=1000;
# cycle_no=1;

# number of multiple imputations;
mi_no=50;

# Vectors and matrices holding the parameter estimates;

BD_mean_vec=BD_mean_var_vec=CC_mean_vec=CC_mean_var_vec=rep(NA, cycle_no);

# multiple imputation;
# For the marginal mean estimates;
# bugs mean the joint modeling imputation using bugs; the JM method;
bugs_mean_mat=bugs_mean_var_mat=matrix(NA, cycle_no, mi_no);

# glm is the general location modeling imputation; the GLOM method;
glm_mean_mat=glm_mean_var_mat=matrix(NA, cycle_no, mi_no);

# For the logistic regression slope coefficient;

BD_slope_vec=BD_slope_var_vec=rep(NA, cycle_no);
CC_slope_vec=CC_slope_var_vec=rep(NA, cycle_no);
bugs_slope_mat=bugs_slope_var_mat=glm_slope_mat=glm_slope_var_mat=matrix(NA, cycle_no, mi_no);


# begin the simulation;

for (cycle in 1:cycle_no)

{

set.seed(cycle);

x1_x2=mvrnorm(n = rowobs, mu=mu, Sigma=Sigma);

x1=x1_x2[,1];
x2=x1_x2[,2];

# generate complete data
# logistic outcome;
y=rbinom(rowobs,size=1,prob=exp(beta0+beta1*x1+beta2*x2)/(1+exp(beta0+beta1*x1+beta2*x2)));

# apply the missing cases to complete data;
# missing completely at random;

missing_indi=cbind(rbinom(rowobs,size=1,prob=0.15), rbinom(rowobs,size=1, prob=0.15), rbinom(rowobs, size=1, prob=0.15));

y_miss=y;
x1_miss=x1;
x2_miss=x2;

# y_miss[miss_indi==1]=NA;
# miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];


y_miss[missing_indi[,1]==1]=NA;
x1_miss[missing_indi[,2]==1]=NA;
x2_miss[missing_indi[,3]==1]=NA;

y_obs_no=length(y_miss[!is.na(y_miss)]);
y_mis_no=rowobs-y_obs_no;
y_obs=y[missing_indi[,1]==0];

# Before-deletion analysis;

# sample mean;

BD_mean_vec[cycle]=mean(y);
BD_mean_var_vec[cycle]=var(y)/rowobs;

# logistic regression coefficient for the slope of x1;
BD_logistic=summary(glm(y~x1+x2, family=binomial(link="logit")));
BD_slope_vec[cycle]=BD_logistic$coeff[2,1];
BD_slope_var_vec[cycle]=BD_logistic$coeff[2,2]^2;

# complete-case analysis;
# mean estimands;
CC_mean=mean(y_miss, na.rm="T");
CC_mean_vec[cycle]=CC_mean;
CC_mean_var_vec[cycle]=var(y_miss, na.rm="T")/y_obs_no;

# logistic regression coefficient for the slope of x1;
CC_logistic=summary(glm(y_miss~x1_miss+x2_miss, family=binomial(link="logit")));
CC_slope_vec[cycle]=CC_logistic$coeff[2,1];
CC_slope_var_vec[cycle]=CC_logistic$coeff[2,2]^2;


# winbugs imputation;

# apply missing cases to the data;
M=rowobs;

process_stage2_data_simulate=list(M=M, y_miss=y_miss, x1_miss=x1_miss, x2_miss=x2_miss);
process_stage2_init= function(){list(beta0=-1, beta1=0.5, beta2=0.5, mu=0, alpha0=0, alpha1=0, tau=1, psi=1)};


# If we want to draw the posterior samples then output the model parameters;
# process_stage2_parameters=c("mu", "beta0", "beta1", "beta2", "mu", "alpha0", "alpha1", "tau", "psi", "y", "x1", "x2");

# Otherwise we only need to output the imputed values for the variables;
process_stage2_parameters=c("y_impute", "x1_impute", "x2_impute");


gibbs_no=10000;

# imputation Gibbs;

process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file,
    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=gibbs_no/(2*mi_no), bugs.seed=cycle,
    bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);

attach.bugs(process_stage2.sim);
rm(process_stage2.sim);

y_completed_bugs=y_miss;
x1_completed_bugs=x1_miss;
x2_completed_bugs=x2_miss;

for (i in 1:mi_no)
{
y_completed_bugs[missing_indi[,1]==1]=y_impute[i,];
x1_completed_bugs[missing_indi[,2]==1]=x1_impute[i,];
x2_completed_bugs[missing_indi[,3]==1]=x2_impute[i,];


# marginal means;
bugs_mean=mean(y_completed_bugs);
bugs_var=var(y_completed_bugs)/rowobs;

bugs_mean_mat[cycle,i]=bugs_mean;
bugs_mean_var_mat[cycle,i]=bugs_var;

# logistic regression coefficient for x1 slope;
bugs_logistic=summary(glm(y_completed_bugs~x1_completed_bugs+x2_completed_bugs, family=binomial(link="logit")));
bugs_slope_mat[cycle,i]=bugs_logistic$coeff[2,1];
bugs_slope_var_mat[cycle,i]=bugs_logistic$coeff[2,2]^2;

}

# general location modeling imputation;
# make 0/1 variable to 1/2 factor;
y_miss_recode=y_miss+1;
ori_data=cbind(y_miss_recode, x1_miss, x2_miss);
s=prelim.mix(ori_data,1);
thetahat=em.mix(s);

for (i in 1:mi_no)
{
rngseed(i); # set random number generator seed
newtheta=da.mix(s, thetahat, steps=200, showits=TRUE) # take 200 steps
# getparam.mix(s, newtheta, corr=TRUE);
imp_data=imp.mix(s, newtheta);

y_completed_glm=imp_data[,1]-1;
x1_completed_glm=imp_data[,2];
x2_completed_glm=imp_data[,3];


# multiple imputation analysis;

# marginal means;
glm_mean=mean(y_completed_glm);
glm_var=var(y_completed_glm)/rowobs;

glm_mean_mat[cycle,i]=glm_mean;
glm_mean_var_mat[cycle,i]=glm_var;

# logistic regression coefficient for x1 slope;
glm_logistic=summary(glm(y_completed_glm~x1_completed_glm+x2_completed_glm, family=binomial(link="logit")));
glm_slope_mat[cycle,i]=glm_logistic$coeff[2,1];
glm_slope_var_mat[cycle,i]=glm_logistic$coeff[2,2]^2;



}

cat("the cycle is", cycle, "\n");

}

# the end of mulitple imputation;
#############################################################################
# Simulation results;
# Table 6.5 
# marginal mean estimand;

# marginal mean of Y;
true_mean=mean(BD_mean_vec);
true_mean;

# before-deletion;
BD=performance(BD_mean_vec, BD_mean_var_vec, true_mean);
BD;

# complete-case analysis;
CC=performance(CC_mean_vec, CC_mean_var_vec, true_mean);
CC;

# Multiple imputation analysis;
# JM/bugs imputation;
bugs_MI=mi_performance(bugs_mean_mat, bugs_mean_var_mat, true_mean, rowobs);
bugs_MI;

# general location modeling imputation;
glm_MI=mi_performance(glm_mean_mat, glm_mean_var_mat, true_mean, rowobs);
glm_MI;

##########################################################################
# slope estimand;

true_slope=mean(BD_slope_vec);
true_slope;

# before-deletion;
BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;

# complete-case analysis;
CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;

# JM/bugs imputation;
bugs_MI=mi_performance(bugs_slope_mat, bugs_slope_var_mat, true_slope, rowobs);
bugs_MI;

# general location modeling imputation;
glm_MI=mi_performance(glm_slope_mat, glm_slope_var_mat, true_slope, rowobs);
glm_MI;

##########################################################################
# Make trace plot of posterior draws from the last simulation replicate
# For Fig. 6.6
process_stage2_parameters=c("beta0", "beta1", "beta2", "mu", "alpha0", "alpha1", "tau", "psi");

process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file,
    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=1, bugs.seed=cycle,
    bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);

attach.bugs(process_stage2.sim);
rm(process_stage2.sim);


plot(seq(gibbs_no/2+1, gibbs_no, 1), beta0, type="l", xlab="iteration", main="", ylab="beta0");
plot(seq(gibbs_no/2+1, gibbs_no, 1), beta1, type="l", xlab="iteration", main="", ylab="beta1");
plot(seq(gibbs_no/2+1, gibbs_no, 1), beta2, type="l", xlab="iteration", main="", ylab="beta2");
plot(seq(gibbs_no/2+1, gibbs_no, 1), mu, type="l", xlab="iteration", main="", ylab="mu");
plot(seq(gibbs_no/2+1, gibbs_no, 1), alpha0, type="l", xlab="iteration", main="", ylab="alpha0");
plot(seq(gibbs_no/2+1, gibbs_no, 1), alpha1, type="l", xlab="iteration", main="", ylab="alpha1");
plot(seq(gibbs_no/2+1, gibbs_no, 1), tau, type="l", xlab="iteration", main="", ylab="tau");
plot(seq(gibbs_no/2+1, gibbs_no, 1), psi, type="l", xlab="iteration", main="", ylab="psi");

