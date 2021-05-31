# Example 6.8

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


# read the WinBUGS code;
model.file <- system.file(package="R2WinBUGS", "model", "logitmodel_missinginteractionforR.txt")
# Let's take a look:
file.show(model.file)

# logistic regression coefficients;
beta0=-1.5;
beta1=1/3;
beta2=1/6;

# mean and variance of x;
mu=2;
sigma=1;

# complete-data sample size;
rowobs=2000;

# number of simulations;
cycle_no=1000;

# number of imputations;
mi_no=50;

# matrices holding the logistic regression coefficient estimates;

BD_slope_vec=BD_slope_var_vec=rep(NA, cycle_no);
CC_slope_vec=CC_slope_var_vec=rep(NA, cycle_no);
bugs_slope_mat=bugs_slope_var_mat=pas_slope_mat=pas_slope_var_mat=jav_slope_mat=jav_slope_var_mat=matrix(NA, cycle_no, mi_no);


# begin the simulation;

for (cycle in 1:cycle_no)
{
set.seed(cycle);

x=mu+sqrt(sigma)*rnorm(rowobs);

y=rbinom(rowobs,size=1,prob=exp(beta0+beta1*x+beta2*x^2)/(1+exp(beta0+beta1*x+beta2*x^2)));


# apply the missing cases to complete data;
# missing at random;

alpha0=2;

alpha1=-2;

missing_indi=rbinom(rowobs,size=1,prob=1/(1+exp(alpha0+alpha1*y)));


y_miss=y;
x_miss=x;

# y_miss[miss_indi==1]=NA;
# miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];

x_miss[missing_indi==1]=NA;

# before deletion analysis;
# logistic regression coefficient for the slope of x^2;
BD_logistic=summary(glm(y~x+I(x^2), family=binomial(link="logit")));
BD_slope_vec[cycle]=BD_logistic$coeff[3,1];
BD_slope_var_vec[cycle]=BD_logistic$coeff[3,2]^2;

# logistic regression coefficient for the slope of x^2;
CC_logistic=summary(glm(y_miss~x_miss+I(x_miss^2), family=binomial(link="logit")));
CC_slope_vec[cycle]=CC_logistic$coeff[3,1];
CC_slope_var_vec[cycle]=CC_logistic$coeff[3,2]^2;


# winbugs imputation;

# apply missing cases to the data;
M=rowobs;

process_stage2_data_simulate=list("M","y_miss","x_miss");
process_stage2_init= function(){list(beta0=-1.5, beta1=1/3, beta2=1/6, mu=2, tau=1)};
# process_stage2_parameters=c("mu", "beta0", "beta1", "beta2", "mu", "alpha0", "alpha1", "tau", "psi", "y", "x1", "x2");
process_stage2_parameters=c("x_impute");


gibbs_no=10000;

process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file,
    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=gibbs_no/(2*mi_no), bugs.seed=cycle,
    bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);

#process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file="C:/Users/WDQ7/reliability_test_forR.txt",
#    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=1, bugs.seed=197789,
#    bugs.directory="C:/Program Files (x86)/OpenBUGS/OpenBUGS321", summary.only=FALSE, debug=FALSE,
#    working.directory=NULL, clearWD=TRUE);


attach.bugs(process_stage2.sim);
rm(process_stage2.sim);

y_completed_bugs=y;
x_completed_bugs=x_miss;

for (i in 1:mi_no)
{
x_completed_bugs[missing_indi==1]=x_impute[i,];

# logistic regression coefficient for x-square slope;
bugs_logistic=summary(glm(y_completed_bugs~x_completed_bugs+I(x_completed_bugs^2), family=binomial(link="logit")));
bugs_slope_mat[cycle,i]=bugs_logistic$coeff[3,1];
bugs_slope_var_mat[cycle,i]=bugs_logistic$coeff[3,2]^2;

}

# passive imputation;
y_completed_passive=y;
x_completed_passive=x_miss;
for (i in 1:mi_no)
{
x_imputed_passive=mice.impute.norm(x_miss, ry=as.logical(1-missing_indi), seed=i, x=y);
x_completed_passive[missing_indi==1]=x_imputed_passive;

# logistic regression coefficient for x-square slope;
pas_logistic=summary(glm(y_completed_passive~x_completed_passive+I(x_completed_passive^2), family=binomial(link="logit")));
pas_slope_mat[cycle,i]=pas_logistic$coeff[3,1];
pas_slope_var_mat[cycle,i]=pas_logistic$coeff[3,2]^2;
}

# just another variable imputation
y_miss_recode=y_miss+1;
x2_miss=x_miss^2;
ori_data=cbind(y_miss_recode, x_miss, x2_miss);
s=prelim.mix(ori_data,1);
thetahat=em.mix(s);

y_completed_jav=y;
for (i in 1:mi_no)
{
rngseed(i); # set random number generator seed
newtheta=da.mix(s, thetahat, steps=200, showits=TRUE) # take 200 steps
# getparam.mix(s, newtheta, corr=TRUE);
imp_data=imp.mix(s, newtheta);
x_completed_jav=imp_data[,2];
x2_completed_jav=imp_data[,3];


# multiple imputation analysis;

# logistic regression coefficient for x1 slope;
jav_logistic=summary(glm(y_completed_jav~x_completed_jav+x2_completed_jav, family=binomial(link="logit")));
jav_slope_mat[cycle,i]=jav_logistic$coeff[3,1];
jav_slope_var_mat[cycle,i]=jav_logistic$coeff[3,2]^2;


}

cat("the cycle is", cycle, "\n");

}

# the end of mulitple imputation;

##########################################################################
# Simulation results
# slope estimand;
# Table 6.6
# complete-data analysis;

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

# passive imputation;
pas_MI=mi_performance(pas_slope_mat, pas_slope_var_mat, true_slope, rowobs);
pas_MI;

# Just another variable imputation;
jav_MI=mi_performance(jav_slope_mat, jav_slope_var_mat, true_slope, rowobs);
jav_MI;


