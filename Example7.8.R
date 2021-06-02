# Example 7.8

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


# WinBUGS code for JM
model0.file <- system.file(package="R2WinBUGS", "model", "logitmodel_missinginteraction_JMforR.txt")
# Let's take a look:
file.show(model0.file)

# WinBUGS code for FCS (a)
model1.file <- system.file(package="R2WinBUGS", "model", "logitmodel_missinginteraction_fcsforR.txt")
# Let's take a look:
file.show(model1.file)

# WinBUGS code for FCS (b)
model2.file <- system.file(package="R2WinBUGS", "model", "logitmodel_missinginteraction_fcs_mix_forR.txt")
# Let's take a look:
file.show(model2.file)


# regression coefficients;

beta0=-1.5;
beta1=1/3;
beta2=1/6;

# mean and variance of x;
mu=2;
sigma=1;

# complete-data sample size;

rowobs=1000;

# number of simulations;

cycle_no=1000;

# number of imputations;
mi_no=50;

# Vectors and matrices holding the parameter estimates;

BD_slope_vec=BD_slope_var_vec=rep(NA, cycle_no);
CC_slope_vec=CC_slope_var_vec=rep(NA, cycle_no);

JM_slope_mat=JM_slope_var_mat=FCS_a_slope_mat=FCS_a_slope_var_mat=FCS_b_slope_mat=FCS_b_slope_var_mat=matrix(NA, cycle_no, mi_no);


# begin the simulation;

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu+sqrt(sigma)*rnorm(rowobs);

y=rbinom(rowobs,size=1,prob=exp(beta0+beta1*x+beta2*x^2)/(1+exp(beta0+beta1*x+beta2*x^2)));


# apply the missing cases to complete data;

missing_indi=cbind(rbinom(rowobs,size=1,prob=0.20), rbinom(rowobs,size=1, prob=0.20));

missing_indi_x=missing_indi[,1];
missing_indi_y=missing_indi[,2];

y_miss=y;
x_miss=x;

x_miss[missing_indi_x==1]=NA;
y_miss[missing_indi_y==1]=NA;

# before deletion analysis;

# logistic regression coefficient for the slope of x^2;
BD_logistic=summary(glm(y~x+I(x^2), family=binomial(link="logit")));
BD_slope_vec[cycle]=BD_logistic$coeff[3,1];
BD_slope_var_vec[cycle]=BD_logistic$coeff[3,2]^2;

# logistic regression coefficient for the slope of x^2;
CC_logistic=summary(glm(y_miss~x_miss+I(x_miss^2), family=binomial(link="logit")));
CC_slope_vec[cycle]=CC_logistic$coeff[3,1];
CC_slope_var_vec[cycle]=CC_logistic$coeff[3,2]^2;


# JM winbugs imputation;


M=rowobs;

gibbs_no=10000;

process_stage2_data_simulate_JM=list("M","y_miss", "x_miss");
process_stage2_init_JM= function(){list(beta0=-1.5, beta1=1/3, beta2=1/6, mu=2, tau=1)};
process_stage2_parameters_JM=c("x_impute", "y_impute");

process_stage2.sim.JM <- bugs(data=process_stage2_data_simulate_JM, inits=process_stage2_init_JM, parameters=process_stage2_parameters_JM, model.file=model0.file,
    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=gibbs_no/(2*mi_no), bugs.seed=cycle,
    bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);


attach.bugs(process_stage2.sim.JM);
rm(process_stage2.sim.JM);

y_completed_bugs_JM=y_miss;
x_completed_bugs_JM=x_miss;

for (i in 1:mi_no)
{
x_completed_bugs_JM[missing_indi_x==1]=x_impute[i,];
y_completed_bugs_JM[missing_indi_y==1]=y_impute[i,];

# logistic regression coefficient for x-square slope;
JM_logistic=summary(glm(y_completed_bugs_JM~x_completed_bugs_JM+I(x_completed_bugs_JM^2), family=binomial(link="logit")));
JM_slope_mat[cycle,i]=JM_logistic$coeff[3,1];
JM_slope_var_mat[cycle,i]=JM_logistic$coeff[3,2]^2;

}

# FCS (a) imputation assuming equal variance;

process_stage2_data_simulate_FCS=list("M","y_miss", "x_miss");
process_stage2_init_FCS= function(){list(beta0=-1.5, beta1=1/3, beta2=1/6, alpha0=2, alpha1=0, tau=1)};
process_stage2_parameters_FCS=c("x_impute", "y_impute");

process_stage2.sim.FCS <- bugs(data=process_stage2_data_simulate_FCS, inits=process_stage2_init_FCS, parameters=process_stage2_parameters_FCS, model.file=model1.file,
    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=gibbs_no/(2*mi_no), bugs.seed=cycle,
    bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);


attach.bugs(process_stage2.sim.FCS);
rm(process_stage2.sim.FCS);

y_completed_FCS_a=y_miss;
x_completed_FCS_a=x_miss;

for (i in 1:mi_no)
{
x_completed_FCS_a[missing_indi_x==1]=x_impute[i,];
y_completed_FCS_a[missing_indi_y==1]=y_impute[i,];


# logistic regression coefficient for x-square slope;
FCS_a_logistic=summary(glm(y_completed_FCS_a~x_completed_FCS_a+I(x_completed_FCS_a^2), family=binomial(link="logit")));
FCS_a_slope_mat[cycle,i]=FCS_a_logistic$coeff[3,1];
FCS_a_slope_var_mat[cycle,i]=FCS_a_logistic$coeff[3,2]^2;

}

# FCS (b) imputation assuming unequal variances;

process_stage2_data_simulate_FCS_mix=list("M","y_miss", "x_miss");
process_stage2_init_FCS_mix= function(){list(beta0=-1.5, beta1=1/3, beta2=1/6, mu=c(2,2), tau=c(1,1))};
process_stage2_parameters_FCS_mix=c("x_impute", "y_impute");

process_stage2.sim.FCS.mix <- bugs(data=process_stage2_data_simulate_FCS_mix, inits=process_stage2_init_FCS_mix, parameters=process_stage2_parameters_FCS_mix, model.file=model2.file,
    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=gibbs_no/(2*mi_no), bugs.seed=cycle,
    bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, 
    working.directory=NULL, clearWD=TRUE, debug=FALSE);


attach.bugs(process_stage2.sim.FCS.mix);
rm(process_stage2.sim.FCS.mix);

y_completed_FCS_b=y_miss;
x_completed_FCS_b=x_miss;
for (i in 1:mi_no)
{
x_completed_FCS_b[missing_indi_x==1]=x_impute[i,];
y_completed_FCS_b[missing_indi_y==1]=y_impute[i,];

# multiple imputation analysis;

# logistic regression coefficient for x1 slope;
FCS_b_logistic=summary(glm(y_completed_FCS_b~x_completed_FCS_b+I(x_completed_FCS_b^2), family=binomial(link="logit")));
FCS_b_slope_mat[cycle,i]=FCS_b_logistic$coeff[3,1];
FCS_b_slope_var_mat[cycle,i]=FCS_b_logistic$coeff[3,2]^2;

}

cat("the cycle is", cycle, "\n");

}

# the end of mulitple imputation;

##########################################################################
# Simulation results
# Table 7.5
# slope estimand;

true_slope=mean(BD_slope_vec);

# before-deletion;
BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;

# complete-case analysis;
CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;

# JM imputation;
JM_MI=mi_performance(JM_slope_mat, JM_slope_var_mat, true_slope, rowobs);
JM_MI;

# FCS (a) imputation;
FCS_a_MI=mi_performance(FCS_a_slope_mat, FCS_a_slope_var_mat, true_slope, rowobs);
FCS_a_MI;

# FCS (b) imputation;
FCS_b_MI=mi_performance(FCS_b_slope_mat, FCS_b_slope_var_mat, true_slope, rowobs);
FCS_b_MI;

# generate the logistic outcome;
# beta0=-1.5;
# beta1=1/3;
# beta2=1/6;

# generate the univariate normal x variable;
# mu=2;
# sigma=1;

# rowobs=1000000;
# x=mu+sqrt(sigma)*rnorm(rowobs);

# y=rbinom(rowobs,size=1,prob=exp(beta0+beta1*x+beta2*x^2)/(1+exp(beta0+beta1*x+beta2*x^2)));


