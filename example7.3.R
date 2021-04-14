library(R2WinBUGS);
library(MASS);
library(mice);
library(norm);
library(mix);
# library(HI);
options(digits=4);
rm(list=ls());


# setwd('\\\\cdc.gov\\private\\L728\\book\\example8.4')

# generate the logistic outcome;
beta0=-1.5;
beta1=1/3;
beta2=1/6;

# generate the univariate normal x variable;
mu=2;
sigma=1;


rowobs=2000;
cycle_no=1000;
mi_no=50;

# matrices holding the parameter estimates;

BD_slope_vec=BD_slope_var_vec=rep(NA, cycle_no);
CC_slope_vec=CC_slope_var_vec=rep(NA, cycle_no);
bugs_slope_mat=bugs_slope_var_mat=pas_slope_mat=pas_slope_var_mat=jav_slope_mat=jav_slope_var_mat=matrix(NA, cycle_no, mi_no);


# generate data for example 8.5
for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu+sqrt(sigma)*rnorm(rowobs);

y2=rbinom(rowobs,size=1,prob=exp(beta0+beta1*x+beta2*x^2)/(1+exp(beta0+beta1*x+beta2*x^2)));


# apply the missing cases to complete data;
# missing at random;

alpha0=2;

alpha1=-2;

missing_indi=rbinom(rowobs,size=1,prob=1/(1+exp(alpha0+alpha1*y2)));


y2_miss=y2;
x_miss=x;

# y2_miss[miss_indi==1]=NA;
# miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];

x_miss[missing_indi==1]=NA;

# before deletion analysis;
# compelete data analysis;

# logistic regression coefficient for the slope of x^2;
BD_logistic=summary(glm(y2~x+I(x^2), family=binomial(link="logit")));
BD_slope_vec[cycle]=BD_logistic$coeff[3,1];
BD_slope_var_vec[cycle]=BD_logistic$coeff[3,2]^2;

# logistic regression coefficient for the slope of x^2;
CC_logistic=summary(glm(y2_miss~x_miss+I(x_miss^2), family=binomial(link="logit")));
CC_slope_vec[cycle]=CC_logistic$coeff[3,1];
CC_slope_var_vec[cycle]=CC_logistic$coeff[3,2]^2;


# winbugs imputation;

# apply missing cases to the data;
M=rowobs;

process_stage2_data_simulate=list("M","y2", "x_miss");
process_stage2_init= function(){list(beta0=-1.5, beta1=1/3, beta2=1/6, mu=2, tau=1)};
# process_stage2_parameters=c("mu", "beta0", "beta1", "beta2", "mu", "alpha0", "alpha1", "tau", "psi", "y2", "x1", "x2");
process_stage2_parameters=c("x_impute");


gibbs_no=10000;


process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file="C:/Users/WDQ7/book_program/logitmodel_missinginteractionforR.txt",
    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=gibbs_no/(2*mi_no), bugs.seed=cycle,
    bugs.directory="C:/Users/WDQ7/Desktop/WinBUGS14", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);





#process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file="C:/Users/WDQ7/reliability_test_forR.txt",
#    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=1, bugs.seed=197789,
#    bugs.directory="C:/Program Files (x86)/OpenBUGS/OpenBUGS321", summary.only=FALSE, debug=FALSE,
#    working.directory=NULL, clearWD=TRUE);


attach.bugs(process_stage2.sim);
rm(process_stage2.sim);

y2_completed_bugs=y2;
x_completed_bugs=x_miss;

for (i in 1:mi_no)
{
x_completed_bugs[missing_indi==1]=x_impute[i,];

# logistic regression coefficient for x-square slope;
bugs_logistic=summary(glm(y2_completed_bugs~x_completed_bugs+I(x_completed_bugs^2), family=binomial(link="logit")));
bugs_slope_mat[cycle,i]=bugs_logistic$coeff[3,1];
bugs_slope_var_mat[cycle,i]=bugs_logistic$coeff[3,2]^2;

}

# passive imputation;
y2_completed_passive=y2;
x_completed_passive=x_miss;
for (i in 1:mi_no)
{
x_imputed_passive=mice.impute.norm(x_miss, ry=as.logical(1-missing_indi), seed=i, x=y2);
x_completed_passive[missing_indi==1]=x_imputed_passive;

# logistic regression coefficient for x-square slope;
pas_logistic=summary(glm(y2_completed_passive~x_completed_passive+I(x_completed_passive^2), family=binomial(link="logit")));
pas_slope_mat[cycle,i]=pas_logistic$coeff[3,1];
pas_slope_var_mat[cycle,i]=pas_logistic$coeff[3,2]^2;
}

# just another variable imputation
y2_miss_recode=y2_miss+1;
x2_miss=x_miss^2;
ori_data=cbind(y2_miss_recode, x_miss, x2_miss);
s=prelim.mix(ori_data,1);
thetahat=em.mix(s);

y2_completed_jav=y2;
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
jav_logistic=summary(glm(y2_completed_jav~x_completed_jav+x2_completed_jav, family=binomial(link="logit")));
jav_slope_mat[cycle,i]=jav_logistic$coeff[3,1];
jav_slope_var_mat[cycle,i]=jav_logistic$coeff[3,2]^2;



}

cat("the cycle is", cycle, "\n");

}

# the end of mulitple imputation;

##########################################################################
# slope estimand;

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

bugs_mi_mean=rep(NA, cycle_no);
bugs_mi_var=rep(NA, cycle_no);
bugs_mi_df=rep(NA, cycle_no);
bugs_mi_f=rep(NA, cycle_no);

pas_mi_mean=rep(NA, cycle_no);
pas_mi_var=rep(NA, cycle_no);
pas_mi_df=rep(NA, cycle_no);
pas_mi_f=rep(NA, cycle_no);

jav_mi_mean=rep(NA, cycle_no);
jav_mi_var=rep(NA, cycle_no);
jav_mi_df=rep(NA, cycle_no);
jav_mi_f=rep(NA, cycle_no);



for (i in 1:cycle_no)
{
bugs_summary=pool.scalar(bugs_slope_mat[i,], bugs_slope_var_mat[i,], n=rowobs);
bugs_mi_mean[i]=bugs_summary$qbar;
bugs_mi_var[i]=bugs_summary$t;
bugs_mi_df[i]=bugs_summary$df;
bugs_mi_f[i]=bugs_summary$f;

pas_summary=pool.scalar(pas_slope_mat[i,], pas_slope_var_mat[i,], n=rowobs);
pas_mi_mean[i]=pas_summary$qbar;
pas_mi_var[i]=pas_summary$t;
pas_mi_df[i]=pas_summary$df;
pas_mi_f[i]=pas_summary$f;


jav_summary=pool.scalar(jav_slope_mat[i,], jav_slope_var_mat[i,], n=rowobs);
jav_mi_mean[i]=jav_summary$qbar;
jav_mi_var[i]=jav_summary$t;
jav_mi_df[i]=jav_summary$df;
jav_mi_f[i]=jav_summary$f;


}



# For slope;

# proper imputation;

bugs_mi_bias=mean(bugs_mi_mean)-true_slope;
bugs_mi_bias;

bugs_mi_bias/true_slope;

bugs_mi_mse=mean((bugs_mi_mean-true_slope)^2);
bugs_mi_mse;

# coverage;
bugs_slope_low95_vec=bugs_mi_mean-qt(.975, bugs_mi_df)*sqrt(bugs_mi_var);
bugs_slope_up95_vec=bugs_mi_mean+qt(.975, bugs_mi_df)*sqrt(bugs_mi_var);

bugs_slope_length=mean(bugs_slope_up95_vec-bugs_slope_low95_vec);
bugs_slope_length;


bugs_slope_coverage=(bugs_slope_low95_vec < true_slope)*(bugs_slope_up95_vec > true_slope);

mean(bugs_slope_coverage);

# variance estimates;
mean(sqrt(bugs_mi_var));
sqrt(var(bugs_mi_mean));

# passive imputation;

pas_mi_bias=mean(pas_mi_mean)-true_slope;
pas_mi_bias;

pas_mi_bias/true_slope;

pas_mi_mse=mean((pas_mi_mean-true_slope)^2);
pas_mi_mse;

# coverage;
pas_slope_low95_vec=pas_mi_mean-qt(.975, pas_mi_df)*sqrt(pas_mi_var);
pas_slope_up95_vec=pas_mi_mean+qt(.975, pas_mi_df)*sqrt(pas_mi_var);

pas_slope_length=mean(pas_slope_up95_vec-pas_slope_low95_vec);
pas_slope_length;

pas_slope_coverage=(pas_slope_low95_vec < true_slope)*(pas_slope_up95_vec > true_slope);

mean(pas_slope_coverage);

# variance estimates;
mean(sqrt(pas_mi_var));
sqrt(var(pas_mi_mean));

# JAV

jav_mi_bias=mean(jav_mi_mean)-true_slope;
jav_mi_bias;

jav_mi_bias/true_slope;

jav_mi_mse=mean((jav_mi_mean-true_slope)^2);
jav_mi_mse;

# coverage;
jav_slope_low95_vec=jav_mi_mean-qt(.975, jav_mi_df)*sqrt(jav_mi_var);
jav_slope_up95_vec=jav_mi_mean+qt(.975, jav_mi_df)*sqrt(jav_mi_var);

jav_slope_length=mean(jav_slope_up95_vec-jav_slope_low95_vec);
jav_slope_length;

jav_slope_coverage=(jav_slope_low95_vec < true_slope)*(jav_slope_up95_vec > true_slope);

mean(jav_slope_coverage);

# variance estimates;
mean(sqrt(jav_mi_var));
sqrt(var(jav_mi_mean));
