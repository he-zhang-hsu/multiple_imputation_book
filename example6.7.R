library(R2WinBUGS);
library(MASS);
library(mice);
library(norm);
library(mix);
# library(HI);
options(digits=4);
rm(list=ls());


# setwd('\\\\cdc.gov\\private\\L728\\book\\example8.4')
# get the winbugs model file;

model.file <- system.file(package="R2WinBUGS", "model", "logitmodel_missing_bivariate_forR.txt")
# Let's take a look:
file.show(model.file)


# generate the logistic outcome;
beta0=-1;
beta1=0.5;
beta2=0.5;

# generate bivariate normal x variables;
mu=rep(0,2);
Sigma=matrix(c(1,0.5,0.5,1), nrow=2,ncol=2);


rowobs=1000;
# cycle_no=1000;
cycle_no=1;
mi_no=50;

# matrices holding the parameter estimates;
# complete-data inferences;
y2_mean_vec=y2_var_vec=cc_mean_vec=cc_var_vec=rep(NA, cycle_no);

bugs_mean_mat=bugs_var_mat=glm_mean_mat=glm_var_mat=matrix(NA, cycle_no, mi_no);

BD_slope_vec=BD_slope_var_vec=rep(NA, cycle_no);
CC_slope_vec=CC_slope_var_vec=rep(NA, cycle_no);
bugs_slope_mat=bugs_slope_var_mat=glm_slope_mat=glm_slope_var_mat=matrix(NA, cycle_no, mi_no);


# generate data for example 8.4
for (cycle in 1:cycle_no)

{
set.seed(cycle);

x1_x2=mvrnorm(n = rowobs, mu=mu, Sigma=Sigma);

x1=x1_x2[,1];
x2=x1_x2[,2];

# logistic outcome;
y2=rbinom(rowobs,size=1,prob=exp(beta0+beta1*x1+beta2*x2)/(1+exp(beta0+beta1*x1+beta2*x2)));

# apply the missing cases to complete data;
# missing completely at random;

missing_indi=cbind(rbinom(rowobs,size=1,prob=0.15), rbinom(rowobs,size=1, prob=0.15), rbinom(rowobs, size=1, prob=0.15));

y2_miss=y2;
x1_miss=x1;
x2_miss=x2;

# y2_miss[miss_indi==1]=NA;
# miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];


y2_miss[missing_indi[,1]==1]=NA;
x1_miss[missing_indi[,2]==1]=NA;
x2_miss[missing_indi[,3]==1]=NA;

y2_obs_no=length(y2_miss[!is.na(y2_miss)]);
y2_mis_no=rowobs-y2_obs_no;
y2_obs=y2[missing_indi[,1]==0];

# before deletion analysis;
# compelete data analysis;
# sample mean;
y2_mean=mean(y2);
y2_var=var(y2)/rowobs;
y2_mean_vec[cycle]=y2_mean;
y2_var_vec[cycle]=y2_var;

# logistic regression coefficient for the slope of x1;
BD_logistic=summary(glm(y2~x1+x2, family=binomial(link="logit")));
BD_slope_vec[cycle]=BD_logistic$coeff[2,1];
BD_slope_var_vec[cycle]=BD_logistic$coeff[2,2]^2;

# complete-case analysis;
# mean estimands;
cc_mean=mean(y2_miss, na.rm="T");
cc_mean_vec[cycle]=cc_mean;
cc_var_vec[cycle]=var(y2_miss, na.rm="T")/y2_obs_no;

# logistic regression coefficient for the slope of x1;
CC_logistic=summary(glm(y2_miss~x1_miss+x2_miss, family=binomial(link="logit")));
CC_slope_vec[cycle]=CC_logistic$coeff[2,1];
CC_slope_var_vec[cycle]=CC_logistic$coeff[2,2]^2;


# winbugs imputation;

# apply missing cases to the data;
M=rowobs;

process_stage2_data_simulate=list(M=M, y2_miss=y2_miss, x1_miss=x1_miss, x2_miss=x2_miss);
process_stage2_init= function(){list(beta0=-1, beta1=0.5, beta2=0.5, mu=0, alpha0=0, alpha1=0, tau=1, psi=1)};
# process_stage2_parameters=c("mu", "beta0", "beta1", "beta2", "mu", "alpha0", "alpha1", "tau", "psi", "y2", "x1", "x2");
process_stage2_parameters=c("y2_impute", "x1_impute", "x2_impute");


gibbs_no=10000;

# imputation Gibbs;

process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file,
    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=gibbs_no/(2*mi_no), bugs.seed=cycle,
    bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);

attach.bugs(process_stage2.sim);
rm(process_stage2.sim);


# parameter draw Gibbs;

#process_stage2_parameters=c("beta0", "beta1", "beta2", "mu", "alpha0", "alpha1", "tau", "psi");

#process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file,
#    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=1, bugs.seed=cycle,
#    bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
#    working.directory=NULL, clearWD=TRUE);

#plot(seq(gibbs_no/2+1, gibbs_no, 1), beta0, type="l", xlab="iteration", main="", ylab="beta0");
#plot(seq(gibbs_no/2+1, gibbs_no, 1), beta1, type="l", xlab="iteration", main="", ylab="beta1");
#plot(seq(gibbs_no/2+1, gibbs_no, 1), beta2, type="l", xlab="iteration", main="", ylab="beta2");
#plot(seq(gibbs_no/2+1, gibbs_no, 1), mu, type="l", xlab="iteration", main="", ylab="mu");
#plot(seq(gibbs_no/2+1, gibbs_no, 1), alpha0, type="l", xlab="iteration", main="", ylab="alpha0");
#plot(seq(gibbs_no/2+1, gibbs_no, 1), alpha1, type="l", xlab="iteration", main="", ylab="alpha1");
#plot(seq(gibbs_no/2+1, gibbs_no, 1), tau, type="l", xlab="iteration", main="", ylab="tau");
#plot(seq(gibbs_no/2+1, gibbs_no, 1), psi, type="l", xlab="iteration", main="", ylab="psi");


attach.bugs(process_stage2.sim);
rm(process_stage2.sim);











y2_completed_bugs=y2_miss;
x1_completed_bugs=x1_miss;
x2_completed_bugs=x2_miss;

for (i in 1:mi_no)
{
y2_completed_bugs[missing_indi[,1]==1]=y2_impute[i,];
x1_completed_bugs[missing_indi[,2]==1]=x1_impute[i,];
x2_completed_bugs[missing_indi[,3]==1]=x2_impute[i,];


# marginal means;
bugs_mean=mean(y2_completed_bugs);
bugs_var=var(y2_completed_bugs)/rowobs;

bugs_mean_mat[cycle,i]=bugs_mean;
bugs_var_mat[cycle,i]=bugs_var;

# logistic regression coefficient for x1 slope;
bugs_logistic=summary(glm(y2_completed_bugs~x1_completed_bugs+x2_completed_bugs, family=binomial(link="logit")));
bugs_slope_mat[cycle,i]=bugs_logistic$coeff[2,1];
bugs_slope_var_mat[cycle,i]=bugs_logistic$coeff[2,2]^2;

}

# general location modeling imputation;
y2_miss_recode=y2_miss+1;
ori_data=cbind(y2_miss_recode, x1_miss, x2_miss);
s=prelim.mix(ori_data,1);
thetahat=em.mix(s);

for (i in 1:mi_no)
{
rngseed(i); # set random number generator seed
newtheta=da.mix(s, thetahat, steps=200, showits=TRUE) # take 200 steps
# getparam.mix(s, newtheta, corr=TRUE);
imp_data=imp.mix(s, newtheta);

y2_completed_glm=imp_data[,1]-1;
x1_completed_glm=imp_data[,2];
x2_completed_glm=imp_data[,3];


# multiple imputation analysis;

# marginal means;
glm_mean=mean(y2_completed_glm);
glm_var=var(y2_completed_glm)/rowobs;

glm_mean_mat[cycle,i]=glm_mean;
glm_var_mat[cycle,i]=glm_var;

# logistic regression coefficient for x1 slope;
glm_logistic=summary(glm(y2_completed_glm~x1_completed_glm+x2_completed_glm, family=binomial(link="logit")));
glm_slope_mat[cycle,i]=glm_logistic$coeff[2,1];
glm_slope_var_mat[cycle,i]=glm_logistic$coeff[2,2]^2;



}

cat("the cycle is", cycle, "\n");

}

# the end of mulitple imputation;
#############################################################################
# marginal mean estimand;

true=mean(y2_mean_vec);
true;

mean_bias=mean(y2_mean_vec-true);
mean_bias;

mean_bias/true;

true_mse=mean((y2_mean_vec-true)^2);
true_mse;

# lower and upper 95% CI;
y2_mean_low95_vec=y2_mean_vec-1.96*sqrt(y2_var_vec);
y2_mean_up95_vec=y2_mean_vec+1.96*sqrt(y2_var_vec);

mean_length=mean(y2_mean_up95_vec-y2_mean_low95_vec);
mean_length;

coverage=(y2_mean_low95_vec < true)*(y2_mean_up95_vec > true);

mean(coverage);



# variance estimates;
mean(sqrt(y2_var_vec));
sqrt(var(y2_mean_vec));


# complete case analysis;
cc_true=mean(cc_mean_vec);
cc_bias=cc_true-true;
cc_bias;

cc_bias/true;

cc_mse=mean((cc_mean_vec-true)^2);
cc_mse;

# lower and upper 95% CI;
cc_mean_low95_vec=cc_mean_vec-1.96*sqrt(cc_var_vec);
cc_mean_up95_vec=cc_mean_vec+1.96*sqrt(cc_var_vec);

mean_cc_length=mean(cc_mean_up95_vec-cc_mean_low95_vec);
mean_cc_length;

cc_coverage=(cc_mean_low95_vec < true)*(cc_mean_up95_vec > true);

mean(cc_coverage);

# variance estimates;
mean(sqrt(cc_var_vec));
sqrt(var(cc_mean_vec));


# multiple imputation estimates;
# bugs imputation

bugs_mi_mean=rep(NA, cycle_no);
bugs_mi_var=rep(NA, cycle_no);
bugs_mi_df=rep(NA, cycle_no);
bugs_mi_f=rep(NA, cycle_no);

# glm imputation;

glm_mi_mean=rep(NA, cycle_no);
glm_mi_var=rep(NA, cycle_no);
glm_mi_df=rep(NA, cycle_no);
glm_mi_f=rep(NA, cycle_no);

for (i in 1:cycle_no)
{
bugs_summary=pool.scalar(bugs_mean_mat[i,], bugs_var_mat[i,], n=rowobs);
bugs_mi_mean[i]=bugs_summary$qbar;
bugs_mi_var[i]=bugs_summary$t;
bugs_mi_df[i]=bugs_summary$df;
bugs_mi_f[i]=bugs_summary$f;

glm_summary=pool.scalar(glm_mean_mat[i,], glm_var_mat[i,], n=rowobs);
glm_mi_mean[i]=glm_summary$qbar;
glm_mi_var[i]=glm_summary$t;
glm_mi_df[i]=glm_summary$df;
glm_mi_f[i]=glm_summary$f;



}

# bugs imputation;
bugs_mi_true=mean(bugs_mi_mean);
bugs_mi_bias=bugs_mi_true-true;

bugs_mi_bias;

bugs_mi_bias/true;

bugs_mi_mse=mean((bugs_mi_mean-true)^2);
bugs_mi_mse;

# coverage;
bugs_mean_low95_vec=bugs_mi_mean-qt(.975, bugs_mi_df)*sqrt(bugs_mi_var);
bugs_mean_up95_vec=bugs_mi_mean+qt(.975, bugs_mi_df)*sqrt(bugs_mi_var);

mean_bugs_length=mean(bugs_mean_up95_vec-bugs_mean_low95_vec);
mean_bugs_length;

bugs_mi_coverage=(bugs_mean_low95_vec < true)*(bugs_mean_up95_vec > true);
mean(bugs_mi_coverage);

# variance estimates;
mean(sqrt(bugs_mi_var));
sqrt(var(bugs_mi_mean));


# glm imputation;
glm_mi_true=mean(glm_mi_mean);
glm_mi_bias=glm_mi_true-true;

glm_mi_bias;
glm_mi_bias/true;

glm_mi_mse=mean((glm_mi_mean-true)^2);
glm_mi_mse;

# coverage;
glm_mean_low95_vec=glm_mi_mean-qt(.975, glm_mi_df)*sqrt(glm_mi_var);
glm_mean_up95_vec=glm_mi_mean+qt(.975, glm_mi_df)*sqrt(glm_mi_var);

# trans_mean_low95_vec=trans_mean_low95_vec[!is.na(trans_mean_low95_vec)];
# trans_mean_up95_vec=trans_mean_up95_vec[!is.na(trans_mean_up95_vec)];


mean_glm_length=mean(glm_mean_up95_vec-glm_mean_low95_vec);
mean_glm_length;

glm_mi_coverage=(glm_mean_low95_vec < true)*(glm_mean_up95_vec > true);
mean(glm_mi_coverage);

# variance estimates;
mean(sqrt(glm_mi_var));
sqrt(var(glm_mi_mean));


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

glm_mi_mean=rep(NA, cycle_no);
glm_mi_var=rep(NA, cycle_no);
glm_mi_df=rep(NA, cycle_no);
glm_mi_f=rep(NA, cycle_no);


for (i in 1:cycle_no)
{
bugs_summary=pool.scalar(bugs_slope_mat[i,], bugs_slope_var_mat[i,], n=rowobs);
bugs_mi_mean[i]=bugs_summary$qbar;
bugs_mi_var[i]=bugs_summary$t;
bugs_mi_df[i]=bugs_summary$df;
bugs_mi_f[i]=bugs_summary$f;

glm_summary=pool.scalar(glm_slope_mat[i,], glm_slope_var_mat[i,], n=rowobs);
glm_mi_mean[i]=glm_summary$qbar;
glm_mi_var[i]=glm_summary$t;
glm_mi_df[i]=glm_summary$df;
glm_mi_f[i]=glm_summary$f;


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

glm_mi_bias=mean(glm_mi_mean)-true_slope;
glm_mi_bias;

glm_mi_bias/true_slope;

glm_mi_mse=mean((glm_mi_mean-true_slope)^2);
glm_mi_mse;

# coverage;
glm_slope_low95_vec=glm_mi_mean-qt(.975, glm_mi_df)*sqrt(glm_mi_var);
glm_slope_up95_vec=glm_mi_mean+qt(.975, glm_mi_df)*sqrt(glm_mi_var);

glm_slope_length=mean(glm_slope_up95_vec-glm_slope_low95_vec);
glm_slope_length;

glm_slope_coverage=(glm_slope_low95_vec < true_slope)*(glm_slope_up95_vec > true_slope);

mean(glm_slope_coverage);

# variance estimates;
mean(sqrt(glm_mi_var));
sqrt(var(glm_mi_mean));
