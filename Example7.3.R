# Example 7.3

rm(list=ls());
# library(car);
# library(mvtnorm);
library(mice);
library(norm);
library(R2WinBUGS);
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


# WinBUGS code using cut() to perform
# FCS imputation;

model2.file <- system.file(package="R2WinBUGS", "model", "normal_fcs_cut.txt")
# Let's take a look:
file.show(model2.file)



# set up the random seed;
set.seed(197789);


# complete-data sample size;

rowobs=1000;

# mean and variance of x;
mu_x=1;
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

# Vectors and matrices holding the parameter estimates;
# xyslope means regressing x on y;

BD_slope_vec=BD_slope_var_vec=BD_xyslope_vec=BD_xyslope_var_vec=rep(NA, cycle_no);

CC_mean_vec=CC_mean_var_vec=CC_slope_vec=CC_slope_var_vec=CC_xyslope_vec=CC_xyslope_var_vec=rep(NA, cycle_no);

# Multiple imputation
# MI for JM imputation;
# MICE for the FCS imputation by R MICE
# BUGS for the FCS imputation by WinBUGS

MI_slope_mat=MI_slope_var_mat=MI_xyslope_mat=MI_xyslope_var_mat=matrix(NA, cycle_no, mi_no);

MICE_slope_mat=MICE_slope_var_mat=MICE_xyslope_mat=MICE_xyslope_var_mat=matrix(NA, cycle_no, mi_no);

BUGS_slope_mat=BUGS_slope_var_mat=BUGS_xyslope_mat=BUGS_xyslope_var_mat=matrix(NA, cycle_no, mi_no);


# begin the simulation;

for (cycle in 1:cycle_no)

{

set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the y variable;
y=beta0+beta1*x+error;

# MCAR;
# pose missing data for both x and y;

missing_indi=cbind(rbinom(rowobs,size=1,prob=0.20), rbinom(rowobs,size=1, prob=0.20));

missing_indi_x=missing_indi[,1];
missing_indi_y=missing_indi[,2];

y_miss=y;
x_miss=x;

x_miss[missing_indi_x==1]=NA;
y_miss[missing_indi_y==1]=NA;


# regression analysis;
BD=summary(lm(y~x));

BD_slope_vec[cycle]=BD$coeff[2,1];
BD_slope_var_vec[cycle]=BD$coeff[2,2]^2;

# reverse regression analysis;
# regressing x on y;
BD_xy=summary(lm(x~y));

BD_xyslope_vec[cycle]=BD_xy$coeff[2,1];
BD_xyslope_var_vec[cycle]=BD_xy$coeff[2,2]^2;


# complete-case analysis;

# regression analysis;
CC=summary(lm(y_miss~x_miss));

CC_slope_vec[cycle]=CC$coeff[2,1];
CC_slope_var_vec[cycle]=CC$coeff[2,2]^2;

# reverse regression analysis;
CC_xy=summary(lm(x_miss~y_miss));

CC_xyslope_vec[cycle]=CC_xy$coeff[2,1];
CC_xyslope_var_vec[cycle]=CC_xy$coeff[2,2]^2;

# multiple imputation analysis;
# JM imputation;
# by R norm

ori_data=cbind(y_miss, x_miss);
s=prelim.norm(ori_data);
thetahat=em.norm(s);

for (i in 1:mi_no)
{
rngseed(i); # set random number generator seed
newtheta=da.norm(s, thetahat, steps=200, showits=TRUE) # take 200 steps
imp_data=imp.norm(s, newtheta);

y_completed_norm=imp_data[,1];
x_completed_norm=imp_data[,2];

# regression analysis;

IMP=summary(lm(y_completed_norm~x_completed_norm));

MI_slope_mat[cycle,i]=IMP$coeff[2,1];
MI_slope_var_mat[cycle,i]=IMP$coeff[2,2]^2;

# reverse regression analysis;

IMP_xy=summary(lm(x_completed_norm~y_completed_norm));

MI_xyslope_mat[cycle,i]=IMP_xy$coeff[2,1];
MI_xyslope_var_mat[cycle,i]=IMP_xy$coeff[2,2]^2;

}

# FCS imputation by R MICE

ori_data_frame=as.data.frame(ori_data);
y_completed_mice=y_miss;
x_completed_mice=x_miss;
for (i in 1:mi_no)
{

 u=mice(data=ori_data_frame, m=1, method=c('norm', 'norm'), seed=i);
 y_imputed_mice=u$imp$y_miss[[1]];
 y_completed_mice[missing_indi_y==1]=y_imputed_mice;
 
x_imputed_mice=u$imp$x_miss[[1]];
x_completed_mice[missing_indi_x==1]=x_imputed_mice;

# regression analysis;

IMP=summary(lm(y_completed_mice~x_completed_mice));

MICE_slope_mat[cycle,i]=IMP$coeff[2,1];
MICE_slope_var_mat[cycle,i]=IMP$coeff[2,2]^2;

# reverse regression analysis;

IMP_xy=summary(lm(x_completed_mice~y_completed_mice));

MICE_xyslope_mat[cycle,i]=IMP_xy$coeff[2,1];
MICE_xyslope_var_mat[cycle,i]=IMP_xy$coeff[2,2]^2;


}

# bugs imputation for FCS
M=rowobs;

gibbs_no=10000;


process_stage2_data_simulate_JM=list("M","y_miss", "x_miss");
process_stage2_init_JM=function(){list(beta0=-2, beta1=1, alpha0=1.5, alpha1=0.5, tau_x=2, tau_y=1)};
# process_stage2_init_JM=function(){list(beta0=-2, beta1=1, mu=1, tau_x=2, tau_y=1)};

process_stage2_parameters_JM=c("y_impute", "x_impute");

process_stage2.sim.JM <- bugs(data=process_stage2_data_simulate_JM, inits=process_stage2_init_JM, parameters=process_stage2_parameters_JM, model.file=model2.file,
    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=gibbs_no/(2*mi_no), bugs.seed=cycle,
    bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE,
    working.directory=NULL, clearWD=TRUE, debug=FALSE);

attach.bugs(process_stage2.sim.JM);
rm(process_stage2.sim.JM);



y_completed_bugs=y_miss;
x_completed_bugs=x_miss;

for (i in 1:mi_no)
{
y_completed_bugs[missing_indi_y==1]=y_impute[i,];
x_completed_bugs[missing_indi_x==1]=x_impute[i,];

# regression analysis;

IMP=summary(lm(y_completed_bugs~x_completed_bugs));

BUGS_slope_mat[cycle,i]=IMP$coeff[2,1];
BUGS_slope_var_mat[cycle,i]=IMP$coeff[2,2]^2;

# reverse regression analysis;

IMP_xy=summary(lm(x_completed_bugs~y_completed_bugs));

BUGS_xyslope_mat[cycle,i]=IMP_xy$coeff[2,1];
BUGS_xyslope_var_mat[cycle,i]=IMP_xy$coeff[2,2]^2;



}



cat("the cycle is", cycle, "\n");




}

# the end of multiple imputation;


##########################################################################
# Simulation results
# Table 7.3
# Slope of regressing y on x;

true_slope=mean(BD_slope_vec);
true_slope;

# before-deletion;
BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;

# complete-case analysis;
CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;

# JM imputation;
MI=mi_performance(MI_slope_mat, MI_slope_var_mat, true_slope, rowobs);
MI;

# FCS imputation by R MICE;
MICE_MI=mi_performance(MICE_slope_mat, MICE_slope_var_mat, true_slope, rowobs);
MICE_MI;


# FCS imputation by WinBUGS;
BUGS_MI=mi_performance(BUGS_slope_mat, BUGS_slope_var_mat, true_slope, rowobs);
BUGS_MI;

#######################################################################
# reverse slope;

# complete-data analysis;

true_xyslope=mean(BD_xyslope_vec);
true_xyslope;

# before-deletion;
BD=performance(BD_xyslope_vec, BD_xyslope_var_vec, true_xyslope);
BD;

# complete-case analysis;
CC=performance(CC_xyslope_vec, CC_xyslope_var_vec, true_xyslope);
CC;

# JM imputation;
MI=mi_performance(MI_xyslope_mat, MI_xyslope_var_mat, true_xyslope, rowobs);
MI;

# FCS imputation by R MICE;
MICE_MI=mi_performance(MICE_xyslope_mat, MICE_xyslope_var_mat, true_xyslope, rowobs);
MICE_MI;


# FCS imputation by WinBUGS;
BUGS_MI=mi_performance(BUGS_xyslope_mat, BUGS_xyslope_var_mat, true_xyslope, rowobs);
BUGS_MI;






