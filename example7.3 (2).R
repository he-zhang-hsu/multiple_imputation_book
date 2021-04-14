################## data generate ##############
# Program: boxcox_impute_expolore_simulation; ############
# generate data, and try out different missing data methods;
# run some simulations to assess the biasedness of various methods;

# library(car);
# library(mvtnorm);
library(mice);
library(norm);
library(R2WinBUGS);
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



model0.file <- system.file(package="R2WinBUGS", "model", "normal_fcs.txt")
# Let's take a look:
file.show(model0.file)

model1.file <- system.file(package="R2WinBUGS", "model", "normal_jm.txt")
# Let's take a look:
file.show(model1.file)


model2.file <- system.file(package="R2WinBUGS", "model", "normal_fcs_cut.txt")
# Let's take a look:
file.show(model2.file)



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


cycle_no=1000;
mi_no=50;

# matrices holding the parameter estimates;
# complete-data inferences;
# xyslope means regressing x on y;
BD_mean_vec=BD_mean_var_vec=BD_slope_vec=BD_slope_var_vec=BD_xyslope_vec=BD_xyslope_var_vec=rep(NA, cycle_no);
CC_mean_vec=CC_mean_var_vec=CC_slope_vec=CC_slope_var_vec=CC_xyslope_vec=CC_xyslope_var_vec=rep(NA, cycle_no);
SI_mean_vec=SI_mean_var_vec=SI_slope_vec=SI_slope_var_vec=SI_xyslope_vec=SI_xyslope_var_vec=rep(NA, cycle_no);
INDI_xyslope_vec=INDI_xyslope_var_vec=rep(NA, cycle_no);

MI_mean_mat=MI_mean_var_mat=MI_slope_mat=MI_slope_var_mat=MI_xyslope_mat=MI_xyslope_var_mat=matrix(NA, cycle_no, mi_no);
MICE_mean_mat=MICE_mean_var_mat=MICE_slope_mat=MICE_slope_var_mat=MICE_xyslope_mat=MICE_xyslope_var_mat=matrix(NA, cycle_no, mi_no);
BUGS_mean_mat=BUGS_mean_var_mat=BUGS_slope_mat=BUGS_slope_var_mat=BUGS_xyslope_mat=BUGS_xyslope_var_mat=matrix(NA, cycle_no, mi_no);



for (cycle in 1:cycle_no)

{

set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# error distribution
# could be normal, or could be other distributions;
error=sqrt(var_error)*rnorm(rowobs);

# generate the latent variable;
y=beta0+beta1*x+error;

# MCAR;

missing_indi=cbind(rbinom(rowobs,size=1,prob=0.20), rbinom(rowobs,size=1, prob=0.20));

missing_indi_x=missing_indi[,1];
missing_indi_y=missing_indi[,2];

y_miss=y;
x_miss=x;

# y2_miss[miss_indi==1]=NA;
# miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];

x_miss[missing_indi_x==1]=NA;
y_miss[missing_indi_y==1]=NA;

# write.table(cbind(y_miss, x_miss), row.name=FALSE, file="C:/Users/Guanghui He/Personal/Yulei/research/MIbook/chapter7_multivariate_MICE/program/xytest.txt");

# summary(lm(x~y));

# set up the missing data as MCAR on x;
# miss_indi=runif(n=rowobs)<0.40;

# set up the missing data as MAR on x;
# alpha0=-1.4;
# alpha1=1;
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x)/(1+exp(alpha0+alpha1*x)));


# obs_indi=1-miss_indi;

# y_miss=y;
# y_miss[miss_indi==1]=NA;
# miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
# obs_no=length(y_miss[!is.na(y_miss)]);
# mis_no=rowobs-obs_no;
# y_obs=y[miss_indi==0];
# x_obs=x[miss_indi==0];
# x_mis=x[miss_indi==1];

# compelete data analysis;
# marginal mean;
# BD_mean_vec[cycle]=mean(y);
# BD_mean_var_vec[cycle]=var(y)/rowobs;


# regression analysis;
BD=summary(lm(y~x));

BD_slope_vec[cycle]=BD$coeff[2,1];
BD_slope_var_vec[cycle]=BD$coeff[2,2]^2;

# reverse regression analysis;
BD_xy=summary(lm(x~y));

BD_xyslope_vec[cycle]=BD_xy$coeff[2,1];
BD_xyslope_var_vec[cycle]=BD_xy$coeff[2,2]^2;


# different missing data methods;

# complete-case analysis;

# marginal mean;
# CC_mean_vec[cycle]=mean(y_obs);
# CC_mean_var_vec[cycle]=var(y_obs)/obs_no;


# regression analysis;
CC=summary(lm(y_miss~x_miss));

CC_slope_vec[cycle]=CC$coeff[2,1];
CC_slope_var_vec[cycle]=CC$coeff[2,2]^2;

# reverse regression analysis;
CC_xy=summary(lm(x_miss~y_miss));

CC_xyslope_vec[cycle]=CC_xy$coeff[2,1];
CC_xyslope_var_vec[cycle]=CC_xy$coeff[2,2]^2;





# multiple imputation analysis;
# joint model imputation;


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

# MICE imputation;
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

# plot the distribution of the data;

# check the performance of the estimates;
# population quantity

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



# multiple imputation;

slope_mi_mean=rep(NA, cycle_no);
slope_mi_var=rep(NA, cycle_no);
slope_mi_df=rep(NA, cycle_no);
slope_mi_f=rep(NA, cycle_no);

slope_mice_mean=rep(NA, cycle_no);
slope_mice_var=rep(NA, cycle_no);
slope_mice_df=rep(NA, cycle_no);
slope_mice_f=rep(NA, cycle_no);

slope_bugs_mean=rep(NA, cycle_no);
slope_bugs_var=rep(NA, cycle_no);
slope_bugs_df=rep(NA, cycle_no);
slope_bugs_f=rep(NA, cycle_no);




for (i in 1:cycle_no)
{
slope_summary=pool.scalar(MI_slope_mat[i,], MI_slope_var_mat[i,], n=rowobs);
slope_mi_mean[i]=slope_summary$qbar;
slope_mi_var[i]=slope_summary$t;
slope_mi_df[i]=slope_summary$df;
slope_mi_f[i]=slope_summary$f;


slope_summary=pool.scalar(MICE_slope_mat[i,], MICE_slope_var_mat[i,], n=rowobs);
slope_mice_mean[i]=slope_summary$qbar;
slope_mice_var[i]=slope_summary$t;
slope_mice_df[i]=slope_summary$df;
slope_mice_f[i]=slope_summary$f;

slope_summary=pool.scalar(BUGS_slope_mat[i,], BUGS_slope_var_mat[i,], n=rowobs);
slope_bugs_mean[i]=slope_summary$qbar;
slope_bugs_var[i]=slope_summary$t;
slope_bugs_df[i]=slope_summary$df;
slope_bugs_f[i]=slope_summary$f;



}


# JM imputation;

slope_mi_bias=mean(slope_mi_mean)-true_slope;
slope_mi_bias;
slope_mi_bias/true_slope;

slope_mi_mse=mean((slope_mi_mean-true_slope)^2);
slope_mi_mse;

# coverage;
MI_slope_low95_vec=slope_mi_mean-qt(.975, slope_mi_df)*sqrt(slope_mi_var);
MI_slope_up95_vec=slope_mi_mean+qt(.975, slope_mi_df)*sqrt(slope_mi_var);

MI_slope_length=mean(MI_slope_up95_vec-MI_slope_low95_vec);
MI_slope_length;


MI_slope_coverage=(MI_slope_low95_vec < true_slope)*(MI_slope_up95_vec > true_slope);

mean(MI_slope_coverage);

# variance estimates;
mean(sqrt(slope_mi_var));
sqrt(var(slope_mi_mean));


# MICE imputation;

slope_mice_bias=mean(slope_mice_mean)-true_slope;
slope_mice_bias;
slope_mice_bias/true_slope;

slope_mice_mse=mean((slope_mice_mean-true_slope)^2);
slope_mice_mse;

# coverage;
MICE_slope_low95_vec=slope_mice_mean-qt(.975, slope_mice_df)*sqrt(slope_mice_var);
MICE_slope_up95_vec=slope_mice_mean+qt(.975, slope_mice_df)*sqrt(slope_mice_var);

MICE_slope_length=mean(MICE_slope_up95_vec-MICE_slope_low95_vec);
MICE_slope_length;


MICE_slope_coverage=(MICE_slope_low95_vec < true_slope)*(MICE_slope_up95_vec > true_slope);

mean(MICE_slope_coverage);

# variance estimates;
mean(sqrt(slope_mice_var));
sqrt(var(slope_mice_mean));


# BUGS imputation;

slope_bugs_bias=mean(slope_bugs_mean)-true_slope;
slope_bugs_bias;
slope_bugs_bias/true_slope;

slope_bugs_mse=mean((slope_bugs_mean-true_slope)^2);
slope_bugs_mse;

# coverage;
BUGS_slope_low95_vec=slope_bugs_mean-qt(.975, slope_bugs_df)*sqrt(slope_bugs_var);
BUGS_slope_up95_vec=slope_bugs_mean+qt(.975, slope_bugs_df)*sqrt(slope_bugs_var);

BUGS_slope_length=mean(BUGS_slope_up95_vec-BUGS_slope_low95_vec);
BUGS_slope_length;


BUGS_slope_coverage=(BUGS_slope_low95_vec < true_slope)*(BUGS_slope_up95_vec > true_slope);

mean(BUGS_slope_coverage);

# variance estimates;
mean(sqrt(slope_mice_var));
sqrt(var(slope_mice_mean));


#######################################################################
# reverse slope;

# complete-data analysis;

true_xyslope=mean(BD_xyslope_vec);
true_xyslope;


# performance;
BD_xyslope_mse=mean((BD_xyslope_vec-true_xyslope)^2);
BD_xyslope_mse;

# lower and upper 95% CI;
BD_xyslope_low95_vec=BD_xyslope_vec-1.96*sqrt(BD_xyslope_var_vec);
BD_xyslope_up95_vec=BD_xyslope_vec+1.96*sqrt(BD_xyslope_var_vec);

BD_xyslope_length=mean(BD_xyslope_up95_vec-BD_xyslope_low95_vec);
BD_xyslope_length;

BD_xyslope_cov=(BD_xyslope_low95_vec < true_xyslope)*(BD_xyslope_up95_vec > true_xyslope);
mean(BD_xyslope_cov);

# variance estimates;
mean(sqrt(BD_xyslope_var_vec));
sqrt(var(BD_xyslope_vec));


# complete-case analysis;
CC_xyslope_bias=mean(CC_xyslope_vec)-true_xyslope;
CC_xyslope_bias;

CC_xyslope_bias/true_xyslope;

CC_xyslope_mse=mean((CC_xyslope_vec-true_xyslope)^2);
CC_xyslope_mse;

# lower and upper 95% CI;
CC_xyslope_low95_vec=CC_xyslope_vec-1.96*sqrt(CC_xyslope_var_vec);
CC_xyslope_up95_vec=CC_xyslope_vec+1.96*sqrt(CC_xyslope_var_vec);

CC_xyslope_length=mean(CC_xyslope_up95_vec-CC_xyslope_low95_vec);
CC_xyslope_length;

CC_xyslope_cov=(CC_xyslope_low95_vec < true_xyslope)*(CC_xyslope_up95_vec > true_xyslope);

mean(CC_xyslope_cov);

# variance estimates;
mean(sqrt(CC_xyslope_var_vec));
sqrt(var(CC_xyslope_vec));
# multiple imputation estimates;

xyslope_mi_mean=rep(NA, cycle_no);
xyslope_mi_var=rep(NA, cycle_no);
xyslope_mi_df=rep(NA, cycle_no);
xyslope_mi_f=rep(NA, cycle_no);

xyslope_mice_mean=rep(NA, cycle_no);
xyslope_mice_var=rep(NA, cycle_no);
xyslope_mice_df=rep(NA, cycle_no);
xyslope_mice_f=rep(NA, cycle_no);

xyslope_bugs_mean=rep(NA, cycle_no);
xyslope_bugs_var=rep(NA, cycle_no);
xyslope_bugs_df=rep(NA, cycle_no);
xyslope_bugs_f=rep(NA, cycle_no);



for (i in 1:cycle_no)
{
xyslope_summary=pool.scalar(MI_xyslope_mat[i,], MI_xyslope_var_mat[i,], n=rowobs);
xyslope_mi_mean[i]=xyslope_summary$qbar;
xyslope_mi_var[i]=xyslope_summary$t;
xyslope_mi_df[i]=xyslope_summary$df;
xyslope_mi_f[i]=xyslope_summary$f;

xyslope_summary=pool.scalar(MICE_xyslope_mat[i,], MICE_xyslope_var_mat[i,], n=rowobs);
xyslope_mice_mean[i]=xyslope_summary$qbar;
xyslope_mice_var[i]=xyslope_summary$t;
xyslope_mice_df[i]=xyslope_summary$df;
xyslope_mice_f[i]=xyslope_summary$f;

xyslope_summary=pool.scalar(BUGS_xyslope_mat[i,], BUGS_xyslope_var_mat[i,], n=rowobs);
xyslope_bugs_mean[i]=xyslope_summary$qbar;
xyslope_bugs_var[i]=xyslope_summary$t;
xyslope_bugs_df[i]=xyslope_summary$df;
xyslope_bugs_f[i]=xyslope_summary$f;




}



# normal imputation;

xyslope_mi_bias=mean(xyslope_mi_mean)-true_xyslope;
xyslope_mi_bias;
xyslope_mi_bias/true_xyslope;

xyslope_mi_mse=mean((xyslope_mi_mean-true_xyslope)^2);
xyslope_mi_mse;

# coverage;
MI_xyslope_low95_vec=xyslope_mi_mean-qt(.975, xyslope_mi_df)*sqrt(xyslope_mi_var);
MI_xyslope_up95_vec=xyslope_mi_mean+qt(.975, xyslope_mi_df)*sqrt(xyslope_mi_var);

MI_xyslope_length=mean(MI_xyslope_up95_vec-MI_xyslope_low95_vec);
MI_xyslope_length;


MI_xyslope_coverage=(MI_xyslope_low95_vec < true_xyslope)*(MI_xyslope_up95_vec > true_xyslope);

mean(MI_xyslope_coverage);

# variance estimates;
mean(sqrt(xyslope_mi_var));
sqrt(var(xyslope_mi_mean));


# MICE imputation;

xyslope_mice_bias=mean(xyslope_mice_mean)-true_xyslope;
xyslope_mice_bias;
xyslope_mice_bias/true_xyslope;

xyslope_mice_mse=mean((xyslope_mice_mean-true_xyslope)^2);
xyslope_mice_mse;

# coverage;
MICE_xyslope_low95_vec=xyslope_mice_mean-qt(.975, xyslope_mice_df)*sqrt(xyslope_mice_var);
MICE_xyslope_up95_vec=xyslope_mice_mean+qt(.975, xyslope_mice_df)*sqrt(xyslope_mice_var);

MICE_xyslope_length=mean(MICE_xyslope_up95_vec-MICE_xyslope_low95_vec);
MICE_xyslope_length;


MICE_xyslope_coverage=(MICE_xyslope_low95_vec < true_xyslope)*(MICE_xyslope_up95_vec > true_xyslope);

mean(MICE_xyslope_coverage);

# variance estimates;
mean(sqrt(xyslope_mice_var));
sqrt(var(xyslope_mice_mean));

# BUGS imputation;


xyslope_bugs_bias=mean(xyslope_bugs_mean)-true_xyslope;
xyslope_bugs_bias;
xyslope_bugs_bias/true_xyslope;

xyslope_bugs_mse=mean((xyslope_bugs_mean-true_xyslope)^2);
xyslope_bugs_mse;

# coverage;
BUGS_xyslope_low95_vec=xyslope_bugs_mean-qt(.975, xyslope_bugs_df)*sqrt(xyslope_bugs_var);
BUGS_xyslope_up95_vec=xyslope_bugs_mean+qt(.975, xyslope_bugs_df)*sqrt(xyslope_bugs_var);

BUGS_xyslope_length=mean(BUGS_xyslope_up95_vec-BUGS_xyslope_low95_vec);
BUGS_xyslope_length;


BUGS_xyslope_coverage=(BUGS_xyslope_low95_vec < true_xyslope)*(BUGS_xyslope_up95_vec > true_xyslope);

mean(BUGS_xyslope_coverage);

# variance estimates;
mean(sqrt(xyslope_bugs_var));
sqrt(var(xyslope_bugs_mean));





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


