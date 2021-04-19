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
beta2=1;
beta3=1;

# error variance;
var_error=1;


cycle_no=1000;
mi_no=50;

# matrices holding the parameter estimates;
# complete-data inferences;
# xyslope means regressing x on y;
BD_mean_vec=BD_mean_var_vec=CC_mean_vec=CC_mean_var_vec=p_mis_vec=rep(NA, cycle_no);
MI_mean_x1_mat=MI_mean_x1_var_mat=MI_mean_x1plusx2_mat=MI_mean_x1plusx2_var_mat=MI_mean_x1x2_mat=MI_mean_x1x2_var_mat=matrix(NA, cycle_no, mi_no);
MI_rsquare_x1_mat=MI_rsquare_x1plusx2_mat=MI_rsquare_x1x2_mat=matrix(NA, cycle_no, mi_no);




for (cycle in 1:cycle_no)

{
set.seed(cycle);

x1=mu_x+sqrt(var_x)*rnorm(rowobs);
x2=mu_x+sqrt(var_x)*rnorm(rowobs);


# error distribution
# could be normal, or could be other distributions;
error=sqrt(var_error)*rnorm(rowobs);

# generate the latent variable;
y=beta0+beta1*x1+beta2*x2+beta3*x1*x2+error;

# set up the missing data as MCAR on x;
miss_indi=runif(n=rowobs)<0.40;

# set up the missing data as MAR on x;
# alpha0=-1.4;
# alpha1=0.3;
# alpha2=0.7;
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x1+alpha2*x2)/(1+exp(alpha0+alpha1*x1+alpha2*x2)));


y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
# x_obs=x[miss_indi==0];
# x_mis=x[miss_indi==1];
p_mis_vec[cycle]=mis_no/rowobs;

# compelete data analysis;
# marginal mean;
BD_mean_vec[cycle]=mean(y);
BD_mean_var_vec[cycle]=var(y)/rowobs;


# regression analysis;
# BD=summary(lm(y~x));

# BD_slope_vec[cycle]=BD$coeff[2,1];
# BD_slope_var_vec[cycle]=BD$coeff[2,2]^2;

# reverse regression analysis;
# BD_xy=summary(lm(x~y));

# BD_xyslope_vec[cycle]=BD_xy$coeff[2,1];
# BD_xyslope_var_vec[cycle]=BD_xy$coeff[2,2]^2;


# different missing data methods;

# complete-case analysis;

# marginal mean;
CC_mean_vec[cycle]=mean(y_obs);
CC_mean_var_vec[cycle]=var(y_obs)/obs_no;


# regression analysis;
# CC=summary(lm(y_obs~x_obs));

# CC_slope_vec[cycle]=CC$coeff[2,1];
# CC_slope_var_vec[cycle]=CC$coeff[2,2]^2;

# reverse regression analysis;
# CC_xy=summary(lm(x_obs~y_obs));

# CC_xyslope_vec[cycle]=CC_xy$coeff[2,1];
# CC_xyslope_var_vec[cycle]=CC_xy$coeff[2,2]^2;



# now impute the missing y2's to get estimates of the marginal;
y_completed_x1=y_completed_x1plusx2=y_completed_x1x2=y_miss;

for (i in 1:mi_no)
{
# imputation model only includes x1;
y_imputed_x1=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x1);

y_completed_x1[miss_seq]=y_imputed_x1;

# completed-data analysis;

# marginal mean;
MI_mean_x1_mat[cycle,i]=mean(y_completed_x1);
MI_mean_x1_var_mat[cycle,i]=var(y_completed_x1)/rowobs;
MI_rsquare_x1_mat[cycle,i]=summary(lm(y_completed_x1~x1))$r.squared;

# imputation model includes both x1 and x2;
y_imputed_x1plusx2=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x1,x2));

y_completed_x1plusx2[miss_seq]=y_imputed_x1plusx2;

# completed-data analysis;

# marginal mean;
MI_mean_x1plusx2_mat[cycle,i]=mean(y_completed_x1plusx2);
MI_mean_x1plusx2_var_mat[cycle,i]=var(y_completed_x1plusx2)/rowobs;
MI_rsquare_x1plusx2_mat[cycle,i]=summary(lm(y_completed_x1plusx2~x1+x2))$r.squared;

# imputation model includes both x1 and x2 and their interactions;
y_imputed_x1x2=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x1,x2, x1*x2));

y_completed_x1x2[miss_seq]=y_imputed_x1x2;

# completed-data analysis;

# marginal mean;
MI_mean_x1x2_mat[cycle,i]=mean(y_completed_x1x2);
MI_mean_x1x2_var_mat[cycle,i]=var(y_completed_x1x2)/rowobs;
MI_rsquare_x1x2_mat[cycle,i]=summary(lm(y_completed_x1x2~x1*x2))$r.squared;



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

# performance
BD_mean_mse=mean((BD_mean_vec-true_mean)^2);
BD_mean_mse;

# lower and upper 95% CI;
BD_mean_low95_vec=BD_mean_vec-1.96*sqrt(BD_mean_var_vec);
BD_mean_up95_vec=BD_mean_vec+1.96*sqrt(BD_mean_var_vec);

BD_mean_length=mean(BD_mean_up95_vec-BD_mean_low95_vec);
BD_mean_length;

BD_mean_cov=(BD_mean_low95_vec < true_mean)*(BD_mean_up95_vec > true_mean);

mean(BD_mean_cov);

# variance estimates;
mean(BD_mean_var_vec);
var(BD_mean_vec);


# complete case analysis;

CC_mean_bias=mean(CC_mean_vec)-true_mean;
CC_mean_bias;

CC_mean_bias/true_mean;


CC_mean_mse=mean((CC_mean_vec-true_mean)^2);
CC_mean_mse;

# lower and upper 95% CI;
CC_mean_low95_vec=CC_mean_vec-1.96*sqrt(CC_mean_var_vec);
CC_mean_up95_vec=CC_mean_vec+1.96*sqrt(CC_mean_var_vec);

CC_mean_length=mean(CC_mean_up95_vec-CC_mean_low95_vec);
CC_mean_length;

CC_mean_cov=(CC_mean_low95_vec < true_mean)*(CC_mean_up95_vec > true_mean);

mean(CC_mean_cov);

# variance estimates;
mean(CC_mean_var_vec);
var(CC_mean_vec);


# multiple imputation estimates;
mean_mi_x1_mean=rep(NA, cycle_no);
mean_mi_x1_var=rep(NA, cycle_no);
mean_mi_x1_df=rep(NA, cycle_no);
mean_mi_x1_f=rep(NA, cycle_no);

# test=mi.inference(as.list(MI_mean_x1_mat[1,]), as.list(sqrt(MI_mean_x1_var_mat[i,])));


for (i in 1:cycle_no)
{
mean_summary=pool.scalar(MI_mean_x1_mat[i,], MI_mean_x1_var_mat[i,], n=rowobs);
mean_mi_x1_mean[i]=mean_summary$qbar;
mean_mi_x1_var[i]=mean_summary$t;
mean_mi_x1_df[i]=mean_summary$df;
mean_mi_x1_f[i]=mean_summary$fmi;


}


# For marginal means;
# proper imputation;

mean_mi_x1_bias=mean(mean_mi_x1_mean)-true_mean;
mean_mi_x1_bias;

mean_mi_x1_bias/true_mean;

mean_mi_x1_mse=mean((mean_mi_x1_mean-true_mean)^2);
mean_mi_x1_mse;

# coverage;
MI_x1_mean_low95_vec=mean_mi_x1_mean-qt(.975, mean_mi_x1_df)*sqrt(mean_mi_x1_var);
MI_x1_mean_up95_vec=mean_mi_x1_mean+qt(.975, mean_mi_x1_df)*sqrt(mean_mi_x1_var);

MI_x1_mean_length=mean(MI_x1_mean_up95_vec-MI_x1_mean_low95_vec);

MI_x1_mean_length;

MI_x1_mean_coverage=(MI_x1_mean_low95_vec < true_mean)*(MI_x1_mean_up95_vec > true_mean);
mean(MI_x1_mean_coverage);

# variance estimates;
mean(sqrt(mean_mi_x1_var));
sqrt(var(mean_mi_x1_mean));

mean_mi_x1_rsquare=rowMeans(MI_rsquare_x1_mat);

mean(mean_mi_x1_rsquare);

mean_mi_x1_f_approx=p_mis_vec*(1-mean_mi_x1_rsquare)/(1-p_mis_vec*mean_mi_x1_rsquare);

mean(mean_mi_x1_f);
mean(mean_mi_x1_f_approx);

plot(mean_mi_x1_f, mean_mi_x1_f_approx);


# multiple imputation estimates;
mean_mi_x1plusx2_mean=rep(NA, cycle_no);
mean_mi_x1plusx2_var=rep(NA, cycle_no);
mean_mi_x1plusx2_df=rep(NA, cycle_no);
mean_mi_x1plusx2_f=rep(NA, cycle_no);

# test=mi.inference(as.list(MI_mean_x1_mat[1,]), as.list(sqrt(MI_mean_x1_var_mat[i,])));


for (i in 1:cycle_no)
{
mean_summary=pool.scalar(MI_mean_x1plusx2_mat[i,], MI_mean_x1plusx2_var_mat[i,], n=rowobs);
mean_mi_x1plusx2_mean[i]=mean_summary$qbar;
mean_mi_x1plusx2_var[i]=mean_summary$t;
mean_mi_x1plusx2_df[i]=mean_summary$df;
mean_mi_x1plusx2_f[i]=mean_summary$fmi;


}


# For marginal means;
# proper imputation;

mean_mi_x1plusx2_bias=mean(mean_mi_x1plusx2_mean)-true_mean;
mean_mi_x1plusx2_bias;

mean_mi_x1plusx2_bias/true_mean;

mean_mi_x1plusx2_mse=mean((mean_mi_x1plusx2_mean-true_mean)^2);
mean_mi_x1plusx2_mse;

# coverage;
MI_x1plusx2_mean_low95_vec=mean_mi_x1plusx2_mean-qt(.975, mean_mi_x1plusx2_df)*sqrt(mean_mi_x1plusx2_var);
MI_x1plusx2_mean_up95_vec=mean_mi_x1plusx2_mean+qt(.975, mean_mi_x1plusx2_df)*sqrt(mean_mi_x1plusx2_var);

MI_x1plusx2_mean_length=mean(MI_x1plusx2_mean_up95_vec-MI_x1plusx2_mean_low95_vec);

MI_x1plusx2_mean_length;

MI_x1plusx2_mean_coverage=(MI_x1plusx2_mean_low95_vec < true_mean)*(MI_x1plusx2_mean_up95_vec > true_mean);
mean(MI_x1plusx2_mean_coverage);

# variance estimates;
mean(sqrt(mean_mi_x1plusx2_var));
sqrt(var(mean_mi_x1plusx2_mean));

mean_mi_x1plusx2_rsquare=rowMeans(MI_rsquare_x1plusx2_mat);

mean(mean_mi_x1plusx2_rsquare);


mean_mi_x1plusx2_f_approx=p_mis_vec*(1-mean_mi_x1plusx2_rsquare)/(1-p_mis_vec*mean_mi_x1plusx2_rsquare);

mean(mean_mi_x1plusx2_f);
mean(mean_mi_x1plusx2_f_approx);

plot(mean_mi_x1plusx2_f, mean_mi_x1plusx2_f_approx);


# multiple imputation estimates;
mean_mi_x1x2_mean=rep(NA, cycle_no);
mean_mi_x1x2_var=rep(NA, cycle_no);
mean_mi_x1x2_df=rep(NA, cycle_no);
mean_mi_x1x2_f=rep(NA, cycle_no);

# test=mi.inference(as.list(MI_mean_x1_mat[1,]), as.list(sqrt(MI_mean_x1_var_mat[i,])));


for (i in 1:cycle_no)
{
mean_summary=pool.scalar(MI_mean_x1x2_mat[i,], MI_mean_x1x2_var_mat[i,], n=rowobs);
mean_mi_x1x2_mean[i]=mean_summary$qbar;
mean_mi_x1x2_var[i]=mean_summary$t;
mean_mi_x1x2_df[i]=mean_summary$df;
mean_mi_x1x2_f[i]=mean_summary$fmi;


}


# For marginal means;
# proper imputation;

mean_mi_x1x2_bias=mean(mean_mi_x1x2_mean)-true_mean;
mean_mi_x1x2_bias;

mean_mi_x1x2_bias/true_mean;

mean_mi_x1x2_mse=mean((mean_mi_x1x2_mean-true_mean)^2);
mean_mi_x1x2_mse;

# coverage;
MI_x1x2_mean_low95_vec=mean_mi_x1x2_mean-qt(.975, mean_mi_x1x2_df)*sqrt(mean_mi_x1x2_var);
MI_x1x2_mean_up95_vec=mean_mi_x1x2_mean+qt(.975, mean_mi_x1x2_df)*sqrt(mean_mi_x1x2_var);

MI_x1x2_mean_length=mean(MI_x1x2_mean_up95_vec-MI_x1x2_mean_low95_vec);

MI_x1x2_mean_length;

MI_x1x2_mean_coverage=(MI_x1x2_mean_low95_vec < true_mean)*(MI_x1x2_mean_up95_vec > true_mean);
mean(MI_x1x2_mean_coverage);

# variance estimates;
mean(sqrt(mean_mi_x1x2_var));
sqrt(var(mean_mi_x1x2_mean));

mean_mi_x1x2_rsquare=rowMeans(MI_rsquare_x1x2_mat);

mean(mean_mi_x1x2_rsquare);


mean_mi_x1x2_f_approx=p_mis_vec*(1-mean_mi_x1x2_rsquare)/(1-p_mis_vec*mean_mi_x1x2_rsquare);

mean(mean_mi_x1x2_f);
mean(mean_mi_x1x2_f_approx);

plot(mean_mi_x1x2_f, mean_mi_x1x2_f_approx);



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
mean(BD_slope_var_vec);
var(BD_slope_vec);


# complete-case analysis;
CC_slope_bias=mean(CC_slope_vec)-true_slope;
CC_slope_bias;

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
mean(CC_slope_var_vec);
var(CC_slope_vec);


# multiple imputation estimates;

slope_mi_mean=rep(NA, cycle_no);
slope_mi_var=rep(NA, cycle_no);
slope_mi_df=rep(NA, cycle_no);
slope_mi_f=rep(NA, cycle_no);

slope_IM_mi_mean=rep(NA, cycle_no);
slope_IM_mi_var=rep(NA, cycle_no);
slope_IM_mi_df=rep(NA, cycle_no);
slope_IM_mi_f=rep(NA, cycle_no);



for (i in 1:cycle_no)
{
slope_summary=pool.scalar(MI_slope_mat[i,], MI_slope_var_mat[i,]);
slope_mi_mean[i]=slope_summary$qbar;
slope_mi_var[i]=slope_summary$t;
slope_mi_df[i]=slope_summary$df;
slope_mi_f[i]=slope_summary$f;

IM_slope_summary=pool.scalar(IM_MI_slope_mat[i,], IM_MI_slope_var_mat[i,]);
slope_IM_mi_mean[i]=IM_slope_summary$qbar;
slope_IM_mi_var[i]=IM_slope_summary$t;
slope_IM_mi_df[i]=IM_slope_summary$df;
slope_IM_mi_f[i]=IM_slope_summary$f;



}



# For slope;

# proper imputation;

slope_mi_bias=mean(slope_mi_mean)-true_slope;
slope_mi_bias;

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
mean(slope_mi_var);
var(slope_mi_mean);

# improper imputation;

slope_IM_mi_bias=mean(slope_IM_mi_mean)-true_slope;
slope_IM_mi_bias;

slope_IM_mi_mse=mean((slope_IM_mi_mean-true_slope)^2);
slope_IM_mi_mse;

# coverage;
IM_MI_slope_low95_vec=slope_IM_mi_mean-qt(.975, slope_IM_mi_df)*sqrt(slope_IM_mi_var);
IM_MI_slope_up95_vec=slope_IM_mi_mean+qt(.975, slope_IM_mi_df)*sqrt(slope_IM_mi_var);

IM_MI_slope_length=mean(IM_MI_slope_up95_vec-IM_MI_slope_low95_vec);
IM_MI_slope_length;


IM_MI_slope_coverage=(IM_MI_slope_low95_vec < true_slope)*(IM_MI_slope_up95_vec > true_slope);

mean(IM_MI_slope_coverage);

# variance estimates;
mean(slope_IM_mi_var);
var(slope_IM_mi_mean);

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
mean(BD_xyslope_var_vec);
var(BD_xyslope_vec);


# complete-case analysis;
CC_xyslope_bias=mean(CC_xyslope_vec)-true_xyslope;
CC_xyslope_bias;

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
mean(CC_xyslope_var_vec);
var(CC_xyslope_vec);


# multiple imputation estimates;

xyslope_mi_mean=rep(NA, cycle_no);
xyslope_mi_var=rep(NA, cycle_no);
xyslope_mi_df=rep(NA, cycle_no);
xyslope_mi_f=rep(NA, cycle_no);

xyslope_IM_mi_mean=rep(NA, cycle_no);
xyslope_IM_mi_var=rep(NA, cycle_no);
xyslope_IM_mi_df=rep(NA, cycle_no);
xyslope_IM_mi_f=rep(NA, cycle_no);



for (i in 1:cycle_no)
{
xyslope_summary=pool.scalar(MI_xyslope_mat[i,], MI_xyslope_var_mat[i,]);
xyslope_mi_mean[i]=xyslope_summary$qbar;
xyslope_mi_var[i]=xyslope_summary$t;
xyslope_mi_df[i]=xyslope_summary$df;
xyslope_mi_f[i]=xyslope_summary$f;

IM_xyslope_summary=pool.scalar(IM_MI_xyslope_mat[i,], IM_MI_xyslope_var_mat[i,]);
xyslope_IM_mi_mean[i]=IM_xyslope_summary$qbar;
xyslope_IM_mi_var[i]=IM_xyslope_summary$t;
xyslope_IM_mi_df[i]=IM_xyslope_summary$df;
xyslope_IM_mi_f[i]=IM_xyslope_summary$f;



}



# proper imputation;

xyslope_mi_bias=mean(xyslope_mi_mean)-true_xyslope;
xyslope_mi_bias;

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
mean(xyslope_mi_var);
var(xyslope_mi_mean);

# improper imputation;

xyslope_IM_mi_bias=mean(xyslope_IM_mi_mean)-true_xyslope;
xyslope_IM_mi_bias;

xyslope_IM_mi_mse=mean((xyslope_IM_mi_mean-true_xyslope)^2);
xyslope_IM_mi_mse;

# coverage;
IM_MI_xyslope_low95_vec=xyslope_IM_mi_mean-qt(.975, xyslope_IM_mi_df)*sqrt(xyslope_IM_mi_var);
IM_MI_xyslope_up95_vec=xyslope_IM_mi_mean+qt(.975, xyslope_IM_mi_df)*sqrt(xyslope_IM_mi_var);

IM_MI_xyslope_length=mean(IM_MI_xyslope_up95_vec-IM_MI_xyslope_low95_vec);
IM_MI_xyslope_length;


IM_MI_xyslope_coverage=(IM_MI_xyslope_low95_vec < true_xyslope)*(IM_MI_xyslope_up95_vec > true_xyslope);

mean(IM_MI_xyslope_coverage);

# variance estimates;
mean(xyslope_IM_mi_var);
var(xyslope_IM_mi_mean);


