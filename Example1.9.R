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


x=mu_x+sqrt(var_x)*rnorm(rowobs);

# error distribution
# could be normal, or could be other distributions;
error=sqrt(var_error)*rnorm(rowobs);

# generate the latent variable;
y=beta0+beta1*x+error;

# set up the missing data as MCAR on x;
# miss_indi=runif(n=rowobs)<0.40;

# before deletion;
bd_mean=mean(y);

bd_mean;

bd_se=sd(y)/sqrt(rowobs);

bd_se;

bd_fit=lm(y~x);

plot(x,y, main="Before Deletion");
abline(bd_fit);


# set up the missing data as MCAR on y;
set.seed(197789);

miss_indi=rbinom(n=rowobs, size=1, prob=0.4);

y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
x_obs=x[miss_indi==0];
x_mis=x[miss_indi==1];

# MCAR;

mcar_mean=mean(y_obs, na.rm=T);

mcar_mean;

mcar_se=sd(y_obs)/sqrt(obs_no);

mcar_se;

mcar_fit=lm(y_miss~x);

plot(x,y_miss, ylab="y_mcar", main="MCAR");
abline(mcar_fit);

mcar_logistic=summary(glm(miss_indi~x, family=binomial(link="logit")));

mcar_logistic;

y_mcar=y_miss;

# MAR;
set.seed(197789);

alpha0=-1.4;
alpha1=1;
miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x)/(1+exp(alpha0+alpha1*x)));


y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
x_obs=x[miss_indi==0];
x_mis=x[miss_indi==1];

mar_mean=mean(y_obs, na.rm=T);

mar_mean;

mar_se=sd(y_obs)/sqrt(obs_no);

mar_se;

mar_fit=lm(y_miss~x);

plot(x,y_miss, ylab="y_mar", main="MAR");
abline(mar_fit);

mar_logistic=summary(glm(miss_indi~x, family=binomial(link="logit")));

mar_logistic;



y_mar=y_miss;

# MNAR;
set.seed(197789);

alpha0=0.2;
alpha1=0.2;
alpha2=1;
miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x+alpha2*y)/(1+exp(alpha0+alpha1*x+alpha2*y)));


y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
x_obs=x[miss_indi==0];
x_mis=x[miss_indi==1];
mis_no;

mnar_mean=mean(y_obs, na.rm=T);

mnar_mean;

mnar_se=sd(y_obs)/sqrt(obs_no);

mnar_se;

mnar_fit=lm(y_miss~x);

mnar_logistic=summary(glm(miss_indi~x, family=binomial(link="logit")));

mnar_logistic;



plot(x,y_miss, ylab="y_mnar", main="MNAR");
abline(mnar_fit);



y_mnar=y_miss;

# box plots for all the graphs;

y_data=cbind(y, y_mcar, y_mar, y_mnar);

boxplot(y_data);

# compelete data analysis;
# marginal mean;
BD_mean_vec[cycle]=mean(y);
BD_mean_var_vec[cycle]=var(y)/rowobs;


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
CC_mean_vec[cycle]=mean(y_obs);
CC_mean_var_vec[cycle]=var(y_obs)/obs_no;


# regression analysis;
CC=summary(lm(y_obs~x_obs));

CC_slope_vec[cycle]=CC$coeff[2,1];
CC_slope_var_vec[cycle]=CC$coeff[2,2]^2;

# reverse regression analysis;
CC_xy=summary(lm(x_obs~y_obs));

CC_xyslope_vec[cycle]=CC_xy$coeff[2,1];
CC_xyslope_var_vec[cycle]=CC_xy$coeff[2,2]^2;



# now impute the missing y2's to get estimates of the marginal;
y_completed=y_completed_IM=y_miss;

for (i in 1:mi_no)
{
# proper model imputation;
y_imputed=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);


y_completed[miss_seq]=y_imputed;

# completed-data analysis;

# marginal mean;
MI_mean_mat[cycle,i]=mean(y_completed);
MI_mean_var_mat[cycle,i]=var(y_completed)/rowobs;

# regression analysis;

IMP=summary(lm(y_completed~x));

MI_slope_mat[cycle,i]=IMP$coeff[2,1];
MI_slope_var_mat[cycle,i]=IMP$coeff[2,2]^2;

# reverse regression analysis;

IMP_xy=summary(lm(x~y_completed));

MI_xyslope_mat[cycle,i]=IMP_xy$coeff[2,1];
MI_xyslope_var_mat[cycle,i]=IMP_xy$coeff[2,2]^2;

y_imputed_IM=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=rep(1,rowobs));


y_completed_IM[miss_seq]=y_imputed_IM;

# completed-data analysis;

# marginal mean;
IM_MI_mean_mat[cycle,i]=mean(y_completed_IM);
IM_MI_mean_var_mat[cycle,i]=var(y_completed_IM)/rowobs;

# regression analysis;

IM_IMP=summary(lm(y_completed_IM~x));

IM_MI_slope_mat[cycle,i]=IM_IMP$coeff[2,1];
IM_MI_slope_var_mat[cycle,i]=IM_IMP$coeff[2,2]^2;

# reverse regression analysis;

IM_IMP_xy=summary(lm(x~y_completed_IM));

IM_MI_xyslope_mat[cycle,i]=IM_IMP_xy$coeff[2,1];
IM_MI_xyslope_var_mat[cycle,i]=IM_IMP_xy$coeff[2,2]^2;



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
mean_mi_mean=rep(NA, cycle_no);
mean_mi_var=rep(NA, cycle_no);
mean_mi_df=rep(NA, cycle_no);
mean_mi_f=rep(NA, cycle_no);

mean_IM_mi_mean=rep(NA, cycle_no);
mean_IM_mi_var=rep(NA, cycle_no);
mean_IM_mi_df=rep(NA, cycle_no);
mean_IM_mi_f=rep(NA, cycle_no);



for (i in 1:cycle_no)
{
mean_summary=pool.scalar(MI_mean_mat[i,], MI_mean_var_mat[i,]);
mean_mi_mean[i]=mean_summary$qbar;
mean_mi_var[i]=mean_summary$t;
mean_mi_df[i]=mean_summary$df;
mean_mi_f[i]=mean_summary$f;

IM_mean_summary=pool.scalar(IM_MI_mean_mat[i,], IM_MI_mean_var_mat[i,]);
mean_IM_mi_mean[i]=IM_mean_summary$qbar;
mean_IM_mi_var[i]=IM_mean_summary$t;
mean_IM_mi_df[i]=IM_mean_summary$df;
mean_IM_mi_f[i]=IM_mean_summary$f;


}


# For marginal means;
# proper imputation;

mean_mi_bias=mean(mean_mi_mean)-true_mean;
mean_mi_bias;

mean_mi_mse=mean((mean_mi_mean-true_mean)^2);
mean_mi_mse;

# coverage;
MI_mean_low95_vec=mean_mi_mean-qt(.975, mean_mi_df)*sqrt(mean_mi_var);
MI_mean_up95_vec=mean_mi_mean+qt(.975, mean_mi_df)*sqrt(mean_mi_var);

MI_mean_length=mean(MI_mean_up95_vec-MI_mean_low95_vec);

MI_mean_length;

MI_mean_coverage=(MI_mean_low95_vec < true_mean)*(MI_mean_up95_vec > true_mean);
mean(MI_mean_coverage);

# variance estimates;
mean(mean_mi_var);
var(mean_mi_mean);

# improper imputation;

mean_IM_mi_bias=mean(mean_IM_mi_mean)-true_mean;
mean_IM_mi_bias;

mean_IM_mi_mse=mean((mean_IM_mi_mean-true_mean)^2);
mean_IM_mi_mse;

# coverage;
IM_MI_mean_low95_vec=mean_IM_mi_mean-qt(.975, mean_IM_mi_df)*sqrt(mean_IM_mi_var);
IM_MI_mean_up95_vec=mean_IM_mi_mean+qt(.975, mean_IM_mi_df)*sqrt(mean_IM_mi_var);

IM_MI_mean_length=mean(IM_MI_mean_up95_vec-IM_MI_mean_low95_vec);

IM_MI_mean_length;

IM_MI_mean_coverage=(IM_MI_mean_low95_vec < true_mean)*(IM_MI_mean_up95_vec > true_mean);
mean(IM_MI_mean_coverage);

# variance estimates;
mean(mean_IM_mi_var);
var(mean_IM_mi_mean);


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


