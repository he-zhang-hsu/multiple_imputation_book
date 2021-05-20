rm(list=ls());
library(mice);
library(norm);
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


rowobs=1000;

# set up the random seed;
set.seed(197789);

# generate data;
# meand and variance of the X covariate;
mu_x=1;
var_x=1;

# regression parameters; 
beta0=-2;
beta1=1;

# error variance;
var_error=1;


# number of simulations;
cycle_no=1000;

# matrices holding the parameter estimates;
# complete-data inferences;
# xyslope means regressing x on y;
# BD for before-deletion;
BD_mean_vec=BD_mean_var_vec=BD_slope_vec=BD_slope_var_vec=BD_xyslope_vec=BD_xyslope_var_vec=rep(NA, cycle_no);

# CC for complete-case analysis;
CC_mean_vec=CC_mean_var_vec=CC_slope_vec=CC_slope_var_vec=CC_xyslope_vec=CC_xyslope_var_vec=rep(NA, cycle_no);

# RP for regression prediction imputation;
RP_mean_vec=RP_mean_var_vec=RP_slope_vec=RP_slope_var_vec=RP_xyslope_vec=RP_xyslope_var_vec=rep(NA, cycle_no);

# INDI for missing data indicator method;
INDI_xyslope_vec=INDI_xyslope_var_vec=rep(NA, cycle_no);

for (cycle in 1:cycle_no)

{

set.seed(cycle);

# X is the covariate;
x=mu_x+sqrt(var_x)*rnorm(rowobs);

# normal error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the outcome variable Y;
y=beta0+beta1*x+error;

# set up the missing data as MCAR on x;
# These results are in Table 2.1;
# miss_indi=runif(n=rowobs)<0.40;

# set up the missing data as MAR on x;
# These results are in Table 2.2;
alpha0=-1.4;
alpha1=1;
miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x)/(1+exp(alpha0+alpha1*x)));


obs_indi=1-miss_indi;

y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
x_obs=x[miss_indi==0];
x_mis=x[miss_indi==1];

# compelete data analysis;
# marginal mean of Y;
BD_mean_vec[cycle]=mean(y);
BD_mean_var_vec[cycle]=var(y)/rowobs;


# regression analysis;
# regressing Y on X;
BD=summary(lm(y~x));

BD_slope_vec[cycle]=BD$coeff[2,1];
BD_slope_var_vec[cycle]=BD$coeff[2,2]^2;

# reverse regression analysis;
# regressing X on Y;
BD_xy=summary(lm(x~y));

BD_xyslope_vec[cycle]=BD_xy$coeff[2,1];
BD_xyslope_var_vec[cycle]=BD_xy$coeff[2,2]^2;


# different missing data methods;

# complete-case analysis;

# marginal mean;
CC_mean_vec[cycle]=mean(y_obs);
CC_mean_var_vec[cycle]=var(y_obs)/obs_no;


# regression analysis;
# regressing Y on X;
CC=summary(lm(y_obs~x_obs));

CC_slope_vec[cycle]=CC$coeff[2,1];
CC_slope_var_vec[cycle]=CC$coeff[2,2]^2;

# reverse regression analysis;
# regressing X on Y;
CC_xy=summary(lm(x_obs~y_obs));

CC_xyslope_vec[cycle]=CC_xy$coeff[2,1];
CC_xyslope_var_vec[cycle]=CC_xy$coeff[2,2]^2;

# missing data indicator method;
y_miss_0=y_miss;
y_miss_0[miss_indi==1]=0;

INDI_xy=summary(lm(x ~ obs_indi+y_miss_0));
INDI_xyslope_vec[cycle]=INDI_xy$coeff[3,1];
INDI_xyslope_var_vec[cycle]=INDI_xy$coeff[3,2]^2;

# now impute the missing Y's to get estimates of the marginal;
y_completed=y_miss;

# single imputation method using regression prediction;
# using R mice;
y_imputed=mice.impute.norm.predict(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);


y_completed[miss_seq]=y_imputed;

# completed-data analysis;

# marginal mean of Y;
RP_mean_vec[cycle]=mean(y_completed);
RP_mean_var_vec[cycle]=var(y_completed)/rowobs;

# regression analysis;
# regressing Y on X;
IMP=summary(lm(y_completed~x));

RP_slope_vec[cycle]=IMP$coeff[2,1];
RP_slope_var_vec[cycle]=IMP$coeff[2,2]^2;

# reverse regression analysis;
# regressing X on Y;
IMP_xy=summary(lm(x~y_completed));

RP_xyslope_vec[cycle]=IMP_xy$coeff[2,1];
RP_xyslope_var_vec[cycle]=IMP_xy$coeff[2,2]^2;


cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;


#############################################################################
# Simulation estimates;

# marginal mean of Y;
true_mean=mean(BD_mean_vec);
true_mean;

# before-deletion;
BD=performance(BD_mean_vec, BD_mean_var_vec, true_mean);
BD;

# complete-case analysis;
CC=performance(CC_mean_vec, CC_mean_var_vec, true_mean);
CC;

# RP imputation;
RP=performance(RP_mean_vec, RP_mean_var_vec, true_mean);
RP;


##########################################################################
# slope estimand for regessing Y on X;

true_slope=mean(BD_slope_vec);
true_slope;

# before-deletion;
BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;

# complete-case analysis;
CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;

# RP imputation;
RP=performance(RP_slope_vec, RP_slope_var_vec, true_slope);
RP;



#######################################################################
# slope estimand for regressing X on Y;

true_xyslope=mean(BD_xyslope_vec);
true_xyslope;

# before-deletion; 
BD=performance(BD_xyslope_vec, BD_xyslope_var_vec, true_xyslope);
BD;

# complete-case analysis;

CC=performance(CC_xyslope_vec, CC_xyslope_var_vec, true_xyslope);
CC;

# RP imputation;
RP=performance(RP_xyslope_vec, RP_xyslope_var_vec, true_xyslope);
RP;

# missingness indicator method;
INDI=performance(INDI_xyslope_vec, INDI_xyslope_var_vec, true_xyslope);
INDI;


###################################
# plot for Fig. 2.2 using the data from last simulation replicate;
# using only 250 data points;

y_plot_obs=y;
y_plot_miss=y;
y_plot_obs[miss_indi==1]=NA;
y_plot_miss[miss_indi==0]=NA;

y_plot=cbind(y_plot_obs, y_plot_miss);

# plot the observed and missing points;
# scatter plot
matplot(x[1:250], y_plot[1:250,], pch=c(1,16), main="Scatter plot before deletion",
 xlab="x", ylab="y");
# matplot(y_plot[1:250,], x[1:250], pch=c(1,16), xlab="y", ylab="x");

# histogram
hist(y_plot_miss, xlab="y_miss", main="Histogram of deleted values");

# plot the observed and imputed points from RP;

y_plot_imputed=y_completed;
y_plot_imputed[miss_indi==0]=NA;
y_plot_impute=cbind(y_plot_obs, y_plot_imputed);


# scatter plot;
matplot(x[1:250], y_plot_impute[1:250,], pch=c(1,17), xlab="x", ylab="y",
   main="Scatter plot after regression prediction");

# matplot(y_plot_impute[1:250,], x[1:250], pch=c(1,17), xlab="y", ylab="x");

# histogram;
hist(y_plot_imputed, xlab="y_predicted", main="Histogram of predicted values");
