# Example 12.12

rm(list=ls());
# library(car);
# library(mvtnorm);
library(mice);
library(norm);
# library(HI);
options(digits=4);

# set up the random seed;
set.seed(197789);

# Complete-data sample size;
rowobs=1000;

# Meand and variance of x;
mu_x=1;
var_x=1;

# regression parameters;
# Scenario I;
# regression parameters 
beta0=-2;
beta1=1;

# error variance;
var_error=1;

# number of simulations;
cycle_no=1000;

# number of multiple imputations;
mi_no=50;

# matrices holding the predictions or imputations;
SI_pred_vec=rep(NA, cycle_no);

MI_pred_mat=matrix(NA, cycle_no, mi_no);

# begin the simulation;

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# normal error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the y variable;
# Scenario I;
y=beta0+beta1*x+error;


# set up the missing data as MCAR on x;
miss_indi=runif(n=rowobs)<0.40;

# set up the missing data as MAR on x;
# alpha0=-1.4;
# alpha1=1;
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x)/(1+exp(alpha0+alpha1*x)));


obs_indi=1-miss_indi;

y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
x_obs=x[miss_indi==0];
x_mis=x[miss_indi==1];

# now impute the missing y's to get estimates of the marginal;
y_completed=y_miss;

# single imputation method;
y_imputed=mice.impute.norm.predict(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);


y_completed[miss_seq]=y_imputed;

# calculate the prediction accruacy;
SI_pred_vec[cycle]=mean((y_imputed-y[miss_indi==1])^2);

# multiple imputation analysis;

# now impute the missing y's to get estimates of the marginal;
y_completed_MI=y_miss;

# the matrix holds the multiply imputed data;

y_completed_MI_mat=matrix(NA, nrow=rowobs, ncol=mi_no);

for (i in 1:mi_no)
{

set.seed(i);

y_imputed_MI=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);


y_completed_MI[miss_seq]=y_imputed_MI;


y_completed_MI_mat[,i]=y_completed_MI;

# prediction accuracy for multiple imputations;
MI_pred_mat[cycle,i]=mean((y_imputed_MI-y[miss_indi==1])^2);

}


cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

###########################################################################
# Check the prediction accuracy
# single imputation;
mean(SI_pred_vec);

# multiple imputation;
mean(MI_pred_mat);



