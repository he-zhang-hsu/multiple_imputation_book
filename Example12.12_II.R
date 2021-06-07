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

# complete-data sample size;
rowobs=1000;

# mean and variance of x;
mu_x=1;
var_x=1;

# regression parameters; 
beta0=-2;
beta1=1;
beta2=1;

# error variance;
var_error=1;

# number of simulations;
cycle_no=1000;

# number of multiple imputations;
mi_no=50;

# For prediction error;
# MI for linear model imputation;
MI_pred_mat=matrix(NA, cycle_no, mi_no);


# IM_MI for quadratic model imputation;
IM_pred_mat=matrix(NA, cycle_no, mi_no);


# begin the simulation;

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# normal error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the y variable;
y=beta0+beta1*x+beta2*x^2+error;

# set up the missing data as MCAR on x;
miss_indi=runif(n=rowobs)<0.40;

# set up the missing data as MAR on x;
# alpha0=-2;
# alpha1=1;
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x)/(1+exp(alpha0+alpha1*x)));


y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
x_obs=x[miss_indi==0];
x_mis=x[miss_indi==1];

# now impute the missing y's to get estimates of the marginal;
y_completed=y_completed_IM=y_miss;

for (i in 1:mi_no)
{
set.seed(i);
# linear normal model imputation;
y_imputed=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);


y_completed[miss_seq]=y_imputed;

# prediction error;

MI_pred_mat[cycle,i]=mean((y_imputed-y[miss_indi==1])^2);



# imputation method including both x and x^2;
y_imputed_IM=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x,x^2));

y_completed_IM[miss_seq]=y_imputed_IM;

# prediction error;

IM_pred_mat[cycle,i]=mean((y_imputed_IM-y[miss_indi==1])^2);




}

cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;



# prediction error;

mean(MI_pred_mat);
mean(IM_pred_mat);
