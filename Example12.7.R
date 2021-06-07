# Example 12.7

rm(list=ls());
library(mice);
library(norm);
# library(HI);
options(digits=4);

# number of complete-data sample size;
rowobs=1000;

# number of simulations;
cycle_no=1000;

# mean and variance of y;
mu_y=1;
var_y=1;


# generate data;
set.seed(197789);

y_matrix=matrix(rnorm(rowobs*cycle_no,mean=mu_y, sd=sqrt(var_y)), rowobs, cycle_no);
ymiss_matrix=y_matrix;
ymiss_matrix[1:rowobs/2,]=NA;


# number of multiple imputations;
mi_no_vector=c(2,5,10,20,50,100,200,1000);

experiment_no=length(mi_no_vector);

mean_estimate_mat=fmi_estimate_mat=matrix(NA, cycle_no, experiment_no);


# set.seed(197789);

for (l in 1:experiment_no)
{

set.seed(2*l);

mi_no=mi_no_vector[l];

MI_mean_mat=MI_mean_var_mat=matrix(NA, cycle_no, mi_no);


for (cycle in 1:cycle_no)

{


y_miss=ymiss_matrix[,cycle];

miss_indi=is.na(y_miss);

miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y_miss[miss_indi==0];

# generate imputed data;

# obtain observed-data sufficient statistics;
mu_hat=mean(y_obs);
s_square=1/(obs_no-1)*crossprod(y_obs-mu_hat);

# multiple imputation;
y_completed=y_miss;

for (i in 1:mi_no)
{

# set.seed(2*i);

# draw sigma_square;
sigma_square=1/(rgamma(1, shape=obs_no/2-1/2, scale=2/((obs_no-1)*s_square)));
mu=rnorm(1, mu_hat, sd=sqrt(sigma_square/obs_no));

y_imputed=rnorm(mis_no, mu, sd=sqrt(sigma_square));


y_completed[miss_seq]=y_imputed;


# completed-data analysis;

# marginal mean;
MI_mean_mat[cycle,i]=mean(y_completed);
MI_mean_var_mat[cycle,i]=var(y_completed)/rowobs;

}

cat("the cycle is", cycle, "\n");



}

# multiple imputation estimates;

mean_mi_mean=rep(NA, cycle_no);
mean_mi_var=rep(NA, cycle_no);
mean_mi_df=rep(NA, cycle_no);
mean_mi_f=rep(NA, cycle_no);
mean_mi_bw=rep(NA, cycle_no);


for (k in 1:cycle_no)
{
mean_summary=pool.scalar(MI_mean_mat[k,], MI_mean_var_mat[k,], n=rowobs);
mean_mi_mean[k]=mean_summary$qbar;
mean_mi_var[k]=mean_summary$t;
mean_mi_df[k]=mean_summary$df;
mean_mi_f[k]=mean_summary$f;
mean_mi_bw[k]=mean_summary$b;

}


# For marginal means;

mean_estimate_mat[,l]=mean_mi_mean;

fmi_estimate_mat[,l]=mean_mi_f;

cat("the number of imputations is", mi_no, "\n");


}

# colMeans(fmi_estimate_mat);
# apply(fmi_estimate_mat, 2, sd);
# plot(log(mi_no_vector),colMeans(fmi_estimate_mat));


# Fig. 12.4

boxplot(fmi_estimate_mat, names=mi_no_vector, xlab="number of imputations", ylab="FMI estimates");
abline(h=0.5, col="red");

# colMeans(mean_estimate_mat);

# boxplot(mean_estimate_mat);

# apply(mean_estimate_mat, 2, sd);
