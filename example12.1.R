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

# error variance;
var_error=1;


cycle_no=1;
mi_no=50;



for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# error distribution
# could be normal, or could be other distributions;
error=sqrt(var_error)*rnorm(rowobs);

# generate the latent variable;
y=beta0+beta1*x+error;

# set up the missing data as MCAR on x;
# miss_indi=runif(n=rowobs)<0.40;

# set up the missing data as MAR on x;
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


# now impute the missing y2's to get estimates of the marginal;
y_completed=y_completed_IM=y_miss;

for (i in 1:mi_no)
{
# proper model imputation;
y_imputed=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);


y_completed[miss_seq]=y_imputed;



}

cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

# plot the distribution of the data;

# MCAR ;
# hist(y_obs, freq=FALSE, main="MCAR");
# hist(y_imputed, freq=FALSE, xlab="y_imp", main="MCAR");

# plot(x_obs, y_obs, xlim=c(-2,5), main="MCAR");
# abline(lm(y_obs~x_obs));
# lm(y_obs ~ x_obs);

# plot(x_mis, y_imputed, xlim=c(-2,5), pch=17, col="red", ylab="y_imp", main="MCAR");
# abline(lm(y_imputed ~ x_mis), col="red");
# lm(y_imputed ~ x_mis);




# MAR;

hist(y_obs, freq=FALSE, main="MAR");
hist(y_imputed, freq=FALSE, xlab="y_imp", main="MAR");


plot(x_obs, y_obs, xlim=c(-2,5), main="MAR");
abline(lm(y_obs~x_obs));
lm(y_obs ~ x_obs);

plot(x_mis, y_imputed, xlim=c(-2,5), pch=17, col="red", ylab="y_imp", main="MAR");
abline(lm(y_imputed ~ x_mis), col="red");
lm(y_imputed ~ x_mis);





mean(y_obs);
var(y_obs);
mean(y_imputed);
var(y_imputed);

t.test(y_obs, y_imputed);
ks.test(y_obs, y_imputed);

