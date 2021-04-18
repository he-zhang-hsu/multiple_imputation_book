################## data generate ##############
# Program: boxcox_impute_expolore_simulation; ############
# generate data, and try out different missing data methods;
# run some simulations to assess the biasedness of various methods;

# library(car);
library(mice);
library(norm);
# adaptive metroplis rejection sampling;
library(HI);
# include the library MASS for generating from multivariate normal distribution;
library(MASS);
# include the library msm and bayesm for using the function of drawing truncated normal
# distribution;
library(msm);
library(bayesm);
library(mvtnorm);
library(sn);
library(arm);
# library(MCMCpack);


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

bootstrap=function(n)
{
 rowobs=n;  
 donor_vector=ceiling(runif(rowobs)*rowobs);
 donor_vector;
}

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
estimand=as.data.frame(cbind(bias,rbias,mse,mean_length,mean_cov,sd,se))
}








rowobs=1000;

# set up the random seed;
set.seed(197789);

# generate data;
# distribution of x, the covariate;
# could be normal, or could be other distributions;
# fixed at mean 10, and variance 1;
mu_x=1;
var_x=1/4;

# regression parameters;
# both beta0 and beta1 are fixed at 
beta0=-2;
# beta1=1; # beta1=1 creates around 30% of 1's;
beta1=2; # beta1=2 creates around 50% of 1's;
# beta1=0.5;
#beta1=0; # creates around 12% of 1's;
# beta1=-1; # beta1=-1 creates around 7% of 1's;
# beta1=-2; # beta1=-2 creates around 6.7% of 1's;
# beta1=-3; # beta1=-3 creates around 7% of 1's;
# beta1=-4; # beta1=-4 creates around 1.57% of 1's;
# beta1=-5; # beta1=-5 creates 1.25% of 1's;
# beta1=-6;



cycle_no=1000;
mi_no=80;
sens=0.9;
spec=0.8;

# matrices holding the parameter estimates;
# complete-data inferences;
BD_mean_vec=BD_mean_var_vec=rep(NA, cycle_no);
CC_mean_vec=CC_mean_var_vec=rep(NA, cycle_no);
SI_mean_vec=SI_mean_var_vec=rep(NA, cycle_no);
mle_mean_mat=mle_mean_var_mat=matrix(NA, cycle_no, mi_no);

BD_slope_vec=BD_slope_var_vec=rep(NA, cycle_no);
CC_slope_vec=CC_slope_var_vec=rep(NA, cycle_no);
SI_slope_vec=SI_slope_var_vec=rep(NA, cycle_no);
mle_slope_mat=mle_slope_var_mat=matrix(NA, cycle_no, mi_no);




for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);
# x=runif(rowobs, min=0, max=2);

# generate the logistic regression outcome y2;
y2_true=rbinom(n=rowobs, size=1, prob=exp(beta0+beta1*x)/(1+exp(beta0+beta1*x)));

# assign the misclassification;
y2=rep(NA, rowobs);
alpha0=log(1/spec-1);
alpha1=log(1/(1-sens)-1)-log(1/spec-1);
y2=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*y2_true)/(1+exp(alpha0+alpha1*y2_true)));

# set up the missing data as MCAR on x;

miss_indi=runif(n=rowobs)<0.80;
y2_miss=y2_true;
y2_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y2_miss[!is.na(y2_miss)]);
mis_no=rowobs-obs_no;
y2_obs=y2_true[miss_indi==0];
x_obs=x[miss_indi==0];
x_miss=x[miss_indi==1];

# observed data is cbind(x, y2_miss);

# compelete data analysis;
# sample mean;
y2_mean=mean(y2_true);
y2_var=var(y2_true)/rowobs;
BD_mean_vec[cycle]=y2_mean;
BD_mean_var_vec[cycle]=y2_var;

# logistic regression coefficient for slope;
BD_logistic=summary(glm(y2_true~x, family=binomial(link="logit")));
BD_slope_vec[cycle]=BD_logistic$coeff[2,1];
BD_slope_var_vec[cycle]=BD_logistic$coeff[2,2]^2;


# different missing data methods;

# complete-case analysis;
# use the validation sample
# mean estimands;
cc_mean=mean(y2_miss, na.rm="T");
CC_mean_vec[cycle]=cc_mean;
CC_mean_var_vec[cycle]=var(y2_miss, na.rm="T")/obs_no;

# logistic regression coefficient for slope;
CC_logistic=summary(glm(y2_obs~x_obs, family=binomial(link="logit")));
CC_slope_vec[cycle]=CC_logistic$coeff[2,1];
CC_slope_var_vec[cycle]=CC_logistic$coeff[2,2]^2;

# use the misreported version;
SI_mean_vec[cycle]=mean(y2);
SI_mean_var_vec[cycle]=var(y2)/rowobs;

# logistic regression coefficient for slope;
SI_logistic=summary(glm(y2~x, family=binomial(link="logit")));
SI_slope_vec[cycle]=SI_logistic$coeff[2,1];
SI_slope_var_vec[cycle]=SI_logistic$coeff[2,2]^2;




# now impute the missing y2's to get estimates of the marginal;
y2_completed_mle=y2_miss;
x_new=cbind(x, y2);

for (i in 1:mi_no)
{

set.seed(i);

# bivariate model imputation; 
# Rubin's imputation based on MLE;
y2_imputed_mle=mice.impute.logreg(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=x_new);

y2_completed_mle[miss_seq]=y2_imputed_mle;

# marginal means;
mle_mean=mean(y2_completed_mle);
mle_var=var(y2_completed_mle)/rowobs;

mle_mean_mat[cycle,i]=mle_mean;
mle_mean_var_mat[cycle,i]=mle_var;

# logistic regression coefficient for slope;
mle_logistic=summary(glm(y2_completed_mle~x, family=binomial(link="logit")));
mle_slope_mat[cycle,i]=mle_logistic$coeff[2,1];
mle_slope_var_mat[cycle,i]=mle_logistic$coeff[2,2]^2;




}

cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

# plot the distribution of the data;

# check the performance of the estimates;
# population quantity
# mean estimand;
# true=pnorm((beta0+beta1*mu_x-cut_off)/(beta1^2*var_x+var_error));

#############################################################################
# marginal mean estimand;

true=mean(BD_mean_vec);
true;

true;


BD=performance(BD_mean_vec, BD_mean_var_vec, true);
BD;
CC=performance(CC_mean_vec, CC_mean_var_vec, true);
CC;
SI=performance(SI_mean_vec, SI_mean_var_vec, true);
SI;

MI=mi_performance(mle_mean_mat, mle_mean_var_mat, true, rowobs);
MI;





##########################################################################
# slope estimand;

# complete-data analysis;

true_slope=mean(BD_slope_vec);
true_slope;

BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;
CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;
SI=performance(SI_slope_vec, SI_slope_var_vec, true_slope);
SI;

MI=mi_performance(mle_slope_mat, mle_slope_var_mat, true_slope, rowobs);
MI;


