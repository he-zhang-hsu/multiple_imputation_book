# Example 4.4

rm(list=ls());

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


# function bootstrap() provided the sample index of a nonparameteric sample;
bootstrap=function(n)
{
 rowobs=n;  
# donor_vector=ceiling(runif(rowobs)*rowobs);
donor_vector=sample(seq(1,n,1), replace=TRUE); 
donor_vector;
}

# function mice.impute.lda.proper() imputes the binary data
# using discriminant analysis model and is from R MICE

mice.impute.lda.proper <- function(y, ry, x, wy = NULL, ...) { 
if (is.null(wy)) wy <- !ry 
fy <- as.factor(y) 
nc <- length(levels(fy)) 
    
   #   SvB June 2009 - take bootstrap sample of training data 
      idx <- sample((1:length(y))[ry], size=sum(ry), replace=TRUE) 
      x[ry,] <- x[idx,] 
      y[ry] <- y[idx] 
   #   end bootstrap 
fy <- as.factor(y) 
nc <- length(levels(fy)) 

   fit <- lda(x, fy, subset = ry) 
   post <- predict(fit, x[wy, , drop = FALSE])$posterior 
   un <- rep(runif(sum(wy)), each = nc) 
   idx <- 1 + apply(un > apply(post, 1, cumsum), 2, sum) 
   return(levels(fy)[idx]) 
 } 


# the complete-data sample size in the simulation;

rowobs=1000;

# set up the random seed;
set.seed(197789);

# number of simulations;
cycle_no=1000;

# number of multiple imputations;

mi_no=50;

# generate complete data using
# discriminant analysis model;
p=0.5;
# p=0.1;

# c=1;
# c=2;
# c=3;
c=4;

# Vectors and matrices holding the parameter estimates;
# For marginal means;
# Before-deletion analysis;
BD_mean_vec=BD_var_vec=rep(NA, cycle_no);

# Complete-case analysis;
CC_mean_vec=CC_var_vec=rep(NA, cycle_no);

# For multiple imputation;
# Use MLE logistic regression imputation;
mle_mean_mat=mle_var_mat=matrix(NA, cycle_no, mi_no);

# Use discriminant analysis imputation;
lda_mean_mat=lda_var_mat=matrix(NA, cycle_no, mi_no);

# For logistic regression slope coefficient;
# Before-deletion analysis;
BD_slope_vec=BD_slope_var_vec=rep(NA, cycle_no);

# Complete-case analysis;
CC_slope_vec=CC_slope_var_vec=rep(NA, cycle_no);

# For multiple imputation;
# Use MLE logistic regression imputation;
mle_slope_mat=mle_slope_var_mat=matrix(NA, cycle_no, mi_no);

# Use discriminant analysis imputation;
lda_slope_mat=lda_slope_var_mat=matrix(NA, cycle_no, mi_no);


# Begin the simulation;

for (cycle in 1:cycle_no)

{
set.seed(cycle);

# generate the model using the discriminant analysis model;

y=rbinom(n=rowobs, size=1, prob=p);

error=rnorm(rowobs);

# generate the discriminant x;
x=1+c*y+error;

# 40% missing;
# MCAR;
miss_indi=runif(n=rowobs)<0.40;

y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
x_obs=x[miss_indi==0];
x_miss=x[miss_indi==1];

# observed data is cbind(x, y_miss);

# compelete data analysis;
# sample mean;
BD_mean=mean(y);
BD_var=var(y)/rowobs;
BD_mean_vec[cycle]=BD_mean;
BD_var_vec[cycle]=BD_var;

# logistic regression coefficient for slope;
BD_logistic=summary(glm(y~x, family=binomial(link="logit")));
BD_slope_vec[cycle]=BD_logistic$coeff[2,1];
BD_slope_var_vec[cycle]=BD_logistic$coeff[2,2]^2;

# different missing data methods;

# complete-case analysis;
# mean estimands;
CC_mean=mean(y_miss, na.rm="T");
CC_mean_vec[cycle]=CC_mean;
CC_var_vec[cycle]=var(y_miss, na.rm="T")/obs_no;

# logistic regression coefficient for slope;
CC_logistic=summary(glm(y_obs~x_obs, family=binomial(link="logit")));
CC_slope_vec[cycle]=CC_logistic$coeff[2,1];
CC_slope_var_vec[cycle]=CC_logistic$coeff[2,2]^2;

# now impute the missing y's to get estimates of the marginal;
y_completed_mle=y_completed_lda=y_miss;

for (i in 1:mi_no)
{
set.seed(i);
 
# MLE logistic regression imputation;
y_imputed_mle=mice.impute.logreg(y_miss, ry=as.logical(1-miss_indi), seed=i, x=as.matrix(x));

# discriminant imputation based;

y_imputed_lda=mice.impute.lda.proper(y_miss, ry=as.logical(1-miss_indi), seed=i, x=as.matrix(x));

y_completed_mle[miss_seq]=y_imputed_mle;
y_completed_lda[miss_seq]=as.numeric(y_imputed_lda);


# For MLE logistic regression imputation;
# marginal means;
mle_mean=mean(y_completed_mle);
mle_var=var(y_completed_mle)/rowobs;

mle_mean_mat[cycle,i]=mle_mean;
mle_var_mat[cycle,i]=mle_var;

# logistic regression coefficient for slope;
mle_logistic=summary(glm(y_completed_mle~x, family=binomial(link="logit")));
mle_slope_mat[cycle,i]=mle_logistic$coeff[2,1];
mle_slope_var_mat[cycle,i]=mle_logistic$coeff[2,2]^2;


# For discriminant analysis imputation;
lda_mean=mean(y_completed_lda);
lda_mean_var=var(y_completed_lda)/rowobs;

lda_mean_mat[cycle,i]=lda_mean;
lda_var_mat[cycle,i]=lda_mean_var;

# logistic regression coefficient for slope;
lda_logistic=summary(glm(y_completed_lda~x, family=binomial(link="logit")));
lda_slope_mat[cycle,i]=lda_logistic$coeff[2,1];
lda_slope_var_mat[cycle,i]=lda_logistic$coeff[2,2]^2;



}

cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

#############################################################################
# Simulation results;
# marginal mean estimand;
true_mean=mean(BD_mean_vec);
true_mean;

# before-deletion;
BD=performance(BD_mean_vec, BD_var_vec, true_mean);
BD;

# complete-case analysis;
CC=performance(CC_mean_vec, CC_var_vec, true_mean);
CC;

# Multiple imputation analysis;
# MLE logistic regression imputation;

mle_MI=mi_performance(mle_mean_mat, mle_var_mat, true_mean, rowobs);
mle_MI;

# discriminant analysis imputation;
lda_MI=mi_performance(lda_mean_mat, lda_var_mat, true_mean, rowobs);
lda_MI;



##########################################################################
# slope estimand;
# Table 4.6
# complete-data analysis;

true_slope=mean(BD_slope_vec);
true_slope;

# before-deletion;
BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;

# complete-case analysis;
CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;

# Multiple imputation analysis;
# MLE logistic regression imputation;

mle_MI=mi_performance(mle_slope_mat, mle_slope_var_mat, true_slope, rowobs);
mle_MI;

# discriminant analysis imputation;
lda_MI=mi_performance(lda_slope_mat, lda_slope_var_mat, true_slope, rowobs);
lda_MI;






