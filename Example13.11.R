# Example 13.11

rm(list=ls());
library(sampleSelection)
library(R2WinBUGS);
library(mice);
library(norm);
library(MASS);
options(digits=4);

# function parameter_draw() draws parameters from a Bayesian linear regression model;

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



# read the data;

a1c=read.csv(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter13_univariate_nonignorable\\program\\A1c.csv");

rowobs=nrow(a1c);

attach(a1c);

# BMI has two missing values;
sum(is.na(male));
sum(is.na(white));

# impute the two missing values in BMI;
BMI_mean=mean(BMI, na.rm=T);
BMI_sd=sd(BMI, na.rm=T);

set.seed(197789);
BMI[61]=rnorm(1, mean=BMI_mean, sd=BMI_sd);

set.seed(197552);
BMI[81]=rnorm(1, mean=BMI_mean, sd=BMI_sd);

# log transformation of HgA1c;

ln_hga1c=HgA1c.o=log(HgA1c)

y_miss=ln_hga1c;

# prepare the dataset;

obs_no=rowobs-sum(ma1c);
mis_no=sum(ma1c);


a1c_obs=cbind(Diabetes, Age, BMI, white, y_miss)[ma1c==0,];
x_mis=cbind(Diabetes, Age, BMI, white)[ma1c==1,];
y_completed_sens=y_miss;

# beta_draw=parameter_draw(a1c_obs);

# MNAR imputation sensitivity analysis;
# the number of imputations;
mi_no=50;

# sensitivity parameter vector;

c_beta_vec=seq(-1,1,0.1);

# Results from Table 13.6

for (l in 1:length(c_beta_vec))
{
c_beta=c_beta_vec[l];

MI_mean_vec=MI_mean_var_vec=rep(NA, mi_no);

for (i in 1:mi_no)
{
 set.seed(i);
 # paramter draws from observed data;

  proper_draw=parameter_draw(a1c_obs);
  beta_draw_ori=proper_draw[1:5];
  beta_diabetes_draw=beta_draw_ori[2];
  beta_diabetes_new=beta_diabetes_draw*(1+c_beta);
  beta_draw=beta_draw_ori;
  beta_draw[2]=beta_diabetes_new;
  sigma2_draw=proper_draw[6];

  y_imputed_sens=cbind(rep(1,mis_no), x_mis)%*%beta_draw+sqrt(sigma2_draw)*rnorm(mis_no);
  y_completed_sens[ma1c==1]=y_imputed_sens;




  # estimate the mean;
  MI_mean_vec[i]=mean(exp(y_completed_sens));
  MI_mean_var_vec[i]=var(exp(y_completed_sens))/rowobs;
# MI_mean_vec[i]=mean(y_completed_sens);
# MI_mean_var_vec[i]=var(y_completed_sens)/rowobs;



}
 
# MI_mean_vec;

# combine the multiple imputation estimates;
mean_summary=pool.scalar(MI_mean_vec, MI_mean_var_vec, n=rowobs);
mean_mi_mean=mean_summary$qbar;
mean_mi_var=mean_summary$t;
mean_mi_df=mean_summary$df;
mean_mi_f=mean_summary$f;

MI_mean_low95=mean_mi_mean-qt(.975, mean_mi_df)*sqrt(mean_mi_var);
MI_mean_up95=mean_mi_mean+qt(.975, mean_mi_df)*sqrt(mean_mi_var);

cat("the c_beta is", c_beta, "\n");

# mean_mi_mean;
# MI_mean_low95;
# MI_mean_up95;

cat("the mean_mi_mean is", mean_mi_mean, "\n");
cat("the lower bound is", MI_mean_low95, "\n");
cat("the upper bound is", MI_mean_up95, "\n");


}

