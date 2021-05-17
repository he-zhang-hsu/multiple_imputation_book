# library(car);
# library(mvtnorm);
library(mice);
library(norm);
# library(HI);
options(digits=4);
rm(list=ls());

# complete-data sample size is 1000;

rowobs=1000;

# set up the random seed;
set.seed(197789);

# generate data;
# distribution of x, the covariate, is from N(1,1);
mu_x=1;
var_x=1;
x=mu_x+sqrt(var_x)*rnorm(rowobs);


# error variance;
var_error=1;

# normal error distribution
error=sqrt(var_error)*rnorm(rowobs);


# regression parameters;
beta0=-2;
beta1=1;

# generate the outcome variable;
y=beta0+beta1*x+error;

# before deletion;
# Table 1.10
bd_mean=mean(y);

bd_mean;

bd_se=sd(y)/sqrt(rowobs);

bd_se;

bd_fit=lm(y~x);

# Fig. 1.3
plot(x,y, main="Before Deletion");
abline(bd_fit);


# Missing data on y based on MCAR;
set.seed(197789);

# missingness indicator;
miss_indi=rbinom(n=rowobs, size=1, prob=0.4);

y_miss=y;

# assign missing values to y;

y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
x_obs=x[miss_indi==0];
x_mis=x[miss_indi==1];

# Table 1.10

mcar_mean=mean(y_obs, na.rm=T);

mcar_mean;

mcar_se=sd(y_obs)/sqrt(obs_no);

mcar_se;

mcar_fit=lm(y_miss~x);

mcar_logistic=summary(glm(miss_indi~x, family=binomial(link="logit")));

mcar_logistic;

# Fig. 1.3
plot(x,y_miss, ylab="y_mcar", ylim=c(-6, 4), main="MCAR");
abline(mcar_fit);

y_mcar=y_miss;

# Missing data on y based on MAR;

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

# Table 1.10

mar_mean=mean(y_obs, na.rm=T);

mar_mean;

mar_se=sd(y_obs)/sqrt(obs_no);

mar_se;

mar_fit=lm(y_miss~x);

mar_logistic=summary(glm(miss_indi~x, family=binomial(link="logit")));

mar_logistic;


# Fig. 1.3

plot(x,y_miss, ylab="y_mar", ylim=c(-6, 4), main="MAR");
abline(mar_fit);

y_mar=y_miss;

# Missing data on y based on MNAR;

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

# Table 1.10
mnar_mean=mean(y_obs, na.rm=T);

mnar_mean;

mnar_se=sd(y_obs)/sqrt(obs_no);

mnar_se;

mnar_fit=lm(y_miss~x);

mnar_logistic=summary(glm(miss_indi~x, family=binomial(link="logit")));

mnar_logistic;


# Fig. 1.3

plot(x,y_miss, ylab="y_mnar", ylim=c(-6,4), main="MNAR");
abline(mnar_fit);

y_mnar=y_miss;

# box plots for y's in all scenarios;
# Fig. 1.2

y_data=cbind(y, y_mcar, y_mar, y_mnar);

boxplot(y_data);
