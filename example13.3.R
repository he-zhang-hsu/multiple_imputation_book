# selection model simulation
# library(car);
library(mvtnorm);
library(MASS);

# sample size cannot go higher!
set.seed(197789);

rowobs=10000;
mu_y=0;
mu_x=0;
# rho=0;
# rho=0.2;
# rho=0.4
# rho=0.6
rho=0.8;
sigma2_y=1;
Cov_matrix=matrix(c(sigma2_y, sqrt(sigma2_y)*rho, sqrt(sigma2_y)*rho,1), nrow=2, ncol=2);
mu=c(mu_y, mu_x);


# generate bivariate normal distribution;
data=mvrnorm(n = rowobs, mu=mu, Sigma=Cov_matrix);

# extract data;
y=data[,1];

# hist(y, xlab="rho=0", main="");
# abline(v=mean(y), col="red");
# qqnorm(y, xlab="rho=0", main="");
# qqline(y);



x=data[,2];

y_obs=y[x>0];
obs_no=length(y_obs);

# hist(y_obs, xlab="y_obs, rho=0");
# hist(y_obs, xlab="rho=0", main="");
# abline(v=mean(y_obs), col="red");
# qqnorm(y_obs, xlab="rho=0", main="");
# qqline(y_obs);

# hist(y_obs, xlab="rho=0.2", main="");
# abline(v=mean(y_obs), col="red");
# qqnorm(y_obs, xlab="rho=0.2", main="");
# qqline(y_obs);

# hist(y_obs, xlab="rho=0.4", main="");
# abline(v=mean(y_obs), col="red");
# qqnorm(y_obs, xlab="rho=0.4", main="");
# qqline(y_obs);

# hist(y_obs, xlab="rho=0.6", main="");
# abline(v=mean(y_obs), col="red");
# qqnorm(y_obs, xlab="rho=0.6", main="");
# qqline(y_obs);

hist(y_obs, xlab="rho=0.8", main="");
abline(v=mean(y_obs), col="red");
qqnorm(y_obs, xlab="rho=0.8", main="");
qqline(y_obs);

mean(y_obs);


# breaks=seq(-5,9,0.25);

# par(mfrow=c(2,2));
# For plot in Chapter 7 example 5;
# red for the original y;
hist(y, main="", xlab="Y_obs: tilted bars; rho=0.8")
# blue for the observed y;
hist(y_obs, add=T, density=25, angle=60 )

# hist(x, breaks=seq(-4,11,0.5));


y_miss=y[x<0];
miss_no=length(y_miss);

# miss_no is better less than obs_no;

summary(y_obs)-summary(y_miss);

par(mfrow=c(1,1));

par(mfrow=c(2,1));
hist(y);
hist(y_obs);
u=density(y);
plot(u);
v=density(y_obs);
plot(v);

densityPlot(y);


# impute the missing values using the bottom 100*(1-rho) percentile of the observed data;

y_obs_sort=sort(y_obs);
y_miss_impute=y_obs_sort[1:((1-rho)*length(y_obs_sort))];



summary(y_miss_impute)-summary(y_miss);


par(mfrow=c(2,1));
hist(y_miss_impute);
hist(y_miss);



# rejection sampling for y_miss;
n_sample=100000;
y_miss_impute_reject_normal=rnorm(n_sample, mu_y, sd=sqrt(sigma2_y));
sample_selection=(runif(n_sample)<pnorm((-mu_x-rho*(y_miss_impute_reject_normal-mu_y)/sqrt(sigma2_y))/(sqrt(1-rho^2))));
y_miss_impute_reject_final=(y_miss_impute_reject_normal[sample_selection])[1:miss_no];
 
summary(y_miss_impute_reject_final)-summary(y_miss);

par(mfrow=c(2,1));
hist(y_miss_impute_reject_final);
hist(y_miss);

# importance sampling for y_miss;
y_miss_impute_importance_normal=y_miss_impute_reject_normal;
prob_importance=pnorm((-mu_x-rho*(y_miss_impute_reject_normal-mu_y)/sqrt(sigma2_y))/(sqrt(1-rho^2)))/pnorm(-mu_x);
y_miss_impute_importance_final=sample(y_miss_impute_importance_normal, size=miss_no, replace = FALSE, prob = prob_importance/sum(prob_importance));

summary(y_miss_impute_importance_final)-summary(y_miss);

par(mfrow=c(2,1));
hist(y_miss_impute_importance_final);
hist(y_miss);

# importance sampling for y_miss using y_obs;
obs_prob_importance=pnorm((-mu_x-rho*(y_obs-mu_y)/sqrt(sigma2_y))/(sqrt(1-rho^2)))/(1-pnorm((-mu_x-rho*(y_obs-mu_y)/sqrt(sigma2_y))/(sqrt(1-rho^2))));
y_miss_impute_obs_importance=sample(y_obs, size=miss_no, replace = TRUE, prob = obs_prob_importance/sum(obs_prob_importance));

summary(y_miss_impute_obs_importance)-summary(y_miss);

par(mfrow=c(2,1));
hist(y_miss_impute_obs_importance);
hist(y_miss);

# sensitivity analysis;

# set rho;
rho_new=-0.95;
mu_x_new=qnorm(obs_no/(obs_no+miss_no));
# lambda_new is the inverse mills ratio;
lambda_new=dnorm(mu_x_new)/pnorm(mu_x_new);
sigma2_y_new=var(y_obs)/(1-rho_new^2*lambda_new*(lambda_new+mu_x_new));
mu_y_new=mean(y_obs)-rho_new*lambda_new*sqrt(sigma2_y_new);

# importance sampling for y_miss using y_obs;
# use the selection model parameters;
obs_prob_importance=pnorm((-mu_x_new-rho_new*(y_obs-mu_y_new)/sqrt(sigma2_y_new))/(sqrt(1-rho_new^2)))/(1-pnorm((-mu_x_new-rho_new*(y_obs-mu_y_new)/sqrt(sigma2_y_new))/(sqrt(1-rho_new^2))));
y_miss_impute_obs_importance=sample(y_obs, size=miss_no, replace = TRUE, prob = obs_prob_importance/sum(obs_prob_importance));

summary(y_miss_impute_obs_importance)-summary(y_miss);

par(mfrow=c(2,1));
hist(y_miss_impute_obs_importance);
hist(y_miss);
