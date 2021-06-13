# Example 13.3

# selection model simulation
# library(car);
library(mvtnorm);
library(MASS);

# set the random seed;
set.seed(197789);

rowobs=10000;
mu_y=0;
mu_x=0;
sigma2_y=1;

# correlation coefficient in the Heckman model;
rho=0;
# rho=0.2;
# rho=0.4
# rho=0.6
# rho=0.8;

Cov_matrix=matrix(c(sigma2_y, sqrt(sigma2_y)*rho, sqrt(sigma2_y)*rho,1), nrow=2, ncol=2);
mu=c(mu_y, mu_x);

# generate complete data from a bivariate normal distribution;
data=mvrnorm(n = rowobs, mu=mu, Sigma=Cov_matrix);

# extract data;
y=data[,1];

x=data[,2];

# observed values;
y_obs=y[x>0];
obs_no=length(y_obs);

# Fig. 13.1 1st row;
hist(y_obs, xlab="rho=0", main="");
abline(v=mean(y_obs), col="red");
qqnorm(y_obs, xlab="rho=0", main="");
qqline(y_obs);

# Fig. 13.1 2nd row;
# hist(y_obs, xlab="rho=0.2", main="");
# abline(v=mean(y_obs), col="red");
# qqnorm(y_obs, xlab="rho=0.2", main="");
# qqline(y_obs);

# Fig. 13.1 3rd row;
# hist(y_obs, xlab="rho=0.4", main="");
# abline(v=mean(y_obs), col="red");
# qqnorm(y_obs, xlab="rho=0.4", main="");
# qqline(y_obs);

# Fig. 13.1 4th row;
# hist(y_obs, xlab="rho=0.6", main="");
# abline(v=mean(y_obs), col="red");
# qqnorm(y_obs, xlab="rho=0.6", main="");
# qqline(y_obs);

# Fig. 13.1 5th row;
# hist(y_obs, xlab="rho=0.8", main="");
# abline(v=mean(y_obs), col="red");
# qqnorm(y_obs, xlab="rho=0.8", main="");
# qqline(y_obs);


# Fig. 13.4

# choose rho=0.8
hist(y, main="", xlab="Y_obs: tilted bars; rho=0.8")
hist(y_obs, add=T, density=25, angle=60 )

