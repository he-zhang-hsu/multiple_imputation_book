# example 3.4

library(mice);


# generate 100 samples from N(0,1);

set.seed(197789);

n=1000;

y=rnorm(n,1);

# remove the last 500 cases;

n_obs=500;

y_obs=y[1:n_obs];

n_mis=n-n_obs;

# complete case analysis;
CC_mean=mean(y_obs);
CC_var=var(y_obs)/n_obs;

CC_mean;
CC_var;

# assign missing data;
miss_indi=c(rep(0, n_obs), rep(1, n_obs));

y_miss=y;
y_miss[miss_indi==1]=NA;


# generate imputed data;

# obtain observed-data sufficient statistics;
mu_hat=CC_mean;
s_square=1/(n_obs-1)*crossprod(y_obs-mu_hat);

# multiple imputation;

# M=2;

# M=5;

# M=20;

# M=50;

# M=200;

# M=1000;

M=10000;

MI_mean_vec=rep(NA, M);
MI_mean_var_vec=rep(NA, M);

for (i in 1:M)
{
set.seed(i);

# draw sigma_square;
sigma_square=1/(rgamma(1, shape=n_obs/2-1/2, scale=2/((n_obs-1)*s_square)));
mu=rnorm(1, mu_hat, sd=sqrt(sigma_square/n_obs));

y_imputed=rnorm(n_mis, mu, sd=sqrt(sigma_square));

y_completed=c(y_obs, y_imputed);

MI_mean_vec[i]=mean(y_completed);
MI_mean_var_vec[i]=var(y_completed)/n;



}

# pool the summary results;

mean_summary=pool.scalar(MI_mean_vec, MI_mean_var_vec, n=n);

mean_summary;

hist(mean_summary$u);

hist(mean_summary$qhat);

