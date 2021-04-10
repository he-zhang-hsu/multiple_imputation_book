# example 3.1

# generate 100 samples from N(0,1);

n=1000;

y=rnorm(n,1);

# remove the last 25 cases;

n_obs=500;

y_obs=y[1:n_obs];

n_mis=n-n_obs;

sample=100000;

# generate the posterior distribution of the mean of y_obs;

y_obs_bar=mean(y_obs);

set.seed(197789);

mu_obs=rnorm(sample, mean=y_obs_bar, sd=sqrt(1/n_obs));

hist(mu_obs);
mean(mu_obs);
var(mu_obs);


# generate imputed data;
M=50;

mu_mi_matrix=matrix(NA, nrow=sample/M, ncol=M);

for (i in 1:M)
{
set.seed(i);
mu_draw=rnorm(1, mean=y_obs_bar, sd=sqrt(1/n_obs));
y_impute=rnorm(n_mis, mean=mu_draw, sd=sqrt(1));
y_completed=c(y_obs, y_impute);

mu_mi_matrix[,i]=rnorm(sample/M, mean=mean(y_completed), sd=sqrt(1/n));
if (i==1) mu_mi_matrix_1=rnorm(sample, mean=mean(y_completed), sd=sqrt(1/n));

}

hist(mu_obs);

hist(mu_mi_matrix[,1]);

par(mfrow=c(1,3))
hist(mu_mi_matrix, xlab="", main="", xlim=c(.75, 1.2));
hist(mu_obs, xlab="", main="", xlim=c(.75, 1.2));
hist(mu_mi_matrix_1, xlab="", main="", xlim=c(.75, 1.2));

par(mfrow=c(1,1));

mean(mu_obs);
var(mu_obs);
mean(mu_mi_matrix);
var(as.vector(mu_mi_matrix));

mean(mu_mi_matrix_1);
var(mu_mi_matrix_1);


summary(mu_obs);
summary(as.vector(mu_mi_matrix));
summary(mu_mi_matrix_1);


