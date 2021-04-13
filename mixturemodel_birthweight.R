################## data generate ##############
# Program: boxcox_impute_expolore_simulation; ############
# generate data, and try out different missing data methods;
# run some simulations to assess the biasedness of various methods;

# library(car);
# library(mvtnorm);
library(mice);
library(norm);
library(R2WinBUGS);
# library(HI);
options(digits=4);
rm(list=ls());

model0.file <- system.file(package="R2WinBUGS", "model", "mixturemodel_complete.txt")
# Let's take a look:
file.show(model0.file)

mode20.file <- system.file(package="R2WinBUGS", "model", "mixturemodel_impute.txt")
# Let's take a look:
file.show(mode20.file)




data=read.table(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter6_multivariate_jointmodeling\\program\\week30.txt", col.names=c("week", "dbirwt", "dplural", "sex"), na.strings=".");


rowobs=nrow(data);

hist(data$dbirwt, xlab="birthweight", main="");

y=sort(data$dbirwt);

mean(y);
sd(y);
summary(y);


T=rep(NA, rowobs);
T[1]=1;
T[rowobs]=2;

# output the data;
# varname_cov=c("y[]", "T[]");
# write.table(cbind(y, T), file=paste("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter6_multivariate_jointmodeling\\program\\birthweightforbugs.dat",sep=""), row.name=FALSE, col.name=varname_cov);

# fit mixture models to the data;

N=rowobs;
alpha=c(1,1);

process_stage2_data_simulate=list("y","T","N","alpha");


process_stage2_init=function(){list(lambda = c(1500, NA),theta=1500,sigma1=400,sigma2=500,P=c(0.5, 0.5))};

# process_stage2_parameters=c("mu.beta", "sigma", "Y_impute");
# process_stage2_parameters=c("mu.beta", "sigma", "Y");
process_stage2_parameters=c("mean1", "mean2", "sigma1", "sigma2", "P1");

gibbs_no=20000;


process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file=model0.file,
    n.chains=1, n.thin=10, n.iter=gibbs_no, n.burnin=gibbs_no/2, bugs.seed=197789, bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);



attach.bugs(process_stage2.sim);
rm(process_stage2.sim);

mean(mean1);
sd(mean1);
mean(mean2);
sd(mean2);

plot(mean1, type="l");

mean(sigma1);
sd(sigma1);
mean(sigma2);
sd(sigma2);

plot(sigma1, type="l");


mean(sigma1^2);
sd(sigma1^2);
mean(sigma2^2);
sd(sigma2^2);

mean(P1);
sd(P1);


# assing missing values to Y;
set.seed(19770809);
miss_indi=runif(n=rowobs)<0.30;
y_miss=y;
y_miss[miss_indi==1]=NA;

hist(y_miss, xlab="observed birthweight", main="");

# write.table(cbind(y_miss, T), file=paste("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter6_multivariate_jointmodeling\\program\\birthweightforbugs_miss.dat",sep=""), row.name=FALSE, col.name=varname_cov);

# fit mixture models for incomplete data;

N=rowobs;
alpha=c(1,1);

process_stage2_data_simulate=list("y_miss","T","N","alpha");


process_stage2_init=function(){list(lambda = c(1500, NA),theta=1500,sigma1=400,sigma2=500,P=c(0.5, 0.5))};

# process_stage2_parameters=c("mu.beta", "sigma", "Y_impute");
# process_stage2_parameters=c("mu.beta", "sigma", "Y");
process_stage2_parameters=c("mean1", "mean2", "sigma1", "sigma2", "P1");

gibbs_no=20000;


process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file=mode20.file,
    n.chains=1, n.thin=10, n.iter=gibbs_no, n.burnin=gibbs_no/2, bugs.seed=197789, bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);

attach.bugs(process_stage2.sim);
rm(process_stage2.sim);



mean(mean1);
sd(mean1);
mean(mean2);
sd(mean2);

plot(mean1, type="l");

mean(sigma1);
sd(sigma1);
mean(sigma2);
sd(sigma2);

plot(sigma1, type="l");


mean(sigma1^2);
sd(sigma1^2);
mean(sigma2^2);
sd(sigma2^2);

mean(P1);
sd(P1);


# multiple imputation using mixture models;

N=rowobs;
alpha=c(1,1);

process_stage2_data_simulate=list("y_miss","T","N","alpha");


process_stage2_init=function(){list(lambda = c(1500, NA),theta=1500,sigma1=400,sigma2=500,P=c(0.5, 0.5))};

# process_stage2_parameters=c("mu.beta", "sigma", "Y_impute");
# process_stage2_parameters=c("mu.beta", "sigma", "Y");
process_stage2_parameters=c("y_impute");

gibbs_no=20000;
mi_no=20;

process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file=mode20.file,
    n.chains=1, n.thin=gibbs_no/(2*mi_no), n.iter=gibbs_no, n.burnin=gibbs_no/2, bugs.seed=197789, bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);

attach.bugs(process_stage2.sim);
rm(process_stage2.sim);

# plot the imputed data using normal mixture distribution;

hist(y_impute[1,], xlab="imputed birthweight: mixture model", main="");

# multiple imputation using the normal model;
# obtain observed-data sufficient statistics;
y_obs=y_miss[miss_indi==0];
n_obs=length(y_obs);
n_mis=rowobs-n_obs;

mu_hat=mean(y_obs);

s_square=1/(n_obs-1)*crossprod(y_obs-mu_hat);

for (i in 1:mi_no)
{
set.seed(i);

# draw sigma_square;
sigma_square=1/(rgamma(1, shape=n_obs/2-1/2, scale=2/((n_obs-1)*s_square)));
mu=rnorm(1, mu_hat, sd=sqrt(sigma_square/n_obs));

y_imputed=rnorm(n_mis, mu, sd=sqrt(sigma_square));

y_completed=c(y_obs, y_imputed);

}

# plot the imputed data using the normal distribution;

hist(y_imputed, xlab="imputed birthweight: normal model", main="");





