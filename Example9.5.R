# Example 9.5

rm(list=ls());
library(R2WinBUGS);
library(MASS);
library(mice);
library(norm);
library(mix);
# library(HI);
options(digits=4);

# WinBUGS syntax

# Imputation program assuming correlated random intercepts and slopes;

model0.file <- system.file(package="R2WinBUGS", "model", "logisticmixed_imputeforR.txt")
# Let's take a look:
file.show(model0.file)

# Imputation program assuming only random intercepts;

model1.file <- system.file(package="R2WinBUGS", "model", "logisticmixed_ranint_only_imputeforR.txt")
# Let's take a look:
file.show(model1.file)

# Imputation program assuming independent random intercepts and slopes;

model2.file <- system.file(package="R2WinBUGS", "model", "logisticmixed_ranint_slope_imputeforR.txt")
# Let's take a look:
file.show(model2.file)

# Imputation program assuming common intercepts and slopes;

model3.file <- system.file(package="R2WinBUGS", "model", "logisticmixed_fixed_imputeforR.txt")
# Let's take a look:
file.show(model3.file)



# read the data;

oridata=scan(file=paste("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter9_longitudinal\\program\\psid_uni_uni_hr_2172_log.dat"), na.strings=".");

orimatrix=matrix(oridata, nrow=2172, ncol=34, byrow=TRUE);

# oriy concatenate the response vector;
oriy=orimatrix[,13:31];

oriyvec=as.vector(t(oriy));

summary(oriyvec);

# create binary variable;

oriy_binary=1*(oriy>=2);
mean(oriy_binary, na.rm=T);
sum(is.na(oriy_binary))/length(oriy_binary);

# create the data;

N=2172;
T=19;
x = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19);
xbar=10;
Y=oriy_binary;

process_stage2_data_simulate=list("Y","x","xbar", "T","N");

# correlated random intercept and slope;

# create the initial values;

process_stage2_init=function(){list(beta.c  = 0.197, alpha.c = 0.046, slope=0, sigma.alpha=2.187, sigma.beta.alpha=0.298)};

process_stage2_parameters=c("alpha.c", "alpha", "beta.c", "beta", "sigma.alpha", "sigma.beta.alpha", "slope");

# process_stage2_data_simulate=list("M","y2_miss", "x1_miss", "x2_miss");
# process_stage2_init= function(){list(beta0=-3, beta1=1.5, beta2=1.5, mu=0, alpha0=0, alpha1=0, tau=1, psi=1)};
# process_stage2_parameters=c("mu", "beta0", "beta1", "beta2", "mu", "alpha0", "alpha1", "tau", "psi", "y2", "x1", "x2");


gibbs_no=10000;


process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file=model0.file,
    n.chains=1, n.thin=100, n.iter=gibbs_no, n.burnin=gibbs_no/2, bugs.seed=197789, bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);

names(process_stage2.sim);

# Table 9.4;
# dic value;

process_stage2.sim$pD;
process_stage2.sim$DIC;


# dim(process_stage2.sim$sims.array);

# process_stage2.sim$sims.array[,,"alpha[1]"];

# process_stage2.sim$sims.array[,,1:10];


attach.bugs(process_stage2.sim);

# posterior means;
alpha.mean=colMeans(alpha);
alpha.c.mean=mean(alpha.c);
beta.mean=colMeans(beta);
beta.c.mean=mean(beta.c);

alpha.mean;
alpha.c.mean;
beta.mean;
beta.c.mean;
alpha.c.mean-beta.c.mean*xbar;

# plot the posterior means of the intercept;
alpha_original=alpha.mean+alpha.c.mean-beta.c.mean*xbar-beta.mean*xbar;

# Fig. 9.5
# Top left
hist(alpha_original, xlab="Posterior Means of Random Intercepts at Baseline", main="");

summary(alpha_original);

beta_original=beta.mean+beta.c.mean;

# Top right
hist(beta_original, xlab="Posterior Means of Random Slopes", main="");

summary(beta_original);

# Bottom left
plot(alpha_original, beta_original, xlab="Posterior Means of Intercepts", ylab="Posterior Means of Slopes", main="");


# plot trajectories for some families;
coef=cbind(alpha_original, beta_original);
time_mat=rbind(1, x);

data_time=coef%*%time_mat;

# Bottom right;
matplot(x+1978, t(data_time[1:25,]), col=rep("black", 25), lty=1, lwd=rep(1,25), type="l", xlab="time", ylab="Individual Trajectory at Logit Scale");

# work on the variance components;
var_alpha=sigma.alpha^2;
var_beta=slope^2*var_alpha+sigma.beta.alpha^2;
rho=sqrt(slope^2*var_alpha/var_beta);
cov_alpha_beta=rho*sqrt(var_alpha*var_beta);
mean(var_alpha);
mean(var_beta);
mean(rho);
mean(cov_alpha_beta);


var_alpha_original=var_alpha-2*xbar*rho*sqrt(var_alpha*var_beta)+xbar^2*var_beta;
var_beta_original=var_beta;
cov_alpha_beta_original=rho*sqrt(var_alpha*var_beta)-xbar*var_beta;
rho_alpha_beta_original=cov_alpha_beta_original/sqrt(var_alpha_original*var_beta_original);

mean(var_alpha_original);
mean(var_beta_original);
mean(rho_alpha_beta_original);


# process_stage2.sim;


# random intercept model;
# create the initial values;

process_stage2_int_init=function(){list(beta.c=0.197, alpha.c=0.046, sigma.alpha=2.187)};

process_stage2_int_parameters=c("alpha.c", "alpha", "beta.c", "sigma.alpha");


gibbs_no=10000;


process_stage2_int.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_int_init, parameters=process_stage2_int_parameters, model.file=model1.file,
    n.chains=1, n.thin=100, n.iter=gibbs_no, n.burnin=gibbs_no/2, bugs.seed=197789, bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);

names(process_stage2_int.sim);

# dic value;
# Table 9.4;
process_stage2_int.sim$pD;
process_stage2_int.sim$DIC;

# random intercept and slope model, assuming independence
# create the initial values;

process_stage2_int_slope_init=function(){list(beta.c=0.197, alpha.c=0.046, sigma.alpha=2.187, sigma.beta.alpha=0.298)};

process_stage2_int_slope_parameters=c("alpha.c", "alpha", "beta.c", "sigma.alpha", "sigma.beta.alpha");

gibbs_no=10000;

process_stage2_int_slope.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_int_slope_init, parameters=process_stage2_int_slope_parameters, model.file=model2.file,
    n.chains=1, n.thin=100, n.iter=gibbs_no, n.burnin=gibbs_no/2, bugs.seed=197789, bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);

names(process_stage2_int_slope.sim);

# dic value;
# Table 9.4;
process_stage2_int_slope.sim$pD;
process_stage2_int_slope.sim$DIC;


# common intercept and slope model;

# create the initial values;

process_stage2_fixed_init=function(){list(beta.c=0.197, alpha.c=0.046)};

process_stage2_fixed_parameters=c("alpha.c", "beta.c");

gibbs_no=10000;


process_stage2_fixed.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_fixed_init, parameters=process_stage2_fixed_parameters, model.file=model3.file,
    n.chains=1, n.thin=100, n.iter=gibbs_no, n.burnin=gibbs_no/2, bugs.seed=197789, bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);

names(process_stage2_fixed.sim);

# dic value;
# Table 9.4;
process_stage2_fixed.sim$pD;
process_stage2_fixed.sim$DIC;




