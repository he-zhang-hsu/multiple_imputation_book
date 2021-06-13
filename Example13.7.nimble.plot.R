# Example 13.7 plot

rm(list=ls());
library(coda)
library(lattice)
library(gridExtra)
library(MASS)
library(Matrix)
library(mice)
library(mitools)
library(survival)
library(nimble);
library(FKF);
library(MASS);
library(pryr);


# Set up data

# nimble;

# set up the data;
set.seed(197789);
rowobs=1000;
M=rowobs;
mu_x=1;
var_x=1;
x=mu_x+sqrt(var_x)*rnorm(rowobs);



# regression parameters; 
beta0=-2;
beta1=1;

# error variance;
var_error=1;

# normal error distribution
error=sqrt(var_error)*rnorm(rowobs);

# generate the outcome variable;
y=beta0+beta1*x+error;

# set up the missing data as MNAR on Y only;
# alpha0=-1.4;
# alpha1=-1;
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*y)/(1+exp(alpha0+alpha1*y)));

# set up the missing data as MNAR on both Y and X;
alpha0=-1.4;
alpha1=-0.75;
alpha2=0.25;
miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*y+alpha2*x)/(1+exp(alpha0+alpha1*y+alpha2*x)));

m=miss_indi;

y_miss=y;
y_miss[miss_indi==1]=NA;

# model for mnar only related to Y;

selection_impute=nimbleCode(
# MODEL   
 {			
		# complete-data likelihood;
			for( i in 1 : M ) {
			# outcome model;
				y_miss[i] ~ dnorm(mu[i], tau);
				mu[i] <- beta0+beta1*x[i];
			# response model
			     m[i] ~ dbern(p[i]);
				logit(p[i]) <-alpha0+alpha1*y_miss[i];	
                    y_impute[i] <- y_miss[i];
         		}
			# prior for complete-model parameters;
			beta0 ~ dnorm(0.0, 1.0E-3);
			beta1 ~ dnorm(0.0, 1.0E-3);
			alpha0 ~ dnorm(0.0, 1.0E-3);
			alpha1 ~ dnorm(0.0, 1.0E-3);
			tau ~ dgamma(1.0E-3, 1.0E-3);
			
}
)

############################
# model for mnar related to both Y and X;
selection_impute=nimbleCode(
# MODEL   
# this model estimates the parameters of a logistic model with a mismeasured predictor
# for missing values;
# MODEL   
{			
		# complete-data likelihood;
			for( i in 1 : M ) {
			# outcome model;
				y_miss[i] ~ dnorm(mu[i], tau);
				mu[i] <- beta0+beta1*x[i];	
			# response model
			     m[i] ~ dbern(p[i]);
				logit(p[i]) <-alpha0+alpha1*x[i]+alpha2*y_miss[i];	
                    y_impute[i] <- y_miss[i];
         		}
			# prior for complete-model parameters;
			beta0 ~ dnorm(0.0, 1.0E-3);
			beta1 ~ dnorm(0.0, 1.0E-3);
		#	alpha0 ~ dlogis(0,1);
	          alpha0 ~ dnorm(0.0, 1.0E-3);	
			alpha1 ~ dnorm(0.0, 1.0E-3);
			alpha2 ~ dnorm(0.0, 1.0E-3);
			tau ~ dgamma(1.0E-3, 1.0E-3);			
}

)



selection.constants=list(M=M);

selection.data=list(y_miss=y_miss, m=m, x=x);

# for mmnar model relates only to Y;
# selection.inits=list(beta0=0, beta1=0, tau=1, alpha0=0, alpha1=0);

# for mnar model relates to both Y and X;
selection.inits=list(beta0=0, beta1=0, tau=1, alpha0=0, alpha1=0, alpha2=0);


# the approach compiling the data;

selection_impute_model=nimbleModel(selection_impute, constants=selection.constants, 
data=selection.data, inits=selection.inits);

# selection_impute_model$getNodeNames();
selection.configure=configureMCMC(selection_impute_model);

# for mnar model relates only to Y;
# selection.configure$addMonitors('alpha0', 'alpha1', 'tau', 'beta0', 'beta1');

# for mnar model relates to both Y and X;
selection.configure$addMonitors('alpha0', 'alpha1', 'alpha2', 'tau', 'beta0', 'beta1');

selection.configure;


selection.mcmc=buildMCMC(selection.configure);

# compile the Nimble project;

impute=compileNimble(selection_impute_model, selection.mcmc);

niter=20000;
set.seed(0);

# run MCMC;
impute$selection.mcmc$run(niter);

# MCMC results;
results=as.matrix(impute$selection.mcmc$mvSamples);
head(results);

# results[,"alpha0"];
# results[,"alpha1"];
# results[,"beta0"];
# results[,"beta1"];

results[,"alpha0"];
results[,"alpha1"];
results[,"alpha2"];
results[,"beta0"];
results[,"beta1"];


# posterior summary;
# summary(results[(niter/2):niter,c("beta0", "beta1", "alpha0", "alpha1")]);

summary(results[(niter/2):niter,c("beta0", "beta1", "alpha0", "alpha1", "alpha2")]);

# Fig. 13.2
# plot(seq(niter/2,niter,1), type="l", results[(niter/2):niter, "beta0"], xlab="iteration", ylab="beta0");

plot(seq(niter/2,niter,1), type="l", results[(niter/2):niter,"beta1"], xlab="iteration", ylab="beta1");

plot(seq(niter/2,niter,1), type="l", results[(niter/2):niter,"alpha0"], xlab="iteration", ylab="gamma0");

plot(seq(niter/2,niter,1), type="l", results[(niter/2):niter,"alpha1"], xlab="iteration", ylab="gamma1");

plot(seq(niter/2,niter,1), type="l", results[(niter/2):niter,"alpha2"], xlab="iteration", ylab="gamma2");

