# Example 8.6

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
# getNimbleOptions();
nimbleOptions(clearNimbleFunctionsAfterCompiling=TRUE);

# nmlnclb generate survival data;

nm1nc1b <- function( n, beta, rho, a,g, c,tao )
{
  ok <- 0
  while( ok == 0 )
    {
### generating the data
### the covariates
zz=mvrnorm(n,mu=c(0,0),Sigma=matrix(c(1,rho,rho,1),2,2))
zm=zz[,1]
zc=zz[,2]
z <- cbind(zm, zc)
### 
risk <- exp( z %*% beta )
u1 <- runif(n)
t <- c * ( -log(1-u1) / risk )^(1/tao)
u2 <- runif(n)
ctime <- c * ( -log(1-u2) )

case <- rep(0, n)
case[ t <= ctime ] <- 1

x <- t
x[case==0] <- ctime[case==0]

# temp <- quantile( x, 0.95 )
# case[ x >= temp ] <- 0

###  assign the conditional selection probability
sele.prob <- 1/(1 + exp(a + g*case ))
u1 <- runif( n )
sele.ind <- rep(0, n)
sele.ind[ u1 < sele.prob ] <- 1 
   ok <- ifelse( min(table(sele.ind, case)) >= 5, 1, 0 ) 
  }
data <- data.frame( x, zm, zc, case, sele.ind ) 

return( data )
}


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
estimand=as.data.frame(cbind(bias,rbias,mse,mean_length,mean_cov,sd,se))
}


#######################################################
# simulation;

cycle_no=100;
mi_no=50;
niter=5000;
burnin=niter/2;


# matrices holding the parameter estimates;
# BD means using true y;
BD_m_slope_vec=BD_m_slope_var_vec=BD_c_slope_vec=BD_c_slope_var_vec=rep(NA, cycle_no);
# CC means using y from validation data;
CC_m_slope_vec=CC_m_slope_var_vec=CC_c_slope_vec=CC_c_slope_var_vec=rep(NA, cycle_no);
# MI means using multiply imputed y;
MI_m_slope_mat=MI_m_slope_var_mat=MI_c_slope_mat=MI_c_slope_var_mat=matrix(NA, cycle_no, mi_no);

# population parameters;

n=100;
shape=1;
beta=c( 2, 2 );
rho=0.5;


# R nimble code;
# initialize the object;
# generate survival data;

# nimble;


for (cycle in 1:cycle_no)

{
set.seed(2*cycle);

# generate survival data;
# cycle=1;

data=nm1nc1b(n, beta, rho, 2,-4, 1,shape);

observetime=data$x;
zm=data$zm
zc=data$zc
case=data$case
sele.ind=data$sele.ind
# assign missing data;
zm[sele.ind==0]=NA;

failtime=c(sort(observetime[case==1]), max(observetime));


# compelete data (using true y) analysis;

# before deletion;
zz=cbind(data$zm, data$zc);

coxfit1=coxph(Surv(observetime, case) ~ zz)
beta.zm=coxfit1$coef[1];
beta.zc=coxfit1$coef[2];
se.zm=sqrt(coxfit1$var[1,1]);
se.zc=sqrt(coxfit1$var[2,2]);

BD_m_slope_vec[cycle]=beta.zm;
BD_m_slope_var_vec[cycle]=se.zm^2;

BD_c_slope_vec[cycle]=beta.zc;
BD_c_slope_var_vec[cycle]=se.zc^2;

# different missing data methods;

# complete-case analysis;
zz_mis=cbind(zm, zc);

coxfit2=coxph(Surv(observetime, case) ~ zz_mis)
beta.zm.mis=coxfit2$coef[1];
beta.zc.mis=coxfit2$coef[2];
se.zm.mis=sqrt(coxfit2$var[1,1])
se.zc.mis=sqrt(coxfit2$var[2,2]);

CC_m_slope_vec[cycle]=beta.zm.mis;
CC_m_slope_var_vec[cycle]=se.zm.mis^2;

CC_c_slope_vec[cycle]=beta.zc.mis;
CC_c_slope_var_vec[cycle]=se.zc.mis^2;


# multiple imputation analysis;

# setup the data for BUGS;
# setup the data
N=n;
T=length(failtime)-1;
Y=matrix(NA, nrow=N, ncol=T);
dN=matrix(NA, nrow=N, ncol=T);
eps=1E-10;

# Set up data
for(i in 1:N) 
{
for(j in 1:T) 
{
# risk set = 1 if obs.t >= t
Y[i,j] <- as.numeric((observetime[i] - failtime[j] + eps)>=0)

	# counting process jump = 1 if obs.t in [ t[j], t[j+1] )
	#                      i.e. if t[j] <= obs.t < t[j+1]
dN[i, j] <- Y[i,j] * as.numeric((failtime[j + 1] - observetime[i]-eps)>=0)*case[i]
}
}
# Y
# dN
# nimble;

survival_impute=nimbleCode(
{
#model
# {
	# Model 
		for(j in 1:T) {
#            beta0[j] ~ dnorm(0, 0.001);
			for(i in 1:N) {
				dN[i, j]   ~ dpois(Idt[i, j])              # Likelihood
		        Idt[i, j] <- Y[i, j] * exp(betaM * zm[i]+betaC*zc[i]) * dL0[j] 	# Intensity 
			}     
			dL0[j] ~ dgamma(mu[j], c)
			mu[j] <- dL0.star[j] * c    # prior mean hazard

		}
		for (i in 1:N)
		{
		zm[i] ~ dnorm(mu_zm[i], tau);
	     mu_zm[i]<-alpha0+alpha1*zc[i];
	}
		c <- 0.001
		r <- 0.1 
		for (j in 1 : T) 
{  dL0.star[j] <- r * (failtime[j + 1] - failtime[j])  } 
		betaM ~ dnorm(0.0,0.001)
		betaC ~ dnorm(0.0, 0.001)              
 	   alpha0 ~ dnorm(0.0, 0.001);
		alpha1 ~ dnorm(0.0, 0.001);
		 tau ~ dgamma(.001, .001);
		sigma <- sqrt(1/tau);
#	}
}
)


survival.constants=list(N=n, T=length(failtime)-1, eps=1.0E-10);

survival.data=list(observetime=observetime, Y=Y, dN=dN, case=case, failtime=failtime, zm=zm, zc=zc);

survival.inits=list(betaM = 2, betaC=2, alpha0=0, alpha1=0.5, tau=1, zm=rnorm(n), dL0=rep(1,length(failtime)-1));



# compiling;
survival_impute_model=nimbleModel(survival_impute, constants=survival.constants, 
data=survival.data, inits=survival.inits);

survival.configure=configureMCMC(survival_impute_model);

# survival.configure;

survival.configure$resetMonitors();

survival.configure$addMonitors('zm');

survival.configure$setThin(burnin/mi_no);

survival.mcmc=buildMCMC(survival.configure);

survival_impute_compiled=compileNimble(survival_impute_model, survival.mcmc);

# runMCMC(survival_impute_model, niter=100);

survival_impute_compiled$survival.mcmc$run(niter);

results=as.matrix(survival_impute_compiled$survival.mcmc$mvSamples);
no_mi_2=nrow(results);
# throw away the burn-in samples;
results_final=results[(no_mi_2/2+1):no_mi_2,];
# head(results);

# release memory?


# nimble:::clearCompiled(survival_impute_model);

rm(survival_impute_model)
rm(survival_impute_compiled);
rm(Y);
rm(dN);

# no compiling


# survival_impute.mcmc=nimbleMCMC(code=survival_impute, data=survival.data, 
# constants=survival.constants, nchains=1, inits=survival.inits, niter=niter, nburnin=niter/2, 
# summary=TRUE, thin=(niter-burnin)/mi_no, WAIC=FALSE, monitors=c('zm'));

# completed data;
# results=as.data.frame(survival_impute.mcmc$samples);

# multiple imputation analysis;

for (i in 1:mi_no)
{
zm_completed=results_final[i,];
zz_completed=cbind(zm_completed, zc);

# completed-data analysis;
coxfit3=coxph(Surv(observetime, case) ~ zz_completed)

beta.zm.impute=coxfit3$coef[1];
beta.zc.impute=coxfit3$coef[2];
se.zm.impute=sqrt(coxfit3$var[1,1]);
se.zc.impute=sqrt(coxfit3$var[2,2]);

MI_m_slope_mat[cycle,i]=beta.zm.impute;
MI_m_slope_var_mat[cycle,i]=se.zm.impute^2;


MI_c_slope_mat[cycle,i]=beta.zc.impute;
MI_c_slope_var_mat[cycle,i]=se.zc.impute^2;


}



cat("the cycle is", cycle, "\n");


}

# assess the performance;

##########################################################################
# Simulation results;
# Table 8.8

# complete-data analysis;

true_m_slope=mean(BD_m_slope_vec);
true_m_slope;

BD=performance(BD_m_slope_vec, BD_m_slope_var_vec, true_m_slope);
BD;
CC=performance(CC_m_slope_vec, CC_m_slope_var_vec, true_m_slope);
CC;
MI=mi_performance(MI_m_slope_mat, MI_m_slope_var_mat, true_m_slope, n);
MI;



true_c_slope=mean(BD_c_slope_vec);
true_c_slope;

BD=performance(BD_c_slope_vec, BD_c_slope_var_vec, true_c_slope);
BD;
CC=performance(CC_c_slope_vec, CC_c_slope_var_vec, true_c_slope);
CC;
MI=mi_performance(MI_c_slope_mat, MI_c_slope_var_mat, true_c_slope, n);
MI;


