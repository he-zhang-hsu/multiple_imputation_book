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


improper_imputation=function(data)
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
 # beta=mvrnorm(n=1, mu=beta_hat, Sigma=sigma_square*inv_v_beta);

 parameter=c(beta_hat,sigma_square);
 parameter;
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


# model.file <- system.file(package="R2WinBUGS", "model", "mnar_logit_yonly.txt")
# Let's take a look:
# file.show(model.file)

# model1.file <- system.file(package="R2WinBUGS", "model", "mnar_logit_yandx.txt")
# Let's take a look:
# file.show(model1.file)

model2.file <- system.file(package="R2WinBUGS", "model", "mnar_logit_yandx_new.txt")
# Let's take a look:
file.show(model2.file)




rowobs=1000;

# set up the random seed;
set.seed(197789);

# generate data;
# distribution of x, the covariate;
# could be normal, or could be other distributions;
# fixed at mean 10, and variance 1;
mu_x=1;
var_x=1;

# auxiliary variable z;
mu_z=1;
var_z=3;

# regression parameters;
# both beta0, beta1 and beta2 are fixed at 
beta0=-2;
# beta1=0;
beta1=1;

# error variance;
var_error=1;


cycle_no=1000;
mi_no=50;

# matrices holding the parameter estimates;
# complete-data inferences;
# xyslope means regressing x on y;
BD_mean_vec=BD_mean_var_vec=BD_slope_vec=BD_slope_var_vec=BD_xyslope_vec=BD_xyslope_var_vec=rep(NA, cycle_no);
CC_mean_vec=CC_mean_var_vec=CC_slope_vec=CC_slope_var_vec=CC_xyslope_vec=CC_xyslope_var_vec=rep(NA, cycle_no);
SI_mean_vec=SI_mean_var_vec=SI_slope_vec=SI_slope_var_vec=SI_xyslope_vec=SI_xyslope_var_vec=rep(NA, cycle_no);
INDI_xyslope_vec=INDI_xyslope_var_vec=rep(NA, cycle_no);

# inclusive imputation;
MI_mean_mat=MI_mean_var_mat=MI_slope_mat=MI_slope_var_mat=MI_xyslope_mat=MI_xyslope_var_mat=matrix(NA, cycle_no, mi_no);
# exclusive imputation;
EX_mean_mat=EX_mean_var_mat=EX_slope_mat=EX_slope_var_mat=EX_xyslope_mat=EX_xyslope_var_mat=matrix(NA, cycle_no, mi_no);
# BUGS imputation;
BUGS_mean_mat=BUGS_mean_var_mat=BUGS_slope_mat=BUGS_slope_var_mat=BUGS_xyslope_mat=BUGS_xyslope_var_mat=matrix(NA, cycle_no, mi_no);

# record the logistic regression coefficient of missingness indicator on x,y,z;
logit_coef_mat=matrix(NA, cycle_no, 3);

cycle=1;

for (cycle in 1:cycle_no)

{
set.seed(cycle);

x=mu_x+sqrt(var_x)*rnorm(rowobs);

# error distribution
# could be normal, or could be other distributions;
error=sqrt(var_error)*rnorm(rowobs);

# generate the outcome variable;
y=beta0+beta1*x+error;


# generaite auxiliary variable z;
# not related to y nor related to the missingness mechanism;
# z=mu_z+sqrt(var_z)*rnorm(rowobs);

# related to Y but not related to the missingness mechanism;
error_z=sqrt(var_error)*rnorm(rowobs);

z=2+y+error_z;


# Z not related to Y but related to missingness mechanism;

# z=mu_z+sqrt(var_z)*rnorm(rowobs);

# alpha0=-1.7;
# alpha1=1
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*z)/(1+exp(alpha0+alpha1*z)));

# Z related to Y and also related to missingness mechainsim;
# error_z=sqrt(var_error)*rnorm(rowobs);
# z=2+y+error_z;
# alpha0=-1.7;
# alpha1=1
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*z)/(1+exp(alpha0+alpha1*z)));



# set up the missing data as MCAR on x;
# miss_indi=runif(n=rowobs)<0.40;

# set up the missing data as MAR on x;
# alpha0=-1.4;
# alpha1=1;
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*x)/(1+exp(alpha0+alpha1*x)));


# set up the missing data as MNAR on Y only;
# alpha0=-1.4;
# alpha1=-1;
# miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*y)/(1+exp(alpha0+alpha1*y)));

# set up the missing data as MNAR depends on both Y and X only;
alpha0=-1.4;
alpha1=-0.75;
alpha2=0.25;
miss_indi=rbinom(n=rowobs, size=1, prob=exp(alpha0+alpha1*y+alpha2*x)/(1+exp(alpha0+alpha1*y+alpha2*x)));



propensity_logistic=summary(glm(miss_indi~x+y+z, family=binomial(link="logit")));
logit_coef_mat[cycle,]=propensity_logistic$coeff[2:4,1];



obs_indi=1-miss_indi;

y_miss=y;
y_miss[miss_indi==1]=NA;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=length(y_miss[!is.na(y_miss)]);
mis_no=rowobs-obs_no;
y_obs=y[miss_indi==0];
x_obs=x[miss_indi==0];
x_mis=x[miss_indi==1];

# output the data;
# varname_cov=c("m[]", "y_miss[]", "x[]");

# write.table(cbind(miss_indi, y_miss, x), file=paste("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter13_univariate_nonignorable\\program\\example13.7.test1.dat",sep=""), row.name=FALSE, col.name=varname_cov);



# compelete data analysis;
# marginal mean;
BD_mean_vec[cycle]=mean(y);
BD_mean_var_vec[cycle]=var(y)/rowobs;


# regression analysis;
BD=summary(lm(y~x));

BD_slope_vec[cycle]=BD$coeff[2,1];
BD_slope_var_vec[cycle]=BD$coeff[2,2]^2;

# reverse regression analysis;
BD_xy=summary(lm(x~y));

BD_xyslope_vec[cycle]=BD_xy$coeff[2,1];
BD_xyslope_var_vec[cycle]=BD_xy$coeff[2,2]^2;


# different missing data methods;

# complete-case analysis;

# marginal mean;
CC_mean_vec[cycle]=mean(y_obs);
CC_mean_var_vec[cycle]=var(y_obs)/obs_no;


# regression analysis;
CC=summary(lm(y_obs~x_obs));

CC_slope_vec[cycle]=CC$coeff[2,1];
CC_slope_var_vec[cycle]=CC$coeff[2,2]^2;

# reverse regression analysis;
CC_xy=summary(lm(x_obs~y_obs));

CC_xyslope_vec[cycle]=CC_xy$coeff[2,1];
CC_xyslope_var_vec[cycle]=CC_xy$coeff[2,2]^2;

# missing data indicator method;
y_miss_0=y_miss;
y_miss_0[miss_indi==1]=0;

INDI_xy=summary(lm(x ~ obs_indi+y_miss_0));
INDI_xyslope_vec[cycle]=INDI_xy$coeff[3,1];
INDI_xyslope_var_vec[cycle]=INDI_xy$coeff[3,2]^2;



# now impute the missing y's to get estimates of the marginal;
y_completed=y_miss;

# single imputation method;
y_imputed=mice.impute.norm.predict(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);


y_completed[miss_seq]=y_imputed;

# completed-data analysis;

# marginal mean;
SI_mean_vec[cycle]=mean(y_completed);
SI_mean_var_vec[cycle]=var(y_completed)/rowobs;

# regression analysis;

IMP=summary(lm(y_completed~x));

SI_slope_vec[cycle]=IMP$coeff[2,1];
SI_slope_var_vec[cycle]=IMP$coeff[2,2]^2;

# reverse regression analysis;

IMP_xy=summary(lm(x~y_completed));

SI_xyslope_vec[cycle]=IMP_xy$coeff[2,1];
SI_xyslope_var_vec[cycle]=IMP_xy$coeff[2,2]^2;


# multiple imputation analysis;

# now impute the missing y's to get estimates of the marginal;
y_completed_MI=y_completed_EX=y_miss;

# the matrix holds the multiply imputed data;

y_completed_MI_mat=y_completed_EX_mat=matrix(NA, nrow=rowobs, ncol=mi_no);

for (i in 1:mi_no)
{
# inclusive imputation;

set.seed(i);

y_imputed_MI=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);


y_completed_MI[miss_seq]=y_imputed_MI;

y_completed_MI_mat[,i]=y_completed_MI;

# completed-data analysis;

# marginal mean;
MI_mean_mat[cycle,i]=mean(y_completed_MI);
MI_mean_var_mat[cycle,i]=var(y_completed_MI)/rowobs;

# regression analysis;

IMP=summary(lm(y_completed_MI~x));

MI_slope_mat[cycle,i]=IMP$coeff[2,1];
MI_slope_var_mat[cycle,i]=IMP$coeff[2,2]^2;

# reverse regression analysis;

IMP_xy=summary(lm(x~y_completed_MI));

MI_xyslope_mat[cycle,i]=IMP_xy$coeff[2,1];
MI_xyslope_var_mat[cycle,i]=IMP_xy$coeff[2,2]^2;

# exclusive imputation;


y_imputed_EX=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x,z));


y_completed_EX[miss_seq]=y_imputed_EX;

y_completed_EX_mat[,i]=y_completed_EX;

# completed-data analysis;

# marginal mean;
EX_mean_mat[cycle,i]=mean(y_completed_EX);
EX_mean_var_mat[cycle,i]=var(y_completed_EX)/rowobs;

# regression analysis;

IMP=summary(lm(y_completed_EX~x));

EX_slope_mat[cycle,i]=IMP$coeff[2,1];
EX_slope_var_mat[cycle,i]=IMP$coeff[2,2]^2;

# reverse regression analysis;

IMP_xy=summary(lm(x~y_completed_EX));

EX_xyslope_mat[cycle,i]=IMP_xy$coeff[2,1];
EX_xyslope_var_mat[cycle,i]=IMP_xy$coeff[2,2]^2;





}


# BUGS imputation;


# apply missing cases to the data;
M=rowobs;
gibbs_no=20000;


# for model, response only depends on y;

# process_stage2_data_simulate=list(M=M, y_miss=y_miss, x=x, m=miss_indi);
# process_stage2_init= function(){list(beta0=0, beta1=0, tau=1, alpha0=0, alpha1=0)};
# process_stage2_parameters=c("y_impute");

# process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file,
#    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=gibbs_no/(2*mi_no), bugs.seed=cycle,
#    bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
#    working.directory=NULL, clearWD=TRUE);


# for model1, response depends on both y and x;


process_stage2_data_simulate=list(M=M, y_miss=y_miss, x=x, m=miss_indi);
process_stage2_init= function(){list(beta0=-2, beta1=1, beta2=0, sigma=1, alpha0=0, alpha1=0.5)};
process_stage2_parameters=c("y_impute");

process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model2.file,
    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=gibbs_no/(2*mi_no), bugs.seed=cycle,
    bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);

attach.bugs(process_stage2.sim);
rm(process_stage2.sim);

y_completed_BUGS=y_miss;

for (i in 1:mi_no)

{

y_completed_BUGS[miss_indi==1]=y_impute[i,];

# completed-data analysis;

# marginal mean;
BUGS_mean_mat[cycle,i]=mean(y_completed_BUGS);
BUGS_mean_var_mat[cycle,i]=var(y_completed_BUGS)/rowobs;

# regression analysis;

IMP=summary(lm(y_completed_BUGS~x));

BUGS_slope_mat[cycle,i]=IMP$coeff[2,1];
BUGS_slope_var_mat[cycle,i]=IMP$coeff[2,2]^2;

# reverse regression analysis;

IMP_xy=summary(lm(x~y_completed_BUGS));

BUGS_xyslope_mat[cycle,i]=IMP_xy$coeff[2,1];
BUGS_xyslope_var_mat[cycle,i]=IMP_xy$coeff[2,2]^2;




}



cat("the cycle is", cycle, "\n");


}

# the end of multiple imputation;

# plot the distribution of the data;

# check the performance of the estimates;
# population quantity

#############################################################################
# marginal mean estimand;
true_mean=mean(BD_mean_vec);
true_mean;


BD=performance(BD_mean_vec, BD_mean_var_vec, true_mean);
BD;
CC=performance(CC_mean_vec, CC_mean_var_vec, true_mean);
CC;

MI=mi_performance(MI_mean_mat, MI_mean_var_mat, true_mean, rowobs);
MI;

EX=mi_performance(EX_mean_mat, EX_mean_var_mat, true_mean, rowobs);
EX;


BUGS=mi_performance(BUGS_mean_mat, BUGS_mean_var_mat, true_mean, rowobs);
BUGS;

hist(rowMeans(BUGS_mean_mat));

hist(BUGS_mean_mat[1,]);


##########################################################################
# slope estimand;

# complete-data analysis;

true_slope=mean(BD_slope_vec);
true_slope;

BD=performance(BD_slope_vec, BD_slope_var_vec, true_slope);
BD;
CC=performance(CC_slope_vec, CC_slope_var_vec, true_slope);
CC;


MI=mi_performance(MI_slope_mat, MI_slope_var_mat, true_slope, rowobs);
MI;

EX=mi_performance(EX_slope_mat, EX_slope_var_mat, true_slope, rowobs);
EX;


BUGS=mi_performance(BUGS_slope_mat, BUGS_slope_var_mat, true_slope, rowobs);
BUGS;



#######################################################################
# reverse slope;

# complete-data analysis;

true_xyslope=mean(BD_xyslope_vec);
true_xyslope;

BD=performance(BD_xyslope_vec, BD_xyslope_var_vec, true_xyslope);
BD;
CC=performance(CC_xyslope_vec, CC_xyslope_var_vec, true_xyslope);
CC;


MI=mi_performance(MI_xyslope_mat, MI_xyslope_var_mat, true_xyslope, rowobs);
MI;

EX=mi_performance(EX_xyslope_mat, EX_xyslope_var_mat, true_xyslope, rowobs);
EX;

BUGS=mi_performance(BUGS_xyslope_mat, BUGS_xyslope_var_mat, true_xyslope, rowobs);
BUGS;

