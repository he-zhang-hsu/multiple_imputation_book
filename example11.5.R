library(R2WinBUGS);
library(MASS);
library(mice);
library(norm);
library(mix);
# library(HI);
options(digits=4);
rm(list=ls());


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





# setwd('\\\\cdc.gov\\private\\L728\\book\\example8.4')
# get the winbugs model file;

model.file <- system.file(package="R2WinBUGS", "model", "logitmodel_mismeasure_bivariate_forR.txt")
# Let's take a look:
file.show(model.file)


# generate the logistic outcome;
beta0=-1;
beta1=0.5;
beta2=0.5;

# generate bivariate normal x variables;
mu=rep(0,2);
Sigma=matrix(c(1,0.5,0.5,1), nrow=2,ncol=2);


rowobs=1000;
cycle_no=1000;
mi_no=80;


# x1_x2=mvrnorm(n = rowobs, mu=mu, Sigma=Sigma);

# x1=x1_x2[,1];
# x2=x1_x2[,2];


# logistic outcome;
# y2=rbinom(rowobs,size=1,prob=exp(beta0+beta1*x1+beta2*x2)/(1+exp(beta0+beta1*x1+beta2*x2)));

c0=-0.2;
c1=1.2;
c2=1;

# y is observed after applying some error;

# varname_cov=c("y2[]", "x1[]", "x1_o[]","x2[]");

# write.table(cbind(y2, x1, x1_o, x2), file=paste("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter11_measurementerror\\program\\example11.5.test.dat",sep=""), row.name=FALSE, col.name=varname_cov);




# matrices holding the parameter estimates;
# complete-data inferences;
BD_mean_vec=BD_mean_var_vec=CC_mean_vec=CC_mean_var_vec=SI_mean_vec=SI_mean_var_vec=rep(NA, cycle_no);

bugs_mean_mat=bugs_mean_var_mat=glm_mean_mat=glm_mean_var_mat=matrix(NA, cycle_no, mi_no);

BD_slope1_vec=BD_slope1_var_vec=rep(NA, cycle_no);
CC_slope1_vec=CC_slope1_var_vec=rep(NA, cycle_no);
SI_slope1_vec=SI_slope1_var_vec=rep(NA, cycle_no);
bugs_slope1_mat=bugs_slope1_var_mat=glm_slope1_mat=glm_slope1_var_mat=matrix(NA, cycle_no, mi_no);


BD_slope2_vec=BD_slope2_var_vec=rep(NA, cycle_no);
CC_slope2_vec=CC_slope2_var_vec=rep(NA, cycle_no);
SI_slope2_vec=SI_slope2_var_vec=rep(NA, cycle_no);
bugs_slope2_mat=bugs_slope2_var_mat=glm_slope2_mat=glm_slope2_var_mat=matrix(NA, cycle_no, mi_no);




# generate data for example 8.4
for (cycle in 1:cycle_no)

{
set.seed(cycle);

x1_x2=mvrnorm(n = rowobs, mu=mu, Sigma=Sigma);

x1_BD=x1=x1_x2[,1];
x2=x1_x2[,2];

# logistic outcome;
y2=rbinom(rowobs,size=1,prob=exp(beta0+beta1*x1+beta2*x2)/(1+exp(beta0+beta1*x1+beta2*x2)));

# apply the mismeasurement to x1;
x1_o=c0+c1*x1+c2*rnorm(rowobs);

# missing completely at random;
# set up the missing data as MCAR on x;
miss_indi=runif(n=rowobs)<0.80;
x1[miss_indi==1]=NA;
mis_no=sum(miss_indi);
obs_no=rowobs-mis_no;

# before deletion analysis;
# compelete data analysis;
# sample mean;
BD_mean_vec[cycle]=mean(x1_BD);
BD_mean_var_vec[cycle]=var(x1_BD)/rowobs;

# logistic regression coefficient for the slope of x1;
BD_logistic=summary(glm(y2~x1_BD+x2, family=binomial(link="logit")));
BD_slope1_vec[cycle]=BD_logistic$coeff[2,1];
BD_slope1_var_vec[cycle]=BD_logistic$coeff[2,2]^2;

BD_slope2_vec[cycle]=BD_logistic$coeff[3,1];
BD_slope2_var_vec[cycle]=BD_logistic$coeff[3,2]^2;



# complete-case analysis using the validation study;
# mean estimands;
CC_mean_vec[cycle]=mean(x1, na.rm="T");
CC_mean_var_vec[cycle]=var(x1, na.rm="T")/obs_no;

# logistic regression coefficient for the slope of x1;
CC_logistic=summary(glm(y2~x1+x2, family=binomial(link="logit")));
CC_slope1_vec[cycle]=CC_logistic$coeff[2,1];
CC_slope1_var_vec[cycle]=CC_logistic$coeff[2,2]^2;

CC_slope2_vec[cycle]=CC_logistic$coeff[3,1];
CC_slope2_var_vec[cycle]=CC_logistic$coeff[3,2]^2;


# using the wrong x1;
# mean estimands;
SI_mean_vec[cycle]=mean(x1_o);
SI_mean_var_vec[cycle]=var(x1_o)/rowobs;

# logistic regression coefficient for the slope of x1;
SI_logistic=summary(glm(y2~x1_o+x2, family=binomial(link="logit")));
SI_slope1_vec[cycle]=SI_logistic$coeff[2,1];
SI_slope1_var_vec[cycle]=SI_logistic$coeff[2,2]^2;

SI_slope2_vec[cycle]=SI_logistic$coeff[3,1];
SI_slope2_var_vec[cycle]=SI_logistic$coeff[3,2]^2;



# winbugs imputation;

# apply missing cases to the data;
M=rowobs;

process_stage2_data_simulate=list(M=M, y2=y2, x1=x1, x2=x2, x1_o=x1_o);
process_stage2_init= function(){list(beta0=-1, beta1=0.5, beta2=0.5, mu=0, tau=1, alpha0_12=0, alpha1_12=0, alpha0_m1=0, alpha1_m1=0, psi_12=1, psi_m1=1)};
# process_stage2_parameters=c("mu", "beta0", "beta1", "beta2", "mu", "alpha0", "alpha1", "tau", "psi", "y2", "x1", "x2");
process_stage2_parameters=c("x1_impute");


gibbs_no=20000;


process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file,
    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=gibbs_no/(2*mi_no), bugs.seed=cycle,
    bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);





#process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file="C:/Users/WDQ7/reliability_test_forR.txt",
#    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=1, bugs.seed=197789,
#    bugs.directory="C:/Program Files (x86)/OpenBUGS/OpenBUGS321", summary.only=FALSE, debug=FALSE,
#    working.directory=NULL, clearWD=TRUE);


attach.bugs(process_stage2.sim);
rm(process_stage2.sim);

x1_completed_bugs=x1;

for (i in 1:mi_no)

{
x1_completed_bugs[miss_indi==1]=x1_impute[i,];

# marginal means;

bugs_mean_mat[cycle,i]=mean(x1_completed_bugs);
bugs_mean_var_mat[cycle,i]=var(x1_completed_bugs)/rowobs;

# logistic regression coefficient for x1 slope;
bugs_logistic=summary(glm(y2~x1_completed_bugs+x2, family=binomial(link="logit")));
bugs_slope1_mat[cycle,i]=bugs_logistic$coeff[2,1];
bugs_slope1_var_mat[cycle,i]=bugs_logistic$coeff[2,2]^2;

bugs_slope2_mat[cycle,i]=bugs_logistic$coeff[3,1];
bugs_slope2_var_mat[cycle,i]=bugs_logistic$coeff[3,2]^2;


}

# linear model imputation;
x1_completed_glm=x1;
x_new=cbind(x1_o, y2, x2);

for (i in 1:mi_no)
{

# proper model imputation;
set.seed(i);

x1_imputed_glm=mice.impute.norm(x1, ry=as.logical(1-miss_indi), seed=i, x=x_new);
x1_completed_glm[miss_indi==1]=x1_imputed_glm;


# multiple imputation analysis;

# marginal means;
glm_mean_mat[cycle,i]=mean(x1_completed_glm);;
glm_mean_var_mat[cycle,i]=var(x1_completed_glm)/rowobs;

# logistic regression coefficient for x1 slope;
glm_logistic=summary(glm(y2~x1_completed_glm+x2, family=binomial(link="logit")));
glm_slope1_mat[cycle,i]=glm_logistic$coeff[2,1];
glm_slope1_var_mat[cycle,i]=glm_logistic$coeff[2,2]^2;

glm_slope2_mat[cycle,i]=glm_logistic$coeff[3,1];
glm_slope2_var_mat[cycle,i]=glm_logistic$coeff[3,2]^2;



}

cat("the cycle is", cycle, "\n");

}

# the end of mulitple imputation;
#############################################################################
# marginal mean estimand;

true_mean=mean(BD_mean_vec);
true_mean;


BD=performance(BD_mean_vec, BD_mean_var_vec, true_mean);
BD;
CC=performance(CC_mean_vec, CC_mean_var_vec, true_mean);
CC;
SI=performance(SI_mean_vec, SI_mean_var_vec, true_mean);
SI;

bugs=mi_performance(bugs_mean_mat, bugs_mean_var_mat, true_mean, rowobs);
bugs;



glm=mi_performance(glm_mean_mat, glm_mean_var_mat, true_mean, rowobs);
glm;




##########################################################################
# slope estimand;

# complete-data analysis;

true_slope1=mean(BD_slope1_vec);
true_slope1;


BD=performance(BD_slope1_vec, BD_slope1_var_vec, true_slope1);
BD;
CC=performance(CC_slope1_vec, CC_slope1_var_vec, true_slope1);
CC;
SI=performance(SI_slope1_vec, SI_slope1_var_vec, true_slope1);
SI;

bugs=mi_performance(bugs_slope1_mat, bugs_slope1_var_mat, true_slope1, rowobs);
bugs;

glm=mi_performance(glm_slope1_mat, glm_slope1_var_mat, true_slope1, rowobs);
glm;

true_slope2=mean(BD_slope2_vec);
true_slope2;


BD=performance(BD_slope2_vec, BD_slope2_var_vec, true_slope2);
BD;
CC=performance(CC_slope2_vec, CC_slope2_var_vec, true_slope2);
CC;
SI=performance(SI_slope2_vec, SI_slope2_var_vec, true_slope2);
SI;

bugs=mi_performance(bugs_slope2_mat, bugs_slope2_var_mat, true_slope2, rowobs);
bugs;

glm=mi_performance(glm_slope2_mat, glm_slope2_var_mat, true_slope2, rowobs);
glm;

