# Example 6.2 
# Multiple imputation using trivariate t-distribution;

rm(list=ls());

# library(car);
# library(mvtnorm);
library(mice);
library(norm);
library(MASS);
library(miscF);
library(tmvtnorm);
library(TruncatedNormal);
library(mvtnorm);
library(BRugs);
library(crch);

options(digits=4);

# read the data;
data=read.table(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter6_multivariate_jointmodeling\\program\\colorado_ami_chf_pb.txt", col.names=c("ami", "chf", "cap"), na.strings=".");

# Logit transformation;
trans_data=log(data/(1-data));

total_no=nrow(data);
rowobs=total_no;

missindi=cbind(is.na(data$ami),is.na(data$chf),is.na(data$cap));
obsindi=1-missindi;
obs_no=total_no-colSums(missindi);
# cc_no is the number of complete-cases;
cc_no=sum(obsindi[,1]*obsindi[,2]*obsindi[,3]);

hosp_performance=as.matrix(data);

trans_hosp_performance=as.matrix(trans_data);

# figure out the missingness pattern;
A=trans_hosp_performance;

# The following code prepares missingness pattern matrices

# scenario 1;
# xmis_ymis_zmis is the matrix 
# where all three variables are missing;
# n_xmis_ymis_zmis is the number of cases
# where all three variables are missing;
xmis_ymis_zmis=A[is.na(A[,1])==TRUE & is.na(A[,2])==TRUE & is.na(A[,3])==TRUE,];
n_xmis_ymis_zmis=nrow(xmis_ymis_zmis);

# scenario 2;

xobs_yobs_zobs=A[is.na(A[,1])==FALSE & is.na(A[,2])==FALSE & is.na(A[,3])==FALSE,];
n_xobs_yobs_zobs=nrow(xobs_yobs_zobs);

# scenario 3;
xobs_ymis_zobs=A[is.na(A[,1])==FALSE & is.na(A[,2])==TRUE & is.na(A[,3])==FALSE,];
n_xobs_ymis_zobs=nrow(xobs_ymis_zobs);

xobs_ymis_zobs_base=xobs_ymis_zobs[,c(1,3)];
xobs_ymis_zobs_base;

# scenario 4;

xobs_yobs_zmis=A[is.na(A[,1])==FALSE & is.na(A[,2])==FALSE & is.na(A[,3])==TRUE,];
n_xobs_yobs_zmis=nrow(xobs_yobs_zmis);

xobs_yobs_zmis_base=xobs_yobs_zmis[,c(1,2)];
xobs_yobs_zmis_base;

# scenario 5;
xmis_yobs_zobs=A[is.na(A[,1])==TRUE & is.na(A[,2])==FALSE & is.na(A[,3])==FALSE,];
n_xmis_yobs_zobs=nrow(xmis_yobs_zobs);

xmis_yobs_zobs_base=xmis_yobs_zobs[,c(2,3)];

xmis_yobs_zobs_base;

# scenario 6 (no cases);
xobs_ymis_zmis=A[is.na(A[,1])==FALSE & is.na(A[,2])==TRUE & is.na(A[,3])==TRUE,];
n_xobs_ymis_zmis=nrow(xobs_ymis_zmis);

xobs_ymis_zmis_base=xobs_ymis_zmis[,1];

xobs_ymis_zmis_base;

# scenario 7;
xmis_yobs_zmis=A[is.na(A[,1])==TRUE & is.na(A[,2])==FALSE & is.na(A[,3])==TRUE,];
n_xmis_yobs_zmis=nrow(xmis_yobs_zmis);

xmis_yobs_zmis_base=xmis_yobs_zmis[,2];
xmis_yobs_zmis_base;

# scenario 8;
xmis_ymis_zobs=A[is.na(A[,1])==TRUE & is.na(A[,2])==TRUE & is.na(A[,3])==FALSE,];
n_xmis_ymis_zobs=nrow(xmis_ymis_zobs);

xmis_ymis_zobs_base=xmis_ymis_zobs[,3];
xmis_ymis_zobs_base;

n_xobs_yobs_zobs+n_xobs_ymis_zobs+n_xobs_yobs_zmis+n_xmis_yobs_zobs+n_xobs_ymis_zmis+
n_xmis_yobs_zmis+n_xmis_ymis_zobs+n_xmis_ymis_zmis;


# try fit the trivariate t-model for the observed cases;
# First obtain the initial values for the mean vector
# and covariance matrix using multivariate normal model;
trans_s=prelim.norm(trans_hosp_performance);

# obtain the initial value for thetahat;
trans_thetahat=em.norm(trans_s)

para_initial=getparam.norm(trans_s, trans_thetahat, corr=FALSE);

mu_initial=para_initial$mu;

Sigma_initial=para_initial$sigma;

# set the initial value
# for the degrees of freedom as 30
nu_initial=30;

# Fit the tri-variate t model to complete cases;
result_cc=mvt.mcmc(xobs_yobs_zobs, prior.lower.v=1, prior.upper.v=200, initial.v=nu_initial, 
initial.Sigma=Sigma_initial) 

# Mean vector
Mu_cc=colMeans(result_cc$Mu.save)
Mu_cc
# Scale matrix 
Sigma_cc=apply(result_cc$Sigma.save, c(1,2), mean)
Sigma_cc 
# Degrees of freedom
nu_cc=mean(result_cc$v.save)
nu_cc

Sigma_cc*nu_cc/(nu_cc-2);

# For example, plot the 
# Posterior samples of the 
# degrees of freedom
hist(result_cc$v.save);
plot(result_cc$v.save);

# initial imputation;
Impute_cc=matrix(rep(Mu_cc, rowobs), nrow=rowobs, byrow=T)+rmvt(n=rowobs, sigma = Sigma_cc, df = nu_cc);

temp_impute=A;

temp_impute[is.na(A)==TRUE]=Impute_cc[is.na(A)==TRUE];

# Data augmentation;

iter_no=1000;
parameter_draw_matrix=matrix(NA, iter_no, 10);
mcmc_burnin=100;
mcmc_iter=100;

impute_data_array=array(data=NA, dim=c(50, rowobs, 3));

# begin the iteration
# This algorithm takes a very long time

for (i in 1:iter_no)
{

# Estimate the mean and covariance matrix only using data with at least one observed element;

temp_impute_obs=temp_impute[is.na(A[,1])==FALSE | is.na(A[,2])==FALSE | is.na(A[,3])==FALSE,];

# parameter_draw=mvt.mcmc(temp_impute_obs, prior.lower.v=1, prior.upper.v=200, initial.v=nu_cc, initial.Sigma=Sigma_cc, nmcmc=mcmc_iter, nburn=mcmc_burnin);

# The Posterior step of DA;
parameter_draw=mvt.mcmc(temp_impute_obs, prior.lower.v=1, prior.upper.v=200, nmcmc=mcmc_iter, nburn=mcmc_burnin);

Mu_draw=parameter_draw$Mu.save[mcmc_iter,];
Sigma_draw=parameter_draw$Sigma.save[,,mcmc_iter];
v_draw=parameter_draw$v.save[mcmc_iter];

parameter_draw_matrix[i, 1:3]=Mu_draw;
parameter_draw_matrix[i, 4:9]=
c(Sigma_draw[1,1], Sigma_draw[1,2], Sigma_draw[1,3], Sigma_draw[2,2], Sigma_draw[2,3],
Sigma_draw[3,3]);
parameter_draw_matrix[i,10]=v_draw;

# obtain the conditional parameters;
# Disect the trivariate-t model
# to conditinal t-distributions;

# scenario 3 ymis, xobs, zobs;

ymis_xobs_zobs_predbase=cbind(1, t(t(xobs_ymis_zobs_base)-Mu_draw[c(1,3)]));

ymis_xobs_zobs_slope=c(Mu_draw[2], matrix(c(Sigma_draw[2,1], Sigma_draw[2,3]),1,2)%*%solve(matrix(c(Sigma_draw[1,1], Sigma_draw[1,3], Sigma_draw[3,1], Sigma_draw[3,3]),2,2)));

mu_ymis_xobs_zobs =ymis_xobs_zobs_predbase %*% ymis_xobs_zobs_slope;
 
sigma2_ymis_xobs_zobs=Sigma_draw[2,2]-matrix(c(Sigma_draw[2,1], Sigma_draw[2,3]),1,2)%*%solve(matrix(c(Sigma_draw[1,1], Sigma_draw[1,3], Sigma_draw[3,1], Sigma_draw[3,3]),2,2))%*%c(Sigma_draw[2,1], Sigma_draw[2,3]);

v_ymis_xobs_zobs=v_draw+2;

d_ymis_xobs_zobs=diag(t(t(xobs_ymis_zobs_base)-Mu_draw[c(1,3)])%*%solve(matrix(c(Sigma_draw[1,1], Sigma_draw[1,3], Sigma_draw[3,1], Sigma_draw[3,3]),2,2))%*%(t(xobs_ymis_zobs_base)-Mu_draw[c(1,3)]));


# scenario 4 xobs yobs zmis;

zmis_xobs_yobs_predbase=cbind(1, t(t(xobs_yobs_zmis_base)-Mu_draw[c(1,2)]));

zmis_xobs_yobs_slope=c(Mu_draw[3], matrix(c(Sigma_draw[3,1], Sigma_draw[3,2]),1,2)%*%solve(matrix(c(Sigma_draw[1,1], Sigma_draw[1,2], Sigma_draw[2,1], Sigma_draw[2,2]),2,2)));

mu_zmis_xobs_yobs =zmis_xobs_yobs_predbase %*% zmis_xobs_yobs_slope;
 
sigma2_zmis_xobs_yobs=Sigma_draw[3,3]-matrix(c(Sigma_draw[3,1], Sigma_draw[3,2]),1,2)%*%solve(matrix(c(Sigma_draw[1,1], Sigma_draw[1,2], Sigma_draw[2,1], Sigma_draw[2,2]),2,2))%*%c(Sigma_draw[3,1], Sigma_draw[3,2]);

v_zmis_xobs_yobs=v_draw+2;

d_zmis_xobs_yobs=diag(t(t(xobs_yobs_zmis_base)-Mu_draw[c(1,2)])%*%solve(matrix(c(Sigma_draw[1,1], Sigma_draw[1,2], Sigma_draw[2,1], Sigma_draw[2,2]),2,2))%*%(t(xobs_yobs_zmis_base)-Mu_draw[c(1,2)]));

# scenario 5 xmis yobs zobs;

xmis_yobs_zobs_predbase=cbind(1, t(t(xmis_yobs_zobs_base)-Mu_draw[c(2,3)]));

xmis_yobs_zobs_slope=c(Mu_draw[1], matrix(c(Sigma_draw[1,2], Sigma_draw[1,3]),1,2)%*%solve(matrix(c(Sigma_draw[2,2], Sigma_draw[2,3], Sigma_draw[3,2], Sigma_draw[3,3]),2,2)));

mu_xmis_yobs_zobs =xmis_yobs_zobs_predbase %*% xmis_yobs_zobs_slope;
 
sigma2_xmis_yobs_zobs=Sigma_draw[1,1]-matrix(c(Sigma_draw[1,2], Sigma_draw[1,3]),1,2)%*%solve(matrix(c(Sigma_draw[2,2], Sigma_draw[2,3], Sigma_draw[3,2], Sigma_draw[3,3]),2,2))%*%c(Sigma_draw[1,2], Sigma_draw[1,3]);

v_xmis_yobs_zobs=v_draw+2;

d_xmis_yobs_zobs=diag(t(t(xmis_yobs_zobs_base)-Mu_draw[c(2,3)])%*%solve(matrix(c(Sigma_draw[2,2], Sigma_draw[2,3], Sigma_draw[3,2], Sigma_draw[3,3]),2,2))%*%(t(xmis_yobs_zobs_base)-Mu_draw[c(2,3)]));

# no scenario 6;

# scenario 7 xmis, yobs, zmis;
mu_xmis_zmis_yobs=matrix(rep(Mu_draw[c(1,3)], n_xmis_yobs_zmis), nrow=n_xmis_yobs_zmis, byrow=T)+(xmis_yobs_zmis_base-Mu_draw[2])%*%t(c(Sigma_draw[1,2], Sigma_draw[3,2]))/Sigma_draw[2,2];

sigma2_xmis_zmis_yobs=matrix(c(Sigma_draw[1,1], Sigma_draw[1,3], Sigma_draw[3,1], Sigma_draw[3,3]), nrow=2)-c(Sigma_draw[1,2], Sigma_draw[3,2])%*%t(c(Sigma_draw[1,2], Sigma_draw[3,2]))/Sigma_draw[2,2]

v_xmis_zmis_yobs=v_draw+1;

d_xmis_zmis_yobs=(xmis_yobs_zmis_base-Mu_draw[2])^2/Sigma_draw[2,2];

# scenario 8 xmis ymis zobs;
mu_xmis_ymis_zobs=matrix(rep(Mu_draw[c(1,2)], n_xmis_ymis_zobs), nrow=n_xmis_ymis_zobs, byrow=T)+(xmis_ymis_zobs_base-Mu_draw[3])%*%t(c(Sigma_draw[1,3], Sigma_draw[2,3]))/Sigma_draw[3,3];

sigma2_xmis_ymis_zobs=matrix(c(Sigma_draw[1,1], Sigma_draw[1,2], Sigma_draw[2,1], Sigma_draw[2,2]), nrow=2)-c(Sigma_draw[1,3], Sigma_draw[2,3])%*%t(c(Sigma_draw[1,3], Sigma_draw[2,3]))/Sigma_draw[3,3]

v_xmis_ymis_zobs=v_draw+1;

d_xmis_ymis_zobs=(xmis_ymis_zobs_base-Mu_draw[3])^2/Sigma_draw[3,3];


# imputation;
# The I-step of the DA algorithm;

# scenario 1, xmis, ymis, zmis;
ximpute_yimpute_zimpute=matrix(rep(Mu_draw, n_xmis_ymis_zmis), nrow=n_xmis_ymis_zmis, byrow=T)+rmvt(n=n_xmis_ymis_zmis, sigma = Sigma_draw, df = v_draw);

# scenario 3, ymis, xobs, zobs;
xobs_yimpute_zobs=mu_ymis_xobs_zobs+sqrt(sigma2_ymis_xobs_zobs*(v_draw+d_ymis_xobs_zobs)/v_ymis_xobs_zobs)*rt(n=n_xobs_ymis_zobs, df=v_ymis_xobs_zobs);

# scenario 4 xobs yobs zmis;
xobs_yobs_zimpute=mu_zmis_xobs_yobs+sqrt(sigma2_zmis_xobs_yobs*(v_draw+d_zmis_xobs_yobs)/v_zmis_xobs_yobs)*rt(n=n_xobs_yobs_zmis, df=v_zmis_xobs_yobs);

# scenario 5 xmis yobs zobs;
ximpute_yobs_zobs=mu_xmis_yobs_zobs+sqrt(sigma2_xmis_yobs_zobs*(v_draw+d_xmis_yobs_zobs)/v_xmis_yobs_zobs)*rt(n=n_xmis_yobs_zobs, df=v_xmis_yobs_zobs);

# no scenario 6;

# scenario 7 xmis, yobs, zmis;
ximpute_yobs_zimpute=mu_xmis_zmis_yobs+sqrt(v_draw+d_xmis_zmis_yobs/v_xmis_zmis_yobs)*rmvt(n=n_xmis_yobs_zmis, sigma = sigma2_xmis_zmis_yobs, df = v_xmis_zmis_yobs);

# scenario 8 xmis, ymis, zobs;

ximpute_yimpute_zobs=mu_xmis_ymis_zobs+sqrt(v_draw+d_xmis_ymis_zobs/v_xmis_ymis_zobs)*rmvt(n=n_xmis_ymis_zobs, sigma = sigma2_xmis_ymis_zobs, df= v_xmis_ymis_zobs);

# scenario 1;
temp_impute[is.na(A[,1])==TRUE & is.na(A[,2])==TRUE & is.na(A[,3])==TRUE, ]=ximpute_yimpute_zimpute;

# scenario 3, ymis, xobs, zobs;
temp_impute[is.na(A[,1])==FALSE & is.na(A[,2])==TRUE & is.na(A[,3])==FALSE, 2]=xobs_yimpute_zobs;

# scenario 4 xobs yobs zmis;

temp_impute[is.na(A[,1])==FALSE & is.na(A[,2])==FALSE & is.na(A[,3])==TRUE, 3]=xobs_yobs_zimpute;


# scenario 5 xmis yobs zobs;
temp_impute[is.na(A[,1])==TRUE & is.na(A[,2])==FALSE & is.na(A[,3])==FALSE, 1]=ximpute_yobs_zobs;



# no scenario 6;

# scenario 7 xmis, yobs, zmis;

temp_impute[is.na(A[,1])==TRUE & is.na(A[,2])==FALSE & is.na(A[,3])==TRUE, c(1,3)]=ximpute_yobs_zimpute;

# scenario 8 xmis, ymis, zobs;

temp_impute[is.na(A[,1])==TRUE & is.na(A[,2])==TRUE & is.na(A[,3])==FALSE, c(1,2)]=ximpute_yimpute_zobs;

# output the dataset;

if (i%%20==0) impute_data_array[i/20,,]=temp_impute;

  
cat("the iteration is", i, "\n");
cat("the parameter draw is", parameter_draw_matrix[i,], "\n");
}

# Fig. 6.4 

plot(parameter_draw_matrix[,10], type="l", xlab="iteration", ylab="t degrees of freedom");

# Posterior means of model parameters
colMeans(parameter_draw_matrix);

hist(parameter_draw_matrix[,10], main="", xlab="t degrees of freedom");

# Multiple imputation estimates
# Table 6.1 
# TMI method;

mi_no=50;
trans_est_mat=trans_var_mat=matrix(NA, nrow=14, ncol=mi_no);

for(i in 1:mi_no)
{

trans_hosp_imp=exp(impute_data_array[i,,])/(1+exp(impute_data_array[i,,]));
# completed-data estimates;
trans_est_mat[1:3,i]=colMeans(trans_hosp_imp);
trans_var_mat[1:3,i]=diag(var(trans_hosp_imp))/rowobs;
trans_est_mat[4,i]=mean(trans_hosp_imp[,1]>0.9);
trans_var_mat[4,i]=var(trans_hosp_imp[,1]>0.9)/rowobs;
trans_est_mat[5,i]=mean(trans_hosp_imp[,1]>0.95);
trans_var_mat[5,i]=var(trans_hosp_imp[,1]>0.95)/rowobs;
trans_est_mat[6,i]=mean(trans_hosp_imp[,1]>0.99);
trans_var_mat[6,i]=var(trans_hosp_imp[,1]>0.99)/rowobs;
trans_est_mat[7,i]=mean(trans_hosp_imp[,2]>0.9);
trans_var_mat[7,i]=var(trans_hosp_imp[,2]>0.9)/rowobs;
trans_est_mat[8,i]=mean(trans_hosp_imp[,2]>0.95);
trans_var_mat[8,i]=var(trans_hosp_imp[,2]>0.95)/rowobs;
trans_est_mat[9,i]=mean(trans_hosp_imp[,2]>0.99);
trans_var_mat[9,i]=var(trans_hosp_imp[,2]>0.99)/rowobs;
trans_est_mat[10,i]=mean(trans_hosp_imp[,3]>0.9);
trans_var_mat[10,i]=var(trans_hosp_imp[,3]>0.9)/rowobs;
trans_est_mat[11,i]=mean(trans_hosp_imp[,3]>0.95);
trans_var_mat[11,i]=var(trans_hosp_imp[,3]>0.95)/rowobs;
trans_est_mat[12,i]=mean(trans_hosp_imp[,3]>0.99);
trans_var_mat[12,i]=var(trans_hosp_imp[,3]>0.99)/rowobs;
trans_est_mat[13,i]=mean(trans_hosp_imp[,1]>0.9 & trans_hosp_imp[,2]> 0.9 & trans_hosp_imp[,3]>0.9);
trans_var_mat[13,i]=var(trans_hosp_imp[,1]>0.9 & trans_hosp_imp[,2]> 0.9 & trans_hosp_imp[,3]>0.9)/rowobs;
trans_est_mat[14,i]=mean(trans_hosp_imp[,1]>0.95 & trans_hosp_imp[,2]> 0.95 & trans_hosp_imp[,3]>0.95);
trans_var_mat[14,i]=var(trans_hosp_imp[,1]>0.95 & trans_hosp_imp[,2]> 0.95 & trans_hosp_imp[,3]>0.95)/rowobs;

}

# summarize the inference;
trans_summary_matrix=matrix(NA, 14,3);
for (i in 1:14)
{
marginal_summary=pool.scalar(trans_est_mat[i,], trans_var_mat[i,], n=rowobs);
marginal_mi_mean=marginal_summary$qbar;
marginal_mi_var=marginal_summary$t;
marginal_mi_df=marginal_summary$df;
marginal_mi_f=marginal_summary$f;

trans_summary_matrix[i,1]=marginal_mi_mean;
trans_summary_matrix[i,2]=marginal_mi_mean-1.96*sqrt(marginal_mi_var);
trans_summary_matrix[i,3]=marginal_mi_mean+1.96*sqrt(marginal_mi_var);
}

trans_summary_matrix;
