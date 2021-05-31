# Example 6.2

rm(list=ls());
# library(car);
# library(mvtnorm);
library(mice);
library(norm);
# library(HI);
options(digits=4);

# read the data

data=read.table(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter6_multivariate_jointmodeling\\program\\colorado_ami_chf_pb.txt", col.names=c("ami", "chf", "cap"), na.strings=".");


# trans_data=asin(sqrt(data));
# Logit transformation;
trans_data=log(data/(1-data));

# Some descriptive statistics
summary(data$ami);
summary(data$chf);
summary(data$cap);

total_no=nrow(data);
rowobs=total_no;

# missingness rates;

sum(is.na(data$ami))/total_no;
sum(is.na(data$chf))/total_no;
sum(is.na(data$cap))/total_no;

missindi=cbind(is.na(data$ami),is.na(data$chf),is.na(data$cap));
obsindi=1-missindi;
obs_no=total_no-colSums(missindi);
# cc_no is the number of complete-cases;
cc_no=sum(obsindi[,1]*obsindi[,2]*obsindi[,3]);

hosp_performance=as.matrix(data);

# plots of Example 6.2;
# Fig. 6.1 Top row
# Histograms on the original scale;
hist(hosp_performance[,1], xlab="AMI", main="");
hist(hosp_performance[,2], xlab="CHF", main="");
hist(hosp_performance[,3], xlab="CAP", main="");

# plots on the logit transformation;
trans_hosp_performance=as.matrix(trans_data);

# Fig. 6.1 middle and bottom row;
# Histograms and QQ plots on the transformed scale;
hist(trans_hosp_performance[,1], xlab="Logit of AMI", main="");
hist(trans_hosp_performance[,2], xlab="Logit of CHF", main="");
hist(trans_hosp_performance[,3], xlab="Logit of CAP", main="");

qqnorm(trans_hosp_performance[,1]);
qqline(trans_hosp_performance[,1]);

qqnorm(trans_hosp_performance[,2]);
qqline(trans_hosp_performance[,2]);

qqnorm(trans_hosp_performance[,3]);
qqline(trans_hosp_performance[,3]);


# Fig. 6.2 left panel
# Scatter plots on the original scale
plot(data, main="Scatter Plots of AMI, CHF, and CAP");

# Fig. 6.2 right panel;
# Scatter plots on the transformed scale 
plot(trans_data, main="Scatter Plots of Logit AMI, CHF, and CAP");


# Using complete-case analysis;
# Results in Table 6.1
cc_summary_matrix=matrix(NA, nrow=14, ncol=3);
cc_summary_matrix[1:3,1]=colMeans(hosp_performance, na.rm=T);
cc_summary_matrix[1:3,2]=cc_summary_matrix[1:3,1]-1.96*sqrt(diag(var(hosp_performance, na.rm=T))/obs_no);
cc_summary_matrix[1:3,3]=cc_summary_matrix[1:3,1]+1.96*sqrt(diag(var(hosp_performance, na.rm=T))/obs_no);

cc_summary_matrix[4,1]=mean(hosp_performance[,1]>0.9,na.rm=T);
cc_summary_matrix[4,2]=cc_summary_matrix[4,1]-1.96*sqrt(var(hosp_performance[,1]>0.9,na.rm=T)/obs_no[1]);
cc_summary_matrix[4,3]=cc_summary_matrix[4,1]+1.96*sqrt(var(hosp_performance[,1]>0.9,na.rm=T)/obs_no[1]);

cc_summary_matrix[5,1]=mean(hosp_performance[,1]>0.95,na.rm=T);
cc_summary_matrix[5,2]=cc_summary_matrix[5,1]-1.96*sqrt(var(hosp_performance[,1]>0.95,na.rm=T)/obs_no[1]);
cc_summary_matrix[5,3]=cc_summary_matrix[5,1]+1.96*sqrt(var(hosp_performance[,1]>0.95,na.rm=T)/obs_no[1]);

cc_summary_matrix[6,1]=mean(hosp_performance[,1]>0.99,na.rm=T);
cc_summary_matrix[6,2]=cc_summary_matrix[6,1]-1.96*sqrt(var(hosp_performance[,1]>0.99,na.rm=T)/obs_no[1]);
cc_summary_matrix[6,3]=cc_summary_matrix[6,1]+1.96*sqrt(var(hosp_performance[,1]>0.99,na.rm=T)/obs_no[1]);

cc_summary_matrix[7,1]=mean(hosp_performance[,2]>0.9,na.rm=T);
cc_summary_matrix[7,2]=cc_summary_matrix[7,1]-1.96*sqrt(var(hosp_performance[,2]>0.9,na.rm=T)/obs_no[2]);
cc_summary_matrix[7,3]=cc_summary_matrix[7,1]+1.96*sqrt(var(hosp_performance[,2]>0.9,na.rm=T)/obs_no[2]);

cc_summary_matrix[8,1]=mean(hosp_performance[,2]>0.95,na.rm=T);
cc_summary_matrix[8,2]=cc_summary_matrix[8,1]-1.96*sqrt(var(hosp_performance[,2]>0.95,na.rm=T)/obs_no[2]);
cc_summary_matrix[8,3]=cc_summary_matrix[8,1]+1.96*sqrt(var(hosp_performance[,2]>0.95,na.rm=T)/obs_no[2]);

cc_summary_matrix[9,1]=mean(hosp_performance[,2]>0.99,na.rm=T);
cc_summary_matrix[9,2]=cc_summary_matrix[9,1]-1.96*sqrt(var(hosp_performance[,2]>0.99,na.rm=T)/obs_no[2]);
cc_summary_matrix[9,3]=cc_summary_matrix[9,1]+1.96*sqrt(var(hosp_performance[,2]>0.99,na.rm=T)/obs_no[2]);

cc_summary_matrix[10,1]=mean(hosp_performance[,3]>0.9,na.rm=T);
cc_summary_matrix[10,2]=cc_summary_matrix[10,1]-1.96*sqrt(var(hosp_performance[,3]>0.9,na.rm=T)/obs_no[3]);
cc_summary_matrix[10,3]=cc_summary_matrix[10,1]+1.96*sqrt(var(hosp_performance[,3]>0.9,na.rm=T)/obs_no[3]);

cc_summary_matrix[11,1]=mean(hosp_performance[,3]>0.95,na.rm=T);
cc_summary_matrix[11,2]=cc_summary_matrix[11,1]-1.96*sqrt(var(hosp_performance[,3]>0.95,na.rm=T)/obs_no[3]);
cc_summary_matrix[11,3]=cc_summary_matrix[11,1]+1.96*sqrt(var(hosp_performance[,3]>0.95,na.rm=T)/obs_no[3]);

cc_summary_matrix[12,1]=mean(hosp_performance[,3]>0.99,na.rm=T);
cc_summary_matrix[12,2]=cc_summary_matrix[12,1]-1.96*sqrt(var(hosp_performance[,3]>0.99,na.rm=T)/obs_no[3]);
cc_summary_matrix[12,3]=cc_summary_matrix[12,1]+1.96*sqrt(var(hosp_performance[,3]>0.99,na.rm=T)/obs_no[3]);

cc_summary_matrix[13,1]=mean(hosp_performance[,1]>0.9 & hosp_performance[,2]> 0.9 & hosp_performance[,3]>0.9, na.rm=T);
cc_summary_matrix[13,2]=cc_summary_matrix[13,1]-1.96*sqrt(var(hosp_performance[,1]>0.9 & hosp_performance[,2]> 0.9 & hosp_performance[,3]>0.9, na.rm=T)/cc_no);
cc_summary_matrix[13,3]=cc_summary_matrix[13,1]+1.96*sqrt(var(hosp_performance[,1]>0.9 & hosp_performance[,2]> 0.9 & hosp_performance[,3]>0.9, na.rm=T)/cc_no);

cc_summary_matrix[14,1]=mean(hosp_performance[,1]>0.95 & hosp_performance[,2]> 0.95 & hosp_performance[,3]>0.95, na.rm=T);
cc_summary_matrix[14,2]=cc_summary_matrix[14,1]-1.96*sqrt(var(hosp_performance[,1]>0.95 & hosp_performance[,2]> 0.95 & hosp_performance[,3]>0.95, na.rm=T)/cc_no);
cc_summary_matrix[14,3]=cc_summary_matrix[14,1]+1.96*sqrt(var(hosp_performance[,1]>0.95 & hosp_performance[,2]> 0.95 & hosp_performance[,3]>0.95, na.rm=T)/cc_no);

cc_summary_matrix;

# use multivariate normal imputation;
# on the original scale;
# initialize using norm;
s=prelim.norm(hosp_performance);

# obtain the initial value for thetahat;
thetahat=em.norm(s)

mi_no=50;
norm_est_mat=norm_var_mat=matrix(NA, nrow=14, ncol=mi_no);

for(i in 1:mi_no)
{
rngseed(i);

theta=da.norm(s,thetahat,steps=200,showits=T);

# imputation;
hosp_norm_imp=imp.norm(s,theta,hosp_performance); 

# rounding the values;
hosp_norm_imp[hosp_norm_imp>=1]=1;
hosp_norm_imp[hosp_norm_imp<=0]=0;

# completed-data estimates;
norm_est_mat[1:3,i]=colMeans(hosp_norm_imp);
norm_var_mat[1:3,i]=diag(var(hosp_norm_imp))/rowobs;

norm_est_mat[4,i]=mean(hosp_norm_imp[,1]>0.9);
norm_var_mat[4,i]=var(hosp_norm_imp[,1]>0.9)/rowobs;
norm_est_mat[5,i]=mean(hosp_norm_imp[,1]>0.95);
norm_var_mat[5,i]=var(hosp_norm_imp[,1]>0.95)/rowobs;
norm_est_mat[6,i]=mean(hosp_norm_imp[,1]>0.99);
norm_var_mat[6,i]=var(hosp_norm_imp[,1]>0.99)/rowobs;

norm_est_mat[7,i]=mean(hosp_norm_imp[,2]>0.9);
norm_var_mat[7,i]=var(hosp_norm_imp[,2]>0.9)/rowobs;
norm_est_mat[8,i]=mean(hosp_norm_imp[,2]>0.95);
norm_var_mat[8,i]=var(hosp_norm_imp[,2]>0.95)/rowobs;
norm_est_mat[9,i]=mean(hosp_norm_imp[,2]>0.99);
norm_var_mat[9,i]=var(hosp_norm_imp[,2]>0.99)/rowobs;

norm_est_mat[10,i]=mean(hosp_norm_imp[,3]>0.9);
norm_var_mat[10,i]=var(hosp_norm_imp[,3]>0.9)/rowobs;
norm_est_mat[11,i]=mean(hosp_norm_imp[,3]>0.95);
norm_var_mat[11,i]=var(hosp_norm_imp[,3]>0.95)/rowobs;
norm_est_mat[12,i]=mean(hosp_norm_imp[,3]>0.99);
norm_var_mat[12,i]=var(hosp_norm_imp[,3]>0.99)/rowobs;

norm_est_mat[13,i]=mean(hosp_norm_imp[,1]>0.9 & hosp_norm_imp[,2]> 0.9 & hosp_norm_imp[,3]>0.9);
norm_var_mat[13,i]=var(hosp_norm_imp[,1]>0.9 & hosp_norm_imp[,2]> 0.9 & hosp_norm_imp[,3]>0.9)/rowobs;
norm_est_mat[14,i]=mean(hosp_norm_imp[,1]>0.95 & hosp_norm_imp[,2]> 0.95 & hosp_norm_imp[,3]>0.95);
norm_var_mat[14,i]=var(hosp_norm_imp[,1]>0.95 & hosp_norm_imp[,2]> 0.95 & hosp_norm_imp[,3]>0.95)/rowobs;

}

# get the posterior inference;

rngseed(197552);

theta0=thetahat;
getparam.norm(s, thetahat);

iter_no=1000;
parameter=matrix(NA, nrow=iter_no, ncol=9);

for (iter in 1:iter_no)
{
theta=da.norm(s,theta0, steps=1, showits=T);
parameter_iter=getparam.norm(s,theta, corr=TRUE);
parameter[iter,1:3]=parameter_iter$mu;
# parameter[iter,4:6]=parameter_iter$sigma[1,];
# parameter[iter,7:8]=parameter_iter$sigma[2,2:3];
# parameter[iter,9]=parameter_iter$sigma[3,3];
parameter[iter,4:6]=parameter_iter$sdv;
parameter[iter,7:8]=parameter_iter$r[1,2:3];
parameter[iter,9]=parameter_iter$r[2,3];
theta0=theta;
}

# Posterior means of model paraemters;
colMeans(parameter);

# summarize the inference;
# Results in Table 6.1
# NMI method;
norm_summary_matrix=matrix(NA, 14,3);
for (i in 1:14)
{
marginal_summary=pool.scalar(norm_est_mat[i,], norm_var_mat[i,], n=rowobs);
marginal_mi_mean=marginal_summary$qbar;
marginal_mi_var=marginal_summary$t;
marginal_mi_df=marginal_summary$df;
marginal_mi_f=marginal_summary$f;

norm_summary_matrix[i,1]=marginal_mi_mean;
norm_summary_matrix[i,2]=marginal_mi_mean-1.96*sqrt(marginal_mi_var);
norm_summary_matrix[i,3]=marginal_mi_mean+1.96*sqrt(marginal_mi_var);
}

norm_summary_matrix;


# multivariate normal imputation
# on transformed scale;
trans_s=prelim.norm(trans_hosp_performance);

# obtain the initial value for thetahat;
trans_thetahat=em.norm(trans_s)

mi_no=50;
trans_est_mat=trans_var_mat=matrix(NA, nrow=14, ncol=mi_no);

for(i in 1:mi_no)
{
rngseed(i);

trans_theta=da.norm(trans_s,trans_thetahat,steps=200,showits=T);

# imputation;
trans_hosp_imp_first=imp.norm(trans_s,trans_theta,trans_hosp_performance); 

# transforming back;
# trans_hosp_imp=(sin(trans_hosp_imp_first))^2;
trans_hosp_imp=exp(trans_hosp_imp_first)/(1+exp(trans_hosp_imp_first));
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
# Results in Table 6.1
# LMI method
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

# obtain the MCMC chains of the posterior samples of the parameters;
theta0=trans_thetahat;
getparam.norm(trans_s, trans_thetahat, corr=TRUE);

iter_no=2000;
parameter=matrix(NA, nrow=iter_no, ncol=9);

for (iter in 1:iter_no)
{
theta=da.norm(trans_s,theta0, steps=1, showits=T);
parameter_iter=getparam.norm(trans_s, theta, corr=TRUE);
parameter[iter,1:3]=parameter_iter$mu;
# parameter[iter,4:6]=parameter_iter$sigma[1,];
# parameter[iter,7:8]=parameter_iter$sigma[2,2:3];
# parameter[iter,9]=parameter_iter$sigma[3,3];
parameter[iter,4:6]=parameter_iter$sdv;
parameter[iter,7:8]=parameter_iter$r[1,2:3];
parameter[iter,9]=parameter_iter$r[2,3];
theta=theta0;
}

# Posterior means;
colMeans(parameter);

# make plots of the MCMC chain;
# Fig. 6.3

plot(seq(iter_no/2+1, iter_no, 1), parameter[(iter_no/2+1):iter_no,1], type="l", xlab="iteration", main="", ylab="mu_ami");
plot(seq(iter_no/2+1, iter_no, 1), parameter[(iter_no/2+1):iter_no,2], type="l", xlab="iteration", main="", ylab="mu_chf");
plot(seq(iter_no/2+1, iter_no, 1), parameter[(iter_no/2+1):iter_no,3], type="l", xlab="iteration", main="", ylab="mu_cap");
plot(seq(iter_no/2+1, iter_no, 1), parameter[(iter_no/2+1):iter_no,4], type="l", xlab="iteration", main="", ylab="sigma_ami");
plot(seq(iter_no/2+1, iter_no, 1), parameter[(iter_no/2+1):iter_no,5], type="l", xlab="iteration", main="", ylab="sigma_chf");
plot(seq(iter_no/2+1, iter_no, 1), parameter[(iter_no/2+1):iter_no,6], type="l", xlab="iteration", main="", ylab="sigma_cap");
plot(seq(iter_no/2+1, iter_no, 1), parameter[(iter_no/2+1):iter_no,7], type="l", xlab="iteration", main="", ylab="rho_ami_chf");
plot(seq(iter_no/2+1, iter_no, 1), parameter[(iter_no/2+1):iter_no,8], type="l", xlab="iteration", main="", ylab="rho_ami_cap");
plot(seq(iter_no/2+1, iter_no, 1), parameter[(iter_no/2+1):iter_no,9], type="l", xlab="iteration", main="", ylab="rho_chf_cap");
