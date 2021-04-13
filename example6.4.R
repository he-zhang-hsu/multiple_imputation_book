################## data generate ##############
# Program: boxcox_impute_expolore_simulation; ############
# generate data, and try out different missing data methods;
# run some simulations to assess the biasedness of various methods;

library(car);
library(mvtnorm);
library(mice);
library(norm);
library(mix);
# library(HI);
options(digits=4);
rm(list=ls());

data=read.table(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter6_multivariate_jointmodeling\\program\\sbirth.txt", col.names=c("RESTATUS","PLDEL3","REGNRES","CITRSPOP","METRORES",	"CNTRSPOP",	"DMAGE", "MRACE3",
				  "MEDUC6",		"DMAR",	"MPLBIRR",	"ADEQUACY",
				  "NLBNL",	"NLBND",	"NOTERM",	"DLIVORD",
				  "DTOTORD",	"MPRE5",	"DFAGE",	"FRACE4",	"NPREVIS",
				  "DGESTAT",	"CSEX",	"DBIRWT",	"DPLURAL",
				  "FMAPS", "DELMETH5",	"CIGAR",	"DRINK",
				  "WTGAIN", "CONGNTL",	"MEDICALRISK", "OBSTETRIC",
				  "NEWBORN", "LABCOMP",	"m",	"gestat_orig"), na.strings=".");

rowobs=nrow(data);


dgestat=data$DGESTAT
mean(dgestat, na.rm=TRUE);
sd(dgestat, na.rm=TRUE);
table(dgestat,useNA="always");
summary(dgestat);

dbirwt=data$DBIRWT;

mean(dbirwt, na.rm=TRUE);
sd(dbirwt, na.rm=TRUE);
summary(dbirwt);


frace=data$FRACE;

table(frace, useNA="always");

mrace=data$MRACE;

table(mrace, useNA="always");


meduc=data$MEDUC;
meduc_new=meduc;


# summary statistics;
summary(dgestat);
summary(dbirwt);
summary(frace);
table(frace, useNA="always");
table(mrace, useNA="always");
table(meduc, useNA="always");



# collapse the two less than high school into one group;
meduc_new[meduc_new==1 | meduc_new==2]=2;
meduc_new=meduc_new-1;

# dmar=data$DMAR;


# make some plots;
hist(dgestat, xlab="DGESTAT", main="");
qqnorm(dgestat);
qqline(dgestat);

hist(dbirwt, xlab="DBIRWT", main="");
qqnorm(dbirwt);
qqline(dbirwt);

lambda_com=powerTransform(dgestat~sqrt(dbirwt));

plot(dgestat, sqrt(dbirwt));

dbirwt_sqrt=sqrt(dbirwt);

hist(dbirwt_sqrt, xlab="Square root of DBIRWT", main="");

plot(dbirwt_sqrt, dgestat, xlab="Square root of DBIRWT", ylab="DGESTAT");


# form the new dataset;
ori_data=cbind(frace, mrace, meduc_new, dgestat, dbirwt_sqrt);

# ori_data=cbind(frace, mrace, dmar, dgestat, dbirwt);



# apply 10% missing data to mrace, dmar, and dbirwt;
set.seed(197789);

miss_indi=matrix(rbinom(n=3*rowobs, size=1, prob=0.1), nrow=rowobs, ncol=3);

# set some of the ori_data as missing values;
ori_data[,c(2,3,5)][miss_indi==1]=NA;

# the CC missingness rates;
# the CC missingness rates;


meduc_new=ori_data[,3];
mrace_new=ori_data[,2];
frace_new=ori_data[,1];
dgestat_new=ori_data[,4];
dbirwt_new=ori_data[,5]^2;


mrace_dummy1=mrace_dummy2=rep(NA,rowobs);
mrace_dummy1[mrace_new==1]=0;
mrace_dummy2[mrace_new==1]=0;
mrace_dummy1[mrace_new==2]=1;
mrace_dummy2[mrace_new==2]=0;
mrace_dummy1[mrace_new==3]=0;
mrace_dummy2[mrace_new==3]=1;

frace_dummy1=frace_dummy2=rep(NA,rowobs);
frace_dummy1[frace_new==1]=0;
frace_dummy2[frace_new==1]=0;
frace_dummy1[frace_new==2]=1;
frace_dummy2[frace_new==2]=0;
frace_dummy1[frace_new==3]=0;
frace_dummy2[frace_new==3]=1;

meduc_dummy1=meduc_dummy2=meduc_dummy3=rep(NA,rowobs);
meduc_dummy1[meduc_new==1]=0;
meduc_dummy2[meduc_new==1]=0;
meduc_dummy3[meduc_new==1]=0;


meduc_dummy1[meduc_new==2]=1;
meduc_dummy2[meduc_new==2]=0;
meduc_dummy3[meduc_new==2]=0;


meduc_dummy1[meduc_new==3]=0;
meduc_dummy2[meduc_new==3]=1;
meduc_dummy3[meduc_new==3]=0;


meduc_dummy1[meduc_new==4]=0;
meduc_dummy2[meduc_new==4]=0;
meduc_dummy3[meduc_new==4]=1;




# run a regression analysis;
cc_main=lm(dgestat_new ~ dbirwt_new+I(dbirwt_new*dbirwt_new)+mrace_dummy1+mrace_dummy2+frace_dummy1+frace_dummy2+meduc_dummy1+meduc_dummy2+meduc_dummy3);
summary(cc_main);
cc_summary_matrix=matrix(NA, 10, 3);

cc_summary_matrix[,1]=summary(cc_main)$coef[,1];
cc_summary_matrix[,2]=cc_summary_matrix[,1]-1.96*summary(cc_main)$coef[,2];
cc_summary_matrix[,3]=cc_summary_matrix[,1]+1.96*summary(cc_main)$coef[,2];


# general location model imputation;
s=prelim.mix(ori_data,3);
thetahat=em.mix(s);
getparam.mix(s, thetahat, corr=TRUE);

mi_no=50;
sat_est_mat=sat_var_mat=matrix(NA, nrow=10, ncol=mi_no);


for (i in 1:mi_no)
{
rngseed(i); # set random number generator seed
newtheta=da.mix(s, thetahat, steps=200, showits=TRUE) # take 200 steps
# getparam.mix(s, newtheta, corr=TRUE);
imp_data=imp.mix(s, newtheta);

# set up the dummies;
meduc_new=imp_data[,3];
mrace_new=imp_data[,2];
frace_new=imp_data[,1];
dgestat_new=imp_data[,4];
dbirwt_new=imp_data[,5]^2;


mrace_dummy1=mrace_dummy2=rep(NA,rowobs);
mrace_dummy1[mrace_new==1]=0;
mrace_dummy2[mrace_new==1]=0;
mrace_dummy1[mrace_new==2]=1;
mrace_dummy2[mrace_new==2]=0;
mrace_dummy1[mrace_new==3]=0;
mrace_dummy2[mrace_new==3]=1;

frace_dummy1=frace_dummy2=rep(NA,rowobs);
frace_dummy1[frace_new==1]=0;
frace_dummy2[frace_new==1]=0;
frace_dummy1[frace_new==2]=1;
frace_dummy2[frace_new==2]=0;
frace_dummy1[frace_new==3]=0;
frace_dummy2[frace_new==3]=1;

meduc_dummy1=meduc_dummy2=meduc_dummy3=rep(NA,rowobs);
meduc_dummy1[meduc_new==1]=0;
meduc_dummy2[meduc_new==1]=0;
meduc_dummy3[meduc_new==1]=0;


meduc_dummy1[meduc_new==2]=1;
meduc_dummy2[meduc_new==2]=0;
meduc_dummy3[meduc_new==2]=0;


meduc_dummy1[meduc_new==3]=0;
meduc_dummy2[meduc_new==3]=1;
meduc_dummy3[meduc_new==3]=0;


meduc_dummy1[meduc_new==4]=0;
meduc_dummy2[meduc_new==4]=0;
meduc_dummy3[meduc_new==4]=1;

# run a regression analysis;
# run a regression analysis;
impute_main=lm(dgestat_new ~ dbirwt_new+I(dbirwt_new*dbirwt_new)+mrace_dummy1+mrace_dummy2+frace_dummy1+frace_dummy2+meduc_dummy1+meduc_dummy2+meduc_dummy3);
sat_est_mat[,i]=summary(impute_main)$coef[,1];
sat_var_mat[,i]=summary(impute_main)$coef[,2]^2;

}

# summarize the inference;
sat_summary_matrix=matrix(NA,10,4);
for (i in 1:10)
{
sat_summary=pool.scalar(sat_est_mat[i,], sat_var_mat[i,], n=rowobs);
sat_mi_mean=sat_summary$qbar;
sat_mi_var=sat_summary$t;
sat_mi_df=sat_summary$df;
sat_mi_f=sat_summary$f;

sat_summary_matrix[i,1]=sat_mi_mean;
sat_summary_matrix[i,2]=sat_mi_mean-qt(.975, sat_mi_df)*sqrt(sat_mi_var);
sat_summary_matrix[i,3]=sat_mi_mean+qt(.975, sat_mi_df)*sqrt(sat_mi_var);
sat_summary_matrix[i,4]=sqrt(sat_mi_var);
}

sat_summary_matrix;


# imputation just using the linear term of dbirwt;

ori_data=cbind(frace, mrace, meduc_new, dgestat, dbirwt);


# set some of the ori_data as missing values;
ori_data[,c(2,3,5)][miss_indi==1]=NA;


# general location model imputation;
s=prelim.mix(ori_data,3);
thetahat=em.mix(s);
getparam.mix(s, thetahat, corr=TRUE);

mi_no=50;
sat_est_mat=sat_var_mat=matrix(NA, nrow=10, ncol=mi_no);


for (i in 1:mi_no)
{
rngseed(i); # set random number generator seed
newtheta=da.mix(s, thetahat, steps=200, showits=TRUE) # take 200 steps
# getparam.mix(s, newtheta, corr=TRUE);
imp_data=imp.mix(s, newtheta);

# set up the dummies;
meduc_new=imp_data[,3];
mrace_new=imp_data[,2];
frace_new=imp_data[,1];
dgestat_new=imp_data[,4];
# linear term;
dbirwt_new=imp_data[,5];


mrace_dummy1=mrace_dummy2=rep(NA,rowobs);
mrace_dummy1[mrace_new==1]=0;
mrace_dummy2[mrace_new==1]=0;
mrace_dummy1[mrace_new==2]=1;
mrace_dummy2[mrace_new==2]=0;
mrace_dummy1[mrace_new==3]=0;
mrace_dummy2[mrace_new==3]=1;

frace_dummy1=frace_dummy2=rep(NA,rowobs);
frace_dummy1[frace_new==1]=0;
frace_dummy2[frace_new==1]=0;
frace_dummy1[frace_new==2]=1;
frace_dummy2[frace_new==2]=0;
frace_dummy1[frace_new==3]=0;
frace_dummy2[frace_new==3]=1;

meduc_dummy1=meduc_dummy2=meduc_dummy3=rep(NA,rowobs);
meduc_dummy1[meduc_new==1]=0;
meduc_dummy2[meduc_new==1]=0;
meduc_dummy3[meduc_new==1]=0;


meduc_dummy1[meduc_new==2]=1;
meduc_dummy2[meduc_new==2]=0;
meduc_dummy3[meduc_new==2]=0;


meduc_dummy1[meduc_new==3]=0;
meduc_dummy2[meduc_new==3]=1;
meduc_dummy3[meduc_new==3]=0;


meduc_dummy1[meduc_new==4]=0;
meduc_dummy2[meduc_new==4]=0;
meduc_dummy3[meduc_new==4]=1;

# run a regression analysis;
# run a regression analysis;
impute_main=lm(dgestat_new ~ dbirwt_new+I(dbirwt_new*dbirwt_new)+mrace_dummy1+mrace_dummy2+frace_dummy1+frace_dummy2+meduc_dummy1+meduc_dummy2+meduc_dummy3);
sat_est_mat[,i]=summary(impute_main)$coef[,1];
sat_var_mat[,i]=summary(impute_main)$coef[,2]^2;

}

# summarize the inference;
sat_summary_matrix=matrix(NA,10,4);
for (i in 1:10)
{
sat_summary=pool.scalar(sat_est_mat[i,], sat_var_mat[i,], n=rowobs);
sat_mi_mean=sat_summary$qbar;
sat_mi_var=sat_summary$t;
sat_mi_df=sat_summary$df;
sat_mi_f=sat_summary$f;

sat_summary_matrix[i,1]=sat_mi_mean;
sat_summary_matrix[i,2]=sat_mi_mean-qt(.975, sat_mi_df)*sqrt(sat_mi_var);
sat_summary_matrix[i,3]=sat_mi_mean+qt(.975, sat_mi_df)*sqrt(sat_mi_var);
sat_summary_matrix[i,4]=sqrt(sat_mi_var);
}

sat_summary_matrix;

# latent-variable model imputation using R jomo;

# form the new dataset;
ori_data=cbind(frace, mrace, meduc_new, dgestat, dbirwt_sqrt);



# set some of the ori_data as missing values;
ori_data[,c(2,3,5)][miss_indi==1]=NA;



ori_data=as.data.frame(ori_data);

ori_data$frace=as.factor(ori_data$frace);
ori_data$mrace=as.factor(ori_data$mrace);
ori_data$meduc_new=as.factor(ori_data$meduc_new);


Y.con=ori_data[,c("dgestat","dbirwt_sqrt")];
 
Y.cat=ori_data[,c("frace", "mrace", "meduc_new"), drop=FALSE];
 
Y.numcat=c(3,3,4);

# X=data.frame(rep(1,300),sldata[,c("sex")]) 
# colnames(X)<-c("const", "sex") 
# beta.start<-matrix(0,2,5) 
# l1cov.start<-diag(1,5) 
# l1cov.prior=diag(1,5); 
nburn=as.integer(2000); 
nbetween=as.integer(200); 
nimp=mi_no;
im
imp<-jomo1mix(Y.con=Y.con,Y.cat=Y.cat,Y.numcat=Y.numcat,X=NULL,beta.start=NULL,l1cov.start=NULL,
l1cov.prior=NULL,nburn=nburn,nbetween=nbetween,nimp=nimp, output=1, out.iter=10);


prob_est_mat=prob_var_mat=matrix(NA, nrow=10, ncol=mi_no);
imp_data=matrix(NA, rowobs,5);

for (i in 1:mi_no)
{


ori_data_imp=imp[imp$Imputation==i,1:5];
imp_data[,1]=as.numeric(ori_data_imp[,3]);
imp_data[,2]=as.numeric(ori_data_imp[,4]);
imp_data[,3]=as.numeric(ori_data_imp[,5]);
imp_data[,4]=ori_data_imp[,1];
imp_data[,5]=ori_data_imp[,2];





# set up the dummies;
meduc_new=imp_data[,3];
mrace_new=imp_data[,2];
frace_new=imp_data[,1];
dgestat_new=imp_data[,4];
# quadratic term;
dbirwt_new=imp_data[,5]^2;


mrace_dummy1=mrace_dummy2=rep(NA,rowobs);
mrace_dummy1[mrace_new==1]=0;
mrace_dummy2[mrace_new==1]=0;
mrace_dummy1[mrace_new==2]=1;
mrace_dummy2[mrace_new==2]=0;
mrace_dummy1[mrace_new==3]=0;
mrace_dummy2[mrace_new==3]=1;

frace_dummy1=frace_dummy2=rep(NA,rowobs);
frace_dummy1[frace_new==1]=0;
frace_dummy2[frace_new==1]=0;
frace_dummy1[frace_new==2]=1;
frace_dummy2[frace_new==2]=0;
frace_dummy1[frace_new==3]=0;
frace_dummy2[frace_new==3]=1;

meduc_dummy1=meduc_dummy2=meduc_dummy3=rep(NA,rowobs);
meduc_dummy1[meduc_new==1]=0;
meduc_dummy2[meduc_new==1]=0;
meduc_dummy3[meduc_new==1]=0;


meduc_dummy1[meduc_new==2]=1;
meduc_dummy2[meduc_new==2]=0;
meduc_dummy3[meduc_new==2]=0;


meduc_dummy1[meduc_new==3]=0;
meduc_dummy2[meduc_new==3]=1;
meduc_dummy3[meduc_new==3]=0;


meduc_dummy1[meduc_new==4]=0;
meduc_dummy2[meduc_new==4]=0;
meduc_dummy3[meduc_new==4]=1;

# run a regression analysis;
# run a regression analysis;
impute_main=lm(dgestat_new ~ dbirwt_new+I(dbirwt_new*dbirwt_new)+mrace_dummy1+mrace_dummy2+frace_dummy1+frace_dummy2+meduc_dummy1+meduc_dummy2+meduc_dummy3);
prob_est_mat[,i]=summary(impute_main)$coef[,1];
prob_var_mat[,i]=summary(impute_main)$coef[,2]^2;

}

# summarize the inference;
prob_summary_matrix=matrix(NA,10,4);
for (i in 1:10)
{
prob_summary=pool.scalar(prob_est_mat[i,], prob_var_mat[i,], n=rowobs);
prob_mi_mean=prob_summary$qbar;
prob_mi_var=prob_summary$t;
prob_mi_df=prob_summary$df;
prob_mi_f=prob_summary$f;

prob_summary_matrix[i,1]=prob_mi_mean;
prob_summary_matrix[i,2]=prob_mi_mean-qt(.975, prob_mi_df)*sqrt(prob_mi_var);
prob_summary_matrix[i,3]=prob_mi_mean+qt(.975, prob_mi_df)*sqrt(prob_mi_var);
prob_summary_matrix[i,4]=sqrt(prob_mi_var);
}

prob_summary_matrix;







# posterior estiamtes of the bivariate probit models;
imp<-jomo1mix.MCMCchain(Y.con=Y.con,Y.cat=Y.cat,Y.numcat=Y.numcat,nburn=5000)

plot(c(1:5000),imp$collectbeta[1,1,1:5000],type="l")
#Or similarly we can check the convergence of any element of omega:
plot(c(1:5000),imp$collectomega[2,2,1:5000],type="l")







