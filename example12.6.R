library(openxlsx)
library(mice);
library(norm);
# adaptive metroplis rejection sampling;
library(HI);
# include the library MASS for generating from multivariate normal distribution;
library(MASS);
# include the library msm and bayesm for using the function of drawing truncated normal
# distribution;
library(msm);
library(bayesm);
library(mvtnorm);
library(sn);
library(arm);
# library(MCMCpack);
library(cat);
library(jomo);


# read the data set from Zhou et al. 2017;
cdat<-read.xlsx("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter4_univariate_model_parametric\\program\\paper1_dataset.xlsx", sheet=1, startRow=1, colNames=T);



# fetch the outcome variable;
table(cdat$carercvd);

carercvd=as.numeric(cdat$carercvd);
carercvd_obs=carercvd[!is.na(carercvd)];

rowobs=length(carercvd);
mis_no=sum(is.na(carercvd));
obs_no=rowobs-mis_no;


# create a binary outcome variable of carercvd;
carercvd_23=carercvd;
# collapse the category 2 and 3;
carercvd_23[carercvd_23>=2]=2;

# recode carercvd_23 into a binary indicator;
carercvd_23=carercvd_23-1;

# 1 very satisfied and 0 somewhat satisfied or not at all satisfied;

carercvd_23=1-carercvd_23;

# fetch the covariates;
# gender;
sex=(cdat$sex);

# general health;
genhlth=cdat$genhlth;

# education level;
X_educag=cdat$X_educag;

# having health coverage;

hlthpln1=cdat$hlthpln1;

# having delayed health care;
delaym=cdat$delaym;

# coverage=1 yes, 0 no;
coverage=2-hlthpln1;


selfhealth=genhlth;
selfhealth[selfhealth<=3]=1;
selfhealth[selfhealth>3]=0;

table(selfhealth);
table(coverage);
table(delaym);


# creat the orginal dataset;
# carercvd (1: very satisfied; 0: somewhat or not satisfied;
# selfhealth (1: excellent/very good/good; 0: fair/poor);
# coverage (1: yes; 0: no);
# delaym (1:yes; 0: no);

ori_data=cbind(carercvd_23, selfhealth, coverage, delaym);

# creat 10% missing data for selfhealth, coverage, and delaym;
set.seed(197789);

miss_indi=matrix(rbinom(n=3*rowobs, size=1, prob=0.1), nrow=rowobs, ncol=3);

# set some of the ori_data as missing values;
ori_data[,2:4][miss_indi==1]=NA;

# complete-case analysis;
cc_est_mat=cc_var_mat=rep(NA,9);
ori_data_1=ori_data[!is.na(ori_data[,1]),1];
ori_data_cc=ori_data[complete.cases(ori_data),];
ori_data_cc_111=ori_data_cc[ori_data_cc[,2]==1 & ori_data_cc[,3]==1 & ori_data_cc[,4]==1,];
ori_data_cc_110=ori_data_cc[ori_data_cc[,2]==1 & ori_data_cc[,3]==1 & ori_data_cc[,4]==0,];
ori_data_cc_101=ori_data_cc[ori_data_cc[,2]==1 & ori_data_cc[,3]==0 & ori_data_cc[,4]==1,];
ori_data_cc_100=ori_data_cc[ori_data_cc[,2]==1 & ori_data_cc[,3]==0 & ori_data_cc[,4]==0,];
ori_data_cc_011=ori_data_cc[ori_data_cc[,2]==0 & ori_data_cc[,3]==1 & ori_data_cc[,4]==1,];
ori_data_cc_010=ori_data_cc[ori_data_cc[,2]==0 & ori_data_cc[,3]==1 & ori_data_cc[,4]==0,];
ori_data_cc_001=ori_data_cc[ori_data_cc[,2]==0 & ori_data_cc[,3]==0 & ori_data_cc[,4]==1,];
ori_data_cc_000=ori_data_cc[ori_data_cc[,2]==0 & ori_data_cc[,3]==0 & ori_data_cc[,4]==0,];

nrow(ori_data_cc)/rowobs;




cc_est_mat[1]=mean(ori_data_1);
cc_var_mat[1]=var(ori_data_1)/length(ori_data_1);

cc_est_mat[2]=mean(ori_data_cc_111[,1]);
cc_var_mat[2]=var(ori_data_cc_111[,1])/nrow(ori_data_cc_111);

cc_est_mat[3]=mean(ori_data_cc_110[,1]);
cc_var_mat[3]=var(ori_data_cc_110[,1])/nrow(ori_data_cc_110);

cc_est_mat[4]=mean(ori_data_cc_101[,1]);
cc_var_mat[4]=var(ori_data_cc_101[,1])/nrow(ori_data_cc_101);

cc_est_mat[5]=mean(ori_data_cc_100[,1]);
cc_var_mat[5]=var(ori_data_cc_100[,1])/nrow(ori_data_cc_100);

cc_est_mat[6]=mean(ori_data_cc_011[,1]);
cc_var_mat[6]=var(ori_data_cc_011[,1])/nrow(ori_data_cc_011);

cc_est_mat[7]=mean(ori_data_cc_010[,1]);
cc_var_mat[7]=var(ori_data_cc_010[,1])/nrow(ori_data_cc_010);

cc_est_mat[8]=mean(ori_data_cc_001[,1]);
cc_var_mat[8]=var(ori_data_cc_001[,1])/nrow(ori_data_cc_001);

cc_est_mat[9]=mean(ori_data_cc_000[,1]);
cc_var_mat[9]=var(ori_data_cc_000[,1])/nrow(ori_data_cc_000);

cc_summary_matrix=matrix(NA,9,3);

cc_summary_matrix[,1]=cc_est_mat;
cc_summary_matrix[,2]=cc_est_mat-1.96*sqrt(cc_var_mat);
cc_summary_matrix[,3]=cc_est_mat+1.96*sqrt(cc_var_mat);

cc_summary_matrix;


# prepare for the imputation;
# library cat takes 1/2 factors;

ori_data_new=ori_data+1;

ori_data_concat=rbind(ori_data_new, matrix(NA, rowobs, 4));

# using cat to do the imputation;
setup=prelim.cat(ori_data_concat);

# saturated model;
# test_margins=c(1,2,3);


theta_saturated=em.cat(setup);

mi_no=1000;
sat_est_mat=sat_rep_mat=matrix(NA, nrow=9, ncol=mi_no);


for (i in 1:mi_no)
{

rngseed(i);

theta_s=da.cat(setup,theta_saturated,steps=200,showits=T);

# imputation;
ori_data_new_imp=imp.cat(setup,theta_s); 
ori_data_imp=ori_data_new_imp-1;

ori_data_mi=ori_data_imp[1:rowobs,];
ori_data_rep=ori_data_imp[(rowobs+1):(2*rowobs),];

# for imputation;

# marginal mean of carercvd;
sat_est_mat[1,i]=mean(ori_data_mi[,1]);

# means of carercvd across 8 groups;
ori_data_mi_111=ori_data_mi[ori_data_mi[,2]==1 & ori_data_mi[,3]==1 & ori_data_mi[,4]==1,];
ori_data_mi_110=ori_data_mi[ori_data_mi[,2]==1 & ori_data_mi[,3]==1 & ori_data_mi[,4]==0,];
ori_data_mi_101=ori_data_mi[ori_data_mi[,2]==1 & ori_data_mi[,3]==0 & ori_data_mi[,4]==1,];
ori_data_mi_100=ori_data_mi[ori_data_mi[,2]==1 & ori_data_mi[,3]==0 & ori_data_mi[,4]==0,];
ori_data_mi_011=ori_data_mi[ori_data_mi[,2]==0 & ori_data_mi[,3]==1 & ori_data_mi[,4]==1,];
ori_data_mi_010=ori_data_mi[ori_data_mi[,2]==0 & ori_data_mi[,3]==1 & ori_data_mi[,4]==0,];
ori_data_mi_001=ori_data_mi[ori_data_mi[,2]==0 & ori_data_mi[,3]==0 & ori_data_mi[,4]==1,];
ori_data_mi_000=ori_data_mi[ori_data_mi[,2]==0 & ori_data_mi[,3]==0 & ori_data_mi[,4]==0,];

sat_est_mat[2,i]=mean(ori_data_mi_111[,1]);
sat_est_mat[3,i]=mean(ori_data_mi_110[,1]);
sat_est_mat[4,i]=mean(ori_data_mi_101[,1]);
sat_est_mat[5,i]=mean(ori_data_mi_100[,1]);
sat_est_mat[6,i]=mean(ori_data_mi_011[,1]);
sat_est_mat[7,i]=mean(ori_data_mi_010[,1]);
sat_est_mat[8,i]=mean(ori_data_mi_001[,1]);
sat_est_mat[9,i]=mean(ori_data_mi_000[,1]);

# for replicates;

# marginal mean of carercvd;
sat_rep_mat[1,i]=mean(ori_data_rep[,1]);

# means of carercvd across 8 groups;
ori_data_rep_111=ori_data_rep[ori_data_rep[,2]==1 & ori_data_rep[,3]==1 & ori_data_rep[,4]==1,];
ori_data_rep_110=ori_data_rep[ori_data_rep[,2]==1 & ori_data_rep[,3]==1 & ori_data_rep[,4]==0,];
ori_data_rep_101=ori_data_rep[ori_data_rep[,2]==1 & ori_data_rep[,3]==0 & ori_data_rep[,4]==1,];
ori_data_rep_100=ori_data_rep[ori_data_rep[,2]==1 & ori_data_rep[,3]==0 & ori_data_rep[,4]==0,];
ori_data_rep_011=ori_data_rep[ori_data_rep[,2]==0 & ori_data_rep[,3]==1 & ori_data_rep[,4]==1,];
ori_data_rep_010=ori_data_rep[ori_data_rep[,2]==0 & ori_data_rep[,3]==1 & ori_data_rep[,4]==0,];
ori_data_rep_001=ori_data_rep[ori_data_rep[,2]==0 & ori_data_rep[,3]==0 & ori_data_rep[,4]==1,];
ori_data_rep_000=ori_data_rep[ori_data_rep[,2]==0 & ori_data_rep[,3]==0 & ori_data_rep[,4]==0,];

sat_rep_mat[2,i]=mean(ori_data_rep_111[,1]);
sat_rep_mat[3,i]=mean(ori_data_rep_110[,1]);
sat_rep_mat[4,i]=mean(ori_data_rep_101[,1]);
sat_rep_mat[5,i]=mean(ori_data_rep_100[,1]);
sat_rep_mat[6,i]=mean(ori_data_rep_011[,1]);
sat_rep_mat[7,i]=mean(ori_data_rep_010[,1]);
sat_rep_mat[8,i]=mean(ori_data_rep_001[,1]);
sat_rep_mat[9,i]=mean(ori_data_rep_000[,1]);



}


# summarize the differences;
sat_diff_mat=sat_rep_mat-sat_est_mat;

# posterior predictive P-values;
rowMeans(sat_diff_mat>0);

rowMeans(sat_rep_mat);
rowMeans(sat_est_mat);



# multiple imputation using main effects only;
main_margins=c(1,0,2,0,3,0,4);
theta_margin=ecm.cat(setup, margins=main_margins);



main_est_mat=main_rep_mat=matrix(NA, nrow=9, ncol=mi_no);


for (i in 1:mi_no)
{

rngseed(i);

theta_m=dabipf(setup,main_margins, theta_margin,steps=200,showits=T);

# imputation;
ori_data_new_imp=imp.cat(setup,theta_m); 
ori_data_imp=ori_data_new_imp-1;

ori_data_mi=ori_data_imp[1:rowobs,];
ori_data_rep=ori_data_imp[(rowobs+1):(2*rowobs),];

# for imputation;

# marginal mean of carercvd;
main_est_mat[1,i]=mean(ori_data_mi[,1]);

# means of carercvd across 8 groups;
ori_data_mi_111=ori_data_mi[ori_data_mi[,2]==1 & ori_data_mi[,3]==1 & ori_data_mi[,4]==1,];
ori_data_mi_110=ori_data_mi[ori_data_mi[,2]==1 & ori_data_mi[,3]==1 & ori_data_mi[,4]==0,];
ori_data_mi_101=ori_data_mi[ori_data_mi[,2]==1 & ori_data_mi[,3]==0 & ori_data_mi[,4]==1,];
ori_data_mi_100=ori_data_mi[ori_data_mi[,2]==1 & ori_data_mi[,3]==0 & ori_data_mi[,4]==0,];
ori_data_mi_011=ori_data_mi[ori_data_mi[,2]==0 & ori_data_mi[,3]==1 & ori_data_mi[,4]==1,];
ori_data_mi_010=ori_data_mi[ori_data_mi[,2]==0 & ori_data_mi[,3]==1 & ori_data_mi[,4]==0,];
ori_data_mi_001=ori_data_mi[ori_data_mi[,2]==0 & ori_data_mi[,3]==0 & ori_data_mi[,4]==1,];
ori_data_mi_000=ori_data_mi[ori_data_mi[,2]==0 & ori_data_mi[,3]==0 & ori_data_mi[,4]==0,];

main_est_mat[2,i]=mean(ori_data_mi_111[,1]);
main_est_mat[3,i]=mean(ori_data_mi_110[,1]);
main_est_mat[4,i]=mean(ori_data_mi_101[,1]);
main_est_mat[5,i]=mean(ori_data_mi_100[,1]);
main_est_mat[6,i]=mean(ori_data_mi_011[,1]);
main_est_mat[7,i]=mean(ori_data_mi_010[,1]);
main_est_mat[8,i]=mean(ori_data_mi_001[,1]);
main_est_mat[9,i]=mean(ori_data_mi_000[,1]);

# for replicates;

# marginal mean of carercvd;
main_rep_mat[1,i]=mean(ori_data_rep[,1]);

# means of carercvd across 8 groups;
ori_data_rep_111=ori_data_rep[ori_data_rep[,2]==1 & ori_data_rep[,3]==1 & ori_data_rep[,4]==1,];
ori_data_rep_110=ori_data_rep[ori_data_rep[,2]==1 & ori_data_rep[,3]==1 & ori_data_rep[,4]==0,];
ori_data_rep_101=ori_data_rep[ori_data_rep[,2]==1 & ori_data_rep[,3]==0 & ori_data_rep[,4]==1,];
ori_data_rep_100=ori_data_rep[ori_data_rep[,2]==1 & ori_data_rep[,3]==0 & ori_data_rep[,4]==0,];
ori_data_rep_011=ori_data_rep[ori_data_rep[,2]==0 & ori_data_rep[,3]==1 & ori_data_rep[,4]==1,];
ori_data_rep_010=ori_data_rep[ori_data_rep[,2]==0 & ori_data_rep[,3]==1 & ori_data_rep[,4]==0,];
ori_data_rep_001=ori_data_rep[ori_data_rep[,2]==0 & ori_data_rep[,3]==0 & ori_data_rep[,4]==1,];
ori_data_rep_000=ori_data_rep[ori_data_rep[,2]==0 & ori_data_rep[,3]==0 & ori_data_rep[,4]==0,];

main_rep_mat[2,i]=mean(ori_data_rep_111[,1]);
main_rep_mat[3,i]=mean(ori_data_rep_110[,1]);
main_rep_mat[4,i]=mean(ori_data_rep_101[,1]);
main_rep_mat[5,i]=mean(ori_data_rep_100[,1]);
main_rep_mat[6,i]=mean(ori_data_rep_011[,1]);
main_rep_mat[7,i]=mean(ori_data_rep_010[,1]);
main_rep_mat[8,i]=mean(ori_data_rep_001[,1]);
main_rep_mat[9,i]=mean(ori_data_rep_000[,1]);





}

# summarize the differences;
main_diff_mat=main_rep_mat-main_est_mat;

# posterior predictive P-values;
rowMeans(main_diff_mat>0);

rowMeans(main_rep_mat);
rowMeans(main_est_mat);



# multiple imputation using bivariate probit models;


prob_est_mat=prob_rep_mat=matrix(NA, nrow=9, ncol=mi_no);

# data processing;

ori_data_concat=rbind(ori_data, matrix(NA, rowobs, 4));

ori_data_factor=as.data.frame(ori_data_concat);

ori_data_factor$carercvd_23<-as.factor(ori_data_factor$carercvd_23);
ori_data_factor$selfhealth<-as.factor(ori_data_factor$selfhealth);
ori_data_factor$coverage<-as.factor(ori_data_factor$coverage);
ori_data_factor$delaym<-as.factor(ori_data_factor$delaym);

nburn=1000;

imp=jomo1cat(ori_data_factor, rep(2,4), X=NULL, beta.start=NULL, l1cov.start=NULL, l1cov.prior=NULL, nburn=nburn, nbetween=200, nimp=mi_no,output=1, out.iter=10)


for (i in 1:mi_no)
{

# fetch the imputation;

ori_data_imp=imp[imp$Imputation==i,1:4];
ori_data_imp[,1]=as.numeric(ori_data_imp[,1]);
ori_data_imp[,2]=as.numeric(ori_data_imp[,2]);
ori_data_imp[,3]=as.numeric(ori_data_imp[,3]);
ori_data_imp[,4]=as.numeric(ori_data_imp[,4]);

ori_data_imp=ori_data_imp-1;

ori_data_mi=ori_data_imp[1:rowobs,];
ori_data_rep=ori_data_imp[(rowobs+1):(2*rowobs),];

# for imputation;

# marginal mean of carercvd;
prob_est_mat[1,i]=mean(ori_data_mi[,1]);

# means of carercvd across 8 groups;
ori_data_mi_111=ori_data_mi[ori_data_mi[,2]==1 & ori_data_mi[,3]==1 & ori_data_mi[,4]==1,];
ori_data_mi_110=ori_data_mi[ori_data_mi[,2]==1 & ori_data_mi[,3]==1 & ori_data_mi[,4]==0,];
ori_data_mi_101=ori_data_mi[ori_data_mi[,2]==1 & ori_data_mi[,3]==0 & ori_data_mi[,4]==1,];
ori_data_mi_100=ori_data_mi[ori_data_mi[,2]==1 & ori_data_mi[,3]==0 & ori_data_mi[,4]==0,];
ori_data_mi_011=ori_data_mi[ori_data_mi[,2]==0 & ori_data_mi[,3]==1 & ori_data_mi[,4]==1,];
ori_data_mi_010=ori_data_mi[ori_data_mi[,2]==0 & ori_data_mi[,3]==1 & ori_data_mi[,4]==0,];
ori_data_mi_001=ori_data_mi[ori_data_mi[,2]==0 & ori_data_mi[,3]==0 & ori_data_mi[,4]==1,];
ori_data_mi_000=ori_data_mi[ori_data_mi[,2]==0 & ori_data_mi[,3]==0 & ori_data_mi[,4]==0,];

prob_est_mat[2,i]=mean(ori_data_mi_111[,1]);
prob_est_mat[3,i]=mean(ori_data_mi_110[,1]);
prob_est_mat[4,i]=mean(ori_data_mi_101[,1]);
prob_est_mat[5,i]=mean(ori_data_mi_100[,1]);
prob_est_mat[6,i]=mean(ori_data_mi_011[,1]);
prob_est_mat[7,i]=mean(ori_data_mi_010[,1]);
prob_est_mat[8,i]=mean(ori_data_mi_001[,1]);
prob_est_mat[9,i]=mean(ori_data_mi_000[,1]);

# for replicates;

# marginal mean of carercvd;
prob_rep_mat[1,i]=mean(ori_data_rep[,1]);

# means of carercvd across 8 groups;
ori_data_rep_111=ori_data_rep[ori_data_rep[,2]==1 & ori_data_rep[,3]==1 & ori_data_rep[,4]==1,];
ori_data_rep_110=ori_data_rep[ori_data_rep[,2]==1 & ori_data_rep[,3]==1 & ori_data_rep[,4]==0,];
ori_data_rep_101=ori_data_rep[ori_data_rep[,2]==1 & ori_data_rep[,3]==0 & ori_data_rep[,4]==1,];
ori_data_rep_100=ori_data_rep[ori_data_rep[,2]==1 & ori_data_rep[,3]==0 & ori_data_rep[,4]==0,];
ori_data_rep_011=ori_data_rep[ori_data_rep[,2]==0 & ori_data_rep[,3]==1 & ori_data_rep[,4]==1,];
ori_data_rep_010=ori_data_rep[ori_data_rep[,2]==0 & ori_data_rep[,3]==1 & ori_data_rep[,4]==0,];
ori_data_rep_001=ori_data_rep[ori_data_rep[,2]==0 & ori_data_rep[,3]==0 & ori_data_rep[,4]==1,];
ori_data_rep_000=ori_data_rep[ori_data_rep[,2]==0 & ori_data_rep[,3]==0 & ori_data_rep[,4]==0,];

prob_rep_mat[2,i]=mean(ori_data_rep_111[,1]);
prob_rep_mat[3,i]=mean(ori_data_rep_110[,1]);
prob_rep_mat[4,i]=mean(ori_data_rep_101[,1]);
prob_rep_mat[5,i]=mean(ori_data_rep_100[,1]);
prob_rep_mat[6,i]=mean(ori_data_rep_011[,1]);
prob_rep_mat[7,i]=mean(ori_data_rep_010[,1]);
prob_rep_mat[8,i]=mean(ori_data_rep_001[,1]);
prob_rep_mat[9,i]=mean(ori_data_rep_000[,1]);





}



# summarize the differences;
prob_diff_mat=prob_rep_mat-prob_est_mat;

rowMeans(prob_rep_mat);
rowMeans(prob_est_mat);

# posterior predictive P-values;
rowMeans(prob_diff_mat>0);



# posterior estiamtes of the bivariate probit models;
imp<-jomo1cat.MCMCchain(ori_data_factor, rep(2,4),nburn=nburn)
