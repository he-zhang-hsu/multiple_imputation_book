# Example 4.7

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
library(nnet);
library(brant);


# read the data set from Zhou et al. 2017;
cdat<-read.xlsx("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter4_univariate_model_parametric\\program\\paper1_dataset.xlsx", sheet=1, startRow=1, colNames=T);

# 1430 cases;
nrow(cdat);

# fetch the outcome variable;
table(cdat$carercvd);

carercvd=as.numeric(cdat$carercvd);
carercvd_obs=carercvd[!is.na(carercvd)];

rowobs=length(carercvd);
mis_no=sum(is.na(carercvd));
obs_no=rowobs-mis_no;

# rowobs is the sample size of the dataset;
rowobs;
# mis_no is the number of missing cases in carercvd;
mis_no;
# obs_no is the number of observed cases in carercvd;
obs_no;

# create a binary outcome variable of carercvd;

carercvd_23=carercvd;
# collapse the category 2 and 3;
carercvd_23[carercvd_23>=2]=2;

# recode carercvd_23 into a binary indicator;
carercvd_23=carercvd_23-1;

# miss_indi is the missingness indicator;
miss_indi=is.na(carercvd_23);
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];


# fetch the covariates;
# gender;
sex=(cdat$sex);

# table(sex);

# general health;
genhlth=cdat$genhlth;

# education level;
X_educag=cdat$X_educag;

# table(X_educag);
# table(genhlth);


# having health coverage;

hlthpln1=cdat$hlthpln1;

# table(hlthpln1);

# having delayed health care;
delaym=cdat$delaym;

# table(delaym);

# age group;
age_g=cdat$X_age_g;

# table(age_g);

# race groups;

race_g=cdat$X_race_g1;

# table(race_g);

# income group;

incomeg=cdat$X_incomg;

# table(incomeg);

# employ1;

employ1=cdat$employ1;

# table(employ1);

# table(cdat);



# descriptive statistics of the dataset
# Table 4.10;
# fully observed covariates;
covariates=cbind(sex, genhlth, age_g, X_educag, hlthpln1, delaym);

table(sex)/rowobs;
table(genhlth)/rowobs;
table(age_g)/rowobs;
table(X_educag)/rowobs;
table(hlthpln1)/rowobs;
table(delaym)/rowobs;
table(carercvd)/obs_no;


# propensity model analysis;
# running a logistic regression model for the missingness indicator
# not missing completely at random;
# Table 4.11;
covariates_miss=covariates[miss_seq,];
covariates_obs=covariates[-miss_seq,];
propensity=glm(miss_indi~covariates,family=binomial(link="logit"));
summary(propensity);

table(miss_indi, delaym);

# logistic regression analysis for carercvd_23 using complete cases;
# logistic_cc=glm(carercvd_23~covariates,family=binomial(link="logit"));
# summary(logistic_cc);

# multinomial logistic regression model for carercvd;
# Table 4.12;

# Treat all covariates as continuous
generallogistic_cc=multinom(carercvd~covariates);
summary(generallogistic_cc);


# generallogistic_cc=multinom(carercvd~sex+as.factor(genhlth), data=cdat);
# summary(generallogistic_cc);


# proportional odds model;
# cdat$carercvd_factor <- factor(cdat$carercvd);

# Run the ordinal logistic regression model
# ppods <- polr(carercvd_factor ~ sex+genhlth+age_g+X_educag+hlthpln1+delaym, data=cdat);
# brant(ppods);

# colMeans(covariates_obs);
# colMeans(covariates_miss);


# multiply impute the polytomous variable;
# 50 imputations;
mi_no=50;

y2_completed_glogit=carercvd;
y2_miss=carercvd;

# p1/p2/p3 are the proportion estimates of the 3 categories of carercvd;
p1_mean_mat=p1_var_mat=p2_mean_mat=p2_var_mat=p3_mean_mat=p3_var_mat=rep(NA, mi_no);


# multinomial logistic regression imputation;

for (i in 1:mi_no)
{
set.seed(i);
y2_imputed_mle=mice.impute.polyreg(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=covariates);
y2_completed_glogit[miss_seq]=as.numeric(y2_imputed_mle);
p1=sum(y2_completed_glogit==1)/rowobs;
var_p1=p1*(1-p1)/rowobs;
p2=sum(y2_completed_glogit==2)/rowobs;
var_p2=p2*(1-p2)/rowobs;
p3=sum(y2_completed_glogit==3)/rowobs;
var_p3=p3*(1-p3)/rowobs;

p1_mean_mat[i]=p1;
p1_var_mat[i]=var_p1;
p2_mean_mat[i]=p2;
p2_var_mat[i]=var_p2;
p3_mean_mat[i]=p3;
p3_var_mat[i]=var_p3;

}


# Results in Table 4.13;
# complete case analysis for the 3 proportions;

p1_obs=sum(carercvd_obs==1)/obs_no;
se_p1_obs=sqrt(p1_obs*(1-p1_obs)/obs_no);
p2_obs=sum(carercvd_obs==2)/obs_no;
se_p2_obs=sqrt(p2_obs*(1-p2_obs)/obs_no);
p3_obs=sum(carercvd_obs==3)/obs_no;
se_p3_obs=sqrt(p3_obs*(1-p3_obs)/obs_no);

p1_obs;
se_p1_obs;
p2_obs;
se_p2_obs;
p3_obs;
se_p3_obs;

# combine the multiple imputation estimates;
p1_summary=pool.scalar(p1_mean_mat, p1_var_mat, n=rowobs);
p1_mi_mean=p1_summary$qbar;

p1_mi_mean;

p1_mi_var=p1_summary$t;
sqrt(p1_mi_var);

p1_summary$fmi;

p2_summary=pool.scalar(p2_mean_mat, p2_var_mat, n=rowobs);
p2_mi_mean=p2_summary$qbar;

p2_mi_mean;

p2_mi_var=p2_summary$t;
sqrt(p2_mi_var);

p2_summary$fmi;

p3_summary=pool.scalar(p3_mean_mat, p3_var_mat, n=rowobs);
p3_mi_mean=p3_summary$qbar;

p3_mi_mean;

p3_mi_var=p3_summary$t;
sqrt(p3_mi_var);

p3_summary$fmi;

# proportional odds model imputation;
# test=mice.impute.polr(as.factor(y2_miss), ry=as.logical(1-miss_indi), seed=197789, x=covariates);
# test=mice.impute.polr(test_data$V1, ry=as.logical(1-miss_indi), seed=197789, x=test_data$V2);
# u=polr(as.factor(carercvd)~sex, data=cdat);

y2_miss=carercvd;
test_data=as.data.frame(cbind(as.factor(y2_miss), covariates));
test_data$V1=as.factor(test_data$V1);
u=mice(data=test_data, m=50, method=c('polr'), seed=197789);

y2_completed_prop=carercvd;

p1_mean_mat=p1_var_mat=p2_mean_mat=p2_var_mat=p3_mean_mat=p3_var_mat=rep(NA, mi_no);

for (i in 1:mi_no)
{
y2_completed_prop[miss_seq]=as.numeric(u$imp$V1[,i]);
p1=sum(y2_completed_prop==1)/rowobs;
var_p1=p1*(1-p1)/rowobs;
p2=sum(y2_completed_prop==2)/rowobs;
var_p2=p2*(1-p2)/rowobs;
p3=sum(y2_completed_prop==3)/rowobs;
var_p3=p3*(1-p3)/rowobs;

p1_mean_mat[i]=p1;
p1_var_mat[i]=var_p1;
p2_mean_mat[i]=p2;
p2_var_mat[i]=var_p2;
p3_mean_mat[i]=p3;
p3_var_mat[i]=var_p3;

}


# combine the multiple imputation estimates;
# proportional odds imputation;
# The results are similar to those from the 
# multinomial logistic regression imputation;
p1_summary=pool.scalar(p1_mean_mat, p1_var_mat, n=rowobs);
p1_mi_mean=p1_summary$qbar;

p1_mi_mean;

p1_mi_var=p1_summary$t;
sqrt(p1_mi_var);

p1_summary$fmi;

p2_summary=pool.scalar(p2_mean_mat, p2_var_mat, n=rowobs);
p2_mi_mean=p2_summary$qbar;

p2_mi_mean;

p2_mi_var=p2_summary$t;
sqrt(p2_mi_var);

p2_summary$fmi;

p3_summary=pool.scalar(p3_mean_mat, p3_var_mat, n=rowobs);
p3_mi_mean=p3_summary$qbar;

p3_mi_mean;

p3_mi_var=p3_summary$t;
sqrt(p3_mi_var);

p3_summary$fmi;



