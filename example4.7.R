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

rowobs;
mis_no;
obs_no;

# create a binary outcome variable of carercvd;

carercvd_23=carercvd;
# collapse the category 2 and 3;
carercvd_23[carercvd_23>=2]=2;

# recode carercvd_23 into a binary indicator;
carercvd_23=carercvd_23-1;

miss_indi=is.na(carercvd_23);
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];


# fetch the covariates;
# gender;
sex=(cdat$sex);

table(sex);

# general health;
genhlth=cdat$genhlth;

# education level;
X_educag=cdat$X_educag;

table(X_educag);


table(genhlth);


# having health coverage;

hlthpln1=cdat$hlthpln1;

table(hlthpln1);


# having delayed health care;
delaym=cdat$delaym;

table(delaym);

# age group;
age_g=cdat$X_age_g;

table(age_g);

# race groups;

race_g=cdat$X_race_g1;

table(race_g);


# income group;

incomeg=cdat$X_incomg;

table(incomeg);

# employ1;

employ1=cdat$employ1;

table(employ1);

table(cdat);

obs_no;




# fully observed covariates;
covariates=cbind(sex, genhlth, age_g, X_educag, hlthpln1, delaym);

table(sex)/rowobs;
table(genhlth)/rowobs;
table(age_g)/rowobs;
table(X_educag)/rowobs;
table(hlthpln1)/rowobs;
table(delaym)/rowobs;
table(carercvd)/obs_no;

# descriptive statistics of covariates;

covariates_miss=covariates[miss_seq,];
covariates_obs=covariates[-miss_seq,];


# propensity model analysis;
# not missing completely at random;
propensity=glm(miss_indi~covariates,family=binomial(link="logit"));
summary(propensity);

table(miss_indi, delaym);

# logistic regression analysis for carercvd_23 using complete cases;
logistic_cc=glm(carercvd_23~covariates,family=binomial(link="logit"));
summary(logistic_cc);

# multinomial logistic regression model;

generallogistic_cc=multinom(carercvd~covariates);

summary(generallogistic_cc);

generallogistic_cc=multinom(carercvd~sex+as.factor(genhlth), data=cdat);

summary(generallogistic_cc);


# proportional odds model;

cdat$carercvd_factor <- factor(cdat$carercvd);

# Run the ordinal logistic regression model
ppods <- polr(carercvd_factor ~ sex+genhlth+age_g+X_educag+hlthpln1+delaym, data=cdat);

brant(ppods);

colMeans(covariates_obs);
colMeans(covariates_miss);


# impute polytomous variable;

mi_no=50;

y2_completed_glogit=carercvd;
y2_miss=carercvd;

p1_mean_mat=p1_var_mat=p2_mean_mat=p2_var_mat=p3_mean_mat=p3_var_mat=rep(NA, mi_no);


# general logit model;


for (i in 1:mi_no)
{
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
# mle imputation;
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

# proportional odds model;
# test=mice.impute.polr(as.factor(y2_miss), ry=as.logical(1-miss_indi), seed=197789, x=covariates);
# test=mice.impute.polr(test_data$V1, ry=as.logical(1-miss_indi), seed=197789, x=test_data$V2);
# u=polr(as.factor(carercvd)~sex, data=cdat);

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





# logistic regression imputation for carercvd_23;
# mi_no=50;

# marginal_mean_mat=marginal_var_mat=obs_mean_mat=obs_var_mat=trans_mean_mat=trans_var_mat=rep(NA, mi_no);


# y2_completed_mle=y2_completed_boot=y2_completed_bayes=carercvd_23;

# y2_miss=carercvd_23;

# for (i in 1:mi_no)
# {


# Rubin's mle imputation;
# y2_imputed_mle=mice.impute.logreg(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=covariates);

# bootstrap imputation;
# y2_imputed_boot=mice.impute.logreg.boot(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=covariates);


# imputation model based on exact MCMC algorithm;
# diffuse prior;
# coef_s1x_posterior=bayesglm(y2_miss~covariates, family=binomial(link="logit"), prior.scale=Inf, prior.df=Inf);  
# coef_s1x_draw=as.vector(sim(coef_s1x_posterior, n.sims=1)@coef);
# prob_predict=exp(cbind(rep(1,mis_no),covariates_miss)%*%coef_s1x_draw)/(1+exp(cbind(rep(1,mis_no),covariates_miss)%*%coef_s1x_draw));
# y2_imputed_bayes=rbinom(n=mis_no, size=1, prob=prob_predict);


# y2_completed_mle[miss_seq]=y2_imputed_mle;
# y2_completed_boot[miss_seq]=y2_imputed_boot;
# y2_completed_bayes[miss_seq]=y2_imputed_bayes;


# marginal means;
# obs_mean=mean(y2_completed_mle);
# obs_var=var(y2_completed_mle)/rowobs;

# obs_mean_mat[i]=obs_mean;
# obs_var_mat[i]=obs_var;

# logistic regression coefficient for slope;
# mle_logistic=summary(glm(y2_completed_mle~x, family=binomial(link="logit")));
# mle_slope_mat[cycle,i]=mle_logistic$coeff[2,1];
# mle_slope_var_mat[cycle,i]=mle_logistic$coeff[2,2]^2;



# trans_mean=mean(y2_completed_boot);
# trans_mean_var=var(y2_completed_boot)/rowobs;

trans_mean_mat[i]=trans_mean;
trans_var_mat[i]=trans_mean_var;



# logistic regression coefficient for slope;
# boot_logistic=summary(glm(y2_completed_boot~x, family=binomial(link="logit")));
# boot_slope_mat[cycle,i]=boot_logistic$coeff[2,1];
# boot_slope_var_mat[cycle,i]=boot_logistic$coeff[2,2]^2;


marginal_mean=mean(y2_completed_bayes);
marginal_mean_var=var(y2_completed_bayes)/rowobs;

marginal_mean_mat[i]=marginal_mean;
marginal_var_mat[i]=marginal_mean_var;

# logistic regression coefficient for slope;
# bayes_logistic=summary(glm(y2_completed_bayes~x, family=binomial(link="logit")));
# bayes_slope_mat[cycle,i]=bayes_logistic$coeff[2,1];
# bayes_slope_var_mat[cycle,i]=bayes_logistic$coeff[2,2]^2;


}

# complete-case analysis;
# point estimates;
mean(carercvd_23, na.rm=T);
# standard error;
sqrt(var(carercvd_23, na.rm=T)/obs_no);

# multiple imputation analysis;
# mle imputation;
obs_summary=pool.scalar(obs_mean_mat, obs_var_mat);
obs_mi_mean=obs_summary$qbar;

obs_mi_mean;

obs_mi_var=obs_summary$t;
sqrt(obs_mi_var);

obs_summary$fmi;


# bootstrap imputation;
trans_summary=pool.scalar(trans_mean_mat, trans_var_mat);
trans_mi_mean=trans_summary$qbar;

trans_mi_mean;

trans_mi_var=trans_summary$t;

sqrt(trans_mi_var);

trans_summary$fmi;

# bayes imputation;

marginal_summary=pool.scalar(marginal_mean_mat, marginal_var_mat);
marginal_mi_mean=marginal_summary$qbar;

marginal_mi_mean;

marginal_mi_var=marginal_summary$t;

sqrt(marginal_mi_var);

marginal_summary$fmi;


# collapse category 1 and 2
# create a binary outcome variable of carercvd;

carercvd_12=carercvd;
# collapse the category 1 and 2;
carercvd_12[carercvd_12<=2]=2;

# recode carercvd_12 into a binary indicator;
carercvd_12=carercvd_12-2;

# logistic regression analysis for carercvd_12 using complete cases;
logistic_cc=glm(carercvd_12~covariates,family=binomial(link="logit"));
summary(logistic_cc);


# logistic regression imputation for carercvd_12;
mi_no=50;

marginal_mean_mat=marginal_var_mat=obs_mean_mat=obs_var_mat=trans_mean_mat=trans_var_mat=rep(NA, mi_no);


y2_completed_mle=y2_completed_boot=y2_completed_bayes=carercvd_12;

y2_miss=carercvd_12;

for (i in 1:mi_no)
{


# Rubin's mle imputation;
y2_imputed_mle=mice.impute.logreg(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=covariates);

# bootstrap imputation;
y2_imputed_boot=mice.impute.logreg.boot(y2_miss, ry=as.logical(1-miss_indi), seed=i, x=covariates);


# imputation model based on exact MCMC algorithm;
# diffuse prior;
coef_s1x_posterior=bayesglm(y2_miss~covariates, family=binomial(link="logit"), prior.scale=Inf, prior.df=Inf);  
coef_s1x_draw=as.vector(sim(coef_s1x_posterior, n.sims=1)@coef);
prob_predict=exp(cbind(rep(1,mis_no),covariates_miss)%*%coef_s1x_draw)/(1+exp(cbind(rep(1,mis_no),covariates_miss)%*%coef_s1x_draw));
y2_imputed_bayes=rbinom(n=mis_no, size=1, prob=prob_predict);


y2_completed_mle[miss_seq]=y2_imputed_mle;
y2_completed_boot[miss_seq]=y2_imputed_boot;
y2_completed_bayes[miss_seq]=y2_imputed_bayes;


# marginal means;
obs_mean=mean(y2_completed_mle);
obs_var=var(y2_completed_mle)/rowobs;

obs_mean_mat[i]=obs_mean;
obs_var_mat[i]=obs_var;

# logistic regression coefficient for slope;
# mle_logistic=summary(glm(y2_completed_mle~x, family=binomial(link="logit")));
# mle_slope_mat[cycle,i]=mle_logistic$coeff[2,1];
# mle_slope_var_mat[cycle,i]=mle_logistic$coeff[2,2]^2;



trans_mean=mean(y2_completed_boot);
trans_mean_var=var(y2_completed_boot)/rowobs;

trans_mean_mat[i]=trans_mean;
trans_var_mat[i]=trans_mean_var;



# logistic regression coefficient for slope;
# boot_logistic=summary(glm(y2_completed_boot~x, family=binomial(link="logit")));
# boot_slope_mat[cycle,i]=boot_logistic$coeff[2,1];
# boot_slope_var_mat[cycle,i]=boot_logistic$coeff[2,2]^2;


marginal_mean=mean(y2_completed_bayes);
marginal_mean_var=var(y2_completed_bayes)/rowobs;

marginal_mean_mat[i]=marginal_mean;
marginal_var_mat[i]=marginal_mean_var;

# logistic regression coefficient for slope;
# bayes_logistic=summary(glm(y2_completed_bayes~x, family=binomial(link="logit")));
# bayes_slope_mat[cycle,i]=bayes_logistic$coeff[2,1];
# bayes_slope_var_mat[cycle,i]=bayes_logistic$coeff[2,2]^2;


}

# complete-case analysis;
# point estimates;
mean(carercvd_12, na.rm=T);
# standard error;
sqrt(var(carercvd_12, na.rm=T)/obs_no);

# multiple imputation analysis;
# mle imputation;
obs_summary=pool.scalar(obs_mean_mat, obs_var_mat);
obs_mi_mean=obs_summary$qbar;

obs_mi_mean;

obs_mi_var=obs_summary$t;
sqrt(obs_mi_var);

obs_summary$fmi;

# bootstrap imputation;
trans_summary=pool.scalar(trans_mean_mat, trans_var_mat);
trans_mi_mean=trans_summary$qbar;

trans_mi_mean;

trans_mi_var=trans_summary$t;

sqrt(trans_mi_var);

# bayes imputation;

trans_summary$fmi;

marginal_summary=pool.scalar(marginal_mean_mat, marginal_var_mat);
marginal_mi_mean=marginal_summary$qbar;

marginal_mi_mean;

marginal_mi_var=marginal_summary$t;

sqrt(marginal_mi_var);

marginal_summary$fmi;
