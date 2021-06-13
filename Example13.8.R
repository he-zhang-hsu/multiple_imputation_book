# Example 13.8

rm(list=ls());
library(sampleSelection)
library(R2WinBUGS);
library(mice);
library(norm);
options(digits=4);

# read the data;

a1c=read.csv(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter13_univariate_nonignorable\\program\\A1c.csv");

# sample size;
rowobs=nrow(a1c);

attach(a1c);

# right-skewed distribution;
# Fig. 13.3 top left;
# HgA1c is same as HbA1c;
hist(HgA1c, main="", xlab="HbA1c");

mean(HgA1c[Diabetes==1], na.rm=T);
mean(HgA1c[Diabetes==0], na.rm=T);

# Fig. 13.3 top right
hist(log(HgA1c), main="", xlab="log(HbA1c)");

# ma1c is the missingness indicator;

mean(ma1c);

sum(is.na(Age));
sum(is.na(Diabetes));
sum(is.na(BMI));

# BMI has two missing values;
sum(is.na(male));
sum(is.na(white));

# impute the two missing values in BMI using some simple procedures;
BMI_mean=mean(BMI, na.rm=T);
BMI_sd=sd(BMI, na.rm=T);

set.seed(197789);
BMI[61]=rnorm(1, mean=BMI_mean, sd=BMI_sd);

set.seed(197552);
BMI[81]=rnorm(1, mean=BMI_mean, sd=BMI_sd);

# the association between ma1c and male;
table(ma1c, male);


# descriptive statistics of covariates;
mean(Age);
sd(Age);

mean(Diabetes);

mean(BMI);
sd(BMI);

mean(white);

# Table 13.3
# propensity score analysis;
propensity=glm(1-ma1c~Age+Diabetes+BMI+white, family=binomial(link="probit"))

summary(propensity);

z=rep(0,length(Age));
z[is.na(HgA1c)==FALSE]=1;

# log transformation of HgA1c;

HgA1c.o=log(HgA1c)

# linear regression model on the transformed scale;
linear=lm(HgA1c.o ~ Age+Diabetes+BMI+white);
summary(linear);

# linear regression model on the original scale;
# linear_untran=lm(HgA1c ~ Age+Diabetes+BMI+white);
# summary(linear_untran);

# Heckman selection model fitting;
HgA1c.o[z==0]=0;

selection_fit=selection(z~Age+Diabetes+BMI+white,HgA1c.o~Age+BMI+Diabetes+white)

summary(selection_fit);

# Direct imputation under the selection model;

ln_hga1c=HgA1c.o;
ln_hga1c[ma1c==1]=NA;

mi_no=50;

# impute missing ln_hga1c using the selection model;
MI_mean_mat=MI_mean_var_mat=BUGS_mean_mat=BUGS_mean_var_mat=rep(NA, mi_no);

# WinBUGS file;

model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c.txt")
# Let's take a look:
file.show(model.file)

gibbs_no=20000;
mi_no=50;

process_stage2_data_simulate=list(n=rowobs, ma1c=ma1c, ln_hga1c=ln_hga1c, Age=Age, Diabetes=Diabetes, BMI=BMI, white=white);
process_stage2_init= function(){list(a = c(1.742, -0.003, 0.005, 0.253, -0.005), b = c(0.216, -0.005, 0.741, 0.016, 0.064), theta=0.15, sigma=0.16)};
process_stage2_parameters=c("y_impute");

process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file,
    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=gibbs_no/(2*mi_no), bugs.seed=197789,
    bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);

attach.bugs(process_stage2.sim);
rm(process_stage2.sim);

y_miss=ln_hga1c;

y_completed_BUGS=y_completed_MAR=y_miss;

y_completed_BUGS_mat=y_completed_MAR_mat=matrix(NA, nrow=rowobs, ncol=mi_no);

# combine the multiple imputation results;

for (i in 1:mi_no)

{

y_completed_BUGS[ma1c==1]=y_impute[i,];

y_completed_BUGS_mat[,i]=y_completed_BUGS;

BUGS_mean_mat[i]=mean(exp(y_completed_BUGS));
BUGS_mean_var_mat[i]=var(exp(y_completed_BUGS))/rowobs;

}

# Fig. 13.3 bottom left;
hist(y_completed_BUGS, xlab="Completed log(HbA1c) by selection model", main="");

# MAR imputation;

for (i in 1:mi_no)

{

y_imputed_MI=mice.impute.norm(y_miss, ry=as.logical(1-ma1c), seed=i, x=cbind(Age, Diabetes, BMI, white));

y_completed_MAR[ma1c==1]=y_imputed_MI;

y_completed_MAR_mat[,i]=y_completed_MAR;

MI_mean_mat[i]=mean(exp(y_completed_MAR));
MI_mean_var_mat[i]=var(exp(y_completed_MAR))/rowobs;



}

# Fig. 13.3 bottom right;
hist(y_completed_MAR, xlab="Completed log(HbA1c) by MAR", main="");

# combine the estimates;
BUGS_summary=pool.scalar(BUGS_mean_mat, BUGS_mean_var_mat, n=rowobs);
MI_summary=pool.scalar(MI_mean_mat, MI_mean_var_mat, n=rowobs);

BUGS_summary;
MI_summary;

# Table 13.4;

BUGS_summary$qbar
BUGS_summary$qbar-qt(.975, BUGS_summary$df)*sqrt(BUGS_summary$t);
BUGS_summary$qbar+qt(.975, BUGS_summary$df)*sqrt(BUGS_summary$t);
BUGS_summary$fmi;

MI_summary$qbar;
MI_summary$qbar-qt(.975, MI_summary$df)*sqrt(MI_summary$t);
MI_summary$qbar+qt(.975, MI_summary$df)*sqrt(MI_summary$t);
MI_summary$fmi;

# complete-case mean;
obs_no=rowobs-sum(ma1c);

CC_mean=mean(HgA1c, na.rm=T);

CC_mean-1.96*sd(HgA1c, na.rm=T)/sqrt(obs_no);
CC_mean+1.96*sd(HgA1c, na.rm=T)/sqrt(obs_no);

# Bayesian estimates for the model parameters;
# Table 13.3

process_stage2_data_simulate=list(n=rowobs, ma1c=ma1c, ln_hga1c=ln_hga1c, Age=Age, Diabetes=Diabetes, BMI=BMI, white=white);
process_stage2_init= function(){list(a = c(1.742, -0.003, 0.005, 0.253, -0.005), b = c(0.216, -0.005, 0.741, 0.016, 0.064), theta=0.15, sigma=0.16)};
process_stage2_parameters=c("a.cons", "a.Age", "a.Diabetes", "a.BMI", "a.white","b.cons","b.Age","b.Diabetes","b.BMI","b.white","sigma","rho");

burn_in=5000;
total=burn_in+100*1000;

process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file,
    n.chains=1, n.iter=total, n.burnin=burn_in, n.thin=(total-burn_in)/1000, bugs.seed=197789,
    bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);

attach.bugs(process_stage2.sim);
rm(process_stage2.sim);

mean(a.cons);
sd(a.cons);
mean(a.Age);
sd(a.Age);
mean(a.Diabetes);
sd(a.Diabetes);
mean(a.BMI);
sd(a.BMI);
mean(a.white);
sd(a.white);

mean(b.cons);
sd(b.cons);
mean(b.Age);
sd(b.Age);
mean(b.Diabetes);
sd(b.Diabetes);
mean(b.BMI);
sd(b.BMI);
mean(b.white);
sd(b.white);

mean(sigma);
sd(sigma);
mean(rho);
sd(rho);



