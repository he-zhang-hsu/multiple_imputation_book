# Example 13.12
rm(list=ls());

library(sampleSelection)
library(R2WinBUGS);
library(mice);
library(norm);
options(digits=4);

# read the data;
a1c=read.csv(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter13_univariate_nonignorable\\program\\A1c.csv");

rowobs=nrow(a1c);

attach(a1c);

# BMI has two missing values;
sum(is.na(male));
sum(is.na(white));

# impute the two missing values in BMI;
BMI_mean=mean(BMI, na.rm=T);
BMI_sd=sd(BMI, na.rm=T);

set.seed(197789);
BMI[61]=rnorm(1, mean=BMI_mean, sd=BMI_sd);

set.seed(197552);
BMI[81]=rnorm(1, mean=BMI_mean, sd=BMI_sd);


z=rep(0,length(Age));
z[is.na(HgA1c)==FALSE]=1;

# log transformation of HgA1c;

HgA1c.o=log(HgA1c)


HgA1c.o[z==0]=0;

ln_hga1c=HgA1c.o;
ln_hga1c[ma1c==1]=NA;


# impute missing ln_hga1c using the selection model
# by setting the correlation coefficient to be different values;

mi_no=50;
MI_mean_mat=MI_mean_var_mat=BUGS_mean_mat=BUGS_mean_var_mat=rep(NA, mi_no);

# WinBUGS files;
# rho is set as -0.9;
model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_ne9.txt")
# rho is set as -0.7;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_ne7.txt")
# rho is set as -0.5;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_ne5.txt")
# rho is set as -0.3;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_ne3.txt")
# rho is set as -0.1;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_ne1.txt")
# rho is set as 0.1;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_po1.txt")
# rho is set as 0.3;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_po3.txt")
# rho is set as 0.5;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_po5.txt")

# rho is set as 0.7;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_po7.txt")


# rho is set as 0.9;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_po9.txt")


# Let's take a look:
file.show(model.file)

gibbs_no=20000;

process_stage2_data_simulate=list(n=rowobs, ma1c=ma1c, ln_hga1c=ln_hga1c, Age=Age, Diabetes=Diabetes, BMI=BMI, white=white);
process_stage2_init= function(){list(a = c(1.742, -0.003, 0.005, 0.253, -0.005), b = c(0.216, -0.005, 0.741, 0.016, 0.064), sigma=0.16)};
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



for (i in 1:mi_no)

{

y_completed_BUGS[ma1c==1]=y_impute[i,];

y_completed_BUGS_mat[,i]=y_completed_BUGS;

BUGS_mean_mat[i]=mean(exp(y_completed_BUGS));
BUGS_mean_var_mat[i]=var(exp(y_completed_BUGS))/rowobs;

}

# Fig. 13.5 with the selected rho
hist(y_completed_BUGS, xlab="Completed log(HgA1c) by selection model, rho=-.9", main="");

# combine the estimates for the mean;
# Table 13.7 with the selected rho
BUGS_summary=pool.scalar(BUGS_mean_mat, BUGS_mean_var_mat, n=rowobs);
BUGS_summary;
BUGS_summary$qbar
BUGS_summary$qbar-qt(.975, BUGS_summary$df)*sqrt(BUGS_summary$t);
BUGS_summary$qbar+qt(.975, BUGS_summary$df)*sqrt(BUGS_summary$t);
BUGS_summary$fmi;


