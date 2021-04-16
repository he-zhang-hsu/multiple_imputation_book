time.begin=date();
library(pan);
library(MASS);


# begin the main program;

# begin the replication;

# for data set with age0-12 yrs old children, the subject no. is 2172;
# for data set with age3-12 yrs old children, the subject no. is 1793.

# import the missingdata set;
# oridata=scan(file=paste("//cdc.gov/private/L728/wdq7/book/chapter 11/psid_uni_uni_hr_2172_log.dat"), na.strings=".");

oridata=scan(file=paste("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter9_longitudinal\\program\\psid_uni_uni_hr_2172_log.dat"), na.strings=".");



orimatrix=matrix(oridata, nrow=2172, ncol=34, byrow=TRUE);

# oriy concatenate the response vector;
oriy=orimatrix[,13:31];

# plot the distribution of oriy;

hist(oriy, xlab="LOGINR", main="");
qqnorm(oriy, main="QQ plot of LOGINR");
qqline(oriy);


oriy_mean=colMeans(oriy, na.rm=T);
yplot=rbind(oriy[1:25,], oriy_mean);

# plot the 1st 25 subjects and the mean response;

time=seq(1,19,1);
year=seq(from=1979, to=1997, by=1);


matplot(year, t(yplot), col=c(rep("black", 25), "red"), lty=1, lwd=c(rep(1,25),3),type="l", ylab="LOGINR");

# plot the imputations;

# imputation.plot;

# plot the sample trajectory plot;
sample=scan(file=("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter9_longitudinal\\program\\individual.dat"), na.strings=".");

sample_matrix=matrix(sample, nrow=19, ncol=20);
matplot(x=year, y=sample_matrix, type="l", lty=rep(1,19), col=rep(1,19), ylab="LOGINR", main="Trajectory plot");


part_combine_impu=scan(file=("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter9_longitudinal\\program\\part_combine_impu.dat"), what=list(subject=0, impu=0, individual=0, imputed=0, observed=0, year=0), na.strings=".");

combine=cbind(part_combine_impu$subject, part_combine_impu$impu, part_combine_impu$individual,
part_combine_impu$imputed, part_combine_impu$observed, part_combine_impu$year);

data_11=combine[1:19,];

# plot 11;
 matplot(x=data_11[,6], y=data_11[,3:5], ylim=c(0,4), lty=c(1,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");


data_12=combine[20:38,];

# plot 12;
 matplot(x=data_12[,6], y=data_12[,3:5], ylim=c(0,4), lty=c(1,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");

data_13=combine[39:57,];

# plot 13;
 matplot(x=data_13[,6], y=data_13[,3:5], ylim=c(0,4), lty=c(1,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");


data_21=combine[58:76,];

# plot 21;
 matplot(x=data_21[,6], y=data_21[,3:5], ylim=c(0,4), lty=c(1,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");

data_22=combine[77:95,];

# plot 22;
matplot(x=data_22[,6], y=data_22[,3:5], ylim=c(0,4), lty=c(1,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");

data_23=combine[96:114,];

# plot 23;
matplot(x=data_23[,6], y=data_23[,3:5], ylim=c(0,4), lty=c(1,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");

data_31=combine[115:133,];

# plot 31;
matplot(x=data_31[,6], y=data_31[,3:5], ylim=c(0,4), lty=c(1,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");

data_32=combine[134:152,];

# plot 32;
matplot(x=data_32[,6], y=data_32[,3:5], ylim=c(0,4), lty=c(1,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");

data_33=combine[153:171,];

# plot 33;
matplot(x=data_33[,6], y=data_33[,3:5], ylim=c(0,4), lty=c(1,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");






oriyvec=as.vector(t(oriy));

# subj is the vector that function pan needs;
subj=rep(1:2172, each=19);

# xt is the time covariate;
xt=cbind(rep(1,19),seq(from=1, to=19, by=1));

# zerot is the zero covariate;
zerot=matrix(0, nrow=20, ncol=2);

# pred is the predictor covaraite matrix that function pan needs;
pred=matrix(0, nrow=2172*19, ncol=9);

for(i in 1:2172)
{
 if (orimatrix[i,32]==0 & orimatrix[i,33]==0 & orimatrix[i,34]==0 )
 {
 pred[((i-1)*19+1):(i*19),1:2]=xt;
 }
 if (orimatrix[i,32]==0 & orimatrix[i,33]==0 & orimatrix[i,34]==1 )
 {
 pred[((i-1)*19+1):(i*19),1:2]=xt;
 pred[((i-1)*19+1):(i*19),7]=rep(1,19);
 }
 if (orimatrix[i,32]==0 & orimatrix[i,33]==1 & orimatrix[i,34]==0 )
 {
 pred[((i-1)*19+1):(i*19),1:2]=xt;
 pred[((i-1)*19+1):(i*19),5:6]=xt;
 }
 if (orimatrix[i,32]==0 & orimatrix[i,33]==1 & orimatrix[i,34]==1 )
 {
 pred[((i-1)*19+1):(i*19),1:2]=xt;
 pred[((i-1)*19+1):(i*19),5:6]=xt;
 pred[((i-1)*19+1):(i*19),7]=rep(1,19);
 }
 if (orimatrix[i,32]==1 & orimatrix[i,33]==0 & orimatrix[i,34]==0 )
 {
 pred[((i-1)*19+1):(i*19),1:2]=xt;
 pred[((i-1)*19+1):(i*19),3:4]=xt;
 }
 if (orimatrix[i,32]==1 & orimatrix[i,33]==0 & orimatrix[i,34]==1 )
 {
 pred[((i-1)*19+1):(i*19),1:2]=xt;
 pred[((i-1)*19+1):(i*19),3:4]=xt;
 pred[((i-1)*19+1):(i*19),7]=rep(1,19);
 }
 if (orimatrix[i,32]==1 & orimatrix[i,33]==1 & orimatrix[i,34]==0 )
 {
 pred[((i-1)*19+1):(i*19),1:2]=xt;
 pred[((i-1)*19+1):(i*19),3:4]=xt;
 pred[((i-1)*19+1):(i*19),5:6]=xt;
 }
 if (orimatrix[i,32]==1 & orimatrix[i,33]==1 & orimatrix[i,34]==1 )
 {
 pred[((i-1)*19+1):(i*19),1:2]=xt;
 pred[((i-1)*19+1):(i*19),3:4]=xt;
 pred[((i-1)*19+1):(i*19),5:6]=xt;
 pred[((i-1)*19+1):(i*19),7]=rep(1,19);
 }
 pred[((i-1)*19+1):(i*19),8:9]=xt;
 
}

# xcol indicates the fixed effects covariate;
xcol=1:7;

# zcol indicates the random effects covariate;
zcol=8:9;

# prior is the hyperparameters that for sigma and random effects;

prior=list(a=1, Binv=1, c=5, Dinv=5*diag(2));

# begin imputation;
# the first imputation;

for (M in 1:5)
{

result=pan(y=oriyvec, subj=subj, pred=pred, xcol=xcol, zcol=zcol, prior=prior, seed=M, iter=3000);

imp=result$y;
impmatrix1=matrix(imp, nrow=2172, ncol=19, byrow=TRUE);
impmatrix=cbind(orimatrix[,1:12], impmatrix1, orimatrix[,32:34]);

# output the data set;
# determine the setnumber;

set.number=M;
options(digits=3);
# write.table(impmatrix, file=paste("h:\\files_to_copy\\yulei_thesis\\chapter5\\data\\psid\\imputed_data\\psid_uni_uni_2172_log_lin_impu", set.number, ".dat",
# sep=""), row.name=FALSE, col.name=FALSE, na=".");
#cat ("export the", M, "th imputed data set for missing data set", "\n");
# end of multiple imputation;
}

# some descriptive statistics;
time=seq(1,19,1);

matplot(time, t(oriy[1:5,]));

summary(oriyvec);

oriy_binary=1*(oriy>=2);
oriy_binary_vec=as.vector(t(oriy_binary));

# output the dataset;
# write.table(oriy_binary_vec, file=paste("//cdc.gov/private/L728/wdq7/book/chapter 11/psid_uni_uni_hr_2172_log_binary.dat",sep=""), row.name=FALSE, col.name=FALSE);

varname=c("Y[,1]", "Y[,2]", "Y[,3]","Y[,4]","Y[,5]","Y[,6]","Y[,7]","Y[,8]","Y[,9]","Y[,10]","Y[,11]","Y[,12]","Y[,13]","Y[,14]","Y[,15]","Y[,16]","Y[,17]","Y[,18]","Y[,19]");
 
# write.table(oriy_binary, file=paste("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter11_multileveldata\\program\\psid_uni_uni_hr_2172_log_binary.dat",sep=""), row.name=FALSE, col.name=varname);

# write.table(oriy_binary[1:50,], file=paste("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter11_multileveldata\\program\\psid_uni_uni_hr_2172_log_binary_part.dat",sep=""), row.name=FALSE, col.name=varname);

write.table(oriy, file=paste("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter11_multileveldata\\program\\psid_uni_uni_hr_2172_log.dat",sep=""), row.name=FALSE, col.name=varname);

# output the dataset with baseline covariates;

varname_cov=c("Y[,1]", "Y[,2]", "Y[,3]","Y[,4]","Y[,5]","Y[,6]","Y[,7]","Y[,8]","Y[,9]","Y[,10]","Y[,11]","Y[,12]","Y[,13]","Y[,14]","Y[,15]","Y[,16]","Y[,17]","Y[,18]","Y[,19]", "health[]", "edu[]", "race[]");

write.table(orimatrix[,13:34], file=paste("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter11_multileveldata\\program\\psid_uni_uni_hr_2172_log_cov.dat",sep=""), row.name=FALSE, col.name=varname_cov);

# assign some covariate observations as missing;
cov_matrix=orimatrix[,32:34];
health=orimatrix[,32];
edu=orimatrix[,33];
race=orimatrix[,34];

mean(race);
glm(edu ~ race);
glm(health ~ edu + race + edu*race);

rowobs=nrow(orimatrix);

# creat 10% missing data for selfhealth, coverage, and delaym;
set.seed(197789);

miss_indi=matrix(rbinom(n=3*rowobs, size=1, prob=0.1), nrow=rowobs, ncol=3);

# set some of the ori_data as missing values;
cov_matrix[miss_indi==1]=NA;

write.table(cbind(oriy, cov_matrix), file=paste("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter11_multileveldata\\program\\psid_uni_uni_hr_2172_log_cov_missing.dat",sep=""), row.name=FALSE, col.name=varname_cov);



oriy_binary_vec=as.vector(t(oriy_binary[1:50,]));
subj_binary=rep(1:50, each=19);
time_binary=rep(seq(1,19,1)-10, 50);

binary_data=cbind(subj_binary, oriy_binary_vec, time_binary);
write.table(binary_data, file=paste("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter11_multileveldata\\program\\test_binary_part.dat",sep=""), row.name=FALSE, col.name=c("subj", "Y", "Time"));


test=rep(0.125,100);

time.end=date();
