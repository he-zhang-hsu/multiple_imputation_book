# Example 9.4


time.begin=date();
library(pan);
library(MASS);



# import the missingdata set;

oridata=scan(file=paste("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter9_longitudinal\\program\\psid_uni_uni_hr_2172_log.dat"), na.strings=".");

orimatrix=matrix(oridata, nrow=2172, ncol=34, byrow=TRUE);

# oriy concatenate the response vector;
oriy=orimatrix[,13:31];

# plot the distribution of oriy;
# oriy is already log transformed;
# Fig. 9.2

hist(oriy, xlab="LOGINR", main="");

qqnorm(oriy, main="QQ plot of LOGINR");
qqline(oriy);

# plot the 1st 25 subjects and the mean response;
# Fig. 9.3
oriy_mean=colMeans(oriy, na.rm=T);
yplot=rbind(oriy[1:25,], oriy_mean);

time=seq(1,19,1);
year=seq(from=1979, to=1997, by=1);

matplot(year, t(yplot), col=c(rep("black", 25), "red"), lty=1, lwd=c(rep(1,25),3),type="l", ylab="LOGINR");

# plot the imputations;

# Fig. 9.4
part_combine_impu=scan(file=("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter9_longitudinal\\program\\part_combine_impu.dat"), what=list(subject=0, impu=0, individual=0, imputed=0, observed=0, year=0), na.strings=".");

combine=cbind(part_combine_impu$subject, part_combine_impu$impu, part_combine_impu$individual,
part_combine_impu$imputed, part_combine_impu$observed, part_combine_impu$year);

# Top row;
data_11=combine[1:19,];

# plot 11;
 matplot(x=data_11[,6], y=data_11[,3:5], ylim=c(0,4), lty=c(2,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");


data_12=combine[20:38,];

# plot 12;
 matplot(x=data_12[,6], y=data_12[,3:5], ylim=c(0,4), lty=c(2,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");

data_13=combine[39:57,];

# plot 13;
 matplot(x=data_13[,6], y=data_13[,3:5], ylim=c(0,4), lty=c(2,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");

# Middle row;
data_21=combine[58:76,];

# plot 21;
 matplot(x=data_21[,6], y=data_21[,3:5], ylim=c(0,4), lty=c(2,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");

data_22=combine[77:95,];

# plot 22;
matplot(x=data_22[,6], y=data_22[,3:5], ylim=c(0,4), lty=c(2,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");

data_23=combine[96:114,];

# plot 23;
matplot(x=data_23[,6], y=data_23[,3:5], ylim=c(0,4), lty=c(2,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");

# Bottom row;
data_31=combine[115:133,];

# plot 31;
matplot(x=data_31[,6], y=data_31[,3:5], ylim=c(0,4), lty=c(2,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");

data_32=combine[134:152,];

# plot 32;
matplot(x=data_32[,6], y=data_32[,3:5], ylim=c(0,4), lty=c(2,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");

data_33=combine[153:171,];

# plot 33;
matplot(x=data_33[,6], y=data_33[,3:5], ylim=c(0,4), lty=c(2,100,100),type="lpp", pch=c(1,17,1), col=c("blue","red","black"), lwd=c(1,3,3), xlab="", ylab="", main="");


