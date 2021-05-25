
rm(list=ls());
library(mice);
library(norm);
# library(HI);
options(digits=4);

# Example 4.2;
# Read the gestational age data from Example 4.2

data=read.table(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter4_univariate_model_parametric\\program\\dgestat.txt", col.names=c("dgestat", "dbirwt", "dmage", "dfage", "fmaps", "wtgain", "nprevis", "cigar", "drink"), na.strings=".");

# fetch the outcome;
y=data$dgestat;
dbirwt=data$dbirwt;

miss_indi=is.na(y);

# missingness rate;
mean(miss_indi);
total_no=nrow(data);
mis_no=sum(miss_indi);
obs_no=total_no-mis_no;

y_miss=y;
miss_seq=cbind(seq(1,length(y),1), miss_indi)[miss_indi==1,1];

# descriptive statistics;

# summary(data);

# Table 4.3;
colMeans(data, na.rm=T);

sd(data$dgestat,na.rm=T);
sd(data$dbirwt, na.rm=T);
sd(data$dmage, na.rm=T);
sd(data$dfage, na.rm=T);
sd(data$fmaps, na.rm=T);
sd(data$wtgain, na.rm=T);
sd(data$nprevis, na.rm=T);
sd(data$cigar, na.rm=T);
sd(data$drink, na.rm=T);


# dbirtw_w39=dbirwt[!is.na(y) & (y==39)];
# hist(dbirtw_w39);



# fetch the predictors;
x=data[,2:9];

# Multiple imputation for missing gestational age using normal linear model;
mi_no=20;

y_completed=y_miss;

y_imputed_mat=matrix(NA, mis_no, mi_no);

for (i in 1:mi_no)
{
# proper model imputation;
set.seed(i);
y_imputed=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);
y_completed[miss_seq]=y_imputed;
y_imputed_mat[,i]=y_imputed;
}

# Histogram plots for example 4.2
# Fig. 4.3

# Plot observed values;
hist(y_miss, xlim=c(19, 48), xlab="DGESTAT", main="Observed", breaks=seq(19,48,1))

# Plot the imputed missing values from the last imputation;
hist(y_imputed_mat[,20], xlim=c(19,48), xlab="DGESTAT", breaks=seq(19,48,1), main="Normal Imputation", add=F)

