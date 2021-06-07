# Example 12.11

rm(list=ls());
# library(car);
# library(mvtnorm);
library(mice);
library(norm);
# library(HI);
options(digits=4);

# read the data
data=read.table(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter4_univariate_model_parametric\\program\\dgestat.txt", col.names=c("dgestat", "dbirwt", "dmage", "dfage", "fmaps", "wtgain", "nprevis", "cigar", "drink"), na.strings=".");

# fetch the outcome;
y=data$dgestat;
miss_indi=is.na(y);

# missingness rate;
mean(miss_indi);
rowobs=total_no=nrow(data);
mis_no=sum(miss_indi);
obs_no=total_no-mis_no;

y_miss=y;
miss_seq=cbind(seq(1,length(y),1), miss_indi)[miss_indi==1,1];

# descriptive statistics;

summary(data);
table(data$dgestat);

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



# descriptive statistics;
plot(data$dbirwt, data$dgestat, xlab="DBIRWT", ylab="DGESTAT", main="Observed");
# plot(data$dgestat, log(data$dbirwt), xlab="DBIRWT", ylab="DGESTAT", main="Observed");
plot(data$dmage, data$dgestat);
plot(data$dfage, data$dgestat);
plot(data$fmaps, data$dgestat);
plot(data$wtgain, data$dgestat);
plot(data$nprevis, data$dgestat);
plot(data$cigar, data$dgestat);
plot(data$drink, data$dgestat);

dbirwt=data$dbirwt;
dmage=data$dmage;
dfage=data$dfage;
fmaps=data$fmaps;
wtgain=data$wtgain;
nprevis=data$nprevis;
cigar=data$cigar;
drink=data$drink;

# fetch the predictors;
# with birth weight
x=data[,2:9];

# without birth weight;
x_nodbirwt=data[,3:9];


# Model I
# normal imputation without DBIRWT;

mi_no=100;

# matrix holding mean and variance estimates of multiple imputation;
MI_mean_mat=MI_mean_var_mat=rep(NA, mi_no);

y_completed_nm=y_miss;

for (i in 1:mi_no)
{
set.seed(i);
y_imputed_nm=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x_nodbirwt);
y_completed_nm[miss_seq]=y_imputed_nm;
MI_mean_mat[i]=mean(y_completed_nm);
MI_mean_var_mat[i]=var(y_completed_nm)/rowobs;
}

mean_summary=pool.scalar(MI_mean_mat, MI_mean_var_mat, n=rowobs);

# Table 12.6
mean_summary$qbar;
sqrt(mean_summary$t);
mean_summary$fmi;

# model II
# normal imputation with DBIRWT;

mi_no=100;


# matrix holding mean and variance estimates of multiple imputation;
MI_mean_mat=MI_mean_var_mat=rep(NA, mi_no);


y_completed_nm=y_miss;

for (i in 1:mi_no)
{
# proper model imputation;
set.seed(i);
y_imputed_nm=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);
y_completed_nm[miss_seq]=y_imputed_nm;
MI_mean_mat[i]=mean(y_completed_nm);
MI_mean_var_mat[i]=var(y_completed_nm)/rowobs;
}

# Table 12.6
mean_summary=pool.scalar(MI_mean_mat, MI_mean_var_mat, n=rowobs);
mean_summary$qbar;
sqrt(mean_summary$t);
mean_summary$fmi;


# Model III;
# qudratic normal imputation;
dbirwt_std=(dbirwt-mean(dbirwt))/sd(dbirwt);
dbirwt_std_square=dbirwt_std^2;

x_new=cbind(dbirwt_std, dbirwt_std_square, data[,3:9]);

mi_no=100;

y_completed_nm_square=y_miss;


# matrix holding mean and variance estimates of multiple imputation;
MI_mean_mat=MI_mean_var_mat=rep(NA, mi_no);

for (i in 1:mi_no)
{
set.seed(i);
y_imputed_nm_square=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x_new);
y_completed_nm_square[miss_seq]=y_imputed_nm_square;
MI_mean_mat[i]=mean(y_completed_nm_square);
MI_mean_var_mat[i]=var(y_completed_nm_square)/rowobs;

}

mean_summary=pool.scalar(MI_mean_mat, MI_mean_var_mat, n=rowobs);

# Table 12.6
mean_summary$qbar;
sqrt(mean_summary$t);
mean_summary$fmi;



# Model IV;
# quadratic PMM imputation;

mi_no=100;

y_completed_pmm_square=y_miss;


# matrix holding mean and variance estimates of multiple imputation;
MI_mean_mat=MI_mean_var_mat=rep(NA, mi_no);



for (i in 1:mi_no)
{
# proper model imputation;
set.seed(i);
y_imputed_pmm_square=mice.impute.pmm(y_miss, ry=as.logical(1-miss_indi), donors=5, seed=i, x=x_new);
y_completed_pmm_square[miss_seq]=y_imputed_pmm_square;
MI_mean_mat[i]=mean(y_completed_pmm_square);
MI_mean_var_mat[i]=var(y_completed_pmm_square)/rowobs;


}

mean_summary=pool.scalar(MI_mean_mat, MI_mean_var_mat, n=rowobs);

# Table 12.6
mean_summary$qbar;
sqrt(mean_summary$t);
mean_summary$fmi;

