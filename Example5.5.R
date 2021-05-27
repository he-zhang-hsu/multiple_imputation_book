# Example 5.5

rm(list=ls());

# library(car);
# library(mvtnorm);
library(mice);
library(norm);
# library(HI);
options(digits=4);


data=read.table(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter4_univariate_model_parametric\\program\\dgestat.txt", col.names=c("dgestat", "dbirwt", "dmage", "dfage", "fmaps", "wtgain", "nprevis", "cigar", "drink"), na.strings=".");

# fetch the outcome;
y=data$dgestat;
miss_indi=is.na(y);

# missingness rate;
mean(miss_indi);
total_no=nrow(data);
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

# Scatter plots;
# Fig. 5.10, the top panel;

plot(data$dbirwt, data$dgestat, xlab="DBIRWT", ylab="DGESTAT", main="Observed");
# plot(data$dgestat, log(data$dbirwt), xlab="DBIRWT", ylab="DGESTAT", main="Observed");
# plot(data$dmage, data$dgestat);
# plot(data$dfage, data$dgestat);
# plot(data$fmaps, data$dgestat);
# plot(data$wtgain, data$dgestat);
# plot(data$nprevis, data$dgestat);
# plot(data$cigar, data$dgestat);
# plot(data$drink, data$dgestat);

dbirwt=data$dbirwt;
dmage=data$dmage;
dfage=data$dfage;
fmaps=data$fmaps;
wtgain=data$wtgain;
nprevis=data$nprevis;
cigar=data$cigar;
drink=data$drink;


# fetch the predictors for imputation;
x=data[,2:9];


# run the propensity score analysis;
# dgestat_ps=glm(miss_indi ~ dbirwt+dmage+dfage+fmaps+wtgain+nprevis+cigar+drink,  family=binomial(link="logit"));

# summary(dgestat_ps);

# some exploratory analyses;

# dgestat_cen=y-mean(y, na.rm=T);

# dbirwt_pred=lm(dbirwt ~ dgestat_cen + I(dgestat_cen*dgestat_cen) + dmage + dfage + fmaps + wtgain + nprevis + cigar + drink);
# summary(dbirwt_pred);

# model 1
# normal model imputation;

mi_no=20;

y_completed_nm=y_miss;

y_imputed_mat_nm=matrix(NA, mis_no, mi_no);

for (i in 1:mi_no)
{
# normal model imputation;
set.seed(i);
y_imputed_nm=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);
y_completed_nm[miss_seq]=y_imputed_nm;
y_imputed_mat_nm[,i]=y_imputed_nm;
}

# plot the scatter plot of the imputed values from the 5th imputation;
# Fig. 5.10, second row left
plot(dbirwt[miss_seq], y_imputed_mat_nm[,20], pch=17, col="red", xlab="DBIRWT", ylab="DGESTAT", main="Normal Imputation");


# model 2
# PMM imputation;

mi_no=20;

y_completed_pmm=y_miss;

y_imputed_mat_pmm=matrix(NA, mis_no, mi_no);

for (i in 1:mi_no)
{
# proper model imputation;
set.seed(i);
y_imputed_pmm=mice.impute.pmm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);
y_completed_pmm[miss_seq]=y_imputed_pmm;
y_imputed_mat_pmm[,i]=y_imputed_pmm;
}

# plot the scatter plot of the imputed values from the 5th imputation;
# Fig. 5.10, second row right
plot(dbirwt[miss_seq], y_imputed_mat_pmm[,20], pch=17, col="red", xlab="DBIRWT", ylab="DGESTAT", main="PMM Imputation");

# model 3;
# normal model imputation including the quadratic term of dbirwt;
# standardize dbirwt;
dbirwt_std=(dbirwt-mean(dbirwt))/sd(dbirwt);
dbirwt_std_square=dbirwt_std^2;

x_new=cbind(dbirwt_std, dbirwt_std_square, data[,3:9]);

mi_no=20;

y_completed_nm_square=y_miss;

y_imputed_mat_nm_square=matrix(NA, mis_no, mi_no);

for (i in 1:mi_no)
{
# proper model imputation;
set.seed(i);
y_imputed_nm_square=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x_new);
y_completed_nm_square[miss_seq]=y_imputed_nm_square;
y_imputed_mat_nm_square[,i]=y_imputed_nm_square;
}

# plot the scatter plot of the imputed values from the 5th imputation;
# Fig. 5.10, third row left
plot(dbirwt[miss_seq], y_imputed_mat_nm_square[,20], pch=17, col="red", xlab="DBIRWT", ylab="DGESTAT", main="Normal Quadratic Imputation");

# model 4;
# PMM imputation including the quadratic term of dbirwt;

mi_no=20;

y_completed_pmm_square=y_miss;

y_imputed_mat_pmm_square=matrix(NA, mis_no, mi_no);

for (i in 1:mi_no)
{
# proper model imputation;
set.seed(i);
y_imputed_pmm_square=mice.impute.pmm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x_new);
y_completed_pmm_square[miss_seq]=y_imputed_pmm_square;
y_imputed_mat_pmm_square[,i]=y_imputed_pmm_square;
}

# plot the scatter plot of the imputed values from the 5th imputation;
# Fig. 5.10, third row right
plot(dbirwt[miss_seq], y_imputed_mat_pmm_square[,20], pch=17, col="red", xlab="DBIRWT", ylab="DGESTAT", main="PMM Quadratic Imputation");



# Histograms;

# Fig. 5.9
# Top row: obseved values
hist(y_miss, xlim=c(19, 48), xlab="DGESTAT", main="Observed", breaks=seq(19,48,1))

# Second row left: normal model imputation
hist(y_imputed_mat_nm[,20], xlim=c(19,48), xlab="DGESTAT", breaks=seq(19,48,1), main="Normal Imputation", add=F)

# Second row right: PMM imputation
hist(y_imputed_mat_pmm[,20], xlim=c(19,48), breaks=seq(19,48,1), xlab="DGESTAT", main="PMM Imputation", add=F)

# Third row left: normal model imputation including the quadratic term of dbirwt
hist(y_imputed_mat_nm_square[,20], xlim=c(19,48), breaks=seq(19,48,1), xlab="DGESTAT", main="Normal Quadratic Imputation", add=F)

# Third row right: PMM imputation including the quadratic term of dbirwt
hist(y_imputed_mat_pmm_square[,20], xlim=c(19,48), breaks=seq(19,48,1), xlab="DGESTAT", main="PMM Quadratic Imputation", add=F)


