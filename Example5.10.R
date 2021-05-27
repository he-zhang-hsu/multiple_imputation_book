# Example 5.10

rm(list=ls());
# library(car);
# library(mvtnorm);
library(mice);
library(norm);
# library(HI);
options(digits=4);

# read the gestational age dataset;
data=read.table(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter4_univariate_model_parametric\\program\\dgestat.txt", col.names=c("dgestat", "dbirwt", "dmage", "dfage", "fmaps", "wtgain", "nprevis", "cigar", "drink"), na.strings=".");

# fetch the outcome;
y=data$dgestat;
miss_indi=is.na(y);
dgestat=y;

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



# descriptive statistics;
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

# fetch the predictors;
x=data[,2:9];


# run the propensity score analysis;
dgestat_ps=glm(1-miss_indi ~ dbirwt+dmage+dfage+fmaps+wtgain+nprevis+cigar+drink,  family=binomial(link="logit"));

summary(dgestat_ps);

predicted=predict(dgestat_ps, x, type="response");
weight_ori=1/predicted;

# top code the propensity weights;  
weight_ori[weight_ori >1.5]=1.5;
# weight=(weight_ori-mean(weight_ori))/sd(weight_ori);

weight=weight_ori;
summary(weight_ori);


# for plots in example 5.10
# Fig. 5.11 left panel;
plot(predicted, dgestat, xlab="response propensity score");

# Fig. 5.11 right panel;
plot(weight, dgestat, xlab="response weight");



# some exploratory analyses;

# regressing birthweight on other variables;
dgestat_cen=y-mean(y, na.rm=T);
dgestat=y;

dbirwt_pred=lm(dbirwt ~ dgestat_cen + I(dgestat_cen*dgestat_cen) + dmage + dfage + fmaps + wtgain + nprevis + cigar + drink);
summary(dbirwt_pred);

# regressing gestational age on other variables;

dgestat_pred_linear=lm(dgestat ~ dbirwt + dmage + dfage + fmaps + wtgain + nprevis + cigar + drink);
dgestat_pred_quad=lm(dgestat ~ dbirwt + I(dbirwt*dbirwt)+ dmage + dfage + fmaps + wtgain + nprevis + cigar + drink);

summary(dgestat_pred_quad);


# Different multiple imputation methods
# model I
# normal imputation; NMI

mi_no=20;

y_completed_nm=y_miss;

y_imputed_mat_nm=matrix(NA, mis_no, mi_no);

for (i in 1:mi_no)
{
# proper model imputation;
set.seed(i);
y_imputed_nm=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);
y_completed_nm[miss_seq]=y_imputed_nm;
y_imputed_mat_nm[,i]=y_imputed_nm;
}

# plot the histograms of the observed data and 5 imputed data;
# plot(dbirwt[miss_seq], y_imputed_mat_nm[,20], pch=17, col="red", xlab="DBIRWT", ylab="DGESTAT", main="Normal Imputation");

# R^2 statistics
dgestat_nm_quad=lm(y_completed_nm ~ dbirwt + I(dbirwt*dbirwt)+ dmage + dfage + fmaps + wtgain + nprevis + cigar + drink);

summary(dgestat_nm_quad);


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

plot(dbirwt[miss_seq], y_imputed_mat_pmm[,20], pch=17, col="red", xlab="DBIRWT", ylab="DGESTAT", main="PMM Imputation");

# model 3;
# qudratic normal imputation; NMI-QUAD
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

# plot(dbirwt[miss_seq], y_imputed_mat_nm_square[,20], pch=17, col="red", xlab="DBIRWT", ylab="DGESTAT", main="Normal Quadratic Imputation");

# R^2 statistics
dgestat_nm_square_quad=lm(y_completed_nm_square ~ dbirwt + I(dbirwt*dbirwt)+ dmage + dfage + fmaps + wtgain + nprevis + cigar + drink);
summary(dgestat_nm_square_quad);

# model 4;
# quadratic PMM imputation;

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

plot(dbirwt[miss_seq], y_imputed_mat_pmm_square[,20], pch=17, col="red", xlab="DBIRWT", ylab="DGESTAT", main="PMM Quadratic Imputation");

# model 5;
# PMI: linear imputation using weights;
mi_no=20;

y_completed_weight=y_miss;

y_imputed_mat_weight=matrix(NA, mis_no, mi_no);

for (i in 1:mi_no)
{
# proper model imputation;
set.seed(2*i);
y_imputed_weight=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x,weight));
y_completed_weight[miss_seq]=y_imputed_weight;
y_imputed_mat_weight[,i]=y_imputed_weight;
}

# for plot in example 5.10.
# Fig. 5.12 left panel
plot(dbirwt[miss_seq], y_imputed_mat_weight[,20], pch=17, col="red", xlab="DBIRWT", ylab="DGESTAT", main="Propensity Weight Imputation");

# R^2 statistics
dgestat_weight_quad=lm(y_completed_weight ~ dbirwt + I(dbirwt*dbirwt)+ dmage + dfage + fmaps + wtgain + nprevis + cigar + drink);
summary(dgestat_weight_quad);

# model 6 pmm with covariates and weight;

mi_no=20;

y_completed_weight_pmm=y_miss;

y_imputed_mat_weight_pmm=matrix(NA, mis_no, mi_no);

for (i in 1:mi_no)
{
# proper model imputation;
set.seed(i);
y_imputed_weight_pmm=mice.impute.pmm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=cbind(x,weight));
y_completed_weight_pmm[miss_seq]=y_imputed_weight_pmm;
y_imputed_mat_weight_pmm[,i]=y_imputed_weight_pmm;
}


plot(dbirwt[miss_seq], y_imputed_mat_weight_pmm[,20], pch=17, col="red", xlab="DBIRWT", ylab="DGESTAT", main="Propensity Weight PMM");

dgestat_weight_pmm_quad=lm(y_completed_weight_pmm ~ dbirwt + I(dbirwt*dbirwt)+ dmage + dfage + fmaps + wtgain + nprevis + cigar + drink);

summary(dgestat_weight_pmm_quad);



# model 7 with weight interact with other covariates;
# PMI-INT
x_new=cbind(x, weight*x);
# x_new=cbind(x, predicted*x);

mi_no=20;

y_completed_weight_int=y_miss;

y_imputed_mat_weight_int=matrix(NA, mis_no, mi_no);

for (i in 1:mi_no)
{
# proper model imputation;
set.seed(2*i);
y_imputed_weight_int=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x_new);
y_completed_weight_int[miss_seq]=y_imputed_weight_int;
y_imputed_mat_weight_int[,i]=y_imputed_weight_int;
}


# for plot in example 5.10.
# Fig. 5.12 right panel
plot(dbirwt[miss_seq], y_imputed_mat_weight_int[,20], pch=17, col="red", xlab="DBIRWT", ylab="DGESTAT", main="Propensity Weight Interaction");

# R^2 statistics
dgestat_weight_int_quad=lm(y_completed_weight_int ~ dbirwt + I(dbirwt*dbirwt)+ dmage + dfage + fmaps + wtgain + nprevis + cigar + drink);

summary(dgestat_weight_int_quad);



