################## data generate ##############
# Program: boxcox_impute_expolore_simulation; ############
# generate data, and try out different missing data methods;
# run some simulations to assess the biasedness of various methods;

# library(car);
# library(mvtnorm);
library(mice);
library(norm);
# library(HI);
options(digits=4);
rm(list=ls());

# Example 12.10;

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






# make some plots;

# fetch the predictors;
x=data[,2:9];
x_nodbirwt=data[,3:9];


# run the propensity score analysis;
dgestat_ps=glm(miss_indi ~ dbirwt+dmage+dfage+fmaps+wtgain+nprevis+cigar+drink,  family=binomial(link="logit"));

summary(dgestat_ps);

# some exploratory analyses;

dgestat_cen=y-mean(y, na.rm=T);

dbirwt_pred=lm(dbirwt ~ dgestat_cen + I(dgestat_cen*dgestat_cen) + dmage + dfage + fmaps + wtgain + nprevis + cigar + drink);
summary(dbirwt_pred);

# model 0
# normal imputation without DBIRWT;

mi_no=100;

# matrix holding mean and variance estimates of multiple imputation;
MI_mean_mat=MI_mean_var_mat=rep(NA, mi_no);



y_completed_nm=y_miss;

for (i in 1:mi_no)
{
# proper model imputation;
set.seed(i);
y_imputed_nm=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x_nodbirwt);
y_completed_nm[miss_seq]=y_imputed_nm;
MI_mean_mat[i]=mean(y_completed_nm);
MI_mean_var_mat[i]=var(y_completed_nm)/rowobs;
}

mean_summary=pool.scalar(MI_mean_mat, MI_mean_var_mat, n=rowobs);
mean_summary$qbar;
sqrt(mean_summary$t);
mean_summary$fmi;

#mean_mi_x1_mean[i]=mean_summary$qbar;
#mean_mi_x1_var[i]=mean_summary$t;
#mean_mi_x1_df[i]=mean_summary$df;
#mean_mi_x1_f[i]=mean_summary$fmi;




# model I
# normal imputation;

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

mean_summary=pool.scalar(MI_mean_mat, MI_mean_var_mat, n=rowobs);
mean_summary$qbar;
sqrt(mean_summary$t);
mean_summary$fmi;



# model 2
# PMM imputation;

# mi_no=20;

# y_completed_pmm=y_miss;

# y_imputed_mat_pmm=matrix(NA, mis_no, mi_no);

# for (i in 1:mi_no)
# {
# proper model imputation;
# set.seed(i);
# y_imputed_pmm=mice.impute.pmm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x);
# y_completed_pmm[miss_seq]=y_imputed_pmm;
# y_imputed_mat_pmm[,i]=y_imputed_pmm;
# }

# plot(dbirwt[miss_seq], y_imputed_mat_pmm[,20], pch=17, col="red", xlab="DBIRWT", ylab="DGESTAT", main="PMM Imputation");

# model 3;
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
# proper model imputation;
set.seed(i);
y_imputed_nm_square=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x_new);
y_completed_nm_square[miss_seq]=y_imputed_nm_square;
MI_mean_mat[i]=mean(y_completed_nm_square);
MI_mean_var_mat[i]=var(y_completed_nm_square)/rowobs;

}

mean_summary=pool.scalar(MI_mean_mat, MI_mean_var_mat, n=rowobs);
mean_summary$qbar;
sqrt(mean_summary$t);
mean_summary$fmi;



# model 4;
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
mean_summary$qbar;
sqrt(mean_summary$t);
mean_summary$fmi;



# red for observed cases and blue for imputed cases;
# par(mfrow=c(2,1));

par(mfrow=c(1,1));


hist(y_miss, xlim=c(19, 48), xlab="DGESTAT", main="Observed", breaks=seq(19,48,1))


hist(y_imputed_mat_nm[,20], xlim=c(19,48), xlab="DGESTAT", breaks=seq(19,48,1), main="Normal Imputation", add=F)


hist(y_imputed_mat_pmm[,20], xlim=c(19,48), breaks=seq(19,48,1), xlab="DGESTAT", main="PMM Imputation", add=F)


hist(y_imputed_mat_nm_square[,20], xlim=c(19,48), breaks=seq(19,48,1), xlab="DGESTAT", main="Normal Quadratic Imputation", add=F)

hist(y_imputed_mat_pmm_square[,20], xlim=c(19,48), breaks=seq(19,48,1), xlab="DGESTAT", main="PMM Quadratic Imputation", add=F)

plot(ecdf(y_miss));

plot(ecdf(y_imputed_mat_nm[,20]));

plot(ecdf(y_imputed_mat_pmm[,20]));

plot(ecdf(y_imputed_mat_nm_square[,20]));

plot(ecdf(y_imputed_mat_pmm_square[,20]));



# box()

summary(as.vector(y_imputed_mat));


# h2<-rnorm(1000,4)
# h1<-rnorm(1000,6)

# Histogram Grey Color
# hist(h1, col=rgb(0.1,0.1,0.1,0.5),xlim=c(0,10), ylim=c(0,200))
# hist(h2, , add=T)
# box()

# Histogram Colored (blue and red)
# hist(h1, col=rgb(1,0,0,0.5),xlim=c(0,10), ylim=c(0,200), main=???¡ì?¨¨Overlapping Histogram???¡ì?¨¨, xlab=???¡ì?¨¨Variable???¡ì?¨¨)
# hist(h2, col=rgb(0,0,1,0.5), add=T)
# box()



#Random numbers
# h2<-rnorm(1000,4)
# h1<-rnorm(800,6)

# combined on density scales;


# Histogram Grey Color
# hist(h1,col=rgb(0.1,0.1,0.1,0.5), freq=FALSE)
# hist(h2,col=rgb(0.8,0.8,0.8,0.5), freq=FALSE, add=T)
# box()

# Histogram Colored (blue and red)

# hist(h1, col=rgb(1,0,0,0.5), freq=FALSE)
# hist(h2, col=rgb(0,0,1,0.5), freq=FALSE, add=T)
# box()

