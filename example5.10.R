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

# Example 5.9;

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


# run the propensity score analysis;
dgestat_ps=glm(1-miss_indi ~ dbirwt+dmage+dfage+fmaps+wtgain+nprevis+cigar+drink,  family=binomial(link="logit"));

summary(dgestat_ps);

predicted=predict(dgestat_ps, x, type="response");
weight_ori=1/predicted;  
weight_ori[weight_ori >1.5]=1.5;
# weight=(weight_ori-mean(weight_ori))/sd(weight_ori);

weight=weight_ori;
summary(weight_ori);


# for plots in example 5.10
plot(predicted, dgestat, xlab="response propensity score");

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


# model I
# normal imputation;

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
plot(dbirwt[miss_seq], y_imputed_mat_nm[,20], pch=17, col="red", xlab="DBIRWT", ylab="DGESTAT", main="Normal Imputation");

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
# qudratic normal imputation;
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

plot(dbirwt[miss_seq], y_imputed_mat_nm_square[,20], pch=17, col="red", xlab="DBIRWT", ylab="DGESTAT", main="Normal Quadratic Imputation");

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
# linear imputation using weights;
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

plot(dbirwt[miss_seq], y_imputed_mat_weight[,20], pch=17, col="red", xlab="DBIRWT", ylab="DGESTAT", main="Propensity Weight Imputation");

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
x_new=cbind(x, weight*x);
# x_new=cbind(x, predicted*x);

mi_no=20;

y_completed_weight_int=y_miss;

y_imputed_mat_weight_int=matrix(NA, mis_no, mi_no);

for (i in 1:mi_no)
{
# proper model imputation;
set.seed(i*2);
y_imputed_weight_int=mice.impute.norm(y_miss, ry=as.logical(1-miss_indi), seed=i, x=x_new);
y_completed_weight_int[miss_seq]=y_imputed_weight_int;
y_imputed_mat_weight_int[,i]=y_imputed_weight_int;
}


# for plot in example 5.10.

plot(dbirwt[miss_seq], y_imputed_mat_weight_int[,20], pch=17, col="red", xlab="DBIRWT", ylab="DGESTAT", main="Propensity Weight Interaction");

dgestat_weight_int_quad=lm(y_completed_weight_int ~ dbirwt + I(dbirwt*dbirwt)+ dmage + dfage + fmaps + wtgain + nprevis + cigar + drink);

summary(dgestat_weight_int_quad);




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

