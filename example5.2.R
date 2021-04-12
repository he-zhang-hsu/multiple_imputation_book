rm(list=ls());



# income_weight=read.table(file="//cdc.gov/private/L728/wdq7/NHIS_imputation/data_for_book/rawincome.dat", na.strings=".");

income_weight=read.table(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter1_missingdataproblem\\program\\rawincome.dat", na.strings=".");

income=income_weight[,1];
weight=income_weight[,2];

# the full sample size of NHIS 2016;
rowobs=length(income);

# the number of missing cases;

mis_no=sum(is.na(income));

# the missingess frequency;
mis_no/rowobs;


lambda_com=powerTransform(income+0.001);

lambda_com;


quantile(income, 0.95, na.rm=T);




summary(income);

hist(income^(1/3));
hist(income);


# top code the income value;


income[income >206000]=206000;

hist(income, main="Histogram of Observed Family Total Income", xlab="Family Total Income in $", freq=TRUE);

# saved as example1.1.hist.ps;

qqnorm(income, main="Normal Q-Q Plot of Observed Family Total Income");
qqline(income);


# saved as example1.1.qq.ps;


# plot the transformed scale income;


income[income >206000]=206000;

# save as example 5.2 plot;
hist(income^(1/3), main="Histogram of Transformed Income", xlab="Transfromed Income");

qqnorm(income^(1/3), main="QQ Plot of Transformed Income");
qqline(income^(1/3), main="QQ Plot of Transformed Income");


# read imputation of income;
# this part should go to imputation diagnostics;

income_imputed=read.table(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter5_univariate_model_robust\\program\\income2016_1.txt", na.strings=".");

income=income_imputed[,1];
missingindicator=income_imputed[,2];

# make two plots;
# the imputation tends to impute lower income people;
par(mfrow=c(1,2));
hist(income[missingindicator==1]^(1/3), freq=FALSE);
hist(income[missingindicator==0]^(1/3), freq=FALSE);

qqnorm(income[missingindicator==1]^(1/3));
qqline(income[missingindicator==1]^(1/3));
qqnorm(income[missingindicator==0]^(1/3));
qqline(income[missingindicator==0]^(1/3));

# compare the distribution of the observed and imputed data;
t.test(income[missingindicator==1], income[missingindicator==0]);





hist(weight, main="Distribution of Family-level Survey Weight", freq=TRUE);
summary(weight);
weighted.mean(income, w=weight, na.rm=T);


