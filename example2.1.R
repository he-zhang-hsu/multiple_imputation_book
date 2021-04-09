rm(list=ls());
library(boot);
library(Hmisc);


# sample.size=scan(file="//cdc.gov/private/L728/wdq7/old_projects/hospital_classification/data/sample_size.dat", na.strings=".");
# bmxarmc=scan(file="C:\\Users\\yhsc\\Personal\\Yulei\\research\\MIbook\\chapter2_statisticalbackground\\examples\\bmxarmc.txt", na.strings=".");

# data=read.table(file="C:\\Users\\yhsc\\Personal\\Yulei\\research\\MIbook\\chapter2_statisticalbackground\\examples\\bmxarmc.txt", col.names=c("weight", "bmx"), na.strings=".");

# Numbers of 1's and 0's for the RAND I health insurance coverage example

data=rep(0, 2304);
data[1:2152]=1;


meanfunction=function(x, indices) {mean(x[indices])}

percentilefunction5=function(x, indices) {quantile(x[indices], 0.05)}

percentilefunction95=function(x, indices) {quantile(x[indices], 0.95)}

set.seed(197789);

# bootstrap standard errors;
initial_result=boot(data=data, statistic=meanfunction, R=10000);
CI=boot.ci(initial_result);

sqrt(var(initial_result$t));


# plot the bootstrap estimate distribution
hist(initial_result$t, xlab="Proportion estimate", main="Bootstrap samples");

qqnorm(initial_result$t, main="Normal QQ plot");
qqline(initial_result$t);



# simulate the posterior estimate for Example 2.2;

Y=2152;
n=2304;
a=Y+1;
b=n-Y+1;

set.seed(197789);

x=rbeta(10000, a, b);

mean(x);
var(x);
quantile(x,.025);
quantile(x,.975);


hist(x, main="Posterior samples", xlab="Theta");
qqnorm(x, main="Normal QQ plot");
qqline(x);

