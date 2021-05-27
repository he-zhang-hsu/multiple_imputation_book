# Example 5.2

rm(list=ls());



# read the file;

income=scan(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter1_missingdataproblem\\program\\income.dat");

# Fig. 5.2;
hist(income^(1/3), main="Histogram of Transformed Income", xlab="Transfromed Income");

qqnorm(income^(1/3), main="QQ Plot of Transformed Income");
qqline(income^(1/3), main="QQ Plot of Transformed Income");
