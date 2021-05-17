rm(list=ls());



# income_weight=read.table(file="//cdc.gov/private/L728/wdq7/NHIS_imputation/data_for_book/rawincome.dat", na.strings=".");

# income_weight=read.table(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter1_missingdataproblem\\program\\rawincome.dat", na.strings=".");


#  income=income_weight[,1];
# weight=income_weight[,2];


# top coding
# income[income >206000]=206000;

# income_file=write.table(income, file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter1_missingdataproblem\\program\\income.dat", 
# row.names=FALSE, col.names=FALSE, na="NA");

# read the file;

income=scan(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter1_missingdataproblem\\program\\income.dat");


# the full sample size;
rowobs=length(income);

# the number of missing cases;

mis_no=sum(is.na(income));

# the missingess frequency;
mis_no/rowobs;


# lambda_com=powerTransform(income+0.001);

# lambda_com;

# summary(income);

# hist(income^(1/3));
# hist(income);


# The following plots are for Fig. 1.1 in the book

hist(income, main="Histogram of observed income values", xlab="Total family income in $", freq=TRUE);

# saved as example1.1.hist.ps;

qqnorm(income, main="QQ plot of observed income values");
qqline(income);

# saved as example1.1.qq.ps;



