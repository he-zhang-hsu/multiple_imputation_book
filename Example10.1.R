# Example 10.1


# Read the data;

# The final survey weights from RANDS 1 data;
rands_weight=read.table(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter10_survey\\program\\rands1_finalweight.txt", na.strings=".");

# Fig. 10.1

hist(rands_weight$V1, main="", xlab="Survey Weights of RANDS 1");
