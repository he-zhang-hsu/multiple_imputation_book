# Example 11.11

# Generate four binary variables
# r1-r3 and col
# These variables are going to be used in a SAS file

library(bindata);


cor.mat<-matrix (c(1, -0.25, -0.15, 0.1, -0.25, 1.00, 0.35, 0.2,  -0.15, 0.35, 1.00, 0.3, 0.1, 0.2, 0.3,1), nrow=4, ncol=4, byrow=TRUE)

data<-rmvbin(1000000, margprob=c(0.8, 0.3, 0.2, 0.3), bincorr=cor.mat)

summary(data)

cor(data)

write.csv(data,file="\\\\cdc.gov\\private\\M132\\jnu5\\QSDMB\\Incidence Model\\Multiple Imputation\\Rick_tran_cate\\program\\Simulation\\Data Generation\\data.csv")
