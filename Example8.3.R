# Example 8.3

rm(list=ls());
library(R2WinBUGS);
library(MASS);
library(mice);
library(norm);
library(mix);
# library(HI);
options(digits=4);

# WinBUGS syntax;

model0.file <- system.file(package="R2WinBUGS", "model", "weibull_imputeforR.txt")
# Let's take a look:
file.show(model0.file)


# read the data;


# create the working data;
M=4;
N=20;


t = matrix(c(12, 1, 21, 25, 11, 26, 27, 30, 13, 12, 21, 20, 23, 25, 23, 29, 35, NA, 31, 36,  
32, 27, 23, 12, 18, NA, NA, 38, 29, 30, NA, 32, NA, NA, NA, NA, 25, 30, 37, 27,  
22, 26, NA, 28, 19, 15, 12, 35, 35, 10, 22, 18, NA, 12, NA, NA, 31, 24, 37, 29,  
27, 18, 22, 13, 18, 29, 28, NA, 16, 22, 26, 19, NA, NA, 17, 28, 26, 12, 17, 26), M, N, byrow=T);

t.cen =matrix(c( 0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 40,  0,  0,  
0,   0,  0,  0,  0, 40, 40,  0,  0,  0, 40,  0, 40, 40, 40, 40,  0,  0,  0,  0,  
			  0,  0, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0, 24,  0, 40, 40,  0,  0,  0,  0,  
0,  0,  0,  0,  0,  0,  0, 20,  0,  0,  0,  0, 29, 10,  0,  0,  0,  0,  0,  0), M, N, byrow=T);

process_stage2_data_simulate=list("M","N","t", "t.cen");


# create the initial values;

process_stage2_init=function(){
list(beta=rep(0,4),r=1)
};

process_stage2_parameters=c("t_impute");

gibbs_no=10000;

process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file=model0.file,
    n.chains=1, n.thin=100, n.iter=gibbs_no, n.burnin=gibbs_no/2, bugs.seed=197552, bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);

names(process_stage2.sim);

attach.bugs(process_stage2.sim);

# the 1st dimention is the number of iterations (50);

# the 1st 15 cases of the 2nd and 3rd dimension correspond to 
# the censored cases;
# Table 8.4

# apply(t_impute, c(2,3), mean);
# t_impute[,1,1]-t_impute[,4,1];
summary(t_impute[,1,1]);
summary(t_impute[,1,2]);
summary(t_impute[,1,3]);
summary(t_impute[,1,4]);
summary(t_impute[,1,5]);
summary(t_impute[,1,6]);
summary(t_impute[,1,7]);
summary(t_impute[,1,8]);
summary(t_impute[,1,9]);
summary(t_impute[,1,10]);
summary(t_impute[,2,1]);
summary(t_impute[,2,2]);
summary(t_impute[,2,3]);
summary(t_impute[,2,4]);
summary(t_impute[,2,5]);





