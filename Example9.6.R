# Example 9.6

rm(list=ls());
library(R2WinBUGS);
library(MASS);
library(mice);
library(norm);
library(mix);
# library(HI);
options(digits=4);

# WinBUGS syntax;

model0.file <- system.file(package="R2WinBUGS", "model", "linearmixed_covmissing_imputeforR.txt")
# Let's take a look:
file.show(model0.file)


# read the data;

oridata=scan(file=paste("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter9_longitudinal\\program\\psid_uni_uni_hr_2172_log.dat"), na.strings=".");

orimatrix=matrix(oridata, nrow=2172, ncol=34, byrow=TRUE);

# oriy concatenate the response vector;
oriy=orimatrix[,13:31];

oriyvec=as.vector(t(oriy));

summary(oriyvec);
   

# fetch the baseline covariates;

# assign some covariate observations as missing;
health=orimatrix[,32];
edu=orimatrix[,33];
race=orimatrix[,34];

mean(race);
mean(health);
mean(edu);
# glm(edu ~ race);
# glm(health ~ edu + race + edu*race);

# sample size;

rowobs=nrow(orimatrix);

# creat 10% missing data for health, edu, and race;
set.seed(197789);

miss_indi=matrix(rbinom(n=3*rowobs, size=1, prob=0.1), nrow=rowobs, ncol=3);
# set some of the ori_data as missing values;

health[miss_indi[,1]==1]=NA;
edu[miss_indi[,2]==1]=NA;
race[miss_indi[,3]==1]=NA;

# create the working data;

N=2172;
T=19;
x = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19);
xbar=10;
Y=oriy;

process_stage2_data_simulate=list("Y","x","xbar", "health", "edu", "race", "T","N");

# random intercept and slope;

# create the initial values;

process_stage2_init=function(){
list(beta.c=0.037, alpha.c=1.463, alpha.health=0.093, alpha.edu=0.38, alpha.race=0.550, slope=0.023, sigma.alpha=0.448, 
sigma.beta.alpha=0.038, sigma.error=0.368, mu_race=0.065, alpha.edu.race=-0.739, beta.edu.race=0.87, 
alpha.health.edu.race=1.05, beta.health.edu=0.36, beta.health.race=0.81, beta.health.edu.race=0.137,
alpha=rep(0,rowobs), beta=rep(0,rowobs))
};

process_stage2_parameters=c("beta.c", "alpha.c", "alpha.health", "alpha.edu", "alpha.race", "slope", "sigma.alpha", 
"sigma.beta.alpha", "sigma.error", "mu_race", "alpha.edu.race", "beta.edu.race", "alpha.health.edu.race", 
"beta.health.edu", "beta.health.race", "beta.health.edu.race");

gibbs_no=20000;

process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file=model0.file,
    n.chains=1, n.thin=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, bugs.seed=197552, bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);

names(process_stage2.sim);

# dic value;

# process_stage2.sim$pD;
# process_stage2.sim$DIC;

attach.bugs(process_stage2.sim);

# posterior means;
# Table 9.5

mean(beta.c);
sd(beta.c);
quantile(beta.c, c(.025, .975));

mean(alpha.c);
sd(alpha.c);
quantile(alpha.c, c(.025, .975));

mean(alpha.health);
sd(alpha.health);
quantile(alpha.health, c(.025, .975));

mean(alpha.edu);
sd(alpha.edu);
quantile(alpha.edu, c(.025, .975));

mean(alpha.race);
sd(alpha.race);
quantile(alpha.race, c(.025, .975));

mean(slope);
sd(slope);
quantile(slope, c(.025, .975));

mean(sigma.alpha);
sd(sigma.alpha);
quantile(sigma.alpha, c(.025, .975));

mean(sigma.beta.alpha);
sd(sigma.beta.alpha);
quantile(sigma.beta.alpha, c(.025, .975));

mean(sigma.error);
sd(sigma.error);
quantile(sigma.error, c(.025, .975));

mean(mu_race);
sd(mu_race);
quantile(mu_race, c(.025, .975));

mean(alpha.edu.race);
sd(alpha.edu.race);
quantile(alpha.edu.race, c(.025, .975));

mean(beta.edu.race);
sd(beta.edu.race);
quantile(beta.edu.race, c(.025, .975));

mean(alpha.health.edu.race);
sd(alpha.health.edu.race);
quantile(alpha.health.edu.race, c(.025, .975));

mean(beta.health.edu);
sd(beta.health.edu);
quantile(beta.health.edu, c(.025, .975));

mean(beta.health.race);
sd(beta.health.race);
quantile(beta.health.race, c(.025, .975));

mean(beta.health.edu.race);
sd(beta.health.edu.race);
quantile(beta.health.edu.race, c(.025, .975));


