rm(list=ls());

library(R2WinBUGS);


# setwd('\\\\cdc.gov\\private\\L728\\book\\example8.4')

# An example model file is given in:
model.file <- system.file(package="R2WinBUGS", "model", "schools.txt")
# Let's take a look:
file.show(model.file)

# Some example data (see ?schools for details):
data(schools)
schools

J <- nrow(schools)
y <- schools$estimate
sigma.y <- schools$sd
data <- list(J=J, y=y, sigma.y=sigma.y)
inits <- function(){
    list(theta=rnorm(J, 0, 100), mu.theta=rnorm(1, 0, 100),
         sigma.theta=runif(1, 0, 100))
}
## or alternatively something like:
# inits <- list(
#   list(theta=rnorm(J, 0, 90), mu.theta=rnorm(1, 0, 90),
#        sigma.theta=runif(1, 0, 90)),
#   list(theta=rnorm(J, 0, 100), mu.theta=rnorm(1, 0, 100),
#        sigma.theta=runif(1, 0, 100))
#   list(theta=rnorm(J, 0, 110), mu.theta=rnorm(1, 0, 110),
#        sigma.theta=runif(1, 0, 110)))

parameters <- c("theta", "mu.theta", "sigma.theta")

## Not run: 
## You may need to edit "bugs.directory",
## also you need write access in the working directory:
schools.sim <- bugs(data, inits, parameters, model.file,
    n.chains=3, n.iter=5000,
  bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", debug=TRUE);

#    bugs.directory="c:/Program Files/WinBUGS14/")
print(schools.sim)
plot(schools.sim)

## End(Not run)


# generate data for example 8.4


model.file <- system.file(package="R2WinBUGS", "model", "logitmodel_missing_forR.txt")
# Let's take a look:
file.show(model.file)

set.seed(197789);
x=rnorm(1000);
beta0=-1;
beta1=1;
y2=rbinom(1000,size=1,prob=exp(beta0+beta1*x)/(1+exp(beta0+beta1*x)));
data=cbind(y2,x);

# apply missing cases to the data;

missing_indi=cbind(rbinom(1000,size=1,prob=0.25), rbinom(1000,size=1, prob=0.25));

y2[missing_indi[,1]==1]=NA;
x[missing_indi[,2]==1]=NA;

missdata=cbind(y2,x);

M=1000;






# process_stage2_data_simulate=list("M", "y2", "x");
process_stage2_data_simulate=list(M=M, y2=y2, x=x);

process_stage2_init= function(){list(beta0=-1, beta1=1, mu=0, tau=1)};
process_stage2_parameters=c("mu", "beta0", "beta1", "y2", "x");

gibbs_no=20000;


process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file,
    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=1, bugs.seed=197789,
    bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", debug=TRUE);



#process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file="C:/Users/WDQ7/reliability_test_forR.txt",
#    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=1, bugs.seed=197789,
#    bugs.directory="C:/Program Files (x86)/OpenBUGS/OpenBUGS321", summary.only=FALSE, debug=FALSE,
#    working.directory=NULL, clearWD=TRUE);


attach.bugs(process_stage2.sim);
rm(process_stage2.sim);

nrow(y2);
nrow(x);
ncol(y2);
ncol(x);

median(mu0);
median(tau2);
median(sigma2);
median(sqrt(tau2/sigma2));
median(K);

# output the data;


hist(process_data$r2/process_data$r1, breaks=seq(0, 1.0, 0.02), xlab="Aspirin raw rate", main="Disribution of Hospital-level Aspirin Raw Rates");
abline(v=0.90);

process_data=process_data[process_data$r1 < 30,];
# process_data=process_data[process_data$r1 >= 30,];
mu_actual=mean(process_data$r2/process_data$r1);
tau_actual=sd(process_data$r2/process_data$r1);
tau_square_actual=tau_actual^2;
sigma_square_actual=mu_actual-mu_actual^2-tau_actual^2;
K_actual=sqrt(tau_square_actual)/sqrt(sigma_square_actual);

(1-0.96^2)/0.96^2*K^2;

# number of hospitals;
nrow(process_data);
summary(process_data$r1);
hist(process_data$r1, xlab="sample size", main="");
hist(1/process_data$r1);
mean(1/process_data$r1);

# K is the ratio between between-provider SD and within-provider SD 
# K_vec=seq(0.2,5.0,0.05);
K_vec=seq(0.2, 2, 0.02);
length_K=length(K_vec);
rho_mle_vec=rep(0, length_K);
rho_eb_vec=rep(0, length_K);
# rho_prob1_vec fixes the probability treshold first;
rho_prob1_vec=rep(0, length_K);

# rho_prob2_vec fixes the cut-off point first;
rho_prob2_vec=rep(0, length_K);

for (i in 1:length_K)
{
K=K_vec[i];
rho_mle=K/sqrt(K^2+mean(1/process_data$r1));
rho_eb=sum(K/(K^2+1/process_data$r1))/sqrt(nrow(process_data)*sum(1/(K^2+1/process_data$r1)));

# prob is the probability threshold for the prob1 method;
prob=0.9;

# cut is the cut-off threshold for the prob2 method;
cut=0.9;

numerator=sum(K^3/(K^2+1/process_data$r1));
denominator1=nrow(process_data)*sum(K^4/(K^2+1/process_data$r1));
denominator2=qnorm(prob)^2*(nrow(process_data)*sum(K^2/(process_data$r1*K^2+1))-sum(sqrt(K^2/(process_data$r1*K^2+1)))^2);
rho_prob1=numerator/sqrt(denominator1+denominator2);

numerator2=sum(sqrt(process_data$r1*K^2/(K^2+1/process_data$r1)));
denominator2_1=nrow(process_data)*sum(process_data$r1);
denominator2_2=qnorm(cut)^2*(nrow(process_data)*sum(process_data$r1+1/K^2)-sum(sqrt(process_data$r1+1/K^2))^2);
rho_prob2=numerator2/sqrt(denominator2_1+denominator2_2);



rho_mle_vec[i]=rho_mle;
rho_eb_vec[i]=rho_eb;
rho_prob1_vec[i]=rho_prob1;
rho_prob2_vec[i]=rho_prob2;

}

rho_mle_vec_all=rho_mle_vec;
rho_eb_vec_all=rho_eb_vec;
rho_prob1_vec_all=rho_prob1_vec;
rho_prob2_vec_all=rho_prob2_vec;

rho_mle_vec_lt30=rho_mle_vec;
rho_eb_vec_lt30=rho_eb_vec;
rho_prob_vec_lt30=rho_prob_vec;

rho_mle_vec_ge30=rho_mle_vec;
rho_eb_vec_ge30=rho_eb_vec;
rho_prob_vec_ge30=rho_prob_vec;


matplot(K_vec, cbind(rho_mle_vec_all^2, rho_eb_vec_all^2), xlab="K=SD(BW)/SD(WI)", ylab="RELIABILITY", main="1:DIR, 2:SHR");



matplot(K_vec, cbind(rho_mle_vec_all^2, rho_eb_vec_all^2, rho_prob1_vec_all^2), xlab="K=SD(BW)/SD(WI)", ylab="RELIABILITY", main="1:DIR, 2:SHR, 3:PROB for all sample");

matplot(K_vec, cbind(rho_mle_vec_all, rho_eb_vec_all, rho_prob_vec_all, rho_mle_vec_lt30, rho_eb_vec_lt30, rho_prob_vec_lt30, rho_mle_vec_ge30, rho_eb_vec_ge30, rho_prob_vec_ge30), xlab="K=SD(BW)/SD(WI)", ylab="CORR", main="1/2/3: MLE/BMEAN/BPROB (all), 4/5/6: MLE/BMEAN/BPROB (<30), 7/8/9 MLE/BMEAN/BROB (>=30)");


matplot(K_vec, cbind(rho_mle_vec_lt30, rho_eb_vec_lt30, rho_prob_vec_lt30), xlab="SD-BW/SD_WI", ylab="CORR", main="1:MLE, 2:BMEAN, 3:BPROB for sample size lt 30");
matplot(K_vec, cbind(rho_mle_vec_ge30, rho_eb_vec_ge30, rho_prob_vec_ge30), xlab="SD-BW/SD_WI", ylab="CORR", main="1:MLE, 2:BMEAN, 3:BPROB for sample size ge 30", ylim=c(0.5, 1.0));


# calculate sensitivies of three approaches;
top_prob=0.9;
sens_mle_vec=rep(0, length_K);
sens_eb_vec=rep(0, length_K);
sens_prob_vec=rep(0, length_K);

for (i in 1:length_K)
{
corr_mle=matrix(c(1,rho_mle_vec[i],rho_mle_vec[i],1),2,2);
corr_eb=matrix(c(1,rho_eb_vec[i],rho_eb_vec[i],1),2,2);
corr_prob=matrix(c(1,rho_prob_vec[i],rho_prob_vec[i],1),2,2);

sens_mle_vec[i]=pmvnorm(lower=rep(qnorm(top_prob),2), upper=Inf, corr=corr_mle)/(1-top_prob);
sens_eb_vec[i]=pmvnorm(lower=rep(qnorm(top_prob),2), upper=Inf, corr=corr_eb)/(1-top_prob);
sens_prob_vec[i]=pmvnorm(lower=rep(qnorm(top_prob),2), upper=Inf, corr=corr_prob)/(1-top_prob);

}


matplot(K_vec, cbind(sens_mle_vec, sens_eb_vec, sens_prob_vec), xlab="SD-BW/SD_WI", ylab="Sensitivity", main="1:MLE, 2:BMEAN, 3:BPROB for classifying top 10%");


matplot(K_vec, cbind(sens_mle_vec, sens_eb_vec, sens_prob_vec), xlab="SD-BW/SD_WI", ylab="Sensitivity", main="1:MLE, 2:BMEAN, 3:BPROB for classifying top 10% and sample size gt 30");


# calculate the sample size cut-off point for MLE;
C_vec=seq(0.80,0.99, 0.01);
length_C=length(C_vec);
n_min_mat=matrix(0, nrow=length_C, ncol=length_K);

for (i in 1:length_C)
{
 for (j in 1:length_K)
{
n_min=(sqrt(1+C_vec[i]^2/(1-C_vec[i]^2))-1)/(2*K_vec[j]^2);
n_min_mat[i,j]=n_min;
}
}

# calculat the relationship between correlation coefficient and sensitivity
# calculate sensitivies of three approaches;
top_prob=0.5;
# top_prob=0.9;
rho_vec=seq(0, 0.99, 0.01);
length_rho=length(rho_vec);

sens_vec=rep(0, length_rho);

for (i in 1:length_rho)
{
corr_mat=matrix(c(1,rho_vec[i],rho_vec[i],1),2,2);
sens_vec[i]=pmvnorm(lower=rep(qnorm(top_prob),2), upper=Inf, corr=corr_mat)/(1-top_prob);
}

plot(rho_vec, sens_vec, xlab="correlation coefficient", ylab="sensitivity", main="sentivity for classifying top 10%");

# calculat the average of the inverse of the sample size;
hospital_size=sort(process_data$r1);
hospital_size_mean_vec=rep(0, length(hospital_size));
for (i in 1:length(hospital_size))
{
 hospital_size_part=hospital_size[i:length(hospital_size)]
 hospital_size_mean_vec[i]=mean(1/hospital_size_part);
}

# keep the cases whose numerator is greater than 30;
process_data=process_data[process_data$r1 >=30,];
nrow(process_data);
summary(process_data$r1);
hist(process_data$r1);
mean(process_data$r2/process_data$r1);
sd(process_data$r2/process_data$r1);
hist(process_data$r2/process_data$r1, xlab="performance rate", main="CAP5: Initial antibiotics received");
abline(v=0.85);

# for (i in 1:nrow(process_data))
# {
# if (process_data$r2[i]==0) process_data$r2[i]=0.5;
# if (process_data$r1[i]==process_data$r2[i]) 
# process_data$r1[i]=process_data$r1[i]+0.5;
# }


# simulate some data;
sample_size=hospital_size[1:1000];
K=0.5;
tau_square=1;
sigma_square=tau_square/K^2;
u_true=rnorm(1000, sd=sqrt(tau_square));
y_bar=u_true+rnorm(1000, sd=sqrt(sigma_square/sample_size));
cor(u_true, y_bar);
B=tau_square/(tau_square+sigma_square/sample_size);
y_eb=B*y_bar+(1-B)*(mean(u_true));
cor(u_true, y_eb);

rho_mle=K/sqrt(K^2+mean(1/sample_size));
rho_eb=sum(K/(K^2+1/sample_size))/sqrt(length(sample_size)*sum(1/(K^2+1/sample_size)));

# output the dataset;
write.table(cbind(sample_size, y_bar), file ="H:/cdrive_copy/eligibility_profiling_analysis/research_notes/samplesize/test.dat", row.name=FALSE, col.name=FALSE, sep=" ");

# calculate the rho for each hospital;
q=0.9;
prob_class=0.5;

process_data=process_data[process_data$r1 < 30,];

sample_size=process_data$r1;
N=length(sample_size);

mu=mean(process_data$r2/process_data$r1);
sigma=sd(process_data$r2/process_data$r1);
sigma_square=sigma^2;
tau_square=mu-mu^2-sigma^2;
K=sqrt(tau_square)/sigma;


true_mean=mu;
true_rho_vec=sqrt(K^2/(K^2+1/sample_size));



marginal_var_raw=tau_square+sigma_square/sample_size;
marginal_var_mean=tau_square^2/marginal_var_raw;
marginal_var_prob=marginal_var_mean;

raw_mixture=norMix(rep(true_mean, N), sig2=marginal_var_raw);
quantile_raw_mixture=qnorMix(p=q, obj=raw_mixture);

mean_mixture=norMix(rep(true_mean, N), sig2=marginal_var_mean);
quantile_mean_mixture=qnorMix(p=q, obj=mean_mixture);

prob_mixture=norMix(true_mean-qnorm(prob_class)*sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size), sig2=marginal_var_prob);
quantile_prob_mixture=qnorMix(p=q, obj=prob_mixture);

cut_off_raw=(quantile_raw_mixture-true_mean)/sqrt(marginal_var_raw);
cut_off_mean=(quantile_mean_mixture-true_mean)/sqrt(marginal_var_mean);
cut_off_prob=(quantile_prob_mixture-true_mean+qnorm(prob_class)*sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size))/sqrt(marginal_var_prob);

raw_mat=cbind(cut_off_raw, true_rho_vec);
mean_mat=cbind(cut_off_mean, true_rho_vec);
prob_mat=cbind(cut_off_prob, true_rho_vec);


sens_mixture_raw_vec=apply(raw_mat, 1, sens_formula_vec, q=q);
sens_mixture_mean_vec=apply(mean_mat, 1, sens_formula_vec, q=q);
sens_mixture_prob_vec=apply(prob_mat, 1, sens_formula_vec, q=q);

matplot(sample_size, cbind(sens_mixture_raw_vec, sens_mixture_mean_vec, sens_mixture_prob_vec));

mean(sens_mixture_raw_vec);
mean(sens_mixture_mean_vec);
mean(sens_mixture_prob_vec);


# calculat the sensitivies on the samples size when K is changing;
# assume mu=0, tau_square=1, and change K;
# calculate the rho for each hospital;
q=0.9;
# q=0.75;
prob_class=0.9;
# prob_class=0.9;
# prob_class=0.7;
# prob_class=0.6;
# prob_class=0.46;

# process_data=process_data[process_data$r1 < 30,];
# process_data$r1=c(3,50);


# sample_size=process_data$r1;
N=length(sample_size);

sample_repeat=table(sort(sample_size));
begin_index=rep(1, length(sample_repeat));
end_index=rep(0, length(sample_repeat));

for (i in 2:length(sample_repeat))
{
 begin_index[i]=begin_index[i-1]+sample_repeat[i-1];
}

for (i in 1:length(sample_repeat))
{
 end_index[i]=begin_index[i]+sample_repeat[i]-1;
}
sample_size_distinct=sort(sample_size)[begin_index];


q=0.9;
# mu_actual=0.7;
# tau_square_actual=0.04;
# sigma_square_actual=0.26;
# K_actual=0.4;
mu_actual=-3.48;
tau_square_actual=0.29;
sigma_square_actual=2.31;
K_actual=0.36;

mu=mu_actual;
tau_square=tau_square_actual;
# we could choose mu=0 and tau_square=1;
K_vec=seq(0.2,2.0,0.02);
length_K=length(K_vec);
sens_raw_vec=rep(0, length_K);
sens_mean_vec=rep(0, length_K);
sens_prob_vec=rep(0, length_K);
sens_prob2_vec=rep(0, length_K);


spec_raw_vec=rep(0, length_K);
spec_mean_vec=rep(0, length_K);
spec_prob_vec=rep(0, length_K);
spec_prob2_vec=rep(0, length_K);

# set q=0, 1/3*q, and 2/3*q, and q;
cut_prob=mu+qnorm(q)*sqrt(tau_square);
# cut_prob=mu;

for (i in 1:length_K)
{
true_mean=mu;
K=K_vec[i];
sigma_square=tau_square/K^2;

true_rho_vec=sqrt(K^2/(K^2+1/sample_size));
marginal_var_raw=tau_square+sigma_square/sample_size;
marginal_var_mean=tau_square^2/marginal_var_raw;
marginal_var_prob=marginal_var_mean;
marginal_var_prob2=marginal_var_mean;

raw_mixture=norMix(rep(true_mean, N), sig2=marginal_var_raw);
quantile_raw_mixture=qnorMix(p=q, obj=raw_mixture);

mean_mixture=norMix(rep(true_mean, N), sig2=marginal_var_mean);
quantile_mean_mixture=qnorMix(p=q, obj=mean_mixture);

prob_mixture=norMix(true_mean-qnorm(prob_class)*sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size), sig2=marginal_var_prob);
quantile_prob_mixture=qnorMix(p=q, obj=prob_mixture);


prob2_mixture=norMix((true_mean-cut_prob)/sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size), sig2=sample_size*K^2);
quantile_prob2_mixture=qnorMix(p=q, obj=prob2_mixture);


cut_off_raw=(quantile_raw_mixture-true_mean)/sqrt(marginal_var_raw);
cut_off_mean=(quantile_mean_mixture-true_mean)/sqrt(marginal_var_mean);
cut_off_prob=(quantile_prob_mixture-true_mean+qnorm(prob_class)*sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size))/sqrt(marginal_var_prob);
cut_off_prob2=(quantile_prob2_mixture*sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size)-true_mean+cut_prob)/sqrt(marginal_var_prob2);


raw_mat=cbind(cut_off_raw, true_rho_vec);
mean_mat=cbind(cut_off_mean, true_rho_vec);
prob_mat=cbind(cut_off_prob, true_rho_vec);
prob2_mat=cbind(cut_off_prob2, true_rho_vec);


sens_mixture_raw_vec=apply(raw_mat, 1, sens_formula_vec, q=q);
sens_mixture_mean_vec=apply(mean_mat, 1, sens_formula_vec, q=q);
sens_mixture_prob_vec=apply(prob_mat, 1, sens_formula_vec, q=q);
sens_mixture_prob2_vec=apply(prob2_mat, 1, sens_formula_vec, q=q);

spec_mixture_raw_vec=apply(raw_mat, 1, spec_formula_vec, q=q);
spec_mixture_mean_vec=apply(mean_mat, 1, spec_formula_vec, q=q);
spec_mixture_prob_vec=apply(prob_mat, 1, spec_formula_vec, q=q);
spec_mixture_prob2_vec=apply(prob2_mat, 1, spec_formula_vec, q=q);



# matplot(sample_size, cbind(sens_mixture_raw_vec, sens_mixture_mean_vec, sens_mixture_prob_vec));

sens_raw_vec[i]=mean(sens_mixture_raw_vec);
sens_mean_vec[i]=mean(sens_mixture_mean_vec);
sens_prob_vec[i]=mean(sens_mixture_prob_vec);
sens_prob2_vec[i]=mean(sens_mixture_prob2_vec);

spec_raw_vec[i]=mean(spec_mixture_raw_vec);
spec_mean_vec[i]=mean(spec_mixture_mean_vec);
spec_prob_vec[i]=mean(spec_mixture_prob_vec);
spec_prob2_vec[i]=mean(spec_mixture_prob2_vec);

}

# the following is Fig. 3 top panel;
# pdf("picture.pdf",height=5,width=5)


lends = c("DIR", "SHR");
matplot(K_vec, cbind(sens_raw_vec, sens_mean_vec), xlab="K=SD(BW)/SD(WI)", ylab="Sensitivity for classifying top 10%", type="p", pch=c(0,1));
points(cbind(1.75, c(.63, .66)), pch=c(0,1), col= 1:2, cex = 1.25);
text(cbind(1.9, c(.63, .66)), lends, col= 1:2, cex = 1);

matplot(matrix(1:12, 4), type="c", lty=1, lwd=10, lend=lends)
text(cbind(2.5, 2*c(1,3,5)-.4), lends, col= 1:3, cex = 1.5)
lends <- c("round","butt","square")
matplot(matrix(1:12, 4), type="c", lty=1, lwd=10, lend=lends)
text(cbind(2.5, 2*c(1,3,5)-.4), lends, col= 1:3, cex = 1.5)


dev.off()

matplot(K_vec, cbind(sens_raw_vec, sens_mean_vec, sens_prob_vec, sens_prob2_vec), main="1/2/3: DIR/SHR/PROB", xlab="K=SD(BW)/SD(WI)", ylab="Sensitivity for classifying top 10%", type="l");

plot(K_vec, sens_prob_vec-sens_mean_vec);

matplot(K_vec, cbind(spec_raw_vec, spec_mean_vec, spec_prob_vec, spec_prob2_vec), main="1/2/3: MLE/BMEAN/BPROB", xlab="K=SD(BW)/SD(WI)", ylab="Specificity for classifying top 10%");

cbind(sens_raw_vec, sens_mean_vec, sens_prob_vec, sens_prob2_vec);

# calculate sensitivity/specificity for the shrinkage and posterior probability methods;
# for the posterior probability methods using different thresholds;
q=0.9;


# mu_actual=0.7;
# tau_square_actual=0.04;
# sigma_square_actual=0.26;
# K_actual=0.4;
mu_actual=-3.48;
tau_square_actual=0.29;
sigma_square_actual=2.31;
K_actual=0.36;

mu=mu_actual;
tau_square=tau_square_actual;



# we could choose mu=0 and tau_square=1;
K_vec=seq(0.2,2.0,0.02);
length_K=length(K_vec);

sens_mean_vec=rep(0, length_K);
spec_mean_vec=rep(0, length_K);


for (i in 1:length_K)
{
true_mean=mu;
K=K_vec[i];
sigma_square=tau_square/K^2;

true_rho_vec=sqrt(K^2/(K^2+1/sample_size));

marginal_var_raw=tau_square+sigma_square/sample_size;
marginal_var_mean=tau_square^2/marginal_var_raw;

mean_mixture=norMix(rep(true_mean, N), sig2=marginal_var_mean);
quantile_mean_mixture=qnorMix(p=q, obj=mean_mixture);

cut_off_mean=(quantile_mean_mixture-true_mean)/sqrt(marginal_var_mean);
mean_mat=cbind(cut_off_mean, true_rho_vec);
sens_mixture_mean_vec=apply(mean_mat, 1, sens_formula_vec, q=q);
spec_mixture_mean_vec=apply(mean_mat, 1, spec_formula_vec, q=q);

sens_mean_vec[i]=mean(sens_mixture_mean_vec);
spec_mean_vec[i]=mean(spec_mixture_mean_vec);

}


sens_prob_vec=rep(0, length_K);
sens_prob2_vec=rep(0, length_K);
spec_prob_vec=rep(0, length_K);
spec_prob2_vec=rep(0, length_K);

sens_prob_mat=matrix(0, length_K, 5);
sens_prob2_mat=matrix(0, length_K, 5);
spec_prob_mat=matrix(0, length_K, 5);
spec_prob2_mat=matrix(0, length_K, 5);


prob_class_vec=c(0.6, 0.7, 0.8, 0.9, 0.95);
q_vec=c(0.25, 0.5, 0.75, 0.9, 0.95);


for (s in 1:5)
{

cut_prob=mu+qnorm(q_vec[s])*sqrt(tau_square);
prob_class=prob_class_vec[s];

for (i in 1:length_K)
{
true_mean=mu;



K=K_vec[i];
sigma_square=tau_square/K^2;


true_rho_vec=sqrt(K^2/(K^2+1/sample_size));
marginal_var_raw=tau_square+sigma_square/sample_size;
marginal_var_mean=tau_square^2/marginal_var_raw;
marginal_var_prob=marginal_var_mean;
marginal_var_prob2=marginal_var_mean;

prob_mixture=norMix(true_mean-qnorm(prob_class)*sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size), sig2=marginal_var_prob);
quantile_prob_mixture=qnorMix(p=q, obj=prob_mixture);

prob2_mixture=norMix((true_mean-cut_prob)/sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size), sig2=sample_size*K^2);
quantile_prob2_mixture=qnorMix(p=q, obj=prob2_mixture);


cut_off_prob=(quantile_prob_mixture-true_mean+qnorm(prob_class)*sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size))/sqrt(marginal_var_prob);
cut_off_prob2=(quantile_prob2_mixture*sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size)-true_mean+cut_prob)/sqrt(marginal_var_prob2);


prob_mat=cbind(cut_off_prob, true_rho_vec);
prob2_mat=cbind(cut_off_prob2, true_rho_vec);


sens_mixture_prob_vec=apply(prob_mat, 1, sens_formula_vec, q=q);
sens_mixture_prob2_vec=apply(prob2_mat, 1, sens_formula_vec, q=q);

spec_mixture_prob_vec=apply(prob_mat, 1, spec_formula_vec, q=q);
spec_mixture_prob2_vec=apply(prob2_mat, 1, spec_formula_vec, q=q);



sens_prob_vec[i]=mean(sens_mixture_prob_vec);
sens_prob2_vec[i]=mean(sens_mixture_prob2_vec);

spec_prob_vec[i]=mean(spec_mixture_prob_vec);
spec_prob2_vec[i]=mean(spec_mixture_prob2_vec);

}

sens_prob_mat[,s]=sens_prob_vec;
sens_prob2_mat[,s]=sens_prob2_vec;
spec_prob_mat[,s]=spec_prob_vec;
spec_prob2_mat[,s]=spec_prob2_vec;

}

matplot(K_vec, cbind(sens_mean_vec, sens_prob_mat[,5]), main="1/2/3: MLE/BMEAN/BPROB", xlab="K=SD(BW)/SD(WI)", ylab="Specificity for classifying top 10%");

# Fig. 3 bottom left
matplot(K_vec, sens_prob_mat-sens_mean_vec, xlab="K=SD(BW)/SD(WI)", ylab="Sensitivity difference between PROB1 and SHR", type="p", pch=rep(2,4));
abline(h=0);
lends= c("PROB1,Pprob=.6", "PROB1,Pprob=.7","PROB1,Pprob=.8","PROB1,Pprob=.9","PROB1,Pprob=.95");

points(cbind(1.45, c(-0.011, -0.0105, -0.010, -0.0095, -0.009)), pch=rep(2,5), col= 1:5, cex = 1.25);
text(cbind(1.75, c(-0.011, -0.0105, -0.010, -0.0095, -0.009)), lends, col= 1:5, cex = 1);


matplot(K_vec, spec_prob_mat-spec_mean_vec);

# Fig. 3 bottom right;
matplot(K_vec, sens_prob2_mat-sens_mean_vec, xlab="K=SD(BW)/SD(WI)", ylab="Sensitivity difference between PROB2 and SHR", type="p", pch=rep(3,5));
abline(h=0);
lends= c("PROB2,s=.25", "PROB2,s=.5","PROB2,s=.75","PROB2,s=.9","PROB2,s=.95");

points(cbind(.25, c(-0.17, -0.16, -0.15, -0.14, -0.13)), pch=rep(3,5), col= 1:5, cex = 1.25);
text(cbind(.5, c(-0.17, -0.16, -0.15, -0.14, -0.13)), lends, col= 1:5, cex = 1);


matplot(K_vec, spec_prob2_mat-spec_mean_vec);

sens_prob_mat;


# evaluat the average sensitivity and contribution from each hospital of three methods;
q=0.9;
prob_class=0.9;
N=length(sample_size);

# mu_actual=0.7;
# tau_square_actual=0.04;
# sigma_square_actual=0.26;
mu_actual=-3.48;
tau_square_actual=0.29;
sigma_square_actual=2.31;


K_actual=sqrt(tau_square_actual/sigma_square_actual);
mu=mu_actual;
tau_square=tau_square_actual;


cut_prob=mu+qnorm(q)*sqrt(tau_square);

true_mean=mu_actual;
K=K_actual;
sigma_square=sigma_square_actual;
tau_square=sigma_square*K^2;

true_rho_vec=sqrt(K^2/(K^2+1/sample_size));
marginal_var_raw=tau_square+sigma_square/sample_size;
marginal_var_mean=tau_square^2/marginal_var_raw;
marginal_var_prob=marginal_var_mean;
marginal_var_prob2=marginal_var_mean;


raw_mixture=norMix(rep(true_mean, N), sig2=marginal_var_raw);
quantile_raw_mixture=qnorMix(p=q, obj=raw_mixture);

mean_mixture=norMix(rep(true_mean, N), sig2=marginal_var_mean);
quantile_mean_mixture=qnorMix(p=q, obj=mean_mixture);

prob_mixture=norMix(true_mean-qnorm(prob_class)*sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size), sig2=marginal_var_prob);
quantile_prob_mixture=qnorMix(p=q, obj=prob_mixture);

prob2_mixture=norMix((true_mean-cut_prob)/sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size), sig2=sample_size*K^2);
quantile_prob2_mixture=qnorMix(p=q, obj=prob2_mixture);



cut_off_raw=(quantile_raw_mixture-true_mean)/sqrt(marginal_var_raw);
cut_off_mean=(quantile_mean_mixture-true_mean)/sqrt(marginal_var_mean);
cut_off_prob=(quantile_prob_mixture-true_mean+qnorm(prob_class)*sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size))/sqrt(marginal_var_prob);
cut_off_prob2=(quantile_prob2_mixture*sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size)-true_mean+cut_prob)/sqrt(marginal_var_prob2);



raw_mat=cbind(cut_off_raw, true_rho_vec);
mean_mat=cbind(cut_off_mean, true_rho_vec);
prob_mat=cbind(cut_off_prob, true_rho_vec);
prob2_mat=cbind(cut_off_prob2, true_rho_vec);


sens_mixture_raw_vec=apply(raw_mat, 1, sens_formula_vec, q=q);
sens_mixture_mean_vec=apply(mean_mat, 1, sens_formula_vec, q=q);
sens_mixture_prob_vec=apply(prob_mat, 1, sens_formula_vec, q=q);
sens_mixture_prob2_vec=apply(prob2_mat, 1, sens_formula_vec, q=q);


spec_mixture_raw_vec=apply(raw_mat, 1, spec_formula_vec, q=q);
spec_mixture_mean_vec=apply(mean_mat, 1, spec_formula_vec, q=q);
spec_mixture_prob_vec=apply(prob_mat, 1, spec_formula_vec, q=q);
spec_mixture_prob2_vec=apply(prob2_mat, 1, spec_formula_vec, q=q);

# Fig. 4;

lends= c("DIR", "SHR","PROB1","PROB2");


matplot(sample_size, cbind(sens_mixture_raw_vec, sens_mixture_mean_vec, sens_mixture_prob_vec, sens_mixture_prob2_vec), main="", xlab="sample size", ylab="Sample size-specific contribution to sensitivity", type="p", pch=c(0,1,2,3), lwd=1);
points(cbind(180, c(0, .07, .14, .21)), pch=c(0,1,2,3), col= 1:4, cex = 1.25);
text(cbind(200, c(0, .07, .14, .21)), lends, col= 1:4, cex = 1);



matplot(sample_size, cbind(spec_mixture_raw_vec, spec_mixture_mean_vec, spec_mixture_prob_vec, spec_mixture_prob2_vec), main="", xlab="sample size", ylab="Sample size-specific contribution to specificity", type="p", pch=c(0,1,2,3));
points(cbind(180, c(0.7, .72, .74, .76)), pch=c(0,1,2,3), col= 1:4, cex = 1.25);
text(cbind(200, c(0.7, .72, .74, .76)), lends, col= 1:4, cex = 1);


sens_matrix=cbind(sample_size, sens_mixture_raw_vec, sens_mixture_mean_vec, sens_mixture_prob_vec, sens_mixture_prob2_vec);

sens_matrix[order(sample_size),];

plot(sample_size, sens_mixture_prob_vec);

mean(sens_mixture_raw_vec);
mean(sens_mixture_mean_vec);
mean(sens_mixture_prob_vec);
mean(sens_mixture_prob2_vec);

mean(spec_mixture_raw_vec);
mean(spec_mixture_mean_vec);
mean(spec_mixture_prob_vec);
mean(spec_mixture_prob2_vec);



# plot the relative probability of identifying hospitals with sample size n;

relative_prob_raw=1/N*pnorm(-cut_off_raw)/(1-q);
relative_prob_mean=1/N*pnorm(-cut_off_mean)/(1-q);
relative_prob_prob=1/N*pnorm(-cut_off_prob)/(1-q);

relative_prob_raw_sort=relative_prob_raw[order(sample_size)];
relative_prob_raw_true=relative_prob_raw_sort[begin_index]*table(sample_size);

relative_prob_mean_sort=relative_prob_mean[order(sample_size)];
relative_prob_mean_true=relative_prob_mean_sort[begin_index]*table(sample_size);

relative_prob_prob_sort=relative_prob_prob[order(sample_size)];
relative_prob_prob_true=relative_prob_prob_sort[begin_index]*table(sample_size);

matplot(sample_size, cbind(relative_prob_raw, relative_prob_mean, relative_prob_prob));

matplot(sample_size_distinct, cbind(relative_prob_raw_true, relative_prob_mean_true, relative_prob_prob_true), xlab="distinct sample size", ylab="probability", main="1/2/3: MLE/BMEAN/BPROB");

# calculate the selection probability into the bottom tier;
relative_prob_bot_raw=1/N*pnorm(cut_off_raw)/q;
relative_prob_bot_mean=1/N*pnorm(cut_off_mean)/q;
relative_prob_bot_prob=1/N*pnorm(cut_off_prob)/q;

relative_prob_raw_sort=relative_prob_bot_raw[order(sample_size)];
relative_prob_raw_true=relative_prob_raw_sort[begin_index]*table(sample_size);

relative_prob_mean_sort=relative_prob_bot_mean[order(sample_size)];
relative_prob_mean_true=relative_prob_mean_sort[begin_index]*table(sample_size);

relative_prob_prob_sort=relative_prob_bot_prob[order(sample_size)];
relative_prob_prob_true=relative_prob_prob_sort[begin_index]*table(sample_size);

matplot(sample_size_distinct, cbind(relative_prob_raw_true, relative_prob_mean_true, relative_prob_prob_true), xlab="distinct sample size", ylab="probability", main="1/2/3: MLE/BMEAN/BPROB");

# Evaluate the sensivity/specificity from all 4 methods when q changes;
# evaluat the average sensitivity and contribution from each hospital of three methods;
# mu_actual=0.7;
# tau_square_actual=0.04;
# sigma_square_actual=0.26;
mu_actual=-3.48;
tau_square_actual=0.29;
sigma_square_actual=2.31;

K_actual=sqrt(tau_square_actual/sigma_square_actual);
mu=mu_actual;
tau_square=tau_square_actual;
N=length(sample_size);


q_vec=seq(0.5, 0.9, 0.01);
# q_vec=seq(0.001, 0.999, .001);
prob_class=0.9;
cut_prob_vec=mu_actual+qnorm(q_vec)*sqrt(tau_square);

sens_mixture_raw=rep(0, length(q_vec));
sens_mixture_mean=rep(0, length(q_vec));
sens_mixture_prob=rep(0, length(q_vec));
sens_mixture_prob2=rep(0, length(q_vec));

spec_mixture_raw=rep(0, length(q_vec));
spec_mixture_mean=rep(0, length(q_vec));
spec_mixture_prob=rep(0, length(q_vec));
spec_mixture_prob2=rep(0, length(q_vec));



true_mean=mu_actual;
K=K_actual;
sigma_square=sigma_square_actual;
tau_square=sigma_square*K^2;

true_rho_vec=sqrt(K^2/(K^2+1/sample_size));
marginal_var_raw=tau_square+sigma_square/sample_size;
marginal_var_mean=tau_square^2/marginal_var_raw;
marginal_var_prob=marginal_var_mean;
marginal_var_prob2=marginal_var_mean;

for (i in 1:length(q_vec))
{

raw_mixture=norMix(rep(true_mean, N), sig2=marginal_var_raw);
quantile_raw_mixture=qnorMix(p=q_vec[i], obj=raw_mixture);

mean_mixture=norMix(rep(true_mean, N), sig2=marginal_var_mean);
quantile_mean_mixture=qnorMix(p=q_vec[i], obj=mean_mixture);

prob_mixture=norMix(true_mean-qnorm(prob_class)*sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size), sig2=marginal_var_prob);
quantile_prob_mixture=qnorMix(p=q_vec[i], obj=prob_mixture);

prob2_mixture=norMix((true_mean-cut_prob_vec[i])/sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size), sig2=sample_size*K^2);
quantile_prob2_mixture=qnorMix(p=q_vec[i], obj=prob2_mixture);



cut_off_raw=(quantile_raw_mixture-true_mean)/sqrt(marginal_var_raw);
cut_off_mean=(quantile_mean_mixture-true_mean)/sqrt(marginal_var_mean);
cut_off_prob=(quantile_prob_mixture-true_mean+qnorm(prob_class)*sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size))/sqrt(marginal_var_prob);
cut_off_prob2=(quantile_prob2_mixture*sqrt(K^2/(K^2+1/sample_size)*sigma_square/sample_size)-true_mean+cut_prob_vec[i])/sqrt(marginal_var_prob2);



raw_mat=cbind(cut_off_raw, true_rho_vec);
mean_mat=cbind(cut_off_mean, true_rho_vec);
prob_mat=cbind(cut_off_prob, true_rho_vec);
prob2_mat=cbind(cut_off_prob2, true_rho_vec);


sens_mixture_raw_vec=apply(raw_mat, 1, sens_formula_vec, q=q_vec[i]);
sens_mixture_mean_vec=apply(mean_mat, 1, sens_formula_vec, q=q_vec[i]);
sens_mixture_prob_vec=apply(prob_mat, 1, sens_formula_vec, q=q_vec[i]);
sens_mixture_prob2_vec=apply(prob2_mat, 1, sens_formula_vec, q=q_vec[i]);


spec_mixture_raw_vec=apply(raw_mat, 1, spec_formula_vec, q=q_vec[i]);
spec_mixture_mean_vec=apply(mean_mat, 1, spec_formula_vec, q=q_vec[i]);
spec_mixture_prob_vec=apply(prob_mat, 1, spec_formula_vec, q=q_vec[i]);
spec_mixture_prob2_vec=apply(prob2_mat, 1, spec_formula_vec, q=q_vec[i]);


sens_mixture_raw[i]=mean(sens_mixture_raw_vec);
sens_mixture_mean[i]=mean(sens_mixture_mean_vec);
sens_mixture_prob[i]=mean(sens_mixture_prob_vec);
sens_mixture_prob2[i]=mean(sens_mixture_prob2_vec);

spec_mixture_raw[i]=mean(spec_mixture_raw_vec);
spec_mixture_mean[i]=mean(spec_mixture_mean_vec);
spec_mixture_prob[i]=mean(spec_mixture_prob_vec);
spec_mixture_prob2[i]=mean(spec_mixture_prob2_vec);


}


# Fig. 5;
lends= c("DIR", "SHR","PROB1","PROB2");


matplot(q_vec, cbind(sens_mixture_raw, sens_mixture_mean, sens_mixture_prob, sens_mixture_prob2), xlab="cut-off point", ylab="sensitivity", pch=c(0,1,2,3) );

points(cbind(0.5, c(0.76, .77, .78, .79)), pch=c(0,1,2,3), col= 1:4, cex = 1.25);
text(cbind(0.55, c(0.76, .77, .78, .79)), lends, col= 1:4, cex = 1);



matplot(q_vec, cbind(spec_mixture_raw, spec_mixture_mean, spec_mixture_prob, spec_mixture_prob2), xlab="cut-off point", ylab="specificity", pch=c(0,1,2,3));
points(cbind(0.83, c(.890, .897, .904, .911)), pch=c(0,1,2,3), col= 1:4, cex = 1.25);
text(cbind(0.88, c(.890, .897, .904, .911)), lends, col= 1:4, cex = 1);


# caclulate the total correct classification rate;
pc_mixture_raw=(1-q_vec)*sens_mixture_raw+q_vec*spec_mixture_raw;
pc_mixture_mean=(1-q_vec)*sens_mixture_mean+q_vec*spec_mixture_mean;
pc_mixture_prob=(1-q_vec)*sens_mixture_prob+q_vec*spec_mixture_prob;
pc_mixture_prob2=(1-q_vec)*sens_mixture_prob2+q_vec*spec_mixture_prob2;



plot(1-spec_mixture_raw, sens_mixture_raw);

plot(1-spec_mixture_mean, sens_mixture_mean);

matplot(cbind(1-spec_mixture_raw, 1-spec_mixture_mean, 1-spec_mixture_prob, 1-spec_mixture_prob2), cbind(sens_mixture_raw, sens_mixture_mean, sens_mixture_prob, sens_mixture_prob2));


# caclulate the sensitivity/specificity for hospitals when removing sample size one by one;
# sample_size=process_data$r1;
N=length(sample_size);
q=0.9;
prob_class=0.9;

sample_repeat=table(sort(sample_size));
begin_index=rep(1, length(sample_repeat));
end_index=rep(0, length(sample_repeat));

for (i in 2:length(sample_repeat))
{
 begin_index[i]=begin_index[i-1]+sample_repeat[i-1];
}

for (i in 1:length(sample_repeat))
{
 end_index[i]=begin_index[i]+sample_repeat[i]-1;
}
sample_size_distinct=sort(sample_size)[begin_index];

sample_size_sort=sort(sample_size);


# estimate the model parameters by removing the sample one by one;
distinct_sample_size=length(sample_size_distinct);

y_bar_sort=y_bar[order(sample_size)];
mu0_vec=rep(0, distinct_sample_size);
tau2_vec=rep(0, distinct_sample_size);
sigma2_vec=rep(0, distinct_sample_size);
K_vec=rep(0, distinct_sample_size);
cut_prob_vec=rep(0, distinct_sample_size);

for (i in 1:distinct_sample_size)
{
M=length(sample_size_sort[begin_index[i]:length(sample_size_sort)]);
y_bar=y_bar_sort[begin_index[i]:length(sample_size_sort)];
sample_size=sample_size_sort[begin_index[i]:length(sample_size_sort)];
process_stage2_data_simulate=list("M", "sample_size", "y_bar");
process_stage2_init= function(){list(prec1 = 1, prec0=1, mu0 = 0)};
process_stage2_parameters=c("mu0", "sigma2", "tau2", "K");

gibbs_no=20000;


process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file="C:/Users/WDQ7/reliability_test_forR.txt",
    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=1, bugs.seed=197789,
    bugs.directory="C:/Users/WDQ7/Desktop/WinBUGS14", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);



attach.bugs(process_stage2.sim);
rm(process_stage2.sim);

mu0_vec[i]=median(mu0);
tau2_vec[i]=median(tau2);
sigma2_vec[i]=median(sigma2);
K_vec[i]=median(K);
# cut_prob_vec[i]=mu0_vec[i]+qnorm(q)*sqrt(tau2_vec[i]);
cut_prob_vec[i]=-mu0_vec[i]+qnorm(q)*sqrt(tau2_vec[i]);


}


sens_losingsample_raw=rep(0, length(sample_size_distinct));
sens_losingsample_mean=rep(0, length(sample_size_distinct));
sens_losingsample_prob=rep(0, length(sample_size_distinct));
sens_losingsample_prob2=rep(0, length(sample_size_distinct));

spec_losingsample_raw=rep(0, length(sample_size_distinct));
spec_losingsample_mean=rep(0, length(sample_size_distinct));
spec_losingsample_prob=rep(0, length(sample_size_distinct));
spec_losingsample_prob2=rep(0, length(sample_size_distinct));


# evaluat the average sensitivity and contribution from each hospital of three methods;
for (i in 1:distinct_sample_size)
{
# true_mean=mu0_vec[i];
true_mean=-mu0_vec[i];
# K=K_vec[i];
sigma_square=sigma2_vec[i];
tau_square=tau2_vec[i];
K=sqrt(tau_square)/sqrt(sigma_square);
N=length(sample_size_sort[begin_index[i]:length(sample_size_sort)]);
cut_prob=cut_prob_vec[i];

true_rho_vec=sqrt(K^2/(K^2+1/sample_size_sort[begin_index[i]:length(sample_size_sort)]));
marginal_var_raw=tau_square+sigma_square/sample_size_sort[begin_index[i]:length(sample_size_sort)];
marginal_var_mean=tau_square^2/marginal_var_raw;
marginal_var_prob=marginal_var_prob2=marginal_var_mean;

raw_mixture=norMix(rep(true_mean, N), sig2=marginal_var_raw);
quantile_raw_mixture=qnorMix(p=q, obj=raw_mixture);

mean_mixture=norMix(rep(true_mean, N), sig2=marginal_var_mean);
quantile_mean_mixture=qnorMix(p=q, obj=mean_mixture);

prob_mixture=norMix(true_mean-qnorm(prob_class)*sqrt(K^2/(K^2+1/sample_size_sort[begin_index[i]:length(sample_size_sort)])*sigma_square/sample_size_sort[begin_index[i]:length(sample_size_sort)]), sig2=marginal_var_prob);
quantile_prob_mixture=qnorMix(p=q, obj=prob_mixture);

prob2_mixture=norMix((true_mean-cut_prob)/sqrt(K^2/(K^2+1/sample_size_sort[begin_index[i]:length(sample_size_sort)])*sigma_square/sample_size_sort[begin_index[i]:length(sample_size_sort)]), sig2=sample_size_sort[begin_index[i]:length(sample_size_sort)]*K^2);

quantile_prob2_mixture=qnorMix(p=q, obj=prob2_mixture);



cut_off_raw=(quantile_raw_mixture-true_mean)/sqrt(marginal_var_raw);
cut_off_mean=(quantile_mean_mixture-true_mean)/sqrt(marginal_var_mean);
cut_off_prob=(quantile_prob_mixture-true_mean+qnorm(prob_class)*sqrt(K^2/(K^2+1/sample_size_sort[begin_index[i]:length(sample_size_sort)])*sigma_square/sample_size_sort[begin_index[i]:length(sample_size_sort)]))/sqrt(marginal_var_prob);
cut_off_prob2=(quantile_prob2_mixture*sqrt(K^2/(K^2+1/sample_size_sort[begin_index[i]:length(sample_size_sort)])*sigma_square/sample_size_sort[begin_index[i]:length(sample_size_sort)])-true_mean+cut_prob)/sqrt(marginal_var_prob2);


raw_mat=cbind(cut_off_raw, true_rho_vec);
mean_mat=cbind(cut_off_mean, true_rho_vec);
prob_mat=cbind(cut_off_prob, true_rho_vec);
prob2_mat=cbind(cut_off_prob2, true_rho_vec);


sens_mixture_raw_vec=apply(raw_mat, 1, sens_formula_vec, q=q);
sens_mixture_mean_vec=apply(mean_mat, 1, sens_formula_vec, q=q);
sens_mixture_prob_vec=apply(prob_mat, 1, sens_formula_vec, q=q);
sens_mixture_prob2_vec=apply(prob2_mat, 1, sens_formula_vec, q=q);


spec_mixture_raw_vec=apply(raw_mat, 1, spec_formula_vec, q=q);
spec_mixture_mean_vec=apply(mean_mat, 1, spec_formula_vec, q=q);
spec_mixture_prob_vec=apply(prob_mat, 1, spec_formula_vec, q=q);
spec_mixture_prob2_vec=apply(prob2_mat, 1, spec_formula_vec, q=q);

sens_losingsample_raw[i]=mean(sens_mixture_raw_vec);

sens_losingsample_mean[i]=mean(sens_mixture_mean_vec);

sens_losingsample_prob[i]=mean(sens_mixture_prob_vec);

sens_losingsample_prob2[i]=mean(sens_mixture_prob2_vec);

spec_losingsample_raw[i]=mean(spec_mixture_raw_vec);

spec_losingsample_mean[i]=mean(spec_mixture_mean_vec);

spec_losingsample_prob[i]=mean(spec_mixture_prob_vec);

spec_losingsample_prob2[i]=mean(spec_mixture_prob2_vec);




}

matplot(sample_size_distinct[1:30], cbind(sens_losingsample_raw, sens_losingsample_mean, sens_losingsample_prob, sens_losingsample_prob2)[1:30,], xlab="Sample size excluded", ylab="Sensitivity", main="1/2/3: MLE/BMEAN/BPROB");

matplot(sample_size_distinct[1:30], cbind(spec_losingsample_raw, spec_losingsample_mean, spec_losingsample_prob)[1:30,], xlab="Sample size excluded", ylab="Specificity", main="1/2/3: MLE/BMEAN/BPROB");

# Fig. 6
plot(sample_size_distinct[1:30], sens_losingsample_prob2[1:30], xlab="Sample size excluded", ylab="Sensitivity Calculated from PROB2 Method", main="", pch=3);
plot(sample_size_distinct[1:30], spec_losingsample_mean[1:30], xlab="Sample size excluded", ylab="Specificity", main="");

lends= c("PROB2");


plot(sample_size_distinct[1:18], sens_losingsample_prob2[1:18], xlab="Sample size excluded", ylab="Sensitivity Calculated from PROB2 Method", main="", pch=3, type="b");

points(cbind(25, .73), pch=3, col= 1, cex = 1.25);
text(cbind(27, .73), lends, col= 1, cex = 1);


plot(sample_size_distinct[1:30], sens_losingsample_prob2[1:30], xlab="Sample size excluded", ylab="Sensitivity", main="");
plot(sample_size_distinct[1:30], sens_losingsample_mean[1:30], xlab="Sample size excluded", ylab="Sensitivity", main="");



# generate simulated samples for illustration Fig. 1;
set.seed(197789);
mu_simulate=rnorm(329, mean=3.48, sd=sqrt(0.29));
ybar_simulate=mu_simulate+rnorm(329, sd=sqrt(2.31/sample_size));
mu_density=density(mu_simulate, bw=0.3);
dy=density(x=mu_simulate, n=1000, bw=0.3);
plot(x=dy$x, y=dy$y, type="l", main="g=0.9, h=-0.2");
dz=density(x=ybar_simulate, n=1000, bw=0.30);
matplot(x=cbind(dy$x, dz$x), y=cbind(dy$y, dz$y), type="l", xlab="performance", ylab="density");
abline(v=quantile(mu_simulate, 0.9), col="red", lty=2);
abline(v=quantile(ybar_simulate, 0.9));

lends= c("Sample Average", "True Quality");
# lines(cbind(5, c(0.15, .2)), lty=c(1,2), col=1:2 , cex = 1.25);
lines(cbind(c(4.75,5.25), 0.15), lty=1, col=1, cex=1.25);
lines(cbind(c(4.75,5.25), 0.2), lty=2, col=2, cex=1.25);

text(cbind(5.9, c(0.15, .2)), lends, col= 1:2, cex = 0.75);


plot(mu_simulate, ybar_simulate, xlab="true qualiy", ylab="sample average");
abline(h=quantile(ybar_simulate, 0.9), v=quantile(mu_simulate, 0.9));
sens_simulate=sum(mu_simulate > quantile(mu_simulate, 0.9) & ybar_simulate > quantile(ybar_simulate, 0.9))/(329*(1-0.9));
quantile(mu_simulate, 0.9);
quantile(ybar_simulate, 0.9);

spec_simulate=sum(mu_simulate < quantile(mu_simulate, 0.9) & ybar_simulate < quantile(ybar_simulate, 0.9))/(3304*0.9);


# calculate the variations of the sensitivity/specificity under the model;
# for wait time, we need to reverse the sign;
y_bar=-process$log_waittime;

process_stage2_data_simulate=list("M", "sample_size", "y_bar");
process_stage2_init= function(){list(prec1 = 1, prec0=1, mu0 = 0)};
process_stage2_parameters=c("mu0", "sigma2", "tau2", "K");

gibbs_no=20000;


process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file="C:/Users/WDQ7/reliability_test_forR.txt",
    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=1, bugs.seed=197789,
    bugs.directory="C:/Users/WDQ7/Desktop/WinBUGS14", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);



attach.bugs(process_stage2.sim);
rm(process_stage2.sim);

median(mu0);
median(tau2);
median(sigma2);
median(sqrt(tau2/sigma2));
median(K);

prob_class=0.9;
q=0.9;

N=length(sample_size);


mu0=mu0[9001:10000];
tau2=tau2[9001:10000];
sigma2=sigma2[9001:10000];
K=K[9001:10000];

# mu0=0.7;
# tau2=0.0416;
# sigma2=0.26;
# K=0.4;
sens_mixture_raw=rep(0, length(mu0));
sens_mixture_mean=rep(0, length(mu0));
sens_mixture_prob=rep(0, length(mu0));
sens_mixture_prob2=rep(0, length(mu0));

spec_mixture_raw=rep(0, length(mu0));
spec_mixture_mean=rep(0, length(mu0));
spec_mixture_prob=rep(0, length(mu0));
spec_mixture_prob2=rep(0, length(mu0));



for (i in 1:length(mu0))
{

true_rho_vec=sqrt(K[i]^2/(K[i]^2+1/sample_size));
marginal_var_raw=tau2[i]+sigma2[i]/sample_size;
marginal_var_mean=tau2[i]^2/marginal_var_raw;
marginal_var_prob=marginal_var_mean;
marginal_var_prob2=marginal_var_mean;

cut_prob=mu0[i]+qnorm(q)*sqrt(tau2[i]);


raw_mixture=norMix(rep(mu0[i], N), sig2=marginal_var_raw);
quantile_raw_mixture=qnorMix(p=q, obj=raw_mixture);

mean_mixture=norMix(rep(mu0[i], N), sig2=marginal_var_mean);
quantile_mean_mixture=qnorMix(p=q, obj=mean_mixture);

prob_mixture=norMix(mu0[i]-qnorm(prob_class)*sqrt(K[i]^2/(K[i]^2+1/sample_size)*sigma2[i]/sample_size), sig2=marginal_var_prob);
quantile_prob_mixture=qnorMix(p=q, obj=prob_mixture);

prob2_mixture=norMix((mu0[i]-cut_prob)/sqrt(K[i]^2/(K[i]^2+1/sample_size)*sigma2[i]/sample_size), sig2=sample_size*K[i]^2);
quantile_prob2_mixture=qnorMix(p=q, obj=prob2_mixture);



cut_off_raw=(quantile_raw_mixture-mu0[i])/sqrt(marginal_var_raw);
cut_off_mean=(quantile_mean_mixture-mu0[i])/sqrt(marginal_var_mean);
cut_off_prob=(quantile_prob_mixture-mu0[i]+qnorm(prob_class)*sqrt(K[i]^2/(K[i]^2+1/sample_size)*sigma2[i]/sample_size))/sqrt(marginal_var_prob);
cut_off_prob2=(quantile_prob2_mixture*sqrt(K[i]^2/(K[i]^2+1/sample_size)*sigma2[i]/sample_size)-mu0[i]+cut_prob)/sqrt(marginal_var_prob2);



raw_mat=cbind(cut_off_raw, true_rho_vec);
mean_mat=cbind(cut_off_mean, true_rho_vec);
prob_mat=cbind(cut_off_prob, true_rho_vec);
prob2_mat=cbind(cut_off_prob2, true_rho_vec);


sens_mixture_raw_vec=apply(raw_mat, 1, sens_formula_vec, q=q);
sens_mixture_mean_vec=apply(mean_mat, 1, sens_formula_vec, q=q);
sens_mixture_prob_vec=apply(prob_mat, 1, sens_formula_vec, q=q);
sens_mixture_prob2_vec=apply(prob2_mat, 1, sens_formula_vec, q=q);


spec_mixture_raw_vec=apply(raw_mat, 1, spec_formula_vec, q=q);
spec_mixture_mean_vec=apply(mean_mat, 1, spec_formula_vec, q=q);
spec_mixture_prob_vec=apply(prob_mat, 1, spec_formula_vec, q=q);
spec_mixture_prob2_vec=apply(prob2_mat, 1, spec_formula_vec, q=q);


sens_mixture_raw[i]=mean(sens_mixture_raw_vec);
sens_mixture_mean[i]=mean(sens_mixture_mean_vec);
sens_mixture_prob[i]=mean(sens_mixture_prob_vec);
sens_mixture_prob2[i]=mean(sens_mixture_prob2_vec);

spec_mixture_raw[i]=mean(spec_mixture_raw_vec);
spec_mixture_mean[i]=mean(spec_mixture_mean_vec);
spec_mixture_prob[i]=mean(spec_mixture_prob_vec);
spec_mixture_prob2[i]=mean(spec_mixture_prob2_vec);

}

quantile(sens_mixture_raw, c(0.025, 0.5, 0.975));
quantile(sens_mixture_mean, c(0.025, 0.5, 0.975));
quantile(sens_mixture_prob, c(0.025, 0.5, 0.975));
quantile(sens_mixture_prob2, c(0.025, 0.5, 0.975));

quantile(spec_mixture_raw, c(0.025, 0.5, 0.975));
quantile(spec_mixture_mean, c(0.025, 0.5, 0.975));
quantile(spec_mixture_prob, c(0.025, 0.5, 0.975));
quantile(spec_mixture_prob2, c(0.025, 0.5, 0.975));


# study the sensitivity/specificity of DIR and SHR methods using external threshold;
# mu_actual=0.7;
# tau_square_actual=0.04;
# sigma_square_actual=0.26;
# c1_vec=seq(0.50,0.60,0.01);

mu_actual=3.48;
tau_square_actual=0.29;
sigma_square_actual=2.31;


c1_vec=seq(2.40, 4.56, 0.05);
length_c1=length(c1_vec);
K_actual=sqrt(tau_square_actual/sigma_square_actual);
K=K_actual;
true_rho_vec=sqrt(K^2/(K^2+1/sample_size));
marginal_var_raw=tau_square_actual+sigma_square_actual/sample_size;
marginal_var_mean=tau_square_actual^2/marginal_var_raw;

sens_mixture_raw=rep(0, length_c1);
sens_mixture_mean=rep(0, length_c1);
spec_mixture_raw=rep(0, length_c1);
spec_mixture_mean=rep(0, length_c1);

for (i in 1:length_c1)
{
cut_off_raw=(c1_vec[i]-mu_actual)/sqrt(marginal_var_raw);
cut_off_mean=(c1_vec[i]-mu_actual)/sqrt(marginal_var_mean);

q=pnorm((c1_vec[i]-mu_actual)/sqrt(tau_square_actual));
raw_mat=cbind(cut_off_raw, true_rho_vec);
mean_mat=cbind(cut_off_mean, true_rho_vec);

sens_mixture_raw_vec=apply(raw_mat, 1, sens_formula_vec, q=q);
sens_mixture_mean_vec=apply(mean_mat, 1, sens_formula_vec, q=q);
spec_mixture_raw_vec=apply(raw_mat, 1, spec_formula_vec, q=q);
spec_mixture_mean_vec=apply(mean_mat, 1, spec_formula_vec, q=q);

pcc_mixture_raw_vec=sens_mixture_raw_vec*(1-pnorm((c1_vec[i]-mu_actual)/sqrt(tau_square_actual)))+spec_mixture_raw_vec*pnorm((c1_vec[i]-mu_actual)/sqrt(tau_square_actual));
pcc_mixture_mean_vec=sens_mixture_mean_vec*(1-pnorm((c1_vec[i]-mu_actual)/sqrt(tau_square_actual)))+spec_mixture_mean_vec*pnorm((c1_vec[i]-mu_actual)/sqrt(tau_square_actual));



sens_mixture_raw[i]=mean(sens_mixture_raw_vec);
sens_mixture_mean[i]=mean(sens_mixture_mean_vec);
spec_mixture_raw[i]=mean(spec_mixture_raw_vec);
spec_mixture_mean[i]=mean(spec_mixture_mean_vec);


}

pcc_mixture_raw=sens_mixture_raw*(1-pnorm((c1_vec-mu_actual)/sqrt(tau_square_actual)))+spec_mixture_raw*pnorm((c1_vec-mu_actual)/sqrt(tau_square_actual));
pcc_mixture_mean=sens_mixture_mean*(1-pnorm((c1_vec-mu_actual)/sqrt(tau_square_actual)))+spec_mixture_mean*pnorm((c1_vec-mu_actual)/sqrt(tau_square_actual));

lends = c("DIR", "SHR");

matplot((c1_vec-mu_actual)/sqrt(tau_square_actual), cbind(sens_mixture_raw, sens_mixture_mean), xlab="Standardized difference between the threshold and mean", ylab="Sensitivity", pch=c(0,1));
points(cbind(-1.8, c(.60, .64)), pch=c(0,1), col= 1:2, cex = 1.25);
text(cbind(-1.5, c(.60, .64)), lends, col= 1:2, cex = 1);


matplot((c1_vec-mu_actual)/sqrt(tau_square_actual), cbind(spec_mixture_raw, spec_mixture_mean), xlab="Standardized difference between the threshold and mean", ylab="Specificity", pch=c(0,1));
points(cbind(1.5, c(.60, .64)), pch=c(0,1), col= 1:2, cex = 1.25);
text(cbind(1.8, c(.60, .64)), lends, col= 1:2, cex = 1);


matplot((c1_vec-mu_actual)/sqrt(tau_square_actual), cbind(pcc_mixture_raw, pcc_mixture_mean), xlab="Standardized difference between the threshold and mean", ylab="Probablity of correct classification", pch=c(0,1));
points(cbind(1.5, c(.895, .905)), pch=c(0,1), col= 1:2, cex = 1.25);
text(cbind(1.8, c(.895, .905)), lends, col= 1:2, cex = 1);

matplot(cbind(1-spec_mixture_raw, 1-spec_mixture_mean), cbind(sens_mixture_raw, sens_mixture_mean));

plot(1-spec_mixture_raw, sens_mixture_raw);
plot(1-spec_mixture_mean, sens_mixture_mean);



# calculate the reliability measure for binary data;
mu=log(0.7/0.3);
tau_vec=rep(0,100);
k_vec=rep(0,100);
for (i in 1:100)
{
tau=i/10;
tau_vec[i]=tau;
k=tau/sqrt((exp(mu)+exp(-mu))*(1+tau^2/2)+2);
k_vec[i]=k;
}

plot(tau_vec, k_vec);

mu=1.101;
tau=sqrt(1.583);
tau=2;
k=tau/sqrt((exp(mu)+exp(-mu))*(1+tau^2/2)+2);
k=0.7;

tau=k*sqrt((exp(mu)+exp(-mu)+2)/(1-k^2/2*(exp(mu)+exp(-mu))));


