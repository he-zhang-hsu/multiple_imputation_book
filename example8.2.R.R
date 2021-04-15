#data.gen<-function (N = 200, dep = 1, b0 = 4, b1 = -2, b2 = 0.5, b3 = -2, b4 = 2, b5 = 2, r0 = 3, r1 = -3, r2 = 0.5, r3 = -2, r4 = 1.5,r5 = 2., t.50 = 0.835, t.35 = 0.978) 
# {
rm(list=ls());
library(survival);   
library(MASS);
library(rms);
library(mice);
library(norm);

parameter_draw=function(data)
{
 # obtain all the computational components;
 p=ncol(data);
 n_total=nrow(data);
 x=cbind(rep(1,n_total), data[,1:(p-1)]);
 y=data[,p];
 # perform least-square calculation;
 beta_hat=solve(t(x)%*%x)%*%t(x)%*%y;
 s_square=1/(n_total-p)*crossprod(y-x%*%beta_hat);
 inv_v_beta=solve(t(x)%*%x);
 
 # draw sigma_square;
 # set.seed(19760421);
 sigma_square=1/(rgamma(1, shape=n_total/2-p/2, scale=2/((n_total-p)*s_square)));
 
 # draw beta;
 beta=mvrnorm(n=1, mu=beta_hat, Sigma=sigma_square*inv_v_beta);

 parameter=c(beta,sigma_square);
 parameter;
}


# First obtain the median of the complete survival time;
rowobs=1000;
b0 = -1;
b1 = 1;

r0 = -1;
r1 = 1;
var_error=1/9;

cycle_no=1000;
mi_no=50;
gibbs_no=200;

T0_median_vec=rep(NA,cycle_no);
for (cycle in 1:cycle_no)

{
set.seed(cycle);

# generate complete survival time, censoring time, and observed time;

X1 <- runif(rowobs);
# u.t is the error for the survival time;
u.t <- rnorm(rowobs, 0, sd=sqrt(var_error));

# T0 is the survival time;
T0 <- exp(b0+b1 * X1 + u.t);

T0_median_vec[cycle]=median(T0);
}

# t.50 is taken as the population median of the survival time;
t.50=mean(T0_median_vec);

t.50;

# matrices holding the parameter estimates;
# complete-data inferences;
# xyslope means regressing x on y;
BD_mean_vec=BD_mean_var_vec=BD_slope_vec=BD_slope_var_vec=BD_xyslope_vec=BD_xyslope_var_vec=rep(NA, cycle_no);
CC_mean_vec=CC_mean_var_vec=CC_slope_vec=CC_slope_var_vec=CC_xyslope_vec=CC_xyslope_var_vec=rep(NA, cycle_no);
MI_mean_mat=MI_mean_var_mat=MI_slope_mat=MI_slope_var_mat=MI_xyslope_mat=MI_xyslope_var_mat=matrix(NA, cycle_no, mi_no);
# IM_MI indicates improper MI
IM_MI_mean_mat=IM_MI_mean_var_mat=IM_MI_slope_mat=IM_MI_slope_var_mat=IM_MI_xyslope_mat=IM_MI_xyslope_var_mat=matrix(NA, cycle_no, mi_no);


for (cycle in 1:cycle_no)

{
set.seed(cycle);

# generate complete survival time, censoring time, and observed time;

X1 <- runif(rowobs);
# u.t is the error for the survival time;
u.t <- rnorm(rowobs, 0, sd=sqrt(var_error));

# T0 is the survival time;
T0 <- exp(b0+b1 * X1 + u.t);

# u.c is the error for the censoring time;
#	if(dep == 0); censoring completely at random;
# T.c is the censoring time;
# T.c <- rexp(rowobs, 1);

# if (dep==1) censoring at random;
u.c= rnorm(rowobs, 0, sd=sqrt(var_error));
T.c <- exp(r0+r1 * X1 +u.c);


	
# correlation between survival and censoring time;
cortc <- cor(T0, T.c);
# cortc
# create the censoring;
# T.o is the observed time;
T.o <- rep(0, rowobs)
I.d <- rep(0, rowobs)
I.d.t <- rep(1, rowobs)
for(k in 1:rowobs) {
	T.o[k] <- min(T0[k], T.c[k])
	if(T.o[k] == T0[k])
	I.d[k] <- 1
	}

# cbind(T0, T.c, T.o, I.d);
# mean(I.d);


# before-deletion time analysis for KM estimate;
# test for survival;

# before-deletion analysis;
KM.fo <- survfit(Surv(T0, I.d.t) ~ 1)
T50.fo <- max(KM.fo$time[KM.fo$time <= t.50])
St50.fo <- KM.fo$surv[KM.fo$time == T50.fo]
se50.fo <- KM.fo$surv[KM.fo$time == T50.fo] * KM.fo$std[KM.fo$time == T50.fo]
# low50.fo <- St50.fo - 1.96 * se50.fo
# high50.fo <- St50.fo + 1.96 * se50.fo

BD_mean_vec[cycle]=St50.fo;
BD_mean_var_vec[cycle]=se50.fo^2;



# observed data analysis;
KM.po <- survfit(Surv(T.o, I.d) ~ 1)
T50.po <- max(KM.po$time[KM.po$time <= t.50])
St50.po <- KM.po$surv[KM.po$time == T50.po]
se50.po <- KM.po$surv[KM.po$time == T50.po] * KM.po$std[KM.po$time == T50.po]	
# low50.po <- St50.po - 1.96 * se50.po
# high50.po <- St50.po + 1.96 * se50.po

CC_mean_vec[cycle]=St50.po;
CC_mean_var_vec[cycle]=se50.po^2;


# imputation programs;

# take the log-transformation;
y=log(T0);
y_miss=log(T.o);
x=cbind(rep(1,rowobs), X1);
p=ncol(x);
# missing data are the ones who are censored;
miss_indi=1-I.d;
miss_seq=cbind(seq(1,rowobs,1), miss_indi)[miss_indi==1,1];
obs_no=sum(I.d);
mis_no=rowobs-obs_no;
x_mis=x[miss_indi==1,];



# imputation;
for (i in 1:mi_no)
{
# initialize the completed data;
y_completed=y_miss;
beta_matrix=matrix(NA, gibbs_no,p);
sigma_square_vector=rep(NA, gibbs_no);

set.seed(i);

for (j in 1:gibbs_no)
{
# perform least-square calculation;
 beta_hat=solve(t(x)%*%x)%*%t(x)%*%y_completed;
 s_square=1/(rowobs-p)*crossprod(y_completed-x%*%beta_hat);
 inv_v_beta=solve(t(x)%*%x);
 
 # draw sigma_square;
 # set.seed(19760421);
 sigma_square=1/(rgamma(1, shape=rowobs/2-p/2, scale=2/((rowobs-p)*s_square)));
 sigma=sqrt(sigma_square);
 sigma_square_vector[j]=sigma_square;
 # draw beta;
 beta=mvrnorm(n=1, mu=beta_hat, Sigma=sigma_square*inv_v_beta);
 beta_matrix[j,]=beta;
 # impute the missing data;

 y_impute=x_mis%*%beta+sigma*qnorm(runif(mis_no)*(1-pnorm((y_miss[miss_indi==1]-x_mis%*%beta)/sigma))+pnorm((y_miss[miss_indi==1]-x_mis%*%beta)/sigma));

 y_completed[miss_indi==1]=y_impute;

}

# multiple imputation analysis;
T.completed=exp(y_completed);
KM.mi <- survfit(Surv(T.completed, I.d.t) ~ 1)
T50.mi <- max(KM.mi$time[KM.mi$time <= t.50])
St50.mi <- KM.mi$surv[KM.mi$time == T50.mi]
se50.mi <- KM.mi$surv[KM.mi$time == T50.mi] * KM.mi$std[KM.mi$time == T50.mi]	
# low50.mi <- St50.mi - 1.96 * se50.mi
# high50.mi <- St50.mi + 1.96 * se50.mi

MI_mean_mat[cycle,i]=St50.mi;
MI_mean_var_mat[cycle,i]=se50.mi^2;

# end of multiple imputation;
}

# diagnostics;
# naive estimate;
# summary(lm(y_miss ~ X1));

# summary(sigma_square_vector);
# summary(beta_matrix[,1]);
# summary(beta_matrix[,2]);
# hist(sigma_square_vector);
# hist(beta_matrix[,1]);
# hist(beta_matrix[,2]);

# plot(sigma_square_vector);
# plot(beta_matrix[,1]);
# plot(beta_matrix[,2]);
# cbind(y_miss, y_completed, miss_indi, I.d);

# imputation data analysis;
# buckly James estimator;
# f_fo <- bj(Surv(T0, I.d.t) ~ X1, x=TRUE, y=TRUE)

# f_po <- bj(Surv(T.o, I.d) ~ X1, x=TRUE, y=TRUE)

# f_mi <- bj(Surv(T.completed, I.d.t) ~ X1, x=TRUE, y=TRUE)

cat("the cycle is", cycle, "\n");

# end of cycle;

}


#############################################################################
# marginal mean estimand;
true_mean=mean(BD_mean_vec);
true_mean;

# performance
BD_mean_mse=mean((BD_mean_vec-true_mean)^2);
BD_mean_mse;

# lower and upper 95% CI;
BD_mean_low95_vec=BD_mean_vec-1.96*sqrt(BD_mean_var_vec);
BD_mean_up95_vec=BD_mean_vec+1.96*sqrt(BD_mean_var_vec);

BD_mean_length=mean(BD_mean_up95_vec-BD_mean_low95_vec);
BD_mean_length;

BD_mean_cov=(BD_mean_low95_vec < true_mean)*(BD_mean_up95_vec > true_mean);

mean(BD_mean_cov);

# variance estimates;
mean(sqrt(BD_mean_var_vec));
sqrt(var(BD_mean_vec));


# complete case analysis;

CC_mean_bias=mean(CC_mean_vec)-true_mean;
CC_mean_bias;
CC_mean_bias/true_mean;

CC_mean_mse=mean((CC_mean_vec-true_mean)^2);
CC_mean_mse;

# lower and upper 95% CI;
CC_mean_low95_vec=CC_mean_vec-1.96*sqrt(CC_mean_var_vec);
CC_mean_up95_vec=CC_mean_vec+1.96*sqrt(CC_mean_var_vec);

CC_mean_length=mean(CC_mean_up95_vec-CC_mean_low95_vec);
CC_mean_length;

CC_mean_cov=(CC_mean_low95_vec < true_mean)*(CC_mean_up95_vec > true_mean);

mean(CC_mean_cov);

# variance estimates;
mean(sqrt(CC_mean_var_vec));
sqrt(var(CC_mean_vec));


# multiple imputation estimates;
mean_mi_mean=rep(NA, cycle_no);
mean_mi_var=rep(NA, cycle_no);
mean_mi_df=rep(NA, cycle_no);
mean_mi_f=rep(NA, cycle_no);



for (i in 1:cycle_no)
{
mean_summary=pool.scalar(MI_mean_mat[i,], MI_mean_var_mat[i,], n=rowobs);
mean_mi_mean[i]=mean_summary$qbar;
mean_mi_var[i]=mean_summary$t;
mean_mi_df[i]=mean_summary$df;
mean_mi_f[i]=mean_summary$f;

# mean_mi_df[i]=(mi_no-1)*(1+mean_summary$ubar/((1+1/mi_no)*mean_summary$b))^2;
# rm=(1+1/mi_no)*mean_summary$b/mean_summary$ubar;
# mean_mi_f[i]=(rm+2/(mean_mi_df[i]+3))/(rm+1);
}


# For marginal means;
# proper imputation;

mean_mi_bias=mean(mean_mi_mean)-true_mean;
mean_mi_bias;
mean_mi_bias/true_mean;

mean_mi_mse=mean((mean_mi_mean-true_mean)^2);
mean_mi_mse;

# coverage;
MI_mean_low95_vec=mean_mi_mean-qt(.975, mean_mi_df)*sqrt(mean_mi_var);
MI_mean_up95_vec=mean_mi_mean+qt(.975, mean_mi_df)*sqrt(mean_mi_var);

MI_mean_length=mean(MI_mean_up95_vec-MI_mean_low95_vec);

MI_mean_length;

MI_mean_coverage=(MI_mean_low95_vec < true_mean)*(MI_mean_up95_vec > true_mean);
mean(MI_mean_coverage);

# variance estimates;
mean(sqrt(mean_mi_var));
sqrt(var(mean_mi_mean));

