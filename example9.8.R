library(R2WinBUGS);
library(MASS);
library(mice);
library(norm);
library(mix);
# library(HI);
options(digits=4);
rm(list=ls());
# winbugs imputation;

# apply missing cases to the data;
# data
# process_stage2_data_simulate=list(t = structure(.Data = 
#		    c(12,    1, 21, 25, 11, 26, 27, 30, 13, 12, 21, 20, 23, 25, 23, 29, 35, NA, 31, 36,  
#			 32, 27, 23, 12, 18, NA, NA, 38, 29, 30, NA, 32, NA, NA, NA, NA, 25, 30, 37, 27,  
#			 22, 26, NA, 28, 19, 15, 12, 35, 35, 10, 22, 18, NA, 12, NA, NA, 31, 24, 37, 29,  
#			 27, 18, 22, 13, 18, 29, 28, NA, 16, 22, 26, 19, NA, NA, 17, 28, 26, 12, 17, 26),
#			.Dim = c(4, 20)),
#		t.cen = structure(.Data = 
#		    c(  0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 40,  0,  0,  
#			    0,   0,  0,  0,  0, 40, 40,  0,  0,  0, 40,  0, 40, 40, 40, 40,  0,  0,  0,  0,  
#			    0,  0, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0, 24,  0, 40, 40,  0,  0,  0,  0,  
#			    0,  0,  0,  0,  0,  0,  0, 20,  0,  0,  0,  0, 29, 10,  0,  0,  0,  0,  0,  0),
#			.Dim = c(4, 20)),
#		M = 4, N = 20);


# t = structure(.Data = 
#	    c(12,    1, 21, 25, 11, 26, 27, 30, 13, 12, 21, 20, 23, 25, 23, 29, 35, NA, 31, 36,  
#			 32, 27, 23, 12, 18, NA, NA, 38, 29, 30, NA, 32, NA, NA, NA, NA, 25, 30, 37, 27,  
#			 22, 26, NA, 28, 19, 15, 12, 35, 35, 10, 22, 18, NA, 12, NA, NA, 31, 24, 37, 29,  
#			 27, 18, 22, 13, 18, 29, 28, NA, 16, 22, 26, 19, NA, NA, 17, 28, 26, 12, 17, 26),
#			.Dim = c(4, 20));

#t.cen = structure(.Data = 
#		    c(  0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 40,  0,  0,  
#			    0,   0,  0,  0,  0, 40, 40,  0,  0,  0, 40,  0, 40, 40, 40, 40,  0,  0,  0,  0,  
#			    0,  0, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0, 24,  0, 40, 40,  0,  0,  0,  0,  
#			    0,  0,  0,  0,  0,  0,  0, 20,  0,  0,  0,  0, 29, 10,  0,  0,  0,  0,  0,  0),
#			.Dim = c(4, 20));

model0.file <- system.file(package="R2WinBUGS", "model", "linearmixed_imputeforR.txt")
# Let's take a look:
file.show(model0.file)



rowobs=N=30;
T=5;
x=c(8.0, 15.0, 22.0, 29.0, 36.0);
Y=matrix(c(151, 199, 246, 283, 320,145, 199, 249, 293, 354,147, 214, 263, 312, 328,155, 200, 237, 272, 297,135, 188, 230, 280, 323,159, 210, 252, 298, NA,141, 189, 231, 275, NA,159, 201, 248, 297, NA,177, 236, 285, 350, NA,134, 182, 220, 260, NA,160, 208, 261, NA, NA,143, 188, 220, NA, NA,154, 200, 244, NA, NA,171, 221, 270, NA, NA,163, 216, 242, NA, NA,160, 207, 248, NA, NA,142, 187, 234, NA, NA,156, 203, 243, NA, NA,157, 212, 259, NA, NA,152, 203, 246, NA, NA,154, 205, NA, NA, NA,139, 190, NA, NA, NA,146, 191, NA, NA, NA,157, 211, NA, NA, NA,132, 185, NA, NA, NA,160, NA, NA, NA, NA,169, NA, NA, NA, NA,157, NA, NA, NA, NA,137, NA, NA, NA, NA,153, NA, NA, NA, NA),nrow=30, ncol=5, byrow=T);
mean=c(0,0);
Omega=matrix(c(200,0,0,0.2), nrow=2, ncol=2);
prec = matrix(c(10^(-6),0,0,10^(-6)), nrow=2, ncol=2);

# plot the data;

matplot(x,t(Y), pch=19, xlab="day", ylab="weight", type="b");


process_stage2_init=function(){list(mu.beta=c(0,0), tauC=1)};

# process_stage2_parameters=c("mu.beta", "sigma", "Y_impute");
# process_stage2_parameters=c("mu.beta", "sigma", "Y");
process_stage2_parameters=c("Y", "mu.beta", "sigma", "R");


process_stage2_data_simulate=list("Y","x", "T", "N", "Omega", "mean", "prec");


# process_stage2_data_simulate=list("M","y2_miss", "x1_miss", "x2_miss");
# process_stage2_init= function(){list(beta0=-3, beta1=1.5, beta2=1.5, mu=0, alpha0=0, alpha1=0, tau=1, psi=1)};
# process_stage2_parameters=c("mu", "beta0", "beta1", "beta2", "mu", "alpha0", "alpha1", "tau", "psi", "y2", "x1", "x2");


gibbs_no=10000;


process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file=model0.file,
    n.chains=1, n.thin=100, n.iter=gibbs_no, n.burnin=gibbs_no/2, bugs.seed=197789, bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);

# fetch the object;

# the following program is in example 11.9
names(process_stage2.sim);
dim(process_stage2.sim$sims.array);



# fetch the imputation;

Y_6_5_impute=process_stage2.sim$sims.array[,,"Y[6,5]"];
Y_7_5_impute=process_stage2.sim$sims.array[,,"Y[7,5]"];
Y_8_5_impute=process_stage2.sim$sims.array[,,"Y[8,5]"];
Y_9_5_impute=process_stage2.sim$sims.array[,,"Y[9,5]"];
Y_10_5_impute=process_stage2.sim$sims.array[,,"Y[10,5]"];
Y_11_5_impute=process_stage2.sim$sims.array[,,"Y[11,5]"];
Y_12_5_impute=process_stage2.sim$sims.array[,,"Y[12,5]"];
Y_13_5_impute=process_stage2.sim$sims.array[,,"Y[13,5]"];
Y_14_5_impute=process_stage2.sim$sims.array[,,"Y[14,5]"];
Y_15_5_impute=process_stage2.sim$sims.array[,,"Y[15,5]"];
Y_16_5_impute=process_stage2.sim$sims.array[,,"Y[16,5]"];
Y_17_5_impute=process_stage2.sim$sims.array[,,"Y[17,5]"];
Y_18_5_impute=process_stage2.sim$sims.array[,,"Y[18,5]"];
Y_19_5_impute=process_stage2.sim$sims.array[,,"Y[19,5]"];
Y_20_5_impute=process_stage2.sim$sims.array[,,"Y[20,5]"];
Y_21_5_impute=process_stage2.sim$sims.array[,,"Y[21,5]"];
Y_22_5_impute=process_stage2.sim$sims.array[,,"Y[22,5]"];
Y_23_5_impute=process_stage2.sim$sims.array[,,"Y[23,5]"];
Y_24_5_impute=process_stage2.sim$sims.array[,,"Y[24,5]"];
Y_25_5_impute=process_stage2.sim$sims.array[,,"Y[25,5]"];
Y_26_5_impute=process_stage2.sim$sims.array[,,"Y[26,5]"];
Y_27_5_impute=process_stage2.sim$sims.array[,,"Y[27,5]"];
Y_28_5_impute=process_stage2.sim$sims.array[,,"Y[28,5]"];
Y_29_5_impute=process_stage2.sim$sims.array[,,"Y[29,5]"];
Y_30_5_impute=process_stage2.sim$sims.array[,,"Y[30,5]"];

Y_11_4_impute=process_stage2.sim$sims.array[,,"Y[11,4]"];
Y_12_4_impute=process_stage2.sim$sims.array[,,"Y[12,4]"];
Y_13_4_impute=process_stage2.sim$sims.array[,,"Y[13,4]"];
Y_14_4_impute=process_stage2.sim$sims.array[,,"Y[14,4]"];
Y_15_4_impute=process_stage2.sim$sims.array[,,"Y[15,4]"];
Y_16_4_impute=process_stage2.sim$sims.array[,,"Y[16,4]"];
Y_17_4_impute=process_stage2.sim$sims.array[,,"Y[17,4]"];
Y_18_4_impute=process_stage2.sim$sims.array[,,"Y[18,4]"];
Y_19_4_impute=process_stage2.sim$sims.array[,,"Y[19,4]"];
Y_20_4_impute=process_stage2.sim$sims.array[,,"Y[20,4]"];
Y_21_4_impute=process_stage2.sim$sims.array[,,"Y[21,4]"];
Y_22_4_impute=process_stage2.sim$sims.array[,,"Y[22,4]"];
Y_23_4_impute=process_stage2.sim$sims.array[,,"Y[23,4]"];
Y_24_4_impute=process_stage2.sim$sims.array[,,"Y[24,4]"];
Y_25_4_impute=process_stage2.sim$sims.array[,,"Y[25,4]"];
Y_26_4_impute=process_stage2.sim$sims.array[,,"Y[26,4]"];
Y_27_4_impute=process_stage2.sim$sims.array[,,"Y[27,4]"];
Y_28_4_impute=process_stage2.sim$sims.array[,,"Y[28,4]"];
Y_29_4_impute=process_stage2.sim$sims.array[,,"Y[29,4]"];
Y_30_4_impute=process_stage2.sim$sims.array[,,"Y[30,4]"];


Y_21_3_impute=process_stage2.sim$sims.array[,,"Y[21,3]"];
Y_22_3_impute=process_stage2.sim$sims.array[,,"Y[22,3]"];
Y_23_3_impute=process_stage2.sim$sims.array[,,"Y[23,3]"];
Y_24_3_impute=process_stage2.sim$sims.array[,,"Y[24,3]"];
Y_25_3_impute=process_stage2.sim$sims.array[,,"Y[25,3]"];
Y_26_3_impute=process_stage2.sim$sims.array[,,"Y[26,3]"];
Y_27_3_impute=process_stage2.sim$sims.array[,,"Y[27,3]"];
Y_28_3_impute=process_stage2.sim$sims.array[,,"Y[28,3]"];
Y_29_3_impute=process_stage2.sim$sims.array[,,"Y[29,3]"];
Y_30_3_impute=process_stage2.sim$sims.array[,,"Y[30,3]"];


Y_26_2_impute=process_stage2.sim$sims.array[,,"Y[26,2]"];
Y_27_2_impute=process_stage2.sim$sims.array[,,"Y[27,2]"];
Y_28_2_impute=process_stage2.sim$sims.array[,,"Y[28,2]"];
Y_29_2_impute=process_stage2.sim$sims.array[,,"Y[29,2]"];
Y_30_2_impute=process_stage2.sim$sims.array[,,"Y[30,2]"];

summary(Y_26_2_impute);
summary(Y_26_3_impute);
summary(Y_26_4_impute);
summary(Y_26_5_impute);

# obtain the model estimates;
#attach.bugs(process_stage2.sim);
#colMeans(mu.beta);
#var(mu.beta);
#mean(sigma^2);
#var(sigma^2);
# Omega.beta=apply(R,1,solve);
#rowMeans(Omega.beta);


# store estimates;
mi_no=50;


miximpute2_mean_vec=miximpute2_var_vec=miximpute3_mean_vec=miximpute3_var_vec=miximpute4_mean_vec=miximpute4_var_vec=miximpute5_mean_vec=miximpute5_var_vec=rep(NA, mi_no);
nmimpute2_mean_vec=nmimpute2_var_vec=nmimpute3_mean_vec=nmimpute3_var_vec=nmimpute4_mean_vec=nmimpute4_var_vec=nmimpute5_mean_vec=nmimpute5_var_vec=rep(NA, mi_no);

rowobs=nrow(Y);
obs_wave2=rowobs-sum(is.na(Y[,2]));
obs_wave3=rowobs-sum(is.na(Y[,3]));
obs_wave4=rowobs-sum(is.na(Y[,4]));
obs_wave5=rowobs-sum(is.na(Y[,5]));


# complete-case analysis;
 
cc2_mean=mean(Y[,2], na.rm="T");
cc2_var=var(Y[,2], na.rm="T")/obs_wave2;
cc3_mean=mean(Y[,3], na.rm="T");
cc3_var=var(Y[,3], na.rm="T")/obs_wave3;
cc4_mean=mean(Y[,4], na.rm="T");
cc4_var=var(Y[,4], na.rm="T")/obs_wave4;
cc5_mean=mean(Y[,5], na.rm="T");
cc5_var=var(Y[,5], na.rm="T")/obs_wave5;



# construct multiple imputation estimates for the means of waves 2-5;
y_miss=y_completed=Y;

for (i in 1:mi_no)
{
y_completed[6,5]=Y_6_5_impute[i];
y_completed[7,5]=Y_7_5_impute[i];
y_completed[8,5]=Y_8_5_impute[i];
y_completed[9,5]=Y_9_5_impute[i];
y_completed[10,5]=Y_10_5_impute[i];
y_completed[11,5]=Y_11_5_impute[i];
y_completed[12,5]=Y_12_5_impute[i];
y_completed[13,5]=Y_13_5_impute[i];
y_completed[14,5]=Y_14_5_impute[i];
y_completed[15,5]=Y_15_5_impute[i];
y_completed[16,5]=Y_16_5_impute[i];
y_completed[17,5]=Y_17_5_impute[i];
y_completed[18,5]=Y_18_5_impute[i];
y_completed[19,5]=Y_19_5_impute[i];
y_completed[20,5]=Y_20_5_impute[i];
y_completed[21,5]=Y_21_5_impute[i];
y_completed[22,5]=Y_22_5_impute[i];
y_completed[23,5]=Y_23_5_impute[i];
y_completed[24,5]=Y_24_5_impute[i];
y_completed[25,5]=Y_25_5_impute[i];
y_completed[26,5]=Y_26_5_impute[i];
y_completed[27,5]=Y_27_5_impute[i];
y_completed[28,5]=Y_28_5_impute[i];
y_completed[29,5]=Y_29_5_impute[i];
y_completed[30,5]=Y_30_5_impute[i];

y_completed[11,4]=Y_11_4_impute[i];
y_completed[12,4]=Y_12_4_impute[i];
y_completed[13,4]=Y_13_4_impute[i];
y_completed[14,4]=Y_14_4_impute[i];
y_completed[15,4]=Y_15_4_impute[i];
y_completed[16,4]=Y_16_4_impute[i];
y_completed[17,4]=Y_17_4_impute[i];
y_completed[18,4]=Y_18_4_impute[i];
y_completed[19,4]=Y_19_4_impute[i];
y_completed[20,4]=Y_20_4_impute[i];
y_completed[21,4]=Y_21_4_impute[i];
y_completed[22,4]=Y_22_4_impute[i];
y_completed[23,4]=Y_23_4_impute[i];
y_completed[24,4]=Y_24_4_impute[i];
y_completed[25,4]=Y_25_4_impute[i];
y_completed[26,4]=Y_26_4_impute[i];
y_completed[27,4]=Y_27_4_impute[i];
y_completed[28,4]=Y_28_4_impute[i];
y_completed[29,4]=Y_29_4_impute[i];
y_completed[30,4]=Y_30_4_impute[i];

y_completed[21,3]=Y_21_3_impute[i];
y_completed[22,3]=Y_22_3_impute[i];
y_completed[23,3]=Y_23_3_impute[i];
y_completed[24,3]=Y_24_3_impute[i];
y_completed[25,3]=Y_25_3_impute[i];
y_completed[26,3]=Y_26_3_impute[i];
y_completed[27,3]=Y_27_3_impute[i];
y_completed[28,3]=Y_28_3_impute[i];
y_completed[29,3]=Y_29_3_impute[i];
y_completed[30,3]=Y_30_3_impute[i];

y_completed[26,2]=Y_26_2_impute[i];
y_completed[27,2]=Y_27_2_impute[i];
y_completed[28,2]=Y_28_2_impute[i];
y_completed[29,2]=Y_29_2_impute[i];
y_completed[30,2]=Y_30_2_impute[i];


# multiple imputation estimates;
miximpute2_mean_vec[i]=mean(y_completed[,2]);
miximpute2_var_vec[i]=var(y_completed[,2])/rowobs;
miximpute3_mean_vec[i]=mean(y_completed[,3]);
miximpute3_var_vec[i]=var(y_completed[,3])/rowobs;
miximpute4_mean_vec[i]=mean(y_completed[,4]);
miximpute4_var_vec[i]=var(y_completed[,4])/rowobs;
miximpute5_mean_vec[i]=mean(y_completed[,5]);
miximpute5_var_vec[i]=var(y_completed[,5])/rowobs;

}


# multivariate normal imputation;
# y_miss_data=cbind(y_miss[,1], y_miss[,2], y_miss[,3], y_miss[,4], y_miss[,5]);

# s=prelim.norm(y_miss_data);

# thetahat=em.norm(s);

# rngseed(1);

# theta=da.norm(s,thetahat,steps=200, showits=T);

for (i in 1:mi_no)
{
# rngseed(i);



# imputation;
y_completed_norm=mice(y_miss, seed=i, m=1, maxit=100);
y_completed_norm_final=complete(y_completed_norm);



# completed data estimates;

# multiple imputation estimates;
nmimpute2_mean_vec[i]=mean(y_completed_norm_final[,2]);
nmimpute2_var_vec[i]=var(y_completed_norm_final[,2])/rowobs;
nmimpute3_mean_vec[i]=mean(y_completed_norm_final[,3]);
nmimpute3_var_vec[i]=var(y_completed_norm_final[,3])/rowobs;
nmimpute4_mean_vec[i]=mean(y_completed_norm_final[,4]);
nmimpute4_var_vec[i]=var(y_completed_norm_final[,4])/rowobs;
nmimpute5_mean_vec[i]=mean(y_completed_norm_final[,5]);
nmimpute5_var_vec[i]=var(y_completed_norm_final[,5])/rowobs;



}

# combine the multiple imputation estimates;
miximpute2_summary=pool.scalar(miximpute2_mean_vec, miximpute2_var_vec, n=rowobs);
miximpute2_mi_mean=miximpute2_summary$qbar;
miximpute2_mi_var=miximpute2_summary$t;
miximpute2_mi_df=miximpute2_summary$df;
miximpute2_mi_f=miximpute2_summary$f;

miximpute2_mi_mean;
sqrt(miximpute2_mi_var);

nmimpute2_summary=pool.scalar(nmimpute2_mean_vec, nmimpute2_var_vec);
nmimpute2_mi_mean=nmimpute2_summary$qbar;
nmimpute2_mi_var=nmimpute2_summary$t;
nmimpute2_mi_df=nmimpute2_summary$df;
nmimpute2_mi_f=nmimpute2_summary$f;

nmimpute2_mi_mean;
sqrt(nmimpute2_mi_var);


miximpute3_summary=pool.scalar(miximpute3_mean_vec, miximpute3_var_vec);
miximpute3_mi_mean=miximpute3_summary$qbar;
miximpute3_mi_var=miximpute3_summary$t;
miximpute3_mi_df=miximpute3_summary$df;
miximpute3_mi_f=miximpute3_summary$f;

miximpute3_mi_mean;
sqrt(miximpute3_mi_var);

nmimpute3_summary=pool.scalar(nmimpute3_mean_vec, nmimpute3_var_vec);
nmimpute3_mi_mean=nmimpute3_summary$qbar;
nmimpute3_mi_var=nmimpute3_summary$t;
nmimpute3_mi_df=nmimpute3_summary$df;
nmimpute3_mi_f=nmimpute3_summary$f;

nmimpute3_mi_mean;
sqrt(nmimpute3_mi_var);

miximpute4_summary=pool.scalar(miximpute4_mean_vec, miximpute4_var_vec);
miximpute4_mi_mean=miximpute4_summary$qbar;
miximpute4_mi_var=miximpute4_summary$t;
miximpute4_mi_df=miximpute4_summary$df;
miximpute4_mi_f=miximpute4_summary$f;

miximpute4_mi_mean;
sqrt(miximpute4_mi_var);

nmimpute4_summary=pool.scalar(nmimpute4_mean_vec, nmimpute4_var_vec);
nmimpute4_mi_mean=nmimpute4_summary$qbar;
nmimpute4_mi_var=nmimpute4_summary$t;
nmimpute4_mi_df=nmimpute4_summary$df;
nmimpute4_mi_f=nmimpute4_summary$f;

nmimpute4_mi_mean;
sqrt(nmimpute4_mi_var);


miximpute5_summary=pool.scalar(miximpute5_mean_vec, miximpute5_var_vec);
miximpute5_mi_mean=miximpute5_summary$qbar;
miximpute5_mi_var=miximpute5_summary$t;
miximpute5_mi_df=miximpute5_summary$df;
miximpute5_mi_f=miximpute5_summary$f;

miximpute5_mi_mean;
sqrt(miximpute5_mi_var);

nmimpute5_summary=pool.scalar(nmimpute5_mean_vec, nmimpute5_var_vec);
nmimpute5_mi_mean=nmimpute5_summary$qbar;
nmimpute5_mi_var=nmimpute5_summary$t;
nmimpute5_mi_df=nmimpute5_summary$df;
nmimpute5_mi_f=nmimpute5_summary$f;

nmimpute5_mi_mean;
sqrt(nmimpute5_mi_var);


cc2_mean;
sqrt(cc2_var);
cc3_mean;
sqrt(cc3_var);
cc4_mean;
sqrt(cc4_var);
cc5_mean;
sqrt(cc5_var);


# plot(seq(1,50,1), Y_6_5_impute);



# attach.bugs(process_stage2.sim);

process_stage2.sim;

#process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file="C:/Users/WDQ7/reliability_test_forR.txt",
#    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=1, bugs.seed=197789,
#    bugs.directory="C:/Program Files (x86)/OpenBUGS/OpenBUGS321", summary.only=FALSE, debug=FALSE,
#    working.directory=NULL, clearWD=TRUE);


# For logistic random-effects model

Y_complete=c(151, 199, 246, 283, 320, 145, 199, 249, 293, 354, 147, 214, 263, 312, 328, 155, 200, 237, 272, 297, 135, 188, 230, 280, 323,159, 210, 252, 298, 331,141, 189, 231, 275, 305,159, 201, 248, 297, 338,177, 236, 285, 350, 376,
							 134, 182, 220, 260, 296,
							 160, 208, 261, 313, 352,
							 143, 188, 220, 273, 314,
							 154, 200, 244, 289, 325,
							 171, 221, 270, 326, 358,
							 163, 216, 242, 281, 312,
							 160, 207, 248, 288, 324,
							 142, 187, 234, 280, 316,
							 156, 203, 243, 283, 317,
							 157, 212, 259, 307, 336,
							 152, 203, 246, 286, 321,
							 154, 205, 253, 298, 334,
							 139, 190, 225, 267, 302,
							 146, 191, 229, 272, 302,
							 157, 211, 250, 285, 323,
							 132, 185, 237, 286, 331,
							 160, 207, 257, 303, 345,
							 169, 216, 261, 295, 333,
							 157, 205, 248, 289, 316,
							 137, 180, 219, 258, 291,
							 153, 200, 244, 286, 324);

summary(Y_complete);

Y_binary=1*(Y_complete>=250);

Y_binary=1*(Y_complete>=300);

Y_binary=1*(Y_complete>=200);

