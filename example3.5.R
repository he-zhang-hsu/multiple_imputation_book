# example 3.5

library(mice);

# the mean estimate of 5 NHIS 2016 imputed family income;
mean_vec=c(66044, 66107, 65999, 66091, 66044);
se_vec=c(690.58, 701.29, 684.69, 690.83, 696.10);
n=40875;

mean_summary=pool.scalar(mean_vec, se_vec^2, n=n);

mean_summary;

sqrt(mean_summary$t);

mean_summary$qbar-qt(.975, mean_summary$df)*sqrt(mean_summary$t);
mean_summary$qbar+qt(.975, mean_summary$df)*sqrt(mean_summary$t);

lambda=1.2*mean_summary$b/mean_summary$t

lambda;

low95_vec=mi_mean-qt(.975, mi_df)*sqrt(mi_var);
up95_vec=mi_mean+qt(.975, mi_df)*sqrt(mi_var);
