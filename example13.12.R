#Data analysis
library(sampleSelection)
library(R2WinBUGS);
library(mice);
library(norm);
options(digits=4);
rm(list=ls());


# read the data;

a1c=read.csv(file="C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter13_univariate_nonignorable\\program\\A1c.csv");

rowobs=nrow(a1c);

attach(a1c);

# right-skewed distribution;

hist(HgA1c, main="");
mean(HgA1c[Diabetes==1], na.rm=T);
mean(HgA1c[Diabetes==0], na.rm=T);


hist(log(HgA1c), main="");

# ma1c is the missingness indicator;

mean(ma1c);

sum(is.na(Age));
sum(is.na(Diabetes));
sum(is.na(BMI));
# BMI has two missing values;
sum(is.na(male));
sum(is.na(white));

# impute the two missing values in BMI;
BMI_mean=mean(BMI, na.rm=T);
BMI_sd=sd(BMI, na.rm=T);

set.seed(197789);
BMI[61]=rnorm(1, mean=BMI_mean, sd=BMI_sd);

set.seed(197552);
BMI[81]=rnorm(1, mean=BMI_mean, sd=BMI_sd);

# the association between ma1c and male;
table(ma1c, male);


# descriptive statistics of covariates;
mean(Age);
sd(Age);

mean(Diabetes);

mean(BMI);
sd(BMI);

mean(white);


# propensity score analysis;
propensity=glm(1-ma1c~Age+Diabetes+BMI+white, family=binomial(link="probit"))

summary(propensity);

z=rep(0,length(Age));
z[is.na(HgA1c)==FALSE]=1;

# log transformation of HgA1c;

HgA1c.o=log(HgA1c)

linear=lm(HgA1c.o ~ Age+Diabetes+BMI+white);
summary(linear);

linear_untran=lm(HgA1c ~ Age+Diabetes+BMI+white);
summary(linear_untran);

HgA1c.o[z==0]=0;

selection_fit=selection(z~Age+Diabetes+BMI+white,HgA1c.o~Age+BMI+Diabetes+white)


# output the dataset;

ln_hga1c=HgA1c.o;
ln_hga1c[ma1c==1]=NA;

# varname_cov=c("ma1c[]", "ln_hga1c[]", "Age[]", "Diabetes[]", "BMI[]", "white[]");

# write.table(cbind(ma1c, ln_hga1c, Age, Diabetes, BMI, white), file=paste("C:\\Users\\Guanghui He\\Personal\\Yulei\\research\\MIbook\\chapter13_univariate_nonignorable\\program\\a1c.dat",sep=""), row.name=FALSE, col.name=varname_cov);

# impute missing ln_hga1c using the selection model;
mi_no=50;
MI_mean_mat=MI_mean_var_mat=BUGS_mean_mat=BUGS_mean_var_mat=rep(NA, mi_no);

# rho is set as -0.9;
model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_ne9.txt")
# rho is set as -0.7;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_ne7.txt")
# rho is set as -0.5;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_ne5.txt")
# rho is set as -0.3;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_ne3.txt")
# rho is set as -0.1;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_ne1.txt")
# rho is set as 0.1;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_po1.txt")
# rho is set as 0.3;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_po3.txt")
# rho is set as 0.5;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_po5.txt")

# rho is set as 0.7;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_po7.txt")


# rho is set as 0.9;
# model.file <- system.file(package="R2WinBUGS", "model", "selectionmodel_a1c_po9.txt")






# Let's take a look:
file.show(model.file)

gibbs_no=20000;

process_stage2_data_simulate=list(n=rowobs, ma1c=ma1c, ln_hga1c=ln_hga1c, Age=Age, Diabetes=Diabetes, BMI=BMI, white=white);
process_stage2_init= function(){list(a = c(1.742, -0.003, 0.005, 0.253, -0.005), b = c(0.216, -0.005, 0.741, 0.016, 0.064), sigma=0.16)};
process_stage2_parameters=c("y_impute");

process_stage2.sim <- bugs(data=process_stage2_data_simulate, inits=process_stage2_init, parameters=process_stage2_parameters, model.file,
    n.chains=1, n.iter=gibbs_no, n.burnin=gibbs_no/2, n.thin=gibbs_no/(2*mi_no), bugs.seed=197789,
    bugs.directory="C:/Users/Guanghui He/Personal/Yulei/research/winbugs14_full_patched/WinBUGS14/", summary.only=FALSE, debug=FALSE,
    working.directory=NULL, clearWD=TRUE);

attach.bugs(process_stage2.sim);
rm(process_stage2.sim);

y_miss=ln_hga1c;

y_completed_BUGS=y_completed_MAR=y_miss;

y_completed_BUGS_mat=y_completed_MAR_mat=matrix(NA, nrow=rowobs, ncol=mi_no);



for (i in 1:mi_no)

{

y_completed_BUGS[ma1c==1]=y_impute[i,];

y_completed_BUGS_mat[,i]=y_completed_BUGS;

BUGS_mean_mat[i]=mean(exp(y_completed_BUGS));
BUGS_mean_var_mat[i]=var(exp(y_completed_BUGS))/rowobs;

}

hist(y_completed_BUGS, xlab="Completed log(HgA1c) by selection model, rho=-.9", main="");

# combine the estimates;
BUGS_summary=pool.scalar(BUGS_mean_mat, BUGS_mean_var_mat, n=rowobs);
BUGS_summary;
BUGS_summary$qbar
BUGS_summary$qbar-qt(.975, BUGS_summary$df)*sqrt(BUGS_summary$t);
BUGS_summary$qbar+qt(.975, BUGS_summary$df)*sqrt(BUGS_summary$t);
BUGS_summary$fmi;


# complete-case mean;
obs_no=rowobs-sum(ma1c);

CC_mean=mean(HgA1c, na.rm=T);

CC_mean-1.96*sd(HgA1c, na.rm=T)/sqrt(obs_no);
CC_mean+1.96*sd(HgA1c, na.rm=T)/sqrt(obs_no);



##############################################################################
function (rho,K=5) 
{
   library(sampleSelection)
   CC=CC.se=CC.cr=0
   MI=MI.sd=MI.se=MI.cr=0
   a1c=read.csv("A1c.csv")
   attach(a1c)
   n=length(HgA1c)
   y=HgA1c
   z=rep(0,n)
   z[is.na(y)==FALSE]=1
   y.o=y
   y.o[z==0]=0
   x=Diabetes
   mrate=(1-mean(z))
   print(mrate)
#CC
   CC=mean(y[z==1])
   CC.se=sqrt(var(y[z==1])/length(y[z==1]))
   CC.l=CC-1.96*CC.se
   CC.u=CC+1.96*CC.se

#MI
   dat.MI=matrix(nrow = n, ncol =K)
   dat.MI[,1:K]=y
   boot=matrix(sample(c(1:n), n*K, replace=T), nrow = n, ncol = K) 
   for(k in 1:K){
      boot.k=boot[,k]
      y.b=y[boot.k]
      z.b=z[boot.k]
      x.b=x[boot.k]   
      y.b1=y.b[z.b==1 & x.b==1]
      y.b0=y.b[z.b==1 & x.b==0]
      y.b1.s=sort(y.b1)
      y.b0.s=sort(y.b0)
      n1.o=length(y.b1)
      n0.o=length(y.b0)
      n1.s=ceiling(n1.o*(1-abs(rho)))
      n0.s=ceiling(n0.o*(1-abs(rho)))
      if(rho>=0){
         MI1.set=y.b1.s[1:n1.s]
         MI0.set=y.b0.s[1:n0.s]
      }
      if(rho<0){
         MI1.set=y.b1.s[(n1.o-n1.s+1):n1.o]
         MI0.set=y.b0.s[(n0.o-n0.s+1):n0.o]
      }
      for(i in 1:n){
          if(z[i]==0 & x[i]==1)
             dat.MI[i,k]=sample(MI1.set,1)
          if(z[i]==0 & x[i]==0)
             dat.MI[i,k]=sample(MI0.set,1)
      }
   }
   mi.est=rep(0,3)
   for(k in 1:K){
      mi.est[1]=mi.est[1]+mean(dat.MI[,k])/K
      mi.est[2]=mi.est[2]+(mean(dat.MI[,k]))^2
      mi.est[3]=mi.est[3]+var(dat.MI[,k])/n/K
   }
   MI=mi.est[1]
   B.mi=(mi.est[2]-K*mi.est[1]^2)/(K-1)
   MI.se=sqrt(mi.est[3]+((K+1)*B.mi)/K) 
   v.mi=(K-1)*(1+mi.est[3]/((K+1)*B.mi/K))^2
   t.crit=qt(0.975,v.mi)
   MI.l=mi.est[1]-t.crit*MI.se
   MI.u=mi.est[1]+t.crit*MI.se

   round(rbind(c(CC,CC.se,CC.l,CC.u),c(MI,MI.se,MI.l,MI.u)),3)
}


sen.MNAR=function (rep,n,rho,K=5) 
{
	library(mvtnorm)
	FO=FO.sd=FO.se=FO.cr=0
	CC=CC.sd=CC.se=CC.cr=0
	MI=MI.sd=MI.se=MI.cr=0
	mrate=0
	for(try in 1:rep){
		dat=rmvnorm(n,mean=c(0,0),sigma=matrix(c(1,rho,rho,1),ncol=2))
		y=dat[,1]
		x=dat[,2]
		z=rep(0,n)
		z[x>0]=1
		mrate=mrate+(1-mean(z))/rep
		#FO
		u.fo=mean(y)
		FO=FO+u.fo/rep
		FO.sd=FO.sd+u.fo^2
		fo.se=sqrt(var(y)/n)
		FO.se=FO.se+fo.se/rep
		fo.l=u.fo-1.96*fo.se
		fo.u=u.fo+1.96*fo.se
		if(fo.l<=0 & fo.u>=0)
		  FO.cr=FO.cr+1/rep
		#CC
		u.cc=mean(y[z==1])
		CC=CC+u.cc/rep
		CC.sd=CC.sd+u.cc^2
		cc.se=sqrt(var(y[z==1])/length(y[z==1]))
		CC.se=CC.se+cc.se/rep
		cc.l=u.cc-1.96*cc.se
		cc.u=u.cc+1.96*cc.se
		if(cc.l<=0 & cc.u>=0)
		  CC.cr=CC.cr+1/rep
		#MI
		dat.MI=matrix(nrow = n, ncol =K)
		dat.MI[,1:K]=y
		boot=matrix(sample(c(1:n), n*K, replace=T), nrow = n, ncol = K) 
		for(k in 1:K){
			boot.k=boot[,k]
			y.b=y[boot.k]
			z.b=z[boot.k]
			y.b.z=y.b[z.b==1]
			y.b.z.s=sort(y.b.z)
			n.o=length(z.b[z.b==1])
		    n.s=ceiling(n.o*(1-rho))
		    MI.set=y.b.z.s[1:n.s]
            for(i in 1:n){
            	if(z[i]==0)
                	dat.MI[i,k]=sample(MI.set,1)
            }			
		}
		mi.est=rep(0,3)
		for(k in 1:K){
			mi.est[1]=mi.est[1]+mean(dat.MI[,k])/K
			mi.est[2]=mi.est[2]+(mean(dat.MI[,k]))^2
			mi.est[3]=mi.est[3]+var(dat.MI[,k])/n/K
		}
		MI=MI+mi.est[1]/rep
		MI.sd=MI.sd+mi.est[1]^2
	    B.mi=(mi.est[2]-K*mi.est[1]^2)/(K-1)
        se.mi=sqrt(mi.est[3]+((K+1)*B.mi)/K)
        MI.se=MI.se+se.mi/rep
        v.mi=(K-1)*(1+mi.est[3]/((K+1)*B.mi/K))^2
        t.crit=qt(0.975,v.mi)
        l.mi=mi.est[1]-t.crit*se.mi
        u.mi=mi.est[1]+t.crit*se.mi
        if(l.mi<=0 & u.mi>=0)
           MI.cr=MI.cr+1/rep
	}
	FO.sd=sqrt((FO.sd-rep*FO^2)/(rep-1))
	CC.sd=sqrt((CC.sd-rep*CC^2)/(rep-1))
	MI.sd=sqrt((MI.sd-rep*MI^2)/(rep-1))
	print(round(c(rep,n,rho,K,mrate),3))
	round(rbind(c(FO,FO.sd,FO.se,FO.cr),c(CC,CC.sd,CC.se,CC.cr),c(MI,MI.sd,MI.se,MI.cr)),3)
}
