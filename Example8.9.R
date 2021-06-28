# Example 8.9

library(survival)
library(gtsummary)
library(readr)
library(NNMIS)
library(dplyr)

# read the data;
BCA=read_csv("SEERS4BCA.csv")
BCA$race=rep(NA,length(BCA[,1]))
BCA$race[BCA$white==1]="White"
BCA$race[BCA$black==1]="Black"
BCA$race[BCA$other==1]="Other"
BCA$her2=factor(BCA$her2,levels=c(0,1),labels=c("Negative","Positive"))
BCA$surgery=1-BCA$nosurgery
BCA$surgery=factor(BCA$surgery,levels=c(0,1),labels=c("No","Yes"))
BCA$rad=1-BCA$norad
BCA$rad=factor(BCA$rad,levels=c(0,1),labels=c("No","Yes"))

# Summary of patient characteristics 
# Table 8.11
BCA.s=BCA %>% select(age,race,her2,surgery,rad)
tbl_summary(
BCA.s,
label=list(age~"Age",race~"Race",her2~"HER2",surgery~"Surgery",rad~"Radiation"),
statistic=list(all_continuous()~"{mean}±{sd}"),
digits = list(all_categorical() ~ c(0, 2),all_continuous() ~ c(2, 2)),
missing = "ifany", # don't list missing data separately
missing_text=c("Missing")
) %>%
modify_header(label = "**Variable**") %>% # update the column header
bold_labels()


#PO
fit.cc=coxph(Surv(time,status)~her2+age+black+other+surgery+rad,data=BCA)

#impute missing covariate using NNMIS package
BCA.mi=NNMIS(BCA$her2,xa=cbind(BCA$age,BCA$black,BCA$other,BCA$surgery,BCA$rad),xb=cbind(BCA$age,BCA$black,BCA$other,BCA$surgery,BCA$rad),time=BCA$time,event=BCA$status)

#fit Cox regression based on the imputed datasets using coxph.pool in NNMIS
fit.nnmi=coxph.pool(BCA.mi, BCA$time, BCA$status,cbind(BCA$age,BCA$black,BCA$other,BCA$surgery,BCA$rad))

est=round(cbind(summary(fit.cc)$coefficients[,c(1,3,5)],fit.nnmi[,c(1,3,6)]),4)
row.names(est)=c("HER2","Age","Black","Other","Surgery","Radiation")
colnames(est)=c("coef.cc","se.cc","p.cc","coef.nnmi","se.nnmi","p.nnmi")
est

# Table 8.13 is on the hazards ratio



