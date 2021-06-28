# Example 8.5

library(survival)
library(readr)

# functions;
MITIND=function(dat, simulation, choice, N.MI=10, n.s=10, w1=0.8, w2=0.2)
{
    #This function imputes event times for censored observations and then derives KM estimates based on Hsu et al. (2006)
    #dat: 1st column: observed time; 2nd column: censoring indicator; 3rd-7th column: covariates
    #simulation: 0: data analysis; 1: simulation 
    #choice: choice for the covariates included in the two working models; 
    #       1: all 5 covariates included in both working models
    #       2: all 5 covariates included in the event time working model and the first 3 covariates included in the censoring time working model
    #       3: the first 3 covariates included in the event time working model and all 5 covariates included in the censoring time working model
    #       4: the first 3 covariates included in both working models
    #N.MI: # of imputes for each censored observation
    #n.s: size of the nearest neighborhood
    
    library(survival)
    library(TruncatedDistributions)
    library(truncdist)
    X1=as.matrix(dat[,3])
    X2=as.matrix(dat[,4])
    X3=as.matrix(dat[,5])
    X4=as.matrix(dat[,6])
    X5=as.matrix(dat[,7])
    T.o=as.matrix(dat[,1])
    I.d=as.matrix(dat[,2])
    crate=(1-mean(I.d))
    N=length(T.o)  
    
    #NNMI
    c1=N.MI
    c2=N.MI+1
    c3=2*N.MI
    #initialization of imputed data matrix 
    dat.NNMI = matrix(nrow = N, ncol = 2 * N.MI)
    dat.NNMI[, (1:c1)] = T.o
    dat.NNMI[, (c2:c3)] = I.d
    dat.PMI=dat.NNMI
    ID=seq(1,N,1)
    ID.c=ID[I.d==0]      #all IDs with censored observations
    N.c=length(ID.c)     #number of censored observations
    
    for(j in 1:N.MI) {
        boot.j=sample(ID,N,replace=T) #observations selected for bootstrap sample
        I.Bd=I.d[boot.j]
        T.Bo=T.o[boot.j]
        X1.B=X1[boot.j]
        X2.B=X2[boot.j]
        X3.B=X3[boot.j]
        X4.B=X4[boot.j]
        X5.B=X5[boot.j]
        #fit Weibull regression for PMI
        if(choice==1 | choice==2)
            fit.w=survreg(Surv(T.Bo, I.Bd) ~ X1.B + X2.B + X3.B + X4.B + X5.B, dist="weibull",scale=0)
        if(choice==3 | choice==4)
            fit.w=survreg(Surv(T.Bo, I.Bd) ~ X1.B + X2.B + X3.B, dist="weibull",scale=0)
        a.w=1/summary(fit.w)$scale
        
        #fit the two working models for each bootstrap sample
        if(choice==1){
            Xt=as.matrix(cbind(X1, X2, X3, X4,X5))
            Xc=as.matrix(cbind(X1, X2, X3, X4,X5))
            Xt.B=as.matrix(cbind(X1.B, X2.B, X3.B, X4.B,X5.B))
            Xc.B=as.matrix(cbind(X1.B, X2.B, X3.B, X4.B,X5.B))
            cox.Bt=coxph(Surv(T.Bo, I.Bd) ~ X1.B + X2.B + X3.B + X4.B + X5.B)
            cox.Bc=coxph(Surv(T.Bo, (1 - I.Bd)) ~ X1.B + X2.B + X3.B + X4.B + X5.B)
        }
        if(choice==2){
            Xt=as.matrix(cbind(X1, X2, X3, X4,X5))
            Xc=as.matrix(cbind(X1, X2, X3))
            Xt.B=as.matrix(cbind(X1.B, X2.B, X3.B, X4.B,X5.B))
            Xc.B=as.matrix(cbind(X1.B, X2.B, X3.B))
            cox.Bt=coxph(Surv(T.Bo, I.Bd) ~ X1.B + X2.B + X3.B + X4.B + X5.B)
            cox.Bc=coxph(Surv(T.Bo, (1 - I.Bd)) ~ X1.B + X2.B + X3.B)
        }
        if(choice==3){
            Xt=as.matrix(cbind(X1, X2, X3))
            Xc=as.matrix(cbind(X1, X2, X3, X4,X5))
            Xt.B=as.matrix(cbind(X1.B, X2.B, X3.B))
            Xc.B=as.matrix(cbind(X1.B, X2.B, X3.B, X4.B,X5.B))
            cox.Bt=coxph(Surv(T.Bo, I.Bd) ~ X1.B + X2.B + X3.B)
            cox.Bc=coxph(Surv(T.Bo, (1 - I.Bd)) ~ X1.B + X2.B + X3.B + X4.B + X5.B)
        }
        if(choice==4){
            Xt=as.matrix(cbind(X1, X2, X3))
            Xc=as.matrix(cbind(X1, X2, X3))
            Xt.B=as.matrix(cbind(X1.B, X2.B, X3.B))
            Xc.B=as.matrix(cbind(X1.B, X2.B, X3.B))
            cox.Bt=coxph(Surv(T.Bo, I.Bd) ~ X1.B + X2.B + X3.B)
            cox.Bc=coxph(Surv(T.Bo, (1 - I.Bd)) ~ X1.B + X2.B + X3.B)
        }
        coef.Bt=cox.Bt$coef
        coef.Bc=cox.Bc$coef
        
        #risk scores
        score.Bt=Xt.B %*% coef.Bt
        score.Bc=Xc.B %*% coef.Bc
        #standardized risk scores
        sscore.Bt=(score.Bt - mean(score.Bt))/sqrt(as.numeric(var(score.Bt)))
        sscore.Bc=(score.Bc - mean(score.Bc))/sqrt(as.numeric(var(score.Bc)))
        for(i in 1:N.c) {    #only impute for censored observations
            ID.i=ID.c[i]
            #PMI: imputation based on Weibull distribution
            lti.w=as.numeric(c(1,Xt[ID.i,  ]) %*% as.matrix(summary(fit.w)$coefficients))
            lambdai.w=exp(-lti.w*a.w)
            bi.w=(1/lambdai.w)^(1/a.w)
            r.w=rtweibull(1,shape=a.w,scale=bi.w,T.o[ID.i],Inf)
            dat.PMI[ID.i,j]=r.w
            dat.PMI[ID.i,N.MI+j]=1
            
            #NNMI: imputation based on Kaplan-Meier imputation in Hsu et al. (2006)
            #risk score for the ith censored observation of the original data
            score.ti=Xt[ID.i,  ] %*% coef.Bt
            score.ci=Xc[ID.i,  ] %*% coef.Bc
            #standardized risk score for the ith censored observation of the original data
            sscore.ti=(score.ti - mean(score.Bt))/sqrt(as.numeric(var(score.Bt)))
            sscore.ci=(score.ci - mean(score.Bc))/sqrt(as.numeric(var(score.Bc)))
            #distance between the ith censored observation and the observations in the bootstrap sample    
            D=(w1 * (sscore.Bt - as.numeric(sscore.ti))^2 + w2 * (sscore.Bc - as.numeric(sscore.ci))^2)^0.5
            T.Bso=T.Bo[T.Bo > T.o[ID.i]]
            I.Bsd=I.Bd[T.Bo > T.o[ID.i]]
            D.s=D[T.Bo > T.o[ID.i]]
            T.Bo.mi=T.Bso[order(rank(D.s))]
            I.Bd.mi=I.Bsd[order(rank(D.s))]
            n.mi= min(n.s, length(T.Bo.mi))
            if(n.mi > 0) {
                dat.mi=cbind(T.Bo.mi[1:n.mi], I.Bd.mi[1:n.mi])
                if(is.na(sum(dat.mi[, 1])) == F) {
                    BNNMI=survfit(Surv(dat.mi[, 1], dat.mi[, 2])~1)
                    F.BNNMI=1 - BNNMI$surv
                    u.im1 <- runif(1)
                    if(u.im1 <= max(F.BNNMI) & sum(dat.mi[, 2]) >0) {
                        dat.NNMI[ID.i, j]=min(BNNMI$time[F.BNNMI == min(F.BNNMI[F.BNNMI >= u.im1])])
                        dat.NNMI[ID.i, (N.MI + j)]= 1
                    }
                    if(u.im1 > max(F.BNNMI) & sum(dat.mi[, 2]) >0) {
                        dat.NNMI[ID.i, j]= max(BNNMI$time)
                        dat.NNMI[ID.i, (N.MI + j)]= 0
                    }
                    if(sum(dat.mi[, 2]) == 0) {
                        dat.NNMI[ID.i, j]= max(dat.mi[, 1])
                        dat.NNMI[ID.i, (N.MI + j)]= 0
                    }
                }
            }
        }
    }
    Tpmi.u=sort(unique(as.numeric(dat.PMI[,1:c1])))
    Npmi.u=length(Tpmi.u)
    pmi=matrix(0,nrow=Npmi.u,ncol=3)
    To.u=sort(unique(as.numeric(dat.NNMI[,1:c1])))
    No.u=length(To.u)
    nnmi=matrix(0,nrow=No.u,ncol=3)
    for(j in 1:c1) {
        KM.pmi=survfit(Surv(dat.PMI[, j], dat.PMI[, (c1 + j)])~1)
        for(i in 1:Npmi.u){
            if(Tpmi.u[i]<min(KM.pmi$time)){
                pmi[i,1]=pmi[i,1] + 1.0/c1
                pmi[i,2]=pmi[i,2] + 0/c1
                pmi[i,3]=pmi[i,3] + 1.0^2
            }
            if(Tpmi.u[i]>=min(KM.pmi$time)){
                ti.pmi=max(KM.pmi$time[KM.pmi$time <= Tpmi.u[i]])
                pmi[i,1]=pmi[i,1] + KM.pmi$surv[KM.pmi$time == ti.pmi]/c1
                pmi[i,2]=pmi[i,2] + (KM.pmi$surv[KM.pmi$time == ti.pmi] * KM.pmi$std.err[KM.pmi$time == ti.pmi])^2/c1
                pmi[i,3]=pmi[i,3] + (KM.pmi$surv[KM.pmi$time == ti.pmi])^2
            }
        }
        KM.nnmi=survfit(Surv(dat.NNMI[, j], dat.NNMI[, (c1 + j)])~1)
        for(i in 1:No.u){
            if(To.u[i]<min(KM.nnmi$time)){
                nnmi[i,1]=nnmi[i,1] + 1.0/c1
                nnmi[i,2]=nnmi[i,2] + 0/c1
                nnmi[i,3]=nnmi[i,3] + 1.0^2
            }
            if(To.u[i]>=min(KM.nnmi$time)){
                ti.nnmi=max(KM.nnmi$time[KM.nnmi$time <= To.u[i]])
                nnmi[i,1]=nnmi[i,1] + KM.nnmi$surv[KM.nnmi$time == ti.nnmi]/c1
                nnmi[i,2]=nnmi[i,2] + (KM.nnmi$surv[KM.nnmi$time == ti.nnmi] * KM.nnmi$std.err[KM.nnmi$time == ti.nnmi])^2/c1
                nnmi[i,3]=nnmi[i,3] + (KM.nnmi$surv[KM.nnmi$time == ti.nnmi])^2
            }
        }
    }
    #PMI
    PMI=matrix(0,nrow=Npmi.u,ncol=4)
    PMI[,1]=Tpmi.u                                  #time point
    for(i in 1:Npmi.u){
        PMI[i,2]=pmi[i,1]                           #S(t) estimate
        U.pmi=pmi[i,2]                              #within a MI dataset variation
        B.pmi=(pmi[i,3]-c1*pmi[i,1]^2)/(c1-1)       #between MI datasets variation
        var50.pmi=U.pmi + (c2 *B.pmi)/c1
        PMI[i,3]=sqrt(var50.pmi)                    #standard error 
        PMI[i,4]=(c1-1)*(1+(U.pmi*c1/c2/B.pmi))^2   #degrees of freedom for t distribution
    }
    
    #NNMI
    NNMI=matrix(0,nrow=No.u,ncol=4)
    NNMI[,1]=To.u                                    #time point
    for(i in 1:No.u){
        NNMI[i,2]=nnmi[i,1]                           #S(t) estimate
        U.nnmi=nnmi[i,2]                              #within a MI dataset variation
        B.nnmi=(nnmi[i,3]-c1*nnmi[i,1]^2)/(c1-1)       #between MI datasets variation
        var50.nnmi=U.nnmi + (c2 *B.nnmi)/c1
        NNMI[i,3]=sqrt(var50.nnmi)                    #standard error
        NNMI[i,4]=(c1-1)*(1+(U.nnmi*c1/c2/B.nnmi))^2   #degrees of freedom for t distribution
    }
    #time: time 
    #surv: KM S(t) estimates 
    #se: standard error 
    #df: degrees of freedom for t distribution 
    output=list("PMI.time"=PMI[,1],"PMI.surv"=PMI[,2],"PMI.se"=PMI[,3],"PMI.df"=PMI[,4],
                "NNMI.time"=NNMI[,1],"NNMI.surv"=NNMI[,2],"NNMI.se"=NNMI[,3],"NNMI.df"=NNMI[,4])
    return(output)
}



# read the data;
pca=read_csv("pca.csv")

#Working PH model for event time
PH.e=coxph(Surv(time,status)~age+log(psa)+gscore+stage+dose,data=pca)

#Working PH model for censoring time
PH.c=coxph(Surv(time,1-status)~age+log(psa)+gscore+stage+dose,data=pca)

#put the estimation results side by side for the two working models
# Table 8.6
print(cbind(summary(PH.e)$coefficients[,c(1,3,5)],summary(PH.c)$coefficients[,c(1,3,5)]))

#S(t) estimates at 5 and 10 years
#PO
KM.po=survfit(Surv(pca$time,pca$status)~1)
t5.po = max(KM.po$time[KM.po$time <= 5])
PO5 = round(KM.po$surv[KM.po$time == t5.po],4)
PO5.se = round(KM.po$surv[KM.po$time == t5.po] * KM.po$std.err[KM.po$time == t5.po],4)
t10.po = max(KM.po$time[KM.po$time <= 10])
PO10 = round(KM.po$surv[KM.po$time == t10.po],4)
PO10.se = round(KM.po$surv[KM.po$time == t10.po] * KM.po$std.err[KM.po$time == t10.po],4)
PO.est=c(PO5,PO5.se,PO10,PO10.se)

#MI
KM.mi=MITIND(pca,0,1)
#PMI
t5.pmi = max(KM.mi$PMI.time[KM.mi$PMI.time <= 5])
PMI5 = round(KM.mi$PMI.surv[KM.mi$PMI.time == t5.pmi],4)
PMI5.se = round(KM.mi$PMI.se[KM.mi$PMI.time == t5.pmi],4)
t10.pmi = max(KM.mi$PMI.time[KM.mi$PMI.time <= 10])
PMI10 = round(KM.mi$PMI.surv[KM.mi$PMI.time == t10.pmi],4)
PMI10.se = round(KM.mi$PMI.se[KM.mi$PMI.time == t10.pmi],4)
PMI.est=c(PMI5,PMI5.se,PMI10,PMI10.se)

#NNMI
t5.nnmi = max(KM.mi$NNMI.time[KM.mi$NNMI.time <= 5])
NNMI5 = round(KM.mi$NNMI.surv[KM.mi$NNMI.time == t5.nnmi],4)
NNMI5.se = round(KM.mi$NNMI.se[KM.mi$NNMI.time == t5.nnmi],4)
t10.nnmi = max(KM.mi$NNMI.time[KM.mi$NNMI.time <= 10])
NNMI10 = round(KM.mi$NNMI.surv[KM.mi$NNMI.time == t10.nnmi],4)
NNMI10.se = round(KM.mi$NNMI.se[KM.mi$NNMI.time == t10.nnmi],4)
NNMI.est=c(NNMI5,NNMI5.se,NNMI10,NNMI10.se)
St.est=rbind(PO.est,PMI.est,NNMI.est)
row.names(St.est)=c("PO","PMI","NNMI")
colnames(St.est)=c("S(t)","se(5)","S(10)","se(10)")

# Table 8.7
print(St.est)

#Kaplan-Meier curves
# Fig. 8.3
plot(KM.po,xlab="Time (years)",ylab="Recurrence-free probability",conf.int=FALSE, mark.time=FALSE,ylim=c(0.65,1),xlim=c(0,15),lty=1,lwd=2,col="blue")
lines(KM.mi$PMI.time,KM.mi$PMI.surv,type="s",lty=5,col="green",lwd=2)
lines(KM.mi$NNMI.time,KM.mi$NNMI.surv,type="s",lty=3,col="red",lwd=2)
legend(8,1,c("PO","PMI","NNMI"),lty=c(1,5,3),col=c("blue","green","red"),lwd=2)



