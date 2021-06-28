# Example 8.4

# install.packages("TruncatedDistributions", repos="http://R-Forge.R-project.org")

library(survival);
library(TruncatedDistributions)
# library(truncdist);

# functions;

dataind=function (N = 100, b0 = 4, b1 = -2, b2 = 0.5, b3 = -2, b4 = 2, b5 = 2, r0 = 3, r1 = -3, r2 = 0.5, r3 = -2, r4 = 1.5, r5 = 2)
{
    #N: sample size
    #b0-b5: coefficents for the hazard function of event time, h(t)=t^b0*exp(b1*X1+b2*X2+b3*X3+b4*X4+b5*X5)
    #r0-r5: coefficients for the hazard function of censoring time, h(c)=c^r0*exp(r1*X1+r2*X2+r3*X3+r4*X4+r5*X5)
    
    u.t=runif(N)
    u.c=runif(N)
    X1=runif(N)
    X2=runif(N)
    X3=runif(N)
    X4=runif(N)
    X5=runif(N)
    #generate event times based on h(t)=t^b0*exp(b1*X1+b2*X2+b3*X3+b4*X4+b5*X5)
    T.e=(-logb(u.t)*(b0+1)/exp(b1 * X1 + b2 * X2 + b3 * X3 + b4 * X4 + b5 * X5))^(1/(b0+1))
    #generate censoring times based on h(c)=c^r0*exp(c1*X1+c2*X2+c3*X3+c4*X4+c5*X5)
    T.c=(- logb(u.c)*(r0+1)/exp(r1 * X1 + r2 * X2 + r3 * X3 + r4 * X4 + r5 * X5))^(1/(r0+1))
    T.o=apply(cbind(T.e,T.c),1,min)
    I.d=as.numeric(T.o==T.e)
    crate=1-mean(I.d)
    data=cbind(T.o,I.d,X1,X2,X3,X4,X5,T.e)
    return(list("data"=data,"t50"=median(T.e),"crate"=crate))
}

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


simulation=function(rep=500,N=200,b0 = 4, b1 = -2, b2 = 0.5, b3 = -2, b4 = 2, b5 = 2, r0 = 3, r1 = -3, r2 = 0.5, r3 = -2, r4 = 1.5, r5 = 2,choice=1,N.MI=10,n.s=10,w1=0.8,w2=0.2)
{
    library(survival)
    t50=dataind(50000)$t50
    FO=FO.ss=FO.se=FO.cr=0
    PO=PO.ss=PO.se=PO.cr=0
    PMI=PMI.ss=PMI.se=PMI.cr=0
    NNMI=NNMI.ss=NNMI.se=NNMI.cr=0
    crate=0
    for(try in 1:rep){
        print(try)
        set.seed(try)
        sample=dataind(N,b0,b1,b2,b3,b4,b5,r0,r1,r2,r3,r4,r5)
        crate=crate+sample$crate/rep
        data=sample$data
        
        #FO
        KM.fo=survfit(Surv(data[,8],rep(1,N))~1)
        t50.fo=max(KM.fo$time[KM.fo$time <= t50])
        FO=FO+KM.fo$surv[KM.fo$time == t50.fo]/rep
        FO.se=FO.se+(KM.fo$surv[KM.fo$time == t50.fo] * KM.fo$std.err[KM.fo$time == t50.fo])^2/rep
        FO.ss=FO.ss+(KM.fo$surv[KM.fo$time == t50.fo])^2
        low.fo=KM.fo$surv[KM.fo$time == t50.fo]-1.96*KM.fo$surv[KM.fo$time == t50.fo] * KM.fo$std.err[KM.fo$time == t50.fo]
        high.fo=KM.fo$surv[KM.fo$time == t50.fo]+1.96*KM.fo$surv[KM.fo$time == t50.fo] * KM.fo$std.err[KM.fo$time == t50.fo]
        if(low.fo<=0.50 & high.fo>=0.50)
            FO.cr=FO.cr+100/rep
        #PO
        KM.po=survfit(Surv(data[,1],data[,2])~1)
        t50.po=max(KM.po$time[KM.po$time <= t50])
        PO=PO+KM.po$surv[KM.po$time == t50.po]/rep
        PO.se=PO.se+(KM.po$surv[KM.po$time == t50.po] * KM.po$std.err[KM.po$time == t50.po])^2/rep
        PO.ss=PO.ss+(KM.po$surv[KM.po$time == t50.po])^2
        low.po=KM.po$surv[KM.po$time == t50.po]-1.96*KM.po$surv[KM.po$time == t50.po] * KM.po$std.err[KM.po$time == t50.po]
        high.po=KM.po$surv[KM.po$time == t50.po]+1.96*KM.po$surv[KM.po$time == t50.po] * KM.po$std.err[KM.po$time == t50.po]
        if(low.po<=0.50 & high.po>=0.50)
            PO.cr=PO.cr+100/rep
        #MI
        KM.mi=MITIND(data,simulation=1,choice)
        #print(KM.mi)
        #PMI
        t50.pmi = max(KM.mi$PMI.time[KM.mi$PMI.time <= t50])
        PMI =PMI+KM.mi$PMI.surv[KM.mi$PMI.time == t50.pmi]/rep
        PMI.se =PMI.se+KM.mi$PMI.se[KM.mi$PMI.time == t50.pmi]/rep
        PMI.ss=PMI.ss+(KM.mi$PMI.surv[KM.mi$PMI.time == t50.pmi])^2
        low.pmi=KM.mi$PMI.surv[KM.mi$PMI.time == t50.pmi]-qt(0.975,KM.mi$PMI.df[KM.mi$PMI.time == t50.pmi])*KM.mi$PMI.se[KM.mi$PMI.time == t50.pmi]
        high.pmi=KM.mi$PMI.surv[KM.mi$PMI.time == t50.pmi]+qt(0.975,KM.mi$PMI.df[KM.mi$PMI.time == t50.pmi])*KM.mi$PMI.se[KM.mi$PMI.time == t50.pmi]
        if(low.pmi<=0.50 & high.pmi>=0.50)
            PMI.cr=PMI.cr+100/rep
        #NNMI
        t50.nnmi = max(KM.mi$NNMI.time[KM.mi$NNMI.time <= t50])
        NNMI =NNMI+KM.mi$NNMI.surv[KM.mi$NNMI.time == t50.nnmi]/rep
        NNMI.se =NNMI.se+KM.mi$NNMI.se[KM.mi$NNMI.time == t50.nnmi]/rep
        NNMI.ss=NNMI.ss+(KM.mi$NNMI.surv[KM.mi$NNMI.time == t50.nnmi])^2
        low.nnmi=KM.mi$NNMI.surv[KM.mi$NNMI.time == t50.nnmi]-qt(0.975,KM.mi$NNMI.df[KM.mi$NNMI.time == t50.nnmi])*KM.mi$NNMI.se[KM.mi$NNMI.time == t50.nnmi]
        high.nnmi=KM.mi$NNMI.surv[KM.mi$NNMI.time == t50.nnmi]+qt(0.975,KM.mi$NNMI.df[KM.mi$NNMI.time == t50.nnmi])*KM.mi$NNMI.se[KM.mi$NNMI.time == t50.nnmi]
        if(low.nnmi<=0.50 & high.nnmi>=0.50)
            NNMI.cr=NNMI.cr+100/rep
    }
    FO.sd=sqrt((FO.ss-rep*FO^2)/(rep-1))
    PO.sd=sqrt((PO.ss-rep*PO^2)/(rep-1))
    PMI.sd=sqrt((PMI.ss-rep*PMI^2)/(rep-1))
    NNMI.sd=sqrt((NNMI.ss-rep*NNMI^2)/(rep-1))
    output=matrix(nrow=4,ncol=5)
    output[1,]=c(round(FO,3),round(100*(FO-FO)/FO,1),round(FO.sd,3),round(FO.se,3),round(FO.cr,1))
    output[2,]=c(round(PO,3),round(100*(PO-FO)/FO,1),round(PO.sd,3),round(PO.se,3),round(PO.cr,1))
    output[3,]=c(round(PMI,3),round(100*(PMI-FO)/FO,1),round(PMI.sd,3),round(PMI.se,3),round(PMI.cr,1))
    output[4,]=c(round(NNMI,3),round(100*(NNMI-FO)/FO,1),round(NNMI.sd,3),round(NNMI.se,3),round(NNMI.cr,1))
    colnames(output)=c("Estimate","RBIAS(%)","SD","SE","COV(%)")
    row.names(output)=c("FO","PO","PMI","NNMI")
    print("Scenario:")
    print(choice)
    print("Censoring rate (%):")
    print(round(100*crate,2))
    output
}

simulation(rep=500,N=200,b0 = 4, b1 = -2, b2 = 0.5, b3 = -2, b4 = 2, b5 = 2, r0 = 3, r1 = -3, r2 = 0.5, r3 = -2, r4 = 1.5, r5 = 2,choice=1,N.MI=10,n.s=10,w1=0.8,w2=0.2)

# ouput for Table 8.4
