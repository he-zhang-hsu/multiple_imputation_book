
/*libname a "\\cdc.gov\private\M132\jnu5\QSDMB\Incidence Model\Multiple Imputation\Rick_tran_cate\data";
options noxwait nofmterr; 


*Before the real simulation, get the marginal probability and correlations;

data MI_M;
	set a.MI_M;
	if r5=1 or R07=1 or R08=1 or R09=1 or R10=1 then r5_new=1;
	else r5_new=0;
	if r5=. and R07=. and R08=. and R09=. and R10=. then r5_new=.;
run;

proc print data=MI_M (obs=10);
	var r1 r3 r5 R07 R08 R09 R10 r5_new race1;
run;

proc freq data=MI_M; 
	tables r1 r3 r5_new r5 R07 R08 R09 R10 race1;
run;

proc corr data=MI_M;
	var r1 r3 r5_new race1;
run;



*look at the missing proportion of each category, assume we do not have missing for the covariate;

proc freq data=MI_M; tables r1 r3 r5_new race1;run;

*/

/*read in the simulation data;

data b.data4; set data4;run;

data data; set b.data4;run;

*add the simulation number and subject ID to it, generate the id first but itselves;

data id;
	do simu=1 to 1000;
		do id=1 to 999;*if we put 1000 then we have 100100 observations;
		output;
		end;
		output;
	end;
run;


*merge it back to the data;
data b.simu_data4(rename=(V2=r1 V3=r3 V4=r5 VAR5=co1));
	set data; set id;
run;

proc means data=b.simu_data4; 
	var r1 r3 r5 co1;
run;

proc corr data=b.simu_data4;
var r1 r3 r5 co1;
run;


proc print data=b.simu_data4 (obs=10); run;

*/




*MAR;
%macro simul(nsim); 

libname b "\\cdc.gov\private\M132\jnu5\QSDMB\Incidence Model\Multiple Imputation\Rick_tran_cate\program\Simulation\Data Generation";

libname a "\\cdc.gov\private\M132\jnu5\QSDMB\Incidence Model\Multiple Imputation\Rick_tran_cate\program\Simulation\result\rick\set1\missing7";


data truth;
	set b.simu_rick1;
	if r1=1 and r3=1 then xmode=1;
	else if r1=1 and r3=0 then xmode=2;
	else if r1=0 and r3=1 then xmode=3;
	else if r1=0 and r3=0 and r5=1 then xmode=4;
	else if r1=0 and r3=0 and r5=0 then xmode=9;
	else xmode=.;
run;

proc freq  data=truth; tables xmode/nocum ;
	ods output Freq.Table1.OneWayFreqs=truth2;
run;


data a.truth; set truth2; truth=percent/100;keep xmode truth;run;


*start simulation;

ods listing close;                                                                                                                    
ods select none;                                                                                                                        
filename junk dummy;                                                                                                                    
proc printto  log=junk; run;                                                                                        
                                                                                                                                        
%do simu=1 %to &nsim;                                                                                                                    
                                                                                                                                        
proc datasets library=work kill;run;     

data simu; 
	set b.simu_rick1;
	if simu=&simu;
run;

                                                                                                                                        
data samplesize;                                                                                                                        
      seedsample1=&simu;                                                                                                                 
      seedsample2=2*&simu;
 	  seedsample3=3*&simu;
      call symput("seedsample1",trim(left(put(seedsample1,8.))));                                                                       
      call symput("seedsample2",trim(left(put(seedsample2,8.))));  
	  call symput("seedsample3",trim(left(put(seedsample3,8.))));  
run;                                                                                                                                   

                                                                                                                          
*The next, lets set 20% MCAR missing to r1 and 40% MCAR missing to R02;

*set up missing for ;  
data sample1; set simu;
	if co1=1 then do;
         bin1=ranbin(&seedsample1.,1,0.4);
         OUTPUT;
      END;

	  if co1=0 then do;
		bin2=ranbin(&seedsample2.,1,0.8);
		bin3=ranbin(&seedsample2.,1,0.7);
         OUTPUT;
      END;
run;

data MI_M3; set sample1; 
	if bin1=1 then r1=r1; 
	if bin1=0 then r1=.;
	if bin2=1 then r3=r3;
	if bin2=0 then r3=.;
	if bin3=1 then r5=r5;
	if bin3=0 then r5=.;
run;
               


*based on the missing data on r1 and r3, set missing data for xmode;

data MI_M4;
	set MI_M3;
	if r1=1 and r3=1 then Trans=1;
	else if r1=1 and r3=0 then Trans=2;
	else if r1=0 and r3=1 then Trans=3;
	else if r1=0 and r3=0 and r5=1 then Trans=4;
	else if r1=0 and r3=0 and r5=0 then Trans=9;
	else Trans=.;
	drop xmode;
run;

data MI_M4(rename=(trans=xmode)); set MI_M4;run;

proc freq data=mi_m4; tables xmode;run;



*Method 1.1, impute categorical variable;
* JAV imputation;
proc mi data=MI_M4
     out=xtemp  /* output dataset containing imputation results */
     nimpute=5 /* number of imputations, change if necessary */
     seed=1234 /* same seed will generate same results */
    noprint;   /*suppress printing results in output window */
  class xmode ;*race1-race3 region1-region9 rmsa1-rmsa2 hivcat1-hivcat2 origin1-origin2 facility1-facility5 delaygrp1-delaygrp2;
  fcs discrim(xmode= r1 r3 r5 co1);
var  r1 r3 r5 co1 xmode ;
run;



  *convert xmode to dummy variables;
  data xtemp2;
  	set xtemp;
	if xmode=1 then xmode1=1;else xmode1=0;
	if xmode=2 then xmode2=1;else xmode2=0; 
	if xmode=3 then xmode3=1;else xmode3=0;
	if xmode=4 then xmode4=1;else xmode4=0;
	if xmode=9 then xmode9=1;else xmode9=0;
run;



proc means data=xtemp2 ;
var xmode1 xmode2 xmode3 xmode4 xmode9;
output out=xout1 mean=mean_xmode1-mean_xmode4 mean_xmode9 
stderr=se_xmode1-se_xmode4 se_xmode9;
by _imputation_;
run;

proc mianalyze data=xout1;
modeleffects mean_xmode1 mean_xmode2 mean_xmode3 mean_xmode4 mean_xmode9 ;
stderr se_xmode1 se_xmode2 se_xmode3 se_xmode4 se_xmode9;
ods output Mianalyze.ParameterEstimates=method1;
run;


data method1_new; 
	set method1;
	if parm="mean_xmode1" then xmode=1;
	if parm="mean_xmode2" then xmode=2;
	if parm="mean_xmode3" then xmode=3;
	if parm="mean_xmode4" then xmode=4;
	if parm="mean_xmode9" then xmode=9;
	simu=&simu.;
	keep xmode Estimate StdErr LCLMean UCLMean simu;
run;


proc append base=a.method1_final data=method1_new;run;

*Method 1.2, impute categorical variable;
* active imputation;
proc mi data=MI_M4
     out=xtemp_new  /* output dataset containing imputation results */
     nimpute=5 /* number of imputations, change if necessary */
     seed=1234 /* same seed will generate same results */
    noprint;   /*suppress printing results in output window */
  class xmode ;*race1-race3 region1-region9 rmsa1-rmsa2 hivcat1-hivcat2 origin1-origin2 facility1-facility5 delaygrp1-delaygrp2;
  fcs discrim(xmode= co1);
var  co1 xmode ;
run;



  *convert xmode to dummy variables;
  data xtemp2_new;
  	set xtemp_new;
	if xmode=1 then xmode1=1;else xmode1=0;
	if xmode=2 then xmode2=1;else xmode2=0; 
	if xmode=3 then xmode3=1;else xmode3=0;
	if xmode=4 then xmode4=1;else xmode4=0;
	if xmode=9 then xmode9=1;else xmode9=0;
run;



proc means data=xtemp2_new ;
var xmode1 xmode2 xmode3 xmode4 xmode9;
output out=xout1_new mean=mean_xmode1-mean_xmode4 mean_xmode9 
stderr=se_xmode1-se_xmode4 se_xmode9;
by _imputation_;
run;

proc mianalyze data=xout1_new;
modeleffects mean_xmode1 mean_xmode2 mean_xmode3 mean_xmode4 mean_xmode9 ;
stderr se_xmode1 se_xmode2 se_xmode3 se_xmode4 se_xmode9;
ods output Mianalyze.ParameterEstimates=method12;
run;


data method12_new; 
	set method12;
	if parm="mean_xmode1" then xmode=1;
	if parm="mean_xmode2" then xmode=2;
	if parm="mean_xmode3" then xmode=3;
	if parm="mean_xmode4" then xmode=4;
	if parm="mean_xmode9" then xmode=9;
	simu=&simu.;
	keep xmode Estimate StdErr LCLMean UCLMean simu;
run;


proc append base=a.method12_final data=method12_new;run;

*Method2.11, impute r1, r3 and r5 simutaneously, no interaction;
* Passive II;
data Method2_1;set MI_M4;run;

proc mi data=Method2_1 out=xtemp3  nimpute=5;
class r1 r3 r5;*r1 r3 r5;
*fcs logistic(xmode= r1 r3 r5 co1 / details);
fcs logistic(r1= r3 r5 co1 / details); 
fcs logistic(r3= r1 r5 co1 / details); 
fcs logistic(r5= r1 r3 co1 / details); 

/* fcs discrim(r1= r3 r5 co1);
fcs discrim(r3= r1 r5 co1);
fcs discrim(r5= r1 r3 co1); */
var co1 r1 r3 r5;
run;


data xtemp4;
	set xtemp3;
	if r1=1 and r3=1 then Trans=1;
	else if r1=1 and r3=0 then Trans=2;
	else if r1=0 and r3=1 then Trans=3;
	else if r1=0 and r3=0 and r5=1 then Trans=4;
	else if r1=0 and r3=0 and r5=0 then Trans=9;
	else Trans=.;
	drop xmode;
run;

data xtemp4(rename=(trans=xmode)); set xtemp4;run;

  *convert xmode to dummy variables;
  data xtemp5;
  	set xtemp4;
	if xmode=1 then xmode1=1;else xmode1=0;
	if xmode=2 then xmode2=1;else xmode2=0; 
	if xmode=3 then xmode3=1;else xmode3=0;
	if xmode=4 then xmode4=1;else xmode4=0;
	if xmode=9 then xmode9=1;else xmode9=0;
run;


proc means data=xtemp5;
var xmode1 xmode2 xmode3 xmode4 xmode9;
output out=xout2 mean=mean_xmode1-mean_xmode4 mean_xmode9 
stderr=se_xmode1-se_xmode4 se_xmode9;
by _imputation_;
run;

proc mianalyze data=xout2;
modeleffects mean_xmode1 mean_xmode2 mean_xmode3 mean_xmode4 mean_xmode9 ;
stderr se_xmode1 se_xmode2 se_xmode3 se_xmode4 se_xmode9;
ods output Mianalyze.ParameterEstimates=method21;
run;


data method21_new; 
	set method21;
	if parm="mean_xmode1" then xmode=1;
	if parm="mean_xmode2" then xmode=2;
	if parm="mean_xmode3" then xmode=3;
	if parm="mean_xmode4" then xmode=4;
	if parm="mean_xmode9" then xmode=9;
	simu=&simu.;
	keep xmode Estimate StdErr LCLMean UCLMean simu;
run;


proc append base=a.method21_final data=method21_new;run;

*Method 2.12, have intercation terms;
* Passive III;
proc mi data=Method2_1  nimpute=5 out=xtemp31;
class r1 r3 r5;
fcs logistic(r1= r3 r5 co1 r3*r5 r3*co1 r5*co1 r3*r5*co1 / details); 
fcs logistic(r3= r1 r5 co1 r1*r5 r1*co1 r5*co1 r1*r5*co1 / details);
fcs logistic(r5= r1 r3 co1 r1*r3 r1*co1 r3*co1 r1*r3*co1 / details); 
var r1 r3 r5 co1;
run;



data xtemp41;
	set xtemp31;
	if r1=1 and r3=1 then Trans=1;
	else if r1=1 and r3=0 then Trans=2;
	else if r1=0 and r3=1 then Trans=3;
	else if r1=0 and r3=0 and r5=1 then Trans=4;
	else if r1=0 and r3=0 and r5=0 then Trans=9;
	else Trans=.;
	drop xmode;
run;

data xtemp41(rename=(trans=xmode)); set xtemp41;run;

  *convert xmode to dummy variables;
  data xtemp51;
  	set xtemp41;
	if xmode=1 then xmode1=1;else xmode1=0;
	if xmode=2 then xmode2=1;else xmode2=0; 
	if xmode=3 then xmode3=1;else xmode3=0;
	if xmode=4 then xmode4=1;else xmode4=0;
	if xmode=9 then xmode9=1;else xmode9=0;
run;


proc means data=xtemp51;
var xmode1 xmode2 xmode3 xmode4 xmode9;
output out=xout212 mean=mean_xmode1-mean_xmode4 mean_xmode9 
stderr=se_xmode1-se_xmode4 se_xmode9;
by _imputation_;
run;

proc mianalyze data=xout212;
modeleffects mean_xmode1 mean_xmode2 mean_xmode3 mean_xmode4 mean_xmode9 ;
stderr se_xmode1 se_xmode2 se_xmode3 se_xmode4 se_xmode9;
ods output Mianalyze.ParameterEstimates=method212;
run;


data method212_new; 
	set method212;
	if parm="mean_xmode1" then xmode=1;
	if parm="mean_xmode2" then xmode=2;
	if parm="mean_xmode3" then xmode=3;
	if parm="mean_xmode4" then xmode=4;
	if parm="mean_xmode9" then xmode=9;
	simu=&simu.;
	keep xmode Estimate StdErr LCLMean UCLMean simu;
run;


proc append base=a.method212_final data=method212_new;run;


*Method2_2, follow the hiarachy of the rules;
* Passive I;
data Method2_2;set MI_M4;run;


*impute r1 and r3 first;
proc mi data=MI_M4
     out=xtemp6  /* output dataset containing imputation results */
     nimpute=5 /* number of imputations, change if necessary */
     seed=1234 /* same seed will generate same results */
    noprint;   /*suppress printing results in output window */
  class r1 r3 co1;*race1-race3 region1-region9 rmsa1-rmsa2 hivcat1-hivcat2 origin1-origin2 facility1-facility5 delaygrp1-delaygrp2;
 fcs logistic(r1  = r3 co1 r3*co1);
  fcs logistic(r3 = r1 co1 r1*co1);
var  r1 r3 co1 ;
run;


*imputation1;
data xtemp6_1; set xtemp6; if _imputation_=1;run;

data xtemp6_11; set xtemp6_1;if r1~=0 or r3~=0;run;

*r1=r3=0 then impute r5;
data xtemp6_12; set xtemp6_1; if r1=0 and r3=0;run;

proc mi data=xtemp6_12
     out=xtemp6_12_new  /* output dataset containing imputation results */
     nimpute=1 /* number of imputations, change if necessary */
     seed=1234 /* same seed will generate same results */
    noprint;   /*suppress printing results in output window */
  class r5;*race1-race3 region1-region9 rmsa1-rmsa2 hivcat1-hivcat2 origin1-origin2 facility1-facility5 delaygrp1-delaygrp2;
  fcs logistic(r5= co1);
var  co1 r5;
run;

*merge this dataset with the part which does not need to impute;
data xtemp6_1_new; set xtemp6_12_new xtemp6_11;drop xmode;run;


*imputation2;
data xtemp6_2; set xtemp6; if _imputation_=2;run;

data xtemp6_21; set xtemp6_2;if r1~=0 or r3~=0;run;

*r1=r3=0 then impute r5;
data xtemp6_22; set xtemp6_2; if r1=0 and r3=0;run;

proc mi data=xtemp6_22
     out=xtemp6_22_new  /* output dataset containing imputation results */
     nimpute=1 /* number of imputations, change if necessary */
     seed=1234 /* same seed will generate same results */
    noprint;   /*suppress printing results in output window */
  class r5;*race1-race3 region1-region9 rmsa1-rmsa2 hivcat1-hivcat2 origin1-origin2 facility1-facility5 delaygrp1-delaygrp2;
  fcs logistic(r5= co1);
var  co1 r5;
run;

*merge this dataset with the part which does not need to impute;
data xtemp6_2_new; set xtemp6_22_new xtemp6_21;drop xmode;run;



*imputation3;
data xtemp6_3; set xtemp6; if _imputation_=3;run;

data xtemp6_31; set xtemp6_3;if r1~=0 or r3~=0;run;

*r1=r3=0 then impute r5;
data xtemp6_32; set xtemp6_3; if r1=0 and r3=0;run;

proc mi data=xtemp6_32
     out=xtemp6_32_new  /* output dataset containing imputation results */
     nimpute=1 /* number of imputations, change if necessary */
     seed=1234 /* same seed will generate same results */
    noprint;   /*suppress printing results in output window */
  class r5;*race1-race3 region1-region9 rmsa1-rmsa2 hivcat1-hivcat2 origin1-origin2 facility1-facility5 delaygrp1-delaygrp2;
  fcs logistic(r5= co1);
var  co1 r5;
run;

*merge this dataset with the part which does not need to impute;
data xtemp6_3_new; set xtemp6_32_new xtemp6_31;drop xmode;run;



*imputation4;
data xtemp6_4; set xtemp6; if _imputation_=4;run;

data xtemp6_41; set xtemp6_4;if r1~=0 or r3~=0;run;

*r1=r3=0 then impute r5;
data xtemp6_42; set xtemp6_4; if r1=0 and r3=0;run;

proc mi data=xtemp6_42
     out=xtemp6_42_new  /* output dataset containing imputation results */
     nimpute=1 /* number of imputations, change if necessary */
     seed=1234 /* same seed will generate same results */
    noprint;   /*suppress printing results in output window */
  class r5;*race1-race3 region1-region9 rmsa1-rmsa2 hivcat1-hivcat2 origin1-origin2 facility1-facility5 delaygrp1-delaygrp2;
  fcs logistic(r5= co1);
var  co1 r5;
run;

*merge this dataset with the part which does not need to impute;
data xtemp6_4_new; set xtemp6_42_new xtemp6_41;drop xmode;run;



*imputation5;
data xtemp6_5; set xtemp6; if _imputation_=5;run;

data xtemp6_51; set xtemp6_5;if r1~=0 or r3~=0;run;

*r1=r3=0 then impute r5;
data xtemp6_52; set xtemp6_5; if r1=0 and r3=0;run;

proc mi data=xtemp6_52
     out=xtemp6_52_new  /* output dataset containing imputation results */
     nimpute=1 /* number of imputations, change if necessary */
     seed=1234 /* same seed will generate same results */
    noprint;   /*suppress printing results in output window */
  class r5;*race1-race3 region1-region9 rmsa1-rmsa2 hivcat1-hivcat2 origin1-origin2 facility1-facility5 delaygrp1-delaygrp2;
  fcs logistic(r5= co1);
var  co1 r5;
run;

*merge this dataset with the part which does not need to impute;
data xtemp6_5_new; set xtemp6_52_new xtemp6_51;drop xmode;run;


*merge all 5 imputations;

data xtemp7; set xtemp6_1_new xtemp6_2_new xtemp6_3_new xtemp6_4_new xtemp6_5_new;run;


*look at the performance of r1, r3 and r5;
proc means data=xtemp7;
var r1 r3 r5;
output out=xout_r2 mean=mean_r1 mean_r3 mean_r5 
stderr=se_r1 se_r3 se_r5;
by _imputation_;
run;

proc mianalyze data=xout_r2;
modeleffects mean_r1 mean_r3 mean_r5  ;
stderr se_r1 se_r3 se_r5;
ods output Mianalyze.ParameterEstimates=method22_r;
run;

data method22_r;
set method22_r;
simu=&simu.;
*keep r1 r3 r5 Estimate StdErr LCLMean UCLMean simu;
run;

proc append base=a.method22_r data=method22_r;run;

data xtemp8;
	set xtemp7;
	if r1=1 and r3=1 then xmode=1;
	else if r1=1 and r3=0 then xmode=2;
	else if r1=0 and r3=1 then xmode=3;
	else if r1=0 and r3=0 and r5=1 then xmode=4;
	else if r1=0 and r3=0 and r5=0 then xmode=9;
	else xmode=.;
run;

  *convert xmode to dummy variables;
  data xtemp9;
  	set xtemp8;
	if xmode=1 then xmode1=1;else xmode1=0;
	if xmode=2 then xmode2=1;else xmode2=0; 
	if xmode=3 then xmode3=1;else xmode3=0;
	if xmode=4 then xmode4=1;else xmode4=0;
	if xmode=9 then xmode9=1;else xmode9=0;
run;


proc means data=xtemp9;
var xmode1 xmode2 xmode3 xmode4 xmode9;
output out=xout3 mean=mean_xmode1-mean_xmode4 mean_xmode9 
stderr=se_xmode1-se_xmode4 se_xmode9;
by _imputation_;
run;

proc mianalyze data=xout3;
modeleffects mean_xmode1 mean_xmode2 mean_xmode3 mean_xmode4 mean_xmode9 ;
stderr se_xmode1 se_xmode2 se_xmode3 se_xmode4 se_xmode9;
ods output Mianalyze.ParameterEstimates=method22;
run;


data method22_new; 
	set method22;
	if parm="mean_xmode1" then xmode=1;
	if parm="mean_xmode2" then xmode=2;
	if parm="mean_xmode3" then xmode=3;
	if parm="mean_xmode4" then xmode=4;
	if parm="mean_xmode9" then xmode=9;
	simu=&simu.;
	keep xmode Estimate StdErr LCLMean UCLMean simu;
run;


proc append base=a.method22_final data=method22_new;run;

%end;

/* Table 11.16 */

*summarize method 1;
* JAV

proc sort data=a.truth; by xmode;run;
proc sort data=a.method1_final;by xmode;run;

data a.together1; merge a.truth a.method1_final; 
	by xmode;
bias=estimate-truth;
if LCLMean<=truth<= UCLMean then cover=1;
else cover=0; 
run;

*summarize method 12;
* active imputation;

proc sort data=a.truth; by xmode;run;
proc sort data=a.method12_final;by xmode;run;

data a.together12; merge a.truth a.method12_final; 
	by xmode;
bias=estimate-truth;
if LCLMean<=truth<= UCLMean then cover=1;
else cover=0; 
run;


*summarize method 21;
* Passive II
proc sort data=a.truth; by xmode;run;
proc sort data=a.method21_final;by xmode;run;

data a.together21; merge a.truth a.method21_final; 
	by xmode;
bias=estimate-truth;
if LCLMean<=truth<= UCLMean then cover=1;
else cover=0; 
run;


*summarize method 212;
* Passive III;
proc sort data=a.truth; by xmode;run;
proc sort data=a.method212_final;by xmode;run;

data a.together212; merge a.truth a.method212_final; 
	by xmode;
bias=estimate-truth;
if LCLMean<=truth<= UCLMean then cover=1;
else cover=0; 
run;


*summarize method 22;
* Passive I;
proc sort data=a.truth; by xmode;run;
proc sort data=a.method22_final;by xmode;run;

data a.together22; merge a.truth a.method22_final; 
	by xmode;
bias=estimate-truth;
if LCLMean<=truth<= UCLMean then cover=1;
else cover=0; 
run;

%mend;

%simul(nsim=1000);















