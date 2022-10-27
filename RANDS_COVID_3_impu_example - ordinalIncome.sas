
/* the data RANDS_COVID_3.txt can be downloaded at https://www.cdc.gov/nchs/rands/data.htm */

libname d "C:\Users\wdq7\OneDrive - CDC\+My_Large_Workspace\machine_learning\Project\Imputation_software";
data d.RANDS_COVID_3;
infile 'C:\Users\wdq7\OneDrive - CDC\+My_Large_Workspace\machine_learning\Project\Imputation_software\RANDS_COVID_3.txt';
input     AGE 2. +1
    ALT_TELMED 2. +1
    ALT_TELMED_TOTALTIME 3. +1
    ANXFREQ 2. +1
    ANXFREQ_TOTALTIME 3. +1
    ANXLEVEL 2. +1
    ANXLEVEL_TOTALTIME 3. +1
    ANXMED 2. +1
    ANXMED_TOTALTIME 3. +1
    ASEV 2. +1
    ASTILL 2. +1
    BURDEN1 2. +1
    BURDEN2 2. +1
    CANEV 2. +1
    CHDEV 2. +1
    CHLEV 2. +1
    COGMEMDIFF 2. +1
    COMDIFF 2. +1
    COPDEV 2. +1
    COVIDNOCAR_A 2. +1
    COVIDNOCAR_B 2. +1
    COVIDNOCAR_C 2. +1
    COVIDNOCAR_D 2. +1
    COVIDNOCAR_E 2. +1
    COVIDNOCAR_F 2. +1
    COVIDNOCAR_G 2. +1
    COVIDNOCAR_H 2. +1
    COVIDNOCAR_I 2. +1
    COVID_HES 2. +1
    COVID_IMPACT1 2. +1
    COVID_IMPACT2 2. +1
    COVID_IMPACT3 2. +1
    COVID_IMPACT4 2. +1
    COVID_NOWK 2. +1
    COVID_S6A_A 2. +1
    COVID_S6A_B 2. +1
    COVID_S6A_C 2. +1
    COVID_S6A_D 2. +1
    COVID_S6A_E 2. +1
    COVID_S6A_F 2. +1
    COVID_S6A_G 2. +1
    COVID_S6A_H 2. +1
    COVID_S6A_I 2. +1
    COVID_S6B_A 2. +1
    COVID_S6B_B 2. +1
    COVID_S6B_C 2. +1
    COVID_S6B_E 2. +1
    COVID_S6B_G 2. +1
    COVID_S6B_H 2. +1
    COVID_S6B_I 2. +1
    COVID_S6C_A 2. +1
    COVID_S6C_B 2. +1
    COVID_S6C_C 2. +1
    COVID_S6C_D 2. +1
    COVID_S6C_E 2. +1
    COVID_S6C_F 2. +1
    COVID_S6C_G 2. +1
    COVID_S6C_H 2. +1
    COVID_S6C_I 2. +1
    CaseID 5. +1
    DEPFREQ 2. +1
    DEPFREQ_TOTALTIME 3. +1
    DEPLEVEL 2. +1
    DEPLEVEL_TOTALTIME 3. +1
    DEPMED 2. +1
    DEPMED_TOTALTIME 3. +1
    DIBEV 2. +1
    DIFF 2. +1
    DNGCARE 2. +1
    EDUC 1. +1
    EMPLASTWK 2. +1
    EMPLOY 1. +1
    GAD7_A 2. +1
    GAD7_B 2. +1
    GAD7_TOTALTIME 3. +1
    GENDER 1. +1
    HEARINGDF 2. +1
    HHSIZE 1. +1
    HICOV 2. +1
    HOME_TYPE 1. +1
    HOUSING 1. +1
    HYPEV 2. +1
    INCOME 2. +1
    INTERNET 1. +1
    LIFESAT1 2. +1
    LIFESAT1_TOTALTIME 3. +1
    LIFESAT2 2. +1
    LIFESAT2_TOTALTIME 3. +1
    MARITAL 1. +1
    MODE_PREF $4. +1
    NOCARTYP_A 2. +1
    NOCARTYP_B 2. +1
    NOCARTYP_C 2. +1
    NOCARTYP_D 2. +1
    NOCARTYP_E 2. +1
    NOCARTYP_F 2. +1
    NOCARTYP_G 2. +1
    NOCARTYP_H 2. +1
    NOCARTYP_I 2. +1
    PHONESERVICE 1. +1
    PHQ_A 2. +1
    PHQ_B 2. +1
    PHQ_TOTALTIME 4. +1
    PHSTAT 2. +1
    PLIFESAT_TOTALTIME 4. +1
    PREDIB 2. +1
    PROBE_ANX_1 1. +1
    PROBE_ANX_2 1. +1
    PROBE_ANX_3 1. +1
    PROBE_ANX_4 1. +1
    PROBE_ANX_5 1. +1
    PROBE_ANX_6 1. +1
    PROBE_ANX_7 1. +1
    PROBE_ANX_DK 1. +1
    PROBE_ANX_REF 1. +1
    PROBE_DEP_1 1. +1
    PROBE_DEP_2 1. +1
    PROBE_DEP_3 1. +1
    PROBE_DEP_4 1. +1
    PROBE_DEP_5 1. +1
    PROBE_DEP_6 1. +1
    PROBE_DEP_7 1. +1
    PROBE_DEP_DK 1. +1
    PROBE_DEP_REF 1. +1
    PROBE_DIST_TIMER_TOTALTIME 4. +1
    PROBE_GETCAR 2. +1
    PROBE_LIFESAT_1 1. +1
    PROBE_LIFESAT_2 1. +1
    PROBE_LIFESAT_3 1. +1
    PROBE_LIFESAT_4 1. +1
    PROBE_LIFESAT_5 1. +1
    PROBE_LIFESAT_DK 1. +1
    PROBE_LIFESAT_REF 1. +1
    PROBE_MASKTYPE_A 2. +1
    PROBE_MASKTYPE_B 2. +1
    PROBE_MASKTYPE_C 2. +1
    PROBE_MASKTYPE_D 2. +1
    PROBE_MASKTYPE_E 2. +1
    PROBE_MASKTYPE_F 2. +1
    PROBE_MASKUSE 2. +1
    PROBE_NOCAR_A 2. +1
    PROBE_NOCAR_B 2. +1
    PROBE_NOCAR_C 2. +1
    PROBE_NOCAR_D 2. +1
    PROBE_NOCAR_E 2. +1
    PROBE_REL6_TIMER_TOTALTIME 4. +1
    PROBE_TELMEDALT1 2. +1
    PROBE_TELMEDALT2 2. +1
    PROBE_TELMEDALT3 2. +1
    PROBE_TELMEDORIG_1 1. +1
    PROBE_TELMEDORIG_2 1. +1
    PROBE_TELMEDORIG_3 1. +1
    PROBE_TELMEDORIG_DK 1. +1
    PROBE_TELMEDORIG_REF 1. +1
    PROBE_VAX_1 1. +1
    PROBE_VAX_10 1. +1
    PROBE_VAX_2 1. +1
    PROBE_VAX_3 1. +1
    PROBE_VAX_4 1. +1
    PROBE_VAX_5 1. +1
    PROBE_VAX_6 1. +1
    PROBE_VAX_7 1. +1
    PROBE_VAX_8 1. +1
    PROBE_VAX_9 1. +1
    PROBE_VAX_DK 1. +1
    PROBE_VAX_REF 1. +1
    PROBE_VAX_TIMER_TOTALTIME 4. +1
    P_LIFEOPEN 1. +1
    P_LIFESAT 1. +1
    P_PROBEANXDEP 1. +1
    P_PROBETEL 1. +1
    P_TELMEDEXP 1. +1
    P_VAXOPEN 1. +1
    RACETHNICITY 1. +1
    REGION4 1. +1
    REL1 2. +1
    REL2 2. +1
    REL3 2. +1
    REL4 2. +1
    REL5 2. +1
    REL6 2. +1
    REL7 2. +1
    REL_TOTALTIME 5. +1
    SMKEV 2. +1
    SMKNOW 2. +1
    SUM_GAD7 1. +1
    SUM_PHQ 1. +1
    SURV_MODE 1. +1
    S_VPSU 1. +1
    S_VSTRAT 2. +1
    TELMED 2. +1
    TELMEDALT1_TOTALTIME 4. +1
    TELMEDALT2_TOTALTIME 3. +1
    TELMEDALT3_TOTALTIME 3. +1
    TELMEDNEW 2. +1
    TELMEDORIG_TOTALTIME 4. +1
    TELMEDUSE 2. +1
    TELMED_TOTALTIME 3. +1
    UPPSLFCR 2. +1
    USPLKIND 2. +1
    USUALPL 2. +1
    VAX_BELIEF_A 2. +1
    VAX_BELIEF_B 2. +1
    VAX_BELIEF_C 2. +1
    VAX_BELIEF_D 2. +1
    VAX_BELIEF_E 2. +1
    VAX_BELIEF_F 2. +1
    VAX_BELIEF_G 2. +1
    VAX_COVID 2. +1
    VAX_HERD 2. +1
    VAX_HES 2. +1
    VAX_KNOW 2. +1
    VAX_MD 2. +1
    VAX_PLAN 2. +1
    VAX_RISK 2. +1
    VAX_SIDE 2. +1
    VAX_TOTALTIME 5. +1
    VISIONDF 2. +1
    WEIGHT 8.5 +1
    WEIGHT_CALIBRATED 8.5 +1
    duration 3. ;
* optional - create formats with separate program (RANDS_COVID_3_formats.sas) first, uncomment code below to associate formats ;
/*
format P_LIFESAT P_LIFESAT. ;
format P_LIFEOPEN P_LIFEOPEN. ;
format P_VAXOPEN P_VAXOPEN. ;
format P_TELMEDEXP P_TELMEDEXP. ;
format P_PROBETEL P_PROBETEL. ;
format P_PROBEANXDEP P_PROBEANXDEP. ;
format PHSTAT PHSTAT. ;
format LIFESAT1 LIFESAT1x. ;
format LIFESAT2 LIFESAT2x. ;
format PROBE_LIFESAT_1 PROBE_LIFESAT_1x. ;
format PROBE_LIFESAT_2 PROBE_LIFESAT_2x. ;
format PROBE_LIFESAT_3 PROBE_LIFESAT_3x. ;
format PROBE_LIFESAT_4 PROBE_LIFESAT_4x. ;
format PROBE_LIFESAT_5 PROBE_LIFESAT_5x. ;
format PROBE_LIFESAT_DK PROBE_LIFESAT_DK. ;
format PROBE_LIFESAT_REF PROBE_LIFESAT_REF. ;
format HYPEV HYPEV. ;
format CHLEV CHLEV. ;
format CHDEV CHDEV. ;
format ASEV ASEV. ;
format COPDEV COPDEV. ;
format CANEV CANEV. ;
format ASTILL ASTILL. ;
format PREDIB PREDIB. ;
format DIBEV DIBEV. ;
format SMKEV SMKEV. ;
format SMKNOW SMKNOW. ;
format EMPLASTWK EMPLASTWK. ;
format COVID_NOWK COVID_NOWK. ;
format HICOV HICOV. ;
format DNGCARE DNGCARE. ;
format USUALPL USUALPL. ;
format USPLKIND USPLKIND. ;
format TELMED TELMED. ;
format ALT_TELMED ALT_TELMED. ;
format PROBE_TELMEDORIG_1 PROBE_TELMEDORIG_1x. ;
format PROBE_TELMEDORIG_2 PROBE_TELMEDORIG_2x. ;
format PROBE_TELMEDORIG_3 PROBE_TELMEDORIG_3x. ;
format PROBE_TELMEDORIG_DK PROBE_TELMEDORIG_DK. ;
format PROBE_TELMEDORIG_REF PROBE_TELMEDORIG_REF. ;
format PROBE_TELMEDALT1 PROBE_TELMEDALT1x. ;
format PROBE_TELMEDALT2 PROBE_TELMEDALT2x. ;
format PROBE_TELMEDALT3 PROBE_TELMEDALT3x. ;
format TELMEDUSE TELMEDUSE. ;
format TELMEDNEW TELMEDNEW. ;
format NOCARTYP_A NOCARTYP_A. ;
format NOCARTYP_B NOCARTYP_B. ;
format NOCARTYP_C NOCARTYP_C. ;
format NOCARTYP_D NOCARTYP_D. ;
format NOCARTYP_E NOCARTYP_E. ;
format NOCARTYP_F NOCARTYP_F. ;
format NOCARTYP_G NOCARTYP_G. ;
format NOCARTYP_H NOCARTYP_H. ;
format NOCARTYP_I NOCARTYP_I. ;
format COVIDNOCAR_A COVIDNOCAR_A. ;
format COVIDNOCAR_B COVIDNOCAR_B. ;
format COVIDNOCAR_C COVIDNOCAR_C. ;
format COVIDNOCAR_D COVIDNOCAR_D. ;
format COVIDNOCAR_E COVIDNOCAR_E. ;
format COVIDNOCAR_F COVIDNOCAR_F. ;
format COVIDNOCAR_G COVIDNOCAR_G. ;
format COVIDNOCAR_H COVIDNOCAR_H. ;
format COVIDNOCAR_I COVIDNOCAR_I. ;
format PROBE_NOCAR_A PROBE_NOCAR_A. ;
format PROBE_NOCAR_B PROBE_NOCAR_B. ;
format PROBE_NOCAR_C PROBE_NOCAR_C. ;
format PROBE_NOCAR_D PROBE_NOCAR_D. ;
format PROBE_NOCAR_E PROBE_NOCAR_E. ;
format PROBE_GETCAR PROBE_GETCAR. ;
format VAX_HES VAX_HES. ;
format VAX_SIDE VAX_SIDE. ;
format VAX_KNOW VAX_KNOW. ;
format VAX_MD VAX_MD. ;
format VAX_RISK VAX_RISK. ;
format VAX_HERD VAX_HERD. ;
format COVID_HES COVID_HES. ;
format VAX_COVID VAX_COVID. ;
format VAX_PLAN VAX_PLAN. ;
format PROBE_VAX_1 PROBE_VAX_1x. ;
format PROBE_VAX_2 PROBE_VAX_2x. ;
format PROBE_VAX_3 PROBE_VAX_3x. ;
format PROBE_VAX_4 PROBE_VAX_4x. ;
format PROBE_VAX_5 PROBE_VAX_5x. ;
format PROBE_VAX_6 PROBE_VAX_6x. ;
format PROBE_VAX_7 PROBE_VAX_7x. ;
format PROBE_VAX_8 PROBE_VAX_8x. ;
format PROBE_VAX_9 PROBE_VAX_9x. ;
format PROBE_VAX_10 PROBE_VAX_10x. ;
format PROBE_VAX_DK PROBE_VAX_DK. ;
format PROBE_VAX_REF PROBE_VAX_REF. ;
format VAX_BELIEF_A VAX_BELIEF_A. ;
format VAX_BELIEF_B VAX_BELIEF_B. ;
format VAX_BELIEF_C VAX_BELIEF_C. ;
format VAX_BELIEF_D VAX_BELIEF_D. ;
format VAX_BELIEF_E VAX_BELIEF_E. ;
format VAX_BELIEF_F VAX_BELIEF_F. ;
format VAX_BELIEF_G VAX_BELIEF_G. ;
format VAX_TOTALTIME VAX_TOTALTIME. ;
format COVID_S6A_A COVID_S6A_A. ;
format COVID_S6A_B COVID_S6A_B. ;
format COVID_S6A_C COVID_S6A_C. ;
format COVID_S6A_D COVID_S6A_D. ;
format COVID_S6A_E COVID_S6A_E. ;
format COVID_S6A_F COVID_S6A_F. ;
format COVID_S6A_G COVID_S6A_G. ;
format COVID_S6A_H COVID_S6A_H. ;
format COVID_S6A_I COVID_S6A_I. ;
format COVID_S6B_A COVID_S6B_A. ;
format COVID_S6B_B COVID_S6B_B. ;
format COVID_S6B_C COVID_S6B_C. ;
format COVID_S6B_E COVID_S6B_E. ;
format COVID_S6B_G COVID_S6B_G. ;
format COVID_S6B_H COVID_S6B_H. ;
format COVID_S6B_I COVID_S6B_I. ;
format COVID_S6C_A COVID_S6C_A. ;
format COVID_S6C_B COVID_S6C_B. ;
format COVID_S6C_C COVID_S6C_C. ;
format COVID_S6C_D COVID_S6C_D. ;
format COVID_S6C_E COVID_S6C_E. ;
format COVID_S6C_F COVID_S6C_F. ;
format COVID_S6C_G COVID_S6C_G. ;
format COVID_S6C_H COVID_S6C_H. ;
format COVID_S6C_I COVID_S6C_I. ;
format PROBE_MASKTYPE_A PROBE_MASKTYPE_A. ;
format PROBE_MASKTYPE_B PROBE_MASKTYPE_B. ;
format PROBE_MASKTYPE_C PROBE_MASKTYPE_C. ;
format PROBE_MASKTYPE_D PROBE_MASKTYPE_D. ;
format PROBE_MASKTYPE_E PROBE_MASKTYPE_E. ;
format PROBE_MASKTYPE_F PROBE_MASKTYPE_F. ;
format PROBE_MASKUSE PROBE_MASKUSE. ;
format COVID_IMPACT1 COVID_IMPACT1x. ;
format COVID_IMPACT2 COVID_IMPACT2x. ;
format COVID_IMPACT3 COVID_IMPACT3x. ;
format COVID_IMPACT4 COVID_IMPACT4x. ;
format GAD7_A GAD7_A. ;
format GAD7_B GAD7_B. ;
format PHQ_A PHQ_A. ;
format PHQ_B PHQ_B. ;
format VISIONDF VISIONDF. ;
format HEARINGDF HEARINGDF. ;
format DIFF DIFF. ;
format COMDIFF COMDIFF. ;
format COGMEMDIFF COGMEMDIFF. ;
format UPPSLFCR UPPSLFCR. ;
format ANXFREQ ANXFREQ. ;
format ANXMED ANXMED. ;
format ANXLEVEL ANXLEVEL. ;
format DEPFREQ DEPFREQ. ;
format DEPMED DEPMED. ;
format DEPLEVEL DEPLEVEL. ;
format PROBE_ANX_1 PROBE_ANX_1x. ;
format PROBE_ANX_2 PROBE_ANX_2x. ;
format PROBE_ANX_3 PROBE_ANX_3x. ;
format PROBE_ANX_4 PROBE_ANX_4x. ;
format PROBE_ANX_5 PROBE_ANX_5x. ;
format PROBE_ANX_6 PROBE_ANX_6x. ;
format PROBE_ANX_7 PROBE_ANX_7x. ;
format PROBE_ANX_DK PROBE_ANX_DK. ;
format PROBE_ANX_REF PROBE_ANX_REF. ;
format PROBE_DEP_1 PROBE_DEP_1x. ;
format PROBE_DEP_2 PROBE_DEP_2x. ;
format PROBE_DEP_3 PROBE_DEP_3x. ;
format PROBE_DEP_4 PROBE_DEP_4x. ;
format PROBE_DEP_5 PROBE_DEP_5x. ;
format PROBE_DEP_6 PROBE_DEP_6x. ;
format PROBE_DEP_7 PROBE_DEP_7x. ;
format PROBE_DEP_DK PROBE_DEP_DK. ;
format PROBE_DEP_REF PROBE_DEP_REF. ;
format REL1 REL1x. ;
format REL2 REL2x. ;
format REL3 REL3x. ;
format REL4 REL4x. ;
format REL5 REL5x. ;
format REL6 REL6x. ;
format REL7 REL7x. ;
format REL_TOTALTIME REL_TOTALTIME. ;
format BURDEN1 BURDEN1x. ;
format BURDEN2 BURDEN2x. ;
format SURV_MODE SURV_MODE. ;
format GENDER GENDER. ;
format RACETHNICITY RACETHNICITY. ;
format EDUC EDUC. ;
format MARITAL MARITAL. ;
format EMPLOY EMPLOY. ;
format INCOME INCOME. ;
format REGION4 REGION4x. ;
format INTERNET INTERNET. ;
format HOUSING HOUSING. ;
format HOME_TYPE HOME_TYPE. ;
format PHONESERVICE PHONESERVICE. ;
*/
run;

/* check the original data */

proc contents data = d.RANDS_COVID_3;
run;

/* extract a subset variable */
data rands_covid3_ori;
  set d.RANDS_COVID_3;
    keep CaseID S_VPSU      S_VSTRAT     WEIGHT     	WEIGHT_CALIBRATED 
	AGE				EDUC        EMPLOY     	 GENDER 
    HHSIZE      	HICOV       HOME_TYPE    HOUSING    	HYPEV 
    INCOME      	INTERNET    MARITAL      PHONESERVICE   PHSTAT 
    RACETHNICITY    REGION4     EMPLASTWK  SMKEV        SMKNOW         
    USUALPL ;
	/* change code 98 to missingness */
	if HICOV=98 then HICOV=.;
	if HYPEV=98 then HYPEV=.;
	if PHSTAT=98 then PHSTAT=.;
	if EMPLASTWK=98 then EMPLASTWK=.;
	if SMKEV=98 then SMKEV=.;
	if SMKNOW=98 then SMKNOW=.;
	if USUALPL=98 then USUALPL=.;
run;


/* the missingness rate */

proc freq data=rands_covid3_ori;
  tables /* CaseID */ S_VPSU      S_VSTRAT /* WEIGHT WEIGHT_CALIBRATED */      
	AGE				EDUC        EMPLOY     	 GENDER 
    HHSIZE      	HICOV       HOME_TYPE    HOUSING    	HYPEV 
    INCOME      	INTERNET    MARITAL      PHONESERVICE   PHSTAT 
    RACETHNICITY    REGION4     EMPLASTWK  SMKEV        SMKNOW         
    USUALPL / missing;
run;

proc univariate data = rands_covid3_ori normal plot;
  var weight_calibrated;
  histogram weight_calibrated;
run;


/* strata and psu */
proc freq data = rands_covid3_ori;
  tables S_VSTRAT S_VPSU;
run;

/* create a new variable that has the frequency of S_VSTRAT */
proc freq data = rands_covid3_ori;
  tables S_VSTRAT;
  ods output OneWayFreqs = strata_freq (keep = S_VSTRAT Frequency);
run;



/* merge the strata freq to the original data */
proc sort data = rands_covid3_ori;
  by S_VSTRAT;
run;

proc sort data = strata_freq;
  by S_VSTRAT;
run;

data rands_covid3_ori;
  merge rands_covid3_ori strata_freq;
  by S_VSTRAT;
  S_VSTRAT_COMBINE = S_VSTRAT;
  if Frequency < 10 then S_VSTRAT_COMBINE = 99;
  drop Frequency;
run;

proc freq data = rands_covid3_ori;
  tables S_VSTRAT S_VSTRAT_COMBINE / missing;
run;

/* focus on the income and marital */
/* Income range from 1 to 16 */
proc freq data= rands_covid3_ori;
  tables INCOME MARITAL / missing;
run;

/* exploratory analysis of income */
proc univariate data = rands_covid3_ori normal plot;
  var INCOME;
  histogram INCOME;
run;

/* survey analysis of income */
/* average 10.38, 95% CI is (10.14, 10.62) */
proc surveymeans data=rands_covid3_ori;
weight WEIGHT_CALIBRATED;
strata S_VSTRAT;
cluster S_VPSU;
var INCOME;
run;

/* The relationship between income and weights */
/* a slight positive corrleation between calibrated weights and income */
proc corr data = rands_covid3_ori;
  var INCOME WEIGHT_CALIBRATED;
  * var INCOME WEIGHT;
run;

proc plot data = rands_covid3_ori;
  plot INCOME * WEIGHT_CALIBRATED = "*";
run;

/* sampling strata */
proc sort data = rands_covid3_ori;
  by S_VSTRAT;
run;

proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by S_VSTRAT;
run;

/* age variable */
proc means data = rands_covid3_ori;
   var AGE;
   weight WEIGHT_CALIBRATED;
run;

/* higher age with lower income */
proc corr data = rands_covid3_ori;
  var INCOME AGE;
 * weight WEIGHT_CALIBRATED;
run;

proc plot data = rands_covid3_ori;
  plot INCOME * AGE = "*";
run;

/* gender */
/* male have higher income than female */
proc sort data = rands_covid3_ori;
  by GENDER;
run;

proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by GENDER;
run;

/* education */
/* higher education has higher income */
proc sort data = rands_covid3_ori;
  by EDUC;
run;

proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by EDUC;
run;

/* race */
/* White > other > Hispanic > Black */
proc sort data = rands_covid3_ori;
  by RACETHNICITY;
run;

proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by RACETHNICITY;
run;

/* census region */
/* Northeast > West > Midwest > South */
proc sort data = rands_covid3_ori;
  by REGION4;
run;

proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by REGION4;
run;

/* current employment status */
/* 1 > 2 > 7 >5 > 3 > 4 > 6 */
proc sort data = rands_covid3_ori;
  by EMPLOY;
run;

proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by EMPLOY;
run;

/* Employment last week */
/* 1 > 2 */
proc sort data = rands_covid3_ori;
  by EMPLASTWK;
run;

proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by EMPLASTWK;
run;

/* household size */
/* 5 > 4 > 2 >3 > 6 > 1 */
proc sort data=rands_covid3_ori;
  by HHSIZE;
run;

proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by HHSIZE;
run;

/* Health insurance coverage */
/* Yes > No */
proc sort data = rands_covid3_ori;
  by HICOV;
run;

proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by HICOV;
run;

/* Home type */
/* 1 > 2 > 3 > 4 */
proc sort data = rands_covid3_ori;
  by HOME_TYPE;
run;

proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by HOME_TYPE;
run;

/* Housing */
/* 1 > 2 > 3 */
proc sort data = rands_covid3_ori;
  by HOUSING;
run;

proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by HOUSING;
run;

/* Hypertension */
/* Yes < No */
proc sort data = rands_covid3_ori;
  by HYPEV;
 run;

proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by HYPEV;
run;

/* Internet */
/* No < Yes */
proc sort data = rands_covid3_ori;
  by INTERNET;
run;

proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by INTERNET;
run;

/* phone service */
/* 2 > 3 > 5 > 4 > 1 */
proc sort data = rands_covid3_ori;
  by PHONESERVICE;
run;

proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by PHONESERVICE;
run;

/* self-rated health status */
/* 2 > 1 > 3 > 4 > 5 */
/* 1 and 2 vs. 3 - 5 */
proc sort data = rands_covid3_ori;
  by PHSTAT;
run;

/* ever smoking */
/* Yes < No */
proc sort data = rands_covid3_ori;
  by SMKEV;
run;

proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by SMKEV;
run;


proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by PHSTAT;
run;

/* usual place of care */
/* 1 > 3 > 2 */
proc sort data = rands_covid3_ori;
  by USUALPL;
run;

proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by USUALPL;
run;

/* marital status */
/* 1 > 6 > 5 > 3 > 2 > 4 */
proc sort data = rands_covid3_ori;
  by MARITAL;
run;

proc means data = rands_covid3_ori;
  var INCOME;
  weight WEIGHT_CALIBRATED;
  by MARITAL;
run;



/* create a new binary marital status variable */
data rands_covid3_new;
  set rands_covid3_ori;
  if MARITAL = 1 or MARITAL = 6 then MARITAL_NEW = 1;
  else MARITAL_NEW =0;
run;

proc freq data = rands_covid3_new;
  tables MARITAL MARITAL_NEW / missing;
run;

/* using a regression to find the predictors of income */
proc surveyreg data = rands_covid3_new;
  weight WEIGHT_CALIBRATED;
  strata S_VSTRAT;
  cluster S_VPSU;
  /* treat HHSIZE as a continuous variable */
  class EDUC GENDER INTERNET MARITAL_NEW
         EMPLOY HICOV HOME_TYPE HOUSING HYPEV 
         PHONESERVICE PHSTAT RACETHNICITY REGION4 EMPLASTWK SMKEV MARITAL_NEW USUALPL;
  model INCOME = AGE EDUC GENDER INTERNET MARITAL_NEW 	
		EMPLOY HHSIZE HICOV HOME_TYPE HOUSING HYPEV 
        PHONESERVICE PHSTAT RACETHNICITY REGION4 EMPLASTWK SMKEV USUALPL / solution;
run;

/* exploratory analysis for marital status (marital_new) */
/* overall mean */
/* the mean is 0.613 (0.592, 0.634) */
proc surveymeans data = rands_covid3_new;
  var MARITAL_NEW;
  strata S_VSTRAT;
  cluster S_VPSU;
  weight WEIGHT_CALIBRATED;
run;

/* relationship with other variables */
proc surveylogistic data = rands_covid3_new;
  weight WEIGHT_CALIBRATED;
  strata S_VSTRAT;
  cluster S_VPSU;
  class EDUC GENDER INTERNET MARITAL_NEW
         EMPLOY HICOV HOME_TYPE HOUSING HYPEV 
         PHONESERVICE PHSTAT RACETHNICITY REGION4 EMPLASTWK SMKEV USUALPL;
  model MARITAL_NEW (event = '1') = INCOME AGE EDUC GENDER INTERNET 	
		EMPLOY HHSIZE HICOV HOME_TYPE HOUSING HYPEV 
        PHONESERVICE PHSTAT RACETHNICITY REGION4 EMPLASTWK SMKEV USUALPL;
run;


/* generate missing data for marital status and INCOME */
data rands_covid_new;
 * call streaminit(197789);
  set rands_covid3_new;
  p_miss_INCOME = exp(-2+0.5*EDUC-0.5*GENDER-0.01*AGE+0.5*INTERNET)/(1+exp(-2+0.5*EDUC-0.5*GENDER-0.01*AGE+0.5*INTERNET));
  rnumber_INCOME = ranuni(20110411);
  If rnumber_INCOME < p_miss_INCOME then R_miss_INCOME =1; 
  else R_miss_INCOME=0;  
  
  /* If R_miss_INCOME = 1 then INCOME=.; */

  /* generate missing data for marital status depending on age, education, household size */
  p_miss_MARITAL = exp(-2-0.5*EDUC+0.75*HHSIZE+0.01*AGE-0.75*INTERNET)/(1+exp(-2-0.5*EDUC+0.75*HHSIZE+0.01*AGE-0.75*INTERNET));

  rnumber_MARITAL = ranuni(20140201);
  If rnumber_MARITAL < p_miss_MARITAL then R_miss_MARITAL =1; 
  else R_miss_MARITAL=0;  

 /* concatenate the strata and psu together */
  str_psu= S_VPSU||S_VSTRAT;

run;
/*
proc freq data = rands_covid_new;
  tables str_psu / missing;
run;
*/
proc freq data = rands_covid_new;
  table R_miss_INCOME R_miss_MARITAL / missing;
run;

proc means data = rands_covid_new;
  var p_miss_INCOME p_miss_MARITAL;
run;


/* the observed income means is 10.17 (9.91, 10.43) */
proc sort data = rands_covid_new;
  by R_miss_INCOME;
run;

proc surveymeans data = rands_covid_new;
  var INCOME;
  by R_miss_INCOME;
  weight WEIGHT_CALIBRATED;
  strata S_VSTRAT;
  cluster S_VPSU;
run;

/* the observed marital status is 0.589 (0.565, 0.614) */ 
proc sort data = rands_covid_new;
  by R_miss_MARITAL;
run;

proc surveymeans data = rands_covid_new;
  var MARITAL_NEW;
  by R_miss_MARITAL;
  weight WEIGHT_CALIBRATED;
  strata S_VSTRAT;
  cluster S_VPSU;
run;

/* set up the missing data */
data rands_covid_missing;
  set rands_covid_new;
  If R_miss_INCOME = 1 then INCOME=.;
  If R_miss_MARITAL = 1 then MARITAL_NEW=.;
run;

proc freq data = rands_covid_missing;
  tables INCOME MARITAL_NEW / missing;
run;


/* The propensity score model for the nonresponse */

/* using survey logistic */
proc surveylogistic data = rands_covid_missing;
  weight WEIGHT_CALIBRATED;
  strata S_VSTRAT      ;
  cluster S_VPSU;
  model R_miss_income (event = "1") = EDUC GENDER AGE INTERNET;
run;

proc surveylogistic data = rands_covid_missing;
  weight WEIGHT_CALIBRATED;
  strata S_VSTRAT      ;
  cluster S_VPSU;
  model R_miss_MARITAL (event = "1") = EDUC HHSIZE AGE INTERNET;
run;


/* The chi-square test for the nonresponse indicator */
/*
proc surveyfreq data=rands_covid_missing;
weight WEIGHT_CALIBRATED;
strata S_VSTRAT      ;
cluster S_VPSU;
table R_miss_income*(EDUC GENDER INTERNET)/chisq;
run;
*/

/* set up the imputation model for the income variable */
/* the simplest model includes age educ, gender, internet, marital_new, hhsize and combined strata */
proc mi data =rands_covid_missing seed =197789 out= income_impute  nimpute =5 
    min = 1 . . . . . . . . 
    max = 16 . . . . . . . .   
;
   class EDUC GENDER INTERNET MARITAL_NEW S_VSTRAT_COMBINE ; 
   fcs nbiter=20 reg (INCOME/details) logistic (MARITAL_NEW / details likelihood=augment) ; 
   * fcs nbiter=20 regpmm (INCOME/details) logistic (MARITAL_NEW / details likelihood=augment);
   * fcs nbiter=20 reg (INCOME/details) discrim ( MARITAL_NEW   /classeffects =include details ); 
   * fcs nbiter=20 regpmm (INCOME/details) discrim ( MARITAL_NEW   /classeffects =include details );  
   var INCOME AGE WEIGHT_CALIBRATED EDUC GENDER INTERNET MARITAL_NEW HHSIZE	
     S_VSTRAT_COMBINE;  
run;

/* the complicated model includes more variables */
proc mi data =rands_covid_missing seed =197789 out= income_impute  nimpute =20 
  min = 1 . . . . . . . . . . . . . . . . . . . . 
  max = 16 . . . . . . . . . . . . . . . . . . . .   
;
   class EDUC GENDER INTERNET MARITAL_NEW S_VSTRAT_COMBINE
         EMPLOY HICOV HOME_TYPE HOUSING HYPEV 
         PHONESERVICE PHSTAT RACETHNICITY REGION4 EMPLASTWK SMKEV USUALPL; 
   * fcs nbiter=20 reg (INCOME/details) logistic (MARITAL_NEW / details likelihood = augment);
   * fcs nbiter=20 regpmm (INCOME/details) logistic (MARITAL_NEW /details likelihood = augment);
   * fcs nbiter=20 reg (INCOME/details) discrim ( MARITAL_NEW   /classeffects =include details ); 
     fcs nbiter=20 regpmm (INCOME/details) discrim ( MARITAL_NEW   /classeffects =include details ); 
var  INCOME AGE WEIGHT_CALIBRATED EDUC GENDER INTERNET MARITAL_NEW	
     	S_VSTRAT_COMBINE
		EMPLOY HHSIZE HICOV HOME_TYPE HOUSING HYPEV 
        PHONESERVICE PHSTAT RACETHNICITY REGION4 EMPLASTWK SMKEV USUALPL;  
run;


/* however some of the imputed values are out of bound */
proc univariate data = income_impute normal plot;
  var INCOME;
run;

/* post processing it */
data income_impute;
  set income_impute;
  INCOME_new=round(INCOME);
  if INCOME_NEW < 1 then INCOME_NEW = 1;
  if INCOME_NEW > 16 then INCOME_NEW = 16;
run;

/*
proc univariate data = income_impute normal plot;
  var INCOME INCOME_NEW;
  histogram INCOME INCOME_NEW;
run; 

proc freq data = income_impute;
  tables INCOME INCOME_NEW / missing;
run;

*/
/*income MIanalyze*/
proc surveymeans data=income_impute;
weight WEIGHT_CALIBRATED;
strata S_VSTRAT;
cluster S_VPSU;
var INCOME INCOME_NEW MARITAL_NEW;
by   _imputation_;;
ods output  Statistics = mean_income_imp;
run;
/*
proc print data = mean_income_imp;
run;
*/
/* for the original imputed values */
/* the simple model reg + logistic results is 10.33 (10.06, 10.59) */
/* the simple model reg + discriminant is 10.36 (10.11, 10.60) */
/* the simple model regpmm + logistic results is 10.34 (10.11, 10.58) */
/* the simple model regpmm + discrimiant is 10.31 (10.06, 10.56) */

/* the complicated model reg + logistic results is 10.40 (10.16, 10.63) */
/* the complicated model regpmm + logistic results is 10.37 (10.13, 10.61) */
/* the complicated model reg + discriminant results is 10.39 (10.15, 10.62) */
/* the complicated model regpmm + discriminant results is 10.37 (10.14, 10.61) */


/* adding constraint for the income, simple normal model + logistic results is 10.18 (9.94, 10.41) */
/* adding constraint for the income, simple normal model + discriminant anylsis is 10.22 (9.99, 10.45) */
/* adding constraint for the income, complicated normal model + logistic results is 10.26 (10.03, 10.49) */
/* adding constraint for the income, complicated normal model + discriminant results is 10.28 (10.05, 10.52) */
proc MIanalyze data =mean_income_imp edf=88;
modeleffects mean;
StdErr  StdErr;
where VarName = 'INCOME';
ods output ParameterEstimates=MI_results_income;
run;

/* for the rounded values */
/* the simple model reg + logistic results is 10.28 (10.03, 10.53) */
/* the simple model reg + discrimiant results is 10.31 (10.07, 10.55) */


/* the complicated reg + logistic results is 10.35 (10.12, 10.58) */
/* the complicated reg + discriminant results is 10.34 (10.11, 10.57) */

/* adding constraint for the income, simple normal model + logistic results is 10.18 (9.94, 10.41) */
/* adding constraint for the income, simple normal model + discriminant analysis is 10.22 (9.99, 10.45) */
/* adding constraint for the income, complicated normal model + logistic results is 10.26 (10.03, 10.49) */
/* adding constraint for the income, complicated normal model + discriminant results is 10.28 (10.05, 10.52) */
proc MIanalyze data =mean_income_imp edf=88;
modeleffects mean;
StdErr  StdErr;
where VarName = 'INCOME_new';
ods output ParameterEstimates=MI_results_income_round;
run;

/* for the marital status */
/* the simple model reg + logistic results is 0.635 (0.611, 0.658) */
/* the simple model reg + discriminant results is 0.638 (0.615, 0.661) */
/* the simple model regpmm + logistic results is 0.639 (0.615, 0.662) */
/* the simple model regpmm + discriminant results is 0.636 (0.614, 0.658) */

/* the complicated model reg + logistic results is 0.620 (0.598, 0.642) */
/* the complicated model regpmm + logistic results is 0.619 (0.594, 0.645) */
/* the complicated model reg + discriminant is 0.621 (0.598, 0.645) */
/* the complicated model regpmm + discriminant is 0.620 (0.594, 0.646) */

/* adding constraint for the income, simple normal model + logistic results is 0.638 (0.609, 0.667) */
/* adding constraint for the income, simple normal model + discriminant anaylsis 0.638 (0.615, 0.660) */
/* adding constraint for the income, complicated normal model + logistic results is 0.618 (0.595, 0.642) */
/* adding constraint for the income, complicated normal model + discrimiant results is 0.619 (0.596, 0.643) */
proc MIanalyze data =mean_income_imp edf=88;
modeleffects mean;
StdErr  StdErr;
where VarName = 'MARITAL_NEW';
ods output ParameterEstimates=MI_results_marital;
run;


/* a small test of data synthesis */
data synthesis_test;
  set rands_covid3_ori;
  ID = 1;
  keep INCOME AGE WEIGHT_CALIBRATED EDUC GENDER INTERNET MARITAL HHSIZE	
     S_VSTRAT_COMBINE ID;
run;

data synthesis_missing;
  set synthesis_test;
  /* set all variables as missing except the weight and strata */
  INCOME=.;
  AGE=.;
  EDUC=.;
  GENDER=.;
  INTERNET=.;
  MARITAL=.;
  HHSIZE=.;
  ID = 0;
run;

/* concatenate the two datasets */
data synthesis_final;
  set synthesis_test synthesis_missing;
run;

proc freq data = synthesis_final;
  tables ID / missing;
run;

/* multiple imputation for data synthesis */
proc mi data =synthesis_final seed =197789 out= synthesis_final_impute  nimpute =50  ;
   class EDUC GENDER INTERNET MARITAL S_VSTRAT_COMBINE ; 
    * fcs nbiter=20 reg (INCOME/details) logistic (MARITAL_NEW / details likelihood=augment) ; 
     fcs nbiter=100 reg (INCOME AGE HHSIZE /details) logistic (EDUC GENDER INTERNET MARITAL / details likelihood=augment);
   * fcs nbiter=20 reg (INCOME/details) discrim ( MARITAL_NEW   /classeffects =include details ); 
   * fcs nbiter=20 regpmm (INCOME/details) discrim ( MARITAL_NEW   /classeffects =include details );  
   var INCOME AGE WEIGHT_CALIBRATED EDUC GENDER INTERNET MARITAL HHSIZE	
     S_VSTRAT_COMBINE;  
run;

/* data analysis */
proc surveymeans data=synthesis_final_impute;
weight WEIGHT_CALIBRATED;
strata S_VSTRAT_COMBINE;
var INCOME AGE EDUC GENDER INTERNET MARITAL HHSIZE;
where ID = 0;
by   _imputation_;
ods output  Statistics = result_synthesis;
run;


proc print data = result_synthesis;
run;

proc MIanalyze data = result_synthesis edf=88;
modeleffects mean;
StdErr  StdErr;
where VarName = 'INCOME';
ods output ParameterEstimates= synthesis_results_final;
run;

/* original data results */
proc surveymeans data = synthesis_test;
  strata S_VSTRAT_COMBINE;
  weight WEIGHT_CALIBRATED;
  var INCOME;
run;

proc MIanalyze data = result_synthesis edf=88;
modeleffects mean;
StdErr  StdErr;
where VarName = 'HHSIZE';
ods output ParameterEstimates= synthesis_results_final;
run;

/* original data results */
proc surveymeans data = synthesis_test;
  strata S_VSTRAT_COMBINE;
  weight WEIGHT_CALIBRATED;
  var HHSIZE;
run;

