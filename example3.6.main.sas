*********************************************************************
Comining CHISQ TEST 
*********************************************************************;

* USER NOTE: REPLACE CORRECT PATH AND PARAMETERS BEFORE EXECUTING THE SAS PROGRAM;

%LET YEAR = 2016;       *** provide survey year ***;

  *** path to store the imputed income SAS datasets ***;
* LIBNAME  NHIS   "C:\Guangyu\Yulei\nhis2016";
* LIBNAME  LIBRARY   "C:\Guangyu\Yulei\nhis2016";
* DEFINE VARIABLE VALUES FOR REPORTS;

*  USE THE STATEMENT "PROC FORMAT LIBRARY=LIBRARY"
     TO PERMANENTLY STORE THE FORMAT DEFINITIONS;

*  USE THE STATEMENT "PROC FORMAT" IF YOU DO NOT WISH
      TO PERMANENTLY STORE THE FORMATS.;

**PROC FORMAT LIBRARY=LIBRARY;

LIBNAME  NHIS   "\\cdc.gov\private\A728\wdq7\imputation_diagnostics\data";
LIBNAME  LIBRARY   "\\cdc.gov\private\A728\wdq7\imputation_diagnostics\data";


PROC FORMAT;

   VALUE INC000X
      25				= "25 Income Imputation"
   ;

   VALUE INC001X
   	  0					 = "0 Reported"
	  1					 = "1 Imputed; no information"
	  2					 = "2 Imputed; reported in categories"
   ;

   VALUE INC002X
      0					= "0 Not top-coded"
	  1					= "1 Top-coded"
   ;

   VALUE INC003X
      0					= "0 Not imputed"
	  1					= "1 Imputed"
   ;

   VALUE INC004X
      1					= "1 Employed"
	  2					= "2 Not employed"
   ;

   VALUE INC005X
      0                  = "0 Reported"
      1                  = "1 Imputed"
   ;
RUN;

%macro allimp;
%do IMPNUM = 1 %to 5;

  *** path to the imputed income ASCII datasets ***;
* FILENAME  ASCIIDAT  "C:\Guangyu\Yulei\nhis2016\INCMIMP&IMPNUM..DAT";

FILENAME  ASCIIDAT  "\\cdc.gov\private\A728\wdq7\imputation_diagnostics\data\INCMIMP&IMPNUM..DAT";


DATA NHIS.INCMIMP&IMPNUM;   *** CREATE A SAS DATA SET ***;

   INFILE ASCIIDAT PAD LRECL=44;

   * DEFINE LENGTH OF ALL VARIABLES;

   LENGTH
		RECTYPE 3		SRVY_YR 4		HHX $6			FMX $2
		FPX $2			IMPNUM 3		FAMINCF2 3		TCINCM_F 3
		FAMINCI2 8		POVRATI3 8		EMPLOY_F 3		EMPLOY_I 3
		ERNYR_F 3		TCEARN_F 3		ERNYR_I2 8
		;

   * INPUT ALL VARIABLES;
   * IMPLIED DECIMAL IS FORMALLY PLACED IN THE APPROPRIATE LOCATION FOR THE VARIABLE POVRATI3;

   INPUT
		RECTYPE 	1-2			SRVY_YR 	3-6		
		HHX 		7-12		FMX 		13-14
		FPX 		15-16		IMPNUM 		17 		
		FAMINCF2	18			TCINCM_F 	19
		FAMINCI2 	20-25		POVRATI3 	26-34 .3
		EMPLOY_F 	35 			EMPLOY_I 	36 			
		ERNYR_F 	37 			TCEARN_F 	38
		ERNYR_I2 	39-44
		;

   * DEFINE VARIABLE LABELS;

   LABEL
		RECTYPE		= "File identifier type"
		SRVY_YR		= "Year of National Health Interview Survey"
		HHX			= "HH identifier"
		FMX			= "Family identifier"
		FPX			= "Person number identifier"
		IMPNUM		= "Imputation number"
		FAMINCF2	= "Family income imputation flag"
		TCINCM_F	= "Family income/poverty ratio top-coded flag"
		FAMINCI2	= "Top-coded family income"
		POVRATI3	= "Ratio of family income to poverty threshold"
		EMPLOY_F	= "Employment status imputation flag"
		EMPLOY_I	= "Person's employment status"
		ERNYR_F		= "Person's earnings imputation flag"
		TCEARN_F	= "Person's earnings top-coded flag"
		ERNYR_I2	= "Person's total earnings last year (top-coded)"
		;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   FORMAT
		RECTYPE		INC000X.
		FAMINCF2	INC001X.
		TCINCM_F	INC002X.
		EMPLOY_F	INC003X.
		EMPLOY_I	INC004X.
		ERNYR_F		INC005X.
		TCEARN_F	INC002X.
		;
RUN;

PROC CONTENTS DATA=NHIS.INCMIMP&IMPNUM;
   TITLE1 "CONTENTS OF THE &YEAR NHIS IMPUTED INCOME FILE, DATASET &IMPNUM";
RUN;
PROC FREQ DATA=NHIS.INCMIMP&IMPNUM;

   TABLES   RECTYPE		SRVY_YR		IMPNUM		FAMINCF2
			TCINCM_F	EMPLOY_F	EMPLOY_I	ERNYR_F
			TCEARN_F
			;
   TITLE1 "FREQUENCY REPORT FOR &YEAR NHIS IMPUTED INCOME FILE, DATASET &IMPNUM";
   TITLE2 '(UNWEIGHTED)';
PROC MEANS DATA=NHIS.INCMIMP&IMPNUM;

	VAR		FAMINCI2	POVRATI3	ERNYR_I2
			;
   TITLE1 "MEANS FOR &YEAR NHIS IMPUTED INCOME FILE, DATASET &IMPNUM";
   TITLE2 '(UNWEIGHTED)';

* USER NOTE: TO SEE UNFORMATTED VALUES IN PROCEDURES, ADD THE
             STATEMENT: FORMAT _ALL_;
RUN;
%end;
%mend allimp;
%allimp;


%macro incomedata;
%do i =1 %to 5;
data income&i;
set NHIS.INCMIMP&i;
if FAMINCF2 =0 then income_ori=FAMINCI2;
else if FAMINCF2  in ( 1 2 ) then income_ori=.;
HHX_FMX=HHX||FMX;
if income_ori= . then income_M=1;
else income_M=0;

if 0<=FAMINCI2<35000 then income_grp_imp=1;
else if 35000<=FAMINCI2<75000 then income_grp_imp=2;
else if 75000<=FAMINCI2<100000 then income_grp_imp=3;
else if 100000<=FAMINCI2  then income_grp_imp=4;

keep HHX FMX  HHX_FMX income_ori income_M FAMINCI2 FAMINCF2 income_grp_imp IMPNUM;
run;


proc sort data=income&i out=income&i nodupkey;
by HHX_FMX;
run;
proc sort data=income&i;
by HHX FMX ;
run;
%end;
%mend;
%incomedata;


data FAMILYXX;
set NHIS.FAMILYXX;

if FM_STRCP =99 then FM_STRCP=.;
if FWICYN in ( 7 8 9) then FWICYN=.;
if FSNAP in ( 7 8 9) then FSNAP=.;
if FGAH in ( 7 8 9) then FGAH=.;
if Houseown in ( 7 8 9) then Houseown=.;
if FMEDBILL in ( 7 8 9) then FMEDBILL=.;
if CURWRKN in ( 7 8 9) then CURWRKN=.;
if TELCELN in ( 7 8 9) then TELCELN=.;
if FSRUNOUT in (7 8 9) then FSRUNOUT=.;
if FSALYN in (7 8 9) then FSALYN=.;
if FDIVDYN in (7 8 9) then FDIVDYN=.;
if FM_EDUC1 in (01 02 03 04) then do;
FM_EDUC1_Recode= 1 ;/*high school and less*/
educ_1_dummy1=1;
educ_1_dummy2=0;
end;

else if FM_EDUC1 in (05 06 07) then do;
FM_EDUC1_Recode= 2 ;/*some college*/
educ_1_dummy1=0;
educ_1_dummy2=1;
end;

else if FM_EDUC1 in (08 09) then do;
FM_EDUC1_Recode= 3 ; /*college and graduate*/
educ_1_dummy1=0;
educ_1_dummy2=0;
end;

else if FM_EDUC1 in (97 98 99) then do;
FM_EDUC1_Recode=.;
educ_1_dummy1=.;
educ_1_dummy2=.;
end;

Familyind=1;
keep HHX FMX FM_EDUC1_Recode Familyind WTFA_FAM;
run;


proc sort data=FAMILYXX ;
by HHX FMX ;
run;

data household;
set NHIS.HOUSEHLD;
keep HHX PSTRAT PPSU;
run;

proc sort data=household;
by HHX;
run;


%macro imputed_all;
%do i=1 %to 5;

data family_w_income&i;
merge income&i FAMILYXX;
by HHX FMX ;
run;

data family_w_income&i;
merge family_w_income&i household;
by HHX ;
if Familyind=1;
run;

%end;
%mend;
%imputed_all;

Data all_income;
set family_w_income1 -family_w_income5;
run;


proc sort data=all_income;
by IMPNUM;
run;

/* wald chi-square */
/*
proc surveyfreq data=all_income;
table income_grp_imp*FM_EDUC1_Recode/row WCHISQ ;
weight WTFA_FAM;
strata PSTRAT;
cluster PPSU;
by IMPNUM;
ods output CrossTabs=out_percent;
run;
*/

/* rao scott chi-square default */
/* For Table 3.5 */
proc surveyfreq data=all_income;
tables income_grp_imp * FM_EDUC1_Recode /chisq;
weight WTFA_FAM;
cluster PPSU;
strata PSTRAT;
by IMPNUM;
* ods output CrossTabls=out_percent1;
run;

/* rao scott chi-square default */
/*
proc surveyfreq data=all_income;
tables income_grp_imp * FM_EDUC1_Recode /chisq (firstorder);
weight WTFA_FAM;
cluster PPSU;
strata PSTRAT;
by IMPNUM;
* ods output CrossTabls=out_percent1;
run;
*/

/*
proc surveyfreq data=all_income;
tables income_grp_imp * FM_EDUC1_Recode /chisq (modified);
weight WTFA_FAM;
cluster PPSU;
strata PSTRAT;
by IMPNUM;
* ods output CrossTabls=out_percent1;
run;
*/

/* rao scott chi-square second order */
/*
proc surveyfreq data=all_income;
tables income_grp_imp * FM_EDUC1_Recode /chisq (secondorder);
weight WTFA_FAM;
cluster PPSU;
strata PSTRAT;
by IMPNUM;
* ods output CrossTabls=out_percent1;
run;
*/

/* SAS macro combining chi-square statistics from multiply imputed datasets */
/* we use the default Rao Scott Chi-square test statistics for survey data */
/* For Table 3.5 */
%macro combchi(df=,chi=);
proc iml;
  df=&df;
  g2={&chi};
  m=ncol(g2);
  g=sqrt(g2);
  mg2=sum(g2)/m;
  r=(1+1/m)*(ssq(g)-(sum(g)**2)/m)/(m-1);
  f=(mg2/df - r*(m-1)/(m+1))/(1+r);
  ddf=(m-1)*(1+1/r)**2/df**(3/m);
  p=1-probf(f,df,ddf);
  print f df ddf;
  print p;
run;
%mend combchi;
%combchi(df=6, chi=8837.57 9192.08 9103.73 9211.00 8815.56);

