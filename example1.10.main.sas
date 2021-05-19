*********************************************************************
 July 20, 2017

 THIS IS AN EXAMPLE OF A SAS PROGRAM THAT CREATES A SAS
 FILE FROM THE PUBLIC USE IMPUTED INCOME ASCII FILES

 THIS IS STORED IN INCMIMP.SAS
*********************************************************************;

* USER NOTE: REPLACE CORRECT PATH AND PARAMETERS BEFORE EXECUTING THE SAS PROGRAM;

%LET YEAR = 2016;       *** provide survey year ***;

  *** path to store the imputed income SAS datasets ***;
LIBNAME  NHIS   "\\cdc.gov\private\A728\wdq7\imputation_diagnostics\data";
* LIBNAME  LIBRARY   "C:\Guangyu\Yulei\nhis2016";
* DEFINE VARIABLE VALUES FOR REPORTS;

*  USE THE STATEMENT "PROC FORMAT LIBRARY=LIBRARY"
     TO PERMANENTLY STORE THE FORMAT DEFINITIONS;

*  USE THE STATEMENT "PROC FORMAT" IF YOU DO NOT WISH
      TO PERMANENTLY STORE THE FORMATS.;

**PROC FORMAT LIBRARY=LIBRARY;
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


data income1;
 set NHIS.INCMIMP1;
* set NHIS.INCMIMP2;
* set NHIS.INCMIMP3;
* set NHIS.INCMIMP4;
* set NHIS.INCMIMP5;
if FAMINCF2 =0 then income_ori=FAMINCI2;
else if FAMINCF2  in ( 1 2 ) then income_ori=.;
HHX_FMX=HHX||FMX;
if income_ori= . then income_M=1;
else income_M=0;

keep HHX FMX  HHX_FMX income_ori income_M FAMINCI2 FAMINCF2;
run;


proc sort data=income1 out=income1 nodupkey;
by HHX_FMX;
run;
proc sort data=income1;
by HHX FMX ;
run;

/*
proc contents data=income1;
run;

proc univariate data=income1 normal plot ;
  var faminci2;
run;

proc freq data=income1;
  table famincf2 income_M / missing;
run;
*/

/* save dataset for imputation 1 */
/*
data _null_;
  set income1;
  file '\\cdc.gov\private\A728\wdq7\imputation_diagnostics\data\income2016_1.txt';
  put faminci2 income_M;
run;
*/

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
run;

proc sort data=FAMILYXX ;
by HHX FMX ;
run;

data family_w_income1;
merge income1 FAMILYXX;
by HHX FMX ;
run;

data household;
set NHIS.HOUSEHLD;
keep HHX PSTRAT PPSU;
run;

proc sort data=household;
by HHX;
run;

data family_w_income;
merge family_w_income1 household;
by HHX ;
if Familyind=1;
run;


/*
proc contents data=family_w_income;
run;
*/

/* histogram of income */
/*
proc univariate data=family_w_income normal plot;
  var income_ori;
  histogram income_ori;
run;
*/

**********************************
Chapter 1, Table 1.9, missing data pattern
********************************;
proc mi data=family_w_income nimpute=0;
var income_ori CURWRKN FDIVDYN FM_EDUC1_Recode ;
ods select misspattern;
run;


**********************************
Chapter 1, Table 1.11, missing data pattern
********************************;
/*
proc surveymeans data=family_w_income;
var income_ori;
class incgrp5;
strata PSTRAT;
cluster PPSU;
weight wtfa_fam;
run;
*/

proc freq data=family_w_income;
table incgrp5;
weight wtfa_fam;
run;

proc freq data=family_w_income;
table incgrp5*income_M;
weight wtfa_fam;
run;

proc sort data=family_w_income;
  by FM_EDUC1_Recode;
run;

/*
proc surveymeans data=family_w_income;
var income_ori;
by FM_EDUC1_Recode;
strata PSTRAT;
cluster PPSU;
weight wtfa_fam;
run;
*/

proc sort data=family_w_income;
  by FSALYN;
run;

/*
proc surveymeans data=family_w_income;
var income_ori;
by FSALYN;
strata PSTRAT;
cluster PPSU;
weight wtfa_fam;
run;
*/

**********************************
Chapter 1, Table 1.12, missing data pattern
********************************;
/*
proc means data=family_w_income;
var income_ori;
class FM_EDUC1_Recode;
run;
*/

proc freq data=family_w_income;
table FM_EDUC1_Recode/missing;
weight wtfa_fam;
run;

proc freq data=family_w_income;
table FM_EDUC1_Recode*income_M/missing;
weight wtfa_fam;
run;

**********************************
Chapter 1, Table 1.13,
********************************;


proc freq data=family_w_income;
table FM_EDUC1_Recode*incgrp5*income_M/missing;
weight wtfa_fam;
run;

/*
proc freq data=family_w_income;
table fsalyn*income_M/missing;
weight wtfa_fam;
run;
*/

**********************************
Chapter 1, Table 1.14, 
********************************;

proc surveylogistic data=family_w_income ;
class incgrp5;
model income_M (descending) =incgrp5;
where incgrp5 in ( 1 2 3 4);
strata PSTRAT;
cluster PPSU;
weight wtfa_fam;
run;


proc surveylogistic data=family_w_income ;
class incgrp5 FM_EDUC1_Recode;
model income_M (descending) =incgrp5 FM_EDUC1_Recode;
where incgrp5 in ( 1 2 3 4);
strata PSTRAT;
cluster PPSU;
weight wtfa_fam;
run;

proc surveylogistic data=family_w_income ;
class incgrp5 FM_EDUC1_Recode FSALYN;
model income_M (descending) =incgrp5 FM_EDUC1_Recode FSALYN;
where incgrp5 in ( 1 2 3 4);
strata PSTRAT;
cluster PPSU;
weight wtfa_fam;
run;

