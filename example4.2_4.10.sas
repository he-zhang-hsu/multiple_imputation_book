libname t "\\cdc.gov\private\A728\wdq7\imputation_diagnostics\data";

/* read the full data */
data sbirth_output_missing_withmatch;
  set t.sbirth_output_missing_withmatch;
run;


/* save the dataset */

data _null_;
  set sbirth_output_missing_withmatch;
  file '\\cdc.gov\private\A728\wdq7\imputation_diagnostics\data\dgestat.txt';
  put dgestat dbirwt dmage  dfage  fmaps  wtgain  nprevis  cigar  drink;
run;


/* some descriptive statistics */
/*
proc freq data=sbirth_output_missing_withmatch;
  tables dgestat m dbirwt dmage dfage fmaps wtgain nprevis cigar drink /missing;
run;
*/


/* Example 4.2, the predictor includes dbirwt */
/* multiple imputation using normal linear regression model with 20 imputations */
proc mi data=sbirth_output_missing_withmatch out=sbirth_mi_model2 seed=197789 nimpute=20;
 var  dbirwt dmage  dfage fmaps wtgain  nprevis cigar  drink dgestat;
 monotone reg(dgestat= dbirwt dmage  dfage fmaps wtgain  nprevis cigar  drink/ details);
 
run;

/*
data test;
  set sbirth_output_missing_withmatch;
  dbirwt2=dbirwt**2;
run;

proc mi data=test out=sbirth_mi_model2 seed=197789 nimpute=50;
 var  dbirwt dbirwt2 dmage  dfage fmaps wtgain  nprevis cigar  drink dgestat;
 monotone reg(dgestat= dbirwt dbirwt2 dmage  dfage fmaps wtgain  nprevis cigar  drink/ details);
 
run;
*/

/* Results for Table 4.4 */
proc means data=sbirth_output_missing_withmatch mean stderr nmiss;
  var dgestat;
title 'complete-case analysis of means of dgestat';
run;



/* combining results of regressing dgestat on other variables */
proc reg data=sbirth_mi_model2 outest=outreg covout noprint;
   model dgestat= dbirwt dmage  dfage fmaps wtgain  nprevis cigar  drink;
   by _Imputation_;
run;

proc mianalyze data=outreg;
   modeleffects Intercept dbirwt dmage  dfage fmaps wtgain  nprevis cigar  drink ;
title 'regression analysis using dgestat as the outcome, MI';
run;

/* complete-case analysis results of regressing dgestat on other variables */

proc reg data=sbirth_output_missing_withmatch;
  model dgestat=dbirwt dmage  dfage fmaps wtgain  nprevis cigar  drink;
title 'regression analysis using dgestat as the outcome, CC';
run;

/* combining results of regressing dbirwt on dgestat and other variables */
proc reg data=sbirth_mi_model2 outest=outreg_rev covout noprint;
   model dbirwt= dgestat dmage  dfage fmaps wtgain  nprevis cigar  drink;
   by _Imputation_;
run;

proc mianalyze data=outreg_rev;
   modeleffects Intercept dgestat dmage  dfage fmaps wtgain  nprevis cigar  drink ;
title 'regression analysis using dgestat as the predictor, MI';
run;

/* complete-case analysis of regressing dbirwt on dgestat and other variables */
proc reg data=sbirth_output_missing_withmatch;
  model dbirwt=dgestat dmage  dfage fmaps wtgain  nprevis cigar  drink;
title 'regression analysis using dgestat as the predictor, CC';
run;

/* Exampe 4.10, the predictor of imputation excludes dbirwt */
/* Table 4.15 */
proc mi data=sbirth_output_missing_withmatch out=sbirth_mi_model2 seed=197789 nimpute=20;
 var dbirwt dmage  dfage fmaps wtgain  nprevis cigar  drink dgestat;
 monotone reg(dgestat= dmage  dfage fmaps wtgain  nprevis cigar  drink/ details); 
run;


/* combining results of regressing dbirwt on dgestat and other variables */
proc reg data=sbirth_mi_model2 outest=outreg_rev covout noprint;
   model dbirwt= dgestat dmage  dfage fmaps wtgain  nprevis cigar  drink;
   by _Imputation_;
run;

proc mianalyze data=outreg_rev;
   modeleffects Intercept dgestat dmage  dfage fmaps wtgain  nprevis cigar  drink ;
title 'regression analysis using dgestat as the predictor, MI';
run;
