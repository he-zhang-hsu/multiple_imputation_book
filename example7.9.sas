libname t "\\cdc.gov\private\A728\wdq7\imputation_diagnostics\data";


/*Descriptive analysis*/

/* Tables 7.7-7.9 */
proc freq data=t.sbirth;
  tables RESTATUS PLDEL3 REGNRES CITRSPOP METRORES CNTRSPOP MRACE3
         MEDUC6 DMAR MPLBIRR ADEQUACY MPRE5 FRACE4 CSEX DPLURAL DELMETH5 
         NLBNL NLBND NOTERM MEDICALRISK OBSTETRIC CONGNTL NEWBORN LABCOMP NPREVIS CIGAR DRINK
         / missing nocum nocol nopercent;
run;

proc means data=t.sbirth;
     var DBIRWT DMAGE DFAGE FMAPS WTGAIN DGESTAT NPREVIS CIGAR DRINK CONGNTL	MEDICALRISK	OBSTETRIC
				  NEWBORN		LABCOMP;
run;

/* prepare a working dataset */
data sbirth_missing;
  set t.sbirth;
  keep RESTATUS PLDEL3 REGNRES CITRSPOP METRORES CNTRSPOP MRACE3
         MEDUC6 DMAR MPLBIRR ADEQUACY MPRE5 FRACE4 CSEX DPLURAL DELMETH5 
         NLBNL NLBND NOTERM MEDICALRISK OBSTETRIC CONGNTL NEWBORN LABCOMP
		DBIRWT DMAGE DFAGE FMAPS WTGAIN DGESTAT NPREVIS CIGAR DRINK;
run;

proc contents data=sbirth_missing;
run;


/* FCS imputation */
/* this is for Example 7.9 of the book */
proc mi data=sbirth_missing seed=197789 out=sbirth_impute nimpute=1 ;
   class RESTATUS PLDEL3 REGNRES CITRSPOP METRORES CNTRSPOP MRACE3
         MEDUC6 DMAR MPLBIRR ADEQUACY MPRE5 FRACE4 CSEX DPLURAL DELMETH5;

fcs nbiter=20 logistic (pldel3 / details likelihood=augment);
fcs nbiter=20 logistic (meduc6 / details link=glogit likelihood=augment);
fcs nbiter=20 logistic (mplbirr / details  likelihood=augment);
fcs nbiter=20 logistic (adequacy / details link=glogit likelihood=augment);
fcs nbiter=20 logistic (mpre5 / details link=glogit likelihood=augment);
fcs nbiter=20 logistic (frace4 / details link=glogit likelihood=augment);
fcs nbiter=20 logistic (delmeth5 / details link=glogit likelihood=augment);
fcs nbiter=20 regpmm (nlbnl / details);
fcs nbiter=20 regpmm (nlbnd / details );
fcs nbiter=20 regpmm (noterm / details );
fcs nbiter=20 regpmm (medicalrisk / details);
fcs nbiter=20 regpmm (obstetric / details);
fcs nbiter=20 regpmm (congntl / details);
fcs nbiter=20 regpmm (newborn / details);
fcs nbiter=20 regpmm (labcomp / details);
fcs nbiter=20 regpmm (dbirwt / details);
fcs nbiter=20 regpmm (dfage / details);
fcs nbiter=20 regpmm (fmaps / details);
fcs nbiter=20 regpmm (wtgain / details);
fcs nbiter=20 regpmm (nprevis / details);
fcs nbiter=20 regpmm (cigar / details);
fcs nbiter=20 regpmm (drink / details);
fcs nbiter=20 regpmm (dgestat / details);


   var RESTATUS PLDEL3 REGNRES CITRSPOP METRORES CNTRSPOP MRACE3
         MEDUC6 DMAR MPLBIRR ADEQUACY MPRE5 FRACE4 CSEX DPLURAL DELMETH5 
         NLBNL NLBND NOTERM MEDICALRISK OBSTETRIC CONGNTL NEWBORN LABCOMP
		DBIRWT DMAGE DFAGE FMAPS WTGAIN DGESTAT NPREVIS CIGAR DRINK;

run;

proc contents data=sbirth_impute;
run;

/* There are a few numbers which take decimals really close to 0 so we round them */
data sbirth_impute_new;
  set sbirth_impute;
  NLBNL=round(NLBNL);
  NLBND=round(NLBND);
  NOTERM=round(NOTERM);
  MEDICALRISK=round(MEDICALRISK);
  OBSTETRIC=round(OBSTETRIC);
  CONGNTL=round(CONGNTL);
  NEWBORN=round(NEWBORN);
  LABCOMP=round(LABCOMP);
  DBIRWT=round(DBIRWT);
  DMAGE=round(DMAGE);
  DFAGE=round(DFAGE);
  FMAPS=round(FMAPS);
  WTGAIN=round(WTGAIN);
  DGESTAT=round(DGESTAT);
  NPREVIS=round(NPREVIS);
  CIGAR=round(CIGAR);
  DRINK=round(DRINK);
run;

/* Frequency of categorical variables for completed data */
proc freq data=sbirth_impute_new;
  tables RESTATUS PLDEL3 REGNRES CITRSPOP METRORES CNTRSPOP MRACE3
         MEDUC6 DMAR MPLBIRR ADEQUACY MPRE5 FRACE4 CSEX DPLURAL DELMETH5 
         NLBNL NLBND NOTERM MEDICALRISK OBSTETRIC CONGNTL NEWBORN LABCOMP
         / missing nocum nocol nopercent;
run;

/* Mean and SD of continuous variables for completed data */
proc means data=sbirth_impute_new;
     var DBIRWT DMAGE DFAGE FMAPS WTGAIN DGESTAT NPREVIS CIGAR DRINK;
run;

/* the end of multiple imputation book example */
