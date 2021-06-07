libname t "\\cdc.gov\private\A728\wdq7\imputation_diagnostics\data";

/************************************************ 
	PSmatching.sas   adapted from
	
	Paper 185-2007  SAS Global Forum 2007
	Local and Global Optimal Propensity Score Matching
	Marcelo Coca-Perraillon
	Health Care Policy Department, Harvard Medical School, Boston, MA
	
	-------------------------------
 Treatment and Control observations must be in separate datasets such that 
    Control data includes: idC =  subject_id, pscoreC = propensity score
    Treatment data includes: idT, pscoreT
    id must be numeric	
    
	method = NN (nearest neighbor), caliper, or radius 
	
	caliper value = max for matching
	
	replacement = yes/no  whether controls can be matched to more than one case
	
	out = output data set name
	
	example call:
	
	 %PSMatching(datatreatment= T, datacontrol= C, method= NN,
	  numberofcontrols= 1, caliper=, replacement= no, out= matches);
	  
 Output format:
           Id                  Matched
       Selected     PScore       To       PScore
Obs     Control    Control    TreatID     Treat
  1      18628     0.39192     16143     0.39192
  2      18505     0.23029     16158     0.23002
  3      15589     0.29260     16112     0.29260

All other variables discarded.   Reformat for merge on subject_id with original data:

  data pairs;
    set  matches;
	subject_id = IdSelectedControl; pscore = PScoreControl; pair = _N_;
	output;
	subject_id = MatchedToTreatID; pscore = PScoreTreat; pair = _N_;
	output;
	keep subject_id pscore pair;


************************************************/

%macro PSMatching(datatreatment=, datacontrol=, method=, numberofcontrols=, caliper=,
 replacement=, out=);

/* Create copies of the treated units if N > 1 */;
 data _Treatment0(drop= i);
  set &datatreatment;
  do i= 1 to &numberofcontrols;
  RandomNumber= ranuni(12345);
output;
end;
run;
/* Randomly sort both datasets */
proc sort data= _Treatment0 out= _Treatment(drop= RandomNumber);
by RandomNumber;
run;
data _Control0;
set &datacontrol;
RandomNumber= ranuni(45678);
run;
proc sort data= _Control0 out= _Control(drop= RandomNumber);
by RandomNumber;
run;

 data Matched(keep = IdSelectedControl PScoreControl MatchedToTreatID PScoreTreat);
  length pscoreC 8;
  length idC 8;
/* Load Control dataset into the hash object */
  if _N_= 1 then do;
declare hash h(dataset: "_Control", ordered: 'no');
declare hiter iter('h');
h.defineKey('idC');
h.defineData('pscoreC', 'idC');
h.defineDone();
call missing(idC, pscoreC);
end;
/* Open the treatment */
set _Treatment;
%if %upcase(&method) ~= RADIUS %then %do;
retain BestDistance 99;
%end;
/* Iterate over the hash */
rc= iter.first();
if (rc=0) then BestDistance= 99;
do while (rc = 0);
/* Caliper */
%if %upcase(&method) = CALIPER %then %do;
if (pscoreT - &caliper) <= pscoreC <= (pscoreT + &caliper) then do;
ScoreDistance = abs(pscoreT - pscoreC);
if ScoreDistance < BestDistance then do;
BestDistance = ScoreDistance;
IdSelectedControl = idC;
PScoreControl =  pscoreC;
MatchedToTreatID = idT;
PScoreTreat = pscoreT;
end;
end;
%end;
/* NN */
%if %upcase(&method) = NN %then %do;
ScoreDistance = abs(pscoreT - pscoreC);
if ScoreDistance < BestDistance then do;
BestDistance = ScoreDistance;
IdSelectedControl = idC;
PScoreControl =  pscoreC;
MatchedToTreatID = idT;
PScoreTreat = pscoreT;
end;
%end;

%if %upcase(&method) = NN or %upcase(&method) = CALIPER %then %do;
rc = iter.next();
/* Output the best control and remove it */
if (rc ~= 0) and BestDistance ~=99 then do;
output;
%if %upcase(&replacement) = NO %then %do;
rc1 = h.remove(key: IdSelectedControl);
%end;
end;
%end;
/* Radius */
%if %upcase(&method) = RADIUS %then %do;
if (pscoreT - &caliper) <= pscoreC <= (pscoreT + &caliper) then do;
IdSelectedControl = idC;
PScoreControl =  pscoreC;
MatchedToTreatID = idT;
PScoreTreat = pscoreT;
output;
end;
rc = iter.next();
%end;
end;
run;
/* Delete temporary tables. Quote for debugging */
proc datasets;
delete _:(gennum=all);

run;
 data &out;
   set Matched;
 run;
%mend PSMatching;

/* This macro performs bivariate matching */
%macro PS_cov_group_matching(numbermi=&numbermi, var=&var, data=&data, imputevar=&imputevar);

%if %sysfunc(exist(a)) %then
        %do;
        proc datasets nolist;
          delete a b matches pair;
        run;
        quit;
    %end;

 
%do i=1 %to &numbermi; 

data a;
  set &data;
 
  where &imputevar=&i;

run;

data a_missing;
  set a;
  idT=subject_id;
  pscoreT=presponse;
  where m=1 and &var = 1;
  keep idT pscoreT;
run;

data a_observed;
  set a;
  idC=subject_id;
  pscoreC=presponse;
  where m=0 and &var =1;
  keep idC pscoreC;
run;

%PSMatching(datatreatment= a_missing, datacontrol= a_observed, method= NN,
	  numberofcontrols= 1, caliper=, replacement= no, out= matches);


  
data pairs;
    set  matches;
	subject_id = IdSelectedControl; pscore = PScoreControl; covar_pair = _N_;
	output;
	subject_id = MatchedToTreatID; pscore = PScoreTreat; covar_pair = _N_;
	output;
	keep subject_id pscore covar_pair;
run;



proc sort data=a;
  by subject_id;
run;

proc sort data=pairs;
  by subject_id;
run;

data a;
  merge a pairs;
  by subject_id;
run;

proc append base=b data=a force;

run;

%end;

data &data._withmatch;
  set b;
run;


%mend PS_cov_group_matching;


/******************************************************************************/
/* Program : stddiff.sas
/* Purpose : SAS macro to calculate the Standardized Difference
/* Usage : %stddiff(inds = Studydata, groupvar = dex,
/* numvars = age bmi/r glucose,
/* charvars = female surgtype,
/* stdfmt = 8.5,
/* outds = std_result);
/*******************************************************************************/
/* NOTE: All binary variables must be coded as 0 and 1 in the dataset
/* PARAMETERS:
/* inds:       input dataset
/* groupvar:   a binary variable, must be coded as 0 and 1
/* numvars:    a list of continuous variables.
/*             "/r" denotes to use the rank-based mean and SD to calculate Stddiff
/* charvars:   a list of categorical variables. If a variable is a binary categorical variable,
/*             it must be coded as 0 and 1 since we use the level = 0 as the reference level.
/* wtvar:      a weight variable.
/* stdfmt = 8.5 the format of Standardized Difference
/* outds output result dataset
/*********************************************************************************/

/*options  symbolgen mlogic mprint; */
 
%macro 	stddiff( inds = , 
				groupvar =, 
				numvars = , 
				charvars = ,  
				wtvar = ,
   				stdfmt = 8.4,
				outds = stddiff_result ); 

/* create a table to store stddiff */
proc sql; 
   create table &outds.  
       (VarName char(32), 
   		Stddiff char (10)
       ); 
quit; 

/* delete records if the group variable is missing */

data base_data; 
 	set &inds.; 
 	where &GroupVar. ne .; 
run; 

/* remove leading or tailing blanks */
%let groupvar = %sysfunc(strip(&GroupVar.)); 

					/****************************************/
					/* part 1: compare continuous variables */
					/****************************************/

%if %length(&numvars.) > 0 %then %do; 

/* remove multiple blanks and get the total number of continuous variables */
	%let numvar = %sysfunc(compbl(&numvars.)); 
	%let numvar = %sysfunc(strip(&numvar.)); 
	%let n_convar = %sysfunc(countc(&numvar.,' ')); 
	%let n_convar = %eval(&n_convar. + 1); 

/* summarize variables one-by-one */
	%do ii = 1 %to &n_convar.; 
    	%let convar = %sysfunc(scan(&numvar.,&ii.,' ')); 

    /* if requires rank-based mean and std for skewed variables */
		%if %index(&convar., /r) > 0 %then %do; 
    		%let convar = %sysfunc(scan(&convar.,1,'/')); 
    		%let convar = %sysfunc(strip(&convar.)); 

    		data temp_1; 
     			set base_data (keep = &groupvar. &convar. &wtvar.); 
    		run; 

    /* rank a variable */
    		proc rank data=temp_1 out=temp_2; 
          		var &convar.; 
        		ranks rank_&convar.; 
    		run; 

    /* get ranked-mean and sd */

			proc means data = temp_2;
				class &groupvar.;
				var rank_&convar.;
				weight &wtvar.;
				output out = temp_3 mean = _mean_  std = _std_;
			run;

			data  temp_3;
				set temp_3;
				where _type_ = 1;
			run;

			proc sort data = temp_3;
				by &groupvar.;
			run;
	     	%end; 
   
	/* for normal-distributed variable */

		%else %do; 
	    	%let convar = %sysfunc(strip(&convar.)); 
	    	data temp_1; 
	    	 	set base_data (keep = &groupvar. &convar. &wtvar.); 
	    	run; 
	    	data temp_2; 
	     		set temp_1; 
	    	run; 

	    /* get mean and sd */

			proc means data = temp_2;
				class &groupvar.;
				var &convar.;
				weight &wtvar.;
				output out = temp_3 mean = _mean_  std = _std_;
			run;

			data  temp_3;
				set temp_3;
				where _type_ = 1;
			run;

			proc sort data = temp_3;
				by &groupvar.;
			run;

    		%end; 

/* calculate stddiff */   
	   proc sql; 
	    	create table temp_4 as  
	    		select (a._mean_ - b._mean_)/ 
				sqrt((a._std_**2 + b._std_**2)/2) as d 
	    		from temp_3(where = (&groupvar = 1)) as a, 
	      			 temp_3(where = (&groupvar = 0)) as b; 
	   quit; 
	     
	   data temp_5; 
	   		set temp_4; 
	        stddiff = compress(put(d,&stdfmt.)); 
	        keep stddiff; 
	    run; 

	/* insert into std table */
	   proc sql noprint; 
	    	select stddiff into: std_value from temp_5; 
	    	insert into &outds.  values("&convar.", "&std_value."); 
	   quit; 

	/* delete temporary data sets */

	   proc datasets lib = work nodetails nolist; 
	    delete  temp_1 - temp_5; 
	   quit; 
   	   %end;  
%end; 

		/**********************************************/
		/* part 2: compare categorical variables      */
		/**********************************************/

%if %length(&charvars.) > 0 %then %do; 
	%let n_charvar = %sysfunc(countw(&charvars.)); 

/* get column percents for each levels of the variable by the group */
	%do jj = 1 %to &n_charvar.; 
   		%let char_var = %scan(&charvars., &jj.); 
   		%let char_var = %sysfunc(strip(&char_var.)); 
		  data temp_1; 
		   	set base_data (keep = &groupvar. &char_var. &wtvar.); 
		  run; 
    
		  proc sql; 
		   	create table temp_2 as 
		   	select distinct &char_var. as &char_var.
		   	from temp_1
			where &char_var. is not missing; 
		  quit; 

		  proc sql noprint; 
		   	select count(*) into :_mylevel_ from temp_2; 
		  quit; 

   		%let _mylevel_ = %sysfunc(strip(&_mylevel_.)); 

		  data temp_3; 
		   	set temp_2; 
		   		do &groupvar. = 0,1 ; 
		    	output; 
		   		end; 
		  run;

		ods output CrossTabFreqs = temp_4; 
  		proc freq data = temp_1; 
   			table &char_var. * &groupvar.; 
			%if %length(&wtvar.) > 0 %then %do;
				weight &wtvar.;
				%end;
  		run; 

	  	proc sql; 
	   		create table  temp_5 as 
	   		select a.*, b.ColPercent 
	   		from temp_3 as a 
	   		left join temp_4 as b 
	   		on 	a.&groupvar. = b.&groupvar. and  
	     		a.&char_var. = b.&char_var.; 
	  	quit; 

	  	data temp_6; 
	   		set temp_5; 
	   		if ColPercent = . then ColPercent = 0; 
	  	run; 

  		proc sort data = temp_6 out = catfreq; 
   			by &groupvar. &char_var.; 
  		run; 
  
  		proc datasets lib = work nodetails nolist; 
   			delete  temp_1 - temp_6; 
  		quit; 

/* if a categorical variable only has one level: 0 or 1 */
/* stddiff = 0 */
		%if &_mylevel_. = 1 %then %do; 
  			proc sql noprint; 
   				insert into &outds.  values("&char_var.", "0"); 
  			quit; 
   		%end; 

/* if a categorical variable  has two level: 0 and 1 */
/* it is a binary variable, using two sample proportation formula */
		%else %if &_mylevel_. = 2 %then %do; 

  			data temp_7; 
   				set catfreq; 
  				where &char_var. = 1; 
   				ColPercent = ColPercent/100; 
  			run; 

  			proc sql; 
   				create table temp_8 as  
   				select (a.ColPercent - b.ColPercent)/(sqrt((a.ColPercent*(1- 
						a.ColPercent) +  
     					b.ColPercent*(1-b.ColPercent))/2)) as d 
   				from temp_7(where = (&groupvar = 1)) as a, 
     		    	 temp_7(where = (&groupvar = 0)) as b; 
 		 	quit; 
     
  			data temp_9; 
          		set temp_8; 
          		stddiff = compress(put(d,&stdfmt.)); 
                keep stddiff; 
            run; 

  			proc sql noprint; 
   				select  stddiff into: std_value from temp_9; 
   					insert into &outds.  values("&char_var.", "&std_value."); 
  			quit; 

  			proc datasets lib = work nodetails nolist; 
   				delete  temp_7 temp_8 temp_9; 
  			quit; 
    	%end; 
/* if a categorical variable  has more than two level such as a, b and c */
		%else %if &_mylevel_. > 2 %then %do; 
   			%let _k_ = %eval(&_mylevel_. - 1); 
   			%let _k_ = %sysfunc(strip(&_k_.)); 
  			data temp_7; 
   				set catfreq; 
  				by &groupvar.; 
   				if last.&groupvar. then delete; 
   				ColPercent = ColPercent/100; 
  			run; 

  			proc sql noprint; 
   				select ColPercent into :tlist separated by ' '  
				from temp_7 where &groupvar. = 1; 

   				select ColPercent into :clist separated by ' '  
				from temp_7 where &groupvar. = 0; 
  			quit; 

/* vector T, C and T-C */
  			data t_1; 
   				array t{*}  t1- t&_k_.   (&tlist.); 
   				array c{*}  c1- c&_k_.   (&clist.); 
   				array tc{*} tc1 - tc&_k_. ; 
   				do i = 1 to dim(t); 
    				tc{i} = t{i} - c{i}; 
   				end; 
   			drop i; 
  			run; 

/* each column has one element of a S covariance matrix (k x k) */

			%let _dm = ; 
			%let _dm = %eval(&_k_.*&_k_.); 
  			data covdata; 
   				array t{*}  t1- t&_k_.  (&tlist.); 
   				array c{*}  c1- c&_k_.   (&clist.); 
   				array cv{&_k_.,&_k_.} x1 -x&_dm.; 
   				do i = 1 to &_k_.; 
    				do j = 1 to &_k_.; 
     					if i = j then do; 
      						cv{i,j} = 0.5*(t{i}*(1-t{i}) + c{i}*(1-c{i})); 
      						end; 
     					else do; 
      						cv{i,j} = -0.5 * (t[i] * t[j] + c[i] * c[j]); 
      						end; 
    					if cv{&_k_.,&_k_.] ne . then output; 
    				end; 
  				end; 
  			run; 

  			proc transpose data = covdata(keep = x1 -x&_dm.) out = covdata_1; 
  			run; 

  			data covdata_2; 
   				set covdata_1; 
   				retain id gp 1; 
   				if mod(_n_ - 1,&_k_.) = 0 then gp = gp + 1; 
  			run; 

		  	proc sort data = covdata_2 ; 
		   		by gp id; 
		  	run;   

			data covdata_3; 
		   		set covdata_2; 
		   		by gp id; 
		   		retain lp; 
		   		if first.gp then lp = 0; 
		   		lp = lp+1; 
		  	run; 

/* transpose to a S variance-covariance matrix format */
           
		  	data covdata_4; 
		   		set covdata_3; 
		   		retain y1-y&_k_.; 
		   		array cy{1:&_k_.} y1-y&_k_.; 
		   		by gp id; 
		   		if first.gp then do; 
		    		do k = 1 to &_k_.; 
		     			cy{k} = .; 
				    end; 
		   		end; 
		   		cy{lp} = col1; 
		   		if last.gp then output; 
		   		keep y:; 
		  	run; 

/* get inverse of S matrix */
		  data A_1; 
		   set covdata_4; 
		   array _I{*} I1-I&_k_.; 
		   do j=1 to &_k_.; 
		    if j=_n_ then _I[j]=1;  
		    else _I[j]=0; 
		   end; 
		   drop j; 
		  run; 

/* solve the inverse of the matrix */

  %macro inv; 
    	%do j=1 %to &_k_.; 
    		proc orthoreg data=A_1 outest=A_inv_&j.(keep=y1-y&_k_.) 
     			noprint singular=1E-16; 
     			model I&j=y1-y&_k_. /noint; 
    		run; 
    		quit; 
    	%end; 

   		data A_inverse; 
    		set %do j=1 %to &_k_.; 
     		A_inv_&j 
     	%end;; 
   		run; 
  %mend; 
  %inv; 

  		proc transpose data=A_inverse out=A_inverse_t; 
  		run; 

   /* calculate the mahalanobis distance */
  		data t_2; 
   			set A_inverse_t; 
   			array t{*}  t1- t&_k_.  (&tlist.); 
   			array c{*}  c1- c&_k_.  (&clist.); 
   			i = _n_; 
   			trt = t{i}; 
   			ctl = c{i}; 
   			tc = t{i} - c{i}; 
  		run; 
 
		data t_3; 
   			set t_2; 
   			array aa{&_k_.} col1 - col&_k_.; 
   			array bb{&_k_.} bb1- bb&_k_.; 
   			do i = 1 to &_k_.; 
    			bb{i} = aa{i}*tc; 
   			end; 
  		run; 

  		proc summary data = t_3 ; 
   			var bb1-bb&_k_.; 
   			output out = t_4 sum =; 
  		run; 

  		data t_5; 
   			merge t_1 t_4; 
   			array d1{*} tc1- tc&_k_. ; 
   			array d2{*} bb1-bb&_k_.; 
   			array d3{*} y1-y&_k_.; 
   			do i = 1 to &_k_.; 
   				d3{i} = d1{i}*d2{i}; 
   			end; 
   			d = sqrt(sum(of y1-y&_k_.)); 
   			stddiff = compress(put(d,&stdfmt.));      
   			keep stddiff; 
  		run; 

  		proc sql noprint; 
   			select  stddiff into: std_value from t_5; 
   			insert into &outds.  values("&char_var.", "&std_value."); 
  		quit; 
   
  		proc datasets lib = work nodetails nolist; 
   			delete  covdata covdata_1 covdata_2 covdata_3 covdata_4 
      		A_1 A_inverse A_inverse_t t_1 t_2 t_3 t_4 t_5
     		A_inv_:; 
  		quit; 
	   %end; 
	%end; 
%end; 

proc datasets lib = work nodetails nolist; 
  delete Catfreq  Base_data temp_7; 
quit; 

proc print data = &outds.; 
	title 'Calculated Standardized Difference';
run; 

title;

%mend; 

/* This macro evaluate the standard differences for multiply imputed data */

%macro mi_stddiff (numbermi=&numbermi, var=&var, data=&data, groupvar=&groupvar, pairvar=&pairvar, imputevar=&imputevar, cutoff=&cutoff);

%if %sysfunc(exist(stddiff_result)) %then
        %do;
        proc datasets nolist;
          delete stddiff_result stddiff_result_final a;
        run;
        quit;
    %end;

 
%do i=1 %to &numbermi; 

data a;
  set &data;
  where &imputevar=&i and &pairvar > 0;
run;

%stddiff( inds = a, 
               groupvar =&groupvar, 
               numvars = &var, 
               charvars = ,  
               stdfmt = 8.4,
               outds = stddiff_result );


proc append base=stddiff_result_final data=stddiff_result force;

run;
    

 
%end;

data stddiff_result_final;
  set stddiff_result_final;
  _imputation_=_n_;
  stddiff_new=stddiff*1.0;
  if abs(stddiff_new) > &cutoff then stddiff_indi=1;
  else stddiff_indi=0;
run;

proc means data=stddiff_result_final noprint;
  var stddiff_new;
  output out=b mean(stddiff_new)=stddiff_mean mean(stddiff_indi)=stddiff_indi_mean; 
run;

data stddiff_result_&data.;
  if _N_=1 then set b (keep=stddiff_mean stddiff_indi_mean);
  set stddiff_result_final;
run;

proc print data=stddiff_result_&data.;
title "standardized differences for imputed data &data.";
run;

%mend mi_stddiff;

/* This macro evaluates he variance ratios of multiply imputed data */

%macro mi_varratio (numbermi=&numbermi, var=&var, data=&data, groupvar=&groupvar, pairvar=&pairvar, imputevar=&imputevar, cutoff1=&cutoff1, cutoff2=&cutoff2);
%if %sysfunc(exist(varratio_result)) %then
        %do;
        proc datasets nolist;
          delete varratio_result varratio_result_final a b c;
        run;
        quit;
    %end;

 
%do i=1 %to &numbermi; 

data a;
  set &data;
  where &imputevar=&i and &pairvar > 0;
run;

proc means data = a noprint;
				class &groupvar;
				var &var;
				output out = b var = _var_;
			run;

data  b;
	set b;
	where _type_ = 1;
run;

proc sort data = b;
	by &groupvar;
	run;
			
 proc sql; 
	  create table varratio_result as  
	    		select 	a._var_ / b._var_ as varratio 
	    		from b (where = (&groupvar = 1)) as a, 
	      			 b (where = (&groupvar = 0)) as b; 
	   quit; 


proc append base=varratio_result_final data=varratio_result force;
run;
    

%end;

data varratio_result_final;
  set varratio_result_final;
  _imputation_=_n_;
  if varratio > &cutoff1 or varratio < &cutoff2 then varratio_indi=1;
  else varratio_indi=0;
run;

proc means data=varratio_result_final noprint;
  var varratio;
  output out=c mean(varratio)=varratio_mean mean(varratio_indi)=varratio_indi_mean; 
run;

data varratio_result_&data.;
  if _N_=1 then set c (keep=varratio_mean varratio_indi_mean);
  set varratio_result_final;
run;

proc print data=varratio_result_&data.;
title "variance ratio for imputed data &data.";
run;

%mend mi_varratio;

/* This macro subgroup the dataset into 5 groups */
%macro subset(data=&data, classvar=&classvar);

%do i=1 %to 5;
  data &data._s&i;
    set &data;
    where &classvar=&i;
  run;
%end;

%mend subset;

/* This macro performs bivariate matching */
/* create 5 groups based on the quantiles of the outcome variable */
%macro PSoutcomematching(numbermi=&numbermi, data=, groupvar=, propensityvar=, outcomevar=, imputevar=);  

%if %sysfunc(exist(&data)) %then
        %do;
        proc datasets nolist;
          delete &data._withmatch;
        run;
        quit;
    %end;


%do i=1 %to &numbermi; 

data temp;
  set &data;
  where &imputevar=&i;
  sub_id=_n_;
run;

proc sort data=temp;
  by decending &groupvar;
run;


data temp_missing;
  set temp;
  idT=sub_id;
  pscoreT=&propensityvar;
  outcomeT=&outcomevar;
  where &groupvar=1;
  group=&groupvar;
  keep idT pscoreT outcomeT group;
run;
/*
proc print data=sbirth_mi_model1_missing (obs=100);
run;
*/
data temp_observed;
  set temp;
  idC=sub_id;
  pscoreC=&propensityvar;
  group=&groupvar;
  outcomeC=&outcomevar;
  where &groupvar=0;
  keep idC pscoreC outcomeC group;
run;
/*
proc print data=sbirth_mi_model1_observed (obs=100);
run;
*/
data _TreatControl;
set temp_missing (rename=(idT= id pscoreT= pscore outcomeT= &outcomevar))
temp_observed (rename=(idC= id pscoreC= pscore outcomeC= &outcomevar));
run;
/*
proc contents data=_TreatControl;
run;
*/
proc princomp data= _TreatControl std out= _TreatControlDC;
var pscore &outcomevar;
run;

data _Treatment0(rename=(prin1= pscoreT prin2= outcomeT id= idT))
_Control0(rename=(prin1= pscoreC prin2= outcomeC id= idC));
set _TreatControlDC;
RandomNumber= ranuni(123456);
if group= 1 then output _Treatment0;
else if group= 0 then output _Control0;
run;
/*
proc print data=_Treatment0 (obs=100);
run;

proc print data=_Control0 (obs=100);
run;
*/
proc sort data= _Treatment0 out= _Treatment;
by RandomNumber;
run;

proc sort data= _Control0 out= _Control;
by RandomNumber;
run;

data MatchedMH(keep = IdSelectedControl MatchedToTreatID);
length idC pscoreC outcomeC 8;
if _N_= 1 then do;
declare hash h(dataset: "_Control", ordered: 'no');
declare hiter iter('h');
h.defineKey('idC');
h.defineData('pscoreC', 'idC', 'outcomeC');
h.defineDone();
call missing(pscoreC, idC, outcomeC);
end;
set _Treatment;
retain BestDistance 99;
rc=iter.first();
if (rc=0) then BestDistance= 99;
do while (rc = 0);
* Euclidean distance;
ScoreDistance = sqrt((pscoreT-pscoreC)**2 + (outcomeC-outcomeT)**2);
if ScoreDistance < BestDistance then do;
BestDistance = ScoreDistance;
IdSelectedControl = idC;
MatchedToTreatID = idT;
end;
rc = iter.next();
if (rc ~= 0) then do;
output;
rc1 = h.remove(key: IdSelectedControl);
end;
end;
run;

/* create 5 groups based on the quantiles of the outcome variable */
data MatchedMH;
  set MatchedMH;
  idC=IdSelectedControl;
run;

proc sort data=MatchedMH;
 by idC;
run;

proc sort data=temp_observed;
  by idC;
run;

data MatchedMH_final;
  merge MatchedMH (in=in_match) temp_observed;
  by idC;
  if in_match > 0;
  drop idC group;
run;

proc univariate data=MatchedMH_final noprint;
  var outcomeC;
   output out=Pctls pctlpts  = 20 40 50 60 80
                    pctlpre  = outcomeC;

run;

data MatchedMH_final;
  if _N_=1 then set pctls;
  set MatchedMH_final;
  if outcomeC <= outcomeC20 then outcome_match_class=1;
  else if outcomeC <= outcomeC40 then outcome_match_class=2;
  else if outcomeC <= outcomeC60 then outcome_match_class=3;
  else if outcomeC <= outcomeC80 then outcome_match_class=4;
  else outcome_match_class=5;
run;

data outcome_pairs;
    set  MatchedMH_final;
	sub_id = IdSelectedControl; outcome_pair = _N_;
	output;
	sub_id = MatchedToTreatID;  outcome_pair = _N_;
	output;
	keep sub_id  outcome_pair outcome_match_class;
run;


proc sort data=outcome_pairs;
  by sub_id;
run;

proc sort data=temp;
  by sub_id;
run;


data temp_withmatch;
  merge temp outcome_pairs;
  by sub_id;
run;

proc append base=&data._withmatch data=temp_withmatch force;

run;
    

%end;

%mend PSoutcomematching;

/* This macro performe bivariate matching */
/* create 5 groups based on the quantiles of the propensity variable */
%macro PSpropensitymatching(numbermi=&numbermi, data=, groupvar=, propensityvar=, outcomevar=, imputevar=);  

%if %sysfunc(exist(&data)) %then
        %do;
        proc datasets nolist;
          delete &data._withmatch;
        run;
        quit;
    %end;


%do i=1 %to &numbermi; 

data temp;
  set &data;
  where &imputevar=&i;
  sub_id=_n_;
run;

proc sort data=temp;
  by decending &groupvar;
run;


data temp_missing;
  set temp;
  idT=sub_id;
  pscoreT=&propensityvar;
  outcomeT=&outcomevar;
  where &groupvar=1;
  group=&groupvar;
  keep idT pscoreT outcomeT group;
run;
/*
proc print data=sbirth_mi_model1_missing (obs=100);
run;
*/
data temp_observed;
  set temp;
  idC=sub_id;
  pscoreC=&propensityvar;
  group=&groupvar;
  outcomeC=&outcomevar;
  where &groupvar=0;
  keep idC pscoreC outcomeC group;
run;
/*
proc print data=sbirth_mi_model1_observed (obs=100);
run;
*/
data _TreatControl;
set temp_missing (rename=(idT= id pscoreT= pscore outcomeT= &outcomevar))
temp_observed (rename=(idC= id pscoreC= pscore outcomeC= &outcomevar));
run;
/*
proc contents data=_TreatControl;
run;
*/
proc princomp data= _TreatControl std out= _TreatControlDC;
var pscore &outcomevar;
run;

data _Treatment0(rename=(prin1= pscoreT prin2= outcomeT id= idT))
_Control0(rename=(prin1= pscoreC prin2= outcomeC id= idC));
set _TreatControlDC;
RandomNumber= ranuni(123456);
if group= 1 then output _Treatment0;
else if group= 0 then output _Control0;
run;
/*
proc print data=_Treatment0 (obs=100);
run;

proc print data=_Control0 (obs=100);
run;
*/
proc sort data= _Treatment0 out= _Treatment;
by RandomNumber;
run;

proc sort data= _Control0 out= _Control;
by RandomNumber;
run;

data MatchedMH(keep = IdSelectedControl MatchedToTreatID);
length idC pscoreC outcomeC 8;
if _N_= 1 then do;
declare hash h(dataset: "_Control", ordered: 'no');
declare hiter iter('h');
h.defineKey('idC');
h.defineData('pscoreC', 'idC', 'outcomeC');
h.defineDone();
call missing(pscoreC, idC, outcomeC);
end;
set _Treatment;
retain BestDistance 99;
rc=iter.first();
if (rc=0) then BestDistance= 99;
do while (rc = 0);
* Euclidean distance;
ScoreDistance = sqrt((pscoreT-pscoreC)**2 + (outcomeC-outcomeT)**2);
if ScoreDistance < BestDistance then do;
BestDistance = ScoreDistance;
IdSelectedControl = idC;
MatchedToTreatID = idT;
end;
rc = iter.next();
if (rc ~= 0) then do;
output;
rc1 = h.remove(key: IdSelectedControl);
end;
end;
run;

/* create 5 groups based on the quantiles of the propensity variable */
data MatchedMH;
  set MatchedMH;
  idC=IdSelectedControl;
run;

proc sort data=MatchedMH;
 by idC;
run;

proc sort data=temp_observed;
  by idC;
run;

data MatchedMH_final;
  merge MatchedMH (in=in_match) temp_observed;
  by idC;
  if in_match > 0;
  drop idC group;
run;

proc univariate data=MatchedMH_final noprint;
  var pscoreC;
   output out=Pctls pctlpts  = 20 40 50 60 80
                    pctlpre  = pscoreC;

run;

data MatchedMH_final;
  if _N_=1 then set pctls;
  set MatchedMH_final;
  if pscoreC <= pscoreC20 then pscore_match_class=1;
  else if pscoreC <= pscoreC40 then pscore_match_class=2;
  else if pscoreC <= pscoreC60 then pscore_match_class=3;
  else if pscoreC <= pscoreC80 then pscore_match_class=4;
  else pscore_match_class=5;
run;

data propensity_pairs;
    set  MatchedMH_final;
	sub_id = IdSelectedControl; propensity_pair = _N_;
	output;
	sub_id = MatchedToTreatID;  propensity_pair = _N_;
	output;
	keep sub_id  propensity_pair pscore_match_class;
run;


proc sort data=propensity_pairs;
  by sub_id;
run;

proc sort data=temp;
  by sub_id;
run;


data temp_withmatch;
  merge temp propensity_pairs;
  by sub_id;
run;

proc append base=&data._withmatch data=temp_withmatch force;

run;
    

%end;


%mend PSpropensitymatching;

/* subgroup means of dgestat by variables */
%macro groupmeans(group=&group);

proc sort data=t.sbirth;
  by &group;
run;

proc means data=t.sbirth mean stderr n;
  var dgestat;
  by &group;
run;

%mend groupmeans;

%macro continuousgroupmeans(var=&var);

proc sort data=sbirth;
  by m;
run;

proc means data=sbirth mean stderr n;
  var &var;
  by m;
run;

%mend continuousgroupmeans;


%macro grouplogistic(var=&var);

proc logistic data=sbirth;
  class &var;
  model m=&var;
run;

%mend grouplogistic;

/* Data processing */
/*Descriptive analysis*/
proc freq data=t.sbirth;
  tables RESTATUS PLDEL3 REGNRES CITRSPOP METRORES CNTRSPOP MRACE3
         MEDUC6 DMAR MPLBIRR ADEQUACY MPRE5 FRACE4 CSEX DPLURAL DELMETH5 
         NLBNL NLBND NOTERM MEDICALRISK OBSTETRIC CONGNTL NEWBORN LABCOMP
         / missing nocum nocol nopercent;
run;


/* run this code */
/* collapse some of the categorical variables */
data sbirth;
  set t.sbirth;

  log_dgestat=log(dgestat);
  sqrt_dgestat=sqrt(dgestat);
  cubic_dgestat=dgestat**(4/3);

  dbirwt_square=dbirwt*dbirwt;

  if dbirwt >= 3345 then dbirwt_half=1;
  else dbirwt_half=0;

  if restatus=1 then restatus_new=0;
  else if restatus in (2,3) then restatus_new=1;
  else restatus_new=.;

  pldel3_new=pldel3-1;

  regnres_new=regnres;
  if regnres=0 then regnres_new=.;

  if citrspop in (0,1,2,3) then citrspop_new=0;
  else if citrspop=9 then citrspop_new=1;
  else citrspop_new=.;

  metrores_new=metrores-1;

  if cntrspop in (0,1,2,3) then cntrspop_new=0;
  else if cntrspop=9 then cntrspop_new=1;
  else cntrspop_new=.;

  mrace3_new=mrace3;

  meduc6_new=meduc6;

  dmar_new=dmar-1;

  mplbirr_new=mplbirr-1;

  adequacy_new=adequacy;

  mpre5_new=mpre5;

  frace4_new=frace4;

  csex_new=csex-1;

  if dplural=1 then dplural_new=0;
  else if dplural in (2,3,4,5) then dplural_new=1;

  if delmeth5 in (1,2) then delmeth5_new=0;
  else if delmeth5 in (3,4) then delmeth5_new=1;
  else delmeth5_new=.;

  nlbnl_new=nlbnl;
  if nlbnl >=3 then nlbnl_new=3;

  nlbnd_new=nlbnd;
  if nlbnd  >=1 then nlbnd_new=1;

  noterm_new=noterm;
  if noterm >=2 then noterm_new=2;

  medicalrisk_new=medicalrisk;
  if medicalrisk >=2 then medicalrisk_new=2;

  obstetric_new=obstetric;
  if obstetric >=4 then obstetric_new=4;

  congntl_new=congntl;
  if congntl >=1 then congntl_new=1;

  newborn_new=newborn;
  if newborn >=1 then newborn_new=1;

  labcomp_new=labcomp;
  if labcomp >=2 then labcomp_new=2;


  /* create binary dummies */
  restatus_new_dummy=restatus_new;
  pldel3_new_dummy=pldel3_new;
  if regnres_new=1 then do; regnres_new_dummy1=1; regnres_new_dummy2=0; regnres_new_dummy3=0; end;
  if regnres_new=2 then do; regnres_new_dummy1=0; regnres_new_dummy2=1; regnres_new_dummy3=0; end;
  if regnres_new=3 then do; regnres_new_dummy1=0; regnres_new_dummy2=0; regnres_new_dummy3=1; end;
  if regnres_new=4 then do; regnres_new_dummy1=-1; regnres_new_dummy2=-1; regnres_new_dummy3=-1; end;

  citrspop_new_dummy=citrspop_new;

  metrores_new_dummy=metrores_new;

  cntrspop_new_dummy=cntrspop_new;

  if mrace3_new=1 then do; mrace3_new_dummy1=1; mrace3_new_dummy2=0; end;
  if mrace3_new=2 then do; mrace3_new_dummy1=0; mrace3_new_dummy2=1; end;
  if mrace3_new=3 then do; mrace3_new_dummy1=-1; mrace3_new_dummy2=-1; end;

  if meduc6_new=1 then do; meduc6_new_dummy1=1; meduc6_new_dummy2=0; meduc6_new_dummy3=0; meduc6_new_dummy4=0; end;
  if meduc6_new=2 then do; meduc6_new_dummy1=0; meduc6_new_dummy2=1; meduc6_new_dummy3=0; meduc6_new_dummy4=0; end;
  if meduc6_new=3 then do; meduc6_new_dummy1=0; meduc6_new_dummy2=0; meduc6_new_dummy3=1; meduc6_new_dummy4=0; end;
  if meduc6_new=4 then do; meduc6_new_dummy1=0; meduc6_new_dummy2=0; meduc6_new_dummy3=0; meduc6_new_dummy4=1; end;
  if meduc6_new=5 then do; meduc6_new_dummy1=-1; meduc6_new_dummy2=-1; meduc6_new_dummy3=-1; meduc6_new_dummy4=-1; end;

  dmar_new_dummy=dmar_new;

  mplbirr_new_dummy=mplbirr_new;

  if adequacy_new=1 then do; adequacy_new_dummy1=1; adequacy_new_dummy2=0; end;
  if adequacy_new=2 then do; adequacy_new_dummy1=0; adequacy_new_dummy2=1; end;
  if adequacy_new=3 then do; adequacy_new_dummy1=-1; adequacy_new_dummy2=-1; end;

  if mpre5_new=1 then do; mpre5_new_dummy1=1; mpre5_new_dummy2=0; mpre5_new_dummy3=0; end;
  if mpre5_new=2 then do; mpre5_new_dummy1=0; mpre5_new_dummy2=1; mpre5_new_dummy3=0; end;
  if mpre5_new=3 then do; mpre5_new_dummy1=0; mpre5_new_dummy2=0; mpre5_new_dummy3=1; end;
  if mpre5_new=4 then do; mpre5_new_dummy1=-1; mpre5_new_dummy2=-1; mpre5_new_dummy3=-1; end;

  if frace4_new=1 then do; frace4_new_dummy1=1; frace4_new_dummy2=0; end;
  if frace4_new=2 then do; frace4_new_dummy1=0; frace4_new_dummy2=1; end;
  if frace4_new=3 then do; frace4_new_dummy1=-1; frace4_new_dummy2=-1; end;

  csex_new_dummy=csex_new;

  dplural_new_dummy=dplural_new;

  delmeth5_new_dummy=delmeth5_new;

  if nlbnl_new=0 then do; nlbnl_new_dummy1=1; nlbnl_new_dummy2=0; nlbnl_new_dummy3=0; end;
  if nlbnl_new=1 then do; nlbnl_new_dummy1=0; nlbnl_new_dummy2=1; nlbnl_new_dummy3=0; end;
  if nlbnl_new=2 then do; nlbnl_new_dummy1=0; nlbnl_new_dummy2=0; nlbnl_new_dummy3=1; end;
  if nlbnl_new=3 then do; nlbnl_new_dummy1=-1; nlbnl_new_dummy2=-1; nlbnl_new_dummy3=-1; end;

  nlbnd_new_dummy=nlbnd_new;

  if noterm_new=0 then do; noterm_new_dummy1=1; noterm_new_dummy2=0; end;
  if noterm_new=1 then do; noterm_new_dummy1=0; noterm_new_dummy2=1; end;
  if noterm_new=2 then do; noterm_new_dummy1=-1; noterm_new_dummy2=-1; end;

  if medicalrisk_new=0 then do; medicalrisk_new_dummy1=1; medicalrisk_new_dummy2=0; end;
  if medicalrisk_new=1 then do; medicalrisk_new_dummy1=0; medicalrisk_new_dummy2=1; end;
  if medicalrisk_new=2 then do; medicalrisk_new_dummy1=-1; medicalrisk_new_dummy2=-1; end;

  if obstetric_new=0 then do; obstetric_new_dummy1=1; obstetric_new_dummy2=0; obstetric_new_dummy3=0; obstetric_new_dummy4=0; end;
  if obstetric_new=1 then do; obstetric_new_dummy1=0; obstetric_new_dummy2=1; obstetric_new_dummy3=0; obstetric_new_dummy4=0; end;
  if obstetric_new=2 then do; obstetric_new_dummy1=0; obstetric_new_dummy2=0; obstetric_new_dummy3=1; obstetric_new_dummy4=0; end;
  if obstetric_new=3 then do; obstetric_new_dummy1=0; obstetric_new_dummy2=0; obstetric_new_dummy3=0; obstetric_new_dummy4=1; end;
  if obstetric_new=4 then do; obstetric_new_dummy1=-1; obstetric_new_dummy2=-1; obstetric_new_dummy3=-1; obstetric_new_dummy4=-1; end;

  congntl_new_dummy=congntl_new;

  newborn_new_dummy=newborn_new;

  if labcomp_new=0 then do; labcomp_new_dummy1=1; labcomp_new_dummy2=0; end;
  if labcomp_new=1 then do; labcomp_new_dummy1=0; labcomp_new_dummy2=1; end;
  if labcomp_new=2 then do; labcomp_new_dummy1=-1; labcomp_new_dummy2=-1; end;

run;


/* run this code */
/* fit a logistic model */
/* and obtain the propensity score */
proc logistic data=sbirth;
  class regnres_new mrace3_new meduc6_new adequacy_new mpre5_new frace4_new nlbnl_new noterm_new medicalrisk_new
        obstetric_new labcomp_new;
  model m= restatus_new pldel3_new regnres_new citrspop_new metrores_new cntrspop_new mrace3_new meduc6_new
  dmar_new mplbirr_new adequacy_new mpre5_new frace4_new csex_new dplural_new delmeth5_new nlbnl_new nlbnd_new
  noterm_new medicalrisk_new obstetric_new congntl_new newborn_new labcomp_new

dbirwt dbirwt_square dmage dfage fmaps wtgain  nprevis cigar  drink 

citrspop_new*regnres_new meduc6_new * mpre5_new dmar_new * mpre5_new 

dmage*cigar / rsquare
;

output out=sbirth_output predicted=presponse;



 roc 'full' restatus_new pldel3_new regnres_new citrspop_new metrores_new cntrspop_new mrace3_new meduc6_new
  dmar_new mplbirr_new adequacy_new mpre5_new frace4_new csex_new dplural_new delmeth5_new nlbnl_new nlbnd_new
  noterm_new medicalrisk_new obstetric_new congntl_new newborn_new labcomp_new

dbirwt dmage  dfage fmaps wtgain  nprevis cigar  drink ;

 roc 'interaction' restatus_new pldel3_new regnres_new citrspop_new metrores_new cntrspop_new mrace3_new meduc6_new
  dmar_new mplbirr_new adequacy_new mpre5_new frace4_new csex_new dplural_new delmeth5_new nlbnl_new nlbnd_new
  noterm_new medicalrisk_new obstetric_new congntl_new newborn_new labcomp_new

dbirwt dmage  dfage fmaps wtgain  nprevis cigar  drink 

citrspop_new*regnres_new meduc6_new * mpre5_new dmar_new * mpre5_new 

dmage*cigar 
;

run;

/* run this code */
/* This is the estimated propensity */
data sbirth_output;
  set sbirth_output;
  propensity=log(presponse/(1-presponse));
run;

/* run this code */
proc univariate data=sbirth_output noprint;
  var presponse;
   output out=Pctls pctlpts  = 20 40 50 60 80
                    pctlpre  = presponse;

run;

proc sort data=sbirth_output;
  by m;
run;


/* run this code */
/* take the subset data for imputation */
/* where propensity score is not missing */
data sbirth_output_missing;
  set sbirth_output;
  where presponse ~=.;
  if presponse <= 0.79073 then presponse_class=1;
  else if presponse <= 0.83150 then presponse_class=2;
  else if presponse <= 0.85962 then presponse_class=3;
  else if presponse <= 0.88542 then presponse_class=4;
  else presponse_class=5;
  if presponse <= 0.84634 then presponse_half=1;
  else presponse_half=2;
  subject_id=_n_;
  if presponse_class=1 then do;
  presponse_dummy1=1; presponse_dummy2=0; presponse_dummy3=0; presponse_dummy4=0;
  end;
  if presponse_class=2 then do;
  presponse_dummy1=0; presponse_dummy2=1; presponse_dummy3=0; presponse_dummy4=0;
  end;
  if presponse_class=3 then do;
  presponse_dummy1=0; presponse_dummy2=0; presponse_dummy3=1; presponse_dummy4=0;
  end;
  if presponse_class=4 then do;
  presponse_dummy1=0; presponse_dummy2=0; presponse_dummy3=0; presponse_dummy4=1;
  end;
  if presponse_class=5 then do;
  presponse_dummy1=-1; presponse_dummy2=-1; presponse_dummy3=-1; presponse_dummy4=-1;
  end;
  presponse_m_dummy1=presponse_dummy1*m;
  presponse_m_dummy2=presponse_dummy2*m;
  presponse_m_dummy3=presponse_dummy3*m;
  presponse_m_dummy4=presponse_dummy4*m;

  if presponse_half=1 then presponse_half_dummy=1;
  else presponse_half_dummy=0;
  presponse_half_m_dummy=presponse_half_dummy*m;
run;


 /*  Control data includes: idC =  subject_id, pscoreC = propensity score
    Treatment data includes: idT, pscoreT
 */
/* run this code */
data sbirth_output_missing_missing;
  set sbirth_output_missing;
  idT=subject_id;
  pscoreT=presponse;
  where m=1;
  keep idT pscoreT;
run;

data sbirth_output_missing_observed;
  set sbirth_output_missing;
  idC=subject_id;
  pscoreC=presponse;
  where m=0;
  keep idC pscoreC;
run;

/* perform propensity mathcing before multiple imputation */
%PSMatching(datatreatment= sbirth_output_missing_missing, datacontrol= sbirth_output_missing_observed, method= NN,
	  numberofcontrols= 1, caliper=, replacement= no, out= matches);


proc print data=matches;
run;

proc univariate data=matches noprint;
  var pscorecontrol;
   output out=Pctls pctlpts  = 20 40 50 60 80
                    pctlpre  = pscorecontrol;

run;

proc print data=pctls;
run;

data matches;
  set matches;
  if pscorecontrol <= 0.76382 then p_match_class=1;
  else if pscorecontrol <= 0.80891 then p_match_class=2;
  else if pscorecontrol <= 0.84055 then p_match_class=3;
  else if pscorecontrol <= 0.87017 then p_match_class=4;
  else p_match_class=5;
run;

  
data pairs;
    set  matches;
	subject_id = IdSelectedControl; pscore = PScoreControl; pair = _N_;
	output;
	subject_id = MatchedToTreatID; pscore = PScoreTreat; pair = _N_;
	output;
	keep subject_id pscore pair p_match_class;
run;


proc sort data=pairs;
  by subject_id;
run;

proc sort data=sbirth_output_missing;
  by subject_id;
run;

data sbirth_output_missing_withmatch;
  merge sbirth_output_missing pairs;
  by subject_id;
run;

proc sort data=sbirth_output_missing_withmatch;
  by pair;
run;

proc print data=sbirth_output_missing_withmatch (obs=100);
  where pair > 0;
  var m pair subject_id presponse pscore p_match_class dbirwt_half;
run;

data final;
  set sbirth_output_missing_withmatch;
  keep m pair subject_id presponse pscore p_match_class 
  dbirwt dbirwt_square dmage  dfage fmaps wtgain  nprevis cigar  drink dgestat;
run;

data t.sbirth_output_missing_withmatch;
  * set sbirth_output_missing_withmatch;
  set final;
run;

proc contents data=t.sbirth_output_missing_withmatch;
run;

/* sbirth_output_missing_withmatch  is the working dataset for imputation diagnostics */

/* Diagnostic results for 4 imputation models  */
/* Table 12.2 */
/* model 1, including all continous covariates except dbirwt */

data sbirth_output_missing_withmatch;
  set t.sbirth_output_missing_withmatch;
run;


/* run this code */
proc mi data=sbirth_output_missing_withmatch out=sbirth_mi_model1  seed=197789 nimpute=100;
 var   dbirwt  dmage  dfage fmaps wtgain  nprevis cigar  drink dgestat;
 monotone reg(dgestat= /* dbirwt */ dmage  dfage  fmaps  wtgain  nprevis  cigar  drink/ details);
 
run;


/* evaluate the standardized differences */


/* run this code */
/* evaluate the standardized difference of the outcome */

/* grouping the dataset into 5 subgroups by the propensity score */ 
%subset(data=sbirth_mi_model1, classvar=p_match_class)

/* overall differences */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model1, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)
/*
proc print data=stddiff_result_final;
run;
*/
proc means data=stddiff_result_final;
run;

/* by propensity groups */
/* group 1 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model1_s1, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 2 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model1_s2, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)


proc means data=stddiff_result_final;
run;

/* group 3 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model1_s3, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 4 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model1_s4, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 5 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model1_s5, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)


proc means data=stddiff_result_final;
run;


/* Variance ratio */
/* overall */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model1, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.07, cutoff2=0.94)

proc means data=varratio_result_final;
run;

/* group 1 */

%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model1_s1, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)
proc print data=VARRATIO_RESULT_FINAL;
run;

proc means data=varratio_result_final;
run;

/* group 2 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model1_s2, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 3 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model1_s3, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 4 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model1_s4, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 5 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model1_s5, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;



/* strata based on the outcome */
%PSoutcomematching(numbermi=100, data=sbirth_mi_model1, groupvar=m, propensityvar=presponse, outcomevar=dgestat, imputevar=_imputation_);  

%subset(data=sbirth_mi_model1_withmatch, classvar=outcome_match_class)

/* overall differences of dbirwt */

%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 1 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s1, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 2 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s2, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 3 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s3, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)


proc means data=stddiff_result_final;
run;

/* group 4 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s4, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 5 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s5, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)


proc means data=stddiff_result_final;
run;


/* Variance ratio */

/* Overall */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.07, cutoff2=0.94)

proc means data=varratio_result_final;
run;

/* group 1 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s1, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.13, cutoff2=0.89)

proc means data=varratio_result_final;
run;

/* group 2 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s2, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.16, cutoff2=0.86)

proc means data=varratio_result_final;
run;

/* group 3 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s3, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.16, cutoff2=0.86)

proc means data=varratio_result_final;
run;

/* group 4 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s4, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.17, cutoff2=0.85)

proc means data=varratio_result_final;
run;

/* group 5 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s5, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.16, cutoff2=0.86)

proc means data=varratio_result_final;
run;

/* strata based on the propensity score */
%PSpropensitymatching(numbermi=100, data=sbirth_mi_model1, groupvar=m, propensityvar=presponse, outcomevar=dgestat, imputevar=_imputation_);  

%subset(data=sbirth_mi_model1_withmatch, classvar=pscore_match_class)

/* overall 
%mi_stddiff(numbermi=100, var= dbirwt, data=sbirth_mi_model1_withmatch, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;
*/

/* group 1 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s1, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)


proc means data=stddiff_result_final;
run;

/* group 2 */

%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s2, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 3 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s3, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 4 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s4, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 5 */

%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s5, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;


/* Variance ratio */
/* Overall 

%mi_varratio(numbermi=100, var= dbirwt, data=sbirth_mi_model1_withmatch, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.07, cutoff2=0.94)

proc means data=varratio_result_final;
run;

*/

/* group 1 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s1, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 2 */
%mi_varratio(numbermi=100, var=dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s2, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 3 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s3, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 4 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s4, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 5 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model1_withmatch_s5, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* the end of the model 1 */


/* model 2, including all continous covariates */
proc mi data=sbirth_output_missing_withmatch out=sbirth_mi_model2 seed=197789 nimpute=100;
 var  dbirwt dmage  dfage fmaps wtgain  nprevis cigar  drink dgestat;
 monotone reg(dgestat= dbirwt dmage  dfage fmaps wtgain  nprevis cigar  drink/ details);
 
run;

/* evaluate the standardized difference */
/* grouping the dataset into 5 subgroups by the propensity score */
/* for the outcome */
%subset(data=sbirth_mi_model2, classvar=p_match_class)

/* Overall differences */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model2, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 1 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model2_s1, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 2 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model2_s2, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 3 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model2_s3, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 4 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model2_s4, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 5 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model2_s5, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;


/* Variance ratio */
/* Overall */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model2, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.07, cutoff2=0.94)

proc means data=varratio_result_final;
run;

/* group 1 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model2_s1, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 2 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model2_s2, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 3 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model2_s3, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 4 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model2_s4, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 5 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model2_s5, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* for the covariates */
/* stratified by the outcome */
%PSoutcomematching(numbermi=100, data=sbirth_mi_model2, groupvar=m, propensityvar=presponse, outcomevar=dgestat, imputevar=_imputation_);  

%subset(data=sbirth_mi_model2_withmatch, classvar=outcome_match_class)

/* overall */
%mi_stddiff(numbermi=100, var=dbirwt /* fmaps */, data=sbirth_mi_model2_withmatch, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)


proc means data=stddiff_result_final;
run;

/* group 1 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model2_withmatch_s1, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 2 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model2_withmatch_s2, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 3 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model2_withmatch_s3, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 4 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model2_withmatch_s4, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 5 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model2_withmatch_s5, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;


/* Variance ratio */
/* Overall */

%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model2_withmatch, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.07, cutoff2=0.94)
proc print data=VARRATIO_RESULT_FINAL;
run;

proc means data=varratio_result_final;
run;

/* group 1 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model2_withmatch_s1, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.13, cutoff2=0.89)

proc means data=varratio_result_final;
run;

/* group 2 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model2_withmatch_s2, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.16, cutoff2=0.86)

proc means data=varratio_result_final;
run;

/* group 3 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model2_withmatch_s3, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.16, cutoff2=0.86)

proc means data=varratio_result_final;
run;

/* group 4 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model2_withmatch_s4, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.17, cutoff2=0.85)

proc means data=varratio_result_final;
run;

/* group 5 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model2_withmatch_s5, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.16, cutoff2=0.86)

proc means data=varratio_result_final;
run;


/* strata based on the propensity score */
%PSpropensitymatching(numbermi=100, data=sbirth_mi_model2, groupvar=m, propensityvar=presponse, outcomevar=dgestat, imputevar=_imputation_);  

%subset(data=sbirth_mi_model2_withmatch, classvar=pscore_match_class)

/*
%mi_stddiff(numbermi=100, var= dbirwt, data=sbirth_mi_model2_withmatch, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc print data=stddiff_result_final;
run;

proc means data=stddiff_result_final;
run;
*/

/* group 1 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model2_withmatch_s1, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 2 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model2_withmatch_s2, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 3 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model2_withmatch_s3, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 4 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model2_withmatch_s4, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 5 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model2_withmatch_s5, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;


/* Variance ratio */
/*
%mi_varratio(numbermi=100, var= dbirwt, data=sbirth_mi_model2_withmatch, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.07, cutoff2=0.94)
proc print data=VARRATIO_RESULT_FINAL;
run;

proc means data=varratio_result_final;
run;
*/

/* group 1 */
%mi_varratio(numbermi=100, var=dbirwt, data=sbirth_mi_model2_withmatch_s1, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 2 */
%mi_varratio(numbermi=100, var= dbirwt, data=sbirth_mi_model2_withmatch_s2, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 3 */
%mi_varratio(numbermi=100, var= dbirwt, data=sbirth_mi_model2_withmatch_s3, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 4 */
%mi_varratio(numbermi=100, var= dbirwt, data=sbirth_mi_model2_withmatch_s4, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 5 */
%mi_varratio(numbermi=100, var=dbirwt, data=sbirth_mi_model2_withmatch_s5, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* the end of model 2 */

/* model 3, including all continous covariates plus dbirwt square */
proc mi data=sbirth_output_missing_withmatch out=sbirth_mi_model3 seed=197789 nimpute=100;
 var  dbirwt dbirwt_square dmage  dfage fmaps wtgain  nprevis cigar  drink dgestat;
 monotone reg(dgestat= dbirwt dbirwt_square dmage  dfage fmaps wtgain  nprevis cigar  drink/ details);
 
run;


/* evaluate the standardized difference */
/* for the outcome */
/* grouping the dataset into 5 subgroups by the propensity score */

%subset(data=sbirth_mi_model3, classvar=p_match_class)

/* overall */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model3, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 1 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model3_s1, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)


proc means data=stddiff_result_final;
run;

/* group 2 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model3_s2, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 3 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model3_s3, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 4 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model3_s4, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 5 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model3_s5, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* Overall variance ratio */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model3, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.07, cutoff2=0.94)

proc means data=varratio_result_final;
run;

/* group 1 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model3_s1, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 2 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model3_s2, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 3 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model3_s3, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 4 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model3_s4, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 5 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model3_s5, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;


/* for the covariates */
/* stratified by the outcome */
%PSoutcomematching(numbermi=100, data=sbirth_mi_model3, groupvar=m, propensityvar=presponse, outcomevar=dgestat, imputevar=_imputation_);  

%subset(data=sbirth_mi_model3_withmatch, classvar=outcome_match_class)

/* overall difference */
%mi_stddiff(numbermi=100, var=dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)


proc means data=stddiff_result_final;
run;

/* group 1 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s1, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 2 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s2, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)


proc means data=stddiff_result_final;
run;

/* group 3 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s3, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 4 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s4, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 5 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s5, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;


/* Variance ratio */
/* Overall */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.07, cutoff2=0.94)

proc means data=varratio_result_final;
run;

/* group 1 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s1, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.13, cutoff2=0.89)

proc means data=varratio_result_final;
run;

/* group 2 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s2, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.16, cutoff2=0.86)

proc means data=varratio_result_final;
run;

/* group 3 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s3, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.16, cutoff2=0.86)

proc means data=varratio_result_final;
run;

/* group 4 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s4, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.17, cutoff2=0.85)

proc means data=varratio_result_final;
run;

/* group 5 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s5, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.16, cutoff2=0.86)

proc means data=varratio_result_final;
run;

/* strata based on the propensity score */
%PSpropensitymatching(numbermi=100, data=sbirth_mi_model3, groupvar=m, propensityvar=presponse, outcomevar=dgestat, imputevar=_imputation_);  


%subset(data=sbirth_mi_model3_withmatch, classvar=pscore_match_class)

/*
%mi_stddiff(numbermi=100, var= dbirwt, data=sbirth_mi_model3_withmatch, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc print data=stddiff_result_final;
run;

proc means data=stddiff_result_final;
run;
*/

/* group 1 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s1, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)


proc means data=stddiff_result_final;
run;

/* group 2 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s2, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 3 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s3, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 4 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s4, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 5 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s5, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;


/* Variance ratio */
/*

%mi_varratio(numbermi=100, var= dbirwt, data=sbirth_mi_model3_withmatch, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.07, cutoff2=0.94)
proc print data=VARRATIO_RESULT_FINAL;
run;

proc means data=varratio_result_final;
run;
*/

/* group 1 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s1, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 2 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s2, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 3 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s3, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 4 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s4, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 5 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model3_withmatch_s5, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* the end of model 3 */

/* model 4, predictive mean matching */

proc mi data=sbirth_output_missing_withmatch out=sbirth_mi_model4 seed=197789 nimpute=100;
 var  dbirwt dbirwt_square dmage  dfage fmaps wtgain  nprevis cigar  drink dgestat;
 monotone /* reg */  regpmm (dgestat= dbirwt dbirwt_square dmage  dfage fmaps wtgain  nprevis cigar  drink/ details k=10);
 
run;


/* evaluate the standardized difference */
/* for the outcome */
/* grouping the dataset into 5 subgroups by the propensity score */

%subset(data=sbirth_mi_model4, classvar=p_match_class)


%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model4, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)


proc means data=stddiff_result_final;
run;

/* group 1 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model4_s1, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)


proc means data=stddiff_result_final;
run;

/* group 2 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model4_s2, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 3 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model4_s3, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 4 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model4_s4, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 5 */
%mi_stddiff(numbermi=100, var=dgestat, data=sbirth_mi_model4_s5, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;


/* Variance ratio */
/* Overall */

%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model4, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.07, cutoff2=0.94)

proc means data=varratio_result_final;
run;

/* group 1 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model4_s1, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 2 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model4_s2, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 3 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model4_s3, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 4 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model4_s4, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 5 */
%mi_varratio(numbermi=100, var=dgestat, data=sbirth_mi_model4_s5, groupvar=m, pairvar=pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;


/* for the covariates */
/* stratified by the outcome */
%PSoutcomematching(numbermi=100, data=sbirth_mi_model4, groupvar=m, propensityvar=presponse, outcomevar=dgestat, imputevar=_imputation_);  

%subset(data=sbirth_mi_model4_withmatch, classvar=outcome_match_class)

/* for some impuations, because of the categorical values, we cannot separte into 5 categories */

/* overall difference */

%mi_stddiff(numbermi=100, var=dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 1 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s1, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 2 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s2, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 3 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s3, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 4 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s4, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 5 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s5, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;



/* Variance ratio */
/* Overall */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.07, cutoff2=0.94)

proc means data=varratio_result_final;
run;

/* group 1 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s1, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.13, cutoff2=0.89)


proc means data=varratio_result_final;
run;

/* group 2 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s2, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.16, cutoff2=0.86)

proc means data=varratio_result_final;
run;

/* group 3 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s3, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.16, cutoff2=0.86)

proc means data=varratio_result_final;
run;

/* group 4 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s4, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.17, cutoff2=0.85)

proc means data=varratio_result_final;
run;

/* group 5 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s5, groupvar=m, pairvar=outcome_pair, imputevar=_imputation_, cutoff1=1.16, cutoff2=0.86)

proc means data=varratio_result_final;
run;

/* strata based on the propensity score */
%PSpropensitymatching(numbermi=100, data=sbirth_mi_model4, groupvar=m, propensityvar=presponse, outcomevar=dgestat, imputevar=_imputation_);  


%subset(data=sbirth_mi_model4_withmatch, classvar=pscore_match_class)

/*
%mi_stddiff(numbermi=100, var= dbirwt, data=sbirth_mi_model4_withmatch, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc print data=stddiff_result_final;
run;

proc means data=stddiff_result_final;
run;
*/

/* group 1 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s1, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)


proc means data=stddiff_result_final;
run;

/* group 2 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s2, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)


proc means data=stddiff_result_final;
run;

/* group 3 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s3, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)


proc means data=stddiff_result_final;
run;

/* group 4 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s4, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)

proc means data=stddiff_result_final;
run;

/* group 5 */
%mi_stddiff(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s5, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff=0.1)


proc means data=stddiff_result_final;
run;


/* Variance ratio */
/*

%mi_varratio(numbermi=100, var= dbirwt, data=sbirth_mi_model4_withmatch, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.07, cutoff2=.94)
proc print data=VARRATIO_RESULT_FINAL;
run;

proc means data=varratio_result_final;
run;
*/

/* group 1 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s1, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 2 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s2, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 3 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s3, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 4 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s4, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* group 5 */
%mi_varratio(numbermi=100, var= dbirwt /* fmaps */, data=sbirth_mi_model4_withmatch_s5, groupvar=m, pairvar=propensity_pair, imputevar=_imputation_, cutoff1=1.15, cutoff2=0.87)

proc means data=varratio_result_final;
run;

/* the end of model 4 */
