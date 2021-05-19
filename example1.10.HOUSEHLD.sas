*********************************************************************
 MAY 30, 2017  4:46 PM
 
 This is an example of a SAS program that creates a SAS
 file from the 2016 NHIS Public Use HOUSEHLD.DAT ASCII file
 
 This is stored in HOUSEHLD.SAS
*********************************************************************;

* USER NOTE: PLACE NEXT STATEMENT IN SUBSEQUENT PROGRAMS;
LIBNAME  NHIS     "\\cdc.gov\private\A728\wdq7\imputation_diagnostics\data";

* USER NOTE: PLACE NEXT STATEMENT IN SUBSEQUENT PROGRAMS
             IF YOU ALLOW PROGRAM TO PERMANENTLY STORE FORMATS;
LIBNAME  LIBRARY  "\\cdc.gov\private\A728\wdq7\imputation_diagnostics\data";

FILENAME ASCIIDAT '\\cdc.gov\private\A728\wdq7\imputation_diagnostics\data\HOUSEHLD.dat';

* DEFINE VARIABLE VALUES FOR REPORTS;

*  USE THE STATEMENT "PROC FORMAT LIBRARY=LIBRARY"
     TO PERMANENTLY STORE THE FORMAT DEFINITIONS;

*  USE THE STATEMENT "PROC FORMAT" IF YOU DO NOT WISH
      TO PERMANENTLY STORE THE FORMATS;

PROC FORMAT LIBRARY=LIBRARY;
*PROC FORMAT;

   VALUE $GROUPC
      ' '< - HIGH   = "Range of Values"
   ;
   VALUE GROUPN
      LOW - HIGH   = "Range of Values"
   ;
   VALUE HHP001X
      10                 = "10 Household"
      20                 = "20 Person"
      25                 = "25 Income Imputation"
      30                 = "30 Sample Adult"
      38                 = "38 Functioning and Disability"
      40                 = "40 Sample Child"
      60                 = "60 Family"
      63                 = "63 Family Disability Questions"
      65                 = "65 Paradata"
      70                 = "70 Injury/Poisoning Episode"
      75                 = "75 Injury/Poisoning Verbatim"
   ;
   VALUE HHP002X
      .                   = '.'
      OTHER              = "Survey Year"
   ;
   VALUE HHP004X
      1                  = "1 Quarter 1"
      2                  = "2 Quarter 2"
      3                  = "3 Quarter 3"
      4                  = "4 Quarter 4"
   ;
   VALUE HHP005X
      01                 = "01 January"
      02                 = "02 February"
      03                 = "03 March"
      04                 = "04 April"
      05                 = "05 May"
      06                 = "06 June"
      07                 = "07 July"
      08                 = "08 August"
      09                 = "09 September"
      10                 = "10 October"
      11                 = "11 November"
      12                 = "12 December"
   ;
   VALUE HHP008X
      01                 = "01 House, apartment, flat, condo"
      02                 = "02 HU in nontransient hotel, motel"
      03                 = "03 HU-permanent in transient hotel, motel"
      04                 = "04 HU in rooming house"
      05                 = "05 Mobile home/trailer w/no permanent rooms added"
      06                 = "06 Mobile home/trailer w/1+ permanent rooms added"
      07                 = "07 HU not specified above"
      08                 = "08 Quarters not HU in room or board house"
      09                 = "09 Unit not permanent-transient hotel, motel"
      10                 = "10 Unoccupied site for mobile home/trailer/tent"
      11                 = "11 Student quarters in college dormitory"
      12                 = "12 Group quarter unit not specified above"
      98                 = "98 Not ascertained"
   ;
   VALUE HHP009X
      1                  = "1 Refused"
      2                  = "2 No one home - repeated calls"
      3                  = "3 Temporarily absent"
      4                  = "4 Language problem"
      5                  = "5 Other"
   ;
   VALUE HHP015X
      1                  = "1 Northeast"
      2                  = "2 Midwest"
      3                  = "3 South"
      4                  = "4 West"
   ;

DATA NHIS.HOUSEHLD;
   * CREATE A SAS DATA SET;
   INFILE ASCIIDAT PAD LRECL=47;

   * DEFINE LENGTH OF ALL VARIABLES;

   LENGTH

      /* IDN LENGTHS */

      RECTYPE    3   SRVY_YR    4   HHX      $ 6   INTV_QRT   3
      INTV_MON   3   WTIA_HH    8   WTFA_HH    8

      /* COV LENGTHS */

      LIVQRT     3

      /* MHH LENGTHS */

      NON_INTV   3   ACPT_FAM   3   REJ_FAM    3   ACPT_PER   3
      REJ_PER    3   ACPTCHLD   3

      /* UCF LENGTHS */

      REGION     3   PSTRAT     4   PPSU       4 ;

   * INPUT ALL VARIABLES;

   INPUT

      /* IDN LOCATIONS */

      RECTYPE       1 -   2    SRVY_YR       3 -   6
      HHX      $    7 -  12    INTV_QRT     13 -  13
      INTV_MON     14 -  15    WTIA_HH      16 -  21 .1
      WTFA_HH      22 -  27

      /* COV LOCATIONS */

      LIVQRT       28 -  29

      /* MHH LOCATIONS */

      NON_INTV     30 -  30    ACPT_FAM     31 -  32
      REJ_FAM      33 -  34    ACPT_PER     35 -  36
      REJ_PER      37 -  38    ACPTCHLD     39 -  40

      /* UCF LOCATIONS */

      REGION       41 -  41    PSTRAT       42 -  44
      PPSU         45 -  47;

   * DEFINE VARIABLE LABELS;

   LABEL

      /* IDN LABELS */

      RECTYPE    ="File type identifier"
      SRVY_YR    ="Year of National Health Interview Survey"
      HHX        ="Household Number"
      INTV_QRT   ="Interview Quarter"
      INTV_MON   ="Interview Month"
      WTIA_HH    ="Weight - Interim Annual"
      WTFA_HH    ="Weight - Final Annual"

      /* COV LABELS */

      LIVQRT     ="Type of Living Quarters"

      /* MHH LABELS */

      NON_INTV   ="Category of type A non-response"
      ACPT_FAM   ="Number of families in HH responding"
      REJ_FAM    ="Number of families in HH not responding"
      ACPT_PER   ="Number of persons in HH responding"
      REJ_PER    ="Number of persons in HH not responding"
      ACPTCHLD   ="Number of children in HH that responded"

      /* UCF LABELS */

      REGION     ="Region"
      PSTRAT     ="Pseudo-stratum for public use file variance estimation"
      PPSU       ="Pseudo-PSU for public use file variance estimation"
   ;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   *    (PUT ASTERISK (*) BEFORE WORD "FORMAT"
        IN NEXT STATEMENT TO PREVENT FORMAT
        ASSOCIATIONS BEING STORED WITH DATA SET);
   FORMAT

      /* IDN FORMAT ASSOCIATIONS */

      RECTYPE   HHP001X.   SRVY_YR   HHP002X.   HHX       $GROUPC.
      INTV_QRT  HHP004X.   INTV_MON  HHP005X.   WTIA_HH   GROUPN.
      WTFA_HH   GROUPN.

      /* COV FORMAT ASSOCIATIONS */

      LIVQRT    HHP008X.

      /* MHH FORMAT ASSOCIATIONS */

      NON_INTV  HHP009X.

      /* UCF FORMAT ASSOCIATIONS */

      REGION    HHP015X.;
RUN;

PROC CONTENTS DATA=NHIS.HOUSEHLD;
   TITLE1 'CONTENTS OF THE 2016 NHIS Household FILE';

PROC FREQ DATA=NHIS.HOUSEHLD NOTITLE;
   TABLES RECTYPE/LIST MISSING;
   TITLE1 'FREQUENCY REPORT FOR 2016 NHIS Household FILE';
   TITLE2 '(WEIGHTED)';
   WEIGHT WTFA_HH;

PROC FREQ DATA=NHIS.HOUSEHLD NOTITLE;
   TABLES RECTYPE/LIST MISSING;
   TITLE1 'FREQUENCY REPORT FOR 2016 NHIS Household FILE';
   TITLE2 '(UNWEIGHTED)';

* USER NOTE: TO SEE UNFORMATTED VALUES IN PROCEDURES, ADD THE
             STATEMENT: FORMAT _ALL_;
RUN;
