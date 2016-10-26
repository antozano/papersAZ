*****************************************************************************************;
*                                                                                       *;
* Macros to implement Lumley and Levy case-crossover approach (Environmetrics 2000)     *;
* Call the creatLumleyHazard macro to execute.                                          *;
*                                                                                       *;
* Greg Wellenius (gwelleni@hsph.harvard.edu)                                            *;
* July 14, 2003     
*revised Antonella Zanobetti
*
*;
*****************************************************************************************;
libname a "C:\Antonella\MORTALITA'";

ods listing close;*to suppress the output printing;
options nonotes nosource nosource2 ; *suppresses LOG WINDOW printing;


proc contents data=a.cataniaexp;run;



/** matching every third day **/

/* Macro def for creating the hazard data set */
%macro makecontrol (daynum=);
  data control; set cases;
    length case time 3;
    case=0; time=2;
    if day ne &daynum; /* take all days in the stratum except the index day */
    matchday=&daynum;  /* keep track of the index day that created this data set */
    matchdow=weekday(date); /* in case you want to stratify on day of week too */
	
	*extract every third day for controls;
	test=day-&daynum;
	test2=mod(test,3);
	if test2=0;
	drop test test2;run;

  data control;
    merge control (in=a drop=day) exposures (in=b drop=date);
    by stratum matchday;
    if a=1 AND b=1; 

  proc append data=control base=hazard;
%mend;* makecontrol;


%macro createLumleyHazard (poll,poll1,causa);
/** create cases and exposure data **/

   /* All exposure data (date, temp, PM, rhum, etc) should go into this data set */;
   data exposures; set a.&poll;
        day=day(data);
        month=month(data);
        yearex=year(data);
		date=mdy(month,day,yearex);
		dow=weekday(date);
		if date ne .;
     keep date dow day month yearex
  TMAXC TMEDIAC TMINC UMIDITApct VENTOMEDIAkmh VISIBILITAkm pm10 temp1 temp2 temp3 temp01 temp02 temp03 pm101 pm102 pm103 pm1001 pm1002 pm1003;
	run;
     
/* All case data (date, # of events, patient id, etc) should go into this data set */
data cases; set a.&poll1;
        day=day(dod);
        month=month(dod);
        year=year(dod);
		date=mdy(month,day,year);
		if date ne .;
		if &causa=1;
keep date id sesso eta cvd allresp tumore mi ihd chf stroke copd pneum resp;
run;

proc sort data=cases; by date;run;
proc sort data=exposures; by date;run;

   /* For air pollution stuff: change the dateto be your first date */
  %let startDate='01Jan2006'D;
   %let dateInterval='MONTH';
   %let numDays=31;


   data exposures; set exposures;
     if date >= &startDate;
     stratum = intck(&dateInterval,&startDate, date); /* Number each stratum */
     day = datdif(intnx(&dateInterval,&startDate, stratum),date,'act/act')+1; /* Number each day within each stratum */
     matchday = day;
     dow = weekday(date);


   data cases; set cases;
     if date >= &startDate;
     stratum = intck(&dateInterval,&startDate, date); /* Number each stratum */
     day = datdif(intnx(&dateInterval,&startDate, stratum),date,'act/act')+1; /* Number each day within each stratum */
     matchday = day;
     matchdow = weekday(date);

   /*proc means data=cases;*/
   proc sort data=exposures;
     by date;
   proc sort data=cases;
     by date;

   /* Start off the new data set with the cases and their actual exposure info */
   data hazard;
     merge cases (in=a) exposures (in=b);
     by date;
     if a=1 AND b=1; 
     length case time 3; /*keeps the data set smaller*/
     case=1; time=1;
     matchday=day;

   /* proc means data=hazard; */
   proc sort data=cases;
     by stratum matchday;
   proc sort data=exposures;
     by stratum matchday;

   %do i=1 %to &numDays;
     %makecontrol(daynum=&i);
   %end;

/** add in the data below all other variables or interaction you want to create */;
data a.day_&poll1&causa;set hazard;
if dow=1 then wd1=1; else wd1=0;
if dow=2 then wd2=1; else wd2=0;
if dow=3 then wd3=1; else wd3=0;
if dow=4 then wd4=1; else wd4=0;
if dow=5 then wd5=1; else wd5=0;
if dow=6 then wd6=1; else wd6=0;
run;
%mend createLumleyHazard;
%createLumleyHazard  (cataniaexp,catania,cvd);
%createLumleyHazard  (cataniaexp,catania,stroke);
%createLumleyHazard  (cataniaexp,catania,mi);
%createLumleyHazard  (cataniaexp,catania,ihd);
%createLumleyHazard  (cataniaexp,catania,chf);
%createLumleyHazard  (cataniaexp,catania,tumore);
%createLumleyHazard  (cataniaexp,catania,copd);
%createLumleyHazard  (cataniaexp,catania,pneum);
%createLumleyHazard  (cataniaexp,catania,allresp);
%createLumleyHazard  (cataniaexp,catania,resp);



/*** this is for individual data ***/

%MACRO city(haz,n,pol,out1);

/** regression  **/
proc sort data=a.&haz; by id date;run;
	ods output ParameterEstimates=mypar2(keep= parameter Estimate  StdErr);
	proc phreg data=a.&haz nosummary ;
  	model time*case(0) = &pol  tmediac wd1-wd6/ties=discrete;
    	strata  id date; 
		run;

	data est0; set mypar2;
	reg=&n;
	tval=Estimate/StdErr;
	retain recount;
	recount+1;
	if recount<3;
	run;	
	proc append data=est0 base=a.catania_&out1;run; 
%mend;

%city (Day_cataniacvd,1,pm10,cvdexp);
%city (Day_cataniacvd,2,pm101,cvdexp);
%city (Day_cataniacvd,3,pm102,cvdexp);
%city (Day_cataniacvd,4,pm103,cvdexp);
%city (Day_cataniacvd,5,pm1001,cvdexp);
%city (Day_cataniacvd,6,pm1002,cvdexp);
%city (Day_cataniacvd,7,pm1003,cvdexp);


%city (Day_cataniastroke,1,pm10,strokeexp);
%city (Day_cataniastroke,2,pm101,strokeexp);
%city (Day_cataniastroke,3,pm102,strokeexp);
%city (Day_cataniastroke,4,pm103,strokeexp);
%city (Day_cataniastroke,5,pm1001,strokeexp);
%city (Day_cataniastroke,6,pm1002,strokeexp);
%city (Day_cataniastroke,7,pm1003,strokeexp);

%city (Day_cataniami,1,pm10,miexp);
%city (Day_cataniami,2,pm101,miexp);
%city (Day_cataniami,3,pm102,miexp);
%city (Day_cataniami,4,pm103,miexp);
%city (Day_cataniami,5,pm1001,miexp);
%city (Day_cataniami,6,pm1002,miexp);
%city (Day_cataniami,7,pm1003,miexp);


%city (Day_cataniatumore,1,pm10,tumoreexp);
%city (Day_cataniatumore,2,pm101,tumoreexp);
%city (Day_cataniatumore,3,pm102,tumoreexp);
%city (Day_cataniatumore,4,pm103,tumoreexp);
%city (Day_cataniatumore,5,pm1001,tumoreexp);
%city (Day_cataniatumore,6,pm1002,tumoreexp);
%city (Day_cataniatumore,7,pm1003,tumoreexp);


%city (Day_cataniaallresp,1,pm10,allrespexp);
%city (Day_cataniaallresp,2,pm101,allrespexp);
%city (Day_cataniaallresp,3,pm102,allrespexp);
%city (Day_cataniaallresp,4,pm103,allrespexp);
%city (Day_cataniaallresp,5,pm1001,allrespexp);
%city (Day_cataniaallresp,6,pm1002,allrespexp);
%city (Day_cataniaallresp,7,pm1003,allrespexp);


%city (Day_cataniaresp,1,pm10,respexp);
%city (Day_cataniaresp,2,pm101,respexp);
%city (Day_cataniaresp,3,pm102,respexp);
%city (Day_cataniaresp,4,pm103,respexp);
%city (Day_cataniaresp,5,pm1001,respexp);
%city (Day_cataniaresp,6,pm1002,respexp);
%city (Day_cataniaresp,7,pm1003,respexp);


%city (Day_cataniacopd,1,pm10,copdexp);
%city (Day_cataniacopd,2,pm101,copdexp);
%city (Day_cataniacopd,3,pm102,copdexp);
%city (Day_cataniacopd,4,pm103,copdexp);
%city (Day_cataniacopd,5,pm1001,copdexp);
%city (Day_cataniacopd,6,pm1002,copdexp);
%city (Day_cataniacopd,7,pm1003,copdexp);


%city (Day_cataniapneum,1,pm10,pneumexp);
%city (Day_cataniapneum,2,pm101,pneumexp);
%city (Day_cataniapneum,3,pm102,pneumexp);
%city (Day_cataniapneum,4,pm103,pneumexp);
%city (Day_cataniapneum,5,pm1001,pneumexp);
%city (Day_cataniapneum,6,pm1002,pneumexp);
%city (Day_cataniapneum,7,pm1003,pneumexp);


%city (Day_cataniaihd,1,pm10,ihdexp);
%city (Day_cataniaihd,2,pm101,ihdexp);
%city (Day_cataniaihd,3,pm102,ihdexp);
%city (Day_cataniaihd,4,pm103,ihdexp);
%city (Day_cataniaihd,5,pm1001,ihdexp);
%city (Day_cataniaihd,6,pm1002,ihdexp);
%city (Day_cataniaihd,7,pm1003,ihdexp);

%city (Day_cataniachf,1,pm10,chfexp);
%city (Day_cataniachf,2,pm101,chfexp);
%city (Day_cataniachf,3,pm102,chfexp);
%city (Day_cataniachf,4,pm103,chfexp);
%city (Day_cataniachf,5,pm1001,chfexp);
%city (Day_cataniachf,6,pm1002,chfexp);
%city (Day_cataniachf,7,pm1003,chfexp);


data cvd;set a.Catania_cvdexp;
outcome='cvd';
run;
data chf;set a.Catania_chfexp;
outcome='chf';
run;
data ihd;set a.Catania_ihdexp;
outcome='ihd';
run;
data tumore;set a.Catania_tumoreexp;
outcome='tumore';
run;
data copd;set a.Catania_copdexp;
outcome='copd';
run;
data pneum;set a.Catania_pneumexp;
outcome='pneum';
run;
data allresp;set a.Catania_allrespexp;
outcome='allresp';
run;
data resp;set a.Catania_respexp;
outcome='resp';
run;
data stroke;set a.Catania_strokeexp;
outcome='stroke';
run;
data mi;set a.Catania_miexp;
outcome='mi';
run;
data a.rescatapm10; set tumore  mi stroke cvd ihd chf resp allresp copd pneum;run;



%MACRO city(haz,n,pol,out1);

/** regression  **/
proc sort data=a.&haz; by id date;run;
	ods output ParameterEstimates=mypar2(keep= parameter Estimate  StdErr);
	proc phreg data=a.&haz nosummary ;
  	model time*case(0) = &pol pm10 wd1-wd6/ties=discrete;
    	strata  id date; 
		run;

	data est0; set mypar2;
	reg=&n;
	tval=Estimate/StdErr;
	retain recount;
	recount+1;
	if recount<3;
	run;	
	proc append data=est0 base=a.catania_&out1;run; 
%mend;


%city (Day_cataniastroke,1,tmediac,stroketemp);
%city (Day_cataniastroke,2,temp1,stroketemp);
%city (Day_cataniastroke,3,temp2,stroketemp);
%city (Day_cataniastroke,4,temp3,stroketemp);
%city (Day_cataniastroke,5,temp01,stroketemp);
%city (Day_cataniastroke,6,temp02,stroketemp);
%city (Day_cataniastroke,7,temp03,stroketemp);


%city (Day_cataniastroke,8,tmediac,stroketemp);
%city (Day_cataniastroke,9,temp1,stroketemp);
%city (Day_cataniastroke,10,temp2,stroketemp);
%city (Day_cataniastroke,11,temp3,stroketemp);
%city (Day_cataniastroke,12,temp01,stroketemp);
%city (Day_cataniastroke,13,temp02,stroketemp);
%city (Day_cataniastroke,14,temp03,stroketemp);
