#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <time.h>
#include <string.h>


#define NX 720           /* 360 or 720 */
#define NY 360           /* 180 or 360 */
#define NEWNT 102         /* 101/102, 61/62,61,61,61,61,56,56  */
#define NCCODE 192       /* number of countries  139 or 192 */
#define NZ 31            /* Number of z dims tracked 16 or 31 (z is distance from managed cells) */
#define MAXZ 21           /* 11 or 21 Maximum z before we get tired, and spread remaining harvest over all forested cells with z >= this value */ 
#define BEST_CASE 1               /* 1= best case, 0= other case */
#define BEST_CASE_MIN_FLOWS_T4 0    /* 1= best_sage or best_hyde or */
#define BEST_CASE_MIN_FLOWS_T5 0 
                                    /* best_sage3=same as best_sage except set l0 or other instead of l1 */
                                    /* best_hyde3=same as best_hyde except set l0 or other instead of l1 */

                                    /* 0= other case or */
                                    /* best_sage2=same as best_sage except set t1 or other instead of t4 */
                                    /* best_hyde2=same as best_hyde except set t1 or other instead of t4 */

#define COUNTRY_TEST 0              /* 0=global, 1=country */  
#define COUNTRY_TEST_GCODE 11        /* 221=US, 34=CANADA, 28=BRAZIL, CAMBODIA=35, AUS=11 */
#define COUNTRY_TEST_INDEX 6         /* 126=US, 20=CANADA, 16=BRAZIL, CAMBODIA=21, AUS=6 */
#define CONTERMINOUS 0                /* 1=CONTERMINOUS US ONLY,0=ALL US */
  
#define TIME_CHECK 0                 /* 1=Restricted # of time steps */
#define TC 2                           /* # of time steps desired */    

#define TOTAL_HARVEST_SWITCH 1
#define SECONDARY_HARVEST_SWITCH 1
#define VIRGIN_HARVEST_SWITCH 1
#define FORCE_HARVEST_SWITCH 1
#define CPAVG 1
#define ORIGINAL_DATA 0               /* 1=print out original states not using
                                         transitions, 0=print out states being
                                         updated with transitions */

#define TB2BB 1.0         /* total biomass to bole biomass ratio, value                                                     of 2.0 assumes the WH numbers are bole biomass only 
                              and we need to cut twice as much biomass  */

#define PROB_FNAME "phbio.average.7states.txt"
#define NPHB 50            /* filelength of probability of harvest given biomass */

#define SMART_FLOW_BUG_PRINT 1
#define STATE_BUG_PRINT 1
#define FLOW_BUG_PRINT 1

#define ZEROVALUE 0.000000000000000000000
#define ZEROVALUE_CHECK 0.000000001


#define NREG 21
#define HEADLEN 6
#define CLASS_NUM 6      /*c,p,o,i,w,u*/
#define FLEN 14
#define PATH1 "/data/" 

#define CAT01 "backup/projects/glm/inputs/hyde_1.0/1deg_grids/"
#define CAT0 "backup/projects/glm/inputs/sage/1deg_grids/"
#define CAT1 "backup/projects/glm/inputs/hyde_1.0/1deg_grids/nodata/"
#define CAT2 "backup/projects/glm/inputs/hyde_1.0/1deg_grids/future/"
#define CAT3 "backup/projects/glm/code/public_inputs/other/"
#define CAT4 "backup/projects/glm/inputs/hyde_3.0/1deg_grids/"
#define CAT5 "backup/projects/glm/code/public_inputs/hyde_3.0/half_deg_grids/"
#define CAT6 "backup/projects/glm/code/public_inputs/hyde_3.0/half_deg_grids/nodata/"

int yrstart;
int future_rate_option, nodata_option, logging_option, res_option, hist_option;
char trun[10], finput_dir[150];


char fnames[CLASS_NUM][FLEN] = { 
{"gcrop.files."},
{"gpast.files."},		
{"gothr.files."},
{"gicew.files."},
{"gwatr.files."},
{"gurbn.files."},
};

struct new_data {
  int zdis;
  double vbh, sbh, sma, smb, flowvs, sbh2, vbh2, sbh3;
  double c,p,v,i,w,s,u;  
  double flowcp, flowpc, flowpv, flowvp, flowvc, flowcv;
  double flowsp,flowps,flowsc,flowcs; 
  double flowcu, flowpu, flowvu, flowsu, flowuc, flowup, flowus; 
} newdata[NY][NX][NEWNT];  

struct czt_data {
  double vb;                  /* units MgC */
} cztdata[NCCODE][NZ][NEWNT];

struct ct_data {
  double wh, whr, stbh, vwh;  /* units MgC */
  double cavg, pavg, cgarea, pgarea, cpavg, cpgarea;
  double wh_at_zmax, wh_lessthan_zmax;
  double smb, vnfb, sarea, smb_nf, sarea_nf;
  double converted_forest_land;
  double converted_forest_land_total;
  double flowvc_prime,flowsc_prime,flowcs_prime,flowvp_prime,flowsp_prime,flowps_prime;
  double flowvu_prime,flowsu_prime,flowus_prime;
  int pcount;
} ctdata[NCCODE][NEWNT];

struct c_data {
  double area, areaf, areanf;
  int ccode, newccode, continent_code, rcode;
  int smart_flow_option, converted_forest_land_option;
  int zdis_option, adjust_smart_flow_option, harvest_option;
} cdata[NCCODE];

struct r_data {
  int rcode, continent_code;
} rdata[NREG];

struct other_data {
  int usflag, gcode, newgcode, gcontinent_code, fnf, rcode, newrcode, shiftcult;
  double vba, vnppa;
} dstatic[NY][NX];


char infilelist[CLASS_NUM][NEWNT][120];
char header[HEADLEN][25];
double lat[NY], garea[NY][NX], lon[NX], phb[NPHB], check_total[NEWNT], check_total3[NEWNT];
int year[NEWNT];
char cname[NCCODE][50], rname[NREG][50];

   
struct sum_check {
  double c,p,v,i,w,s,u;
  double vf,vnf,sf,snf;
  double smaf,smanf,vbh,sbh,vbh2,sbh2,sbh3,wh_unmet;
  double total, wh, cavg, pavg, cpavg;  
  double flowcp, flowpc, flowpv, flowvp, flowvc, flowcv;
  double flowsp,flowps,flowsc,flowcs;
  double flowcu,flowpu,flowvu,flowsu,flowuc,flowup,flowus;
  double flowsbh, flowsbh2, flowvbh, flowvbh2, flowsbh3;
  double converted_forest_land, converted_forest_land_total;
  double flowvc_prime,flowsc_prime,flowcs_prime,flowvp_prime,flowsp_prime,flowps_prime;
  double flowvu_prime,flowsu_prime,flowus_prime;
} byclass_sum[NEWNT];



struct check_data {
  double c,p,v,s,u;
} chtest[NY][NX][NEWNT];  


void get_flisting();
void initialize();
void initialize_checker(int country_code);
void option_settings();
void read_data_pastruns();
void read_other_data();
void cellinfo();
void read_country_names();
void read_country_codes();
void read_regional_codes();
void read_regional_names();
void read_continent_codes();
void update_vb(int it);
void transitions();
void secondary_harvest(int it, int i);
void virgin_harvest(int it, int i);
void force_harvest2(int it, int i);
void country_cpavg(int it, int i);
void smart_flow(int k, int m, int it);
void alternative_smart_flow1(int k, int m, int it);
void alternative_smart_flow2(int k, int m, int it);
void adjust_smart_flow1(int k, int m, int it, int i);
void adjust_smart_flow2(int k, int m, int it, int i);
void adjust_smart_flow3(int k, int m, int it, int i);
void update_states(int it, int zmax);
void zdis_calc(int it);
float prob_harv(float biomass);
void country_timeseries_checker(int country_code, char country_name[50]);
void loop_call_for_country_final_stats();
void country_final_stats(int country_code, int istart);
void global_timeseries_checker(int country_code, char country_name[50]);
void continent_timeseries_checker(int continent_code, char continent_name[50]);
void regional_timeseries_checker(int regional_code, char regional_name[50]);
void output_updated_states();
void output_updated_states2();
void output_updated_states3();
void compute_country_flows(int i, int it);
void output_lu();
void country_primeflow_print();

/*****************************************************************************/
int main(int ac, char *av[]){

  FILE *outfile;
  int i, country_code, continent_code, regional_code;
  char country_name[50], continent_name[50], regional_name[50];


  /* arguments brought in from shell script are set to global variables and
     used throughout */

  setbuf(stdout, NULL);
 
  strcpy(trun,av[1]);
  future_rate_option = atoi(av[2]);
  nodata_option=atoi(av[3]);
  logging_option=atoi(av[4]);
  strcpy(finput_dir,av[5]);
  printf("finput_dir %s\n", finput_dir); 
  res_option= atoi(av[6]);
  hist_option = atoi(av[7]);

  printf("res_option = %d\n",res_option);

  printf("\n\nPROGRAM: 1700-2000 GLM\n\n");


  option_settings();

  if(hist_option==1){
    if(strcmp("tone",trun) == 0) {yrstart=1700;}
    else if(strcmp("ttwo",trun) == 0) {yrstart=1761;}
    else if(strcmp("tthree",trun) == 0) {yrstart=1822;}
    else if(strcmp("tfour",trun) == 0) {yrstart=1883;}
    else if(strcmp("tfive",trun) == 0) {yrstart=1944;}
    else if(strcmp("tsix",trun) == 0) {yrstart=2005;}
    else {yrstart=2050;}
  }
  else if(hist_option==2){
    if(strcmp("tone",trun) == 0) {yrstart=1500;}
    else if(strcmp("ttwo",trun) == 0) {yrstart=1601;}
    else if(strcmp("tthree",trun) == 0) {yrstart=1702;}
    else if(strcmp("tfour",trun) == 0) {yrstart=1803;}
    else if(strcmp("tfive",trun) == 0) {yrstart=1904;}
    else if(strcmp("tsix",trun) == 0) {yrstart=2005;}
    else {yrstart=2050;}
  }
  else if(hist_option==3){
    if(strcmp("tone",trun) == 0) {yrstart=1850;}
    else if(strcmp("ttwo",trun) == 0) {yrstart=1881;}
    else if(strcmp("tthree",trun) == 0) {yrstart=1912;}
    else if(strcmp("tfour",trun) == 0) {yrstart=1943;}
    else if(strcmp("tfive",trun) == 0) {yrstart=1974;}
    else if(strcmp("tsix",trun) == 0) {yrstart=2005;}
    else {yrstart=2050;}
  }


#if 1
  initialize(); 
  read_data_pastruns(); 
#endif


#if 1

  cellinfo();

#if 1
  read_country_names();
  read_country_codes();
  read_regional_names();
  read_regional_codes();
  read_continent_codes();
  read_other_data();
#endif
#if 1  
  transitions(); 
#endif

#if 1
  if(strcmp("tone",trun) == 0) {
    outfile=fopen("check_total.txt","w");
  }
  else{
    outfile=fopen("check_total.txt","a");
  }
  for (i=1;i<NEWNT;i++){
    fprintf(outfile,"total at time %d is %lf\n",i,check_total[i]);
  }
  fclose(outfile);
#endif

#if 1
  if(strcmp("tone",trun) == 0) {
    outfile=fopen("check_total3.txt","w");
  }
  else{
    outfile=fopen("check_total3.txt","a");
  }
  for (i=1;i<NEWNT;i++){
    fprintf(outfile,"total at time %d is %lf\n",i,check_total3[i]);
  }
  fclose(outfile);
#endif

#if 1
  for (i=0;i<NREG;i++){
    if(rdata[i].rcode > 0 && rdata[i].rcode < 18){
      strcpy(regional_name,rname[i]);
      regional_code=rdata[i].rcode;
      regional_timeseries_checker(regional_code,regional_name);
    }
  }
#endif


#if 1

  if(res_option==1){
  strcpy(country_name,"aus");
  country_code=11;
  country_timeseries_checker(country_code,country_name);

  strcpy(country_name,"us");
  country_code=221;
  country_timeseries_checker(country_code,country_name);

  strcpy(country_name,"ak");
  country_code=260;
  country_timeseries_checker(country_code,country_name);

  strcpy(country_name,"br");
  country_code=28;
  country_timeseries_checker(country_code,country_name);

  strcpy(country_name,"cn");
  country_code=34;
  country_timeseries_checker(country_code,country_name); 
  }
  if(res_option==2){
  strcpy(country_name,"aus");
  country_code=36;
  country_timeseries_checker(country_code,country_name);

  strcpy(country_name,"us");
  country_code=840;
  country_timeseries_checker(country_code,country_name);

  strcpy(country_name,"ak");
  country_code=900;
  country_timeseries_checker(country_code,country_name);

  strcpy(country_name,"br");
  country_code=76;
  country_timeseries_checker(country_code,country_name);

  strcpy(country_name,"cn");
  country_code=156;
  country_timeseries_checker(country_code,country_name); 
  }
  strcpy(country_name,"global");
#if COUNTRY_TEST
  country_code=COUNTRY_TEST_GCODE;
#else
  country_code=1000;
#endif
  global_timeseries_checker(country_code,country_name); 


#endif


#if 1
  strcpy(continent_name,"na");
  continent_code=1;
  continent_timeseries_checker(continent_code,continent_name);

  strcpy(continent_name,"sa");
  continent_code=2;
  continent_timeseries_checker(continent_code,continent_name);

  strcpy(continent_name,"eu");
  continent_code=3;
  continent_timeseries_checker(continent_code,continent_name);

  strcpy(continent_name,"asia");
  continent_code=4;
  continent_timeseries_checker(continent_code,continent_name);

  strcpy(continent_name,"africa");
  continent_code=5;
  continent_timeseries_checker(continent_code,continent_name);

  strcpy(continent_name,"aus");
  continent_code=6;
  continent_timeseries_checker(continent_code,continent_name);

  loop_call_for_country_final_stats();

#endif


#if 1
  output_updated_states();
#endif
#if 1
  output_updated_states2();
#endif
#if 1
  output_updated_states3();
#endif

#if 0
  output_lu(); 
#endif

  country_primeflow_print();

#endif
 
  return;
}

/***********************************************************************/
void option_settings(){

  FILE *infile;
  char new_path[90], tmp_path0[90], tmp_path1[90], fname[90];
  int i;


  strcpy(fname,"clist.txt");

  infile=fopen(fname,"r");

  for (i=0;i<NCCODE;i++){

    fscanf(infile,"%d %d %d %d %d\n",&cdata[i].smart_flow_option,&cdata[i].converted_forest_land_option,&cdata[i].zdis_option,&cdata[i].adjust_smart_flow_option,&cdata[i].harvest_option);

  }
  fclose(infile);


}
/***********************************************************************/
void initialize(){

  int i, iz, k, m, it;
  
 

  year[0]=yrstart;
  for (it=1;it<NEWNT;it++) year[it]=year[0]+it;


  for (it=0;it<NEWNT;it++){
    
    for (i=0;i<NCCODE;i++) {


      cdata[i].area=ZEROVALUE;
      cdata[i].areaf=ZEROVALUE;
      cdata[i].areanf=ZEROVALUE;

      ctdata[i][it].wh=ZEROVALUE;
      ctdata[i][it].whr=ZEROVALUE;
      ctdata[i][it].stbh=ZEROVALUE;
      ctdata[i][it].vwh=ZEROVALUE;
      
      ctdata[i][it].wh_lessthan_zmax=ZEROVALUE;
      ctdata[i][it].wh_at_zmax=ZEROVALUE;
      ctdata[i][it].vnfb=ZEROVALUE; 
      ctdata[i][it].smb=ZEROVALUE;
      ctdata[i][it].sarea=ZEROVALUE;
      ctdata[i][it].sarea_nf=ZEROVALUE;
      ctdata[i][it].smb_nf=ZEROVALUE;
      ctdata[i][it].pavg=ZEROVALUE;  
      ctdata[i][it].cavg=ZEROVALUE;  
      ctdata[i][it].cpavg=ZEROVALUE;  
      ctdata[i][it].cgarea=ZEROVALUE;
      ctdata[i][it].pgarea=ZEROVALUE;
      ctdata[i][it].cpgarea=ZEROVALUE;
      ctdata[i][it].pcount=0;
      ctdata[i][it].converted_forest_land=ZEROVALUE;
      ctdata[i][it].converted_forest_land_total=ZEROVALUE;
      ctdata[i][it].flowvc_prime=ZEROVALUE;
      ctdata[i][it].flowsc_prime=ZEROVALUE;
      ctdata[i][it].flowcs_prime=ZEROVALUE;
      ctdata[i][it].flowvp_prime=ZEROVALUE;
      ctdata[i][it].flowsp_prime=ZEROVALUE;
      ctdata[i][it].flowps_prime=ZEROVALUE;
      ctdata[i][it].flowus_prime=ZEROVALUE;
      ctdata[i][it].flowsu_prime=ZEROVALUE;
      ctdata[i][it].flowvu_prime=ZEROVALUE;

      for (iz=0;iz<NZ;iz++){
	cztdata[i][iz][it].vb=ZEROVALUE;
      }
    }
    
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){


	newdata[k][m][it].c=ZEROVALUE;
	newdata[k][m][it].p=ZEROVALUE;
	newdata[k][m][it].v=ZEROVALUE;
	newdata[k][m][it].s=ZEROVALUE;
	newdata[k][m][it].i=ZEROVALUE;
	newdata[k][m][it].w=ZEROVALUE;
	newdata[k][m][it].u=ZEROVALUE;
	newdata[k][m][it].flowcp=ZEROVALUE;
	newdata[k][m][it].flowpc=ZEROVALUE;
	newdata[k][m][it].flowpv=ZEROVALUE;
	newdata[k][m][it].flowvp=ZEROVALUE;
	newdata[k][m][it].flowvc=ZEROVALUE;
	newdata[k][m][it].flowcv=ZEROVALUE;
	newdata[k][m][it].flowsc=ZEROVALUE;
	newdata[k][m][it].flowcs=ZEROVALUE;
	newdata[k][m][it].flowsp=ZEROVALUE;
	newdata[k][m][it].flowps=ZEROVALUE;
	newdata[k][m][it].flowvs=ZEROVALUE;
	newdata[k][m][it].flowcu=ZEROVALUE;
	newdata[k][m][it].flowpu=ZEROVALUE;
	newdata[k][m][it].flowsu=ZEROVALUE;
	newdata[k][m][it].flowvu=ZEROVALUE;
	newdata[k][m][it].flowuc=ZEROVALUE;
	newdata[k][m][it].flowup=ZEROVALUE;
	newdata[k][m][it].flowus=ZEROVALUE;
	newdata[k][m][it].zdis=NZ-1;
        newdata[k][m][it].sma=ZEROVALUE;
        newdata[k][m][it].smb=ZEROVALUE;
        newdata[k][m][it].sbh2=ZEROVALUE;
        newdata[k][m][it].sbh3=ZEROVALUE;
        newdata[k][m][it].vbh2=ZEROVALUE;
	newdata[k][m][it].vbh=ZEROVALUE;
	newdata[k][m][it].sbh=ZEROVALUE;    
	dstatic[k][m].gcode=0;
        dstatic[k][m].newgcode=9999; 
	dstatic[k][m].gcontinent_code=0;
	dstatic[k][m].fnf=0;
	dstatic[k][m].vba=ZEROVALUE;
	dstatic[k][m].vnppa=ZEROVALUE;
	dstatic[k][m].rcode=0;
	dstatic[k][m].newrcode=0;
	dstatic[k][m].shiftcult=0;

      }
    }
  }
  
  return;
}
     
/***********************************************************************/
void get_flisting(){

  FILE *infile;
  int i, j, cnum;
  char new_path[90], in_cat[90]; 
  char tmp_path0[90], tmp_path1[90], tmp_path2[90], tmp_path3[90], tmp_path4[90];
  char class_fname[90], fname[90];
  
 
  printf("\ngetting file listings...\n");


  strcpy(new_path,PATH1);


  for (i=0;i<CLASS_NUM;i++){
    
    if(nodata_option == 6) strcpy(in_cat,CAT6);
    if(nodata_option == 5) strcpy(in_cat,CAT5); 
    if(nodata_option == 4) strcpy(in_cat,CAT4);     
    if(nodata_option == 3) strcpy(in_cat,CAT0);  
    if(nodata_option == 2) strcpy(in_cat,CAT1);  
    if(nodata_option == 1) strcpy(in_cat,CAT01); 
  
    strcpy(tmp_path0,new_path);
    strcat(tmp_path0,in_cat);
    strcpy(tmp_path1,tmp_path0);
    
    if ((nodata_option ==5) || (nodata_option == 4) || (nodata_option == 6)){
      if (hist_option == 1) strcpy(tmp_path4,"1700-2005/urban/");
      if (hist_option == 2) strcpy(tmp_path4,"1500-2005/urban/");
      if (hist_option == 3) strcpy(tmp_path4,"1850-2005/urban/");

      strcat(tmp_path1,tmp_path4);
    }

    strcat(tmp_path1,fnames[i]);
    strcpy(class_fname,tmp_path1);
    strcat(class_fname,trun);

    infile=fopen(class_fname,"r");
    
    for (j=0;j<NEWNT;j++) {
      

      fscanf(infile,"%s\n",fname);
      
      strcpy(tmp_path2,new_path);
      strcat(tmp_path2,in_cat);
      strcpy(tmp_path3,tmp_path2);
if (nodata_option ==6){
      if (hist_option == 1) strcpy(tmp_path4,"1700-2005/");
      if (hist_option == 2) strcpy(tmp_path4,"1500-2005/");
      if (hist_option == 3) strcpy(tmp_path4,"1850-2005/");

      strcat(tmp_path3,tmp_path4);
    }
      strcat(tmp_path3,"urban/");
      strcat(tmp_path3,fname);
      
      strcpy(infilelist[i][j],tmp_path3);
      
      
    }
    fclose(infile);

  }

  return;
}
/*******************************************************************/

void read_data_pastruns(){

  FILE *infile;
  int ic, i, k, m, it, cnum;
  char newpath[90], tmppath[90], ftag1[90], ftag2[90], fin[90];
  

  get_flisting(); 

  printf("finputdir: %s\n",finput_dir); 
  printf("\nreading data for runs prior/up to 2000...\n");

  
  for (ic=0;ic<CLASS_NUM;ic++){   

 
    for (it=0;it<NEWNT;it++){ 
     
    infile=fopen(infilelist[ic][it],"r");
 
#if 1
    for (i=0;i<HEADLEN;i++) fgets(header[i],25,infile);
#endif
  
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){   
	  if(ic == 0) fscanf(infile,"%lf\n",&newdata[k][m][it].c);
	  if(ic == 1) fscanf(infile,"%lf\n",&newdata[k][m][it].p);
	  if(ic == 2) fscanf(infile,"%lf\n",&newdata[k][m][it].v);
	  if(ic == 3) fscanf(infile,"%lf\n",&newdata[k][m][it].i);
	  if(ic == 4) fscanf(infile,"%lf\n",&newdata[k][m][it].w);  
	  if(ic == 5) fscanf(infile,"%lf\n",&newdata[k][m][it].u);  
      } 
    }
    fclose(infile);
    
    }
  }


  if(strcmp("tone",trun) != 0){
    
    printf("%s\n",trun);

    strcpy(newpath,PATH1);
    strcpy(tmppath,newpath);
    strcat(tmppath,finput_dir);
    strcat(tmppath,"updated_states/");
 
    cnum=7; 
    for (ic=0;ic<cnum;ic++){
      
       if(ic == 0) strcpy(ftag1,"gothr.");
       if(ic == 1) strcpy(ftag1,"gsecd."); 
       if(ic == 2) strcpy(ftag1,"gcrop.");
       if(ic == 3) strcpy(ftag1,"gpast."); 
       if(ic == 4) strcpy(ftag1,"gssmb."); 
       if(ic == 5) strcpy(ftag1,"gssma."); 
       if(ic == 6) strcpy(ftag1,"gurbn."); 

       if(hist_option==1){
	 if(strcmp("ttwo",trun) == 0) strcat(ftag1,"1761.txt");
	 if(strcmp("tthree",trun) == 0) strcat(ftag1,"1822.txt"); 
	 if(strcmp("tfour",trun) == 0) strcat(ftag1,"1883.txt");
	 if(strcmp("tfive",trun) == 0) strcat(ftag1,"1944.txt");
       }
       else if(hist_option==2){
	 if(strcmp("ttwo",trun) == 0) strcat(ftag1,"1601.txt");
	 if(strcmp("tthree",trun) == 0) strcat(ftag1,"1702.txt"); 
	 if(strcmp("tfour",trun) == 0) strcat(ftag1,"1803.txt");
	 if(strcmp("tfive",trun) == 0) strcat(ftag1,"1904.txt");
       }
       else if(hist_option==3){
	 if(strcmp("ttwo",trun) == 0) strcat(ftag1,"1881.txt");
	 if(strcmp("tthree",trun) == 0) strcat(ftag1,"1912.txt"); 
	 if(strcmp("tfour",trun) == 0) strcat(ftag1,"1943.txt");
	 if(strcmp("tfive",trun) == 0) strcat(ftag1,"1974.txt");
       }

       strcpy(ftag2,ftag1);
       strcpy(fin,tmppath);
       strcat(fin,ftag2);
    
   
       infile=fopen(fin,"r");
       
#if 1
       for (i=0;i<HEADLEN;i++) fgets(header[i],25,infile);
#endif


       for (k=0;k<NY;k++){
	 for (m=0;m<NX;m++){
	   
	   if(ic == 0){ 
	     fscanf(infile,"%lf\n",&newdata[k][m][0].v);
	   }
	   if(ic == 1){ 
	     fscanf(infile,"%lf\n",&newdata[k][m][0].s);
	   }
	   if(ic == 2){ 
	     fscanf(infile,"%lf\n",&newdata[k][m][0].c);
	   }
	   if(ic == 3){ 
	     fscanf(infile,"%lf\n",&newdata[k][m][0].p);
	   } 
	   if(ic == 4){
	     fscanf(infile,"%lf\n",&newdata[k][m][0].smb); 
	   } 
	   if(ic == 5) fscanf(infile,"%lf\n",&newdata[k][m][0].sma); 
	   if(ic == 6){
	     fscanf(infile,"%lf\n",&newdata[k][m][0].u);
	   }
      	 }
       }
       fclose(infile);
       
    }
  }




  
  return;
}
/********************************************************************/
void cellinfo(){

  FILE *infile;
  char new_path[90], tmp_path0[90], tmp_path1[90], tmp_path4[90], new_path2[90], in_cat[90];
  int k, m, it;
  float er=6370., pi=3.141592654;
  float nwloncorner=-180., nwlatcorner=90., dtocenter, dstep;
  double summer=ZEROVALUE, summer2=ZEROVALUE;

 
  
  if (res_option==1){
    dtocenter=-0.5;
    dstep=-1.0;
  }
  else if (res_option==2){
    dtocenter=-0.25;
    dstep=-0.5;
  }

  lat[0]=nwlatcorner+dtocenter;
  for (k=1;k<NY;k++){
    lat[k]=lat[k-1]+dstep;
  }
  
  
  if (res_option==1){
    dtocenter=0.5;
    dstep=1.0;
  }
  else if (res_option==2){
    dtocenter=0.25;
    dstep=0.5;
  }    
 
  lon[0]=nwloncorner+dtocenter;
  for (m=1;m<NX;m++){
     lon[m]=lon[m-1]+dstep;
  }
  
  
   strcpy(new_path,PATH1);
  
   strcpy(in_cat,CAT3); 
   strcpy(tmp_path0,new_path);
   strcat(tmp_path0,in_cat);
   strcpy(tmp_path1,tmp_path0);
   strcpy(tmp_path4,"cellarea/");
   strcat(tmp_path1,tmp_path4);
   strcpy(new_path2,tmp_path1);

   if (res_option==2){
     strcat(new_path2,"cellarea_halfdeg.txt");
   } 

  infile=fopen(new_path2,"r");
  
  
  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){ 
      
    fscanf(infile,"%lf\n",&garea[k][m]);
    
    garea[k][m]=garea[k][m]*1000.00*1000.00;   /* convert cell area from km2 to m2 */
    
    }
  } 
  fclose(infile);
  
 

  
  return;
}



/********************************************************************/
void read_country_codes(){

  FILE *cinfile;
  int i, k, m, it;
  char new_path[90], tmp_path0[90], tmp_path1[90];
  char gfname[90], ffname[90], outfname[90], scfname[90];


  printf("\nreading country codes...\n");
  
  
  strcpy(new_path,PATH1);
  strcpy(tmp_path0,new_path);
  strcat(tmp_path0,CAT3);
  strcpy(tmp_path1,tmp_path0);
  if ((nodata_option == 6) || (nodata_option == 5)) {
    strcat(tmp_path1,"ccodes/ccodes.txt.sort2wh");
  }
  else {
    strcat(tmp_path1,"ccodes/add_ak/ccodes.txt.sort2wh_139countries");
  }
  strcpy(ffname,tmp_path1);
  
  cinfile=fopen(ffname,"r");
  for (i=0;i<NCCODE;i++) fscanf(cinfile,"%d\n",&cdata[i].ccode);
  fclose(cinfile);
  
  if (res_option==2){
    strcpy(new_path,PATH1);
    strcpy(tmp_path0,new_path);
    strcat(tmp_path0,CAT3);
    strcpy(tmp_path1,tmp_path0);
    strcat(tmp_path1,"shift_cult/shiftcult_map_halfdeg.txt");
    strcpy(scfname,tmp_path1);
  
    cinfile=fopen(scfname,"r");
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){
	fscanf(cinfile,"%d\n",&dstatic[k][m].shiftcult);
      }
    }
    fclose(cinfile);
  }

  strcpy(new_path,PATH1);
  strcpy(tmp_path0,new_path);
  strcat(tmp_path0,CAT3);
  strcpy(tmp_path1,tmp_path0);


  if (res_option==1){
    if ((nodata_option == 6) || (nodata_option == 5)){
      strcat(tmp_path1,"ccodes/add_ak/gcodes_ak_conform.asc");
    }
    else {
      strcat(tmp_path1,"ccodes/add_ak/gcodes_ak_conform.asc");
    }
  }
  else if (res_option==2){
    strcat(tmp_path1,"ccodes/ccodes_half_deg.txt");
  } 

  strcpy(gfname,tmp_path1);

  cinfile=fopen(gfname,"r");


  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){
      fscanf(cinfile,"%d\n",&dstatic[k][m].gcode);
    }
  }
  fclose(cinfile);
  
  
  
  for (i=0;i<NCCODE;i++) {
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){
	
	if(dstatic[k][m].gcode == cdata[i].ccode) dstatic[k][m].newgcode=i;
	
      }
    }
    cdata[i].newccode=i;
  }

  
 
  return;
}
/********************************************************************/
void read_regional_codes(){

  FILE *cinfile;
  int i, k, m, it;
  int ccodetmp[NCCODE], rcodetmp[NCCODE];
  char new_path[90], tmp_path0[90], tmp_path1[90];
  char gfname[90], ffname[90], outfname[90];



  printf("\nreading regional codes...\n");
  

  strcpy(new_path,PATH1);
  strcpy(tmp_path0,new_path);
  strcat(tmp_path0,CAT3);
  strcpy(tmp_path1,tmp_path0);
  strcat(tmp_path1,"wood_harvest/codes.txt");
  strcpy(ffname,tmp_path1);
  
  cinfile=fopen(ffname,"r");
  for (i=0;i<NREG;i++) fscanf(cinfile,"%d\n",&rdata[i].rcode);
  fclose(cinfile);
  

  
  strcpy(new_path,PATH1);
  strcpy(tmp_path0,new_path);
  strcat(tmp_path0,CAT3);
  strcpy(tmp_path1,tmp_path0);
  strcat(tmp_path1,"wood_harvest/continent_codes.txt");
  strcpy(ffname,tmp_path1);
  
  cinfile=fopen(ffname,"r");
  for (i=0;i<NREG;i++) fscanf(cinfile,"%d\n",&rdata[i].continent_code);
  fclose(cinfile);
  



  strcpy(new_path,PATH1);
  strcpy(tmp_path0,new_path);
  strcat(tmp_path0,CAT3);
  strcpy(tmp_path1,tmp_path0);
  strcat(tmp_path1,"wood_harvest/codes2glm.txt");
  strcpy(ffname,tmp_path1);
  

  cinfile=fopen(ffname,"r");
  for (i=0;i<NCCODE;i++) fscanf(cinfile,"%d %d\n",&rcodetmp[i],&ccodetmp[i]);
  fclose(cinfile);
  


  for (i=0;i<NCCODE;i++) {
   for (k=0;k<NCCODE;k++) {
	
     if(cdata[i].ccode == ccodetmp[k]) cdata[i].rcode=rcodetmp[k];
     
   }
  }
 
  

  for (i=0;i<NCCODE;i++) {
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){
	
	if(dstatic[k][m].gcode == cdata[i].ccode) dstatic[k][m].rcode=cdata[i].rcode;

	
      }
    }
  }
 


  for (i=0;i<NREG;i++) {
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){
	
	if(dstatic[k][m].rcode == rdata[i].rcode) dstatic[k][m].newrcode=i;
	
      }
    }
  }

  
 
  return;
}

/********************************************************************/
void read_continent_codes(){

  FILE *cinfile;
  int i, k, m, it;
  char new_path[90], tmp_path0[90], tmp_path1[90];
  char gfname[90], ffname[90], outfname[90];



  printf("\nreading continent codes...\n");
  
  
  strcpy(new_path,PATH1);
  strcpy(tmp_path0,new_path);
  strcat(tmp_path0,CAT3);
  strcpy(tmp_path1,tmp_path0);
  strcat(tmp_path1,"ccodes/continent.codes.txt.sort2wh");
  strcpy(ffname,tmp_path1);
  
  cinfile=fopen(ffname,"r");
  for (i=0;i<NCCODE;i++) fscanf(cinfile,"%d\n",&cdata[i].continent_code);
  fclose(cinfile);
  
  
  strcpy(new_path,PATH1);
  strcpy(tmp_path0,new_path);
  strcat(tmp_path0,CAT3);
  strcpy(tmp_path1,tmp_path0);

  if (res_option==1){
    strcat(tmp_path1,"ccodes/add_ak/gcodes_continent_conform.asc");
  }
  else if (res_option==2){
    strcat(tmp_path1,"ccodes/gcodes_continent_half_deg_DUMMY.asc");
  } 
  
  strcpy(gfname,tmp_path1);

  cinfile=fopen(gfname,"r");


  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){
      fscanf(cinfile,"%d\n",&dstatic[k][m].gcontinent_code);
    }
  }
  fclose(cinfile);
  
  return;
}
/********************************************************************/
void read_other_data(){

  FILE *infile, *testfile;
  int i, k, m, it;
  double sumtest=ZEROVALUE, summer, wh_multiplier=1.0;
  double tester=ZEROVALUE; 
  char new_path[90], tmp_path0[90], tmp_path1[90], ffname[90];


  printf("\nreading other data...\n");

  /* static biomass grid, initial units=kgC/m2 */
  
  strcpy(new_path,PATH1);
  strcpy(tmp_path0,new_path);
  strcat(tmp_path0,CAT3);
  strcpy(tmp_path1,tmp_path0);

  if (res_option==1){
    strcat(tmp_path1,"miami_biomass_v3/miami_1deg_conform.txt");
  }
  else if (res_option==2){
    strcat(tmp_path1,"miami_biomass_v3/miami_halfdeg_conform.txt");
  } 

  
  strcpy(ffname,tmp_path1);
  
 infile=fopen(ffname,"r");

  

  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){

      fscanf(infile,"%lf\n",&dstatic[k][m].vba);

      /* putting carbon in terms of aboveground w/ 0.75 */
       dstatic[k][m].vba=dstatic[k][m].vba*0.75; 

      /* to not divide by zero for flowvbh */
       if(dstatic[k][m].vba < 0.01) dstatic[k][m].vba=0.01; 

      /* resetting forest definition to MLU */
      if(dstatic[k][m].vba >= 2.0){
	dstatic[k][m].fnf=1;
      }
      else{
	dstatic[k][m].fnf=0;
      }

    }
  }
  fclose(infile);
  printf("finished reading vba\n");


  /* grid of vnpp from miami model in columnar format, initial units=kgC/m2/yr */

  strcpy(new_path,PATH1);
  strcpy(tmp_path0,new_path);
  strcat(tmp_path0,CAT3);
  strcpy(tmp_path1,tmp_path0);

  if (res_option==1){
    strcat(tmp_path1,"miami_npp/miami.dat.in_conform");
  }
  else if (res_option==2){
    strcat(tmp_path1,"miami_npp/miami.half_deg.in_conform");
  } 


  strcpy(ffname,tmp_path1);
  
  infile=fopen(ffname,"r");
  
  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){
      fscanf(infile,"%lf\n",&dstatic[k][m].vnppa);
    }
  }
  fclose(infile);
  
  printf("finished reading npp\n");

  /* wood harvest by country thru time, initial units=MgC */
  
  strcpy(new_path,PATH1);
  strcpy(tmp_path0,new_path);
  strcat(tmp_path0,CAT3);
  strcpy(tmp_path1,tmp_path0);


  if(hist_option==1){
    if ((nodata_option==5) || (nodata_option==6)){
      if(logging_option == 4) {
	strcat(tmp_path1,"wood_harvest/nodata/1700-2005/woodharvin.txt.in.nodata.");
      }
      else if(logging_option == 1) {
	strcat(tmp_path1,"wood_harvest/1700-2005/woodharvin.txt.in.");
      }
      else if(logging_option == 2){
	strcat(tmp_path1,"wood_harvest/1700-2005/woodharvin.txt.in.");
	wh_multiplier=1.5;
      }
      else if(logging_option == 3){
	strcat(tmp_path1,"wood_harvest/1700-2005/woodharvin.txt.in."); 
	wh_multiplier=0.5;
      }
      else if(logging_option == 0){
	strcat(tmp_path1,"wood_harvest/1700-2005/woodharvin.txt.in."); 
	wh_multiplier=0;
      }
    }
    else {
      if(logging_option == 4) {
	strcat(tmp_path1,"wood_harvest/nodata/woodharvin.txt.in.nodata.");
      }
      else if(logging_option == 1) {
	strcat(tmp_path1,"wood_harvest/1700-2000/Old/woodharvin.txt.in.");
      }
      else if(logging_option == 2){
	strcat(tmp_path1,"wood_harvest/1700-2000/Old/woodharvin.txt.in.");
	wh_multiplier=1.5;
      }
      else if(logging_option == 3){
	strcat(tmp_path1,"wood_harvest/1700-2000/Old/woodharvin.txt.in."); 
	wh_multiplier=0.5;
      }
      else if(logging_option == 0){
	strcat(tmp_path1,"wood_harvest/1700-2000/Old/woodharvin.txt.in."); 
	wh_multiplier=0;
      }
    }
  }
  else if(hist_option==2){
    if(logging_option == 4) {
      strcat(tmp_path1,"wood_harvest/nodata/1500-2005/woodharvin.txt.in.nodata.");
    }
    else if(logging_option == 1) {
      strcat(tmp_path1,"wood_harvest/1500-2005/woodharvin.txt.in.");
    }
    else if(logging_option == 2){
      strcat(tmp_path1,"wood_harvest/1500-2005/woodharvin.txt.in.");
      wh_multiplier=1.5;
    }
    else if(logging_option == 3){
      strcat(tmp_path1,"wood_harvest/1500-2005/woodharvin.txt.in."); 
      wh_multiplier=0.5;
    }
    else if(logging_option == 0){
      strcat(tmp_path1,"wood_harvest/1500-2005/woodharvin.txt.in."); 
      wh_multiplier=0;
    }
  }
  else if(hist_option==3){
    if(logging_option == 4) {
      strcat(tmp_path1,"wood_harvest/nodata/1850-2005/woodharvin.txt.in.");
    }
    else if(logging_option == 1) {
      strcat(tmp_path1,"wood_harvest/1850-2005/woodharvin.txt.in.");
    }
    else if(logging_option == 2){
      strcat(tmp_path1,"wood_harvest/1850-2000/woodharvin.txt.in.");
      wh_multiplier=1.5;
    }
    else if(logging_option == 3){
      strcat(tmp_path1,"wood_harvest/1850-2000/woodharvin.txt.in."); 
      wh_multiplier=0.5;
    }
    else if(logging_option == 0){
      strcat(tmp_path1,"wood_harvest/1850-2005/woodharvin.txt.in."); 
      wh_multiplier=0;
    }
  }

  strcpy(ffname,tmp_path1);
  strcat(ffname,trun);
  

  infile=fopen(ffname,"r");


  
  for (it=0;it<NEWNT;it++){
    for (i=0;i<NCCODE;i++) {
      
      fscanf(infile,"%lf\n",&ctdata[i][it].wh);     
      ctdata[i][it].wh*=TB2BB; 
      ctdata[i][it].wh*=wh_multiplier;

      if(logging_option == 0) ctdata[i][it].wh=ZEROVALUE;

     }
   }
  fclose(infile);



  /* read probability of harvest biomass for secondary */

  strcpy(ffname,PROB_FNAME);
  infile=fopen(ffname,"r");

  for (i=0;i<NPHB;i++) fscanf(infile,"%*s %lf\n",&phb[i]);      
  fclose(infile);



  return;
}


/********************************************************************/
void read_country_names(){

  FILE *cinfile;
  int i;
  char new_path[90], tmp_path0[90], tmp_path1[90], ffname[90];
  

  printf("\nreading country names...\n");
 
  strcpy(new_path,PATH1);
  strcpy(tmp_path0,new_path);
  strcat(tmp_path0,CAT3);
  strcpy(tmp_path1,tmp_path0);

  if ((nodata_option == 6) || (nodata_option == 5)){
    strcat(tmp_path1,"ccodes/cnames.txt.sort2wh");
  }
  else {
    strcat(tmp_path1,"ccodes/add_ak/cnames.txt.sort2wh_139countries");
  }
  strcpy(ffname,tmp_path1);
  
  cinfile=fopen(ffname,"r");
  for (i=0;i<NCCODE;i++) {
    fscanf(cinfile,"%s\n",&cname[i]);
  }
  fclose(cinfile);

  return;
}
/********************************************************************/
void read_regional_names(){

  FILE *cinfile;
  int i;
  char new_path[90], tmp_path0[90], tmp_path1[90], ffname[90];
  

  printf("\nreading regional names...\n");
 
  strcpy(new_path,PATH1);
  strcpy(tmp_path0,new_path);
  strcat(tmp_path0,CAT3);
  strcpy(tmp_path1,tmp_path0);
  strcat(tmp_path1,"wood_harvest/names.txt");
  strcpy(ffname,tmp_path1);
  
  cinfile=fopen(ffname,"r");
  for (i=0;i<NREG;i++) {
    fscanf(cinfile,"%s\n",&rname[i]);
    printf("%s\n",rname[i]);

  }
  fclose(cinfile);

 

  return;

}

/********************************************************************/
void country_primeflow_print(){

  FILE *testfile;  
  int i, it;
 

  testfile=fopen("country.primeflow.txt","w");


  for (it=0;it<NEWNT;it++){
    for (i=0;i<NCCODE;i++){
     
      fprintf(testfile,"cname %s vc_prime %lf sc_prime %lf cs_prime %lf vp_prime %lf sp_prime %lf ps_prime %lf us_prime %lf su_prime %lf vu_prime %lf\n",cname[i],ctdata[i][it].flowvc_prime,ctdata[i][it].flowsc_prime,ctdata[i][it].flowcs_prime,ctdata[i][it].flowvp_prime,ctdata[i][it].flowsp_prime,ctdata[i][it].flowps_prime,ctdata[i][it].flowus_prime,ctdata[i][it].flowsu_prime,ctdata[i][it].flowvu_prime);
      
    }
  }

  fclose(testfile);
  
  return;

}
	    
	

/********************************************************************/
/** newdata= yearly interpolated data where it=0 is 1700
    transitions are computed here where it=1 is transition from 1700 to 1701
**/

void transitions(){
  
  char outstat[2];
  int iz, i, ic, k, m, it, ih, zmax=0, izz;
  double tester1, tester2, tester3, tester4;
  


  printf("\ncomputing transitions...\n");



  /* wh checker files written to below */

  if(strcmp("tone",trun) == 0) {
    strcpy(outstat,"w");
  }
  else {
    strcpy(outstat,"a");
  }



#if TIME_CHECK
  
  for (it=0;it<TC;it++){
#else
    
    for (it=0;it<NEWNT-1;it++){

  check_total3[it]=0;
#endif
            
   
      for (k=0;k<NY;k++){
	for (m=0;m<NX;m++){


#if COUNTRY_TEST
 
	  if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){

#else	   
	    if(dstatic[k][m].gcode > 0){
#endif	      
	      i=dstatic[k][m].newgcode;

#if BEST_CASE
	      if((dstatic[k][m].gcontinent_code == 3) || (dstatic[k][m].gcontinent_code == 4)){
		
		/* minimum flows, secondary priority, Eurasia only */ 
		
		alternative_smart_flow2(k,m,it);
		
#if BEST_CASE_MIN_FLOWS_T4    /* abandonment in the tropical Eurasia, secondary priority */
		if((lat[k] <= 23.5) && (lat[k] >= -23.5)) adjust_smart_flow2(k,m,it,i);
#endif

#if BEST_CASE_MIN_FLOWS_T5   /*abandonment in shiftcult map */
		if(dstatic[k][m].shiftcult==1) adjust_smart_flow2(k,m,it,i);
#endif
		
	      }
	      else{
		
		/* minimum flows everywhere with abandonment in the tropics, primary priority */
		
		alternative_smart_flow1(k,m,it);
		
#if BEST_CASE_MIN_FLOWS_T4    /* abandonment in the tropical Eurasia, primary priority */
		if((lat[k] <= 23.5) && (lat[k] >= -23.5)) adjust_smart_flow1(k,m,it,i);
#endif
#if BEST_CASE_MIN_FLOWS_T5   /*abandonment in shiftcult map */
		if(dstatic[k][m].shiftcult==1) adjust_smart_flow1(k,m,it,i);
#endif
		
	      }


#else

	      if(cdata[i].smart_flow_option == 1){
	
		if(cdata[i].adjust_smart_flow_option == 1) {
		  alternative_smart_flow1(k,m,it);
		  
		}
		else if(cdata[i].adjust_smart_flow_option == 2) {
		  alternative_smart_flow1(k,m,it);	     
		  adjust_smart_flow1(k,m,it,i);
		  
		  
		}
		else if(cdata[i].adjust_smart_flow_option == 3){
		  alternative_smart_flow1(k,m,it);
		  
		  if((lon[m]>-106.0) && (lon[m]<-71.0) && (lat[k]<20.0) && (lat[k]>5.0)) adjust_smart_flow3(k,m,it,i);
		  if((lon[m]>-81.0) && (lon[m]<-60.0) && (lat[k]<5.0) && (lat[k]>-25.0)) adjust_smart_flow3(k,m,it,i);
		  if((lon[m]>-60.0) && (lon[m]<-40.0) && (lat[k]<5.0) && (lat[k]>-10.0)) adjust_smart_flow3(k,m,it,i);
		  if((lon[m]>-20.0) && (lon[m]<55.0) && (lat[k]<15.0) && (lat[k]>-15.0)) adjust_smart_flow3(k,m,it,i);
		  if((lon[m]>90.0) && (lon[m]<105.0) && (lat[k]<30.0) && (lat[k]>-11.0)) adjust_smart_flow3(k,m,it,i);
		  if((lon[m]>105.0) && (lon[m]<160.0) && (lat[k]<20.0) && (lat[k]>-11.0)) adjust_smart_flow3(k,m,it,i);
		  
		}
		else if(cdata[i].adjust_smart_flow_option ==4){  

		  
		  alternative_smart_flow1(k,m,it);
		  
		  /* additional abandonment in tropics */
		  if((lat[k] <= 23.5) && (lat[k] >= -23.5)) adjust_smart_flow1(k,m,it,i);
		  
		}
		else if (cdata[i].adjust_smart_flow_option == 5){
		  alternative_smart_flow1(k,m,it);
	    if(dstatic[k][m].gcode>0){
	      check_total3[it]+=(newdata[k][m][it].v+newdata[k][m][it].s+newdata[k][m][it].c+newdata[k][m][it].p+newdata[k][m][it].u+newdata[k][m][it].i+newdata[k][m][it].w)*garea[k][m];
	    }
		  if(dstatic[k][m].shiftcult==1) adjust_smart_flow1(k,m,it,i);

		}
		
		
	      } /*end smart_flow if */
	      
	      else if(cdata[i].smart_flow_option == 2){
		
		if(cdata[i].adjust_smart_flow_option == 1){
		  alternative_smart_flow2(k,m,it);
		  
		}
		else if(cdata[i].adjust_smart_flow_option == 2) {	 
		  alternative_smart_flow2(k,m,it);
		  adjust_smart_flow2(k,m,it,i);
		  
		}
		else if(cdata[i].adjust_smart_flow_option == 3) {	
		  alternative_smart_flow2(k,m,it);
		  
		  if((lon[m]>-106.0) && (lon[m]<-71.0) && (lat[k]<20.0) && (lat[k]>5.0)) adjust_smart_flow3(k,m,it,i);
		  if((lon[m]>-81.0) && (lon[m]<-60.0) && (lat[k]<5.0) && (lat[k]>-25.0)) adjust_smart_flow3(k,m,it,i);
		  if((lon[m]>-60.0) && (lon[m]<-40.0) && (lat[k]<5.0) && (lat[k]>-10.0)) adjust_smart_flow3(k,m,it,i);
		  if((lon[m]>-20.0) && (lon[m]<55.0) && (lat[k]<15.0) && (lat[k]>-15.0)) adjust_smart_flow3(k,m,it,i);
		  if((lon[m]>90.0) && (lon[m]<105.0) && (lat[k]<30.0) && (lat[k]>-11.0)) adjust_smart_flow3(k,m,it,i);
		  if((lon[m]>105.0) && (lon[m]<160.0) && (lat[k]<20.0) && (lat[k]>-11.0)) adjust_smart_flow3(k,m,it,i);
		  
		}
		else if(cdata[i].adjust_smart_flow_option == 4) {
		  alternative_smart_flow2(k,m,it);
		  
		  /* additional abandonment in tropics */
		  if((lat[k] <= 23.5) && (lat[k] >= -23.5)) adjust_smart_flow2(k,m,it,i);
		}

		else if(cdata[i].adjust_smart_flow_option == 5) {
		  alternative_smart_flow2(k,m,it);

	    if(dstatic[k][m].gcode>0){
	      check_total3[it]+=(newdata[k][m][it].v+newdata[k][m][it].s+newdata[k][m][it].c+newdata[k][m][it].p+newdata[k][m][it].u+newdata[k][m][it].i+newdata[k][m][it].w)*garea[k][m];
	    }

		  if(dstatic[k][m].shiftcult==1) adjust_smart_flow2(k,m,it,i);
		}
		
		
	      } /*end smart_flow else if */
	      
	      else{ /* cdata[i].smart_flow_option == 3 */
		smart_flow(k,m,it);
		/* need an adjusted smart flow for abandonment */
		
	      } /*end smart_flow else if */
	      
	      
#endif /* end of BEST_CASE if */
	      
	      
	      
	    } /* end of gcode if */
	    
	    
	    if(dstatic[k][m].fnf == 1){
	      
	      ctdata[i][it].converted_forest_land+=(((newdata[k][m][it].flowvp+newdata[k][m][it].flowvc+newdata[k][m][it].flowvu)*dstatic[k][m].vba)+((newdata[k][m][it].flowsp+newdata[k][m][it].flowsc+newdata[k][m][it].flowsu)*newdata[k][m][it].smb))*garea[k][m]/1000.0; 	
	      
	    }
	    
	    
	    
	  
	    
	  }  /* end of m */
	}  /* end of k */
	
	
	
    
#if TOTAL_HARVEST_SWITCH

    printf("beginning harvest, time= %d\n",year[it]); 

    zdis_calc(it);

    update_vb(it);


/* country loop begin */

#if COUNTRY_TEST
   for(i=COUNTRY_TEST_INDEX;i<COUNTRY_TEST_INDEX+1;i++){ 
#else
      for (i=0;i<NCCODE;i++){
#endif

	
	ctdata[i][it].stbh=ZEROVALUE;

	

	ctdata[i][it].whr=ctdata[i][it].wh;     /* units = MgC */ 




	/* priority for wood clearing */

#if BEST_CASE


	if((cdata[i].continent_code == 3) || (cdata[i].continent_code == 4)){

	  /* if Eurasia, for h1 do not count clearing in harvest, but track its total */
          /* if Eurasia, for h3 count clearing in harvest and track its total */


	  if(nodata_option == 1){

	    ctdata[i][it].converted_forest_land_total+=ctdata[i][it].converted_forest_land;
	    ctdata[i][it].converted_forest_land=ZEROVALUE;

	  }
	  if(nodata_option == 3){

	    if(ctdata[i][it].converted_forest_land <= ctdata[i][it].wh){
	      ctdata[i][it].converted_forest_land_total+=ctdata[i][it].converted_forest_land;
	      ctdata[i][it].whr-=ctdata[i][it].converted_forest_land;
	    }
	    else{
	      ctdata[i][it].whr=ZEROVALUE;
	      ctdata[i][it].converted_forest_land_total+=ctdata[i][it].converted_forest_land;
	      ctdata[i][it].converted_forest_land=ctdata[i][it].wh;
	    }

	  }
	  if(nodata_option == 4){

	    ctdata[i][it].converted_forest_land_total+=ctdata[i][it].converted_forest_land;
	    ctdata[i][it].converted_forest_land=ZEROVALUE;
	  }
	  if((nodata_option == 6) || (nodata_option == 5)){
	    ctdata[i][it].converted_forest_land_total+=ctdata[i][it].converted_forest_land; 
	    ctdata[i][it].converted_forest_land=ZEROVALUE; 
	  }

	}
	else{

	  /* if not Eurasia, do not count clearing in harvest, but track its total */

	  ctdata[i][it].converted_forest_land_total+=ctdata[i][it].converted_forest_land;
	  ctdata[i][it].converted_forest_land=ZEROVALUE;

	}


#else


	/* checks for counting converted land in wood harvest */

	if(cdata[i].converted_forest_land_option == 2){
	  if(ctdata[i][it].converted_forest_land <= ctdata[i][it].wh){
            ctdata[i][it].converted_forest_land_total+=ctdata[i][it].converted_forest_land;
	    ctdata[i][it].whr-=ctdata[i][it].converted_forest_land;
	  }
	  else{
	    ctdata[i][it].whr=ZEROVALUE;
	    ctdata[i][it].converted_forest_land_total+=ctdata[i][it].converted_forest_land;
    	    ctdata[i][it].converted_forest_land=ctdata[i][it].wh;
	  }
	}
	else{
	  ctdata[i][it].converted_forest_land_total+=ctdata[i][it].converted_forest_land;
	  ctdata[i][it].converted_forest_land=ZEROVALUE;
	}

#endif



	



	/* priority for wood harvest */

#if BEST_CASE


	/* if best case, Eurasia does secondary priority for wood, non-Eurasia does primay */

	if((cdata[i].continent_code == 3) || (cdata[i].continent_code == 4)){	


#if SECONDARY_HARVEST_SWITCH 
	  secondary_harvest(it,i);
#endif 	  	  

	  ctdata[i][it].vwh = ctdata[i][it].whr;

#if VIRGIN_HARVEST_SWITCH	        
	  virgin_harvest(it,i);
#endif

	}
	else{	 


	  ctdata[i][it].vwh = ctdata[i][it].whr;
	  
#if VIRGIN_HARVEST_SWITCH	        
	  virgin_harvest(it,i);
#endif
	  
#if SECONDARY_HARVEST_SWITCH 
	  secondary_harvest(it,i);
#endif 

	}


#else




	  /********************************/
	if(cdata[i].harvest_option == 2){

#if SECONDARY_HARVEST_SWITCH 
	  secondary_harvest(it,i);

#endif 	  	  
	  
	  ctdata[i][it].vwh = ctdata[i][it].whr;
	  
#if VIRGIN_HARVEST_SWITCH	        
	  virgin_harvest(it,i);

#endif
	  
	}
	else{

	  ctdata[i][it].vwh = ctdata[i][it].whr;

#if VIRGIN_HARVEST_SWITCH	        
	  virgin_harvest(it,i);

#endif

#if SECONDARY_HARVEST_SWITCH 
	  secondary_harvest(it,i);

#endif 

	}
	  /********************************/

#endif  /* end if for BEST_CASE if */	






#if FORCE_HARVEST_SWITCH	
      	if(ctdata[i][it].whr > 0.01) force_harvest2(it,i);
#endif



    
	 } /* end country */

	
    
#endif   /* end total_harvest_switch */
     
      
 
	 update_states(it,zmax); 



#if CPAVG
	 for(i=0;i<NCCODE;i++) country_cpavg(it,i);
#endif



      
    }  /* end time */
    
    
       
 
    return;
    
  }

/***********************************************************************/
void secondary_harvest(int it, int i){

int k, m;
double sbh;


  /* calculate secondary biomass harvest for each grid cell (sbh), and
     track the accumulating amount of sbh at the country level (stbh)  */



        for (k=0;k<NY;k++){
          for (m=0;m<NX;m++){


#if COUNTRY_TEST

            if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){

#else
              if(dstatic[k][m].gcode > 0){
#endif


                if(dstatic[k][m].gcode == cdata[i].ccode){



                   if(dstatic[k][m].fnf == 1){

                    if(ctdata[i][it].whr > ZEROVALUE){

                      sbh=newdata[k][m][it].smb*
                        prob_harv(newdata[k][m][it].smb)*
                        newdata[k][m][it].s*garea[k][m];


                      if(sbh <= ctdata[i][it].whr*1000.){
                        newdata[k][m][it].sbh=sbh;
                      }
                      else{
                        newdata[k][m][it].sbh=ctdata[i][it].whr*1000.;
                      }

		      

                    } /*end of whr */


                      ctdata[i][it].stbh+=newdata[k][m][it].sbh/1000.;


                      ctdata[i][it].whr-=newdata[k][m][it].sbh/1000.;


                  } /* end of fnf */
                }  /* end of gcode */
              }


            } /* end of m */
          } /* end of k */

         return;

}

/***********************************************************************/
void virgin_harvest(int it, int i){
                      
  FILE *testfile;
  int k, m, iz, izz, zmax=0, im, j;
  double total_avail;
  char outstat[2];



  if((strcmp("tone",trun) == 0) && (it==0)) {
    strcpy(outstat,"w");
  }
  else {
    strcpy(outstat,"a");
  }
  testfile=fopen("zmax.test",outstat);


    /* first determine zmax; the maximum # cells away from the focal cell (the agricultural
       cell) needed to go to attain wood harvest demand.
       zmax=0; have enough biomass in focal (agricultural) cell 
       zmax=1; need to go to the next adjacent cell 
       zmax=2; etc...zmax=10 or MAXZ (MAXZ=11; maximum z before we get tired, and then 
       we spread remaining harvest over all forested cells with iz >= this value */




        if(ctdata[i][it].vwh < cztdata[i][0][it].vb)
          zmax=0;

        else if(ctdata[i][it].vwh < cztdata[i][0][it].vb+cztdata[i][1][it].vb){
          zmax=1;
        }
        else if(ctdata[i][it].vwh < cztdata[i][0][it].vb+cztdata[i][1][it].vb
                +cztdata[i][2][it].vb){
          zmax=2;
        }
        else if(ctdata[i][it].vwh < cztdata[i][0][it].vb+cztdata[i][1][it].vb
                +cztdata[i][2][it].vb
                +cztdata[i][3][it].vb){
          zmax=3;
        }
        else if(ctdata[i][it].vwh < cztdata[i][0][it].vb+cztdata[i][1][it].vb
                +cztdata[i][2][it].vb
                +cztdata[i][3][it].vb+cztdata[i][4][it].vb){
          zmax=4;
        }
        else if(ctdata[i][it].vwh < cztdata[i][0][it].vb+cztdata[i][1][it].vb
                +cztdata[i][2][it].vb
                +cztdata[i][3][it].vb+cztdata[i][4][it].vb+cztdata[i][5][it].vb){
          zmax=5;
        }
        else if(ctdata[i][it].vwh < cztdata[i][0][it].vb+cztdata[i][1][it].vb
                +cztdata[i][2][it].vb
                +cztdata[i][3][it].vb+cztdata[i][4][it].vb+cztdata[i][5][it].vb+
                cztdata[i][6][it].vb){
          zmax=6;
        }
        else if(ctdata[i][it].vwh < cztdata[i][0][it].vb+cztdata[i][1][it].vb
                +cztdata[i][2][it].vb
                +cztdata[i][3][it].vb+cztdata[i][4][it].vb+cztdata[i][5][it].vb+
                cztdata[i][6][it].vb+cztdata[i][7][it].vb){
          zmax=7;
        }
        else if(ctdata[i][it].vwh < cztdata[i][0][it].vb+cztdata[i][1][it].vb
                +cztdata[i][2][it].vb
                +cztdata[i][3][it].vb+cztdata[i][4][it].vb+cztdata[i][5][it].vb+
                cztdata[i][6][it].vb+cztdata[i][7][it].vb+cztdata[i][8][it].vb){
          zmax=8;
        }
        else if(ctdata[i][it].vwh < cztdata[i][0][it].vb+cztdata[i][1][it].vb
                +cztdata[i][2][it].vb
                +cztdata[i][3][it].vb+cztdata[i][4][it].vb+cztdata[i][5][it].vb+
                cztdata[i][6][it].vb+cztdata[i][7][it].vb+cztdata[i][8][it].vb
                +cztdata[i][9][it].vb){
          zmax=9;
        }
        else if(ctdata[i][it].vwh < cztdata[i][0][it].vb+cztdata[i][1][it].vb
                +cztdata[i][2][it].vb
                +cztdata[i][3][it].vb+cztdata[i][4][it].vb+cztdata[i][5][it].vb+
                cztdata[i][6][it].vb+cztdata[i][7][it].vb+cztdata[i][8][it].vb
                +cztdata[i][9][it].vb+cztdata[i][10][it].vb){
          zmax=10;
        }
        else {
          zmax=MAXZ;
        }




        if(zmax == 0){
	  ctdata[i][it].wh_at_zmax=ctdata[i][it].whr;
	}
        else{
	  for(j=0;j<zmax;j++) ctdata[i][it].wh_lessthan_zmax+=cztdata[i][j][it].vb;
	  ctdata[i][it].wh_at_zmax=ctdata[i][it].whr-ctdata[i][it].wh_lessthan_zmax;
	}
       

	
        /***************/
        /* determine virgin biomass harvest (vbh) in units of kgC and the flow
           of virgin to secondary (flowvs) */


        for (k=0;k<NY;k++){
          for (m=0;m<NX;m++){


#if COUNTRY_TEST

          if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){
#else
            if(dstatic[k][m].gcode > 0){
#endif


              if(dstatic[k][m].gcode == cdata[i].ccode) {

                if(ctdata[i][it].whr > ZEROVALUE) {


		  iz=newdata[k][m][it].zdis;


                   if(dstatic[k][m].fnf == 1){



                    /* determine how to take biomass for harvest demand */

                    if(zmax<MAXZ){ /*business as usual*/


                      if(iz < zmax){ /* take as much as possible to fulfill harvest demand */

                        newdata[k][m][it].flowvs=(newdata[k][m][it].v-newdata[k][m][it].flowvc-newdata[k][m][it].flowvp-newdata[k][m][it].flowvu);
			
                        newdata[k][m][it].vbh=newdata[k][m][it].flowvs*garea[k][m]*dstatic[k][m].vba;


                      }
                      else if(iz == zmax){

                        /*if iz=zmax; take what is needed proportionally out of all cells from 
                          focal cell out iz=zmax cells away */



                        if(cztdata[i][iz][it].vb > ZEROVALUE){

                          newdata[k][m][it].flowvs=(newdata[k][m][it].v-
                              newdata[k][m][it].flowvc-newdata[k][m][it].flowvp-newdata[k][m][it].flowvu)*
                              ctdata[i][it].wh_at_zmax/cztdata[i][iz][it].vb;

			 

                          newdata[k][m][it].vbh=newdata[k][m][it].flowvs*garea[k][m]*
                                                 dstatic[k][m].vba;
                        }
                      }





                    } /* end of business as usual, zmax < MAXZ */

                    else{ /* zmax>=MAXZ */


                      if(iz < MAXZ){ /*zmax=MAXZ*/
			/*most dry countries with reported wood harvest and not enough forest*/
			/* take everything possible */

                        newdata[k][m][it].flowvs=(newdata[k][m][it].v-
                           newdata[k][m][it].flowvc-newdata[k][m][it].flowvp-newdata[k][m][it].flowvu);

			

                        newdata[k][m][it].vbh=newdata[k][m][it].flowvs*garea[k][m]*
                          dstatic[k][m].vba;


                      }
                      else if(iz<=NZ){ /*from MAXZ=zmax to NZ*/

			/* rare, reported wood harvest, but have to go more than MAXZ to get it*/
			/* take what is needed proportionally out of each cell */

                        total_avail=ZEROVALUE;

                        for(izz=MAXZ;izz<NZ;izz++)
                          total_avail+=cztdata[i][izz][it].vb;

                        if(total_avail > ZEROVALUE){
     
                           if(ctdata[i][it].wh_at_zmax <= total_avail){

			     newdata[k][m][it].flowvs=(newdata[k][m][it].v-
						       newdata[k][m][it].flowvc-newdata[k][m][it].flowvp-newdata[k][m][it].flowvu)*ctdata[i][it].wh_at_zmax/total_avail;
			   }
			   else{
			     newdata[k][m][it].flowvs=(newdata[k][m][it].v-
					       newdata[k][m][it].flowvc-newdata[k][m][it].flowvp-newdata[k][m][it].flowvu);

			   }
			}

			/*safe, but still may not get enough harvest*/

			

                        newdata[k][m][it].vbh=newdata[k][m][it].flowvs*garea[k][m]*
                          dstatic[k][m].vba;


                        fprintf(testfile,"potential zmax trouble spot, %d %s zmax %d iz %d flowvs %lf wh_at_zmax %lf total %lf\n",year[it],cname[i],zmax,iz,newdata[k][m][it].flowvs,ctdata[i][it].wh_at_zmax,total_avail);  

                      }
                      else {
                        printf("iz %d zmax %d\n",iz,zmax);
                      }

                    }  /* end of else iz < MAXZ */



		   }/* end of fnf*/

                }/* end of whr*/


                ctdata[i][it].whr-=newdata[k][m][it].vbh/1000.;
		

              }/*gcode*/
            }


          }  /* end of m */
        } /* end of k */

	    fclose(testfile);

        return;

}

/***********************************************************************/
void force_harvest2(int it, int i){

  FILE *testfile;
  int k, m, flag=0;
  double whr_orig=ZEROVALUE, flowvs_add=ZEROVALUE, flowss_add=ZEROVALUE, wh_from_sbh_only=ZEROVALUE;
  double tester1=ZEROVALUE, tester2=ZEROVALUE, tester3=ZEROVALUE, tester4=ZEROVALUE, tester5=ZEROVALUE;
  char outstat[2];

  
   
 
 
/* hold original amount of whr for spreading proportionally */

 whr_orig=ctdata[i][it].whr;



/* compute country level smb, smb_nf, vnfb, s area for force harvest */


 for (k=0;k<NY;k++){
   for (m=0;m<NX;m++){
     
     
#if COUNTRY_TEST
     if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){
       
#else
       if(dstatic[k][m].gcode > 0){
#endif
	 
	 if(dstatic[k][m].gcode == cdata[i].ccode){
	   
	 	     
	   if(dstatic[k][m].fnf == 0){

	     /* country level virgin nonforested biomass */
	     ctdata[i][it].vnfb+=(newdata[k][m][it].v-newdata[k][m][it].flowvc-newdata[k][m][it].flowvp-newdata[k][m][it].flowvu)*garea[k][m]*dstatic[k][m].vba/1000.;
	    
	     /* country level secondary nonforested biomass */
             ctdata[i][it].smb_nf+=(newdata[k][m][it].s-newdata[k][m][it].flowsc-newdata[k][m][it].flowsp-newdata[k][m][it].flowsu)*garea[k][m]*newdata[k][m][it].smb/1000.;

	     /* country level secondary nonforested area */
	       ctdata[i][it].sarea_nf+=(newdata[k][m][it].s+newdata[k][m][it].flowps+newdata[k][m][it].flowcs+newdata[k][m][it].flowvs+newdata[k][m][it].flowus-newdata[k][m][it].flowsc-newdata[k][m][it].flowsp-newdata[k][m][it].flowsu)*garea[k][m]; 
	     
 
	   }
           else {

	     /* country level secondary forest biomass from sbh only */
	     wh_from_sbh_only+=newdata[k][m][it].sbh/1000.;

	     /* country level secondary forest biomass */
	     ctdata[i][it].smb+=(newdata[k][m][it].s-newdata[k][m][it].flowsc-newdata[k][m][it].flowsp-newdata[k][m][it].flowsu)*((garea[k][m]*newdata[k][m][it].smb-newdata[k][m][it].sbh)/1000.); 

	     /* country level secondary forested area */
	        ctdata[i][it].sarea+=(newdata[k][m][it].s+newdata[k][m][it].flowps+newdata[k][m][it].flowcs+newdata[k][m][it].flowvs+newdata[k][m][it].flowus-newdata[k][m][it].flowsc-newdata[k][m][it].flowsp-newdata[k][m][it].flowsu)*garea[k][m]; 
	

	   }
	
   
	 }  /* end of gcode */
       }
       
     } /* end of m */
   } /* end of k */
   




   /* 4 cases for getting remaining wood harvest (whr) */


   if(ctdata[i][it].smb >= ctdata[i][it].whr){

     flag=1;  /* case 1, taking needed younger secondary forest proportionally to satisfy whr */                                 

     for (k=0;k<NY;k++){
       for (m=0;m<NX;m++){
	 
	 
#if COUNTRY_TEST
	 
	 if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){
	   
#else
	   if(dstatic[k][m].gcode > 0){
#endif
	     
	     if(dstatic[k][m].gcode == cdata[i].ccode){
	       
	       
	       if(dstatic[k][m].fnf == 1){
		 
		 if(ctdata[i][it].whr > ZEROVALUE){
		   
		   
		   newdata[k][m][it].sbh2=whr_orig/ctdata[i][it].smb*(newdata[k][m][it].s-newdata[k][m][it].flowsc-newdata[k][m][it].flowsp-newdata[k][m][it].flowsu)*(garea[k][m]*newdata[k][m][it].smb-newdata[k][m][it].sbh);
		   
		   ctdata[i][it].whr-=newdata[k][m][it].sbh2/1000.;
		   
		   
		 } /*end of whr */
	       } /*end of fnf */
	     }  /* end of gcode */
	   }
	   
	 } /* end of m */
       } /* end of k */

       
     }/* end case 1*/
     

     
     else if((cdata[i].smart_flow_option == 1) && ((ctdata[i][it].smb+ctdata[i][it].vnfb) >= ctdata[i][it].whr)){    

       flag=2;  /* case 2, primary priority, not enough young secondary forest to satisfy demand, 
		   therefore taking all young secondary forest and needed virgin nonforest proportionally 
		   to satisfy whr */

  

       for (k=0;k<NY;k++){
	 for (m=0;m<NX;m++){
	   
	   
#if COUNTRY_TEST
	   
	   if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){
	     
#else
	     if(dstatic[k][m].gcode > 0){
#endif
	       
	       if(dstatic[k][m].gcode == cdata[i].ccode){
		 
		 
		 if(ctdata[i][it].whr > ZEROVALUE){
		   
		   if(dstatic[k][m].fnf == 1) {

		     newdata[k][m][it].sbh2=(newdata[k][m][it].smb*garea[k][m]-newdata[k][m][it].sbh)*(newdata[k][m][it].s-newdata[k][m][it].flowsc-newdata[k][m][it].flowsp-newdata[k][m][it].flowsu);

		     ctdata[i][it].whr-=newdata[k][m][it].sbh2/1000.;

		   } /* end of fnf*/		   
		 } /*end of whr*/
	       }  /* end of gcode */
	     }
	     
	   } /* end of m */
	 } /* end of k */



	 /* take the whr proportionally out of vnfb; note needed to consider how much was taken
            out of secondary forest (ctdata[i][it].smb) */


	 for (k=0;k<NY;k++){
	   for (m=0;m<NX;m++){
	     
	     
#if COUNTRY_TEST
	     
	   if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){
	     
#else
	     if(dstatic[k][m].gcode > 0){
#endif
	       
	       if(dstatic[k][m].gcode == cdata[i].ccode){
		 
		 
		 if(ctdata[i][it].whr > ZEROVALUE){
		   
		   if(dstatic[k][m].fnf == 0) {
		     
		    
                        flowvs_add=(newdata[k][m][it].v-newdata[k][m][it].flowvc-newdata[k][m][it].flowvp-newdata[k][m][it].flowvu)*(whr_orig-ctdata[i][it].smb)/ctdata[i][it].vnfb;

		     if(flowvs_add <= -ZEROVALUE) {flowvs_add = ZEROVALUE;}
		     
		     newdata[k][m][it].flowvs+=flowvs_add;
		     
		     newdata[k][m][it].vbh2=flowvs_add*dstatic[k][m].vba*garea[k][m];
		     
		     ctdata[i][it].whr-=newdata[k][m][it].vbh2/1000.;
		     
		   } /*end of fnf */
		   
		  

		   
		 } /*end of whr*/
		 
	       }  /* end of gcode */
	     }
	     
	   } /* end of m */
	 } /* end of k */
	 
	 
       } /*end case 2 */
	 

	 
       else if((cdata[i].smart_flow_option == 2) && ((ctdata[i][it].smb+ctdata[i][it].smb_nf) >= ctdata[i][it].whr)){  

       flag=3;  /* case 3, secondary priority, taking all young secondary forest and needed secondary 
                   nonforest proportionally to satisfy whr */


     
       for (k=0;k<NY;k++){
	 for (m=0;m<NX;m++){
	   
	   
#if COUNTRY_TEST
	   
	   if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){
	     
#else
	     if(dstatic[k][m].gcode > 0){
#endif
	       
	       if(dstatic[k][m].gcode == cdata[i].ccode){
		 
		 
		 if(ctdata[i][it].whr > ZEROVALUE){
		   
		   if(dstatic[k][m].fnf == 1) {

		     newdata[k][m][it].sbh2=(newdata[k][m][it].smb*garea[k][m]-newdata[k][m][it].sbh)*(newdata[k][m][it].s-newdata[k][m][it].flowsc-newdata[k][m][it].flowsp-newdata[k][m][it].flowsu);

		     ctdata[i][it].whr-=newdata[k][m][it].sbh2/1000.;


		   } /* end of fnf*/		   
		 } /*end of whr*/
	       }  /* end of gcode */
	     }
	     
	   } /* end of m */
	 } /* end of k */



	 /* take the whr proportionally out of secondary nonforest; note needed to consider how much was taken
            out of secondary forest (ctdata[i][it].smb) */


	 for (k=0;k<NY;k++){
	   for (m=0;m<NX;m++){
	     
	     
#if COUNTRY_TEST
	     
	   if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){
	     
#else
	     if(dstatic[k][m].gcode > 0){
#endif
	       
	       if(dstatic[k][m].gcode == cdata[i].ccode){
		 
		 
		 if(ctdata[i][it].whr > ZEROVALUE){
		   
		   if(dstatic[k][m].fnf == 0) {
		     
		    
                        flowss_add=(newdata[k][m][it].s-newdata[k][m][it].flowsc-newdata[k][m][it].flowsp-newdata[k][m][it].flowsu)*(whr_orig-ctdata[i][it].smb)/ctdata[i][it].smb_nf;
		 
		     if(flowss_add <= -ZEROVALUE) {flowss_add = ZEROVALUE;}

		     newdata[k][m][it].sbh3=flowss_add*newdata[k][m][it].smb*garea[k][m];
		     
		     ctdata[i][it].whr-=newdata[k][m][it].sbh3/1000.;


		   } /*end of fnf */
		   
		  

		   
		 } /*end of whr*/
		 
	       }  /* end of gcode */
	     }
	     
	   } /* end of m */
	 } /* end of k */
	 
	 
       } /*end case 3 */


	 
       else{
	 
  
         flag=4;   /* case 4, taking all young secondary forest, then depending on priority;
                      primary priority take all vnfb and any additional amount needed from smb_nf
                      secondary priority, take all smb_nf and any additional amount needed from vnfb */


	 for (k=0;k<NY;k++){
	   for (m=0;m<NX;m++){

	     
#if COUNTRY_TEST
	     
	     if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){
	       
#else
	       if(dstatic[k][m].gcode > 0){
#endif
		 
		 if(dstatic[k][m].gcode == cdata[i].ccode){
		   
		   if(dstatic[k][m].fnf == 1){

		     /* take all young secondary forest */
		       
		     newdata[k][m][it].sbh2=(newdata[k][m][it].smb*garea[k][m]-newdata[k][m][it].sbh)*(newdata[k][m][it].s-newdata[k][m][it].flowsc-newdata[k][m][it].flowsp-newdata[k][m][it].flowsu);

		     ctdata[i][it].whr-=newdata[k][m][it].sbh2/1000.;

		   }
		   else{

		     if(cdata[i].smart_flow_option == 1){

		       /* primary priority, first take all vnfb */


		       flowvs_add=newdata[k][m][it].v-newdata[k][m][it].flowvc-newdata[k][m][it].flowvp-newdata[k][m][it].flowvu;
		       
		       if(flowvs_add <= -ZEROVALUE) {flowvs_add = ZEROVALUE;}
		       
		       newdata[k][m][it].flowvs+=flowvs_add;
		     
		       newdata[k][m][it].vbh2=flowvs_add*dstatic[k][m].vba*garea[k][m];
		     
		       ctdata[i][it].whr-=newdata[k][m][it].vbh2/1000.;

  
		       /* primary priority, not enough smb and vnfb, determined how much can be taken from smb_nf */

		       /* first try taking it proportionally */

		      
                       if(ctdata[i][it].smb_nf >= (whr_orig-ctdata[i][it].smb-ctdata[i][it].vnfb)){
		    
			 
                            flowss_add=(newdata[k][m][it].s-newdata[k][m][it].flowsc-newdata[k][m][it].flowsp-newdata[k][m][it].flowsu)*(whr_orig-ctdata[i][it].smb-ctdata[i][it].vnfb)/ctdata[i][it].smb_nf;
			   
			 if(flowss_add <= -ZEROVALUE) {flowss_add = ZEROVALUE;}
			   
			 newdata[k][m][it].sbh3=flowss_add*newdata[k][m][it].smb*garea[k][m];
			   
			 ctdata[i][it].whr-=newdata[k][m][it].sbh3/1000.;	

			
			 			 
		       }
		       else{ /* not enough in smb_nf, therefore take all remaining smb_nf, the rest is unmet */

			 flowss_add=newdata[k][m][it].s-newdata[k][m][it].flowsc-newdata[k][m][it].flowsp-newdata[k][m][it].flowsu;
			   
			 if(flowss_add <= -ZEROVALUE) {flowss_add = ZEROVALUE;}
			   
			 newdata[k][m][it].sbh3=flowss_add*newdata[k][m][it].smb*garea[k][m];

			 ctdata[i][it].whr-=newdata[k][m][it].sbh3/1000.;

		       }


		     } 
		     else { /* if smart_flow_option == 2 */

		       /* secondary priority, first take all smb_nf */

			 flowss_add=newdata[k][m][it].s-newdata[k][m][it].flowsc-newdata[k][m][it].flowsp-newdata[k][m][it].flowsu;
			   
			 if(flowss_add <= -ZEROVALUE) {flowss_add = ZEROVALUE;}
			   
			 newdata[k][m][it].sbh3=flowss_add*newdata[k][m][it].smb*garea[k][m];

			 ctdata[i][it].whr-=newdata[k][m][it].sbh3/1000.;


		       /* secondary priority, not enough smb and smb_nf, determined how much can be 
                          taken from vnfb */

		       /* first try taking it proportionally */

			 
			 if(ctdata[i][it].vnfb >= (whr_orig-ctdata[i][it].smb-ctdata[i][it].smb_nf)){

			 

flowvs_add=(newdata[k][m][it].v-newdata[k][m][it].flowvc-newdata[k][m][it].flowvp-newdata[k][m][it].flowvu)*(whr_orig-ctdata[i][it].smb-ctdata[i][it].smb_nf)/ctdata[i][it].vnfb;

			 if(flowvs_add <= -ZEROVALUE) {flowvs_add = ZEROVALUE;}

			 newdata[k][m][it].flowvs+=flowvs_add;
		     
			 newdata[k][m][it].vbh2=flowvs_add*dstatic[k][m].vba*garea[k][m];

		     
			 ctdata[i][it].whr-=newdata[k][m][it].vbh2/1000.;
	   
		       }
		       else{ /* not enough in vnfb, therefore take all remaining vnfb, the rest is unmet */


			 flowvs_add=newdata[k][m][it].v-newdata[k][m][it].flowvc-newdata[k][m][it].flowvp-newdata[k][m][it].flowvu;
		       
			 if(flowvs_add <= -ZEROVALUE) {flowvs_add = ZEROVALUE;}
		       
			 newdata[k][m][it].flowvs+=flowvs_add;
		     
			 newdata[k][m][it].vbh2=flowvs_add*dstatic[k][m].vba*garea[k][m];
		     
			 ctdata[i][it].whr-=newdata[k][m][it].vbh2/1000.;

  		       }


		     } /* ends else of smart_flow if */
		   
		   } /* ends else of fnf */
		 

		 }  /* end of gcode */
	       }
	       
	     } /* end of m */
	   } /* end of k */
	   
	  	   
	 } /* end case 4*/
	 

	 for (k=0;k<NY;k++){
	   for (m=0;m<NX;m++){
	     
#if COUNTRY_TEST
	     
	     if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){
	       
#else
	       if(dstatic[k][m].gcode > 0){
#endif
		 
		 if(dstatic[k][m].gcode == cdata[i].ccode){


		   tester1+=newdata[k][m][it].vbh/1000.;
		   tester2+=newdata[k][m][it].vbh2/1000.;
		   tester3+=newdata[k][m][it].sbh/1000.;
		   tester4+=newdata[k][m][it].sbh2/1000.;
		   tester5+=newdata[k][m][it].sbh3/1000.;

		   
		 }
	       }
	     }
	   }



	   if((strcmp("tone",trun) == 0) && (it==0)) {
	     strcpy(outstat,"w");
	   }
	   else {
	     strcpy(outstat,"a");
	   }
	   testfile=fopen("wh.latest",outstat);
	   
       
	   fprintf(testfile,"%d %s flag %d fao# %lf wh_needed %lf vbh %lf vbh2 %lf sbh %lf sbh2 %lf sbh3 %lf whtot %lf smb %lf smb_nf %lf vnfb %lf\n",year[it],cname[i],flag,ctdata[i][it].wh,whr_orig,tester1,tester2,tester3,tester4,tester5,(tester1+tester2+tester3+tester4+tester5),ctdata[i][it].smb,ctdata[i][it].smb_nf,ctdata[i][it].vnfb);
	   
	   
	   fclose(testfile);
	   

	 return;
	 
       }

/***********************************************************************/

float prob_harv(float biomass){
  
  /* test function for development */
  
  int i, flagger=0;
  float p=ZEROVALUE;
  



  for (i=0;i<NPHB-1;i++){

    if(biomass > (i*0.5) && biomass <= ((i+1)*0.5)) {

       p=phb[i];
       flagger=1;
    }
    if(i==0 && biomass <= (i*1.0)) flagger=1;

  }
  if(flagger == 0) printf("BUG: prob_harv .not. set, b: %f\n",biomass);


  return(p);
  
}
/***********************************************************************/
void update_vb(int it){
  
  int iz, k, m, i;
  
  
  /* determine virgin biomass at the country level (vb); vb is the virgin
     available for a particular country for particular iz; therefore given
     a country code and a value of iz, vb is known */ 
  
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){


#if COUNTRY_TEST
                      
            if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){
                      
#else
              if(dstatic[k][m].gcode > 0){
#endif


  
        i=dstatic[k][m].newgcode;
	
	  
	
	  for (iz=0;iz<NZ;iz++){
	    
	    if(newdata[k][m][it].zdis == iz){



	      if(dstatic[k][m].fnf == 1){
 
		
		/* vb in units = MgC */
		cztdata[i][iz][it].vb+=
		  (newdata[k][m][it].v-newdata[k][m][it].flowvc-newdata[k][m][it].flowvp-newdata[k][m][it].flowvu)*garea[k][m]*dstatic[k][m].vba/1000.0;

              
		
	      }
	    }
	  }

	}
	  
	
      } /* end of m */
    } /* end of k */
    
   


  return;
}

/********************************************************************/
void zdis_calc(int it){
  
  FILE *testfile;
  int i;
  int iz, k, m, zcheck, zvalue, flagger; 
  

  
  /* determine focal cell for zdis, zdis is initialized to NZ earlier for all times*/ 
   
  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){
       if(dstatic[k][m].gcode > 0){
      i=dstatic[k][m].newgcode;

      
      if(cdata[i].zdis_option == 1){

  	/* setting zdis=0 for land use cells and masking out ocean */
      
      
	 if(newdata[k][m][it].c > ZEROVALUE || newdata[k][m][it].p > ZEROVALUE ||
	 newdata[k][m][it].s > ZEROVALUE || newdata[k][m][it].u > ZEROVALUE) {
      

	
        newdata[k][m][it].zdis = 0;
	
      }  
      
      }

      else{
       newdata[k][m][it].zdis = 0;
      }
       }
    }
  }
  
  
  /* check all other cells relative to focal cell */
  
  
  for (iz=0;iz<NZ;iz++){
    
    zcheck=iz;
    zvalue=iz+1;
    
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){
	
#if 0
#if COUNTRY_TEST
	
	if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){   
#else
	  if(dstatic[k][m].gcode > 0){
#endif	
#endif	    
	    if((k-1 > -1) && (m-1 > -1)) {
	      if(newdata[k][m][it].zdis == zcheck) {
		if(newdata[k-1][m-1][it].zdis > zcheck) {
 		  newdata[k-1][m-1][it].zdis = zvalue;
		}
	      }
	    }
	    
	    if((k-1 > -1) && (m == m)) {
	      if(newdata[k][m][it].zdis == zcheck) {
		if(newdata[k-1][m][it].zdis > zcheck) {
                 
		  newdata[k-1][m][it].zdis = zvalue;
		}
	      }
	    }
	
	    if((k-1 > -1) && (m+1 < NX)) {
	      if(newdata[k][m][it].zdis == zcheck) {
		if(newdata[k-1][m+1][it].zdis > zcheck) {
		  newdata[k-1][m+1][it].zdis = zvalue;
		}
	      }
	    }
	    
	    if((k == k) && (m-1 > -1)) {
	      if(newdata[k][m][it].zdis == zcheck) {
		if(newdata[k][m-1][it].zdis > zcheck) {
		  newdata[k][m-1][it].zdis = zvalue;
		}
	      }
	    }
	    
	    if((k == k) && (m+1 < NX)) {
	      if(newdata[k][m][it].zdis == zcheck) {
		if(newdata[k][m+1][it].zdis > zcheck) {
		  newdata[k][m+1][it].zdis = zvalue;
		}
	      }
	    }
	    
	    if((k+1 < NY) && (m-1 > -1)) {
	      if(newdata[k][m][it].zdis == zcheck) {
		if(newdata[k+1][m-1][it].zdis > zcheck) {
		  newdata[k+1][m-1][it].zdis = zvalue;
		}
	      }
	    }
	    
	    if((k+1 < NY) && (m == m)) {
	      if(newdata[k][m][it].zdis == zcheck) {
		if(newdata[k+1][m][it].zdis > zcheck) {
		  newdata[k+1][m][it].zdis = zvalue;
		}
	      }
	    }
	    
	    if((k+1 < NY) && (m+1 < NX)){ 
	      if(newdata[k][m][it].zdis == zcheck) {
		if(newdata[k+1][m+1][it].zdis > zcheck) {
		  newdata[k+1][m+1][it].zdis = zvalue;
		}
	      }      
	    }
	    
	   

	}   /* end of m */
      }    /* end of k */



    }     /* end of iz */
    


	
  return;
}

/***********************************************************************/
void country_cpavg(int it, int i){
	  
  int k, m, im, cellcounter=0;
  double garea_sum=ZEROVALUE;


  
  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){
      
      
#if COUNTRY_TEST      
	      if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){
#else
		if(dstatic[k][m].gcode > 0){
#endif
		  
		  
		  if(dstatic[k][m].gcode == cdata[i].ccode) {
		    
		    ctdata[i][it].cavg+=newdata[k][m][it].c*garea[k][m]; 
                    ctdata[i][it].pavg+=newdata[k][m][it].p*garea[k][m];

		    ctdata[i][it].cpavg+=(newdata[k][m][it].c+newdata[k][m][it].p)*garea[k][m];

		    if((newdata[k][m][it].c > ZEROVALUE) || (newdata[k][m][it].p > ZEROVALUE)) ctdata[i][it].cpgarea+=garea[k][m];


                    if(newdata[k][m][it].c > ZEROVALUE) ctdata[i][it].cgarea+=garea[k][m];
                    if(newdata[k][m][it].p > ZEROVALUE) ctdata[i][it].pgarea+=garea[k][m];
		    if(newdata[k][m][it].p > ZEROVALUE) ctdata[i][it].pcount++;



		  }/*gcode*/
		}
		  
		  
	      }  /* end of m */
    } /* end of k */

    if(cdata[i].ccode > 0) {
      if(ctdata[i][it].cgarea > ZEROVALUE) ctdata[i][it].cavg=ctdata[i][it].cavg/ctdata[i][it].cgarea;
      if(ctdata[i][it].pgarea > ZEROVALUE) ctdata[i][it].pavg=ctdata[i][it].pavg/ctdata[i][it].pgarea;
      if(ctdata[i][it].cpgarea > ZEROVALUE) ctdata[i][it].cpavg=ctdata[i][it].cpavg/ctdata[i][it].cpgarea;
    }


      
    return;
	      
  }


/********************************************************************/
void update_states(int it, int zmax){

  int i,k,m, iz;
  float sma_area_gained, sma_area_lost, sma_area_notlost;
  float sum_test;
  double gsum=ZEROVALUE, tester=ZEROVALUE;

  check_total[it]=0;


  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){

#if COUNTRY_TEST
	  
     if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){    
		
#else
	    if(dstatic[k][m].gcode > 0){
#endif	
      

#if FLOW_BUG_PRINT



    	/* checking for negative zero flows */

	if(newdata[k][m][it].flowcp < ZEROVALUE) {printf("Update states flow Bug: flowcp<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowcp);newdata[k][m][it].flowcp = ZEROVALUE;}
	if(newdata[k][m][it].flowpc < ZEROVALUE) {printf("Update states flow Bug: flowpc<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowpc);newdata[k][m][it].flowpc = ZEROVALUE;}
	if(newdata[k][m][it].flowpv < ZEROVALUE) {printf("Update states flow Bug: flowpv<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowpv);newdata[k][m][it].flowpv = ZEROVALUE;}
	if(newdata[k][m][it].flowvp < ZEROVALUE) {printf("Update states flow Bug: flowvp<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowvp);newdata[k][m][it].flowvp = ZEROVALUE;}
	if(newdata[k][m][it].flowvc < ZEROVALUE) {printf("Update states flow Bug: flowvc<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowvc);newdata[k][m][it].flowvc = ZEROVALUE;}
	if(newdata[k][m][it].flowcv < ZEROVALUE) {printf("Update states flow Bug: flowcv<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowcv);newdata[k][m][it].flowcv = ZEROVALUE;}
	if(newdata[k][m][it].flowsp < ZEROVALUE) {printf("Update states flow Bug: flowsp<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowsp);newdata[k][m][it].flowsp = ZEROVALUE;}
	if(newdata[k][m][it].flowps < ZEROVALUE) {printf("Update states flow Bug: flowps<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowps);newdata[k][m][it].flowps = ZEROVALUE;}
	if(newdata[k][m][it].flowsc < ZEROVALUE) {printf("Update states flow Bug: flowsc<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowsc);newdata[k][m][it].flowsc = ZEROVALUE;}
	if(newdata[k][m][it].flowcs < ZEROVALUE) {printf("Update states flow Bug: flowcs<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowcs);newdata[k][m][it].flowcs = ZEROVALUE;}
	if(newdata[k][m][it].flowvs < ZEROVALUE) {printf("Update states flow Bug: flowvs<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowvs);newdata[k][m][it].flowvs = ZEROVALUE;}
	if(newdata[k][m][it].flowcu < ZEROVALUE) {printf("Update states flow Bug: flowcu<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowcu);newdata[k][m][it].flowcu = ZEROVALUE;}
	if(newdata[k][m][it].flowpu < ZEROVALUE) {printf("Update states flow Bug: flowpu<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowpu);newdata[k][m][it].flowpu = ZEROVALUE;}
	if(newdata[k][m][it].flowvu < ZEROVALUE) {printf("Update states flow Bug: flowvu<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowvu);newdata[k][m][it].flowvu = ZEROVALUE;}
	if(newdata[k][m][it].flowsu < ZEROVALUE) {printf("Update states flow Bug: flowsu<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowsu);newdata[k][m][it].flowsu = ZEROVALUE;}
	if(newdata[k][m][it].flowuc < ZEROVALUE) {printf("Update states flow Bug: flowuc<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowuc);newdata[k][m][it].flowuc = ZEROVALUE;}
	if(newdata[k][m][it].flowup < ZEROVALUE) {printf("Update states flow Bug: flowup<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowup);newdata[k][m][it].flowup = ZEROVALUE;}
	if(newdata[k][m][it].flowus < ZEROVALUE) {printf("Update states flow Bug: flowus<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowus);newdata[k][m][it].flowus = ZEROVALUE;}
	


	if(newdata[k][m][it].vbh < ZEROVALUE) {printf("Harvest Bug: vbh<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].vbh);newdata[k][m][it].vbh = ZEROVALUE;}
	if(newdata[k][m][it].sbh < ZEROVALUE) {printf("Harvest Bug: sbh<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].sbh);newdata[k][m][it].sbh = ZEROVALUE;}
	if(newdata[k][m][it].vbh2 < ZEROVALUE) {printf("Harvest Bug: vbh2<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].vbh2);newdata[k][m][it].vbh2 = ZEROVALUE;}
	if(newdata[k][m][it].sbh2 < ZEROVALUE) {printf("Harvest Bug: sbh2<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].sbh2);newdata[k][m][it].sbh2 = ZEROVALUE;}


#else	

	if(newdata[k][m][it].flowcp < ZEROVALUE) {newdata[k][m][it].flowcp = ZEROVALUE;}
	if(newdata[k][m][it].flowpc < ZEROVALUE) {newdata[k][m][it].flowpc = ZEROVALUE;}
	if(newdata[k][m][it].flowpv < ZEROVALUE) {newdata[k][m][it].flowpv = ZEROVALUE;}
	if(newdata[k][m][it].flowvp < ZEROVALUE) {newdata[k][m][it].flowvp = ZEROVALUE;}
	if(newdata[k][m][it].flowvc < ZEROVALUE) {newdata[k][m][it].flowvc = ZEROVALUE;}
	if(newdata[k][m][it].flowcv < ZEROVALUE) {newdata[k][m][it].flowcv = ZEROVALUE;}
	if(newdata[k][m][it].flowsp < ZEROVALUE) {newdata[k][m][it].flowsp = ZEROVALUE;}
	if(newdata[k][m][it].flowps < ZEROVALUE) {newdata[k][m][it].flowps = ZEROVALUE;}
	if(newdata[k][m][it].flowsc < ZEROVALUE) {newdata[k][m][it].flowsc = ZEROVALUE;}
	if(newdata[k][m][it].flowcs < ZEROVALUE) {newdata[k][m][it].flowcs = ZEROVALUE;}
	if(newdata[k][m][it].flowvs < ZEROVALUE) {newdata[k][m][it].flowvs = ZEROVALUE;}

	if(newdata[k][m][it].flowcu < ZEROVALUE) {newdata[k][m][it].flowcu = ZEROVALUE;}
	if(newdata[k][m][it].flowpu < ZEROVALUE) {newdata[k][m][it].flowpu = ZEROVALUE;}
	if(newdata[k][m][it].flowvu < ZEROVALUE) {newdata[k][m][it].flowvu = ZEROVALUE;}
	if(newdata[k][m][it].flowsu < ZEROVALUE) {newdata[k][m][it].flowsu = ZEROVALUE;}
	if(newdata[k][m][it].flowuc < ZEROVALUE) {newdata[k][m][it].flowuc = ZEROVALUE;}
	if(newdata[k][m][it].flowup < ZEROVALUE) {newdata[k][m][it].flowup = ZEROVALUE;}
	if(newdata[k][m][it].flowus < ZEROVALUE) {newdata[k][m][it].flowus = ZEROVALUE;}

#endif


	/*updating states*/
	

	newdata[k][m][it+1].v=newdata[k][m][it].v+newdata[k][m][it].flowpv+
	  newdata[k][m][it].flowcv-newdata[k][m][it].flowvp-
	  newdata[k][m][it].flowvc-newdata[k][m][it].flowvs-newdata[k][m][it].flowvu;
	
	newdata[k][m][it+1].s=newdata[k][m][it].s+newdata[k][m][it].flowps+
	  newdata[k][m][it].flowvs+newdata[k][m][it].flowcs-
	  newdata[k][m][it].flowsp-newdata[k][m][it].flowsc-newdata[k][m][it].flowsu+newdata[k][m][it].flowus;
	   

	newdata[k][m][it+1].c=newdata[k][m][it].c +
	  newdata[k][m][it].flowpc + newdata[k][m][it].flowvc +
	  newdata[k][m][it].flowsc - newdata[k][m][it].flowcp -
	  newdata[k][m][it].flowcv - newdata[k][m][it].flowcs -
          newdata[k][m][it].flowcu + newdata[k][m][it].flowuc;
	
	newdata[k][m][it+1].p=newdata[k][m][it].p +
	  newdata[k][m][it].flowcp + newdata[k][m][it].flowvp +
	  newdata[k][m][it].flowsp - newdata[k][m][it].flowpc -
	  newdata[k][m][it].flowpv - newdata[k][m][it].flowps -
          newdata[k][m][it].flowpu + newdata[k][m][it].flowup;

	newdata[k][m][it+1].u=newdata[k][m][it].u +
	  newdata[k][m][it].flowcu + newdata[k][m][it].flowvu +
	  newdata[k][m][it].flowsu + newdata[k][m][it].flowpu - 
          newdata[k][m][it].flowuc - newdata[k][m][it].flowus - 
          newdata[k][m][it].flowup;


#if STATE_BUG_PRINT
      
	if(newdata[k][m][it+1].v < ZEROVALUE) { 
           printf("Update states Bug: v<0 at t+1, resetting to 0 by force. k=%d m=%d it=%d it+1=%d v(t)=%lf v(t+1)=%lf pv %lf cv %lf vp %lf vc %lf vs %lf\n",k,m,it,it+1,newdata[k][m][it].v,newdata[k][m][it+1].v,newdata[k][m][it].flowpv,newdata[k][m][it].flowcv,newdata[k][m][it].flowvp,newdata[k][m][it].flowvc,newdata[k][m][it].flowvs);                
           newdata[k][m][it+1].v = ZEROVALUE;
        }
	if(newdata[k][m][it+1].s < ZEROVALUE) {         
           printf("Update states Bug: s<0 at t+1, resetting to 0 by force. k=%d m=%d it=%d it+1=%d s(t)=%lf s(t+1)=%lf cs %lf ps %lf vs %lf sc %lf sp %lf\n",k,m,it,it+1,newdata[k][m][it].s,newdata[k][m][it+1].s,newdata[k][m][it].flowcs,newdata[k][m][it].flowps,newdata[k][m][it].flowvs,newdata[k][m][it].flowsc,newdata[k][m][it].flowsp);
           newdata[k][m][it+1].s =ZEROVALUE;  
        }
	if(newdata[k][m][it+1].c < ZEROVALUE) {   
           printf("Update states Bug: c<0 at t+1, resetting to 0 by force. k=%d m=%d it=%d it+1=%d c(t)=%lf c(t+1)=%lf pc %lf vc %lf sc %lf cp %lf cv %lf cs %lf\n",k,m,it,it+1,newdata[k][m][it].c,newdata[k][m][it+1].c,newdata[k][m][it].flowpc,newdata[k][m][it].flowvc,newdata[k][m][it].flowsc,newdata[k][m][it].flowcp,newdata[k][m][it].flowcv,newdata[k][m][it].flowcs);
           newdata[k][m][it+1].c = ZEROVALUE; 
        }
	if(newdata[k][m][it+1].p < ZEROVALUE) {
           printf("Update states Bug: p<0 at t+1, resetting to 0 by force. k=%d m=%d it=%d it+1=%d p(t)=%lf p(t+1)=%lf cp %lf vp %lf sp %lf pc %lf pv %lf ps %lf\n",k,m,it,it+1,newdata[k][m][it].p,newdata[k][m][it+1].p,newdata[k][m][it].flowcp,newdata[k][m][it].flowvp,newdata[k][m][it].flowsp,newdata[k][m][it].flowpc,newdata[k][m][it].flowpv,newdata[k][m][it].flowps);
           newdata[k][m][it+1].p = ZEROVALUE; 
        }
	if(newdata[k][m][it+1].u < ZEROVALUE) {   
           printf("Update states Bug: u<0 at t+1, resetting to 0 by force. k=%d m=%d it=%d it+1=%d u(t)=%lf u(t+1)=%lf pu %lf vu %lf su %lf up %lf cu %lf us %lf\n",k,m,it,it+1,newdata[k][m][it].u,newdata[k][m][it+1].u,newdata[k][m][it].flowpu,newdata[k][m][it].flowvu,newdata[k][m][it].flowsu,newdata[k][m][it].flowup,newdata[k][m][it].flowcu,newdata[k][m][it].flowus);
           newdata[k][m][it+1].u = ZEROVALUE; 
        }

#else

	if(newdata[k][m][it+1].v < ZEROVALUE) newdata[k][m][it+1].v = ZEROVALUE;  
	if(newdata[k][m][it+1].s < ZEROVALUE) newdata[k][m][it+1].s = ZEROVALUE;
	if(newdata[k][m][it+1].c < ZEROVALUE) newdata[k][m][it+1].c = ZEROVALUE;
	if(newdata[k][m][it+1].p < ZEROVALUE) newdata[k][m][it+1].p = ZEROVALUE; 
	if(newdata[k][m][it+1].u < ZEROVALUE) newdata[k][m][it+1].u = ZEROVALUE; 
	
#endif



	sum_test=ZEROVALUE;
	
	sum_test=newdata[k][m][it+1].c+newdata[k][m][it+1].p+
                 newdata[k][m][it+1].v+newdata[k][m][it+1].s+
                 newdata[k][m][it+1].i+newdata[k][m][it+1].w+newdata[k][m][it+1].u;
	
	
	
        if(newdata[k][m][it].smb > ZEROVALUE_CHECK){

	  sma_area_gained=newdata[k][m][it].flowvs+newdata[k][m][it].flowcs+newdata[k][m][it].flowps+newdata[k][m][it].flowus+(newdata[k][m][it].sbh+newdata[k][m][it].sbh2+newdata[k][m][it].sbh3)/newdata[k][m][it].smb/garea[k][m];
	  
	  sma_area_lost=newdata[k][m][it].flowsp+newdata[k][m][it].flowsc+newdata[k][m][it].flowsu+
	    (newdata[k][m][it].sbh+newdata[k][m][it].sbh2+newdata[k][m][it].sbh3)/newdata[k][m][it].smb/garea[k][m];

	} 
        else {
	 
	  sma_area_gained=newdata[k][m][it].flowvs+newdata[k][m][it].flowcs+newdata[k][m][it].flowps+newdata[k][m][it].flowus;
	  sma_area_lost=newdata[k][m][it].flowsp+newdata[k][m][it].flowsc+newdata[k][m][it].flowsu;
	}


	sma_area_notlost=newdata[k][m][it].s-sma_area_lost;
	

 
	if((newdata[k][m][it].s+sma_area_gained-sma_area_lost) > ZEROVALUE_CHECK){
	  
	 


	  newdata[k][m][it+1].sma=(sma_area_notlost*(newdata[k][m][it].sma+1.0)+
	               sma_area_gained*1.0)/ 
                      (newdata[k][m][it].s+sma_area_gained-sma_area_lost);


	}
	else{
	 
	  newdata[k][m][it+1].sma=1.0;
	}



	
	if(dstatic[k][m].vba > ZEROVALUE){
	  
/* 0.38 is wood fraction of NPP, 0.75 is the aboveground fraction of NPP */

	

	 newdata[k][m][it+1].smb=dstatic[k][m].vba*(1.0000000-exp(-(dstatic[k][m].vnppa*0.75*0.38*newdata[k][m][it].sma)/dstatic[k][m].vba));   


	}
	else{
	  newdata[k][m][it+1].smb=ZEROVALUE;
	}
	
	
	if(newdata[k][m][it+1].sma < ZEROVALUE) {
	   printf("Bug: sma < 0, being reset to zero by force. k=%d m=%d it+1=%d sma(t) %10.15lf sma(t+1) %10.15lf sma_area_notlost %lf s(t) %lf sma_area_lost %lf sma_area_gained %lf smb(t) %lf\n",k,m,it+1,newdata[k][m][it].sma,newdata[k][m][it+1].sma,sma_area_notlost,newdata[k][m][it].s,sma_area_lost,sma_area_gained,newdata[k][m][it].smb);
	  
	   newdata[k][m][it+1].sma = 1.0;
	}
	
	
	   
	    }

	    if(dstatic[k][m].gcode>0){
	      check_total[it]+=(newdata[k][m][it].v+newdata[k][m][it].s+newdata[k][m][it].c+newdata[k][m][it].p+newdata[k][m][it].u+newdata[k][m][it].i+newdata[k][m][it].w)*garea[k][m];
	    }

         } /* end of m */
  } /* end of k */

   

  return;
}





/********************************************************************/
void output_lu(){

  FILE *outfile;
  int i, ih, k, m, it, itstart;
  char outstat[2];
  char latbase[4], lonbase[4];
  char alatshort[4], alatmedium[5], alatlong[6];
  char alonshorter[4], alonshort[5], alonmedium[6], alonlong[7];
  char fname[170], fout_trans[170];
  double flowsbh, flowvbh, flowsbh2, flowvbh2, flowsbh3;

  
  strcpy(lonbase,"lon");
  strcpy(latbase,"lat");
  
  
  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){
      
      strcpy(fname,"");
      strcpy(fout_trans,"");
      
      strcat(fname,PATH1);
      strcat(fname,finput_dir);
      strcat(fname,"lu/");
      strcat(fname,latbase);

      if(res_option==1) {
     
	if(lat[k] < -9.5 && lat[k] > -100.5) { 
	  sprintf(alatlong,"%5.1f",lat[k]);
	  strcat(fname,alatlong);
	}
	if(lat[k] > 9.5 && lat[k] < 100.5) {
	  sprintf(alatmedium,"%4.1f",lat[k]);
	  strcat(fname,alatmedium);
	}
	if(lat[k] > -10.5 && lat[k] < 0.) {
	  sprintf(alatmedium,"%4.1f",lat[k]);
	  strcat(fname,alatmedium);
	}
	if(lat[k] < 10.5 && lat[k] > 0.) {
	  sprintf(alatshort,"%3.1f",lat[k]);
	  strcat(fname,alatshort);
	} 
      
	strcat(fname,lonbase);
      
	if(lon[m] < -99.5) { 
	  sprintf(alonlong,"%6.1f",lon[m]);
	  strcat(fname,alonlong);
	}
	if(lon[m] > 99.5) {
	  sprintf(alonmedium,"%5.1f",lon[m]);
	  strcat(fname,alonmedium);
	}
	if(lon[m] > -100.5 && lon[m] < -9.5) {
	  sprintf(alonmedium,"%5.1f",lon[m]);
	  strcat(fname,alonmedium);
	}
	if(lon[m] < 100.5 && lon[m] > 9.5) {
	  sprintf(alonshort,"%4.1f",lon[m]);
	  strcat(fname,alonshort);
	}
	if(lon[m] > -10.5 && lon[m] < 0.) {
	  sprintf(alonshort,"%4.1f",lon[m]);
	  strcat(fname,alonshort);
	}
	if(lon[m] < 10.5 && lon[m] > 0.) {
	  sprintf(alonshorter,"%3.1f",lon[m]);
	  strcat(fname,alonshorter);
	}
      }
      if(res_option==2) {
	if(lat[k] < -9.75 && lat[k] > -100.25) { 
	  sprintf(alatlong,"%6.2f",lat[k]);
	  strcat(fname,alatlong);
	}
	if(lat[k] > 9.75 && lat[k] < 100.25) {
	  sprintf(alatmedium,"%5.2f",lat[k]);
	  strcat(fname,alatmedium);
	}
	if(lat[k] > -10.25 && lat[k] < 0.) {
	  sprintf(alatmedium,"%5.2f",lat[k]);
	  strcat(fname,alatmedium);
	}
	if(lat[k] < 10.25 && lat[k] > 0.) {
	  sprintf(alatshort,"%4.2f",lat[k]);
	  strcat(fname,alatshort);
	} 
      
	strcat(fname,lonbase);
      
	if(lon[m] < -99.75) { 
	  sprintf(alonlong,"%7.2f",lon[m]);
	  strcat(fname,alonlong);
	}
	if(lon[m] > 99.75) {
	  sprintf(alonmedium,"%6.2f",lon[m]);
	  strcat(fname,alonmedium);
	}
	if(lon[m] > -100.25 && lon[m] < -9.75) {
	  sprintf(alonmedium,"%6.2f",lon[m]);
	  strcat(fname,alonmedium);
	}
	if(lon[m] < 100.25 && lon[m] > 9.75) {
	  sprintf(alonshort,"%5.2f",lon[m]);
	  strcat(fname,alonshort);
	}
	if(lon[m] > -10.25 && lon[m] < 0.) {
	  sprintf(alonshort,"%5.2f",lon[m]);
	  strcat(fname,alonshort);
	}
	if(lon[m] < 10.25 && lon[m] > 0.) {
	  sprintf(alonshorter,"%4.2f",lon[m]);
	  strcat(fname,alonshorter);
	}
      }
      
      strcat(fname,".lu");
      
      for (i=0;i<sizeof(fname);i++) 
	if(fname[i] != ' ') fout_trans[i] = fname[i];  


#if COUNTRY_TEST
	  
	  if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){    	
#else
	    if(dstatic[k][m].gcode > 0){
#endif	
      
	
	if(strcmp("tone",trun) == 0) {
          strcpy(outstat,"w");
	}
	else {
          strcpy(outstat,"a");
	}

	outfile=fopen(fout_trans,outstat);
	
	if(strcmp("tone",trun) == 0) 
          fprintf(outfile,"year cp pc pv vp vc cv sc cs sp ps vs cu pu vu su uc up us sbh f_sbh vbh f_vbh sbh2 f_sbh2 vbh2 f_vbh2 sbh3 f_sbh3\n");
	
	
	for (it=0;it<NEWNT-1;it++){


        if(newdata[k][m][it].flowcp <= -ZEROVALUE) {newdata[k][m][it].flowcp = ZEROVALUE;}
        if(newdata[k][m][it].flowpc <= -ZEROVALUE) {newdata[k][m][it].flowpc = ZEROVALUE;}
        if(newdata[k][m][it].flowpv <= -ZEROVALUE) {newdata[k][m][it].flowpv = ZEROVALUE;}        
        if(newdata[k][m][it].flowvp <= -ZEROVALUE) {newdata[k][m][it].flowvp = ZEROVALUE;}        
        if(newdata[k][m][it].flowvc <= -ZEROVALUE) {newdata[k][m][it].flowvc = ZEROVALUE;}        
        if(newdata[k][m][it].flowcv <= -ZEROVALUE) {newdata[k][m][it].flowcv = ZEROVALUE;}        
        if(newdata[k][m][it].flowsp <= -ZEROVALUE) {newdata[k][m][it].flowsp = ZEROVALUE;}        
        if(newdata[k][m][it].flowps <= -ZEROVALUE) {newdata[k][m][it].flowps = ZEROVALUE;}        
        if(newdata[k][m][it].flowsc <= -ZEROVALUE) {newdata[k][m][it].flowsc = ZEROVALUE;}        
        if(newdata[k][m][it].flowcs <= -ZEROVALUE) {newdata[k][m][it].flowcs = ZEROVALUE;}
        if(newdata[k][m][it].flowvs <= -ZEROVALUE) {newdata[k][m][it].flowvs = ZEROVALUE;}
        if(newdata[k][m][it].flowcu <= -ZEROVALUE) {newdata[k][m][it].flowcu = ZEROVALUE;}
        if(newdata[k][m][it].flowpu <= -ZEROVALUE) {newdata[k][m][it].flowpu = ZEROVALUE;}
        if(newdata[k][m][it].flowvu <= -ZEROVALUE) {newdata[k][m][it].flowvu = ZEROVALUE;}
        if(newdata[k][m][it].flowsu <= -ZEROVALUE) {newdata[k][m][it].flowsu = ZEROVALUE;}
        if(newdata[k][m][it].flowuc <= -ZEROVALUE) {newdata[k][m][it].flowuc = ZEROVALUE;}
        if(newdata[k][m][it].flowup <= -ZEROVALUE) {newdata[k][m][it].flowup = ZEROVALUE;}
        if(newdata[k][m][it].flowus <= -ZEROVALUE) {newdata[k][m][it].flowus = ZEROVALUE;}

	flowsbh=flowsbh2=flowvbh=flowvbh2=flowsbh3=ZEROVALUE;
	
	if(newdata[k][m][it].smb > ZEROVALUE) flowsbh=newdata[k][m][it].sbh/newdata[k][m][it].smb/garea[k][m];
	if(dstatic[k][m].vba > ZEROVALUE) flowvbh=newdata[k][m][it].vbh/dstatic[k][m].vba/garea[k][m];
	if(newdata[k][m][it].smb > ZEROVALUE) flowsbh2=newdata[k][m][it].sbh2/newdata[k][m][it].smb/garea[k][m];
	if(dstatic[k][m].vba > ZEROVALUE) flowvbh2=newdata[k][m][it].vbh2/dstatic[k][m].vba/garea[k][m];
	if(newdata[k][m][it].smb > ZEROVALUE) flowsbh3=newdata[k][m][it].sbh3/newdata[k][m][it].smb/garea[k][m];



          fprintf(outfile,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                  year[it],
		  newdata[k][m][it].flowcp,
		  newdata[k][m][it].flowpc,
		  newdata[k][m][it].flowpv,
		  newdata[k][m][it].flowvp,
		  newdata[k][m][it].flowvc,
		  newdata[k][m][it].flowcv,
		  newdata[k][m][it].flowsc,
		  newdata[k][m][it].flowcs,
		  newdata[k][m][it].flowsp,
		  newdata[k][m][it].flowps,
                  newdata[k][m][it].flowvs,
                  newdata[k][m][it].flowcu,
                  newdata[k][m][it].flowpu,
                  newdata[k][m][it].flowvu,
                  newdata[k][m][it].flowsu,
                  newdata[k][m][it].flowuc,
                  newdata[k][m][it].flowup,
                  newdata[k][m][it].flowus,
		  newdata[k][m][it].sbh/1000,
		  flowsbh,
                  newdata[k][m][it].vbh/1000, 
		  flowvbh,
		  newdata[k][m][it].sbh2/1000,
		  flowsbh2,
                  newdata[k][m][it].vbh2/1000,
		  flowvbh2,
		  newdata[k][m][it].sbh3/1000,
		  flowsbh3);
	}
       fclose(outfile);
      }
      
    }
  }
  return;
}

/********************************************************************/
void output_updated_states(){
  
  FILE *outfile; 
  int ih, ic, k, m, it, cnum, nit;
  char newpath[190], newpath1[190], fn[30], ayear[10], fout[190];
  
  printf("output updated states \n");
  strcpy(newpath,PATH1);
  strcat(newpath,finput_dir);
  strcat(newpath,"updated_states/");


  cnum=8;
  
  
  for (ic=0;ic<cnum;ic++){

     for (it=0;it<NEWNT;it++){ 


      nit=it;
      strcpy(newpath1,newpath);
      sprintf(ayear,"%d.txt",year[nit]);

      if(ic == 0) strcpy(fn,"gothr.");
      if(ic == 1) strcpy(fn,"gcrop.");
      if(ic == 2) strcpy(fn,"gpast.");
      if(ic == 3) strcpy(fn,"gsecd.");
      if(ic == 4) strcpy(fn,"gssmb.");
      if(ic == 5) strcpy(fn,"gssma.");
      if(ic == 6) strcpy(fn,"gsumm.");
      if(ic == 7) strcpy(fn,"gurbn.");

      strcat(fn,ayear);
      strcat(newpath1,fn);
      strcpy(fout,newpath1); 

      outfile=fopen(fout,"w");

      for (ih=0;ih<HEADLEN;ih++) fprintf(outfile,"%s",header[ih]);
      
 
      for (k=0;k<NY;k++){
	for (m=0;m<NX;m++){
	  

	    if(ic == 0) fprintf(outfile,"%lf ",newdata[k][m][nit].v);
	    if(ic == 1) fprintf(outfile,"%lf ",newdata[k][m][nit].c);
	    if(ic == 2) fprintf(outfile,"%lf ",newdata[k][m][nit].p);	  
	    if(ic == 3) fprintf(outfile,"%lf ",newdata[k][m][nit].s);
	    if(ic == 4) fprintf(outfile,"%lf ",newdata[k][m][nit].smb);
	    if(ic == 5) fprintf(outfile,"%lf ",newdata[k][m][nit].sma);
	    if(ic == 6) fprintf(outfile,"%lf ",newdata[k][m][nit].c+newdata[k][m][nit].p+newdata[k][m][nit].s+newdata[k][m][nit].v+newdata[k][m][nit].i+newdata[k][m][nit].w+newdata[k][m][nit].u); 
	    if(ic == 7) fprintf(outfile,"%lf ",newdata[k][m][nit].u);

	}
        fprintf(outfile,"\n");

      }
      fclose(outfile);
    }
  }



  return;
}

/********************************************************************/
void output_updated_states2(){
  
  FILE *outfile; 
  int ih, ic, k, m, it, cnum, nit;
  char newpath[190], newpath1[190], fn[30], ayear[10], fout[190];
  double  flowcp,flowpc,flowpv,flowvp,flowvc,flowcv,flowsc,flowcs,flowsp,flowps;
  double flowsbh,flowvbh,flowsbh2,flowvbh2,flowsbh3;

  
  printf("output updated states 2\n");
  strcpy(newpath,PATH1);
  strcat(newpath,finput_dir);
  strcat(newpath,"updated_states/");


  /* cnum=18; */
  cnum=6;
  
  for (ic=0;ic<cnum;ic++){

  
    
     for (it=0;it<NEWNT;it++){ 




      nit=it;



  
      strcpy(newpath1,newpath);
      sprintf(ayear,"%d.txt",year[nit]);


      if(ic == 0) strcpy(fn,"gsbh1.");
      if(ic == 1) strcpy(fn,"gvbh1.");
      if(ic == 2) strcpy(fn,"gsbh2.");
      if(ic == 3) strcpy(fn,"gvbh2.");
      if(ic == 4) strcpy(fn,"gsbh3.");
      if(ic == 5) strcpy(fn,"gzdis.");
      if(ic == 6) strcpy(fn,"gvvvf.");
      if(ic == 7) strcpy(fn,"gsssf.");
      if(ic == 8) strcpy(fn,"gvvnf.");
      if(ic == 9) strcpy(fn,"gssnf.");
      if(ic == 10) strcpy(fn,"gallt.");
      if(ic == 11) strcpy(fn,"glunl.");
      if(ic == 12) strcpy(fn,"gluwl.");
      if(ic == 13) strcpy(fn,"ginwl.");
      if(ic == 14) strcpy(fn,"ginnl.");
      if(ic == 15) strcpy(fn,"gvhar.");
      if(ic == 16) strcpy(fn,"gshar.");
      if(ic == 17) strcpy(fn,"ginet.");

      strcat(fn,ayear);
      strcat(newpath1,fn);
      strcpy(fout,newpath1); 

      outfile=fopen(fout,"w");

      for (ih=0;ih<HEADLEN;ih++) fprintf(outfile,"%s",header[ih]);
      
 
      for (k=0;k<NY;k++){
	for (m=0;m<NX;m++){
	  

	  flowcp=newdata[k][m][nit].flowcp*garea[k][m]/1000.0/1000.0;
	  flowpc=newdata[k][m][nit].flowpc*garea[k][m]/1000.0/1000.0;
	  flowpv=newdata[k][m][nit].flowpv*garea[k][m]/1000.0/1000.0;
	  flowvp=newdata[k][m][nit].flowvp*garea[k][m]/1000.0/1000.0;
	  flowvc=newdata[k][m][nit].flowvc*garea[k][m]/1000.0/1000.0;
	  flowcv=newdata[k][m][nit].flowcv*garea[k][m]/1000.0/1000.0;
	  flowsc=newdata[k][m][nit].flowsc*garea[k][m]/1000.0/1000.0;
	  flowcs=newdata[k][m][nit].flowcs*garea[k][m]/1000.0/1000.0;
	  flowsp=newdata[k][m][nit].flowsp*garea[k][m]/1000.0/1000.0;
	  flowps=newdata[k][m][nit].flowps*garea[k][m]/1000.0/1000.0;
          flowsbh=ZEROVALUE;
	  if(newdata[k][m][nit].smb > ZEROVALUE) flowsbh=newdata[k][m][nit].sbh/newdata[k][m][nit].smb/1000.0/1000.0;
          flowvbh=ZEROVALUE;
	  if(dstatic[k][m].vba > ZEROVALUE) flowvbh=newdata[k][m][nit].vbh/dstatic[k][m].vba/1000.0/1000.0;
          flowsbh2=ZEROVALUE;
	  if(newdata[k][m][nit].smb > ZEROVALUE) flowsbh2=newdata[k][m][nit].sbh2/newdata[k][m][nit].smb/1000.0/1000.0;
          flowvbh2=ZEROVALUE;
	  if(dstatic[k][m].vba > ZEROVALUE) flowvbh2=newdata[k][m][nit].vbh2/dstatic[k][m].vba/1000.0/1000.0;

	  flowsbh3=ZEROVALUE;
	  if(newdata[k][m][nit].smb > ZEROVALUE) flowsbh3=newdata[k][m][nit].sbh3/newdata[k][m][nit].smb/garea[k][m];

	  if(ic == 0) fprintf(outfile,"%lf ",newdata[k][m][nit].sbh);
	  if(ic == 1) fprintf(outfile,"%lf ",newdata[k][m][nit].vbh);
	  if(ic == 2) fprintf(outfile,"%lf ",newdata[k][m][nit].sbh2);
	  if(ic == 3) fprintf(outfile,"%lf ",newdata[k][m][nit].vbh2);
	  if(ic == 4) fprintf(outfile,"%lf ",newdata[k][m][nit].sbh3);
	  if(ic == 5) fprintf(outfile,"%d ",newdata[k][m][nit].zdis);           
	  if(ic == 6){
	    if(dstatic[k][m].fnf == 1){
	      fprintf(outfile,"%lf ",newdata[k][m][nit].v);
	    }
	    else{
	      fprintf(outfile,"0.00000000 ");
	    }
	  }
	  if(ic == 7){
	    if(dstatic[k][m].fnf == 1){
	      fprintf(outfile,"%lf ",newdata[k][m][nit].s);
	    }
	    else{
	      fprintf(outfile,"0.00000000 ");
	    }
	  }
	  if(ic == 8){
	    if(dstatic[k][m].fnf == 0){
	      fprintf(outfile,"%lf ",newdata[k][m][nit].v);
	    }
	    else{
	      fprintf(outfile,"0.00000000 ");
	    }
	  }
	  if(ic == 9){
	    if(dstatic[k][m].fnf == 0){
	      fprintf(outfile,"%lf ",newdata[k][m][nit].s);
	    }
	    else{
	      fprintf(outfile,"0.00000000 ");
	    }
	  }
	  if(ic == 10) fprintf(outfile,"%lf ",(flowcp+flowpc+flowvp+flowvc+flowsc+flowsp)+(flowsbh+flowvbh+flowsbh2+flowvbh2)+(flowps+flowcs+flowpv+flowcv));

	  if(ic == 11) fprintf(outfile,"%lf ",(flowcp+flowpc+flowvp+flowvc+flowsc+flowsp)-(flowps+flowcs+flowpv+flowcv));

	  if(ic == 12) fprintf(outfile,"%lf ",(flowcp+flowpc+flowvp+flowvc+flowsc+flowsp)+(flowsbh+flowvbh+flowsbh2+flowvbh2)-(flowps+flowcs+flowpv+flowcv));
	  if(ic == 13) fprintf(outfile,"%lf ",(flowcp+flowpc+flowvp+flowvc+flowsc+flowsp)+(flowsbh+flowvbh+flowsbh2+flowvbh2));
	  if(ic == 14) fprintf(outfile,"%lf ",(flowcp+flowpc+flowvp+flowvc+flowsc+flowsp));
	  if(ic == 15) fprintf(outfile,"%lf ",flowvbh+flowvbh2);
	  if(ic == 16) fprintf(outfile,"%lf ",flowsbh+flowsbh2+flowsbh3);
	  if(ic == 17) fprintf(outfile,"%lf ",(flowvp+flowvc+flowsc+flowsp)-(flowps+flowcs+flowpv+flowcv)+(flowvbh+flowvbh2));


	}
        fprintf(outfile,"\n");

      }
      fclose(outfile);
    }
  }



  return;
}

/********************************************************************/
void output_updated_states3(){
  
  FILE *outfile; 
  int ih, ic, k, m, it, cnum, nit;
  char newpath[190], newpath1[190], fn[30], ayear[10], fout[190];
  double  flowcp,flowpc,flowpv,flowvp,flowvc,flowcv,flowsc,flowcs,flowsp,flowps;
  double flowsbh,flowvbh,flowsbh2,flowvbh2,flowsbh3;
  double flowcu,flowpu,flowvu,flowsu,flowuc,flowup,flowus;
  
  printf("output updated states 3\n");
  strcpy(newpath,PATH1);
  strcat(newpath,finput_dir);
  strcat(newpath,"updated_states/");


  cnum=20;
  
  
  for (ic=0;ic<cnum;ic++){

 
    
     for (it=0;it<NEWNT;it++){ 




      nit=it;


      strcpy(newpath1,newpath);
      sprintf(ayear,"%d.txt",year[nit]);

      if(ic == 0) strcpy(fn,"gflcp.");
      if(ic == 1) strcpy(fn,"gflpc.");
      if(ic == 2) strcpy(fn,"gflcs.");
      if(ic == 3) strcpy(fn,"gflsc.");
      if(ic == 4) strcpy(fn,"gflps.");
      if(ic == 5) strcpy(fn,"gflsp.");
      if(ic == 6) strcpy(fn,"gflvp.");
      if(ic == 7) strcpy(fn,"gflvc.");
      if(ic == 8) strcpy(fn,"gfvh1.");
      if(ic == 9) strcpy(fn,"gfvh2.");
      if(ic == 10) strcpy(fn,"gfsh1.");
      if(ic == 11) strcpy(fn,"gfsh2.");
      if(ic == 12) strcpy(fn,"gfsh3.");
      if(ic == 13) strcpy(fn,"gflcu.");
      if(ic == 14) strcpy(fn,"gflpu.");
      if(ic == 15) strcpy(fn,"gflvu.");
      if(ic == 16) strcpy(fn,"gflsu.");
      if(ic == 17) strcpy(fn,"gfluc.");
      if(ic == 18) strcpy(fn,"gflup.");
      if(ic == 19) strcpy(fn,"gflus.");

      strcat(fn,ayear);
      strcat(newpath1,fn);
      strcpy(fout,newpath1); 

      outfile=fopen(fout,"w");

      for (ih=0;ih<HEADLEN;ih++) fprintf(outfile,"%s",header[ih]);
      
 
      for (k=0;k<NY;k++){
	for (m=0;m<NX;m++){
	  

	  flowcp=newdata[k][m][nit].flowcp;
	  flowpc=newdata[k][m][nit].flowpc;
	  flowpv=newdata[k][m][nit].flowpv;
	  flowvp=newdata[k][m][nit].flowvp;
	  flowvc=newdata[k][m][nit].flowvc;
	  flowcv=newdata[k][m][nit].flowcv;
	  flowsc=newdata[k][m][nit].flowsc;
	  flowcs=newdata[k][m][nit].flowcs;
	  flowsp=newdata[k][m][nit].flowsp;
	  flowps=newdata[k][m][nit].flowps;

	  flowcu=newdata[k][m][nit].flowcu;
	  flowpu=newdata[k][m][nit].flowpu;
	  flowsu=newdata[k][m][nit].flowsu;
	  flowvu=newdata[k][m][nit].flowvu;
	  flowuc=newdata[k][m][nit].flowuc;
	  flowup=newdata[k][m][nit].flowup;
	  flowus=newdata[k][m][nit].flowus;

          flowsbh=ZEROVALUE;
	  if(newdata[k][m][nit].smb > ZEROVALUE) flowsbh=newdata[k][m][nit].sbh/newdata[k][m][nit].smb/garea[k][m];
          flowvbh=ZEROVALUE;
	  if(dstatic[k][m].vba > ZEROVALUE) flowvbh=newdata[k][m][nit].vbh/dstatic[k][m].vba/garea[k][m];
          flowsbh2=ZEROVALUE;
	  if(newdata[k][m][nit].smb > ZEROVALUE) flowsbh2=newdata[k][m][nit].sbh2/newdata[k][m][nit].smb/garea[k][m];
          flowvbh2=ZEROVALUE;
	  if(dstatic[k][m].vba > ZEROVALUE) flowvbh2=newdata[k][m][nit].vbh2/dstatic[k][m].vba/garea[k][m];
	  flowsbh3=ZEROVALUE;
	  if(newdata[k][m][nit].smb > ZEROVALUE) flowsbh3=newdata[k][m][nit].sbh3/newdata[k][m][nit].smb/garea[k][m];


	  if(ic == 0) fprintf(outfile,"%lf ",flowcp);
	  if(ic == 1) fprintf(outfile,"%lf ",flowpc);
	  if(ic == 2) fprintf(outfile,"%lf ",flowcs);
	  if(ic == 3) fprintf(outfile,"%lf ",flowsc);
	  if(ic == 4) fprintf(outfile,"%lf ",flowps);
	  if(ic == 5) fprintf(outfile,"%lf ",flowsp);
	  if(ic == 6) fprintf(outfile,"%lf ",flowvp);
	  if(ic == 7) fprintf(outfile,"%lf ",flowvc);
	  if(ic == 8) fprintf(outfile,"%lf ",flowvbh);
	  if(ic == 9) fprintf(outfile,"%lf ",flowvbh2);
	  if(ic == 10) fprintf(outfile,"%lf ",flowsbh);
	  if(ic == 11) fprintf(outfile,"%lf ",flowsbh2);
	  if(ic == 12) fprintf(outfile,"%lf ",flowsbh3);
	  if(ic == 13) fprintf(outfile,"%lf ",flowcu);
	  if(ic == 14) fprintf(outfile,"%lf ",flowpu);
	  if(ic == 15) fprintf(outfile,"%lf ",flowvu);
	  if(ic == 16) fprintf(outfile,"%lf ",flowsu);
	  if(ic == 17) fprintf(outfile,"%lf ",flowuc);
	  if(ic == 18) fprintf(outfile,"%lf ",flowup);
	  if(ic == 19) fprintf(outfile,"%lf ",flowus);

	}
        fprintf(outfile,"\n");

      }
      fclose(outfile);
    }
  }



  return;
}
/*******************************************************************/

void smart_flow(int k, int m, int it){

  double vstotal=ZEROVALUE;
  double delta_c=ZEROVALUE, delta_p=ZEROVALUE, delta_o=ZEROVALUE, delta_u=ZEROVALUE, unmet_du=ZEROVALUE;
  

  vstotal=newdata[k][m][it].v+newdata[k][m][it].s; 
  
   /* Urban flows first */
   delta_u = newdata[k][m][it+1].u - newdata[k][m][it].u;

   if(fabs(delta_u)<ZEROVALUE_CHECK){
     /*do nothing because delta_u is zero */ 
   }
   else if (delta_u<=ZEROVALUE){
     if((newdata[k][m][it+1].v + newdata[k][m][it+1].s) < -delta_u){
       newdata[k][m][it].flowus=newdata[k][m][it+1].v + newdata[k][m][it+1].s;
       unmet_du = -delta_u - newdata[k][m][it+1].v - newdata[k][m][it+1].s;
       newdata[k][m][it].flowuc=unmet_du*newdata[k][m][it].c/(newdata[k][m][it].c + newdata[k][m][it+1].p);
       newdata[k][m][it].flowup=unmet_du*newdata[k][m][it].p/(newdata[k][m][it].c + newdata[k][m][it+1].p);
     }
     else newdata[k][m][it].flowus=-delta_u;
   }
   else if (delta_u>ZEROVALUE){
     if((newdata[k][m][it].c+newdata[k][m][it].p+newdata[k][m][it].s)>=delta_u){
       newdata[k][m][it].flowcu=newdata[k][m][it].c/(newdata[k][m][it].c+newdata[k][m][it].p+newdata[k][m][it].s)*delta_u;
       newdata[k][m][it].flowpu=newdata[k][m][it].p/(newdata[k][m][it].c+newdata[k][m][it].p+newdata[k][m][it].s)*delta_u;
       newdata[k][m][it].flowsu=newdata[k][m][it].s/(newdata[k][m][it].c+newdata[k][m][it].p+newdata[k][m][it].s)*delta_u;
     }
     else{
       newdata[k][m][it].flowcu=newdata[k][m][it].c;
       newdata[k][m][it].flowpu=newdata[k][m][it].p;
       newdata[k][m][it].flowsu=newdata[k][m][it].s;
       newdata[k][m][it].flowvu=delta_u-newdata[k][m][it].c-newdata[k][m][it].p-newdata[k][m][it].s;
     }
   }

  delta_c = newdata[k][m][it+1].c - newdata[k][m][it].c;
  delta_p = newdata[k][m][it+1].p - newdata[k][m][it].p;
  delta_o = newdata[k][m][it+1].v + newdata[k][m][it+1].s -newdata[k][m][it].v-newdata[k][m][it].s;

   /*amend deltas */
  delta_c = delta_c+newdata[k][m][it].flowcu-newdata[k][m][it].flowus;
  delta_p = delta_p+newdata[k][m][it].flowpu-newdata[k][m][it].flowup;
  delta_o = delta_o-newdata[k][m][it].flowus+newdata[k][m][it].flowsu+newdata[k][m][it].flowvu;
 

  /*three cases*/
  if((fabs(delta_c)<ZEROVALUE_CHECK) && (fabs(delta_p)<ZEROVALUE_CHECK) && (fabs(delta_o)<ZEROVALUE_CHECK)){
    /*all zero*/
    
  }
  /*two negative*/
  else if((delta_c<=ZEROVALUE) && (delta_p<=ZEROVALUE) && (delta_o>ZEROVALUE)){
    
    newdata[k][m][it].flowps=-delta_p;
    newdata[k][m][it].flowcs=-delta_c;    
  }
  else if((delta_c<=ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o<=ZEROVALUE)){
    
    if(vstotal > ZEROVALUE){
      newdata[k][m][it].flowvp=-delta_o*newdata[k][m][it].v/vstotal;
      newdata[k][m][it].flowsp=-delta_o*(1-newdata[k][m][it].v/vstotal);
    }
    newdata[k][m][it].flowcp=-delta_c;
    
  }
  else if((delta_c>ZEROVALUE)&& (delta_p<=ZEROVALUE) && (delta_o<=ZEROVALUE)){
    
    newdata[k][m][it].flowpc=-delta_p; 
    if(vstotal > ZEROVALUE){
      newdata[k][m][it].flowvc=-delta_o*newdata[k][m][it].v/vstotal;
      newdata[k][m][it].flowsc=-delta_o*(1-newdata[k][m][it].v/vstotal);
    }
    
  } 
  /*two positive*/
  else if((delta_c<=ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o>ZEROVALUE)){

    newdata[k][m][it].flowcs=delta_o;
    newdata[k][m][it].flowcp=delta_p;
  }
  else if((delta_c>ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o<=ZEROVALUE)){

    if(vstotal > ZEROVALUE){
      newdata[k][m][it].flowvp=delta_p*newdata[k][m][it].v/vstotal;
      newdata[k][m][it].flowsp=delta_p*(1-newdata[k][m][it].v/vstotal);
      newdata[k][m][it].flowvc=delta_c*newdata[k][m][it].v/vstotal;
      newdata[k][m][it].flowsc=delta_c*(1-newdata[k][m][it].v/vstotal);
    }
    
  }
  else if((delta_c>ZEROVALUE)&& (delta_p<=ZEROVALUE) && (delta_o>ZEROVALUE)){

    newdata[k][m][it].flowps=delta_o;
    newdata[k][m][it].flowpc=delta_c;
    
  } 



# if SMART_FLOW_BUG_PRINT
  
    	/* checking for negative zero flows */

	if(newdata[k][m][it].flowcp < ZEROVALUE) {printf("Smart flow Bug: flowcp<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowcp);newdata[k][m][it].flowcp = ZEROVALUE;}
	if(newdata[k][m][it].flowpc < ZEROVALUE) {printf("Smart flow Bug: flowpc<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowpc);newdata[k][m][it].flowpc = ZEROVALUE;}
	if(newdata[k][m][it].flowpv < ZEROVALUE) {printf("Smart flow Bug: flowpv<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowpv);newdata[k][m][it].flowpv = ZEROVALUE;}
	if(newdata[k][m][it].flowvp < ZEROVALUE) {printf("Smart flow Bug: flowvp<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowvp);newdata[k][m][it].flowvp = ZEROVALUE;}
	if(newdata[k][m][it].flowvc < ZEROVALUE) {printf("Smart flow Bug: flowvc<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowvc);newdata[k][m][it].flowvc = ZEROVALUE;}
	if(newdata[k][m][it].flowcv < ZEROVALUE) {printf("Smart flow Bug: flowcv<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowcv);newdata[k][m][it].flowcv = ZEROVALUE;}
	if(newdata[k][m][it].flowsp < ZEROVALUE) {printf("Smart flow Bug: flowsp<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowsp);newdata[k][m][it].flowsp = ZEROVALUE;}
	if(newdata[k][m][it].flowps < ZEROVALUE) {printf("Smart flow Bug: flowps<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowps);newdata[k][m][it].flowps = ZEROVALUE;}
	if(newdata[k][m][it].flowsc < ZEROVALUE) {printf("Smart flow Bug: flowsc<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowsc);newdata[k][m][it].flowsc = ZEROVALUE;}
	if(newdata[k][m][it].flowcs < ZEROVALUE) {printf("Smart flow Bug: flowcs<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowcs);newdata[k][m][it].flowcs = ZEROVALUE;}
	if(newdata[k][m][it].flowvs < ZEROVALUE) {printf("Smart flow Bug: flowvs<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowvs);newdata[k][m][it].flowvs = ZEROVALUE;}
	if(newdata[k][m][it].flowcu < ZEROVALUE) {printf("Smart flow Bug: flowcu<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowcu);newdata[k][m][it].flowcu = ZEROVALUE;}
	if(newdata[k][m][it].flowpu < ZEROVALUE) {printf("Smart flow Bug: flowpu<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowpu);newdata[k][m][it].flowpu = ZEROVALUE;}
	if(newdata[k][m][it].flowvu < ZEROVALUE) {printf("Smart flow Bug: flowvu<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowvu);newdata[k][m][it].flowvu = ZEROVALUE;}
	if(newdata[k][m][it].flowsu < ZEROVALUE) {printf("Smart flow Bug: flowsu<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowsu);newdata[k][m][it].flowsu = ZEROVALUE;}
	if(newdata[k][m][it].flowuc < ZEROVALUE) {printf("Smart flow Bug: flowuc<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowuc);newdata[k][m][it].flowuc = ZEROVALUE;}
	if(newdata[k][m][it].flowup < ZEROVALUE) {printf("Smart flow Bug: flowup<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowup);newdata[k][m][it].flowup = ZEROVALUE;}
	if(newdata[k][m][it].flowus < ZEROVALUE) {printf("Smart flow Bug: flowus<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowus);newdata[k][m][it].flowus = ZEROVALUE;}


#endif
    

  return;
  
}
 

/********************************************************/
 void alternative_smart_flow1(int k, int m, int it){

   double vtotal=ZEROVALUE, stotal=ZEROVALUE, vstotal=ZEROVALUE;
   double p_other=ZEROVALUE, c_other=ZEROVALUE, total_demand=ZEROVALUE;
   double delta_c=ZEROVALUE, delta_p=ZEROVALUE, delta_o=ZEROVALUE, delta_u=ZEROVALUE, unmet_du=ZEROVALUE;



   /* Urban flows first */
   delta_u = newdata[k][m][it+1].u - newdata[k][m][it].u;

   if(fabs(delta_u)<=ZEROVALUE_CHECK){
     /*do nothing because delta_u is zero */ 
   }
   else if (delta_u<=ZEROVALUE){
     newdata[k][m][it].flowuc=newdata[k][m][it+1].c/(newdata[k][m][it+1].c+newdata[k][m][it+1].p+newdata[k][m][it+1].s+newdata[k][m][it+1].v)*(-delta_u);
     newdata[k][m][it].flowup=newdata[k][m][it+1].p/(newdata[k][m][it+1].c+newdata[k][m][it+1].p+newdata[k][m][it+1].s+newdata[k][m][it+1].v)*(-delta_u);
     newdata[k][m][it].flowus=(newdata[k][m][it+1].s+newdata[k][m][it+1].v)/(newdata[k][m][it+1].c+newdata[k][m][it+1].p+newdata[k][m][it+1].s+newdata[k][m][it+1].v)*(-delta_u);
   }
   else if (delta_u>ZEROVALUE){
     if((newdata[k][m][it].c+newdata[k][m][it].p+newdata[k][m][it].s)>=delta_u){
       newdata[k][m][it].flowcu=newdata[k][m][it].c/(newdata[k][m][it].c+newdata[k][m][it].p+newdata[k][m][it].s)*delta_u;
       newdata[k][m][it].flowpu=newdata[k][m][it].p/(newdata[k][m][it].c+newdata[k][m][it].p+newdata[k][m][it].s)*delta_u;
       newdata[k][m][it].flowsu=newdata[k][m][it].s/(newdata[k][m][it].c+newdata[k][m][it].p+newdata[k][m][it].s)*delta_u;
     }
     else{
       newdata[k][m][it].flowcu=newdata[k][m][it].c;
       newdata[k][m][it].flowpu=newdata[k][m][it].p;
       newdata[k][m][it].flowsu=newdata[k][m][it].s;
       newdata[k][m][it].flowvu=delta_u-newdata[k][m][it].c-newdata[k][m][it].p-newdata[k][m][it].s;
     }
   }


  delta_c = newdata[k][m][it+1].c - newdata[k][m][it].c;
  delta_p = newdata[k][m][it+1].p - newdata[k][m][it].p;
  delta_o = newdata[k][m][it+1].v + newdata[k][m][it+1].s -newdata[k][m][it].v-newdata[k][m][it].s;

  /*amend deltas */
  delta_c = delta_c+newdata[k][m][it].flowcu-newdata[k][m][it].flowuc;
  delta_p = delta_p+newdata[k][m][it].flowpu-newdata[k][m][it].flowup;
  delta_o = delta_o-newdata[k][m][it].flowus+newdata[k][m][it].flowsu+newdata[k][m][it].flowvu;

    
  /*three cases*/
  if((fabs(delta_c)<ZEROVALUE_CHECK) && (fabs(delta_p)<ZEROVALUE_CHECK) && (fabs(delta_o)<ZEROVALUE_CHECK)){
    /*all zero*/
    
  }
  /*two negative*/
  else if((delta_c<=ZEROVALUE) && (delta_p<=ZEROVALUE) && (delta_o>ZEROVALUE)){
    
    newdata[k][m][it].flowps=-delta_p;
    newdata[k][m][it].flowcs=-delta_c;    
  }
  else if((delta_c<=ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o<=ZEROVALUE)){

    newdata[k][m][it].flowcp=-delta_c; 


    if((newdata[k][m][it].v-newdata[k][m][it].flowvu) >= -delta_o) {

      newdata[k][m][it].flowvp=-delta_o;
    }
    else {
      newdata[k][m][it].flowvp=newdata[k][m][it].v-newdata[k][m][it].flowvu;
      newdata[k][m][it].flowsp=-delta_o-(newdata[k][m][it].v-newdata[k][m][it].flowvu);
    }

    
  }
  else if((delta_c>ZEROVALUE)&& (delta_p<=ZEROVALUE) && (delta_o<=ZEROVALUE)){
    
    newdata[k][m][it].flowpc=-delta_p; 


    if((newdata[k][m][it].v-newdata[k][m][it].flowvu) >= -delta_o) {

      newdata[k][m][it].flowvc=-delta_o;
    }
    else {
      newdata[k][m][it].flowvc=newdata[k][m][it].v-newdata[k][m][it].flowvu;
      newdata[k][m][it].flowsc=-delta_o-(newdata[k][m][it].v-newdata[k][m][it].flowvu);
    }
    

  } 
  /*two positive*/
  else if((delta_c<=ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o>ZEROVALUE)){

    newdata[k][m][it].flowcs=delta_o;
    newdata[k][m][it].flowcp=delta_p;
  }
  else if((delta_c>ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o<=ZEROVALUE)){

 
    if((newdata[k][m][it].v-newdata[k][m][it].flowvu) >= -delta_o) {

      newdata[k][m][it].flowvc=delta_c;
      newdata[k][m][it].flowvp=delta_p;

    }
    else {


      newdata[k][m][it].flowvc=(newdata[k][m][it].v-newdata[k][m][it].flowvu)*delta_c/(delta_p+delta_c);
      newdata[k][m][it].flowvp=(newdata[k][m][it].v-newdata[k][m][it].flowvu)*delta_p/(delta_p+delta_c);
    
      newdata[k][m][it].flowsc=delta_c-newdata[k][m][it].flowvc;
      newdata[k][m][it].flowsp=delta_p-newdata[k][m][it].flowvp;
    }
    

  }
  else if((delta_c>ZEROVALUE)&& (delta_p<=ZEROVALUE) && (delta_o>ZEROVALUE)){

    newdata[k][m][it].flowps=delta_o;
    newdata[k][m][it].flowpc=delta_c;
    
  } 

  


# if SMART_FLOW_BUG_PRINT


    	/* checking for negative zero flows */

	if(newdata[k][m][it].flowcp < ZEROVALUE) {printf("Smart flow Bug: flowcp<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowcp);newdata[k][m][it].flowcp = ZEROVALUE;}
	if(newdata[k][m][it].flowpc < ZEROVALUE) {printf("Smart flow Bug: flowpc<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowpc);newdata[k][m][it].flowpc = ZEROVALUE;}
	if(newdata[k][m][it].flowpv < ZEROVALUE) {printf("Smart flow Bug: flowpv<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowpv);newdata[k][m][it].flowpv = ZEROVALUE;}
	if(newdata[k][m][it].flowvp < ZEROVALUE) {printf("Smart flow Bug: flowvp<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowvp);newdata[k][m][it].flowvp = ZEROVALUE;}
	if(newdata[k][m][it].flowvc < ZEROVALUE) {printf("Smart flow Bug: flowvc<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowvc);newdata[k][m][it].flowvc = ZEROVALUE;}
	if(newdata[k][m][it].flowcv < ZEROVALUE) {printf("Smart flow Bug: flowcv<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowcv);newdata[k][m][it].flowcv = ZEROVALUE;}
	if(newdata[k][m][it].flowsp < ZEROVALUE) {printf("Smart flow Bug: flowsp<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowsp);newdata[k][m][it].flowsp = ZEROVALUE;}
	if(newdata[k][m][it].flowps < ZEROVALUE) {printf("Smart flow Bug: flowps<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowps);newdata[k][m][it].flowps = ZEROVALUE;}
	if(newdata[k][m][it].flowsc < ZEROVALUE) {printf("Smart flow Bug: flowsc<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowsc);newdata[k][m][it].flowsc = ZEROVALUE;}
	if(newdata[k][m][it].flowcs < ZEROVALUE) {printf("Smart flow Bug: flowcs<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowcs);newdata[k][m][it].flowcs = ZEROVALUE;}
	if(newdata[k][m][it].flowvs < ZEROVALUE) {printf("Smart flow Bug: flowvs<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowvs);newdata[k][m][it].flowvs = ZEROVALUE;}
	if(newdata[k][m][it].flowcu < ZEROVALUE) {printf("Smart flow Bug: flowcu<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowcu);newdata[k][m][it].flowcu = ZEROVALUE;}
	if(newdata[k][m][it].flowpu < ZEROVALUE) {printf("Smart flow Bug: flowpu<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowpu);newdata[k][m][it].flowpu = ZEROVALUE;}
	if(newdata[k][m][it].flowvu < ZEROVALUE) {printf("Smart flow Bug: flowvu<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowvu);newdata[k][m][it].flowvu = ZEROVALUE;}
	if(newdata[k][m][it].flowsu < ZEROVALUE) {printf("Smart flow Bug: flowsu<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowsu);newdata[k][m][it].flowsu = ZEROVALUE;}
	if(newdata[k][m][it].flowuc < ZEROVALUE) {printf("Smart flow Bug: flowuc<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowuc);newdata[k][m][it].flowuc = ZEROVALUE;}
	if(newdata[k][m][it].flowup < ZEROVALUE) {printf("Smart flow Bug: flowup<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowup);newdata[k][m][it].flowup = ZEROVALUE;}
	if(newdata[k][m][it].flowus < ZEROVALUE) {printf("Smart flow Bug: flowus<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowus);newdata[k][m][it].flowus = ZEROVALUE;}


#endif
    


  return;
  
}


/********************************************************************/
void alternative_smart_flow2(int k, int m, int it){

  double vtotal=ZEROVALUE, stotal=ZEROVALUE, vstotal=ZEROVALUE;
  double p_other=ZEROVALUE, c_other=ZEROVALUE, total_demand=ZEROVALUE;
  double delta_c=ZEROVALUE, delta_p=ZEROVALUE, delta_o=ZEROVALUE, delta_u=ZEROVALUE, unmet_du=ZEROVALUE;



   /* Urban flows first */
   delta_u = newdata[k][m][it+1].u - newdata[k][m][it].u;

   if(fabs(delta_u)<ZEROVALUE_CHECK){
     /*do nothing because delta_u is zero */ 
   }
   else if (delta_u<=ZEROVALUE){
     newdata[k][m][it].flowuc=newdata[k][m][it+1].c/(newdata[k][m][it+1].c+newdata[k][m][it+1].p+newdata[k][m][it+1].s+newdata[k][m][it+1].v)*(-delta_u);
     newdata[k][m][it].flowup=newdata[k][m][it+1].p/(newdata[k][m][it+1].c+newdata[k][m][it+1].p+newdata[k][m][it+1].s+newdata[k][m][it+1].v)*(-delta_u);
     newdata[k][m][it].flowus=(newdata[k][m][it+1].s+newdata[k][m][it+1].v)/(newdata[k][m][it+1].c+newdata[k][m][it+1].p+newdata[k][m][it+1].s+newdata[k][m][it+1].v)*(-delta_u);
   }
   else if (delta_u>ZEROVALUE){
     if((newdata[k][m][it].c+newdata[k][m][it].p+newdata[k][m][it].s)>=delta_u){
       newdata[k][m][it].flowcu=newdata[k][m][it].c/(newdata[k][m][it].c+newdata[k][m][it].p+newdata[k][m][it].s)*delta_u;
       newdata[k][m][it].flowpu=newdata[k][m][it].p/(newdata[k][m][it].c+newdata[k][m][it].p+newdata[k][m][it].s)*delta_u;
       newdata[k][m][it].flowsu=newdata[k][m][it].s/(newdata[k][m][it].c+newdata[k][m][it].p+newdata[k][m][it].s)*delta_u;
     }
     else{
       newdata[k][m][it].flowcu=newdata[k][m][it].c;
       newdata[k][m][it].flowpu=newdata[k][m][it].p;
       newdata[k][m][it].flowsu=newdata[k][m][it].s;
       newdata[k][m][it].flowvu=delta_u-newdata[k][m][it].c-newdata[k][m][it].p-newdata[k][m][it].s;
     }
   }


  delta_c = newdata[k][m][it+1].c - newdata[k][m][it].c;
  delta_p = newdata[k][m][it+1].p - newdata[k][m][it].p;
  delta_o = newdata[k][m][it+1].v + newdata[k][m][it+1].s -newdata[k][m][it].v-newdata[k][m][it].s;

  /*amend deltas */
  delta_c = delta_c+newdata[k][m][it].flowcu-newdata[k][m][it].flowuc;
  delta_p = delta_p+newdata[k][m][it].flowpu-newdata[k][m][it].flowup;
  delta_o = delta_o-newdata[k][m][it].flowus+newdata[k][m][it].flowsu+newdata[k][m][it].flowvu;

  
  /*three cases*/
  if((fabs(delta_c)<ZEROVALUE_CHECK) && (fabs(delta_p)<ZEROVALUE_CHECK) && (fabs(delta_o)<ZEROVALUE_CHECK)){
    /*all zero*/
    
  }
  /*two negative*/
  else if((delta_c<=ZEROVALUE) && (delta_p<=ZEROVALUE) && (delta_o>ZEROVALUE)){
    
    newdata[k][m][it].flowps=-delta_p;
    newdata[k][m][it].flowcs=-delta_c;    
  }
  else if((delta_c<=ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o<=ZEROVALUE)){

    newdata[k][m][it].flowcp=-delta_c; 


    if((newdata[k][m][it].s-newdata[k][m][it].flowsu) >= -delta_o) {

      newdata[k][m][it].flowsp=-delta_o;
    }
    else {
      newdata[k][m][it].flowsp=newdata[k][m][it].s-newdata[k][m][it].flowsu;
      newdata[k][m][it].flowvp=-delta_o-(newdata[k][m][it].s-newdata[k][m][it].flowsu);
    }

    
  }
  else if((delta_c>ZEROVALUE)&& (delta_p<=ZEROVALUE) && (delta_o<=ZEROVALUE)){
    
    newdata[k][m][it].flowpc=-delta_p; 


    if((newdata[k][m][it].s-newdata[k][m][it].flowsu) >= -delta_o) {

      newdata[k][m][it].flowsc=-delta_o;
    }
    else {
      newdata[k][m][it].flowsc=newdata[k][m][it].s-newdata[k][m][it].flowsu;
      newdata[k][m][it].flowvc=-delta_o-(newdata[k][m][it].s-newdata[k][m][it].flowsu);
    }
    

  } 
  /*two positive*/
  else if((delta_c<=ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o>ZEROVALUE)){

    newdata[k][m][it].flowcs=delta_o;
    newdata[k][m][it].flowcp=delta_p;
  }
  else if((delta_c>ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o<=ZEROVALUE)){

 
    if((newdata[k][m][it].s-newdata[k][m][it].flowsu) >= -delta_o) {

      newdata[k][m][it].flowsc=delta_c;
      newdata[k][m][it].flowsp=delta_p;

    }
    else {

      newdata[k][m][it].flowsc=(newdata[k][m][it].s-newdata[k][m][it].flowsu)*delta_c/(delta_p+delta_c);
      newdata[k][m][it].flowsp=(newdata[k][m][it].s-newdata[k][m][it].flowsu)*delta_p/(delta_p+delta_c);
    
      newdata[k][m][it].flowvc=delta_c-newdata[k][m][it].flowsc;
      newdata[k][m][it].flowvp=delta_p-newdata[k][m][it].flowsp;
    }
    

  }
  else if((delta_c>ZEROVALUE)&& (delta_p<=ZEROVALUE) && (delta_o>ZEROVALUE)){

    newdata[k][m][it].flowps=delta_o;
    newdata[k][m][it].flowpc=delta_c;
    
  } 



# if SMART_FLOW_BUG_PRINT

 	if((newdata[k][m][it].flowvc < ZEROVALUE) || (newdata[k][m][it].flowvp < ZEROVALUE) || (newdata[k][m][it].flowsc < ZEROVALUE) || (newdata[k][m][it].flowsp < ZEROVALUE)) printf("v(t) %lf s(t) %lf c(t) %lf p(t) %lf v(t+1) %lf s(t+1) %lf c(t+1) %lf p(t+1) %lf k %d m %d vc %lf vp %lf sc %lf sp %lf delta_p %lf delta_o %lf delta_c %lf\n",newdata[k][m][it].v,newdata[k][m][it].s,newdata[k][m][it].c,newdata[k][m][it].p,newdata[k][m][it+1].v,newdata[k][m][it+1].s,newdata[k][m][it+1].c,newdata[k][m][it+1].p,k,m,newdata[k][m][it].flowvc,newdata[k][m][it].flowvp,newdata[k][m][it].flowsc,newdata[k][m][it].flowsp,delta_p,delta_o,delta_c); 
  
    	/* checking for negative zero flows */

	if(newdata[k][m][it].flowcp < ZEROVALUE) {printf("Smart flow Bug: flowcp<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[k][m][it].flowcp);newdata[k][m][it].flowcp = ZEROVALUE;}
	if(newdata[k][m][it].flowpc < ZEROVALUE) {printf("Smart flow Bug: flowpc<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[k][m][it].flowpc);newdata[k][m][it].flowpc = ZEROVALUE;}
	if(newdata[k][m][it].flowpv < ZEROVALUE) {printf("Smart flow Bug: flowpv<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[k][m][it].flowpv);newdata[k][m][it].flowpv = ZEROVALUE;}
	if(newdata[k][m][it].flowvp < ZEROVALUE) {printf("Smart flow Bug: flowvp<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[k][m][it].flowvp);newdata[k][m][it].flowvp = ZEROVALUE;}
	if(newdata[k][m][it].flowvc < ZEROVALUE) {printf("Smart flow Bug: flowvc<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[k][m][it].flowvc);newdata[k][m][it].flowvc = ZEROVALUE;}
	if(newdata[k][m][it].flowcv < ZEROVALUE) {printf("Smart flow Bug: flowcv<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[k][m][it].flowcv);newdata[k][m][it].flowcv = ZEROVALUE;}
	if(newdata[k][m][it].flowsp < ZEROVALUE) {printf("Smart flow Bug: flowsp<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[k][m][it].flowsp);newdata[k][m][it].flowsp = ZEROVALUE;}
	if(newdata[k][m][it].flowps < ZEROVALUE) {printf("Smart flow Bug: flowps<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[k][m][it].flowps);newdata[k][m][it].flowps = ZEROVALUE;}
	if(newdata[k][m][it].flowsc < ZEROVALUE) {printf("Smart flow Bug: flowsc<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[k][m][it].flowsc);newdata[k][m][it].flowsc = ZEROVALUE;}
	if(newdata[k][m][it].flowcs < ZEROVALUE) {printf("Smart flow Bug: flowcs<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[k][m][it].flowcs);newdata[k][m][it].flowcs = ZEROVALUE;}
	if(newdata[k][m][it].flowvs < ZEROVALUE) {printf("Smart flow Bug: flowvs<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[k][m][it].flowvs);newdata[k][m][it].flowvs = ZEROVALUE;}
	if(newdata[k][m][it].flowcu < ZEROVALUE) {printf("Smart flow Bug: flowcu<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowcu);newdata[k][m][it].flowcu = ZEROVALUE;}
	if(newdata[k][m][it].flowpu < ZEROVALUE) {printf("Smart flow Bug: flowpu<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowpu);newdata[k][m][it].flowpu = ZEROVALUE;}
	if(newdata[k][m][it].flowvu < ZEROVALUE) {printf("Smart flow Bug: flowvu<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowvu);newdata[k][m][it].flowvu = ZEROVALUE;}
	if(newdata[k][m][it].flowsu < ZEROVALUE) {printf("Smart flow Bug: flowsu<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowsu);newdata[k][m][it].flowsu = ZEROVALUE;}
	if(newdata[k][m][it].flowuc < ZEROVALUE) {printf("Smart flow Bug: flowuc<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowuc);newdata[k][m][it].flowuc = ZEROVALUE;}
	if(newdata[k][m][it].flowup < ZEROVALUE) {printf("Smart flow Bug: flowup<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowup);newdata[k][m][it].flowup = ZEROVALUE;}
	if(newdata[k][m][it].flowus < ZEROVALUE) {printf("Smart flow Bug: flowus<0, resetting to 0 by force %10.15lf\n",newdata[k][m][it].flowus);newdata[k][m][it].flowus = ZEROVALUE;}

#endif
    

  return;
  
}

/********************************************************************/
void adjust_smart_flow1(int k, int m, int it, int i){

  double vtotal=ZEROVALUE, stotal=ZEROVALUE, vstotal=ZEROVALUE;
  double flowcs_prime=ZEROVALUE, flowps_prime=ZEROVALUE;
  double flowvc_prime=ZEROVALUE, flowvp_prime=ZEROVALUE;
  double flowsc_prime=ZEROVALUE, flowsp_prime=ZEROVALUE;
  double time_k, crop, past;




#if 1
  vtotal=newdata[k][m][it].v-newdata[k][m][it].flowvc-newdata[k][m][it].flowvp-newdata[k][m][it].flowvu; 
 
  stotal=newdata[k][m][it].s-newdata[k][m][it].flowsc-newdata[k][m][it].flowsp-newdata[k][m][it].flowsu;

  vstotal=vtotal+stotal;

  crop=newdata[k][m][it].c-newdata[k][m][it].flowcs-newdata[k][m][it].flowcp-newdata[k][m][it].flowcu;

  past=newdata[k][m][it].p-newdata[k][m][it].flowps-newdata[k][m][it].flowpc-newdata[k][m][it].flowpu;
#endif
 
  if((crop+past) > ZEROVALUE){


    /* abandonment of crop of pasture */


    /* best case method, desired abandonment = 1/15 per year */
    /* determine time_k */
    time_k=1.0/15.0;
    flowcs_prime=crop*time_k; 
    flowps_prime=past*time_k;

  
  
    if((flowcs_prime+flowps_prime) <= vstotal){

      if((flowcs_prime+flowps_prime) <= vtotal) {

	flowvc_prime=flowcs_prime;
	flowvp_prime=flowps_prime;
    
      }
      else { /* not enough v, but enough v+s */

	flowvc_prime=flowcs_prime/(flowcs_prime+flowps_prime)*vtotal;
	flowvp_prime=flowps_prime/(flowcs_prime+flowps_prime)*vtotal;
	
	flowsc_prime=flowcs_prime-flowvc_prime;
	flowsp_prime=flowps_prime-flowvp_prime;
	
      }

    }
    else{ /* not enough v+s, take all v+s (everything) that is possible */

      flowvc_prime=flowcs_prime/(flowcs_prime+flowps_prime)*vtotal;
      flowvp_prime=flowps_prime/(flowcs_prime+flowps_prime)*vtotal;
      
      flowsc_prime=flowcs_prime/(flowcs_prime+flowps_prime)*stotal;
      flowsp_prime=flowps_prime/(flowcs_prime+flowps_prime)*stotal;
      

      /* re-define desired abandonment (flowcs_prime and flowps_prime) to be what is available */

      flowcs_prime=flowvc_prime+flowsc_prime; 
      flowps_prime=flowvp_prime+flowsp_prime;

    }

 
 
    ctdata[i][it].flowvc_prime+=flowvc_prime*garea[k][m]/1000./1000.;
    ctdata[i][it].flowvp_prime+=flowvp_prime*garea[k][m]/1000./1000.;
    ctdata[i][it].flowsp_prime+=flowsp_prime*garea[k][m]/1000./1000.;
    ctdata[i][it].flowsc_prime+=flowsc_prime*garea[k][m]/1000./1000.;
    ctdata[i][it].flowps_prime+=flowps_prime*garea[k][m]/1000./1000.;
    ctdata[i][it].flowcs_prime+=flowcs_prime*garea[k][m]/1000./1000.;  

    newdata[k][m][it].flowvp+=flowvp_prime;
    newdata[k][m][it].flowvc+=flowvc_prime;
    newdata[k][m][it].flowsp+=flowsp_prime;
    newdata[k][m][it].flowsc+=flowsc_prime;
    newdata[k][m][it].flowps+=flowps_prime;
    newdata[k][m][it].flowcs+=flowcs_prime;



  }    


  return;
  
}
/********************************************************************/
void adjust_smart_flow2(int k, int m, int it, int i){

  double vtotal=ZEROVALUE, stotal=ZEROVALUE, vstotal=ZEROVALUE;
  double flowcs_prime=ZEROVALUE, flowps_prime=ZEROVALUE;
  double flowvc_prime=ZEROVALUE, flowvp_prime=ZEROVALUE;
  double flowsc_prime=ZEROVALUE, flowsp_prime=ZEROVALUE;
  double time_k, crop, past;



#if 1
  vtotal=newdata[k][m][it].v-newdata[k][m][it].flowvc-newdata[k][m][it].flowvp-newdata[k][m][it].flowvu; 
 
  stotal=newdata[k][m][it].s-newdata[k][m][it].flowsc-newdata[k][m][it].flowsp-newdata[k][m][it].flowsu;

  vstotal=vtotal+stotal;

  crop=newdata[k][m][it].c-newdata[k][m][it].flowcs-newdata[k][m][it].flowcp-newdata[k][m][it].flowcu;

  past=newdata[k][m][it].p-newdata[k][m][it].flowps-newdata[k][m][it].flowpc-newdata[k][m][it].flowpu;
#endif
 
  if((crop+past) > ZEROVALUE){

    /* abandonment of crop of pasture */



    /* best case method, desired abandonment = 1/15 per year */
    /* determine time_k */
    time_k=1.0/15.0;
    flowcs_prime=crop*time_k; 
    flowps_prime=past*time_k;

  
    if((flowcs_prime+flowps_prime) <= vstotal){

      if((flowcs_prime+flowps_prime) <= stotal) {

	flowsc_prime=flowcs_prime;
	flowsp_prime=flowps_prime;
    
      }
      else { /* not enough s, but enough v+s */

	flowsc_prime=flowcs_prime/(flowcs_prime+flowps_prime)*stotal;
	flowsp_prime=flowps_prime/(flowcs_prime+flowps_prime)*stotal;
	
	flowvc_prime=flowcs_prime-flowsc_prime;
	flowvp_prime=flowps_prime-flowsp_prime;
	
      }

    }
    else{ /* not enough v+s, take all v+s (everything) that is possible */

      flowvc_prime=flowcs_prime/(flowcs_prime+flowps_prime)*vtotal;
      flowvp_prime=flowps_prime/(flowcs_prime+flowps_prime)*vtotal;
      
      flowsc_prime=flowcs_prime/(flowcs_prime+flowps_prime)*stotal;
      flowsp_prime=flowps_prime/(flowcs_prime+flowps_prime)*stotal;     

      /* re-define desired abandonment (flowcs_prime and flowps_prime) to be what is available */

      flowcs_prime=flowvc_prime+flowsc_prime; 
      flowps_prime=flowvp_prime+flowsp_prime;

      
    }

    ctdata[i][it].flowvc_prime+=flowvc_prime*garea[k][m]/1000./1000.;
    ctdata[i][it].flowvp_prime+=flowvp_prime*garea[k][m]/1000./1000.;
    ctdata[i][it].flowsp_prime+=flowsp_prime*garea[k][m]/1000./1000.;
    ctdata[i][it].flowsc_prime+=flowsc_prime*garea[k][m]/1000./1000.;
    ctdata[i][it].flowps_prime+=flowps_prime*garea[k][m]/1000./1000.;
    ctdata[i][it].flowcs_prime+=flowcs_prime*garea[k][m]/1000./1000.; 

    newdata[k][m][it].flowvp+=flowvp_prime;
    newdata[k][m][it].flowvc+=flowvc_prime;
    newdata[k][m][it].flowsp+=flowsp_prime;
    newdata[k][m][it].flowsc+=flowsc_prime;
    newdata[k][m][it].flowps+=flowps_prime;
    newdata[k][m][it].flowcs+=flowcs_prime;




  }    



  return;
  
}


/********************************************************************/
void adjust_smart_flow3(int k, int m, int it, int i){

  FILE *outfile;
  double vtotal=ZEROVALUE, stotal=ZEROVALUE, vstotal=ZEROVALUE;
  double flowcs_prime=ZEROVALUE, flowps_prime=ZEROVALUE;
  double flowvc_prime=ZEROVALUE, flowvp_prime=ZEROVALUE;
  double flowsc_prime=ZEROVALUE, flowsp_prime=ZEROVALUE;
  double time_k, crop=ZEROVALUE, past=ZEROVALUE, time_fallow_max=ZEROVALUE, time_fallow_min=ZEROVALUE;
  double area_secondary_peryear=ZEROVALUE, area_secondary_fallowed=ZEROVALUE;
  double time_fallow=ZEROVALUE;




  /* secondary is the priority, only abandoning crop for now */


  vtotal=newdata[k][m][it].v-newdata[k][m][it].flowvc-newdata[k][m][it].flowvp; 
    
  stotal=newdata[k][m][it].s+newdata[k][m][it].flowcs+newdata[k][m][it].flowps-newdata[k][m][it].flowsc-newdata[k][m][it].flowsp;

  vstotal=vtotal+stotal;

  crop=newdata[k][m][it].c+newdata[k][m][it].flowvc+newdata[k][m][it].flowsc+newdata[k][m][it].flowpc-newdata[k][m][it].flowcs-newdata[k][m][it].flowcp;

  past=newdata[k][m][it].p+newdata[k][m][it].flowvp+newdata[k][m][it].flowsp+newdata[k][m][it].flowcp-newdata[k][m][it].flowps-newdata[k][m][it].flowpc;


  if( (vstotal > ZEROVALUE) && ((crop > ZEROVALUE) && (crop < 0.10)) ){

    /* time before abandoning cultivated land */
    time_k=4.0;
    
    /* fraction v+s required to replace crop abandoned */
    
    if(newdata[k][m][it].sma < 1.0) {
      time_fallow=time_k;
    }
    else{
      time_fallow=vstotal/(crop/time_k+stotal/(2.0*newdata[k][m][it].sma-1.0));
    }
    
    /* setting a max, allow us to return to fallowed secondary every 20 years */
    time_fallow_max=20.0;
    
    /* setting a min, allow us to return to fallowed secondary every time_k years */
    time_fallow_min=time_k;

    
    if(time_fallow > time_fallow_max){
      
      time_fallow=time_fallow_max;
      
      if(newdata[k][m][it].sma < 1.0) {
	area_secondary_fallowed=ZEROVALUE;
      }
      else{
	area_secondary_fallowed=stotal*(1.0-time_fallow/(2.0*newdata[k][m][it].sma-1.0));
      }
      
    }
    else if(time_fallow > time_k){
      
      if(newdata[k][m][it].sma < 1.0) {
	area_secondary_fallowed=ZEROVALUE;
      }
      else{
	area_secondary_fallowed=stotal*(1.0-time_fallow/(2.0*newdata[k][m][it].sma-1.0));
      }
      
    }
    else{ /*less than 4.0, slow down abandonment */
      
      time_fallow=time_k;    
      
      if(newdata[k][m][it].sma < 1.0) {
	area_secondary_fallowed=ZEROVALUE;
      }
      else{
	area_secondary_fallowed=stotal*(1.0-time_fallow/(2.0*newdata[k][m][it].sma-1.0));
      }
      
    }
    
    if(area_secondary_fallowed < ZEROVALUE) area_secondary_fallowed=ZEROVALUE;
    
    
    if((crop/time_k) < (vtotal+area_secondary_fallowed)){
      
      flowcs_prime=crop/time_k;
    }
    else{
      
      flowcs_prime=vtotal+area_secondary_fallowed;
    }
    
    if(area_secondary_fallowed < flowcs_prime){

      flowsc_prime=area_secondary_fallowed;
    }
    else{
      
      flowsc_prime=flowcs_prime;
    }
    
    flowvc_prime=flowcs_prime-flowsc_prime;
    
    
    
    
    ctdata[i][it].flowvc_prime+=flowvc_prime*garea[k][m]/1000.0/1000.0;
    ctdata[i][it].flowsc_prime+=flowsc_prime*garea[k][m]/1000.0/1000.0;
    ctdata[i][it].flowcs_prime+=flowcs_prime*garea[k][m]/1000.0/1000.0;
    
    newdata[k][m][it].flowvc+=flowvc_prime;
    newdata[k][m][it].flowsc+=flowsc_prime;
    newdata[k][m][it].flowcs+=flowcs_prime;


   
  } /* vstotal > zerovalue endif */

  return;
  
}
/********************************************************************/

void initialize_checker(int country_code){

  int k, m, it, i;


  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){


      for (it=0;it<NEWNT;it++){
	
	chtest[k][m][it].c=ZEROVALUE;
	chtest[k][m][it].p=ZEROVALUE;
	chtest[k][m][it].v=ZEROVALUE;
	chtest[k][m][it].s=ZEROVALUE;
	chtest[k][m][it].u=ZEROVALUE;
      
      }

     

      dstatic[k][m].usflag=1;

      if((CONTERMINOUS == 1) && (country_code == 221)){


	if(lat[k] > 50.) dstatic[k][m].usflag=0;
	if((lat[k] < 25.) && (lon[m] > -60.)) dstatic[k][m].usflag=0;
	if(lon[m] < -130.) dstatic[k][m].usflag=0;

	  
      }

    }
  }



  for (it=0;it<NEWNT;it++){

    byclass_sum[it].c=ZEROVALUE;
    byclass_sum[it].p=ZEROVALUE;
    byclass_sum[it].v=ZEROVALUE;
    byclass_sum[it].i=ZEROVALUE;
    byclass_sum[it].w=ZEROVALUE;
    byclass_sum[it].s=ZEROVALUE;
    byclass_sum[it].u=ZEROVALUE;
    byclass_sum[it].total=ZEROVALUE;
    byclass_sum[it].vf=ZEROVALUE;
    byclass_sum[it].vnf=ZEROVALUE;
    byclass_sum[it].sf=ZEROVALUE;
    byclass_sum[it].snf=ZEROVALUE;
    byclass_sum[it].smaf=ZEROVALUE;
    byclass_sum[it].smanf=ZEROVALUE;
    byclass_sum[it].vbh=ZEROVALUE;
    byclass_sum[it].sbh=ZEROVALUE;
    byclass_sum[it].vbh2=ZEROVALUE;
    byclass_sum[it].sbh2=ZEROVALUE;
    byclass_sum[it].sbh3=ZEROVALUE;
    byclass_sum[it].wh_unmet=ZEROVALUE;
    byclass_sum[it].wh=ZEROVALUE;
    byclass_sum[it].cavg=ZEROVALUE;
    byclass_sum[it].pavg=ZEROVALUE;
    byclass_sum[it].cpavg=ZEROVALUE;
    byclass_sum[it].flowcp=ZEROVALUE;
    byclass_sum[it].flowpc=ZEROVALUE;
    byclass_sum[it].flowpv=ZEROVALUE;
    byclass_sum[it].flowvp=ZEROVALUE;
    byclass_sum[it].flowvc=ZEROVALUE;
    byclass_sum[it].flowcv=ZEROVALUE;
    byclass_sum[it].flowsp=ZEROVALUE;
    byclass_sum[it].flowps=ZEROVALUE;
    byclass_sum[it].flowsc=ZEROVALUE;
    byclass_sum[it].flowcs=ZEROVALUE;
    byclass_sum[it].flowcu=ZEROVALUE;
    byclass_sum[it].flowpu=ZEROVALUE;
    byclass_sum[it].flowvu=ZEROVALUE;
    byclass_sum[it].flowsu=ZEROVALUE;
    byclass_sum[it].flowuc=ZEROVALUE;
    byclass_sum[it].flowup=ZEROVALUE;
    byclass_sum[it].flowus=ZEROVALUE;
    byclass_sum[it].flowsbh=ZEROVALUE;
    byclass_sum[it].flowsbh2=ZEROVALUE;
    byclass_sum[it].flowsbh3=ZEROVALUE;
    byclass_sum[it].flowvbh=ZEROVALUE;
    byclass_sum[it].flowvbh2=ZEROVALUE;
    byclass_sum[it].converted_forest_land=ZEROVALUE;
    byclass_sum[it].converted_forest_land_total=ZEROVALUE;
    byclass_sum[it].flowcs_prime=ZEROVALUE;
    byclass_sum[it].flowvc_prime=ZEROVALUE;
    byclass_sum[it].flowsc_prime=ZEROVALUE;
    byclass_sum[it].flowvp_prime=ZEROVALUE;
    byclass_sum[it].flowsp_prime=ZEROVALUE;
    byclass_sum[it].flowps_prime=ZEROVALUE;
    byclass_sum[it].flowsu_prime=ZEROVALUE;
    byclass_sum[it].flowvu_prime=ZEROVALUE;
    byclass_sum[it].flowus_prime=ZEROVALUE;
  }

  return;
}
/********************************************************************/

void global_timeseries_checker(int country_code, char country_name[50]){

  FILE *outfile, *outfile2, *outfile3;
  int i, ic, k, m, it;
  double garea_sum=ZEROVALUE, land_garea_sum=ZEROVALUE, total_garea_sum=ZEROVALUE;
  double nf_s_area=ZEROVALUE, fnf_s_area=ZEROVALUE, check_total2[NEWNT];
  char outstat[2], ffname[90];
 
  printf("printing global tester...\n"); 

  initialize_checker(country_code);
 

  strcpy(ffname,country_name);  
  strcat(ffname,".test");



  if(strcmp("tone",trun) == 0) {
    outfile=fopen(ffname,"w");
    outfile2=fopen("global.primeflow.txt","w");
    outfile3=fopen("check_diff.txt","w");
  }
  else{
    outfile=fopen(ffname,"a");
    outfile2=fopen("global.primeflow.txt","a");
    outfile3=fopen("check_diff.txt","a");
  }


   for (it=0;it<NEWNT-1;it++){
     check_total2[it]=0;
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){   
	
#if COUNTRY_TEST
	
	if((dstatic[k][m].gcode==country_code) && (dstatic[k][m].usflag==1)){   	
#else
	   if(dstatic[k][m].gcode > 0){ 
#endif		
	   
	    if(it == 0) {
	      chtest[k][m][it].c = newdata[k][m][it].c;
	      chtest[k][m][it].p = newdata[k][m][it].p;
	      chtest[k][m][it].v = newdata[k][m][it].v;
	      chtest[k][m][it].s = newdata[k][m][it].s;
	      chtest[k][m][it].u = newdata[k][m][it].u;
	    }
	    
	    chtest[k][m][it+1].c=chtest[k][m][it].c +
	      newdata[k][m][it].flowpc + newdata[k][m][it].flowvc + 
	      newdata[k][m][it].flowsc - 
              newdata[k][m][it].flowcp - newdata[k][m][it].flowcv -
              newdata[k][m][it].flowcs - newdata[k][m][it].flowcu + newdata[k][m][it].flowuc; 
	    
	    chtest[k][m][it+1].p=chtest[k][m][it].p +
	      newdata[k][m][it].flowcp + newdata[k][m][it].flowvp + 
	      newdata[k][m][it].flowsp - 
              newdata[k][m][it].flowpc - newdata[k][m][it].flowpv - 
              newdata[k][m][it].flowps - newdata[k][m][it].flowpu + newdata[k][m][it].flowup;
	    
	    chtest[k][m][it+1].v=chtest[k][m][it].v +
	      newdata[k][m][it].flowpv + newdata[k][m][it].flowcv -
	      newdata[k][m][it].flowvp - newdata[k][m][it].flowvc -   
	      newdata[k][m][it].flowvs - newdata[k][m][it].flowvu;
	    
	    chtest[k][m][it+1].s=chtest[k][m][it].s +
	      newdata[k][m][it].flowps + newdata[k][m][it].flowvs + 
	      newdata[k][m][it].flowcs -
	      newdata[k][m][it].flowsp - newdata[k][m][it].flowsc -
              newdata[k][m][it].flowsu + newdata[k][m][it].flowus;

	    chtest[k][m][it+1].u=chtest[k][m][it].u +
	      newdata[k][m][it].flowcu + newdata[k][m][it].flowvu + 
              newdata[k][m][it].flowsu + newdata[k][m][it].flowpu - 
	      newdata[k][m][it].flowuc - newdata[k][m][it].flowus -
	      newdata[k][m][it].flowup;  

	  

check_total2[it]+=(newdata[k][m][it].v+newdata[k][m][it].s+newdata[k][m][it].c+newdata[k][m][it].p+newdata[k][m][it].u+newdata[k][m][it].i+newdata[k][m][it].w)*garea[k][m];
	    


#if 1
	  } /* end of country_test */
#endif

	} /* end of m */
      } /* end of m */
    }  /* end of it */





  for (it=0;it<NEWNT-1;it++){
    

    for (i=0;i<NCCODE;i++) {

#if COUNTRY_TEST
	
      if(cdata[i].ccode == country_code) {
	byclass_sum[it].wh=ctdata[i][it].wh;
	byclass_sum[it].converted_forest_land=ctdata[i][it].converted_forest_land;
	byclass_sum[it].converted_forest_land_total=ctdata[i][it].converted_forest_land_total;
	byclass_sum[it].flowcs_prime+=ctdata[i][it].flowcs_prime;
	byclass_sum[it].flowvc_prime+=ctdata[i][it].flowvc_prime;
	byclass_sum[it].flowsc_prime+=ctdata[i][it].flowsc_prime;
	byclass_sum[it].flowvp_prime+=ctdata[i][it].flowvp_prime;
	byclass_sum[it].flowsp_prime+=ctdata[i][it].flowsp_prime;
	byclass_sum[it].flowps_prime+=ctdata[i][it].flowps_prime;
	byclass_sum[it].flowus_prime+=ctdata[i][it].flowus_prime;
	byclass_sum[it].flowsu_prime+=ctdata[i][it].flowsu_prime;
	byclass_sum[it].flowvu_prime+=ctdata[i][it].flowvu_prime;
      }
#else
      byclass_sum[it].converted_forest_land+=ctdata[i][it].converted_forest_land;
      byclass_sum[it].converted_forest_land_total+=ctdata[i][it].converted_forest_land_total;
      byclass_sum[it].wh+=ctdata[i][it].wh;

      if(cdata[i].ccode > 0){
    
	byclass_sum[it].flowcs_prime+=ctdata[i][it].flowcs_prime;
	byclass_sum[it].flowvc_prime+=ctdata[i][it].flowvc_prime;
	byclass_sum[it].flowsc_prime+=ctdata[i][it].flowsc_prime;
	byclass_sum[it].flowvp_prime+=ctdata[i][it].flowvp_prime;
	byclass_sum[it].flowsp_prime+=ctdata[i][it].flowsp_prime;
	byclass_sum[it].flowps_prime+=ctdata[i][it].flowps_prime;
	byclass_sum[it].flowus_prime+=ctdata[i][it].flowus_prime;
	byclass_sum[it].flowsu_prime+=ctdata[i][it].flowsu_prime;
	byclass_sum[it].flowvu_prime+=ctdata[i][it].flowvu_prime;
      }
#endif
   
    } /* end of i */

 

    fnf_s_area=ZEROVALUE;
    nf_s_area=ZEROVALUE;
    
    

    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){		

	i=dstatic[k][m].newgcode;

#if COUNTRY_TEST
	
	if((dstatic[k][m].gcode==country_code) && (dstatic[k][m].usflag==1)){   	
#else
	    if(dstatic[k][m].gcode > 0){ 
#endif		

	     
          
	    if(dstatic[k][m].fnf == 1){
		
	      byclass_sum[it].vf+=garea[k][m]*chtest[k][m][it].v;
	      byclass_sum[it].sf+=garea[k][m]*chtest[k][m][it].s;
	      byclass_sum[it].smaf+=chtest[k][m][it].s*garea[k][m]*newdata[k][m][it].sma; 
	      fnf_s_area += chtest[k][m][it].s*garea[k][m];
	    }
	    else{
	      byclass_sum[it].vnf+=garea[k][m]*chtest[k][m][it].v;
	      byclass_sum[it].snf+=garea[k][m]*chtest[k][m][it].s;
	      byclass_sum[it].smanf+=chtest[k][m][it].s*garea[k][m]*newdata[k][m][it].sma;
	      nf_s_area += chtest[k][m][it].s*garea[k][m];
	    }
	    
	    byclass_sum[it].vbh+=newdata[k][m][it].vbh;
	    byclass_sum[it].sbh+=newdata[k][m][it].sbh;
	    byclass_sum[it].vbh2+=newdata[k][m][it].vbh2;
	    byclass_sum[it].sbh2+=newdata[k][m][it].sbh2;
	    byclass_sum[it].sbh3+=newdata[k][m][it].sbh3;
	    byclass_sum[it].wh_unmet+=ZEROVALUE; 
	    
	    
	    byclass_sum[it].c+=garea[k][m]*chtest[k][m][it].c; 
	    byclass_sum[it].p+=garea[k][m]*chtest[k][m][it].p;
	    byclass_sum[it].v+=garea[k][m]*chtest[k][m][it].v; 
	    byclass_sum[it].i+=garea[k][m]*newdata[k][m][it].i;  
	    byclass_sum[it].w+=garea[k][m]*newdata[k][m][it].w;  
	    byclass_sum[it].s+=garea[k][m]*chtest[k][m][it].s; 
	    byclass_sum[it].u+=garea[k][m]*chtest[k][m][it].u; 

	    
	    byclass_sum[it].flowcp+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcp;
	    byclass_sum[it].flowpc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowpc;
	    byclass_sum[it].flowpv+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowpv;
	    byclass_sum[it].flowvp+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowvp;
	    byclass_sum[it].flowvc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowvc;
	    byclass_sum[it].flowcv+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcv;
	    byclass_sum[it].flowsp+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowsp;
	    byclass_sum[it].flowps+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowps;
	    byclass_sum[it].flowsc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowsc;
	    byclass_sum[it].flowcs+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcs;
	    byclass_sum[it].flowcu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcu;
	    byclass_sum[it].flowpu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowpu;
	    byclass_sum[it].flowvu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowvu;
	    byclass_sum[it].flowsu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowsu;
	    byclass_sum[it].flowuc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowuc;
	    byclass_sum[it].flowup+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowup;
	    byclass_sum[it].flowus+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowus;


	    if(dstatic[k][m].fnf == 1){

	      byclass_sum[it].flowvbh+=newdata[k][m][it].vbh/dstatic[k][m].vba/1000.0/1000.0;
	      
	      if(newdata[k][m][it].smb > ZEROVALUE){
		
		byclass_sum[it].flowsbh+=newdata[k][m][it].sbh/newdata[k][m][it].smb/1000.0/1000.0;
		byclass_sum[it].flowsbh2+=newdata[k][m][it].sbh2/newdata[k][m][it].smb/1000.0/1000.0;
	      }
	      
	    }
	    else{
	      if(dstatic[k][m].vba > ZEROVALUE) byclass_sum[it].flowvbh2+=newdata[k][m][it].vbh2/dstatic[k][m].vba/1000.0/1000.0; 
	      if(newdata[k][m][it].smb > ZEROVALUE) byclass_sum[it].flowsbh3+=newdata[k][m][it].sbh3/newdata[k][m][it].smb/1000.0/1000.0;
	    }
	    
	    
	    if(it == 0 ) land_garea_sum+=(chtest[k][m][it].v+chtest[k][m][it].s+chtest[k][m][it].p+chtest[k][m][it].c+chtest[k][m][it].u)*garea[k][m]; 
	    
	    if(it == 0) garea_sum+=garea[k][m];
	      
#if 1
	}
#endif

	  if(it==0) total_garea_sum+=garea[k][m];

	  
	} /* end of m*/
      } /* end of k*/
    
     
      byclass_sum[it].total=byclass_sum[it].c+byclass_sum[it].p+
	byclass_sum[it].v+byclass_sum[it].i+
	byclass_sum[it].w+byclass_sum[it].s+byclass_sum[it].u;  
      

      if(fnf_s_area > ZEROVALUE) byclass_sum[it].smaf=byclass_sum[it].smaf/fnf_s_area;
    
      if(nf_s_area > ZEROVALUE) byclass_sum[it].smanf=byclass_sum[it].smanf/nf_s_area;
    
    

      if((it == 0) && (strcmp("tone",trun) == 0)) fprintf(outfile,"yr c p v i w s u fsum vf vnf sf snf smaf smanf vbh sbh vbh2 sbh2 sbh3 sum wh clearing_amount_total clearing_amount_used wh-sum-clearing_amount_ifused(unmet) flowcp flowpc flowpv flowvp flowvc flowcv flowsp flowps flowsc flowcs flowcu flowpu flowvu flowsu flowuc flowup flowus flowvbh flowsbh flowvbh2 flowsbh2 flowsbh3 land_asum_km2 land_iw_asum_km2 total_asum_km2\n");
      
      fprintf(outfile,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	      year[it],
	      byclass_sum[it].c/garea_sum,
	      byclass_sum[it].p/garea_sum,
	      byclass_sum[it].v/garea_sum,
	      byclass_sum[it].i/garea_sum,
	      byclass_sum[it].w/garea_sum,
	      byclass_sum[it].s/garea_sum,
	      byclass_sum[it].u/garea_sum,
	      byclass_sum[it].total/garea_sum,
	      byclass_sum[it].vf/garea_sum,
	      byclass_sum[it].vnf/garea_sum,
	      byclass_sum[it].sf/garea_sum,
	      byclass_sum[it].snf/garea_sum,
	      byclass_sum[it].smaf,
	      byclass_sum[it].smanf,
	      byclass_sum[it].vbh/1000.,
	      byclass_sum[it].sbh/1000.,
	      byclass_sum[it].vbh2/1000.,
	      byclass_sum[it].sbh2/1000.,
	      byclass_sum[it].sbh3/1000.,
	      (byclass_sum[it].vbh+byclass_sum[it].sbh+byclass_sum[it].vbh2+byclass_sum[it].sbh2+byclass_sum[it].sbh3)/1000.,
              byclass_sum[it].wh,
	      byclass_sum[it].converted_forest_land_total,
	      byclass_sum[it].converted_forest_land,
	      byclass_sum[it].wh-((byclass_sum[it].vbh+byclass_sum[it].sbh+byclass_sum[it].vbh2+byclass_sum[it].sbh2+byclass_sum[it].sbh3)/1000.+byclass_sum[it].converted_forest_land),
	      byclass_sum[it].flowcp,
	      byclass_sum[it].flowpc,
	      byclass_sum[it].flowpv,
	      byclass_sum[it].flowvp,
	      byclass_sum[it].flowvc,
	      byclass_sum[it].flowcv,
	      byclass_sum[it].flowsp,
	      byclass_sum[it].flowps,
	      byclass_sum[it].flowsc,
	      byclass_sum[it].flowcs,
	      byclass_sum[it].flowcu,
	      byclass_sum[it].flowpu,
	      byclass_sum[it].flowvu,
	      byclass_sum[it].flowsu,
	      byclass_sum[it].flowuc,
	      byclass_sum[it].flowup,
	      byclass_sum[it].flowus,
	      byclass_sum[it].flowvbh,
	      byclass_sum[it].flowsbh,
	      byclass_sum[it].flowvbh2,
	      byclass_sum[it].flowsbh2,
	      byclass_sum[it].flowsbh3,
	      land_garea_sum/1000.0/1000.0,
	      garea_sum/1000.0/1000.0,
	      total_garea_sum/1000.0/1000.0); 


      if((it == 0) && (strcmp("tone",trun) == 0)) fprintf(outfile2,"yr cs vc sc vp sp ps vu su us\n");
      
      fprintf(outfile2,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	      year[it],
	      byclass_sum[it].flowcs_prime,
	      byclass_sum[it].flowvc_prime,
	      byclass_sum[it].flowsc_prime,
	      byclass_sum[it].flowvp_prime,
	      byclass_sum[it].flowsp_prime,
	      byclass_sum[it].flowps_prime,
	      byclass_sum[it].flowvu_prime,
	      byclass_sum[it].flowsu_prime,
	      byclass_sum[it].flowus_prime);

           fprintf(outfile3,"diffrence in totals: %lf\n",check_total2[it]-check_total[it]); 

    } /* end of it */

    fclose(outfile);
    fclose(outfile2);
    fclose(outfile3);


    return;

     
  }

/********************************************************************/

void country_timeseries_checker(int country_code, char country_name[50]){

  FILE *outfile;
  int i, ic, k, m, it;
  double garea_sum=ZEROVALUE, land_garea_sum=ZEROVALUE, total_garea_sum=ZEROVALUE;
  double nf_s_area=ZEROVALUE, fnf_s_area=ZEROVALUE;
  char outstat[2], ffname[90];
  double tester1, tester2, tester3, tester4;



  printf("printing country tester...\n");  

  initialize_checker(country_code);

  printf("code %d\n",country_code);


  strcpy(ffname,country_name);  
  strcat(ffname,".country.test");


  if(strcmp("tone",trun) == 0) 
    outfile=fopen(ffname,"w");
  else
    outfile=fopen(ffname,"a");


  for (it=0;it<NEWNT-1;it++){
    
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){   
	
	
	if((dstatic[k][m].gcode==country_code) && (dstatic[k][m].usflag==1)){   	
	  
	  if(it == 0) {
	    chtest[k][m][it].c = newdata[k][m][it].c;
	    chtest[k][m][it].p = newdata[k][m][it].p;
	    chtest[k][m][it].v = newdata[k][m][it].v;
	    chtest[k][m][it].s = newdata[k][m][it].s;
	    chtest[k][m][it].u = newdata[k][m][it].u;
	  }
	  
	  chtest[k][m][it+1].c=chtest[k][m][it].c +
	    newdata[k][m][it].flowpc + newdata[k][m][it].flowvc + 
	    newdata[k][m][it].flowsc - newdata[k][m][it].flowcp - 
	    newdata[k][m][it].flowcv - newdata[k][m][it].flowcs +
	    newdata[k][m][it].flowuc - newdata[k][m][it].flowcu; 
	  
	  chtest[k][m][it+1].p=chtest[k][m][it].p +
	    newdata[k][m][it].flowcp + newdata[k][m][it].flowvp + 
	    newdata[k][m][it].flowsp - newdata[k][m][it].flowpc - 
	    newdata[k][m][it].flowpv - newdata[k][m][it].flowps +
	    newdata[k][m][it].flowup - newdata[k][m][it].flowpu;
	    
	  chtest[k][m][it+1].v=chtest[k][m][it].v +
	    newdata[k][m][it].flowpv + newdata[k][m][it].flowcv -
	    newdata[k][m][it].flowvp - newdata[k][m][it].flowvc -   
	    newdata[k][m][it].flowvs - newdata[k][m][it].flowvu;
	    
	  chtest[k][m][it+1].s=chtest[k][m][it].s +
	    newdata[k][m][it].flowvs + newdata[k][m][it].flowps + 
	    newdata[k][m][it].flowcs -
	    newdata[k][m][it].flowsp - newdata[k][m][it].flowsc +
	    newdata[k][m][it].flowus - newdata[k][m][it].flowsu;  

	  chtest[k][m][it+1].u=chtest[k][m][it].u +
	    newdata[k][m][it].flowcu + newdata[k][m][it].flowvu + 
	    newdata[k][m][it].flowsu + newdata[k][m][it].flowpu - 
	    newdata[k][m][it].flowuc - newdata[k][m][it].flowup -
	    newdata[k][m][it].flowus;
	  
	 
	} /* end of country_test */
      } /* end of m */
    } /* end of m */
  }  /* end of it */
  
 

  for (it=0;it<NEWNT-1;it++){
    

    for (i=0;i<NCCODE;i++) {

      if(cdata[i].ccode == country_code) {
	byclass_sum[it].wh=ctdata[i][it].wh;
	byclass_sum[it].converted_forest_land=ctdata[i][it].converted_forest_land;
	byclass_sum[it].converted_forest_land_total=ctdata[i][it].converted_forest_land_total;
     
      }
      
    } /* end of i */

 
    fnf_s_area=ZEROVALUE;
    nf_s_area=ZEROVALUE;
    
    

    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){		

	if((dstatic[k][m].gcode==country_code) && (dstatic[k][m].usflag==1)){   	

	
	  if(dstatic[k][m].fnf == 1){
		
	    byclass_sum[it].vf+=garea[k][m]*chtest[k][m][it].v;
	    byclass_sum[it].sf+=garea[k][m]*chtest[k][m][it].s;
	    byclass_sum[it].smaf+=chtest[k][m][it].s*garea[k][m]*newdata[k][m][it].sma;
	    fnf_s_area += chtest[k][m][it].s*garea[k][m];
	  }
	  else{
	    byclass_sum[it].vnf+=garea[k][m]*chtest[k][m][it].v;
	    byclass_sum[it].snf+=garea[k][m]*chtest[k][m][it].s;
	    byclass_sum[it].smanf+=chtest[k][m][it].s*garea[k][m]*newdata[k][m][it].sma;
	    nf_s_area += chtest[k][m][it].s*garea[k][m];
	  }
	    
	  byclass_sum[it].vbh+=newdata[k][m][it].vbh/1000.;
	  byclass_sum[it].sbh+=newdata[k][m][it].sbh/1000.;
	  byclass_sum[it].vbh2+=newdata[k][m][it].vbh2/1000.;
	  byclass_sum[it].sbh2+=newdata[k][m][it].sbh2/1000.;
	  byclass_sum[it].sbh3+=newdata[k][m][it].sbh3/1000.;
	  byclass_sum[it].wh_unmet+=ZEROVALUE;  
	  
	  
	  byclass_sum[it].c+=garea[k][m]*chtest[k][m][it].c; 
	  byclass_sum[it].p+=garea[k][m]*chtest[k][m][it].p;
	  byclass_sum[it].v+=garea[k][m]*chtest[k][m][it].v; 
	  byclass_sum[it].i+=garea[k][m]*newdata[k][m][it].i;  
	  byclass_sum[it].w+=garea[k][m]*newdata[k][m][it].w;  
	  byclass_sum[it].s+=garea[k][m]*chtest[k][m][it].s; 
	  byclass_sum[it].u+=garea[k][m]*chtest[k][m][it].u; 
	  
	  byclass_sum[it].flowcp+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcp;
	  byclass_sum[it].flowpc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowpc;
	  byclass_sum[it].flowpv+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowpv;
	  byclass_sum[it].flowvp+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowvp;
	  byclass_sum[it].flowvc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowvc;
	  byclass_sum[it].flowcv+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcv;
	  byclass_sum[it].flowsp+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowsp;
	  byclass_sum[it].flowps+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowps;
	  byclass_sum[it].flowsc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowsc;
	  byclass_sum[it].flowcs+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcs;
	  byclass_sum[it].flowcu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcu;
	  byclass_sum[it].flowpu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowpu;
	  byclass_sum[it].flowvu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowvu;
	  byclass_sum[it].flowsu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowsu;
	  byclass_sum[it].flowuc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowuc;
	  byclass_sum[it].flowup+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowup;
	  byclass_sum[it].flowus+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowus;

          
	  
	  if(dstatic[k][m].fnf == 1){

	    if(dstatic[k][m].vba > ZEROVALUE) byclass_sum[it].flowvbh+=newdata[k][m][it].vbh/dstatic[k][m].vba/1000.0/1000.0;

	    if(newdata[k][m][it].smb > ZEROVALUE){
	     
	      if(newdata[k][m][it].smb > ZEROVALUE) byclass_sum[it].flowsbh+=newdata[k][m][it].sbh/newdata[k][m][it].smb/1000.0/1000.0;
	      if(newdata[k][m][it].smb > ZEROVALUE) byclass_sum[it].flowsbh2+=newdata[k][m][it].sbh2/newdata[k][m][it].smb/1000.0/1000.0;

	      if(newdata[k][m][it].smb > ZEROVALUE) byclass_sum[it].flowsbh3+=newdata[k][m][it].sbh3/newdata[k][m][it].smb/1000.0/1000.0;



	    }
	    
	  }
	  else{
	    if(dstatic[k][m].vba > ZEROVALUE) byclass_sum[it].flowvbh2+=newdata[k][m][it].vbh2/dstatic[k][m].vba/1000.0/1000.0;
	  } 
	  	    
	  
	  
	  if(it == 0 ) land_garea_sum+=(chtest[k][m][it].v+chtest[k][m][it].s+chtest[k][m][it].p+chtest[k][m][it].c+chtest[k][m][it].u)*garea[k][m]; 
	  
	  if(it == 0) garea_sum+=garea[k][m];
	      

	} /* end of country_test */

        if(it==0) total_garea_sum+=garea[k][m];


      } /* end of m*/
    } /* end of k*/
    
  
    byclass_sum[it].total+=byclass_sum[it].c+byclass_sum[it].p+
      byclass_sum[it].v+byclass_sum[it].i+
      byclass_sum[it].w+byclass_sum[it].s+byclass_sum[it].u;  
      
      
    if(fnf_s_area > ZEROVALUE) byclass_sum[it].smaf=byclass_sum[it].smaf/fnf_s_area;
    
    if(nf_s_area > ZEROVALUE) byclass_sum[it].smanf=byclass_sum[it].smanf/nf_s_area;
    


    if((it == 0) && (strcmp("tone",trun) == 0)) fprintf(outfile,"yr c p v i w s u fsum vf vnf sf snf smaf smanf vbh sbh vbh2 sbh2 sbh3 sum wh clearing_amount_total clearing_amount_used wh-sum-clearing_amount_ifused(unmet) flowcp flowpc flowpv flowvp flowvc flowcv flowsp flowps flowsc flowcs flowcu flowpu flowvu flowsu flowuc flowup flowus flowvbh flowsbh flowvbh2 flowsbh2 flowsbh3 land_asum_km2 land_iw_asum_km2 total_garea_sum_km2\n");
      
    fprintf(outfile,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	    year[it],
	    byclass_sum[it].c/garea_sum,
	    byclass_sum[it].p/garea_sum,
	    byclass_sum[it].v/garea_sum,
	    byclass_sum[it].i/garea_sum,
	    byclass_sum[it].w/garea_sum,
	    byclass_sum[it].s/garea_sum,
	    byclass_sum[it].u/garea_sum,
	    byclass_sum[it].total/garea_sum,
	    byclass_sum[it].vf/garea_sum,
	    byclass_sum[it].vnf/garea_sum,
	    byclass_sum[it].sf/garea_sum,
	    byclass_sum[it].snf/garea_sum,
	    byclass_sum[it].smaf,
	    byclass_sum[it].smanf,
	    byclass_sum[it].vbh,
	    byclass_sum[it].sbh,
	    byclass_sum[it].vbh2,
	    byclass_sum[it].sbh2,
	    byclass_sum[it].sbh3,
	    byclass_sum[it].vbh+byclass_sum[it].sbh+byclass_sum[it].vbh2+byclass_sum[it].sbh2+byclass_sum[it].sbh3,
	    byclass_sum[it].wh,
	    byclass_sum[it].converted_forest_land_total,
	    byclass_sum[it].converted_forest_land,
	    byclass_sum[it].wh-(byclass_sum[it].vbh+byclass_sum[it].sbh+byclass_sum[it].vbh2+byclass_sum[it].sbh2+byclass_sum[it].sbh3+byclass_sum[it].converted_forest_land),
	    byclass_sum[it].flowcp,
	    byclass_sum[it].flowpc,
   	    byclass_sum[it].flowpv,
	    byclass_sum[it].flowvp,
	    byclass_sum[it].flowvc,
	    byclass_sum[it].flowcv,
	    byclass_sum[it].flowsp,
	    byclass_sum[it].flowps,
	    byclass_sum[it].flowsc,
	    byclass_sum[it].flowcs,
	    byclass_sum[it].flowcu,
	    byclass_sum[it].flowpu,
	    byclass_sum[it].flowvu,
	    byclass_sum[it].flowsu,
	    byclass_sum[it].flowuc,
	    byclass_sum[it].flowup,
	    byclass_sum[it].flowus,
	    byclass_sum[it].flowvbh,
	    byclass_sum[it].flowsbh,
	    byclass_sum[it].flowvbh2,
	    byclass_sum[it].flowsbh2,
	    byclass_sum[it].flowsbh3,
	    land_garea_sum/1000.0/1000.0,
	    garea_sum/1000.0/1000.0,
	    total_garea_sum/1000.0/1000.0 ); 
      
  } /* end of it */

  fclose(outfile);

  

  return;
  
}

/********************************************************************/

void continent_timeseries_checker(int continent_code, char continent_name[50]){

  FILE *outfile, *outfile2;
  int i, ic, k, m, it;
  double garea_sum=ZEROVALUE, land_garea_sum=ZEROVALUE, total_garea_sum=ZEROVALUE;
  double nf_s_area=ZEROVALUE, fnf_s_area=ZEROVALUE;
  char outstat[2], ffname[90], ffname2[90];
  double tester1, tester2, tester3, tester4;


 
  printf("printing continent tester...\n");



  initialize_checker(continent_code);


  strcpy(ffname,continent_name);  
  strcat(ffname,".continent.test");

  strcpy(ffname2,continent_name);  
  strcat(ffname2,".primeflow.txt");




  if(strcmp("tone",trun) == 0){ 
    outfile=fopen(ffname,"w");
    outfile2=fopen(ffname2,"w");
  }
  else{
    outfile=fopen(ffname,"a");
    outfile2=fopen(ffname2,"a");
  }
    


  for (it=0;it<NEWNT-1;it++){
    
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){   
	
	
	if(dstatic[k][m].gcontinent_code==continent_code){   	
	  
	  if(it == 0) {
	    chtest[k][m][it].c = newdata[k][m][it].c;
	    chtest[k][m][it].p = newdata[k][m][it].p;
	    chtest[k][m][it].v = newdata[k][m][it].v;
	    chtest[k][m][it].s = newdata[k][m][it].s;
	    chtest[k][m][it].u = newdata[k][m][it].u;
	  }
	  
	  chtest[k][m][it+1].c=chtest[k][m][it].c +
	    newdata[k][m][it].flowpc + newdata[k][m][it].flowvc + 
	    newdata[k][m][it].flowsc - newdata[k][m][it].flowcp - 
	    newdata[k][m][it].flowcv - newdata[k][m][it].flowcs +
	    newdata[k][m][it].flowuc - newdata[k][m][it].flowcu; 
	  
	  chtest[k][m][it+1].p=chtest[k][m][it].p +
	    newdata[k][m][it].flowcp + newdata[k][m][it].flowvp + 
	    newdata[k][m][it].flowsp - newdata[k][m][it].flowpc - 
	    newdata[k][m][it].flowpv - newdata[k][m][it].flowps +
	    newdata[k][m][it].flowup - newdata[k][m][it].flowpu;
	    
	  chtest[k][m][it+1].v=chtest[k][m][it].v +
	    newdata[k][m][it].flowpv + newdata[k][m][it].flowcv -
	    newdata[k][m][it].flowvp - newdata[k][m][it].flowvc -   
	    newdata[k][m][it].flowvs - newdata[k][m][it].flowvu;
	    
	  chtest[k][m][it+1].s=chtest[k][m][it].s +
	    newdata[k][m][it].flowvs + newdata[k][m][it].flowps + 
	    newdata[k][m][it].flowcs -
	    newdata[k][m][it].flowsp - newdata[k][m][it].flowsc +
	    newdata[k][m][it].flowus - newdata[k][m][it].flowsu; 

	  chtest[k][m][it+1].u=chtest[k][m][it].u +
	    newdata[k][m][it].flowcu + newdata[k][m][it].flowpu + 
	    newdata[k][m][it].flowsu + newdata[k][m][it].flowvu - 
	    newdata[k][m][it].flowuc - newdata[k][m][it].flowup -
	    newdata[k][m][it].flowus;
 
	  
	} /* end of continent_test */
      } /* end of m */
    } /* end of m */
  }  /* end of it */
  

  for (it=0;it<NEWNT-1;it++){
    

    for (i=0;i<NCCODE;i++) {

      if(cdata[i].continent_code == continent_code) {
	byclass_sum[it].wh+=ctdata[i][it].wh;
	byclass_sum[it].converted_forest_land+=ctdata[i][it].converted_forest_land;
	byclass_sum[it].converted_forest_land_total+=ctdata[i][it].converted_forest_land_total;

	byclass_sum[it].flowcs_prime+=ctdata[i][it].flowcs_prime;
	byclass_sum[it].flowvc_prime+=ctdata[i][it].flowvc_prime;
	byclass_sum[it].flowsc_prime+=ctdata[i][it].flowsc_prime;
	byclass_sum[it].flowvp_prime+=ctdata[i][it].flowvp_prime;
	byclass_sum[it].flowsp_prime+=ctdata[i][it].flowsp_prime;
	byclass_sum[it].flowps_prime+=ctdata[i][it].flowps_prime;
	byclass_sum[it].flowus_prime+=ctdata[i][it].flowus_prime;
	byclass_sum[it].flowsu_prime+=ctdata[i][it].flowsu_prime;
	byclass_sum[it].flowvu_prime+=ctdata[i][it].flowvu_prime;
      }
      
    } /* end of i */

 
    fnf_s_area=ZEROVALUE;
    nf_s_area=ZEROVALUE;
    
    

    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){		

	if(dstatic[k][m].gcontinent_code==continent_code){   	

	
	  if(dstatic[k][m].fnf == 1){
		
	    byclass_sum[it].vf+=garea[k][m]*chtest[k][m][it].v;
	    byclass_sum[it].sf+=garea[k][m]*chtest[k][m][it].s;
	    byclass_sum[it].smaf+=chtest[k][m][it].s*garea[k][m]*newdata[k][m][it].sma;
	    fnf_s_area += chtest[k][m][it].s*garea[k][m];
	  }
	  else{
	    byclass_sum[it].vnf+=garea[k][m]*chtest[k][m][it].v;
	    byclass_sum[it].snf+=garea[k][m]*chtest[k][m][it].s;
	    byclass_sum[it].smanf+=chtest[k][m][it].s*garea[k][m]*newdata[k][m][it].sma;
	    nf_s_area += chtest[k][m][it].s*garea[k][m];
	  }
	    
	  byclass_sum[it].vbh+=newdata[k][m][it].vbh/1000.;
	  byclass_sum[it].sbh+=newdata[k][m][it].sbh/1000.;
	  byclass_sum[it].vbh2+=newdata[k][m][it].vbh2/1000.;
	  byclass_sum[it].sbh2+=newdata[k][m][it].sbh2/1000.;
	  byclass_sum[it].sbh3+=newdata[k][m][it].sbh3/1000.;
	  byclass_sum[it].wh_unmet+=ZEROVALUE;  
	  
	  
	  byclass_sum[it].c+=garea[k][m]*chtest[k][m][it].c; 
	  byclass_sum[it].p+=garea[k][m]*chtest[k][m][it].p;
	  byclass_sum[it].v+=garea[k][m]*chtest[k][m][it].v; 
	  byclass_sum[it].i+=garea[k][m]*newdata[k][m][it].i;  
	  byclass_sum[it].w+=garea[k][m]*newdata[k][m][it].w;  
	  byclass_sum[it].s+=garea[k][m]*chtest[k][m][it].s; 
	  byclass_sum[it].u+=garea[k][m]*chtest[k][m][it].u; 
	  
	  byclass_sum[it].flowcp+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcp;
	  byclass_sum[it].flowpc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowpc;
	  byclass_sum[it].flowpv+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowpv;
	  byclass_sum[it].flowvp+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowvp;
	  byclass_sum[it].flowvc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowvc;
	  byclass_sum[it].flowcv+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcv;
	  byclass_sum[it].flowsp+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowsp;
	  byclass_sum[it].flowps+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowps;
	  byclass_sum[it].flowsc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowsc;
	  byclass_sum[it].flowcs+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcs;
	  byclass_sum[it].flowcu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcu;
	  byclass_sum[it].flowpu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowpu;
	  byclass_sum[it].flowvu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowvu;
	  byclass_sum[it].flowsu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowsu;
	  byclass_sum[it].flowuc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowuc;
	  byclass_sum[it].flowup+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowup;
	  byclass_sum[it].flowus+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowus;

          
	  
	  if(dstatic[k][m].fnf == 1){

	    byclass_sum[it].flowvbh+=newdata[k][m][it].vbh/dstatic[k][m].vba/1000.0/1000.0;

	    if(newdata[k][m][it].smb > ZEROVALUE){
	     
	      byclass_sum[it].flowsbh+=newdata[k][m][it].sbh/newdata[k][m][it].smb/1000.0/1000.0;
	      byclass_sum[it].flowsbh2+=newdata[k][m][it].sbh2/newdata[k][m][it].smb/1000.0/1000.0;
	      byclass_sum[it].flowsbh3+=newdata[k][m][it].sbh3/newdata[k][m][it].smb/1000.0/1000.0;
	    }
	    
	  }
	  else{
	    if(dstatic[k][m].vba > ZEROVALUE) byclass_sum[it].flowvbh2+=newdata[k][m][it].vbh2/dstatic[k][m].vba/1000.0/1000.0;
	  } 
	  	    
       
	  
	  if(it == 0 ) land_garea_sum+=(chtest[k][m][it].v+chtest[k][m][it].s+chtest[k][m][it].p+chtest[k][m][it].c+chtest[k][m][it].u)*garea[k][m]; 
	  
	  if(it == 0) garea_sum+=garea[k][m];
	      

	} /* end of continent_test */

        if(it==0) total_garea_sum+=garea[k][m];


      } /* end of m*/
    } /* end of k*/
    
    
    byclass_sum[it].total+=byclass_sum[it].c+byclass_sum[it].p+
      byclass_sum[it].v+byclass_sum[it].i+
      byclass_sum[it].w+byclass_sum[it].s+byclass_sum[it].u;  
      
      
    if(fnf_s_area > ZEROVALUE) byclass_sum[it].smaf=byclass_sum[it].smaf/fnf_s_area;
    
    if(nf_s_area > ZEROVALUE) byclass_sum[it].smanf=byclass_sum[it].smanf/nf_s_area;    
    
 

    if((it == 0) && (strcmp("tone",trun) == 0)) fprintf(outfile,"yr c p v i w s u fsum vf vnf sf snf smaf smanf vbh sbh vbh2 sbh2 sbh3 sum wh clearing_amount_total clearing_amount_used wh-sum-clearing_amount_ifused(unmet) flowcp flowpc flowpv flowvp flowvc flowcv flowsp flowps flowsc flowcs flowcu flowpu flowvu flowsu flowuc flowup flowus flowvbh flowsbh flowvbh2 flowsbh2 flowsbh3 land_asum_km2 land_iw_asum_km2 total_garea_sum_km2\n");
      
     fprintf(outfile,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	    year[it],
	    byclass_sum[it].c/garea_sum,
	    byclass_sum[it].p/garea_sum,
	    byclass_sum[it].v/garea_sum,
	    byclass_sum[it].i/garea_sum,
	    byclass_sum[it].w/garea_sum,
	    byclass_sum[it].s/garea_sum,
	    byclass_sum[it].u/garea_sum,
	    byclass_sum[it].total/garea_sum,
	    byclass_sum[it].vf/garea_sum,
	    byclass_sum[it].vnf/garea_sum,
	    byclass_sum[it].sf/garea_sum,
	    byclass_sum[it].snf/garea_sum,
	    byclass_sum[it].smaf,
	    byclass_sum[it].smanf,
	    byclass_sum[it].vbh,
	    byclass_sum[it].sbh,
	    byclass_sum[it].vbh2,
	    byclass_sum[it].sbh2,
	     byclass_sum[it].sbh3,
	    byclass_sum[it].vbh+byclass_sum[it].sbh+byclass_sum[it].vbh2+byclass_sum[it].sbh2+byclass_sum[it].sbh3,
	     byclass_sum[it].wh,
	     byclass_sum[it].converted_forest_land_total,
	     byclass_sum[it].converted_forest_land,
	    byclass_sum[it].wh-(byclass_sum[it].vbh+byclass_sum[it].sbh+byclass_sum[it].vbh2+byclass_sum[it].sbh2+byclass_sum[it].sbh3+byclass_sum[it].converted_forest_land),
	     byclass_sum[it].flowcp,
	    byclass_sum[it].flowpc,
   	    byclass_sum[it].flowpv,
	    byclass_sum[it].flowvp,
	    byclass_sum[it].flowvc,
	    byclass_sum[it].flowcv,
	    byclass_sum[it].flowsp,
	    byclass_sum[it].flowps,
	    byclass_sum[it].flowsc,
	    byclass_sum[it].flowcs,
	    byclass_sum[it].flowcu,
	    byclass_sum[it].flowpu,
	    byclass_sum[it].flowvu,
	    byclass_sum[it].flowsu,
	    byclass_sum[it].flowuc,
	    byclass_sum[it].flowup,
	    byclass_sum[it].flowus,
	    byclass_sum[it].flowvbh,
	    byclass_sum[it].flowsbh,
	    byclass_sum[it].flowvbh2,
	    byclass_sum[it].flowsbh2,
	     byclass_sum[it].flowsbh3,
	    land_garea_sum/1000.0/1000.0,
	     garea_sum/1000.0/1000.0,
	     total_garea_sum/1000.0/1000.0 ); 


     if((it == 0) && (strcmp("tone",trun) == 0)) fprintf(outfile2,"yr cs vc sc vp sp ps su us vu\n");
      
     fprintf(outfile2,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	     year[it],
	     byclass_sum[it].flowcs_prime,
	     byclass_sum[it].flowvc_prime,
	     byclass_sum[it].flowsc_prime,
	     byclass_sum[it].flowvp_prime,
	     byclass_sum[it].flowsp_prime,
	     byclass_sum[it].flowps_prime,
	     byclass_sum[it].flowsu_prime,
	     byclass_sum[it].flowus_prime,
	     byclass_sum[it].flowvu_prime);   

 
  } /* end of it */

  fclose(outfile);
  fclose(outfile2);


  return;
  
}

/********************************************************************/

void regional_timeseries_checker(int regional_code, char regional_name[50]){

  FILE *outfile;
  int i, ic, k, m, it;
  double garea_sum=ZEROVALUE, land_garea_sum=ZEROVALUE, total_garea_sum=ZEROVALUE;
  double nf_s_area=ZEROVALUE, fnf_s_area=ZEROVALUE;
  char outstat[2], ffname[90];
  double tester1, tester2, tester3, tester4;



  printf("printing regional tester...\n");  

  initialize_checker(regional_code);

  printf("code %d\n",regional_code);


  strcpy(ffname,regional_name);  
  strcat(ffname,".regional.test");


  if(strcmp("tone",trun) == 0) 
    outfile=fopen(ffname,"w");
  else
    outfile=fopen(ffname,"a");


  for (it=0;it<NEWNT-1;it++){
    
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){   
	
	
	if((dstatic[k][m].rcode==regional_code) && (dstatic[k][m].usflag==1)){   	
	  
	  if(it == 0) {
	    chtest[k][m][it].c = newdata[k][m][it].c;
	    chtest[k][m][it].p = newdata[k][m][it].p;
	    chtest[k][m][it].v = newdata[k][m][it].v;
	    chtest[k][m][it].s = newdata[k][m][it].s;
	    chtest[k][m][it].u = newdata[k][m][it].u;
	  }
	  
	  chtest[k][m][it+1].c=chtest[k][m][it].c +
	    newdata[k][m][it].flowpc + newdata[k][m][it].flowvc + 
	    newdata[k][m][it].flowsc - newdata[k][m][it].flowcp - 
	    newdata[k][m][it].flowcv - newdata[k][m][it].flowcs +
	    newdata[k][m][it].flowuc - newdata[k][m][it].flowcu; 
	  
	  chtest[k][m][it+1].p=chtest[k][m][it].p +
	    newdata[k][m][it].flowcp + newdata[k][m][it].flowvp + 
	    newdata[k][m][it].flowsp - newdata[k][m][it].flowpc - 
	    newdata[k][m][it].flowpv - newdata[k][m][it].flowps +
	    newdata[k][m][it].flowup - newdata[k][m][it].flowpu;
	    
	  chtest[k][m][it+1].v=chtest[k][m][it].v +
	    newdata[k][m][it].flowpv + newdata[k][m][it].flowcv -
	    newdata[k][m][it].flowvp - newdata[k][m][it].flowvc -   
	    newdata[k][m][it].flowvs - newdata[k][m][it].flowvu;
	    
	  chtest[k][m][it+1].s=chtest[k][m][it].s +
	    newdata[k][m][it].flowvs + newdata[k][m][it].flowps + 
	    newdata[k][m][it].flowcs -
	    newdata[k][m][it].flowsp - newdata[k][m][it].flowsc +
	    newdata[k][m][it].flowus - newdata[k][m][it].flowsu;  

	  chtest[k][m][it+1].u=chtest[k][m][it].u +
	    newdata[k][m][it].flowpu + newdata[k][m][it].flowvu + 
	    newdata[k][m][it].flowsu + newdata[k][m][it].flowcu - 
	    newdata[k][m][it].flowus -
	    newdata[k][m][it].flowuc - newdata[k][m][it].flowup; 
	  
	 
	} /* end of country_test */
      } /* end of m */
    } /* end of m */
  }  /* end of it */
  
 

  for (it=0;it<NEWNT-1;it++){
    

    for (i=0;i<NCCODE;i++) {

      
      if(cdata[i].rcode == regional_code) {

	byclass_sum[it].wh+=ctdata[i][it].wh;
	byclass_sum[it].converted_forest_land+=ctdata[i][it].converted_forest_land;
	byclass_sum[it].converted_forest_land_total+=ctdata[i][it].converted_forest_land_total;

      }
      
    } /* end of i */

 
    fnf_s_area=ZEROVALUE;
    nf_s_area=ZEROVALUE;
    
    

    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){		

	if((dstatic[k][m].rcode==regional_code) && (dstatic[k][m].usflag==1)){   	

	
	  if(dstatic[k][m].fnf == 1){
		
	    byclass_sum[it].vf+=garea[k][m]*chtest[k][m][it].v;
	    byclass_sum[it].sf+=garea[k][m]*chtest[k][m][it].s;
	    byclass_sum[it].smaf+=chtest[k][m][it].s*garea[k][m]*newdata[k][m][it].sma;
	    fnf_s_area += chtest[k][m][it].s*garea[k][m];
	  }
	  else{
	    byclass_sum[it].vnf+=garea[k][m]*chtest[k][m][it].v;
	    byclass_sum[it].snf+=garea[k][m]*chtest[k][m][it].s;
	    byclass_sum[it].smanf+=chtest[k][m][it].s*garea[k][m]*newdata[k][m][it].sma;
	    nf_s_area += chtest[k][m][it].s*garea[k][m];
	  }
	    
	  byclass_sum[it].vbh+=newdata[k][m][it].vbh/1000.;
	  byclass_sum[it].sbh+=newdata[k][m][it].sbh/1000.;
	  byclass_sum[it].vbh2+=newdata[k][m][it].vbh2/1000.;
	  byclass_sum[it].sbh2+=newdata[k][m][it].sbh2/1000.;
	  byclass_sum[it].sbh3+=newdata[k][m][it].sbh3/1000.;
	  byclass_sum[it].wh_unmet+=ZEROVALUE; 
	  
	  
	  byclass_sum[it].c+=garea[k][m]*chtest[k][m][it].c; 
	  byclass_sum[it].p+=garea[k][m]*chtest[k][m][it].p;
	  byclass_sum[it].v+=garea[k][m]*chtest[k][m][it].v; 
	  byclass_sum[it].i+=garea[k][m]*newdata[k][m][it].i;  
	  byclass_sum[it].w+=garea[k][m]*newdata[k][m][it].w;  
	  byclass_sum[it].s+=garea[k][m]*chtest[k][m][it].s; 
	  byclass_sum[it].u+=garea[k][m]*chtest[k][m][it].u; 
	  
	  byclass_sum[it].flowcp+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcp;
	  byclass_sum[it].flowpc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowpc;
	  byclass_sum[it].flowpv+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowpv;
	  byclass_sum[it].flowvp+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowvp;
	  byclass_sum[it].flowvc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowvc;
	  byclass_sum[it].flowcv+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcv;
	  byclass_sum[it].flowsp+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowsp;
	  byclass_sum[it].flowps+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowps;
	  byclass_sum[it].flowsc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowsc;
	  byclass_sum[it].flowcs+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcs;
	  byclass_sum[it].flowcu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcu;
	  byclass_sum[it].flowpu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowpu;
	  byclass_sum[it].flowvu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowvu;
	  byclass_sum[it].flowsu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowsu;
	  byclass_sum[it].flowuc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowuc;
	  byclass_sum[it].flowup+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowup;
	  byclass_sum[it].flowus+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowus;

          
	  
	  if(dstatic[k][m].fnf == 1){

	    if(dstatic[k][m].vba > ZEROVALUE) byclass_sum[it].flowvbh+=newdata[k][m][it].vbh/dstatic[k][m].vba/1000.0/1000.0;

	    if(newdata[k][m][it].smb > ZEROVALUE){
	     
	      if(newdata[k][m][it].smb > ZEROVALUE) byclass_sum[it].flowsbh+=newdata[k][m][it].sbh/newdata[k][m][it].smb/1000.0/1000.0;
	      if(newdata[k][m][it].smb > ZEROVALUE) byclass_sum[it].flowsbh2+=newdata[k][m][it].sbh2/newdata[k][m][it].smb/1000.0/1000.0;

	      if(newdata[k][m][it].smb > ZEROVALUE) byclass_sum[it].flowsbh3+=newdata[k][m][it].sbh3/newdata[k][m][it].smb/1000.0/1000.0;



	    }
	    
	  }
	  else{
	    if(dstatic[k][m].vba > ZEROVALUE) byclass_sum[it].flowvbh2+=newdata[k][m][it].vbh2/dstatic[k][m].vba/1000.0/1000.0;
	  } 
	  	    
	  
	  
	  if(it == 0 ) land_garea_sum+=(chtest[k][m][it].v+chtest[k][m][it].s+chtest[k][m][it].p+chtest[k][m][it].c+chtest[k][m][it].u)*garea[k][m]; 
	  
	  if(it == 0) garea_sum+=garea[k][m];
	      

	} /* end of country_test */

        if(it==0) total_garea_sum+=garea[k][m];


      } /* end of m*/
    } /* end of k*/
    
  
    byclass_sum[it].total+=byclass_sum[it].c+byclass_sum[it].p+
      byclass_sum[it].v+byclass_sum[it].i+
      byclass_sum[it].w+byclass_sum[it].s+byclass_sum[it].u;  
      
      
    if(fnf_s_area > ZEROVALUE) byclass_sum[it].smaf=byclass_sum[it].smaf/fnf_s_area;
    
    if(nf_s_area > ZEROVALUE) byclass_sum[it].smanf=byclass_sum[it].smanf/nf_s_area;
    


    if((it == 0) && (strcmp("tone",trun) == 0)) fprintf(outfile,"yr c p v i w s u fsum vf vnf sf snf smaf smanf vbh sbh vbh2 sbh2 sbh3 sum wh clearing_amount_total clearing_amount_used wh-sum-clearing_amount_ifused(unmet) flowcp flowpc flowpv flowvp flowvc flowcv flowsp flowps flowsc flowcs flowcu flowpu flowvu flowsu flowuc flowup flowus flowvbh flowsbh flowvbh2 flowsbh2 flowsbh3 land_asum_km2 land_iw_asum_km2 total_garea_sum_km2\n");
      
    fprintf(outfile,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	    year[it],
	    byclass_sum[it].c/garea_sum,
	    byclass_sum[it].p/garea_sum,
	    byclass_sum[it].v/garea_sum,
	    byclass_sum[it].i/garea_sum,
	    byclass_sum[it].w/garea_sum,
	    byclass_sum[it].s/garea_sum,
	    byclass_sum[it].u/garea_sum,
	    byclass_sum[it].total/garea_sum,
	    byclass_sum[it].vf/garea_sum,
	    byclass_sum[it].vnf/garea_sum,
	    byclass_sum[it].sf/garea_sum,
	    byclass_sum[it].snf/garea_sum,
	    byclass_sum[it].smaf,
	    byclass_sum[it].smanf,
	    byclass_sum[it].vbh,
	    byclass_sum[it].sbh,
	    byclass_sum[it].vbh2,
	    byclass_sum[it].sbh2,
	    byclass_sum[it].sbh3,
	    byclass_sum[it].vbh+byclass_sum[it].sbh+byclass_sum[it].vbh2+byclass_sum[it].sbh2+byclass_sum[it].sbh3,
	    byclass_sum[it].wh,
	    byclass_sum[it].converted_forest_land_total,
	    byclass_sum[it].converted_forest_land,
	    byclass_sum[it].wh-(byclass_sum[it].vbh+byclass_sum[it].sbh+byclass_sum[it].vbh2+byclass_sum[it].sbh2+byclass_sum[it].sbh3+byclass_sum[it].converted_forest_land),
	    byclass_sum[it].flowcp,
	    byclass_sum[it].flowpc,
   	    byclass_sum[it].flowpv,
	    byclass_sum[it].flowvp,
	    byclass_sum[it].flowvc,
	    byclass_sum[it].flowcv,
	    byclass_sum[it].flowsp,
	    byclass_sum[it].flowps,
	    byclass_sum[it].flowsc,
	    byclass_sum[it].flowcs,
	    byclass_sum[it].flowcu,
	    byclass_sum[it].flowpu,
	    byclass_sum[it].flowvu,
	    byclass_sum[it].flowsu,
	    byclass_sum[it].flowuc,
	    byclass_sum[it].flowup,
	    byclass_sum[it].flowus,
	    byclass_sum[it].flowvbh,
	    byclass_sum[it].flowsbh,
	    byclass_sum[it].flowvbh2,
	    byclass_sum[it].flowsbh2,
	    byclass_sum[it].flowsbh3,
	    land_garea_sum/1000.0/1000.0,
	    garea_sum/1000.0/1000.0,
	    total_garea_sum/1000.0/1000.0 ); 
      
  } /* end of it */

  fclose(outfile);

  

  return;
  
}

/********************************************************************/

void loop_call_for_country_final_stats(){

  int i, country_code;


  for (i=0;i<NCCODE;i++){

     country_code=cdata[i].ccode;
     country_final_stats(country_code,i);


  } /* end of i */

  return;
}

/********************************************************************/

void country_final_stats(int country_code, int istart){

  FILE *outfile;
  int i, ic, k, m, it;
  double garea_sum=ZEROVALUE, land_garea_sum=ZEROVALUE, total_garea_sum=ZEROVALUE;
  double nf_s_area=ZEROVALUE, fnf_s_area=ZEROVALUE;
  char outstat[2], ffname[90];
  double tester1, tester2, tester3, tester4;



  initialize_checker(country_code);


  if((strcmp("tone",trun) == 0) && (istart == 0))
    outfile=fopen("country.final.stats.txt","w");
  else
    outfile=fopen("country.final.stats.txt","a");  



  for (it=0;it<NEWNT-1;it++){
    
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){   
	
	
	if((dstatic[k][m].gcode==country_code) && (dstatic[k][m].usflag==1)){   	
	  
	  if(it == 0) {
	    chtest[k][m][it].c = newdata[k][m][it].c;
	    chtest[k][m][it].p = newdata[k][m][it].p;
	    chtest[k][m][it].v = newdata[k][m][it].v;
	    chtest[k][m][it].s = newdata[k][m][it].s;
	    chtest[k][m][it].u = newdata[k][m][it].u;
	  }
	  
	  chtest[k][m][it+1].c=chtest[k][m][it].c +
	    newdata[k][m][it].flowpc + newdata[k][m][it].flowvc + 
	    newdata[k][m][it].flowsc - newdata[k][m][it].flowcp - 
	    newdata[k][m][it].flowcv - newdata[k][m][it].flowcs +
	    newdata[k][m][it].flowuc - newdata[k][m][it].flowcu; 
	  
	  chtest[k][m][it+1].p=chtest[k][m][it].p +
	    newdata[k][m][it].flowcp + newdata[k][m][it].flowvp + 
	    newdata[k][m][it].flowsp - newdata[k][m][it].flowpc - 
	    newdata[k][m][it].flowpv - newdata[k][m][it].flowps +
	    newdata[k][m][it].flowup - newdata[k][m][it].flowpu;
	    
	  chtest[k][m][it+1].v=chtest[k][m][it].v +
	    newdata[k][m][it].flowpv + newdata[k][m][it].flowcv -
	    newdata[k][m][it].flowvp - newdata[k][m][it].flowvc -   
	    newdata[k][m][it].flowvs - newdata[k][m][it].flowvu;
	    
	  chtest[k][m][it+1].s=chtest[k][m][it].s +
	    newdata[k][m][it].flowvs + newdata[k][m][it].flowps + 
	    newdata[k][m][it].flowcs -
	    newdata[k][m][it].flowsp - newdata[k][m][it].flowsc +
	    newdata[k][m][it].flowus - newdata[k][m][it].flowsu;  
	 
	  chtest[k][m][it+1].u=chtest[k][m][it].u +
	    newdata[k][m][it].flowcu + newdata[k][m][it].flowpu + 
	    newdata[k][m][it].flowvu + newdata[k][m][it].flowsu - 
	    newdata[k][m][it].flowuc - newdata[k][m][it].flowup -
	    newdata[k][m][it].flowus; 
	 

	} /* end of country_test */
      } /* end of m */
    } /* end of m */
  }  /* end of it */
  
 

  for (it=0;it<NEWNT-1;it++){
    

    for (i=0;i<NCCODE;i++) {

      if(cdata[i].ccode == country_code) {
	byclass_sum[it].wh=ctdata[i][it].wh;
	byclass_sum[it].converted_forest_land=ctdata[i][it].converted_forest_land;
	byclass_sum[it].converted_forest_land_total=ctdata[i][it].converted_forest_land_total;
      
      }
      
    } /* end of i */

 
    fnf_s_area=ZEROVALUE;
    nf_s_area=ZEROVALUE;
    
    

    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){		

	if((dstatic[k][m].gcode==country_code) && (dstatic[k][m].usflag==1)){   	

	
	  if(dstatic[k][m].fnf == 1){
		
	    byclass_sum[it].vf+=garea[k][m]*chtest[k][m][it].v;
	    byclass_sum[it].sf+=garea[k][m]*chtest[k][m][it].s;
	    byclass_sum[it].smaf+=chtest[k][m][it].s*garea[k][m]*newdata[k][m][it].sma;
	    fnf_s_area += chtest[k][m][it].s*garea[k][m];
	  }
	  else{
	    byclass_sum[it].vnf+=garea[k][m]*chtest[k][m][it].v;
	    byclass_sum[it].snf+=garea[k][m]*chtest[k][m][it].s;
	    byclass_sum[it].smanf+=chtest[k][m][it].s*garea[k][m]*newdata[k][m][it].sma;
	    nf_s_area += chtest[k][m][it].s*garea[k][m];
	  }
	    
	  byclass_sum[it].vbh+=newdata[k][m][it].vbh/1000.;
	  byclass_sum[it].sbh+=newdata[k][m][it].sbh/1000.;
	  byclass_sum[it].vbh2+=newdata[k][m][it].vbh2/1000.;
	  byclass_sum[it].sbh2+=newdata[k][m][it].sbh2/1000.;
	  byclass_sum[it].sbh3+=newdata[k][m][it].sbh3/1000.;
	  byclass_sum[it].wh_unmet+=ZEROVALUE;   
	  
	  
	  byclass_sum[it].c+=garea[k][m]*chtest[k][m][it].c; 
	  byclass_sum[it].p+=garea[k][m]*chtest[k][m][it].p;
	  byclass_sum[it].v+=garea[k][m]*chtest[k][m][it].v; 
	  byclass_sum[it].i+=garea[k][m]*newdata[k][m][it].i;  
	  byclass_sum[it].w+=garea[k][m]*newdata[k][m][it].w;  
	  byclass_sum[it].s+=garea[k][m]*chtest[k][m][it].s; 
	  byclass_sum[it].u+=garea[k][m]*chtest[k][m][it].u; 
	  
	  byclass_sum[it].flowcp+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcp;
	  byclass_sum[it].flowpc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowpc;
	  byclass_sum[it].flowpv+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowpv;
	  byclass_sum[it].flowvp+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowvp;
	  byclass_sum[it].flowvc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowvc;
	  byclass_sum[it].flowcv+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcv;
	  byclass_sum[it].flowsp+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowsp;
	  byclass_sum[it].flowps+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowps;
	  byclass_sum[it].flowsc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowsc;
	  byclass_sum[it].flowcs+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcs;
	  byclass_sum[it].flowcu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowcu;
	  byclass_sum[it].flowpu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowpu;
	  byclass_sum[it].flowvu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowvu;
	  byclass_sum[it].flowsu+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowsu;
	  byclass_sum[it].flowuc+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowuc;
	  byclass_sum[it].flowup+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowup;
	  byclass_sum[it].flowus+=garea[k][m]/1000.0/1000.0*newdata[k][m][it].flowus;

          
	  
	  if(dstatic[k][m].fnf == 1){

	    if(dstatic[k][m].vba > ZEROVALUE) byclass_sum[it].flowvbh+=newdata[k][m][it].vbh/dstatic[k][m].vba/1000.0/1000.0;

	    if(newdata[k][m][it].smb > ZEROVALUE){
	     
	      if(newdata[k][m][it].smb > ZEROVALUE) byclass_sum[it].flowsbh+=newdata[k][m][it].sbh/newdata[k][m][it].smb/1000.0/1000.0;
	      if(newdata[k][m][it].smb > ZEROVALUE) byclass_sum[it].flowsbh2+=newdata[k][m][it].sbh2/newdata[k][m][it].smb/1000.0/1000.0;
	      if(newdata[k][m][it].smb > ZEROVALUE) byclass_sum[it].flowsbh3+=newdata[k][m][it].sbh3/newdata[k][m][it].smb/1000.0/1000.0;
	    }
	    
	  }
	  else{
	    if(dstatic[k][m].vba > ZEROVALUE) byclass_sum[it].flowvbh2+=newdata[k][m][it].vbh2/dstatic[k][m].vba/1000.0/1000.0;
	  } 
	  	    
	  
	  
	  if(it == 0 ) land_garea_sum+=(chtest[k][m][it].v+chtest[k][m][it].s+chtest[k][m][it].p+chtest[k][m][it].c+chtest[k][m][it].u)*garea[k][m]; 
	  
	  if(it == 0) garea_sum+=garea[k][m];
	      

	} /* end of country_test */

        if(it==0) total_garea_sum+=garea[k][m];


      } /* end of m*/
    } /* end of k*/
    
  
    byclass_sum[it].total+=byclass_sum[it].c+byclass_sum[it].p+
      byclass_sum[it].v+byclass_sum[it].i+
      byclass_sum[it].w+byclass_sum[it].s+byclass_sum[it].u;  
      
      
    if(fnf_s_area > ZEROVALUE) byclass_sum[it].smaf=byclass_sum[it].smaf/fnf_s_area;
    
    if(nf_s_area > ZEROVALUE) byclass_sum[it].smanf=byclass_sum[it].smanf/nf_s_area;
    
    if(hist_option==1){

      if((strcmp("tone",trun) == 0) && (istart == 0) && (year[it] == 1700)) fprintf(outfile,"yr cname c p v i w s u fsum vf vnf sf snf smaf smanf vbh sbh vbh2 sbh2 sbh3 sum wh clearing_amount_total clearing_amount_used wh-sum-clearing_amount_ifused(unmet) flowcp flowpc flowpv flowvp flowvc flowcv flowsp flowps flowsc flowcs flowcu flowpu flowvu flowsu flowuc flowup flowus flowvbh flowsbh flowvbh2 flowsbh3 flowsbh2 land_asum_km2 land_iw_asum_km2 total_garea_sum_km2\n");
    }
    else if(hist_option==2){
      if((strcmp("tone",trun) == 0) && (istart == 0) && (year[it] == 1500)) fprintf(outfile,"yr cname c p v i w s u fsum vf vnf sf snf smaf smanf vbh sbh vbh2 sbh2 sbh3 sum wh clearing_amount_total clearing_amount_used wh-sum-clearing_amount_ifused(unmet) flowcp flowpc flowpv flowvp flowvc flowcv flowsp flowps flowsc flowcs flowcu flowpu flowvu flowsu flowuc flowup flowus flowvbh flowsbh flowvbh2 flowsbh3 flowsbh2 land_asum_km2 land_iw_asum_km2 total_garea_sum_km2\n");
    }
    else if(hist_option==3){
      if((strcmp("tone",trun) == 0) && (istart == 0) && (year[it] == 1850)) fprintf(outfile,"yr cname c p v i w s u fsum vf vnf sf snf smaf smanf vbh sbh vbh2 sbh2 sbh3 sum wh clearing_amount_total clearing_amount_used wh-sum-clearing_amount_ifused(unmet) flowcp flowpc flowpv flowvp flowvc flowcv flowsp flowps flowsc flowcs flowcu flowpu flowvu flowsu flowuc flowup flowus flowvbh flowsbh flowvbh2 flowsbh3 flowsbh2 land_asum_km2 land_iw_asum_km2 total_garea_sum_km2\n");
    }
      
    fprintf(outfile,"%d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	    year[it],cname[istart],
	    byclass_sum[it].c/garea_sum,
	    byclass_sum[it].p/garea_sum,
	    byclass_sum[it].v/garea_sum,
	    byclass_sum[it].i/garea_sum,
	    byclass_sum[it].w/garea_sum,
	    byclass_sum[it].s/garea_sum,
	    byclass_sum[it].u/garea_sum,
	    byclass_sum[it].total/garea_sum,
	    byclass_sum[it].vf/garea_sum,
	    byclass_sum[it].vnf/garea_sum,
	    byclass_sum[it].sf/garea_sum,
	    byclass_sum[it].snf/garea_sum,
	    byclass_sum[it].smanf,
	    byclass_sum[it].vbh,
	    byclass_sum[it].sbh,
	    byclass_sum[it].vbh2,
	    byclass_sum[it].sbh2,
	    byclass_sum[it].sbh3,
	    byclass_sum[it].vbh+byclass_sum[it].sbh+byclass_sum[it].vbh2+byclass_sum[it].sbh2+byclass_sum[it].sbh3,
	    byclass_sum[it].wh,
	    byclass_sum[it].converted_forest_land_total,
	    byclass_sum[it].converted_forest_land,
	    byclass_sum[it].wh-(byclass_sum[it].vbh+byclass_sum[it].sbh+byclass_sum[it].vbh2+byclass_sum[it].sbh2+byclass_sum[it].sbh3+byclass_sum[it].converted_forest_land),
	    byclass_sum[it].flowcp,
	    byclass_sum[it].flowpc,
   	    byclass_sum[it].flowpv,
	    byclass_sum[it].flowvp,
	    byclass_sum[it].flowvc,
	    byclass_sum[it].flowcv,
	    byclass_sum[it].flowsp,
	    byclass_sum[it].flowps,
	    byclass_sum[it].flowsc,
	    byclass_sum[it].flowcs,
	    byclass_sum[it].flowcu,
	    byclass_sum[it].flowpu,
	    byclass_sum[it].flowvu,
	    byclass_sum[it].flowsu,
	    byclass_sum[it].flowuc,
	    byclass_sum[it].flowup,
	    byclass_sum[it].flowus,
	    byclass_sum[it].flowvbh,
	    byclass_sum[it].flowsbh,
	    byclass_sum[it].flowvbh2,
	    byclass_sum[it].flowsbh2,
	    byclass_sum[it].flowsbh3,
	    land_garea_sum/1000.0/1000.0,
	    garea_sum/1000.0/1000.0,
	    total_garea_sum/1000.0/1000.0 ); 
      
  } /* end of it */

  fclose(outfile);

  

  return;
  
}

/********************************************************************/
