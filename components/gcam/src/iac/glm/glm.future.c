// NOTES: currently opening and closing netcdf files every timestep could
//        keep open for length of run.
//  check header for future 27 vs 25
//  hist vs future nreg 21 vs 14

#include "dictionary.h"
#include "iniparser.h"
#include <math.h>
#include <netcdf.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#define CLASS_NUM 5      /*c,p,o,i,w*/
#define FLEN 14
#define HEADLEN 6
#define NCCODE 192       /* number of countries 192 */
#define NEWNT 2         /* 102, 62, or 32 */   //jt add
#define NPHB 50   /* filelength of probability of harvest given biomass */
#define NREG 14
#define NX 720           /* 360 or 720 */
#define NY 360           /* 180 or 360 */
#define NZ 31            /* Number of z dims tracked 31 (z is distance from managed cells) */
#define ZEROVALUE 0.000000000000000000000
#define ZEROVALUE_CHECK 0.000000001

//#define NEWNT 51         /*51     101, 61,61,61,61,61,56,56  */
#define MAXZ_21 1
#define DMAX 31

  
#define ZEROVALUE 0.000000000000000000000
#define ZEROVALUE_CHECK 0.000000001

char casestr[256];
char config_file[256];
char restart_file[256];
char* CAT3; 
char* CAT5; 
char* CAT6; 
char* PATH1;
char* PROB_FNAME;       /* "phbio.average.7states.txt"*/
char* casename;
char* ccodes_file;
char* ccodes_map;
char* cellinfo_file;
char* cnames_file;
char* contcodes_file;
char* finput_dir;
char* future_crop_path;
char* future_icew_path;
char* future_initial_state;
char* future_othr_path;
char* future_past_path;
char* future_secd_path;
char* future_watr_path;
char* gcodes_cont_map;
char* hyde_crop_path;
char* hyde_icew_path;
char* hyde_othr_path;
char* hyde_past_path;
char* hyde_watr_path;
char* miami_biomass_file_vba;
char* miami_biomass_file_vnppa;
char* outstat;
char* rcp_woodharvest_file;
char* rcp_woodharvest_nodata_file;
char* regcodes_map;
char* regnames_file;
char* restart_filename;
char* runtype;
char* s;
char* shiftcult_map;
char* static_regions_file;
char* static_vba_file;
char* trun;
char* whcodes2glm_map;
char* whcodes_file;
char* whcontcodes_file;
char* woodharvest_file;
char* woodharvest_nodata_file;
double TB2BB; /* total biomass to bole biomass ratio, value of 2.0 */
              /* assumes the WH numbers are bole biomass only */
	      /* and we need to cut twice as much biomass  */
int base_year;  /* netCDF file id for state */
int coterminous;  /* netCDF file id for state */
int initialrun;  /* netCDF file id for state */
int input_netcdf_files;  /* netCDF file id for state */
int input_text_files;  /* netCDF file id for state */
int ncid_state;  /* netCDF file id for state */
int print_file_headers;
int output_netcdf;  /* netCDF file id for state */
int output_text;  /* netCDF file id for state */
int region_test,region_test_gcode,region_test_index;
int future_rate_option, nodata_option, logging_option, res_option, hist_option;
int start_year,stop_year,restart,netcdf;  //jtadd
int yrstart,trun_years;

//global params read in from config file
int MAXZ;
int BEST_CASE;
int BEST_CASE_MIN_FLOWS_T5;
int BEST_CASE_MIN_FLOWS_T4;
int TOTAL_HARVEST_SWITCH;
int SECONDARY_HARVEST_SWITCH;
int VIRGIN_HARVEST_SWITCH;
int FORCE_HARVEST_SWITCH;
int CPAVG;
int smart_flow_option, converted_forest_land_option;
int zdis_option, adjust_smart_flow_option, harvest_option;
int initial_run;
int SMART_FLOW_BUG_PRINT;
int STATE_BUG_PRINT;
int STATE_PRINT;
int FLOW_BUG_PRINT;
int REGION_TEST_TMP;           /* 0=global, 1=country */  
int REGION_TEST_GCODE;         /* 221=US, 34=CANADA, 28=BRAZIL, CAMBODIA=35, AUS=11 */
int REGION_TEST_INDEX;         /* 126=US, 20=CANADA, 16=BRAZIL, CAMBODIA=21, AUS=6 */
int CONTERMINOUS;              /* 1=CONTERMINOUS US ONLY,0=ALL US */

int yrstart;
int future_rate_option, nodata_option, logging_option, suit_option, res_option, newnt_option, hist_option, gridded_wh,use_urban;
//jtchar trun[10], finput_dir[150], hist_data[150];
char hist_data[150];
double popsum_agriculture=ZEROVALUE;
double wh_ratio[NCCODE];
double hist_reg_wh[NREG];


struct new_data {
  int zdis[NY][NX];
  double vbh[NY][NX], sbh[NY][NX], sma[NY][NX], smb[NY][NX], flowvs[NY][NX], sbh2[NY][NX], vbh2[NY][NX], sbh3[NY][NX];
  double c[NY][NX],p[NY][NX],v[NY][NX],i[NY][NX],w[NY][NX],s[NY][NX],u[NY][NX];  
  double flowcp[NY][NX], flowpc[NY][NX], flowpv[NY][NX], flowvp[NY][NX], flowvc[NY][NX], flowcv[NY][NX];
  double flowsp[NY][NX],flowps[NY][NX],flowsc[NY][NX],flowcs[NY][NX], whd[NY][NX], wh_unmet[NY][NX], converted_forest_land[NY][NX], converted_forest_land_total[NY][NX];  
  double flowcu[NY][NX],flowpu[NY][NX],flowsu[NY][NX],flowvu[NY][NX],flowus[NY][NX],flowuc[NY][NX],flowup[NY][NX],flowvbh[NY][NX],flowvbh2[NY][NX],flowsbh[NY][NX],flowsbh2[NY][NX],flowsbh3[NY][NX];
} newdata[NEWNT];  

struct new_data_scale {
  double vbh[NY][NX], sbh[NY][NX], sbh2[NY][NX], vbh2[NY][NX], sbh3[NY][NX];
} newdata_scale;  

struct rzt_data {
  double vb;                  /* units MgC */
} rztdata[NREG][NZ];

struct czt_data {
  double vb;                  /* units MgC */
} cztdata[NCCODE][NZ];

struct rt_data {
  double wh, whr, stbh, vwh;  /* units MgC */
  double wh_at_zmax, wh_lessthan_zmax;
  double smb, vnfb, sarea, smb_nf, sarea_nf, fh_sbh2, fh_sbh3, fh_vbh2;
  double converted_forest_land;
  double converted_forest_land_total;
  double flowvc_prime,flowsc_prime,flowcs_prime,flowvp_prime,flowsp_prime,flowps_prime;
  double flowvu_prime,flowsu_prime, flowus_prime;
  double reg_avail, reg_unmet;
} rtdata[NREG];

struct ct_data {
  double wh, whr, stbh, vwh, wh1, orig_whd;  /* units MgC */
  double wh_at_zmax, wh_lessthan_zmax;
  double smb, vnfb, sarea, smb_nf, sarea_nf, fh_sbh2, fh_sbh3, fh_vbh2, predict_b;
  double converted_forest_land;
  double converted_forest_land_total;
  double flowvc_prime,flowsc_prime,flowcs_prime,flowvp_prime,flowsp_prime,flowps_prime;
} ctdata[NCCODE];

struct c_data {
  int ccode, newccode, continent_code, rcode;
  int smart_flow_option, converted_forest_land_option;
  int zdis_option, harvest_option;
  double wh2005;
} cdata[NCCODE];

struct r_data {
  int rcode, continent_code;  
  int smart_flow_option, converted_forest_land_option;
  int zdis_option, adjust_smart_flow_option, harvest_option;
} rdata[NREG];

struct other_data {
  int usflag, gcode, newgcode, continent_code, fnf, rcode, newrcode, shiftcult;
  double vba, vnppa;
} dstatic[NY][NX];


char header[HEADLEN][27];
double lat[NY], garea[NY][NX], lon[NX], phb[NPHB];
/* int year[NEWNT]; */
char cname[NCCODE][50], rname[NREG][50];

   
struct sum_check {
  double c,p,v,i,w,s,u;
  double vf,vnf,sf,snf;
  double smaf,smanf,vbh,sbh,vbh2,sbh2,sbh3,wh_unmet;
  double total, wh;  
  double flowcp, flowpc, flowpv, flowvp, flowvc, flowcv;
  double flowcu,flowpu,flowvu,flowsu,flowuc,flowup,flowus;
  double flowsp,flowps,flowsc,flowcs;
  double flowsbh, flowsbh2, flowvbh, flowvbh2, flowsbh3;
  double converted_forest_land, converted_forest_land_total;
  double flowvc_prime,flowsc_prime,flowcs_prime,flowvp_prime,flowsp_prime,flowps_prime;
  double flowvu_prime,flowsu_prime,flowus_prime;
} byclass_sum;



struct check_data {
  double c,p,v,s,u;
} chtest[NY][NX][NEWNT];  


void adjust_smart_flow1(int k, int m, int it, int i);
void adjust_smart_flow2(int k, int m, int it, int i);
void adjust_smart_flow3(int k, int m, int it, int i);
void alternative_smart_flow1(int k, int m, int it);
void alternative_smart_flow2(int k, int m, int it);
void cellinfo();
void check_err(const int stat, const int line, const char *file);
void compute_country_flows(int i, int it);
void continent_timeseries_checker(int continent_code, char continent_name[50], int curryear);
void country_final_stats(int country_code, int istart,int curryear);
void country_primeflow_print(int curryear);
void country_timeseries_checker(int country_code, char country_name[50]);
void country_unmet_wh(int it, int i);
void create_restart_file(void);
void finalize(void);
void force_harvest(int it, int i, int curryear);
void force_harvest2(int it, int i, int curryear);
void force_harvest3(int it, int i, int curryear);
void global_timeseries_checker(int regional_code, char regional_name[50],int curryear);
void harvest_gridded_data(int it);
void initialize();
void initialize_checker(int regional_code);
void initialize_woodharvest_country_ratios(int baseyear);
void loop_call_for_country_final_stats(int curryear);
int  mkdir2 (char *str1, char *str2, int mode);
int  open2 (char *str1, char *str2, int flags, int mode);
void option_settings(int smart_flow_option,
		     int converted_forest_land_option,
		     int zdis_option,
		     int adjust_smart_flow_option,
		     int harvest_option);
void output_lu(int curryear);
void output_lu_nc(int curryear);
void output_updated_states(int curryear);
void output_updated_states2(int curryear);
void output_updated_states3(int curryear);
void output_updated_states_nc(int curryear);
void predict_available_biomass(int it, int i);
float prob_harv(float biomass);
void read_continent_codes();
void read_country_codes();
void read_country_names();
/* void read_data_futureruns(int curryear); */
void read_data_water_ice_nc(int curryear);
void read_data_futureruns_nc(int curryear);
void read_other_data();
void read_restart(void);
void read_regional_codes();
void read_regional_names();
/* void read_woodharvest_data(int curryear); */
void read_woodharvest_data_nc(int curryear);
void regional_timeseries_checker(int regional_code, char regional_name[50], int curryear);
void run_model();
void secondary_harvest(int it, int i);
void secondary_harvest2(int it, int i);
void secondary_harvest3(int it, int i);
void smart_flow(int k, int m, int it);
char* strcatjt(char *dest, const char *src);
void step(int *curryear);;
void stepglm_ccsm(int *year,double *glmi,int *glmi_fdim1, int *glmi_fdim2,double *glmi_wh,int *glmi_wh_fdim1, double *glmo, int *glmo_fdim1, int *glmo_fdim2);
void transitions(int curryear);
void update_states(int it, int zmax);
void update_vb(int it);
void update_vb2(int it);
void update_vb3(int it);
void virgin_harvest(int it, int i,int curryear);
void virgin_harvest2(int it, int i,int curryear);
void virgin_harvest3(int it, int i,int curryear);
void zdis_calc(int it);
void zdis_calc2(int it);
inline double fround(double n, unsigned d);
/*****************************************************************************/
#if (defined MAIN)
int main(int argc, char *argv[]){

#include "iniparser.h"

  setbuf(stdout, NULL);

  printf("\n");
  if( argc != 1 )
    {
      printf("Usage Error:\tIncorrect number of arguments\n\n");
      printf("usage: %s <file_name>\n", argv[0]);
      return;
    }

  initialize();

  run_model();

  finalize();

  return;
}
#endif
/*****************************************************************************/
void option_settings(int smart_flow_option,
		     int converted_forest_land_option,
		     int zdis_option,
		     int adjust_smart_flow_option,
		     int harvest_option) {
  int i;

  printf("\noption_settings...\n");

  for (i=0;i<NREG;i++){
    rdata[i].smart_flow_option=smart_flow_option;
    rdata[i].converted_forest_land_option=converted_forest_land_option;
    rdata[i].zdis_option=zdis_option;
    rdata[i].adjust_smart_flow_option=adjust_smart_flow_option;
    rdata[i].harvest_option=harvest_option;
  }

  if(rdata[i].zdis_option == 1) suit_option=1;
  if(rdata[i].zdis_option == 2) suit_option=1;

}
/***********************************************************************/
void initialize(){

  dictionary* ini ;
  int status;
  int i, iz, k, m, it, j;
  char runtypebuf[50];
  char new_path[90];
  int varidrst_past,
    varidrst_othr,
    varidrst_urbn,
    varidrst_secd,
    varidrst_suit,
    varidrst_ice,
    varidrst_water,
    varidrst_crop,
    varidrst_sma,
    varidrst_smb,
    varidrst_pasture;                          /* netCDF ID */
  size_t count[] = {1, NY, NX};
  size_t start[] = {0, 0, 0}; /* start at first value */
  
  printf("\n\nPROGRAM: 2005-2100 GLM\n\n");

  /* Load up parser file */
  ini = iniparser_load("glm.fut.conf");
  iniparser_dump(ini, stderr);

  /* Get output from config file */
  s          = iniparser_getstring(ini, "output directory:output dir", NULL);
  printf("main output directory:     [%s]\n", s ? s : "UNDEF");

  /* Get casename from config file */
  casename   = iniparser_getstring(ini, "control:case name", NULL);
  printf("casename:     [%s]\n", casename ? casename : "UNDEF");

  /* create output/lu/tester/updated_states dirs */
  status     =  mkdir(s,0755);
  status     =  mkdir2(s,"/lu",0755);
  status     =  mkdir2(s,"/tester",0755);
  status     =  mkdir2(s,"/updated_states",0755);
  status     =  mkdir2(s,"/casename",0755);

  /* parse control options from config file */

  future_rate_option   = iniparser_getint(ini, "control:future rate option", -999);
  nodata_option        = iniparser_getint(ini, "control:nodata_option", -999);
  logging_option       = iniparser_getint(ini, "control:logging_option", -999);
  suit_option          = iniparser_getint(ini, "control:suit_option", -999);
  start_year           = iniparser_getint(ini, "control:start year", -999);
  stop_year            = iniparser_getint(ini, "control:stop year", -999);
  stop_year            = iniparser_getint(ini, "control:stop year", -999);
  res_option           = iniparser_getint(ini, "control:res option", -999);
  smart_flow_option    = iniparser_getint(ini, "control:smart_flow_option", -999);
  converted_forest_land_option=iniparser_getint(ini, "control:converted_forest_land_option", -999);
  zdis_option          = iniparser_getint(ini, "control:zdis_option", -999);
  adjust_smart_flow_option = iniparser_getint(ini, "control:adjust_smart_flow_option",-999);
  harvest_option       = iniparser_getint(ini, "control:harvest_option",-999);
  use_urban            = iniparser_getint(ini, "control:use_urban",-999);


//global params read in from config file
  MAXZ                 = iniparser_getint(ini,"control:maxz",-999);
  BEST_CASE            = iniparser_getint(ini,"control:best_case",-999);
  BEST_CASE_MIN_FLOWS_T5= iniparser_getint(ini,"control:best_case_min_flows_t5",-999);
  TOTAL_HARVEST_SWITCH = iniparser_getint(ini,"control:total_harvest_switch",-999);
  SECONDARY_HARVEST_SWITCH= iniparser_getint(ini,"control:secondary_harvest_switch",-999);
  VIRGIN_HARVEST_SWITCH= iniparser_getint(ini,"control:virgin_harvest_switch",-999);
  FORCE_HARVEST_SWITCH = iniparser_getint(ini,"control:force_harvest_switch",-999);
  CPAVG                = iniparser_getint(ini,"control:cpavg",-999);
  BEST_CASE_MIN_FLOWS_T4= iniparser_getint(ini,"control:best_case_min_flows_t4",-999);
  REGION_TEST_TMP      = iniparser_getint(ini,"control:region_test_tmp",-999);
  REGION_TEST_GCODE    = iniparser_getint(ini,"control:region_test_gcode",-999);
  REGION_TEST_INDEX    = iniparser_getint(ini,"control:region_test_index",-999);
  CONTERMINOUS          = iniparser_getint(ini,"control:conterminous",-999);


  /* total biomass to bole biomass ratio, value of 2.0 assumes the */
  /*WH numbers are bole biomass only and we need to cut twice as much biomass  */
  TB2BB                = iniparser_getdouble(ini,"control:tb2bb",1.0);   
  PROB_FNAME           = iniparser_getstring(ini,"control:phbio_filename",NULL);
  runtype              = iniparser_getstring(ini,"control:runtype",NULL);
  SMART_FLOW_BUG_PRINT = iniparser_getint(ini,"debug options:smart_flow_bug_print",-999);
  STATE_BUG_PRINT      = iniparser_getint(ini,"debug options:state_bug_print",-999);
  STATE_PRINT          = iniparser_getint(ini,"debug options:state_print",-999);
  FLOW_BUG_PRINT       = iniparser_getint(ini,"debug options:flow_bug_print",-999);
  PATH1                = iniparser_getstring(ini,"control:top level glm dir",NULL);  
  output_text          = iniparser_getint(ini,"control:output_text",-999);
  output_netcdf        = iniparser_getint(ini,"control:output_netcdf",-999);
  input_text_files     = iniparser_getint(ini,"control:input_text_files",-999);
  input_netcdf_files   = iniparser_getint(ini,"control:input_netcdf_files",-999);

  hyde_crop_path       = iniparser_getstring(ini,"hyde_datasets:hyde_crop_path",NULL);
  hyde_othr_path       = iniparser_getstring(ini,"hyde_datasets:hyde_othr_path",NULL);
  hyde_past_path       = iniparser_getstring(ini,"hyde_datasets:hyde_past_path",NULL);
  hyde_icew_path       = iniparser_getstring(ini,"hyde_datasets:hyde_icew_path",NULL);
  hyde_watr_path       = iniparser_getstring(ini,"hyde_datasets:hyde_watr_path",NULL);
  woodharvest_file     = iniparser_getstring(ini,"woodharvest_datasets:woodharvest_file",NULL);
  woodharvest_nodata_file = iniparser_getstring(ini,"woodharvest_datasets:woodharvest_nodata_file",NULL);
  cellinfo_file           = iniparser_getstring(ini,"other_datasets:cellinfo_file",NULL);
  ccodes_file             = iniparser_getstring(ini,"other_datasets:ccodes_file",NULL);
  ccodes_map              = iniparser_getstring(ini,"other_datasets:ccodes_map",NULL);
  cnames_file             = iniparser_getstring(ini,"other_datasets:cnames_file",NULL);
  regnames_file           = iniparser_getstring(ini,"other_datasets:regnames_file",NULL);
  contcodes_file          = iniparser_getstring(ini,"other_datasets:contcodes_file",NULL);
  shiftcult_map           = iniparser_getstring(ini,"other_datasets:shiftcult_map",NULL);
  whcodes_file            = iniparser_getstring(ini,"other_datasets:whcodes_file",NULL);
  whcontcodes_file        = iniparser_getstring(ini,"other_datasets:whcontcodes_file",NULL);
  whcodes2glm_map         = iniparser_getstring(ini,"other_datasets:whcodes2glm_map",NULL);
  regcodes_map            = iniparser_getstring(ini,"other_datasets:regcodes_map",NULL);
  gcodes_cont_map         = iniparser_getstring(ini,"other_datasets:gcodes_cont_map",NULL);
  miami_biomass_file_vba = iniparser_getstring(ini,"other_datasets:miami_biomass_file_vba",NULL);
  miami_biomass_file_vnppa = iniparser_getstring(ini,"other_datasets:miami_biomass_file_vnppa",NULL);
  future_crop_path       = iniparser_getstring(ini,"future_datasets:future_crop_path",NULL);
  future_othr_path       = iniparser_getstring(ini,"future_datasets:future_othr_path",NULL);
  future_past_path       = iniparser_getstring(ini,"future_datasets:future_past_path",NULL);
  future_icew_path       = iniparser_getstring(ini,"future_datasets:future_icew_path",NULL);
  future_initial_state   = iniparser_getstring(ini,"future_datasets:future_initial_state",NULL);
  future_watr_path       = iniparser_getstring(ini,"future_datasets:future_watr_path",NULL);
  future_secd_path       = iniparser_getstring(ini,"future_datasets:future_secd_path",NULL);

  static_regions_file    = iniparser_getstring(ini,"output_files:static_region_file",NULL);
  woodharvest_file       = iniparser_getstring(ini,"woodharvest_datasets:woodharvest_file",NULL);
  rcp_woodharvest_file   = iniparser_getstring(ini,"woodharvest_datasets:rcp_woodharvest_file",NULL);

  //needed for original way of running glm.c
  hist_option          = iniparser_getint(ini,"old control:hist_option",-999);
  trun                 = iniparser_getstring(ini,"old control:trun",NULL);
  finput_dir           = iniparser_getstring(ini,"old control:foutput_dir",NULL);
  CAT3                 = iniparser_getstring(ini,"old control:cat3",NULL);
  CAT5                 = iniparser_getstring(ini,"old control:cat5",NULL);
  CAT6                 = iniparser_getstring(ini,"old control:cat6",NULL);

  strcpy(runtypebuf,runtype);
  //jtfix  if(strcmp("initial",lowercase_string(runtypebuf)) == 0) {initialrun=1;}
  if(strcmp("initial",runtype) == 0) {initialrun=1;}

  option_settings(smart_flow_option,
		  converted_forest_land_option,
		  zdis_option,
		  adjust_smart_flow_option,
		  harvest_option);


  
  if(strcmp("tone",trun) == 0) {yrstart=1700;}
  else if(strcmp("ttwo",trun) == 0) {yrstart=1760;}
  else if(strcmp("tthree",trun) == 0) {yrstart=1820;}
  else if(strcmp("tfour",trun) == 0) {yrstart=1880;}
  else if(strcmp("tfive",trun) == 0) {yrstart=1940;}
  else if(strcmp("tsix",trun) == 0) {yrstart=2005;}
  else {yrstart=2050;}


/*   year[0]=yrstart; */
/*   for (it=1;it<NEWNT;it++){ */
/*       year[it]=year[0]+it; */
/*   } */


  for (i=0;i<NREG;i++) {

    hist_reg_wh[i]=ZEROVALUE;
    wh_ratio[i]=ZEROVALUE;
    
    rtdata[i].wh=ZEROVALUE;
    rtdata[i].whr=ZEROVALUE;
    rtdata[i].stbh=ZEROVALUE;
    rtdata[i].vwh=ZEROVALUE;
    
    rtdata[i].wh_lessthan_zmax=ZEROVALUE;
    rtdata[i].wh_at_zmax=ZEROVALUE;
    rtdata[i].vnfb=ZEROVALUE; 
    rtdata[i].smb=ZEROVALUE;
    rtdata[i].fh_sbh2=ZEROVALUE;
    rtdata[i].fh_sbh3=ZEROVALUE;
    rtdata[i].fh_vbh2=ZEROVALUE;
    rtdata[i].sarea=ZEROVALUE;
    rtdata[i].sarea_nf=ZEROVALUE;
    rtdata[i].smb_nf=ZEROVALUE;
    rtdata[i].converted_forest_land=ZEROVALUE;
    rtdata[i].converted_forest_land_total=ZEROVALUE;
    rtdata[i].flowvc_prime=ZEROVALUE;
    rtdata[i].flowsc_prime=ZEROVALUE;
    rtdata[i].flowcs_prime=ZEROVALUE;
    rtdata[i].flowvp_prime=ZEROVALUE;
    rtdata[i].flowsp_prime=ZEROVALUE;
    rtdata[i].flowps_prime=ZEROVALUE;
    rtdata[i].reg_unmet=ZEROVALUE;
    rtdata[i].reg_avail=ZEROVALUE;
    
    for (iz=0;iz<NZ;iz++){
      rztdata[i][iz].vb=ZEROVALUE;
    }
  }
  for (i=0;i<NCCODE;i++) {
    ctdata[i].wh=ZEROVALUE;
    ctdata[i].whr=ZEROVALUE;
    ctdata[i].stbh=ZEROVALUE;
    ctdata[i].vwh=ZEROVALUE;
    ctdata[i].wh1=ZEROVALUE;
    ctdata[i].orig_whd=ZEROVALUE;
    
    ctdata[i].wh_lessthan_zmax=ZEROVALUE;
    ctdata[i].wh_at_zmax=ZEROVALUE;
    ctdata[i].vnfb=ZEROVALUE; 
    ctdata[i].smb=ZEROVALUE;
    ctdata[i].fh_sbh2=ZEROVALUE;
    ctdata[i].fh_sbh3=ZEROVALUE;
    ctdata[i].fh_vbh2=ZEROVALUE;
    ctdata[i].sarea=ZEROVALUE;
    ctdata[i].sarea_nf=ZEROVALUE;
    ctdata[i].smb_nf=ZEROVALUE;
    ctdata[i].converted_forest_land=ZEROVALUE;
    ctdata[i].converted_forest_land_total=ZEROVALUE;
    ctdata[i].predict_b=ZEROVALUE;
    ctdata[i].flowvc_prime=ZEROVALUE;
    ctdata[i].flowsc_prime=ZEROVALUE;
    ctdata[i].flowcs_prime=ZEROVALUE;
    ctdata[i].flowvp_prime=ZEROVALUE;
    ctdata[i].flowsp_prime=ZEROVALUE;
    ctdata[i].flowps_prime=ZEROVALUE;
    
    for (iz=0;iz<NZ;iz++){
      cztdata[i][iz].vb=ZEROVALUE;
    }
  }
    
  for (it=0;it<NEWNT;it++){
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){
	

	newdata[it].c[k][m]=ZEROVALUE;
	newdata[it].p[k][m]=ZEROVALUE;
	newdata[it].v[k][m]=ZEROVALUE;
	newdata[it].s[k][m]=ZEROVALUE;
	newdata[it].i[k][m]=ZEROVALUE;
	newdata[it].w[k][m]=ZEROVALUE;
	newdata[it].u[k][m]=ZEROVALUE;
	newdata[it].flowcp[k][m]=ZEROVALUE;
	newdata[it].flowpc[k][m]=ZEROVALUE;
	newdata[it].flowpv[k][m]=ZEROVALUE;
	newdata[it].flowvp[k][m]=ZEROVALUE;
	newdata[it].flowvc[k][m]=ZEROVALUE;
	newdata[it].flowcv[k][m]=ZEROVALUE;
	newdata[it].flowsc[k][m]=ZEROVALUE;
	newdata[it].flowcs[k][m]=ZEROVALUE;
	newdata[it].flowsp[k][m]=ZEROVALUE;
	newdata[it].flowps[k][m]=ZEROVALUE;
	newdata[it].flowvs[k][m]=ZEROVALUE;
	newdata[it].flowcu[k][m]=ZEROVALUE;
	newdata[it].flowpu[k][m]=ZEROVALUE;
	newdata[it].flowsu[k][m]=ZEROVALUE;
	newdata[it].flowvu[k][m]=ZEROVALUE;
	newdata[it].flowuc[k][m]=ZEROVALUE;
	newdata[it].flowup[k][m]=ZEROVALUE;
	newdata[it].flowus[k][m]=ZEROVALUE;
	newdata[it].whd[k][m]=ZEROVALUE;
	newdata[it].wh_unmet[k][m]=ZEROVALUE;
	newdata[it].converted_forest_land[k][m]=ZEROVALUE;
	newdata[it].converted_forest_land_total[k][m]=ZEROVALUE;
        newdata[it].sma[k][m]=ZEROVALUE;
        newdata[it].smb[k][m]=ZEROVALUE;
        newdata[it].sbh2[k][m]=ZEROVALUE;
        newdata[it].sbh3[k][m]=ZEROVALUE;
        newdata[it].vbh2[k][m]=ZEROVALUE;
	newdata[it].vbh[k][m]=ZEROVALUE;
	newdata[it].sbh[k][m]=ZEROVALUE;  
 	dstatic[k][m].gcode=0;
	dstatic[k][m].rcode=0;
        dstatic[k][m].newgcode=9999; 
        dstatic[k][m].newrcode=0; 
	dstatic[k][m].continent_code=0;
	dstatic[k][m].fnf=0;
	dstatic[k][m].vba=ZEROVALUE;
	dstatic[k][m].vnppa=ZEROVALUE;
	dstatic[k][m].shiftcult=0;

      }
    }
  }
  
  cellinfo();
  read_country_names();
  read_country_codes();
  read_regional_codes();
  read_regional_names();
  read_continent_codes();
  read_other_data();
  initialize_woodharvest_country_ratios(2005);
  if(!initialrun) {
    read_restart(); 
  }else{
    //    strcpy(new_path,PATH1); 
    strcpy(new_path,future_initial_state); 
    printf("read in initial data... %s\n",new_path);
    start[0]=start[1]=start[2]=0;
    status = nc_open(new_path, NC_NOWRITE, &ncid_state);
    check_err(status,__LINE__,__FILE__);
    
    status = nc_inq_varid (ncid_state, "gothr", &varidrst_othr);
    check_err(status,__LINE__,__FILE__);
    status = nc_get_vara_double(ncid_state, varidrst_othr, start, count, &newdata[0].v[0][0]);
    check_err(status,__LINE__,__FILE__);
    
    if (use_urban) {
      status = nc_inq_varid (ncid_state, "gurbn", &varidrst_urbn);
      check_err(status,__LINE__,__FILE__);
      status = nc_get_vara_double(ncid_state, varidrst_urbn, start, count, &newdata[0].u[0][0]);
      check_err(status,__LINE__,__FILE__);
    }
    
    status = nc_inq_varid (ncid_state, "gsecd", &varidrst_secd);
    check_err(status,__LINE__,__FILE__);
    status = nc_get_vara_double(ncid_state, varidrst_secd, start, count, &newdata[0].s[0][0]);
    check_err(status,__LINE__,__FILE__);
    
    status = nc_inq_varid (ncid_state, "gcrop", &varidrst_crop);
    check_err(status,__LINE__,__FILE__);
    status = nc_get_vara_double(ncid_state, varidrst_crop, start, count, &newdata[0].c[0][0]);
    check_err(status,__LINE__,__FILE__);
    
    status = nc_inq_varid (ncid_state, "gpast", &varidrst_past);
    check_err(status,__LINE__,__FILE__);
    status = nc_get_vara_double(ncid_state, varidrst_past, start, count, &newdata[0].p[0][0]);
    check_err(status,__LINE__,__FILE__);
    
    status = nc_inq_varid (ncid_state, "gssmb", &varidrst_smb);
    check_err(status,__LINE__,__FILE__);
    status = nc_get_vara_double(ncid_state, varidrst_smb, start, count, &newdata[0].smb[0][0]);
    check_err(status,__LINE__,__FILE__);
    
    status = nc_inq_varid (ncid_state, "gssma", &varidrst_sma);
    check_err(status,__LINE__,__FILE__);
    status = nc_get_vara_double(ncid_state, varidrst_sma, start, count, &newdata[0].sma[0][0]);
    check_err(status,__LINE__,__FILE__);
    
    status = nc_close(ncid_state);
  }
  
  return;
}
/*****************************************************************************/
void run_model(){

  // start time loop
  int curryear;
  for (curryear=start_year;curryear<stop_year;curryear++){
    step(&curryear);
  } 
  return;
}
/*****************************************************************************/
void stepglm_ccsm(int *year,double *glmi,int *glmi_fdim1, int *glmi_fdim2,double *glmi_wh,int *glmi_wh_fdim1, double *glmo, int *glmo_fdim1, int *glmo_fdim2){
  int k,m,i,iz;
  int country_code, continent_code, regional_code,curryear;
  char country_name[50], continent_name[50], regional_name[50];

  curryear=*year;
  printf("curryear=%ld start_year=%d stop_year=%d\n",curryear,start_year,stop_year);
  if (curryear >=start_year && curryear<=stop_year){

    if (initialrun && curryear==start_year) {
      outstat="w";
      print_file_headers=1;
    }else{
      outstat="a";
      print_file_headers=0;
    }

    read_data_water_ice_nc(curryear); 

    //    read_data_futureruns_nc(curryear); 
    //    unpack glmi and glmi_wh into appropriate arrays
    i=-1;
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){
	i+=1;
	newdata[1].v[k][m]=glmi[i];
	i+=1;
	newdata[1].c[k][m]=glmi[i];
	i+=1;
	newdata[1].p[k][m]=glmi[i];
      }
    }
    for (i=0;i<NREG;i++){
      rtdata[i].wh==glmi_wh[i];
    }

    transitions(curryear); 
    
#if 0
    for (i=0;i<NREG;i++){
      
      if(rdata[i].rcode > 0 && rdata[i].rcode < 19){
	strcpy(regional_name,rname[i]);
	regional_code=rdata[i].rcode;
	regional_timeseries_checker(regional_code,regional_name,curryear);
      }
    }
#endif

#if 1

    strcpy(regional_name,"global");
#if REGION_TEST
    regional_code=REGION_TEST_GCODE;
#else
    regional_code=1000;
#endif
    global_timeseries_checker(regional_code,regional_name,curryear); 
    
#endif


#if 0
    strcpy(continent_name,"na");
    continent_code=1;
    continent_timeseries_checker(continent_code,continent_name,curryear);
    
    strcpy(continent_name,"sa");
    continent_code=2;
    continent_timeseries_checker(continent_code,continent_name,curryear);
    
    strcpy(continent_name,"eu");
    continent_code=3;
    continent_timeseries_checker(continent_code,continent_name,curryear);
    
    strcpy(continent_name,"asia");
    continent_code=4;
    continent_timeseries_checker(continent_code,continent_name,curryear);
    
    strcpy(continent_name,"africa");
    continent_code=5;
    continent_timeseries_checker(continent_code,continent_name,curryear);
    
    strcpy(continent_name,"aus");
    continent_code=6;
    continent_timeseries_checker(continent_code,continent_name,curryear);
#endif

#if 0
    loop_call_for_country_final_stats(curryear);
#endif
    
    if (output_text) {
      output_updated_states(curryear);
      output_updated_states2(curryear);
      if (use_urban) output_updated_states3(curryear);
      output_lu(curryear); 
    }
    if (output_netcdf) {
      output_updated_states_nc(curryear);
      output_lu_nc(curryear);
    }
    
    country_primeflow_print(curryear);

  } 

  // Fill glmo output array

  int km=-1;
  for ( k=0;k<NY;k++){
    for ( m=0;m<NX;m++){
      km=km+1;
      //      printf("k=%d m=%d\n",k,m);
      //      if (newdata[0].c[k][m] != 0.) printf("filling glmo1.c[%ld] with %f\n kk=%d mm=%d",km,newdata[0].c[k][m],k,m);
      glmo[km] = newdata[1].c[k][m];
      km=km+1;
      glmo[km] = newdata[1].p[k][m];
      km=km+1;
      glmo[km] = newdata[1].v[k][m];
      km=km+1;
      glmo[km] = newdata[1].s[k][m];
      km=km+1;
      glmo[km] = newdata[1].flowvbh[k][m];
      km=km+1;
      glmo[km] = newdata[1].flowvbh2[k][m];
      km=km+1;
      //      if (newdata[1].flowsbh[k][m] != 0.) printf("filling glmo7.flowsbh[%ld] with %f kk=%d mm=%d",km,newdata[1].flowsbh[k][m],k,m);
      glmo[km] = newdata[1].flowsbh[k][m];
      km=km+1;
      glmo[km] = newdata[1].flowsbh2[k][m];
      km=km+1;
      glmo[km] = newdata[1].flowsbh3[k][m];
    }
  }

  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){   
      newdata[0].c[k][m]=newdata[1].c[k][m];
      newdata[0].p[k][m]=newdata[1].p[k][m];
      newdata[0].v[k][m]=newdata[1].v[k][m];
      //	    newdata[0].i[k][m]=newdata[1].i[k][m];
      //	    newdata[0].w[k][m]=newdata[1].w[k][m];
      newdata[0].s[k][m]=newdata[1].s[k][m];
      if (use_urban) newdata[0].u[k][m]=newdata[1].u[k][m];
      newdata[0].sma[k][m]=newdata[1].sma[k][m];
      newdata[0].smb[k][m]=newdata[1].smb[k][m];
      newdata[1].s[k][m]=0.;
      newdata[1].sma[k][m]=0.;
      newdata[1].smb[k][m]=0.;
      
      newdata[0].flowcp[k][m]=ZEROVALUE;
      newdata[0].flowpc[k][m]=ZEROVALUE;
      newdata[0].flowpv[k][m]=ZEROVALUE;
      newdata[0].flowvp[k][m]=ZEROVALUE;
      newdata[0].flowvc[k][m]=ZEROVALUE;
      newdata[0].flowcv[k][m]=ZEROVALUE;
      newdata[0].flowsc[k][m]=ZEROVALUE;
      newdata[0].flowcs[k][m]=ZEROVALUE;
      newdata[0].flowsp[k][m]=ZEROVALUE;
      newdata[0].flowps[k][m]=ZEROVALUE;
      newdata[0].flowvs[k][m]=ZEROVALUE;
      newdata[0].flowcu[k][m]=ZEROVALUE;
      newdata[0].flowpu[k][m]=ZEROVALUE;
      newdata[0].flowsu[k][m]=ZEROVALUE;
      newdata[0].flowvu[k][m]=ZEROVALUE;
      newdata[0].flowuc[k][m]=ZEROVALUE;
      newdata[0].flowup[k][m]=ZEROVALUE;
      newdata[0].flowus[k][m]=ZEROVALUE;
      newdata[0].sbh2[k][m]=ZEROVALUE;
      newdata[0].sbh3[k][m]=ZEROVALUE;
      newdata[0].vbh2[k][m]=ZEROVALUE;
      newdata[0].vbh[k][m]=ZEROVALUE;
      newdata[0].sbh[k][m]=ZEROVALUE; 
    } 
  }
  
  for (i=0;i<NREG;i++) {
    
    hist_reg_wh[i]=ZEROVALUE;
    wh_ratio[i]=ZEROVALUE;
    
    rtdata[i].wh=ZEROVALUE;
    rtdata[i].whr=ZEROVALUE;
    rtdata[i].stbh=ZEROVALUE;
    rtdata[i].vwh=ZEROVALUE;
    
    rtdata[i].wh_lessthan_zmax=ZEROVALUE;
    rtdata[i].wh_at_zmax=ZEROVALUE;
    rtdata[i].vnfb=ZEROVALUE; 
    rtdata[i].smb=ZEROVALUE;
    rtdata[i].fh_sbh2=ZEROVALUE;
    rtdata[i].fh_sbh3=ZEROVALUE;
    rtdata[i].fh_vbh2=ZEROVALUE;
    rtdata[i].sarea=ZEROVALUE;
    rtdata[i].sarea_nf=ZEROVALUE;
    rtdata[i].smb_nf=ZEROVALUE;
    rtdata[i].converted_forest_land=ZEROVALUE;
    rtdata[i].converted_forest_land_total=ZEROVALUE;
    rtdata[i].flowvc_prime=ZEROVALUE;
    rtdata[i].flowsc_prime=ZEROVALUE;
    rtdata[i].flowcs_prime=ZEROVALUE;
    rtdata[i].flowvp_prime=ZEROVALUE;
    rtdata[i].flowsp_prime=ZEROVALUE;
    rtdata[i].flowps_prime=ZEROVALUE;
    rtdata[i].flowus_prime=ZEROVALUE;
    rtdata[i].flowsu_prime=ZEROVALUE;
    rtdata[i].flowvu_prime=ZEROVALUE;
    rtdata[i].reg_unmet=ZEROVALUE;
    rtdata[i].reg_avail=ZEROVALUE;
    
    for (iz=0;iz<NZ;iz++){
      rztdata[i][iz].vb=ZEROVALUE;
    }
  }
  for (i=0;i<NCCODE;i++) {
    
    
    ctdata[i].wh=ZEROVALUE;
    ctdata[i].whr=ZEROVALUE;
    ctdata[i].stbh=ZEROVALUE;
    ctdata[i].vwh=ZEROVALUE;
    ctdata[i].wh1=ZEROVALUE;
    ctdata[i].orig_whd=ZEROVALUE;
    
    ctdata[i].wh_lessthan_zmax=ZEROVALUE;
    ctdata[i].wh_at_zmax=ZEROVALUE;
    ctdata[i].vnfb=ZEROVALUE; 
    ctdata[i].smb=ZEROVALUE;
    ctdata[i].fh_sbh2=ZEROVALUE;
    ctdata[i].fh_sbh3=ZEROVALUE;
    ctdata[i].fh_vbh2=ZEROVALUE;
    ctdata[i].sarea=ZEROVALUE;
    ctdata[i].sarea_nf=ZEROVALUE;
    ctdata[i].smb_nf=ZEROVALUE;
    ctdata[i].converted_forest_land=ZEROVALUE;
    ctdata[i].converted_forest_land_total=ZEROVALUE;
    ctdata[i].predict_b=ZEROVALUE;
    ctdata[i].flowvc_prime=ZEROVALUE;
    ctdata[i].flowsc_prime=ZEROVALUE;
    ctdata[i].flowcs_prime=ZEROVALUE;
    ctdata[i].flowvp_prime=ZEROVALUE;
    ctdata[i].flowsp_prime=ZEROVALUE;
    ctdata[i].flowps_prime=ZEROVALUE;
    
    for (iz=0;iz<NZ;iz++){
      cztdata[i][iz].vb=ZEROVALUE;
    }
  }
  return;
}

/*****************************************************************************/
void step(int *year){

  int country_code, continent_code, regional_code,curryear;
  char country_name[50], continent_name[50], regional_name[50];
  curryear=*year;
  printf("curryear=%ld start_year=%d stop_year=%d\n",curryear,start_year,stop_year);
  if (curryear >=start_year && curryear<=stop_year){

    if (initialrun && curryear==start_year) {
      outstat="w";
      print_file_headers=1;
    }else{
      outstat="a";
      print_file_headers=0;
    }

    read_data_futureruns_nc(curryear); 

    transitions(curryear); 
    
#if 0
    for (i=0;i<NREG;i++){
      
      if(rdata[i].rcode > 0 && rdata[i].rcode < 19){
	strcpy(regional_name,rname[i]);
	regional_code=rdata[i].rcode;
	regional_timeseries_checker(regional_code,regional_name,curryear);
      }
    }
#endif

#if 1

    strcpy(regional_name,"global");
#if REGION_TEST
    regional_code=REGION_TEST_GCODE;
#else
    regional_code=1000;
#endif
    global_timeseries_checker(regional_code,regional_name,curryear); 
    
#endif


#if 0
    strcpy(continent_name,"na");
    continent_code=1;
    continent_timeseries_checker(continent_code,continent_name,curryear);
    
    strcpy(continent_name,"sa");
    continent_code=2;
    continent_timeseries_checker(continent_code,continent_name,curryear);
    
    strcpy(continent_name,"eu");
    continent_code=3;
    continent_timeseries_checker(continent_code,continent_name,curryear);
    
    strcpy(continent_name,"asia");
    continent_code=4;
    continent_timeseries_checker(continent_code,continent_name,curryear);
    
    strcpy(continent_name,"africa");
    continent_code=5;
    continent_timeseries_checker(continent_code,continent_name,curryear);
    
    strcpy(continent_name,"aus");
    continent_code=6;
    continent_timeseries_checker(continent_code,continent_name,curryear);
#endif

#if 0
    loop_call_for_country_final_stats(curryear);
#endif
    
    if (output_text) {
      output_updated_states(curryear);
      output_updated_states2(curryear);
      if (use_urban) output_updated_states3(curryear);
      output_lu(curryear); 
    }
    if (output_netcdf) {
      output_updated_states_nc(curryear);
      output_lu_nc(curryear);
    }
    
    country_primeflow_print(curryear);

  } 
  return;
}

/***********************************************************************/
void finalize(void){
 if (output_netcdf) output_updated_states_nc(stop_year);
  create_restart_file();
  // closed here to save having to open and read up through dataset each
  // timestep - just used for backward compatability with tone,two... way
  // of doing things.
  return;
}
/***********************************************************************/
/********************************************************************/
void read_data_futureruns_nc(int curryear){

  char new_path[90];
  double first_year;
  int ic, i, k, m, it, cnum,j,iz;
  int it1;
  int ncid_past,
    ncid_othr,
    ncid_urbn,
    ncid_secd,
    ncid_suit,
    ncid_ice,
    ncid_water,
    ncid_crop,
    ncid_pasture;                          /* netCDF ID */
  int varid_past,
    varid_othr,
    varid_urbn,
    varid_secd,
    varid_suit,
    varid_ice,
    varid_water,
    varid_crop,
    varid_time,
    varid_pasture;                          /* netCDF ID */
  int status;                       /* error status */
  int varidrst_past,
    varidrst_othr,
    varidrst_urbn,
    varidrst_secd,
    varidrst_suit,
    varidrst_ice,
    varidrst_water,
    varidrst_crop,
    varidrst_sma,
    varidrst_smb,
    varidrst_pasture;                          /* netCDF ID */
  size_t count1[] = {1};
  size_t count[] = {1, NY, NX};
  size_t start1[] = {0}; /* start at first value */
  size_t start[] = {0, 0, 0}; /* start at first value */

  printf("\nreading data for future runs, 2005 and beyond...\n");
  strcpy(new_path,PATH1); 
  strcat(new_path,future_crop_path); 

  status = nc_open(new_path, NC_NOWRITE, &ncid_crop);
  check_err(status,__LINE__,__FILE__);
  status = nc_inq_varid (ncid_crop, "cropland", &varid_crop);
  check_err(status,__LINE__,__FILE__);
  status = nc_inq_varid (ncid_crop, "time", &varid_time);
  check_err(status,__LINE__,__FILE__);
  status = nc_get_vara_double(ncid_crop, varid_time, start1, count1, &first_year);
  check_err(status,__LINE__,__FILE__);

  strcpy(new_path,PATH1); 
  strcat(new_path,future_past_path); 

  status = nc_open(new_path, NC_NOWRITE, &ncid_past);
  check_err(status,__LINE__,__FILE__);
  status = nc_inq_varid (ncid_past, "pasture", &varid_past);
  check_err(status,__LINE__,__FILE__);

  strcpy(new_path,PATH1); 
  strcat(new_path,future_othr_path); 

  status = nc_open(new_path, NC_NOWRITE, &ncid_othr);
  check_err(status,__LINE__,__FILE__);
  status = nc_inq_varid (ncid_othr, "primaryland", &varid_othr);
  check_err(status,__LINE__,__FILE__);

  strcpy(new_path,PATH1); 
  strcat(new_path,future_icew_path); 

  status = nc_open(new_path, NC_NOWRITE, &ncid_ice);
  check_err(status,__LINE__,__FILE__);
  status = nc_inq_varid (ncid_ice, "ice", &varid_ice);
  check_err(status,__LINE__,__FILE__);

  strcpy(new_path,PATH1); 
  strcat(new_path,future_watr_path); 

  status = nc_open(new_path, NC_NOWRITE, &ncid_water);
  check_err(status,__LINE__,__FILE__);
  status = nc_inq_varid (ncid_water, "water", &varid_water);
  check_err(status,__LINE__,__FILE__);

  //  status = nc_open(future_suit_path, NC_NOWRITE, &ncid_suit);
  //  check_err(status,__LINE__,__FILE__);
  //  status = nc_inq_varid (ncid_suit, "suit", &varid_suit);
  //  check_err(status,__LINE__,__FILE__);

  //For each class type: if first time through prime the pipe with two years
  //otherwise copy the updated values from position 1 to position 0 and put
  //new data into position 1.  The copy should be replaced by an index that 
  // just alternates between 0 and 1 at the end of every timestep. 

  /* read in constructed future states */
  if (curryear==start_year) {
    for (it = 0;it<2;it++) {
      //start reading at current year index
      start[0]=curryear-first_year+it; 
      count[0]=1;
      status = nc_get_vara_double(ncid_crop, varid_crop, start, count, &newdata[it].c[0][0]);
      check_err(status,__LINE__,__FILE__);
      status = nc_get_vara_double(ncid_past, varid_past, start, count, &newdata[it].p[0][0]);
      check_err(status,__LINE__,__FILE__);
      status = nc_get_vara_double(ncid_othr, varid_othr, start, count, &newdata[it].v[0][0]);
      check_err(status,__LINE__,__FILE__);
      if (use_urban) {
	status = nc_get_vara_double(ncid_urbn, varid_urbn, start, count, &newdata[it].u[0][0]);
	check_err(status,__LINE__,__FILE__);
      }
      //for water and ice only read in first (only) time slice 
      start[0]=0;
      status = nc_get_vara_double(ncid_ice, varid_ice, start, count, &newdata[it].i[0][0]);
      check_err(status,__LINE__,__FILE__);
      status = nc_get_vara_double(ncid_water, varid_water, start, count, &newdata[it].w[0][0]);
      check_err(status,__LINE__,__FILE__);
      for (k=0;k<NY;k++){
	for (m=0;m<NX;m++){   
	  newdata[it].c[k][m]=fround(newdata[it].c[k][m],6);
	  newdata[it].p[k][m]=fround(newdata[it].p[k][m],6);
	  newdata[it].v[k][m]=fround(newdata[it].v[k][m],6);
	  if (use_urban) newdata[it].u[k][m]=fround(newdata[it].u[k][m],6);
	  newdata[0].i[k][m]=fround(newdata[0].i[k][m],6);
	  newdata[0].w[k][m]=fround(newdata[0].w[k][m],6);

      /* suitability based on distance, need zdis algorithm and zdis initialized to NZ */
	  if(suit_option == 1){
	    newdata[it].zdis[k][m]=NZ-1;
	  }
	}
      }
    }
  }else {
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){   
	newdata[0].c[k][m]=newdata[1].c[k][m];
	newdata[0].p[k][m]=newdata[1].p[k][m];
	newdata[0].v[k][m]=newdata[1].v[k][m];
	//	    newdata[0].i[k][m]=newdata[1].i[k][m];
	//	    newdata[0].w[k][m]=newdata[1].w[k][m];
	newdata[0].s[k][m]=newdata[1].s[k][m];
	if (use_urban) newdata[0].u[k][m]=newdata[1].u[k][m];
	newdata[0].sma[k][m]=newdata[1].sma[k][m];
	newdata[0].smb[k][m]=newdata[1].smb[k][m];
	newdata[1].s[k][m]=0.;
	newdata[1].sma[k][m]=0.;
	newdata[1].smb[k][m]=0.;
	
	newdata[0].flowcp[k][m]=ZEROVALUE;
	newdata[0].flowpc[k][m]=ZEROVALUE;
	newdata[0].flowpv[k][m]=ZEROVALUE;
	newdata[0].flowvp[k][m]=ZEROVALUE;
	newdata[0].flowvc[k][m]=ZEROVALUE;
	newdata[0].flowcv[k][m]=ZEROVALUE;
	newdata[0].flowsc[k][m]=ZEROVALUE;
	newdata[0].flowcs[k][m]=ZEROVALUE;
	newdata[0].flowsp[k][m]=ZEROVALUE;
	newdata[0].flowps[k][m]=ZEROVALUE;
	newdata[0].flowvs[k][m]=ZEROVALUE;
	newdata[0].flowcu[k][m]=ZEROVALUE;
	newdata[0].flowpu[k][m]=ZEROVALUE;
	newdata[0].flowsu[k][m]=ZEROVALUE;
	newdata[0].flowvu[k][m]=ZEROVALUE;
	newdata[0].flowuc[k][m]=ZEROVALUE;
	newdata[0].flowup[k][m]=ZEROVALUE;
	newdata[0].flowus[k][m]=ZEROVALUE;
	newdata[0].sbh2[k][m]=ZEROVALUE;
	newdata[0].sbh3[k][m]=ZEROVALUE;
	newdata[0].vbh2[k][m]=ZEROVALUE;
	newdata[0].vbh[k][m]=ZEROVALUE;
	newdata[0].sbh[k][m]=ZEROVALUE; 
      } 
    }
    start[0]=curryear-first_year+1;
    count[0]=1;
    status = nc_get_vara_double(ncid_crop, varid_crop, start, count, &newdata[1].c[0][0]);
    check_err(status,__LINE__,__FILE__);
    status = nc_get_vara_double(ncid_past, varid_past, start, count, &newdata[1].p[0][0]);
    check_err(status,__LINE__,__FILE__);
    status = nc_get_vara_double(ncid_othr, varid_othr, start, count, &newdata[1].v[0][0]);
    check_err(status,__LINE__,__FILE__);
    if (use_urban) {
      status = nc_get_vara_double(ncid_urbn, varid_urbn, start, count, &newdata[1].u[0][0]);
      check_err(status,__LINE__,__FILE__);
    }
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){   
	newdata[1].c[k][m]=fround(newdata[1].c[k][m],6);
	newdata[1].p[k][m]=fround(newdata[1].p[k][m],6);
	newdata[1].v[k][m]=fround(newdata[1].v[k][m],6);
	if (use_urban) newdata[1].u[k][m]=fround(newdata[1].u[k][m],6);
      }
    }
    //	 don't need to get dat for ice and water always same
  }
    
  //  if this is the 1st step of run then overwrite 
  //  time[0] with future initial state, time[1] filled in from boundary data above.
  
  //  if(!initialrun && start_year==curryear){
  if(start_year==curryear){
    //    strcpy(new_path,PATH1); 
    strcpy(new_path,future_initial_state); 
    printf("newdata values in read_data_past... %s\n",new_path);
    start[0]=start[1]=start[2]=0;
    status = nc_open(new_path, NC_NOWRITE, &ncid_state);
    check_err(status,__LINE__,__FILE__);
    
    status = nc_inq_varid (ncid_state, "gothr", &varidrst_othr);
    check_err(status,__LINE__,__FILE__);
    status = nc_get_vara_double(ncid_state, varidrst_othr, start, count, &newdata[0].v[0][0]);
    check_err(status,__LINE__,__FILE__);

    if (use_urban) {
      status = nc_inq_varid (ncid_state, "gurbn", &varidrst_urbn);
      check_err(status,__LINE__,__FILE__);
      status = nc_get_vara_double(ncid_state, varidrst_urbn, start, count, &newdata[0].u[0][0]);
      check_err(status,__LINE__,__FILE__);
    }
    
    status = nc_inq_varid (ncid_state, "gsecd", &varidrst_secd);
    check_err(status,__LINE__,__FILE__);
    status = nc_get_vara_double(ncid_state, varidrst_secd, start, count, &newdata[0].s[0][0]);
    check_err(status,__LINE__,__FILE__);

    status = nc_inq_varid (ncid_state, "gcrop", &varidrst_crop);
    check_err(status,__LINE__,__FILE__);
    status = nc_get_vara_double(ncid_state, varidrst_crop, start, count, &newdata[0].c[0][0]);
    check_err(status,__LINE__,__FILE__);

    status = nc_inq_varid (ncid_state, "gpast", &varidrst_past);
    check_err(status,__LINE__,__FILE__);
    status = nc_get_vara_double(ncid_state, varidrst_past, start, count, &newdata[0].p[0][0]);
    check_err(status,__LINE__,__FILE__);

    status = nc_inq_varid (ncid_state, "gssmb", &varidrst_smb);
    check_err(status,__LINE__,__FILE__);
    status = nc_get_vara_double(ncid_state, varidrst_smb, start, count, &newdata[0].smb[0][0]);
    check_err(status,__LINE__,__FILE__);

    status = nc_inq_varid (ncid_state, "gssma", &varidrst_sma);
    check_err(status,__LINE__,__FILE__);
    status = nc_get_vara_double(ncid_state, varidrst_sma, start, count, &newdata[0].sma[0][0]);
    check_err(status,__LINE__,__FILE__);

    status = nc_close(ncid_state);
  }
  
  if (0) {
  
  for (it1=0;it1<2;it1++){
    printf("newdata values in read_data_past...\n");
    printf("it1=%ld\n",it1);
    for (j=0;j<NY;j++){
      for (i=0;i<NX;i++){   
	if (newdata[it1].c[j][i] != 0.) printf("curryear,it1, j,i %ld  %ld %ld %ld %lf\n",curryear,it1,j,i,newdata[it1].c[j][i]);
	/*     printf("%lf\n",newdata[it1].p[j][i]); */
	/*     printf("%lf\n",newdata[it1].v[j][i]); */
	/*     printf("%lf\n",newdata[it1].i[j][i]); */
	/*     printf("%lf\n",newdata[it1].w[j][i]); */
	/*     printf("%lf\n",newdata[it1].s[j][i]); */
	/*     printf("%lf\n",newdata[it1].sma[j][i]); */
	/*     printf("%lf\n",newdata[it1].smb[j][i]); */
      }
    }
  }
  }
status = nc_close(ncid_past);
check_err(status,__LINE__,__FILE__);
status = nc_close(ncid_ice);
check_err(status,__LINE__,__FILE__);
status = nc_close(ncid_water);
check_err(status,__LINE__,__FILE__);
status = nc_close(ncid_othr);
check_err(status,__LINE__,__FILE__);
if (use_urban) {
  status = nc_close(ncid_othr);
  check_err(status,__LINE__,__FILE__);
}
status = nc_close(ncid_crop);
check_err(status,__LINE__,__FILE__);

for (i=0;i<NREG;i++) {
  
  hist_reg_wh[i]=ZEROVALUE;
  wh_ratio[i]=ZEROVALUE;
  
  rtdata[i].wh=ZEROVALUE;
  rtdata[i].whr=ZEROVALUE;
  rtdata[i].stbh=ZEROVALUE;
  rtdata[i].vwh=ZEROVALUE;
  
  rtdata[i].wh_lessthan_zmax=ZEROVALUE;
  rtdata[i].wh_at_zmax=ZEROVALUE;
  rtdata[i].vnfb=ZEROVALUE; 
  rtdata[i].smb=ZEROVALUE;
  rtdata[i].fh_sbh2=ZEROVALUE;
  rtdata[i].fh_sbh3=ZEROVALUE;
  rtdata[i].fh_vbh2=ZEROVALUE;
  rtdata[i].sarea=ZEROVALUE;
  rtdata[i].sarea_nf=ZEROVALUE;
  rtdata[i].smb_nf=ZEROVALUE;
  rtdata[i].converted_forest_land=ZEROVALUE;
  rtdata[i].converted_forest_land_total=ZEROVALUE;
  rtdata[i].flowvc_prime=ZEROVALUE;
  rtdata[i].flowsc_prime=ZEROVALUE;
  rtdata[i].flowcs_prime=ZEROVALUE;
  rtdata[i].flowvp_prime=ZEROVALUE;
  rtdata[i].flowsp_prime=ZEROVALUE;
  rtdata[i].flowps_prime=ZEROVALUE;
  rtdata[i].flowus_prime=ZEROVALUE;
  rtdata[i].flowsu_prime=ZEROVALUE;
  rtdata[i].flowvu_prime=ZEROVALUE;
  rtdata[i].reg_unmet=ZEROVALUE;
  rtdata[i].reg_avail=ZEROVALUE;
  
  for (iz=0;iz<NZ;iz++){
    rztdata[i][iz].vb=ZEROVALUE;
  }
 }
for (i=0;i<NCCODE;i++) {
  
  
  ctdata[i].wh=ZEROVALUE;
  ctdata[i].whr=ZEROVALUE;
  ctdata[i].stbh=ZEROVALUE;
  ctdata[i].vwh=ZEROVALUE;
  ctdata[i].wh1=ZEROVALUE;
  ctdata[i].orig_whd=ZEROVALUE;
  
  ctdata[i].wh_lessthan_zmax=ZEROVALUE;
  ctdata[i].wh_at_zmax=ZEROVALUE;
  ctdata[i].vnfb=ZEROVALUE; 
  ctdata[i].smb=ZEROVALUE;
  ctdata[i].fh_sbh2=ZEROVALUE;
  ctdata[i].fh_sbh3=ZEROVALUE;
  ctdata[i].fh_vbh2=ZEROVALUE;
  ctdata[i].sarea=ZEROVALUE;
  ctdata[i].sarea_nf=ZEROVALUE;
  ctdata[i].smb_nf=ZEROVALUE;
  ctdata[i].converted_forest_land=ZEROVALUE;
  ctdata[i].converted_forest_land_total=ZEROVALUE;
  ctdata[i].predict_b=ZEROVALUE;
  ctdata[i].flowvc_prime=ZEROVALUE;
  ctdata[i].flowsc_prime=ZEROVALUE;
  ctdata[i].flowcs_prime=ZEROVALUE;
  ctdata[i].flowvp_prime=ZEROVALUE;
  ctdata[i].flowsp_prime=ZEROVALUE;
  ctdata[i].flowps_prime=ZEROVALUE;
  
  for (iz=0;iz<NZ;iz++){
    cztdata[i][iz].vb=ZEROVALUE;
  }
 }
    
return;
}
/********************************************************************/
void read_data_water_ice_nc(int curryear){

  char new_path[90];
  double first_year;
  int ic, i, k, m, it, cnum,j,iz;
  int it1;
  int ncid_ice,
    ncid_water;

  int varid_ice,
    varid_water;

  int status;                       /* error status */
  int varidrst_ice,
    varidrst_water;

  size_t count[] = {1, NY, NX};
  size_t start[] = {0, 0, 0}; /* start at first value */

  printf("\nreading data for water and ice...\n");

  strcpy(new_path,PATH1); 
  strcat(new_path,future_icew_path); 
  
  status = nc_open(new_path, NC_NOWRITE, &ncid_ice);
  check_err(status,__LINE__,__FILE__);
  status = nc_inq_varid (ncid_ice, "ice", &varid_ice);
  check_err(status,__LINE__,__FILE__);

  strcpy(new_path,PATH1); 
  strcat(new_path,future_watr_path); 

  status = nc_open(new_path, NC_NOWRITE, &ncid_water);
  check_err(status,__LINE__,__FILE__);
  status = nc_inq_varid (ncid_water, "water", &varid_water);
  check_err(status,__LINE__,__FILE__);

  /* read in constructed future states */
  if (curryear==start_year) {
    for (it = 0;it<2;it++) {
      //start reading at current year index
      count[0]=1;
      start[0]=0;
      status = nc_get_vara_double(ncid_ice, varid_ice, start, count, &newdata[it].i[0][0]);
      check_err(status,__LINE__,__FILE__);
      status = nc_get_vara_double(ncid_water, varid_water, start, count, &newdata[it].w[0][0]);
      check_err(status,__LINE__,__FILE__);
      for (k=0;k<NY;k++){
	for (m=0;m<NX;m++){   
	  newdata[0].i[k][m]=fround(newdata[0].i[k][m],6);
	  newdata[0].w[k][m]=fround(newdata[0].w[k][m],6);
	}
      }
    }
  }
  
  status = nc_close(ncid_ice);
  check_err(status,__LINE__,__FILE__);
  status = nc_close(ncid_water);
  check_err(status,__LINE__,__FILE__);
  
  return;
}
/********************************************************************/
void
check_err(const int stat, const int line, const char *file) {
    if (stat != NC_NOERR) {
           (void) fprintf(stderr, "line %d of %s: %s\n", line, file, 
                          nc_strerror(stat));
        exit(1);
    }
}
/********************************************************************/
void cellinfo(){

  FILE *infile;
  //  char new_path[90], tmp_path0[90], tmp_path1[90], tmp_path4[90], new_path2[90], in_cat[90];
  char new_path[190];
  int k, m, it;
  float er=6370., pi=3.141592654;
  float nwloncorner=-180., nwlatcorner=90., dtocenter, dstep;
  double summer=ZEROVALUE, summer2=ZEROVALUE;
 
  if (res_option==2){
    dtocenter=-0.25;
    dstep=-0.5;
  }

  lat[0]=nwlatcorner+dtocenter;
  for (k=1;k<NY;k++){
    lat[k]=lat[k-1]+dstep;
  }
  
 
  if (res_option==2){
    dtocenter=0.25;
    dstep=0.5;
  }    
 
  lon[0]=nwloncorner+dtocenter;
  for (m=1;m<NX;m++){
     lon[m]=lon[m-1]+dstep;
  }
  
  
/*    strcpy(new_path,PATH1); */
/*    strcpy(in_cat,CAT3);  */
/*    strcpy(tmp_path0,new_path); */
/*    strcat(tmp_path0,in_cat); */
/*    strcpy(tmp_path1,tmp_path0); */
/*    strcpy(tmp_path4,"cellarea/"); */
/*    strcat(tmp_path1,tmp_path4); */
/*    strcpy(new_path2,tmp_path1); */

/*    if (res_option==2){ */
/*      strcat(new_path2,"cellarea_halfdeg.txt"); */
/*    }  */
  strcpy(new_path,PATH1);
  strcat(new_path,cellinfo_file);

  infile=fopen(new_path,"r");
  
  
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
  
  
  /* strcpy(new_path,PATH1); */
  /* strcpy(tmp_path0,new_path); */
  /* strcat(tmp_path0,CAT3); */
  /* strcpy(tmp_path1,tmp_path0); */

  /* if ((nodata_option==5)||(nodata_option==6)) { */
  /*   strcat(tmp_path1,"ccodes/ccodes.txt.sort2wh"); */
  /* } */
  /* else { */
  /*   strcat(tmp_path1,"ccodes/ccodes.txt.sort2wh_139countries"); */
  /* } */
 
  /* strcpy(ffname,tmp_path1); */
  strcpy(ffname,PATH1);
  strcat(ffname,ccodes_file);
  
  cinfile=fopen(ffname,"r");
  for (i=0;i<NCCODE;i++) fscanf(cinfile,"%d\n",&cdata[i].ccode);
  fclose(cinfile);
  

  if (res_option==2){
    /* strcpy(new_path,PATH1); */
    /* strcpy(tmp_path0,new_path); */
    /* strcat(tmp_path0,CAT3); */
    /* strcpy(tmp_path1,tmp_path0); */
    /* strcat(tmp_path1,"shift_cult/shiftcult_map_halfdeg.txt"); */
    /* strcpy(scfname,tmp_path1); */
  strcpy(scfname,PATH1);
  strcat(scfname,shiftcult_map);
  
    cinfile=fopen(scfname,"r");
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){
	fscanf(cinfile,"%d\n",&dstatic[k][m].shiftcult);
      }
    }
    fclose(cinfile);
  }

  
  /* strcpy(new_path,PATH1); */
  /* strcpy(tmp_path0,new_path); */
  /* strcat(tmp_path0,CAT3); */
  /* strcpy(tmp_path1,tmp_path0); */


  /* if (res_option==1){ */
  /*   if ((nodata_option==5)||(nodata_option==6)){ */
  /*     strcat(tmp_path1,"ccodes/gcodes_ak_conform.asc"); */
  /*   } */
  /*   else { */
  /*     strcat(tmp_path1,"ccodes/gcodes_ak_conform.asc"); */
  /*   } */
  /* } */
  /* else if (res_option==2){ */
  /*   strcat(tmp_path1,"ccodes/ccodes_half_deg.txt"); */
  /* }  */


  /* strcpy(gfname,tmp_path1); */

  strcpy(tmp_path1,PATH1);
  strcat(tmp_path1,ccodes_map);
  cinfile=fopen(tmp_path1,"r");  


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
  }




  
 
  return;
}
/********************************************************************/
void read_regional_codes(){

  FILE *cinfile;
  FILE *outfile; 
  int i, k, m, it;
  int tmpcode[NCCODE];
  char new_path[90], tmp_path0[90], tmp_path1[90], newpath1[90], newpath[90];
  char gfname[90], ffname[90], outfname[90], fn[90],fout[90];

  printf("\nreading regional codes...\n");

 /* strcpy(new_path,PATH1); */
 /*  strcpy(tmp_path0,new_path); */
 /*  strcat(tmp_path0,CAT3); */
 /*  strcpy(tmp_path1,tmp_path0); */
 /*  if(res_option==1) strcat(tmp_path1,"wood_harvest/codes.txt"); */
 /*   else if(res_option==2) strcat(tmp_path1,"wood_harvest/codes_halfdeg_minicam.txt");  */
 /*  /\* else if(res_option==2) strcat(tmp_path1,"wood_harvest/codes_halfdeg.txt"); *\/ */

 /*  strcpy(ffname,tmp_path1); */
  
  strcpy(new_path,PATH1);
  strcat(new_path,whcodes_file);

  cinfile=fopen(new_path,"r");
  for (i=0;i<NREG;i++) fscanf(cinfile,"%d\n",&rdata[i].rcode);
  fclose(cinfile);

 /* strcpy(new_path,PATH1); */
 /*  strcpy(tmp_path0,new_path); */
 /*  strcat(tmp_path0,CAT3); */
 /*  strcpy(tmp_path1,tmp_path0); */

 /*    strcat(tmp_path1,"wood_harvest/continent_codes_minicam_test.txt");   */
 /*      /\*  strcat(tmp_path1,"wood_harvest/continent_codes24.txt"); *\/ */

 /*  strcpy(ffname,tmp_path1); */
  
  strcpy(new_path,PATH1);
  strcat(new_path,whcontcodes_file);

  cinfile=fopen(new_path,"r");
  for (i=0;i<NREG;i++) fscanf(cinfile,"%d\n",&rdata[i].continent_code);
  fclose(cinfile);


  /* strcpy(new_path,PATH1); */
  /* strcpy(tmp_path0,new_path); */
  /* strcat(tmp_path0,CAT3); */
  /* strcpy(tmp_path1,tmp_path0); */
  /* if(res_option==1) strcat(tmp_path1,"wood_harvest/codes2glm.txt"); */
  /* else if(res_option==2) strcat(tmp_path1,"wood_harvest/codes2glm_minicam.txt");  */
  /*     /\*   else if(res_option==2) strcat(tmp_path1,"wood_harvest/codes2glm_halfdeg_new3.txt");  *\/  */
  /* strcpy(ffname,tmp_path1); */


  strcpy(new_path,PATH1);
  strcat(new_path,whcodes2glm_map);
  

  cinfile=fopen(new_path,"r");
  for (i=0;i<NCCODE;i++) fscanf(cinfile,"%d %d\n",&cdata[i].rcode,&tmpcode[i]);
  fclose(cinfile);
  

  
  for (i=0;i<NCCODE;i++) {
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){
	
	if(dstatic[k][m].gcode == tmpcode[i]) dstatic[k][m].rcode=cdata[i].rcode;

	
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
  if (use_urban) {
    outfile=fopen(static_regions_file,"w");
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){
	fprintf(outfile,"%d ",dstatic[k][m].rcode);
      }
      fprintf(outfile,"\n");
    }
    fclose(outfile);
  }
 
  for (i=0;i<NCCODE;i++) {
    cdata[i].smart_flow_option = rdata[cdata[i].rcode-1].smart_flow_option;
    cdata[i].zdis_option = rdata[cdata[i].rcode-1].zdis_option;
    cdata[i].converted_forest_land_option = rdata[cdata[i].rcode-1].converted_forest_land_option;
    cdata[i].harvest_option = rdata[cdata[i].rcode-1].harvest_option;

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
  
  
  /* strcpy(new_path,PATH1); */
  /* strcpy(tmp_path0,new_path); */
  /* strcat(tmp_path0,CAT3); */
  /* strcpy(tmp_path1,tmp_path0); */
  /* strcat(tmp_path1,"ccodes/continent.codes.txt.sort2wh"); */
  /* strcpy(ffname,tmp_path1); */
  
  strcpy(ffname,PATH1); 
  strcat(ffname,contcodes_file); 
  cinfile=fopen(ffname,"r");
  for (i=0;i<NCCODE;i++) fscanf(cinfile,"%d\n",&cdata[i].continent_code);
  fclose(cinfile);
  
  
  strcpy(new_path,PATH1);
  /* strcpy(tmp_path0,new_path); */
  /* strcat(tmp_path0,CAT3); */
  /* strcpy(tmp_path1,tmp_path0); */

  /* if (res_option==1){ */
  /*   strcat(tmp_path1,"ccodes/gcodes_continent_conform.asc"); */
  /* } */
  /* else if (res_option==2){ */
  /*   strcat(tmp_path1,"ccodes/gcodes_continent_half_deg_DUMMY.asc"); */
  /* }  */

  /* strcpy(gfname,tmp_path1); */
  strcpy(new_path,PATH1);
  strcat(new_path,gcodes_cont_map);

  cinfile=fopen(new_path,"r");


  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){
      fscanf(cinfile,"%d\n",&dstatic[k][m].continent_code);
    }
  }
  fclose(cinfile);
  
  return;
}
/********************************************************************/
void read_other_data(){

  FILE *infile, *testfile, *outfile;
  int i, k, m, it, hist_t, reg_num;
/*   double sumtest=ZEROVALUE, summer, wh_multiplier=1.0; */
  double wh_multiplier=1.0;
/*   double tester=ZEROVALUE, junk;  */
  double junk; 
/*   char new_path[130], tmp_path0[130], tmp_path1[130], ffname[130], newpath[90], newpath1[90], fn[90], fout[90], tmppath[130], ftag1[7], ftag2[90],fin[90],ayear[30]; */
  char new_path[130];
/*   char new_path[130], tmp_path0[130], tmp_path1[130], ffname[130], newpath[90], newpath1[90], fn[90], fout[90], tmppath[130], ftag1[7], ftag2[90],fin[90],ayear[30]; */


  printf("\nreading other data...\n");

  /* static biomass grid, initial units=kgC/m2 */
  
/*   strcpy(new_path,PATH1); */
/*   strcpy(tmp_path0,new_path); */
/*   strcat(tmp_path0,CAT3); */
/*   strcpy(tmp_path1,tmp_path0);  */

/*   if (res_option==1){ */
/*     strcat(tmp_path1,"miami_biomass_v3/miami_1deg_conform.txt"); */
/*   } */
/*   else if (res_option==2){ */
/*     strcat(tmp_path1,"miami_biomass_v3/miami_halfdeg_conform.txt"); */
/*   }  */
  
/*   strcpy(ffname,tmp_path1); */
  
/*  infile=fopen(ffname,"r"); */
  strcpy(new_path,PATH1);
  strcat(new_path,miami_biomass_file_vba);

  infile=fopen(new_path,"r");


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


  if (use_urban) {
    outfile=fopen(static_vba_file,"w");
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){
	fprintf(outfile,"%d ",dstatic[k][m].vba);
      }
      fprintf(outfile,"\n");
    }
    fclose(outfile);
  }



  /* grid of vnpp from miami model in columnar format, initial units=kgC/m2/yr */

/*   strcpy(new_path,PATH1); */
/*   strcpy(tmp_path0,new_path); */
/*   strcat(tmp_path0,CAT3); */
/*   strcpy(tmp_path1,tmp_path0); */

/*   if (res_option==1){ */
/*     strcat(tmp_path1,"miami_npp/miami.dat.in_conform"); */
/*   } */
/*   else if (res_option==2){ */
/*     strcat(tmp_path1,"miami_npp/miami.half_deg.in_conform"); */
/*   }  */

/*   strcpy(ffname,tmp_path1); */
  
/*   infile=fopen(ffname,"r"); */
  strcpy(new_path,PATH1);
  strcat(new_path,miami_biomass_file_vnppa);

  infile=fopen(new_path,"r");
 
  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){
      fscanf(infile,"%lf\n",&dstatic[k][m].vnppa);
    }
  }
  fclose(infile);

  /* read probability of harvest biomass for secondary */
  
  strcpy(new_path,PATH1);
  strcat(new_path,PROB_FNAME);


  infile=fopen(new_path,"r");

/*   strcpy(ffname,PROB_FNAME); */
/*   infile=fopen(ffname,"r"); */
  
  for (i=0;i<NPHB;i++) fscanf(infile,"%*s %lf\n",&phb[i]);      
  fclose(infile);
  

  return;
}


/********************************************************************/
/* void read_woodharvest_data(int curryear){ */

/*   FILE *infile; */
/*   int i, k, m, it,yr,baseyear,hist_t,reg_num; */
/*   double wh_multiplier=1.0,junk; */
/*   char new_path[90], tmp_path0[90], tmp_path1[90], ffname[90], newpath[90], newpath1[90], tmppath[130], ftag1[7], ftag2[90],fin[90],ayear[30]; */



/*   if(gridded_wh ==0){ */

/*   /\* wood harvest by country thru time, initial units=MgC *\/ */
  
/*     strcpy(new_path,PATH1); */
/*     strcpy(tmp_path0,new_path); */
/*     strcat(tmp_path0,CAT3); */
/*     strcpy(tmp_path1,tmp_path0); */
    
/*     if((strcmp("tsix",trun) == 0) || (strcmp("tseven",trun) == 0)){ */
      
/*       if(future_rate_option == 1) {  */
/* 	strcat(tmp_path1,"wood_harvest/future/image/a1b/woodharvest_image.txt.warpf90.a1b."); */
/*       } */
/*       else if(future_rate_option == 2) {  */
/* 	strcat(tmp_path1,"wood_harvest/future/image/b1/woodharvest_image.txt.warpf90.b1."); */
/*       } */
/*       else if(future_rate_option == 3) { */
/* 	if((logging_option==1)||(logging_option==0)){ */
/* 	  strcat(tmp_path1,"wood_harvest/RCP/rcp_wh_minicam4."); */
/* 	} */
/* 	else if(logging_option==4){ */
/* 	  strcat(tmp_path1,"wood_harvest/RCP/rcp_wh_minicam_nd."); */
/* 	}  */
/*       } */
/*       else if(future_rate_option == 4) { */
/* 	if((logging_option==1)||(logging_option==0)){ */
/* 	  strcat(tmp_path1,"wood_harvest/RCP/rcp_wh_aim4."); */
/* 	} */
/* 	else if(logging_option==4){ */
/* 	  strcat(tmp_path1,"wood_harvest/RCP/rcp_wh_aim_nd."); */
/* 	} */
	
/*       } */
/*       else if(future_rate_option == 5) { */
/* 	if((logging_option==1)||(logging_option==0)){ */
/* 	  strcat(tmp_path1,"wood_harvest/RCP/rcp_wh_image4."); */
/* 	} */
/* 	else if(logging_option==4){ */
/* 	  strcat(tmp_path1,"wood_harvest/RCP/rcp_wh_image_nd."); */
/* 	} */
	
/*       } */
/*       else if(future_rate_option == 6) { */
/* 	if((logging_option==1)||(logging_option==0)){ */
/* 	  strcat(tmp_path1,"wood_harvest/RCP/rcp_wh_iiasa4."); */
/* 	} */
/* 	else if(logging_option==4){ */
/* 	  strcat(tmp_path1,"wood_harvest/RCP/rcp_wh_iiasa_nd."); */
/* 	} */
	
/*       } */
      
/*     } */
/*     else{ */
      
/*       if(logging_option == 4) { */
/* 	strcat(tmp_path1,"wood_harvest/nodata/woodharvin.txt.in.nodata."); */
/*       } */
/*       else if(logging_option == 1) { */
/* 	strcat(tmp_path1,"wood_harvest/1700-2000/woodharvin.txt.in."); */
/*       } */
/*       else if(logging_option == 2){ */
/* 	strcat(tmp_path1,"wood_harvest/1700-2000/woodharvin.txt.in."); */
/* 	wh_multiplier=1.5; */
/*       } */
/*       else{  /\* logging_option == 3 or 0 *\/ */
/* 	strcat(tmp_path1,"wood_harvest/1700-2000/woodharvin.txt.in.");  */
/* 	wh_multiplier=0.5; */
/*       } */
      
/*     } */
/*     strcpy(ffname,tmp_path1); */
/*     strcat(ffname,trun); */
    
/*     infile=fopen(ffname,"r"); */
/*     for (it=0;it<NEWNT;it++){ */
/*       for (i=0;i<NREG;i++) { */
	
/* 	fscanf(infile,"%lf\n",&rtdata[i].wh);     */
	
/* 	if(logging_option == 0) rtdata[i].wh=ZEROVALUE; */
	
/*       } */
/*     } */
/*     fclose(infile); */
    
/*   } */
  
/*   /\* read in future gridded WH data *\/ */
/*   if(gridded_wh==1){ */
/*     if((logging_option==1)|(logging_option==4)){ */
/*       strcpy(newpath,PATH1); */
      
/*       for (it=0;it<NEWNT;it++){ */
	
/* 	strcpy(tmppath,newpath); */
	
/* 	if(res_option==2) strcat(tmppath,CAT6); */
	
/* 	if(logging_option==4) strcat(tmppath,"nodata/"); */
	
/* 	strcpy(ftag1,"gfwhd."); */
	
/* 	sprintf(ayear,"%d.txt",curryear); */
/* 	strcat(ftag1,ayear);       */
/* 	strcpy(ftag2,ftag1); */
/* 	strcpy(fin,tmppath); */
/* 	strcat(fin,ftag2); */
     	
/* 	infile=fopen(fin,"r"); */
	
/* #if 1 */
/* 	for (i=0;i<HEADLEN;i++) fgets(header[i],27,infile); */
	
/* #endif */
	
	
/* 	for (k=0;k<NY;k++){ */
/* 	  for (m=0;m<NX;m++){ */
	    
/* 	    fscanf(infile,"%lf\n",&newdata[it].whd[k][m]); */
	  
/* 	  } */
/* 	} */
/* 	fclose(infile); */
	
/*       } /\* end it *\/ */

/*     } /\* logging_option if *\/ */
/*   } /\* gridded_wh if *\/ */
  
  
/*   if(gridded_wh == 0){ */
/*     /\* compute the ratio of each country's WH in 2005 to the region's WH in 2005 *\/ */
/*     strcpy(new_path,PATH1); */
/*     strcpy(tmp_path0,new_path); */
/*     strcat(tmp_path0,CAT3); */
/*     strcpy(tmp_path1,tmp_path0); */
    
/*     if(hist_option==1){ */
/*       hist_t = 62; */
/*       if(logging_option == 1) { */
/* 	strcat(tmp_path1,"wood_harvest/1700-2005/woodharvin.txt.in."); */
/*       } */
/*       else if(logging_option == 4){ */
/* 	strcat(tmp_path1,"wood_harvest/nodata/1700-2005/woodharvin.txt.in.nodata.");  */
/*       } */
/*       else if(logging_option == 0){ */
/* 	strcat(tmp_path1,"wood_harvest/1700-2005/woodharvin.txt.in.");  */
/* 	wh_multiplier=0; */
/*       } */
/*     }  */
/*     else if(hist_option==2){ */
/*       hist_t = 102; */
/*       if(logging_option == 1) { */
/* 	strcat(tmp_path1,"wood_harvest/1500-2005/woodharvin.txt.in."); */
/*       } */
/*       else if(logging_option == 4){ */
/* 	strcat(tmp_path1,"wood_harvest/nodata/1500-2005/woodharvin.txt.in.nodata.");  */
/*       } */
/*       else if(logging_option == 0){ */
/* 	strcat(tmp_path1,"wood_harvest/1500-2005/woodharvin.txt.in.");  */
/* 	wh_multiplier=0; */
/*       } */
      
/*     } */
/*     else if ((hist_option==3) || (hist_option==4 && use_urban)){ */
/*       hist_t = 32; */
/*       if(logging_option == 1) { */
/* 	strcat(tmp_path1,"wood_harvest/1850-2005/woodharvin.txt.in."); */
/*       } */
/*       else if(logging_option == 4){ */
/* 	strcat(tmp_path1,"wood_harvest/nodata/1850-2005/woodharvin.txt.in.");  */
/*       } */
/*       else if(logging_option == 0){ */
/* 	strcat(tmp_path1,"wood_harvest/1850-2005/woodharvin.txt.in.");  */
/* 	wh_multiplier=0; */
/*       } */
/*     } */
    
    
/*     strcpy(ffname,tmp_path1); */
/*     strcat(ffname,"tfive"); */
    
/*     infile=fopen(ffname,"r"); */
    
/*     for (it=0;it<hist_t;it++){ */
/*       for (i=0;i<NCCODE;i++){ */
/* 	if(it==hist_t-1){  */
/* 	  fscanf(infile,"%lf\n",&cdata[i].wh2005);      */
/* 	  cdata[i].wh2005*=TB2BB;  */
/* 	  if(logging_option == 0) cdata[i].wh2005=ZEROVALUE; */
/* 	} */
/* 	else{ */
/* 	  fscanf(infile,"%lf\n",&junk);  */
/* 	} */
/*       } */
      
/*     } */
/*     fclose(infile); */
    
    
/*     for (i=0;i<NCCODE;i++){ */
/*       reg_num = cdata[i].rcode; */
/*       hist_reg_wh[reg_num]+=cdata[i].wh2005; */
/*     } */
    
    
/*     for (i=0;i<NCCODE;i++){ */
/*       if(hist_reg_wh[cdata[i].rcode]<ZEROVALUE_CHECK){ */
/* 	wh_ratio[i] = ZEROVALUE; */
/*       } */
/*       else{ */
/* 	wh_ratio[i] = cdata[i].wh2005/hist_reg_wh[cdata[i].rcode]; */
/*       } */
      
/*     } */

/*   } /\* end gridded wh if *\/  */
  
/*   return; */
/* } */


/********************************************************************/
void read_woodharvest_data_nc(int curryear){
  //need to just malloc ctdatawh after getting country code size
  double first_year,rtdatawh[NREG];
  size_t rnum;
  int wh_multiplier=1.0;
  int ncid,i,region_code;            /* netCDF ID */
  int status;                       /* error status */
  int varid_time,varid_wh,dimid_region;        /* variable ID */
  size_t count1[] = {1};
  size_t start1[] = {0}; /* start at first value */
  size_t count2[] = {1,0};
  size_t start2[] = {0,0}; /* start at first value */

  printf("\n reading woodharvest data for year %d...\n",curryear);

  /* wood harvest projections */

  if(logging_option == 4) {
    status = nc_open(rcp_woodharvest_nodata_file, NC_NOWRITE, &ncid);
  }else{
    status = nc_open(rcp_woodharvest_file, NC_NOWRITE, &ncid);
  }
  check_err(status,__LINE__,__FILE__);

  status = nc_inq_varid (ncid, "time", &varid_time);
  check_err(status,__LINE__,__FILE__);

  status = nc_get_vara_double(ncid, varid_time, start1, count1, &first_year);
  check_err(status,__LINE__,__FILE__);

  status = nc_inq_dimid (ncid, "region_code", &dimid_region);
  check_err(status,__LINE__,__FILE__);
  status = nc_inq_dimlen(ncid, dimid_region, &rnum);
  check_err(status,__LINE__,__FILE__);


  status = nc_inq_varid (ncid, "woodharvest", &varid_wh);
  check_err(status,__LINE__,__FILE__);

  start2[0]=curryear-first_year;
  start2[1]=0;
  count2[0]=1;
  count2[1]=NREG;
  status = nc_get_vara_double(ncid, varid_wh, start2, count2, &rtdatawh);
  check_err(status,__LINE__,__FILE__);

  for (i=0;i<NREG;i++) {
    if(logging_option == 0) {
      rtdata[i].wh=ZEROVALUE;
    }else{
      rtdata[i].wh=rtdatawh[i];
    }
  }
  status = nc_close(ncid);
  check_err(status,__LINE__,__FILE__);

  return;
}
/********************************************************************/
void initialize_woodharvest_country_ratios(int baseyear){
  //need to just malloc ctdatawh after getting country code size
  double first_year,ctdatawh[NCCODE];
  size_t cnum;
  int reg_num; 
  int wh_multiplier=1.0;
  int ncid,i;                       /* netCDF ID */
  int status;                       /* error status */
  int varid_time,varid_wh,dimid_country;        /* variable ID */
  size_t count1[] = {1};
  size_t start1[] = {0}; /* start at first value */
  size_t count2[] = {1,0};
  size_t start2[] = {0,0}; /* start at first value */

  printf("\n reading woodharvest data for year %d...\n",baseyear);

  /* wood harvest by country thru time, initial units=MgC */

  if(logging_option == 4) {
    status = nc_open(woodharvest_nodata_file, NC_NOWRITE, &ncid);
  }else{
    status = nc_open(woodharvest_file, NC_NOWRITE, &ncid);
  }
  check_err(status,__LINE__,__FILE__);

  status = nc_inq_varid (ncid, "time", &varid_time);
  check_err(status,__LINE__,__FILE__);

  status = nc_get_vara_double(ncid, varid_time, start1, count1, &first_year);
  check_err(status,__LINE__,__FILE__);

  status = nc_inq_varid (ncid, "woodharvest", &varid_wh);
  check_err(status,__LINE__,__FILE__);

  start2[0]=baseyear-first_year;
  start2[1]=0;
  count2[0]=1;
  count2[1]=cnum;
  status = nc_get_vara_double(ncid, varid_wh, start2, count2, &ctdatawh[0]);
  check_err(status,__LINE__,__FILE__);

  if(logging_option == 0) wh_multiplier=0;

  for (i=0;i<NCCODE;i++) {
    ctdata[i].wh=ctdatawh[i]*TB2BB*wh_multiplier;
  }

  status = nc_close(ncid);
  check_err(status,__LINE__,__FILE__);

    
  for (i=0;i<NCCODE;i++){
    reg_num = cdata[i].rcode;
    hist_reg_wh[reg_num]+=cdata[i].wh2005;
  }
    
    
  for (i=0;i<NCCODE;i++){
    if(hist_reg_wh[cdata[i].rcode]<ZEROVALUE_CHECK){
      wh_ratio[i] = ZEROVALUE;
    }
    else{
      wh_ratio[i] = cdata[i].wh2005/hist_reg_wh[cdata[i].rcode];
    }
    
  }
  return;
}
/********************************************************************/
 
void read_country_names(){

  FILE *cinfile;
  int i;
  char new_path[90], tmp_path0[90], tmp_path1[90], ffname[90];
  

  printf("\nreading country names...\n");
 
  /* strcpy(new_path,PATH1); */
  /* strcpy(tmp_path0,new_path); */
  /* strcat(tmp_path0,CAT3); */
  /* strcpy(tmp_path1,tmp_path0); */
  /* strcat(tmp_path1,"ccodes/cnames.txt.sort2wh"); */
  /* strcpy(ffname,tmp_path1); */
  
  strcpy(new_path,PATH1);
  strcat(new_path,cnames_file);

  cinfile=fopen(new_path,"r");
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
 
  /* strcpy(new_path,PATH1); */
  /* strcpy(tmp_path0,new_path); */
  /* strcat(tmp_path0,CAT3); */
  /* strcpy(tmp_path1,tmp_path0); */
  /* /\* strcat(tmp_path1,"wood_harvest/names24.txt"); *\/ */
  /*     strcat(tmp_path1,"wood_harvest/names_minicam.txt");  */
  /*    strcpy(ffname,tmp_path1); */
  
  strcpy(new_path,PATH1);
  strcat(new_path,regnames_file);

  cinfile=fopen(new_path,"r");
  for (i=0;i<NREG;i++) {
    fscanf(cinfile,"%s\n",&rname[i]);
  }
  fclose(cinfile);

 

  return;

}


/********************************************************************/
void country_primeflow_print(int curryear){

  FILE *testfile;  
  int i, it;
 

  testfile=fopen("regional.primeflow.txt","w");


    for (i=0;i<NREG;i++){
     
      if (use_urban) {
      fprintf(testfile,"rname %s vc_prime %lf sc_prime %lf cs_prime %lf vp_prime %lf sp_prime %lf ps_prime %lf\n",rname[i],rtdata[i].flowvc_prime,rtdata[i].flowsc_prime,rtdata[i].flowcs_prime,rtdata[i].flowvp_prime,rtdata[i].flowsp_prime,rtdata[i].flowps_prime);
      }else{
	fprintf(testfile,"rname %s vc_prime %lf sc_prime %lf cs_prime %lf vp_prime %lf sp_prime %lf ps_prime %lf us_prime %lf su_prime %lf vu_prime %lf\n",rname[i],rtdata[i].flowvc_prime,rtdata[i].flowsc_prime,rtdata[i].flowcs_prime,rtdata[i].flowvp_prime,rtdata[i].flowsp_prime,rtdata[i].flowps_prime,rtdata[i].flowus_prime,rtdata[i].flowsu_prime,rtdata[i].flowvu_prime);
      }
    }

  fclose(testfile);
  
 

  return;

}
	    
	

/********************************************************************/
/** newdata= yearly interpolated data where it=0 is 1700
    transitions are computed here where it=1 is transition from 1700 to 1701
**/

void transitions(int curryear){
  
  char newpath6[190], newpath7[190], fn[30], fout[190];
  int iz, i, ic, k, m, it, ih, zmax=0, izz;
  double tester1, tester2, tester3, tester4;
  

  printf("\ncomputing transitions...\n");



  /* wh checker files written to below */

  it=0;  
      for (k=0;k<NY;k++){
	for (m=0;m<NX;m++){

 
	    if(dstatic[k][m].rcode > 0){


	      i=dstatic[k][m].newrcode;
	   
	      
	      if (BEST_CASE){
	      
	      if((dstatic[k][m].continent_code == 3) || (dstatic[k][m].continent_code == 4)){
		
		/* minimum flows, secondary priority, Eurasia only */ 
		
		alternative_smart_flow2(k,m,it);
		
		if (BEST_CASE_MIN_FLOWS_T4) {    /* abandonment in the tropical Eurasia, secondary priority */
		if((lat[k] <= 23.5) && (lat[k] >= -23.5)) adjust_smart_flow2(k,m,it,i);
		}
		if (BEST_CASE_MIN_FLOWS_T5){   /*abandonment in shiftcult map */
		if(dstatic[k][m].shiftcult==1) adjust_smart_flow2(k,m,it,i);
		}
	      }
	      else{
		
		/* minimum flows everywhere with abandonment in the tropics, primary priority */
		
		  alternative_smart_flow1(k,m,it);
		  
		  if (BEST_CASE_MIN_FLOWS_T4){    /* abandonment in tropical non Eurasia, primary priority */
		  if((lat[k] <= 23.5) && (lat[k] >= -23.5)) adjust_smart_flow1(k,m,it,i);
		  }
		  if (BEST_CASE_MIN_FLOWS_T5){   /*abandonment in shiftcult map */
		if(dstatic[k][m].shiftcult==1) adjust_smart_flow1(k,m,it,i);
		  }		  
	      }

	      
	      }else{
 
		if(rdata[i].smart_flow_option == 1){
		  
		  
		  if(rdata[i].adjust_smart_flow_option == 1) {
		    
		    alternative_smart_flow1(k,m,it);
		    
		  }
		  else if(rdata[i].adjust_smart_flow_option == 2) {
		    
		    alternative_smart_flow1(k,m,it);	     
		    adjust_smart_flow1(k,m,it,i);
		    
		    
		  }
		  else if(rdata[i].adjust_smart_flow_option == 3 && use_urban){
		    
		    alternative_smart_flow1(k,m,it);
		    
		    if((lon[m]>-106.0) && (lon[m]<-71.0) && (lat[k]<20.0) && (lat[k]>5.0)) adjust_smart_flow3(k,m,it,i);
		    if((lon[m]>-81.0) && (lon[m]<-60.0) && (lat[k]<5.0) && (lat[k]>-25.0)) adjust_smart_flow3(k,m,it,i);
		    if((lon[m]>-60.0) && (lon[m]<-40.0) && (lat[k]<5.0) && (lat[k]>-10.0)) adjust_smart_flow3(k,m,it,i);
		    if((lon[m]>-20.0) && (lon[m]<55.0) && (lat[k]<15.0) && (lat[k]>-15.0)) adjust_smart_flow3(k,m,it,i);
		    if((lon[m]>90.0) && (lon[m]<105.0) && (lat[k]<30.0) && (lat[k]>-11.0)) adjust_smart_flow3(k,m,it,i);
		    if((lon[m]>105.0) && (lon[m]<160.0) && (lat[k]<20.0) && (lat[k]>-11.0)) adjust_smart_flow3(k,m,it,i);
		  
		  }
		 
		  else if(rdata[i].adjust_smart_flow_option == 4){  

		    alternative_smart_flow1(k,m,it);

		    /* additional abandonment in tropics */
		    if((lat[k] <= 23.5) && (lat[k] >= -23.5)) adjust_smart_flow1(k,m,it,i);

		  }
		else if (rdata[i].adjust_smart_flow_option == 5){
		  alternative_smart_flow1(k,m,it);
		  if(dstatic[k][m].shiftcult==1) adjust_smart_flow1(k,m,it,i);

		}
		


		} /*end smart_flow if */

		else if(rdata[i].smart_flow_option == 2){


		  if(rdata[i].adjust_smart_flow_option == 1){

		    alternative_smart_flow2(k,m,it);

		  }
		  else if(rdata[i].adjust_smart_flow_option == 2) {	 

		    alternative_smart_flow2(k,m,it);
		    adjust_smart_flow2(k,m,it,i);

		  }
		  
		  else if(rdata[i].adjust_smart_flow_option == 3 && use_urban) {	

		    alternative_smart_flow2(k,m,it);
    
		    if((lon[m]>-106.0) && (lon[m]<-71.0) && (lat[k]<20.0) && (lat[k]>5.0)) adjust_smart_flow3(k,m,it,i);
		    if((lon[m]>-81.0) && (lon[m]<-60.0) && (lat[k]<5.0) && (lat[k]>-25.0)) adjust_smart_flow3(k,m,it,i);
		    if((lon[m]>-60.0) && (lon[m]<-40.0) && (lat[k]<5.0) && (lat[k]>-10.0)) adjust_smart_flow3(k,m,it,i);
		    if((lon[m]>-20.0) && (lon[m]<55.0) && (lat[k]<15.0) && (lat[k]>-15.0)) adjust_smart_flow3(k,m,it,i);
		    if((lon[m]>90.0) && (lon[m]<105.0) && (lat[k]<30.0) && (lat[k]>-11.0)) adjust_smart_flow3(k,m,it,i);
		    if((lon[m]>105.0) && (lon[m]<160.0) && (lat[k]<20.0) && (lat[k]>-11.0)) adjust_smart_flow3(k,m,it,i);

		  }
		  else {   /* rdata[i].adjust_smart_flow_option == 4 */

		    alternative_smart_flow2(k,m,it);

		    /* additional abandonment in tropics */
		    if((lat[k] <= 23.5) && (lat[k] >= -23.5)) adjust_smart_flow2(k,m,it,i);
		  }


		} /*end smart_flow else if */

		else{ /* rdata[i].smart_flow_option == 3 */

		  smart_flow(k,m,it);
		  /* need an adjusted smart flow for abandonment */
				
		} /*end smart_flow else if */


	      } /* end of BEST_CASE if */


	
	    } /* end of rcode if */

		
	    if(dstatic[k][m].gcode > 0){


	      i=dstatic[k][m].newgcode;
              if(dstatic[k][m].fnf == 1){
		
		  ctdata[i].converted_forest_land+=(((newdata[it].flowvp[k][m]+newdata[it].flowvc[k][m])*dstatic[k][m].vba)+((newdata[it].flowsp[k][m]+newdata[it].flowsc[k][m])*newdata[it].smb[k][m]))*garea[k][m]/1000.0; 	
		
		}
	    } /* end gcode if */
            
            if(dstatic[k][m].fnf == 1){
              newdata[it].converted_forest_land[k][m]=(((newdata[it].flowvp[k][m]+newdata[it].flowvc[k][m])*dstatic[k][m].vba)+((newdata[it].flowsp[k][m]+newdata[it].flowsc[k][m])*newdata[it].smb[k][m]))*garea[k][m]/1000.0;
            }
	
	  }  /* end of m */
	}  /* end of k */
	
	
	

      if(gridded_wh==0){
	if (TOTAL_HARVEST_SWITCH){

    printf("beginning harvest, time= %d\n",curryear); 

    if(suit_option == 1) zdis_calc(it);

    update_vb2(it);



    for (i=0;i<NCCODE;i++){
      predict_available_biomass(it,i);
      
    }

    /* divide regional WH between countries according to country ratios */
    for (i=0;i<NCCODE;i++){
      ctdata[i].wh = wh_ratio[i]*rtdata[cdata[i].rcode-1].wh;
     
    }


/* regional loop begin */

#if REGION_TEST
   for(i=REGION_TEST_INDEX;i<REGION_TEST_INDEX+1;i++){ 
#else
      for (i=0;i<NCCODE;i++){
#endif
	
	ctdata[i].stbh=ZEROVALUE;
	ctdata[i].whr=ctdata[i].wh;     /* units = MgC */ 


	/* priority for wood clearing */

	if (BEST_CASE){


	if((cdata[i].continent_code == 3) || (cdata[i].continent_code == 4)){

	  /* if Eurasia, for h1 do not count clearing in harvest, but track its total */
          /* if Eurasia, for h3 count clearing in harvest and track its total */


	  if(nodata_option == 1){

	    ctdata[i].converted_forest_land_total+=ctdata[i].converted_forest_land;
	    ctdata[i].converted_forest_land=ZEROVALUE;

	  }
	  if(nodata_option == 3){

	    if(ctdata[i].converted_forest_land <= ctdata[i].wh){
	      ctdata[i].converted_forest_land_total+=ctdata[i].converted_forest_land;
	      ctdata[i].whr-=ctdata[i].converted_forest_land;
	    }
	    else{
	      ctdata[i].whr=ZEROVALUE;
	      ctdata[i].converted_forest_land_total+=ctdata[i].converted_forest_land;
	      ctdata[i].converted_forest_land=ctdata[i].wh;
	    }

	  }
          if(nodata_option == 4){

	    ctdata[i].converted_forest_land_total+=ctdata[i].converted_forest_land;
	    ctdata[i].converted_forest_land=ZEROVALUE;
	  }
	  if((nodata_option == 5)||(nodata_option == 6)){

	    ctdata[i].converted_forest_land_total+=ctdata[i].converted_forest_land;
	    ctdata[i].converted_forest_land=ZEROVALUE;
	  }

	}
	
        else{

	  /* if not Eurasia, do not count clearing in harvest, but track its total */

	  ctdata[i].converted_forest_land_total+=ctdata[i].converted_forest_land;
	  ctdata[i].converted_forest_land=ZEROVALUE;

	}


	}else{


	/* checks for counting converted land in wood harvest */

	if(cdata[i].converted_forest_land_option == 2){
	  if(ctdata[i].converted_forest_land <= ctdata[i].wh){
            ctdata[i].converted_forest_land_total+=ctdata[i].converted_forest_land;
	    ctdata[i].whr-=ctdata[i].converted_forest_land;
	  }
	  else{
	    ctdata[i].whr=ZEROVALUE;
	    ctdata[i].converted_forest_land_total+=ctdata[i].converted_forest_land;
    	    ctdata[i].converted_forest_land=ctdata[i].wh;
	  }
	}
	else{
	  ctdata[i].converted_forest_land_total+=ctdata[i].converted_forest_land;
	  ctdata[i].converted_forest_land=ZEROVALUE;
	}

	}


	/* priority for wood harvest */

      } /* end country loop */

      /* predict unmet WH in each country and track regional unmet/available WH predictions */
      for (i=0;i<NCCODE;i++){
	if(ctdata[i].whr <= ctdata[i].predict_b){
          ctdata[i].predict_b -= ctdata[i].whr;
          rtdata[cdata[i].rcode-1].reg_avail+=ctdata[i].predict_b;
	}
	else{
	  rtdata[cdata[i].rcode-1].reg_unmet+=(ctdata[i].whr-ctdata[i].predict_b);
	  ctdata[i].whr = ctdata[i].predict_b;
	  ctdata[i].predict_b = 0;
	}
      }

      /* reassign regional unmet to countries with available land */
      for (i=0;i<NCCODE;i++){
	if(rtdata[cdata[i].rcode-1].reg_unmet<rtdata[cdata[i].rcode-1].reg_avail){
	  ctdata[i].whr+=ctdata[i].predict_b/rtdata[cdata[i].rcode-1].reg_avail*rtdata[cdata[i].rcode-1].reg_unmet;
	}
        else{
          if (ctdata[i].predict_b != 0.){
            printf("tcxglm1: %f %f %f\n",ctdata[i].whr,ctdata[i].predict_b,rtdata[cdata[i].rcode-1].reg_avail);
  	    ctdata[i].whr+=ctdata[i].predict_b/rtdata[cdata[i].rcode-1].reg_avail*rtdata[cdata[i].rcode-1].reg_avail;
	  }
        }
      }
   } /* end gridded_wh if  */
    
   if(gridded_wh==1){
     harvest_gridded_data(it);
     for (i=0;i<NCCODE;i++){
       country_unmet_wh(it,i);
      
     }
    if(suit_option == 1) zdis_calc(it);

    update_vb3(it);

   } /* end gridded_wh if */

#if REGION_TEST
   for(i=REGION_TEST_INDEX;i<REGION_TEST_INDEX+1;i++){ 
#else
      for (i=0;i<NCCODE;i++){
#endif

	if (BEST_CASE){


	/* if best case, Eurasia does secondary priority for wood, non-Eurasia does primay */

	if((cdata[i].continent_code == 3) || (cdata[i].continent_code == 4)){	



	  if (SECONDARY_HARVEST_SWITCH ){
	  if(gridded_wh==0){
	    secondary_harvest2(it,i);
	  }
	  else{
	    secondary_harvest3(it,i);
	  }

	  
	  } 	  	  

	  ctdata[i].vwh = ctdata[i].whr;

	  if (VIRGIN_HARVEST_SWITCH ){
	  if(gridded_wh==0){
	    virgin_harvest2(it,i,curryear);
	  }
	  else{
	    virgin_harvest3(it,i,curryear);
	  }

	  }

	}
	else{	 


	  ctdata[i].vwh = ctdata[i].whr;
	  
	  if (VIRGIN_HARVEST_SWITCH){
	  if(gridded_wh==0){
	    virgin_harvest2(it,i,curryear);
	  }
	  else{
	    virgin_harvest3(it,i,curryear);
	  }

	 
	  }
	 
 
	  if (SECONDARY_HARVEST_SWITCH ){
	  if(gridded_wh==0){
	    secondary_harvest2(it,i);
	  }
	  else{
	    secondary_harvest3(it,i);
	  }

	  
	  }

	}


	}else{



	/********************************/
	if(cdata[i].harvest_option == 2){

	  if (SECONDARY_HARVEST_SWITCH ){
	  if(gridded_wh==0){
	    secondary_harvest2(it,i);
	  }
	  else{
	    secondary_harvest3(it,i);
	  }

	  
	  }  
	  
	  ctdata[i].vwh = ctdata[i].whr;
	  
	  if (VIRGIN_HARVEST_SWITCH ){
	  if(gridded_wh==0){
	    virgin_harvest2(it,i,curryear);
	  }
	  else{
	    virgin_harvest3(it,i,curryear);
	  }

	  
	}
	  
	}
	else{

	  ctdata[i].vwh = ctdata[i].whr;

	  if (VIRGIN_HARVEST_SWITCH ){
	 	 
	  if(gridded_wh==0){
	    virgin_harvest2(it,i,curryear);
	  }
	  else{
	    virgin_harvest3(it,i,curryear);
	  }
	  
	  }
	  


	  if (SECONDARY_HARVEST_SWITCH ){
	   if(gridded_wh==0){
	    secondary_harvest2(it,i);
	  }
	  else{
	    secondary_harvest3(it,i);
	  }

	  
	  } 
	 
	}
	  /********************************/

	}  /* end if for BEST_CASE if */	


	if (FORCE_HARVEST_SWITCH ){
       

	  if(ctdata[i].whr > 0.01){
	    if(gridded_wh==0){
	      force_harvest2(it,i,curryear);
	    }
	    else{
	      force_harvest3(it,i,curryear);
	    }
	  }
	
	}



    
      } /* end regional */
      
	
    
   }   /* end total_harvest_switch */


      update_states(it,zmax); 


      
      //jt   }  /* end time */
    
    
       
 
   return;
    
}

/***********************************************************************/
void secondary_harvest(int it, int i){

  int k, m;
  double sbh, sbhtest=ZEROVALUE; 



  /* calculate secondary biomass harvest for each grid cell (sbh), and
     track the accumulating amount of sbh at the country level (stbh)  */



  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){


#if COUNTRY_TEST

      if(dstatic[k][m].rcode==REGION_TEST_GCODE){

#else
	if(dstatic[k][m].rcode > 0){
#endif


	  if(dstatic[k][m].rcode == rdata[i].rcode){



	    if(dstatic[k][m].fnf == 1){

	      if(rtdata[i].whr > ZEROVALUE){

		sbh=newdata[it].smb[k][m]*
		  prob_harv(newdata[it].smb[k][m])*
		  newdata[it].s[k][m]*garea[k][m];


		if(sbh <= rtdata[i].whr*1000.){
		  newdata[it].sbh[k][m]=sbh;
		}
		else{
		  newdata[it].sbh[k][m]=rtdata[i].whr*1000.;
		}
		
 if(dstatic[k][m].rcode == 3) sbhtest+=newdata[it].sbh[k][m]/1000.;
 
	      rtdata[i].stbh+=newdata[it].sbh[k][m]/1000.;
	      
	      
	      rtdata[i].whr-=newdata[it].sbh[k][m]/1000.;

		
	      } /*end of whr */

	      
	      
	    } /* end of fnf */
	  }  /* end of gcode */
	}
	
	
      } /* end of m */
    } /* end of k */

    
    

    return;

  }

/***********************************************************************/
  void secondary_harvest2(int it, int i){

  int k, m;
  double sbh, sbhtest=ZEROVALUE; 



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

	      if(ctdata[i].whr > ZEROVALUE){

		sbh=newdata[it].smb[k][m]*
		  prob_harv(newdata[it].smb[k][m])*
		  newdata[it].s[k][m]*garea[k][m];


		if(sbh <= ctdata[i].whr*1000.){
		  newdata[it].sbh[k][m]=sbh;
		}
		else{
		  newdata[it].sbh[k][m]=ctdata[i].whr*1000.;
		}
		
 if(dstatic[k][m].rcode == 3) sbhtest+=newdata[it].sbh[k][m]/1000.;
 
	      ctdata[i].stbh+=newdata[it].sbh[k][m]/1000.;
	      
	      
	      ctdata[i].whr-=newdata[it].sbh[k][m]/1000.;

		
	      } /*end of whr */

	      
	    } /* end of fnf */
	  }  /* end of gcode */
	}
	
	
      } /* end of m */
    } /* end of k */

    
   

    return;

  }

/***********************************************************************/
void virgin_harvest(int it, int i,int curryear){
                      
  FILE *testfile;
  int k, m, iz, izz, zmax=0, im, j;
  double total_avail, vbhsummer=ZEROVALUE;

  testfile=fopen("zmax.test",outstat);

  
  /* first determine zmax; the maximum # cells away from the focal cell (the agricultural
     cell) needed to go to attain wood harvest demand.
     zmax=0; have enough biomass in focal (agricultural) cell 
     zmax=1; need to go to the next adjacent cell 
     zmax=2; etc...zmax=10 or MAXZ (MAXZ=11; maximum z before we get tired, and then 
     we spread remaining harvest over all forested cells with iz >= this value */
  
 
  if(rtdata[i].vwh < rztdata[i][0].vb)
    zmax=0;

  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb){
    zmax=1;
    
  }
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb){
    zmax=2;
    
  }
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +rztdata[i][3].vb){
    zmax=3;
  }
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +rztdata[i][3].vb+rztdata[i][4].vb){
    zmax=4;
  }
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +rztdata[i][3].vb+rztdata[i][4].vb+rztdata[i][5].vb){
    zmax=5;
  }
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +rztdata[i][3].vb+rztdata[i][4].vb+rztdata[i][5].vb+
	  rztdata[i][6].vb){
    zmax=6;
  }
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +rztdata[i][3].vb+rztdata[i][4].vb+rztdata[i][5].vb+
	  rztdata[i][6].vb+rztdata[i][7].vb){
    zmax=7;
  }
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +                +rztdata[i][3].vb+rztdata[i][4].vb+rztdata[i][5].vb+
	  rztdata[i][6].vb+rztdata[i][7].vb+rztdata[i][8].vb){
    zmax=8;
  }
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +rztdata[i][3].vb+rztdata[i][4].vb+rztdata[i][5].vb+
	  rztdata[i][6].vb+rztdata[i][7].vb+rztdata[i][8].vb
	  +rztdata[i][9].vb){
    zmax=9;
  }
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +rztdata[i][3].vb+rztdata[i][4].vb+rztdata[i][5].vb+
	  rztdata[i][6].vb+rztdata[i][7].vb+rztdata[i][8].vb
	  +rztdata[i][9].vb+rztdata[i][10].vb){
    zmax=10;
  }

#if MAXZ_21	        
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +rztdata[i][3].vb+rztdata[i][4].vb+rztdata[i][5].vb+
	  rztdata[i][6].vb+rztdata[i][7].vb+rztdata[i][8].vb
	  +rztdata[i][9].vb+rztdata[i][10].vb+rztdata[i][11].vb){
    zmax=11;
  }
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +rztdata[i][3].vb+rztdata[i][4].vb+rztdata[i][5].vb+
	  rztdata[i][6].vb+rztdata[i][7].vb+rztdata[i][8].vb
	  +rztdata[i][9].vb+rztdata[i][10].vb+rztdata[i][11].vb+rztdata[i][12].vb){
    zmax=12;
  }	 
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +rztdata[i][3].vb+rztdata[i][4].vb+rztdata[i][5].vb+
	  rztdata[i][6].vb+rztdata[i][7].vb+rztdata[i][8].vb
	  +rztdata[i][9].vb+rztdata[i][10].vb+rztdata[i][11].vb+rztdata[i][12].vb+rztdata[i][13].vb){
    zmax=13;
  }
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +rztdata[i][3].vb+rztdata[i][4].vb+rztdata[i][5].vb+
	  rztdata[i][6].vb+rztdata[i][7].vb+rztdata[i][8].vb
	  +rztdata[i][9].vb+rztdata[i][10].vb+rztdata[i][11].vb+rztdata[i][12].vb+rztdata[i][13].vb+rztdata[i][14].vb){
    zmax=14;
  }
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +rztdata[i][3].vb+rztdata[i][4].vb+rztdata[i][5].vb+
	  rztdata[i][6].vb+rztdata[i][7].vb+rztdata[i][8].vb
	  +rztdata[i][9].vb+rztdata[i][10].vb+rztdata[i][11].vb+rztdata[i][12].vb+rztdata[i][13].vb+rztdata[i][14].vb+rztdata[i][15].vb){
    zmax=15;
  }
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +rztdata[i][3].vb+rztdata[i][4].vb+rztdata[i][5].vb+
	  rztdata[i][6].vb+rztdata[i][7].vb+rztdata[i][8].vb
	  +rztdata[i][9].vb+rztdata[i][10].vb+rztdata[i][11].vb+rztdata[i][12].vb+rztdata[i][13].vb+rztdata[i][14].vb+rztdata[i][15].vb+rztdata[i][16].vb){
    zmax=16;
  }
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +rztdata[i][3].vb+rztdata[i][4].vb+rztdata[i][5].vb+
	  rztdata[i][6].vb+rztdata[i][7].vb+rztdata[i][8].vb
	  +rztdata[i][9].vb+rztdata[i][10].vb+rztdata[i][11].vb+rztdata[i][12].vb+rztdata[i][13].vb+rztdata[i][14].vb+rztdata[i][15].vb+rztdata[i][16].vb+rztdata[i][17].vb){
    zmax=17;
  }	
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +rztdata[i][3].vb+rztdata[i][4].vb+rztdata[i][5].vb+
	  rztdata[i][6].vb+rztdata[i][7].vb+rztdata[i][8].vb
	  +rztdata[i][9].vb+rztdata[i][10].vb+rztdata[i][11].vb+rztdata[i][12].vb+rztdata[i][13].vb+rztdata[i][14].vb+rztdata[i][15].vb+rztdata[i][16].vb+rztdata[i][17].vb+rztdata[i][18].vb){
    zmax=18;
  }
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +rztdata[i][3].vb+rztdata[i][4].vb+rztdata[i][5].vb+
	  rztdata[i][6].vb+rztdata[i][7].vb+rztdata[i][8].vb
	  +rztdata[i][9].vb+rztdata[i][10].vb+rztdata[i][11].vb+rztdata[i][12].vb+rztdata[i][13].vb+rztdata[i][14].vb+rztdata[i][15].vb+rztdata[i][16].vb+rztdata[i][17].vb+rztdata[i][18].vb+rztdata[i][19].vb){
    zmax=19;
  }
  else if(rtdata[i].vwh < rztdata[i][0].vb+rztdata[i][1].vb
	  +rztdata[i][2].vb
	  +rztdata[i][3].vb+rztdata[i][4].vb+rztdata[i][5].vb+
	  rztdata[i][6].vb+rztdata[i][7].vb+rztdata[i][8].vb
	  +rztdata[i][9].vb+rztdata[i][10].vb+rztdata[i][11].vb+rztdata[i][12].vb+rztdata[i][13].vb+rztdata[i][14].vb+rztdata[i][15].vb+rztdata[i][16].vb+rztdata[i][17].vb+rztdata[i][18].vb+rztdata[i][19].vb+rztdata[i][20].vb){
    zmax=20;
  }	 
#endif

  else {
    zmax=MAXZ; 
  }
  
  


  if(zmax == 0){
    rtdata[i].wh_at_zmax=rtdata[i].whr;
  }
  else{
    for(j=0;j<zmax;j++) rtdata[i].wh_lessthan_zmax+=rztdata[i][j].vb;
    rtdata[i].wh_at_zmax=rtdata[i].whr-rtdata[i].wh_lessthan_zmax;
  }

  

  
  /***************/
  /* determine virgin biomass harvest (vbh) in units of kgC and the flow
     of virgin to secondary (flowvs) */
  
 
  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){


#if REGION_TEST

      if(dstatic[k][m].rcode==REGION_TEST_GCODE){
#else
	if(dstatic[k][m].rcode > 0){
#endif


	  if(dstatic[k][m].rcode == rdata[i].rcode) {
	 
	    if(rtdata[i].whr > ZEROVALUE) {


	      iz=newdata[it].zdis[k][m];





	      if(dstatic[k][m].fnf == 1){



		/* determine how to take biomass for harvest demand */

		if(zmax<MAXZ){ /*business as usual*/


		  if(iz < zmax){ /* take as much as possible to fulfill harvest demand */
		      
		    newdata[it].flowvs[k][m]=newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m];
		    
		    newdata[it].vbh[k][m]=newdata[it].flowvs[k][m]*garea[k][m]*dstatic[k][m].vba;
		       
		    
		  }
		  else if(iz == zmax){


		    if(rztdata[i][iz].vb > ZEROVALUE){
		      
		      newdata[it].flowvs[k][m]=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*rtdata[i].wh_at_zmax/rztdata[i][iz].vb;
		      
		      
		      newdata[it].vbh[k][m]=newdata[it].flowvs[k][m]*garea[k][m]*dstatic[k][m].vba;
		    }	  
	

		  }

		  if(it ==0 && dstatic[k][m].rcode == 3) {

		    vbhsummer+=newdata[it].vbh[k][m]/1000;

		  
		  }



		} /* end of business as usual, zmax < MAXZ */

		else{ 
		  if(iz < MAXZ){ /*zmax=MAXZ*/
		    /*most dry countries with reported wood harvest and not enough forest*/
		    /* take everything possible */
		    
		    newdata[it].flowvs[k][m]=(newdata[it].v[k][m]-
					      newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m]);
		    
		    
		    newdata[it].vbh[k][m]=newdata[it].flowvs[k][m]*garea[k][m]*
		      dstatic[k][m].vba;
		    
		    
		  }
		  else if(iz<=NZ){ /*from MAXZ=zmax to NZ*/
		    
		    /* rare, reported wood harvest, but have to go more than MAXZ to get it*/
		    /* take what is needed proportionally out of each cell */

		    total_avail=ZEROVALUE;
		    
		    for(izz=MAXZ;izz<NZ;izz++)
		      total_avail+=rztdata[i][izz].vb;
		    
		    if(total_avail > ZEROVALUE){
		      
		      if(rtdata[i].wh_at_zmax <= total_avail){
			
			newdata[it].flowvs[k][m]=(newdata[it].v[k][m]-
						  newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*
			  rtdata[i].wh_at_zmax/total_avail;
		      }
		      else{
			newdata[it].flowvs[k][m]=(newdata[it].v[k][m]-
						  newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m]);
			
		      }
		    }
		    
		    /*safe, but still may not get enough harvest*/
		    
		   
		    newdata[it].vbh[k][m]=newdata[it].flowvs[k][m]*garea[k][m]*
		      dstatic[k][m].vba;
		    
		    
		    fprintf(testfile,"potential zmax trouble spot, %d %s zmax %d iz %d flowvs %lf wh_at_zmax %lf total %lf\n",curryear,cname[i],zmax,iz,newdata[it].flowvs[k][m],rtdata[i].wh_at_zmax,total_avail);  
		    

		  }
		  else {
		    printf("iz %d zmax %d\n",iz,zmax);
		  }
		  
		}  /* end of else iz < MAXZ */



	      }/* end of fnf*/

rtdata[i].whr-=newdata[it].vbh[k][m]/1000.;

	    }/* end of whr*/


	   

	  }/*gcode*/
	}


      }  /* end of m */
    } /* end of k */


    fclose(testfile);

 
    return;

  }
/***********************************************************************/
  void virgin_harvest2(int it, int i,int curryear){
                      
  FILE *testfile;
  int k, m, iz, izz, zmax=0, im, j;
  double total_avail, vbhsummer=ZEROVALUE;
  
  /* first determine zmax; the maximum # cells away from the focal cell (the agricultural
     cell) needed to go to attain wood harvest demand.
     zmax=0; have enough biomass in focal (agricultural) cell 
     zmax=1; need to go to the next adjacent cell 
     zmax=2; etc...zmax=10 or MAXZ (MAXZ=11; maximum z before we get tired, and then 
     we spread remaining harvest over all forested cells with iz >= this value */
  
  testfile=fopen("zmax.test",outstat);
 

  if(ctdata[i].vwh < cztdata[i][0].vb)
    zmax=0;
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb){
    zmax=1;
    
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb){
    zmax=2;
   
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb){
    zmax=3;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb){
    zmax=4;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb){
    zmax=5;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb){
    zmax=6;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb){
    zmax=7;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +                +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb){
    zmax=8;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb){
    zmax=9;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb){
    zmax=10;
  }

#if MAXZ_21	        
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb){
    zmax=11;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb){
    zmax=12;
  }	 
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb+cztdata[i][13].vb){
    zmax=13;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb+cztdata[i][13].vb+cztdata[i][14].vb){
    zmax=14;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb+cztdata[i][13].vb+cztdata[i][14].vb+cztdata[i][15].vb){
    zmax=15;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb+cztdata[i][13].vb+cztdata[i][14].vb+cztdata[i][15].vb+cztdata[i][16].vb){
    zmax=16;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb+cztdata[i][13].vb+cztdata[i][14].vb+cztdata[i][15].vb+cztdata[i][16].vb+cztdata[i][17].vb){
    zmax=17;
  }	
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb+cztdata[i][13].vb+cztdata[i][14].vb+cztdata[i][15].vb+cztdata[i][16].vb+cztdata[i][17].vb+cztdata[i][18].vb){
    zmax=18;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb+cztdata[i][13].vb+cztdata[i][14].vb+cztdata[i][15].vb+cztdata[i][16].vb+cztdata[i][17].vb+cztdata[i][18].vb+cztdata[i][19].vb){
    zmax=19;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb+cztdata[i][13].vb+cztdata[i][14].vb+cztdata[i][15].vb+cztdata[i][16].vb+cztdata[i][17].vb+cztdata[i][18].vb+cztdata[i][19].vb+cztdata[i][20].vb){
    zmax=20;
  }	 
#endif

  else {
    zmax=MAXZ; 
  }

  if(zmax == 0){
    ctdata[i].wh_at_zmax=ctdata[i].whr;
  }
  else{
    for(j=0;j<zmax;j++) ctdata[i].wh_lessthan_zmax+=cztdata[i][j].vb;
    ctdata[i].wh_at_zmax=ctdata[i].whr-ctdata[i].wh_lessthan_zmax;
  }

  

  
  /***************/
  /* determine virgin biomass harvest (vbh) in units of kgC and the flow
 [k][m]    of virgin to secondary (flowvs) */
  
 
  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){


#if COUNTRY_TEST

      if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){
#else
	if(dstatic[k][m].gcode > 0){
#endif


	  if(dstatic[k][m].gcode == cdata[i].ccode) {
	 
	    if(ctdata[i].whr > ZEROVALUE) {


	      iz=newdata[it].zdis[k][m];





	      if(dstatic[k][m].fnf == 1){



		/* determine how to take biomass for harvest demand */

		if(zmax<MAXZ){ /*business as usual*/


		  if(iz < zmax){ /* take as much as possible to fulfill harvest demand */
		      
		    newdata[it].flowvs[k][m]=newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m];
		    
		    newdata[it].vbh[k][m]=newdata[it].flowvs[k][m]*garea[k][m]*dstatic[k][m].vba;
		       
		    
		  }
		  else if(iz == zmax){


		    if(cztdata[i][iz].vb > ZEROVALUE){
		      
		      newdata[it].flowvs[k][m]=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*ctdata[i].wh_at_zmax/cztdata[i][iz].vb;
		      
		      
		      newdata[it].vbh[k][m]=newdata[it].flowvs[k][m]*garea[k][m]*dstatic[k][m].vba;
		    }	  
	

		  }

		  if(it ==0 && dstatic[k][m].rcode == 3) {

		    vbhsummer+=newdata[it].vbh[k][m]/1000;

		   

		  }



		} /* end of business as usual, zmax < MAXZ */

		else{ 

		  if(iz < MAXZ){ /*zmax=MAXZ*/
		    /*most dry countries with reported wood harvest and not enough forest*/
		    /* take everything possible */
		    
		    newdata[it].flowvs[k][m]=(newdata[it].v[k][m]-
					      newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m]);
		    
		    
		    newdata[it].vbh[k][m]=newdata[it].flowvs[k][m]*garea[k][m]*
		      dstatic[k][m].vba;
		    
		    
		  }
		  else if(iz<=NZ){ /*from MAXZ=zmax to NZ*/
		    
		    /* rare, reported wood harvest, but have to go more than MAXZ to get it*/
		    /* take what is needed proportionally out of each cell */

		    total_avail=ZEROVALUE;
		    
		    for(izz=MAXZ;izz<NZ;izz++)
		      total_avail+=cztdata[i][izz].vb;
		    
		    if(total_avail > ZEROVALUE){
		      
		      if(ctdata[i].wh_at_zmax <= total_avail){
			
			newdata[it].flowvs[k][m]=(newdata[it].v[k][m]-
						  newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*
			  ctdata[i].wh_at_zmax/total_avail;
		      }
		      else{
			newdata[it].flowvs[k][m]=(newdata[it].v[k][m]-
						  newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m]);
			
		      }
		    }
		    
		    /*safe, but still may not get enough harvest*/
		    
		   
		    
		    newdata[it].vbh[k][m]=newdata[it].flowvs[k][m]*garea[k][m]*
		      dstatic[k][m].vba;
		    
		    
		    fprintf(testfile,"potential zmax trouble spot, %d %s zmax %d iz %d flowvs %lf wh_at_zmax %lf total %lf\n",curryear,cname[i],zmax,iz,newdata[it].flowvs[k][m],ctdata[i].wh_at_zmax,total_avail);  
		    

		  }
		  else {
		    printf("iz %d zmax %d\n",iz,zmax);
		  }
		  
		}  /* end of else iz < MAXZ */



	      }/* end of fnf*/

ctdata[i].whr-=newdata[it].vbh[k][m]/1000.;

	    }/* end of whr*/


	  

	  }/*gcode*/
	}


      }  /* end of m */
    } /* end of k */

    fclose(testfile);


    return;

  }
/***********************************************************************/
  void force_harvest(int it, int i, int curryear){

  FILE *testfile;
  int k, m, flag=0;
  double whr_orig=ZEROVALUE, flowvs_add=ZEROVALUE, flowss_add=ZEROVALUE, wh_from_sbh_only=ZEROVALUE;
  double tester1=ZEROVALUE, tester2=ZEROVALUE, tester3=ZEROVALUE, tester4=ZEROVALUE, tester5=ZEROVALUE;
 
 
  /* hold original amount of whr for spreading proportionally */

  whr_orig=rtdata[i].whr;



  /* compute country level smb, smb_nf, vnfb, s area for force harvest */
  
  
  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){
      
      
#if REGION_TEST
      
      if(dstatic[k][m].rcode==REGION_TEST_GCODE){
	   
#else
	if(dstatic[k][m].rcode > 0){
#endif
	     
	  if(dstatic[k][m].rcode == rdata[i].rcode){	   
	 	     
	    if(dstatic[k][m].fnf == 0){
	     
	      /* country level virgin nonforested biomass */
	      rtdata[i].vnfb+=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*garea[k][m]*dstatic[k][m].vba/1000.;
		 
	      /* country level secondary nonforested biomass */
	      rtdata[i].smb_nf+=(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*garea[k][m]*newdata[it].smb[k][m]/1000.;

	      /* country level secondary nonforested area */
	      rtdata[i].sarea_nf+=(newdata[it].s[k][m]+newdata[it].flowps[k][m]+newdata[it].flowcs[k][m]+newdata[it].flowvs[k][m]+newdata[it].flowus[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*garea[k][m];
 
	    }
	    else {

	      /* country level secondary forest biomass from sbh only */
	      wh_from_sbh_only+=newdata[it].sbh[k][m]/1000.;
	     
	      /* country level secondary forest biomass */
	      rtdata[i].smb+=(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*((garea[k][m]*newdata[it].smb[k][m]-newdata[it].sbh[k][m])/1000.); 

	     /* country level secondary forested area */
	      rtdata[i].sarea+=(newdata[it].s[k][m]+newdata[it].flowps[k][m]+newdata[it].flowcs[k][m]+newdata[it].flowvs[k][m]+newdata[it].flowus[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*garea[k][m];
	      
	    }
	
   
	  }  /* end of gcode */
	}
       
      } /* end of m */
    } /* end of k */
   

    /* 4 cases for getting remaining wood harvest (whr) */


    if(rtdata[i].smb >= rtdata[i].whr){


      flag=1;  /* case 1, taking needed younger secondary forest proportionally to satisfy whr */                                 

      for (k=0;k<NY;k++){
	for (m=0;m<NX;m++){
	 
#if REGION_TEST
	 
	  if(dstatic[k][m].rcode==REGION_TEST_GCODE){
	   
#else
	    if(dstatic[k][m].rcode > 0){
#endif
	     
	      if(dstatic[k][m].rcode == rdata[i].rcode){

	       
		if(dstatic[k][m].fnf == 1){
		 
		  if(rtdata[i].whr > ZEROVALUE){
		   
		   
		    newdata[it].sbh2[k][m]=whr_orig/rtdata[i].smb*(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*(garea[k][m]*newdata[it].smb[k][m]-newdata[it].sbh[k][m]);
		   
		    rtdata[i].whr-=newdata[it].sbh2[k][m]/1000.;
		   
		   
		  } /*end of whr */
		} /*end of fnf */
	      }  /* end of gcode */
	    }
	   
	  } /* end of m */
	} /* end of k */

       
      }/* end case 1*/
     

     
      else if((rdata[i].smart_flow_option == 1) && ((rtdata[i].smb+rtdata[i].vnfb) >= rtdata[i].whr)){    

	flag=2;  /* case 2, primary priority, not enough young secondary forest to satisfy demand, 
		   therefore taking all young secondary forest and needed virgin nonforest proportionally 
		   to satisfy whr */

  

	for (k=0;k<NY;k++){
	  for (m=0;m<NX;m++){
	   
#if REGION_TEST
	 
	    if(dstatic[k][m].rcode==REGION_TEST_GCODE){
	   
#else
	      if(dstatic[k][m].rcode > 0){
#endif
	     
		if(dstatic[k][m].rcode == rdata[i].rcode){

		 
		  if(rtdata[i].whr > ZEROVALUE){
		   
		    if(dstatic[k][m].fnf == 1) {

		      newdata[it].sbh2[k][m]=(newdata[it].smb[k][m]*garea[k][m]-newdata[it].sbh[k][m])*(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m]);

		      rtdata[i].whr-=newdata[it].sbh2[k][m]/1000.;
rtdata[i].fh_sbh2+=newdata[it].sbh2[k][m]/1000.;
		    } /* end of fnf*/		   
		  } /*end of whr*/
		}  /* end of gcode */
	      }
	     
	    } /* end of m */
	  } /* end of k */



	 /* take the whr proportionally out of vnfb; note needed to consider how much was taken
            out of secondary forest (rtdata[i].smb) */


	 for (k=0;k<NY;k++){
	   for (m=0;m<NX;m++){
	     
#if REGION_TEST
	     
	     if(dstatic[k][m].rcode==REGION_TEST_GCODE){
	   
#else
	       if(dstatic[k][m].rcode > 0){
#endif
	     
		 if(dstatic[k][m].rcode == rdata[i].rcode){

		 
		   if(rtdata[i].whr > ZEROVALUE){
		   
		     if(dstatic[k][m].fnf == 0) {
		     
	     	       flowvs_add=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*(whr_orig-rtdata[i].smb)/rtdata[i].vnfb; 
		     

		       if(flowvs_add <= -ZEROVALUE) {flowvs_add = ZEROVALUE;}
		     
		       newdata[it].flowvs[k][m]+=flowvs_add;
		     
		       newdata[it].vbh2[k][m]=flowvs_add*dstatic[k][m].vba*garea[k][m];
		     
		       rtdata[i].whr-=newdata[it].vbh2[k][m]/1000.;
		  rtdata[i].fh_vbh2+=newdata[it].vbh2[k][m]/1000.;   
		     } /*end of fnf */
		   
		    

		   
		   } /*end of whr*/
		 
		 }  /* end of gcode */
	       }
	     
	     } /* end of m */
	   } /* end of k */

	   
	 } /*end case 2 */
	 

	 
	 else if((rdata[i].smart_flow_option == 2) && ((rtdata[i].smb+rtdata[i].smb_nf) >= rtdata[i].whr)){  

	   flag=3;  /* case 3, secondary priority, taking all young secondary forest and needed secondary 
                   nonforest proportionally to satisfy whr */


     
	   for (k=0;k<NY;k++){
	     for (m=0;m<NX;m++){
	   
#if REGION_TEST
	 
	       if(dstatic[k][m].rcode==REGION_TEST_GCODE){
	   
#else
		 if(dstatic[k][m].rcode > 0){
#endif
	     
		   if(dstatic[k][m].rcode == rdata[i].rcode){

		 
		     if(rtdata[i].whr > ZEROVALUE){
		   
		       if(dstatic[k][m].fnf == 1) {

			 newdata[it].sbh2[k][m]=(newdata[it].smb[k][m]*garea[k][m]-newdata[it].sbh[k][m])*(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m]);

			 rtdata[i].whr-=newdata[it].sbh2[k][m]/1000.;

		       } /* end of fnf*/		   
		     } /*end of whr*/
		   }  /* end of gcode */
		 }
	     
	       } /* end of m */
	     } /* end of k */


	     /* take the whr proportionally out of secondary nonforest; note needed to consider how much was taken
		out of secondary forest (rtdata[i].smb) */


	     for (k=0;k<NY;k++){
	       for (m=0;m<NX;m++){
	     
#if REGION_TEST
	     
		 if(dstatic[k][m].rcode==REGION_TEST_GCODE){
	   
#else
		   if(dstatic[k][m].rcode > 0){
#endif
	     
		     if(dstatic[k][m].rcode == rdata[i].rcode){
		       
		       if(rtdata[i].whr > ZEROVALUE){
		   
			 if(dstatic[k][m].fnf == 0) {
		     
			    flowss_add=(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*(whr_orig-rtdata[i].smb)/rtdata[i].smb_nf;
			   
			   if(flowss_add <= -ZEROVALUE) {flowss_add = ZEROVALUE;}

			   newdata[it].sbh3[k][m]=flowss_add*newdata[it].smb[k][m]*garea[k][m];
		     
			   rtdata[i].whr-=newdata[it].sbh3[k][m]/1000.;


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

	    
#if REGION_TEST
	 
		   if(dstatic[k][m].rcode==REGION_TEST_GCODE){
	   
#else
		     if(dstatic[k][m].rcode > 0){
#endif
	     
		       if(dstatic[k][m].rcode == rdata[i].rcode){

		   
			 if(dstatic[k][m].fnf == 1){

			   /* take all young secondary forest */
		       
			   newdata[it].sbh2[k][m]=(newdata[it].smb[k][m]*garea[k][m]-newdata[it].sbh[k][m])*(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m]);

			   rtdata[i].whr-=newdata[it].sbh2[k][m]/1000.;
rtdata[i].fh_sbh2+=newdata[it].sbh2[k][m]/1000.;
			 }
		       }
		     }
		   }
		 }

	
	       for (k=0;k<NY;k++){
		 for (m=0;m<NX;m++){

	    
#if REGION_TEST
	 
		   if(dstatic[k][m].rcode==REGION_TEST_GCODE){
	   
#else
		     if(dstatic[k][m].rcode > 0){
#endif
	     
		       if(dstatic[k][m].rcode == rdata[i].rcode){

		   
			 if(dstatic[k][m].fnf == 0){

			   if(rdata[i].smart_flow_option == 1){

			     /* primary priority, first take all vnfb */


			     flowvs_add=newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m];
		       
			     if(flowvs_add <= -ZEROVALUE) {flowvs_add = ZEROVALUE;}
		       
			     newdata[it].flowvs[k][m]+=flowvs_add;
		     
			     newdata[it].vbh2[k][m]=flowvs_add*dstatic[k][m].vba*garea[k][m];
		     
			     rtdata[i].whr-=newdata[it].vbh2[k][m]/1000.;
rtdata[i].fh_vbh2+=newdata[it].vbh2[k][m]/1000.;
  


			     /* primary priority, not enough smb and vnfb, determined how much can be taken from smb_nf */

			     /* first try taking it proportionally */

 if(rtdata[i].smb_nf >= (whr_orig-rtdata[i].smb-rtdata[i].vnfb)){

     
		    
			         flowss_add=(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*(whr_orig-rtdata[i].smb-rtdata[i].vnfb)/rtdata[i].smb_nf; 
			      
			   
			       if(flowss_add <= -ZEROVALUE) {flowss_add = ZEROVALUE;}
			   
			       newdata[it].sbh3[k][m]=flowss_add*newdata[it].smb[k][m]*garea[k][m];
			   
			           rtdata[i].whr-=newdata[it].sbh3[k][m]/1000.;   
			        rtdata[i].fh_sbh3+=newdata[it].sbh3[k][m]/1000.; 
			       
			     }
			     else{ /* not enough in smb_nf, therefore take all remaining smb_nf, the rest is unmet */
			 
			       flowss_add=newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m];
			   
			       if(flowss_add <= -ZEROVALUE) {flowss_add = ZEROVALUE;}
			   
			       newdata[it].sbh3[k][m]=flowss_add*newdata[it].smb[k][m]*garea[k][m];

			        rtdata[i].whr-=newdata[it].sbh3[k][m]/1000.; 
			        rtdata[i].fh_sbh3+=newdata[it].sbh3[k][m]/1000.; 
			     }


			   } 
			   else { /* if smart_flow_option == 2 */
		       
			     /* secondary priority, first take all smb_nf */

			     flowss_add=newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m];
			   
			     if(flowss_add <= -ZEROVALUE) {flowss_add = ZEROVALUE;}
			   
			     newdata[it].sbh3[k][m]=flowss_add*newdata[it].smb[k][m]*garea[k][m];

			     rtdata[i].whr-=newdata[it].sbh3[k][m]/1000.;
rtdata[i].fh_sbh3+=newdata[it].sbh3[k][m]/1000.;


			     /* secondary priority, not enough smb and smb_nf, determined how much can be 
                          taken from vnfb */

			     /* first try taking it proportionally */

 if(rtdata[i].vnfb >= (whr_orig-rtdata[i].smb-rtdata[i].smb_nf)){
     
			         flowvs_add=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*(whr_orig-rtdata[i].smb-rtdata[i].smb_nf)/rtdata[i].vnfb;
			      

			       if(flowvs_add <= -ZEROVALUE) {flowvs_add = ZEROVALUE;}

			       newdata[it].flowvs[k][m]+=flowvs_add;
		     
			       newdata[it].vbh2[k][m]=flowvs_add*dstatic[k][m].vba*garea[k][m];
		     
			       rtdata[i].whr-=newdata[it].vbh2[k][m]/1000.;
	   rtdata[i].fh_vbh2+=newdata[it].vbh2[k][m]/1000.;
			     }
			     else{ /* not enough in vnfb, therefore take all remaining vnfb, the rest is unmet */


			       flowvs_add=newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m];
		       
			       if(flowvs_add <= -ZEROVALUE) {flowvs_add = ZEROVALUE;}
		       
			       newdata[it].flowvs[k][m]+=flowvs_add;
		     
			       newdata[it].vbh2[k][m]=flowvs_add*dstatic[k][m].vba*garea[k][m];
		     
			       rtdata[i].whr-=newdata[it].vbh2[k][m]/1000.;
rtdata[i].fh_vbh2+=newdata[it].vbh2[k][m]/1000.;
			     }


			   } /* ends else of smart_flow i[k][m]f */
		   
			 } /* ends else of fnf */
		 

		       }  /* end of gcode */
		     }
	       
		   } /* end of m */
		 } /* end of k */
		 
	       } /* end case 4*/
	 


	       for (k=0;k<NY;k++){
		 for (m=0;m<NX;m++){

#if REGION_TEST
	 
		   if(dstatic[k][m].rcode==REGION_TEST_GCODE){
		     
#else
		     if(dstatic[k][m].rcode > 0){
#endif
	     
		       if(dstatic[k][m].rcode == rdata[i].rcode){

		 
			   tester1+=newdata[it].vbh[k][m]/1000.;
			   tester2+=newdata[it].vbh2[k][m]/1000.;
			   tester3+=newdata[it].sbh[k][m]/1000.;
			   tester4+=newdata[it].sbh2[k][m]/1000.;
			   tester5+=newdata[it].sbh3[k][m]/1000.;

		   
		       }
		     }
		   }
		 }



		 testfile=fopen("wh.latest",outstat);
	   
		 
		 fprintf(testfile,"%d %s flag %d fao# %lf wh_needed %lf vbh %lf vbh2 %lf sbh %lf sbh2 %lf sbh3 %lf whtot %lf smb %lf smb_nf %lf vnfb %lf\n",curryear,cname[i],flag,rtdata[i].wh,whr_orig,tester1,tester2,tester3,tester4,tester5,(tester1+tester2+tester3+tester4+tester5),rtdata[i].smb,rtdata[i].smb_nf,rtdata[i].vnfb);
		   
		 
		 fclose(testfile);
	   

		 return;
	 
	       }

/***********************************************************************/
void force_harvest2(int it, int i, int curryear){

  FILE *testfile;
  int k, m, flag=0;
  double whr_orig=ZEROVALUE, flowvs_add=ZEROVALUE, flowss_add=ZEROVALUE, wh_from_sbh_only=ZEROVALUE;
  double tester1=ZEROVALUE, tester2=ZEROVALUE, tester3=ZEROVALUE, tester4=ZEROVALUE, tester5=ZEROVALUE;
 
  /* hold original amount of whr for spreading proportionally */

  whr_orig=ctdata[i].whr;



  /* compute country level smb, smb_nf, vnfb, s area for force harvest */
  
  
  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){
      
      
#if REGION_TEST
      
      if(dstatic[k][m].gcode==REGION_TEST_GCODE){
	   
#else
	if(dstatic[k][m].gcode > 0){
#endif
	     
	  if(dstatic[k][m].gcode == cdata[i].ccode){	   
	 	     
	    if(dstatic[k][m].fnf == 0){
	     
	      /* country level virgin nonforested biomass */
	      ctdata[i].vnfb+=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*garea[k][m]*dstatic[k][m].vba/1000.;
		 
	      /* country level secondary nonforested biomass */
	      ctdata[i].smb_nf+=(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*garea[k][m]*newdata[it].smb[k][m]/1000.;

	      /* country level secondary nonforested area */
	      ctdata[i].sarea_nf+=(newdata[it].s[k][m]+newdata[it].flowps[k][m]+newdata[it].flowcs[k][m]+newdata[it].flowvs[k][m]+newdata[it].flowus[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*garea[k][m];
 
	    }
	    else {

	      /* country level secondary forest biomass from sbh only */
	      wh_from_sbh_only+=newdata[it].sbh[k][m]/1000.;
	     
	      /* country level secondary forest biomass */
	      ctdata[i].smb+=(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*((garea[k][m]*newdata[it].smb[k][m]-newdata[it].sbh[k][m])/1000.); 

	     /* country level secondary forested area */
	      ctdata[i].sarea+=(newdata[it].s[k][m]+newdata[it].flowps[k][m]+newdata[it].flowcs[k][m]+newdata[it].flowvs[k][m]+newdata[it].flowus[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*garea[k][m];
	      
	    }
	
   
	  }  /* end of gcode */
	}
       
      } /* end of m */
    } /* end of k */
   




    /* 4 cases for getting remaining wood harvest (whr) */


    if(ctdata[i].smb >= ctdata[i].whr){


      flag=1;  /* case 1, taking needed younger secondary forest proportionally to satisfy whr */                                 

      for (k=0;k<NY;k++){
	for (m=0;m<NX;m++){
	 
#if REGION_TEST
	 
	  if(dstatic[k][m].gcode==REGION_TEST_GCODE){
	   
#else
	    if(dstatic[k][m].gcode > 0){
#endif
	     
	      if(dstatic[k][m].gcode == cdata[i].ccode){

	       
		if(dstatic[k][m].fnf == 1){
		 
		  if(ctdata[i].whr > ZEROVALUE){
		   
		   
		    newdata[it].sbh2[k][m]=whr_orig/ctdata[i].smb*(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*(garea[k][m]*newdata[it].smb[k][m]-newdata[it].sbh[k][m]);
		   
		    ctdata[i].whr-=newdata[it].sbh2[k][m]/1000.;
		   
		   
		  } /*end of whr */
		} /*end of fnf */
	      }  /* end of gcode */
	    }
	   
	  } /* end of m */
	} /* end of k */

       
      }/* end case 1*/
     

     
      else if((cdata[i].smart_flow_option == 1) && ((ctdata[i].smb+ctdata[i].vnfb) >= ctdata[i].whr)){    

	flag=2;  /* case 2, primary priority, not enough young secondary forest to satisfy demand, 
		   therefore taking all young secondary forest and needed virgin nonforest proportionally 
		   to satisfy whr */

  

	for (k=0;k<NY;k++){
	  for (m=0;m<NX;m++){
	   
#if REGION_TEST
	 
	    if(dstatic[k][m].gcode==REGION_TEST_GCODE){
	   
#else
	      if(dstatic[k][m].gcode > 0){
#endif
	     
		if(dstatic[k][m].gcode == cdata[i].ccode){

		 
		  if(ctdata[i].whr > ZEROVALUE){
		   
		    if(dstatic[k][m].fnf == 1) {

		      newdata[it].sbh2[k][m]=(newdata[it].smb[k][m]*garea[k][m]-newdata[it].sbh[k][m])*(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m]);

		      ctdata[i].whr-=newdata[it].sbh2[k][m]/1000.;
ctdata[i].fh_sbh2+=newdata[it].sbh2[k][m]/1000.;
		    } /* end of fnf*/		   
		  } /*end of whr*/
		}  /* end of gcode */
	      }
	     
	    } /* end of m */
	  } /* end of k */



	 /* take the whr proportionally out of vnfb; note needed to consider how much was taken
            out of secondary forest (ctdata[i].smb) */


	 for (k=0;k<NY;k++){
	   for (m=0;m<NX;m++){
	     
#if REGION_TEST
	     
	     if(dstatic[k][m].gcode==REGION_TEST_GCODE){
	   
#else
	       if(dstatic[k][m].gcode > 0){
#endif
	     
		 if(dstatic[k][m].gcode == cdata[i].ccode){

		 
		   if(ctdata[i].whr > ZEROVALUE){
		   
		     if(dstatic[k][m].fnf == 0) {
		     
	     	       flowvs_add=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*(whr_orig-ctdata[i].smb)/ctdata[i].vnfb; 
		      

		       if(flowvs_add <= -ZEROVALUE) {flowvs_add = ZEROVALUE;}
		     
		       newdata[it].flowvs[k][m]+=flowvs_add;
		     
		       newdata[it].vbh2[k][m]=flowvs_add*dstatic[k][m].vba*garea[k][m];
		     
		       ctdata[i].whr-=newdata[it].vbh2[k][m]/1000.;
		  ctdata[i].fh_vbh2+=newdata[it].vbh2[k][m]/1000.;   
		     } /*end of fnf */
		   
		    

		   
		   } /*end of whr*/
		 
		 }  /* end of gcode */
	       }
	     
	     } /* end of m */
	   } /* end of k */

	  
	 } /*end case 2 */
	 

	 
	 else if((cdata[i].smart_flow_option == 2) && ((ctdata[i].smb+ctdata[i].smb_nf) >= ctdata[i].whr)){  

	   flag=3;  /* case 3, secondary priority, taking all young secondary forest and needed secondary 
                   nonforest proportionally to satisfy whr */


     
	   for (k=0;k<NY;k++){
	     for (m=0;m<NX;m++){
	   
#if REGION_TEST
	 
	       if(dstatic[k][m].gcode==REGION_TEST_GCODE){
	   
#else
		 if(dstatic[k][m].gcode > 0){
#endif
	     
		   if(dstatic[k][m].gcode == cdata[i].ccode){

		 
		     if(ctdata[i].whr > ZEROVALUE){
		   
		       if(dstatic[k][m].fnf == 1) {

			 newdata[it].sbh2[k][m]=(newdata[it].smb[k][m]*garea[k][m]-newdata[it].sbh[k][m])*(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m]);

			 ctdata[i].whr-=newdata[it].sbh2[k][m]/1000.;

		       } /* end of fnf*/		   
		     } /*end of whr*/
		   }  /* end of gcode */
		 }
	     
	       } /* end of m */
	     } /* end of k */


	     /* take the whr proportionally out of secondary nonforest; note needed to consider how much was taken
		out of secondary forest (ctdata[i].smb) */


	     for (k=0;k<NY;k++){
	       for (m=0;m<NX;m++){
	     
#if REGION_TEST
	     
		 if(dstatic[k][m].gcode==REGION_TEST_GCODE){
	   
#else
		   if(dstatic[k][m].gcode > 0){
#endif
	     
		     if(dstatic[k][m].gcode == cdata[i].ccode){
		       
		       if(ctdata[i].whr > ZEROVALUE){
		   
			 if(dstatic[k][m].fnf == 0) {
		     
			    flowss_add=(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*(whr_orig-ctdata[i].smb)/ctdata[i].smb_nf;
		
		 
			   if(flowss_add <= -ZEROVALUE) {flowss_add = ZEROVALUE;}

			   newdata[it].sbh3[k][m]=flowss_add*newdata[it].smb[k][m]*garea[k][m];
		     
			   ctdata[i].whr-=newdata[it].sbh3[k][m]/1000.;


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

	    
#if REGION_TEST
	 
		   if(dstatic[k][m].gcode==REGION_TEST_GCODE){
	   
#else
		     if(dstatic[k][m].gcode > 0){
#endif
	     
		       if(dstatic[k][m].gcode == cdata[i].ccode){

		   
			 if(dstatic[k][m].fnf == 1){

			   /* take all young secondary forest */
		       
			   newdata[it].sbh2[k][m]=(newdata[it].smb[k][m]*garea[k][m]-newdata[it].sbh[k][m])*(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m]);

			   ctdata[i].whr-=newdata[it].sbh2[k][m]/1000.;
ctdata[i].fh_sbh2+=newdata[it].sbh2[k][m]/1000.;
			 }
		       }
		     }
		   }
		 }

		

	       for (k=0;k<NY;k++){
		 for (m=0;m<NX;m++){

	    
#if REGION_TEST
	 
		   if(dstatic[k][m].gcode==REGION_TEST_GCODE){
	   
#else
		     if(dstatic[k][m].gcode > 0){
#endif
	     
		       if(dstatic[k][m].gcode == cdata[i].ccode){

		   
			 if(dstatic[k][m].fnf == 0){

			   if(cdata[i].smart_flow_option == 1){

			     /* primary priority, first take all vnfb */


			     flowvs_add=newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m];
		       
			     if(flowvs_add <= -ZEROVALUE) {flowvs_add = ZEROVALUE;}
		       
			     newdata[it].flowvs[k][m]+=flowvs_add;
		     
			     newdata[it].vbh2[k][m]=flowvs_add*dstatic[k][m].vba*garea[k][m];
		     
			     ctdata[i].whr-=newdata[it].vbh2[k][m]/1000.;
ctdata[i].fh_vbh2+=newdata[it].vbh2[k][m]/1000.;
  


			     /* primary priority, not enough smb and vnfb, determined how much can be taken from smb_nf */

			     /* first try taking it proportionally */

 if(ctdata[i].smb_nf >= (whr_orig-ctdata[i].smb-ctdata[i].vnfb)){

		    
			         flowss_add=(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*(whr_orig-ctdata[i].smb-ctdata[i].vnfb)/ctdata[i].smb_nf; 
			      
			   
			       if(flowss_add <= -ZEROVALUE) {flowss_add = ZEROVALUE;}
			   
			       newdata[it].sbh3[k][m]=flowss_add*newdata[it].smb[k][m]*garea[k][m];
			   
			           ctdata[i].whr-=newdata[it].sbh3[k][m]/1000.;   
			        ctdata[i].fh_sbh3+=newdata[it].sbh3[k][m]/1000.; 
			    
			 			 
			     }
			     else{ /* not enough in smb_nf, therefore take all remaining smb_nf, the rest is unmet */
			 
			       flowss_add=newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m];
			   
			       if(flowss_add <= -ZEROVALUE) {flowss_add = ZEROVALUE;}
			   
			       newdata[it].sbh3[k][m]=flowss_add*newdata[it].smb[k][m]*garea[k][m];

			        ctdata[i].whr-=newdata[it].sbh3[k][m]/1000.; 
			        ctdata[i].fh_sbh3+=newdata[it].sbh3[k][m]/1000.; 
			     }


			   } 
			   else { /* if smart_flow_option == 2 */
		       
			     /* secondary priority, first take all smb_nf */

			     flowss_add=newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m];
			   
			     if(flowss_add <= -ZEROVALUE) {flowss_add = ZEROVALUE;}
			   
			     newdata[it].sbh3[k][m]=flowss_add*newdata[it].smb[k][m]*garea[k][m];

			     ctdata[i].whr-=newdata[it].sbh3[k][m]/1000.;
ctdata[i].fh_sbh3+=newdata[it].sbh3[k][m]/1000.;


			     /* secondary priority, not enough smb and smb_nf, determined how much can be 
                          taken from vnfb */

			     /* first try taking it proportionally */

 if(ctdata[i].vnfb >= (whr_orig-ctdata[i].smb-ctdata[i].smb_nf)){
		    
			         flowvs_add=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*(whr_orig-ctdata[i].smb-ctdata[i].smb_nf)/ctdata[i].vnfb;
			       

			       if(flowvs_add <= -ZEROVALUE) {flowvs_add = ZEROVALUE;}

			       newdata[it].flowvs[k][m]+=flowvs_add;
		     
			       newdata[it].vbh2[k][m]=flowvs_add*dstatic[k][m].vba*garea[k][m];
		     
			       ctdata[i].whr-=newdata[it].vbh2[k][m]/1000.;
	   ctdata[i].fh_vbh2+=newdata[it].vbh2[k][m]/1000.;
			     }
			     else{ /* not enough in vnfb, therefore take all remaining vnfb, the rest is unmet */


			       flowvs_add=newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m];
		       
			       if(flowvs_add <= -ZEROVALUE) {flowvs_add = ZEROVALUE;}
		       
			       newdata[it].flowvs[k][m]+=flowvs_add;
		     
			       newdata[it].vbh2[k][m]=flowvs_add*dstatic[k][m].vba*garea[k][m];
		     
			       ctdata[i].whr-=newdata[it].vbh2[k][m]/1000.;
ctdata[i].fh_vbh2+=newdata[it].vbh2[k][m]/1000.;
			     }


			   } /* ends else of smart_flow i[k][m]f */
		   
			 } /* ends else of fnf */
		 

		       }  /* end of gcode */
		     }
	       
		   } /* end of m */
		 } /* end of k */
		
	  	   
	       } /* end case 4*/
	 


	       for (k=0;k<NY;k++){
		 for (m=0;m<NX;m++){

#if REGION_TEST
	 
		   if(dstatic[k][m].gcode==REGION_TEST_GCODE){
		     
#else
		     if(dstatic[k][m].gcode > 0){
#endif
	     
		       if(dstatic[k][m].gcode == cdata[i].ccode){

		 
			   tester1+=newdata[it].vbh[k][m]/1000.;
			   tester2+=newdata[it].vbh2[k][m]/1000.;
			   tester3+=newdata[it].sbh[k][m]/1000.;
			   tester4+=newdata[it].sbh2[k][m]/1000.;
			   tester5+=newdata[it].sbh3[k][m]/1000.;

		   
		       }
		     }
		   }
		 }

		 testfile=fopen("wh.latest",outstat);
	   
		 
		 fprintf(testfile,"%d %s flag %d fao# %lf wh_needed %lf vbh %lf vbh2 %lf sbh %lf sbh2 %lf sbh3 %lf whtot %lf smb %lf smb_nf %lf vnfb %lf\n",curryear,cname[i],flag,ctdata[i].wh,whr_orig,tester1,tester2,tester3,tester4,tester5,(tester1+tester2+tester3+tester4+tester5),ctdata[i].smb,ctdata[i].smb_nf,ctdata[i].vnfb);
		   
		 
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
      
#if REGION_TEST
	
      if(dstatic[k][m].rcode==REGION_TEST_GCODE){
	   
#else
	if(dstatic[k][m].rcode > 0){
#endif
	     
	  i=dstatic[k][m].newrcode;
	  iz=newdata[it].zdis[k][m];
  

	  if(dstatic[k][m].fnf == 1){
		
		
	    /* vb in units = MgC */
	    rztdata[i][iz].vb+=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*garea[k][m]*dstatic[k][m].vba/1000.0;
		
	  }

	}  
      } /* end of m */
    } /* end of k */
    
    
    return;
  }

/********************************************************************/
void update_vb2(int it){
  
  int iz, k, m, i;
  
  
  /* determine virgin biomass at the country level (vb); vb is the virgin
     available for a particular country for particular iz; therefore given
     a country code and a value of iz, vb is known */ 
  
  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){
      
#if REGION_TEST
	
      if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){
	   
#else
	if(dstatic[k][m].gcode > 0){
#endif
	     
	  i=dstatic[k][m].newgcode;
	  iz=newdata[it].zdis[k][m];
  

	  if(dstatic[k][m].fnf == 1){
		
		
	    /* vb in units = MgC */
	    cztdata[i][iz].vb+=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*garea[k][m]*dstatic[k][m].vba/1000.0;
		
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
      if(dstatic[k][m].rcode > 0){
      i=dstatic[k][m].newrcode;

      
      if(rdata[i].zdis_option == 1){

  	/* setting zdis=0 for land use cells and masking out ocean */
      
      
	if(newdata[it].c[k][m] > ZEROVALUE || newdata[it].p[k][m] > ZEROVALUE ||
	   newdata[it].s[k][m] > ZEROVALUE || newdata[it].u[k][m] > ZEROVALUE) {
	
	  newdata[it].zdis[k][m] = 0;
	
	}  
	

      }

      else{
	newdata[it].zdis[k][m]  = 0;
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
	

	    if((k-1 > -1) && (m-1 > -1)) {
	      if(newdata[it].zdis[k][m] == zcheck) {
		if(newdata[it].zdis[k-1][m-1] > zcheck) {
 		  newdata[it].zdis[k-1][m-1] = zvalue;
		}
	      }
	    }
	    
	    if((k-1 > -1) && (m == m)) {
	      if(newdata[it].zdis[k][m] == zcheck) {
		if(newdata[it].zdis[k-1][m] > zcheck) {
                 
		  newdata[it].zdis[k-1][m] = zvalue;
		}
	      }
	    }
	
	    if((k-1 > -1) && (m+1 < NX)) {
	      if(newdata[it].zdis[k][m] == zcheck) {
		if(newdata[it].zdis[k-1][m+1] > zcheck) {
		  newdata[it].zdis[k-1][m+1] = zvalue;
		}
	      }
	    }
	    
	    if((k == k) && (m-1 > -1)) {
	      if(newdata[it].zdis[k][m] == zcheck) {
		if(newdata[it].zdis[k][m-1] > zcheck) {
		  newdata[it].zdis[k][m-1] = zvalue;
		}
	      }
	    }
	    
	    if((k == k) && (m+1 < NX)) {
	      if(newdata[it].zdis[k][m] == zcheck) {
		if(newdata[it].zdis[k][m+1] > zcheck) {
		  newdata[it].zdis[k][m+1] = zvalue;
		}
	      }
	    }
	    
	    if((k+1 < NY) && (m-1 > -1)) {
	      if(newdata[it].zdis[k][m] == zcheck) {
		if(newdata[it].zdis[k+1][m-1] > zcheck) {
		  newdata[it].zdis[k+1][m-1] = zvalue;
		}
	      }
	    }
	    
	    if((k+1 < NY) && (m == m)) {
	      if(newdata[it].zdis[k][m] == zcheck) {
		if(newdata[it].zdis[k+1][m] > zcheck) {
		  newdata[it].zdis[k+1][m] = zvalue;
		}
	      }
	    }
	    
	    if((k+1 < NY) && (m+1 < NX)){ 
	      if(newdata[it].zdis[k][m] == zcheck) {
		if(newdata[it].zdis[k+1][m+1] > zcheck) {
		  newdata[it].zdis[k+1][m+1] = zvalue;
		}
	      }      
	    }
	    
	    /* }   end if gcode */

	}   /* end of m */
      }    /* end of k */


    }     /* end of iz */

	
  return;
}

/***********************************************************************/
void zdis_calc2(int it){
  
  FILE *testfile;
  int i;
  int iz, k, m, zcheck, zvalue, flagger; 
  

  
  /* determine focal cell for zdis, zdis is initialized to NZ earlier for all times*/ 
   
  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){
      
      i=dstatic[k][m].newgcode;

      
      if(cdata[i].zdis_option == 1){

  	/* setting zdis=0 for land use cells and masking out ocean */
      
      
	if(newdata[it].c[k][m] > ZEROVALUE || newdata[it].p[k][m] > ZEROVALUE ||
	   newdata[it].s[k][m] > ZEROVALUE) {
	
	  newdata[it].zdis[k][m]  = 0;
	
	}  


      }

      else{
	newdata[it].zdis[k][m]  = 0;
      }
      
    }
  }
  
  
  /* check all other cells relative to focal cell */
  
  
  for (iz=0;iz<NZ;iz++){
    
    zcheck=iz;
    zvalue=iz+1;
    
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){
	

	    if((k-1 > -1) && (m-1 > -1)) {
	      if(newdata[it].zdis[k-1][m] == zcheck) {
		if(newdata[it].zdis[k-1][m-1] > zcheck) {
 		  newdata[it].zdis[k-1][m-1] = zvalue;
		}
	      }
	    }
	    
	    if((k-1 > -1) && (m == m)) {
	      if(newdata[it].zdis[k][m] == zcheck) {
		if(newdata[it].zdis[k-1][m] > zcheck) {
                 
		  newdata[it].zdis[k-1][m] = zvalue;
		}
	      }
	    }
	
	    if((k-1 > -1) && (m+1 < NX)) {
	      if(newdata[it].zdis[k][m] == zcheck) {
		if(newdata[it].zdis[k-1][m+1] > zcheck) {
		  newdata[it].zdis[k-1][m+1] = zvalue;
		}
	      }
	    }
	    
	    if((k == k) && (m-1 > -1)) {
	      if(newdata[it].zdis[k][m] == zcheck) {
		if(newdata[it].zdis[k][m-1] > zcheck) {
		  newdata[it].zdis[k][m-1] = zvalue;
		}
	      }
	    }
	    
	    if((k == k) && (m+1 < NX)) {
	      if(newdata[it].zdis[k][m] == zcheck) {
		if(newdata[it].zdis[k][m+1] > zcheck) {
		  newdata[it].zdis[k][m+1] = zvalue;
		}
	      }
	    }
	    
	    if((k+1 < NY) && (m-1 > -1)) {
	      if(newdata[it].zdis[k][m] == zcheck) {
		if(newdata[it].zdis[k+1][m-1] > zcheck) {
		  newdata[it].zdis[k+1][m-1]= zvalue;
		}
	      }
	    }
	    
	    if((k+1 < NY) && (m == m)) {
	      if(newdata[it].zdis[k][m] == zcheck) {
		if(newdata[it].zdis[k+1][m] > zcheck) {
		  newdata[it].zdis[k+1][m] = zvalue;
		}
	      }
	    }
	    
	    if((k+1 < NY) && (m+1 < NX)){ 
	      if(newdata[it].zdis[k][m] == zcheck) {
		if(newdata[it].zdis[k+1][m+1] > zcheck) {
		  newdata[it].zdis[k+1][m+1] = zvalue;
		}
	      }      
	    }
	    
	    /* }   end if gcode */

	}   /* end of m */
      }    /* end of k */     


    }     /* end of iz */
    
	
  return;
}

/***********************************************************************/
void update_states(int it, int zmax){

  int i,k,m, iz;
  float sma_area_gained, sma_area_lost, sma_area_notlost;
  float sum_test;
  double gsum=ZEROVALUE, tester=ZEROVALUE;



  
  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){


#if REGION_TEST
	 
	 if(dstatic[k][m].rcode==REGION_TEST_GCODE){
	   
#else
	   if(dstatic[k][m].rcode > 0){
#endif
	     
	       

	     if (FLOW_BUG_PRINT){


    	/* checking for negative zero flows */

	if(newdata[it].flowcp[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowcp<0, resetting to 0 by force %10.15lf\n",newdata[it].flowcp[k][m]);newdata[it].flowcp[k][m] = ZEROVALUE;}
	if(newdata[it].flowpc[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowpc<0, resetting to 0 by force %10.15lf\n",newdata[it].flowpc[k][m]);newdata[it].flowpc[k][m] = ZEROVALUE;}
	if(newdata[it].flowpv[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowpv<0, resetting to 0 by force %10.15lf\n",newdata[it].flowpv[k][m]);newdata[it].flowpv[k][m] = ZEROVALUE;}
	if(newdata[it].flowvp[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowvp<0, resetting to 0 by force %10.15lf\n",newdata[it].flowvp[k][m]);newdata[it].flowvp[k][m] = ZEROVALUE;}
	if(newdata[it].flowvc[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowvc<0, resetting to 0 by force %10.15lf\n",newdata[it].flowvc[k][m]);newdata[it].flowvc[k][m] = ZEROVALUE;}
	if(newdata[it].flowcv[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowcv<0, resetting to 0 by force %10.15lf\n",newdata[it].flowcv[k][m]);newdata[it].flowcv[k][m] = ZEROVALUE;}
	if(newdata[it].flowsp[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowsp<0, resetting to 0 by force %10.15lf\n",newdata[it].flowsp[k][m]);newdata[it].flowsp[k][m] = ZEROVALUE;}
	if(newdata[it].flowps[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowps<0, resetting to 0 by force %10.15lf\n",newdata[it].flowps[k][m]);newdata[it].flowps[k][m] = ZEROVALUE;}
	if(newdata[it].flowsc[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowsc<0, resetting to 0 by force %10.15lf\n",newdata[it].flowsc[k][m]);newdata[it].flowsc[k][m] = ZEROVALUE;}
	if(newdata[it].flowcs[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowcs<0, resetting to 0 by force %10.15lf\n",newdata[it].flowcs[k][m]);newdata[it].flowcs[k][m] = ZEROVALUE;}
	if(newdata[it].flowvs[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowvs<0, resetting to 0 by force %10.15lf\n",newdata[it].flowvs[k][m]);newdata[it].flowvs[k][m] = ZEROVALUE;}
	if(newdata[it].flowcu[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowcu<0, resetting to 0 by force %10.15lf\n",newdata[it].flowcu[k][m]);newdata[it].flowcu[k][m] = ZEROVALUE;}
	if(newdata[it].flowpu[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowpu<0, resetting to 0 by force %10.15lf\n",newdata[it].flowpu[k][m]);newdata[it].flowpu[k][m] = ZEROVALUE;}
	if(newdata[it].flowvu[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowvu<0, resetting to 0 by force %10.15lf\n",newdata[it].flowvu[k][m]);newdata[it].flowvu[k][m] = ZEROVALUE;}
	if(newdata[it].flowsu[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowsu<0, resetting to 0 by force %10.15lf\n",newdata[it].flowsu[k][m]);newdata[it].flowsu[k][m] = ZEROVALUE;}
	if(newdata[it].flowuc[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowuc<0, resetting to 0 by force %10.15lf\n",newdata[it].flowuc[k][m]);newdata[it].flowuc[k][m] = ZEROVALUE;}
	if(newdata[it].flowup[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowup<0, resetting to 0 by force %10.15lf\n",newdata[it].flowup[k][m]);newdata[it].flowup[k][m] = ZEROVALUE;}
	if(newdata[it].flowus[k][m] < ZEROVALUE) {printf("Update states flow Bug: flowus<0, resetting to 0 by force %10.15lf\n",newdata[it].flowus[k][m]);newdata[it].flowus[k][m] = ZEROVALUE;}
	


	if(newdata[it].vbh[k][m] < ZEROVALUE) {printf("Harvest Bug: vbh<0, resetting to 0 by force %10.15lf k: %d m: %d\n",newdata[it].vbh[k][m],k,m);newdata[it].vbh[k][m] = ZEROVALUE;}
	if(newdata[it].sbh[k][m] < ZEROVALUE) {printf("Harvest Bug: sbh<0, resetting to 0 by force %10.15lf\n",newdata[it].sbh[k][m]);newdata[it].sbh[k][m] = ZEROVALUE;}
	if(newdata[it].vbh2[k][m] < ZEROVALUE) {printf("Harvest Bug: vbh2<0, resetting to 0 by force %10.15lf\n",newdata[it].vbh2[k][m]);newdata[it].vbh2[k][m] = ZEROVALUE;}
	if(newdata[it].sbh2[k][m] < ZEROVALUE) {printf("Harvest Bug: sbh2<0, resetting to 0 by force %10.15lf\n",newdata[it].sbh2[k][m]);newdata[it].sbh2[k][m] = ZEROVALUE;}



	     }else{

	if(newdata[it].flowcp[k][m] < ZEROVALUE) {newdata[it].flowcp[k][m] = ZEROVALUE;}
	if(newdata[it].flowpc[k][m] < ZEROVALUE) {newdata[it].flowpc[k][m] = ZEROVALUE;}
	if(newdata[it].flowpv[k][m] < ZEROVALUE) {newdata[it].flowpv[k][m] = ZEROVALUE;}
	if(newdata[it].flowvp[k][m] < ZEROVALUE) {newdata[it].flowvp[k][m] = ZEROVALUE;}
	if(newdata[it].flowvc[k][m] < ZEROVALUE) {newdata[it].flowvc[k][m] = ZEROVALUE;}
	if(newdata[it].flowcv[k][m] < ZEROVALUE) {newdata[it].flowcv[k][m] = ZEROVALUE;}
	if(newdata[it].flowsp[k][m] < ZEROVALUE) {newdata[it].flowsp[k][m] = ZEROVALUE;}
	if(newdata[it].flowps[k][m] < ZEROVALUE) {newdata[it].flowps[k][m] = ZEROVALUE;}
	if(newdata[it].flowsc[k][m] < ZEROVALUE) {newdata[it].flowsc[k][m] = ZEROVALUE;}
	if(newdata[it].flowcs[k][m] < ZEROVALUE) {newdata[it].flowcs[k][m] = ZEROVALUE;}

	if(newdata[it].flowvs[k][m] < ZEROVALUE) {newdata[it].flowvs[k][m] = ZEROVALUE;}
	if(newdata[it].flowcu[k][m] < ZEROVALUE) {newdata[it].flowcu[k][m] = ZEROVALUE;}
	if(newdata[it].flowpu[k][m] < ZEROVALUE) {newdata[it].flowpu[k][m] = ZEROVALUE;}
	if(newdata[it].flowvu[k][m] < ZEROVALUE) {newdata[it].flowvu[k][m] = ZEROVALUE;}
	if(newdata[it].flowsu[k][m] < ZEROVALUE) {newdata[it].flowsu[k][m] = ZEROVALUE;}
	if(newdata[it].flowuc[k][m] < ZEROVALUE) {newdata[it].flowuc[k][m] = ZEROVALUE;}
	if(newdata[it].flowup[k][m] < ZEROVALUE) {newdata[it].flowup[k][m] = ZEROVALUE;}
	if(newdata[it].flowus[k][m] < ZEROVALUE) {newdata[it].flowus[k][m] = ZEROVALUE;}

	     }


	/*updating states*/
	

	newdata[it+1].v[k][m]=newdata[it].v[k][m]+newdata[it].flowpv[k][m]+
	  newdata[it].flowcv[k][m]-newdata[it].flowvp[k][m]-
	  newdata[it].flowvc[k][m]-newdata[it].flowvs[k][m]-newdata[it].flowvu[k][m];
	
	newdata[it+1].s[k][m]=newdata[it].s[k][m]+newdata[it].flowps[k][m]+
	  newdata[it].flowvs[k][m]+newdata[it].flowcs[k][m]-
	  newdata[it].flowsp[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsu[k][m]+newdata[it].flowus[k][m];
	   

	newdata[it+1].c[k][m]=newdata[it].c[k][m] +
	  newdata[it].flowpc[k][m] + newdata[it].flowvc[k][m] +
	  newdata[it].flowsc[k][m] - newdata[it].flowcp[k][m] -
	  newdata[it].flowcv[k][m] - newdata[it].flowcs[k][m] -
	  newdata[it].flowcu[k][m] + newdata[it].flowuc[k][m];
	
	newdata[it+1].p[k][m]=newdata[it].p[k][m] +
	  newdata[it].flowcp[k][m] + newdata[it].flowvp[k][m] +
	  newdata[it].flowsp[k][m] - newdata[it].flowpc[k][m] -
	  newdata[it].flowpv[k][m] - newdata[it].flowps[k][m] -
	  newdata[it].flowpu[k][m] + newdata[it].flowup[k][m];

	newdata[it+1].u[k][m]=newdata[it].u[k][m] +
	  newdata[it].flowcu[k][m] + newdata[it].flowvu[k][m] +
	  newdata[it].flowsu[k][m] + newdata[it].flowpu[k][m] - 
          newdata[it].flowuc[k][m] - newdata[it].flowus[k][m] - 
          newdata[it].flowup[k][m];

	if (STATE_BUG_PRINT){
      
	if(newdata[it+1].v[k][m] < ZEROVALUE) { 
           printf("Update states Bug: v<0 at t+1, resetting to 0 by force. k=%d m=%d it=%d it+1=%d v(t)=%lf v(t+1)=%lf pv %lf cv %lf vp %lf vc %lf vs %lf\n",k,m,it,it+1,newdata[it].v[k][m],newdata[it+1].v[k][m],newdata[it].flowpv[k][m],newdata[it].flowcv[k][m],newdata[it].flowvp[k][m],newdata[it].flowvc[k][m],newdata[it].flowvs[k][m]);                
           newdata[it+1].v[k][m] = ZEROVALUE;
        }
	if(newdata[it+1].s[k][m] < ZEROVALUE) {         
           printf("Update states Bug: s<0 at t+1, resetting to 0 by force. k=%d m=%d it=%d it+1=%d s(t)=%lf s(t+1)=%lf cs %lf ps %lf vs %lf sc %lf sp %lf\n",k,m,it,it+1,newdata[it].s[k][m],newdata[it+1].s[k][m],newdata[it].flowcs[k][m],newdata[it].flowps[k][m],newdata[it].flowvs[k][m],newdata[it].flowsc[k][m],newdata[it].flowsp[k][m]);
           newdata[it+1].s[k][m] =ZEROVALUE;  
        }
	if(newdata[it+1].c[k][m] < ZEROVALUE) {   
           printf("Update states Bug: c<0 at t+1, resetting to 0 by force. k=%d m=%d it=%d it+1=%d c(t)=%lf c(t+1)=%lf pc %lf vc %lf sc %lf cp %lf cv %lf cs %lf\n",k,m,it,it+1,newdata[it].c[k][m],newdata[it+1].c[k][m],newdata[it].flowpc[k][m],newdata[it].flowvc[k][m],newdata[it].flowsc[k][m],newdata[it].flowcp[k][m],newdata[it].flowcv[k][m],newdata[it].flowcs[k][m]);
           newdata[it+1].c[k][m] = ZEROVALUE; 
        }
	if(newdata[it+1].p[k][m] < ZEROVALUE) {
           printf("Update states Bug: p<0 at t+1, resetting to 0 by force. k=%d m=%d it=%d it+1=%d p(t)=%lf p(t+1)=%lf cp %lf vp %lf sp %lf pc %lf pv %lf ps %lf\n",k,m,it,it+1,newdata[it].p[k][m],newdata[it+1].p[k][m],newdata[it].flowcp[k][m],newdata[it].flowvp[k][m],newdata[it].flowsp[k][m],newdata[it].flowpc[k][m],newdata[it].flowpv[k][m],newdata[it].flowps[k][m]);
           newdata[it+1].p[k][m] = ZEROVALUE; 
	}
	if(newdata[it+1].u[k][m] < ZEROVALUE) {
           printf("Update states Bug: u<0 at t+1, resetting to 0 by force. k=%d\
 m=%d it=%d it+1=%d u(t)=%lf u(t+1)=%lf pu %lf vu %lf su %lf up %lf cu %lf us %\
lf\n",k,m,it,it+1,newdata[it].u[k][m],newdata[it+1].u[k][m],newdata[it].flowpu[k][m],newdata[it].flowvu[k][m],newdata[it].flowsu[k][m],newdata[it].flowup[k][m],newdata[it].flowcu[k][m],newdata[it].flowus[k][m]);
           newdata[it+1].u[k][m] = ZEROVALUE; 

        }


	}else{

	if(newdata[it+1].v[k][m] < ZEROVALUE) newdata[it+1].v[k][m] = ZEROVALUE;  
	if(newdata[it+1].s[k][m] < ZEROVALUE) newdata[it+1].s[k][m] = ZEROVALUE;
	if(newdata[it+1].c[k][m] < ZEROVALUE) newdata[it+1].c[k][m] = ZEROVALUE;
	if(newdata[it+1].p[k][m] < ZEROVALUE) newdata[it+1].p[k][m] = ZEROVALUE; 
	if(newdata[it+1].u[k][m] < ZEROVALUE) newdata[it+1].u[k][m] = ZEROVALUE; 
	
	}



	sum_test=ZEROVALUE;
	
	sum_test=newdata[it+1].c[k][m]+newdata[it+1].p[k][m]+
                 newdata[it+1].v[k][m]+newdata[it+1].s[k][m]+
                 newdata[it+1].i[k][m]+newdata[it+1].w[k][m]+newdata[it+1].u[k][m];



        if(newdata[it].smb[k][m] > ZEROVALUE_CHECK){

	  sma_area_gained=newdata[it].flowvs[k][m]+newdata[it].flowcs[k][m]+newdata[it].flowps[k][m]+newdata[it].flowus[k][m]+(newdata[it].sbh[k][m]+newdata[it].sbh2[k][m]+newdata[it].sbh3[k][m])/newdata[it].smb[k][m]/garea[k][m];
	  
	  sma_area_lost=newdata[it].flowsp[k][m]+newdata[it].flowsc[k][m]+newdata[it].flowsu[k][m]+
	    (newdata[it].sbh[k][m]+newdata[it].sbh2[k][m]+newdata[it].sbh3[k][m])/newdata[it].smb[k][m]/garea[k][m];

	} 
        else {
	  
	  sma_area_gained=newdata[it].flowvs[k][m]+newdata[it].flowcs[k][m]+newdata[it].flowps[k][m]+newdata[it].flowus[k][m];
	  sma_area_lost=newdata[it].flowsp[k][m]+newdata[it].flowsc[k][m]+newdata[it].flowsu[k][m];
	}


	sma_area_notlost=newdata[it].s[k][m]-sma_area_lost;
	

 
	if((newdata[it].s[k][m]+sma_area_gained-sma_area_lost) > ZEROVALUE_CHECK){
	  

	  newdata[it+1].sma[k][m]=(sma_area_notlost*(newdata[it].sma[k][m]+1.0)+
	               sma_area_gained*1.0)/ 
                      (newdata[it].s[k][m]+sma_area_gained-sma_area_lost);


	}
	else{
	 
	  newdata[it+1].sma[k][m]=1.0;
	}



	
	if(dstatic[k][m].vba > ZEROVALUE){
	  
	  /* 0.38 is wood fraction of NPP, 0.75 is the aboveground fraction of NPP */


	  newdata[it+1].smb[k][m]=dstatic[k][m].vba*(1.0000000-exp(-(dstatic[k][m].vnppa*0.75*0.38*newdata[it].sma[k][m])/dstatic[k][m].vba));   


	}
	else{
	  newdata[it+1].smb[k][m]=ZEROVALUE;
	}
	
	
	if(newdata[it+1].sma[k][m] < ZEROVALUE) {
	  printf("Bug: sma < 0, being reset to zero by force. k=%d m=%d it+1=%d sma(t) %10.15lf sma(t+1) %10.15lf sma_area_notlost %lf s(t) %lf sma_area_lost %lf sma_area_gained %lf smb(t) %lf\n",k,m,it+1,newdata[it].sma[k][m],newdata[it+1].sma[k][m],sma_area_notlost,newdata[it].s[k][m],sma_area_lost,sma_area_gained,newdata[it].smb[k][m]);
	   
	   newdata[it+1].sma[k][m] = 1.0;
	}
	
	 
	   }

         } /* end of m */
    } /* end of k */



  return;
}





/********************************************************************/
void output_lu(int curryear){

  FILE *outfile;
  int i, ih, k, m, itstart;
  char latbase[4], lonbase[4];
  //  char alatshort[4], alatmedium[5], alatlong[6];
  //  char alonshorter[4], alonshort[5], alonmedium[6], alonlong[7];
  char alatshort[5], alatmedium[6], alatlong[7];
  char alonshorter[5], alonshort[6], alonmedium[7], alonlong[8];
  char fname[170], fout_trans[170];
  double flowsbh, flowvbh, flowsbh2, flowvbh2, flowsbh3;

  
  strcpy(lonbase,"lon");
  strcpy(latbase,"lat");
  
  
  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){
      
      strcpy(fname,"");
      /*      strcpy(fout_trans,""); */
      
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
      
#if REGION_TEST
	  
	  if(dstatic[k][m].rcode==REGION_TEST_GCODE){    	
#else
	    if(dstatic[k][m].rcode > 0){
#endif	
      
	
	outfile=fopen(fname,outstat);
	
	if(print_file_headers) {
	  if (use_urban) {
	    fprintf(outfile,"year cp pc pv vp vc cv sc cs sp ps vs cu pu vu su uc up us sbh f_sbh vbh f_vbh sbh2 f_sbh2 vbh2 f_vbh2 sbh3 f_sbh3\n");
	  }else{
          fprintf(outfile,"year cp pc pv vp vc cv sc cs sp ps vs sbh f_sbh vbh f_vbh sbh2 f_sbh2 vbh2 f_vbh2 sbh3 f_sbh3\n");
	  }
	}

        if(newdata[0].flowcp[k][m] <= -ZEROVALUE) {newdata[0].flowcp[k][m] = ZEROVALUE;}
        if(newdata[0].flowpc[k][m] <= -ZEROVALUE) {newdata[0].flowpc[k][m] = ZEROVALUE;}
        if(newdata[0].flowpv[k][m] <= -ZEROVALUE) {newdata[0].flowpv[k][m] = ZEROVALUE;}        
        if(newdata[0].flowvp[k][m] <= -ZEROVALUE) {newdata[0].flowvp[k][m] = ZEROVALUE;}        
        if(newdata[0].flowvc[k][m] <= -ZEROVALUE) {newdata[0].flowvc[k][m] = ZEROVALUE;}        
        if(newdata[0].flowcv[k][m] <= -ZEROVALUE) {newdata[0].flowcv[k][m] = ZEROVALUE;}        
        if(newdata[0].flowsp[k][m] <= -ZEROVALUE) {newdata[0].flowsp[k][m] = ZEROVALUE;}        
        if(newdata[0].flowps[k][m] <= -ZEROVALUE) {newdata[0].flowps[k][m] = ZEROVALUE;}        
        if(newdata[0].flowsc[k][m] <= -ZEROVALUE) {newdata[0].flowsc[k][m] = ZEROVALUE;}        
        if(newdata[0].flowcs[k][m] <= -ZEROVALUE) {newdata[0].flowcs[k][m] = ZEROVALUE;}
        if(newdata[0].flowvs[k][m] <= -ZEROVALUE) {newdata[0].flowvs[k][m] = ZEROVALUE;}
        if(newdata[0].flowcu[k][m] <= -ZEROVALUE) {newdata[0].flowcu[k][m] = ZEROVALUE;}
        if(newdata[0].flowpu[k][m] <= -ZEROVALUE) {newdata[0].flowpu[k][m] = ZEROVALUE;}
        if(newdata[0].flowvu[k][m] <= -ZEROVALUE) {newdata[0].flowvu[k][m] = ZEROVALUE;}
        if(newdata[0].flowsu[k][m] <= -ZEROVALUE) {newdata[0].flowsu[k][m] = ZEROVALUE;}
        if(newdata[0].flowuc[k][m] <= -ZEROVALUE) {newdata[0].flowuc[k][m] = ZEROVALUE;}
        if(newdata[0].flowup[k][m] <= -ZEROVALUE) {newdata[0].flowup[k][m] = ZEROVALUE;}
        if(newdata[0].flowus[k][m] <= -ZEROVALUE) {newdata[0].flowus[k][m] = ZEROVALUE;}

	flowsbh=flowsbh2=flowvbh=flowvbh2=flowsbh3=ZEROVALUE;
	
	if(newdata[0].smb[k][m] > ZEROVALUE) flowsbh=newdata[0].sbh[k][m]/newdata[0].smb[k][m]/garea[k][m];
	if(dstatic[k][m].vba > ZEROVALUE) flowvbh=newdata[0].vbh[k][m]/dstatic[k][m].vba/garea[k][m];
	if(newdata[0].smb[k][m] > ZEROVALUE) flowsbh2=newdata[0].sbh2[k][m]/newdata[0].smb[k][m]/garea[k][m];
	if(dstatic[k][m].vba > ZEROVALUE) flowvbh2=newdata[0].vbh2[k][m]/dstatic[k][m].vba/garea[k][m];
	if(newdata[0].smb[k][m] > ZEROVALUE) flowsbh3=newdata[0].sbh3[k][m]/newdata[0].smb[k][m]/garea[k][m];

        if (use_urban) {

          fprintf(outfile,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                  curryear,
		  newdata[0].flowcp[k][m],
		  newdata[0].flowpc[k][m],
		  newdata[0].flowpv[k][m],
		  newdata[0].flowvp[k][m],
		  newdata[0].flowvc[k][m],
		  newdata[0].flowcv[k][m],
		  newdata[0].flowsc[k][m],
		  newdata[0].flowcs[k][m],
		  newdata[0].flowsp[k][m],
		  newdata[0].flowps[k][m],
                  newdata[0].flowvs[k][m],
                  newdata[0].flowcu[k][m],
                  newdata[0].flowpu[k][m],
                  newdata[0].flowvu[k][m],
                  newdata[0].flowsu[k][m],
                  newdata[0].flowuc[k][m],
                  newdata[0].flowup[k][m],
                  newdata[0].flowus[k][m],
		  newdata[0].sbh[k][m]/1000,
		  flowsbh,
                  newdata[0].vbh[k][m]/1000, 
		  flowvbh,
		  newdata[0].sbh2[k][m]/1000,
		  flowsbh2,
                  newdata[0].vbh2[k][m]/1000,
		  flowvbh2,
                  newdata[0].sbh3[k][m]/1000,
		  flowsbh3);
	}else{
      
          fprintf(outfile,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                  curryear,
		  newdata[0].flowcp[k][m],
		  newdata[0].flowpc[k][m],
		  newdata[0].flowpv[k][m],
		  newdata[0].flowvp[k][m],
		  newdata[0].flowvc[k][m],
		  newdata[0].flowcv[k][m],
		  newdata[0].flowsc[k][m],
		  newdata[0].flowcs[k][m],
		  newdata[0].flowsp[k][m],
		  newdata[0].flowps[k][m],
                  newdata[0].flowvs[k][m],
		  newdata[0].sbh[k][m]/1000,
		  flowsbh,
                  newdata[0].vbh[k][m]/1000, 
		  flowvbh,
		  newdata[0].sbh2[k][m]/1000,
		  flowsbh2,
                  newdata[0].vbh2[k][m]/1000,
		  flowvbh2,
                  newdata[0].sbh3[k][m]/1000,
		  flowsbh3);
	}
       fclose(outfile);
      }
      
    }
  }
  return;
}
/********************************************************************/
void output_updated_states_nc(int curryear){
  int stat,k,m,timeidx;
/*   char ayear[7], state_file[256]; */
  char state_file[256];
  static long start[] = {0, 0, 0};   /* start at first value */
  static long count[] = {1, NY, NX};


    /* dimension ids */

  int lat_dim;
  int lon_dim;
  int time_dim;

    /* dimension lengths */

  size_t lat_len = NY;
  size_t lon_len = NX;  
  
  size_t time_len = NC_UNLIMITED;

    /* Variable Ids */
  int cell_area_id;
  int gcrop_id;
  int gothr_id;
  int gurbn_id;
  int gpast_id;
  int gsecd_id;
  int gssma_id;
  int gssmb_id;
  int gsumm_id;
  int gflcp_id;
  int gflcs_id;
  int gflpc_id;
  int gflps_id;
  int gflsc_id;
  int gflsp_id;
  int gflvc_id;
  int gflvp_id;
  int gfsh1_id;
  int gfsh2_id;
  int gfsh3_id;
  int gfvh1_id;
  int gfvh2_id;
  int gsbh1_id;
  int gsbh2_id;
  int gsbh3_id;
  int gvbh1_id;
  int gvbh2_id;
  int gzdis_id;
  int lat_id;
  int lon_id;
  int time_id;
  
    /* rank (number of dimensions) for each variable */
#   define RANK_cell_area 2
#   define RANK_gcrop 3
#   define RANK_gothr 3
#   define RANK_gurbn 3
#   define RANK_gpast 3
#   define RANK_gsecd 3
#   define RANK_gssma 3
#   define RANK_gssmb 3
#   define RANK_gsumm 3
#   define RANK_gflcp 3
#   define RANK_gflcs 3
#   define RANK_gflpc 3
#   define RANK_gflps 3
#   define RANK_gflsc 3
#   define RANK_gflsp 3
#   define RANK_gflvc 3
#   define RANK_gflvp 3
#   define RANK_gfsh1 3
#   define RANK_gfsh2 3
#   define RANK_gfsh3 3
#   define RANK_gfvh1 3
#   define RANK_gfvh2 3
#   define RANK_gsbh1 3
#   define RANK_gsbh2 3
#   define RANK_gsbh3 3
#   define RANK_gvbh1 3
#   define RANK_gvbh2 3
#   define RANK_gzdis 3
#   define RANK_lat 1
#   define RANK_lon 1
#   define RANK_time 1

    /* variable shapes */
    int cell_area_dims[RANK_cell_area];
    int gcrop_dims[RANK_gcrop];
    int gothr_dims[RANK_gothr];
    int gurbn_dims[RANK_gurbn];
    int gpast_dims[RANK_gpast];
    int gsecd_dims[RANK_gsecd];
    int gssma_dims[RANK_gssma];
    int gssmb_dims[RANK_gssmb];
    int gsumm_dims[RANK_gsumm];
    int gflcp_dims[RANK_gflcp];
    int gflcs_dims[RANK_gflcs];
    int gflpc_dims[RANK_gflpc];
    int gflps_dims[RANK_gflps];
    int gflsc_dims[RANK_gflsc];
    int gflsp_dims[RANK_gflsp];
    int gflvc_dims[RANK_gflvc];
    int gflvp_dims[RANK_gflvp];
    int gfsh1_dims[RANK_gfsh1];
    int gfsh2_dims[RANK_gfsh2];
    int gfsh3_dims[RANK_gfsh3];
    int gfvh1_dims[RANK_gfvh1];
    int gfvh2_dims[RANK_gfvh2];
    int gsbh1_dims[RANK_gsbh1];
    int gsbh2_dims[RANK_gsbh2];
    int gsbh3_dims[RANK_gsbh3];
    int gvbh1_dims[RANK_gvbh1];
    int gvbh2_dims[RANK_gvbh2];
    int gzdis_dims[RANK_gzdis];
    int lat_dims[RANK_lat];
    int lon_dims[RANK_lon];
    int time_dims[RANK_time];
    char curryearc[8];

    printf("output updated state to netcdf file \n");

    // update flows needed for output to clm - these are not needed for restart
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){
	//	printf("k is %ld m is %ld",k,m);
	newdata[0].flowsbh[k][m]=(newdata[0].smb[k][m] > ZEROVALUE) ? newdata[0].sbh[k][m]/newdata[0].smb[k][m]/garea[k][m] : ZEROVALUE;
	//	if (newdata[0].smb[k][m] > ZEROVALUE) {
	//	  printf("sbh smb garea flowsbh %f %f %f %f \n", newdata[0].sbh[k][m],newdata[0].smb[k][m],garea[k][m], newdata[0].sbh[k][m]/newdata[0].smb[k][m]/garea[k][m]);
	//	    }
	newdata[0].flowvbh[k][m]=(dstatic[k][m].vba > ZEROVALUE) ? newdata[0].vbh[k][m]/dstatic[k][m].vba/garea[k][m] : ZEROVALUE;
	newdata[0].flowsbh2[k][m]=(newdata[0].smb[k][m] > ZEROVALUE) ? newdata[0].sbh2[k][m]/newdata[0].smb[k][m]/garea[k][m] : ZEROVALUE;
	newdata[0].flowvbh2[k][m]=(dstatic[k][m].vba > ZEROVALUE) ? newdata[0].vbh2[k][m]/dstatic[k][m].vba/garea[k][m] : ZEROVALUE;
	newdata[0].flowsbh3[k][m]=(newdata[0].smb[k][m] > ZEROVALUE) ? newdata[0].flowsbh3[k][m]=newdata[0].sbh3[k][m]/newdata[0].smb[k][m]/garea[k][m] : ZEROVALUE;
      }
    }

    sprintf(curryearc, "%d", curryear);
    //jt fix the casename
    strcpy(state_file,"test.glm.state.");
    strcat(state_file,strcat(curryearc,".nc"));
    /* enter define mode */
    stat = nc_create(state_file, NC_CLOBBER, &ncid_state);
    check_err(stat,__LINE__,__FILE__);

    /* define dimensions */
    stat = nc_def_dim(ncid_state, "lat", lat_len, &lat_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid_state, "lon", lon_len, &lon_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid_state, "time", time_len, &time_dim);
    check_err(stat,__LINE__,__FILE__);

    /* define variables */

    cell_area_dims[0] = lat_dim;
    cell_area_dims[1] = lon_dim;
    stat = nc_def_var(ncid_state, "cell_area", NC_FLOAT, RANK_cell_area, cell_area_dims, &cell_area_id);
    check_err(stat,__LINE__,__FILE__);

    gcrop_dims[0] = time_dim;
    gcrop_dims[1] = lat_dim;
    gcrop_dims[2] = lon_dim;
    stat = nc_def_var(ncid_state, "gcrop", NC_FLOAT, RANK_gcrop, gcrop_dims, &gcrop_id);
    check_err(stat,__LINE__,__FILE__);

    gothr_dims[0] = time_dim;
    gothr_dims[1] = lat_dim;
    gothr_dims[2] = lon_dim;
    stat = nc_def_var(ncid_state, "gothr", NC_FLOAT, RANK_gothr, gothr_dims, &gothr_id);
    check_err(stat,__LINE__,__FILE__);

    gurbn_dims[0] = time_dim;
    gurbn_dims[1] = lat_dim;
    gurbn_dims[2] = lon_dim;
    stat = nc_def_var(ncid_state, "gurbn", NC_FLOAT, RANK_gurbn, gurbn_dims, &gurbn_id);
    check_err(stat,__LINE__,__FILE__);

    gpast_dims[0] = time_dim;
    gpast_dims[1] = lat_dim;
    gpast_dims[2] = lon_dim;
    stat = nc_def_var(ncid_state, "gpast", NC_FLOAT, RANK_gpast, gpast_dims, &gpast_id);
    check_err(stat,__LINE__,__FILE__);

    gsecd_dims[0] = time_dim;
    gsecd_dims[1] = lat_dim;
    gsecd_dims[2] = lon_dim;
    stat = nc_def_var(ncid_state, "gsecd", NC_FLOAT, RANK_gsecd, gsecd_dims, &gsecd_id);
    check_err(stat,__LINE__,__FILE__);

    gssma_dims[0] = time_dim;
    gssma_dims[1] = lat_dim;
    gssma_dims[2] = lon_dim;
    stat = nc_def_var(ncid_state, "gssma", NC_FLOAT, RANK_gssma, gssma_dims, &gssma_id);
    check_err(stat,__LINE__,__FILE__);

    gssmb_dims[0] = time_dim;
    gssmb_dims[1] = lat_dim;
    gssmb_dims[2] = lon_dim;
    stat = nc_def_var(ncid_state, "gssmb", NC_FLOAT, RANK_gssmb, gssmb_dims, &gssmb_id);
    check_err(stat,__LINE__,__FILE__);

    gsumm_dims[0] = time_dim;
    gsumm_dims[1] = lat_dim;
    gsumm_dims[2] = lon_dim;
    stat = nc_def_var(ncid_state, "gsumm", NC_FLOAT, RANK_gsumm, gsumm_dims, &gsumm_id);
    check_err(stat,__LINE__,__FILE__);

/*     gflcp_dims[0] = time_dim; */
/*     gflcp_dims[1] = lat_dim; */
/*     gflcp_dims[2] = lon_dim; */
/*     stat = nc_def_var(ncid_state, "gflcp", NC_FLOAT, RANK_gflcp, gflcp_dims, &gflcp_id); */
/*     check_err(stat,__LINE__,__FILE__); */

/*     gflcs_dims[0] = time_dim; */
/*     gflcs_dims[1] = lat_dim; */
/*     gflcs_dims[2] = lon_dim; */
/*     stat = nc_def_var(ncid_state, "gflcs", NC_FLOAT, RANK_gflcs, gflcs_dims, &gflcs_id); */
/*     check_err(stat,__LINE__,__FILE__); */

/*     gflpc_dims[0] = time_dim; */
/*     gflpc_dims[1] = lat_dim; */
/*     gflpc_dims[2] = lon_dim; */
/*     stat = nc_def_var(ncid_state, "gflpc", NC_FLOAT, RANK_gflpc, gflpc_dims, &gflpc_id); */
/*     check_err(stat,__LINE__,__FILE__); */

/*     gflps_dims[0] = time_dim; */
/*     gflps_dims[1] = lat_dim; */
/*     gflps_dims[2] = lon_dim; */
/*     stat = nc_def_var(ncid_state, "gflps", NC_FLOAT, RANK_gflps, gflps_dims, &gflps_id); */
/*     check_err(stat,__LINE__,__FILE__); */

/*     gflsc_dims[0] = time_dim; */
/*     gflsc_dims[1] = lat_dim; */
/*     gflsc_dims[2] = lon_dim; */
/*     stat = nc_def_var(ncid_state, "gflsc", NC_FLOAT, RANK_gflsc, gflsc_dims, &gflsc_id); */
/*     check_err(stat,__LINE__,__FILE__); */

/*     gflsp_dims[0] = time_dim; */
/*     gflsp_dims[1] = lat_dim; */
/*     gflsp_dims[2] = lon_dim; */
/*     stat = nc_def_var(ncid_state, "gflsp", NC_FLOAT, RANK_gflsp, gflsp_dims, &gflsp_id); */
/*     check_err(stat,__LINE__,__FILE__); */

/*     gflvc_dims[0] = time_dim; */
/*     gflvc_dims[1] = lat_dim; */
/*     gflvc_dims[2] = lon_dim; */
/*     stat = nc_def_var(ncid_state, "gflvc", NC_FLOAT, RANK_gflvc, gflvc_dims, &gflvc_id); */
/*     check_err(stat,__LINE__,__FILE__); */

/*     gflvp_dims[0] = time_dim; */
/*     gflvp_dims[1] = lat_dim; */
/*     gflvp_dims[2] = lon_dim; */
/*     stat = nc_def_var(ncid_state, "gflvp", NC_FLOAT, RANK_gflvp, gflvp_dims, &gflvp_id); */
/*     check_err(stat,__LINE__,__FILE__); */

    gfsh1_dims[0] = time_dim;
    gfsh1_dims[1] = lat_dim;
    gfsh1_dims[2] = lon_dim;
    stat = nc_def_var(ncid_state, "gfsh1", NC_FLOAT, RANK_gfsh1, gfsh1_dims, &gfsh1_id);
    check_err(stat,__LINE__,__FILE__);

    gfsh2_dims[0] = time_dim;
    gfsh2_dims[1] = lat_dim;
    gfsh2_dims[2] = lon_dim;
    stat = nc_def_var(ncid_state, "gfsh2", NC_FLOAT, RANK_gfsh2, gfsh2_dims, &gfsh2_id);
    check_err(stat,__LINE__,__FILE__);

    gfsh3_dims[0] = time_dim;
    gfsh3_dims[1] = lat_dim;
    gfsh3_dims[2] = lon_dim;
    stat = nc_def_var(ncid_state, "gfsh3", NC_FLOAT, RANK_gfsh3, gfsh3_dims, &gfsh3_id);
    check_err(stat,__LINE__,__FILE__);

    gfvh1_dims[0] = time_dim;
    gfvh1_dims[1] = lat_dim;
    gfvh1_dims[2] = lon_dim;
    stat = nc_def_var(ncid_state, "gfvh1", NC_FLOAT, RANK_gfvh1, gfvh1_dims, &gfvh1_id);
    check_err(stat,__LINE__,__FILE__);

    gfvh2_dims[0] = time_dim;
    gfvh2_dims[1] = lat_dim;
    gfvh2_dims[2] = lon_dim;
    stat = nc_def_var(ncid_state, "gfvh2", NC_FLOAT, RANK_gfvh2, gfvh2_dims, &gfvh2_id);
    check_err(stat,__LINE__,__FILE__);

/*     gsbh1_dims[0] = time_dim; */
/*     gsbh1_dims[1] = lat_dim; */
/*     gsbh1_dims[2] = lon_dim; */
/*     stat = nc_def_var(ncid_state, "gsbh1", NC_FLOAT, RANK_gsbh1, gsbh1_dims, &gsbh1_id); */
/*     check_err(stat,__LINE__,__FILE__); */

/*     gsbh2_dims[0] = time_dim; */
/*     gsbh2_dims[1] = lat_dim; */
/*     gsbh2_dims[2] = lon_dim; */
/*     stat = nc_def_var(ncid_state, "gsbh2", NC_FLOAT, RANK_gsbh2, gsbh2_dims, &gsbh2_id); */
/*     check_err(stat,__LINE__,__FILE__); */

/*     gsbh3_dims[0] = time_dim; */
/*     gsbh3_dims[1] = lat_dim; */
/*     gsbh3_dims[2] = lon_dim; */
/*     stat = nc_def_var(ncid_state, "gsbh3", NC_FLOAT, RANK_gsbh3, gsbh3_dims, &gsbh3_id); */
/*     check_err(stat,__LINE__,__FILE__); */

/*     gvbh1_dims[0] = time_dim; */
/*     gvbh1_dims[1] = lat_dim; */
/*     gvbh1_dims[2] = lon_dim; */
/*     stat = nc_def_var(ncid_state, "gvbh1", NC_FLOAT, RANK_gvbh1, gvbh1_dims, &gvbh1_id); */
/*     check_err(stat,__LINE__,__FILE__); */

/*     gvbh2_dims[0] = time_dim; */
/*     gvbh2_dims[1] = lat_dim; */
/*     gvbh2_dims[2] = lon_dim; */
/*     stat = nc_def_var(ncid_state, "gvbh2", NC_FLOAT, RANK_gvbh2, gvbh2_dims, &gvbh2_id); */
/*     check_err(stat,__LINE__,__FILE__); */

/*     gzdis_dims[0] = time_dim; */
/*     gzdis_dims[1] = lat_dim; */
/*     gzdis_dims[2] = lon_dim; */
/*     stat = nc_def_var(ncid_state, "gzdis", NC_FLOAT, RANK_gzdis, gzdis_dims, &gzdis_id); */
/*     check_err(stat,__LINE__,__FILE__); */

    lat_dims[0] = lat_dim;
    stat = nc_def_var(ncid_state, "lat", NC_DOUBLE, RANK_lat, lat_dims, &lat_id);
    check_err(stat,__LINE__,__FILE__);

    lon_dims[0] = lon_dim;
    stat = nc_def_var(ncid_state, "lon", NC_DOUBLE, RANK_lon, lon_dims, &lon_id);
    check_err(stat,__LINE__,__FILE__);

    time_dims[0] = time_dim;
    stat = nc_def_var(ncid_state, "time", NC_DOUBLE, RANK_time, time_dims, &time_id);
    check_err(stat,__LINE__,__FILE__);

    /* assign global attributes */
    { /* CDI */
    stat = nc_put_att_text(ncid_state, NC_GLOBAL, "CDI", 36, "Climate Data Interface version 1.4.1");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* Conventions */
    stat = nc_put_att_text(ncid_state, NC_GLOBAL, "Conventions", 6, "CF-1.0");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* history */
    stat = nc_put_att_text(ncid_state, NC_GLOBAL, "history", 1, " ");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* CDO */
    stat = nc_put_att_text(ncid_state, NC_GLOBAL, "CDO", 67, "Climate Data Operators version 1.4.1 (http://www.mpimet.mpg.de/cdo)");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* NCO */
    stat = nc_put_att_text(ncid_state, NC_GLOBAL, "NCO", 5, "4.0.5");
    check_err(stat,__LINE__,__FILE__);
    }

    /* assign per-variable attributes */
    { /* long_name */
    stat = nc_put_att_text(ncid_state, cell_area_id, "long_name", 17, "area of grid cell");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* standard_name */
    stat = nc_put_att_text(ncid_state, cell_area_id, "standard_name", 4, "area");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid_state, cell_area_id, "units", 3, "m^2");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid_state, gcrop_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid_state, gcrop_id, "long_name", 20, "fraction of cropland");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid_state, gothr_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid_state, gothr_id, "long_name", 24, "fraction of primary land");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid_state, gurbn_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid_state, gurbn_id, "long_name", 24, "fraction of urban land");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid_state, gpast_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid_state, gpast_id, "long_name", 19, "fraction of pasture");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid_state, gsecd_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid_state, gsecd_id, "long_name", 26, "fraction of secondary_land");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid_state, gssma_id, "units", 4, "year");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid_state, gssma_id, "long_name", 38, "mean age of secondary land in gridcell");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid_state, gssmb_id, "units", 7, "kgC/m^2");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid_state, gssmb_id, "long_name", 50, "mean biomass density of secondary land in gridcell");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid_state, gsumm_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid_state, gsumm_id, "long_name", 63, "sum of crop,pasture,primary,secondary,urban,ice,water fractions");
    check_err(stat,__LINE__,__FILE__);
    }
/*     { /\* units *\/ */
/*     stat = nc_put_att_text(ncid_state, gflcp_id, "units", 8, "fraction"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* long_name *\/ */
/*     stat = nc_put_att_text(ncid_state, gflcp_id, "long_name", 44, "fraction of cropland transitioned to pasture"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* units *\/ */
/*     stat = nc_put_att_text(ncid_state, gflcs_id, "units", 8, "fraction"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* long_name *\/ */
/*     stat = nc_put_att_text(ncid_state, gflcs_id, "long_name", 51, "fraction of cropland transitioned to secondary land"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* units *\/ */
/*     stat = nc_put_att_text(ncid_state, gflpc_id, "units", 8, "fraction"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* long_name *\/ */
/*     stat = nc_put_att_text(ncid_state, gflpc_id, "long_name", 44, "fraction of pasture transitioned to cropland"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* units *\/ */
/*     stat = nc_put_att_text(ncid_state, gflps_id, "units", 8, "fraction"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* long_name *\/ */
/*     stat = nc_put_att_text(ncid_state, gflps_id, "long_name", 50, "fraction of pasture transitioned to secondary land"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* units *\/ */
/*     stat = nc_put_att_text(ncid_state, gflsc_id, "units", 8, "fraction"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* long_name *\/ */
/*     stat = nc_put_att_text(ncid_state, gflsc_id, "long_name", 46, "fraction of secondary transitioned to cropland"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* units *\/ */
/*     stat = nc_put_att_text(ncid_state, gflsp_id, "units", 8, "fraction"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* long_name *\/ */
/*     stat = nc_put_att_text(ncid_state, gflsp_id, "long_name", 45, "fraction of secondary transitioned to pasture"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* units *\/ */
/*     stat = nc_put_att_text(ncid_state, gflvc_id, "units", 8, "fraction"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* long_name *\/ */
/*     stat = nc_put_att_text(ncid_state, gflvc_id, "long_name", 44, "fraction of primary transitioned to cropland"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* units *\/ */
/*     stat = nc_put_att_text(ncid_state, gflvp_id, "units", 8, "fraction"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* long_name *\/ */
/*     stat = nc_put_att_text(ncid_state, gflvp_id, "long_name", 43, "fraction of primary transitioned to pasture"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
    { /* units */
    stat = nc_put_att_text(ncid_state, gfsh1_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid_state, gfsh1_id, "long_name", 78, "gridcell fraction that had wood harvested from mature secondary forested land ");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid_state, gfsh2_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid_state, gfsh2_id, "long_name", 77, "gridcell fraction that had wood harvested from young secondary forested land ");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid_state, gfsh3_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid_state, gfsh3_id, "long_name", 80, "gridcell fraction that had wood harvested from young secondary nonforested land ");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid_state, gfvh1_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid_state, gfvh1_id, "long_name", 69, "gridcell fraction that had wood harvested from primary forested land ");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid_state, gfvh2_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid_state, gfvh2_id, "long_name", 72, "gridcell fraction that had wood harvested from primary nonforested land ");
    check_err(stat,__LINE__,__FILE__);
    }
/*     { /\* units *\/ */
/*     stat = nc_put_att_text(ncid_state, gsbh1_id, "units", 3, "kgC"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* long_name *\/ */
/*     stat = nc_put_att_text(ncid_state, gsbh1_id, "long_name", 60, "mature secondary forest biomass harvested from each gridcell"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* units *\/ */
/*     stat = nc_put_att_text(ncid_state, gsbh2_id, "units", 3, "kgC"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* long_name *\/ */
/*     stat = nc_put_att_text(ncid_state, gsbh2_id, "long_name", 59, "young secondary forest biomass harvested from each gridcell"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* units *\/ */
/*     stat = nc_put_att_text(ncid_state, gsbh3_id, "units", 3, "kgC"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* long_name *\/ */
/*     stat = nc_put_att_text(ncid_state, gsbh3_id, "long_name", 58, "secondary non-forest biomass harvested from each gridcell "); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* units *\/ */
/*     stat = nc_put_att_text(ncid_state, gvbh1_id, "units", 3, "kgC"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* long_name *\/ */
/*     stat = nc_put_att_text(ncid_state, gvbh1_id, "long_name", 52, "primary forest biomass harvested from each gridcell "); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* units *\/ */
/*     stat = nc_put_att_text(ncid_state, gvbh2_id, "units", 3, "kgC"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* long_name *\/ */
/*     stat = nc_put_att_text(ncid_state, gvbh2_id, "long_name", 57, "primary non-forest biomass harvested from each gridcell  "); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* units *\/ */
/*     stat = nc_put_att_text(ncid_state, gzdis_id, "units", 8, "fraction"); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
/*     { /\* long_name *\/ */
/*     stat = nc_put_att_text(ncid_state, gzdis_id, "long_name", 1, " "); */
/*     check_err(stat,__LINE__,__FILE__); */
/*     } */
    { /* long_name */
    stat = nc_put_att_text(ncid_state, lat_id, "long_name", 8, "latitude");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid_state, lat_id, "units", 13, "degrees_north");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* standard_name */
    stat = nc_put_att_text(ncid_state, lat_id, "standard_name", 8, "latitude");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid_state, lon_id, "long_name", 9, "longitude");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid_state, lon_id, "units", 12, "degrees_east");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* standard_name */
    stat = nc_put_att_text(ncid_state, lon_id, "standard_name", 9, "longitude");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid_state, time_id, "units", 13, "year as %Y.%f");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* calendar */
    stat = nc_put_att_text(ncid_state, time_id, "calendar", 19, "proleptic_gregorian");
    check_err(stat,__LINE__,__FILE__);
    }


    /* leave define mode */
    stat = nc_enddef (ncid_state);
    check_err(stat,__LINE__,__FILE__);

    /* at end of run put out time 1 of array */
    timeidx=0;
    if (curryear == stop_year) timeidx=1;
    
    /* assign variable data */
    stat = nc_put_vara_double(ncid_state, gothr_id, start, count, &newdata[timeidx].v[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid_state, gurbn_id, start, count, &newdata[timeidx].u[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid_state, gcrop_id, start, count, &newdata[timeidx].c[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid_state, gpast_id, start, count, &newdata[timeidx].p[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid_state, gsecd_id, start, count, &newdata[timeidx].s[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid_state, gssmb_id, start, count, &newdata[timeidx].smb[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid_state, gssma_id, start, count, &newdata[timeidx].sma[0][0]);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_put_vara_double(ncid_state, gfvh1_id, start, count, &newdata[timeidx].flowvbh[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid_state, gfvh2_id, start, count, &newdata[timeidx].flowvbh2[0][0]);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_put_vara_double(ncid_state, gfsh1_id, start, count, &newdata[timeidx].flowsbh[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid_state, gfsh2_id, start, count, &newdata[timeidx].flowsbh2[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid_state, gfsh3_id, start, count, &newdata[timeidx].flowsbh3[0][0]);
    check_err(stat,__LINE__,__FILE__);

/*     for (k=0;k<NY;k++){ */
/*       for (m=0;m<NX;m++){ */
/* 	sumdata[k][m]=newdata[0].c[k][m]+newdata[0].p[k][m]+newdata[0].s[k][m]+newdata[0].v[k][m]+newdata[0].i[k][m]+newdata[0].w[k][m];  */
/*       } */
/*     }       */
/*     stat = nc_put_vara_double(ncid_state, gsumm_id, start, count, &sumdata); */
/*     check_err(stat,__LINE__,__FILE__); */
    
    stat = nc_close(ncid_state);
    check_err(stat,__LINE__,__FILE__);

    return;
}

/********************************************************************/
void output_lu_nc(int curryear){
  int stat,k,m,ncid;
  char curryearc[8], lu_file[256];
  static long start[] = {0, 0, 0};   /* start at first value */
  static long count[] = {1, NY, NX};
  double f_sbh[NY][NX],f_sbh2[NY][NX],f_vbh[NY][NX],f_vbh2[NY][NX],f_sbh3[NY][NX],curryeard;

    /* dimension ids */
    int lat_dim;
    int lon_dim;
    int time_dim;

    /* dimension lengths */
    size_t lat_len = 360;
    size_t lon_len = 720;
    size_t time_len = NC_UNLIMITED;

    /* variable ids */
    int cell_area_id;
    int cp_id;
    int pc_id;
    int pv_id;
    int vp_id;
    int vc_id;
    int cv_id;
    int sc_id;
    int cs_id;
    int sp_id;
    int ps_id;
    int vs_id;
    int sbh_id;
    int vbh_id;
    int sbh2_id;
    int vbh2_id;
    int sbh3_id;
    int f_sbh3_id;
    int f_sbh2_id;
    int f_sbh_id;
    int f_vbh2_id;
    int f_vbh_id;
    int lat_id;
    int lon_id;
    int time_id;

    /* rank (number of dimensions) for each variable */
#   define RANK_cell_area 2
#   define RANK_cp 3
#   define RANK_pc 3
#   define RANK_pv 3
#   define RANK_vp 3
#   define RANK_vc 3
#   define RANK_cv 3
#   define RANK_sc 3
#   define RANK_cs 3
#   define RANK_sp 3
#   define RANK_ps 3
#   define RANK_vs 3
#   define RANK_sbh 3
#   define RANK_vbh 3
#   define RANK_sbh2 3
#   define RANK_vbh2 3
#   define RANK_sbh3 3
#   define RANK_f_sbh3 3
#   define RANK_f_sbh2 3
#   define RANK_f_sbh 3
#   define RANK_f_vbh2 3
#   define RANK_f_vbh 3
#   define RANK_lat 1
#   define RANK_lon 1
#   define RANK_time 1

    /* variable shapes */
    int cell_area_dims[RANK_cell_area];
    int cp_dims[RANK_cp];
    int pc_dims[RANK_pc];
    int pv_dims[RANK_pv];
    int vp_dims[RANK_vp];
    int vc_dims[RANK_vc];
    int cv_dims[RANK_cv];
    int sc_dims[RANK_sc];
    int cs_dims[RANK_cs];
    int sp_dims[RANK_sp];
    int ps_dims[RANK_ps];
    int vs_dims[RANK_vs];
    int sbh_dims[RANK_sbh];
    int vbh_dims[RANK_vbh];
    int sbh2_dims[RANK_sbh2];
    int vbh2_dims[RANK_vbh2];
    int sbh3_dims[RANK_sbh3];
    int f_sbh3_dims[RANK_f_sbh3];
    int f_sbh2_dims[RANK_f_sbh2];
    int f_sbh_dims[RANK_f_sbh];
    int f_vbh2_dims[RANK_f_vbh2];
    int f_vbh_dims[RANK_f_vbh];
    int lat_dims[RANK_lat];
    int lon_dims[RANK_lon];
    int time_dims[RANK_time];

    printf("output land change netcdf file \n");
 
    sprintf(curryearc, "%d", curryear);
    //jt fix the casename
    strcpy(lu_file,"test.glm.lu.");
    strcat(lu_file,strcat(curryearc,".nc"));

    /* enter define mode */
    stat = nc_create(lu_file, NC_CLOBBER, &ncid);
    check_err(stat,__LINE__,__FILE__);

    /* define dimensions */
    stat = nc_def_dim(ncid, "lat", lat_len, &lat_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "lon", lon_len, &lon_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "time", time_len, &time_dim);
    check_err(stat,__LINE__,__FILE__);

    /* define variables */

    cell_area_dims[0] = lat_dim;
    cell_area_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "cell_area", NC_FLOAT, RANK_cell_area, cell_area_dims, &cell_area_id);
    check_err(stat,__LINE__,__FILE__);

    cp_dims[0] = time_dim;
    cp_dims[1] = lat_dim;
    cp_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "cp", NC_FLOAT, RANK_cp, cp_dims, &cp_id);
    check_err(stat,__LINE__,__FILE__);

    pc_dims[0] = time_dim;
    pc_dims[1] = lat_dim;
    pc_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "pc", NC_FLOAT, RANK_pc, pc_dims, &pc_id);
    check_err(stat,__LINE__,__FILE__);

    pv_dims[0] = time_dim;
    pv_dims[1] = lat_dim;
    pv_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "pv", NC_FLOAT, RANK_pv, pv_dims, &pv_id);
    check_err(stat,__LINE__,__FILE__);

    vp_dims[0] = time_dim;
    vp_dims[1] = lat_dim;
    vp_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "vp", NC_FLOAT, RANK_vp, vp_dims, &vp_id);
    check_err(stat,__LINE__,__FILE__);

    vc_dims[0] = time_dim;
    vc_dims[1] = lat_dim;
    vc_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "vc", NC_FLOAT, RANK_vc, vc_dims, &vc_id);
    check_err(stat,__LINE__,__FILE__);

    cv_dims[0] = time_dim;
    cv_dims[1] = lat_dim;
    cv_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "cv", NC_FLOAT, RANK_cv, cv_dims, &cv_id);
    check_err(stat,__LINE__,__FILE__);

    sc_dims[0] = time_dim;
    sc_dims[1] = lat_dim;
    sc_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "sc", NC_FLOAT, RANK_sc, sc_dims, &sc_id);
    check_err(stat,__LINE__,__FILE__);

    cs_dims[0] = time_dim;
    cs_dims[1] = lat_dim;
    cs_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "cs", NC_FLOAT, RANK_cs, cs_dims, &cs_id);
    check_err(stat,__LINE__,__FILE__);

    sp_dims[0] = time_dim;
    sp_dims[1] = lat_dim;
    sp_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "sp", NC_FLOAT, RANK_sp, sp_dims, &sp_id);
    check_err(stat,__LINE__,__FILE__);

    ps_dims[0] = time_dim;
    ps_dims[1] = lat_dim;
    ps_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "ps", NC_FLOAT, RANK_ps, ps_dims, &ps_id);
    check_err(stat,__LINE__,__FILE__);

    vs_dims[0] = time_dim;
    vs_dims[1] = lat_dim;
    vs_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "vs", NC_FLOAT, RANK_vs, vs_dims, &vs_id);
    check_err(stat,__LINE__,__FILE__);

    sbh_dims[0] = time_dim;
    sbh_dims[1] = lat_dim;
    sbh_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "sbh", NC_FLOAT, RANK_sbh, sbh_dims, &sbh_id);
    check_err(stat,__LINE__,__FILE__);

    vbh_dims[0] = time_dim;
    vbh_dims[1] = lat_dim;
    vbh_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "vbh", NC_FLOAT, RANK_vbh, vbh_dims, &vbh_id);
    check_err(stat,__LINE__,__FILE__);

    sbh2_dims[0] = time_dim;
    sbh2_dims[1] = lat_dim;
    sbh2_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "sbh2", NC_FLOAT, RANK_sbh2, sbh2_dims, &sbh2_id);
    check_err(stat,__LINE__,__FILE__);

    vbh2_dims[0] = time_dim;
    vbh2_dims[1] = lat_dim;
    vbh2_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "vbh2", NC_FLOAT, RANK_vbh2, vbh2_dims, &vbh2_id);
    check_err(stat,__LINE__,__FILE__);

    sbh3_dims[0] = time_dim;
    sbh3_dims[1] = lat_dim;
    sbh3_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "sbh3", NC_FLOAT, RANK_sbh3, sbh3_dims, &sbh3_id);
    check_err(stat,__LINE__,__FILE__);

    f_sbh3_dims[0] = time_dim;
    f_sbh3_dims[1] = lat_dim;
    f_sbh3_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "f_sbh3", NC_FLOAT, RANK_f_sbh3, f_sbh3_dims, &f_sbh3_id);
    check_err(stat,__LINE__,__FILE__);

    f_sbh2_dims[0] = time_dim;
    f_sbh2_dims[1] = lat_dim;
    f_sbh2_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "f_sbh2", NC_FLOAT, RANK_f_sbh2, f_sbh2_dims, &f_sbh2_id);
    check_err(stat,__LINE__,__FILE__);

    f_sbh_dims[0] = time_dim;
    f_sbh_dims[1] = lat_dim;
    f_sbh_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "f_sbh", NC_FLOAT, RANK_f_sbh, f_sbh_dims, &f_sbh_id);
    check_err(stat,__LINE__,__FILE__);

    f_vbh2_dims[0] = time_dim;
    f_vbh2_dims[1] = lat_dim;
    f_vbh2_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "f_vbh2", NC_FLOAT, RANK_f_vbh2, f_vbh2_dims, &f_vbh2_id);
    check_err(stat,__LINE__,__FILE__);

    f_vbh_dims[0] = time_dim;
    f_vbh_dims[1] = lat_dim;
    f_vbh_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "f_vbh", NC_FLOAT, RANK_f_vbh, f_vbh_dims, &f_vbh_id);
    check_err(stat,__LINE__,__FILE__);

    lat_dims[0] = lat_dim;
    stat = nc_def_var(ncid, "lat", NC_FLOAT, RANK_lat, lat_dims, &lat_id);
    check_err(stat,__LINE__,__FILE__);

    lon_dims[0] = lon_dim;
    stat = nc_def_var(ncid, "lon", NC_FLOAT, RANK_lon, lon_dims, &lon_id);
    check_err(stat,__LINE__,__FILE__);

    time_dims[0] = time_dim;
    stat = nc_def_var(ncid, "time", NC_FLOAT, RANK_time, time_dims, &time_id);
    check_err(stat,__LINE__,__FILE__);

    /* assign global attributes */
    { /* CDI */
    stat = nc_put_att_text(ncid, NC_GLOBAL, "CDI", 36, "Climate Data Interface version 1.4.1");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* Conventions */
    stat = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "CF-1.0");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* history */
    stat = nc_put_att_text(ncid, NC_GLOBAL, "history", 1, " ");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* CDO */
    stat = nc_put_att_text(ncid, NC_GLOBAL, "CDO", 67, "Climate Data Operators version 1.4.1 (http://www.mpimet.mpg.de/cdo)");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* NCO */
    stat = nc_put_att_text(ncid, NC_GLOBAL, "NCO", 5, "4.0.5");
    check_err(stat,__LINE__,__FILE__);
    }


    /* assign per-variable attributes */
    { /* long_name */
    stat = nc_put_att_text(ncid, cell_area_id, "long_name", 17, "area of grid cell");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* standard_name */
    stat = nc_put_att_text(ncid, cell_area_id, "standard_name", 4, "area");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, cell_area_id, "units", 3, "m^2");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, cp_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, cp_id, "long_name", 31, "transition from crop to pasture");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, pc_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, pc_id, "long_name", 31, "transition from pasture to crop");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, pv_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, pv_id, "long_name", 39, "transition from pasture to primary land");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, vp_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, vp_id, "long_name", 39, "transition from primary land to pasture");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, vc_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, vc_id, "long_name", 36, "transition from primary land to crop");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, cv_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, cv_id, "long_name", 36, "transition from crop to primary land");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, sc_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, sc_id, "long_name", 38, "transition from secondary land to crop");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, cs_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, cs_id, "long_name", 38, "transition from crop to secondary land");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, sp_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, sp_id, "long_name", 41, "transition from secondary land to pasture");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, ps_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, ps_id, "long_name", 46, "transition from pasture land to secondary land");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, vs_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, vs_id, "long_name", 46, "transition from primary land to secondary land");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, sbh_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, sbh_id, "long_name", 7, "unknown");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, vbh_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, vbh_id, "long_name", 7, "unknown");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, sbh2_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, sbh2_id, "long_name", 7, "unknown");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, vbh2_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, vbh2_id, "long_name", 7, "unknown");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, sbh3_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, sbh3_id, "long_name", 7, "unknown");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, f_sbh3_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, f_sbh3_id, "long_name", 7, "unknown");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, f_sbh2_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, f_sbh2_id, "long_name", 7, "unknown");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, f_sbh_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, f_sbh_id, "long_name", 7, "unknown");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, f_vbh2_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, f_vbh2_id, "long_name", 7, "unknown");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, f_vbh_id, "units", 8, "fraction");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, f_vbh_id, "long_name", 7, "unknown");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, lat_id, "long_name", 8, "latitude");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, lat_id, "units", 13, "degrees_north");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* standard_name */
    stat = nc_put_att_text(ncid, lat_id, "standard_name", 8, "latitude");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* long_name */
    stat = nc_put_att_text(ncid, lon_id, "long_name", 9, "longitude");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, lon_id, "units", 12, "degrees_east");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* standard_name */
    stat = nc_put_att_text(ncid, lon_id, "standard_name", 9, "longitude");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* units */
    stat = nc_put_att_text(ncid, time_id, "units", 13, "year as %Y.%f");
    check_err(stat,__LINE__,__FILE__);
    }
    { /* calendar */
    stat = nc_put_att_text(ncid, time_id, "calendar", 19, "proleptic_gregorian");
    check_err(stat,__LINE__,__FILE__);
    }


    /* leave define mode */
    stat = nc_enddef (ncid);
    check_err(stat,__LINE__,__FILE__);

    /* assign variable data */
    //jt  this is output only if(dstatic[k][m].gcode > 0)  need to
    //jt  check if rest of values should be set to missing.


    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){

        f_sbh[k][m]=f_sbh2[k][m]=f_vbh[k][m]=f_vbh2[k][m]=f_sbh3[k][m]=ZEROVALUE;

        if(newdata[0].flowcp[k][m] <= -ZEROVALUE) {newdata[0].flowcp[k][m] = ZEROVALUE;}
        if(newdata[0].flowpc[k][m] <= -ZEROVALUE) {newdata[0].flowpc[k][m] = ZEROVALUE;}
        if(newdata[0].flowpv[k][m] <= -ZEROVALUE) {newdata[0].flowpv[k][m] = ZEROVALUE;}        
        if(newdata[0].flowvp[k][m] <= -ZEROVALUE) {newdata[0].flowvp[k][m] = ZEROVALUE;}        
        if(newdata[0].flowvc[k][m] <= -ZEROVALUE) {newdata[0].flowvc[k][m] = ZEROVALUE;}        
        if(newdata[0].flowcv[k][m] <= -ZEROVALUE) {newdata[0].flowcv[k][m] = ZEROVALUE;}        
        if(newdata[0].flowsp[k][m] <= -ZEROVALUE) {newdata[0].flowsp[k][m] = ZEROVALUE;}        
        if(newdata[0].flowps[k][m] <= -ZEROVALUE) {newdata[0].flowps[k][m] = ZEROVALUE;}        
        if(newdata[0].flowsc[k][m] <= -ZEROVALUE) {newdata[0].flowsc[k][m] = ZEROVALUE;}        
        if(newdata[0].flowcs[k][m] <= -ZEROVALUE) {newdata[0].flowcs[k][m] = ZEROVALUE;}
        if(newdata[0].flowvs[k][m] <= -ZEROVALUE) {newdata[0].flowvs[k][m] = ZEROVALUE;}


	if(newdata[0].smb[k][m] > ZEROVALUE) f_sbh[k][m]=newdata[0].sbh[k][m]/newdata[0].smb[k][m]/garea[k][m];
	if(dstatic[k][m].vba > ZEROVALUE) f_vbh[k][m]=newdata[0].vbh[k][m]/dstatic[k][m].vba/garea[k][m];
	if(newdata[0].smb[k][m] > ZEROVALUE) f_sbh2[k][m]=newdata[0].sbh2[k][m]/newdata[0].smb[k][m]/garea[k][m];
	if(dstatic[k][m].vba > ZEROVALUE) f_vbh2[k][m]=newdata[0].vbh2[k][m]/dstatic[k][m].vba/garea[k][m];
	if(newdata[0].smb[k][m] > ZEROVALUE) f_sbh3[k][m]=newdata[0].sbh3[k][m]/newdata[0].smb[k][m]/garea[k][m];
    
    
	newdata_scale.sbh[k][m]=newdata[0].sbh[k][m]/1000;
	newdata_scale.vbh[k][m]=newdata[0].vbh[k][m]/1000;
	newdata_scale.sbh2[k][m]=newdata[0].sbh2[k][m]/1000;
	newdata_scale.vbh2[k][m]=newdata[0].vbh2[k][m]/1000;
	newdata_scale.sbh3[k][m]=newdata[0].sbh3[k][m]/1000;

      }
    }
    curryeard=curryear;
    stat = nc_put_var_double(ncid, time_id, &curryeard);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, cp_id, start, count, &newdata[0].flowcp[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, pc_id, start, count, &newdata[0].flowpc[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, pv_id, start, count, &newdata[0].flowpv[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, vp_id, start, count, &newdata[0].flowvp[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, vc_id, start, count, &newdata[0].flowvc[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, cv_id, start, count, &newdata[0].flowcv[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, sc_id, start, count, &newdata[0].flowsc[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, cs_id, start, count, &newdata[0].flowcs[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, sp_id, start, count, &newdata[0].flowsp[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, ps_id, start, count, &newdata[0].flowps[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, vs_id, start, count, &newdata[0].flowvs[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, sbh_id, start, count, &newdata_scale.sbh[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, f_sbh_id, start, count, &f_sbh[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, vbh_id, start, count, &newdata_scale.vbh[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, f_vbh_id, start, count, &f_vbh[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, sbh2_id, start, count, &newdata_scale.sbh2[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, f_sbh2_id, start, count, &f_sbh2[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, vbh2_id, start, count, &newdata_scale.vbh2[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, f_vbh2_id, start, count, &f_vbh2[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, sbh3_id, start, count, &newdata_scale.sbh3[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_vara_double(ncid, f_sbh3_id, start, count, &f_sbh3[0][0]);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_close(ncid);
    check_err(stat,__LINE__,__FILE__);

    return;
}

/********************************************************************/
void output_updated_states(int curryear){
  
  FILE *outfile; 
  int ih, ic, k, m, it, cnum, nit;
  char newpath[190], newpath1[190], fn[30], ayear[10], fout[190];
  
  
  strcpy(newpath,PATH1);
  strcat(newpath,finput_dir);
  strcat(newpath,"updated_states/");


  cnum=7;
  if (use_urban) cnum=8;
  
  
  for (ic=0;ic<cnum;ic++){

     for (it=0;it<NEWNT;it++){ 


      nit=it;
      strcpy(newpath1,newpath);
      sprintf(ayear,"%d.txt",curryear);

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

      for (ih=0;ih<HEADLEN;ih++) {
	fprintf(outfile,"%s",header[ih]);
	

      }
      
 
      for (k=0;k<NY;k++){
	for (m=0;m<NX;m++){
	  

	    if(ic == 0) fprintf(outfile,"%lf ",newdata[nit].v[k][m]);
	    if(ic == 1) fprintf(outfile,"%lf ",newdata[nit].c[k][m]);
	    if(ic == 2) fprintf(outfile,"%lf ",newdata[nit].p[k][m]);	  
	    if(ic == 3) fprintf(outfile,"%lf ",newdata[nit].s[k][m]);
	    if(ic == 4) fprintf(outfile,"%lf ",newdata[nit].smb[k][m]);
	    if(ic == 5) fprintf(outfile,"%lf ",newdata[nit].sma[k][m]);
	    if (use_urban) {
	      if(ic == 6) fprintf(outfile,"%lf ",newdata[nit].c[k][m]+newdata[nit].p[k][m]+newdata[nit].s[k][m]+newdata[nit].v[k][m]+newdata[nit].i[k][m]+newdata[nit].w[k][m]+newdata[nit].u[k][m]); 
	      if(ic == 7) fprintf(outfile,"%lf ",newdata[nit].u[k][m]);
	    }else{
	    if(ic == 6) fprintf(outfile,"%lf ",newdata[nit].c[k][m]+newdata[nit].p[k][m]+newdata[nit].s[k][m]+newdata[nit].v[k][m]+newdata[nit].i[k][m]+newdata[nit].w[k][m]); 
	    }
	}
        fprintf(outfile,"\n");

      }
      fclose(outfile);
    }
  }



  return;
}

/********************************************************************/
void output_updated_states2(int curryear){
  
  FILE *outfile; 
  int ih, ic, k, m, it, cnum, nit;
  char newpath[190], newpath1[190], fn[30], ayear[10], fout[190];
  double  flowcp,flowpc,flowpv,flowvp,flowvc,flowcv,flowsc,flowcs,flowsp,flowps;
  double flowsbh,flowvbh,flowsbh2,flowvbh2,flowsbh3;

  
  
  strcpy(newpath,PATH1);
  strcat(newpath,finput_dir);
  strcat(newpath,"updated_states/");

  cnum=6;
  
  for (ic=0;ic<cnum;ic++){
    
     for (it=0;it<NEWNT;it++){ 

      nit=it;

      strcpy(newpath1,newpath);
      sprintf(ayear,"%d.txt",curryear);


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
	  

	  flowcp=newdata[nit].flowcp[k][m]*garea[k][m]/1000.0/1000.0;
	  flowpc=newdata[nit].flowpc[k][m]*garea[k][m]/1000.0/1000.0;
	  flowpv=newdata[nit].flowpv[k][m]*garea[k][m]/1000.0/1000.0;
	  flowvp=newdata[nit].flowvp[k][m]*garea[k][m]/1000.0/1000.0;
	  flowvc=newdata[nit].flowvc[k][m]*garea[k][m]/1000.0/1000.0;
	  flowcv=newdata[nit].flowcv[k][m]*garea[k][m]/1000.0/1000.0;
	  flowsc=newdata[nit].flowsc[k][m]*garea[k][m]/1000.0/1000.0;
	  flowcs=newdata[nit].flowcs[k][m]*garea[k][m]/1000.0/1000.0;
	  flowsp=newdata[nit].flowsp[k][m]*garea[k][m]/1000.0/1000.0;
	  flowps=newdata[nit].flowps[k][m]*garea[k][m]/1000.0/1000.0;
          flowsbh=ZEROVALUE;
	  if(newdata[nit].smb[k][m] > ZEROVALUE) flowsbh=newdata[nit].sbh[k][m]/newdata[nit].smb[k][m]/1000.0/1000.0;
          flowvbh=ZEROVALUE;
	  if(dstatic[k][m].vba > ZEROVALUE) flowvbh=newdata[nit].vbh[k][m]/dstatic[k][m].vba/1000.0/1000.0;
          flowsbh2=ZEROVALUE;
	  if(newdata[nit].smb[k][m] > ZEROVALUE) flowsbh2=newdata[nit].sbh2[k][m]/newdata[nit].smb[k][m]/1000.0/1000.0;
          flowvbh2=ZEROVALUE;
	  if(dstatic[k][m].vba > ZEROVALUE) flowvbh2=newdata[nit].vbh2[k][m]/dstatic[k][m].vba/1000.0/1000.0;

	  flowsbh3=ZEROVALUE;
	  if(newdata[nit].smb[k][m] > ZEROVALUE) flowsbh3=newdata[nit].sbh3[k][m]/newdata[nit].smb[k][m]/garea[k][m];

	  if(ic == 0) fprintf(outfile,"%lf ",newdata[nit].sbh[k][m]);
	  if(ic == 1) fprintf(outfile,"%lf ",newdata[nit].vbh[k][m]);
	  if(ic == 2) fprintf(outfile,"%lf ",newdata[nit].sbh2[k][m]);
	  if(ic == 3) fprintf(outfile,"%lf ",newdata[nit].vbh2[k][m]);
	  if(ic == 4) fprintf(outfile,"%lf ",newdata[nit].sbh3[k][m]);
	  if(ic == 5) fprintf(outfile,"%d ",newdata[nit].zdis[k][m]);           
	  if(ic == 6){
	    if(dstatic[k][m].fnf == 1){
	      fprintf(outfile,"%lf ",newdata[nit].v[k][m]);
	    }
	    else{
	      fprintf(outfile,"0.00000000 ");
	    }
	  }
	  if(ic == 7){
	    if(dstatic[k][m].fnf == 1){
	      fprintf(outfile,"%lf ",newdata[nit].s[k][m]);
	    }
	    else{
	      fprintf(outfile,"0.00000000 ");
	    }
	  }
	  if(ic == 8){
	    if(dstatic[k][m].fnf == 0){
	      fprintf(outfile,"%lf ",newdata[nit].v[k][m]);
	    }
	    else{
	      fprintf(outfile,"0.00000000 ");
	    }
	  }
	  if(ic == 9){
	    if(dstatic[k][m].fnf == 0){
	      fprintf(outfile,"%lf ",newdata[nit].s[k][m]);
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
void output_updated_states3(int curryear){
  
  FILE *outfile; 
  int ih, ic, k, m, it, cnum, nit;
  char newpath[190], newpath1[190], fn[30], ayear[10], fout[190];
  double  flowcp,flowpc,flowpv,flowvp,flowvc,flowcv,flowsc,flowcs,flowsp,flowps;
  double flowsbh,flowvbh,flowsbh2,flowvbh2,flowsbh3;
  double flowcu,flowpu,flowvu,flowsu,flowuc,flowup,flowus;

  
  
  strcpy(newpath,PATH1);
  strcat(newpath,finput_dir);
  strcat(newpath,"updated_states/");


  cnum=13;
  if (use_urban) cnum=20;
  
  
  for (ic=0;ic<cnum;ic++){

     for (it=0;it<NEWNT;it++){ 

      nit=it;


      strcpy(newpath1,newpath);
      sprintf(ayear,"%d.txt",curryear);

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
      if (use_urban) {
	if(ic == 13) strcpy(fn,"gflcu.");
	if(ic == 14) strcpy(fn,"gflpu.");
	if(ic == 15) strcpy(fn,"gflvu.");
	if(ic == 16) strcpy(fn,"gflsu.");
	if(ic == 17) strcpy(fn,"gfluc.");
	if(ic == 18) strcpy(fn,"gflup.");
	if(ic == 19) strcpy(fn,"gflus.");
      }
      strcat(fn,ayear);
      strcat(newpath1,fn);
      strcpy(fout,newpath1); 

      outfile=fopen(fout,"w");

      for (ih=0;ih<HEADLEN;ih++) fprintf(outfile,"%s",header[ih]);
      
 
      for (k=0;k<NY;k++){
	for (m=0;m<NX;m++){
	  

	  flowcp=newdata[nit].flowcp[k][m];
	  flowpc=newdata[nit].flowpc[k][m];
	  flowpv=newdata[nit].flowpv[k][m];
	  flowvp=newdata[nit].flowvp[k][m];
	  flowvc=newdata[nit].flowvc[k][m];
	  flowcv=newdata[nit].flowcv[k][m];
	  flowsc=newdata[nit].flowsc[k][m];
	  flowcs=newdata[nit].flowcs[k][m];
	  flowsp=newdata[nit].flowsp[k][m];
	  flowps=newdata[nit].flowps[k][m];
	  flowcu=newdata[nit].flowcu[k][m];
	  flowpu=newdata[nit].flowpu[k][m];
	  flowsu=newdata[nit].flowsu[k][m];
	  flowvu=newdata[nit].flowvu[k][m];
	  flowuc=newdata[nit].flowuc[k][m];
	  flowup=newdata[nit].flowup[k][m];
	  flowus=newdata[nit].flowus[k][m];
          flowsbh=ZEROVALUE;
	  if(newdata[nit].smb[k][m] > ZEROVALUE) flowsbh=newdata[nit].sbh[k][m]/newdata[nit].smb[k][m]/garea[k][m];
          flowvbh=ZEROVALUE;
	  if(dstatic[k][m].vba > ZEROVALUE) flowvbh=newdata[nit].vbh[k][m]/dstatic[k][m].vba/garea[k][m];
          flowsbh2=ZEROVALUE;
	  if(newdata[nit].smb[k][m] > ZEROVALUE) flowsbh2=newdata[nit].sbh2[k][m]/newdata[nit].smb[k][m]/garea[k][m];
          flowvbh2=ZEROVALUE;
	  if(dstatic[k][m].vba > ZEROVALUE) flowvbh2=newdata[nit].vbh2[k][m]/dstatic[k][m].vba/garea[k][m];
	  flowsbh3=ZEROVALUE;
	  if(newdata[nit].smb[k][m] > ZEROVALUE) flowsbh3=newdata[nit].sbh3[k][m]/newdata[nit].smb[k][m]/garea[k][m];


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
          if(use_urban) {
	    if(ic == 13) fprintf(outfile,"%lf ",flowcu);
	    if(ic == 14) fprintf(outfile,"%lf ",flowpu);
	    if(ic == 15) fprintf(outfile,"%lf ",flowvu);
	    if(ic == 16) fprintf(outfile,"%lf ",flowsu);
	    if(ic == 17) fprintf(outfile,"%lf ",flowuc);
	    if(ic == 18) fprintf(outfile,"%lf ",flowup);
	    if(ic == 19) fprintf(outfile,"%lf ",flowus);
	  }
	}
        fprintf(outfile,"\n");

      }
      fclose(outfile);
    }
  }



  return;
}

/********************************************************/

void smart_flow(int k, int m, int it){

  double vstotal=ZEROVALUE;
  double delta_c=ZEROVALUE, delta_p=ZEROVALUE, delta_o=ZEROVALUE, delta_u=ZEROVALUE, unmet_du=ZEROVALUE;
  

  vstotal=newdata[it].v[k][m]+newdata[it].s[k][m]; 
  
  if (use_urban) {

/* Urban flows first */
   delta_u = newdata[it+1].u[k][m] - newdata[it].u[k][m];

   if(fabs(delta_u)<ZEROVALUE_CHECK){
     /*do nothing because delta_u is zero */ 
   }
   else if (delta_u<=ZEROVALUE){
     if((newdata[it+1].v[k][m] + newdata[it+1].s[k][m]) < -delta_u){
       newdata[it].flowus[k][m]=newdata[it+1].v[k][m] + newdata[it+1].s[k][m];
       unmet_du = -delta_u - newdata[it+1].v[k][m] - newdata[it+1].s[k][m];
       newdata[it].flowuc[k][m]=unmet_du*newdata[it].c[k][m]/(newdata[it].c[k][m] + newdata[it+1].p[k][m]);
       newdata[it].flowup[k][m]=unmet_du*newdata[it].p[k][m]/(newdata[it].c[k][m] + newdata[it+1].p[k][m]);
     }
     else newdata[it].flowus[k][m]=-delta_u;
   }
   else if (delta_u>ZEROVALUE){
     if((newdata[it].c[k][m]+newdata[it].p[k][m]+newdata[it].s[k][m])>=delta_u){
       newdata[it].flowcu[k][m]=newdata[it].c[k][m]/(newdata[it].c[k][m]+newdata[it].p[k][m]+newdata[it].s[k][m])*delta_u;
       newdata[it].flowpu[k][m]=newdata[it].p[k][m]/(newdata[it].c[k][m]+newdata[it].p[k][m]+newdata[it].s[k][m])*delta_u;
       newdata[it].flowsu[k][m]=newdata[it].s[k][m]/(newdata[it].c[k][m]+newdata[it].p[k][m]+newdata[it].s[k][m])*delta_u;
     }
     else{
       newdata[it].flowcu[k][m]=newdata[it].c[k][m];
       newdata[it].flowpu[k][m]=newdata[it].p[k][m];
       newdata[it].flowsu[k][m]=newdata[it].s[k][m];
       newdata[it].flowvu[k][m]=delta_u-newdata[it].c[k][m]-newdata[it].p[k][m]-newdata[it].s[k][m];
     }
   }
  }

  //  delta_c = newdata[it+1].c[k][m] - newdata[it].c[k][m];
  //  delta_p = newdata[it+1].p[k][m] - newdata[it].p[k][m];
  //  delta_o = newdata[it+1].v[k][m] + newdata[it+1].s[k][m] -newdata[it].v[k][m]-newdata[it].s[k][m];

  delta_c = delta_c+newdata[it].flowcu[k][m]-newdata[it].flowus[k][m];
  delta_p = delta_p+newdata[it].flowpu[k][m]-newdata[it].flowup[k][m];
  delta_o = delta_o-newdata[it].flowus[k][m]+newdata[it].flowsu[k][m]+newdata[it].flowvu[k][m];
  
  /*three cases*/
  if((fabs(delta_c)<ZEROVALUE_CHECK) && (fabs(delta_p)<ZEROVALUE_CHECK) && (fabs(delta_o)<ZEROVALUE_CHECK)){
    /*all zero*/
    
  }
  /*two negative*/
  else if((delta_c<=ZEROVALUE) && (delta_p<=ZEROVALUE) && (delta_o>ZEROVALUE)){
    
    newdata[it].flowps[k][m]=-delta_p;
    newdata[it].flowcs[k][m]=-delta_c;    
  }
  else if((delta_c<=ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o<=ZEROVALUE)){
    
    if(vstotal > ZEROVALUE){
      newdata[it].flowvp[k][m]=-delta_o*newdata[it].v[k][m]/vstotal;
      newdata[it].flowsp[k][m]=-delta_o*(1-newdata[it].v[k][m]/vstotal);
    }
    newdata[it].flowcp[k][m]=-delta_c;
    
  }
  else if((delta_c>ZEROVALUE)&& (delta_p<=ZEROVALUE) && (delta_o<=ZEROVALUE)){
    
    newdata[it].flowpc[k][m]=-delta_p; 
    if(vstotal > ZEROVALUE){
      newdata[it].flowvc[k][m]=-delta_o*newdata[it].v[k][m]/vstotal;
      newdata[it].flowsc[k][m]=-delta_o*(1-newdata[it].v[k][m]/vstotal);
    }
    
  } 
  /*two positive*/
  else if((delta_c<=ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o>ZEROVALUE)){

    newdata[it].flowcs[k][m]=delta_o;
    newdata[it].flowcp[k][m]=delta_p;
  }
  else if((delta_c>ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o<=ZEROVALUE)){

    if(vstotal > ZEROVALUE){
      newdata[it].flowvp[k][m]=delta_p*newdata[it].v[k][m]/vstotal;
      newdata[it].flowsp[k][m]=delta_p*(1-newdata[it].v[k][m]/vstotal);
      newdata[it].flowvc[k][m]=delta_c*newdata[it].v[k][m]/vstotal;
      newdata[it].flowsc[k][m]=delta_c*(1-newdata[it].v[k][m]/vstotal);
    }
    
  }
  else if((delta_c>ZEROVALUE)&& (delta_p<=ZEROVALUE) && (delta_o>ZEROVALUE)){

    newdata[it].flowps[k][m]=delta_o;
    newdata[it].flowpc[k][m]=delta_c;
    
  }  


  if (SMART_FLOW_BUG_PRINT){
  
    	/* checking for negative zero flows */

	if(newdata[it].flowcp[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowcp<0, resetting to 0 by force %10.15lf\n",newdata[it].flowcp[k][m]);newdata[it].flowcp[k][m] = ZEROVALUE;}
	if(newdata[it].flowpc[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowpc<0, resetting to 0 by force %10.15lf\n",newdata[it].flowpc[k][m]);newdata[it].flowpc[k][m] = ZEROVALUE;}
	if(newdata[it].flowpv[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowpv<0, resetting to 0 by force %10.15lf\n",newdata[it].flowpv[k][m]);newdata[it].flowpv[k][m] = ZEROVALUE;}
	if(newdata[it].flowvp[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowvp<0, resetting to 0 by force %10.15lf\n",newdata[it].flowvp[k][m]);newdata[it].flowvp[k][m] = ZEROVALUE;}
	if(newdata[it].flowvc[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowvc<0, resetting to 0 by force %10.15lf\n",newdata[it].flowvc[k][m]);newdata[it].flowvc[k][m] = ZEROVALUE;}
	if(newdata[it].flowcv[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowcv<0, resetting to 0 by force %10.15lf\n",newdata[it].flowcv[k][m]);newdata[it].flowcv[k][m] = ZEROVALUE;}
	if(newdata[it].flowsp[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowsp<0, resetting to 0 by force %10.15lf\n",newdata[it].flowsp[k][m]);newdata[it].flowsp[k][m] = ZEROVALUE;}
	if(newdata[it].flowps[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowps<0, resetting to 0 by force %10.15lf\n",newdata[it].flowps[k][m]);newdata[it].flowps[k][m] = ZEROVALUE;}
	if(newdata[it].flowsc[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowsc<0, resetting to 0 by force %10.15lf\n",newdata[it].flowsc[k][m]);newdata[it].flowsc[k][m] = ZEROVALUE;}
	if(newdata[it].flowcs[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowcs<0, resetting to 0 by force %10.15lf\n",newdata[it].flowcs[k][m]);newdata[it].flowcs[k][m] = ZEROVALUE;}
	if(newdata[it].flowvs[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowvs<0, resetting to 0 by force %10.15lf\n",newdata[it].flowvs[k][m]);newdata[it].flowvs[k][m] = ZEROVALUE;}
	if(newdata[it].flowcu[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowcu<0, resetting to 0 by force %10.15lf\n",newdata[it].flowcu[k][m]);newdata[it].flowcu[k][m] = ZEROVALUE;}
	if(newdata[it].flowpu[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowpu<0, resetting to 0 by force %10.15lf\n",newdata[it].flowpu[k][m]);newdata[it].flowpu[k][m] = ZEROVALUE;}
	if(newdata[it].flowvu[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowvu<0, resetting to 0 by force %10.15lf\n",newdata[it].flowvu[k][m]);newdata[it].flowvu[k][m] = ZEROVALUE;}
	if(newdata[it].flowsu[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowsu<0, resetting to 0 by force %10.15lf\n",newdata[it].flowsu[k][m]);newdata[it].flowsu[k][m] = ZEROVALUE;}
	if(newdata[it].flowuc[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowuc<0, resetting to 0 by force %10.15lf\n",newdata[it].flowuc[k][m]);newdata[it].flowuc[k][m] = ZEROVALUE;}
	if(newdata[it].flowup[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowup<0, resetting to 0 by force %10.15lf\n",newdata[it].flowup[k][m]);newdata[it].flowup[k][m] = ZEROVALUE;}
	if(newdata[it].flowus[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowus<0, resetting to 0 by force %10.15lf\n",newdata[it].flowus[k][m]);newdata[it].flowus[k][m] = ZEROVALUE;}



  }
    

  return;
  
}
 

/********************************************************/
void alternative_smart_flow1(int k, int m, int it){

  double vtotal=ZEROVALUE, stotal=ZEROVALUE, vstotal=ZEROVALUE;
  double p_other=ZEROVALUE, c_other=ZEROVALUE, total_demand=ZEROVALUE;
  double delta_c=ZEROVALUE, delta_p=ZEROVALUE, delta_o=ZEROVALUE, delta_u=ZEROVALUE, unmet_du=ZEROVALUE;

  if (use_urban) {
   /* Urban flows first */
   delta_u = newdata[it+1].u[k][m] - newdata[it].u[k][m];

   if(fabs(delta_u)<=ZEROVALUE_CHECK){
     /*do nothing because delta_u is zero */ 
   }
   else if (delta_u<=ZEROVALUE){
     newdata[it].flowuc[k][m]=newdata[it+1].c[k][m]/(newdata[it+1].c[k][m]+newdata[it+1].p[k][m]+newdata[it+1].s[k][m]+newdata[it+1].v[k][m])*(-delta_u);
     newdata[it].flowup[k][m]=newdata[it+1].p[k][m]/(newdata[it+1].c[k][m]+newdata[it+1].p[k][m]+newdata[it+1].s[k][m]+newdata[it+1].v[k][m])*(-delta_u);
     newdata[it].flowus[k][m]=(newdata[it+1].s[k][m]+newdata[it+1].v[k][m])/(newdata[it+1].c[k][m]+newdata[it+1].p[k][m]+newdata[it+1].s[k][m]+newdata[it+1].v[k][m])*(-delta_u);
   }
   else if (delta_u>ZEROVALUE){
     if((newdata[it].c[k][m]+newdata[it].p[k][m]+newdata[it].s[k][m])>=delta_u){
       newdata[it].flowcu[k][m]=newdata[it].c[k][m]/(newdata[it].c[k][m]+newdata[it].p[k][m]+newdata[it].s[k][m])*delta_u;
       newdata[it].flowpu[k][m]=newdata[it].p[k][m]/(newdata[it].c[k][m]+newdata[it].p[k][m]+newdata[it].s[k][m])*delta_u;
       newdata[it].flowsu[k][m]=newdata[it].s[k][m]/(newdata[it].c[k][m]+newdata[it].p[k][m]+newdata[it].s[k][m])*delta_u;
     }
     else{
       newdata[it].flowcu[k][m]=newdata[it].c[k][m];
       newdata[it].flowpu[k][m]=newdata[it].p[k][m];
       newdata[it].flowsu[k][m]=newdata[it].s[k][m];
       newdata[it].flowvu[k][m]=delta_u-newdata[it].c[k][m]-newdata[it].p[k][m]-newdata[it].s[k][m];
     }
   }
  
  }

  vtotal=newdata[it].v[k][m];     
  stotal=newdata[it].s[k][m];
  vstotal=vtotal+stotal;

  delta_c = newdata[it+1].c[k][m] - newdata[it].c[k][m];
  delta_p = newdata[it+1].p[k][m] - newdata[it].p[k][m];
  delta_o = newdata[it+1].v[k][m] + newdata[it+1].s[k][m] -newdata[it].v[k][m]-newdata[it].s[k][m];
 
      /*amend deltas */
  delta_c = delta_c+newdata[it].flowcu[k][m]-newdata[it].flowuc[k][m];
  delta_p = delta_p+newdata[it].flowpu[k][m]-newdata[it].flowup[k][m];
  delta_o = delta_o-newdata[it].flowus[k][m]+newdata[it].flowsu[k][m]+newdata[it].flowvu[k][m];
    
  /*three cases*/
  if((fabs(delta_c)<ZEROVALUE_CHECK) && (fabs(delta_p)<ZEROVALUE_CHECK) && (fabs(delta_o)<ZEROVALUE_CHECK)){
    /*all zero*/
    
  }
  /*two negative*/
  else if((delta_c<=ZEROVALUE) && (delta_p<=ZEROVALUE) && (delta_o>ZEROVALUE)){
    
    newdata[it].flowps[k][m]=-delta_p;
    newdata[it].flowcs[k][m]=-delta_c;    
  }
  else if((delta_c<=ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o<=ZEROVALUE)){

    newdata[it].flowcp[k][m]=-delta_c; 


    if((newdata[it].v[k][m] - newdata[it].flowvu[k][m]) >= -delta_o) {

      newdata[it].flowvp[k][m]=-delta_o;
    }
    else {
      newdata[it].flowvp[k][m]=newdata[it].v[k][m] - newdata[it].flowvu[k][m];
      newdata[it].flowsp[k][m]=-delta_o-(newdata[it].v[k][m] - newdata[it].flowvu[k][m]);
    }

    
  }
  else if((delta_c>ZEROVALUE)&& (delta_p<=ZEROVALUE) && (delta_o<=ZEROVALUE)){
    
    newdata[it].flowpc[k][m]=-delta_p; 


    if((newdata[it].v[k][m] - newdata[it].flowvu[k][m]) >= -delta_o) {

      newdata[it].flowvc[k][m]=-delta_o;
    }
    else {
      newdata[it].flowvc[k][m]=newdata[it].v[k][m];
      newdata[it].flowsc[k][m]=-delta_o-newdata[it].v[k][m];
    }
    

  } 
  /*two positive*/
  else if((delta_c<=ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o>ZEROVALUE)){

    newdata[it].flowcs[k][m]=delta_o;
    newdata[it].flowcp[k][m]=delta_p;
  }
  else if((delta_c>ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o<=ZEROVALUE)){

 
    if((newdata[it].v[k][m] - newdata[it].flowvu[k][m]) >= -delta_o) {

      newdata[it].flowvc[k][m]=delta_c;
      newdata[it].flowvp[k][m]=delta_p;

    }
    else {


      newdata[it].flowvc[k][m]=(newdata[it].v[k][m] - newdata[it].flowvu[k][m])*delta_c/(delta_p+delta_c);
      newdata[it].flowvp[k][m]=(newdata[it].v[k][m] - newdata[it].flowvu[k][m])*delta_p/(delta_p+delta_c);
    
      newdata[it].flowsc[k][m]=delta_c-newdata[it].flowvc[k][m];
      newdata[it].flowsp[k][m]=delta_p-newdata[it].flowvp[k][m];
    }
    

  }
  else if((delta_c>ZEROVALUE)&& (delta_p<=ZEROVALUE) && (delta_o>ZEROVALUE)){

    newdata[it].flowps[k][m]=delta_o;
    newdata[it].flowpc[k][m]=delta_c;
    
  } 
  


  if (SMART_FLOW_BUG_PRINT){


    	/* checking for negative zero flows */

	if(newdata[it].flowcp[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowcp<0, resetting to 0 by force %10.15lf\n",newdata[it].flowcp[k][m]);newdata[it].flowcp[k][m] = ZEROVALUE;}
	if(newdata[it].flowpc[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowpc<0, resetting to 0 by force %10.15lf\n",newdata[it].flowpc[k][m]);newdata[it].flowpc[k][m] = ZEROVALUE;}
	if(newdata[it].flowpv[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowpv<0, resetting to 0 by force %10.15lf\n",newdata[it].flowpv[k][m]);newdata[it].flowpv[k][m] = ZEROVALUE;}
	if(newdata[it].flowvp[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowvp<0, resetting to 0 by force %10.15lf\n",newdata[it].flowvp[k][m]);newdata[it].flowvp[k][m] = ZEROVALUE;}
	if(newdata[it].flowvc[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowvc<0, resetting to 0 by force %10.15lf\n",newdata[it].flowvc[k][m]);newdata[it].flowvc[k][m] = ZEROVALUE;}
	if(newdata[it].flowcv[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowcv<0, resetting to 0 by force %10.15lf\n",newdata[it].flowcv[k][m]);newdata[it].flowcv[k][m] = ZEROVALUE;}
	if(newdata[it].flowsp[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowsp<0, resetting to 0 by force %10.15lf\n",newdata[it].flowsp[k][m]);newdata[it].flowsp[k][m] = ZEROVALUE;}
	if(newdata[it].flowps[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowps<0, resetting to 0 by force %10.15lf\n",newdata[it].flowps[k][m]);newdata[it].flowps[k][m] = ZEROVALUE;}
	if(newdata[it].flowsc[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowsc<0, resetting to 0 by force %10.15lf\n",newdata[it].flowsc[k][m]);newdata[it].flowsc[k][m] = ZEROVALUE;}
	if(newdata[it].flowcs[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowcs<0, resetting to 0 by force %10.15lf\n",newdata[it].flowcs[k][m]);newdata[it].flowcs[k][m] = ZEROVALUE;}
	if(newdata[it].flowvs[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowvs<0, resetting to 0 by force %10.15lf\n",newdata[it].flowvs[k][m]);newdata[it].flowvs[k][m] = ZEROVALUE;}
	if(newdata[it].flowcu[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowcu<0, resetting to 0 by force %10.15lf\n",newdata[it].flowcu[k][m]);newdata[it].flowcu[k][m] = ZEROVALUE;}
	if(newdata[it].flowpu[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowpu<0, resetting to 0 by force %10.15lf\n",newdata[it].flowpu[k][m]);newdata[it].flowpu[k][m] = ZEROVALUE;}
	if(newdata[it].flowvu[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowvu<0, resetting to 0 by force %10.15lf\n",newdata[it].flowvu[k][m]);newdata[it].flowvu[k][m] = ZEROVALUE;}
	if(newdata[it].flowsu[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowsu<0, resetting to 0 by force %10.15lf\n",newdata[it].flowsu[k][m]);newdata[it].flowsu[k][m] = ZEROVALUE;}
	if(newdata[it].flowuc[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowuc<0, resetting to 0 by force %10.15lf\n",newdata[it].flowuc[k][m]);newdata[it].flowuc[k][m] = ZEROVALUE;}
	if(newdata[it].flowup[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowup<0, resetting to 0 by force %10.15lf\n",newdata[it].flowup[k][m]);newdata[it].flowup[k][m] = ZEROVALUE;}
	if(newdata[it].flowus[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowus<0, resetting to 0 by force %10.15lf\n",newdata[it].flowus[k][m]);newdata[it].flowus[k][m] = ZEROVALUE;}

  }
    


  return;
  
}


/********************************************************************/
void alternative_smart_flow2(int k, int m, int it){

  double p_other=ZEROVALUE, c_other=ZEROVALUE, total_demand=ZEROVALUE;
  double delta_c=ZEROVALUE, delta_p=ZEROVALUE, delta_o=ZEROVALUE, delta_u=ZEROVALUE, unmet_du=ZEROVALUE;

  
  if (use_urban) {
   /* Urban flows first */
   delta_u = newdata[it+1].u[k][m] - newdata[it].u[k][m];

   if(fabs(delta_u)<ZEROVALUE_CHECK){
     /*do nothing because delta_u is zero */ 
   }
   else if (delta_u<=ZEROVALUE){
     newdata[it].flowuc[k][m]=newdata[it+1].c[k][m]/(newdata[it+1].c[k][m]+newdata[it+1].p[k][m]+newdata[it+1].s[k][m]+newdata[it+1].v[k][m])*(-delta_u);
     newdata[it].flowup[k][m]=newdata[it+1].p[k][m]/(newdata[it+1].c[k][m]+newdata[it+1].p[k][m]+newdata[it+1].s[k][m]+newdata[it+1].v[k][m])*(-delta_u);
     newdata[it].flowus[k][m]=(newdata[it+1].s[k][m]+newdata[it+1].v[k][m])/(newdata[it+1].c[k][m]+newdata[it+1].p[k][m]+newdata[it+1].s[k][m]+newdata[it+1].v[k][m])*(-delta_u);
   }
   else if (delta_u>ZEROVALUE){
     if((newdata[it].c[k][m]+newdata[it].p[k][m]+newdata[it].s[k][m])>=delta_u){
       newdata[it].flowcu[k][m]=newdata[it].c[k][m]/(newdata[it].c[k][m]+newdata[it].p[k][m]+newdata[it].s[k][m])*delta_u;
       newdata[it].flowpu[k][m]=newdata[it].p[k][m]/(newdata[it].c[k][m]+newdata[it].p[k][m]+newdata[it].s[k][m])*delta_u;
       newdata[it].flowsu[k][m]=newdata[it].s[k][m]/(newdata[it].c[k][m]+newdata[it].p[k][m]+newdata[it].s[k][m])*delta_u;
     }
     else{
       newdata[it].flowcu[k][m]=newdata[it].c[k][m];
       newdata[it].flowpu[k][m]=newdata[it].p[k][m];
       newdata[it].flowsu[k][m]=newdata[it].s[k][m];
       newdata[it].flowvu[k][m]=delta_u-newdata[it].c[k][m]-newdata[it].p[k][m]-newdata[it].s[k][m];
     }
   }
  }

  delta_c = newdata[it+1].c[k][m] - newdata[it].c[k][m];
  delta_p = newdata[it+1].p[k][m] - newdata[it].p[k][m];
  delta_o = newdata[it+1].v[k][m] + newdata[it+1].s[k][m] -newdata[it].v[k][m]-newdata[it].s[k][m];

  /*amend deltas */
  delta_c = delta_c+newdata[it].flowcu[k][m]-newdata[it].flowuc[k][m];
  delta_p = delta_p+newdata[it].flowpu[k][m]-newdata[it].flowup[k][m];
  delta_o = delta_o-newdata[it].flowus[k][m]+newdata[it].flowsu[k][m]+newdata[it].flowvu[k][m];
  
  /*three cases*/
  if((fabs(delta_c)<ZEROVALUE_CHECK) && (fabs(delta_p)<ZEROVALUE_CHECK) && (fabs(delta_o)<ZEROVALUE_CHECK)){
    /*all zero*/
    
  }
  /*two negative*/
  else if((delta_c<=ZEROVALUE) && (delta_p<=ZEROVALUE) && (delta_o>ZEROVALUE)){
    
    newdata[it].flowps[k][m]=-delta_p;
    newdata[it].flowcs[k][m]=-delta_c;    
  }
  else if((delta_c<=ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o<=ZEROVALUE)){

    newdata[it].flowcp[k][m]=-delta_c; 


    if((newdata[it].s[k][m] - newdata[it].flowsu[k][m]) >= -delta_o) {

      newdata[it].flowsp[k][m]=-delta_o;
    }
    else {
      newdata[it].flowsp[k][m]=(newdata[it].s[k][m] - newdata[it].flowsu[k][m]);
      newdata[it].flowvp[k][m]=-delta_o-(newdata[it].s[k][m] - newdata[it].flowsu[k][m]);
    }

    
  }
  else if((delta_c>ZEROVALUE)&& (delta_p<=ZEROVALUE) && (delta_o<=ZEROVALUE)){
    
    newdata[it].flowpc[k][m]=-delta_p; 


    if((newdata[it].s[k][m] - newdata[it].flowsu[k][m]) >= -delta_o) {

      newdata[it].flowsc[k][m]=-delta_o;
    }
    else {
      newdata[it].flowsc[k][m]=(newdata[it].s[k][m] - newdata[it].flowsu[k][m]);
      newdata[it].flowvc[k][m]=-delta_o-(newdata[it].s[k][m] - newdata[it].flowsu[k][m]);
    }
    

  } 
  /*two positive*/
  else if((delta_c<=ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o>ZEROVALUE)){

    newdata[it].flowcs[k][m]=delta_o;
    newdata[it].flowcp[k][m]=delta_p;
  }
  else if((delta_c>ZEROVALUE) && (delta_p>ZEROVALUE) && (delta_o<=ZEROVALUE)){

 
    if((newdata[it].s[k][m] - newdata[it].flowsu[k][m]) >= -delta_o) {

      newdata[it].flowsc[k][m]=delta_c;
      newdata[it].flowsp[k][m]=delta_p;

    }
    else {

      newdata[it].flowsc[k][m]=(newdata[it].s[k][m] - newdata[it].flowsu[k][m])*delta_c/(delta_p+delta_c);
      newdata[it].flowsp[k][m]=(newdata[it].s[k][m] - newdata[it].flowsu[k][m])*delta_p/(delta_p+delta_c);
    
      newdata[it].flowvc[k][m]=delta_c-newdata[it].flowsc[k][m];
      newdata[it].flowvp[k][m]=delta_p-newdata[it].flowsp[k][m];
    }
    

  }
  else if((delta_c>ZEROVALUE)&& (delta_p<=ZEROVALUE) && (delta_o>ZEROVALUE)){

    newdata[it].flowps[k][m]=delta_o;
    newdata[it].flowpc[k][m]=delta_c;
    
  } 
  

  if (SMART_FLOW_BUG_PRINT){

 	if((newdata[it].flowvc[k][m] < ZEROVALUE) || (newdata[it].flowvp[k][m] < ZEROVALUE) || (newdata[it].flowsc[k][m] < ZEROVALUE) || (newdata[it].flowsp[k][m] < ZEROVALUE)) printf("v(t) %lf s(t) %lf c(t) %lf p(t) %lf v(t+1) %lf s(t+1) %lf c(t+1) %lf p(t+1) %lf k %d m %d vc %lf vp %lf sc %lf sp %lf delta_p %lf delta_o %lf delta_c %lf\n",newdata[it].v[k][m],newdata[it].s[k][m],newdata[it].c[k][m],newdata[it].p[k][m],newdata[it+1].v[k][m],newdata[it+1].s[k][m],newdata[it+1].c[k][m],newdata[it+1].p[k][m],k,m,newdata[it].flowvc[k][m],newdata[it].flowvp[k][m],newdata[it].flowsc[k][m],newdata[it].flowsp[k][m],delta_p,delta_o,delta_c); 
  
    	/* checking for negative zero flows */

	if(newdata[it].flowcp[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowcp<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[it].flowcp[k][m]);newdata[it].flowcp[k][m] = ZEROVALUE;}
	if(newdata[it].flowpc[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowpc<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[it].flowpc[k][m]);newdata[it].flowpc[k][m] = ZEROVALUE;}
	if(newdata[it].flowpv[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowpv<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[it].flowpv[k][m]);newdata[it].flowpv[k][m] = ZEROVALUE;}
	if(newdata[it].flowvp[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowvp<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[it].flowvp[k][m]);newdata[it].flowvp[k][m] = ZEROVALUE;}
	if(newdata[it].flowvc[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowvc<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[it].flowvc[k][m]);newdata[it].flowvc[k][m] = ZEROVALUE;}
	if(newdata[it].flowcv[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowcv<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[it].flowcv[k][m]);newdata[it].flowcv[k][m] = ZEROVALUE;}
	if(newdata[it].flowsp[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowsp<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[it].flowsp[k][m]);newdata[it].flowsp[k][m] = ZEROVALUE;}
	if(newdata[it].flowps[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowps<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[it].flowps[k][m]);newdata[it].flowps[k][m] = ZEROVALUE;}
	if(newdata[it].flowsc[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowsc<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[it].flowsc[k][m]);newdata[it].flowsc[k][m] = ZEROVALUE;}
	if(newdata[it].flowcs[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowcs<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[it].flowcs[k][m]);newdata[it].flowcs[k][m] = ZEROVALUE;}
	if(newdata[it].flowvs[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowvs<0, resetting to 0 by force k %d m %d %10.15lf\n",k,m,newdata[it].flowvs[k][m]);newdata[it].flowvs[k][m] = ZEROVALUE;}
	if(newdata[it].flowcu[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowcu<0, resetting to 0 by force %10.15lf\n",newdata[it].flowcu);newdata[it].flowcu[k][m] = ZEROVALUE;}
	if(newdata[it].flowpu[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowpu<0, resetting to 0 by force %10.15lf\n",newdata[it].flowpu);newdata[it].flowpu[k][m] = ZEROVALUE;}
	if(newdata[it].flowvu[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowvu<0, resetting to 0 by force %10.15lf\n",newdata[it].flowvu);newdata[it].flowvu[k][m] = ZEROVALUE;}
	if(newdata[it].flowsu[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowsu<0, resetting to 0 by force %10.15lf\n",newdata[it].flowsu);newdata[it].flowsu[k][m] = ZEROVALUE;}
	if(newdata[it].flowuc[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowuc<0, resetting to 0 by force %10.15lf\n",newdata[it].flowuc);newdata[it].flowuc[k][m] = ZEROVALUE;}
	if(newdata[it].flowup[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowup<0, resetting to 0 by force %10.15lf\n",newdata[it].flowup);newdata[it].flowup[k][m] = ZEROVALUE;}
	if(newdata[it].flowus[k][m] < ZEROVALUE) {printf("Smart flow Bug: flowus<0, resetting to 0 by force %10.15lf\n",newdata[it].flowus);newdata[it].flowus[k][m] = ZEROVALUE;}


  }
   

  return;
  
}

/********************************************************************/
void adjust_smart_flow1(int k, int m, int it, int i){

  double vtotal=ZEROVALUE, stotal=ZEROVALUE, vstotal=ZEROVALUE;
  double flowcs_prime=ZEROVALUE, flowps_prime=ZEROVALUE;
  double flowvc_prime=ZEROVALUE, flowvp_prime=ZEROVALUE;
  double flowsc_prime=ZEROVALUE, flowsp_prime=ZEROVALUE;
  double time_k, crop, past;


  vtotal=newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m]; 
 
  stotal=newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m];

  vstotal=vtotal+stotal;

  crop=newdata[it].c[k][m]-newdata[it].flowcs[k][m]-newdata[it].flowcp[k][m]-newdata[it].flowcu[k][m];

  past=newdata[it].p[k][m]-newdata[it].flowps[k][m]-newdata[it].flowpc[k][m]-newdata[it].flowpu[k][m];

 
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


 
    rtdata[i].flowvc_prime+=flowvc_prime*garea[k][m]/1000./1000.;
    rtdata[i].flowvp_prime+=flowvp_prime*garea[k][m]/1000./1000.;
    rtdata[i].flowsp_prime+=flowsp_prime*garea[k][m]/1000./1000.;
    rtdata[i].flowsc_prime+=flowsc_prime*garea[k][m]/1000./1000.;
    rtdata[i].flowps_prime+=flowps_prime*garea[k][m]/1000./1000.;
    rtdata[i].flowcs_prime+=flowcs_prime*garea[k][m]/1000./1000.;  

    newdata[it].flowvp[k][m]+=flowvp_prime;
    newdata[it].flowvc[k][m]+=flowvc_prime;
    newdata[it].flowsp[k][m]+=flowsp_prime;
    newdata[it].flowsc[k][m]+=flowsc_prime;
    newdata[it].flowps[k][m]+=flowps_prime;
    newdata[it].flowcs[k][m]+=flowcs_prime;




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
  vtotal=newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m]; 
 
  stotal=newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m];

  vstotal=vtotal+stotal;

  crop=newdata[it].c[k][m]-newdata[it].flowcs[k][m]-newdata[it].flowcp[k][m]-newdata[it].flowcu[k][m];

  past=newdata[it].p[k][m]-newdata[it].flowps[k][m]-newdata[it].flowpc[k][m]-newdata[it].flowpu[k][m];
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

    rtdata[i].flowvc_prime+=flowvc_prime*garea[k][m]/1000./1000.;
    rtdata[i].flowvp_prime+=flowvp_prime*garea[k][m]/1000./1000.;
    rtdata[i].flowsp_prime+=flowsp_prime*garea[k][m]/1000./1000.;
    rtdata[i].flowsc_prime+=flowsc_prime*garea[k][m]/1000./1000.;
    rtdata[i].flowps_prime+=flowps_prime*garea[k][m]/1000./1000.;
    rtdata[i].flowcs_prime+=flowcs_prime*garea[k][m]/1000./1000.; 

    newdata[it].flowvp[k][m]+=flowvp_prime;
    newdata[it].flowvc[k][m]+=flowvc_prime;
    newdata[it].flowsp[k][m]+=flowsp_prime;
    newdata[it].flowsc[k][m]+=flowsc_prime;
    newdata[it].flowps[k][m]+=flowps_prime;
    newdata[it].flowcs[k][m]+=flowcs_prime;




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


  vtotal=newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m]; 
    
  stotal=newdata[it].s[k][m]+newdata[it].flowcs[k][m]+newdata[it].flowps[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m];

  vstotal=vtotal+stotal;

  crop=newdata[it].c[k][m]+newdata[it].flowvc[k][m]+newdata[it].flowsc[k][m]+newdata[it].flowpc[k][m]-newdata[it].flowcs[k][m]-newdata[it].flowcp[k][m]-newdata[it].flowcu[k][m];

  past=newdata[it].p[k][m]+newdata[it].flowvp[k][m]+newdata[it].flowsp[k][m]+newdata[it].flowcp[k][m]-newdata[it].flowps[k][m]-newdata[it].flowpc[k][m]-newdata[it].flowpu[k][m];


  
  /* FROLKING METHOD MARCH 23, 2004, Ensures secondary is fallowed before it is re-cultivated */

  if( (vstotal > ZEROVALUE) && ((crop > ZEROVALUE) && (crop < 0.10)) ){

    /* time before abandoning cultivated land */
    time_k=4.0;
    
    /* fraction v+s required to replace crop abandoned */
    
    if(newdata[it].sma[k][m] < 1.0) {
      time_fallow=time_k;
    }
    else{
      time_fallow=vstotal/(crop/time_k+stotal/(2.0*newdata[it].sma[k][m]-1.0));
    }
    
    /* setting a max, allow us to return to fallowed secondary every 20 years */
    time_fallow_max=20.0;
    
    /* setting a min, allow us to return to fallowed secondary every time_k years */
    time_fallow_min=time_k;

    
    if(time_fallow > time_fallow_max){
      
      time_fallow=time_fallow_max;
      
      if(newdata[it].sma[k][m] < 1.0) {
	area_secondary_fallowed=ZEROVALUE;
      }
      else{
	area_secondary_fallowed=stotal*(1.0-time_fallow/(2.0*newdata[it].sma[k][m]-1.0));
      }
      
    }
    else if(time_fallow > time_k){
      
      if(newdata[it].sma[k][m] < 1.0) {
	area_secondary_fallowed=ZEROVALUE;
      }
      else{
	area_secondary_fallowed=stotal*(1.0-time_fallow/(2.0*newdata[it].sma[k][m]-1.0));
      }
      
    }
    else{ /*less than 4.0, slow down abandonment */
      
      time_fallow=time_k;    
      
      if(newdata[it].sma[k][m] < 1.0) {
	area_secondary_fallowed=ZEROVALUE;
      }
      else{
	area_secondary_fallowed=stotal*(1.0-time_fallow/(2.0*newdata[it].sma[k][m]-1.0));
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
    
    
    
    
    rtdata[i].flowvc_prime+=flowvc_prime*garea[k][m]/1000.0/1000.0;
    rtdata[i].flowsc_prime+=flowsc_prime*garea[k][m]/1000.0/1000.0;
    rtdata[i].flowcs_prime+=flowcs_prime*garea[k][m]/1000.0/1000.0;
    
    newdata[it].flowvc[k][m]+=flowvc_prime;
    newdata[it].flowsc[k][m]+=flowsc_prime;
    newdata[it].flowcs[k][m]+=flowcs_prime;



  } /* vstotal > zerovalue endif */

  return;
  
}

/********************************************************************/

void initialize_checker(int regional_code){

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

      if((CONTERMINOUS == 1) && (regional_code == 221)){

	if(lat[k] > 50.) dstatic[k][m].usflag=0;
	if((lat[k] > -60.)) dstatic[k][m].usflag=0;
	if(lon[m] < -130.) dstatic[k][m].usflag=0;
	  
      }

    }
  }



  for (it=0;it<NEWNT;it++){

    byclass_sum.c=ZEROVALUE;
    byclass_sum.p=ZEROVALUE;
    byclass_sum.v=ZEROVALUE;
    byclass_sum.i=ZEROVALUE;
    byclass_sum.w=ZEROVALUE;
    byclass_sum.s=ZEROVALUE;
    byclass_sum.u=ZEROVALUE;
    byclass_sum.total=ZEROVALUE;
    byclass_sum.vf=ZEROVALUE;
    byclass_sum.vnf=ZEROVALUE;
    byclass_sum.sf=ZEROVALUE;
    byclass_sum.snf=ZEROVALUE;
    byclass_sum.smaf=ZEROVALUE;
    byclass_sum.smanf=ZEROVALUE;
    byclass_sum.vbh=ZEROVALUE;
    byclass_sum.sbh=ZEROVALUE;
    byclass_sum.vbh2=ZEROVALUE;
    byclass_sum.sbh2=ZEROVALUE;
    byclass_sum.sbh3=ZEROVALUE;
    byclass_sum.wh_unmet=ZEROVALUE;
    byclass_sum.wh=ZEROVALUE;
    byclass_sum.flowcp=ZEROVALUE;
    byclass_sum.flowpc=ZEROVALUE;
    byclass_sum.flowpv=ZEROVALUE;
    byclass_sum.flowvp=ZEROVALUE;
    byclass_sum.flowvc=ZEROVALUE;
    byclass_sum.flowcv=ZEROVALUE;
    byclass_sum.flowsp=ZEROVALUE;
    byclass_sum.flowps=ZEROVALUE;
    byclass_sum.flowsc=ZEROVALUE;
    byclass_sum.flowcs=ZEROVALUE;
    byclass_sum.flowcu=ZEROVALUE;
    byclass_sum.flowpu=ZEROVALUE;
    byclass_sum.flowvu=ZEROVALUE;
    byclass_sum.flowsu=ZEROVALUE;
    byclass_sum.flowuc=ZEROVALUE;
    byclass_sum.flowup=ZEROVALUE;
    byclass_sum.flowus=ZEROVALUE;
    byclass_sum.flowsbh=ZEROVALUE;
    byclass_sum.flowsbh2=ZEROVALUE;
    byclass_sum.flowsbh3=ZEROVALUE;
    byclass_sum.flowvbh=ZEROVALUE;
    byclass_sum.flowvbh2=ZEROVALUE;
    byclass_sum.converted_forest_land=ZEROVALUE;
    byclass_sum.converted_forest_land_total=ZEROVALUE;
    byclass_sum.flowcs_prime=ZEROVALUE;
    byclass_sum.flowvc_prime=ZEROVALUE;
    byclass_sum.flowsc_prime=ZEROVALUE;
    byclass_sum.flowvp_prime=ZEROVALUE;
    byclass_sum.flowsp_prime=ZEROVALUE;
    byclass_sum.flowps_prime=ZEROVALUE;
    byclass_sum.flowsu_prime=ZEROVALUE;
    byclass_sum.flowvu_prime=ZEROVALUE;
    byclass_sum.flowus_prime=ZEROVALUE;
  }

  return;
}
/********************************************************************/

 void global_timeseries_checker(int regional_code, char regional_name[50], int curryear){

  FILE *outfile, *outfile2;
  int i, ic, k, m, it;
  double garea_sum=ZEROVALUE, land_garea_sum=ZEROVALUE, total_garea_sum=ZEROVALUE;
  double nf_s_area=ZEROVALUE, fnf_s_area=ZEROVALUE;
  char ffname[90];


 
  printf("printing global tester...\n"); 

  initialize_checker(regional_code);
 

  strcpy(ffname,regional_name);  
  strcat(ffname,".test");

  outfile=fopen(ffname,outstat);
  outfile2=fopen("global.primeflow.txt",outstat);

  it=0;
    
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){   
	
#if REGION_TEST
	
	if((dstatic[k][m].rcode==regional_code) && (dstatic[k][m].usflag==1)){   	
#else
	   if(dstatic[k][m].gcode > 0){ 
#endif		
	   
	    if(it == 0) {
	      chtest[k][m][it].c = newdata[it].c[k][m];
	      chtest[k][m][it].p = newdata[it].p[k][m];
	      chtest[k][m][it].v = newdata[it].v[k][m];
	      chtest[k][m][it].s = newdata[it].s[k][m];
	      chtest[k][m][it].u = newdata[it].u[k][m];
	    }
	   
	    chtest[k][m][it+1].c=chtest[k][m][it].c +
	      newdata[it].flowpc[k][m] + newdata[it].flowvc[k][m] + 
	      newdata[it].flowsc[k][m] - newdata[it].flowcp[k][m] - 
	      newdata[it].flowcv[k][m] - newdata[it].flowcs[k][m] - newdata[it].flowcu[k][m] + newdata[it].flowuc[k][m]; 
	    
	    chtest[k][m][it+1].p=chtest[k][m][it].p +
	      newdata[it].flowcp[k][m] + newdata[it].flowvp[k][m] + 
	      newdata[it].flowsp[k][m] - newdata[it].flowpc[k][m] - 
	      newdata[it].flowpv[k][m] - newdata[it].flowps[k][m] - newdata[it].flowpu[k][m] + newdata[it].flowup[k][m];
	    
	    chtest[k][m][it+1].v=chtest[k][m][it].v +
	      newdata[it].flowpv[k][m] + newdata[it].flowcv[k][m] -
	      newdata[it].flowvp[k][m] - newdata[it].flowvc[k][m] -   
	      newdata[it].flowvs[k][m] - newdata[it].flowvu[k][m];
	    
	    chtest[k][m][it+1].s=chtest[k][m][it].s +
	      newdata[it].flowvs[k][m] + newdata[it].flowps[k][m] + 
	      newdata[it].flowcs[k][m] -
	      newdata[it].flowsp[k][m] - newdata[it].flowsc[k][m] -  
	      newdata[it].flowsu[k][m] + newdata[it].flowus[k][m];  
	    
	    chtest[k][m][it+1].u=chtest[k][m][it].u +
	      newdata[it].flowcu[k][m] + newdata[it].flowvu[k][m] + 
              newdata[it].flowsu[k][m] + newdata[it].flowpu[k][m] - 
	      newdata[it].flowuc[k][m] - newdata[it].flowus[k][m] -
	      newdata[it].flowup[k][m];  
	    
	    


	  } /* end of region_test/rcode */

	} /* end of m */
      } /* end of m */

    
    for (i=0;i<NREG;i++) {

#if REGION_TEST
	
      if(rdata[i].rcode == regional_code) {

	byclass_sum.wh=rtdata[i].wh;
	byclass_sum.converted_forest_land=rtdata[i].converted_forest_land;
	byclass_sum.converted_forest_land_total=rtdata[i].converted_forest_land_total;
	byclass_sum.flowcs_prime+=rtdata[i].flowcs_prime;
	byclass_sum.flowvc_prime+=rtdata[i].flowvc_prime;
	byclass_sum.flowsc_prime+=rtdata[i].flowsc_prime;
	byclass_sum.flowvp_prime+=rtdata[i].flowvp_prime;
	byclass_sum.flowsp_prime+=rtdata[i].flowsp_prime;
	byclass_sum.flowps_prime+=rtdata[i].flowps_prime;
	byclass_sum.flowus_prime+=rtdata[i].flowus_prime;
	byclass_sum.flowsu_prime+=rtdata[i].flowsu_prime;
	byclass_sum.flowvu_prime+=rtdata[i].flowvu_prime;

      }
#else

      byclass_sum.wh+=rtdata[i].wh;

      if(rdata[i].rcode > 0){

	byclass_sum.flowcs_prime+=rtdata[i].flowcs_prime;
	byclass_sum.flowvc_prime+=rtdata[i].flowvc_prime;
	byclass_sum.flowsc_prime+=rtdata[i].flowsc_prime;
	byclass_sum.flowvp_prime+=rtdata[i].flowvp_prime;
	byclass_sum.flowsp_prime+=rtdata[i].flowsp_prime;
	byclass_sum.flowps_prime+=rtdata[i].flowps_prime;
	byclass_sum.flowus_prime+=rtdata[i].flowus_prime;
	byclass_sum.flowsu_prime+=rtdata[i].flowsu_prime;
	byclass_sum.flowvu_prime+=rtdata[i].flowvu_prime;

      }
#endif
   
    } /* end of i */

    for (i=0;i<NCCODE;i++) {
      byclass_sum.converted_forest_land+=ctdata[i].converted_forest_land;
      byclass_sum.converted_forest_land_total+=ctdata[i].converted_forest_land_total;
    }


    fnf_s_area=ZEROVALUE;
    nf_s_area=ZEROVALUE;
    
    if(gridded_wh==1){
      byclass_sum.wh=ZEROVALUE;
      byclass_sum.converted_forest_land=ZEROVALUE;
      byclass_sum.converted_forest_land_total=ZEROVALUE;
    }
    
    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){		

#if REGION_TEST
	
	if((dstatic[k][m].rcode==regional_code) && (dstatic[k][m].usflag==1)){   	
#else
	    if(dstatic[k][m].gcode > 0){ 
#endif		
	      if(gridded_wh==1){
	      byclass_sum.wh+=newdata[it].whd[k][m];
	      byclass_sum.converted_forest_land+=newdata[it].converted_forest_land[k][m];
	      byclass_sum.converted_forest_land_total+=newdata[it].converted_forest_land_total[k][m];
	      }
	     	     
	    if(dstatic[k][m].fnf == 1){
		
	      byclass_sum.vf+=garea[k][m]*chtest[k][m][it].v;
	      byclass_sum.sf+=garea[k][m]*chtest[k][m][it].s;
	      byclass_sum.smaf+=chtest[k][m][it].s*garea[k][m]*newdata[it].sma[k][m]; 
	      fnf_s_area += chtest[k][m][it].s*garea[k][m];
	    }
	    else{
	      byclass_sum.vnf+=garea[k][m]*chtest[k][m][it].v;
	      byclass_sum.snf+=garea[k][m]*chtest[k][m][it].s;
	      byclass_sum.smanf+=chtest[k][m][it].s*garea[k][m]*newdata[it].sma[k][m];
	      nf_s_area += chtest[k][m][it].s*garea[k][m];
	    }
	    
	    byclass_sum.vbh+=newdata[it].vbh[k][m];
	    byclass_sum.sbh+=newdata[it].sbh[k][m];
	    byclass_sum.vbh2+=newdata[it].vbh2[k][m];
	    byclass_sum.sbh2+=newdata[it].sbh2[k][m];
	    byclass_sum.sbh3+=newdata[it].sbh3[k][m];
	    byclass_sum.wh_unmet+=ZEROVALUE; 
	    
	    
	    byclass_sum.c+=garea[k][m]*chtest[k][m][it].c; 
	    byclass_sum.p+=garea[k][m]*chtest[k][m][it].p;
	    byclass_sum.v+=garea[k][m]*chtest[k][m][it].v; 
	    byclass_sum.i+=garea[k][m]*newdata[it].i[k][m];  
	    byclass_sum.w+=garea[k][m]*newdata[it].w[k][m];  
	    byclass_sum.s+=garea[k][m]*chtest[k][m][it].s; 
	    byclass_sum.u+=garea[k][m]*chtest[k][m][it].u; 

	   

	    byclass_sum.flowcp+=garea[k][m]/1000.0/1000.0*newdata[it].flowcp[k][m];
	    byclass_sum.flowpc+=garea[k][m]/1000.0/1000.0*newdata[it].flowpc[k][m];
	    byclass_sum.flowpv+=garea[k][m]/1000.0/1000.0*newdata[it].flowpv[k][m];
	    byclass_sum.flowvp+=garea[k][m]/1000.0/1000.0*newdata[it].flowvp[k][m];
	    byclass_sum.flowvc+=garea[k][m]/1000.0/1000.0*newdata[it].flowvc[k][m];
	    byclass_sum.flowcv+=garea[k][m]/1000.0/1000.0*newdata[it].flowcv[k][m];
	    byclass_sum.flowsp+=garea[k][m]/1000.0/1000.0*newdata[it].flowsp[k][m];
	    byclass_sum.flowps+=garea[k][m]/1000.0/1000.0*newdata[it].flowps[k][m];
	    byclass_sum.flowsc+=garea[k][m]/1000.0/1000.0*newdata[it].flowsc[k][m];
	    byclass_sum.flowcs+=garea[k][m]/1000.0/1000.0*newdata[it].flowcs[k][m];
	    byclass_sum.flowcu+=garea[k][m]/1000.0/1000.0*newdata[it].flowcu[k][m];
	    byclass_sum.flowpu+=garea[k][m]/1000.0/1000.0*newdata[it].flowpu[k][m];
	    byclass_sum.flowvu+=garea[k][m]/1000.0/1000.0*newdata[it].flowvu[k][m];
	    byclass_sum.flowsu+=garea[k][m]/1000.0/1000.0*newdata[it].flowsu[k][m];
	    byclass_sum.flowuc+=garea[k][m]/1000.0/1000.0*newdata[it].flowuc[k][m];
	    byclass_sum.flowup+=garea[k][m]/1000.0/1000.0*newdata[it].flowup[k][m];
	    byclass_sum.flowus+=garea[k][m]/1000.0/1000.0*newdata[it].flowus[k][m];


	    if(dstatic[k][m].fnf == 1){

	      byclass_sum.flowvbh+=newdata[it].vbh[k][m]/dstatic[k][m].vba/1000.0/1000.0;
	      
	      if(newdata[it].smb[k][m] > ZEROVALUE){
		
		byclass_sum.flowsbh+=newdata[it].sbh[k][m]/newdata[it].smb[k][m]/1000.0/1000.0;
		byclass_sum.flowsbh2+=newdata[it].sbh2[k][m]/newdata[it].smb[k][m]/1000.0/1000.0;
	      }
	      
	    }
	    else{
	      if(dstatic[k][m].vba > ZEROVALUE) byclass_sum.flowvbh2+=newdata[it].vbh2[k][m]/dstatic[k][m].vba/1000.0/1000.0; 
	      if(newdata[it].smb[k][m] > ZEROVALUE) byclass_sum.flowsbh3+=newdata[it].sbh3[k][m]/newdata[it].smb[k][m]/1000.0/1000.0;
	    }
	    

	    
	    if(it == 0) land_garea_sum+=(chtest[k][m][it].v+chtest[k][m][it].s+chtest[k][m][it].p+chtest[k][m][it].c+chtest[k][m][it].u)*garea[k][m]; 
	    
	    if(it == 0) garea_sum+=garea[k][m];
	      

	    } /* gcode */

	    if(it==0) total_garea_sum+=garea[k][m];

	  
	} /* end of m*/
      } /* end of k*/
   
      


      byclass_sum.total+=byclass_sum.c+byclass_sum.p+
	byclass_sum.v+byclass_sum.i+
	byclass_sum.w+byclass_sum.s+byclass_sum.u;  

     
      
      if(fnf_s_area > ZEROVALUE) byclass_sum.smaf=byclass_sum.smaf/fnf_s_area;
    
      if(nf_s_area > ZEROVALUE) byclass_sum.smanf=byclass_sum.smanf/nf_s_area;
    
 
      if (use_urban) {

	if(print_file_headers) fprintf(outfile,"yr c p v i w s u fsum vf vnf sf snf smaf smanf vbh sbh vbh2 sbh2 sbh3 sum wh clearing_amount_total clearing_amount_used wh-sum-clearing_amount_ifused(unmet) flowcp flowpc flowpv flowvp flowvc flowcv flowsp flowps flowsc flowcs flowcu flowpu flowvu flowsu flowuc flowup flowus flowvbh flowsbh flowvbh2 flowsbh2 flowsbh3 land_asum_km2 land_iw_asum_km2 total_asum_km2\n");
      

      fprintf(outfile,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	      curryear,
	      byclass_sum.c/garea_sum,
	      byclass_sum.p/garea_sum,
	      byclass_sum.v/garea_sum,
	      byclass_sum.i/garea_sum,
	      byclass_sum.w/garea_sum,
	      byclass_sum.s/garea_sum,
	      byclass_sum.u/garea_sum,
	      byclass_sum.total/garea_sum,
	      byclass_sum.vf/garea_sum,
	      byclass_sum.vnf/garea_sum,
	      byclass_sum.sf/garea_sum,
	      byclass_sum.snf/garea_sum,
	      byclass_sum.smaf,
	      byclass_sum.smanf,
	      byclass_sum.vbh/1000.,
	      byclass_sum.sbh/1000.,
	      byclass_sum.vbh2/1000.,
	      byclass_sum.sbh2/1000.,
	      byclass_sum.sbh3/1000.,
	      (byclass_sum.vbh+byclass_sum.sbh+byclass_sum.vbh2+byclass_sum.sbh2+byclass_sum.sbh3)/1000.,
              byclass_sum.wh,
	      byclass_sum.converted_forest_land_total,
	      byclass_sum.converted_forest_land,
	      byclass_sum.wh-((byclass_sum.vbh+byclass_sum.sbh+byclass_sum.vbh2+byclass_sum.sbh2+byclass_sum.sbh3)/1000.+byclass_sum.converted_forest_land),
	      byclass_sum.flowcp,
	      byclass_sum.flowpc,
	      byclass_sum.flowpv,
	      byclass_sum.flowvp,
	      byclass_sum.flowvc,
	      byclass_sum.flowcv,
	      byclass_sum.flowsp,
	      byclass_sum.flowps,
	      byclass_sum.flowsc,
	      byclass_sum.flowcs,
	      byclass_sum.flowcu,
	      byclass_sum.flowpu,
	      byclass_sum.flowvu,
	      byclass_sum.flowsu,
	      byclass_sum.flowuc,
	      byclass_sum.flowup,
	      byclass_sum.flowus,
	      byclass_sum.flowvbh,
	      byclass_sum.flowsbh,
	      byclass_sum.flowvbh2,
	      byclass_sum.flowsbh2,
	      byclass_sum.flowsbh3,
	      land_garea_sum/1000.0/1000.0,
	      garea_sum/1000.0/1000.0,
	      total_garea_sum/1000.0/1000.0); 

      if(print_file_headers) fprintf(outfile2,"yr cs vc sc vp sp ps vu su us\n");
      
      fprintf(outfile2,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	      curryear,
	      byclass_sum.flowcs_prime,
	      byclass_sum.flowvc_prime,
	      byclass_sum.flowsc_prime,
	      byclass_sum.flowvp_prime,
	      byclass_sum.flowsp_prime,
	      byclass_sum.flowps_prime,
	      byclass_sum.flowvu_prime,
	      byclass_sum.flowsu_prime,
	      byclass_sum.flowus_prime);

      }else{

        if(print_file_headers) fprintf(outfile,"yr c p v i w s fsum vf vnf sf snf smaf smanf vbh sbh vbh2 sbh2 sbh3 sum wh clearing_amount_total clearing_amount_used wh-sum-clearing_amount_ifused(unmet) flowcp[k][m] flowpc flowpv flowvp flowvc flowcv flowsp flowps flowsc flowcs flowvbh flowsbh flowvbh2 flowsbh2 flowsbh3 land_asum_km2 land_iw_asum_km2 total_asum_km2\n");
      

      fprintf(outfile,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	      curryear,
	      byclass_sum.c/garea_sum,
	      byclass_sum.p/garea_sum,
	      byclass_sum.v/garea_sum,
	      byclass_sum.i/garea_sum,
	      byclass_sum.w/garea_sum,
	      byclass_sum.s/garea_sum,
	      byclass_sum.total/garea_sum,
	      byclass_sum.vf/garea_sum,
	      byclass_sum.vnf/garea_sum,
	      byclass_sum.sf/garea_sum,
	      byclass_sum.snf/garea_sum,
	      byclass_sum.smaf,
	      byclass_sum.smanf,
	      byclass_sum.vbh/1000.,
	      byclass_sum.sbh/1000.,
	      byclass_sum.vbh2/1000.,
	      byclass_sum.sbh2/1000.,
	      byclass_sum.sbh3/1000.,
	      (byclass_sum.vbh+byclass_sum.sbh+byclass_sum.vbh2+byclass_sum.sbh2+byclass_sum.sbh3)/1000.,
              byclass_sum.wh,
	      byclass_sum.converted_forest_land_total,
	      byclass_sum.converted_forest_land,
	      byclass_sum.wh-((byclass_sum.vbh+byclass_sum.sbh+byclass_sum.vbh2+byclass_sum.sbh2+byclass_sum.sbh3)/1000.+byclass_sum.converted_forest_land),
	      byclass_sum.flowcp,
	      byclass_sum.flowpc,
	      byclass_sum.flowpv,
	      byclass_sum.flowvp,
	      byclass_sum.flowvc,
	      byclass_sum.flowcv,
	      byclass_sum.flowsp,
	      byclass_sum.flowps,
	      byclass_sum.flowsc,
	      byclass_sum.flowcs,
	      byclass_sum.flowvbh,
	      byclass_sum.flowsbh,
	      byclass_sum.flowvbh2,
	      byclass_sum.flowsbh2,
	      byclass_sum.flowsbh3,
	      land_garea_sum/1000.0/1000.0,
	      garea_sum/1000.0/1000.0,
	      total_garea_sum/1000.0/1000.0); 



      if(print_file_headers) fprintf(outfile2,"yr cs vc sc vp sp ps\n");
      
      fprintf(outfile2,"%d %lf %lf %lf %lf %lf %lf\n",
	      curryear,
	      byclass_sum.flowcs_prime,
	      byclass_sum.flowvc_prime,
	      byclass_sum.flowsc_prime,
	      byclass_sum.flowvp_prime,
	      byclass_sum.flowsp_prime,
	      byclass_sum.flowps_prime);
    } /* end of else urban */

    fclose(outfile);
    fclose(outfile2);


    return;

     
  }

/********************************************************************/

  void regional_timeseries_checker(int regional_code, char regional_name[50],int curryear){

  FILE *outfile;
  int i, ic, k, m, it;
  double garea_sum=ZEROVALUE, land_garea_sum=ZEROVALUE, total_garea_sum=ZEROVALUE;
  double nf_s_area=ZEROVALUE, fnf_s_area=ZEROVALUE;
  char ffname[90];
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
	    chtest[k][m][it].c = newdata[it].c[k][m];
	    chtest[k][m][it].p = newdata[it].p[k][m];
	    chtest[k][m][it].v = newdata[it].v[k][m];
	    chtest[k][m][it].s = newdata[it].s[k][m];
	    chtest[k][m][it].u = newdata[it].u[k][m];
	  }
	  
	  chtest[k][m][it+1].c=chtest[k][m][it].c +
	    newdata[it].flowpc[k][m] + newdata[it].flowvc[k][m] + 
	    newdata[it].flowsc[k][m] - newdata[it].flowcp[k][m] - 
	    newdata[it].flowcv[k][m] - newdata[it].flowcs[k][m] + 
	    newdata[it].flowuc[k][m] - newdata[it].flowcu[k][m]; 
	  
	  chtest[k][m][it+1].p=chtest[k][m][it].p +
	    newdata[it].flowcp[k][m] + newdata[it].flowvp[k][m] + 
	    newdata[it].flowsp[k][m] - newdata[it].flowpc[k][m] - 
	    newdata[it].flowpv[k][m] - newdata[it].flowps[k][m] +
	    newdata[it].flowup[k][m] - newdata[it].flowpu[k][m];
	    
	  chtest[k][m][it+1].v=chtest[k][m][it].v +
	    newdata[it].flowpv[k][m] + newdata[it].flowcv[k][m] -
	    newdata[it].flowvp[k][m] - newdata[it].flowvc[k][m] -   
	    newdata[it].flowvs[k][m] - newdata[it].flowvu[k][m];
	    
	  chtest[k][m][it+1].s=chtest[k][m][it].s +
	    newdata[it].flowvs[k][m] + newdata[it].flowps[k][m] + 
	    newdata[it].flowcs[k][m] -
	    newdata[it].flowsp[k][m] - newdata[it].flowsc[k][m] +  
	    newdata[it].flowus[k][m] - newdata[it].flowsu[k][m];

	  chtest[k][m][it+1].u=chtest[k][m][it].u +
	    newdata[it].flowcu[k][m] + newdata[it].flowvu[k][m] + 
	    newdata[it].flowsu[k][m] + newdata[it].flowpu[k][m] - 
	    newdata[it].flowuc[k][m] - newdata[it].flowup[k][m] -
	    newdata[it].flowus[k][m];  
	  
	 
	} /* end of country_test */
      } /* end of m */
    } /* end of m */
  }  /* end of it */
  
 

  for (it=0;it<NEWNT-1;it++){
    

    for (i=0;i<NREG;i++) {

      if(rdata[i].rcode == regional_code) {
	byclass_sum.wh=rtdata[i].wh;
	byclass_sum.converted_forest_land=rtdata[i].converted_forest_land;
	byclass_sum.converted_forest_land_total=rtdata[i].converted_forest_land_total;

      }
      
    } /* end of i */

 
    fnf_s_area=ZEROVALUE;
    nf_s_area=ZEROVALUE;
    
    

    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){		

	if((dstatic[k][m].rcode==regional_code) && (dstatic[k][m].usflag==1)){   	

	
	  if(dstatic[k][m].fnf == 1){
		
	    byclass_sum.vf+=garea[k][m]*chtest[k][m][it].v;
	    byclass_sum.sf+=garea[k][m]*chtest[k][m][it].s;
	    byclass_sum.smaf+=chtest[k][m][it].s*garea[k][m]*newdata[it].sma[k][m];
	    fnf_s_area += chtest[k][m][it].s*garea[k][m];
	  }
	  else{
	    byclass_sum.vnf+=garea[k][m]*chtest[k][m][it].v;
	    byclass_sum.snf+=garea[k][m]*chtest[k][m][it].s;
	    byclass_sum.smanf+=chtest[k][m][it].s*garea[k][m]*newdata[it].sma[k][m];
	    nf_s_area += chtest[k][m][it].s*garea[k][m];
	  }
	    
	  byclass_sum.vbh+=newdata[it].vbh[k][m]/1000.;
	  byclass_sum.sbh+=newdata[it].sbh[k][m]/1000.;
	  byclass_sum.vbh2+=newdata[it].vbh2[k][m]/1000.;
	  byclass_sum.sbh2+=newdata[it].sbh2[k][m]/1000.;
	  byclass_sum.sbh3+=newdata[it].sbh3[k][m]/1000.;
	  byclass_sum.wh_unmet+=ZEROVALUE; 
	  
	  
	  byclass_sum.c+=garea[k][m]*chtest[k][m][it].c; 
	  byclass_sum.p+=garea[k][m]*chtest[k][m][it].p;
	  byclass_sum.v+=garea[k][m]*chtest[k][m][it].v; 
	  byclass_sum.i+=garea[k][m]*newdata[it].i[k][m];  
	  byclass_sum.w+=garea[k][m]*newdata[it].w[k][m];  
	  byclass_sum.s+=garea[k][m]*chtest[k][m][it].s; 
	  byclass_sum.u+=garea[k][m]*chtest[k][m][it].u; 
	  
	  byclass_sum.flowcp+=garea[k][m]/1000.0/1000.0*newdata[it].flowcp[k][m];
	  byclass_sum.flowpc+=garea[k][m]/1000.0/1000.0*newdata[it].flowpc[k][m];
	  byclass_sum.flowpv+=garea[k][m]/1000.0/1000.0*newdata[it].flowpv[k][m];
	  byclass_sum.flowvp+=garea[k][m]/1000.0/1000.0*newdata[it].flowvp[k][m];
	  byclass_sum.flowvc+=garea[k][m]/1000.0/1000.0*newdata[it].flowvc[k][m];
	  byclass_sum.flowcv+=garea[k][m]/1000.0/1000.0*newdata[it].flowcv[k][m];
	  byclass_sum.flowsp+=garea[k][m]/1000.0/1000.0*newdata[it].flowsp[k][m];
	  byclass_sum.flowps+=garea[k][m]/1000.0/1000.0*newdata[it].flowps[k][m];
	  byclass_sum.flowsc+=garea[k][m]/1000.0/1000.0*newdata[it].flowsc[k][m];
	  byclass_sum.flowcs+=garea[k][m]/1000.0/1000.0*newdata[it].flowcs[k][m];
	  byclass_sum.flowcu+=garea[k][m]/1000.0/1000.0*newdata[it].flowcu[k][m];
	  byclass_sum.flowpu+=garea[k][m]/1000.0/1000.0*newdata[it].flowpu[k][m];
	  byclass_sum.flowvu+=garea[k][m]/1000.0/1000.0*newdata[it].flowvu[k][m];
	  byclass_sum.flowsu+=garea[k][m]/1000.0/1000.0*newdata[it].flowsu[k][m];
	  byclass_sum.flowuc+=garea[k][m]/1000.0/1000.0*newdata[it].flowuc[k][m];
	  byclass_sum.flowup+=garea[k][m]/1000.0/1000.0*newdata[it].flowup[k][m];
	  byclass_sum.flowus+=garea[k][m]/1000.0/1000.0*newdata[it].flowus[k][m];

          
	  
	  if(dstatic[k][m].fnf == 1){

	    if(dstatic[k][m].vba > ZEROVALUE) byclass_sum.flowvbh+=newdata[it].vbh[k][m]/dstatic[k][m].vba/1000.0/1000.0;

	    if(newdata[it].smb[k][m] > ZEROVALUE){
	     
	      if(newdata[it].smb[k][m] > ZEROVALUE) byclass_sum.flowsbh+=newdata[it].sbh[k][m]/newdata[it].smb[k][m]/1000.0/1000.0;
	      if(newdata[it].smb[k][m] > ZEROVALUE) byclass_sum.flowsbh2+=newdata[it].sbh2[k][m]/newdata[it].smb[k][m]/1000.0/1000.0;

	      if(newdata[it].smb[k][m] > ZEROVALUE) byclass_sum.flowsbh3+=newdata[it].sbh3[k][m]/newdata[it].smb[k][m]/1000.0/1000.0;



	    }
	    
	  }
	  else{
	    if(dstatic[k][m].vba > ZEROVALUE) byclass_sum.flowvbh2+=newdata[it].vbh2[k][m]/dstatic[k][m].vba/1000.0/1000.0;
	  } 
	  	    
	  
	  
	  if(it == 0 ) land_garea_sum+=(chtest[k][m][it].v+chtest[k][m][it].s+chtest[k][m][it].p+chtest[k][m][it].c+chtest[k][m][it].u)*garea[k][m]; 
	  
	  if(it == 0) garea_sum+=garea[k][m];
	      

	} /* end of country_test */

        if(it==0) total_garea_sum+=garea[k][m];


      } /* end of m*/
    } /* end of k*/
    
  
    byclass_sum.total+=byclass_sum.c+byclass_sum.p+
      byclass_sum.v+byclass_sum.i+
      byclass_sum.w+byclass_sum.s+byclass_sum.u;  
      
      
    if(fnf_s_area > ZEROVALUE) byclass_sum.smaf=byclass_sum.smaf/fnf_s_area;
    
    if(nf_s_area > ZEROVALUE) byclass_sum.smanf=byclass_sum.smanf/nf_s_area;
    
    if (use_urban) {

    if(print_file_headers) fprintf(outfile,"yr c p v i w s u fsum vf vnf sf snf smaf smanf vbh sbh vbh2 sbh2 sbh3 sum wh clearing_amount_total clearing_amount_used wh-sum-clearing_amount_ifused(unmet) flowcp flowpc flowpv flowvp flowvc flowcv flowsp flowps flowsc flowcs flowcu flowpu flowvu flowsu flowuc flowup flowus flowvbh flowsbh flowvbh2 flowsbh2 flowsbh3 land_asum_km2 land_iw_asum_km2 total_garea_sum_km2\n");
      
    fprintf(outfile,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	    curryear,
	    byclass_sum.c/garea_sum,
	    byclass_sum.p/garea_sum,
	    byclass_sum.v/garea_sum,
	    byclass_sum.i/garea_sum,
	    byclass_sum.w/garea_sum,
	    byclass_sum.s/garea_sum,
	    byclass_sum.u/garea_sum,
	    byclass_sum.total/garea_sum,
	    byclass_sum.vf/garea_sum,
	    byclass_sum.vnf/garea_sum,
	    byclass_sum.sf/garea_sum,
	    byclass_sum.snf/garea_sum,
	    byclass_sum.smaf,
	    byclass_sum.smanf,
	    byclass_sum.vbh,
	    byclass_sum.sbh,
	    byclass_sum.vbh2,
	    byclass_sum.sbh2,
	    byclass_sum.sbh3,
	    byclass_sum.vbh+byclass_sum.sbh+byclass_sum.vbh2+byclass_sum.sbh2+byclass_sum.sbh3,
	    byclass_sum.wh,
	    byclass_sum.converted_forest_land_total,
	    byclass_sum.converted_forest_land,
	    byclass_sum.wh-(byclass_sum.vbh+byclass_sum.sbh+byclass_sum.vbh2+byclass_sum.sbh2+byclass_sum.sbh3+byclass_sum.converted_forest_land),
	    byclass_sum.flowcp,
	    byclass_sum.flowpc,
   	    byclass_sum.flowpv,
	    byclass_sum.flowvp,
	    byclass_sum.flowvc,
	    byclass_sum.flowcv,
	    byclass_sum.flowsp,
	    byclass_sum.flowps,
	    byclass_sum.flowsc,
	    byclass_sum.flowcs,
	    byclass_sum.flowcu,
	    byclass_sum.flowpu,
	    byclass_sum.flowvu,
	    byclass_sum.flowsu,
	    byclass_sum.flowuc,
	    byclass_sum.flowup,
	    byclass_sum.flowus,
	    byclass_sum.flowvbh,
	    byclass_sum.flowsbh,
	    byclass_sum.flowvbh2,
	    byclass_sum.flowsbh2,
	    byclass_sum.flowsbh3,
	    land_garea_sum/1000.0/1000.0,
	    garea_sum/1000.0/1000.0,
	    total_garea_sum/1000.0/1000.0 ); 
    }else{

    if(print_file_headers) fprintf(outfile,"yr c p v i w s fsum vf vnf sf snf smaf smanf vbh sbh vbh2 sbh2 sbh3 sum wh clearing_amount_total clearing_amount_used wh-sum-clearing_amount_ifused(unmet) flowcp[k][m] flowpc flowpv flowvp flowvc flowcv flowsp flowps flowsc flowcs flowvbh flowsbh flowvbh2 flowsbh2 flowsbh3 land_asum_km2 land_iw_asum_km2 total_garea_sum_km2\n");
      
    fprintf(outfile,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	    curryear,
	    byclass_sum.c/garea_sum,
	    byclass_sum.p/garea_sum,
	    byclass_sum.v/garea_sum,
	    byclass_sum.i/garea_sum,
	    byclass_sum.w/garea_sum,
	    byclass_sum.s/garea_sum,
	    byclass_sum.total/garea_sum,
	    byclass_sum.vf/garea_sum,
	    byclass_sum.vnf/garea_sum,
	    byclass_sum.sf/garea_sum,
	    byclass_sum.snf/garea_sum,
	    byclass_sum.smaf,
	    byclass_sum.smanf,
	    byclass_sum.vbh,
	    byclass_sum.sbh,
	    byclass_sum.vbh2,
	    byclass_sum.sbh2,
	    byclass_sum.sbh3,
	    byclass_sum.vbh+byclass_sum.sbh+byclass_sum.vbh2+byclass_sum.sbh2+byclass_sum.sbh3,
	    byclass_sum.wh,
	    byclass_sum.converted_forest_land_total,
	    byclass_sum.converted_forest_land,
	    byclass_sum.wh-(byclass_sum.vbh+byclass_sum.sbh+byclass_sum.vbh2+byclass_sum.sbh2+byclass_sum.sbh3+byclass_sum.converted_forest_land),
	    byclass_sum.flowcp,
	    byclass_sum.flowpc,
   	    byclass_sum.flowpv,
	    byclass_sum.flowvp,
	    byclass_sum.flowvc,
	    byclass_sum.flowcv,
	    byclass_sum.flowsp,
	    byclass_sum.flowps,
	    byclass_sum.flowsc,
	    byclass_sum.flowcs,
	    byclass_sum.flowvbh,
	    byclass_sum.flowsbh,
	    byclass_sum.flowvbh2,
	    byclass_sum.flowsbh2,
	    byclass_sum.flowsbh3,
	    land_garea_sum/1000.0/1000.0,
	    garea_sum/1000.0/1000.0,
	    total_garea_sum/1000.0/1000.0 ); 
    }      
  } /* end of it */

  fclose(outfile);

  

  return;
  
}


/*****************************************************************/
  void continent_timeseries_checker(int continent_code, char continent_name[50], int curryear){

  FILE *outfile, *outfile2;
  int i, ic, k, m, it;
  double garea_sum=ZEROVALUE, land_garea_sum=ZEROVALUE, total_garea_sum=ZEROVALUE;
  double nf_s_area=ZEROVALUE, fnf_s_area=ZEROVALUE;
  char ffname[90], ffname2[90];
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
	
	
	if(dstatic[k][m].continent_code==continent_code){   	
	  
	  if(it == 0) {
	    chtest[k][m][it].c = newdata[it].c[k][m];
	    chtest[k][m][it].p = newdata[it].p[k][m];
	    chtest[k][m][it].v = newdata[it].v[k][m];
	    chtest[k][m][it].s = newdata[it].s[k][m];
	    chtest[k][m][it].u = newdata[it].u[k][m];
	  }
	  
	  chtest[k][m][it+1].c=chtest[k][m][it].c +
	    newdata[it].flowpc[k][m] + newdata[it].flowvc[k][m] + 
	    newdata[it].flowsc[k][m] - newdata[it].flowcp[k][m] - 
	    newdata[it].flowcv[k][m] - newdata[it].flowcs[k][m] + 
	    newdata[it].flowuc[k][m] - newdata[it].flowcu[k][m]; 
	  
	  chtest[k][m][it+1].p=chtest[k][m][it].p +
	    newdata[it].flowcp[k][m] + newdata[it].flowvp[k][m] + 
	    newdata[it].flowsp[k][m] - newdata[it].flowpc[k][m] - 
	    newdata[it].flowpv[k][m] - newdata[it].flowps[k][m] +
	    newdata[it].flowup[k][m] - newdata[it].flowpu[k][m];
	    
	  chtest[k][m][it+1].v=chtest[k][m][it].v +
	    newdata[it].flowpv[k][m] + newdata[it].flowcv[k][m] -
	    newdata[it].flowvp[k][m] - newdata[it].flowvc[k][m] -   
	    newdata[it].flowvs[k][m] - newdata[it].flowvu[k][m];
	    
	  chtest[k][m][it+1].s=chtest[k][m][it].s +
	    newdata[it].flowvs[k][m] + newdata[it].flowps[k][m] + 
	    newdata[it].flowcs[k][m] -
	    newdata[it].flowsp[k][m] - newdata[it].flowsc[k][m] +  
            newdata[it].flowus[k][m] - newdata[it].flowsu[k][m];

	  chtest[k][m][it+1].u=chtest[k][m][it].u +
	    newdata[it].flowcu[k][m] + newdata[it].flowpu[k][m] + 
	    newdata[it].flowsu[k][m] + newdata[it].flowvu[k][m] - 
	    newdata[it].flowuc[k][m] - newdata[it].flowup[k][m] -
	    newdata[it].flowus[k][m];  
	  
	} /* end of continent_test */
      } /* end of m */
    } /* end of m */
  }  /* end of it */
  

  for (it=0;it<NEWNT-1;it++){
    

    for (i=0;i<NREG;i++) {

      if(rdata[i].continent_code == continent_code) {
	byclass_sum.wh+=rtdata[i].wh;
	byclass_sum.converted_forest_land+=rtdata[i].converted_forest_land;
	byclass_sum.converted_forest_land_total+=rtdata[i].converted_forest_land_total;

	byclass_sum.flowcs_prime+=rtdata[i].flowcs_prime;
	byclass_sum.flowvc_prime+=rtdata[i].flowvc_prime;
	byclass_sum.flowsc_prime+=rtdata[i].flowsc_prime;
	byclass_sum.flowvp_prime+=rtdata[i].flowvp_prime;
	byclass_sum.flowsp_prime+=rtdata[i].flowsp_prime;
	byclass_sum.flowps_prime+=rtdata[i].flowps_prime;
	byclass_sum.flowus_prime+=rtdata[i].flowus_prime;
	byclass_sum.flowsu_prime+=rtdata[i].flowsu_prime;
	byclass_sum.flowvu_prime+=rtdata[i].flowvu_prime;
      }
      
    } /* end of i */

 
    fnf_s_area=ZEROVALUE;
    nf_s_area=ZEROVALUE;
    
    

    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){		

	if(dstatic[k][m].continent_code==continent_code){   	

	
	  if(dstatic[k][m].fnf == 1){
		
	    byclass_sum.vf+=garea[k][m]*chtest[k][m][it].v;
	    byclass_sum.sf+=garea[k][m]*chtest[k][m][it].s;
	    byclass_sum.smaf+=chtest[k][m][it].s*garea[k][m]*newdata[it].sma[k][m];
	    fnf_s_area += chtest[k][m][it].s*garea[k][m];
	  }
	  else{
	    byclass_sum.vnf+=garea[k][m]*chtest[k][m][it].v;
	    byclass_sum.snf+=garea[k][m]*chtest[k][m][it].s;
	    byclass_sum.smanf+=chtest[k][m][it].s*garea[k][m]*newdata[it].sma[k][m];
	    nf_s_area += chtest[k][m][it].s*garea[k][m];
	  }
	    
	  byclass_sum.vbh+=newdata[it].vbh[k][m]/1000.;
	  byclass_sum.sbh+=newdata[it].sbh[k][m]/1000.;
	  byclass_sum.vbh2+=newdata[it].vbh2[k][m]/1000.;
	  byclass_sum.sbh2+=newdata[it].sbh2[k][m]/1000.;
	  byclass_sum.sbh3+=newdata[it].sbh3[k][m]/1000.;
	  byclass_sum.wh_unmet+=ZEROVALUE; 
	  
	  
	  byclass_sum.c+=garea[k][m]*chtest[k][m][it].c; 
	  byclass_sum.p+=garea[k][m]*chtest[k][m][it].p;
	  byclass_sum.v+=garea[k][m]*chtest[k][m][it].v; 
	  byclass_sum.i+=garea[k][m]*newdata[it].i[k][m];  
	  byclass_sum.w+=garea[k][m]*newdata[it].w[k][m];  
	  byclass_sum.s+=garea[k][m]*chtest[k][m][it].s; 
	  byclass_sum.u+=garea[k][m]*chtest[k][m][it].u; 
	  
	  byclass_sum.flowcp+=garea[k][m]/1000.0/1000.0*newdata[it].flowcp[k][m];
	  byclass_sum.flowpc+=garea[k][m]/1000.0/1000.0*newdata[it].flowpc[k][m];
	  byclass_sum.flowpv+=garea[k][m]/1000.0/1000.0*newdata[it].flowpv[k][m];
	  byclass_sum.flowvp+=garea[k][m]/1000.0/1000.0*newdata[it].flowvp[k][m];
	  byclass_sum.flowvc+=garea[k][m]/1000.0/1000.0*newdata[it].flowvc[k][m];
	  byclass_sum.flowcv+=garea[k][m]/1000.0/1000.0*newdata[it].flowcv[k][m];
	  byclass_sum.flowsp+=garea[k][m]/1000.0/1000.0*newdata[it].flowsp[k][m];
	  byclass_sum.flowps+=garea[k][m]/1000.0/1000.0*newdata[it].flowps[k][m];
	  byclass_sum.flowsc+=garea[k][m]/1000.0/1000.0*newdata[it].flowsc[k][m];
	  byclass_sum.flowcs+=garea[k][m]/1000.0/1000.0*newdata[it].flowcs[k][m];
	  byclass_sum.flowcu+=garea[k][m]/1000.0/1000.0*newdata[it].flowcu[k][m];
	  byclass_sum.flowpu+=garea[k][m]/1000.0/1000.0*newdata[it].flowpu[k][m];
	  byclass_sum.flowvu+=garea[k][m]/1000.0/1000.0*newdata[it].flowvu[k][m];
	  byclass_sum.flowsu+=garea[k][m]/1000.0/1000.0*newdata[it].flowsu[k][m];
	  byclass_sum.flowuc+=garea[k][m]/1000.0/1000.0*newdata[it].flowuc[k][m];
	  byclass_sum.flowup+=garea[k][m]/1000.0/1000.0*newdata[it].flowup[k][m];
	  byclass_sum.flowus+=garea[k][m]/1000.0/1000.0*newdata[it].flowus[k][m];

          
	  
	  if(dstatic[k][m].fnf == 1){

	    byclass_sum.flowvbh+=newdata[it].vbh[k][m]/dstatic[k][m].vba/1000.0/1000.0;

	    if(newdata[it].smb[k][m] > ZEROVALUE){
	     
	      byclass_sum.flowsbh+=newdata[it].sbh[k][m]/newdata[it].smb[k][m]/1000.0/1000.0;
	      byclass_sum.flowsbh2+=newdata[it].sbh2[k][m]/newdata[it].smb[k][m]/1000.0/1000.0;
	      byclass_sum.flowsbh3+=newdata[it].sbh3[k][m]/newdata[it].smb[k][m]/1000.0/1000.0;
	    }
	    
	  }
	  else{
	    if(dstatic[k][m].vba > ZEROVALUE) byclass_sum.flowvbh2+=newdata[it].vbh2[k][m]/dstatic[k][m].vba/1000.0/1000.0;
	  } 
	  	    
	  
	  if(it == 0 ) land_garea_sum+=(chtest[k][m][it].v+chtest[k][m][it].s+chtest[k][m][it].p+chtest[k][m][it].c+chtest[k][m][it].u)*garea[k][m]; 
	  
	  if(it == 0) garea_sum+=garea[k][m];
	      

	} /* end of continent_test */

        if(it==0) total_garea_sum+=garea[k][m];


      } /* end of m*/
    } /* end of k*/
    
    
    byclass_sum.total+=byclass_sum.c+byclass_sum.p+
      byclass_sum.v+byclass_sum.i+
      byclass_sum.w+byclass_sum.s+byclass_sum.u;  
      
      
    if(fnf_s_area > ZEROVALUE) byclass_sum.smaf=byclass_sum.smaf/fnf_s_area;
    
    if(nf_s_area > ZEROVALUE) byclass_sum.smanf=byclass_sum.smanf/nf_s_area;    
    
    if (use_urban) {

    if(print_file_headers) fprintf(outfile,"yr c p v i w s u fsum vf vnf sf snf smaf smanf vbh sbh vbh2 sbh2 sbh3 sum wh clearing_amount_total clearing_amount_used wh-sum-clearing_amount_ifused(unmet) flowcp flowpc flowpv flowvp flowvc flowcv flowsp flowps flowsc flowcs flowcu flowpu flowvu flowsu flowuc flowup flowus flowvbh flowsbh flowvbh2 flowsbh2 flowsbh3 land_asum_km2 land_iw_asum_km2 total_garea_sum_km2\n");
      
     fprintf(outfile,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	    curryear,
	    byclass_sum.c/garea_sum,
	    byclass_sum.p/garea_sum,
	    byclass_sum.v/garea_sum,
	    byclass_sum.i/garea_sum,
	    byclass_sum.w/garea_sum,
	    byclass_sum.s/garea_sum,
	    byclass_sum.u/garea_sum,
	    byclass_sum.total/garea_sum,
	    byclass_sum.vf/garea_sum,
	    byclass_sum.vnf/garea_sum,
	    byclass_sum.sf/garea_sum,
	    byclass_sum.snf/garea_sum,
	    byclass_sum.smaf,
	    byclass_sum.smanf,
	    byclass_sum.vbh,
	    byclass_sum.sbh,
	    byclass_sum.vbh2,
	    byclass_sum.sbh2,
	     byclass_sum.sbh3,
	    byclass_sum.vbh+byclass_sum.sbh+byclass_sum.vbh2+byclass_sum.sbh2+byclass_sum.sbh3,
	     byclass_sum.wh,
	     byclass_sum.converted_forest_land_total,
	     byclass_sum.converted_forest_land,
	    byclass_sum.wh-(byclass_sum.vbh+byclass_sum.sbh+byclass_sum.vbh2+byclass_sum.sbh2+byclass_sum.sbh3+byclass_sum.converted_forest_land),
	     byclass_sum.flowcp,
	    byclass_sum.flowpc,
   	    byclass_sum.flowpv,
	    byclass_sum.flowvp,
	    byclass_sum.flowvc,
	    byclass_sum.flowcv,
	    byclass_sum.flowsp,
	    byclass_sum.flowps,
	    byclass_sum.flowsc,
	    byclass_sum.flowcs,
	    byclass_sum.flowcu,
	    byclass_sum.flowpu,
	    byclass_sum.flowvu,
	    byclass_sum.flowsu,
	    byclass_sum.flowuc,
	    byclass_sum.flowup,
	    byclass_sum.flowus,
	    byclass_sum.flowvbh,
	    byclass_sum.flowsbh,
	    byclass_sum.flowvbh2,
	    byclass_sum.flowsbh2,
	     byclass_sum.flowsbh3,
	    land_garea_sum/1000.0/1000.0,
	     garea_sum/1000.0/1000.0,
	     total_garea_sum/1000.0/1000.0 ); 


     if(print_file_headers) fprintf(outfile2,"yr cs vc sc vp sp ps su us vu\n");
      
     fprintf(outfile2,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	     curryear,
	     byclass_sum.flowcs_prime,
	     byclass_sum.flowvc_prime,
	     byclass_sum.flowsc_prime,
	     byclass_sum.flowvp_prime,
	     byclass_sum.flowsp_prime,
	     byclass_sum.flowps_prime,
	     byclass_sum.flowsu_prime,
	     byclass_sum.flowus_prime,
	     byclass_sum.flowvu_prime);   


    }else{

    if(print_file_headers) fprintf(outfile,"yr c p v i w s fsum vf vnf sf snf smaf smanf vbh sbh vbh2 sbh2 sbh3 sum wh clearing_amount_total clearing_amount_used wh-sum-clearing_amount_ifused(unmet) flowcp[k][m] flowpc flowpv flowvp flowvc flowcv flowsp flowps flowsc flowcs flowvbh flowsbh flowvbh2 flowsbh2 flowsbh3 land_asum_km2 land_iw_asum_km2 total_garea_sum_km2\n");
      
     fprintf(outfile,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	    curryear,
	    byclass_sum.c/garea_sum,
	    byclass_sum.p/garea_sum,
	    byclass_sum.v/garea_sum,
	    byclass_sum.i/garea_sum,
	    byclass_sum.w/garea_sum,
	    byclass_sum.s/garea_sum,
	    byclass_sum.total/garea_sum,
	    byclass_sum.vf/garea_sum,
	    byclass_sum.vnf/garea_sum,
	    byclass_sum.sf/garea_sum,
	    byclass_sum.snf/garea_sum,
	    byclass_sum.smaf,
	    byclass_sum.smanf,
	    byclass_sum.vbh,
	    byclass_sum.sbh,
	    byclass_sum.vbh2,
	    byclass_sum.sbh2,
	     byclass_sum.sbh3,
	    byclass_sum.vbh+byclass_sum.sbh+byclass_sum.vbh2+byclass_sum.sbh2+byclass_sum.sbh3,
	     byclass_sum.wh,
	     byclass_sum.converted_forest_land_total,
	     byclass_sum.converted_forest_land,
	    byclass_sum.wh-(byclass_sum.vbh+byclass_sum.sbh+byclass_sum.vbh2+byclass_sum.sbh2+byclass_sum.sbh3+byclass_sum.converted_forest_land),
	     byclass_sum.flowcp,
	    byclass_sum.flowpc,
   	    byclass_sum.flowpv,
	    byclass_sum.flowvp,
	    byclass_sum.flowvc,
	    byclass_sum.flowcv,
	    byclass_sum.flowsp,
	    byclass_sum.flowps,
	    byclass_sum.flowsc,
	    byclass_sum.flowcs,
	    byclass_sum.flowvbh,
	    byclass_sum.flowsbh,
	    byclass_sum.flowvbh2,
	    byclass_sum.flowsbh2,
	     byclass_sum.flowsbh3,
	    land_garea_sum/1000.0/1000.0,
	     garea_sum/1000.0/1000.0,
	     total_garea_sum/1000.0/1000.0 ); 


     if((it == 0) && (strcmp("tone",trun) == 0)) fprintf(outfile2,"yr cs vc sc vp sp ps\n");
      
     fprintf(outfile2,"%d %lf %lf %lf %lf %lf %lf\n",
	     curryear,
	     byclass_sum.flowcs_prime,
	     byclass_sum.flowvc_prime,
	     byclass_sum.flowsc_prime,
	     byclass_sum.flowvp_prime,
	     byclass_sum.flowsp_prime,
	     byclass_sum.flowps_prime);   

    } 
  } /* end of it */

  fclose(outfile);
  fclose(outfile2);


  return;
  
}
/********************************************************************/

void loop_call_for_country_final_stats(int curryear){

  int i, country_code;


  for (i=0;i<NCCODE;i++){

     country_code=cdata[i].ccode;
     country_final_stats(country_code,i,curryear);


  } /* end of i */

  return;
}

/********************************************************************/

 void country_final_stats(int country_code, int istart, int curryear){

  FILE *outfile;
  int i, ic, k, m, it;
  double garea_sum=ZEROVALUE, land_garea_sum=ZEROVALUE, total_garea_sum=ZEROVALUE;
  double nf_s_area=ZEROVALUE, fnf_s_area=ZEROVALUE;
  char ffname[90];
  double tester1, tester2, tester3, tester4;



  initialize_checker(country_code);


  if((strcmp("tone",trun) == 0) && (istart == 0))
    outfile=fopen("country.final.stats.txt","w");
  else
    outfile=fopen("country.final.stats.txt","a");  

  it=0;

    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){   
	
	
	if((dstatic[k][m].gcode==country_code) && (dstatic[k][m].usflag==1)){   	
	  
	  if(it == 0) {
	    chtest[k][m][it].c = newdata[it].c[k][m];
	    chtest[k][m][it].p = newdata[it].p[k][m];
	    chtest[k][m][it].v = newdata[it].v[k][m];
	    chtest[k][m][it].s = newdata[it].s[k][m];
	    chtest[k][m][it].u = newdata[it].u[k][m];
	  }
	  
	  chtest[k][m][it+1].c=chtest[k][m][it].c +
	    newdata[it].flowpc[k][m] + newdata[it].flowvc[k][m] + 
	    newdata[it].flowsc[k][m] - newdata[it].flowcp[k][m] - 
	    newdata[it].flowcv[k][m] - newdata[it].flowcs[k][m] + 
	    newdata[it].flowuc[k][m] - newdata[it].flowcu[k][m]; 
	  
	  chtest[k][m][it+1].p=chtest[k][m][it].p +
	    newdata[it].flowcp[k][m] + newdata[it].flowvp[k][m] + 
	    newdata[it].flowsp[k][m] - newdata[it].flowpc[k][m] - 
	    newdata[it].flowpv[k][m] - newdata[it].flowps[k][m] +
	    newdata[it].flowup[k][m] - newdata[it].flowpu[k][m];
	    
	  chtest[k][m][it+1].v=chtest[k][m][it].v +
	    newdata[it].flowpv[k][m] + newdata[it].flowcv[k][m] -
	    newdata[it].flowvp[k][m] - newdata[it].flowvc[k][m] -   
	    newdata[it].flowvs[k][m] - newdata[it].flowvu[k][m];
	    
	  chtest[k][m][it+1].s=chtest[k][m][it].s +
	    newdata[it].flowvs[k][m] + newdata[it].flowps[k][m] + 
	    newdata[it].flowcs[k][m] -
	    newdata[it].flowsp[k][m] - newdata[it].flowsc[k][m] +
	    newdata[it].flowus[k][m] - newdata[it].flowsu[k][m];

	  chtest[k][m][it+1].u=chtest[k][m][it].u +
	    newdata[it].flowcu[k][m] + newdata[it].flowpu[k][m] + 
	    newdata[it].flowvu[k][m] + newdata[it].flowsu[k][m] - 
	    newdata[it].flowuc[k][m] - newdata[it].flowup[k][m] -
	    newdata[it].flowus[k][m];   
	  
	 
	} /* end of country_test */
      } /* end of m */
    } /* end of m */

    it=0;
    for (i=0;i<NCCODE;i++) {

      if(cdata[i].ccode == country_code) {
	byclass_sum.wh=rtdata[i].wh;
	byclass_sum.converted_forest_land=rtdata[i].converted_forest_land;
	byclass_sum.converted_forest_land_total=rtdata[i].converted_forest_land_total;

      }
      
    } /* end of i */

 
    fnf_s_area=ZEROVALUE;
    nf_s_area=ZEROVALUE;
    
    

    for (k=0;k<NY;k++){
      for (m=0;m<NX;m++){		

	if((dstatic[k][m].gcode==country_code) && (dstatic[k][m].usflag==1)){   	

	
	  if(dstatic[k][m].fnf == 1){
		
	    byclass_sum.vf+=garea[k][m]*chtest[k][m][it].v;
	    byclass_sum.sf+=garea[k][m]*chtest[k][m][it].s;
	    byclass_sum.smaf+=chtest[k][m][it].s*garea[k][m]*newdata[it].sma[k][m];
	    fnf_s_area += chtest[k][m][it].s*garea[k][m];
	  }
	  else{
	    byclass_sum.vnf+=garea[k][m]*chtest[k][m][it].v;
	    byclass_sum.snf+=garea[k][m]*chtest[k][m][it].s;
	    byclass_sum.smanf+=chtest[k][m][it].s*garea[k][m]*newdata[it].sma[k][m];
	    nf_s_area += chtest[k][m][it].s*garea[k][m];
	  }
	    
	  byclass_sum.vbh+=newdata[it].vbh[k][m]/1000.;
	  byclass_sum.sbh+=newdata[it].sbh[k][m]/1000.;
	  byclass_sum.vbh2+=newdata[it].vbh2[k][m]/1000.;
	  byclass_sum.sbh2+=newdata[it].sbh2[k][m]/1000.;
	  byclass_sum.sbh3+=newdata[it].sbh3[k][m]/1000.;
	  byclass_sum.wh_unmet+=ZEROVALUE;  
	  
	  
	  byclass_sum.c+=garea[k][m]*chtest[k][m][it].c; 
	  byclass_sum.p+=garea[k][m]*chtest[k][m][it].p;
	  byclass_sum.v+=garea[k][m]*chtest[k][m][it].v; 
	  byclass_sum.i+=garea[k][m]*newdata[it].i[k][m];  
	  byclass_sum.w+=garea[k][m]*newdata[it].w[k][m];  
	  byclass_sum.s+=garea[k][m]*chtest[k][m][it].s; 
	  byclass_sum.u+=garea[k][m]*chtest[k][m][it].u; 
	  
	  byclass_sum.flowcp+=garea[k][m]/1000.0/1000.0*newdata[it].flowcp[k][m];
	  byclass_sum.flowpc+=garea[k][m]/1000.0/1000.0*newdata[it].flowpc[k][m];
	  byclass_sum.flowpv+=garea[k][m]/1000.0/1000.0*newdata[it].flowpv[k][m];
	  byclass_sum.flowvp+=garea[k][m]/1000.0/1000.0*newdata[it].flowvp[k][m];
	  byclass_sum.flowvc+=garea[k][m]/1000.0/1000.0*newdata[it].flowvc[k][m];
	  byclass_sum.flowcv+=garea[k][m]/1000.0/1000.0*newdata[it].flowcv[k][m];
	  byclass_sum.flowsp+=garea[k][m]/1000.0/1000.0*newdata[it].flowsp[k][m];
	  byclass_sum.flowps+=garea[k][m]/1000.0/1000.0*newdata[it].flowps[k][m];
	  byclass_sum.flowsc+=garea[k][m]/1000.0/1000.0*newdata[it].flowsc[k][m];
	  byclass_sum.flowcs+=garea[k][m]/1000.0/1000.0*newdata[it].flowcs[k][m];
	  byclass_sum.flowcu+=garea[k][m]/1000.0/1000.0*newdata[it].flowcu[k][m];
	  byclass_sum.flowpu+=garea[k][m]/1000.0/1000.0*newdata[it].flowpu[k][m];
	  byclass_sum.flowvu+=garea[k][m]/1000.0/1000.0*newdata[it].flowvu[k][m];
	  byclass_sum.flowsu+=garea[k][m]/1000.0/1000.0*newdata[it].flowsu[k][m];
	  byclass_sum.flowuc+=garea[k][m]/1000.0/1000.0*newdata[it].flowuc[k][m];
	  byclass_sum.flowup+=garea[k][m]/1000.0/1000.0*newdata[it].flowup[k][m];
	  byclass_sum.flowus+=garea[k][m]/1000.0/1000.0*newdata[it].flowus[k][m];

          
	  
	  if(dstatic[k][m].fnf == 1){

	    if(dstatic[k][m].vba > ZEROVALUE) byclass_sum.flowvbh+=newdata[it].vbh[k][m]/dstatic[k][m].vba/1000.0/1000.0;

	    if(newdata[it].smb[k][m] > ZEROVALUE){
	     
	      if(newdata[it].smb[k][m] > ZEROVALUE) byclass_sum.flowsbh+=newdata[it].sbh[k][m]/newdata[it].smb[k][m]/1000.0/1000.0;
	      if(newdata[it].smb[k][m] > ZEROVALUE) byclass_sum.flowsbh2+=newdata[it].sbh2[k][m]/newdata[it].smb[k][m]/1000.0/1000.0;
	      if(newdata[it].smb[k][m] > ZEROVALUE) byclass_sum.flowsbh3+=newdata[it].sbh3[k][m]/newdata[it].smb[k][m]/1000.0/1000.0;
	    }
	    
	  }
	  else{
	    if(dstatic[k][m].vba > ZEROVALUE) byclass_sum.flowvbh2+=newdata[it].vbh2[k][m]/dstatic[k][m].vba/1000.0/1000.0;
	  } 
	  	    
	  
	  
	  if(it == 0 ) land_garea_sum+=(chtest[k][m][it].v+chtest[k][m][it].s+chtest[k][m][it].p+chtest[k][m][it].c+chtest[k][m][it].u)*garea[k][m]; 
	  
	  if(it == 0) garea_sum+=garea[k][m];
	      

	} /* end of country_test */

        if(it==0) total_garea_sum+=garea[k][m];


      } /* end of m*/
    } /* end of k*/
    
  
    byclass_sum.total+=byclass_sum.c+byclass_sum.p+
      byclass_sum.v+byclass_sum.i+
      byclass_sum.w+byclass_sum.s+byclass_sum.u;  
      
      
    if(fnf_s_area > ZEROVALUE) byclass_sum.smaf=byclass_sum.smaf/fnf_s_area;
    
    if(nf_s_area > ZEROVALUE) byclass_sum.smanf=byclass_sum.smanf/nf_s_area;
    

    if (use_urban) {

      if(print_file_headers) fprintf(outfile,"yr cname c p v i w s u fsum vf vnf sf snf smaf smanf vbh sbh vbh2 sbh2 sbh3 sum wh clearing_amount_total clearing_amount_used wh-sum-clearing_amount_ifused(unmet) flowcp flowpc flowpv flowvp flowvc flowcv flowsp flowps flowsc flowcs flowcu flowpu flowvu flowsu flowuc flowup flowus flowvbh flowsbh flowvbh2 flowsbh3 flowsbh2 land_asum_km2 land_iw_asum_km2 total_garea_sum_km2\n");
      
    fprintf(outfile,"%d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	    curryear,cname[istart],
	    byclass_sum.c/garea_sum,
	    byclass_sum.p/garea_sum,
	    byclass_sum.v/garea_sum,
	    byclass_sum.i/garea_sum,
	    byclass_sum.w/garea_sum,
	    byclass_sum.s/garea_sum,
	    byclass_sum.u/garea_sum,
	    byclass_sum.total/garea_sum,
	    byclass_sum.vf/garea_sum,
	    byclass_sum.vnf/garea_sum,
	    byclass_sum.sf/garea_sum,
	    byclass_sum.snf/garea_sum,
	    byclass_sum.smaf,
	    byclass_sum.smanf,
	    byclass_sum.vbh,
	    byclass_sum.sbh,
	    byclass_sum.vbh2,
	    byclass_sum.sbh2,
	    byclass_sum.sbh3,
	    byclass_sum.vbh+byclass_sum.sbh+byclass_sum.vbh2+byclass_sum.sbh2+byclass_sum.sbh3,
	    byclass_sum.wh,
	    byclass_sum.converted_forest_land_total,
	    byclass_sum.converted_forest_land,
	    byclass_sum.wh-(byclass_sum.vbh+byclass_sum.sbh+byclass_sum.vbh2+byclass_sum.sbh2+byclass_sum.sbh3+byclass_sum.converted_forest_land),
	    byclass_sum.flowcp,
	    byclass_sum.flowpc,
   	    byclass_sum.flowpv,
	    byclass_sum.flowvp,
	    byclass_sum.flowvc,
	    byclass_sum.flowcv,
	    byclass_sum.flowsp,
	    byclass_sum.flowps,
	    byclass_sum.flowsc,
	    byclass_sum.flowcs,
	    byclass_sum.flowcu,
	    byclass_sum.flowpu,
	    byclass_sum.flowvu,
	    byclass_sum.flowsu,
	    byclass_sum.flowuc,
	    byclass_sum.flowup,
	    byclass_sum.flowus,
	    byclass_sum.flowvbh,
	    byclass_sum.flowsbh,
	    byclass_sum.flowvbh2,
	    byclass_sum.flowsbh2,
	    byclass_sum.flowsbh3,
	    land_garea_sum/1000.0/1000.0,
	    garea_sum/1000.0/1000.0,
	    total_garea_sum/1000.0/1000.0 ); 

    }else{

      if(print_file_headers) fprintf(outfile,"yr cname c p v i w s fsum vf vnf sf snf smaf smanf vbh sbh vbh2 sbh2 sbh3 sum wh clearing_amount_total clearing_amount_used wh-sum-clearing_amount_ifused(unmet) flowcp[k][m] flowpc flowpv flowvp flowvc flowcv flowsp flowps flowsc flowcs flowvbh flowsbh flowvbh2 flowsbh3 flowsbh2 land_asum_km2 land_iw_asum_km2 total_garea_sum_km2\n");
      
    fprintf(outfile,"%d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	    curryear,cname[istart],
	    byclass_sum.c/garea_sum,
	    byclass_sum.p/garea_sum,
	    byclass_sum.v/garea_sum,
	    byclass_sum.i/garea_sum,
	    byclass_sum.w/garea_sum,
	    byclass_sum.s/garea_sum,
	    byclass_sum.total/garea_sum,
	    byclass_sum.vf/garea_sum,
	    byclass_sum.vnf/garea_sum,
	    byclass_sum.sf/garea_sum,
	    byclass_sum.snf/garea_sum,
	    byclass_sum.smaf,
	    byclass_sum.smanf,
	    byclass_sum.vbh,
	    byclass_sum.sbh,
	    byclass_sum.vbh2,
	    byclass_sum.sbh2,
	    byclass_sum.sbh3,
	    byclass_sum.vbh+byclass_sum.sbh+byclass_sum.vbh2+byclass_sum.sbh2+byclass_sum.sbh3,
	    byclass_sum.wh,
	    byclass_sum.converted_forest_land_total,
	    byclass_sum.converted_forest_land,
	    byclass_sum.wh-(byclass_sum.vbh+byclass_sum.sbh+byclass_sum.vbh2+byclass_sum.sbh2+byclass_sum.sbh3+byclass_sum.converted_forest_land),
	    byclass_sum.flowcp,
	    byclass_sum.flowpc,
   	    byclass_sum.flowpv,
	    byclass_sum.flowvp,
	    byclass_sum.flowvc,
	    byclass_sum.flowcv,
	    byclass_sum.flowsp,
	    byclass_sum.flowps,
	    byclass_sum.flowsc,
	    byclass_sum.flowcs,
	    byclass_sum.flowvbh,
	    byclass_sum.flowsbh,
	    byclass_sum.flowvbh2,
	    byclass_sum.flowsbh2,
	    byclass_sum.flowsbh3,
	    land_garea_sum/1000.0/1000.0,
	    garea_sum/1000.0/1000.0,
	    total_garea_sum/1000.0/1000.0 ); 
    }      

  fclose(outfile);

  

  return;
  
}

/********************************************************************/
 void predict_available_biomass(int it,int i){
  int k, m, iz;
  double vf=ZEROVALUE, vnf=ZEROVALUE, sf=ZEROVALUE, snf=ZEROVALUE;

  
  for (iz=0;iz<NZ;iz++){
    vf+=cztdata[i][iz].vb;
  }


  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){


#if COUNTRY_TEST

      if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){

#else
	if(dstatic[k][m].gcode > 0){
#endif


	  if(dstatic[k][m].gcode == cdata[i].ccode){


	    if(dstatic[k][m].fnf == 1){

	      sf+=newdata[it].smb[k][m]*(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*garea[k][m]/1000.;
	      
	    } 
	    else{

	      snf+=newdata[it].smb[k][m]*(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*garea[k][m]/1000.;
		vnf+=dstatic[k][m].vba*(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*garea[k][m]/1000.;

	    } /* end of fnf */
	  
           }  /* end of gcode */
	}
	
      } /* end of m */
    } /* end of k */

    ctdata[i].predict_b = vf+vnf+sf+snf;

  return;

}

/********************************************************************/
  void harvest_gridded_data(int it){
    int k, m, i;
    double vf=ZEROVALUE, vnf=ZEROVALUE, sf=ZEROVALUE, snf=ZEROVALUE;
    double virgin_biomass, secondary_biomass, whr;


  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){

#if COUNTRY_TEST

      if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){

#else
	if(dstatic[k][m].gcode > 0){
#endif
	  i=dstatic[k][m].newgcode;

          virgin_biomass = dstatic[k][m].vba*(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*garea[k][m]/1000.;
          secondary_biomass = newdata[it].smb[k][m]*(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*garea[k][m]/1000.;
	  
          whr = newdata[it].whd[k][m];

	  if (BEST_CASE){

	if((cdata[i].continent_code == 3) || (cdata[i].continent_code == 4)){

	  /* if Eurasia, for h1 do not count clearing in harvest, but track its total */
          /* if Eurasia, for h3 count clearing in harvest and track its total */

	  if((nodata_option == 5)||(nodata_option == 6)){

	    newdata[it].converted_forest_land_total[k][m]=newdata[it].converted_forest_land[k][m];
	    newdata[it].converted_forest_land[k][m]=ZEROVALUE;
	  }

	}
	
        else{

	  /* if not Eurasia, do not count clearing in harvest, but track its total */

	    newdata[it].converted_forest_land_total[k][m]=newdata[it].converted_forest_land[k][m];
	    newdata[it].converted_forest_land[k][m]=ZEROVALUE;
	}
	  }else{
	/* checks for counting converted land in wood harvest */

	if(cdata[i].converted_forest_land_option == 2){
	  if(newdata[it].converted_forest_land[k][m] <= newdata[it].whd[k][m]){
            newdata[it].converted_forest_land_total[k][m]=newdata[it].converted_forest_land[k][m];
	    whr-=newdata[it].converted_forest_land[k][m];
	  }
	  else{
	    whr=ZEROVALUE;
	    newdata[it].converted_forest_land_total[k][m]=newdata[it].converted_forest_land[k][m];
    	    newdata[it].converted_forest_land[k][m]=newdata[it].whd[k][m];
	  }
	}
	else{
	    newdata[it].converted_forest_land_total[k][m]=newdata[it].converted_forest_land[k][m];
	    newdata[it].converted_forest_land[k][m]=ZEROVALUE;
	}
	  }

	if(whr>ZEROVALUE){

	  if (BEST_CASE){
	  /* if best case, Eurasia does secondary priority for wood, non-Eurasia does primay */
	  if(isnan(newdata[it].flowvs[k][m])){
	   
	  }
	  if((cdata[i].continent_code == 3) || (cdata[i].continent_code == 4)){	
	    if(secondary_biomass>=whr){
	      if(dstatic[k][m].fnf==1){
		newdata[it].sbh[k][m]= whr*1000;
	      }
	      else{
		newdata[it].sbh3[k][m]=whr*1000;		     
	      }
	      whr=ZEROVALUE;
	    }
	    else{
	      if(dstatic[k][m].fnf==1){
		newdata[it].sbh[k][m]=secondary_biomass*1000;
	      }
	      else{
		newdata[it].sbh3[k][m]=secondary_biomass*1000;		     
	      }
	      whr=whr-secondary_biomass;
	    }

	    if(virgin_biomass>whr){
	      if(dstatic[k][m].fnf==1){
		newdata[it].vbh[k][m]= whr*1000;
	      }
	      else{
		newdata[it].vbh2[k][m]=whr*1000;		     
	      }
	      newdata[it].flowvs[k][m]+=((newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*whr/virgin_biomass);

	 

	      whr=ZEROVALUE;
	  
	    }
	    else{
	      if(dstatic[k][m].fnf==1){
		newdata[it].vbh[k][m]=virgin_biomass*1000;
	      }
	      else{
		newdata[it].vbh2[k][m]=virgin_biomass*1000;		     
	      }
	      newdata[it].flowvs[k][m]+=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m]);

	      

             whr=whr-virgin_biomass;
	      
	    }

	  }
	  else{	 /* non-Eurasia */
	    if(virgin_biomass>=whr){
	      if(dstatic[k][m].fnf==1){
		newdata[it].vbh[k][m]= whr*1000;
	      }
	      else{
		newdata[it].vbh2[k][m]=whr*1000;		     
	      }
	      newdata[it].flowvs[k][m]+=((newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*whr/virgin_biomass);
	      
	 

whr=ZEROVALUE;
	      
	    }
	    else{
	      if(dstatic[k][m].fnf==1){
		newdata[it].vbh[k][m]=virgin_biomass*1000;
	      }
	      else{
		newdata[it].vbh2[k][m]=virgin_biomass*1000;		     
	      }
	      newdata[it].flowvs[k][m]+=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m]);
	  
	     
	 

  whr=whr-virgin_biomass;

	    }
	    if(secondary_biomass>=whr){
	      if(dstatic[k][m].fnf==1){
		newdata[it].sbh[k][m]= whr*1000;
	      }
	      else{
		newdata[it].sbh3[k][m]=whr*1000;		     
	      }
	      whr=ZEROVALUE;
	    }
	    else{
	      if(dstatic[k][m].fnf==1){
		newdata[it].sbh[k][m]=secondary_biomass*1000;
	      }
	      else{
		newdata[it].sbh3[k][m]=secondary_biomass*1000;		     
	      }
	      whr=whr-secondary_biomass;
	    }
	 
	  }


	  }else{

	  if(cdata[i].harvest_option == 2){

	    if(secondary_biomass>=whr){
	      if(dstatic[k][m].fnf==1){
		newdata[it].sbh[k][m]= whr*1000;
	      }
	      else{
		newdata[it].sbh3[k][m]=whr*1000;		     
	      }
	      whr=ZEROVALUE;
	    }
	    else{
	      if(dstatic[k][m].fnf==1){
		newdata[it].sbh[k][m]=secondary_biomass*1000;
	      }
	      else{
		newdata[it].sbh3[k][m]=secondary_biomass*1000;		     
	      }
	      whr=whr-secondary_biomass;
	    }

	    if((virgin_biomass>=whr) && (whr > ZEROVALUE)){
	      if(dstatic[k][m].fnf==1){
		newdata[it].vbh[k][m]= whr*1000;
	      }
	      else{
		newdata[it].vbh2[k][m]=whr*1000;		     
	      }
	      newdata[it].flowvs[k][m]+=((newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*whr/virgin_biomass);
	      whr=ZEROVALUE;
	    }
	    else if((virgin_biomass<whr) && (whr > ZEROVALUE)){
	      if(dstatic[k][m].fnf==1){
		newdata[it].vbh[k][m]=virgin_biomass*1000;
	      }
	      else{
		newdata[it].vbh2[k][m]=virgin_biomass*1000;		     
	      }
	      newdata[it].flowvs[k][m]+=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m]);
	      whr=whr-virgin_biomass;
	    }
	  
	  }
	  else{
	    if(virgin_biomass>=whr){
	      if(dstatic[k][m].fnf==1){
		newdata[it].vbh[k][m]= whr*1000;
	      }
	      else{
		newdata[it].vbh2[k][m]=whr*1000;		     
	      }
	      newdata[it].flowvs[k][m]+=((newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m])*whr/virgin_biomass);
	      whr=ZEROVALUE;
	    }
	    else{
	      if(dstatic[k][m].fnf==1){
		newdata[it].vbh[k][m]=virgin_biomass*1000;
	      }
	      else{
		newdata[it].vbh2[k][m]=virgin_biomass*1000;		     
	      }
	      newdata[it].flowvs[k][m]+=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvu[k][m]);
	      whr=whr-virgin_biomass;
	    }
	    if((secondary_biomass>=whr)&&(whr>ZEROVALUE)){
	      if(dstatic[k][m].fnf==1){
		newdata[it].sbh[k][m]= whr*1000;
	      }
	      else{
		newdata[it].sbh3[k][m]=whr*1000;		     
	      }
	      whr=ZEROVALUE;
	    }
	    else if((secondary_biomass<whr)&&(whr>ZEROVALUE)){
	      if(dstatic[k][m].fnf==1){
		newdata[it].sbh[k][m]=secondary_biomass*1000;
	      }
	      else{
		newdata[it].sbh3[k][m]=secondary_biomass*1000;		     
	      }
	      whr=whr-secondary_biomass;
	    }

	  }	  

	  }  /* end if for BEST_CASE if */	
	 
	  newdata[it].wh_unmet[k][m] = whr;

	} /* end of whd */
           }  /* end of gcode */
      } /* end of m */
    } /* end of k */


  return;

}

/********************************************************************/
void country_unmet_wh(int it, int i){
  int k, m;
  double vf=ZEROVALUE, vnf=ZEROVALUE, sf=ZEROVALUE, snf=ZEROVALUE;


  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){


#if COUNTRY_TEST

      if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){

#else
	if(dstatic[k][m].gcode > 0){
#endif


	  if(dstatic[k][m].gcode == cdata[i].ccode){
	    ctdata[i].whr+=newdata[it].wh_unmet[k][m];
	    ctdata[i].wh1+=(newdata[it].sbh[k][m]+newdata[it].sbh2[k][m]+newdata[it].sbh3[k][m]+newdata[it].vbh[k][m]+newdata[it].vbh2[k][m])/1000;
	    ctdata[i].orig_whd+=newdata[it].whd[k][m];
	  }
	  
           }  /* end of gcode */
	
	
      } /* end of m */
    } /* end of k */


  return;

}

/********************************************************************/
void update_vb3(int it){
  
  int iz, k, m, i;
  
  
  /* determine virgin biomass at the country level (vb); vb is the virgin
     available for a particular country for particular iz; therefore given
     a country code and a value of iz, vb is known */ 
  
  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){
      
#if REGION_TEST
	
      if(dstatic[k][m].gcode==COUNTRY_TEST_GCODE){
	   
#else
	if(dstatic[k][m].gcode > 0){
#endif
	     
	  i=dstatic[k][m].newgcode;
	  iz=newdata[it].zdis[k][m] ;
  

	  if(dstatic[k][m].fnf == 1){
		
		
	    /* vb in units = MgC */
	    cztdata[i][iz].vb+=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvs[k][m]-newdata[it].flowvu[k][m])*garea[k][m]*dstatic[k][m].vba/1000.0;
		
	  }


	}  
      } /* end of m */
    } /* end of k */
    
    
    return;
  }

/********************************************************************/
  void secondary_harvest3(int it, int i){

  int k, m;
  double sbh, sbhtest=ZEROVALUE; 



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

	      if(ctdata[i].whr > ZEROVALUE){

		sbh=(newdata[it].smb[k][m]-(newdata[it].sbh[k][m]+newdata[it].sbh2[k][m])/garea[k][m])*
		  prob_harv((newdata[it].smb[k][m]-(newdata[it].sbh[k][m]+newdata[it].sbh2[k][m])/garea[k][m]))*
		  newdata[it].s[k][m]*garea[k][m];


		if(sbh <= ctdata[i].whr*1000.){
		  newdata[it].sbh[k][m]+=sbh;
if(dstatic[k][m].rcode == 3) sbhtest+=sbh/1000.;
 
	      ctdata[i].stbh+=sbh/1000.;
	      
	      
	      ctdata[i].whr-=sbh/1000.;
	      
		}
		else{
		  newdata[it].sbh[k][m]+=ctdata[i].whr*1000.;
if(dstatic[k][m].rcode == 3) sbhtest+=ctdata[i].whr;
 
	      ctdata[i].stbh+=ctdata[i].whr;
	      
	      
	      ctdata[i].whr-=ctdata[i].whr;
		}
		
		sbh=0;
		
	      } /*end of whr */
	   
	      
	    } /* end of fnf */
	  }  /* end of gcode */
	}
	
	
      } /* end of m */
    } /* end of k */


    return;

  }

/***********************************************************************/
  void virgin_harvest3(int it, int i,int curryear){
                      
  FILE *testfile;
  int k, m, iz, izz, zmax=0, im, j;
  double total_avail, vbhsummer=ZEROVALUE, flowvs;

 
  testfile=fopen("zmax.test",outstat);

  
  /* first determine zmax; the maximum # cells away from the focal cell (the agricultural
     cell) needed to go to attain wood harvest demand.
     zmax=0; have enough biomass in focal (agricultural) cell 
     zmax=1; need to go to the next adjacent cell 
     zmax=2; etc...zmax=10 or MAXZ (MAXZ=11; maximum z before we get tired, and then 
     we spread remaining harvest over all forested cells with iz >= this value */
  
 
  if(ctdata[i].vwh < cztdata[i][0].vb)
    zmax=0;
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb){
    zmax=1;

  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb){
    zmax=2;

  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb){
    zmax=3;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb){
    zmax=4;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb){
    zmax=5;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb){
    zmax=6;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb){
    zmax=7;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +                +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb){
    zmax=8;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb){
    zmax=9;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb){
    zmax=10;
  }

#if MAXZ_21	        
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb){
    zmax=11;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb){
    zmax=12;
  }	 
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb+cztdata[i][13].vb){
    zmax=13;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb+cztdata[i][13].vb+cztdata[i][14].vb){
    zmax=14;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb+cztdata[i][13].vb+cztdata[i][14].vb+cztdata[i][15].vb){
    zmax=15;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb+cztdata[i][13].vb+cztdata[i][14].vb+cztdata[i][15].vb+cztdata[i][16].vb){
    zmax=16;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb+cztdata[i][13].vb+cztdata[i][14].vb+cztdata[i][15].vb+cztdata[i][16].vb+cztdata[i][17].vb){
    zmax=17;
  }	
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb+cztdata[i][13].vb+cztdata[i][14].vb+cztdata[i][15].vb+cztdata[i][16].vb+cztdata[i][17].vb+cztdata[i][18].vb){
    zmax=18;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb+cztdata[i][13].vb+cztdata[i][14].vb+cztdata[i][15].vb+cztdata[i][16].vb+cztdata[i][17].vb+cztdata[i][18].vb+cztdata[i][19].vb){
    zmax=19;
  }
  else if(ctdata[i].vwh < cztdata[i][0].vb+cztdata[i][1].vb
	  +cztdata[i][2].vb
	  +cztdata[i][3].vb+cztdata[i][4].vb+cztdata[i][5].vb+
	  cztdata[i][6].vb+cztdata[i][7].vb+cztdata[i][8].vb
	  +cztdata[i][9].vb+cztdata[i][10].vb+cztdata[i][11].vb+cztdata[i][12].vb+cztdata[i][13].vb+cztdata[i][14].vb+cztdata[i][15].vb+cztdata[i][16].vb+cztdata[i][17].vb+cztdata[i][18].vb+cztdata[i][19].vb+cztdata[i][20].vb){
    zmax=20;
  }	 
#endif

  else {
    zmax=MAXZ;
  }

  if(zmax == 0){
    ctdata[i].wh_at_zmax=ctdata[i].whr;
  }
  else{
    for(j=0;j<zmax;j++) ctdata[i].wh_lessthan_zmax+=cztdata[i][j].vb;
    ctdata[i].wh_at_zmax=ctdata[i].whr-ctdata[i].wh_lessthan_zmax;
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
	 
	    if(ctdata[i].whr > ZEROVALUE) {


	      iz=newdata[it].zdis[k][m] ;





	      if(dstatic[k][m].fnf == 1){



		/* determine how to take biomass for harvest demand */

		if(zmax<MAXZ){ /*business as usual*/


		  if(iz < zmax){ /* take as much as possible to fulfill harvest demand */
		      
		    flowvs=newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvs[k][m]-newdata[it].flowvu[k][m]; 
		    
		    newdata[it].vbh[k][m]+=flowvs*garea[k][m]*dstatic[k][m].vba;                            
		    newdata[it].flowvs[k][m]+=flowvs;
		       
		    
		  }
		  else if(iz == zmax){


		    if(cztdata[i][iz].vb > ZEROVALUE){
		      
		     
		      flowvs=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvs[k][m]-newdata[it].flowvu[k][m])*ctdata[i].wh_at_zmax/cztdata[i][iz].vb;
		      
		      newdata[it].vbh[k][m]+=flowvs*garea[k][m]*dstatic[k][m].vba;
                      newdata[it].flowvs[k][m]+=flowvs;
		    }	  
	

		  }

		  if(it ==0 && dstatic[k][m].rcode == 3) {

		    vbhsummer+=newdata[it].vbh[k][m]/1000;


		  }



		} /* end of business as usual, zmax < MAXZ */

		else{ 

		  if(iz < MAXZ){ /*zmax=MAXZ*/
		    /*most dry countries with reported wood harvest and not enough forest*/
		    /* take everything possible */
		    
		   
		    flowvs=(newdata[it].v[k][m]-
		       newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvs[k][m]-newdata[it].flowvu[k][m]);
		    
		    newdata[it].vbh[k][m]+=flowvs*garea[k][m]*
		      dstatic[k][m].vba;
		     newdata[it].flowvs[k][m]+=flowvs;
		    
		  }
		  else if(iz<=NZ){ /*from MAXZ=zmax to NZ*/
		    
		    /* rare, reported wood harvest, but have to go more than MAXZ to get it*/
		    /* take what is needed proportionally out of each cell */

		    total_avail=ZEROVALUE;
		    
		    for(izz=MAXZ;izz<NZ;izz++)
		      total_avail+=cztdata[i][izz].vb;
		    
		    if(total_avail > ZEROVALUE){
		      
		      if(ctdata[i].wh_at_zmax <= total_avail){
			
			
                         flowvs=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvs[k][m]-newdata[it].flowvu[k][m])*ctdata[i].wh_at_zmax/total_avail;
		      }
		      else{
			
			flowvs=(newdata[it].v[k][m]-
			  newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvs[k][m]-newdata[it].flowvu[k][m]);
		      }
		    }
		    
		    /*safe, but still may not get enough harvest*/
		   
		    newdata[it].vbh[k][m]+=flowvs*garea[k][m]*
		      dstatic[k][m].vba;
		    newdata[it].flowvs[k][m]+=flowvs;
		    
		    fprintf(testfile,"potential zmax trouble spot, %d %s zmax %d iz %d flowvs %lf wh_at_zmax %lf total %lf\n",curryear,cname[i],zmax,iz,newdata[it].flowvs[k][m],ctdata[i].wh_at_zmax,total_avail);  
		    

		  }
		  else {
		    printf("iz %d zmax %d\n",iz,zmax);
		  }
		  
		}  /* end of else iz < MAXZ */



	      }/* end of fnf*/

	      ctdata[i].whr-=(flowvs*garea[k][m]*dstatic[k][m].vba)/1000.;
	      flowvs=0;
	    }/* end of whr*/


	  }/*gcode*/
	}


      }  /* end of m */
    } /* end of k */

    fclose(testfile);

 
    return;

  }
/***********************************************************************/
void force_harvest3(int it, int i, int curryear){

  FILE *testfile;
  int k, m, flag=0;
  double whr_orig=ZEROVALUE, flowvs_add=ZEROVALUE, flowss_add=ZEROVALUE, wh_from_sbh_only=ZEROVALUE, sbh2, sbh3, vbh2;
  double tester1=ZEROVALUE, tester2=ZEROVALUE, tester3=ZEROVALUE, tester4=ZEROVALUE, tester5=ZEROVALUE;

  
   
 
 
  /* hold original amount of whr for spreading proportionally */

  whr_orig=ctdata[i].whr;



  /* compute country level smb, smb_nf, vnfb, s area for force harvest */
  
  
  for (k=0;k<NY;k++){
    for (m=0;m<NX;m++){
      
      
#if REGION_TEST
      
      if(dstatic[k][m].gcode==REGION_TEST_GCODE){
	   
#else
	if(dstatic[k][m].gcode > 0){
#endif
	     
	  if(dstatic[k][m].gcode == cdata[i].ccode){	   
	 	     
	    if(dstatic[k][m].fnf == 0){
	     
	      /* country level virgin nonforested biomass */
	      ctdata[i].vnfb+=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvs[k][m]-newdata[it].flowvu[k][m])*garea[k][m]*dstatic[k][m].vba/1000.;
		 
	      /* country level secondary nonforested biomass */
	      ctdata[i].smb_nf+=(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowvu[k][m])*garea[k][m]*(newdata[it].smb[k][m]-newdata[it].sbh3[k][m]/garea[k][m])/1000.;

	      /* country level secondary nonforested area */
	      ctdata[i].sarea_nf+=(newdata[it].s[k][m]+newdata[it].flowps[k][m]+newdata[it].flowcs[k][m]+newdata[it].flowvs[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m]+newdata[it].flowus[k][m])*garea[k][m];
 
	    }
	    else {

	      /* country level secondary forest biomass from sbh only */
	      wh_from_sbh_only+=newdata[it].sbh[k][m]/1000.;
	     
	      /* country level secondary forest biomass */
	      ctdata[i].smb+=(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*((garea[k][m]*newdata[it].smb[k][m]-newdata[it].sbh[k][m]-newdata[it].sbh2[k][m])/1000.); 

	     /* country level secondary forested area */
	      ctdata[i].sarea+=(newdata[it].s[k][m]+newdata[it].flowps[k][m]+newdata[it].flowcs[k][m]+newdata[it].flowvs[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m]+newdata[it].flowus[k][m])*garea[k][m];
	      
	    }
	
   
	  }  /* end of gcode */
	}
       
      } /* end of m */
    } /* end of k */
   

    /* 4 cases for getting remaining wood harvest (whr) */


    if(ctdata[i].smb >= ctdata[i].whr){


      flag=1;  /* case 1, taking needed younger secondary forest proportionally to satisfy whr */                                 

      for (k=0;k<NY;k++){
	for (m=0;m<NX;m++){
	 
#if REGION_TEST
	 
	  if(dstatic[k][m].gcode==REGION_TEST_GCODE){
	   
#else
	    if(dstatic[k][m].gcode > 0){
#endif
	     
	      if(dstatic[k][m].gcode == cdata[i].ccode){

	       
		if(dstatic[k][m].fnf == 1){
		 
		  if(ctdata[i].whr > ZEROVALUE){
		   
		   
		   sbh2=whr_orig/ctdata[i].smb*(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*(garea[k][m]*newdata[it].smb[k][m]-newdata[it].sbh[k][m]-newdata[it].sbh2[k][m]);
		   
		    ctdata[i].whr-=sbh2/1000.;
		    newdata[it].sbh2[k][m]+=sbh2;
		   
		  } /*end of whr */
		} /*end of fnf */
	      }  /* end of gcode */
	    }
	   
	  } /* end of m */
	} /* end of k */

       
      }/* end case 1*/
     

     
      else if((cdata[i].smart_flow_option == 1) && ((ctdata[i].smb+ctdata[i].vnfb) >= ctdata[i].whr)){    

	flag=2;  /* case 2, primary priority, not enough young secondary forest to satisfy demand, 
		   therefore taking all young secondary forest and needed virgin nonforest proportionally 
		   to satisfy whr */

  

	for (k=0;k<NY;k++){
	  for (m=0;m<NX;m++){
	   
#if REGION_TEST
	 
	    if(dstatic[k][m].gcode==REGION_TEST_GCODE){
	   
#else
	      if(dstatic[k][m].gcode > 0){
#endif
	     
		if(dstatic[k][m].gcode == cdata[i].ccode){

		 
		  if(ctdata[i].whr > ZEROVALUE){
		   
		    if(dstatic[k][m].fnf == 1) {

		      sbh2=(newdata[it].smb[k][m]*garea[k][m]-newdata[it].sbh[k][m]-newdata[it].sbh2[k][m])*(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m]);

		      ctdata[i].whr-=sbh2/1000.;
		      newdata[it].sbh2[k][m]+=sbh2;
ctdata[i].fh_sbh2+=sbh2/1000.;
		    } /* end of fnf*/		   
		  } /*end of whr*/
		}  /* end of gcode */
	      }
	     
	    } /* end of m */
	  } /* end of k */



	 /* take the whr proportionally out of vnfb; note needed to consider how much was taken
            out of secondary forest (ctdata[i].smb) */


	 for (k=0;k<NY;k++){
	   for (m=0;m<NX;m++){
	     
#if REGION_TEST
	     
	     if(dstatic[k][m].gcode==REGION_TEST_GCODE){
	   
#else
	       if(dstatic[k][m].gcode > 0){
#endif
	     
		 if(dstatic[k][m].gcode == cdata[i].ccode){

		 
		   if(ctdata[i].whr > ZEROVALUE){
		   
		     if(dstatic[k][m].fnf == 0) {
		     
	     	       flowvs_add=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvs[k][m]-newdata[it].flowvu[k][m])*(whr_orig-ctdata[i].smb)/ctdata[i].vnfb; 
		      
		       if(flowvs_add <= -ZEROVALUE) {flowvs_add = ZEROVALUE;}
		     
		       newdata[it].flowvs[k][m]+=flowvs_add;
		     
		       newdata[it].vbh2[k][m]+=flowvs_add*dstatic[k][m].vba*garea[k][m];
		     
		       ctdata[i].whr-=flowvs_add*dstatic[k][m].vba*garea[k][m]/1000.;
		  ctdata[i].fh_vbh2+=flowvs_add*dstatic[k][m].vba*garea[k][m]/1000.;   
		     } /*end of fnf */
		   
  
		   } /*end of whr*/
		 
		 }  /* end of gcode */
	       }
	     
	     } /* end of m */
	   } /* end of k */


	 } /*end case 2 */
	 

	 
	 else if((cdata[i].smart_flow_option == 2) && ((ctdata[i].smb+ctdata[i].smb_nf) >= ctdata[i].whr)){  

	   flag=3;  /* case 3, secondary priority, taking all young secondary forest and needed secondary 
                   nonforest proportionally to satisfy whr */


     
	   for (k=0;k<NY;k++){
	     for (m=0;m<NX;m++){
	   
#if REGION_TEST
	 
	       if(dstatic[k][m].gcode==REGION_TEST_GCODE){
	   
#else
		 if(dstatic[k][m].gcode > 0){
#endif
	     
		   if(dstatic[k][m].gcode == cdata[i].ccode){

		 
		     if(ctdata[i].whr > ZEROVALUE){
		   
		       if(dstatic[k][m].fnf == 1) {

			 sbh2=(newdata[it].smb[k][m]*garea[k][m]-newdata[it].sbh[k][m]-newdata[it].sbh2[k][m])*(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m]);

			 ctdata[i].whr-=sbh2/1000.;
			 newdata[it].sbh2[k][m]+=sbh2;
		       } /* end of fnf*/		   
		     } /*end of whr*/
		   }  /* end of gcode */
		 }
	     
	       } /* end of m */
	     } /* end of k */


	     /* take the whr proportionally out of secondary nonforest; note needed to consider how much was taken
		out of secondary forest (ctdata[i].smb) */


	     for (k=0;k<NY;k++){
	       for (m=0;m<NX;m++){
	     
#if REGION_TEST
	     
		 if(dstatic[k][m].gcode==REGION_TEST_GCODE){
	   
#else
		   if(dstatic[k][m].gcode > 0){
#endif
	     
		     if(dstatic[k][m].gcode == cdata[i].ccode){
		       
		       if(ctdata[i].whr > ZEROVALUE){
		   
			 if(dstatic[k][m].fnf == 0) {
		     
			    flowss_add=(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*(whr_orig-ctdata[i].smb)/ctdata[i].smb_nf;
			   
		 
			   if(flowss_add <= -ZEROVALUE) {flowss_add = ZEROVALUE;}

			   sbh3=flowss_add*(newdata[it].smb[k][m]-newdata[it].sbh3[k][m]/garea[k][m])*garea[k][m];
		     
			   ctdata[i].whr-=sbh3/1000.;
			   newdata[it].sbh3[k][m]+=sbh3;

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

	    
#if REGION_TEST
	 
		   if(dstatic[k][m].gcode==REGION_TEST_GCODE){
	   
#else
		     if(dstatic[k][m].gcode > 0){
#endif
	     
		       if(dstatic[k][m].gcode == cdata[i].ccode){

		   
			 if(dstatic[k][m].fnf == 1){

			   /* take all young secondary forest */
		       
			   sbh2=(newdata[it].smb[k][m]*garea[k][m]-newdata[it].sbh[k][m]-newdata[it].sbh2[k][m])*(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m]);

			   ctdata[i].whr-=sbh2/1000.;
ctdata[i].fh_sbh2+=sbh2/1000.;
 newdata[it].sbh2[k][m]+=sbh2;
			 }
		       }
		     }
		   }
		 }

		

	       for (k=0;k<NY;k++){
		 for (m=0;m<NX;m++){

	    
#if REGION_TEST
	 
		   if(dstatic[k][m].gcode==REGION_TEST_GCODE){
	   
#else
		     if(dstatic[k][m].gcode > 0){
#endif
	     
		       if(dstatic[k][m].gcode == cdata[i].ccode){

		   
			 if(dstatic[k][m].fnf == 0){

			   if(cdata[i].smart_flow_option == 1){

			     /* primary priority, first take all vnfb */


			     flowvs_add=newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvs[k][m]-newdata[it].flowvu[k][m];
		       
			     if(flowvs_add <= -ZEROVALUE) {flowvs_add = ZEROVALUE;}
		       
			     newdata[it].flowvs[k][m]+=flowvs_add;
		     
			     newdata[it].vbh2[k][m]+=flowvs_add*dstatic[k][m].vba*garea[k][m];
		     
			     ctdata[i].whr-=flowvs_add*dstatic[k][m].vba*garea[k][m]/1000.;
ctdata[i].fh_vbh2+=flowvs_add*dstatic[k][m].vba*garea[k][m]/1000.;
  


			     /* primary priority, not enough smb and vnfb, determined how much can be taken from smb_nf */

			     /* first try taking it proportionally */

 if(ctdata[i].smb_nf >= (whr_orig-ctdata[i].smb-ctdata[i].vnfb)){

    
			         flowss_add=(newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m])*(whr_orig-ctdata[i].smb-ctdata[i].vnfb)/ctdata[i].smb_nf; 
			      
			   
			       if(flowss_add <= -ZEROVALUE) {flowss_add = ZEROVALUE;}
			   
			       sbh3=flowss_add*(newdata[it].smb[k][m]-newdata[it].sbh3[k][m]/garea[k][m])*garea[k][m];
			   
			           ctdata[i].whr-sbh3/1000.;   
			        ctdata[i].fh_sbh3+=sbh3/1000.; 
				newdata[it].sbh3[k][m]+=sbh3;
			       
			 			 
			     }
			     else{ /* not enough in smb_nf, therefore take all remaining smb_nf, the rest is unmet */
			 
			       flowss_add=newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m];
			   
			       if(flowss_add <= -ZEROVALUE) {flowss_add = ZEROVALUE;}
			   
			       sbh3=flowss_add*(newdata[it].smb[k][m]-newdata[it].sbh3[k][m]/garea[k][m])*garea[k][m];

			        ctdata[i].whr-=sbh3/1000.; 
			        ctdata[i].fh_sbh3+=sbh3/1000.; 
				newdata[it].sbh3[k][m]+=sbh3;
}


			   } 
			   else { /* if smart_flow_option == 2 */
		       
			     /* secondary priority, first take all smb_nf */

			     flowss_add=newdata[it].s[k][m]-newdata[it].flowsc[k][m]-newdata[it].flowsp[k][m]-newdata[it].flowsu[k][m];
			   
			     if(flowss_add <= -ZEROVALUE) {flowss_add = ZEROVALUE;}
			   
			     sbh3=flowss_add*(newdata[it].smb[k][m]-newdata[it].sbh3[k][m]/garea[k][m])*garea[k][m];

			     ctdata[i].whr-=sbh3/1000.;
ctdata[i].fh_sbh3+=sbh3/1000.;
 newdata[it].sbh3[k][m]+=sbh3;

			     /* secondary priority, not enough smb and smb_nf, determined how much can be 
                          taken from vnfb */

			     /* first try taking it proportionally */

 if(ctdata[i].vnfb >= (whr_orig-ctdata[i].smb-ctdata[i].smb_nf)){
     
		    
			         flowvs_add=(newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvs[k][m]-newdata[it].flowvu[k][m])*(whr_orig-ctdata[i].smb-ctdata[i].smb_nf)/ctdata[i].vnfb;
			      

			       if(flowvs_add <= -ZEROVALUE) {flowvs_add = ZEROVALUE;}

			       newdata[it].flowvs[k][m]+=flowvs_add;
		     
			       newdata[it].vbh2[k][m]+=flowvs_add*dstatic[k][m].vba*garea[k][m];
		     
			       ctdata[i].whr-=flowvs_add*dstatic[k][m].vba*garea[k][m]/1000.;
	   ctdata[i].fh_vbh2+=flowvs_add*dstatic[k][m].vba*garea[k][m]/1000.;
			     }
			     else{ /* not enough in vnfb, therefore take all remaining vnfb, the rest is unmet */


			       flowvs_add=newdata[it].v[k][m]-newdata[it].flowvc[k][m]-newdata[it].flowvp[k][m]-newdata[it].flowvs[k][m]-newdata[it].flowvu[k][m];
		       
			       if(flowvs_add <= -ZEROVALUE) {flowvs_add = ZEROVALUE;}
		       
			       newdata[it].flowvs[k][m]+=flowvs_add;
		     
			       newdata[it].vbh2[k][m]+=flowvs_add*dstatic[k][m].vba*garea[k][m];
		     
			       ctdata[i].whr-=flowvs_add*dstatic[k][m].vba*garea[k][m]/1000.;
ctdata[i].fh_vbh2+=flowvs_add*dstatic[k][m].vba*garea[k][m]/1000.;
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

#if REGION_TEST
	 
		   if(dstatic[k][m].gcode==REGION_TEST_GCODE){
		     
#else
		     if(dstatic[k][m].gcode > 0){
#endif
	     
		       if(dstatic[k][m].gcode == cdata[i].ccode){

		 
			   tester1+=newdata[it].vbh[k][m]/1000.;
			   tester2+=newdata[it].vbh2[k][m]/1000.;
			   tester3+=newdata[it].sbh[k][m]/1000.;
			   tester4+=newdata[it].sbh2[k][m]/1000.;
			   tester5+=newdata[it].sbh3[k][m]/1000.;

		   
		       }
		     }
		   }
		 }

		 testfile=fopen("wh.latest",outstat);
	   
		 fprintf(testfile,"%d %s flag %d fao# %lf wh_needed %lf vbh %lf vbh2 %lf sbh %lf sbh2 %lf sbh3 %lf whtot %lf smb %lf smb_nf %lf vnfb %lf\n",curryear,cname[i],flag,ctdata[i].wh,whr_orig,tester1,tester2,tester3,tester4,tester5,(tester1+tester2+tester3+tester4+tester5),ctdata[i].smb,ctdata[i].smb_nf,ctdata[i].vnfb);
		   
		 
		 fclose(testfile);
	   

		 return;
	 
	       }

/***********************************************************************/
/* round number n to d decimal points */

inline double fround(double n, unsigned d)
{
  return floor(n * pow(10., d) + .5) / pow(10., d);
}
/********************************************************************/
     int
       mkdir2 (char *str1, char *str2, int mode)
     {
       char *name = (char *) alloca (strlen (str1) + strlen (str2) + 1);
       strcat (strcpy (name, str1), str2);
       return mkdir(name, mode);
     }
/********************************************************************/
char * strcatjt(char *dest, const char *src) {
  if (strlen(src) + 1 >
      sizeof(dest) - strlen(dest))
    err(1, "dest would be truncated");
  (void)strncat(dest, src,sizeof(dest) - strlen(dest) - 1);
  return dest;
}
/********************************************************************/
int open2 (char *str1, char *str2, int flags, int mode)
{
  char *name = (char *) alloca (strlen (str1) + strlen (str2) + 1);
  strcpy (name, str1);
  strcat (name, str2);
  return open (name, flags, mode);
}
/********************************************************************/
void create_restart_file(void)
{
	FILE	*	ini ;
	int status;
        char start_yearc[4],stop_yearc[4];
        sprintf(start_yearc, "%d", start_year);
        sprintf(stop_yearc, "%d", stop_year);
        strcpy(casestr,casename);
        status = unlink("restart_glm.ini");
       	ini = fopen("restart_glm.ini", "w");
	strcpy(restart_file,strcat(casestr,".glm.state."));
	strcat(restart_file,strcat(stop_yearc,".nc"));
	fprintf(ini,
		"#\n"
		"# This is a restart file for glm \n"
		"#\n"
		"\n"
		"[restart]\n"
		"\n"
		"restart_file        = %s \n"
		"restart_year        = %d \n"
		"\n"
		"\n",restart_file,stop_year);
   
	fclose(ini);
}
/********************************************************************/
void read_restart(void)
{
	dictionary*	ini_rest;
	int status;

        ini_rest         = iniparser_load("restart_glm.ini");
        restart_filename = iniparser_getstring(ini_rest, "restart:restart_file", NULL);
        start_year       = iniparser_getint(ini_rest, "restart:restart_year", 1500);

        if (start_year >= stop_year) {
	  fprintf(stderr, "glm:read_restart: start_year of restart run < stop_year\n");
	  exit(1);
	}
}
/*-------------------------------------------------------------------------*/
/**
   @file    iniparser.c
   @author  N. Devillard
   @date    Sep 2007
   @version 3.0
   @brief   Parser for ini files.
*/
/*--------------------------------------------------------------------------*/
/*
    $Id: iniparser.c,v 2.18 2008-01-03 18:35:39 ndevilla Exp $
    $Revision: 2.18 $
    $Date: 2008-01-03 18:35:39 $
*/
/*---------------------------- Includes ------------------------------------*/
#include <ctype.h>
#include "iniparser.h"

/*---------------------------- Defines -------------------------------------*/
#define ASCIILINESZ         (1024)
#define INI_INVALID_KEY     ((char*)-1)

/*---------------------------------------------------------------------------
                        Private to this module
 ---------------------------------------------------------------------------*/
/**
 * This enum stores the status for each parsed line (internal use only).
 */
typedef enum _line_status_ {
    LINE_UNPROCESSED,
    LINE_ERROR,
    LINE_EMPTY,
    LINE_COMMENT,
    LINE_SECTION,
    LINE_VALUE
} line_status ;

/*-------------------------------------------------------------------------*/
/**
  @brief	Convert a string to lowercase.
  @param	s	String to convert.
  @return	ptr to statically allocated string.

  This function returns a pointer to a statically allocated string
  containing a lowercased version of the input string. Do not free
  or modify the returned string! Since the returned string is statically
  allocated, it will be modified at each function call (not re-entrant).
 */
/*--------------------------------------------------------------------------*/
static char * strlwc(const char * s)
{
    static char l[ASCIILINESZ+1];
    int i ;

    if (s==NULL) return NULL ;
    memset(l, 0, ASCIILINESZ+1);
    i=0 ;
    while (s[i] && i<ASCIILINESZ) {
        l[i] = (char)tolower((int)s[i]);
        i++ ;
    }
    l[ASCIILINESZ]=(char)0;
    return l ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Remove blanks at the beginning and the end of a string.
  @param	s	String to parse.
  @return	ptr to statically allocated string.

  This function returns a pointer to a statically allocated string,
  which is identical to the input string, except that all blank
  characters at the end and the beg. of the string have been removed.
  Do not free or modify the returned string! Since the returned string
  is statically allocated, it will be modified at each function call
  (not re-entrant).
 */
/*--------------------------------------------------------------------------*/
static char * strstrip(char * s)
{
    static char l[ASCIILINESZ+1];
	char * last ;
	
    if (s==NULL) return NULL ;
    
	while (isspace((int)*s) && *s) s++;
	memset(l, 0, ASCIILINESZ+1);
	strcpy(l, s);
	last = l + strlen(l);
	while (last > l) {
		if (!isspace((int)*(last-1)))
			break ;
		last -- ;
	}
	*last = (char)0;
	return (char*)l ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Get number of sections in a dictionary
  @param    d   Dictionary to examine
  @return   int Number of sections found in dictionary

  This function returns the number of sections found in a dictionary.
  The test to recognize sections is done on the string stored in the
  dictionary: a section name is given as "section" whereas a key is
  stored as "section:key", thus the test looks for entries that do not
  contain a colon.

  This clearly fails in the case a section name contains a colon, but
  this should simply be avoided.

  This function returns -1 in case of error.
 */
/*--------------------------------------------------------------------------*/
int iniparser_getnsec(dictionary * d)
{
    int i ;
    int nsec ;

    if (d==NULL) return -1 ;
    nsec=0 ;
    for (i=0 ; i<d->size ; i++) {
        if (d->key[i]==NULL)
            continue ;
        if (strchr(d->key[i], ':')==NULL) {
            nsec ++ ;
        }
    }
    return nsec ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Get name for section n in a dictionary.
  @param    d   Dictionary to examine
  @param    n   Section number (from 0 to nsec-1).
  @return   Pointer to char string

  This function locates the n-th section in a dictionary and returns
  its name as a pointer to a string statically allocated inside the
  dictionary. Do not free or modify the returned string!

  This function returns NULL in case of error.
 */
/*--------------------------------------------------------------------------*/
char * iniparser_getsecname(dictionary * d, int n)
{
    int i ;
    int foundsec ;

    if (d==NULL || n<0) return NULL ;
    foundsec=0 ;
    for (i=0 ; i<d->size ; i++) {
        if (d->key[i]==NULL)
            continue ;
        if (strchr(d->key[i], ':')==NULL) {
            foundsec++ ;
            if (foundsec>n)
                break ;
        }
    }
    if (foundsec<=n) {
        return NULL ;
    }
    return d->key[i] ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Dump a dictionary to an opened file pointer.
  @param    d   Dictionary to dump.
  @param    f   Opened file pointer to dump to.
  @return   void

  This function prints out the contents of a dictionary, one element by
  line, onto the provided file pointer. It is OK to specify @c stderr
  or @c stdout as output files. This function is meant for debugging
  purposes mostly.
 */
/*--------------------------------------------------------------------------*/
void iniparser_dump(dictionary * d, FILE * f)
{
    int     i ;

    if (d==NULL || f==NULL) return ;
    for (i=0 ; i<d->size ; i++) {
        if (d->key[i]==NULL)
            continue ;
        if (d->val[i]!=NULL) {
            fprintf(f, "[%s]=[%s]\n", d->key[i], d->val[i]);
        } else {
            fprintf(f, "[%s]=UNDEF\n", d->key[i]);
        }
    }
    return ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Save a dictionary to a loadable ini file
  @param    d   Dictionary to dump
  @param    f   Opened file pointer to dump to
  @return   void

  This function dumps a given dictionary into a loadable ini file.
  It is Ok to specify @c stderr or @c stdout as output files.
 */
/*--------------------------------------------------------------------------*/
void iniparser_dump_ini(dictionary * d, FILE * f)
{
    int     i, j ;
    char    keym[ASCIILINESZ+1];
    int     nsec ;
    char *  secname ;
    int     seclen ;

    if (d==NULL || f==NULL) return ;

    nsec = iniparser_getnsec(d);
    if (nsec<1) {
        /* No section in file: dump all keys as they are */
        for (i=0 ; i<d->size ; i++) {
            if (d->key[i]==NULL)
                continue ;
            fprintf(f, "%s = %s\n", d->key[i], d->val[i]);
        }
        return ;
    }
    for (i=0 ; i<nsec ; i++) {
        secname = iniparser_getsecname(d, i) ;
        seclen  = (int)strlen(secname);
        fprintf(f, "\n[%s]\n", secname);
        sprintf(keym, "%s:", secname);
        for (j=0 ; j<d->size ; j++) {
            if (d->key[j]==NULL)
                continue ;
            if (!strncmp(d->key[j], keym, seclen+1)) {
                fprintf(f,
                        "%-30s = %s\n",
                        d->key[j]+seclen+1,
                        d->val[j] ? d->val[j] : "");
            }
        }
    }
    fprintf(f, "\n");
    return ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Get the string associated to a key
  @param    d       Dictionary to search
  @param    key     Key string to look for
  @param    def     Default value to return if key not found.
  @return   pointer to statically allocated character string

  This function queries a dictionary for a key. A key as read from an
  ini file is given as "section:key". If the key cannot be found,
  the pointer passed as 'def' is returned.
  The returned char pointer is pointing to a string allocated in
  the dictionary, do not free or modify it.
 */
/*--------------------------------------------------------------------------*/
char * iniparser_getstring(dictionary * d, const char * key, char * def)
{
    char * lc_key ;
    char * sval ;

    if (d==NULL || key==NULL)
        return def ;

    lc_key = strlwc(key);
    sval = dictionary_get(d, lc_key, def);
    return sval ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Get the string associated to a key, convert to an int
  @param    d Dictionary to search
  @param    key Key string to look for
  @param    notfound Value to return in case of error
  @return   integer

  This function queries a dictionary for a key. A key as read from an
  ini file is given as "section:key". If the key cannot be found,
  the notfound value is returned.

  Supported values for integers include the usual C notation
  so decimal, octal (starting with 0) and hexadecimal (starting with 0x)
  are supported. Examples:

  "42"      ->  42
  "042"     ->  34 (octal -> decimal)
  "0x42"    ->  66 (hexa  -> decimal)

  Warning: the conversion may overflow in various ways. Conversion is
  totally outsourced to strtol(), see the associated man page for overflow
  handling.

  Credits: Thanks to A. Becker for suggesting strtol()
 */
/*--------------------------------------------------------------------------*/
int iniparser_getint(dictionary * d, const char * key, int notfound)
{
    char    *   str ;

    str = iniparser_getstring(d, key, INI_INVALID_KEY);
    if (str==INI_INVALID_KEY) return notfound ;
    return (int)strtol(str, NULL, 0);
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Get the string associated to a key, convert to a double
  @param    d Dictionary to search
  @param    key Key string to look for
  @param    notfound Value to return in case of error
  @return   double

  This function queries a dictionary for a key. A key as read from an
  ini file is given as "section:key". If the key cannot be found,
  the notfound value is returned.
 */
/*--------------------------------------------------------------------------*/
double iniparser_getdouble(dictionary * d, char * key, double notfound)
{
    char    *   str ;

    str = iniparser_getstring(d, key, INI_INVALID_KEY);
    if (str==INI_INVALID_KEY) return notfound ;
    return atof(str);
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Get the string associated to a key, convert to a boolean
  @param    d Dictionary to search
  @param    key Key string to look for
  @param    notfound Value to return in case of error
  @return   integer

  This function queries a dictionary for a key. A key as read from an
  ini file is given as "section:key". If the key cannot be found,
  the notfound value is returned.

  A true boolean is found if one of the following is matched:

  - A string starting with 'y'
  - A string starting with 'Y'
  - A string starting with 't'
  - A string starting with 'T'
  - A string starting with '1'

  A false boolean is found if one of the following is matched:

  - A string starting with 'n'
  - A string starting with 'N'
  - A string starting with 'f'
  - A string starting with 'F'
  - A string starting with '0'

  The notfound value returned if no boolean is identified, does not
  necessarily have to be 0 or 1.
 */
/*--------------------------------------------------------------------------*/
int iniparser_getboolean(dictionary * d, const char * key, int notfound)
{
    char    *   c ;
    int         ret ;

    c = iniparser_getstring(d, key, INI_INVALID_KEY);
    if (c==INI_INVALID_KEY) return notfound ;
    if (c[0]=='y' || c[0]=='Y' || c[0]=='1' || c[0]=='t' || c[0]=='T') {
        ret = 1 ;
    } else if (c[0]=='n' || c[0]=='N' || c[0]=='0' || c[0]=='f' || c[0]=='F') {
        ret = 0 ;
    } else {
        ret = notfound ;
    }
    return ret;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Finds out if a given entry exists in a dictionary
  @param    ini     Dictionary to search
  @param    entry   Name of the entry to look for
  @return   integer 1 if entry exists, 0 otherwise

  Finds out if a given entry exists in the dictionary. Since sections
  are stored as keys with NULL associated values, this is the only way
  of querying for the presence of sections in a dictionary.
 */
/*--------------------------------------------------------------------------*/
int iniparser_find_entry(
    dictionary  *   ini,
    char        *   entry
)
{
    int found=0 ;
    if (iniparser_getstring(ini, entry, INI_INVALID_KEY)!=INI_INVALID_KEY) {
        found = 1 ;
    }
    return found ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Set an entry in a dictionary.
  @param    ini     Dictionary to modify.
  @param    entry   Entry to modify (entry name)
  @param    val     New value to associate to the entry.
  @return   int 0 if Ok, -1 otherwise.

  If the given entry can be found in the dictionary, it is modified to
  contain the provided value. If it cannot be found, -1 is returned.
  It is Ok to set val to NULL.
 */
/*--------------------------------------------------------------------------*/
int iniparser_set(dictionary * ini, char * entry, char * val)
{
    return dictionary_set(ini, strlwc(entry), val) ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Delete an entry in a dictionary
  @param    ini     Dictionary to modify
  @param    entry   Entry to delete (entry name)
  @return   void

  If the given entry can be found, it is deleted from the dictionary.
 */
/*--------------------------------------------------------------------------*/
void iniparser_unset(dictionary * ini, char * entry)
{
    dictionary_unset(ini, strlwc(entry));
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Load a single line from an INI file
  @param    input_line  Input line, may be concatenated multi-line input
  @param    section     Output space to store section
  @param    key         Output space to store key
  @param    value       Output space to store value
  @return   line_status value
 */
/*--------------------------------------------------------------------------*/
static line_status iniparser_line(
    char * input_line,
    char * section,
    char * key,
    char * value)
{   
    line_status sta ;
    char        line[ASCIILINESZ+1];
    int         len ;

    strcpy(line, strstrip(input_line));
    len = (int)strlen(line);

    sta = LINE_UNPROCESSED ;
    if (len<1) {
        /* Empty line */
        sta = LINE_EMPTY ;
    } else if (line[0]=='#') {
        /* Comment line */
        sta = LINE_COMMENT ; 
    } else if (line[0]=='[' && line[len-1]==']') {
        /* Section name */
        sscanf(line, "[%[^]]", section);
        strcpy(section, strstrip(section));
        strcpy(section, strlwc(section));
        sta = LINE_SECTION ;
    } else if (sscanf (line, "%[^=] = \"%[^\"]\"", key, value) == 2
           ||  sscanf (line, "%[^=] = '%[^\']'",   key, value) == 2
           ||  sscanf (line, "%[^=] = %[^;#]",     key, value) == 2) {
        /* Usual key=value, with or without comments */
        strcpy(key, strstrip(key));
        strcpy(key, strlwc(key));
        strcpy(value, strstrip(value));
        /*
         * sscanf cannot handle '' or "" as empty values
         * this is done here
         */
        if (!strcmp(value, "\"\"") || (!strcmp(value, "''"))) {
            value[0]=0 ;
        }
        sta = LINE_VALUE ;
    } else if (sscanf(line, "%[^=] = %[;#]", key, value)==2
           ||  sscanf(line, "%[^=] %[=]", key, value) == 2) {
        /*
         * Special cases:
         * key=
         * key=;
         * key=#
         */
        strcpy(key, strstrip(key));
        strcpy(key, strlwc(key));
        value[0]=0 ;
        sta = LINE_VALUE ;
    } else {
        /* Generate syntax error */
        sta = LINE_ERROR ;
    }
    return sta ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Parse an ini file and return an allocated dictionary object
  @param    ininame Name of the ini file to read.
  @return   Pointer to newly allocated dictionary

  This is the parser for ini files. This function is called, providing
  the name of the file to be read. It returns a dictionary object that
  should not be accessed directly, but through accessor functions
  instead.

  The returned dictionary must be freed using iniparser_freedict().
 */
/*--------------------------------------------------------------------------*/
dictionary * iniparser_load(const char * ininame)
{
    FILE * in ;

    char line    [ASCIILINESZ+1] ;
    char section [ASCIILINESZ+1] ;
    char key     [ASCIILINESZ+1] ;
    char tmp     [ASCIILINESZ+1] ;
    char val     [ASCIILINESZ+1] ;

    int  last=0 ;
    int  len ;
    int  lineno=0 ;
    int  errs=0;

    dictionary * dict ;

    if ((in=fopen(ininame, "r"))==NULL) {
        fprintf(stderr, "iniparser: cannot open %s\n", ininame);
        return NULL ;
    }

    dict = dictionary_new(0) ;
    if (!dict) {
        fclose(in);
        return NULL ;
    }

    memset(line,    0, ASCIILINESZ);
    memset(section, 0, ASCIILINESZ);
    memset(key,     0, ASCIILINESZ);
    memset(val,     0, ASCIILINESZ);
    last=0 ;

    while (fgets(line+last, ASCIILINESZ-last, in)!=NULL) {
        lineno++ ;
        len = (int)strlen(line)-1;
        /* Safety check against buffer overflows */
        if (line[len]!='\n') {
            fprintf(stderr,
                    "iniparser: input line too long in %s (%d)\n",
                    ininame,
                    lineno);
            dictionary_del(dict);
            fclose(in);
            return NULL ;
        }
        /* Get rid of \n and spaces at end of line */
        while ((len>=0) &&
                ((line[len]=='\n') || (isspace(line[len])))) {
            line[len]=0 ;
            len-- ;
        }
        /* Detect multi-line */
        if (line[len]=='\\') {
            /* Multi-line value */
            last=len ;
            continue ;
        } else {
            last=0 ;
        }
        switch (iniparser_line(line, section, key, val)) {
            case LINE_EMPTY:
            case LINE_COMMENT:
            break ;

            case LINE_SECTION:
            errs = dictionary_set(dict, section, NULL);
            break ;

            case LINE_VALUE:
            sprintf(tmp, "%s:%s", section, key);
            errs = dictionary_set(dict, tmp, val) ;
            break ;

            case LINE_ERROR:
            fprintf(stderr, "iniparser: syntax error in %s (%d):\n",
                    ininame,
                    lineno);
            fprintf(stderr, "-> %s\n", line);
            errs++ ;
            break;

            default:
            break ;
        }
        memset(line, 0, ASCIILINESZ);
        last=0;
        if (errs<0) {
            fprintf(stderr, "iniparser: memory allocation failure\n");
            break ;
        }
    }
    if (errs) {
        dictionary_del(dict);
        dict = NULL ;
    }
    fclose(in);
    return dict ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Free all memory associated to an ini dictionary
  @param    d Dictionary to free
  @return   void

  Free all memory associated to an ini dictionary.
  It is mandatory to call this function before the dictionary object
  gets out of the current context.
 */
/*--------------------------------------------------------------------------*/
void iniparser_freedict(dictionary * d)
{
    dictionary_del(d);
}

/* vim: set ts=4 et sw=4 tw=75 */
/********************************************************************/
/*-------------------------------------------------------------------------*/
/**
   @file	dictionary.c
   @author	N. Devillard
   @date	Sep 2007
   @version	$Revision: 1.27 $
   @brief	Implements a dictionary for string variables.

   This module implements a simple dictionary object, i.e. a list
   of string/string associations. This object is useful to store e.g.
   informations retrieved from a configuration file (ini files).
*/
/*--------------------------------------------------------------------------*/

/*
	$Id: dictionary.c,v 1.27 2007-11-23 21:39:18 ndevilla Exp $
	$Revision: 1.27 $
*/
/*---------------------------------------------------------------------------
   								Includes
 ---------------------------------------------------------------------------*/
#include "dictionary.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/** Maximum value size for integers and doubles. */
#define MAXVALSZ	1024

/** Minimal allocated number of entries in a dictionary */
#define DICTMINSZ	128

/** Invalid key token */
#define DICT_INVALID_KEY    ((char*)-1)

/*---------------------------------------------------------------------------
  							Private functions
 ---------------------------------------------------------------------------*/

/* Doubles the allocated size associated to a pointer */
/* 'size' is the current allocated size. */
static void * mem_double(void * ptr, int size)
{
    void * newptr ;
 
    newptr = calloc(2*size, 1);
    if (newptr==NULL) {
        return NULL ;
    }
    memcpy(newptr, ptr, size);
    free(ptr);
    return newptr ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Duplicate a string
  @param    s String to duplicate
  @return   Pointer to a newly allocated string, to be freed with free()

  This is a replacement for strdup(). This implementation is provided
  for systems that do not have it.
 */
/*--------------------------------------------------------------------------*/
static char * xstrdup(char * s)
{
    char * t ;
    if (!s)
        return NULL ;
    t = malloc(strlen(s)+1) ;
    if (t) {
        strcpy(t,s);
    }
    return t ;
}

/*---------------------------------------------------------------------------
  							Function codes
 ---------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/**
  @brief	Compute the hash key for a string.
  @param	key		Character string to use for key.
  @return	1 unsigned int on at least 32 bits.

  This hash function has been taken from an Article in Dr Dobbs Journal.
  This is normally a collision-free function, distributing keys evenly.
  The key is stored anyway in the struct so that collision can be avoided
  by comparing the key itself in last resort.
 */
/*--------------------------------------------------------------------------*/
unsigned dictionary_hash(char * key)
{
	int			len ;
	unsigned	hash ;
	int			i ;

	len = strlen(key);
	for (hash=0, i=0 ; i<len ; i++) {
		hash += (unsigned)key[i] ;
		hash += (hash<<10);
		hash ^= (hash>>6) ;
	}
	hash += (hash <<3);
	hash ^= (hash >>11);
	hash += (hash <<15);
	return hash ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Create a new dictionary object.
  @param	size	Optional initial size of the dictionary.
  @return	1 newly allocated dictionary objet.

  This function allocates a new dictionary object of given size and returns
  it. If you do not know in advance (roughly) the number of entries in the
  dictionary, give size=0.
 */
/*--------------------------------------------------------------------------*/
dictionary * dictionary_new(int size)
{
	dictionary	*	d ;

	/* If no size was specified, allocate space for DICTMINSZ */
	if (size<DICTMINSZ) size=DICTMINSZ ;

	if (!(d = (dictionary *)calloc(1, sizeof(dictionary)))) {
		return NULL;
	}
	d->size = size ;
	d->val  = (char **)calloc(size, sizeof(char*));
	d->key  = (char **)calloc(size, sizeof(char*));
	d->hash = (unsigned int *)calloc(size, sizeof(unsigned));
	return d ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Delete a dictionary object
  @param	d	dictionary object to deallocate.
  @return	void

  Deallocate a dictionary object and all memory associated to it.
 */
/*--------------------------------------------------------------------------*/
void dictionary_del(dictionary * d)
{
	int		i ;

	if (d==NULL) return ;
	for (i=0 ; i<d->size ; i++) {
		if (d->key[i]!=NULL)
			free(d->key[i]);
		if (d->val[i]!=NULL)
			free(d->val[i]);
	}
	free(d->val);
	free(d->key);
	free(d->hash);
	free(d);
	return ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Get a value from a dictionary.
  @param	d		dictionary object to search.
  @param	key		Key to look for in the dictionary.
  @param    def     Default value to return if key not found.
  @return	1 pointer to internally allocated character string.

  This function locates a key in a dictionary and returns a pointer to its
  value, or the passed 'def' pointer if no such key can be found in
  dictionary. The returned character pointer points to data internal to the
  dictionary object, you should not try to free it or modify it.
 */
/*--------------------------------------------------------------------------*/
char * dictionary_get(dictionary * d, char * key, char * def)
{
	unsigned	hash ;
	int			i ;

	hash = dictionary_hash(key);
	for (i=0 ; i<d->size ; i++) {
        if (d->key[i]==NULL)
            continue ;
        /* Compare hash */
		if (hash==d->hash[i]) {
            /* Compare string, to avoid hash collisions */
            if (!strcmp(key, d->key[i])) {
				return d->val[i] ;
			}
		}
	}
	return def ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Set a value in a dictionary.
  @param    d       dictionary object to modify.
  @param    key     Key to modify or add.
  @param    val     Value to add.
  @return   int     0 if Ok, anything else otherwise

  If the given key is found in the dictionary, the associated value is
  replaced by the provided one. If the key cannot be found in the
  dictionary, it is added to it.

  It is Ok to provide a NULL value for val, but NULL values for the dictionary
  or the key are considered as errors: the function will return immediately
  in such a case.

  Notice that if you dictionary_set a variable to NULL, a call to
  dictionary_get will return a NULL value: the variable will be found, and
  its value (NULL) is returned. In other words, setting the variable
  content to NULL is equivalent to deleting the variable from the
  dictionary. It is not possible (in this implementation) to have a key in
  the dictionary without value.

  This function returns non-zero in case of failure.
 */
/*--------------------------------------------------------------------------*/
int dictionary_set(dictionary * d, char * key, char * val)
{
	int			i ;
	unsigned	hash ;

	if (d==NULL || key==NULL) return -1 ;
	
	/* Compute hash for this key */
	hash = dictionary_hash(key) ;
	/* Find if value is already in dictionary */
	if (d->n>0) {
		for (i=0 ; i<d->size ; i++) {
            if (d->key[i]==NULL)
                continue ;
			if (hash==d->hash[i]) { /* Same hash value */
				if (!strcmp(key, d->key[i])) {	 /* Same key */
					/* Found a value: modify and return */
					if (d->val[i]!=NULL)
						free(d->val[i]);
                    d->val[i] = val ? xstrdup(val) : NULL ;
                    /* Value has been modified: return */
					return 0 ;
				}
			}
		}
	}
	/* Add a new value */
	/* See if dictionary needs to grow */
	if (d->n==d->size) {

		/* Reached maximum size: reallocate dictionary */
		d->val  = (char **)mem_double(d->val,  d->size * sizeof(char*)) ;
		d->key  = (char **)mem_double(d->key,  d->size * sizeof(char*)) ;
		d->hash = (unsigned int *)mem_double(d->hash, d->size * sizeof(unsigned)) ;
        if ((d->val==NULL) || (d->key==NULL) || (d->hash==NULL)) {
            /* Cannot grow dictionary */
            return -1 ;
        }
		/* Double size */
		d->size *= 2 ;
	}

    /* Insert key in the first empty slot */
    for (i=0 ; i<d->size ; i++) {
        if (d->key[i]==NULL) {
            /* Add key here */
            break ;
        }
    }
	/* Copy key */
	d->key[i]  = xstrdup(key);
    d->val[i]  = val ? xstrdup(val) : NULL ;
	d->hash[i] = hash;
	d->n ++ ;
	return 0 ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Delete a key in a dictionary
  @param	d		dictionary object to modify.
  @param	key		Key to remove.
  @return   void

  This function deletes a key in a dictionary. Nothing is done if the
  key cannot be found.
 */
/*--------------------------------------------------------------------------*/
void dictionary_unset(dictionary * d, char * key)
{
	unsigned	hash ;
	int			i ;

	if (key == NULL) {
		return;
	}

	hash = dictionary_hash(key);
	for (i=0 ; i<d->size ; i++) {
        if (d->key[i]==NULL)
            continue ;
        /* Compare hash */
		if (hash==d->hash[i]) {
            /* Compare string, to avoid hash collisions */
            if (!strcmp(key, d->key[i])) {
                /* Found key */
                break ;
			}
		}
	}
    if (i>=d->size)
        /* Key not found */
        return ;

    free(d->key[i]);
    d->key[i] = NULL ;
    if (d->val[i]!=NULL) {
        free(d->val[i]);
        d->val[i] = NULL ;
    }
    d->hash[i] = 0 ;
    d->n -- ;
    return ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Dump a dictionary to an opened file pointer.
  @param	d	Dictionary to dump
  @param	f	Opened file pointer.
  @return	void

  Dumps a dictionary onto an opened file pointer. Key pairs are printed out
  as @c [Key]=[Value], one per line. It is Ok to provide stdout or stderr as
  output file pointers.
 */
/*--------------------------------------------------------------------------*/
void dictionary_dump(dictionary * d, FILE * out)
{
	int		i ;

	if (d==NULL || out==NULL) return ;
	if (d->n<1) {
		fprintf(out, "empty dictionary\n");
		return ;
	}
	for (i=0 ; i<d->size ; i++) {
        if (d->key[i]) {
            fprintf(out, "%20s\t[%s]\n",
                    d->key[i],
                    d->val[i] ? d->val[i] : "UNDEF");
        }
	}
	return ;
}


/* Test code */
#ifdef TESTDIC
#define NVALS 20000
int main(int argc, char *argv[])
{
	dictionary	*	d ;
	char	*	val ;
	int			i ;
	char		cval[90] ;

	/* Allocate dictionary */
	printf("allocating...\n");
	d = dictionary_new(0);
	
	/* Set values in dictionary */
	printf("setting %d values...\n", NVALS);
	for (i=0 ; i<NVALS ; i++) {
		sprintf(cval, "%04d", i);
		dictionary_set(d, cval, "salut");
	}
	printf("getting %d values...\n", NVALS);
	for (i=0 ; i<NVALS ; i++) {
		sprintf(cval, "%04d", i);
		val = dictionary_get(d, cval, DICT_INVALID_KEY);
		if (val==DICT_INVALID_KEY) {
			printf("cannot get value for key [%s]\n", cval);
		}
	}
    printf("unsetting %d values...\n", NVALS);
	for (i=0 ; i<NVALS ; i++) {
		sprintf(cval, "%04d", i);
		dictionary_unset(d, cval);
	}
    if (d->n != 0) {
        printf("error deleting values\n");
    }
	printf("deallocating...\n");
	dictionary_del(d);
	return 0 ;
}
#endif

/* vim: set ts=4 et sw=4 tw=75 */
/********************************************************************/
