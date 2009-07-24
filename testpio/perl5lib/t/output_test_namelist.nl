&cam_inparm
 absems_data		= '/fs/cgd/csm/inputdata/atm/cam/rad/abs_ems_factors_fastvx.c030508.nc'
 bnd_topo		= '/fs/cgd/csm/inputdata/atm/cam/topo/USGS-gtopo30_10x15_remap_c050520.nc'
 complex1		= (-5.2345d+23,+5.98765E-100), (5.2345,5.98765E-33)
 diag_cnst_conv_tend		= 'q_only'
 dtime		= 1800
 ncdata		= '/fs/cgd/csm/inputdata/atm/cam/inic/fv/cami_0000-01-01_10x15_L26_c030918.nc'
 phys_alltoall		= 13
 scenario_prognostic_sulfur		= 'RAMPED'
/
&clm_inparm
 fatmgrid		= '/fs/cgd/csm/inputdata/lnd/clm2/griddata/griddata_10x15_USGS_070110.nc'
 fatmlndfrc		= '/fs/cgd/csm/inputdata/lnd/clm2/griddata/fracdata_10x15_USGS_070110.nc'
 fpftcon		= '/fs/cgd/csm/inputdata/lnd/clm2/pftdata/pft-physiology.c070207'
 fsurdat		= '/fs/cgd/csm/inputdata/lnd/clm2/surfdata/surfdata_10x15_USGS_070110.nc'
 hist_fincl1		= 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy',
         'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy',
         'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy',
         'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy',
         'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy',
         'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy',
         'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy',
         'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy', 'thingy',
         'thingy', 'thingy', 'thingy'
/
&generic_namelist
 stop_option		= 'nyears'
/
&seq_infodata_inparm
 case_name		= 'case_name_in_quotes'
 orb_iyear_ad		= 1990
 start_type		= 'startup'
/
&seq_timemgr_inparm
 atm_cpl_dt		= 3600
 ice_cpl_dt		= 3600
 lnd_cpl_dt		= 3600
 ocn_cpl_dt		= 3600
 restart_option		= 'monthly'
 start_ymd		= 20000101
 stop_n		= 1
 stop_option		= 'ndays'
/
#!--------------------------------------------------------------------------------------------------------------------------
#! output_test_namelist.nl.tmp:: Add note to end
#!--------------------------------------------------------------------------------------------------------------------------
