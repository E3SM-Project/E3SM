&cam_inparm
 absems_data = '/fs/cgd/csm/inputdata/atm/cam/rad/abs_ems_factors_fastvx.c030508.nc'
 bnd_topo = '/fs/cgd/csm/inputdata/atm/cam/topo/USGS-gtopo30_10x15_remap_c050520.nc'
 complex1 = (-5.2345d+23,+5.98765E-100), (5.2345,5.98765E-33)
 diag_cnst_conv_tend = 'q_only'
 dtime = 1800
 ncdata = '/fs/cgd/csm/inputdata/atm/cam/inic/fv/cami_0000-01-01_10x15_L26_c030918.nc'
 phys_alltoall = 13
 scenario_prognostic_sulfur = 'RAMPED'
/
&clm_inparm
 fatmgrid = '/fs/cgd/csm/inputdata/lnd/clm2/griddata/griddata_10x15_USGS_070110.nc'
 fatmlndfrc = '/fs/cgd/csm/inputdata/lnd/clm2/griddata/fracdata_10x15_USGS_070110.nc'
 fpftcon = '/fs/cgd/csm/inputdata/lnd/clm2/pftdata/pft-physiology.c070207'
 fsurdat = '/fs/cgd/csm/inputdata/lnd/clm2/surfdata/surfdata_10x15_USGS_070110.nc'
 hist_fincl1 = 'thingfirst', 'thing0', 'thing1', 'thing2', 'thing3', 'thing4', 'thing5', 'thing6', 'thing7',
         'thing8', 'thing9', 'thing10', 'thing11', 'thing12', 'thing13', 'thing14', 'thing15', 'thing16',
         'thing17', 'thing18', 'thing19', 'thing20', 'thing21', 'thing22', 'thing23', 'thing24', 'thing25',
         'thing26', 'thing27', 'thing28', 'thing29', 'thing30', 'thing31', 'thing32', 'thing33', 'thing34',
         'thing35', 'thing36', 'thing37', 'thing38', 'thing39', 'thing40', 'thing41', 'thing42', 'thing43',
         'thing44', 'thing45', 'thing46', 'thing47', 'thing48', 'thing49', 'thing50', 'thing51', 'thing52',
         'thing53', 'thing54', 'thing55', 'thing56', 'thing57', 'thing58', 'thing59', 'thing60', 'thing61',
         'thing62', 'thing63', 'thing64', 'thing65', 'thing66', 'thing67', 'thing68', 'thing69', 'thing70',
         'thing71', 'thing72', 'thing73', 'thing74', 'thing75', 'thing76', 'thing77', 'thing78', 'thing79',
         'thing80', 'thing81', 'thing82', 'thing83', 'thing84', 'thing85', 'thing86', 'thing87', 'thing88',
         'thing89', 'thing90', 'thing91', 'thing92', 'thing93', 'thing94', 'thing95', 'thing96', 'thing97',
         'thing98'
 hist_fincl2 = ''
/
&generic_namelist
 stop_option = 'nyears'
/
&seq_infodata_inparm
 case_name = 'case_name_in_quotes'
 orb_eccen = 1.e-12
 orb_iyear_ad = 1990
 start_type = 'startup'
/
&seq_timemgr_inparm
 atm_cpl_dt = 3600
 ice_cpl_dt = 3600
 lnd_cpl_dt = 3600
 ocn_cpl_dt = 3600
 restart_option = 'monthly'
 start_ymd = 20000101
 stop_n = 1
 stop_option = 'ndays'
/
#!--------------------------------------------------------------------------------------------------------------------------
#! output_test_namelist.nl.tmp:: Add note to end
#!--------------------------------------------------------------------------------------------------------------------------
