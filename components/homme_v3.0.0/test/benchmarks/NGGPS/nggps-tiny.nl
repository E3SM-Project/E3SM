&ctl_nl
NThreads=-1   ! use maximum threads available
partmethod    = 4
topology      = "cube"
test_case     = "jw_baroclinic"
u_perturb = 1
rotate_grid = 0
ne=4
qsize = 10
! tstep=40  360 steps = 4h  (Bechmark reports 2h time)
nmax = 360
statefreq=360
restartfreq   = 43200
restartfile   = "./R0001"
theta_hydrostatic_mode=.true.
runtype       = 0
tstep=40
rsplit=3
qsplit = 1
tstep_type = 5
integration   = "explicit"
nu=4.5e17
nu_div=4.5e17
nu_p=4.5e17
nu_q=4.5e17
nu_s=4.5e17
nu_top = 0e5
se_ftype     = 0
limiter_option = 9
vert_remap_q_alg = 1
hypervis_order = 2
hypervis_subcycle=1
hypervis_subcycle_q=1
/
&vert_nl
vfile_mid = './sabm-128.ascii'
vfile_int = './sabi-128.ascii'
/

&prof_inparm
profile_outpe_num = 100
profile_single_file		= .true.
/

&analysis_nl
! to compare with EUL ref solution:
! interp_nlat = 32
! interp_nlon = 64
! interp_gridtype=2
 
 output_timeunits=1,1
 output_frequency=0,0
 output_start_time=0,0
 output_end_time=30000,30000
 output_varnames1='ps','zeta','dp3d'
 output_varnames2='Q','Q2','Q3','Q4','Q5'
 io_stride=8
 output_type = 'netcdf' 
/

