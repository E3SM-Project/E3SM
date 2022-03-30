&ctl_nl
NThreads=-1
partmethod    = 4
topology      = "cube"
test_case     = "jw_baroclinic"
u_perturb = 1
rotate_grid = 0
ne=256
qsize = 10
! run 4h.   (Bechmark reports 2h time)
nmax = 360
!ndays=2
statefreq=200
disable_diagnostics = .true.
!disable_diagnostics = .false.
runtype       = 0
theta_hydrostatic_mode=.false.
theta_advect_form = 1
tstep=40
rsplit=8
qsplit = 1
tstep_type = 7
integration   = "explicit"
nu=1.6e12
nu_div=1.6e12
nu_p=1.6e12
nu_q=1.6e12
nu_s=1.6e12
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
 output_varnames1='ps','zeta'
 output_varnames2='Q','Q2','Q3','Q4','Q5'
 io_stride=8
 output_type = 'pnetcdf' 
/

