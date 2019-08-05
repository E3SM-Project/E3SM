&ctl_nl
NThreads=-1
partmethod    = 4
topology      = "cube"
test_case     = "jw_baroclinic"
u_perturb = 1
rotate_grid = 0
ne=1024
qsize = 10
! tstep=9  200 steps = 30min  (Bechmark reports 30min time)
! tstep=10 180 steps = 30min  (Bechmark reports 30min time)
nmax = 180
statefreq=90
disable_diagnostics = .true.
theta_hydrostatic_mode=.false.
theta_advect_form = 1
runtype       = 0
tstep=10
rsplit=8
qsplit = 1
tstep_type = 7
integration   = "explicit"
nu=2.5e10
nu_div=2.5e10
nu_p=2.5e10
nu_q=2.5e10
nu_s=2.5e10
nu_top = 0e5
se_ftype     = 0
limiter_option = 9
vert_remap_q_alg = 1
hypervis_order = 2
hypervis_subcycle=1
hypervis_subcycle_q=1
/
&vert_nl
vform         = "ccm"
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
 output_type = 'pnetcdf64' 
/

