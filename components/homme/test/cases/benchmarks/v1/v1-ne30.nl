&ctl_nl
NThreads=-1
partmethod    = 4
topology      = "cube"
test_case     = "jw_baroclinic"
u_perturb = 1
rotate_grid = 0
ne=30
ndays=1
qsize = 40
statefreq=360
disable_diagnostics = .true.
restartfreq   = 43200
restartfile   = "./R0001"
runtype       = 0
mesh_file='/dev/null'
tstep=300      ! ne30: 300  ne120: 75
rsplit=3       ! ne30: 3   ne120:  2
qsplit = 1
tstep_type = 5
integration   = "explicit"
theta_hydrostatic_mode=.false.
nu=1e15
nu_div=2.5e15
nu_p=1e15
nu_q=1e15
nu_s=1e15
nu_top = 2.5e5
se_ftype     = 0
limiter_option = 8
vert_remap_q_alg = 1
hypervis_scaling=0
hypervis_order = 2
hypervis_subcycle=3    ! ne30: 3  ne120: 4
/
&vert_nl
vfile_mid = './acme-72m.ascii'
vfile_int = './acme-72i.ascii'
/

&prof_inparm
profile_outpe_num = 100
profile_single_file		= .true.
/

&analysis_nl
! disabled
 output_timeunits=1,1
 output_frequency=0,0
 output_start_time=0,0
 output_end_time=30000,30000
 output_varnames1='ps','zeta','T','geo'
 output_varnames2='Q','Q2','Q3','Q4','Q5'
 io_stride=8
 output_type = 'netcdf' 
/

