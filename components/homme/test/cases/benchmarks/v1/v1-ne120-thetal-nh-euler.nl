&ctl_nl
!-------------------------------------EULER
rsplit=2                             !NOTE THAT RSPLIT HAS TO DIVIDE # OF TSTEPS/NMAX
qsplit=1                             !NOTE THAT QSPLIT HAS TO DIVIDE # OF TSTEPS/NMAX
se_ftype=0
!-------------------------------------HYDRO/NON
!theta_hydrostatic_mode=.true.      
!tstep_type    = 5                   
theta_hydrostatic_mode = .false.   
tstep_type    = 7                  
!-------------------------------------RES,RUN,OUTPUT
ne=120
qsize = 40
ndays=1
statefreq=999999999
disable_diagnostics = .true.
tstep=75      ! ne30: 300  ne120: 75
!-------------------------------------HV
nu=1e13
nu_div=2.5e13
nu_p=1e13
nu_q=1e13
nu_s=1e13
nu_top = 2.5e5
hypervis_scaling=0
hypervis_order = 2
hypervis_subcycle=4
!-------------------------------------UNLIKELY TO CHANGE
NThreads=-1
partmethod    = 4
topology      = "cube"
test_case     = "jw_baroclinic"
u_perturb = 1
rotate_grid = 0
restartfreq   = 999999999
restartfile   = "./R0001"
runtype       = 0
mesh_file='/dev/null'
integration   = "explicit"
limiter_option = 9 ! this is diff from what's for preqx
vert_remap_q_alg = 1
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

