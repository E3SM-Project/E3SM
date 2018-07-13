&ctl_nl
theta_hydrostatic_mode = .true.
dcmip4_moist  = 0
dcmip4_X      = 1.0
NThreads=1
partmethod    = 4
topology      = "cube"
test_case     = "dcmip2012_test4"
u_perturb = 1
rotate_grid = 0
ne=30
qsize = 0
ndays = 15
statefreq=288
runtype       = 0
mesh_file='/dev/null'
tstep=300
rsplit=3
qsplit = 1
tstep_type = 5
integration   = "explicit"
nu=1e15
nu_div=1e15
nu_p=1e15
nu_q=1e15
nu_s=1e15
nu_top = 0
se_ftype     = 0
limiter_option = 9
vert_remap_q_alg = 0
hypervis_scaling=0
hypervis_order = 2
hypervis_subcycle = 1
hypervis_subcycle = 1
/
&vert_nl
vform         = "ccm"
vfile_mid = '../vcoord/camm-30.ascii'
vfile_int = '../vcoord/cami-30.ascii'
/

&prof_inparm
profile_outpe_num = 100
profile_single_file             = .true.
/

&analysis_nl
! to compare with EUL ref solution:
! interp_nlat = 512
! interp_nlon = 1024
 interp_gridtype=2

 output_timeunits=1              ! 1- days, 2 hours, 0 - tsteps
 output_frequency=1
 output_start_time=7
 output_end_time=3000
 output_varnames1='ps','zeta','u','v','T'
 num_io_procs      = 16
 output_type = 'netcdf'
 output_prefix = 'preqx-X1-'
/





