&ctl_nl
theta_hydrostatic_mode = .false.
dcmip4_moist  = 0
dcmip4_X      = 100.0
vert_num_threads = 1
NThreads=1
partmethod    = 4
topology      = "cube"
test_case     = "dcmip2012_test4"
u_perturb = 1
rotate_grid = 0
ne=11
qsize = 0
nmax = 38880
statefreq=1080
restartfreq =   38880
restartfile   = "./restart/R000038880"
runtype       = 0
mesh_file='/dev/null'
tstep=.2
rsplit=0
qsplit = 1
tstep_type = 5
integration   = "explicit"
nu=2e10
nu_p=2e10
nu_q=2e10
nu_s=2e10
nu_top = 0
se_ftype     = 0
limiter_option = 8
vert_remap_q_alg = 1
hypervis_scaling=0
hypervis_order = 2
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

 output_timeunits=1              ! 1- days, 2 hours, 0 - tsteps
 output_frequency=1
 output_start_time=9
 output_end_time=3000
 output_varnames1='ps','zeta','u','v','T'
 num_io_procs      = 16
 output_type = 'netcdf'
 output_prefix = 'hydro-X100-erun1-'
/





