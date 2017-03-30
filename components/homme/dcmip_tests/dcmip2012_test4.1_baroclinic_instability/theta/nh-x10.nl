&ctl_nl
theta_hydrostatic_mode = .false.
dcmip4_moist  = 0
dcmip4_X      = 10.0
vert_num_threads = 1
NThreads=1
partmethod    = 4
topology      = "cube"
test_case     = "dcmip2012_test4"
u_perturb = 1
rotate_grid = 0
ne=30
qsize = 0
nmax = 648000
statefreq=6000
runtype       = 0
mesh_file='/dev/null'
tstep=0.2
rsplit=8
qsplit = 1
tstep_type = 5
energy_fixer  = -1
integration   = "explicit"
smooth        = 0
nu=1e12
nu_div=1e12
nu_p=1e12
nu_q=1e12
nu_s=1e12
nu_top = 0
se_ftype     = 0
limiter_option = 8
vert_remap_q_alg = 1
hypervis_scaling=0
hypervis_order = 2
hypervis_subcycle=3
hypervis_subcycle=3
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

 output_timeunits=0              ! 1 is days, 2- hours, 0 - tsteps
 output_frequency=43200
 output_start_time=302400
 output_end_time=300000000
 output_varnames1='ps','zeta','u','v','T'
 num_io_procs      = 16
 output_type = 'netcdf'
 output_prefix= 'nonhydro-X10-'
/





