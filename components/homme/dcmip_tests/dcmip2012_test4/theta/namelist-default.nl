&ctl_nl
theta_hydrostatic_mode = .true.
dcmip4_moist  = 0
dcmip4_X      = 1000.0
vert_num_threads = 1
NThreads=1
partmethod    = 4
topology      = "cube"
test_case     = "dcmip2012_test4"
u_perturb = 1
rotate_grid = 0
ne=30
qsize = 1
nmax = 6480
ndays          = 0
statefreq=60
runtype       = 0
mesh_file='/dev/null'
tstep=0.2
rsplit=0
qsplit = 1
tstep_type = 5
energy_fixer  = -1
integration   = "explicit"
smooth        = 0
nu=1e6
nu_div=1e6
nu_p=1e6
nu_q=1e6
nu_s=1e6
nu_top = 0
se_ftype     = 0
limiter_option = 8
vert_remap_q_alg = 1
hypervis_scaling=0
hypervis_order = 2
hypervis_subcycle=3
hypervis_subcycle=3
/
&solver_nl
precon_method = "identity"
maxits        = 500
tol           = 1.e-9
/
&filter_nl
filter_type   = "taylor"
transfer_type = "bv"
filter_freq   = 0
filter_mu     = 0.04D0
p_bv          = 12.0D0
s_bv          = .666666666666666666D0
wght_fm       = 0.10D0
kcut_fm       = 2
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

 output_timeunits=0              ! 1 is days, 2 is hours?
 output_frequency=432
 output_start_time=3024
 output_end_time=300000000
 output_varnames1='ps','zeta','u','v','T'
 io_stride=8
 output_type = 'netcdf'
/





