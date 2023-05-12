&ctl_nl
NThreads          = -1
partmethod        = 4
topology          = "cube"
test_case         = "jw_baroclinic"
u_perturb         = 1
rotate_grid       = 0
ne                = 1024
qsize             = 10
nmax              = 4096
disable_diagnostics = .true.
statefreq         = 99999999
restartfreq       = 43200
restartfile       = "./R0001"
runtype           = 0
mesh_file         = '/dev/null'
tstep             = 10.0
integration       = "explicit"
smooth            = 0
nu                = 2.50e+10
nu_div            = 2.50e+10
nu_p              = 2.50e+10
nu_q              = 2.50e+10
nu_s              =-1
nu_top            = 0
se_ftype          = 0
limiter_option    = 9
vert_remap_q_alg  = 10
hypervis_scaling  = 0
hypervis_order    = 2
hypervis_subcycle = 1
hypervis_subcycle_tom  = 0
theta_hydrostatic_mode = false
theta_advect_form = 1
tstep_type        = 10
dt_remap_factor     = 2
dt_tracer_factor    = 8
hypervis_subcycle_q = 8
transport_alg       = 12
semi_lagrange_hv_q  = 1
semi_lagrange_cdr_check = .false.
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
vfile_mid = '/gpfs/alpine/scratch/ambradl/cli115/vcoord/sabm-128.ascii'
vfile_int = '/gpfs/alpine/scratch/ambradl/cli115/vcoord/sabi-128.ascii'
/

&prof_inparm
profile_outpe_num   = 100
profile_single_file = .true.
/

!  timunits: 0= steps, 1=days, 2=hours
&analysis_nl
 interp_gridtype   = 2
 output_timeunits  = 1,1
 output_frequency  = 0,0
 output_start_time = 0,0
 output_end_time   = 30000,30000
 output_varnames1  = 'ps','zeta','u','v','T'
 output_varnames2  = 'Q','Q2','Q3','Q4'
 io_stride         = 8
 output_type       = 'netcdf'
/
