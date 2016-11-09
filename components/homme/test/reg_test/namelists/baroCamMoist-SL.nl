&ctl_nl
NThreads                     = 1
partmethod                   = 4
topology                     = "cube"
test_case                    = "asp_baroclinic"
rotate_grid                  = 0
ne                           = 2
qsize                        = 4
tstep_type                   = 5 
ndays                        = 1
statefreq                    = 36
restartfreq                  = 43200
restartfile                  = "./R0001"
runtype                      = 0
tstep                        = 1200
rsplit                       = 3
qsplit                       = 3
use_semi_lagrange_transport  = .true.
cubed_sphere_map = 2
integration                  = "explicit"
nu                           = 5e16
nu_s                         = -1  ! use same value as nu
nu_q                         = 5e16  
nu_p                         = 5e16
nu_div                       = -1
limiter_option               = 8
energy_fixer                 = -1
hypervis_order               = 2
hypervis_subcycle            = 4
u_perturb                    = 1
vert_remap_q_alg = 1
! disable_diagnostics          = .true.
moisture = 'notdry'
/


&filter_nl
filter_freq   = 0
/

&vert_nl
vform         = "ccm"
vfile_mid     = "vcoord/camm-26.ascii"
vfile_int     = "vcoord/cami-26.ascii"
/

&prof_inparm
profile_outpe_num   = 100
profile_single_file	= .true.
/

&analysis_nl
output_prefix     = "camBaroMoist-SL-"
interp_gridtype   = 2
output_timeunits  = 1,1
output_frequency  = 1,1
output_start_time = 0,0
output_end_time   = 30,30
output_varnames1  = 'zeta', 'u', 'v', 'ps', 'dp3d'
output_varnames2  = 'Q', 'Q2', 'Q3', 'Q4','phys_lat','phys_lon'
io_stride         = 8
output_type       = 'netcdf' 
/

