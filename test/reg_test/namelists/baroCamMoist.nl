&ctl_nl
NThreads                     = 1
partmethod                   = 4
topology                     = "cube"
test_case                    = "asp_baroclinic"
rotate_grid                  = 0
ne                           = 6
qsize                        = 25
tstep_type                   = 5 
ndays                        = 3
statefreq                    = 45
restartfreq                  = 43200
restartfile                  = "./R0001"
runtype                      = 0
tstep                        = 480
rsplit                       = 3
qsplit                       = 1
psurf_vis                    = 0  
integration                  = "explicit"
smooth                       = 0
nu                           = 5e16
nu_s                         = -1  ! use same value as nu
nu_q                         = 5e16  
nu_p                         = 5e16
nu_div                       = -1
npdg=0
limiter_option               = 8
energy_fixer                 = -1
hypervis_order               = 2
hypervis_subcycle            = 4
u_perturb                    = 1
vert_remap_q_alg = 1
tracer_advection_formulation = 1
! disable_diagnostics          = .true.
moisture = 'notdry'
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
vfile_mid     = "vcoord/camm-26.ascii"
vfile_int     = "vcoord/cami-26.ascii"
/

&prof_inparm
profile_outpe_num   = 100
profile_single_file	= .true.
/

&analysis_nl
output_prefix     = "camBaroMoist-"
interp_gridtype   = 2
output_timeunits  = 1,1
output_frequency  = 3,3
output_start_time = 0,0
output_end_time   = 30,30
output_varnames1  = 'zeta', 'u', 'v', 'ps', 'dp3d'
output_varnames2  = 'Q', 'Q2', 'Q3', 'Q4','phys_lat','phys_lon'
io_stride         = 8
output_type       = 'netcdf' 
/

