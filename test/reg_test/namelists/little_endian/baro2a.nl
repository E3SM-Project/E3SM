&ctl_nl
NThreads                     = 1
partmethod                   = 4
topology                     = "cube"
test_case                    = "asp_baroclinic"
rotate_grid                  = 0
ne                           = 15
qsize                        = 4
tracer_advection_formulation = 1
tstep_type                   = 0
compute_mean_flux            = 0  
ndays                        = 9
statefreq                    = 480
accumfreq                    = -1
accumstart                   = 200
accumstop                    = 1200
restartfreq                  = 43200
restartfile                  = "./R0001"
runtype                      = 0
tstep                        = 180
qsplit                       = 1
rk_stage_user                = 0
psurf_vis                    = 0  
integration                  = "explicit"
smooth                       = 0.005
nu                           = 1e16
nu_s                         = -1        ! use same value as nu
nu_q                         = 1e16    
nu_p                         = 0
limiter_option               = 8 
energy_fixer                 = -1
hypervis_order               = 2
hypervis_subcycle            = -1
u_perturb                    = 1
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
vfile_mid     = "vcoord/camm-26.fbin.littleendian"
vfile_int     = "vcoord/cami-26.fbin.littleendian"
/

&prof_inparm
profile_outpe_num   = 100
profile_single_file = .true.
/

&analysis_nl
output_prefix     = "baro2a-"
interp_gridtype   = 2
output_timeunits  = 1,1
output_frequency  = 9,9
output_start_time = 0,0
output_end_time   = 30,30
output_varnames1  = 'u', 'v', 'ps', 'T', 'zeta'
output_varnames2  = 'Q', 'Q2', 'Q3', 'Q4'
io_stride         = 8
output_type       = 'netcdf' 
/

