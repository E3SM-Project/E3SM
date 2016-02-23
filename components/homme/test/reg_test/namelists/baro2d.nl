&ctl_nl
NThreads                     = 1
partmethod                   = 4
topology                     = "cube"
test_case                    = "asp_baroclinic"
rotate_grid                  = 0
ne                           = 15
qsize                        = 4
ntrac                        = 0
tstep_type                   = 1 
ndays                        = 9
statefreq                    = 576
restartfreq                  = 43200
restartfile                  = "./R0001"
runtype                      = 0
tstep                        = 150
qsplit                       = 4
rk_stage_user                = 3
integration                  = "explicit"
nu                           = 1e16
nu_div                       = 2.5e16
nu_s                         = -1        ! use same value as nu
nu_q                         = 1e16    
nu_p                         = 1e16
limiter_option               = 8 
energy_fixer                 = -1
hypervis_order               = 2
hypervis_subcycle            = 2
hypervis_subcycle_q          = 1
u_perturb                    = 1
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
output_prefix     = "baro2d-"
interp_gridtype   = 2
output_timeunits  = 1,1,1
output_frequency  = 9,9,9
output_start_time = 0,0,0
output_end_time   = 30,30,30
output_varnames1  = 'u', 'v', 'ps', 'T', 'zeta'
output_varnames2  = 'Q', 'Q2', 'Q3', 'Q4'
output_varnames3  = 'phys_lat','phys_lon'
io_stride         = 8
output_type       = 'netcdf' 
/

