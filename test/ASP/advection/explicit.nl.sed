&ctl_nl
NThreads      = 8
partmethod    = 4
topology      = "cube"
test_case     = "asp_tracer"
ne            = NE
qsize         = 2
tracer_advection_formulation = 1
ndays          = 12
rotate_grid   = 0
statefreq     = SFREQ
accumfreq     = -1
accumstart    = 200
accumstop     = 1200
restartfreq   = 43200
restartfile   = "./R0001"
runtype       = 0
tstep         = TSTEP
integration   = "explicit"
smooth        = 0.00               ! default = 0.05
nu            = NU1
nu_s          = NU2
nu_p          = 0
limiter_option = 4
hypervis_order = 2
hypervis_subcycle= 1
prescribed_wind = 1
energy_fixer = -1
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
filter_freq_advection   = 0
filter_mu_advection   = 0.00
p_bv          = 12.0D0
s_bv          = .80
wght_fm       = 0.10D0
kcut_fm       = 2
/
&vert_nl
vform         = "ccm"
vfile_mid     = "vcoord/aspL60_mid.ascii"
vfile_int     = "vcoord/aspL60_int.ascii"
/
&analysis_XX
 interp_gridtype=1
 output_timeunits=0
 output_frequency=0
/
&analysis_nl
 interp_gridtype=1
! every 6h
! output_timeunits=2
! output_frequency=6
 output_timeunits=1
 output_frequency=2
 output_start_time=0
 output_end_time=28800
 output_varnames1='Q','Q2'
 interp_type=0
 io_stride=8
 output_type = 'netcdf'
/
