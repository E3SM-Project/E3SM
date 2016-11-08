&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "jw_baroclinic"
ne            = 6
qsize = 1
nmax          = 864   
statefreq     = 60
restartfreq =   864   
restartfile   = "./R0001"
runtype       = 0
tstep         = 900 
tstep_type    = 5
integration   = "explicit"
smooth        = 0.0000
energy_fixer  = 0
nu            = 4e16
nu_s          = 4e16
nu_p          = 0
nu_top        = 0
hypervis_order = 2
hypervis_subcycle = 1
u_perturb      = 1
initial_total_mass = 0
limiter_option = 8
/
&filter_nl
filter_freq   = 0
/
&vert_nl
vform         = "ccm"
vfile_mid     = "vcoord/camm-26.ascii"
vfile_int     = "vcoord/cami-26.ascii"
/
&analysis_nl
output_prefix    = "baro2c-run1-"
output_timeunits = 1
output_frequency = 3
output_varnames1 = 'u', 'v', 'ps', 'T'
interp_type      = 0          ! native high order
output_type      = 'netcdf'
io_stride        = 8
/



