&ctl_nl
nthreads      = 1
partmethod    = 4
topology      = "plane"
geometry      = "plane"
cubed_sphere_map = 2
test_case     = "planar_dbl_vrtx"
ne_x            = 20
ne_y            = 20
nmax          = 250
statefreq     = 250
restartfreq   = -1
runtype       = 0
tstep         = 100.0
integration   = "explicit"
smooth        = .0025
hypervis_order = 2
hypervis_subcycle = 1
hypervis_scaling=3.0
nu=0.216784
/
&filter_nl
transfer_type = "bv"
filter_type   = "taylor"
filter_freq   = 0
filter_mu     = 0.05D0
p_bv          = 12.0D0
s_bv          = .666666666666666666D0
kcut_fm       = 2
wght_fm       = 0.10D0
/
&analysis_nl
output_timeunits=0,
output_frequency=250,
output_dir        = "./movies/"               ! destination dir for netcdf file
interp_type       = 0                         ! 0=native grid, 1=bilinear
output_type       ='netcdf'                   ! netcdf or pnetcdf
/
