&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "jw_baroclinic"
ne            = NE
qsize = 0
nmax          = 1
statefreq     = 1
accumfreq     = -1
accumstart    = 200
accumstop     = 1200
restartfreq   = 43200
restartfile   = "./R0001"
runtype       = -1
tstep         = TSTEP
integration   = "explicit"
smooth        = 0.005           
nu            = 0
nu_s          = 0
nu_p          = 0
hypervis_order = 2
hypervis_subcycle = 1
u_perturb      = 1
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
/
&vert_nl
vform         = "ccm"
vfile_mid     = "/home/mataylo/codes/homme/test/vcoord/sabm-20.fbin.littleendian"
vfile_int     = "/home/mataylo/codes/homme/test/vcoord/sabi-20.fbin.littleendian"
/
&analysis_nl
 output_timeunits=1
 output_frequency=0
 infilenames='h1-tavg.nc'
 output_varnames1='psxx'
! num_io_procs = 1 
 io_stride = 8
 output_type = 'netcdf'
 interp_gridtype = 2 
 interp_type = 0
 interp_nlat=256
 interp_nlon=512
/

