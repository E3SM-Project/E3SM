&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "held_suarez0"
ne            = NE
meshfile      = "not used"
qsize = 0
nmax          = 0
statefreq     = 1
restartfreq   = 43200
restartfile   = "./R0001"
runtype       = 0
tstep         = 0.0001
tstep_type   = 1
qsplit=1
integration   = "explicit"
nu=0
nu_s          = 0
nu_p          = 0
smooth_phis_numcycle = 0
smooth_phis_nudt = 0
hypervis_scaling = 0 
hypervis_order = 2
hypervis_subcycle= 1
hypervis_subcycle_q= -1
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
vfile_mid     = "../vcoord/aspL20_mid.isotherm.ascii"
vfile_int     = "../vcoord/aspL20_int.isotherm.ascii"
/
&analysis_nl
 output_timeunits=1
 output_frequency=1
 infilenames='h1-tavg.nc'
 output_varnames1='geos'
 output_type='netcdf'
! num_io_procs = 1
 io_stride = 16
/

&prof_inparm
profile_outpe_num = 8
profile_single_file		= .true.
/
