&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "jw_baroclinic"
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
smooth        = 0
nu=0
nu_s          = 0
nu_p          = 0
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
 output_varnames1='area','corners','cv_lat','cv_lon'
! output_varnames1='area','corners','hypervis','cv_lat','cv_lon'
! output_type='netcdf'  
! output_type='pnetcdf'  
! output_type='pnetcdf64'  ! needed for ne1024
 output_type='netcdf4p'  ! needed for ne1024
! num_io_procs = 1
 io_stride = 8
/

&prof_inparm
profile_outpe_num = 8
profile_single_file		= .true.
/
