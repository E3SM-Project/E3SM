&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "held_suarez0"
ne            = 4
mesh_file = 'none'    
qsize = 0
nmax          = 1
statefreq     = 1
restartfreq   = 43200
restartfile   = "./R0001"
runtype = 0
tstep         = 0.0001
tstep_type   = 1
qsplit=1
integration   = "explicit"
smooth        = 0
nu=0
nu_s          = 0
nu_p          = 0
hypervis_order = 2
hypervis_subcycle=1
hypervis_subcycle_q=1
u_perturb      = 1
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
 output_prefix     = "template1-"
 output_timeunits=1
 output_frequency = 1
 infilenames='h0-tavg.nc'
 output_varnames1='area','corners','hypervis','cv_lat','cv_lon','phys_lat','phys_lon','phys_cv_lat','phys_cv_lon'
! removing 'phys_area' since this is not BFB accross different runs 
 output_type='netcdf'
! num_io_procs = 1
 io_stride = 8
/

&prof_inparm
profile_outpe_num = 8
profile_single_file		= .true.
/
