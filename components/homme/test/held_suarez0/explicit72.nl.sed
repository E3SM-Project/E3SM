&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "held_suarez0"
ne            = NE
ndays         = 400
statefreq     = SFREQ
theta_hydrostatic_mode = .true.
tstep_type    = 5
!theta_hydrostatic_mode = .false.
!tstep_type    = 9
qsize         = 1
theta_advect_form = 1
pgrad_correction=1
hv_ref_profiles=2
hv_theta_correction=1
limiter_option = 9
dt_remap_factor = 1
dt_tracer_factor = 1
vert_remap_q_alg = 10
restartfreq   =  1
restartfile   = "restart/R0001"
restartdir    = "./restart/"
runtype       = RUNTYPE
tstep         = TSTEP
integration   = "explicit"
nu            = NU1
nu_top = 2.5e5  ! default 2.5e5    HSV1 1.5ok.  2.0 bad
                ! timesplit version.  5e5 works.  10e5 crashes.  
hypervis_scaling = 3  ! 0 for constant coeff HV
hypervis_order = 2
hypervis_subcycle = 2
hypervis_subcycle_tom = 1
/
&vert_nl
vfile_mid     = "../vcoord/acme-72m.ascii"
vfile_int     = "../vcoord/acme-72i.ascii"
/
&analysis_nl
infilenames=''
output_timeunits=1,2,2    ! 0=tsteps, 1=days, 2=hours, 3=seconds
output_frequency=1,0,0    ! 0 to disable.
output_start_time=200,14640,0
output_end_time=500,14641,9999999
output_varnames1='ps','omega'
!output_varnames1='u','v','T','ps','geos','omega'
! debug output
output_varnames2='u','v','T','zeta','div','ps','geo','dp3d','geos','Th','omega'
! output3: hourly data for 20 days  
output_varnames3='geos','omega','zeta','ps','div'
io_stride = 32
/

&prof_inparm
profile_outpe_num = 100
profile_single_file		= .true.
/
