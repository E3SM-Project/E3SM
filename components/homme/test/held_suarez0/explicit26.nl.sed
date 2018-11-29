&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "held_suarez0"
ne            = NE
ndays         = 400
statefreq     = SFREQ
!cubed_sphere_map = 0
theta_hydrostatic_mode = .true.
tstep_type    = 5
!theta_hydrostatic_mode = .false.
!tstep_type    = 7
!vert_remap_q_alg=1
qsize         = 1
limiter_option = 9
rsplit        = 3
restartfreq   =  1
restartfile   = "restart/R0001"
restartdir    = "./restart/"
!restartdir    = "bglockless:./restart/"
runtype       = RUNTYPE
tstep         = TSTEP
integration   = "explicit"
nu            = NU1
nu_s          = NU2
nu_p          = NUP
!nu_top = 2.5e5
hypervis_order = 2
hypervis_subcycle = 1
/
&vert_nl
vform         = "ccm"
vfile_mid     = "../vcoord/camm-26.ascii"
vfile_int     = "../vcoord/cami-26.ascii"
/
&analysis_nl
infilenames=''
output_timeunits=1,0,2    ! 1=days, 2=hours, 3=seconds
output_frequency=1,1,1
output_start_time=600,2090,4800
output_end_time=30000,2110,5280
output_varnames1='u','v','T','zeta','div','ps'
!output_varnames1='u','v','T','zeta','ps','Q','DIFFT'
! debug output
!output_varnames2='u','v','T','zeta','div','ps','geo','dp3d','geos','Th'
! output3: hourly data for 20 days  
output_varnames3='geos','omega','zeta','ps','div'
io_stride = 32
/

&prof_inparm
profile_outpe_num = 100
profile_single_file		= .true.
/








