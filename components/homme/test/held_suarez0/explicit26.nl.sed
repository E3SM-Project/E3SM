&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "held_suarez0"
ne            = NE
ndays         = 400
statefreq     = SFREQ
!theta_hydrostatic_mode = .true.
!tstep_type    = 5
theta_hydrostatic_mode = .false.
tstep_type    = 7
qsize         = 1
theta_advect_form = 1
limiter_option = 9
rsplit        = 6
restartfreq   =  1
restartfile   = "restart/R0001"
restartdir    = "./restart/"
runtype       = RUNTYPE
tstep         = TSTEP
integration   = "explicit"
nu            = NU1
nu_s          = NU1
nu_p          = NU1
nu_q          = NU1
nu_div        = NU2
nu_top = 2.5e5
hypervis_order = 2
hypervis_subcycle = 2
se_ftype=0
/
&vert_nl
vform         = "ccm"
vfile_mid     = "../vcoord/camm-26.ascii"
vfile_int     = "../vcoord/cami-26.ascii"
/
&analysis_nl
infilenames=''
output_timeunits=1,0,2    ! 1=days, 2=hours, 3=seconds
output_frequency=1,0,0    ! 0 to disable
output_start_time=600,0,0
output_end_time=30000,999999999,0
output_varnames1='u','v','T','zeta','div','ps','geos','omega'
!output_varnames1='u','v','T','zeta','ps','Q','DIFFT'
! debug output
output_varnames2='u','v','T','zeta','div','ps','geo','dp3d','geos','Th'
! output3: hourly data for 20 days  
output_varnames3='geos','omega','zeta','ps','div'
io_stride = 32
/

&prof_inparm
profile_outpe_num = 100
profile_single_file		= .true.
/
