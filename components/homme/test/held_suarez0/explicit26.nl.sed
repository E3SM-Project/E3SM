&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
theta_hydrostatic_mode = .true.
test_case     = "held_suarez0"
ne            = NE
ndays         = 400
statefreq     = SFREQ
tstep_type    = 5
qsize         = 1
limiter_option = 9
rsplit        = 1
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
output_timeunits=1,1,2    ! 1=days, 2=hours, 3=seconds
output_frequency=1,1,1
output_start_time=600,600,4800
output_end_time=30000,30000,5280
output_varnames1='u','v','T','zeta','div','ps'
!output_varnames1='u','v','T','zeta','ps','Q','DIFFT'
!output_varnames2='DIFFT','DIFFU','DIFFV','CONVU','CONVV','FU','FV'
! output3: hourly data for 20 days  
output_varnames3='geos','omega','zeta','ps','div'
io_stride = 32
/

&prof_inparm
profile_outpe_num = 100
profile_single_file		= .true.
/








