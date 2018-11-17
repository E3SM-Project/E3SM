&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "held_suarez0"
ne            = NE
ndays         = 400
statefreq     = SFREQ
tstep_type    = 5
qsize         = 1
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
output_start_time=600,600,100
output_end_time=30000,30000,110
output_varnames1='u','v','T','zeta','ps'
!output_varnames1='u','v','T','zeta','ps','Q','DIFFT'
!output_varnames2='DIFFT','DIFFU','DIFFV','CONVU','CONVV','FU','FV'
output_varnames3='geos','u','v','zeta','ps','div'
io_stride = 32
/

&prof_inparm
profile_outpe_num = 100
profile_single_file		= .true.
/








