!=======================================================!
! 	1 day = 1 * 24 * 3600 = 86400 sec		!
! 	nmax  = ndays * 86400 / tstep = 30
! 	12 days at 30.0 stepsize: nmax= 34560		!
!=======================================================!
&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "swtc2"
ndays         = 1
statefreq = 80
tasknum       = 0
restartfreq   = -1
restartfile   = "./restart/R000000050"
runtype       = 0
tstep = 90
integration   = "explicit"
rk_stage_user = 0
LFTfreq       = 1
smooth        = 0.0
nu= 1e-9
nu_s= 1e-9
hypervis_order = 2
hypervis_subcycle = 2
hypervis_power = 0
mesh_file = 'mountain_10_x2.g'
hypervis_scaling = 3.2
/
&filter_nl
filter_freq   = 0
/
&analysis_nl
!=======================================================!
!  currently up to 5 streams are allowed		!
!  output_stream_count=1				!
!							!
!  timunits: 0= steps, 1=days, 2=hours			!
!  output_timeunits=1,2 				!
!  output_start_time=0,1176				!			
!  output_end_time=-1,-1				!
!  output_frequency=1,1 				!
!  output_dir ="./movies/"				!
!							!
!  allowed variables: 'ps   ','geop ','u    ','v    ',	!
!                     'latp ','lonp ','latv ','lonv ',	!
!                     'elem ','Time ' 			!
!							!
!  output_varnames1-5					!
!=======================================================!
output_start_time= 0
output_end_time  = -1
output_timeunits = 1
output_frequency = 1
interp_nlat=181
interp_nlon=360
output_varnames1 = 'geop','zeta'
io_stride = 8
output_type = 'netcdf'
/







