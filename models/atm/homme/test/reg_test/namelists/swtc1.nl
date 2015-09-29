!=======================================================!
! 	1 day = 1 * 24 * 3600 = 86400 sec		!
! 	nmax  = ndays = 12
! 	12 days at 120.0 stepsize: nmax= 8640 		!
!=======================================================!
&ctl_nl
NThreads          = 1
partmethod        = 4
topology          = "cube"
test_case         = 'swtc1'   
ne                = 10
ndays             = 12
statefreq         = 180
tasknum           = 0
restartfreq       = -1
restartfile       = "./restart/R000000050"
runtype           = 0
tstep             = 480
integration       = "runge_kutta"
rk_stage_user     = 3
LFTfreq           = 0
smooth            = 0
limiter_option    = 4
nu                = 0
nu_s              = 0
hypervis_order    = 2
hypervis_subcycle = 1
/
&solver_nl
precon_method     = "block_jacobi"
maxits            = 100
tol               = 1.e-12
/
&filter_nl
transfer_type     = "bv"
filter_type       = "taylor"
filter_freq       = 0
filter_mu         = 0.005
p_bv              = 12.0D0
s_bv              = .80
wght_fm           = 0.10D0
kcut_fm           = 2
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
output_start_time = 0
output_end_time   = 12
output_frequency  = 12
output_timeunits  = 1
output_varnames1  = 'u', 'v', 'ps', 'geop'
interp_nlat       = 128
interp_nlon       = 256
output_type       = 'netcdf'
/

