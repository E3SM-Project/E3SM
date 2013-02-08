!=======================================================!
! 	1 day = 1 * 24 * 3600 = 86400 sec		!
! 	nmax  = ndays * 86400 / tstep 			!
! 	1 days at 10.0 stepsize: nmax= 8640 		!
!		  15.0 stepsize: nmax= 4320		!
!                 30.0 stepsize: nmax= 2880		!
!=======================================================!
!	1 hour= 1 * 3600 = 3600 sec			!
!	1 hr  at 10.0 stepsize:  remapfreq= 360  	!
!	1 hr  at 15.0 stepsize:  remapfreq= 240  	!
!	1 hr  at 30.0 stepsize:  remapfreq= 120  	!
!=======================================================!
&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
!test_case     = "heldsuarez"
test_case    = "jw_bcl"
ne            = 20
ndays         = 30
remapfreq     = 16
statefreq     = 8640
accumfreq     = 1080
accumstart    = 200
accumstop     = 400
tasknum       = 0
restartfreq   =-1
restartfile   = "./restart/"
runtype       = 0
tstep         = 10.0D0
integration   = "explicit"
smooth        = 0.005
!nu            = 7.0e5
nu            = 1.0e5
/
&solver_nl
precon_method = "block_jacobi"
maxits        = 100
tol           = 1.e-12
/
&filter_nl
transfer_type = "bv"
filter_type   = "taylor"
filter_freq   = 1
filter_mu     = 0.05D0
p_bv          = 12.0D0
s_bv          = 0.666666666666666666D0
wght_fm       = 0.10D0
kcut_fm       = 2
/
&vert_nl
vform         = "ccm"
vfile_mid     = "vcoord/camm-26.fbin.littleendian"
vfile_int     = "vcoord/cami-26.fbin.littleendian"
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
!                     'elem ','Time ','T    ','p3d  ', 	!
!		      'dgs  '				!
!  'dgs' includes 'ps','T850','zeta850' 		!
!							!
!  output_varnames1-5					!
!=======================================================!
output_frequency = 1
output_timeunits = 1
output_varnames1 ='ps','T','zeta'
output_type='netcdf'
io_stride=8
/
&dg_nl
/
