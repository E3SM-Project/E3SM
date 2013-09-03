!===================================================================================!
! TEST INPUT for CWFV                                                              !
! Christoph Erath                                                                   !
! 	1 day = 1 * 24 * 3600 = 86400 sec		                                            !
! 	nmax  = ndays * 86400 / tstep 			                                            !
! 	12 days (1036800 s) 		                                                        !
!===================================================================================!
&ctl_nl
NThreads      = 1
partmethod    = 4
test_case     = "cwfv_boomerang"
ne            = 120 !4 ! number of elements is ne*ne on each face, number must be >=2, o
                  !ne element per face is not allowed
test_cfldep   = .TRUE.
ndays	        = 12
ntrac         = 2
tstep         = 75   ! use factor 2 of 600 to reach the end time
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
!  output_dir ="./results/" is default				!
!							!
!  allowed variables: 'ps   ','geop ', 'c    ', u    ','v    ',	!
!                     'latp ','lonp ','latv ','lonv ',	!
!                     'elem ','Time ' 			!
!							!
!  output_varnames1-5					!
!=======================================================!    
output_start_time = 0
output_end_time   = 288
output_frequency  = 6
output_timeunits  = 2
output_varnames1  = 'geop'
! interp_nlat       = 128
! interp_nlon       = 256
output_dir ="/glade/home/erath/hommelynx/test/cwfvtest/results/"
output_type='netcdf'
/


! this will be used by homme only if namelist file is 'input.nl'
&prof_inparm
profile_single_file = .true.
profile_outpe_num = 5000
/




!END of FILE- an extra line is important!