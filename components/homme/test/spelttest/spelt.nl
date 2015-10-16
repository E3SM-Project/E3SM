!===================================================================================!
! TEST INPUT for SPELT                                                              !
! Christoph Erath                                                                   !
! 	1 day = 1 * 24 * 3600 = 86400 sec		                                            !
! 	nmax  = ndays * 86400 / tstep 			                                            !
! 	12 days (1036800 s) 		                                                        !
!===================================================================================!
&ctl_nl
NThreads      = 1
partmethod    = 4
test_case     = "spelt_boomerang"
ne            = 20 !52 !52 !4 ! number of elements is ne*ne on each face, number must be >=2, o
                  !ne element per face is not allowed
test_cfldep   = .FALSE.
ndays	        = 12
ntrac         = 4
!nmax = 2
tstep         = 600 !300 !600   ! use factor 2 of 600 to reach the end time
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
output_end_time   = 100000
output_frequency  = 2
output_timeunits  = 1
output_varnames1  = 'geop','C', 'C2', 'C3'
interp_nlat       = 256
interp_nlon       = 512
output_type='netcdf'
/


! this will be used by homme only if namelist file is 'input.nl'
&prof_inparm
profile_single_file = .true.
profile_outpe_num = 0
/




!END of FILE- an extra line is important!