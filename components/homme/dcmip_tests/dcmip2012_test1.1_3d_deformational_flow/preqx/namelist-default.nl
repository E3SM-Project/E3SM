!
! namelist for dcmip2012 test1-1: 3d deformational flow
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2012_test1_1"       ! test identifier
  ne                = 30                        ! number of elements per cube face
  qsize             = 4                         ! num tracer fields
  ndays             = 12                         ! num simulation days: 0 = use nmax steps
  statefreq         = 200                       ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 = new run
  tstep             = 33                        ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 1                         ! 1 => default method
  smooth            = 0                         ! timestep smooting
  nu                = 1.585e13                  ! hyperviscosity minimum L2 error at 10^13.2
  nu_s              = 1.585e13
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  prescribed_wind   = 1
/
&filter_nl/
&solver_nl
  precon_method     = "identity"
  maxits            = 50
  tol               = 1.e-7
/
&vert_nl
  vform             = "ccm"                     ! vertical coordinate type "ccm"=hybrid pressure/terrain
  vanalytic         = 1                         ! set vcoords in initialization routine
  vtop              = 0.2549944                 ! vertical coordinate at top of atm 254.9 hPa (12km)
/
&analysis_nl
  output_dir        = "./movies/"              ! destination dir for netcdf file
  output_timeunits  = 1,                        ! 1=days, 2=hours, 0=timesteps
  output_frequency  = 1,                        ! output every 4 hours
  output_varnames1  ='u','geo','Q','Q2','Q3','Q4' ! variables to write to file
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16         
  interp_nlat       = 256
  interp_nlon       = 512
  interp_gridtype   = 2                         ! gauss grid
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
