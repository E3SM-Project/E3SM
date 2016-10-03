!
! namelist for dcmip2012 test 3-1 nonhydrostatic gravity waves
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2012_test3"         ! test identifier
  ne                = 7                         ! number of elements per cube face
  qsize             = 1                         ! num tracer fields
  ndays             = 0                         ! num simulation days: 0 = use nmax steps
  nmax              = 7200                      ! total number of steps: 7200 = 3600s / tstep=0.5s
  statefreq         = 20                        ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 = new run
  tstep             = 0.5                       ! largest timestep
  integration       = 'explicit'
  smooth            = 0                         ! timestep smooting (nonzero smoothing breaks this test)
  nu                = 8.0e6                     ! reduced earth hyperviz =1e15 / 500^3
  nu_s              = 8.0e6
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  rearth            = 50969.76                  ! scaled earth radius = a/125.0
  omega             = 0.0                       ! earth angular speed = 0.0
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
/
&analysis_nl
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 0,                        ! 1=days, 2=hours, 0=timesteps
  output_frequency  = 200,                      ! 100 sec / 0.5 sec per step
  output_varnames1  ='T','ps','u','v','omega'   ! variables to write to file
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  io_stride         = 8
  interp_nlat       = 65
  interp_nlon       = 128
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
