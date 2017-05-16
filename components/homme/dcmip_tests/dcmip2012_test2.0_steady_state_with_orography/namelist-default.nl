!
! namelist for dcmip2012 test2-0: steady-state atmosphere with orography
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2012_test2_0"       ! test identifier
  ne                = 30                        ! number of elements per cube face
  qsize             = 1                         ! num tracer fields
  ndays             = 6                         ! num simulation days: 0 = use nmax steps
  statefreq         = 200                       ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 = new run
  tstep             = 20                        ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 1                         ! 1 => default method
  smooth            = 0                         ! timestep smooting
  nu                = 3e14                      ! hyperviscosity
  nu_s              = 3e14
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
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
  vtop              = 2.05e-1                   ! vertical coordinate at top of atm (z=12000m)
/
&analysis_nl
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 2,                        ! 1=days, 2=hours, 0=timesteps
  output_frequency  = 4,                        ! output every 4 hours
  output_varnames1  ='T','ps','u','v','omega'   ! variables to write to file
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16         
  interp_nlat       = 256
  interp_nlon       = 512
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
