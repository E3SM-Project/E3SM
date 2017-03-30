!
! namelist for dcmip2012 test2-0: steady-state atmosphere with orography
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2012_test2_0"       ! test identifier
  theta_hydrostatic_mode = .true.
  ne                = 30                        ! number of elements per cube face
  qsize             = 0                         ! num tracer fields
  ndays             = 6                         ! num simulation days: 0 = use nmax steps
  statefreq         = 288                       ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 = new run
  tstep             = 300                       ! largest timestep in seconds
  rsplit            = 3
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 5                         ! 1 => default method
  nu                = 1e15                      ! hyperviscosity
  nu_s              = 1e15
  nu_p              = 1e15
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  omega             = 0.0                       ! earth angular speed = 0.0
/
&vert_nl
  vform             = "ccm"                     ! vertical coordinate type "ccm"=hybrid pressure/terrain
  vanalytic         = 1                         ! set vcoords in initialization routine
  vtop              = 2.05e-1                   ! vertical coordinate at top of atm (z=12000m)
/
&analysis_nl
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 2,                        ! 1=days, 2=hours, 0=timesteps
  output_frequency  = 12,                       ! output every 12 hours
  output_varnames1  ='T','ps','u','v','omega','geo'   ! variables to write to file
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16         
  interp_nlat       = 91
  interp_nlon       = 360
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
