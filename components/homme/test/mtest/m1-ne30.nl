! M1 tests as in KSP2015, runs with 60 levels.
!_______________________________________________________________________
&ctl_nl
  theta_hydrostatic_mode = .false.
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "mtest1"       ! test identifier
  ne                = 30                        ! number of elements per cube face
  qsize             = 0                         ! num tracer fields
  ndays             = 0                         ! num simulation days: 0 => use nmax steps
  nmax              = 72000                     ! 7200s / 0.1s per step = 72000 steps
  statefreq         = 1000                      ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 0.1                       ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 5                         ! 1 => default method
  nu                = 2.2e8                       ! reduced planet hyperviscosity hv/500^3
  nu_s              = 2.2e8
  nu_p              = 2.2e8
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  rearth            = 38248.0                   ! reduced planet radius rearth = a/500.0
  omega             = 0.0                       ! earth angular speed = 0.0
/
&vert_nl
  vform             = "ccm"                     ! vertical coordinate type "ccm"=hybrid pressure/terrain
  vanalytic         = 1                         ! set vcoords in initialization routine
  vtop              = 3.2818e-2                 ! vertical coordinate at top of atm (z=30km)
/
&analysis_nl
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 0,                        ! 1=days, 2=hours, 0=timesteps
  output_frequency  = 72000,                     ! 100s /0.1s = 1000 steps between outputs
  output_varnames1  ='T','ps','u','v','omega'   ! variables to write to file
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16         
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
