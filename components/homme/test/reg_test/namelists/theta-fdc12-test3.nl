&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2012_test3"       ! test identifier
  theta_hydrostatic_mode = .false.
  rsplit            = 1
  ne                = ${HOMME_TEST_NE}
  qsize             = 1                         ! num tracer fields
  nmax              = 180                       ! 7200s / 0.4s per step = 18000 steps
  statefreq         = 180                       ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 1.0                       ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 9                        ! 1 => default method
  nu                = 5.0e8                     ! reduced planet hyperviscosity hv/500^3
  nu_p              = 5.0e8
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  rearth            = 50969.0                   ! reduced planet radius rearth = a/500.0
  omega             = 0.0                       ! earth angular speed = 0.0
  limiter_option    = 9
  vert_remap_q_alg  = 1
/
&vert_nl
  vform             = "ccm"                     ! vertical coordinate type "ccm"=hybrid pressure/terrain
  vanalytic         = 1                         ! set vcoords in initialization routine
  vtop              = 3.2818e-2                 ! vertical coordinate at top of atm (z=30km)
/
&analysis_nl
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 0,                        ! 0=timesteps, 1=days, 2=hours, 3=seconds
  output_frequency  = 180,                      ! 720 seconds (10+1 outputs)
  output_varnames1  ='T','ps','u','v','omega','geo_i'   ! variables to write to file
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  !num_io_procs      = 16
  io_stride         = 8
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file   = .true.
/
