!
! namelist for dcmip2012 test2-1: nonhydro mountain waves without shear
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2012_test2_1"       ! test identifier
  theta_hydrostatic_mode = .false.
  rsplit = 0
  ne                = 20                        ! number of elements per cube face
  qsize             = 1                         ! num tracer fields
  ndays             = 0                         ! num simulation days: 0 => use nmax steps
  nmax              = 14400                     ! 7200s / 0.1s per step = 72000 steps
  statefreq         = 30                      ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 0.5                       ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 5                         ! 1 => default method
  nu                = 3.2e7                       ! reduced planet hyperviscosity hv/500^3
  nu_s              = 3.2e7
  nu_p              = 3.2e7
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  rearth            = 12752.0                   ! reduced planet radius rearth = a/500.0
  omega             = 0.0                       ! earth angular speed = 0.0
  dcmip2_x_ueq      = 20.0                      ! windspeed at equator  (m/s)
  dcmip2_x_h0       = 250.0                     ! mountain height       (m)
  dcmip2_x_d        = 5000.0                    ! mountain half-width   (m)
  dcmip2_x_xi       = 4000.0                    ! mountain wavelength   (m)
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
  vtop              = 3.2818e-2                 ! vertical coordinate at top of atm (z=30km)
/
&analysis_nl
  output_dir        = "./movies/"              ! destination dir for netcdf file
  output_timeunits  = 0,                        ! 1=days, 2=hours, 0=timesteps
  output_frequency  = 1000,                     ! 100s /0.1s = 1000 steps between outputs
  output_varnames1  ='T','ps','u','v','omega'   ! variables to write to file
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16         
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
