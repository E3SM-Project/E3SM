!
! theta: namelist for dcmip2016 test 3: supercell storm (small planet X=120)
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2016_test3"         ! test identifier
  ne                = 60                        ! number of elements per cube face
  qsize             = 4                         ! num tracer fields
  nmax              = 28800                     ! 7200s(120min)/tstep
  statefreq         = 100                       ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 0.25                      ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 5
  rsplit            = 0
  qsplit            = 10
  nu                = 6.3e7                     ! default= 1e15/(120)^3 *(ne30/ne60)**3.2
  nu_s              = 6.3e7
  nu_p              = 0
  nu_q              = 0
  nu_top            = 0                         ! 2.5e5/(120)^(1)
  limiter_option    = 4
  dcmip16_mu        = 500.0d0                   ! additional uniform viscosity
  dcmip16_mu_s      = 1500.0d0
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  rearth            = 53133                     ! 6.376E6  / 120
  omega             = 0
  se_ftype          = 0
  moisture          = 'wet'
  theta_hydrostatic_mode = .false.
/
&vert_nl
  vform             = "ccm"                     ! vertical coordinate type "ccm"=hybrid pressure/terrain
  vanalytic         = 1                         ! set vcoords in initialization routine
/
&analysis_nl
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 3                         ! 0=timesteps, 1=days, 2=hours, 3=seconds
  output_frequency  = 60                       ! 300 seconds
  output_varnames1  ='geo','w','Q','Q2','Q3','precl'   ! variables to write to file
  interp_nlon       = 720
  interp_nlat       = 362
  interp_gridtype   = 2                         ! gauss grid
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  interp_lon0       = -180.0                    ! shift lon range to [-180,+180)
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16         
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
