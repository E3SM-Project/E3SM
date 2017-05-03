!
! namelist for dcmip2012 test2-1: nonhydro mountain waves without shear
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2016_test2"         ! test identifier
  theta_hydrostatic_mode = .true.
  rsplit            = 3
  ne                = 7                         ! number of elements per cube face
  qsize             = 3                         ! num tracer fields
  ndays             = 10
  statefreq         = 300                       ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 60                        ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 5                         ! 1 => default method
  nu                = 1e13
  nu_s              = 1e13
  nu_p              = 1e13
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
/
&vert_nl
  vform             = "ccm"                     ! vertical coordinate type "ccm"=hybrid pressure/terrain
  vanalytic         = 1                         ! set vcoords in initialization routine
/
&analysis_nl
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 2,                        ! 0=timesteps, 1=days, 2=hours, 3=seconds
  output_frequency  = 4,                        
  output_varnames1  ='T','Th','p','pnh','ps','u','v','w','geo','Q','Q2','Q3'   ! variables to write to file
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16
  interp_nlon       = 360
  interp_nlat       = 181
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
