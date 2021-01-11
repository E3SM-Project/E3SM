!
! namelist for planar nonhydrostatic gravity waves
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "plane"                   ! mesh type: planar
  geometry          = "plane"                   ! mesh type: planar
  test_case         = "planar_nonhydro_gravity_wave"         		! test identifier
  theta_hydrostatic_mode = .true.
  transport_alg     = 0            
  theta_advect_form = 1
  vert_remap_q_alg  = 10
  moisture          = 'dry'
  ne_x              = 33                        ! number of elements in x-dir
  ne_y              = 33                        ! number of elements in y-dir
  qsize             = 0                         ! num tracer fields
  nmax              = 200                       ! total number of steps: 600s / tstep
  statefreq         = 200                       ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 = new run
  qsplit            = -1                        ! timesteps set via se_tstep, dt_remap_factor, dt_tracer_factor
  rsplit            = -1
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 5                         ! 1 => default method
  tstep             = 3.0                       ! dynamics timestep
  dt_remap_factor   = 1				! remap every 1 time steps
  dt_tracer_factor  = 1                         ! tracers run at dynamics time step
  hypervis_scaling  = 3.0
  nu                = 0.216784
  nu_top            = 0.0             
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  hypervis_subcycle_tom = 1
  hypervis_subcycle_q = 1
/
&vert_nl
  vform             = "ccm"                     ! vertical coordinate type "ccm"=hybrid pressure/terrain
  vanalytic         = 1                         ! set vcoords in initialization routine
  vtop              = 2.73919e-1                ! vertical coordinate at top of atm (z=10000m)
/
&analysis_nl
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 0,                        ! 1=days, 2=hours, 0=timesteps
  output_frequency  = 200,                      ! steps
  output_varnames1  ='T','Th','ps','u','v','omega'   ! variables to write to file
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 1         
  interp_nlat       = 128
  interp_nlon       = 256
  interp_gridtype   = 2                         ! gauss grid
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/


