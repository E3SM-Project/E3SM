!
! namelist for planar nonhydrostatic gravity waves
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "plane"                   ! mesh type: planar
  geometry          = "plane"                   ! mesh type: planar
  test_case         = "planar_rising_bubble"         		! test identifier
  theta_hydrostatic_mode = .false.
!  theta_hydrostatic_mode = .true.
  transport_alg     = 0            
  theta_advect_form = 1
  vert_remap_q_alg  = 10
  moisture          = 'moist'
  ne_x              = 40 !68                       ! number of elements in x-dir
  ne_y              = 4                       ! number of elements in y-dir
  qsize             = 3                         ! num tracer fields
  nmax              = 14000                     ! total number of steps: 600s / tstep
  statefreq         = 100                       ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 = new run
  qsplit            = -1                        ! timesteps set via se_tstep, dt_remap_factor, dt_tracer_factor
  rsplit            = -1
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 9                         ! 1 => default method
!  tstep_type        = 5                         ! 1 => default method
  tstep             = 0.02 !! 20x20x60lev 0.1                       ! dynamics timestep
  dt_remap_factor   = 1				! remap every 1 time steps
  dt_tracer_factor  = 1                        ! tracers run at dynamics time step
  hypervis_order    = 2 
  hypervis_scaling  = 3.0
!  hypervis_scaling  = 0.0
  nu                = 0.216784  !dry ran with 0.01
!!!0.02 !!!!default EAM is 3.4e-8 for earth circ=40000
  nu_top            = 0.0             
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  hypervis_subcycle_tom = 1
  hypervis_subcycle_q = 1
  se_ftype=2
  limiter_option    = 9
  planar_slice=.true.
  lx = 20000.0
  ly = 1000.0
  sx = -10000.0
  sy = -500.0
  bubble_zcenter=2000.0
  bubble_ztop=20000.0
  bubble_xyradius=3000.0
  bubble_zradius =1000.0
  bubble_cosine=.true.
!  bubble_cosine=.false.
!  bubble_moist=.false.
  bubble_moist=.true.
  bubble_T0=300.0
  bubble_dT=2.0
  bubble_moist_dq=0.2
/
&vert_nl
  vanalytic         = 1                         ! set vcoords in initialization routine
  vtop              = 2.73919e-1                ! vertical coordinate at top of atm (z=10000m)
/
&analysis_nl
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 0,                        ! 1=days, 2=hours, 0=timesteps
  output_frequency  = 200,                        ! steps
  output_varnames1  ='T','Th','ps','geo','Q1','Q2','Q3'   ! variables to write to file
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


