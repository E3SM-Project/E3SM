! Prescribed-wind planar transport with spatially and temporally varying ps.

&ctl_nl
  nthreads          = 1
  partmethod        = 4
  topology          = "plane"
  geometry          = "plane"
  test_case         = "planar_transport_a"
  prescribed_wind   = 1
  theta_hydrostatic_mode = .true.
  vert_remap_q_alg  = 10
  moisture          = 'moist'
  qsize             = 3
  statefreq         = 1000
  restartfreq       = -1
  runtype           = 0
  integration       = 'explicit'
  tstep_type        = 5
  ne_x              = 128
  ne_y              = 5
  planar_slice      = .true.
  ndays             = 12
  hypervis_order    = 2
  hypervis_scaling  = 3.2
  nu                = 0.01
  nu_top            = 0.0
  nu_p              = 42 ! to satisfy kokkos exe
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hypervis subcycling
  hypervis_subcycle_tom = 1
  limiter_option    = 9
  tstep             = 1800
  transport_alg       = 12
  dt_remap_factor     = 0
  dt_tracer_factor    = 6
  hypervis_subcycle_q = 6
  semi_lagrange_hv_q  = 1
  semi_lagrange_trajectory_nsubstep  = 4        ! >0: Enhanced trajectory method
  semi_lagrange_trajectory_nvelocity = 3        ! >2: Extra velocity snapshots
  semi_lagrange_halo = 3
  semi_lagrange_nearest_point_lev = 0
/
&vert_nl
  vanalytic         = 1                         ! set vcoords in initialization routine
  vtop              = 2.73919e-1                ! vertical coordinate at top of atm (z=10000m)
/
&analysis_nl
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 1,                        ! 1=days, 2=hours, 0=timesteps
  output_frequency  = 6,
  output_varnames1  = 'Q','Q2','Q3','ps'
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       = 'netcdf'                  ! netcdf or pnetcdf
  num_io_procs      = 1
  interp_nlat       = 128
  interp_nlon       = 256
  interp_gridtype   = 2                         ! gauss grid
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
