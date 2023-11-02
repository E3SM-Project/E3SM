&ctl_nl
  nthreads          = 1
  partmethod        = 4
  topology          = "plane"
  geometry          = "plane"
  test_case         = "planar_rising_bubble"
  theta_hydrostatic_mode = .false.
  theta_advect_form = 1
  vert_remap_q_alg  = 10
  moisture          = 'moist'
  qsize             = 3
  statefreq         = 1000
  restartfreq       = -1
  runtype           = 0
  qsplit            = -1
  rsplit            = -1
  integration       = 'explicit'
  tstep_type        = 9
  ne_x              = 32
  ne_y              = 5
  nmax              = 300 ! 1800 is good for an actual simulation
  tstep             = 0.8
  hypervis_order    = 2
  hypervis_scaling  = 3.2
  nu                = 0.01
  nu_top            = 0.0
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hypervis subcycling
  hypervis_subcycle_tom = 1
  se_ftype=2
  limiter_option    = 9
  planar_slice=.true.
  lx = 20000.0
  ly = 1000.0
  sx = -10000.0
  sy = -500.0
  bubble_zcenter=2000.0
  bubble_ztop=20000.0
  bubble_xyradius=5000.0
  bubble_zradius =1000.0
  bubble_cosine=.true.
  bubble_moist=.true.
  bubble_T0=300.0
  bubble_dT=4.0
  bubble_moist_drh=1.0
  bubble_prec_type=1
  dt_remap_factor     = 5
  transport_alg       = 12
  dt_tracer_factor    = 10
  hypervis_subcycle_q = 10
  semi_lagrange_hv_q  = 1
  semi_lagrange_nearest_point_lev = 256
/
&vert_nl
  vanalytic         = 1                         ! set vcoords in initialization routine
  vtop              = 2.73919e-1                ! vertical coordinate at top of atm (z=10000m)
/
&analysis_nl
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 0,                        ! 1=days, 2=hours, 0=timesteps
  output_frequency  = 300,
  output_varnames1  ='T','Th','ps','geo','Q','u','v' !,'Q2','Q3'
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
