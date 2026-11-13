&ctl_nl
  nthreads          = -1
  partmethod        = 4
  topology          = "cube"
  test_case         = 'dcmip2012_test1_3e_conv'
  prescribed_wind   = 1
  qsize             = 4
  ndays             = 1
  statefreq         = 240
  restartfreq       = -1
  runtype           = 0
  ne                = 20
  integration       = 'explicit'
  tstep_type        = 1
  smooth            = 0
  nu                = 1.585e13 ! nu values are irrelevant
  nu_s              = 1.585e13
  nu_p              = 42 ! to satisfy kokkos exe
  se_ftype          = -1
  limiter_option    = 9
  hypervis_order    = 2
  hypervis_subcycle = 1
  moisture          = 'dry'
  theta_hydrostatic_mode = .true.
  dcmip16_prec_type = 1
  dcmip16_pbl_type  = -1
  transport_alg                     = 12
  semi_lagrange_cdr_alg             = 3
  semi_lagrange_cdr_check           = .false.
  semi_lagrange_hv_q                = 0
  semi_lagrange_nearest_point_lev   = 0
  semi_lagrange_halo                = 2
  dt_remap_factor                    = 0
  dt_tracer_factor                   = 4
  tstep                              = 200.0
  semi_lagrange_trajectory_nsubstep  = 2
  semi_lagrange_trajectory_nvelocity = 3
  semi_lagrange_diagnostics = 1
  hypervis_subcycle_q = 0
  limiter_option    = 9
  vert_remap_q_alg  = 10
/
&vert_nl
  vanalytic         = 1
  vtop              = 0.2549944
/
&analysis_nl
  output_dir        = "./movies/"
  output_timeunits  = 2,                    ! 1=days, 2=hours, 0=timesteps
  output_frequency  = 2,
  output_varnames1  = 'ps','Q','u','v'
  interp_type       = 0
  output_type       = 'netcdf'
  num_io_procs      = 16
  interp_nlon       = 180
  interp_nlat       = 91
  interp_gridtype   = 2
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
