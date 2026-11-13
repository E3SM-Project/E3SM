&ctl_nl
  nthreads          = -1                        ! use OMP_NUM_THREADS
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2012_test1_3a_conv" ! test identifier
  prescribed_wind   = 1
  ne                = 5                         ! number of elements per cube face
  qsize             = 4                         ! num tracer fields
  ndays             = 12
  statefreq         = 10                        ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 28800                     ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 1 
  dt_remap_factor   = 1
  dt_tracer_factor  = 2
  smooth = 0
  nu                = 1.585e13
  nu_s              = 1.585e13
  se_ftype          = -1
  limiter_option    = 9
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  moisture          = 'dry'
  theta_hydrostatic_mode = .true.
  dcmip16_prec_type = 1                         ! 0=kessler physics
  dcmip16_pbl_type  = -1                        ! 0=reed-jablonowski pbl, -1 = none
  transport_alg     = 12
  semi_lagrange_cdr_alg   = 20
  semi_lagrange_cdr_check = .true.
  semi_lagrange_nearest_point_lev = 256
  limiter_option    = 9
  hypervis_subcycle_q = 2
  vert_remap_q_alg   = 10
  semi_lagrange_hv_q = 1
/
&vert_nl
  vanalytic         = 1
  vtop              = 0.2549944
/
&analysis_nl
!  output_prefix     = "PREFIX"
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 1,                        ! 1=days, 2=hours, 0=timesteps
  output_frequency  = 6,
  output_varnames1  ='u','Q','Q2','Q3','Q4'     ! variables to write to file
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16
  interp_nlon       = 180
  interp_nlat       = 91
  interp_gridtype   = 2
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
