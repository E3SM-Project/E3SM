!
! preqx: namelist for dcmip2016 test1: moist baroclininc wave
!_______________________________________________________________________
&ctl_nl
  nthreads          = -1                        ! use OMP_NUM_THREADS
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2016_test1_pg2"     ! test identifier
  ne                = 4                         ! number of elements per cube face
  qsize             = 6                         ! num tracer fields
  ndays             = 6
  statefreq         = 50                        ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 900                       ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 5 
  dt_tracer_factor  = 6
  dt_remap_factor   = 2
  nu                = 6e17                      ! default= 1e15*(ne30/ne8)**3.2 = 6.9e16
  nu_s              = 6e17
  nu_p              = 6e17  
  nu_top            = 0                         ! default = 2.5e5
  limiter_option    = 9
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  moisture          = 'wet'
  theta_hydrostatic_mode = .true.
  dcmip16_prec_type = 1                         ! 0=kessler physics
  dcmip16_pbl_type  = -1                        ! 0=reed-jablonowski pbl, -1 = none
  transport_alg     = 12
  semi_lagrange_cdr_alg   = 2
  semi_lagrange_cdr_check = .true.
  semi_lagrange_nearest_point_lev = 100
  vert_remap_q_alg   = 10
  nu_q = 0
/
&vert_nl
  vfile_mid         = "vcoord/camm-30.ascii"
  vfile_int         = "vcoord/cami-30.ascii"
/
&analysis_nl
!  output_prefix     = "dc1SUFFIX-"
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 1,                        ! 0=timesteps, 1=days, 2=hours, 3=seconds
  output_frequency  = 3
  output_varnames1  ='T','ps','Q','Q4','Q5','precl'  ! variables to write to file
  interp_type       = 1                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16
  interp_nlon       = 180
  interp_nlat       = 91
  interp_gridtype   = 1
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
