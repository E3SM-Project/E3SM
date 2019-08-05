!
! preqx: namelist for dcmip2016 test1: moist baroclininc wave
!_______________________________________________________________________
&ctl_nl
  nthreads          = -1                        ! use OMP_NUM_THREADS
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2016_test1"         ! test identifier
  ne                = 8                         ! number of elements per cube face
  qsize             = 6                         ! num tracer fields
  nmax              = 3
  statefreq         = 1                         ! number of steps between screen dumps
  restartfreq       = 0                         ! don't write restart files if < 0
  restartfile   = "./restart/R000000648"
  se_ftype = -1
  runtype           = 1                         ! 0 => new run
  tstep             = 1200                       ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 7 
  rsplit            = 0
  qsplit            = 1
  nu                = 0e16                      ! default= 1e15*(ne30/ne8)**3.2 = 6.9e16
  nu_s              = 0e16
  nu_p              = 0e16  
  nu_top            = 0                         ! default = 2.5e5
  limiter_option    = 9
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  moisture          = 'wet'
  theta_hydrostatic_mode = .false.
  dcmip16_prec_type = 1                         ! 0=kessler physics
  dcmip16_pbl_type  = -1                        ! 0=reed-jablonowski pbl, -1 = none
/
&vert_nl
  vform             = "ccm"
  vfile_mid         = "../vcoord/camm-30.ascii"
  vfile_int         = "../vcoord/cami-30.ascii"
/
&analysis_nl
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 2,                        ! 0=timesteps, 1=days, 2=hours, 3=seconds
  output_frequency  = 24                        ! every N hours
  output_varnames1  ='T','ps','Q','precl','zeta'   ! variables to write to file
  interp_type       = 1                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16
  interp_gridtype   = 1
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
