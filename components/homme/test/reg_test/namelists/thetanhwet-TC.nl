!
! theta: namelist for dcmip2016 test2: tropical cyclone
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2016_test2"         ! test identifier
  ne                = 8                         ! number of elements per cube face
  qsize             = 4                         ! num tracer fields
  ndays             = 1
  statefreq         = 24                        ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 450                      ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 7
  rsplit            = 2
  qsplit            = 2
  nu                = 3e16                      ! default= 1e15*(ne30/ne8)**3.2 = 6.9e16
  nu_s              = 3e16
  nu_p              = 3e16  
  nu_top            = 2.5e5                     ! default = 2.5e5
  limiter_option    = 8
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  moisture          = 'wet'
  theta_hydrostatic_mode = .false.
  dcmip16_prec_type = 0                         ! 0=kessler,     1= reed-jablonowski
  dcmip16_pbl_type  = 0                         ! 0=basic pbl,   1= bryan pbl
/
&vert_nl
  vform         = "ccm"
  vfile_mid     = "./vcoord/camm-30.ascii"
  vfile_int     = "./vcoord/cami-30.ascii"
/
&analysis_nl
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 2,                        ! 0=timesteps, 1=days, 2=hours, 3=seconds
  output_frequency  = 24                        ! every 3 hours
  output_varnames1  ='T','w','ps','Q','Q2'      ! variables to write to file
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16
  interp_nlon       = 180
  interp_nlat       = 91
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
