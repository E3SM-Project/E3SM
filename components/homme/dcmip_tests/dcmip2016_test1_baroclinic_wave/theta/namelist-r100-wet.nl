!
! theta: namelist for dcmip2016 test1: moist baroclinic wave
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2016_test1"         ! test identifier
  ne                = 30                        ! number of elements per cube face
  qsize             = 5                         ! num tracer fields
  ndays             = 30
  statefreq         = 10                        ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 80                        ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 7 
  rsplit            = 1
  qsplit            = 1
  nu                = 1e15                      ! default= 1e15*(ne30/ne30)**3.2 = 1e15
  nu_s              = 1e15
  nu_p              = 1e15
  nu_top            = 0                         ! default = 2.5e5
  limiter_option = 9
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  moisture          = 'wet'
  theta_hydrostatic_mode = .false.
  dcmip16_prec_type = 0                         ! 0=kessler physics
  dcmip16_pbl_type  = -1                        ! 0=reed-jablonowski pbl, -1 = none
/
&vert_nl
  vform             = "ccm"
  vfile_mid         = "../vcoord/camm-30.ascii"
  vfile_int         = "../vcoord/cami-30.ascii"
/
&analysis_nl
  output_prefix     = "r100-wet-"
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 2,                        ! 0=timesteps, 1=days, 2=hours, 3=seconds
  output_frequency  = 6                         ! every 6 hours
  output_varnames1  ='T','ps','pnh','geo','u','v','w','omega','Th','Q','Q2','Q3','Q4','Q5','rho','precl','zeta'   ! variables to write to file
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16
  interp_nlon       = 360
  interp_nlat       = 181
  interp_gridtype   = 1
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
