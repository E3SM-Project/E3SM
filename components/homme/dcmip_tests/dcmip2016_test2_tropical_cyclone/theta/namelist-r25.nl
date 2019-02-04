!
! theta: namelist for dcmip2016 test2: tropical cyclone
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2016_test2"         ! test identifier
  ne                = 120                       ! number of elements per cube face (1/4dg)
  qsize             = 4                         ! num tracer fields
  ndays             = 10
  statefreq         = 10                        ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 25                        ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 5
  rsplit            = 3
  qsplit            = 1
  nu                = 1.2e13                    ! default= 1e15*(ne30/ne60)**3.2 = 1.1e14
  nu_s              = 1.2e13
  nu_p              = 1.2e13
  nu_top            = 2.5e5                     ! default = 2.5e5
  limiter_option = 9
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  moisture          = 'wet'
  theta_hydrostatic_mode = .true.
  dcmip16_prec_type = 0                         ! 0=kessler prec 1= reed prec
  dcmip16_pbl_type  = 0                         ! 0=reed pbl     1= bryan pbl
/
&vert_nl
  vform         = "ccm"
  vfile_mid     = "../vcoord/camm-30.ascii"
  vfile_int     = "../vcoord/cami-30.ascii"
/
&analysis_nl
  output_prefix     = "r25-prec0-pbl0-"             ! which prec & pbl type?
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 2,                        ! 0=timesteps, 1=days, 2=hours, 3=seconds
  output_frequency  = 6                         ! every 6 hours
  output_varnames1  ='T','ps','pnh','geo','u','v','w','omega','Th','Q','Q2','Q3','precl'   ! variables to write to file
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16
  interp_nlon       = 1440
  interp_nlat       = 721
  interp_gridtype   = 1
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
