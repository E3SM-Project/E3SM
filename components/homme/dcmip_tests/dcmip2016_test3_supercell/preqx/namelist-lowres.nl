!
! preqx: namelist for dcmip2016 test 3: supercell storm (small planet X=120)
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2016_test3"         ! test identifier
  ne                = 7                         ! number of elements per cube face
  qsize             = 4                         ! num tracer fields: qv,qc,qr
  nmax              = 36000                     ! 7200s(120min)/tstep
  statefreq         = 10                        ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 0.2                       ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 5
  rsplit            = 1
  qsplit            = 1
  nu                = 4.57e10                   ! 1e15/(120)^3 *(ne30/ne7)
  nu_s              = 4.57e10
  nu_p              = 4.57e10
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  rearth            = 53133                     ! 6.376E6  / 120
  omega             = 0
  se_ftype          = 0
  moisture          = 'dry'
/
&vert_nl
  vform             = "ccm"                     ! vertical coordinate type "ccm"=hybrid pressure/terrain
  vanalytic         = 1                         ! set vcoords in initialization routine
  vtop              = 5e-2                      ! vertical coordinate at top of atm (z=20km)
/
&analysis_nl
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 0   !3,                   ! 0=timesteps, 1=days, 2=hours, 3=seconds
  output_frequency  = 100 !900,                 ! 900 seconds (15 minutes)
  output_varnames1  ='T','p','ps','pnh','geo','u','v','omega','Th','Q','Q2','Q3','Q4'   ! variables to write to file
  interp_nlon       = 360
  interp_nlat       = 181
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16         
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
