!
! namelist for development testing of theta model
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "jw_baroclinic"           ! test identifier
  ne                = ${HOMME_TEST_NE}          ! number of elements per cube face
  qsize             = ${HOMME_TEST_QSIZE}       ! num tracer fields
  ndays             = ${HOMME_TEST_NDAYS}
  statefreq         = 24                        ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 450                      ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 5
  rsplit            = ${HOMME_TEST_RSPLIT}
  qsplit            = ${HOMME_TEST_QSPLIT}
  nu                = 3e16                      ! default= 1e15*(ne30/ne8)**3.2 = 6.9e16
  nu_s              = 3e16
  nu_p              = 3e16  
  nu_top            = 0e5                     ! default = 2.5e5
  limiter_option    = 8
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  moisture          = '${HOMME_TEST_MOISTURE}'
  dcmip16_prec_type = 0                         ! 0=kessler,     1= reed-jablonowski
  dcmip16_pbl_type  = 0                         ! 0=basic pbl,   1= bryan pbl
  vert_remap_q_alg  = 1
  theta_hydrostatic_mode = .${HOMME_TEST_HYDROSTATIC_MODE}.
/
&vert_nl
  vfile_mid     = "./vcoord/${HOMME_TEST_VCOORD_MID_FILE}"
  vfile_int     = "./vcoord/${HOMME_TEST_VCOORD_INT_FILE}"
/
&analysis_nl
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 2,                        ! 0=timesteps, 1=days, 2=hours, 3=seconds
  output_frequency  = 24                        ! every 3 hours
  output_varnames1  ='T','w','ps','Q','Q4'      ! variables to write to file
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
