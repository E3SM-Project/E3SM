#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#if !defined(CAM) && !defined(SCREAM)
#include "homme_git_sha.h"
#endif

module namelist_mod

  use kinds,      only: real_kind, iulog
  use params_mod, only: recursive, sfcurve, SPHERE_COORDS, Z2_NO_TASK_MAPPING
  use cube_mod,   only: rotate_grid
#ifdef CAM
  use dyn_grid,   only: fv_nphys
#endif
  use physical_constants, only: rearth, rrearth, omega
#if (defined MODEL_THETA_L && defined ARKODE)
  use arkode_mod, only: rel_tol, abs_tol, calc_nonlinear_stats, use_column_solver
#endif
use physical_constants, only : rearth, rrearth, DD_PI

use physical_constants, only : scale_factor, scale_factor_inv, domain_size, laplacian_rigid_factor
use physical_constants, only : Sx, Sy, Lx, Ly, dx, dy, dx_ref, dy_ref

  use control_mod, only : &
    MAX_STRING_LEN,&
    MAX_FILE_LEN,  &
    partmethod,    &       ! Mesh partitioning method (METIS)
    coord_transform_method,    &       !how to represent the coordinates.
    z2_map_method,    &       !zoltan2 how to perform mapping (network-topology aware)
    topology,      &       ! Mesh topology
    geometry,      &       ! Mesh geometry
    test_case,     &       ! test case
    planar_slice,     &
    numnodes,      &
    sub_case,      &
    statefreq,     &       ! number of steps per printstate call
    restartfreq,   &
    restartfile,   &       ! name of the restart file for INPUT
    restartdir,    &       ! name of the restart directory for OUTPUT
    runtype,       &
    integration,   &       ! integration method
    theta_hydrostatic_mode,       &   
    transport_alg , &      ! SE Eulerian, classical SL, cell-integrated SL
    semi_lagrange_cdr_alg, &     ! see control_mod for semi_lagrange_* descriptions
    semi_lagrange_cdr_check, &
    semi_lagrange_hv_q, &
    semi_lagrange_nearest_point_lev, &
    tstep_type,    &
    cubed_sphere_map, &
    qsplit,        &
    rsplit,        &
    dt_tracer_factor, &
    dt_remap_factor, &
    rk_stage_user, &
    LFTfreq,       &
    prescribed_wind, &
    ftype,         &
    limiter_option,&
    nu,            &
    nu_s,          &
    nu_q,          &
    nu_div,        &
    nu_p,          &
    nu_top,        &
    tom_sponge_start, &
    dcmip16_mu,     &
    dcmip16_mu_s,   &
    dcmip16_mu_q,   &
    dcmip16_prec_type, &
    dcmip16_pbl_type,&
    interp_lon0,    &
    hypervis_scaling,   &  ! use tensor HV instead of scalar coefficient
    disable_diagnostics, & ! use to disable diagnostics for timing reasons
    hypervis_order,       &
    hypervis_subcycle,    &
    hypervis_subcycle_tom,&
    hypervis_subcycle_q,  &
    smooth_phis_numcycle, &
    smooth_phis_p2filt, &
    smooth_phis_nudt,     &
    initial_total_mass,   & ! set > 0 to set the initial_total_mass
    u_perturb,     &        ! J&W baroclinic test perturbation size
    moisture,      &
    use_moisture,      &
    vfile_mid,     &
    vfile_int,     &
    vanalytic,     &
    vtop,          &
    precon_method, &
    maxits,        &
    tol,           &
    debug_level,   &
    theta_advect_form,   &
    vtheta_thresh,   &
    dp3d_thresh,   &
    pgrad_correction,    &
    hv_ref_profiles,     &
    hv_theta_correction, &
    hv_theta_thresh, &
    vert_remap_q_alg, &
    vert_remap_u_alg, &
    se_fv_phys_remap_alg, &
    timestep_make_subcycle_parameters_consistent


!PLANAR setup
#if !defined(CAM) && !defined(SCREAM)
  use control_mod, only:              &
    set_planar_defaults,&
    bubble_T0, &
    bubble_dT, &
    bubble_xycenter, &
    bubble_zcenter, &
    bubble_ztop, &
    bubble_xyradius, &
    bubble_zradius, &
    bubble_cosine, &
    bubble_moist, &
    bubble_moist_drh, &
    bubble_rh_background, &
    bubble_prec_type, &
    case_planar_bubble
#endif


! control parameters for dcmip stand-alone tests
#if !defined(CAM) && !defined(SCREAM)
  use control_mod, only:              &
    pertlim,                          &
    dcmip2_0_h0,                      &
    dcmip2_0_rm,                      &
    dcmip2_0_zetam,                   &
    dcmip2_x_ueq,                     &
    dcmip2_x_h0,                      &
    dcmip2_x_d,                       &
    dcmip2_x_xi,                      &
    dcmip4_moist,                     &
    dcmip4_X
#endif

  use thread_mod,     only: nthreads, omp_set_num_threads, omp_get_max_threads, vthreads
  use dimensions_mod, only: ne, ne_x, ne_y, np, nnodes, nmpi_per_node, npart, qsize, qsize_d, set_mesh_dimensions
#if defined(CAM) || defined(SCREAM)
  use time_mod,       only: tstep, nsplit, smooth
#else
  use time_mod,       only: tstep, ndays,nmax, nendstep,secpday, smooth, secphr, nsplit
#endif
  use parallel_mod,   only: parallel_t,  iam, abortmp, &
       partitionfornodes, mpireal_t, mpilogical_t, mpiinteger_t, mpichar_t
  use cg_mod,         only: cg_no_debug

#if !defined(CAM) && !defined(SCREAM)
  use interpolate_mod, only : vector_uvars, vector_vvars, max_vecvars, interpolate_analysis, replace_vec_by_vordiv
  use common_io_mod, only : &
       output_prefix,       &
       output_type,         &
       output_timeunits,    &
       output_start_time,   &
       output_end_time,     &
       output_frequency,    &
       output_dir,          &
       output_varnames1,    &
       output_varnames2,    &
       output_varnames3,    &
       output_varnames4,    &
       output_varnames5,    &
       max_output_variables,&
       max_output_streams,  &
       num_io_procs,        &
       io_stride,           &
       varname_len,         &
       infilenames,         &
       tool,                &
       MAX_INFILES
  use physical_constants, only: omega
  use common_movie_mod,   only : setvarnames
#endif

  use interpolate_mod, only : set_interp_parameter, get_interp_parameter


  !=======================================================================================================!
  ! This module should contain no global data and should only be used where readnl is called

  implicit none
  private


  public :: readnl

 contains

  ! ============================================
  ! readnl:
  ! read in the namelists...
  ! ============================================

#if defined(CAM) || defined(SCREAM)
  subroutine readnl(par, NLFileName)
    use shr_file_mod,      only: getunit=>shr_file_getUnit, freeunit=>shr_file_freeUnit
#ifndef HOMME_WITHOUT_PIOLIBRARY
    use mesh_mod, only : MeshOpen
#endif
    character(len=*), intent(in) :: NLFilename  ! namelist filename
#else
  subroutine readnl(par)
#ifndef HOMME_WITHOUT_PIOLIBRARY
    use mesh_mod, only : MeshOpen
#endif
#endif
    type (parallel_t), intent(in) ::  par
    character(len=MAX_FILE_LEN) :: mesh_file
    integer :: se_ftype, se_limiter_option
    integer :: se_nsplit
    integer :: interp_nlat, interp_nlon, interp_gridtype, interp_type
    integer :: i, ii, j
    integer  :: ierr
    character(len=80) :: errstr, arg
    real(kind=real_kind) :: dt_max, se_tstep
#if defined(CAM) || defined(SCREAM)
    character(len=MAX_STRING_LEN) :: se_topology
    integer :: se_partmethod
    integer :: se_ne
    integer :: unitn
    character(len=*), parameter ::  subname = "homme:namelist_mod"
#endif
    ! ============================================
    ! Namelists
    ! ============================================

    namelist /ctl_nl/ PARTMETHOD,                &         ! mesh partitioning method
                      COORD_TRANSFORM_METHOD,    &         ! Zoltan2 coordinate transformation method.
                      Z2_MAP_METHOD,             &         ! Zoltan2 processor mapping (network-topology aware) method.
                      TOPOLOGY,                  &         ! mesh topology
                      GEOMETRY,                  &         ! mesh geometry
#if defined(CAM) || defined(SCREAM)
      se_partmethod,     &
      se_topology,       &
      se_ne,             &
      se_limiter_option, &
#else
      qsize,             &         ! number of SE tracers
      nthreads,          &         ! number of threads per process
      limiter_option,    &
      smooth,            &         ! timestep Filter
      pertlim,           &         ! temperature initial perturbation
      omega,                   &   ! scaled rotation rate
      rearth,                  &   ! scaled earth radius
#endif
      COORD_TRANSFORM_METHOD, &
      Z2_MAP_METHOD,  &
      vthreads,      &             ! number of vertical/column threads per horizontal thread
      npart,         &
      numnodes,      &
      ne,            &             ! element resolution factor
      ne_x,            &             ! element resolution factor in x-dir for planar
      ne_y,            &             ! element resolution factor in y-dir for planar
      statefreq,     &             ! number of steps per printstate call
      integration,   &             ! integration method
      theta_hydrostatic_mode,       &   
      transport_alg , &      ! SE Eulerian, classical SL, cell-integrated SL
      semi_lagrange_cdr_alg, &
      semi_lagrange_cdr_check, &
      semi_lagrange_hv_q, &
      semi_lagrange_nearest_point_lev, &
      tstep_type,    &
      cubed_sphere_map, &
      qsplit,        &
      rsplit,        &
      dt_tracer_factor, &
      dt_remap_factor, &
      rk_stage_user, &
      LFTfreq,       &
      disable_diagnostics, &
      prescribed_wind, &
      se_ftype,        &           ! forcing type
      nu,            &
      nu_s,          &
      nu_q,          &
      nu_div,        &
      nu_p,          &
      nu_top,        &
      tom_sponge_start, &
      dcmip16_mu,     &
      dcmip16_mu_s,   &
      dcmip16_mu_q,   &
      dcmip16_prec_type,&
      dcmip16_pbl_type,&
      hypervis_order,    &
      hypervis_subcycle, &
      hypervis_subcycle_tom, &
      hypervis_subcycle_q, &
      hypervis_scaling, &
      smooth_phis_numcycle, &
      smooth_phis_p2filt, &
      smooth_phis_nudt, &
      initial_total_mass, &
      u_perturb,     &
      rotate_grid,   &
      mesh_file,     &               ! Name of mesh file
      theta_advect_form,     &
      vtheta_thresh,         &
      dp3d_thresh,         &
      pgrad_correction,      &
      hv_ref_profiles,       &
      hv_theta_correction,   &
      hv_theta_thresh,   &
      vert_remap_q_alg, &
      vert_remap_u_alg, &
      se_fv_phys_remap_alg


#if defined(CAM) || defined(SCREAM)
    namelist  /ctl_nl/ SE_NSPLIT,  &                ! number of dynamics steps per physics timestep
      se_tstep
#else
    namelist /ctl_nl/test_case,       &             ! test case idenitfier
      planar_slice,     &             ! const y-dir
      sub_case,        &             ! generic test case parameter
      nmax,            &             ! number of steps
      ndays,           &             ! number of days to simulate
      restartfreq,     &
      restartfile,     &             ! name of the restart file for INPUT
      restartdir,      &             ! name of the restart directory for OUTPUT
      runtype,         &
      tstep,           &             ! dynamics time step
      moisture
    ! control parameters for dcmip stand-alone tests
    namelist /ctl_nl/     &
      dcmip2_0_h0,        & !dcmip2-0 mountain height           (meters)
      dcmip2_0_rm,        & !dcmip2-0 mountain range radius     (radians)
      dcmip2_0_zetam,     & !dcmip2-0 mountain range half width (radians)
      dcmip2_x_ueq,       & !dcmip2-x windspeed at equator      (m/s)
      dcmip2_x_h0,        & !dcmip2-x mountain height           (m)
      dcmip2_x_d,         & !dcmip2-x mountain half-width       (m)
      dcmip2_x_xi,        & !dcmip2-x mountain wavelength       (m)
      dcmip4_moist,       & !dcmip4   moist, 0 or 1
      dcmip4_X              !dcmip4   scaling factor, nondim
!PLANAR
    namelist /ctl_nl/ &
      lx, ly, &
      sx, sy, &
      bubble_T0, &
      bubble_dT, &
      bubble_xycenter, &
      bubble_zcenter, &
      bubble_ztop, &
      bubble_xyradius, &
      bubble_zradius, &
      bubble_cosine, &
      bubble_moist, &
      bubble_moist_drh, &
      bubble_rh_background, &
      bubble_prec_type
    namelist /vert_nl/        &
      vfile_mid,          &
      vfile_int,          &
      vanalytic,          & ! use analytically generated vertical levels
      vtop                  ! top coordinate level. used when vanaltic=1
    namelist /analysis_nl/    &
      output_prefix,       &
      output_timeunits,    &
      output_start_time,   &
      output_end_time,     &
      output_frequency,    &
      output_dir,          &
      output_type,         &
      io_stride,           &
      num_io_procs,        &
      infilenames,         &
      replace_vec_by_vordiv, &
      vector_uvars,        &
      vector_vvars,        &
      output_varnames1,    &
      output_varnames2,    &
      output_varnames3,    &
      output_varnames4,    &
      output_varnames5,    &
      interp_nlat,         &
      interp_nlon,         &
      interp_lon0,         &
      interp_gridtype,     &
      interp_type,         &
      interpolate_analysis, &
      tool
#endif
! ^ ifndef CAM


    namelist /solver_nl/precon_method, &
      maxits,        &
      tol,           &
      debug_level

#if (defined MODEL_THETA_L && defined ARKODE)
    namelist /arkode_nl/ &
      rel_tol, &
      abs_tol, &
      calc_nonlinear_stats, &
      use_column_solver
#endif


    ! ==========================
    ! Set the default partmethod
    ! ==========================
    PARTMETHOD    = SFCURVE
    COORD_TRANSFORM_METHOD = SPHERE_COORDS
    Z2_MAP_METHOD = Z2_NO_TASK_MAPPING
    npart         = 1
    se_tstep=-1
#if !defined(CAM) && !defined(SCREAM)
    ndays         = 0
    nmax          = 12
    nthreads = 1
    se_ftype = ftype   ! MNL: For non-CAM runs, ftype=0 in control_mod
    nsplit = 1
    pertlim = 0.0_real_kind
#endif
    sub_case      = 1
    numnodes      = -1
    restartfreq   = -100
    restartdir    = "./restart/"
    runtype       = 0
    statefreq     = 1
    integration   = "explicit"
    moisture      = "dry"
    nu_top=0
    initial_total_mass=0
    mesh_file='none'
    ne              = 0
    ne_x              = 0
    ne_y              = 0
    transport_alg = 0
    semi_lagrange_cdr_alg = 3
    semi_lagrange_cdr_check = .false.
    semi_lagrange_hv_q = 1
    semi_lagrange_nearest_point_lev = 256
    disable_diagnostics = .false.
    se_fv_phys_remap_alg = 1
    planar_slice = .false.

    theta_hydrostatic_mode = .true.    ! for preqx, this must be .true.
#if ( defined MODEL_THETA_C || defined MODEL_THETA_L ) 
    theta_hydrostatic_mode = .false.   ! default NH
#endif


    ! =======================
    ! Read namelist variables
    ! =======================

!   write(iulog,*) "masterproc=",par%masterproc

    if (par%masterproc) then
       write(iulog,*)"reading ctl namelist..."
#if defined(CAM) || defined(SCREAM)
       unitn=getunit()
       open( unitn, file=trim(nlfilename), status='old' )
       ierr = 1
       do while ( ierr /= 0 )
          read (unitn,ctl_nl,iostat=ierr)
          if (ierr < 0) then
            write(6,*) 'ierr =',ierr
             call abortmp( subname//':: namelist read returns an'// &
                  ' end of file or end of record condition' )
          end if
       end do
       close( unitn )
       call freeunit( unitn )
#elif defined(OSF1) || defined(_NAMELIST_FROM_FILE)
       open(unit=7,file="input.nl",status="OLD")
       read(unit=7,nml=ctl_nl)
#else
      print *, 'Reading namelist ctl_nl from standard input'
      read(*,nml=ctl_nl)
#endif
#ifndef _USEMETIS
      ! override METIS options to SFCURVE
      if (partmethod>=0 .and. partmethod<=3) partmethod=SFCURVE
#endif
       ! ========================
       ! if this is a restart run
       ! ========================
       if(runtype .eq. 1) then
          write(iulog,*)"readnl: restartfile = ",restartfile
       else if(runtype < 0) then
          write(iulog,*)'readnl: runtype=', runtype,' interpolation mode '
       endif



       ! ================================================
       ! if ndays is defined in the namelist, then
       ! moviefreq and restartfreq are interpreted to be in units of days.
       ! Both must be converted to numbers of steps.
       ! ================================================
#if !defined(CAM) && !defined(SCREAM)
       if (tstep <= 0) then
          call abortmp('tstep must be > 0')
       end if
       if (ndays>0) then
          nmax = ndays * (secpday/tstep)
          if (restartfreq>0) restartfreq=restartfreq*(secpday/tstep)
       end if
       nEndStep = nmax
#endif

       if((integration .ne. "explicit").and.(integration .ne. "runge_kutta").and. &
                    (integration .ne. "full_imp")) then
          call abortmp('integration must be explicit, full_imp, or runge_kutta')
       end if

       if (integration == "full_imp") then
          if (tstep_type<10) then
             ! namelist did not set a valid tstep_type. pick one:
             tstep_type=11   ! backward euler
             !tstep_type=12  ! BDF2 with BE bootstrap
          endif
       endif

#if !defined(CAM) && !defined(SCREAM)

       !checks before the next NL
       !check on ne
       if (topology == "plane" .and. ne /= 0) then
          call abortmp('cannot set ne for planar topology, use ne_x and ne_y instead')
       end if

       !setting default PLANAR values if they are not set in ctl_nl
       call set_planar_defaults()

       !checks on planar tests
       if (test_case(1:7)=="planar_") then

       !check on Lx, Ly
       if(lx .le. 0.d0 .or. ly .le. 0.d0) then
          call abortmp('for planar tests, planar_lx and planar_ly should be >0')
       endif

       endif


#ifdef _PRIM
       write(iulog,*) "reading physics namelist..."
       if (test_case      == "held_suarez0"   .or. &
           test_case(1:10)== "baroclinic"     .or. &
           test_case(1:13)== "jw_baroclinic"  .or. &
           test_case(1:5) == "dcmip"          .or. &
           test_case(1:5) == "mtest"          .or. &
           test_case      == "planar_hydro_gravity_wave"            .or. &
           test_case      == "planar_nonhydro_gravity_wave"           .or. &
           test_case      == "planar_hydro_mtn_wave"            .or. &
           test_case      == "planar_nonhydro_mtn_wave"           .or. &
           test_case      == "planar_schar_mtn_wave"            .or. &
           test_case      == "planar_rising_bubble"             .or. &
           test_case      == "planar_density_current"             .or. &
           test_case      == "planar_baroclinic_instab"             .or. &
           test_case      == "planar_moist_rising_bubble"            .or. &
           test_case      == "planar_moist_density_current"            .or. &
           test_case      == "planar_moist_baroclinic_instab"            .or. &
           test_case      == "planar_tropical_cyclone"             .or. &
           test_case      == "planar_supercell"             .or. &
           test_case(1:4) == "asp_")  then
         write(iulog,*) "reading vertical namelist..."

#if defined(OSF1) || defined(_NAMELIST_FROM_FILE)
         read(unit=7,nml=vert_nl)
#else
         read(*,nml=vert_nl)
#endif
         vfile_mid  = trim(adjustl(vfile_mid))
         vfile_int  = trim(adjustl(vfile_int))

         write(iulog,*) '  vfile_mid=',trim(vfile_mid)
         write(iulog,*) '  vfile_int=',trim(vfile_int)
         write(iulog,*) '  vanalytic=',vanalytic
         if(vanalytic==1) then
         write(iulog,*) '  vtop=',vtop
         endif

       end if
#endif


!      Default interpolation grid  (0 = auto compute based on ne,nv)  interpolation is off by default
       interp_nlat =  0
       interp_nlon = 0
       interp_lon0 = 0
       interp_gridtype = 2
       interp_type = 0
       replace_vec_by_vordiv(:)=.false.
       vector_uvars(:)=''
       vector_vvars(:)=''
       vector_uvars(1:12) = (/ &
            'U         ','UBOT      ','U200      ','U250      ','U850      ','FU        ',&
            'CONVU     ','DIFFU     ','UTGWORO   ','UFLX      ','MET_U     ','MET_U_tend' /)
       vector_vvars(1:12) = (/ &
            'V         ','VBOT      ','V200      ','V250      ','V850      ','FV        ',&
            'CONVV     ','DIFFV     ','VTGWORO   ','VFLX      ','MET_V     ','MET_V_tend' /)
       infilenames(:) = ''
       output_prefix = ""
       output_start_time=0
       output_end_time =0
       output_frequency=0
       output_timeunits=0
! Default is to write all variables at time 0 and time nendstep.
       output_end_time(1) = -1
       output_frequency(1) = nendstep
       output_dir = "./movies/"
       output_varnames1='     '
       output_varnames2='     '
       output_varnames3='     '
       output_varnames4='     '
       output_varnames5='     '
       io_stride=0
       num_io_procs=0
       output_type = 'netcdf' ! Change by MNL
!     output_type = 'pnetcdf'
       tool = 'none'


#ifndef HOMME_WITHOUT_PIOLIBRARY
       write(iulog,*)"reading analysis namelist..."
#if defined(OSF1) || defined(_NAMELIST_FROM_FILE)
       read(unit=7,nml=analysis_nl)
#else
       read(*,nml=analysis_nl)
#endif

      if (io_stride .eq.0 .and. num_io_procs .eq.0) then
         ! user did not set anything
         io_stride=1
         num_io_procs = par%nprocs
      else if (io_stride.eq.0) then  ! user set num_io_procs
         io_stride = max(1,par%nprocs/num_io_procs)
      else if (num_io_procs .eq. 0) then   ! user set io_stride
         num_io_procs=max(1,par%nprocs/io_stride)
      else ! user set both parameters
      endif
      ! sanity check
      if(num_io_procs*io_stride>par%nprocs) then
         if (io_stride > par%nprocs) io_stride=par%nprocs
         num_io_procs=par%nprocs/io_stride
      end if

      if(output_varnames1(1).eq.'     ') then
         if( runtype>=0) then
            call setvarnames(output_varnames1)
         else ! interpolation mode
            output_varnames1='all'
         end if
      end if
       do i=1,max_output_streams
          if(output_frequency(i)>0 .and. output_end_time(i)==0) output_end_time(i)=-1
          if(output_timeunits(i).eq.1) then  ! per_day
             output_frequency(i) = output_frequency(i)*(secpday/tstep)
             output_start_time(i)= output_start_time(i)*(secpday/tstep)
             if (output_end_time(i)>=0) &
                  output_end_time(i)  = output_end_time(i)*(secpday/tstep)
          else if(output_timeunits(i).eq.2) then  ! per_hour
             output_frequency(i) = output_frequency(i)*(secphr/tstep)
             output_start_time(i)= output_start_time(i)*(secphr/tstep)
             if (output_end_time(i)>=0) &
                  output_end_time(i)  = output_end_time(i)*(secphr/tstep)
          else if(output_timeunits(i).eq.3) then  ! per_seconds
             output_frequency(i) = output_frequency(i)/tstep
             output_start_time(i)= output_start_time(i)/tstep
             if (output_end_time(i)>=0) &
                  output_end_time(i)  = output_end_time(i)/tstep
          end if
          if(output_end_time(i)<0) output_end_time(i)=nEndStep
          if ( output_start_time(i) > output_end_time(i) ) output_frequency(i)=0
       end do
#endif

#if (defined MODEL_THETA_L && defined ARKODE)
       write(iulog,*)"reading arkode namelist..."
#if defined(OSF1) || defined(_NAMELIST_FROM_FILE)
       read(unit=7,nml=arkode_nl)
#else
       read(*,nml=arkode_nl)
#endif
#endif

!=======================================================================================================!

#if defined(OSF1) || defined(_NAMELIST_FROM_FILE)
       close(unit=7)
#endif
#endif
! ^ if !defined(CAM) && !defined(SCREAM)
       ierr = timestep_make_subcycle_parameters_consistent(par, rsplit, qsplit, &
            dt_remap_factor, dt_tracer_factor)

#if defined(CAM) || defined(SCREAM)
       limiter_option=se_limiter_option
       partmethod = se_partmethod
       ne         = se_ne
       topology   = se_topology
       qsize      = qsize_d
       nsplit     = se_nsplit
       tstep      = se_tstep
       if (tstep > 0) then
          if (par%masterproc .and. nsplit > 0) then
             write(iulog,'(a,i3,a)') &
                  'se_tstep and se_nsplit were specified; changing se_nsplit from ', &
                  nsplit, ' to -1.'
          end if
          nsplit = -1
       end if
#endif
    end if


    call MPI_barrier(par%comm,ierr)

    npart  = par%nprocs

    ! Broadcast namelist variables to all MPI processes

    call MPI_bcast(Z2_MAP_METHOD ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(COORD_TRANSFORM_METHOD ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(PARTMETHOD ,     1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(TOPOLOGY,        MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
    call MPI_bcast(geometry,        MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
    call MPI_bcast(test_case,       MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
    call MPI_bcast(ne,              1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(ne_x,              1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(ne_y,              1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(qsize,           1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(sub_case,        1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(statefreq,       1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(restartfreq,     1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(runtype,         1,MPIinteger_t,par%root,par%comm,ierr)

#if !defined(CAM) && !defined(SCREAM)
    if(test_case == "dcmip2012_test4") then
       rearth = rearth/dcmip4_X
       omega = omega*dcmip4_X
    endif

    call MPI_bcast(pertlim,         1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(nmax,            1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(NTHREADS,        1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(ndays,           1, MPIinteger_t, par%root,par%comm,ierr)
    nEndStep = nmax

    call MPI_bcast(omega,           1, MPIreal_t,    par%root,par%comm,ierr)
    call MPI_bcast(rearth,          1, MPIreal_t,    par%root,par%comm,ierr)
    rrearth = 1.0d0/rearth

    call MPI_bcast(dcmip2_0_h0,     1, MPIreal_t,    par%root,par%comm,ierr)
    call MPI_bcast(dcmip2_0_rm,     1, MPIreal_t,    par%root,par%comm,ierr)
    call MPI_bcast(dcmip2_0_zetam,  1, MPIreal_t,    par%root,par%comm,ierr)
    call MPI_bcast(dcmip2_x_ueq,    1, MPIreal_t,    par%root,par%comm,ierr)
    call MPI_bcast(dcmip2_x_h0,     1, MPIreal_t,    par%root,par%comm,ierr)
    call MPI_bcast(dcmip2_x_d,      1, MPIreal_t,    par%root,par%comm,ierr)
    call MPI_bcast(dcmip2_x_xi,     1, MPIreal_t,    par%root,par%comm,ierr)
    call MPI_bcast(dcmip4_moist,    1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(dcmip4_X,        1, MPIreal_t,    par%root,par%comm,ierr)
#endif
    call MPI_bcast(vthreads  ,      1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(smooth,          1, MPIreal_t,    par%root,par%comm,ierr)
    call MPI_bcast(NSPLIT,          1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(tstep,           1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(limiter_option,  1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(se_ftype,        1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(theta_advect_form,1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(vtheta_thresh,    1, MPIreal_t, par%root,par%comm,ierr)
    call MPI_bcast(dp3d_thresh,    1, MPIreal_t, par%root,par%comm,ierr)
    call MPI_bcast(pgrad_correction,   1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(hv_ref_profiles,    1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(hv_theta_correction,1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(hv_theta_thresh,1, MPIreal_t, par%root,par%comm,ierr)
    call MPI_bcast(vert_remap_q_alg,1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(vert_remap_u_alg,1, MPIinteger_t, par%root,par%comm,ierr)

    call MPI_bcast(nu,              1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(nu_s,            1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(nu_q,            1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(nu_div,          1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(nu_p,            1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(nu_top,          1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(tom_sponge_start,1, MPIreal_t   , par%root,par%comm,ierr)

    call MPI_bcast(dcmip16_mu,      1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(dcmip16_mu_s,    1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(dcmip16_mu_q,    1, MPIreal_t   , par%root,par%comm,ierr)

    call MPI_bcast(dcmip16_prec_type, 1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(dcmip16_pbl_type , 1, MPIinteger_t, par%root,par%comm,ierr)

    call MPI_bcast(disable_diagnostics,1,MPIlogical_t,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_order,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_scaling,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_subcycle,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_subcycle_tom,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_subcycle_q,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(smooth_phis_numcycle,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(smooth_phis_p2filt,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(smooth_phis_nudt,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(initial_total_mass ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(u_perturb     ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(rotate_grid   ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(integration,MAX_STRING_LEN,MPIChar_t ,par%root,par%comm,ierr)
#ifndef HOMME_WITHOUT_PIOLIBRARY
    call MPI_bcast(mesh_file,MAX_FILE_LEN,MPIChar_t ,par%root,par%comm,ierr)
#endif

!PLANAR
    call MPI_bcast(lx     ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(ly     ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(sx     ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(sy     ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(planar_slice ,1,MPIlogical_t,par%root,par%comm,ierr)

!PLANAR
#if !defined(CAM) && !defined (SCREAM)
    call MPI_bcast(bubble_T0 ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(bubble_dT ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(bubble_xycenter,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(bubble_zcenter ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(bubble_ztop ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(bubble_xyradius ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(bubble_zradius ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(bubble_cosine ,1,MPIlogical_t,par%root,par%comm,ierr)
    call MPI_bcast(bubble_moist ,1,MPIlogical_t,par%root,par%comm,ierr)
    call MPI_bcast(bubble_moist_drh ,1,MPIreal_t,par%root,par%comm,ierr)
    call MPI_bcast(bubble_rh_background ,1,MPIreal_t,par%root,par%comm,ierr)
    call MPI_bcast(bubble_prec_type, 1, MPIinteger_t, par%root,par%comm,ierr)

    call MPI_bcast(case_planar_bubble,1,MPIlogical_t,par%root,par%comm,ierr)
#endif

    call MPI_bcast(theta_hydrostatic_mode ,1,MPIlogical_t,par%root,par%comm,ierr)
    call MPI_bcast(transport_alg ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(semi_lagrange_cdr_alg ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(semi_lagrange_cdr_check ,1,MPIlogical_t,par%root,par%comm,ierr)
    call MPI_bcast(semi_lagrange_hv_q ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(semi_lagrange_nearest_point_lev ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(tstep_type,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(cubed_sphere_map,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(qsplit,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(rsplit,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(dt_tracer_factor,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(dt_remap_factor,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(rk_stage_user,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(LFTfreq,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(prescribed_wind,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(moisture,MAX_STRING_LEN,MPIChar_t ,par%root,par%comm,ierr)
    call MPI_bcast(se_fv_phys_remap_alg,1,MPIinteger_t ,par%root,par%comm,ierr)

    call MPI_bcast(restartfile,MAX_STRING_LEN,MPIChar_t ,par%root,par%comm,ierr)
    call MPI_bcast(restartdir,MAX_STRING_LEN,MPIChar_t ,par%root,par%comm,ierr)

    if (integration == "full_imp") then
       call MPI_bcast(precon_method,MAX_STRING_LEN,MPIChar_t,par%root,par%comm,ierr)
       call MPI_bcast(maxits     ,1,MPIinteger_t,par%root,par%comm,ierr)
       call MPI_bcast(tol        ,1,MPIreal_t   ,par%root,par%comm,ierr)
    end if

    call MPI_bcast(vfile_mid, MAX_STRING_LEN, MPIChar_t   , par%root, par%comm,ierr)
    call MPI_bcast(vfile_int, MAX_STRING_LEN, MPIChar_t   , par%root, par%comm,ierr)
    call MPI_bcast(vanalytic, 1,              MPIinteger_t, par%root, par%comm,ierr)
    call MPI_bcast(vtop     , 1,              MPIreal_t   , par%root, par%comm,ierr)

#if !defined(CAM) && !defined(SCREAM)
#ifndef HOMME_WITHOUT_PIOLIBRARY
    call MPI_bcast(output_prefix,MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
    call MPI_bcast(output_timeunits ,max_output_streams,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(output_start_time ,max_output_streams,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(output_end_time ,max_output_streams,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(output_frequency ,max_output_streams,MPIinteger_t,par%root,par%comm,ierr)

    call MPI_bcast(output_dir   ,MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
    call MPI_bcast(output_varnames1 ,varname_len*max_output_variables,MPIChar_t,par%root,par%comm,ierr)
    call MPI_bcast(output_varnames2 ,varname_len*max_output_variables,MPIChar_t,par%root,par%comm,ierr)
    call MPI_bcast(output_varnames3 ,varname_len*max_output_variables,MPIChar_t,par%root,par%comm,ierr)
    call MPI_bcast(output_varnames4 ,varname_len*max_output_variables,MPIChar_t,par%root,par%comm,ierr)
    call MPI_bcast(output_varnames5 ,varname_len*max_output_variables,MPIChar_t,par%root,par%comm,ierr)
    call MPI_bcast(io_stride , 1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(num_io_procs , 1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(output_type , 9,MPIChar_t,par%root,par%comm,ierr)
    call MPI_bcast(infilenames ,160*MAX_INFILES ,MPIChar_t,par%root,par%comm,ierr)
    call MPI_bcast(tool,MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
#endif

#if (defined MODEL_THETA_L && defined ARKODE)
    call MPI_bcast(rel_tol, 1, MPIreal_t, par%root, par%comm, ierr)
    call MPI_bcast(abs_tol, 1, MPIreal_t, par%root, par%comm, ierr)
    call MPI_bcast(calc_nonlinear_stats, 1, MPIlogical_t, par%root, par%comm, ierr)
    call MPI_bcast(use_column_solver, 1, MPIlogical_t, par%root, par%comm, ierr)
#endif

    ! should we assume Q(:,:,:,1) has water vapor:
    use_moisture = ( moisture /= "dry") 
    if (qsize<1) use_moisture = .false.  

    if(par%masterproc) print *, "use moisture:", use_moisture


    ! use maximum available:
    if (NThreads == -1) then
#if defined(HORIZ_OPENMP) || defined (COLUMN_OPENMP)
       NThreads = omp_get_max_threads()
#else
       NThreads = 1  
#endif       
    endif
    
    ! sanity check on thread count
    ! HOMME will run if if nthreads > max, but gptl will print out GB of warnings.
    if (NThreads > omp_get_max_threads()) then
       if(par%masterproc) print *, "Main:NThreads=",NThreads
       if(par%masterproc) print *, 'omp_get_max_threads() = ',OMP_get_max_threads()
       if(par%masterproc) print *, 'requested threads exceeds OMP_get_max_threads()'
       call abortmp('stopping')
    endif
    call omp_set_num_threads(NThreads)

    ! if user sets hypervis_subcycle=-1, then use automatic formula
    if (hypervis_subcycle==-1) then
       if (np==4 .and. topology == "cube") then
          ! 1.25d23 worked out by testing, for nv==4
          ! a little confusing:
          ! u,v:  nu and hypervis_subcycle
          ! T:    nu_s and hypervis_subcycle
          ! Q:    nu and hypervis_subcycle_q
          dt_max = 1.25d23/(nu*ne**4.0)
          hypervis_subcycle_q = ceiling( tstep/dt_max )
          hypervis_subcycle   = ceiling( tstep/dt_max )
       else
          call abortmp('hypervis_subcycle auto determine only supported for nv==4 and topology==cube')
       endif
    endif
#endif
    ! set defautl for dynamics remap
    if (vert_remap_u_alg == -2) vert_remap_u_alg = vert_remap_q_alg
    
    ! more thread error checks:  
#ifdef HORIZ_OPENMP
    if(par%masterproc) write(iulog,*)'-DHORIZ_OPENMP enabled'
#else
    if(par%masterproc) write(iulog,*)'-DHORIZ_OPENMP disabled'
#endif
#ifdef COLUMN_OPENMP
    if(par%masterproc) write(iulog,*)'-DCOLUMN_OPENMP enabled'
#else
    if(par%masterproc) write(iulog,*)'-DCOLUMN_OPENMP disabled'
#endif

if (topology == "plane" .and. mesh_file /= "none") then
  call abortmp("RRM grids not yet supported for plane")
end if

    if (ne /=0 .or. ne_x /=0 .or. ne_y /=0) then
    if (mesh_file /= "none" .and. mesh_file /= "/dev/null") then
      write (*,*) "namelist_mod: mesh_file:",trim(mesh_file), &
                  " and ne/ne_x/ne_y:",ne,ne_x,ne_y," are both specified in the input file."
      write (*,*) "Specify one or the other, but not both."
      call abortmp("Do not specify ne (or ne_x, ne_y) if using a mesh file input.")
    end if
    end if
    if (par%masterproc) write (iulog,*) "Mesh File:", trim(mesh_file)
    if (ne.eq.0 .and. ne_x .eq. 0 .and. ne_y .eq. 0) then
#ifndef HOMME_WITHOUT_PIOLIBRARY
       call set_mesh_dimensions()
       if (par%masterproc) write (iulog,*) "Opening Mesh File:", trim(mesh_file)
       call MeshOpen(mesh_file, par)
#else
       call abortmp("Build is without PIO library, mesh runs (ne=0) are not supported.")
#endif
    end if
    ! set map
    if (cubed_sphere_map<0) then
#if ( defined MODEL_THETA_C || defined MODEL_THETA_L ) 
       cubed_sphere_map=2  ! theta model default = element local
#else
       cubed_sphere_map=0  ! default is equi-angle gnomonic
#endif
    endif
    if (ne.eq.0 .and. ne_x .eq. 0 .and. ne_y .eq. 0) cubed_sphere_map=2  ! must use element_local for var-res grids
    if (par%masterproc) write (iulog,*) "Reference element projection: cubed_sphere_map=",cubed_sphere_map

! set geometric factors
! Ideally this stuff moves into a separate sub routine
! Along with all the other things that are test case specific (rearth/Omega scaling for small-earth experiments, etc.)
    if (geometry == "plane") then
      scale_factor = 1.0D0
      scale_factor_inv = 1.0D0
      laplacian_rigid_factor = 0.0D0 !this eliminates the correction to ensure the Laplacian doesn't damp rigid motion

! makes the y-direction cells identical in size to the x-dir cells
! this is important for hyperviscosity, etc.
! Also adjusts Sy so y-dir domain is centered at 0
      if (planar_slice .eqv. .true.) then
        Ly = (Lx / ne_x) * ne_y
        Sy = -Ly/2.0D0
      end if

      domain_size = Lx * Ly
      dx = Lx/ne_x
      dy = Ly/ne_y
      dx_ref = 1.0D0/ne_x
      dy_ref = 1.0D0/ne_y

    else if (geometry == "sphere") then
      scale_factor = rearth
      scale_factor_inv = rrearth
      domain_size = 4.0D0*DD_PI
      laplacian_rigid_factor = rrearth
    end if ! if plane

    if (topology == "plane" .and. hypervis_scaling==0) then
      call abortmp("Error: planar grids require the use of tensor HV")
    end if

    ftype = se_ftype

#ifdef _PRIM
    if (limiter_option==8 .or. limiter_option==84 .or. limiter_option == 9) then
       if (hypervis_subcycle_q/=1 .and. transport_alg == 0) then
          call abortmp('limiter 8,84,9 require hypervis_subcycle_q=1')
       endif
    endif
    if (transport_alg == 0 .and. dt_remap_factor > 0 .and. dt_remap_factor < dt_tracer_factor) then
       call abortmp('Only SL transport supports vertical remap time step < tracer time step.')
    end if
#endif

#if !defined(CAM) && !defined(SCREAM)
!standalone homme does not support ftype=1 (cause it is identical to ftype=0).
!also, standalone ftype=0 is the same as standalone ftype=2.
    if ((ftype == 0).or.(ftype == 2).or.(ftype == 3).or.(ftype == 4).or.(ftype == -1)) then
    else
       call abortmp('Standalone homme supports only se_ftype=-1,0,2,3,4')
    endif
#endif

    if((prescribed_wind/=0).and.(prescribed_wind/=1))then
          call abortmp('prescribed_wind should be either 0 or 1')
    endif

!=======================================================================================================!
#if defined(CAM) || defined(SCREAM)
    nmpi_per_node=1
#endif
#if !defined(CAM) && !defined(SCREAM)
    call MPI_bcast(interpolate_analysis, 7,MPIlogical_t,par%root,par%comm,ierr)
    call MPI_bcast(interp_nlat , 1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(interp_nlon , 1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(interp_lon0 , 1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(interp_gridtype , 1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(interp_type , 1,MPIinteger_t,par%root,par%comm,ierr)

    call MPI_bcast(replace_vec_by_vordiv ,MAX_VECVARS ,MPIlogical_t,par%root,par%comm,ierr)
    call MPI_bcast(vector_uvars ,10*MAX_VECVARS ,MPIChar_t,par%root,par%comm,ierr)
    call MPI_bcast(vector_vvars ,10*MAX_VECVARS ,MPIChar_t,par%root,par%comm,ierr)

    call set_interp_parameter('gridtype',interp_gridtype)
    call set_interp_parameter("itype",interp_type)
    if(any(interpolate_analysis)) then
! FIX: THIS MUST CHANGE FOR PLANAR INTERPOLATION, SINCE WE DONT HAVE NE ANYMORE...
       if (interp_nlon==0 .or. interp_nlat==0) then
          ! compute interpolation grid based on number of points around equator
          call set_interp_parameter('auto',4*ne*(np-1))
          interp_nlat = get_interp_parameter('nlat')
          interp_nlon = get_interp_parameter('nlon')
       else
          call set_interp_parameter('nlat',interp_nlat)
          call set_interp_parameter('nlon',interp_nlon)
       endif
    endif
#endif
! ^ ifndef CAM

    ! some default diffusion coefficiets
    if(nu_s<0)    nu_s  = nu
    if(nu_q<0)    nu_q  = nu
    if(nu_div<0)  nu_div= nu
    if(nu_p<0) then                                                                           
       if (rsplit==0) then                                                                    
          nu_p=0  ! eulerian code traditionally run with nu_p=0                               
       else                                                                                   
          nu_p=nu                                                                             
       endif                                                                                  
    endif 
    if(dcmip16_mu_s<0)    dcmip16_mu_s  = dcmip16_mu
    if(dcmip16_mu_q<0)    dcmip16_mu_q  = dcmip16_mu_s

    nnodes = npart/nmpi_per_node
    if(numnodes > 0 ) then
        nnodes = numnodes
        nmpi_per_node = npart/nnodes
    endif

    ! ====================================================================
    !  Do not perform node level partitioning if you are only on one node
    ! ====================================================================
    ! PartitionForNodes=.FALSE.
    if((nnodes .eq. 1) .and. PartitionForNodes) PartitionForNodes=.FALSE.

    if (par%masterproc) then
       write(iulog,*)"done reading namelist..."

       write(iulog,*)"readnl: topology      = ",TRIM( TOPOLOGY )
       write(iulog,*)"readnl: geometry      = ",TRIM( geometry )
#if !defined(CAM) && !defined(SCREAM)
       write(iulog,*)"readnl: test_case     = ",TRIM(test_case)
       write(iulog,*)"readnl: omega         = ",omega
       write(iulog,*)"readnl: sub_case      = ",sub_case
       write(iulog,*)"readnl: ndays         = ",ndays
       write(iulog,*)"readnl: nmax          = ",nmax
       write(iulog,*)"readnl: pertlim       = ",pertlim

       write(iulog,*)"readnl: qsize,qsize_d = ",qsize,qsize_d
       if (qsize>qsize_d) then
          call abortmp('user specified qsize > qsize_d parameter in dimensions_mod.F90')
       endif
       write(iulog,*)"readnl: NThreads      = ",NTHREADS
#endif

if (topology == "cube") then
       write(iulog,*)"readnl: ne,np         = ",ne,np
else if (topology == "plane") then
  write(iulog,*)"readnl: ne_x,ne_y,np         = ",ne_x,ne_y,np
end if
       write(iulog,*)"readnl: partmethod    = ",PARTMETHOD
       write(iulog,*)"readnl: COORD_TRANSFORM_METHOD    = ",COORD_TRANSFORM_METHOD
       write(iulog,*)"readnl: Z2_MAP_METHOD    = ",Z2_MAP_METHOD

       write(iulog,*)'readnl: nmpi_per_node = ',nmpi_per_node
       write(iulog,*)"readnl: vthreads      = ",vthreads
       write(iulog,*)'readnl: nnodes        = ',nnodes
       write(iulog,*)'readnl: npart         = ',npart

       print *
       write(iulog,*)"readnl: integration   = ",trim(integration)
       if (integration == "explicit" ) then
          write(iulog,*)"readnl: LF-trapazoidal freq= ",LFTfreq
       endif
       if (integration == "runge_kutta"  ) then
          write(iulog,*)"readnl: rk_stage_user   = ",rk_stage_user
       endif
       write(iulog,*)"readnl: theta_hydrostatic_mode = ",theta_hydrostatic_mode
       write(iulog,*)"readnl: transport_alg   = ",transport_alg
       write(iulog,*)"readnl: semi_lagrange_cdr_alg   = ",semi_lagrange_cdr_alg
       write(iulog,*)"readnl: semi_lagrange_cdr_check   = ",semi_lagrange_cdr_check
       write(iulog,*)"readnl: semi_lagrange_hv_q   = ",semi_lagrange_hv_q
       write(iulog,*)"readnl: semi_lagrange_nearest_point_lev   = ",semi_lagrange_nearest_point_lev
       write(iulog,*)"readnl: tstep_type    = ",tstep_type
       write(iulog,*)"readnl: theta_advect_form = ",theta_advect_form
       write(iulog,*)"readnl: vtheta_thresh     = ",vtheta_thresh
       write(iulog,*)"readnl: dp3d_thresh     = ",dp3d_thresh
       write(iulog,*)"readnl: pgrad_correction  = ",pgrad_correction
       write(iulog,*)"readnl: hv_ref_profiles   = ",hv_ref_profiles
       write(iulog,*)"readnl: hv_theta_correction= ",hv_theta_correction
       write(iulog,*)"readnl: hv_theta_thresh   = ",hv_theta_thresh
       if (hv_ref_profiles==0 .and. hv_theta_correction==1) then
          call abortmp("hv_theta_correction=1 requires hv_ref_profiles=1 or 2")
       endif
       
       write(iulog,*)"readnl: vert_remap_q_alg  = ",vert_remap_q_alg
       write(iulog,*)"readnl: vert_remap_u_alg  = ",vert_remap_u_alg
#if defined(CAM) || defined(SCREAM)
       write(iulog,*)"readnl: se_nsplit         = ", NSPLIT
       write(iulog,*)"readnl: se_tstep         = ", tstep
       write(iulog,*)"readnl: se_ftype          = ",ftype
       write(iulog,*)"readnl: se_limiter_option = ",limiter_option
#else
       write(iulog,*)"readnl: tstep          = ",tstep
       write(iulog,*)"readnl: ftype          = ",ftype
       write(iulog,*)"readnl: limiter_option = ",limiter_option
#endif
       write(iulog,*)"readnl: qsplit (deprecated) = ",qsplit
       write(iulog,*)"readnl: vertical remap frequency rsplit (0=disabled) (deprecated): ",rsplit
       write(iulog,*)"readnl: dt_tracer_factor = ",dt_tracer_factor
       write(iulog,*)"readnl: vertical remap frequency dt_remap_factor (0=disabled): ",dt_remap_factor

       write(iulog,*)"readnl: runtype       = ",runtype
       write(iulog,*)"readnl: se_fv_phys_remap_alg = ",se_fv_phys_remap_alg

       if(hypervis_scaling /=0)then
          write(iulog,*)"Tensor hyperviscosity:  hypervis_scaling=",hypervis_scaling
       else
          write(iulog,*)"Constant (hyper)viscosity used."
       endif

       write(iulog,*)"hypervis_subcycle     = ",hypervis_subcycle
       write(iulog,*)"hypervis_subcycle_tom = ",hypervis_subcycle_tom
       write(iulog,*)"hypervis_subcycle_q   = ",hypervis_subcycle_q
       write(iulog,'(a,2e9.2)')"viscosity:  nu (vor/div) = ",nu,nu_div
       write(iulog,'(a,2e9.2)')"viscosity:  nu_s      = ",nu_s
       write(iulog,'(a,2e9.2)')"viscosity:  nu_q      = ",nu_q
       write(iulog,'(a,2e9.2)')"viscosity:  nu_p      = ",nu_p
       write(iulog,'(a,2e9.2)')"viscosity:  nu_top      = ",nu_top
       write(iulog,'(a,2e9.2)')"viscosity:  tom_sponge_start  = ",tom_sponge_start

       if(dcmip16_mu/=0)  write(iulog,'(a,2e9.2)')"1st order viscosity:  dcmip16_mu   = ",dcmip16_mu
       if(dcmip16_mu_s/=0)write(iulog,'(a,2e9.2)')"1st order viscosity:  dcmip16_mu_s = ",dcmip16_mu_s
       if(dcmip16_mu_q/=0)write(iulog,'(a,2e9.2)')"1st order viscosity:  dcmip16_mu_q = ",dcmip16_mu_q

       if(initial_total_mass>0) then
          write(iulog,*) "initial_total_mass = ",initial_total_mass
       end if

#if !defined(CAM) && !defined(SCREAM)
#ifndef HOMME_WITHOUT_PIOLIBRARY
       write(iulog,*)"  analysis: output_prefix = ",TRIM(output_prefix)
       write(iulog,*)"  analysis: io_stride = ",io_stride
       write(iulog,*)"  analysis: num_io_procs = ",num_io_procs

       do i=1,max_output_streams
          if(output_frequency(i) .gt. 0) then
             write(iulog,*)"  analysis stream     :",i
             write(iulog,*)"  analysis : start_time", output_start_time(i)
             write(iulog,*)"  analysis : end_time  ", output_end_time(i)
             write(iulog,*)"  analysis : frequency ", output_frequency(i)
             select case (i)
             case (1)
                write(*,'(10(a,2x))')"   analysis : variables ",(trim(output_varnames1(j)),j=1,max_output_variables)
             case (2)
                write(*,'(10(a,2x))')"   analysis : variables ",(trim(output_varnames2(j)),j=1,max_output_variables)
             case (3)
                write(*,'(10(a,2x))')"   analysis : variables ",(trim(output_varnames3(j)),j=1,max_output_variables)
             case (4)
                write(*,'(10(a,2x))')"   analysis : variables ",(trim(output_varnames4(j)),j=1,max_output_variables)
             case (5)
                write(*,'(10(a,2x))')"   analysis : variables ",(trim(output_varnames5(j)),j=1,max_output_variables)
             end select
          end if
       end do
#endif

#if (defined MODEL_THETA_L && defined ARKODE)
       write(iulog,*)""
       write(iulog,*)"arkode: rel_tol = ",rel_tol
       write(iulog,*)"arkode: abs_tol = ",abs_tol
       write(iulog,*)"arkode: calc_nonlinear_stats = ",calc_nonlinear_stats
       write(iulog,*)"arkode: use_column_solver = ",use_column_solver
#endif

       ! display physical constants for HOMME stand alone simulations
       write(iulog,*)""
       write(iulog,*)"physconst: omega  = ",omega
       write(iulog,*)"physconst: rearth = ",rearth
       write(iulog,*)"physconst: rrearth= ",rrearth
       write(iulog,*)""
#endif

#if !defined(CAM) && !defined(SCREAM)
       write(iulog,*)" analysis interpolation = ", interpolate_analysis

       if(any(interpolate_analysis)) then
          write(iulog,*)" analysis interp nlat = ",interp_nlat
          write(iulog,*)" analysis interp nlon = ",interp_nlon
          write(iulog,*)" analysis interp lon0 = ",interp_lon0
          write(iulog,*)" analysis interp gridtype = ",interp_gridtype
          write(iulog,*)" analysis interpolation type = ",interp_type
       end if
#endif
! ^ ifndef CAM

#if !defined(CAM) && !defined(SCREAM)
      write(iulog,*)"HOMME SHA = ", HOMME_SHA1
#endif

      call print_clear_message()

!=======================================================================================================!
    endif

  end subroutine readnl

  subroutine print_clear_message()
    ! Those familiar with Homme can deduce dycore and transport type from
    ! options. But we should provide friendly message for others.
#if defined MODEL_THETA_L
    write(iulog,*) 'Running dycore: theta-l'
#elif defined _PRIM
    write(iulog,*) 'Running dycore: preqx'
#endif
    if (transport_alg == 0) then
       write(iulog,*) 'Running tracer transport method: horizonatally Eulerian'
    else
       write(iulog,*) 'Running tracer transport method: semi-Lagrangian'
    end if
  end subroutine print_clear_message

end module namelist_mod
