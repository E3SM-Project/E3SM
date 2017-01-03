#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module namelist_mod

  use kinds,      only: real_kind, iulog
  use params_mod, only: recursive, sfcurve
  use cube_mod,   only: rotate_grid
  use physical_constants, only: rearth, rrearth, omega

  use control_mod, only : &
    MAX_STRING_LEN,&
    MAX_FILE_LEN,  &
    partmethod,    &       ! Mesh partitioning method (METIS)
    topology,      &       ! Mesh topology
    test_case,     &       ! test case
    uselapi,       &
    multilevel,    &
    numnodes,      &
    sub_case,      &
    tasknum,       &       ! used dg model in AIX machine
    remapfreq,     &       ! number of steps per remapping call
    remap_type,    &       ! selected remapping option
    statefreq,     &       ! number of steps per printstate call
    restartfreq,   &
    restartfile,   &       ! name of the restart file for INPUT
    restartdir,    &       ! name of the restart directory for OUTPUT
    runtype,       &
    integration,   &       ! integration method
    theta_hydrostatic_mode,       &   
    use_semi_lagrange_transport , &   ! conservation or non-conservation formulaton
    use_semi_lagrange_transport_local_conservation , &   ! local conservation vs. global 
    tstep_type,    &
    cubed_sphere_map, &
    qsplit,        &
    rsplit,        &
    rk_stage_user, &
    LFTfreq,       &
    prescribed_wind, &
    ftype,         &
    energy_fixer,  &
    limiter_option,&
    fine_ne,       &
    max_hypervis_courant, &
    nu,            &
    nu_s,          &
    nu_q,          &
    nu_div,        &
    nu_p,          &
    nu_top,        &
    hypervis_scaling,   &  ! use tensor HV instead of scalar coefficient
    disable_diagnostics, & ! use to disable diagnostics for timing reasons
    psurf_vis,    &
    hypervis_order,       &
    hypervis_power,       &
    hypervis_subcycle,    &
    hypervis_subcycle_q,  &
    smooth_phis_numcycle, &
    smooth_sgh_numcycle,  &
    smooth_phis_nudt,     &
    initial_total_mass,   & ! set > 0 to set the initial_total_mass
    u_perturb,     &        ! J&W baroclinic test perturbation size
    moisture,      &
    vform,         &
    vfile_mid,     &
    vfile_int,     &
    vanalytic,     &
    vtop,          &
    precon_method, &
    maxits,        &
    tol,           &
    debug_level,   &
    vert_remap_q_alg

#ifndef CAM
  use control_mod, only:              &
    pertlim,                          &
    dcmip2_0_h0,                      &
    dcmip2_0_rm,                      &
    dcmip2_0_zetam,                   &
    dcmip2_x_ueq,                     &
    dcmip2_x_h0,                      &
    dcmip2_x_d,                       &
    dcmip2_x_xi
#endif

  use thread_mod,     only: nthreads, nthreads_accel, omp_set_num_threads, omp_get_max_threads, vert_num_threads, vthreads
  use dimensions_mod, only: ne, np, nnodes, nmpi_per_node, npart, qsize, qsize_d, set_mesh_dimensions
#ifdef CAM
  use time_mod,       only: nsplit, smooth, phys_tscale
#else
  use time_mod,       only: tstep, ndays,nmax, nendstep,secpday, smooth, secphr, nsplit, phys_tscale
#endif
  use parallel_mod,   only: parallel_t,  iam, abortmp, &
       partitionfornodes, useframes, mpireal_t, mpilogical_t, mpiinteger_t, mpichar_t
  use cg_mod,         only: cg_no_debug

#ifndef CAM
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

#ifdef CAM
  subroutine readnl(par, NLFileName)
    use units, only : getunit, freeunit
    use mesh_mod, only : MeshOpen
    character(len=*), intent(in) :: NLFilename  ! namelist filename
#else
  subroutine readnl(par)
    use mesh_mod, only : MeshOpen
#endif
    type (parallel_t), intent(in) ::  par
    character(len=MAX_FILE_LEN) :: mesh_file
    integer :: se_ftype, se_limiter_option
    integer :: se_phys_tscale, se_nsplit
    integer :: interp_nlat, interp_nlon, interp_gridtype, interp_type
    integer :: i, ii, j
    integer  :: ierr
    character(len=80) :: errstr, arg
    real(kind=real_kind) :: dt_max
#ifdef CAM
    character(len=MAX_STRING_LEN) :: se_topology
    integer :: se_partmethod
    integer :: se_ne
    integer :: unitn
    character(len=*), parameter ::  subname = "homme:namelist_mod"
    ! These items are only here to keep readnl from crashing. Remove when possible
    integer :: se_fv_nphys
    character(len=80)  :: se_write_phys_grid
    character(len=256) :: se_phys_grid_file
#endif
    ! ============================================
    ! Namelists
    ! ============================================

    namelist /ctl_nl/ PARTMETHOD,       &         ! mesh partitioning method
                      TOPOLOGY,         &         ! mesh topology
#ifdef CAM
      se_partmethod,     &
      se_topology,       &
      se_ne,             &
      se_limiter_option, &
      vthreads,          &         ! number of vertical/column threads per horizontal thread
#else
      qsize,             &         ! number of SE tracers
      nthreads,          &         ! number of threads per process
      vert_num_threads,  &         ! number of threads per process
      nthreads_accel,    &         ! number of threads per an accelerator process
      limiter_option,    &
      smooth,            &         ! timestep Filter
      omega,             &
      pertlim,           &         ! temperature initial perturbation
      omega,                   &   ! scaled rotation rate
      rearth,                  &   ! scaled earth radius
#endif
      npart,         &
      uselapi,       &
      multilevel,    &
      useframes,     &
      numnodes,      &
      ne,            &             ! element resolution factor
      tasknum,       &
      remapfreq,     &             ! number of steps per remapping call
      remap_type,    &             ! selected remapping option
      statefreq,     &             ! number of steps per printstate call
      integration,   &             ! integration method
      theta_hydrostatic_mode,       &   
      use_semi_lagrange_transport , &
      use_semi_lagrange_transport_local_conservation , &
      tstep_type,    &
      cubed_sphere_map, &
      qsplit,        &
      rsplit,        &
      rk_stage_user, &
      LFTfreq,       &
      disable_diagnostics, &
      prescribed_wind, &
      se_ftype,        &           ! forcing type
      energy_fixer,    &
      fine_ne,         &
      max_hypervis_courant, &
      nu,            &
      nu_s,          &
      nu_q,          &
      nu_div,        &
      nu_p,          &
      nu_top,        &
      psurf_vis,     &
      hypervis_order,    &
      hypervis_power,    &
      hypervis_subcycle, &
      hypervis_subcycle_q, &
      hypervis_scaling, &
      smooth_phis_numcycle, &
      smooth_sgh_numcycle, &
      smooth_phis_nudt, &
      initial_total_mass, &
      u_perturb,     &
      rotate_grid,   &
      mesh_file,     &               ! Name of mesh file
      vert_remap_q_alg


#ifdef CAM
    namelist  /ctl_nl/ SE_NSPLIT,  &                ! number of dynamics steps per physics timestep
      se_phys_tscale
#else
    namelist /ctl_nl/test_case,       &             ! test case idenitfier
      sub_case,        &             ! generic test case parameter
      nmax,            &             ! number of steps
      ndays,           &             ! number of days to simulate
      restartfreq,     &
      restartfile,     &             ! name of the restart file for INPUT
      restartdir,      &             ! name of the restart directory for OUTPUT
      runtype,         &
      tstep,           &             ! tracer time step
      moisture
    ! control parameters for dcmip stand-alone tests
    namelist /ctl_nl/     &
      dcmip2_0_h0,        & !dcmip2-0 mountain height           (meters)
      dcmip2_0_rm,        & !dcmip2-0 mountain range radius     (radians)
      dcmip2_0_zetam,     & !dcmip2-0 mountain range half width (radians)
      dcmip2_x_ueq,       & !dcmip2-x windspeed at equator      (m/s)
      dcmip2_x_h0,        & !dcmip2-x mountain height           (m)
      dcmip2_x_d,         & !dcmip2-x mountain half-width       (m)
      dcmip2_x_xi           !dcmip2-x mountain wavelength       (m)
    namelist /vert_nl/        &
      vform,              &
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
      interp_gridtype,     &
      interp_type,         &
      interpolate_analysis
#endif
! ^ ifndef CAM


    namelist /solver_nl/precon_method, &
      maxits,        &
      tol,           &
      debug_level


    ! ==========================
    ! Set the default partmethod
    ! ==========================

    PARTMETHOD    = RECURSIVE
    npart         = 1
    useframes     = 0
    multilevel    = 1
    uselapi       = .TRUE.
#ifdef CAM
    ! set all CAM defaults
    ! CAM requires forward-in-time, subcycled dynamics
    ! RK2 3 stage tracers, sign-preserving conservative
    tstep_type              = 1      ! forward-in-time RK methods
    qsplit=4; rk_stage_user=3
    se_limiter_option=4
    se_ftype = 2
    se_partmethod = -1
    se_ne       = -1
    se_topology = 'none'
    se_phys_tscale=0
    se_nsplit = 1
    qsize = qsize_d
    vthreads = 1
#else
    ndays         = 0
    nmax          = 12
    nthreads = 1
    vert_num_threads = 1
    nthreads_accel = -1
    se_ftype = ftype   ! MNL: For non-CAM runs, ftype=0 in control_mod
    phys_tscale=0
    nsplit = 1
    pertlim = 0.0_real_kind
#endif
    sub_case      = 1
    numnodes      = -1
    restartfreq   = -100
    restartdir    = "./restart/"
    runtype       = 0
    statefreq     = 1
    remapfreq     = 240
    remap_type    = "parabolic"
    tasknum       =-1
    integration   = "explicit"
    moisture      = "dry"
    nu_top=0
    initial_total_mass=0
    mesh_file='none'
    ne              = 0
    use_semi_lagrange_transport   = .false.
    use_semi_lagrange_transport_local_conservation   = .false.
    disable_diagnostics = .false.


    ! =======================
    ! Read namelist variables
    ! =======================

!   write(iulog,*) "masterproc=",par%masterproc

    if (par%masterproc) then

       write(iulog,*)"reading ctl namelist..."
#if defined(CAM)
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
      !=================================
      ! override the selected partition
      ! method and set it to SFCURVE
      !=================================
      PARTMETHOD = SFCURVE
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
#ifndef CAM
       if (tstep <= 0) then
          call abortmp('tstep must be > 0')
       end if
       if (ndays .gt. 0) then
          nmax = ndays * (secpday/tstep)
          restartfreq  = restartfreq*(secpday/tstep)
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

#ifndef CAM

#ifdef _PRIM
       write(iulog,*) "reading physics namelist..."
       if (test_case      == "held_suarez0"   .or. &
           test_case(1:10)== "baroclinic"     .or. &
           test_case(1:13)== "jw_baroclinic"  .or. &
           test_case(1:5) == "dcmip"          .or. &
           test_case(1:4) == "asp_")  then
         write(iulog,*) "reading vertical namelist..."

#if defined(OSF1) || defined(_NAMELIST_FROM_FILE)
         read(unit=7,nml=vert_nl)
#else
         read(*,nml=vert_nl)
#endif
         vform      = trim(adjustl(vform))
         vfile_mid  = trim(adjustl(vfile_mid))
         vfile_int  = trim(adjustl(vfile_int))

         write(iulog,*) '  vform =',trim(vform)
         write(iulog,*) '  vfile_mid=',trim(vfile_mid)
         write(iulog,*) '  vfile_int=',trim(vfile_int)
         write(iulog,*) '  vanalytic=',vanalytic
         if(vanalytic==1) then
         write(iulog,*) '  vtop=',vtop
         endif

       end if
#endif

!      Default interpolation grid  (0 = auto compute based on ne,nv)  interpolation is off by default
#ifdef PIO_INTERP
       interpolate_analysis=.true.
#else
       interpolate_analysis=.false.
#endif
       interp_nlat =  0
       interp_nlon = 0
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
             output_end_time(i)  = output_end_time(i)*(secpday/tstep)
          else if(output_timeunits(i).eq.2) then  ! per_hour
             output_frequency(i) = output_frequency(i)*(secphr/tstep)
             output_start_time(i)= output_start_time(i)*(secphr/tstep)
             output_end_time(i)  = output_end_time(i)*(secphr/tstep)
          end if
          if(output_end_time(i)<0) then
             output_end_time(i)=nEndStep
          endif
       end do

!=======================================================================================================!

#if defined(OSF1) || defined(_NAMELIST_FROM_FILE)
       close(unit=7)
#endif
#endif
! ^ ifndef CAM
    end if

#ifdef CAM
    if(se_partmethod /= -1) partmethod = se_partmethod
    if(se_ne /= -1) ne = se_ne
    if(se_topology .ne. 'none') topology = se_topology
#endif

    call MPI_barrier(par%comm,ierr)

    npart  = par%nprocs

    ! Broadcast namelist variables to all MPI processes

    call MPI_bcast(PARTMETHOD ,     1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(TOPOLOGY,        MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
    call MPI_bcast(test_case,       MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
    call MPI_bcast(tasknum,         1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(ne,              1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(qsize,           1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(sub_case,        1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(remapfreq,       1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(remap_type, MAX_STRING_LEN, MPIChar_t, par%root, par%comm, ierr)
    call MPI_bcast(statefreq,       1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(restartfreq,     1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(multilevel,      1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(useframes,       1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(runtype,         1,MPIinteger_t,par%root,par%comm,ierr)

#ifdef CAM
    phys_tscale     = se_phys_tscale
    limiter_option  = se_limiter_option
    nsplit          = se_nsplit
    call MPI_bcast(vthreads  ,      1, MPIinteger_t, par%root,par%comm,ierr)
#else
    call MPI_bcast(omega,           1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(pertlim,         1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(tstep,           1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(nmax,            1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(NTHREADS,        1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(vert_num_threads,1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(nthreads_accel,  1, MPIinteger_t, par%root,par%comm,ierr)
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
#endif
    call MPI_bcast(smooth,          1, MPIreal_t,    par%root,par%comm,ierr)
    call MPI_bcast(phys_tscale,     1, MPIreal_t,    par%root,par%comm,ierr)
    call MPI_bcast(NSPLIT,          1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(limiter_option,  1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(se_ftype,        1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(vert_remap_q_alg,1, MPIinteger_t, par%root,par%comm,ierr)

    call MPI_bcast(fine_ne,         1, MPIinteger_t, par%root,par%comm,ierr)
    call MPI_bcast(max_hypervis_courant,1,MPIreal_t, par%root,par%comm,ierr)
    call MPI_bcast(nu,              1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(nu_s,            1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(nu_q,            1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(nu_div,          1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(nu_p,            1, MPIreal_t   , par%root,par%comm,ierr)
    call MPI_bcast(nu_top,          1, MPIreal_t   , par%root,par%comm,ierr)

    call MPI_bcast(disable_diagnostics,1,MPIlogical_t,par%root,par%comm,ierr)
    call MPI_bcast(psurf_vis,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_order,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_power,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_scaling,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_subcycle,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_subcycle_q,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(smooth_phis_numcycle,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(smooth_sgh_numcycle,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(smooth_phis_nudt,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(initial_total_mass ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(u_perturb     ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(rotate_grid   ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(integration,MAX_STRING_LEN,MPIChar_t ,par%root,par%comm,ierr)
    call MPI_bcast(mesh_file,MAX_FILE_LEN,MPIChar_t ,par%root,par%comm,ierr)
    call MPI_bcast(theta_hydrostatic_mode ,1,MPIlogical_t,par%root,par%comm,ierr)
    call MPI_bcast(use_semi_lagrange_transport ,1,MPIlogical_t,par%root,par%comm,ierr)
    call MPI_bcast(use_semi_lagrange_transport_local_conservation ,1,MPIlogical_t,par%root,par%comm,ierr)
    call MPI_bcast(tstep_type,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(cubed_sphere_map,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(qsplit,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(rsplit,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(rk_stage_user,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(LFTfreq,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(prescribed_wind,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(moisture,MAX_STRING_LEN,MPIChar_t ,par%root,par%comm,ierr)

    call MPI_bcast(restartfile,MAX_STRING_LEN,MPIChar_t ,par%root,par%comm,ierr)
    call MPI_bcast(restartdir,MAX_STRING_LEN,MPIChar_t ,par%root,par%comm,ierr)

    call MPI_bcast(uselapi,1,MPIlogical_t,par%root,par%comm,ierr)

    if (integration == "full_imp") then
       call MPI_bcast(precon_method,MAX_STRING_LEN,MPIChar_t,par%root,par%comm,ierr)
       call MPI_bcast(maxits     ,1,MPIinteger_t,par%root,par%comm,ierr)
       call MPI_bcast(tol        ,1,MPIreal_t   ,par%root,par%comm,ierr)
    end if

    call MPI_bcast(vform    , MAX_STRING_LEN, MPIChar_t   , par%root, par%comm,ierr)
    call MPI_bcast(vfile_mid, MAX_STRING_LEN, MPIChar_t   , par%root, par%comm,ierr)
    call MPI_bcast(vfile_int, MAX_STRING_LEN, MPIChar_t   , par%root, par%comm,ierr)
    call MPI_bcast(vanalytic, 1,              MPIinteger_t, par%root, par%comm,ierr)
    call MPI_bcast(vtop     , 1,              MPIreal_t   , par%root, par%comm,ierr)

#ifndef CAM
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

#ifdef IS_ACCELERATOR
    if (nthreads_accel > 0) then
        nthreads = nthreads_accel
    end if
#endif

    ! sanity check on thread count
    ! HOMME will run if if nthreads > max, but gptl will print out GB of warnings.
    if (NThreads*vert_num_threads > omp_get_max_threads()) then
       if(par%masterproc) write(iulog,*) "Main:NThreads=",NThreads
       if(par%masterproc) print *,'omp_get_max_threads() = ',OMP_get_max_threads()
       if(par%masterproc) print *,'requested threads exceeds OMP_get_max_threads()'
       call abortmp('stopping')
    endif
    call omp_set_num_threads(NThreads*vert_num_threads)


    ! if user sets hypervis_subcycle=-1, then use automatic formula
    if (hypervis_subcycle==-1) then
       if (np==4) then
          ! 1.25d23 worked out by testing, for nv==4
          ! a little confusing:
          ! u,v:  nu and hypervis_subcycle
          ! T:    nu_s and hypervis_subcycle
          ! Q:    nu and hypervis_subcycle_q
          dt_max = 1.25d23/(nu*ne**4.0)
          hypervis_subcycle_q = ceiling( tstep/dt_max )
          hypervis_subcycle   = ceiling( tstep/dt_max )
       else
          call abortmp('hypervis_subcycle auto determine only supported for nv==4')
       endif
    endif
#endif

    if (ne /=0) then
    if (mesh_file /= "none" .and. mesh_file /= "/dev/null") then
      write (*,*) "namelist_mod: mesh_file:",trim(mesh_file), &
                  " and ne:",ne," are both sepcified in the input file."
      write (*,*) "Specify one or the other, but not both."
      call abortmp("Do not specify ne if using a mesh file input.")
    end if
    end if
    if (par%masterproc) write (iulog,*) "Mesh File:", trim(mesh_file)
    if (ne.eq.0) then
       if (par%masterproc) write (iulog,*) "Opening Mesh File:", trim(mesh_file)
      call set_mesh_dimensions()
      call MeshOpen(mesh_file, par)
    end if
    ! set map
    if (cubed_sphere_map<0) then
       cubed_sphere_map=0  ! default is equi-angle gnomonic
       if (ne.eq.0) cubed_sphere_map=2  ! element_local for var-res grids
    endif
    if (par%masterproc) write (iulog,*) "Reference element projection: cubed_sphere_map=",cubed_sphere_map

!logic around different hyperviscosity options
    if (hypervis_power /= 0) then
      if (hypervis_scaling /= 0) then
        print *,'Both hypervis_power and hypervis_scaling are nonzero.'
        print *,'(1) Set hypervis_power=1, hypervis_scaling=0 for HV based on an element area.'
        print *,'(2) Set hypervis_power=0 and hypervis_scaling=1 for HV based on a tensor.'
        print *,'(3) Set hypervis_power=0 and hypervis_scaling=0 for constant HV.'
          call abortmp("Error: hypervis_power>0 and hypervis_scaling>0")
      endif
    endif


    ftype = se_ftype

#ifdef _PRIM
    rk_stage_user=3  ! 3d PRIM code only supports 3 stage RK tracer advection
    ! CHECK phys timescale, requires se_ftype=0 (pure tendencies for forcing)
    if (phys_tscale/=0) then
       if (ftype>0) call abortmp('user specified se_phys_tscale requires se_ftype<=0')
    endif
    if (limiter_option==8 .or. limiter_option==84) then
       if (hypervis_subcycle_q/=1) then
          call abortmp('limiter 8,84 requires hypervis_subcycle_q=1')
       endif
    endif
#endif

    if((prescribed_wind/=0).and.(prescribed_wind/=1))then
          call abortmp('prescribed_wind should be either 0 or 1')
    endif

    
    if (use_semi_lagrange_transport .and. rsplit == 0) then
       call abortmp('The semi-Lagrange Transport option requires 0 < rsplit')
    end if

!=======================================================================================================!
#ifdef CAM
    nmpi_per_node=1
#endif
#ifndef CAM
    call MPI_bcast(interpolate_analysis, 7,MPIlogical_t,par%root,par%comm,ierr)
    call MPI_bcast(interp_nlat , 1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(interp_nlon , 1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(interp_gridtype , 1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(interp_type , 1,MPIinteger_t,par%root,par%comm,ierr)

    call MPI_bcast(replace_vec_by_vordiv ,MAX_VECVARS ,MPIlogical_t,par%root,par%comm,ierr)
    call MPI_bcast(vector_uvars ,10*MAX_VECVARS ,MPIChar_t,par%root,par%comm,ierr)
    call MPI_bcast(vector_vvars ,10*MAX_VECVARS ,MPIChar_t,par%root,par%comm,ierr)

    call set_interp_parameter('gridtype',interp_gridtype)
    call set_interp_parameter("itype",interp_type)
    if(any(interpolate_analysis)) then
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
    if(nu_s<0) nu_s=nu
    if(nu_q<0) nu_q=nu
    if(nu_div<0) nu_div=nu


    if (multilevel <= 0) then
      nmpi_per_node = 1
    endif

    nnodes = npart/nmpi_per_node

    if(numnodes > 0 .and. multilevel .eq. 1) then
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
#ifndef CAM
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
       write(iulog,*)"readnl: vert_num_threads = ",vert_num_threads
       write(iulog,*)"readnl: nthreads_accel = ",nthreads_accel
#endif

       write(iulog,*)"readnl: ne,np         = ",NE,np
       write(iulog,*)"readnl: partmethod    = ",PARTMETHOD
       write(iulog,*)'readnl: nmpi_per_node = ',nmpi_per_node
       write(iulog,*)"readnl: vthreads      = ",vthreads
       write(iulog,*)'readnl: multilevel    = ',multilevel
       write(iulog,*)'readnl: useframes     = ',useframes
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
       write(iulog,*)"readnl: use_semi_lagrange_transport   = ",use_semi_lagrange_transport
       write(iulog,*)"readnl: use_semi_lagrange_transport_local_conservation=",use_semi_lagrange_transport_local_conservation
       write(iulog,*)"readnl: tstep_type    = ",tstep_type
       write(iulog,*)"readnl: vert_remap_q_alg  = ",vert_remap_q_alg
#ifdef CAM
       write(iulog,*)"readnl: se_nsplit         = ", NSPLIT
       write(iulog,*)"readnl: se_ftype          = ",ftype
       write(iulog,*)"readnl: se_limiter_option = ",limiter_option
#else
       write(iulog,*)"readnl: tstep          = ",tstep
       write(iulog,*)"readnl: ftype          = ",ftype
       write(iulog,*)"readnl: limiter_option = ",limiter_option
#endif
       write(iulog,*)"readnl: qsplit        = ",qsplit
       write(iulog,*)"readnl: vertical remap frequency rsplit (0=disabled): ",rsplit

       write(iulog,*)"readnl: runtype       = ",runtype

       if (hypervis_power /= 0)then
          write(iulog,*)"Variable scalar hyperviscosity: hypervis_power=",hypervis_power
          write(iulog,*)"max_hypervis_courant = ", max_hypervis_courant
          write(iulog,*)"Equivalent ne in fine region = ", fine_ne
       elseif(hypervis_scaling /=0)then
          write(iulog,*)"Tensor hyperviscosity:  hypervis_scaling=",hypervis_scaling
       else
          write(iulog,*)"Constant (hyper)viscosity used."
       endif

       write(iulog,*)"hypervis_subcycle, hypervis_subcycle_q = ",&
            hypervis_subcycle,hypervis_subcycle_q
       !write(iulog,*)"psurf_vis: ",psurf_vis
       write(iulog,'(a,2e9.2)')"viscosity:  nu (vor/div) = ",nu,nu_div
       write(iulog,'(a,2e9.2)')"viscosity:  nu_s      = ",nu_s
       write(iulog,'(a,2e9.2)')"viscosity:  nu_q      = ",nu_q
       write(iulog,'(a,2e9.2)')"viscosity:  nu_p      = ",nu_p
       write(iulog,'(a,2e9.2)')"viscosity:  nu_top      = ",nu_top
       write(iulog,*)"PHIS smoothing:  ",smooth_phis_numcycle,smooth_phis_nudt
       write(iulog,*)"SGH  smoothing:  ",smooth_sgh_numcycle

       if(initial_total_mass>0) then
          write(iulog,*) "initial_total_mass = ",initial_total_mass
       end if

#ifndef CAM
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

       ! display physical constants for HOMME stand alone simulations
       write(iulog,*)""
       write(iulog,*)"physconst: omega  = ",omega
       write(iulog,*)"physconst: rearth = ",rearth
       write(iulog,*)"physconst: rrearth= ",rrearth
       write(iulog,*)""
#endif

#ifndef CAM
       write(iulog,*)" analysis interpolation = ", interpolate_analysis

       if(any(interpolate_analysis)) then
          write(iulog,*)" analysis interp nlat = ",interp_nlat
          write(iulog,*)" analysis interp nlon = ",interp_nlon
          write(iulog,*)" analysis interp gridtype = ",interp_gridtype
          write(iulog,*)" analysis interpolation type = ",interp_type
       end if
#endif
! ^ ifndef CAM

!=======================================================================================================!
    endif

  end subroutine readnl

end module namelist_mod
