module homme_params_mod
  use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_double, c_bool, C_NULL_CHAR

  !**********************************
  ! TODO: there's a lot of params that can be removed for scream,
  !       but that requires using namelists different from those
  !       stored in homme. Until then, add a lot of stuff to
  !       the ctl_nl namelist declaration here, which we will
  !       just ignore. When you start using scream namelist files,
  !       you should revisit this file, and purge all the stuff
  !       that is no longer needed
  !***********************************
  use kinds,        only: real_kind, iulog
  use parallel_mod, only: abortmp

  use physical_constants, only: &
    rearth,                     &
    rrearth,                    &
    omega

  use dimensions_mod, only: &
    ne,                     &
    np,                     &
    npart,                  &
    qsize,                  &   ! number of SE tracers
    qsize_d                     ! Compile time upper bound for qsize

  use cube_mod, only: rotate_grid
  use time_mod, only: tstep, nEndStep, secpday, ndays, nmax, &
                      smooth ! Unused in SCREAM
  use thread_mod,   only: nthreads, vthreads      ! Unused in SCREAM
  use parallel_mod, only: mpireal_t, mpilogical_t, mpiinteger_t, mpichar_t

  use control_mod, only:    &
    MAX_STRING_LEN,         &
    MAX_FILE_LEN,           &
    cubed_sphere_map,       &
    partmethod,             &    ! Mesh partitioning method (METIS)
    coord_transform_method, &
    vert_remap_q_alg,       &
    theta_advect_form,      &
    theta_hydrostatic_mode, &   
    topology,               &    ! Mesh topology
    tstep_type,             &
    ftype,                  & ! Unused in SCREAM
    dt_tracer_factor,       &
    dt_remap_factor,        &
    qsplit,                 &
    rsplit,                 &
    limiter_option,         &
    moisture,               &
    use_moisture,           &
    runtype,                & ! Unused in SCREAM
    max_hypervis_courant,   &
    nu,                     &
    nu_s,                   &
    nu_q,                   &
    nu_div,                 &
    nu_p,                   &
    nu_top,                 &
    hypervis_order,         &
    hypervis_scaling,       &
    hypervis_power,         &
    hypervis_subcycle,      &
    hypervis_subcycle_tom,  &
    hypervis_subcycle_q,    &
    transport_alg,          &
    disable_diagnostics,    &
    test_case,              &
    dcmip16_mu,             &
    dcmip16_mu_s,           &
    dcmip16_mu_q,           &
    dcmip16_prec_type,      &
    dcmip16_pbl_type,       &
    integration,            & ! Unused in SCREAM
    restartfreq,            & ! Unused in SCREAM
    restartfile,            & ! Unused in SCREAM
    u_perturb,              &
    se_fv_phys_remap_alg,   &
    statefreq,              & ! Unused in SCREAM
    vform,                  &
    vfile_mid,              &
    vfile_int,              &
    timestep_make_subcycle_parameters_consistent

  public :: init_params_f90
  public :: get_homme_int_param_f90
  public :: get_homme_real_param_f90
  public :: get_homme_bool_param_f90

  public :: set_homme_int_param_f90
  public :: set_homme_real_param_f90
  public :: set_homme_bool_param_f90
contains

  subroutine init_params_f90 (nl_fname_c) bind(c)
    use iso_c_binding,     only: c_ptr, c_f_pointer
    use shr_file_mod,      only: getunit=>shr_file_getUnit, freeunit=>shr_file_freeUnit
    use homme_context_mod, only: is_parallel_inited, is_params_inited, &
                                 par
    use params_mod, only: SFCURVE, SPHERE_COORDS

    ! Inputs
    type (c_ptr), intent(in) :: nl_fname_c

    ! Locals
    character(len=256), pointer :: nl_fname
    integer :: ierr, unitn, se_ftype
    real(kind=real_kind) :: dt_max
    character(len=MAX_FILE_LEN) :: mesh_file ! Unused in SCREAM

    !-----------------------------!
    !   Declare valid namelists   !
    !-----------------------------!

    namelist /ctl_nl/test_case, &
      u_perturb,                &
      partmethod,               &         ! mesh partitioning method
      coord_transform_method,   &
      vert_remap_q_alg,         &
      theta_advect_form,        &
      theta_hydrostatic_mode,   &   
      ne,                       &
      ndays,                    &
      nmax,                     &
      rotate_grid,              &
      topology,                 &         ! Mesh topology
      tstep_type,               &
      tstep,                    &
      qsplit,                   &
      qsize,                    &         ! number of SE tracers
      rsplit,                   &
      omega,                    &
      rearth,                   &
      limiter_option,           &
      nu,                       &
      nu_s,                     &
      nu_q,                     &
      nu_div,                   &
      nu_p,                     &
      nu_top,                   &
      hypervis_order,           &
      hypervis_scaling,         &
      hypervis_subcycle,        &
      hypervis_subcycle_tom,    &
      moisture,                 & ! Unused in SCREAM
      vthreads,                 & ! Unused in SCREAM
      nthreads,                 & ! Unused in SCREAM
      statefreq,                & ! Unused in SCREAM
      restartfreq,              & ! Unused in SCREAM
      restartfile,              & ! Unused in SCREAM
      runtype,                  & ! Unused in SCREAM
      mesh_file,                & ! Unused in SCREAM
      integration,              & ! Unused in SCREAM
      smooth,                   & ! Unused in SCREAM
      se_ftype                    ! Unused in SCREAM

    namelist /vert_nl/    &
      vform,              &
      vfile_mid,          &
      vfile_int

    ! Sanity check (we need to bcast options to all ranks)
    if (.not. is_parallel_inited) then
      call abortmp ("Error! 'init_parallel_f90' was not called yet.\n")
    endif

    !-----------------------------!
    !        Set defaults         !
    !-----------------------------!

    tstep = -1
    ndays =  0
    nmax  = 12
    moisture = 'dry'
    partmethod = SFCURVE
    coord_transform_method = SPHERE_COORDS
    transport_alg = 0
    runtype = 0
    statefreq = 99999

    !-----------------------------!
    !     Parse namelist file     !
    !-----------------------------!

    ! Open namelist file
    unitn = getunit()
    call c_f_pointer(nl_fname_c,nl_fname)
    open( unitn, file=trim(nl_fname), status='old' )

    ! Parse all namelist sections
    read (unit=unitn,nml=ctl_nl,iostat=ierr)
    if (ierr < 0) then
      write(6,*) 'ierr =',ierr
      call abortmp( '[scream] init_params_f90 :: namelist read returns an'// &
             ' end of file or end of record condition' )
    end if
    read(unit=unitn,nml=vert_nl)
    vform      = trim(adjustl(vform))
    vfile_mid  = trim(adjustl(vfile_mid))
    vfile_int  = trim(adjustl(vfile_int))

    ! Close namelist file
    close( unitn )
    call freeunit( unitn )

    !-----------------------------!
    !      Broadcast options      !
    !-----------------------------!

    ! Geometry params
    call MPI_bcast(partmethod,      1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(ne,              1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(cubed_sphere_map,1, MPIinteger_t ,par%root, par%comm, ierr)
    call MPI_bcast(topology, MAX_STRING_LEN, MPIChar_t, par%root, par%comm, ierr)

    ! Time-stepping params
    call MPI_bcast(nmax,       1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(ndays,      1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(tstep_type, 1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(qsplit,     1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(rsplit,     1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(tstep,      1, MPIreal_t,    par%root, par%comm, ierr)

    ! Algorithmic params
    call MPI_bcast(limiter_option,         1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(theta_advect_form,      1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(vert_remap_q_alg,       1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(theta_hydrostatic_mode, 1, MPIlogical_t, par%root, par%comm, ierr)
    call MPI_bcast(transport_alg ,         1, MPIinteger_t, par%root, par%comm, ierr)

    ! Physical params
    call MPI_bcast(omega,  1, MPIreal_t, par%root, par%comm, ierr)
    call MPI_bcast(rearth, 1, MPIreal_t, par%root, par%comm, ierr)

    ! Hyperviscosity params
    call MPI_bcast(nu,     1, MPIreal_t, par%root, par%comm, ierr)
    call MPI_bcast(nu_s,   1, MPIreal_t, par%root, par%comm, ierr)
    call MPI_bcast(nu_q,   1, MPIreal_t, par%root, par%comm, ierr)
    call MPI_bcast(nu_div, 1, MPIreal_t, par%root, par%comm, ierr)
    call MPI_bcast(nu_p,   1, MPIreal_t, par%root, par%comm, ierr)
    call MPI_bcast(nu_top, 1, MPIreal_t, par%root, par%comm, ierr)

    call MPI_bcast(hypervis_order,        1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(hypervis_power,        1, MPIreal_t,    par%root, par%comm, ierr)
    call MPI_bcast(hypervis_scaling,      1, MPIreal_t,    par%root, par%comm, ierr)
    call MPI_bcast(hypervis_subcycle,     1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(hypervis_subcycle_tom, 1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(hypervis_subcycle_q,   1, MPIinteger_t, par%root, par%comm, ierr)

    ! Case config params
    call MPI_bcast(disable_diagnostics, 1, MPIlogical_t, par%root, par%comm, ierr)
    call MPI_bcast(qsize,               1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(statefreq,           1, MPIinteger_t, par%root, par%comm, ierr)

    call MPI_bcast(test_case, MAX_STRING_LEN, MPIChar_t, par%root, par%comm, ierr)
    call MPI_bcast(moisture,  MAX_STRING_LEN, MPIChar_t, par%root, par%comm, ierr)

    call MPI_bcast(dcmip16_mu,        1, MPIreal_t   , par%root, par%comm, ierr)
    call MPI_bcast(dcmip16_mu_s,      1, MPIreal_t   , par%root, par%comm, ierr)
    call MPI_bcast(dcmip16_mu_q,      1, MPIreal_t   , par%root, par%comm, ierr)
    call MPI_bcast(dcmip16_prec_type, 1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(dcmip16_pbl_type , 1, MPIinteger_t, par%root, par%comm, ierr)

    ! Vertical coord params
    call MPI_bcast(vform,     MAX_STRING_LEN, MPIChar_t, par%root, par%comm, ierr)
    call MPI_bcast(vfile_mid, MAX_STRING_LEN, MPIChar_t, par%root, par%comm, ierr)
    call MPI_bcast(vfile_int, MAX_STRING_LEN, MPIChar_t, par%root, par%comm, ierr)

    !-----------------------------!
    !   Check/deduce parameters   !
    !-----------------------------!


    ierr = timestep_make_subcycle_parameters_consistent(par, rsplit, qsplit, dt_remap_factor, dt_tracer_factor)

    ftype = 0
    npart = par%nprocs

    use_moisture = ( moisture /= "dry") 

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

    if (ndays>0) then
      nmax = INT(ndays * (secpday/tstep))
    end if
    nEndStep = nmax

    ! set map
    if (cubed_sphere_map<0) then
       cubed_sphere_map=2  ! theta model default = element local
    endif

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

    if (limiter_option==8 .or. limiter_option==84 .or. limiter_option == 9) then
       if (hypervis_subcycle_q/=1 .and. transport_alg == 0) then
          call abortmp('limiter 8,84,9 require hypervis_subcycle_q=1')
       endif
    endif
    if (transport_alg == 0 .and. dt_remap_factor > 0 .and. dt_remap_factor < dt_tracer_factor) then
       call abortmp('Only SL transport supports vertical remap time step < tracer time step.')
    end if

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

    if (qsize>qsize_d) then
      print *, "qsize, qsize_d:", qsize, qsize_d
      call abortmp('user specified qsize > qsize_d parameter in dimensions_mod.F90')
    endif

    !-----------------------------!
    !       Print parameters      !
    !-----------------------------!

    if (par%masterproc) then
       write(iulog,*)"done reading namelist..."

       write(iulog,*)"readnl: topology      = ",TRIM( TOPOLOGY )

       write(iulog,*)"readnl: test_case     = ",TRIM(test_case)
       write(iulog,*)"readnl: omega         = ",omega
       write(iulog,*)"readnl: ndays         = ",ndays
       write(iulog,*)"readnl: nmax          = ",nmax

       write(iulog,*)"readnl: qsize,qsize_d = ",qsize,qsize_d

       write(iulog,*)"readnl: ne,np         = ",NE,np
       write(iulog,*)"readnl: partmethod    = ",PARTMETHOD
       write(iulog,*)"readnl: COORD_TRANSFORM_METHOD    = ",COORD_TRANSFORM_METHOD

       write(iulog,*)"readnl: theta_hydrostatic_mode = ",theta_hydrostatic_mode
       write(iulog,*)"readnl: transport_alg   = ",transport_alg
       write(iulog,*)"readnl: tstep_type    = ",tstep_type
       write(iulog,*)"readnl: theta_advect_form = ",theta_advect_form
       write(iulog,*)"readnl: vert_remap_q_alg  = ",vert_remap_q_alg
       write(iulog,*)"readnl: tstep          = ",tstep
       write(iulog,*)"readnl: ftype          = ",ftype
       write(iulog,*)"readnl: limiter_option = ",limiter_option
       write(iulog,*)"readnl: dt_tracer_factor = ",dt_tracer_factor
       write(iulog,*)"readnl: vertical remap frequency dt_remap_factor (0=disabled): ",dt_remap_factor

       write(iulog,*)"readnl: se_fv_phys_remap_alg = ",se_fv_phys_remap_alg

       if (hypervis_power /= 0)then
          write(iulog,*)"Variable scalar hyperviscosity: hypervis_power=",hypervis_power
          write(iulog,*)"max_hypervis_courant = ", max_hypervis_courant
       elseif(hypervis_scaling /=0)then
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

      write(iulog,*) '  vform =',trim(vform)
      write(iulog,*) '  vfile_mid=',trim(vfile_mid)
      write(iulog,*) '  vfile_int=',trim(vfile_int)

      ! display physical constants for HOMME stand alone simulations
      write(iulog,*)""
      write(iulog,*)"physconst: omega  = ",omega
      write(iulog,*)"physconst: rearth = ",rearth
      write(iulog,*)"physconst: rrearth= ",rrearth
      write(iulog,*)""
    endif
    
    is_params_inited = .true.
  end subroutine init_params_f90

  function get_homme_int_param_f90 (param_name_c) result(param_value) bind(c)
    use dimensions_mod,    only: nelemd, qsize, nlev, np
    use control_mod,       only: ftype
    use time_mod,          only: nmax
    use homme_context_mod, only: elem
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: param_name_c
    integer (kind=c_int) :: param_value
    !
    ! Local(s)
    !
    character(len=256), pointer :: param_name
    integer :: len
    integer :: dims(4)

    call c_f_pointer(param_name_c,param_name)
    len = index(param_name, C_NULL_CHAR) -1
    select case(param_name(1:len))
      case("ftype")
        param_value = ftype
      case("qsize")
        param_value = qsize
      case("nmax")
        param_value = nmax
      case("num momentum forcings")
        dims = SHAPE(elem(1)%derived%FM)
        param_value = dims(3)
      case default
        call abortmp ("[get_homme_int_param_f90] Error! Unrecognized parameter name.")
    end select 

  end function get_homme_int_param_f90

  function get_homme_real_param_f90 (param_name_c) result(param_value) bind(c)
    use control_mod,    only: nu, nu_div, nu_p, nu_q, nu_s, hypervis_scaling
    use time_mod,       only: tstep
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: param_name_c
    real (kind=c_double) :: param_value
    !
    ! Local(s)
    !
    character(len=256), pointer :: param_name
    integer :: len

    call c_f_pointer(param_name_c,param_name)
    len = index(param_name, C_NULL_CHAR) -1
    select case(param_name(1:len))
      case("nu")
        param_value = nu
      case("nu_div")
        param_value = nu_div
      case("nu_p")
        param_value = nu_p
      case("nu_q")
        param_value = nu_q
      case("nu_s")
        param_value = nu_s
      case("hypervis_scaling")
        param_value = hypervis_scaling
      case("dt")
        param_value = tstep
      case default
        call abortmp ("[get_homme_real_param_f90] Error! Unrecognized parameter name.")
    end select 

  end function get_homme_real_param_f90

  function get_homme_bool_param_f90 (param_name_c) result(param_value) bind(c)
    use control_mod,    only: moisture
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: param_name_c
    logical (kind=c_bool) :: param_value
    !
    ! Local(s)
    !
    character(len=256), pointer :: param_name
    integer :: len

    call c_f_pointer(param_name_c,param_name)
    len = index(param_name, C_NULL_CHAR) -1
    select case(param_name(1:len))
      case("moisture")
        if (moisture .eq. 'dry') then
          param_value = .false.
        else
          param_value = .true.
        endif
      case default
        call abortmp ("[get_homme_bool_param_f90] Error! Unrecognized parameter name.")
    end select 

  end function get_homme_bool_param_f90

  subroutine set_homme_int_param_f90 (param_name_c, param_value) bind(c)
    use dimensions_mod,    only: qsize, nlev, np
    use control_mod,       only: ftype, use_moisture
    use time_mod,          only: nmax
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: param_name_c
    integer (kind=c_int), intent(in) :: param_value
    !
    ! Local(s)
    !
    character(len=256), pointer :: param_name
    integer :: len

    call c_f_pointer(param_name_c,param_name)
    len = index(param_name, C_NULL_CHAR) -1
    select case(param_name(1:len))
      case("ftype")
        ftype = param_value
      case("qsize")
        qsize = param_value
        if (qsize<1) use_moisture = .false.
      case("nmax")
        nmax = param_value
      case default
        call abortmp ("[set_homme_int_param_f90] Error! Unrecognized parameter name.")
    end select 

  end subroutine set_homme_int_param_f90

  subroutine set_homme_real_param_f90 (param_name_c, param_value) bind(c)
    use control_mod,    only: nu, nu_div, nu_p, nu_q, nu_s, hypervis_scaling
    use time_mod,       only: tstep
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: param_name_c
    real (kind=c_double), intent(in) :: param_value
    !
    ! Local(s)
    !
    character(len=256), pointer :: param_name
    integer :: len

    call c_f_pointer(param_name_c,param_name)
    len = index(param_name, C_NULL_CHAR) -1
    select case(param_name(1:len))
      case("nu")
        nu = param_value
      case("nu_div")
        nu_div = param_value
      case("nu_p")
        nu_p = param_value
      case("nu_q")
        nu_q = param_value
      case("nu_s")
        nu_s = param_value
      case("hypervis_scaling")
        hypervis_scaling = param_value
      case("dt")
        tstep = param_value
      case default
        call abortmp ("[set_homme_real_param_f90] Error! Unrecognized parameter name.")
    end select 

  end subroutine set_homme_real_param_f90

  subroutine set_homme_bool_param_f90 (param_name_c, param_value) bind(c)
    use control_mod,    only: moisture
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: param_name_c
    logical (kind=c_bool), intent(in) :: param_value
    !
    ! Local(s)
    !
    character(len=256), pointer :: param_name
    integer :: len

    call c_f_pointer(param_name_c,param_name)
    len = index(param_name, C_NULL_CHAR) -1
    select case(param_name(1:len))
      case("moisture")
        if (param_value) then
          moisture = 'dry'
        else
          moisture = 'notdry'
        endif
      case default
        call abortmp ("[set_homme_bool_param_f90] Error! Unrecognized parameter name.")
    end select 

  end subroutine set_homme_bool_param_f90

end module homme_params_mod
