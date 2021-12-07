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
    qsize_d,                &
    npart

  use cube_mod, only: rotate_grid
  use time_mod, only: tstep, nsplit
  use thread_mod,   only: nthreads, vthreads      ! Unused in SCREAM
  use parallel_mod, only: mpireal_t, mpilogical_t, mpiinteger_t, mpichar_t
  use physical_constants, only : scale_factor, scale_factor_inv, domain_size, laplacian_rigid_factor, DD_PI

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
    geometry,               &         ! Mesh geometry
    tstep_type,             &
    ftype,                  &
    dt_tracer_factor,       &
    dt_remap_factor,        &
    rsplit,                 &
    qsplit,                 &
    limiter_option,         &
    moisture,               &
    use_moisture,           &
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
    u_perturb,              &
    se_fv_phys_remap_alg,   &
    runtype,                & ! Unused in SCREAM
    statefreq,              & ! Unused in SCREAM
    vfile_mid,              &
    vfile_int,              &
    timestep_make_subcycle_parameters_consistent

  public :: init_params_f90
  public :: get_homme_int_param_f90
  public :: get_homme_real_param_f90
  public :: get_homme_bool_param_f90
  public :: get_homme_nsplit_f90

  public :: set_homme_int_param_f90
  public :: set_homme_real_param_f90
  public :: set_homme_bool_param_f90

contains

  subroutine init_params_f90 (nl_fname_c) bind(c)
    use iso_c_binding,     only: c_ptr, c_f_pointer, C_NULL_CHAR
    use shr_file_mod,      only: getunit=>shr_file_getUnit, freeunit=>shr_file_freeUnit
    use homme_context_mod, only: is_parallel_inited, is_params_inited, &
                                 par
    use params_mod, only: SFCURVE, SPHERE_COORDS

    ! Inputs
    type (c_ptr), intent(in) :: nl_fname_c

    ! Locals
    character(len=256), pointer :: nl_fname_ptr
    character(len=256) :: nl_fname
    integer :: ierr, unitn, se_ftype, str_len
    real(kind=real_kind) :: dt_max
    character(len=MAX_FILE_LEN) :: mesh_file ! Unused in SCREAM

    !-----------------------------!
    !   Declare valid namelists   !
    !-----------------------------!

    namelist /ctl_nl/test_case, &
      u_perturb,                &
      partmethod,               &         ! mesh partitioning method
      cubed_sphere_map,         &
      vert_remap_q_alg,         &
      theta_advect_form,        &
      theta_hydrostatic_mode,   &   
      ne,                       &
      tstep,                    &
      rotate_grid,              &
      topology,                 &         ! Mesh topology
      tstep_type,               &
      dt_tracer_factor,         &
      dt_remap_factor,          &
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
      statefreq,                & ! Unused in SCREAM
      mesh_file,                & ! Unused in SCREAM
      se_ftype

    namelist /vert_nl/    &
      vfile_mid,          &
      vfile_int

    ! Sanity check (we need to bcast options to all ranks)
    if (.not. is_parallel_inited) then
      call abortmp ("Error! 'init_parallel_f90' was not called yet.\n")
    endif

    !-----------------------------!
    !        Set defaults         !
    !-----------------------------!

    nsplit = -1
    tstep = 0
    moisture = 'dry'
    cubed_sphere_map = 2
    partmethod = SFCURVE
    coord_transform_method = SPHERE_COORDS
    transport_alg = 0
    runtype = 0
    statefreq = 99999
    geometry = "sphere"
    se_ftype = ftype

    !-----------------------------!
    !     Parse namelist file     !
    !-----------------------------!

    ! Open namelist file
    unitn = getunit()
    call c_f_pointer(nl_fname_c,nl_fname_ptr)
    str_len = index(nl_fname_ptr, C_NULL_CHAR) - 1
    nl_fname = trim(nl_fname_ptr(1:str_len))
    open( unitn, file=nl_fname, status='old' )

    print *, "Reading namelist options from file:", nl_fname

    ! Parse all namelist sections
    read (unit=unitn,nml=ctl_nl,iostat=ierr)
    if (ierr < 0) then
      write(6,*) 'ierr =',ierr
      call abortmp( '[scream] init_params_f90 :: namelist read returns an'// &
             ' end of file or end of record condition' )
    end if
    read(unit=unitn,nml=vert_nl)
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
    call MPI_bcast(geometry, MAX_STRING_LEN, MPIChar_t, par%root, par%comm, ierr)

    ! Time-stepping params
    call MPI_bcast(tstep,            1, MPIreal_t,    par%root, par%comm, ierr)
    call MPI_bcast(tstep_type,       1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(dt_tracer_factor, 1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(dt_remap_factor,  1, MPIinteger_t, par%root, par%comm, ierr)
    call MPI_bcast(se_ftype,         1, MPIinteger_t, par%root, par%comm, ierr)

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
    call MPI_bcast(statefreq,           1, MPIinteger_t, par%root, par%comm, ierr)

    call MPI_bcast(test_case, MAX_STRING_LEN, MPIChar_t, par%root, par%comm, ierr)
    call MPI_bcast(moisture,  MAX_STRING_LEN, MPIChar_t, par%root, par%comm, ierr)

    ! Vertical coord params
    call MPI_bcast(vfile_mid, MAX_STRING_LEN, MPIChar_t, par%root, par%comm, ierr)
    call MPI_bcast(vfile_int, MAX_STRING_LEN, MPIChar_t, par%root, par%comm, ierr)

    !-----------------------------!
    !   Check/deduce parameters   !
    !-----------------------------!

    ierr = timestep_make_subcycle_parameters_consistent(par, rsplit, qsplit, dt_remap_factor, dt_tracer_factor)

    ftype = se_ftype
    npart = par%nprocs

    use_moisture = ( moisture /= "dry") 

    ! set map
    if (cubed_sphere_map<0) then
       cubed_sphere_map=2  ! theta model default = element local
    endif

    if (geometry=="sphere") then
      scale_factor           = rearth
      scale_factor_inv       = rrearth
      domain_size            = 4.0D0*DD_PI
      laplacian_rigid_factor = rrearth
    else
      call abortmp("Error: scream only supports 'sphere' geometry, for now.")
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

    !-----------------------------!
    !       Print parameters      !
    !-----------------------------!

    if (par%masterproc) then
       write(iulog,*)"done reading namelist..."

       write(iulog,*)"homme namelist: topology      = ",TRIM( TOPOLOGY )

       write(iulog,*)"homme namelist: test_case     = ",TRIM(test_case)
       write(iulog,*)"homme namelist: omega         = ",omega

       write(iulog,*)"homme namelist: qsize_d = ",qsize_d

       write(iulog,*)"homme namelist: ne,np         = ",ne,np
       write(iulog,*)"homme namelist: partmethod    = ",PARTMETHOD
       write(iulog,*)"homme namelist: cubed_sphere_map = ",cubed_sphere_map
       write(iulog,*)"homme namelist: COORD_TRANSFORM_METHOD    = ",COORD_TRANSFORM_METHOD

       write(iulog,*)"homme namelist: theta_hydrostatic_mode = ",theta_hydrostatic_mode
       write(iulog,*)"homme namelist: transport_alg   = ",transport_alg
       write(iulog,*)"homme namelist: theta_advect_form = ",theta_advect_form
       write(iulog,*)"homme namelist: vert_remap_q_alg  = ",vert_remap_q_alg
       write(iulog,*)"homme namelist: tstep          = ",tstep
       write(iulog,*)"homme namelist: tstep_type    = ",tstep_type
       write(iulog,*)"homme namelist: ftype          = ",ftype
       write(iulog,*)"homme namelist: limiter_option = ",limiter_option
       write(iulog,*)"homme namelist: dt_tracer_factor = ",dt_tracer_factor
       write(iulog,*)"homme namelist: vertical remap frequency dt_remap_factor (0=disabled): ",dt_remap_factor

       write(iulog,*)"homme namelist: se_fv_phys_remap_alg = ",se_fv_phys_remap_alg

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
    use dimensions_mod,    only: qsize, nlev, np
    use control_mod,       only: ftype
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

  function get_homme_nsplit_f90 (atm_dt) result(nsplit_out) bind(c)
    use control_mod, only: compute_nsplit=>timestep_make_eam_parameters_consistent
    use homme_context_mod, only: par
    !
    ! Input(s)
    !
    integer (kind=c_int) :: atm_dt
    !
    ! Local(s)
    !
    integer (kind=c_int) :: nsplit_out
    integer :: ierr, nstep_factor

    nsplit_out = -1
    ierr = compute_nsplit(par, dt_remap_factor, dt_tracer_factor, nsplit_out, nstep_factor, tstep, atm_dt)
    if (ierr .ne. 0) then
      call abortmp ("[get_homme_nsplit_f90] Error! Something went wrong in timestep_make_eam_parameters_consistent.")
    endif
    if (nsplit .eq. -1) then
      nsplit = nsplit_out
    else
      ! For now, do not allow atm_dt to change throughout the simulation.
      ! This might actually be safe, and a potentially useful feature,
      ! so feel free to add it
      if (nsplit .ne. nsplit_out) then
        call abortmp ("[get_homme_nsplit_f90]\n" // &
                      "  Error! nsplit was already computed, but had a different value.\n" // &
                      "  We currently do not allow dt to change during the simulation.\n")
      endif
    endif

  end function get_homme_nsplit_f90


  subroutine set_homme_int_param_f90 (param_name_c, param_value) bind(c)
    use dimensions_mod,    only: qsize, nlev
    use control_mod,       only: ftype, use_moisture
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
      case("ne")
        ne = param_value
      case("ftype")
        ftype = param_value
      case("qsize")
        qsize = param_value
        if (qsize>qsize_d) then
          print *, "qsize, qsize_d:", qsize, qsize_d
          call abortmp('user specified qsize > qsize_d parameter in dimensions_mod.F90')
        endif
        if (qsize<1) use_moisture = .false.
      case default
        call abortmp ("[set_homme_int_param_f90] Error! Unrecognized parameter name.")
    end select 

  end subroutine set_homme_int_param_f90

  subroutine set_homme_real_param_f90 (param_name_c, param_value) bind(c)
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: param_name_c
    real (kind=c_double), intent(in) :: param_value
    !
    ! Local(s)
    !
    character(len=256), pointer :: param_name
    character(len=256) :: param_val_str

    call c_f_pointer(param_name_c,param_name)
    write (param_val_str,'(F10.15)') param_value

    call abortmp ("[set_homme_real_param_f90] Error! This method does not currently have a valid use.\n" // &
                  "  Attempt to set " // param_name // " = " // param_val_str)

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
