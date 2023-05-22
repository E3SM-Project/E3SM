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
    ftype,                  &
    dt_tracer_factor,       &
    dt_remap_factor

  implicit none

  public :: init_params_f90
  public :: get_homme_int_param_f90
  public :: get_homme_real_param_f90
  public :: get_homme_bool_param_f90
  public :: get_homme_nsplit_f90

  public :: set_homme_int_param_f90
  public :: set_homme_real_param_f90
  public :: set_homme_bool_param_f90

  logical :: nsplit_inited = .false.

contains

  subroutine init_params_f90 (nl_fname_c) bind(c)
    use iso_c_binding,     only: c_ptr, c_f_pointer, C_NULL_CHAR
    use namelist_mod,      only: readnl
    use shr_file_mod,      only: getunit=>shr_file_getUnit, freeunit=>shr_file_freeUnit
    use homme_context_mod, only: is_params_inited, par

    ! Inputs
    type (c_ptr), intent(in) :: nl_fname_c

    ! Locals
    character(len=256), pointer :: nl_fname_ptr
    character(len=256) :: nl_fname
    integer :: str_len

    call c_f_pointer(nl_fname_c,nl_fname_ptr)
    str_len = index(nl_fname_ptr, C_NULL_CHAR) - 1
    nl_fname = trim(nl_fname_ptr(1:str_len))

    call readnl(par,nl_fname)
    
    is_params_inited = .true.
  end subroutine init_params_f90

  function get_homme_int_param_f90 (param_name_c) result(param_value) bind(c)
    use dimensions_mod,    only: qsize, nlev, np, ne, nelem
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
      case("ne")
        param_value = ne
      case("nelem")
        param_value = nelem
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
    if (nsplit_inited) then
      ! For now, do not allow atm_dt to change throughout the simulation.
      ! This might actually be safe, and a potentially useful feature,
      ! so feel free to add it
      if (nsplit .ne. nsplit_out) then
        call abortmp ("[get_homme_nsplit_f90]\n" // &
                      "  Error! nsplit was already computed, but had a different value.\n" // &
                      "  We currently do not allow dt to change during the simulation.\n")
      endif
    else
      nsplit = nsplit_out
      nsplit_inited = .true.
    endif

  end function get_homme_nsplit_f90


  subroutine set_homme_int_param_f90 (param_name_c, param_value) bind(c)
    use dimensions_mod,    only: qsize
    use control_mod,       only: ftype, use_moisture
    use homme_context_mod, only: tl
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
      case("num_steps")
        tl%nstep = param_value
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
    use control_mod,    only: moisture, use_moisture
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
          moisture = 'notdry'
          use_moisture = .true.
        else
          moisture = 'dry'
          use_moisture = .false.
        endif
      case default
        call abortmp ("[set_homme_bool_param_f90] Error! Unrecognized parameter name.")
    end select 

  end subroutine set_homme_bool_param_f90

end module homme_params_mod
