#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module scream_homme_interface_mod
  use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_double, c_bool, C_NULL_CHAR

  ! Import homme types to be used for module-level variables
  ! TODO: some of these are only needed cause homme's f90 interfaces expect them,
  !       but they are discarded and never used in scream (e.g., domain1d).
  !       We need to find a way to eliminate these useless dependencies
  ! NOTE: perhaps we should call the init_homme_f90_structures subroutine
  !       directly from atm_mct_init
  use domain_mod,    only: domain1d_t
  use element_mod,   only: element_t
  use parallel_mod,  only: parallel_t, abortmp
  use time_mod,      only: timelevel_t
  use hybvcoord_mod, only: hvcoord_t, hvcoord_init
  use hybrid_mod,    only: hybrid_t, hybrid_create
  use perf_mod,      only: t_initf, t_finalizef, t_prf, t_startf, t_stopf

  implicit none
  private 

  type (element_t), pointer :: elem(:) => null()
  type (domain1d_t), pointer :: dom_mt(:) => null()
  type (parallel_t)  :: par
  type (timelevel_t) :: tl
  type (hybrid_t)    :: hybrid
  type (hvcoord_t)   :: hvcoord

  integer :: nets, nete, next_output_step

  logical :: is_half_inited = .false.
  logical :: is_inited = .false.

  public :: init_homme1_f90
  public :: init_homme2_f90
  public :: was_init_homme1_called_f90
  public :: was_init_homme2_called_f90
  public :: finalize_homme_f90
  public :: get_elem_cols_gids_f90
  public :: get_unique_cols_f90
  public :: get_num_owned_columns_f90

  public :: get_homme_int_param_value_f90
  public :: get_homme_real_param_value_f90
  public :: get_homme_bool_param_value_f90
contains

  subroutine init_homme1_f90 (f_comm) bind(c)
    use parallel_mod,     only: initmp_from_par
    use prim_driver_mod,  only: prim_init1, prim_create_c_data_structures
    use control_mod,      only: vfile_mid, vfile_int
    use common_io_mod,    only: infilenames
    use interpolate_driver_mod, only : pio_read_phis
    !
    ! Input(s)
    !
    integer (kind=c_int), intent(in) :: f_comm
    !
    ! Local(s)
    !
    integer :: ierr

    if (is_half_inited) then
      call abortmp ("Error! init_homme1_f90 was already called.\n")
    endif

    call init_parallel(f_comm)
    call initmp_from_par(par)

    if(par%masterproc) print *,"Primitive Equation Init1..."
    call t_initf('input.nl',LogPrint=par%masterproc, &
                 Mpicom=par%comm, MasterTask=par%masterproc)
    call t_startf('Total')

    if(par%masterproc) print *,"Entering prim_init1"

    call prim_init1(elem,par,dom_mt,tl)

    ! In c++ the threading is handled inside the c++ code.
    ! Hence, we only have 1 entry in dom_mt
    nets = dom_mt(0)%start
    nete = dom_mt(0)%end

    ! Fake hybrid, since kokkos handles hall parallelism in hommexx
    hybrid = hybrid_create(par,0,1)

    ! Init hvcoord
    hvcoord = hvcoord_init(vfile_mid, vfile_int, .true., hybrid%masterthread, ierr)

    is_half_inited = .true.

    if (infilenames(1)/='') call pio_read_phis(elem,hybrid%par)

    if(par%masterproc) print *,"Entering prim_create_c_data_structures"

    call prim_create_c_data_structures (tl, hvcoord)
  end subroutine init_homme1_f90

  subroutine init_homme2_f90 () bind(c)
    use prim_driver_mod,  only: prim_init_kokkos_states, prim_init_kokkos_functors
    use prim_driver_base, only: prim_init2
    use control_mod,      only: runtype
    use common_io_mod,    only: output_dir
#ifdef VERTICAL_INTERPOLATION
    use netcdf_interp_mod, only: netcdf_interp_init, netcdf_interp_write, netcdf_interp_finish
#endif
#ifdef PIO_INTERP
    use interp_movie_mod,  only: interp_movie_output, interp_movie_finish, interp_movie_init
#else
    use prim_movie_mod,    only: prim_movie_output, prim_movie_finish, prim_movie_init
#endif
    use common_movie_mod,  only: nextOutputStep
    !
    ! Local(s)
    !
    integer :: ierr

    if(par%masterproc) print *,"Entering init_homme2_f90"

    if (.not. is_half_inited) then
      call abortmp ("Error! init_homme1_f90 was not called yet.\n")
    endif

    if (is_inited) then
      call abortmp ("Error! init_homme2_f90 was already called.\n")
    endif

    call prim_init2(elem, hybrid, nets, nete, tl, hvcoord)

    ! Initialize Elements and Diagnostics
    call prim_init_kokkos_states (elem)

    ! Initialize the functors, including their boundary exchanges
    call prim_init_kokkos_functors ()

    ! Here we get sure the directory specified
    ! in the input namelist file in the 
    ! variable 'output_dir' does exist.
    ! this avoids a abort deep within the PIO 
    ! library (SIGABRT:signal 6) which in most
    ! architectures produces a core dump.
    if (par%masterproc) then 
       open(unit=447,file=trim(output_dir) // "/output_dir_test",iostat=ierr)
       if ( ierr==0 ) then
          print *,'Directory ',trim(output_dir), ' does exist: initialing IO'
          close(447)
       else
          print *,'Error creating file in directory ',trim(output_dir)
          call abortmp("Please be sure the directory exist or specify 'output_dir' in the namelist.")
       end if
    endif
    
    if(par%masterproc) print *,"I/O init..."
    ! initialize history files.  filename constructed with restart time
    ! so we have to do this after ReadRestart in prim_init2 above
#ifdef VERTICAL_INTERPOLATION
    call netcdf_interp_init(elem, hybrid, hvcoord)
#elif defined PIO_INTERP
    call interp_movie_init( elem, par,  hvcoord, tl )
#else
    call prim_movie_init( elem, par, hvcoord, tl )
#endif

    ! output initial state for NEW runs (not restarts or branch runs)
    if (runtype == 0 ) then
      if(par%masterproc) print *,"Output of initial state..."
#ifdef VERTICAL_INTERPOLATION
      call netcdf_interp_write(elem, tl, hybrid, hvcoord)
#elif defined(PIO_INTERP)
      call interp_movie_output(elem, tl, par, 0d0, hvcoord=hvcoord)
#else
      call prim_movie_output(elem, tl, hvcoord, par)
#endif
    endif

    next_output_step = nextOutputStep(tl)

    if(par%masterproc) print *,"Entering main timestepping loop"

    is_inited = .true.
  end subroutine init_homme2_f90

  function was_init_homme1_called_f90 () result(is_it) bind(c)
    logical(kind=c_bool) :: is_it
    is_it = is_inited
  end function was_init_homme1_called_f90

  function was_init_homme2_called_f90 () result(is_it) bind(c)
    logical(kind=c_bool) :: is_it
    is_it = is_half_inited
  end function was_init_homme2_called_f90

  subroutine run_homme_f90 (dt) bind(c)
    use control_mod,     only: restartfreq
    use dimensions_mod,  only: nelemd
    use prim_driver_mod, only: prim_run_subcycle
    use time_mod,        only: tstep
#ifdef VERTICAL_INTERPOLATION
    use netcdf_interp_mod, only: netcdf_interp_write
#elif defined PIO_INTERP
    use interp_movie_mod,  only: interp_movie_output
#else
    use prim_movie_mod,    only: prim_movie_output
#endif
    use common_movie_mod,  only: nextOutputStep
    use restart_io_mod,    only: writerestart
    !
    ! Input(s)
    !
    real (kind=c_double), intent(in) :: dt

    if (.not. is_inited) then
      call abortmp ("Error! Homme was not initialized yet (or was already finalized).\n")
    endif

    ! Set dt in the time mod
    tstep = dt

    if (par%masterproc) print *, "HOMME step: ", tl%nstep

    call prim_run_subcycle(elem,hybrid,nets,nete,dt,.false.,tl,hvcoord,1)

    if (tl%nstep .ge. next_output_step) then
#ifdef VERTICAL_INTERPOLATION
      call netcdf_interp_write(elem, tl, hybrid, hvcoord)
#elif defined PIO_INTERP
      call interp_movie_output(elem, tl, par, 0d0,hvcoord=hvcoord)
#else
      call prim_movie_output(elem, tl, hvcoord, par)
#endif

      ! ============================================================
      ! Write restart files if required 
      ! ============================================================
      if((restartfreq > 0) .and. (MODULO(tl%nstep,restartfreq) ==0)) then 
        call WriteRestart(elem, 0,1,nelemd,tl)
      endif

      next_output_step = nextOutputStep(tl)
    endif

  end subroutine run_homme_f90

  subroutine finalize_homme_f90 () bind(c)
    use prim_driver_mod,   only: prim_finalize
#ifdef VERTICAL_INTERPOLATION
    use netcdf_interp_mod, only: netcdf_interp_finish
#elif defined(PIO_INTERP)
    use interp_movie_mod,  only: interp_movie_finish
#else
    use prim_movie_mod,    only: prim_movie_finish
#endif

    if (.not. is_inited) then
      call abortmp ("Error! Homme was already finalized.\n")
    endif

#ifdef VERTICAL_INTERPOLATION
    call netcdf_interp_finish()
#elif defined PIO_INTERP
    call interp_movie_finish()
#else
    call prim_movie_finish()
#endif

    ! Deallocate the pointers
    deallocate (elem)
    deallocate (dom_mt)

    call t_stopf('Total')

    ! Finalize and write the timings
    if(par%masterproc) print *,"writing timing data"
    call t_prf('HommeTime', par%comm)
    if(par%masterproc) print *,"calling t_finalizef"
    call t_finalizef()

    is_inited = .false.
    is_half_inited = .false.

  end subroutine finalize_homme_f90

  subroutine init_parallel (f_comm)
    !
    ! Input(s)
    !
    integer, intent(in) :: f_comm
    !
    ! Local(s)
    !
    integer :: ierr

    ! Initialize parallel structure
    par%comm = f_comm
    par%root = 0
    par%dynproc = .true.

    call MPI_comm_rank(par%comm,par%rank,ierr)
    call MPI_comm_size(par%comm,par%nprocs,ierr)

    par%masterproc = .false.
    if (par%rank .eq. par%root) par%masterproc = .true.
  end subroutine init_parallel

  subroutine get_elem_cols_gids_f90 (gids_ptr) bind(c)
    use dimensions_mod, only: np,nelemd
    use kinds,          only: long_kind
    !
    ! Input(s)
    !
    type (c_ptr),         intent(in) :: gids_ptr
    !
    ! Local(s)
    !
    integer (kind=long_kind), pointer :: gids(:)
    integer :: i,j,ie

    if (.not. is_half_inited) then
      call abortmp ("Error! init_homme1_f90 was not called yet.\n")
    endif

    call c_f_pointer(gids_ptr, gids, [nelemd*np*np])
    do ie=1,nelemd
      do j=1,np
        do i=1,np
          gids(j*np+i) = elem(ie)%gdofP(i,j)
        enddo
      enddo
    enddo
  end subroutine get_elem_cols_gids_f90

  subroutine get_unique_cols_f90 (gids_ptr, elgp_ptr) bind(c)
    use dimensions_mod, only: np, nelemd
    use element_mod,    only: index_t
    use kinds,          only: long_kind
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: gids_ptr, elgp_ptr
    !
    ! Local(s)
    !
    type (index_t) :: idxP
    integer (kind=long_kind), pointer :: gids(:)
    integer (kind=c_int), pointer :: elgp(:,:)
    integer :: i,ip,jp,ie,ncols,icol

    if (.not. is_half_inited) then
      call abortmp ("Error! init_homme1_f90 was not called yet.\n")
    endif

    ncols = get_num_owned_columns_f90()
    call c_f_pointer(gids_ptr, gids, [ncols])
    call c_f_pointer(elgp_ptr, elgp, [3,ncols])

    icol = 1
    do ie=1,nelemd
      idxP = elem(ie)%idxP
      ncols = idxP%NumUniquePts
      do i=1,ncols
        ip = idxP%ia(i)
        jp = idxP%ja(i)
        elgp(1,icol) = ie
        elgp(2,icol) = jp
        elgp(3,icol) = ip

        gids(icol) = elem(ie)%gdofP(ip,jp)
      enddo
    enddo
  end subroutine get_unique_cols_f90

  function get_num_owned_columns_f90 () result (num_cols) bind(c)
    use dimensions_mod, only: nelemd
    !
    ! Local(s)
    !
    integer (kind=c_int) :: num_cols
    integer :: ie

    if (.not. is_half_inited) then
      call abortmp ("Error! init_homme1_f90 was not called yet.\n")
    endif

    num_cols = 0
    do ie=1,nelemd
      num_cols = num_cols + elem(ie)%idxP%NumUniquePts
    enddo
  end function get_num_owned_columns_f90

  function get_homme_int_param_value_f90 (param_name_c) result(param_value) bind(c)
    use dimensions_mod, only: nelemd, qsize
    use control_mod,    only: ftype
    use time_mod,       only: nmax
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

    if (.not. is_half_inited) then
      call abortmp ("Error! init_homme1_f90 was not called yet.\n")
    endif

    call c_f_pointer(param_name_c,param_name)
    len = index(param_name, C_NULL_CHAR) -1
    select case(param_name(1:len))
      case("ftype")
        param_value = ftype
      case("nelemd")
        param_value = nelemd
      case("qsize")
        param_value = qsize
      case("nmax")
        param_value = nmax
      case("num momentum forcings")
        dims = SHAPE(elem(1)%derived%FM)
        param_value = dims(3)
      case default
        call abortmp ("[get_homme_int_param_value_f90] Error! Unrecognized parameter name.")
    end select 

  end function get_homme_int_param_value_f90

  function get_homme_real_param_value_f90 (param_name_c) result(param_value) bind(c)
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

    if (.not. is_half_inited) then
      call abortmp ("Error! init_homme1_f90 was not called yet.\n")
    endif

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
        call abortmp ("[get_homme_real_param_value_f90] Error! Unrecognized parameter name.")
    end select 

  end function get_homme_real_param_value_f90

  function get_homme_bool_param_value_f90 (param_name_c) result(param_value) bind(c)
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

    if (.not. is_half_inited) then
      call abortmp ("Error! init_homme1_f90 was not called yet.\n")
    endif

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
        call abortmp ("[get_homme_bool_param_value_f90] Error! Unrecognized parameter name.")
    end select 

  end function get_homme_bool_param_value_f90

end module scream_homme_interface_mod
