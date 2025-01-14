module homme_context_mod

  use iso_c_binding, only: c_bool, c_int

  ! This module contains homme's "global state", which other module
  ! can alter as they see fit

  ! Import homme types to be used for module-level variables
  use domain_mod,          only: domain1d_t
  use derivative_mod_base, only: derivative_t
  use element_mod,         only: element_t
  use parallel_mod,        only: parallel_t, abortmp
  use time_mod,            only: timelevel_t
  use hybvcoord_mod,       only: hvcoord_t
  use hybrid_mod,          only: hybrid_t
  use prim_driver_base,    only: deriv => deriv1

  implicit none
  private

  ! Making deriv (which is an alias to deriv1 in prim_driver_base) available
  public :: deriv

  type (element_t),    pointer, public :: elem(:) => null()
  type (domain1d_t),   pointer, public :: dom_mt(:) => null()
  type (parallel_t)           , public :: par
  type (timelevel_t)          , public :: tl
  type (hybrid_t)             , public :: hybrid
  type (hvcoord_t)            , public :: hvcoord

  integer, public :: npes        = 1
  integer, public :: iam         = 0
  logical, public :: masterproc  = .false.
  character(len=256), public :: homme_log_fname = ""
  logical :: homme_log_set = .false.

  logical, public :: is_parallel_inited          = .false.
  logical, public :: is_params_inited            = .false.
  logical, public :: is_geometry_inited          = .false.
  logical, public :: is_data_structures_inited   = .false.
  logical, public :: is_model_inited             = .false.
  logical, public :: is_hommexx_functors_inited  = .false.

  ! Functions callable from C, mostly to detect whether some parts were already inited
  public :: init_parallel_f90
  public :: is_parallel_inited_f90
  public :: is_params_inited_f90
  public :: is_geometry_inited_f90
  public :: is_data_structures_inited_f90
  public :: is_model_inited_f90
  public :: is_hommexx_functors_inited_f90
  public :: close_homme_log

contains

  subroutine set_homme_log_file_name_f90(c_str) bind(c)
    use shr_file_mod, only: shr_file_getUnit
    use iso_c_binding, only: C_NULL_CHAR, c_ptr, c_f_pointer
    use kinds, only: iulog
    type (c_ptr), intent(in) :: c_str
    !
    ! Local(s)
    !
    character(len=256), pointer :: full_name
    character(len=256) :: path, fname
    integer :: len, slash, ierr

    call c_f_pointer(c_str,full_name)
    len = index(full_name, C_NULL_CHAR) -1
    if (len>0) then
      ! Search last slash in the (trimmed) full name
      slash = index(full_name(1:len),'/',back=.true.)

      ! Note: if there's no slash (relative filename),
      ! then slash=0, and path is the empty string.
      ! Otherwise, path ends with the slash
      path = full_name(1:slash)
      fname = full_name(slash+1:len)

      homme_log_fname = trim(path)//"homme_"//fname

      iulog = shr_file_getunit()
      if (masterproc) then
        ! Create the homme log file on root rank...
        open (unit=iulog,file=trim(homme_log_fname),status='REPLACE', &
              action='WRITE', access='SEQUENTIAL', position="append")
        write(iulog,*) " ---- HOMME LOG FILE ----"
        flush(iulog)
      endif
      call mpi_barrier(par%comm,ierr)
      if (.not. masterproc) then
        ! ... and open it on all other ranks
        open (unit=iulog,file=trim(homme_log_fname),status='OLD', &
              action='WRITE', access='SEQUENTIAL', position="append")
      endif

      homme_log_set = .true.
    endif
  end subroutine set_homme_log_file_name_f90

  subroutine close_homme_log ()
    use shr_file_mod, only: shr_file_freeUnit
    use kinds, only: iulog

    if (homme_log_set) then
      close (iulog)
      call shr_file_freeUnit(iulog)
      homme_log_fname = ""
      homme_log_set = .false.
    endif
  end subroutine close_homme_log

  subroutine init_parallel_f90 (f_comm) bind(c)
    use parallel_mod,   only: initmp_from_par, abortmp
    use dimensions_mod, only: npart
    use hybrid_mod,     only: hybrid_create
    interface
      subroutine reset_cxx_comm (f_comm) bind(c)
        use iso_c_binding, only: c_int
        !
        ! Inputs
        !
        integer(kind=c_int), intent(in) :: f_comm
      end subroutine reset_cxx_comm
    end interface
    !
    ! Input(s)
    !
    integer (kind=c_int), intent(in) :: f_comm
    !
    ! Local(s)
    !
    integer :: ierr

    if (is_parallel_inited) then
      call abortmp ("Error! 'homme_init_parallel' was already called.\n")
    endif

    ! Initialize parallel structure
    par%comm = f_comm
    par%root = 0
    par%dynproc = .true.

    call MPI_comm_rank(par%comm,par%rank,ierr)
    call MPI_comm_size(par%comm,par%nprocs,ierr)

    par%masterproc = .false.
    if (par%rank .eq. par%root) par%masterproc = .true.

    call initmp_from_par(par)

    ! No horizontal threading in f90 when using C++
    hybrid = hybrid_create(par,0,1)

    ! Set number of mesh partitions equal to the number of ranks
    npart = par%nprocs

    ! Init the comm in Hommexx
    call reset_cxx_comm (f_comm)

    npes = par%nprocs
    iam  = par%rank
    masterproc = par%masterproc

    is_parallel_inited = .true.
  end subroutine init_parallel_f90

  function is_parallel_inited_f90 () result(inited) bind(c)
    logical (kind=c_bool) :: inited

    inited = LOGICAL(is_parallel_inited,kind=c_bool)
  end function is_parallel_inited_f90

  function is_params_inited_f90 () result(inited) bind(c)
    logical (kind=c_bool) :: inited

    inited = LOGICAL(is_params_inited,kind=c_bool)
  end function is_params_inited_f90

  function is_geometry_inited_f90 () result(inited) bind(c)
    logical (kind=c_bool) :: inited

    inited = LOGICAL(is_geometry_inited,kind=c_bool)
  end function is_geometry_inited_f90

  function is_data_structures_inited_f90 () result(inited) bind(c)
    logical (kind=c_bool) :: inited

    inited = LOGICAL(is_data_structures_inited,kind=c_bool)
  end function is_data_structures_inited_f90

  function is_model_inited_f90 () result(inited) bind(c)
    logical (kind=c_bool) :: inited

    inited = LOGICAL(is_model_inited,kind=c_bool)
  end function is_model_inited_f90

  function is_hommexx_functors_inited_f90 () result(inited) bind(c)
    logical (kind=c_bool) :: inited

    inited = LOGICAL(is_hommexx_functors_inited,kind=c_bool)
  end function is_hommexx_functors_inited_f90

end module homme_context_mod
