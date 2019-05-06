module inic_analytic_utils

  !-----------------------------------------------------------------------
  !
  ! Purpose: Set analytic initial conditions based on input coordinates
  !
  !
  !-----------------------------------------------------------------------
  use cam_logfile,         only: iulog
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use cam_abortutils,      only: endrun
  use shr_sys_mod,         only: shr_sys_flush

  implicit none
  private

  ! Public interfaces
  public :: analytic_ic_readnl   ! Read dyn_test_nl namelist
  public :: analytic_ic_active   ! .true. if analytic IC should be set
  public :: analytic_ic_is_moist ! .true. if IC are moist

  ! Private module variables
  integer,                   parameter         :: scheme_len = 32
  logical                                      :: moist = .false.

  ! Protected resource
  character(len=scheme_len), public, protected :: analytic_ic_type = 'none'

!==============================================================================
CONTAINS
!==============================================================================

  logical function analytic_ic_active()
    analytic_ic_active = (trim(analytic_ic_type) /= 'none')
  end function analytic_ic_active

  logical function analytic_ic_is_moist()
    analytic_ic_is_moist = moist
  end function analytic_ic_is_moist

  subroutine analytic_ic_readnl(nlfile)

    use namelist_utils, only: find_group_name
    use units,          only: getunit, freeunit
    use spmd_utils,     only: masterproc, masterprocid, mpicom, mpi_character, mpi_logical
    use shr_string_mod, only: shr_string_toLower

    ! Dummy argument
    character(len=*), intent(in)   :: nlfile  ! filepath of namelist input file

    !
    ! Local variables
    integer                        :: unitn, ierr
    logical                        :: nl_not_found
    character(len=128)             :: msg
    character(len=*), parameter    :: subname = 'ANALYTIC_IC_READNL'

#ifdef ANALYTIC_IC
    ! History namelist items
    namelist /analytic_ic_nl/ analytic_ic_type

    if (masterproc) then
      unitn = getunit()
      open(unitn, file=trim(nlfile), status='old')
      call find_group_name(unitn, 'analytic_ic_nl', status=ierr)
      if (ierr == 0) then
        nl_not_found = .false.
        write(iulog, *) 'Read in analytic_ic_nl namelist from: ',trim(nlfile)
        read(unitn, analytic_ic_nl, iostat=ierr)
        if (ierr /= 0) then
          write(msg, '(a,i0)')                                             &
               ': ERROR reading namelist, analytic_ic_nl, iostat = ', ierr
          call endrun(subname//trim(msg))
        end if
      else
        nl_not_found = .true.
      end if
      close(unitn)
      call freeunit(unitn)

      analytic_ic_type = shr_string_toLower(analytic_ic_type)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(analytic_ic_type, len(analytic_ic_type), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(nl_not_found, 1, mpi_logical, masterprocid, mpicom, ierr)

    if (nl_not_found) then
       ! If analytic IC functionality is turned on (via a configure switch), then
       ! build-namelist supplies the namelist group.  If not found then nothing
       ! to do.
       return
    else
       select case(trim(analytic_ic_type))
       case('held_suarez_1994')
          msg = 'Dynamics state will be set to Held-Suarez (1994) initial conditions.'
       case('baroclinic_wave')
          moist = .true.
          msg = 'Dynamics state will be set to a baroclinic wave initial condition.'
       case('dry_baroclinic_wave')
          moist = .false.
          msg = 'Dynamics state will be set to a dry baroclinic wave initial condition.'
       case('none')
          msg = subname//': ERROR: analytic_ic_type must be set'
          write(iulog, *) msg
          call endrun(msg)
       case default
          msg = subname//': ERROR: analytic_ic_type not recognized: '//trim(analytic_ic_type)
          write(iulog, *) msg
          call endrun(msg)
       end select

    end if

    ! Write out initial condition scheme info
    if (masterproc) then
       write(iulog, *) msg
    end if
#else
    analytic_ic_type = 'none'
    moist = .false.
#endif

  end subroutine analytic_ic_readnl

end module inic_analytic_utils
