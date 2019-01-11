module dwav_shr_mod

  ! !USES:

  use shr_kind_mod   , only : IN=>SHR_KIND_IN, R8=>SHR_KIND_R8
  use shr_kind_mod   , only : CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_file_mod   , only : shr_file_getunit, shr_file_freeunit
  use shr_sys_mod    , only : shr_sys_flush, shr_sys_abort
  use shr_strdata_mod, only : shr_strdata_type, shr_strdata_readnml
  use shr_mpi_mod    , only : shr_mpi_bcast

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: dwav_shr_read_namelists

  !--------------------------------------------------------------------------
  ! Public data
  !--------------------------------------------------------------------------

  ! input namelist variables
  character(CL) , public :: restfilm              ! model restart file namelist
  character(CL) , public :: restfils              ! stream restart file namelist
  logical       , public :: force_prognostic_true ! if true set prognostic true

  ! variables obtained from namelist read
  character(CL) , public :: rest_file             ! restart filename
  character(CL) , public :: rest_file_strm        ! restart filename for streams
  character(CL) , public :: datamode              ! mode
  character(len=*), public, parameter :: nullstr = 'undefined'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine dwav_shr_read_namelists(filename, mpicom, my_task, master_task, &
       logunit, SDWAV, wav_present, wav_prognostic)

    ! !DESCRIPTION: Read in dwav namelists
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*)       , intent(in)    :: filename          ! input namelist filename
    integer(IN)            , intent(in)    :: mpicom            ! mpi communicator
    integer(IN)            , intent(in)    :: my_task           ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task       ! task number of master task
    integer(IN)            , intent(in)    :: logunit           ! logging unit number
    type(shr_strdata_type) , intent(inout) :: SDWAV
    logical                , intent(out)   :: wav_present       ! flag
    logical                , intent(out)   :: wav_prognostic    ! flag

    !--- local variables ---
    integer(IN)   :: nunit       ! unit number
    integer(IN)   :: ierr        ! error code
    character(CL) :: decomp      ! decomp strategy - not used for NUOPC - but still needed in namelist for now

    !--- formats ---
    character(*), parameter :: F00   = "('(dwav_comp_init) ',8a)"
    character(*), parameter :: F0L   = "('(dwav_comp_init) ',a, l2)"
    character(*), parameter :: F01   = "('(dwav_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(dwav_comp_init) ',a,4es13.6)"
    character(*), parameter :: F06   = "('(dwav_comp_init) ',a,5l3)"
    character(*), parameter :: subName = "(shr_dwav_read_namelists) "
    !-------------------------------------------------------------------------------

    !----- define namelist -----
    namelist / dwav_nml / decomp, &
         restfilm, restfils, force_prognostic_true

    decomp     = "1d"
    restfilm   = trim(nullstr)
    restfils   = trim(nullstr)
    force_prognostic_true = .false.

    if (my_task == master_task) then
       nunit = shr_file_getUnit() ! get unused unit number
       open (nunit,file=trim(filename),status="old",action="read")
       read (nunit,nml=dwav_nml,iostat=ierr)
       close(nunit)
       call shr_file_freeUnit(nunit)
       if (ierr > 0) then
          write(logunit,F01) 'ERROR: reading input namelist, '//trim(filename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(filename))
       end if
       write(logunit,F00)' decomp     = ',trim(decomp)
       write(logunit,F00)' restfilm   = ',trim(restfilm)
       write(logunit,F00)' restfils   = ',trim(restfils)
       write(logunit,F0L)' force_prognostic_true = ',force_prognostic_true
    endif

    call shr_mpi_bcast(decomp  ,mpicom,'decomp')
    call shr_mpi_bcast(restfilm,mpicom,'restfilm')
    call shr_mpi_bcast(restfils,mpicom,'restfils')
    call shr_mpi_bcast(force_prognostic_true,mpicom,'force_prognostic_true')

    rest_file = trim(restfilm)
    rest_file_strm = trim(restfils)

    !----------------------------------------------------------------------------
    ! Read dshr namelist
    !----------------------------------------------------------------------------

    call shr_strdata_readnml(SDWAV,trim(filename),mpicom=mpicom)

    !----------------------------------------------------------------------------
    ! Determine and validate datamode
    !----------------------------------------------------------------------------

    datamode = trim(SDWAV%dataMode)

    if ( trim(datamode) == 'NULL' .or. &
         trim(datamode) == 'COPYALL') then
       if (my_task == master_task) then
          write(logunit,F00) 'dwav datamode = ',trim(datamode)
       end if
    else
       write(logunit,F00) ' ERROR illegal dwav datamode = ',trim(datamode)
       call shr_sys_abort()
    end if

    !----------------------------------------------------------------------------
    ! Determine present and prognostic flags
    !----------------------------------------------------------------------------

    wav_present    = .false.
    wav_prognostic = .false.
    if (force_prognostic_true) then
       wav_present    = .true.
       wav_prognostic = .true.
    endif
    if (trim(datamode) /= 'NULL') then
       wav_present = .true.
    end if

  end subroutine dwav_shr_read_namelists

end module dwav_shr_mod
