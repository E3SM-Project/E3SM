module drof_shr_mod

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

  public :: drof_shr_read_namelists

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

  subroutine drof_shr_read_namelists(filename, mpicom, my_task, master_task, &
       logunit, SDROF, rof_present, rof_prognostic, rofice_present, flood_present)

    ! !DESCRIPTION: Read in drof namelists
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*)       , intent(in)    :: filename          ! input namelist filename
    integer(IN)            , intent(in)    :: mpicom            ! mpi communicator
    integer(IN)            , intent(in)    :: my_task           ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task       ! task number of master task
    integer(IN)            , intent(in)    :: logunit           ! logging unit number
    type(shr_strdata_type) , intent(inout) :: SDROF
    logical                , intent(out)   :: rof_present       ! flag
    logical                , intent(out)   :: rof_prognostic    ! flag
    logical, optional      , intent(out)   :: rofice_present    ! flag
    logical, optional      , intent(out)   :: flood_present     ! flag

    !--- local variables ---
    integer(IN)   :: nunit       ! unit number
    integer(IN)   :: ierr        ! error code
    character(CL) :: decomp      ! decomp strategy - not used for NUOPC - but still needed in namelist for now

    !--- formats ---
    character(*), parameter :: F00   = "('(drof_comp_init) ',8a)"
    character(*), parameter :: F0L   = "('(drof_comp_init) ',a, l2)"
    character(*), parameter :: F01   = "('(drof_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(drof_comp_init) ',a,4es13.6)"
    character(*), parameter :: F06   = "('(drof_comp_init) ',a,5l3)"
    character(*), parameter :: subName = "(shr_drof_read_namelists) "
    !-------------------------------------------------------------------------------

    !----- define namelist -----
    namelist / drof_nml / decomp, &
         restfilm, restfils, force_prognostic_true

    !----------------------------------------------------------------------------
    ! Read drof_in
    !----------------------------------------------------------------------------

    restfilm   = trim(nullstr)
    restfils   = trim(nullstr)
    force_prognostic_true = .false.
    if (my_task == master_task) then
       nunit = shr_file_getUnit() ! get unused unit number
       open (nunit,file=trim(filename),status="old",action="read")
       read (nunit,nml=drof_nml,iostat=ierr)
       close(nunit)
       call shr_file_freeUnit(nunit)
       if (ierr > 0) then
          write(logunit,F01) 'ERROR: reading input namelist, '//trim(filename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(filename))
       end if
       write(logunit,F00)' restfilm   = ',trim(restfilm)
       write(logunit,F00)' restfils   = ',trim(restfils)
       write(logunit,F0L)' force_prognostic_true = ',force_prognostic_true
    endif
    call shr_mpi_bcast(restfilm,mpicom,'restfilm')
    call shr_mpi_bcast(restfils,mpicom,'restfils')
    call shr_mpi_bcast(force_prognostic_true,mpicom,'force_prognostic_true')

    rest_file = trim(restfilm)
    rest_file_strm = trim(restfils)

    !----------------------------------------------------------------------------
    ! Read dshr namelist
    !----------------------------------------------------------------------------

    call shr_strdata_readnml(SDROF,trim(filename),mpicom=mpicom)

    datamode = trim(SDROF%dataMode)

    ! Validate mode

    if (trim(datamode) == 'NULL' .or. &
        trim(datamode) == 'COPYALL') then
        if (my_task == master_task) then
           write(logunit,F00) 'drof datamode = ',trim(datamode)
        end if
    else
      write(logunit,F00) ' ERROR illegal drof datamode = ',trim(datamode)
      call shr_sys_abort()
   end if

    !----------------------------------------------------------------------------
    ! Determine present and prognostic flags
    !----------------------------------------------------------------------------

    rof_present = .true.
    if (trim(datamode) == 'NULL') then
       rof_present = .false.
    end if

    rof_prognostic = .false.
    if (force_prognostic_true) then
       rof_prognostic = .true.
    end if

    if (present(rofice_present)) then
       rofice_present = .false.
       if (trim(datamode) /= 'NULL') then
          rofice_present = .true.
       end if
    end if

    if (present(flood_present)) then 
       flood_present  = .false.
    end if

  end subroutine drof_shr_read_namelists

end module drof_shr_mod
