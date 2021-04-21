module dlnd_shr_mod

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

  public :: dlnd_shr_read_namelists

  !--------------------------------------------------------------------------
  ! Public data
  !--------------------------------------------------------------------------

  ! input namelist variables
  character(CL) , public :: decomp                ! decomp strategy
  character(CL) , public :: restfilm              ! model restart file namelist
  character(CL) , public :: restfils              ! stream restart file namelist
  logical       , public :: force_prognostic_true ! if true set prognostic true
  character(CL) , public :: domain_fracname       ! name of fraction field on first stream file ('null' implies this is unused)

  ! variables obtained from namelist read
  character(CL) , public :: rest_file             ! restart filename
  character(CL) , public :: rest_file_strm        ! restart filename for streams
  character(CL) , public :: datamode              ! mode
  character(len=*), public, parameter :: nullstr = 'undefined'
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine dlnd_shr_read_namelists(mpicom, my_task, master_task, &
       inst_index, inst_suffix, inst_name, &
       logunit, shrlogunit, SDLND, lnd_present, lnd_prognostic)

    ! !DESCRIPTION: Read in dlnd namelists
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    integer(IN)            , intent(in)    :: mpicom            ! mpi communicator
    integer(IN)            , intent(in)    :: my_task           ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task       ! task number of master task
    integer                , intent(in)    :: inst_index        ! number of current instance (ie. 1)
    character(len=16)      , intent(in)    :: inst_suffix       ! char string associated with instance
    character(len=16)      , intent(in)    :: inst_name         ! fullname of current instance (ie. "lnd_0001")
    integer(IN)            , intent(in)    :: logunit           ! logging unit number
    integer(IN)            , intent(in)    :: shrlogunit        ! original log unit and level
    type(shr_strdata_type) , intent(inout) :: SDLND
    logical                , intent(out)   :: lnd_present       ! flag
    logical                , intent(out)   :: lnd_prognostic    ! flag

    !--- local variables ---
    character(CL) :: fileName    ! generic file name
    integer(IN)   :: nunit       ! unit number
    integer(IN)   :: ierr        ! error code

    !--- formats ---
    character(*), parameter :: F00   = "('(dlnd_comp_init) ',8a)"
    character(*), parameter :: F0L   = "('(dlnd_comp_init) ',a, l2)"
    character(*), parameter :: F01   = "('(dlnd_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(dlnd_comp_init) ',a,4es13.6)"
    character(*), parameter :: F06   = "('(dlnd_comp_init) ',a,5l3)"
    character(*), parameter :: subName = "(shr_dlnd_read_namelists) "
    !-------------------------------------------------------------------------------

    !----- define namelist -----
    namelist / dlnd_nml / &
         decomp, restfilm, restfils, force_prognostic_true, domain_fracname

    !----------------------------------------------------------------------------
    ! Determine input filenamname
    !----------------------------------------------------------------------------

    filename = "dlnd_in"//trim(inst_suffix)

    !----------------------------------------------------------------------------
    ! Read dlnd_in
    !----------------------------------------------------------------------------

    filename   = "dlnd_in"//trim(inst_suffix)
    decomp     = "1d"
    restfilm   = trim(nullstr)
    restfils   = trim(nullstr)
    force_prognostic_true = .false.
    domain_fracname = trim(nullstr)
    if (my_task == master_task) then
       nunit = shr_file_getUnit() ! get unused unit number
       open (nunit,file=trim(filename),status="old",action="read")
       read (nunit,nml=dlnd_nml,iostat=ierr)
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
       write(logunit,F00)' domain_fracname = ',trim(domain_fracname)
    endif
    call shr_mpi_bcast(decomp  ,mpicom,'decomp')
    call shr_mpi_bcast(restfilm,mpicom,'restfilm')
    call shr_mpi_bcast(restfils,mpicom,'restfils')
    call shr_mpi_bcast(force_prognostic_true,mpicom,'force_prognostic_true')
    call shr_mpi_bcast(domain_fracname,mpicom,'domain_fracname')

    rest_file = trim(restfilm)
    rest_file_strm = trim(restfils)

    !----------------------------------------------------------------------------
    ! Read dshr namelist
    !----------------------------------------------------------------------------

    call shr_strdata_readnml(SDLND,trim(filename),mpicom=mpicom)

    !----------------------------------------------------------------------------
    ! Determine and validate datamode
    !----------------------------------------------------------------------------

    datamode = trim(SDLND%dataMode)

    if (trim(datamode) == 'NULL' .or. &
         trim(datamode) == 'COPYALL') then
       if (my_task == master_task) then
          write(logunit,F00) 'dlnd datamode = ',trim(datamode)
       end if
    else
       write(logunit,F00) ' ERROR illegal dlnd datamode = ',trim(datamode)
       call shr_sys_abort()
    end if

    !----------------------------------------------------------------------------
    ! Determine present and prognostic flags
    !----------------------------------------------------------------------------

    lnd_present    = .false.
    lnd_prognostic = .false.
    if (force_prognostic_true) then
       lnd_present    = .true.
       lnd_prognostic = .true.
    endif
    if (trim(datamode) /= 'NULL') then
       lnd_present = .true.
    end if

  end subroutine dlnd_shr_read_namelists

end module dlnd_shr_mod
