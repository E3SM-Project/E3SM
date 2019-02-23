module docn_shr_mod

  ! !USES:

  use shr_kind_mod   , only : IN=>SHR_KIND_IN, R8=>SHR_KIND_R8
  use shr_kind_mod   , only : CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_file_mod   , only : shr_file_getunit, shr_file_freeunit
  use shr_sys_mod    , only : shr_sys_abort
  use shr_strdata_mod, only : shr_strdata_type, shr_strdata_readnml
  use shr_mpi_mod    , only : shr_mpi_bcast

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: docn_shr_read_namelists

  !--------------------------------------------------------------------------
  ! Public data
  !--------------------------------------------------------------------------

  ! Note that model decomp will now come from reading in the mesh directly

  ! input namelist variables
  character(CL) , public :: restfilm              ! model restart file namelist
  character(CL) , public :: restfils              ! stream restart file namelist
  logical       , public :: force_prognostic_true ! if true set prognostic true

  ! variables obtained from namelist read
  character(CL) , public :: rest_file             ! restart filename
  character(CL) , public :: rest_file_strm        ! restart filename for streams
  character(CL) , public :: datamode              ! mode
  integer(IN)   , public :: aquap_option
  character(len=*), public, parameter :: nullstr = 'undefined'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine docn_shr_read_namelists(filename, mpicom, my_task, master_task, &
       logunit, SDOCN, ocn_present, ocn_prognostic, ocnrof_prognostic)

    ! !DESCRIPTION: Read in docn namelists
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*)       , intent(in)    :: filename          ! input namelist filename
    integer(IN)            , intent(in)    :: mpicom            ! mpi communicator
    integer(IN)            , intent(in)    :: my_task           ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task       ! task number of master task
    integer(IN)            , intent(in)    :: logunit           ! logging unit number
    type(shr_strdata_type) , intent(inout) :: SDOCN
    logical                , intent(out)   :: ocn_present       ! flag
    logical                , intent(out)   :: ocn_prognostic    ! flag
    logical                , intent(out)   :: ocnrof_prognostic ! flag

    !--- local variables ---
    integer(IN)   :: nunit       ! unit number
    integer(IN)   :: ierr        ! error code
    character(CL) :: decomp      ! decomp strategy - not used for NUOPC - but still needed in namelist for now

    !--- formats ---
    character(*), parameter :: F00   = "('(docn_comp_init) ',8a)"
    character(*), parameter :: F0L   = "('(docn_comp_init) ',a, l2)"
    character(*), parameter :: F01   = "('(docn_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(docn_comp_init) ',a,4es13.6)"
    character(*), parameter :: F06   = "('(docn_comp_init) ',a,5l3)"
    character(*), parameter :: subName = "(shr_docn_read_namelists) "
    !-------------------------------------------------------------------------------

    !----- define namelist -----
    namelist / docn_nml / decomp, &
         restfilm, restfils, force_prognostic_true

    !----------------------------------------------------------------------------
    ! Read docn_in
    !----------------------------------------------------------------------------

    restfilm   = trim(nullstr)
    restfils   = trim(nullstr)
    force_prognostic_true = .false.
    if (my_task == master_task) then
       nunit = shr_file_getUnit() ! get unused unit number
       open (nunit,file=trim(filename),status="old",action="read")
       read (nunit,nml=docn_nml,iostat=ierr)
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

    call shr_strdata_readnml(SDOCN,trim(filename),mpicom=mpicom)

    datamode = trim(SDOCN%dataMode)

    ! Special logic for prescribed aquaplanet

    if (datamode(1:9) == 'SST_AQUAP' .and. trim(datamode) /= 'SST_AQUAPFILE') then
       ! First determine the prescribed aquaplanet option
       if (len_trim(datamode) == 10) then
          read(datamode(10:10),'(i1)') aquap_option
       else if (len_trim(datamode) == 11) then
          read(datamode(10:11),'(i2)') aquap_option
       end if
       ! Now remove the index from the datamode value, to have a generic setting
       ! for use below
       datamode = "SST_AQUAPANAL"
    end if

    ! Validate mode

    if ( trim(datamode) == 'NULL'          .or. &
         trim(datamode) == 'SSTDATA'       .or. &
         trim(datamode) == 'SST_AQUAPANAL' .or. &
         trim(datamode) == 'SST_AQUAPFILE' .or. &
         trim(datamode) == 'COPYALL'       .or. &
         trim(datamode) == 'IAF'           .or. &
         trim(datamode) == 'SOM'           .or. &
         trim(datamode) == 'SOM_AQUAP') then
       if (my_task == master_task) then
          write(logunit,F00) ' docn datamode = ',trim(datamode)
       end if
    else
       write(logunit,F00) ' ERROR illegal docn datamode = ',trim(datamode)
       call shr_sys_abort()
    endif

    !----------------------------------------------------------------------------
    ! Determine present and prognostic flag
    !----------------------------------------------------------------------------

    ocn_present = .true.
    if (trim(datamode) == 'NULL') then
       ocn_present = .false.
    end if

    ocn_prognostic = .false.
    if (force_prognostic_true) then
       ocn_prognostic  = .true.
    end if
    if (trim(datamode) == 'IAF') then
       ocn_prognostic = .true.
    end if
    if (trim(datamode) == 'SOM' .or. trim(datamode) == 'SOM_AQUAP') then
       ocn_prognostic = .true.
    endif

    ocnrof_prognostic = .false.
    if (force_prognostic_true .or. (trim(datamode) == 'IAF')) then
       ocnrof_prognostic = .true.
    end if
       
  end subroutine docn_shr_read_namelists

end module docn_shr_mod
