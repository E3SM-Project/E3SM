module dice_shr_mod

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

  public :: dice_shr_read_namelists

  !--------------------------------------------------------------------------
  ! Public data
  !--------------------------------------------------------------------------

  ! Note that model decomp will now come from reading in the mesh directly

  ! input namelist variables
  character(CL) , public :: restfilm              ! model restart file namelist
  character(CL) , public :: restfils              ! stream restart file namelist
  real(R8)      , public :: flux_swpf             ! short-wave penatration factor
  real(R8)      , public :: flux_Qmin             ! bound on melt rate
  logical       , public :: flux_Qacc             ! activates water accumulation/melt wrt Q
  real(R8)      , public :: flux_Qacc0            ! initial water accumulation value

  ! variables obtained from namelist read
  character(CL) , public :: rest_file             ! restart filename
  character(CL) , public :: rest_file_strm        ! restart filename for streams
  character(CL) , public :: datamode              ! mode
  character(len=*), public, parameter :: nullstr = 'undefined'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine dice_shr_read_namelists(filename, mpicom, my_task, master_task, &
       logunit, SDICE, ice_present, ice_prognostic)

    ! !DESCRIPTION: Read in dice namelists
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*)       , intent(in)    :: filename       ! input namelist filename
    integer(IN)            , intent(in)    :: mpicom         ! mpi communicator
    integer(IN)            , intent(in)    :: my_task        ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task    ! task number of master task
    integer(IN)            , intent(in)    :: logunit        ! logging unit number
    type(shr_strdata_type) , intent(inout) :: SDICE
    logical                , intent(out)   :: ice_present    ! flag
    logical                , intent(out)   :: ice_prognostic ! flag

    !--- local variables ---
    integer(IN)   :: nunit       ! unit number
    integer(IN)   :: ierr        ! error code
    character(CL) :: decomp      ! decomp strategy - not used for NUOPC - but still needed in namelist for now

    !--- formats ---
    character(*), parameter :: F00   = "('(dice_comp_init) ',8a)"
    character(*), parameter :: F01   = "('(dice_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(dice_comp_init) ',a,4es13.6)"
    character(*), parameter :: F06   = "('(dice_comp_init) ',a,5l3)"
    character(*), parameter :: subName = "(shr_dice_read_namelists) "
    !-------------------------------------------------------------------------------

    !----- define namelist -----
    namelist / dice_nml / decomp, &
         flux_swpf, flux_Qmin, flux_Qacc, flux_Qacc0, restfilm, restfils

    !----------------------------------------------------------------------------
    ! Read dice_in
    !----------------------------------------------------------------------------

    flux_swpf  =     0.0_R8  ! no penetration
    flux_Qmin  =  -300.0_R8  ! kg/s/m^2
    flux_Qacc  = .false.     ! no accumulation
    flux_Qacc0 =     0.0_R8  ! no water
    restfilm   = trim(nullstr)
    restfils   = trim(nullstr)

    if (my_task == master_task) then
       nunit = shr_file_getUnit() ! get unused unit number
       open (nunit,file=trim(filename),status="old",action="read")
       read (nunit,nml=dice_nml,iostat=ierr)
       close(nunit)
       call shr_file_freeUnit(nunit)
       if (ierr > 0) then
          write(logunit,F01) 'ERROR: reading input namelist, '//trim(filename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(filename))
       end if
       write(logunit,F02)' flux_swpf  = ',flux_swpf
       write(logunit,F02)' flux_Qmin  = ',flux_Qmin
       write(logunit,F06)' flux_Qacc  = ',flux_Qacc
       write(logunit,F02)' flux_Qacc0 = ',flux_Qacc0
       write(logunit,F00)' restfilm   = ',trim(restfilm)
       write(logunit,F00)' restfils   = ',trim(restfils)
    endif

    call shr_mpi_bcast(flux_swpf ,mpicom,'flux_swpf')
    call shr_mpi_bcast(flux_Qmin ,mpicom,'flux_Qmin')
    call shr_mpi_bcast(flux_Qacc ,mpicom,'flux_Qacc')
    call shr_mpi_bcast(flux_Qacc0,mpicom,'flux_Qacc0')
    call shr_mpi_bcast(restfilm,mpicom,'restfilm')
    call shr_mpi_bcast(restfils,mpicom,'restfils')

    rest_file = trim(restfilm)
    rest_file_strm = trim(restfils)

    !----------------------------------------------------------------------------
    ! Read dshr namelist
    !----------------------------------------------------------------------------

    call shr_strdata_readnml(SDICE,trim(filename),mpicom=mpicom)

    ! Validate mode

    datamode = trim(SDICE%dataMode)
    if (trim(datamode) == 'NULL' .or. &
        trim(datamode) == 'SSTDATA' .or. &
        trim(datamode) == 'COPYALL') then
       if (my_task == master_task) then
          write(logunit,F00) ' dice datamode = ',trim(datamode)
       end if
    else
       write(logunit,F00) ' ERROR illegal dice datamode = ',trim(datamode)
       call shr_sys_abort()
    endif

    !----------------------------------------------------------------------------
    ! Determine present and prognostic flag
    !----------------------------------------------------------------------------

    if (trim(datamode) == 'NULL') then
       ice_present = .false.
       ice_prognostic = .false.
    else
       ice_present = .true.
       ice_prognostic = .true.
    end if

  end subroutine dice_shr_read_namelists

end module dice_shr_mod
