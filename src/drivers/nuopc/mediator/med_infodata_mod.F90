module med_infodata_mod

  ! !DESCRIPTION: A module to get, put, and store some standard scalar data

  ! !USES:

  use med_constants_mod , only: CL, R8
  use seq_comm_mct      , only: num_inst_atm, num_inst_lnd, num_inst_rof
  use seq_comm_mct      , only: num_inst_ocn, num_inst_ice, num_inst_glc, num_inst_wav
  use esmFlds           , only: ncomps

  implicit none
  private  ! default private

  ! !PUBLIC TYPES:

  public :: med_infodata_type

  ! !PUBLIC MEMBER FUNCTIONS

  public :: med_infodata_init1           ! Initialize before clocks are initialized
  public :: med_infodata_init2           ! Init after clocks are initialized
  public :: med_infodata_GetData         ! Get values from infodata object
  public :: med_infodata_CopyStateToInfodata
  public :: med_infodata_CopyInfodataToState

  ! !PUBLIC DATA MEMBERS:
  public :: med_infodata                  ! instance of infodata datatype

  ! Type to hold pause/resume signaling information
  type seq_pause_resume_type
     private
     character(CL) :: atm_resume(num_inst_atm) = ' ' ! atm resume file
     character(CL) :: lnd_resume(num_inst_lnd) = ' ' ! lnd resume file
     character(CL) :: ice_resume(num_inst_ice) = ' ' ! ice resume file
     character(CL) :: ocn_resume(num_inst_ocn) = ' ' ! ocn resume file
     character(CL) :: glc_resume(num_inst_glc) = ' ' ! glc resume file
     character(CL) :: rof_resume(num_inst_rof) = ' ' ! rof resume file
     character(CL) :: wav_resume(num_inst_wav) = ' ' ! wav resume file
     character(CL) :: cpl_resume = ' '               ! cpl resume file
  end type seq_pause_resume_type

  ! InputInfo derived type
  type med_infodata_type
     private
     !--- set via components and held fixed after initialization ---
     integer :: nx(ncomps) = -1              ! global nx
     integer :: ny(ncomps) = -1              ! global ny
     logical :: rofice_present = .false.     ! does rof have iceberg coupling on
     logical :: rof_prognostic = .false.     ! does rof component need input data
     logical :: flood_present = .false.      ! does rof have flooding on
     logical :: iceberg_prognostic = .false. ! does the ice model support icebergs
     logical :: glclnd_present = .false.     ! does glc have land coupling fields on
     logical :: glcocn_present = .false.     ! does glc have ocean runoff on
     logical :: glcice_present = .false.     ! does glc have iceberg coupling on
     logical :: glc_coupled_fluxes = .false. ! does glc send fluxes to other components
                                             ! (only relevant if glc_present is .true.)

     !--- set via components and may be time varying ---
     real(R8)       :: nextsw_cday = -1.0_R8 ! calendar of next atm shortwave
     real(R8)       :: precip_fact =  1.0_R8 ! precip factor
     type(seq_pause_resume_type), pointer :: pause_resume => NULL()

     !--- set by driver and may be time varying
     logical                 :: glc_valid_input = .true. ! is valid accumulated data being sent to prognostic glc
  end type med_infodata_type

  type (med_infodata_type), target :: med_infodata ! single instance for cpl and all comps

  ! used/reused in module

  character(*),parameter :: u_FILE_u = &
    __FILE__

!===============================================================================
CONTAINS
!===============================================================================

  subroutine med_infodata_init1(infodata)

    ! !DESCRIPTION:
    ! Initialize pause_resume

    ! !INPUT/OUTPUT PARAMETERS:
    type(med_infodata_type), intent(INOUT) :: infodata  ! infodata object

    !----- local -----
    character(len=*),    parameter :: subname = '(med_infodata_Init1) '
    !-------------------------------------------------------------------------------

    if (associated(infodata%pause_resume)) then
       deallocate(infodata%pause_resume)
    end if
    nullify(infodata%pause_resume)

  end subroutine med_infodata_init1

  !===============================================================================
  subroutine med_infodata_init2(infodata)

    ! !DESCRIPTION: re-initialize pause-resume depending on the time manager setup

    ! !INPUT/OUTPUT PARAMETERS:
    type(med_infodata_type), intent(INOUT) :: infodata  ! infodata object
    !----------------------------------------------------------

    ! TODO: this must be re-implemented for nuopc
    !| If pause/resume is active, initialize the resume data
    allocate(infodata%pause_resume)

  end subroutine med_infodata_init2

  !================================================================================
  subroutine med_infodata_CopyStateToInfodata(State, infodata, type, vm, rc)
    use ESMF                  , only : ESMF_State, ESMF_Field, ESMF_StateItem_Flag
    use ESMF                  , only : ESMF_StateGet, ESMF_FieldGet, ESMF_LogWrite
    use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LOGMSG_INFO
    use ESMF                  , only : ESMF_STATEITEM_NOTFOUND, operator(==)
    use ESMF                  , only : ESMF_VMBroadCast, ESMF_VM, ESMF_VMGet
    use esmFlds               , only : compname
    use shr_nuopc_scalars_mod , only : flds_scalar_num, flds_scalar_name
    use shr_nuopc_scalars_mod , only : flds_scalar_index_nx, flds_scalar_index_ny
    use shr_nuopc_scalars_mod , only : flds_scalar_index_nextsw_cday
    use shr_nuopc_scalars_mod , only : flds_scalar_index_flood_present
    use shr_nuopc_scalars_mod , only : flds_scalar_index_rofice_present
    use shr_nuopc_scalars_mod , only : flds_scalar_index_precip_fact
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_chkErr

    ! ----------------------------------------------
    ! Copy scalar data from State to local data on root then broadcast data
    ! to all PETs in component.
    ! ----------------------------------------------

    type(ESMF_State),        intent(in)     :: State
    type(med_infodata_type), intent(inout)  :: infodata
    character(len=*),        intent(in)     :: type
    type(ESMF_VM)                           :: vm
    integer,                 intent(inout)  :: rc

    ! local variables
    integer                         :: n
    integer                         :: mytask, ierr, len
    type(ESMF_Field)                :: field
    type(ESMF_StateItem_Flag)       :: itemType
    real(R8), pointer               :: farrayptr(:,:)
    real(R8)                        :: data(flds_scalar_num)
    character(len=32)               :: ntype
    integer                         :: dbrc
    character(len=1024)             :: msgString
    character(len=*), parameter     :: subname='(med_infodata_CopyStateToInfodata)'
    !----------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateGet(State, itemName=trim(flds_scalar_name), itemType=itemType, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (itemType == ESMF_STATEITEM_NOTFOUND) then
       call ESMF_LogWrite(trim(subname)//": "//trim(flds_scalar_name)//" not found", ESMF_LOGMSG_INFO, &
            line=__LINE__, file=u_FILE_u, rc=dbrc)
    else
      call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=field, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      if (mytask == 0) then
        call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (size(data) < flds_scalar_num .or. size(farrayptr) < flds_scalar_num) then
          call ESMF_LogWrite(trim(subname)//": ERROR on data size", ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif
        data(1:flds_scalar_num) = farrayptr(1,1:flds_scalar_num)
      endif

      call ESMF_VMBroadCast(vm, data, flds_scalar_num, 0, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      do n = 1,ncomps
         ntype = trim(compname(n))//'2cpli'
         if (trim(type) == trim(ntype)) then
            infodata%nx(n) = nint(data(flds_scalar_index_nx))
            infodata%ny(n) = nint(data(flds_scalar_index_ny))
            write(msgString,'(2i8,2l4)') nint(data(flds_scalar_index_nx)),nint(data(flds_scalar_index_ny))
            call ESMF_LogWrite(trim(subname)//":"//trim(type)//":"//trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
         endif
      enddo

      if (type == 'atm2cpli') then
        infodata%nextsw_cday = data(flds_scalar_index_nextsw_cday)

      elseif (type == 'ocn2cpli') then
        infodata%precip_fact=data(flds_scalar_index_precip_fact)

      elseif (type == 'ice2cpli') then
        ! nothing

      elseif (type == 'lnd2cpli') then
        ! nothing

      elseif (type == 'rof2cpli') then
        infodata%flood_present=(nint(data(flds_scalar_index_flood_present)) /= 0)
        infodata%rofice_present=(nint(data(flds_scalar_index_rofice_present)) /= 0)

      elseif (type == 'wav2cpli') then
        ! nothing

      elseif (type == 'glc2cpli') then
        ! nothing

      elseif (type == 'atm2cpl') then
         infodata%nextsw_cday=data(flds_scalar_index_nextsw_cday)

      elseif (type == 'ocn2cpl') then
         infodata%precip_fact=data(flds_scalar_index_precip_fact)

      elseif (type == 'ice2cpl') then
        ! nothing

      elseif (type == 'lnd2cpl') then
        ! nothing

      elseif (type == 'rof2cpl') then
        ! nothing

      elseif (type == 'wav2cpl') then
        ! nothing

      elseif (type == 'glc2cpl') then
        ! nothing

      else
         call ESMF_LogWrite(trim(subname)//": ERROR in type = "//trim(type), &
              ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u, rc=dbrc)
        rc = ESMF_FAILURE
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif

    endif

  end subroutine med_infodata_CopyStateToInfodata

  !================================================================================
  subroutine med_infodata_CopyInfodataToState(infodata, State, type, mytask, rc)
    use ESMF                  , only : ESMF_State, ESMF_StateGet, ESMF_Field, ESMF_StateItem_Flag, ESMF_FieldGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_STATEITEM_NOTFOUND
    use ESMF                  , only : operator(==), ESMF_FAILURE
    use shr_nuopc_scalars_mod , only : flds_scalar_num, flds_scalar_name
    use shr_nuopc_scalars_mod , only : flds_scalar_index_nx, flds_scalar_index_ny
    use shr_nuopc_scalars_mod , only : flds_scalar_index_nextsw_cday
    use shr_nuopc_scalars_mod , only : flds_scalar_index_precip_fact
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_chkErr

    ! ----------------------------------------------
    ! Copy local scalar data into State, root only,
    ! but called on all PETs in component
    ! ----------------------------------------------

    type(med_infodata_type),intent(in):: infodata
    type(ESMF_State),  intent(inout)  :: State
    character(len=*),  intent(in)     :: type
    integer         ,  intent(in)     :: mytask
    integer,           intent(inout)  :: rc

    ! local variables
    type(ESMF_Field)            :: field
    type(ESMF_StateItem_Flag)   :: ItemType
    real(R8), pointer :: farrayptr(:,:)
    real(R8)          :: nextsw_cday, precip_fact
    integer :: dbrc
    character(len=*), parameter :: subname='(med_infodata_CopyInfodataToState)'
    !----------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(flds_scalar_name), itemType=itemType, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (itemType == ESMF_STATEITEM_NOTFOUND) then

       call ESMF_LogWrite(trim(subname)//": "//trim(flds_scalar_name)//" not found", &
            ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u, rc=dbrc)

    else

      call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=field, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      if (mytask == 0) then
        call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        if (size(farrayptr) < flds_scalar_num) then
           call ESMF_LogWrite(trim(subname)//": ERROR on data size", &
                ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif
        farrayptr(1,flds_scalar_index_nextsw_cday) = infodata%nextsw_cday
        farrayptr(1,flds_scalar_index_precip_fact) = infodata%precip_fact
      endif

    endif

  end subroutine med_infodata_CopyInfodataToState

  !===============================================================================
  subroutine med_infodata_GetData( infodata, ncomp, flux_epbal, flux_epbalfact, nx, ny)

    ! Get values out of the infodata object.
    use med_constants_mod     , only : CL, IN
    use med_internalstate_mod , only : logunit, loglevel
    use shr_sys_mod           , only : shr_sys_abort

    ! !INPUT/OUTPUT PARAMETERS:
    type(med_infodata_type) , intent(IN)  :: infodata       ! Input CCSM structure
    integer(IN),   optional , intent(IN)  :: ncomp          ! Component ID
    character(CL), optional , intent(IN)  :: flux_epbal     ! selects E,P,R adjustment technique
    real(R8),      optional , intent(OUT) :: flux_epbalfact ! adjusted precip factor
    integer(IN),   optional , intent(OUT) :: nx             ! nx
    integer(IN),   optional , intent(OUT) :: ny             ! ny

    !----- local -----
    character(len=*), parameter :: subname = '(med_infodata_GetData) '
    !-------------------------------------------------------------------------------

    if ( present(flux_epbalfact)) then
       if (.not. present(flux_epbal)) then
          call shr_sys_abort(subname // "Must provide flux_epbal as an input argument to determine infodata%precip_fact")
       end if

       flux_epbalfact = 1.0_R8
       if (trim(flux_epbal) == 'ocn') then
          flux_epbalfact = infodata%precip_fact
          if (flux_epbalfact <= 0.0_R8) then
             if (loglevel > 0) then
                write(logunit,'(2a,e16.6)') trim(subname),' WARNING: factor from ocn = ',flux_epbalfact
                write(logunit,'(2a)') trim(subname),' WARNING: resetting flux_epbalfact to 1.0'
             end if
             flux_epbalfact = 1.0_R8
          end if
       end if
    endif

    if (present(nx)) then
       if (.not.present(ncomp)) then
          call shr_sys_abort(subname // " Must provide nx")
       endif
       nx = infodata%nx(ncomp)
    endif

    if (present(ny)) then
       if (.not.present(ncomp)) then
          call shr_sys_abort(subname // "Must provide ny")
       endif
       ny = infodata%ny(ncomp)
    endif

  end subroutine med_infodata_GetData

end module med_infodata_mod
