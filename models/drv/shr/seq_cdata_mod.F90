module seq_cdata_mod

  ! Common data

  use shr_kind_mod     , only: r8=> shr_kind_r8
  use shr_sys_mod      , only: shr_sys_flush
  use shr_sys_mod      , only: shr_sys_abort
  use seq_infodata_mod , only: seq_infodata_type

  use mct_mod
  use seq_comm_mct

  implicit none
  save
  private

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public seq_cdata_init
  public seq_cdata_clean
  public seq_cdata_setptrs
  public seq_cdata_info

!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------
  type seq_cdata
    !--- in general, this just groups together related data via pointers

    !--- actual memory
    character(len=16) :: name        ! user defined name
    character(len=8)  :: set         ! flag related to state of type
    integer           :: ID          ! component id
    integer           :: mpicom      ! mpi communicator

    !--- pointers
    type(mct_gGrid)  ,pointer :: dom => null()        ! domain info
    type(mct_gsMap)  ,pointer :: gsMap => null()      ! decomp info
    type(seq_infodata_type),pointer :: infodata => null() ! Input init object
  end type seq_cdata
  public seq_cdata

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------
  character(*),parameter :: isset = 'set'
  character(*),parameter :: unset = 'NOT'

!==============================================================================
contains
!==============================================================================

  subroutine seq_cdata_init(cdata,ID,dom,gsMap,infodata,name)

    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata)     ,intent(inout)       :: cdata      ! initialized
    integer             ,intent(in)          :: ID         ! component id
    type(mct_gGrid)     ,intent(in),target   :: dom        ! domain
    type(mct_gsMap)     ,intent(in),target   :: gsMap      ! decomp
    type(seq_infodata_type) &
                        ,intent(in),target   :: infodata   ! INIT object
    character(len=*)    ,intent(in),optional :: name       ! user defined name
    !
    ! Local variables
    !
    character(*),parameter :: subName = '(seq_cdata_init) '
    integer :: mpicom     ! mpi communicator
    logical :: iamroot    ! iamroot
    !-----------------------------------------------------------------------

    call seq_comm_setptrs(ID,mpicom=mpicom,iamroot=iamroot)

    if (trim(cdata%set) == trim(isset)) then
       if (iamroot) write(logunit,*) trim(subName), &
         ': WARNING cdata is already set, cleaning'
       call seq_cdata_clean(cdata)
    endif

    if (present(name)) then
      cdata%name   =  name
    else
      cdata%name   =  'undefined'
    endif

    cdata%set      =  isset
    cdata%ID       =  ID
    cdata%mpicom   =  mpicom
    cdata%dom      => dom
    cdata%gsMap    => gsMap
    cdata%infodata => infodata

  end subroutine seq_cdata_init

!==============================================================================

  subroutine seq_cdata_clean(cdata)

    implicit none
    !-----------------------------------------------------------------------

    ! Arguments
    type(seq_cdata),intent(inout)     :: cdata    ! initialized

    ! Local variables
    character(*),parameter :: subName = '(seq_cdata_clean) '

    !-----------------------------------------------------------------------

    if (trim(cdata%set) == trim(isset)) then
      cdata%name     = 'cleaned'
      nullify(cdata%dom)
      nullify(cdata%gsMap)
      nullify(cdata%infodata)
    endif

    cdata%set = unset
    cdata%ID = -1
    cdata%mpicom = -1

  end subroutine seq_cdata_clean

!==============================================================================

  subroutine seq_cdata_setptrs(cdata,ID,mpicom,dom,gsMap,infodata,name)

    implicit none
    !-----------------------------------------------------------------------

    ! Arguments
    type(seq_cdata),intent(in) :: cdata      ! input

    integer          ,optional         :: ID         ! component id
    integer          ,optional         :: mpicom     ! mpi comm
    type(mct_gGrid)  ,optional,pointer :: dom        ! domain
    type(mct_gsMap)  ,optional,pointer :: gsMap      ! decomp
    type(seq_infodata_type) &
                     ,optional,pointer :: infodata   ! INIT object
    character(len=*) ,optional         :: name       ! name

    ! Local variables
    character(*),parameter :: subName = '(seq_cdata_setptrs) '
    logical :: iamroot

    !-----------------------------------------------------------------------

    call seq_comm_setptrs(cdata%ID,iamroot=iamroot)

    if (trim(cdata%set) == trim(unset)) then
       if (iamroot) write(logunit,*) trim(subName), &
         ': WARNING cdata is not set, return'
       return
    endif

    if (present(name)) then
      name = cdata%name
    endif
    if (present(ID)) then
      ID = cdata%ID
    endif
    if (present(mpicom)) then
      mpicom = cdata%mpicom
    endif
    if (present(dom)) then
      dom => cdata%dom
    endif
    if (present(gsMap)) then
      gsMap => cdata%gsMap
    endif
    if (present(infodata)) then
      infodata => cdata%infodata
    endif

  end subroutine seq_cdata_setptrs

!==============================================================================

  subroutine seq_cdata_info(cdata)

    implicit none
    !-----------------------------------------------------------------------

    ! Arguments
    type(seq_cdata), intent(in) :: cdata

    ! Local variables
    character(*),parameter :: subName = '(seq_cdata_info) '
    logical :: iamroot

    !-----------------------------------------------------------------------

    iamroot = .false.
    if (trim(cdata%set) == trim(isset)) then
      call seq_comm_setptrs(cdata%ID,iamroot=iamroot)
    else
      write(logunit,*) subName,':',' uninitialized'
    endif

    if (iamroot) then
      write(logunit,*) ' '
      write(logunit,*) trim(subName),':'
      write(logunit,*) '  name        = ',cdata%name
      write(logunit,*) '  ID          = ',cdata%ID
      write(logunit,*) '  mpicom      = ',cdata%mpicom
      !--tc these lead to death in mct if dims or lsize are 0
      ! write(logunit,*) '  dom dims    = ',mct_gGrid_dims(cdata%dom)
      ! write(logunit,*) '  dom lsize   = ',mct_gGrid_lsize(cdata%dom)
      ! --mvr these lead to seg fault/core dump on tempest (irix64)
      ! write(logunit,*) '  dom IList   = ',mct_avect_exportIList2c(cdata%dom%data)
      ! write(logunit,*) '  dom RList   = ',mct_avect_exportRList2c(cdata%dom%data)
      write(logunit,*) '  gsMap gsize = ',mct_gsMap_gsize(cdata%gsMap)
      write(logunit,*) '  gsMap lsize = ',mct_gsMap_lsize(cdata%gsMap,cdata%mpicom)
      write(logunit,*) ' '
    endif

  end subroutine seq_cdata_info

end module seq_cdata_mod
