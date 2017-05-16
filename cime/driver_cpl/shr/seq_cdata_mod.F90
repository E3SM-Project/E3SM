module seq_cdata_mod

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
  public :: seq_cdata_setptrs
  public :: seq_cdata_init    ! only used by xxx_comp_esm.F90 for data models

  !--------------------------------------------------------------------------
  ! Public data
  !--------------------------------------------------------------------------
  ! in general, this type just groups together related data via pointers

  type seq_cdata
     character(len=16)                :: name               ! user defined name
     integer                          :: ID                 ! component id
     integer                          :: mpicom             ! mpi communicator
     type(mct_gGrid)         ,pointer :: dom => null()      ! domain info
     type(mct_gsMap)         ,pointer :: gsMap => null()    ! decomp info
     type(seq_infodata_type) ,pointer :: infodata => null() ! Input init object
  end type seq_cdata

  public seq_cdata

!==============================================================================
contains
!==============================================================================

  subroutine seq_cdata_setptrs(cdata, ID, mpicom, dom, gsMap, infodata, name)

    !-----------------------------------------------------------------------
    !
    ! Arguments
    type(seq_cdata)         ,intent(in)       :: cdata      ! input
    integer                 ,optional         :: ID         ! component id
    integer                 ,optional         :: mpicom     ! mpi comm
    type(mct_gGrid)         ,optional,pointer :: dom        ! domain
    type(mct_gsMap)         ,optional,pointer :: gsMap      ! decomp
    type(seq_infodata_type) ,optional,pointer :: infodata   ! INIT object
    character(len=*)        ,optional         :: name       ! name
    !
    ! Local variables
    character(*),parameter :: subName = '(seq_cdata_setptrs) '
    !-----------------------------------------------------------------------

    if (present(name     )) name     =  cdata%name
    if (present(ID       )) ID       =  cdata%ID
    if (present(mpicom   )) mpicom   =  cdata%mpicom
    if (present(dom      )) dom      => cdata%dom
    if (present(gsMap    )) gsMap    => cdata%gsMap
    if (present(infodata )) infodata => cdata%infodata

  end subroutine seq_cdata_setptrs

  !===============================================================================

  subroutine seq_cdata_init(cdata,ID,dom,gsMap,infodata,name)

    !-----------------------------------------------------------------------
    ! Description
    ! This is here only for backwards compatibility with current data model
    ! xxx_comp_esmf.F90 interfaces
    !
    ! Arguments
    implicit none
    type(seq_cdata)         ,intent(inout)       :: cdata      ! initialized
    integer                 ,intent(in)          :: ID         ! component id
    type(mct_gGrid)         ,intent(in),target   :: dom        ! domain
    type(mct_gsMap)         ,intent(in),target   :: gsMap      ! decomp
    type(seq_infodata_type) ,intent(in),target   :: infodata   ! INIT object
    character(len=*)        ,intent(in),optional :: name       ! user defined name
    !
    ! Local variables
    !
    integer :: mpicom     ! mpi communicator
    character(*),parameter :: subName = '(seq_cdata_init) '
    logical :: iamroot    ! iamroot
    !-----------------------------------------------------------------------

    call seq_comm_setptrs(ID, mpicom=mpicom, iamroot=iamroot)

    if (present(name)) then
      cdata%name   =  name
    else
      cdata%name   =  'undefined'
    endif
    cdata%ID       =  ID
    cdata%mpicom   =  mpicom
    cdata%dom      => dom
    cdata%gsMap    => gsMap
    cdata%infodata => infodata

  end subroutine seq_cdata_init

end module seq_cdata_mod
