!BOP ===========================================================================
!
! !MODULE: shr_moab_mod -- moab-based utilities
!
! !INTERFACE: ------------------------------------------------------------------
module shr_moab_mod

 use shr_sys_mod,      only: shr_sys_abort
 use seq_comm_mct,     only: logunit
 use seq_comm_mct,     only: atm_pg_active, mbaxid,mb_scm_land,mblxid
 use shr_kind_mod,     only: CXX => shr_kind_CXX
 use shr_kind_mod    , only: R8 => SHR_KIND_R8
 use mct_mod,          only: mct_aVect
 use iso_c_binding
 
 implicit none
 private

 public :: mbGetnCells
 public :: mbGetCellTagVals
 public :: mbSetCellTagVals
 public :: moab_set_tag_from_av
 public :: moab_set_av_from_tag

contains

!===============================================================================
! !BOP =========================================================================
!
! !IROUTINE: mbGetnCells -- return number of cells in input moab app
!
! !DESCRIPTION:
!
!     Return the number of cells in the mesh associated with the input moab app id
!
! !REVISION HISTORY:
!     2024-Jan-02 - R. Jacob - initial version
! !INTERFACE: ------------------------------------------------------------------
 integer function mbGetnCells(moabid)
   use iMOAB , only : iMOAB_GetMeshInfo

   ! Input arguments
   integer, intent(in):: moabid

   ! Local variables
   integer :: nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3) ! for moab info
   integer :: ierr

   character(*), parameter   :: subname = '(mbGetnCells) '
!-----------------------------------------------------------------------
!

   ierr  = iMOAB_GetMeshInfo ( moabid, nvert, nvise, nbl, nsurf, nvisBC );

   if (ierr .ne. 0) then
       write(logunit,*) subname,' error in getting info '
       call shr_sys_abort(subname//' error in getting info ')
   endif

   ! only case on coupler side that we actually want number of vertices is when we use spectral atm
   if ((.not. atm_pg_active .and. (mbaxid .eq. moabid)) .or. &
       (mb_scm_land .and. (mblxid .eq. moabid))) then
    mbGetnCells = nvert(1)
   else
    mbGetnCells = nvise(1)
   endif

  end function mbGetnCells

!===============================================================================
! !BOP =========================================================================
!
! !IROUTINE: mbGetCellTagVals -- get values from input tagname and mpid and return in input array
!
! !DESCRIPTION:
!
!      Get values from input tagname and mpid and return in input array
!
! !REVISION HISTORY:
!     2024-Jan-05 - R. Jacob - initial version
! !INTERFACE: ------------------------------------------------------------------

  subroutine mbGetCellTagVals(mbid, intag,inarray,nMax)
    use iMOAB , only : iMOAB_GetDoubleTagStorage

    integer,intent(in)                   :: mbid
    character(len=*),intent(in)          :: intag
    real(R8)   ,intent(inout)            :: inarray(nMax) 
    integer,intent(in)                   :: nMax

! Local variables
    integer :: ierr, ent_type
    character(CXX)           :: tagname
    character(*), parameter   :: subname = '(mbGetCellTagVals) '
!-----------------------------------------------------------------------
!
    if ((.not. atm_pg_active .and. (mbaxid .eq. mbid)) .or. &
       (mb_scm_land .and. (mblxid .eq. mbid))) then
      ent_type = 0
    else
      ent_type = 1
    endif
    
    tagname = trim(intag)//C_NULL_CHAR
    ierr = iMOAB_GetDoubleTagStorage (mbid, tagname, nMax , ent_type, inarray)
    if (ierr .ne. 0) then
        write(logunit,*) subname,' error in getting tag storage '
        call shr_sys_abort(subname//' error in getting tag storage ')
    endif


  end subroutine
!===============================================================================
! !BOP =========================================================================
!
! !IROUTINE: mbSetCellTagVals -- Set values to input tagname and mpid using input array
!
! !DESCRIPTION:
!
!      Set values to input tagname and mpid using input array
!
! !REVISION HISTORY:
!     2024-Jan-05 - R. Jacob - initial version
! !INTERFACE: ------------------------------------------------------------------

  subroutine mbSetCellTagVals(mbid, intag,inarray,nMax)
    use iMOAB , only : iMOAB_SetDoubleTagStorage

    integer,intent(in)                   :: mbid
    character(len=*),intent(in)          :: intag
    real(R8)   ,intent(inout)            :: inarray(nMax) 
    integer,intent(in)                   :: nMax

! Local variables
    integer :: ierr, ent_type
    character(CXX)           :: tagname
    character(*), parameter   :: subname = '(mbSetCellTagVals) '
!-----------------------------------------------------------------------
!

    if ((.not. atm_pg_active .and. (mbaxid .eq. mbid)) .or. &
       (mb_scm_land .and. (mblxid .eq. mbid))) then
      ent_type = 0
    else
      ent_type = 1
    endif
    tagname = trim(intag)//C_NULL_CHAR
    ierr = iMOAB_SetDoubleTagStorage (mbid, tagname, nMax , ent_type, inarray)
    if (ierr .ne. 0) then
        write(logunit,*) subname,' error in setting tag storage '
        call shr_sys_abort(subname//' error in setting tag storage ')
    endif

  end subroutine
!===============================================================================
! !BOP =========================================================================
!
! !IROUTINE: moab_set_tag_from_av -- set moab tag values from an attribute vector
!
! !DESCRIPTION:
!
!      Set field method for data models, to set moab tags from data fields in AVs
!
! !REVISION HISTORY:
!     2024-Jan-05 - R. Jacob - initial version
!     2025-Jan-18 - Moved from seq_flds_mod to shr_moab_mod
! !INTERFACE: ------------------------------------------------------------------

  subroutine moab_set_tag_from_av(tagname, avx, index, mbapid, dataarr, lsize)
    use iMOAB,        only: iMOAB_SetDoubleTagStorage

    character(len=*), intent(in) :: tagname
    type(mct_aVect), intent(in) :: avx
    integer, intent(in) :: index
    integer, intent(in) :: mbapid !  moab app id
    real(R8), intent(inout) :: dataarr(:)
    integer, intent(in) :: lsize

    ! Local variables
    integer :: ierr
    character(*), parameter :: subname = '(moab_set_tag_from_av) '
!-----------------------------------------------------------------------
!
    dataarr(:) = avx%rAttr(index, :)
    ierr = iMOAB_SetDoubleTagStorage ( mbapid, tagname, lsize, &
                                        0, & ! data on vertices
                                        dataarr )
    if (ierr > 0 )  &
        call shr_sys_abort(subname//'Error: fail to set tag values for '//tagname)

  end subroutine moab_set_tag_from_av
!===============================================================================
! !BOP =========================================================================
!
! !IROUTINE: moab_set_av_from_tag -- set attribute vector values from a moab tag
!
! !DESCRIPTION:
!
!      Get field values from moab tags and set them into an attribute vector
!
! !REVISION HISTORY:
!     2025-Jan-18 - R. Jacob - initial version
! !INTERFACE: ------------------------------------------------------------------

  subroutine moab_set_av_from_tag(tagname, avx, index, mbapid, dataarr, lsize)
    use iMOAB,        only: iMOAB_GetDoubleTagStorage

    character(len=*), intent(in) :: tagname
    type(mct_aVect), intent(inout) :: avx
    integer, intent(in) :: index
    integer, intent(in) :: mbapid !  moab app id
    real(R8), intent(inout) :: dataarr(:)
    integer, intent(in) :: lsize

    ! Local variables
    integer :: ierr
    character(*), parameter :: subname = '(moab_set_av_from_tag) '
!-----------------------------------------------------------------------
!
    ierr = iMOAB_GetDoubleTagStorage ( mbapid, tagname, lsize, &
                                        0, & ! data on vertices
                                        dataarr )
    if (ierr > 0 )  &
        call shr_sys_abort(subname//'Error: fail to get tag values for '//tagname)
    avx%rAttr(index, :) = dataarr(:)

  end subroutine moab_set_av_from_tag
end module shr_moab_mod
