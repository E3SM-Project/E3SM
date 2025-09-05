!BOP ===========================================================================
!
! !MODULE: shr_moab_mod -- moab-based utilities
!
! !INTERFACE: ------------------------------------------------------------------
module shr_moab_mod

 use shr_sys_mod,      only: shr_sys_abort
 use seq_comm_mct,     only: logunit
 use shr_kind_mod,     only: CXX => shr_kind_CXX
 use shr_kind_mod    , only: R8 => SHR_KIND_R8
 use iso_c_binding
 
 implicit none
 private

 public :: mbGetnCells 
 public :: mbGetCellTagVals
 public :: mbSetCellTagVals

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

   mbGetnCells = nvise(1)

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

    ent_type=1
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

    ent_type=1
    tagname = trim(intag)//C_NULL_CHAR
    ierr = iMOAB_SetDoubleTagStorage (mbid, tagname, nMax , ent_type, inarray)
    if (ierr .ne. 0) then
        write(logunit,*) subname,' error in setting tag storage '
        call shr_sys_abort(subname//' error in setting tag storage ')
    endif

  end subroutine
end module shr_moab_mod
