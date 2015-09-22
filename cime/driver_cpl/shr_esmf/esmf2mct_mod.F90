module esmf2mct_mod
#ifdef USE_ESMF_LIB

use ESMF
use esmfshr_util_mod, only: esmfshr_util_ArrayCopy
use mct_mod
use esmfshr_util_mod
use shr_kind_mod, only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8
use shr_sys_mod , only: shr_sys_abort
use shr_log_mod, only: shr_log_unit, shr_log_level

implicit none

! 
! Author: Fei Liu
!
! This module implements methods that work with ESMF_Array and MCT attribute
! vectors (attrVect) and global segmentation map (gsMap).
! Another layer of utils take these 
! array related methods and put them into interface blocks. user using
! the higher level utils module see/use names such as esmf2mct_init, etc
! Array, AttrVect, and gsMap APIs

public esmf2mct_init
public esmf2mct_copy

private
interface esmf2mct_init
    module procedure esmf2mct_init_AvsGSmap_Arrays
    module procedure esmf2mct_init_Av_Array
    module procedure esmf2mct_init_GGrid_Array
    module procedure esmf2mct_init_GSmap_Gindex
    module procedure esmf2mct_init_GSmap_Distgrid
    module procedure esmf2mct_init_GSmap_Array
    module procedure esmf2mct_init_Gindex_Distgrid
    module procedure esmf2mct_init_Gindex_Array
end interface

interface esmf2mct_copy
    module procedure esmf2mct_copy_a
    module procedure esmf2mct_copy_alist
end interface

contains
!--------------------------------------------------------------------------

! initialize an AttrVect from ESMF_Array
subroutine esmf2mct_init_Av_Array(array, attrVect, rc)

    implicit none

    ! inout parameters
    type(ESMF_Array), intent(inout)         :: array
    type(mct_aVect), intent(out)            :: attrVect
    integer, intent(out), optional          :: rc

    ! local variables
    integer                                 :: localrc, dsize
    character(len=8096)                     :: mct_names
    character(len=*),parameter :: subname = 'esmf2mct_init_Av_Array'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    call ESMF_AttributeGet(array, name="mct_names", value=mct_names, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    ! construct AttrVect
    ! compute dsize for mct - PE local 
    call esmfshr_util_ArrayGetSize(array, lsize2=dsize)

    call mct_aVect_init(aV=attrVect, rList=trim(mct_names), lsize=dsize)

end subroutine esmf2mct_init_Av_Array

!--------------------------------------------------------------------------

! initialize a GGrid from ESMF_Array
subroutine esmf2mct_init_GGrid_Array(array, ggrid, rc)

    implicit none

    ! inout parameters
    type(ESMF_Array), intent(inout)         :: array
    type(mct_ggrid), intent(out)            :: ggrid
    integer, intent(out), optional          :: rc

    ! local variables
    integer                                 :: localrc, dsize, lsize
    character(len=8096)                     :: mct_names
    integer, pointer                        :: gindex(:)
    type(ESMF_distgrid)                     :: distgrid
    character(len=*),parameter :: subname = 'esmf2mct_init_GGrid_Array'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    call ESMF_AttributeGet(array, name="mct_names", value=mct_names, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    ! construct ggrid
    ! compute dsize for mct - PE local 
    call esmfshr_util_ArrayGetSize(array, lsize2=dsize)

    if (len_trim(mct_names) == 0 .and. dsize == 0) then
       return   ! empty domain
    endif

    call mct_gGrid_init( GGrid=ggrid, &
         CoordChars=trim(mct_names), lsize=dsize)
    call mct_aVect_zero(ggrid%data)

    call ESMF_ArrayGet(array, distgrid=distgrid, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    call ESMF_DistgridGet(distgrid,localDE=0,elementCount=lsize,rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    if (dsize /= lsize) then
       write(shr_log_unit,*) trim(subname),' ERROR: inconsistent sizes ',dsize,lsize
       call shr_sys_abort(trim(subname))
    endif

    allocate(gindex(lsize))
    call ESMF_DistgridGet(distgrid,localDE=0,seqIndexList=gindex,rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    call mct_gGrid_importIAttr(ggrid,'GlobGridNum',gindex,lsize)

    deallocate(gindex)

end subroutine esmf2mct_init_GGrid_Array

!--------------------------------------------------------------------------

! initialize a gsMap from an array
subroutine esmf2mct_init_GSmap_Array(array, COMPID, gsmap, mpicom, gsize, rc)

    implicit none

    ! inout parameters
    type(ESMF_Array), intent(in)            :: array
    integer, intent(in)                     :: COMPID
    type(mct_gsMap), intent(out)            :: gsmap
    integer, intent(in)                     :: mpicom
    integer, intent(in), optional           :: gsize
    integer, intent(out), optional          :: rc

    ! local variables 
    type(ESMF_DistGrid)                     :: distgrid
    integer                                 :: localrc
    character(len=*),parameter :: subname = 'esmf2mct_init_GSmap_Array'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    call ESMF_ArrayGet(array, distgrid=distgrid, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    if (present(gsize)) then
       call esmf2mct_init_GSmap_Distgrid(distgrid, COMPID, gsmap, mpicom, gsize, rc=localrc)
    else
       call esmf2mct_init_GSmap_Distgrid(distgrid, COMPID, gsmap, mpicom, rc=localrc)
    endif
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

end subroutine esmf2mct_init_GSmap_Array

!--------------------------------------------------------------------------

! initialize a gsMap from an array
subroutine esmf2mct_init_Gindex_Array(array, gindex, rc)

    implicit none

    ! inout parameters
    type(ESMF_Array), intent(in)            :: array
    integer, intent(inout), pointer         :: gindex(:)
    integer, intent(out), optional          :: rc

    ! local variables 
    type(ESMF_DistGrid)                     :: distgrid
    integer                                 :: localrc
    character(len=*),parameter :: subname = 'esmf2mct_init_Gindex_Array'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    call ESMF_ArrayGet(array, distgrid=distgrid, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    call esmf2mct_init_Gindex_Distgrid(distgrid,gindex,rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

end subroutine esmf2mct_init_Gindex_Array

!--------------------------------------------------------------------------

! initialize a gsMap from a distgrid
subroutine esmf2mct_init_GSmap_Distgrid(distgrid, COMPID, gsmap, mpicom, gsize, rc)

    implicit none

    ! inout parameters
    type(ESMF_DistGrid), intent(in)         :: distgrid
    integer, intent(in)                     :: COMPID
    type(mct_gsMap), intent(out)            :: gsmap
    integer, intent(in)                     :: mpicom
    integer, intent(in), optional           :: gsize
    integer, intent(out), optional          :: rc

    ! local variables 
    type(ESMF_VM)                           :: vm
    integer                                 :: dsize, localrc
    integer, allocatable                    :: gindex(:)
    character(len=*),parameter :: subname = 'esmf2mct_init_GSmap_Distgrid'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    ! query index list
    call ESMF_DistGridGet(distgrid, localDe=0, elementCount=dsize, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    allocate(gindex(dsize), stat=localrc)
    if(localrc /= 0) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    call ESMF_DistGridGet(distgrid, localDe=0, seqIndexList=gindex, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    ! finally, construct the gsMap
    if(present(gsize)) then
        call mct_gsMap_init(gsMap, gindex, mpicom, COMPID, dsize, gsize)
    else
        call mct_gsMap_init(gsMap, gindex, mpicom, COMPID, dsize)
    endif
    deallocate(gindex)

end subroutine esmf2mct_init_GSmap_Distgrid

!--------------------------------------------------------------------------

! initialize a gsMap from a distgrid
subroutine esmf2mct_init_Gindex_Distgrid(distgrid, gindex, rc)

    implicit none

    ! inout parameters
    type(ESMF_DistGrid), intent(in)         :: distgrid
    integer, intent(inout), pointer         :: gindex(:)
    integer, intent(out), optional          :: rc

    ! local variables 
    integer                                 :: dsize, localrc
    character(len=*),parameter :: subname = 'esmf2mct_init_Gindex_Distgrid'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    ! query index list
    call ESMF_DistGridGet(distgrid, localDe=0, elementCount=dsize, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    allocate(gindex(dsize), stat=localrc)
    if(localrc /= 0) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    call ESMF_DistGridGet(distgrid, localDe=0, seqIndexList=gindex, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

end subroutine esmf2mct_init_Gindex_Distgrid

!--------------------------------------------------------------------------

! initialize a gsMap from a gindex
subroutine esmf2mct_init_GSmap_Gindex(gindex, COMPID, gsmap, mpicom, gsize, rc)

    implicit none

    ! inout parameters
    integer, intent(in)                     :: gindex(:)
    integer, intent(in)                     :: COMPID
    type(mct_gsMap), intent(out)            :: gsmap
    integer, intent(in)                     :: mpicom
    integer, intent(in), optional           :: gsize
    integer, intent(out), optional          :: rc

    ! local variables 
    type(ESMF_VM)                           :: vm
    integer                                 :: dsize, localrc
    character(len=*),parameter :: subname = 'esmf2mct_init_GSmap_Gindex'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    ! construct the gsMap
    dsize = size(gindex, 1)
    if(present(gsize)) then
        call mct_gsMap_init(gsMap, gindex, mpicom, COMPID, dsize, gsize)
    else
        call mct_gsMap_init(gsMap, gindex, mpicom, COMPID, dsize)
    endif

end subroutine esmf2mct_init_GSmap_Gindex

!--------------------------------------------------------------------------

! copy data from ESMF_Array to attrVect
subroutine esmf2mct_copy_a(array, attrVect, rc)

    implicit none

    ! inout parameters
    type(ESMF_Array), intent(in)            :: array
    type(mct_aVect), intent(inout)          :: attrVect
    integer, intent(out), optional          :: rc

    ! local variables
    integer                                 :: a_dsize, a_nfields, localrc
    real(ESMF_KIND_R8), pointer             :: a_fptr(:,:)
    integer                                 :: av_dsize, av_nfields
    integer                                 :: i, j, a_off1, a_off2
    character(len=*),parameter :: subname = 'esmf2mct_copy_a'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    call ESMF_ArrayGet(array, localDe=0, farrayPtr=a_fptr, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

! tcraig - this call seems to change answers, why!!
!    call esmfshr_util_ArrayGetSize(array, lsize1=a_nfields, lsize2=a_dsize)
    a_dsize = ubound(a_fptr, 2)-lbound(a_fptr, 2)+1
    a_nfields = ubound(a_fptr, 1)-lbound(a_fptr, 1)+1

    av_dsize = mct_aVect_lsize(aV=attrVect) 
    av_nfields = mct_aVect_nRAttr(aV=attrVect)

    if (a_dsize == 0 .or. a_nfields == 0) then
       ! nothing to copy, skip size check too
       return
    endif

    if(a_dsize /= av_dsize .or. a_nfields /= av_nfields) then
        write(shr_log_unit,*) trim(subname),' ERROR in sizes ',a_dsize,av_dsize,a_nfields,av_nfields
        write(shr_log_unit,*) trim(subname),' ERROR av nflds = ',trim(mct_aVect_exportRlist2c(attrVect))
        call shr_sys_abort(trim(subname))
    endif

    ! attrVect%rattr = a_fptr
    a_off2 = lbound(a_fptr, 2)-1
    a_off1 = lbound(a_fptr, 1)-1

    do j = 1, a_dsize
    do i = 1, a_nfields
       attrVect%rattr(i, j) = a_fptr(a_off1+i, a_off2+j)
    enddo
    enddo

end subroutine esmf2mct_copy_a
 
!--------------------------------------------------------------------------

! copy data from ESMF_Array to attrVect
subroutine esmf2mct_copy_alist(array, attrVect, list, rc)

    ! inout parameters
    type(ESMF_Array), intent(inout)         :: array
    type(mct_aVect), intent(inout)          :: attrVect
    character(len=*), intent(in)            :: list
    integer, intent(out), optional          :: rc

    ! local variables
    type(ESMF_Array)                        :: dstArray
    type(ESMF_DistGrid)                     :: distgrid
    integer                                 :: localrc

    if(present(rc)) rc = ESMF_SUCCESS

    call ESMF_ArrayGet(array, distgrid=distgrid, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
   
    !--------------------------------------------------------------------------
    ! use the attribute vector Fortran array pointer to create a temporary Array;
    ! this avoids unnecessary memory copy and code is more efficient
    !--------------------------------------------------------------------------
    dstArray = ESMF_ArrayCreate(farrayPtr=attrVect%rattr, distgrid=distgrid, distgridToArrayMap=(/2/), &
        rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(dstArray, name="mct_names", &
        value=trim(mct_aVect_exportRList2c(attrVect)), rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    ! do a field name based copy between incongruent src and dst Array pairs
    call esmfshr_util_ArrayCopy(array, dstArray, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayDestroy(dstArray, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

end subroutine esmf2mct_copy_alist

!--------------------------------------------------------------------------

! convert ESMF_Array to Attribute Vector and Global Segmentation Map
subroutine esmf2mct_init_AvsGSmap_Arrays(importArray, exportArray, COMPID, mpicom, importAttrVect, exportAttrVect, gsmap, rc)

    implicit none

    ! inout parameters
    type(ESMF_Array), intent(inout)         :: importArray
    type(ESMF_Array), intent(inout)         :: exportArray
    integer, intent(in)                     :: COMPID
    integer, intent(in)                     :: mpicom
    type(mct_aVect), intent(out)            :: importAttrVect
    type(mct_aVect), intent(out)            :: exportAttrVect
    type(mct_gsMap), intent(out)            :: gsmap
    integer, intent(out), optional          :: rc

    ! internal variables
    integer                                 :: localrc
    type(ESMF_DistGrid)                     :: i_distgrid, e_distgrid
    character(len=*),parameter :: subname = 'esmf2mct_init_AvsGSmap_Arrays'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    ! query array for distgrid
    call ESMF_ArrayGet(importArray, distgrid=i_distgrid, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(exportArray, distgrid=e_distgrid, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

!    if(.not. ESMF_DistGridMatch(i_distgrid, e_distgrid)) then
!         call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
!    endif
    if (ESMF_DistGridMatch(i_distgrid, e_distgrid) /= ESMF_DISTGRIDMATCH_EXACT) then
         call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    endif

    ! initializations
    call esmf2mct_init_GSmap_Distgrid(i_distgrid, COMPID, gsmap, mpicom, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    call esmf2mct_init_Av_Array(importArray, importAttrVect, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    call esmf2mct_init_Av_Array(exportArray, exportAttrVect, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

end subroutine esmf2mct_init_AvsGSmap_Arrays

!--------------------------------------------------------------------------

#endif
end module
