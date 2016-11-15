module esmfshr_util_mod
#ifdef USE_ESMF_LIB

use ESMF
use shr_kind_mod, only : SHR_KIND_IN
use shr_string_mod, only : shr_string_ListGetNum
use shr_string_mod, only : shr_string_ListGetIndexF
use shr_string_mod, only : shr_string_ListGetName
use shr_sys_mod, only : shr_sys_flush, shr_sys_abort
use shr_log_mod, only : shr_log_unit, shr_log_level

implicit none

!
! Author: Fei Liu
!
! public methods from this module, another layer of utils take these 
! array related methods and put them into interface blocks. user using
! the higher level utils module see/use names such as esmf2mct_init, etc
!
public esmfshr_util_ArrayCopy
public esmfshr_util_ArraySum
public esmfshr_util_ArrayZero
public esmfshr_util_ArrayAvg
public esmfshr_util_ArrayGetIndex
public esmfshr_util_ArrayGetName
public esmfshr_util_ArrayGetSize
public esmfshr_util_ArrayPutField
public esmfshr_util_ArrayGetField
public esmfshr_util_DistgridCreate
public esmfshr_util_StateArrayDestroy
public esmfshr_util_StateADistgridDestroy
public esmfshr_util_CheckRC

private

contains

!--------------------------------------------------------------------------

! query the n-th index position of a string from colon delimited strings
! it's not the position of the rattr.
function esmfshr_util_ArrayGetIndex(array, rattr, meta_attname, rc)

    implicit none

    ! inout parameters
    type(ESMF_Array), intent(inout)         :: array
    character(len=*), intent(in)            :: rattr
    character(len=*), intent(in), optional  :: meta_attname
    integer, intent(out), optional          :: rc

    ! return
    integer                                 :: esmfshr_util_ArrayGetIndex

    ! local
    integer                                 :: attr_len, shift, localrc
    character(len=9)                        :: meta_str = 'mct_names'
    character(len=8196)                     :: rattrs
    logical                                 :: found, last
    character(len=*),parameter :: subname = 'esmfshr_util_ArrayGetIndex'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    esmfshr_util_ArrayGetIndex = 0
    found = .false.
    last = .false.

    attr_len = len_trim(rattr)

    if(present(meta_attname)) then
        call ESMF_AttributeGet(array, name=meta_attname, value=rattrs, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(array, name=meta_str, value=rattrs, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    endif

    esmfshr_util_ArrayGetIndex = shr_string_listGetIndexF(trim(rattrs),trim(rattr))

end function esmfshr_util_ArrayGetIndex

!--------------------------------------------------------------------------

! query the n-th index position of a string from colon delimited strings
! it's not the position of the rattr.
subroutine esmfshr_util_ArrayGetName(array, kfld, name, meta_attname, rc)

    implicit none

    ! inout parameters
    type(ESMF_Array), intent(inout)         :: array
    integer, intent(in)                     :: kfld
    character(len=*), intent(out)           :: name
    character(len=*), intent(in), optional  :: meta_attname
    integer, intent(out), optional          :: rc

    ! local
    integer                                 :: localrc
    character(len=9)                        :: meta_str = 'mct_names'
    character(len=8196)                     :: rattrs
    character(len=*),parameter :: subname = 'esmfshr_util_ArrayGetName'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    name = ''

    if(present(meta_attname)) then
        call ESMF_AttributeGet(array, name=meta_attname, value=rattrs, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(array, name=meta_str, value=rattrs, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    endif

    call shr_string_listGetName(trim(rattrs),kfld,name)

end subroutine esmfshr_util_ArrayGetName

!--------------------------------------------------------------------------

! query the n-th index position of a string from colon delimited strings
! it's not the position of the rattr.
subroutine esmfshr_util_ArrayGetSize(array, lsize1, lsize2, rc)

    implicit none

    ! inout parameters
    type(ESMF_Array), intent(in)            :: array
    integer, intent(out),optional           :: lsize1
    integer, intent(out),optional           :: lsize2
    integer, intent(out),optional           :: rc

    ! local
    integer  :: s1, s2, lrc
    real(ESMF_KIND_R8), pointer :: fptr(:, :)
    character(len=*),parameter :: subname = 'esmfshr_util_ArrayGetSize'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    call ESMF_ArrayGet(array, localDe=0, farrayPtr = fptr, rc=lrc)
    if(lrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=lrc, endflag=ESMF_END_ABORT)
!    call ESMF_ArrayGet(array, dim=1, localDe=0, indexCount=s1)
!    call ESMF_ArrayGet(array, dim=2, localDe=0, indexCount=s2)

    s1 = ubound(fptr,1) - lbound(fptr,1) + 1
    s2 = ubound(fptr,2) - lbound(fptr,2) + 1

    if (present(lsize1)) then
       lsize1 = s1
    endif

    if (present(lsize2)) then
       lsize2 = s2
    endif

end subroutine esmfshr_util_ArrayGetSize

!--------------------------------------------------------------------------

! put a field into the array
subroutine esmfshr_util_ArrayPutField(array, fname, buf, rc)

    implicit none

    ! inout parameters
    type(ESMF_Array), intent(inout)   :: array
    character(len=*), intent(in)      :: fname
    real*8, pointer, intent(in)       :: buf(:)
    integer, intent(out), optional    :: rc

    ! local
    integer                           :: localrc,sa,sb,lba2,uba2
    integer                           :: nfld
    real(ESMF_KIND_R8), pointer       :: fptr(:,:)
    character(len=*),parameter :: subname = 'esmfshr_util_ArrayPutField'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    sa   = -999
    sb   = -999
    lba2 = -999
    uba2 = -999
    nfld = -999

    nfld = esmfshr_util_ArrayGetIndex(array,fname)

    if (nfld < 1) then
       write(shr_log_unit,*) subname,' ERROR: fld not found on array ',trim(fname)
!!       call shr_sys_abort(subname)
    endif

    sb = size(buf)
    call esmfshr_util_ArrayGetSize(array, lsize2=sa)

    call ESMF_ArrayGet(array, localDe=0, farrayPtr=fptr, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    lba2 = lbound(fptr,2)
    uba2 = ubound(fptr,2)

    if (sb /= sa .or. sa /= uba2-lba2+1 ) then
       write(shr_log_unit,*) subname,' ERROR: sb /= sa',sb,sa,trim(fname),lba2,uba2
       call shr_sys_abort(subname)
    endif

    fptr(nfld,lba2:uba2) = buf(1:sa)

end subroutine esmfshr_util_ArrayPutField
!--------------------------------------------------------------------------

! put a field into the array
subroutine esmfshr_util_ArrayGetField(array, fname, buf, rc)

    implicit none

    ! inout parameters
    type(ESMF_Array), intent(inout)   :: array
    character(len=*), intent(in)      :: fname
    real*8, intent(inout)             :: buf(:)
    integer, intent(out), optional    :: rc

    ! local
    integer                           :: localrc,sa,sb,lba2,uba2
    integer                           :: nfld
    real(ESMF_KIND_R8), pointer       :: fptr(:,:)
    character(len=*),parameter :: subname = 'esmfshr_util_ArrayGetField'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    sa   = -999
    sb   = -999
    lba2 = -999
    uba2 = -999
    nfld = -999

    nfld = esmfshr_util_ArrayGetIndex(array,fname)
    if (nfld < 1) then
       write(shr_log_unit,*) subname,' ERROR: fname not on array ',trim(fname)
!!       call shr_sys_abort(subname)
    endif

    sb = size(buf)
    call esmfshr_util_ArrayGetSize(array, lsize2=sa)

    call ESMF_ArrayGet(array, localDe=0, farrayPtr=fptr, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    lba2 = lbound(fptr,2)
    uba2 = ubound(fptr,2)

    if (sa /= sb .or. sa /= uba2-lba2+1) then
       write(shr_log_unit,*) subname,' ERROR: sa /= sb',sb,sa,trim(fname),lba2,uba2
       call shr_sys_abort(subname)
    endif

    buf = 0.0
    buf(1:sa) = fptr(nfld,lba2:uba2)

end subroutine esmfshr_util_ArrayGetField
!--------------------------------------------------------------------------

! put a field into the array
subroutine esmfshr_util_StateArrayDestroy(state, array_name, rc)

    implicit none

    ! inout parameters
    type(ESMF_State), intent(inout)   :: state
    character(len=*), intent(in)      :: array_name
    integer, intent(out), optional    :: rc

    ! local
    type(ESMF_Array)    :: array
    type(ESMF_DistGrid) :: distgrid
    integer             :: lrc
    character(len=*),parameter :: subname = 'esmfshr_util_StateArrayDestroy'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemName=trim(array_name), array=array, rc=lrc)
    call esmfshr_util_CheckRC(lrc,trim(subname)//':StateGet array')

    call ESMF_ArrayGet(array, distgrid=distgrid, rc=lrc)
    call esmfshr_util_CheckRC(lrc,trim(subname)//':ArrayGet distgrid')

! not safe in esmf
!    call ESMF_DistgridDestroy(distgrid, rc=lrc)
!    call esmfshr_util_CheckRC(lrc,trim(subname)//':DistgridDestroy distgrid')

    call ESMF_ArrayDestroy(array, rc=lrc)
    call esmfshr_util_CheckRC(lrc,trim(subname)//':ArrayDestroy array')

end subroutine esmfshr_util_StateArrayDestroy

!--------------------------------------------------------------------------

! initialize a distgrid
subroutine esmfshr_util_DistgridCreate(distgrid, gsize,nxg,nyg,mpicom,decomp,rc)

    implicit none

    ! inout parameters
    type(ESMF_distgrid), intent(inout) :: distgrid
    integer(SHR_KIND_IN),intent(out)   :: gsize
    integer(SHR_KIND_IN),intent(in)    :: nxg
    integer(SHR_KIND_IN),intent(in)    :: nyg
    integer(SHR_KIND_IN),intent(in)    :: mpicom
    character(len=*), intent(in)       :: decomp
    integer, intent(out), optional     :: rc

    ! local

    integer(SHR_KIND_IN) :: k,n,npes,ierr,i, my_task
    integer(SHR_KIND_IN), pointer :: start(:)     ! for gsmap initialization
    integer(SHR_KIND_IN), pointer :: length(:)    ! for gsmap initialization
    integer(SHR_KIND_IN), pointer :: pe_loc(:)    ! for gsmap initialization
    integer(SHR_KIND_IN), pointer :: arbIndex(:)  ! for gsmap initialization
    integer(SHR_KIND_IN), pointer :: deBlockList(:,:,:)  ! for gsmap initialization
    character(*), parameter :: subname = '(distgrid_create) '
    character(*), parameter :: F00   = "('(distgrid_create) ',8a)"
    character(*), parameter :: F01   = "('(distgrid_create) ',a,7i8)"
    character(*), parameter :: F02   = "('(distgrid_create) ',a,4i8)"

    ! ---------------------------------------------
    rc = ESMF_SUCCESS

    call mpi_comm_size(mpicom,npes,ierr)
    call mpi_comm_rank(mpicom,my_task,ierr)
    allocate(start(npes),length(npes),pe_loc(npes))

    start = 0
    length = 0
    gsize = nxg*nyg
    if (gsize > 0) then
       allocate(deBlockList(1,2,npes))
       do n = 1,npes
          if (trim(decomp) == '1d') then
             length(n)  = gsize/npes
             if (n <= mod(gsize,npes)) length(n) = length(n) + 1
          elseif (trim(decomp) == 'root') then
             length = 0
             length(1) = gsize
          else
             call shr_sys_abort(subname//' ERROR decomp not allowed')
          endif
          if (n == 1) then
             start(n) = 1
          else
             start(n) = start(n-1) + length(n-1)
          endif
          pe_loc(n) = n-1
!debug          if (my_task == 0) then
!             write(shr_log_unit,F01) 'gsmap: ',my_task,k,n,start(n),length(n),pe_loc(n), gsize
!             call shr_sys_flush(shr_log_unit)
!          endif
          deBlockList(:,1,n) = (/start(n)/)
          deBlockList(:,2,n) = (/start(n)+length(n)-1/)
       enddo
       distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/gsize/), deBlockList=deBlockList, rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       deallocate(deBlockList)
    else
       distgrid = ESMF_DistGridCreate(minIndex=(/0/), maxIndex=(/-1/), rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    endif

    deallocate(start,length,pe_loc)

end subroutine esmfshr_util_DistgridCreate

!--------------------------------------------------------------------------

! put a field into the array
subroutine esmfshr_util_StateADistgridDestroy(state, array_name, rc)

    implicit none

    ! inout parameters
    type(ESMF_State), intent(inout)   :: state
    character(len=*), intent(in)      :: array_name
    integer, intent(out), optional    :: rc

    ! local
    type(ESMF_Array)    :: array
    type(ESMF_DistGrid) :: distgrid
    integer             :: lrc
    character(len=*),parameter :: subname = 'esmfshr_util_StateADistgridDestroy'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemName=trim(array_name), array=array, rc=lrc)
    call esmfshr_util_CheckRC(lrc,trim(subname)//':StateGet array')

    call ESMF_ArrayGet(array, distgrid=distgrid, rc=lrc)
    call esmfshr_util_CheckRC(lrc,trim(subname)//':ArrayGet distgrid')

    call ESMF_DistgridDestroy(distgrid, rc=lrc)
    call esmfshr_util_CheckRC(lrc,trim(subname)//':DistgridDestroy distgrid')

end subroutine esmfshr_util_StateADistgridDestroy
!--------------------------------------------------------------------------

! put a field into the array
subroutine esmfshr_util_CheckRC(rc,string,abort)

    implicit none

    ! inout parameters
    integer         , intent(in)           :: rc
    character(len=*), intent(in), optional :: string
    logical         , intent(in), optional :: abort

    ! local
    character(len=1024) :: lstring
    integer             :: lrc
    logical             :: labort
    character(len=*),parameter :: subname = 'esmfshr_util_CheckRC'
    !---------------------------------------

    lstring = ''
    if (present(string)) then
       lstring = trim(string)
    endif

    labort = .true.
    if (present(abort)) then
       labort = abort
    endif
    
    if (rc /= ESMF_SUCCESS) then
       write(shr_log_unit,*) trim(subname),' ERROR: ',trim(lstring)
       if (labort) then
          call shr_sys_abort(trim(subname)//' ERROR: '//trim(lstring))
          call ESMF_Finalize(rc=lrc, endflag=ESMF_END_ABORT)
       endif
    endif

end subroutine esmfshr_util_CheckRC

!--------------------------------------------------------------------------

! copy matching data from src Array to dst Array
subroutine esmfshr_util_ArrayCopy(srcArray, dstArray, rc)

    implicit none

    ! arguments
    type(ESMF_Array), intent(inout)         :: srcArray
    type(ESMF_Array), intent(inout)         :: dstArray
    integer, intent(out),optional           :: rc

    ! local
    real(ESMF_KIND_R8), pointer             :: srcFptr(:,:), dstFptr(:,:)
    integer                                 :: srcLB(2), dstLB(2), srcUB(2), dstUB(2)
    integer                                 :: i, j, localrc
    integer                                 :: src_nfields, dst_nfields, im
    integer, allocatable                    :: src2dst_map(:)
    integer, allocatable                    :: s1(:),s2(:),d1(:),d2(:)
    integer, allocatable                    :: srcidx(:), dstidx(:)
    integer                                 :: ssize,dsize, map_nfields
    character(len=8196)                     :: src_attname, dst_attname
    character(len=128)                      :: src_att, dst_att
    character(len=*),parameter :: subname = 'esmfshr_util_ArrayCopy'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    call ESMF_ArrayGet(srcArray, localDe=0, farrayPtr=srcFptr, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    call ESMF_ArrayGet(dstArray, localDe=0, farrayPtr=dstFptr, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(srcArray, name="mct_names", value=src_attname, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeGet(dstArray, name="mct_names", value=dst_attname, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    src_nfields = shr_string_ListGetNum(src_attname)
    dst_nfields = shr_string_ListGetNum(dst_attname)

    if (src_nfields < 1 .or. dst_nfields < 1) then
       return
    endif
    
    allocate(src2dst_map(src_nfields), stat=localrc)
    allocate(s1(src_nfields),s2(src_nfields),stat=localrc)
    allocate(d1(dst_nfields),d2(dst_nfields),stat=localrc)
    if(localrc /= 0) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    ssize = len_trim(src_attname)
    s1(1) = 1
    s2(src_nfields) = ssize
    do i = 2,src_nfields
       s1(i) = s1(i-1) + index(src_attname(s1(i-1):ssize),':')
       s2(i-1) = s1(i) - 2
    enddo

    dsize = len_trim(dst_attname)
    d1(1) = 1
    d2(dst_nfields) = dsize
    do i = 2,dst_nfields
       d1(i) = d1(i-1) + index(dst_attname(d1(i-1):dsize),':')
       d2(i-1) = d1(i) - 2
    enddo

    src2dst_map = -1
    ! build the map in O(n*n) time
    do i = 1, src_nfields
        src_att = src_attname(s1(i):s2(i))
        do j = 1, dst_nfields
            dst_att = dst_attname(d1(j):d2(j))
            if(trim(src_att) == trim(dst_att)) then
                if(src2dst_map(i) == -1) then
                    src2dst_map(i) = j
                    exit    ! early exit 
                else
                    ! cannot have a single key maps to multiple values, abort
                    call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
                endif
            endif
        enddo
    enddo

    deallocate(s1,s2,d1,d2)

    srcLB = lbound(srcFptr)
    dstLB = lbound(dstFptr)
    srcUB = ubound(srcFptr)
    dstUB = ubound(dstFptr)
    if(srcLB(2) /= dstLB(2) .or. srcUB(2) /= dstUB(2)) then
        ! src and dst Array must have identical shape, otherwise abort
        call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    endif
        
    ! construct index arrays to speed up copy in linear time
    map_nfields = 0
    do i = srcLB(1), srcUB(1)
        if(src2dst_map(i) /= -1) map_nfields = map_nfields + 1
    enddo
    allocate(srcidx(map_nfields), dstidx(map_nfields), stat=localrc)
    if(localrc /= 0) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    j = 1
    do i = srcLB(1), srcUB(1)
        if(src2dst_map(i) /= -1) then
            srcidx(j) = i
            dstidx(j) = src2dst_map(i-srcLB(1)+1) + dstLB(1) - 1
            j = j + 1
        endif
    enddo
        
    ! do data copy, this is faster than reverse loop order and block copy
    do j = srcLB(2), srcUB(2)
        do i = 1, map_nfields
               dstFptr(dstidx(i), j) = srcFptr(srcidx(i), j)
        enddo
    enddo

    deallocate(src2dst_map, srcidx, dstidx)

end subroutine esmfshr_util_ArrayCopy

!--------------------------------------------------------------------------

! add matching data values from src Array to dst Array
subroutine esmfshr_util_ArraySum(srcArray, dstArray, rc)

    implicit none

    ! arguments
    type(ESMF_Array), intent(inout)         :: srcArray
    type(ESMF_Array), intent(inout)         :: dstArray
    integer, intent(out),optional           :: rc
    ! local
    real(ESMF_KIND_R8), pointer             :: srcFptr(:,:), dstFptr(:,:)
    integer                                 :: srcLB(2), dstLB(2), srcUB(2), dstUB(2)
    integer                                 :: i, j, localrc
    character(len=*),parameter :: subname = 'esmfshr_util_ArraySum'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    call ESMF_ArrayGet(srcArray, localDe=0, farrayPtr=srcFptr, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    call ESMF_ArrayGet(dstArray, localDe=0, farrayPtr=dstFptr, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    srcLB = lbound(srcFptr)
    dstLB = lbound(dstFptr)
    srcUB = ubound(srcFptr)
    dstUB = ubound(dstFptr)
    if(srcLB(1) /= dstLB(1) .or. srcLB(2) /= dstLB(2) .or. &
       srcUB(1) /= dstUB(1) .or. srcUB(2) /= dstUB(2)) then
        call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    endif
        
    do j = srcLB(2), srcUB(2)
        do i = srcLB(1), srcUB(1)
            dstFptr(i, j) = dstFptr(i, j) + srcFptr(i, j)
        enddo
    enddo

end subroutine esmfshr_util_ArraySum

!--------------------------------------------------------------------------

! zero the values of a esmf Array
subroutine esmfshr_util_ArrayZero(array, rc)

    implicit none

    ! arguments
    type(ESMF_Array), intent(inout)         :: array
    integer, intent(out),optional           :: rc
    ! local
    real(ESMF_KIND_R8), pointer             :: fptr(:, :)
    integer                                 :: localrc
    character(len=*),parameter :: subname = 'esmfshr_util_ArrayZero'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    call ESMF_ArrayGet(array, localDe=0, farrayPtr=fptr, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    fptr(:, :) = 0.0_ESMF_KIND_R8

end subroutine esmfshr_util_ArrayZero

!--------------------------------------------------------------------------

! Average the values of a esmf Array
subroutine esmfshr_util_ArrayAvg(array, count, rc)

    implicit none

    ! arguments
    type(ESMF_Array), intent(inout)         :: array
    integer, intent(in)                     :: count
    integer, intent(out),optional           :: rc
    ! local
    real(ESMF_KIND_R8), pointer             :: fptr(:, :)
    integer                                 :: localrc
    character(len=*),parameter :: subname = 'esmfshr_util_ArrayAvg'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    call ESMF_ArrayGet(array, localDe=0, farrayPtr=fptr, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

    fptr(:, :) = fptr(:, :)/count

end subroutine esmfshr_util_ArrayAvg

#endif
end module
