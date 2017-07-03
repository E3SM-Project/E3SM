module Lookup_Table_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  type, abstract, public :: lookup_table_base_type
    PetscInt :: dim
    PetscInt :: dims(3)
    PetscReal, pointer :: data(:)    
    class(lookup_table_axis_type), pointer :: axis1
  contains
    procedure(LookupTableEvaluateDummy), deferred, public :: Sample 
  end type lookup_table_base_type
  
  type, public, extends(lookup_table_base_type) :: lookup_table_uniform_type
    class(lookup_table_axis_type), pointer :: axis2
    class(lookup_table_axis_type), pointer :: axis3
  contains
    procedure, public :: Sample => LookupTableEvaluateUniform
  end type lookup_table_uniform_type
  
  type, public, extends(lookup_table_base_type) :: lookup_table_general_type
    class(lookup_table_axis2_general_type), pointer :: axis2
  contains
    procedure, public :: Sample => LookupTableEvaluateGeneral
  end type lookup_table_general_type
  
  type, public :: lookup_table_axis_type
    PetscInt :: itype
    PetscInt :: saved_index
    PetscReal, pointer :: values(:)
  end type lookup_table_axis_type
  
  type, public, extends(lookup_table_axis_type) :: lookup_table_axis2_general_type
    PetscInt :: saved_index2
  end type lookup_table_axis2_general_type
  
  abstract interface
    function LookupTableEvaluateDummy(this,lookup1,lookup2,lookup3)
      import lookup_table_base_type
      implicit none
      class(lookup_table_base_type) :: this
      PetscReal :: lookup1
      PetscReal, optional :: lookup2
      PetscReal, optional :: lookup3
      PetscReal :: LookupTableEvaluateDummy
    end function LookupTableEvaluateDummy
  end interface
  
  interface LookupTableTest
    module procedure LookupTableTest1D
    module procedure LookupTableTest2D
  end interface

  interface LookupTableDestroy
    module procedure LookupTableUniformDestroy
    module procedure LookupTableGeneralDestroy
  end interface

  public :: LookupTableCreateUniform, &
            LookupTableCreateGeneral, &
            LookupTableDestroy, &
            LookupTableTest
  
contains

! ************************************************************************** !

subroutine LookupTableBaseInit(lookup_table)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 

  implicit none
  
  class(lookup_table_base_type) :: lookup_table
  
  lookup_table%dim = 0
  lookup_table%dims = 0
  nullify(lookup_table%data)
  nullify(lookup_table%axis1)

end subroutine LookupTableBaseInit

! ************************************************************************** !

function LookupTableCreateUniform(dim)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 

  implicit none
  
  PetscInt :: dim
  
  class(lookup_table_uniform_type), pointer :: LookupTableCreateUniform

  class(lookup_table_uniform_type), pointer :: lookup_table
  
  allocate(lookup_table)
  call LookupTableBaseInit(lookup_table)
  lookup_table%dim = dim
  nullify(lookup_table%axis2)
  nullify(lookup_table%axis3)
  allocate(lookup_table%axis1)
  call LookupTableAxisInit(lookup_table%axis1)
  if (dim > 1) then
    allocate(lookup_table%axis2)
    call LookupTableAxisInit(lookup_table%axis2)
  endif
  if (dim > 2) then
    allocate(lookup_table%axis3)
    call LookupTableAxisInit(lookup_table%axis3)
  endif
  
  LookupTableCreateUniform => lookup_table

end function LookupTableCreateUniform

! ************************************************************************** !

function LookupTableCreateGeneral(dim)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 

  implicit none
  
  PetscInt :: dim
  
  class(lookup_table_general_type), pointer :: LookupTableCreateGeneral

  class(lookup_table_general_type), pointer :: lookup_table
  
  allocate(lookup_table)
  call LookupTableBaseInit(lookup_table)
  lookup_table%dim = dim
  nullify(lookup_table%axis2)
  allocate(lookup_table%axis1)
  call LookupTableAxisInit(lookup_table%axis1)
  if (dim > 1) then
    allocate(lookup_table%axis2)
    call LookupTableAxisInit(lookup_table%axis2)
    lookup_table%axis2%saved_index2 = 1
  endif
  
  LookupTableCreateGeneral => lookup_table

end function LookupTableCreateGeneral

! ************************************************************************** !

subroutine LookupTableAxisInit(axis)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 

  implicit none
  
  class(lookup_table_axis_type) :: axis
  
  axis%itype = 0
  axis%saved_index = 1
  nullify(axis%values)
  
end subroutine LookupTableAxisInit

! ************************************************************************** !

function LookupTableEvaluateUniform(this,lookup1,lookup2,lookup3)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  implicit none
  
  class(lookup_table_uniform_type) :: this
  PetscReal :: lookup1
  PetscReal, optional :: lookup2
  PetscReal, optional :: lookup3 
  
  PetscReal :: LookupTableEvaluateUniform

  call LookupTableIndexUniform(this,lookup1,lookup2,lookup3)
  if (present(lookup3)) then
!    call LookupTableInterpolate3DUniform(this,lookup1,lookup2,lookup3,LookupTableEvaluateUniform)
  else if (present(lookup2)) then
    call LookupTableInterpolate2DUniform(this,lookup1,lookup2,LookupTableEvaluateUniform)
  else
    call LookupTableInterpolate1D(this,lookup1,LookupTableEvaluateUniform)
  endif
  
end function LookupTableEvaluateUniform


! ************************************************************************** !

function LookupTableEvaluateGeneral(this,lookup1,lookup2,lookup3)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  implicit none
  
  class(lookup_table_general_type) :: this
  PetscReal :: lookup1
  PetscReal, optional :: lookup2
  PetscReal, optional :: lookup3 
  
  PetscReal :: LookupTableEvaluateGeneral

  call LookupTableIndexGeneral(this,lookup1,lookup2,lookup3)
  if (present(lookup3)) then
!    call LookupTableInterpolate3DGeneral(this,lookup1,lookup2,lookup3,LookupTableEvaluateGeneral)
  else if (present(lookup2)) then
    call LookupTableInterpolate2DGeneral(this,lookup1,lookup2,LookupTableEvaluateGeneral)
  else
    call LookupTableInterpolate1D(this,lookup1,LookupTableEvaluateGeneral)
  endif
  
end function LookupTableEvaluateGeneral

! ************************************************************************** !

subroutine LookupTableIndexUniform(this,lookup1,lookup2,lookup3)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  implicit none
  
  class(lookup_table_uniform_type) :: this
  PetscReal :: lookup1
  PetscReal, optional :: lookup2
  PetscReal, optional :: lookup3 
  
  call LookupTableAxisIndexUniform(this%axis1,lookup1)
  if (associated(this%axis2)) then
    call LookupTableAxisIndexUniform(this%axis2,lookup2)
  endif
  if (associated(this%axis3)) then
    call LookupTableAxisIndexUniform(this%axis3,lookup3)
  endif

end subroutine LookupTableIndexUniform

! ************************************************************************** !

subroutine LookupTableIndexGeneral(this,lookup1,lookup2,lookup3)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  implicit none
  
  class(lookup_table_general_type) :: this
  PetscReal :: lookup1
  PetscReal, optional :: lookup2
  PetscReal, optional :: lookup3 
  
  PetscInt :: ja, jb
  PetscInt :: istart, iend
  class(lookup_table_axis2_general_type), pointer :: axis2
  
  ! axis 1 corresponds to the j dim below
  call LookupTableAxisIndexGeneral(lookup1,this%axis1%values, &
                                   this%axis1%saved_index)
  
  if (.not. associated(this%axis2)) return
  
  axis2 => this%axis2
 
  ja = this%axis1%saved_index
  if (ja > 0) then
    jb = max(min(ja+1,this%dims(1)),1)
  else
    ja = 1
    jb = 1
  endif
  
  iend = ja*this%dims(2)
  istart = iend - this%dims(2) + 1
  call LookupTableAxisIndexGeneral(lookup2,axis2%values(istart:iend), &
                                   axis2%saved_index)
  if (ja /= jb) then
    iend = jb*this%dims(2)
    istart = iend - this%dims(2) + 1
    call LookupTableAxisIndexGeneral(lookup2,axis2%values(istart:iend), &
                                     axis2%saved_index2)
  else 
    axis2%saved_index2 = axis2%saved_index
  endif

end subroutine LookupTableIndexGeneral

! ************************************************************************** !

subroutine LookupTableAxisIndexUniform(this,lookup1)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  implicit none
  
  class(lookup_table_axis_type) :: this
  PetscReal :: lookup1
  
  PetscInt :: i1
  PetscInt :: size1
  PetscReal :: begin1
  
  size1 = size(this%values)
  begin1 = this%values(1)
  i1 = int((lookup1 - begin1) / (this%values(size1) - begin1) * (size1-1) + 1)
  ! truncate i1 to zero indicating the value is below the range specified
  i1 = max(min(i1,size1),0)
  this%saved_index = i1

end subroutine LookupTableAxisIndexUniform

! ************************************************************************** !

subroutine LookupTableAxisIndexGeneral(lookup1,values,saved_index)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  implicit none
  
  PetscReal :: lookup1
  PetscInt :: saved_index
  PetscReal :: values(:)
  
  PetscInt :: i1
  PetscInt :: j1
  PetscInt :: mid1
  PetscInt :: sizei
 
  sizei = size(values)
  if (lookup1 < values(1)) then
    saved_index = 0
    return
  else if (lookup1 > values(sizei)) then
    saved_index = sizei
    return
  endif
  i1 = max(min(saved_index,sizei-1),1)
  if (lookup1 > values(i1+1) .or. &
      lookup1 < values(i1)) then
    ! move either up or down array
    if (lookup1 > values(i1+1)) then
      i1 = i1+1
      j1 = sizei
    else 
      j1 = i1
      i1 = 1
    endif
    do
      mid1 = (j1+i1) / 2
      if (lookup1 > values(mid1)) then
        i1 = mid1
      else
        j1 = mid1
      endif
      if (j1-i1 <= 1) exit
    enddo
  endif
  saved_index = i1

end subroutine LookupTableAxisIndexGeneral

! ************************************************************************** !

subroutine LookupTableInterpolate1D(this,lookup1,result)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  use Utility_module, only : Interpolate
  
  implicit none
  
  class(lookup_table_base_type) :: this
  PetscReal :: lookup1
  PetscReal :: result
  
  PetscInt :: i1, j1
  
  i1 = this%axis1%saved_index
  if (i1 > 0) then
    j1 = max(min(i1+1,this%dims(1)),1)
    call Interpolate(this%axis1%values(j1),this%axis1%values(i1),lookup1, &
                     this%data(j1),this%data(i1),result)
  else ! catch values below axis range
    result = this%data(1)
  endif
  
end subroutine LookupTableInterpolate1D

! ************************************************************************** !

subroutine LookupTableInterpolate2DUniform(this,lookup1,lookup2,result)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  use Utility_module, only : Interpolate, InterpolateBilinear

  implicit none
  
  class(lookup_table_uniform_type) :: this
  PetscReal :: lookup1
  PetscReal :: lookup2
  PetscReal :: result
  
  PetscInt :: i1, i2, j1, j2
  PetscReal :: x1, x2, y1, y2, z1, z2, z3, z4
  PetscInt :: sizei, sizej

  !  x1,y2,z3 ------ x2,y2,z4
  !     |               |
  !     |               |
  !     |   x,y         |
  !     |               |
  !  x1,y1,z1 ------ x2,y1,z2                   
  
  result = UNINITIALIZED_DOUBLE
  sizei = this%dims(1)
  sizej = this%dims(2)
  i1 = this%axis1%saved_index
  j1 = this%axis2%saved_index
  ! index axes
  if (i1 > 0) then
    i2 = max(min(i1+1,sizei),1)
  else
    i1 = 1
    i2 = 1
  endif
  if (j1 > 0) then
    j2 = max(min(j1+1,sizej),1)
  else
    j1 = 1
    j2 = 1
  endif
  if (i2 == i1) then
    if (j2 == j1) then
      ! corner of domain
      result = this%data(i1+(j1-1)*sizei)
    else
      y1 = this%axis2%values(j1)
      y2 = this%axis2%values(j2)
      z1 = this%data(i1+(j1-1)*sizei)
      z3 = this%data(i2+(j2-1)*sizei)
      call Interpolate(y2,y1,lookup2,z3,z1,result)
    endif
  else if (j2 == j1) then
    x1 = this%axis1%values(i1)
    x2 = this%axis1%values(i2)
    z1 = this%data(i1+(j1-1)*sizei)
    z2 = this%data(i2+(j2-1)*sizei)
    call Interpolate(x2,x1,lookup1,z2,z1,result)
  else
    x1 = this%axis1%values(i1)
    x2 = this%axis1%values(i2)
    y1 = this%axis2%values(j1)
    y2 = this%axis2%values(j2)
    z1 = this%data(i1+(j1-1)*sizei)
    z2 = this%data(i2+(j1-1)*sizei)
    z3 = this%data(i1+(j2-1)*sizei)
    z4 = this%data(i2+(j2-1)*sizei)
    result = InterpolateBilinear(lookup1,lookup2,x1,x2,y1,y2,z1,z2,z3,z4)
  endif
                   
end subroutine LookupTableInterpolate2DUniform

! ************************************************************************** !

subroutine LookupTableInterpolate2DGeneral(this,lookup1,lookup2,result)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  use Utility_module, only : Interpolate

  implicit none
  
  class(lookup_table_general_type) :: this
  PetscInt :: i1a
  PetscInt :: i1b
  PetscInt :: ja
  PetscReal :: lookup1
  PetscReal :: lookup2
  PetscReal :: result
  
  PetscInt :: i2a, i2b
  PetscInt :: jb
  PetscInt :: ii1, ii2, iia, iib
  PetscReal :: x1, x2, z1, z2, xa, xb
  PetscReal :: interp_a, interp_b
  PetscInt :: sizei, sizej

  
  !         x2,y2,z4
  !           /|
  !          / |
  !         /  |
  !        /   |
  !       /    |
  !      /     |
  !  x1,y2,z3  |
  !     |      |
  !     |  x,y |
  !     |      |
  !  x1,y1,z1 -x2,y1,z2                   
 
  result = UNINITIALIZED_DOUBLE
  sizei = this%dims(2)
  sizej = this%dims(1)
  ! index axes
  ja = this%axis1%saved_index
  i1a = this%axis2%saved_index
  i1b = this%axis2%saved_index2
  if (ja > 0) then
    jb = max(min(ja+1,sizej),1)
  else
    ja = 1
    jb = 1
  endif
  if (i1a > 0) then
    i2a = max(min(i1a+1,sizei),1)
  else
    i1a = 1
    i2a = 1
  endif
  if (i1b > 0) then
    i2b = max(min(i1b+1,sizei),1)
  else
    i1b = 1
    i2b = 1
  endif
  if (jb == ja) then
    ! only use ja/i1a/i2a
    if (i2a == i1a) then
      ! corner of domain
      result = this%data(i1a+(ja-1)*sizei)
    else
      ii1 = i1a+(ja-1)*sizei
      ii2 = i2a+(ja-1)*sizei
      x1 = this%axis2%values(ii1)
      x2 = this%axis2%values(ii2)
      z1 = this%data(ii1)
      z2 = this%data(ii2)
      call Interpolate(x2,x1,lookup2,z2,z1,result)
    endif
  else
    ! ja / i*a
    if (i2a == i1a) then
      interp_a = this%data(i1a+(ja-1)*sizei)
    else
      iia = i1a+(ja-1)*sizei
      iib = i2a+(ja-1)*sizei
      x1 = this%axis2%values(iia)
      x2 = this%axis2%values(iib)
      z1 = this%data(iia)
      z2 = this%data(iib)
      call Interpolate(x2,x1,lookup2,z2,z1,interp_a)
    endif
    ! jb / i*b
    if (i2b == i1b) then
      interp_b = this%data(i1b+(jb-1)*sizei)
    else
      iia = i1b+(jb-1)*sizei
      iib = i2b+(jb-1)*sizei
      x1 = this%axis2%values(iia)
      x2 = this%axis2%values(iib)
      z1 = this%data(iia)
      z2 = this%data(iib)
      call Interpolate(x2,x1,lookup2,z2,z1,interp_b)
    endif
    xa = this%axis1%values(ja)
    xb = this%axis1%values(jb)
    call Interpolate(xb,xa,lookup1,interp_b,interp_a,result)
  endif
                   
end subroutine LookupTableInterpolate2DGeneral

! ************************************************************************** !

subroutine LookupTableTest1D(lookup_table,lookup,desired_result)
  ! 
  ! Deallocates any allocated pointers in axis
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  use Utility_module
  
  implicit none
  
  class(lookup_table_base_type) :: lookup_table
  PetscReal :: lookup
  PetscReal :: desired_result
  
  PetscReal :: result
  
100 format(3(f12.6),l)

  result = lookup_table%Sample(lookup)
  write(*,100) lookup, desired_result, result, Equal(result,desired_result)

end subroutine LookupTableTest1D

! ************************************************************************** !

subroutine LookupTableTest2D(lookup_table,lookup1,lookup2,desired_result)
  ! 
  ! Deallocates any allocated pointers in axis
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  use Utility_module
  
  implicit none
  
  class(lookup_table_base_type) :: lookup_table
  PetscReal :: lookup1
  PetscReal :: lookup2
  PetscReal :: desired_result
  
  PetscReal :: result
  
100 format(4(f12.6),l)

  result = lookup_table%Sample(lookup1,lookup2)
  write(*,100) lookup1, lookup2, desired_result, result, &
    Equal(result,desired_result)

end subroutine LookupTableTest2D

! ************************************************************************** !

subroutine LookupTableAxisDestroy(axis)
  ! 
  ! Deallocates any allocated pointers in axis
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  use Utility_module
  
  implicit none
  
  class(lookup_table_axis_type), pointer :: axis
  
  if (.not.associated(axis)) return
  
  call DeallocateArray(axis%values)
  deallocate(axis)
  nullify(axis)

end subroutine LookupTableAxisDestroy

! ************************************************************************** !

subroutine LookupTableUniformDestroy(lookup_table)
  ! 
  ! Deallocates any allocated pointers in lookup table
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  use Utility_module
  
  implicit none
  
  class(lookup_table_uniform_type), pointer :: lookup_table
  
  call LookupTableAxisDestroy(lookup_table%axis1)
  call LookupTableAxisDestroy(lookup_table%axis2)
  call LookupTableAxisDestroy(lookup_table%axis3)
  call DeallocateArray(lookup_table%data)
  deallocate(lookup_table)
  nullify(lookup_table)

end subroutine LookupTableUniformDestroy

! ************************************************************************** !

subroutine LookupTableGeneralDestroy(lookup_table)
  ! 
  ! Deallocates any allocated pointers in lookup table
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  use Utility_module
  
  implicit none
  
  class(lookup_table_general_type), pointer :: lookup_table
  
  if (.not.associated(lookup_table)) return

  call LookupTableAxisDestroy(lookup_table%axis1)
  ! axis2 is different type
  if (associated(lookup_table%axis2)) then
    call DeallocateArray(lookup_table%axis2%values)
    deallocate(lookup_table%axis2)
    nullify(lookup_table%axis2)
  endif
  call DeallocateArray(lookup_table%data)
  deallocate(lookup_table)
  nullify(lookup_table)

end subroutine LookupTableGeneralDestroy

end module Lookup_Table_module
