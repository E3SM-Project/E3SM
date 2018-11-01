module Lookup_Table_module
  
#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module

  implicit none
  
  private

  ! Variable extrapolation types
  PetscInt, parameter, public :: VAR_EXTRAP_CONST_VAL = 1
  PetscInt, parameter, public :: VAR_EXTRAP_CONST_GRAD = 2
  ! Variable interpolation types
  PetscInt, parameter :: VAR_INTERP_LINEAR = 1
  PetscInt, parameter :: VAR_INTERP_X_LINLOG = 2
  
  type, abstract, public :: lookup_table_base_type
    PetscInt :: dim
    PetscInt :: dims(3)
    PetscReal, pointer :: data(:)
    PetscReal, pointer :: var_data(:,:)
    type(lookup_table_var_ptr_type), pointer :: var_array(:)
    type(lookup_table_var_list_type), pointer :: vars    
    class(lookup_table_axis_type), pointer :: axis1
  contains
    procedure(LookupTableEvaluateDummy), deferred, public :: Sample
    procedure(LookupTableValAndGradDummy),deferred,public :: SampleAndGradient
    procedure, public :: LookupTableVarConvFactors  
    procedure, public :: LookupTableVarsInit
    procedure, public :: VarPointAndUnitConv
    procedure, public :: SetupConstGradExtrap
    procedure, public :: LookupTableVarInitGradients
    procedure, public :: SetupVarLinLogInterp
  end type lookup_table_base_type
  
  type, public, extends(lookup_table_base_type) :: lookup_table_uniform_type
    class(lookup_table_axis_type), pointer :: axis2
    class(lookup_table_axis_type), pointer :: axis3
  contains
    procedure, public :: Sample => LookupTableEvaluateUniform
    procedure, public :: SampleAndGradient => ValAndGradUniform
  end type lookup_table_uniform_type
  
  type, public, extends(lookup_table_base_type) :: lookup_table_general_type
    class(lookup_table_axis2_general_type), pointer :: axis2
  contains
    procedure, public :: Sample => LookupTableEvaluateGeneral
    procedure, public :: SampleAndGradient => ValAndGradGeneral
  end type lookup_table_general_type
  
  type, public :: lookup_table_axis_type
    PetscInt :: itype
    PetscInt :: saved_index
    PetscReal, pointer :: values(:)
  end type lookup_table_axis_type
  
  type, public, extends(lookup_table_axis_type) :: lookup_table_axis2_general_type
    PetscInt :: saved_index2
  end type lookup_table_axis2_general_type

  type, public :: lookup_table_var_type
    PetscInt :: id
    PetscInt :: iname
    PetscInt :: data_idx
    PetscInt :: extrapolation_itype
    PetscInt :: interp_type
    character(len=MAXWORDLENGTH) :: internal_units
    character(len=MAXWORDLENGTH) :: user_units
    PetscReal :: conversion_factor
    PetscReal, pointer :: data(:)
    PetscReal :: sample
    PetscReal, pointer :: sample_grad(:)
    type(lookup_table_var_type), pointer :: next
  end type  lookup_table_var_type

  type, public :: lookup_table_var_ptr_type
    type(lookup_table_var_type), pointer :: ptr
  end type lookup_table_var_ptr_type

  type, public :: lookup_table_var_list_type
    PetscInt :: num_lookup_table_vars
    type(lookup_table_var_type), pointer :: first
    type(lookup_table_var_type), pointer :: last
    !type(lookup_table_var_type), pointer :: array(:)
  end type lookup_table_var_list_type
  
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
    
    subroutine LookupTableValAndGradDummy(this,var_iname,lookup1,lookup2,lookup3)
      import lookup_table_base_type
      implicit none
      class(lookup_table_base_type) :: this
      PetscInt, intent(in) :: var_iname
      PetscReal :: lookup1
      PetscReal, optional :: lookup2
      PetscReal, optional :: lookup3
    end subroutine LookupTableValAndGradDummy
    
  end interface
  
  interface LookupTableTest
    module procedure LookupTableTest1D
    module procedure LookupTableTest2D
  end interface

  interface LookupTableCreateGeneral
    module procedure LookupTableCreateGeneralDim
    module procedure LookupTableCreateGeneralNoDim  
  end interface

  interface LookupTableDestroy
    module procedure LookupTableUniformDestroy
    module procedure LookupTableGeneralDestroy
  end interface

  public :: LookupTableCreateUniform, &
            LookupTableCreateGeneral, &
            LookupTableAxisInit, &
            LookupTableDestroy, &
            LookupTableTest, &
            LookupTableVarInitList, &
            CreateLookupTableVar, &
            LookupTableVarAddToList, &
            LookupTableVarListDestroy            
  
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
  nullify(lookup_table%var_data)
  nullify(lookup_table%var_array)
  nullify(lookup_table%vars)
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

function LookupTableCreateGeneralDim(dim)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14 - 
  ! Modified by Paolo Orsini: change of function name to add new constructor
  ! Date: 05/02/18
  ! 

  implicit none
  
  PetscInt :: dim
  
  class(lookup_table_general_type), pointer :: LookupTableCreateGeneralDim

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
  
  LookupTableCreateGeneralDim => lookup_table

end function LookupTableCreateGeneralDim

! ************************************************************************** !

function LookupTableCreateGeneralNoDim()
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/02/18
  ! 
  ! contructor before knowing the lookup table dims

  implicit none
  
  class(lookup_table_general_type), pointer :: LookupTableCreateGeneralNoDim

  class(lookup_table_general_type), pointer :: lookup_table
  
  allocate(lookup_table)
  call LookupTableBaseInit(lookup_table)
  nullify(lookup_table%axis2)
  
  LookupTableCreateGeneralNoDim => lookup_table

end function LookupTableCreateGeneralNoDim

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

subroutine ValAndGradUniform(this,var_iname,lookup1,lookup2,lookup3)
  ! 
  ! Computes value and gradient for given coordinates
  !
  ! Author: Paolo Orsini
  ! Date: 05/11/18
  ! 
  implicit none
  
  class(lookup_table_uniform_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscReal :: lookup1
  PetscReal, optional :: lookup2
  PetscReal, optional :: lookup3

  call LookupTableIndexUniform(this,lookup1,lookup2,lookup3)
  if (present(lookup3)) then
    !3D interpolation and gradient computation not yet supported
  else if (present(lookup2)) then
    call InterpExtrapGrad2DUniform(this,var_iname,lookup1,lookup2)
  else
    call InterpExtrapGrad1D(this,var_iname,lookup1)
  endif
  
end subroutine ValAndGradUniform

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

subroutine ValAndGradGeneral(this,var_iname,lookup1,lookup2,lookup3)
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/21/18
  ! 
  implicit none
  
  class(lookup_table_general_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscReal :: lookup1
  PetscReal, optional :: lookup2
  PetscReal, optional :: lookup3  
  
  call LookupTableIndexGeneral(this,lookup1,lookup2,lookup3)
  
  if (present(lookup3)) then
    ! 3D general lookup table not supported
  else if (present(lookup2)) then
    call InterpExtrapGradGeneral2D(this,var_iname,lookup1,lookup2)
  else
    call InterpExtrapGrad1D(this,var_iname,lookup1)
  end if    
  
  
end subroutine ValAndGradGeneral

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

subroutine InterpExtrapGrad1D(this,var_iname,lookup1)
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/12/18
  ! 
  use Utility_module, only : Interpolate, GradientLinear
  
  implicit none
  
  class(lookup_table_base_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscReal, intent(in) :: lookup1
  
  PetscReal :: x1, x2, x1_grad, x2_grad
  PetscInt :: i1, i2, i1_grad, i2_grad
  PetscInt :: var_idx, sizei
  
  i1 = this%axis1%saved_index
  sizei = this%dims(1)
  
  if (i1 > 0) then
    i2 = max(min(i1+1,this%dims(1)),1)
  else
    i1=1
    i2=1
  end if    
  
  if (i1 == i2) then !end points
    if (i1 == 1) then
      i1_grad = 1
      i2_grad = 2
    else if (i1 == sizei) then
      i1_grad = sizei - 1
      i2_grad = sizei
    end if
    !x1 = this%axis1%values(i1)
    x1_grad = this%axis1%values(i1_grad)
    x2_grad = this%axis1%values(i2_grad)
    !x_fract_grad = (lookup1 - x1_grad) / ( x2_grad - x1_grad )
    if (var_iname == -1 ) then
      do var_idx = 1,size(this%var_array(:))
        if ( associated(this%var_array(var_idx)%ptr) ) then
          call Var1DExtrapolate(this,var_idx,i1,i1_grad,i2_grad, &
                                      x1_grad,x2_grad,lookup1)
        end if     
      end do
    else 
      call Var1DExtrapolate(this,var_iname,i1,i1_grad,i2_grad, &
                                  x1_grad,x2_grad,lookup1)
    end if
  else ! away from end points
    x1 = this%axis1%values(i1)
    x2 = this%axis1%values(i2)
    if (var_iname == -1 ) then
      do var_idx = 1,size(this%var_array(:))
        if ( associated(this%var_array(var_idx)%ptr) ) then
          call Var1DInterpolate(this,var_idx,i1,i2,x1,x2,lookup1)
        end if     
      end do
    else
      call Var1DInterpolate(this,var_iname,i1,i2,x1,x2,lookup1)
    end if    
  end if  
  
end subroutine InterpExtrapGrad1D

! ************************************************************************** !

subroutine Var1DExtrapolate(this,var_iname,i1,i1_grad,i2_grad, &
                            x1_grad,x2_grad,lookup)
  ! 
  ! Author: Paolo Orsini
  ! Date: 02/06/18
  !   
  implicit none

  class(lookup_table_base_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscInt, intent(in) :: i1,i1_grad,i2_grad
  PetscReal, intent(in) :: x1_grad,x2_grad
  PetscReal, intent(in) :: lookup
          
  if ( this%var_array(var_iname)%ptr%extrapolation_itype == &          
       VAR_EXTRAP_CONST_GRAD ) then
    call Var1DInterpolate(this,var_iname,i1_grad,i2_grad, &
                                   x1_grad,x2_grad,lookup)
  else if( this%var_array(var_iname)%ptr%extrapolation_itype == &
         VAR_EXTRAP_CONST_VAL) then
    this%var_array(var_iname)%ptr%sample = &
                         this%var_array(var_iname)%ptr%data(i1)
    this%var_array(var_iname)%ptr%sample_grad(1) = 0.0d0
      
  end if

end subroutine Var1DExtrapolate

! ************************************************************************** !

subroutine Var1DInterpolate(this,var_iname,i1,i2,x1,x2,lookup)
  ! 
  ! Author: Paolo Orsini
  ! Date: 02/06/18
  ! 
  use Utility_module, only : Interpolate, GradientLinear
  
  implicit none
    
  class(lookup_table_base_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscInt, intent(in) :: i1,i2
  PetscReal, intent(in) :: x1,x2
  PetscReal, intent(in) :: lookup

  PetscReal :: val1, val2

  val1 = this%var_array(var_iname)%ptr%data(i1)
  val2 = this%var_array(var_iname)%ptr%data(i2)
  if (this%var_array(var_iname)%ptr%interp_type  == VAR_INTERP_X_LINLOG) then
     val1 = dlog(val1)
     val2 = dlog(val2)
  end if   
  call GradientLinear(x2,x1,val2,val1, &
                      this%var_array(var_iname)%ptr%sample_grad(1))
  call Interpolate(x2,x1,lookup,val2,val1, &
                   this%var_array(var_iname)%ptr%sample)
  if (this%var_array(var_iname)%ptr%interp_type  == VAR_INTERP_X_LINLOG) then
    this%var_array(var_iname)%ptr%sample = &
        dexp(this%var_array(var_iname)%ptr%sample)
    this%var_array(var_iname)%ptr%sample_grad(1) = &
        this%var_array(var_iname)%ptr%sample_grad(1) * &
        this%var_array(var_iname)%ptr%sample
  end if                 
  
end subroutine Var1DInterpolate

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

subroutine InterpExtrapGrad2DUniform(this,var_iname,lookup1,lookup2)
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/25/18
  ! 

  implicit none
  
  class(lookup_table_uniform_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscReal :: lookup1
  PetscReal :: lookup2
  class(lookup_table_axis_type), pointer :: axis1 => null()
  class(lookup_table_axis_type), pointer :: axis2 => null()

  PetscInt :: i1, i2, j1, j2
  PetscInt :: i1_grad, i2_grad, j1_grad, j2_grad
  PetscReal :: x1, x2
  PetscReal :: x1_grad, x2_grad
  PetscReal :: x_frac
  PetscInt :: sizei, sizej
  PetscInt :: var_idx

  !  x1,y2,z3 ------ x2,y2,z4
  !     |               |
  !     |               |
  !     |   x,y         |
  !     |               |
  !  x1,y1,z1 ------ x2,y1,z2


  axis1 => this%axis1
  axis2 => this%axis2
  
  sizei = this%dims(1)
  sizej = this%dims(2)
  i1 = axis1%saved_index
  j1 = axis2%saved_index
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
  
  if (j2 == j1) then
    if ( j1 == 1 ) then
      j1_grad = 1
      j2_grad = 2
    else if ( j1 == sizej ) then
      j1_grad = sizej -1
      j2_grad = sizej
    end if
  end if

  !compute i indices for interpolation/extrapolation
  if (i2 == i1) then
    !perform extrapolation in the x-direction
    if ( i1 == 1) then
      i1_grad = 1
      i2_grad = 1 + 1
    else if ( i1 == sizei ) then
      i1_grad = sizei - 1
      i2_grad = sizei
    end if
    !extrapolation in the x-direction
    x1_grad = axis1%values(i1_grad)
    x2_grad = axis1%values(i2_grad)
    x_frac = (lookup1 - x1_grad) / ( x2_grad - x1_grad )    
    if (var_iname == -1 ) then
      do var_idx = 1,size(this%var_array(:))
        if ( associated(this%var_array(var_idx)%ptr) ) then
          call VarUniformXExtrapolate(this,var_idx,i1,i1_grad,i2_grad,sizei, &
                                      j1,j1_grad,j2,j2_grad,&
                                      x1_grad,x2_grad,x_frac,lookup1,lookup2)        
        end if  
      end do
    else
      call VarUniformXExtrapolate(this,var_iname,i1,i1_grad,i2_grad,sizei, &
                                  j1,j1_grad,j2,j2_grad,&
                                  x1_grad,x2_grad,x_frac,lookup1,lookup2)      
    end if    
  else !away from end point in the x-direction
    !inteprolation in the x-direction
    x1 = axis1%values(i1)
    x2 = axis1%values(i2)
    x_frac = (lookup1 - x1) / ( x2 - x1 )    
    if (var_iname == -1 ) then
      do var_idx = 1,size(this%var_array(:))
        if ( associated(this%var_array(var_idx)%ptr) ) then    
          call VarUniformXInterpolate(this,var_idx,j1,j1_grad,j2,j2_grad, &
                                     i1,i2,sizei,x1,x2,x_frac,lookup1,lookup2)
        end if
      end do
    else
      call VarUniformXInterpolate(this,var_iname,j1,j1_grad,j2,j2_grad, &
                                 i1,i2,sizei,x1,x2,x_frac,lookup1,lookup2)      
    end if                                   
  end if
  
  nullify(axis1)
  nullify(axis2)
                   
end subroutine InterpExtrapGrad2DUniform

! ************************************************************************** !

subroutine VarUniformXExtrapolate(this,var_iname,i_col,i1_grad,i2_grad,sizei, &
                                  j1,j1_grad,j2,j2_grad,&
                                  x1_grad,x2_grad,x_frac_grad,lookup1,lookup2)
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/25/18
  !     

  implicit none
  
  class(lookup_table_uniform_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscInt, intent(in) :: i_col,i1_grad,i2_grad,sizei
  PetscInt, intent(in) :: j1,j1_grad,j2,j2_grad
  PetscReal, intent(in) :: x1_grad,x2_grad,x_frac_grad,lookup1,lookup2
  
                           
  if (this%var_array(var_iname)%ptr%extrapolation_itype == &
      VAR_EXTRAP_CONST_GRAD ) then
    ! extrapolate  
    call VarUniformXInterpolate(this,var_iname,j1,j1_grad,j2,j2_grad, &
                                i1_grad,i2_grad,sizei,x1_grad,x2_grad, &
                                x_frac_grad,lookup1,lookup2)
  else if ( this%var_array(var_iname)%ptr%extrapolation_itype == & 
            VAR_EXTRAP_CONST_VAL ) then
    !value and gradient in the edge at x = const
    call UniformYValAndGrad(this,var_iname,j1,j2,j1_grad,j2_grad,i_col,sizei, &
                            lookup2,this%var_array(var_iname)%ptr%sample, &              
                            this%var_array(var_iname)%ptr%sample_grad(2))
    this%var_array(var_iname)%ptr%sample_grad(1) = 0.0d0        
    if (this%var_array(var_iname)%ptr%interp_type == &
                                 VAR_INTERP_X_LINLOG ) then
      this%var_array(var_iname)%ptr%sample = &
           dexp(this%var_array(var_iname)%ptr%sample)
      this%var_array(var_iname)%ptr%sample_grad(2) = &
            this%var_array(var_iname)%ptr%sample_grad(2) * &
            this%var_array(var_iname)%ptr%sample
    end if        
  end if

end subroutine VarUniformXExtrapolate

! ************************************************************************** !

subroutine VarUniformXInterpolate(this,var_iname,j1,j1_grad,j2,j2_grad, &
                                 i1,i2,sizei,x1,x2,x_frac,lookup1,lookup2)
! 
! Author: Paolo Orsini
! Date: 05/17/18
!     
  use Utility_module, only : Interpolate, GradientLinear

  implicit none

  class(lookup_table_uniform_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscInt, intent(in) :: j1,j1_grad
  PetscInt, intent(in) :: j2,j2_grad
  PetscInt, intent(in) :: i1,i2,sizei
  PetscReal, intent(in) :: x1,x2,x_frac,lookup1,lookup2

  PetscReal :: val_i1,grad_i1,val_i2,grad_i2

  !left
  call UniformYValAndGrad(this,var_iname,j1,j2,j1_grad,j2_grad,i1, &
                                        sizei,lookup2,val_i1,grad_i1)
  !right
  call UniformYValAndGrad(this,var_iname,j1,j2,j1_grad,j2_grad,i2, &
                                          sizei,lookup2,val_i2,grad_i2)
  call GradientLinear(x2,x1,val_i2,val_i1, &
                      this%var_array(var_iname)%ptr%sample_grad(1))
  call Interpolate(x2,x1,lookup1,val_i2,val_i1, &
                   this%var_array(var_iname)%ptr%sample)
  this%var_array(var_iname)%ptr%sample_grad(2) = &
                       grad_i1 * (1.0 - x_frac) + grad_i2 * x_frac
  if (this%var_array(var_iname)%ptr%interp_type == &
                                      VAR_INTERP_X_LINLOG ) then
    this%var_array(var_iname)%ptr%sample = &
       dexp(this%var_array(var_iname)%ptr%sample)
    this%var_array(var_iname)%ptr%sample_grad(1:2) =  &
       this%var_array(var_iname)%ptr%sample_grad(1:2) * &
       this%var_array(var_iname)%ptr%sample
  end if

end subroutine VarUniformXInterpolate

! ************************************************************************** !

subroutine UniformYValAndGrad(this,var_iname,j1,j2,j1_grad,j2_grad,i_col, &
                              sizei,lookup2,val,grad_val)                              
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/25/18
  !   
  ! Given a data set of point in y this routine computes either:
  ! 1) z and dz_dy
  !    OR
  ! 2) ln(z) and d(ln(z))/dy
  !
  ! depending on the interpolation method chosen for the variable
  !
  use Utility_module, only : Interpolate, GradientLinear

  implicit none
  
  class(lookup_table_uniform_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscInt, intent(in) :: j1,j2
  PetscInt, intent(in) :: j1_grad,j2_grad
  PetscInt, intent(in) :: i_col, sizei
  PetscReal,intent(in) :: lookup2
  PetscReal, intent(out) :: val, grad_val
  
  PetscReal :: y1, y2, z1, z2
  PetscReal :: y1_grad, y2_grad
  
  
  if (j2 == j1) then
    if ( this%var_array(var_iname)%ptr%extrapolation_itype == &
         VAR_EXTRAP_CONST_GRAD ) then
      y1_grad = this%axis2%values(j1_grad)
      y2_grad = this%axis2%values(j2_grad)
      z1 = this%var_array(var_iname)%ptr%data(i_col+(j1_grad-1)*sizei)
      z2 = this%var_array(var_iname)%ptr%data(i_col+(j2_grad-1)*sizei)
      call GradientLinear(y2_grad,y1_grad,z2,z1,grad_val)
      call Interpolate(y2_grad,y1_grad,lookup2,z2,z1,val)
       if ( this%var_array(var_iname)%ptr%interp_type == &
            VAR_INTERP_X_LINLOG ) then
          grad_val = ( 1.0 / val ) * grad_val
          val = dlog(val)
       end if    
    else if ( this%var_array(var_iname)%ptr%extrapolation_itype == &
              VAR_EXTRAP_CONST_VAL ) then
      val = this%var_array(var_iname)%ptr%data(i_col+(j1-1)*sizei)
      if ( this%var_array(var_iname)%ptr%interp_type == &
           VAR_INTERP_X_LINLOG ) then
        val = dlog(val)
      end if              
      grad_val = 0.0
    end if  
  else !away from end points
    y1 = this%axis2%values(j1)
    y2 = this%axis2%values(j2)
    z1 = this%var_array(var_iname)%ptr%data(i_col+(j1-1)*sizei)
    z2 = this%var_array(var_iname)%ptr%data(i_col+(j2-1)*sizei)
    call Interpolate(y2,y1,lookup2,z2,z1,val)
    call GradientLinear(y2,y1,z2,z1,grad_val)
     if ( this%var_array(var_iname)%ptr%interp_type == &
          VAR_INTERP_X_LINLOG ) then
       grad_val = 1.0 / val * grad_val
       val = dlog(val)
     end if    
  endif
  
end subroutine UniformYValAndGrad

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

subroutine InterpExtrapGradGeneral2D(this,var_iname,lookup1,lookup2)
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/17/18
  ! 
  use Utility_module, only : Interpolate, GradientLinear

  implicit none
  
  class(lookup_table_general_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscReal, intent(in) :: lookup1
  PetscReal, intent(in) :: lookup2
  
  PetscInt :: var_idx
  PetscInt :: ja
  PetscInt :: i1a, i1b
  PetscInt :: jb
  PetscInt :: i2a, i2b
  PetscReal :: xa, xb
  PetscInt :: sizei, sizej
  PetscInt :: ja_grad, jb_grad
  PetscInt :: i1a_grad, i2a_grad, i1b_grad, i2b_grad
  PetscReal :: xa_grad, xb_grad, x_frac

  
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
  !
  !    ja     jb
  ! note that j refer to x
 
  ! x used for both i and j to indicate all operations are 1D
  sizei = this%dims(2)  ! i referes to y 
  sizej = this%dims(1)  ! j refers to x
  x_frac = 0.0
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
  
  if ( i1a == i2a ) then
    if ( i1a == 1 ) then
      i1a_grad = 1
      i2a_grad = 2
    else if (i1a == sizei) then
      i1a_grad = sizei - 1
      i2a_grad = sizei
    end if  
  end if  

  if ( i1b == i2b ) then
    if ( i1b == 1 ) then
      i1b_grad = 1
      i2b_grad = 2
    else if (i1b == sizei) then
      i1b_grad = sizei - 1
      i2b_grad = sizei
    end if  
  end if
  
  if (jb == ja) then
  !left and right table edges (included the corners that are extrapolated)
    if ( ja == 1 ) then
      ja_grad = 1
      jb_grad = 2
      !new look up to find i1b_grad, i2b_grad 
      call GeneralFindIndicesAxis2(this,jb_grad,lookup2,i1b,i2b, &
                                   i1b_grad,i2b_grad)
    else if ( ja == sizej ) then
      ja_grad = sizej - 1
      jb_grad = sizej
      !new look up to find i1a_grad, i2a_grad 
      call GeneralFindIndicesAxis2(this,ja_grad,lookup2,i1b,i2b, &
                                   i1a_grad,i2a_grad)
    end if
    xa_grad = this%axis1%values(ja_grad)
    xb_grad = this%axis1%values(jb_grad)
    x_frac = (lookup1 - xa_grad) / (xb_grad - xa_grad)
    !interpolation/extrapolation left and right columns             
    if (var_iname == -1 ) then
      do var_idx = 1,size(this%var_array(:))
        if ( associated(this%var_array(var_idx)%ptr) ) then
          call VarGeneralXExtrapolate(this,var_idx,ja,i1a,i2a,sizei, &
                                      ja_grad,i1a_grad,i2a_grad, &
                                      jb_grad,i1b_grad,i2b_grad, &
                                      xa_grad,xb_grad,x_frac,lookup1,lookup2)
        end if
      end do                             
    else
      ! to call function var_extrapolate function
      call VarGeneralXExtrapolate(this,var_iname,ja,i1a,i2a,sizei, &
                                        ja_grad,i1a_grad,i2a_grad, &
                                        jb_grad,i1b_grad,i2b_grad, &
                                        xa_grad,xb_grad,x_frac,lookup1,lookup2) 
    end if
  else !away from left/right table edges - no extrapolation in x directions
    xa = this%axis1%values(ja)
    xb = this%axis1%values(jb)
    x_frac = ( lookup1 - xa ) / ( xb - xa )
    if (var_iname == -1 ) then
      do var_idx = 1,size(this%var_array(:))
        if ( associated(this%var_array(var_idx)%ptr) ) then
          call VarGeneralXInterpolate(this,var_idx,ja,i1a,i2a, &
                     i1a_grad,i2a_grad,sizei,jb,i1b,i2b,i1b_grad,i2b_grad, &
                     xa,xb,x_frac,lookup1,lookup2) 
        end if
      end do          
    else
      call VarGeneralXInterpolate(this,var_iname,ja,i1a,i2a, &
                 i1a_grad,i2a_grad,sizei,jb,i1b,i2b,i1b_grad,i2b_grad, &
                 xa,xb,x_frac,lookup1,lookup2)
    end if
  end if        
                   
end subroutine InterpExtrapGradGeneral2D

! ************************************************************************** !

subroutine VarGeneralXExtrapolate(this,var_iname,ja,i1a,i2a,sizei, &
                                  ja_grad,i1a_grad,i2a_grad, &
                                  jb_grad,i1b_grad,i2b_grad, &
                                  xa_grad,xb_grad,x_frac,lookup1,lookup2)
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/17/18
  !     
  use Utility_module, only : Interpolate, GradientLinear

  implicit none
  
  class(lookup_table_general_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscInt, intent(in) :: ja,i1a,i2a,sizei
  PetscInt, intent(in) :: ja_grad,i1a_grad,i2a_grad
  PetscInt, intent(in) :: jb_grad,i1b_grad,i2b_grad
  PetscReal, intent(in) :: xa_grad,xb_grad,x_frac,lookup1,lookup2
  
  !PetsReal :: val_a,val_b,grad_a,grad_b
                           
  if (this%var_array(var_iname)%ptr%extrapolation_itype == &
      VAR_EXTRAP_CONST_GRAD ) then

    call VarGeneralXInterpolate(this,var_iname,ja_grad,i1a,i2a,&
                  i1a_grad,i2a_grad,sizei,jb_grad,i1a,i2a,i1b_grad,i2b_grad, &
                              xa_grad,xb_grad,x_frac,lookup1,lookup2)
  else if ( this%var_array(var_iname)%ptr%extrapolation_itype == & 
            VAR_EXTRAP_CONST_VAL ) then
    !value and gradient in the edge at x = const
    call GeneralYValAndGrad(this,var_iname,ja,i1a,i2a,i1a_grad,i2a_grad, &
                             sizei,lookup2, &
                             this%var_array(var_iname)%ptr%sample, &
                             this%var_array(var_iname)%ptr%sample_grad(2))
    this%var_array(var_iname)%ptr%sample_grad(1) = 0.0d0                         
    if (this%var_array(var_iname)%ptr%interp_type == &
                                 VAR_INTERP_X_LINLOG ) then
      this%var_array(var_iname)%ptr%sample = &
           dexp(this%var_array(var_iname)%ptr%sample)
      this%var_array(var_iname)%ptr%sample_grad(2) = &
            this%var_array(var_iname)%ptr%sample_grad(2) * &
            this%var_array(var_iname)%ptr%sample
    end if        
  end if

end subroutine VarGeneralXExtrapolate

! ************************************************************************** !

subroutine VarGeneralXInterpolate(this,var_iname,ja,i1a,i2a,i1a_grad,i2a_grad, &
                                 sizei,jb,i1b,i2b,i1b_grad,i2b_grad,&
                                 xa,xb,x_frac,lookup1,lookup2)
! 
! Author: Paolo Orsini
! Date: 05/17/18
!     
  use Utility_module, only : Interpolate, GradientLinear

  implicit none

  class(lookup_table_general_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscInt, intent(in) :: ja,i1a,i2a,i1a_grad,i2a_grad,sizei
  PetscInt, intent(in) :: jb,i1b,i2b,i1b_grad,i2b_grad
  PetscReal, intent(in) :: xa,xb,x_frac,lookup1,lookup2

  PetscReal :: val_a,grad_a,val_b,grad_b

  !left
  call GeneralYValAndGrad(this,var_iname,ja,i1a,i2a,i1a_grad,i2a_grad, &
                          sizei,lookup2,val_a,grad_a)
  !right  
  call GeneralYValAndGrad(this,var_iname,jb,i1b,i2b,i1b_grad,i2b_grad, &
                          sizei,lookup2,val_b,grad_b)
  call GradientLinear(xb,xa,val_b,val_a, &
                      this%var_array(var_iname)%ptr%sample_grad(1))
  call Interpolate(xb,xa,lookup1,val_b,val_a, &
                   this%var_array(var_iname)%ptr%sample)
  this%var_array(var_iname)%ptr%sample_grad(2) = &
                       grad_a * (1.0 - x_frac) + grad_b * x_frac
  if (this%var_array(var_iname)%ptr%interp_type == &
                                      VAR_INTERP_X_LINLOG ) then
    this%var_array(var_iname)%ptr%sample = &
       dexp(this%var_array(var_iname)%ptr%sample)
    this%var_array(var_iname)%ptr%sample_grad(1:2) =  &
       this%var_array(var_iname)%ptr%sample_grad(1:2) * &
       this%var_array(var_iname)%ptr%sample
  end if

end subroutine VarGeneralXInterpolate

! ************************************************************************** !

subroutine GeneralYValAndGrad(this,var_iname,j_col,i1,i2,i1_grad,i2_grad, &
                              sizei,lookup,val,grad_val)
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/17/18
  !   
  ! Given a data set of point in y this routine computes either:
  ! 1) z and dz_dy
  !    OR
  ! 2) ln(z) and d(ln(z))/dy
  !
  ! depending on the interpolation method chosen for the variable
  !
  use Utility_module, only : Interpolate, GradientLinear

  implicit none
  
  class(lookup_table_general_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscInt, intent(in) :: j_col
  PetscInt, intent(in) :: i1,i2
  PetscInt, intent(in) :: i1_grad, i2_grad
  PetscInt, intent(in) :: sizei
  PetscReal,intent(in) :: lookup
  PetscReal, intent(out) :: val, grad_val
  
  !PetscInt :: i1_grad, i2_grad
  PetscInt :: iia, iib
  PetscReal :: x1, x2, z1, z2
  PetscReal :: x1_grad, x2_grad

  
  if (i2 == i1) then
    if ( this%var_array(var_iname)%ptr%extrapolation_itype == &
         VAR_EXTRAP_CONST_GRAD ) then
      iia = i1+(j_col-1)*sizei
      x1 = this%axis2%values(iia)
      iia = i1_grad+(j_col-1)*sizei
      iib = i2_grad+(j_col-1)*sizei
      x1_grad = this%axis2%values(iia)
      x2_grad = this%axis2%values(iib)
      z1 = this%var_array(var_iname)%ptr%data(iia)
      z2 = this%var_array(var_iname)%ptr%data(iib)
      !log numerical derivatives for testing
      ! if ( this%var_array(var_iname)%ptr%interp_type == &
      !      VAR_INTERP_X_LINLOG ) then
      !   z1 = dlog(z1)
      !   z2 = dlog(z2)
      ! end if     
      call GradientLinear(x2_grad,x1_grad,z2,z1,grad_val)
      call Interpolate(x2_grad,x1_grad,lookup,z2,z1,val)
      !log analytical derivatives
      if ( this%var_array(var_iname)%ptr%interp_type == &
           VAR_INTERP_X_LINLOG ) then
        grad_val = ( 1.0 / val) * grad_val
        val = dlog(val)
      end if           
    else if ( this%var_array(var_iname)%ptr%extrapolation_itype == &
              VAR_EXTRAP_CONST_VAL ) then
      val = this%var_array(var_iname)%ptr%data(i1+(j_col-1)*sizei)
      if ( this%var_array(var_iname)%ptr%interp_type == &
           VAR_INTERP_X_LINLOG ) then
        val = dlog(val)   
      end if           
      grad_val = 0.0
    end if  
  else !away from end points
    iia = i1+(j_col-1)*sizei
    iib = i2+(j_col-1)*sizei
    x1 = this%axis2%values(iia)
    x2 = this%axis2%values(iib)
    z1 = this%var_array(var_iname)%ptr%data(iia)
    z2 = this%var_array(var_iname)%ptr%data(iib)
    !log numerical derivatives for testing
    ! if ( this%var_array(var_iname)%ptr%interp_type == &
    !      VAR_INTERP_X_LINLOG ) then
    !   z1 = dlog(z1)
    !   z2 = dlog(z2)
    ! end if         
    call Interpolate(x2,x1,lookup,z2,z1,val)
    call GradientLinear(x2,x1,z2,z1,grad_val)
    if ( this%var_array(var_iname)%ptr%interp_type == &
          VAR_INTERP_X_LINLOG ) then
      grad_val = 1.0 / val * grad_val 
      val = dlog(val)
    end if      
  endif
  
end subroutine GeneralYValAndGrad

! ************************************************************************** !

subroutine GeneralFindIndicesAxis2(this,j_col,lookup,i1,i2,i1_grad,i2_grad)
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/21/18
  !   
  implicit none  
  
  class(lookup_table_general_type) :: this
  PetscInt, intent(in)  :: j_col
  PetscReal, intent(in)  :: lookup
  PetscInt, intent(out) :: i1
  PetscInt, intent(out) :: i2
  PetscInt, intent(out) :: i1_grad
  PetscInt, intent(out) :: i2_grad

  PetscInt :: iend, istart, sizei
  
  sizei = this%dims(2)
  iend = j_col*this%dims(2)
  istart = iend - this%dims(2) + 1
  call LookupTableAxisIndexGeneral(lookup,this%axis2%values(istart:iend),i1)
  if (i1 > 0) then
    i2 = max(min(i1+1,sizei),1)
  else
    i1 = 1
    i2 = 1
  end if
  
  if ( i1 == i2 ) then
    if ( i1 == 1 ) then
      i1_grad = 1
      i2_grad = 2
    else if ( i1 == sizei) then
      i1_grad = sizei - 1
      i2_grad = sizei
    end if    
  end if  
  
end subroutine GeneralFindIndicesAxis2

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

subroutine LookupTableBaseDestroy(lookup_table)
  !
  ! Deallocates any allocated pointers in lookup table base type
  !
  ! Author: Paolo Orsini
  ! Date: 12/01/17
  !
  use Utility_module

  implicit none

  PetscInt :: i

  class(lookup_table_base_type) :: lookup_table

  call DeallocateArray(lookup_table%data)
  call DeallocateArray(lookup_table%var_data)
  call LookupTableAxisDestroy(lookup_table%axis1)

  if (associated(lookup_table%var_array)) then
    do i =1,size(lookup_table%var_array(:))
      nullify(lookup_table%var_array(i)%ptr)
    end do
    deallocate(lookup_table%var_array)
    nullify(lookup_table%var_array)
  end if

  call LookupTableVarListDestroy(lookup_table%vars)

end subroutine LookupTableBaseDestroy

! ************************************************************************** !

subroutine LookupTableUniformDestroy(lookup_table)
  ! 
  ! Deallocates any allocated pointers in lookup table
  ! 
  ! Author: Glenn Hammond, Paolo Orsini
  ! Date: 10/15/14, 04/18/2018
  ! 
  use Utility_module
  
  implicit none
  
  class(lookup_table_uniform_type), pointer :: lookup_table
  
  call LookupTableBaseDestroy(lookup_table)
  call LookupTableAxisDestroy(lookup_table%axis2)
  call LookupTableAxisDestroy(lookup_table%axis3)
  deallocate(lookup_table)
  nullify(lookup_table)

end subroutine LookupTableUniformDestroy

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

  call LookupTableBaseDestroy(lookup_table)
  ! axis2 is different type
  if (associated(lookup_table%axis2)) then
    call DeallocateArray(lookup_table%axis2%values)
    deallocate(lookup_table%axis2)
    nullify(lookup_table%axis2)
  endif
  deallocate(lookup_table)
  nullify(lookup_table)

end subroutine LookupTableGeneralDestroy

! ************************************************************************** !
! ** lookup table variable procedures
! ************************************************************************** !

subroutine LookupTableVarsInit(this,n_var_max)
  !
  ! Initializes a lookup table variables
  !
  ! Author: Paolo Orsini
  ! Date: 05/04/18
  !

  implicit none

  class(lookup_table_base_type) :: this
  PetscInt :: n_var_max
  
  PetscInt :: i_var
  
  allocate(this%vars)
  call LookupTableVarInitList(this%vars)
  allocate(this%var_array(n_var_max))
  do i_var = 1, n_var_max
    nullify(this%var_array(i_var)%ptr)
  end do

end subroutine LookupTableVarsInit

! ************************************************************************** !

subroutine LookupTableVarInitList(list)
  !
  ! Initializes a lookup table var list
  !
  ! Author: Paolo Orsini
  ! Date: 12/01/17
  !

  implicit none

  type(lookup_table_var_list_type) :: list

  nullify(list%first)
  nullify(list%last)
  list%num_lookup_table_vars = 0

end subroutine LookupTableVarInitList

! ************************************************************************** !

subroutine LookupTableVarAddToList(new_var,list)
  !
  ! Adds a new lookup table var to a lookup table vars list
  !
  ! Author: Paolo Orsini
  ! Date: 12/01/2017
  !

  implicit none

  type(lookup_table_var_type), pointer :: new_var
  type(lookup_table_var_list_type) :: list

  list%num_lookup_table_vars = list%num_lookup_table_vars + 1
  new_var%id = list%num_lookup_table_vars
  if (.not.associated(list%first)) list%first => new_var
  if (associated(list%last)) list%last%next => new_var
  list%last => new_var

end subroutine LookupTableVarAddToList

! ************************************************************************** !

function CreateLookupTableVar(var_iname,internal_units,user_units,data_idx)
  !
  ! Creates a new lookup table var
  !
  ! Author: Paolo Orsini
  ! Date: 12/01/2017 - modfied 04/17/2018
  !

  implicit none

  type(lookup_table_var_type), pointer :: CreateLookupTableVar
  PetscInt, intent(in) :: var_iname
  character(len=MAXWORDLENGTH), intent(in) :: internal_units
  character(len=MAXWORDLENGTH), intent(in) :: user_units
  PetscInt, intent(in) :: data_idx

  allocate(CreateLookupTableVar)

  CreateLookupTableVar%id = UNINITIALIZED_INTEGER
  CreateLookupTableVar%iname = var_iname
  CreateLookupTableVar%data_idx = data_idx
  CreateLookupTableVar%extrapolation_itype = UNINITIALIZED_INTEGER
  CreateLookupTableVar%interp_type = VAR_INTERP_LINEAR
  CreateLookupTableVar%internal_units = trim(internal_units)
  CreateLookupTableVar%user_units = trim(user_units)
  CreateLookupTableVar%conversion_factor = UNINITIALIZED_DOUBLE
  nullify(CreateLookupTableVar%data)
  CreateLookupTableVar%sample = UNINITIALIZED_DOUBLE
  nullify(CreateLookupTableVar%sample_grad)
  nullify(CreateLookupTableVar%next)

end function CreateLookupTableVar

! ************************************************************************** !

subroutine LookupTableVarInitGradients(this,option)
  !
  ! allocate lookup variable gradients
  !
  ! Author: Paolo Orsini
  ! Date: 06/04/2018
  !
  
  use Option_module

  implicit none
  
  class(lookup_table_base_type) :: this
  type(option_type) :: option

  PetscInt :: prop_idx

  if ( this%dim > 0 ) then
    do prop_idx = 1,size(this%var_array(:))
      if ( associated(this%var_array(prop_idx)%ptr) ) then
        allocate(this%var_array(prop_idx)%ptr%sample_grad(this%dim))
        this%var_array(prop_idx)%ptr%sample_grad = UNINITIALIZED_DOUBLE
      end if  
    end do
  else
    option%io_buffer = "LookupTableVarInitGradients: cannot initialise " // &
                        "var gradient before defining table dims"
    call printErrMsg(option)    
  end if  

end subroutine LookupTableVarInitGradients

! ************************************************************************** !

subroutine LookupTableVarConvFactors(this,option)
  !
  ! Compute conversion varctors for all lookup table variables
  !
  ! Author: Paolo Orsini
  ! Date: 12/01/2017
  !

  use Option_module
  use Units_module

  implicit none

  class(lookup_table_base_type) :: this
  type(option_type) :: option

  type(lookup_table_var_type), pointer :: var

  var => this%vars%first

  do
    if (.not.associated(var)) exit
      var%conversion_factor = &
          UnitsConvertToInternal(var%user_units,var%internal_units,option)
      var => var%next
  enddo

end subroutine LookupTableVarConvFactors

! ************************************************************************** !

subroutine VarPointAndUnitConv(this,option)
  !
  ! Points variable data arrays to their table columns given the locations
  ! Convert variable data units given the convrsion factors 
  ! previously computed
  !
  ! Author: Paolo Orsini
  ! Date: 05/04/2018
  !

  use Option_module
  use Units_module

  implicit none

  class(lookup_table_base_type) :: this
  type(option_type) :: option

  PetscInt :: prop_idx, data_idx

  do prop_idx = 1,size(this%var_array(:))
    if ( associated(this%var_array(prop_idx)%ptr) ) then
      data_idx = this%var_array(prop_idx)%ptr%data_idx
      this%var_array(prop_idx)%ptr%data => this%var_data(data_idx,:)
      this%var_data(data_idx,:) = this%var_data(data_idx,:) * &
                              this%var_array(prop_idx)%ptr%conversion_factor
    end if
  end do

end subroutine VarPointAndUnitConv

! ************************************************************************** !

subroutine SetupConstGradExtrap(this,option)
  !
  ! Points variable data arrays to their table columns given the locations
  ! Convert variable data units given the convrsion factors 
  ! previously computed
  !
  ! Author: Paolo Orsini
  ! Date: 05/04/2018
  !

  use Option_module

  implicit none

  class(lookup_table_base_type) :: this
  type(option_type) :: option

  PetscInt :: prop_idx
  
  do prop_idx = 1,size(this%var_array(:))
    if ( associated(this%var_array(prop_idx)%ptr) ) then
      this%var_array(prop_idx)%ptr%extrapolation_itype = VAR_EXTRAP_CONST_GRAD
    end if
  end do  
  
end subroutine SetupConstGradExtrap
    
! ************************************************************************** !    

subroutine SetupVarLinLogInterp(this,var_iname,option)
  ! 
  ! Author: Paolo Orsini
  ! Date: 06/06/18
  ! 
  
  use Option_module
  
  implicit none
  
  class(lookup_table_base_type) :: this
  PetscInt, intent(in) :: var_iname
  type(option_type) :: option
  
  if ( associated(this%var_array) ) then
    if ( associated(this%var_array(var_iname)%ptr) ) then
      this%var_array(var_iname)%ptr%interp_type = VAR_INTERP_X_LINLOG
    else
      option%io_buffer = "SetupVarLinLogInterp: cannot setup " // &
            "LinLog inteprolation method for a var not defined as lookupvar"
      call printErrMsg(option)      
    end if  
  end if

end subroutine SetupVarLinLogInterp

! ************************************************************************** !

subroutine LookupTableVarListDestroy(lookup_table_vars_list)
  !
  ! Deallocates a list of lookup table vars
  !
  ! Author: Paolo Orsini
  ! Date: 12/01/2017
  !

  implicit none

  type(lookup_table_var_list_type), pointer :: lookup_table_vars_list

  type(lookup_table_var_type), pointer :: var, prev_var

  if (.not.associated(lookup_table_vars_list)) return

  var => lookup_table_vars_list%first
  do
    if (.not.associated(var)) exit
    prev_var => var
    var => var%next
    call LookupTableVarDestroy(prev_var)
  enddo

  lookup_table_vars_list%num_lookup_table_vars = 0
  nullify(lookup_table_vars_list%first)
  nullify(lookup_table_vars_list%last)

  deallocate(lookup_table_vars_list)
  nullify(lookup_table_vars_list)

end subroutine LookupTableVarListDestroy

! ************************************************************************** !

subroutine LookupTableVarDestroy(lookup_table_var)
  !
  ! Deallocates all members of a lookup table var
  !
  ! Author: Paolo Orsini
  ! Date: 12/01/2017 - mod 04/17/2018
  !
  use Utility_module

  implicit none

  type(lookup_table_var_type) :: lookup_table_var

  if (associated(lookup_table_var%data)) then
    !data is only a pointer to a slice of var_data
    nullify(lookup_table_var%data)
  end if
  call DeallocateArray(lookup_table_var%sample_grad)
  nullify(lookup_table_var%next)

end subroutine LookupTableVarDestroy

! ************************************************************************** !


end module Lookup_Table_module
