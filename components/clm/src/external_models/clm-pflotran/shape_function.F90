module Shape_Function_module

  use Gauss_module
  use Grid_Unstructured_Cell_module
  use PFLOTRAN_Constants_module

  implicit none

#include "petsc/finclude/petscsys.h"

  private
  
  type, public :: shapefunction_type
    PetscInt :: EleType               ! element type
    PetscReal, pointer :: zeta(:)     ! coordinates of point in reference element
    PetscReal, pointer :: N(:)        ! shape function for all nodes evaluated at zeta (N is a row vector)
    PetscReal, pointer :: DN(:,:)     ! derivatives of shape function with respect to zeta
    PetscReal, pointer :: coord(:,:)  ! local coordinates of the vertices in the reference element
  end type shapefunction_type
    
  public :: ShapeFunctionInitialize, ShapeFunctionCalculate, &
            ShapeFunctionDestroy
  
  contains

! ************************************************************************** !

subroutine ShapeFunctionInitialize(shapefunction)
  ! 
  ! Allocate memory for shapefunction type
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  ! 
 
  type(shapefunction_type) :: shapefunction
  PetscReal, pointer :: coord(:,:)
  
  select case(shapefunction%EleType)
    case(LINE_TYPE)
      allocate(shapefunction%N(TWO_INTEGER))
      allocate(shapefunction%DN(TWO_INTEGER,ONE_INTEGER)) 
      allocate(shapefunction%zeta(ONE_INTEGER))   
      allocate(shapefunction%coord(TWO_INTEGER,ONE_INTEGER))
    case(QUAD_TYPE)
      allocate(shapefunction%N(FOUR_INTEGER))
      allocate(shapefunction%DN(FOUR_INTEGER,TWO_INTEGER))
      allocate(shapefunction%zeta(TWO_INTEGER))
      allocate(shapefunction%coord(FOUR_INTEGER,TWO_INTEGER))
    case(WEDGE_TYPE)
      allocate(shapefunction%N(SIX_INTEGER))
      allocate(shapefunction%DN(SIX_INTEGER,THREE_INTEGER))
      allocate(shapefunction%zeta(THREE_INTEGER))  
      allocate(shapefunction%coord(SIX_INTEGER,THREE_INTEGER))
    case(TET_TYPE)
      allocate(shapefunction%N(FOUR_INTEGER))
      allocate(shapefunction%DN(FOUR_INTEGER,THREE_INTEGER))
      allocate(shapefunction%zeta(THREE_INTEGER)) 
      allocate(shapefunction%coord(FOUR_INTEGER,THREE_INTEGER))
    case(PYR_TYPE)
      allocate(shapefunction%N(FIVE_INTEGER))
      allocate(shapefunction%DN(FIVE_INTEGER,THREE_INTEGER))
      allocate(shapefunction%zeta(THREE_INTEGER)) 
      allocate(shapefunction%coord(FIVE_INTEGER,THREE_INTEGER))
    case(HEX_TYPE)
      allocate(shapefunction%N(EIGHT_INTEGER))
      allocate(shapefunction%DN(EIGHT_INTEGER,THREE_INTEGER))
      allocate(shapefunction%zeta(THREE_INTEGER))
      allocate(shapefunction%coord(EIGHT_INTEGER,THREE_INTEGER))
    case default
      print *, 'Error: Invalid EleType. Enter an EleType from L2,' // &
               ' T3, Q4, B8 only.'       
  end select   
  
  shapefunction%zeta = 0.d0
  shapefunction%N  = 0.d0
  shapefunction%DN = 0.d0
  shapefunction%coord = 0.d0
  
  coord => shapefunction%coord  
  
  select case(shapefunction%EleType)
    case(LINE_TYPE)
      coord(1,1) = -1.d0
      coord(2,1) = 1.d0
    case(QUAD_TYPE)
      coord(1,1) = -1.d0
      coord(1,2) = -1.d0
      coord(2,1) = 1.d0
      coord(2,2) = -1.d0
      coord(3,1) = 1.d0
      coord(3,2) = 1.d0
      coord(4,1) = -1.d0
      coord(4,2) = 1.d0
    case(WEDGE_TYPE)
      coord(1,:) = -(/0.d0,1.d0,1.d0/)
      coord(2,:) = -(/-1.d0,-1.d0,1.d0/)
      coord(3,:) = -(/1.d0,-1.d0,1.d0/)
      coord(4,:) = -(/0.d0,1.d0,-1.d0/)
      coord(5,:) = -(/-1.d0,-1.d0,-1.d0/)
      coord(6,:) = -(/1.d0,-1.d0,-1.d0/)
    case(TET_TYPE)
      coord(1,:) = -(/0.d0,1.d0,1.d0/)
      coord(2,:) = -(/-1.d0,-1.d0,1.d0/)
      coord(3,:) = -(/1.d0,-1.d0,1.d0/)
      coord(4,:) = -(/0.d0,0.d0,-1.d0/)
    case(PYR_TYPE)
      coord(1,:) = -(/1.d0,1.d0,1.d0/)
      coord(2,:) = -(/-1.d0,1.d0,1.d0/)
      coord(3,:) = -(/-1.d0,-1.d0,1.d0/)
      coord(4,:) = -(/1.d0,-1.d0,1.d0/)
      coord(5,:) = -(/0.d0,0.d0,-1.d0/)
    case(HEX_TYPE)
      coord(1,:) = -(/1.d0,1.d0,1.d0/)
      coord(2,:) = -(/-1.d0,1.d0,1.d0/)
      coord(3,:) = -(/-1.d0,-1.d0,1.d0/)
      coord(4,:) = -(/1.d0,-1.d0,1.d0/)
      coord(5,:) = -(/1.d0,1.d0,-1.d0/)
      coord(6,:) = -(/-1.d0,1.d0,-1.d0/)
      coord(7,:) = -(/-1.d0,-1.d0,-1.d0/)
      coord(8,:) = -(/1.d0,-1.d0,-1.d0/)
    case default
      print *, 'Error: Invalid EleType. Enter an EleType from L2,' // &
               ' T3, Q4, B8 only.'       
  end select     
   
end subroutine ShapeFunctionInitialize

! ************************************************************************** !

subroutine ShapeFunctionCalculate(shapefunction)
  ! 
  ! Subroutine provides shape functions and its derivatives
  ! at a given spatial point (in the reference element) 'zeta' for various finite
  ! elements.
  ! Input variables
  ! EleType: element type
  ! L2: one-dimensional element -1 <= x <= +1
  ! Q4: four node quadrilateral element -1 <= x, y <= +1
  ! Q9: nine node quadrilateral element
  ! T3: three node right-angled triangle with h = b = 1
  ! T6: six node right-angled triangle
  ! B8: eight node right-angled tetrahedron element
  ! zeta: coordinates of a point in the reference element
  ! Output variables
  ! N: shape functions for all nodes evaluated at zeta (N is a row
  ! vector!)
  ! DN: derivatives of shape functions with respect to zeta
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  ! 

  type(shapefunction_type) :: shapefunction
  PetscReal, pointer :: zeta(:)
  PetscReal, pointer :: N(:)
  PetscReal, pointer :: DN(:,:)
  PetscInt :: i
    
  N => shapefunction%N
  DN => shapefunction%DN
  zeta => shapefunction%zeta
  
  do i = 1, size(zeta)
    if (zeta(i) < -1.d0 .or. zeta(i) > 1.d0) then
      print *, 'Error: Enter values between -1 and 1 for zeta'
      stop
    endif
  enddo
  
  select case(shapefunction%EleType)
    case(LINE_TYPE)
      N(1) = 0.5d0*(1.d0 - zeta(1))
      N(2) = 0.5d0*(1.d0 + zeta(1))
      DN(1,1) = -0.5d0
      DN(2,1) = 0.5d0
    case(QUAD_TYPE)
      N(1) = 1.d0/4.d0*(1.d0 - zeta(1))*(1.d0 - zeta(2))
      N(2) = 1.d0/4.d0*(1.d0 + zeta(1))*(1.d0 - zeta(2))
      N(3) = 1.d0/4.d0*(1.d0 + zeta(1))*(1.d0 + zeta(2))
      N(4) = 1.d0/4.d0*(1.d0 - zeta(1))*(1.d0 + zeta(2))
      DN(1,1) = -1.d0/4.d0*(1.d0 - zeta(2))
      DN(1,2) = -1.d0/4.d0*(1.d0 - zeta(1))
      DN(2,1) = 1.d0/4.d0*(1.d0 - zeta(2)) 
      DN(2,2) = -1.d0/4.d0*(1.d0 + zeta(1))
      DN(3,1) = 1.d0/4.d0*(1.d0 + zeta(2))  
      DN(3,2) = 1.d0/4.d0*(1.d0 + zeta(1))
      DN(4,1) = -1.d0/4.d0*(1.d0 + zeta(2))  
      DN(4,2) = 1.d0/4.d0*(1.d0 - zeta(1))
    case(WEDGE_TYPE)
      N(1) = 1.d0/4.d0*(1.d0 - zeta(2))*(1.d0 - zeta(3))
      N(2) = 1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 + zeta(2))*(1.d0 - zeta(3))
      N(3) = 1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 + zeta(2))*(1.d0 - zeta(3))
      N(4) = 1.d0/4.d0*(1.d0 - zeta(2))*(1.d0 + zeta(3))
      N(5) = 1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 + zeta(2))*(1.d0 + zeta(3))
      N(6) = 1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 + zeta(2))*(1.d0 + zeta(3))
         
      DN(1,:) = (/0.d0, &
                  1.d0/4.d0*(-1.d0)*(1.d0 - zeta(3)), &
                  1.d0/4.d0*(1.d0 - zeta(2))*(-1.d0)/) 
      DN(2,:) = (/1.d0/8.d0*(+1.d0)*(1.d0 + zeta(2))*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(+1.d0)*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 + zeta(2))*(-1.d0)/)
      DN(3,:) = (/1.d0/8.d0*(-1.d0)*(1.d0 + zeta(2))*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(+1.d0)*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 + zeta(2))*(-1.d0)/)
      DN(4,:) = (/0.d0, &
                  1.d0/4.d0*(-1.d0)*(1.d0 + zeta(3)), &
                  1.d0/4.d0*(1.d0 - zeta(2))*(+1.d0)/)
      DN(5,:) = (/1.d0/8.d0*(+1.d0)*(1.d0 + zeta(2))*(1.d0 + zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(+1.d0)*(1.d0 + zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 + zeta(2))*(+1.d0)/)
      DN(6,:) = (/1.d0/8.d0*(-1.d0)*(1.d0 + zeta(2))*(1.d0 + zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(+1.d0)*(1.d0 + zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 + zeta(2))*(+1.d0)/) 
      
    case(TET_TYPE)
      N(1) = 1.d0/4.d0*(1.d0 - zeta(2))*(1.d0 - zeta(3))
      N(2) = 1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 + zeta(2))*(1.d0 - zeta(3))
      N(3) = 1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 + zeta(2))*(1.d0 - zeta(3))
      N(4) = 1.d0/2.d0*(1.d0 + zeta(3))
         
      DN(1,:) = (/0.d0, &
                  1.d0/4.d0*(-1.d0)*(1.d0 - zeta(3)), &
                  1.d0/4.d0*(1.d0 - zeta(2))*(-1.d0)/) 
      DN(2,:) = (/1.d0/8.d0*(+1.d0)*(1.d0 + zeta(2))*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(+1.d0)*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 + zeta(2))*(-1.d0)/)
      DN(3,:) = (/1.d0/8.d0*(-1.d0)*(1.d0 + zeta(2))*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(+1.d0)*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 + zeta(2))*(-1.d0)/)
      DN(4,:) = (/0.d0,0.d0,1.d0/2.d0/)
      
    case(PYR_TYPE)
      N(1) = 1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 - zeta(2))*(1.d0 - zeta(3))
      N(2) = 1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 - zeta(2))*(1.d0 - zeta(3))
      N(3) = 1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 + zeta(2))*(1.d0 - zeta(3))
      N(4) = 1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 + zeta(2))*(1.d0 - zeta(3))
      N(5) = 1.d0/2.d0*(1.d0 + zeta(3))
         
      DN(1,:) = (/1.d0/8.d0*(-1.d0)*(1.d0 - zeta(2))*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(-1.d0)*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 - zeta(2))*(-1.d0)/) 
      DN(2,:) = (/1.d0/8.d0*(+1.d0)*(1.d0 - zeta(2))*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(-1.d0)*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 - zeta(2))*(-1.d0)/)
      DN(3,:) = (/1.d0/8.d0*(+1.d0)*(1.d0 + zeta(2))*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(+1.d0)*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 + zeta(2))*(-1.d0)/)
      DN(4,:) = (/1.d0/8.d0*(-1.d0)*(1.d0 + zeta(2))*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(+1.d0)*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 + zeta(2))*(-1.d0)/)
      DN(5,:) = (/0.d0,0.d0,1.d0/2.d0/)
      
    case(HEX_TYPE)
      N(1) = 1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 - zeta(2))*(1.d0 - zeta(3))
      N(2) = 1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 - zeta(2))*(1.d0 - zeta(3))
      N(3) = 1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 + zeta(2))*(1.d0 - zeta(3))
      N(4) = 1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 + zeta(2))*(1.d0 - zeta(3))
      N(5) = 1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 - zeta(2))*(1.d0 + zeta(3))
      N(6) = 1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 - zeta(2))*(1.d0 + zeta(3))
      N(7) = 1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 + zeta(2))*(1.d0 + zeta(3))
      N(8) = 1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 + zeta(2))*(1.d0 + zeta(3))
         
      DN(1,:) = (/1.d0/8.d0*(-1.d0)*(1.d0 - zeta(2))*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(-1.d0)*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 - zeta(2))*(-1.d0)/) 
      DN(2,:) = (/1.d0/8.d0*(+1.d0)*(1.d0 - zeta(2))*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(-1.d0)*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 - zeta(2))*(-1.d0)/)
      DN(3,:) = (/1.d0/8.d0*(+1.d0)*(1.d0 + zeta(2))*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(+1.d0)*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 + zeta(2))*(-1.d0)/)
      DN(4,:) = (/1.d0/8.d0*(-1.d0)*(1.d0 + zeta(2))*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(+1.d0)*(1.d0 - zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 + zeta(2))*(-1.d0)/)
      DN(5,:) = (/1.d0/8.d0*(-1.d0)*(1.d0 - zeta(2))*(1.d0 + zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(-1.d0)*(1.d0 + zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 - zeta(2))*(+1.d0)/)
      DN(6,:) = (/1.d0/8.d0*(+1.d0)*(1.d0 - zeta(2))*(1.d0 + zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(-1.d0)*(1.d0 + zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 - zeta(2))*(+1.d0)/)
      DN(7,:) = (/1.d0/8.d0*(+1.d0)*(1.d0 + zeta(2))*(1.d0 + zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(+1.d0)*(1.d0 + zeta(3)), &
                  1.d0/8.d0*(1.d0 + zeta(1))*(1.d0 + zeta(2))*(+1.d0)/)
      DN(8,:) = (/1.d0/8.d0*(-1.d0)*(1.d0 + zeta(2))*(1.d0 + zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(+1.d0)*(1.d0 + zeta(3)), &
                  1.d0/8.d0*(1.d0 - zeta(1))*(1.d0 + zeta(2))*(+1.d0)/)
    case default
      print *, 'Error: Invalid EleType. Enter an EleType from L2,' // &
               ' T3, Q4, B8 only.'
  end select
      
end subroutine ShapeFunctionCalculate

! ************************************************************************** !

subroutine ShapeFunctionDestroy(shapefunction)
  ! 
  ! Dellocate memory for shapefunction type
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  ! 

  type(shapefunction_type) :: shapefunction
  
  deallocate(shapefunction%N)
  nullify(shapefunction%N)
  deallocate(shapefunction%DN)
  nullify(shapefunction%DN)
  deallocate(shapefunction%zeta)
  nullify(shapefunction%zeta)
  deallocate(shapefunction%coord)
  nullify(shapefunction%coord)
  
end subroutine ShapeFunctionDestroy
     
end module Shape_Function_module

