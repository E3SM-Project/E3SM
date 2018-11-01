module Geometry_module
 
#include "petsc/finclude/petscsys.h"
   use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: point3d_type
    PetscInt :: id
    PetscReal :: x
    PetscReal :: y
    PetscReal :: z
  end type point3d_type
  
  type, public :: plane_type
    PetscReal :: A
    PetscReal :: B
    PetscReal :: C
    PetscReal :: D
  end type plane_type
  
  type, public :: polygonal_volume_type
    type(point3d_type), pointer :: xy_coordinates(:)
    type(point3d_type), pointer :: xz_coordinates(:)
    type(point3d_type), pointer :: yz_coordinates(:)
  end type polygonal_volume_type

  interface GeometryPointInPolygon
    module procedure GeometryPointInPolygon1
    module procedure GeometryPointInPolygon2
  end interface GeometryPointInPolygon
  
  interface GeometryComputePlaneWithPoints
    module procedure GeometryComputePlaneWithPoints1
    module procedure GeometryComputePlaneWithPoints2
  end interface GeometryComputePlaneWithPoints
  
  public :: GeometryCreatePolygonalVolume, &
            GeometryReadCoordinates, &
            GeometryReadCoordinate, &
            GeometryPointInPolygonalVolume, &
            GeometryCopyCoordinates, &
            GeometryDestroyPolygonalVolume, &
            GeometryPointInPolygon, &
            GeometryComputePlaneWithPoints, &
            GeomComputePlaneWithGradients, &
            GeometryProjectPointOntoPlane, &
            GeometryGetPlaneIntercept, &
            GeometryGetPlaneZIntercept, &
            GeomGetPlaneGradientinXandY, &
            GeomComputeDistanceFromPlane
!           12345678901234567890123456789012
contains

! ************************************************************************** !

function GeometryCreatePolygonalVolume()
  ! 
  ! Creates a polygonal volume.  I.e. a volume
  ! defined by 3 polygons, on each plane of the
  ! principle coordinate system x,y,z
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/12
  ! 

  implicit none
  
  type(polygonal_volume_type), pointer :: GeometryCreatePolygonalVolume
  
  type(polygonal_volume_type), pointer :: polygonal_volume
  
  allocate(polygonal_volume)
  nullify(polygonal_volume%xy_coordinates)
  nullify(polygonal_volume%xz_coordinates)
  nullify(polygonal_volume%yz_coordinates)
  
  GeometryCreatePolygonalVolume => polygonal_volume

end function GeometryCreatePolygonalVolume

! ************************************************************************** !

subroutine GeometryReadCoordinates(input,option,region_name,coordinates)
  ! 
  ! Reads a list of coordinates
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/12
  ! 

  use Input_Aux_module
  use Option_module

  implicit none
  
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: region_name
  type(point3d_type), pointer :: coordinates(:)
  
  PetscInt :: icount
  PetscInt, parameter :: max_num_coordinates = 100
  type(point3d_type) :: temp_coordinates(max_num_coordinates)

  icount = 0
  do
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'REGION')
    if (InputCheckExit(input,option)) exit              
    icount = icount + 1
    if (icount > max_num_coordinates) then
      write(option%io_buffer, &
            '(''Number of coordinates in region '',a, &
          &'' exceeds limit of '',i3)') region_name, max_num_coordinates
      option%io_buffer = trim(option%io_buffer) // &
        ' Increase size of PetscInt, parameter :: max_num_coordinates ' // &
        ' in RegionReadCoordinates()'
      call printErrMsg(option)
    endif
    call GeometryReadCoordinate(input,option,temp_coordinates(icount),'REGION')
  enddo
  allocate(coordinates(icount))
  do icount = 1, size(coordinates)
    coordinates(icount)%x = temp_coordinates(icount)%x
    coordinates(icount)%y = temp_coordinates(icount)%y
    coordinates(icount)%z = temp_coordinates(icount)%z
  enddo

end subroutine GeometryReadCoordinates

! ************************************************************************** !

subroutine GeometryReadCoordinate(input,option,coordinate,error_string)
  ! 
  ! Reads a list of coordinates
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/12
  ! 

  use Input_Aux_module
  use Option_module

  implicit none
  
  type(input_type) :: input
  type(option_type) :: option
  character(len=*) :: error_string
  type(point3d_type) :: coordinate
  
  call InputReadDouble(input,option,coordinate%x) 
  call InputErrorMsg(input,option,'x-coordinate',error_string)
  call InputReadDouble(input,option,coordinate%y)
  call InputErrorMsg(input,option,'y-coordinate',error_string)
  call InputReadDouble(input,option,coordinate%z)
  call InputErrorMsg(input,option,'z-coordinate',error_string)

end subroutine GeometryReadCoordinate

! ************************************************************************** !

function GeometryPointInPolygonalVolume(x,y,z,polygonal_volume,option)
  ! 
  ! Determines whether a point in xyz space is
  ! within a polygonal volume defined by polygons
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/12
  ! 
 
  use Option_module

  implicit none
  
  PetscReal :: x, y, z
  type(polygonal_volume_type) :: polygonal_volume
  type(option_type) :: option
  
  PetscBool :: xy, xz, yz
  PetscBool :: GeometryPointInPolygonalVolume
    
  GeometryPointInPolygonalVolume = PETSC_FALSE
  
  ! if a polygon is not defined in a particular direction, it is assumed that
  ! the point is within the "infinite" polygon
  
  ! XY plane
  if (associated(polygonal_volume%xy_coordinates)) then
    if (.not.GeometryPointInPolygon(x,y,z,Z_DIRECTION, &
                                  polygonal_volume%xy_coordinates)) then
      return
    endif
  endif

  ! XZ plane
  if (associated(polygonal_volume%xz_coordinates)) then
    if (.not.GeometryPointInPolygon(x,y,z,Y_DIRECTION, &
                                  polygonal_volume%xz_coordinates)) then
      return
    endif
  endif

  ! YZ plane
  if (associated(polygonal_volume%yz_coordinates)) then
    if (.not.GeometryPointInPolygon(x,y,z,X_DIRECTION, &
                                  polygonal_volume%yz_coordinates)) then
      return
    endif
  endif
  
  ! if point is not within any of polygons above, the function will return
  ! prior to this point.
  GeometryPointInPolygonalVolume = PETSC_TRUE
  
end function GeometryPointInPolygonalVolume

! ************************************************************************** !

function GeometryPointInPolygon1(x,y,z,axis,coordinates)
  ! 
  ! Determines whether a point in xyz space is within
  ! a 2d polygon based on coordinate object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/12
  ! 

  implicit none
  
  PetscReal :: x, y, z
  PetscInt :: axis
  type(point3d_type) :: coordinates(:)
  
  PetscInt, parameter :: max_num_coordinates = 100
  PetscReal :: xx, yy
  PetscInt :: i, num_coordinates
  PetscReal :: xx_array(max_num_coordinates), yy_array(max_num_coordinates)
  
  PetscBool :: GeometryPointInPolygon1
  
  GeometryPointInPolygon1 = PETSC_FALSE
  
  num_coordinates = size(coordinates)
  select case(axis)
    case(Z_DIRECTION)
      xx = x
      yy = y
      do i = 1, num_coordinates
        xx_array(i) = coordinates(i)%x
        yy_array(i) = coordinates(i)%y
      enddo
    case(Y_DIRECTION)
      xx = x
      yy = z
      do i = 1, num_coordinates
        xx_array(i) = coordinates(i)%x
        yy_array(i) = coordinates(i)%z
      enddo
    case(X_DIRECTION)
      xx = y
      yy = z
      do i = 1, num_coordinates
        xx_array(i) = coordinates(i)%y
        yy_array(i) = coordinates(i)%z
      enddo
  end select
  if (num_coordinates == 2) then
    GeometryPointInPolygon1 = &
      GeometryPointInRectangle(xx,yy,xx_array,yy_array)
  else
    GeometryPointInPolygon1 = &
      GeometryPointInPolygon2(xx,yy,xx_array,yy_array,num_coordinates)
  endif
  
end function GeometryPointInPolygon1

! ************************************************************************** !

function GeometryPointInPolygon2(x,y,x_array,y_array,num_coordinates)
  ! 
  ! Determines whether a point in xy space is within
  ! a polygon
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/12
  ! 

  implicit none
  
  PetscReal :: x, y
  PetscReal :: x_array(:), y_array(:)
  PetscInt :: num_coordinates
  
  PetscBool :: GeometryPointInPolygon2
  PetscInt :: i, j
  
  GeometryPointInPolygon2 = PETSC_FALSE

  j = 1
  do i = 1, num_coordinates
    j = i + 1
    if (i == num_coordinates) j = 1
    if ((y_array(i) < y .and. y_array(j) >= y) .or. &
        (y_array(j) < y .and. y_array(i) >= y)) then
      if ((x_array(i) + &
           (y-y_array(i))/(y_array(j)-y_array(i))*(x_array(j)-x_array(i))) &
          < x) then
        GeometryPointInPolygon2 = .not.GeometryPointInPolygon2
      endif
    endif    
  enddo

end function GeometryPointInPolygon2

! ************************************************************************** !

function GeometryPointInRectangle(x,y,x_array,y_array)
  ! 
  ! Determines whether a point in xy space is within
  ! a box
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/12
  ! 

  implicit none
  
  PetscReal :: x, y
  PetscReal :: x_array(:), y_array(:)

  PetscBool :: GeometryPointInRectangle
  
  PetscReal :: x1, y1
  PetscReal :: x2, y2
  
  GeometryPointInRectangle = PETSC_FALSE
  
  ! only using first 2 values in each array
  
  x1 = minval(x_array(1:2))
  y1 = minval(y_array(1:2))
  x2 = maxval(x_array(1:2))
  y2 = maxval(y_array(1:2))
  
  
  if (x >= x1 .and. x <= x2 .and. y >= y1 .and. y <= y2) &
    GeometryPointInRectangle = PETSC_TRUE

end function GeometryPointInRectangle

! ************************************************************************** !

subroutine GeometryCopyCoordinates(coordinates_in,coordinates_out)
  ! 
  ! Deallocates a polygonal volume object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/12
  ! 

  implicit none
  
  type(point3d_type) :: coordinates_in(:)
  type(point3d_type), pointer :: coordinates_out(:)
  
  PetscInt :: num_coordinates
  PetscInt :: i
  
  num_coordinates = size(coordinates_in)
  allocate(coordinates_out(num_coordinates))
  do i = 1, num_coordinates
   coordinates_out(i)%x = coordinates_in(i)%x 
   coordinates_out(i)%y = coordinates_in(i)%y
   coordinates_out(i)%z = coordinates_in(i)%z 
  enddo
 
end subroutine GeometryCopyCoordinates

! ************************************************************************** !

subroutine GeometryDestroyPolygonalVolume(polygonal_volume)
  ! 
  ! Deallocates a polygonal volume object
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/09
  ! 

  implicit none
  
  type(polygonal_volume_type), pointer :: polygonal_volume
  
  if (.not.associated(polygonal_volume)) return
  
  if (associated(polygonal_volume%xy_coordinates)) &
    deallocate(polygonal_volume%xy_coordinates)
  nullify(polygonal_volume%xy_coordinates)
  if (associated(polygonal_volume%xz_coordinates)) &
    deallocate(polygonal_volume%xz_coordinates)
  nullify(polygonal_volume%xz_coordinates)
  if (associated(polygonal_volume%yz_coordinates)) &
    deallocate(polygonal_volume%yz_coordinates)
  nullify(polygonal_volume%yz_coordinates)
  
  deallocate(polygonal_volume)
  nullify(polygonal_volume)
  
end subroutine GeometryDestroyPolygonalVolume

! ************************************************************************** !

subroutine GeometryComputePlaneWithPoints1(plane,x1,y1,z1,x2,y2,z2,x3,y3,z3)
  ! 
  ! Calculates the plane defined by a point and gradients in x and y
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/09, 02/01/17
  ! 

  implicit none
  
  type(plane_type) :: plane
  PetscReal :: x1,y1,z1
  PetscReal :: x2,y2,z2
  PetscReal :: x3,y3,z3

  PetscReal :: x12, y12, z12
  PetscReal :: x13, y13, z13

#if 0
  ! more math here
  plane%A = y1*(z2-z3)+y2*(z3-z1)+y3*(z1-z2)
  plane%B = z1*(x2-x3)+z2*(x3-x1)+z3*(x1-x2)
  plane%C = x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)
  plane%D = -1.d0*(x1*(y2*z3-y3*z2)+x2*(y3*z1-y1*z3)+x3*(y1*z2-y2*z1))
#else
  ! less math here
  x12 = x2-x1
  y12 = y2-y1
  z12 = z2-z1
  x13 = x3-x1
  y13 = y3-y1
  z13 = z3-z1
  plane%A = y12*z13-z12*y13
  plane%B = z12*x13-x12*z13
  plane%C = x12*y13-y12*x13
  plane%D = -1.d0*(plane%A*x1+plane%B*y1+plane%C*z1)
#endif
  
end subroutine GeometryComputePlaneWithPoints1

! ************************************************************************** !

subroutine GeometryComputePlaneWithPoints2(plane,point1,point2,point3)
  ! 
  ! Calculates the plane defined by three points
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/09
  ! 

  implicit none
  
  type(plane_type) :: plane
  type(point3d_type) :: point1, point2, point3
  
  call GeometryComputePlaneWithPoints1(plane, &
                                       point1%x,point1%y,point1%z, &
                                       point2%x,point2%y,point2%z, &
                                       point3%x,point3%y,point3%z)

end subroutine GeometryComputePlaneWithPoints2

! ************************************************************************** !

subroutine GeomComputePlaneWithGradients(plane,x,y,z,dz_dx,dz_dy)
  ! 
  ! Calculates the plane defined by a point and gradients in x and y
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/09
  ! 

  implicit none
  
  type(plane_type) :: plane
  PetscReal :: x
  PetscReal :: y
  PetscReal :: z
  PetscReal :: dz_dx
  PetscReal :: dz_dy
  
  PetscReal :: x1, y1, z1, x2, y2, z2, x3, y3, z3
  
  x1 = x
  y1 = y
  z1 = z
  x2 = x + 1.d0
  y2 = y
  z2 = z + dz_dx
  x3 = x
  y3 = y + 1.d0
  z3 = z + dz_dy
  
  call GeometryComputePlaneWithPoints(plane,x1,y1,z1,x2,y2,z2,x3,y3,z3)

end subroutine GeomComputePlaneWithGradients

! ************************************************************************** !

subroutine GeometryProjectPointOntoPlane(plane,point,intercept)
  ! 
  ! Calculates the intercept of a point with a plane
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/22/11
  ! 

  implicit none
  
  type(plane_type) :: plane
  type(point3d_type) :: point
  type(point3d_type) :: intercept

  PetscReal :: scalar
  
  ! plane equation:
  !   A*x + B*y + C*z + D = 0

  scalar = (plane%A*point%x + plane%B*point%y + plane%C*point%z + plane%D) / &
           (plane%A*plane%A + plane%B*plane%B + plane%C*plane%C)
  
  intercept%x = point%x - plane%A * scalar
  intercept%y = point%y - plane%B * scalar
  intercept%z = point%z - plane%C * scalar

end subroutine GeometryProjectPointOntoPlane

! ************************************************************************** !

subroutine GeometryGetPlaneIntercept(plane,point1,point2,intercept)
  ! 
  ! Calculates the intercept of a line with a plane
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/09
  ! 

  implicit none
  
  type(plane_type) :: plane
  type(point3d_type) :: point1, point2
  type(point3d_type) :: intercept

  PetscReal :: x1,y1,z1
  PetscReal :: x2,y2,z2
  PetscReal :: u
    
  x1 = point1%x
  y1 = point1%y
  z1 = point1%z
  x2 = point2%x
  y2 = point2%y
  z2 = point2%z
 
 
  u = (plane%A*x1 + plane%B*y1 + plane%C*z1 + plane%D) / &
      (plane%A*(x1-x2) + plane%B*(y1-y2) + plane%C*(z1-z2))

  intercept%x = point1%x + u*(point2%x-point1%x)
  intercept%y = point1%y + u*(point2%y-point1%y)
  intercept%z = point1%z + u*(point2%z-point1%z)

end subroutine GeometryGetPlaneIntercept

! ************************************************************************** !

function GeometryGetPlaneZIntercept(plane,x,y)
  ! 
  ! Calculates the intercept of a line with an x,y coordinate
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/01/17
  ! 

  implicit none
  
  type(plane_type) :: plane
  PetscReal :: x
  PetscReal :: y

  PetscReal :: GeometryGetPlaneZIntercept

  ! Ax + By + Cz + D = 0
  ! z = (Ax + By + D)/-C
    
  GeometryGetPlaneZIntercept = -1.d0 * &
                               (x * plane%A + y * plane%B + plane%D) / &
                               plane%C

end function GeometryGetPlaneZIntercept

! ************************************************************************** !

subroutine GeomGetPlaneGradientinXandY(plane,dz_dx,dz_dy)
  ! 
  ! Calculates the intercept of a line with a plane
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/09
  ! 

  implicit none
  
  type(plane_type) :: plane
  PetscReal :: dz_dx
  PetscReal :: dz_dy
    
  dz_dx = -1.d0*plane%A/plane%C
  dz_dy = -1.d0*plane%B/plane%C

end subroutine GeomGetPlaneGradientinXandY

! ************************************************************************** !

function GeomComputeDistanceFromPlane(plane,point)
  ! 
  ! Calculates the distance of a point from a plane
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/11
  ! 

  implicit none
  
  type(plane_type) :: plane
  type(point3d_type) :: point
  
  PetscReal :: GeomComputeDistanceFromPlane

  GeomComputeDistanceFromPlane = &
    (plane%A*point%x + plane%B*point%y + plane%C*point%z + plane%D) / &
    sqrt(plane%A*plane%A+plane%B*plane%B+plane%C*plane%C)
  
end function GeomComputeDistanceFromPlane

end module Geometry_module
