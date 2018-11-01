module Connection_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private
  type, public :: connection_set_type
    PetscInt :: id
    PetscInt :: itype                  ! connection type (boundary, internal, source sink
    PetscInt :: num_connections
    PetscInt :: offset
    PetscInt, pointer :: local(:)      ! 1 if connection is local, 0 if connection is ghosted
    PetscInt, pointer :: id_up(:)      ! list of ids of upwind cells
    PetscInt, pointer :: id_dn(:)      ! list of ids of downwind cells
    PetscInt, pointer :: id_up2(:)     ! list of ids of 2nd upwind cells
    PetscInt, pointer :: id_dn2(:)     ! list of ids of 2nd downwind cells
    PetscReal, pointer :: dist(:,:)    ! list of distance vectors, size(-1:3,num_connections) where
                                       !   -1 = fraction upwind
                                       !   0 = magnitude of distance 
                                       !   1-3 = components of unit vector
    PetscReal, pointer :: intercp(:,:) ! x,y,z location of intercept between the line connecting
                                       ! upwind and downwind cells with the face shared by the cells
    PetscReal, pointer :: area(:)      ! list of areas of faces normal to distance vectors
    PetscReal, pointer :: cntr(:,:)    ! coordinates (1:3, num_connections) of the mass center of the face
    PetscInt, pointer :: face_id(:)    ! list of ids of faces (in local order)
    type(connection_set_type), pointer :: next
  end type connection_set_type


  ! pointer data structure required for making an array of region pointers in F90
  type, public :: connection_set_ptr_type
    type(connection_set_type), pointer :: ptr           ! pointer to the connection_set_type
  end type connection_set_ptr_type 
  
  type, public :: connection_set_list_type
    PetscInt :: num_connection_objects
    type(connection_set_type), pointer :: first
    type(connection_set_type), pointer :: last
    type(connection_set_ptr_type), pointer :: array(:)
  end type connection_set_list_type
  
  public :: ConnectionCreate, &
            ConnectionAddToList, &
            ConnectionGetNumberInList, &
            ConnectionInitList, &
            ConnectionCalculateDistances, &
            ConnectionDestroyList, &
            ConnectionDestroy
  
contains

! ************************************************************************** !

function ConnectionCreate(num_connections,connection_itype)
  ! 
  ! Allocates and initializes a new connection
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/07
  ! 

  implicit none
  
  PetscInt :: num_connections
  PetscInt :: connection_itype
  
  type(connection_set_type), pointer :: ConnectionCreate

  type(connection_set_type), pointer :: connection

  allocate(connection)
  connection%id = 0
  connection%itype = connection_itype
  connection%offset = 0
  connection%num_connections = num_connections
  nullify(connection%local)
  nullify(connection%id_up)
  nullify(connection%id_dn)
  nullify(connection%id_up2)
  nullify(connection%id_dn2)
  nullify(connection%face_id)
  nullify(connection%dist)
  nullify(connection%intercp)
  nullify(connection%area)
  nullify(connection%cntr)
  select case(connection_itype)
    case(INTERNAL_CONNECTION_TYPE)
      allocate(connection%id_up(num_connections))
      allocate(connection%id_dn(num_connections))
      allocate(connection%dist(-1:3,num_connections))
      allocate(connection%intercp(1:3,num_connections))
      allocate(connection%area(num_connections))
      allocate(connection%face_id(num_connections))
      connection%id_up = 0
      connection%id_dn = 0
      connection%face_id = 0
      connection%dist = 0.d0
      connection%intercp = 0.d0
      connection%area = 0.d0
    case(BOUNDARY_CONNECTION_TYPE)
      allocate(connection%id_dn(num_connections))
      allocate(connection%dist(-1:3,num_connections))
      allocate(connection%intercp(1:3,num_connections))
      allocate(connection%area(num_connections))
      allocate(connection%face_id(num_connections))
      connection%id_dn = 0
      connection%dist = 0.d0
      connection%intercp = 0.d0
      connection%area = 0.d0
    case(SRC_SINK_CONNECTION_TYPE,INITIAL_CONNECTION_TYPE)
      allocate(connection%id_dn(num_connections))
      connection%id_dn = 0
  end select
  nullify(connection%next)
  
  ConnectionCreate => connection

end function ConnectionCreate

! ************************************************************************** !

function ConnectionGetNumberInList(list)
  ! 
  ! Returns the number of connections in a list
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/19/07
  ! 

  implicit none
  
  type(connection_set_list_type) :: list

  PetscInt :: ConnectionGetNumberInList
  type(connection_set_type), pointer :: cur_connection_set
  
  ConnectionGetNumberInList = 0
  cur_connection_set => list%first
  do
    if (.not.associated(cur_connection_set)) exit
    ConnectionGetNumberInList = ConnectionGetNumberInList + &
                                cur_connection_set%num_connections
    cur_connection_set => cur_connection_set%next
  enddo

end function ConnectionGetNumberInList

! ************************************************************************** !

subroutine ConnectionInitList(list)
  ! 
  ! InitConnectionModule: Initializes module variables, lists, arrays.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/07
  ! 

  implicit none

  type(connection_set_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_connection_objects = 0

end subroutine ConnectionInitList

! ************************************************************************** !

subroutine ConnectionAddToList(new_connection_set,list)
  ! 
  ! Adds a new connection of the module global list of
  ! connections
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/07
  ! 

  implicit none
  
  type(connection_set_type), pointer :: new_connection_set
  type(connection_set_list_type) :: list
  
  list%num_connection_objects = list%num_connection_objects + 1
  new_connection_set%id = list%num_connection_objects
  if (.not.associated(list%first)) list%first => new_connection_set
  if (associated(list%last)) list%last%next => new_connection_set
  list%last => new_connection_set
  
end subroutine ConnectionAddToList

! ************************************************************************** !

subroutine ConnectionConvertListToArray(list)
  ! 
  ! Creates an array of pointers to the
  ! connections in the connection list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/07
  ! 

  implicit none
  
  type(connection_set_list_type) :: list
    
  type(connection_set_type), pointer :: cur_connection_set
  
  
  allocate(list%array(list%num_connection_objects))
  
  cur_connection_set => list%first
  do 
    if (.not.associated(cur_connection_set)) exit
    list%array(cur_connection_set%id)%ptr => cur_connection_set
    cur_connection_set => cur_connection_set%next
  enddo

end subroutine ConnectionConvertListToArray

! ************************************************************************** !

subroutine ConnectionCalculateDistances(dist,gravity,distance_upwind, &
                                        distance_downwind,distance_gravity, &
                                        upwind_weight)
  ! 
  ! Calculates the various distances and weights used in a flux calculation.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  implicit none
  
  PetscReal, intent(in) :: dist(-1:3)
  PetscReal, intent(in) :: gravity(3)
  
  PetscReal, intent(out) :: distance_upwind
  PetscReal, intent(out) :: distance_downwind
  PetscReal, intent(out) :: distance_gravity
  PetscReal, intent(out) :: upwind_weight
  
  ! dist(-1) = scalar - fraction upwind
  ! dist(0) = scalar - magnitude of distance
  ! gravity = vector(3)
  ! dist(1:3) = vector(3) - unit vector
  distance_gravity = dist(0) * &                  ! distance_gravity = dx*g*n
                     dot_product(gravity,dist(1:3))
  distance_upwind = dist(0)*dist(-1)
  distance_downwind = dist(0)-distance_upwind ! should avoid truncation error
  ! upweight could be calculated as 1.d0-fraction_upwind
  ! however, this introduces ever so slight error causing pflow-overhaul not
  ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
  upwind_weight = distance_downwind/(distance_upwind+distance_downwind)
  
end subroutine ConnectionCalculateDistances

! ************************************************************************** !

subroutine ConnectionDestroy(connection)
  ! 
  ! Deallocates a connection
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 
  use Utility_module, only : DeallocateArray
  
  implicit none
  
  type(connection_set_type), pointer :: connection
  
  if (.not.associated(connection)) return
  
  call DeallocateArray(connection%local)
  call DeallocateArray(connection%id_up)
  call DeallocateArray(connection%id_dn)
  call DeallocateArray(connection%id_up2)
  call DeallocateArray(connection%id_dn2)
  call DeallocateArray(connection%face_id)
  call DeallocateArray(connection%dist)
  call DeallocateArray(connection%intercp)
  call DeallocateArray(connection%area)
  call DeallocateArray(connection%cntr)
  
  nullify(connection%next)
  
  deallocate(connection)
  nullify(connection)

end subroutine ConnectionDestroy

! ************************************************************************** !

subroutine ConnectionDestroyList(list)
  ! 
  ! Deallocates the module global list and array of regions
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/07
  ! 

  implicit none
  
  type(connection_set_list_type), pointer :: list
    
  type(connection_set_type), pointer :: cur_connection_set, prev_connection_set
  
  if (.not.associated(list)) return
  
  if (associated(list%array)) deallocate(list%array)
  nullify(list%array)
  
  cur_connection_set => list%first
  do 
    if (.not.associated(cur_connection_set)) exit
    prev_connection_set => cur_connection_set
    cur_connection_set => cur_connection_set%next
    call ConnectionDestroy(prev_connection_set)
  enddo
  
  nullify(list%first)
  nullify(list%last)
  list%num_connection_objects = 0
  
  deallocate(list)
  nullify(list)

end subroutine ConnectionDestroyList

end module Connection_module
