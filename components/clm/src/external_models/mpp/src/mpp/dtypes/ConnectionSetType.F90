module ConnectionSetType

#ifdef USE_PETSC_LIB

  use ArrayDimThree, only : array_dim3_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private
#include <petsc/finclude/petsc.h>

  type, public :: connection_type
     private
     PetscInt                       :: type         ! horizontal or vertical
     PetscInt                       :: id_up        ! IDs of upwind cells [-]
     PetscInt                       :: id_dn        ! IDs of downwind cells [-]
     PetscReal                      :: area         ! area normal to connection [m^2]

     ! Distances are computed as : downwind - upwind
     PetscReal                      :: dist_up      ! Magnitude of distance from centroid of upwind cell to centroid of cell face [m]
     PetscReal                      :: dist_dn      ! Magnitude of distance from centroid of downwind cell to centroid of cell face [m]

     type(array_dim3_type), pointer :: dist_unitvec ! Unit vector [-]
   contains
     procedure, public :: Init            => ConnInit
     procedure, public :: SetType         => ConnSetType
     procedure, public :: SetIDUp         => ConnSetIDUp
     procedure, public :: SetIDDn         => ConnSetIDDn
     procedure, public :: SetArea         => ConnSetArea
     procedure, public :: SetDistUp       => ConnSetDistUp
     procedure, public :: SetDistDn       => ConnSetDistDn
     procedure, public :: SetDistUnitVec  => ConnSetDistUnitVec
     procedure, public :: GetType         => ConnGetType
     procedure, public :: GetIDUp         => ConnGetIDUp
     procedure, public :: GetIDDn         => ConnGetIDDn
     procedure, public :: GetArea         => ConnGetArea
     procedure, public :: GetDistUp       => ConnGetDistUp
     procedure, public :: GetDistDn       => ConnGetDistDn
     procedure, public :: GetDistUnitVec  => ConnGetDistUnitVec
     procedure, public :: GetDistUnitVecX => ConnGetDistUnitVecX
     procedure, public :: GetDistUnitVecY => ConnGetDistUnitVecY
     procedure, public :: GetDistUnitVecZ => ConnGetDistUnitVecZ
  end type connection_type

  type, public :: connection_set_type

     PetscInt                           :: id              ! identifier
     PetscInt                           :: num_connections ! total num. of connections
     type(connection_type)    , pointer :: conn(:)         ! information about all connections
     type(connection_set_type), pointer :: next

  end type connection_set_type

  type, public :: connection_set_list_type
     PetscInt                           :: num_connection_list
     type(connection_set_type), pointer :: first
     type(connection_set_type), pointer :: last
  end type connection_set_list_type

  public :: ConnectionSetNew
  public :: ConnectionSetDestroy
  public :: ConnectionSetListInit
  public :: ConnectionSetListAddSet
  public :: ConnectionSetListClean
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine ConnInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize a connection
    !
    ! !USES:
    use MultiPhysicsProbConstants   , only : CONN_VERTICAL
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this

    this%type    = CONN_VERTICAL
    this%id_up   = 0
    this%id_dn   = 0
    this%area    = 0.d0
    this%dist_up = 0.d0
    this%dist_dn = 0.d0
    allocate(this%dist_unitvec)
    
  end subroutine ConnInit

  !------------------------------------------------------------------------
  subroutine ConnSetType(this, val)
    !
    ! !DESCRIPTION:
    ! Set type of connection
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this
    PetscInt               :: val

    this%type = val
    
  end subroutine ConnSetType

  !------------------------------------------------------------------------
  subroutine ConnSetIDUp(this, val)
    !
    ! !DESCRIPTION:
    ! Set ID of upwind cell
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this
    PetscInt               :: val

    this%id_up = val
    
  end subroutine ConnSetIDUp

  !------------------------------------------------------------------------
  subroutine ConnSetIDDn(this, val)
    !
    ! !DESCRIPTION:
    ! Set ID of downwind cell
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this
    PetscInt               :: val

    this%id_dn = val
    
  end subroutine ConnSetIDDn

  !------------------------------------------------------------------------
  subroutine ConnSetArea(this, val)
    !
    ! !DESCRIPTION:
    ! Set area of connection
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this
    PetscReal              :: val

    this%area = val
    
  end subroutine ConnSetArea

  !------------------------------------------------------------------------
  subroutine ConnSetDistUp(this, val)
    !
    ! !DESCRIPTION:
    ! Set upwind distance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this
    PetscReal              :: val

    this%dist_up = val
    
  end subroutine ConnSetDistUp

  !------------------------------------------------------------------------
  subroutine ConnSetDistDn(this, val)
    !
    ! !DESCRIPTION:
    ! Set downwind distance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this
    PetscReal              :: val

    this%dist_dn = val
    
  end subroutine ConnSetDistDn

  !------------------------------------------------------------------------
  subroutine ConnSetDistUnitVec(this, u_x, u_y, u_z)
    !
    ! !DESCRIPTION:
    ! Set downwind distance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this
    PetscReal              :: u_x, u_y, u_z

    this%dist_unitvec%arr(1) = u_x
    this%dist_unitvec%arr(2) = u_y
    this%dist_unitvec%arr(3) = u_z
    
  end subroutine ConnSetDistUnitVec

  !------------------------------------------------------------------------
  function ConnGetIDUp(this)
    !
    ! !DESCRIPTION:
    ! Return ID of upwind cell
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this
    !
    PetscInt               :: ConnGetIDUp

    ConnGetIDUp = this%id_up
    
  end function ConnGetIDUp

  !------------------------------------------------------------------------
  function ConnGetIDDn(this)
    !
    ! !DESCRIPTION:
    ! Return ID of downwind cell
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this
    !
    PetscInt               :: ConnGetIDDn

    ConnGetIDDn = this%id_dn
    
  end function ConnGetIDDn

  !------------------------------------------------------------------------
  function ConnGetType(this)
    !
    ! !DESCRIPTION:
    ! Return connection type
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this
    !
    PetscInt               :: ConnGetType

    ConnGetType = this%type
    
  end function ConnGetType

  !------------------------------------------------------------------------
  function ConnGetArea(this)
    !
    ! !DESCRIPTION:
    ! Return area of connection
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this
    !
    PetscReal              :: ConnGetArea

    ConnGetArea = this%area
    
  end function ConnGetArea

  !------------------------------------------------------------------------
  function ConnGetDistUp(this)
    !
    ! !DESCRIPTION:
    ! Return upwind distance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this
    !
    PetscReal              :: ConnGetDistUp

    ConnGetDistUp = this%dist_up
    
  end function ConnGetDistUp

  !------------------------------------------------------------------------
  function ConnGetDistDn(this)
    !
    ! !DESCRIPTION:
    ! Return downwind distance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this
    !
    PetscReal              :: ConnGetDistDn

    ConnGetDistDn = this%dist_dn
    
  end function ConnGetDistDn

  !------------------------------------------------------------------------
  function ConnGetDistUnitVec(this)
    !
    ! !DESCRIPTION:
    ! Return all components of unit distance vector
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this
    !
    PetscReal              :: ConnGetDistUnitVec(3)

    ConnGetDistUnitVec(1) = this%GetDistUnitVecX()
    ConnGetDistUnitVec(2) = this%GetDistUnitVecY()
    ConnGetDistUnitVec(3) = this%GetDistUnitVecZ()
    
  end function ConnGetDistUnitVec

  !------------------------------------------------------------------------
  function ConnGetDistUnitVecX(this)
    !
    ! !DESCRIPTION:
    ! Return X component of unit vector
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this
    !
    PetscReal              :: ConnGetDistUnitVecX

    ConnGetDistUnitVecX = this%dist_unitvec%arr(1)
    
  end function ConnGetDistUnitVecX

  !------------------------------------------------------------------------
  function ConnGetDistUnitVecY(this)
    !
    ! !DESCRIPTION:
    ! Return Y component of unit vector
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this
    !
    PetscReal              :: ConnGetDistUnitVecY

    ConnGetDistUnitVecY = this%dist_unitvec%arr(2)
    
  end function ConnGetDistUnitVecY

  !------------------------------------------------------------------------
  function ConnGetDistUnitVecZ(this)
    !
    ! !DESCRIPTION:
    ! Return Z component of unit vector
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(connection_type) :: this
    !
    PetscReal              :: ConnGetDistUnitVecZ

    ConnGetDistUnitVecZ = this%dist_unitvec%arr(3)
    
  end function ConnGetDistUnitVecZ

  !------------------------------------------------------------------------
  function ConnectionSetNew(num_connections)
    !
    ! !DESCRIPTION:
    ! Return a new connection set
    !
    ! !USES:
    use MultiPhysicsProbConstants   , only : CONN_VERTICAL
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscInt, intent(in)              :: num_connections
    !
    ! !LOCAL VARIABLES:
    type(connection_set_type),pointer :: ConnectionSetNew
    type(connection_set_type),pointer :: conn_set
    PetscInt                          :: iconn

    allocate(conn_set)

    conn_set%id                = 0
    conn_set%num_connections   = num_connections

    allocate(conn_set%conn(num_connections));
    do iconn = 1, num_connections
       call conn_set%conn(iconn)%Init()
    end do

    nullify(conn_set%next)

    ConnectionSetNew => conn_set

  end function ConnectionSetNew

  !------------------------------------------------------------------------
  subroutine ConnectionSetListInit(list)
    !
    ! !DESCRIPTION:
    ! Initialize a connection set list
    !
    implicit none
    !
    ! !ARGUMENTS
    type(connection_set_list_type) :: list

    list%num_connection_list   = 0

    nullify(list%first)
    nullify(list%last )

  end subroutine ConnectionSetListInit

  !------------------------------------------------------------------------
  subroutine ConnectionSetListAddSet(list, new_conn_set)
    !
    ! !DESCRIPTION:
    ! Add a connection set to a connection set list.
    !
    implicit none
    !
    ! !ARGUMENTS
    type(connection_set_list_type)    :: list
    type(connection_set_type),pointer :: new_conn_set

    list%num_connection_list = list%num_connection_list + 1
    new_conn_set%id          = list%num_connection_list

    if (.not.associated(list%first)) then
       list%first => new_conn_set
    endif

    if (associated(list%last)) then
       list%last%next => new_conn_set
    endif

    list%last => new_conn_set

  end subroutine ConnectionSetListAddSet

  !------------------------------------------------------------------------
  subroutine ConnectionSetDestroy(conn_set)
    !
    ! !DESCRIPTION:
    ! Release all allocated memory
    !
    implicit none
    !
    ! !ARGUMENTS
    type(connection_set_type),pointer :: conn_set

    if (.not.(associated(conn_set))) return

    deallocate(conn_set%conn)
    nullify(conn_set%conn)

    nullify(conn_set%next)
    deallocate(conn_set)
    nullify(conn_set)

  end subroutine ConnectionSetDestroy

  !------------------------------------------------------------------------
  subroutine ConnectionSetListClean(list)
    !
    ! !DESCRIPTION:
    ! Release all allocated memory
    !
    implicit none
    !
    ! !ARGUMENTS
    type(connection_set_list_type), intent(inout) :: list
    !
    ! !LOCAL VARIABLES:
    type(connection_set_type), pointer            :: curr_conn_set, prev_conn_set

    curr_conn_set => list%first
    do
       if (.not.associated(curr_conn_set)) exit
       prev_conn_set => curr_conn_set
       curr_conn_set => curr_conn_set%next
       call ConnectionSetDestroy(prev_conn_set)
    enddo

    call ConnectionSetListInit(list)

  end subroutine ConnectionSetListClean

#endif

end module ConnectionSetType
