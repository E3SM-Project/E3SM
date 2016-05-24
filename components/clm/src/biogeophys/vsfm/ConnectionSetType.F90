#ifdef USE_PETSC_LIB


module ConnectionSetType

  use ArrayDimThree, only : array_dim3_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private
#include "finclude/petscsys.h"

  type, public :: connection_set_type

     PetscInt                                     :: id              ! identifier

     PetscInt                                     :: num_connections ! total num. of connections
     PetscInt, pointer                            :: id_up(:)        ! IDs of upwind cells [-]
     PetscInt, pointer                            :: id_dn(:)        ! IDs of downwind cells [-]

     PetscReal, pointer                           :: area(:)         ! area normal to connection [m^2]

                                                                     ! Distances are computed as : downwind - upwind
     PetscReal, pointer                           :: dist_up(:)      ! Magnitude of distance from centroid of
                                                                     ! upwind cell to centroid of cell face [m]
     PetscReal, pointer                           :: dist_dn(:)      ! Magnitude of distance from centroid of
                                                                     ! downwind cell to centroid of cell face [m]
     type(array_dim3_type), dimension(:), pointer :: dist_unitvec    ! Unit vector [-]

     type(connection_set_type), pointer           :: next

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
  function ConnectionSetNew(num_connections)
    !
    ! !DESCRIPTION:
    ! Return a new connection set
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    implicit none
    !
    ! !ARGUMENTS
    PetscInt, intent(in)              :: num_connections
    !
    ! !LOCAL VARIABLES:
    type(connection_set_type),pointer :: ConnectionSetNew
    type(connection_set_type),pointer :: conn_set

    allocate(conn_set)

    conn_set%id                = 0
    conn_set%num_connections   = num_connections

    allocate(conn_set%id_up       (num_connections)); conn_set%id_up  (:) = 0
    allocate(conn_set%id_dn       (num_connections)); conn_set%id_dn  (:) = 0
    allocate(conn_set%area        (num_connections)); conn_set%area   (:) = nan
    allocate(conn_set%dist_up     (num_connections)); conn_set%dist_up(:) = nan
    allocate(conn_set%dist_dn     (num_connections)); conn_set%dist_dn(:) = nan
    allocate(conn_set%dist_unitvec(num_connections));

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

    if (associated(conn_set%id_up        )) deallocate(conn_set%id_up)
    if (associated(conn_set%id_dn        )) deallocate(conn_set%id_dn)
    if (associated(conn_set%area         )) deallocate(conn_set%area)
    if (associated(conn_set%dist_up      )) deallocate(conn_set%dist_up)
    if (associated(conn_set%dist_dn      )) deallocate(conn_set%dist_dn)
    if (associated(conn_set%dist_unitvec )) deallocate(conn_set%dist_unitvec)

    nullify(conn_set%id_up        )
    nullify(conn_set%id_up        )
    nullify(conn_set%area         )
    nullify(conn_set%dist_up      )
    nullify(conn_set%dist_dn      )
    nullify(conn_set%dist_unitvec )
    nullify(conn_set%next         )

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

end module ConnectionSetType

#endif
