#ifdef USE_PETSC_LIB


module MultiPhysicsProbBaseType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Base object for multi-physics problem
  !-----------------------------------------------------------------------
  ! !USES:
  use MeshType, only : mesh_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include "finclude/petscsys.h"

  type, public :: multiphysicsprob_base_type
     character(len =256)      :: name        ! name of the multi-physics problem (MPP)
     PetscInt                 :: id          ! identifier for the MPP
     class(mesh_type),pointer :: meshes(:)   ! meshes associated with the MPP
     PetscInt                 :: nmesh       ! number of meshes in the MPP
     PetscInt                 :: solver_type ! identifier for the type of PETSc solver
   contains
     procedure, public :: Init
     procedure, public :: Clean
  end type multiphysicsprob_base_type

  public :: MPPBaseInit
  public :: MMPBaseClean

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine Init(this)
    !
    ! !DESCRIPTION:
    ! Wrapper subroutine for initialization
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this

    call MPPBaseInit(this)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine MPPBaseInit(this)
    !
    ! !DESCRIPTION:
    ! Initializes a MPP object
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this

    this%name         = ""
    this%id           = 0
    this%nmesh        = 0
    this%solver_type  = 0

    nullify(this%meshes)

  end subroutine MPPBaseInit

  !------------------------------------------------------------------------
  subroutine Clean(this)
    !
    ! !DESCRIPTION:
    ! Wrapper subroutine for freeing up memory
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this

    call MMPBaseClean(this)

  end subroutine Clean

  !------------------------------------------------------------------------
  subroutine MMPBaseClean(this)
    !
    ! !DESCRIPTION:
    ! Free up memory associated with meshes
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    PetscInt                          :: imesh

    do imesh = 1, this%nmesh
       call this%meshes(imesh)%Clean()
    enddo

  end subroutine MMPBaseClean

end module MultiPhysicsProbBaseType

#endif
