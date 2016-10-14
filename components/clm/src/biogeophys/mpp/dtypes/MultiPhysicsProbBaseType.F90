
module MultiPhysicsProbBaseType

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Base object for multi-physics problem
  !-----------------------------------------------------------------------
  ! !USES:
  use MeshType        , only : mesh_type
  use mpp_varctl      , only : iulog
  use mpp_abortutils  , only : endrun
  use mpp_shr_log_mod , only : errMsg => shr_log_errMsg
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
     procedure, public :: SetName                    => MPPSetName
     procedure, public :: SetID                      => MPPSetID
     procedure, public :: SetNumMeshes               => MPPSetNumMeshes
     procedure, public :: MeshSetName                => MPPMeshSetName
     procedure, public :: MeshSetOrientation         => MPPMeshSetOrientation
     procedure, public :: MeshSetID                  => MPPMeshSetID
     procedure, public :: MeshSetDimensions          => MPPMeshSetDimensions
     procedure, public :: MeshSetGeometricAttributes => MPPMeshSetGeometricAttributes
     procedure, public :: MeshSetGridCellFilter      => MPPMeshSetGridCellFilter
     procedure, public :: MeshComputeVolume          => MPPMeshComputeVolume
     procedure, public :: MeshSetConnectionSet       => MPPMeshSetConnectionSet
  end type multiphysicsprob_base_type

  public :: MPPBaseInit
  public :: MMPBaseClean
  public :: MPPSetName
  public :: MPPSetNumMeshes

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
  subroutine MPPSetName(this, name)
    !
    ! !DESCRIPTION:
    ! Set name of MPP
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    character(len =*)                 :: name

    this%name = trim(name)

  end subroutine MPPSetName

  !------------------------------------------------------------------------
  subroutine MPPSetNumMeshes(this, nmesh)
    !
    ! !DESCRIPTION:
    ! Set number of meshes
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    PetscInt                          :: nmesh
    !
    ! !LOCAL VARIABLES:
    PetscInt                          :: imesh

    if (nmesh <= 0) then
       write(iulog,*) 'Attempting to set invalid value for number of meshes.' // &
            ' nmesh = ',nmesh
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (this%nmesh /= 0) then
       write(iulog,*) 'Number of meshes has already been set.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    this%nmesh = nmesh
    allocate(this%meshes(this%nmesh))

    do imesh = 1, nmesh
       call this%meshes(imesh)%Init()
    enddo

  end subroutine MPPSetNumMeshes

  !------------------------------------------------------------------------
  subroutine MPPSetID(this, id)
    !
    ! !DESCRIPTION:
    ! Set ID for MPP
    !
    use MultiPhysicsProbConstants, only : MPP_VSFM_SNES_CLM
    use MultiPhysicsProbConstants, only : MPP_THERMAL_TBASED_KSP_CLM
    use MultiPhysicsProbConstants, only : MPP_THERMAL_EBASED_SNES_CLM
    use MultiPhysicsProbConstants, only : MPP_TH_SNES_CLM
    use MultiPhysicsProbConstants, only : PETSC_SNES
    use MultiPhysicsProbConstants, only : PETSC_KSP
    use MultiPhysicsProbConstants, only : SOE_RE_ODE
    use MultiPhysicsProbConstants, only : SOE_THERMAL_TBASED
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    PetscInt                          :: id

    select case (id)
    case (MPP_VSFM_SNES_CLM, MPP_THERMAL_EBASED_SNES_CLM, MPP_TH_SNES_CLM)
       this%id                    = id
       this%solver_type           = PETSC_SNES

    case (MPP_THERMAL_TBASED_KSP_CLM)
       this%id                    = id
       this%solver_type           = PETSC_KSP

    case default
       write(iulog,*) 'Attempting to set unsupported MPP ID. id = ', id
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine MPPSetID

  !------------------------------------------------------------------------
  subroutine MPPMeshSetName(this, imesh, mesh_name)
    !
    ! !DESCRIPTION:
    ! Set name of mesh
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    PetscInt                          :: imesh
    character(len =*)                 :: mesh_name

    call CheckMeshIndex(this, imesh)
    call this%meshes(imesh)%SetName(mesh_name)

  end subroutine MPPMeshSetName

  !------------------------------------------------------------------------
  subroutine MPPMeshSetOrientation(this, imesh, orientation)
    !
    ! !DESCRIPTION:
    ! Set mesh orientation
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    PetscInt                          :: imesh
    PetscInt                          :: orientation

    call CheckMeshIndex(this, imesh)
    call this%meshes(imesh)%SetOrientation(orientation)

  end subroutine MPPMeshSetOrientation

  !------------------------------------------------------------------------
  subroutine MPPMeshSetID(this, imesh, id)
    !
    ! !DESCRIPTION:
    ! Set mesh ID
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    PetscInt                          :: imesh
    PetscInt                          :: id

    call CheckMeshIndex(this, imesh)
    call this%meshes(imesh)%SetID(id)

  end subroutine MPPMeshSetID

  !------------------------------------------------------------------------
  subroutine MPPMeshSetDimensions(this, imesh, ncells_local, ncells_ghost, nlev)
    !
    ! !DESCRIPTION:
    ! Set mesh dimension
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    PetscInt                          :: imesh
    PetscInt                          :: ncells_local
    PetscInt                          :: ncells_ghost
    PetscInt                          :: nlev

    call CheckMeshIndex(this, imesh)
    call this%meshes(imesh)%SetDimensions(ncells_local, ncells_ghost, nlev)

  end subroutine MPPMeshSetDimensions

  !------------------------------------------------------------------------
  subroutine MPPMeshSetGeometricAttributes(this, imesh, var_id, values)
    !
    ! !DESCRIPTION:
    ! Set mesh geometric attributes
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    PetscInt                          :: imesh
    PetscInt                          :: var_id
    PetscReal, pointer                :: values(:)

    call CheckMeshIndex(this, imesh)
    call this%meshes(imesh)%SetGeometricAttributes(var_id, values)

  end subroutine MPPMeshSetGeometricAttributes

  !------------------------------------------------------------------------
  subroutine MPPMeshSetGridCellFilter(this, imesh, values)
    !
    ! !DESCRIPTION:
    ! Set grid filter
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    PetscInt                          :: imesh
    PetscInt, pointer                 :: values(:)

    call CheckMeshIndex(this, imesh)
    call this%meshes(imesh)%SetGridCellFilter( values)

  end subroutine MPPMeshSetGridCellFilter

  !------------------------------------------------------------------------
  subroutine MPPMeshComputeVolume(this, imesh)
    !
    ! !DESCRIPTION:
    ! Compute volume for cells of a mesh
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    PetscInt                          :: imesh

    call CheckMeshIndex(this, imesh)
    call this%meshes(imesh)%ComputeVolume()

  end subroutine MPPMeshComputeVolume

  !------------------------------------------------------------------------
  subroutine MPPMeshSetConnectionSet(this, imesh, conn_type, nconn, id_up, id_dn, &
       dist_up, dist_dn, area)
    !
    ! !DESCRIPTION:
    ! Set connection set for a mesh
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    PetscInt                          :: imesh
    PetscInt                          :: conn_type
    PetscInt                          :: nconn
    PetscInt, pointer                 :: id_up(:)
    PetscInt, pointer                 :: id_dn(:)
    PetscReal, pointer                :: dist_up(:)
    PetscReal, pointer                :: dist_dn(:)
    PetscReal, pointer                :: area(:)

    call CheckMeshIndex(this, imesh)
    call this%meshes(imesh)%SetConnectionSet(conn_type, nconn, id_up, id_dn, &
         dist_up, dist_dn, area)

  end subroutine MPPMeshSetConnectionSet

  !------------------------------------------------------------------------
  subroutine CheckMeshIndex(this, imesh)
    !
    ! !DESCRIPTION:
    ! Check mesh index is valid
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    PetscInt                          :: imesh

    if (imesh <= 0) then
       write(iulog,*) 'imesh is an invalid non-positive number.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (imesh > this%nmesh) then
       write(iulog,*) 'imesh (=', imesh, ') > this%nmesh (=',this%nmesh,')'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

  end subroutine CheckMeshIndex

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

#endif

end module MultiPhysicsProbBaseType
