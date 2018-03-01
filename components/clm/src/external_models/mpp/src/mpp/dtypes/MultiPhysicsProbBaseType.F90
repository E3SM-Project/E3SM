
module MultiPhysicsProbBaseType

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Base object for multi-physics problem
  !-----------------------------------------------------------------------

#include <petsc/finclude/petsc.h>
  use petscsys
  use petscvec
  use petscmat
  use petscts
  use petscsnes
  use petscdm
  use petscdmda
  !
  ! !USES:
  use MeshType        , only : mesh_type
  use mpp_varctl      , only : iulog
  use mpp_abortutils  , only : endrun
  use mpp_shr_log_mod , only : errMsg => shr_log_errMsg
  use SystemOfEquationsBaseType
  use SystemOfEquationsBasePointerType   , only : sysofeqns_base_pointer_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: multiphysicsprob_base_type
     character(len =256)                         :: name          ! name of the multi-physics problem (MPP)
     PetscInt                                    :: id            ! identifier for the MPP
     class(mesh_type)                  , pointer :: meshes(:)     ! meshes associated with the MPP
     PetscInt                                    :: nmesh         ! number of meshes in the MPP
     PetscInt                                    :: solver_type   ! identifier for the type of PETSc solver
     class(sysofeqns_base_type)        , pointer :: soe           ! system-of-equations for the MPP
     type(sysofeqns_base_pointer_type) , pointer :: soe_ptr ! workaround to send user defined context to PETSc
   contains
     procedure, public :: Init
     procedure, public :: Clean
     procedure, public :: SetName                       => MPPSetName
     procedure, public :: SetID                         => MPPSetID
     procedure, public :: SetNumMeshes                  => MPPSetNumMeshes
     procedure, public :: MeshSetName                   => MPPMeshSetName
     procedure, public :: MeshSetOrientation            => MPPMeshSetOrientation
     procedure, public :: MeshSetID                     => MPPMeshSetID
     procedure, public :: MeshSetDimensions             => MPPMeshSetDimensions
     procedure, public :: MeshSetGeometricAttributes    => MPPMeshSetGeometricAttributes
     procedure, public :: MeshSetGridCellFilter         => MPPMeshSetGridCellFilter
     procedure, public :: MeshComputeVolume             => MPPMeshComputeVolume
     procedure, public :: MeshSetConnectionSet          => MPPMeshSetConnectionSet
     procedure, public :: AddGovEqn                     => MPPAddGovEqn
     procedure, public :: AddGovEqnWithMeshRank         => MPPAddGovEqnWithMeshRank
     procedure, public :: SetMeshesOfGoveqns            => MPPSetMeshesOfGoveqns
     procedure, public :: SetMeshesOfGoveqnsByMeshRank  => MPPSetMeshesOfGoveqnsByMeshRank
     procedure, public :: SetMPIRank                    => MPPSetMPIRank
     procedure, public :: GetMPIRank                    => MPPGetMPIRank
     procedure, public :: GovEqnUpdateBCConnectionSet   => MPPGovEqnUpdateBCConnectionSet
     procedure, public :: GovEqnSetCouplingVars         => MPPGovEqnSetCouplingVars
     procedure, public :: GovEqnSetInternalCouplingVars => MPPGovEqnSetInternalCouplingVars
     procedure, public :: GovEqnSetBothCouplingVars     => MPPGovEqnSetBothCouplingVars
     procedure, public :: GovEqnAddCouplingCondition    => MPPGovEqnAddCouplingCondition
  end type multiphysicsprob_base_type

  public :: MPPBaseInit
  public :: MMPBaseClean
  public :: MPPSetName
  public :: MPPSetNumMeshes
  public :: MPPSetupProblem

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

    nullify(this%meshes  )
    nullify(this%soe     )
    nullify(this%soe_ptr )

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
       dist_up, dist_dn, area, itype)
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
    PetscInt, pointer                 :: itype(:)

    call CheckMeshIndex(this, imesh)
    call this%meshes(imesh)%SetConnectionSet(conn_type, nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype)

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

  !------------------------------------------------------------------------
  subroutine MPPAddGovEqn(this, geq_type, name, mesh_itype)
    !
    ! !DESCRIPTION:
    ! Adds a governing equation to the MPP
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    PetscInt                          :: geq_type
    character(len =*)                 :: name
    PetscInt                          :: mesh_itype

    call this%soe%AddGovEqn(geq_type, name, mesh_itype)

  end subroutine MPPAddGovEqn

  !------------------------------------------------------------------------
  subroutine MPPAddGovEqnWithMeshRank(this, geq_type, name, mesh_rank)
    !
    ! !DESCRIPTION:
    ! Adds a governing equation to the MPP
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    PetscInt                          :: geq_type
    character(len =*)                 :: name
    PetscInt                          :: mesh_rank

    call this%soe%AddGovEqnWithMeshRank(geq_type, name, mesh_rank)

  end subroutine MPPAddGovEqnWithMeshRank

  !------------------------------------------------------------------------
  subroutine MPPSetMeshesOfGoveqns(this)
    !
    ! !DESCRIPTION:
    ! Set association of governing equations and meshes
    !
    use GoverningEquationBaseType, only : goveqn_base_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this

    call this%soe%SetMeshesOfGoveqns(this%meshes, this%nmesh)

  end subroutine MPPSetMeshesOfGoveqns

  !------------------------------------------------------------------------
  subroutine MPPSetMeshesOfGoveqnsByMeshRank(this)
    !
    ! !DESCRIPTION:
    ! Set association of governing equations and meshes
    !
    use GoverningEquationBaseType, only : goveqn_base_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this

    call this%soe%SetMeshesOfGoveqnsByMeshRank(this%meshes, this%nmesh)

  end subroutine MPPSetMeshesOfGoveqnsByMeshRank

  !------------------------------------------------------------------------
  subroutine MPPSetMPIRank(this, rank)
    !
    ! !DESCRIPTION:
    ! Sets MPI rank
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    PetscInt                          :: rank

    if (associated(this%soe)) then
       this%soe%mpi_rank = rank
    endif

  end subroutine MPPSetMPIRank

  !------------------------------------------------------------------------
  subroutine MPPGetMPIRank(this, rank)
    !
    ! !DESCRIPTION:
    ! Returns MPI rank
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    PetscInt, intent(out)             :: rank

    rank = -1
    if (associated(this%soe)) then
       rank = this%soe%mpi_rank
    endif

  end subroutine MPPGetMPIRank

  !------------------------------------------------------------------------
  subroutine MPPGovEqnUpdateBCConnectionSet(this, igoveqn, icond, &
       var_type, nval, values)
    !
    ! !DESCRIPTION:
    ! For a boundary condition of a given governing equation, update distance
    ! for a downstream cell.
    !
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use GoverningEquationBaseType , only : goveqn_base_type
    use MultiPhysicsProbConstants , only : VAR_DIST_DN
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type)   :: this
    PetscInt                            :: igoveqn
    PetscInt                            :: icond
    PetscInt                            :: nval
    PetscInt                            :: var_type
    PetscReal                 , pointer :: values (:)
    !
    class(goveqn_base_type)   , pointer :: cur_goveq
    type(condition_type)      , pointer :: cur_cond
    type(connection_set_type) , pointer :: cur_conn_set
    PetscInt                            :: ii
    PetscInt                            :: iconn
    PetscInt                            :: bc_idx
    PetscBool                           :: bc_found

    if (igoveqn > this%soe%ngoveqns) then
       write(iulog,*) 'Attempting to access governing equation that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq => this%soe%goveqns
    do ii = 1, igoveqn-1
       cur_goveq => cur_goveq%next
    enddo

    bc_found = PETSC_FALSE
    bc_idx = 0
    cur_cond => cur_goveq%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit

       bc_idx = bc_idx + 1
       if (bc_idx == icond) then
          bc_found = PETSC_TRUE

          cur_conn_set => cur_cond%conn_set
          if (nval /= cur_conn_set%num_connections) then
             write(iulog,*) 'Number of values to update connections ' // &
                  'do not match number of connections.'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          endif

          do iconn = 1, cur_conn_set%num_connections

             select case(var_type)
             case (VAR_DIST_DN)
                call cur_conn_set%conn(iconn)%SetDistDn(values(iconn))
             case default
                write(iulog,*) 'Unknown variable type'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end select
          enddo

          exit

       end if

       cur_cond => cur_cond%next
    enddo

    if (.not.bc_found) then
       write(iulog,*) 'Failed to find icond = ',icond,' in the boundary condition list.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

  end subroutine MPPGovEqnUpdateBCConnectionSet

  !------------------------------------------------------------------------
  subroutine MPPGovEqnSetCouplingVars(this, igoveqn, nvars, &
       var_ids, goveqn_ids)
    !
    ! !DESCRIPTION:
    ! In order to couple the given governing equation, add:
    ! - ids of variables needed, and
    ! - ids of governing equations from which variables are needed.
    ! needed for coupling
    ! 
    ! !USES:
    use ConditionType             , only : condition_type
    use GoverningEquationBaseType , only : goveqn_base_type
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use CouplingVariableType      , only : coupling_variable_type
    use CouplingVariableType      , only : CouplingVariableCreate
    use CouplingVariableType      , only : CouplingVariableListAddCouplingVar
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type)      :: this
    PetscInt                               :: igoveqn
    PetscInt                               :: nvars
    PetscInt                     , pointer :: var_ids(:)
    PetscInt                     , pointer :: goveqn_ids(:)
    !
    class(goveqn_base_type)      , pointer :: cur_goveq_1
    class(goveqn_base_type)      , pointer :: cur_goveq_2
    type(condition_type)         , pointer :: cur_cond_1
    type(condition_type)         , pointer :: cur_cond_2
    type(coupling_variable_type) , pointer :: cpl_var
    PetscInt                               :: ii
    PetscInt                               :: ieqn
    PetscInt                               :: ivar
    PetscInt                               :: bc_idx_1
    PetscInt                               :: bc_idx_2
    PetscInt                               :: bc_offset_1
    PetscBool                              :: bc_found

    if (igoveqn > this%soe%ngoveqns) then
       write(iulog,*) 'Attempting to set coupling vars for governing ' // &
            'equation that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq_1 => this%soe%goveqns
    do ii = 1, igoveqn-1
       cur_goveq_1 => cur_goveq_1%next
    end do

    do ivar = 1, nvars

       if (goveqn_ids(ivar) > this%soe%ngoveqns) then
          write(iulog,*) 'Attempting to set coupling vars to a governing ' // &
               'equation that is not in the list'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       bc_found    = PETSC_FALSE
       bc_idx_1    = 1
       bc_offset_1 = 0

       cur_cond_1 => cur_goveq_1%boundary_conditions%first
       do
          if (.not.associated(cur_cond_1)) exit

          ! Is this the appropriate BC?
          if (cur_cond_1%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then
             do ieqn = 1, cur_cond_1%num_other_goveqs
                if (cur_cond_1%list_id_of_other_goveqs(ieqn) == goveqn_ids(ivar) ) then
                   bc_found = PETSC_TRUE
                   exit
                endif
             enddo
          endif

          if (bc_found) exit

          bc_idx_1    = bc_idx_1    + 1
          bc_offset_1 = bc_offset_1 + cur_cond_1%conn_set%num_connections

          cur_cond_1 => cur_cond_1%next
       enddo

       if (.not.bc_found) then
          write(iulog,*)'For goveqn%name = ',trim(cur_goveq_1%name) // &
               ', no coupling boundary condition found to copule it with ' // &
               'equation_number = ', goveqn_ids(ivar)
       endif

       cur_goveq_2 => this%soe%goveqns
       do ii = 1, goveqn_ids(ivar)-1
          cur_goveq_2 => cur_goveq_2%next
       enddo

       bc_found    = PETSC_FALSE
       bc_idx_2    = 1

       cur_cond_2 => cur_goveq_2%boundary_conditions%first
       do
          if (.not.associated(cur_cond_2)) exit

          ! Is this the appropriate BC?
          if (cur_cond_2%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then
             do ieqn = 1, cur_cond_2%num_other_goveqs
                if (cur_cond_2%list_id_of_other_goveqs(ieqn) == igoveqn ) then
                   bc_found = PETSC_TRUE
                   exit
                endif
             enddo
          endif

          if (bc_found) exit

          bc_idx_2    = bc_idx_2    + 1

          cur_cond_2 => cur_cond_2%next
       enddo

       if (.not.bc_found) then
          write(iulog,*)'For goveqn%name = ',trim(cur_goveq_2%name) // &
               ', no coupling boundary condition found to copule it with ' // &
               'equation_number = ', bc_idx_2
       endif

       cpl_var => CouplingVariableCreate()

       cpl_var%variable_type                     = var_ids(ivar)
       cpl_var%num_cells                         = cur_cond_1%conn_set%num_connections
       cpl_var%rank_of_coupling_goveqn           = goveqn_ids(ivar)
       cpl_var%variable_is_bc_in_coupling_goveqn = PETSC_TRUE
       cpl_var%offset_of_bc_in_current_goveqn    = bc_offset_1
       cpl_var%rank_of_bc_in_current_goveqn      = bc_idx_1
       cpl_var%rank_of_bc_in_coupling_goveqn     = bc_idx_2

       call CouplingVariableListAddCouplingVar(cur_goveq_1%coupling_vars, cpl_var)

    enddo

  end subroutine MPPGovEqnSetCouplingVars
  
  !------------------------------------------------------------------------
  subroutine MPPGovEqnSetInternalCouplingVars(this, igoveqn, nvars, &
       var_ids, goveqn_ids)
    !
    ! !DESCRIPTION:
    ! In order to couple the given governing equation, add:
    ! - ids of variables needed, and
    ! - ids of governing equations from which variables are needed.
    ! needed for coupling
    ! 
    ! !USES:
    use ConditionType             , only : condition_type
    use GoverningEquationBaseType , only : goveqn_base_type
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use CouplingVariableType      , only : coupling_variable_type
    use CouplingVariableType      , only : CouplingVariableCreate
    use CouplingVariableType      , only : CouplingVariableListAddCouplingVar
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type)      :: this
    PetscInt                               :: igoveqn
    PetscInt                               :: nvars
    PetscInt                     , pointer :: var_ids(:)
    PetscInt                     , pointer :: goveqn_ids(:)
    !
    class(goveqn_base_type)      , pointer :: cur_goveq_1
    class(goveqn_base_type)      , pointer :: cur_goveq_2
    type(condition_type)         , pointer :: cur_cond_1
    type(condition_type)         , pointer :: cur_cond_2
    type(coupling_variable_type) , pointer :: cpl_var
    PetscInt                               :: ii
    PetscInt                               :: ieqn
    PetscInt                               :: ivar
    PetscInt                               :: bc_idx_1
    PetscInt                               :: bc_idx_2
    PetscInt                               :: bc_offset_1
    PetscBool                              :: bc_found

    if (igoveqn > this%soe%ngoveqns) then
       write(iulog,*) 'Attempting to set coupling vars for governing ' // &
            'equation that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq_1 => this%soe%goveqns
    do ii = 1, igoveqn-1
       cur_goveq_1 => cur_goveq_1%next
    end do

    do ivar = 1, nvars

       if (goveqn_ids(ivar) > this%soe%ngoveqns) then
          write(iulog,*) 'Attempting to set coupling vars to a governing ' // &
               'equation that is not in the list'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif
       
       cpl_var => CouplingVariableCreate()

       cpl_var%variable_type                     = var_ids(ivar)
       cpl_var%rank_of_coupling_goveqn           = goveqn_ids(ivar)
       cpl_var%variable_is_bc_in_coupling_goveqn = PETSC_FALSE

       call CouplingVariableListAddCouplingVar(cur_goveq_1%coupling_vars, cpl_var)

    enddo

  end subroutine MPPGovEqnSetInternalCouplingVars

  !------------------------------------------------------------------------
  subroutine MPPGovEqnSetBothCouplingVars(this, igoveqn, nvars, &
       var_ids, goveqn_ids, is_bc)
    !
    ! !DESCRIPTION:
    ! In order to couple the given governing equation, add:
    ! - ids of variables needed, and
    ! - ids of governing equations from which variables are needed.
    ! needed for coupling
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use GoverningEquationBaseType , only : goveqn_base_type
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use CouplingVariableType      , only : coupling_variable_type
    use CouplingVariableType      , only : CouplingVariableCreate
    use CouplingVariableType      , only : CouplingVariableListAddCouplingVar
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type)      :: this
    PetscInt                               :: igoveqn
    PetscInt                               :: nvars
    PetscInt                     , pointer :: var_ids(:)
    PetscInt                     , pointer :: goveqn_ids(:)
    PetscInt                     , pointer :: is_bc(:)
    !
    class(goveqn_base_type)      , pointer :: cur_goveq_1
    class(goveqn_base_type)      , pointer :: cur_goveq_2
    type(condition_type)         , pointer :: cur_cond_1
    type(condition_type)         , pointer :: cur_cond_2
    type(coupling_variable_type) , pointer :: cpl_var
    PetscInt                               :: ii
    PetscInt                               :: ieqn
    PetscInt                               :: ivar
    PetscInt                               :: bc_idx_1
    PetscInt                               :: bc_idx_2
    PetscInt                               :: bc_offset_1
    PetscBool                              :: bc_found

    if (igoveqn > this%soe%ngoveqns) then
       write(iulog,*) 'Attempting to set coupling vars for governing ' // &
            'equation that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq_1 => this%soe%goveqns
    do ii = 1, igoveqn-1
       cur_goveq_1 => cur_goveq_1%next
    end do

    do ivar = 1, nvars

       if (goveqn_ids(ivar) > this%soe%ngoveqns) then
          write(iulog,*) 'Attempting to set coupling vars to a governing ' // &
               'equation that is not in the list'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       if (is_bc(ivar) == 1) then
          bc_found    = PETSC_FALSE
          bc_idx_1    = 1
          bc_offset_1 = 0

          cur_cond_1 => cur_goveq_1%boundary_conditions%first
          do
             if (.not.associated(cur_cond_1)) exit

             ! Is this the appropriate BC?
             if (cur_cond_1%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then
                do ieqn = 1, cur_cond_1%num_other_goveqs
                   if (cur_cond_1%list_id_of_other_goveqs(ieqn) == goveqn_ids(ivar) ) then
                      bc_found = PETSC_TRUE
                      exit
                   endif
                enddo
             endif

             if (bc_found) exit

             bc_idx_1    = bc_idx_1    + 1
             bc_offset_1 = bc_offset_1 + cur_cond_1%conn_set%num_connections

             cur_cond_1 => cur_cond_1%next
          enddo

          if (.not.bc_found) then
             write(iulog,*)'For goveqn%name = ',trim(cur_goveq_1%name) // &
                  ', no coupling boundary condition found to copule it with ' // &
                  'equation_number = ', goveqn_ids(ivar)
          endif

          cur_goveq_2 => this%soe%goveqns
          do ii = 1, goveqn_ids(ivar)-1
             cur_goveq_2 => cur_goveq_2%next
          enddo

          bc_found    = PETSC_FALSE
          bc_idx_2    = 1

          cur_cond_2 => cur_goveq_2%boundary_conditions%first
          do
             if (.not.associated(cur_cond_2)) exit

             ! Is this the appropriate BC?
             if (cur_cond_2%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then
                do ieqn = 1, cur_cond_2%num_other_goveqs
                   if (cur_cond_2%list_id_of_other_goveqs(ieqn) == igoveqn ) then
                      bc_found = PETSC_TRUE
                      exit
                   endif
                enddo
             endif

             if (bc_found) exit

             bc_idx_2    = bc_idx_2    + 1

             cur_cond_2 => cur_cond_2%next
          enddo

          if (.not.bc_found) then
             write(iulog,*)'For goveqn%name = ',trim(cur_goveq_2%name) // &
                  ', no coupling boundary condition found to copule it with ' // &
                  'equation_number = ', bc_idx_2
          endif

          cpl_var => CouplingVariableCreate()

          cpl_var%variable_type                     = var_ids(ivar)
          cpl_var%num_cells                         = cur_cond_1%conn_set%num_connections
          cpl_var%rank_of_coupling_goveqn           = goveqn_ids(ivar)
          cpl_var%variable_is_bc_in_coupling_goveqn = PETSC_TRUE
          cpl_var%offset_of_bc_in_current_goveqn    = bc_offset_1
          cpl_var%rank_of_bc_in_current_goveqn      = bc_idx_1
          cpl_var%rank_of_bc_in_coupling_goveqn     = bc_idx_2

          call CouplingVariableListAddCouplingVar(cur_goveq_1%coupling_vars, cpl_var)

    else

          cpl_var => CouplingVariableCreate()

          cpl_var%variable_type                     = var_ids(ivar)
          cpl_var%rank_of_coupling_goveqn           = goveqn_ids(ivar)
          cpl_var%variable_is_bc_in_coupling_goveqn = PETSC_FALSE

          call CouplingVariableListAddCouplingVar(cur_goveq_1%coupling_vars, cpl_var)

       endif

    enddo

  end subroutine MPPGovEqnSetBothCouplingVars

  !------------------------------------------------------------------------
  subroutine MPPGovEqnAddCouplingCondition(this, ieqn_1, ieqn_2, &
       iregion_1, iregion_2)
    !
    ! !DESCRIPTION:
    ! Adds a boundary condition to couple ieqn_1 and ieqn_2
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type) :: this
    PetscInt                          :: ieqn_1
    PetscInt                          :: ieqn_2
    PetscInt                          :: iregion_1
    PetscInt                          :: iregion_2
    !
    character(len=256)                :: name
    PetscInt                          :: num_other_goveqs
    PetscInt, pointer                 :: id_of_other_goveqs(:)

    num_other_goveqs = 1
    allocate(id_of_other_goveqs(num_other_goveqs))

    write(name,*) ieqn_2
    name = 'BC_for_coupling_with_equation_' // trim(adjustl(name))
    id_of_other_goveqs(1) = ieqn_2
    call this%soe%AddCouplingBCsInGovEqn(ieqn_1, &
         name, '[K]', iregion_1, num_other_goveqs, id_of_other_goveqs)

    write(name,*) ieqn_1
    name = 'BC_for_coupling_with_equation_' // trim(adjustl(name))
    id_of_other_goveqs(1) = ieqn_1
    call this%soe%AddCouplingBCsInGovEqn(ieqn_2,  &
         name, '[K]', iregion_2, num_other_goveqs, id_of_other_goveqs)

    deallocate(id_of_other_goveqs)

  end subroutine MPPGovEqnAddCouplingCondition

  !------------------------------------------------------------------------
  subroutine MPPSetupProblem(this)
    !
    ! !DESCRIPTION:
    ! Sets the PETSc SNES problem
    !
    use GoverningEquationBaseType , only : goveqn_base_type
    use MultiPhysicsProbConstants , only : PETSC_SNES
    use MultiPhysicsProbConstants , only : PETSC_KSP
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type)      :: this

    select case(this%solver_type)
    case (PETSC_SNES)
       call MPPSetupProblemSNES(this)
    case (PETSC_KSP)
       call MPPSetupProblemKSP(this)
    case default
       write(iulog,*) 'VSFMMPPSetup: Unknown this%solver_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    call this%soe%CreateVectorsForGovEqn()

    this%soe%solver%petsc_solver_type = this%solver_type
    
  end subroutine MPPSetupProblem

  !------------------------------------------------------------------------
  subroutine MPPSetupProblemSNES(this)
    !
    ! !DESCRIPTION:
    ! Sets the PETSc SNES problem
    !
    use GoverningEquationBaseType        , only : goveqn_base_type
    use SystemOfEquationsBasePointerType , only : SOEResidual
    use SystemOfEquationsBasePointerType , only : SOEJacobian
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type)      :: this
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type)    , pointer   :: cur_goveq
    class(sysofeqns_base_type) , pointer   :: base_soe
    PetscInt                               :: size
    PetscInt                               :: igoveq
    PetscErrorCode                         :: ierr
    DM                         , pointer   :: dms(:)
    PetscReal                  , parameter :: atol    = PETSC_DEFAULT_REAL
    PetscReal                  , parameter :: rtol    = PETSC_DEFAULT_REAL
    PetscReal                  , parameter :: stol    = 1.d-10
    PetscInt                   , parameter :: max_it  = PETSC_DEFAULT_INTEGER
    PetscInt                   , parameter :: max_f   = PETSC_DEFAULT_INTEGER
    character(len=256)                     :: name

    base_soe         => this%soe
    this%soe_ptr%ptr => this%soe
    
    ! Create PETSc DM for each governing equation

    allocate(dms(base_soe%ngoveqns))

    igoveq = 0
    cur_goveq => base_soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       igoveq = igoveq + 1
       size   = cur_goveq%mesh%ncells_local

       call DMDACreate1d(PETSC_COMM_SELF, &
            DM_BOUNDARY_NONE, size, 1, 1, &
            PETSC_NULL_INTEGER, dms(igoveq), ierr);
       CHKERRQ(ierr)

       write(name,*) igoveq
       name = 'fgoveq_' // trim(adjustl(name))
       call DMSetOptionsPrefix(dms(igoveq), name, ierr); CHKERRQ(ierr)

       call DMSetFromOptions(dms(igoveq), ierr); CHKERRQ(ierr)
       call DMSetUp(         dms(igoveq) , ierr); CHKERRQ(ierr)

       write(name,*) igoveq
       name = 'goveq_' // trim(adjustl(name))
       call DMDASetFieldName(dms(igoveq), 0, name, ierr); CHKERRQ(ierr)

       cur_goveq => cur_goveq%next
    enddo

    ! DM-Composite approach

    ! Create DMComposite: temperature
    call DMCompositeCreate(PETSC_COMM_SELF, base_soe%solver%dm, ierr); CHKERRQ(ierr)
    call DMSetOptionsPrefix(base_soe%solver%dm, "temperature_", ierr); CHKERRQ(ierr)

    ! Add DMs to DMComposite
    do igoveq = 1, base_soe%ngoveqns
       call DMCompositeAddDM(base_soe%solver%dm, dms(igoveq), ierr); CHKERRQ(ierr)
    enddo

    ! Setup DM
    call DMSetUp(base_soe%solver%dm, ierr); CHKERRQ(ierr)

    ! Create matrix
    call DMCreateMatrix    (base_soe%solver%dm   , base_soe%solver%Amat, ierr); CHKERRQ(ierr)

    call MatSetOption      (base_soe%solver%Amat , MAT_NEW_NONZERO_LOCATION_ERR , &
         PETSC_FALSE, ierr); CHKERRQ(ierr)
    call MatSetOption      (base_soe%solver%Amat , MAT_NEW_NONZERO_ALLOCATION_ERR, &
         PETSC_FALSE, ierr); CHKERRQ(ierr)

    call MatSetFromOptions (base_soe%solver%Amat , ierr); CHKERRQ(ierr)

    call DMCreateMatrix     (base_soe%solver%dm , base_soe%solver%jac, ierr ); CHKERRQ(ierr)
    call MatSetOption       (base_soe%solver%jac, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)
    call MatSetOption       (base_soe%solver%jac, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)

    ! Create vectors
    call DMCreateGlobalVector(base_soe%solver%dm, base_soe%solver%soln         , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(base_soe%solver%dm, base_soe%solver%rhs          , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(base_soe%solver%dm, base_soe%solver%soln_prev    , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(base_soe%solver%dm, base_soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(base_soe%solver%dm, base_soe%solver%res          , ierr); CHKERRQ(ierr)

    ! Initialize vectors
    call VecZeroEntries(base_soe%solver%soln          , ierr); CHKERRQ(ierr)
    call VecZeroEntries(base_soe%solver%soln_prev     , ierr); CHKERRQ(ierr)
    call VecZeroEntries(base_soe%solver%soln_prev_clm , ierr); CHKERRQ(ierr)
    call VecZeroEntries(base_soe%solver%res          , ierr); CHKERRQ(ierr)

    ! Create SNES
    call SNESCreate             (PETSC_COMM_SELF , base_soe%solver%snes, ierr); CHKERRQ(ierr)
    !call SNESSetOptionsPrefix   (base_soe%solver%snes  , "temperature_", ierr); CHKERRQ(ierr)

    call SNESSetTolerances(base_soe%solver%snes, atol, rtol, stol, &
                           max_it, max_f, ierr); CHKERRQ(ierr)

    call SNESSetFunction(base_soe%solver%snes, base_soe%solver%res, SOEResidual, &
         this%soe_ptr, ierr); CHKERRQ(ierr)

    call SNESSetJacobian(base_soe%solver%snes, base_soe%solver%jac, base_soe%solver%jac,     &
         SOEJacobian, this%soe_ptr, ierr); CHKERRQ(ierr)

    call SNESSetFromOptions(base_soe%solver%snes, ierr); CHKERRQ(ierr)

    ! Cleanup
    do igoveq = 1, base_soe%ngoveqns
       call DMDestroy(dms(igoveq), ierr); CHKERRQ(ierr)
    enddo
    deallocate(dms)

  end subroutine MPPSetupProblemSNES

  !------------------------------------------------------------------------
  subroutine MPPSetupProblemKSP(this)
    !
    ! !DESCRIPTION:
    ! Sets the PETSc KSP problem
    !
    use GoverningEquationBaseType , only : goveqn_base_type
    use SystemOfEquationsBasePointerType, only : SOEComputeRHS, SOEComputeOperators
    !
    implicit none
    !
    ! !ARGUMENTS
    class(multiphysicsprob_base_type)      :: this
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type)       , pointer :: cur_goveq
    class(sysofeqns_base_type)    , pointer :: base_soe
    PetscInt                                :: size
    PetscInt                                :: igoveq
    PetscErrorCode                          :: ierr
    DM                            , pointer :: dms(:)
    character(len=256)                      :: name

    base_soe         => this%soe
    this%soe_ptr%ptr => this%soe

    ! Create PETSc DM for each governing equation

    allocate(dms(base_soe%ngoveqns))

    igoveq = 0
    cur_goveq => base_soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       igoveq = igoveq + 1
       size   = cur_goveq%mesh%ncells_local

       call DMDACreate1d(PETSC_COMM_SELF, &
            DM_BOUNDARY_NONE, size, 1, 0, &
            PETSC_NULL_INTEGER, dms(igoveq), ierr);
       CHKERRQ(ierr)

       write(name,*) igoveq
       name = 'fgoveq_' // trim(adjustl(name))
       call DMSetOptionsPrefix(dms(igoveq), name, ierr); CHKERRQ(ierr)

       call DMSetFromOptions(dms(igoveq), ierr); CHKERRQ(ierr)
       call DMSetUp         (dms(igoveq), ierr); CHKERRQ(ierr)

       write(name,*) igoveq
       name = 'goveq_' // trim(adjustl(name))
       call DMDASetFieldName(dms(igoveq), 0, name, ierr); CHKERRQ(ierr)

       cur_goveq => cur_goveq%next
    enddo

    ! DM-Composite approach

    ! Create DMComposite: temperature
    call DMCompositeCreate(PETSC_COMM_SELF, base_soe%solver%dm, ierr); CHKERRQ(ierr)
    call DMSetOptionsPrefix(base_soe%solver%dm, "temperature_", ierr); CHKERRQ(ierr)

    ! Add DMs to DMComposite
    do igoveq = 1, base_soe%ngoveqns
       call DMCompositeAddDM(base_soe%solver%dm, dms(igoveq), ierr); CHKERRQ(ierr)
    enddo

    ! Setup DM
    call DMSetUp(base_soe%solver%dm, ierr); CHKERRQ(ierr)

    ! Create matrix
    call DMCreateMatrix    (base_soe%solver%dm   , base_soe%solver%Amat, ierr); CHKERRQ(ierr)

    call MatSetOption      (base_soe%solver%Amat , MAT_NEW_NONZERO_LOCATION_ERR , &
         PETSC_FALSE, ierr); CHKERRQ(ierr)
    call MatSetOption      (base_soe%solver%Amat , MAT_NEW_NONZERO_ALLOCATION_ERR, &
         PETSC_FALSE, ierr); CHKERRQ(ierr)

    call MatSetFromOptions (base_soe%solver%Amat , ierr); CHKERRQ(ierr)

    ! Create vectors
    call DMCreateGlobalVector(base_soe%solver%dm, base_soe%solver%soln         , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(base_soe%solver%dm, base_soe%solver%rhs          , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(base_soe%solver%dm, base_soe%solver%soln_prev    , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(base_soe%solver%dm, base_soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

    ! Initialize vectors
    call VecZeroEntries(base_soe%solver%soln          , ierr); CHKERRQ(ierr)
    call VecZeroEntries(base_soe%solver%rhs           ,  ierr); CHKERRQ(ierr)
    call VecZeroEntries(base_soe%solver%soln_prev     ,  ierr); CHKERRQ(ierr)
    call VecZeroEntries(base_soe%solver%soln_prev_clm ,  ierr); CHKERRQ(ierr)

    ! Create KSP
    call KSPCreate              (PETSC_COMM_SELF , base_soe%solver%ksp, ierr); CHKERRQ(ierr)
    call KSPSetOptionsPrefix    (base_soe%solver%ksp   , "temperature_", ierr); CHKERRQ(ierr)

    call KSPSetComputeRHS       (base_soe%solver%ksp   , SOEComputeRHS      , &
         this%soe_ptr, ierr); CHKERRQ(ierr)
    call KSPSetComputeOperators (base_soe%solver%ksp   , SOEComputeOperators, &
         this%soe_ptr, ierr); CHKERRQ(ierr)

    call KSPSetFromOptions      (base_soe%solver%ksp   , ierr); CHKERRQ(ierr)

    ! Cleanup
    do igoveq = 1, base_soe%ngoveqns
       call DMDestroy(dms(igoveq), ierr); CHKERRQ(ierr)
    enddo
    deallocate(dms)

  end subroutine MPPSetupProblemKSP

#endif

end module MultiPhysicsProbBaseType
