
module SystemOfEquationsVSFMType

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Object for the Variable Saturate Flow Model (VSFM) system-of-equations
  !
  ! !USES:
  use mpp_varctl                     , only : iulog
  use mpp_abortutils                 , only : endrun
  use mpp_shr_log_mod                , only : errMsg => shr_log_errMsg
  use SystemOfEquationsVSFMAuxType   , only : sysofeqns_vsfm_auxvar_type
  use SystemOfEquationsBaseType      , only : sysofeqns_base_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscts.h"
#include "finclude/petscts.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscdmda.h"
#include "finclude/petscdmda.h90"
#include "finclude/petscviewer.h"

  type, public, extends(sysofeqns_base_type) :: sysofeqns_vsfm_type
     type (sysofeqns_vsfm_auxvar_type), pointer :: aux_vars_in(:)            ! Internal state.
     type (sysofeqns_vsfm_auxvar_type), pointer :: aux_vars_bc(:)            ! Boundary conditions.
     type (sysofeqns_vsfm_auxvar_type), pointer :: aux_vars_ss(:)            ! Source-sink.
     PetscInt, pointer                          :: soe_auxvars_bc_offset (:) ! Cummulative sum of number of control volumes associated with each boundary condition.
     PetscInt, pointer                          :: soe_auxvars_ss_offset (:) ! Cummulative sum of number of control volumes associated with each source-sink condition.
     PetscInt, pointer                          :: soe_auxvars_bc_ncells (:) ! Number of control volumes associated with each boundary condition.
     PetscInt, pointer                          :: soe_auxvars_ss_ncells (:) ! Number of control volumes associated with each source-sink condition.
     PetscInt                                   :: num_auxvars_in            ! Number of auxvars associated with internal state.
     PetscInt                                   :: num_auxvars_in_local      ! Number of auxvars associated with internal state.
     PetscInt                                   :: num_auxvars_bc            ! Number of auxvars associated with boundary condition.
     PetscInt                                   :: num_auxvars_ss            ! Number of auxvars associated with source-sink condition.
   contains
     procedure, public :: Init                   => VSFMSOEInit
     procedure, public :: Setup                  => VSFMSOESetup
     procedure, public :: Residual               => VSFMSOEResidual
     procedure, public :: Jacobian               => VSFMJacobian
     procedure, public :: PreSolve               => VSFMSOEPreSolve
     procedure, public :: SetDataFromCLM         => VSFMSOESetDataFromCLM
     procedure, public :: GetDataForCLM          => VSFMSOEGetDataForCLM
     procedure, public :: PostSolve              => VSFMSOEPostSolve
     procedure, public :: PostStepDT             => VSFMSPostStepDT
     procedure, public :: PreStepDT              => VSFMSPreStepDT
     procedure, public :: GetConditionNames      => VSFMSGetConditionNames
     procedure, public :: SetDataFromCLMForGhost => VSFMSOESetDataFromCLMForGhost
     procedure, public :: ComputeLateralFlux     => VSFMComputeLateralFlux
     procedure, public :: AddGovEqn              => VSFMAddGovEqn
     procedure, public :: CreateVectorsForGovEqn => VSFMCreateVectorsForGovEqn
  end type sysofeqns_vsfm_type

  public :: VSFMSOESetAuxVars,         &
            VSFMSOEUpdateConnections,  &
            VSFMSOEUpdateBCConnections

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine VSFMSOEInit(this)
    !
    ! !DESCRIPTION:
    ! Initializes module variables and data structures
    !
    ! !USES:
    use SystemOfEquationsBaseType, only : SOEBaseInit
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type) :: this

    call SOEBaseInit(this)

    this%num_auxvars_in         = 0
    this%num_auxvars_in_local   = 0
    this%num_auxvars_bc         = 0
    this%num_auxvars_ss         = 0

    nullify(this%aux_vars_in           )
    nullify(this%aux_vars_bc           )
    nullify(this%aux_vars_ss           )

    nullify(this%soe_auxvars_bc_offset )
    nullify(this%soe_auxvars_ss_offset )
    nullify(this%soe_auxvars_bc_ncells )
    nullify(this%soe_auxvars_ss_ncells )

  end subroutine VSFMSOEInit

  !------------------------------------------------------------------------
  subroutine VSFMSOESetup(this, mpi_rank, discretization_type, mpp_type, soe_type, meshes, nmesh)
    !
    ! !DESCRIPTION:
    ! Sets up SoE for the VSFM
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : SOE_RE_ODE
    use MultiPhysicsProbConstants, only : MPP_VSFM_SNES_CLM
    use MeshType, only                  : mesh_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type) :: this
    PetscInt                   :: mpi_rank
    integer, intent(in)                               :: discretization_type
    PetscInt                   :: mpp_type
    PetscInt                   :: soe_type
    class(mesh_type), pointer  :: meshes(:)
    PetscInt                   :: nmesh

    call this%Init()

    this%mpi_rank = mpi_rank

    select case (soe_type)
    case (SOE_RE_ODE)
       this%itype = SOE_RE_ODE
       select case(mpp_type)
       case (MPP_VSFM_SNES_CLM)
          call VSFMSOERichardsEqnODESetup(this, discretization_type, meshes, nmesh)
       case default
          write(iulog,*) 'VSFMSOESetup: Unknown mpp_type'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
    case default
       write(iulog,*) 'VSFMSOESetup: Unknown soe_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine VSFMSOESetup

  !------------------------------------------------------------------------
  subroutine VSFMSOERichardsEqnODESetup(vsfm_soe, discretization_type,  meshes, nmesh)
    !
    ! !DESCRIPTION:
    ! Sets up SoE for the VSFM that uses PETSc SNES.
    !
    ! !USES:
    use MeshType                        , only : mesh_type
    use MeshType                        , only : MeshCreateConnectionSet
    use ConditionType                   , only : condition_type
    use ConditionType                   , only : ConditionNew
    use ConditionType                   , only : ConditionListAddCondition
    use MultiPhysicsProbConstants       , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants       , only : SOIL_CELLS
    use MultiPhysicsProbConstants       , only : COND_MASS_RATE
    use MultiPhysicsProbConstants       , only : COND_MASS_FLUX
    use MultiPhysicsProbConstants       , only : COND_BC
    use MultiPhysicsProbConstants       , only : COND_SS
    use MultiPhysicsProbConstants       , only : COND_NULL
    use MultiPhysicsProbConstants       , only : COND_SEEPAGE_BC
    use MultiPhysicsProbConstants       , only : PRESSURE_REF
    use MultiPhysicsProbConstants       , only : DISCRETIZATION_THREE_DIM
    use MultiPhysicsProbConstants       , only : DISCRETIZATION_VERTICAL_WITH_SS
    use GoveqnRichardsODEPressureType   , only : goveqn_richards_ode_pressure_type
    use EOSWaterMod                     , only : DENSITY_TGDPB01
    use SystemOfEquationsBaseType       , only : SOESetMeshesOfGoveqns
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)                        :: vsfm_soe
    integer, intent(in)                               :: discretization_type
    class(mesh_type), pointer                         :: meshes(:)
    type(condition_type),pointer                      :: bc
    type(condition_type),pointer                      :: ss
    PetscInt                                          :: nmesh
    !
    ! !LOCAL VARIABLES:
    class (goveqn_richards_ode_pressure_type),pointer :: goveq_richards_pres
    PetscInt                                          :: nvars
    PetscInt                                          :: iauxvar
    PetscInt                                          :: icond
    PetscInt                                          :: offset
    PetscInt                                          :: num_bc
    PetscInt                                          :: num_ss
    PetscInt                                          :: total_ncells_for_bc
    PetscInt                                          :: total_ncells_for_ss
    PetscInt, pointer                                 :: ncells_for_bc(:)
    PetscInt, pointer                                 :: ncells_for_ss(:)
    PetscErrorCode                                    :: ierr

    vsfm_soe%name     = "SOEs for VSFM using ODE approach"
    vsfm_soe%ngoveqns = 1

    ! Add governing-equations to the system-of-equations
    allocate(goveq_richards_pres)
    call goveq_richards_pres%Setup()
    goveq_richards_pres%id_in_list = 1

    vsfm_soe%goveqns      => goveq_richards_pres

    ! Assign mesh to each governing equation
    call SOESetMeshesOfGoveqns(vsfm_soe, meshes, nmesh)

    ! Set BCs/SSs for each governing equation
    ss               => ConditionNew()
    ss%name          = 'Infiltration_Flux'
    ss%units         = 'kg/s'
    ss%itype         = COND_MASS_RATE
    ss%region_itype  = SOIL_TOP_CELLS
    allocate(ss%conn_set)
    call MeshCreateConnectionSet(  &
         goveq_richards_pres%mesh, &
         ss%region_itype,          &
         ss%conn_set,              &
         ss%ncells)
    allocate(ss%value(ss%ncells))
    ss%value(:) = 0.d0
    call ConditionListAddCondition(goveq_richards_pres%source_sinks, ss)
    nullify(ss)

    ss              => ConditionNew()
    ss%name         = 'Evapotranspiration_Flux'
    ss%units        = 'kg/s'
    ss%itype        = COND_MASS_RATE
    ss%region_itype = SOIL_CELLS
    allocate(ss%conn_set)
    call MeshCreateConnectionSet(  &
         goveq_richards_pres%mesh, &
         ss%region_itype,          &
         ss%conn_set,              &
         ss%ncells)
    allocate(ss%value(ss%ncells))
    ss%value(:) = 0.d0
    call ConditionListAddCondition(goveq_richards_pres%source_sinks, ss)
    nullify(ss)

    ss               => ConditionNew()
    ss%name          = 'Dew_Flux'
    ss%units         = 'kg/s'
    ss%itype         = COND_MASS_RATE
    ss%region_itype  = SOIL_TOP_CELLS
    allocate(ss%conn_set)
    call MeshCreateConnectionSet(  &
         goveq_richards_pres%mesh, &
         ss%region_itype,          &
         ss%conn_set,              &
         ss%ncells)
    allocate(ss%value(ss%ncells))
    ss%value(:) = 0.d0
    call ConditionListAddCondition(goveq_richards_pres%source_sinks, ss)
    nullify(ss)

    ss               => ConditionNew()
    ss%name          = 'Drainage_Flux'
    ss%units         = 'kg/s'
    ss%itype         = COND_MASS_RATE
    ss%region_itype  = SOIL_CELLS
    allocate(ss%conn_set)
    call MeshCreateConnectionSet(  &
         goveq_richards_pres%mesh, &
         ss%region_itype,          &
         ss%conn_set,              &
         ss%ncells)
    allocate(ss%value(ss%ncells))
    ss%value(:) = 0.d0
    call ConditionListAddCondition(goveq_richards_pres%source_sinks, ss)
    nullify(ss)

    ss               => ConditionNew()
    ss%name          = 'Snow_Disappearance_Flux'
    ss%units         = 'kg/s'
    ss%itype         = COND_MASS_RATE
    ss%region_itype  = SOIL_TOP_CELLS
    allocate(ss%conn_set)
    call MeshCreateConnectionSet(  &
         goveq_richards_pres%mesh, &
         ss%region_itype,          &
         ss%conn_set,              &
         ss%ncells)
    allocate(ss%value(ss%ncells))
    ss%value(:) = 0.d0
    call ConditionListAddCondition(goveq_richards_pres%source_sinks, ss)
    nullify(ss)

    ss               => ConditionNew()
    ss%name          = 'Sublimation_Flux'
    ss%units         = 'kg/s'
    ss%itype         = COND_MASS_RATE
    ss%region_itype  = SOIL_TOP_CELLS
    allocate(ss%conn_set)
    call MeshCreateConnectionSet(  &
         goveq_richards_pres%mesh, &
         ss%region_itype,          &
         ss%conn_set,              &
         ss%ncells)
    allocate(ss%value(ss%ncells))
    ss%value(:) = 0.d0
    call ConditionListAddCondition(goveq_richards_pres%source_sinks, ss)
    nullify(ss)

    if (discretization_type == DISCRETIZATION_THREE_DIM        .or. &
        discretization_type == DISCRETIZATION_VERTICAL_WITH_SS ) then
       ss               => ConditionNew()
       ss%name          = 'Lateral_flux'
       ss%units         = 'kg/s'
       ss%itype         = COND_MASS_RATE
       ss%region_itype  = SOIL_CELLS
       allocate(ss%conn_set)
       call MeshCreateConnectionSet(  &
            goveq_richards_pres%mesh, &
            ss%region_itype,          &
            ss%conn_set,              &
            ss%ncells)
       allocate(ss%value(ss%ncells))
       ss%value(:) = 0.d0
       call ConditionListAddCondition(goveq_richards_pres%source_sinks, ss)
       nullify(ss)

       bc              => ConditionNew()
       bc%name         = 'Seepage_bc'
       bc%units        = 'Pa'
       bc%itype        = COND_SEEPAGE_BC
       bc%region_itype = SOIL_TOP_CELLS
       allocate(bc%conn_set)
       call MeshCreateConnectionSet(  &
            goveq_richards_pres%mesh, &
            bc%region_itype,          &
            bc%conn_set,              &
            bc%ncells)
       allocate(bc%value(bc%ncells))
       bc%value(:) = PRESSURE_REF
       call ConditionListAddCondition(goveq_richards_pres%boundary_conditions, bc)
       nullify(bc)

    endif

    ! Allocate memory for aux vars

    call goveq_richards_pres%AllocateAuxVars()
    !call goveq_richards_pres%SetDensityType(DENSITY_TGDPB01)

    vsfm_soe%num_auxvars_in       = goveq_richards_pres%mesh%ncells_all
    vsfm_soe%num_auxvars_in_local = goveq_richards_pres%mesh%ncells_local

    allocate(vsfm_soe%aux_vars_in(vsfm_soe%num_auxvars_in))

    do iauxvar = 1,vsfm_soe%num_auxvars_in
       call vsfm_soe%aux_vars_in(iauxvar)%Init()
       vsfm_soe%aux_vars_in(iauxvar)%is_in      = PETSC_TRUE
       vsfm_soe%aux_vars_in(iauxvar)%goveqn_id  = 1
    enddo

    call goveq_richards_pres%NumCellsInConditions( &
         COND_BC,                                  & ! Boundary condition ID
         COND_NULL,                                & ! ID of any condition to be excluded
         num_bc,                                   & ! Num. of BCs
         ncells_for_bc                             & ! Num. of control volumes in each BC
         )

    call goveq_richards_pres%NumCellsInConditions( &
         COND_SS,                                  &  ! Source-sink condition ID
         COND_NULL,                                &  ! ID of any condition to be excluded
         num_ss,                                   &  ! Num. of source-sink conditions
         ncells_for_ss                             &  ! Num. of control volumes in each source-sink condition
         )

    allocate(vsfm_soe%soe_auxvars_bc_offset(num_bc))
    allocate(vsfm_soe%soe_auxvars_ss_offset(num_ss))
    allocate(vsfm_soe%soe_auxvars_bc_ncells(num_bc))
    allocate(vsfm_soe%soe_auxvars_ss_ncells(num_ss))

    ! Find total number of control volumes associated with all BCs
    total_ncells_for_bc = 0
    do icond = 1, num_bc
       vsfm_soe%soe_auxvars_bc_offset(icond) = total_ncells_for_bc
       vsfm_soe%soe_auxvars_bc_ncells(icond) = ncells_for_bc(icond)

       total_ncells_for_bc = total_ncells_for_bc + ncells_for_bc(icond)
    enddo

    ! Allocate memory for soe-auxvars associated with BC
    vsfm_soe%num_auxvars_bc = total_ncells_for_bc
    allocate(vsfm_soe%aux_vars_bc(total_ncells_for_bc))

    ! Initialize memory for soe-auxvars associated with BC
    offset = 0
    do icond = 1, num_bc
       do iauxvar = 1, ncells_for_bc(icond)
          call vsfm_soe%aux_vars_bc(iauxvar+offset)%Init()

          vsfm_soe%aux_vars_bc(iauxvar+offset)%is_bc        = PETSC_TRUE
          vsfm_soe%aux_vars_bc(iauxvar+offset)%goveqn_id    = 1
          vsfm_soe%aux_vars_bc(iauxvar+offset)%condition_id = icond
       enddo
       offset = offset + ncells_for_bc(icond)
    enddo

    ! Find total number of control volumes associated with all SS
    total_ncells_for_ss = 0
    do icond = 1, num_ss
       vsfm_soe%soe_auxvars_ss_offset(icond) = total_ncells_for_ss
       vsfm_soe%soe_auxvars_ss_ncells(icond) = ncells_for_ss(icond)

       total_ncells_for_ss = total_ncells_for_ss + ncells_for_ss(icond)
    enddo

    ! Allocate memory for soe-auxvars associated with BC
    vsfm_soe%num_auxvars_ss = total_ncells_for_ss
    allocate(vsfm_soe%aux_vars_ss(total_ncells_for_ss))

    ! Initialize memory for soe-auxvars associated with SS
    offset = 0
    do icond = 1, num_ss
       do iauxvar = 1, ncells_for_ss(icond)
          call vsfm_soe%aux_vars_ss(iauxvar+offset)%Init()

          vsfm_soe%aux_vars_ss(iauxvar+offset)%is_ss        = PETSC_TRUE
          vsfm_soe%aux_vars_ss(iauxvar+offset)%goveqn_id    = 1
          vsfm_soe%aux_vars_ss(iauxvar+offset)%condition_id = icond
       enddo
       offset = offset + ncells_for_ss(icond)
    enddo

    ! Set the variables needed by a given governing equation from
    ! other governing equation

    nvars = 0
    goveq_richards_pres%nvars_needed_from_other_goveqns = nvars

    ! Create PETSc vector to store accumulation term at the begining
    ! of the timestep

    call VecCreateSeq(PETSC_COMM_SELF, goveq_richards_pres%mesh%ncells_local, &
         goveq_richards_pres%accum_prev, ierr)
    CHKERRQ(ierr)

    if (associated(ncells_for_bc)) deallocate(ncells_for_bc)
    if (associated(ncells_for_ss)) deallocate(ncells_for_ss)

  end subroutine VSFMSOERichardsEqnODESetup

  !------------------------------------------------------------------------
  subroutine VSFMSOEResidual(this, snes, X, F, ierr)
    !
    ! !DESCRIPTION:
    ! Performs residual function evaluation for the VSFM
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use MultiPhysicsProbConstants     , only : AUXVAR_INTERNAL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)      :: this
    SNES                            :: snes
    Vec                             :: X
    Vec                             :: F
    PetscErrorCode                  :: ierr
    !
    ! !LOCAL VARIABLES:
    PetscInt                        :: dm_id
    PetscInt                        :: nDM
    PetscInt                        :: offset
    DM, pointer                     :: dms(:)
    Vec, pointer                    :: X_subvecs(:)
    Vec, pointer                    :: F_subvecs(:)
    class(goveqn_base_type),pointer :: cur_goveq
    class(goveqn_base_type),pointer :: cur_goveq_1
    class(goveqn_base_type),pointer :: cur_goveq_2
    PetscInt                        :: row, col
    PetscViewer                     :: viewer
    character(len=256)              :: string

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(this%dm, nDM, ierr); CHKERRQ(ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(this%dm, dms, ierr); CHKERRQ(ierr)

    ! Allocate vectors for individual GEs
    allocate(X_subvecs(    nDM))
    allocate(F_subvecs(    nDM))

    ! Get vectors (X,F) for individual GEs
    call DMCompositeGetAccessArray(this%dm, X, nDM, PETSC_NULL_INTEGER, X_subvecs, &
         ierr); CHKERRQ(ierr)
    call DMCompositeGetAccessArray(this%dm, F, nDM, PETSC_NULL_INTEGER, F_subvecs, &
         ierr); CHKERRQ(ierr)

    ! 1) {X}  ---> sim_aux()
    call VSFMSOEUpdateAuxVarsODE(this, X)

    ! 2.1) GE ---> GetFromSimAux()
    ! Get pointers to governing-equations
    offset = 0
    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)
          class is (goveqn_richards_ode_pressure_type)
             call cur_goveq%GetFromSOEAuxVarsIntrn(this%aux_vars_in, offset)
       end select

       call cur_goveq%UpdateAuxVarsIntrn()
       offset = offset + cur_goveq%mesh%ncells_local

       cur_goveq => cur_goveq%next
    enddo

    ! 3  ) GE_1 <---> GE_2 exchange AuxVars()
    do row = 1,nDM
       do col = row+1,nDM
          call this%SetPointerToIthGovEqn(row, cur_goveq_1)
          call this%SetPointerToIthGovEqn(col, cur_goveq_2)
          call VSFMSOEGovEqnExchangeAuxVars(cur_goveq_1, cur_goveq_2)
          call VSFMSOEGovEqnExchangeAuxVars(cur_goveq_2, cur_goveq_1)
       enddo
    enddo

    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       call cur_goveq%UpdateAuxVars()
       cur_goveq => cur_goveq%next
    enddo

    ! Call Residual
    dm_id = 0
    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       dm_id = dm_id + 1

       call VecZeroEntries(F_subvecs(dm_id), ierr); CHKERRQ(ierr)

       call cur_goveq%Residual(X_subvecs(dm_id), &
            F_subvecs(dm_id),                    &
            ierr); CHKERRQ(ierr)

       cur_goveq => cur_goveq%next
    enddo

    ! Restore vectors (u,udot,F) for individual GEs
    call DMCompositeRestoreAccessArray(this%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)
    call DMCompositeRestoreAccessArray(this%dm, F, nDM, PETSC_NULL_INTEGER, &
         F_subvecs, ierr); CHKERRQ(ierr)

    ! Free memory
    deallocate(dms)
    deallocate(X_subvecs)
    deallocate(F_subvecs)

  end subroutine VSFMSOEResidual

  !------------------------------------------------------------------------
  subroutine VSFMJacobian(this, snes, X, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Computes jacobian for the VSFM
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use MultiPhysicsProbConstants     , only : AUXVAR_INTERNAL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)      :: this
    SNES                            :: snes
    Vec                             :: X
    Mat                             :: A
    Mat                             :: B
    PetscErrorCode                  :: ierr
    !
    ! !LOCAL VARIABLES:
    PetscInt                        :: row
    PetscInt                        :: col
    PetscInt                        :: nDM
    PetscInt                        :: icell

    IS,pointer                      :: is(:)
    DM, pointer                     :: dms(:)
    Vec, pointer                    :: X_subvecs(:)
    Mat, pointer                    :: B_submats(:,:)
    class(goveqn_base_type),pointer :: cur_goveq_1
    class(goveqn_base_type),pointer :: cur_goveq_2
    PetscViewer                     :: viewer
    character(len=256)              :: string

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(this%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(this%dm, dms, ierr); CHKERRQ(ierr)

    ! Allocate vectors for individual GEs
    allocate(X_subvecs(    nDM))

    ! Get vectors (X) for individual GEs
    call DMCompositeGetAccessArray(this%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)

    ! Initialize the matrix
    call MatZeroEntries(B, ierr); CHKERRQ(ierr)

    ! Get submatrices
    allocate(is(nDM))
    allocate(B_submats(nDM,nDM))
    call DMCompositeGetLocalISs(this%dm, is, ierr); CHKERRQ(ierr)
    do row = 1,nDM
       do col = 1,nDM
          call MatGetLocalSubMatrix(B, is(row), is(col), B_submats(row,col), &
               ierr); CHKERRQ(ierr)
       enddo
    enddo

    ! Jacobian and JacobianOffDiag
    row = 0
    cur_goveq_1 => this%goveqns
    do
       if (.not.associated(cur_goveq_1)) exit

       row = row + 1

       call cur_goveq_1%Jacobian(X_subvecs(row), &
            B_submats(row,row),                  &
            B_submats(row,row),                  &
            ierr); CHKERRQ(ierr)

       cur_goveq_2 => cur_goveq_1%next
       col = row
       do
          if (.not.associated(cur_goveq_2)) exit

          col = col + 1

          ! J = dF_1/dx_2
          call cur_goveq_1%JacobianOffDiag(X_subvecs(row),        &
                                           X_subvecs(col),        &
                                           B_submats(row,col),    &
                                           B_submats(row,col),    &
                                           cur_goveq_2%id,        &
                                           cur_goveq_2%id_in_list,&
                                           ierr); CHKERRQ(ierr)

          ! J = dF_2/dx_1
          call cur_goveq_2%JacobianOffDiag(X_subvecs(col),        &
                                           X_subvecs(row),        &
                                           B_submats(col,row),    &
                                           B_submats(col,row),    &
                                           cur_goveq_1%id,        &
                                           cur_goveq_1%id_in_list,&
                                           ierr); CHKERRQ(ierr)

          cur_goveq_2 => cur_goveq_2%next
       enddo

       cur_goveq_1 => cur_goveq_1%next
    enddo

    ! Restore vectors (X) for individual GEs
    call DMCompositeRestoreAccessArray(this%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)

    ! Restore submatrices
    do row = 1,nDM
       do col = 1,nDM
          call MatRestoreLocalSubMatrix(B, is(row), is(col), B_submats(row,col), &
               ierr); CHKERRQ(ierr)
       enddo
    enddo

    ! Destroy IS
    do row = 1,nDM
       call ISDestroy(is(row), ierr); CHKERRQ(ierr)
    enddo

    ! Assemble matrix
    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    call MatAssemblyEnd(  B, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    if ( A /= B) then
       call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
       call MatAssemblyEnd(  A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    endif

    ! Free memory
    deallocate(dms       )
    deallocate(X_subvecs )
    deallocate(is        )
    deallocate(B_submats )

  end subroutine VSFMJacobian

  !------------------------------------------------------------------------
  subroutine VSFMSOEUpdateAuxVarsODE(vsfm_soe, X)
    !
    ! !DESCRIPTION:
    ! Updates the SoE vars for the discretized ODE based on the input
    ! vector X
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : SOE_RE_ODE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)    :: vsfm_soe
    Vec                           :: X
    !
    ! !LOCAL VARIABLES:

    select case (vsfm_soe%itype)
    case (SOE_RE_ODE)
       call VSFMSOEUpdateAuxVarsRichEqnODE(vsfm_soe, X)
    case default
       write(iulog,*) 'VSFMSOEUpdateAuxVars: unknown vsfm_soe%itype'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine VSFMSOEUpdateAuxVarsODE

  !------------------------------------------------------------------------
  subroutine VSFMSOEUpdateAuxVarsRichEqnODE(vsfm_soe, X)
    !
    ! !DESCRIPTION:
    ! Updates the SoE vars for the discretized ODE of Richards equation
    ! based on the input vector X
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type) :: vsfm_soe
    Vec                        :: X
    !
    ! !LOCAL VARIABLES:
    PetscInt                   :: dm_id
    PetscInt                   :: nDM
    DM, pointer                :: dms(:)
    Vec, pointer               :: X_subvecs(:)
    PetscInt                   :: size
    PetscInt                   :: offset
    PetscErrorCode             :: ierr

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(vsfm_soe%dm, nDM, ierr); CHKERRQ(ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(vsfm_soe%dm, dms, ierr); CHKERRQ(ierr)

    ! Allocate vectors for individual GEs
    allocate(X_subvecs(nDM))

    ! Get vectors (X) for individual GEs
    call DMCompositeGetAccessArray(vsfm_soe%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)

    ! Update the SoE auxvars
    offset = 0
    do dm_id = 1, nDM
       call VSFMSOESetAuxVars(vsfm_soe, AUXVAR_INTERNAL, VAR_PRESSURE, &
            X_subvecs(dm_id), offset)
       call VecGetSize(X_subvecs(dm_id), size, ierr); CHKERRQ(ierr)
       offset = offset + size
    enddo

    ! Restore vectors (u,udot,F) for individual GEs
    call DMCompositeRestoreAccessArray(vsfm_soe%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)

    ! Free memory
    deallocate(dms)
    deallocate(X_subvecs)

  end subroutine VSFMSOEUpdateAuxVarsRichEqnODE

  !------------------------------------------------------------------------
  subroutine VSFMSOESetAuxVars(vsfm_soe, auxvar_type, var_type, &
       var_vec, offset)
    !
    ! !DESCRIPTION:
    ! Set values in SoE auxvars.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)                 :: vsfm_soe
    PetscInt                                   :: auxvar_type
    PetscInt, intent(in)                       :: var_type
    Vec                                        :: var_vec
    !
    ! !LOCAL VARIABLES:
    PetscReal, pointer                         :: var_p(:)
    type (sysofeqns_vsfm_auxvar_type), pointer :: avars(:)
    PetscInt                                   :: nauxvar
    PetscInt                                   :: nvar
    PetscInt, optional                         :: offset
    PetscInt                                   :: iauxvar
    PetscInt                                   :: iauxvar_off
    PetscErrorCode                             :: ierr

    select case(auxvar_type)
    case (AUXVAR_INTERNAL)
       avars => vsfm_soe%aux_vars_in
    case default
       write(iulog,*) 'VSFMSOESetAuxVars: auxvar_type not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    if (present(offset)) then
       iauxvar_off = offset
    else
       iauxvar_off = 0
    endif

    nauxvar = size(avars)

    call VecGetLocalSize(var_vec, nvar, ierr); CHKERRQ(ierr)

    if (nvar+iauxvar_off > nauxvar) then
       write(iulog,*) 'VSFMSOESetAuxVars: nvar+iauxvar_off > nauxvar.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    call VecGetArrayF90(var_vec, var_p, ierr); CHKERRQ(ierr)
    do iauxvar = 1, nvar
       call avars(iauxvar + iauxvar_off)%SetValue(var_type, var_p(iauxvar))
    enddo

    call VecRestoreArrayF90(var_vec, var_p, ierr); CHKERRQ(ierr)

  end subroutine VSFMSOESetAuxVars

  !------------------------------------------------------------------------
  subroutine VSFMSOEBasePrintInfo(this)
    !
    ! !DESCRIPTION:
    ! Display information about the VSFM problem
    !
    ! !USES:
    use GoverningEquationBaseType       , only : goveqn_base_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)      :: this
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type),pointer :: cur_goveqn

    write(iulog,*)'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(iulog,*)'  SystemOfEqns_name  : ',trim(this%name)
    write(iulog,*)'  SystemOfEqns_itype : ',this%itype
    write(iulog,*)''

    cur_goveqn => this%goveqns
    do
       if (.not.associated(cur_goveqn)) exit
       call cur_goveqn%PrintInfo()
       cur_goveqn => cur_goveqn%next
    enddo

    write(iulog,*)'  No. of aux_vars_in : ',size(this%aux_vars_in)
    write(iulog,*)''
    write(iulog,*)'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  end subroutine VSFMSOEBasePrintInfo

  !------------------------------------------------------------------------
  subroutine VSFMSOEPreSolve(this)
    !
    ! !DESCRIPTION:
    ! Peform operations before PETSc solver are called.
    ! - For ODE based Richards equation VSFM solver, this subroutine calls
    !   indiviaul governing equations to compute the accumulation term
    !   that depends on the solution from previous time step.
    !
    ! !USES:
    use MultiPhysicsProbConstants     , only : SOE_RE_ODE
    use MultiPhysicsProbConstants     , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants     , only : AUXVAR_BC
    use MultiPhysicsProbConstants     , only : AUXVAR_SS
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type) :: this
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type),pointer :: cur_goveq
    PetscInt :: offset

    select case (this%itype)
    case(SOE_RE_ODE)

       ! 1) {soln_prev}  ---> sim_aux()
       call VSFMSOEUpdateAuxVarsODE(this, this%soln_prev)

       ! 2) GE ---> GetFromSimAux()
       offset = 0
       cur_goveq => this%goveqns
       do
          if (.not.associated(cur_goveq)) exit
          select type(cur_goveq)
             class is (goveqn_richards_ode_pressure_type)

             call cur_goveq%GetFromSOEAuxVarsIntrn(this%aux_vars_in, offset)
             call cur_goveq%GetFromSOEAuxVarsBC(this%aux_vars_bc)
             call cur_goveq%GetFromSOEAuxVarsSS(this%aux_vars_ss)

             call cur_goveq%UpdateAuxVars()

             call cur_goveq%PreSolve()

             offset = offset + cur_goveq%mesh%ncells_local

          end select
          cur_goveq => cur_goveq%next
       enddo

    case default
       write(iulog,*) 'VSFMSOESetup: Unknown soe_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine VSFMSOEPreSolve

  !------------------------------------------------------------------------
  subroutine VSFMSOEPostSolve(this)
    !
    ! !DESCRIPTION:
    ! Peform operations after a successful call to the PETSc solver.
    !
    ! !USES:
    use MultiPhysicsProbConstants     , only : SOE_RE_ODE
    use MultiPhysicsProbConstants     , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants     , only : AUXVAR_BC
    use MultiPhysicsProbConstants     , only : AUXVAR_SS
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)      :: this
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type),pointer :: cur_goveq
    PetscErrorCode                  :: ierr

    call VecCopy(this%soln, this%soln_prev,ierr); CHKERRQ(ierr)

    select case (this%itype)
    case(SOE_RE_ODE)

       cur_goveq => this%goveqns
       do
          if (.not.associated(cur_goveq)) exit
          select type(cur_goveq)
             class is (goveqn_richards_ode_pressure_type)

             call cur_goveq%SetDataInSOEAuxVar(AUXVAR_INTERNAL, this%aux_vars_in)
             call cur_goveq%SetDataInSOEAuxVar(AUXVAR_BC      , this%aux_vars_bc)

          end select
          cur_goveq => cur_goveq%next
       enddo

    case default
       write(iulog,*) 'VSFMSOESetup: Unknown soe_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine VSFMSOEPostSolve

  !------------------------------------------------------------------------
  subroutine VSFMSOESetDataFromCLM(this, soe_auxvar_type, var_type, &
       soe_auxvar_id, data_1d)
    !
    ! !DESCRIPTION:
    ! Used by CLM to set values of boundary conditions and source-sink
    ! terms for the VSFM solver.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : SOE_RE_ODE
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants, only : AUXVAR_BC
    use MultiPhysicsProbConstants, only : AUXVAR_SS
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)                 :: this
    PetscInt, intent(in)                       :: var_type
    PetscInt                                   :: soe_auxvar_type
    PetscInt                                   :: soe_auxvar_id
    PetscReal                                  :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                   :: iauxvar
    PetscInt                                   :: iauxvar_off
    PetscInt                                   :: nauxvar
    type (sysofeqns_vsfm_auxvar_type), pointer :: auxvars(:)

    select case(soe_auxvar_type)
    case(AUXVAR_INTERNAL)
       auxvars      => this%aux_vars_in
       iauxvar_off  = 0
       nauxvar      = this%num_auxvars_in
    case(AUXVAR_BC)
       auxvars      => this%aux_vars_bc
       iauxvar_off  = this%soe_auxvars_bc_offset(soe_auxvar_id)
       nauxvar      = this%soe_auxvars_bc_ncells(soe_auxvar_id)
    case(AUXVAR_SS)
       auxvars      => this%aux_vars_ss
       iauxvar_off  = this%soe_auxvars_ss_offset(soe_auxvar_id)
       nauxvar      = this%soe_auxvars_ss_ncells(soe_auxvar_id)
    case default
       write(iulog,*) 'VSFMSOESetDataFromCLM: Unknown soe_auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    if (size(data_1d) > nauxvar) then
       write(iulog,*) 'VSFMSOESetDataFromCLM: size(data_1d) > nauxvar'
       write(iulog,*) 'size(data_1d) = ',size(data_1d)
       write(iulog,*) 'nauxvar       = ', nauxvar
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    do iauxvar = 1, size(data_1d)
       call auxvars(iauxvar + iauxvar_off)%SetValue(var_type, data_1d(iauxvar))
    enddo

  end subroutine VSFMSOESetDataFromCLM

  !------------------------------------------------------------------------
  subroutine VSFMSOESetDataFromCLMForGhost(this, soe_auxvar_type, var_type, &
       soe_auxvar_id, data_1d)
    !
    ! !DESCRIPTION:
    ! Used by CLM to set values of boundary conditions and source-sink
    ! terms for the VSFM solver.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : SOE_RE_ODE
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants, only : AUXVAR_BC
    use MultiPhysicsProbConstants, only : AUXVAR_SS
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)                 :: this
    PetscInt, intent(in)                       :: var_type
    PetscInt                                   :: soe_auxvar_type
    PetscInt                                   :: soe_auxvar_id
    PetscReal                                  :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                   :: iauxvar
    PetscInt                                   :: iauxvar_beg
    PetscInt                                   :: iauxvar_end
    PetscInt                                   :: nauxvar
    type (sysofeqns_vsfm_auxvar_type), pointer :: auxvars(:)

    select case(soe_auxvar_type)
    case(AUXVAR_INTERNAL)
       auxvars      => this%aux_vars_in
       nauxvar      = this%num_auxvars_in
       iauxvar_beg  = this%num_auxvars_in_local + 1
       iauxvar_end  = nauxvar
    case default
       write(iulog,*) 'VSFMSOESetDataFromCLM: Unknown soe_auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    if (size(data_1d) > nauxvar) then
       write(iulog,*) 'VSFMSOESetDataFromCLMForGhost: size(data_1d) > nauxvar'
       write(iulog,*) 'size(data_1d) = ',size(data_1d)
       write(iulog,*) 'nauxvar       = ', nauxvar
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    do iauxvar = iauxvar_beg, iauxvar_end
       call auxvars(iauxvar)%SetValue(var_type, data_1d(iauxvar))
    enddo

  end subroutine VSFMSOESetDataFromCLMForGhost

  !------------------------------------------------------------------------
  subroutine VSFMSOEGetDataForCLM(this, soe_auxvar_type, var_type, &
       soe_auxvar_id, data_1d)
    !
    ! !DESCRIPTION:
    ! Used by CLM to extracted values from the VSFM solver
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : SOE_RE_ODE
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants, only : AUXVAR_BC
    use MultiPhysicsProbConstants, only : AUXVAR_SS
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)                 :: this
    PetscInt, intent(in)                       :: var_type
    PetscInt                                   :: soe_auxvar_type
    PetscInt                                   :: soe_auxvar_id
    PetscReal                                  :: data_1d(:)
    PetscInt                                   :: nsize
    !
    ! !LOCAL VARIABLES:
    PetscInt                                   :: iauxvar
    PetscInt                                   :: iauxvar_off
    PetscInt                                   :: nauxvar
    PetscReal                                  :: var_value
    type (sysofeqns_vsfm_auxvar_type), pointer :: auxvars(:)

    select case(soe_auxvar_type)
    case(AUXVAR_INTERNAL)
       auxvars      => this%aux_vars_in
       iauxvar_off  = 0
       nauxvar      = this%num_auxvars_in
    case(AUXVAR_SS)
       auxvars      => this%aux_vars_ss
       iauxvar_off  = this%soe_auxvars_ss_offset(soe_auxvar_id)
       nauxvar      = this%soe_auxvars_ss_ncells(soe_auxvar_id)
    case(AUXVAR_BC)
       auxvars      => this%aux_vars_bc
       iauxvar_off  = this%soe_auxvars_bc_offset(soe_auxvar_id)
       nauxvar      = this%soe_auxvars_bc_ncells(soe_auxvar_id)
    case default
       write(iulog,*) 'VSFMSOEGetDataForCLM: Unknown soe_auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    if (size(data_1d) > nauxvar) then
       write(iulog,*) 'VSFMSOEGetDataForCLM: size(data_1d) > nauxvar'
       write(iulog,*) 'size(data_1d) = ',size(data_1d)
       write(iulog,*) 'nauxvar       = ', nauxvar
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    do iauxvar = 1, size(data_1d)
       call auxvars(iauxvar + iauxvar_off)%GetValue(var_type, var_value)
       data_1d(iauxvar) = var_value
    enddo

  end subroutine VSFMSOEGetDataForCLM

  !------------------------------------------------------------------------
  subroutine VSFMComputeLateralFlux(this, dt)
    !
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use MultiPhysicsProbConstants     , only : AUXVAR_SS
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)    :: this
    PetscReal                     :: dt
    !
    class(goveqn_base_type),pointer :: cur_goveq

    ! Set timestep
    call this%SetDtime(dt)

    ! Do any pre-solve operations
    call this%PreSolve()

    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)
          class is (goveqn_richards_ode_pressure_type)
             call cur_goveq%ComputeLateralFlux()
             call cur_goveq%SetDataInSOEAuxVar(AUXVAR_SS, this%aux_vars_ss)
       end select

       cur_goveq => cur_goveq%next
    enddo

  end subroutine VSFMComputeLateralFlux

  !------------------------------------------------------------------------
  subroutine VSFMSPreStepDT(this)
    !
    ! !DESCRIPTION:
    ! This subroutines copies solution vector at previous CLM time step
    ! before StepDT is called.
    !
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)      :: this
    !
    class(goveqn_base_type),pointer :: cur_goveq
    PetscErrorCode                  :: ierr
    PetscInt                        :: iauxvar

    call VecCopy(this%soln_prev_clm, this%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(this%soln_prev_clm, this%soln     , ierr); CHKERRQ(ierr)

    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)
          class is (goveqn_richards_ode_pressure_type)
          call cur_goveq%PreStepDT()
       end select
       cur_goveq => cur_goveq%next
    enddo

  end subroutine VSFMSPreStepDT

  !------------------------------------------------------------------------
  subroutine VSFMSPostStepDT(this)
    !
    ! !DESCRIPTION:
    ! This subroutines make a copy of solution vector post StepDT is
    ! called.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type) :: this
    PetscErrorCode             :: ierr

    call VecCopy(this%soln_prev, this%soln_prev_clm, ierr); CHKERRQ(ierr)

  end subroutine VSFMSPostStepDT

  !------------------------------------------------------------------------
  subroutine VSFMSGetConditionNames(this, cond_type, cond_type_to_exclude, &
       num_conds, cond_names)
    !
    ! !DESCRIPTION:
    ! Returns the total number and names of conditions (eg. boundary condition
    ! or source-sink) present
    !
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)      :: this
    PetscInt, intent(in)            :: cond_type
    PetscInt, intent(in)            :: cond_type_to_exclude
    character (len=256), pointer    :: cond_names(:)
    PetscInt                        :: num_conds
    !
    class(goveqn_base_type),pointer :: cur_goveq
    PetscInt                        :: nn
    PetscErrorCode                  :: ierr
    PetscInt                        :: num_conds_tmp
    character (len=256), pointer    :: cond_names_tmp(:)

    num_conds = 0
    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)
          class is (goveqn_richards_ode_pressure_type)
             call cur_goveq%NumConditions(cond_type, cond_type_to_exclude, num_conds_tmp)
          class default
             write(iulog,*) 'VSFMSGetConditionNames: Supported for only goveqn_richards_ode_pressure_type'
             call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
       num_conds = num_conds + num_conds_tmp
       cur_goveq => cur_goveq%next
    enddo

    allocate(cond_names(num_conds))

    num_conds = 0
    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)
          class is (goveqn_richards_ode_pressure_type)
          call cur_goveq%GetConditionNames(cond_type, cond_type_to_exclude, num_conds_tmp, cond_names_tmp)
          if (num_conds_tmp > 0) then
             do nn = 1, num_conds_tmp
                num_conds = num_conds + 1
                cond_names(num_conds) = cond_names_tmp(nn)
             enddo
             deallocate(cond_names_tmp)
          endif
       end select
       cur_goveq => cur_goveq%next
    enddo

  end subroutine VSFMSGetConditionNames

  !------------------------------------------------------------------------
  subroutine VSFMAddGovEqn(this, geq_type, name, mesh_itype)
    !
    ! !DESCRIPTION:
    ! Adds a governing equation to system-of-equations
    !
    ! !USES:
    use SystemOfEquationsBaseType     , only : SOEBaseInit
    use GoverningEquationBaseType     , only : goveqn_base_type
    use MultiPhysicsProbConstants     , only : GE_RE
    use MultiPhysicsProbConstants     , only : MESH_CLM_SOIL_COL
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type) :: this
    PetscInt                   :: geq_type
    character(len=*)           :: name
    PetscInt                   :: mesh_itype
    !
    ! !LOCAL VARIABLES:
    class (goveqn_richards_ode_pressure_type) , pointer :: goveq_richards
    class(goveqn_base_type),pointer                     :: cur_goveqn
    integer                                             :: igoveqn

    cur_goveqn => this%goveqns

    do igoveqn = 1, this%ngoveqns - 1
       cur_goveqn => cur_goveqn%next
    enddo

    this%ngoveqns = this%ngoveqns + 1

    select case(geq_type)
       case (GE_RE)

          allocate(goveq_richards)
          call goveq_richards%Setup()

          goveq_richards%name        = trim(name)
          goveq_richards%id_in_list  = this%ngoveqns
          goveq_richards%mesh_itype  = mesh_itype

          if (this%ngoveqns == 1) then
             this%goveqns => goveq_richards
          else
             cur_goveqn%next => goveq_richards
          endif

       case default
          write(iulog,*) 'Unknown governing equation type'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select

     end subroutine VSFMAddGovEqn

  !------------------------------------------------------------------------
  subroutine VSFMCreateVectorsForGovEqn(this)
    !
    ! !DESCRIPTION:
    ! Creates vectors required by each governing equation
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)      :: this
    !
    class(goveqn_base_type),pointer :: cur_goveq

    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          call cur_goveq%CreateVectors()
       class default
          write(iulog,*) 'Unsupported cur_goveq type'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
       cur_goveq => cur_goveq%next
    enddo

  end subroutine VSFMCreateVectorsForGovEqn

  !------------------------------------------------------------------------

  subroutine VSFMSOEGovEqnExchangeAuxVars(cur_goveq_1, cur_goveq_2)
    !
    ! !DESCRIPTION:
    !
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use ConnectionSetType             , only : connection_set_type
    use ConditionType                 , only : condition_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type),pointer :: cur_goveq_1
    class(goveqn_base_type),pointer :: cur_goveq_2
    !
    ! !LOCAL VARIABLES:
    type(connection_set_type), pointer                  :: cur_conn_set_2
    type(condition_type),pointer                        :: cur_cond_2
    PetscInt                                            :: idx
    PetscInt                                            :: iauxvar
    PetscInt                                            :: ivar
    PetscInt                                            :: var_type
    PetscInt                                            :: ghosted_id
    PetscInt                                            :: bc_idx
    PetscInt                                            :: bc_offset
    PetscInt                                            :: bc_auxvar_idx_of_other_goveqn
    PetscReal                                           :: var_value
    PetscBool                                           :: bc_found
    PetscBool                                           :: bc_type

    do ivar = 1,cur_goveq_1%nvars_needed_from_other_goveqns
       if (cur_goveq_1%ids_of_other_goveqns(ivar) == &
           cur_goveq_2%id_in_list) then

          var_type                      = cur_goveq_1%var_ids_needed_from_other_goveqns(ivar)
          bc_type                       = cur_goveq_1%is_bc_auxvar_type(ivar)
          bc_offset                     = cur_goveq_1%bc_auxvar_offset(ivar)
          bc_auxvar_idx_of_other_goveqn = cur_goveq_1%bc_auxvar_idx_of_other_goveqn(ivar)

          !GB: Presently assuming that auxvars for ALL cells is needed from
          !    the other goveqn. Additionally, it is assumed that the retrieved
          !    data is needed for 'aux_vars_in' and not for boundary-conditions or
          !    source-sinks
          if (.not.bc_type) then
             do ghosted_id = 1,cur_goveq_1%mesh%ncells_local

                select type(cur_goveq_2)
                class is (goveqn_richards_ode_pressure_type)
                   call cur_goveq_2%aux_vars_in(ghosted_id)%GetValue(var_type, var_value)
                class default
                   write(iulog,*)'VSFMSOEGovEqnExchangeAuxVars: Unknown class'
                   call endrun(msg=errMsg(__FILE__, __LINE__))
                end select

                select type(cur_goveq_1)
                class is (goveqn_richards_ode_pressure_type)
                   if (.not.bc_type) then
                      call cur_goveq_1%aux_vars_in(ghosted_id)%SetValue(var_type, var_value)
                   else
                      call cur_goveq_1%aux_vars_bc(ghosted_id+bc_offset)%SetValue(var_type, var_value)
                   endif

                class default
                   write(iulog,*)'VSFMSOEGovEqnExchangeAuxVars: Unknown class'
                   call endrun(msg=errMsg(__FILE__, __LINE__))
                end select

             enddo

           else

              bc_idx = 1
              bc_found = PETSC_FALSE
              cur_cond_2 => cur_goveq_2%boundary_conditions%first
              do
                 if (.not.associated(cur_cond_2)) exit
                 cur_conn_set_2 => cur_cond_2%conn_set
                 if (bc_idx == bc_auxvar_idx_of_other_goveqn) then
                    bc_found = PETSC_TRUE
                    exit
                 endif

                 bc_idx = bc_idx + 1
                 cur_cond_2 => cur_cond_2%next
              enddo

              if (.not.bc_found) then
                 write(iulog,*) 'VSFMSOEGovEqnExchangeAuxVars: BC not found'
                 call endrun(msg=errMsg(__FILE__, __LINE__))
              endif

              if (cur_conn_set_2%num_connections /= &
                  cur_goveq_1%bc_auxvar_ncells(ivar) ) then
                 write(iulog,*) 'VSFMSOEGovEqnExchangeAuxVars: Number of '
                 call endrun(msg=errMsg(__FILE__, __LINE__))
              endif

              do iauxvar = 1, cur_goveq_1%bc_auxvar_ncells(ivar)

                 ! Get value from cur_goveq_2
                 select type(cur_goveq_2)
                 class is (goveqn_richards_ode_pressure_type)
                    idx = cur_conn_set_2%id_dn(iauxvar)
                    call cur_goveq_2%aux_vars_in(idx)%GetValue(var_type, var_value)
                 class default
                    write(iulog,*)'VSFMSOEGovEqnExchangeAuxVars: Unknown class'
                    call endrun(msg=errMsg(__FILE__, __LINE__))
                 end select

                 idx = iauxvar + bc_offset

                 ! Set value from cur_goveq_1
                 select type(cur_goveq_1)
                 class is (goveqn_richards_ode_pressure_type)
                    call cur_goveq_1%aux_vars_bc(idx)%SetValue(var_type, var_value)
                 class default
                    write(iulog,*)'VSFMSOEGovEqnExchangeAuxVars: Unknown class'
                    call endrun(msg=errMsg(__FILE__, __LINE__))
                 end select
              enddo
           endif ! if (.not.bc_type)

        endif ! (
     enddo

  end subroutine VSFMSOEGovEqnExchangeAuxVars

  !------------------------------------------------------------------------
  subroutine VSFMSOEUpdateConnections(this, mpp_id)
    !
    ! !DESCRIPTION:
    !
    use GoverningEquationBaseType     , only : goveqn_base_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)    :: this
    PetscInt, intent(in)          :: mpp_id
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type),pointer :: cur_goveq_1
    class(goveqn_base_type),pointer :: cur_goveq_2
    PetscInt :: row, col
    PetscBool :: check

    check = PETSC_TRUE

    do row = 1, this%ngoveqns
       do col = row+1, this%ngoveqns
          call this%SetPointerToIthGovEqn(row, cur_goveq_1)
          call this%SetPointerToIthGovEqn(col, cur_goveq_2)
          call VSFMSOEUpdateBCConnections(cur_goveq_1, cur_goveq_2, mpp_id)
       enddo
    enddo

  end subroutine VSFMSOEUpdateConnections

  !------------------------------------------------------------------------

  subroutine VSFMSOEUpdateBCConnections(cur_goveq_1, cur_goveq_2, mpp_id)
    !
    ! !DESCRIPTION:
    !
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use ConnectionSetType, only             : connection_set_type
    use ConditionType, only                 : condition_type
    use RichardsODEPressureAuxType
    use SaturationFunction
    use MultiPhysicsProbConstants , only : MPP_VSFM_SNES_CLM
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type),pointer :: cur_goveq_1
    class(goveqn_base_type),pointer :: cur_goveq_2
    PetscInt, intent(in)            :: mpp_id
    !
    ! !LOCAL VARIABLES:
    type(connection_set_type), pointer                  :: cur_conn_set_1
    type(connection_set_type), pointer                  :: cur_conn_set_2
    type(condition_type),pointer                        :: cur_cond_1
    type(condition_type),pointer                        :: cur_cond_2
    type (rich_ode_pres_auxvar_type), pointer           :: aux_vars_bc_1(:)
    type (rich_ode_pres_auxvar_type), pointer           :: aux_vars_bc_2(:)
    type (rich_ode_pres_auxvar_type)                    :: tmp_aux_var_bc_1
    type (rich_ode_pres_auxvar_type)                    :: tmp_aux_var_bc_2
    PetscInt                                            :: iauxvar
    PetscInt                                            :: ivar
    PetscInt                                            :: sum_conn_1
    PetscInt                                            :: sum_conn_2
    PetscInt                                            :: bc_idx
    PetscInt                                            :: bc_auxvar_idx
    PetscInt                                            :: bc_auxvar_idx_of_other_goveqn
    PetscInt                                            :: id_dn_1
    PetscInt                                            :: id_dn_2
    PetscReal                                           :: x_up, y_up, z_up
    PetscReal                                           :: x_dn, y_dn, z_dn
    PetscReal                                           :: x_half, y_half, z_half
    PetscReal                                           :: dx, dy, dz
    PetscReal                                           :: dist, dist_up, dist_dn
    PetscReal                                           :: area_1, area_2
    PetscBool                                           :: bc_found
    PetscBool                                           :: bc_type

    PetscReal, dimension(3) :: vec3Swap
    PetscReal :: valSwap
    type(saturation_params_type) :: satParamsSwap

    select type(cur_goveq_1)
    class is (goveqn_richards_ode_pressure_type)
        aux_vars_bc_1 => cur_goveq_1%aux_vars_bc
    class default
        write(iulog,*)'VSFMSOEUpdateBCConnections: Unknown class'
        call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    select type(cur_goveq_2)
    class is (goveqn_richards_ode_pressure_type)
        aux_vars_bc_2 => cur_goveq_2%aux_vars_bc
    class default
        write(iulog,*)'VSFMSOEUpdateBCConnections: Unknown class'
        call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    do ivar = 1,cur_goveq_1%nvars_needed_from_other_goveqns
       if (cur_goveq_1%ids_of_other_goveqns(ivar) == &
           cur_goveq_2%id_in_list) then

          bc_type                        = cur_goveq_1%is_bc_auxvar_type(ivar)
          bc_auxvar_idx                  = cur_goveq_1%bc_auxvar_idx(ivar)
          bc_auxvar_idx_of_other_goveqn  = cur_goveq_1%bc_auxvar_idx_of_other_goveqn(ivar)

          if (bc_type) then

             bc_idx = 1
             sum_conn_1 = 0
             bc_found = PETSC_FALSE
             cur_cond_1 => cur_goveq_1%boundary_conditions%first
             do
                if (.not.associated(cur_cond_1)) exit
                cur_conn_set_1 => cur_cond_1%conn_set
                if (bc_idx == bc_auxvar_idx) then
                   bc_found = PETSC_TRUE
                   exit
                endif

                bc_idx = bc_idx + 1
                sum_conn_1 = sum_conn_1 + cur_conn_set_1%num_connections
                cur_cond_1 => cur_cond_1%next
             enddo

             if (.not.bc_found) then
                write(iulog,*) 'VSFMSOEGovEqnExchangeAuxVars: BC not found'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             bc_idx = 1
             sum_conn_2 = 0
             bc_found = PETSC_FALSE
             cur_cond_2 => cur_goveq_2%boundary_conditions%first
             do
                if (.not.associated(cur_cond_2)) exit
                cur_conn_set_2 => cur_cond_2%conn_set
                if (bc_idx == bc_auxvar_idx_of_other_goveqn) then
                   bc_found = PETSC_TRUE
                   exit
                endif

                bc_idx = bc_idx + 1
                sum_conn_2 = sum_conn_2 + cur_conn_set_2%num_connections
                cur_cond_2 => cur_cond_2%next
             enddo

             if (.not.bc_found) then
                write(iulog,*) 'VSFMSOEGovEqnExchangeAuxVars: BC not found'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             if (cur_conn_set_1%num_connections /= &
                 cur_conn_set_2%num_connections ) then
                write(iulog,*) 'VSFMSOEGovEqnExchangeAuxVars: Number of cells in BC of eq-1 and eq-2 are not same.'
                write(iulog,*) 'Eq-1'
                write(iulog,*) '   Name      : ',trim(cur_goveq_1%name)
                write(iulog,*) '   BC name   : ',trim(cur_cond_1%name)
                write(iulog,*) '   BC ncells : ',cur_conn_set_1%num_connections
                write(iulog,*) 'Eq-2'
                write(iulog,*) '   Name      : ',trim(cur_goveq_2%name)
                write(iulog,*) '   BC name   : ',trim(cur_cond_2%name)
                write(iulog,*) '   BC ncells : ',cur_conn_set_2%num_connections
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             if (cur_goveq_2%id_in_list > cur_goveq_1%id_in_list) then
                cur_cond_2%swap_order = PETSC_TRUE
             else
                cur_cond_1%swap_order = PETSC_TRUE
             endif

             do iauxvar = 1, cur_conn_set_1%num_connections


               !
               ! Eq-1
               !              unit_vec
               !            <----------
               !        "dn"           "up"
               !    _____________ _____________
               !   |             |             |
               !   |             |             |
               !   |      x      o      x      |
               !   |             |             |
               !   |   mesh-eq1  |   mesh-eq2  |
               !   |_____________|_____________|
               !
               ! Eq-2
               !        "up"           "dn"
               !            ---------->
               !              unit_vec
               !

                id_dn_1 = cur_conn_set_1%id_dn(iauxvar)
                id_dn_2 = cur_conn_set_2%id_dn(iauxvar)

                x_dn = cur_goveq_1%mesh%x(id_dn_1)
                y_dn = cur_goveq_1%mesh%y(id_dn_1)
                z_dn = cur_goveq_1%mesh%z(id_dn_1)

                x_up = cur_goveq_2%mesh%x(id_dn_2)
                y_up = cur_goveq_2%mesh%y(id_dn_2)
                z_up = cur_goveq_2%mesh%z(id_dn_2)

                dx = -(x_up - x_dn)
                dy = -(y_up - y_dn)
                dz = -(z_up - z_dn)

                dist = (dx**2.d0 + dy**2.d0 + dz**2.d0)**0.5d0

                cur_conn_set_1%dist_unitvec(iauxvar)%arr(1) = dx/dist
                cur_conn_set_1%dist_unitvec(iauxvar)%arr(2) = dy/dist
                cur_conn_set_1%dist_unitvec(iauxvar)%arr(3) = dz/dist

                x_half = (x_dn + x_up)/2.d0
                y_half = (y_dn + y_up)/2.d0
                z_half = (z_dn + z_up)/2.d0

                cur_conn_set_1%id_up(iauxvar)   = id_dn_2

                ! unit_vector for eq1 = -unit_vector for eq2
                cur_conn_set_2%dist_unitvec(iauxvar)%arr = &
                    -cur_conn_set_1%dist_unitvec(iauxvar)%arr

                ! dist_up eq2 = dist_dn eq1
                ! dist_up eq1 = dist_dn eq2
                cur_conn_set_2%dist_up(iauxvar) = cur_conn_set_1%dist_dn(iauxvar)
                cur_conn_set_1%dist_up(iauxvar) = cur_conn_set_2%dist_dn(iauxvar)

                cur_conn_set_2%id_up(iauxvar)   = id_dn_1

                sum_conn_1 = sum_conn_1 + 1
                sum_conn_2 = sum_conn_2 + 1

                call tmp_aux_var_bc_1%Copy(aux_vars_bc_1(sum_conn_1))
                call tmp_aux_var_bc_2%Copy(aux_vars_bc_2(sum_conn_2))

                call aux_vars_bc_1(sum_conn_1)%Copy(tmp_aux_var_bc_2)
                call aux_vars_bc_2(sum_conn_2)%Copy(tmp_aux_var_bc_1)

             enddo

          endif
       endif
    enddo

  end subroutine VSFMSOEUpdateBCConnections

#endif

end module SystemOfEquationsVSFMType
