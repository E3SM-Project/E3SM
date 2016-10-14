
module MultiPhysicsProbVSFM

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Object for variably saturated flow model
  !-----------------------------------------------------------------------

  ! !USES:
  use mpp_varctl                         , only : iulog
  use mpp_abortutils                     , only : endrun
  use mpp_shr_log_mod                    , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbBaseType           , only : multiphysicsprob_base_type
  use SystemOfEquationsVSFMType          , only : sysofeqns_vsfm_type
  use SystemOfEquationsBasePointerType   , only : sysofeqns_base_pointer_type

  implicit none
  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscts.h"
#include "finclude/petscts.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscsnes.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscdmda.h"
#include "finclude/petscdmda.h90"
#include "finclude/petscviewer.h"

  type, public, extends(multiphysicsprob_base_type) :: mpp_vsfm_type
     class(sysofeqns_vsfm_type),pointer          :: sysofeqns
     type(sysofeqns_base_pointer_type), pointer  :: sysofeqns_ptr
   contains
     procedure, public :: Init                         => VSFMMPPInit
     procedure, public :: Clean                        => VSFMMPPClean
     procedure, public :: Setup                        => VSFMMPPSetup
     procedure, public :: Restart                      => VSFMMPPRestart
     procedure, public :: UpdateSysOfEqnsAuxVars       => VSFMMPPUpdateSysOfEqnsAuxVars
     procedure, public :: SetMPIRank                   => VSFMMPPSetMPIRank
     procedure, public :: SetMeshesOfGoveqns           => VSFMMPPSetMeshesOfGoveqns
     procedure, public :: AddGovEqn                    => VSFMMPPAddGovEqn
     procedure, public :: GovEqnAddCondition           => VSFMMPPGovEqnAddCondition
     procedure, public :: AllocateAuxVars              => VSFMMPPAllocateAuxVars
     procedure, public :: SetupProblem                 => VSFMMPPSetupProblem
     procedure, public :: GovEqnUpdateConditionConnSet => VSFMMPPGovEqnUpdateConditionConnSet
     procedure, public :: GovEqnSetCouplingVars        => VSFMMPPGovEqnSetCouplingVars

  end type mpp_vsfm_type

  type(mpp_vsfm_type), public, target :: vsfm_mpp
  public :: VSFMMPPSetSoils

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine VSFMMPPInit(this)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbBaseType , only : MPPBaseInit
    use MultiPhysicsProbConstants, only : MPP_VSFM_SNES_CLM
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type) :: this

    call MPPBaseInit(this)

    allocate(this%sysofeqns)
    call this%sysofeqns%Init()

    allocate(this%sysofeqns_ptr)
    nullify(this%sysofeqns_ptr%ptr)

  end subroutine VSFMMPPInit

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetMPIRank(this, rank)
    !
    ! !DESCRIPTION:
    ! Sets MPI rank
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)    :: this
    PetscInt                :: rank

    if (associated(this%sysofeqns)) then
       this%sysofeqns%mpi_rank = rank
    endif

  end subroutine VSFMMPPSetMPIRank

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetup(this, begg, endg, begc, endc, mpi_rank, &
       ugrid, grc_landunit_indices, lun_coli, lun_colf,           &
       discretization_type, ncols_ghost, filter_vsfmc,            &
       xc_col, yc_col, zc_col, z, zi, dz,                         &
       area_col, grid_owner, col_itype,                           &
       watsat, hksat, bsw, sucsat, eff_porosity, residual_sat,    &
       zwt, vsfm_satfunc_type, density_type)
    !
    ! !DESCRIPTION:
    ! Sets up the Variably Saturated Flow Model (VSFM) - Multi-Phyiscs Problem 
    ! (MPP):
    !
    ! - Creates the mesh,
    ! - Sets up the system of equation (SoE) in the VSFM-MPP,
    ! - Sets up the PETSc SNES, and
    ! - Initializes the VSFM-MPP.
    !
    ! !USES
    use MultiPhysicsProbConstants, only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants, only : PETSC_SNES
    use MultiPhysicsProbConstants, only : MPP_VSFM_SNES_CLM
    use MultiPhysicsProbConstants, only : SOE_RE_ODE
    use UnstructuredGridType     , only : ugrid_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)           :: this
    integer, intent(in)            :: begg,endg
    integer, intent(in)            :: begc,endc
    PetscInt, intent(in)           :: mpi_rank
    type(ugrid_type), pointer      :: ugrid
    integer, intent(in)            :: grc_landunit_indices(:,:)
    integer, intent(in)            :: lun_coli(:)
    integer, intent(in)            :: lun_colf(:)
    integer, intent(in)            :: discretization_type
    integer, intent(in)            :: ncols_ghost          ! number of ghost/halo columns
    PetscInt, pointer, intent(in)  :: filter_vsfmc(:) ! column filter for soil points
    PetscReal, pointer, intent(in) :: xc_col(:)
    PetscReal, pointer, intent(in) :: yc_col(:)
    PetscReal, pointer, intent(in) :: zc_col(:)
    PetscReal, pointer, intent(in) :: z(:,:)
    PetscReal, pointer, intent(in) :: zi(:,:)
    PetscReal, pointer, intent(in) :: dz(:,:)
    PetscReal, pointer, intent(in) :: area_col(:)
    PetscInt, pointer, intent(in)  :: grid_owner(:)
    PetscInt, pointer, intent(in)  :: col_itype(:)
    PetscReal, intent(in), pointer :: watsat(:,:)
    PetscReal, intent(in), pointer :: hksat(:,:)
    PetscReal, intent(in), pointer :: bsw(:,:)
    PetscReal, intent(in), pointer :: sucsat(:,:)
    PetscReal, intent(in), pointer :: eff_porosity(:,:)
    PetscReal, intent(in), pointer :: residual_sat(:,:)
    PetscReal, intent(in), pointer :: zwt(:)
    character(len=32), intent(in)  :: vsfm_satfunc_type
    PetscInt, intent(in)           :: density_type
    !
    ! !LOCAL VARIABLES:
    PetscInt                       :: soe_type

    ! Initialize
    call this%Init()

    ! Create all meshes needed by various governing equations
    this%nmesh = 1
    allocate(this%meshes(this%nmesh))

    this%id                  = MPP_VSFM_SNES_CLM
    this%solver_type         = PETSC_SNES
    soe_type                 = SOE_RE_ODE

    call this%meshes(1)%CreateFromCLMCols(begg, endg, begc, endc, &
         ugrid, grc_landunit_indices, lun_coli, lun_colf,         &
         discretization_type, ncols_ghost,                        &
         xc_col, yc_col, zc_col, zi, dz,                          &
         area_col, grid_owner, col_itype,                         &
         filter_vsfmc)

    ! Setup the system-of-equations
    allocate(this%sysofeqns)
    call this%sysofeqns%Setup(mpi_rank, discretization_type, this%id, &
         soe_type, this%meshes, this%nmesh)

    this%sysofeqns%solver_type = this%solver_type

    ! Setup the PETSc
    select case(this%solver_type)
    case (PETSC_SNES)
      call VSFMMPPSetupPetscSNESSetup(this)
    case default
       write(iulog,*) 'VSFMMPPSetup: Unknown this%solver_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    ! Initliaze the VSFM-MPP
    call VSFMMPPInitialize(this, begc, endc, ncols_ghost, &
         zc_col, filter_vsfmc, watsat, hksat, bsw, sucsat, &
         eff_porosity, residual_sat, zwt, vsfm_satfunc_type, density_type)

  end subroutine VSFMMPPSetup

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetupPetscSNESSetup(vsfm_mpp)
    !
    ! !DESCRIPTION:
    ! Sets up the PETSc SNES solver for the VSFM
    !
    use SystemOfEquationsBasePointerType , only : SOEResidual
    use SystemOfEquationsBasePointerType , only : SOEJacobian
    use SystemOfEquationsVSFMType        , only : VSFMSOESetAuxVars
    use GoverningEquationBaseType        , only : goveqn_base_type
    use GoveqnRichardsODEPressureType    , only : goveqn_richards_ode_pressure_type
    use MultiPhysicsProbConstants        , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants        , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants        , only : MPP_VSFM_SNES_CLM
    use mpp_abortutils                   , only : endrun
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                              :: vsfm_mpp
    !
    ! !LOCAL VARIABLES:
    PetscInt                                          :: size
    DM                                                :: dm_p                ! DM object for pressure equation
    PetscErrorCode                                    :: ierr
    class(goveqn_base_type),pointer                   :: cur_goveq
    class (goveqn_richards_ode_pressure_type),pointer :: goveq_richards_pres
    class(sysofeqns_vsfm_type),pointer                :: vsfm_soe
    PetscReal, parameter                              :: atol    = PETSC_DEFAULT_REAL
    PetscReal, parameter                              :: rtol    = PETSC_DEFAULT_REAL
    PetscReal, parameter                              :: stol    = 1.d-10
    PetscInt, parameter                               :: max_it  = PETSC_DEFAULT_INTEGER
    PetscInt, parameter                               :: max_f   = PETSC_DEFAULT_INTEGER

    vsfm_soe => vsfm_mpp%sysofeqns
    vsfm_mpp%sysofeqns_ptr%ptr => vsfm_mpp%sysofeqns

    ! Get pointers to governing-equations

    cur_goveq => vsfm_soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          goveq_richards_pres => cur_goveq
       end select

       cur_goveq => cur_goveq%next
    enddo

    ! DM-Composite approach

    ! Create PETSc DM for pressure-equation
    size = goveq_richards_pres%mesh%ncells_local

    select case(vsfm_mpp%id)
    case (MPP_VSFM_SNES_CLM)
       call DMDACreate1d(PETSC_COMM_SELF,     &
                         DM_BOUNDARY_NONE,    &
                         size,                &
                         1,                   &
                         1,                   &
                         PETSC_NULL_INTEGER,  &
                         dm_p,                &
                         ierr); CHKERRQ(ierr)
       case default
          write(iulog,*)'VSFMMPPSetupPetscSNESSetup: Unknown vsfm_mpp%id'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    call DMSetOptionsPrefix (dm_p , "fp_" , ierr               ); CHKERRQ(ierr)
    call DMDASetFieldName   (dm_p , 0     , "pressure_" , ierr ); CHKERRQ(ierr)
    call DMSetFromOptions   (dm_p , ierr                       ); CHKERRQ(ierr)

    ! Create DMComposite: pressure
    call DMCompositeCreate  (PETSC_COMM_SELF , vsfm_soe%dm , ierr ); CHKERRQ(ierr)
    call DMSetOptionsPrefix (vsfm_soe%dm     , "pressure_" , ierr ); CHKERRQ(ierr)
    call DMCompositeAddDM   (vsfm_soe%dm     , dm_p        , ierr ); CHKERRQ(ierr)
    call DMDestroy          (dm_p            , ierr               ); CHKERRQ(ierr)
    call DMSetUp            (vsfm_soe%dm     , ierr               ); CHKERRQ(ierr)

    call DMCreateMatrix     (vsfm_soe%dm     , vsfm_soe%jac, ierr ); CHKERRQ(ierr)
    call MatSetOption       (vsfm_soe%jac, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)
    call MatSetOption       (vsfm_soe%jac, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)

    ! Create solution vector
    call DMCreateGlobalVector(vsfm_soe%dm , vsfm_soe%soln          , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(vsfm_soe%dm , vsfm_soe%soln_prev     , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(vsfm_soe%dm , vsfm_soe%soln_prev_clm , ierr); CHKERRQ(ierr)

    call VecZeroEntries(vsfm_soe%soln         , ierr); CHKERRQ(ierr)
    call VecZeroEntries(vsfm_soe%soln_prev    , ierr); CHKERRQ(ierr)
    call VecZeroEntries(vsfm_soe%soln_prev_clm, ierr); CHKERRQ(ierr)

    ! SNES
    call SNESCreate(PETSC_COMM_SELF, vsfm_soe%snes, ierr); CHKERRQ(ierr)
    call SNESSetTolerances(vsfm_soe%snes, atol, rtol, stol, &
                           max_it, max_f, ierr); CHKERRQ(ierr)

    call SNESSetFunction(vsfm_soe%snes, PETSC_NULL_OBJECT, SOEResidual, &
                         vsfm_mpp%sysofeqns_ptr, ierr); CHKERRQ(ierr)
    call SNESSetJacobian(vsfm_soe%snes, vsfm_soe%jac, vsfm_soe%jac,     &
                         SOEJacobian, vsfm_mpp%sysofeqns_ptr, ierr); CHKERRQ(ierr)

    call SNESSetFromOptions(vsfm_soe%snes, ierr); CHKERRQ(ierr)

    ! Get pointers to governing-equations
    call vsfm_soe%CreateVectorsForGovEqn()

  end subroutine VSFMMPPSetupPetscSNESSetup

  !------------------------------------------------------------------------
  subroutine VSFMMPPInitialize(vsfm_mpp, begc, endc, ncols_ghost, &
       zc_col, filter_vsfmc, watsat, hksat, bsw, sucsat, eff_porosity, &
       residual_sat, zwt, vsfm_satfunc_type, density_type)

    !
    ! !DESCRIPTION:
    ! Initialize the VSFM by:
    !
    ! - Setting soil properties from CLM,
    ! - Setting initial conditions from CLM,
    ! - Performing pre- and post-solve for the VSFM.
    !
    use SystemOfEquationsBasePointerType , only : SOEResidual
    use SystemOfEquationsBasePointerType , only : SOEJacobian
    use SystemOfEquationsVSFMType        , only : VSFMSOESetAuxVars
    use GoverningEquationBaseType        , only : goveqn_base_type
    use GoveqnRichardsODEPressureType    , only : goveqn_richards_ode_pressure_type
    use MultiPhysicsProbConstants        , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants        , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants        , only : MPP_VSFM_SNES_CLM
    use mpp_abortutils                   , only : endrun
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                              :: vsfm_mpp
    integer, intent(in)                               :: begc,endc
    integer, intent(in)                               :: ncols_ghost          ! number of ghost/halo columns
    PetscReal, pointer, intent(in)                    :: zc_col(:)
    integer, intent(in)                               :: filter_vsfmc(:)
    PetscReal, intent(in), pointer                    :: watsat(:,:)
    PetscReal, intent(in), pointer                    :: hksat(:,:)
    PetscReal, intent(in), pointer                    :: bsw(:,:)
    PetscReal, intent(in), pointer                    :: sucsat(:,:)
    PetscReal, intent(in), pointer                    :: eff_porosity(:,:)
    PetscReal, intent(in), pointer                    :: residual_sat(:,:)
    PetscReal, intent(in), pointer                    :: zwt(:)
    character(len=32), intent(in)                     :: vsfm_satfunc_type
    PetscInt, intent(in)                              :: density_type
    !
    ! !LOCAL VARIABLES:
    PetscInt                                          :: size
    class(goveqn_base_type),pointer                   :: cur_goveq
    class (goveqn_richards_ode_pressure_type),pointer :: goveq_richards_pres
    class(sysofeqns_vsfm_type),pointer                :: vsfm_soe
    Vec                                               :: variable_vec
    PetscReal, pointer                                :: vec_p(:)
    PetscErrorCode                                    :: ierr

    vsfm_soe => vsfm_mpp%sysofeqns

    ! Set initial coniditions
    call VSFMMPPSetSoils(vsfm_mpp, begc, endc, ncols_ghost, filter_vsfmc, &
         watsat, hksat, bsw, sucsat, eff_porosity, residual_sat,          &
         vsfm_satfunc_type, density_type)

    call VSFMMPPSetICs(  vsfm_mpp, begc, endc, zc_col, filter_vsfmc, zwt)

    ! Get pointers to governing-equations
    cur_goveq => vsfm_soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          goveq_richards_pres => cur_goveq
       end select
       cur_goveq => cur_goveq%next
    enddo

    ! Create PETSc DM for pressure-equation
    size = goveq_richards_pres%mesh%ncells_local

    call VecCreateSeq(PETSC_COMM_SELF,size,variable_vec,ierr); CHKERRQ(ierr)
    call VecGetArrayF90(variable_vec, vec_p, ierr); CHKERRQ(ierr)
    vec_p = 273.15d0 + 25.0d0
    call VecRestoreArrayF90(variable_vec, vec_p, ierr); CHKERRQ(ierr)

    call VSFMSOESetAuxVars(vsfm_soe, AUXVAR_INTERNAL, VAR_TEMPERATURE, variable_vec)
    call VecDestroy(variable_vec,ierr); CHKERRQ(ierr)

    ! PreSolve: Allows saturation value to be computed based on ICs and stored
    !           in GE auxvar
    call vsfm_soe%SetDtime(1.d0)
    call vsfm_soe%PreSolve()

    ! PostSolve: Allows saturation value stored in GE auxvar to be copied into
    !            SoE auxvar
    call vsfm_soe%PostSolve()

  end subroutine VSFMMPPInitialize

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetICs(vsfm_mpp, begc, endc, zc_col, filter_vsfmc, zwt)
    !
    ! !DESCRIPTION:
    ! Sets inital conditions for VSFM solver
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : MPP_VSFM_SNES_CLM
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                 :: vsfm_mpp
    integer, intent(in)                  :: begc,endc
    PetscReal, pointer, intent(in)       :: zc_col(:)
    integer, intent(in)                  :: filter_vsfmc(:)
    PetscReal, intent(in), pointer       :: zwt(:)

    select case(vsfm_mpp%id)
    case (MPP_VSFM_SNES_CLM)
       call VSFMMPPSetICsSNESCLM(vsfm_mpp, begc, endc, zc_col, filter_vsfmc, zwt)
    case default
       write(iulog,*) 'VSFMMPPSetICs: Unknown mpp_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine VSFMMPPSetICs

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetSoils(vsfm_mpp, begc, endc, ncols_ghost, filter_vsfmc, &
       watsat, hksat, bsw, sucsat, eff_porosity, residual_sat, &
       vsfm_satfunc_type, density_type)
    !
    ! !DESCRIPTION:
    ! Sets soil properties for VSFM solver
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : MPP_VSFM_SNES_CLM
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)           :: vsfm_mpp
    integer, intent(in)            :: begc,endc
    integer, intent(in)            :: ncols_ghost
    integer, intent(in)            :: filter_vsfmc(:) ! column filter for soil points
    PetscReal, intent(in), pointer :: watsat(:,:)
    PetscReal, intent(in), pointer :: hksat(:,:)
    PetscReal, intent(in), pointer :: bsw(:,:)
    PetscReal, intent(in), pointer :: sucsat(:,:)
    PetscReal, intent(in), pointer :: eff_porosity(:,:)
    PetscReal, intent(in), pointer :: residual_sat(:,:)
    character(len=32), intent(in)  :: vsfm_satfunc_type
    PetscInt                       :: density_type

    select case(vsfm_mpp%id)
    case (MPP_VSFM_SNES_CLM)
       call VSFMMPPSetSoilsCLM(vsfm_mpp, begc, endc, ncols_ghost, filter_vsfmc, &
            watsat, hksat, bsw, sucsat, eff_porosity, residual_sat, &
            vsfm_satfunc_type, density_type)
    case default
       write(iulog,*) 'VSFMMPPSetSoils: Unknown mpp_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine VSFMMPPSetSoils

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetSoilsCLM(vsfm_mpp, begc, endc, ncols_ghost, filter_vsfmc, &
       watsat, hksat, bsw, sucsat, eff_porosity, residual_sat,                   &
       vsfm_satfunc_type, density_type)
    !
    ! !DESCRIPTION:
    ! Sets soil properties for VSFM solver from CLM
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use RichardsODEPressureAuxType    , only : rich_ode_pres_auxvar_type
    use SaturationFunction            , only : SatFunc_Set_BC
    use SaturationFunction            , only : SatFunc_Set_SBC_bz2
    use SaturationFunction            , only : SatFunc_Set_SBC_bz3
    use SaturationFunction            , only : SatFunc_Set_VG
    use PorosityFunctionMod           , only : PorosityFunctionSetConstantModel
    use ConditionType                 , only : condition_type
    use ConnectionSetType             , only : connection_set_type
    use mpp_varpar                    , only : nlevgrnd
    use mpp_varcon                    , only : grav, denh2o
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                              :: vsfm_mpp
    integer, intent(in)                               :: begc,endc
    integer, intent(in)                               :: ncols_ghost
    integer, intent(in)                               :: filter_vsfmc(:)
    PetscReal, intent(in), pointer                    :: watsat(:,:)
    PetscReal, intent(in), pointer                    :: hksat(:,:)
    PetscReal, intent(in), pointer                    :: bsw(:,:)
    PetscReal, intent(in), pointer                    :: sucsat(:,:)
    PetscReal, intent(in), pointer                    :: eff_porosity(:,:)
    PetscReal, intent(in), pointer                    :: residual_sat(:,:)
    character(len=32), intent(in)                     :: vsfm_satfunc_type
    PetscInt                                          :: density_type
    !
    ! !LOCAL VARIABLES:
    class (goveqn_richards_ode_pressure_type),pointer :: goveq_richards_ode_pres
    class(sysofeqns_vsfm_type),pointer                :: vsfm_soe
    class(goveqn_base_type),pointer                   :: cur_goveq
    type (rich_ode_pres_auxvar_type), pointer         :: ode_aux_vars_in(:)
    type (rich_ode_pres_auxvar_type), pointer         :: ode_aux_vars_bc(:)
    type (rich_ode_pres_auxvar_type), pointer         :: ode_aux_vars_ss(:)
    type(condition_type),pointer                      :: cur_cond
    type(connection_set_type), pointer                :: cur_conn_set
    PetscInt                                          :: ghosted_id
    PetscInt                                          :: sum_conn
    PetscInt                                          :: iconn
    PetscInt                                          :: ncells
    PetscInt                                          :: j,c,g,l
    PetscInt                                          :: icell
    PetscReal                                         :: perm
    PetscReal                                         :: por
    PetscReal                                         :: sat_res
    PetscReal                                         :: alpha
    PetscReal                                         :: lambda
    PetscReal, parameter                              :: vish2o = 0.001002d0    ! [N s/m^2] @ 20 degC
    PetscInt                                          :: first_active_hydro_col_id
    PetscInt                                          :: col_id
    PetscBool                                         :: dae_eqn
    PetscInt :: ierr

    first_active_hydro_col_id = -1

    vsfm_soe => vsfm_mpp%sysofeqns

    cur_goveq => vsfm_soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          goveq_richards_ode_pres => cur_goveq
          dae_eqn = PETSC_FALSE
       end select

       cur_goveq => cur_goveq%next
    enddo

    call goveq_richards_ode_pres%SetDensityType(density_type)

    ! Find the first hydrologically active soil column id
    first_active_hydro_col_id = -1
    do c = begc, endc
       if (filter_vsfmc(c) == 1) then
          if (first_active_hydro_col_id == -1) then
             first_active_hydro_col_id = c
             exit
          endif
       endif
    enddo

    if (first_active_hydro_col_id == -1) then
       write(iulog,*)'No active soil hydrology column found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (dae_eqn) then
       write(iulog,*)'DAE formulation of Richards equation is not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    else

       ode_aux_vars_in => goveq_richards_ode_pres%aux_vars_in

       ! In first pass, set hydraulic properties for soil columns that are
       ! active in CLM

       do c = begc, endc + ncols_ghost
          do j = 1, nlevgrnd

             icell = (c - begc)*nlevgrnd + j
             if (filter_vsfmc(c) == 1 ) then
                ! Columns on which hydrology is performed
                col_id = c
             else
                col_id = first_active_hydro_col_id
             endif

             ! perm = hydraulic-conductivity * viscosity / ( density * gravity )
             ![m^2]  [mm/sec]   * [Ns/m^2] /([kg/m^3] * [m/s^2]) * [m/mm]
             perm  = hksat(col_id,j) * vish2o   /(denh2o   * grav   ) * 0.001d0

             ! alpha [1/Pa]; while sucsat [mm of H20]
             ! [Pa] = [mm of H20] * 0.001 [m/mm] * 1000 [kg/m^3] * 9.81 [m/sec^2]
             alpha = 1.d0/(sucsat(col_id,j)*grav)

             ! lambda = 1/bsw
             lambda = 1.d0/bsw(col_id,j)

             sat_res = residual_sat(col_id,j)
             por = watsat(col_id,j)

             ode_aux_vars_in(icell)%perm(1:3) = perm
             ode_aux_vars_in(icell)%por       = por

             call PorosityFunctionSetConstantModel(ode_aux_vars_in(icell)%porParams, &
                  ode_aux_vars_in(icell)%por)

             if (vsfm_satfunc_type == 'brooks_corey') then
                call SatFunc_Set_BC(ode_aux_vars_in(icell)%satParams,      &
                                    sat_res,                               &
                                    alpha,                                 &
                                    lambda)
             elseif (vsfm_satfunc_type == 'smooth_brooks_corey_bz2') then
                call SatFunc_Set_SBC_bz2(ode_aux_vars_in(icell)%satParams, &
                                         sat_res,                          &
                                         alpha,                            &
                                         lambda,                           &
                                         -0.9d0/alpha)
             elseif (vsfm_satfunc_type == 'smooth_brooks_corey_bz3') then
                call SatFunc_Set_SBC_bz3(ode_aux_vars_in(icell)%satParams, &
                                         sat_res,                          &
                                         alpha,                            &
                                         lambda,                           &
                                         -0.9d0/alpha)
             elseif (vsfm_satfunc_type == 'van_genuchten') then
                call SatFunc_Set_VG(ode_aux_vars_in(icell)%satParams,      &
                                    sat_res,                               &
                                    alpha,                                 &
                                    lambda)
             else
                call endrun(msg='ERROR:: Unknown vsfm_satfunc_type = '//vsfm_satfunc_type//&
                  errMsg(__FILE__, __LINE__))
             endif

          enddo
       enddo

       ! Set soil properties for boundary-condition auxvars
       sum_conn = 0
       ode_aux_vars_bc => goveq_richards_ode_pres%aux_vars_bc
       cur_cond        => goveq_richards_ode_pres%boundary_conditions%first

       do
          if (.not.associated(cur_cond)) exit
          cur_conn_set => cur_cond%conn_set

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             ghosted_id = cur_conn_set%id_dn(iconn)

             ode_aux_vars_bc(sum_conn)%perm(:)       = ode_aux_vars_in(ghosted_id)%perm(:)
             ode_aux_vars_bc(sum_conn)%por           = ode_aux_vars_in(ghosted_id)%por
             ode_aux_vars_bc(sum_conn)%satParams     = ode_aux_vars_in(ghosted_id)%satParams
             ode_aux_vars_bc(sum_conn)%porParams     = ode_aux_vars_in(ghosted_id)%porParams

             ode_aux_vars_bc(sum_conn)%pressure_prev = 3.5355d3

             call ode_aux_vars_bc(sum_conn)%satParams%Copy(ode_aux_vars_in(ghosted_id)%satParams)

          enddo
          cur_cond => cur_cond%next
       enddo

       ! Set soil properties for source-sink auxvars
       sum_conn = 0

       ode_aux_vars_ss => goveq_richards_ode_pres%aux_vars_ss
       cur_cond        => goveq_richards_ode_pres%source_sinks%first

       do
          if (.not.associated(cur_cond)) exit
          cur_conn_set => cur_cond%conn_set

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             ghosted_id = cur_conn_set%id_dn(iconn)

             ode_aux_vars_ss(sum_conn)%perm(:)       = ode_aux_vars_in(ghosted_id)%perm(:)
             ode_aux_vars_ss(sum_conn)%por           = ode_aux_vars_in(ghosted_id)%por
             ode_aux_vars_ss(sum_conn)%satParams     = ode_aux_vars_in(ghosted_id)%satParams
             ode_aux_vars_ss(sum_conn)%porParams     = ode_aux_vars_in(ghosted_id)%porParams

             ode_aux_vars_ss(sum_conn)%pressure_prev = 3.5355d3

          enddo
          cur_cond => cur_cond%next
       enddo

    endif

  end subroutine VSFMMPPSetSoilsCLM

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetICsSNESCLM(vsfm_mpp, begc, endc, zc_col, filter_vsfmc, zwt)
    !
    ! !DESCRIPTION:
    ! Sets inital conditions for VSFM solver from CLM
    !
    ! !USES:
    use SystemOfEquationsBasePointerType , only : SOEIFunction, SOEIJacobian
    use GoverningEquationBaseType        , only : goveqn_base_type
    use GoveqnRichardsODEPressureType    , only : goveqn_richards_ode_pressure_type
    use mpp_varpar                       , only : nlevgrnd
    use MultiPhysicsProbConstants        , only : GRAVITY_CONSTANT
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                              :: vsfm_mpp
    integer, intent(in)                               :: begc,endc
    PetscReal, pointer, intent(in)                    :: zc_col(:)
    integer, intent(in)                               :: filter_vsfmc(:)
    PetscReal, intent(in), pointer                    :: zwt(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                          :: size
    PetscErrorCode                                    :: ierr
    class(goveqn_base_type),pointer                   :: cur_goveq
    class (goveqn_richards_ode_pressure_type),pointer :: goveq_richards_ode_pres
    class(sysofeqns_vsfm_type),pointer                :: vsfm_soe
    PetscInt                                          :: nDM
    DM, pointer                                       :: dms(:)
    Vec, pointer                                      :: soln_subvecs(:)
    PetscReal, pointer                                :: press_p(:)
    PetscInt                                          :: ghosted_id
    PetscReal                                         :: initial_sat, initial_pressure
    PetscInt                                          :: j,c,g,l
    PetscInt                                          :: first_active_hydro_col_id
    PetscInt                                          :: icell
    PetscInt                                          :: col_id

    first_active_hydro_col_id = -1
    vsfm_soe => vsfm_mpp%sysofeqns

    cur_goveq => vsfm_soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          goveq_richards_ode_pres => cur_goveq
          exit
       end select
       cur_goveq => cur_goveq%next
    enddo

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(vsfm_soe%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(vsfm_soe%dm, dms, ierr)

    ! Allocate vectors for individual GEs
    allocate(soln_subvecs(nDM))

    ! Get solution vectors for individual GEs
    call DMCompositeGetAccessArray(vsfm_soe%dm, vsfm_soe%soln, nDM, &
                                   PETSC_NULL_INTEGER, soln_subvecs, ierr)

    call VecGetArrayF90(soln_subvecs(1), press_p, ierr)

    first_active_hydro_col_id = -1
    do c = begc, endc
       if (filter_vsfmc(c) == 1) then
          if (first_active_hydro_col_id == -1) then
             first_active_hydro_col_id = c
             exit
          endif
       endif
    enddo

    if (first_active_hydro_col_id == -1) then
       write(iulog,*)'No active soil hydrology column found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    icell = 0
    do c = begc, endc
       do j = 1, nlevgrnd

          icell = (c - begc)*nlevgrnd + j
          if (filter_vsfmc(c) == 1) then

             ! Columns on which hydrology is performed
             col_id = c
          else
             col_id = first_active_hydro_col_id
          endif

          press_p(icell) = 101325.d0 + &
                           997.16d0*GRAVITY_CONSTANT * &
                           (-zwt(col_id) - (goveq_richards_ode_pres%mesh%z(icell) - zc_col(c)) )
       enddo
    enddo

    call VecRestoreArrayF90(soln_subvecs(1), press_p, ierr)

    ! Restore solution vectors for individual GEs
    call DMCompositeRestoreAccessArray(vsfm_soe%dm, vsfm_soe%soln, nDM, &
                                       PETSC_NULL_INTEGER, soln_subvecs, ierr)

    call VecCopy(vsfm_soe%soln, vsfm_soe%soln_prev,     ierr); CHKERRQ(ierr)
    call VecCopy(vsfm_soe%soln, vsfm_soe%soln_prev_clm, ierr); CHKERRQ(ierr)

    ! Free memory
    deallocate(dms)
    deallocate(soln_subvecs)

  end subroutine VSFMMPPSetICsSNESCLM

  !------------------------------------------------------------------------
  subroutine VSFMMPPUpdateSysOfEqnsAuxVars(this, gov_eqn_id, auxvar_type, &
                                           var_id, condition_id, num_values, &
                                           values)
    !
    ! !DESCRIPTION:
    ! Updates SoE auxvars associated with:
    ! - Internal connections,
    ! - Boundary conditions, and
    ! - Source-sink terms.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants, only : AUXVAR_BC
    use MultiPhysicsProbConstants, only : AUXVAR_SS
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)               :: this
    PetscInt, intent(in)               :: gov_eqn_id
    PetscInt, intent(in)               :: auxvar_type
    PetscInt, intent(in)               :: var_id
    PetscInt, intent(in)               :: condition_id
    PetscInt, intent(in)               :: num_values
    PetscReal,intent(in)               :: values(:)
    ! !LOCAL VARIABLES:
    class(sysofeqns_vsfm_type),pointer :: vsfm_soe
    PetscInt                           :: iauxvar
    PetscInt                           :: counter

    vsfm_soe => this%sysofeqns

    select case(auxvar_type)
    case (AUXVAR_INTERNAL)
       write(iulog,*)'VSFMMPPUpdateSysOfEqnsAuxVars: Add code for AUXVAR_INTERNAL'
       call endrun(msg=errMsg(__FILE__, __LINE__))

    case (AUXVAR_BC)
       counter = 0
       do iauxvar = 1, vsfm_soe%num_auxvars_bc
          if (vsfm_soe%aux_vars_bc(iauxvar)%is_bc                         .and. &
              vsfm_soe%aux_vars_bc(iauxvar)%goveqn_id == gov_eqn_id       .and. &
              vsfm_soe%aux_vars_bc(iauxvar)%condition_id == condition_id        &
              ) then

             counter = counter + 1
             vsfm_soe%aux_vars_bc(iauxvar)%condition_value = values(counter)

          end if
       enddo

       if (counter /= num_values) then
          write(iulog,*)'VSFMMPPUpdateSysOfEqnsAuxVars: All BC values not assigned to ' // &
             'vsfm_soe%aux_vars_bc'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

    case (AUXVAR_SS)
       counter = 0
       do iauxvar = 1, vsfm_soe%num_auxvars_ss
          if (vsfm_soe%aux_vars_ss(iauxvar)%is_ss                         .and. &
              vsfm_soe%aux_vars_ss(iauxvar)%goveqn_id == gov_eqn_id       .and. &
              vsfm_soe%aux_vars_ss(iauxvar)%condition_id == condition_id        &
              ) then

             counter = counter + 1
             vsfm_soe%aux_vars_ss(iauxvar)%condition_value = values(counter)

          end if
       enddo

       if (counter /= num_values) then
          write(iulog,*)'VSFMMPPUpdateSysOfEqnsAuxVars: All BC values not assigned to ' // &
             'vsfm_soe%aux_vars_ss'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

    case default
       write(iulog,*) 'VSFMMPPUpdateSysOfEqnsAuxVars: Unknown auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine VSFMMPPUpdateSysOfEqnsAuxVars

  !------------------------------------------------------------------------
  subroutine VSFMMPPRestart(this, data_1d)
    !
    ! !DESCRIPTION:
    ! Updates PETSc vectors that store the solution at previous time step 
    ! (%soln_prev) and the first guess at next time step (%soln).
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                              :: this
    PetscReal                                         :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                          :: size
    PetscErrorCode                                    :: ierr
    class(goveqn_base_type),pointer                   :: cur_goveq
    class (goveqn_richards_ode_pressure_type),pointer :: goveq_richards_pres
    class(sysofeqns_vsfm_type),pointer                :: vsfm_soe
    PetscInt                                          :: nDM
    DM, pointer                                       :: dms(:)
    Vec, pointer                                      :: soln_subvecs(:)
    PetscReal, pointer                                :: press_p(:)
    PetscInt                                          :: local_id

    vsfm_soe => this%sysofeqns

    cur_goveq => vsfm_soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          goveq_richards_pres => cur_goveq
          exit
       end select
       cur_goveq => cur_goveq%next
    enddo

    if (size(data_1d) /= goveq_richards_pres%mesh%ncells_local) then
       write(iulog,*) 'VSFMMPPRestart: size(data_1d) /= goveq_richards_pres%mesh%ncells_local'
       write(iulog,*) 'size(data_1d)                    = ',size(data_1d)
       write(iulog,*) 'goveq_richards_pres%mesh%ncells_local  = ', goveq_richards_pres%mesh%ncells_local
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(vsfm_soe%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(vsfm_soe%dm, dms, ierr)

    ! Allocate vectors for individual GEs
    allocate(soln_subvecs(nDM))

    ! Get solution vectors for individual GEs
    call DMCompositeGetAccessArray(vsfm_soe%dm, vsfm_soe%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    call VecGetArrayF90(soln_subvecs(1), press_p, ierr)
    do local_id = 1, goveq_richards_pres%mesh%ncells_local
       press_p(local_id) = data_1d(local_id)
    enddo
    call VecRestoreArrayF90(soln_subvecs(1), press_p, ierr)

    ! Restore solution vectors for individual GEs
    call DMCompositeRestoreAccessArray(vsfm_soe%dm, vsfm_soe%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    call VecCopy(vsfm_soe%soln, vsfm_soe%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(vsfm_soe%soln, vsfm_soe%soln_prev_clm, ierr); CHKERRQ(ierr)

    ! Free up memory
    deallocate(dms)
    deallocate(soln_subvecs)

  end subroutine VSFMMPPRestart

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetMeshesOfGoveqns(this)
    !
    ! !DESCRIPTION:
    ! Set association of governing equations and meshes
    !
    use GoverningEquationBaseType, only : goveqn_base_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type) :: this

    call this%sysofeqns%SetMeshesOfGoveqns(this%meshes, this%nmesh)

  end subroutine VSFMMPPSetMeshesOfGoveqns

  !------------------------------------------------------------------------
  subroutine VSFMMPPAddGovEqn(this, geq_type, name, mesh_itype)
    !
    ! !DESCRIPTION:
    ! Adds a governing equation to the MPP
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type) :: this
    PetscInt             :: geq_type
    character(len =*)    :: name
    PetscInt             :: mesh_itype

    call this%sysofeqns%AddGovEqn(geq_type, name, mesh_itype)

  end subroutine VSFMMPPAddGovEqn


  !------------------------------------------------------------------------
  subroutine VSFMMPPGovEqnAddCondition(this, igoveqn, ss_or_bc_type, name, unit, &
       cond_type, region_type, id_of_other_goveq, conn_set)
    !
    ! !DESCRIPTION:
    ! Adds a boundary/source-sink condition to a governing equation
    !
    use GoverningEquationBaseType, only : goveqn_base_type
    use ConnectionSetType        , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)              :: this
    PetscInt                          :: igoveqn
    PetscInt                          :: ss_or_bc_type
    character(len =*)                 :: name
    character(len =*)                 :: unit
    PetscInt                          :: cond_type
    PetscInt                          :: region_type
    PetscInt, optional                :: id_of_other_goveq
    type(connection_set_type),pointer, optional :: conn_set
    !
    class(goveqn_base_type),pointer   :: cur_goveq
    class(goveqn_base_type),pointer   :: other_goveq
    PetscInt                          :: ii

    if (igoveqn > this%sysofeqns%ngoveqns) then
       write(iulog,*) 'Attempting to add condition for governing equation ' // &
            'that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq => this%sysofeqns%goveqns
    do ii = 1, igoveqn-1
       cur_goveq => cur_goveq%next
    enddo

    if (.not.present(id_of_other_goveq)) then
       if (.not.present(conn_set)) then
          call cur_goveq%AddCondition(ss_or_bc_type, name, unit, &
               cond_type, region_type)
       else
          call cur_goveq%AddCondition(ss_or_bc_type, name, unit, &
               cond_type, region_type, conn_set=conn_set)
       endif
    else

       other_goveq => this%sysofeqns%goveqns
       do ii = 1,id_of_other_goveq-1
          other_goveq => other_goveq%next
       enddo

       if (.not.present(conn_set)) then
          call cur_goveq%AddCondition(ss_or_bc_type, name, unit, &
               cond_type, region_type, id_of_other_goveq, other_goveq%id )
       else
          call cur_goveq%AddCondition(ss_or_bc_type, name, unit, &
               cond_type, region_type, id_of_other_goveq=id_of_other_goveq, &
               itype_of_other_goveq = other_goveq%id, conn_set=conn_set)
       endif
    endif

  end subroutine VSFMMPPGovEqnAddCondition

  !------------------------------------------------------------------------
  subroutine VSFMMPPGovEqnUpdateConditionConnSet(this, igoveqn, icond, &
       ss_or_bc_type, nconn,  conn_id_up, conn_id_dn, &
       conn_dist_up, conn_dist_dn,  conn_unitvec, conn_area)
    !
    ! !DESCRIPTION:
    ! For a given governing equation, updates connection set of a
    ! boundary/source-sink condition
    !
    use GoverningEquationBaseType, only : goveqn_base_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)            :: this
    PetscInt                        :: igoveqn
    PetscInt                        :: icond
    PetscInt                        :: ss_or_bc_type
    PetscInt                        :: nconn
    PetscInt, pointer               :: conn_id_up(:)   !
    PetscInt, pointer               :: conn_id_dn(:)   !
    PetscReal, pointer              :: conn_dist_up(:) !
    PetscReal, pointer              :: conn_dist_dn(:) !
    PetscReal, pointer              :: conn_area(:)    !
    PetscReal, pointer              :: conn_type(:)    !
    PetscReal, pointer              :: conn_unitvec(:,:)
    !
    class(goveqn_base_type),pointer :: cur_goveq
    class(goveqn_base_type),pointer :: other_goveq
    PetscInt                        :: ii

    if (igoveqn > this%sysofeqns%ngoveqns) then
       write(iulog,*) 'Attempting to add condition for governing equation ' // &
            'that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq => this%sysofeqns%goveqns
    do ii = 1, igoveqn-1
       cur_goveq => cur_goveq%next
    enddo

    call cur_goveq%UpdateConditionConnSet(icond, ss_or_bc_type, nconn, &
         conn_id_up, conn_id_dn, conn_dist_up, conn_dist_dn,    &
         conn_unitvec, conn_area)

  end subroutine VSFMMPPGovEqnUpdateConditionConnSet

  !------------------------------------------------------------------------
  subroutine VSFMMPPAllocateAuxVars(this)
    !
    ! !DESCRIPTION:
    ! Allocates auxvars for governing equations and system-of-governing-eqns
    !
    use SystemOfEquationsBaseType           , only : sysofeqns_base_type
    use GoverningEquationBaseType           , only : goveqn_base_type
    use GoveqnRichardsODEPressureType       , only : goveqn_richards_ode_pressure_type
    use MultiPhysicsProbConstants           , only : COND_BC
    use MultiPhysicsProbConstants           , only : COND_SS
    use MultiPhysicsProbConstants           , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                   :: this
    !
    class(sysofeqns_base_type), pointer    :: soe_base
    class(sysofeqns_vsfm_type), pointer :: soe
    class(goveqn_base_type), pointer       :: cur_goveq
    PetscInt                               :: igoveqn
    PetscInt                               :: num_bc
    PetscInt                               :: num_ss
    PetscInt                               :: icond
    PetscInt                               :: iauxvar
    PetscInt                               :: iauxvar_beg, iauxvar_end
    PetscInt                               :: iauxvar_beg_bc, iauxvar_end_bc
    PetscInt                               :: iauxvar_beg_ss, iauxvar_end_ss
    PetscInt                               :: count_bc, count_ss
    PetscInt                               :: offset_bc, offset_ss
    PetscInt, pointer                      :: ncells_for_bc(:)
    PetscInt, pointer                      :: ncells_for_ss(:)
    PetscInt, pointer                      :: offsets_bc(:)
    PetscInt, pointer                      :: offsets_ss(:)

    soe_base => this%sysofeqns

    select type(soe_base)
    class is(sysofeqns_vsfm_type)
       soe => this%sysofeqns
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    !
    ! Pass-1: Determine total number of BCs (excluding BCs
    !         needed for coupling various governing equations)
    !         and SSs for all governing equations
    !
    igoveqn = 0
    cur_goveq => soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          call cur_goveq%AllocateAuxVars()
          call cur_goveq%GetNumCellsInConditions(COND_BC, &
               COND_DIRICHLET_FRM_OTR_GOVEQ, num_bc, ncells_for_bc)
          call cur_goveq%GetNumCellsInConditions(COND_SS, -9999, &
               num_ss, ncells_for_ss)

       end select

       igoveqn = igoveqn + 1

       soe%num_auxvars_in = soe%num_auxvars_in + &
            cur_goveq%mesh%ncells_all

       do icond = 1, num_bc
          soe%num_auxvars_bc = soe%num_auxvars_bc + &
               ncells_for_bc(icond)
       enddo

       do icond = 1, num_ss
          soe%num_auxvars_ss = soe%num_auxvars_ss + &
               ncells_for_ss(icond)
       enddo

       if (num_bc > 0) deallocate(ncells_for_bc)
       if (num_ss > 0) deallocate(ncells_for_ss)

       cur_goveq => cur_goveq%next
    enddo

    ! Allocate memory
    allocate(soe%aux_vars_in           (soe%num_auxvars_in))
    allocate(soe%aux_vars_bc           (soe%num_auxvars_bc))
    allocate(soe%aux_vars_ss           (soe%num_auxvars_ss))

    allocate(soe%soe_auxvars_bc_offset (soe%num_auxvars_bc))
    allocate(soe%soe_auxvars_ss_offset (soe%num_auxvars_bc))
    allocate(soe%soe_auxvars_bc_ncells (soe%num_auxvars_bc))
    allocate(soe%soe_auxvars_ss_ncells (soe%num_auxvars_ss))

    igoveqn        = 0
    iauxvar_beg    = 0
    iauxvar_end    = 0
    iauxvar_beg_bc = 0
    iauxvar_end_bc = 0
    iauxvar_beg_ss = 0
    iauxvar_end_ss = 0
    count_bc       = 0
    count_ss       = 0
    offset_bc      = 0
    offset_ss      = 0

    !
    ! Pass-2: Set values for auxvars of SoE
    !
    cur_goveq => soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          call cur_goveq%GetNumCellsInConditions(COND_BC, &
               COND_DIRICHLET_FRM_OTR_GOVEQ, num_bc, ncells_for_bc)
          call cur_goveq%GetNumCellsInConditions(COND_SS, -9999, &
               num_ss, ncells_for_ss)

       end select

       igoveqn = igoveqn + 1

       iauxvar_beg = iauxvar_end + 1
       iauxvar_end = iauxvar_end + cur_goveq%mesh%ncells_all

       do iauxvar = iauxvar_beg, iauxvar_end
          call soe%aux_vars_in(iauxvar)%Init()

          soe%aux_vars_in(iauxvar)%is_in     = PETSC_TRUE
          soe%aux_vars_in(iauxvar)%goveqn_id = igoveqn
       enddo

       allocate(offsets_bc(num_bc))
       allocate(offsets_ss(num_ss))

       allocate(soe%soe_auxvars_bc_offset(num_bc))
       allocate(soe%soe_auxvars_ss_offset(num_ss))
       allocate(soe%soe_auxvars_bc_ncells(num_bc))
       allocate(soe%soe_auxvars_ss_ncells(num_ss))

       do icond = 1, num_bc
          count_bc = count_bc + 1

          soe%soe_auxvars_bc_offset(count_bc) = offset_bc
          offsets_bc(icond)                   = offset_bc
          soe%soe_auxvars_bc_ncells(count_bc) = ncells_for_bc(icond)
          offset_bc                           = offset_bc + ncells_for_bc(icond)

          iauxvar_beg_bc = iauxvar_end_bc + 1
          iauxvar_end_bc = iauxvar_end_bc + ncells_for_bc(icond)

          do iauxvar = iauxvar_beg_bc, iauxvar_end_bc
             call soe%aux_vars_bc(iauxvar)%Init()

             soe%aux_vars_bc(iauxvar)%is_bc        = PETSC_TRUE
             soe%aux_vars_bc(iauxvar)%goveqn_id    = igoveqn
             soe%aux_vars_bc(iauxvar)%condition_id = icond
          enddo
       enddo

       do icond = 1, num_ss
          count_ss = count_ss + 1
          soe%soe_auxvars_ss_offset(count_ss) = offset_ss
          offsets_ss(icond)                   = offset_ss
          soe%soe_auxvars_ss_ncells(count_ss) = ncells_for_ss(icond)
          offset_ss                           = offset_ss + ncells_for_ss(icond)

          iauxvar_beg_ss = iauxvar_end_ss + 1
          iauxvar_end_ss = iauxvar_end_ss + ncells_for_ss(icond)

          do iauxvar = iauxvar_beg_ss, iauxvar_end_ss
             call soe%aux_vars_ss(iauxvar)%Init()

             soe%aux_vars_ss(iauxvar)%is_ss        = PETSC_TRUE
             soe%aux_vars_ss(iauxvar)%goveqn_id    = igoveqn
             soe%aux_vars_ss(iauxvar)%condition_id = icond
          enddo
       enddo

       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          call cur_goveq%SetSOEAuxVarOffsets(num_bc, offsets_bc, num_ss, offsets_ss)
       end select

       if (num_bc > 0) deallocate(ncells_for_bc)
       if (num_ss > 0) deallocate(ncells_for_ss)

       deallocate(offsets_bc)
       deallocate(offsets_ss)

       cur_goveq => cur_goveq%next
    enddo

  end subroutine VSFMMPPAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine VSFMMPPGovEqnSetCouplingVars(this, igoveqn, nvars, &
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
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)              :: this
    PetscInt                          :: igoveqn
    PetscInt                          :: nvars
    PetscInt, pointer                 :: var_ids(:)
    PetscInt, pointer                 :: goveqn_ids(:)
    !
    class(goveqn_base_type) , pointer :: cur_goveq_1
    class(goveqn_base_type) , pointer :: cur_goveq_2
    type(condition_type)    , pointer :: cur_cond_1
    type(condition_type)    , pointer :: cur_cond_2
    PetscInt                          :: ii
    PetscInt                          :: ieqn
    PetscInt                          :: ivar
    PetscInt                          :: bc_idx_1
    PetscInt                          :: bc_idx_2
    PetscInt                          :: bc_offset_1
    PetscBool                         :: bc_found

    if (igoveqn > this%sysofeqns%ngoveqns) then
       write(iulog,*) 'Attempting to set coupling vars for governing ' // &
            'equation that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq_1 => this%sysofeqns%goveqns
    do ii = 1, igoveqn-1
       cur_goveq_1 => cur_goveq_1%next
    end do

    call cur_goveq_1%AllocVarsFromOtherGEs(nvars)

    do ivar = 1, nvars

       if (goveqn_ids(ivar) > this%sysofeqns%ngoveqns) then
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

       cur_goveq_2 => this%sysofeqns%goveqns
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

       cur_goveq_1%var_ids_needed_from_other_goveqns (ivar) = var_ids(ivar)
       cur_goveq_1%ids_of_other_goveqns              (ivar) = goveqn_ids(ivar)
       cur_goveq_1%is_bc_auxvar_type                 (ivar) = PETSC_TRUE
       cur_goveq_1%bc_auxvar_offset                  (ivar) = bc_offset_1
       cur_goveq_1%bc_auxvar_ncells                  (ivar) = cur_cond_1%conn_set%num_connections
       cur_goveq_1%bc_auxvar_idx                     (ivar) = bc_idx_1
       cur_goveq_1%bc_auxvar_idx_of_other_goveqn     (ivar) = bc_idx_2

    enddo

  end subroutine VSFMMPPGovEqnSetCouplingVars

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetupProblem(this)
    !
    ! !DESCRIPTION:
    ! Sets up the PETSc related data structure
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : SOE_RE_ODE
    use MultiPhysicsProbConstants , only : PETSC_SNES
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type) :: this

    select case(this%solver_type)
    case (PETSC_SNES)
      call VSFMMPPSetupPetscSNESSetup(this)
    case default
       write(iulog,*) 'VSFMMPPSetup: Unknown this%solver_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    this%sysofeqns%solver_type = this%solver_type
    this%sysofeqns%itype       = SOE_RE_ODE

  end subroutine VSFMMPPSetupProblem

  !------------------------------------------------------------------------
  subroutine VSFMMPPClean(this)
    !
    ! !DESCRIPTION:
    ! Release memory allocated to the object
    !
    ! !USES:
    use MultiPhysicsProbBaseType , only : MMPBaseClean
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type) :: this

    call MMPBaseClean(this)

  end subroutine VSFMMPPClean

#endif

end module MultiPhysicsProbVSFM
