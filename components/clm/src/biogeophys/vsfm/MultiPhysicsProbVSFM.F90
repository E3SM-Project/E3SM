#ifdef USE_PETSC_LIB


module MultiPhysicsProbVSFM
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Object for variably saturated flow model
  !-----------------------------------------------------------------------

  ! !USES:
  use SoilStateType                      , only : soilstate_type
  use WaterstateType                     , only : waterstate_type
  use SoilHydrologyType                  , only : soilhydrology_type
  use clm_varctl                         , only : iulog
  use abortutils                         , only : endrun
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
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
     procedure, public :: Init                   => VSFMMPPInit
     procedure, public :: Clean                  => VSFMMPPClean
     procedure, public :: Setup                  => VSFMMPPSetup
     procedure, public :: Restart                => VSFMMPPRestart
     procedure, public :: UpdateSysOfEqnsAuxVars => VSFMMPPUpdateSysOfEqnsAuxVars
  end type mpp_vsfm_type

  type(mpp_vsfm_type), public, target :: vsfm_mpp

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

    this%name   = 'Variably-Saturated-Flow-Model'
    this%id     = MPP_VSFM_SNES_CLM
    this%nmesh  = 1

    nullify(this%sysofeqns)
    allocate(this%sysofeqns_ptr)
    nullify(this%sysofeqns_ptr%ptr)

  end subroutine VSFMMPPInit

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetup(this, begc, endc, num_hydrologyc, filter_hydrologyc, &
       soilstate_vars, waterstate_vars, soilhydrology_vars)
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
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                 :: this
    integer, intent(in)                  :: begc,endc
    integer, intent(in)                  :: num_hydrologyc       ! number of column soil points in column filter
    integer, intent(in)                  :: filter_hydrologyc(:) ! column filter for soil points
    type(soilstate_type) , intent(in)    :: soilstate_vars
    type(waterstate_type)  , intent(in)  :: waterstate_vars
    type(soilhydrology_type), intent(in) :: soilhydrology_vars
    !
    ! !LOCAL VARIABLES:
    PetscInt                             :: soe_type

    ! Initialize
    call this%Init()

    ! Create all meshes needed by various governing equations
    this%nmesh = 1
    allocate(this%meshes(this%nmesh))

    this%id           = MPP_VSFM_SNES_CLM
    this%solver_type  = PETSC_SNES
    soe_type          = SOE_RE_ODE
    call this%meshes(1)%Create(MESH_CLM_SOIL_COL, begc, endc, &
                               waterstate_vars, soilhydrology_vars)

    ! Setup the system-of-equations
    allocate(this%sysofeqns)
    call this%sysofeqns%Setup(this%id, soe_type, this%meshes, this%nmesh)

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
    call VSFMMPPInitialize(this, begc, endc, soilstate_vars, &
             waterstate_vars, soilhydrology_vars)

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
    use abortutils                       , only : endrun
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
    size = goveq_richards_pres%mesh%ncells

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

  end subroutine VSFMMPPSetupPetscSNESSetup

  !------------------------------------------------------------------------
  subroutine VSFMMPPInitialize(vsfm_mpp, begc, endc, &
       soilstate_vars, waterstate_vars, soilhydrology_vars)

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
    use abortutils                       , only : endrun
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                              :: vsfm_mpp
    integer, intent(in)                               :: begc,endc
    type(soilstate_type) , intent(in)                 :: soilstate_vars
    type(waterstate_type), intent(in)                 :: waterstate_vars
    type(soilhydrology_type), intent(in)              :: soilhydrology_vars
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
    call VSFMMPPSetSoils(vsfm_mpp, begc, endc, soilstate_vars     )
    call VSFMMPPSetICs(  vsfm_mpp, begc, endc, soilhydrology_vars )

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
    size = goveq_richards_pres%mesh%ncells

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
  subroutine VSFMMPPSetICs(vsfm_mpp, begc, endc, soilhydrology_vars)
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
    type(soilhydrology_type), intent(in) :: soilhydrology_vars

    select case(vsfm_mpp%id)
    case (MPP_VSFM_SNES_CLM)
       call VSFMMPPSetICsSNESCLM(vsfm_mpp, begc, endc, soilhydrology_vars)
    case default
       write(iulog,*) 'VSFMMPPSetICs: Unknown mpp_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine VSFMMPPSetICs

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetSoils(vsfm_mpp, begc, endc, soilstate_vars)
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
    class(mpp_vsfm_type)              :: vsfm_mpp
    integer, intent(in)               :: begc,endc
    type(soilstate_type) , intent(in) :: soilstate_vars

    select case(vsfm_mpp%id)
    case (MPP_VSFM_SNES_CLM)
       call VSFMMPPSetSoilsCLM(vsfm_mpp, begc, endc, soilstate_vars)
    case default
       write(iulog,*) 'VSFMMPPSetSoils: Unknown mpp_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine VSFMMPPSetSoils

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetSoilsCLM(vsfm_mpp, begc, endc, soilstate_vars)
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
    use clm_varpar                    , only : nlevgrnd
    use ColumnType                    , only : col
    use LandunitType                  , only : lun
    use landunit_varcon               , only : istcrop, istsoil
    use column_varcon                 , only : icol_road_perv
    use clm_varcon                    , only : grav, denh2o
    use clm_varctl                    , only : vsfm_satfunc_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                              :: vsfm_mpp
    integer, intent(in)                               :: begc,endc
    type(soilstate_type) , intent(in)                 :: soilstate_vars
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

    associate(&
         smpmin       =>    soilstate_vars%smpmin_col       , & ! Input:  [real(r8) (:)   ]  restriction for min of soil potential (mm)
         watsat       =>    soilstate_vars%watsat_col       , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         hksat        =>    soilstate_vars%hksat_col        , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
         bsw          =>    soilstate_vars%bsw_col          , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
         sucsat       =>    soilstate_vars%sucsat_col       , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
         eff_porosity =>    soilstate_vars%eff_porosity_col , & ! Input:  [real(r8) (:,:) ]  effective porosity = porosity - vol_ice
         rootr_col    =>    soilstate_vars%rootr_col        , & ! Input:  [real(r8) (:,:) ]  effective fraction of roots in each soil layer
         smp_l        =>    soilstate_vars%smp_l_col        , & ! Input:  [real(r8) (:,:) ]  soil matrix potential [mm]
         hk_l         =>    soilstate_vars%hk_l_col         , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity (mm/s)
         rootr_pft    =>    soilstate_vars%rootr_patch        & ! Input:  [real(r8) (:,:) ]  effective fraction of roots in each soil layer
         )

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

    ! Find the first hydrologically active soil column id
    first_active_hydro_col_id = -1
    do c = begc, endc
       l = col%landunit(c)

       if (col%active(c) .and. &
           (lun%itype(l) == istsoil .or. col%itype(c) == icol_road_perv .or. &
            lun%itype(l) == istcrop)) then
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

       !icell = 0
       do c = begc, endc
          do j = 1, nlevgrnd

             !icell = icell + 1
             icell = (c - begc)*nlevgrnd + j
             g = col%gridcell(c)
             l = col%landunit(c)

             if (col%active(c) .and. &
                 (lun%itype(l) == istsoil .or. col%itype(c) == icol_road_perv .or. &
                  lun%itype(l) == istcrop)) then
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

             sat_res = 0.d0
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

    end associate

  end subroutine VSFMMPPSetSoilsCLM

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetICsSNESCLM(vsfm_mpp, begc, endc, soilhydrology_vars)
    !
    ! !DESCRIPTION:
    ! Sets inital conditions for VSFM solver from CLM
    !
    ! !USES:
    use SystemOfEquationsBasePointerType , only : SOEIFunction, SOEIJacobian
    use GoverningEquationBaseType        , only : goveqn_base_type
    use GoveqnRichardsODEPressureType    , only : goveqn_richards_ode_pressure_type
    use clm_varpar                       , only : nlevgrnd, nlevsoi
    use ColumnType                       , only : col
    use LandunitType                     , only : lun
    use column_varcon                    , only : icol_road_perv, icol_road_imperv
    use landunit_varcon                  , only : istice, istwet, istice_mec, istcrop, istsoil
    use clm_varcon                       , only : grav, denh2o
    use MultiPhysicsProbConstants        , only : GRAVITY_CONSTANT
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                              :: vsfm_mpp
    integer, intent(in)                               :: begc,endc
    type(soilhydrology_type), intent(in)              :: soilhydrology_vars
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

    associate(                                 &
         zwt_col => soilhydrology_vars%zwt_col &
         )

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
       l = col%landunit(c)

       if (col%active(c) .and. &
           (lun%itype(l) == istsoil .or. col%itype(c) == icol_road_perv .or. &
            lun%itype(l) == istcrop)) then
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
          g = col%gridcell(c)
          l = col%landunit(c)

          if (col%active(c) .and. &
              (lun%itype(l) == istsoil .or. col%itype(c) == icol_road_perv .or. &
                lun%itype(l) == istcrop)) then
             ! Columns on which hydrology is performed
             col_id = c
          else
             col_id = first_active_hydro_col_id
          endif

          press_p(icell) = 101325.d0 + &
                           997.16d0*GRAVITY_CONSTANT*(-zwt_col(col_id) - goveq_richards_ode_pres%mesh%z(icell))
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

    end associate

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

    if (size(data_1d) /= goveq_richards_pres%mesh%ncells) then
       write(iulog,*) 'VSFMMPPRestart: size(data_1d) /= goveq_richards_pres%mesh%ncells'
       write(iulog,*) 'size(data_1d)                    = ',size(data_1d)
       write(iulog,*) 'goveq_richards_pres%mesh%ncells  = ', goveq_richards_pres%mesh%ncells
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
    do local_id = 1, goveq_richards_pres%mesh%ncells
       press_p(local_id) = data_1d(local_id)
    enddo
    call VecRestoreArrayF90(soln_subvecs(1), press_p, ierr)

    ! Restore solution vectors for individual GEs
    call DMCompositeRestoreAccessArray(vsfm_soe%dm, vsfm_soe%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    call VecCopy(vsfm_soe%soln, vsfm_soe%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(vsfm_soe%soln, vsfm_soe%soln_prev_clm, ierr); CHKERRQ(ierr)

  end subroutine VSFMMPPRestart

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

end module MultiPhysicsProbVSFM

#endif
