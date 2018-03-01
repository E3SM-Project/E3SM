
module MultiPhysicsProbVSFM

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Object for variably saturated flow model
  !-----------------------------------------------------------------------

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                         , only : iulog
  use mpp_abortutils                     , only : endrun
  use mpp_shr_log_mod                    , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbBaseType           , only : multiphysicsprob_base_type
  use SystemOfEquationsVSFMType          , only : sysofeqns_vsfm_type
  use SystemOfEquationsBasePointerType   , only : sysofeqns_base_pointer_type
  use SystemOfEquationsBaseType          , only : sysofeqns_base_type
  use petscsys
  use petscvec
  use petscmat
  use petscts
  use petscsnes
  use petscdm
  use petscdmda

  implicit none
  private


  type, public, extends(multiphysicsprob_base_type) :: mpp_vsfm_type
   contains
     procedure, public :: Init            => VSFMMPPInit
     procedure, public :: Restart         => VSFMMPPRestart
     procedure, public :: AllocateAuxVars => VSFMMPPAllocateAuxVars
     procedure, public :: SetupProblem    => VSFMMPPSetupProblem
  end type mpp_vsfm_type

  type(mpp_vsfm_type), public, target :: vsfm_mpp
  public :: VSFMMPPSetSoils
  public :: VSFMMPPSetDensityType
  public :: VSFMMPPSetSoilPorosity
  public :: VSFMMPPSetSoilPermeability
  public :: VSFMMPPSetRelativePermeability
  public :: VSFMMPPSetSaturationFunction
  public :: VSFMMPPSetAuxVarConnRealValue
  public :: VSFMMPPSetAuxVarConnIntValue
  public :: VSFMMPPSetSourceSinkAuxVarRealValue
  public :: VSFMMPPSetSaturationFunctionAuxVarConn

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
    !
    class(sysofeqns_vsfm_type), pointer :: sysofeqns

    call MPPBaseInit(this)

    allocate(sysofeqns)
    call sysofeqns%Init()

    this%soe => sysofeqns

    allocate(this%soe_ptr)
    nullify(this%soe_ptr%ptr)

  end subroutine VSFMMPPInit

  !------------------------------------------------------------------------
  subroutine VSFMMPPInitialize(vsfm_mpp, begc, endc, ncols_ghost, &
       zc_col, filter_vsfmc, watsat, hksat, bsw, sucsat,          &
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
    integer, intent(in), pointer                      :: filter_vsfmc(:)
    PetscReal, intent(in), pointer                    :: watsat(:,:)
    PetscReal, intent(in), pointer                    :: hksat(:,:)
    PetscReal, intent(in), pointer                    :: bsw(:,:)
    PetscReal, intent(in), pointer                    :: sucsat(:,:)
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
    class(sysofeqns_base_type),pointer                :: base_soe
    Vec                                               :: variable_vec
    PetscReal, pointer                                :: vec_p(:)
    PetscErrorCode                                    :: ierr

    base_soe => vsfm_mpp%soe
    
    select type(base_soe)
    class is (sysofeqns_vsfm_type)
       vsfm_soe => base_soe
    end select

    ! Set initial coniditions
    call VSFMMPPSetSoils(vsfm_mpp, begc, endc, ncols_ghost, filter_vsfmc, &
         watsat, hksat, bsw, sucsat, residual_sat,          &
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
       watsat, hksat, bsw, sucsat, residual_sat, &
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
    integer, intent(in), pointer   :: filter_vsfmc(:) ! column filter for soil points
    PetscReal, intent(in), pointer :: watsat(:,:)
    PetscReal, intent(in), pointer :: hksat(:,:)
    PetscReal, intent(in), pointer :: bsw(:,:)
    PetscReal, intent(in), pointer :: sucsat(:,:)
    PetscReal, intent(in), pointer :: residual_sat(:,:)
    character(len=32), intent(in)  :: vsfm_satfunc_type
    PetscInt                       :: density_type

    select case(vsfm_mpp%id)
    case (MPP_VSFM_SNES_CLM)
       call VSFMMPPSetSoilsCLM(vsfm_mpp, begc, endc, ncols_ghost, filter_vsfmc, &
            watsat, hksat, bsw, sucsat, residual_sat, &
            vsfm_satfunc_type, density_type)
    case default
       write(iulog,*) 'VSFMMPPSetSoils: Unknown mpp_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine VSFMMPPSetSoils

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetSoilsCLM(vsfm_mpp, begc, endc, ncols_ghost, filter_vsfmc, &
       watsat, hksat, bsw, sucsat, residual_sat,                   &
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
    integer, intent(in), pointer                      :: filter_vsfmc(:)
    PetscReal, intent(in), pointer                    :: watsat(:,:)
    PetscReal, intent(in), pointer                    :: hksat(:,:)
    PetscReal, intent(in), pointer                    :: bsw(:,:)
    PetscReal, intent(in), pointer                    :: sucsat(:,:)
    PetscReal, intent(in), pointer                    :: residual_sat(:,:)
    character(len=32), intent(in)                     :: vsfm_satfunc_type
    PetscInt                                          :: density_type
    !
    ! !LOCAL VARIABLES:
    class (goveqn_richards_ode_pressure_type),pointer :: goveq_richards_ode_pres
    class(sysofeqns_vsfm_type),pointer                :: vsfm_soe
    class(sysofeqns_base_type),pointer                :: base_soe
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

    first_active_hydro_col_id = -1

    base_soe => vsfm_mpp%soe
    
    select type(base_soe)
    class is (sysofeqns_vsfm_type)
       vsfm_soe => base_soe
    end select
    
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
             ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

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
             ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

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
    class(sysofeqns_base_type),pointer                :: base_soe
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

    base_soe => vsfm_mpp%soe
    
    select type(base_soe)
    class is (sysofeqns_vsfm_type)
       vsfm_soe => base_soe
    end select
    
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
    call DMCompositeGetNumberDM(base_soe%solver%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(base_soe%solver%dm, dms, ierr)

    ! Allocate vectors for individual GEs
    allocate(soln_subvecs(nDM))

    ! Get solution vectors for individual GEs
    call DMCompositeGetAccessArray(base_soe%solver%dm, base_soe%solver%soln, nDM, &
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
    call DMCompositeRestoreAccessArray(base_soe%solver%dm, base_soe%solver%soln, nDM, &
                                       PETSC_NULL_INTEGER, soln_subvecs, ierr)

    call VecCopy(base_soe%solver%soln, base_soe%solver%soln_prev,     ierr); CHKERRQ(ierr)
    call VecCopy(base_soe%solver%soln, base_soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

    ! Free memory
    deallocate(dms)
    deallocate(soln_subvecs)

  end subroutine VSFMMPPSetICsSNESCLM

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
    class(sysofeqns_base_type),pointer                :: base_soe
    PetscInt                                          :: nDM
    DM, pointer                                       :: dms(:)
    Vec, pointer                                      :: soln_subvecs(:)
    PetscReal, pointer                                :: press_p(:)
    PetscInt                                          :: local_id

    base_soe => this%soe
    
    select type(base_soe)
    class is (sysofeqns_vsfm_type)
       vsfm_soe => base_soe
    class default
       write(iulog,*)'SoE is not of sysofeqns_vsfm_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select
    
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
    call DMCompositeGetNumberDM(base_soe%solver%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(base_soe%solver%dm, dms, ierr)

    ! Allocate vectors for individual GEs
    allocate(soln_subvecs(nDM))

    ! Get solution vectors for individual GEs
    call DMCompositeGetAccessArray(base_soe%solver%dm, base_soe%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    call VecGetArrayF90(soln_subvecs(1), press_p, ierr)
    do local_id = 1, goveq_richards_pres%mesh%ncells_local
       press_p(local_id) = data_1d(local_id)
    enddo
    call VecRestoreArrayF90(soln_subvecs(1), press_p, ierr)

    ! Restore solution vectors for individual GEs
    call DMCompositeRestoreAccessArray(base_soe%solver%dm, base_soe%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    call VecCopy(base_soe%solver%soln, base_soe%solver%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(base_soe%solver%soln, base_soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

    ! Free up memory
    deallocate(dms)
    deallocate(soln_subvecs)

  end subroutine VSFMMPPRestart

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
    use MultiPhysicsProbConstants           , only : COND_NULL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                 :: this
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_vsfm_type) , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: igoveqn
    PetscInt                             :: num_bc
    PetscInt                             :: num_ss
    PetscInt                             :: nconn_in
    PetscInt                             :: icond
    PetscInt                             :: iauxvar
    PetscInt                             :: iauxvar_beg, iauxvar_end
    PetscInt                             :: iauxvar_beg_bc, iauxvar_end_bc
    PetscInt                             :: iauxvar_beg_ss, iauxvar_end_ss
    PetscInt                             :: iauxvar_beg_conn_in, iauxvar_end_conn_in
    PetscInt                             :: cond_itype_to_exclude
    PetscInt                             :: count_bc, count_ss
    PetscInt                             :: offset_bc, offset_ss
    PetscInt                             :: cell_id
    PetscInt                             :: num_cells_bc
    PetscInt                             :: num_cells_ss
    PetscInt                   , pointer :: ncells_for_bc(:)
    PetscInt                   , pointer :: ncells_for_ss(:)
    PetscInt                   , pointer :: offsets_bc(:)
    PetscInt                   , pointer :: offsets_ss(:)
    PetscInt                   , pointer :: cell_id_dn_in_bc(:)
    PetscInt                   , pointer :: cell_id_dn_in_ss(:)
    PetscBool                  , pointer :: grid_cell_active(:)

    base_soe => this%soe

    select type(base_soe)
    class is(sysofeqns_vsfm_type)
       soe => base_soe
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
    soe%num_auxvars_conn_in = 0
    cur_goveq => soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          call cur_goveq%AllocateAuxVars()
          call cur_goveq%GetNumInternalConnections(nconn_in)
       end select

       cond_itype_to_exclude = COND_DIRICHLET_FRM_OTR_GOVEQ
       call cur_goveq%GetNCellsInCondsExcptCondItype(COND_BC, &
            cond_itype_to_exclude, num_bc, ncells_for_bc)

       cond_itype_to_exclude = COND_NULL
       call cur_goveq%GetNCellsInCondsExcptCondItype(COND_SS, &
            cond_itype_to_exclude, num_ss, ncells_for_ss)

       soe%num_auxvars_conn_in = soe%num_auxvars_conn_in + nconn_in
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
    allocate(soe%aux_vars_conn_in      (soe%num_auxvars_conn_in))

    allocate(soe%soe_auxvars_bc_offset (soe%num_auxvars_bc))
    allocate(soe%soe_auxvars_ss_offset (soe%num_auxvars_ss))
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
    iauxvar_beg_conn_in    = 0
    iauxvar_end_conn_in    = 0

    !
    ! Pass-2: Set values for auxvars of SoE
    !
    cur_goveq => soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       ! Boundary condition
       cond_itype_to_exclude = COND_DIRICHLET_FRM_OTR_GOVEQ
       call cur_goveq%GetNCellsInCondsExcptCondItype(COND_BC, &
            cond_itype_to_exclude, num_bc, ncells_for_bc)
       call cur_goveq%GetConnIDDnForCondsExcptCondItype(COND_BC, &
            cond_itype_to_exclude, num_cells_bc, cell_id_dn_in_bc)

       ! Source-sinks
       cond_itype_to_exclude = COND_NULL
       call cur_goveq%GetNCellsInCondsExcptCondItype(COND_SS, &
            cond_itype_to_exclude, num_ss, ncells_for_ss)
       call cur_goveq%GetConnIDDnForCondsExcptCondItype(COND_SS, &
            cond_itype_to_exclude, num_cells_ss, cell_id_dn_in_ss)

       call cur_goveq%GetMeshGridCellIsActive(grid_cell_active)

       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          call cur_goveq%GetNumInternalConnections(nconn_in)

       end select

       igoveqn = igoveqn + 1

       iauxvar_beg = iauxvar_end + 1
       iauxvar_end = iauxvar_end + cur_goveq%mesh%ncells_all

       ! AuxVars associated with grid cells
       cell_id     = 0
       do iauxvar = iauxvar_beg, iauxvar_end

          cell_id = cell_id + 1
          call soe%aux_vars_in(iauxvar)%Init()

          soe%aux_vars_in(iauxvar)%is_in     = PETSC_TRUE
          soe%aux_vars_in(iauxvar)%goveqn_id = igoveqn
          soe%aux_vars_in(iauxvar)%is_active = grid_cell_active(cell_id)
       enddo

       ! AuxVars associated with internal connections
       iauxvar_beg_conn_in = iauxvar_end_conn_in + 1
       iauxvar_end_conn_in = iauxvar_end_conn_in + nconn_in

       do iauxvar = iauxvar_beg_conn_in, iauxvar_end_conn_in
          call soe%aux_vars_conn_in(iauxvar)%Init()

          soe%aux_vars_conn_in(iauxvar)%is_in = PETSC_TRUE
          soe%aux_vars_conn_in(iauxvar)%goveqn_id = igoveqn

       end do

       allocate(offsets_bc(num_bc))
       allocate(offsets_ss(num_ss))

       do icond = 1, num_bc
          count_bc = count_bc + 1

          soe%soe_auxvars_bc_offset(count_bc) = offset_bc
          offsets_bc(icond)                   = offset_bc
          soe%soe_auxvars_bc_ncells(count_bc) = ncells_for_bc(icond)
          offset_bc                           = offset_bc + ncells_for_bc(icond)

          iauxvar_beg_bc = iauxvar_end_bc + 1
          iauxvar_end_bc = iauxvar_end_bc + ncells_for_bc(icond)

          do iauxvar = iauxvar_beg_bc, iauxvar_end_bc

             cell_id = cell_id_dn_in_bc(iauxvar)

             call soe%aux_vars_bc(iauxvar)%Init()

             soe%aux_vars_bc(iauxvar)%is_bc        = PETSC_TRUE
             soe%aux_vars_bc(iauxvar)%goveqn_id    = igoveqn
             soe%aux_vars_bc(iauxvar)%condition_id = icond
             soe%aux_vars_bc(iauxvar)%is_active    = grid_cell_active(cell_id)
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

             cell_id = cell_id_dn_in_ss(iauxvar)

             call soe%aux_vars_ss(iauxvar)%Init()

             soe%aux_vars_ss(iauxvar)%is_ss        = PETSC_TRUE
             soe%aux_vars_ss(iauxvar)%goveqn_id    = igoveqn
             soe%aux_vars_ss(iauxvar)%condition_id = icond
             soe%aux_vars_ss(iauxvar)%is_active    = grid_cell_active(cell_id)
          enddo
       enddo

       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          call cur_goveq%SetSOEAuxVarOffsets(num_bc, offsets_bc, num_ss, offsets_ss)
       end select

       if (num_bc > 0) then
          deallocate(ncells_for_bc    )
          deallocate(offsets_bc       )
          deallocate(cell_id_dn_in_bc )
       end if
       if (num_ss > 0) then
          deallocate(ncells_for_ss    )
          deallocate(offsets_ss       )
          deallocate(cell_id_dn_in_ss )
       end if

       deallocate(grid_cell_active)

       cur_goveq => cur_goveq%next
    enddo

    do iauxvar = 1, soe%num_auxvars_conn_in
       call soe%aux_vars_conn_in(iauxvar)%Init()
    end do

  end subroutine VSFMMPPAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetupProblem(this, dyn_linesearch)
    !
    ! !DESCRIPTION:
    ! Sets up the PETSc related data structure
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : SOE_RE_ODE
    use MultiPhysicsProbBaseType  , only : MPPSetupProblem
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type) :: this
    PetscBool, optional  :: dyn_linesearch

    call MPPSetupProblem(this)

    this%soe%itype = SOE_RE_ODE

    if (present(dyn_linesearch)) then
       call this%soe%solver%SetUseDynLineSearch(dyn_linesearch)
    endif

  end subroutine VSFMMPPSetupProblem

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetAuxVarConnRealValue(this, igoveqn, auxvar_type, &
       var_type, var_value)
    !
    ! !DESCRIPTION:
    ! Set real value of conn-auxvars
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type) , intent(inout) :: this
    PetscInt             , intent(in)    :: igoveqn
    PetscInt             , intent(in)    :: auxvar_type
    PetscInt             , intent(in)    :: var_type
    PetscReal, pointer   , intent(in)    :: var_value(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                             :: ii
    class(goveqn_base_type) , pointer    :: cur_goveq

    if (igoveqn > this%soe%ngoveqns) then
       write(iulog,*) 'Attempting to access governing equation ' // &
            'that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq => this%soe%goveqns
    do ii = 1, igoveqn-1
       cur_goveq => cur_goveq%next
    enddo

    select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          call cur_goveq%SetAuxVarConnRealValue(auxvar_type, var_type, var_value)

       class default
          write(iulog,*) 'VSFMMPPSetAuxVarConnRealValue only supports ' // &
               'goveqn_richards_ode_pressure_type'
          call endrun(msg=errMsg(__FILE__,__LINE__))
    end select

  end subroutine VSFMMPPSetAuxVarConnRealValue

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetAuxVarConnIntValue(this, igoveqn, auxvar_type, &
       var_type, var_value)
    !
    ! !DESCRIPTION:
    ! Set integer value of conn-auxvars
    !
    ! !USES:
    use MultiPhysicsProbConstants     , only : AUXVAR_BC
    use MultiPhysicsProbConstants     , only : AUXVAR_INTERNAL
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type) , intent(inout) :: this
    PetscInt             , intent(in)    :: igoveqn
    PetscInt             , intent(in)    :: auxvar_type
    PetscInt             , intent(in)    :: var_type
    PetscInt, pointer    , intent(in)    :: var_value(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                             :: ii
    class(goveqn_base_type) , pointer    :: cur_goveq

    if (igoveqn > this%soe%ngoveqns) then
       write(iulog,*) 'Attempting to access governing equation ' // &
            'that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq => this%soe%goveqns
    do ii = 1, igoveqn-1
       cur_goveq => cur_goveq%next
    enddo

    select type(cur_goveq)
    class is (goveqn_richards_ode_pressure_type)
       call cur_goveq%SetAuxVarConnIntValue(auxvar_type, var_type, var_value)

    class default
       write(iulog,*) 'VSFMMPPSetAuxVarConnRealValue only supports ' // &
            'goveqn_richards_ode_pressure_type'
       call endrun(msg=errMsg(__FILE__,__LINE__))
    end select


  end subroutine VSFMMPPSetAuxVarConnIntValue

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetSoilPermeability(this, igoveqn, perm_x, perm_y, perm_z)
    !
    ! !DESCRIPTION:
    ! Set soil permeability values
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use RichardsODEPressureAuxType    , only : rich_ode_pres_auxvar_type
    use ConditionType                 , only : condition_type
    use ConnectionSetType             , only : connection_set_type
    use MultiPhysicsProbConstants     , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                      , intent(inout)       :: this
    PetscInt                                  , intent(in)          :: igoveqn
    PetscReal                                 , pointer, intent(in) :: perm_x(:)
    PetscReal                                 , pointer, intent(in) :: perm_y(:)
    PetscReal                                 , pointer, intent(in) :: perm_z(:)
    !
    ! !LOCAL VARIABLES:
    class (goveqn_richards_ode_pressure_type) , pointer             :: goveq_richards_ode_pres
    class(sysofeqns_vsfm_type)                , pointer             :: vsfm_soe
    class(goveqn_base_type)                   , pointer             :: cur_goveq
    class(sysofeqns_base_type)                , pointer             :: base_soe
    type (rich_ode_pres_auxvar_type)          , pointer             :: ode_aux_vars_in(:)
    type (rich_ode_pres_auxvar_type)          , pointer             :: ode_aux_vars_bc(:)
    type (rich_ode_pres_auxvar_type)          , pointer             :: ode_aux_vars_ss(:)
    PetscInt                                                        :: icell
    PetscInt                                                        :: ii
    type(condition_type)                      , pointer             :: cur_cond
    type(connection_set_type)                 , pointer             :: cur_conn_set
    PetscInt                                                        :: ghosted_id
    PetscInt                                                        :: sum_conn
    PetscInt                                                        :: iconn

    base_soe => this%soe
    
    select type(base_soe)
    class is (sysofeqns_vsfm_type)
       vsfm_soe => base_soe
    end select

    if (igoveqn > this%soe%ngoveqns) then
       write(iulog,*) 'Attempting to add condition for governing equation ' // &
            'that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq => this%soe%goveqns
    do ii = 1, igoveqn-1
       cur_goveq => cur_goveq%next
    enddo

    select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          goveq_richards_ode_pres => cur_goveq
       class default
          write(iulog,*)'Only goveqn_richards_ode_pressure_type supported'
          call endrun(msg=errMsg(__FILE__,__LINE__))
    end select

    if (goveq_richards_ode_pres%mesh%ncells_local /= size(perm_x)) then
       write(iulog,*) 'No. of values for soil permeability is not equal to no. of grid cells.'
       write(iulog,*) 'No. of soil porosity values = ',size(perm_x)
       write(iulog,*) 'No. of grid cells           = ',goveq_richards_ode_pres%mesh%ncells_local
    endif

    ode_aux_vars_in => goveq_richards_ode_pres%aux_vars_in

    ! Set soil properties for internal auxvars
    do icell = 1, size(perm_x)
       ode_aux_vars_in(icell)%perm(1) = perm_x(icell)
       ode_aux_vars_in(icell)%perm(2) = perm_y(icell)
       ode_aux_vars_in(icell)%perm(3) = perm_z(icell)
    enddo

    ! Set soil properties for boundary-condition auxvars
    sum_conn = 0
    ode_aux_vars_bc => goveq_richards_ode_pres%aux_vars_bc
    cur_cond        => goveq_richards_ode_pres%boundary_conditions%first

    do
       if (.not.associated(cur_cond)) exit
       if (cur_cond%itype /= COND_DIRICHLET_FRM_OTR_GOVEQ) then
          cur_conn_set => cur_cond%conn_set

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

             ode_aux_vars_bc(sum_conn)%perm(:) = ode_aux_vars_in(ghosted_id)%perm(:)

          enddo
       else
          cur_conn_set => cur_cond%conn_set

          ghosted_id = 1
          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             ode_aux_vars_bc(sum_conn)%perm(:) = ode_aux_vars_in(ghosted_id)%perm(:)
          enddo
          
       endif
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
          ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

          ode_aux_vars_ss(sum_conn)%perm (:) = ode_aux_vars_in(ghosted_id)%perm(:)

       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine VSFMMPPSetSoilPermeability

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetRelativePermeability(this, igoveqn, relperm_type, param_1, param_2)
    !
    ! !DESCRIPTION:
    ! Set soil permeability values
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use RichardsODEPressureAuxType    , only : rich_ode_pres_auxvar_type
    use ConditionType                 , only : condition_type
    use ConnectionSetType             , only : connection_set_type
    use MultiPhysicsProbConstants     , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use SaturationFunction            , only : RELPERM_FUNC_MUALEM
    use SaturationFunction            , only : RELPERM_FUNC_WEIBULL
    use SaturationFunction            , only : RELPERM_FUNC_CAMPBELL
    use SaturationFunction            , only : SatFunc_Set_Weibull_RelPerm
    use SaturationFunction            , only : SatFunc_Set_Campbell_RelPerm
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                      , intent(inout)       :: this
    PetscInt                                  , intent(in)          :: igoveqn
    PetscInt                                  , pointer, intent(in) :: relperm_type(:)
    PetscReal                                 , pointer, intent(in) :: param_1(:)
    PetscReal                                 , pointer, intent(in) :: param_2(:)
    !
    ! !LOCAL VARIABLES:
    class (goveqn_richards_ode_pressure_type) , pointer             :: goveq_richards_ode_pres
    class(goveqn_base_type)                   , pointer             :: cur_goveq
    type (rich_ode_pres_auxvar_type)          , pointer             :: ode_aux_vars_in(:)
    type (rich_ode_pres_auxvar_type)          , pointer             :: ode_aux_vars_bc(:)
    type (rich_ode_pres_auxvar_type)          , pointer             :: ode_aux_vars_ss(:)
    PetscInt                                                        :: icell
    PetscInt                                                        :: igov
    type(condition_type)                      , pointer             :: cur_cond
    type(connection_set_type)                 , pointer             :: cur_conn_set
    PetscInt                                                        :: ghosted_id
    PetscInt                                                        :: sum_conn
    PetscInt                                                        :: iconn

    if (igoveqn > this%soe%ngoveqns) then
       write(iulog,*) 'Attempting to add condition for governing equation ' // &
            'that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq => this%soe%goveqns
    do igov = 1, igoveqn-1
       cur_goveq => cur_goveq%next
    enddo

    select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          goveq_richards_ode_pres => cur_goveq
       class default
          write(iulog,*)'Only goveqn_richards_ode_pressure_type supported'
          call endrun(msg=errMsg(__FILE__,__LINE__))
    end select

    if (goveq_richards_ode_pres%mesh%ncells_local /= size(relperm_type)) then
       write(iulog,*) 'No. of values for relative permeability is not equal to no. of grid cells.'
       write(iulog,*) 'No. of rel. perm. values = ',size(relperm_type)
       write(iulog,*) 'No. of grid cells        = ',goveq_richards_ode_pres%mesh%ncells_local
    endif

    ode_aux_vars_in => goveq_richards_ode_pres%aux_vars_in

    ! Set soil properties for internal auxvars
    do icell = 1, size(param_1)

       select case (relperm_type(icell))
       case (RELPERM_FUNC_MUALEM)

       case (RELPERM_FUNC_WEIBULL)
          call SatFunc_Set_Weibull_RelPerm(      &
               ode_aux_vars_in(icell)%satParams, &
               param_1(icell),                   &
               param_2(icell))

       case (RELPERM_FUNC_CAMPBELL)
          call SatFunc_Set_Campbell_RelPerm(              &
               ode_aux_vars_in(icell)%satParams, &
               param_1(icell),                   &
               param_2(icell))

       case default
          call endrun(msg='ERROR:: Unknown relperm_type ' // &
               errMsg(__FILE__, __LINE__))
       end select
    end do

    ! Set soil properties for boundary-condition auxvars
    sum_conn = 0
    ode_aux_vars_bc => goveq_richards_ode_pres%aux_vars_bc
    cur_cond        => goveq_richards_ode_pres%boundary_conditions%first

    do
       if (.not.associated(cur_cond)) exit
       if (cur_cond%itype /= COND_DIRICHLET_FRM_OTR_GOVEQ) then
          cur_conn_set => cur_cond%conn_set

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

             call ode_aux_vars_bc(sum_conn)%satParams%Copy(ode_aux_vars_in(ghosted_id)%satParams)

          enddo
       else
          cur_conn_set => cur_cond%conn_set

          ghosted_id = 1
          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             call ode_aux_vars_bc(sum_conn)%satParams%Copy(ode_aux_vars_in(ghosted_id)%satParams)
          enddo

       endif
       cur_cond => cur_cond%next
    enddo

    ! Set satParams for source-sink auxvars
    sum_conn = 0

    ode_aux_vars_ss => goveq_richards_ode_pres%aux_vars_ss
    cur_cond        => goveq_richards_ode_pres%source_sinks%first

    do
       if (.not.associated(cur_cond)) exit
       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1
          ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()
          call ode_aux_vars_ss(sum_conn)%satParams%Copy(ode_aux_vars_in(ghosted_id)%satParams)

       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine VSFMMPPSetRelativePermeability

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetDensityType(this, igoveqn, density_type)
    !
    ! !DESCRIPTION:
    ! Set density type
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use RichardsODEPressureAuxType    , only : rich_ode_pres_auxvar_type
    use ConditionType                 , only : condition_type
    use ConnectionSetType             , only : connection_set_type
    use PorosityFunctionMod           , only : PorosityFunctionSetConstantModel
    use MultiPhysicsProbConstants     , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                      , intent(inout)       :: this
    PetscInt                                  , intent(in)          :: igoveqn
    PetscInt                                  , intent(in)          :: density_type
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type)                   , pointer             :: cur_goveq
    PetscInt                                                        :: ii

    if (igoveqn > this%soe%ngoveqns) then
       write(iulog,*) 'Attempting to add condition for governing equation ' // &
            'that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq => this%soe%goveqns
    do ii = 1, igoveqn-1
       cur_goveq => cur_goveq%next
    enddo

    select type(cur_goveq)
    class is (goveqn_richards_ode_pressure_type)
       call cur_goveq%SetDensityType(density_type)
    class default
       write(iulog,*)'Only goveqn_richards_ode_pressure_type supported'
       call endrun(msg=errMsg(__FILE__,__LINE__))
    end select

  end subroutine VSFMMPPSetDensityType

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetSoilPorosity(this, igoveqn, por)
    !
    ! !DESCRIPTION:
    ! Set soil porosity value
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use RichardsODEPressureAuxType    , only : rich_ode_pres_auxvar_type
    use ConditionType                 , only : condition_type
    use ConnectionSetType             , only : connection_set_type
    use PorosityFunctionMod           , only : PorosityFunctionSetConstantModel
    use MultiPhysicsProbConstants     , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                      , intent(inout)       :: this
    PetscInt                                  , intent(in)          :: igoveqn
    PetscReal                                 , pointer, intent(in) :: por(:)
    !
    ! !LOCAL VARIABLES:
    class (goveqn_richards_ode_pressure_type) , pointer             :: goveq_richards_ode_pres
    class(sysofeqns_vsfm_type)                , pointer             :: vsfm_soe
    class(sysofeqns_base_type)                , pointer             :: base_soe
    class(goveqn_base_type)                   , pointer             :: cur_goveq
    type (rich_ode_pres_auxvar_type)          , pointer             :: ode_aux_vars_in(:)
    type (rich_ode_pres_auxvar_type)          , pointer             :: ode_aux_vars_bc(:)
    type (rich_ode_pres_auxvar_type)          , pointer             :: ode_aux_vars_ss(:)
    PetscInt                                                        :: icell
    PetscInt                                                        :: ii
    type(condition_type)                      , pointer             :: cur_cond
    type(connection_set_type)                 , pointer             :: cur_conn_set
    PetscInt                                                        :: ghosted_id
    PetscInt                                                        :: sum_conn
    PetscInt                                                        :: iconn

    base_soe => this%soe
    
    select type(base_soe)
    class is (sysofeqns_vsfm_type)
       vsfm_soe => base_soe
    end select

    if (igoveqn > this%soe%ngoveqns) then
       write(iulog,*) 'Attempting to add condition for governing equation ' // &
            'that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq => this%soe%goveqns
    do ii = 1, igoveqn-1
       cur_goveq => cur_goveq%next
    enddo

    select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          goveq_richards_ode_pres => cur_goveq
       class default
          write(iulog,*)'Only goveqn_richards_ode_pressure_type supported'
          call endrun(msg=errMsg(__FILE__,__LINE__))
    end select

    if (goveq_richards_ode_pres%mesh%ncells_local /= size(por)) then
       write(iulog,*) 'No. of values for soil porosity is not equal to no. of grid cells.'
       write(iulog,*) 'No. of soil porosity values = ',size(por)
       write(iulog,*) 'No. of grid cells           = ',goveq_richards_ode_pres%mesh%ncells_local
    endif

    ode_aux_vars_in => goveq_richards_ode_pres%aux_vars_in

    ! Set soil properties for internal auxvars
    do icell = 1, size(por)
       ode_aux_vars_in(icell)%por = por(icell)
       call PorosityFunctionSetConstantModel(ode_aux_vars_in(icell)%porParams, &
            ode_aux_vars_in(icell)%por)
    enddo

    ! Set soil properties for boundary-condition auxvars
    sum_conn = 0
    ode_aux_vars_bc => goveq_richards_ode_pres%aux_vars_bc
    cur_cond        => goveq_richards_ode_pres%boundary_conditions%first

    do
       if (.not.associated(cur_cond)) exit
       if (cur_cond%itype /= COND_DIRICHLET_FRM_OTR_GOVEQ) then

          cur_conn_set => cur_cond%conn_set

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

             ode_aux_vars_bc(sum_conn)%por       = ode_aux_vars_in(ghosted_id)%por
             ode_aux_vars_bc(sum_conn)%porParams = ode_aux_vars_in(ghosted_id)%porParams

          enddo
       else
          cur_conn_set => cur_cond%conn_set

          ghosted_id = 1
          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1

             ode_aux_vars_bc(sum_conn)%por       = ode_aux_vars_in(ghosted_id)%por
             ode_aux_vars_bc(sum_conn)%porParams = ode_aux_vars_in(ghosted_id)%porParams

          enddo
       end if
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
          ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

          ode_aux_vars_ss(sum_conn)%por       = ode_aux_vars_in(ghosted_id)%por
          ode_aux_vars_ss(sum_conn)%porParams = ode_aux_vars_in(ghosted_id)%porParams

       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine VSFMMPPSetSoilPorosity

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetSaturationFunction(this, igoveqn, vsfm_satfunc_type, &
       alpha, lambda, sat_res)
    !
    ! !DESCRIPTION:
    ! Set saturation function
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use RichardsODEPressureAuxType    , only : rich_ode_pres_auxvar_type
    use ConditionType                 , only : condition_type
    use ConnectionSetType             , only : connection_set_type
    use SaturationFunction            , only : SatFunc_Set_BC
    use SaturationFunction            , only : SatFunc_Set_SBC_bz2
    use SaturationFunction            , only : SatFunc_Set_SBC_bz3
    use SaturationFunction            , only : SatFunc_Set_VG
    use SaturationFunction            , only : SatFunc_Set_FETCH2
    use SaturationFunction            , only : SatFunc_Set_Chuang
    use SaturationFunction            , only : SAT_FUNC_VAN_GENUCHTEN
    use SaturationFunction            , only : SAT_FUNC_BROOKS_COREY
    use SaturationFunction            , only : SAT_FUNC_SMOOTHED_BROOKS_COREY_BZ2
    use SaturationFunction            , only : SAT_FUNC_SMOOTHED_BROOKS_COREY_BZ3
    use SaturationFunction            , only : SAT_FUNC_FETCH2
    use SaturationFunction            , only : SAT_FUNC_CHUANG
    use MultiPhysicsProbConstants     , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                      , intent(inout)       :: this
    PetscInt                                  , intent(in)          :: igoveqn
    PetscInt                                  , pointer, intent(in) :: vsfm_satfunc_type(:)
    PetscReal                                 , pointer, intent(in) :: alpha(:)
    PetscReal                                 , pointer, intent(in) :: lambda(:)
    PetscReal                                 , pointer, intent(in) :: sat_res(:)
    !
    ! !LOCAL VARIABLES:
    class (goveqn_richards_ode_pressure_type) , pointer             :: goveq_richards_ode_pres
    class(sysofeqns_vsfm_type)                , pointer             :: vsfm_soe
    class(sysofeqns_base_type)                , pointer             :: base_soe
    class(goveqn_base_type)                   , pointer             :: cur_goveq
    type (rich_ode_pres_auxvar_type)          , pointer             :: ode_aux_vars_in(:)
    type (rich_ode_pres_auxvar_type)          , pointer             :: ode_aux_vars_bc(:)
    type (rich_ode_pres_auxvar_type)          , pointer             :: ode_aux_vars_ss(:)
    PetscInt                                                        :: icell
    PetscInt                                                        :: ii
    type(condition_type)                      , pointer             :: cur_cond
    type(connection_set_type)                 , pointer             :: cur_conn_set
    PetscInt                                                        :: ghosted_id
    PetscInt                                                        :: sum_conn
    PetscInt                                                        :: iconn

    base_soe => this%soe
    
    select type(base_soe)
    class is (sysofeqns_vsfm_type)
       vsfm_soe => base_soe
    end select

    if (igoveqn > this%soe%ngoveqns) then
       write(iulog,*) 'Attempting to add condition for governing equation ' // &
            'that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq => this%soe%goveqns
    do ii = 1, igoveqn-1
       cur_goveq => cur_goveq%next
    enddo

    select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          goveq_richards_ode_pres => cur_goveq
       class default
          write(iulog,*)'Only goveqn_richards_ode_pressure_type supported'
          call endrun(msg=errMsg(__FILE__,__LINE__))
    end select

    if (goveq_richards_ode_pres%mesh%ncells_local /= size(alpha)) then
       write(iulog,*) 'No. of values for saturation function is not equal to no. of grid cells.'
       write(iulog,*) 'No. of soil sat. func. values = ',size(alpha)
       write(iulog,*) 'No. of grid cells             = ',goveq_richards_ode_pres%mesh%ncells_local
    endif

    ode_aux_vars_in => goveq_richards_ode_pres%aux_vars_in

    ! Set soil properties for internal auxvars

    do icell = 1, size(alpha)

       select case (vsfm_satfunc_type(icell))
       case (SAT_FUNC_BROOKS_COREY)
          call SatFunc_Set_BC(                   &
               ode_aux_vars_in(icell)%satParams, &
               sat_res(icell),                   &
               alpha(icell),                     &
               lambda(icell))

       case (SAT_FUNC_SMOOTHED_BROOKS_COREY_BZ2)
          call SatFunc_Set_SBC_bz2(              &
               ode_aux_vars_in(icell)%satParams, &
               sat_res(icell),                   &
               alpha(icell),                     &
               lambda(icell),                    &
               -0.9d0/alpha(icell))

       case (SAT_FUNC_SMOOTHED_BROOKS_COREY_BZ3)
          call SatFunc_Set_SBC_bz3(              &
               ode_aux_vars_in(icell)%satParams, &
               sat_res(icell),                   &
               alpha(icell),                     &
               lambda(icell),                    &
               -0.9d0/alpha(icell))

       case (SAT_FUNC_VAN_GENUCHTEN)
          call SatFunc_Set_VG(                   &
               ode_aux_vars_in(icell)%satParams, &
               sat_res(icell),                   &
               alpha(icell),                     &
               lambda(icell))

       case (SAT_FUNC_FETCH2)
          call SatFunc_Set_FETCH2(               &
               ode_aux_vars_in(icell)%satParams, &
               alpha(icell),                     &
               lambda(icell))

       case (SAT_FUNC_CHUANG)
          call SatFunc_Set_Chuang(               &
               ode_aux_vars_in(icell)%satParams, &
               alpha(icell),                     &
               lambda(icell))

       case default
          call endrun(msg='ERROR:: Unknown vsfm_satfunc_type ' // &
               errMsg(__FILE__, __LINE__))
       end select
    end do

    ! Set soil properties for boundary-condition auxvars
    sum_conn = 0
    ode_aux_vars_bc => goveq_richards_ode_pres%aux_vars_bc
    cur_cond        => goveq_richards_ode_pres%boundary_conditions%first

    do
       if (.not.associated(cur_cond)) exit
       if (cur_cond%itype /= COND_DIRICHLET_FRM_OTR_GOVEQ) then
          cur_conn_set => cur_cond%conn_set

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

             call ode_aux_vars_bc(sum_conn)%satParams%Copy(ode_aux_vars_in(ghosted_id)%satParams)

          enddo
       else
          cur_conn_set => cur_cond%conn_set

          ghosted_id = 1
          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1

             call ode_aux_vars_bc(sum_conn)%satParams%Copy(ode_aux_vars_in(ghosted_id)%satParams)

          enddo
       endif
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
          ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

          call ode_aux_vars_ss(sum_conn)%satParams%Copy(ode_aux_vars_in(ghosted_id)%satParams)

       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine VSFMMPPSetSaturationFunction

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetSourceSinkAuxVarRealValue(this, igoveqn, &
       var_type, var_value)
    !
    ! !DESCRIPTION:
    ! Set real values for source-sink auxvars
    !
    ! !USES:
    use MultiPhysicsProbConstants     , only : AUXVAR_BC
    use MultiPhysicsProbConstants     , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants     , only : VAR_POT_MASS_SINK_PRESSURE
    use MultiPhysicsProbConstants     , only : VAR_POT_MASS_SINK_EXPONENT
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use RichardsODEPressureAuxType    , only : rich_ode_pres_auxvar_type
    use ConditionType                 , only : condition_type
    use ConnectionSetType             , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                      , intent(inout)          :: this
    PetscInt                                  , intent(in)             :: igoveqn
    PetscInt                                  , intent(in)             :: var_type
    PetscReal                                 , pointer   , intent(in) :: var_value(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                                           :: ii
    class(goveqn_base_type)                   , pointer                :: cur_goveq
    class (goveqn_richards_ode_pressure_type) , pointer                :: goveq_richards_ode_pres
    type (rich_ode_pres_auxvar_type)          , pointer                :: ode_aux_vars_ss(:)
    type(condition_type)                      , pointer                :: cur_cond
    type(connection_set_type)                 , pointer                :: cur_conn_set
    PetscInt                                                           :: sum_conn
    PetscInt                                                           :: iconn

    if (igoveqn > this%soe%ngoveqns) then
       write(iulog,*) 'Attempting to access governing equation ' // &
            'that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq => this%soe%goveqns
    do ii = 1, igoveqn-1
       cur_goveq => cur_goveq%next
    enddo

    select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          goveq_richards_ode_pres => cur_goveq

       class default
          write(iulog,*) 'VSFMMPPSetSourceSinkAuxVarRealValue only supports ' // &
               'goveqn_richards_ode_pressure_type'
          call endrun(msg=errMsg(__FILE__,__LINE__))
    end select

    ! Set soil properties for source-sink auxvars
    sum_conn = 0

    ode_aux_vars_ss => goveq_richards_ode_pres%aux_vars_ss
    cur_cond        => goveq_richards_ode_pres%source_sinks%first

    do
       if (.not.associated(cur_cond)) exit
       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1

          select case(var_type)
              case (VAR_POT_MASS_SINK_PRESSURE)
                 ode_aux_vars_ss(sum_conn)%pot_mass_sink_pressure = var_value(sum_conn)
              case (VAR_POT_MASS_SINK_EXPONENT)
                 ode_aux_vars_ss(sum_conn)%pot_mass_sink_exponent = var_value(sum_conn)
              case default
                 write(iulog,*) 'VSFMMPPSetSourceSinkAuxVarRealValue: Unknown var_type'
                 call endrun(msg=errMsg(__FILE__,__LINE__))
              end select

       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine VSFMMPPSetSourceSinkAuxVarRealValue

  !------------------------------------------------------------------------
  subroutine VSFMMPPSetSaturationFunctionAuxVarConn(this, igoveqn, &
       auxvar_conn_type, set_upwind_auxvar, satfunc_itype, param_1, param_2, param_3)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use GoverningEquationBaseType      , only : goveqn_base_type
    use GoveqnRichardsODEPressureType  , only : goveqn_richards_ode_pressure_type
    use RichardsODEPressureConnAuxType , only : rich_ode_pres_conn_auxvar_type
    use ConditionType                  , only : condition_type
    use ConnectionSetType              , only : connection_set_type
    use SaturationFunction             , only : SatFunc_Set_Campbell_RelPerm
    use SaturationFunction             , only : SatFunc_Set_Weibull_RelPerm
    use SaturationFunction             , only : SatFunc_Set_Chuang
    use SaturationFunction             , only : SatFunc_Set_VG
    use SaturationFunction             , only : RELPERM_FUNC_CAMPBELL
    use SaturationFunction             , only : RELPERM_FUNC_WEIBULL
    use SaturationFunction             , only : RELPERM_FUNC_MUALEM
    use SaturationFunction             , only : SAT_FUNC_CHUANG
    use SaturationFunction             , only : SAT_FUNC_VAN_GENUCHTEN
    use MultiPhysicsProbConstants      , only : AUXVAR_CONN_BC
    use MultiPhysicsProbConstants      , only : AUXVAR_CONN_INTERNAL
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    class(mpp_vsfm_type)                      , intent(inout)       :: this
    PetscInt                                  , intent(in)          :: igoveqn
    PetscInt                                  , intent(in)          :: auxvar_conn_type
    PetscBool                                 , intent(in)          :: set_upwind_auxvar(:)
    PetscInt                                  , pointer, intent(in) :: satfunc_itype(:)
    PetscReal                                 , pointer, intent(in) :: param_1(:)
    PetscReal                                 , pointer, intent(in) :: param_2(:)
    PetscReal                      , optional , pointer, intent(in) :: param_3(:)
    !
    ! !LOCAL VARIABLES:
    class (goveqn_richards_ode_pressure_type) , pointer             :: goveq_richards_ode_pres
    class(sysofeqns_vsfm_type)                , pointer             :: vsfm_soe
    class(sysofeqns_base_type)                , pointer             :: base_soe
    class(goveqn_base_type)                   , pointer             :: cur_goveq
    type (rich_ode_pres_conn_auxvar_type)     , pointer             :: conn_aux_vars(:)
    type(condition_type)                      , pointer             :: cur_cond
    type(connection_set_type)                 , pointer             :: cur_conn_set
    PetscInt                                                        :: ii
    PetscInt                                                        :: ghosted_id
    PetscInt                                                        :: sum_conn
    PetscInt                                                        :: iconn

    base_soe => this%soe
    
    select type(base_soe)
    class is (sysofeqns_vsfm_type)
       vsfm_soe => base_soe
    end select

    if (igoveqn > this%soe%ngoveqns) then
       write(iulog,*) 'Attempting to add condition for governing equation ' // &
            'that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq => this%soe%goveqns
    do ii = 1, igoveqn-1
       cur_goveq => cur_goveq%next
    enddo

    select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          goveq_richards_ode_pres => cur_goveq
       class default
          write(iulog,*)'Only goveqn_richards_ode_pressure_type supported'
          call endrun(msg=errMsg(__FILE__,__LINE__))
    end select


    select case(auxvar_conn_type)
    case (AUXVAR_CONN_INTERNAL)
       conn_aux_vars => goveq_richards_ode_pres%aux_vars_conn_in
       cur_conn_set  => goveq_richards_ode_pres%mesh%intrn_conn_set_list%first
       sum_conn = 0
       do
          if (.not.associated(cur_conn_set)) exit

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1

             if (sum_conn > size(param_1)) then
                write(iulog,*) 'No. of values for saturation function is not equal to no. connections.'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end if

             select case(satfunc_itype(sum_conn))
             case (0)
                ! Do nothing

             case (RELPERM_FUNC_MUALEM)
                if (set_upwind_auxvar(sum_conn)) then
                   conn_aux_vars(sum_conn)%satParams_up%relperm_func_type = RELPERM_FUNC_MUALEM
                   conn_aux_vars(sum_conn)%satParams_up%alpha   = param_1(sum_conn)
                   conn_aux_vars(sum_conn)%satParams_up%vg_m    = param_2(sum_conn)
                   conn_aux_vars(sum_conn)%satParams_up%sat_res = param_3(sum_conn)
                   conn_aux_vars(sum_conn)%satParams_up%vg_n    = 1.d0/(1.d0 - param_2(sum_conn))
                else
                   conn_aux_vars(sum_conn)%satParams_dn%relperm_func_type = RELPERM_FUNC_MUALEM
                   conn_aux_vars(sum_conn)%satParams_dn%alpha             = param_1(sum_conn)
                   conn_aux_vars(sum_conn)%satParams_dn%vg_m              = param_2(sum_conn)
                   conn_aux_vars(sum_conn)%satParams_dn%sat_res           = param_3(sum_conn)
                   conn_aux_vars(sum_conn)%satParams_dn%vg_n              = 1.d0/(1.d0 - param_2(sum_conn))
                end if
             case (RELPERM_FUNC_CAMPBELL)

                if (set_upwind_auxvar(sum_conn)) then
                   call SatFunc_Set_Campbell_RelPerm(conn_aux_vars(sum_conn)%satParams_up, &
                        param_1(sum_conn), param_2(sum_conn))
                else
                   call SatFunc_Set_Campbell_RelPerm(conn_aux_vars(sum_conn)%satParams_dn, &
                        param_1(sum_conn), param_2(sum_conn))
                endif

             case (RELPERM_FUNC_WEIBULL)

                if (set_upwind_auxvar(sum_conn)) then
                   call SatFunc_Set_Weibull_RelPerm(conn_aux_vars(sum_conn)%satParams_up, &
                        param_1(sum_conn), param_2(sum_conn))
                else
                   call SatFunc_Set_Weibull_RelPerm(conn_aux_vars(sum_conn)%satParams_dn, &
                        param_1(sum_conn), param_2(sum_conn))
                endif

             case (SAT_FUNC_CHUANG)

                if (set_upwind_auxvar(sum_conn)) then
                   call SatFunc_Set_Chuang(conn_aux_vars(sum_conn)%satParams_up, &
                        param_1(sum_conn), param_2(sum_conn))
                else
                   call SatFunc_Set_Chuang(conn_aux_vars(sum_conn)%satParams_dn, &
                        param_1(sum_conn), param_2(sum_conn))
                endif

             case (SAT_FUNC_VAN_GENUCHTEN)
                if (set_upwind_auxvar(sum_conn)) then
                   call SatFunc_Set_VG(conn_aux_vars(sum_conn)%satParams_up, &
                        param_3(sum_conn), param_2(sum_conn), param_1(sum_conn))
                else
                   call SatFunc_Set_VG(conn_aux_vars(sum_conn)%satParams_dn, &
                        param_3(sum_conn), param_2(sum_conn), param_1(sum_conn))
                endif

             case default
                write(iulog,*)'Unknown satfunc_itype_variable ', satfunc_itype(sum_conn)
                call endrun(msg=errMsg(__FILE__,__LINE__))
             end select

          end do

          cur_conn_set => cur_conn_set%next

       end do

    case(AUXVAR_CONN_BC)

       ! Set soil properties for boundary-condition auxvars
       sum_conn = 0
       conn_aux_vars => goveq_richards_ode_pres%aux_vars_conn_bc
       cur_cond      => goveq_richards_ode_pres%boundary_conditions%first

       do
          if (.not.associated(cur_cond)) exit
          cur_conn_set => cur_cond%conn_set

          do iconn = 1, cur_conn_set%num_connections
             sum_conn   = sum_conn + 1
             ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

             if (sum_conn> size(param_1)) then
                write(iulog,*) 'No. of values for saturation function is not equal to no. connections.'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             select case(satfunc_itype(sum_conn))
             case (0)
                ! Do nothing

             case (RELPERM_FUNC_MUALEM)
                if (set_upwind_auxvar(sum_conn)) then
                   conn_aux_vars(sum_conn)%satParams_up%relperm_func_type = RELPERM_FUNC_MUALEM
                   conn_aux_vars(sum_conn)%satParams_up%alpha   = param_1(sum_conn)
                   conn_aux_vars(sum_conn)%satParams_up%vg_m    = param_2(sum_conn)
                   conn_aux_vars(sum_conn)%satParams_up%sat_res = param_3(sum_conn)
                   conn_aux_vars(sum_conn)%satParams_up%vg_n    = 1.d0/(1.d0 - param_2(sum_conn))
                else
                   conn_aux_vars(sum_conn)%satParams_dn%relperm_func_type = RELPERM_FUNC_MUALEM
                   conn_aux_vars(sum_conn)%satParams_dn%alpha             = param_1(sum_conn)
                   conn_aux_vars(sum_conn)%satParams_dn%vg_m              = param_2(sum_conn)
                   conn_aux_vars(sum_conn)%satParams_dn%sat_res           = param_3(sum_conn)
                   conn_aux_vars(sum_conn)%satParams_dn%vg_n              = 1.d0/(1.d0 - param_2(sum_conn))
                end if
             case (RELPERM_FUNC_CAMPBELL)

                if (set_upwind_auxvar(sum_conn)) then
                   call SatFunc_Set_Campbell_RelPerm(conn_aux_vars(sum_conn)%satParams_up, &
                        param_1(sum_conn), param_2(sum_conn))
                else
                   call SatFunc_Set_Campbell_RelPerm(conn_aux_vars(sum_conn)%satParams_dn, &
                        param_1(sum_conn), param_2(sum_conn))
                endif

             case (RELPERM_FUNC_WEIBULL)

                if (set_upwind_auxvar(sum_conn)) then
                   call SatFunc_Set_Weibull_RelPerm(conn_aux_vars(sum_conn)%satParams_up, &
                        param_1(sum_conn), param_2(sum_conn))
                else
                   call SatFunc_Set_Weibull_RelPerm(conn_aux_vars(sum_conn)%satParams_dn, &
                        param_1(sum_conn), param_2(sum_conn))
                endif

             case (SAT_FUNC_CHUANG)

                if (set_upwind_auxvar(sum_conn)) then
                   call SatFunc_Set_Chuang(conn_aux_vars(sum_conn)%satParams_up, &
                        param_1(sum_conn), param_2(sum_conn))
                else
                   call SatFunc_Set_Chuang(conn_aux_vars(sum_conn)%satParams_dn, &
                        param_1(sum_conn), param_2(sum_conn))
                endif

             case (SAT_FUNC_VAN_GENUCHTEN)
                if (set_upwind_auxvar(sum_conn)) then
                   call SatFunc_Set_VG(conn_aux_vars(sum_conn)%satParams_up, &
                        param_3(sum_conn), param_2(sum_conn), param_1(sum_conn))
                else
                   call SatFunc_Set_VG(conn_aux_vars(sum_conn)%satParams_dn, &
                        param_3(sum_conn), param_2(sum_conn), param_1(sum_conn))
                endif

             case default
                write(iulog,*)'Unknown satfunc_itype_variable ', satfunc_itype(sum_conn)
                call endrun(msg=errMsg(__FILE__,__LINE__))
             end select

          enddo
          cur_cond => cur_cond%next
       enddo

    case default
       write(iulog,*)'Only supports AUXVAR_CONN_BC type'
       call endrun(msg=errMsg(__FILE__,__LINE__))

    end select

  end subroutine VSFMMPPSetSaturationFunctionAuxVarConn

#endif

end module MultiPhysicsProbVSFM
