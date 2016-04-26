
module MultiPhysicsProbThermal

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Object for thermal model
  !-----------------------------------------------------------------------

  ! !USES:
  use SoilStateType                      , only : soilstate_type
  use WaterstateType                     , only : waterstate_type
  use SoilHydrologyType                  , only : soilhydrology_type
  use clm_varctl                         , only : iulog
  use abortutils                         , only : endrun
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbBaseType           , only : multiphysicsprob_base_type
  use SystemOfEquationsThermalType       , only : sysofeqns_thermal_type
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

  type, public, extends(multiphysicsprob_base_type) :: mpp_thermal_type
     class(sysofeqns_thermal_type),pointer          :: sysofeqns
     type(sysofeqns_base_pointer_type), pointer  :: sysofeqns_ptr
   contains
     procedure, public :: Init                   => ThermalMPPInit
     procedure, public :: Setup                  => ThermalMPPSetup
  end type mpp_thermal_type

  type(mpp_thermal_type), public, target :: thermal_mpp

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine ThermalMPPInit(this)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbBaseType , only : MPPBaseInit
    use MultiPhysicsProbConstants, only : MPP_THERMAL_TBASED_KSP_CLM
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this

    call MPPBaseInit(this)

    this%name   = 'Snow + Standing water + Soil Thermal model using temperature'
    this%id     = MPP_THERMAL_TBASED_KSP_CLM
    this%nmesh  = 1

    nullify(this%sysofeqns)
    allocate(this%sysofeqns_ptr)
    nullify(this%sysofeqns_ptr%ptr)

  end subroutine ThermalMPPInit

  !------------------------------------------------------------------------
  subroutine ThermalMPPSetup(this, begg, endg, begc, endc, ncols_ghost, filter_hydrologyc, &
       xc_col, yc_col, zc_col, area_col, grid_owner, &
       soilstate_vars, waterstate_vars, soilhydrology_vars)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants, only : MESH_CLM_THERMAL_SOIL_COL
    use MultiPhysicsProbConstants, only : MESH_CLM_SNOW_COL
    use MultiPhysicsProbConstants, only : MESH_CLM_SSW_COL
    use MultiPhysicsProbConstants, only : PETSC_KSP
    use MultiPhysicsProbConstants, only : MPP_THERMAL_TBASED_KSP_CLM
    use MultiPhysicsProbConstants, only : SOE_THERMAL_TBASED
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this
    integer, intent(in)                  :: begg,endg
    integer, intent(in)                  :: begc,endc
    integer, intent(in)                  :: ncols_ghost          ! number of ghost/halo columns
    integer, intent(in)                  :: filter_hydrologyc(:) ! column filter for soil points
    PetscReal, pointer, intent(in)       :: xc_col(:)
    PetscReal, pointer, intent(in)       :: yc_col(:)
    PetscReal, pointer, intent(in)       :: zc_col(:)
    PetscReal, pointer, intent(in)       :: area_col(:)
    PetscInt, pointer, intent(in)        :: grid_owner(:)
    type(soilstate_type) , intent(in)    :: soilstate_vars
    type(waterstate_type)  , intent(in)  :: waterstate_vars
    type(soilhydrology_type), intent(in) :: soilhydrology_vars
    !
    ! !LOCAL VARIABLES:
    PetscInt                             :: soe_type

    ! Initialize
    call this%Init()

    ! Create all meshes needed by various governing equations
    this%nmesh = 3
    allocate(this%meshes(this%nmesh))

    this%id           = MPP_THERMAL_TBASED_KSP_CLM
    this%solver_type  = PETSC_KSP
    soe_type          = SOE_THERMAL_TBASED
    
    call this%meshes(1)%Create(MESH_CLM_SNOW_COL, begg, endg, begc, endc, &
         -1, ncols_ghost, &
         xc_col, yc_col, zc_col, area_col, grid_owner, &
         waterstate_vars, soilhydrology_vars)

    call this%meshes(2)%Create(MESH_CLM_SSW_COL, begg, endg, begc, endc, &
         -1, ncols_ghost, &
         xc_col, yc_col, zc_col, area_col, grid_owner, &
         waterstate_vars, soilhydrology_vars)
 
    call this%meshes(3)%Create(MESH_CLM_THERMAL_SOIL_COL, begg, endg, begc, endc, &
         -1, ncols_ghost, &
         xc_col, yc_col, zc_col, area_col, grid_owner, &
         waterstate_vars, soilhydrology_vars)

   ! Setup the system-of-equations
    allocate(this%sysofeqns)
    call this%sysofeqns%Setup(this%id, soe_type, this%meshes, this%nmesh)

    this%sysofeqns%solver_type = this%solver_type

    ! Setup the PETSc
    select case(this%solver_type)
    case (PETSC_KSP)
      call ThermalMPPKSPSetup(this, begc, endc, soilstate_vars)
    case default
       write(iulog,*) 'VSFMMPPSetup: Unknown this%solver_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermalMPPSetup

  !------------------------------------------------------------------------
  subroutine ThermalMPPKSPSetup(therm_mpp, begc, endc, soilstate_vars)
    !
    ! !DESCRIPTION:
    !
    use GoverningEquationBaseType
    use GoveqnThermalKSPTemperatureSnowType
    use GoveqnThermalKSPTemperatureSSWType
    use GoveqnThermalKSPTemperatureSoilType
    use SystemOfEquationsBasePointerType
    use abortutils, only : endrun
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type)                             :: therm_mpp
    integer, intent(in)                                 :: begc,endc
    type(soilstate_type) , intent(in)                   :: soilstate_vars
    !
    ! !LOCAL VARIABLES:
    PetscInt :: size
    PetscErrorCode                                        :: ierr
    DM                                                    :: dm_snow
    DM                                                    :: dm_sh2o
    DM                                                    :: dm_soil
    class(goveqn_base_type),pointer                       :: cur_goveq
    class (goveqn_thermal_ksp_temp_snow_type),pointer     :: goveq_snow
    class (goveqn_thermal_ksp_temp_ssw_type),pointer      :: goveq_sh2o
    class (goveqn_thermal_ksp_temp_soil_type),pointer     :: goveq_soil
    class(sysofeqns_thermal_type),pointer                 :: therm_soe

    therm_soe => therm_mpp%sysofeqns
    therm_mpp%sysofeqns_ptr%ptr => therm_mpp%sysofeqns

    ! Get pointers to governing-equations

    cur_goveq => therm_soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_thermal_ksp_temp_snow_type)
          goveq_snow => cur_goveq

       class is (goveqn_thermal_ksp_temp_ssw_type)
          goveq_sh2o => cur_goveq

       class is (goveqn_thermal_ksp_temp_soil_type)
          goveq_soil => cur_goveq
       end select

       cur_goveq => cur_goveq%next
    enddo

    ! DM-Composite approach


    ! Create PETSc DM for each governing equation
    size = goveq_snow%mesh%ncells_local
    call DMDACreate1d(PETSC_COMM_SELF,     &
                      DM_BOUNDARY_NONE,    &
                      size,                &
                      1,                   &
                      0,                   &
                      PETSC_NULL_INTEGER,  &
                      dm_snow,             &
                      ierr); CHKERRQ(ierr)
    call DMSetOptionsPrefix(dm_snow, "fsnow_", ierr); CHKERRQ(ierr)
    call DMDASetFieldName(dm_snow, 0, "snow_", ierr); CHKERRQ(ierr)
    call DMSetFromOptions(dm_snow, ierr); CHKERRQ(ierr)

    size = goveq_sh2o%mesh%ncells_local
    call DMDACreate1d(PETSC_COMM_SELF,     &
                      DM_BOUNDARY_NONE,    &
                      size,                &
                      1,                   &
                      0,                   &
                      PETSC_NULL_INTEGER,  &
                      dm_sh2o,             &
                      ierr); CHKERRQ(ierr)
    call DMSetOptionsPrefix(dm_sh2o, "fsh2o_", ierr); CHKERRQ(ierr)
    call DMDASetFieldName(dm_sh2o, 0, "sh2o_", ierr); CHKERRQ(ierr)
    call DMSetFromOptions(dm_sh2o, ierr); CHKERRQ(ierr)
    

    size = goveq_soil%mesh%ncells_local
    call DMDACreate1d(PETSC_COMM_SELF,     &
                      DM_BOUNDARY_NONE,    &
                      size,                &
                      1,                   &
                      0,                   &
                      PETSC_NULL_INTEGER,  &
                      dm_soil,             &
                      ierr); CHKERRQ(ierr)
    call DMSetOptionsPrefix(dm_soil, "fsoil_", ierr); CHKERRQ(ierr)
    call DMDASetFieldName(dm_soil, 0, "soil_", ierr); CHKERRQ(ierr)
    call DMSetFromOptions(dm_soil, ierr); CHKERRQ(ierr)


    ! Create DMComposite: temperature
    call DMCompositeCreate(PETSC_COMM_SELF, therm_soe%dm, ierr); CHKERRQ(ierr)
    call DMSetOptionsPrefix(therm_soe%dm, "temperature_", ierr); CHKERRQ(ierr)

    call DMCompositeAddDM(therm_soe%dm, dm_snow, ierr); CHKERRQ(ierr)
    call DMCompositeAddDM(therm_soe%dm, dm_sh2o, ierr); CHKERRQ(ierr)
    call DMCompositeAddDM(therm_soe%dm, dm_soil, ierr); CHKERRQ(ierr)

    call DMDestroy(dm_snow, ierr); CHKERRQ(ierr)
    call DMDestroy(dm_sh2o, ierr); CHKERRQ(ierr)
    call DMDestroy(dm_soil, ierr); CHKERRQ(ierr)

    call DMSetUp(therm_soe%dm, ierr); CHKERRQ(ierr)

    call DMCreateMatrix(therm_soe%dm, therm_soe%Amat, ierr); CHKERRQ(ierr)
    call MatSetOption(therm_soe%Amat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)
    call MatSetOption(therm_soe%Amat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)
    call MatSetFromOptions(therm_soe%Amat, ierr); CHKERRQ(ierr)

    ! Create solution vector
    call DMCreateGlobalVector(therm_soe%dm, therm_soe%soln         , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(therm_soe%dm, therm_soe%rhs          , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(therm_soe%dm, therm_soe%soln_prev    , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(therm_soe%dm, therm_soe%soln_prev_clm, ierr); CHKERRQ(ierr)

    call VecZeroEntries(therm_soe%soln, ierr); CHKERRQ(ierr)
    call VecZeroEntries(therm_soe%rhs,  ierr); CHKERRQ(ierr)
    call VecZeroEntries(therm_soe%soln_prev,  ierr); CHKERRQ(ierr)
    call VecZeroEntries(therm_soe%soln_prev_clm,  ierr); CHKERRQ(ierr)

    ! KSP
    call KSPCreate(PETSC_COMM_SELF, therm_soe%ksp, ierr); CHKERRQ(ierr)
    call KSPSetOptionsPrefix(therm_soe%ksp, "temperature_", ierr); CHKERRQ(ierr)
    call KSPSetComputeRHS(therm_soe%ksp, SOEComputeRHS, therm_mpp%sysofeqns_ptr, ierr); CHKERRQ(ierr)
    call KSPSetComputeOperators(therm_soe%ksp, SOEComputeOperators, therm_mpp%sysofeqns_ptr, ierr); CHKERRQ(ierr)
    call KSPSetFromOptions(therm_soe%ksp, ierr); CHKERRQ(ierr)

    ! Set soil thermal properties
    call MPPThermalSetSoils(therm_mpp, begc, endc, soilstate_vars)

  end subroutine ThermalMPPKSPSetup

  !------------------------------------------------------------------------
  subroutine MPPThermalSetSnow(therm_mpp, begc, endc)
    !
    ! !DESCRIPTION:
    !
    use GoverningEquationBaseType
    use GoveqnThermalKSPTemperatureSnowType
    use ThermalKSPTemperatureSnowAuxType
    use clm_varpar      , only : nlevgrnd, nlevsoi, nlevsno
    use ColumnType      , only : col
    use LandunitType    , only   : lun
    use landunit_varcon , only   : istsoil
    use clm_varcon      , only   : grav, denh2o
    use landunit_varcon , only : istcrop, istsoil
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type)                             :: therm_mpp
    integer              , intent(in)                   :: begc,endc
    !
    ! !LOCAL VARIABLES:
    class (goveqn_thermal_ksp_temp_snow_type) , pointer :: goveq_soil
    class(sysofeqns_thermal_type)             , pointer :: therm_soe
    class(goveqn_base_type)                   , pointer :: cur_goveq
    type (therm_ksp_temp_snow_auxvar_type)    , pointer :: aux_vars_in(:)
    PetscInt                                            :: j,c,g,l
    PetscInt                                            :: icell
    PetscBool                                           :: found

    therm_soe => therm_mpp%sysofeqns

    found = PETSC_FALSE
    cur_goveq => therm_soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_thermal_ksp_temp_snow_type)
          goveq_soil => cur_goveq
          found = PETSC_TRUE
       end select

       cur_goveq => cur_goveq%next
    enddo

    if (.not.found) then
       write(iulog,*)'Thermal governing equation for snow NOT found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    aux_vars_in => goveq_soil%aux_vars_in

    ! Set thermal properties

    icell = 0
    do c = begc, endc

       ! Soil layers
       do j = -nlevsno+1, 0

         icell = icell + 1

         g = col%gridcell(c)
         l = col%landunit(c)

         if (col%active(c) .and. .not.lun%lakpoi(l) .and. .not.lun%urbpoi(l)) then
            aux_vars_in(icell)%is_active = PETSC_TRUE
         else
            aux_vars_in(icell)%is_active = PETSC_FALSE
         endif

      enddo
    enddo

  end subroutine MPPThermalSetSnow

  !------------------------------------------------------------------------
  subroutine MPPThermalSetSoils(therm_mpp, begc, endc, soilstate_vars)
    !
    ! !DESCRIPTION:
    !
    use GoverningEquationBaseType
    use GoveqnThermalKSPTemperatureSoilType
    use ThermalKSPTemperatureSoilAuxType
    use clm_varpar      , only : nlevgrnd, nlevsoi, nlevsno
    use ColumnType      , only : col
    use LandunitType    , only   : lun
    use landunit_varcon , only   : istsoil
    use clm_varcon      , only   : grav, denh2o
    use landunit_varcon , only : istcrop, istsoil
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type)                             :: therm_mpp
    integer              , intent(in)                   :: begc,endc
    type(soilstate_type) , intent(in)                   :: soilstate_vars
    !
    ! !LOCAL VARIABLES:
    class (goveqn_thermal_ksp_temp_soil_type) , pointer :: goveq_soil
    class(sysofeqns_thermal_type)             , pointer :: therm_soe
    class(goveqn_base_type)                   , pointer :: cur_goveq
    type (therm_ksp_temp_soil_auxvar_type)    , pointer :: aux_vars_in(:)
    PetscInt                                            :: j,c,g,l
    PetscInt                                            :: icell
    PetscInt                                            :: col_id
    PetscInt                                            :: first_active_col_id
    PetscBool                                           :: found

    associate(&
         watsat       =>    soilstate_vars%watsat_col        , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         csol         =>    soilstate_vars%csol_col          , & ! Input:  [real(r8) (:,:) ]  heat capacity, soil solids (J/m**3/Kelvin)
         tkmg         =>    soilstate_vars%tkmg_col          , & ! Input:  [real(r8) (:,:) ]  thermal conductivity, soil minerals  [W/m-K]
         tkdry        =>    soilstate_vars%tkdry_col           & ! Input:  [real(r8) (:,:) ]  thermal conductivity, dry soil (W/m/Kelvin)
         )

    therm_soe => therm_mpp%sysofeqns

    found = PETSC_FALSE
    cur_goveq => therm_soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_thermal_ksp_temp_soil_type)
          goveq_soil => cur_goveq
          found = PETSC_TRUE
       end select

       cur_goveq => cur_goveq%next
    enddo

    if (.not.found) then
       write(iulog,*)'Thermal governing equation for soil NOT found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    aux_vars_in => goveq_soil%aux_vars_in

    first_active_col_id = -1
    do c = begc, endc
       l = col%landunit(c)

       if (col%active(c) .and. .not.lun%lakpoi(l) .and. .not.lun%urbpoi(l)) then
          if (first_active_col_id == -1) then
             first_active_col_id = c
             exit
          endif
       endif
    enddo

    if (first_active_col_id == -1) then
       write(iulog,*)'No active soil column found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Set thermal properties

    icell = 0
    do c = begc, endc

       ! Soil layers
      do j = 1, nlevgrnd

         icell = icell + 1

         g = col%gridcell(c)
         l = col%landunit(c)

         if (col%active(c) .and. .not.lun%lakpoi(l) .and. .not.lun%urbpoi(l)) then
            col_id                       = c
            aux_vars_in(icell)%is_active = PETSC_TRUE
         else
            col_id                       = first_active_col_id 
            aux_vars_in(icell)%is_active = PETSC_FALSE
         endif

         if (j > nlevsoi) then
            aux_vars_in(icell)%is_soil_shallow = PETSC_FALSE
         else
            aux_vars_in(icell)%is_soil_shallow = PETSC_TRUE
         endif         

         aux_vars_in(icell)%itype                 = lun%itype(l)
         aux_vars_in(icell)%por                   = watsat(col_id,j)
         aux_vars_in(icell)%therm_cond_minerals   = tkmg(col_id,j)
         aux_vars_in(icell)%therm_cond_dry        = tkdry(col_id,j)
         aux_vars_in(icell)%heat_cap_minerals_puv = csol(col_id,j)

      enddo
    enddo

    end associate
  end subroutine MPPThermalSetSoils

#endif

end module MultiPhysicsProbThermal
