module EOS_Gas_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use EOSData_module
  use Utility_module, only : Equal
  
  implicit none

  private
  
  ! module variables
  PetscReal :: fmw_gas           !kg/Kmol
  PetscReal :: constant_density
  PetscReal :: constant_enthalpy
  PetscReal :: constant_viscosity
  PetscReal :: constant_henry
  !density at reference P and T
  PetscReal :: reference_density_kg
  ! specific to RKS EOS
  PetscBool :: rks_use_hydrogen
  PetscBool :: rks_use_effect_critical_props
  PetscBool :: rks_use_cubic_root_solution
  PetscReal :: rks_Tc
  PetscReal :: rks_Pc
  PetscReal :: rks_acentric
  PetscReal :: rks_coeff_a
  PetscReal :: rks_coeff_b
  
  ! exponential
  PetscReal :: exponent_reference_density
  PetscReal :: exponent_reference_pressure
  PetscReal :: exponent_gas_compressibility

  ! This is the offset added to temperature [C] used to calculate the energy
  ! equation of state.  Swithing between 0. and 273.15 greatly changes results.
#if defined(MATCH_TOUGH2)
  PetscReal, parameter :: T_energy_offset = 0.d0
#else
  PetscReal, parameter :: T_energy_offset = 273.15d0
#endif

  ! EOS databases
  class(eos_database_type), pointer :: eos_dbase

  ! PVT tables - eos_tables
  class(eos_table_type), pointer :: pvt_table

  ! In order to support generic EOS subroutines, we need the following:
  ! 1. An interface declaration that defines the argument list (best to have 
  !    "Dummy" appended.
  ! 2. A procedure pointer that is initially set to null.  This pointer is
  !    pointed to the appropriate subroutine later on (e.g. EOSGasInit())
  ! 3. An interface for derivative/non-derivative versions

  ! procedure pointer declarations
  procedure(EOSGasViscosityDummy), pointer :: EOSGasViscosityPtr => null()
  procedure(EOSGasDensityEnergyDummy), pointer :: &
    EOSGasDensityEnergyPtr => null()
  procedure(EOSGasDensityDummy), pointer :: EOSGasDensityPtr => null()
  procedure(EOSGasEnergyDummy), pointer :: EOSGasEnergyPtr => null()
  procedure(EOSGasHenryDummy), pointer :: EOSGasHenryPtr => null()
  
  ! interface blocks
  interface
    subroutine EOSGasViscosityDummy(T, P_comp, P_gas, Rho_comp, V_mix, &
                                    calculate_derivative, dV_dT, dV_dPcomp, &
                                    dV_dPgas, dV_dRhocomp, ierr, table_idxs)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P_comp   ! air pressure [Pa]
      PetscReal, intent(in) :: P_gas    ! gas pressure [Pa]
      PetscReal, intent(in) :: Rho_comp ! air density [C]
      PetscReal, intent(out) :: V_mix   ! mixture viscosity
      PetscBool, intent(in) :: calculate_derivative
      PetscReal, intent(out) :: dV_dT       ! derivative wrt temperature
      PetscReal, intent(out) :: dV_dPcomp   ! derivative wrt component pressure
      PetscReal, intent(out) :: dV_dPgas    ! derivative wrt gas pressure
      PetscReal, intent(out) :: dV_dRhocomp ! derivative wrt component density      
      PetscErrorCode, intent(out) :: ierr
      PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
    end subroutine EOSGasViscosityDummy
    subroutine EOSGasDensityDummy(T,P,Rho_gas,dRho_dT,dRho_dP,ierr,table_idxs)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P        ! pressure [Pa]
      PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
      PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
      PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressuret
      PetscErrorCode, intent(out) :: ierr
      PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
    end subroutine EOSGasDensityDummy
    subroutine EOSGasEnergyDummy(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P        ! pressure [Pa]
      PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
      PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
      PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
      PetscReal, intent(out) :: U       ! internal energy [J/kmol]
      PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
      PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSGasEnergyDummy
    subroutine EOSGasDensityEnergyDummy(T,P,Rho_gas,dRho_dT,dRho_dP, &
                                        H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr, &
                                        table_idxs)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P        ! pressure [Pa]
      PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
      PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
      PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressuret
      PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
      PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
      PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
      PetscReal, intent(out) :: U       ! internal energy [J/kmol]
      PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
      PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
      PetscErrorCode, intent(out) :: ierr
      PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
    end subroutine EOSGasDensityEnergyDummy
    subroutine EOSGasHenryDummy(T,Psat,Hc,calculate_derivative, &
                                Psat_p,Psat_T,Hc_P,Hc_T)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: Psat     ! saturation pressure
      PetscReal, intent(out) :: Hc      ! Henry's constant
      PetscBool, intent(in) :: calculate_derivative
      PetscReal, intent(in) :: Psat_P   ! derivative Psat wrt pressure
      PetscReal, intent(in) :: Psat_T   ! derivative Psat wrt temperature
      PetscReal, intent(out) :: Hc_P    ! derivative Henry's constant wrt pressure
      PetscReal, intent(out) :: Hc_T    ! derivative Henry's constant wrt temperature
    end subroutine EOSGasHenryDummy
  end interface
  
  ! interfaces for derivative/non-derivative versions that are visible outside
  ! the module.
  interface EOSGasViscosity
    procedure EOSGasViscosityNoDerive
    procedure EOSGasViscosityDerive
  end interface
  interface EOSGasDensity
    procedure EOSGasDensityNoDerive
    procedure EOSGasDensityDerive
  end interface
  interface EOSGasEnergy
    procedure EOSGasEnergyNoDerive
    procedure EOSGasEnergyDerive
  end interface
  interface EOSGasDensityEnergy
    procedure EOSGasDenEnthNoDerive
    procedure EOSGasDenEnthDerive
  end interface
  interface EOSGasHenry
    procedure EOSGasHenryNoDerive
    procedure EOSGasHenryDerive
  end interface

  ! the "public" definition that makes subroutines visible outside.
  public :: EOSGasInit, &
            EOSGasVerify, &
            EOSGasViscosity, &
            EOSGasDensity, &
            EOSGasEnergy, &
            EOSGasDensityEnergy, &
            EOSGasHenry, &
            EOSGasInputRecord, &
            EOSGasTest
            
  public :: EOSGasSetDensityIdeal, &
            EOSGasSetEnergyIdeal, &
            EOSGasSetDensityRKS, &
            EOSGasSetDensityPRMethane, &
            EOSGasSetEnergyIdealMethane, &
            EOSGasSetFMWConstant, &
            EOSGasGetFMW, & 
            EOSGasSetDensityConstant, &
            EOSGasSetEnergyConstant, &
            EOSGasSetViscosityConstant, &
            EOSGasSetHenry, &
            EOSGasSetHenryConstant, &
            EOSGasSetEOSDBase, &
            EOSGasSetReferenceDensity, &
            EOSGasGetReferenceDensity, &
            EOSGasSetPVDG, &
            EOSGasTableProcess, &
            EOSGasDBaseDestroy
 
  contains

! ************************************************************************** !

subroutine EOSGasInit()

  implicit none

  fmw_gas = UNINITIALIZED_DOUBLE  
  constant_density = UNINITIALIZED_DOUBLE
  constant_viscosity = UNINITIALIZED_DOUBLE
  constant_enthalpy = UNINITIALIZED_DOUBLE
  constant_henry = UNINITIALIZED_DOUBLE

  ! exponential
  exponent_reference_density = UNINITIALIZED_DOUBLE
  exponent_reference_pressure = UNINITIALIZED_DOUBLE
  exponent_gas_compressibility = UNINITIALIZED_DOUBLE

  ! for RKS gas density
  rks_use_hydrogen = PETSC_FALSE
  rks_use_effect_critical_props = PETSC_FALSE
  rks_Tc = UNINITIALIZED_DOUBLE
  rks_Pc = UNINITIALIZED_DOUBLE
  rks_acentric = UNINITIALIZED_DOUBLE
  rks_coeff_a = UNINITIALIZED_DOUBLE
  rks_coeff_b = UNINITIALIZED_DOUBLE

  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasDensityPtr => EOSGasDensityIdeal
  EOSGasEnergyPtr => EOSGasEnergyIdeal
  EOSGasViscosityPtr => EOSGasViscosity1
  EOSGasHenryPtr => EOSGasHenry_air
  
end subroutine EOSGasInit

! ************************************************************************** !

subroutine EOSGasVerify(ierr,error_string)

  implicit none
  
  PetscErrorCode, intent(out) :: ierr
  character(len=MAXSTRINGLENGTH), intent(out) :: error_string
  
  ierr = 0
  
  error_string = ''
  if ((associated(EOSGasDensityPtr,EOSGasDensityIdeal) .and. &
        Initialized(constant_density)) .or. &
      (associated(EOSGasDensityPtr,EOSGasDensityRKS) .and. &
        Initialized(constant_density)) .or. &
      (associated(EOSGasEnergyPtr,EOSGasEnergyIdeal) .and. &
        Initialized(constant_enthalpy)) &
     ) then
    ierr = 1
  endif

  if (associated(EOSGasDensityPtr,EOSGasDensityConstant) .and. &
      Uninitialized(constant_density)) then
    error_string = trim(error_string) // &
      ' CONSTANT density not set.'
    ierr = 1
  endif
  
  if (associated(EOSGasEnergyPtr,EOSGasEnergyConstant) .and. &
      Uninitialized(constant_enthalpy)) then
    error_string = trim(error_string) // &
      ' CONSTANT enthalpy not set.'
    ierr = 1
  endif
  
  if ((associated(EOSGasViscosityPtr, &
                  EOSGasViscosityConstant) .and. &
       Uninitialized(constant_viscosity)) .or. &
      (associated(EOSGasViscosityPtr, &
                  EOSGasViscosity1) .and. &
       Initialized(constant_viscosity))) then
    ierr = 1
  endif
  
  if (associated(EOSGasHenryPtr, &
                 EOSGasHenryConstant) .and. &
      Uninitialized(constant_henry)) then
    error_string = trim(error_string) // " Henry's constant not set"
    ierr = 1
  endif
  
  if (associated(EOSGasDensityPtr,EOSGasDensityRKS)) then
    if (rks_use_hydrogen) then
      ! Assign default hydrogen parameters if not assigned
      if (Uninitialized(rks_Tc)) then
        rks_Tc = 41.67d0
        ierr = 5
        error_string = trim(error_string) // " Tc"
      endif
      if (Uninitialized(rks_Pc)) then
        rks_Pc = 2.1029d6
        ierr = 5
        error_string = trim(error_string) // " Pc"
      endif
      if (Uninitialized(rks_coeff_a)) then
        rks_coeff_a = 0.42747d0
        ierr = 5
        error_string = trim(error_string) // " omega_a"
      endif
      if (Uninitialized(rks_coeff_b)) then
        rks_coeff_b = 0.08664d0
        ierr = 5
        error_string = trim(error_string) // " omega_b"
      endif
    else
      if (Uninitialized(rks_Tc) .or. Uninitialized(rks_coeff_a) .or. &
        Uninitialized(rks_Pc) .or. Uninitialized(rks_coeff_b) .or. &
        Uninitialized(rks_acentric)) then
          error_string = trim(error_string) // &
          " RKS parameters not set for non-hydrogen gas"
          ierr = 1
      endif
    endif
  endif

  if (Uninitialized(fmw_gas)) then
    fmw_gas = FMWAIR
    !ierr = 6
    !error_string = trim(error_string) // " FMWAIR"
  end if

  if ( associated(pvt_table) ) then
    if(Uninitialized(reference_density_kg)) then
      error_string = trim(error_string) // &
      'A reference (e.g. Surface) density must be specified ' // &
      'using either REFERENCE_DENSITY, SURFACE_DENSITY or STANDARD_DENSITY '
      ierr = 1
    end if
  end if
      
end subroutine EOSGasVerify

! ************************************************************************** !

subroutine EOSGasSetDensityIdeal()

  implicit none
  
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasDensityPtr => EOSGasDensityIdeal
  
end subroutine EOSGasSetDensityIdeal

! ************************************************************************** !

subroutine EOSGasSetDensityRKS(use_hydrogen,use_effective_properties, &
                               use_cubic_root_solution, &
                               tc,pc,acen,omegaa,omegab)

  implicit none

  PetscBool :: use_hydrogen
  PetscBool :: use_effective_properties
  PetscBool :: use_cubic_root_solution
  PetscReal :: tc
  PetscReal :: pc
  PetscReal :: acen
  PetscReal :: omegaa
  PetscReal :: omegab
  
  rks_use_hydrogen = use_hydrogen
  rks_use_effect_critical_props = use_effective_properties
  rks_use_cubic_root_solution = use_cubic_root_solution
  rks_Tc = tc
  rks_Pc = pc
  rks_acentric = acen
  rks_coeff_a = omegaa
  rks_coeff_b = omegab
  
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasDensityPtr => EOSGasDensityRKS
  
end subroutine EOSGasSetDensityRKS

! ************************************************************************** !

subroutine EOSGasSetDensityPRMethane()

  implicit none
  
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasDensityPtr => EOSGasDensityPRMethane
  
end subroutine EOSGasSetDensityPRMethane

! ************************************************************************** !

subroutine EOSGasSetEnergyIdeal()

  implicit none
  
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasEnergyPtr => EOSGasEnergyIdeal
  
end subroutine EOSGasSetEnergyIdeal

! ************************************************************************** !

subroutine EOSGasSetEnergyIdealMethane()

  implicit none
  
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasEnergyPtr => EOSGasEnergyIdealMethane
  
end subroutine EOSGasSetEnergyIdealMethane

! ************************************************************************** !

subroutine EOSGasSetFMWConstant(input_fmw_gas)

  implicit none

  PetscReal :: input_fmw_gas
  
  fmw_gas = input_fmw_gas

end subroutine EOSGasSetFMWConstant

! ************************************************************************** !

function EOSGasGetFMW()

  implicit none
  
  PetscReal :: EOSGasGetFMW
  
  EOSGasGetFMW = fmw_gas
  
end function EOSGasGetFMW

! ************************************************************************** !

subroutine EOSGasSetReferenceDensity(input_ref_density_kg)

  implicit none

  PetscReal :: input_ref_density_kg

  reference_density_kg = input_ref_density_kg

end subroutine EOSGasSetReferenceDensity

! ************************************************************************** !

function EOSGasGetReferenceDensity()

  implicit none

  PetscReal :: EOSGasGetReferenceDensity

  EOSGasGetReferenceDensity = reference_density_kg

end function EOSGasGetReferenceDensity

! ************************************************************************** !

subroutine EOSGasSetDensityConstant(density)

  implicit none
  
  PetscReal :: density
  
  constant_density = density  
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasDensityPtr => EOSGasDensityConstant
  
end subroutine EOSGasSetDensityConstant

! ************************************************************************** !

subroutine EOSGasSetEnergyConstant(enthalpy)

  implicit none
  
  PetscReal :: enthalpy
  
  constant_enthalpy = enthalpy  
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasEnergyPtr => EOSGasEnergyConstant
  
end subroutine EOSGasSetEnergyConstant

! ************************************************************************** !

subroutine EOSGasSetViscosityConstant(viscosity)

  implicit none
  
  PetscReal :: viscosity
  
  constant_viscosity = viscosity  
  EOSGasViscosityPtr => EOSGasViscosityConstant
  
end subroutine EOSGasSetViscosityConstant

! ************************************************************************** !

subroutine EOSGasSetHenry()

  implicit none
  
  EOSGasHenryPtr => EOSGasHenry_air
  
end subroutine EOSGasSetHenry

! ************************************************************************** !

subroutine EOSGasSetHenryConstant(henrys_constant)

  implicit none
  
  PetscReal :: henrys_constant
  
  constant_henry = henrys_constant
  EOSGasHenryPtr => EOSGasHenryConstant
  
end subroutine EOSGasSetHenryConstant

! ************************************************************************** !

subroutine EOSGasViscosityNoDerive(T, P_comp, P_gas, Rho_comp, V_mix, ierr, &
                                   table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P_comp   ! air pressure [Pa]
  PetscReal, intent(in) :: P_gas    ! gas pressure [Pa]
  PetscReal, intent(in) :: Rho_comp ! air density [C]
  PetscReal, intent(out) :: V_mix   ! mixture viscosity
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
  
  PetscReal :: dum1, dum2, dum3, dum4
  
  call EOSGasViscosityPtr(T, P_comp, P_gas, Rho_comp, V_mix, PETSC_FALSE, &
                          dum1,dum2,dum3,dum4,ierr,table_idxs)
  
end subroutine EOSGasViscosityNoDerive

! ************************************************************************** !

subroutine EOSGasViscosityDerive(T, P_comp, P_gas, Rho_comp, &
                                 dRho_dT, dRho_dPcomp, dRho_dPgas, &
                                 dPcomp_dT, dPcomp_dPgas, &
                                 V_mix, dV_dT, dV_dPcomp, dV_dPgas, ierr, &
                                 table_idxs)
  implicit none

  PetscReal, intent(in) :: T             ! temperature [C]
  PetscReal, intent(in) :: P_comp        ! air pressure [Pa]
  PetscReal, intent(in) :: P_gas         ! gas pressure [Pa]
  PetscReal, intent(in) :: Rho_comp      ! air density [C]
  PetscReal, intent(in) :: dRho_dT       ! derivative of density wrt temperature
  PetscReal, intent(in) :: dRho_dPcomp   ! derivative of density wrt gas comp pressure
  PetscReal, intent(in) :: dRho_dPgas    ! derivative of density wrt gas pressure
  PetscReal, intent(in) :: dPcomp_dT     ! derivative of gas comp pressure wrt temperature
  PetscReal, intent(in) :: dPcomp_dPgas  ! derivative of gas comp pressure wrt gas pressure
  PetscReal, intent(out) :: V_mix        ! mixture viscosity
  PetscReal, intent(out) :: dV_dT        ! derivative gas viscosity wrt temperature
  PetscReal, intent(out) :: dV_dPcomp    ! derivative gas viscosity wrt gas comp press
  PetscReal, intent(out) :: dV_dPgas     ! derivative gas viscosity wrt gas press
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
  
  !geh: at very low temperatures, the derivative wrt Rhocomp is very sensitive to
  !     the perturbation.  Need a value as large as 1.d-3 at 2C to match analtyical.
  PetscReal, parameter :: pert_tol = 1.d-8
  PetscReal :: pert
    
  PetscReal :: T_pert
  PetscReal :: P_comp_pert
  PetscReal :: P_gas_pert
  PetscReal :: Rho_comp_pert
  PetscReal :: V_mix_pert
  PetscReal :: dV_dRhocomp
  PetscReal :: dV_dT_, dV_dPcomp_, dV_dPgas_, dV_dRhocomp_
  PetscReal :: dum1, dum2, dum3, dum4
  
  dV_dT = 0.d0
  dV_dPgas = 0.d0
  dV_dPcomp = 0.d0
  dV_dRhocomp = 0.d0
  
!#define NUMERICAL_DERIVATIVE_VISCOSITY
#define PARTIALS
    
  ! We have to calcualte the derivative numerically
  call EOSGasViscosityPtr(T, P_comp, P_gas, Rho_comp, V_mix, PETSC_TRUE, &
                          dV_dT_, dV_dPcomp_, dV_dPgas_, dV_dRhocomp_, ierr, &
                          table_idxs)
#if defined(PARTIALS)
  dV_dT_ = dV_dT_ + dV_dRhocomp_*dRho_dT + dV_dPcomp_*dPcomp_dT
  dV_dPgas_ = dV_dPgas_ + dV_dPcomp_*dPcomp_dPgas + dV_dRhocomp_*dRho_dPgas
  ! put dv_dPcomp last to avoid double counting for partials above
  dV_dPcomp_ = dV_dPcomp_ + dV_dRhocomp_*dRho_dPcomp
#endif

#if defined(NUMERICAL_DERIVATIVE_VISCOSITY)                          
  ! temperature
  pert = pert_tol * T
  T_pert = T + pert
  call EOSGasViscosityPtr(T_pert, P_comp, P_gas, Rho_comp, V_mix_pert, &
                          PETSC_FALSE,dum1,dum2,dum3,dum4,ierr,table_idxs)
  dV_dT = dV_dT + (V_mix_pert - V_mix)/pert
  ! gas component pressure
  pert = pert_tol * P_comp
  P_comp_pert = P_comp + pert
  call EOSGasViscosityPtr(T, P_comp_pert, P_gas, Rho_comp, V_mix_pert, &
                          PETSC_FALSE,dum1,dum2,dum3,dum4,ierr,table_idxs)
  dV_dPcomp = (V_mix_pert - V_mix)/pert
  ! gas pressure
  pert = pert_tol * P_gas
  P_gas_pert = P_gas + pert
  call EOSGasViscosityPtr(T, P_comp, P_gas_pert, Rho_comp, V_mix_pert, &
                          PETSC_FALSE,dum1,dum2,dum3,dum4,ierr,table_idxs)
  dV_dPgas = (V_mix_pert - V_mix)/pert
  ! component density
  pert = pert_tol * Rho_comp
  Rho_comp_pert = Rho_comp + pert
  call EOSGasViscosityPtr(T, P_comp, P_gas, Rho_comp_pert, V_mix_pert, &
                          PETSC_FALSE,dum1,dum2,dum3,dum4,ierr,table_idxs)
  dV_dRhocomp = (V_mix_pert - V_mix)/pert
#if defined(PARTIALS)
  dV_dT = dV_dT + dV_dRhocomp * dRho_dT + dV_dPcomp * dPcomp_dT
  dV_dPgas = dV_dPgas + dV_dPcomp * dPcomp_dPgas + dV_dRhocomp * dRho_dPgas
  dV_dPcomp = dV_dPcomp + dV_dRhocomp * dRho_dPcomp
#endif
#else
  dV_dPcomp = dV_dPcomp_
  dV_dT = dV_dT_
  dV_dPgas = dV_dPgas_
#endif

end subroutine EOSGasViscosityDerive

! ************************************************************************** !

subroutine EOSGasViscosity1(T, P_comp, P_gas, Rho_comp, V_mix, &
                            calculate_derivative, dV_dT, dV_dPcomp, &
                            dV_dPgas, dV_dRhocomp, ierr, table_idxs)

  implicit none

  PetscReal, intent(in) :: T            ! temperature [C]
  PetscReal, intent(in) :: P_comp       ! air pressure [Pa]
  PetscReal, intent(in) :: P_gas        ! gas pressure [Pa]
  PetscReal, intent(in) :: Rho_comp     ! air density [C]
  PetscReal, intent(out) :: V_mix       ! mixture viscosity
  PetscBool, intent(in) :: calculate_derivative
  PetscReal, intent(out) :: dV_dT       ! derivative wrt temperature
  PetscReal, intent(out) :: dV_dPcomp   ! derivative wrt component pressure
  PetscReal, intent(out) :: dV_dPgas    ! derivative wrt gas pressure
  PetscReal, intent(out) :: dV_dRhocomp ! derivative wrt component density
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
  
  !geh: copied from gas_eos_mod.F90
  !
  ! REFERENCES
  ! THIS ROUTINE IS LARGELY ADAPTED FROM THE TOUGH CODE.
  ! this routine computes the viscosity of vapor-air mixtures.
  ! it uses a modified version of a formulation based on kinetic
  ! gas theory, as given by j.o. hirschfelder, c.f. curtiss, and
  ! r.b. bird, molecular theory of gases and liquids, john wiley
  ! & sons, 1954, pp. 528-530.
  ! the modification made to the hirschfelder et al. expressions is
  ! that for vapor viscosity accurate (empirical) values are used,
  ! rather than the first order expression of kinetic theory.
  ! the formulation matches experimental data on viscosities of
  ! vapor-air mixtures in the temperature range from 100 to 150
  ! deg. c, for all compositions, to better than 4%.
  ! 
  !PetscReal, intent(in) :: t     ! [C]
  !PetscReal, intent(in) :: p_air ! [Pa]
  !PetscReal, intent(in) :: p_gas ! [Pa]
  !PetscReal, intent(in) :: d_air ! [kmol/m^3]
  !PetscReal, intent(out) :: visg ! [Pa-s]

  PetscReal ::  fair,fwat,cair,cwat
  PetscReal :: p_air, d_air, visg

  data  fair,   fwat,    cair,  cwat &
        /97.d0, 363.d0, 3.617d0, 2.655d0/
 
  PetscReal :: fmix,cmix,d,xga,xg1,tk,trd1,trd3,ome1,ome3,ard,fmw3,vis1, &
               v1,vs,vis2,vis3,z1,g,h,e,z2,z3
               
  PetscReal :: dtrd1_dT, dtrd3_dT
  PetscReal :: dome1_dtrd1, dome1_dT
  PetscReal :: dome3_dtrd3, dome3_dT
  PetscReal :: dard_dtrd3, dard_dT
  PetscReal :: dvis1_dtrd1, dvis1_dome1, dvis1_dT
  PetscReal :: dvis2_dvs, dvis2_dT, dvis2_dRhocomp
  PetscReal :: dvis3_dtrd3, dvis3_dome3, dvis3_dT
  PetscReal :: dv1_dT
  PetscReal :: dvs_dT, dvs_dd, dvs_dRhocomp
  PetscReal :: dd_dRhocomp
  PetscReal :: dxga_dPcomp, dxga_dPgas
  PetscReal :: dxg1_dPcomp, dxg1_dPgas
  PetscReal :: dg_dxga, dg_dPcomp, dg_dPgas
  PetscReal :: de_dxga, de_dxg1, de_dvis1, de_dvis2, de_dvis3, &
               de_dT, de_dRhocomp, de_dPcomp, de_dPgas
  PetscReal :: dh_dxg1, dh_dPcomp, dh_dPgas
  PetscReal :: dz1_dxga, dz1_dxg1, dz1_dvis1, dz1_dvis2, dz1_dvis3, &
               dz1_dT, dz1_dRhocomp, dz1_dPcomp, dz1_dPgas
  PetscReal :: dz2_dard, dz2_de, dz2_dg, dz2_dh, dz2_dvis1, dz2_dvis2, &
               dz2_dvis3, dz2_dT, dz2_dRhocomp, dz2_dPcomp, dz2_dPgas
  PetscReal :: dz3_dard, dz3_de, dz3_dg, dz3_dh, dz3_dvis1, dz3_dvis2, &
               dz3_dvis3, dz3_dT, dz3_dRhocomp, dz3_dPcomp, dz3_dPgas, &
               dz3_dxga, dz3_dxg1
  PetscReal :: dvisg_dz1, dvisg_dz2, dvisg_dz3
  PetscReal :: ard_0point6, one_over_vis1sq, one_over_vis2sq, z1pz2
  PetscReal :: one_over_trd3
  PetscReal :: tempreal
  
!c======================================================================

  p_air = P_comp
  d_air = Rho_comp
          
  fmix = sqrt (fair*fwat)
  cmix = (cair+cwat)*0.5d0

!      do k = 1,nb
!       if (iphas(k).eq.2 .or. iphas(k).eq.0) then

      d = d_air * FMWAIR   
      xga = p_air / P_gas
      xg1 = 1.D0 - xga
      tk  = t + 273.15d0
      trd1 = tk/fair
      trd3 = tk/fmix
      one_over_trd3 = 1.d0/trd3
      ome1 = 1.188d0/trd1-0.051d0
      ome3 = (1.480d0-0.412d0*log(trd3))*one_over_trd3
      ard  = 1.095d0*one_over_trd3
      fmw3 = 2.d0*FMWAIR*FMWH2O/(FMWAIR+FMWH2O)
      vis1 = 266.93d-7*sqrt(FMWAIR*trd1*fair)/(cair*cair*ome1*trd1)
      v1 = .407d0*t +80.4d0
      if (t <= 350.d0) then
        vs = 1.d-7*(v1-d*(1858.d0-5.9d0*t )*1.d-3)
      else
!cpcl .      vs = 1.d-7*(v1 + 0.353d0*d + 676.5d-6*d**2 + 102.1d-9*d**3)
        vs = 1.d-7*(v1 + (0.353d0 + (676.5d-6 + 102.1d-9*d)*d)*d)
      endif
      vis2 = 10.d0*vs
      vis3 = 266.93d-7*sqrt(fmw3*trd3*fmix)/(cmix*cmix*ome3*trd3)
      z1 = xga*xga/vis1+2.d0*xg1*xga/vis3+xg1*xg1/vis2
      g = xga*xga*FMWAIR/FMWH2O
      h = xg1*xg1*FMWH2O/FMWAIR
      e = (2.d0*xga*xg1*FMWAIR*FMWH2O/fmw3**2.d0)*vis3/(vis1*vis2)
      z2 = 0.6d0*ard*(g/vis1+e+h/vis2)
      z3 = 0.6d0*ard*(g+e*(vis1+vis2)-2.d0*xga*xg1+h)
      visg  = (1.d0+z3)/(z1+z2)*.1d0 
      
  V_mix = visg
  
  if (calculate_derivative) then

    if (xg1 < 1.d-100) then
      print *, 'Derivatives in EOSGasViscosity1() may not be calculated &
        &when the water vapor pressure is essentially zero.'
      stop
    endif
  
    ! air density (kg)
    ! d = d_air * FMWAIR   
    dd_dRhocomp = FMWAIR
  
    ! xga
    ! xga = p_air / P_gas
    dxga_dPcomp = 1.d0 / P_gas
    dxga_dPgas = -1.d0 * p_air / (P_gas*P_gas)
  
    ! xg1
    ! xg1 = 1.D0 - xga
    dxg1_dPcomp = -1.d0*dxga_dPcomp
    dxg1_dPgas = -1.d0*dxga_dPgas
  
    ! trd3 wrt t
    ! trd1 = tk/fair
    dtrd1_dT = 1.d0/fair

    ! trd3 wrt t
    ! trd3 = tk/fmix
    dtrd3_dT = 1.d0/fmix
  
    ! ome1 wrt trd1
    ! ome1 = 1.188d0/trd1-0.051d0
    dome1_dtrd1 = -1.188d0/(trd1*trd1)
    dome1_dT = dome1_dtrd1*dtrd1_dT
  
    ! ome3 wrt trd3
    ! ome3 = (1.480d0-0.412d0*log(trd3))/trd3
    dome3_dtrd3 = (-0.412d0*one_over_trd3 - ome3)*one_over_trd3 ! dlnx / dx = 1/x dx
    dome3_dT = dome3_dtrd3*dtrd3_dT
  
    ! ard wrt trd3
    ! ard = 1.095d0/trd3
    dard_dtrd3 = -1.095d0*one_over_trd3*one_over_trd3
    dard_dT = dard_dtrd3*dtrd3_dT
  
    ! vis1 wrt trd1, dome1
    ! vis1 = 266.93d-7*sqrt(FMWAIR*trd1*fair)/(cair*cair*ome1*trd1)
    tempreal = FMWAIR*trd1*fair
    dvis1_dtrd1 = (0.5d0*vis1 - vis1)/trd1
    dvis1_dome1 = -1.d0*vis1/ome1
    dvis1_dT = dvis1_dtrd1*dtrd1_dT + &
               dvis1_dome1*dome1_dT
  
    ! v1 wrt t
    ! v1 = .407d0*t +80.4d0
    dv1_dT = 0.407d0
  
    ! vs wrt v1, d, t
    if (t .le.350.d0) then
      ! vs = 1.d-7*(v1-d*(1858.d0-5.9d0*t )*1.d-3)
      dvs_dT = 1.d-7*(dv1_dT-d*(-5.9d0)*1.d-3)
      dvs_dd = 1.d-7*(-1.d0*(1858.d0-5.9d0*t )*1.d-3)
      dvs_dRhocomp = dvs_dd*dd_dRhocomp
    else
      ! vs = 1.d-7*(v1 + (0.353d0 + (676.5d-6 + 102.1d-9*d)*d)*d)
      dvs_dT = 1.d-7*(dv1_dT)
      dvs_dd = 1.d-7*(0.353d0 + 2.d0*676.5d-6*d + 3.d0*102.1d-9*d*d)
      dvs_dRhocomp = dvs_dd*dd_dRhocomp
    endif  
  
    ! vis2 wrt vs
    ! vis2 = 10.d0*vs
    dvis2_dvs = 10.d0
    dvis2_dT = dvis2_dvs*dvs_dT
    dvis2_dRhocomp = dvis2_dvs*dvs_dRhocomp
        
    ! vis3 wrt trd3, dome3
    ! vis3 = 266.93d-7*sqrt(fmw3*trd3*fmix)/(cmix*cmix*ome3*trd3)
    tempreal = fmw3*trd3*fmix
    dvis3_dtrd3 = (0.5d0*vis3 - vis3)*one_over_trd3
    dvis3_dome3 = -1.d0*vis3/ome3      
    dvis3_dT = dvis3_dtrd3*dtrd3_dT + dvis3_dome3*dome3_dT
        
    one_over_vis1sq = 1.d0/(vis1*vis1)
    one_over_vis2sq = 1.d0/(vis2*vis2)
    ! z1 wrt vis#,xg#
    ! z1 = xga*xga/vis1+2.d0*xg1*xga/vis3+xg1*xg1/vis2
    dz1_dxga = 2.d0*xga/vis1 + 2.d0*xg1/vis3
    dz1_dxg1 = 2.d0*xga/vis3 + 2.d0*xg1/vis2
    dz1_dvis1 = -1.d0*xga*xga*one_over_vis1sq
    dz1_dvis2 = -1.d0*xg1*xg1*one_over_vis2sq
    dz1_dvis3 = -2.d0*xg1*xga/(vis3*vis3)
    dz1_dT = dz1_dvis1*dvis1_dT + dz1_dvis2*dvis2_dT + dz1_dvis3*dvis3_dT
    dz1_dRhocomp = dz1_dvis2*dvis2_dRhocomp
    dz1_dPcomp = dz1_dxga*dxga_dPcomp + dz1_dxg1*dxg1_dPcomp
    dz1_dPgas = dz1_dxga*dxga_dPgas + dz1_dxg1*dxg1_dPgas
  
    ! g wrt xga
    ! g = xga*xga*FMWAIR/FMWH2O
    dg_dxga = 2.d0*g/xga
    dg_dPcomp = dg_dxga*dxga_dPcomp
    dg_dPgas = dg_dxga*dxga_dPgas

    ! h wrt xg1
    ! h = xg1*xg1*FMWH2O/FMWAIR
    dh_dxg1 = 2.d0*h/xg1
    dh_dPcomp = dh_dxg1*dxg1_dPcomp
    dh_dPgas = dh_dxg1*dxg1_dPgas
  
    ! e wrt xg#, vis#
    ! e = (2.d0*xga*xg1*FMWAIR*FMWH2O/fmw3**2.d0)*vis3/(vis1*vis2)
    de_dxga = e/xga
    de_dxg1 = e/xg1
    de_dvis1 = -1.d0*e/vis1
    de_dvis2 = -1.d0*e/vis2
    de_dvis3 = e/vis3
    de_dT = de_dvis1*dvis1_dT + de_dvis2*dvis2_dT + de_dvis3*dvis3_dT
    de_dRhocomp = de_dvis2*dvis2_dRhocomp
    de_dPcomp = de_dxga*dxga_dPcomp + de_dxg1*dxg1_dPcomp
    de_dPgas = de_dxga*dxga_dPgas + de_dxg1*dxg1_dPgas
      
    ard_0point6 = 0.6d0*ard
    ! z2 wrt ard,e,g,h,vis#,xg#
    ! z2 = 0.6d0*ard*(g/vis1+e+h/vis2)
    dz2_dard = 0.6d0*(g/vis1+e+h/vis2)
    dz2_dg = ard_0point6*(1.d0/vis1)
    dz2_dvis1 = -ard_0point6*(g*one_over_vis1sq)
    dz2_de = ard_0point6
    dz2_dh = ard_0point6*(1.d0/vis2)
    dz2_dvis2 = -ard_0point6*(h*one_over_vis2sq)
    dz2_dT = dz2_dard*dard_dT + dz2_dvis1*dvis1_dT + &
             dz2_de*de_dT + dz2_dvis2*dvis2_dT
    dz2_dRhocomp = dz2_de*de_dRhocomp + dz2_dvis2*dvis2_dRhocomp
    dz2_dPcomp = dz2_dg*dg_dPcomp + dz2_de*de_dPcomp + dz2_dh*dh_dPcomp
    dz2_dPgas = dz2_dg*dg_dPgas + dz2_de*de_dPgas + dz2_dh*dh_dPgas
  
    ! z3 wrt ard,e,g,h,vis#,xg#
    ! z3 = 0.6d0*ard*(g+e*(vis1+vis2)-2.d0*xga*xg1+h)
    dz3_dard = 0.6d0*(g+e*(vis1+vis2)-2.d0*xga*xg1+h)
    dz3_dg = ard_0point6
    dz3_de = ard_0point6*((vis1+vis2))
    dz3_dvis1 = ard_0point6*(e)
    dz3_dvis2 = ard_0point6*(e)
    dz3_dxga = ard_0point6*(-2.d0*xg1)
    dz3_dxg1 = ard_0point6*(-2.d0*xga)
    dz3_dh = ard_0point6
        
    dz3_dT = dz3_dard*dard_dT + dz3_de*de_dT + dz3_dvis1*dvis1_dT + dz3_dvis2*dvis2_dT
    dz3_dRhocomp = dz3_de*de_dRhocomp + dz3_dvis2*dvis2_dRhocomp
    dz3_dPcomp = dz3_dxga*dxga_dPcomp + dz3_dxg1*dxg1_dPcomp + dz3_dg*dg_dPcomp + dz3_de*de_dPcomp + dz3_dh*dh_dPcomp
    dz3_dPgas = dz3_dxga*dxga_dPgas + dz3_dxg1*dxg1_dPgas + dz3_dg*dg_dPgas + dz3_de*de_dPgas + dz3_dh*dh_dPgas
  
    ! visg wrt z#
    ! visg  = (1.d0+z3)/(z1+z2)*.1d0 
    z1pz2 = z1+z2
    dvisg_dz3 = 1.d0/z1pz2*.1d0 
    dvisg_dz1 = -1.d0*(1.d0+z3)/(z1pz2*z1pz2)*.1d0 
    dvisg_dz2 = -1.d0*(1.d0+z3)/(z1pz2*z1pz2)*.1d0 
    
    dV_dT = dvisg_dz1*dz1_dT + dvisg_dz2*dz2_dT + dvisg_dz3*dz3_dT
    dV_dRhocomp = dvisg_dz1*dz1_dRhocomp + dvisg_dz2*dz2_dRhocomp + dvisg_dz3*dz3_dRhocomp
    dV_dPcomp = dvisg_dz1*dz1_dPcomp + dvisg_dz2*dz2_dPcomp + dvisg_dz3*dz3_dPcomp
    dV_dPgas = dvisg_dz1*dz1_dPgas + dvisg_dz2*dz2_dPgas + dvisg_dz3*dz3_dPgas
  endif
  
end subroutine EOSGasViscosity1

! ************************************************************************** !

subroutine EOSGasViscosityConstant(T, P_comp, P_gas, Rho_comp, V_mix, &
                                   calculate_derivative, dV_dT, dV_dPcomp, &
                                   dV_dPgas, dV_dRhocomp, ierr, table_idxs)
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P_comp   ! air pressure [Pa]
  PetscReal, intent(in) :: P_gas    ! gas pressure [Pa]
  PetscReal, intent(in) :: Rho_comp ! air density [C]
  PetscReal, intent(out) :: V_mix   ! mixture viscosity
  PetscBool, intent(in) :: calculate_derivative
  PetscReal, intent(out) :: dV_dT       ! derivative wrt temperature
  PetscReal, intent(out) :: dV_dPcomp   ! derivative wrt component pressure
  PetscReal, intent(out) :: dV_dPgas    ! derivative wrt gas pressure
  PetscReal, intent(out) :: dV_dRhocomp ! derivative wrt component density    
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  V_mix = constant_viscosity
  
  dV_dT = 0.d0
  dV_dPcomp = 0.d0
  dV_dPgas = 0.d0
  dV_dRhocomp = 0.d0
  
end subroutine EOSGasViscosityConstant

! ************************************************************************** !

subroutine EOSGasViscosityEOSDBase(T, P_comp, P_gas, Rho_comp, V_mix, &
                                   calculate_derivative, dV_dT, dV_dPcomp, &
                                   dV_dPgas, dV_dRhocomp, ierr, table_idxs)
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P_comp   ! air pressure [Pa]
  PetscReal, intent(in) :: P_gas    ! gas pressure [Pa]
  PetscReal, intent(in) :: Rho_comp ! air density [C]
  PetscReal, intent(out) :: V_mix   ! mixture viscosity [Pa-s]
  PetscBool, intent(in) :: calculate_derivative
  PetscReal, intent(out) :: dV_dT       ! derivative wrt db temperature [Pa-s/C]
  PetscReal, intent(out) :: dV_dPcomp   ! derivative wrt component pressure [Pa-s/Pa]
  PetscReal, intent(out) :: dV_dPgas    ! derivative wrt gas pressure [Pa-s/Pa]
  PetscReal, intent(out) :: dV_dRhocomp ! derivative wrt component density [Pa-s/kmol/m3]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: NaN

  !ierr initialised in EOSEOSProp 
  !call eos_dbase%EOSProp(T,P_gas,EOS_VISCOSITY,V_mix,ierr)
  call eos_dbase%EOSPropGrad(T,P_gas,EOS_VISCOSITY,V_mix,dV_dT,dV_dPgas,ierr)

  dV_dPcomp = 0.0d0
  dV_dRhocomp = 0.0d0

  ! initialize to derivative to NaN so that not mistakenly used.
  ! NaN = 0.0
  ! NaN = 1.0/NaN
  ! NaN = 0.0d0*NaN
  ! 
  ! dV_dT = NaN
  ! dV_dPcomp = NaN
  ! dV_dPgas = NaN
  ! dV_dRhocomp = NaN
  ! 
  ! if (calculate_derivative) then
  !   ! not yet implemented
  !   ierr = 99 !error 99 points out that deriv are asked but not available yet. 
  !   print*, "EOSGasViscosityEOSDBase - Viscosity derivatives not supported"
  !   stop
  ! end if
  
end subroutine EOSGasViscosityEOSDBase

! ************************************************************************** !

subroutine EOSGasViscosityTable(T, P_comp, P_gas, Rho_comp, V_mix, &
                                   calculate_derivative, dV_dT, dV_dPcomp, &
                                   dV_dPgas, dV_dRhocomp, ierr, table_idxs)
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P_comp   ! air pressure [Pa]
  PetscReal, intent(in) :: P_gas    ! gas pressure [Pa]
  PetscReal, intent(in) :: Rho_comp ! air density [C]
  PetscReal, intent(out) :: V_mix   ! mixture viscosity [Pa-s]
  PetscBool, intent(in) :: calculate_derivative
  PetscReal, intent(out) :: dV_dT       ! derivative wrt db temperature [Pa-s/C]
  PetscReal, intent(out) :: dV_dPcomp   ! derivative wrt component pressure [Pa-s/Pa]
  PetscReal, intent(out) :: dV_dPgas    ! derivative wrt gas pressure [Pa-s/Pa]
  PetscReal, intent(out) :: dV_dRhocomp ! derivative wrt component density [Pa-s/kmol/m3]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: NaN

  !call pvt_table%EOSProp(T,P_gas,EOS_VISCOSITY,V_mix,table_idxs,ierr)
  
  call pvt_table%EOSPropGrad(T,P_gas,EOS_VISCOSITY,V_mix,dV_dT,dV_dPgas, &
                             table_idxs,ierr)

  dV_dPcomp = 0.0d0
  dV_dRhocomp = 0.0d0
  ! initialize to derivative to NaN so that not mistakenly used.
  ! NaN = 0.0
  ! NaN = 1.0/NaN
  ! NaN = 0.0d0*NaN
  ! 
  ! dV_dT = NaN
  ! dV_dPcomp = NaN
  ! dV_dPgas = NaN
  ! dV_dRhocomp = NaN
  ! 
  ! if (calculate_derivative) then
  !   ! not yet implemented
  !   ierr = 99 !error 99 points out that deriv are asked but not available yet.
  !   print*, "EOSGasViscosityTable - Viscosity derivatives not supported"
  !   stop
  ! end if

end subroutine EOSGasViscosityTable

! ************************************************************************** !

subroutine EOSGasDensityNoDerive(T,P,Rho_gas,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
  
  PetscReal :: dum1, dum2
  
  ! derivatives are so cheap, just compute them
  call EOSGasDensityPtr(T,P,Rho_gas,dum1,dum2,ierr,table_idxs)
  
end subroutine EOSGasDensityNoDerive

! ************************************************************************** !

subroutine EOSGasDensityDerive(T,P,Rho_gas,dRho_dT,dRho_dP,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
  
  call EOSGasDensityPtr(T,P,Rho_gas,dRho_dT,dRho_dP,ierr,table_idxs)
                          
end subroutine EOSGasDensityDerive

! ************************************************************************** !

subroutine EOSGasEnergyNoDerive(T,P,H,U,ierr)

  implicit none
  
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dum1, dum2, dum3, dum4
  
  call EOSGasEnergyPtr(T,P,H,dum1,dum2,U,dum3,dum4,ierr)
  
end subroutine EOSGasEnergyNoDerive

! ************************************************************************** !

subroutine EOSGasEnergyDerive(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
    
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr
  
  call EOSGasEnergyPtr(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
  
end subroutine EOSGasEnergyDerive

! ************************************************************************** !

subroutine EOSGasDenEnthNoDerive(T,P,Rho_gas,H,U,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
  
  PetscReal :: dum1, dum2, dum3, dum4, dum5, dum6
  
  call EOSGasDensityEnergyPtr(T,P,Rho_gas,dum1,dum2, &
                              H,dum3,dum4,U,dum5,dum6,ierr,table_idxs)
  
end subroutine EOSGasDenEnthNoDerive

! ************************************************************************** !

subroutine EOSGasDenEnthDerive(T,P,Rho_gas,dRho_dT,dRho_dP, &
                               H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressuret
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
  
  call EOSGasDensityEnergyPtr(T,P,Rho_gas,dRho_dT,dRho_dP, &
                              H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr,table_idxs)
  
end subroutine EOSGasDenEnthDerive

! ************************************************************************** !

subroutine EOSGasDensityEnergyGeneral(T,P,Rho_gas,dRho_dT,dRho_dP, &
                                      H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr,&
                                      table_idxs)
                                      
  implicit none
  
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressuret
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
  
  call EOSGasDensityPtr(T,P,Rho_gas,dRho_dT,dRho_dP,ierr,table_idxs)
  call EOSGasEnergyPtr(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)  
  
end subroutine EOSGasDensityEnergyGeneral

! ************************************************************************** !

subroutine EOSGasDensityIdeal(T,P,Rho_gas,dRho_dT,dRho_dP,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal  T_kelvin

  T_kelvin = T + 273.15d0
  Rho_gas = P / T_kelvin / IDEAL_GAS_CONSTANT * 1.d-3 ! mol/m^3 -> kmol/m^3

  dRho_dP =  Rho_gas / P
  dRho_dT = -Rho_gas / T_kelvin

end subroutine EOSGasDensityIdeal

! ************************************************************************** !

subroutine EOSGasDensityRKS(T,P,Rho_gas,dRho_dT,dRho_dP,ierr,table_idxs)
! Redlich-Kwong-Soave (RKS) equation of state used in BRAGFLO.  See
! Soave, Giorgio, 1972, "Equilibrium constants from a modified Redlich-Kwong
! equation of state", Chem. Eng. Sci., V27, pp 1197-1203.
! current version is for hydrogen only.
!
! Author: Heeho Park
! Date: 05/08/14
!
  use PFLOTRAN_Constants_module, only : UNINITIALIZED_DOUBLE

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
  
  PetscReal :: T_kelvin, RT, alpha
  PetscReal :: Tc_eff, Pc_eff
  PetscReal :: coeff_alpha
  
  ! Newton method approach
  PetscReal :: a_Newton, b_Newton , a_RT, p_RT
  PetscReal :: b2, V, f, dfdV, dVd
  PetscInt :: i
  PetscReal :: coef(4)
  PetscInt, parameter :: maxit = 50

  PetscReal, parameter :: PMIN = 0.1d0  ! from BRAGFLO line 6580
  PetscReal :: Pg

  ! analytical cubic root approach
  PetscReal :: A_cubic, B_cubic
  PetscReal :: c1, c0, xn, yn, del2, del1, h, theta, z1

  !for hydrogen
  ! American Petroleum Institute,
  ! "Technical Data Book-Petroleum Refining" (1977)
  
  dRho_dT = UNINITIALIZED_DOUBLE
  dRho_dP = UNINITIALIZED_DOUBLE
  
  T_kelvin = T + 273.15d0
  RT = IDEAL_GAS_CONSTANT * T_kelvin

  if (rks_use_effect_critical_props) then
    ! Effective critical hydrogen temperature and pressure from BRAGFLO
    ! Taken from BRAGFLO DENGZ()
    Tc_eff = rks_Tc/(1.d0+21.8d-3/(T_kelvin*2.01588d-3))
    Pc_eff = rks_Pc/(1.d0+44.2d-3/(T_kelvin*2.01588d-3))
  else
    Tc_eff = rks_Tc
    Pc_eff = rks_Pc
  endif
  
  if (rks_use_hydrogen) then
  ! geh: commenting out Peter's suggestion that Heeho implemented in favor of
  !      what is implemented in BRAGFLO.
#if 0
    ! suggested by Peter Lichtner
    ! this function honors alpha(Tc)=1 while closely fitting Graboski curve from
    ! "A Modified Soave Equation of State for Phase Equilibrium
    !  Calculations. 3. Systems Containing Hydrogen"
    alpha = exp(0.340d0*(1.d0-T_kelvin/Tc_eff))
#else
    ! from BRAGFLO DENGZ()
    alpha = 1.202d0*exp(-0.30288d0*(T_kelvin/Tc_eff))
#endif
  else
    coeff_alpha = 0.48508d0 + rks_acentric*(1.55171d0 - 0.15613*rks_acentric)
    alpha = (1.d0 + coeff_alpha*(1.d0 - sqrt(T_Kelvin/Tc_eff)))**2.d0
  end if

  if (rks_use_cubic_root_solution) then ! analtyical cubic root soltion

!    Pg = P
    Pg = max(PMIN,P)
  
    A_cubic = rks_coeff_a*alpha*(Pg/Pc_eff)/(T_kelvin/Tc_eff)**2
    B_cubic = rks_coeff_b*(Pg/Pc_eff)/(T_kelvin/Tc_eff)

    ! Solve the cubic eos for compressibility z = pv/(RT)
    ! z^3 - z^2 + (A_cubic-B_cubic-B_cubic^2)*z - (A_cubic*B_cubic) = 0
    ! cubic(a,b,c,d) = (1, -1, c1, c0)
    c1 = A_cubic-B_cubic-B_cubic*B_cubic
    c0 = -A_cubic*B_cubic

    ! Analytical solution to cubic root
    ! See Abramowitz and Stegun pg17 (requires imaginary numbers)
    ! See Cubic function from Wikipedia, and 
    ! Nickalls 1993,  A new approach to solving the cubic Cardans solution 
    ! revealed you get three-real-roots if yn^2 < hn^2, which is the case here
    ! the first root, z1, is the correct one for gas

    xn = 1.0d0/3.0d0
    yn = -2.0d0/27.0d0 + c1/3.0d0 + c0
    del2 = (1.0d0 - 3.0d0*c1)/9.0d0
    del1 = sqrt(del2)
    h = 2.0d0 * del2 * del1

    theta = acos(-yn/h)/3.0d0
    z1 = xn + 2.0d0*del1*cos(theta)
    ! z2 = xn + 2.0d0*del1*cos(theta + 4.0d0/3.0d0*pi)
    ! z3 = xn + 2.0d0*del1*cos(theta + 2.0d0/3.0d0*pi)


    Rho_gas = Pg/(z1*RT) * 1.0d-3 ! mol/m^3 -> kmol/m^3
    ! BRAGFLO: if pressure is less than zero, 
    ! the density is set to 1/10 the density at 1e4 Pa
    Rho_gas = max(Rho_gas, 1.d0/RT) ! kmol/m^3

    if (Pg < 0.d0) then
      print *
      print *, 'EOSGasDensityRKS received a bad pressure: ', Pg
      print *, 'EOSGasDensityRKS temperature: ', T
      print *, 'EOSGasDensityRKS calculated a density: ', Pg/(z1*RT)*1.0d-3
      print *, 'EOSGasDensityRKS reset density to: ', Rho_gas
      print *
    endif

  else ! Newton's method

    Pg = max(PMIN,P)

    a_Newton = rks_coeff_a * alpha * (IDEAL_GAS_CONSTANT * Tc_eff)**2.d0 / &
               Pc_eff
    b_Newton = rks_coeff_b * IDEAL_GAS_CONSTANT * Tc_eff / Pc_eff
  
    a_RT = a_Newton / RT
    P_RT = Pg / RT
    coef(4) = Pg / RT
    coef(3) = 1.0d0
    coef(2) = a_RT - b_Newton - P_RT*b_Newton**2.d0
    coef(1) = a_RT*b_Newton
    V = RT/Pg  ! initial guess
  
    !Newton iteration to find a Volume of gas
    do i = 1, maxit
      f = V*(V*(coef(4)*V - coef(3)) + coef(2)) - coef(1)
      dfdV = V*(3.0d0*coef(4)*V - 2.0d0*coef(3)) + coef(2)
      if (.not. Equal(dfdV,0.d0)) then
        dVd = F/dfdV
        V = V - dVd
        if (.not. Equal(V,0.d0)) then
          if (abs(dVd/V) .lt. 1.d-10) then
            Rho_gas = 1.d0/V * 1.d-3 ! mol/m^3 -> kmol/m^3
            exit
          end if
        else
          print *, 'Error: RKS gas density cannot be solved'
        end if 
      else
        print *, 'Error: Zero Slope in Equation of State'
      end if
    end do

    if (Rho_gas < 0.d0) then
      print *
      print *, 'EOSGasDensityRKS calculated a negative density: ', Rho_gas
      print *, 'Stopping....'
      print *
      stop
    endif

  endif

end subroutine EOSGasDensityRKS

! ************************************************************************** !

subroutine EOSGasDensityPRMethane(T,P,Rho_gas,dRho_dT,dRho_dP,ierr,table_idxs)
! Peng-Robinson (PR) equation of state. This is the uncorrected version 
! sufficient for methane. For other higher hydrocarbons one might need
! the corrected version published in 1978.
! Reference: Peng, D. Y., and Robinson, D. B. (1976).
! A New Two-Constant Equation of State Industrial and Engineering Chemistry: 
! Fundamentals 15: 59–64 DOI:10.1021/i160057a011
! current version is for methane only.
! if omega is greater than 0.49 then alpha needs to be modified below
! Satish Karra, LANL
! 09/30/2014

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
  
  PetscReal :: T_kelvin, RT, alpha,  T_r, Z
  PetscReal :: f, dfdZ, dZd
  PetscReal :: a, b
  PetscInt :: i
  PetscReal, parameter :: coeff_a = 0.45724
  PetscReal, parameter :: coeff_b = 0.0778
  
  ! for methane critical temp, pressure 
  PetscReal, parameter :: Tc = 191.15    ! in K
  PetscReal, parameter :: Pc = 4.641e6   ! in Pa
  PetscReal, parameter :: omega = 0.0115 ! acentric factor 
  
  !solver
  PetscReal :: coef(4)
  PetscInt, parameter :: maxit = 50
  
  T_kelvin = T + 273.15d0
  RT = IDEAL_GAS_CONSTANT*T_kelvin
  T_r = T_kelvin/Tc
  alpha = (1 + (0.3744 + 1.5422*omega - 0.26992*omega**2)*(1-T_r**(0.5)))**2
  
  a = coeff_a * alpha * (IDEAL_GAS_CONSTANT * Tc)**2 / Pc
  b = coeff_b * IDEAL_GAS_CONSTANT * Tc / Pc

  A = alpha * a / (RT)**2
  B = b * P / RT 

  coef(1) = 1.0d0
  coef(2) = B - 1.0d0
  coef(3) = A - 2 * B - 3 * B**2
  coef(4) = B**3 +  B**2 - A * B 
  Z = 1.0d0  ! initial guess for compressibility
  
  !Newton iteration to find a Volume of gas
  do i = 1, maxit
    f = coef(1) * Z**3 + coef(2) * Z**2 + coef(3) * Z + coef(4)
    dfdZ = coef(1) * 3 * Z**2  + coef(2) * 2 * Z + coef(3)
    if (.not. Equal(dfdZ,0.d0)) then
      dZd = f/dfdZ
      Z = Z - dZd
      if (.not. Equal(Z,0.d0)) then
        if (abs(dZd/Z) .lt. 1.d-10) then
          Rho_gas = 1.d0/Z * P/RT * 1d-3 ! mol/m^3 -> kmol/m^3
          exit
        end if
      else
        print *, 'Error: PR76 gas density cannot be solved'
      end if 
    else
      print *, 'Error: Zero Slope in PR76 Equation of State'
    end if
  end do
  
end subroutine EOSGasDensityPRMethane

! ************************************************************************** !

subroutine EOSGasDensityEOSDBase(T, P, Rho_gas, dRho_dT, dRho_dP, ierr, &
                                 table_idxs)
  !
  ! Author: Paolo Orsini
  ! Date: 07/27/17
  ! 
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt table temperature [kmol/m^3/Pa]
  PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt table pressure [kmol/m^3/C]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: NaN

  !ierr initialised in EOSEOSProp 
  !call eos_dbase%EOSProp(T,P,EOS_DENSITY,Rho_gas,ierr)
  call eos_dbase%EOSPropGrad(T,P,EOS_DENSITY,Rho_gas,dRho_dT,dRho_dP,ierr)

  !PO todo: conversion when loaidng database to do this operation only once
           ! conversion to molar density
           ! kg/m3 * kmol/kg  = kmol/m3
  Rho_gas = Rho_gas / fmw_gas ! kmol/m^3
  dRho_dT = dRho_dT / fmw_gas
  dRho_dP = dRho_dP / fmw_gas

  ! ! initialize to derivative to NaN so that not mistakenly used.
  ! NaN = 0.0
  ! NaN = 1.0/NaN
  ! NaN = 0.0d0*NaN
  ! dRho_dT = NaN
  ! dRho_dP = NaN

end subroutine EOSGasDensityEOSDBase

! ************************************************************************** !

subroutine EOSGasDensityTable(T, P, Rho_gas, dRho_dT, dRho_dP, ierr, &
                                 table_idxs)
  !
  ! Author: Paolo Orsini
  ! Date: 10/30/17
  !

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! oil density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt table temperature [kmol/m^3/C]
  PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt table pressure [kmol/m^3/Pa]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: NaN

  !Rho from pvt table is in kmol/m3
  !call pvt_table%EOSProp(T,P,EOS_DENSITY,Rho_gas,table_idxs,ierr)

  call pvt_table%EOSPropGrad(T,P,EOS_DENSITY,Rho_gas,dRho_dT,dRho_dP, &
                             table_idxs,ierr)

  ! initialize to derivative to NaN so that not mistakenly used.
  ! NaN = 0.0
  ! NaN = 1.0/NaN
  ! NaN = 0.0d0*NaN
  ! dRho_dT = NaN
  ! dRho_dP = NaN

end subroutine EOSGasDensityTable

! ************************************************************************** !

subroutine EOSGasFugacity(T,P,Rho_gas,phi,ierr)
! current version is for hydrogen only. See
! Soave, Giorgio, 1972, "Equilibrium constants from a modified Redlich-Kwong
! equation of state", Chem. Eng. Sci., V27, pp 1197-1203.

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(in) :: Rho_gas  ! gas density [kmol/m^3]
  PetscReal, intent(out) :: phi     ! fugacity coefficient
  PetscErrorCode, intent(out) :: ierr

  !for hydrogen
  ! American Petroleum Institute,
  ! "Technical Data Book-Petroleum Refining" (1977)
  PetscReal, parameter :: Tc = 41.67d0   ! critical temperature
  PetscReal, parameter :: Pc = 2.1029d6  ! critical pressure

  PetscReal :: T_kelvin, Z, V, Tr, LHS
  PetscReal :: A, B, alpha

  V = 1.d0 / Rho_gas * 1.d-3
  T_kelvin = T + 273.15d0
  Tr = T_kelvin/Tc  ! reduced temperature
  alpha = EXP(0.340d0*(1.d0-Tr))  ! alpha(T)

  A = 0.42747 * alpha * (P/Pc) / (T/Tc)**2
  B = 0.08664 * (P/Pc) / (T/Tc)
  Z = P*V / (IDEAL_GAS_CONSTANT*T)
  
  ! Fugacity Coefficient
  LHS = Z - 1.d0 - LOG(Z-B) - A/B * LOG((Z+B)/Z)
  phi = EXP(LHS)  ! dimensionless (Pa/Pa)
  
end subroutine EOSGasFugacity

! ************************************************************************** !

subroutine EOSGasEnergyIdealMethane(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
  
! Based on EOSGasEnergyIdeal. Modified for Methane
! Satish Karra, LANL
! 09/30/2014
  
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr

  ! Cpg units: J/mol-K
  PetscReal, parameter :: Cp_methane = 35.69 ! at 298.15 K and 1 bar http://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Units=SI&Mask=1
  PetscReal :: T_energy
  PetscReal :: T_k

  ierr = 0
  ! T_energy is either T or T + 273.15
  ! do not change below
  T_energy = T + T_energy_offset
  T_k = T + 273.15d0
  H = Cp_methane * T_energy * 1.d3  ! J/mol -> J/kmol
  U = H - IDEAL_GAS_CONSTANT * T_k * 1.d3 ! J/mol -> J/kmol

  dH_dP = 0.d0
  dH_dT = Cp_methane * 1.d3
  dU_dP = 0.d0
  dU_dT = (Cp_methane - IDEAL_GAS_CONSTANT) * 1.d3
    
end subroutine EOSGasEnergyIdealMethane

! ************************************************************************** !

subroutine EOSGasEnergyIdeal(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
    
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr

  ! Cv_air units: J/mol-K
  PetscReal, parameter :: Cv_air = 20.85 ! heat capacity wiki [J/mol-K]
  PetscReal :: T_energy
  PetscReal :: T_k

  ! T_energy is either T or T + 273.15
  ! do not change below
  T_energy = T + T_energy_offset
  T_k = T + 273.15d0
  U = Cv_air * T_energy * 1.d3  ! J/mol -> J/kmol
  H = U + IDEAL_GAS_CONSTANT * T_k * 1.d3 ! J/mol -> J/kmol

  dU_dP = 0.d0
  dU_dT = Cv_air * 1.d3
  dH_dP = 0.d0
  dH_dT = dU_dT + IDEAL_GAS_CONSTANT * 1.d3
    
end subroutine EOSGasEnergyIdeal

! ************************************************************************** !

subroutine EOSGasEnergyEOSDBase(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
    
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature [J/kmol/C]
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure [J/kmol/Pa]
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature [J/kmol/C]
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure [J/kmol/Pa]
  PetscErrorCode, intent(out) :: ierr

  PetscReal :: NaN

  !ierr initialised in EOSProp - PO: should do only one lookup here
  !call eos_dbase%EOSProp(T,P,EOS_ENTHALPY,H,ierr)
  !call eos_dbase%EOSProp(T,P,EOS_INTERNAL_ENERGY,U,ierr)

  call eos_dbase%EOSPropGrad(T,P,EOS_ENTHALPY,H,dH_dT,dH_dP,ierr)
  call eos_dbase%EOSPropGrad(T,P,EOS_INTERNAL_ENERGY,U,dU_dT,dU_dP,ierr)
  
  !PO todo: conversion when loaidng database to do this operation only once
  ! conversion to molar energy
  ! J/kg * kg/Kmol = J/Kmol
  H = H  * fmw_gas
  dH_dT = dH_dT * fmw_gas
  dH_dP = dH_dP * fmw_gas
  
  U = U  * fmw_gas
  dU_dT = dU_dT * fmw_gas
  dU_dP = dU_dP * fmw_gas

  ! NaN = 0.0
  ! NaN = 1.0/NaN
  ! NaN = 0.0d0*NaN
  ! 
  ! dU_dP = NaN
  ! dU_dT = NaN
  ! dH_dP = NaN
  ! dH_dT = NaN
    
end subroutine EOSGasEnergyEOSDBase

! ************************************************************************** !

subroutine EOSGasDensityConstant(T,P,Rho_gas,dRho_dT,dRho_dP,ierr,table_idxs)

  implicit none
  
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
  
  ierr = 0
  Rho_gas = constant_density ! kmol/m^3

  dRho_dT = 0.d0
  dRho_dP = 0.d0

end subroutine EOSGasDensityConstant

! ************************************************************************** !

subroutine EOSGasEnergyConstant(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
    
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: T_energy

  ierr = 0
  H = constant_enthalpy ! J/kmol
  ! T_energy is either T or T + 273.15
  ! do not change below
  T_energy = T + T_energy_offset
  U = H - IDEAL_GAS_CONSTANT * T_energy * 1.d3 ! J/mol -> J/kmol
  
  dH_dP = 0.d0
  dH_dT = 0.d0
  dU_dP = 0.d0
  dU_dT = -1.d0 * IDEAL_GAS_CONSTANT * 1.d3
  
  print *, 'Calculation of internal gas energy not set up properly in ' // &
    'EOSGasEnergyConstant.'
  stop  
  
end subroutine EOSGasEnergyConstant

! ************************************************************************** !

subroutine EOSGasHenryNoDerive(T,Psat,Hc)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: Psat     ! saturation pressure
  PetscReal, intent(out) :: Hc      ! Henry's constant
  
  PetscReal :: dum1, dum2, dum3, dum4
  
  call EOSGasHenryPtr(T,Psat,Hc,PETSC_FALSE,dum1,dum2,dum3,dum4)
                          
end subroutine EOSGasHenryNoDerive

! ************************************************************************** !

subroutine EOSGasHenryDerive(T,Psat,Psat_P,Psat_T,Hc,Hc_P,Hc_T)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: Psat     ! saturation pressure
  PetscReal, intent(in) :: Psat_P   ! derivative Psat wrt pressure
  PetscReal, intent(in) :: Psat_T   ! derivative Psat wrt temperature
  PetscReal, intent(out) :: Hc      ! Henry's constant
  PetscReal, intent(out) :: Hc_P    ! derivative Henry's constant wrt pressure
  PetscReal, intent(out) :: Hc_T    ! derivative Henry's constant wrt temperature
  
!  PetscReal :: pert, T_pert, Hc_pert, dum1, dum2
  
  call EOSGasHenryPtr(T,Psat,Hc,PETSC_TRUE,Psat_P,Psat_T,Hc_P,Hc_T)
!  pert = 1.d-8 * T
!  T_pert = T + pert
!  call EOSGasHenryPtr(T_pert,Psat,Hc_pert,PETSC_FALSE,Psat_P,Psat_T,dum1,dum2)
!  T_pert = (Hc_pert - Hc)/pert
                          
end subroutine EOSGasHenryDerive

#if 0
!TODO(geh): remove if derivative version matching
! ************************************************************************** !

subroutine EOSGasHenry_air_noderiv(tc,ps,Henry)
! 
!   Calculates Henry's constant as a function of temperature [C], 
!   and saturation pressure [Pa].  
!
!   Fernandez-Prini, F., J. Alvarez and A. Harvey (2003) Henry's Constants and
!   Vapor-Liquid Distribution Constants for Gaseous Solutes in H2O and D2O at
!   High Temperatures, J. Phys. Chem. Ref. Data, Vol. 32, No. 2, Equation 15
!   with coefficients A,B,C from Table 3 for N2(g)

    implicit none
    
    PetscReal, intent(in) :: tc
    PetscReal, intent(in) :: ps
    PetscReal, intent(out) :: Henry

    PetscReal  Tr,tao,tmp,t
    PetscReal, parameter :: a=-9.67578d0, b=4.72162d0, c=11.70585d0
    PetscReal, parameter :: Tcl=647.096d0 ! H2O critical temp(K) from IAPWS(1995b)

    t = tc+273.15D0
    Tr = t/Tcl
    tao = 1.D0-Tr
    tmp = a/Tr + b * tao**0.355d0/Tr + c * (Tr**(-0.41d0)) * exp(tao)
    Henry = exp(tmp)*ps

end subroutine EOSGasHenry_air_noderiv
#endif
! ************************************************************************** !

subroutine EOSGasHenry_air(T,Psat,Hc,calculate_derivative, &
                           Psat_P,Psat_T,Hc_P,Hc_T)
! 
!   Calculates Henry's constant as a function of temperature [C], 
!   and saturation pressure [Pa].  
!
!   Fernandez-Prini, F., J. Alvarez and A. Harvey (2003) Henry's Constants and
!   Vapor-Liquid Distribution Constants for Gaseous Solutes in H2O and D2O at
!   High Temperatures, J. Phys. Chem. Ref. Data, Vol. 32, No. 2, Equation 15
!   with coefficients A,B,C from Table 3 for N2(g)

    implicit none
    
    PetscReal, intent(in) :: T        ! temperature [C]
    PetscReal, intent(in) :: Psat     ! saturation pressure
    PetscReal, intent(out) :: Hc      ! Henry's constant
    PetscBool, intent(in) :: calculate_derivative
    PetscReal, intent(in) :: Psat_P   ! derivative Psat wrt pressure
    PetscReal, intent(in) :: Psat_T   ! derivative Psat wrt temperature
    PetscReal, intent(out) :: Hc_P    ! derivative Henry's constant wrt pressure
    PetscReal, intent(out) :: Hc_T    ! derivative Henry's constant wrt temperature

    PetscReal :: Tr,tau,TK
    PetscReal :: tmp, tmpA, tmpB, tmpC
    PetscReal :: dTr_dT, dtau_dT, dtmp_dTr, dtmp_dtau
    PetscReal, parameter :: a=-9.67578d0, b=4.72162d0, c=11.70585d0 ! for N2
    PetscReal, parameter :: Tcw=647.096d0 ! H2O critical temp from IAPWS(1995b)

    TK = T+273.15D0
    Tr = TK/Tcw
    tau = 1.D0-Tr
    !tmp = a/Tr + b * tau**0.355d0/Tr + c * (Tr**(-0.41d0)) * exp(tau)
    tmpA = a/Tr
    tmpB = b*(tau**0.355d0)/Tr
    tmpC = c*(Tr**(-0.41d0))*exp(tau)
    tmp = tmpA+tmpB+tmpC
    Hc = exp(tmp)*Psat
    if (calculate_derivative) then
      dTr_dT = 1.d0/Tcw
      dtau_dT = -1.d0*dTr_dT
      dtmp_dTr = (-tmpA-tmpB-0.41d0*tmpC)/Tr
      dtmp_dtau = 0.355d0*tmpB/tau+tmpC
      Hc_T = Hc*((dtmp_dTr*dTr_dT+dtmp_dtau*dtau_dT)+Psat_T/Psat)
!      tmp = ((-a/Tr+b*(-0.355d0*tau**(-0.645d0)-tau**0.355d0/Tr))/Tr - &
!             c*exp(tau)*(tau**(-0.41d0))*(0.41d0/Tr-1.d0))/Tcw
!      Hc_T = Hc*(tmp+Psat_T/Psat)
      Hc_P = Psat_P*Hc/Psat
    endif

end subroutine EOSGasHenry_air

! ************************************************************************** !

subroutine EOSGasHenryConstant(T,Psat,Hc,calculate_derivative, &
                               Psat_P,Psat_T,Hc_P,Hc_T)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: Psat     ! saturation pressure
  PetscReal, intent(out) :: Hc      ! Henry's constant
  PetscBool, intent(in) :: calculate_derivative
  PetscReal, intent(in) :: Psat_P   ! derivative Psat wrt pressure
  PetscReal, intent(in) :: Psat_T   ! derivative Psat wrt temperature
  PetscReal, intent(out) :: Hc_P    ! derivative Henry's constant wrt pressure
  PetscReal, intent(out) :: Hc_T    ! derivative Henry's constant wrt temperature

  Hc = constant_henry
  Hc_P = 0.d0
  Hc_T = 0.d0
    
end subroutine EOSGasHenryConstant

! **************************************************************************** !

subroutine EOSGasSetEOSDBase(filename,option)

  use Option_module

  implicit none

  character(len=MAXWORDLENGTH) :: filename
  type(option_type) :: option

  eos_dbase => EOSDatabaseCreate(filename,'gas_database')
  call eos_dbase%Read(option)

  !set property function pointers
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasDensityPtr => EOSGasDensityEOSDBase
  EOSGasEnergyPtr => EOSGasEnergyEOSDBase
  EOSGasViscosityPtr => EOSGasViscosityEOSDBase

end subroutine EOSGasSetEOSDBase

! **************************************************************************** !

subroutine EOSGasSetPVDG(input,option)
  !
  ! Author: Paolo Orsini
  ! Date: 10/30/17
  !
  ! Set up a PVDG table

  use Option_module
  use Input_Aux_module
  use Lookup_Table_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  type(lookup_table_var_type), pointer :: db_var => null()
  character(len=MAXWORDLENGTH) :: internal_units, user_units
  PetscInt :: data_idx

  pvt_table => EOSTableCreate('PVDG',option)

  pvt_table%num_prop = 2
  
  ! units initially assing default values - overwritten by units specified 
  ! in the table input 
  internal_units = '' !assign default value by SetDefaultInternalUnits
  user_units = ''     !assign default value by SetMetricUnits
  
  !adding FVF 
  data_idx = 1 !position of FVF in the table (after pressure)
  db_var => CreateLookupTableVar(EOS_FVF,internal_units,user_units,data_idx)
  call pvt_table%AddEOSProp(db_var,option)
  nullify(db_var)

  !adding VISCOSITY 
  data_idx = 2 !position of FVF in the table (after pressure)
  db_var => CreateLookupTableVar(EOS_VISCOSITY,internal_units,user_units, &
                                 data_idx)
  call pvt_table%AddEOSProp(db_var,option)
  nullify(db_var)

  !set Default internal must be called before Set Metric
  call pvt_table%SetDefaultInternalUnits(option)
  call pvt_table%SetMetricUnits(option)

  call pvt_table%Read(input,option)

  call EOSTableAddToList(pvt_table,eos_table_list)

  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasDensityPtr => EOSGasDensityTable
  EOSGasViscosityPtr => EOSGasViscosityTable

end subroutine EOSGasSetPVDG

! ************************************************************************** !

subroutine EOSGasTableProcess(option)
  !
  ! Author: Paolo Orsini
  ! Date: 10/30/17
  !
  ! Processes gas pvt table - once the entire input deck has been read

  use Option_module

  implicit none

  type(option_type) :: option

  if (.not.associated(pvt_table)) return

  select case(pvt_table%name)
  case("PVDG")
      call pvt_table%ConvertFVFtoMolarDensity(fmw_gas,reference_density_kg)
  end select

end subroutine EOSGasTableProcess

! ************************************************************************** !

subroutine EOSGasInputRecord()
  ! 
  ! Prints ingested equation of state information to the input record file.
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/04/2016
  ! 
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: word1, word2
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'GAS'
  
  ! gas density [kg/m^3]
  if (associated(EOSGasDensityEnergyPtr,EOSGasDensityEnergyGeneral) .and. &
      associated(EOSGasDensityPtr,EOSGasDensityConstant)) then
    write(id,'(a29)',advance='no') 'gas density: '
    write(word1,*) constant_density
    write(id,'(a)') 'constant, ' // adjustl(trim(word1)) // ' kg/m^3'
  endif
  if (associated(EOSGasDensityEnergyPtr,EOSGasDensityEnergyGeneral) .and. &
      associated(EOSGasDensityPtr,EOSGasDensityRKS)) then
    write(id,'(a29)',advance='no') 'gas density: '
    write(id,'(a)') 'rsk'
    write(id,'(a29)',advance='no') 'critical temperature: '
    write(word1,*) rks_Tc
    write(id,'(a)') adjustl(trim(word1)) // ' C'
    write(id,'(a29)',advance='no') 'critical pressure: '
    write(word1,*) rks_Pc
    write(id,'(a)') adjustl(trim(word1)) // ' Pa'
    if (rks_use_cubic_root_solution) then
      write(id,'(a29)',advance='no') 'use analytical cubic root solution: True'
    else
      write(id,'(a29)',advance='no') 'use analytical cubic root solution: False'
    endif
    if (.not.rks_use_hydrogen) then
      write(id,'(a29)',advance='no') 'acentric factor: '
      write(word1,*) rks_acentric
      write(id,'(a)') adjustl(trim(word1)) 
    endif
    write(id,'(a29)',advance='no') 'omega A: '
    write(word1,*) rks_coeff_a
    write(id,'(a)') adjustl(trim(word1))
    write(id,'(a29)',advance='no') 'omega B: '
    write(word1,*) rks_coeff_b
    write(id,'(a)') adjustl(trim(word1))
  endif
  if (associated(EOSGasDensityEnergyPtr,EOSGasDensityEnergyGeneral) .and. &
      associated(EOSGasDensityPtr,EOSGasDensityPRMethane)) then
    write(id,'(a29)',advance='no') 'gas density: '
    write(id,'(a)') 'pr methane'
  endif
  if (associated(EOSGasDensityEnergyPtr,EOSGasDensityEnergyGeneral) .and. &
      associated(EOSGasDensityPtr,EOSGasDensityIdeal)) then
    write(id,'(a29)',advance='no') 'gas density: '
    write(id,'(a)') 'default, ideal'
  endif
  
  ! gas viscosity [Pa-s]
  if (associated(EOSGasViscosityPtr,EOSGasViscosityConstant)) then
    write(id,'(a29)',advance='no') 'gas viscosity: '
    write(word1,*) constant_viscosity
    write(id,'(a)') 'constant, ' // trim(word1) // ' Pa-sec'
  endif
  if (associated(EOSGasViscosityPtr,EOSGasViscosity1)) then
    write(id,'(a29)',advance='no') 'gas viscosity: '
    write(id,'(a)') 'default'
  endif
  
  ! gas enthalpy [J/kmol]
  if (associated(EOSGasDensityEnergyPtr,EOSGasDensityEnergyGeneral) .and. &
      associated(EOSGasEnergyPtr,EOSGasEnergyConstant)) then
    write(id,'(a29)',advance='no') 'gas enthalpy: '
    write(word1,*) constant_enthalpy
    write(id,'(a)') 'constant, ' // trim(word1) // 'J/kmol'
  endif
  if (associated(EOSGasDensityEnergyPtr,EOSGasDensityEnergyGeneral) .and. &
      associated(EOSGasEnergyPtr,EOSGasEnergyIdeal)) then
    write(id,'(a29)',advance='no') 'gas enthalpy: '
    write(id,'(a)') 'default, ideal'
  endif
  if (associated(EOSGasDensityEnergyPtr,EOSGasDensityEnergyGeneral) .and. &
      associated(EOSGasEnergyPtr,EOSGasEnergyIdealMethane)) then
    write(id,'(a29)',advance='no') 'gas enthalpy: '
    write(id,'(a)') 'ideal methane'
  endif
  
  ! henry's constant
  if (associated(EOSGasHenryPtr,EOSGasHenryConstant)) then
    write(id,'(a29)',advance='no') "henry's constant: "
    write(word1,*) constant_henry
    write(id,'(a)') 'constant, ' // trim(word1) 
  endif
  if (associated(EOSGasHenryPtr,EOSGasHenry_air)) then
    write(id,'(a29)',advance='no') "henry's constant: "
    write(id,'(a)') 'default, air'
  endif
  
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  
end subroutine EOSGasInputRecord

! ************************************************************************** !

subroutine EOSGasTest(temp_low,temp_high,pres_low,pres_high, &
                        ntemp,npres,uniform_temp,uniform_pres,filename)

  use EOS_Water_module, only : EOSWaterSaturationPressure

  implicit none

  PetscReal :: temp_low
  PetscReal :: temp_high
  PetscReal :: pres_low
  PetscReal :: pres_high
  PetscInt :: npres
  PetscInt :: ntemp
  PetscBool :: uniform_temp
  PetscBool :: uniform_pres
  character(len=MAXWORDLENGTH) :: filename

  PetscReal, allocatable :: temp(:)
  PetscReal, allocatable :: pres(:)
  PetscReal, allocatable :: density_kg(:,:)
  PetscReal, allocatable :: enthalpy(:,:)
  PetscReal, allocatable :: internal_energy(:,:)
  PetscReal, allocatable :: viscosity(:,:)
  PetscReal, allocatable :: saturation_pressure_array(:)
  PetscReal, allocatable :: henry(:)
  PetscReal :: dum1, dum2, dum3, dum4
  PetscInt :: itemp, ipres
  PetscReal :: ln_low, ln_high
  PetscReal :: saturation_pressure
  PetscReal :: air_pressure
  PetscReal :: NaN
  character(len=MAXWORDLENGTH) :: eos_density_name
  character(len=MAXWORDLENGTH) :: eos_energy_name
  character(len=MAXWORDLENGTH) :: eos_viscosity_name
  character(len=MAXWORDLENGTH) :: eos_saturation_pressure_name
  character(len=MAXWORDLENGTH) :: eos_henry_name
  character(len=MAXSTRINGLENGTH) :: header, string

  PetscErrorCode :: ierr

  NaN = 0.d0
  NaN = 1.d0/NaN
  NaN = 0.d0*NaN

  allocate(temp(ntemp))
  temp = UNINITIALIZED_DOUBLE
  allocate(pres(npres))
  pres = UNINITIALIZED_DOUBLE
  allocate(density_kg(npres,ntemp))
  density_kg = UNINITIALIZED_DOUBLE
  allocate(viscosity(npres,ntemp))
  viscosity = UNINITIALIZED_DOUBLE
  allocate(enthalpy(npres,ntemp))
  enthalpy = UNINITIALIZED_DOUBLE
  allocate(internal_energy(npres,ntemp))
  internal_energy = UNINITIALIZED_DOUBLE
  allocate(saturation_pressure_array(ntemp))
  saturation_pressure_array = UNINITIALIZED_DOUBLE
  allocate(henry(ntemp))
  henry = UNINITIALIZED_DOUBLE

  if (uniform_pres) then
    do ipres = 1, npres
      pres(ipres) = &
        (pres_high-pres_low)/max(dble(npres-1),1.d0) * (ipres-1) + pres_low
    enddo
  else
    ln_high = log(pres_high)
    ln_low = log(pres_low)
    do ipres = 1, npres
      pres(ipres) = &
        exp((ln_high-ln_low)/max(dble(npres-1),1.d0) * (ipres-1) + ln_low)
    enddo
  endif

  if (uniform_temp) then
    do itemp = 1, ntemp
      temp(itemp) = &
        (temp_high-temp_low)/max(dble(ntemp-1),1.d0) * (itemp-1) + temp_low
    enddo
  else
    ln_high = log(temp_high)
    ln_low = log(temp_low)
    do itemp = 1, ntemp
      temp(itemp) = &
        exp((ln_high-ln_low)/max(dble(ntemp-1),1.d0) * (itemp-1) + ln_low)
    enddo
  endif

  ! density
  if (associated(EOSGasDensityPtr,EOSGasDensityConstant)) then
    eos_density_name = 'Constant'
  else if (associated(EOSGasDensityPtr,EOSGasDensityIdeal)) then
    eos_density_name = 'Ideal Gas Law'
  else if (associated(EOSGasDensityPtr,EOSGasDensityRKS)) then
    eos_density_name = 'Redlich-Kwong-Soave'
  else if (associated(EOSGasDensityPtr,EOSGasDensityPRMethane)) then
    eos_density_name = 'Peng-Robinson Methane'
  else 
    eos_density_name = 'Unknown'
  endif

  ! energy
  if (associated(EOSGasEnergyPtr,EOSGasEnergyConstant)) then
    eos_energy_name = 'Constant'
  else if (associated(EOSGasEnergyPtr,EOSGasEnergyIdeal)) then
    eos_energy_name = 'Ideal Gas Law'
  else if (associated(EOSGasEnergyPtr,EOSGasEnergyIdealMethane)) then
    eos_energy_name = 'Ideal Gas Law - Methane'
  else
    eos_energy_name = 'Unknown'
  endif

  ! viscosity
  if (associated(EOSGasViscosityPtr,EOSGasViscosityConstant)) then
    eos_viscosity_name = 'Constant'
  else if (associated(EOSGasViscosityPtr,EOSGasViscosity1)) then
    eos_viscosity_name = 'Default'
  else
    eos_viscosity_name = 'Unknown'
  endif

  ! Henry's constant
  if (associated(EOSGasHenryPtr,EOSGasHenryConstant)) then
    eos_henry_name = 'Constant'
  else if (associated(EOSGasHenryPtr,EOSGasHenry_air)) then
    eos_henry_name = 'Default'
  else
    eos_henry_name = 'Unknown'
  endif

  ! saturation pressure
  eos_saturation_pressure_name = 'Unknown'

  do itemp = 1, ntemp
    ! Need saturation pressure to calculate air partial pressure
    call EOSWaterSaturationPressure(temp(itemp),saturation_pressure,ierr)
    saturation_pressure_array(itemp) = saturation_pressure
    call EOSGasHenry(temp(itemp),saturation_pressure,henry(itemp))
    do ipres = 1, npres
      call EOSGasDensityPtr(temp(itemp),pres(ipres), &
                            density_kg(ipres,itemp), &
                            dum1,dum2,ierr)
      call EOSGasEnergyPtr(temp(itemp),pres(ipres),&
                           enthalpy(ipres,itemp),dum1,dum2, &
                           internal_energy(ipres,itemp),dum3,dum4,ierr)
      air_pressure = pres(ipres) - saturation_pressure
      if (air_pressure > 0.d0) then
        call EOSGasViscosityPtr(temp(itemp),air_pressure,pres(ipres), &
                                density_kg(ipres,itemp), &
                                viscosity(ipres,itemp),PETSC_FALSE, &
                                dum1,dum2,dum3,dum4,ierr)
      else
        viscosity(ipres,itemp) = NaN
      endif
    enddo
  enddo

100 format(100es16.8)
  if (len_trim(filename) == 0) then
    string = 'eos_gas_test.txt'
  else
    string = filename
  endif
  open(unit=IUNIT_TEMP,file=string)
  header = 'T[C], P[Pa], &
    &Density (' // trim(eos_density_name) // ') [kmol/m^3], &
    &Enthalpy (' // trim(eos_energy_name) // ') [J/kmol], &
    &Internal Energy (' // trim(eos_energy_name) // ') [J/kmol], &
    &Viscosity (' // trim(eos_viscosity_name) // ') [Pa-s], &
    &Saturation Pressure (' // trim(eos_saturation_pressure_name) // ") [Pa], &
    &Henry's Constant (" // trim(eos_henry_name) // ') [-]'
  write(IUNIT_TEMP,'(a)') trim(header)
  write(IUNIT_TEMP,'(100i9)') ntemp, npres
  do itemp = 1, ntemp
    do ipres = 1, npres
      write(IUNIT_TEMP,100) temp(itemp), pres(ipres), &
            density_kg(ipres,itemp), enthalpy(ipres,itemp), &
            internal_energy(ipres,itemp), viscosity(ipres,itemp), &
            saturation_pressure_array(itemp), henry(itemp)
    enddo
  enddo
  close(IUNIT_TEMP)

  deallocate(temp)
  deallocate(pres)
  deallocate(density_kg)
  deallocate(enthalpy)
  deallocate(internal_energy)
  deallocate(viscosity)
  deallocate(saturation_pressure_array)
  deallocate(henry)

end subroutine EOSGasTest

! ************************************************************************** !

subroutine EOSGasDBaseDestroy()

  implicit none

  !destroy databases
  call EOSDatabaseDestroy(eos_dbase)

  nullify(pvt_table)

end subroutine EOSGasDBaseDestroy

! ************************************************************************** !

end module EOS_Gas_module
