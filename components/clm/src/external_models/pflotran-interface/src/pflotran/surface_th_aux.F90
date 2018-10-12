module Surface_TH_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public :: Surface_TH_auxvar_type
    PetscReal :: h        ! enthalpy -- not currently used
    PetscReal :: u        ! internal energy -- not currently used
    PetscReal :: pc       ! pressure change -- not currently used
    PetscReal :: Cw       ! Specific heat capacity of surface water
    PetscReal :: Ci       ! Specific heat capacity of surface ice
    PetscReal :: Cwi      ! Weighted average of Cw and Ci
      ! RTM: Note that I believe we can simple make Cw and Ci Fortran 
      ! parameters, but I'm keeping things as is for now.
    PetscReal :: k_therm  ! Thermal conductivity of surface water
    PetscReal :: unfrozen_fraction ! Proportion of unfrozen surface water.
    PetscReal :: den_water_kg  ! Density [kg/m^3] of liquid water ONLY.
      ! Note that we currently use the den_kg(1) field of the global aux var type 
      ! to store the density of the liquid/ice mixture.  We also need to track
      ! the liquid density because it is required in in the advective term of 
      ! the temperature equation.  
  end type Surface_TH_auxvar_type

  type, public :: Surface_TH_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(Surface_TH_auxvar_type), pointer :: auxvars(:)
    type(Surface_TH_auxvar_type), pointer :: auxvars_bc(:)
    type(Surface_TH_auxvar_type), pointer :: auxvars_ss(:)
  end type Surface_TH_type

  public :: SurfaceTHAuxCreate, &
            SurfaceTHAuxDestroy, &
            SurfaceTHAuxVarCompute, &
            SurfaceTHAuxVarComputeUnfrozen, &
            SurfaceTHAuxVarInit, &
            SurfaceTHAuxVarCopy

contains

! ************************************************************************** !

function SurfaceTHAuxCreate(option)
  ! 
  ! This routine creates an empty object
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module

  implicit none
  
  type(option_type) :: option
  type(Surface_TH_type), pointer :: SurfaceTHAuxCreate
  
  type(Surface_TH_type), pointer :: aux

  allocate(aux)
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  aux%n_zero_rows = 0
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  SurfaceTHAuxCreate => aux
  
end function SurfaceTHAuxCreate

! ************************************************************************** !

subroutine SurfaceTHAuxVarInit(auxvar,option)
  ! 
  ! This routine initilizes an object.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 

  use Option_module
  use PFLOTRAN_Constants_module, only : UNINITIALIZED_DOUBLE

  implicit none
  
  type(Surface_TH_auxvar_type) :: auxvar
  type(option_type) :: option

  PetscReal :: uninit_value
  uninit_value = UNINITIALIZED_DOUBLE
  
  auxvar%h   = uninit_value
  auxvar%u   = uninit_value
  auxvar%pc  = uninit_value
  auxvar%Cw  = 4.188d3     ! [J/kg/K]
  !auxvar%Ci = 2.050d3     ! [J/kg/K]
  auxvar%Ci  = 4.188d3     ! [J/kg/K]
  auxvar%Cwi = 4.188d3     ! [J/kg/K]
  auxvar%k_therm = 0.57d0 ! [J/s/m/K]
  auxvar%unfrozen_fraction = 1.d0
  auxvar%den_water_kg = 1.d3  ! [kg/m^3]

end subroutine SurfaceTHAuxVarInit

! ************************************************************************** !

subroutine SurfaceTHAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! This routine makes a copy of an object.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 

  use Option_module

  implicit none
  
  type(Surface_TH_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%h = auxvar%h
  auxvar2%u = auxvar%u
  auxvar2%pc = auxvar%pc
  auxvar2%Cw = auxvar%Cw
  auxvar2%Ci = auxvar%Ci
  auxvar2%Cwi = auxvar%Cwi
  auxvar2%k_therm = auxvar%k_therm
  auxvar2%unfrozen_fraction = auxvar%unfrozen_fraction
  auxvar2%den_water_kg = auxvar%den_water_kg

end subroutine SurfaceTHAuxVarCopy

! ************************************************************************** !

subroutine SurfaceTHAuxVarCompute(xx,auxvar,global_auxvar, &
                                  option)
  ! 
  ! This routine computes values for auxiliary variables.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use Surface_Global_Aux_module
  
  use EOS_Water_module  
  use PFLOTRAN_Constants_module, only : DUMMY_VALUE,MIN_SURFACE_WATER_HEIGHT

  implicit none

  type(option_type) :: option
  PetscReal :: xx(option%nflowdof)
  type(Surface_TH_auxvar_type) :: auxvar
  type(surface_global_auxvar_type) :: global_auxvar
  PetscReal :: por, perm

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: di_kg, dwi_kg
    ! Densities of ice and water-ice mixture, respectively.
  PetscReal :: unfrozen_fraction
  PetscReal :: kr, ds_dp, dkr_dp
  PetscReal :: dvis_dt, dvis_dp, dvis_dpsat
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  PetscReal :: dpw_dp
  PetscReal :: dpsat_dt
  PetscReal :: k_therm_w, k_therm_i
  
  global_auxvar%den_kg(1) = 0.d0

  auxvar%h = 0.d0
  auxvar%u = 0.d0
  kr = 0.d0
 
  global_auxvar%head(1) = xx(1)
  !global_auxvar%temp = xx(2)
    ! RTM: Why is the above commented out?  Is one of these internal 
    ! energy instead of temperature?
 

!***************  Liquid phase properties **************************

  pw = option%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0

  if (global_auxvar%head(1) < MIN_SURFACE_WATER_HEIGHT) then
    global_auxvar%is_dry = PETSC_TRUE
    call EOSWaterDensity(0.0d0,pw,dw_kg,dw_mol,ierr)
    call EOSWaterEnthalpy(0.0d0,pw,hw,ierr)
  else
    global_auxvar%is_dry = PETSC_FALSE
    call EOSWaterDensity(global_auxvar%temp,pw,dw_kg,dw_mol,ierr)
    call EOSWaterEnthalpy(global_auxvar%temp,pw,hw,ierr)
  endif

  ! J/kmol -> whatever units
  hw = hw * option%scale
  
  global_auxvar%den_kg(1) = dw_kg
  !di_kg = 917.d0 ![kg/m^3]
  di_kg = dw_kg

  k_therm_w = 0.57d0 ! [J/s/m/K]
  k_therm_i = 2.18d0 ! [J/s/m/K]
    ! RTM: Same warning for thermal conductivities--these should be computed.
  
  ! RTM: These are being set but we are not using them now.  We should get rid 
  ! of them if we settle on not using an enthalpy formulation for the energy 
  ! equation.
  auxvar%h = hw
  auxvar%u = auxvar%h - pw / dw_mol * option%scale

  ! Compute unfrozen fraction, and then compute the weighted averages of 
  ! density, specific heat capacity, thermal conductivity
  unfrozen_fraction = SurfaceTHAuxVarComputeUnfrozen(global_auxvar%temp)
  auxvar%unfrozen_fraction = unfrozen_fraction
  global_auxvar%den_kg(1) = unfrozen_fraction * dw_kg + (1.d0 - unfrozen_fraction) * di_kg
  auxvar%Cwi = unfrozen_fraction * auxvar%Cw + (1.d0 - unfrozen_fraction) * auxvar%Ci
  auxvar%k_therm = unfrozen_fraction * k_therm_w + (1.d0 - unfrozen_fraction) * k_therm_i
  
  
end subroutine SurfaceTHAuxVarCompute

! ************************************************************************** !

function SurfaceTHAuxVarComputeUnfrozen(temp)
  ! 
  ! This returns the unfrozen fraction, given a temperature.
  ! RTM: This is currently a public function, but is only called within
  ! Surface_TH_Aux_module, so we may want to make it private.
  ! We may also want to move this into its own module if we add some
  ! more sophisticated ways of computing this quantity.
  ! 

  implicit none

  PetscReal :: SurfaceTHAuxVarComputeUnfrozen

  PetscReal :: temp

  ! For now, we simply say that everything is unfrozen if the temperature is 
  ! above the nominal freezing point, and everything is frozen otherwise.
  ! We may want to turn this step function into a steep sigmoidal one for 
  ! numerical reasons.
  ! At some point, we may also want to consider a more detailed model of the 
  ! freezing and thawing process that might partition the frozen and unfrozen 
  ! phases according to a more physical mechanism (maybe using state transition
  ! theory?).
  if (temp > 0.d0) then
    SurfaceTHAuxVarComputeUnfrozen = 1.d0
  else
    SurfaceTHAuxVarComputeUnfrozen = 0.d0
  endif

end function SurfaceTHAuxVarComputeUnfrozen

! ************************************************************************** !

subroutine SurfaceTHAuxDestroy(aux)
  ! 
  ! This routine deallocates an object.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 

  implicit none

  type(Surface_TH_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  if (associated(aux%auxvars)) deallocate(aux%auxvars)
  nullify(aux%auxvars)
  if (associated(aux%auxvars_bc)) deallocate(aux%auxvars_bc)
  nullify(aux%auxvars_bc)
  if (associated(aux%auxvars_ss)) deallocate(aux%auxvars_ss)
  nullify(aux%auxvars_ss)
  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)

  deallocate(aux)
  nullify(aux)  

  end subroutine SurfaceTHAuxDestroy

end module Surface_TH_Aux_module
