module Init_Subsurface_Tran_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private


  public :: InitSubsurfTranSetupRealization
  
contains

! ************************************************************************** !

subroutine InitSubsurfTranSetupRealization(realization)
  ! 
  ! Initializes material property data structres and assign them to the domain.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  ! 
  use Realization_Subsurface_class
  use Option_module
  
  use Reactive_Transport_module
  use Global_module
  use Condition_Control_module
  use Variables_module
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  
  option => realization%option
  
  call RTSetup(realization)

  ! initialize densities and saturations
  if (option%nflowdof == 0) then
    call GlobalSetAuxVarScalar(realization,option%reference_pressure, &
                               LIQUID_PRESSURE)
    call GlobalSetAuxVarScalar(realization,option%reference_temperature, &
                               TEMPERATURE)
    call GlobalSetAuxVarScalar(realization,option%reference_saturation, &
                               LIQUID_SATURATION)
    call GlobalSetAuxVarScalar(realization, &
                               option%reference_density(option%liquid_phase), &
                               LIQUID_DENSITY)
    if (option%transport%nphase > 1) then
      call GlobalSetAuxVarScalar(realization,1.d0-option%reference_saturation, &
                                 GAS_SATURATION)
      call GlobalSetAuxVarScalar(realization, &
                                 option%reference_density(option%gas_phase), &
                                 GAS_DENSITY)
    endif
  else
    call GlobalUpdateAuxVars(realization,TIME_T,0.d0)
    call GlobalWeightAuxVars(realization,0.d0)
  endif

  ! initial concentrations must be assigned after densities are set !!!
  call CondControlAssignTranInitCond(realization)
  
end subroutine InitSubsurfTranSetupRealization

end module Init_Subsurface_Tran_module
