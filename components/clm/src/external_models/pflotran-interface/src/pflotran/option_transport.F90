module Option_Transport_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: transport_option_type 
  
    PetscInt :: rt_idof
    PetscInt :: reactive_transport_coupling
    PetscInt :: tvd_flux_limiter
    PetscBool :: store_fluxes
    PetscReal :: tran_weight_t0, tran_weight_t1
    
    PetscReal :: inf_rel_update_tol
    PetscReal :: inf_scaled_res_tol
  
    PetscBool :: jumpstart_kinetic_sorption
    PetscBool :: no_checkpoint_kinetic_sorption
    PetscBool :: no_restart_kinetic_sorption
    PetscBool :: no_restart_mineral_vol_frac
    PetscBool :: only_vertical_tran
    PetscBool :: numerical_derivatives

    PetscInt :: nphase
        
  end type transport_option_type
  
  public :: OptionTransportCreate, &
            OptionTransportInitAll, &
            OptionTransportInitRealization, &
            OptionTransportDestroy

contains

! ************************************************************************** !

function OptionTransportCreate()
  ! 
  ! Allocates and initializes a new OptionTransport object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  type(transport_option_type), pointer :: OptionTransportCreate
  
  type(transport_option_type), pointer :: option
  
  allocate(option)

  ! DO NOT initialize members of the option type here.  One must decide 
  ! whether the member needs initialization once for all stochastic 
  ! simulations or initialization for every realization (e.g. within multiple 
  ! stochastic simulations).  This is done in OptionTransportInitAll() and
  ! OptionTransportInitRealization()
  call OptionTransportInitAll(option)
  OptionTransportCreate => option
  
end function OptionTransportCreate

! ************************************************************************** !

subroutine OptionTransportInitAll(option)
  ! 
  ! Initializes all option variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  type(transport_option_type) :: option
  
  ! These variables should only be initialized once at the beginning of a
  ! PFLOTRAN run (regardless of whether stochastic)
  
 
  call OptionTransportInitRealization(option)

end subroutine OptionTransportInitAll

! ************************************************************************** !

subroutine OptionTransportInitRealization(option)
  ! 
  ! Initializes option variables specific to a single
  ! realization
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  type(transport_option_type) :: option
  
  ! These variables should be initialized once at the beginning of every 
  ! PFLOTRAN realization or simulation of a single realization
    
  option%tvd_flux_limiter = 1
  option%rt_idof = UNINITIALIZED_INTEGER
  option%store_fluxes = PETSC_FALSE
  
  option%reactive_transport_coupling = GLOBAL_IMPLICIT
  option%numerical_derivatives = PETSC_FALSE
  
  option%jumpstart_kinetic_sorption = PETSC_FALSE
  option%no_checkpoint_kinetic_sorption = PETSC_FALSE
  option%no_restart_kinetic_sorption = PETSC_FALSE
  option%no_restart_mineral_vol_frac = PETSC_FALSE
  
  option%tran_weight_t0 = 0.d0
  option%tran_weight_t1 = 0.d0

  option%inf_rel_update_tol = UNINITIALIZED_DOUBLE
  option%inf_scaled_res_tol = UNINITIALIZED_DOUBLE 

  option%nphase = 1

  option%only_vertical_tran = PETSC_FALSE
  
end subroutine OptionTransportInitRealization

! ************************************************************************** !

subroutine OptionTransportDestroy(option)
  ! 
  ! Deallocates an option
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  ! 

  implicit none
  
  type(transport_option_type), pointer :: option
  
  if (.not.associated(option)) return
  ! all kinds of stuff needs to be added here.

  ! all the below should be placed somewhere other than option.F90
  
  deallocate(option)
  nullify(option)
  
end subroutine OptionTransportDestroy

end module Option_Transport_module
