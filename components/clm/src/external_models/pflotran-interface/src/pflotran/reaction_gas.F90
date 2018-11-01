module Reaction_Gas_module

  use Reaction_Aux_module
  use Reactive_Transport_Aux_module  
  use Global_Aux_module
  use Reaction_Gas_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
 
  private

#include "petsc/finclude/petscsys.h"

  
  public :: RGasRead, &
            RTotalGas

contains

! ************************************************************************** !

subroutine RGasRead(gas_species_list,gas_type,error_msg,input,option)
  ! 
  ! Reads immobile species
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/02/13/ 08/01/16
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none
  
  type(gas_species_type), pointer :: gas_species_list
  PetscInt :: gas_type
  character(len=MAXSTRINGLENGTH) :: error_msg
  type(input_type), pointer :: input
  type(option_type) :: option
  
  type(gas_species_type), pointer :: new_gas_species, &
                                     prev_gas_species

  ! since both active and passive gases are in the same list, skip to the
  ! end of the list if it exists.
  if (associated(gas_species_list)) then
    prev_gas_species => gas_species_list
    do
      if (.not.associated(prev_gas_species%next)) exit
      prev_gas_species => prev_gas_species%next
    enddo
  else
    nullify(prev_gas_species)
  endif
  ! read in new gases
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    new_gas_species => GasSpeciesCreate()
    call InputReadWord(input,option,new_gas_species%name,PETSC_TRUE)  
    call InputErrorMsg(input,option,'keyword',error_msg)    
    new_gas_species%itype = gas_type
    if (associated(prev_gas_species)) then
      prev_gas_species%next => new_gas_species
      new_gas_species%id = prev_gas_species%id + 1
    else
      gas_species_list => new_gas_species
      new_gas_species%id = 1
    endif
    prev_gas_species => new_gas_species
    nullify(new_gas_species)
  enddo                                          
                                          
end subroutine RGasRead

! ************************************************************************** !

subroutine RTotalGas(rt_auxvar,global_auxvar,reaction,option)
  ! 
  ! Computes the total component concentrations and derivative with
  ! respect to free-ion
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/01/16
  ! 

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt, parameter :: iphase = 2
  PetscInt :: i, j, igas, icomp, jcomp, ncomp
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: lnQK, tempreal
  PetscReal :: RT
  PetscReal :: gas_concentration
  type(gas_type), pointer :: gas
  
  rt_auxvar%total(:,iphase) = 0.d0 !debugging 
  
  gas => reaction%gas
  ! units of ideal gas constant = J/mol-K = kPa-L/mol-K
  ! units of RT = Pa-L/mol
  RT = IDEAL_GAS_CONSTANT*(global_auxvar%temp+273.15d0)*1.d3
  
  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
  
  ! initialize derivatives
  rt_auxvar%aqueous%dtotal(:,:,iphase) = 0.d0
   
  do igas = 1, gas%nactive_gas ! for each secondary species
    ! compute secondary species concentration
    lnQK = -gas%acteqlogK(igas)*LOG_TO_LN

    ! activity of water
    if (gas%acteqh2oid(igas) > 0) then
      lnQK = lnQK + gas%acteqh2ostoich(igas)*rt_auxvar%ln_act_h2o
    endif

    ncomp = gas%acteqspecid(0,igas)
    do i = 1, ncomp
      icomp = gas%acteqspecid(i,igas)
      lnQK = lnQK + gas%acteqstoich(i,igas)*ln_act(icomp)
    enddo
    ! units = bars
    rt_auxvar%gas_pp(igas) = exp(lnQK)
    ! unit = mol/L gas
    gas_concentration = rt_auxvar%gas_pp(igas) * 1.d5 / RT

    ! add contribution to primary totals
    ! units of total = mol/L gas
    do i = 1, ncomp
      icomp = gas%acteqspecid(i,igas)
      rt_auxvar%total(icomp,iphase) = rt_auxvar%total(icomp,iphase) + &
                                      gas%acteqstoich(i,igas)* &
                                      gas_concentration
    enddo
    
    ! add contribution to derivatives of total with respect to free
    ! units of dtotal = kg water / L gas
    do j = 1, ncomp
      jcomp = gas%acteqspecid(j,igas)
      tempreal = gas%acteqstoich(j,igas)*exp(lnQK-ln_conc(jcomp))*1.d5/RT
      do i = 1, ncomp
        icomp = gas%acteqspecid(i,igas)
        rt_auxvar%aqueous%dtotal(icomp,jcomp,iphase) = &
          rt_auxvar%aqueous%dtotal(icomp,jcomp,iphase) + &
          gas%acteqstoich(i,igas)*tempreal
      enddo
    enddo
  enddo

end subroutine RTotalGas

! ************************************************************************** !

  
end module Reaction_Gas_module
