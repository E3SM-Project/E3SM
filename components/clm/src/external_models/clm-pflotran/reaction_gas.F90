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
            RTotalGas, &
            RTotalCO2

contains

! ************************************************************************** !

subroutine RGasRead(gas_species_list,gas_type,error_msg,input,option)
  ! 
  ! Reads immobile species
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/02/13/ 08/01/16
  ! 

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
  
  if (option%nphase < 2 .or. option%iflowmode /= G_MODE) return
  
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

subroutine RTotalCO2(rt_auxvar,global_auxvar,reaction,option)
  ! 
  ! Computes the total component concentrations and derivative with
  ! respect to free-ion for CO2 modes; this is legacy cod3
  ! 
  ! Author: Glenn Hammond, but originally by Chuan Lu
  ! Date: 08/01/16
  ! 

  use Option_module
  use EOS_Water_module
  use co2eos_module, only: Henry_duan_sun

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option

  PetscErrorCode :: ierr
  PetscInt :: iphase
  PetscInt :: icomp
  PetscReal :: tempreal
  PetscReal :: dg,dddt,dddp,fg,dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,&
               yco2,pco2,sat_pressure,lngamco2
  PetscReal :: den_kg_per_L
  PetscReal :: den
  PetscReal :: lnQK
  PetscReal :: m_cl, m_na, muco2, xmass
  PetscReal :: pressure, temperature, xphico2
  PetscInt :: ieqgas

! *********** Add SC phase and gas contributions ***********************  
  ! CO2-specific
  iphase = 2

  if (iphase > option%nphase) return 
  rt_auxvar%total(:,iphase) = 0.D0
  rt_auxvar%aqueous%dtotal(:,:,iphase) = 0.D0

  if (associated(global_auxvar%xmass)) xmass = global_auxvar%xmass(iphase)
  den_kg_per_L = global_auxvar%den_kg(iphase)*xmass*1.d-3

  if (global_auxvar%sat(iphase) > 1.D-20) then
    do ieqgas = 1, reaction%gas%npassive_gas ! all gas phase species are secondary

      pressure = global_auxvar%pres(2)
      temperature = global_auxvar%temp
      xphico2 = global_auxvar%fugacoeff(1)
      den = global_auxvar%den(2)
 
      call EOSWaterSaturationPressure(temperature, sat_pressure, ierr)
      pco2 = pressure - sat_pressure
!     call co2_span_wagner(pressure*1.D-6,temperature+273.15D0,dg,dddt,dddp,fg, &
!              dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)
!
!            fg = fg*1D6
!            xphico2 = fg / pco2
!            global_auxvar%fugacoeff(1) = xphico2


      if (abs(reaction%species_idx%co2_gas_id) == ieqgas ) then

        if (reaction%species_idx%na_ion_id /= 0 .and. reaction%species_idx%cl_ion_id /= 0) then
          m_na = rt_auxvar%pri_molal(reaction%species_idx%na_ion_id)
          m_cl = rt_auxvar%pri_molal(reaction%species_idx%cl_ion_id)
          call Henry_duan_sun(temperature,pressure*1D-5,muco2, &
                lngamco2,m_na,m_cl)
        else
          call Henry_duan_sun(temperature,pressure*1D-5,muco2, &
                lngamco2,option%m_nacl,option%m_nacl)
        endif
        !lnQk = - log(muco2) 
        lnQk = - log(muco2)-lngamco2
           
      else   
        lngamco2 = 0.d0
        lnQK = -reaction%gas%acteqlogK(ieqgas)*LOG_TO_LN
      endif 
          
      if (reaction%gas%acteqh2oid(ieqgas) > 0) then
        lnQK = lnQK + reaction%gas%acteqh2ostoich(ieqgas)*rt_auxvar%ln_act_h2o
      endif
   
   ! contribute to %total          
   !     do i = 1, ncomp
   ! removed loop over species, suppose only one primary species is related
      icomp = reaction%gas%acteqspecid(1,ieqgas)
      pressure = pressure * 1.D-5
        
!     rt_auxvar%gas_pp(ieqgas) = &
!         exp(lnQK+lngamco2)*rt_auxvar%pri_molal(icomp) &
!         /(IDEAL_GAS_CONSTANT*1.d-2*(temperature+273.15D0)*xphico2)

!     This form includes factor Z in pV = ZRT for nonideal gas
      rt_auxvar%gas_pp(ieqgas) = &
          exp(lnQK)*rt_auxvar%pri_act_coef(icomp)*rt_auxvar%pri_molal(icomp)* &
          den/pressure/xphico2

      rt_auxvar%total(icomp,iphase) = rt_auxvar%total(icomp,iphase) + &
          reaction%gas%acteqstoich(1,ieqgas)* &
          rt_auxvar%gas_pp(ieqgas)

!       print *,'RTotal: ',icomp,ieqgas,pressure, temperature, xphico2, &
!         global_auxvar%sat(iphase),rt_auxvar%gas_pp(ieqgas), &
!         rt_auxvar%pri_act_coef(icomp)*exp(lnQK)*rt_auxvar%pri_molal(icomp) &
!         /pressure/xphico2*den


   ! contribute to %dtotal
   !      tempreal = exp(lnQK+lngamco2)/pressure/xphico2*den
!     tempreal = rt_auxvar%pri_act_coef(icomp)*exp(lnQK) &
!         /pressure/xphico2*den
      tempreal = rt_auxvar%gas_pp(ieqgas)/rt_auxvar%pri_molal(icomp)
      rt_auxvar%aqueous%dtotal(icomp,icomp,iphase) = &
          rt_auxvar%aqueous%dtotal(icomp,icomp,iphase) + &
          reaction%gas%acteqstoich(1,ieqgas)*tempreal
    enddo
  ! rt_auxvar%total(:,iphase) = rt_auxvar%total(:,iphase)!*den_kg_per_L
  ! units of dtotal = kg water/L water
  ! rt_auxvar%dtotal(:, :,iphase) = rt_auxvar%dtotal(:,:,iphase)!*den_kg_per_L
  endif

end subroutine RTotalCO2
  
end module Reaction_Gas_module
