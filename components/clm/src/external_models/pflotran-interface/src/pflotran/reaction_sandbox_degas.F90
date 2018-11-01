module Reaction_Sandbox_degas_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Sandbox_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module
  
  implicit none
  
  private
  
  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_degas_type
    PetscInt  :: ispec_co2a, ispec_n2oa, ispec_n2a
    PetscInt  :: ispec_co2g, ispec_n2og, ispec_n2g
    PetscInt  :: ispec_proton
    PetscInt  :: ispec_himm
    PetscReal :: k_kinetic_co2, k_kinetic_n2o, k_kinetic_n2
    PetscReal :: k_kinetic_h
    PetscBool :: b_fixph
    PetscReal :: fixph

  contains
    procedure, public :: ReadInput => degasRead
    procedure, public :: Setup => degasSetup
    procedure, public :: Evaluate => degasReact
    procedure, public :: Destroy => degasDestroy
  end type reaction_sandbox_degas_type

  public :: degasCreate

contains

! ************************************************************************** !
!
! degasCreate: Allocates degas reaction object.
!
! ************************************************************************** !
function degasCreate()

  implicit none
  
  class(reaction_sandbox_degas_type), pointer :: degasCreate

! 4. Add code to allocate object and initialized all variables to zero and
!    nullify all pointers. E.g.
  allocate(degasCreate)
  degasCreate%ispec_co2a = 0
  degasCreate%ispec_co2g = 0
  degasCreate%ispec_n2oa = 0
  degasCreate%ispec_n2og = 0
  degasCreate%ispec_n2a = 0
  degasCreate%ispec_n2g = 0
  degasCreate%ispec_proton = 0
  degasCreate%k_kinetic_co2 = 1.d-5
  degasCreate%k_kinetic_n2o = 1.d-5
  degasCreate%k_kinetic_n2 = 1.d-5
  degasCreate%k_kinetic_h = 1.d-5
  degasCreate%fixph = 6.5d0
  degasCreate%b_fixph = PETSC_FALSE
  nullify(degasCreate%next)  
      
end function degasCreate

! ************************************************************************** !
!
! degasRead: Reads input deck for degas reaction parameters (if any)
!
! ************************************************************************** !
subroutine degasRead(this,input,option)

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_degas_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,DEGAS')
    call StringToUpper(word)   

    select case(trim(word))
      case('KINETIC_CONSTANT_CO2')
         call InputReadDouble(input,option,this%k_kinetic_co2)
         call InputErrorMsg(input,option,'CO2 degas kinetic rate constant', &
                     'CHEMISTRY,REACTION_SANDBOX,DEGAS,REACTION')
      case('KINETIC_CONSTANT_N2O')
         call InputReadDouble(input,option,this%k_kinetic_n2o)
         call InputErrorMsg(input,option,'N2O degas kinetic rate constant', &
                     'CHEMISTRY,REACTION_SANDBOX,DEGAS,REACTION')
      case('KINETIC_CONSTANT_N2')
         call InputReadDouble(input,option,this%k_kinetic_n2)
         call InputErrorMsg(input,option,'N2 degas kinetic rate constant', &
                     'CHEMISTRY,REACTION_SANDBOX,DEGAS,REACTION')
      case('KINETIC_CONSTANT_H+')
         call InputReadDouble(input,option,this%k_kinetic_h)
         call InputErrorMsg(input,option,'H+ fix pH  kinetic rate constant', &
                     'CHEMISTRY,REACTION_SANDBOX,DEGAS,REACTION')
      case('FIXPH')
         call InputReadDouble(input,option,this%fixph)
         call InputErrorMsg(input,option,'fix ph', &
                     'CHEMISTRY,REACTION_SANDBOX,DEGAS,REACTION')
         this%b_fixph = PETSC_TRUE
      case default
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DEGAS,' // &
            'REACTION keyword: ' // trim(word) // ' not recognized.'
          call printErrMsg(option)
    end select
  enddo
  
end subroutine degasRead

! ************************************************************************** !
!
! degasSetup: Sets up the degas reaction either with parameters either
!                read from the input deck or hardwired.
! author: Guoping Tang
! date: 04/03/2014
!
! ************************************************************************** !
subroutine degasSetup(this,reaction,option)

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Option_module
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName

  implicit none
  
  class(reaction_sandbox_degas_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word

  word = 'CO2(aq)'
  this%ispec_co2a = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'N2O(aq)'
  this%ispec_n2oa = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'N2(aq)'
  this%ispec_n2a = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  ! (TODO) currently, gas CO2 related process not ready, so it's assumed as immobile species.
  word = 'CO2imm'
  this%ispec_co2g = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
  word = 'N2Oimm'
  this%ispec_n2og = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)

  word = 'N2imm'
  this%ispec_n2g = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)

  if(this%b_fixph) then
     word = 'H+'
     this%ispec_proton = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

     word = 'Himm'
     this%ispec_himm = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
   
     if(this%ispec_proton < 0) then
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DEGAS,' // &
            'H+ is not defined even though pH needs to be fixed!'
        call printErrMsg(option)
     endif
 
     if(this%ispec_himm < 0) then
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DEGAS,' // &
            'Himm is not defined even though pH needs to be fixed!'
        call printErrMsg(option)
     endif
  endif

end subroutine degasSetup

! ************************************************************************** !
!
! degasReact: Evaluates reaction storing residual and/or Jacobian
!
! ************************************************************************** !
subroutine degasReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class, only : material_auxvar_type


  
  implicit none


  class(reaction_sandbox_degas_type) :: this  
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscErrorCode :: ierr

  PetscInt, parameter :: iphase = 1
  PetscReal, parameter :: H2O_kg_mol = 18.01534d-3 ! kg mol-1 for pure water
  PetscReal, parameter :: rgas = 8.3144621d0       ! m3 Pa K-1 mol-1

  PetscInt :: ires_co2a, ires_co2g
  PetscInt :: ires_n2oa, ires_n2og
  PetscInt :: ires_n2a, ires_n2g
  PetscInt :: ires_proton, ires_himm

  PetscReal :: c_co2_aq, c_n2o_aq, c_n2_aq         ! actual gas solution conc.: mole/L
  PetscReal :: c_co2_eq, c_n2o_eq, c_n2_eq         ! gas solubility: mole/L
  PetscReal :: tc
  PetscReal :: volume, porosity, lsat, isat, total_sal
  PetscReal :: air_press, air_vol, air_molar
  PetscReal :: co2_p, co2_molar, xmole_co2, xmass_co2   ! in air, or in solution (x)
  PetscReal :: co2_rho, co2_fg, co2_xphi, co2_henry, co2_poyn
  PetscReal :: n2o_p, n2o_molar, xmole_n2o, xmass_n2o   ! in air, or in solution (x)
  PetscReal :: n2o_fg, n2o_xphi, n2o_henry
  PetscReal :: n2_p, n2_molar, xmole_n2, xmass_n2       ! in air, or in solution (x)
  PetscReal :: n2_fg, n2_xphi, n2_henry
  PetscReal :: c_h, c_h_fix
  PetscReal :: temp_real, rate, drate
  PetscReal :: convert_molal_to_molar
  PetscReal :: xmass
 
  PetscInt :: ires
!-------------------------------------------------------------------------------------

  xmass = 1.d0

  if (associated(global_auxvar%xmass)) xmass = global_auxvar%xmass(iphase)
  
  if (reaction%initialize_with_molality) then
    convert_molal_to_molar = global_auxvar%den_kg(iphase)*xmass/1000.d0
  else
    convert_molal_to_molar = 1.d0
  endif

  !default values for calculating gas solubility
  tc = option%reference_temperature
  air_press = option%reference_pressure
  lsat = 0.50d0  ! 50% saturation assumed as default
  isat = 0.d0

  if (option%iflowmode == TH_MODE) then

      air_press = max(air_press, global_auxvar%pres(1))      ! total (air)gas pressure: water pressure if over atm. press., otherwise atm. press.
      lsat = global_auxvar%sat(1)
      if (option%iflowmode == TH_MODE) then
         tc = global_auxvar%temp
         !isat = th_auxvar%sat_ice       ! (TODO) not yet figure out how to point to 'th_auxvar'
      endif
#ifdef CLM_PFLOTRAN
  elseif (option%ntrandof.gt.0 ) then
      air_press = max(air_press, global_auxvar%pres(1))      ! total (air)gas pressure: water pressure if over atm. press., otherwise atm. press.
      lsat = global_auxvar%sat(1)
      tc = global_auxvar%temp
  else
      option%io_buffer='reaction_sandbox_degas ' // &
                 'not supported for the modes applied for.'
      call printErrMsg(option)
#endif
  endif

  porosity = material_auxvar%porosity
  volume = material_auxvar%volume
  !air_vol = max(0.01d0, porosity * (1.d0-lsat-isat))               ! min. 0.01 to avoid math. issue (temporarily off - TODO)
  air_vol = 1.d0                                                    ! atm. air volume fraction (no adjusting of soil air volume)
  air_molar = air_press/rgas/(tc+273.15d0)                          ! molAir/m3

!
!------------------------------------------------------------------------------------
! co2(aq) <==> co2(g)
!
  if (this%ispec_co2a > 0 .and. this%ispec_co2g > 0) then

    ires_co2a = this%ispec_co2a
    ires_co2g = this%ispec_co2g + reaction%offset_immobile

    c_co2_aq = rt_auxvar%total(this%ispec_co2a, iphase)

    co2_p = 350.0d-6 * option%reference_pressure
#ifdef CLM_PFLOTRAN
    ! resetting 'co2g' from CLM after adjusting
    if (this%ispec_co2g > 0) then
       co2_molar = rt_auxvar%immobile(this%ispec_co2g)/air_vol          ! molCO2/m3 bulk soil --> mol/m3 air space
       co2_p = co2_molar/air_molar*air_press                            ! mole fraction --> Pa
    endif
#endif

    ! the 'Duan's equations appear generates 'infinity' at low CO2 pressure (probably due to for deep earth geochemistry?)
    !temp_real = max(min(tc,65.d0), 1.d-20)                                   ! 'duanco2' only functions from 0 - 65oC
    !call duanco2(temp_real,co2_p, co2_rho, co2_fg, co2_xphi)     ! only need 'co2_xphi' (fugacity coefficient) for the following call
    !call Henry_CO2_noderiv(xmole_co2,xmass_co2,temp_real,co2_p,co2_xphi,co2_henry,co2_poyn)   ! 'xmolco2': mol fraction; 'xmco2': mass fraction (CO2:CO2+H2O)

    ! Weiss (1974)
    temp_real = max(min(tc,40.d0), -1.d0)                         ! Weiss (1974) functions from -1 - 40oC
    total_sal = 1.d-20                                            ! no accounting for salinity
    call weiss_co2(temp_real, air_press, total_sal, co2_p, &
         xmole_co2, xmass_co2, co2_henry, co2_fg, co2_xphi)
    c_co2_eq =  xmole_co2/H2O_kg_mol                              ! moleCO2/mole solution -> moleCO2/kg (L) water solution

    temp_real = volume * 1000.0d0 * porosity * lsat                      ! kgH2O
    !rate = this%k_kinetic_co2 * (c_co2_aq/c_co2_eq-1.0d0) * temp_real   ! moles/second: moles/kgW/s * (-) * kgW
    rate = this%k_kinetic_co2 * (c_co2_aq - c_co2_eq) * temp_real        ! moles/second: 1/s * moles/kgW * kgW

    ! degas occurs if over-saturated, or gas dissolves if high gas conc.
    if(abs(rate) > 1.0d-20) then

      Residual(ires_co2a) = Residual(ires_co2a) + rate
      Residual(ires_co2g) = Residual(ires_co2g) - rate

      if (compute_derivative) then
        !drate = this%k_kinetic_co2 /c_co2_eq * temp_real
        drate = this%k_kinetic_co2 * temp_real                    ! ????? -- this seems incorrect (TODO)

        Jacobian(ires_co2a,ires_co2a) = Jacobian(ires_co2a,ires_co2a) + drate * &
                rt_auxvar%aqueous%dtotal(this%ispec_co2a,this%ispec_co2a,iphase)

        Jacobian(ires_co2g,ires_co2a) = Jacobian(ires_co2g,ires_co2a) - drate
      endif

    endif

  endif ! end of 'if (this%ispec_co2a > 0 .and. this%ispec_co2g > 0)' block

!
!------------------------------------------------------------------------------------
! n2o(aq) <==> n2o(g)
!
  if (this%ispec_n2oa > 0 .and. this%ispec_n2og > 0) then

    ires_n2oa = this%ispec_n2oa
    ires_n2og = this%ispec_n2og + reaction%offset_immobile

    c_n2o_aq = rt_auxvar%total(this%ispec_n2oa, iphase)

    n2o_p = 310.0d-9 * option%reference_pressure                        ! default (310 ppbv N2O in atm. in about 1990s)
#ifdef CLM_PFLOTRAN
    ! resetting 'n2og' from CLM after adjusting via 'N2Oimm'
    if (this%ispec_n2og > 0) then
       n2o_molar = rt_auxvar%immobile(this%ispec_n2og)/air_vol          ! molN2O/m3 bulk soil --> mol/m3 air space
       n2o_p = n2o_molar/air_molar*air_press                            ! mole fraction --> Pa
    endif
#endif

    temp_real = max(min(tc,40.d0), 1.d-20)
    total_sal = 1.0d-20
    call weiss_price_n2o (temp_real, air_press, total_sal, n2o_p, &
         xmole_n2o, xmass_n2o, n2o_henry, n2o_fg, n2o_xphi)
    c_n2o_eq =  xmole_n2o/H2O_kg_mol                                    ! moleN2O/mole solution -> moleN2O/kg (L) water solution

    temp_real = volume * 1000.0d0 * porosity * lsat                     ! kgH2O (L)
    !rate = this%k_kinetic_n2 * (c_n2o_aq/c_n2o_eq-1.0d0) * temp_real   ! moles/second: moles/kgW/s * (-) * kgW
    rate = this%k_kinetic_n2o * (c_n2o_aq - c_n2o_eq) * temp_real       ! moles/second: 1/s * moles/kgW * kgW

    ! degas occurs if over-saturated, or gas dissolves if high gas conc.
    if(abs(rate) > 1.0d-20) then

      Residual(ires_n2oa) = Residual(ires_n2oa) + rate
      Residual(ires_n2og) = Residual(ires_n2og) - rate

      if (compute_derivative) then
        !drate = this%k_kinetic_n2o /c_n2o_eq * temp_real
        drate = this%k_kinetic_n2o * temp_real

        Jacobian(ires_n2oa,ires_n2oa) = Jacobian(ires_n2oa,ires_n2oa) + drate * &
            rt_auxvar%aqueous%dtotal(this%ispec_n2oa,this%ispec_n2oa,iphase)

        Jacobian(ires_n2og,ires_n2oa) = Jacobian(ires_n2og,ires_n2oa) - drate
      endif

    endif

  endif ! end of 'if (this%ispec_n2oa > 0 .and. this%ispec_n2og > 0)' block

!
!------------------------------------------------------------------------------------
! n2(aq) <==> n2(g)
!
  if (this%ispec_n2a > 0 .and. this%ispec_n2g > 0) then
    ires_n2a = this%ispec_n2a
    ires_n2g = this%ispec_n2g + reaction%offset_immobile

    c_n2_aq = rt_auxvar%total(this%ispec_n2a, iphase)                   ! M (mole-N/L)

    n2_p = 0.78084d0 * option%reference_pressure                        ! default
#ifdef CLM_PFLOTRAN
    ! resetting 'n2g' from CLM after adjusting via 'N2imm'
    if (this%ispec_n2g > 0) then
       n2_molar = rt_auxvar%immobile(this%ispec_n2g)/air_vol          ! molN2-N/m3 bulk soil --> mol/m3 air space
       n2_p = n2_molar/air_molar*air_press                            ! mole fraction --> Pa
    endif
#endif

    temp_real = max(min(tc,40.d0), -2.d0)
    total_sal = 1.0d-20
    call weiss_n2(temp_real, air_press, total_sal, n2_p, &
         xmole_n2, xmass_n2, n2_henry)
    c_n2_eq =  xmole_n2/H2O_kg_mol                            ! moleN2/mole solution -> moleN2/kg (L) water solution

    temp_real = volume * porosity * lsat * 1.d3                       ! kgH2O (L) <== m3
    !rate = this%k_kinetic_n2 * (c_n2_aq/c_n2_eq-1.0d0) * temp_real   ! moles/second: moles/kgW/s * (-) * kgW
    rate = this%k_kinetic_n2 * (c_n2_aq - c_n2_eq) * temp_real        ! moles/second: 1/s * moles/kgW * kgW

    ! degas occurs if over-saturated, or gas dissolves if high gas conc.
    if(abs(rate) > 1.0d-20) then

      Residual(ires_n2a) = Residual(ires_n2a) + rate
      Residual(ires_n2g) = Residual(ires_n2g) - rate

      if (compute_derivative) then
        !drate = this%k_kinetic_n2 /c_n2_eq * temp_real
        drate = this%k_kinetic_n2 * temp_real

        Jacobian(ires_n2a,ires_n2a) = Jacobian(ires_n2a,ires_n2a) + drate * &
            rt_auxvar%aqueous%dtotal(this%ispec_n2a,this%ispec_n2a,iphase)

        Jacobian(ires_n2g,ires_n2a) = Jacobian(ires_n2g,ires_n2a) - drate

      endif
    endif

  endif ! end of 'if(this%ispec_n2a > 0 .and. this%ispec_n2g > 0)' block

!
!------------------------------------------------------------------------------------------------
!
  ! the following is for setting pH at a fixed value from input
  ! (TODO) if PFLOTRAN H+/OH- related processes are in place, this should be removed
  if(this%b_fixph) then
    ires_proton = this%ispec_proton
    ires_himm = this%ispec_himm + reaction%offset_immobile

    c_h = rt_auxvar%pri_molal(this%ispec_proton) * convert_molal_to_molar

    c_h_fix = 10.0d0 ** (-1.0d0 * this%fixph) / rt_auxvar%pri_act_coef(this%ispec_proton)

    temp_real = volume * 1000.0d0 * porosity * lsat        ! kgH2O
    if(PETSC_FALSE) then
      rate = this%k_kinetic_h * (c_h/c_h_fix - 1.0d0) * temp_real
    else
      rate = this%k_kinetic_h * (c_h - c_h_fix) * temp_real
    endif  
     
    if(abs(rate) > 1.0d-20) then
       Residual(ires_proton) = Residual(ires_proton) + rate
       Residual(ires_himm) = Residual(ires_himm) - rate

       if (compute_derivative) then
          if (PETSC_FALSE) then
             drate = this%k_kinetic_h * convert_molal_to_molar /c_h_fix * temp_real
          else
             drate = this%k_kinetic_h * convert_molal_to_molar * temp_real
          endif
          Jacobian(ires_proton,ires_proton) = Jacobian(ires_proton,ires_proton) + drate
          Jacobian(ires_himm,ires_proton) = Jacobian(ires_proton,ires_himm) - drate
       endif

    endif

  endif 

#ifdef DEBUG
  do ires=1, reaction%ncomp
    temp_real = Residual(ires)

    if (abs(temp_real) > huge(temp_real)) then
      write(option%fid_out, *) 'infinity of Residual matrix checking at ires=', ires
      write(option%fid_out, *) 'Reaction Sandbox: DEGAS'

      option%io_buffer = 'checking infinity of Residuals matrix @ degasReact '
      call printErrMsg(option)
    endif

    if (temp_real /= temp_real) then
      write(option%fid_out, *) 'NaN of Residual matrix checking at ires=', ires
      write(option%fid_out, *) 'Reaction Sandbox: DEGAS'
      option%io_buffer = ' checking NaN of Residuals matrix  @ degasReact '
      call printErrMsg(option)
    endif

  enddo
#endif

end subroutine degasReact

! ************************************************************************** !
!
! degasDestroy: Destroys allocatable or pointer objects created in this 
!                  module
!
! ************************************************************************** !
subroutine degasDestroy(this)

  implicit none
  
  class(reaction_sandbox_degas_type) :: this  

end subroutine degasDestroy

! ************************************************************************** !

  subroutine weiss_co2 (tt, tp, ts, pco2, xmole, xmass, kh, fg, phi)
  !
  ! Weiss R. F., 1974. Carbon dioxide in water and seawater: the solubility of
  ! a non-ideal gas. Marine chemistry, 2(1974)203-215
  !
  ! Input: tt   [C]        temperature
  !        tp   [Pa]       total air pressure
  !        ts   [%o]       total salinity (parts per thousands)
  !        pco2 [Pa]       CO2 partial pressure
  ! Output: fg  [Pa]       CO2 fugacity
  !         phi [-]        CO2 fugacity coefficient
  !         kh  [Pa]       CO2 Henry's Law Constant
  !         xmole  [-]     mole fraction of CO2 solubility
  !         xmass  [-]     mass fraction of CO2 solubility

      implicit none

#include "petsc/finclude/petscsys.h"

      PetscReal, intent(in) :: tt, tp, ts, pco2
      PetscReal, intent(out):: xmole, xmass, kh, fg, phi

      PetscReal, parameter :: atm  = 1.01325d5        ! Pa of 1 atm
      PetscReal, parameter :: rgas = 0.08205601       ! L atm mol-1 K-1
      PetscReal, parameter :: xmwco2  = 44.0095d-3    ! kg mol-1
      PetscReal, parameter :: xmwh2o  = 18.01534d-3   ! kg mol-1

      PetscReal :: tk, tk2, tk3, tk_100k
      PetscReal :: p_rt
      PetscReal :: x1,x2
      PetscReal :: k0, epsilon, bt
      PetscReal :: vbar, cco2

      PetscReal :: a1,a2,a3,b1,b2,b3

      ! empirical parameters for temperature effect (for K0 in moles/L/atm)
      data a1    /-58.0931d0/
      data a2    / 90.5069d0/
      data a3    / 22.2940d0/
      ! empirical parameters for salt-water effect
      data b1    / 0.027766d0/
      data b2    /-0.025888d0/
      data b3    / 0.0050578d0/

      ! variable values transformation
      tk = tt + 273.15d0
      tk2= tk*tk
      tk3= tk2*tk
      tk_100k = tk/100.d0
      p_rt = (tp/atm)/rgas/tk                       ! mol L-1
      p_rt = p_rt/1000.d0                           ! mol cm-3  (this conversion needed? - the original paper did say)
      x1 = pco2/tp                                  ! mole fractions of binary mixture of co2-air
      x2 = 1.d0 - x1

      ! fugacity of co2 gas
      epsilon = 57.7d0 - 0.118d0*tk                                 ! eq.(11): cm3/mol
      bt = -1636.75d0+12.0408d0*tk-3.27957d-2*tk2+3.16528d-5*tk3    ! eq.(6): cm3/mol
      fg = pco2*dexp((bt+2.d0*x2*x2*epsilon)*p_rt)                  ! eq.(9): Pa (pco2=xl*tp)

      ! fugacity coefficient
      phi = fg / pco2

      ! K0 in Weiss(1974)'s paper
      k0 = a1+a2/tk_100k+a3*dlog(tk_100k) + &
           ts*(b1+b2*tk_100k+b3*tk_100k*tk_100k)    ! eq. (12): ln(mol/L/atm)
      k0 = dexp(k0)                                 ! mol/L/atm
      k0 = k0/atm                                   ! mol/L/pa

      ! solubility
      vbar = dexp((1.0d0-tp/atm)*30.0d-3/rgas/tk)    ! 30.0 is a general value in cm3/mol for converting dissolved co2 molar to vol in solution (eq.(5) the 3rd term)
      cco2 = k0*fg*vbar                              ! eq. (5): mol/L
      xmole = cco2/(1.d0/xmwh2o)                     ! mole fraction: assuming solution volume is all water (1L=1kgH2O)
      xmass = xmole*xmwco2/(xmole*xmwco2+(1.d0-xmole)*xmwh2o)    ! mass fraction

      ! Henry's Law constant: kh = pco2/[xmole]
      kh = 1.d0/k0                      ! Pa L mol-1: pco2 in pa, solubility in mol/L, with formula: [C]=Pc/kh
      kh = kh*(1.0d0/xmwh2o)            ! Pa (pa mol mol-1): pn2o in pa, solubility in mole fraction (1Lsolution = 1kgH2O)

      return

    end subroutine weiss_co2

! ************************************************************************** !

  subroutine weiss_price_n2o (tt, tp, ts, pn2o, xmole, xmass, kh, fg, phi)
  !
  ! Weiss and Price, 1980. Nitrous oxide solubility in water and seawater. Marine
  ! chemistry, 8(1980)347-359
  !
  ! Input: tt   [C]        temperature
  !        tp   [Pa]       total air pressure
  !        ts   [%o]       total salinity (parts per thousands)
  !        pn2o [Pa]       N2O partial pressure
  ! Output: fg  [Pa]       N2O fugacity
  !         phi [-]        N2O fugacity coefficient
  !         kh  [Pa]       N2O Henry's Law Constant
  !         xmole  [-]     mole fraction of N2O solubility
  !         xmass  [-]     mass fraction of N2O solubility

      implicit none

#include "petsc/finclude/petscsys.h"

      PetscReal, intent(in) :: tt, tp, ts, pn2o
      PetscReal, intent(out):: xmole, xmass, kh, fg, phi

      PetscReal, parameter :: atm  = 1.01325d5        ! Pa of 1 atm
      PetscReal, parameter :: rgas = 0.08205601       ! L atm mol-1 K-1
      PetscReal, parameter :: xmwn2o  = 44.01287d-3   ! kg mol-1
      PetscReal, parameter :: xmwh2o  = 18.01534d-3   ! kg mol-1

      PetscReal :: tk, tk2, tk_100k
      PetscReal :: p_rt
      PetscReal :: x1,x2
      PetscReal :: k0, epsilon, bt
      PetscReal :: vbar, cn2o

      PetscReal :: a1,a2,a3,b1,b2,b3

      ! empirical parameters for temperature effect
      data a1    /-62.7076d0/
      data a2    / 97.3066d0/
      data a3    / 24.1406d0/
      ! empirical parameters for salt-water effect
      data b1    /-0.058420d0/
      data b2    / 0.033193d0/
      data b3    /-0.0051313d0/

      ! variable values transformation
      tk = tt + 273.15d0
      tk2= tk*tk
      tk_100k = tk/100.d0
      p_rt = (tp/atm)/rgas/tk                       ! mol L-1
      p_rt = p_rt/1000.d0                           ! mol cm-3  (this conversion needed? - the original paper did say)
      x1 = pn2o/tp                                  ! mole fractions of binary mixture of n2o-air
      x2 = 1.d0 - x1

      ! fugacity of n2o gas
      epsilon = 65.0d0 - 0.1338d0*tk                ! eq.(6): cm3/mol
      bt = -905.95d0+4.1685d0*tk-0.0052734d0*tk2    ! eq.(4): cm3/mol
      fg = pn2o*dexp((bt+2.d0*x2*x2*epsilon)*p_rt)  ! eq.(5): Pa (so, pn2o is moist-air based)

      ! fugacity coefficient
      phi = fg / pn2o

      ! K0 in Weiss and Price (1980)'s paper
      k0 = a1+a2/tk_100k+a3*dlog(tk_100k) + &
           ts*(b1+b2*tk_100k+b3*tk_100k*tk_100k)    ! eq. (12): ln(mol/L/atm)
      k0 = dexp(k0)                                 ! mol/L/atm
      k0 = k0/atm                                   ! mol/L/pa

      ! solubility
      vbar = dexp((1.0d0-tp/atm)*32.3d-3/rgas/tk)    ! 32.3 is in cm3/mol for converting dissolved n2o molar to vol in solution (unit inconsistency??)
      cn2o = k0*fg*vbar                              ! eq. (1): mol/L
      xmole = cn2o/(1.d0/xmwh2o)                     ! mole fraction: assuming solution volume is all water (1L=1kgH2O)
      xmass = xmole*xmwn2o/(xmole*xmwn2o+(1.d0-xmole)*xmwh2o)    ! mass fraction

      ! Henry's Law constant: kh = pn2o/[xmole]
      kh = 1.d0/k0                      ! Pa L mol-1: pn2o in pa, solubility in mol/L
      kh = kh*(1.0d0/xmwh2o)            ! Pa (pa mol mol-1): pn2o in pa, solubility in mole fraction (1Lsolution = 1kgH2O)

      return

    end subroutine weiss_price_n2o

! ************************************************************************** !

  subroutine weiss_n2 (tt, tp, ts, pn2, xmole, xmass, kh)
  !
  ! Weiss, R. F., 1970. The solubility of nitrogen, oxygen and argon in water and seawater.
  ! Deep-Sea Research, 17(1970)721-735.
  !
  ! Input: tt   [C]        temperature (-2 ~ 40 oC)
  !        tp   [Pa]       total air pressure
  !        ts   [%o]       total salinity (parts per thousands, per mil) ( 0 - 40)
  !        pn2  [Pa]       N2 partial pressure
  ! Output:
  !         kh  [Pa]       N2 Henry's Law Constant
  !         xmole  [-]     mole fraction of N2 solubility
  !         xmass  [-]     mass fraction of N2 solubility

      implicit none

#include "petsc/finclude/petscsys.h"

      PetscReal, intent(in) :: tt, tp, ts, pn2
      PetscReal, intent(out):: xmole, xmass, kh

      PetscReal, parameter :: atmn2= 0.78084d0         ! mole fraction of N2 in atm
      PetscReal, parameter :: atm  = 1.01325d5         ! Pa of 1 atm
      PetscReal, parameter :: rgas = 0.08205601       ! L atm mol-1 K-1
      PetscReal, parameter :: xmwn2   = 28.01344d-3   ! kg mol-1
      PetscReal, parameter :: xmwh2o  = 18.01534d-3   ! kg mol-1

      PetscReal :: tk, tk_100k
      PetscReal :: k0, phi
      PetscReal :: cn2

      PetscReal :: a1,a2,a3,a4,b1,b2,b3

      ! empirical parameters for temperature effect, under moist air at 1 atm.
      data a1    /-172.4965d0/
      data a2    / 248.4262d0/
      data a3    / 143.3483d0/
      data a4    /-21.7120d0/
      ! empirical parameters for salt-water effect, under moist air at 1 atm
      data b1    /-0.049781d0/
      data b2    /-0.025018d0/
      data b3    /-0.0034861d0/

      ! variable values transformation
      tk = tt + 273.15d0
      tk_100k = tk/100.d0

      ! solubility at STP (0.101325 MPa and 298.15 K) in Weiss (1970)'s paper
      k0 = a1+a2/tk_100k+a3*dlog(tk_100k) + a4*tk_100k + &
           ts*(b1+b2*tk_100k+b3*tk_100k*tk_100k)    ! eq. (4): ln(ml-n2(STP)/L-H2O)
      k0 = dexp(k0)                                 ! ml-N2(STP)/L-H2O
      k0 = (k0*1.d-3)/rgas/298.15d0                 ! mol-N2(STP)/L-H2O

      ! Henry's Law constant at STP (and atm. pN2) : kh = pn2/[xmole]
      kh = (atmn2*atm)/k0               ! Pa L mol-1: pn2 in pa, solubility in mol/L
      kh = kh*(1.0d0/xmwh2o)            ! Pa (pa mol mol-1): pn2 in pa, solubility in mole fraction (1Lsolution = 1kgH2O)

      ! fugacity coefficient(not yet - TODO)
      phi = 1.d0

      ! tp adjusting (not yet - TODO)


      ! solubility at pN2 as input, assuming ideal gas and constant Kh
      cn2 = (pn2*phi)/kh                            ! mol/L
      xmole = cn2/(1.d0/xmwh2o)                     ! mole fraction: assuming solution volume is all water (1L=1kgH2O)
      xmass = xmole*xmwn2/(xmole*xmwn2+(1.d0-xmole)*xmwh2o)    ! mass fraction

      return

    end subroutine weiss_n2

end module Reaction_Sandbox_degas_class
