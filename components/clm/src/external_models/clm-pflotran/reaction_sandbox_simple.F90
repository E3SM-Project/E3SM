module Reaction_Sandbox_Simple_class

  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_simple_type
    ! Aqueous species
    PetscInt :: species_Aaq_id
    PetscInt :: species_Baq_id
    PetscInt :: species_Caq_id
    PetscInt :: species_Daq_id
    PetscInt :: species_Eaq_id
    PetscInt :: species_Faq_id
    ! Immobile species (e.g. biomass)
    PetscInt :: species_Xim_id
    PetscInt :: species_Yim_id
  contains
    procedure, public :: Setup => SimpleSetup
    procedure, public :: Evaluate => SimpleReact
  end type reaction_sandbox_simple_type

  public :: SimpleCreate

contains

! ************************************************************************** !

function SimpleCreate()
  ! 
  ! Allocates simple reaction object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/03/15

  implicit none
  
  class(reaction_sandbox_simple_type), pointer :: SimpleCreate

  allocate(SimpleCreate)
  nullify(SimpleCreate%next)  
      
end function SimpleCreate

! ************************************************************************** !

subroutine SimpleSetup(this,reaction,option)
  ! 
  ! Sets up the simple reaction with hardwired parameters
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/03/15

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  use Option_module

  implicit none
  
  class(reaction_sandbox_simple_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word

  ! Aqueous species
  word = 'Aaq'
  this%species_Aaq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Baq'
  this%species_Baq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Caq'
  this%species_Caq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Daq'
  this%species_Daq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Eaq'
  this%species_Eaq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Faq'
  this%species_Faq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)

  ! Immobile species
  word = 'Xim'
  this%species_Xim_id = &
    GetImmobileSpeciesIDFromName(word,reaction%immobile,option)
  word = 'Yim'
  this%species_Yim_id = &
    GetImmobileSpeciesIDFromName(word,reaction%immobile,option)

      
end subroutine SimpleSetup

! ************************************************************************** !

subroutine SimpleReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/03/15
  ! 

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  
  implicit none
  
  class(reaction_sandbox_simple_type) :: this  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: volume                 ! m^3 bulk
  PetscReal :: porosity               ! m^3 pore space / m^3 bulk
  PetscReal :: liquid_saturation      ! m^3 water / m^3 pore space
  PetscReal :: L_water                ! L water
  
  PetscReal :: Aaq, Baq, Caq, Daq, Eaq, Faq  ! mol/L
  PetscReal :: Xim, Yim  ! mol/m^3
  PetscReal :: Rate,Rate1,Rate2
  PetscReal :: RateA, RateB, RateC, RateD, RateE, RateF, RateX, RateY  ! mol/sec
  PetscReal :: stoichA, stoichB, stoichC, stoichD, stoichE, stoichF
  PetscReal :: stoichX, stoichY
  PetscReal :: k, kr, k1, k2  ! units are problem specific
  PetscReal :: K_Aaq, K_Baq ! [mol/L]
  
  porosity = material_auxvar%porosity
  liquid_saturation = global_auxvar%sat(iphase)
  volume = material_auxvar%volume
  L_water = porosity*liquid_saturation*volume*1.d3  ! 1.d3 converts m^3 water -> L water
  
  Aaq = rt_auxvar%total(this%species_Aaq_id,iphase)
  Baq = rt_auxvar%total(this%species_Baq_id,iphase)
  Caq = rt_auxvar%total(this%species_Caq_id,iphase)
  Daq = rt_auxvar%total(this%species_Daq_id,iphase)
  Eaq = rt_auxvar%total(this%species_Eaq_id,iphase)
  Faq = rt_auxvar%total(this%species_Faq_id,iphase)
  
  Xim = rt_auxvar%immobile(this%species_Xim_id)
  Yim = rt_auxvar%immobile(this%species_Yim_id)

  ! initialize all rates to zero
  Rate = 0.d0

  RateA = 0.d0
  RateB = 0.d0
  RateC = 0.d0
  RateD = 0.d0
  RateE = 0.d0
  RateF = 0.d0
  RateX = 0.d0
  RateY = 0.d0

  ! stoichiometries
  ! reactants have negative stoichiometry
  ! products have positive stoichiometry
  stoichA = 0.d0
  stoichB = 0.d0
  stoichC = 0.d0
  stoichD = 0.d0
  stoichE = 0.d0
  stoichF = 0.d0
  stoichX = 0.d0
  stoichY = 0.d0
  
  ! kinetic rate constants
  k = 0.d0
  kr = 0.d0

  ! Monod half-saturation constants
  K_Aaq = 0.d0
  K_Baq = 0.d0

! PCL: A -> B -> C
! kinetic rate constants
! non-stationary state
! k1 = 1.d-1
! k2 = 5.d-2

! stationary state
! k1 = 1.d-1
! k2 = 1.d-0

! Rate1 = k1 * Aaq * L_water
! Rate2 = k2 * Baq * L_water

! RateA = -Rate1
! RateB =  Rate1 - Rate2
! RateC =  Rate2

! END PCL

  ! zero-order (A -> C)
  !k = 0.d0 ! WARNING: Too high a rate can result in negative concentrations
            !          which are non-physical and the code will not converge.
  !stoichA = -1.d0
  !stoichC = 1.d0
  !Rate = k * L_water
  !RateA = stoichA * Rate
  !RateC = stoichC * Rate
  
  ! first-order (A -> C)
  !k = 1.d-10
  !stoichA = -1.d0
  !stoichC = 1.d0
  !Rate = k * Aaq * L_water
  !RateA = stoichA * Rate
  !RateC = stoichC * Rate
  
  ! second-order (A + B -> C)
  !k = 1.d-10
  !stoichA = -1.d0
  !stoichB = -1.d0
  !stoichC = 1.d0
  !Rate = k * Aaq * Baq * L_water
  !RateA = stoichA * Rate
  !RateB = stoichB * Rate
  !RateC = stoichC * Rate
  
  ! Monod (A -> C)
  !k = 1.d-12
  !K_Aaq = 5.d-4
  !stoichA = -1.d0
  !stoichC = 1.d0
  !Rate = k * Aaq / (K_Aaq + Aaq) * L_water
  !RateA = stoichA * Rate
  !RateC = stoichC * Rate
  
  ! multiplicative Monod w/biomass
  ! A + 2B -> C
  !k = 1.d-7
  !K_Aaq = 5.d-4
  !K_Baq = 5.d-4
  !stoichA = -1.d0
  !stoichB = -2.d0
  !stoichC = 1.d0
  !Rate = k * Xim * Aaq / (K_Aaq + Aaq) * Baq / (K_Baq + Baq) * L_water
  !RateA = stoichA * Rate
  !RateB = stoichB * Rate
  !RateC = stoichC * Rate
  
  ! first-order forward - reverse (A <-> C)
  !k = 1.d-6
  !kr = 1.d-7
  !stoichA = -1.d0
  !stoichC = 1.d0
  !Rate = (k * Aaq - kr * Caq) * L_water
  !RateA = stoichA * Rate
  !RateC = stoichC * Rate
  
  ! NOTES
  ! 1. Always subtract contribution from residual
  ! 2. Units of residual are moles/second  
  Residual(this%species_Aaq_id) = Residual(this%species_Aaq_id) - RateA
  Residual(this%species_Baq_id) = Residual(this%species_Baq_id) - RateB
  Residual(this%species_Caq_id) = Residual(this%species_Caq_id) - RateC
  Residual(this%species_Daq_id) = Residual(this%species_Daq_id) - RateD
  Residual(this%species_Eaq_id) = Residual(this%species_Eaq_id) - RateE
  Residual(this%species_Faq_id) = Residual(this%species_Faq_id) - RateF
  Residual(this%species_Xim_id) = Residual(this%species_Xim_id) - RateX
  Residual(this%species_Yim_id) = Residual(this%species_Yim_id) - RateY
  
end subroutine SimpleReact

end module Reaction_Sandbox_Simple_class
