module Reaction_Microbial_Aux_module
  
  use Reaction_Database_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  PetscInt, parameter, public :: INHIBITION_THRESHOLD = 1
  PetscInt, parameter, public :: INHIBITION_THERMODYNAMIC = 2
  PetscInt, parameter, public :: INHIBITION_MONOD = 3
  PetscInt, parameter, public :: INHIBITION_INVERSE_MONOD = 4
  
  type, public :: microbial_rxn_type
    PetscInt :: id
    PetscInt :: itype
    character(len=MAXSTRINGLENGTH) :: reaction
    PetscReal :: rate_constant
    PetscReal :: activation_energy
    PetscBool :: print_me
    type(database_rxn_type), pointer :: dbaserxn    
    type(monod_type), pointer :: monod
    type(inhibition_type), pointer :: inhibition
    type(microbial_biomass_type), pointer :: biomass
    type(microbial_rxn_type), pointer :: next
  end type microbial_rxn_type
  
  type, public :: monod_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: species_name
    PetscReal :: half_saturation_constant
    PetscReal :: threshold_concentration
    type(monod_type), pointer :: next
  end type monod_type

  type, public :: inhibition_type
    PetscInt :: id
    PetscInt :: itype
    character(len=MAXWORDLENGTH) :: species_name
    PetscReal :: inhibition_constant
    PetscReal :: inhibition_constant2
    type(inhibition_type), pointer :: next
  end type inhibition_type

  type, public :: microbial_biomass_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: species_name
    PetscReal :: yield
  end type microbial_biomass_type
  
  type, public :: microbial_type

    PetscInt :: nrxn
    
    type(microbial_rxn_type), pointer :: microbial_rxn_list

    ! microbial reactions
    PetscReal, pointer :: rate_constant(:)
    PetscReal, pointer :: activation_energy(:)
    PetscReal, pointer :: stoich(:,:)
    PetscInt, pointer :: specid(:,:)
    PetscInt, pointer :: biomassid(:)
    PetscReal, pointer :: biomass_yield(:)
    PetscInt, pointer :: monodid(:,:)
    PetscInt, pointer :: inhibitionid(:,:)
    PetscInt, pointer :: monod_specid(:)
    PetscReal, pointer :: monod_K(:)
    PetscReal, pointer :: monod_Cth(:)
    PetscInt, pointer :: inhibition_type(:)
    PetscInt, pointer :: inhibition_specid(:)
    PetscReal, pointer :: inhibition_C(:)
    PetscReal, pointer :: inhibition_C2(:)
    
  end type microbial_type

  public :: MicrobialCreate, &
            MicrobialRxnCreate, &
            MicrobialMonodCreate, &
            MicrobialInhibitionCreate, &
            MicrobialBiomassCreate, &
            MicrobialGetMonodCount, &
            MicrobialGetInhibitionCount, &
            MicrobialGetBiomassCount, &
            MicrobialRxnDestroy, &
            MicrobialBiomassDestroy, &
            MicrobialDestroy
             
contains

! ************************************************************************** !

function MicrobialCreate()
  ! 
  ! Allocate and initialize microbial object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/12
  ! 

  implicit none
  
  type(microbial_type), pointer :: MicrobialCreate
  
  type(microbial_type), pointer :: microbial

  allocate(microbial)  
    
  nullify(microbial%microbial_rxn_list)
    
  microbial%nrxn = 0

  nullify(microbial%rate_constant)
  nullify(microbial%activation_energy)
  nullify(microbial%stoich)
  nullify(microbial%specid)
  nullify(microbial%biomassid)
  nullify(microbial%biomass_yield)
  nullify(microbial%monodid)
  nullify(microbial%inhibitionid)
  nullify(microbial%monod_specid)
  nullify(microbial%monod_K)
  nullify(microbial%monod_Cth)
  nullify(microbial%inhibition_type)
  nullify(microbial%inhibition_specid)
  nullify(microbial%inhibition_C)
  nullify(microbial%inhibition_C2)
  
  MicrobialCreate => microbial
  
end function MicrobialCreate

! ************************************************************************** !

function MicrobialRxnCreate()
  ! 
  ! Allocate and initialize a microbial object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/12
  ! 

  implicit none
  
  type(microbial_rxn_type), pointer :: MicrobialRxnCreate
  
  type(microbial_rxn_type), pointer :: microbial_rxn

  allocate(microbial_rxn)  
  microbial_rxn%id = 0
  microbial_rxn%itype = 0
  microbial_rxn%reaction = ''
  microbial_rxn%rate_constant = 0.d0
  microbial_rxn%activation_energy = 0.d0
  microbial_rxn%print_me = PETSC_FALSE
  nullify(microbial_rxn%biomass)
  nullify(microbial_rxn%dbaserxn)
  nullify(microbial_rxn%monod)
  nullify(microbial_rxn%inhibition)
  nullify(microbial_rxn%next)
  
  MicrobialRxnCreate => microbial_rxn
  
end function MicrobialRxnCreate

! ************************************************************************** !

function MicrobialMonodCreate()
  ! 
  ! Allocate and initialize a microbial monod object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/12
  ! 

  implicit none
  
  type(monod_type), pointer :: MicrobialMonodCreate
  
  type(monod_type), pointer :: monod

  allocate(monod)  
  monod%id = 0
  monod%species_name = ''
  monod%half_saturation_constant = 0.d0
  monod%threshold_concentration = 0.d0
  nullify(monod%next)
  
  MicrobialMonodCreate => monod
  
end function MicrobialMonodCreate

! ************************************************************************** !

function MicrobialInhibitionCreate()
  ! 
  ! Allocate and initialize a microbial inhibition
  ! object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/12
  ! 

  implicit none
  
  type(inhibition_type), pointer :: MicrobialInhibitionCreate
  
  type(inhibition_type), pointer :: inhibition
  
  allocate(inhibition)  
  inhibition%id = 0
  inhibition%itype = 0
  inhibition%species_name = ''
  inhibition%inhibition_constant = UNINITIALIZED_DOUBLE
  inhibition%inhibition_constant2 = 0.d0
  nullify(inhibition%next)
  
  MicrobialInhibitionCreate => inhibition
  
end function MicrobialInhibitionCreate

! ************************************************************************** !

function MicrobialBiomassCreate()
  ! 
  ! Allocate and initialize a microbial biomass object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/02/13
  ! 

  implicit none
  
  type(microbial_biomass_type), pointer :: MicrobialBiomassCreate
  
  type(microbial_biomass_type), pointer :: biomass

  allocate(biomass)  
  biomass%id = 0
  biomass%species_name = ''
  biomass%yield = 0.d0
  
  MicrobialBiomassCreate => biomass
  
end function MicrobialBiomassCreate

! ************************************************************************** !

function MicrobialGetMonodCount(microbial_rxn)
  ! 
  ! Counts number of monod expressions in
  ! microbial reaction
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/12
  ! 

  implicit none
  
  type(microbial_rxn_type) :: microbial_rxn
  
  PetscInt :: MicrobialGetMonodCount
  
  type(monod_type), pointer :: cur_monod
  PetscInt :: icount
  
  icount = 0
  cur_monod => microbial_rxn%monod
  do
    if (.not.associated(cur_monod)) exit
    icount = icount + 1
    cur_monod => cur_monod%next
  enddo
  
  MicrobialGetMonodCount = icount
  
end function MicrobialGetMonodCount

! ************************************************************************** !

function MicrobialGetInhibitionCount(microbial_rxn)
  ! 
  ! Counts number of inhibiton expressions in
  ! microbial reaction
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/12
  ! 

  implicit none
  
  type(microbial_rxn_type) :: microbial_rxn
  
  PetscInt :: MicrobialGetInhibitionCount
  
  type(inhibition_type), pointer :: cur_inhibition
  PetscInt :: icount
  
  icount = 0
  cur_inhibition => microbial_rxn%inhibition
  do
    if (.not.associated(cur_inhibition)) exit
    icount = icount + 1
    cur_inhibition => cur_inhibition%next
  enddo
  
  MicrobialGetInhibitionCount = icount
  
end function MicrobialGetInhibitionCount

! ************************************************************************** !

function MicrobialGetBiomassCount(microbial)
  ! 
  ! Returns the number of biomass species
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/02/13
  ! 

  implicit none
  
  PetscInt :: MicrobialGetBiomassCount
  type(microbial_type) :: microbial

  type(microbial_rxn_type), pointer :: microbial_rxn

  MicrobialGetBiomassCount = 0
  microbial_rxn => microbial%microbial_rxn_list
  do
    if (.not.associated(microbial_rxn%biomass)) exit
    MicrobialGetBiomassCount = MicrobialGetBiomassCount + 1
    microbial_rxn => microbial_rxn%next
  enddo

end function MicrobialGetBiomassCount

! ************************************************************************** !

subroutine MicrobialRxnDestroy(microbial)
  ! 
  ! Deallocates a microbial rxn object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/12
  ! 

  implicit none
    
  type(microbial_rxn_type), pointer :: microbial

  call DatabaseRxnDestroy(microbial%dbaserxn)
  call MicrobialMonodDestroy(microbial%monod)
  call MicrobialInhibitionDestroy(microbial%inhibition)
  call MicrobialBiomassDestroy(microbial%biomass)

  deallocate(microbial)  
  nullify(microbial)

end subroutine MicrobialRxnDestroy

! ************************************************************************** !

recursive subroutine MicrobialMonodDestroy(monod)
  ! 
  ! Deallocates a microbial monod object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/12
  ! 

  implicit none
    
  type(monod_type), pointer :: monod

  if (.not.associated(monod)) return
  
  call MicrobialMonodDestroy(monod%next)
  
  deallocate(monod)
  nullify(monod)
  
end subroutine MicrobialMonodDestroy

! ************************************************************************** !

recursive subroutine MicrobialInhibitionDestroy(inhibition)
  ! 
  ! Deallocates a microbial inhibition object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/12
  ! 

  implicit none
    
  type(inhibition_type), pointer :: inhibition

  if (.not. associated(inhibition)) return

  call MicrobialInhibitionDestroy(inhibition%next)
  
  deallocate(inhibition)
  nullify(inhibition)

end subroutine MicrobialInhibitionDestroy

! ************************************************************************** !

subroutine MicrobialBiomassDestroy(biomass)
  ! 
  ! Deallocates a microbial biomass object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/02/13
  ! 

  implicit none
    
  type(microbial_biomass_type), pointer :: biomass

  if (.not.associated(biomass)) return
  
  deallocate(biomass)
  nullify(biomass)
  
end subroutine MicrobialBiomassDestroy

! ************************************************************************** !

subroutine MicrobialDestroy(microbial)
  ! 
  ! Deallocates a microbial object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/29/08
  ! 

  use Utility_module, only: DeallocateArray
  
  implicit none

  type(microbial_type), pointer :: microbial
  
  type(microbial_rxn_type), pointer :: cur_microbial, prev_microbial

  if (.not.associated(microbial)) return
  
  ! microbial reactions
  cur_microbial => microbial%microbial_rxn_list
  do
    if (.not.associated(cur_microbial)) exit
    prev_microbial => cur_microbial
    cur_microbial => cur_microbial%next
    call MicrobialRxnDestroy(prev_microbial)
  enddo    
  nullify(microbial%microbial_rxn_list)
  
  call DeallocateArray(microbial%rate_constant)
  call DeallocateArray(microbial%activation_energy)
  call DeallocateArray(microbial%stoich)
  call DeallocateArray(microbial%specid)
  call DeallocateArray(microbial%biomassid)
  call DeallocateArray(microbial%biomass_yield)
  call DeallocateArray(microbial%monodid)
  call DeallocateArray(microbial%inhibitionid)
  call DeallocateArray(microbial%monod_specid)
  call DeallocateArray(microbial%monod_K)
  call DeallocateArray(microbial%monod_Cth)
  call DeallocateArray(microbial%inhibition_type)
  call DeallocateArray(microbial%inhibition_specid)
  call DeallocateArray(microbial%inhibition_C)
  call DeallocateArray(microbial%inhibition_C2)
  
  deallocate(microbial)
  nullify(microbial)

end subroutine MicrobialDestroy

end module Reaction_Microbial_Aux_module
