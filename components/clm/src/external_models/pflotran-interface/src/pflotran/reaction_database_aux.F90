module Reaction_Database_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  type, public :: database_rxn_type
    PetscInt :: nspec
    character(len=MAXWORDLENGTH), pointer :: spec_name(:)
    PetscReal, pointer :: stoich(:)
    PetscInt, pointer :: spec_ids(:)
    PetscReal, pointer :: logK(:)
    PetscReal, pointer :: logKCoeff_hpt(:)
  end type database_rxn_type

  public :: BasisAlignSpeciesInRxn, &
            BasisSubSpeciesInGasOrSecRxn, &
            BasisSubSpeciesinMineralRxn, &
            DatabaseRxnCreate, &
            DatabaseRxnCreateFromRxnString, &
            DatabaseCheckLegitimateLogKs, &
            DatabaseRxnDestroy
            
contains

! ************************************************************************** !

function DatabaseRxnCreate()
  ! 
  ! Allocate and initialize an equilibrium reaction
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/01/08
  ! 

  implicit none
    
  type(database_rxn_type), pointer :: DatabaseRxnCreate

  type(database_rxn_type), pointer :: dbaserxn

  allocate(dbaserxn)
  dbaserxn%nspec = 0
  nullify(dbaserxn%spec_name)
  nullify(dbaserxn%stoich)
  nullify(dbaserxn%spec_ids)
  nullify(dbaserxn%logK)
  
  DatabaseRxnCreate => dbaserxn
  
end function DatabaseRxnCreate

! ************************************************************************** !

function DatabaseRxnCreateFromRxnString(reaction_string, &
                                        naqcomp, aq_offset, &
                                        primary_aq_species_names, &
                                        nimcomp, im_offset, &
                                        primary_im_species_names, &
                                        consider_immobile_species,&
                                        option)
  ! 
  ! Creates a database reaction given a
  ! reaction string
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/12
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use String_module
  use Input_Aux_module
  
  implicit none
  
  character(len=MAXSTRINGLENGTH) :: reaction_string
  PetscInt :: naqcomp ! mobile aqueoues species
  PetscInt :: aq_offset ! offset for aqueous species
  character(len=MAXWORDLENGTH) :: primary_aq_species_names(naqcomp)
  PetscInt :: nimcomp ! immobile primary speces (e.g. biomass)
  PetscInt :: im_offset ! offset for aqueous species
  character(len=MAXWORDLENGTH), pointer :: primary_im_species_names(:)
  PetscBool :: consider_immobile_species
  type(option_type) :: option
    
  type(database_rxn_type), pointer :: DatabaseRxnCreateFromRxnString

  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=MAXWORDLENGTH) :: word, word2
  PetscInt :: icount
  PetscInt :: midpoint
  PetscInt :: i, j, idum
  PetscReal :: value
  PetscBool :: negative_flag
  PetscBool :: found
  PetscErrorCode :: ierr
  type(database_rxn_type), pointer :: dbaserxn
  
  
  dbaserxn => DatabaseRxnCreate()

  icount = 0
  ! Be sure to copy as words are removed when read.  Need to full string for 
  ! later below
  string = reaction_string
  do
    ierr = 0
    call InputReadWord(string,word,PETSC_TRUE,ierr)
    if (InputError(ierr)) exit

    select case(word)
      case('+')
      case('-')
      case('=','<=>','<->','->','=>')
      case default
      ! try reading as double precision
      string2 = word
      if (.not.StringStartsWithAlpha(string2) .and. &
          StringIntegerDoubleOrWord(string2) /= STRING_IS_A_WORD) then
        ! the word is the stoichiometry value
      else
        ! check water
        word2 = 'H2O'
        if (.not.StringCompareIgnoreCase(word,word2)) then
          ! the word is the species name
          icount = icount + 1
        endif
      endif
    end select

  enddo
      
  ! load species into database format
  dbaserxn%nspec = icount
  allocate(dbaserxn%spec_name(icount))
  dbaserxn%spec_name = ''
  allocate(dbaserxn%stoich(icount))
  dbaserxn%stoich = UNINITIALIZED_DOUBLE
  allocate(dbaserxn%spec_ids(icount))
  dbaserxn%spec_ids = 0

  string = reaction_string
  icount = 1
  ! midpoint points to the first product species, as in
  ! reactant1 + reactant2 <-> product1 + product2
  midpoint = 0
  negative_flag = PETSC_FALSE
  do
    !geh: This conditional ensures that if water is at the end of
    !     the reaction expression, it is skipped.
    if (icount > dbaserxn%nspec) exit
        
    ierr = 0
    call InputReadWord(string,word,PETSC_TRUE,ierr)
    if (InputError(ierr)) exit

    select case(word)
      case('+')
      case('-')
        ! toggle negative flag
        if (negative_flag) then
          negative_flag = PETSC_FALSE
        else
          negative_flag = PETSC_TRUE
        endif
      case('=','<=>','<->','->','=>')
        midpoint = icount
      case default
        ! try reading as double precision
        string2 = word
        if (.not.StringStartsWithAlpha(string2) .and. &
            StringIntegerDoubleOrWord(string2) /= STRING_IS_A_WORD) then
          ! negate if a product
          call InputReadDouble(string2,option,value,ierr)
          if (ierr /= 0) then
            option%io_buffer = 'Keyword "' // trim(word) // &
               '" not recognized in reaction string "' // &
               trim(reaction_string) // '".'
            call printErrMsg(option)
          endif
          ! negate if negative stoichiometry
          if (negative_flag) value = -1.0*value
          dbaserxn%stoich(icount) = value
        else
          dbaserxn%spec_name(icount) = word
          if (negative_flag .and. &
              (dbaserxn%stoich(icount) + 999.d0) < 1.d-10) then
            dbaserxn%stoich(icount) = -1.d0
          endif

          ! set the primary aqueous species id
          found = PETSC_FALSE
          do i = 1, naqcomp
            if (StringCompare(word,primary_aq_species_names(i), &
                              MAXWORDLENGTH)) then
              dbaserxn%spec_ids(icount) = i + aq_offset
              found = PETSC_TRUE
              exit      
            endif
          enddo
          ! set the primary immobile species id
          if (.not.found .and. consider_immobile_species) then
            do i = 1, nimcomp
              if (StringCompare(word,primary_im_species_names(i), &
                                MAXWORDLENGTH)) then
                dbaserxn%spec_ids(icount) = i + im_offset
                found = PETSC_TRUE
                exit      
              endif
            enddo
          endif
          
          ! check water
          word2 = 'H2O'
          if (StringCompareIgnoreCase(word,word2)) then
            ! set stoichiometry back to uninitialized
            dbaserxn%stoich(icount) = UNINITIALIZED_DOUBLE
            ! don't increment icount
          else if (.not.found) then
            if (consider_immobile_species) then
              option%io_buffer = 'Species ' // trim(word) // &
                        ' in reaction not found among primary species list.'
            else
              option%io_buffer = 'Species ' // trim(word) // &
                ' in reaction not found among primary aqueous species list.'
            endif
            call printErrMsg(option)     
          else
            icount = icount + 1
          endif
        endif
        negative_flag = PETSC_FALSE
    end select
  enddo
      
  ! if no stoichiometry specified, default = 1.
  do i = 1, dbaserxn%nspec
    if ((dbaserxn%stoich(i) + 999.d0) < 1.d-10) dbaserxn%stoich(i) = 1.d0
  enddo
  if (midpoint > 0) then
    ! negate stoichiometries after midpoint
    do i = midpoint, dbaserxn%nspec
      dbaserxn%stoich(i) = -1.d0*dbaserxn%stoich(i)
    enddo
  endif
  ! now negate all stoichiometries to have - for reactants; + for products
  do i = 1, dbaserxn%nspec
    dbaserxn%stoich(i) = -1.d0*dbaserxn%stoich(i)
  enddo
  ! reorder species ids in ascending order
  do i = 1, dbaserxn%nspec
    do j = i+1, dbaserxn%nspec
      if (dbaserxn%spec_ids(i) > dbaserxn%spec_ids(j)) then
        ! swap ids
        idum = dbaserxn%spec_ids(j)
        dbaserxn%spec_ids(j) = dbaserxn%spec_ids(i)
        dbaserxn%spec_ids(i) = idum
        ! swap stoichiometry
        value = dbaserxn%stoich(j)
        dbaserxn%stoich(j) = dbaserxn%stoich(i)
        dbaserxn%stoich(i) = value
        ! swap names
        word = dbaserxn%spec_name(j)
        dbaserxn%spec_name(j) = dbaserxn%spec_name(i)
        dbaserxn%spec_name(i) = word
      endif
    enddo
  enddo
  
  DatabaseRxnCreateFromRxnString => dbaserxn
  
end function DatabaseRxnCreateFromRxnString

! ************************************************************************** !

subroutine BasisAlignSpeciesInRxn(num_basis_species,basis_names, &
                                  num_rxn_species,rxn_species_names, &
                                  rxn_stoich,rxn_species_ids,species_name, &
                                  option)
  ! 
  ! Aligns the ordering of species in reaction with
  ! the current basis
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/07/08
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use String_module
  
  implicit none
  
  PetscInt :: num_basis_species
  character(len=MAXWORDLENGTH) :: basis_names(num_basis_species), species_name
  PetscInt :: num_rxn_species
  character(len=MAXWORDLENGTH) :: rxn_species_names(num_rxn_species)
  PetscReal :: rxn_stoich(num_rxn_species)
  PetscInt :: rxn_species_ids(num_rxn_species)
  type(option_type) :: option

  PetscInt :: i_rxn_species
  PetscInt :: i_basis_species
  PetscReal :: stoich_new(num_basis_species)
  PetscBool :: found
  
  stoich_new = 0.d0
  do i_rxn_species = 1, num_rxn_species
    found = PETSC_FALSE
    do i_basis_species = 1, num_basis_species
      if (StringCompare(rxn_species_names(i_rxn_species), &
                            basis_names(i_basis_species), &
                            MAXWORDLENGTH)) then
        stoich_new(i_basis_species) = rxn_stoich(i_rxn_species)
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      option%io_buffer = trim(rxn_species_names(i_rxn_species)) // &
               ' not found in basis (BasisAlignSpeciesInRxn) for species ' // &
               trim(species_name)
      call printErrMsg(option)
    endif
  enddo
  
  ! zero everthing out
  rxn_species_names = ''
  rxn_stoich = 0.d0
  rxn_species_ids = 0

  ! fill in
  i_rxn_species = 0
  do i_basis_species = 1, num_basis_species
    if (dabs(stoich_new(i_basis_species)) > 1.d-40) then
      i_rxn_species = i_rxn_species + 1
      rxn_species_names(i_rxn_species) = basis_names(i_basis_species)
      rxn_stoich(i_rxn_species) = stoich_new(i_basis_species)
      rxn_species_ids(i_rxn_species) = i_basis_species
    endif
  enddo
  
  if (i_rxn_species /= num_rxn_species) then
    write(option%io_buffer,*) &
                   'Number of reaction species does not match original:', &
                    i_rxn_species, num_rxn_species
    call printErrMsg(option)
  endif

end subroutine BasisAlignSpeciesInRxn 

! ************************************************************************** !

subroutine BasisSubSpeciesInGasOrSecRxn(name1,dbaserxn1,dbaserxn2,scale)
  ! 
  ! Swaps out a chemical species in a chemical
  ! reaction, replacing it with the species in a
  ! secondary reaction (swaps 1 into 2)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/06/08
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use String_module

  implicit none
  
  character(len=MAXWORDLENGTH) :: name1
  type(database_rxn_type) :: dbaserxn1
  type(database_rxn_type) :: dbaserxn2
  PetscReal :: scale
  
  PetscInt :: i, j, tempcount, prevcount
  character(len=MAXWORDLENGTH) :: tempnames(20)
  PetscReal :: tempstoich(20)
  PetscBool :: found

  tempnames = ''
  tempstoich = 0.d0

  ! load species in reaction other than species 1 into new arrays
  scale = 1.d0
  tempcount = 0
  do i=1,dbaserxn2%nspec
    if (.not.StringCompare(name1, &
                             dbaserxn2%spec_name(i), &
                             MAXWORDLENGTH)) then
      tempcount = tempcount + 1
      tempnames(tempcount) = dbaserxn2%spec_name(i)
      tempstoich(tempcount) = dbaserxn2%stoich(i)
    else
      scale = dbaserxn2%stoich(i)
    endif
  enddo
  
  ! search for duplicate species and add stoichs or add new species
  ! if not duplicated
  do j=1,dbaserxn1%nspec
    found = PETSC_FALSE
    do i=1,tempcount
      if (StringCompare(tempnames(i), &
                          dbaserxn1%spec_name(j), &
                          MAXWORDLENGTH)) then
        tempstoich(i) = tempstoich(i) + scale*dbaserxn1%stoich(j)
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      tempcount = tempcount + 1
      tempnames(tempcount) = dbaserxn1%spec_name(j)
      tempstoich(tempcount) = scale*dbaserxn1%stoich(j)
    endif
  enddo
  
  ! deallocate arrays
  deallocate(dbaserxn2%spec_name)
  deallocate(dbaserxn2%stoich)
  
  ! check for zero stoichiometries due to cancelation
  prevcount = tempcount
  tempcount = 0
  do i=1,prevcount
    if (dabs(tempstoich(i)) > 1.d-10) then
      tempcount = tempcount + 1
      tempnames(tempcount) = tempnames(i)
      tempstoich(tempcount) = tempstoich(i)
    endif
  enddo
  
  tempnames(tempcount+1:) = ''
  tempstoich(tempcount+1:) = 0.d0
  
  ! reallocate
  allocate(dbaserxn2%spec_name(tempcount))
  allocate(dbaserxn2%stoich(tempcount))
  
  ! fill arrays in dbaserxn
  dbaserxn2%nspec = tempcount
  do i=1,tempcount
    dbaserxn2%spec_name(i) = tempnames(i)
    dbaserxn2%stoich(i) = tempstoich(i)
  enddo
  
  dbaserxn2%logK = dbaserxn2%logK + scale*dbaserxn1%logK

end subroutine BasisSubSpeciesInGasOrSecRxn

! ************************************************************************** !

subroutine BasisSubSpeciesInMineralRxn(name,sec_dbaserxn,mnrl_dbaserxn,scale)
  ! 
  ! Swaps out a chemical species in a chemical
  ! reaction, replacing it with the species in a
  ! secondary reaction (swaps 1 into 2)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/06/08
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use String_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: name
  type(database_rxn_type) :: sec_dbaserxn
  type(database_rxn_type) :: mnrl_dbaserxn
  PetscReal :: scale
  
  PetscInt :: i, j, tempcount, prevcount
  character(len=MAXWORDLENGTH) :: tempnames(20)
  PetscReal :: tempstoich(20)
  PetscBool :: found

  tempnames = ''
  tempstoich = 0.d0
  
  ! load species in reaction other than species 1 into new arrays
  scale = 1.d0
  tempcount = 0
  do i=1,mnrl_dbaserxn%nspec
    if (.not.StringCompare(name, &
                           mnrl_dbaserxn%spec_name(i), &
                           MAXWORDLENGTH)) then
      tempcount = tempcount + 1
      tempnames(tempcount) = mnrl_dbaserxn%spec_name(i)
      tempstoich(tempcount) = mnrl_dbaserxn%stoich(i)
    else
      scale = mnrl_dbaserxn%stoich(i)
    endif
  enddo
  
  ! search for duplicate species and add stoichs or add new species
  ! if not duplicated
  do j=1,sec_dbaserxn%nspec
    found = PETSC_FALSE
    do i=1,mnrl_dbaserxn%nspec
      if (StringCompare(tempnames(i), &
                        sec_dbaserxn%spec_name(j), &
                        MAXWORDLENGTH)) then
        tempstoich(i) = tempstoich(i) + scale*sec_dbaserxn%stoich(j)
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      tempcount = tempcount + 1
      tempnames(tempcount) = sec_dbaserxn%spec_name(j)
      tempstoich(tempcount) = scale*sec_dbaserxn%stoich(j)
    endif
  enddo
  
  ! deallocate arrays
  deallocate(mnrl_dbaserxn%spec_name)
  deallocate(mnrl_dbaserxn%stoich)
  
  ! check for zero stoichiometries due to cancelation
  prevcount = tempcount
  tempcount = 0
  do i=1,prevcount
    if (dabs(tempstoich(i)) > 1.d-10) then
      tempcount = tempcount + 1
      tempnames(tempcount) = tempnames(i)
      tempstoich(tempcount) = tempstoich(i)
    endif
  enddo

  tempnames(tempcount+1:) = ''
  tempstoich(tempcount+1:) = 0.d0
    
  ! reallocate
  allocate(mnrl_dbaserxn%spec_name(tempcount))
  allocate(mnrl_dbaserxn%stoich(tempcount))
  
  ! fill arrays in dbaserxn
  mnrl_dbaserxn%nspec = tempcount
  do i=1,tempcount
    mnrl_dbaserxn%spec_name(i) = tempnames(i)
    mnrl_dbaserxn%stoich(i) = tempstoich(i)
  enddo
  
  mnrl_dbaserxn%logK = mnrl_dbaserxn%logK + scale*sec_dbaserxn%logK

end subroutine BasisSubSpeciesInMineralRxn

! ************************************************************************** !

function DatabaseCheckLegitimateLogKs(dbaserxn,species_name,temperatures, &
                                      option)
  ! 
  ! Checks whether legitimate log Ks exist for
  ! all database temperatures if running
  ! non-isothermal
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/07/13
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  
  implicit none
    
  type(database_rxn_type), pointer :: dbaserxn
  character(len=MAXWORDLENGTH) :: species_name
  PetscReal :: temperatures(:)
  type(option_type) :: option

  PetscBool :: DatabaseCheckLegitimateLogKs
  
  PetscInt :: itemp
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  
  DatabaseCheckLegitimateLogKs = PETSC_TRUE

  if (.not.associated(dbaserxn) .or. option%use_isothermal) return
  
  string = ''
  do itemp = 1, size(dbaserxn%logK)
    if (dabs(dbaserxn%logK(itemp) - 500.) < 1.d-10) then
      write(word,'(f5.1)') temperatures(itemp)
      string = trim(string) // ' ' // word
      DatabaseCheckLegitimateLogKs = PETSC_FALSE
    endif
  enddo
  
  if (.not.DatabaseCheckLegitimateLogKs) then
    option%io_buffer = 'Undefined log Ks for temperatures (' // &
                       trim(adjustl(string)) // ') for species "' // &
                       trim(species_name) // '" in database.'
    call printWrnMsg(option)
  endif
  
end function DatabaseCheckLegitimateLogKs

! ************************************************************************** !

subroutine DatabaseRxnDestroy(dbaserxn)
  ! 
  ! Deallocates a database reaction
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/29/08
  ! 
  use Utility_module, only : DeallocateArray
  
  implicit none
    
  type(database_rxn_type), pointer :: dbaserxn

  if (.not.associated(dbaserxn)) return
  
  if (associated(dbaserxn%spec_name)) deallocate(dbaserxn%spec_name)
  nullify(dbaserxn%spec_name)
  call DeallocateArray(dbaserxn%spec_ids)
  call DeallocateArray(dbaserxn%stoich)
  call DeallocateArray(dbaserxn%logK)

  deallocate(dbaserxn)  
  nullify(dbaserxn)

end subroutine DatabaseRxnDestroy

end module Reaction_Database_Aux_module
