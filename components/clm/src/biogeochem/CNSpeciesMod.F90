module CNSpeciesMod

  !-----------------------------------------------------------------------
  ! Module holding information about different species available in the CN code (C, C13,
  ! C14, N).
  !
  !
  ! NOTE(wjs, 2016-06-05) Eventually I could imagine having a cn_species base class, with
  ! derived classes for each species type - so a NUTRIENT_SPECIES_c class, a NUTRIENT_SPECIES_c13
  ! class, a NUTRIENT_SPECIES_c14 class and a NUTRIENT_SPECIES_n class. These would contain methods
  ! to handle calculations specific to each species type. For example, there could be a
  ! carbon_multiplier method that returns the species-specific multiplier that you would
  ! apply to a variable in units of gC/m2 to give you g[this species]/m2 (this would
  ! depend on pft type).
  !
  ! Basically, anywhere where there is code that has a conditional based on the constants
  ! defined here, we could replace that with polymorphism using a cn_species class.
  !
  ! Eventually I think it would make sense to make this contain an instance of
  ! species_base_type (i.e., the class used to determine history & restart field names),
  ! with forwarding methods. So then (e.g.) a cn_products_type object would just contain a
  ! cn_species object (which in turn would contain a species_metadata [or whatever we call
  ! it] object).

  implicit none
  private

  integer, parameter, public :: NUTRIENT_SPECIES_C12 = 1
  integer, parameter, public :: NUTRIENT_SPECIES_C13 = 2
  integer, parameter, public :: NUTRIENT_SPECIES_C14 = 3
  integer, parameter, public :: NUTRIENT_SPECIES_N   = 4
  integer, parameter, public :: NUTRIENT_SPECIES_P   = 5

  public :: species_from_string      ! convert a string representation to one of the constants defined here
  public :: species_name_from_string ! convert a string representation to nutrient name
  public :: species_history_name_suffix_from_string

contains

  !-----------------------------------------------------------------------
  function species_from_string(species_string) result(species)
    !
    ! !DESCRIPTION:
    ! Convert a string representation to one of the constants defined here
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer :: species  ! function result
    character(len=*), intent(in) :: species_string  ! string representation of species (should be lowercase)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'species_from_string'
    !-----------------------------------------------------------------------

    select case (species_string)
    case ('c12')
       species = NUTRIENT_SPECIES_C12
    case ('c13')
       species = NUTRIENT_SPECIES_C13
    case ('c14')
       species = NUTRIENT_SPECIES_C14
    case ('n')
       species = NUTRIENT_SPECIES_N
    case ('p')
       species = NUTRIENT_SPECIES_P
    end select

  end function species_from_string

  !-----------------------------------------------------------------------
  function species_name_from_string(species_string) result(species_name)
    !
    ! !DESCRIPTION:
    ! Convert a string representation to species name
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: species_string  ! string representation of species (should be lowercase)
    !
    ! !LOCAL VARIABLES:
    character(len=3)             :: species_name

    character(len=*), parameter :: subname = 'species_name_from_string'
    !-----------------------------------------------------------------------

    select case (species_string)
    case ('c12')
       species_name = 'C'
    case ('c13')
       species_name = 'C13'
    case ('c14')
       species_name = 'C14'
    case ('n')
       species_name = 'N'
    case ('p')
       species_name = 'P'
    end select

  end function species_name_from_string

  !-----------------------------------------------------------------------
  function species_history_name_suffix_from_string(species_string) result(name_suffix)
    !
    ! !DESCRIPTION:
    ! Convert a string representation to species name
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: species_string  ! string representation of species (should be lowercase)
    !
    ! !LOCAL VARIABLES:
    character(len=3)             :: name_suffix

    character(len=*), parameter :: subname = 'species_name_from_string'
    !-----------------------------------------------------------------------

    select case (species_string)
    case ('c12')
       name_suffix = ''
    case ('c13')
       name_suffix = 'C13_'
    case ('c14')
       name_suffix = 'C14_'
    case ('n')
       name_suffix = ''
    case ('p')
       name_suffix = ''
    end select

  end function species_history_name_suffix_from_string

end module CNSpeciesMod
