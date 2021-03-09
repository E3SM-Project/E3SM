module SpeciesMod

  !-----------------------------------------------------------------------
  ! Module holding information about different species available in the CN code (C, C13,
  ! C14, N).
  !
  !
  ! NOTE(wjs, 2016-06-05) Eventually I could imagine having a cn_species base class, with
  ! derived classes for each species type - so a cn_species_c class, a cn_species_c13
  ! class, a cn_species_c14 class and a cn_species_n class. These would contain methods
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

  integer, parameter, public :: CN_SPECIES_C12 = 1
  integer, parameter, public :: CN_SPECIES_C13 = 2
  integer, parameter, public :: CN_SPECIES_C14 = 3
  integer, parameter, public :: CN_SPECIES_N   = 4
  integer, parameter, public :: CN_SPECIES_P   = 5

  public :: species_from_string  ! convert a string representation to one of the constants defined here

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
       species = CN_SPECIES_C12
    case ('c13')
       species = CN_SPECIES_C13
    case ('c14')
       species = CN_SPECIES_C14
    case ('n')
       species = CN_SPECIES_N
    case ('p')
       species = CN_SPECIES_P
    end select

  end function species_from_string


end module SpeciesMod
