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
  public :: nutrient_name
  public :: nutrient_units
  public :: nutrient_text
  public :: nutrient_restart_name
  public :: by_species
  public :: hist_default
  public :: is_carbon            ! whether a species is a carbon isotope

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


  pure logical function is_carbon(species)
    ! Whether species is one of the carbon isotopes
    integer, intent(in) :: species
    is_carbon = species == CN_SPECIES_C12 .or. species == CN_SPECIES_C13 .or. species == CN_SPECIES_C14
  end function is_carbon

  pure function species_label(species) result(label)
    ! C, C13, C14, N or P
    integer, intent(in) :: species
    character(len=:), allocatable :: label
    select case (species)
    case (CN_SPECIES_C13)
       label = 'C13'
    case (CN_SPECIES_C14)
       label = 'C14'
    case (CN_SPECIES_N)
       label = 'N'
    case (CN_SPECIES_P)
       label = 'P'
    case default
       label = 'C'
    end select
  end function species_label

  pure function fill_template(template, with) result(text)
    ! template with each '#' replaced by `with`
    character(len=*), intent(in) :: template, with
    character(len=:), allocatable :: text
    integer :: i
    text = ''
    do i = 1, len(template)
       if (template(i:i) == '#') then
          text = text // with
       else
          text = text // template(i:i)
       end if
    end do
  end function fill_template

  pure function nutrient_name(species, template) result(name)
    ! history field name: '#' is the element letter, carbon isotopes prefixed (C13_LEAF#)
    integer, intent(in) :: species
    character(len=*), intent(in) :: template
    character(len=:), allocatable :: name
    if (species == CN_SPECIES_C13 .or. species == CN_SPECIES_C14) then
       name = species_label(species) // '_' // fill_template(template, 'C')
    else
       name = fill_template(template, species_label(species))
    end if
  end function nutrient_name

  pure function nutrient_units(species, template) result(units)
    ! units: '#' is the element (or isotope) label (g#/m^2 -> gC13/m^2)
    integer, intent(in) :: species
    character(len=*), intent(in) :: template
    character(len=:), allocatable :: units
    units = fill_template(template, species_label(species))
  end function nutrient_units

  pure function nutrient_text(species, template) result(text)
    ! long name: '#' is the element letter, carbon isotopes prefixed (C13 leaf C)
    integer, intent(in) :: species
    character(len=*), intent(in) :: template
    character(len=:), allocatable :: text
    if (species == CN_SPECIES_C13 .or. species == CN_SPECIES_C14) then
       text = species_label(species) // ' ' // fill_template(template, 'C')
    else
       text = fill_template(template, species_label(species))
    end if
  end function nutrient_text

  pure function nutrient_restart_name(species, template) result(name)
    ! restart variable name: '#' is the element letter in lowercase, carbon isotopes suffixed (leaf#_storage -> leafc_storage_13)
    integer, intent(in) :: species
    character(len=*), intent(in) :: template
    character(len=:), allocatable :: name
    select case (species)
    case (CN_SPECIES_C13)
       name = fill_template(template, 'c') // '_13'
    case (CN_SPECIES_C14)
       name = fill_template(template, 'c') // '_14'
    case (CN_SPECIES_N)
       name = fill_template(template, 'n')
    case (CN_SPECIES_P)
       name = fill_template(template, 'p')
    case default
       name = fill_template(template, 'c')
    end select
  end function nutrient_restart_name

  pure function by_species(species, c12, c13, c14, n, p) result(text)
    ! the text given for species
    integer, intent(in) :: species
    character(len=*), intent(in) :: c12, c13, c14, n, p
    character(len=:), allocatable :: text
    select case (species)
    case (CN_SPECIES_C13)
       text = c13
    case (CN_SPECIES_C14)
       text = c14
    case (CN_SPECIES_N)
       text = n
    case (CN_SPECIES_P)
       text = p
    case default
       text = c12
    end select
  end function by_species

  pure function hist_default(active) result(default)
    ! a history field's default= argument
    logical, intent(in) :: active
    character(len=:), allocatable :: default
    if (active) then
       default = 'active'
    else
       default = 'inactive'
    end if
  end function hist_default

end module SpeciesMod
