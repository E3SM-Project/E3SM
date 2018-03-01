module BeTR_LandunitType
  use bshr_kind_mod   , only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type, public :: betr_landunit_type

     ! topological mapping functionality
     integer , pointer :: itype     (:) => null() ! landunit type
     logical , pointer :: ifspecial (:) => null() ! true=>landunit is not vegetated
  end type betr_landunit_type

  type(betr_landunit_type), public :: betr_lun  !geomorphological landunits

end module BeTR_LandunitType
