module BeTR_pftvarconType

implicit none
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type, public :: betr_pftvarcon_type
   integer :: nc3_arctic_grass       !value for C3 arctic grass
   integer :: nc3_nonarctic_grass    !value for C3 non-arctic grass
   integer :: nc4_grass              !value for C4 grass
   integer :: noveg                  !value for not vegetated
  contains
    procedure, public :: is_grass_patch
  end type betr_pftvarcon_type

  type(betr_pftvarcon_type), public :: betr_pftvarcon

contains
  function is_grass_patch(this, vtype)result(yesno)

  implicit none
  class(betr_pftvarcon_type), intent(in) :: this
  integer, intent(in) :: vtype

  logical :: yesno

  yesno = (vtype == this%nc3_arctic_grass) .or. &
          (vtype == this%nc3_nonarctic_grass) .or. &
          (vtype == this%nc4_grass)

  end function is_grass_patch
!------------------------------------------------------------------------
  function is_veg_patch(this, vtype)result(yesno)

  implicit none
  class(betr_pftvarcon_type), intent(in) :: this
  integer, intent(in) :: vtype
  logical :: yesno

  yesno = (vtype /=  this%noveg)
  end function is_veg_patch

end module BeTR_pftvarconType
