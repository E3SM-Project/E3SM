module BeTR_PatchType
  use bshr_kind_mod, only: r8 => shr_kind_r8
implicit none
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type, public :: betr_patch_type
   real(r8), pointer :: wtcol(:)  => null()  !weight relative to a column
   integer,  pointer :: column(:) => null() !column index
   integer,  pointer :: itype(:) => null()  !patch vegetation
   real(r8), pointer :: crop(:)  => null()  !crop pft: 0. = not crop, 1. = crop pft
   integer :: npfts
  contains
    procedure, public :: init
  end type betr_patch_type

  public :: create_betr_patch_type
contains
  function create_betr_patch_type()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(betr_patch_type), pointer :: create_betr_patch_type
    class(betr_patch_type), pointer :: pft

    allocate(pft)
    create_betr_patch_type => pft

  end function create_betr_patch_type
  !-------------------------------------------------------------------------------
  subroutine Init(this, bounds)
  use betr_varcon         , only : ispval => bispval
  use bshr_infnan_mod     , only : nan => shr_infnan_nan, assignment(=)
  use betr_decompMod      , only : betr_bounds_type
  implicit none
  class(betr_patch_type), intent(inout) :: this
  type(betr_bounds_type), intent(in) :: bounds
  integer :: begp, endp
  integer :: lbj, ubj

  begp = bounds%begp; endp = bounds%endp
  allocate(this%wtcol(begp:endp)); this%wtcol(:) = nan
  allocate(this%column(begp:endp)); this%column(:) = 1
  allocate(this%itype(begp:endp)); this%itype(:) = ispval
  allocate(this%crop(begp:endp)); this%crop(:) = nan
  end subroutine Init
end module BeTR_PatchType
