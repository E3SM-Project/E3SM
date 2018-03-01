module BeTR_ColumnType
  !

  ! Column type for data transfer between lsm and BeTR
  ! !PUBLIC TYPES:
  use bshr_kind_mod   , only : r8 => shr_kind_r8
  implicit none
  save
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type, public :: betr_column_type
     ! g/l/c/p hierarchy, local g/l/c/p cells only
   integer , pointer :: snl      (:)   => null() ! number of snow layers
   real(r8), pointer :: dz       (:,:) => null() ! layer thickness (m)  (-nlevsno+1:nlevgrnd)
   real(r8), pointer :: zi       (:,:) => null() ! interface level below a "z" level (m) (-nlevsno+0:nlevgrnd)
   real(r8), pointer :: z        (:,:) => null() ! interface level below a "z" level (m) (-nlevsno+0:nlevgrnd)
   integer , pointer :: pfti     (:)   => null() ! beginning pft index for each column
   integer , pointer :: pftf     (:)   => null() ! ending pft index for each column
   integer , pointer :: npfts    (:)   => null() ! number of patches for each column
  contains
    procedure, public :: Init
  end type betr_column_type

  public :: create_betr_column_type
contains
  function create_betr_column_type()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(betr_column_type), pointer :: create_betr_column_type
    class(betr_column_type), pointer :: col

    allocate(col)
    create_betr_column_type => col

  end function create_betr_column_type
  !-------------------------------------------------------------------------------
  subroutine Init(this, bounds)
  use betr_varcon, only : ispval => bispval
  use betr_decompMod , only : betr_bounds_type
  implicit none
  class(betr_column_type), intent(inout) :: this
  type(betr_bounds_type), intent(in) :: bounds
  integer :: begc, endc
  integer :: lbj, ubj

  lbj = bounds%lbj; ubj = bounds%ubj
  begc = bounds%begc; endc = bounds%endc

  allocate(this%snl(begc:endc))          ; this%snl(:) = ispval
  allocate(this%dz(begc:endc, lbj:ubj))  ;
  allocate(this%zi(begc:endc,lbj-1:ubj)) ;
  allocate(this%z(begc:endc,lbj:ubj))    ;
  allocate(this%npfts(begc:endc))
  allocate(this%pfti(begc:endc))
  allocate(this%pftf(begc:endc))
  end subroutine Init
end module BeTR_ColumnType
