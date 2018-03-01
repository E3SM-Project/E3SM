module ColumnType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Column data type allocation and initialization
  ! --------------------------------------------------------
  ! column types can have values of
  ! --------------------------------------------------------
  !   1  => (istsoil)          soil (vegetated or bare soil)
  !   2  => (istcrop)          crop (only for crop configuration)
  !   3  => (istice)           land ice
  !   4  => (istice_mec)       land ice (multiple elevation classes)
  !   5  => (istdlak)          deep lake
  !   6  => (istwet)           wetland
  !   71 => (icol_roof)        urban roof
  !   72 => (icol_sunwall)     urban sunwall
  !   73 => (icol_shadewall)   urban shadewall
  !   74 => (icol_road_imperv) urban impervious road
  !   75 => (icol_road_perv)   urban pervious road
  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varcon     , only : ispval
  implicit none
  save
  private

  !----------------------------------------------------
  ! column data type
  !----------------------------------------------------

  type, public :: column_type
    integer , pointer :: snl(:)            => null()    !number of snow layers
    real(r8), pointer :: zi(:,:)       => null()        !interface level below a "z" level (m) (-nlevsno+0:nlevgrnd)
    real(r8), pointer :: dz(:,:)        => null()       !layer thickness (m)  (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: z(:,:)       => null()         !layer depth (m) (-nlevsno+1:nlevgrnd)

    integer , pointer :: pfti                 (:)  => null() ! beginning pft index for each column
    integer , pointer :: pftf                 (:)  => null() ! ending pft index for each column
    real(r8), pointer :: wtgcell              (:)  => null() ! weight (relative to gridcell)
    integer , pointer :: gridcell             (:)  => null() ! index into gridcell level quantities

     logical , pointer :: active               (:) => null()  ! true=>do computations on this column
     integer , pointer :: landunit             (:) => null()  ! index into landunit level quantities
     integer , pointer :: itype                (:) => null()  ! column type
     real(r8), pointer :: wtlunit              (:) => null()  ! weight (relative to landunit)
     integer , pointer :: npfts                (:)  => null() ! number of patches for each column
  contains
    procedure          :: Init
    procedure, private :: InitAllocate
  end type column_type
  type(column_type), public, target :: col !column data structure (soil/snow/canopy columns)
  contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(column_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !ARGUMENTS:
    class(column_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    integer :: lbj,  ubj
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    lbj  = bounds%lbj; ubj = bounds%ubj

    allocate(this%snl(begc:endc))          ; this%snl(:) = ispval
    allocate(this%dz(begc:endc, lbj:ubj))  ;
    allocate(this%zi(begc:endc,lbj-1:ubj)) ;
    allocate(this%z(begc:endc,lbj:ubj))    ;
    allocate(this%landunit(begc:endc))     ; this%landunit(:) = ispval
    allocate(this%gridcell(begc:endc))     ; this%gridcell(:) = ispval
    allocate(this%npfts(begc:endc))
    allocate(this%pfti(begc:endc))
    allocate(this%pftf(begc:endc))
  end subroutine InitAllocate
end module ColumnType
