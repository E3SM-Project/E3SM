Module SoilHydrologyType

! DESCRIPTION
! derived data for soilhydrology

 use shr_kind_mod   , only : r8 => shr_kind_r8
 use decompMod      , only : bounds_type
implicit none

  type, public :: soilhydrology_type

  real(r8), pointer :: fracice_col       (:,:)  => null() ! col fractional impermeability (-)
  real(r8), pointer :: zwts_col           (:)  => null()  ! the shallower between zwt_perch and zwt
  real(r8), pointer :: qcharge_col      (:)  => null()   ! bottom of soil col flux, (mm/s)  
  contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
  end type soilhydrology_type

  contains


  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(soilhydrology_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    ! !ARGUMENTS:
    class(soilhydrology_type) :: this
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

    allocate(this%fracice_col       (begc:endc,lbj:ubj))        ; this%fracice_col       (:,:)   = nan
    allocate(this%zwts_col           (begc:endc))                ; this%zwts_col         (:)     = nan
    allocate(this%qcharge_col      (begc:endc))                 ; this%qcharge_col      (:)     = nan
  end subroutine InitAllocate
end Module SoilHydrologyType
