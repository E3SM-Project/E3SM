module SoilStateType
  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varcon      , only : spval
  implicit none
  save
  private

  !----------------------------------------------------
  ! column physical state variables structure
  !----------------------------------------------------
  type, public :: soilstate_type
    real(r8), pointer :: bsw_col(:,:)        => null()      !Clapp and Hornberger "b" (nlevgrnd)
    real(r8), pointer :: watsat_col(:,:)      => null()     !volumetric soil water at saturation (porosity) (nlevgrnd)
    real(r8), pointer :: eff_porosity_col(:,:) => null()    !effective porosity = porosity - vol_ice (nlevgrnd)
    real(r8), pointer :: soilpsi_col          (:,:) => null() ! col soil water potential in each soil layer (MPa) (CN)
    real(r8), pointer :: cellorg_col          (:,:)=> null() ! col organic matter for gridcell containing column (1:nlevsoi)
    real(r8), pointer :: cellclay_col         (:,:) => null()! clay value for gridcell containing column (1:nlevsoi)
    real(r8), pointer :: cellsand_col         (:,:) => null()! sand value for gridcell containing column (1:nlevsoi)
    real(r8), pointer :: bd_col               (:,:)=> null() ! col bulk density of dry soil material [kg/m^3] (CN)
    real(r8), pointer :: watfc_col            (:,:)=> null() ! col volumetric soil water at field capacity (nlevsoi)
    real(r8), pointer :: sucsat_col           (:,:)=> null() ! col minimum soil suction (mm) (nlevgrnd)
    real(r8), pointer :: rootfr_patch         (:,:) => null()! patch fraction of roots in each soil layer (nlevgrnd)
  contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
  end type soilstate_type

  contains


  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(soilstate_type) :: this
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
    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    integer :: lbj,  ubj
    !------------------------------------------------------------------------

    begp = bounds%begp; endp=bounds%endp
    begc = bounds%begc; endc= bounds%endc
    lbj  = bounds%lbj; ubj = bounds%ubj
    allocate(this%bsw_col(begc:endc, lbj:ubj))         ; this%bsw_col(:,:) = nan
    allocate(this%watsat_col(begc:endc, lbj:ubj))      ; this%watsat_col(:,:) = nan
    allocate(this%eff_porosity_col(begc:endc, lbj:ubj)); this%eff_porosity_col(:,:) = nan
    allocate(this%soilpsi_col          (begc:endc,lbj:ubj))            ; this%soilpsi_col          (:,:) = nan
    allocate(this%cellorg_col          (begc:endc,lbj:ubj))            ; this%cellorg_col          (:,:) = nan
    allocate(this%cellclay_col         (begc:endc,lbj:ubj))            ; this%cellclay_col         (:,:) = nan
    allocate(this%cellsand_col         (begc:endc,lbj:ubj))            ; this%cellsand_col         (:,:) = nan
    allocate(this%bd_col               (begc:endc,lbj:ubj))            ; this%bd_col               (:,:) = nan
    allocate(this%watfc_col            (begc:endc,lbj:ubj))            ; this%watfc_col            (:,:) = nan
    allocate(this%sucsat_col           (begc:endc,lbj:ubj))            ; this%sucsat_col           (:,:) = spval
    allocate(this%rootfr_patch         (begp:endp,lbj:ubj))            ; this%rootfr_patch         (:,:) = nan
  end subroutine InitAllocate
end module SoilStateType
