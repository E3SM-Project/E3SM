module WaterfluxType

  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun

  implicit none
  save
  private

  !----------------------------------------------------
  ! column water flux variables structure
  !----------------------------------------------------
  type, public :: waterflux_type
    real(r8), pointer :: qflx_adv_col(:,:)    => null()   !advection velocity from one layer to another, (0:nlevgrnd), positive downward
    real(r8), pointer :: qflx_infl_col(:)	    => null()   !infiltration (mm H2O /s)
    real(r8), pointer :: qflx_surf_col(:)	    => null()   !surface runoff (mm H2O /s)
    real(r8), pointer :: h2oliq_vol_tendency(:,:)   => null()      !temporal change of water during the solution of soil water movement
    real(r8), pointer :: qflx_gross_evap_soil_col (:) => null()  ! col gross infiltration from soil, this satisfies the relationship qflx_infl_col = qflx_gross_infl_soil_col-qflx_gross_evap_soil_col
    real(r8), pointer :: qflx_gross_infl_soil_col (:) => null()  ! col gross infiltration, before considering the evaporation, mm/s
    real(r8), pointer :: qflx_rootsoi_col         (:,:)=> null() ! col root and soil water exchange [mm H2O/s] [+ into root]
    real(r8), pointer :: qflx_drain_vr_col        (:,:)=> null() ! col liquid water losted as drainage (m /time step)
    real(r8), pointer :: qflx_totdrain_col        (:) => null()  ! col total liquid water drainage  (m/time step), updated in betr
    real(r8), pointer :: qflx_dew_grnd_col        (:)=> null()   ! col ground surface dew formation (mm H2O /s) [+] (+ = to atm); usually eflx_bot >= 0)
    real(r8), pointer :: qflx_dew_snow_col        (:) => null()  ! col surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_sub_snow_vol_col    (:)=> null()
    real(r8), pointer :: qflx_sub_snow_col        (:)=> null()   ! col sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_h2osfc2topsoi_col   (:) => null()  ! col liquid water coming from surface standing water top soil (mm H2O/s)
    real(r8), pointer :: qflx_snow2topsoi_col     (:) => null()  ! col liquid water coming from residual snow to topsoil (mm H2O/s)
    real(r8), pointer :: qflx_tran_veg_patch      (:)=> null()
    real(r8), pointer :: qflx_rootsoi_patch       (:,:) => null()! pft root and soil water exchange [mm H2O/s] [+ into atmosphere]
    real(r8), pointer :: qflx_rootsoi_frac_patch  (:,:) => null()

  contains
    procedure          :: Init
    procedure, private :: InitAllocate
  end type waterflux_type

  contains
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(waterflux_type) :: this
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
    class(waterflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    integer :: lbj,  ubj
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begp = bounds%begp; endp= bounds%endp
    lbj  = bounds%lbj; ubj = bounds%ubj

    allocate(this%qflx_adv_col(begc:endc, lbj-1:ubj))
    allocate(this%qflx_infl_col(begc:endc))
    allocate(this%qflx_surf_col(begc:endc))
    allocate(this%qflx_gross_evap_soil_col (begc:endc))              ; this%qflx_gross_evap_soil_col (:)   = nan
    allocate(this%qflx_gross_infl_soil_col (begc:endc))              ; this%qflx_gross_infl_soil_col (:)   = nan
    allocate(this%qflx_rootsoi_col         (begc:endc,lbj:ubj))      ; this%qflx_rootsoi_col         (:,:) = nan
    allocate(this%qflx_drain_vr_col        (begc:endc,lbj:ubj))      ; this%qflx_drain_vr_col        (:,:) = nan
    allocate(this%qflx_dew_grnd_col        (begc:endc))              ; this%qflx_dew_grnd_col        (:)   = nan
    allocate(this%qflx_dew_snow_col        (begc:endc))              ; this%qflx_dew_snow_col        (:)   = nan
    allocate(this%qflx_sub_snow_vol_col    (begc:endc))              ; this%qflx_sub_snow_vol_col    (:)   = 0._r8
    allocate(this%qflx_sub_snow_col        (begc:endc))              ; this%qflx_sub_snow_col        (:)   = 0.0_r8
    allocate(this%qflx_snow2topsoi_col     (begc:endc))              ; this%qflx_snow2topsoi_col     (:)   = nan
    allocate(this%qflx_h2osfc2topsoi_col   (begc:endc))              ; this%qflx_h2osfc2topsoi_col   (:)   = nan
    allocate(this%qflx_tran_veg_patch      (begp:endp))              ; this%qflx_tran_veg_patch      (:)   = nan
    allocate(this%qflx_rootsoi_patch       (begp:endp,lbj:ubj))      ; this%qflx_rootsoi_patch       (:,:) = nan
    allocate( this%qflx_totdrain_col       (begc:endc))              ; this%qflx_totdrain_col        (:)   = nan
    allocate(this%qflx_rootsoi_frac_patch  (begp:endp,lbj:ubj))      ; this%qflx_rootsoi_frac_patch(:,:) = nan
  end subroutine InitAllocate




end module WaterfluxType
