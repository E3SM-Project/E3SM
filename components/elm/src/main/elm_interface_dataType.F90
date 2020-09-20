module elm_interface_dataType

!=================================================================================================
! ALM Thermal(T)-Hydrology (H) & BioGeoChemistry (BGC) Interface: Data Type (Variables)
! created: 8/25/2015
! update: 9/16/2016, 2/2/2017, May-2017, June-2017
!=================================================================================================
  ! USES:
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_infnan_mod        , only : nan => shr_infnan_nan, assignment(=)

  use elm_interface_thType
  use elm_interface_bgcType

  implicit none

  private

!-------------------------------------------------------------------------------------------------
  type, public :: elm_interface_data_type

     ! col dimension
     real(r8), pointer :: z                                         (:,:)   ! layer depth (m) (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: zi                                        (:,:)   ! layer depth (m) (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: dz                                        (:,:)   ! layer thickness (m)  (-nlevsno+1:nlevgrnd)

     ! soilstate_vars:
     real(r8), pointer :: bd_col                                    (:,:)   ! col bulk density of dry soil material [kg/m^3] (CN)
     real(r8), pointer :: hksat_col                                 (:,:)   ! col hydraulic conductivity at saturation (mm H2O /s)
     real(r8), pointer :: bsw_col                                   (:,:)   ! col Clapp and Hornberger "b" (nlevgrnd)
     real(r8), pointer :: watsat_col                                (:,:)   ! col volumetric soil water at saturation (porosity)
     real(r8), pointer :: watmin_col                                (:,:)   ! col minimum volumetric soil water (nlevsoi)
     real(r8), pointer :: sucsat_col                                (:,:)   ! col minimum soil suction (mm) (nlevgrnd)
     real(r8), pointer :: sucmin_col                                (:,:)   ! col minimum allowable soil liquid suction pressure (mm) [Note: sucmin_col is a negative value, while sucsat_col is a positive quantity]
     real(r8), pointer :: watfc_col                                 (:,:)   ! col volumetric soil water at field capacity (nlevsoi)
     real(r8), pointer :: porosity_col                              (:,:)   ! col soil porisity (1-bulk_density/soil_density) (VIC)
     real(r8), pointer :: eff_porosity_col                          (:,:)   ! col effective porosity = porosity - vol_ice (nlevgrnd)
     real(r8), pointer :: cellorg_col                               (:,:)   ! col organic matter for gridcell containing column (1:nlevsoi)
     real(r8), pointer :: rootfr_col                                (:,:)   ! col fraction of roots in each soil layer (nlevgrnd)

     real(r8), pointer :: tkfrz_col                                 (:,:)   ! col thermal conductivity, frozen soil  [W/m-K] (new) (nlevgrnd)
     real(r8), pointer :: tkdry_col                                 (:,:)   ! col thermal conductivity, dry soil (W/m/Kelvin) (nlevgrnd)
     real(r8), pointer :: tkwet_col                                 (:,:)   ! col thermal conductivity, saturated soil [W/m-K] (new) (nlevgrnd)
     real(r8), pointer :: csol_col                                  (:,:)   ! col heat capacity, soil solids (J/m**3/Kelvin) (nlevgrnd)

     ! thermal-hydrology:
     type(elm_interface_th_datatype) :: th

     ! biogeochemistry:
     type(elm_interface_bgc_datatype):: bgc

     !
  contains
     procedure , public  :: Init
     procedure , private :: InitAllocate
  end type elm_interface_data_type
!-------------------------------------------------------------------------------------------------
  

contains

!-------------------------------------------------------------------------------------------------
  subroutine Init(this, bounds)
     use decompMod               , only : bounds_type
     class(elm_interface_data_type) :: this
     type(bounds_type), intent(in)  :: bounds

     call this%InitAllocate (bounds)

     call this%th%Init(bounds)

     call this%bgc%Init(bounds)

  end subroutine Init
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------

  subroutine InitAllocate(this, bounds)
    ! USES
    use clm_varpar            , only : nlevsno, nlevgrnd
    use clm_varcon            , only : spval
    use decompMod             , only : bounds_type

    ! ARGUMENTS:
    real(r8) :: ival  = 0.0_r8  ! initial value
    class(elm_interface_data_type)      :: this
    type(bounds_type), intent(in)       :: bounds

    ! LOCAL VARIABLES:
    integer  :: begg, endg
    integer  :: begc, endc
    integer  :: begp, endp
    !------------------------------------------------------------------------
    begg = bounds%begg; endg= bounds%endg
    begc = bounds%begc; endc= bounds%endc
    begp = bounds%begp; endp= bounds%endp

    ! col:
    allocate(this%z                     (begc:endc,-nlevsno+1:nlevgrnd))      ; this%z                    (:,:) = nan
    allocate(this%zi                    (begc:endc,-nlevsno+0:nlevgrnd))      ; this%zi                   (:,:) = nan
    allocate(this%dz                    (begc:endc,-nlevsno+1:nlevgrnd))      ; this%dz                   (:,:) = nan

    ! soilstate_vars:
    allocate(this%bd_col                (begc:endc,1:nlevgrnd))               ; this%bd_col               (:,:) = nan
    allocate(this%hksat_col             (begc:endc,1:nlevgrnd))               ; this%hksat_col            (:,:) = spval
    allocate(this%bsw_col               (begc:endc,1:nlevgrnd))               ; this%bsw_col              (:,:) = nan
    allocate(this%watsat_col            (begc:endc,1:nlevgrnd))               ; this%watsat_col           (:,:) = nan
    allocate(this%watmin_col            (begc:endc,1:nlevgrnd))               ; this%watmin_col           (:,:) = nan
    allocate(this%sucsat_col            (begc:endc,1:nlevgrnd))               ; this%sucsat_col           (:,:) = spval
    allocate(this%sucmin_col            (begc:endc,1:nlevgrnd))               ; this%sucmin_col           (:,:) = spval
    allocate(this%watfc_col             (begc:endc,1:nlevgrnd))               ; this%watfc_col            (:,:) = nan
    allocate(this%porosity_col          (begc:endc,1:nlevgrnd))               ; this%porosity_col         (:,:) = spval
    allocate(this%eff_porosity_col      (begc:endc,1:nlevgrnd))               ; this%eff_porosity_col     (:,:) = spval
    allocate(this%cellorg_col           (begc:endc,1:nlevgrnd))               ; this%cellorg_col          (:,:) = nan
    allocate(this%rootfr_col            (begc:endc,1:nlevgrnd))               ; this%rootfr_col           (:,:) = nan

    allocate(this%tkwet_col             (begc:endc,1:nlevgrnd))               ; this%tkwet_col            (:,:) = nan
    allocate(this%tkdry_col             (begc:endc,1:nlevgrnd))               ; this%tkdry_col            (:,:) = nan
    allocate(this%tkfrz_col             (begc:endc,1:nlevgrnd))               ; this%tkfrz_col            (:,:) = nan
    allocate(this%csol_col              (begc:endc,1:nlevgrnd))               ; this%csol_col             (:,:) = nan

  end subroutine InitAllocate
!-------------------------------------------------------------------------------------------------

end module elm_interface_dataType
