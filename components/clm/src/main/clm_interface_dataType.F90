module clm_interface_dataType

!=================================================================================================
! ALM Thermal(T)-Hydrology (H) & BioGeoChemistry (BGC) Interface: Data Type (Variables)
! created: 8/25/2015
! update: 9/16/2016, 2/2/2017, May-2017, June-2017
! update: 5/13/2019 (all are IMPLICITLY column-based coupling,
!        esp. after ELM v2 data-structure modification @3/12/2019, commit 1bf22e32d)
!=================================================================================================
  ! USES:
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_infnan_mod        , only : nan => shr_infnan_nan, assignment(=)

  use clm_interface_thType
  use clm_interface_bgcType

  implicit none

  private

!-------------------------------------------------------------------------------------------------
  type, public :: clm_interface_data_type

     ! coupling level (grid/column/patch):
     integer :: cpl_level = 2                                ! coupling level: 1 - grid, 2 - column (default), 3 - patch

     ! col dimension
     real(r8), pointer :: z                                     (:,:)   ! layer depth (m) (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: zi                                    (:,:)   ! layer depth (m) (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: dz                                    (:,:)   ! layer thickness (m)  (-nlevsno+1:nlevgrnd)
     integer,  pointer :: snl                                   (:)     ! number of snow layers (negative)

     ! soilstate_vars:
     real(r8), pointer :: bd                                    (:,:)   ! bulk density of dry soil material [kg/m^3] (CN)
     real(r8), pointer :: hksat                                 (:,:)   ! hydraulic conductivity at saturation (mm H2O /s)
     real(r8), pointer :: bsw                                   (:,:)   ! Clapp and Hornberger "b" (nlevgrnd)
     real(r8), pointer :: watsat                                (:,:)   ! volumetric soil water at saturation (porosity)
     real(r8), pointer :: watmin                                (:,:)   ! minimum volumetric soil water (nlevsoi)
     real(r8), pointer :: sucsat                                (:,:)   ! minimum soil suction (mm) (nlevgrnd)
     real(r8), pointer :: sucmin                                (:,:)   ! minimum allowable soil liquid suction pressure (mm) [Note: sucmin is a negative value, while sucsat is a positive quantity]
     real(r8), pointer :: watfc                                 (:,:)   ! volumetric soil water at field capacity (nlevsoi)
     real(r8), pointer :: porosity                              (:,:)   ! soil porisity (1-bulk_density/soil_density) (VIC)
     real(r8), pointer :: eff_porosity                          (:,:)   ! effective porosity = porosity - vol_ice (nlevgrnd)
     real(r8), pointer :: cellorg                               (:,:)   ! organic matter for gridcell containing column (1:nlevsoi)
     real(r8), pointer :: rootfr                                (:,:)   ! fraction of roots in each soil layer (nlevgrnd)

     real(r8), pointer :: tkfrz                                 (:,:)   ! thermal conductivity, frozen soil  [W/m-K] (new) (nlevgrnd)
     real(r8), pointer :: tkdry                                 (:,:)   ! thermal conductivity, dry soil (W/m/Kelvin) (nlevgrnd)
     real(r8), pointer :: tkwet                                 (:,:)   ! thermal conductivity, saturated soil [W/m-K] (new) (nlevgrnd)
     real(r8), pointer :: csol                                  (:,:)   ! heat capacity, soil solids (J/m**3/Kelvin) (nlevgrnd)

     ! thermal-hydrology:
     type(clm_interface_th_datatype) :: th

     ! biogeochemistry:
     type(clm_interface_bgc_datatype):: bgc

     !
  contains
     procedure , public  :: Init
     procedure , private :: InitAllocate
  end type clm_interface_data_type
!-------------------------------------------------------------------------------------------------
  

contains

!-------------------------------------------------------------------------------------------------
  subroutine Init(this, bounds)
     use decompMod               , only : bounds_type
     class(clm_interface_data_type) :: this
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
    class(clm_interface_data_type)      :: this
    type(bounds_type), intent(in)       :: bounds

    ! LOCAL VARIABLES:
    integer  :: begc, endc
    !------------------------------------------------------------------------
    begc = bounds%begc; endc= bounds%endc
    if (this%cpl_level==1) then
      begc = bounds%begg; endc= bounds%endg
    elseif (this%cpl_level==3) then
      begc = bounds%begp; endc= bounds%endp
    endif

    ! col:
    allocate(this%z                     (begc:endc,-nlevsno+1:nlevgrnd))      ; this%z                    (:,:) = nan
    allocate(this%zi                    (begc:endc,-nlevsno+0:nlevgrnd))      ; this%zi                   (:,:) = nan
    allocate(this%dz                    (begc:endc,-nlevsno+1:nlevgrnd))      ; this%dz                   (:,:) = nan
    allocate(this%snl                   (begc:endc))                          ; this%snl                  (:)   = 0

    ! soilstate_vars:
    allocate(this%bd                (begc:endc,1:nlevgrnd))               ; this%bd               (:,:) = nan
    allocate(this%hksat             (begc:endc,1:nlevgrnd))               ; this%hksat            (:,:) = spval
    allocate(this%bsw               (begc:endc,1:nlevgrnd))               ; this%bsw              (:,:) = nan
    allocate(this%watsat            (begc:endc,1:nlevgrnd))               ; this%watsat           (:,:) = nan
    allocate(this%watmin            (begc:endc,1:nlevgrnd))               ; this%watmin           (:,:) = nan
    allocate(this%sucsat            (begc:endc,1:nlevgrnd))               ; this%sucsat           (:,:) = spval
    allocate(this%sucmin            (begc:endc,1:nlevgrnd))               ; this%sucmin           (:,:) = spval
    allocate(this%watfc             (begc:endc,1:nlevgrnd))               ; this%watfc            (:,:) = nan
    allocate(this%porosity          (begc:endc,1:nlevgrnd))               ; this%porosity         (:,:) = spval
    allocate(this%eff_porosity      (begc:endc,1:nlevgrnd))               ; this%eff_porosity     (:,:) = spval
    allocate(this%cellorg           (begc:endc,1:nlevgrnd))               ; this%cellorg          (:,:) = nan
    allocate(this%rootfr            (begc:endc,1:nlevgrnd))               ; this%rootfr           (:,:) = nan

    allocate(this%tkwet             (begc:endc,1:nlevgrnd))               ; this%tkwet            (:,:) = nan
    allocate(this%tkdry             (begc:endc,1:nlevgrnd))               ; this%tkdry            (:,:) = nan
    allocate(this%tkfrz             (begc:endc,1:nlevgrnd))               ; this%tkfrz            (:,:) = nan
    allocate(this%csol              (begc:endc,1:nlevgrnd))               ; this%csol             (:,:) = nan

  end subroutine InitAllocate
!-------------------------------------------------------------------------------------------------

end module clm_interface_dataType
