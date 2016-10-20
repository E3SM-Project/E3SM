module ColumnEnergyFluxType

  #include "shr_assert.h"

  !------------------------------------------------------------------------------
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use clm_varcon     , only : spval
  use decompMod      , only : bounds_type
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use PatchType      , only : pft                
  !-------------------------------------------------------------------------------
  implicit none
  save
  private





  type, public :: soilcol_energy_flux
   real(r8), pointer :: eflx_snomelt_col        (:)   ! col snow melt heat flux (W/m**2)
   real(r8), pointer :: eflx_snomelt_r_col      (:)   ! col rural snow melt heat flux (W/m**2)
   real(r8), pointer :: eflx_snomelt_u_col      (:)   ! col urban snow melt heat flux (W/m**2)
   real(r8), pointer :: eflx_bot_col            (:)   ! col heat flux from beneath the soil or ice column (W/m**2)
   real(r8), pointer :: eflx_fgr12_col          (:)   ! col ground heat flux between soil layers 1 and 2 (W/m**2)
   real(r8), pointer :: eflx_fgr_col            (:,:) ! col (rural) soil downward heat flux (W/m2) (1:nlevgrnd)  (pos upward; usually eflx_bot >= 0)
   real(r8), pointer :: eflx_building_heat_col  (:)   ! col heat flux from urban building interior to urban walls, roof (W/m**2)
   real(r8), pointer :: eflx_urban_ac_col       (:)   ! col urban air conditioning flux (W/m**2)
   real(r8), pointer :: eflx_urban_heat_col     (:)   ! col urban heating flux (W/m**2)
   real(r8), pointer :: errsoi_col              (:)   ! soil/lake energy conservation error   (W/m**2)

   ! Latent heat
   real(r8), pointer :: htvp_col                (:)   ! latent heat of vapor of water (or sublimation) [j/kg]

   ! Balance Checks
   real(r8), pointer :: errseb_col              (:)   ! surface energy conservation error     (W/m**2)
   real(r8), pointer :: errsol_col              (:)   ! solar radiation conservation error    (W/m**2)
   real(r8), pointer :: errlon_col              (:)   ! longwave radiation conservation error (W/m**2)

  contains
    procedure, public :: Init => init_col_ef
    procedure, public :: InitAcclocate => initallocate_col_ef
    procedure, public :: Clean => clean_col_ef

end type soilcol_energy_flux


! Begin Subroutines
! --------------------------------------------------------------------

subroutine init_col_ef(this, bounds, t_grnd_col)

    class(energyflux_type)         :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: t_grnd_col( bounds%begc: )

    SHR_ASSERT_ALL((ubound(t_grnd_col) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))

    call this%InitAllocate ( bounds )
    call this%InitHistory ( bounds )
    call this%InitCold ( bounds, t_grnd_col ) 

end subroutine init_col_ef




subroutine initallocate_col_ef(this, begc, endc)
  class(soilcol_energy_flux) :: this
  integer, intent(in) :: begc   ! beginning soil column index
  integer, intent(in) :: endc   ! ending soil column index    

  allocate( this%eflx_bot_col            (begc:endc))             ; this%eflx_bot_col            (:)   = nan
  allocate( this%eflx_snomelt_col        (begc:endc))             ; this%eflx_snomelt_col        (:)   = nan
  allocate( this%eflx_snomelt_r_col      (begc:endc))             ; this%eflx_snomelt_r_col      (:)   = nan
  allocate( this%eflx_snomelt_u_col      (begc:endc))             ; this%eflx_snomelt_u_col      (:)   = nan
  allocate( this%eflx_fgr12_col          (begc:endc))             ; this%eflx_fgr12_col          (:)   = nan
  allocate( this%eflx_fgr_col            (begc:endc, 1:nlevgrnd)) ; this%eflx_fgr_col            (:,:) = nan
  allocate( this%eflx_building_heat_col  (begc:endc))             ; this%eflx_building_heat_col  (:)   = nan
  allocate( this%eflx_urban_ac_col       (begc:endc))             ; this%eflx_urban_ac_col       (:)   = nan
  allocate( this%eflx_urban_heat_col     (begc:endc))             ; this%eflx_urban_heat_col     (:)   = nan
  allocate( this%htvp_col                (begc:endc))             ; this%htvp_col                (:)   = nan
  allocate( this%errsoi_col              (begc:endc))             ; this%errsoi_col              (:)   = nan
  allocate( this%errseb_col              (begc:endc))             ; this%errseb_col              (:)   = nan
  allocate( this%errsol_col              (begc:endc))             ; this%errsol_col              (:)   = nan
  allocate( this%errlon_col              (begc:endc))             ; this%errlon_col              (:)   = nan

end subroutine init_col_ef

subroutine clean_col_ef(this)
    class(soilcol_energy_flux) :: this

    deallocate( this%eflx_bot_col)
    deallocate( this%eflx_snomelt_col )
    deallocate( this%eflx_snomelt_r_col )
    deallocate( this%eflx_snomelt_u_col )
    deallocate( this%eflx_fgr12_col )
    deallocate( this%eflx_fgr_col )
    deallocate( this%eflx_building_heat_col )
    deallocate( this%eflx_urban_ac_col )
    deallocate( this%eflx_urban_heat_col )
    deallocate( this%htvp_col )
    deallocate( this%errsoi_col )
    deallocate( this%errseb_col )
    deallocate( this%errsol_col )
    deallocate( this%errlon_col )

end subroutine clean_col_ef

! --------------------------------------------------------------------
! End Subroutines

end module ColumnEnergyFluxType 
