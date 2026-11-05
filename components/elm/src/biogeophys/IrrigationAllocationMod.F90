module IrrigationAllocationMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Pure helper routines for two-way (tw_irr) irrigation coupling with MOSART:
  ! (1) IrrigationWindowCloses: whether the river request window, of length
  !     river_dt [s], closes at the current driver time.
  ! (2) IrrigationDemandWindow: advance the patch request pipeline. The
  !     accumulated request [mm] becomes a pending mean rate [mm/s] when the
  !     window closes; the lagged request is the one the current supply answers.
  ! (3) IrrigationDemandFlux: split the gridcell surface supply, and the
  !     optional extra groundwater for the deficit, among patches in proportion
  !     to their lagged request. Each patch is capped at its own request and
  !     there is no division by patch area.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: IrrigationDemandFlux    !split the gridcell surface supply into patches
  public :: IrrigationDemandWindow  !advance the patch request pipeline
  public :: IrrigationWindowCloses  !test whether the river request window should close
  !-----------------------------------------------------------------------  
contains
  !-----------------------------------------------------------------------
  pure logical function IrrigationWindowCloses(driver_days, river_dt, river_phase)
    !
    ! !DESCRIPTION:
    ! Return true if the irrigation request window, which is aligned with the
    ! river coupling interval, closes at the current driver time.
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: driver_days   ! current driver clock time [days]
    integer, intent(in) :: river_dt       ! river coupling interval, i.e. request window length [s]
    integer, intent(in) :: river_phase    ! offset of the river run within river_dt [s]
    !-----------------------------------------------------------------------

    ! seq_timemgr_EClockGetData(curr_time=...) returns DAYS, not seconds.
    ! Reduce before converting to integer so century-scale runs do not overflow.
    IrrigationWindowCloses = modulo(nint(modulo(driver_days*86400._r8, &
         real(river_dt,r8)))-river_phase,river_dt) == 0
  end function IrrigationWindowCloses
  !-----------------------------------------------------------------------
  pure subroutine IrrigationDemandWindow(accum, pending, lagged, close_window, river_dt)
    !
    ! !DESCRIPTION:
    ! Advance the patch request pipeline by one land coupling step. The pending
    ! request becomes the lagged request that the current supply answers. When
    ! the window closes, the accumulated request becomes the new pending mean
    ! rate and the accumulator is reset.
    !
    ! !ARGUMENTS:
    real(r8), intent(inout) :: accum(:)    ! patch surface request accumulated in the open window [mm]; reset when it closes
    real(r8), intent(inout) :: pending(:)  ! patch request sent to the river, mean rate over the last closed window [mm/s]
    real(r8), intent(inout) :: lagged(:)   ! patch request that the current river supply answers [mm/s]
    logical, intent(in) :: close_window    ! true if the request window closes at this driver time
    integer, intent(in) :: river_dt        ! river coupling interval, i.e. request window length [s]
    !-----------------------------------------------------------------------

    lagged = pending
    if (close_window) then
      pending = accum/real(river_dt,r8)
      accum = 0._r8
    end if
  end subroutine IrrigationDemandWindow
  !-----------------------------------------------------------------------
  pure subroutine IrrigationDemandFlux(request, grid_request, supply, deficit, nonlake, extra_gw, surface, groundwater)
    !
    ! !DESCRIPTION:
    ! Allocate the gridcell river supply, and optionally extra groundwater for
    ! the unmet deficit, to one patch in proportion to its share of the gridcell
    ! request. The allocation never exceeds the patch request.
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: request       ! lagged surface request of this patch [mm/s per patch area]
    real(r8), intent(in)  :: grid_request  ! patch-weighted sum of request over the gridcell [mm/s per gridcell area]
    real(r8), intent(in)  :: supply        ! surface water supplied by the river [mm/s per gridcell area]
    real(r8), intent(in)  :: deficit       ! request unmet by the river supply [mm/s per gridcell area]
    real(r8), intent(in)  :: nonlake       ! summed gridcell weight of non-lake patches [-]
    logical,  intent(in)  :: extra_gw      ! true to meet the deficit with extra groundwater pumping
    real(r8), intent(out) :: surface       ! surface water irrigation applied to this patch [mm/s per patch area]
    real(r8), intent(out) :: groundwater   ! extra groundwater irrigation applied to this patch [mm/s per patch area]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: surface_volume             ! supply allocated to irrigation, capped at grid_request [mm/s per gridcell area]
    real(r8) :: groundwater_volume         ! extra groundwater allocated, capped at the unmet request [mm/s per gridcell area]
    !-----------------------------------------------------------------------

    surface = 0._r8
    groundwater = 0._r8
    if (request <= 0._r8 .or. grid_request <= 0._r8) return

    surface_volume = min(max(0._r8, supply)*nonlake, grid_request)
    
    ! Divide volumes before multiplying by patch demand. No inverse patch area.
    surface = (surface_volume/grid_request)*request
    if (extra_gw) then
      groundwater_volume = min(max(0._r8, deficit)*nonlake, max(0._r8, grid_request-surface_volume))
      groundwater = (groundwater_volume/grid_request)*request
    end if
  end subroutine IrrigationDemandFlux
end module IrrigationAllocationMod
