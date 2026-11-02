module IrrigationAllocationMod
  use shr_kind_mod, only : r8 => shr_kind_r8
  implicit none
  private
  public :: IrrigationDemandFlux, IrrigationDemandWindow, IrrigationWindowCloses
contains
  pure logical function IrrigationWindowCloses(driver_days, river_dt, river_phase)
    real(r8), intent(in) :: driver_days
    integer, intent(in) :: river_dt, river_phase
    ! seq_timemgr_EClockGetData(curr_time=...) returns DAYS, not seconds.
    ! Reduce before converting to integer so century-scale runs do not overflow.
    IrrigationWindowCloses = modulo(nint(modulo(driver_days*86400._r8, &
         real(river_dt,r8)))-river_phase,river_dt) == 0
  end function IrrigationWindowCloses

  pure subroutine IrrigationDemandWindow(accum, pending, lagged, close_window, river_dt)
    real(r8), intent(inout) :: accum(:), pending(:), lagged(:)
    logical, intent(in) :: close_window
    integer, intent(in) :: river_dt
    lagged = pending
    if (close_window) then
      pending = accum/real(river_dt,r8)
      accum = 0._r8
    end if
  end subroutine IrrigationDemandWindow

  pure subroutine IrrigationDemandFlux(request, grid_request, supply, deficit, nonlake, extra_gw, surface, groundwater)
    real(r8), intent(in) :: request, grid_request, supply, deficit, nonlake
    logical, intent(in) :: extra_gw
    real(r8), intent(out) :: surface, groundwater
    real(r8) :: surface_volume, groundwater_volume
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
