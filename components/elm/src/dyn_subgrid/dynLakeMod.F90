module dynLakeMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Dynamic lake fraction driven by MOSART-Lake (WP-A of the E3SM dynamic-lake
  ! coupling; heat-off, water-only).
  !
  ! Ownership split: MOSART-Lake owns the lake WATER (storage, area, inflow/outflow);
  ! ELM owns the surface ENERGY balance over the lake fraction. Each day this
  ! module reads the lake surface area MOSART sent through the coupler
  ! (atm2lnd_vars%lake_{r,t}_Asur_grc, m2, same grid) and re-splits the lake
  ! landunit between its two columns:
  !
  !   lake column    (col_pp%is_lake)  weight r      = wet part of the footprint
  !   lakebed column (col_pp%is_soil)  weight 1 - r  = exposed lakebed (soil physics)
  !
  ! The lake LANDUNIT weight (PCT_LAKE on the surface dataset) never changes: it
  ! is the maximum footprint the lake can occupy, and
  !
  !   r = (Asur / A_land) / wt_landunit,   A_land = gridcell area * land fraction
  !
  ! Two inconsistencies between the MOSART lake data and the ELM surface data
  ! are WARNED about (counted, first few printed), never aborted on:
  !   (1) Asur / A_land > 1          the MOSART lake is larger than the gridcell
  !   (2) r > 1                      the MOSART lake exceeds the PCT_LAKE footprint
  ! r is clamped to [0,1] in both cases. With consistently generated datasets
  ! (one lake mask -> PCT_LAKE and MOSART lake parameters) neither fires.
  !
  ! Requires dyn_lake=.true. (driver), create_lakebed_column=.true. (ELM),
  ! lnd and rof on the same grid (checked in prep_lnd_mod). Water/heat of the
  ! weight change is handled by the dynamic-subgrid conservation machinery
  ! (dyn_hwcontent_init/final -> qflx_liq_dynbal, eflx_dynbal); the lake water
  ! body itself is excluded there because MOSART owns it.
  !
  ! History: FLAKE_DYN (r * wt_landunit, gridcell wet-lake fraction) and
  ! LAKE_ASUR_RATIO (r before clamping) on atm2lnd_vars.
  !---------------------------------------------------------------------------

#include "shr_assert.h"

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use abortutils      , only : endrun
  use elm_varctl      , only : iulog, use_dyn_lake, create_lakebed_column, create_twolakes_per_gridcell
  use elm_varcon      , only : spval, ispval
  use landunit_varcon , only : istdlak
  use decompMod       , only : bounds_type
  use GridcellType    , only : grc_pp
  use TopounitType    , only : top_pp
  use LandunitType    , only : lun_pp
  use ColumnType      , only : col_pp
  use domainMod       , only : ldomain
  use atm2lndType     , only : atm2lnd_type
  use elm_time_manager, only : is_beg_curr_day, get_nstep
  use spmdMod         , only : masterproc

  implicit none
  private

  public :: dynlake_driver

  integer , parameter :: max_print   = 10       ! per warning type, per call, per task
  ! Per-cell: has this cell's split been set from a valid MOSART area since (re)start? Until it
  ! has, the split is applied at the first step a valid area arrives (not only at the day
  ! boundary), so the fully-wet initial split lasts one coupling interval, not a day+.
  logical, allocatable, save :: split_initialized(:)
  ! The coupler hands ELM ZEROS (not spval) for rof fields until MOSART's first export
  ! (first rof coupling interval after a cold start or a restart). Updating from those zeros
  ! would collapse every lake to the lakebed on day 1. MOSART therefore exports
  ! Sr_lake_valid = 1 with its lake fields; a cell is updated only once that flag has arrived.
  ! Per cell and data-carried, so no MPI collective is needed (PE-layout independent).
  real(r8), parameter :: ratio_tol   = 1.e-6_r8 ! tolerance before a ">1" counts as a warning

contains

  !-----------------------------------------------------------------------
  subroutine dynlake_driver (bounds, atm2lnd_vars)
    !
    ! !DESCRIPTION:
    ! Once per day (first time step of the day), set the lake / lakebed column
    ! weights of every lake landunit from the MOSART-Lake surface area.
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds
    type(atm2lnd_type), intent(inout) :: atm2lnd_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: g, t, l, c
    real(r8) :: asur_r, asur_t, asur        ! MOSART lake surface area (m2)
    real(r8) :: a_land                      ! gridcell land area (m2)
    real(r8) :: f_lake                      ! Asur / A_land
    real(r8) :: wt_l                        ! lake landunit weight in topounit (max footprint)
    real(r8) :: r                           ! wet share of the footprint
    real(r8) :: flake_g                     ! gridcell wet-lake fraction (diagnostic)
    integer  :: n_upd, n_over_cell, n_over_foot, n_nolun, n_nomos, n_wait
    character(len=*), parameter :: subname = 'dynlake_driver'
    !-----------------------------------------------------------------------

    if (.not. use_dyn_lake) return

    if (.not. create_lakebed_column) then
       call endrun(msg=subname//': use_dyn_lake requires create_lakebed_column=.true.'//errMsg(__FILE__, __LINE__))
    end if
    if (create_twolakes_per_gridcell) then
       call endrun(msg=subname//': create_twolakes_per_gridcell is not yet driven by MOSART-Lake; set it .false.'// &
            errMsg(__FILE__, __LINE__))
    end if

    if (.not. allocated(split_initialized)) then
       allocate(split_initialized(bounds%begg:bounds%endg)); split_initialized(:) = .false.
    end if

    ! Daily cadence, like the other dynamic-landunit sources (the rof->lnd exchange
    ! itself is every coupling interval; a landunit weight does not need sub-daily updates) --
    ! except for cells whose split has not yet been set from a valid area (cold start / restart)
    if (.not. is_beg_curr_day() .and. all(split_initialized)) return

    n_upd = 0; n_over_cell = 0; n_over_foot = 0; n_nolun = 0; n_nomos = 0; n_wait = 0

    do g = bounds%begg, bounds%endg

       asur_r = atm2lnd_vars%lake_r_Asur_grc(g)
       asur_t = atm2lnd_vars%lake_t_Asur_grc(g)

       ! Not received yet (spval, or the coupler's zero buffer before MOSART's first export):
       ! keep the current split (also the first day after a restart)
       if (asur_r >= 0.5_r8*spval .or. asur_t >= 0.5_r8*spval) cycle
       if (atm2lnd_vars%lake_valid_grc(g) < 0.5_r8) then
          n_wait = n_wait + 1
          cycle
       end if
       if (.not. is_beg_curr_day() .and. split_initialized(g)) cycle   ! mid-day: only first-time splits
       split_initialized(g) = .true.

       asur   = max(asur_r, 0._r8) + max(asur_t, 0._r8)
       a_land = ldomain%area(g) * 1.e6_r8 * ldomain%frac(g)     ! km2 -> m2, land part only
       if (a_land <= 0._r8) cycle

       f_lake = asur / a_land

       ! Warning (1): MOSART lake larger than the gridcell land area
       if (f_lake > 1._r8 + ratio_tol) then
          n_over_cell = n_over_cell + 1
          if (n_over_cell <= max_print) then
             write(iulog,'(a,2f9.3,a,es11.3,a,es11.3,a,f8.3)') subname// &
                  ' WARNING lake area > gridcell land area at (lat,lon)=', &
                  grc_pp%latdeg(g), grc_pp%londeg(g), '  Asur[m2]=', asur, '  A_land[m2]=', a_land, &
                  '  ratio=', f_lake
          end if
       end if

       flake_g = 0._r8

       do t = grc_pp%topi(g), grc_pp%topf(g)

          l = top_pp%landunit_indices(istdlak, t)
          if (l == ispval) then
             ! MOSART has a lake here but ELM has no lake landunit (PCT_LAKE = 0): nothing to move
             if (asur > 0._r8) then
                n_nolun = n_nolun + 1
                if (n_nolun <= max_print) then
                   write(iulog,'(a,2f9.3,a,es11.3)') subname// &
                        ' WARNING MOSART lake but no ELM lake landunit (PCT_LAKE=0) at (lat,lon)=', &
                        grc_pp%latdeg(g), grc_pp%londeg(g), '  Asur[m2]=', asur
                end if
             end if
             atm2lnd_vars%lake_asur_ratio_grc(g) = 0._r8
             cycle
          end if

          wt_l = lun_pp%wttopounit(l)
          if (wt_l <= 0._r8) then
             atm2lnd_vars%lake_asur_ratio_grc(g) = 0._r8
             cycle
          end if

          ! Warning (3): ELM has a lake landunit here but MOSART reports no lake at all
          ! (area and volume both zero after the first exchange): the footprint becomes lakebed
          if (asur <= 0._r8 .and. atm2lnd_vars%lake_r_Vtot_grc(g) <= 0._r8 .and. &
              atm2lnd_vars%lake_t_Vtot_grc(g) <= 0._r8) then
             n_nomos = n_nomos + 1
             if (n_nomos <= max_print) then
                write(iulog,'(a,2f9.3,a,f8.4)') subname// &
                     ' WARNING ELM lake landunit but no MOSART lake at (lat,lon)=', &
                     grc_pp%latdeg(g), grc_pp%londeg(g), '  PCT_LAKE/100=', wt_l
             end if
          end if

          r = f_lake / wt_l
          atm2lnd_vars%lake_asur_ratio_grc(g) = r        ! unclamped, so the inconsistency is visible in history

          ! Warning (2): MOSART lake exceeds the PCT_LAKE (maximum) footprint
          if (r > 1._r8 + ratio_tol) then
             n_over_foot = n_over_foot + 1
             if (n_over_foot <= max_print) then
                write(iulog,'(a,2f9.3,a,f8.4,a,f8.4,a,f8.3)') subname// &
                     ' WARNING lake area exceeds PCT_LAKE footprint at (lat,lon)=', &
                     grc_pp%latdeg(g), grc_pp%londeg(g), '  Asur/A_land=', f_lake, '  PCT_LAKE/100=', wt_l, &
                     '  ratio=', r
             end if
          end if
          r = max(0._r8, min(1._r8, r))

          ! Re-split the lake landunit: lake column = wet, lakebed soil column = exposed
          do c = lun_pp%coli(l), lun_pp%colf(l)
             if (col_pp%is_lake(c)) then
                col_pp%wtlunit(c) = r
             else if (col_pp%is_soil(c)) then
                col_pp%wtlunit(c) = 1._r8 - r
             end if
          end do

          flake_g = flake_g + r * wt_l * top_pp%wtgcell(t)
          n_upd = n_upd + 1
       end do

       atm2lnd_vars%flake_dyn_grc(g) = flake_g
    end do

    ! Daily summary for THIS task (no collective on the physics path; every task writes to
    ! its own iulog, the master's reaches lnd.log). Counts > 0 flag WP-C dataset inconsistency.
    if (n_over_cell > 0 .or. n_over_foot > 0 .or. n_nolun > 0 .or. n_nomos > 0 .or. &
        (masterproc .and. n_wait > 0)) then
       write(iulog,'(a,i10,a,i8,a,i8,a,i8,a,i8,a,i8,a,i8)') subname//' nstep=', get_nstep(), &
            '  [this task] lake landunits updated=', n_upd, '  waiting for MOSART=', n_wait, &
            '  WARN lake>gridcell=', n_over_cell, '  WARN lake>PCT_LAKE=', n_over_foot, &
            '  WARN MOSART lake w/o ELM landunit=', n_nolun, '  WARN ELM lake w/o MOSART lake=', n_nomos
    end if

  end subroutine dynlake_driver

end module dynLakeMod
