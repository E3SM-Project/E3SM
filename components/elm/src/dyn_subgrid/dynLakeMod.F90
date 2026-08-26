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
  use spmdMod         , only : masterproc, mpicom, MPI_INTEGER, MPI_SUM

  implicit none
  private

  public :: dynlake_driver

  integer , parameter :: max_print   = 10       ! per warning type, per call, per task
  ! Latch: the coupler hands ELM ZEROS (not spval) for rof fields until MOSART's first export
  ! (first rof coupling interval after a cold start or a restart). Updating from those zeros
  ! would collapse every lake to the lakebed on day 1, so no update is made on this task
  ! until a nonzero lake area/volume has been seen at least once.
  logical, save :: lake_fields_received = .false.
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
    integer  :: n_upd, n_over_cell, n_over_foot, n_nolun
    integer  :: nloc(4), nglob(4), ier          ! task / global warning counters
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

    ! Daily cadence, like the other dynamic-landunit sources (the rof->lnd exchange
    ! itself is every coupling interval; a landunit weight does not need sub-daily updates)
    if (.not. is_beg_curr_day()) return

    if (.not. lake_fields_received) then
       do g = bounds%begg, bounds%endg
          asur_r = atm2lnd_vars%lake_r_Asur_grc(g); asur_t = atm2lnd_vars%lake_t_Asur_grc(g)
          if (asur_r >= 0.5_r8*spval .or. asur_t >= 0.5_r8*spval) cycle
          if (asur_r > 0._r8 .or. asur_t > 0._r8 .or. &
              atm2lnd_vars%lake_r_Vtot_grc(g) > 0._r8 .or. atm2lnd_vars%lake_t_Vtot_grc(g) > 0._r8) then
             lake_fields_received = .true.
             exit
          end if
       end do
       if (.not. lake_fields_received) then
          write(iulog,'(a,i10,a)') subname//' nstep=', get_nstep(), &
               '  MOSART lake fields not received yet on this task -- lake/lakebed split unchanged'
          return
       end if
    end if

    n_upd = 0; n_over_cell = 0; n_over_foot = 0; n_nolun = 0

    do g = bounds%begg, bounds%endg

       asur_r = atm2lnd_vars%lake_r_Asur_grc(g)
       asur_t = atm2lnd_vars%lake_t_Asur_grc(g)

       ! Not received yet (spval until the first rof->lnd exchange): keep the current split
       if (asur_r >= 0.5_r8*spval .or. asur_t >= 0.5_r8*spval) cycle

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

    ! Global summary from the master task (per-task prints above only reach the master's log)
    nloc = (/ n_upd, n_over_cell, n_over_foot, n_nolun /)
    call mpi_allreduce(nloc, nglob, 4, MPI_INTEGER, MPI_SUM, mpicom, ier)
    if (masterproc .and. (nglob(2) > 0 .or. nglob(3) > 0 .or. nglob(4) > 0)) then
       write(iulog,'(a,i10,a,i8,a,i8,a,i8,a,i8)') subname//' nstep=', get_nstep(), &
            '  lake landunits updated=', nglob(1), &
            '  WARN lake>gridcell=', nglob(2), '  WARN lake>PCT_LAKE=', nglob(3), &
            '  WARN lake w/o landunit=', nglob(4)
    end if

  end subroutine dynlake_driver

end module dynLakeMod
