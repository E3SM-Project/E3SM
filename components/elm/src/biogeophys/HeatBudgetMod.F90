module HeatBudgetMod
  ! !DESCRIPTION:
  ! Heat (energy) budget diagnostic for ELM's lnd.log, mirroring
  ! WaterBudgetMod. Two tables are printed:
  !
  ! - NET HEAT FLUXES: a cross-boundary flux-conservation check
  !   (comparable to the coupler's own NET HEAT BUDGET table in cpl.log).
  ! - HEAT STATES: a beginning/end-of-timestep snapshot of column heat
  !   content (col_es%hc_soisno, snow+soil+lake) plus the soil/lake
  !   solver's own numerical-closure residual (col_ef%errsoi), both
  !   already validated every timestep by SoilFluxesMod/LakeTemperatureMod/
  !   BalanceCheckMod -- see HeatBudget_Run's header for what this
  !   table can and cannot be used to conclude, and HEAT_BUDGET_PLAN.md's
  !   Phase 2 section for the full reasoning.
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use shr_sys_mod       , only : shr_sys_abort
  use decompMod         , only : bounds_type
  use abortutils        , only : endrun
  use elm_varctl        , only : iulog
  use atm2lndType       , only : atm2lnd_type
  use lnd2atmType       , only : lnd2atm_type
  use spmdMod           , only : masterproc
  use elm_varcon        , only : spval
  use ColumnDataType    , only : col_es, col_ef, col_ws
  use VegetationDataType, only : veg_ef
  use GridcellDataType  , only : grc_es, grc_ef
  use timeinfoMod

  implicit none
  save
  private

  public :: HeatBudget_Reset
  public :: HeatBudget_Run
  public :: HeatBudget_Accum
  public :: HeatBudget_Print
  public :: HeatBudget_Restart
  public :: HeatBudget_SetBeginningStates
  public :: HeatBudget_SetEndingStates

  !--- F for flux ---
  ! Named/ordered to match the coupler's own NET HEAT BUDGET rows
  ! (hnetsw/hlwdn/hlwup/hlatvap/hsen) for direct cross-validation against
  ! cpl.log's 'lnd' column.

  integer, parameter :: f_netsw = 1
  integer, parameter :: f_lwdn  = 2
  integer, parameter :: f_lwup  = 3
  integer, parameter :: f_lh    = 4
  integer, parameter :: f_sh    = 5

  integer, parameter, public :: f_size = f_sh

  character(len=12),parameter :: fname(f_size) = &
       (/&
       '      netsw ', &
       '       lwdn ', &
       '       lwup ', &
       '        lh  ', &
       '        sh  '  &
       /)

  !--- P for period ---

  integer, parameter :: p_inst = 1
  integer, parameter :: p_day  = 2
  integer, parameter :: p_mon  = 3
  integer, parameter :: p_ann  = 4
  integer, parameter :: p_inf  = 5

  integer, parameter, public :: p_size = p_inf

  character(len=8),parameter :: pname(p_size) = &
       (/'    inst','   daily',' monthly','  annual','all_time' /)

  real(r8) :: budg_fluxL(f_size, p_size) ! local sum, valid on all pes
  real(r8) :: budg_fluxG(f_size, p_size) ! global sum, valid only on root pe
  real(r8) :: budg_fluxN(f_size, p_size) ! counter, valid only on root pe

  !--- G for the into-ground flux table ---
  ! Rows that (up to errsoi_col's own tiny residual) sum to exactly the
  ! rate of change of col_es%hc_soisno -- i.e. this table's *SUM*, unlike
  ! NET HEAT FLUXES' *SUM*, is directly comparable to HEAT STATES'
  ! *NET CHANGE*. Derived from the exact identity
  !   d(cv*T)/dt = cv*(dT/dt) + T*(dcv/dt)
  ! applied to SoilFluxesMod.F90's own errsoi_patch formula (whose
  ! ΔT/fact terms are exactly the cv*(dT/dt) part, evaluated with cv held
  ! at whatever mass was current for that call):
  !   d(hc_soisno)/dt = eflx_soil_grnd - xmf - xmf_h2osfc
  !        - frac_h2osfc*(t_h2osfc-t_h2osfc_bef)*(c_h2osfc/dtime)
  !        + eflx_h2osfc_to_snow + eflx_building_heat
  !        + tssbef*(cv-cv_bef)/dtime - errsoi_col
  ! grnd and phase use only already-persistent veg_ef/col_ef/col_es/col_ws
  ! fields. masschg (the T*(dcv/dt) term above) needed a small, purely
  ! additive change to SoilTemperatureMod.F90/LakeTemperatureMod.F90: a new
  ! col_es%cv_bef snapshot of layer heat capacity from the previous call,
  ! and col_ef%eflx_hc_masschg computed from it and the existing tssbef.
  ! It captures how hc_soisno changes purely from water mass changing
  ! (infiltration, drainage, ET) while referenced to absolute Kelvin
  ! temperature -- not itself a boundary flux, but the dominant missing
  ! piece once phase is accounted for. LakeTemperatureMod folds its own
  ! analogous correction into eflx_sh_tot instead of keeping it separate,
  ! and masschg does not cover the lake-water layers themselves (cv_lake),
  ! so g_phase and g_masschg both under-represent lake columns specifically
  ! -- any residual for a lake-containing gridcell shows up against
  ! errsoi_col, not as a bug in this table.

  integer, parameter :: g_grnd    = 1
  integer, parameter :: g_phase   = 2
  integer, parameter :: g_masschg = 3

  integer, parameter, public :: g_size = g_masschg

  character(len=12),parameter :: gname(g_size) = &
       (/&
       '        grnd', &
       '       phase', &
       '     masschg'  &
       /)

  real(r8) :: budg_gfluxL(g_size, p_size) ! local sum, valid on all pes
  real(r8) :: budg_gfluxG(g_size, p_size) ! global sum, valid only on root pe
  real(r8) :: budg_gfluxN(g_size, p_size) ! counter, valid only on root pe

  !--- S for state ---
  ! Column heat content (col_es%hc_soisno = snow+soil+lake heat content,
  ! MJ/m2) aggregated to the gridcell, split into two reservoirs the same
  ! way col_es%hc_soi already splits it internally: "Soil" (hc_soi, j>=1
  ! layers only) and "Snow+Lake" (hc_soisno-hc_soi, derived at print time --
  ! snow layers for non-lake columns, snow-on-lake plus the lake water
  ! itself for lake columns; the model has no field that separates those
  ! two for lake columns specifically). Plus the soil/lake solver's own
  ! column-level numerical-closure residual (col_ef%errsoi, W/m2),
  ! aggregated the same way. Canopy heat storage is not tracked anywhere
  ! in the model, so there is deliberately no canopy reservoir here
  ! (unlike a fabricated always-zero one).

  integer, parameter :: s_h_beg     = 1
  integer, parameter :: s_h_end     = 2
  integer, parameter :: s_hsoi_beg  = 3
  integer, parameter :: s_hsoi_end  = 4
  integer, parameter :: s_h_errsoi  = 5

  integer, parameter, public :: s_size = s_h_errsoi

  character(len=12),parameter :: sname(s_size) = &
       (/&
       'total_hc_beg', &
       'total_hc_end', &
       ' soil_hc_beg', &
       ' soil_hc_end', &
       '  errsoi_col'  &
       /)

  real(r8) :: budg_stateL(s_size, p_size) ! local sum, valid on all pes
  real(r8), public :: budg_stateG(s_size, p_size) ! global sum, valid only on root pe

  !----- formats -----
  character(*),parameter :: FA0= "('    ',12x,(3x,a10,2x),' | ',(3x,a10,2x))"
  character(*),parameter :: FF = "('    ',a12,f15.8,' | ',f18.2)"
  character(*),parameter :: HS0= "('    ',12x,3(a18),' | ',(a18))"
  character(*),parameter :: HS = "('    ',a12,2(f18.2),18x,' | ',(f18.2))"
  character(*),parameter :: HS2= "('    ',a12,18x,f18.2,18x,' | ',f18.2)"
  character(*),parameter :: HS3= "('    ',a12,3(f18.2),' | ',(f18.2))"

contains

  !-----------------------------------------------------------------------
  subroutine HeatBudget_Reset(mode)
    !
    use elm_time_manager, only : get_curr_date, get_prev_date
    !
    implicit none
    !
    character(len=*), intent(in),optional :: mode
    !
    integer :: year, mon, day, sec
    integer :: ip
    character(*),parameter :: subName = '(HeatBudget_Reset) '

    if (.not.present(mode)) then
       call get_curr_date(year, mon, day, sec)
       call get_prev_date(year, mon, day, sec)

       do ip = 1,p_size
          if (ip == p_inst) then
             budg_fluxL(:,ip)  = 0.0_r8
             budg_fluxG(:,ip)  = 0.0_r8
             budg_fluxN(:,ip)  = 0.0_r8
             budg_stateL(:,ip) = 0.0_r8
             budg_stateG(:,ip) = 0.0_r8
             budg_gfluxL(:,ip) = 0.0_r8
             budg_gfluxG(:,ip) = 0.0_r8
             budg_gfluxN(:,ip) = 0.0_r8
          endif
          if (ip==p_day .and. sec==0) then
             budg_fluxL(:,ip)  = 0.0_r8
             budg_fluxG(:,ip)  = 0.0_r8
             budg_fluxN(:,ip)  = 0.0_r8
             budg_stateL(:,ip) = 0.0_r8
             budg_stateG(:,ip) = 0.0_r8
             budg_gfluxL(:,ip) = 0.0_r8
             budg_gfluxG(:,ip) = 0.0_r8
             budg_gfluxN(:,ip) = 0.0_r8
          endif
          if (ip==p_mon .and. day==1 .and. sec==0) then
             budg_fluxL(:,ip)  = 0.0_r8
             budg_fluxG(:,ip)  = 0.0_r8
             budg_fluxN(:,ip)  = 0.0_r8
             budg_stateL(:,ip) = 0.0_r8
             budg_stateG(:,ip) = 0.0_r8
             budg_gfluxL(:,ip) = 0.0_r8
             budg_gfluxG(:,ip) = 0.0_r8
             budg_gfluxN(:,ip) = 0.0_r8
          endif
          if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
             budg_fluxL(:,ip)  = 0.0_r8
             budg_fluxG(:,ip)  = 0.0_r8
             budg_fluxN(:,ip)  = 0.0_r8
             budg_stateL(:,ip) = 0.0_r8
             budg_stateG(:,ip) = 0.0_r8
             budg_gfluxL(:,ip) = 0.0_r8
             budg_gfluxG(:,ip) = 0.0_r8
             budg_gfluxN(:,ip) = 0.0_r8
          endif
       enddo

    else

       if (trim(mode) == 'inst') then
          budg_fluxL  (:,p_inst)   = 0.0_r8
          budg_fluxG  (:,p_inst)   = 0.0_r8
          budg_fluxN  (:,p_inst)   = 0.0_r8
          budg_stateL (:,p_inst)   = 0.0_r8
          budg_stateG (:,p_inst)   = 0.0_r8
          budg_gfluxL (:,p_inst)   = 0.0_r8
          budg_gfluxG (:,p_inst)   = 0.0_r8
          budg_gfluxN (:,p_inst)   = 0.0_r8
       elseif (trim(mode) == 'day') then
          budg_fluxL  (:,p_day)    = 0.0_r8
          budg_fluxG  (:,p_day)    = 0.0_r8
          budg_fluxN  (:,p_day)    = 0.0_r8
          budg_stateL (:,p_day)    = 0.0_r8
          budg_stateG (:,p_day)    = 0.0_r8
          budg_gfluxL (:,p_day)    = 0.0_r8
          budg_gfluxG (:,p_day)    = 0.0_r8
          budg_gfluxN (:,p_day)    = 0.0_r8
       elseif (trim(mode) == 'mon') then
          budg_fluxL  (:,p_mon)    = 0.0_r8
          budg_fluxG  (:,p_mon)    = 0.0_r8
          budg_fluxN  (:,p_mon)    = 0.0_r8
          budg_stateL (:,p_mon)    = 0.0_r8
          budg_stateG (:,p_mon)    = 0.0_r8
          budg_gfluxL (:,p_mon)    = 0.0_r8
          budg_gfluxG (:,p_mon)    = 0.0_r8
          budg_gfluxN (:,p_mon)    = 0.0_r8
       elseif (trim(mode) == 'ann') then
          budg_fluxL  (:,p_ann)    = 0.0_r8
          budg_fluxG  (:,p_ann)    = 0.0_r8
          budg_fluxN  (:,p_ann)    = 0.0_r8
          budg_stateL (:,p_ann)    = 0.0_r8
          budg_stateG (:,p_ann)    = 0.0_r8
          budg_gfluxL (:,p_ann)    = 0.0_r8
          budg_gfluxG (:,p_ann)    = 0.0_r8
          budg_gfluxN (:,p_ann)    = 0.0_r8
       elseif (trim(mode) == 'inf') then
          budg_fluxL  (:,p_inf)    = 0.0_r8
          budg_fluxG  (:,p_inf)    = 0.0_r8
          budg_fluxN  (:,p_inf)    = 0.0_r8
          budg_stateL (:,p_inf)    = 0.0_r8
          budg_stateG (:,p_inf)    = 0.0_r8
          budg_gfluxL (:,p_inf)    = 0.0_r8
          budg_gfluxG (:,p_inf)    = 0.0_r8
          budg_gfluxN (:,p_inf)    = 0.0_r8
       elseif (trim(mode) == 'all') then
          budg_fluxL  (:,:)        = 0.0_r8
          budg_fluxG  (:,:)        = 0.0_r8
          budg_fluxN  (:,:)        = 0.0_r8
          budg_stateL (:,:)        = 0.0_r8
          budg_stateG (:,:)        = 0.0_r8
          budg_gfluxL (:,:)        = 0.0_r8
          budg_gfluxG (:,:)        = 0.0_r8
          budg_gfluxN (:,:)        = 0.0_r8
       else
          call shr_sys_abort(subname//' ERROR in mode '//trim(mode))
       endif
    endif

  end subroutine HeatBudget_Reset

  !-----------------------------------------------------------------------
  subroutine HeatBudget_Accum()
    !
    use elm_time_manager, only : get_curr_date, get_prev_date, get_nstep
    !
    implicit none
    !
    integer                :: ip, nf
    integer                :: year_prev, month_prev, day_prev, sec_prev
    integer                :: year_curr, month_curr, day_curr, sec_curr
    character(*),parameter :: subName = '(HeatBudget_Accum)'
    logical                :: update_state_beg, update_state_end

    call get_prev_date(year_prev, month_prev, day_prev, sec_prev)
    call get_curr_date(year_curr, month_curr, day_curr, sec_curr)

    do ip = p_inst+1, p_size
       budg_fluxL(:,ip)  = budg_fluxL(:,ip)  + budg_fluxL(:,p_inst)
       budg_gfluxL(:,ip) = budg_gfluxL(:,ip) + budg_gfluxL(:,p_inst)
       update_state_beg = .false.
       update_state_end = .false.

       select case (ip)
       case (p_day)
          if (sec_prev == 0) update_state_beg = .true.
          if (sec_curr == 0) update_state_end = .true.
       case (p_mon)
          if (sec_prev == 0 .and. day_prev == 1) update_state_beg = .true.
          if (sec_curr == 0 .and. day_curr == 1) update_state_end = .true.
       case (p_ann)
          if (sec_prev == 0 .and. day_prev == 1 .and. month_prev == 1) update_state_beg = .true.
          if (sec_curr == 0 .and. day_curr == 1 .and. month_curr == 1) update_state_end = .true.
       case (p_inf)
          update_state_end = .true.
       end select
       if (get_nstep() == 1) update_state_beg = .true.

       if (update_state_beg) then
          nf = s_h_beg    ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
          nf = s_hsoi_beg ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
       endif

       if (update_state_end) then
          nf = s_h_end    ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
          nf = s_hsoi_end ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
       endif
       nf = s_h_errsoi ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + budg_stateL(nf, p_inst)
    end do
    budg_fluxN(:,:)  = budg_fluxN(:,:)  + 1._r8
    budg_gfluxN(:,:) = budg_gfluxN(:,:) + 1._r8

  end subroutine HeatBudget_Accum

  !-----------------------------------------------------------------------
  subroutine HeatBudget_SetBeginningStates(bounds)
    !
    ! !DESCRIPTION:
    ! Snapshot the grid-level column heat content (snow+soil+lake,
    ! col_es%hc_soisno) at the beginning of the time step, before
    ! dynSubgrid_driver reweights columns for this step. This mirrors
    ! BalanceCheckMod's BeginGridWaterBalance, which captures grc_ws%begwb
    ! at the same point in elm_driver.F90 for exactly the same reason:
    ! taking the snapshot any later would double-count (or miss) the
    ! land-cover reweighting jump that eflx_dynbal already corrects for
    ! in the flux table.
    !
    ! col_es%hc_soisno is a diagnostic (not a prognostic state): it holds
    ! spval until the first call to SoilTemperature/LakeTemperature this
    ! run, so grc_es%beg_hc will legitimately be spval on the very first
    ! timestep (and over any all-urban gridcell, since urban columns are
    ! never included in hc_soisno). HeatBudget_Run skips accumulating the
    ! state terms wherever that happens.
    !
    use subgridAveMod, only : c2g
    !
    implicit none
    !
    type(bounds_type), intent(in) :: bounds

    call c2g(bounds, col_es%hc_soisno(bounds%begc:bounds%endc), &
         grc_es%beg_hc(bounds%begg:bounds%endg), &
         c2l_scale_type='unity', l2g_scale_type='unity')

    call c2g(bounds, col_es%hc_soi(bounds%begc:bounds%endc), &
         grc_es%beg_hc_soi(bounds%begg:bounds%endg), &
         c2l_scale_type='unity', l2g_scale_type='unity')

  end subroutine HeatBudget_SetBeginningStates

  !-----------------------------------------------------------------------
  subroutine HeatBudget_SetEndingStates(bounds)
    !
    ! !DESCRIPTION:
    ! Snapshot the grid-level column heat content (col_es%hc_soisno) and
    ! the soil/lake solver's own numerical-closure residual
    ! (col_ef%errsoi) at the end of the time step, once
    ! SoilTemperature/SoilFluxes/LakeTemperature have all run for this
    ! step. Called per-clump alongside GridBalanceCheck/
    ! WaterBudget_SetEndingMonthlyStates, matching where those grid-level
    ! c2g aggregations already happen in elm_driver.F90 -- HeatBudget_Run
    ! itself (called once, on bounds_proc, right before lnd2atm_vars-based
    ! fluxes are gathered) only reads the results, mirroring
    ! WaterBudget_Run's own read-only use of grc_ws%endwb.
    !
    ! Also aggregate the two terms behind the "into ground" flux table:
    ! grc_ef%heat_into_grnd (p2g of veg_ef%eflx_soil_grnd) and
    ! grc_ef%heat_phase_corr (c2g of a per-column combination of the
    ! phase-change and h2osfc/building-heat terms from
    ! SoilFluxesMod.F90's own errsoi_patch formula -- see HeatBudget_Run's
    ! header for the exact expression and its provenance).
    !
    use subgridAveMod, only : c2g, p2g
    !
    implicit none
    !
    type(bounds_type), intent(in) :: bounds
    !
    integer  :: c
    real(r8) :: phase_col(bounds%begc:bounds%endc)

    call c2g(bounds, col_es%hc_soisno(bounds%begc:bounds%endc), &
         grc_es%end_hc(bounds%begg:bounds%endg), &
         c2l_scale_type='unity', l2g_scale_type='unity')

    call c2g(bounds, col_es%hc_soi(bounds%begc:bounds%endc), &
         grc_es%end_hc_soi(bounds%begg:bounds%endg), &
         c2l_scale_type='unity', l2g_scale_type='unity')

    call c2g(bounds, col_ef%errsoi(bounds%begc:bounds%endc), &
         grc_es%errsoi(bounds%begg:bounds%endg), &
         c2l_scale_type='unity', l2g_scale_type='unity')

    call p2g(bounds, veg_ef%eflx_soil_grnd(bounds%begp:bounds%endp), &
         grc_ef%heat_into_grnd(bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type='unity', l2g_scale_type='unity')

    associate(                                       &
         xmf          => col_ef%xmf                , &
         xmf_h2osfc   => col_ef%xmf_h2osfc         , &
         h2osfc2snow  => col_ef%eflx_h2osfc_to_snow, &
         bldg_heat    => col_ef%eflx_building_heat , &
         frac_h2osfc  => col_ws%frac_h2osfc        , &
         t_h2osfc     => col_es%t_h2osfc           , &
         t_h2osfc_bef => col_es%t_h2osfc_bef       , &
         c_h2osfc     => col_es%c_h2osfc             &
         )

      do c = bounds%begc, bounds%endc
         if (xmf(c) == spval .or. xmf_h2osfc(c) == spval .or. h2osfc2snow(c) == spval &
              .or. bldg_heat(c) == spval .or. frac_h2osfc(c) == spval &
              .or. t_h2osfc(c) == spval .or. t_h2osfc_bef(c) == spval .or. c_h2osfc(c) == spval) then
            phase_col(c) = spval
         else
            phase_col(c) = -xmf(c) - xmf_h2osfc(c) &
                 - frac_h2osfc(c)*(t_h2osfc(c)-t_h2osfc_bef(c))*(c_h2osfc(c)/dtime_mod) &
                 + h2osfc2snow(c) + bldg_heat(c)
         end if
      end do

    end associate

    call c2g(bounds, phase_col(bounds%begc:bounds%endc), &
         grc_ef%heat_phase_corr(bounds%begg:bounds%endg), &
         c2l_scale_type='unity', l2g_scale_type='unity')

    call c2g(bounds, col_ef%eflx_hc_masschg(bounds%begc:bounds%endc), &
         grc_ef%heat_masschg(bounds%begg:bounds%endg), &
         c2l_scale_type='unity', l2g_scale_type='unity')

  end subroutine HeatBudget_SetEndingStates

  !-----------------------------------------------------------------------
  subroutine HeatBudget_Run(bounds, atm2lnd_vars, lnd2atm_vars)
    !
    ! !DESCRIPTION:
    ! Accumulate the instantaneous, area-weighted net heat fluxes
    ! exchanged between ELM and the atmosphere. Mirrors WaterBudget_Run's
    ! gridcell loop and area weighting exactly.
    !
    ! Also accumulate the grid-level heat-content state (col_es%hc_soisno,
    ! now holding this step's ending value since SoilTemperature/
    ! SoilFluxes/LakeTemperature have already run) and the soil/lake
    ! solver's own numerical-closure residual (col_ef%errsoi). Note that
    ! (end_hc - beg_hc) - (net flux table * dt) is NOT a clean numerical
    ! truncation term the way water's errh2o is: it also contains real,
    ! currently-unmodeled advected heat carried by precipitation,
    ! snowmelt, and runoff (see HEAT_BUDGET_PLAN.md, Phase 2 issue #3).
    ! errsoi_col itself, by contrast, IS a validated, tightly-bounded
    ! (<1e-5 W/m2, or the model aborts) internal closure term -- it is
    ! printed here for that reason, not as a stand-in for the former.
    !
    use domainMod, only : ldomain
    use elm_varcon, only : re
    !
    implicit none

    type(bounds_type)        , intent(in) :: bounds
    type(atm2lnd_type)       , intent(in) :: atm2lnd_vars
    type(lnd2atm_type)       , intent(in) :: lnd2atm_vars

    integer  :: g, nf, ip
    real(r8) :: af, one_over_re2

    associate(                                                              &
         forc_lwdn   => atm2lnd_vars%forc_lwrad_not_downscaled_grc        , &
         netsw       => lnd2atm_vars%fsa_grc                              , &
         lwup        => lnd2atm_vars%eflx_lwrad_out_grc                   , &
         eflx_sh_tot => lnd2atm_vars%eflx_sh_tot_grc                      , &
         eflx_lh_tot => lnd2atm_vars%eflx_lh_tot_grc                      , &
         beg_hc_grc  => grc_es%beg_hc                                     , &
         end_hc_grc  => grc_es%end_hc                                     , &
         beg_hcsoi_grc => grc_es%beg_hc_soi                               , &
         end_hcsoi_grc => grc_es%end_hc_soi                               , &
         errsoi_grc  => grc_es%errsoi                                     , &
         grnd_grc    => grc_ef%heat_into_grnd                             , &
         phase_grc   => grc_ef%heat_phase_corr                            , &
         masschg_grc => grc_ef%heat_masschg                                 &
         )

      ip = p_inst

      budg_stateL(:,ip) = 0._r8
      budg_gfluxL(:,ip) = 0._r8
      one_over_re2 = 1._r8/(re**2._r8)

      do g = bounds%begg, bounds%endg

         af   = (ldomain%area(g) * one_over_re2) * & ! area (converting km**2 to radians**2)
                ldomain%frac(g)                      ! land fraction

         nf = f_netsw; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) + netsw(g)*af
         nf = f_lwdn ; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) + forc_lwdn(g)*af
         nf = f_lwup ; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) - lwup(g)*af
         nf = f_lh   ; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) - eflx_lh_tot(g)*af
         nf = f_sh   ; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) - eflx_sh_tot(g)*af

         if (beg_hc_grc(g) /= spval .and. end_hc_grc(g) /= spval) then
            nf = s_h_beg ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + beg_hc_grc(g)*af
            nf = s_h_end ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + end_hc_grc(g)*af
         end if
         if (beg_hcsoi_grc(g) /= spval .and. end_hcsoi_grc(g) /= spval) then
            nf = s_hsoi_beg ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + beg_hcsoi_grc(g)*af
            nf = s_hsoi_end ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + end_hcsoi_grc(g)*af
         end if
         if (errsoi_grc(g) /= spval) then
            nf = s_h_errsoi ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + errsoi_grc(g)*af
         end if

         if (grnd_grc(g) /= spval) then
            nf = g_grnd ; budg_gfluxL(nf,ip) = budg_gfluxL(nf,ip) + grnd_grc(g)*af
         end if
         if (phase_grc(g) /= spval) then
            nf = g_phase ; budg_gfluxL(nf,ip) = budg_gfluxL(nf,ip) + phase_grc(g)*af
         end if
         if (masschg_grc(g) /= spval) then
            nf = g_masschg ; budg_gfluxL(nf,ip) = budg_gfluxL(nf,ip) + masschg_grc(g)*af
         end if
      end do

    end associate

  end subroutine HeatBudget_Run

  !-----------------------------------------------------------------------
  subroutine HeatBudget_Sum0()
    !
    use spmdMod    , only : mpicom
    use shr_mpi_mod, only : shr_mpi_sum
    !
    implicit none
    !
    real(r8)               :: budg_fluxGtmp(f_size,p_size) ! temporary sum
    real(r8)               :: budg_stateGtmp(s_size,p_size) ! temporary sum
    real(r8)               :: budg_gfluxGtmp(g_size,p_size) ! temporary sum
    character(*),parameter :: subName = '(HeatBudget_Sum0)'

    budg_fluxGtmp = 0._r8
    budg_stateGtmp = 0._r8
    budg_gfluxGtmp = 0._r8

    call shr_mpi_sum(budg_fluxL, budg_fluxGtmp, mpicom, subName, all=.true. )
    call shr_mpi_sum(budg_stateL, budg_stateGtmp, mpicom, subName, all=.true. )
    call shr_mpi_sum(budg_gfluxL, budg_gfluxGtmp, mpicom, subName, all=.true. )

    budg_fluxG  = budg_fluxG + budg_fluxGtmp
    budg_stateG = budg_stateGtmp
    budg_gfluxG = budg_gfluxG + budg_gfluxGtmp

    budg_fluxL            = 0._r8 ! reset all fluxes
    budg_stateL(:,p_inst) = 0._r8 ! only reset instantaneous states
    budg_gfluxL            = 0._r8 ! reset all fluxes

  end subroutine HeatBudget_Sum0

  !-----------------------------------------------------------------------
  subroutine HeatBudget_Print(budg_print_inst,  budg_print_daily,  budg_print_month,  &
       budg_print_ann,  budg_print_ltann,  budg_print_ltend)
    !
    use elm_time_manager, only : get_curr_date, get_prev_date, get_nstep, get_step_size
    use shr_const_mod   , only : shr_const_pi
    !
    implicit none
    !
    integer , intent(in) :: budg_print_inst
    integer , intent(in) :: budg_print_daily
    integer , intent(in) :: budg_print_month
    integer , intent(in) :: budg_print_ann
    integer , intent(in) :: budg_print_ltann
    integer , intent(in) :: budg_print_ltend
    !
    ! !LOCAL VARIABLES:
    integer :: f,ip ! data array indicies
    integer :: plev        ! print level
    integer :: year, mon, day, sec
    integer :: cdate
    logical :: sumdone
    real(r8) :: unit_conversion
    real(r8) :: budg_fluxGpr (f_size,p_size) ! values to print, scaled and such
    real(r8) :: budg_gfluxGpr (g_size,p_size) ! values to print, scaled and such

    sumdone = .false.

    if (get_nstep() <= 1) then
       call get_prev_date(year, mon, day, sec);
    else
       call get_curr_date(year, mon, day, sec);
    end if

    cdate = year*10000 + mon*100 + day

    do ip = 1,p_size
       plev = 0
       if (ip == p_inst) then
          plev = max(plev,budg_print_inst)
       endif
       if (ip==p_day .and. sec==0) then
          plev = max(plev,budg_print_daily)
       endif
       if (ip==p_mon .and. day==1 .and. sec==0) then
          plev = max(plev,budg_print_month)
       endif
       if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
          plev = max(plev,budg_print_ann)
       endif
       if (ip==p_inf .and. mon==1 .and. day==1 .and. sec==0) then
          plev = max(plev,budg_print_ltann)
       endif

       if (plev > 0) then
          unit_conversion = 1.d0/(4.0_r8*shr_const_pi)*1.0e6_r8
          if (.not.sumdone) then
             sumdone = .true.
             call HeatBudget_Sum0()
             budg_fluxGpr = budg_fluxG
             budg_fluxGpr = budg_fluxGpr*unit_conversion
             budg_fluxGpr = budg_fluxGpr/budg_fluxN
             budg_gfluxGpr = budg_gfluxG
             budg_gfluxGpr = budg_gfluxGpr*unit_conversion
             budg_gfluxGpr = budg_gfluxGpr/budg_gfluxN
          end if

          if (ip == p_day .and. get_nstep() == 1) cycle
          if (ip == p_mon .and. get_nstep() == 1) cycle
          if (ip == p_ann .and. get_nstep() == 1) cycle
          if (ip == p_inf .and. get_nstep() == 1) cycle

          if (masterproc) then
             write(iulog,*)''
             write(iulog,*)'NET HEAT FLUXES : period ',trim(pname(ip)),': date = ',cdate,sec
             write(iulog,FA0)'  Time  ','  Time    '
             write(iulog,FA0)'averaged','integrated'
             write(iulog,FA0)'W/m2*1e6','MJ/m2*1e6'
             write(iulog,'(32("-"),"|",20("-"))')
             do f = 1, f_size
                write(iulog,FF)fname(f),budg_fluxGpr(f,ip),budg_fluxG(f,ip)*unit_conversion*get_step_size()/1.e6_r8
             end do
             write(iulog,'(32("-"),"|",20("-"))')
             write(iulog,FF)'   *SUM*', &
                  sum(budg_fluxGpr(:,ip)), sum(budg_fluxG(:,ip))*unit_conversion*get_step_size()/1.e6_r8
             write(iulog,'(32("-"),"|",20("-"))')

             write(iulog,*)''
             write(iulog,*)'NET HEAT FLUXES (into ground) : period ',trim(pname(ip)),': date = ',cdate,sec
             write(iulog,FA0)'  Time  ','  Time    '
             write(iulog,FA0)'averaged','integrated'
             write(iulog,FA0)'W/m2*1e6','MJ/m2*1e6'
             write(iulog,'(32("-"),"|",20("-"))')
             do f = 1, g_size
                write(iulog,FF)gname(f),budg_gfluxGpr(f,ip),budg_gfluxG(f,ip)*unit_conversion*get_step_size()/1.e6_r8
             end do
             write(iulog,'(32("-"),"|",20("-"))')
             write(iulog,FF)'   *SUM*', &
                  sum(budg_gfluxGpr(:,ip)), sum(budg_gfluxG(:,ip))*unit_conversion*get_step_size()/1.e6_r8
             write(iulog,'(32("-"),"|",20("-"))')
             write(iulog,*)'  (integrated column is directly comparable to the *NET CHANGE* row of HEAT STATES below)'

             write(iulog,*)''
             write(iulog,*)'HEAT STATES (MJ/m2*1e6): period ',trim(pname(ip)),': date = ',cdate,sec
             write(iulog,HS0) &
                  '       Soil       ', &
                  '     Snow+Lake    ', &
                  '   Grid-level Err ', &
                  '       TOTAL      '
             write(iulog,'(71("-"),"|",20("-"))')
             write(iulog,HS) '         beg', &
                  budg_stateG(s_hsoi_beg,ip)*unit_conversion, &
                  (budg_stateG(s_h_beg,ip)-budg_stateG(s_hsoi_beg,ip))*unit_conversion, &
                  budg_stateG(s_h_beg,ip)*unit_conversion
             write(iulog,HS) '         end', &
                  budg_stateG(s_hsoi_end,ip)*unit_conversion, &
                  (budg_stateG(s_h_end,ip)-budg_stateG(s_hsoi_end,ip))*unit_conversion, &
                  budg_stateG(s_h_end,ip)*unit_conversion
             write(iulog,HS3)'*NET CHANGE*', &
                  (budg_stateG(s_hsoi_end,ip)-budg_stateG(s_hsoi_beg,ip))*unit_conversion, &
                  ((budg_stateG(s_h_end,ip)-budg_stateG(s_hsoi_end,ip)) &
                  -(budg_stateG(s_h_beg,ip)-budg_stateG(s_hsoi_beg,ip)))*unit_conversion, &
                  budg_stateG(s_h_errsoi,ip)*unit_conversion/budg_fluxN(1,ip), &
                  (budg_stateG(s_h_end,ip)-budg_stateG(s_h_beg,ip))*unit_conversion
             write(iulog,'(71("-"),"|",20("-"))')
             write(iulog,HS2)'   *SUM*    ', &
                  (budg_stateG(s_h_end,ip)-budg_stateG(s_h_beg,ip))*unit_conversion &
                  - budg_stateG(s_h_errsoi,ip)*unit_conversion/budg_fluxN(1,ip), &
                  (budg_stateG(s_h_end,ip)-budg_stateG(s_h_beg,ip))*unit_conversion
             write(iulog,'(71("-"),"|",20("-"))')
          end if
       end if
    end do

  end subroutine HeatBudget_Print

  !-----------------------------------------------------------------------
  subroutine HeatBudget_Restart(bounds, ncid, flag)
    !
    use ncdio_pio, only : file_desc_t, ncd_io, ncd_double, ncd_int
    use ncdio_pio, only : ncd_defvar
    !
    implicit none
    !
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    character(len=*),parameter :: subname = 'HeatBudget_Restart'

    select case (trim(flag))
    case ('define')
       call HeatBudget_Restart_Define(bounds, ncid)
    case ('write')
       call HeatBudget_Restart_Write(bounds, ncid, flag)
    case ('read')
       call HeatBudget_Restart_Read(bounds, ncid, flag)
    case default
       write(iulog,*) trim(subname),' ERROR: unknown flag = ',flag
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine HeatBudget_Restart

  !-----------------------------------------------------------------------
  subroutine HeatBudget_Restart_Define(bounds, ncid)
    !
    use ncdio_pio, only : file_desc_t, ncd_io, ncd_double, ncd_defvar
    !
    implicit none
    !
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id

    call ncd_defvar(varname='heat_budg_fluxG', xtype=ncd_double, &
         dim1name='heat_budg_flux', &
         long_name='heat_budg_fluxG', units='W/m2', ncid=ncid)

    call ncd_defvar(varname='heat_budg_fluxN', xtype=ncd_double, &
         dim1name='heat_budg_flux', &
         long_name='heat_budg_fluxN', units='-', ncid=ncid)

    call ncd_defvar(varname='heat_budg_stateG', xtype=ncd_double, &
         dim1name='heat_budg_state', &
         long_name='heat_budg_stateG', units='MJ/m2 or W/m2', ncid=ncid)

    call ncd_defvar(varname='heat_budg_gfluxG', xtype=ncd_double, &
         dim1name='heat_budg_gflux', &
         long_name='heat_budg_gfluxG', units='W/m2', ncid=ncid)

    call ncd_defvar(varname='heat_budg_gfluxN', xtype=ncd_double, &
         dim1name='heat_budg_gflux', &
         long_name='heat_budg_gfluxN', units='-', ncid=ncid)

  end subroutine HeatBudget_Restart_Define

  !-----------------------------------------------------------------------
  subroutine HeatBudget_Restart_Write(bounds, ncid, flag)
    !
    use ncdio_pio   , only : file_desc_t, ncd_io, ncd_double, ncd_int
    use ncdio_pio   , only : ncd_defvar
    use spmdMod     , only : mpicom
    use shr_mpi_mod , only : shr_mpi_sum
    !
    implicit none
    !
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    real(r8) :: budg_fluxGtmp(f_size,p_size) ! temporary sum
    real(r8) :: budg_stateGtmp(s_size,p_size) ! temporary sum
    real(r8) :: budg_gfluxGtmp(g_size,p_size) ! temporary sum
    real(r8) :: budg_fluxG_1D (f_size*p_size)
    real(r8) :: budg_fluxN_1D (f_size*p_size)
    real(r8) :: budg_stateG_1D(s_size*p_size)
    real(r8) :: budg_gfluxG_1D(g_size*p_size)
    real(r8) :: budg_gfluxN_1D(g_size*p_size)
    integer  :: f, s, p, count
    character(*),parameter :: subName = '(HeatBudget_Restart_Write) '

    budg_fluxGtmp = 0._r8
    budg_stateGtmp = 0._r8
    budg_gfluxGtmp = 0._r8

    call shr_mpi_sum(budg_fluxL, budg_fluxGtmp, mpicom, subName, all=.true.)
    call shr_mpi_sum(budg_stateL, budg_stateGtmp, mpicom, subName, all=.true.)
    call shr_mpi_sum(budg_gfluxL, budg_gfluxGtmp, mpicom, subName, all=.true.)

    ! Copy data from 2D into 1D array
    count = 0
    do f = 1, f_size
       do p = 1, p_size
          count = count + 1
          budg_fluxG_1D(count) = budg_fluxG(f,p) + budg_fluxGtmp(f,p)
          budg_fluxN_1D(count) = budg_fluxN(f,p)
       end do
    end do

    ! Copy data from 2D into 1D array
    count = 0
    do s = 1, s_size
       do p = 1, p_size
          count = count + 1
          budg_stateG_1D(count) = budg_stateGtmp(s,p)
       end do
    end do

    ! Copy data from 2D into 1D array
    count = 0
    do f = 1, g_size
       do p = 1, p_size
          count = count + 1
          budg_gfluxG_1D(count) = budg_gfluxG(f,p) + budg_gfluxGtmp(f,p)
          budg_gfluxN_1D(count) = budg_gfluxN(f,p)
       end do
    end do

    call ncd_io(flag=flag, varname='heat_budg_fluxG', data=budg_fluxG_1D, ncid=ncid)
    call ncd_io(flag=flag, varname='heat_budg_fluxN', data=budg_fluxN_1D, ncid=ncid)
    call ncd_io(flag=flag, varname='heat_budg_stateG', data=budg_stateG_1D, ncid=ncid)
    call ncd_io(flag=flag, varname='heat_budg_gfluxG', data=budg_gfluxG_1D, ncid=ncid)
    call ncd_io(flag=flag, varname='heat_budg_gfluxN', data=budg_gfluxN_1D, ncid=ncid)

  end subroutine HeatBudget_Restart_Write

  !-----------------------------------------------------------------------
  subroutine HeatBudget_Restart_Read(bounds, ncid, flag)
    !
    use ncdio_pio, only : file_desc_t, ncd_io, ncd_double, ncd_int
    use ncdio_pio, only : ncd_defvar
    !
    implicit none
    !
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    real(r8) :: budg_fluxG_1D (f_size*p_size)
    real(r8) :: budg_fluxN_1D (f_size*p_size)
    real(r8) :: budg_stateG_1D(s_size*p_size)
    real(r8) :: budg_gfluxG_1D(g_size*p_size)
    real(r8) :: budg_gfluxN_1D(g_size*p_size)
    integer  :: f, s, p, count
    logical  :: readvar_fluxG, readvar_fluxN, readvar_stateG, readvar_gfluxG, readvar_gfluxN

    ! Every restart file that predates this change will lack these
    ! variables; guard with readvar and zero-initialize the accumulators
    ! rather than aborting or reading garbage, so restarting from an older
    ! case works.
    budg_fluxG_1D = 0._r8
    budg_fluxN_1D = 0._r8
    budg_stateG_1D = 0._r8
    budg_gfluxG_1D = 0._r8
    budg_gfluxN_1D = 0._r8

    call ncd_io(flag=flag, varname='heat_budg_fluxG', data=budg_fluxG_1D, ncid=ncid, readvar=readvar_fluxG)
    call ncd_io(flag=flag, varname='heat_budg_fluxN', data=budg_fluxN_1D, ncid=ncid, readvar=readvar_fluxN)
    call ncd_io(flag=flag, varname='heat_budg_stateG', data=budg_stateG_1D, ncid=ncid, readvar=readvar_stateG)
    call ncd_io(flag=flag, varname='heat_budg_gfluxG', data=budg_gfluxG_1D, ncid=ncid, readvar=readvar_gfluxG)
    call ncd_io(flag=flag, varname='heat_budg_gfluxN', data=budg_gfluxN_1D, ncid=ncid, readvar=readvar_gfluxN)

    if (.not. readvar_fluxG) budg_fluxG_1D = 0._r8
    if (.not. readvar_fluxN) budg_fluxN_1D = 0._r8
    if (.not. readvar_stateG) budg_stateG_1D = 0._r8
    if (.not. readvar_gfluxG) budg_gfluxG_1D = 0._r8
    if (.not. readvar_gfluxN) budg_gfluxN_1D = 0._r8

    ! Copy data from 1D into 2D array
    count = 0
    do f = 1, f_size
       do p = 1, p_size
          count = count + 1
          budg_fluxG(f,p) = budg_fluxG_1D(count)
          budg_fluxN(f,p) = budg_fluxN_1D(count)
       end do
    end do

    ! Copy data from 1D into 2D array
    count = 0
    do f = 1, g_size
       do p = 1, p_size
          count = count + 1
          budg_gfluxG(f,p) = budg_gfluxG_1D(count)
          budg_gfluxN(f,p) = budg_gfluxN_1D(count)
       end do
    end do

    ! Copy data from 1D into 2D array. Only masterproc's budg_stateL is
    ! seeded with the restored sum: HeatBudget_Sum0 recomputes budg_stateG
    ! as an MPI sum of budg_stateL across all procs (an overwrite, not an
    ! accumulation, unlike the flux side), so seeding every proc identically
    ! here would multiply the restored value by the number of MPI ranks.
    if (masterproc) then
       count = 0
       do s = 1, s_size
          do p = 1, p_size
             count = count + 1
             budg_stateL(s,p) = budg_stateG_1D(count)
          end do
       end do
    end if

  end subroutine HeatBudget_Restart_Read

end module HeatBudgetMod
