module WaterBudgetMod
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use shr_sys_mod       , only : shr_sys_abort
  use decompMod         , only : bounds_type
  use abortutils        , only : endrun
  use clm_varctl        , only : iulog
  use atm2lndType       , only : atm2lnd_type
  use lnd2atmType       , only : lnd2atm_type
  use WaterstateType    , only : waterstate_type
  use WaterfluxType     , only : waterflux_type
  use spmdMod           , only : masterproc
  use SoilHydrologyType , only : soilhydrology_type

  implicit none
  save
  private

  public :: WaterBudget_Reset
  public :: WaterBudget_Run
  public :: WaterBudget_Accum
  !public :: WaterBudget_Sum_Local
  public :: WaterBudget_Print
  public :: WaterBudget_Restart
  public :: WaterBudget_SetBeginningMonthlyStates
  public :: WaterBudget_SetEndingMonthlyStates

  !--- F for flux ---

  integer, parameter :: f_rain = 1
  integer, parameter :: f_snow = 2
  integer, parameter :: f_evap = 3
  integer, parameter :: f_roff = 4
  integer, parameter :: f_ioff = 5

  integer, parameter, public :: f_size = f_ioff

  character(len=12),parameter :: fname(f_size) = &
       (/&
       '        rain', &
       '        snow', &
       '        evap', &
       '      runoff', &
       '      frzrof'  &
       /)

  !--- S for state ---
  integer, parameter :: s_w_beg     = 1
  integer, parameter :: s_w_end     = 2
  integer, parameter :: s_wcan_beg  = 3
  integer, parameter :: s_wcan_end  = 4
  integer, parameter :: s_wsno_beg  = 5
  integer, parameter :: s_wsno_end  = 6
  integer, parameter :: s_wsfc_beg  = 7
  integer, parameter :: s_wsfc_end  = 8
  integer, parameter :: s_wsliq_beg = 9
  integer, parameter :: s_wsliq_end = 10
  integer, parameter :: s_wsice_beg = 11
  integer, parameter :: s_wsice_end = 12
  integer, parameter :: s_wwa_beg   = 13
  integer, parameter :: s_wwa_end   = 14
  integer, parameter :: s_w_errh2o  = 15

  integer, parameter, public :: s_size = s_w_errh2o

  character(len=12),parameter :: sname(s_size) = &
       (/&
       ' total_w_beg', &
       ' total_w_end', &
       'canopy_w_beg', &
       'canopy_w_end', &
       '  snow_w_beg', &
       '  snow_w_end', &
       '   sfc_w_beg', &
       '   sfc_w_end', &
       'soiliq_w_beg', &
       'soiliq_w_end', &
       'soiice_w_beg', &
       'soiice_w_end', &
       ' aquif_w_beg', &
       ' aquif_w_end', &
       '   error h2o'  &
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

  real(r8) :: budg_fluxL(f_size, p_size)
  real(r8) :: budg_fluxG(f_size, p_size)
  real(r8) :: budg_fluxN(f_size, p_size)

  real(r8) :: budg_stateL(s_size, p_size)
  real(r8), public :: budg_stateG(s_size, p_size)

  logical,save :: first_time = .true.

  !----- formats -----
  character(*),parameter :: FA0= "('    ',12x,(3x,a10,2x),' | ',(3x,a10,2x))"
  character(*),parameter :: FF = "('    ',a12,f15.8,' | ',f18.2)"
  character(*),parameter :: FF2= "('    ',a12,a15,' | ',f18.2)"
  character(*),parameter :: FS = "('    ',a12,6(f18.2),18x,' | ',(f18.2))"
  character(*),parameter :: FS0= "('    ',12x,7(a18),' | ',(a18))"
  character(*),parameter :: FS2= "('    ',a12,54x,f18.2,54x,' | ',f18.2)"
  character(*),parameter :: FS3= "('    ',a12,7(f18.2),' | ',(f18.2))"

contains

  !-----------------------------------------------------------------------
  subroutine WaterBudget_Reset(mode)
    !
    use clm_time_manager, only : get_curr_date, get_prev_date
    !
    implicit none
    !
    character(len=*), intent(in),optional :: mode
    !
    integer :: year, mon, day, sec
    integer :: ip
    character(*),parameter :: subName = '(WaterBudget_Reset) '

    if (.not.present(mode)) then
       call get_curr_date(year, mon, day, sec)
       call get_prev_date(year, mon, day, sec)

       do ip = 1,p_size
          if (ip == p_inst) then
             budg_fluxL(:,ip)  = 0.0_r8
             budg_fluxG(:,ip)  = 0.0_r8
             budg_fluxN(:,ip)  = 0.0_r8
          endif
          if (ip==p_day .and. sec==0) then
             budg_fluxL(:,ip)  = 0.0_r8
             budg_fluxG(:,ip)  = 0.0_r8
             budg_fluxN(:,ip)  = 0.0_r8
          endif
          if (ip==p_mon .and. day==1 .and. sec==0) then
             budg_fluxL(:,ip)  = 0.0_r8
             budg_fluxG(:,ip)  = 0.0_r8
             budg_fluxN(:,ip)  = 0.0_r8
          endif
          if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
             budg_fluxL(:,ip)  = 0.0_r8
             budg_fluxG(:,ip)  = 0.0_r8
             budg_fluxN(:,ip)  = 0.0_r8
          endif
       enddo

    else

       if (trim(mode) == 'inst') then
          budg_fluxL  (:,p_inst)   = 0.0_r8
          budg_fluxG  (:,p_inst)   = 0.0_r8
          budg_fluxN  (:,p_inst)   = 0.0_r8
          budg_stateL (:,p_inst)   = 0.0_r8
          budg_stateG (:,p_inst)   = 0.0_r8
       elseif (trim(mode) == 'day') then
          budg_fluxL  (:,p_day)    = 0.0_r8
          budg_fluxG  (:,p_day)    = 0.0_r8
          budg_fluxN  (:,p_day)    = 0.0_r8
          budg_stateL (:,p_day)    = 0.0_r8
          budg_stateG (:,p_day)    = 0.0_r8
       elseif (trim(mode) == 'mon') then
          budg_fluxL  (:,p_mon)    = 0.0_r8
          budg_fluxG  (:,p_mon)    = 0.0_r8
          budg_fluxN  (:,p_mon)    = 0.0_r8
          budg_stateL (:,p_mon)    = 0.0_r8
          budg_stateG (:,p_mon)    = 0.0_r8
       elseif (trim(mode) == 'ann') then
          budg_fluxL  (:,p_ann)    = 0.0_r8
          budg_fluxG  (:,p_ann)    = 0.0_r8
          budg_fluxN  (:,p_ann)    = 0.0_r8
          budg_stateL (:,p_ann)    = 0.0_r8
          budg_stateG (:,p_ann)    = 0.0_r8
       elseif (trim(mode) == 'inf') then
          budg_fluxL  (:,p_inf)    = 0.0_r8
          budg_fluxG  (:,p_inf)    = 0.0_r8
          budg_fluxN  (:,p_inf)    = 0.0_r8
          budg_stateL (:,p_inf)    = 0.0_r8
          budg_stateG (:,p_inf)    = 0.0_r8
       elseif (trim(mode) == 'all') then
          budg_fluxL  (:,:)        = 0.0_r8
          budg_fluxG  (:,:)        = 0.0_r8
          budg_fluxN  (:,:)        = 0.0_r8
          budg_stateL (:,:)        = 0.0_r8
          budg_stateG (:,:)        = 0.0_r8
       else
          call shr_sys_abort(subname//' ERROR in mode '//trim(mode))
       endif
    endif

  end subroutine WaterBudget_Reset

!-----------------------------------------------------------------------
  subroutine WaterBudget_Accum()
    !
    use clm_time_manager, only : get_curr_date, get_prev_date, get_nstep
    !
    implicit none
    !
    integer                :: ip, nf
    integer                :: year_prev, month_prev, day_prev, sec_prev
    integer                :: year_curr, month_curr, day_curr, sec_curr
    character(*),parameter :: subName = '(WaterBudget_Accum)'
    logical                :: update_state_beg, update_state_end

    call get_prev_date(year_prev, month_prev, day_prev, sec_prev)
    call get_curr_date(year_curr, month_curr, day_curr, sec_curr)

    do ip = p_inst+1, p_size
       budg_fluxL(:,ip) = budg_fluxL(:,ip) + budg_fluxL(:,p_inst)
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
          if (get_nstep() == 1) update_state_beg = .true.
          update_state_end = .true.
       end select

       if (update_state_beg) then
          nf = s_w_beg     ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
          nf = s_wcan_beg  ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
          nf = s_wsno_beg  ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
          nf = s_wsfc_beg  ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
          nf = s_wsliq_beg ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
          nf = s_wsice_beg ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
          nf = s_wwa_beg   ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
       endif

       if (update_state_end) then
          nf = s_w_end     ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
          nf = s_wcan_end  ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
          nf = s_wsno_end  ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
          nf = s_wsfc_end  ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
          nf = s_wsliq_end ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
          nf = s_wsice_end ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
          nf = s_wwa_end   ; budg_stateL(nf,ip) = budg_stateL(nf, p_inst)
       endif
       nf = s_w_errh2o  ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + budg_stateL(nf, p_inst)
    end do
    budg_fluxN(:,:) = budg_fluxN(:,:) + 1._r8
    
  end subroutine WaterBudget_Accum

  !-----------------------------------------------------------------------
  subroutine WaterBudget_Run(bounds, atm2lnd_vars, lnd2atm_vars, waterstate_vars, &
       soilhydrology_vars)
    !
    ! !DESCRIPTION:
    !
    use domainMod, only : ldomain
    use clm_varcon, only : re
    !
    implicit none

    type(bounds_type)        , intent(in) :: bounds
    type(atm2lnd_type)       , intent(in) :: atm2lnd_vars
    type(lnd2atm_type)       , intent(in) :: lnd2atm_vars
    type(waterstate_type)    , intent(in) :: waterstate_vars
    type(soilhydrology_type) , intent(in) :: soilhydrology_vars

    integer  :: g, nf, ip
    real(r8) :: af, one_over_re2

    associate(                                                             &
         forc_rain          => atm2lnd_vars%forc_rain_not_downscaled_grc , &
         forc_snow          => atm2lnd_vars%forc_snow_not_downscaled_grc , &
         qflx_evap_tot      => lnd2atm_vars%qflx_evap_tot_grc            , &
         qflx_rofice        => lnd2atm_vars%qflx_rofice_grc              , &
         qflx_rofliq_qsur   => lnd2atm_vars%qflx_rofliq_qsur_grc         , &
         qflx_rofliq_qsurp  => lnd2atm_vars%qflx_rofliq_qsurp_grc        , &
         qflx_rofliq_qsub   => lnd2atm_vars%qflx_rofliq_qsub_grc         , &
         qflx_rofliq_qsubp  => lnd2atm_vars%qflx_rofliq_qsubp_grc        , &
         qflx_rofliq_qgwl   => lnd2atm_vars%qflx_rofliq_qgwl_grc         , &
         begwb_grc          => waterstate_vars%begwb_grc                 , &
         endwb_grc          => waterstate_vars%endwb_grc                 , &
         beg_wa_grc         => soilhydrology_vars%beg_wa_grc             , &
         beg_h2ocan_grc     => waterstate_vars%beg_h2ocan_grc            , &
         beg_h2osno_grc     => waterstate_vars%beg_h2osno_grc            , &
         beg_h2osfc_grc     => waterstate_vars%beg_h2osfc_grc            , &
         beg_h2osoi_liq_grc => waterstate_vars%beg_h2osoi_liq_grc        , &
         beg_h2osoi_ice_grc => waterstate_vars%beg_h2osoi_ice_grc        , &
         end_wa_grc         => soilhydrology_vars%end_wa_grc             , &
         end_h2ocan_grc     => waterstate_vars%end_h2ocan_grc            , &
         end_h2osno_grc     => waterstate_vars%end_h2osno_grc            , &
         end_h2osfc_grc     => waterstate_vars%end_h2osfc_grc            , &
         end_h2osoi_liq_grc => waterstate_vars%end_h2osoi_liq_grc        , &
         end_h2osoi_ice_grc => waterstate_vars%end_h2osoi_ice_grc        , &
         errh2o_grc         => waterstate_vars%errh2o_grc                  &
         )

      ip = p_inst

      budg_stateL(:,ip) = 0._r8
      one_over_re2 = 1._r8/(re**2._r8)

      do g = bounds%begg, bounds%endg

         af   = (ldomain%area(g) * one_over_re2) * & ! area (converting km**2 to radians**2)
                ldomain%frac(g)                      ! land fraction

         nf = f_rain; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) + forc_rain(g)*af
         nf = f_snow; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) + forc_snow(g)*af
         nf = f_evap; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) - qflx_evap_tot(g)*af
         nf = f_roff; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) &
              - qflx_rofliq_qsur(g)*af - qflx_rofliq_qsurp(g)*af &
              - qflx_rofliq_qsub(g)*af - qflx_rofliq_qsubp(g)*af &
              - qflx_rofliq_qgwl(g)*af
         nf = f_ioff; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) - qflx_rofice(g)*af

         nf = s_w_beg     ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + begwb_grc(g)          *af
         nf = s_w_end     ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + endwb_grc(g)          *af
         nf = s_wcan_beg  ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + beg_h2ocan_grc(g)     *af
         nf = s_wcan_end  ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + end_h2ocan_grc(g)     *af
         nf = s_wsno_beg  ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + beg_h2osno_grc(g)     *af
         nf = s_wsno_end  ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + end_h2osno_grc(g)     *af
         nf = s_wsfc_beg  ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + beg_h2osfc_grc(g)     *af
         nf = s_wsfc_end  ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + end_h2osfc_grc(g)     *af
         nf = s_wsliq_beg ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + beg_h2osoi_liq_grc(g) *af
         nf = s_wsliq_end ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + end_h2osoi_liq_grc(g) *af
         nf = s_wsice_beg ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + beg_h2osoi_ice_grc(g) *af
         nf = s_wsice_end ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + end_h2osoi_ice_grc(g) *af
         nf = s_wwa_beg   ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + beg_wa_grc(g)         *af
         nf = s_wwa_end   ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + end_wa_grc(g)         *af
         nf = s_w_errh2o  ; budg_stateL(nf,ip) = budg_stateL(nf,ip) + errh2o_grc(g)         *af
      end do

    end associate

  end subroutine WaterBudget_Run

!-----------------------------------------------------------------------
  subroutine WaterBudget_Sum0()
    !
    use spmdMod    , only : mpicom
    use shr_mpi_mod, only : shr_mpi_sum
    !
    implicit none
    !
    real(r8)               :: budg_fluxGtmp(f_size,p_size) ! temporary sum
    real(r8)               :: budg_stateGtmp(s_size,p_size) ! temporary sum
    character(*),parameter :: subName = '(WaterBudget_Sum0)'

    budg_fluxGtmp = 0._r8
    budg_stateGtmp = 0._r8

    call shr_mpi_sum(budg_fluxL, budg_fluxGtmp, mpicom, subName)
    call shr_mpi_sum(budg_stateL, budg_stateGtmp, mpicom, subName)

    budg_fluxG  = budg_fluxG + budg_fluxGtmp
    budg_stateG = budg_stateGtmp

    budg_fluxL            = 0._r8 ! reset all fluxes
    budg_stateL(:,p_inst) = 0._r8 ! only reset instantaneous states
    
  end subroutine WaterBudget_Sum0

  !-----------------------------------------------------------------------

  subroutine WaterBudget_Print(budg_print_inst,  budg_print_daily,  budg_print_month,  &
       budg_print_ann,  budg_print_ltann,  budg_print_ltend)
    !
    use clm_time_manager, only : get_curr_date, get_prev_date, get_nstep, get_step_size
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
    integer :: s,f,ic,nf,ip,is ! data array indicies
    integer :: plev        ! print level
    integer :: year, mon, day, sec
    integer :: cdate
    logical :: sumdone
    real(r8) :: unit_conversion
    real(r8) :: budg_fluxGpr (f_size,p_size) ! values to print, scaled and such

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
             call WaterBudget_Sum0()
             budg_fluxGpr = budg_fluxG
             budg_fluxGpr = budg_fluxGpr*unit_conversion
             budg_fluxGpr = budg_fluxGpr/budg_fluxN
          end if

          if (ip == p_day .and. get_nstep() == 1) cycle
          if (ip == p_mon .and. get_nstep() == 1) cycle
          if (ip == p_ann .and. get_nstep() == 1) cycle
          if (ip == p_inf .and. get_nstep() == 1) cycle

          if (masterproc) then
             write(iulog,*)''
             write(iulog,*)'NET WATER FLUXES : period ',trim(pname(ip)),': date = ',cdate,sec
             write(iulog,FA0)'  Time  ','  Time    '
             write(iulog,FA0)'averaged','integrated'
             write(iulog,FA0)'kg/m2s*1e6','kg/m2*1e6'
             write(iulog,'(32("-"),"|",20("-"))')
             do f = 1, f_size
                write(iulog,FF)fname(f),budg_fluxGpr(f,ip),budg_fluxG(f,ip)*unit_conversion*get_step_size()
             end do
             write(iulog,'(32("-"),"|",20("-"))')
             write(iulog,FF)'   *SUM*', &
                  sum(budg_fluxGpr(:,ip)), sum(budg_fluxG(:,ip))*unit_conversion*get_step_size()
             write(iulog,'(32("-"),"|",20("-"))')

             write(iulog,*)''
             write(iulog,*)'WATER STATES (kg/m2*1e6): period ',trim(pname(ip)),': date = ',cdate,sec
             write(iulog,FS0),&
                  '       Canopy  ', &
                  '       Snow    ', &
                  '       SFC     ', &
                  '     Soil Liq  ', &
                  '     Soil Ice  ', &
                  '      Aquifer  ', &
                  ' Grid-level Err', &
                  '       TOTAL   '
             write(iulog,'(143("-"),"|",23("-"))')
             write(iulog,FS), '         beg', &
                  budg_stateG(s_wcan_beg  , ip)*unit_conversion, &
                  budg_stateG(s_wsno_beg  , ip)*unit_conversion, &
                  budg_stateG(s_wsfc_beg  , ip)*unit_conversion, &
                  budg_stateG(s_wsliq_beg , ip)*unit_conversion, &
                  budg_stateG(s_wsice_beg , ip)*unit_conversion, &
                  budg_stateG(s_wwa_beg   , ip)*unit_conversion, &
                  budg_stateG(s_w_beg     , ip)*unit_conversion
             write(iulog,FS), '         end', &
                  budg_stateG(s_wcan_end  , ip)*unit_conversion, &
                  budg_stateG(s_wsno_end  , ip)*unit_conversion, &
                  budg_stateG(s_wsfc_end  , ip)*unit_conversion, &
                  budg_stateG(s_wsliq_end , ip)*unit_conversion, &
                  budg_stateG(s_wsice_end , ip)*unit_conversion, &
                  budg_stateG(s_wwa_end   , ip)*unit_conversion, &
                  budg_stateG(s_w_end     , ip)*unit_conversion
             write(iulog,FS3)'*NET CHANGE*', &
                  (budg_stateG(s_wcan_end  ,ip) - budg_stateG(s_wcan_beg  ,ip))*unit_conversion, &
                  (budg_stateG(s_wsno_end  ,ip) - budg_stateG(s_wsno_beg  ,ip))*unit_conversion, &
                  (budg_stateG(s_wsfc_end  ,ip) - budg_stateG(s_wsfc_beg  ,ip))*unit_conversion, &
                  (budg_stateG(s_wsliq_end ,ip) - budg_stateG(s_wsliq_beg ,ip))*unit_conversion, &
                  (budg_stateG(s_wsice_end ,ip) - budg_stateG(s_wsice_beg ,ip))*unit_conversion, &
                  (budg_stateG(s_wwa_end   ,ip) - budg_stateG(s_wwa_beg   ,ip))*unit_conversion, &
                  budg_stateG(s_w_errh2o ,ip) *unit_conversion, &
                  (budg_stateG(s_w_end     ,ip) - budg_stateG(s_w_beg     ,ip))*unit_conversion
             write(iulog,'(143("-"),"|",23("-"))')
             write(iulog,FS2)'   *SUM*    ', &
                  (budg_stateG(s_wcan_end  ,ip) - budg_stateG(s_wcan_beg  ,ip))*unit_conversion + &
                  (budg_stateG(s_wsno_end  ,ip) - budg_stateG(s_wsno_beg  ,ip))*unit_conversion + &
                  (budg_stateG(s_wsfc_end  ,ip) - budg_stateG(s_wsfc_beg  ,ip))*unit_conversion + &
                  (budg_stateG(s_wsliq_end ,ip) - budg_stateG(s_wsliq_beg ,ip))*unit_conversion + &
                  (budg_stateG(s_wsice_end ,ip) - budg_stateG(s_wsice_beg ,ip))*unit_conversion + &
                  (budg_stateG(s_wwa_end   ,ip) - budg_stateG(s_wwa_beg   ,ip))*unit_conversion + &
                  -budg_stateG(s_w_errh2o ,ip) *unit_conversion, &
                  (budg_stateG(s_w_end     ,ip) - budg_stateG(s_w_beg     ,ip))*unit_conversion
             write(iulog,'(143("-"),"|",23("-"))')
          end if
       end if
    end do

  end subroutine WaterBudget_Print

  !-----------------------------------------------------------------------
  subroutine WaterBudget_Restart(bounds, ncid, flag)
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
    character(len=*),parameter :: subname = 'WaterBudget_Restart'

    select case (trim(flag))
    case ('define')
       call WaterBudget_Restart_Define(bounds, ncid)
    case ('write')
       call WaterBudget_Restart_Write(bounds, ncid, flag)
    case ('read')
       call WaterBudget_Restart_Read(bounds, ncid, flag)
    case default
       write(iulog,*) trim(subname),' ERROR: unknown flag = ',flag
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine WaterBudget_Restart

  !-----------------------------------------------------------------------
  subroutine WaterBudget_Restart_Define(bounds, ncid)
    !
    use ncdio_pio, only : file_desc_t, ncd_io, ncd_double, ncd_defvar
    !
    implicit none
    !
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id

    call ncd_defvar(varname='budg_fluxG', xtype=ncd_double, &
         dim1name='budg_flux', &
         long_name='budg_fluxG', units='mm', ncid=ncid)

    call ncd_defvar(varname='budg_fluxN', xtype=ncd_double, &
         dim1name='budg_flux', &
         long_name='budg_fluxN', units='-', ncid=ncid)

    call ncd_defvar(varname='budg_stateG', xtype=ncd_double, &
         dim1name='budg_state', &
         long_name='budg_stateG', units='mm', ncid=ncid)

  end subroutine WaterBudget_Restart_Define

  !-----------------------------------------------------------------------
  subroutine WaterBudget_Restart_Write(bounds, ncid, flag)
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
    real(r8) :: budg_fluxG_1D (f_size*p_size)
    real(r8) :: budg_fluxN_1D (f_size*p_size)
    real(r8) :: budg_stateG_1D(s_size*p_size)
    integer  :: f, s, p, count
    character(*),parameter :: subName = '(WaterBudget_Restart_Write) '

    call shr_mpi_sum(budg_fluxL, budg_fluxGtmp, mpicom, subName)
    call shr_mpi_sum(budg_stateL, budg_stateGtmp, mpicom, subName)

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

    call ncd_io(flag=flag, varname='budg_fluxG', data=budg_fluxG_1D, ncid=ncid)
    call ncd_io(flag=flag, varname='budg_fluxN', data=budg_fluxN_1D, ncid=ncid)
    call ncd_io(flag=flag, varname='budg_stateG', data=budg_stateG_1D, ncid=ncid)

  end subroutine WaterBudget_Restart_Write

  !-----------------------------------------------------------------------
  subroutine WaterBudget_Restart_Read(bounds, ncid, flag)
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
    integer  :: f, s, p, count

    call ncd_io(flag=flag, varname='budg_fluxG', data=budg_fluxG_1D, ncid=ncid)
    call ncd_io(flag=flag, varname='budg_fluxN', data=budg_fluxN_1D, ncid=ncid)
    call ncd_io(flag=flag, varname='budg_stateG', data=budg_stateG_1D, ncid=ncid)

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
    if (masterproc) then
       count = 0
       do s = 1, s_size
          do p = 1, p_size
             count = count + 1
             budg_stateL(s,p) = budg_stateG_1D(count)
          end do
       end do
    end if

  end subroutine WaterBudget_Restart_Read

   !-----------------------------------------------------------------------
  subroutine WaterBudget_SetBeginningMonthlyStates(bounds, waterstate_vars)
    !
    ! !DESCRIPTION:
    ! Set grid-level water states at the beginning of a month
    !
    ! !USES:
    use subgridAveMod    , only : p2c, c2g
    use clm_varpar       , only : nlevgrnd, nlevsoi, nlevurb
    use clm_varcon       , only : spval
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall 
    use column_varcon    , only : icol_road_perv, icol_road_imperv
    use clm_time_manager , only : get_curr_date, get_prev_date, get_nstep
    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds
    type(waterstate_type)     , intent(inout) :: waterstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: year_prev, month_prev, day_prev, sec_prev
    integer :: year_curr, month_curr, day_curr, sec_curr
    !-----------------------------------------------------------------------

    associate(                                                       & 
         begwb             =>    waterstate_vars%begwb_col         , & ! Output: [real(r8) (:)   ]  water mass begining of the time step
         endwb             =>    waterstate_vars%endwb_col         , & ! Output: [real(r8) (:)   ]  water mass begining of the time step
         tws_month_beg_grc =>    waterstate_vars%tws_month_beg_grc   & ! Output: [real(r8) (:)   ]  grid-level water mass at the begining of a month
         )

      ! Get current and previous dates to determine if a new month started
      call get_prev_date(year_curr, month_curr, day_curr, sec_curr);
      call get_prev_date(year_prev, month_prev, day_prev, sec_prev)

      ! If at the beginning of a simulation, save grid-level TWS based on
      ! 'begwb' from the current time step
      if ( day_curr == 1 .and. sec_curr == 0 .and. get_nstep() <= 1 ) then
         call c2g( bounds, &
              begwb(bounds%begc:bounds%endc), &
              tws_month_beg_grc(bounds%begg:bounds%endg), &
              c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
      endif

      ! If multiple steps into a simulation and the last time step was the
      ! end of a month, save grid-level TWS based on 'endwb' from the last
      ! time step
      if (get_nstep() > 1 .and. day_prev == 1 .and. sec_prev == 0) then
         call c2g( bounds, &
              endwb(bounds%begc:bounds%endc), &
              tws_month_beg_grc(bounds%begg:bounds%endg), &
              c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
      endif

    end associate

  end subroutine WaterBudget_SetBeginningMonthlyStates

   !-----------------------------------------------------------------------
  subroutine WaterBudget_SetEndingMonthlyStates(bounds, waterstate_vars)
    !
    ! !DESCRIPTION:
    ! Set grid-level water states at the end of a month
    !
    ! !USES:
    use subgridAveMod    , only : c2g
    use clm_varcon       , only : spval
    use clm_time_manager , only : get_curr_date, get_nstep
    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds
    type(waterstate_type)     , intent(inout) :: waterstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: year, mon, day, sec
    !-----------------------------------------------------------------------

    associate(                                                       & 
         endwb             =>    waterstate_vars%endwb_col         , & ! Output: [real(r8) (:)   ]  water mass at end of the time step
         tws_month_end_grc =>    waterstate_vars%tws_month_end_grc   & ! Output: [real(r8) (:)   ]  grid-level water mass at the end of a month
         )

      ! If this is the end of a month, save grid-level total water storage
      call get_curr_date(year, mon, day, sec);

      if (get_nstep() >= 1 .and. (day == 1 .and. sec == 0)) then
         call c2g( bounds, &
              endwb(bounds%begc:bounds%endc), &
              tws_month_end_grc(bounds%begg:bounds%endg), &
              c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
      else
         tws_month_end_grc(bounds%begg:bounds%endg) = spval
      end if

    end associate

  end subroutine WaterBudget_SetEndingMonthlyStates

end module WaterBudgetMod
