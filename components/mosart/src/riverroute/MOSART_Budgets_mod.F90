module MOSART_Budgets_mod
! Description: MOSART global water budget disgnostics
! 
! Added by Tian Zhou for E3SM v3, 11/06/2023. 
!-----------------------------------------------------------------------
  ! USES:
  use rof_cpl_indices, only : nt_rtm, nt_nliq, nt_nice
  use RtmVar         , only : iulog
  use RtmSpmd        , only : masterproc
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_sys_mod    , only : shr_sys_abort
  use shr_const_mod  , only : SHR_CONST_REARTH ! earth radius in m
  use shr_const_mod  , only : shr_const_pi

  implicit none
  private

  public MOSART_WaterBudget_Reset
  public MOSART_WaterBudget_Extraction
  public MOSART_WaterBudget_Print
  public MOSART_WaterBudget_Restart_Write
  public MOSART_WaterBudget_Restart_Read
  
    !--- F for flux ---

  integer, parameter :: f_roff = 1
  integer, parameter :: f_rout = 2
  integer, parameter :: f_irri = 3

  integer, parameter, public :: f_size = f_irri

  character(len=16),parameter :: fname(f_size) = &
       (/&
       '     runoff (in)', &
       'streamflow (out)', &
       '           irrig'  &
       /)

   !--- O for "other" term aggregated from flux to state ---

  integer, parameter :: o_othr = 1

  integer, parameter, public :: o_size = o_othr

  !--- S for state ---
  integer, parameter :: s_w_beg         = 1
  integer, parameter :: s_w_end         = 2
  integer, parameter :: s_wchannel_beg  = 3
  integer, parameter :: s_wchannel_end  = 4
  integer, parameter :: s_wsubnet_beg   = 5
  integer, parameter :: s_wsubnet_end   = 6
  integer, parameter :: s_whillslop_beg = 7
  integer, parameter :: s_whillslop_end = 8
  integer, parameter :: s_wres_beg      = 9
  integer, parameter :: s_wres_end      = 10
  integer, parameter :: s_wflood_beg    = 11
  integer, parameter :: s_wflood_end    = 12
  integer, parameter :: s_wother_beg    = 13
  integer, parameter :: s_wother_end    = 14

  integer, parameter, public :: s_size = s_wother_end

  !--- P for period ---

  integer, parameter :: p_inst = 1
  integer, parameter :: p_day  = 2
  integer, parameter :: p_mon  = 3
  integer, parameter :: p_ann  = 4
  integer, parameter :: p_inf  = 5

  integer, parameter, public :: p_size = p_inf

  character(len=8),parameter :: pname(p_size) = &
       (/'    inst','   daily',' monthly','  annual','all_time' /)

  real(r8), public :: rof_budg_fluxG  (f_size, p_size) = 0.0_r8 ! global flux sum (volume/time)
  real(r8), public :: rof_budg_other  (o_size, p_size) = 0.0_r8 ! global "other" term (volume)
  real(r8), public :: rof_budg_stateG (s_size, p_size) = 0.0_r8 ! global state sum (volume)
  real(r8), public :: rof_budg_fluxN  (f_size, p_size) = 0.0_r8 ! counter, valid only on root pe

  real(r8) :: rof_budg_other_rest (o_size,p_size) = 0.0_r8 ! "other" term from restart file

  real(r8), target, public :: rof_budg_fluxG_1D (f_size*p_size) = 0.0_r8
  real(r8), target, public :: rof_budg_stateG_1D(s_size*p_size) = 0.0_r8
  real(r8), target, public :: rof_budg_fluxN_1D (f_size*p_size) = 0.0_r8

  !----- formats -----
  character(*),parameter :: FA0= "('    ',12x,(3x,a10,2x),' | ',(3x,a10,2x))"
  character(*),parameter :: FF = "('',a16,f15.8,' | ',f18.2)"
  character(*),parameter :: FS0= "(' ',12x,6(a16),' | ',(a16))"
  character(*),parameter :: FS2= "(' ',a12,40x,f16.2,40x,' | ',f16.2)"
  character(*),parameter :: FS3= "(' ',a12,6(f16.2),' | ',(f16.2))"

  !----- other variables -----
  real(r8)               :: unit_conversion

contains
  subroutine MOSART_WaterBudget_Reset(mode)
   !
   use RtmTimeManager, only : get_curr_date, get_prev_date
   !
   implicit none
   !
   character(len=*), intent(in),optional :: mode
   !
   integer :: year, mon, day, sec
   integer :: ip
   character(*),parameter :: subName = '(MOSART_WaterBudget_Reset) '

   if (.not.present(mode)) then
      call get_curr_date(year, mon, day, sec)

      do ip = 1,p_size
         if (ip == p_inst) then
            rof_budg_fluxG(:,ip)  = 0.0_r8
            rof_budg_fluxN(:,ip)  = 0.0_r8
         endif
         if (ip==p_day .and. sec==0) then
            rof_budg_fluxG(:,ip)  = 0.0_r8
            rof_budg_fluxN(:,ip)  = 0.0_r8
         endif
         if (ip==p_mon .and. day==1 .and. sec==0) then
            rof_budg_fluxG(:,ip)  = 0.0_r8
            rof_budg_fluxN(:,ip)  = 0.0_r8
         endif
         if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
            rof_budg_fluxG(:,ip)  = 0.0_r8
            rof_budg_fluxN(:,ip)  = 0.0_r8
         endif
      enddo

   else

      if (trim(mode) == 'inst') then
         rof_budg_fluxG  (:,p_inst)   = 0.0_r8
         rof_budg_other  (:,p_inst)   = 0.0_r8
         rof_budg_stateG (:,p_inst)   = 0.0_r8
         rof_budg_fluxN  (:,p_inst)   = 0.0_r8
      elseif (trim(mode) == 'day') then
         rof_budg_fluxG  (:,p_day)    = 0.0_r8
         rof_budg_other  (:,p_day)    = 0.0_r8
         rof_budg_stateG (:,p_day)    = 0.0_r8
         rof_budg_fluxN  (:,p_day)    = 0.0_r8
      elseif (trim(mode) == 'mon') then
         rof_budg_fluxG  (:,p_mon)    = 0.0_r8
         rof_budg_other  (:,p_mon)    = 0.0_r8
         rof_budg_stateG (:,p_mon)    = 0.0_r8
         rof_budg_fluxN  (:,p_mon)    = 0.0_r8
      elseif (trim(mode) == 'ann') then
         rof_budg_fluxG  (:,p_ann)    = 0.0_r8
         rof_budg_other  (:,p_ann)    = 0.0_r8
         rof_budg_stateG (:,p_ann)    = 0.0_r8
         rof_budg_fluxN  (:,p_ann)    = 0.0_r8
      elseif (trim(mode) == 'inf') then
         rof_budg_fluxG  (:,p_inf)    = 0.0_r8
         rof_budg_other  (:,p_inf)    = 0.0_r8
         rof_budg_stateG (:,p_inf)    = 0.0_r8
         rof_budg_fluxN  (:,p_inf)    = 0.0_r8
      elseif (trim(mode) == 'all') then
         rof_budg_fluxG  (:,:)        = 0.0_r8
         rof_budg_other  (:,:)        = 0.0_r8
         rof_budg_stateG (:,:)        = 0.0_r8
         rof_budg_fluxN  (:,:)        = 0.0_r8
      else
         call shr_sys_abort(subname//' ERROR in mode '//trim(mode))
      endif
   endif

  end subroutine MOSART_WaterBudget_Reset
!--------------------------------------------------------------------
  subroutine MOSART_WaterBudget_Extraction(budget_global, budget_terms_total, bv_volt_beg, bv_volt_end, &
   bv_wt_beg, bv_wt_end, bv_wr_beg, bv_wr_end, bv_wh_beg, bv_wh_end, bv_dstor_beg, bv_dstor_end, bv_fp_beg, bv_fp_end, br_sup, &
   budget_in, budget_out, budget_oth)

   use RtmTimeManager, only : get_curr_date, get_prev_date, get_nstep, get_step_size
   implicit none
   !
   integer,  intent(in) :: budget_terms_total
   real(r8), intent(in) :: budget_global(budget_terms_total,nt_rtm)    ! global budget sums
   integer,  intent(in) :: bv_volt_beg  ! = 1  ! initial total volume
   integer,  intent(in) :: bv_volt_end  ! = 2  ! final   total volume
   integer,  intent(in) :: bv_wt_beg    ! = 3  ! initial wt volume
   integer,  intent(in) :: bv_wt_end    ! = 4  ! final   wt volume
   integer,  intent(in) :: bv_wr_beg    ! = 5  ! initial wr volume
   integer,  intent(in) :: bv_wr_end    ! = 6  ! final   wr volume
   integer,  intent(in) :: bv_wh_beg    ! = 7  ! initial wh volume
   integer,  intent(in) :: bv_wh_end    ! = 8  ! final   wh volume
   integer,  intent(in) :: bv_dstor_beg ! = 9  ! initial reservoir storage
   integer,  intent(in) :: bv_dstor_end ! = 10 ! final   reservoir storage
   integer,  intent(in) :: bv_fp_beg    ! = 11 ! Initial water volume over floodplains.
   integer,  intent(in) :: bv_fp_end    ! = 12 ! Final water volume over floodplains.
   integer,  intent(in) :: br_sup       ! = 49 ! supply rate (m3/coupling period)

   real(r8), intent(in) :: budget_in, budget_out, budget_oth

   integer                :: ip, nf
   integer                :: year_prev, month_prev, day_prev, sec_prev
   integer                :: year_curr, month_curr, day_curr, sec_curr
   logical                :: update_state_beg, update_state_end

   character(*),parameter :: subName = '(Extraction) '

   unit_conversion = 1.d0/(4.0_r8*shr_const_pi*SHR_CONST_REARTH**2)*1.0e15_r8

   rof_budg_fluxG(f_roff,p_inst) = budget_in / get_step_size() 
   rof_budg_fluxG(f_rout,p_inst) = -(budget_out - budget_global(br_sup, nt_nliq)) / get_step_size() !subtract the supply part as it has been aggregated to total out
   rof_budg_fluxG(f_irri,p_inst) = -budget_global(br_sup, nt_nliq) / get_step_size() 
   rof_budg_other(o_othr,p_inst) = budget_oth / get_step_size() 

   rof_budg_stateG(s_w_beg, p_inst) =  budget_global(bv_volt_beg,nt_nliq) * unit_conversion ! assign values
   rof_budg_stateG(s_w_end, p_inst) =  budget_global(bv_volt_end,nt_nliq) * unit_conversion ! 

   rof_budg_stateG(s_wchannel_beg, p_inst)  =  budget_global(bv_wr_beg,nt_nliq) * unit_conversion 
   rof_budg_stateG(s_wchannel_end, p_inst)  =  budget_global(bv_wr_end,nt_nliq) * unit_conversion
   rof_budg_stateG(s_wsubnet_beg, p_inst)   =  budget_global(bv_wt_beg,nt_nliq) * unit_conversion 
   rof_budg_stateG(s_wsubnet_end, p_inst)   =  budget_global(bv_wt_end,nt_nliq) * unit_conversion 
   rof_budg_stateG(s_whillslop_beg, p_inst) =  budget_global(bv_wh_beg,nt_nliq) * unit_conversion 
   rof_budg_stateG(s_whillslop_end, p_inst) =  budget_global(bv_wh_end,nt_nliq) * unit_conversion  
   
   rof_budg_stateG(s_wres_beg, p_inst) =  budget_global(bv_dstor_beg,nt_nliq) * unit_conversion
   rof_budg_stateG(s_wres_end, p_inst) =  budget_global(bv_dstor_end,nt_nliq) * unit_conversion
   rof_budg_stateG(s_wflood_beg, p_inst) =  budget_global(bv_fp_beg,nt_nliq) * unit_conversion
   rof_budg_stateG(s_wflood_end, p_inst) =  budget_global(bv_fp_end,nt_nliq) * unit_conversion

   call get_prev_date(year_prev, month_prev, day_prev, sec_prev)
   call get_curr_date(year_curr, month_curr, day_curr, sec_curr)

   do ip = p_inst+1, p_size
      rof_budg_fluxG(:,ip) = rof_budg_fluxG(:,ip) + rof_budg_fluxG(:,p_inst)
      rof_budg_other(:,ip) = rof_budg_other(:,ip) + rof_budg_other(:,p_inst)

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
         nf = s_w_beg        ; rof_budg_stateG(nf,ip) = rof_budg_stateG(nf, p_inst)
         nf = s_wchannel_beg ; rof_budg_stateG(nf,ip) = rof_budg_stateG(nf, p_inst)
         nf = s_wsubnet_beg  ; rof_budg_stateG(nf,ip) = rof_budg_stateG(nf, p_inst)
         nf = s_whillslop_beg; rof_budg_stateG(nf,ip) = rof_budg_stateG(nf, p_inst)
         nf = s_wres_beg     ; rof_budg_stateG(nf,ip) = rof_budg_stateG(nf, p_inst)
         nf = s_wflood_beg   ; rof_budg_stateG(nf,ip) = rof_budg_stateG(nf, p_inst)
         nf = s_wother_beg   ; rof_budg_stateG(nf,ip) = rof_budg_other_rest(o_othr,ip) + (rof_budg_other(o_othr, ip) - rof_budg_other(o_othr,p_inst)) * unit_conversion * get_step_size()
      endif

      if (update_state_end) then
         nf = s_w_end        ; rof_budg_stateG(nf,ip) = rof_budg_stateG(nf, p_inst)
         nf = s_wchannel_end ; rof_budg_stateG(nf,ip) = rof_budg_stateG(nf, p_inst)
         nf = s_wsubnet_end  ; rof_budg_stateG(nf,ip) = rof_budg_stateG(nf, p_inst)
         nf = s_whillslop_end; rof_budg_stateG(nf,ip) = rof_budg_stateG(nf, p_inst)
         nf = s_wres_end     ; rof_budg_stateG(nf,ip) = rof_budg_stateG(nf, p_inst)
         nf = s_wflood_end   ; rof_budg_stateG(nf,ip) = rof_budg_stateG(nf, p_inst)
         nf = s_wother_end   ; rof_budg_stateG(nf,ip) = rof_budg_other_rest(o_othr, ip) + rof_budg_other(o_othr, ip) * unit_conversion * get_step_size() 
      endif

   end do

   rof_budg_fluxN(:,:) = rof_budg_fluxN(:,:) + 1._r8

  end subroutine MOSART_WaterBudget_Extraction
!-----------------------------------------------------------------------
  subroutine MOSART_WaterBudget_Print()
   !
   use RtmTimeManager, only : get_curr_date, get_prev_date, get_nstep, get_step_size
   !
   implicit none
   !
   integer :: budg_print_inst  = 0
   integer :: budg_print_daily = 0
   integer :: budg_print_month = 1
   integer :: budg_print_ann   = 1
   integer :: budg_print_ltann = 1
   integer :: budg_print_ltend = 0
   !
   ! !LOCAL VARIABLES:
   integer :: s,f,ic,nf,ip,is ! data array indicies
   integer :: plev        ! print level
   integer :: year, mon, day, sec
   integer :: cdate
   logical :: sumdone
   real(r8) :: budg_fluxGpr (f_size,p_size) ! values to print, scaled and such

   character(*),parameter :: subName = '(MOSART_WaterBudget_Print) '

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
            unit_conversion = 1.d0/(4.0_r8*shr_const_pi*SHR_CONST_REARTH**2)*1.0e15_r8
         if (.not.sumdone) then
            sumdone = .true.
            budg_fluxGpr = rof_budg_fluxG
            budg_fluxGpr = budg_fluxGpr*unit_conversion
            budg_fluxGpr = budg_fluxGpr/rof_budg_fluxN
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
               write(iulog,FF)fname(f),budg_fluxGpr(f,ip),rof_budg_fluxG(f,ip)*unit_conversion*get_step_size()
            end do
            write(iulog,'(32("-"),"|",20("-"))')
            write(iulog,FF)'   *SUM*', &
                  sum(budg_fluxGpr(:,ip)), sum(rof_budg_fluxG(:,ip))*unit_conversion*get_step_size()
            write(iulog,'(32("-"),"|",20("-"))')
            write(iulog,*)''

            write(iulog,*)'WATER STATES (kg/m2*1e6): period ',trim(pname(ip)),': date = ',cdate,sec
            write(iulog,FS0) &
                  '      channel ', &
                  '   subnetwork ', &
                  '    hillslope ', &
                  '    reservoir ', &
                  '   floodplain ', &
                  '        other ', &
                  '        TOTAL '
            write(iulog,'(110("-"),"|",20("-"))')
            write(iulog,FS3) '         beg', &
                  rof_budg_stateG(s_wchannel_beg, ip), &
                  rof_budg_stateG(s_wsubnet_beg, ip), &
                  rof_budg_stateG(s_whillslop_beg, ip), &
                  rof_budg_stateG(s_wres_beg    , ip), &
                  rof_budg_stateG(s_wflood_beg  , ip), &
               - rof_budg_stateG(s_wother_beg  , ip), & ! print out negative value to reflect the actural calculation
                  rof_budg_stateG(s_w_beg       , ip)- &
                  rof_budg_stateG(s_wother_beg  , ip)+ &
                  rof_budg_stateG(s_wres_beg    , ip) ! s_w_beg doesn't include reservoir storage
            write(iulog,FS3) '         end', &
                  rof_budg_stateG(s_wchannel_end, ip), &
                  rof_budg_stateG(s_wsubnet_end, ip), &
                  rof_budg_stateG(s_whillslop_end, ip), &
                  rof_budg_stateG(s_wres_end    , ip), &
                  rof_budg_stateG(s_wflood_end  , ip), &
               - rof_budg_stateG(s_wother_end  , ip), & ! print out negative value in the table for display only, actual calucation uses positive
                  rof_budg_stateG(s_w_end       , ip)- &
                  rof_budg_stateG(s_wother_end  , ip)+ &
                  rof_budg_stateG(s_wres_end    , ip) ! s_w_end doesn't include reservoir storage
            write(iulog,FS3)'*NET CHANGE*', &
                  (rof_budg_stateG(s_wchannel_end,ip) - rof_budg_stateG(s_wchannel_beg,ip)), &
                  (rof_budg_stateG(s_wsubnet_end,ip) - rof_budg_stateG(s_wsubnet_beg,ip)), &
                  (rof_budg_stateG(s_whillslop_end,ip) - rof_budg_stateG(s_whillslop_beg,ip)), &
                  (rof_budg_stateG(s_wres_end    ,ip) - rof_budg_stateG(s_wres_beg    ,ip)), &
                  (rof_budg_stateG(s_wflood_end  ,ip) - rof_budg_stateG(s_wflood_beg  ,ip)), &
                - (rof_budg_stateG(s_wother_end  ,ip) - rof_budg_stateG(s_wother_beg  ,ip)), &
                  (rof_budg_stateG(s_w_end       ,ip) - rof_budg_stateG(s_wother_end  , ip) + rof_budg_stateG(s_wres_end , ip)) - & 
                  (rof_budg_stateG(s_w_beg       ,ip) - rof_budg_stateG(s_wother_beg  , ip) + rof_budg_stateG(s_wres_beg , ip))
            write(iulog,'(110("-"),"|",20("-"))')
            write(iulog,FS2)'       *SUM*', &
                  (rof_budg_stateG(s_wchannel_end,ip) - rof_budg_stateG(s_wchannel_beg,ip)) + &
                  (rof_budg_stateG(s_wsubnet_end,ip) - rof_budg_stateG(s_wsubnet_beg,ip)) + &
                  (rof_budg_stateG(s_whillslop_end,ip) - rof_budg_stateG(s_whillslop_beg,ip)) + &
                  (rof_budg_stateG(s_wres_end  ,ip) - rof_budg_stateG(s_wres_beg  ,ip)) + &
                  (rof_budg_stateG(s_wflood_end  ,ip) - rof_budg_stateG(s_wflood_beg  ,ip)) - &
                  (rof_budg_stateG(s_wother_end  ,ip) - rof_budg_stateG(s_wother_beg  ,ip)), &
                  (rof_budg_stateG(s_w_end       ,ip) - rof_budg_stateG(s_wother_end  , ip) + rof_budg_stateG(s_wres_end , ip)) - & 
                  (rof_budg_stateG(s_w_beg       ,ip) - rof_budg_stateG(s_wother_beg  , ip) + rof_budg_stateG(s_wres_beg , ip))
            write(iulog,'(110("-"),"|",20("-"))')
         end if
      end if
   end do
  end subroutine MOSART_WaterBudget_Print
!-----------------------------------------------------------------------
  subroutine MOSART_WaterBudget_Restart_Write()
   !
   implicit none
   !
       ! !LOCAL VARIABLES:

   integer  :: f, s, o, p, count
   character(*),parameter :: subName = '(MOSART_WaterBudget_Restart_Write) '

       ! Copy data from 2D into 1D array
   count = 0
   do f = 1, f_size
      do p = 1, p_size
         count = count + 1
         rof_budg_fluxG_1D(count) = rof_budg_fluxG(f,p)
         rof_budg_fluxN_1D(count) = rof_budg_fluxN(f,p)
      end do
   end do

   ! Copy data from 2D into 1D array
   count = 0
   do s = 1, s_size
      do p = 1, p_size
         count = count + 1
         rof_budg_stateG_1D(count) = rof_budg_stateG(s,p)
      end do
   end do

  end subroutine MOSART_WaterBudget_Restart_Write
!-----------------------------------------------------------------------
  subroutine MOSART_WaterBudget_Restart_Read()
   !
   implicit none
   !
   integer  :: f, s, p, o, count
   ! Copy data from 1D into 2D array

   count = 0
   do f = 1, f_size
      do p = 1, p_size
         count = count + 1
         rof_budg_fluxG(f,p) = rof_budg_fluxG_1D(count)
         rof_budg_fluxN(f,p) = rof_budg_fluxN_1D(count)
      end do
   end do

   ! Copy data from 1D into 2D array
   if (masterproc) then
      count = 0
      do s = 1, s_size
         do p = 1, p_size
            count = count + 1
            rof_budg_stateG(s,p) = rof_budg_stateG_1D(count)
         end do
      end do
   end if

   ! save the "other term" from the restart file so that it can
   ! be used as initial values in the water state budget table
   rof_budg_other_rest(o_othr,:) = rof_budg_stateG(s_wother_end,:) 
   ! write(iulog,*) '(TZ)rof_budg_other_rest', rof_budg_other_rest
  end subroutine MOSART_WaterBudget_Restart_Read

end module MOSART_Budgets_mod
