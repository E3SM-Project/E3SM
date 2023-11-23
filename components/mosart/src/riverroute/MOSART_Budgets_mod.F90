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
  
    !--- F for flux ---

  integer, parameter :: f_roff = 1
  integer, parameter :: f_rout = 2
  integer, parameter :: f_irri = 3
  integer, parameter :: f_othr = 4

  integer, parameter, public :: f_size = f_othr

  character(len=16),parameter :: fname(f_size) = &
       (/&
       '     runoff (in)', &
       'streamflow (out)', &
       '           irrig', &
       '           other'  &
       /)

  !--- S for state ---
  integer, parameter :: s_w_beg         = 1
  integer, parameter :: s_w_end         = 2
  integer, parameter :: s_wchannel_beg  = 3
  integer, parameter :: s_wchannel_end  = 4
  integer, parameter :: s_wdam_beg      = 5
  integer, parameter :: s_wdam_end      = 6
  integer, parameter :: s_wflood_beg    = 7
  integer, parameter :: s_wflood_end    = 8
  integer, parameter :: s_wother_beg    = 9
  integer, parameter :: s_wother_end    = 10

  integer, parameter, public :: s_size = s_wother_end

  character(len=13),parameter :: sname(s_size) = &
       (/&
       '  total_w_beg', &
       '  total_w_end', &
       'channel_w_beg', &
       'channel_w_end', &
       '    dam_w_beg', &
       '    dam_w_end', &
       '  flood_w_beg', &
       '  flood_w_end', &
       '  other_w_beg', &
       '  other_w_end'  &
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

  real(r8) :: budg_fluxG(f_size, p_size) ! global flux sum (volume/time)
  real(r8), public :: budg_stateG(s_size, p_size) ! global state sum (volume)
  real(r8) :: budg_fluxN(f_size, p_size) ! counter, valid only on root pe

  !----- formats -----
  character(*),parameter :: FA0= "('    ',12x,(3x,a10,2x),' | ',(3x,a10,2x))"
  character(*),parameter :: FF = "('',a16,f15.8,' | ',f18.2)"
  character(*),parameter :: FS = "('    ',a12,3(f18.2),18x,' | ',(f18.2))"
  character(*),parameter :: FS0= "('    ',12x,4(a18),' | ',(a18))"
  character(*),parameter :: FS2= "('    ',a12,27x,f18.2,27x,' | ',f18.2)"
  character(*),parameter :: FS3= "('    ',a12,4(f18.2),' | ',(f18.2))"

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
            budg_fluxG(:,ip)  = 0.0_r8
            budg_fluxN(:,ip)  = 0.0_r8
         endif
         if (ip==p_day .and. sec==0) then
            budg_fluxG(:,ip)  = 0.0_r8
            budg_fluxN(:,ip)  = 0.0_r8
         endif
         if (ip==p_mon .and. day==1 .and. sec==0) then
            budg_fluxG(:,ip)  = 0.0_r8
            budg_fluxN(:,ip)  = 0.0_r8
         endif
         if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
            budg_fluxG(:,ip)  = 0.0_r8
            budg_fluxN(:,ip)  = 0.0_r8
         endif
      enddo

   else

      if (trim(mode) == 'inst') then
         budg_fluxG  (:,p_inst)   = 0.0_r8
         budg_stateG (:,p_inst)   = 0.0_r8
         budg_fluxN  (:,p_inst)   = 0.0_r8
      elseif (trim(mode) == 'day') then
         budg_fluxG  (:,p_day)    = 0.0_r8
         budg_stateG (:,p_day)    = 0.0_r8
         budg_fluxN  (:,p_day)    = 0.0_r8
      elseif (trim(mode) == 'mon') then
         budg_fluxG  (:,p_mon)    = 0.0_r8
         budg_stateG (:,p_mon)    = 0.0_r8
         budg_fluxN  (:,p_mon)    = 0.0_r8
      elseif (trim(mode) == 'ann') then
         budg_fluxG  (:,p_ann)    = 0.0_r8
         budg_stateG (:,p_ann)    = 0.0_r8
         budg_fluxN  (:,p_ann)    = 0.0_r8
      elseif (trim(mode) == 'inf') then
         budg_fluxG  (:,p_inf)    = 0.0_r8
         budg_stateG (:,p_inf)    = 0.0_r8
         budg_fluxN  (:,p_inf)    = 0.0_r8
      elseif (trim(mode) == 'all') then
         budg_fluxG  (:,:)        = 0.0_r8
         budg_stateG (:,:)        = 0.0_r8
         budg_fluxN  (:,:)        = 0.0_r8
      else
         call shr_sys_abort(subname//' ERROR in mode '//trim(mode))
      endif
   endif

  end subroutine MOSART_WaterBudget_Reset
!--------------------------------------------------------------------
  subroutine MOSART_WaterBudget_Extraction(budget_global, budget_terms_total, bv_volt_beg, bv_volt_end, &
   bv_wt_beg, bv_wt_end, bv_wr_beg, bv_wr_end, bv_wh_beg, bv_wh_end, bv_dstor_beg, bv_dstor_end, bv_fp_beg, bv_fp_end, &
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
   integer,  intent(in) :: bv_dstor_beg ! = 9  ! initial dam storage
   integer,  intent(in) :: bv_dstor_end ! = 10 ! final   dam storage
   integer,  intent(in) :: bv_fp_beg    ! = 11 ! Initial water volume over floodplains.
   integer,  intent(in) :: bv_fp_end    ! = 12 ! Final water volume over floodplains.

   real(r8), intent(in) :: budget_in, budget_out, budget_oth

   integer                :: ip, nf
   integer                :: year_prev, month_prev, day_prev, sec_prev
   integer                :: year_curr, month_curr, day_curr, sec_curr
   logical                :: update_state_beg, update_state_end

   real(r8)               :: budg_fluxGtmp(f_size,p_size) ! temporary sum
   real(r8)               :: budg_stateGtmp(s_size,p_size) ! temporary sum 

   character(*),parameter :: subName = '(Extraction) '

   unit_conversion = 1.d0/(4.0_r8*shr_const_pi*SHR_CONST_REARTH**2)*1.0e15_r8

   budg_fluxG(f_roff,p_inst) = budget_in / get_step_size() 
   budg_fluxG(f_rout,p_inst) = -budget_out / get_step_size() 
   budg_fluxG(f_irri,p_inst) = 0.0_r8 ! need update later 
   budg_fluxG(f_othr,p_inst) = budget_oth / get_step_size() 
     
   budg_stateG(s_w_beg, p_inst) =  budget_global(bv_volt_beg,nt_nliq) ! assign values
   budg_stateG(s_w_end, p_inst) =  budget_global(bv_volt_end,nt_nliq) ! 
   budg_stateG(s_wchannel_beg, p_inst) =  budget_global(bv_wt_beg,nt_nliq) + budget_global(bv_wr_beg,nt_nliq) + budget_global(bv_wh_beg,nt_nliq) 
   budg_stateG(s_wchannel_end, p_inst) =  budget_global(bv_wt_end,nt_nliq) + budget_global(bv_wr_end,nt_nliq) + budget_global(bv_wh_end,nt_nliq) 
   budg_stateG(s_wdam_beg, p_inst) =  budget_global(bv_dstor_beg,nt_nliq) 
   budg_stateG(s_wdam_end, p_inst) =  budget_global(bv_dstor_end,nt_nliq) 
   budg_stateG(s_wflood_beg, p_inst) =  budget_global(bv_fp_beg,nt_nliq)  
   budg_stateG(s_wflood_end, p_inst) =  budget_global(bv_fp_end,nt_nliq) 

   call get_prev_date(year_prev, month_prev, day_prev, sec_prev)
   call get_curr_date(year_curr, month_curr, day_curr, sec_curr)

   do ip = p_inst+1, p_size
      budg_fluxG(:,ip) = budg_fluxG(:,ip) + budg_fluxG(:,p_inst)

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
         nf = s_w_beg        ; budg_stateG(nf,ip) = budg_stateG(nf, p_inst)
         nf = s_wchannel_beg ; budg_stateG(nf,ip) = budg_stateG(nf, p_inst)
         nf = s_wdam_beg     ; budg_stateG(nf,ip) = budg_stateG(nf, p_inst)
         nf = s_wflood_beg   ; budg_stateG(nf,ip) = budg_stateG(nf, p_inst)
         nf = s_wother_beg   ; budg_stateG(nf,ip) = (budg_fluxG(f_othr, ip) - budg_fluxG(f_othr,p_inst))*get_step_size() 
      endif

      if (update_state_end) then
         nf = s_w_end        ; budg_stateG(nf,ip) = budg_stateG(nf, p_inst)
         nf = s_wchannel_end ; budg_stateG(nf,ip) = budg_stateG(nf, p_inst)
         nf = s_wdam_end     ; budg_stateG(nf,ip) = budg_stateG(nf, p_inst)
         nf = s_wflood_end   ; budg_stateG(nf,ip) = budg_stateG(nf, p_inst)
         nf = s_wother_end   ; budg_stateG(nf,ip) = budg_fluxG(f_othr, ip)*get_step_size() 
      endif

   end do
   budg_fluxN(:,:) = budg_fluxN(:,:) + 1._r8

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
          do f = 1, f_size-1 ! Do not include "other" term in the flux budget calcuation
             write(iulog,FF)fname(f),budg_fluxGpr(f,ip),budg_fluxG(f,ip)*unit_conversion*get_step_size()
          end do
          write(iulog,'(32("-"),"|",20("-"))')
          write(iulog,FF)'   *SUM*', &
                budg_fluxGpr(f_roff,ip)+budg_fluxGpr(f_rout,ip)+budg_fluxG(f_irri,ip), &
               (budg_fluxGpr(f_roff,ip)+budg_fluxGpr(f_rout,ip)+budg_fluxG(f_irri,ip))*unit_conversion*get_step_size()
          write(iulog,'(32("-"),"|",20("-"))')

         ! write(iulog,*) trim(subname),'(TZ)budg_fluxG', budg_fluxG(1,3),' fluxN=', budg_fluxN(1,3)

          write(iulog,*)''
          write(iulog,*)'WATER STATES (kg/m2*1e6): period ',trim(pname(ip)),': date = ',cdate,sec
          write(iulog,FS0) &
               '       channel   ', &
               '       reservoir ', &
               '       floodplain', &
               '       other     ', &
               '       TOTAL     '
          write(iulog,'(89("-"),"|",23("-"))')
          write(iulog,FS3) '         beg', &
               budg_stateG(s_wchannel_beg, ip), &
               budg_stateG(s_wdam_beg    , ip), &
               budg_stateG(s_wflood_beg  , ip), &
               budg_stateG(s_wother_beg  , ip), &
               budg_stateG(s_w_beg       , ip)
          write(iulog,FS3) '         end', &
               budg_stateG(s_wchannel_end, ip), &
               budg_stateG(s_wdam_end    , ip), &
               budg_stateG(s_wflood_end  , ip), &
               budg_stateG(s_wother_end  , ip), &
               budg_stateG(s_w_end       , ip)
          write(iulog,FS3)'*NET CHANGE*', &
               (budg_stateG(s_wchannel_end,ip) - budg_stateG(s_wchannel_beg,ip)), &
               (budg_stateG(s_wdam_end    ,ip) - budg_stateG(s_wdam_beg    ,ip)), &
               (budg_stateG(s_wflood_end  ,ip) - budg_stateG(s_wflood_beg  ,ip)), &
               (budg_stateG(s_wother_end  ,ip) - budg_stateG(s_wother_beg  ,ip)), &
               (budg_stateG(s_w_end       ,ip) - budg_stateG(s_w_beg       ,ip))
          write(iulog,'(89("-"),"|",23("-"))')
          write(iulog,FS2)'   *SUM*    ', &
               (budg_stateG(s_wchannel_end  ,ip) - budg_stateG(s_wchannel_beg  ,ip)) + &
               (budg_stateG(s_wdam_end  ,ip) - budg_stateG(s_wdam_beg  ,ip)) + &
               (budg_stateG(s_wflood_end  ,ip) - budg_stateG(s_wflood_beg  ,ip)) - &
               (budg_stateG(s_wother_end  ,ip) - budg_stateG(s_wother_beg  ,ip)), &
               !(budg_stateG(s_w_end     ,ip) - budg_stateG(s_w_beg     ,ip))
               ((budg_stateG(s_wchannel_end  ,ip) - budg_stateG(s_wchannel_beg  ,ip)) + &
               (budg_stateG(s_wdam_end  ,ip) - budg_stateG(s_wdam_beg  ,ip)) + &
               (budg_stateG(s_wflood_end  ,ip) - budg_stateG(s_wflood_beg  ,ip)) - &
               (budg_stateG(s_wother_end  ,ip) - budg_stateG(s_wother_beg  ,ip)))/budg_fluxN(1,ip)
          write(iulog,'(89("-"),"|",23("-"))')
       end if
    end if
 end do

end subroutine MOSART_WaterBudget_Print

end module MOSART_Budgets_mod
