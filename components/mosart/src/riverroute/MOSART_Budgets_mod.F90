module MOSART_Budgets_mod
! Description: MOSART global water budget disgnostics
! 
! Added by Tian Zhou for E3SM v3, 11/06/2023. 
!-----------------------------------------------------------------------
  ! USES:
  use rof_cpl_indices, only : nt_rtm, nt_nliq, nt_nice
  use RtmVar         , only : iulog
  use shr_kind_mod   , only : r8 => shr_kind_r8

  implicit none
  private

  public MOSART_WaterBudget_Extraction

    !--- F for flux ---

  integer, parameter :: f_roff = 1
  integer, parameter :: f_rout = 2
  integer, parameter :: f_ioff = 3
  integer, parameter :: f_iout = 4
  integer, parameter :: f_irri = 5

  integer, parameter, public :: f_size = f_irri

  character(len=12),parameter :: fname(f_size) = &
       (/&
       '      runoff', &
       '     outflow', &
       '      frzrof', &
       '  frzoutflow', &
       '       irrig'  &
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
  integer, parameter :: s_w_errh2o      = 9

  integer, parameter, public :: s_size = s_w_errh2o

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
       '    error h2o'  &
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

  !----- formats -----
  character(*),parameter :: FA0= "('    ',12x,(3x,a10,2x),' | ',(3x,a10,2x))"
  character(*),parameter :: FF = "('    ',a12,f15.8,' | ',f18.2)"
  character(*),parameter :: FF2= "('    ',a12,a15,' | ',f18.2)"
  character(*),parameter :: FS = "('    ',a12,6(f18.2),18x,' | ',(f18.2))"
  character(*),parameter :: FS0= "('    ',12x,7(a18),' | ',(a18))"
  character(*),parameter :: FS2= "('    ',a12,54x,f18.2,54x,' | ',f18.2)"
  character(*),parameter :: FS3= "('    ',a12,7(f18.2),' | ',(f18.2))"

contains
  subroutine MOSART_WaterBudget_Extraction(budget_global, budget_terms_total, bv_volt_i, bv_volt_f, &
   budget_input, budget_output, budget_other)

   use RtmTimeManager, only : get_curr_date, get_nstep
   implicit none
   !
   integer,  intent(in) :: budget_terms_total
   real(r8), intent(in) :: budget_global(budget_terms_total,nt_rtm)    ! global budget sums
   integer,  intent(in) :: bv_volt_i  ! = 1  ! initial total volume
   integer,  intent(in) :: bv_volt_f  ! = 2  ! final   total volume
   real(r8), intent(in) :: budget_input, budget_output, budget_other
   integer :: ip
   character(*),parameter :: subName = '(MOSART_WaterBudget_Extraction) '

   ! budg_fluxG(f_roff,ip) = budget_global(xx,nt_nliq) ! 1 for liquid
   ! budg_fluxG(f_rout,ip) = budget_global(xx,nt_nliq) ! 1 for liquid

   do ip = 1,p_size
     budg_stateG(s_w_beg, ip) =  budget_global(bv_volt_i,nt_nliq) ! 
     budg_stateG(s_w_end, ip) =  budget_global(bv_volt_f,nt_nliq) ! 
   enddo

   write(iulog,'(2a,f22.6  )') trim(subname),'   state begin (TZ)   = ',budg_stateG(s_w_beg,1)
   write(iulog,'(2a,f22.6  )') trim(subname),'   state end   (TZ)   = ',budg_stateG(s_w_end,1)     

  end subroutine MOSART_WaterBudget_Extraction

end module MOSART_Budgets_mod