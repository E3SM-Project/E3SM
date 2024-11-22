!--------------------------------------------------------------------------------
! Manages writing reaction rates to history
!--------------------------------------------------------------------------------
module rate_diags

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_history,  only : fieldname_len
  use cam_history,  only : addfld, horiz_only, add_default
  use cam_history,  only : outfld
  use chem_mods,    only : rxt_tag_cnt, rxt_tag_lst, rxt_tag_map
  use ppgrid,       only : pver
  use mo_constants, only : rgrav
  use phys_control, only : phys_getopts

  implicit none
  private 
  public :: rate_diags_init
  public :: rate_diags_calc

  character(len=fieldname_len) :: rate_names(rxt_tag_cnt)

contains

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine rate_diags_init

    integer :: i, len, pos
    logical  :: history_UCIgaschmbudget_2D ! output 2D gas chemistry tracer concentrations and tendencies
    logical  :: history_UCIgaschmbudget_2D_levels ! output 2D gas chemistry tracer concentrations and tendencies within certain layers

    character(len=64) :: name

     call phys_getopts( history_UCIgaschmbudget_2D_out = history_UCIgaschmbudget_2D, &
                       history_UCIgaschmbudget_2D_levels_out = history_UCIgaschmbudget_2D_levels)

    do i = 1,rxt_tag_cnt
       pos = 0
       pos = index(rxt_tag_lst(i),'tag_')
       if (pos <= 0) pos = index(rxt_tag_lst(i),'usr_')
       if (pos <= 0) pos = index(rxt_tag_lst(i),'cph_')
       if (pos <= 0) pos = index(rxt_tag_lst(i),'ion_')
       if (pos>0) then
          name = 'r_'//trim(rxt_tag_lst(i)(5:))
       else
          name = 'r_'//trim(rxt_tag_lst(i)(1:))
       endif
       len = min(fieldname_len,len_trim(name))
       rate_names(i) = trim(name(1:len))
       call addfld(rate_names(i), (/ 'lev' /),'A', 'molecules/cm3/sec','reaction rate')
    enddo

    if ( history_UCIgaschmbudget_2D ) then
       call addfld('r_lch4_2D', horiz_only, 'A', 'kg/m2/s', 'CH4 vertically integrated reatction rate ')
    endif
    if ( history_UCIgaschmbudget_2D_levels) then 
       call addfld('r_lch4_L1', horiz_only, 'A', 'kg/m2/s', 'CH4 vertically integrated reaction rate from top-of-model to 100 hPa')
       call addfld('r_lch4_L2', horiz_only, 'A', 'kg/m2/s', 'CH4 vertically integrated reaction rate from 100 to 267 hPa')
       call addfld('r_lch4_L3', horiz_only, 'A', 'kg/m2/s', 'CH4 vertically integrated reaction rate from 267 hPa to 856 hPa')
       call addfld('r_lch4_L4', horiz_only, 'A', 'kg/m2/s', 'CH4 vertically integrated reaction rate from 856 hPa to surface')
    endif
    if ( history_UCIgaschmbudget_2D ) then
       call add_default( 'r_lch4_2D', 1, ' ' )
    endif
    if ( history_UCIgaschmbudget_2D_levels ) then
       call add_default( 'r_lch4_L1', 1, ' ' )
       call add_default( 'r_lch4_L2', 1, ' ' )
       call add_default( 'r_lch4_L3', 1, ' ' )
       call add_default( 'r_lch4_L4', 1, ' ' )
    endif
    
  end subroutine rate_diags_init

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine rate_diags_calc( rxt_rates, vmr, m, ncol, lchnk, pver, pdeldry, mbar )

    use mo_rxt_rates_conv, only: set_rates
    use chem_mods,    only : gas_pcnst, rxntot

    real(r8), intent(inout) :: rxt_rates(:,:,:) ! 'molec/cm3/sec'
    real(r8), intent(in)    :: vmr(:,:,:)
    real(r8), intent(in)    :: m(:,:)           ! air density (molecules/cm3)
    integer,  intent(in)    :: ncol, lchnk, pver
    real(r8), intent(in)    :: pdeldry(:,:)
    real(r8), intent(in)    :: mbar(:,:)

    integer :: i, k

    real(r8)  :: rxt_rates_vmr(ncol,pver,max(1,rxntot)) ! 'vmr/sec'
    real(r8)  :: wrk(ncol,pver)
    real(r8)  :: wrk_sum(ncol)
    real(r8)  :: adv_mass_ch4
    logical  :: history_UCIgaschmbudget_2D ! output 2D gas chemistry tracer concentrations and tendencies
    logical  :: history_UCIgaschmbudget_2D_levels ! output 2D gas chemistry tracer concentrations and tendencies within certain layers
    integer  :: gaschmbudget_2D_L1_s ! Start layer of L1 for gas chemistry tracer budget 
    integer  :: gaschmbudget_2D_L1_e ! End layer of L1 for gas chemistry trracer budget
    integer  :: gaschmbudget_2D_L2_s
    integer  :: gaschmbudget_2D_L2_e
    integer  :: gaschmbudget_2D_L3_s
    integer  :: gaschmbudget_2D_L3_e
    integer  :: gaschmbudget_2D_L4_s
    integer  :: gaschmbudget_2D_L4_e

    call phys_getopts( history_UCIgaschmbudget_2D_out = history_UCIgaschmbudget_2D, &
                       history_UCIgaschmbudget_2D_levels_out = history_UCIgaschmbudget_2D_levels, &
                       gaschmbudget_2D_L1_s_out = gaschmbudget_2D_L1_s, &
                       gaschmbudget_2D_L1_e_out = gaschmbudget_2D_L1_e, &
                       gaschmbudget_2D_L2_s_out = gaschmbudget_2D_L2_s, &
                       gaschmbudget_2D_L2_e_out = gaschmbudget_2D_L2_e, &
                       gaschmbudget_2D_L3_s_out = gaschmbudget_2D_L3_s, &
                       gaschmbudget_2D_L3_e_out = gaschmbudget_2D_L3_e, &
                       gaschmbudget_2D_L4_s_out = gaschmbudget_2D_L4_s, &
                       gaschmbudget_2D_L4_e_out = gaschmbudget_2D_L4_e )


    rxt_rates_vmr = 0._r8
    adv_mass_ch4 = 16._r8

    call set_rates( rxt_rates, vmr, ncol )

    rxt_rates_vmr = rxt_rates

    do i = 1, rxt_tag_cnt
       ! convert from vmr/sec to molecules/cm3/sec
       rxt_rates(:ncol,:,rxt_tag_map(i)) = rxt_rates(:ncol,:,rxt_tag_map(i)) *  m(:,:)
       call outfld( rate_names(i), rxt_rates(:ncol,:,rxt_tag_map(i)), ncol, lchnk )

       if (history_UCIgaschmbudget_2D .or. history_UCIgaschmbudget_2D_levels) then

       if (rate_names(i) .eq. 'r_lch4') then
          !kg/m2/sec
          wrk(:ncol,:) = adv_mass_ch4*rxt_rates_vmr(:ncol,:,rxt_tag_map(i))/mbar(:ncol,:) &
                                *pdeldry(:ncol,:)*rgrav

       if (history_UCIgaschmbudget_2D_levels) then
          wrk_sum(:ncol) = 0.0_r8
            do k = gaschmbudget_2D_L1_s, gaschmbudget_2D_L1_e
               wrk_sum(:ncol) = wrk_sum(:ncol) + wrk(:ncol,k)
            enddo
            call outfld( 'r_lch4_L1', wrk_sum(:ncol), ncol ,lchnk )
          wrk_sum(:ncol) = 0.0_r8
            do k = gaschmbudget_2D_L2_s, gaschmbudget_2D_L2_e
               wrk_sum(:ncol) = wrk_sum(:ncol) + wrk(:ncol,k)
            enddo
            call outfld( 'r_lch4_L2', wrk_sum(:ncol), ncol ,lchnk )
          wrk_sum(:ncol) = 0.0_r8
            do k = gaschmbudget_2D_L3_s, gaschmbudget_2D_L3_e
               wrk_sum(:ncol) = wrk_sum(:ncol) + wrk(:ncol,k)
            enddo
            call outfld( 'r_lch4_L3', wrk_sum(:ncol), ncol ,lchnk )
          wrk_sum(:ncol) = 0.0_r8
            do k = gaschmbudget_2D_L4_s, gaschmbudget_2D_L4_e
               wrk_sum(:ncol) = wrk_sum(:ncol) + wrk(:ncol,k)
            enddo
            call outfld( 'r_lch4_L4', wrk_sum(:ncol), ncol ,lchnk )
            do k=2,pver
               wrk(:ncol,1) = wrk(:ncol,1) + wrk(:ncol,k)
            enddo
            call outfld( 'r_lch4_2D', wrk(:ncol,1), ncol ,lchnk )
       elseif (history_UCIgaschmbudget_2D) then          
            do k=2,pver
               wrk(:ncol,1) = wrk(:ncol,1) + wrk(:ncol,k)
            enddo
            call outfld( 'r_lch4_2D', wrk(:ncol,1), ncol ,lchnk )
       endif

       endif

       endif
    enddo

  end subroutine rate_diags_calc

end module rate_diags
