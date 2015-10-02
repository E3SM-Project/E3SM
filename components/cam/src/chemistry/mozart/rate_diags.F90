!--------------------------------------------------------------------------------
! Manages writing reaction rates to history
!--------------------------------------------------------------------------------
module rate_diags

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_history,  only : fieldname_len
  use cam_history,  only : addfld,phys_decomp
  use cam_history,  only : outfld
  use chem_mods,    only : rxt_tag_cnt, rxt_tag_lst, rxt_tag_map
  use ppgrid,       only : pver

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

    character(len=64) :: name

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
       call addfld(rate_names(i), 'molecules/cm3/sec', pver,'A','reaction rate', phys_decomp)
    enddo

  end subroutine rate_diags_init

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine rate_diags_calc( rxt_rates, vmr, m, ncol, lchnk )

    use mo_rxt_rates_conv, only: set_rates

    real(r8), intent(inout) :: rxt_rates(:,:,:) ! 'molec/cm3/sec'
    real(r8), intent(in)    :: vmr(:,:,:)
    real(r8), intent(in)    :: m(:,:)           ! air density (molecules/cm3)
    integer,  intent(in)    :: ncol, lchnk

    integer :: i

    call set_rates( rxt_rates, vmr, ncol )
    
    do i = 1, rxt_tag_cnt
       ! convert from vmr/sec to molecules/cm3/sec
       rxt_rates(:ncol,:,rxt_tag_map(i)) = rxt_rates(:ncol,:,rxt_tag_map(i)) *  m(:,:)
       call outfld( rate_names(i), rxt_rates(:ncol,:,rxt_tag_map(i)), ncol, lchnk )

    enddo
  end subroutine rate_diags_calc

end module rate_diags
