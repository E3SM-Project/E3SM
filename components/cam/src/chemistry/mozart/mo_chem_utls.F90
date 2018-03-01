
module mo_chem_utls

  private
  public :: get_spc_ndx, get_het_ndx, get_extfrc_ndx, get_rxt_ndx, get_inv_ndx

  save

contains

  integer function get_spc_ndx( spc_name )
    !-----------------------------------------------------------------------
    !     ... return overall species index associated with spc_name
    !-----------------------------------------------------------------------

    use chem_mods,     only : gas_pcnst
    use mo_tracname,   only : tracnam => solsym

    implicit none

    !-----------------------------------------------------------------------
    !     ... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: spc_name

    !-----------------------------------------------------------------------
    !     ... local variables
    !-----------------------------------------------------------------------
    integer :: m

    get_spc_ndx = -1
    do m = 1,gas_pcnst
       if( trim( spc_name ) == trim( tracnam(m) ) ) then
          get_spc_ndx = m
          exit
       end if
    end do

  end function get_spc_ndx

  integer function get_inv_ndx( invariant )
    !-----------------------------------------------------------------------
    !     ... return overall external frcing index associated with spc_name
    !-----------------------------------------------------------------------

    use chem_mods,  only : nfs, inv_lst

    implicit none

    !-----------------------------------------------------------------------
    !     ... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: invariant

    !-----------------------------------------------------------------------
    !     ... local variables
    !-----------------------------------------------------------------------
    integer :: m

    get_inv_ndx = -1
    do m = 1,nfs
       if( trim( invariant ) == trim( inv_lst(m) ) ) then
          get_inv_ndx = m
          exit
       end if
    end do

  end function get_inv_ndx

  integer function get_het_ndx( het_name )
    !-----------------------------------------------------------------------
    !     ... return overall het process index associated with spc_name
    !-----------------------------------------------------------------------

    use gas_wetdep_opts,only : gas_wetdep_method, gas_wetdep_list, gas_wetdep_cnt

    implicit none

    !-----------------------------------------------------------------------
    !     ... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: het_name

    !-----------------------------------------------------------------------
    !     ... local variables
    !-----------------------------------------------------------------------
    integer :: m

    get_het_ndx=-1

    do m=1,gas_wetdep_cnt

       if( trim( het_name ) == trim( gas_wetdep_list(m) ) ) then
          get_het_ndx = get_spc_ndx( gas_wetdep_list(m) )
          return
       endif
  
    enddo

  end function get_het_ndx

  integer function get_extfrc_ndx( frc_name )
    !-----------------------------------------------------------------------
    !     ... return overall external frcing index associated with spc_name
    !-----------------------------------------------------------------------

    use chem_mods,  only : extcnt, extfrc_lst

    implicit none

    !-----------------------------------------------------------------------
    !     ... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: frc_name

    !-----------------------------------------------------------------------
    !     ... local variables
    !-----------------------------------------------------------------------
    integer :: m

    get_extfrc_ndx = -1
    if( extcnt > 0 ) then
       do m = 1,max(1,extcnt)
          if( trim( frc_name ) == trim( extfrc_lst(m) ) ) then
             get_extfrc_ndx = m
             exit
          end if
       end do
    end if

  end function get_extfrc_ndx

  integer function get_rxt_ndx( rxt_tag )
    !-----------------------------------------------------------------------
    !     ... return overall external frcing index associated with spc_name
    !-----------------------------------------------------------------------

    use chem_mods,  only : rxt_tag_cnt, rxt_tag_lst, rxt_tag_map

    implicit none

    !-----------------------------------------------------------------------
    !     ... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: rxt_tag

    !-----------------------------------------------------------------------
    !     ... local variables
    !-----------------------------------------------------------------------
    integer :: m

    get_rxt_ndx = -1
    do m = 1,rxt_tag_cnt
       if( trim( rxt_tag ) == trim( rxt_tag_lst(m) ) ) then
          get_rxt_ndx = rxt_tag_map(m)
          exit
       end if
    end do

  end function get_rxt_ndx

end module mo_chem_utls
