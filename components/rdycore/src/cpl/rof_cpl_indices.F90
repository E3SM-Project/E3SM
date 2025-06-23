module rof_cpl_indices

  use shr_sys_mod, only : shr_sys_abort
  implicit none

  save
  private

  public :: rof_cpl_indices_set

!
! !PUBLIC DATA MEMBERS:
!

! driver to rof
  integer, public :: index_x2r_Flrl_rofsur = 0    ! lnd->rof liquid surface runoff forcing from land
  integer, public :: index_x2r_Flrl_rofgwl = 0    ! lnd->rof liquid gwl runoff from land
  integer, public :: index_x2r_Flrl_rofsub = 0    ! lnd->rof liquid subsurface runoff from land
  integer, public :: index_x2r_Flrl_rofdto = 0    ! lnd->rof liquid direct to ocean runoff
  integer, public :: index_x2r_Flrl_rofi   = 0    ! lnd->rof ice runoff forcing from land
  integer, public :: index_x2r_Flrl_demand = 0    ! lnd->rof input total fluxes (<= 0)
  integer, public :: index_x2r_Flrl_Tqsur  = 0    ! lnd->rof Temperature of surface runoff
  integer, public :: index_x2r_Flrl_Tqsub  = 0    ! lnd->rof Temperature of subsurface runoff
  integer, public :: index_x2r_coszen_str  = 0    ! lnd->rof Cosine of Zenith
  integer, public :: index_x2r_So_ssh      = 0    ! ocn->rof ssh from ocean
  integer, public :: nflds_x2r = 0

! rof to driver
  integer, public :: index_r2x_Forr_rofl  = 0     ! rof->ocn liquid runoff to ocean
  integer, public :: index_r2x_Forr_rofi  = 0     ! rof->ocn ice runoff to ocean
  integer, public :: index_r2x_Flrr_flood = 0     ! rof->lnd flood runoff (>fthresh) back to land
  integer, public :: index_r2x_Flrr_volr = 0      ! rof->lnd volr total volume back to land
  integer, public :: index_r2x_Flrr_volrmch = 0   ! rof->lnd volr main channel back to land
  integer, public :: index_r2x_Flrr_supply = 0    ! rof->lnd supply flux for land use
  integer, public :: index_r2x_Flrr_deficit = 0   ! rof->lnd supply deficit
  integer, public :: nflds_r2x = 0

contains

  subroutine rof_cpl_indices_set()
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Set the coupler indices needed by the rdycore model coupler interface.
    !
    ! !USES:
    use seq_flds_mod , only : seq_flds_r2x_fields, seq_flds_x2r_fields
    use mct_mod      , only : mct_aVect, mct_aVect_init, mct_avect_indexra, &
                              mct_aVect_clean, mct_avect_nRattr
    !
    implicit none
    !
    ! !LOCAL VARIABLES:
    type(mct_aVect)   :: avtmp                            ! temporary
    character(len=32) :: subname = 'rof_cpl_indicies_set' ! subroutine name
    !-----------------------------------------------------------------------

    ! x2r
    call mct_aVect_init(avtmp, rList=seq_flds_x2r_fields, lsize=1)

    index_x2r_Flrl_rofsur = mct_avect_indexra(avtmp,'Flrl_rofsur') !'Flrl_rofsur')
    index_x2r_Flrl_rofgwl = mct_avect_indexra(avtmp,'Flrl_rofgwl')
    index_x2r_Flrl_rofsub = mct_avect_indexra(avtmp,'Flrl_rofsub')
    index_x2r_Flrl_rofdto = mct_avect_indexra(avtmp,'Flrl_rofdto',perrwith='quiet')
    index_x2r_Flrl_rofi   = mct_avect_indexra(avtmp,'Flrl_rofi')
    index_x2r_Flrl_demand = mct_avect_indexra(avtmp,'Flrl_demand')
    index_x2r_Flrl_Tqsur  = mct_avect_indexra(avtmp,'Flrl_Tqsur')
    index_x2r_Flrl_Tqsub  = mct_avect_indexra(avtmp,'Flrl_Tqsub')
    index_x2r_So_ssh      = mct_avect_indexra(avtmp,'So_ssh')

    nflds_x2r = mct_avect_nRattr(avtmp)

    call mct_aVect_clean(avtmp)

    ! r2x
    call mct_aVect_init(avtmp, rList=seq_flds_r2x_fields, lsize=1)

    index_r2x_Forr_rofl  = mct_avect_indexra(avtmp,'Forr_rofl')
    index_r2x_Forr_rofi  = mct_avect_indexra(avtmp,'Forr_rofi')
    index_r2x_Flrr_flood = mct_avect_indexra(avtmp,'Flrr_flood')
    index_r2x_Flrr_volr  = mct_avect_indexra(avtmp,'Flrr_volr')
    index_r2x_Flrr_volrmch = mct_avect_indexra(avtmp,'Flrr_volrmch')
    index_r2x_Flrr_supply = mct_avect_indexra(avtmp,'Flrr_supply')
    index_r2x_Flrr_deficit = mct_avect_indexra(avtmp,'Flrr_deficit')

    nflds_r2x = mct_avect_nRattr(avtmp)

    call mct_aVect_clean(avtmp)

  end subroutine rof_cpl_indices_set

end module rof_cpl_indices
