module rof_cpl_indices
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: rof_cpl_indices
!
! !DESCRIPTION:
!    Module containing the indices for the fields passed between ROF and
!    the driver. 
!
! !USES:

  use shr_sys_mod,    only : shr_sys_abort
  implicit none

  SAVE
  private                              ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:

  public :: rof_cpl_indices_set        ! Set the coupler indices

!
! !PUBLIC DATA MEMBERS:
!
  integer, public :: index_x2r_Flrl_rofsur = 0  ! lnd->rof liquid surface runoff forcing from land
  integer, public :: index_x2r_Flrl_rofgwl = 0  ! lnd->rof liquid gwl runoff from land
  integer, public :: index_x2r_Flrl_rofsub = 0  ! lnd->rof liquid subsurface runoff from land
  integer, public :: index_x2r_Flrl_rofdto = 0  ! lnd->rof liquid direct to ocean runoff
  integer, public :: index_x2r_Flrl_rofi  = 0   ! lnd->rof ice runoff forcing from land

  integer, public :: nflds_x2r = 0

  !TODO - nt_rtm and rtm_tracers need to be removed and set by access to the index array
  integer, parameter, public :: nt_rtm = 2    ! number of tracers
  character(len=3), parameter, public :: rtm_tracers(nt_rtm) =  (/'LIQ','ICE'/)

  ! roff to driver (part of land for now) (optional if ROF is off)

  integer, public :: index_r2x_Forr_rofl  = 0   ! rof->ocn liquid runoff to ocean
  integer, public :: index_r2x_Forr_rofi  = 0   ! rof->ocn ice runoff to ocean
  integer, public :: index_r2x_Flrr_flood = 0   ! rof->lnd flood runoff (>fthresh) back to land
  integer, public :: index_r2x_Flrr_volr = 0    ! rof->lnd volr total volume back to land
  integer, public :: index_r2x_Flrr_volrmch = 0 ! rof->lnd volr main channel back to land
  integer, public :: nflds_r2x = 0

!=======================================================================
contains

!=======================================================================


  subroutine rof_cpl_indices_set( )


    !-----------------------------------------------------------------------
    ! !DESCRIPTION: 
    ! Set the coupler indices needed by the rof model coupler interface.
    ! runoff - (rof -> ocn) and (rof->lnd)
    !
    ! !USES:
    use seq_flds_mod  , only: seq_flds_r2x_fields, seq_flds_x2r_fields
    use mct_mod       , only: mct_aVect, mct_aVect_init, mct_avect_indexra, &
                              mct_aVect_clean, mct_avect_nRattr
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !REVISION HISTORY:
    ! Author: Mariana Vertenstein
    !
    ! !LOCAL VARIABLES:
    type(mct_aVect)   :: avtmp      ! temporary av
    character(len=32) :: subname = 'rof_cpl_indices_set'  ! subroutine name
    !-----------------------------------------------------------------------

    ! x2r

    call mct_aVect_init(avtmp, rList=seq_flds_x2r_fields, lsize=1)

    index_x2r_Flrl_rofsur = mct_avect_indexra(avtmp,'Flrl_rofsur') !'Flrl_rofsur')
    index_x2r_Flrl_rofgwl = mct_avect_indexra(avtmp,'Flrl_rofgwl')
    index_x2r_Flrl_rofsub = mct_avect_indexra(avtmp,'Flrl_rofsub')
    index_x2r_Flrl_rofdto = mct_avect_indexra(avtmp,'Flrl_rofdto',perrwith='quiet')
    index_x2r_Flrl_rofi   = mct_avect_indexra(avtmp,'Flrl_rofi')

    nflds_x2r = mct_avect_nRattr(avtmp)

    call mct_aVect_clean(avtmp)

    ! r2x

    call mct_aVect_init(avtmp, rList=seq_flds_r2x_fields, lsize=1)

    index_r2x_Forr_rofl  = mct_avect_indexra(avtmp,'Forr_rofl')
    index_r2x_Forr_rofi  = mct_avect_indexra(avtmp,'Forr_rofi')
    index_r2x_Flrr_flood = mct_avect_indexra(avtmp,'Flrr_flood')
    index_r2x_Flrr_volr  = mct_avect_indexra(avtmp,'Flrr_volr')
    index_r2x_Flrr_volrmch = mct_avect_indexra(avtmp,'Flrr_volrmch')

    nflds_r2x = mct_avect_nRattr(avtmp)

    call mct_aVect_clean(avtmp)

  end subroutine rof_cpl_indices_set

end module rof_cpl_indices


