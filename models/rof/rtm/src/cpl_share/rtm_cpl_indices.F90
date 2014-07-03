module rtm_cpl_indices
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: rtm_cpl_indices
!
! !DESCRIPTION:
!    Module containing the indices for the fields passed between RTM and
!    the driver. 
!
! !USES:

  use shr_sys_mod,    only : shr_sys_abort
  implicit none

  SAVE
  private                              ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:

  public :: rtm_cpl_indices_set        ! Set the coupler indices

!
! !PUBLIC DATA MEMBERS:
!
  integer, public :: index_x2r_Flrl_rofliq  = 0   ! lnd->rtm liquid runoff forcing from land
  integer, public :: index_x2r_Flrl_rofice  = 0   ! lnd->rtm ice runoff forcing from land

  integer, public :: nflds_x2r = 0

  !TODO - nt_rtm and rtm_tracers need to be removed and set by access to the index array
  integer, parameter, public :: nt_rtm = 2    ! number of tracers
  character(len=3), parameter, public :: rtm_tracers(nt_rtm) =  (/'LIQ','ICE'/)

  ! roff to driver (part of land for now) (optional if RTM is off)

  integer, public :: index_r2x_Forr_roff  = 0   ! rtm->ocn liquid runoff to ocean
  integer, public :: index_r2x_Forr_ioff  = 0   ! rtm->ocn ice runoff to ocean
  integer, public :: index_r2x_Flrr_flood = 0   ! rtm->lnd flood runoff (>fthresh) back to land
  integer, public :: index_r2x_Slrr_volr = 0   ! rtm->lnd volr back to land
  integer, public :: nflds_r2x = 0

!=======================================================================
contains

!=======================================================================


  subroutine rtm_cpl_indices_set( )


    !-----------------------------------------------------------------------
    ! !DESCRIPTION: 
    ! Set the coupler indices needed by the rof model coupler interface.
    ! runoff - (rtm -> ocn) and (rtm->lnd)
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
    character(len=32) :: subname = 'rtm_cpl_indices_set'  ! subroutine name
    !-----------------------------------------------------------------------

    ! x2r

    call mct_aVect_init(avtmp, rList=seq_flds_x2r_fields, lsize=1)

    index_x2r_Flrl_rofliq = mct_avect_indexra(avtmp,'Flrl_rofliq')
    index_x2r_Flrl_rofice = mct_avect_indexra(avtmp,'Flrl_rofice')

    nflds_x2r = mct_avect_nRattr(avtmp)

    call mct_aVect_clean(avtmp)

    ! r2x

    call mct_aVect_init(avtmp, rList=seq_flds_r2x_fields, lsize=1)

    index_r2x_Forr_roff  = mct_avect_indexra(avtmp,'Forr_roff')
    index_r2x_Forr_ioff  = mct_avect_indexra(avtmp,'Forr_ioff')
    index_r2x_Flrr_flood = mct_avect_indexra(avtmp,'Flrr_flood')
    index_r2x_Slrr_volr = mct_avect_indexra(avtmp,'Slrr_volr')

    nflds_r2x = mct_avect_nRattr(avtmp)

    call mct_aVect_clean(avtmp)

  end subroutine rtm_cpl_indices_set

end module rtm_cpl_indices


