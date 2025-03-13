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
  integer, public :: index_x2r_Flrl_demand = 0  ! lnd->rof input total fluxes (<= 0)
  integer, public :: index_x2r_Flrl_Tqsur  = 0  ! lnd->rof Temperature of surface runoff
  integer, public :: index_x2r_Flrl_Tqsub  = 0  ! lnd->rof Temperature of subsurface runoff
  integer, public :: index_x2r_Sa_tbot = 0      ! atm->rof air temperature
  integer, public :: index_x2r_Sa_pbot = 0      ! atm->rof surface pressure
  integer, public :: index_x2r_Sa_u    = 0      ! atm->rof zonal velocity
  integer, public :: index_x2r_Sa_v    = 0      ! atm->rof merid velocity
  integer, public :: index_x2r_Sa_shum = 0      ! atm->rof specific humidity
  integer, public :: index_x2r_Faxa_lwdn  = 0   ! atm->rof longwave down flux
  integer, public :: index_x2r_Faxa_swvdr = 0   ! atm->rof shorwave visible direct flux
  integer, public :: index_x2r_Faxa_swvdf = 0   ! atm->rof shorwave visible diffus flux
  integer, public :: index_x2r_Faxa_swndr = 0   ! atm->rof shorwave near-ir direct flux
  integer, public :: index_x2r_Faxa_swndf = 0   ! atm->rof shorwave near-ir diffus flux
  integer, public :: index_x2r_Flrl_rofmud = 0  ! lnd->rof input suspended sediment flux from soil erosion
  integer, public :: index_x2r_Flrl_inundinf = 0! lnd->rof infiltration from floodplain inundation
  integer, public :: nflds_x2r = 0

  integer, public :: index_x2r_coszen_str  = 0   ! lnd->rof Cosine of Zenith
  integer, public :: index_x2r_So_ssh = 0        ! ocn->rof ssh from ocean

  !TODO - nt_rtm and rtm_tracers need to be removed and set by access to the index array
  integer, parameter, public :: nt_rtm = 4    ! number of tracers
  character(len=3), parameter, public :: rtm_tracers(nt_rtm) =  (/'LIQ','ICE','MUD','SAN'/)
  integer, parameter, public :: nt_nliq = 1    ! number of tracers
  integer, parameter, public :: nt_nice = 2    ! number of tracers
  integer, parameter, public :: nt_nmud = 3    ! number of tracers
  integer, parameter, public :: nt_nsan = 4    ! number of tracers

  !Routing methods used for the main-channel
  integer, parameter, public :: KW = 1         ! kinematic wave routing method
  integer, parameter, public :: DW = 2         ! diffusion wave routing method

  ! roff to driver (part of land for now) (optional if ROF is off)

  integer, public :: index_r2x_Forr_rofl  = 0   ! rof->ocn liquid runoff to ocean
  integer, public :: index_r2x_Forr_rofi  = 0   ! rof->ocn ice runoff to ocean
  integer, public :: index_r2x_Forr_rofDIN = 0  ! rof->ocn DIN runoff to ocean
  integer, public :: index_r2x_Forr_rofDIP = 0  ! rof->ocn DIP runoff to ocean
  integer, public :: index_r2x_Forr_rofDON = 0  ! rof->ocn DON runoff to ocean
  integer, public :: index_r2x_Forr_rofDOP = 0  ! rof->ocn DOP runoff to ocean
  integer, public :: index_r2x_Forr_rofDOC = 0  ! rof->ocn DOC runoff to ocean
  integer, public :: index_r2x_Forr_rofPP  = 0  ! rof->ocn PP  runoff to ocean
  integer, public :: index_r2x_Forr_rofDSi = 0  ! rof->ocn DSi runoff to ocean
  integer, public :: index_r2x_Forr_rofPOC = 0  ! rof->ocn POC runoff to ocean
  integer, public :: index_r2x_Forr_rofPN  = 0  ! rof->ocn PN  runoff to ocean
  integer, public :: index_r2x_Forr_rofDIC = 0  ! rof->ocn DIC runoff to ocean
  integer, public :: index_r2x_Forr_rofFe  = 0  ! rof->ocn Fe  runoff to ocean
  integer, public :: index_r2x_Flrr_flood = 0   ! rof->lnd flood runoff (>fthresh) back to land
  integer, public :: index_r2x_Flrr_volr = 0    ! rof->lnd volr total volume back to land
  integer, public :: index_r2x_Flrr_volrmch = 0 ! rof->lnd volr main channel back to land
  integer, public :: index_r2x_Flrr_supply = 0  ! rof->lnd supply flux for land use
  integer, public :: index_r2x_Flrr_deficit = 0 ! rof->lnd supply deficit
  integer, public :: index_r2x_Sr_h2orof      = 0  ! rof->lnd floodplain inundation volume
  integer, public :: index_r2x_Sr_frac_h2orof = 0  ! rof->lnd floodplain inundation fraction
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
    use seq_flds_mod  , only: seq_flds_r2x_fields, seq_flds_x2r_fields, rof_heat, &
                              rof2ocn_nutrients, lnd_rof_two_way, ocn_rof_two_way, &
                              rof_sed
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
    index_x2r_Flrl_demand = mct_avect_indexra(avtmp,'Flrl_demand')
    index_x2r_Flrl_Tqsur  = mct_avect_indexra(avtmp,'Flrl_Tqsur')
    index_x2r_Flrl_Tqsub  = mct_avect_indexra(avtmp,'Flrl_Tqsub')
    if (ocn_rof_two_way) then
      index_x2r_So_ssh      = mct_avect_indexra(avtmp,'So_ssh')
    endif
    if (rof_heat) then
      index_x2r_Sa_tbot     = mct_avect_indexra(avtmp,'Sa_tbot')
      index_x2r_Sa_pbot     = mct_avect_indexra(avtmp,'Sa_pbot')
      index_x2r_Sa_u        = mct_avect_indexra(avtmp,'Sa_u')
      index_x2r_Sa_v        = mct_avect_indexra(avtmp,'Sa_v')
      index_x2r_Sa_shum     = mct_avect_indexra(avtmp,'Sa_shum')
      index_x2r_Faxa_lwdn   = mct_avect_indexra(avtmp,'Faxa_lwdn')
      index_x2r_Faxa_swvdr  = mct_avect_indexra(avtmp,'Faxa_swvdr')
      index_x2r_Faxa_swvdf  = mct_avect_indexra(avtmp,'Faxa_swvdf')
      index_x2r_Faxa_swndr  = mct_avect_indexra(avtmp,'Faxa_swndr')
      index_x2r_Faxa_swndf  = mct_avect_indexra(avtmp,'Faxa_swndf')
    endif

    index_x2r_coszen_str  = mct_avect_indexra(avtmp,'coszen_str')
	if (rof_sed) then
        index_x2r_Flrl_rofmud = mct_avect_indexra(avtmp,'Flrl_rofmud')
	end if
    if (lnd_rof_two_way) then
      index_x2r_Flrl_inundinf =  mct_avect_indexra(avtmp,'Flrl_inundinf')
    endif

    nflds_x2r = mct_avect_nRattr(avtmp)

    call mct_aVect_clean(avtmp)

    ! r2x

    call mct_aVect_init(avtmp, rList=seq_flds_r2x_fields, lsize=1)

    index_r2x_Forr_rofl  = mct_avect_indexra(avtmp,'Forr_rofl')
    index_r2x_Forr_rofi  = mct_avect_indexra(avtmp,'Forr_rofi')
    if (rof2ocn_nutrients) then
       index_r2x_Forr_rofDIN = mct_avect_indexra(avtmp,'Forr_rofDIN')
       index_r2x_Forr_rofDIP = mct_avect_indexra(avtmp,'Forr_rofDIP')
       index_r2x_Forr_rofDON = mct_avect_indexra(avtmp,'Forr_rofDON')
       index_r2x_Forr_rofDOP = mct_avect_indexra(avtmp,'Forr_rofDOP')
       index_r2x_Forr_rofDOC = mct_avect_indexra(avtmp,'Forr_rofDOC')
       index_r2x_Forr_rofPP  = mct_avect_indexra(avtmp,'Forr_rofPP')
       index_r2x_Forr_rofDSi = mct_avect_indexra(avtmp,'Forr_rofDSi')
       index_r2x_Forr_rofPOC = mct_avect_indexra(avtmp,'Forr_rofPOC')
       index_r2x_Forr_rofPN  = mct_avect_indexra(avtmp,'Forr_rofPN')
       index_r2x_Forr_rofDIC = mct_avect_indexra(avtmp,'Forr_rofDIC')
       index_r2x_Forr_rofFe  = mct_avect_indexra(avtmp,'Forr_rofFe')
    endif
    index_r2x_Flrr_flood = mct_avect_indexra(avtmp,'Flrr_flood')
    index_r2x_Flrr_volr  = mct_avect_indexra(avtmp,'Flrr_volr')
    index_r2x_Flrr_volrmch = mct_avect_indexra(avtmp,'Flrr_volrmch')
    index_r2x_Flrr_supply = mct_avect_indexra(avtmp,'Flrr_supply')
    index_r2x_Flrr_deficit = mct_avect_indexra(avtmp,'Flrr_deficit')

    if (lnd_rof_two_way) then
      index_r2x_Sr_h2orof       = mct_avect_indexra(avtmp,'Sr_h2orof')
      index_r2x_Sr_frac_h2orof  = mct_avect_indexra(avtmp,'Sr_frac_h2orof')
    endif
    
    nflds_r2x = mct_avect_nRattr(avtmp)

    call mct_aVect_clean(avtmp)

  end subroutine rof_cpl_indices_set

end module rof_cpl_indices


