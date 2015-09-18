module glc_cpl_indices
  
  use seq_flds_mod
  use mct_mod
!  use glc_constants, only : glc_nec, glc_smb  ! TODO Will these be needed?  If so, need to add MPAS version
  use shr_sys_mod  , only : shr_sys_abort

  implicit none

  SAVE
  public

  integer , parameter, private:: glc_nec_max = 100

  ! Note that, in both the drv -> glc and the glc -> drv fields, index 0 means bare land

! TODO: Comment these variables!

! Jer: ocean coupler links from https://acme-climate.atlassian.net/wiki/pages/viewpage.action?pageId=16253763

  !!! drv -> glc

  ! lnd-derived
  integer, public :: index_x2g_Flgl_qice = 0 !Jer: removed multiple fields, now that these are confined to coupler
  integer, public :: index_x2g_Sl_tsrf   = 0 !Jer: removed multiple fields, now that these are confined to coupler
  
  ! ocn-derived (but potentially coupler-mediated)
  integer, public :: index_x2g_Fogx_qicel = 0 !Subshelf mass flux
  integer, public :: index_x2g_Fogx_qiceh = 0 !Subshelf heat flux Jer:need to rename this everywhere, from So_tglc
  integer, public :: index_x2g_So_phieff = 0

  !!! glc -> drv

  ! to lnd
  integer, public :: index_g2x_Fogg_rofi = 0   ! frozen runoff -> ocn
  integer, public :: index_g2x_Figg_rofi = 0   ! frozen runoff -> ice
  integer, public :: index_g2x_Fogg_rofl = 0   ! liquid runoff -> ocn
  integer, public :: index_g2x_Sg_frac(0:glc_nec_max)   = 0
  integer, public :: index_g2x_Sg_topo(0:glc_nec_max)   = 0
  integer, public :: index_g2x_Flgg_hflx(0:glc_nec_max) = 0

  ! to ocn (but potentially coupler-mediated)

  integer, public :: index_g2x_Sg_tbot = 0 !bottom-layer temperature
  integer, public :: index_g2x_Sg_dztbot = 0 !distance from bottom-layer temperature to ice-ocean interface (to compute T gradient)
  integer, public :: index_g2x_Sg_lithop = 0 !ice sheet lithostatic pressure

  ! common

  integer, public :: index_g2x_Sg_icemask = 0
  integer, public :: index_g2x_Sg_icemask_grounded = 0 !grounded mask
  integer, public :: index_g2x_Sg_icemask_floating = 0 !floating mask
  integer, public :: index_g2x_Sg_icemask_coupled_fluxes = 0
  
contains

  subroutine glc_cpl_indices_set( )

    !-------------------------------------------------------------
    type(mct_aVect)   :: g2x      ! temporary
    type(mct_aVect)   :: x2g      ! temporary
    integer           :: num 
    character(len= 2) :: cnum
    character(len=64) :: name
    !-------------------------------------------------------------

    ! create temporary attribute vectors

    call mct_aVect_init(x2g, rList=seq_flds_x2g_fields, lsize=1)
    call mct_aVect_init(g2x, rList=seq_flds_g2x_fields, lsize=1)

!    glc_nec = 0

    ! glc -> drv

    index_g2x_Fogg_rofi = mct_avect_indexra(g2x,'Fogg_rofi',perrwith='quiet')
    index_g2x_Figg_rofi = mct_avect_indexra(g2x,'Figg_rofi',perrwith='quiet')
    index_g2x_Fogg_rofl = mct_avect_indexra(g2x,'Fogg_rofl',perrwith='quiet')
    
    index_g2x_Sg_tbot =   mct_avect_indexra(g2x,'Sg_tbot',perrwith='quiet')
    index_g2x_Sg_dztbot = mct_avect_indexra(g2x,'Sg_dztbot',perrwith='quiet')
    index_g2x_Sg_lithop = mct_avect_indexra(g2x,'Sg_lithop',perrwith='quiet')

    index_g2x_Sg_icemask =                mct_avect_indexra(g2x,'Sg_icemask',perrwith='quiet')
    index_g2x_Sg_icemask_grounded =       mct_avect_indexra(g2x,'Sg_icemask_grounded',perrwith='quiet')
    index_g2x_Sg_icemask_floating =       mct_avect_indexra(g2x,'Sg_icemask_floating',perrwith='quiet')
    index_g2x_Sg_icemask_coupled_fluxes = mct_avect_indexra(g2x,'Sg_icemask_coupled_fluxes',perrwith='quiet')

    ! drv -> glc

    index_x2g_Fogx_qicel = mct_avect_indexra(x2g,'Fogx_qicel',perrwith='quiet')
    index_x2g_Fogx_qiceh = mct_avect_indexra(x2g,'Fogx_qiceh',perrwith='quiet')
    index_x2g_So_phieff =  mct_avect_indexra(x2g,'So_phieff',perrwith='quiet')



    call mct_aVect_clean(x2g)
    call mct_aVect_clean(g2x)

  end subroutine glc_cpl_indices_set

end module glc_cpl_indices
