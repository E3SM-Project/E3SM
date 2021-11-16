module glc_cpl_indices
  
  use seq_flds_mod
  use mct_mod
  use glc_constants, only : glc_smb
  use shr_sys_mod  , only : shr_sys_abort

  implicit none

  SAVE
  public

  ! drv -> glc

  integer, public :: index_x2g_Sl_tsrf   = 0
  integer, public :: index_x2g_Flgl_qice = 0

  ! glc -> drv

  integer, public :: index_g2x_Fogg_rofi      = 0   ! frozen runoff -> ocn
  integer, public :: index_g2x_Figg_rofi      = 0   ! frozen runoff -> ice
  integer, public :: index_g2x_Fogg_rofl      = 0   ! liquid runoff -> ocn
  integer, public :: index_g2x_Sg_ice_covered = 0
  integer, public :: index_g2x_Sg_topo        = 0
  integer, public :: index_g2x_Flgg_hflx      = 0
  integer, public :: index_g2x_Sg_icemask     = 0
  integer, public :: index_g2x_Sg_icemask_coupled_fluxes = 0

contains

  subroutine glc_cpl_indices_set( )

    !-------------------------------------------------------------
    type(mct_aVect)   :: g2x      ! temporary
    type(mct_aVect)   :: x2g      ! temporary
    !-------------------------------------------------------------

    ! create temporary attribute vectors

    call mct_aVect_init(x2g, rList=seq_flds_x2g_fields, lsize=1)
    call mct_aVect_init(g2x, rList=seq_flds_g2x_fields, lsize=1)

    ! glc -> drv

    index_g2x_Fogg_rofi = mct_avect_indexra(g2x,'Fogg_rofi')
    index_g2x_Figg_rofi = mct_avect_indexra(g2x,'Figg_rofi')
    index_g2x_Fogg_rofl = mct_avect_indexra(g2x,'Fogg_rofl')
    index_g2x_Sg_ice_covered = mct_avect_indexra(g2x,'Sg_ice_covered')
    index_g2x_Sg_topo = mct_avect_indexra(g2x,'Sg_topo')
    index_g2x_Flgg_hflx = mct_avect_indexra(g2x,'Flgg_hflx')
    index_g2x_Sg_icemask = mct_avect_indexra(g2x,'Sg_icemask')
    index_g2x_Sg_icemask_coupled_fluxes = mct_avect_indexra(g2x,'Sg_icemask_coupled_fluxes')

    ! drv -> glc
    index_x2g_Sl_tsrf = mct_avect_indexra(x2g,'Sl_tsrf')
    index_x2g_Flgl_qice = mct_avect_indexra(x2g,'Flgl_qice')

    call mct_aVect_clean(x2g)
    call mct_aVect_clean(g2x)

    ! Set glc_smb
    ! true => get surface mass balance from CLM via coupler (in multiple elev classes)
    ! false => use PDD scheme in GLIMMER
    ! For now, we always use true

    glc_smb = .true.

  end subroutine glc_cpl_indices_set

end module glc_cpl_indices
