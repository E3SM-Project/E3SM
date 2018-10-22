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

  ! drv -> glc

  integer, public :: index_x2g_Flgl_qice = 0 !Jer: removed multiple fields, now that these are confined to coupler
  integer, public :: index_x2g_Fogo_qicel = 0 
  integer, public :: index_x2g_Sl_tsrf   = 0 !Jer: removed multiple fields, now that these are confined to coupler
  integer, public :: index_x2g_Fogo_qiceh = 0 !Jer:need to rename this everywhere...

  ! glc -> drv

  integer, public :: index_g2x_Fogg_rofi = 0   ! frozen runoff -> ocn
  integer, public :: index_g2x_Figg_rofi = 0   ! frozen runoff -> ice
  integer, public :: index_g2x_Fogg_rofl = 0   ! liquid runoff -> ocn
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
    integer           :: num 
    character(len= 2) :: cnum
    character(len=64) :: name
    !-------------------------------------------------------------

    ! create temporary attribute vectors

    call mct_aVect_init(x2g, rList=seq_flds_x2g_fields, lsize=1)
    call mct_aVect_init(g2x, rList=seq_flds_g2x_fields, lsize=1)

    ! glc -> drv

    index_g2x_Fogg_rofi = mct_avect_indexra(g2x,'Fogg_rofi')
    ! solid ice runoff to OCN  (kg/m^2/s = mm H2O/s) (flux: time-average)

    index_g2x_Figg_rofi = mct_avect_indexra(g2x,'Figg_rofi')
    ! solid ice runoff to ICE  (kg/m^2/s = mm H2O/s) (flux: time-average)

    index_g2x_Fogg_rofl = mct_avect_indexra(g2x,'Fogg_rofl')
    ! liquid runoff to OCN  (kg/m^2/s = mm H2O/s) (flux: time-average)

    index_g2x_Sg_ice_covered = mct_avect_indexra(g2x,'Sg_ice_covered')
    ! whether each grid cell is ice-covered (0 or 1) (state: current snapshot)
    ! The downscaling code (lnd -> glc) currently cannot handle fractional values.
    ! because this would require sending SMB for both ice-covered and ice-free
    ! points to a single glc cell. (This could be done, it would just require some
    ! 'rework.)  This mask is used by the coupler to determine where to send
    ! ice sheet landunit SMB vs. bare land landunit SMB.  This mask is used by
    ! CLM to differentiate bare land and ice sheet landunit types.

    index_g2x_Sg_topo = mct_avect_indexra(g2x,'Sg_topo')
    ! surface elevation (m), current snapshot

    index_g2x_Flgg_hflx = mct_avect_indexra(g2x,'Flgg_hflx')
    ! ice surface heat flux (W/m^2, positive down), time average

    index_g2x_Sg_icemask = mct_avect_indexra(g2x,'Sg_icemask')
    ! mask of ice sheet grid coverage (state: current snapshot)
    ! This mask (0 or 1) defines where the ice sheet model is able to handle
    ! surface mass balance, i.e., where there is ice (grounded or floating)
    ! or bare land.  It is calculated as locations where ice is present or
    ! land is present (i.e., excluding open ocean cells on the ice sheet grid).
    ! It is used by CLM (and maybe also the coupler?) to define where SMB
    ! should be sent.

    index_g2x_Sg_icemask_coupled_fluxes = mct_avect_indexra(g2x,'Sg_icemask_coupled_fluxes')
    ! mask of ice sheet grid coverage where we are potentially sending non-zero
    ! fluxes (state: current snapshot)
    ! This mask (0 or 1) indicates where the ice sheet model will actually apply SMB fluxes
    ! from the coupler.  This will be identical to "ice_sheet_grid_mask" for
    ! ice sheet instances that are configured to evolve.  If an ice sheet instance
    ! is configured to not evolve, then this field should be set to 0 everywhere
    ! in that domain.  Where this is 0 the coupler knows to send those fluxes
    ! directly to the runoff model because it is implied that the ice sheet model
    ! will not be applying them there.  This mask exists only to ensure
    ! conservation in configurations that have an ice sheet instance running
    ! but not evolving.  The GLC_TWO_WAY_COUPLING xml variable is what controls
    ! how this should be handled in CLM and the ice sheet model.


    ! drv -> glc

    index_x2g_Sl_tsrf = mct_avect_indexra(x2g,'Sl_tsrf')
    ! Surface temperature (deg C)

    index_x2g_Flgl_qice = mct_avect_indexra(x2g,'Flgl_qice')
    ! surface mass balance (kg/m^2/s)

    call mct_aVect_clean(x2g)
    call mct_aVect_clean(g2x)

  end subroutine glc_cpl_indices_set

end module glc_cpl_indices
