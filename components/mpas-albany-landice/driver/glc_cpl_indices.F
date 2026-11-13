module glc_cpl_indices
  
  use seq_flds_mod
  use mct_mod
!  use glc_constants, only : glc_nec, glc_smb  ! TODO Will these be needed?  If so, need to add MPAS version
  use shr_sys_mod  , only : shr_sys_abort

  implicit none

  SAVE
  public

  integer , parameter, private:: glc_nzoc_max = 100

  ! Note that, in both the drv -> glc and the glc -> drv fields, index 0 means bare land

  integer, public :: index_x2g_Flgl_qice = 0 !Ice sheet surface mass balance
  integer, public :: index_x2g_Sl_tsrf   = 0 !Ice sheet upper surface boundary temperature

  integer, public :: index_x2g_So_tf3d(glc_nzoc_max)   = 0 !Ocean thermal forcing at specified z-levels
  integer, public :: index_x2g_So_tf3d_mask(glc_nzoc_max)   = 0 !mask of ocean thermal forcing at specified z-levels
  integer, public :: index_x2g_Foxo_ismf   = 0 !net ice-shelf melt flux from the ocean, normalized by the floating fraction
  integer, public :: index_x2g_So_tfrz_isf = 0 !ice-ocean interface temperature weighted by the ocean floating fraction
  integer, public :: index_x2g_So_liflfrac = 0 !ocean floating fraction used as the weight for So_tfrz_isf

  integer, public :: index_g2x_Fogg_rofi = 0 ! frozen runoff -> ocn
  integer, public :: index_g2x_Figg_rofi = 0 ! frozen runoff -> ice
  integer, public :: index_g2x_Fogg_rofl = 0 ! liquid runoff -> ocn
  integer, public :: index_g2x_Sg_topo = 0 ! surface elevation
  integer, public :: index_g2x_Flgg_hflx = 0 ! heat flux
  integer, public :: index_g2x_Sg_lithop = 0 !ice sheet lithostatic pressure
  integer, public :: index_g2x_Sg_icemask = 0 !complete grounded/floating ice mask
  integer, public :: index_g2x_Sg_ice_covered = 0 !ice covered ice mask
  integer, public :: index_g2x_Sg_icemask_grounded = 0 !grounded mask
  integer, public :: index_g2x_Sg_icemask_floating = 0 !floating mask
  integer, public :: index_g2x_Sg_below_sea_level_mask = 0 !mask of where the bed is below sea level
  integer, public :: index_g2x_Fogg_hcond = 0 !conductive heat flux into the ocean under floating ice
  integer, public :: index_g2x_Sg_icemask_coupled_fluxes = 0 !ice mask for coupled fluxes

contains

  subroutine glc_cpl_indices_set( )

    use glc_zocnclass_mod

    !-------------------------------------------------------------
    type(mct_aVect)   :: g2x      ! temporary
    type(mct_aVect)   :: x2g      ! temporary

    integer :: glc_nzoc
    integer :: iLev
    character(len=2) :: cnum
    character(len=64) :: varname
    !-------------------------------------------------------------

    glc_nzoc = glc_get_num_zocn_classes()

    ! create temporary attribute vectors

    call mct_aVect_init(x2g, rList=seq_flds_x2g_fields, lsize=1)
    call mct_aVect_init(g2x, rList=seq_flds_g2x_fields, lsize=1)

    !Following block are input boundary conditions for ice sheet model
    index_x2g_Flgl_qice =                 mct_avect_indexra(x2g,'Flgl_qice',perrwith='quiet')
    ! surface mass balance (kg/m^2/s)

    index_x2g_Sl_tsrf =                   mct_avect_indexra(x2g,'Sl_tsrf',perrwith='quiet')
    ! Surface temperature (deg C)

    if (glc_nzoc > 0) then
       do iLev = 1, glc_nzoc
          cnum = glc_zocnclass_as_string(iLev)
          varname = 'So_tf3d' // cnum
          index_x2g_So_tf3d(iLev) = mct_avect_indexra(x2g, trim(varname))

          varname = 'So_tf3d_mask' // cnum
          index_x2g_So_tf3d_mask(iLev) = mct_avect_indexra(x2g, trim(varname))
       enddo
    endif

    index_x2g_Foxo_ismf =                 mct_avect_indexra(x2g,'Foxo_ismf',perrwith='quiet')
    ! net ice-shelf melt flux from the ocean (kg/m^2/s, positive for melt),
    ! normalized by the floating fraction

    index_x2g_So_tfrz_isf =               mct_avect_indexra(x2g,'So_tfrz_isf',perrwith='quiet')
    index_x2g_So_liflfrac =               mct_avect_indexra(x2g,'So_liflfrac',perrwith='quiet')
    ! ice-ocean interface temperature (K) weighted by the ocean floating
    ! fraction, and that weight

    ! ice sheet lithostatic pressure (Pa) for the ocean under ice shelves
    index_g2x_Sg_lithop =                 mct_avect_indexra(g2x,'Sg_lithop',perrwith='quiet')

    !Following block are GLC outputs for other components.

    index_g2x_Fogg_rofi =                 mct_avect_indexra(g2x,'Fogg_rofi',perrwith='quiet')
    ! solid ice runoff to OCN  (kg/m^2/s = mm H2O/s) (flux: time-average)

    index_g2x_Figg_rofi =                 mct_avect_indexra(g2x,'Figg_rofi',perrwith='quiet')
    ! solid ice runoff to ICE  (kg/m^2/s = mm H2O/s) (flux: time-average)

    index_g2x_Fogg_rofl =                 mct_avect_indexra(g2x,'Fogg_rofl',perrwith='quiet')
    ! liquid runoff to OCN  (kg/m^2/s = mm H2O/s) (flux: time-average)

    index_g2x_Sg_topo =                   mct_avect_indexra(g2x,'Sg_topo',perrwith='quiet')
    ! surface elevation (m), current snapshot

    index_g2x_Flgg_hflx =                 mct_avect_indexra(g2x,'Flgg_hflx',perrwith='quiet')
    ! ice surface heat flux (W/m^2, positive down), time average

    index_g2x_Sg_icemask =                mct_avect_indexra(g2x,'Sg_icemask',perrwith='quiet')
    ! mask of ice sheet grid coverage (state: current snapshot)
    ! This mask (0 or 1) defines where the ice sheet model is able to handle
    ! surface mass balance, i.e., where there is ice (grounded or floating)
    ! or bare land.  It is calculated as locations where ice is present or
    ! land is present (i.e., excluding open ocean cells on the ice sheet grid).
    ! It is used by CLM (and maybe also the coupler?) to define where SMB
    ! should be sent.

    index_g2x_Sg_ice_covered = mct_avect_indexra(g2x,'Sg_ice_covered')
    ! whether each grid cell is ice-covered (0 or 1) (state: current snapshot)
    ! The downscaling code (lnd -> glc) currently cannot handle fractional values.
    ! because this would require sending SMB for both ice-covered and ice-free
    ! points to a single glc cell. (This could be done, it would just require some
    ! 'rework.)  This mask is used by the coupler to determine where to send
    ! ice sheet landunit SMB vs. bare land landunit SMB.  This mask is used by
    ! CLM to differentiate bare land and ice sheet landunit types.

    index_g2x_Sg_icemask_grounded =       mct_avect_indexra(g2x,'Sg_icemask_grounded',perrwith='quiet')
    index_g2x_Sg_icemask_floating =       mct_avect_indexra(g2x,'Sg_icemask_floating',perrwith='quiet')
    index_g2x_Sg_below_sea_level_mask =   mct_avect_indexra(g2x,'Sg_below_sea_level_mask',perrwith='quiet')
    index_g2x_Fogg_hcond =                mct_avect_indexra(g2x,'Fogg_hcond',perrwith='quiet')
    index_g2x_Sg_icemask_coupled_fluxes = mct_avect_indexra(g2x,'Sg_icemask_coupled_fluxes',perrwith='quiet')
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

    call mct_aVect_clean(x2g)
    call mct_aVect_clean(g2x)

  end subroutine glc_cpl_indices_set

end module glc_cpl_indices
