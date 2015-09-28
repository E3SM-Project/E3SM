module glc_cpl_indices
  
  use seq_flds_mod
  use mct_mod
  use glc_constants, only : glc_nec, glc_smb  ! TODO Will these be needed?  If so, need to add MPAS version
  use shr_sys_mod  , only : shr_sys_abort

  implicit none

  SAVE
  public

  integer , parameter, private:: glc_nec_max = 100

  ! Note that, in both the drv -> glc and the glc -> drv fields, index 0 means bare land

  integer, public :: index_x2g_Flgl_qice = 0 !Ice sheet surface mass balance
  integer, public :: index_x2g_Sl_tsrf   = 0 !Ice sheet upper surface boundary temperature
  
  integer, public :: index_x2g_So_blt    = 0 !Ice shelf boundary layer ocean temperature
  integer, public :: index_x2g_So_bls    = 0 !Ice shelf boundary layer ocean salinity
  integer, public :: index_x2g_So_htv    = 0 !Ice shelf ocean heat transfer velocity
  integer, public :: index_x2g_So_hsv    = 0 !Ice shelf ocean salinity transfer velocity
  integer, public :: index_x2g_So_phieff = 0 !Ocean effective pressure
  integer, public :: index_x2g_Fogx_qiceli = 0 !Subshelf mass flux
  integer, public :: index_x2g_Fogx_qicehi = 0 !Subshelf heat flux for the ice sheet

  integer, public :: index_g2x_Fogg_rofi = 0 ! frozen runoff -> ocn
  integer, public :: index_g2x_Figg_rofi = 0 ! frozen runoff -> ice
  integer, public :: index_g2x_Fogg_rofl = 0 ! liquid runoff -> ocn
  integer, public :: index_g2x_Sg_frac(0:glc_nec_max)   = 0 !Subgrid elevation fraction
  integer, public :: index_g2x_Sg_topo(0:glc_nec_max)   = 0 !Subgrid elevation
  integer, public :: index_g2x_Flgg_hflx(0:glc_nec_max) = 0 !Subgrid heat flux
  integer, public :: index_g2x_Sg_tbot = 0 !bottom-layer temperature
  integer, public :: index_g2x_Sg_dztbot = 0 !distance from bottom-layer temperature to ice-ocean interface (to compute T gradient)
  integer, public :: index_g2x_Sg_lithop = 0 !ice sheet lithostatic pressure
  integer, public :: index_g2x_Sg_blis = 0 !Boundary layer interface salinity for ocean
  integer, public :: index_g2x_Sg_blit = 0 !Boundary layer interface temperature for ocean
  integer, public :: index_g2x_Fogx_qicelo = 0 !Subshelf heat flux for ocean
  integer, public :: index_g2x_Fogx_qiceho = 0 !Subshelf heat flux for the ice sheet  
  integer, public :: index_g2x_Sg_icemask = 0 !complete grounded/floating ice mask
  integer, public :: index_g2x_Sg_icemask_grounded = 0 !grounded mask
  integer, public :: index_g2x_Sg_icemask_floating = 0 !floating mask
  integer, public :: index_g2x_Sg_icemask_coupled_fluxes = 0 !ice mask for coupled fluxes
  
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

    !Following block are input boundary conditions for ice sheet model
    index_x2g_Flgl_qice =                 mct_avect_indexra(x2g,'Flgl_qice',perrwith='quiet')
    index_x2g_Sl_tsrf =                   mct_avect_indexra(x2g,'Sl_tsrf',perrwith='quiet')
    index_x2g_Fogx_qiceli =               mct_avect_indexra(x2g,'Fogx_qiceli',perrwith='quiet')
    index_x2g_Fogx_qicehi =               mct_avect_indexra(x2g,'Fogx_qicehi',perrwith='quiet')
    index_x2g_So_phieff =                 mct_avect_indexra(x2g,'So_phieff',perrwith='quiet')    

    !Following block of x2g/g2x vectors are used internally within coupler for subshelf melt flux
    !calculations (and so do not have directly-related export-side arrays)
    index_x2g_So_blt =                    mct_avect_indexra(x2g,'So_blt',perrwith='quiet')
    index_x2g_So_bls =                    mct_avect_indexra(x2g,'So_bls',perrwith='quiet')
    index_x2g_So_htv =                    mct_avect_indexra(x2g,'So_htv',perrwith='quiet')
    index_x2g_So_hsv =                    mct_avect_indexra(x2g,'So_hsv',perrwith='quiet')
    index_g2x_Sg_tbot =                   mct_avect_indexra(g2x,'Sg_tbot',perrwith='quiet')
    index_g2x_Sg_dztbot =                 mct_avect_indexra(g2x,'Sg_dztbot',perrwith='quiet')
    index_g2x_Sg_lithop =                 mct_avect_indexra(g2x,'Sg_lithop',perrwith='quiet')
    
    !Following block are outputs for other components.  Note that several fields
    !are calculated internally within the coupler, using the above information as input.
    index_g2x_Fogg_rofi =                 mct_avect_indexra(g2x,'Fogg_rofi',perrwith='quiet')
    index_g2x_Figg_rofi =                 mct_avect_indexra(g2x,'Figg_rofi',perrwith='quiet')
    index_g2x_Fogg_rofl =                 mct_avect_indexra(g2x,'Fogg_rofl',perrwith='quiet')
    index_g2x_Fogx_qicelo =               mct_avect_indexra(g2x,'Fogx_qicelo',perrwith='quiet')
    index_g2x_Fogx_qiceho =               mct_avect_indexra(g2x,'Fogx_qiceho',perrwith='quiet')
    index_g2x_Sg_blis =                   mct_avect_indexra(g2x,'Sg_blis',perrwith='quiet')
    index_g2x_Sg_blit =                   mct_avect_indexra(g2x,'Sg_blit',perrwith='quiet')    
    index_g2x_Sg_icemask =                mct_avect_indexra(g2x,'Sg_icemask',perrwith='quiet')
    index_g2x_Sg_icemask_grounded =       mct_avect_indexra(g2x,'Sg_icemask_grounded',perrwith='quiet')
    index_g2x_Sg_icemask_floating =       mct_avect_indexra(g2x,'Sg_icemask_floating',perrwith='quiet')
    index_g2x_Sg_icemask_coupled_fluxes = mct_avect_indexra(g2x,'Sg_icemask_coupled_fluxes',perrwith='quiet')
    
    do num = 0,glc_nec_max
       write(cnum,'(i2.2)') num
       name = 'Sg_frac' // cnum
       index_g2x_Sg_frac(num)   =         mct_avect_indexra(g2x,trim(name),perrwith='quiet') 
       name = 'Sg_topo' // cnum
       index_g2x_Sg_topo(num)   =         mct_avect_indexra(g2x,trim(name),perrwith='quiet')
       name = 'Flgg_hflx' // cnum
       index_g2x_Flgg_hflx(num) =         mct_avect_indexra(g2x,trim(name),perrwith='quiet')
       
       if ( index_g2x_Sg_frac(num)   == 0 .and. &
            index_g2x_Sg_topo(num)   == 0 .and. &
            index_g2x_Flgg_hflx(num) == 0 ) then
          exit
       end if
       glc_nec = num
    end do
    if (glc_nec == glc_nec_max) then
       write(6,*)'glc_cpl_indices error: glc_nec_max value has been reached ' 
       call shr_sys_abort ('glc_cpl_indices error: glc_nec_cpl cannot equal glc_nec_max')
    end if
    
    call mct_aVect_clean(x2g)
    call mct_aVect_clean(g2x)

  end subroutine glc_cpl_indices_set

end module glc_cpl_indices
