module glc_cpl_indices
  
  use seq_flds_mod
  use mct_mod
  use glc_constants, only : glc_nec, glc_smb
  use shr_sys_mod  , only : shr_sys_abort

  implicit none

  SAVE
  public

  integer , parameter, private:: glc_nec_max = 100

  ! Note that, in both the drv -> glc and the glc -> drv fields, index 0 means bare land

  ! drv -> glc

  integer, public :: index_x2g_Sl_tsrf(0:glc_nec_max)   = 0
  integer, public :: index_x2g_Sl_topo(0:glc_nec_max)   = 0
  integer, public :: index_x2g_Flgl_qice(0:glc_nec_max) = 0

  ! glc -> drv

  integer, public :: index_g2x_Fogg_rofi = 0   ! frozen runoff -> ocn
  integer, public :: index_g2x_Figg_rofi = 0   ! frozen runoff -> ice
  integer, public :: index_g2x_Fogg_rofl = 0   ! liquid runoff -> ocn
  integer, public :: index_g2x_Sg_frac(0:glc_nec_max)   = 0
  integer, public :: index_g2x_Sg_topo(0:glc_nec_max)   = 0
  integer, public :: index_g2x_Flgg_hflx(0:glc_nec_max) = 0
  integer, public :: index_g2x_Sg_icemask = 0
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

    glc_nec = 0

    ! glc -> drv

    index_g2x_Fogg_rofi = mct_avect_indexra(g2x,'Fogg_rofi',perrwith='quiet')
    index_g2x_Figg_rofi = mct_avect_indexra(g2x,'Figg_rofi',perrwith='quiet')
    index_g2x_Fogg_rofl = mct_avect_indexra(g2x,'Fogg_rofl',perrwith='quiet')

    do num = 0,glc_nec_max
       write(cnum,'(i2.2)') num
       name = 'Sg_frac' // cnum
       index_g2x_Sg_frac(num)   = mct_avect_indexra(g2x,trim(name),perrwith='quiet') 
       name = 'Sg_topo' // cnum
       index_g2x_Sg_topo(num)   = mct_avect_indexra(g2x,trim(name),perrwith='quiet')
       name = 'Flgg_hflx' // cnum
       index_g2x_Flgg_hflx(num) = mct_avect_indexra(g2x,trim(name),perrwith='quiet')
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

    index_g2x_Sg_icemask = mct_avect_indexra(g2x,'Sg_icemask',perrwith='quiet')
    index_g2x_Sg_icemask_coupled_fluxes = mct_avect_indexra(g2x,'Sg_icemask_coupled_fluxes',perrwith='quiet')

    ! drv -> glc

    do num = 0,glc_nec
       write(cnum,'(i2.2)') num
       name = 'Sl_tsrf' // cnum
       index_x2g_Sl_tsrf(num)   = mct_avect_indexra(x2g,trim(name))
       name = 'Sl_topo' // cnum
       index_x2g_Sl_topo(num)   = mct_avect_indexra(x2g,trim(name))
       name = 'Flgl_qice' // cnum
       index_x2g_Flgl_qice(num) = mct_avect_indexra(x2g,trim(name))
    end do

    call mct_aVect_clean(x2g)
    call mct_aVect_clean(g2x)

    ! Set glc_smb
    ! true => get surface mass balance from CLM via coupler (in multiple elev classes)
    ! false => use PDD scheme in GLIMMER

    if (glc_nec > 0) then
       glc_smb = .true.
    else
       glc_smb = .false.
    end if

  end subroutine glc_cpl_indices_set

end module glc_cpl_indices
