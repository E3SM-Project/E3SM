module history_defaults
!----------------------------------------------------------------------- 
! 
! Purpose: contains calls to setup default history stuff that has not found
!          a proper home yet. Shouldn't really exist.
!
! Public functions/subroutines:
!   bldfld
! 
! Author: B.A. Boville from code in cam_history.F90
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
  use constituents, only: pcnst, cnst_name
  use ppgrid,       only: pver, pverp
  use pmgrid,       only: plev, plevp
  use dycore,       only: dycore_is

  use cam_history,  only:   addfld, horiz_only, add_default
  implicit none

  PRIVATE

  public :: bldfld

#if ( defined BFB_CAM_SCAM_IOP )
  public :: initialize_iop_history
#endif

CONTAINS


!#######################################################################
  subroutine bldfld ()
!
!----------------------------------------------------------------------- 
! 
! Purpose: 
!
! Build Master Field List of all possible fields in a history file.  Each field has 
! associated with it a "long_name" netcdf attribute that describes what the field is, 
! and a "units" attribute.
! 
! Method: Call a subroutine to add each field
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Local workspace
!
    integer m                     ! Index

!
! Call addfld to add each field to the Master Field List.
!
    call addfld ('SGH',horiz_only,    'I','m','Standard deviation of orography')
    call addfld ('SGH30',horiz_only,    'I','m','Standard deviation of 30s orography')


!jt
!jt Maybe add this to scam specific initialization
!jt

#if ( defined BFB_CAM_SCAM_IOP )
    call addfld ('CLAT1&IC',horiz_only,    'I','none','cos lat for bfb testing', gridname=gauss_grid)
    call add_default ('CLAT1&IC       ',0, 'I')
    call addfld ('CLON1&IC',horiz_only,    'I','none','cos lon for bfb testing', gridname=gauss_grid)
    call add_default ('CLON1&IC       ',0, 'I')
    call addfld ('PHI&IC',horiz_only,    'I','none','lat for bfb testing', gridname=gauss_grid)
    call add_default ('PHI&IC       ',0, 'I')
    call addfld ('LAM&IC',horiz_only,    'I','none','lon for bfb testing', gridname=gauss_grid)
    call add_default ('LAM&IC       ',0, 'I')
#endif

    call addfld ('DQP',(/ 'lev' /), 'A','kg/kg/s','Specific humidity tendency due to precipitation')

  end subroutine bldfld

!#######################################################################
#if ( defined BFB_CAM_SCAM_IOP  )
  subroutine initialize_iop_history()
!
! !DESCRIPTION: 
! !USES:
    use iop
    use ppgrid, only: pver, pverp
    use phys_control,     only: phys_getopts
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer m
!-----------------------------------------------------------------------
    call addfld ('CLAT',horiz_only,    'A','none','cos lat for bfb testing', gridname=gauss_grid)
    call add_default ('CLAT',2,' ')
    call addfld ('q',(/ 'lev' /), 'A','kg/kg','Q for scam',gridname=gauss_grid)
    call add_default ('q',2, ' ')
    call addfld ('u',(/ 'lev' /), 'A','m/s','U for scam',gridname=gauss_grid)
    call add_default ('u',2,' ')
    call addfld ('v',(/ 'lev' /), 'A','m/s','V for scam',gridname=gauss_grid)
    call add_default ('v',2,' ')
    call addfld ('t',(/ 'lev' /), 'A','K','Temperature for scam',gridname=gauss_grid)
    call add_default ('t',2,' ')
    call addfld ('Tg',horiz_only,    'A','K','Surface temperature (radiative) for scam')
    call add_default ('Tg',2,' ')
    call addfld ('Ps',horiz_only, 'A','Pa','Ps for scam',gridname=gauss_grid)
    call add_default ('Ps',2,' ')
    call addfld ('divT3d',(/ 'lev' /), 'A','K','Dynamics Residual for T',gridname=gauss_grid)
    call add_default ('divT3d',2,' ')
    call addfld ('fixmas',horiz_only, 'A','percent','Mass fixer',gridname=gauss_grid)
    call add_default ('fixmas',2,' ')
    call addfld ('beta',horiz_only, 'A','percent','Mass fixer',gridname=gauss_grid)
    call add_default ('beta',2,' ')
    do m=1,pcnst
       call addfld (trim(cnst_name(m))//'_dten',(/ 'lev' /), 'A','kg/kg', &
            trim(cnst_name(m))//' IOP Dynamics Residual for '//trim(cnst_name(m)),gridname=gauss_grid)
       call add_default (trim(cnst_name(m))//'_dten',2,' ')
       call addfld (trim(cnst_name(m))//'_alph',horiz_only, 'A','kg/kg',trim(cnst_name(m))//' alpha constituent fixer',gridname=gauss_grid)
       call add_default (trim(cnst_name(m))//'_alph',2,' ')
       call addfld (trim(cnst_name(m))//'_dqfx',(/ 'lev' /), 'A','kg/kg',trim(cnst_name(m))//' dqfx3 fixer',gridname=gauss_grid)
       call add_default (trim(cnst_name(m))//'_dqfx',2,' ')
    end do
    call addfld ('shflx',horiz_only,    'A','W/m2','Surface sensible heat flux for scam')
    call add_default ('shflx ',2,' ')
    call addfld ('lhflx',horiz_only,    'A','W/m2','Surface latent heat flux for scam')
    call add_default ('lhflx   ',2,' ')
    call addfld ('trefht',horiz_only,    'A','K','Reference height temperature')
    call add_default ('trefht  ',2,' ')
    call addfld ('Tsair',horiz_only,    'A','K','Reference height temperature for scam')
    call add_default ('Tsair  ',2,' ')
    call addfld ('phis',horiz_only,    'I','m2/s2','Surface geopotential for scam')
    call add_default ('phis   ',2,' ')
    call addfld ('Prec',horiz_only,    'A','m/s','Total (convective and large-scale) precipitation rate for scam')
    call add_default ('Prec   ',2,' ')
    call addfld ('omega',(/ 'lev' /), 'A','Pa/s','Vertical velocity (pressure)')
    call add_default ('omega   ',2,' ')

  end subroutine initialize_iop_history
#endif

end module history_defaults
