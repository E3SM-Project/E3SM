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

  use cam_history,  only: phys_decomp, dyn_decomp, addfld, add_default
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
    call addfld ('SGH     ','m       ',1,    'I','Standard deviation of orography',phys_decomp)
    call addfld ('SGH30   ','m       ',1,    'I','Standard deviation of 30s orography',phys_decomp)


!jt
!jt Maybe add this to scam specific initialization
!jt

#if ( defined BFB_CAM_SCAM_IOP )
    call addfld ('CLAT1&IC','none      ',1,    'I','cos lat for bfb testing', dyn_decomp)
    call add_default ('CLAT1&IC       ',0, 'I')
    call addfld ('CLON1&IC','none      ',1,    'I','cos lon for bfb testing', dyn_decomp)
    call add_default ('CLON1&IC       ',0, 'I')
    call addfld ('PHI&IC','none      ',1,    'I','lat for bfb testing', dyn_decomp)
    call add_default ('PHI&IC       ',0, 'I')
    call addfld ('LAM&IC','none      ',1,    'I','lon for bfb testing', dyn_decomp)
    call add_default ('LAM&IC       ',0, 'I')
#endif

    call addfld ('DQP     ','kg/kg/s ',pver, 'A','Specific humidity tendency due to precipitation',phys_decomp)

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
    call addfld ('CLAT ','none      ',1,    'A','cos lat for bfb testing', dyn_decomp)
    call add_default ('CLAT',2,' ')
    call addfld ('q','kg/kg   ',plev, 'A','Q for scam',dyn_decomp)
    call add_default ('q',2, ' ')
    call addfld ('u','m/s     ',plev, 'A','U for scam',dyn_decomp)
    call add_default ('u',2,' ')
    call addfld ('v','m/s     ',plev, 'A','V for scam',dyn_decomp)
    call add_default ('v',2,' ')
    call addfld ('t','K       ',plev, 'A','Temperature for scam',dyn_decomp)
    call add_default ('t',2,' ')
    call addfld ('Tg','K      ',1,    'A','Surface temperature (radiative) for scam',phys_decomp)
    call add_default ('Tg',2,' ')
    call addfld ('Ps','Pa      ',1, 'A','Ps for scam',dyn_decomp)
    call add_default ('Ps',2,' ')
    call addfld ('divT3d','K       ',plev, 'A','Dynamics Residual for T',dyn_decomp)
    call add_default ('divT3d',2,' ')
    call addfld ('fixmas','percent',1, 'A','Mass fixer',dyn_decomp)
    call add_default ('fixmas',2,' ')
    call addfld ('beta','percent  ',1, 'A','Mass fixer',dyn_decomp)
    call add_default ('beta',2,' ')
    do m=1,pcnst
       call addfld (trim(cnst_name(m))//'_dten','kg/kg   ',plev, 'A', &
            trim(cnst_name(m))//' IOP Dynamics Residual for '//trim(cnst_name(m)),dyn_decomp)
       call add_default (trim(cnst_name(m))//'_dten',2,' ')
       call addfld (trim(cnst_name(m))//'_alph','kg/kg   ',1, 'A',trim(cnst_name(m))//' alpha constituent fixer',dyn_decomp)
       call add_default (trim(cnst_name(m))//'_alph',2,' ')
       call addfld (trim(cnst_name(m))//'_dqfx','kg/kg   ',plev, 'A',trim(cnst_name(m))//' dqfx3 fixer',dyn_decomp)
       call add_default (trim(cnst_name(m))//'_dqfx',2,' ')
    end do
    call addfld ('shflx ','W/m2    ',1,    'A','Surface sensible heat flux for scam',phys_decomp)
    call add_default ('shflx ',2,' ')
    call addfld ('lhflx   ','W/m2    ',1,    'A','Surface latent heat flux for scam',phys_decomp)
    call add_default ('lhflx   ',2,' ')
    call addfld ('trefht  ','K       ',1,    'A','Reference height temperature',phys_decomp)
    call add_default ('trefht  ',2,' ')
    call addfld ('Tsair  ','K       ',1,    'A','Reference height temperature for scam',phys_decomp)
    call add_default ('Tsair  ',2,' ')
    call addfld ('phis   ','m2/s2   ',1,    'I','Surface geopotential for scam',phys_decomp)
    call add_default ('phis   ',2,' ')
    call addfld ('Prec   ','m/s     ',1,    'A','Total (convective and large-scale) precipitation rate for scam',phys_decomp)
    call add_default ('Prec   ',2,' ')
    call addfld ('omega   ','Pa/s    ',pver, 'A','Vertical velocity (pressure)',phys_decomp)
    call add_default ('omega   ',2,' ')

  end subroutine initialize_iop_history
#endif

end module history_defaults
