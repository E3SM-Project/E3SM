module tidal_diag 

  !---------------------------------------------------------------------------------
  ! Module to compute fourier coefficients for the diurnal and semidiurnal tide 
  !
  ! Created by: Dan Marsh
  ! Date: 12 May 2008
  !---------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8 => shr_kind_r8
  use ppgrid,        only: pcols, pver

  implicit none

  private

  ! Public interfaces

  public :: tidal_diag_init   ! create coefficient history file variables
  public :: tidal_diag_write  ! calculate and output dignostics
  public :: get_tidal_coeffs

contains

  !===============================================================================

  subroutine  tidal_diag_init()
    !----------------------------------------------------------------------- 
    ! Purpose: create fourier coefficient history file variables
    !-----------------------------------------------------------------------

    use cam_history,        only: addfld, phys_decomp

    call addfld ('T_24_COS','K       ',pver, 'A','Temperature 24hr. cos coeff.',phys_decomp)
    call addfld ('T_24_SIN','K       ',pver, 'A','Temperature 24hr. sin coeff.',phys_decomp)
    call addfld ('T_12_COS','K       ',pver, 'A','Temperature 12hr. cos coeff.',phys_decomp)
    call addfld ('T_12_SIN','K       ',pver, 'A','Temperature 12hr. sin coeff.',phys_decomp)

    call addfld ('U_24_COS','m/s     ',pver, 'A','Zonal wind 24hr. cos coeff.',phys_decomp)
    call addfld ('U_24_SIN','m/s     ',pver, 'A','Zonal wind 24hr. sin coeff.',phys_decomp)
    call addfld ('U_12_COS','m/s     ',pver, 'A','Zonal wind 12hr. cos coeff.',phys_decomp)
    call addfld ('U_12_SIN','m/s     ',pver, 'A','Zonal wind 12hr. sin coeff.',phys_decomp)

    call addfld ('V_24_COS','m/s     ',pver, 'A','Meridional wind 24hr. cos coeff.',phys_decomp)
    call addfld ('V_24_SIN','m/s     ',pver, 'A','Meridional wind 24hr. sin coeff.',phys_decomp)
    call addfld ('V_12_COS','m/s     ',pver, 'A','Meridional wind 12hr. cos coeff.',phys_decomp)
    call addfld ('V_12_SIN','m/s     ',pver, 'A','Meridional wind 12hr. sin coeff.',phys_decomp)

    call addfld ('PS_24_COS','Pa     ',1,    'A','surface pressure 24hr. cos coeff.',phys_decomp)
    call addfld ('PS_24_SIN','Pa     ',1,    'A','surface pressure 24hr. sin coeff.',phys_decomp)
    call addfld ('PS_12_COS','Pa     ',1,    'A','surface pressure 12hr. cos coeff.',phys_decomp)
    call addfld ('PS_12_SIN','Pa     ',1,    'A','surface pressure 12hr. sin coeff.',phys_decomp)

    call addfld ('OMEGA_24_COS','Pa/s',pver, 'A','vertical pressure velocity 24hr. cos coeff.',phys_decomp)
    call addfld ('OMEGA_24_SIN','Pa/s',pver, 'A','vertical pressure velocity 24hr. sin coeff.',phys_decomp)
    call addfld ('OMEGA_12_COS','Pa/s',pver, 'A','vertical pressure velocity 12hr. cos coeff.',phys_decomp)
    call addfld ('OMEGA_12_SIN','Pa/s',pver, 'A','vertical pressure velocity 12hr. sin coeff.',phys_decomp)

    return

  end subroutine tidal_diag_init

  !===============================================================================

  subroutine  tidal_diag_write(state)

    !----------------------------------------------------------------------- 
    ! Purpose: calculate fourier coefficients and save to history files 
    !-----------------------------------------------------------------------
    use cam_history,   only: outfld, hist_fld_active
    use physics_types, only: physics_state

    implicit none

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(physics_state), intent(in) :: state
    !
    !---------------------------Local workspace-----------------------------

    integer  :: lchnk

    real(r8) :: dcoef(4) 
    integer :: ncol

    !-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol = state%ncol

    call get_tidal_coeffs( dcoef )

    if ( hist_fld_active('T_24_COS') .or. hist_fld_active('T_24_SIN') ) then
       call outfld( 'T_24_SIN', state%t(:ncol,:)*dcoef(1), ncol, lchnk )
       call outfld( 'T_24_COS', state%t(:ncol,:)*dcoef(2), ncol, lchnk )
    endif
    if ( hist_fld_active('T_12_COS') .or. hist_fld_active('T_12_SIN') ) then
       call outfld( 'T_12_SIN', state%t(:ncol,:)*dcoef(3), ncol, lchnk )
       call outfld( 'T_12_COS', state%t(:ncol,:)*dcoef(4), ncol, lchnk )
    endif

    if ( hist_fld_active('U_24_COS') .or. hist_fld_active('U_24_SIN') ) then
       call outfld( 'U_24_SIN', state%u(:ncol,:)*dcoef(1), ncol, lchnk )
       call outfld( 'U_24_COS', state%u(:ncol,:)*dcoef(2), ncol, lchnk )
    endif
    if ( hist_fld_active('U_12_COS') .or. hist_fld_active('U_12_SIN') ) then
       call outfld( 'U_12_SIN', state%u(:ncol,:)*dcoef(3), ncol, lchnk )
       call outfld( 'U_12_COS', state%u(:ncol,:)*dcoef(4), ncol, lchnk )
    endif

    if ( hist_fld_active('V_24_COS') .or. hist_fld_active('V_24_SIN') ) then
       call outfld( 'V_24_SIN', state%v(:ncol,:)*dcoef(1), ncol, lchnk )
       call outfld( 'V_24_COS', state%v(:ncol,:)*dcoef(2), ncol, lchnk )
    endif
    if ( hist_fld_active('V_12_COS') .or. hist_fld_active('V_12_SIN') ) then
       call outfld( 'V_12_SIN', state%v(:ncol,:)*dcoef(3), ncol, lchnk )
       call outfld( 'V_12_COS', state%v(:ncol,:)*dcoef(4), ncol, lchnk )
    endif

    if ( hist_fld_active('PS_24_COS') .or. hist_fld_active('PS_24_SIN') ) then
       call outfld( 'PS_24_SIN', state%ps(:ncol)*dcoef(1), ncol, lchnk )
       call outfld( 'PS_24_COS', state%ps(:ncol)*dcoef(2), ncol, lchnk )
    endif
    if ( hist_fld_active('PS_12_COS') .or. hist_fld_active('PS_12_SIN') ) then
       call outfld( 'PS_12_SIN', state%ps(:ncol)*dcoef(3), ncol, lchnk )
       call outfld( 'PS_12_COS', state%ps(:ncol)*dcoef(4), ncol, lchnk )
    endif

    if ( hist_fld_active('OMEGA_24_COS') .or. hist_fld_active('OMEGA_24_SIN') ) then
       call outfld( 'OMEGA_24_SIN', state%omega(:ncol,:)*dcoef(1), ncol, lchnk )
       call outfld( 'OMEGA_24_COS', state%omega(:ncol,:)*dcoef(2), ncol, lchnk )
    endif
    if ( hist_fld_active('OMEGA_12_COS') .or. hist_fld_active('OMEGA_12_SIN') ) then
       call outfld( 'OMEGA_12_SIN', state%omega(:ncol,:)*dcoef(3), ncol, lchnk )
       call outfld( 'OMEGA_12_COS', state%omega(:ncol,:)*dcoef(4), ncol, lchnk )
    endif

    return

  end subroutine tidal_diag_write

  !===============================================================================
  subroutine get_tidal_coeffs( dcoef )

    !----------------------------------------------------------------------- 
    ! Purpose: calculate fourier coefficients
    !-----------------------------------------------------------------------

    use time_manager,  only: get_curr_date               
    use physconst, only: pi, cday

    real(r8), intent(out) :: dcoef(4) 

 !  variables to calculate tidal coeffs
    real(r8), parameter :: pi_x_2 = 2._r8*pi
    real(r8), parameter :: pi_x_4 = 4._r8*pi
    integer  :: year, month
    integer  :: day              ! day of month
    integer  :: tod              ! time of day (seconds past 0Z) 
    real(r8) :: gmtfrac 

 !  calculate multipliers for Fourier transform in time (tidal analysis)
    call get_curr_date(year, month, day, tod)
    gmtfrac = tod / cday

    dcoef(1) = 2._r8*sin(pi_x_2*gmtfrac)
    dcoef(2) = 2._r8*cos(pi_x_2*gmtfrac)
    dcoef(3) = 2._r8*sin(pi_x_4*gmtfrac)
    dcoef(4) = 2._r8*cos(pi_x_4*gmtfrac)

  end subroutine get_tidal_coeffs

end module tidal_diag

