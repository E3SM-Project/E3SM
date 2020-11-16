!-------------------------------------------------------------------------------
! This module uses the Lean solar irradiance data to provide a solar cycle
! scaling factor used in heating rate calculations 
!
! NOTE: this module seems to have been created for RRTMG to update the solar
! irradiance used in the flux calculations without directly modifying the data
! hard-coded into the RRTMG code. This is no longer needed for RRTMGP, but we
! keep this module for compatibility with the calls in physpkg.
!-------------------------------------------------------------------------------
module rad_solar_var

  use shr_kind_mod , only : r8 => shr_kind_r8
  use solar_data,    only : sol_irrad, we, nbins, has_spectrum, sol_tsi
  use solar_data,    only : do_spctrl_scaling
  use cam_abortutils,    only : endrun

  implicit none
  save

  private
  public :: rad_solar_var_init
  public :: get_variability

  real(r8), allocatable :: ref_band_irrad(:)  ! scaling will be relative to ref_band_irrad in each band
  real(r8), allocatable :: irrad(:)           ! solar irradiance at model timestep in each band
  real(r8)              :: tsi_ref            ! total solar irradiance assumed by rrtmg                                                 

  real(r8), allocatable :: radbinmax(:)
  real(r8), allocatable :: radbinmin(:)
  integer :: nradbins
contains

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine rad_solar_var_init( )
    use radconstants,  only : get_number_sw_bands
    use radconstants,  only : get_sw_spectral_boundaries
    use radconstants,  only : get_ref_solar_band_irrad
    use radconstants,  only : get_ref_total_solar_irrad

    integer :: i
    integer :: ierr
    integer :: radmax_loc


    call get_number_sw_bands(nradbins)

    if ( do_spctrl_scaling ) then

       if ( .not.has_spectrum ) then
          call endrun('rad_solar_var_init: solar input file must have irradiance spectrum')
       endif

       allocate (radbinmax(nradbins),stat=ierr)
       if (ierr /= 0) then
          call endrun('rad_solar_var_init: Error allocating space for radbinmax')
       end if

       allocate (radbinmin(nradbins),stat=ierr)
       if (ierr /= 0) then
          call endrun('rad_solar_var_init: Error allocating space for radbinmin')
       end if

       allocate (ref_band_irrad(nradbins), stat=ierr)
       if (ierr /= 0) then
          call endrun('rad_solar_var_init: Error allocating space for ref_band_irrad')
       end if

       allocate (irrad(nradbins), stat=ierr)
       if (ierr /= 0) then
          call endrun('rad_solar_var_init: Error allocating space for irrad')
       end if

       call get_sw_spectral_boundaries(radbinmin, radbinmax, 'nm')

       ! Make sure that the far-IR is included, even if RRTMG does not
       ! extend that far down. 10^5 nm corresponds to a wavenumber of
       ! 100 cm^-1.
       radmax_loc = maxloc(radbinmax,1)
       radbinmax(radmax_loc) = max(100000._r8,radbinmax(radmax_loc))

       ! for rrtmg, reference spectrum from rrtmg
       call get_ref_solar_band_irrad( ref_band_irrad )

    else

       call get_ref_total_solar_irrad(tsi_ref)

    endif

  end subroutine rad_solar_var_init

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine get_variability( sfac )

    real(r8), intent(out) :: sfac(nradbins)       ! scaling factors for CAM heating

    if ( do_spctrl_scaling ) then

      call integrate_spectrum( nbins, nradbins, we, radbinmin, radbinmax, sol_irrad, irrad)

      sfac(:nradbins) = irrad(:nradbins)/ref_band_irrad(:nradbins)

    else

       sfac(:nradbins) = sol_tsi/tsi_ref

    endif

  endsubroutine get_variability

!-------------------------------------------------------------------------------
! private method.........
!-------------------------------------------------------------------------------

  subroutine integrate_spectrum( nsrc, ntrg, src_x, min_trg, max_trg, src, trg )

    use mo_util, only : rebin

    implicit none

    !---------------------------------------------------------------
    !	... dummy arguments
    !---------------------------------------------------------------
    integer,  intent(in)  :: nsrc                  ! dimension source array
    integer,  intent(in)  :: ntrg                  ! dimension target array
    real(r8), intent(in)  :: src_x(nsrc+1)         ! source coordinates
    real(r8), intent(in)  :: max_trg(ntrg)         ! target coordinates
    real(r8), intent(in)  :: min_trg(ntrg)         ! target coordinates
    real(r8), intent(in)  :: src(nsrc)             ! source array
    real(r8), intent(out) :: trg(ntrg)             ! target array
 
    !---------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------
    real(r8) :: trg_x(2), targ(1)         ! target coordinates
    integer  :: i

    do i = 1, ntrg

       trg_x(1) = min_trg(i)
       trg_x(2) = max_trg(i)

       call rebin( nsrc, 1, src_x, trg_x, src(1:nsrc), targ(:) )
       ! W/m2/nm --> W/m2
       trg( i ) = targ(1)*(trg_x(2)-trg_x(1))

    enddo


  end subroutine integrate_spectrum

endmodule rad_solar_var
