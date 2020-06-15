module radconstants

! This module contains constants that are specific to the radiative transfer
! code used in the RRTMG model.
!
! TODO: This all needs to be checked for RRTMGP implementation! Some of this
! might change; band mapping in the shortwave has definitely changed, and this
! module should reflect those changes.
!
! TODO: Should this data be handled in a more robust way? Much of this contains
! explicit mappings to indices, which would probably be better handled with get_
! functions. I.e., get_nswbands() could query the kdist objects in case of
! RRTMGP, and the diag indices could look up the actual bands used in the kdist
! objects as well. On that note, this module should probably go away if
! possible in the future, and we should provide more robust access to the
! radiation interface.

use shr_kind_mod,   only: r8 => shr_kind_r8
use cam_abortutils, only: endrun

implicit none
private
save

! SHORTWAVE DATA

! Number of spectral intervals
integer, parameter, public :: nswbands = 14
integer, parameter, public :: nlwbands = 16

! Wavenumbers of band boundaries
!
! Note: Currently rad_solar_var extends the lowest band down to
! 100 cm^-1 if it is too high to cover the far-IR. Any changes meant
! to affect IR solar variability should take note of this.
!real(r8), dimension(nswbands) :: wavenum_sw_lower, wavenum_sw_upper
real(r8) :: wavenum_sw_lower(nswbands) = & ! in cm^-1
  (/ 820._r8, 2680._r8, 3250._r8, 4000._r8, 4650._r8, 5150._r8, 6150._r8, &
    7700._r8, 8050._r8,12850._r8,16000._r8,22650._r8,29000._r8,38000._r8/)
real(r8) :: wavenum_sw_upper(nswbands) = & ! in cm^-1
  (/2680._r8, 3250._r8, 4000._r8, 4650._r8, 5150._r8, 6150._r8, 7700._r8, &
    8050._r8,12850._r8,16000._r8,22650._r8,29000._r8,38000._r8,50000._r8/)
real(r8) :: wavenum_lw_lower(nlwbands) = &! Longwave spectral band limits (cm-1)
    (/   10._r8,  250._r8, 500._r8,   630._r8,  700._r8,  820._r8,  980._r8, 1080._r8, &
       1180._r8, 1390._r8, 1480._r8, 1800._r8, 2080._r8, 2250._r8, 2390._r8, 2680._r8 /)
real(r8) :: wavenum_lw_upper(nlwbands) = &! Longwave spectral band limits (cm-1)
    (/  250._r8,  500._r8,  630._r8,  700._r8,  820._r8,  980._r8, 1080._r8, 1180._r8, &
       1390._r8, 1480._r8, 1800._r8, 2080._r8, 2250._r8, 2390._r8, 2680._r8, 3250._r8 /)

! Solar irradiance at 1 A.U. in W/m^2 assumed by radiation code
! Rescaled so that sum is precisely 1368.22 and fractional amounts sum to 1.0
! NOTE: This is copied from RRTMG; should be revisited for RRTMGP
real(r8), parameter :: solar_ref_band_irradiance(nswbands) = & 
   (/ &
    12.89_r8,  12.11_r8,  20.3600000000001_r8, 23.73_r8, &
    22.43_r8,  55.63_r8, 102.93_r8, 24.29_r8, &
   345.74_r8, 218.19_r8, 347.20_r8, &
   129.49_r8,  50.15_r8,   3.08_r8 &
   /)

! These are indices to the band for diagnostic output
integer, parameter, public :: idx_sw_diag = 11 ! index to sw visible band
integer, parameter, public :: idx_nir_diag = 9 ! index to sw near infrared (778-1240 nm) band
integer, parameter, public :: idx_uv_diag = 12 ! index to sw uv (345-441 nm) band
integer, parameter, public :: rrtmg_sw_cloudsim_band = 10  ! rrtmg band for .67 micron

! Number of evenly spaced intervals in rh
! The globality of this mesh may not be necessary
! Perhaps it could be specific to the aerosol
! But it is difficult to see how refined it must be
! for lookup.  This value was found to be sufficient
! for Sulfate and probably necessary to resolve the
! high variation near rh = 1.  Alternative methods
! were found to be too slow.
! Optimal approach would be for cam to specify size of aerosol
! based on each aerosol's characteristics.  Radiation 
! should know nothing about hygroscopic growth!
integer, parameter, public :: nrh = 1000  

! LONGWAVE DATA

! These are indices to the band for diagnostic output
integer, parameter, public :: idx_lw_diag = 7 ! index to (H20 window) LW band

integer, parameter, public :: rrtmg_lw_cloudsim_band = 6  ! rrtmg band for 10.5 micron

!These can go away when old camrt disappears
! Index of volc. abs., H2O non-window
integer, public, parameter :: idx_LW_H2O_NONWND=1
! Index of volc. abs., H2O window
integer, public, parameter :: idx_LW_H2O_WINDOW=2
! Index of volc. cnt. abs. 0500--0650 cm-1
integer, public, parameter :: idx_LW_0500_0650=3
! Index of volc. cnt. abs. 0650--0800 cm-1
integer, public, parameter :: idx_LW_0650_0800=4
! Index of volc. cnt. abs. 0800--1000 cm-1
integer, public, parameter :: idx_LW_0800_1000=5
! Index of volc. cnt. abs. 1000--1200 cm-1
integer, public, parameter :: idx_LW_1000_1200=6
! Index of volc. cnt. abs. 1200--2000 cm-1
integer, public, parameter :: idx_LW_1200_2000=7

! GASES TREATED BY RADIATION (line spectrae)

! gasses required by radiation
integer, public, parameter :: gasnamelength = 5
integer, public, parameter :: nradgas = 8
character(len=gasnamelength), public, parameter :: gaslist(nradgas) &
   = (/'H2O  ','O3   ', 'O2   ', 'CO2  ', 'N2O  ', 'CH4  ', 'CFC11', 'CFC12'/)

! what is the minimum mass mixing ratio that can be supported by radiation implementation?
real(r8), public, parameter :: minmmr(nradgas) &
   = epsilon(1._r8)

! Length of "optics type" string specified in optics files.
integer, parameter, public :: ot_length = 32

public :: rad_gas_index

public :: get_number_sw_bands, &
          set_sw_spectral_boundaries, &
          set_lw_spectral_boundaries, &
          get_sw_spectral_boundaries, &
          get_lw_spectral_boundaries, &
          get_sw_spectral_midpoints, &
          get_lw_spectral_midpoints, &
          get_ref_solar_band_irrad, &
          get_ref_total_solar_irrad, &
          get_solar_band_fraction_irrad, &
          get_band_index_sw, &
          get_band_index_lw, &
          test_get_band_index, &
          check_wavenumber_bounds

contains

!------------------------------------------------------------------------------

subroutine get_solar_band_fraction_irrad(fractional_irradiance)
   ! provide Solar Irradiance for each band in RRTMG

   ! fraction of solar irradiance in each band
   real(r8), intent(out) :: fractional_irradiance(1:nswbands)
   real(r8) :: tsi ! total solar irradiance

   tsi = sum(solar_ref_band_irradiance)
   fractional_irradiance = solar_ref_band_irradiance / tsi

end subroutine get_solar_band_fraction_irrad

!------------------------------------------------------------------------------

subroutine get_ref_total_solar_irrad(tsi)
   ! provide Total Solar Irradiance assumed by RRTMG

   real(r8), intent(out) :: tsi

   tsi = sum(solar_ref_band_irradiance)

end subroutine get_ref_total_solar_irrad

!------------------------------------------------------------------------------

subroutine get_ref_solar_band_irrad( band_irrad )

   ! solar irradiance in each band (W/m^2)
   real(r8), intent(out) :: band_irrad(nswbands)
 
   band_irrad = solar_ref_band_irradiance

end subroutine get_ref_solar_band_irrad

!------------------------------------------------------------------------------

subroutine get_number_sw_bands(number_of_bands)

   ! number of solar (shortwave) bands in the rrtmg code
   integer, intent(out) :: number_of_bands

   number_of_bands = nswbands

end subroutine get_number_sw_bands

!------------------------------------------------------------------------------

subroutine set_sw_spectral_boundaries(wavenum_bnds)
   real(r8), intent(in), dimension(:,:) :: wavenum_bnds
   integer :: ibnd
   do ibnd = 1,nswbands
      wavenum_sw_lower(ibnd) = wavenum_bnds(1,ibnd)
      wavenum_sw_upper(ibnd) = wavenum_bnds(2,ibnd)
   end do
end subroutine

!------------------------------------------------------------------------------

subroutine set_lw_spectral_boundaries(wavenum_bnds)
   real(r8), intent(in), dimension(:,:) :: wavenum_bnds
   integer :: ibnd
   do ibnd = 1,nlwbands
      wavenum_lw_lower(ibnd) = wavenum_bnds(1,ibnd)
      wavenum_lw_upper(ibnd) = wavenum_bnds(2,ibnd)
   end do
end subroutine

!------------------------------------------------------------------------------

subroutine get_lw_spectral_boundaries(lower_bounds, upper_bounds, units)
   ! provide spectral boundaries of each longwave band

   real(r8), intent(inout) :: lower_bounds(nlwbands), upper_bounds(nlwbands)
   character(*), intent(in) :: units ! requested units

   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      lower_bounds  = wavenum_lw_lower
      upper_bounds = wavenum_lw_upper
   case('m','meter','meters')
      lower_bounds  = 1.e-2_r8/wavenum_lw_upper
      upper_bounds = 1.e-2_r8/wavenum_lw_lower
   case('nm','nanometer','nanometers')
      lower_bounds  = 1.e7_r8/wavenum_lw_upper
      upper_bounds = 1.e7_r8/wavenum_lw_lower
   case('um','micrometer','micrometers','micron','microns')
      lower_bounds  = 1.e4_r8/wavenum_lw_upper
      upper_bounds = 1.e4_r8/wavenum_lw_lower
   case('cm','centimeter','centimeters')
      lower_bounds  = 1._r8/wavenum_lw_upper
      upper_bounds = 1._r8/wavenum_lw_lower
   case default
      call endrun('get_lw_spectral_boundaries: spectral units not acceptable'//units)
   end select

end subroutine get_lw_spectral_boundaries

!------------------------------------------------------------------------------

subroutine get_sw_spectral_boundaries(lower_bounds, upper_bounds, units)
   ! provide spectral boundaries of each shortwave band

   real(r8), intent(inout) :: lower_bounds(nswbands), upper_bounds(nswbands)
   character(*), intent(in) :: units ! requested units

   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      lower_bounds = wavenum_sw_lower
      upper_bounds = wavenum_sw_upper
   case('m','meter','meters')
      lower_bounds = 1.e-2_r8/wavenum_sw_upper
      upper_bounds = 1.e-2_r8/wavenum_sw_lower
   case('nm','nanometer','nanometers')
      lower_bounds = 1.e7_r8/wavenum_sw_upper
      upper_bounds = 1.e7_r8/wavenum_sw_lower
   case('um','micrometer','micrometers','micron','microns')
      lower_bounds = 1.e4_r8/wavenum_sw_upper
      upper_bounds = 1.e4_r8/wavenum_sw_lower
   case('cm','centimeter','centimeters')
      lower_bounds  = 1._r8/wavenum_sw_upper
      upper_bounds = 1._r8/wavenum_sw_lower
   case default
      call endrun('rad_constants.F90: spectral units not acceptable'//units)
   end select

end subroutine get_sw_spectral_boundaries

!------------------------------------------------------------------------------

subroutine get_sw_spectral_midpoints(band_midpoints, units)
   character(len=*), intent(in) :: units
   real(r8), intent(out) :: band_midpoints(nswbands)
   real(r8), dimension(nswbands) :: lower_bounds, upper_bounds
   integer :: i

   ! Get band limits
   call get_sw_spectral_boundaries(lower_bounds, upper_bounds, units)

   ! Compute band midpoints
   band_midpoints = 0._r8
   do i = 1,nswbands
      band_midpoints(i) = (lower_bounds(i) + upper_bounds(i)) / 2._r8
   end do
end subroutine get_sw_spectral_midpoints

!----------------------------------------------------------------------------

subroutine get_lw_spectral_midpoints(band_midpoints, units)
   character(len=*), intent(in) :: units
   real(r8), intent(out) :: band_midpoints(nlwbands)
   real(r8), dimension(nlwbands) :: lower_bounds, upper_bounds
   integer :: i

   ! Get band limits
   call get_lw_spectral_boundaries(lower_bounds, upper_bounds, units)

   ! Compute band midpoints
   band_midpoints = 0._r8
   do i = 1,nlwbands
      band_midpoints(i) = (lower_bounds(i) + upper_bounds(i)) / 2._r8
   end do
end subroutine get_lw_spectral_midpoints

!----------------------------------------------------------------------------

integer function get_band_index_sw(band, units) result(band_index)
   use assertions, only: assert
   real(r8), intent(in) :: band
   character(len=*), intent(in) :: units
   real(r8) :: band_wavenum
   integer :: idx
   character(len=128) :: err

   ! Spectral bands are in wavenumber, so look for band in wavenumber space
   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      band_wavenum = band
   case('m','meter','meters')
      band_wavenum = 1.e-2_r8 / band
   case('nm','nanometer','nanometers')
      band_wavenum = 1.e7_r8/band
   case('um','micrometer','micrometers','micron','microns')
      band_wavenum = 1.e4_r8/band
   case('cm','centimeter','centimeters')
      band_wavenum  = 1._r8/band
   case default
      call endrun('get_band_index_sw: spectral units not acceptable'//units)
   end select

   ! Look for band
   band_index = -1
   do idx = 1,size(wavenum_sw_lower)
      if (band_wavenum > wavenum_sw_lower(idx) .and. band_wavenum <= wavenum_sw_upper(idx)) then
         band_index = idx
         exit
      end if
   end do

   ! Make sure we found a valid band
   write(err,*) 'get_band_index_sw: index not found for', band, trim(units)
   call assert(band_index > 0, err)
end function get_band_index_sw

!-------------------------------------------------------------------------------

integer function get_band_index_lw(band, units) result(band_index)
   use assertions, only: assert
   real(r8), intent(in) :: band
   character(len=*), intent(in) :: units
   real(r8) :: band_wavenum
   integer :: idx
   character(len=128) :: err

   ! Spectral bands are in wavenumber, so look for band in wavenumber space
   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      band_wavenum = band
   case('m','meter','meters')
      band_wavenum = 1.e-2_r8 / band
   case('nm','nanometer','nanometers')
      band_wavenum = 1.e7_r8/band
   case('um','micrometer','micrometers','micron','microns')
      band_wavenum = 1.e4_r8/band
   case('cm','centimeter','centimeters')
      band_wavenum  = 1._r8/band
   case default
      call endrun('get_band_index_lw: spectral units not acceptable'//units)
   end select

   ! Look for band
   band_index = -1
   do idx = 1,size(wavenum_lw_lower)
      if (band_wavenum > wavenum_lw_lower(idx) .and. band_wavenum <= wavenum_lw_upper(idx)) then
         band_index = idx
         exit
      end if
   end do

   ! Make sure we found a valid band
   write(err,*) 'get_band_index_lw: index not found for', band, trim(units)
   call assert(band_index > 0, err)
end function get_band_index_lw

!-------------------------------------------------------------------------------

! Unit test for get_band_index function
subroutine test_get_band_index()
   use assertions, only: assert
   call assert(get_band_index_sw(0.67_r8, 'micron') == 10, 'get_band_index_sw(0.67 micron) /= 10')
   call assert(get_band_index_lw(10.5_r8, 'micron') == 6 , 'get_band_index_lw(10.5 micron) /= 6' )
end subroutine

!------------------------------------------------------------------------------

integer function rad_gas_index(gasname)

   ! return the index in the gaslist array of the specified gasname

   character(len=*),intent(in) :: gasname
   integer :: igas

   rad_gas_index = -1
   do igas = 1, nradgas
      if (trim(gaslist(igas)).eq.trim(gasname)) then
         rad_gas_index = igas
         return
      endif
   enddo
   call endrun ("rad_gas_index: can not find gas with name "//gasname)
end function rad_gas_index

!------------------------------------------------------------------------------
! Testing routines below here

! Simple check to make sure we are setting the wavenumber bounds correctly. This does not need to
! be run during normal execution, and is included here just to make sure that our setting of
! wavenumber bounds is consistent with what we expect. In general, we are not hard-coding these
! values anymore because with RRTMGP, we allow changing these via the absorption coefficient
! inputfiles.
subroutine check_wavenumber_bounds()
   use assertions, only: assert
   real(r8) :: wavenum_sw_lower_old(nswbands) = & ! in cm^-1
     (/ 820._r8, 2680._r8, 3250._r8, 4000._r8, 4650._r8, 5150._r8, 6150._r8, &
       7700._r8, 8050._r8,12850._r8,16000._r8,22650._r8,29000._r8,38000._r8/)
   real(r8) :: wavenum_sw_upper_old(nswbands) = & ! in cm^-1
     (/2680._r8, 3250._r8, 4000._r8, 4650._r8, 5150._r8, 6150._r8, 7700._r8, &
       8050._r8,12850._r8,16000._r8,22650._r8,29000._r8,38000._r8,50000._r8/)
   real(r8) :: wavenum_lw_lower_old(nlwbands) = &! Longwave spectral band limits (cm-1)
       (/   10._r8,  250._r8, 500._r8,   630._r8,  700._r8,  820._r8,  980._r8, 1080._r8, &
          1180._r8, 1390._r8, 1480._r8, 1800._r8, 2080._r8, 2250._r8, 2390._r8, 2680._r8 /)
   real(r8) :: wavenum_lw_upper_old(nlwbands) = &! Longwave spectral band limits (cm-1)
       (/  250._r8,  500._r8,  630._r8,  700._r8,  820._r8,  980._r8, 1080._r8, 1180._r8, &
          1390._r8, 1480._r8, 1800._r8, 2080._r8, 2250._r8, 2390._r8, 2680._r8, 3250._r8 /)
   call assert(all(wavenum_sw_lower .eq. wavenum_sw_lower_old), 'wavenum_sw_lower')
   call assert(all(wavenum_sw_upper .eq. wavenum_sw_upper_old), 'wavenum_sw_upper')
   call assert(all(wavenum_lw_lower .eq. wavenum_lw_lower_old), 'wavenum_lw_lower')
   call assert(all(wavenum_sw_upper .eq. wavenum_sw_upper_old), 'wavenum_sw_upper')
end subroutine check_wavenumber_bounds

!------------------------------------------------------------------------------

end module radconstants
