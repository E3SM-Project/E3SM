module rad
  use parkind, only : kind_rb ! RRTM expects reals with this kind parameter 
                                   ! (8 byte reals) 
  implicit none

  ! This module holds a bunch of stuff used in the GCSS-CFMIP LES
  !   interface to the RRTM radiation.

  integer, save :: nzrad, nzpatch, npatch_start, npatch_end
  integer, save :: nzrad_old = -1

  real, save :: day_when_patch_tracegases_last_updated = -10.

  ! background sounding, read in from SCAM forcing file.
  integer, parameter :: nzsnd = 120
  real :: psnd(nzsnd) ! pressure sounding read in from SoundingFileName, mb
  real :: tsnd(nzsnd) ! temperature sounding read in from SoundingFileName, K
  real :: qsnd(nzsnd) ! water vapor sounding read in from SoundingFileName, kg/kg

  ! Radiative heating rate (K/s) on model domain  
  real, dimension(:,:,:), allocatable, save :: qrad

  ! Shortwave and longwave heating rates (K/s), fluxes, water paths, and 
  ! effective radii on model domain
  real, dimension(:,:,:), allocatable :: &
       lwUp_3d, lwDown_3d, &
       swUp_3d, swDown_3d, &
       lwHeatingRate_3d, swHeatingRate_3d, &
       lwp_3d, iwp_3d, &
       reliq_3d, reice_3d

  ! surface and top-of-atmosphere (TOA) radiative fluxes
  real, dimension(:,:), allocatable, save :: &
       NetlwUpSurface, &
       NetlwUpSurfaceClearSky, &
       NetlwUpToa, &
       NetlwUpToaClearSky, &
       NetswDownSurface, &
       NetswDownSurfaceClearSky, &
       NetswDownToa, &
       NetswDownToaClearSky, &
       NetswUpToa, &
       insolation_TOA, &
       swDownSurface, &
       lwDownSurface, &
       swnsxy, &
       lwnsxy

  ! input CRM/LES fields to radiation
  logical :: isAllocated_RadInputsOutputs = .false.

  real, dimension(:,:), allocatable :: &
       tabs_slice, & ! absolute temperature, K
       qv_slice, & ! water vapor mass mixing ratio, kg/kg
       qcl_slice, & ! cloud liquid mass mixing ratio, kg/kg
       qci_slice ! cloud ice mass mixing ratio, kg/kg

  real, dimension(:), allocatable :: &
       tg_slice, & ! surface temperature, K
       pres_input, &
       presi_input, &
       latitude_slice, &
       longitude_slice

! Profiles of trace gases
  real, dimension(:), allocatable :: &
       o3_slice, &
       co2_slice, &
       ch4_slice, &
       n2o_slice, &
       o2_slice, &
       cfc11_slice, &
       cfc12_slice, &
       cfc22_slice, &
       ccl4_slice

! Perpetual factor and perpetual cosine solar zenith angle
  real, dimension(:), allocatable :: &
       p_factor_slice, &
       p_coszrs_slice

! fluxes output from RRTM radiation scheme.
  real(kind=kind_rb), dimension(:,:), allocatable :: &
       lwUp , &           ! upward longwave radiative flux (W/m2)
       lwDown , &         ! downward longwave radiative flux (W/m2)
       lwUpClearSky , &   ! clearsky upward longwave radiative flux (W/m2)
       lwDownClearSky , & ! clearsky downward longwave radiative flux (W/m2)
       swUp , &           ! upward shortwave radiative flux (W/m2)
       swDown , &         ! downward shortwave radiative flux (W/m2)
       swUpClearSky , &   ! clearsky upward shortwave radiative flux (W/m2)
       swDownClearSky     ! clearsky downward shortwave radiative flux (W/m2)

! Shortwave and longwave heating rate outputs in K/day
      real (kind=kind_rb), dimension(:,:), allocatable :: &
          swHeatingRate, &
          swHeatingRateClearSky, &
          lwHeatingRate, &
          lwHeatingRateClearSky
  
      real(kind=kind_rb), dimension(:,:), allocatable :: &
          LWP, &        ! liquid water path (g/m2)
          IWP, &        ! ice water path (g/m2)
          liquidRe, &   ! effective radius liquid (microns)
          iceRe         ! effective radius ice (microns)

end module rad
