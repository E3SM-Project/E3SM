module SLakeCon

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: SLakeCon
!
! !DESCRIPTION:
! Module containing constants and parameters for the SLake code (CLM4-LISSS, documented in Subin et al. 2011, JAMES)
!
! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use clm_varpar , only : numrad
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !REVISION HISTORY:
! Created by Zack Subin, 2011
!
!EOP
!-----------------------------------------------------------------------

  integer, private :: i

  !------------------------------------------------------------------
  ! Lake Model Constants
  !------------------------------------------------------------------

  ! Non-tuneable constants for the lake model
  real(r8), parameter :: tdmax = 277._r8   ! temperature of maximum water density (K)
                                           ! This is from Hostetler and Bartlein (1990); more updated sources suggest
                                           ! 277.13 K.

  ! Tuneable constants for the lake model
  real(r8), parameter :: emg_lake = 0.97_r8  ! lake emissivity. This
                                           ! is used for both frozen and unfrozen lakes. This is pulled in from CLM4 and
                                           ! the reference is unclear.

  real(r8) :: alblak(numrad) = &           ! albedo frozen lakes by waveband (1=vis, 2=nir)
                        (/0.60_r8, 0.40_r8/) ! Also unclear what the reference is for this.
  real(r8) :: alblakwi(numrad)            ! albedo of melting lakes due to puddling, open water, or white ice
                                          ! From D. Mironov (2010) Boreal Env. Research
                                          ! To revert albedo of melting lakes to the cold snow-free value, set
                                          ! lake_melt_icealb namelist to 0.60, 0.40 like alblak above.

  real(r8), parameter :: calb = 95.6_r8   ! Coefficient for calculating ice "fraction" for lake surface albedo
                                          ! From D. Mironov (2010) Boreal Env. Research
  real(r8) :: betavis = 0.0_r8            ! The fraction of the visible (e.g. vis not nir from atm) sunlight
                                          ! absorbed in ~1 m of water (the surface layer za_lake).
                                          ! This is roughly the fraction over 700 nm but may depend on the details
                                          ! of atmospheric radiative transfer.
                                          ! As long as NIR = 700 nm and up, this can be zero.
  real(r8), parameter :: z0frzlake = 0.001_r8  ! Momentum Roughness length over frozen lakes without snow  (m)
                                          ! Typical value found in the literature, and consistent with Mironov expressions.
                                          ! See e.g. Morris EM 1989, Andreas EL 1987, Guest & Davidson 1991 (as cited in
                                          ! Vavrus 1996)
  real(r8), parameter :: za_lake = 0.6_r8           ! Base of surface light absorption layer for lakes (m)

  ! For calculating prognostic roughness length
  real(r8), parameter :: cur0    = 0.01_r8  ! min. Charnock parameter
  real(r8), parameter :: cus     = 0.1_r8   ! empirical constant for roughness under smooth flow
  real(r8), parameter :: curm    = 0.1_r8   ! maximum Charnock parameter

  !! The following will be set in initSLake based on namelists.
  real(r8)            :: fcrit              ! critical dimensionless fetch for Charnock parameter.
  real(r8)            :: minz0lake          ! (m) Minimum allowed roughness length for unfrozen lakes.

  ! For calculating enhanced diffusivity
  real(r8), parameter :: n2min = 7.5e-5_r8 ! (s^-2) (yields diffusivity about 6 times km) ! Fang & Stefan 1996

  real(r8)            :: lsadz = 0.03_r8 ! m
  ! Note, this will be adjusted in initSLake if the timestep is not 1800 s.
  ! Lake top numerics can oscillate with 0.01m top layer and 1800 s timestep.
  ! The problem is that the surface flux is fixed during the calculation of the top
  ! layer temperature in the diffusion and not corrected for the tendency of the top layer.
  ! This thickness will be added to all minimum and maximum snow layer thicknesses compared to that used over non-lakes.
  ! Analysis of the CFL condition suggests that the minimum snow layer thickness for 1800 s needs
  ! to be at least ~1.2 cm for the bulk snow values of conductivity and heat capacity
  ! and as much as 2.3 cm for pure ice.
  ! Alternatively, a check could be done in SLakeTemperature in case
  ! t_grnd(c) - t_soisno(c,snl(c)+1) changed sign after the Crank-Nicholson step.
  ! Such an approach, while perhaps allowing additional snow layer resolution, has not been tested.
  ! The approach used over non-lakes is to have a first-order surface flux correction.
  ! We choose not to do that here because t_grnd can vary independently of the top model
  ! layer temperature, while it is fixed to the top layer temperature if tbot > tfrz and
  ! the lake is frozen, or if there is an unstable density gradient in the top unfrozen
  ! lake layer.

  !! The following will be set in initSLake based on namelists.
  real(r8)            :: pudz          ! (m) Optional minimum total ice thickness required to allow lake puddling.
                                       ! Currently used for sensitivity tests only.
  real(r8)            :: depthcrit     ! (m) Depth beneath which to increase mixing. See discussion in Subin et al. 2011
  real(r8)            :: mixfact       ! Mixing increase factor.

  !!!!!!!!!!!
  ! Namelists (some of these have not been extensively tested and are hardwired to default values currently).
  !!!!!!!!!!!

  ! used in SLakeFluxes
  logical,  public :: lake_use_old_fcrit_minz0 = .false. ! true => use old fcrit & minz0 as per Subin et al 2011 form
                                                         ! See initSLakeMod for details. Difference is very small for
                                                         ! small lakes and negligible for large lakes.
                                                         ! Currently hardwired off.

  ! used in SLakeTemperature
  ! Increase mixing by a large factor for deep lakes
  ! Crude but enhanced performance at all 4 deep lakes tested.
  ! See Subin et al 2011 (JAMES) for details
  real(r8), public :: deepmixing_depthcrit = 25._r8     ! (m) minimum lake depth to invoke deepmixing
  real(r8), public :: deepmixing_mixfact   = 10._r8     ! factor to increase mixing by
logical,  public :: lake_no_ed = .false.              ! true => Suppress enhanced diffusion. Small differences.
                                                        ! Currently hardwired .false.
                                                        ! See Subin et al 2011 for details.
                                                        ! Enhanced diffusion is intended for under ice and at large depths.
                                                        ! It is a much smaller change on its own than the "deepmixing"
                                                        ! above, but it increases the effect of deepmixing under ice and for
                                                        ! large depths.


  ! puddling (not extensively tested and currently hardwired off)
  ! used in SLakeTemperature and SurfaceAlbedo
  logical,  public :: lakepuddling = .false.            ! true => suppress convection when greater than minimum amount
                                                        ! of ice is present. This also effectively sets lake_no_melt_icealb.
  real(r8), public :: lake_puddle_thick = 0.2_r8        ! (m) minimum amount of total ice nominal thickness before
                                                        ! convection is suppressed

  ! alblakwi used in SurfaceAlbedo. Will be set by lake_melt_icealb in initSLake.
  real(r8), public :: lake_melt_icealb(numrad) = &      ! Namelist for inputting alblakwi
                      (/ 0.10_r8, 0.10_r8/)

end module SLakeCon
