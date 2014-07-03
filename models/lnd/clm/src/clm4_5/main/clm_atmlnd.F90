module clm_atmlnd

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle atm2lnd, lnd2atm mapping
  !
  ! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar    , only : numrad, ndst, nlevgrnd !ndst = number of dust bins.
  use clm_varcon    , only : rair, grav, cpair, hfus, tfrz
  use clm_varctl    , only : iulog, use_c13, use_cn, use_lch4, use_voc
  use seq_drydep_mod, only : n_drydep, drydep_method, DD_XLND
  use shr_megan_mod , only : shr_megan_mechcomps_n
  use decompMod     , only : bounds_type
  use shr_log_mod   , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

  !----------------------------------------------------
  ! atmosphere -> land variables structure
  !
  ! This structure contains fields that are NOT downscaled - so it's fine for code to
  ! directly use the gridcell-level versions.
  !----------------------------------------------------
  type, public :: atm2lnd_type
     real(r8), pointer :: forc_u(:)        => null() !atm wind speed, east direction (m/s)
     real(r8), pointer :: forc_v(:)        => null() !atm wind speed, north direction (m/s)
     real(r8), pointer :: forc_wind(:)     => null() !atmospheric wind speed
     real(r8), pointer :: forc_hgt(:)      => null() !atmospheric reference height (m)
     real(r8), pointer :: forc_hgt_u(:)    => null() !obs height of wind [m] (new)
     real(r8), pointer :: forc_hgt_t(:)    => null() !obs height of temperature [m] (new)
     real(r8), pointer :: forc_hgt_q(:)    => null() !obs height of humidity [m] (new)
     real(r8), pointer :: forc_vp(:)       => null() !atmospheric vapor pressure (Pa)
     real(r8), pointer :: forc_rh(:)       => null() !atmospheric relative humidity (%)
     real(r8), pointer :: forc_psrf(:)     => null() !surface pressure (Pa)
     real(r8), pointer :: forc_pco2(:)     => null() !CO2 partial pressure (Pa)
     real(r8), pointer :: forc_solad(:,:)  => null() !direct beam radiation (numrad) (vis=forc_sols , nir=forc_soll )
     real(r8), pointer :: forc_solai(:,:)  => null() !diffuse radiation (numrad) (vis=forc_solsd, nir=forc_solld)
     real(r8), pointer :: forc_solar(:)    => null() !incident solar radiation
     real(r8), pointer :: forc_ndep(:)     => null() !nitrogen deposition rate (gN/m2/s)
     real(r8), pointer :: forc_pc13o2(:)   => null() !C13O2 partial pressure (Pa)
     real(r8), pointer :: forc_po2(:)      => null() !O2 partial pressure (Pa)
     real(r8), pointer :: forc_flood(:)    => null() !rof flood (mm/s)
     real(r8), pointer :: volr(:)          => null() !rof volr (m3)
     real(r8), pointer :: forc_aer(:,:)    => null() !aerosol deposition array
     real(r8), pointer :: forc_pch4(:)     => null() !CH4 partial pressure (Pa)
     ! anomaly forcing
     real(r8), pointer ::af_precip(:)      => null() ! anomaly forcing 
     real(r8), pointer ::af_uwind(:)       => null() ! anomaly forcing 
     real(r8), pointer ::af_vwind(:)       => null() ! anomaly forcing 
     real(r8), pointer ::af_tbot(:)        => null() ! anomaly forcing 
     real(r8), pointer ::af_pbot(:)        => null() ! anomaly forcing 
     real(r8), pointer ::af_shum(:)        => null() ! anomaly forcing 
     real(r8), pointer ::af_swdn(:)        => null() ! anomaly forcing 
     real(r8), pointer ::af_lwdn(:)        => null() ! anomaly forcing 
     real(r8), pointer :: bc_precip(:)     => null() ! anomaly forcing - add bias correction
  end type atm2lnd_type

  !----------------------------------------------------
  ! atmosphere -> land variables structure - fields that are downscaled
  ! ----------------------------------------------------
  type, public :: atm2lnd_downscaled_fields_type
     real(r8), pointer :: forc_t(:)        => null() !atmospheric temperature (Kelvin)
     real(r8), pointer :: forc_th(:)       => null() !atm potential temperature (Kelvin)
     real(r8), pointer :: forc_q(:)        => null() !atmospheric specific humidity (kg/kg)
     real(r8), pointer :: forc_pbot(:)     => null() !atmospheric pressure (Pa)
     real(r8), pointer :: forc_rho(:)      => null() !density (kg/m**3)
     real(r8), pointer :: forc_rain(:)     => null() !rain rate [mm/s]
     real(r8), pointer :: forc_snow(:)     => null() !snow rate [mm/s]
     real(r8), pointer :: forc_lwrad(:)    => null() !downwrd IR longwave radiation (W/m**2)
  end type atm2lnd_downscaled_fields_type


  !----------------------------------------------------
  ! land -> atmosphere variables structure
  !----------------------------------------------------
  type, public :: lnd2atm_type
     real(r8), pointer :: t_rad(:)         => null() !radiative temperature (Kelvin)
     real(r8), pointer :: t_ref2m(:)       => null() !2m surface air temperature (Kelvin)
     real(r8), pointer :: q_ref2m(:)       => null() !2m surface specific humidity (kg/kg)
     real(r8), pointer :: u_ref10m(:)      => null() !10m surface wind speed (m/sec)
     real(r8), pointer :: h2osno(:)        => null() !snow water (mm H2O)
     real(r8), pointer :: albd(:,:)        => null() !(numrad) surface albedo (direct)
     real(r8), pointer :: albi(:,:)        => null() !(numrad) surface albedo (diffuse)
     real(r8), pointer :: taux(:)          => null() !wind stress: e-w (kg/m/s**2)
     real(r8), pointer :: tauy(:)          => null() !wind stress: n-s (kg/m/s**2)
     real(r8), pointer :: eflx_lh_tot(:)   => null() !total latent HF (W/m**2)  [+ to atm]
     real(r8), pointer :: eflx_sh_tot(:)   => null() !total sensible HF (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lwrad_out(:)=> null() !IR (longwave) radiation (W/m**2)
     real(r8), pointer :: qflx_evap_tot(:) => null() !qflx_evap_soi + qflx_evap_can + qflx_tran_veg
     real(r8), pointer :: fsa(:)           => null() !solar rad absorbed (total) (W/m**2)
     real(r8), pointer :: nee(:)           => null() !net CO2 flux (kg CO2/m**2/s) [+ to atm]
     real(r8), pointer :: ram1(:)          => null() !aerodynamical resistance (s/m)
     real(r8), pointer :: fv(:)            => null() !friction velocity (m/s) (for dust model)
     real(r8), pointer :: h2osoi_vol(:,:)  => null() !volumetric soil water (0~watsat, m3/m3, nlevgrnd) (for dust model)
     real(r8), pointer :: rofliq(:)        => null() !rof liq forcing
     real(r8), pointer :: rofice(:)        => null() !rof ice forcing
     real(r8), pointer :: flxdst(:,:)      => null() !dust flux (size bins)
     real(r8), pointer :: ddvel(:,:)       => null() !dry deposition velocities
     real(r8), pointer :: flxvoc(:,:)      => null() !VOC flux (size bins)
     real(r8), pointer :: flux_ch4(:)      => null() !net CH4 flux (kg C/m**2/s) [+ to atm]
  end type lnd2atm_type

  ! l2a fields on clm grid
  type(lnd2atm_type),public,target :: clm_l2a

  ! a2l fields on clm grid
  !
  ! This variable contains fields that are not downscaled (so only exist at the grid cell
  ! level); you may use this directly in your code (indexing fields by g).
  type(atm2lnd_type),public,target :: clm_a2l

  
  ! More a2l fields
  ! 
  ! This variable contains fields that are downscaled to the column level. This is the
  ! version of the downscaled fields that should be used by most science code (indexing
  ! fields by c).
  type(atm2lnd_downscaled_fields_type),public,target :: a2l_downscaled_col


  ! DO NOT USE THIS VARIABLE IN YOUR CODE!
  ! 
  ! This variable contains the non-downscaled version of forcing fields. This is currently
  ! used in a handful of places in the code (and needs to be used in lnd_import_export and
  ! the downscaling routines), but in general should NOT be used in new code. Instead use
  ! a2l_downscaled_col, which gives the downscaled versions of these fields.
  !
  ! DO NOT USE THIS VARIABLE IN YOUR CODE!
  type(atm2lnd_downscaled_fields_type),public,target :: a2l_not_downscaled_gcell
  

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: init_atm2lnd              ! Initialize all atmospheric variables required by the land
  public :: init_lnd2atm              ! Initialize all land variables required by the atmosphere
  public :: clm_map2gcell_minimal
  public :: clm_map2gcell
  public :: downscale_forcings        ! Downscale atm forcing fields from gridcell to column

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: init_atm2lnd_type                    ! Initialize atm2lnd derived type variable
  private :: init_atm2lnd_downscaled_fields_type  ! Initialize atm2lnd_downscaled_fields derived type variable
  private :: init_lnd2atm_type                    ! Initialize lnd2atm derived type variable
  private :: downscale_longwave                   ! Downscale longwave radiation from gridcell to column
  private :: build_normalization                  ! Compute normalization factors so that downscaled fields are conservative
  private :: check_downscale_consistency          ! Check consistency of downscaling
  !----------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine init_atm2lnd(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize all atmospheric variables required by the land
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !------------------------------------------------------------------------

    call init_atm2lnd_type(bounds%begg, bounds%endg, clm_a2l)

    call init_atm2lnd_downscaled_fields_type(bounds%begg, bounds%endg, a2l_not_downscaled_gcell)
    call init_atm2lnd_downscaled_fields_type(bounds%begc, bounds%endc, a2l_downscaled_col)

  end subroutine init_atm2lnd

  !------------------------------------------------------------------------
  subroutine init_atm2lnd_type(begg, endg, a2l)
    !
    ! !DESCRIPTION:
    ! Initialize atm2lnd derived type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: begg              ! beginning grid cell index
    integer, intent(in) :: endg              ! ending grid cell index
    type (atm2lnd_type), intent(inout):: a2l
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival   ! initial value
    !------------------------------------------------------------------------

    ! ival = nan      ! causes core dump in map_maparray, tcx fix
    ival = 0.0_r8

    allocate(a2l%forc_u(begg:endg))
    a2l%forc_u(begg:endg)=ival
    allocate(a2l%forc_v(begg:endg))
    a2l%forc_v(begg:endg)=ival
    allocate(a2l%forc_wind(begg:endg))
    a2l%forc_wind(begg:endg)=ival
    allocate(a2l%forc_rh(begg:endg))
    a2l%forc_rh(begg:endg)=ival
    allocate(a2l%forc_hgt(begg:endg))
    a2l%forc_hgt(begg:endg)=ival
    allocate(a2l%forc_hgt_u(begg:endg))
    a2l%forc_hgt_u(begg:endg)=ival
    allocate(a2l%forc_hgt_t(begg:endg))
    a2l%forc_hgt_t(begg:endg)=ival
    allocate(a2l%forc_hgt_q(begg:endg))
    a2l%forc_hgt_q(begg:endg)=ival
    allocate(a2l%forc_vp(begg:endg))
    a2l%forc_vp(begg:endg)=ival
    allocate(a2l%forc_psrf(begg:endg))
    a2l%forc_psrf(begg:endg)=ival
    allocate(a2l%forc_pco2(begg:endg))
    a2l%forc_pco2(begg:endg)=ival
    allocate(a2l%forc_solad(begg:endg,numrad))
    a2l%forc_solad(begg:endg,numrad)=ival
    allocate(a2l%forc_solai(begg:endg,numrad))
    a2l%forc_solai(begg:endg,numrad)=ival
    allocate(a2l%forc_solar(begg:endg))
    a2l%forc_solar(begg:endg)=ival
    allocate(a2l%forc_ndep(begg:endg))
    a2l%forc_ndep(begg:endg)=ival
    allocate(a2l%forc_pc13o2(begg:endg))
    a2l%forc_pc13o2(begg:endg)=ival
    allocate(a2l%forc_po2(begg:endg))
    a2l%forc_po2(begg:endg)=ival
    allocate(a2l%forc_flood(begg:endg))
    a2l%forc_flood(begg:endg)=ival
    allocate(a2l%volr(begg:endg))
    a2l%volr(begg:endg)=ival
    allocate(a2l%forc_aer(begg:endg,14))
    a2l%forc_aer(begg:endg,14)=ival
    allocate(a2l%forc_pch4(begg:endg))
    a2l%forc_pch4(begg:endg)=ival
    ! anomaly forcing
    allocate(a2l%bc_precip(begg:endg))
    a2l%bc_precip(begg:endg) = ival
    allocate(a2l%af_precip(begg:endg))
    a2l%af_precip(begg:endg) = ival
    allocate(a2l%af_uwind(begg:endg))
    a2l%af_uwind(begg:endg) = ival
    allocate(a2l%af_vwind(begg:endg))
    a2l%af_vwind(begg:endg) = ival
    allocate(a2l%af_tbot(begg:endg))
    a2l%af_tbot(begg:endg) = ival
    allocate(a2l%af_pbot(begg:endg))
    a2l%af_pbot(begg:endg) = ival
    allocate(a2l%af_shum(begg:endg))
    a2l%af_shum(begg:endg) = ival
    allocate(a2l%af_swdn(begg:endg))
    a2l%af_swdn(begg:endg) = ival
    allocate(a2l%af_lwdn(begg:endg))
    a2l%af_lwdn(begg:endg) = ival

  end subroutine init_atm2lnd_type

  !------------------------------------------------------------------------
  subroutine init_atm2lnd_downscaled_fields_type(beg, end, a2l)
    !
    ! !DESCRIPTION:
    ! Initialize atm2lnd_downscaled_fields derived type
    !
    ! This can be used to initialize the gridcell-level fields (in which case beg=begg,
    ! end=endg) or the column-level fields (in which case beg=begc, end=endc)
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg              ! beginning index
    integer, intent(in) :: end              ! ending index
    type (atm2lnd_downscaled_fields_type), intent(inout):: a2l
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival   ! initial value
    !------------------------------------------------------------------------

    ! ival = nan      ! causes core dump in map_maparray, tcx fix
    ival = 0.0_r8

    allocate(a2l%forc_t(beg:end))
    a2l%forc_t(beg:end)=ival
    allocate(a2l%forc_q(beg:end))
    a2l%forc_q(beg:end)=ival
    allocate(a2l%forc_pbot(beg:end))
    a2l%forc_pbot(beg:end)=ival
    allocate(a2l%forc_th(beg:end))
    a2l%forc_th(beg:end)=ival
    allocate(a2l%forc_rho(beg:end))
    a2l%forc_rho(beg:end)=ival
    allocate(a2l%forc_lwrad(beg:end))
    a2l%forc_lwrad(beg:end)=ival
    allocate(a2l%forc_rain(beg:end))
    a2l%forc_rain(beg:end)=ival
    allocate(a2l%forc_snow(beg:end))
    a2l%forc_snow(beg:end)=ival
    
  end subroutine init_atm2lnd_downscaled_fields_type


  !------------------------------------------------------------------------
  subroutine init_lnd2atm(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize all land variables required by the atmosphere
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !------------------------------------------------------------------------

    call init_lnd2atm_type(bounds%begg, bounds%endg, clm_l2a)

  end subroutine init_lnd2atm


  !------------------------------------------------------------------------
  subroutine init_lnd2atm_type(begg, endg, l2a)
    !
    ! !DESCRIPTION:
    ! Initialize lnd2atm derived type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: begg              ! beginning grid cell index
    integer, intent(in) :: endg              ! ending grid cell index
    type (lnd2atm_type), intent(inout):: l2a
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival   ! initial value
    !------------------------------------------------------------------------

    ! ival = nan   ! causes core dump in map_maparray, tcx fix
    ival = 0.0_r8

    allocate(l2a%t_rad(begg:endg))
    l2a%t_rad(begg:endg)=ival
    allocate(l2a%t_ref2m(begg:endg))
    l2a%t_ref2m(begg:endg)=ival
    allocate(l2a%q_ref2m(begg:endg))
    l2a%q_ref2m(begg:endg)=ival
    allocate(l2a%u_ref10m(begg:endg))
    l2a%u_ref10m(begg:endg)=ival
    allocate(l2a%h2osno(begg:endg))
    l2a%h2osno(begg:endg)=ival
    allocate(l2a%albd(begg:endg,1:numrad))
    l2a%albd(begg:endg,1:numrad)=ival
    allocate(l2a%albi(begg:endg,1:numrad))
    l2a%albi(begg:endg,1:numrad)=ival
    allocate(l2a%taux(begg:endg))
    l2a%taux(begg:endg)=ival
    allocate(l2a%tauy(begg:endg))
    l2a%tauy(begg:endg)=ival
    allocate(l2a%eflx_lwrad_out(begg:endg))
    l2a%eflx_lwrad_out(begg:endg)=ival
    allocate(l2a%eflx_sh_tot(begg:endg))
    l2a%eflx_sh_tot(begg:endg)=ival
    allocate(l2a%eflx_lh_tot(begg:endg))
    l2a%eflx_lh_tot(begg:endg)=ival
    allocate(l2a%qflx_evap_tot(begg:endg))
    l2a%qflx_evap_tot(begg:endg)=ival
    allocate(l2a%fsa(begg:endg))
    l2a%fsa(begg:endg)=ival
    allocate(l2a%nee(begg:endg))
    l2a%nee(begg:endg)=ival
    allocate(l2a%ram1(begg:endg))
    l2a%ram1(begg:endg)=ival
    allocate(l2a%fv(begg:endg))
    l2a%fv(begg:endg)=ival
    allocate(l2a%h2osoi_vol(begg:endg,1:nlevgrnd))
    l2a%h2osoi_vol(begg:endg,1:nlevgrnd)=ival
    allocate(l2a%rofliq(begg:endg))
    l2a%rofliq(begg:endg)=ival
    allocate(l2a%rofice(begg:endg))
    l2a%rofice(begg:endg)=ival
    allocate(l2a%flxdst(begg:endg,1:ndst))
    l2a%flxdst(begg:endg,1:ndst)=ival
    allocate(l2a%flux_ch4(begg:endg))
    l2a%flux_ch4(begg:endg)=ival
    if (shr_megan_mechcomps_n>0) then
       allocate(l2a%flxvoc(begg:endg,1:shr_megan_mechcomps_n))
       l2a%flxvoc(begg:endg,1:shr_megan_mechcomps_n)=ival
    endif
    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
       allocate(l2a%ddvel(begg:endg,1:n_drydep))
       l2a%ddvel(begg:endg,1:n_drydep)=ival
    end if

  end subroutine init_lnd2atm_type

  !------------------------------------------------------------------------

  subroutine clm_map2gcell_minimal(bounds)
    !
    ! !DESCRIPTION:
    ! Compute l2a component of gridcell derived type. This routine computes
    ! the bare minimum of components necessary to get the first step of a
    ! run started.
    !
    ! !USES:
    use clmtype
    use subgridAveMod
    use clm_varcon  , only : sb
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g             ! index
    real(r8), parameter :: amC   = 12.0_r8 ! Atomic mass number for Carbon
    real(r8), parameter :: amO   = 16.0_r8 ! Atomic mass number for Oxygen
    real(r8), parameter :: amCO2 = amC + 2.0_r8*amO ! Atomic mass number for CO2
    ! The following converts g of C to kg of CO2
    real(r8), parameter :: convertgC2kgCO2 = 1.0e-3_r8 * (amCO2/amC)
    !------------------------------------------------------------------------

    call c2g(bounds, &
         cws%h2osno(bounds%begc:bounds%endc), &
         clm_l2a%h2osno(bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    do g = bounds%begg,bounds%endg
       clm_l2a%h2osno(g) = clm_l2a%h2osno(g)/1000._r8
    end do

    call c2g(bounds, nlevgrnd, &
         cws%h2osoi_vol(bounds%begc:bounds%endc, :), &
         clm_l2a%h2osoi_vol(bounds%begg:bounds%endg, :), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, numrad, &
         pps%albd(bounds%begp:bounds%endp, :), &
         clm_l2a%albd(bounds%begg:bounds%endg, :), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, numrad, &
         pps%albi(bounds%begp:bounds%endp, :), &
         clm_l2a%albi(bounds%begg:bounds%endg, :), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, &
         pef%eflx_lwrad_out(bounds%begp:bounds%endp), &
         clm_l2a%eflx_lwrad_out(bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    do g = bounds%begg,bounds%endg
       clm_l2a%t_rad(g) = sqrt(sqrt(clm_l2a%eflx_lwrad_out(g)/sb))
    end do

  end subroutine clm_map2gcell_minimal

  !------------------------------------------------------------------------

  subroutine clm_map2gcell(bounds)
    !
    ! !DESCRIPTION:
    ! Compute l2a component of gridcell derived type
    !
    ! !USES:
    use clmtype
    use subgridAveMod
    use ch4varcon   , only : ch4offline
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g             ! index
    real(r8), parameter :: amC   = 12.0_r8 ! Atomic mass number for Carbon
    real(r8), parameter :: amO   = 16.0_r8 ! Atomic mass number for Oxygen
    real(r8), parameter :: amCO2 = amC + 2.0_r8*amO ! Atomic mass number for CO2
    ! The following converts g of C to kg of CO2
    real(r8), parameter :: convertgC2kgCO2 = 1.0e-3_r8 * (amCO2/amC)
    !------------------------------------------------------------------------

    ! First, compute the "minimal" set of fields.
    call clm_map2gcell_minimal(bounds)

    call p2g(bounds, &
         pes%t_ref2m(bounds%begp:bounds%endp), &
         clm_l2a%t_ref2m(bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         pes%q_ref2m(bounds%begp:bounds%endp), &
         clm_l2a%q_ref2m(bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         pps%u10_clm(bounds%begp:bounds%endp), &
         clm_l2a%u_ref10m(bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         pmf%taux(bounds%begp:bounds%endp), &
         clm_l2a%taux(bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         pmf%tauy(bounds%begp:bounds%endp), &
         clm_l2a%tauy(bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         pef%eflx_lh_tot(bounds%begp:bounds%endp), &
         clm_l2a%eflx_lh_tot(bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    do g = bounds%begg,bounds%endg
       clm_l2a%eflx_sh_tot(g) = gef%eflx_sh_totg(g)
    end do

    call p2g(bounds, &
         pwf%qflx_evap_tot(bounds%begp:bounds%endp), &
         clm_l2a%qflx_evap_tot(bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, &
         pef%fsa(bounds%begp:bounds%endp), &
         clm_l2a%fsa(bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    if (use_cn) then
       call c2g(bounds, &
            ccf%nee(bounds%begc:bounds%endc), &
            clm_l2a%nee(bounds%begg:bounds%endg), &
            c2l_scale_type= 'unity', l2g_scale_type='unity')
    else
       call p2g(bounds, &
            pcf%fco2(bounds%begp:bounds%endp), &
            clm_l2a%nee(bounds%begg:bounds%endg), &
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
       ! Note that fco2 in is umolC/m2/sec so units need to be changed to gC/m2/sec
       do g = bounds%begg,bounds%endg
          clm_l2a%nee(g) = clm_l2a%nee(g)*12.011e-6_r8
       end do
    end if

    if (use_lch4) then
       if (.not. ch4offline) then
         ! Adjust flux of CO2 by the net conversion of mineralizing C to CH4
          do g = bounds%begg,bounds%endg
             clm_l2a%nee(g) = clm_l2a%nee(g) + gch4%nem(g) ! nem is in g C/m2/sec
             ! nem is calculated in ch4Mod
             ! flux_ch4 is averaged there also.
          end do
       end if
    end if

    call p2g(bounds, &
         pps%fv(bounds%begp:bounds%endp), &
         clm_l2a%fv(bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
       pps%ram1(bounds%begp:bounds%endp), &
       clm_l2a%ram1(bounds%begg:bounds%endg), &
       p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    do g = bounds%begg,bounds%endg
       clm_l2a%rofliq(g) = gwf%qflx_runoffg(g)
       clm_l2a%rofice(g) = gwf%qflx_snwcp_iceg(g)
    end do

    call p2g(bounds, ndst, &
       pdf%flx_mss_vrt_dst(bounds%begp:bounds%endp, :), &
       clm_l2a%flxdst(bounds%begg:bounds%endg, :), &
       p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    if (use_voc) then
       if (shr_megan_mechcomps_n>0) then
          call p2g(bounds, shr_megan_mechcomps_n, &
               pvf%vocflx(bounds%begp:bounds%endp, :), &
               clm_l2a%flxvoc(bounds%begg:bounds%endg, :), &
               p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
       endif
    end if

    if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
       call p2g(bounds, n_drydep, &
           pdd%drydepvel(bounds%begp:bounds%endp, :), &
           clm_l2a%ddvel(bounds%begg:bounds%endg, :), &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
    endif

    ! Convert from gC/m2/s to kgCO2/m2/s
    do g = bounds%begg,bounds%endg
       clm_l2a%nee(g) = clm_l2a%nee(g)*convertgC2kgCO2
    end do

  end subroutine clm_map2gcell

  !-----------------------------------------------------------------------
  subroutine downscale_forcings(bounds,num_do_smb_c,filter_do_smb_c)
    !
    ! !DESCRIPTION:
    ! Downscale atmospheric forcing fields from gridcell to column
    !
    ! Downscaling is done over columns defined by filter_do_smb_c. But we also do direct copies
    ! of gridcell-level forcings into column-level forcings over all other active columns.
    !
    ! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype
    use clm_varcon   , only : rair, cpair, grav, istice_mec, lapse_glcmec, &
         glcmec_rain_snow_threshold
    use domainMod    , only : ldomain
    use QsatMod      , only : Qsat
    use decompMod    , only : bounds_type
    use clm_varctl   , only : iulog, glcmec_downscale_rain_snow_convert
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: num_do_smb_c       ! number of columns in filter_do_smb_c
    integer, intent(in) :: filter_do_smb_c(:) ! filter_do_smb_c giving columns over which downscaling should be done   
    
    ! !LOCAL VARIABLES:
    integer :: g, l, c, fc         ! indices
    integer :: clo, cc

    ! temporaries for topo downscaling
    real(r8) :: hsurf_g,hsurf_c,Hbot
    real(r8) :: zbot_g, tbot_g, pbot_g, thbot_g, qbot_g, qs_g, es_g
    real(r8) :: zbot_c, tbot_c, pbot_c, thbot_c, qbot_c, qs_c, es_c
    real(r8) :: egcm_c, rhos_c
    real(r8) :: dum1,   dum2

    character(len=*), parameter :: subname = 'downscale_forcings'
    !-----------------------------------------------------------------------

    associate(& 
    glc_topo     => cps%glc_topo                     , & ! Input:  [real(r8) (:)]  sfc elevation for glacier_mec column (m)
   
    ! Gridcell-level fields:
    forc_t_g     => a2l_not_downscaled_gcell%forc_t    , & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)        
    forc_th_g    => a2l_not_downscaled_gcell%forc_th   , & ! Input:  [real(r8) (:)]  atmospheric potential temperature (Kelvin)
    forc_q_g     => a2l_not_downscaled_gcell%forc_q    , & ! Input:  [real(r8) (:)]  atmospheric specific humidity (kg/kg)   
    forc_pbot_g  => a2l_not_downscaled_gcell%forc_pbot , & ! Input:  [real(r8) (:)]  atmospheric pressure (Pa)               
    forc_rho_g   => a2l_not_downscaled_gcell%forc_rho  , & ! Input:  [real(r8) (:)]  atmospheric density (kg/m**3)           
    forc_rain_g  => a2l_not_downscaled_gcell%forc_rain , & ! Input:  [real(r8) (:)]  rain rate [mm/s]
    forc_snow_g  => a2l_not_downscaled_gcell%forc_snow , & ! Input:  [real(r8) (:)]  snow rate [mm/s]

    ! Column-level (downscaled) fields:
    forc_t_c     => a2l_downscaled_col%forc_t        , & ! Output: [real(r8) (:)]  atmospheric temperature (Kelvin)        
    forc_th_c    => a2l_downscaled_col%forc_th       , & ! Output: [real(r8) (:)]  atmospheric potential temperature (Kelvin)
    forc_q_c     => a2l_downscaled_col%forc_q        , & ! Output: [real(r8) (:)]  atmospheric specific humidity (kg/kg)   
    forc_pbot_c  => a2l_downscaled_col%forc_pbot     , & ! Output: [real(r8) (:)]  atmospheric pressure (Pa)               
    forc_rho_c   => a2l_downscaled_col%forc_rho      , & ! Output: [real(r8) (:)]  atmospheric density (kg/m**3)           
    forc_rain_c  => a2l_downscaled_col%forc_rain     , & ! Output: [real(r8) (:)]  rain rate [mm/s]
    forc_snow_c  => a2l_downscaled_col%forc_snow       & ! Output: [real(r8) (:)]  snow rate [mm/s]
    )

    ! Initialize column forcing (needs to be done for ALL active columns)
    do c = bounds%begc,bounds%endc
       if (col%active(c)) then
          g = col%gridcell(c)

          forc_t_c(c)     = forc_t_g(g)
          forc_th_c(c)    = forc_th_g(g)
          forc_q_c(c)     = forc_q_g(g)
          forc_pbot_c(c)  = forc_pbot_g(g)
          forc_rho_c(c)   = forc_rho_g(g)
          forc_rain_c(c)  = forc_rain_g(g)
          forc_snow_c(c)  = forc_snow_g(g)
       end if
    end do

    ! Downscale forc_t, forc_th, forc_q, forc_pbot, and forc_rho to columns.
    ! For glacier_mec columns the downscaling is based on surface elevation.
    ! For other columns the downscaling is a simple copy (above).
    
    do fc = 1, num_do_smb_c
       c = filter_do_smb_c(fc)
       l = col%landunit(c)
       g = col%gridcell(c)

       ! This is a simple downscaling procedure 
       ! Note that forc_hgt, forc_u, and forc_v are not downscaled.
       
       hsurf_g      = ldomain%topo(g)     ! gridcell sfc elevation
       hsurf_c      = glc_topo(c)         ! column sfc elevation
       tbot_g       = forc_t_g(g)         ! atm sfc temp
       thbot_g      = forc_th_g(g)        ! atm sfc pot temp
       qbot_g       = forc_q_g(g)         ! atm sfc spec humid
       pbot_g       = forc_pbot_g(g)      ! atm sfc pressure
       zbot_g       = clm_a2l%forc_hgt(g) ! atm ref height
       zbot_c = zbot_g
       tbot_c = tbot_g-lapse_glcmec*(hsurf_c-hsurf_g)   ! sfc temp for column
       Hbot   = rair*0.5_r8*(tbot_g+tbot_c)/grav        ! scale ht at avg temp
       pbot_c = pbot_g*exp(-(hsurf_c-hsurf_g)/Hbot)     ! column sfc press

       ! Derivation of potential temperature calculation:
       ! 
       ! The textbook definition would be:
       ! thbot_c = tbot_c * (p0/pbot_c)^(rair/cpair)
       ! 
       ! Note that pressure is related to scale height as:
       ! pbot_c = p0 * exp(-zbot_c/H)
       !
       ! Using Hbot in place of H, we get:
       ! pbot_c = p0 * exp(-zbot_c/Hbot)
       !
       ! Plugging this in to the textbook definition, then manipulating, we get:
       ! thbot_c = tbot_c * (p0/(p0*exp(-zbot_c/Hbot)))^(rair/cpair)
       !         = tbot_c * (1/exp(-zbot_c/Hbot))^(rair/cpair)
       !         = tbot_c * (exp(zbot_c/Hbot))^(rair/cpair)
       !         = tbot_c * exp((zbot_c/Hbot) * (rair/cpair))

       thbot_c= tbot_c*exp((zbot_c/Hbot)*(rair/cpair))  ! pot temp calc

       call Qsat(tbot_g,pbot_g,es_g,dum1,qs_g,dum2)
       call Qsat(tbot_c,pbot_c,es_c,dum1,qs_c,dum2)

       qbot_c = qbot_g*(qs_c/qs_g)
       egcm_c = qbot_c*pbot_c/(0.622+0.378*qbot_c)
       rhos_c = (pbot_c-0.378*egcm_c) / (rair*tbot_c)

       forc_t_c(c)    = tbot_c
       forc_th_c(c)   = thbot_c
       forc_q_c(c)    = qbot_c
       forc_pbot_c(c) = pbot_c
       forc_rho_c(c)  = rhos_c

       ! Optionally, convert rain to snow or vice versa based on tbot_c
       ! Note: This conversion does not conserve energy.
       !       It would be better to compute the net latent energy associated
       !        with the conversion and to apply it as a pseudo-flux.
       if (glcmec_downscale_rain_snow_convert) then
          if (tbot_c > glcmec_rain_snow_threshold) then  ! too warm for snow
             forc_rain_c(c) = forc_rain_c(c) + forc_snow_c(c)
             forc_snow_c(c) = 0._r8
          else                                           ! too cold for rain
             forc_snow_c(c) = forc_rain_c(c) + forc_snow_c(c)
             forc_rain_c(c) = 0._r8
          endif
       endif   ! glcmec_downscale_rain_snow_convert

    end do    

    call downscale_longwave(bounds, num_do_smb_c, filter_do_smb_c)

    call check_downscale_consistency(bounds)

    end associate
    
  end subroutine downscale_forcings

  !-----------------------------------------------------------------------
  subroutine downscale_longwave(bounds, num_do_smb_c, filter_do_smb_c)
    !
    ! !DESCRIPTION:
    ! Downscale longwave radiation from gridcell to column
    ! Must be done AFTER temperature downscaling
    !
    ! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype
    use clm_varcon   , only : istice_mec, lapse_glcmec
    use domainMod    , only : ldomain
    use clm_varctl   , only : iulog, glcmec_downscale_longwave
    use decompMod    , only : bounds_type
    use abortutils   , only : endrun
    use shr_log_mod  , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: num_do_smb_c       ! number of columns in filter_do_smb_c
    integer, intent(in) :: filter_do_smb_c(:) ! filter_do_smb_c giving columns over which downscaling should be done (currently glcmec columns)
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l,g,fc     ! indices
    real(r8) :: hsurf_c      ! column-level elevation (m)
    real(r8) :: hsurf_g      ! gridcell-level elevation (m)

    real(r8), dimension(bounds%begg : bounds%endg) :: sum_lwrad_g    ! weighted sum of column-level lwrad
    real(r8), dimension(bounds%begg : bounds%endg) :: sum_wts_g      ! sum of weights that contribute to sum_lwrad_g
    real(r8), dimension(bounds%begg : bounds%endg) :: lwrad_norm_g   ! normalization factors
    real(r8), dimension(bounds%begg : bounds%endg) :: newsum_lwrad_g ! weighted sum of column-level lwrad after normalization

    character(len=*), parameter :: subname = 'downscale_longwave'
    !-----------------------------------------------------------------------

    associate(&
    glc_topo     => cps%glc_topo                     , & ! Input:  [real(r8) (:)]  sfc elevation for glacier_mec column (m)

    ! Gridcell-level fields:
    forc_t_g     => a2l_not_downscaled_gcell%forc_t    , & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)        
    forc_lwrad_g => a2l_not_downscaled_gcell%forc_lwrad, & ! Input:  [real(r8) (:)]  downward longwave (W/m**2)

    ! Column-level (downscaled) fields:
    forc_t_c     => a2l_downscaled_col%forc_t        , & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)        
    forc_lwrad_c => a2l_downscaled_col%forc_lwrad      & ! Output: [real(r8) (:)]  downward longwave (W/m**2)
    )

    ! Initialize column forcing (needs to be done for ALL active columns)
    do c = bounds%begc, bounds%endc
       if (col%active(c)) then
          g = col%gridcell(c)
          forc_lwrad_c(c) = forc_lwrad_g(g)
       end if
    end do


    ! Optionally, downscale the longwave radiation, conserving energy
    if (glcmec_downscale_longwave) then
       
       ! Initialize variables related to normalization
       do g = bounds%begg, bounds%endg
          sum_lwrad_g(g) = 0._r8
          sum_wts_g(g) = 0._r8
          newsum_lwrad_g(g) = 0._r8
       end do

       ! Do the downscaling
       do fc = 1, num_do_smb_c
          c = filter_do_smb_c(fc)
          l = col%landunit(c)
          g = col%gridcell(c)

          hsurf_g = ldomain%topo(g)
          hsurf_c = glc_topo(c)
	  
          ! Here we assume that deltaLW = (dLW/dT)*(dT/dz)*deltaz
          ! We get dLW/dT = 4*eps*sigma*T^3 = 4*LW/T from the Stefan-Boltzmann law,
          ! evaluated at the mean temp.
          ! We assume the same temperature lapse rate as above.

          forc_lwrad_c(c) = forc_lwrad_g(g) - &
               4.0_r8 * forc_lwrad_g(g)/(0.5_r8*(forc_t_c(c)+forc_t_g(g))) * &
               lapse_glcmec * (hsurf_c - hsurf_g)

          ! Keep track of the gridcell-level weighted sum for later normalization.
          !
          ! This gridcell-level weighted sum just includes points for which we do the
          ! downscaling (e.g., glc_mec points). Thus the contributing weights
          ! generally do not add to 1. So to do the normalization properly, we also
          ! need to keep track of the weights that have contributed to this sum.
          sum_lwrad_g(g) = sum_lwrad_g(g) + col%wtgcell(c)*forc_lwrad_c(c)
          sum_wts_g(g) = sum_wts_g(g) + col%wtgcell(c)
       end do


       ! Normalize forc_lwrad_c(c) to conserve energy

       call build_normalization(orig_field=forc_lwrad_g(bounds%begg:bounds%endg), &
                                sum_field=sum_lwrad_g(bounds%begg:bounds%endg), &
                                sum_wts=sum_wts_g(bounds%begg:bounds%endg), &
                                norms=lwrad_norm_g(bounds%begg:bounds%endg))

       do fc = 1, num_do_smb_c
          c = filter_do_smb_c(fc)
          l = col%landunit(c)
          g = col%gridcell(c)

          forc_lwrad_c(c) = forc_lwrad_c(c) * lwrad_norm_g(g)
          newsum_lwrad_g(g) = newsum_lwrad_g(g) + col%wtgcell(c)*forc_lwrad_c(c)
       end do


       ! Make sure that, after normalization, the grid cell mean is conserved
       
       do g = bounds%begg, bounds%endg
          if (sum_wts_g(g) > 0._r8) then
             if (abs((newsum_lwrad_g(g) / sum_wts_g(g)) - forc_lwrad_g(g)) > 1.e-8_r8) then
                write(iulog,*) 'g, newsum_lwrad_g, sum_wts_g, forc_lwrad_g: ', &
                     g, newsum_lwrad_g(g), sum_wts_g(g), forc_lwrad_g(g)
                call endrun(msg=' ERROR: Energy conservation error downscaling longwave'//&
                     errMsg(__FILE__, __LINE__))
             end if
          end if
       end do

    end if    ! glcmec_downscale_longwave

    end associate

  end subroutine downscale_longwave


  !-----------------------------------------------------------------------
  subroutine build_normalization(orig_field, sum_field, sum_wts, norms)
    !
    ! !DESCRIPTION:
    ! Build an array of normalization factors that can be applied to a downscaled forcing
    ! field, in order to force the mean of the new field to be the same as the mean of
    ! the old field (for conservation).
    !
    ! This allows for the possibility that only a subset of columns are downscaled. Only
    ! the columns that are adjusted should be included in the weighted sum, sum_field;
    ! sum_wts gives the sum of contributing weights on the grid cell level. 

    ! For example, if a grid cell has an original forcing value of 1.0, and contains 4
    ! columns with the following weights on the gridcell, and the following values after
    ! normalization:
    !
    !       col #:    1     2     3     4
    !      weight:  0.1   0.2   0.3   0.4
    ! downscaled?:  yes   yes    no    no
    !       value:  0.9   1.1   1.0   1.0
    !
    ! Then we would have:
    ! orig_field(g) = 1.0
    ! sum_field(g) = 0.1*0.9 + 0.2*1.1 = 0.31
    ! sum_wts(g) = 0.1 + 0.2 = 0.3
    ! norms(g) = 1.0 / (0.31 / 0.3) = 0.9677
    !
    ! The field can then be normalized as:
    !              forc_lwrad_c(c) = forc_lwrad_c(c) * lwrad_norm_g(g)
    !   where lwrad_norm_g is the array of norms computed by this routine

    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: orig_field(:)  ! the original field, at the grid cell level
    real(r8), intent(in)  :: sum_field(:)   ! the new weighted sum across columns (dimensioned by grid cell)
    real(r8), intent(in)  :: sum_wts(:)     ! sum of the weights used to create sum_field (dimensioned by grid cell)
    real(r8), intent(out) :: norms(:)       ! computed normalization factors

    !-----------------------------------------------------------------------

    SHR_ASSERT((size(orig_field) == size(norms)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT((size(sum_field) == size(norms)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT((size(sum_wts) == size(norms)), errMsg(__FILE__, __LINE__))

    where (sum_wts == 0._r8)
       ! Avoid divide by zero; if sum_wts is 0, then the normalization doesn't matter,
       ! because the adjusted values won't affect the grid cell mean.
       norms = 1.0_r8

    elsewhere (sum_field == 0._r8)
       ! Avoid divide by zero; this should only happen if the gridcell-level value is 0,
       ! in which case the normalization doesn't matter
       norms = 1.0_r8

    elsewhere
       ! The standard case
       norms = orig_field / (sum_field / sum_wts)

    end where

  end subroutine build_normalization


  !-----------------------------------------------------------------------
  subroutine check_downscale_consistency(bounds)
    !
    ! !DESCRIPTION:
    ! Check consistency of downscaling
    !
    ! Note that this operates over more than just the filter used for the downscaling,
    ! because it checks some things outside that filter.
    !
    ! !USES:
    use clmtype
    use abortutils , only : endrun
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use clm_varctl , only : iulog
    use decompMod  , only : bounds_type

    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g, l, c    ! indices
    character(len=*), parameter :: subname = 'check_downscale_consistency'
    !-----------------------------------------------------------------------

    associate(&
    ! Gridcell-level fields:
    forc_t_g     => a2l_not_downscaled_gcell%forc_t    , & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)        
    forc_th_g    => a2l_not_downscaled_gcell%forc_th   , & ! Input:  [real(r8) (:)]  atmospheric potential temperature (Kelvin)
    forc_q_g     => a2l_not_downscaled_gcell%forc_q    , & ! Input:  [real(r8) (:)]  atmospheric specific humidity (kg/kg)   
    forc_pbot_g  => a2l_not_downscaled_gcell%forc_pbot , & ! Input:  [real(r8) (:)]  atmospheric pressure (Pa)               
    forc_rho_g   => a2l_not_downscaled_gcell%forc_rho  , & ! Input:  [real(r8) (:)]  atmospheric density (kg/m**3)           
    forc_rain_g  => a2l_not_downscaled_gcell%forc_rain , & ! Input:  [real(r8) (:)]  rain rate [mm/s]
    forc_snow_g  => a2l_not_downscaled_gcell%forc_snow , & ! Input:  [real(r8) (:)]  snow rate [mm/s]
    forc_lwrad_g => a2l_not_downscaled_gcell%forc_lwrad, & ! Input:  [real(r8) (:)]  downward longwave (W/m**2)

    ! Column-level (downscaled) fields:
    forc_t_c     => a2l_downscaled_col%forc_t        , & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)        
    forc_th_c    => a2l_downscaled_col%forc_th       , & ! Input:  [real(r8) (:)]  atmospheric potential temperature (Kelvin)
    forc_q_c     => a2l_downscaled_col%forc_q        , & ! Input:  [real(r8) (:)]  atmospheric specific humidity (kg/kg)   
    forc_pbot_c  => a2l_downscaled_col%forc_pbot     , & ! Input:  [real(r8) (:)]  atmospheric pressure (Pa)               
    forc_rho_c   => a2l_downscaled_col%forc_rho      , & ! Input:  [real(r8) (:)]  atmospheric density (kg/m**3)           
    forc_rain_c  => a2l_downscaled_col%forc_rain     , & ! Input:  [real(r8) (:)]  rain rate [mm/s]
    forc_snow_c  => a2l_downscaled_col%forc_snow     , & ! Input:  [real(r8) (:)]  snow rate [mm/s]
    forc_lwrad_c => a2l_downscaled_col%forc_lwrad      & ! Input:  [real(r8) (:)]  downward longwave (W/m**2)
    )
    
    ! Make sure that, for urban points, the column-level forcing fields are identical to
    ! the gridcell-level forcing fields. This is needed because the urban-specific code
    ! sometimes uses the gridcell-level forcing fields (and it would take a large
    ! refactor to change this to use column-level fields).
    
    do c = bounds%begc, bounds%endc
       if (col%active(c)) then
          l = col%landunit(c)
          g = col%gridcell(c)

          if (lun%urbpoi(l)) then
             if (forc_t_c(c)     /= forc_t_g(g)    .or. &
                 forc_th_c(c)    /= forc_th_g(g)   .or. &
                 forc_q_c(c)     /= forc_q_g(g)    .or. &
                 forc_pbot_c(c)  /= forc_pbot_g(g) .or. &
                 forc_rho_c(c)   /= forc_rho_g(g)  .or. &
                 forc_rain_c(c)  /= forc_rain_g(g) .or. &
                 forc_snow_c(c)  /= forc_snow_g(g) .or. &
                 forc_lwrad_c(c) /= forc_lwrad_g(g)) then
                write(iulog,*) subname//' ERROR: column-level forcing differs from gridcell-level forcing for urban point'
                write(iulog,*) 'c, g = ', c, g
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end if  ! inequal
          end if  ! urbpoi
       end if  ! active
    end do

    end associate
  end subroutine check_downscale_consistency


end module clm_atmlnd
