      MODULE rad_driver
      
! JUNG: To use cime/src/share/util/shr_orb_mod (R8 data type), 
!       instead of using GQ3D_CODE/RRTMG/RAD_RRTM_CFMIP/shr_orb_mod (R4 data type),
!       some modifications are made (02/23/2018). Follow "JUNG"
!
!       In order to remove MPI structure,
!       - Remove rrtmg_lw_read_nc.f90 and rrtmg_sw_read_nc.f90 
!         Instead, use rrtmg_lw_k_g_constants.f90 and rrtmg_sw_k_g_constants.f90 
!       - Remove subroutine tracesini since subroutine TRACE_GAS_INPUT (in trace_gases.f90)
!         is used ("call tracesini" was commented out by TCRAMS in rad_full.f90).  
!         In TRACE_GAS_INPUT, data is specified, instead of reading.
!       - Comment out all (masterproc) writing
! -------------------------------------------------------------------------- 
! Interface to RRTM longwave and shortwave radiation code.
!   Robert Pincus, November 2007
!  
! Modified by Peter Blossey, July 2009.
!   - interfaced to RRTMG LW v4.8 and SW v3.8.
!   - reversed indices in tracesini to match bottom-up RRTM indexing.
!   - fixed issue w/doperpetual caused by eccf=0.
!   - added extra layer in calls to RRTM to improve heating rates at
!        model top.  Required changes to inatm in rrtmg_lw_rad.f90
!        and inatm_sw in rrtmg_sw_rad.f90.
!   - fixed bug that would zero out cloudFrac if LWP>0 and IWP==0.
!   - changed definition of interfaceT(:,1) to SST.
!   - added o3, etc. profiles to restart file.  Only call tracesini if nrestart==0.
!
! Modified by Peter Blossey, August 2009.
! Changes for CFMIP intercomparison:
!   - converted rad_driver into a routine with input/output
!        arguments for portability of implementation.
!
! -------------------------------------------------------------------------- 

      use shr_orb_mod, only: shr_orb_params
      use cam_rad_parameterizations, only : &
                                     computeRe_Liquid, computeRe_Ice, albedo
      use parkind, only : kind_rb, &   ! (8 byte reals) 
                          kind_rm      ! (4 byte reals)  ! added by JUNG
      implicit none
      private

!JUNG    include "mpif.h"

!------------------------------------------------------------------------------
! Public procedures
!------------------------------------------------------------------------------

      public :: rad_driver_rrtm, initialize_radiation, isInitialized_RadDriver
  
!------------------------------------------------------------------------------
! Constants
!------------------------------------------------------------------------------

      real, parameter :: &
          Pi = 3.14159265358979312, &
          scon = 1367.                 ! solar constant 

! Molecular weights (taken from CAM shrc_const_mod.F90 and physconst.F90)
      real, parameter :: &
          mwdry  = 28.966,   & ! molecular weight dry air
          mwco2  = 44.,      & ! molecular weight co2
          mwh2o  = 18.016,   & ! molecular weight h2o
          mwn2o  = 44.,      & ! molecular weight n2o
          mwch4  = 16.,      & ! molecular weight ch4
          mwf11  = 136.,     & ! molecular weight cfc11
          mwf12  = 120.,     & ! molecular weight cfc12
          mwf22  = 86.,      & ! molecular weight cfc22
          mwccl4 = 154.,     & ! molecular weight ccl4
          mwo3   = 48.,      & ! molecular weight ozone
          mwo2   = 31.9988     ! molecular weight oxygen

! mixingRatioMass = mol_weight/mol_weight_air * mixingRatioVolume
  
!------------------------------------------------------------------------------
! Global storage
!------------------------------------------------------------------------------

      logical :: &
          isInitialized_RadDriver = .false. 

!bloss(072009): changed from mass mixing ratios to volume mixing ratios
!               because we are now using rrtmg_lw.nc sounding for trace gases.

! Profiles of gas volume mixing ratios 

      logical :: &
          isallocated_tracegases
 
      integer :: &
          nz_tracegases
      
      real(kind_rb), dimension(:), allocatable, save :: &
          o3, &       ! ozone
          co2, &      ! carbon dioxide
          ch4, &      ! methane
          n2o, &      ! nitrous oxide
          o2, &       ! oxygen
          cfc11, &    ! CFC-11
          cfc12, &    ! CFC-12
          cfc22, &    ! CFC-22
          ccl4        ! carbon tetrachloride
  
      real, dimension(:,:), allocatable, save :: &
          p_factor_xy, &  ! perpetual-sun factor
          p_coszrs_xy     ! cosine zenith angle
   
      public :: p_factor_xy, p_coszrs_xy
      
! Earth orbital characteristics
! Calculated in shr_orb_mod which is called by rad_driver_rrtm

      real, save ::  &
          eccf,  &     ! eccentricity factor (1./earth-sun dist^2)
          eccen, &     ! Earth eccentricity factor (unitless) (typically 0 to 0.1)
          obliq, &     ! Earth obliquity angle (deg) (-90 to +90) (typically 22-26)
          mvelp        ! Earth moving vernal equinox at perhelion (deg)(0 to 360.0)

! Orbital information after processed by orbit_params
      real, save ::  &
          obliqr, &  ! Earth obliquity in radians
          lambm0, &  ! Mean longitude of perihelion at the vernal equinox (radians)
          mvelpp     ! Earth moving vernal equinox longitude of perihelion plus pi (radians)


!==============================================================================
      CONTAINS 
!==============================================================================

      subroutine rad_driver_rrtm(nx,nzm,lat,pres,presi,tabs,qv,qcl,qci,tg, &
                                 o3, co2, ch4, n2o, o2, cfc11, cfc12, cfc22, ccl4, &
                                 dolongwave,doshortwave,doperpetual,doseasons, &
                                 dosolarconstant,solar_constant,zenith_angle, &
                                 day,day0,latitude,longitude, &
                                 p_factor_slice, p_coszrs_slice, &
                                 ocean,ggr,cp, masterproc, &
                                 lwUp,lwDown,lwUpClearSky,lwDownClearSky, &
                                 swUp,swDown,swUpClearSky,swDownClearSky, &
                                 swHeatingRate, &
                                 swHeatingRateClearSky, &
                                 lwHeatingRate, &
                                 lwHeatingRateClearSky, &
                                 coszrs, &
                                 LWP, IWP, liquidRe, iceRe)

! Astronomy module, for computing solar zenith angle
      use shr_orb_mod, only: shr_orb_params, shr_orb_decl
    
! Radiation solvers
      use rrtmg_sw_rad, only : rrtmg_sw
      use rrtmg_lw_rad, only : rrtmg_lw
      use parrrtm,      only : nbndlw ! Number of LW bands
      use parrrsw,      only : nbndsw, naerec ! Number of SW bands
      use parkind,      only : kind_rb, & ! (8 byte reals) 
                               kind_rm    ! (4 byte reals)   ! Added by JUNG 
      
      IMPLICIT NONE
    
!------------------------------------------------------------------------------
! Input arguments
!------------------------------------------------------------------------------

      integer, intent(in) :: &
          nx, &    ! number of columns for which radiation will be computed
          nzm, &   ! number of model levels in each column.
          lat      ! index of coordinate in y-direction (from 1 to ny).
      
      logical, intent(in) :: &
          masterproc, &        ! .true. if MPI rank==0 (or single-processor run)
          dolongwave, &        ! compute longwave radiation
          doshortwave, &       ! compute shortwave radiation
          doperpetual, &       ! use perpetual (diurnally-averaged) insolation
          doseasons, &         ! allow diurnally-varying insolation to vary with time of year
          dosolarconstant, &   ! allow insolation and zenith angle to be specified if doperpetual==.true.
          ocean                ! .true. if ocean surface, .false. if land
      
      real, intent(in) :: &
          solar_constant, &   ! mean insolation if doperpetual==dosolarconstant==.true.
          zenith_angle, &     ! solar zenith angle if doperpetual==dosolarconstant==.true.
          day, &              ! day of year during iyear (0.0 = 00Z Jan 1)
          day0, &             ! starting day of year for run
          ggr, &              ! gravitational acceleration (~9.8 m/s2)
          cp                  ! specific heat of dry air at constant pressure at 273 K in J/kg/K

      real, intent(in), dimension(nx,nzm) ::   &   ! JUNG_localp
          pres     ! pressure (mb) at center of model levels
          
      real, intent(in), dimension(nx,nzm+1) :: &   ! JUNG_localp
          presi    ! pressure (mb) at model interfaces.
      
      real, intent(in), dimension(nx,nzm) :: &
          tabs, &  ! absolute temperature (K) at model levels
          qv, &    ! water vapor mass mixing ratio (kg/kg) 
          qcl, &   ! cloud liquid mass mixing ratio (kg/kg) 
          qci      ! cloud ice mass mixing ratio (kg/kg) 
      
      real, intent(in), dimension(nx) :: &
          tg, &              ! ground (or sea surface) temperature (K)
          latitude, &        ! latitude (degrees)
          longitude, &       ! longitude (degrees)
          p_factor_slice, &  ! perpetual factor
          p_coszrs_slice     ! perpetual cosine zenith angle

      real, intent(in), dimension(nzm) :: &
          o3, &
          co2, &
          ch4, &
          n2o, &
          o2, &
          cfc11, &
          cfc12, &
          cfc22, &
          ccl4

!------------------------------------------------------------------------------
! Output variables
!
! Note that fluxes are located at interfaces and have an extra 
! level.  The top fluxes, e.g. lwup(:,nzm+2), are approximate
! top-of-atmosphere fluxes.  lwup(:,nzm+1) is the top-of-model flux.
!------------------------------------------------------------------------------

    real(kind=kind_rb), intent(out), dimension(nx,nzm+2) :: &
         lwUp, &            ! upward longwave radiative flux (W/m2)
         lwDown, &          ! downward longwave radiative flux (W/m2)
         lwUpClearSky, &    ! clearsky upward longwave radiative flux (W/m2)
         lwDownClearSky, &  ! clearsky downward longwave radiative flux (W/m2)
         swUp, &            ! upward shortwave radiative flux (W/m2)
         swDown, &          ! downward shortwave radiative flux (W/m2)
         swUpClearSky, &    ! clearsky upward shortwave radiative flux (W/m2)
         swDownClearSky     ! clearsky downward shortwave radiative flux (W/m2)

! Heating rate outputs in K/day
      real (kind=kind_rb), intent(out), dimension(nx,nzm+1) :: &
          swHeatingRate, &
          swHeatingRateClearSky, &
          lwHeatingRate, &
          lwHeatingRateClearSky

      real, intent(out) :: &
          coszrs     ! cosine solar zenith angle

      real (kind = kind_rb), intent(out), dimension(nx, nzm+1) ::     &
          LWP, &        ! liquid water path (g/m2)
          IWP, &        ! ice water path (g/m2)
          liquidRe, &   ! effective radius liquid (microns)
          iceRe         ! effective radius ice (microns)
        
!------------------------------------------------------------------------------
! Local variables 
!------------------------------------------------------------------------------

! Input and output variables for RRTMG SW and LW
!   RRTM specifies the kind of real variables in 
!   Only one column dimension is allowed parkind
!   RRTM is indexed from bottom to top

!bloss: add layer to top to improve top-of-model heating rates.

    real (kind = kind_rb), dimension(nx, nzm+1) ::     &
        layerP, &     ! model layer pressure (mb)
        layerT, &     ! model layer temperature (K)
        layerMass, &  ! layer mass is for convenience
        
        cloudFrac, &  ! cloud fraction
        h2ovmr, &     ! volume mixing ratio H2O
        o3vmr, &      ! volume mixing ratio ozone
        co2vmr, &     ! volume mixing ratio CO2
        ch4vmr, &     ! volume mixing ratio CH4
        n2ovmr, &     ! volume mixing ratio N20
        o2vmr, &      ! volume mixing ratio O2
        cfc11vmr, &   ! volume mixing ratio CFC-11
        cfc12vmr, &   ! volume mixing ratio CFC-12
        cfc22vmr, &   ! volume mixing ratio CFC-22
        ccl4vmr       ! volume mixing ratio CCl4

    real(kind = kind_rb), dimension(nx, nzm+2) :: &
        interfaceP, &    ! Interface pressure (mb)
        interfaceT       ! Interface temperature (K)

      real(kind = kind_rb), dimension(nx) :: &
          surfaceT, &
          solarZenithAngleCos

      real, dimension(nx) ::  &
          asdir, &   ! Surface direct shortwave albedo (0.2-0.7 microns)
          asdif, &   ! Surface diffuse shortwave albedo (0.2-0.7 microns)
          aldir, &   ! Surface direct longwave albedo (0.7-5.0 microns)
          aldif      ! Surface diffuse longwave albedo (0.7-5.0 microns)

      real(kind = kind_rb), dimension(nx, nbndlw) :: &
          surfaceEmissivity    ! default = 0.95

      integer :: &
          overlap   ! RRTMG cloud overlap method (default = 1)

      integer :: i,k
      real    :: dayForSW
   
   ! JUNG
   !  real :: delta
      real(kind = kind_rb) :: delta,eccf_tmp

!bloss: add layer to top to improve top-of-model heating rates.
      
      real(kind = kind_rb), dimension(nbndlw, nx, nzm+2) :: &
          dummyCloudPropsLW     ! Cloud properties for LW calculation (default = 0.) 
      
      real(kind = kind_rb), dimension(nbndsw, nx, nzm+2) :: &
          dummyCloudPropsSW, &  ! Cloud properties for SW calculation (default = 0.)
          dummyTauCloudSW       ! In-cloud optical depth for SW calculation (default = 0.)
      
      real(kind = kind_rb), dimension(nx, nzm+1, nbndlw) :: &
          dummyTauAerosolLW     ! Aerosol properties for LW calculation (default  = 0.) 
      
      real(kind = kind_rb), dimension(nx, nzm+1, nbndsw) :: &
          dummyAerosolProps     ! Aerosol properties for SW calculation (default = 0.) 

! Aerosol properties for iaer = 6
      real(kind = kind_rb), dimension(nx, nzm+1, naerec) :: &
          dummyAerosolProps2    ! Aerosol properties for SW calculation (default = 0.)

!==============================================================================
! Initialize some arrays/parameters used by RRTMG

    overlap = 1 ! specifies RRTMG overlap assumption

    surfaceEmissivity = 0.95 ! Default Value For Sea Surface emissivity

! cloud properties are specified through LWP, IWP, liquidRe, iceRe, cloudFrac
    dummyCloudPropsLW = 0.
    dummyCloudPropsSW = 0.
    dummyTauCloudSW = 0.

! no aerosols in this interface
    dummyTauAerosolLW = 0.
    dummyAerosolProps = 0.
    dummyAerosolProps2 = 0.

! initialize coszrs to non-physical value as a placeholder
    coszrs = -2.

! Fill out 2D arrays needed by RRTMG 
!----------------------------------------------- JUNG_localp
!    layerP(:, 1:nzm) = spread(pres (:), dim = 1, ncopies = nx) 
!    layerP(:, nzm+1) = 0.5*spread(presi(nzm+1), dim = 1, ncopies = nx) ! add layer
!    interfaceP(:, 1:nzm+1) = spread(presi(:), dim = 1, ncopies = nx) 
!    interfaceP(:, nzm+2) = MIN(1.e-4_kind_rb,0.25*layerP(1,nzm+1)) ! near-zero pressure at top of extra layer

    layerP(:, 1:nzm) = pres (:, 1:nzm)
    layerP(:, nzm+1) = 0.5*presi(:, nzm+1) ! add layer
    interfaceP(:, 1:nzm+1) = presi(:, 1:nzm+1)
    interfaceP(:, nzm+2) = MIN(1.e-4_kind_rb,0.25*layerP(:,nzm+1)) ! near-zero pressure at top of extra layer
!-----------------------------------------------

! Convert hPa to Pa in layer mass calculation (kg/m2) 
    layerMass(:, 1:nzm+1) = &
         100. * (interfaceP(:,1:nzm+1) - interfaceP(:,2:nzm+2))/ ggr

! Set up for the radiation computation.
    lwHeatingRate(:, :) = 0.; swHeatingRate(:, :) = 0. 
    lwHeatingRateClearSky(:, :) = 0.; swHeatingRateClearSky(:, :) = 0. 

    layerT(:, 1:nzm) = tabs(:, 1:nzm) 
    layerT(:, nzm+1) = 2.*tabs(:, nzm) - tabs(:, nzm-1) ! add a layer at top.
      
! interpolate to find interface temperatures.
    interfaceT(:, 2:nzm+1) = (layerT(:, 1:nzm) + layerT(:, 2:nzm+1)) / 2. 

! Extrapolate temperature at top from lapse rate within the layer
    interfaceT(:, nzm+2) = 2.*layerT(:, nzm+1) - interfaceT(:, nzm+1)

! Use SST as interface temperature of atmosphere at surface.
    interfaceT(:, 1)  = tg(1:nx) !bloss layerT(:, 1)   + (layerT(:, 1)   - interfaceT(:, 2))   

!------------------------------------------------------------------------------
! Compute cloud IWP/LWP and particle sizes - convert from kg to g

    LWP(:, 1:nzm) = qcl(:, 1:nzm) * 1.e3 * layerMass(:, 1:nzm) 
    LWP(:, nzm+1) = 0. ! zero out extra layer

    IWP(:, 1:nzm) = qci(:, 1:nzm) * 1.e3 * layerMass(:, 1:nzm) 
    IWP(:, nzm+1) = 0. ! zero out extra layer

!bloss(072109): Previous implementation (using where/elsewhere)
!   had the undesirable effect of zeroing out cloudFrac
!   where LWP>0 but IWP==0.

    cloudFrac(:, :) = 0.
    liquidRe(:, :) = 0.
    iceRe(:, :) = 0.

    where(LWP(:, :) > 0.)
      liquidRe(:, :) = computeRe_Liquid(real(layerT), merge(0., 1., ocean))
      cloudFrac(:, :) = 1. 
    end where

    where(IWP(:, :) > 0.)
      !bloss: limit within RRTMG bounds on valid effective radii
      iceRe(:, :) = MAX(5., MIN(140., computeRe_Ice(real(layerT)) ) ) 
      cloudFrac(:, :) = 1. 
    end where

!------------------------------------------------------------------------------
! Volume mixing fractions for gases.

!bloss(072009): Note that o3, etc. are now in ppmv and do not need conversions.

    h2ovmr(:, 1:nzm)   = mwdry/mwh2o * qv(:, 1:nzm) 
    h2ovmr(:, nzm+1)   = mwdry/mwh2o * qv(:, nzm) ! extrapolate into extra layer at top
    o3vmr(:, 1:nzm)    = mwdry/mwo3 * spread(o3(:), dim = 1, ncopies = nx) 
    co2vmr(:, 1:nzm)   = mwdry/mwco2 * spread(co2(:), dim = 1, ncopies = nx) 
    ch4vmr(:, 1:nzm)   = mwdry/mwch4 * spread(ch4(:), dim = 1, ncopies = nx) 
    n2ovmr(:, 1:nzm)   = mwdry/mwn2o * spread(n2o(:), dim = 1, ncopies = nx) 
    o2vmr(:, 1:nzm)    = mwdry/mwo2 * spread(o2(:), dim = 1, ncopies = nx) 
    cfc11vmr(:, 1:nzm) = mwdry/mwf11 * spread(cfc11(:), dim = 1, ncopies = nx) 
    cfc12vmr(:, 1:nzm) = mwdry/mwf12 * spread(cfc12(:), dim = 1, ncopies = nx) 
    cfc22vmr(:, 1:nzm) = mwdry/mwf22 * spread(cfc22(:), dim = 1, ncopies = nx) 
    ccl4vmr(:, 1:nzm)  = mwdry/mwccl4 * spread(ccl4(:), dim = 1, ncopies = nx)
    
!tcram: extrapolate into extra layer at top
    o3vmr(:, nzm+1)    = mwdry/mwo3 * spread(o3(nzm), dim = 1, ncopies = nx) 
    co2vmr(:, nzm+1)   = mwdry/mwco2 * spread(co2(nzm), dim = 1, ncopies = nx)
    ch4vmr(:, nzm+1)   = mwdry/mwch4 * spread(ch4(nzm), dim = 1, ncopies = nx)
    n2ovmr(:, nzm+1)   = mwdry/mwn2o * spread(n2o(nzm), dim = 1, ncopies = nx) 
    o2vmr(:, nzm+1)    = mwdry/mwo2 * spread(o2(nzm), dim = 1, ncopies = nx) 
    cfc11vmr(:, nzm+1) = mwdry/mwf11 * spread(cfc11(nzm), dim = 1, ncopies = nx) 
    cfc12vmr(:, nzm+1) = mwdry/mwf12 * spread(cfc12(nzm), dim = 1, ncopies = nx) 
    cfc22vmr(:, nzm+1) = mwdry/mwf22 * spread(cfc22(nzm), dim = 1, ncopies = nx) 
    ccl4vmr(:, nzm+1)  = mwdry/mwccl4 * spread(ccl4(nzm), dim = 1, ncopies = nx)

!------------------------------------------------------------------------------
! Longwave calculation
!------------------------------------------------------------------------------

    if (dolongwave) then
    
!JUNG      if(lat.eq.1.AND.masterproc) print *, "Computing longwave radiation ..." 
      
      surfaceT(:) = tg(1:nx)

      call rrtmg_lw (nx, nzm+1, overlap,                    & 
           layerP, interfaceP, layerT, interfaceT, surfaceT, &
           h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
           cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, surfaceEmissivity,  &
           2, 3, 1, cloudFrac, &
           dummyCloudPropsLW, IWP, LWP, iceRe, liquidRe, &
           dummyTauAerosolLW, &
           lwUp,lwDown, lwHeatingRate, lwUpClearSky, lwDownClearSky, lwHeatingRateClearSky)

    end if !if(dolongwave)

!------------------------------------------------------------------------------
! End longwave calculation
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Shortwave calculation
!------------------------------------------------------------------------------

    if(doshortwave) then

!JUNG        if(lat.eq.1.AND.masterproc) print *, "Computing shortwave radiation ..." 

! Solar insolation depends on several choices

!---------------
      if(doperpetual) then
        if(dosolarconstant) then
        ! Fix solar constant and zenith angle
          solarZenithAngleCos(:) = cos(zenith_angle * pi/180.)
          eccf = solar_constant/scon
        else          
!---------------
         ! perpetual sun (diurnally-averaged insolation w/insolation-weighted coszrs)
          solarZenithAngleCos(:) = p_coszrs_slice(:)
         
         ! Adjust insolation using the eccentricity factor
!tcram          eccf = p_factor/MAX(p_coszrs, EPSILON(p_coszrs))
          eccf = sum(p_factor_slice(:)/solarZenithAngleCos(:)) / real(nx) 
        end if
      else

!---------------
! diurnally-varying insolation

        if(doseasons) then 
          ! The diurnal cycle of insolation will vary
          ! according to time of year of the current day.
          dayForSW = day 
        else
          ! The diurnal cycle of insolation from the calendar
          ! day on which the simulation starts (day0) will be
          ! repeated throughout the simulation.
          dayForSW = float(floor(day0)) + day - float(floor(day))
        end if

! JUNG -----------------------------------------
!       call shr_orb_decl (dayForSW, eccen, mvelpp, lambm0, obliqr, delta, eccf)
        call shr_orb_decl (real(dayForSW,kind_rb),real(eccen,kind_rb),real(mvelpp,kind_rb), &
                           real(lambm0,kind_rb),real(obliqr,kind_rb),delta,eccf_tmp) 
        eccf = real(eccf_tmp,kind_rm)   
! JUNG -----------------------------------------                        

        ! JUNG: Add do-loop for function zenith because "zenith" is not pure anymore.
        !       Change due to the use of shared version of shr_orb_mode.           
        do i=1,nx 
        solarZenithAngleCos(i) =  &
             zenith(dayForSW, pi * latitude(i)/180., pi * longitude(i)/180.)
        enddo
      end if

!---------------
! coszrs is found in params.f90 and used in the isccp simulator
      coszrs = max(0._kind_rb, solarZenithAngleCos(1))

! We only call the shortwave if the sun is above the horizon. 
!   We assume that the domain is small enough that the entire 
!   thing is either lit or shaded

      if(all(solarZenithAngleCos(:) >= tiny(solarZenithAngleCos))) then 

!JUNG        if(lat.eq.1.AND.masterproc) print *, "Let's do some shortwave" 

        call albedo(ocean, real(solarZenithAngleCos(:)), &
             asdir(:), aldir(:), asdif(:), aldif(:))

        call rrtmg_sw(nx, nzm+1, overlap,                     & 
             layerP, interfaceP, layerT, interfaceT, surfaceT, &
             h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,     &
             real(asdir, kind_rb), real(asdif, kind_rb), real(aldir, kind_rb), real(aldif, kind_rb), &
             solarZenithAngleCos, real(eccf, kind_rb), 0, real(scon, kind_rb),   &
             2, 3, 1, cloudFrac, &
             dummyTauCloudSW, dummyCloudPropsSW, dummyCloudPropsSW, dummyCloudPropsSW, &
             IWP, LWP, iceRe, liquidRe,  &
             dummyAerosolProps, dummyAerosolProps, dummyAerosolProps, dummyAerosolProps2, &
             swUp, swDown, swHeatingRate, swUpClearSky, swDownClearSky, swHeatingRateClearSky)

!JUNG        if(lat.eq.1.AND.masterproc) then
!JUNG          if(doshortwave) then 
!JUNG            if(doperpetual) then
!JUNG              write(*,992) coszrs, SUM(swDown(1:nx,nzm+2))/float(nx)
!JUNG            else
!JUNG              write(*,991) coszrs, SUM(swDown(1:nx,nzm+2))/float(nx), eccf
!JUNG            end if
!JUNG            991 format('radiation: diurnally-varying insolation, coszrs = ',F10.7, &
!JUNG                     ' solin = ',f10.4,' eccf = ',f10.7)
!JUNG            992 format('radiation: diurnally-averaged insolation, coszrs = ',F10.7, &
!JUNG                     ' solin = ',f10.4)
!JUNG          end if
!JUNG          write(*,993) asdir(1), aldir(1), asdif(1), aldif(1)
!JUNG993       format('radiation: surface albedos, asdir= ',F10.7, &
!JUNG               ' aldir = ',f10.7,' asdif = ',f10.7,' aldif = ',f10.7)
!JUNG        end if

      else ! if sun is down
        coszrs = 0.
        swUp = 0.
        swDown = 0.
        swUpClearSky = 0.
        swDownClearSky = 0.
        swHeatingRate = 0.
        swHeatingRateClearSky = 0.
      end if ! if sun is up
      
    end if ! if(doshortwave)
    
!------------------------------------------------------------------------------
! End shortwave calculation
!------------------------------------------------------------------------------

  end subroutine rad_driver_rrtm

!==============================================================================
  subroutine initialize_radiation(nx,ny,cp,iyear,day0,latitude,longitude,doperpetual)
!==============================================================================

    use parkind, only: kind_rb,kind_rm     ! kind_rm is added by JUNG
    use rrtmg_sw_init, only: rrtmg_sw_ini
    use rrtmg_lw_init, only: rrtmg_lw_ini
    
    implicit none

! INPUTS
 
    integer, intent(in) :: nx, ny      ! x- and y-grid sizes
    real, intent(in) :: cp             ! specific heat of dry air at constant pressure, J/kg/K
    integer, intent(in) :: iyear       ! year of simulation (for insolation computation)
    real, intent(in) :: day0           ! day of year during iyear (0.0 = 00Z Jan 1) at start of simulaton
    real, intent(in), dimension(:,:) :: latitude   ! latitude (degrees)
    real, intent(in), dimension(:,:) :: longitude  ! longitude (degrees)
    logical, intent(in) :: doperpetual             ! use perpetual (diurnally-averaged) insolation

! local variables
    real(KIND=kind_rb) :: cpdair
    integer :: ierr, i, j
    real :: p_factor, p_coszrs
    
!   JUNG    
    real(kind=kind_rb) :: eccen_tmp,obliq_tmp,mvelp_tmp
    real(kind=kind_rb) :: obliqr_tmp,lambm0_tmp,mvelpp_tmp

    !bloss  subroutine shr_orb_params
    !bloss  inputs:  iyear, log_print
    !bloss  ouptuts: eccen, obliq, mvelp, obliqr, lambm0, mvelpp

! JUNG ------------------------------------------
!   call shr_orb_params(iyear, eccen, obliq, mvelp, obliqr, lambm0, mvelpp, .false.)
    call shr_orb_params(iyear,eccen_tmp,obliq_tmp,mvelp_tmp,obliqr_tmp, &
                        lambm0_tmp,mvelpp_tmp,.false.)
    eccen = real(eccen_tmp,kind_rm)
    obliq = real(obliq_tmp,kind_rm)
    mvelp = real(mvelp_tmp,kind_rm)
    obliqr = real(obliqr_tmp,kind_rm)
    lambm0 = real(lambm0_tmp,kind_rm)
    mvelpp = real(mvelpp_tmp,kind_rm)
! JUNG ------------------------------------------                        
 
    if(doperpetual) then
      ! perpetual sun (no diurnal cycle)
      !   get diurnally-averaged insolation as a factor of solar constant
      !   as well as mean insolation-weighted solar zenith angle.
      
      do 100 j=1,ny
      do 100 i=1,nx
        call perpetual_factor_coszrs(day0, latitude(i,j), longitude(i,j), p_factor, p_coszrs)
        p_factor_xy(i,j) = p_factor
        p_coszrs_xy(i,j) = p_coszrs
 100  continue
      
    end if
 
    cpdair = cp
    call rrtmg_sw_ini(cpdair)
    call rrtmg_lw_ini(cpdair)
    
    isInitialized_RadDriver = .true. 
  end subroutine initialize_radiation

!==============================================================================
! JUNG  elemental real function zenith(calday, clat, clon)
!       Change due to the use of shared version of shr_orb_mode.  
!       Since shr_orb_decl is not pure in the shared version, 
!       thus, zenith cannot be pure.       
  function zenith(calday, clat, clon) result(zenith_result)
!==============================================================================

! ----------------------------------------------------------------------------
! Astronomy-related procedures 
! ----------------------------------------------------------------------------

     use shr_orb_mod, only : shr_orb_decl, shr_orb_cosz
     
     use parkind, only : kind_rb,kind_rm   ! added by JUNG
                          
     implicit none
     real, intent(in) :: calday, & ! Calendar day, including fraction
                         clat,   & ! Current centered latitude (radians)
                         clon      ! Centered longitude (radians)
     real zenith_result
! JUNG ------------------------------------------------------    
!    real :: delta,eccf  ! delta: Solar declination angle in radians
     real(kind=kind_rb) :: delta,eccf,zenith_tmp  ! delta: Solar declination angle in radians

!    call shr_orb_decl (calday, eccen, mvelpp, lambm0, obliqr, delta, eccf)
!    Compute local cosine solar zenith angle
!    zenith = shr_orb_cosz(calday, clat, clon, delta)

     call shr_orb_decl (real(calday,kind_rb),real(eccen,kind_rb),real(mvelpp,kind_rb), &
                        real(lambm0,kind_rb),real(obliqr,kind_rb), &
                        delta,eccf)

!    Compute local cosine solar zenith angle
     zenith_tmp = shr_orb_cosz(real(calday,kind_rb),real(clat,kind_rb), &
                               real(clon,kind_rb),delta)
     zenith_result = real(zenith_tmp,kind_rm)
!------------------------------------------------------------     
     
  end function zenith

!==============================================================================
  subroutine perpetual_factor_coszrs(day, lat, lon, p_factor, p_coszrs)
  
!  estimate the factor to multiply the solar constant
!  so that the sun hanging perpetually right above
!  the head (zenith angle=0) would produce the same
!  total input the TOA as the sun subgect to diurnal cycle.
!  coded by Marat Khairoutdinov, 2004  
!==============================================================================
     use parkind, only : kind_rb,kind_rm   ! added by JUNG
     use shr_orb_mod, only : shr_orb_decl
     implicit none

     real, intent(in)  :: day, lat, lon ! Day (without fraction); centered lat/lon (degrees) 
     real, intent(out) :: p_factor, p_coszrs

! Local:
     
! JUNG --------------------------------------------------------------------
!    real :: delta,eccf    ! delta: Solar declination angle in radians
     real :: eccf          
     real(kind=kind_rb) :: delta,eccf_tmp  ! delta: Solar declination angle in radians
! JUNG --------------------------------------------------------------------     
          
     integer :: n
     real :: tmp, tmp2, dttime, ttime 
     real :: coszrs
     real :: clat, clon   ! latitude and longitude converted to radians
    
     real, parameter :: dtrad = 60. ! default 60 second interval between insolation computations

    tmp = 0.
    tmp2 = 0.

    clat = pi * lat/180.
    clon = pi * lon/180.

    do n = 1,60*24 ! compute insolation for each minute of the day.

      ttime = day+float(n)*dtrad/86400.

! JUNG ------------------------------------------------------------------
!     call shr_orb_decl (ttime,eccen,mvelpp,lambm0,obliqr,delta,eccf)
      call shr_orb_decl (real(ttime,kind_rb),real(eccen,kind_rb),real(mvelpp,kind_rb), &
                         real(lambm0,kind_rb),real(obliqr,kind_rb), &
                         delta,eccf_tmp)
      
      eccf = real(eccf_tmp,kind_rm) 
! JUNG ------------------------------------------------------------------                       
                         
      coszrs = zenith(ttime, clat, clon)
      tmp  = tmp  + max(0., eccf * coszrs)
      tmp2 = tmp2 + max(0., eccf * coszrs) * max(0.,coszrs)

    end do
    
    tmp = tmp/float(60*24) ! average insolation/scon across minutes
    tmp2 = tmp2/float(60*24) ! average insolation*coszrs/scon across minutes
    
    p_factor = tmp ! mean insolation divided by solar constant
    p_coszrs = tmp2/MAX(tmp, EPSILON(tmp)) ! mean insolation-weighted cosine of solar zenith angle
    
!!$    write(*,*) 'eccentricity factor = ', eccf
!!$    write(*,*) 'delta = ', delta
!!$    write(*,*) 'insolation factor = ', p_factor
!!$    write(*,*) 'insolation-weighted coszrs = ', p_coszrs
!!$    write(*,*) 'perpetual insolation = ', p_factor*scon

  end subroutine perpetual_factor_coszrs

!==============================================================================
!bloss
!bloss(081309): Eliminate restart facility for simplicity of interface.
!bloss          Above-model sounding and trace gas concentrations will
!bloss          be re-computed hourly during the run and at any restart.
!bloss          I realize that this makes consistency between runs with
!bloss          and without restarts difficult if the restarts occur within
!bloss          an hour, but this seems a small price to pay in return for
!bloss          a cleaner interface.
!bloss          Direct complaints to pblossey@gmail.com.
!bloss  ! ----------------------------------------------------------------------------
!bloss  !
!bloss  ! Writing and reading binary restart files
!bloss  !
!bloss  ! ----------------------------------------------------------------------------
!bloss  subroutine write_rad(masterproc,RestartFileName)
!bloss    use rad, only: npatch_start, npatch_end, nzpatch, psnd, tsnd, qsnd
!bloss    implicit none    
!bloss    logical, intent(in) :: masterproc
!bloss    character(LEN=250), intent(in) :: RestartFileName
!bloss
!bloss    if(masterproc) then
!bloss
!bloss      print*,'Writting radiation restart file...'
!bloss      open(56, file = trim(RestartFileName), status='unknown',form='unformatted')
!bloss      ! save trace gas soundings that have been interpolated onto model grid
!bloss      !   also, background sounding used to patch in T, q above model top (if necessary)
!bloss      write(56) o3, co2, ch4, n2o, o2, cfc11, cfc12, cfc22, ccl4, &
!bloss           npatch_start, npatch_end, nzpatch, psnd, tsnd, qsnd
!bloss      close(56)
!bloss
!bloss      print *,'Saved radiation restart file.'
!bloss    end if
!bloss
!bloss  end subroutine write_rad
!bloss  ! ----------------------------------------------------------------------------
!bloss  subroutine read_rad(masterproc,RestartFileName)
!bloss    use rad, only: npatch_start, npatch_end, nzpatch, psnd, tsnd, qsnd
!bloss    implicit none
!bloss    logical, intent(in) :: masterproc
!bloss    character(LEN=250), intent(in) :: RestartFileName
!bloss
!bloss    if(masterproc) print*,'Reading radiation restart file...'
!bloss
!bloss    open(56, file = trim(RestartFileName), status='unknown',form='unformatted')
!bloss    ! read trace gas soundings that have been interpolated onto model grid
!bloss    !   also, background sounding used to patch in T, q above model top (if necessary)
!bloss    read(56) o3, co2, ch4, n2o, o2, cfc11, cfc12, cfc22, ccl4, &
!bloss           npatch_start, npatch_end, nzpatch, psnd, tsnd, qsnd
!bloss    close(56)
!bloss      
!bloss  end subroutine read_rad      

end module rad_driver
