!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_short_wave_absorption_jerlov
!
!> \brief MPAS ocean tracer short wave
!> \author Doug Jacobsen
!> \date   12/17/12
!> \details
!>  This module contains the routine for computing
!>  short wave tendencies using Jerlov
!
!-----------------------------------------------------------------------

module ocn_tracer_short_wave_absorption_variable

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_forcing
   use mpas_stream_manager
   use ocn_constants
   use ocn_framework_forcing

   implicit none

   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_tracer_short_wave_absorption_variable_tend, &
             ocn_tracer_short_wave_absorption_variable_init, &
             ocn_get_variable_sw_fraction,                   &
             ocn_get_os00_coeffs,                            &
             ocn_init_shortwave_forcing_ohlmann,             &
             ocn_get_shortWaveData,                          &
             ocn_shortwave_forcing_write_restart

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_short_wave_absorption_jerlov_tend
!
!> \brief   Computes tendency term for surface fluxes
!> \author  Luke Van Roekel
!> \date    11/10/2015
!> \details
!>  This routine computes the tendency for tracers based on surface fluxes.
!>  This computation is now based on spatially variable chlorophyll, cloud fraction, and zenith angle
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_short_wave_absorption_variable_tend(meshPool, swForcingPool, forcingPool, index_temperature, & !{{{
                                  layerThickness, penetrativeTemperatureFlux, penetrativeTemperatureFluxOBL, tend, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------


      type (mpas_pool_type), intent(in) :: &
         meshPool,                         &          !< Input: mesh information
         swForcingPool,                    &          !< Input: chlorophyll, cloud, zenith data
         forcingPool

      real (kind=RKIND), dimension(:), intent(in) :: &
        penetrativeTemperatureFlux !< Input: penetrative temperature flux through the surface

      real (kind=RKIND), dimension(:), intent(out) :: &
        penetrativeTemperatureFluxOBL

      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness !< Input: Layer thicknesses

      integer, intent(in) :: index_temperature

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, k, depLev, nCells
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND), pointer :: config_surface_buoyancy_depth
      real (kind=RKIND) :: depth
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth
      real (kind=RKIND), dimension(:), allocatable :: weights
      real (kind=RKIND), dimension(:), pointer :: chlorophyllA, zenithAngle, clearSkyRadiation
      character (len=StrKIND), pointer :: config_sw_absorption_type
      real (kind=RKIND), dimension(4) :: Avals, Kvals
      real (kind=RKIND) :: cloudRatio ! cloud Ratio = 1 - incident_sfc_sw_radiation/clearSkyRadiation

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
      call mpas_pool_get_config(ocnConfigs, 'config_surface_buoyancy_depth', config_surface_buoyancy_depth)

      allocate(weights(nVertLevels+1))
      weights = 0.0_RKIND
      weights(1) = 1.0_RKIND
      Avals(:)=0.0_RKIND
      Kvals(:)=0.0_RKIND

      call mpas_pool_get_config(ocnConfigs, 'config_sw_absorption_type', config_sw_absorption_type)
      call mpas_pool_get_array(swForcingPool,'chlorophyllData',chlorophyllA)

      call mpas_pool_get_array(swForcingPool,'zenithAngle',zenithAngle)
      call mpas_pool_get_array(swForcingPool,'clearSkyRadiation',clearSkyRadiation)

      nCells = nCellsArray( 3 )

      !$omp do schedule(runtime) private(depth, k, cloudRatio, depLev)
      do iCell = 1, nCells
        depth = 0.0_RKIND
        cloudRatio = min(1.0_RKIND, 1.0_RKIND - penetrativeTemperatureFlux(iCell)/(hflux_factor*(1.0E-15_RKIND + &
                                                     clearSkyRadiation(iCell))))
        cloudRatio = max(0.0_RKIND, cloudRatio)

        call ocn_get_os00_coeffs(chlorophyllA(iCell),zenithAngle(iCell),cloudRatio,Avals, Kvals)

        do k = 1, maxLevelCell(iCell)
          depth = depth + layerThickness(k, iCell)

          call ocn_get_variable_sw_fraction(depth, weights(k+1), Avals, Kvals)
          tend(index_temperature, k, iCell) = tend(index_temperature, k, iCell) + penetrativeTemperatureFlux(iCell) &
                                            * (weights(k) - weights(k+1) )
        end do

        depth = 0.0_RKIND
        do k=1,maxLevelCell(iCell)
           depth = depth + layerThickness(k,iCell)
           if(depth > abs(config_surface_buoyancy_depth)) exit
        enddo

        if(k == maxLevelCell(iCell) .or. k == 1) then
           depLev=2
        else
           depLev=k
        endif
        penetrativeTemperatureFluxOBL(iCell)=penetrativeTemperatureFlux(iCell)*weights(depLev)

      end do
      !$omp end do

      deallocate(weights)

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_short_wave_absorption_variable_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_short_wave_absorption_variable_init
!
!> \brief   Initializes ocean tracer surface flux quantities
!> \author  Luke Van Roekel
!> \date    11/10/15
!> \details
!>  This routine initializes quantities related to surface fluxes in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_short_wave_absorption_variable_init(domain,err)!{{{

   !--------------------------------------------------------------------

      type(domain_type) :: domain

      integer, intent(out) :: err !< Output: error flag

      character (len=StrKIND), pointer :: config_sw_absorption_type

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_sw_absorption_type', config_sw_absorption_type)

      select case ( trim(config_sw_absorption_type) )
      case ('ohlmann00')
         call ocn_init_shortwave_forcing_ohlmann(domain)
      case default
         call mpas_log_write( &
            'Shortwave parameterization type unknown: config_sw_absortion_type=' // trim(config_sw_absorption_type) // &
            ' Options are: jerlov or ohlmann00 or none', &
            MPAS_LOG_CRIT)
      end select


   end subroutine ocn_tracer_short_wave_absorption_variable_init!}}}

!***********************************************************************

!***********************************************************************
!
!  routine ocn_init_shortwave_forcing_ohlmann
!
!> \brief   Initializes forcing group if parameterization of Ohlmann and Siegel (2000) or Ohlmann(2003)
!>          This parameterization only requires chlorophyll-a concentrations, so only add that
!> \author  Luke Van Roekel
!> \date    11/10/15
!> \details
!>  This routine initializes forcing stream for Ohlmann and Siegel (2000) or Ohlman (2003) parameterization
!
!-----------------------------------------------------------------------

   subroutine ocn_init_shortwave_forcing_ohlmann(domain)!{{{

   type(domain_type) :: domain

   logical, pointer :: &
        config_do_restart

   character(len=strKIND) :: &
        forcingIntervalMonthly,  &
        forcingReferenceTimeMonthly

   call MPAS_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)

   forcingIntervalMonthly = "0000-01-00_00:00:00"
   forcingReferenceTimeMonthly = "0000-01-01_00:00:00"

   call MPAS_forcing_init_group( forcingGroupHead,  &
        "shortwave_monthly_observations", &
        domain, &
        '0000-01-01_00:00:00', &
        '0000-01-01_00:00:00', &
        '0001-00-00_00:00:00', &
        config_do_restart)

   call MPAS_forcing_init_field( domain % streamManager, &
        forcingGroupHead, &
        'shortwave_monthly_observations', &
        'chlorophyllData', &
        'shortwave_forcing_data', &
        'shortwave',  &
        'chlorophyllData',  &
        'constant',  &
        forcingReferenceTimeMonthly,  &
        forcingIntervalMonthly)

   call MPAS_forcing_init_field( domain % streamManager, &
        forcingGroupHead, &
        'shortwave_monthly_observations', &
        'clearSkyRadiation', &
        'shortwave_forcing_data', &
        'shortwave',  &
        'clearSkyRadiation',  &
        'constant',  &
        forcingReferenceTimeMonthly,  &
        forcingIntervalMonthly)

   call MPAS_forcing_init_field( domain % streamManager, &
        forcingGroupHead, &
        'shortwave_monthly_observations', &
        'zenithAngle', &
        'shortwave_forcing_data', &
        'shortwave',  &
        'zenithAngle',  &
        'constant',  &
        forcingReferenceTimeMonthly,  &
        forcingIntervalMonthly)

   call MPAS_forcing_init_field_data( forcingGroupHead, &
        'shortwave_monthly_observations', &
        domain % streamManager, &
        config_do_restart, &
        .false.)

   end subroutine ocn_init_shortwave_forcing_ohlmann!}}}

!***********************************************************************

!***********************************************************************
!
!  routine get_shortWaveData
!
!> \brief   retrieve data needed to compute penetration of shortwave radiation
!> \author  Luke Van Roekel
!> \date    11/10/15
!> \details
!>  This routine calls mpas_forcing routines to acquire needed shortwave data and interpolates
!>    between time levels
!
!-----------------------------------------------------------------------

    subroutine ocn_get_shortWaveData( streamManager, &
        domain, &
        simulationClock, &
        firstTimeStep) !{{{

        type (MPAS_streamManager_type), intent(inout) :: streamManager

        type (domain_type) :: domain
        type (MPAS_timeInterval_type) :: timeStepSW
        type (MPAS_clock_type) :: simulationClock

        logical,pointer :: config_use_activeTracers_surface_bulk_forcing
        logical, intent(in) :: firstTimeStep
        character(len=strKind), pointer :: config_sw_absorption_type
        character(len=strKind), pointer :: config_dt
        real(kind=RKIND) :: dt


        call MPAS_pool_get_config(domain%configs, 'config_use_activeTracers_surface_bulk_forcing', &
                                  config_use_activeTracers_surface_bulk_forcing)
        call MPAS_pool_get_config(domain%configs, 'config_sw_absorption_type', config_sw_absorption_type)
        call MPAS_pool_get_config(domain%configs, 'config_dt', config_dt)

        call mpas_set_timeInterval(timeStepSW,timeString=config_dt)
        call mpas_get_timeInterval(timeStepSW,dt=dt)

        if(trim(config_sw_absorption_type) == 'ohlmann00' .and. config_use_activeTracers_surface_bulk_forcing) then
             call MPAS_forcing_get_forcing(forcingGroupHead, &
                  'shortwave_monthly_observations', streamManager, dt)
        endif

    end subroutine ocn_get_shortWaveData!}}}


!***********************************************************************


!***********************************************************************
!
!  routine ocn_get_variable_fractions
!
!> \brief   Computes short wave absorption fractions
!> \author  Luke Van Roekel
!> \date    11/10/2015
!> \details
!>  Computes fraction of solar short-wave flux penetrating to
!>  specified depth due to time and space varying chlorophyll, cloud fraction, and zenith angle
!> based on:
!>     Ohlmann and Siegel (2000), Ohlmann (2003), Manizza et al. (2005)

!
!-----------------------------------------------------------------------
   subroutine ocn_get_variable_sw_fraction(depth, weight, Avals, Kvals)!{{{
!  Note: below 200m the solar penetration gets set to zero,
!     otherwise the limit for the exponent ($+/- 5678$) needs to be
!     taken care of.

      real (kind=RKIND), intent(in) :: depth !< Input: Depth of bottom of cell
      real (kind=RKIND), intent(in),dimension(4) :: Avals  !< Input: spectral partitioning of radiation
      real (kind=RKIND), intent(in),dimension(4) :: Kvals  !< Input: extinction coefficients for different radiation bands
      real (kind=RKIND), intent(out) :: weight !< Output: Weight for Jerlov absorption

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!
      integer :: k

      real (kind=RKIND), parameter :: depth_cutoff = -200.0_RKIND

!-----------------------------------------------------------------------
!
!  compute absorption fraction
!
!-----------------------------------------------------------------------

      if (-depth < depth_cutoff) then
         weight = 0.0_RKIND
      else
         weight=0.0_RKIND
         do k=1,4
            weight = weight + Avals(k)*exp(-depth*Kvals(k))
         enddo
      endif
   end subroutine ocn_get_variable_sw_fraction!}}}

!***********************************************************************



!***********************************************************************
!
!  routine ocn_get_os00_coeffs
!
!> \brief   Computes coefficients for spatially varying penetrating shortwave via Ohlmann and Siegel 2000
!> \author  Luke Van Roekel
!> \date    11/10/2015
!> \details
!>   This will fill in coefficients for the penetrating shortwave parameterization
!>      A1*exp(-K1*depth) + A2*exp(-K2*depth) + A3*exp(-K3*depth) + A4*exp(-K4*depth)
!>   For Ohlmann and Siegel (2000) A4 no approximations are made
!>   Here the IR portion of the spectrum is also decomposed.  Four exponential terms are used
!>   NOTE: other schemes can be easily recovered.  For example, Ohlmann (2003) requires a
!>      coefficient changes only Use equations 6a - 6d and 7
!            Ohlmann, JC, 2003: Ocean Radiant Heating in Climate Models, J.Clim, v16, 1337-1351
!     To recover the Morel (1988) scheme used in Manizza et al. (2005), use equations 1-4 of
!           Manizza, M, C LeQuere, AJ Watson, ET Buitenhuis, 2005: Bio-optical feedbacks among phytoplankton
!              upper ocean physics and sea-ice in a global model, Geophys, Res. Lett
!
!-----------------------------------------------------------------------

  subroutine ocn_get_os00_coeffs(chlorophyllA,zenithAngle,cloudFraction,Avals,Kvals)!{{{

     real(kind=RKIND), intent(in) :: chlorophyllA, zenithAngle, cloudFraction
     real(kind=RKIND), intent(out), dimension(4) :: Avals, Kvals

     if(cloudFraction > 0.1_RKIND) then ! cloudy skies
         Avals(1) = 0.026_RKIND*chlorophyllA + 0.112_RKIND*cloudFraction + 0.366_RKIND
         Avals(2) = -0.009_RKIND*chlorophyllA + 0.034_RKIND*cloudFraction + 0.207_RKIND
         Avals(3) = -0.015_RKIND*chlorophyllA -0.006_RKIND*cloudFraction + 0.188_RKIND
         Avals(4) = -0.003_RKIND*chlorophyllA -0.131_RKIND*cloudFraction + 0.169_RKIND
         Kvals(1) = 0.063_RKIND*chlorophyllA -0.015_RKIND*cloudFraction + 0.082_RKIND
         Kvals(2) = 0.278_RKIND*chlorophyllA -0.562_RKIND*cloudFraction + 1.02_RKIND
         Kvals(3) = 3.91_RKIND*chlorophyllA -12.91_RKIND*cloudFraction + 16.62_RKIND
         Kvals(4) = 16.64_RKIND*chlorophyllA -478.28_RKIND*cloudFraction + 736.56_RKIND
     else ! clear skies
         Avals(1) = 0.033_RKIND*chlorophyllA -0.025_RKIND*zenithAngle + 0.419_RKIND
         Avals(2) = -0.010_RKIND*chlorophyllA -0.007_RKIND*zenithAngle + 0.231_RKIND
         Avals(3) = -0.019_RKIND*chlorophyllA -0.003_RKIND*zenithAngle + 0.195_RKIND
         Avals(4) = -0.006_RKIND*chlorophyllA -0.004_RKIND*zenithAngle + 0.154_RKIND
         Kvals(1) = 0.066_RKIND*chlorophyllA + 0.006_RKIND*zenithAngle + 0.066_RKIND
         Kvals(2) = 0.396_RKIND*chlorophyllA -0.027_RKIND*zenithAngle + 0.866_RKIND
         Kvals(3) = 7.68_RKIND*chlorophyllA -2.49_RKIND*zenithAngle + 17.81_RKIND
         Kvals(4) = 51.27_RKIND*chlorophyllA + 13.14_RKIND*zenithAngle + 665.19_RKIND
     endif

   end subroutine ocn_get_os00_coeffs!}}}

!***********************************************************************


!***********************************************************************
!
!  routine ocn_shortwave_forcing_write_restart
!
!> \brief   writes restart timestamp for SW data to be read in on future restart
!> \author  Luke Van Roekel
!> \date    11/16/2015

!
!-----------------------------------------------------------------------

   subroutine ocn_shortwave_forcing_write_restart(domain)!{{{

      type(domain_type) :: domain

      character(len=strKind), pointer :: config_sw_absorption_type


      call MPAS_pool_get_config(domain % configs, "config_sw_absorption_type", config_sw_absorption_type)

      if( trim(config_sw_absorption_type) == 'ohlmann00'  ) then
          call MPAS_forcing_write_restart_times(forcingGroupHead)
      endif

    end subroutine ocn_shortwave_forcing_write_restart!}}}


end module ocn_tracer_short_wave_absorption_variable


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
