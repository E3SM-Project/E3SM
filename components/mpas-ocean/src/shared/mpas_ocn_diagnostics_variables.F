! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_diagnostics_variables
!
!> \brief MPAS ocean diagnostics variable definitions
!> \author Matt Turner
!> \date   13 May 2020
!> \details
!>  This module contains the definitions of all variables contianed in the
!>  diagnosticsPool, and also gets the variables from the Pool in init.
!
!-----------------------------------------------------------------------

module ocn_diagnostics_variables

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_stream_manager
   use mpas_io_units

   use ocn_config

   implicit none
   public
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   real (kind=RKIND), dimension(:,:), pointer :: zMid
   real (kind=RKIND), dimension(:,:), pointer :: zTop
   real (kind=RKIND), dimension(:,:), pointer :: divergence
   real (kind=RKIND), dimension(:,:), pointer :: circulation
   real (kind=RKIND), dimension(:,:), pointer :: relativeVorticity
   real (kind=RKIND), dimension(:,:), pointer :: relativeVorticityCell
   real (kind=RKIND), dimension(:,:), pointer :: normPlanetVortEdge
   real (kind=RKIND), dimension(:,:), pointer :: normRelVortEdge
   real (kind=RKIND), dimension(:,:), pointer :: normalizedRelativeVorticityCell
   real (kind=RKIND), dimension(:,:), pointer :: density
   real (kind=RKIND), dimension(:,:), pointer :: displaceddensity
   real (kind=RKIND), dimension(:,:), pointer :: potentialdensity
   real (kind=RKIND), dimension(:,:), pointer :: montgomeryPotential
   real (kind=RKIND), dimension(:,:), pointer :: pressure
   real (kind=RKIND), dimension(:,:), pointer :: thermExpCoeff
   real (kind=RKIND), dimension(:,:), pointer :: salineContractCoeff
   real (kind=RKIND), dimension(:,:), pointer :: BruntVaisalaFreqTop
   real (kind=RKIND), dimension(:,:), pointer :: tangentialVelocity
   real (kind=RKIND), dimension(:,:), pointer :: layerThickEdgeFlux
   real (kind=RKIND), dimension(:,:), pointer :: layerThickEdgeMean
   real (kind=RKIND), dimension(:,:), pointer :: kineticEnergyCell
   real (kind=RKIND), dimension(:,:), pointer :: vertVelocityTop
   real (kind=RKIND), dimension(:,:), pointer :: vertTransportVelocityTop
   real (kind=RKIND), dimension(:,:), pointer :: vertGMBolusVelocityTop
   real (kind=RKIND), dimension(:,:), pointer :: normalGMBolusVelocity
   real (kind=RKIND), dimension(:,:), pointer :: normalTransportVelocity
   real (kind=RKIND), dimension(:,:), pointer :: RiTopOfCell
   real (kind=RKIND), dimension(:,:), pointer :: RiTopOfEdge

   real (kind=RKIND), dimension(:), pointer :: gradSSH
   real (kind=RKIND), dimension(:), pointer :: pressureAdjustedSSH

   real (kind=RKIND), dimension(:,:), pointer :: tracersSurfaceValue
   real (kind=RKIND), dimension(:,:), pointer :: tracersSurfaceLayerValue
   real (kind=RKIND), dimension(:),   pointer :: normalVelocitySurfaceLayer
   real (kind=RKIND), dimension(:),   pointer :: indexSurfaceLayerDepth

   real (kind=RKIND), dimension(:), pointer :: sfcFlxAttCoeff
   real (kind=RKIND), dimension(:), pointer :: surfaceFluxAttenuationCoefficientRunoff

   real (kind=RKIND), dimension(:), pointer :: landIceFrictionVelocity
   real (kind=RKIND), dimension(:), pointer :: topDrag
   real (kind=RKIND), dimension(:), pointer :: topDragMag
   real (kind=RKIND), dimension(:), pointer :: topographic_wave_drag
   real (kind=RKIND), dimension(:,:), pointer :: landIceBoundaryLayerTracers
   real (kind=RKIND), dimension(:,:), pointer :: landIceTracerTransferVelocities

   integer, pointer :: indexBLT, indexBLS, indexHeatTrans, indexSaltTrans

   real (kind=RKIND), dimension(:), pointer :: penetrativeTemperatureFluxOBL
   real (kind=RKIND), dimension(:), pointer :: bottomLayerShortwaveTemperatureFlux

   real (kind=RKIND), dimension(:), pointer :: surfaceBuoyancyForcing
   real (kind=RKIND), dimension(:), pointer :: surfaceFrictionVelocity
   real (kind=RKIND), dimension(:), pointer :: surfacePressure

   real (kind=RKIND), dimension(:,:), pointer :: &
      transportVelocityX, transportVelocityY, transportVelocityZ, transportVelocityZonal, &
      transportVelocityMeridional, GMBolusVelocityX, GMBolusVelocityY, &
      GMBolusVelocityZ, GMBolusVelocityZonal, GMBolusVelocityMeridional, gmStreamFuncTopOfEdge, &
      GMStreamFuncX, GMStreamFuncY, GMStreamFuncZ, GMStreamFuncZonal, GMStreamFuncMeridional, &
      gmStreamFuncTopOfCell

   real (kind=RKIND), dimension(:,:), pointer :: rx1Edge
   real (kind=RKIND), dimension(:,:), pointer :: rx1Cell
   real (kind=RKIND), dimension(:), pointer :: rx1MaxEdge
   real (kind=RKIND), dimension(:), pointer :: rx1MaxCell
   real (kind=RKIND), pointer :: globalRx1Max

   integer, dimension(:), pointer :: smoothingMask
   real (kind=RKIND), dimension(:,:), pointer :: verticalStretch

   integer, dimension(:), pointer :: modifyLandIcePressureMask

   real (kind=RKIND), dimension(:,:), pointer :: velocityMeridional
   real (kind=RKIND), dimension(:,:), pointer :: velocityZonal

   real(kind=RKIND), dimension(:,:), pointer :: gmKappaScaling
   real(kind=RKIND), dimension(:,:), pointer :: RediKappaScaling
   real(kind=RKIND), dimension(:,:), pointer :: RediKappaSfcTaper
   real(kind=RKIND), dimension(:,:), pointer :: k33
   real(kind=RKIND), dimension(:), pointer   :: gmBolusKappa
   real(kind=RKIND), dimension(:), pointer   :: cGMphaseSpeed
   real(kind=RKIND), dimension(:,:,:), pointer :: slopeTriadUp, slopeTriadDown

   integer, dimension(:), pointer :: indMLD
   real (kind=RKIND), dimension(:,:), pointer :: velocityX, velocityY, velocityZ

   real (kind=RKIND), dimension(:,:), pointer :: vertAleTransportTop

   real (kind=RKIND), dimension(:,:), pointer :: viscosity
   real (kind=RKIND), dimension(:,:), pointer :: vertViscTopOfEdge
   real (kind=RKIND), dimension(:,:), pointer :: vertViscTopOfCell
   real (kind=RKIND), dimension(:,:), pointer :: wettingVelocity
   type (field2DReal), pointer :: wettingVelocityField

   real (kind=RKIND), dimension(:,:), pointer :: vertDiffTopOfCell
   real (kind=RKIND), dimension(:,:,:), pointer :: vertNonLocalFlux
   integer, dimension(:,:), pointer :: rediLimiterCount
   real (kind=RKIND), dimension(:,:,:), pointer :: activeTracerSurfaceFluxTendency
   real (kind=RKIND), dimension(:,:,:), pointer :: activeTracerNonLocalTendency
   real (kind=RKIND), dimension(:,:,:), pointer :: activeTracerHorMixTendency
   real (kind=RKIND), dimension(:,:), pointer :: temperatureShortWaveTendency
   real (kind=RKIND), dimension(:), pointer :: RediKappa

   ! pointers for tendencies used in diagnostic budget computation
   real (kind=RKIND), dimension(:,:,:), pointer, contiguous :: &
      activeTracerHorizontalAdvectionTendency, &
      activeTracerVerticalAdvectionTendency,   &
      activeTracerHorizontalAdvectionEdgeFlux, &
      activeTracerVerticalAdvectionTopFlux,    &
      activeTracerVertMixTendency

   real (kind=RKIND), pointer :: daysSinceStartOfSim
   character (len=StrKIND), pointer :: xtime, simulationStartTime

   real (kind=RKIND), dimension(:,:), pointer :: bulkRichardsonNumber
   real (kind=RKIND), dimension(:,:), pointer :: bulkRichardsonNumberBuoy
   real (kind=RKIND), dimension(:,:), pointer :: bulkRichardsonNumberShear

   real (kind=RKIND), dimension(:), pointer :: boundaryLayerDepth
   real (kind=RKIND), dimension(:), pointer :: boundaryLayerDepthSmooth
   real (kind=RKIND), dimension(:), pointer :: indexBoundaryLayerDepth
   integer, pointer :: index_vertNonLocalFluxTemp

   real (kind=RKIND), dimension(:), pointer :: salinitySurfaceRestoringTendency

   real (kind=RKIND), dimension(:,:), pointer :: surfaceVelocity
   real (kind=RKIND), dimension(:), pointer :: gradSSHZonal, gradSSHMeridional
   real (kind=RKIND), dimension(:), pointer :: gradSSHX
   real (kind=RKIND), dimension(:), pointer :: gradSSHY
   real (kind=RKIND), dimension(:), pointer :: gradSSHZ

   real (kind=RKIND), dimension(:), pointer :: barotropicForcing
   real (kind=RKIND), dimension(:), pointer :: barotropicCoriolisTerm
   real (kind=RKIND), dimension(:), pointer :: barotropicThicknessFlux

   real (kind=RKIND), dimension(:, :), pointer :: SSHGradient
   integer, pointer :: indexSurfaceVelocityZonal, indexSurfaceVelocityMeridional
   integer, pointer :: indexSSHGradientZonal, indexSSHGradientMeridional

   type (field1DInteger), pointer :: smoothingMaskField
   type (field2DReal), pointer :: verticalStretchField

   real (kind=RKIND), pointer :: globalVerticalStretchMax, globalVerticalStretchMin
   real (kind=RKIND), dimension(:), pointer :: betaEdge
   real (kind=RKIND), dimension(:), pointer :: eddyLength
   real (kind=RKIND), dimension(:), pointer :: eddyTime
   real (kind=RKIND), dimension(:), pointer :: gmResolutionTaper

   ! Semi-implicit Array Pointers
   real (kind=RKIND), dimension(:), pointer :: SIvec_r0,SIvec_r00,SIvec_r1 ,SIvec_rh0,SIvec_rh1,SIvec_ph0,SIvec_ph1
   real (kind=RKIND), dimension(:), pointer :: SIvec_v0,SIvec_s0 ,SIvec_s1 ,SIvec_sh0,SIvec_sh1
   real (kind=RKIND), dimension(:), pointer :: SIvec_t0,SIvec_t1 ,SIvec_q0 ,SIvec_q1 ,SIvec_qh0
   real (kind=RKIND), dimension(:), pointer :: SIvec_w0,SIvec_w1 ,SIvec_wh0,SIvec_wh1,SIvec_y0
   real (kind=RKIND), dimension(:), pointer :: SIvec_z0,SIvec_z1 ,SIvec_zh0,SIvec_zh1

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_diagnostics_variables_init
   public :: ocn_diagnostics_variables_destroy

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_diagnostics_variables_init
!
!> \brief   Initializes variables in the diagnosticsPool
!> \author  Matt Turner
!> \date    13 May 2020
!> \details
!>  This routine initializes all of the pointers associated with
!>  variables contained in the diagnosticsPool.
!
!-----------------------------------------------------------------------
   subroutine ocn_diagnostics_variables_init(domain, jenkinsOn, hollandJenkinsOn, err)!{{{
      type (domain_type), intent(in) :: domain

      logical, intent(in) :: jenkinsOn
      logical, intent(in) :: hollandjenkinsOn

      integer, intent(out) :: err !< Output: Error flag

      type (mpas_pool_type), pointer :: diagnosticsPool

      err = 0

      call mpas_pool_get_subpool(domain % blocklist % structs, 'diagnostics', diagnosticsPool)

      if (config_use_topographic_wave_drag) then
         call mpas_pool_get_array(diagnosticsPool, 'topographic_wave_drag', topographic_wave_drag)
      end if

      if ( trim(config_ocean_run_mode) == 'forward' ) then
         call mpas_pool_get_array(diagnosticsPool, 'wettingVelocity', &
                   wettingVelocity)
         call mpas_pool_get_field(diagnosticsPool, 'wettingVelocity', &
                   wettingVelocityField)
      end if

      if ( config_use_GM ) then
         call mpas_pool_get_array(diagnosticsPool, 'RediKappa', &
                   RediKappa)
         call mpas_pool_get_array(diagnosticsPool, 'RediKappaScaling', &
                   RediKappaScaling)
         call mpas_pool_get_array(diagnosticsPool, 'RediKappaSfcTaper', &
                   RediKappaSfcTaper)
         call mpas_pool_get_array(diagnosticsPool, 'rediLimiterCount', &
                   rediLimiterCount)

         call mpas_pool_get_array(diagnosticsPool, 'gmBolusKappa',  &
                   gmBolusKappa)
         call mpas_pool_get_array(diagnosticsPool, 'gmKappaScaling', &
                   gmKappaScaling)
         call mpas_pool_get_array(diagnosticsPool, 'gmResolutionTaper', &
                   gmResolutionTaper)
      end if

      if ( config_compute_active_tracer_budgets ) then
         call mpas_pool_get_array(diagnosticsPool, &
                  'activeTracerSurfaceFluxTendency', &
                   activeTracerSurfaceFluxTendency)
         call mpas_pool_get_array(diagnosticsPool, &
                  'temperatureShortWaveTendency', &
                   temperatureShortWaveTendency)
         call mpas_pool_get_array(diagnosticsPool, &
                  'activeTracerNonLocalTendency', &
                   activeTracerNonLocalTendency)
         call mpas_pool_get_array(diagnosticsPool, &
                  'activeTracerHorMixTendency', &
                   activeTracerHorMixTendency)
         call mpas_pool_get_array(diagnosticsPool, &
                  'activeTracerVerticalAdvectionTopFlux', &
                   activeTracerVerticalAdvectionTopFlux)
         call mpas_pool_get_array(diagnosticsPool, &
                  'activeTracerHorizontalAdvectionEdgeFlux', &
                   activeTracerHorizontalAdvectionEdgeFlux)
         call mpas_pool_get_array(diagnosticsPool, &
                  'activeTracerHorizontalAdvectionTendency', &
                   activeTracerHorizontalAdvectionTendency)
         call mpas_pool_get_array(diagnosticsPool, &
                  'activeTracerVerticalAdvectionTendency',  &
                   activeTracerVerticalAdvectionTendency)
         call mpas_pool_get_array(diagnosticsPool, &
                  'activeTracerVertMixTendency', &
                   activeTracerVertMixTendency)
      end if

      if ( config_use_activeTracers_surface_restoring ) then
         call mpas_pool_get_array(diagnosticsPool, &
                  'salinitySurfaceRestoringTendency', &
                   salinitySurfaceRestoringTendency)
      end if

      ! Retrieve pointers for configurations where land_ice_flux_mode is enabled
      if ( trim(config_land_ice_flux_mode) /= 'off' ) then
         call mpas_pool_get_array(diagnosticsPool, 'topDrag', topDrag)
         call mpas_pool_get_array(diagnosticsPool, 'topDragMagnitude', topDragMag)
         call mpas_pool_get_dimension(diagnosticsPool, &
                   'index_landIceBoundaryLayerTemperature', indexBLT)
         call mpas_pool_get_dimension(diagnosticsPool, &
                   'index_landIceBoundaryLayerSalinity', indexBLS)
         call mpas_pool_get_dimension(diagnosticsPool, &
                   'index_landIceHeatTransferVelocity', indexHeatTrans)
         call mpas_pool_get_dimension(diagnosticsPool, &
                   'index_landIceSaltTransferVelocity', indexSaltTrans)
         call mpas_pool_get_array(diagnosticsPool, &
                   'landIceBoundaryLayerTracers', landIceBoundaryLayerTracers)
         call mpas_pool_get_array(diagnosticsPool, &
                   'landIceTracerTransferVelocities', &
                    landIceTracerTransferVelocities)
         call mpas_pool_get_array(diagnosticsPool, &
                   'landIceFrictionVelocity', landIceFrictionVelocity)
      end if

      if ( trim(config_ocean_run_mode) == 'init' ) then
         call mpas_pool_get_array(diagnosticsPool, 'rx1Edge', rx1Edge, 1)
         call mpas_pool_get_array(diagnosticsPool, 'rx1Cell', rx1Cell, 1)
         call mpas_pool_get_array(diagnosticsPool, 'rx1MaxEdge', rx1MaxEdge, 1)
         call mpas_pool_get_array(diagnosticsPool, 'rx1InitSmoothingMask', &
                   smoothingMask, 1)
         call mpas_pool_get_array(diagnosticsPool, 'verticalStretch', &
                   verticalStretch, 1)
         call mpas_pool_get_array(diagnosticsPool, &
                  'globalVerticalStretchMax', &
                   globalVerticalStretchMax, 1)
         call mpas_pool_get_array(diagnosticsPool, &
                  'globalVerticalStretchMin', &
                   globalVerticalStretchMin, 1)
      end if

      if (config_use_topographic_wave_drag) then
         call mpas_pool_get_array(diagnosticsPool, &
                  'topographic_wave_drag', &
                   topographic_wave_drag)
      end if

      if ( trim(config_ocean_run_mode) == 'init' .or. &
                config_AM_debugDiagnostics_enable ) then
         call mpas_pool_get_array(diagnosticsPool, 'rx1MaxCell', &
                   rx1MaxCell, 1)
         call mpas_pool_get_array(diagnosticsPool, 'globalRx1Max', &
                   globalRx1Max, 1)
      end if

      if ( trim(config_ocean_run_mode) == 'forward' .or. &
           trim(config_ocean_run_mode) == 'analysis' ) then
         call mpas_pool_get_array(diagnosticsPool, 'k33', k33)

         call mpas_pool_get_array(diagnosticsPool, 'gradSSH', gradSSH)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHX', gradSSHX)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHY', gradSSHY)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHZ', gradSSHZ)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHZonal', &
                   gradSSHZonal)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHMeridional', &
                   gradSSHMeridional)
         call mpas_pool_get_array(diagnosticsPool, 'SSHGradient', SSHGradient)

         call mpas_pool_get_array(diagnosticsPool, 'viscosity', viscosity)
         call mpas_pool_get_array(diagnosticsPool, 'divergence', divergence)
         call mpas_pool_get_array(diagnosticsPool, 'circulation', circulation)

         call mpas_pool_get_array(diagnosticsPool, 'velocityX', velocityX)
         call mpas_pool_get_array(diagnosticsPool, 'velocityY', velocityY)
         call mpas_pool_get_array(diagnosticsPool, 'velocityZ', velocityZ)
         call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', &
                   velocityMeridional)
         call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', &
                   velocityZonal)

         call mpas_pool_get_array(diagnosticsPool, 'RiTopOfCell', RiTopOfCell)
         call mpas_pool_get_array(diagnosticsPool, 'RiTopOfEdge', RiTopOfEdge)

         call mpas_pool_get_array(diagnosticsPool, 'slopeTriadUp', &
                   slopeTriadUp)
         call mpas_pool_get_array(diagnosticsPool, 'slopeTriadDown', &
                   slopeTriadDown)

         ! GM-related fields
         call mpas_pool_get_array(diagnosticsPool, 'cGMphaseSpeed', &
                   cGMphaseSpeed)
         call mpas_pool_get_array(diagnosticsPool, 'gmStreamFuncTopOfCell', &
                   gmStreamFuncTopOfCell)
         call mpas_pool_get_array(diagnosticsPool, 'gmStreamFuncTopOfEdge', &
                   gmStreamFuncTopOfEdge)
         call mpas_pool_get_array(diagnosticsPool, 'GMStreamFuncX', &
                   GMStreamFuncX)
         call mpas_pool_get_array(diagnosticsPool, 'GMStreamFuncY', &
                   GMStreamFuncY)
         call mpas_pool_get_array(diagnosticsPool, 'GMStreamFuncZ', &
                   GMStreamFuncZ)
         call mpas_pool_get_array(diagnosticsPool, 'GMStreamFuncZonal', &
                   GMStreamFuncZonal)
         call mpas_pool_get_array(diagnosticsPool, 'GMStreamFuncMeridional', &
                   GMStreamFuncMeridional)
         call mpas_pool_get_array(diagnosticsPool, 'GMBolusVelocityX', &
                   GMBolusVelocityX)
         call mpas_pool_get_array(diagnosticsPool, 'GMBolusVelocityY', &
                   GMBolusVelocityY)
         call mpas_pool_get_array(diagnosticsPool, 'GMBolusVelocityZ', &
                   GMBolusVelocityZ)
         call mpas_pool_get_array(diagnosticsPool, 'GMBolusVelocityZonal', &
                   GMBolusVelocityZonal)
         call mpas_pool_get_array(diagnosticsPool, 'GMBolusVelocityMeridional', &
                   GMBolusVelocityMeridional)
         call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', &
                   normalGMBolusVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'vertGMBolusVelocityTop', &
                   vertGMBolusVelocityTop)

         call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdgeFlux', &
                   layerThickEdgeFlux)
         call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdgeMean', &
                   layerThickEdgeMean)

         call mpas_pool_get_array(diagnosticsPool, &
                  'normalizedRelativeVorticityEdge', &
                   normRelVortEdge)
         call mpas_pool_get_array(diagnosticsPool, &
                  'normalizedPlanetaryVorticityEdge', &
                   normPlanetVortEdge)
         call mpas_pool_get_array(diagnosticsPool, 'relativeVorticity', &
                   relativeVorticity)
         call mpas_pool_get_array(diagnosticsPool, 'relativeVorticityCell', &
                   relativeVorticityCell)

         call mpas_pool_get_array(diagnosticsPool, 'surfaceVelocity', &
                   surfaceVelocity)

         call mpas_pool_get_array(diagnosticsPool, 'vertVelocityTop', &
                   vertVelocityTop)
         call mpas_pool_get_array(diagnosticsPool, 'vertNonLocalFlux', &
                   vertNonLocalFlux)
         call mpas_pool_get_array(diagnosticsPool, 'vertViscTopOfEdge', &
                   vertViscTopOfEdge)
         call mpas_pool_get_array(diagnosticsPool, 'vertViscTopOfCell', &
                   vertViscTopOfCell)
         call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', &
                   vertDiffTopOfCell)
         call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', &
                   vertAleTransportTop)
         call mpas_pool_get_array(diagnosticsPool, 'vertTransportVelocityTop', &
                   vertTransportVelocityTop)

         call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', &
                   barotropicForcing)
         call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', &
                   kineticEnergyCell)

         call mpas_pool_get_array(diagnosticsPool, 'transportVelocityX', &
                   transportVelocityX)
         call mpas_pool_get_array(diagnosticsPool, 'transportVelocityY', &
                   transportVelocityY)
         call mpas_pool_get_array(diagnosticsPool, 'transportVelocityZ', &
                   transportVelocityZ)
         call mpas_pool_get_array(diagnosticsPool, 'transportVelocityZonal', &
                   transportVelocityZonal)
         call mpas_pool_get_array(diagnosticsPool, &
                  'transportVelocityMeridional', &
                   transportVelocityMeridional)

         call mpas_pool_get_array(diagnosticsPool, 'tangentialVelocity', &
                   tangentialVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'montgomeryPotential', &
                   montgomeryPotential)
         call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop', &
                   BruntVaisalaFreqTop)

         call mpas_pool_get_array(diagnosticsPool, 'bulkRichardsonNumber', &
                   bulkRichardsonNumber)
         call mpas_pool_get_array(diagnosticsPool, 'bulkRichardsonNumberBuoy',&
                   bulkRichardsonNumberBuoy)
         call mpas_pool_get_array(diagnosticsPool, 'bulkRichardsonNumberShear',&
                   bulkRichardsonNumberShear)

         ! set pointers for fields related forcing at ocean surface
         call mpas_pool_get_array(diagnosticsPool, 'surfaceBuoyancyForcing', &
                   surfaceBuoyancyForcing)
         call mpas_pool_get_array(diagnosticsPool, 'surfaceFrictionVelocity', &
                   surfaceFrictionVelocity)

         call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', &
                   barotropicThicknessFlux)

         call mpas_pool_get_array(diagnosticsPool, &
                  'normalTransportVelocity', &
                   normalTransportVelocity)
         call mpas_pool_get_array(diagnosticsPool, &
                  'normalVelocitySurfaceLayer', &
                   normalVelocitySurfaceLayer)
         call mpas_pool_get_array(diagnosticsPool, &
                  'penetrativeTemperatureFluxOBL', &
                   penetrativeTemperatureFluxOBL)
         call mpas_pool_get_array(diagnosticsPool, &
                  'normalizedRelativeVorticityCell', &
                   normalizedRelativeVorticityCell)

         call mpas_pool_get_array(diagnosticsPool, 'tracersSurfaceLayerValue', &
                   tracersSurfaceLayerValue)

         ! Set pointer for index variables
         call mpas_pool_get_dimension(diagnosticsPool, &
                  'index_surfaceVelocityZonal', &
                   indexSurfaceVelocityZonal)
         call mpas_pool_get_dimension(diagnosticsPool, &
                  'index_surfaceVelocityMeridional', &
                   indexSurfaceVelocityMeridional)
         call mpas_pool_get_dimension(diagnosticsPool, &
                  'index_vertNonLocalFluxTemp', &
                   index_vertNonLocalFluxTemp)
         call mpas_pool_get_dimension(diagnosticsPool, &
                  'index_SSHGradientMeridional', &
                   indexSSHGradientMeridional)
         call mpas_pool_get_array(diagnosticsPool, &
                  'indexSurfaceLayerDepth', &
                   indexSurfaceLayerDepth)
         call mpas_pool_get_array(diagnosticsPool, &
                  'indexBoundaryLayerDepth', &
                   indexBoundaryLayerDepth)
         call mpas_pool_get_dimension(diagnosticsPool, &
                  'index_SSHGradientZonal', &
                   indexSSHGradientZonal)
      end if

      call mpas_pool_get_array(diagnosticsPool, &
               "daysSinceStartOfSim", &
                daysSinceStartOfSim)
      call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)
      call mpas_pool_get_array(diagnosticsPool, 'simulationStartTime', &
                simulationStartTime)

      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
      call mpas_pool_get_array(diagnosticsPool, 'zTop', zTop)
      call mpas_pool_get_array(diagnosticsPool, 'density', density)
      call mpas_pool_get_array(diagnosticsPool, 'pressure', pressure)
      call mpas_pool_get_array(diagnosticsPool, 'displacedDensity', &
                displacedDensity)
      call mpas_pool_get_array(diagnosticsPool, 'potentialDensity', &
                potentialDensity)
      call mpas_pool_get_array(diagnosticsPool, &
               'inSituThermalExpansionCoeff', &
                thermExpCoeff)
      call mpas_pool_get_array(diagnosticsPool, &
               'inSituSalineContractionCoeff', &
                salineContractCoeff)
      call mpas_pool_get_array(diagnosticsPool, &
               'pressureAdjustedSSH', &
                pressureAdjustedSSH)
      call mpas_pool_get_array(diagnosticsPool, &
               'tracersSurfaceValue', &
                tracersSurfaceValue)
      call mpas_pool_get_array(diagnosticsPool, &
               'surfaceFluxAttenuationCoefficient', &
                sfcFlxAttCoeff)
      call mpas_pool_get_array(diagnosticsPool, &
               'surfaceFluxAttenuationCoefficientRunoff', &
                surfaceFluxAttenuationCoefficientRunoff)
      call mpas_pool_get_array(diagnosticsPool, &
               'boundaryLayerDepth', &
                boundaryLayerDepth)
      call mpas_pool_get_array(diagnosticsPool, &
               'boundaryLayerDepthSmooth', &
                boundaryLayerDepthSmooth)
      call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm', &
                barotropicCoriolisTerm)
      call mpas_pool_get_array(diagnosticsPool, 'betaEdge', betaEdge)
      call mpas_pool_get_array(diagnosticsPool, 'eddyLength', eddyLength)
      call mpas_pool_get_array(diagnosticsPool, 'eddyTime', eddyTime)
      call mpas_pool_get_array(diagnosticsPool, &
               'bottomLayerShortwaveTemperatureFlux', &
                bottomLayerShortwaveTemperatureFlux)
      call mpas_pool_get_field(diagnosticsPool, 'rx1InitSmoothingMask', &
                smoothingMaskField, 1)
      call mpas_pool_get_field(diagnosticsPool, 'verticalStretch', &
                verticalStretchField, 1)
      call mpas_pool_get_array(diagnosticsPool, &
               'modifyLandIcePressureMask', &
                modifyLandIcePressureMask)
      call mpas_pool_get_array(diagnosticsPool, 'indMLD', indMLD)
      call mpas_pool_get_array(diagnosticsPool, 'surfacePressure', &
                surfacePressure)

      ! Semi-implicit Array Pointer retrievals
      if (trim(config_time_integrator) == 'split_implicit') then
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_r0', SIvec_r0)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_r1', SIvec_r1)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_v0', SIvec_v0)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_w0', SIvec_w0)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_w1', SIvec_w1)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_t0', SIvec_t0)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_q0', SIvec_q0)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_s0', SIvec_s0)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_s1', SIvec_s1)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_t0', SIvec_t0)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_t1', SIvec_t1)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_y0', SIvec_y0)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_z0', SIvec_z0)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_z1', SIvec_z1)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_r00', SIvec_r00)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_rh0', SIvec_rh0)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_rh1', SIvec_rh1)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_wh0', SIvec_wh0)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_wh1', SIvec_wh1)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_ph0', SIvec_ph0)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_ph1', SIvec_ph1)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_qh0', SIvec_qh0)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_sh0', SIvec_sh0)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_sh1', SIvec_sh1)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_zh0', SIvec_zh0)
         call mpas_pool_get_array(diagnosticsPool, 'SIvec_zh1', SIvec_zh1)
      end if

      ! Copy diagnostic variables to accelerator
      if ( trim(config_land_ice_flux_mode) /= 'off' ) then
         !$acc enter data create(                                 &
         !$acc                   topDrag,                         &
         !$acc                   topDragMag,                      &
         !$acc                   indexBLT,                        &
         !$acc                   indexBLS,                        &
         !$acc                   indexHeatTrans,                  &
         !$acc                   indexSaltTrans,                  &
         !$acc                   landIceBoundaryLayerTracers,     &
         !$acc                   landIceTracerTransferVelocities, &
         !$acc                   landIceFrictionVelocity)
      end if
      if (config_use_topographic_wave_drag) then
         !$acc enter data create(topographic_wave_drag)
      end if
      if ( trim(config_ocean_run_mode) == 'init' ) then
         !$acc enter data create(                          &
         !$acc                   rx1Edge,                  &
         !$acc                   rx1Cell,                  &
         !$acc                   rx1MaxEdge,               &
         !$acc                   smoothingMask,            &
         !$acc                   verticalStretch,          &
         !$acc                   globalVerticalStretchMax, &
         !$acc                   globalVerticalStretchMin  &
         !$acc                   )
      end if
      if ( trim(config_ocean_run_mode) == 'init' .or. &
                config_AM_debugDiagnostics_enable ) then
         !$acc enter data create(                          &
         !$acc                   rx1MaxCell,               &
         !$acc                   globalRx1Max)
      end if
      if ( trim(config_ocean_run_mode) == 'forward' .or. &
           trim(config_ocean_run_mode) == 'analysis' ) then
         !$acc enter data create(                                 &
         !$acc                   k33,                             &
         !$acc                   gradSSH,                         &
         !$acc                   gradSSHX,                        &
         !$acc                   gradSSHY,                        &
         !$acc                   gradSSHZ,                        &
         !$acc                   viscosity,                       &
         !$acc                   velocityX,                       &
         !$acc                   velocityY,                       &
         !$acc                   velocityZ,                       &
         !$acc                   divergence,                      &
         !$acc                   RiTopOfCell,                     &
         !$acc                   RiTopOfEdge,                     &
         !$acc                   SSHGradient,                     &
         !$acc                   circulation,                     &
         !$acc                   slopeTriadUp,                    &
         !$acc                   gradSSHZonal,                    &
         !$acc                   cGMphaseSpeed,                   &
         !$acc                   velocityZonal,                   &
         !$acc                   GMStreamFuncX,                   &
         !$acc                   GMStreamFuncY,                   &
         !$acc                   GMStreamFuncZ,                   &
         !$acc                   layerThickEdgeFlux,              &
         !$acc                   layerThickEdgeMean,              &
         !$acc                   slopeTriadDown,                  &
         !$acc                   normRelVortEdge,                 &
         !$acc                   surfaceVelocity,                 &
         !$acc                   vertVelocityTop,                 &
         !$acc                   GMBolusVelocityX,                &
         !$acc                   GMBolusVelocityY,                &
         !$acc                   GMBolusVelocityZ,                &
         !$acc                   vertNonLocalFlux,                &
         !$acc                   vertViscTopOfEdge,               &
         !$acc                   gradSSHMeridional,               &
         !$acc                   barotropicForcing,               &
         !$acc                   vertViscTopOfCell,               &
         !$acc                   vertDiffTopOfCell,               &
         !$acc                   GMStreamFuncZonal,               &
         !$acc                   relativeVorticity,               &
         !$acc                   kineticEnergyCell,               &
         !$acc                   normPlanetVortEdge,              &
         !$acc                   velocityMeridional,              &
         !$acc                   transportVelocityX,              &
         !$acc                   transportVelocityY,              &
         !$acc                   transportVelocityZ,              &
         !$acc                   tangentialVelocity,              &
         !$acc                   montgomeryPotential,             &
         !$acc                   BruntVaisalaFreqTop,             &
         !$acc                   vertAleTransportTop,             &
         !$acc                   GMBolusVelocityZonal,            &
         !$acc                   bulkRichardsonNumber,            &
         !$acc                   gmStreamFuncTopOfCell,           &
         !$acc                   gmStreamFuncTopOfEdge,           &
         !$acc                   indexSSHGradientZonal,           &
         !$acc                   relativeVorticityCell,           &
         !$acc                   normalGMBolusVelocity,           &
         !$acc                   vertGMBolusVelocityTop,          &
         !$acc                   transportVelocityZonal,          &
         !$acc                   GMStreamFuncMeridional,          &
         !$acc                   indexSurfaceLayerDepth,          &
         !$acc                   surfaceBuoyancyForcing,          &
         !$acc                   surfaceFrictionVelocity,         &
         !$acc                   indexBoundaryLayerDepth,         &
         !$acc                   barotropicThicknessFlux,         &
         !$acc                   normalTransportVelocity,         &
         !$acc                   vertTransportVelocityTop,        &
         !$acc                   bulkRichardsonNumberBuoy,        &
         !$acc                   tracersSurfaceLayerValue,        &
         !$acc                   GMBolusVelocityMeridional,       &
         !$acc                   bulkRichardsonNumberShear,       &
         !$acc                   indexSurfaceVelocityZonal,       &
         !$acc                   normalVelocitySurfaceLayer,      &
         !$acc                   index_vertNonLocalFluxTemp,      &
         !$acc                   indexSSHGradientMeridional,      &
         !$acc                   transportVelocityMeridional,     &
         !$acc                   penetrativeTemperatureFluxOBL,   &
         !$acc                   indexSurfaceVelocityMeridional,  &
         !$acc                   normalizedRelativeVorticityCell  &
         !$acc                   )
      end if
      if ( trim(config_ocean_run_mode) == 'forward' ) then
         !$acc enter data create( &
         !$acc                   wettingVelocity              &
         !$acc                   )
      end if
      if ( config_use_GM ) then
         !$acc enter data create(                   &
         !$acc                   RediKappa,         &
         !$acc                   gmBolusKappa,      &
         !$acc                   gmKappaScaling,    &
         !$acc                   RediKappaScaling,  &
         !$acc                   rediLimiterCount,  &
         !$acc                   RediKappaSfcTaper, &
         !$acc                   gmResolutionTaper  &
         !$acc                   )
      end if
      if ( config_compute_active_tracer_budgets ) then
         !$acc enter data create(                                         &
         !$acc                   activeTracerHorMixTendency,              &
         !$acc                   activeTracerVertMixTendency,             &
         !$acc                   temperatureShortWaveTendency,            &
         !$acc                   activeTracerNonLocalTendency,            &
         !$acc                   activeTracerSurfaceFluxTendency,         &
         !$acc                   activeTracerVerticalAdvectionTopFlux,    &
         !$acc                   activeTracerVerticalAdvectionTendency,   &
         !$acc                   activeTracerHorizontalAdvectionTendency, &
         !$acc                   activeTracerHorizontalAdvectionEdgeFlux  &
         !$acc                   )
      end if
      if ( config_use_activeTracers_surface_restoring ) then
         !$acc enter data create(                                  &
         !$acc                   salinitySurfaceRestoringTendency  &
         !$acc                   )
      end if
      if ( trim(config_time_integrator) == 'split_implicit' ) then
         !$acc enter data create(                       &
         !$acc                   SIvec_r0,              &
         !$acc                   SIvec_r00,             &
         !$acc                   SIvec_r1,              &
         !$acc                   SIvec_rh0,             &
         !$acc                   SIvec_rh1,             &
         !$acc                   SIvec_ph0,             &
         !$acc                   SIvec_ph1,             &
         !$acc                   SIvec_v0,              &
         !$acc                   SIvec_v0,              &
         !$acc                   SIvec_s0,              &
         !$acc                   SIvec_s1,              &
         !$acc                   SIvec_sh0,             &
         !$acc                   SIvec_sh1,             &
         !$acc                   SIvec_t0,              &
         !$acc                   SIvec_t1,              &
         !$acc                   SIvec_q0,              &
         !$acc                   SIvec_q1,              &
         !$acc                   SIvec_qh0,             &
         !$acc                   SIvec_w0,              &
         !$acc                   SIvec_w1,              &
         !$acc                   SIvec_wh0,             &
         !$acc                   SIvec_wh1,             &
         !$acc                   SIvec_y0,              &
         !$acc                   SIvec_z0,              &
         !$acc                   SIvec_z1,              &
         !$acc                   SIvec_zh0,             &
         !$acc                   SIvec_zh1,             &
         !$acc                   barotropicCoriolisTerm &
         !$acc                   )
      end if

      !$acc enter data create(                                         &
      !$acc                   zMid,                                    &
      !$acc                   zTop,                                    &
      !$acc                   xtime,                                   &
      !$acc                   indMLD,                                  &
      !$acc                   density,                                 &
      !$acc                   betaEdge,                                &
      !$acc                   eddyTime,                                &
      !$acc                   pressure,                                &
      !$acc                   eddyLength,                              &
      !$acc                   thermExpCoeff,                           &
      !$acc                   sfcFlxAttCoeff,                          &
      !$acc                   surfacePressure,                         &
      !$acc                   potentialDensity,                        &
      !$acc                   displacedDensity,                        &
      !$acc                   boundaryLayerDepth,                      &
      !$acc                   salineContractCoeff,                     &
      !$acc                   pressureAdjustedSSH,                     &
      !$acc                   tracersSurfaceValue,                     &
      !$acc                   daysSinceStartOfSim,                     &
      !$acc                   simulationStartTime,                     &
      !$acc                   boundaryLayerDepthSmooth,                &
      !$acc                   surfaceFluxAttenuationCoefficientRunoff  &
      !$acc                   )

    end subroutine ocn_diagnostics_variables_init!}}}

!*******************************************************************************
!
!  ocn_diagnostics_variables_destroy
!
!> \brief Remove diagnostics variables from device
!> \author Rob Aulwes and Phil Jones
!> \date   14 Jan 2020
!> \details
!> This module removes the diagnostics variables from the device
!
!-------------------------------------------------------------------------------

   subroutine ocn_diagnostics_variables_destroy(err) !{{{

      ! Input variables

      ! Output variables

      integer, intent(out) :: &
         err                   ! returned error flag

      ! Local variables

      !***
      !*** end of preamble, begin code
      !***

      err = 0

      ! Remove data from the device
      if ( trim(config_land_ice_flux_mode) /= 'off' ) then
         !$acc exit data delete(                                  &
         !$acc                   topDrag,                         &
         !$acc                   indexBLT,                        &
         !$acc                   indexBLS,                        &
         !$acc                   topDragMag,                      &
         !$acc                   indexHeatTrans,                  &
         !$acc                   indexSaltTrans,                  &
         !$acc                   landIceBoundaryLayerTracers,     &
         !$acc                   landIceTracerTransferVelocities, &
         !$acc                   landIceFrictionVelocity)
      end if
      if (config_use_topographic_wave_drag) then
         !$acc exit data delete(topographic_wave_drag)
      end if
      if ( trim(config_ocean_run_mode) == 'init' ) then
         !$acc exit data delete(                           &
         !$acc                   rx1Edge,                  &
         !$acc                   rx1Cell,                  &
         !$acc                   rx1MaxEdge,               &
         !$acc                   smoothingMask,            &
         !$acc                   verticalStretch,          &
         !$acc                   globalVerticalStretchMax, &
         !$acc                   globalVerticalStretchMin  &
         !$acc                   )
      end if
      if ( trim(config_ocean_run_mode) == 'init' .or. &
                config_AM_debugDiagnostics_enable ) then
         !$acc exit data delete(                           &
         !$acc                   rx1MaxCell,               &
         !$acc                   globalRx1Max)
      end if
      if ( trim(config_ocean_run_mode) == 'forward' .or. &
           trim(config_ocean_run_mode) == 'analysis' ) then
         !$acc exit data delete(                                  &
         !$acc                   k33,                             &
         !$acc                   gradSSH,                         &
         !$acc                   gradSSHX,                        &
         !$acc                   gradSSHY,                        &
         !$acc                   gradSSHZ,                        &
         !$acc                   viscosity,                       &
         !$acc                   velocityX,                       &
         !$acc                   velocityY,                       &
         !$acc                   velocityZ,                       &
         !$acc                   divergence,                      &
         !$acc                   RiTopOfCell,                     &
         !$acc                   RiTopOfEdge,                     &
         !$acc                   SSHGradient,                     &
         !$acc                   circulation,                     &
         !$acc                   slopeTriadUp,                    &
         !$acc                   gradSSHZonal,                    &
         !$acc                   cGMphaseSpeed,                   &
         !$acc                   velocityZonal,                   &
         !$acc                   GMStreamFuncX,                   &
         !$acc                   GMStreamFuncY,                   &
         !$acc                   GMStreamFuncZ,                   &
         !$acc                   layerThickEdgeFlux,              &
         !$acc                   layerThickEdgeMean,              &
         !$acc                   slopeTriadDown,                  &
         !$acc                   normRelVortEdge,                 &
         !$acc                   surfaceVelocity,                 &
         !$acc                   vertVelocityTop,                 &
         !$acc                   GMBolusVelocityX,                &
         !$acc                   GMBolusVelocityY,                &
         !$acc                   GMBolusVelocityZ,                &
         !$acc                   vertNonLocalFlux,                &
         !$acc                   vertViscTopOfEdge,               &
         !$acc                   gradSSHMeridional,               &
         !$acc                   barotropicForcing,               &
         !$acc                   vertViscTopOfCell,               &
         !$acc                   vertDiffTopOfCell,               &
         !$acc                   GMStreamFuncZonal,               &
         !$acc                   relativeVorticity,               &
         !$acc                   kineticEnergyCell,               &
         !$acc                   normPlanetVortEdge,              &
         !$acc                   velocityMeridional,              &
         !$acc                   transportVelocityX,              &
         !$acc                   transportVelocityY,              &
         !$acc                   transportVelocityZ,              &
         !$acc                   tangentialVelocity,              &
         !$acc                   montgomeryPotential,             &
         !$acc                   BruntVaisalaFreqTop,             &
         !$acc                   vertAleTransportTop,             &
         !$acc                   GMBolusVelocityZonal,            &
         !$acc                   bulkRichardsonNumber,            &
         !$acc                   gmStreamFuncTopOfCell,           &
         !$acc                   gmStreamFuncTopOfEdge,           &
         !$acc                   indexSSHGradientZonal,           &
         !$acc                   relativeVorticityCell,           &
         !$acc                   normalGMBolusVelocity,           &
         !$acc                   vertGMBolusVelocityTop,          &
         !$acc                   transportVelocityZonal,          &
         !$acc                   GMStreamFuncMeridional,          &
         !$acc                   indexSurfaceLayerDepth,          &
         !$acc                   surfaceBuoyancyForcing,          &
         !$acc                   surfaceFrictionVelocity,         &
         !$acc                   indexBoundaryLayerDepth,         &
         !$acc                   barotropicThicknessFlux,         &
         !$acc                   normalTransportVelocity,         &
         !$acc                   vertTransportVelocityTop,        &
         !$acc                   bulkRichardsonNumberBuoy,        &
         !$acc                   tracersSurfaceLayerValue,        &
         !$acc                   GMBolusVelocityMeridional,       &
         !$acc                   bulkRichardsonNumberShear,       &
         !$acc                   indexSurfaceVelocityZonal,       &
         !$acc                   normalVelocitySurfaceLayer,      &
         !$acc                   index_vertNonLocalFluxTemp,      &
         !$acc                   indexSSHGradientMeridional,      &
         !$acc                   transportVelocityMeridional,     &
         !$acc                   penetrativeTemperatureFluxOBL,   &
         !$acc                   indexSurfaceVelocityMeridional,  &
         !$acc                   normalizedRelativeVorticityCell  &
         !$acc                   )
      end if
      if ( trim(config_ocean_run_mode) == 'forward' ) then
         !$acc exit data delete(                              &
         !$acc                   wettingVelocity              &
         !$acc                   )
      end if
      if ( config_use_GM ) then
         !$acc exit data delete(                    &
         !$acc                   RediKappa,         &
         !$acc                   gmBolusKappa,      &
         !$acc                   gmKappaScaling,    &
         !$acc                   RediKappaScaling,  &
         !$acc                   rediLimiterCount,  &
         !$acc                   RediKappaSfcTaper, &
         !$acc                   gmResolutionTaper  &
         !$acc                   )
      end if
      if ( config_compute_active_tracer_budgets ) then
         !$acc exit data delete(                                          &
         !$acc                   activeTracerHorMixTendency,              &
         !$acc                   activeTracerVertMixTendency,             &
         !$acc                   temperatureShortWaveTendency,            &
         !$acc                   activeTracerNonLocalTendency,            &
         !$acc                   activeTracerSurfaceFluxTendency,         &
         !$acc                   activeTracerVerticalAdvectionTopFlux,    &
         !$acc                   activeTracerVerticalAdvectionTendency,   &
         !$acc                   activeTracerHorizontalAdvectionTendency, &
         !$acc                   activeTracerHorizontalAdvectionEdgeFlux  &
         !$acc                   )
      end if
      if ( config_use_activeTracers_surface_restoring ) then
         !$acc exit data delete(                                   &
         !$acc                   salinitySurfaceRestoringTendency  &
         !$acc                   )
      end if
      if ( trim(config_time_integrator) == 'split_implicit' ) then
         !$acc exit data delete(                        &
         !$acc                   SIvec_r0,              &
         !$acc                   SIvec_r00,             &
         !$acc                   SIvec_r1,              &
         !$acc                   SIvec_rh0,             &
         !$acc                   SIvec_rh1,             &
         !$acc                   SIvec_ph0,             &
         !$acc                   SIvec_ph1,             &
         !$acc                   SIvec_v0,              &
         !$acc                   SIvec_v0,              &
         !$acc                   SIvec_s0,              &
         !$acc                   SIvec_s1,              &
         !$acc                   SIvec_sh0,             &
         !$acc                   SIvec_sh1,             &
         !$acc                   SIvec_t0,              &
         !$acc                   SIvec_t1,              &
         !$acc                   SIvec_q0,              &
         !$acc                   SIvec_q1,              &
         !$acc                   SIvec_qh0,             &
         !$acc                   SIvec_w0,              &
         !$acc                   SIvec_w1,              &
         !$acc                   SIvec_wh0,             &
         !$acc                   SIvec_wh1,             &
         !$acc                   SIvec_y0,              &
         !$acc                   SIvec_z0,              &
         !$acc                   SIvec_z1,              &
         !$acc                   SIvec_zh0,             &
         !$acc                   SIvec_zh1,             &
         !$acc                   barotropicCoriolisTerm &
         !$acc                   )
      end if

      !$acc exit data delete(                                          &
      !$acc                   zMid,                                    &
      !$acc                   zTop,                                    &
      !$acc                   xtime,                                   &
      !$acc                   indMLD,                                  &
      !$acc                   density,                                 &
      !$acc                   betaEdge,                                &
      !$acc                   eddyTime,                                &
      !$acc                   pressure,                                &
      !$acc                   eddyLength,                              &
      !$acc                   thermExpCoeff,                           &
      !$acc                   sfcFlxAttCoeff,                          &
      !$acc                   surfacePressure,                         &
      !$acc                   potentialDensity,                        &
      !$acc                   displacedDensity,                        &
      !$acc                   boundaryLayerDepth,                      &
      !$acc                   salineContractCoeff,                     &
      !$acc                   pressureAdjustedSSH,                     &
      !$acc                   tracersSurfaceValue,                     &
      !$acc                   daysSinceStartOfSim,                     &
      !$acc                   simulationStartTime,                     &
      !$acc                   boundaryLayerDepthSmooth,                &
      !$acc                   surfaceFluxAttenuationCoefficientRunoff  &
      !$acc                   )

      ! Nullify pointers
      if ( trim(config_land_ice_flux_mode) /= 'off' ) then
         nullify(topDrag,                         &
                 topDragMag,                      &
                 indexBLT,                        &
                 indexBLS,                        &
                 indexHeatTrans,                  &
                 indexSaltTrans,                  &
                 landIceBoundaryLayerTracers,     &
                 landIceTracerTransferVelocities, &
                 landIceFrictionVelocity)
      end if
      if ( trim(config_ocean_run_mode) == 'init' ) then
         nullify(rx1Edge,                  &
                 rx1Cell,                  &
                 rx1MaxEdge,               &
                 smoothingMask,            &
                 verticalStretch,          &
                 globalVerticalStretchMax, &
                 globalVerticalStretchMin)
      end if
      if ( trim(config_ocean_run_mode) == 'init' .or. &
                config_AM_debugDiagnostics_enable ) then
         nullify(rx1MaxCell,               &
                 globalRx1Max)
      end if
      if ( trim(config_ocean_run_mode) == 'forward' .or. &
           trim(config_ocean_run_mode) == 'analysis' ) then
         nullify(k33,                             &
                 gradSSH,                         &
                 gradSSHX,                        &
                 gradSSHY,                        &
                 gradSSHZ,                        &
                 viscosity,                       &
                 velocityX,                       &
                 velocityY,                       &
                 velocityZ,                       &
                 divergence,                      &
                 RiTopOfCell,                     &
                 RiTopOfEdge,                     &
                 SSHGradient,                     &
                 circulation,                     &
                 slopeTriadUp,                    &
                 gradSSHZonal,                    &
                 cGMphaseSpeed,                   &
                 velocityZonal,                   &
                 GMStreamFuncX,                   &
                 GMStreamFuncY,                   &
                 GMStreamFuncZ,                   &
                 layerThickEdgeFlux,              &
                 layerThickEdgeMean,              &
                 slopeTriadDown,                  &
                 normRelVortEdge,                 &
                 surfaceVelocity,                 &
                 vertVelocityTop,                 &
                 GMBolusVelocityX,                &
                 GMBolusVelocityY,                &
                 GMBolusVelocityZ,                &
                 vertNonLocalFlux,                &
                 vertViscTopOfEdge,               &
                 gradSSHMeridional,               &
                 barotropicForcing,               &
                 vertViscTopOfCell,               &
                 vertDiffTopOfCell,               &
                 GMStreamFuncZonal,               &
                 relativeVorticity,               &
                 kineticEnergyCell,               &
                 normPlanetVortEdge,              &
                 velocityMeridional,              &
                 transportVelocityX,              &
                 transportVelocityY,              &
                 transportVelocityZ,              &
                 tangentialVelocity,              &
                 montgomeryPotential,             &
                 BruntVaisalaFreqTop,             &
                 vertAleTransportTop,             &
                 GMBolusVelocityZonal,            &
                 bulkRichardsonNumber,            &
                 gmStreamFuncTopOfCell,           &
                 gmStreamFuncTopOfEdge,           &
                 indexSSHGradientZonal,           &
                 relativeVorticityCell,           &
                 normalGMBolusVelocity,           &
                 vertGMBolusVelocityTop,          &
                 transportVelocityZonal,          &
                 GMStreamFuncMeridional,          &
                 indexSurfaceLayerDepth,          &
                 surfaceBuoyancyForcing,          &
                 surfaceFrictionVelocity,         &
                 indexBoundaryLayerDepth,         &
                 barotropicThicknessFlux,         &
                 normalTransportVelocity,         &
                 vertTransportVelocityTop,        &
                 bulkRichardsonNumberBuoy,        &
                 tracersSurfaceLayerValue,        &
                 GMBolusVelocityMeridional,       &
                 bulkRichardsonNumberShear,       &
                 indexSurfaceVelocityZonal,       &
                 normalVelocitySurfaceLayer,      &
                 index_vertNonLocalFluxTemp,      &
                 indexSSHGradientMeridional,      &
                 transportVelocityMeridional,     &
                 penetrativeTemperatureFluxOBL,   &
                 indexSurfaceVelocityMeridional,  &
                 normalizedRelativeVorticityCell)
      end if
      if ( trim(config_ocean_run_mode) == 'forward' ) then
         nullify(wettingVelocity)
      end if
      if ( config_use_GM ) then
         nullify(RediKappa,         &
                 gmBolusKappa,      &
                 gmKappaScaling,    &
                 RediKappaScaling,  &
                 rediLimiterCount,  &
                 RediKappaSfcTaper, &
                 gmResolutionTaper)
      end if
      if ( config_compute_active_tracer_budgets ) then
         nullify(activeTracerHorMixTendency,              &
                 activeTracerVertMixTendency,             &
                 temperatureShortWaveTendency,            &
                 activeTracerNonLocalTendency,            &
                 activeTracerSurfaceFluxTendency,         &
                 activeTracerVerticalAdvectionTopFlux,    &
                 activeTracerVerticalAdvectionTendency,   &
                 activeTracerHorizontalAdvectionTendency, &
                 activeTracerHorizontalAdvectionEdgeFlux)
      end if
      if ( config_use_activeTracers_surface_restoring ) then
         nullify(salinitySurfaceRestoringTendency)
      end if
      if ( trim(config_time_integrator) == 'split_implicit' ) then
         nullify(SIvec_r0,              &
                 SIvec_r00,             &
                 SIvec_r1,              &
                 SIvec_rh0,             &
                 SIvec_rh1,             &
                 SIvec_ph0,             &
                 SIvec_ph1,             &
                 SIvec_v0,              &
                 SIvec_v0,              &
                 SIvec_s0,              &
                 SIvec_s1,              &
                 SIvec_sh0,             &
                 SIvec_sh1,             &
                 SIvec_t0,              &
                 SIvec_t1,              &
                 SIvec_q0,              &
                 SIvec_q1,              &
                 SIvec_qh0,             &
                 SIvec_w0,              &
                 SIvec_w1,              &
                 SIvec_wh0,             &
                 SIvec_wh1,             &
                 SIvec_y0,              &
                 SIvec_z0,              &
                 SIvec_z1,              &
                 SIvec_zh0,             &
                 SIvec_zh1,             &
                 barotropicCoriolisTerm)
      end if

      nullify(zMid,                                    &
              zTop,                                    &
              xtime,                                   &
              indMLD,                                  &
              density,                                 &
              betaEdge,                                &
              eddyTime,                                &
              pressure,                                &
              eddyLength,                              &
              thermExpCoeff,                           &
              sfcFlxAttCoeff,                          &
              surfacePressure,                         &
              potentialDensity,                        &
              displacedDensity,                        &
              boundaryLayerDepth,                      &
              salineContractCoeff,                     &
              pressureAdjustedSSH,                     &
              tracersSurfaceValue,                     &
              daysSinceStartOfSim,                     &
              simulationStartTime,                     &
              boundaryLayerDepthSmooth,                &
              surfaceFluxAttenuationCoefficientRunoff)

    end subroutine ocn_diagnostics_variables_destroy!}}}

end module ocn_diagnostics_variables

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
