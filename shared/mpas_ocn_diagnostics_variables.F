! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_diagnostics
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
   real (kind=RKIND), dimension(:,:), pointer :: layerThickEdge
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

   real(kind=RKIND), dimension(:,:), pointer :: kappaGMCell
   real(kind=RKIND), dimension(:,:), pointer :: gmKappaScaling
   real(kind=RKIND), dimension(:,:), pointer :: RediKappaScaling
   real(kind=RKIND), dimension(:,:), pointer :: RediKappaSfcTaper
   real(kind=RKIND), dimension(:,:), pointer :: RediKappaCell
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

   real (kind=RKIND), dimension(:,:), pointer :: surfaceVelocity, avgSurfaceVelocity
   real (kind=RKIND), dimension(:), pointer :: gradSSHZonal, gradSSHMeridional
   real (kind=RKIND), dimension(:), pointer :: gradSSHX
   real (kind=RKIND), dimension(:), pointer :: gradSSHY
   real (kind=RKIND), dimension(:), pointer :: gradSSHZ

   real (kind=RKIND), dimension(:), pointer :: barotropicForcing
   real (kind=RKIND), dimension(:), pointer :: barotropicCoriolisTerm
   real (kind=RKIND), dimension(:), pointer :: barotropicThicknessFlux

   real (kind=RKIND), dimension(:, :), pointer :: edgeAreaFractionOfCell
   real (kind=RKIND), dimension(:, :), pointer :: SSHGradient
   integer, pointer :: indexSurfaceVelocityZonal, indexSurfaceVelocityMeridional
   integer, pointer :: indexSSHGradientZonal, indexSSHGradientMeridional

   type (field1DInteger), pointer :: smoothingMaskField
   type (field2DReal), pointer :: verticalStretchField

   real (kind=RKIND), pointer :: globalVerticalStretchMax, globalVerticalStretchMin
   real (kind=RKIND), dimension(:), pointer :: betaEdge
   real (kind=RKIND), dimension(:), pointer :: eddyLength
   real (kind=RKIND), dimension(:), pointer :: eddyTime
   real (kind=RKIND), dimension(:), pointer :: c_visbeck
   real (kind=RKIND), dimension(:), pointer :: gmResolutionTaper

   ! Semi-implicit Array Pointers
   real (kind=RKIND), dimension(:), pointer :: CGvec_r0,CGvec_r00,CGvec_r1 ,CGvec_rh0,CGvec_rh1,CGvec_ph0,CGvec_ph1
   real (kind=RKIND), dimension(:), pointer :: CGvec_v0,CGvec_v1 ,CGvec_s0 ,CGvec_s1 ,CGvec_sh0,CGvec_sh1
   real (kind=RKIND), dimension(:), pointer :: CGvec_t0,CGvec_t1 ,CGvec_q0 ,CGvec_q1 ,CGvec_qh0
   real (kind=RKIND), dimension(:), pointer :: CGvec_w0,CGvec_w1 ,CGvec_wh0,CGvec_wh1,CGvec_y0
   real (kind=RKIND), dimension(:), pointer :: CGvec_z0,CGvec_z1 ,CGvec_zh0,CGvec_zh1

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_diagnostics_variables_init

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

      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
      call mpas_pool_get_array(diagnosticsPool, 'zTop', zTop)
      call mpas_pool_get_array(diagnosticsPool, 'divergence', divergence)
      call mpas_pool_get_array(diagnosticsPool, 'circulation', circulation)
      call mpas_pool_get_array(diagnosticsPool, 'relativeVorticity', relativeVorticity)
      call mpas_pool_get_array(diagnosticsPool, 'relativeVorticityCell', relativeVorticityCell)
      call mpas_pool_get_array(diagnosticsPool, 'normalizedPlanetaryVorticityEdge', normPlanetVortEdge)
      call mpas_pool_get_array(diagnosticsPool, 'normalizedRelativeVorticityEdge', normRelVortEdge)
      call mpas_pool_get_array(diagnosticsPool, 'normalizedRelativeVorticityCell', normalizedRelativeVorticityCell)
      call mpas_pool_get_array(diagnosticsPool, 'density', density)
      call mpas_pool_get_array(diagnosticsPool, 'displacedDensity', displacedDensity)
      call mpas_pool_get_array(diagnosticsPool, 'potentialDensity', potentialDensity)
      call mpas_pool_get_array(diagnosticsPool, 'montgomeryPotential', montgomeryPotential)
      call mpas_pool_get_array(diagnosticsPool, 'pressure', pressure)
      call mpas_pool_get_array(diagnosticsPool, 'inSituThermalExpansionCoeff',thermExpCoeff)
      call mpas_pool_get_array(diagnosticsPool, 'inSituSalineContractionCoeff', salineContractCoeff)
      call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop', BruntVaisalaFreqTop)
      call mpas_pool_get_array(diagnosticsPool, 'tangentialVelocity', tangentialVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThickEdge)
      call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCell)
      call mpas_pool_get_array(diagnosticsPool, 'vertVelocityTop', vertVelocityTop)
      call mpas_pool_get_array(diagnosticsPool, 'vertTransportVelocityTop', vertTransportVelocityTop)
      call mpas_pool_get_array(diagnosticsPool, 'vertGMBolusVelocityTop', vertGMBolusVelocityTop)
      call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'gradSSH', gradSSH)
      call mpas_pool_get_array(diagnosticsPool, 'gradSSHX', gradSSHX)
      call mpas_pool_get_array(diagnosticsPool, 'gradSSHY', gradSSHY)
      call mpas_pool_get_array(diagnosticsPool, 'gradSSHZ', gradSSHZ)
      call mpas_pool_get_array(diagnosticsPool, 'RiTopOfCell', RiTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'RiTopOfEdge', RiTopOfEdge)
      call mpas_pool_get_array(diagnosticsPool, 'pressureAdjustedSSH', pressureAdjustedSSH)

      call mpas_pool_get_array(diagnosticsPool, 'tracersSurfaceValue', tracersSurfaceValue)
      call mpas_pool_get_array(diagnosticsPool, 'tracersSurfaceLayerValue', tracersSurfaceLayerValue)
      call mpas_pool_get_array(diagnosticsPool, 'normalVelocitySurfaceLayer', normalVelocitySurfaceLayer)

      call mpas_pool_get_array(diagnosticsPool, 'surfaceFluxAttenuationCoefficient', sfcFlxAttCoeff)
      call mpas_pool_get_array(diagnosticsPool, 'surfaceFluxAttenuationCoefficientRunoff', surfaceFluxAttenuationCoefficientRunoff)

      call mpas_pool_get_array(diagnosticsPool, 'landIceFrictionVelocity', landIceFrictionVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'topDrag', topDrag)
      call mpas_pool_get_array(diagnosticsPool, 'topDragMagnitude', topDragMag)
      call mpas_pool_get_array(diagnosticsPool, 'penetrativeTemperatureFluxOBL', penetrativeTemperatureFluxOBL)
      call mpas_pool_get_array(diagnosticsPool, 'bottomLayerShortwaveTemperatureFlux', bottomLayerShortwaveTemperatureFlux)
      call mpas_pool_get_array(diagnosticsPool, 'indexSurfaceLayerDepth', indexSurfaceLayerDepth)

      !
      ! set pointers for fields related forcing at ocean surface
      !
      call mpas_pool_get_array(diagnosticsPool, 'surfaceBuoyancyForcing', surfaceBuoyancyForcing)
      call mpas_pool_get_array(diagnosticsPool, 'surfaceFrictionVelocity', surfaceFrictionVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'surfacePressure', surfacePressure)

      call mpas_pool_get_array(diagnosticsPool, 'landIceBoundaryLayerTracers', landIceBoundaryLayerTracers)
      call mpas_pool_get_dimension(diagnosticsPool, 'index_landIceBoundaryLayerTemperature', indexBLT)
      call mpas_pool_get_dimension(diagnosticsPool, 'index_landIceBoundaryLayerSalinity', indexBLS)

      call mpas_pool_get_array(diagnosticsPool, 'transportVelocityX', transportVelocityX)
      call mpas_pool_get_array(diagnosticsPool, 'transportVelocityY', transportVelocityY)
      call mpas_pool_get_array(diagnosticsPool, 'transportVelocityZ', transportVelocityZ)
      call mpas_pool_get_array(diagnosticsPool, 'transportVelocityZonal', transportVelocityZonal)
      call mpas_pool_get_array(diagnosticsPool, 'transportVelocityMeridional', transportVelocityMeridional)

      call mpas_pool_get_array(diagnosticsPool, 'GMBolusVelocityX', GMBolusVelocityX)
      call mpas_pool_get_array(diagnosticsPool, 'GMBolusVelocityY', GMBolusVelocityY)
      call mpas_pool_get_array(diagnosticsPool, 'GMBolusVelocityZ', GMBolusVelocityZ)
      call mpas_pool_get_array(diagnosticsPool, 'GMBolusVelocityZonal', GMBolusVelocityZonal)
      call mpas_pool_get_array(diagnosticsPool, 'GMBolusVelocityMeridional', GMBolusVelocityMeridional)

      call mpas_pool_get_array(diagnosticsPool, 'gmStreamFuncTopOfEdge', gmStreamFuncTopOfEdge)
      call mpas_pool_get_array(diagnosticsPool, 'GMStreamFuncX', GMStreamFuncX)
      call mpas_pool_get_array(diagnosticsPool, 'GMStreamFuncY', GMStreamFuncY)
      call mpas_pool_get_array(diagnosticsPool, 'GMStreamFuncZ', GMStreamFuncZ)
      call mpas_pool_get_array(diagnosticsPool, 'GMStreamFuncZonal', GMStreamFuncZonal)
      call mpas_pool_get_array(diagnosticsPool, 'GMStreamFuncMeridional', GMStreamFuncMeridional)

      if(jenkinsOn .or. hollandJenkinsOn) then
        call mpas_pool_get_array(diagnosticsPool, 'landIceTracerTransferVelocities', landIceTracerTransferVelocities)
        call mpas_pool_get_dimension(diagnosticsPool, 'index_landIceHeatTransferVelocity', indexHeatTrans)
        call mpas_pool_get_dimension(diagnosticsPool, 'index_landIceSaltTransferVelocity', indexSaltTrans)
      end if

      call mpas_pool_get_array(diagnosticsPool, 'rx1Edge', rx1Edge, 1)
      call mpas_pool_get_array(diagnosticsPool, 'rx1Cell', rx1Cell, 1)
      call mpas_pool_get_array(diagnosticsPool, 'rx1MaxEdge', rx1MaxEdge, 1)
      call mpas_pool_get_array(diagnosticsPool, 'rx1MaxCell', rx1MaxCell, 1)

      call mpas_pool_get_array(diagnosticsPool, 'globalRx1Max', globalRx1Max, 1)
      call mpas_pool_get_array(diagnosticsPool, 'globalVerticalStretchMax', globalVerticalStretchMax, 1)
      call mpas_pool_get_array(diagnosticsPool, 'globalVerticalStretchMin', globalVerticalStretchMin, 1)

      call mpas_pool_get_array(diagnosticsPool, 'rx1InitSmoothingMask', smoothingMask, 1)
      call mpas_pool_get_array(diagnosticsPool, 'verticalStretch', verticalStretch, 1)
      call mpas_pool_get_field(diagnosticsPool, 'rx1InitSmoothingMask', smoothingMaskField, 1)
      call mpas_pool_get_field(diagnosticsPool, 'verticalStretch', verticalStretchField, 1)

      call mpas_pool_get_array(diagnosticsPool, 'modifyLandIcePressureMask', modifyLandIcePressureMask)

      call mpas_pool_get_array(diagnosticsPool, 'velocityX', velocityX)
      call mpas_pool_get_array(diagnosticsPool, 'velocityY', velocityY)
      call mpas_pool_get_array(diagnosticsPool, 'velocityZ', velocityZ)
      call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velocityMeridional)
      call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velocityZonal)

      call mpas_pool_get_array(diagnosticsPool, 'indMLD', indMLD)

      call mpas_pool_get_array(diagnosticsPool, 'cGMphaseSpeed', cGMphaseSpeed)
      call mpas_pool_get_array(diagnosticsPool, 'kappaGMCell', kappaGMCell)
      call mpas_pool_get_array(diagnosticsPool, 'gmKappaScaling', gmKappaScaling)
      call mpas_pool_get_array(diagnosticsPool, 'RediKappaScaling', RediKappaScaling)
      call mpas_pool_get_array(diagnosticsPool, 'RediKappaSfcTaper', RediKappaSfcTaper)
      call mpas_pool_get_array(diagnosticsPool, 'slopeTriadUp', slopeTriadUp)
      call mpas_pool_get_array(diagnosticsPool, 'slopeTriadDown', slopeTriadDown)
      call mpas_pool_get_array(diagnosticsPool, 'k33', k33)
      call mpas_pool_get_array(diagnosticsPool, 'gmStreamFuncTopOfCell', gmStreamFuncTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'gmBolusKappa', gmBolusKappa)
      call mpas_pool_get_array(diagnosticsPool, 'RediKappaCell', RediKappaCell)

      call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)

      call mpas_pool_get_array(diagnosticsPool, 'viscosity', viscosity)
      call mpas_pool_get_array(diagnosticsPool, 'wettingVelocity', wettingVelocity)
      call mpas_pool_get_field(diagnosticsPool, 'wettingVelocity', wettingVelocityField)

      !
      ! set pointers for viscosity/diffusivity and initialize to zero
      !
      call mpas_pool_get_array(diagnosticsPool, 'vertViscTopOfEdge', vertViscTopOfEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertViscTopOfCell', vertViscTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', vertDiffTopOfCell)

      call mpas_pool_get_array(diagnosticsPool, 'bulkRichardsonNumber', bulkRichardsonNumber)
      call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth', boundaryLayerDepth)
      call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepthSmooth', boundaryLayerDepthSmooth)
      call mpas_pool_get_array(diagnosticsPool, 'bulkRichardsonNumberBuoy',bulkRichardsonNumberBuoy)
      call mpas_pool_get_array(diagnosticsPool, 'bulkRichardsonNumberShear',bulkRichardsonNumberShear)
      call mpas_pool_get_array(diagnosticsPool, 'indexBoundaryLayerDepth',indexBoundaryLayerDepth)
      call mpas_pool_get_array(diagnosticsPool, 'edgeAreaFractionOfCell', edgeAreaFractionOfCell)

      call mpas_pool_get_array(diagnosticsPool, 'vertNonLocalFlux', vertNonLocalFlux)
      call mpas_pool_get_dimension(diagnosticsPool, 'index_vertNonLocalFluxTemp', index_vertNonLocalFluxTemp)

      call mpas_pool_get_array(diagnosticsPool, 'RediKappa', RediKappa)
      call mpas_pool_get_array(diagnosticsPool, 'rediLimiterCount', rediLimiterCount)
      if (config_compute_active_tracer_budgets) then
         call mpas_pool_get_array(diagnosticsPool,'activeTracerSurfaceFluxTendency',activeTracerSurfaceFluxTendency)
         call mpas_pool_get_array(diagnosticsPool,'temperatureShortWaveTendency',temperatureShortWaveTendency)
         call mpas_pool_get_array(diagnosticsPool,'activeTracerNonLocalTendency',activeTracerNonLocalTendency)
         call mpas_pool_get_array(diagnosticsPool,'activeTracerHorMixTendency',activeTracerHorMixTendency)
      endif

      call mpas_pool_get_array(diagnosticsPool,         &
             'activeTracerVerticalAdvectionTopFlux',    &
              activeTracerVerticalAdvectionTopFlux)
      call mpas_pool_get_array(diagnosticsPool,         &
             'activeTracerHorizontalAdvectionEdgeFlux', &
              activeTracerHorizontalAdvectionEdgeFlux)
      call mpas_pool_get_array(diagnosticsPool,         &
             'activeTracerHorizontalAdvectionTendency', &
              activeTracerHorizontalAdvectionTendency)
      call mpas_pool_get_array(diagnosticsPool,         &
             'activeTracerVerticalAdvectionTendency',   &
              activeTracerVerticalAdvectionTendency)
      call mpas_pool_get_array(diagnosticsPool, &
         'activeTracerVertMixTendency',activeTracerVertMixTendency)

      call mpas_pool_get_array(diagnosticsPool, "daysSinceStartOfSim", daysSinceStartOfSim)
      call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)
      call mpas_pool_get_array(diagnosticsPool, 'simulationStartTime',simulationStartTime)

      call mpas_pool_get_array(diagnosticsPool,'salinitySurfaceRestoringTendency',salinitySurfaceRestoringTendency)

      call mpas_pool_get_array(diagnosticsPool, 'surfaceVelocity', surfaceVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'gradSSHZonal', gradSSHZonal)
      call mpas_pool_get_array(diagnosticsPool, 'gradSSHMeridional', gradSSHMeridional)

      call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
      call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm', barotropicCoriolisTerm)
      call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFlux)

      call mpas_pool_get_array(diagnosticsPool, 'SSHGradient', SSHGradient)

      call mpas_pool_get_dimension(diagnosticsPool, 'index_surfaceVelocityZonal', indexSurfaceVelocityZonal)
      call mpas_pool_get_dimension(diagnosticsPool, 'index_surfaceVelocityMeridional', indexSurfaceVelocityMeridional)
      call mpas_pool_get_dimension(diagnosticsPool, 'index_SSHGradientZonal', indexSSHGradientZonal)
      call mpas_pool_get_dimension(diagnosticsPool, 'index_SSHGradientMeridional', indexSSHGradientMeridional)

      call mpas_pool_get_array(diagnosticsPool, 'betaEdge', betaEdge)
      call mpas_pool_get_array(diagnosticsPool, 'eddyLength', eddyLength)
      call mpas_pool_get_array(diagnosticsPool, 'eddyTime', eddyTime)
      call mpas_pool_get_array(diagnosticsPool, 'c_visbeck', c_visbeck)
      call mpas_pool_get_array(diagnosticsPool, 'gmResolutionTaper', gmResolutionTaper)

      ! Semi-implicit Array Pointer retrievals
      if (trim(config_time_integrator) == 'semi_implicit') then
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0', CGvec_r0)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_r1', CGvec_r1)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_v0', CGvec_v0)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_v1', CGvec_v1)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_w0', CGvec_w0)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_w1', CGvec_w1)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_t0', CGvec_t0)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_q0', CGvec_q0)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_s0', CGvec_s0)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_s1', CGvec_s1)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_t0', CGvec_t0)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_t1', CGvec_t1)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_y0', CGvec_y0)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_z0', CGvec_z0)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_z1', CGvec_z1)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_r00', CGvec_r00)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_rh0', CGvec_rh0)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_rh1', CGvec_rh1)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_wh0', CGvec_wh0)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_wh1', CGvec_wh1)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_ph0', CGvec_ph0)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_ph1', CGvec_ph1)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_qh0', CGvec_qh0)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_sh0', CGvec_sh0)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_sh1', CGvec_sh1)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_zh0', CGvec_zh0)
         call mpas_pool_get_array(diagnosticsPool, 'CGvec_zh1', CGvec_zh1)
      end if

    end subroutine ocn_diagnostics_variables_init!}}}

!***********************************************************************

end module ocn_diagnostics_variables

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
