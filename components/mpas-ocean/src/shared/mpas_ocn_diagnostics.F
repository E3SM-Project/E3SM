! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_diagnostics
!
!> \brief MPAS ocean diagnostics driver
!> \author Mark Petersen
!> \date   23 September 2011
!> \details
!>  This module contains the routines for computing
!>  diagnostic variables, and other quantities such as vertAleTransportTop.
!
!-----------------------------------------------------------------------

module ocn_diagnostics

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_threading
   use mpas_vector_reconstruction
   use mpas_stream_manager
   use mpas_io_units

   use ocn_constants
   use ocn_config
   use ocn_gm
   use ocn_equation_of_state
   use ocn_thick_ale
   use ocn_diagnostics_variables
   use ocn_mesh
   use ocn_surface_land_ice_fluxes
   use ocn_vertical_advection
   use ocn_subgrid

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

   public :: ocn_diagnostic_solve, &
             ocn_diagnostic_solve_layerThicknessEdge, &
             ocn_diagnostic_solve_wctEdge, &
             ocn_relativeVorticity_circulation, &
             ocn_vert_transport_velocity_top, &
             ocn_fuperp, &
             ocn_filter_btr_mode_tend_vel, &
             ocn_reconstruct_eddy_vectors, &
             ocn_compute_mixing_input_fields, &
             ocn_validate_state, &
             ocn_build_log_filename, &
             ocn_diagnostics_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: ke_cell_flag, ke_vertex_flag
   real (kind=RKIND) ::  fCoef
   real (kind=RKIND) ::  landIceTopDragCoeff
   real (kind=RKIND), pointer ::  coef_3rd_order

   ! Methods for computing thickness at edges for flux calculations
   integer :: &
      thickEdgeFluxChoice, &! choice of thickness flux type
      thickEdgeDragChoice   ! choice of thickness drag type
   integer, parameter :: &
      thickEdgeFluxCenter   = 1, &! use mean thickness of cell neighbors
      thickEdgeFluxUpwind   = 2, &! use upwind cell thickness at edge
      thickEdgeFluxConstant = 3   ! use constant thickness in time from
                                  ! RestingThickness, for linear test problems
   integer, parameter :: &
      thickEdgeDragCenter  = 1, &! use mean thickness of cell neighbors
      thickEdgeDragHarMean = 2   ! use harmonic mean cell thickness at edge

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_diagnostic_solve
!
!> \brief   Computes diagnostic variables
!> \author  Mark Petersen
!> \date    23 September 2011
!> \details
!>  This routine computes the diagnostic variables for the ocean
!
!-----------------------------------------------------------------------

   subroutine ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, verticalMeshPool, scratchPool, tracersPool, &
                                   timeLevelIn, full)!{{{

      real (kind=RKIND), intent(in) :: dt !< Input: Time step
      type (mpas_pool_type), intent(in) :: statePool !< Input: State information
      type (mpas_pool_type), intent(in) :: forcingPool !< Input: Forcing information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      type (mpas_pool_type), intent(in) :: verticalMeshPool !< Input: vertical mesh information
      type (mpas_pool_type), intent(in) :: scratchPool !< Input: scratch variables
      type (mpas_pool_type), intent(in) :: tracersPool !< Input: tracer fields
      integer, intent(in), optional :: timeLevelIn !< Input: Time level in state
      logical, intent(in), optional :: full   !< Input: Run all computations?

      integer :: iEdge, iCell, iTracer, k, kmin, kmax
      integer :: err, nCells, nEdges, nVertices, numTracers
      integer :: timeLevel

      real (kind=RKIND) :: coef_3rd_order

      real (kind=RKIND), dimension(:), pointer :: &
        ssh,  seaIcePressure, atmosphericPressure
      real (kind=RKIND), dimension(:,:), pointer :: &
        layerThickness, layerThicknessLag, normalVelocity, restingThickness

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      integer, pointer :: indexTemperaturePtr, indexSalinityPtr
      integer :: indexTemperature, indexSalinity

      logical :: full_compute = .true.

      real (kind=RKIND), dimension(:), pointer :: &
        frazilSurfacePressure, landIcePressure, landIceDraft, landIceFraction

      integer, dimension(:), pointer :: &
        landIceFloatingMask

      call mpas_timer_start('diagnostic solve')

      if ( present(full) ) then
         full_compute = full
      else
         full_compute = .true.
      end if

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_dimension(tracersPool, 'index_temperature', &
                                                 indexTemperaturePtr)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', &
                                                 indexSalinityPtr)
      indexTemperature = indexTemperaturePtr
      indexSalinity    = indexSalinityPtr
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)

      call mpas_pool_get_array(forcingPool, 'seaIcePressure', seaIcePressure)
      call mpas_pool_get_array(forcingPool, 'atmosphericPressure', atmosphericPressure)
      call mpas_pool_get_array(forcingPool, 'frazilSurfacePressure', frazilSurfacePressure)
      call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

      if (landIcePressureOn) then
         call mpas_pool_get_array(forcingPool, 'landIcePressure', landIcePressure)
         call mpas_pool_get_array(forcingPool, 'landIceDraft', landIceDraft)
         call mpas_pool_get_array(forcingPool, 'landIceFloatingMask', landIceFloatingMask)
         call mpas_pool_get_array(forcingPool, 'landIceFraction', landIceFraction)
         if (nCellsAll.gt.0 .AND. landIceFloatingMask(1) == -1) then
            call mpas_log_write('landIceFloatingMask contains the default value which likely indicates that this field is missing in the initial condition file (e.g. because it is meant for an older E3SM version).', &
                                MPAS_LOG_CRIT)
         endif
      end if

      nEdges = nEdgesAll
      nCells = nCellsAll
      nVertices = nVerticesAll

      coef_3rd_order = config_coef_3rd_order

#ifdef MPAS_OPENACC
      !! Listing new inputs required by the first routine where they are used
      !! ocn_diagnostic_solve_layerThicknessEdge :: layerThickness, (not on device)
      !!                                            normalVelocity (not on device)
      !! ocn_relativeVorticity_circulation       ::
      !! ocn_diagnostic_solve_vortVel            :: normalTransportVelocity, (diagnostics array)
      !!                                            normalGMBolusVelocity (diagnostics array)
      !! ocn_diagnostic_solve_kineticEnergy      ::
      !! ocn_diagnostic_solve_vorticity          ::
      !! ocn_diagnostic_solve_surface_pressure   :: atmosphericPressure, (not on device)
      !!                                            seaIcePressure, (not on device)
      !!                                            frazilSurfacePressure, (not on device)
      !!                                            landIcePressure (not on device)
      !! ocn_diagnostic_solve_z_coordinates      ::
      !! ocn_equation_of_state_density           :: activeTracers, (not on device)
      !!                                            tracersSurfaceValue, (diagnostics array)
      !!                                            maxLevelCell, (mesh array)
      !!                                            refBottomDepth  (mesh array)
      !! ocn_diagnostic_solve_pressure           ::
      !! ocn_diagnostic_solve_richardson         :: edgeAreaFractionOfCell (mesh array)

      !! tracersSurfaceValue calc                :: minLevelCell (mesh array)
      !! normalVelocitySurfaceLayer calc         :: minLevelEdgeBot (mesh array)
      !! surfaceFluxAttenuationCoefficientRunoff :: layerThickEdgeFlux (diagnostics array)
      !! surfaceFluxAttenuationCoefficientSubglacialRunoff :: layerThickEdgeFlux (diagnostics array)
      !! ocn_diagnostic_solve_ssh                :: ssh, (not on device)
      !!                                            landIceDraft (not on device)

!     !$acc enter data copyin(layerThickness, normalVelocity)
!     !$acc update device (normalTransportVelocity, &
!     !$acc                normalGMBolusVelocity)

!     !$acc enter data copyin(atmosphericPressure, seaIcePressure)

!     !$acc enter data copyin(activeTracers, ssh)
!     !$acc enter data copyin(activeTracers)
!     !$acc update device(tracersSurfaceValue)

!     if ( associated(frazilSurfacePressure) ) then
!        !$acc enter data copyin(frazilSurfacePressure)
!     endif
!     if (landIcePressureOn) then
!        !$acc enter data copyin(landIcePressure)
!        !$acc enter data copyin(landIceDraft)
!     endif
#endif

      !
      ! Z-coordinates
      ! This section must be placed in the code before computing the density.
      !
      ! inputs : layerThickness
      ! outputs : zMid, zTop, ssh)
      call ocn_diagnostic_solve_z_coordinates(layerThickness, zMid, zTop, ssh)

      ! inputs: layerThickness, normalVelocity
      ! output: layerThickEdgeMean, layerThickEdgeDrag, layerThickEdgeFlux
      call ocn_diagnostic_solve_layerThicknessEdge(normalVelocity, &
                                                   layerThickness, restingThickness, &
                                                   ssh)

      ! inputs: normalVelocity
      ! outputs: relativeVorticity, circulation
      call ocn_relativeVorticity_circulation(relativeVorticity, circulation, normalVelocity, err)

      ! inputs: relativeVorticity, layerThickEdgeFlux, normalVelocity, normalTransportVelocity,
      !         normalGMBolusVelocity
      ! outputs: vertTransportVelocityTop, vertGMBolusVelocityTop, relativeVorticityCell,
      !          divergence, kineticEnergyCell, tangentialVelocity, vertMLEBolusVelocityTop
      ! in and out: vertVelocityTop
      call ocn_diagnostic_solve_vortVel(relativeVorticity, &
             layerThickEdgeFlux, normalVelocity, normalTransportVelocity, &
             vertVelocityTop, vertTransportVelocityTop, relativeVorticityCell, &
             divergence, kineticEnergyCell, tangentialVelocity)

      !compute vertical velocity due to gm bolus velocity if requested
      if (config_use_gm) then
         call ocn_diagnostic_solve_GMvel(layerThickEdgeFlux, normalGMBolusVelocity, vertGMBolusVelocityTop)
      end if

      !compute vertical velocity due to submesoscale eddies (MLE) velocity if requested
      if (config_submesoscale_enable) then
         call ocn_diagnostic_solve_MLEvel(layerThickEdgeFlux, normalMLEVelocity, vertMLEBolusVelocityTop)
      end if

      if ( configVertAdvMethod == vertAdvRemap ) then

         ! Overwrite vertVelocityTop with exact Eulerian solution from remapping
         call mpas_pool_get_array(statePool, 'layerThicknessLag', &
             layerThicknessLag, timeLevel)
         call ocn_diagnostic_solve_vertVel_remap(layerThickness, &
             layerThicknessLag, dt, vertVelocityTop)

         ! Overwrite vertical tracer budget terms with new vertVelocityTop
         if ( config_compute_active_tracer_budgets ) then
            numTracers = size(activeTracers, dim=1)
            do iTracer = 1, numTracers
               !$omp parallel
               !$omp do schedule(runtime) private(k,kmin,kmax)
               do iCell = 1, nCells
                  kmin = minLevelCell(iCell)
                  kmax = maxLevelCell(iCell)

                  activeTracerVerticalAdvectionTopFlux(iTracer, :, iCell) = 0.0_RKIND
                  do k = kmin+1, kmax
                     activeTracerVerticalAdvectionTopFlux(iTracer, k, iCell) = &
                           min(0.0_RKIND,vertVelocityTop(k, iCell))* &
                                         activeTracers(iTracer, k-1, iCell) &
                         + max(0.0_RKIND,vertVelocityTop(k,iCell))* &
                                         activeTracers(iTracer, k, iCell)
                  enddo

                  do k = kmin, kmax-1
                     activeTracerVerticalAdvectionTendency(iTracer, k, iCell) = ( &
                          activeTracerVerticalAdvectionTopFlux(iTracer, k+1, iCell) &
                        - activeTracerVerticalAdvectionTopFlux(iTracer, k, iCell) &
                        ) / layerThickness(k, iCell)
                  enddo
                  activeTracerVerticalAdvectionTendency(iTracer, kmax, iCell) = &
                     - activeTracerVerticalAdvectionTopFlux(iTracer, kmax, iCell) &
                     / layerThickness(kmax, iCell)
               enddo
               !$omp end do
               !$omp end parallel
            enddo
         endif
      endif

      !
      ! Compute kinetic energy
      !
      ! inputs :: normalVelocity
      ! outputs :: kineticEnergyCell
      if ( ke_vertex_flag == 1 ) then
         call ocn_diagnostic_solve_kineticEnergy(normalVelocity, kineticEnergyCell)
      end if

      ! inputs :: normalVelocity, tangentialVelocity, layerThickness, relativeVorticity
      ! outputs :: normRelVortEdge, normPlanetVortEdge,
      !            normalizedRelativeVorticityCell
      call ocn_diagnostic_solve_vorticity(dt, normalVelocity, tangentialVelocity, &
                layerThickness, relativeVorticity, &
                normRelVortEdge, normPlanetVortEdge, &
                normalizedRelativeVorticityCell)

      !
      ! Surface pressure
      ! This section must be placed in the code before computing the density.
      !
      ! inputs: atmosphericPressure, seaIcePressure, frazilSurfacePressure, landIcePressure
      ! outputs: surfacePressure
      call ocn_diagnostic_solve_surface_pressure(forcingPool, atmosphericPressure, &
                seaIcePressure, surfacePressure)


      !
      ! equation of state
      !

      ! Only need densities on 0, 1, and 2 halo cells
      nCells = nCellsAll

      ! inputs: activeTracers, tracersSurfaceValue,maxLevelCell,refBottomDepth
      ! outputs: density, potentialDensity, displacedDensity, thermalExpansionCoeff, salineContractionCoeff
      ! compute in-place density
      call ocn_equation_of_state_density(statePool, meshPool, activeTracers, tracersSurfaceValue, &
                                         nCells, 0, 'relative', density, err, &
                                         thermExpCoeff, salineContractCoeff, &
                                         timeLevelIn=timeLevel)
      ! compute potentialDensity, the density displaced adiabatically to the mid-depth of top layer.
      call ocn_equation_of_state_density(statePool, meshPool, activeTracers, tracersSurfaceValue, &
                                         nCells, 1, 'absolute', potentialDensity, &
                                         err, timeLevelIn=timeLevel)

      ! compute displacedDensity, density displaced adiabatically to the mid-depth one layer deeper.
      ! That is, layer k has been displaced to the depth of layer k+1.
      call ocn_equation_of_state_density(statePool, meshPool, activeTracers, tracersSurfaceValue, &
                                         nCells, 1, 'relative', displacedDensity, &
                                         err, timeLevelIn=timeLevel)

      !
      ! Pressure
      ! This section must be placed in the code after computing the density.
      !
      ! inputs: layerThickness, density, surfacePressure
      ! outputs: montgomeryPotential, pressure
      call ocn_diagnostic_solve_pressure(layerThickness, density, &
                surfacePressure, montgomeryPotential, pressure)

      nCells = nCellsAll

      if ( full_compute ) then

         !
         ! Brunt-Vaisala frequency and gradient Richardson number
         !
         ! inputs: displacedDensity, density, zMid, normalVelocity, edgeAreaFractionOfCell
         ! outputs: BruntVaisalaFreqTop, RiTopOfCell
         call ocn_diagnostic_solve_richardson(displacedDensity, density, zMid, &
                normalVelocity, edgeAreaFractionOfCell, BruntVaisalaFreqTop, RiTopOfCell)

      end if    ! full_compute

      !
      ! extrapolate tracer values to ocean surface
      ! this eventually be a modelled process
      ! at present, just copy k=1 tracer values onto surface values
      ! field will be updated below is better approximations are available

!TDR need to consider how to handle tracersSurfaceValues
#ifdef MPAS_OPENACC
      ! tracersSurfaceValue calc
      ! inputs: activeTracers, minLevelCell
      ! outputs: tracersSurfaceValue
      !$acc kernels present(tracersSurfaceValue, activeTracers, minLevelCell)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
#ifdef MPAS_OPENACC
      !$acc loop independent
#endif
      do iCell = 1, nCells
         do iTracer = 1,size(tracersSurfaceValue,1)
            tracersSurfaceValue(iTracer, iCell) = activeTracers(iTracer, minLevelCell(iCell), iCell)
         enddo
      end do
#ifdef MPAS_OPENACC
      !$acc end kernels
#else
      !$omp end do
#endif

#ifdef MPAS_OPENACC
      ! normalVelocitySurfaceLayer calc
      ! inputs: normalVelocity, minLevelEdgeBot
      ! outputs: normalVelocitySurfaceLayer
      !$acc parallel loop gang vector &
      !$acc          present(normalVelocity, normalVelocitySurfaceLayer, minLevelEdgeBot)
#else
      !$omp do schedule(runtime)
#endif
      do iEdge = 1, nEdges
         normalVelocitySurfaceLayer(iEdge) = normalVelocity(minLevelEdgeBot(iEdge), iEdge)
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      if ( full_compute ) then

         !
         ! average tracer values over the ocean surface layer
         ! the ocean surface layer is generally assumed to be about 0.1 of the boundary layer depth
         !
         ! inputs: layerThickness, activeTracers, layerThicknessEdgeFlux, normalVelocity
         ! outputs: tracersSurfaceLayerValue, indexSurfaceLayerDepth, normalVelocitySurfaceLayer,
         !          surfaceFluxAttenuationCoefficient, surfaceFluxAttenuationCoefficientRunoff
         call ocn_diagnostic_solve_surfaceLayer(layerThickness, activeTracers, &
                layerThickEdgeFlux, normalVelocity, tracersSurfaceLayerValue, &
                indexSurfaceLayerDepth, normalVelocitySurfaceLayer, &
                sfcFlxAttCoeff, surfaceFluxAttenuationCoefficientRunoff)

         if (trim(config_subglacial_runoff_mode) == 'data') then
            ! outputs: surfaceFluxAttenuationCoefficientSubglacialRunoff
            call ocn_diagnostic_solve_surfaceLayer_subglacialRunoff(surfaceFluxAttenuationCoefficientSubglacialRunoff)
         end if

      end if ! full_compute

      !
      !  compute fields needed to compute land-ice fluxes, either in the ocean model or in the coupler
      !
      ! inputs: layerThickness, normalVelocity, landIceFraction, landIceFloatingMask
      ! outputs:
      if (landIcePressureOn) then
         call ocn_compute_land_ice_flux_input_fields(layerThickness, normalVelocity, activeTracers, &
                landIceFraction, landIceFloatingMask, timeLevel, indexTemperature, indexSalinity)
      end if

      if ( full_compute ) then

         !
         ! inputs: ssh, seaIcePressure, landIceDraft
         ! outputs: pressureAdjustedSSH, gradSSH
         if (landIcePressureOn) then
            call ocn_diagnostic_solve_ssh(forcingPool, ssh, seaIcePressure, &
                   pressureAdjustedSSH, gradSSH, landIceDraft=landIceDraft )
        else
            call ocn_diagnostic_solve_ssh(forcingPool, ssh, seaIcePressure, &
                   pressureAdjustedSSH, gradSSH)
         end if

      end if ! full_compute

#ifdef MPAS_OPENACC
      !! Listing outputs by first routine that modifies variable
      !! ocn_diagnostic_solve_layerThicknessEdge :: layerThickEdgeDrag, (diagnostics array)
      !!                                            layerThickEdgeFlux, (diagnostics array)
      !!                                            layerThickEdgeMean (diagnostics array)
      !! ocn_relativeVorticity_circulation       :: relativeVorticity, (diagnostics array)
      !!                                            circulation (diagnostics array)
      !! ocn_diagnostic_solve_vortVel            :: vertTransportVelocityTop, (diagnostics array)
      !!                                            vertGMBolusVelocityTop, (diagnostics array)
      !!                                            relativeVorticityCell, (diagnostics array)
      !!                                            divergence, (diagnostics array)
      !!                                            kineticEnergyCell, (diagnostics array)
      !!                                            tangentialVelocity (diagnostics array)
      !!                                            vertVelocityTop (diagnostics array)
      !! ocn_diagnostic_solve_kineticEnergy      ::
      !! ocn_diagnostic_solve_vorticity          :: normRelVortEdge, (diagnostics array)
      !!                                            normPlanetVortEdge, (diagnostics array)
      !!                                            normalizedRelativeVorticityCell (diagnostics array)
      !! ocn_diagnostic_solve_surface_pressure   :: surfacePressure (diagnostics array)
      !! ocn_diagnostic_solve_z_coordinates      :: zMid, (diagnostics array)
      !!                                            zTop,(diagnostics array)
      !!                                            ssh (not on device)
      !! ocn_equation_of_state_density           :: density, (diagnostics array)
      !!                                            potentialDensity, (diagnostics array)
      !!                                            displacedDensity, (diagnostics array)
      !!                                            thermalExpansionCoeff, (diagnostics array)
      !!                                            salineContractionCoeff (diagnostics array)
      !! ocn_diagnostic_solve_pressure           :: montgomeryPotential, (diagnostics array)
      !!                                            pressure (diagnostics array)
      !! ocn_diagnostic_solve_richardson         :: BruntVaisalaFreqTop, (diagnostics array)
      !!                                            RiTopOfCell (diagnostics array)
      !! tracersSurfaceValue calc                :: tracersSurfaceValue (diagnostics array)
      !! normalVelocitySurfaceLayer calc         :: normalVelocitySurfaceLayer (diagnostics array)
      !! surfaceFluxAttenuationCoefficientRunoff :: tracersSurfaceLayerValue, (diagnostics array)
      !!                                            indexSurfaceLayerDepth, (diagnostics array)
      !!                                            normalVelocitySurfaceLayer, (diagnostics array)
      !!                                            surfaceFluxAttenuationCoefficient, ????
      !!                                            surfaceFluxAttenuationCoefficientRunoff (diagnostics array)
      !! surfaceFluxAttenuationCoefficientSubglacialRunoff :: tracersSurfaceLayerValue, (diagnostics array)
      !!                                            indexSurfaceLayerDepth, (diagnostics array)
      !!                                            normalVelocitySurfaceLayer, (diagnostics array)
      !!                                            surfaceFluxAttenuationCoefficient, ????
      !!                                            surfaceFluxAttenuationCoefficientSubglacialRunoff (diagnostics array)
      !! ocn_diagnostic_solve_ssh                :: pressureAdjustedSSH, (diagnostics array)
      !!                                            gradSSH (diagnostics array)

!     !$acc update host(layerThickEdgeDrag, &
!     !$acc             layerThickEdgeFlux, &
!     !$acc             layerThickEdgeMean)
!     !$acc update host(relativeVorticity, circulation)
!     !$acc update host(vertTransportVelocityTop, &
!     !$acc             vertGMBolusVelocityTop, &
!     !$acc             relativeVorticityCell, &
!     !$acc             divergence, &
!     !$acc             kineticEnergyCell, &
!     !$acc             tangentialVelocity, &
!     !$acc             vertVelocityTop)
!     !$acc update host(normRelVortEdge, normPlanetVortEdge, &
!     !$acc             normalizedRelativeVorticityCell)
!     !$acc update host (surfacePressure)
!     !$acc update host(zMid, zTop)
!     !$acc exit data copyout(ssh)
!     !$acc update host(density, potentialDensity, displacedDensity)
!     !$acc update host(thermExpCoeff,  &
!     !$acc&            salineContractCoeff)
!     !$acc update host(tracersSurfaceValue)

!     !$acc update host(normalVelocitySurfaceLayer)
!     !$acc update host(montgomeryPotential, pressure)

!     if ( full_compute ) then
!        !$acc update host(RiTopOfCell, &
!        !$acc             BruntVaisalaFreqTop)

!        !$acc update host(tracersSurfaceLayerValue, &
!        !$acc             indexSurfaceLayerDepth, &
!        !$acc             normalVelocitySurfaceLayer, &
!        !$acc             sfcFlxAttCoeff, &
!        !$acc             surfaceFluxAttenuationCoefficientRunoff)
!        if (trim(config_subglacial_runoff_mode) == 'data') then
!        !$acc update host(surfaceFluxAttenuationCoefficientSubglacialRunoff)
!        end if

!     end if ! full_compute

!     !$acc exit data delete (atmosphericPressure, seaIcePressure)
!     if ( associated(frazilSurfacePressure) ) then
!        !$acc exit data delete(frazilSurfacePressure)
!     endif
!     if (landIcePressureOn) then
!        !$acc exit data delete(landIcePressure)
!        !$acc exit data delete(landIceDraft)
!     endif

!     !$acc exit data delete (layerThickness, normalVelocity, &
!     !$acc                   activeTracers)
#endif

      call mpas_timer_stop('diagnostic solve')

   end subroutine ocn_diagnostic_solve!}}}

!***********************************************************************
!
!  routine ocn_relativeVorticity_circulation
!
!> \brief   Computes relative vorticity and circulation
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    November 2013
!> \details
!>  Computes relative vorticity and circulation
!
!-----------------------------------------------------------------------

   subroutine ocn_relativeVorticity_circulation(relativeVorticity, circulation, normalVelocity, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         relativeVorticity

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         circulation

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iVertex, iEdge, i, k

      real (kind=RKIND) :: invAreaTri1, r_tmp

      err = 0

#ifdef MPAS_OPENACC
      !$acc parallel loop collapse(2) present(circulation, relativeVorticity)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
      do iVertex = 1, nVerticesAll
         do k = 1, nVertLevels
            circulation(k, iVertex) = 0.0_RKIND
            relativeVorticity(k, iVertex) = 0.0_RKIND
         enddo
      enddo
#ifndef MPAS_OPENACC
      !$omp end do
#endif

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc          present(circulation, relativeVorticity, areaTriangle, edgesOnVertex, &
      !$acc                  maxLevelVertexBot, dcEdge, normalVelocity, edgeSignOnVertex, &
      !$acc                  minLevelVertexTop) &
      !$acc          private(invAreaTri1, iVertex, iEdge, k, r_tmp)
#else
      !$omp do schedule(runtime) private(invAreaTri1, i, iEdge, k, r_tmp)
#endif
      do iVertex = 1, nVerticesAll
         invAreaTri1 = 1.0_RKIND / areaTriangle(iVertex)
         do i = 1, vertexDegree
            iEdge = edgesOnVertex(i, iVertex)
            do k = minLevelVertexTop(iVertex), maxLevelVertexBot(iVertex)
              r_tmp = dcEdge(iEdge) * normalVelocity(k, iEdge)
              circulation(k, iVertex) = circulation(k, iVertex) + edgeSignOnVertex(i, iVertex) * r_tmp
              relativeVorticity(k, iVertex) = relativeVorticity(k, iVertex) + edgeSignOnVertex(i, iVertex) * r_tmp * invAreaTri1
            end do
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

   !--------------------------------------------------------------------

   end subroutine ocn_relativeVorticity_circulation!}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve_wctEdge
!
!> \brief   Computes layer thickness at edge
!> \author  Carolyn Begeman
!> \date    July 2023
!> \details
!>  This routine computes the water column thickness on edges, wctEdge, which
!>  can either be the mean or an upwind value, depending on user input.
!
!-----------------------------------------------------------------------

   subroutine ocn_diagnostic_solve_wctEdge(btrVelocity, deltaSSH, &
                         baroclinicThickness, wctEdge)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
         btrVelocity,       &!< barotropic velocity
         deltaSSH,          &!< ssh change from baroclinic step to barotropic
                             !< substep
         baroclinicThickness !< water column thickness on edges at the last
                             !< baroclinic step (centered or upwinded)

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      ! Outputs are the shared variables from diagnostic_variables:
      real (kind=RKIND), dimension(:), intent(inout) :: &
         wctEdge !< water column thickness on edge

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         iCell, iEdge, &! edge and cell loop indices
         cell1, cell2  ! neighbor cell indices across an edge


      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      ! Compute edge flux based on option set on init
      select case (thickEdgeFluxChoice)

      case (thickEdgeFluxCenter)
         ! Use centered (mean) thickness as flux value

         ! Compute mean layer thickness at edge
#ifdef MPAS_OPENACC
         !$acc enter data copyin(baroclinicThickness, deltaSSH, wctEdge)
         !$acc parallel loop &
         !$acc    present(baroclinicThickness, deltaSSH, wctEdge, cellsOnEdge) &
         !$acc    private(iEdge, cell1, cell2)
#else
         !$omp parallel
         !$omp do schedule(runtime) private (cell1, cell2)
#endif
         do iEdge = 1, nEdgesAll
            if ( maxLevelEdgeTop(iEdge) == 0 ) cycle
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            wctEdge(iEdge) = 0.5_RKIND*(deltaSSH(cell1) + deltaSSH(cell2)) &
                             + baroclinicThickness(iEdge)
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

#ifdef MPAS_OPENACC
         !$acc exit data copyout(wctEdge) delete(baroclinicThickness, deltaSSH)
#endif

      case (thickEdgeFluxUpwind)

         ! Use upwind thickness as the edge flux value
#ifdef MPAS_OPENACC
         !$acc enter data copyin(baroclinicThickness, deltaSSH, btrVelocity, wctEdge)
         !$acc parallel loop &
         !$acc    present(baroclinicThickness, btrVelocity, deltaSSH, cellsOnEdge, maxLevelEdgeTop, wctEdge) &
         !$acc    private(cell1, cell2, iEdge)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(cell1, cell2)
#endif
         do iEdge = 1, nEdgesAll
            if ( maxLevelEdgeTop(iEdge) == 0 ) cycle
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            if (btrVelocity(iEdge) > 0.0_RKIND) then
               wctEdge(iEdge) = deltaSSH(cell1) + baroclinicThickness(iEdge)
            elseif (btrVelocity(iEdge) < 0.0_RKIND) then
               wctEdge(iEdge) = deltaSSH(cell2) + baroclinicThickness(iEdge)
            else
               wctEdge(iEdge) = max(deltaSSH(cell1), deltaSSH(cell2)) &
                                + baroclinicThickness(iEdge)
            end if
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

#ifdef MPAS_OPENACC
         !$acc exit data copyout(wctEdge) delete(baroclinicThickness, deltaSSH, btrVelocity)
#endif

      case default
         ! Should have been caught on init
         call mpas_log_write('Thickness flux option unknown', &
                             MPAS_LOG_CRIT)

      end select

   end subroutine ocn_diagnostic_solve_wctEdge!}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve_layerThicknessEdge
!
!> \brief   Computes layer thickness at edge
!> \author  Matt Turner
!> \date    October 2020
!> \details
!>  This routine computes the diagnostic variables layerThickEdgeDrag,
!>  layerThickEdgeFlux, and layerThickEdgeMean. The Mean is a simple
!>  mean across the edge but the Flux value can either be the mean or
!>  an upwind value, depending on user input.
!
!-----------------------------------------------------------------------

   subroutine ocn_diagnostic_solve_layerThicknessEdge(normalVelocity, &
                                                      layerThickness, restingThickness, &
                                                      ssh)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity,  &!< [in] transport
         layerThickness,  &!< [in] layer thickness at cell center
         restingThickness  !< [in] initial layer thickness at cell center
      real (kind=RKIND), dimension(:), intent(in) :: &
         ssh               !< [in] sea surface height at cell center

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      ! Outputs are the shared variables from diagnostic_variables:
      !real (kind=RKIND), dimension(:,:), intent(out) :: &
      !   layerThickEdgeDrag,    & !< [out] layer thickness on edge for drag
      !   layerThickEdgeMean,    & !< [out] centered layer thickness on edge
      !   layerThickEdgeFlux       !< [out] flux-related thickness on edge

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         iEdge, k,     &! edge and vertical loop indices
         cell1, cell2, &! neighbor cell indices across an edge
         kmin, kmax     ! min,max active vertical levels

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      ! Compute mean layer thickness at edge
#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(layerThickness, layerThickEdgeMean, &
      !$acc            minLevelEdgeBot, maxLevelEdgeTop, cellsOnEdge) &
      !$acc    private(k, kmin, kmax, cell1, cell2)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k, kmin, kmax, cell1, cell2)
#endif
      do iEdge = 1, nEdgesAll
         kmin = minLevelEdgeBot(iEdge)
         kmax = maxLevelEdgeTop(iEdge)
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = 1,nVertLevels
            ! initialize layerThicknessEdgeMean to avoid divide by
            ! zero and NaN problems.
            layerThickEdgeMean(k,iEdge) = -1.0e34_RKIND
         end do
         do k = kmin,kmax
            ! central differenced
            layerThickEdgeMean(k,iEdge) = 0.5_RKIND * &
                                         (layerThickness(k,cell1) + &
                                          layerThickness(k,cell2))
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      ! Compute edge flux based on option set on init
      select case (thickEdgeFluxChoice)

      case (thickEdgeFluxCenter)
         ! Use centered (mean) thickness as flux value

         if ( config_use_subgrid_wetting_drying ) then
#ifdef MPAS_OPENACC
         !$acc update host(ssh)
#endif

            call ocn_subgrid_layerThickEdgeFlux_center(ssh, layerThickEdgeFlux)

#ifdef MPAS_OPENACC
         !$acc update device(layerThickEdgeFlux)
#endif
         else

#ifdef MPAS_OPENACC
            !$acc parallel loop collapse(2) &
            !$acc    present(layerThickEdgeFlux, layerThickEdgeMean)
#else
            !$omp parallel
            !$omp do schedule(runtime) private(k)
#endif
            do iEdge = 1, nEdgesAll
            do k = 1,nVertLevels
               layerThickEdgeFlux(k,iEdge) = &
               layerThickEdgeMean(k,iEdge)
            end do
            end do
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif

         end if  

      case (thickEdgeFluxUpwind)
         ! Use upwind thickness as the edge flux value

         if (config_use_subgrid_wetting_drying) then
#ifdef MPAS_OPENACC
         !$acc update host(ssh, normalVelocity, layerThickness)
#endif

            call ocn_subgrid_layerThickEdgeFlux_upwind(ssh, normalVelocity, layerThickness, layerThickEdgeFlux)

#ifdef MPAS_OPENACC
         !$acc update device(layerThickEdgeFlux)
#endif

         else
#ifdef MPAS_OPENACC
            !$acc parallel loop &
            !$acc    present(normalVelocity, layerThickness, &
            !$acc            minLevelEdgeBot, maxLevelEdgeTop, &
            !$acc            layerThickEdgeFlux, cellsOnEdge) &
            !$acc    private(k, kmin, kmax, cell1, cell2)
#else
            !$omp parallel
            !$omp do schedule(runtime) private(k, kmin, kmax, cell1, cell2)
#endif
            do iEdge = 1, nEdgesAll
               kmin = minLevelEdgeBot(iEdge)
               kmax = maxLevelEdgeTop(iEdge)
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               do k=1,nVertLevels
                  ! initialize layerThicknessEdgeFlux to avoid divide by
                  ! zero and NaN problems.
                  layerThickEdgeFlux(k,iEdge) = -1.0e34_RKIND
               end do
               do k = kmin,kmax
                  if (normalVelocity(k,iEdge) > 0.0_RKIND) then
                     layerThickEdgeFlux(k,iEdge) = layerThickness(k,cell1)
                  elseif (normalVelocity(k,iEdge) < 0.0_RKIND) then
                     layerThickEdgeFlux(k,iEdge) = layerThickness(k,cell2)
                  else
                     layerThickEdgeFlux(k,iEdge) = &
                                         max(layerThickness(k,cell1), &
                                             layerThickness(k,cell2))
                  end if
               end do
            end do
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif

        end if

      case (thickEdgeFluxConstant)
         ! Use linearized version H*u where H is constant in time

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(restingThickness, &
         !$acc            minLevelEdgeBot, maxLevelEdgeTop, &
         !$acc            layerThickEdgeFlux, cellsOnEdge) &
         !$acc    private(k, kmin, kmax, cell1, cell2)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k, kmin, kmax, cell1, cell2)
#endif
         do iEdge = 1, nEdgesAll
            kmin = minLevelEdgeBot(iEdge)
            kmax = maxLevelEdgeTop(iEdge)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1,nVertLevels
               ! initialize layerThicknessEdgeFlux to avoid divide by
               ! zero and NaN problems.
               layerThickEdgeFlux(k,iEdge) = -1.0e34_RKIND
            end do
            do k = kmin,kmax
               layerThickEdgeFlux(k,iEdge) = 0.5_RKIND * &
                                         (restingThickness(k,cell1) + &
                                          restingThickness(k,cell2))
            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      case default
         ! Should have been caught on init
         call mpas_log_write('Thickness flux option unknown', &
                             MPAS_LOG_CRIT)

      end select

      ! Compute edge flux based on option set on init
      select case (thickEdgeDragChoice)

      case (thickEdgeDragCenter)
#ifdef MPAS_OPENACC
         !$acc parallel loop collapse(2) &
         !$acc    present(layerThickEdgeDrag, layerThickEdgeMean)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iEdge = 1, nEdgesAll
         do k = 1,nVertLevels
            layerThickEdgeDrag(k,iEdge) = &
            layerThickEdgeMean(k,iEdge)
         end do
         end do
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif

      case (thickEdgeDragHarMean)

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(layerThickness, layerThickEdgeDrag, &
         !$acc            minLevelEdgeBot, maxLevelEdgeTop, cellsOnEdge) &
         !$acc    private(k, kmin, kmax, cell1, cell2)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k, kmin, kmax, cell1, cell2)
#endif
         do iEdge = 1, nEdgesAll
            kmin = minLevelEdgeBot(iEdge)
            kmax = maxLevelEdgeTop(iEdge)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k = 1,nVertLevels
               ! initialize layerThicknessEdgeMean to avoid divide by
               ! zero and NaN problems.
               layerThickEdgeDrag(k,iEdge) = -1.0e34_RKIND
            end do
            do k = kmin,kmax
               ! central differenced
               layerThickEdgeDrag(k,iEdge) = &
                   2.0_RKIND * layerThickness(k,cell1) * layerThickness(k,cell2)/ &
                   (max(layerThickness(k,cell1)+layerThickness(k,cell2), 1e-9_RKIND))
            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      case default
         ! Should have been caught on init
         call mpas_log_write('Thickness flux option unknown', &
                             MPAS_LOG_CRIT)

      end select


   !--------------------------------------------------------------------

   end subroutine ocn_diagnostic_solve_layerThicknessEdge!}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve_kineticEnergy
!
!> \brief   Computes kinetic energy at cell center
!> \author  Matt Turner
!> \date    October 2020
!> \details
!>  This routine computes the kinetic energy (m^2/s^2) at cell center
!
!-----------------------------------------------------------------------

   subroutine ocn_diagnostic_solve_kineticEnergy(normalVelocity, &
                                                 kineticEnergyCell)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity     !< [in] velocity normal to edge

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         kineticEnergyCell  !< [out] kinetic energy at cell center

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         i, j, k,               &! generic loop indices
         iCell, iEdge, iVertex, &! indices for cell, edge, vertex
         nCells                  ! number of cells and vertices

      real (kind=RKIND) :: &
         rTmp,        &! local tmp for common factors
         invAreaCell1  ! 1/cell area

      ! Scratch Arrays
      real (kind=RKIND), dimension(:,:), allocatable :: &
         kineticEnergyVertex, &! horiz kinetic energy at vertices
         kineticEnergyVertexOnCells

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      nCells = nCellsHalo( 2 )

      allocate(kineticEnergyVertex(nVertLevels, nVerticesAll), &
               kineticEnergyVertexOnCells(nVertLevels, nCells))
      !$acc enter data create(kineticEnergyVertex, &
      !$acc                   kineticEnergyVertexOnCells)

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(kineticEnergyVertex, normalVelocity, &
      !$acc            areaTriangle, dcEdge, dvEdge, edgesOnVertex) &
      !$acc    private(iEdge,rTmp,k,i)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(i, k, iEdge, rTmp)
#endif
      do iVertex = 1, nVerticesAll
         kineticEnergyVertex(:, iVertex) = 0.0_RKIND
         !$acc loop seq
         do i = 1, vertexDegree
            iEdge = edgesOnVertex(i, iVertex)
            rTmp  = dcEdge(iEdge)*dvEdge(iEdge)*0.25_RKIND/ &
                    areaTriangle(iVertex)
            do k = 1, nVertLevels
               kineticEnergyVertex(k,iVertex) = &
               kineticEnergyVertex(k,iVertex) + rTmp* &
                            normalVelocity(k,iEdge)**2
            end do
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(kineticEnergyVertexOnCells, kineticEnergyVertex,&
      !$acc            kiteIndexOnCell, kiteAreasOnVertex, invAreaCell,&
      !$acc            verticesOnCell, nEdgesOnCell) &
      !$acc    private(i,j,k,iVertex,invAreaCell1)
#else
      !$omp do schedule(runtime) &
      !$omp    private(i, j, k, iVertex, invAreaCell1)
#endif
      do iCell = 1, nCells
         kineticEnergyVertexOnCells(:, iCell) = 0.0_RKIND
         invAreaCell1 = invAreaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            j = kiteIndexOnCell(i, iCell)
            iVertex = verticesOnCell(i, iCell)
            !$acc loop seq
            do k = 1, nVertLevels
               kineticEnergyVertexOnCells(k,iCell) = &
               kineticEnergyVertexOnCells(k,iCell) + &
                        kiteAreasOnVertex(j,iVertex)* &
                        kineticEnergyVertex(k,iVertex)*invAreaCell1
            end do
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

      !
      ! Compute kinetic energy in each cell by blending
      ! kineticEnergyCell and kineticEnergyVertexOnCells
      !
#ifdef MPAS_OPENACC
      !$acc parallel loop collapse(2) &
      !$acc    present(kineticEnergyCell, kineticEnergyVertexOnCells)
#else
      !$omp do schedule(runtime) private(k)
#endif
      do iCell = 1, nCells
      do k = 1, nVertLevels
         kineticEnergyCell(k,iCell) = 5.0_RKIND/8.0_RKIND* &
                                      kineticEnergyCell(k,iCell) &
                                    + 3.0_RKIND/8.0_RKIND* &
                                      kineticEnergyVertexOnCells(k,iCell)
      end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !$acc exit data delete(kineticEnergyVertex, &
      !$acc                  kineticEnergyVertexOnCells)
      deallocate(kineticEnergyVertex, &
                 kineticEnergyVertexOnCells)

   !--------------------------------------------------------------------

   end subroutine ocn_diagnostic_solve_kineticEnergy!}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve_vorticity
!
!> \brief   Computes normalized vorticity at various locations
!> \author  Matt Turner
!> \date    October 2020
!> \details
!>  This routine computes a normalized vorticity at centers, edges
!>  based on relative vorticity at cell centers
!
!-----------------------------------------------------------------------

   subroutine ocn_diagnostic_solve_vorticity(dt, &
                         normalVelocity, tangentialVelocity, &
                         layerThickness, relativeVorticity, &
                         normalizedRelativeVorticityEdge, &
                         normalizedPlanetaryVorticityEdge, &
                         normalizedRelativeVorticityCell)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: dt !< [in] Time step

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity,     &!< [in] velocity normal to edge
         tangentialVelocity, &!< [in] velocity tangent to edge
         layerThickness,     &!< [in] layer thickness
         relativeVorticity    !< [in] relative vorticity at cell center

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         normalizedRelativeVorticityEdge,  &!< [out] rel vort at edge
         normalizedPlanetaryVorticityEdge, &!< [out] planet vort edge
         normalizedRelativeVorticityCell    !< [out] rel vort in cell

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         i, j, k,                   &! loop indices
         iCell, iEdge, iVertex,     &! indices for cell, edge, vertex
         nCells, nEdges, nVertices, &! number of cells, edges, vertices
         vertex1, vertex2,          &! neighbor vertex addresses
         cell1, cell2                ! neighbor cell   addresses

      real(kind=RKIND) :: &
         invAreaTri1,          &! 1/triangle area
         invAreaCell1,         &! 1/cell area
         invLength,            &! 1/length
         layerThicknessVertex, &! layer thickness at vertex
         apvm_scale_factor      ! scale factor for APVM form

      ! Local arrays
      ! normalizedPlanetaryVorticityVertex: earth's rotational rate
      !   (Coriolis parameter, f) divided by layer thickness,
      !   defined at vertices units: s^{-1}
      ! normalizedRelativeVorticityVertex: curl of horizontal velocity
      !   divided by layer thickness, defined at vertices units: s^{-1}
      ! vorticityGradientNormalComponent: gradient of vorticity in the
      !   normal direction (positive points from cell1 to cell2)
      !   units: s^{-1} m^{-1}
      ! vorticityGradientTangentialComponent: gradient of vorticity in
      !   the tangent direction (positive points from vertex1 to vertex2)
      !   units: s^{-1} m^{-1}

      real (kind=RKIND), dimension(:,:), allocatable :: &
         normalizedPlanetaryVorticityVertex,  &
         normalizedRelativeVorticityVertex,   &
         vorticityGradientNormalComponent,    &
         vorticityGradientTangentialComponent

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !
      ! Compute normalized relative and planetary vorticity
      !
      allocate(normalizedPlanetaryVorticityVertex(nVertLevels, nVerticesAll), &
               normalizedRelativeVorticityVertex(nVertLevels, nVerticesAll))
      !$acc enter data create(normalizedPlanetaryVorticityVertex, &
      !$acc                   normalizedRelativeVorticityVertex)

      nVertices = nVerticesHalo(2)

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc    present(normalizedRelativeVorticityVertex, &
      !$acc            normalizedPlanetaryVorticityVertex, &
      !$acc            relativeVorticity, fVertex, layerThickness, &
      !$acc            areaTriangle, kiteAreasOnVertex, cellsOnVertex, &
      !$acc            maxLevelVertexBot) &
      !$acc    private(i, k, iCell, invAreaTri1, layerThicknessVertex)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(i, k, iCell, invAreaTri1, layerThicknessVertex)
#endif
      do iVertex = 1, nVertices
         invAreaTri1 = 1.0_RKIND / areaTriangle(iVertex)
         do k = 1, maxLevelVertexBot(iVertex)
            layerThicknessVertex = 0.0_RKIND
            do i = 1, vertexDegree
               iCell = cellsOnVertex(i,iVertex)
               layerThicknessVertex = layerThicknessVertex &
                                    + layerThickness(k,iCell) &
                                    * kiteAreasOnVertex(i,iVertex)
            end do
            layerThicknessVertex = layerThicknessVertex*invAreaTri1
            if (layerThicknessVertex == 0) cycle

            normalizedRelativeVorticityVertex(k,iVertex) = &
                            relativeVorticity(k,iVertex) / &
                            layerThicknessVertex
            normalizedPlanetaryVorticityVertex(k,iVertex) = &
                            fVertex(iVertex)/layerThicknessVertex
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      nEdges = nEdgesHalo(2)

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc    present(normalizedRelativeVorticityEdge, &
      !$acc            normalizedPlanetaryVorticityEdge, &
      !$acc            normalizedRelativeVorticityVertex, &
      !$acc            normalizedPlanetaryVorticityVertex, &
      !&acc            minLevelEdgeTop, maxLevelEdgeBot, &
      !$acc            verticesOnEdge) &
      !$acc    private(k, vertex1, vertex2)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k, vertex1, vertex2)
#endif
      do iEdge = 1, nEdges
         normalizedRelativeVorticityEdge (:,iEdge) = 0.0_RKIND
         normalizedPlanetaryVorticityEdge(:,iEdge) = 0.0_RKIND
         vertex1 = verticesOnEdge(1, iEdge)
         vertex2 = verticesOnEdge(2, iEdge)
         !$acc loop seq
         do k = minLevelEdgeTop(iEdge), maxLevelEdgeBot(iEdge)
            normalizedRelativeVorticityEdge(k,iEdge) = 0.5_RKIND* &
                   (normalizedRelativeVorticityVertex(k,vertex1) &
                  + normalizedRelativeVorticityVertex(k,vertex2))
            normalizedPlanetaryVorticityEdge(k,iEdge) = 0.5_RKIND* &
                   (normalizedPlanetaryVorticityVertex(k,vertex1) &
                  + normalizedPlanetaryVorticityVertex(k,vertex2))
        end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      nCells = nCellsHalo(1)

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc    present(normalizedRelativeVorticityCell, &
      !$acc            normalizedRelativeVorticityVertex, &
      !$acc            invAreaCell, kiteAreasOnVertex, &
      !$acc            nEdgesOnCell, kiteIndexOnCell, &
      !$acc            minLevelCell, maxLevelCell, verticesOnCell) &
      !$acc    private(i, j, k, iVertex, invAreaCell1)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(i, j, k, iVertex, invAreaCell1)
#endif
      do iCell = 1, nCells
         normalizedRelativeVorticityCell(:,iCell) = 0.0_RKIND
         invAreaCell1 = invAreaCell(iCell)

         do i = 1, nEdgesOnCell(iCell)
            j = kiteIndexOnCell(i, iCell)
            iVertex = verticesOnCell(i, iCell)
            !$acc loop seq
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               normalizedRelativeVorticityCell(k,iCell) = &
               normalizedRelativeVorticityCell(k,iCell) + &
                                  kiteAreasOnVertex(j,iVertex)* &
                  normalizedRelativeVorticityVertex(k,iVertex)* &
                  invAreaCell1
            end do
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      nCells = nCellsAll
      nVertices = nVerticesAll
      nEdges = nEdgesAll
      apvm_scale_factor = config_apvm_scale_factor

      ! Diagnostics required for the Anticipated Potential Vorticity Method (apvm).
      if (config_apvm_scale_factor>1e-10) then

         allocate( &
            vorticityGradientNormalComponent    (nVertLevels,nEdges), &
            vorticityGradientTangentialComponent(nVertLevels,nEdges))
         !$acc enter data create(vorticityGradientNormalComponent, &
         !$acc                   vorticityGradientTangentialComponent)

         nEdges = nEdgesHalo( 1 )

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(vorticityGradientNormalComponent, &
         !$acc            vorticityGradientTangentialComponent, &
         !$acc            normalizedRelativeVorticityCell, &
         !$acc            normalizedRelativeVorticityVertex, &
         !$acc            dcEdge, dvEdge, cellsOnEdge, verticesOnedge, &
         !$acc            minLevelEdgeBot, maxLevelEdgeTop, &
         !$acc            minLevelEdgeTop, maxLevelEdgeBot) &
         !$acc    private(k, cell1, cell2, vertex1, vertex2, invLength)
#else
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp    private(k, cell1, cell2, vertex1, vertex2, invLength)
#endif
         do iEdge = 1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            vertex1 = verticesOnedge(1,iEdge)
            vertex2 = verticesOnedge(2,iEdge)

            ! Compute gradient of PV in normal direction

            invLength = 1.0_RKIND/dcEdge(iEdge)
            !$acc loop seq
            do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
               vorticityGradientNormalComponent(k,iEdge) = &
                  (normalizedRelativeVorticityCell(k,cell2) - &
                   normalizedRelativeVorticityCell(k,cell1))*invLength
            enddo

            ! Compute gradient of PV in the tangent direction

            invLength = 1.0_RKIND / dvEdge(iEdge)
            !$acc loop seq
            do k = minLevelEdgeTop(iEdge), maxLevelEdgeBot(iEdge)
              vorticityGradientTangentialComponent(k,iEdge) = &
                 (normalizedRelativeVorticityVertex(k,vertex2) - &
                  normalizedRelativeVorticityVertex(k,vertex1))*invLength
            enddo

         enddo
#ifndef MPAS_OPENACC
         !$omp end do
#endif

         !
         ! Modify PV edge with upstream bias.
         !
#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(normalizedRelativeVorticityEdge,      &
         !$acc            vorticityGradientNormalComponent,     &
         !$acc            vorticityGradientTangentialComponent, &
         !$acc            normalVelocity, tangentialVelocity,   &
         !$acc            minLevelEdgeTop, maxLevelEdgeBot)
#else
         !$omp do schedule(runtime) private(k)
#endif
         do iEdge = 1,nEdges
            !$acc loop seq
            do k = minLevelEdgeTop(iEdge), maxLevelEdgeBot(iEdge)
               normalizedRelativeVorticityEdge(k,iEdge) = &
               normalizedRelativeVorticityEdge(k,iEdge) - &
                    apvm_scale_factor*dt* &
                    (normalVelocity(k,iEdge)* &
                     vorticityGradientNormalComponent(k,iEdge)      &
                   + tangentialVelocity(k,iEdge)* &
                     vorticityGradientTangentialComponent(k,iEdge) )
            end do
         enddo
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

         !$acc exit data delete(vorticityGradientNormalComponent, &
         !$acc                  vorticityGradientTangentialComponent)
         deallocate(vorticityGradientNormalComponent, &
                    vorticityGradientTangentialComponent)

      endif

      !$acc exit data delete(normalizedPlanetaryVorticityVertex, &
      !$acc                  normalizedRelativeVorticityVertex)
      deallocate(normalizedPlanetaryVorticityVertex, &
                 normalizedRelativeVorticityVertex)

   !--------------------------------------------------------------------

   end subroutine ocn_diagnostic_solve_vorticity!}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve_richardson
!
!> \brief   Computes diagnostic variables for gradient richardson
!> \author  Matt Turner
!> \date    October 2020
!> \details
!>  This routine computes the diagnostic variables for gradient richardson
!>   number and brunt viasala
!
!-----------------------------------------------------------------------
   subroutine ocn_diagnostic_solve_richardson(displacedDensity, density, zMid, &
                normalVelocity, edgeAreaFractionOfCell, BruntVaisalaFreqTop, RiTopOfCell)!{{{

      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        displacedDensity
      real (kind=RKIND), dimension(:,:), intent(in) :: &
        density
      real (kind=RKIND), dimension(:,:), intent(in) :: &
        zMid
      real (kind=RKIND), dimension(:,:), intent(in) :: &
        normalVelocity
      real (kind=RKIND), dimension(:,:), intent(in) :: &
        edgeAreaFractionOfCell

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
        BruntVaisalaFreqTop
      real (kind=RKIND), dimension(:,:), intent(out) :: &
        RiTopOfCell

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real(kind=RKIND) :: coef, shear, Ri_edge, dz_edge, shearSquared, delU2, shearMean
      integer :: jCell, iCell, k, i, iEdge
      integer :: nCells
      real(kind=RKIND), dimension(:), allocatable :: RiTopOfCellNorm

      allocate (RiTopOfCellNorm(nVertLevels))
      !
      ! Brunt-Vaisala frequency (this has units of s^{-2})
      !
      coef = -gravity / rho_sw
#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc          present(BruntVaisalaFreqTop, displacedDensity, density, zMid, maxLevelCell, &
      !$acc                  minLevelCell) &
      !$acc          private(k)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iCell = 1, nCellsAll
         BruntVaisalaFreqTop(minLevelCell(iCell),iCell) = 0.0_RKIND
#ifdef MPAS_OPENACC
         !$acc loop seq
#endif
         do k = minLevelCell(iCell)+1, maxLevelCell(iCell)
            BruntVaisalaFreqTop(k,iCell) = coef * (displacedDensity(k-1,iCell) - density(k,iCell)) &
              / (zMid(k-1,iCell) - zMid(k,iCell))
          end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !
      ! Gradient Richardson number
      !

      nCells = nCellsHalo( 2 )
#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc          present(RiTopOfCell, maxLevelCell, nEdgesOnCell, edgesOnCell, normalVelocity, &
      !$acc                  edgeAreaFractionOfCell, zMid, BruntVaisalaFreqTop, minLevelCell) &
      !$acc          private(k, i, iEdge, delU2, shearMean, jCell, shear, dz_edge, Ri_edge, &
      !$acc                  RiTopOfCellNorm)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k, i, iEdge, delU2, shearMean, jCell, shear, &
      !$omp                              dz_edge, Ri_edge, RiTopOfCellNorm)
#endif
      do iCell=1,nCells
         RiTopOfCell(:,iCell) = 100.0_RKIND
         RiTopOfCellNorm(:) = 1.0E-12_RKIND
           do i = 1, nEdgesOnCell(iCell)
             iEdge = edgesOnCell(i, iCell)
            jCell = cellsOnCell(i,iCell)
            do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)
               delU2 = (normalVelocity(k-1,iEdge) - normalVelocity(k,iEdge))**2 + &
                  ( tangentialVelocity(k-1,iEdge) - tangentialVelocity(k,iEdge))**2

               dz_edge = 0.5_RKIND * ( &
                  zMid(k-1,iCell) + zMid(k-1,jCell) - &
                  zMid(k-0,iCell) - zMid(k-0,jCell) )

               shear = delU2 / dz_edge**2

               Ri_edge = max(0.0_RKIND, 0.5_RKIND * ( &
                  BruntVaisalaFreqTop(k,iCell) + &
                  BruntVaisalaFreqTop(k,jCell) ) )/ (shear + 1.0E-12_RKIND)

               RiTopOfCellNorm(k) = RiTopOfCellNorm(k) &
                                    + 0.25_RKIND * dcEdge(iEdge) * dvEdge(iEdge)

               RiTopOfCell(k,iCell) = &
               RiTopOfCell(k,iCell) + 0.25_RKIND * dcEdge(iEdge) * dvEdge(iEdge) * Ri_edge
           enddo
           enddo
         do k = minLevelCell(iCell), maxLevelCell(iCell)
            RiTopOfCell(k,iCell) = RiTopOfCell(k,iCell) / RiTopOfCellNorm(k)
          end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

   deallocate(RiTopOfCellNorm)
   end subroutine ocn_diagnostic_solve_richardson!}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve_surfaceLayer
!
!> \brief   Computes diagnostic variables for surface layer
!> \author  Matt Turner
!> \date    October 2020
!> \details
!>  This routine computes the diagnostic variables for averages of tracer
!>    values at surface layer
!
!-----------------------------------------------------------------------

   subroutine ocn_diagnostic_solve_surfaceLayer(layerThickness, &
                activeTracers, layerThicknessEdgeFlux, normalVelocity, &
                tracersSurfaceLayerValue, indexSurfaceLayerDepth, &
                normalVelocitySurfaceLayer, &
                surfaceFluxAttenuationCoefficient, &
                surfaceFluxAttenuationCoefficientRunoff)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness,         &!< [in] layer thickness
         layerThicknessEdgeFlux, &!< [in] layer thick at edge (flx form)
         normalVelocity           !< [in] velocity

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         activeTracers            !< [in] tracers to be averaged

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         tracersSurfaceLayerValue     !< [out] tracer avg over sfc lyr

      real (kind=RKIND), dimension(:), intent(out) :: &
         indexSurfaceLayerDepth,     &!< [out] layer index at sfc lyr depth
         normalVelocitySurfaceLayer, &!< [out] avg velocity over sfc layer
         surfaceFluxAttenuationCoefficient, &!< [out] sfc flx attenuation
         surfaceFluxAttenuationCoefficientRunoff !< [out] sfc flx attenuation

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         nCells, nEdges,  &! num of cells, edges
         n, nTracers,     &! num of tracers, tracer index
         iCell, iEdge, k, &! loop indices for cell, edge vertical
         kmin, kmax, ksfc,&! vertical indices for active, sfc layers
         cell1, cell2      ! neighbor cell indices across edge

      real(kind=RKIND) :: &
         sfcLayerDepth,   &! surface layer depth
         sumSfcLayer,     &! running sum of depth for sfc lyr calc
         fracSfcLayer,    &! fraction of layer in which sfc layer exists
         rSfcLayer         ! real vertical layer index for sfc lyr

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      ! average tracer values over the ocean surface layer
      ! the ocean surface layer is generally assumed to be about 0.1 of
      !    the boundary layer depth

      ! sfc layer averages currently only used when KPP on
      if (config_use_cvmix_kpp) then

         sfcLayerDepth = config_cvmix_kpp_surface_layer_averaging
         nTracers = size(activeTracers, dim=1)
#ifdef MPAS_OPENACC
         !$acc parallel loop gang vector &
         !$acc    present(tracersSurfaceLayerValue, activeTracers, &
         !$acc            layerThickness, minLevelCell, maxLevelCell, &
         !$acc            indexSurfaceLayerDepth) &
         !$acc    private(k, kmin, kmax, ksfc, &
         !$acc            rSfcLayer, sumSfcLayer, fracSfcLayer)
#else
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp    private(k, kmin, kmax, ksfc, &
         !$omp            rSfcLayer, sumSfcLayer, fracSfcLayer)
#endif
         do iCell=1,nCellsOwned
            kmin = minLevelCell(iCell)
            kmax = maxLevelCell(iCell)
            sumSfcLayer = 0.0_RKIND
            rSfcLayer = kmax
#ifdef MPAS_OPENACC
            !$acc loop seq
#endif
            do k = kmin,kmax
               sumSfcLayer = sumSfcLayer + layerThickness(k,iCell)
               if (sumSfcLayer > sfcLayerDepth) then
                  sumSfcLayer = sumSfcLayer - layerThickness(k,iCell)
                  rSfcLayer = (k-1) + (sfcLayerDepth - sumSfcLayer)/ &
                                       layerThickness(k,iCell)
                  exit
               endif
            end do ! k
            indexSurfaceLayerDepth(iCell) = rSfcLayer
            ksfc = int(rSfcLayer)
            fracSfcLayer = rSfcLayer - real(ksfc,RKIND)

            tracersSurfaceLayerValue(:,iCell) = 0.0_RKIND
#ifdef MPAS_OPENACC
            !$acc loop seq
#endif
            do k=kmin,ksfc
            do n=1,nTracers
               tracersSurfaceLayerValue(n,iCell) = &
               tracersSurfaceLayerValue(n,iCell) + &
                        activeTracers(n,k,iCell)*layerThickness(k,iCell)
            enddo
            enddo

            k = min(ksfc+1, kmax)
            do n=1,nTracers
               tracersSurfaceLayerValue(n,iCell) = &
              (tracersSurfaceLayerValue(n,iCell) + &
                         fracSfcLayer*activeTracers(n,k,iCell)* &
                         layerThickness(k,iCell)) / sfcLayerDepth
            enddo
         enddo ! cell loop
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

         ! average normal velocity values over the ocean surface layer
         ! the ocean surface layer is generally assumed to be about 0.1
         !    of the boundary layer depth

#ifdef MPAS_OPENACC
         !$acc parallel loop gang vector &
         !$acc    present(normalVelocitySurfaceLayer, normalVelocity, &
         !$acc            minLevelEdgeBot, maxLevelEdgeTop, &
         !$acc            cellsOnEdge, layerThicknessEdgeFlux) &
         !$acc    private(k, kmin, kmax, ksfc, cell1, cell2, &
         !$acc            rSfcLayer, sumSfcLayer, fracSfcLayer)
#else
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp    private(k, kmin, kmax, ksfc, cell1, cell2, &
         !$omp            rSfcLayer, sumSfcLayer, fracSfcLayer)
#endif
         do iEdge=1,nEdgesOwned
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            kmin  = minLevelEdgeBot(iEdge)
            kmax  = maxLevelEdgeTop(iEdge)
            sumSfcLayer = 0.0_RKIND
            rSfcLayer   = min(kmin,kmax)

#ifdef MPAS_OPENACC
            !$acc loop seq
#endif
            do k = kmin,kmax
               sumSfcLayer = sumSfcLayer + layerThicknessEdgeFlux(k,iEdge)
               if (sumSfcLayer > sfcLayerDepth) then
                  sumSfcLayer = sumSfcLayer &
                              - layerThicknessEdgeFlux(k,iEdge)
                  rSfcLayer = (k-1) + (sfcLayerDepth - sumSfcLayer)/ &
                                       layerThicknessEdgeFlux(k,iEdge)
                  exit
               endif
            end do ! k
            ksfc = int(rSfcLayer)
            fracSfcLayer = rSfcLayer - real(ksfc,RKIND)

            normalVelocitySurfaceLayer(iEdge) = 0.0_RKIND
#ifdef MPAS_OPENACC
            !$acc loop seq
#endif
            do k = kmin,ksfc
               normalVelocitySurfaceLayer(iEdge) = &
               normalVelocitySurfaceLayer(iEdge) + &
                                       normalVelocity(k,iEdge)* &
                                       layerThicknessEdgeFlux(k,iEdge)
            end do
            k = ksfc+1
            if (k <= kmax) then
               normalVelocitySurfaceLayer(iEdge) = &
              (normalVelocitySurfaceLayer(iEdge) + &
                        fracSfcLayer*normalVelocity(k,iEdge)* &
                        layerThicknessEdgeFlux(k,iEdge)) / sfcLayerDepth
            end if
         enddo ! edge loop
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      endif ! KPP on

      ! compute the attenuation coefficient for surface fluxes

      nCells = nCellsHalo( 1 )

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(surfaceFluxAttenuationCoefficient, &
      !$acc            surfaceFluxAttenuationCoefficientRunoff)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCells
         surfaceFluxAttenuationCoefficient(iCell) = &
                      config_flux_attenuation_coefficient
         surfaceFluxAttenuationCoefficientRunoff(iCell) = &
                      config_flux_attenuation_coefficient_runoff
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

   !--------------------------------------------------------------------

   end subroutine ocn_diagnostic_solve_surfaceLayer !}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve_surfaceLayer_subglacialRunoff
!
!> \brief   Computes diagnostic subglacialRunoff variables for surface layer
!> \author  Irena Vankova
!> \date    July 2024
!> \details
!
!-----------------------------------------------------------------------

   subroutine ocn_diagnostic_solve_surfaceLayer_subglacialRunoff(surfaceFluxAttenuationCoefficientSubglacialRunoff)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(out) :: &
         surfaceFluxAttenuationCoefficientSubglacialRunoff !< [out] sfc flx attenuation

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         nCells,  &! num of cells
         iCell     ! loop indices for cell

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      ! compute the attenuation coefficient for subglacial runoff surface fluxes

      nCells = nCellsHalo( 1 )

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(surfaceFluxAttenuationCoefficientSubglacialRunoff)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCells
            surfaceFluxAttenuationCoefficientSubglacialRunoff(iCell) = &
                      config_flux_attenuation_coefficient_subglacial_runoff
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

   !--------------------------------------------------------------------

   end subroutine ocn_diagnostic_solve_surfaceLayer_subglacialRunoff !}}}

!***********************************************************************


!***********************************************************************
!
!  routine ocn_diagnostic_solve_GMvel
!
!> \brief   Computes vertical velocity from GM bolus velocity
!> \author  Luke Van Roekel
!> \date    August 2022
!> \details
!>  This routine computes the diagnostic vertical velocity from
!>    the normal velocity due to GM bolus velocity
!
!-----------------------------------------------------------------------

   subroutine ocn_diagnostic_solve_GMvel(layerThickEdgeFlux, &
                                         normalGMBolusVelocity, &
                                         vertGMBolusVelocityTop)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickEdgeFlux,     &!< [in] thickness flux at edge
         normalGMBolusVelocity    !< [in] Bolus velocity normal to edge

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         vertGMBolusVelocityTop !< [out] Bolus vertical velocity at top

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         i, k,         &! loop indices for neighbors, vertical levels
         iCell, iEdge, &! addresses for cells, edges
         eoe,          &! address for edges on edge
         kmin, kmax,   &! min, max indices of active vert layers
         edgeSignOnCell_temp ! directional sign across edge

      real(kind=RKIND) :: &
         invAreaCell1, &! 1/cell area
         dvEdge_temp    ! local value of dvEdge

      real (kind=RKIND), dimension(:), allocatable :: &
         div_huGMBolus ! divergence of thick-weighted Bolus velocity

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      ! Return immediately if not needed
      if (.not. config_use_gm) return

      allocate(div_huGMBolus(nVertLevels))
      !$acc enter data create(div_huGMBolus)

      ! Compute vertical velocity from the divergence of horizontal GM
      ! Bolus velocity

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc    present(vertGMBolusVelocityTop, &
      !$acc            normalGMBolusVelocity, layerThickEdgeFlux, &
      !$acc            dvEdge, invAreaCell, edgeSignOnCell, &
      !$acc            nEdgesOnCell, edgesOnCell, &
      !$acc            minLevelCell, maxLevelCell) &
      !$acc    private(div_huGMBolus, i, k, kmin, kmax, iEdge, &
      !$acc            invAreaCell1, edgeSignOnCell_temp, dvEdge_temp)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(div_huGMBolus, i, k, kmin, kmax, iEdge, &
      !$omp            invAreaCell1, edgeSignOnCell_temp, dvEdge_temp)
#endif
      do iCell = 1, nCellsOwned
         div_huGMBolus(:) = 0.0_RKIND
         invAreaCell1 = invAreaCell(iCell)
         kmin = minLevelCell(iCell)
         kmax = maxLevelCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            edgeSignOnCell_temp = edgeSignOnCell(i, iCell)
            dvEdge_temp = dvEdge(iEdge)
            do k = kmin, kmax
               ! Compute vertical velocity from the horizontal GM Bolus velocity
               div_huGMBolus(k) = div_huGMBolus(k) &
                                - layerThickEdgeFlux(k,iEdge)* &
                                  edgeSignOnCell_temp*dvEdge_temp* &
                                  normalGMBolusVelocity(k,iEdge)* &
                                  invAreaCell1
            end do
         end do
         ! Vertical velocity at bottom is zero, initialized above.
         vertGMBolusVelocityTop(1:kmin-1,iCell) = 0.0_RKIND
         vertGMBolusVelocityTop(kmax+1  ,iCell) = 0.0_RKIND
         do k = kmax, 1, -1
            vertGMBolusVelocityTop(k,iCell) = &
            vertGMBolusVelocityTop(k+1,iCell) - div_huGMBolus(k)
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !$acc exit data delete(div_huGMBolus)
      deallocate(div_huGMBolus)

   !--------------------------------------------------------------------

   end subroutine ocn_diagnostic_solve_GMvel!}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve_MLEvel
!
!> \brief   Computes vertical velocity from MLE bolus velocity
!> \author  Luke Van Roekel
!> \date    August 2022
!> \details
!>  This routine computes the diagnostic vertical velocity from
!>    the normal velocity due to MLE (submesoscale eddy) bolus velocity
!
!-----------------------------------------------------------------------

   subroutine ocn_diagnostic_solve_MLEvel(layerThickEdgeFlux, &
                                          normalMLEVelocity, &
                                          vertMLEBolusVelocityTop)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickEdgeFlux,  &!< [in] thickness flux at edge
         normalMLEVelocity     !< [in] MLE velocity normal to edge

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         vertMLEBolusVelocityTop !< [out] MLE velocity at top

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         i, k,         &! loop indices for neighbors, vertical levels
         iCell, iEdge, &! addresses for cells, edges
         eoe,          &! address for edges on edge
         kmin, kmax,   &! min, max indices of active vert layers
         edgeSignOnCell_temp ! directional sign across edge

      real(kind=RKIND) :: &
         invAreaCell1, &! 1/cell area
         dvEdge_temp    ! local value of dvEdge

      real (kind=RKIND), dimension(:), allocatable :: &
         div_huMLEBolus ! divergence of thick-weighted MLE velocity

      !-----------------------------------------------------------------

      ! Return if not needed
      if (.not. config_submesoscale_enable) return

      allocate(div_huMLEBolus(nVertLevels))
      !$acc enter data create(div_huMLEBolus)

      ! Compute vertical velocity from the divergence of horizontal MLE
      ! (submesoscale eddy)  Bolus velocity

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc    present(vertMLEBolusVelocityTop, &
      !$acc            normalMLEVelocity, layerThickEdgeFlux, &
      !$acc            dvEdge, invAreaCell, edgeSignOnCell, &
      !$acc            nEdgesOnCell, edgesOnCell, &
      !$acc            minLevelCell, maxLevelCell) &
      !$acc    private(div_huMLEBolus, i, k, kmin, kmax, iEdge, &
      !$acc            invAreaCell1, edgeSignOnCell_temp, dvEdge_temp)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(div_huMLEBolus, i, k, kmin, kmax, iEdge, &
      !$omp            invAreaCell1, edgeSignOnCell_temp, dvEdge_temp)
#endif
      do iCell = 1, nCellsOwned
         div_huMLEBolus(:) = 0.0_RKIND
         invAreaCell1 = invAreaCell(iCell)
         kmin = minLevelCell(iCell)
         kmax = maxLevelCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            edgeSignOnCell_temp = edgeSignOnCell(i, iCell)
            dvEdge_temp = dvEdge(iEdge)
            do k = kmin,kmax
               div_huMLEBolus(k) = div_huMLEBolus(k) &
                                 - layerThickEdgeFlux(k,iEdge)* &
                                   edgeSignOnCell_temp*dvEdge_temp* &
                                   normalMLEVelocity(k,iEdge)* &
                                   invAreaCell1
            end do
         end do
         ! Vertical velocity at bottom is zero, initialized above.
         vertMLEBolusVelocityTop(1:kmin-1, iCell) = 0.0_RKIND
         vertMLEBolusVelocityTop(kmax+1,   iCell) = 0.0_RKIND
         do k = kmax, 1, -1
            vertMLEBolusVelocityTop(k,iCell) = &
            vertMLEBolusVelocityTop(k+1,iCell) - div_huMLEBolus(k)
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !$acc exit data delete(div_huMLEBolus)
      deallocate(div_huMLEBolus)

   !--------------------------------------------------------------------

   end subroutine ocn_diagnostic_solve_MLEvel!}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve_vortVel
!
!> \brief   Computes several velocity-related fields
!> \author  Matt Turner
!> \date    October 2020
!> \details
!>  This routine a number of fields related to horizontal velocity,
!>  namely kinetic energy, divergence, vorticity, vertical velocity
!>  and tangential velocity
!
!-----------------------------------------------------------------------

   subroutine ocn_diagnostic_solve_vortVel( &
                           relativeVorticity, layerThicknessEdgeFlux, &
                           normalVelocity, normalTransportVelocity, &
                           vertVelocityTop, vertTransportVelocityTop, &
                           relativeVorticityCell, divergence, &
                           kineticEnergyCell, tangentialVelocity)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         relativeVorticity,      &!< [in] relative vorticity in cells
         layerThicknessEdgeFlux, &!< [in] thickness flux at edge
         normalVelocity,         &!< [in] velocity normal to edge
         normalTransportVelocity  !< [in] transport vel normal to edge

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         vertVelocityTop,       &!< [out] vertical velocity top of cell
         vertTransportVelocityTop, &!< [out] vert transport vel at top
         relativeVorticityCell, &!< [out] relative vorticity in cell
         divergence,            &!< [out] divergence
         kineticEnergyCell,     &!< [out] kinetic energy in cell
         tangentialVelocity      !< [out] tangential velocity at edge

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         i, j, k,               &! loop indices
         iCell, iVertex, iEdge, &! indices for cells, vertices, edges
         nCells, nEdges,        &! number of cells, edges
         kmin, kmax,            &! min, max active vert level index
         eoe,                   &! index for edge on edge
         edgeSignOnCell_temp     ! local tmp for sign across edge

      real(kind=RKIND) :: &
         invAreaCell1,      &! 1/cell area
         dcEdge_temp,       &! local val of dcEdge
         dvEdge_temp,       &! local val of dvEdge
         r_tmp,             &! tmp for common factors
         weightsOnEdge_temp  ! local tmp for weights on edge

      real (kind=RKIND), dimension(:), allocatable :: &
         div_hu,         &! local divergence of thick-weighted velocity
         div_huTransport  ! local div of thick-weighted transport vel

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      ! Compute relative vorticity in cells

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc    present(relativeVorticityCell, relativeVorticity, &
      !$acc            invAreaCell, kiteAreasOnVertex, &
      !$acc            nEdgesOnCell, kiteIndexOnCell, verticesOnCell, &
      !$acc            minLevelCell, maxLevelCell) &
      !$acc    private(i, iVertex, j, k, kmin, kmax, invAreaCell1)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(i, iVertex, j, k, kmin, kmax, invAreaCell1)
#endif
      do iCell = 1, nCellsOwned
         relativeVorticityCell(:,iCell) = 0.0_RKIND
         invAreaCell1 = invAreaCell(iCell)
         kmin = minLevelCell(iCell)
         kmax = maxLevelCell(iCell)

         do i = 1, nEdgesOnCell(iCell)
            j = kiteIndexOnCell(i, iCell)
            iVertex = verticesOnCell(i, iCell)
            do k = kmin,kmax
               relativeVorticityCell(k,iCell) = &
               relativeVorticityCell(k,iCell) + &
                        kiteAreasOnVertex(j,iVertex)* &
                        relativeVorticity(k,iVertex)*invAreaCell1
            end do
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      ! Need 0,1,2 halo cells
      nCells = nCellsHalo(2)

      !
      ! Compute divergence, kinetic energy, and vertical velocity
      !
      allocate(div_hu(nVertLevels),div_huTransport(nVertLevels))
      !$acc enter data create(div_hu, div_huTransport)

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc    present(divergence, kineticEnergyCell, &
      !$acc            vertVelocityTop, vertTransportVelocityTop, &
      !$acc            normalVelocity, normalTransportVelocity, &
      !$acc            layerThicknessEdgeFlux, dcEdge, dvEdge, &
      !$acc            invAreaCell, nEdgesOnCell, edgesOnCell, &
      !$acc            edgeSignOnCell, minLevelCell, maxLevelCell) &
      !$acc    private(div_hu, div_huTransport, i, k, kmin, kmax, &
      !$acc            iEdge, invAreaCell1, r_tmp, dcEdge_temp, &
      !$acc            dvEdge_temp, edgeSignOnCell_temp)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(div_hu, div_huTransport, i, k, kmin, kmax, &
      !$omp            iEdge, invAreaCell1, r_tmp, dcEdge_temp, &
      !$omp            dvEdge_temp, edgeSignOnCell_temp)
#endif
      do iCell = 1, nCells
         divergence(:,iCell) = 0.0_RKIND
         kineticEnergyCell(:,iCell) = 0.0_RKIND
         div_hu(:) = 0.0_RKIND
         div_huTransport(:) = 0.0_RKIND
         invAreaCell1 = invAreaCell(iCell)
         kmin = minLevelCell(iCell)
         kmax = maxLevelCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            edgeSignOnCell_temp = edgeSignOnCell(i, iCell)
            dcEdge_temp = dcEdge(iEdge)
            dvEdge_temp = dvEdge(iEdge)
            do k = kmin,kmax
               r_tmp = dvEdge_temp*normalVelocity(k,iEdge)*invAreaCell1

               divergence(k,iCell) = divergence(k,iCell) &
                                   - edgeSignOnCell_temp*r_tmp
               div_hu(k) = div_hu(k) &
                         - layerThicknessEdgeFlux(k,iEdge)* &
                           edgeSignOnCell_temp*r_tmp
               div_huTransport(k) = div_huTransport(k) &
                                  - layerThicknessEdgeFlux(k,iEdge)* &
                                    edgeSignOnCell_temp*dvEdge_temp* &
                                    normalTransportVelocity(k,iEdge)* &
                                    invAreaCell1
               kineticEnergyCell(k,iCell) = kineticEnergyCell(k,iCell) &
                                          + 0.25*r_tmp*dcEdge_temp* &
                                            normalVelocity(k,iEdge)
            end do
         end do
         ! Vertical velocity at bottom is zero, initialized above.
         vertVelocityTop(1:kmin-1,iCell) = 0.0_RKIND
         vertVelocityTop(kmax+1  ,iCell) = 0.0_RKIND
         vertTransportVelocityTop(1:kmin-1,iCell) = 0.0_RKIND
         vertTransportVelocityTop(kmax+1  ,iCell) = 0.0_RKIND
         do k = kmax, 1, -1
            vertVelocityTop(k,iCell) = &
            vertVelocityTop(k+1,iCell) - div_hu(k)
            vertTransportVelocityTop(k,iCell) = &
            vertTransportVelocityTop(k+1,iCell) - div_huTransport(k)
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !$acc exit data delete(div_hu, div_huTransport)
      deallocate(div_hu,div_huTransport)

      ! Compute tangential velocity

      nEdges = nEdgesHalo(2)

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc    present(tangentialVelocity, normalVelocity, &
      !$acc            weightsOnEdge, nEdgesOnEdge, edgesOnEdge, &
      !$acc            minLevelEdgeBot, maxLevelEdgeTop) &
      !$acc    private(i, k, kmin, kmax, eoe, weightsOnEdge_temp)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(i, k, kmin, kmax, eoe, weightsOnEdge_temp)
#endif
      do iEdge = 1, nEdges
         tangentialVelocity(:,iEdge) = 0.0_RKIND
         kmin = minLevelEdgeBot(iEdge)
         kmax = maxLevelEdgeTop(iEdge)
         do i = 1, nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            weightsOnEdge_temp = weightsOnEdge(i, iEdge)
            do k = kmin,kmax
               tangentialVelocity(k,iEdge) = &
               tangentialVelocity(k,iEdge) + weightsOnEdge_temp* &
                                             normalVelocity(k, eoe)
            end do
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

   !--------------------------------------------------------------------

   end subroutine ocn_diagnostic_solve_vortVel!}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve_vertVel_remap
!
!> \brief   Computes diagnostic variables for velocity
!> \author  Carolyn Begeman
!> \date    October 2021
!> \details
!>  This routine computes the Eulerian vertical velocity valid for
!>    Vertical Lagrangian Remapping
!
!-----------------------------------------------------------------------
   subroutine ocn_diagnostic_solve_vertVel_remap(layerThickness, &
             layerThicknessLag, dt, vertVelocityTop)!{{{

      use ocn_mesh

      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness    ! layerThickness after remapping at timeLevelIn
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessLag ! layerThickness before remapping at timeLevelIn
      real (kind=RKIND), intent(in) :: dt

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         vertVelocityTop

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, k
      real (kind=RKIND), dimension(:), allocatable :: &
         heightCell, heightCellLag

      allocate(heightCell(nVertLevels+1))
      allocate(heightCellLag(nVertLevels+1))

      !$omp parallel
      !$omp do schedule(runtime) private(k, heightCell, heightCellLag)
      do iCell = 1, nCellsAll

         heightCell(:) = -1.0_RKIND * bottomDepth(iCell)
         heightCellLag(:) = -1.0_RKIND * bottomDepth(iCell)

         do k = maxLevelCell(iCell), minLevelCell(iCell), -1
            ! reconstruct heights
            heightCellLag(k) = heightCellLag(k + 1) + layerThicknessLag(k, iCell)
            heightCell(k) = heightCell(k + 1) + layerThickness(k, iCell)
         end do

         vertVelocityTop(:, iCell) = 0.0_RKIND

         do k = minLevelCell(iCell), maxLevelCell(iCell)
            ! compute vertical velocity at the top of each Cell
            vertVelocityTop(k,iCell) = (heightCellLag(k) - heightCell(k))/dt
         end do
      end do
      !$omp end do
      !$omp end parallel

      deallocate(heightCell)
      deallocate(heightCellLag)

   end subroutine ocn_diagnostic_solve_vertVel_remap!}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve_surface_pressure
!
!> \brief   Computes diagnostic variables for surface pressure
!> \author  Xylar Asay-Davis
!> \date    Jan 2021
!> \details
!>  This routine computes the diagnostic variable for surface pressure
!>  It must be called in the code after computing the density.
!
!-----------------------------------------------------------------------

   subroutine ocn_diagnostic_solve_surface_pressure(forcingPool, &
               atmosphericPressure, seaIcePressure, surfacePressure)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         forcingPool         !< [in] Forcing information

      real (kind=RKIND), dimension(:), intent(in) :: &
         atmosphericPressure, &!< [in] atmospheric pressure on surface
         seaIcePressure        !< [in] sea-ice pressure on surface

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(out) :: &
         surfacePressure       !< [out] total surface pressure

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: iCell ! cell loop index

      real (kind=RKIND), dimension(:), pointer :: &
         frazilSurfacePressure, &! frazil ice pressure from forcing
         landIcePressure         ! land ice pressure

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      ! Pressure
      !
      call mpas_pool_get_array(forcingPool, 'frazilSurfacePressure', &
                                             frazilSurfacePressure)
      call mpas_pool_get_array(forcingPool, 'landIcePressure', &
                                             landIcePressure)

      ! Compute pressure for generalized coordinates.
      ! Pressure at top surface may be due to atmospheric pressure,
      ! sea ice, frazil ice or the weight of an ice shelf
      ! The latter two are only available when those options are on.

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(surfacePressure, atmosphericPressure, &
      !$acc            seaIcePressure)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCellsAll

         surfacePressure(iCell) = atmosphericPressure(iCell) &
                                + seaIcePressure(iCell)
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

      ! Add frazil component if needed
      if ( associated(frazilSurfacePressure) ) then

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(surfacePressure, frazilSurfacePressure)
#else
         !$omp do schedule(runtime)
#endif
         do iCell = 1, nCellsAll
            surfacePressure(iCell) = surfacePressure(iCell) &
                                   + frazilSurfacePressure(iCell)
         end do
#ifndef MPAS_OPENACC
         !$omp end do
#endif
      endif ! frazil pressure

      ! Add land ice pressure if needed

      if (landIcePressureOn) then
#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(surfacePressure)
         !! need this eventually, but there is currently
         !! an issue, probably with the re-retrieval of pointer
         !!acc    present(surfacePressure, landIcePressure)
#else
         !$omp do schedule(runtime)
#endif
         do iCell = 1, nCellsAll
            surfacePressure(iCell) = surfacePressure(iCell) &
                                   + landIcePressure(iCell)
         end do
#ifndef MPAS_OPENACC
         !$omp end do
#endif
      endif ! land ice pressure
#ifndef MPAS_OPENACC
      !$omp end parallel
#endif

   !--------------------------------------------------------------------

   end subroutine ocn_diagnostic_solve_surface_pressure!}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve_z_coordinates
!
!> \brief   Computes diagnostic variables for mid- and top-layer z, and ssh
!> \author  Xylar Asay-Davis
!> \date    Jan. 2021
!> \details
!>  This routine computes the diagnostic variables for the z-coordinates at the
!>    middle and top of each layer, and the sea-surface height
!
!-----------------------------------------------------------------------
   subroutine ocn_diagnostic_solve_z_coordinates(layerThickness, zMid, zTop, ssh)!{{{

      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         zMid
      real (kind=RKIND), dimension(:,:), intent(out) :: &
         zTop
      real (kind=RKIND), dimension(:), intent(out) :: &
         ssh

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: nCells
      integer :: iCell
      integer :: k

      nCells = nCellsAll

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc          present(maxLevelCell, zMid, layerThickness, zTop, bottomDepth, ssh, &
      !$acc                  minLevelCell) &
      !$acc          private(k)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iCell = 1, nCells

         ! Compute zMid, the z-coordinate of the middle of the layer.
         ! Compute zTop, the z-coordinate of the top of the layer.
         ! Note the negative sign, since bottomDepth is positive
         ! and z-coordinates are negative below the surface.
         k = maxLevelCell(iCell)
         zMid(k:nVertLevels,iCell) = -bottomDepth(iCell) + 0.5_RKIND*layerThickness(k,iCell)
         zTop(k:nVertLevels,iCell) = -bottomDepth(iCell) + layerThickness(k,iCell)

         do k = maxLevelCell(iCell)-1, minLevelCell(iCell), -1
            zMid(k,iCell) = zMid(k+1,iCell)  &
              + 0.5_RKIND*(  layerThickness(k+1,iCell) &
                           + layerThickness(k  ,iCell))
            zTop(k,iCell) = zTop(k+1,iCell)  &
              + layerThickness(k  ,iCell)
         end do

         ! copy zTop(1,iCell) into sea-surface height array
         ssh(iCell) = zTop(minLevelCell(iCell),iCell)

      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      if (config_use_subgrid_wetting_drying) then
#ifdef MPAS_OPENACC
      !$acc update host(layerThickness)
#endif
         do iCell = 1, nCells
               call ocn_subgrid_ssh_lookup(layerThickness(1,iCell), &
                                           subgridWetVolumeCellTable(:,iCell), &
                                           subgridSshCellTableRange(:,iCell), &
                                           bottomDepth(iCell), &
                                           subgridCellBathymetryMin(iCell), & 
                                           ssh(iCell))

               zTop(1,iCell) = ssh(iCell)
               zMid(1,iCell) = -bottomDepth(iCell) + 0.5_RKIND*layerThickness(1,iCell)         
         end do
#ifdef MPAS_OPENACC
      !$acc update device(zTop, zMid, ssh)
#endif
      end if

   end subroutine ocn_diagnostic_solve_z_coordinates!}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve_pressure
!
!> \brief   Computes pressure or Montgomery potential at all levels
!> \author  Matt Turner, Xylar Asay-Davis
!> \date    Jan. 2021
!> \details
!>  This routine computes the pressure or Montgomery potential
!
!-----------------------------------------------------------------------

   subroutine ocn_diagnostic_solve_pressure(layerThickness, density, &
                surfacePressure, montgomeryPotential, pressure)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness, &!< [in] layer thickness
         density          !< [in] density

      real (kind=RKIND), dimension(:), intent(in) :: &
         surfacePressure  !< [in] pressure at surface

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         montgomeryPotential, &!< [out] Montgomery potential at all lvls
         pressure              !< [out] pressure at all levels

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         nCells,    &! number of cells
         iCell, k,  &! cell, vertical loop indices
         kmin, kmax  ! min, max active vertical levels

      real (kind=RKIND), dimension(:), allocatable :: &
         pTop      ! local pressure at top of layer

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !
      ! If needed, compute Montgomery potential
      !   This is the preferred option when layers are isopycnal though
      !     pressure_and_zmid option can be used for isopycnal as well
      !

      nCells = nCellsAll
      if (config_pressure_gradient_type == 'MontgomeryPotential') then

         allocate(pTop(nVertLevels))
         !$acc enter data create(pTop)

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(montgomeryPotential, density, bottomDepth, &
         !$acc            layerThickness, minLevelCell) &
         !$acc    private(pTop, k, kmin)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(pTop, k, kmin)
#endif
         do iCell = 1, nCells

            ! assume atmospheric pressure at surface is zero for now
            pTop(1) = 0.0_RKIND
            ! At top layer it is g*SSH, where SSH may be off by a
            ! constant (ie, bottomDepth can be relative to top or bottom)
            kmin = minLevelCell(iCell)
            montgomeryPotential(kmin,iCell) = gravity &
                   * (bottomDepth(iCell) + &
                      sum(layerThickness(1:nVertLevels,iCell)))

            do k = 2, nVertLevels
               pTop(k) = pTop(k-1) + density(k-1,iCell)*gravity* &
                              layerThickness(k-1,iCell)

               ! from delta M = p delta / density
               montgomeryPotential(k  ,iCell) = &
               montgomeryPotential(k-1,iCell) + pTop(k)* &
                             (1.0_RKIND/density(k  ,iCell) -  &
                              1.0_RKIND/density(k-1,iCell))
            end do

         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

         !$acc exit data delete(pTop)
         deallocate(pTop)

      !
      ! Otherwise, compute pressure for generalized coordinates.
      !

      else

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(pressure, surfacePressure, density, &
         !$acc            layerThickness, minLevelCell, maxLevelCell) &
         !$acc    private(k, kmin, kmax)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k, kmin, kmax)
#endif
         do iCell = 1, nCells
            kmin = minLevelCell(iCell)
            kmax = maxLevelCell(iCell)
            pressure(kmin,iCell) = surfacePressure(iCell) &
                                 + density(kmin,iCell)*gravity* &
                                   0.5_RKIND*layerThickness(kmin,iCell)

            do k = kmin+1, kmax
               pressure(k,iCell) = pressure(k-1,iCell) &
                                 + 0.5_RKIND*gravity* &
                    (  density(k-1,iCell)*layerThickness(k-1,iCell) &
                     + density(k  ,iCell)*layerThickness(k  ,iCell))
            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      endif

   !--------------------------------------------------------------------

   end subroutine ocn_diagnostic_solve_pressure!}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve_ssh
!
!> \brief   Computes diagnostic variables for pressure adjusted ssh and gradSSH
!> \author  Matt Turner
!> \date    October 2020
!> \details
!>  This routine computes the diagnostic variables for pressureAdjustedSSH and gradSSH
!
!-----------------------------------------------------------------------
   subroutine ocn_diagnostic_solve_ssh(forcingPool, ssh, seaIcePressure, &
                pressureAdjustedSSH, gradSSH, landIceDraft)!{{{

      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: forcingPool !< Input: Forcing information
      real (kind=RKIND), dimension(:), intent(in) :: &
         ssh
      real (kind=RKIND), dimension(:), intent(in) :: &
         seaIcePressure
      real (kind=RKIND), dimension(:), optional, intent(in) :: &
         landIceDraft

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(out) :: &
         pressureAdjustedSSH
      real (kind=RKIND), dimension(:), intent(out) :: &
         gradSSH

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: nCells, nEdges
      integer :: iCell, iEdge
      integer :: cell1, cell2

      nCells = nCellsHalo( 1 )
#ifdef MPAS_OPENACC
      !$acc parallel loop present(seaIcePressure, pressureAdjustedSSH, ssh)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCells
         pressureAdjustedSSH(iCell) = ssh(iCell) + ( seaIcePressure(iCell) / ( gravity * rho_sw ) )
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      if (landIcePressureOn) then
#ifdef MPAS_OPENACC
         !! This first directive should be correct, but there may
         !! be a problem with the optional argument. Removing optional
         !! landIceDraft from the present list for now.
         !!acc parallel loop present(pressureAdjustedSSH, landIceDraft)
         !$acc parallel loop present(pressureAdjustedSSH)
#else
         !$omp parallel
         !$omp do schedule(runtime)
#endif
         do iCell = 1, nCells
            ! subtract the land ice draft from the SSH so sea ice doesn't experience tilt
            ! toward land ice
            pressureAdjustedSSH(iCell) = pressureAdjustedSSH(iCell) - landIceDraft(iCell)
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
      end if

      nEdges = nEdgesHalo( 1 )

#ifdef MPAS_OPENACC
      !$acc parallel loop present(pressureAdjustedSSH, dcEdge, cellsOnEdge, edgeMask, gradSSH)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(cell1, cell2)
#endif
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)

         gradSSH(iEdge) = edgeMask(minLevelEdgeBot(iEdge), iEdge) * ( pressureAdjustedSSH(cell2) - pressureAdjustedSSH(cell1) ) / dcEdge(iEdge)
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

   end subroutine ocn_diagnostic_solve_ssh!}}}

!***********************************************************************
!
!  routine ocn_vert_transport_velocity_top
!
!> \brief   Computes vertical transport
!> \author  Mark Petersen
!> \date    August 2013
!> \details
!>  This routine computes the vertical transport through the top of each
!>  cell.
!
!-----------------------------------------------------------------------

   subroutine ocn_vert_transport_velocity_top(verticalMeshPool, &
                     oldLayerThickness, layerThicknessEdgeFlux, &
                     normalVelocity, oldSSH, dt, vertAleTransportTop, &
                     err, newHighFreqThickness)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         verticalMeshPool         !< [in] vertical mesh information

      real (kind=RKIND), dimension(:), intent(in) :: &
         oldSSH                   !< [in] sea surface height at old time

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         oldLayerThickness,      &!< [in] layer thickness at old time
         layerThicknessEdgeFlux, &!< [in] layerThickness at edge
         normalVelocity           !< [in] transport

      !  alters ALE thickness when freq-filtering
      real (kind=RKIND), dimension(:,:), intent(in), optional :: &
         newHighFreqThickness     !< [in] high frequency thickness

      real (kind=RKIND), intent(in) :: &
         dt     !< Input: time step

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         vertAleTransportTop     !< [out] vertical transport top of cell

      integer, intent(out) :: err !< [out] error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         iEdge, iCell, k, i, &! edge, cell, vert, nbr indices
         nCells,             &! number of cells
         kmin, kmax           ! min/max layer indices for active layers

      real (kind=RKIND) :: &
         flux,         &!
         invAreaCell1, &! 1/cell area
         div_hu_btr     ! divergence of h*u for barotropic u

      real (kind=RKIND), dimension(:), allocatable :: &
         projectedSSH ! projectedSSH: projected SSH at a new time

      real (kind=RKIND), dimension(:,:), allocatable :: &
         div_hu,       &! divergence of (thickness*velocity)
         ALE_Thickness  ! ALE thickness at new time

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      err = 0

      if (config_vert_coord_movement == 'impermeable_interfaces' .or. &
          configVertAdvMethod == vertAdvRemap) then
         vertAleTransportTop = 0.0_RKIND
         !$acc update device(vertAleTransportTop)
         return
      end if

      allocate(div_hu(nVertLevels, nCellsAll), &
               projectedSSH(nCellsAll), &
               ALE_Thickness(nVertLevels, nCellsAll))
      !$acc enter data create(projectedSSH, ALE_Thickness, div_hu)

      ! Only need to compute over 0 and 1 halos
      nCells = nCellsHalo( 1 )

      !
      ! thickness-weighted divergence and barotropic divergence
      !
      ! See Ringler et al. (2010) jcp paper, eqn 19, 21, and fig. 3.
#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(projectedSSH, oldSSH, normalVelocity, div_hu, &
      !$acc            layerThicknessEdgeFlux, dvEdge, invAreaCell, &
      !$acc            nEdgesOnCell, edgesOnCell, edgeSignOnCell, &
      !$acc            minLevelEdgeBot, maxLevelEdgeTop) &
      !$acc    private(i, iEdge, k, kmin, kmax, div_hu_btr, &
      !$acc            invAreaCell1, flux)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(i, iEdge, k, kmin, kmax, div_hu_btr, &
      !$omp            invAreaCell1, flux)
#endif
      do iCell = 1, nCells
         div_hu(:,iCell) = 0.0_RKIND
         div_hu_btr      = 0.0_RKIND
         invAreaCell1 = invAreaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            kmin = minLevelEdgeBot(iEdge)
            kmax = maxLevelEdgeTop(iEdge)

            do k = kmin, kmax
               flux = layerThicknessEdgeFlux(k,iEdge)* &
                      normalVelocity(k,iEdge)*dvEdge(iEdge)* &
                      edgeSignOnCell(i,iCell) * invAreaCell1
               div_hu(k,iCell) = div_hu(k,iCell) - flux
               div_hu_btr = div_hu_btr - flux
            end do
         end do
         projectedSSH(iCell) = oldSSH(iCell) - dt*div_hu_btr
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !
      ! Compute desired thickness at new time
      !
      if (present(newHighFreqThickness)) then
        call ocn_ALE_thickness(verticalMeshPool, projectedSSH, &
                               ALE_thickness, err, newHighFreqThickness)
      else
        call ocn_ALE_thickness(verticalMeshPool, projectedSSH, &
                               ALE_thickness, err)
      endif

      !
      ! Vertical transport through layer interfaces
      !
      ! Vertical transport through layer interface at top and
      ! bottom is zero. Here we are using solving the continuity
      ! equation for vertAleTransportTop ($w^t$), and using
      ! ALE_Thickness for thickness at the new time.

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc          present(vertAleTransportTop, ALE_Thickness, &
      !$acc                  oldLayerThickness, div_hu, &
      !$acc                  minLevelCell, maxLevelCell)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iCell = 1,nCells
         vertAleTransportTop(1,iCell) = 0.0_RKIND
         vertAleTransportTop(maxLevelCell(iCell)+1,iCell) = 0.0_RKIND
         do k = maxLevelCell(iCell), minLevelCell(iCell)+1, -1
            vertAleTransportTop(k,iCell) = &
                vertAleTransportTop(k+1,iCell) - div_hu(k,iCell) &
                     - (ALE_Thickness(k,iCell) - &
                        oldLayerThickness(k,iCell))/dt
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !$acc exit data delete(projectedSSH, ALE_Thickness, div_hu)
      deallocate(div_hu, &
                 projectedSSH, &
                 ALE_Thickness)

   !--------------------------------------------------------------------

   end subroutine ocn_vert_transport_velocity_top!}}}

!***********************************************************************
!
!  routine ocn_fuperp
!
!> \brief   Computes f u_perp
!> \author  Mark Petersen
!> \date    23 September 2011
!> \details
!>  This routine computes f u_perp for the ocean
!
!-----------------------------------------------------------------------

   subroutine ocn_fuperp(statePool, meshPool, timeLevelIn)!{{{

      type (mpas_pool_type), intent(inout) :: statePool !< Input/Output: State information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      integer, intent(in), optional :: timeLevelIn !< Input: Input time level for state pool

      integer :: iEdge, cell1, cell2, eoe, j, k
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, normalBaroclinicVelocity

      integer :: timeLevel

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_timer_start("ocn_fuperp")

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocity, timeLevel)

      !DWJ: ADD OMP (Only needed for split explicit)

      !
      ! Put f*normalBaroclinicVelocity^{perp} in u as a work variable
      !
#ifdef MPAS_OPENACC
      !$acc enter data create(normalVelocity) copyin(normalBaroclinicVelocity)

      !$acc parallel loop gang vector &
      !$acc          present(cellsOnEdge, maxLevelEdgeTop, normalVelocity, nEdgesOnEdge, &
      !$acc                  edgesOnEdge, weightsOnEdge, fEdge, normalBaroclinicVelocity, &
      !$acc                  minLevelEdgeBot) &
      !$acc          private(cell1, cell2, k, j, eoe)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(cell1, cell2, k, eoe)
#endif
      do iEdge = 1, nEdgesOwned
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)

            normalVelocity(k,iEdge) = 0.0_RKIND
            do j = 1,nEdgesOnEdge(iEdge)
               eoe = edgesOnEdge(j,iEdge)
               normalVelocity(k,iEdge) = normalVelocity(k,iEdge) + weightsOnEdge(j,iEdge) * normalBaroclinicVelocity(k,eoe) &
                                       * fEdge(eoe)
            end do
         end do
      end do
#ifdef MPAS_OPENACC
      !$acc exit data copyout(normalVelocity) delete(normalBaroclinicVelocity)
#else
      !$omp end do
      !$omp end parallel
#endif

      call mpas_timer_stop("ocn_fuperp")

   end subroutine ocn_fuperp!}}}

!***********************************************************************
!
!  routine ocn_filter_btr_mode_tend_vel
!
!> \brief   ocn_filters barotropic mode out of the velocity tendency
!> \author  Mark Petersen
!> \date    23 September 2011
!> \details
!>  This routine filters barotropic mode out of the velocity tendency.
!
!-----------------------------------------------------------------------
   subroutine ocn_filter_btr_mode_tend_vel(tendPool, statePool, meshPool, timeLevelIn)!{{{

      type (mpas_pool_type), intent(inout) :: tendPool !< Input/Output: Tendency information
      type (mpas_pool_type), intent(in) :: statePool !< Input: State information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      integer, intent(in), optional :: timeLevelIn !< Input: Time level for state pool

      integer :: iEdge, k
      real (kind=RKIND) :: vertSum, normalThicknessFluxSum, thicknessSum
      real (kind=RKIND), dimension(:,:), pointer :: tend_normalVelocity

      integer :: timeLevel

      call mpas_timer_start("ocn_filter_btr_mode_tend_vel")

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_array(tendPool, 'normalVelocity', tend_normalVelocity)

      !$omp parallel
      !$omp do schedule(runtime) private(normalThicknessFluxSum, thicknessSum, k, vertSum)
      do iEdge = 1, nEdgesAll

        ! thicknessSum is initialized outside the loop because on land boundaries
        ! maxLevelEdgeTop=0, but I want to initialize thicknessSum with a
        ! nonzero value to avoid a NaN.
        normalThicknessFluxSum = layerThickEdgeFlux(maxLevelEdgeBot(iEdge),iEdge) * &
                                 tend_normalVelocity(maxLevelEdgeBot(iEdge),iEdge)
        thicknessSum  = layerThickEdgeFlux(maxLevelEdgeBot(iEdge),iEdge)

        do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)
          normalThicknessFluxSum = normalThicknessFluxSum + layerThickEdgeFlux(k,iEdge) * tend_normalVelocity(k,iEdge)
          thicknessSum  =  thicknessSum + layerThickEdgeFlux(k,iEdge)
        enddo

        vertSum = normalThicknessFluxSum / thicknessSum
        do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
          tend_normalVelocity(k,iEdge) = tend_normalVelocity(k,iEdge) - vertSum
        enddo
      enddo ! iEdge
      !$omp end do
      !$omp end parallel

      call mpas_timer_stop("ocn_filter_btr_mode_tend_vel")

   end subroutine ocn_filter_btr_mode_tend_vel!}}}

!***********************************************************************
!
!  routine ocn_compute_mixing_input_fields
!
!> \brief
!>    Compute fields necessary to drive the CVMix KPP and gotm modules
!> \author  Todd Ringler, Luke Van Roekel
!> \date    11 July 2024
!> \details
!>    CVMix/KPP requires the following fields as intent(in):
!>       surfaceBuoyancyForcing
!>       surfaceFrictionVelocity
!>
!
!-----------------------------------------------------------------------

    subroutine ocn_compute_mixing_input_fields(statePool, forcingPool, &
                                            meshPool, timeLevelIn)!{{{

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         statePool,  &!< [in] State information
         meshPool,   &!< [in] Mesh information (for pass to eq state)
         forcingPool  !< [in] Forcing information

      integer, intent(in), optional :: &
         timeLevelIn  !< [in] optional time level for state variables

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer :: &
         iCell, iEdge, i, k, &! loop indices for cell, edge, neighbors
         kmin,               &! topmost active cell index
         nCells,             &! number of cells
         err,                &! local error code
         timeLevel            ! time level for state variables (default 1)

      real (kind=RKIND) :: &
         fracAbsorbed,       &! fraction of sfc flux absorbed
         fracAbsorbedRunoff, &! same for runoff
         fracAbsorbedSubglacialRunoff, &! same for subglacial runoff
         zTop,zBot,          &! temporary variables
         transmissionCoeffTop,transmissionCoeffBot,     &! temporary variables
         sumSurfaceStressSquared ! sum of sfc stress squared

      ! pointers for variable/pool retrievals
      integer, pointer :: &
         indexTempFlux,   &! index of temperature in tracer array
         indexSaltFlux     ! index of salinity    in tracer array

      type (mpas_pool_type), pointer :: &
         tracersSurfaceFluxPool, &! pool of tracer surface fluxes
         tracersPool              ! pool of tracers

      real (kind=RKIND), dimension(:), pointer :: &
         penetrativeTemperatureFlux, &! various sfc flux components
         surfaceThicknessFlux, &
         surfaceThicknessFluxRunoff, &
         surfaceThicknessFluxSubglacialRunoff, &
         rainTemperatureFlux, &
         evapTemperatureFlux, &
         icebergTemperatureFlux, &
         seaIceTemperatureFlux, &
         sfcStressMag

      real (kind=RKIND), dimension(:,:), pointer ::  &
         layerThickness,           &! layer thickness
         normalVelocity,           &! normal velocity
         activeTracersSurfaceFlux, &! sfc flux of active tracers (T,S)
         activeTracersSurfaceFluxRunoff, &! flx of tracers in runoff
         activeTracersSurfaceFluxSubglacialRunoff, &! flx of tracers in subglacial runoff
         nonLocalSurfaceTracerFlux  ! non-local flux of tracers

      real (kind=RKIND), dimension(:,:,:), pointer :: &
         activeTracers  ! array of active tracers

      ! Scratch Arrays
      ! densitySurfaceDisplaced: Density computed by displacing SST
      !    and SSS to every vertical layer within the column
      !    units: kg m^{-3}
      ! thermalExpansionCoeff: Thermal expansion coefficient (alpha),
      !    defined as $-1/\rho d\rho/dT$ (note negative sign)
      !    units: C^{-1}
      ! salineContractionCoeff: Saline contraction coefficient (beta),
      !    defined as $1/\rho d\rho/dS$.  This is also called the
      !    haline contraction coefficient.  units: PSU^{-1}

      real (kind=RKIND), dimension(:,:), allocatable :: &
         densitySurfaceDisplaced, &
         thermalExpansionCoeff,   &
         salineContractionCoeff

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      call mpas_timer_start('Mixing input fields')

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      ! Retrieve variables from pool
      call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux', &
                                               tracersSurfaceFluxPool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      !TODO: should be hard-wired?
      call mpas_pool_get_array(tracersPool, 'activeTracers', &
                                             activeTracers, 1)

      call mpas_pool_get_dimension(tracersSurfaceFluxPool, &
                                  'index_temperatureSurfaceFlux', &
                                   indexTempFlux)
      call mpas_pool_get_dimension(tracersSurfaceFluxPool, &
                                  'index_salinitySurfaceFlux', &
                                   indexSaltFlux)

      call mpas_pool_get_array(statePool, 'normalVelocity', &
                                           normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'layerThickness', &
                                           layerThickness, timeLevel)

      call mpas_pool_get_array(tracersSurfaceFluxPool, &
                              'nonLocalSurfaceTracerFlux', &
                               nonLocalSurfaceTracerFlux)
      call mpas_pool_get_array(tracersSurfaceFluxPool, &
                              'activeTracersSurfaceFlux', &
                               activeTracersSurfaceFlux)
      call mpas_pool_get_array(tracersSurfaceFluxPool, &
                              'activeTracersSurfaceFluxRunoff', &
                               activeTracersSurfaceFluxRunoff)
      call mpas_pool_get_array(tracersSurfaceFluxPool, &
                              'activeTracersSurfaceFluxSubglacialRunoff', &
                               activeTracersSurfaceFluxSubglacialRunoff)

      call mpas_pool_get_array(forcingPool, 'surfaceThicknessFlux', &
                                             surfaceThicknessFlux)
      call mpas_pool_get_array(forcingPool, 'surfaceThicknessFluxRunoff', &
                                             surfaceThicknessFluxRunoff)
      call mpas_pool_get_array(forcingPool, 'surfaceThicknessFluxSubglacialRunoff', &
                                             surfaceThicknessFluxSubglacialRunoff)
      call mpas_pool_get_array(forcingPool, 'penetrativeTemperatureFlux', &
                                             penetrativeTemperatureFlux)
      call mpas_pool_get_array(forcingPool, 'surfaceStressMagnitude', &
                                             sfcStressMag)

      ! allocate scratch space displaced density computation
      ncells = nCellsAll

      allocate(densitySurfaceDisplaced(nVertLevels, nCells), &
               thermalExpansionCoeff(nVertLevels, nCells), &
               salineContractionCoeff(nVertLevels, nCells))
      !$acc enter data create(thermalExpansionCoeff,  &
      !$acc                   salineContractionCoeff, &
      !$acc                   densitySurfaceDisplaced)

      ! Only need to compute over the 0 and 1 halos
      nCells = nCellsHalo( 2 )

      ! compute displaced density defined as density computed
      ! from T, S that have been displaced (adiabatically) from surface

      call ocn_equation_of_state_density(statePool, meshPool, &
                          activeTracers, tracersSurfaceValue, &
                          nCells, 0, 'surfaceDisplaced', &
                          densitySurfaceDisplaced, err, &
                          thermalExpansionCoeff, &
                          salineContractionCoeff, timeLevel)

      !$acc exit data copyout(thermalExpansionCoeff,  &
      !$acc                   salineContractionCoeff, &
      !$acc                   densitySurfaceDisplaced)

      ! compute surface buoyancy forcing based on surface fluxes of
      ! mass, temperature, salinity and frazil (frazil to be added later)
      ! since this computation is confusing, variables, units and sign
      ! convention is repeated here
      ! everything below should be consistent with that specified in
      !   Registry and with the CVMix/KPP documentation:
      !    https://www.dropbox.com/s/6hqgc0rsoa828nf/cvmix_20aug2013.pdf
      !
      ! surfaceThicknessFlux: surface mass flux, m/s, positive into ocean
      ! activeTracersSurfaceFlux(indexTempFlux): non-penetrative
      !    temperature flux, C m/s, positive into ocean
      ! penetrativeTemperatureFlux: penetrative surface temperature
      !    flux at ocean surface, positive into ocean
      ! activeTracersSurfaceFlux(indexSaltFlux): salinity flux, PSU m/s,
      !    positive into ocean
      ! penetrativeTemperatureFluxOBL: penetrative temperature flux
      !    computed at z=OBL, positive down
      !
      ! note: the following fields used the CVMix/KPP computation of
      !    buoyancy forcing are not included here
      !    1. Tm: temperature associated with surfaceThicknessFlux,
      !           C  (here we assume Tm == temperatureSurfaceValue)
      !    2. Sm: salinity associated with surfaceThicknessFlux, PSU
      !           (here we assume Sm == salinitySurfaceValue and
      !            account for salinity flux in
      !            activeTracersSurfaceFlux array)
      !

      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(kmin, fracAbsorbed, fracAbsorbedRunoff, fracAbsorbedSubglacialRunoff, &
      !$omp         zTop, k, zBot, transmissionCoeffTop, transmissionCoeffBot, &
      !$omp         sumSurfaceStressSquared, i, iEdge)
      do iCell = 1, nCells

         kmin = minLevelCell(iCell)
         ! Compute fraction of thickness flux that is in the top model layer
         fracAbsorbed = 1.0_RKIND &
                      - exp( max(-100.0_RKIND, &
                                 -layerThickness(kmin, iCell)/ &
                           config_flux_attenuation_coefficient))
         fracAbsorbedRunoff = 1.0_RKIND &
                      - exp( max(-100.0_RKIND, &
                                 -layerThickness(kmin, iCell)/ &
                           config_flux_attenuation_coefficient_runoff))
         if (trim(config_subglacial_runoff_mode) == 'data') then
             if (config_use_sgr_opt_kpp) then
                if ( trim(config_sgr_flux_vertical_location) == 'top' ) then
                    fracAbsorbedSubglacialRunoff = 1.0_RKIND &
                          - exp( max(-100.0_RKIND, &
                                     -layerThickness(kmin, iCell)/ &
                               config_flux_attenuation_coefficient_subglacial_runoff))
                else if ( trim(config_sgr_flux_vertical_location) == 'uniform' ) then
                   ! calculate total thickness into variable zTop
                   zTop = 0.0_RKIND
                   do k = minLevelCell(iCell), maxLevelCell(iCell)
                      zTop = zTop + layerThickness(k,iCell)
                   end do
                   ! distribute flux evenly throughout water column
                   fracAbsorbedSubglacialRunoff = layerThickness(kmin, iCell) / zTop
                else if ( trim(config_sgr_flux_vertical_location) == 'bottom' ) then
                   zTop = 0.0_RKIND
                   do k = maxLevelCell(iCell), minLevelCell(iCell), -1
                      zBot = zTop - layerThickness(k,iCell)
                      if (k == minLevelCell(iCell)) then
                        transmissionCoeffTop = exp( max(zTop / config_flux_attenuation_coefficient_subglacial_runoff, -100.0_RKIND) )
                        transmissionCoeffBot = exp( max(zBot / config_flux_attenuation_coefficient_subglacial_runoff, -100.0_RKIND) )
                        fracAbsorbedSubglacialRunoff = transmissionCoeffTop - transmissionCoeffBot
                      end if
                      zTop = zBot
                   end do
                end if
             else
                   fracAbsorbedSubglacialRunoff = 0.0_RKIND
             end if
         end if

         ! Store the total tracer flux below in
         ! nonLocalSurfaceTemperatureFlux for use in the CVMix nonlocal
         ! transport code. This includes tracer forcing due to thickness

         nonLocalSurfaceTracerFlux(indexTempFlux, iCell) = &
                   activeTracersSurfaceFlux(indexTempFlux,iCell) &
                 + penetrativeTemperatureFlux(iCell) &
                 - penetrativeTemperatureFluxOBL(iCell)  &
                 - fracAbsorbed*surfaceThicknessFlux(iCell) * &
                   activeTracers(indexTempFlux,kmin,iCell)

         nonLocalSurfaceTracerFlux(indexSaltFlux,iCell) = &
                   activeTracersSurfaceFlux(indexSaltFlux,iCell) &
                 - fracAbsorbed*surfaceThicknessFlux(iCell)* &
                   activeTracers(indexSaltFlux,kmin,iCell) &
                 - fracAbsorbedRunoff*surfaceThicknessFluxRunoff(iCell)* &
                   activeTracers(indexSaltFlux,kmin,iCell)
         if (trim(config_subglacial_runoff_mode) == 'data') then
            nonLocalSurfaceTracerFlux(indexSaltFlux,iCell) = nonLocalSurfaceTracerFlux(indexSaltFlux,iCell) &
                 - fracAbsorbedSubglacialRunoff*surfaceThicknessFluxSubglacialRunoff(iCell)* &
                   activeTracers(indexSaltFlux,kmin,iCell)
         end if

         surfaceBuoyancyForcing(iCell) = &
                  thermalExpansionCoeff(kmin,iCell)* &
                  nonLocalSurfaceTracerFlux(indexTempFlux,iCell) &
                - salineContractionCoeff(kmin,iCell)* &
                  nonLocalSurfaceTracerFlux(indexSaltFlux,iCell)

         ! at this point, surfaceBuoyancyForcing has units of m/s
         ! change into units of m^2/s^3 (which can be thought of as
         ! the flux of buoyancy, units of buoyancy * velocity )
         surfaceBuoyancyForcing(iCell) = surfaceBuoyancyForcing(iCell)* &
                                         gravity

         surfaceFrictionVelocity(iCell) = &
                           sqrt(sfcStressMag(iCell) / rho_sw)

      enddo
      !$omp end do
      !$omp end parallel

      ! deallocate scratch space
      deallocate(thermalExpansionCoeff, &
                 salineContractionCoeff, &
                 densitySurfaceDisplaced)

      call mpas_timer_stop('Mixing input fields')

    !-------------------------------------------------------------------

    end subroutine ocn_compute_mixing_input_fields!}}}

!***********************************************************************
!
!  routine ocn_compute_land_ice_flux_input_fields
!
!> \brief Builds the forcing array for land-ice forcing
!> \author Xylar Asay-Davis
!> \date   09/14/2015
!> \details
!>  This routine builds surface flux arrays related to land-ice forcing.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_land_ice_flux_input_fields(layerThickness, normalVelocity, &
      activeTracers, landIceFraction, landIceFloatingMask, &
      timeLevel, indexTval, indexSval)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel, indexTval, indexSval

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness, normalVelocity

      real (kind=RKIND), dimension(:), intent(in) :: &
         landIceFraction

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         activeTracers

      integer, dimension(:), intent(in) :: &
         landIceFloatingMask

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, iEdge, cell1, cell2, iLevel, i

      real (kind=RKIND) :: blThickness, dz, weightSum, h_nu, Gamma_turb, landIceEdgeFraction, velocityMagnitude

      ! Scratch Arrays
      real (kind=RKIND), dimension(:), allocatable ::  &
         blTempScratch, blSaltScratch

      ! constants for Holland and Jenkins 1999 parameterization of the boundary layer
      real (kind=RKIND), parameter :: &
         Pr = 13.8_RKIND, &             ! the Prandtl number
         Sc = 2432.0_RKIND, &           ! the Schmidt number
         nuSaltWater = 1.95e-6_RKIND, & ! molecular viscosity of sea water (m^2/s)
         kVonKarman = 0.4_RKIND, &      ! the von Karman constant
         xiN = 0.052_RKIND              ! dimensionless planetary boundary layer constant

      integer :: indexBLTval, indexBLSval, indexHeatTransval, indexSaltTransval

      !-----------------------------------------------------------------

      if (( trim(config_land_ice_flux_mode) == "off" ) .or.  &
          ( trim(config_land_ice_flux_mode) == "pressure_only")) then
         ! nothing to do here
         return
      end if

      call mpas_timer_start("land_ice_diagnostic_fields", .false.)

      allocate(blTempScratch(nCellsAll), &
               blSaltScratch(nCellsAll))
      !$acc enter data create(blTempScratch, blSaltScratch)

      ! Compute top drag everywhere there is grounded or floating land ice
#ifdef MPAS_OPENACC
      !$acc enter data copyin(landIceFraction)

      !$acc parallel loop present(cellsOnEdge, kineticEnergyCell, minLevelCell, minLevelEdgeBot) &
      !$acc    present(landIceFraction, topDrag, normalVelocity)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(cell1, cell2, velocityMagnitude, landIceEdgeFraction)
#endif
      do iEdge = 1, nEdgesAll
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)

         ! top drag tau = - CD*|u|*u, where |u| = sqrt(2*KE) = sqrt(KE1 + KE2) from the neighboring cells
         velocityMagnitude = sqrt(kineticEnergyCell(minLevelCell(cell1),cell1) + kineticEnergyCell(minLevelCell(cell2),cell2))
         landIceEdgeFraction = 0.5_RKIND*(landIceFraction(cell1)+landIceFraction(cell2))

         topDrag(iEdge) = - rho_sw * landIceEdgeFraction * landIceTopDragCoeff &
                          * velocityMagnitude * normalVelocity(minLevelEdgeBot(iEdge),iEdge)

      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

      ! compute top drag magnitude at cell centers
#ifdef MPAS_OPENACC
      !$acc parallel loop present(kineticEnergyCell, minLevelCell, topDragMag) &
      !$acc    present(landIceFraction)
#else
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCellsAll
         ! the magnitude of the top drag is CD*u**2 = CD*(2*KE)
         topDragMag(iCell) = rho_sw * landIceFraction(iCell) &
                           * 2.0_RKIND * landIceTopDragCoeff *  kineticEnergyCell(minLevelCell(iCell),iCell)
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      if ( trim(config_land_ice_flux_mode) == "data" ) then
         ! top drag is set, and that's all we need in data (prescribed) mode
         call mpas_timer_stop("land_ice_diagnostic_fields")
         return
      end if

      ! From here on, were in "standalone" or "coupled" mode

      ! compute friction velocity at cell centers
#ifdef MPAS_OPENACC
      !$acc parallel loop present(landIceFrictionVelocity, kineticEnergyCell, minLevelCell)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCellsAll
         ! the friction velocity is the square root of the top drag + variance of tidal velocity
         ! (computed regardless of land-ice coverage)

         landIceFrictionVelocity(iCell) = sqrt(landIceTopDragCoeff * &
                                   (2.0_RKIND * kineticEnergyCell(minLevelCell(iCell),iCell) &
                                   + config_land_ice_flux_tidal_Jourdain_alpha**2 * &
                                   (config_land_ice_flux_tidal_Jourdain_A0 * velocityTidalRMS(iCell) + &
                                   config_land_ice_flux_tidal_Jourdain_U0)**2))
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

      ! average temperature and salinity over horizontal neighbors and the sub-ice-shelf boundary layer
#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(blTempScratch, blSaltScratch, activeTracers, &
      !$acc            minLevelCell, maxLevelCell, layerThickness)
#else
      !$omp do schedule(runtime) private(blThickness, iLevel, dz)
#endif
      do iCell = 1, nCellsAll
         blThickness = 0.0_RKIND
         blTempScratch(iCell) = 0.0_RKIND
         blSaltScratch(iCell) = 0.0_RKIND
         do iLevel = minLevelCell(iCell), maxLevelCell(iCell)
            dz = min(layerThickness(iLevel,iCell),config_land_ice_flux_boundaryLayerThickness-blThickness)
            if(dz <= 0.0_RKIND) exit
            blTempScratch(iCell) = blTempScratch(iCell) + activeTracers(indexTval, iLevel, iCell)*dz
            blSaltScratch(iCell) = blSaltScratch(iCell) + activeTracers(indexSval, iLevel, iCell)*dz
            blThickness = blThickness + dz
         end do
         if(blThickness > 0.0_RKIND) then
           blTempScratch(iCell) = blTempScratch(iCell)/blThickness
           blSaltScratch(iCell) = blSaltScratch(iCell)/blThickness
         end if
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(landIceBoundaryLayerTracers, &
      !$acc            blTempScratch, blSaltScratch, cellMask, &
      !$acc            nEdgesOnCell, cellsOnCell, minLevelCell)
#else
      !$omp do schedule(runtime) private(weightSum, i, cell2)
#endif
      do iCell = 1, nCellsAll
         landIceBoundaryLayerTracers(indexBLT, iCell) = blTempScratch(iCell)
         landIceBoundaryLayerTracers(indexBLS, iCell) = blSaltScratch(iCell)
         if(config_land_ice_flux_boundaryLayerNeighborWeight > 0.0_RKIND) then
            weightSum = 1.0_RKIND
            do i = 1, nEdgesOnCell(iCell)
               cell2 = cellsOnCell(i,iCell)

               landIceBoundaryLayerTracers(indexBLT, iCell) = landIceBoundaryLayerTracers(indexBLT, iCell) &
                 + cellMask(minLevelCell(cell2),cell2)*config_land_ice_flux_boundaryLayerNeighborWeight*blTempScratch(cell2)
               landIceBoundaryLayerTracers(indexBLS, iCell) = landIceBoundaryLayerTracers(indexBLS, iCell) &
                 + cellMask(minLevelCell(cell2),cell2)*config_land_ice_flux_boundaryLayerNeighborWeight*blSaltScratch(cell2)
               weightSum = weightSum + cellMask(minLevelCell(cell2),cell2)*config_land_ice_flux_boundaryLayerNeighborWeight
            end do
            landIceBoundaryLayerTracers(:, iCell) = landIceBoundaryLayerTracers(:, iCell)/weightSum
         end if
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      if(jenkinsOn) then
#ifdef MPAS_OPENACC
         !$acc parallel loop present(landIceTracerTransferVelocities, landIceFrictionVelocity)
#else
         !$omp parallel
         !$omp do schedule(runtime)
#endif
         do iCell = 1, nCellsAll
            ! transfer coefficients from namelist
            landIceTracerTransferVelocities(indexHeatTrans, iCell) = landIceFrictionVelocity(iCell) &
                                             * config_land_ice_flux_jenkins_heat_transfer_coefficient
            landIceTracerTransferVelocities(indexSaltTrans, iCell) = landIceFrictionVelocity(iCell) &
                                             * config_land_ice_flux_jenkins_salt_transfer_coefficient
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
      else if(hollandJenkinsOn) then
#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(fCell, landIceTracerTransferVelocities, &
         !$acc            landIceFrictionVelocity)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(h_nu, Gamma_turb)
#endif
         do iCell = 1, nCellsAll
            ! friction-velocity dependent non-dimensional transfer coefficients from
            ! Holland and Jenkins 1999, (14)-(16) with eta_* = 1
            h_nu = 5.0_RKIND*nuSaltWater/landIceFrictionVelocity(iCell) ! uStar should never be zero because of tidal term

            Gamma_turb = 1.0_RKIND/(2.0_RKIND*xiN) - 1.0_RKIND/kVonKarman
            if(abs(fCell(iCell)) > 0.0_RKIND) then
              Gamma_turb = Gamma_turb + 1.0_RKIND/kVonKarman*log(landIceFrictionVelocity(iCell) &
                *xiN/(abs(fCell(iCell))*h_nu))
            end if

            landIceTracerTransferVelocities(indexHeatTrans, iCell) = 1.0_RKIND/(Gamma_turb + 12.5_RKIND &
                                                                   * Pr**(2.0_RKIND/3.0_RKIND) - 6.0_RKIND)
            landIceTracerTransferVelocities(indexSaltTrans, iCell) = 1.0_RKIND/(Gamma_turb + 12.5_RKIND &
                                                                   * Sc**(2.0_RKIND/3.0_RKIND) - 6.0_RKIND)
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
      end if

      !$acc exit data delete(blTempScratch, blSaltScratch)
      deallocate(blTempScratch, &
                 blSaltScratch)

      ! modify the spatially-varying attenuation coefficient where there is land ice
#ifdef MPAS_OPENACC
      !$acc enter data copyin(landIceFloatingMask)

      !$acc parallel loop present(landIceFloatingMask, sfcFlxAttCoeff)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCellsAll
         if(landIceFloatingMask(iCell) == 1) then
            sfcFlxAttCoeff(iCell) = config_land_ice_flux_attenuation_coefficient
         end if
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      call mpas_timer_stop("land_ice_diagnostic_fields")

   !--------------------------------------------------------------------

   end subroutine ocn_compute_land_ice_flux_input_fields!}}}

!***********************************************************************
!
!  routine ocn_reconstruct_eddy_vectors
!
!> \brief   Computes cell-centered vector diagnostics
!> \author  Mark Petersen
!> \date    May 2014
!> \details
!>  This routine computes cell-centered vector diagnostics
!
!-----------------------------------------------------------------------

   subroutine ocn_reconstruct_eddy_vectors(meshPool) !{{{

      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information

      call mpas_timer_start('reconstruct eddy vecs')

      call mpas_reconstruct(meshPool, normalTransportVelocity,          &
                       transportVelocityX,            &
                       transportVelocityY,            &
                       transportVelocityZ,            &
                       transportVelocityZonal,        &
                       transportVelocityMeridional    &
                      )

      if(config_use_gm) then
         call mpas_reconstruct(meshPool, normalGMBolusVelocity,          &
                       eddyVelocityX,            &
                       eddyVelocityY,            &
                       eddyVelocityZ,            &
                       GMBolusVelocityZonal,        &
                       GMBolusVelocityMeridional    &
                      )

         call mpas_reconstruct(meshPool, gmStreamFuncTopOfEdge,          &
                      GMStreamFuncX,            &
                      GMStreamFuncY,            &
                      GMStreamFuncZ,            &
                      GMStreamFuncZonal,        &
                      GMStreamFuncMeridional    &
                     )
      end if

      if (config_submesoscale_enable) then
         call mpas_reconstruct(meshPool, normalMLEvelocity,          &
                       eddyVelocityX,            &
                       eddyVelocityY,            &
                       eddyVelocityZ,            &
                       mleVelocityZonal,        &
                       mleVelocityMeridional    &
                      )

      end if

      call mpas_timer_stop('reconstruct eddy vecs')

   end subroutine ocn_reconstruct_eddy_vectors!}}}


!***********************************************************************
!
!  routine ocn_validate_state
!
!> \brief   Ocean state validation routine
!> \author  Doug Jacobsen
!> \date    08/11/2016
!> \details
!>  This routine validates that the ocean state is able to continue running
!>  with for the next time step.
!>  If a processor detects that it is unable to continue running, some
!>  diagnostic information is written out about it.
!>  This routine relies on the definition that a NaN does not equal itself
!>  for detecting issues with the state.
!
!-----------------------------------------------------------------------
   subroutine ocn_validate_state(domain, timeLevel)!{{{
      type (domain_type), intent(inout) :: domain
      integer, intent(in), optional :: timeLevel

      integer :: timeLevelLocal

      type (block_type), pointer :: block

      integer, pointer :: nCellsSolve, nEdgesSolve, nVerticesSolve

      type (mpas_pool_type), pointer :: meshPool, statePool, tracersPool
      type (mpas_pool_type), pointer :: forcingPool

      real (kind=RKIND), dimension(:, :), pointer :: layerThickness, normalVelocity
      real (kind=RKIND), dimension(:, :, :), pointer :: activeTracers

      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeBot, maxLevelVertexBot

      real (kind=RKIND), dimension(:), pointer :: real1DArr

      integer :: iCell, iEdge, iTracer, k

      logical :: fatalErrorDetected

      logical :: thickNanFound, velNanFound, tracersNanFound

      real (kind=RKIND) :: minValue, maxValue

      integer :: debugUnit
      character (len=StrKIND) :: debugFilename, fieldName

      if ( present(timeLevel) ) then
         timeLevelLocal = timeLevel
      else
         timeLevelLocal = 1
      end if

      fatalErrorDetected = .false.

      call mpas_new_unit(debugUnit)

      block => domain % blocklist
      do while ( associated(block) )
         thickNanFound = .false.
         velNanFound = .false.
         tracersNanFound = .false.

         debugFilename = ocn_build_log_filename('mpas_ocean_block_stats_', block % blockID)

         ! Get standard fields, for checking errors
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', maxLevelEdgeBot)
         call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)

         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel=timeLevelLocal)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel=timeLevelLocal)
         call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel=timeLevelLocal)

         ! Check for errors in the state fields
         do iCell = 1, nCellsSolve
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               thickNanFound = thickNanFound .or. ( .not. layerThickness(k, iCell) == layerThickness(k, iCell) )
               do iTracer = 1, size(activeTracers, dim=1)
                  tracersNanFound = tracersNanFound .or. &
                                  ( .not. activeTracers(iTracer, k, iCell) == activeTracers(iTracer, k, iCell) )
               end do
            end do
         end do

         do iEdge = 1, nEdgesSolve
            do k = minLevelEdgeTop(iEdge), maxLevelEdgeBot(iEdge)
               velNanFound = velNanFound .or. ( .not. normalVelocity(k, iEdge) == normalVelocity(k, iEdge) )
            end do
         end do

         ! If an error was found, we need to open the file to write data out.
         if ( thickNanFound .or. tracersNanFound .or. velNanFound ) then!{{{
            open(unit=debugUnit, file=debugFilename, form='formatted', status='unknown')

            write(debugUnit, *) 'ERROR: NaN Detected in state see below for which field contained a NaN.'
            write(debugUnit, *) '  -- Statistics information for block fields'

            if ( mpas_stream_mgr_stream_exists(domain % streamManager, 'block_.*') ) then
               call mpas_stream_mgr_block_Write(domain % streamManager, writeBlock=block, streamID='block_.*', &
                  forceWriteNow=.true.)
            end if

            ! Also, write general block information, like lat/lon bounds

            ! Test latCell
            fieldName = 'latCell'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(meshPool, fieldName, real1DArr)
            do iCell = 1, nCellsSolve
               minValue = min( minValue, real1DArr(iCell) )
               maxValue = max( maxValue, real1DArr(iCell) )
            end do
            call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)

            ! Test lonCell
            fieldName = 'lonCell'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(meshPool, fieldName, real1DArr)
            do iCell = 1, nCellsSolve
               minValue = min( minValue, real1DArr(iCell) )
               maxValue = max( maxValue, real1DArr(iCell) )
            end do
            call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)

            ! Test xCell
            fieldName = 'xCell'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(meshPool, fieldName, real1DArr)
            do iCell = 1, nCellsSolve
               minValue = min( minValue, real1DArr(iCell) )
               maxValue = max( maxValue, real1DArr(iCell) )
            end do
            call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)

            ! Test yCell
            fieldName = 'yCell'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(meshPool, fieldName, real1DArr)
            do iCell = 1, nCellsSolve
               minValue = min( minValue, real1DArr(iCell) )
               maxValue = max( maxValue, real1DArr(iCell) )
            end do
            call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)

            ! Test zCell
            fieldName = 'zCell'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(meshPool, fieldName, real1DArr)
            do iCell = 1, nCellsSolve
               minValue = min( minValue, real1DArr(iCell) )
               maxValue = max( maxValue, real1DArr(iCell) )
            end do
            call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)

            ! Test areaCell
            fieldName = 'areaCell'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(meshPool, fieldName, real1DArr)
            do iCell = 1, nCellsSolve
               minValue = min( minValue, real1DArr(iCell) )
               maxValue = max( maxValue, real1DArr(iCell) )
            end do
            call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
         end if!}}}

         ! If there was a thickness NaN found, write out information about fields that affect thickness
         if ( thickNanFound ) then!{{{
            write(debugUnit, *) ''
            write(debugUnit, *) ''
            write(debugUnit, *) 'ERROR: NaN Detected in layerThickness.'
            write(debugUnit, *) '  -- Statistics information for layerThickness fields'

            ! Test seaIcePressure
            fieldName = 'seaIcePressure'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test atmosphericPressure
            fieldName = 'atmosphericPressure'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test rainFlux
            fieldName = 'rainFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test evaporationFlux
            fieldName = 'evaporationFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test snowFlux
            fieldName = 'snowFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test seaIceFreshWaterFlux
            fieldName = 'seaIceFreshWaterFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test iceRunoffFlux
            fieldName = 'iceRunoffFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test riverRunoffFlux
            fieldName = 'riverRunoffFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if
         end if!}}}

         ! If there was a velocity NaN found, write out information about fields that affect velocity
         if ( velNanFound ) then!{{{
            write(debugUnit, *) ''
            write(debugUnit, *) ''
            write(debugUnit, *) 'ERROR: NaN Detected in normalVelocity.'
            write(debugUnit, *) '  -- Statistics information for normalVelocity fields'

            ! Test seaIcePressure
            fieldName = 'seaIcePressure'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test atmosphericPressure
            fieldName = 'atmosphericPressure'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test windStressZonal
            fieldName = 'windStressZonal'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test windStressMeridional
            fieldName = 'windStressMeridional'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test surfaceStressMagnitude
            fieldName = 'surfaceStressMagnitude'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test surfaceStress
            fieldName = 'surfaceStress'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iEdge = 1, nEdgesSolve
                  minValue = min( minValue, real1DArr(iEdge) )
                  maxValue = max( maxValue, real1DArr(iEdge) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test angleEdge
            fieldName = 'angleEdge'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(meshPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iEdge = 1, nEdgesSolve
                  minValue = min( minValue, real1DArr(iEdge) )
                  maxValue = max( maxValue, real1DArr(iEdge) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if
         end if!}}}

         ! If there was a tracers NaN found, write out information about fields that affect tracers
         if ( tracersNanFound ) then!{{{
            write(debugUnit, *) ''
            write(debugUnit, *) ''
            write(debugUnit, *) 'ERROR: NaN Detected in activeTracers.'
            write(debugUnit, *) '  -- Statistics information for activeTracers fields'

            ! Test latentHeatFlux
            fieldName = 'latentHeatFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test sensibleHeatFlux
            fieldName = 'sensibleHeatFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test longWaveHeatFluxUp
            fieldName = 'longWaveHeatFluxUp'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test longWaveHeatFluxDown
            fieldName = 'longWaveHeatFluxDown'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test seaIceHeatFlux
            fieldName = 'seaIceHeatFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test icebergHeatFlux
            fieldName = 'icebergHeatFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test seaIceFreshWaterFlux
            fieldName = 'seaIceFreshWaterFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test icebergFreshWaterFlux
            fieldName = 'icebergFreshWaterFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test rainFlux
            fieldName = 'rainFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test snowFlux
            fieldName = 'snowFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test iceRunoffFlux
            fieldName = 'iceRunoffFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test riverRunoffFlux
            fieldName = 'riverRunoffFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            if (trim(config_subglacial_runoff_mode) == 'data') then
                ! Test subglacialRunoffFlux
                fieldName = 'subglacialRunoffFlux'
                minValue = HUGE(minValue)
                maxValue = -HUGE(maxValue)
                call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
                if ( associated(real1DArr) ) then
                   do iCell = 1, nCellsSolve
                      minValue = min( minValue, real1DArr(iCell) )
                      maxValue = max( maxValue, real1DArr(iCell) )
                   end do
                   call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
                end if
            end if

            ! Test seaIceSalinityFlux
            fieldName = 'seaIceSalinityFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if

            ! Test shortWaveHeatFlux
            fieldName = 'shortWaveHeatFlux'
            minValue = HUGE(minValue)
            maxValue = -HUGE(maxValue)
            call mpas_pool_get_array(forcingPool, fieldName, real1DArr)
            if ( associated(real1DArr) ) then
               do iCell = 1, nCellsSolve
                  minValue = min( minValue, real1DArr(iCell) )
                  maxValue = max( maxValue, real1DArr(iCell) )
               end do
               call ocn_write_field_statistics(debugUnit, fieldName, minValue, maxValue)
            end if
         end if!}}}

         ! If an error was found, we need to close the file now.
         if ( thickNanFound .or. tracersNanFound .or. velNanFound ) then!{{{
            flush(debugUnit)
            close(debugUnit)
            fatalErrorDetected = .true.
         end if!}}}

         block => block % next
      end do

      call mpas_release_unit(debugUnit)

      if ( fatalErrorDetected ) then
         call mpas_log_write( &
            'ERROR: State validation failed. See block stats files for more information.', &
            MPAS_LOG_CRIT)
      end if

   end subroutine ocn_validate_state!}}}

   function ocn_build_log_filename(prefix, identifier) result(filename)!{{{
      character (len=*), intent(in) :: prefix
      integer, intent(in) :: identifier

      character (len=StrKIND) :: filename

      character (len=StrKIND) :: identifierString

      if ( identifier .lt. 10 ) then
         write(identifierString, '(I1)') identifier
      else if ( identifier .lt. 100 ) then
         write(identifierString, '(I2)') identifier
      else if ( identifier .lt. 1000 ) then
         write(identifierString, '(I3)') identifier
      else if ( identifier .lt. 10000 ) then
         write(identifierString, '(I4)') identifier
      else if ( identifier .lt. 100000 ) then
         write(identifierString, '(I5)') identifier
      else if ( identifier .lt. 1000000 ) then
         write(identifierString, '(I6)') identifier
      else if ( identifier .lt. 10000000 ) then
         write(identifierString, '(I7)') identifier
      else
         write(identifierString, '(I99)') identifier
      end if

      filename = trim(prefix) // trim(identifierString)

   end function ocn_build_log_filename!}}}

   subroutine ocn_write_field_statistics(unitNumber, fieldName, minValue, maxValue)
      integer, intent(in) :: unitNumber
      character (len=*), intent(in) :: fieldName
      real (kind=RKIND), intent(in) :: minValue, maxValue

      write(unitNumber, *) '    Field: ', trim(fieldName)
      write(unitNumber, *) '        Min: ', minValue
      write(unitNumber, *) '        Max: ', maxValue

   end subroutine ocn_write_field_statistics

!***********************************************************************
!
!  routine ocn_diagnostics_init
!
!> \brief   Initializes flags used within diagnostics routines.
!> \author  Mark Petersen
!> \date    4 November 2011
!> \details
!>  This routine initializes flags related to quantities computed within
!>  other diagnostics routines.
!
!-----------------------------------------------------------------------
   subroutine ocn_diagnostics_init(domain, err)!{{{
      type(domain_type), intent(in) :: domain
      integer, intent(out) :: err !< Output: Error flag

      err = 0

      if(config_include_KE_vertex) then
         ke_vertex_flag = 1
         ke_cell_flag = 0
      else
         ke_vertex_flag = 0
         ke_cell_flag = 1
      endif

      if ( (trim(config_time_integrator) == 'RK4') &
           .or. (trim(config_time_integrator) == 'LTS') &
           .or. (trim(config_time_integrator) == 'FB_LTS') ) then
         ! For RK4, LTS, or FB_LTS, PV includes f: PV = (eta+f)/h.
         fCoef = 1
      elseif (trim(config_time_integrator) == 'split_explicit' &
        .or.trim(config_time_integrator) == 'unsplit_explicit' &
        .or.trim(config_time_integrator) == 'split_implicit'   &
        .or.trim(config_time_integrator) == 'split_explicit_ab2') then
          ! For split explicit, PV is eta/h because the Coriolis term
          ! is added separately to the momentum tendencies.
          fCoef = 0
      end if

      ! Initialize choice for computing thickness at edges for fluxes
      if (trim(config_thickness_flux_type) == 'centered') then
         thickEdgeFluxChoice = thickEdgeFluxCenter
      elseif (trim(config_thickness_flux_type) == 'constant') then
         thickEdgeFluxChoice = thickEdgeFluxConstant
         if (.not.config_disable_vel_hadv) then
            call mpas_log_write('config_thickness_flux_type set to ' //&
                 & trim(config_thickness_flux_type) // ' but '//&
                 & 'config_disable_vel_hadv is false. '//&
                 & 'Set config_disable_vel_hadv to true for a linear test case. ', &
                 MPAS_LOG_CRIT)
         end if
      elseif (trim(config_thickness_flux_type) == 'upwind') then
         thickEdgeFluxChoice = thickEdgeFluxUpwind
         call mpas_log_write('Thickness flux option set to ' //&
              & trim(config_thickness_flux_type))
      else
         call mpas_log_write('Thickness flux option of ' //&
              & trim(config_thickness_flux_type) // &
              & 'is not known', MPAS_LOG_CRIT)
      end if

      ! Initialize coefficient choice for computing land ice drag
      if (config_use_implicit_top_drag) then
         landIceTopDragCoeff = config_implicit_top_drag_coeff
      else
         landIceTopDragCoeff = config_land_ice_flux_explicit_topDragCoeff
      endif

      if (trim(config_thickness_drag_type) == 'centered') then
         thickEdgeDragChoice = thickEdgeDragCenter
      elseif (trim(config_thickness_drag_type) == 'harmonic-mean') then
         thickEdgeDragChoice = thickEdgeDragHarMean
      else
         call mpas_log_write('Thickness drag option of ' //&
              & trim(config_thickness_drag_type) // &
              & ' is not known', MPAS_LOG_CRIT)
      end if

      call ocn_diagnostics_variables_init(domain, jenkinsOn, &
                                          hollandJenkinsOn, err)

      if ( trim(config_sgr_flux_vertical_location) /= 'top' .and. &
           trim(config_sgr_flux_vertical_location) /= 'uniform' .and. &
           trim(config_sgr_flux_vertical_location) /= 'bottom' ) then
        call mpas_log_write("config_sgr_flux_vertical_location not one of 'bottom', 'uniform', 'top'.", MPAS_LOG_CRIT)
      end if

    end subroutine ocn_diagnostics_init!}}}

!***********************************************************************

end module ocn_diagnostics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
