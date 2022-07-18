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
             ocn_relativeVorticity_circulation, &
             ocn_vert_transport_velocity_top, &
             ocn_fuperp, &
             ocn_filter_btr_mode_tend_vel, &
             ocn_reconstruct_gm_vectors, &
             ocn_compute_kpp_input_fields, &
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
   real (kind=RKIND), pointer ::  coef_3rd_order

   ! Methods for computing thickness at edges for flux calculations
   integer :: thickEdgeFluxChoice  ! choice of thickness flux type
   integer, parameter :: &
      thickEdgeFluxUnknown = 0, &! type unknown or unset
      thickEdgeFluxCenter  = 1, &! use mean thickness of cell neighbors
      thickEdgeFluxUpwind  = 2   ! use upwind cell thickness at edge

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

   subroutine ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, scratchPool, tracersPool, &
                                   timeLevelIn, full)!{{{

      real (kind=RKIND), intent(in) :: dt !< Input: Time step
      type (mpas_pool_type), intent(in) :: statePool !< Input: State information
      type (mpas_pool_type), intent(in) :: forcingPool !< Input: Forcing information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
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
        layerThickness, layerThicknessLag, normalVelocity

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      integer, pointer :: indexTemperature, indexSalinity

      logical :: full_compute = .true.

      real (kind=RKIND), dimension(:), pointer :: &
        frazilSurfacePressure, landIcePressure, landIceDraft, landIceFraction

      integer, dimension(:), pointer :: &
        landIceMask

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

      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)

      call mpas_pool_get_array(forcingPool, 'seaIcePressure', seaIcePressure)
      call mpas_pool_get_array(forcingPool, 'atmosphericPressure', atmosphericPressure)
      call mpas_pool_get_array(forcingPool, 'frazilSurfacePressure', frazilSurfacePressure)

      if (landIcePressureOn) then
         call mpas_pool_get_array(forcingPool, 'landIcePressure', landIcePressure)
         call mpas_pool_get_array(forcingPool, 'landIceDraft', landIceDraft)
         call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)
         call mpas_pool_get_array(forcingPool, 'landIceFraction', landIceFraction)
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

      ! inputs: layerThickness, normalVelocity
      ! output: layerThickEdgeMean, layerThickEdgeFlux
      call ocn_diagnostic_solve_layerThicknessEdge(normalVelocity, &
                                                   layerThickness)

      ! inputs: normalVelocity
      ! outputs: relativeVorticity, circulation
      call ocn_relativeVorticity_circulation(relativeVorticity, circulation, normalVelocity, err)

      ! inputs: relativeVorticity, layerThickEdgeFlux, normalVelocity, normalTransportVelocity, 
      !         normalGMBolusVelocity
      ! outputs: vertTransportVelocityTop, vertGMBolusVelocityTop, relativeVorticityCell, 
      !          divergence, kineticEnergyCell, tangentialVelocity
      ! in and out: vertVelocityTop
      call ocn_diagnostic_solve_vortVel(relativeVorticity, &
             layerThickEdgeFlux, normalVelocity, normalTransportVelocity, &
             normalGMBolusVelocity, vertVelocityTop, vertTransportVelocityTop, &
             vertGMBolusVelocityTop, relativeVorticityCell, divergence, &
             kineticEnergyCell, tangentialVelocity)

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
      ! Z-coordinates
      ! This section must be placed in the code before computing the density.
      !
      ! inputs : layerThickness
      ! outputs : zMid, zTop, ssh)
      call ocn_diagnostic_solve_z_coordinates(layerThickness, zMid, zTop, ssh)

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

      end if ! full_compute

      !
      !  compute fields needed to compute land-ice fluxes, either in the ocean model or in the coupler
      !
      ! inputs: layerThickness, normalVelocity, landIceFraction, landIceMask
      ! outputs: 
      if (landIcePressureOn) then
         call ocn_compute_land_ice_flux_input_fields(layerThickness, normalVelocity, activeTracers, &
                landIceFraction, landIceMask, timeLevel, indexTemperature, indexSalinity)
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
      !! ocn_diagnostic_solve_layerThicknessEdge :: layerThickEdgeFlux, (diagnostics array)
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
      !! ocn_diagnostic_solve_ssh                :: pressureAdjustedSSH, (diagnostics array)
      !!                                            gradSSH (diagnostics array)

!     !$acc update host(layerThickEdgeFlux, layerThickEdgeCenter)
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
!  routine ocn_diagnostic_solve_layerThicknessEdge
!
!> \brief   Computes layer thickness at edge
!> \author  Matt Turner
!> \date    October 2020
!> \details
!>  This routine computes the diagnostic variables layerThickEdgeFlux
!>  and layerThickEdgeMean. The Mean is a simple mean across the edge
!>  but the Flux value can either be the mean or an upwind value,
!>  depending on user input.
!
!-----------------------------------------------------------------------

   subroutine ocn_diagnostic_solve_layerThicknessEdge(normalVelocity, &
                                                      layerThickness)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity,  &!< [in] transport
         layerThickness    !< [in] layer thickness at cell center

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      ! Outputs are the shared variables from diagnostic_variables:
      !real (kind=RKIND), dimension(:,:), intent(out) :: &
      !   layerThickEdgeMean, & !< [out] centered layer thickness on edge
      !   layerThickEdgeFlux    !< [out] flux-related thickness on edge

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

      case (thickEdgeFluxUpwind)
         ! Use upwind thickness as the edge flux value

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
!> \brief   Computes diagnostic variable kineticEnergyCell
!> \author  Matt Turner
!> \date    October 2020
!> \details
!>  This routine computes the diagnostic variable kineticEnergyCell
!
!-----------------------------------------------------------------------
   subroutine ocn_diagnostic_solve_kineticEnergy(normalVelocity, &
                kineticEnergyCell)!{{{

      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity     !< Input: transport

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         kineticEnergyCell

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iVertex, i, iEdge, k, j, nCells, nVertices, iCell
      real (kind=RKIND) :: r_tmp, invAreaCell1

      ! Scratch Arrays
      ! kineticEnergyVertex: kinetic energy of horizontal velocity defined at vertices
      !               units: m^2 s^{-2}
      real (kind=RKIND), dimension(:,:), allocatable :: kineticEnergyVertex
      ! kineticEnergyVertexOnCells: kinetic energy of horizontal velocity defined at vertices
      !                      units: m^2 s^{-2}
      real (kind=RKIND), dimension(:,:), allocatable :: kineticEnergyVertexOnCells
#ifdef MPAS_OPENACC
      !$acc declare device_resident(kineticEnergyVertex, kineticEnergyVertexOnCells)
#endif

      nCells = nCellsHalo( 2 )
      nVertices = nVerticesAll

      allocate(kineticEnergyVertex(nVertLevels, nVertices), &
               kineticEnergyVertexOnCells(nVertLevels, nCells))

#ifdef MPAS_OPENACC
      !$acc parallel loop present(areaTriangle, dcEdge, dvEdge, normalVelocity, edgesOnVertex)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(i, iEdge, r_tmp, k)
#endif
      do iVertex = 1, nVertices
        kineticEnergyVertex(:, iVertex) = 0.0_RKIND
        do i = 1, vertexDegree
          iEdge = edgesOnVertex(i, iVertex)
          r_tmp = dcEdge(iEdge) * dvEdge(iEdge) * 0.25_RKIND / areaTriangle(iVertex)
#ifdef MPAS_OPENACC
          !$acc loop seq
#endif
          do k = 1, nVertLevels
            kineticEnergyVertex(k, iVertex) = kineticEnergyVertex(k, iVertex) + r_tmp * normalVelocity(k, iEdge)**2
          end do
        end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

#ifdef MPAS_OPENACC
      !$acc parallel loop present(kiteIndexOnCell, kiteAreasOnVertex, invAreaCell, verticesOnCell, nEdgesOnCell)
#else
      !$omp do schedule(runtime) private(invAreaCell1, i, j, iVertex, k)
#endif
      do iCell = 1, nCells
        kineticEnergyVertexOnCells(:, iCell) = 0.0_RKIND
        invAreaCell1 = invAreaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          j = kiteIndexOnCell(i, iCell)
          iVertex = verticesOnCell(i, iCell)
#ifdef MPAS_OPENACC
          !$acc loop seq
#endif
          do k = 1, nVertLevels
            kineticEnergyVertexOnCells(k, iCell) = kineticEnergyVertexOnCells(k, iCell) + kiteAreasOnVertex(j, iVertex) &
                                                 * kineticEnergyVertex(k, iVertex) * invAreaCell1
          end do
        end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

      !
      ! Compute kinetic energy in each cell by blending kineticEnergyCell and kineticEnergyVertexOnCells
      !
#ifdef MPAS_OPENACC
      !$acc parallel loop collapse(2) present(kineticEnergyCell)
#else
      !$omp do schedule(runtime) private(k)
#endif
      do iCell = 1, nCells
         do k = 1, nVertLevels
            kineticEnergyCell(k,iCell) = 5.0_RKIND / 8.0_RKIND * kineticEnergyCell(k,iCell) + 3.0_RKIND / 8.0_RKIND &
                                       * kineticEnergyVertexOnCells(k,iCell)
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      deallocate(kineticEnergyVertex, &
                 kineticEnergyVertexOnCells)

   end subroutine ocn_diagnostic_solve_kineticEnergy!}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve_vorticity
!
!> \brief   Computes diagnostic variables for vorticity
!> \author  Matt Turner
!> \date    October 2020
!> \details
!>  This routine computes the diagnostic variables for vorticity
!
!-----------------------------------------------------------------------
   subroutine ocn_diagnostic_solve_vorticity(dt, normalVelocity, tangentialVelocity, &
                layerThickness, relativeVorticity, &
                normalizedRelativeVorticityEdge, normalizedPlanetaryVorticityEdge, &
                normalizedRelativeVorticityCell)!{{{

      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: dt !< Input: Time step
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity     !< Input: transport
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         tangentialVelocity !< Input: transport
      real (kind=RKIND), dimension(:,:), intent(in) :: &
        layerThickness      !< Input: layer thickness
      real (kind=RKIND), dimension(:,:), intent(in) :: &
        relativeVorticity

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
        normalizedRelativeVorticityEdge
      real (kind=RKIND), dimension(:,:), intent(out) :: &
        normalizedPlanetaryVorticityEdge
      real (kind=RKIND), dimension(:,:), intent(out) :: &
        normalizedRelativeVorticityCell

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: nVertices, nEdges, nCells
      integer :: iVertex, k, i, iEdge, iCell, j
      integer :: vertex1, vertex2, cell1, cell2
      real(kind=RKIND) :: invAreaTri1, invAreaCell1, invLength
      real(kind=RKIND) :: layerThicknessVertex
      real(kind=RKIND) :: apvm_scale_factor

      ! normalizedPlanetaryVorticityVertex: earth's rotational rate (Coriolis parameter, f)
      !                                     divided by layer thickness, defined at vertices
      !                              units: s^{-1}
      real (kind=RKIND), dimension(:,:), allocatable :: normalizedPlanetaryVorticityVertex
      ! normalizedRelativeVorticityVertex: curl of horizontal velocity divided by layer thickness,
      !                                    defined at vertices
      !                             units: s^{-1}
      real (kind=RKIND), dimension(:,:), allocatable :: normalizedRelativeVorticityVertex
      ! vorticityGradientNormalComponent: gradient of vorticity in the normal direction
      !                                   (positive points from cell1 to cell2)
      !                            units: s^{-1} m^{-1}
      real (kind=RKIND), dimension(:,:), allocatable :: vorticityGradientNormalComponent
      ! vorticityGradientTangentialComponent: gradient of vorticity in the tangent direction
      !                                       (positive points from vertex1 to vertex2)
      !                                units: s^{-1} m^{-1}
      real (kind=RKIND), dimension(:,:), allocatable :: vorticityGradientTangentialComponent
#ifdef MPAS_OPENACC
      !$acc declare device_resident(normalizedPlanetaryVorticityVertex, &
      !$acc                         normalizedRelativeVorticityVertex, &
      !$acc                         vorticityGradientNormalComponent, &
      !$acc                         vorticityGradientTangentialComponent)
#endif

      !
      ! Compute normalized relative and planetary vorticity
      !
      allocate(normalizedPlanetaryVorticityVertex(nVertLevels, nVerticesAll), &
               normalizedRelativeVorticityVertex(nVertLevels, nVerticesAll))

      nVertices = nVerticesHalo( 2 )

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc          present(areaTriangle, maxLevelVertexBot, &
      !$acc                  layerThickness, kiteAreasOnVertex, &
      !$acc                  relativeVorticity, fVertex, cellsOnVertex) &
      !$acc          private(invAreaTri1, layerThicknessVertex, k, i)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(invAreaTri1, k, layerThicknessVertex, i)
#endif
      do iVertex = 1, nVertices
         invAreaTri1 = 1.0_RKIND / areaTriangle(iVertex)
         do k = 1, maxLevelVertexBot(iVertex)
            layerThicknessVertex = 0.0_RKIND
            do i = 1, vertexDegree
               layerThicknessVertex = layerThicknessVertex + layerThickness(k,cellsOnVertex(i,iVertex)) &
                                    * kiteAreasOnVertex(i,iVertex)
            end do
            layerThicknessVertex = layerThicknessVertex * invAreaTri1
            if (layerThicknessVertex == 0) cycle

            normalizedRelativeVorticityVertex(k,iVertex) = relativeVorticity(k,iVertex) / layerThicknessVertex
            normalizedPlanetaryVorticityVertex(k,iVertex) = fVertex(iVertex) / layerThicknessVertex
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      nEdges = nEdgesHalo( 2 )

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc          present(normalizedRelativeVorticityEdge, &
      !$acc                  normalizedPlanetaryVorticityEdge, &
      !$acc                  verticesOnEdge, maxLevelEdgeBot, &
      !$acc                  minLevelEdgeTop) &
      !$acc          private(vertex1, vertex2, k)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(vertex1, vertex2, k)
#endif
      do iEdge = 1, nEdges
        normalizedRelativeVorticityEdge(:, iEdge) = 0.0_RKIND
        normalizedPlanetaryVorticityEdge(:, iEdge) = 0.0_RKIND
        vertex1 = verticesOnEdge(1, iEdge)
        vertex2 = verticesOnEdge(2, iEdge)
#ifdef MPAS_OPENACC
        !$acc loop seq
#endif
        do k = minLevelEdgeTop(iEdge), maxLevelEdgeBot(iEdge)
              normalizedRelativeVorticityEdge(k, iEdge) = 0.5_RKIND * (normalizedRelativeVorticityVertex(k, vertex1) &
                                                    + normalizedRelativeVorticityVertex(k, vertex2))
              normalizedPlanetaryVorticityEdge(k, iEdge) = 0.5_RKIND * (normalizedPlanetaryVorticityVertex(k, vertex1) &
                                                    + normalizedPlanetaryVorticityVertex(k, vertex2))
        end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      nCells = nCellsHalo ( 1 )

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc          present(normalizedRelativeVorticityCell, maxLevelCell, nEdgesOnCell, &
      !$acc                  areaCell, kiteIndexOnCell, verticesOnCell, kiteAreasOnVertex, &
      !$acc                  invAreaCell, minLevelCell) &
      !$acc          private(invAreaCell1, i, j, iVertex, k)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(invAreaCell1, i, j, iVertex, k)
#endif
      do iCell = 1, nCells
        normalizedRelativeVorticityCell(:, iCell) = 0.0_RKIND
        invAreaCell1 = invAreaCell(iCell)

        do i = 1, nEdgesOnCell(iCell)
          j = kiteIndexOnCell(i, iCell)
          iVertex = verticesOnCell(i, iCell)
#ifdef MPAS_OPENACC
        !$acc loop seq
#endif
          do k = minLevelCell(iCell), maxLevelCell(iCell)
            normalizedRelativeVorticityCell(k, iCell) = normalizedRelativeVorticityCell(k, iCell) &
              + kiteAreasOnVertex(j, iVertex) * normalizedRelativeVorticityVertex(k, iVertex) * invAreaCell1
          end do
        end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      nCells = nCellsAll
      nVertices = nVerticesAll
      nEdges = nCellsAll
      apvm_scale_factor = config_apvm_scale_factor;

      ! Diagnostics required for the Anticipated Potential Vorticity Method (apvm).
      if (config_apvm_scale_factor>1e-10) then

         allocate(vorticityGradientNormalComponent(nVertLevels, nEdges), &
                  vorticityGradientTangentialComponent(nVertLevels, nEdges))

         nEdges = nEdgesHalo( 1 )

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(minLevelEdgeBot, maxLevelEdgeTop, cellsOnEdge, &
         !$acc            minLevelEdgeTop, maxLevelEdgeBot, dcEdge, dvEdge, &
         !$acc            verticesOnedge, normalizedRelativeVorticityCell)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(cell1, cell2, vertex1, vertex2, invLength, k)
#endif
         do iEdge = 1,nEdges
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            vertex1 = verticesOnedge(1, iEdge)
            vertex2 = verticesOnedge(2, iEdge)

            invLength = 1.0_RKIND / dcEdge(iEdge)
            ! Compute gradient of PV in normal direction
            !   ( this computes the gradient for all edges bounding real cells )
#ifdef MPAS_OPENACC
        !$acc loop seq
#endif
            do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
               vorticityGradientNormalComponent(k,iEdge) = &
                  (normalizedRelativeVorticityCell(k,cell2) - normalizedRelativeVorticityCell(k,cell1)) * invLength
            enddo

            invLength = 1.0_RKIND / dvEdge(iEdge)
            ! Compute gradient of PV in the tangent direction
            !   ( this computes the gradient at all edges bounding real cells and distance-1 ghost cells )
#ifdef MPAS_OPENACC
        !$acc loop seq
#endif
            do k = minLevelEdgeTop(iEdge), maxLevelEdgeBot(iEdge)
              vorticityGradientTangentialComponent(k,iEdge) = &
                 (normalizedRelativeVorticityVertex(k,vertex2) - normalizedRelativeVorticityVertex(k,vertex1)) * invLength
            enddo

         enddo
#ifndef MPAS_OPENACC
         !$omp end do
#endif

         !
         ! Modify PV edge with upstream bias.
         !
#ifdef MPAS_OPENACC
         !$acc parallel loop, &
         !$acc    present(minLevelEdgeTop, maxLevelEdgeBot, tangentialVelocity, &
         !$acc            normalVelocity, normalizedRelativeVorticityEdge)
#else
         !$omp do schedule(runtime) private(k)
#endif
         do iEdge = 1,nEdges
#ifdef MPAS_OPENACC
        !$acc loop seq
#endif
            do k = minLevelEdgeTop(iEdge), maxLevelEdgeBot(iEdge)
              normalizedRelativeVorticityEdge(k,iEdge) = normalizedRelativeVorticityEdge(k,iEdge) &
                - apvm_scale_factor * dt * &
                    (  normalVelocity(k,iEdge)     * vorticityGradientNormalComponent(k,iEdge)      &
                     + tangentialVelocity(k,iEdge) * vorticityGradientTangentialComponent(k,iEdge) )
            end do
         enddo
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

         deallocate(vorticityGradientNormalComponent, &
                    vorticityGradientTangentialComponent)

      endif

      deallocate(normalizedPlanetaryVorticityVertex, &
                 normalizedRelativeVorticityVertex)


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

      real(kind=RKIND) :: coef, shearSquared, delU2, shearMean
      integer :: iCell, k, i, iEdge
      integer :: nCells

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
      !$acc          private(k, i, iEdge, delU2, shearSquared, shearMean)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k, shearSquared, i, iEdge, delU2, shearMean)
#endif
      do iCell=1,nCells
         RiTopOfCell(:,iCell) = 100.0_RKIND
         do k = minLevelCell(iCell)+1, maxLevelCell(iCell)
           shearSquared = 0.0_RKIND
           do i = 1, nEdgesOnCell(iCell)
             iEdge = edgesOnCell(i, iCell)
             delU2 = (normalVelocity(k-1,iEdge) - normalVelocity(k,iEdge))**2
             shearSquared = shearSquared + edgeAreaFractionOfCell(i,iCell) * delU2
           enddo
           ! Note that the factor of two is from averaging dot product to cell center on a C-grid
           shearMean = sqrt(2.0_RKIND*shearSquared )
           shearMean = shearMean / (zMid(k-1,iCell) - zMid(k,iCell))
           RiTopOfCell(k,iCell) = BruntVaisalaFreqTop(k,iCell) / (shearMean**2 + 1.0e-10_RKIND)
          end do
          RiTopOfCell(minLevelCell(iCell),iCell) = RiTopOfCell(minLevelCell(iCell)+1,iCell)
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

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
   subroutine ocn_diagnostic_solve_surfaceLayer(layerThickness, activeTracers, &
                layerThicknessEdgeFlux, normalVelocity, tracersSurfaceLayerValue,  &
                indexSurfaceLayerDepth, normalVelocitySurfaceLayer, &
                surfaceFluxAttenuationCoefficient, surfaceFluxAttenuationCoefficientRunoff)!{{{

      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        layerThickness      !< Input: layer thickness
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
        activeTracers
      real (kind=RKIND), dimension(:,:), intent(in) :: &
        layerThicknessEdgeFlux
      real (kind=RKIND), dimension(:,:), intent(in) :: &
        normalVelocity

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         tracersSurfaceLayerValue
      real (kind=RKIND), dimension(:), intent(out) :: &
         indexSurfaceLayerDepth
      real (kind=RKIND), dimension(:), intent(out) :: &
         normalVelocitySurfaceLayer
      real (kind=RKIND), dimension(:), intent(out) :: &
         surfaceFluxAttenuationCoefficient
      real (kind=RKIND), dimension(:), intent(out) :: &
         surfaceFluxAttenuationCoefficientRunoff

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: nCells, nEdges
      integer :: iCell, iEdge
      integer :: k
      real(kind=RKIND) :: rSurfaceLayer
      integer :: cell1, cell2
      real(kind=RKIND) :: surfaceLayerDepth
      real(kind=RKIND) :: sumSurfaceLayer

      !
      ! average tracer values over the ocean surface layer
      ! the ocean surface layer is generally assumed to be about 0.1 of the boundary layer depth
      if(config_use_cvmix_kpp) then

        nCells = nCellsOwned
        surfaceLayerDepth = config_cvmix_kpp_surface_layer_averaging
#ifdef MPAS_OPENACC
        !$acc parallel loop gang vector &
        !$acc          present(tracersSurfaceLayerValue, indexSurfaceLayerDepth, maxLevelCell, &
        !$acc                  layerThickness, activeTracers, minLevelCell) &
        !$acc          private(k, rSurfaceLayer)
#else
        !$omp parallel
        !$omp do schedule(runtime) private(sumSurfaceLayer, k, rSurfaceLayer)
#endif
        do iCell=1,nCells
          sumSurfaceLayer=0.0_RKIND
          tracersSurfaceLayerValue(:,iCell) = 0.0_RKIND
          indexSurfaceLayerDepth(iCell) = -9.e30
#ifdef MPAS_OPENACC
          !$acc loop seq
#endif
          do k = minLevelCell(iCell), maxLevelCell(iCell)
           sumSurfaceLayer = sumSurfaceLayer + layerThickness(k,iCell)
           rSurfaceLayer = maxLevelCell(iCell)
           if(sumSurfaceLayer.gt.surfaceLayerDepth) then
             sumSurfaceLayer = sumSurfaceLayer - layerThickness(k,iCell)
             rSurfaceLayer = int(k-1) + (surfaceLayerDepth-sumSurfaceLayer)/layerThickness(k,iCell)
             indexSurfaceLayerDepth(iCell) = rSurfaceLayer
             exit
           endif
          end do
          tracersSurfaceLayerValue(:, iCell) = 0.0_RKIND
#ifdef MPAS_OPENACC
          !$acc loop seq
#endif
          do k=1,int(rSurfaceLayer)
            tracersSurfaceLayerValue(:,iCell) = tracersSurfaceLayerValue(:,iCell) + activeTracers(:,k,iCell) &
                                              * layerThickness(k,iCell)
          enddo
          k=min( int(rSurfaceLayer)+1, maxLevelCell(iCell) )
          tracersSurfaceLayerValue(:,iCell) = (tracersSurfaceLayerValue(:,iCell) + fraction(rSurfaceLayer) &
                                            * activeTracers(:,k,iCell) * layerThickness(k,iCell)) / surfaceLayerDepth
        enddo
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif

        nEdges = nEdgesOwned

        !
        ! average normal velocity values over the ocean surface layer
        ! the ocean surface layer is generally assumed to be about 0.1 of the boundary layer depth
        !

        surfaceLayerDepth = config_cvmix_kpp_surface_layer_averaging
#ifdef MPAS_OPENACC
        !$acc parallel loop gang vector &
        !$acc          present(normalVelocitySurfaceLayer, cellsOnEdge, maxLevelEdgeTop, &
        !$acc                  layerThicknessEdgeFlux, normalVelocity, minLevelEdgeBot) &
        !$acc          private(cell1, cell2, rSurfaceLayer, sumSurfaceLayer, k)
#else
        !$omp parallel
        !$omp do schedule(runtime) private(cell1, cell2, sumSurfaceLayer, k, rSurfaceLayer)
#endif
        do iEdge=1,nEdges
          normalVelocitySurfaceLayer(iEdge) = 0.0_RKIND
          cell1=cellsOnEdge(1,iEdge)
          cell2=cellsOnEdge(2,iEdge)
          sumSurfaceLayer=0.0_RKIND
          rSurfaceLayer = min(minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge))
#ifdef MPAS_OPENACC
          !$acc loop seq
#endif
          do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
           rSurfaceLayer = k
           sumSurfaceLayer = sumSurfaceLayer + layerThicknessEdgeFlux(k,iEdge)
           if(sumSurfaceLayer.gt.surfaceLayerDepth) then
             sumSurfaceLayer = sumSurfaceLayer - layerThicknessEdgeFlux(k,iEdge)
             rSurfaceLayer = int(k-1) + (surfaceLayerDepth-sumSurfaceLayer)/layerThicknessEdgeFlux(k,iEdge)
             exit
           endif
          end do
          normalVelocitySurfaceLayer(iEdge) = 0.0_RKIND
#ifdef MPAS_OPENACC
          !$acc loop seq
#endif
          do k = minLevelEdgeBot(iEdge),int(rSurfaceLayer)
            normalVelocitySurfaceLayer(iEdge) = normalVelocitySurfaceLayer(iEdge) + normalVelocity(k,iEdge) &
                                              * layerThicknessEdgeFlux(k,iEdge)
          end do
          k=int(rSurfaceLayer)+1
          if(k.le.maxLevelEdgeTop(iEdge)) then
            normalVelocitySurfaceLayer(iEdge) = (normalVelocitySurfaceLayer(iEdge) + fraction(rSurfaceLayer) &
                                              * normalVelocity(k,iEdge) * layerThicknessEdgeFlux(k,iEdge)) / surfaceLayerDepth
          end if
        enddo
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif

      endif

      nCells = nCellsHalo( 1 )

      ! compute the attenuation coefficient for surface fluxes
#ifdef MPAS_OPENACC
      !$acc parallel loop present(surfaceFluxAttenuationCoefficient, surfaceFluxAttenuationCoefficientRunoff)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCells
         surfaceFluxAttenuationCoefficient(iCell) = config_flux_attenuation_coefficient
         surfaceFluxAttenuationCoefficientRunoff(iCell) = config_flux_attenuation_coefficient_runoff
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

   end subroutine ocn_diagnostic_solve_surfaceLayer!}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve_vortVel
!
!> \brief   Computes diagnostic variables for vorticity and velocity
!> \author  Matt Turner
!> \date    October 2020
!> \details
!>  This routine computes the diagnostic variables for vorticity and vertical
!>    velocities
!
!-----------------------------------------------------------------------
   subroutine ocn_diagnostic_solve_vortVel(relativeVorticity, &
                layerThicknessEdgeFlux, normalVelocity, normalTransportVelocity, &
                normalGMBolusVelocity, vertVelocityTop, vertTransportVelocityTop, &
                vertGMBolusVelocityTop, relativeVorticityCell, divergence, &
                kineticEnergyCell, tangentialVelocity)!{{{

      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         relativeVorticity
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdgeFlux
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalTransportVelocity
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalGMBolusVelocity

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         vertVelocityTop
      real (kind=RKIND), dimension(:,:), intent(out) :: &
         vertTransportVelocityTop
      real (kind=RKIND), dimension(:,:), intent(out) :: &
         vertGMBolusVelocityTop
      real (kind=RKIND), dimension(:,:), intent(out) :: &
         relativeVorticityCell
      real (kind=RKIND), dimension(:,:), intent(out) :: &
         divergence
      real (kind=RKIND), dimension(:,:), intent(out) :: &
         kineticEnergyCell
      real (kind=RKIND), dimension(:,:), intent(out) :: &
         tangentialVelocity

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: nCells, nEdges
      integer :: iCell, iVertex, iEdge
      integer :: i, j, k
      integer :: edgeSignOnCell_temp, eoe
      real(kind=RKIND) :: invAreaCell1
      real(kind=RKIND) :: dcEdge_temp, dvEdge_temp, r_tmp, weightsOnEdge_temp

      real (kind=RKIND), dimension(:), allocatable:: div_hu,div_huTransport,div_huGMBolus
#ifdef MPAS_OPENACC
      !$acc declare device_resident(div_hu, div_huTransport, div_huGMBolus)
#endif

      ! Need owned cells for relativeVorticityCell
      nCells = nCellsOwned

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc          present(relativeVorticityCell, areaCell, nEdgesOnCell, &
      !$acc                  kiteIndexOnCell, verticesOnCell, maxLevelCell, &
      !$acc                  relativeVorticity, kiteAreasOnVertex, invAreaCell, &
      !$acc                  minLevelCell) &
      !$acc          private(invAreaCell1, i, j, iVertex, k)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(invAreaCell1, i, j, k, iVertex)
#endif
      do iCell = 1, nCells
        relativeVorticityCell(:,iCell) = 0.0_RKIND
        invAreaCell1 = invAreaCell(iCell)

        do i = 1, nEdgesOnCell(iCell)
          j = kiteIndexOnCell(i, iCell)
          iVertex = verticesOnCell(i, iCell)
          do k = minLevelCell(iCell), maxLevelCell(iCell)
            relativeVorticityCell(k, iCell) = relativeVorticityCell(k, iCell) + kiteAreasOnVertex(j, iVertex) &
                                            * relativeVorticity(k, iVertex) * invAreaCell1
          end do
        end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      ! Need 0,1,2 halo cells
      nCells = nCellsHalo( 2 )

      !
      ! Compute divergence, kinetic energy, and vertical velocity
      !
      allocate(div_hu(nVertLevels),div_huTransport(nVertLevels),div_huGMBolus(nVertLevels))

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc          present(divergence, kineticEnergyCell, &
      !$acc                  areaCell, nEdgesOnCell, edgesOnCell, edgeSignOnCell, dcEdge, dvEdge, &
      !$acc                  maxLevelCell, normalVelocity, divergence, layerThicknessEdgeFlux, &
      !$acc                  normalTransportVelocity, normalGMBolusVelocity, vertVelocityTop, &
      !$acc                  vertTransportVelocityTop, vertGMBolusVelocityTop, invAreaCell, &
      !$acc                  minLevelCell) &
      !$acc          private(invAreaCell1, i, iEdge, edgeSignOnCell_temp, dcEdge_temp, &
      !$acc                  dvEdge_temp, k, r_tmp)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(i, iEdge, invAreaCell1, edgeSignOnCell_temp, dcEdge_temp, dvEdge_temp, &
      !$omp         k, r_tmp, div_hu, div_huTransport, div_huGMBolus)
#endif
      do iCell = 1, nCells
         divergence(:, iCell) = 0.0_RKIND
         kineticEnergyCell(:, iCell) = 0.0_RKIND
         div_hu(:) = 0.0_RKIND
         div_huTransport(:) = 0.0_RKIND
         div_huGMBolus(:) = 0.0_RKIND
         invAreaCell1 = invAreaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            edgeSignOnCell_temp = edgeSignOnCell(i, iCell)
            dcEdge_temp = dcEdge(iEdge)
            dvEdge_temp = dvEdge(iEdge)
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               r_tmp = dvEdge_temp * normalVelocity(k, iEdge) * invAreaCell1

               divergence(k, iCell) = divergence(k, iCell) - edgeSignOnCell_temp * r_tmp
               div_hu(k) = div_hu(k) - layerThicknessEdgeFlux(k, iEdge) * edgeSignOnCell_temp * r_tmp
               kineticEnergyCell(k, iCell) = kineticEnergyCell(k, iCell) &
                                              + 0.25 * r_tmp * dcEdge_temp * normalVelocity(k,iEdge)
               ! Compute vertical velocity from the horizontal total transport
               div_huTransport(k) = div_huTransport(k) &
                                    - layerThicknessEdgeFlux(k, iEdge) * edgeSignOnCell_temp &
                                    * dvEdge_temp * normalTransportVelocity(k, iEdge) * invAreaCell1
               ! Compute vertical velocity from the horizontal GM Bolus velocity
               div_huGMBolus(k) = div_huGMBolus(k) &
                                  - layerThicknessEdgeFlux(k, iEdge) * edgeSignOnCell_temp * dvEdge_temp &
                                  * normalGMBolusVelocity(k, iEdge) * invAreaCell1
            end do
         end do
         ! Vertical velocity at bottom (maxLevelCell(iCell)+1) is zero, initialized above.
         vertVelocityTop(1:minLevelCell(iCell)-1, iCell) = 0.0_RKIND
         vertVelocityTop(maxLevelCell(iCell)+1, iCell) = 0.0_RKIND
         do k = maxLevelCell(iCell), 1, -1
            vertVelocityTop(k,iCell) = vertVelocityTop(k+1,iCell) - div_hu(k)
            vertTransportVelocityTop(k,iCell) = vertTransportVelocityTop(k+1,iCell) - div_huTransport(k)
            vertGMBolusVelocityTop(k,iCell) = vertGMBolusVelocityTop(k+1,iCell) - div_huGMBolus(k)
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      deallocate(div_hu,div_huTransport,div_huGMBolus)

      nEdges = nEdgesHalo( 2 )
#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector &
      !$acc          present(tangentialVelocity, nEdgesOnEdge, edgesOnEdge, weightsOnEdge, &
      !$acc                  maxLevelEdgeTop, normalVelocity, minLevelEdgeBot) &
      !$acc          private(i, eoe, weightsOnEdge_temp, k)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(i, eoe, weightsOnEdge_temp, k)
#endif
      do iEdge = 1, nEdges
         tangentialVelocity(:, iEdge) = 0.0_RKIND
         ! Compute v (tangential) velocities
         do i = 1, nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            weightsOnEdge_temp = weightsOnEdge(i, iEdge)
            do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
               tangentialVelocity(k,iEdge) = tangentialVelocity(k,iEdge) &
                                              + weightsOnEdge_temp * normalVelocity(k, eoe)
            end do
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

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

   end subroutine ocn_diagnostic_solve_z_coordinates!}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve_pressure
!
!> \brief   Computes diagnostic variables for pressure or montgomery potential
!> \author  Matt Turner, Xylar Asay-Davis
!> \date    Jan. 2021
!> \details
!>  This routine computes the diagnostic variables for pressure or montgomery potential
!
!-----------------------------------------------------------------------
   subroutine ocn_diagnostic_solve_pressure(layerThickness, density, &
                surfacePressure, montgomeryPotential, pressure)!{{{

      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         density
      real (kind=RKIND), dimension(:), intent(in) :: &
         surfacePressure

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         montgomeryPotential
      real (kind=RKIND), dimension(:,:), intent(out) :: &
         pressure


      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: nCells
      integer :: iCell
      integer :: k

      real (kind=RKIND), dimension(:), allocatable:: pTop
#ifdef MPAS_OPENACC
      !$acc declare device_resident(pTop)
#endif

      !
      ! Pressure
      ! This section must be placed in the code after computing the density.
      !

      nCells = nCellsAll
      if (config_pressure_gradient_type.eq.'MontgomeryPotential') then

        ! use Montgomery Potential when layers are isopycnal.
        ! However, one may use 'pressure_and_zmid' when layers are isopycnal as well.
        ! Compute pressure at top of each layer, and then Montgomery Potential.

        allocate(pTop(nVertLevels))

#ifdef MPAS_OPENACC
        !$acc parallel loop present(minLevelCell, bottomdepth, montgomeryPotential, &
        !$acc                       density, layerThickness)
#else
        !$omp parallel
        !$omp do schedule(runtime) private(pTop, k)
#endif
        do iCell = 1, nCells

           ! assume atmospheric pressure at the surface is zero for now.
           pTop(1) = 0.0_RKIND
           ! At top layer it is g*SSH, where SSH may be off by a
           ! constant (ie, bottomDepth can be relative to top or bottom)
           montgomeryPotential(minLevelCell(iCell),iCell) = gravity &
              * (bottomDepth(iCell) + sum(layerThickness(1:nVertLevels,iCell)))

           do k = 2, nVertLevels
              pTop(k) = pTop(k-1) + density(k-1,iCell)*gravity* layerThickness(k-1,iCell)

              ! from delta M = p delta / density
              montgomeryPotential(k,iCell) = montgomeryPotential(k-1,iCell) &
                 + pTop(k)*(1.0_RKIND/density(k,iCell) - 1.0_RKIND/density(k-1,iCell))
           end do

        end do
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif

        deallocate(pTop)

      else

#ifdef MPAS_OPENACC
        !$acc parallel loop present(layerThickness, pressure, surfacePressure, density, &
        !$acc                       minLevelCell, maxLevelCell)
#else
        !$omp parallel
        !$omp do schedule(runtime) private(k)
#endif
        do iCell = 1, nCells
           ! Pressure for generalized coordinates.
           pressure(minLevelCell(iCell),iCell) = surfacePressure(iCell) &
             + density(minLevelCell(iCell),iCell)*gravity*0.5_RKIND*layerThickness(minLevelCell(iCell),iCell)

           do k = minLevelCell(iCell)+1, maxLevelCell(iCell)
              pressure(k,iCell) = pressure(k-1,iCell)  &
                + 0.5_RKIND*gravity*(  density(k-1,iCell)*layerThickness(k-1,iCell) &
                                     + density(k  ,iCell)*layerThickness(k  ,iCell))
           end do

        end do
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif

      endif

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

         gradSSH(iEdge) = edgeMask(1, iEdge) * ( pressureAdjustedSSH(cell2) - pressureAdjustedSSH(cell1) ) / dcEdge(iEdge)
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
   subroutine ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, oldLayerThickness, layerThicknessEdgeFlux, &
     normalVelocity, oldSSH, dt, vertAleTransportTop, err, newHighFreqThickness)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool           !< Input: horizonal mesh information

      type (mpas_pool_type), intent(in) :: &
         verticalMeshPool   !< Input: vertical mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         oldLayerThickness    !< Input: layer thickness at old time

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdgeFlux     !< Input: layerThickness interpolated to an edge

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity     !< Input: transport

      real (kind=RKIND), dimension(:), intent(in) :: &
         oldSSH     !< Input: sea surface height at old time

      real (kind=RKIND), dimension(:,:), intent(in), optional :: &
         newHighFreqThickness   !< Input: high frequency thickness.  Alters ALE thickness.

      real (kind=RKIND), intent(in) :: &
         dt     !< Input: time step

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         vertAleTransportTop     !< Output: vertical transport at top of cell

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, iCell, k, i, nCells, kmin, kmax

      real (kind=RKIND) :: flux, invAreaCell1, div_hu_btr
      integer, dimension(:,:), pointer :: edgeSignOnCell

      ! Scratch Arrays
      ! projectedSSH: projected SSH at a new time
      !        units: none
      real (kind=RKIND), dimension(:), allocatable ::  projectedSSH
      ! div_hu: divergence of (thickness*velocity)
      !  units: none
      real (kind=RKIND), dimension(:,:), allocatable :: div_hu
      ! ALE_Thickness: ALE thickness at new time
      !         units: none
      real (kind=RKIND), dimension(:,:), allocatable :: ALE_Thickness
#ifdef MPAS_OPENACC
      !$acc declare device_resident(div_hu)
#endif

      err = 0

      if (config_vert_coord_movement == 'impermeable_interfaces' .or. &
          configVertAdvMethod == vertAdvRemap) then
        vertAleTransportTop=0.0_RKIND
        !$acc update device(vertAleTransportTop)
        return
      end if

      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

      ncells = nCellsAll

      allocate(div_hu(nVertLevels, nCells), &
               projectedSSH(nCells), &
               ALE_Thickness(nVertLevels, nCells))

      ! Only need to compute over 0 and 1 halos
      nCells = nCellsHalo( 1 )

      !
      ! thickness-weighted divergence and barotropic divergence
      !
      ! See Ringler et al. (2010) jcp paper, eqn 19, 21, and fig. 3.
#ifdef MPAS_OPENACC
      !$acc enter data create(projectedSSH, ALE_Thickness) &
      !$acc            copyin(edgeSignOnCell)

      !$acc parallel loop &
      !$acc          present(invAreaCell, nEdgesOnCell, edgesOnCell, maxLevelEdgeTop, &
      !$acc                  layerThicknessEdgeFlux, normalVelocity, dvEdge, edgeSignOnCell, &
      !$acc                  oldSSH, projectedSSH, minLevelEdgeBot) &
      !$acc          private(div_hu_btr, invAreaCell1, i, iEdge, flux, k, kmin, kmax)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(invAreaCell1, i, iEdge, k, flux, div_hu_btr, kmin, kmax)
#endif
      do iCell = 1, nCells
         div_hu(:,iCell) = 0.0_RKIND
         div_hu_btr = 0.0_RKIND
         invAreaCell1 = invAreaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            kmin = minLevelEdgeBot(iEdge)
            kmax = maxLevelEdgeTop(iEdge)

            do k = kmin, kmax
               flux = layerThicknessEdgeFlux(k, iEdge) * normalVelocity(k, iEdge) * dvEdge(iEdge) * edgeSignOnCell(i, iCell) &
                    * invAreaCell1
               div_hu(k,iCell) = div_hu(k,iCell) - flux
               div_hu_btr = div_hu_btr - flux
            end do
         end do
         projectedSSH(iCell) = oldSSH(iCell) - dt*div_hu_btr
      end do
#ifdef MPAS_OPENACC
      !$acc update host(projectedSSH)
      !$acc exit data delete(edgeSignOnCell)
#else
      !$omp end do
      !$omp end parallel
#endif

      !
      ! Compute desired thickness at new time
      !
      if (present(newHighFreqThickness)) then
        call ocn_ALE_thickness(meshPool, verticalMeshPool, projectedSSH, ALE_thickness, err, newHighFreqThickness)
      else
        call ocn_ALE_thickness(meshPool, verticalMeshPool, projectedSSH, ALE_thickness, err)
      endif

      !
      ! Vertical transport through layer interfaces
      !
      ! Vertical transport through layer interface at top and bottom is zero.
      ! Here we are using solving the continuity equation for vertAleTransportTop ($w^t$),
      ! and using ALE_Thickness for thickness at the new time.

#ifdef MPAS_OPENACC
      !$acc update device(ALE_Thickness)

      !$acc parallel loop gang vector &
      !$acc          present(vertAleTransportTop, maxLevelCell, ALE_Thickness, &
      !$acc                  oldLayerThickness, minLevelCell) &
      !$acc          private(k)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iCell = 1,nCells
         vertAleTransportTop(1,iCell) = 0.0_RKIND
         vertAleTransportTop(maxLevelCell(iCell)+1,iCell) = 0.0_RKIND
         do k = maxLevelCell(iCell), minLevelCell(iCell)+1, -1
            vertAleTransportTop(k,iCell) = vertAleTransportTop(k+1,iCell) - div_hu(k,iCell) &
              - (ALE_Thickness(k,iCell) - oldLayerThickness(k,iCell))/dt
         end do
      end do
#ifdef MPAS_OPENACC
      !$acc exit data delete(projectedSSH, ALE_Thickness)
#else
      !$omp end do
      !$omp end parallel
#endif

      deallocate(div_hu, &
                 projectedSSH, &
                 ALE_Thickness)

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
!  routine ocn_compute_KPP_input_fields
!
!> \brief
!>    Compute fields necessary to drive the CVMix KPP module
!> \author  Todd Ringler
!> \date    20 August 2013
!> \details
!>    CVMix/KPP requires the following fields as intent(in):
!>       surfaceBuoyancyForcing
!>       surfaceFrictionVelocity
!>
!
!-----------------------------------------------------------------------

    subroutine ocn_compute_KPP_input_fields(statePool, forcingPool, meshPool, timeLevelIn)!{{{

      type (mpas_pool_type), intent(in) :: statePool !< Input/Output: State information
      type (mpas_pool_type), intent(in) :: forcingPool !< Input: Forcing information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      integer, intent(in), optional :: timeLevelIn

      ! pool pointers
      type (mpas_pool_type), pointer :: tracersSurfaceFluxPool
      type (mpas_pool_type), pointer :: tracersPool

      ! scalars
      integer :: nCells

      ! real pointers
      real (kind=RKIND), dimension(:), pointer :: penetrativeTemperatureFlux, surfaceThicknessFlux, &
           surfaceThicknessFluxRunoff, rainTemperatureFlux, evapTemperatureFlux, &
           icebergTemperatureFlux, seaIceTemperatureFlux

      real (kind=RKIND), dimension(:), pointer :: surfaceStress, surfaceStressMagnitude
      real (kind=RKIND), dimension(:,:), pointer ::  &
           layerThickness, &
           normalVelocity, activeTracersSurfaceFlux, &
           activeTracersSurfaceFluxRunoff, nonLocalSurfaceTracerFlux

      real (kind=RKIND), dimension(:,:,:), pointer :: &
           activeTracers

      ! local
      integer :: iCell, iEdge, i, err, timeLevel
      integer, pointer :: indexTempFlux, indexSaltFlux
      real (kind=RKIND) :: turbulentVelocitySquared, fracAbsorbed, fracAbsorbedRunoff
      real (kind=RKIND) :: sumSurfaceStressSquared

      ! Scratch Arrays
      ! densitySurfaceDisplaced: Density computed by displacing SST and SSS to every vertical
      !                          layer within the column
      !                   units: kg m^{-3}
      real (kind=RKIND), dimension(:,:), allocatable :: densitySurfaceDisplaced
#ifdef MPAS_OPENACC
      !$acc declare device_resident(densitySurfaceDisplaced)
#endif
      ! thermalExpansionCoeff: Thermal expansion coefficient (alpha), defined as $-1/\rho
      !                        d\rho/dT$ (note negative sign)
      !                 units: C^{-1}
      real (kind=RKIND), dimension(:,:), allocatable :: thermalExpansionCoeff
      ! salineContractionCoeff: Saline contraction coefficient (beta), defined as $1/\rho
      !                         d\rho/dS$.  This is also called the haline contraction coefficient.
      !                  units: PSU^{-1}
      real (kind=RKIND), dimension(:,:), allocatable :: salineContractionCoeff

      call mpas_timer_start('KPP input fields')

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux', tracersSurfaceFluxPool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
      ! set the parameter turbulentVelocitySquared
      turbulentVelocitySquared = 0.001_RKIND

      ! set scalar values
      call mpas_pool_get_dimension(tracersSurfaceFluxPool, 'index_temperatureSurfaceFlux', indexTempFlux)
      call mpas_pool_get_dimension(tracersSurfaceFluxPool, 'index_salinitySurfaceFlux', indexSaltFlux)

      ! set pointers into state, mesh, diagnostics and scratch
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)

      call mpas_pool_get_array(tracersSurfaceFluxPool, 'nonLocalSurfaceTracerFlux', nonLocalSurfaceTracerFlux)

      call mpas_pool_get_array(forcingPool, 'surfaceThicknessFlux', surfaceThicknessFlux)
      call mpas_pool_get_array(forcingPool, 'surfaceThicknessFluxRunoff', surfaceThicknessFluxRunoff)
      call mpas_pool_get_array(forcingPool, 'penetrativeTemperatureFlux', penetrativeTemperatureFlux)
      call mpas_pool_get_array(forcingPool, 'surfaceStress', surfaceStress)
      call mpas_pool_get_array(forcingPool, 'surfaceStressMagnitude', surfaceStressMagnitude)
      call mpas_pool_get_array(forcingPool, 'rainTemperatureFlux', rainTemperatureFlux)
      call mpas_pool_get_array(forcingPool, 'evapTemperatureFlux', evapTemperatureFlux)
      call mpas_pool_get_array(forcingPool, 'seaIceTemperatureFlux', seaIceTemperatureFlux)
      call mpas_pool_get_array(forcingPool, 'icebergTemperatureFlux', icebergTemperatureFlux)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
      call mpas_pool_get_array(tracersSurfaceFluxPool, 'activeTracersSurfaceFlux', activeTracersSurfaceFlux)
      call mpas_pool_get_array(tracersSurfaceFluxPool, 'activeTracersSurfaceFluxRunoff', activeTracersSurfaceFluxRunoff)

      ! allocate scratch space displaced density computation
      ncells = nCellsAll

      allocate(densitySurfaceDisplaced(nVertLevels, nCells), &
               thermalExpansionCoeff(nVertLevels, nCells), &
               salineContractionCoeff(nVertLevels, nCells))

      ! Only need to compute over the 0 and 1 halos
      nCells = nCellsHalo( 2 )

#ifdef MPAS_OPENACC
      !$acc enter data create(thermalExpansionCoeff,   &
      !$acc                   salineContractionCoeff)
#endif

      ! compute EOS by displacing SST/SSS to every vertical layer in column

      ! inputs: activeTracers, tracersSurfaceValue,maxLevelCell,refBottomDepth 
      ! outputs: densitySurfaceDisplayed, thermalExpansionCoeff, salineContractionCoeff
      call ocn_equation_of_state_density(statePool, meshPool, activeTracers, tracersSurfaceValue, &
                                         nCells, 0, 'surfaceDisplaced', densitySurfaceDisplaced, &
                                         err, thermalExpansionCoeff, salineContractionCoeff, &
                                         timeLevel)

#ifdef MPAS_OPENACC
      !$acc exit data copyout(thermalExpansionCoeff,   &
      !$acc                   salineContractionCoeff)
#endif

      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(fracAbsorbed, fracAbsorbedRunoff, sumSurfaceStressSquared, i, iEdge)
      do iCell = 1, nCells
         ! compute surface buoyancy forcing based on surface fluxes of mass, temperature, salinity and frazil
         !                                (frazil to be added later)
         ! since this computation is confusing, variables, units and sign convention is repeated here
         ! everything below should be consistent with that specified in Registry
         ! everything below should be consistent with the CVMix/KPP documentation:
         !           https://www.dropbox.com/s/6hqgc0rsoa828nf/cvmix_20aug2013.pdf
         !
         !    surfaceThicknessFlux: surface mass flux, m/s, positive into ocean
         !    activeTracersSurfaceFlux(indexTempFlux): non-penetrative temperature flux, C m/s, positive into ocean
         !    penetrativeTemperatureFlux: penetrative surface temperature flux at ocean surface, positive into ocean
         !    activeTracersSurfaceFlux(indexSaltFlux): salinity flux, PSU m/s, positive into ocean
         !    penetrativeTemperatureFluxOBL: penetrative temperature flux computed at z=OBL, positive down
         !
         ! note: the following fields used the CVMix/KPP computation of buoyancy forcing are not included here
         !    1. Tm: temperature associated with surfaceThicknessFlux, C  (here we assume Tm == temperatureSurfaceValue)
         !    2. Sm: salinity associated with surfaceThicknessFlux, PSU (here we assume Sm == salinitySurfaceValue and account for
         !           salinity flux in activeTracersSurfaceFlux array)
         !

         ! Compute fraction of thickness flux that is in the top model layer
         fracAbsorbed = 1.0_RKIND - exp( max(-layerThickness(minLevelCell(iCell), iCell) / config_flux_attenuation_coefficient, -100.0_RKIND) )
         fracAbsorbedRunoff = 1.0_RKIND - exp( max(-layerThickness(minLevelCell(iCell), iCell) / config_flux_attenuation_coefficient_runoff, &
                              -100.0_RKIND) )
         ! Store the total tracer flux below in nonLocalSurfaceTemperatureFlux for use in the CVMix nonlocal
         ! transport code.  This includes tracer forcing due to thickness

         nonLocalSurfaceTracerFlux(indexTempFlux, iCell) = activeTracersSurfaceFlux(indexTempFlux,iCell) &
                 + penetrativeTemperatureFlux(iCell) - penetrativeTemperatureFluxOBL(iCell)  &
                 - fracAbsorbed * (rainTemperatureFlux(iCell) + evapTemperatureFlux(iCell) + &
                                   seaIceTemperatureFlux(iCell) + icebergTemperatureFlux(iCell)) &
                 - fracAbsorbedRunoff * activeTracersSurfaceFluxRunoff(indexTempFlux, iCell)

         nonLocalSurfaceTracerFlux(indexSaltFlux,iCell) = activeTracersSurfaceFlux(indexSaltFlux,iCell) &
                 - fracAbsorbed * surfaceThicknessFlux(iCell) * activeTracers(indexSaltFlux,minLevelCell(iCell),iCell) &
                 - fracAbsorbedRunoff * surfaceThicknessFluxRunoff(iCell) * activeTracers(indexSaltFlux,minLevelCell(iCell),iCell)

         surfaceBuoyancyForcing(iCell) =  thermalExpansionCoeff (minLevelCell(iCell),iCell) &
                * nonLocalSurfaceTracerFlux(indexTempFlux,iCell) &
                - salineContractionCoeff(minLevelCell(iCell),iCell) * nonLocalSurfaceTracerFlux(indexSaltFlux,iCell)

         ! at this point, surfaceBuoyancyForcing has units of m/s
         ! change into units of m^2/s^3 (which can be thought of as the flux of buoyancy, units of buoyancy * velocity )
         surfaceBuoyancyForcing(iCell) = surfaceBuoyancyForcing(iCell) * gravity

         ! compute magnitude of surface stress
         sumSurfaceStressSquared = 0.0_RKIND
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            sumSurfaceStressSquared = sumSurfaceStressSquared + edgeAreaFractionOfCell(i,iCell) * surfaceStress(iEdge)**2
         enddo
         ! NOTE that the factor of 2 is from averaging dot products to cell centers on a C-grid
         surfaceStressMagnitude(iCell) = sqrt(2.0_RKIND * sumSurfaceStressSquared)
         surfaceFrictionVelocity(iCell) = sqrt(surfaceStressMagnitude(iCell) / rho_sw)

      enddo
      !$omp end do
      !$omp end parallel

#ifdef MPAS_OPENACC
      !$acc exit data delete(thermalExpansionCoeff,   &
      !$acc                  salineContractionCoeff)
#endif

      ! deallocate scratch space
      deallocate(thermalExpansionCoeff, &
                 salineContractionCoeff, &
                 densitySurfaceDisplaced)

      call mpas_timer_stop('KPP input fields')

    end subroutine ocn_compute_KPP_input_fields!}}}


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
      activeTracers, landIceFraction, landIceMask, &
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
         landIceMask

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, iEdge, cell1, cell2, iLevel, i

      real (kind=RKIND) :: blThickness, dz, weightSum, h_nu, Gamma_turb, landIceEdgeFraction, velocityMagnitude

      ! Scratch Arrays
      real (kind=RKIND), dimension(:), allocatable ::  blTempScratch, blSaltScratch
#ifdef MPAS_OPENACC
      !$acc declare device_resident(blTempScratch, blSaltScratch)
#endif

      ! constants for Holland and Jenkins 1999 parameterization of the boundary layer
      real (kind=RKIND), parameter :: &
         Pr = 13.8_RKIND, &             ! the Prandtl number
         Sc = 2432.0_RKIND, &           ! the Schmidt number
         nuSaltWater = 1.95e-6_RKIND, & ! molecular viscosity of sea water (m^2/s)
         kVonKarman = 0.4_RKIND, &      ! the von Karman constant
         xiN = 0.052_RKIND              ! dimensionless planetary boundary layer constant

      integer :: indexBLTval, indexBLSval, indexHeatTransval, indexSaltTransval

      if ( trim(config_land_ice_flux_mode) .ne. 'standalone' .and. trim(config_land_ice_flux_mode) .ne. 'coupled' ) then
         return
      end if

      call mpas_timer_start("land_ice_diagnostic_fields", .false.)

      indexBLTval = indexBLT
      indexBLSval = indexBLS
      indexHeatTransval = indexHeatTrans
      indexSaltTransval = indexSaltTrans

      allocate(blTempScratch(nCellsAll), &
               blSaltScratch(nCellsAll))

      ! Compute top drag
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

         topDrag(iEdge) = - rho_sw * landIceEdgeFraction * config_land_ice_flux_topDragCoeff &
                          * velocityMagnitude * normalVelocity(minLevelEdgeBot(iEdge),iEdge)

      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

      ! compute top drag magnitude and friction velocity at cell centers
#ifdef MPAS_OPENACC
      !$acc parallel loop present(landIceFrictionVelocity, kineticEnergyCell, minLevelCell, topDragMag) &
      !$acc    present(landIceFraction)
#else
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCellsAll
         ! the magnitude of the top drag is CD*u**2 = CD*(2*KE)
         topDragMag(iCell) = rho_sw * landIceFraction(iCell) &
                           * 2.0_RKIND * config_land_ice_flux_topDragCoeff *  kineticEnergyCell(minLevelCell(iCell),iCell)

         ! the friction velocity is the square root of the top drag + variance of tidal velocity
         ! (computed regardless of land-ice coverage)
         landIceFrictionVelocity(iCell) = sqrt(config_land_ice_flux_topDragCoeff * &
                                   (2.0_RKIND * kineticEnergyCell(minLevelCell(iCell),iCell) &
                                   + config_land_ice_flux_rms_tidal_velocity**2))
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

      ! average temperature and salinity over horizontal neighbors and the sub-ice-shelf boundary layer
#ifdef MPAS_OPENACC
      !$acc parallel loop present(minLevelCell, maxLevelCell, layerThickness, activeTracers)
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
      !$acc parallel loop present(cellsOnCell, nEdgesOnCell, minLevelCell) &
      !$acc    present(landIceBoundaryLayerTracers, cellMask)
#else
      !$omp do schedule(runtime) private(weightSum, i, cell2)
#endif
      do iCell = 1, nCellsAll
         landIceBoundaryLayerTracers(indexBLTval, iCell) = blTempScratch(iCell)
         landIceBoundaryLayerTracers(indexBLSval, iCell) = blSaltScratch(iCell)
         if(config_land_ice_flux_boundaryLayerNeighborWeight > 0.0_RKIND) then
            weightSum = 1.0_RKIND
            do i = 1, nEdgesOnCell(iCell)
               cell2 = cellsOnCell(i,iCell)

               landIceBoundaryLayerTracers(indexBLTval, iCell) = landIceBoundaryLayerTracers(indexBLTval, iCell) &
                 + cellMask(minLevelCell(cell2),cell2)*config_land_ice_flux_boundaryLayerNeighborWeight*blTempScratch(cell2)
               landIceBoundaryLayerTracers(indexBLSval, iCell) = landIceBoundaryLayerTracers(indexBLSval, iCell) &
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
            landIceTracerTransferVelocities(indexHeatTransval, iCell) = landIceFrictionVelocity(iCell) &
                                             * config_land_ice_flux_jenkins_heat_transfer_coefficient
            landIceTracerTransferVelocities(indexSaltTransval, iCell) = landIceFrictionVelocity(iCell) &
                                             * config_land_ice_flux_jenkins_salt_transfer_coefficient
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
      else if(hollandJenkinsOn) then
#ifdef MPAS_OPENACC
         !$acc parallel loop present(fCell, landIceTracerTransferVelocities, landIceFrictionVelocity)
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

            landIceTracerTransferVelocities(indexHeatTransval, iCell) = 1.0_RKIND/(Gamma_turb + 12.5_RKIND &
                                                                   * Pr**(2.0_RKIND/3.0_RKIND) - 6.0_RKIND)
            landIceTracerTransferVelocities(indexSaltTransval, iCell) = 1.0_RKIND/(Gamma_turb + 12.5_RKIND &
                                                                   * Sc**(2.0_RKIND/3.0_RKIND) - 6.0_RKIND)
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
      end if

      deallocate(blTempScratch, &
                 blSaltScratch)

      ! modify the spatially-varying attenuation coefficient where there is land ice
#ifdef MPAS_OPENACC
      !$acc enter data copyin(landIceMask)

      !$acc parallel loop present(landIceMask, sfcFlxAttCoeff)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCellsAll
         if(landIceMask(iCell) == 1) then
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
!  routine ocn_reconstruct_gm_vectors
!
!> \brief   Computes cell-centered vector diagnostics
!> \author  Mark Petersen
!> \date    May 2014
!> \details
!>  This routine computes cell-centered vector diagnostics
!
!-----------------------------------------------------------------------

   subroutine ocn_reconstruct_gm_vectors(meshPool) !{{{

      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information

      call mpas_timer_start('reconstruct gm vecs')

      call mpas_reconstruct(meshPool, normalTransportVelocity,          &
                       transportVelocityX,            &
                       transportVelocityY,            &
                       transportVelocityZ,            &
                       transportVelocityZonal,        &
                       transportVelocityMeridional    &
                      )

      call mpas_reconstruct(meshPool, normalGMBolusVelocity,          &
                       GMBolusVelocityX,            &
                       GMBolusVelocityY,            &
                       GMBolusVelocityZ,            &
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

      call mpas_timer_stop('reconstruct gm vecs')

   end subroutine ocn_reconstruct_gm_vectors!}}}


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

      if (trim(config_time_integrator) == 'RK4') then
         ! For RK4, PV includes f: PV = (eta+f)/h.
         fCoef = 1
      elseif (trim(config_time_integrator) == 'split_explicit' &
        .or.trim(config_time_integrator) == 'unsplit_explicit' &
        .or.trim(config_time_integrator) == 'split_implicit') then
          ! For split explicit, PV is eta/h because the Coriolis term
          ! is added separately to the momentum tendencies.
          fCoef = 0
      end if

      ! Initialize choice for computing thickness at edges for fluxes
      if (trim(config_thickness_flux_type) == 'centered') then
         thickEdgeFluxChoice = thickEdgeFluxCenter
      else
         if (config_use_wetting_drying .and. &
             (trim(config_thickness_flux_type) /= 'centered')) then
            if ( trim(config_thickness_flux_type) == 'upwind') then
               thickEdgeFluxChoice = thickEdgeFluxUpwind
            end if
         else
            thickEdgeFluxChoice = thickEdgeFluxUnknown
            call mpas_log_write('Thickness flux option of ' //&
                 & trim(config_thickness_flux_type) // &
                 & 'is not known', MPAS_LOG_CRIT)
         end if
      end if

      call ocn_diagnostics_variables_init(domain, jenkinsOn, &
                                          hollandJenkinsOn, err)

    end subroutine ocn_diagnostics_init!}}}

!***********************************************************************

end module ocn_diagnostics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
