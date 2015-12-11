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

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_timer
   use mpas_threading
   use mpas_vector_reconstruction

   use ocn_constants
   use ocn_gm
   use ocn_equation_of_state
   use ocn_thick_ale
   use ocn_diagnostics_routines

   implicit none
   private
   save

   type (timer_node), pointer :: diagEOSTimer

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
             ocn_vert_transport_velocity_top, &
             ocn_fuperp, &
             ocn_filter_btr_mode_vel, &
             ocn_filter_btr_mode_tend_vel, &
             ocn_reconstruct_gm_vectors, &
             ocn_diagnostics_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: ke_cell_flag, ke_vertex_flag
   real (kind=RKIND) ::  fCoef
   real (kind=RKIND), pointer ::  coef_3rd_order

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

   subroutine ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, timeLevelIn)!{{{

      real (kind=RKIND), intent(in) :: dt !< Input: Time step
      type (mpas_pool_type), intent(in) :: statePool !< Input: State information
      type (mpas_pool_type), intent(in) :: forcingPool !< Input: Forcing information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      type (mpas_pool_type), intent(inout) :: diagnosticsPool  !< Input: diagnostic fields derived from State
      type (mpas_pool_type), intent(in) :: scratchPool !< Input: scratch variables
      type (mpas_pool_type), intent(in) :: tracersPool !< Input: tracer fields
      integer, intent(in), optional :: timeLevelIn !< Input: Time level in state

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j
      integer :: boundaryMask, velMask, err
      integer, pointer  :: nEdgesSolve, nCells, nEdges, nVertices, nVertLevels, vertexDegree

      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, &
        maxLevelCell, maxLevelEdgeTop, maxLevelEdgeBot, &
        maxLevelVertexBot
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, &
        verticesOnEdge, edgesOnEdge, edgesOnVertex,boundaryCell, kiteIndexOnCell, &
        verticesOnCell, edgeSignOnVertex, edgeSignOnCell, edgesOnCell

      real (kind=RKIND) :: d2fdx2_cell1, d2fdx2_cell2, coef_3rd_order, r_tmp, &
        invAreaCell1, invAreaCell2, invAreaTri1, invAreaTri2, invLength, layerThicknessVertex, coef, &
        shearMean, shearSquared, factor, delU2, sumSurfaceLayer, surfaceLayerDepth, rSurfaceLayer

      real (kind=RKIND), dimension(:), allocatable:: pTop, div_hu,div_huTransport,div_huGMBolus

      real (kind=RKIND), dimension(:), pointer :: &
        bottomDepth, fVertex, dvEdge, dcEdge, areaCell, areaTriangle, ssh, seaSurfacePressure, frazilSurfacePressure
      real (kind=RKIND), dimension(:,:), pointer :: &
        weightsOnEdge, kiteAreasOnVertex, layerThicknessEdge, layerThickness, normalVelocity, normalTransportVelocity, normalGMBolusVelocity, tangentialVelocity, pressure,&
        circulation, kineticEnergyCell, montgomeryPotential, vertAleTransportTop, zMid, zTop, divergence, &
        relativeVorticity, relativeVorticityCell, &
        normalizedPlanetaryVorticityEdge, normalizedPlanetaryVorticityVertex, &
        normalizedRelativeVorticityEdge, normalizedRelativeVorticityVertex, normalizedRelativeVorticityCell, &
        density, displacedDensity, potentialDensity, temperature, salinity, kineticEnergyVertex, kineticEnergyVertexOnCells, &
        vertVelocityTop, vertTransportVelocityTop, vertGMBolusVelocityTop, BruntVaisalaFreqTop, &
        vorticityGradientNormalComponent, vorticityGradientTangentialComponent, gradSSH, RiTopOfCell, &
        inSituThermalExpansionCoeff, inSituSalineContractionCoeff

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers, derivTwo
      character :: c1*6

      real (kind=RKIND), dimension(:,:), pointer :: tracersSurfaceValue
      real (kind=RKIND), dimension(:,:), pointer :: tracersSurfaceLayerValue
      real (kind=RKIND), dimension(:),   pointer :: normalVelocitySurfaceLayer
      real (kind=RKIND), dimension(:),   pointer :: indexSurfaceLayerDepth

      type (field2DReal), pointer :: kineticEnergyVertexField, kineticEnergyVertexOnCellsField
      type (field2DReal), pointer :: normalizedRelativeVorticityVertexField, normalizedPlanetaryVorticityVertexField
      type (field2DReal), pointer :: vorticityGradientNormalComponentField, vorticityGradientTangentialComponentField

      integer :: timeLevel
      integer, pointer :: indexTemperature, indexSalinity
      logical, pointer :: config_use_cvmix_kpp
      real (kind=RKIND), pointer :: config_apvm_scale_factor,  config_coef_3rd_order, config_cvmix_kpp_surface_layer_averaging
      character (len=StrKIND), pointer :: config_pressure_gradient_type

      real (kind=RKIND), pointer :: config_flux_attenuation_coefficient
      real (kind=RKIND), dimension(:), pointer :: surfaceFluxAttenuationCoefficient

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_config(ocnConfigs, 'config_apvm_scale_factor', config_apvm_scale_factor)
      call mpas_pool_get_config(ocnConfigs, 'config_pressure_gradient_type', config_pressure_gradient_type)
      call mpas_pool_get_config(ocnConfigs, 'config_coef_3rd_order', config_coef_3rd_order)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_surface_layer_averaging', config_cvmix_kpp_surface_layer_averaging)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)
      call mpas_pool_get_config(ocnConfigs, 'config_flux_attenuation_coefficient', config_flux_attenuation_coefficient)

      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
      call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)

      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
      call mpas_pool_get_array(diagnosticsPool, 'zTop', zTop)
      call mpas_pool_get_array(diagnosticsPool, 'divergence', divergence)
      call mpas_pool_get_array(diagnosticsPool, 'circulation', circulation)
      call mpas_pool_get_array(diagnosticsPool, 'relativeVorticity', relativeVorticity)
      call mpas_pool_get_array(diagnosticsPool, 'relativeVorticityCell', relativeVorticityCell)
      call mpas_pool_get_array(diagnosticsPool, 'normalizedPlanetaryVorticityEdge', normalizedPlanetaryVorticityEdge)
      call mpas_pool_get_array(diagnosticsPool, 'normalizedRelativeVorticityEdge', normalizedRelativeVorticityEdge)
      call mpas_pool_get_array(diagnosticsPool, 'normalizedRelativeVorticityCell', normalizedRelativeVorticityCell)
      call mpas_pool_get_array(diagnosticsPool, 'density', density)
      call mpas_pool_get_array(diagnosticsPool, 'displacedDensity', displacedDensity)
      call mpas_pool_get_array(diagnosticsPool, 'potentialDensity', potentialDensity)
      call mpas_pool_get_array(diagnosticsPool, 'montgomeryPotential', montgomeryPotential)
      call mpas_pool_get_array(diagnosticsPool, 'pressure', pressure)
      call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop', BruntVaisalaFreqTop)
      call mpas_pool_get_array(diagnosticsPool, 'tangentialVelocity', tangentialVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCell)
      call mpas_pool_get_array(diagnosticsPool, 'vertVelocityTop', vertVelocityTop)
      call mpas_pool_get_array(diagnosticsPool, 'vertTransportVelocityTop', vertTransportVelocityTop)
      call mpas_pool_get_array(diagnosticsPool, 'vertGMBolusVelocityTop', vertGMBolusVelocityTop)
      call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'gradSSH', gradSSH)
      call mpas_pool_get_array(diagnosticsPool, 'RiTopOfCell', RiTopOfCell)

      call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnVertex', edgesOnVertex)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(meshPool, 'fVertex', fVertex)
      call mpas_pool_get_array(meshPool, 'derivTwo', derivTwo)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', maxLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)
      call mpas_pool_get_array(meshPool, 'kiteIndexOnCell', kiteIndexOnCell)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'boundaryCell', boundaryCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnVertex', edgeSignOnVertex)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

      call mpas_pool_get_array(forcingPool, 'seaSurfacePressure', seaSurfacePressure)
      call mpas_pool_get_array(forcingPool, 'frazilSurfacePressure', frazilSurfacePressure)
                  
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

      call mpas_pool_get_array(diagnosticsPool, 'tracersSurfaceValue', tracersSurfaceValue)
      call mpas_pool_get_array(diagnosticsPool, 'tracersSurfaceLayerValue', tracersSurfaceLayerValue)
      call mpas_pool_get_array(diagnosticsPool, 'normalVelocitySurfaceLayer', normalVelocitySurfaceLayer)
      call mpas_pool_get_array(diagnosticsPool, 'indexSurfaceLayerDepth', indexSurfaceLayerDepth)

      call mpas_pool_get_array(diagnosticsPool, 'surfaceFluxAttenuationCoefficient', surfaceFluxAttenuationCoefficient)
      !
      ! Compute height on cell edges at velocity locations
      !   Namelist options control the order of accuracy of the reconstructed layerThicknessEdge value
      !

      ! initialize layerThicknessEdge to avoid divide by zero and NaN problems.
      !$omp workshare
      layerThicknessEdge = -1.0e34
      !$omp end workshare
      coef_3rd_order = config_coef_3rd_order 

      !$omp do schedule(runtime) private(cell1, cell2, k)
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = 1, maxLevelEdgeTop(iEdge)
            layerThicknessEdge(k,iEdge) = 0.5 * (layerThickness(k,cell1) + layerThickness(k,cell2))
         end do
      end do
      !$omp end do

      !
      ! set the velocity and height at dummy address
      !    used -1e34 so error clearly occurs if these values are used.
      !

      !$omp workshare
      normalVelocity(:,nEdges+1) = -1e34
      layerThickness(:,nCells+1) = -1e34
      activeTracers(indexTemperature,:,nCells+1) = -1e34
      activeTracers(indexSalinity,:,nCells+1) = -1e34

      divergence(:,:) = 0.0
      vertVelocityTop(:,:)=0.0
      kineticEnergyCell(:,:) = 0.0
      tangentialVelocity(:,:) = 0.0
      !$omp end workshare

      call mpas_threading_barrier()

      call ocn_relativeVorticity_circulation(relativeVorticity, circulation, meshPool, normalVelocity, err)

      call mpas_threading_barrier()
      !$omp workshare
      relativeVorticityCell(:,:) = 0.0
      !$omp end workshare

      !$omp do schedule(runtime) private(invAreaCell1, i, j, k, iVertex)
      do iCell = 1, nCells
        invAreaCell1 = 1.0 / areaCell(iCell)

        do i = 1, nEdgesOnCell(iCell)
          j = kiteIndexOnCell(i, iCell)
          iVertex = verticesOnCell(i, iCell)
          do k = 1, maxLevelCell(iCell)
            relativeVorticityCell(k, iCell) = relativeVorticityCell(k, iCell) + kiteAreasOnVertex(j, iVertex) * relativeVorticity(k, iVertex) * invAreaCell1
          end do
        end do
      end do
      !$omp end do

      !
      ! Compute divergence, kinetic energy, and vertical velocity
      !
      allocate(div_hu(nVertLevels),div_huTransport(nVertLevels),div_huGMBolus(nVertLevels))

      !$omp do schedule(runtime) private(invAreaCell1, iEdge, r_tmp, i, k)
      do iCell = 1, nCells
         div_hu(:) = 0.0
         div_huTransport(:) = 0.0
         div_huGMBolus(:) = 0.0
         invAreaCell1 = 1.0 / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelCell(iCell)
               r_tmp = dvEdge(iEdge) * normalVelocity(k, iEdge) * invAreaCell1

               divergence(k, iCell) = divergence(k, iCell) - edgeSignOnCell(i, iCell) * r_tmp
               div_hu(k)    = div_hu(k) - layerThicknessEdge(k, iEdge) * edgeSignOnCell(i, iCell) * r_tmp 
               kineticEnergyCell(k, iCell) = kineticEnergyCell(k, iCell) + 0.25 * r_tmp * dcEdge(iEdge) * normalVelocity(k,iEdge)

               ! Compute vertical velocity from the horizontal total transport
               div_huTransport(k) = div_huTransport(k) - layerThicknessEdge(k, iEdge) * edgeSignOnCell(i, iCell) * dvEdge(iEdge) * normalTransportVelocity(k, iEdge) * invAreaCell1 
               ! Compute vertical velocity from the horizontal GM Bolus velocity
               div_huGMBolus(k)   = div_huGMBolus(k) - layerThicknessEdge(k, iEdge) * edgeSignOnCell(i, iCell) * dvEdge(iEdge) * normalGMBolusVelocity(k, iEdge) * invAreaCell1
            end do
         end do
         ! Vertical velocity at bottom (maxLevelCell(iCell)+1) is zero, initialized above.
         do k=maxLevelCell(iCell),1,-1
            vertVelocityTop(k,iCell) = vertVelocityTop(k+1,iCell) - div_hu(k)
            vertTransportVelocityTop(k,iCell) = vertTransportVelocityTop(k+1,iCell) - div_huTransport(k)
            vertGMBolusVelocityTop(k,iCell) = vertGMBolusVelocityTop(k+1,iCell) - div_huGMBolus(k)
         end do         
      end do
      !$omp end do

      deallocate(div_hu,div_huTransport,div_huGMBolus)

      !$omp do schedule(runtime) private(eoe, i, k)
      do iEdge = 1, nEdges
         ! Compute v (tangential) velocities
         do i = 1, nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            do k = 1, maxLevelEdgeTop(iEdge) 
               tangentialVelocity(k,iEdge) = tangentialVelocity(k,iEdge) + weightsOnEdge(i,iEdge) * normalVelocity(k, eoe)
            end do
         end do
      end do
      !$omp end do

      !
      ! Compute kinetic energy
      !
      call mpas_pool_get_field(scratchPool, 'kineticEnergyVertex', kineticEnergyVertexField)
      call mpas_pool_get_field(scratchPool, 'kineticEnergyVertexOnCells', kineticEnergyVertexOnCellsField)
      call mpas_allocate_scratch_field(kineticEnergyVertexField, .true.)
      call mpas_allocate_scratch_field(kineticEnergyVertexOnCellsField, .true.)
      call mpas_threading_barrier()

      kineticEnergyVertex         => kineticEnergyVertexField % array
      kineticEnergyVertexOnCells  => kineticEnergyVertexOnCellsField % array

      !$omp workshare
      kineticEnergyVertex(:,:) = 0.0; 
      kineticEnergyVertexOnCells(:,:) = 0.0
      !$omp end workshare
      !$omp do schedule(runtime) private(i, iEdge, r_tmp, k)
      do iVertex = 1, nVertices*ke_vertex_flag
        do i = 1, vertexDegree
          iEdge = edgesOnVertex(i, iVertex)
          r_tmp = dcEdge(iEdge) * dvEdge(iEdge) * 0.25 / areaTriangle(iVertex)
          do k = 1, nVertLevels
            kineticEnergyVertex(k, iVertex) = kineticEnergyVertex(k, iVertex) + r_tmp * normalVelocity(k, iEdge)**2
          end do
        end do
      end do
      !$omp end do

      !$omp do schedule(runtime) private(invAreaCell1, i, j, iVertex, k)
      do iCell = 1, nCells*ke_vertex_flag
        invAreaCell1 = 1.0 / areaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          j = kiteIndexOnCell(i, iCell)
          iVertex = verticesOnCell(i, iCell)
          do k = 1, nVertLevels
            kineticEnergyVertexOnCells(k, iCell) = kineticEnergyVertexOnCells(k, iCell) + kiteAreasOnVertex(j, iVertex) * kineticEnergyVertex(k, iVertex) * invAreaCell1
          end do
        end do
      end do
      !$omp end do
      
      !
      ! Compute kinetic energy in each cell by blending kineticEnergyCell and kineticEnergyVertexOnCells
      !
      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells * ke_vertex_flag
         do k = 1, nVertLevels
            kineticEnergyCell(k,iCell) = 5.0 / 8.0 * kineticEnergyCell(k,iCell) + 3.0 / 8.0 * kineticEnergyVertexOnCells(k,iCell)
         end do
      end do
      !$omp end do

      call mpas_threading_barrier()
      call mpas_deallocate_scratch_field(kineticEnergyVertexField, .true.)
      call mpas_deallocate_scratch_field(kineticEnergyVertexOnCellsField, .true.)

      !
      ! Compute normalized relative and planetary vorticity
      !
      call mpas_pool_get_field(scratchPool, 'normalizedRelativeVorticityVertex', normalizedRelativeVorticityVertexField)
      call mpas_pool_get_field(scratchPool, 'normalizedPlanetaryVorticityVertex', normalizedPlanetaryVorticityVertexField)
      call mpas_allocate_scratch_field(normalizedRelativeVorticityVertexField, .true.)
      call mpas_allocate_scratch_field(normalizedPlanetaryVorticityVertexField, .true.)
      call mpas_threading_barrier()

      normalizedPlanetaryVorticityVertex  => normalizedPlanetaryVorticityVertexField % array
      normalizedRelativeVorticityVertex  => normalizedRelativeVorticityVertexField % array

      !$omp do schedule(runtime) private(invAreaTri1, k, layerThicknessVertex, i)
      do iVertex = 1, nVertices
         invAreaTri1 = 1.0 / areaTriangle(iVertex)
         do k = 1, maxLevelVertexBot(iVertex)
            layerThicknessVertex = 0.0
            do i = 1, vertexDegree
               layerThicknessVertex = layerThicknessVertex + layerThickness(k,cellsOnVertex(i,iVertex)) * kiteAreasOnVertex(i,iVertex)
            end do
            layerThicknessVertex = layerThicknessVertex * invAreaTri1

            normalizedRelativeVorticityVertex(k,iVertex) = relativeVorticity(k,iVertex) / layerThicknessVertex
            normalizedPlanetaryVorticityVertex(k,iVertex) = fVertex(iVertex) / layerThicknessVertex
         end do
      end do
      !$omp end do

      !$omp workshare
      normalizedRelativeVorticityEdge(:,:) = 0.0
      normalizedPlanetaryVorticityEdge(:,:) = 0.0
      !$omp end workshare
      !$omp do schedule(runtime) private(vertex1, vertex2, k)
      do iEdge = 1, nEdges
        vertex1 = verticesOnEdge(1, iEdge)
        vertex2 = verticesOnEdge(2, iEdge)
        do k = 1, maxLevelEdgeBot(iEdge)
          normalizedRelativeVorticityEdge(k, iEdge) = 0.5 * (normalizedRelativeVorticityVertex(k, vertex1) + normalizedRelativeVorticityVertex(k, vertex2))
          normalizedPlanetaryVorticityEdge(k, iEdge) = 0.5 * (normalizedPlanetaryVorticityVertex(k, vertex1) + normalizedPlanetaryVorticityVertex(k, vertex2))
        end do
      end do
      !$omp end do

      !$omp workshare
      normalizedRelativeVorticityCell(:,:) = 0.0
      !$omp end workshare
      !$omp do schedule(runtime) private(invAreaCell1, i, j, iVertex, k)
      do iCell = 1, nCells
        invAreaCell1 = 1.0 / areaCell(iCell)

        do i = 1, nEdgesOnCell(iCell)
          j = kiteIndexOnCell(i, iCell)
          iVertex = verticesOnCell(i, iCell)
          do k = 1, maxLevelCell(iCell)
            normalizedRelativeVorticityCell(k, iCell) = normalizedRelativeVorticityCell(k, iCell) &
              + kiteAreasOnVertex(j, iVertex) * normalizedRelativeVorticityVertex(k, iVertex) * invAreaCell1
          end do
        end do
      end do
      !$omp end do

      ! Diagnostics required for the Anticipated Potential Vorticity Method (apvm).
      if (config_apvm_scale_factor>1e-10) then

         call mpas_pool_get_field(scratchPool, 'vorticityGradientNormalComponent', vorticityGradientNormalComponentField)
         call mpas_pool_get_field(scratchPool, 'vorticityGradientTangentialComponent', vorticityGradientTangentialComponentField)
         call mpas_allocate_scratch_field(vorticityGradientNormalComponentField, .true.)
         call mpas_allocate_scratch_field(vorticityGradientTangentialComponentField, .true.)
         call mpas_threading_barrier()

         vorticityGradientNormalComponent => vorticityGradientNormalComponentField % array
         vorticityGradientTangentialComponent => vorticityGradientTangentialComponentField % array

         !$omp do schedule(runtime) private(cell1, cell2, vertex1, vertex2, invLength, k)
         do iEdge = 1,nEdges
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            vertex1 = verticesOnedge(1, iEdge)
            vertex2 = verticesOnedge(2, iEdge)

            invLength = 1.0 / dcEdge(iEdge)
            ! Compute gradient of PV in normal direction
            !   ( this computes the gradient for all edges bounding real cells )
            do k=1,maxLevelEdgeTop(iEdge)
               vorticityGradientNormalComponent(k,iEdge) = &
                  (normalizedRelativeVorticityCell(k,cell2) - normalizedRelativeVorticityCell(k,cell1)) * invLength
            enddo

            invLength = 1.0 / dvEdge(iEdge)
            ! Compute gradient of PV in the tangent direction
            !   ( this computes the gradient at all edges bounding real cells and distance-1 ghost cells )
            do k = 1,maxLevelEdgeBot(iEdge)
              vorticityGradientTangentialComponent(k,iEdge) = &
                 (normalizedRelativeVorticityVertex(k,vertex2) - normalizedRelativeVorticityVertex(k,vertex1)) * invLength
            enddo

         enddo
         !$omp end do

         !
         ! Modify PV edge with upstream bias.
         !
         !$omp do schedule(runtime) private(k)
         do iEdge = 1,nEdges
            do k = 1,maxLevelEdgeBot(iEdge)
              normalizedRelativeVorticityEdge(k,iEdge) = normalizedRelativeVorticityEdge(k,iEdge) &
                - config_apvm_scale_factor * dt * &
                    (  normalVelocity(k,iEdge)     * vorticityGradientNormalComponent(k,iEdge)      &
                     + tangentialVelocity(k,iEdge) * vorticityGradientTangentialComponent(k,iEdge) )
            enddo
         enddo
         !$omp end do

         call mpas_threading_barrier()
         call mpas_deallocate_scratch_field(vorticityGradientNormalComponentField, .true.)
         call mpas_deallocate_scratch_field(vorticityGradientTangentialComponentField, .true.)

      endif

      call mpas_threading_barrier()
      call mpas_deallocate_scratch_field(normalizedRelativeVorticityVertexField, .true.)
      call mpas_deallocate_scratch_field(normalizedPlanetaryVorticityVertexField, .true.)

      !
      ! equation of state
      !
      call mpas_timer_start("equation of state", .false., diagEOSTimer)
      call mpas_threading_barrier()

      ! compute in-place density
      if (config_pressure_gradient_type.eq.'Jacobian_from_TS') then
         ! only compute EOS derivatives if needed.
         call mpas_pool_get_array(diagnosticsPool, 'inSituThermalExpansionCoeff',inSituThermalExpansionCoeff)
         call mpas_pool_get_array(diagnosticsPool, 'inSituSalineContractionCoeff', inSituSalineContractionCoeff)
         call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, scratchPool, 0, 'relative', density, err, &
              inSituThermalExpansionCoeff, inSituSalineContractionCoeff, timeLevelIn=timeLevel)
      else
         call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, scratchPool, 0, 'relative', density, err, &
              timeLevelIn=timeLevel)
      endif
      call mpas_threading_barrier()

      ! compute potentialDensity, the density displaced adiabatically to the mid-depth of top layer.
      call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, scratchPool, 1, 'absolute', potentialDensity, err, timeLevelIn=timeLevel)

      ! compute displacedDensity, density displaced adiabatically to the mid-depth one layer deeper.  
      ! That is, layer k has been displaced to the depth of layer k+1.
      call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, scratchPool, 1, 'relative', displacedDensity, err, timeLevelIn=timeLevel)
      call mpas_threading_barrier()

      call mpas_timer_stop("equation of state", diagEOSTimer)
      call mpas_threading_barrier()

      !
      ! Pressure
      ! This section must be placed in the code after computing the density.
      !
      if (config_pressure_gradient_type.eq.'MontgomeryPotential') then

        ! use Montgomery Potential when layers are isopycnal.
        ! However, one may use 'pressure_and_zmid' when layers are isopycnal as well.
        ! Compute pressure at top of each layer, and then Montgomery Potential.

        allocate(pTop(nVertLevels))

        !$omp do schedule(runtime) private(k)
        do iCell = 1, nCells

           ! assume atmospheric pressure at the surface is zero for now.
           pTop(1) = 0.0
           ! At top layer it is g*SSH, where SSH may be off by a 
           ! constant (ie, bottomDepth can be relative to top or bottom)
           montgomeryPotential(1,iCell) = gravity &
              * (bottomDepth(iCell) + sum(layerThickness(1:nVertLevels,iCell)))

           do k = 2, nVertLevels
              pTop(k) = pTop(k-1) + density(k-1,iCell)*gravity* layerThickness(k-1,iCell)

              ! from delta M = p delta / density
              montgomeryPotential(k,iCell) = montgomeryPotential(k-1,iCell) &
                 + pTop(k)*(1.0/density(k,iCell) - 1.0/density(k-1,iCell)) 
           end do

        end do
        !$omp end do

        deallocate(pTop)

      else

        !$omp do schedule(runtime) private(k)
        do iCell = 1, nCells
           ! Pressure for generalized coordinates.
           ! Pressure at top surface may be due to atmospheric pressure
           ! or an ice-shelf depression. 
           pressure(1,iCell) = 0.0_RKIND
           if ( associated(frazilSurfacePressure) ) pressure(1,iCell) = pressure(1,iCell) + frazilSurfacePressure(iCell)
           pressure(1,iCell) = pressure(1,iCell) + seaSurfacePressure(iCell)
           pressure(1,iCell) = pressure(1,iCell) + density(1,iCell)*gravity*0.5*layerThickness(1,iCell)

           do k = 2, maxLevelCell(iCell)
              pressure(k,iCell) = pressure(k-1,iCell)  &
                + 0.5*gravity*(  density(k-1,iCell)*layerThickness(k-1,iCell) &
                               + density(k  ,iCell)*layerThickness(k  ,iCell))
           end do

           ! Compute zMid, the z-coordinate of the middle of the layer.
           ! Compute zTop, the z-coordinate of the top of the layer.
           ! Note the negative sign, since bottomDepth is positive
           ! and z-coordinates are negative below the surface.
           k = maxLevelCell(iCell)
           zMid(k:nVertLevels,iCell) = -bottomDepth(iCell) + 0.5*layerThickness(k,iCell)
           zTop(k:nVertLevels,iCell) = -bottomDepth(iCell) +     layerThickness(k,iCell)

           do k = maxLevelCell(iCell)-1, 1, -1
              zMid(k,iCell) = zMid(k+1,iCell)  &
                + 0.5*(  layerThickness(k+1,iCell) &
                       + layerThickness(k  ,iCell))
              zTop(k,iCell) = zTop(k+1,iCell)  &
                       + layerThickness(k  ,iCell)
           end do

           ! copy zTop(1,iCell) into sea-surface height array
           ssh(iCell) = zTop(1,iCell)

        end do
        !$omp end do

      endif

      !
      ! Brunt-Vaisala frequency (this has units of s^{-2})
      !
      coef = -gravity / rho_sw
      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells
         BruntVaisalaFreqTop(1,iCell) = 0.0
         do k = 2, maxLevelCell(iCell)
            BruntVaisalaFreqTop(k,iCell) = coef * (displacedDensity(k-1,iCell) - density(k,iCell)) & 
              / (zMid(k-1,iCell) - zMid(k,iCell))
          end do
      end do
      !$omp end do

      !
      ! Gradient Richardson number
      !
      RiTopOfCell = 100.0
      !$omp do schedule(runtime) private(invAreaCell1, k, shearSquared, i, iEdge, factor, delU2, shearMean)
      do iCell=1,nCells
         invAreaCell1 = 1.0 / areaCell(iCell)
         do k=2,maxLevelCell(iCell)
           shearSquared = 0.0
           do i = 1, nEdgesOnCell(iCell)
             iEdge = edgesOnCell(i, iCell)
             factor = 0.5 * dcEdge(iEdge) * dvEdge(iEdge) * invAreaCell1
             delU2 = (normalVelocity(k-1,iEdge) - normalVelocity(k,iEdge))**2
             shearSquared = shearSquared + factor * delU2
           enddo 
           shearMean = sqrt(shearSquared)
           shearMean = shearMean / (zMid(k-1,iCell) - zMid(k,iCell))
           RiTopOfCell(k,iCell) = BruntVaisalaFreqTop(k,iCell) / (shearMean**2 + 1.0e-10)
          end do
          RiTopOfCell(1,iCell) = RiTopOfCell(2,iCell)
      end do
      !$omp end do

      !
      ! extrapolate tracer values to ocean surface
      ! this eventually be a modelled process
      ! at present, just copy k=1 tracer values onto surface values
      ! field will be updated below is better approximations are available

      !$omp workshare
!TDR need to consider how to handel tracersSurfaceValues
      tracersSurfaceValue(:,:) = activeTracers(:,1,:)
      normalVelocitySurfaceLayer(:) = normalVelocity(1,:)
      !$omp end workshare

      !
      ! average tracer values over the ocean surface layer
      ! the ocean surface layer is generally assumed to be about 0.1 of the boundary layer depth
      if(config_use_cvmix_kpp) then

        !$omp workshare
        tracersSurfaceLayerValue(:,:) = 0.0
        indexSurfaceLayerDepth(:) = -9.e30
        !$omp end workshare
        !$omp do schedule(runtime) private(surfaceLayerDepth, sumSurfaceLayer, k, rSurfaceLayer)
        do iCell=1,nCells
          surfaceLayerDepth = config_cvmix_kpp_surface_layer_averaging
          sumSurfaceLayer=0.0
          do k=1,maxLevelCell(iCell)
           sumSurfaceLayer = sumSurfaceLayer + layerThickness(k,iCell)
           if(sumSurfaceLayer.gt.surfaceLayerDepth) then
             sumSurfaceLayer = sumSurfaceLayer - layerThickness(k,iCell)
             rSurfaceLayer = int(k-1) + (surfaceLayerDepth-sumSurfaceLayer)/layerThickness(k,iCell)
             indexSurfaceLayerDepth(iCell) = rSurfaceLayer
             exit
           endif
          end do
          do k=1,int(rSurfaceLayer)
            tracersSurfaceLayerValue(:,iCell) = tracersSurfaceLayerValue(:,iCell) + activeTracers(:,k,iCell)*layerThickness(k,iCell)
          enddo
          k=int(rSurfaceLayer)+1
          tracersSurfaceLayerValue(:,iCell) = tracersSurfaceLayerValue(:,iCell) + fraction(rSurfaceLayer)*activeTracers(:,k,iCell)*layerThickness(k,iCell)
          tracersSurfaceLayerValue(:,iCell) = tracersSurfaceLayerValue(:,iCell) / surfaceLayerDepth
        enddo
        !$omp end do

        !
        ! average normal velocity values over the ocean surface layer
        ! the ocean surface layer is generally assumed to be about 0.1 of the boundary layer depth
        !
        !$omp workshare
        normalVelocitySurfaceLayer(:) = 0.0_RKIND
        !$omp end workshare
        !$omp do schedule(runtime) private(cell1, cell2, surfaceLayerDepth, sumSurfaceLayer, k, rSurfaceLayer)
        do iEdge=1,nEdges
          cell1=cellsOnEdge(1,iEdge)
          cell2=cellsOnEdge(2,iEdge)
          surfaceLayerDepth = config_cvmix_kpp_surface_layer_averaging
          sumSurfaceLayer=0.0
          rSurfaceLayer = min(1, maxLevelEdgeTop(iEdge))
          do k=1,maxLevelEdgeTop(iEdge)
           rSurfaceLayer = k
           sumSurfaceLayer = sumSurfaceLayer + layerThicknessEdge(k,iEdge)
           if(sumSurfaceLayer.gt.surfaceLayerDepth) then
             sumSurfaceLayer = sumSurfaceLayer - layerThicknessEdge(k,iEdge)
             rSurfaceLayer = int(k-1) + (surfaceLayerDepth-sumSurfaceLayer)/layerThicknessEdge(k,iEdge)
             exit
           endif
          end do
          do k=1,int(rSurfaceLayer)
            normalVelocitySurfaceLayer(iEdge) = normalVelocitySurfaceLayer(iEdge) + normalVelocity(k,iEdge)*layerThicknessEdge(k,iEdge)
          enddo
          k=int(rSurfaceLayer)+1
          if(k.le.maxLevelEdgeTop(iEdge)) then
            normalVelocitySurfaceLayer(iEdge) = normalVelocitySurfaceLayer(iEdge) + fraction(rSurfaceLayer)*normalVelocity(k,iEdge)*layerThicknessEdge(k,iEdge)
            normalVelocitySurfaceLayer(iEdge) = normalVelocitySurfaceLayer(iEdge) / surfaceLayerDepth
          end if
        enddo
        !$omp end do

        !
        !  compute fields used as intent(in) to CVMix/KPP
        call ocn_compute_KPP_input_fields(statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, timeLevel)
      endif

      ! compute the attenuation coefficient for surface fluxes
      !$omp workshare
      surfaceFluxAttenuationCoefficient(:) = config_flux_attenuation_coefficient
      !$omp end workshare

      !
      !  compute fields needed to compute land-ice fluxes, either in the ocean model or in the coupler
      call mpas_timer_start("land_ice_diagnostic_fields", .false.)
      call ocn_compute_land_ice_flux_input_fields(meshPool, statePool, forcingPool, scratchPool, &
                                         diagnosticsPool, timeLevel)
      call mpas_timer_stop("land_ice_diagnostic_fields")

      !$omp do schedule(runtime) private(cell1, cell2)
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)

         gradSSH(1, iEdge) = (ssh(cell2) - ssh(cell1)) / dcEdge(iEdge)
      end do
      !$omp end do

     call mpas_threading_barrier()



   end subroutine ocn_diagnostic_solve!}}}

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
   subroutine ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, oldLayerThickness, layerThicknessEdge, &
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

      type (mpas_pool_type), intent(in) :: scratchPool !< Input: scratch variables

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         oldLayerThickness    !< Input: layer thickness at old time

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdge     !< Input: layerThickness interpolated to an edge

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

      integer :: iEdge, iCell, k, i
      integer, pointer :: nCells, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, &
        maxLevelCell, maxLevelEdgeBot
      integer, dimension(:,:), pointer :: edgesOnCell, edgeSignOnCell

      real (kind=RKIND) :: flux, invAreaCell, div_hu_btr
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real (kind=RKIND), dimension(:), pointer :: &
         projectedSSH       !> projected SSH at new time
      type (field1DReal), pointer :: projectedSSHField
      real (kind=RKIND), dimension(:,:), pointer :: &
         ALE_Thickness, & !> ALE thickness at new time
         div_hu           !> divergence of (thickness*velocity)
      type (field2DReal), pointer :: ALE_ThicknessField, div_huField

      character (len=StrKIND), pointer :: config_vert_coord_movement

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_vert_coord_movement', config_vert_coord_movement)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', maxLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      if (config_vert_coord_movement.eq.'impermeable_interfaces') then
        vertAleTransportTop=0.0
        return
      end if


      call mpas_pool_get_field(scratchPool, 'div_hu', div_huField)
      call mpas_pool_get_field(scratchPool, 'projectedSSH', projectedSSHField)
      call mpas_pool_get_field(scratchPool, 'ALE_Thickness', ALE_ThicknessField)
      call mpas_allocate_scratch_field(div_huField, .true.)
      call mpas_allocate_scratch_field(projectedSSHField, .true.)
      call mpas_allocate_scratch_field(ALE_ThicknessField, .true.)

      call mpas_threading_barrier()

      div_hu => div_huField % array
      projectedSSH => projectedSSHField % array
      ALE_Thickness => ALE_ThicknessField % array

      !
      ! thickness-weighted divergence and barotropic divergence
      !
      ! See Ringler et al. (2010) jcp paper, eqn 19, 21, and fig. 3.
      !$omp do schedule(runtime) private(invAreaCell, i, iEdge, k, flux, div_hu_btr)
      do iCell = 1, nCells
         div_hu(:,iCell) = 0.0_RKIND
         div_hu_btr = 0.0_RKIND
         invAreaCell = 1.0_RKIND / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)

            do k = 1, maxLevelEdgeBot(iEdge)
               flux = layerThicknessEdge(k, iEdge) * normalVelocity(k, iEdge) * dvEdge(iEdge) * edgeSignOnCell(i, iCell) * invAreaCell
               div_hu(k,iCell) = div_hu(k,iCell) - flux
               div_hu_btr = div_hu_btr - flux
            end do
         end do
         projectedSSH(iCell) = oldSSH(iCell) - dt*div_hu_btr
      end do
      !$omp end do

      !
      ! Compute desired thickness at new time
      !
      if (present(newHighFreqThickness)) then
        call ocn_ALE_thickness(meshPool, verticalMeshPool, projectedSSH, ALE_thickness, err, newHighFreqThickness)
      else
        call ocn_ALE_thickness(meshPool, verticalMeshPool, projectedSSH, ALE_thickness, err)
      endif

      call mpas_threading_barrier()

      !
      ! Vertical transport through layer interfaces
      !
      ! Vertical transport through layer interface at top and bottom is zero.
      ! Here we are using solving the continuity equation for vertAleTransportTop ($w^t$),
      ! and using ALE_Thickness for thickness at the new time.

      !$omp do schedule(runtime) private(k)
      do iCell = 1,nCells
         vertAleTransportTop(1,iCell) = 0.0_RKIND
         vertAleTransportTop(maxLevelCell(iCell)+1,iCell) = 0.0_RKIND
         do k = maxLevelCell(iCell),2,-1
            vertAleTransportTop(k,iCell) = vertAleTransportTop(k+1,iCell) - div_hu(k,iCell) &
              - (ALE_Thickness(k,iCell) - oldLayerThickness(k,iCell))/dt
         end do
      end do
      !$omp end do

      call mpas_threading_barrier()
      call mpas_deallocate_scratch_field(div_huField, .true.)
      call mpas_deallocate_scratch_field(projectedSSHField, .true.)
      call mpas_deallocate_scratch_field(ALE_ThicknessField, .true.)

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

      integer :: iEdge, cell1, cell2, eoe, i, j, k
      integer, pointer :: nEdgesSolve
      real (kind=RKIND), dimension(:), pointer :: fEdge
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge, normalVelocity, normalBaroclinicVelocity
      type (dm_info) :: dminfo

      integer, dimension(:), pointer :: maxLevelEdgeTop, nEdgesOnEdge
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnEdge

      integer :: timeLevel

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_timer_start("ocn_fuperp")
      call mpas_threading_barrier()

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocity, timeLevel)

      call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
      call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)

      call mpas_pool_get_array(meshPool, 'fEdge', fEdge)

      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

      !DWJ: ADD OMP (Only needed for split explicit)

      !
      ! Put f*normalBaroclinicVelocity^{perp} in u as a work variable
      !
      !$omp do schedule(runtime) private(cell1, cell2, k, eoe)
      do iEdge = 1, nEdgesSolve
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k = 1, maxLevelEdgeTop(iEdge)

            normalVelocity(k,iEdge) = 0.0
            do j = 1,nEdgesOnEdge(iEdge)
               eoe = edgesOnEdge(j,iEdge)
               normalVelocity(k,iEdge) = normalVelocity(k,iEdge) + weightsOnEdge(j,iEdge) * normalBaroclinicVelocity(k,eoe) * fEdge(eoe) 
            end do
         end do
      end do
      !$omp end do

      call mpas_threading_barrier()
      call mpas_timer_stop("ocn_fuperp")

   end subroutine ocn_fuperp!}}}

!***********************************************************************
!
!  routine ocn_filter_btr_mode_vel
!
!> \brief   filters barotropic mode out of the velocity variable.
!> \author  Mark Petersen
!> \date    23 September 2011
!> \details 
!>  This routine filters barotropic mode out of the velocity variable.
!
!-----------------------------------------------------------------------
   subroutine ocn_filter_btr_mode_vel(statePool, diagnosticsPool, meshPool, timeLevelIn)!{{{

      type (mpas_pool_type), intent(inout) :: statePool !< Input/Output: State information
      type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: Diagnostics information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      integer, intent(in), optional :: timeLevelIn !< Input: Time level for state pool

      integer :: iEdge, k
      integer, pointer :: nEdges
      real (kind=RKIND) :: vertSum, normalThicknessFluxSum, thicknessSum
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge, normalVelocity
      integer, dimension(:), pointer :: maxLevelEdgeTop

      integer :: timeLevel

      call mpas_timer_start("ocn_filter_btr_mode_vel")

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)

      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      !$omp do schedule(runtime) private(normalThicknessFluxSum, thicknessSum, k, vertSum)
      do iEdge = 1, nEdges

        ! thicknessSum is initialized outside the loop because on land boundaries 
        ! maxLevelEdgeTop=0, but I want to initialize thicknessSum with a 
        ! nonzero value to avoid a NaN.
        normalThicknessFluxSum = layerThicknessEdge(1,iEdge) * normalVelocity(1,iEdge)
        thicknessSum  = layerThicknessEdge(1,iEdge)

        do k = 2, maxLevelEdgeTop(iEdge)
          normalThicknessFluxSum = normalThicknessFluxSum + layerThicknessEdge(k,iEdge) * normalVelocity(k,iEdge)
          thicknessSum  =  thicknessSum + layerThicknessEdge(k,iEdge)
        enddo

        vertSum = normalThicknessFluxSum/thicknessSum
        do k = 1, maxLevelEdgeTop(iEdge)
          normalVelocity(k,iEdge) = normalVelocity(k,iEdge) - vertSum
        enddo
      enddo ! iEdge
      !$omp end do

      call mpas_threading_barrier()

      call mpas_timer_stop("ocn_filter_btr_mode_vel")

   end subroutine ocn_filter_btr_mode_vel!}}}

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
   subroutine ocn_filter_btr_mode_tend_vel(tendPool, statePool, diagnosticsPool, meshPool, timeLevelIn)!{{{

      type (mpas_pool_type), intent(inout) :: tendPool !< Input/Output: Tendency information
      type (mpas_pool_type), intent(in) :: statePool !< Input: State information
      type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: Diagnostics information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      integer, intent(in), optional :: timeLevelIn !< Input: Time level for state pool

      integer :: iEdge, k
      integer, pointer :: nEdges
      real (kind=RKIND) :: vertSum, normalThicknessFluxSum, thicknessSum
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge, tend_normalVelocity

      integer, dimension(:), pointer :: maxLevelEdgeTop

      integer :: timeLevel

      call mpas_timer_start("ocn_filter_btr_mode_tend_vel")

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_array(tendPool, 'normalVelocity', tend_normalVelocity)

      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      !$omp do schedule(runtime) private(normalThicknessFluxSum, thicknessSum, vertSum, k)
      do iEdge = 1, nEdges

        ! thicknessSum is initialized outside the loop because on land boundaries 
        ! maxLevelEdgeTop=0, but I want to initialize thicknessSum with a 
        ! nonzero value to avoid a NaN.
        normalThicknessFluxSum = layerThicknessEdge(1,iEdge) * tend_normalVelocity(1,iEdge)
        thicknessSum  = layerThicknessEdge(1,iEdge)

        do k = 2, maxLevelEdgeTop(iEdge)
          normalThicknessFluxSum = normalThicknessFluxSum + layerThicknessEdge(k,iEdge) * tend_normalVelocity(k,iEdge)
          thicknessSum  =  thicknessSum + layerThicknessEdge(k,iEdge)
        enddo

        vertSum = normalThicknessFluxSum / thicknessSum
        do k = 1, maxLevelEdgeTop(iEdge)
          tend_normalVelocity(k,iEdge) = tend_normalVelocity(k,iEdge) - vertSum
        enddo
      enddo ! iEdge
      !$omp end do

      call mpas_timer_stop("ocn_filter_btr_mode_tend_vel")

   end subroutine ocn_filter_btr_mode_tend_vel!}}}

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
   subroutine ocn_diagnostics_init(err)!{{{
      integer, intent(out) :: err !< Output: Error flag

      logical, pointer :: config_include_KE_vertex
      character (len=StrKIND), pointer :: config_time_integrator

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_include_KE_vertex', config_include_KE_vertex)
      call mpas_pool_get_config(ocnConfigs, 'config_time_integrator', config_time_integrator)

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
        .or.trim(config_time_integrator) == 'unsplit_explicit') then
          ! For split explicit, PV is eta/h because the Coriolis term 
          ! is added separately to the momentum tendencies.
          fCoef = 0
      end if

    end subroutine ocn_diagnostics_init!}}}

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
!>       bulkRichardsonNumberBuoy
!>       bulkRichardsonNumberShear
!>
!
!-----------------------------------------------------------------------

    subroutine ocn_compute_KPP_input_fields(statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, timeLevelIn)!{{{

      type (mpas_pool_type), intent(in) :: statePool !< Input/Output: State information
      type (mpas_pool_type), intent(in) :: forcingPool !< Input: Forcing information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      type (mpas_pool_type), intent(inout) :: diagnosticsPool !< Diagnostics information derived from State
      type (mpas_pool_type), intent(in) :: scratchPool !< Input: scratch variables
      integer, intent(in), optional :: timeLevelIn

      ! pool pointers
      type (mpas_pool_type), pointer :: tracersSurfaceFluxPool

      ! scalars
      integer, pointer :: nCells, nVertLevels

      ! integer pointers
      integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell

      ! real pointers
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaCell
      real (kind=RKIND), dimension(:), pointer :: penetrativeTemperatureFlux, surfaceThicknessFlux, &
           surfaceBuoyancyForcing, surfaceFrictionVelocity, penetrativeTemperatureFluxOBL, &
           normalVelocitySurfaceLayer
      real (kind=RKIND), dimension(:), pointer :: surfaceStress, surfaceStressMagnitude
      real (kind=RKIND), dimension(:,:), pointer ::  &
           layerThickness, zMid, zTop, densitySurfaceDisplaced, density, &
           normalVelocity, activeTracersSurfaceFlux, thermalExpansionCoeff, salineContractionCoeff

      real (kind=RKIND), dimension(:), pointer :: &
           indexSurfaceLayerDepth

      real (kind=RKIND), dimension(:,:), pointer ::  & 
           bulkRichardsonNumberBuoy, bulkRichardsonNumberShear

      ! local
      integer :: iCell, iEdge, i, k, err, timeLevel
      integer, pointer :: indexTempFlux, indexSaltFlux
      real (kind=RKIND) :: numerator, denominator, turbulentVelocitySquared
      real (kind=RKIND) :: buoyContribution, shearContribution, factor, deltaVelocitySquared, delU2, invAreaCell

      type (field2DReal), pointer :: densitySurfaceDisplacedField, thermalExpansionCoeffField, salineContractionCoeffField

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux', tracersSurfaceFluxPool)

      ! set the parameter turbulentVelocitySquared
      turbulentVelocitySquared = 0.001_RKIND

      ! set scalar values
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(tracersSurfaceFluxPool, 'index_temperatureSurfaceFlux', indexTempFlux)
      call mpas_pool_get_dimension(tracersSurfaceFluxPool, 'index_salinitySurfaceFlux', indexSaltFlux)

      ! set pointers into state, mesh, diagnostics and scratch
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)

      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
      call mpas_pool_get_array(diagnosticsPool, 'zTop', zTop)
      call mpas_pool_get_array(diagnosticsPool, 'density', density)
      call mpas_pool_get_array(diagnosticsPool, 'surfaceFrictionVelocity', surfaceFrictionVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'penetrativeTemperatureFluxOBL', penetrativeTemperatureFluxOBL)
      call mpas_pool_get_array(diagnosticsPool, 'bulkRichardsonNumberBuoy', bulkRichardsonNumberBuoy)
      call mpas_pool_get_array(diagnosticsPool, 'bulkRichardsonNumberShear', bulkRichardsonNumberShear)
      call mpas_pool_get_array(diagnosticsPool, 'indexSurfaceLayerDepth', indexSurfaceLayerDepth)
      call mpas_pool_get_array(diagnosticsPool, 'surfaceBuoyancyForcing', surfaceBuoyancyForcing)
      call mpas_pool_get_array(diagnosticsPool, 'normalVelocitySurfaceLayer', normalVelocitySurfaceLayer)

      call mpas_pool_get_array(forcingPool, 'surfaceThicknessFlux', surfaceThicknessFlux)
      call mpas_pool_get_array(forcingPool, 'penetrativeTemperatureFlux', penetrativeTemperatureFlux)
      call mpas_pool_get_array(forcingPool, 'surfaceStress', surfaceStress)
      call mpas_pool_get_array(forcingPool, 'surfaceStressMagnitude', surfaceStressMagnitude)

      call mpas_pool_get_array(tracersSurfaceFluxPool, 'activeTracersSurfaceFlux', activeTracersSurfaceFlux)

      ! allocate scratch space displaced density computation
      call mpas_pool_get_field(scratchPool, 'densitySurfaceDisplaced', densitySurfaceDisplacedField)
      call mpas_pool_get_field(scratchPool, 'thermalExpansionCoeff', thermalExpansionCoeffField)
      call mpas_pool_get_field(scratchPool, 'salineContractionCoeff', salineContractionCoeffField)
      call mpas_allocate_scratch_field(densitySurfaceDisplacedField, .true.)
      call mpas_allocate_scratch_field(thermalExpansionCoeffField, .true.)
      call mpas_allocate_scratch_field(salineContractionCoeffField, .true.)
      call mpas_threading_barrier()

      densitySurfaceDisplaced => densitySurfaceDisplacedField % array
      thermalExpansionCoeff => thermalExpansionCoeffField % array
      salineContractionCoeff => salineContractionCoeffField % array

      ! compute EOS by displacing SST/SSS to every vertical layer in column
      call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, scratchPool, 0, 'surfaceDisplaced', densitySurfaceDisplaced, err, &
              thermalExpansionCoeff, salineContractionCoeff, timeLevel)

      !$omp do schedule(runtime) private(invAreaCell, deltaVelocitySquared, i, iEdge, factor, delU2, k, buoyContribution, shearContribution)
      do iCell = 1, nCells
       invAreaCell = 1.0 / areaCell(iCell)

       ! compute surface buoyancy forcing based on surface fluxes of mass, temperature, salinity and frazil (frazil to be added later)
       ! since this computation is confusing, variables, units and sign convention is repeated here
       ! everything below should be consistent with that specified in Registry
       ! everything below should be consistent with the CVMix/KPP documentation: https://www.dropbox.com/s/6hqgc0rsoa828nf/cvmix_20aug2013.pdf
       !
       !    surfaceThicknessFlux: surface mass flux, m/s, positive into ocean
       !    activeTracersSurfaceFlux(indexTempFlux): non-penetrative temperature flux, C m/s, positive into ocean
       !    penetrativeTemperatureFlux: penetrative surface temperature flux at ocean surface, positive into ocean
       !    activeTracersSurfaceFlux(indexSaltFlux): salinity flux, PSU m/s, positive into ocean
       !    penetrativeTemperatureFluxOBL: penetrative temperature flux computed at z=OBL, positive down
       !
       ! note: the following fields used the CVMix/KPP computation of buoyancy forcing are not included here
       !    1. Tm: temperature associated with surfaceThicknessFlux, C  (here we assume Tm == temperatureSurfaceValue)
       !    2. Sm: salinity associated with surfaceThicknessFlux, PSU (here we assume Sm == salinitySurfaceValue and account for salinity flux in activeTracersSurfaceFlux array)
       !
         surfaceBuoyancyForcing(iCell) =  thermalExpansionCoeff (1,iCell) *  &
               (activeTracersSurfaceFlux(indexTempFlux,iCell) + penetrativeTemperatureFlux(iCell) - penetrativeTemperatureFluxOBL(iCell)) &
              - salineContractionCoeff(1,iCell) *  activeTracersSurfaceFlux(indexSaltFlux,iCell)
        
       ! at this point, surfaceBuoyancyForcing has units of m/s 
       ! change into units of m^2/s^3 (which can be thought of as the flux of buoyancy, units of buoyancy * velocity )
         surfaceBuoyancyForcing(iCell) = surfaceBuoyancyForcing(iCell) * gravity

       ! compute magnitude of surface stress
        deltaVelocitySquared = 0.0
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          factor = 0.5 * dcEdge(iEdge) * dvEdge(iEdge) * invAreaCell
          delU2 =  (surfaceStress(iEdge))**2
          deltaVelocitySquared = deltaVelocitySquared + factor * delU2
        enddo
        surfaceStressMagnitude(iCell) = sqrt(deltaVelocitySquared)

       ! compute surface friction velocity
         surfaceFrictionVelocity(iCell) = sqrt(surfaceStressMagnitude(iCell) / rho_sw)


      enddo
      !$omp end do

      call mpas_threading_barrier()

      ! deallocate scratch space
      call mpas_deallocate_scratch_field(densitySurfaceDisplacedField, .true.)
      call mpas_deallocate_scratch_field(thermalExpansionCoeffField, .true.)
      call mpas_deallocate_scratch_field(salineContractionCoeffField, .true.)

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

   subroutine ocn_compute_land_ice_flux_input_fields(meshPool, statePool, &
      forcingPool, scratchPool, diagnosticsPool, timeLevel)!{{{

      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      type (mpas_pool_type), intent(in) :: statePool !< Input: State information
      type (mpas_pool_type), intent(in) :: forcingPool !< Input: Forcing information
      type (mpas_pool_type), intent(in) :: scratchPool !< Input/Output: scratch variables
      type (mpas_pool_type), intent(inout) :: diagnosticsPool !< Input/Output: Diagnostics information

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: tracersPool

      integer :: iCell, iEdge, cell1, cell2, iLevel, i
      integer, pointer :: nCells, nEdges

      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, cellMask

      integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell

      integer, pointer :: indexT, indexS, indexBLT, indexBLS, indexHeatTrans, indexSaltTrans

      character (len=StrKIND), pointer :: config_land_ice_flux_formulation, config_land_ice_flux_mode

      real (kind=RKIND), pointer :: config_land_ice_flux_boundaryLayerThickness, &
                                    config_land_ice_flux_boundaryLayerNeighborWeight, &
                                    config_land_ice_flux_topDragCoeff, &
                                    config_land_ice_flux_rms_tidal_velocity, &
                                    config_land_ice_flux_jenkins_heat_transfer_coefficient, &
                                    config_land_ice_flux_jenkins_salt_transfer_coefficient, &
                                    config_land_ice_flux_attenuation_coefficient

      real (kind=RKIND) :: blThickness, dz, weightSum, h_nu, Gamma_turb, landIceEdgeFraction, velocityMagnitude

      real (kind=RKIND), dimension(:), pointer :: landIceFraction, &
                                                  landIceFrictionVelocity, &
                                                  topDrag, &
                                                  topDragMagnitude, &
                                                  fCell, &
                                                  blTempScratch, blSaltScratch, &
                                                  surfaceFluxAttenuationCoefficient

      real (kind=RKIND), dimension(:,:), pointer :: kineticEnergyCell, layerThickness, normalVelocity, &
                                                    landIceBoundaryLayerTracers, landIceTracerTransferVelocities
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
      type (field1DReal), pointer :: boundaryLayerTemperatureField, boundaryLayerSalinityField

      logical :: jenkinsOn, hollandJenkinsOn

      ! constants for Holland and Jenkins 1999 parameterization of the boundary layer
      real (kind=RKIND), parameter :: &
         Pr = 13.8_RKIND, &             ! the Prandtl number
         Sc = 2432.0_RKIND, &           ! the Schmidt number
         nuSaltWater = 1.95e-6_RKIND, & ! molecular viscosity of sea water (m^2/s)
         kVonKarman = 0.4_RKIND, &      ! the von Karman constant
         xiN = 0.052_RKIND              ! dimensionless planetary boundary layer constant


      call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)
      if ( trim(config_land_ice_flux_mode) == 'off') then
         return
      end if


      jenkinsOn = .false.
      hollandJenkinsOn = .false.
      call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_formulation', config_land_ice_flux_formulation)
      if ( trim(config_land_ice_flux_formulation) == 'Jenkins' ) then
         jenkinsOn = .true.
      else if ( trim(config_land_ice_flux_formulation) == 'HollandJenkins' ) then
         hollandJenkinsOn = .true.
      end if

      call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_topDragCoeff', config_land_ice_flux_topDragCoeff)
      call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_boundaryLayerThickness', config_land_ice_flux_boundaryLayerThickness)
      call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_boundaryLayerNeighborWeight', config_land_ice_flux_boundaryLayerNeighborWeight)
      call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_rms_tidal_velocity', config_land_ice_flux_rms_tidal_velocity)
      call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_attenuation_coefficient', config_land_ice_flux_attenuation_coefficient)

      if(jenkinsOn) then
         call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_jenkins_heat_transfer_coefficient', config_land_ice_flux_jenkins_heat_transfer_coefficient)
         call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_jenkins_salt_transfer_coefficient', config_land_ice_flux_jenkins_salt_transfer_coefficient)
      end if

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellMask', cellMask)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexT)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexS)

      call mpas_pool_get_array(forcingPool, 'landIceFraction', landIceFraction)

      call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCell)

      call mpas_pool_get_array(diagnosticsPool, 'landIceFrictionVelocity', landIceFrictionVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'topDrag', topDrag)
      call mpas_pool_get_array(diagnosticsPool, 'topDragMagnitude', topDragMagnitude)

      call mpas_pool_get_array(diagnosticsPool, 'landIceBoundaryLayerTracers', landIceBoundaryLayerTracers)
      call mpas_pool_get_dimension(diagnosticsPool, 'index_landIceBoundaryLayerTemperature', indexBLT)
      call mpas_pool_get_dimension(diagnosticsPool, 'index_landIceBoundaryLayerSalinity', indexBLS)

      if(jenkinsOn .or. hollandJenkinsOn) then
        call mpas_pool_get_array(diagnosticsPool, 'landIceTracerTransferVelocities', landIceTracerTransferVelocities)
        call mpas_pool_get_dimension(diagnosticsPool, 'index_landIceHeatTransferVelocity', indexHeatTrans)
        call mpas_pool_get_dimension(diagnosticsPool, 'index_landIceSaltTransferVelocity', indexSaltTrans)
      end if
      call mpas_pool_get_array(diagnosticsPool, 'surfaceFluxAttenuationCoefficient', surfaceFluxAttenuationCoefficient)

      call mpas_pool_get_field(scratchPool, 'boundaryLayerTemperatureScratch', boundaryLayerTemperatureField)
      call mpas_pool_get_field(scratchPool, 'boundaryLayerSalinityScratch', boundaryLayerSalinityField)
      call mpas_allocate_scratch_field(boundaryLayerTemperatureField, .true.)
      call mpas_allocate_scratch_field(boundaryLayerSalinityField, .true.)
      call mpas_threading_barrier()
      blTempScratch => boundaryLayerTemperatureField % array
      blSaltScratch => boundaryLayerSalinityField % array

      if(hollandJenkinsOn) then
         call mpas_pool_get_array(meshPool, 'fCell', fCell)
      end if

      ! Compute top drag
      !$omp do schedule(runtime) private(cell1, cell2, velocityMagnitude, landIceEdgeFraction)
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)

         ! top drag tau = - CD*|u|*u, where |u| = sqrt(2*KE) = sqrt(KE1 + KE2) from the neighboring cells
         velocityMagnitude = sqrt(kineticEnergyCell(1,cell1) + kineticEnergyCell(1,cell2))
         landIceEdgeFraction = 0.5_RKIND*(landIceFraction(cell1)+landIceFraction(cell2))

         topDrag(iEdge) = - rho_sw * landIceEdgeFraction * config_land_ice_flux_topDragCoeff &
                          * velocityMagnitude * normalVelocity(1,iEdge)

      end do
      !$omp end do

      ! compute top drag magnitude and friction velocity at cell centers
      !$omp do schedule(runtime)
      do iCell = 1, nCells
         ! the magnitude of the top drag is CD*u**2 = CD*(2*KE)
         topDragMagnitude(iCell) = rho_sw * landIceFraction(iCell) &
                                   * 2.0_RKIND * config_land_ice_flux_topDragCoeff *  kineticEnergyCell(1,iCell)

         ! the friction velocity is the square root of the top drag + variance of tidal velocity (computed regardless of land-ice coverage)
         landIceFrictionVelocity(iCell) = sqrt(config_land_ice_flux_topDragCoeff* (2.0_RKIND * kineticEnergyCell(1,iCell) &
                                   + config_land_ice_flux_rms_tidal_velocity))
      end do
      !$omp end do



      ! average temperature and salinity over horizontal neighbors and the sub-ice-shelf boundary layer
      !$omp do schedule(runtime) private(blThickness, iLevel, dz)
      do iCell = 1, nCells
         blThickness = 0.0_RKIND
         blTempScratch(iCell) = 0.0_RKIND
         blSaltScratch(iCell) = 0.0_RKIND
         do iLevel = 1, maxLevelCell(iCell)
            dz = min(layerThickness(iLevel,iCell),config_land_ice_flux_boundaryLayerThickness-blThickness)
            if(dz <= 0.0_RKIND) exit
            blTempScratch(iCell) = blTempScratch(iCell) + activeTracers(indexT, iLevel, iCell)*dz
            blSaltScratch(iCell) = blSaltScratch(iCell) + activeTracers(indexS, iLevel, iCell)*dz
            blThickness = blThickness + dz
         end do
         if(blThickness > 0.0_RKIND) then
           blTempScratch(iCell) = blTempScratch(iCell)/blThickness
           blSaltScratch(iCell) = blSaltScratch(iCell)/blThickness
         end if
      end do
      !$omp end do

      !$omp do schedule(runtime) private(weightSum, i, cell2)
      do iCell = 1, nCells
         landIceBoundaryLayerTracers(indexBLT, iCell) = blTempScratch(iCell)
         landIceBoundaryLayerTracers(indexBLS, iCell) = blSaltScratch(iCell)
         if(config_land_ice_flux_boundaryLayerNeighborWeight > 0.0_RKIND) then
            weightSum = 1.0_RKIND
            do i = 1, nEdgesOnCell(iCell)
               cell2 = cellsOnCell(i,iCell)

               landIceBoundaryLayerTracers(indexBLT, iCell) = landIceBoundaryLayerTracers(indexBLT, iCell) &
                 + cellMask(1,cell2)*config_land_ice_flux_boundaryLayerNeighborWeight*blTempScratch(cell2)
               landIceBoundaryLayerTracers(indexBLS, iCell) = landIceBoundaryLayerTracers(indexBLS, iCell) &
                 + cellMask(1,cell2)*config_land_ice_flux_boundaryLayerNeighborWeight*blSaltScratch(cell2)
               weightSum = weightSum + cellMask(1,cell2)*config_land_ice_flux_boundaryLayerNeighborWeight
            end do
            landIceBoundaryLayerTracers(:, iCell) = landIceBoundaryLayerTracers(:, iCell)/weightSum
         end if
      end do
      !$omp end do

      if(jenkinsOn) then
         !$omp do schedule(runtime)
         do iCell = 1, nCells
            ! transfer coefficients from namelist
            landIceTracerTransferVelocities(indexHeatTrans, iCell) = landIceFrictionVelocity(iCell)*config_land_ice_flux_jenkins_heat_transfer_coefficient
            landIceTracerTransferVelocities(indexSaltTrans, iCell) = landIceFrictionVelocity(iCell)*config_land_ice_flux_jenkins_salt_transfer_coefficient
         end do
         !$omp end do
      else if(hollandJenkinsOn) then
         !$omp do schedule(runtime) private(h_nu, Gamma_turb)
         do iCell = 1, nCells
            ! friction-velocity dependent non-dimensional transfer coefficients from
            ! Holland and Jenkins 1999, (14)-(16) with eta_* = 1
            h_nu = 5.0_RKIND*nuSaltWater/landIceFrictionVelocity(iCell) ! uStar should never be zero because of tidal term
      
            Gamma_turb = 1.0_RKIND/(2.0_RKIND*xiN) - 1.0_RKIND/kVonKarman
            if(abs(fCell(iCell)) > 0.0_RKIND) then
              Gamma_turb = Gamma_turb + 1.0_RKIND/kVonKarman*log(landIceFrictionVelocity(iCell) &
                *xiN/(abs(fCell(iCell))*h_nu))
            end if
     
            landIceTracerTransferVelocities(indexHeatTrans, iCell) = 1.0_RKIND/(Gamma_turb + 12.5_RKIND*Pr**(2.0_RKIND/3.0_RKIND) - 6.0_RKIND)
            landIceTracerTransferVelocities(indexSaltTrans, iCell) = 1.0_RKIND/(Gamma_turb + 12.5_RKIND*Sc**(2.0_RKIND/3.0_RKIND) - 6.0_RKIND)
         end do
         !$omp end do
      end if

      call mpas_threading_barrier()
      call mpas_deallocate_scratch_field(boundaryLayerTemperatureField, .true.)
      call mpas_deallocate_scratch_field(boundaryLayerSalinityField, .true.)

      ! recompute the spatially-varying attenuation coefficient based on landIceFraction
      !$omp do schedule(runtime)
      do iCell = 1, nCells
         surfaceFluxAttenuationCoefficient(iCell) = landIceFraction(iCell)*config_land_ice_flux_attenuation_coefficient &
           + (1.0_RKIND - landIceFraction(iCell))*surfaceFluxAttenuationCoefficient(iCell)
      end do
      !$omp end do

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

   subroutine ocn_reconstruct_gm_vectors(diagnosticsPool, meshPool) !{{{

      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: Diagnostic information

      real (kind=RKIND), dimension(:,:), pointer :: &
         normalTransportVelocity, transportVelocityX, transportVelocityY, transportVelocityZ, transportVelocityZonal, transportVelocityMeridional, &
         normalGMBolusVelocity, GMBolusVelocityX, GMBolusVelocityY, GMBolusVelocityZ, GMBolusVelocityZonal, GMBolusVelocityMeridional, &
         relativeSlopeTopOfEdge, relativeSlopeTopOfCellX, relativeSlopeTopOfCellY, relativeSlopeTopOfCellZ, relativeSlopeTopOfCellZonal, relativeSlopeTopOfCellMeridional, &
         gmStreamFuncTopOfEdge, GMStreamFuncX, GMStreamFuncY, GMStreamFuncZ, GMStreamFuncZonal, GMStreamFuncMeridional

         call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'transportVelocityX', transportVelocityX)
         call mpas_pool_get_array(diagnosticsPool, 'transportVelocityY', transportVelocityY)
         call mpas_pool_get_array(diagnosticsPool, 'transportVelocityZ', transportVelocityZ)
         call mpas_pool_get_array(diagnosticsPool, 'transportVelocityZonal', transportVelocityZonal)
         call mpas_pool_get_array(diagnosticsPool, 'transportVelocityMeridional', transportVelocityMeridional)

         call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'GMBolusVelocityX', GMBolusVelocityX)
         call mpas_pool_get_array(diagnosticsPool, 'GMBolusVelocityY', GMBolusVelocityY)
         call mpas_pool_get_array(diagnosticsPool, 'GMBolusVelocityZ', GMBolusVelocityZ)
         call mpas_pool_get_array(diagnosticsPool, 'GMBolusVelocityZonal', GMBolusVelocityZonal)
         call mpas_pool_get_array(diagnosticsPool, 'GMBolusVelocityMeridional', GMBolusVelocityMeridional)

         call mpas_pool_get_array(diagnosticsPool, 'relativeSlopeTopOfEdge', relativeSlopeTopOfEdge)
         call mpas_pool_get_array(diagnosticsPool, 'relativeSlopeTopOfCellX', relativeSlopeTopOfCellX)
         call mpas_pool_get_array(diagnosticsPool, 'relativeSlopeTopOfCellY', relativeSlopeTopOfCellY)
         call mpas_pool_get_array(diagnosticsPool, 'relativeSlopeTopOfCellZ', relativeSlopeTopOfCellZ)
         call mpas_pool_get_array(diagnosticsPool, 'relativeSlopeTopOfCellZonal', relativeSlopeTopOfCellZonal)
         call mpas_pool_get_array(diagnosticsPool, 'relativeSlopeTopOfCellMeridional', relativeSlopeTopOfCellMeridional)

         call mpas_pool_get_array(diagnosticsPool, 'gmStreamFuncTopOfEdge', gmStreamFuncTopOfEdge)
         call mpas_pool_get_array(diagnosticsPool, 'GMStreamFuncX', GMStreamFuncX)
         call mpas_pool_get_array(diagnosticsPool, 'GMStreamFuncY', GMStreamFuncY)
         call mpas_pool_get_array(diagnosticsPool, 'GMStreamFuncZ', GMStreamFuncZ)
         call mpas_pool_get_array(diagnosticsPool, 'GMStreamFuncZonal', GMStreamFuncZonal)
         call mpas_pool_get_array(diagnosticsPool, 'GMStreamFuncMeridional', GMStreamFuncMeridional)

         !$omp sections
         !$omp section
         call mpas_reconstruct(meshPool, normalTransportVelocity,          &
                          transportVelocityX,            &
                          transportVelocityY,            &
                          transportVelocityZ,            &
                          transportVelocityZonal,        &
                          transportVelocityMeridional    &
                         )

         !$omp section
         call mpas_reconstruct(meshPool, normalGMBolusVelocity,          &
                          GMBolusVelocityX,            &
                          GMBolusVelocityY,            &
                          GMBolusVelocityZ,            &
                          GMBolusVelocityZonal,        &
                          GMBolusVelocityMeridional    &
                         )

         !$omp section
         call mpas_reconstruct(meshPool, relativeSlopeTopOfEdge,          &
                         relativeSlopeTopOfCellX,            &
                         relativeSlopeTopOfCellY,            &
                         relativeSlopeTopOfCellZ,            &
                         relativeSlopeTopOfCellZonal,        &
                         relativeSlopeTopOfCellMeridional    &
                        )

         !$omp section
         call mpas_reconstruct(meshPool, gmStreamFuncTopOfEdge,          &
                         GMStreamFuncX,            &
                         GMStreamFuncY,            &
                         GMStreamFuncZ,            &
                         GMStreamFuncZonal,        &
                         GMStreamFuncMeridional    &
                        )
         !$omp end sections

   end subroutine ocn_reconstruct_gm_vectors!}}}

!***********************************************************************

end module ocn_diagnostics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
