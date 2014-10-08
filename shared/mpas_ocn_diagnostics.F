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

   use mpas_grid_types
   use mpas_constants
   use mpas_timer
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

   subroutine ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, timeLevelIn)!{{{

      real (kind=RKIND), intent(in) :: dt !< Input: Time step
      type (mpas_pool_type), intent(in) :: statePool !< Input: State information
      type (mpas_pool_type), intent(in) :: forcingPool !< Input: Forcing information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      type (mpas_pool_type), intent(inout) :: diagnosticsPool  !< Input: diagnostic fields derived from State
      type (mpas_pool_type), intent(in) :: scratchPool !< Input: scratch variables
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
        invAreaCell1, invAreaCell2, invAreaTri1, invAreaTri2, invLength, layerThicknessVertex, coef

      real (kind=RKIND), dimension(:), allocatable:: pTop, div_hu,div_huTransport,div_huGMBolus

      real (kind=RKIND), dimension(:), pointer :: &
        bottomDepth, fVertex, dvEdge, dcEdge, areaCell, areaTriangle, ssh, seaSurfacePressure
      real (kind=RKIND), dimension(:,:), pointer :: &
        weightsOnEdge, kiteAreasOnVertex, layerThicknessEdge, layerThickness, normalVelocity, normalTransportVelocity, normalGMBolusVelocity, tangentialVelocity, pressure,&
        circulation, kineticEnergyCell, montgomeryPotential, vertAleTransportTop, zMid, zTop, divergence, &
        relativeVorticity, relativeVorticityCell, &
        normalizedPlanetaryVorticityEdge, normalizedPlanetaryVorticityVertex, &
        normalizedRelativeVorticityEdge, normalizedRelativeVorticityVertex, normalizedRelativeVorticityCell, &
        density, displacedDensity, potentialDensity, temperature, salinity, kineticEnergyVertex, kineticEnergyVertexOnCells, &
        vertVelocityTop, vertTransportVelocityTop, vertGMBolusVelocityTop, BruntVaisalaFreqTop, &
        vorticityGradientNormalComponent, vorticityGradientTangentialComponent, gradSSH, &
        inSituThermalExpansionCoeff, inSituSalineContractionCoeff
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers, derivTwo
      character :: c1*6

      real (kind=RKIND), dimension(:,:), pointer :: tracersSurfaceValue

      type (field2DReal), pointer :: kineticEnergyVertexField, kineticEnergyVertexOnCellsField
      type (field2DReal), pointer :: normalizedRelativeVorticityVertexField, normalizedPlanetaryVorticityVertexField
      type (field2DReal), pointer :: vorticityGradientNormalComponentField, vorticityGradientTangentialComponentField

      integer :: timeLevel
      integer, pointer :: indexTemperature, indexSalinity
      real (kind=RKIND), pointer :: config_density0, config_apvm_scale_factor,  config_coef_3rd_order
      character (len=StrKIND), pointer :: config_pressure_gradient_type

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_config(ocnConfigs, 'config_density0', config_density0)
      call mpas_pool_get_config(ocnConfigs, 'config_apvm_scale_factor', config_apvm_scale_factor)
      call mpas_pool_get_config(ocnConfigs, 'config_pressure_gradient_type', config_pressure_gradient_type)
      call mpas_pool_get_config(ocnConfigs, 'config_coef_3rd_order', config_coef_3rd_order)

      call mpas_pool_get_dimension(statePool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(statePool, 'index_salinity', indexSalinity)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'tracers', tracers, timeLevel)
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
                  
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

      call mpas_pool_get_array(diagnosticsPool, 'tracersSurfaceValue', tracersSurfaceValue)

      !
      ! Compute height on cell edges at velocity locations
      !   Namelist options control the order of accuracy of the reconstructed layerThicknessEdge value
      !

      ! initialize layerThicknessEdge to avoid divide by zero and NaN problems.
      layerThicknessEdge = -1.0e34
      coef_3rd_order = config_coef_3rd_order 

      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = 1, maxLevelEdgeTop(iEdge)
            layerThicknessEdge(k,iEdge) = 0.5 * (layerThickness(k,cell1) + layerThickness(k,cell2))
         end do
      end do

      !
      ! set the velocity and height at dummy address
      !    used -1e34 so error clearly occurs if these values are used.
      !
      normalVelocity(:,nEdges+1) = -1e34
      layerThickness(:,nCells+1) = -1e34
      tracers(indexTemperature,:,nCells+1) = -1e34
      tracers(indexSalinity,:,nCells+1) = -1e34

      divergence(:,:) = 0.0
      vertVelocityTop(:,:)=0.0
      kineticEnergyCell(:,:) = 0.0
      tangentialVelocity(:,:) = 0.0

      call ocn_relativeVorticity_circulation(relativeVorticity, circulation, meshPool, normalVelocity, err)

      relativeVorticityCell(:,:) = 0.0
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

      !
      ! Compute divergence, kinetic energy, and vertical velocity
      !
      allocate(div_hu(nVertLevels),div_huTransport(nVertLevels),div_huGMBolus(nVertLevels))
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
      deallocate(div_hu,div_huTransport,div_huGMBolus)

      do iEdge = 1, nEdges
         ! Compute v (tangential) velocities
         do i = 1, nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            do k = 1, maxLevelEdgeTop(iEdge) 
               tangentialVelocity(k,iEdge) = tangentialVelocity(k,iEdge) + weightsOnEdge(i,iEdge) * normalVelocity(k, eoe)
            end do
         end do
      end do

      !
      ! Compute kinetic energy
      !
      call mpas_pool_get_field(scratchPool, 'kineticEnergyVertex', kineticEnergyVertexField)
      call mpas_pool_get_field(scratchPool, 'kineticEnergyVertexOnCells', kineticEnergyVertexOnCellsField)
      call mpas_allocate_scratch_field(kineticEnergyVertexField, .true.)
      call mpas_allocate_scratch_field(kineticEnergyVertexOnCellsField, .true.)
      kineticEnergyVertex         => kineticEnergyVertexField % array
      kineticEnergyVertexOnCells  => kineticEnergyVertexOnCellsField % array
      kineticEnergyVertex(:,:) = 0.0; 
      kineticEnergyVertexOnCells(:,:) = 0.0
      do iVertex = 1, nVertices*ke_vertex_flag
        do i = 1, vertexDegree
          iEdge = edgesOnVertex(i, iVertex)
          r_tmp = dcEdge(iEdge) * dvEdge(iEdge) * 0.25 / areaTriangle(iVertex)
          do k = 1, nVertLevels
            kineticEnergyVertex(k, iVertex) = kineticEnergyVertex(k, iVertex) + r_tmp * normalVelocity(k, iEdge)**2
          end do
        end do
      end do

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

      !
      ! Compute kinetic energy in each cell by blending kineticEnergyCell and kineticEnergyVertexOnCells
      !
      do iCell = 1, nCells * ke_vertex_flag
         do k = 1, nVertLevels
            kineticEnergyCell(k,iCell) = 5.0 / 8.0 * kineticEnergyCell(k,iCell) + 3.0 / 8.0 * kineticEnergyVertexOnCells(k,iCell)
         end do
      end do

      call mpas_deallocate_scratch_field(kineticEnergyVertexField, .true.)
      call mpas_deallocate_scratch_field(kineticEnergyVertexOnCellsField, .true.)

      !
      ! Compute normalized relative and planetary vorticity
      !
      call mpas_pool_get_field(scratchPool, 'normalizedRelativeVorticityVertex', normalizedRelativeVorticityVertexField)
      call mpas_pool_get_field(scratchPool, 'normalizedPlanetaryVorticityVertex', normalizedPlanetaryVorticityVertexField)
      call mpas_allocate_scratch_field(normalizedRelativeVorticityVertexField, .true.)
      call mpas_allocate_scratch_field(normalizedPlanetaryVorticityVertexField, .true.)
      normalizedPlanetaryVorticityVertex  => normalizedPlanetaryVorticityVertexField % array
      normalizedRelativeVorticityVertex  => normalizedRelativeVorticityVertexField % array
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

      normalizedRelativeVorticityEdge(:,:) = 0.0
      normalizedPlanetaryVorticityEdge(:,:) = 0.0
      do iEdge = 1, nEdges
        vertex1 = verticesOnEdge(1, iEdge)
        vertex2 = verticesOnEdge(2, iEdge)
        do k = 1, maxLevelEdgeBot(iEdge)
          normalizedRelativeVorticityEdge(k, iEdge) = 0.5 * (normalizedRelativeVorticityVertex(k, vertex1) + normalizedRelativeVorticityVertex(k, vertex2))
          normalizedPlanetaryVorticityEdge(k, iEdge) = 0.5 * (normalizedPlanetaryVorticityVertex(k, vertex1) + normalizedPlanetaryVorticityVertex(k, vertex2))
        end do
      end do

      normalizedRelativeVorticityCell(:,:) = 0.0
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

      ! Diagnostics required for the Anticipated Potential Vorticity Method (apvm).
      if (config_apvm_scale_factor>1e-10) then

         call mpas_pool_get_field(scratchPool, 'vorticityGradientNormalComponent', vorticityGradientNormalComponentField)
         call mpas_pool_get_field(scratchPool, 'vorticityGradientTangentialComponent', vorticityGradientTangentialComponentField)
         call mpas_allocate_scratch_field(vorticityGradientNormalComponentField, .true.)
         call mpas_allocate_scratch_field(vorticityGradientTangentialComponentField, .true.)
         vorticityGradientNormalComponent => vorticityGradientNormalComponentField % array
         vorticityGradientTangentialComponent => vorticityGradientTangentialComponentField % array

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

         !
         ! Modify PV edge with upstream bias.
         !
         do iEdge = 1,nEdges
            do k = 1,maxLevelEdgeBot(iEdge)
              normalizedRelativeVorticityEdge(k,iEdge) = normalizedRelativeVorticityEdge(k,iEdge) &
                - config_apvm_scale_factor * dt * &
                    (  normalVelocity(k,iEdge)     * vorticityGradientNormalComponent(k,iEdge)      &
                     + tangentialVelocity(k,iEdge) * vorticityGradientTangentialComponent(k,iEdge) )
            enddo
         enddo
         call mpas_deallocate_scratch_field(vorticityGradientNormalComponentField, .true.)
         call mpas_deallocate_scratch_field(vorticityGradientTangentialComponentField, .true.)

      endif
      call mpas_deallocate_scratch_field(normalizedRelativeVorticityVertexField, .true.)
      call mpas_deallocate_scratch_field(normalizedPlanetaryVorticityVertexField, .true.)

      !
      ! equation of state
      !
      call mpas_timer_start("equation of state", .false., diagEOSTimer)

      ! compute in-place density
      if (config_pressure_gradient_type.eq.'common_level_eos') then
         ! only compute EOS derivatives if needed.
         call mpas_pool_get_array(diagnosticsPool, 'inSituThermalExpansionCoeff',inSituThermalExpansionCoeff)
         call mpas_pool_get_array(diagnosticsPool, 'inSituSalineContractionCoeff', inSituSalineContractionCoeff)
         call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, 0, 'relative', density, err, &
              inSituThermalExpansionCoeff, inSituSalineContractionCoeff, timeLevelIn=timeLevel)
      else
         call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, 0, 'relative', density, err, &
              timeLevelIn=timeLevel)
      endif

      ! compute potentialDensity, the density displaced adiabatically to the mid-depth of top layer.
      call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, 1, 'absolute', potentialDensity, err, timeLevelIn=timeLevel)

      ! compute displacedDensity, density displaced adiabatically to the mid-depth one layer deeper.  
      ! That is, layer k has been displaced to the depth of layer k+1.
      call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, 1, 'relative', displacedDensity, err, timeLevelIn=timeLevel)

      call mpas_timer_stop("equation of state", diagEOSTimer)

      !
      ! Pressure
      ! This section must be placed in the code after computing the density.
      !
      if (config_pressure_gradient_type.eq.'MontgomeryPotential') then

        ! use Montgomery Potential when layers are isopycnal.
        ! However, one may use 'pressure_and_zmid' when layers are isopycnal as well.
        ! Compute pressure at top of each layer, and then Montgomery Potential.
        allocate(pTop(nVertLevels))
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
        deallocate(pTop)

      else

        do iCell = 1, nCells
           ! Pressure for generalized coordinates.
           ! Pressure at top surface may be due to atmospheric pressure
           ! or an ice-shelf depression. 
           pressure(1,iCell) = seaSurfacePressure(iCell) + density(1,iCell)*gravity &
              * 0.5*layerThickness(1,iCell)

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

      endif

      !
      ! Brunt-Vaisala frequency
      !
      coef = -gravity / config_density0
      do iCell = 1, nCells
         BruntVaisalaFreqTop(1,iCell) = 0.0
         do k = 2, maxLevelCell(iCell)
            BruntVaisalaFreqTop(k,iCell) = coef * (displacedDensity(k-1,iCell) - density(k,iCell)) & 
              / (zMid(k-1,iCell) - zMid(k,iCell))
          end do
      end do

      !
      ! extrapolate tracer values to ocean surface
      ! this eventually be a modelled process
      ! at present, just copy k=1 tracer values onto surface values
      tracersSurfaceValue(:,:) = tracers(:,1,:)

      !
      !  compute fields used as intent(in) to CVMix/KPP
      call computeKPPInputFields(statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, timeLevel)

#ifdef MPAS_CESM
      do iEdge = 1, nEdgesSolve
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)

         gradSSH(1, iEdge) = (ssh(cell2) - ssh(cell1)) / dcEdge(iEdge)
      end do
#endif

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
   subroutine ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, oldLayerThickness, layerThicknessEdge, &
     normalVelocity, oldSSH, newHighFreqThickness, dt, vertAleTransportTop, err)!{{{

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
         layerThicknessEdge     !< Input: layerThickness interpolated to an edge

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity     !< Input: transport

      real (kind=RKIND), dimension(:), intent(in) :: &
         oldSSH     !< Input: sea surface height at old time

      real (kind=RKIND), dimension(:,:), intent(in) :: &
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

      real (kind=RKIND) :: flux, invAreaCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real (kind=RKIND), dimension(:), allocatable :: &
         div_hu_btr       !> barotropic divergence of (thickness*velocity)
      real (kind=RKIND), dimension(:,:), allocatable :: &
         ALE_Thickness, & !> ALE thickness at new time
         div_hu           !> divergence of (thickness*velocity)

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

      allocate(div_hu(nVertLevels,nCells), div_hu_btr(nCells), ALE_Thickness(nVertLevels,nCells))

      !
      ! thickness-weighted divergence and barotropic divergence
      !
      ! See Ringler et al. (2010) jcp paper, eqn 19, 21, and fig. 3.
      do iCell = 1, nCells
         div_hu(:,iCell) = 0.0
         div_hu_btr(iCell) = 0.0
         invAreaCell = 1.0 / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)

            do k = 1, maxLevelEdgeBot(iEdge)
               flux = layerThicknessEdge(k, iEdge) * normalVelocity(k, iEdge) * dvEdge(iEdge) * edgeSignOnCell(i, iCell) * invAreaCell
               div_hu(k,iCell) = div_hu(k,iCell) - flux
               div_hu_btr(iCell) = div_hu_btr(iCell) - flux
            end do
         end do

      enddo

      !
      ! Compute desired thickness at new time
      !
      call ocn_ALE_thickness(meshPool, verticalMeshPool, oldSSH, div_hu_btr, newHighFreqThickness, dt, ALE_thickness, err)

      !
      ! Vertical transport through layer interfaces
      !
      ! Vertical transport through layer interface at top and bottom is zero.
      ! Here we are using solving the continuity equation for vertAleTransportTop ($w^t$),
      ! and using ALE_Thickness for thickness at the new time.

      do iCell = 1,nCells
         vertAleTransportTop(1,iCell) = 0.0
         vertAleTransportTop(maxLevelCell(iCell)+1,iCell) = 0.0
         do k = maxLevelCell(iCell),2,-1
            vertAleTransportTop(k,iCell) = vertAleTransportTop(k+1,iCell) - div_hu(k,iCell) &
              - (ALE_Thickness(k,iCell) - oldLayerThickness(k,iCell))/dt
         end do
      end do

      deallocate(div_hu, div_hu_btr, ALE_Thickness)

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

      !
      ! Put f*normalBaroclinicVelocity^{perp} in u as a work variable
      !
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
!  routine computeKPPInputFields
!
!> \brief   
!>    Compute fields necessary to drive the CVMix KPP module
!> \author  Todd Ringler
!> \date    20 August 2013
!> \details
!>    CVMix/KPP requires the following fields as intent(in):
!>       buoyancyForcingOBL
!>       surfaceFrictionVelocity
!>       bulkRichardsonNumber
!>
!
!-----------------------------------------------------------------------

    subroutine computeKPPInputFields(statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, timeLevelIn)!{{{

      type (mpas_pool_type), intent(in) :: statePool !< Input/Output: State information
      type (mpas_pool_type), intent(in) :: forcingPool !< Input: Forcing information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      type (mpas_pool_type), intent(inout) :: diagnosticsPool !< Diagnostics information derived from State
      type (mpas_pool_type), intent(in) :: scratchPool !< Input: scratch variables
      integer, intent(in), optional :: timeLevelIn

      ! scalars
      integer, pointer :: nCells, nVertLevels

      ! integer pointers
      integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell

      ! real pointers
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaCell
      real (kind=RKIND), dimension(:), pointer :: penetrativeTemperatureFlux, surfaceMassFlux, surfaceWindStressMagnitude, &
           buoyancyForcingOBL, surfaceFrictionVelocity, boundaryLayerDepth, penetrativeTemperatureFluxOBL
      real (kind=RKIND), dimension(:,:), pointer ::  &
           layerThickness, zMid, zTop, bulkRichardsonNumber, tracersSurfaceValues, densitySurfaceDisplaced, density, &
           normalVelocity, surfaceTracerFlux, thermalExpansionCoeff, salineContractionCoeff

      ! local
      integer :: iCell, iEdge, i, k, err
      integer, pointer :: indexTempFlux, indexSaltFlux
      integer :: timeLevel
      real (kind=RKIND) :: numerator, denominator, factor, deltaVelocitySquared, turbulentVelocitySquared, delU2, invAreaCell

      type (field2DReal), pointer :: densitySurfaceDisplacedField, thermalExpansionCoeffField, salineContractionCoeffField
      real (kind=RKIND), pointer :: config_density0

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_config(ocnConfigs, 'config_density0', config_density0)

      ! set the parameter turbulentVelocitySquared
      turbulentVelocitySquared = 0.001

      ! set scalar values
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(forcingPool, 'index_surfaceTemperatureFlux', indexTempFlux)
      call mpas_pool_get_dimension(forcingPool, 'index_surfaceSalinityFlux', indexSaltFlux)

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
      call mpas_pool_get_array(diagnosticsPool, 'bulkRichardsonNumber', bulkRichardsonNumber)
      call mpas_pool_get_array(diagnosticsPool, 'tracersSurfaceValue ', tracersSurfaceValues)
      call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth', boundaryLayerDepth)
      call mpas_pool_get_array(diagnosticsPool, 'surfaceFrictionVelocity', surfaceFrictionVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'penetrativeTemperatureFluxOBL', penetrativeTemperatureFluxOBL)
      call mpas_pool_get_array(diagnosticsPool, 'buoyancyForcingOBL', buoyancyForcingOBL)

      call mpas_pool_get_array(forcingPool, 'surfaceMassFlux', surfaceMassFlux)
      call mpas_pool_get_array(forcingPool, 'surfaceTracerFlux', surfaceTracerFlux)
      call mpas_pool_get_array(forcingPool, 'penetrativeTemperatureFlux', penetrativeTemperatureFlux)
      call mpas_pool_get_array(forcingPool, 'surfaceWindStressMagnitude', surfaceWindStressMagnitude)

      ! allocate scratch space displaced density computation
      call mpas_pool_get_field(scratchPool, 'densitySurfaceDisplaced', densitySurfaceDisplacedField)
      call mpas_pool_get_field(scratchPool, 'thermalExpansionCoeff', thermalExpansionCoeffField)
      call mpas_pool_get_field(scratchPool, 'salineContractionCoeff', salineContractionCoeffField)
      call mpas_allocate_scratch_field(densitySurfaceDisplacedField, .true.)
      call mpas_allocate_scratch_field(thermalExpansionCoeffField, .true.)
      call mpas_allocate_scratch_field(salineContractionCoeffField, .true.)
      densitySurfaceDisplaced => densitySurfaceDisplacedField % array
      thermalExpansionCoeff => thermalExpansionCoeffField % array
      salineContractionCoeff => salineContractionCoeffField % array

      ! compute EOS by displacing SST/SSS to every vertical layer in column
      call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, 0, 'surfaceDisplaced', densitySurfaceDisplaced, err, &
              thermalExpansionCoeff, salineContractionCoeff, timeLevel)

      ! set value to out-of-bounds
      bulkRichardsonNumber(:,:) = -1.0e34

      do iCell = 1, nCells
       invAreaCell = 1.0 / areaCell(iCell)

       ! compute surface buoyancy forcing based on surface fluxes of mass, temperature, salinity and frazil (frazil to be added later)
       ! since this computation is confusing, variables, units and sign convention is repeated here
       ! everything below should be consistent with that specified in Registry
       ! everything below should be consistent with the CVMix/KPP documentation: https://www.dropbox.com/s/6hqgc0rsoa828nf/cvmix_20aug2013.pdf
       !
       !    surfaceMassFlux: surface mass flux, m/s, positive into ocean
       !    surfaceTracerFlux(indexTempFlux): non-penetrative temperature flux, C m/s, positive into ocean
       !    penetrativeTemperatureFlux: penetrative surface temperature flux at ocean surface, positive into ocean
       !    surfaceTracerFlux(indexSaltFlux): salinity flux, PSU m/s, positive into ocean
       !    penetrativeTemperatureFluxOBL: penetrative temperature flux computed at z=OBL, positive down
       !
       ! note: the following fields used the CVMix/KPP computation of buoyancy forcing are not included here
       !    1. Tm: temperature associated with surfaceMassFlux, C  (here we assume Tm == temperatureSurfaceValue)
       !    2. Sm: salinity associated with surfaceMassFlux, PSU (here we assume Sm == salinitySurfaceValue and account for salinity flux in surfaceTracerFlux array)
       !
         buoyancyForcingOBL(iCell) =  thermalExpansionCoeff (1,iCell) *  &
               (surfaceTracerFlux(indexTempFlux,iCell) + penetrativeTemperatureFlux(iCell) - penetrativeTemperatureFluxOBL(iCell)) &
              - salineContractionCoeff(1,iCell) *  surfaceTracerFlux(indexSaltFlux,iCell)
        
       ! at this point, buoyancyForcingOBL has units of m/s 
       ! change into units of m/s^3 (which can be thought of as units of buoyancy per second)
         buoyancyForcingOBL(iCell) = buoyancyForcingOBL(iCell) * gravity / max(boundaryLayerDepth(iCell),layerThickness(1,iCell))

       ! compute surface friction velocity
         surfaceFrictionVelocity(iCell) = surfacewindStressMagnitude(iCell) / config_density0

       ! loop over vertical to compute bulk Richardson number
       do k = 1, maxLevelCell(iCell)

        ! find deltaVelocitySquared defined at cell centers
        deltaVelocitySquared = 0.0
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          factor = 0.5 * dcEdge(iEdge) * dvEdge(iEdge) * invAreaCell
          delU2 = (normalVelocity(1,iEdge) - normalVelocity(k,iEdge))**2 
          deltaVelocitySquared = deltaVelocitySquared + factor * delU2
        enddo

        numerator = gravity * (zTop(1,iCell) - zMid(k,iCell)) * (density(k,iCell) - densitySurfaceDisplaced(k,iCell))
        denominator = config_density0 * (deltaVelocitySquared + turbulentVelocitySquared)

        ! compute bulk Richardson number
        bulkRichardsonNumber(k,iCell) = numerator / denominator

       enddo
      enddo

      ! deallocate scratch space
      call mpas_deallocate_scratch_field(densitySurfaceDisplacedField, .true.)
      call mpas_deallocate_scratch_field(thermalExpansionCoeffField, .true.)
      call mpas_deallocate_scratch_field(salineContractionCoeffField, .true.)

    end subroutine computeKPPInputFields!}}}

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

   subroutine ocn_reconstruct_gm_vectors(diagnosticsPool, meshPool)

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

         call mpas_reconstruct(meshPool, relativeSlopeTopOfEdge,          &
                         relativeSlopeTopOfCellX,            &
                         relativeSlopeTopOfCellY,            &
                         relativeSlopeTopOfCellZ,            &
                         relativeSlopeTopOfCellZonal,        &
                         relativeSlopeTopOfCellMeridional    &
                        )

         call mpas_reconstruct(meshPool, gmStreamFuncTopOfEdge,          &
                         GMStreamFuncX,            &
                         GMStreamFuncY,            &
                         GMStreamFuncZ,            &
                         GMStreamFuncZonal,        &
                         GMStreamFuncMeridional    &
                        )

   end subroutine ocn_reconstruct_gm_vectors!}}}

!***********************************************************************

end module ocn_diagnostics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
