! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
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
!>  diagnostic variables, and other quantities such as vertTransportVelocityTop.
!
!-----------------------------------------------------------------------

module ocn_diagnostics

   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_timer

   use ocn_gm
   use ocn_equation_of_state

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
             ocn_diagnostics_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: ke_cell_flag, ke_vertex_flag
   real (kind=RKIND) ::  coef_3rd_order, fCoef

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

   subroutine ocn_diagnostic_solve(dt, s, grid, scratch)!{{{

      implicit none

      real (kind=RKIND), intent(in) :: dt !< Input: Time step
      type (state_type), intent(inout) :: s !< Input/Output: State information
      type (mesh_type), intent(in) :: grid !< Input: Grid information
      type (scratch_type), intent(inout) :: scratch !< Input: scratch variables

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j
      integer :: boundaryMask, velMask, nCells, nEdges, nVertices, nVertLevels, vertexDegree, err

      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, &
        maxLevelCell, maxLevelEdgeTop, maxLevelEdgeBot, &
        maxLevelVertexBot
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, &
        verticesOnEdge, edgesOnEdge, edgesOnVertex,boundaryCell, kiteIndexOnCell, &
        verticesOnCell, edgeSignOnVertex, edgeSignOnCell, edgesOnCell

      real (kind=RKIND) :: d2fdx2_cell1, d2fdx2_cell2, coef_3rd_order, r_tmp, &
        invAreaCell1, invAreaCell2, invAreaTri1, invAreaTri2, invLength, layerThicknessVertex, coef

      real (kind=RKIND), dimension(:), allocatable:: pTop, div_hu

      real (kind=RKIND), dimension(:), pointer :: &
        bottomDepth, fVertex, dvEdge, dcEdge, areaCell, areaTriangle, ssh, seaSurfacePressure
      real (kind=RKIND), dimension(:,:), pointer :: &
        weightsOnEdge, kiteAreasOnVertex, layerThicknessEdge, layerThickness, normalVelocity, tangentialVelocity, pressure,&
        circulation, kineticEnergyCell, montgomeryPotential, vertTransportVelocityTop, zMid, divergence, &
        relativeVorticity, relativeVorticityCell, &
        normalizedPlanetaryVorticityEdge, normalizedPlanetaryVorticityVertex, &
        normalizedRelativeVorticityEdge, normalizedRelativeVorticityVertex, normalizedRelativeVorticityCell, &
        density, displacedDensity, potentialDensity, temperature, salinity, kineticEnergyVertex, kineticEnergyVertexOnCells, uBolusGM, uTransport, &
        vertVelocityTop, BruntVaisalaFreqTop, &
        vorticityGradientNormalComponent, vorticityGradientTangentialComponent
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers, deriv_two
      character :: c1*6

      layerThickness => s % layerThickness % array
      normalVelocity => s % normalVelocity % array
      uTransport  => s % uTransport % array
      uBolusGM    => s % uBolusGM % array
      tangentialVelocity => s % tangentialVelocity % array
      layerThicknessEdge => s % layerThicknessEdge % array
      circulation => s % circulation % array
      divergence  => s % divergence % array
      kineticEnergyCell      => s % kineticEnergyCell % array
      relativeVorticity      => s % relativeVorticity % array
      relativeVorticityCell  => s % relativeVorticityCell % array
      normalizedPlanetaryVorticityEdge   => s % normalizedPlanetaryVorticityEdge % array
      normalizedRelativeVorticityEdge    => s % normalizedRelativeVorticityEdge % array
      normalizedRelativeVorticityCell    => s % normalizedRelativeVorticityCell % array
      density         => s % density % array
      displacedDensity=> s % displacedDensity % array
      potentialDensity=> s % potentialDensity % array
      montgomeryPotential     => s % montgomeryPotential % array
      pressure    => s % pressure % array
      zMid        => s % zMid % array
      ssh         => s % ssh % array
      tracers     => s % tracers % array
      vertVelocityTop => s % vertVelocityTop % array
      BruntVaisalaFreqTop => s % BruntVaisalaFreqTop % array

      weightsOnEdge     => grid % weightsOnEdge % array
      kiteAreasOnVertex => grid % kiteAreasOnVertex % array
      cellsOnEdge       => grid % cellsOnEdge % array
      cellsOnVertex     => grid % cellsOnVertex % array
      verticesOnEdge    => grid % verticesOnEdge % array
      nEdgesOnCell      => grid % nEdgesOnCell % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnCell       => grid % edgesOnCell % array
      edgesOnEdge       => grid % edgesOnEdge % array
      edgesOnVertex     => grid % edgesOnVertex % array
      dcEdge            => grid % dcEdge % array
      dvEdge            => grid % dvEdge % array
      areaCell          => grid % areaCell % array
      areaTriangle      => grid % areaTriangle % array
      bottomDepth       => grid % bottomDepth % array
      fVertex           => grid % fVertex % array
      deriv_two         => grid % deriv_two % array
      maxLevelCell      => grid % maxLevelCell % array
      maxLevelEdgeTop   => grid % maxLevelEdgeTop % array
      maxLevelEdgeBot   => grid % maxLevelEdgeBot % array
      maxLevelVertexBot => grid % maxLevelVertexBot % array
      kiteIndexOnCell => grid % kiteIndexOnCell % array
      verticesOnCell => grid % verticesOnCell % array

      seaSurfacePressure => grid % seaSurfacePressure % array
                  
      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels
      vertexDegree = grid % vertexDegree

      boundaryCell => grid % boundaryCell % array

      edgeSignOnVertex => grid % edgeSignOnVertex % array
      edgeSignOnCell => grid % edgeSignOnCell % array

      !
      ! Compute height on cell edges at velocity locations
      !   Namelist options control the order of accuracy of the reconstructed layerThicknessEdge value
      !

      ! initialize layerThicknessEdge to avoid divide by zero and NaN problems.
      layerThicknessEdge = -1.0e34
      coef_3rd_order = config_coef_3rd_order

      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,maxLevelEdgeTop(iEdge)
            layerThicknessEdge(k,iEdge) = 0.5 * (layerThickness(k,cell1) + layerThickness(k,cell2))
         end do
      end do

      !
      ! set the velocity and height at dummy address
      !    used -1e34 so error clearly occurs if these values are used.
      !
      normalVelocity(:,nEdges+1) = -1e34
      layerThickness(:,nCells+1) = -1e34
      tracers(s % index_temperature,:,nCells+1) = -1e34
      tracers(s % index_salinity,:,nCells+1) = -1e34

      circulation(:,:) = 0.0
      relativeVorticity(:,:) = 0.0
      divergence(:,:) = 0.0
      vertVelocityTop(:,:)=0.0
      kineticEnergyCell(:,:) = 0.0
      tangentialVelocity(:,:) = 0.0
      do iVertex = 1, nVertices
         invAreaTri1 = 1.0 / areaTriangle(iVertex)
         do i = 1, vertexDegree
            iEdge = edgesOnVertex(i, iVertex)
            do k = 1, maxLevelVertexBot(iVertex)
              r_tmp = dcEdge(iEdge) * normalVelocity(k, iEdge)

              circulation(k, iVertex) = circulation(k, iVertex) + edgeSignOnVertex(i, iVertex) * r_tmp 
              relativeVorticity(k, iVertex) = relativeVorticity(k, iVertex) + edgeSignOnVertex(i, iVertex) * r_tmp * invAreaTri1
            end do
         end do
      end do

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

      allocate(div_hu(nVertLevels))
      do iCell = 1, nCells
         div_hu(:) = 0.0
         invAreaCell1 = 1.0 / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelCell(iCell)
               r_tmp = dvEdge(iEdge) * normalVelocity(k, iEdge) * invAreaCell1

               divergence(k, iCell) = divergence(k, iCell) - edgeSignOnCell(i, iCell) * r_tmp
               div_hu(k)    = div_hu(k) - layerThicknessEdge(k, iEdge) * edgeSignOnCell(i, iCell) * r_tmp 
               kineticEnergyCell(k, iCell) = kineticEnergyCell(k, iCell) + 0.25 * r_tmp * dcEdge(iEdge) * normalVelocity(k,iEdge)
            end do
         end do
         ! Vertical velocity at bottom (maxLevelCell(iCell)+1) is zero, initialized above.
         do k=maxLevelCell(iCell),1,-1
            vertVelocityTop(k,iCell) = vertVelocityTop(k+1,iCell) - div_hu(k)
         end do         
      end do
      deallocate(div_hu)

      do iEdge=1,nEdges
         ! Compute v (tangential) velocities
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            do k = 1,maxLevelEdgeTop(iEdge) 
               tangentialVelocity(k,iEdge) = tangentialVelocity(k,iEdge) + weightsOnEdge(i,iEdge) * normalVelocity(k, eoe)
            end do
         end do
      end do

      !
      ! Compute kinetic energy
      !
      call mpas_allocate_scratch_field(scratch % kineticEnergyVertex, .true.)
      call mpas_allocate_scratch_field(scratch % kineticEnergyVertexOnCells, .true.)
      kineticEnergyVertex         => scratch % kineticEnergyVertex % array
      kineticEnergyVertexOnCells  => scratch % kineticEnergyVertexOnCells % array
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
      do iCell=1,nCells*ke_vertex_flag
         do k=1,nVertLevels
            kineticEnergyCell(k,iCell) = 5.0/8.0*kineticEnergyCell(k,iCell) + 3.0/8.0*kineticEnergyVertexOnCells(k,iCell)
         end do
      end do

      call mpas_deallocate_scratch_field(scratch % kineticEnergyVertex, .true.)
      call mpas_deallocate_scratch_field(scratch % kineticEnergyVertexOnCells, .true.)


      !
      ! Compute normalized relative and planetary vorticity
      !
      call mpas_allocate_scratch_field(scratch % normalizedRelativeVorticityVertex, .true.)
      call mpas_allocate_scratch_field(scratch % normalizedPlanetaryVorticityVertex, .true.)
      normalizedPlanetaryVorticityVertex  => scratch % normalizedPlanetaryVorticityVertex % array
      normalizedRelativeVorticityVertex  => scratch % normalizedRelativeVorticityVertex % array
      do iVertex = 1,nVertices
         invAreaTri1 = 1.0 / areaTriangle(iVertex)
         do k=1,maxLevelVertexBot(iVertex)
            layerThicknessVertex = 0.0
            do i=1,vertexDegree
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

         call mpas_allocate_scratch_field(scratch % vorticityGradientNormalComponent, .true.)
         call mpas_allocate_scratch_field(scratch % vorticityGradientTangentialComponent, .true.)
         vorticityGradientNormalComponent => scratch % vorticityGradientNormalComponent % array
         vorticityGradientTangentialComponent => scratch % vorticityGradientTangentialComponent % array
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
         call mpas_deallocate_scratch_field(scratch % vorticityGradientNormalComponent, .true.)
         call mpas_deallocate_scratch_field(scratch % vorticityGradientTangentialComponent, .true.)

      endif
      call mpas_deallocate_scratch_field(scratch % normalizedRelativeVorticityVertex, .true.)
      call mpas_deallocate_scratch_field(scratch % normalizedPlanetaryVorticityVertex, .true.)

      !
      ! equation of state
      !
      ! For an isopycnal model, density should remain constant.
      ! For zlevel, calculate in-situ density
      if (config_vert_coord_movement.ne.'isopycnal') then
         call mpas_timer_start("equation of state", .false., diagEOSTimer)

         ! compute in-place density
         call ocn_equation_of_state_density(s, grid, 0, 'relative', density, err)

         ! compute potentialDensity, the density displaced adiabatically to the mid-depth of top layer.
         call ocn_equation_of_state_density(s, grid, 1, 'absolute', potentialDensity, err)

         ! compute displacedDensity, density displaced adiabatically to the mid-depth one layer deeper.  
         ! That is, layer k has been displaced to the depth of layer k+1.
         call ocn_equation_of_state_density(s, grid, 1, 'relative', displacedDensity, err)

         call mpas_timer_stop("equation of state", diagEOSTimer)
      endif

      !
      ! Pressure
      ! This section must be placed in the code after computing the density.
      !
      if (config_pressure_gradient_type.eq.'MontgomeryPotential') then

        ! For Isopycnal model.
        ! Compute pressure at top of each layer, and then
        ! Montgomery Potential.
        allocate(pTop(nVertLevels))
        do iCell=1,nCells

           ! assume atmospheric pressure at the surface is zero for now.
           pTop(1) = 0.0
           ! For isopycnal mode, p is the Montgomery Potential.
           ! At top layer it is g*SSH, where SSH may be off by a 
           ! constant (ie, bottomDepth can be relative to top or bottom)
           montgomeryPotential(1,iCell) = gravity &
              * (bottomDepth(iCell) + sum(layerThickness(1:nVertLevels,iCell)))

           do k=2,nVertLevels
              pTop(k) = pTop(k-1) + density(k-1,iCell)*gravity* layerThickness(k-1,iCell)

              ! from delta M = p delta / density
              montgomeryPotential(k,iCell) = montgomeryPotential(k-1,iCell) &
                 + pTop(k)*(1.0/density(k,iCell) - 1.0/density(k-1,iCell)) 
           end do

        end do
        deallocate(pTop)

      else

        do iCell=1,nCells
           ! Pressure for generalized coordinates.
           ! Pressure at top surface may be due to atmospheric pressure
           ! or an ice-shelf depression. 
           pressure(1,iCell) = seaSurfacePressure(iCell) + density(1,iCell)*gravity &
              * 0.5*layerThickness(1,iCell)

           do k=2,maxLevelCell(iCell)
              pressure(k,iCell) = pressure(k-1,iCell)  &
                + 0.5*gravity*(  density(k-1,iCell)*layerThickness(k-1,iCell) &
                               + density(k  ,iCell)*layerThickness(k  ,iCell))
           end do

           ! Compute zMid, the z-coordinate of the middle of the layer.
           ! This is used for the density g grad z momentum term.
           ! Note the negative sign, since bottomDepth is positive
           ! and z-coordinates are negative below the surface.
           k = maxLevelCell(iCell)
           zMid(k:nVertLevels,iCell) = -bottomDepth(iCell) + 0.5*layerThickness(k,iCell)

           do k=maxLevelCell(iCell)-1, 1, -1
              zMid(k,iCell) = zMid(k+1,iCell)  &
                + 0.5*(  layerThickness(k+1,iCell) &
                       + layerThickness(k  ,iCell))
           end do

        end do

      endif

      !
      ! Brunt-Vaisala frequency
      !
      coef = -gravity/config_density0
      do iCell=1,nCells
         BruntVaisalaFreqTop(1,iCell) = 0.0
         do k=2,maxLevelCell(iCell)
            BruntVaisalaFreqTop(k,iCell) = coef * (displacedDensity(k-1,iCell) - density(k,iCell)) & 
              / (zMid(k-1,iCell) - zMid(k,iCell))
          end do
      end do

      !
      ! Sea Surface Height
      !
      do iCell=1,nCells
         ! Start at the bottom where we know the depth, and go up.
         ! The bottom depth for this cell is bottomDepth(iCell).
         ! Note the negative sign, since bottomDepth is positive
         ! and z-coordinates are negative below the surface.

         ssh(iCell) = - bottomDepth(iCell) + sum(layerThickness(1:maxLevelCell(iCell),iCell))

      end do

      !
      ! Apply the GM closure as a bolus velocity
      !
      if (config_h_kappa .GE. epsilon(0D0)) then
         call ocn_gm_compute_uBolus(s,grid)
      else
         uBolusGM = 0.0
      end if

   end subroutine ocn_diagnostic_solve!}}}

!***********************************************************************
!
!  routine ocn_vert_transport_velocity_top
!
!> \brief   Computes vertical transport
!> \author  Mark Petersen
!> \date    23 September 2011
!> \details 
!>  This routine computes the vertical transport through the top of each 
!>  cell.
!
!-----------------------------------------------------------------------
   subroutine ocn_vert_transport_velocity_top(grid,layerThickness,layerThicknessEdge,normalVelocity,vertTransportVelocityTop, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness    !< Input: thickness

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdge     !< Input: layerThickness interpolated to an edge

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity     !< Input: transport

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         vertTransportVelocityTop     !< Output: vertical transport at top of cell

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j, cov
      real (kind=RKIND) :: flux, layerThicknessVertex, workpv, density0Inv, thicknessSum, invAreaCell

      integer :: nCells, nEdges, nVertices, nVertLevels, vertexDegree


      real (kind=RKIND), dimension(:), pointer :: &
        dvEdge, areaCell, vertCoordMovementWeights
      real (kind=RKIND), dimension(:), allocatable:: div_hu, h_tend_col
      real (kind=RKIND) :: div_hu_btr

      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, &
        verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex, &
        boundaryEdge, boundaryCell, edgeSignOnCell
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, &
        maxLevelCell, maxLevelEdgeTop, maxLevelEdgeBot, &
        maxLevelVertexBot,  maxLevelVertexTop

      err = 0

      nEdgesOnCell      => grid % nEdgesOnCell % array
      areaCell          => grid % areaCell % array
      cellsOnEdge       => grid % cellsOnEdge % array
      edgesOnCell       => grid % edgesOnCell % array
      edgeSignOnCell    => grid % edgeSignOnCell % array
      maxLevelCell      => grid % maxLevelCell % array
      maxLevelEdgeBot   => grid % maxLevelEdgeBot % array
      dvEdge            => grid % dvEdge % array
      vertCoordMovementWeights => grid % vertCoordMovementWeights % array

      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertLevels = grid % nVertLevels


      if (config_vert_coord_movement.eq.'isopycnal') then
        ! set vertical transport to zero in isopycnal case
        vertTransportVelocityTop=0.0
        return
      end if

      allocate(div_hu(nVertLevels), h_tend_col(nVertLevels))

      !
      ! Compute div(h^{edge} u) for each cell
      ! See Ringler et al. (2010) jcp paper, eqn 19, 21, and fig. 3.
      !

      do iCell=1,nCells
        div_hu(:) = 0.0
        div_hu_btr = 0.0
        thicknessSum = 0.0
        invAreaCell = 1.0 / areaCell(iCell)

        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)

          do k = 1, maxLevelEdgeBot(iEdge)
            flux = layerThicknessEdge(k, iEdge) * normalVelocity(k, iEdge) * dvEdge(iEdge) * edgeSignOnCell(i, iCell) * invAreaCell
            div_hu(k) = div_hu(k) - flux
            div_hu_btr = div_hu_btr - flux
          end do
        end do

        do k = 1, maxLevelCell(iCell)
           h_tend_col(k) = - vertCoordMovementWeights(k) * layerThickness(k, iCell) * div_hu_btr
           thicknessSum = thicknessSum + vertCoordMovementWeights(k) * layerThickness(k, iCell)
        end do

        if(thicknessSum > 0.0) then
           h_tend_col = h_tend_col / thicknessSum
        end if

        ! Vertical transport through layer interface at top and bottom is zero.
        vertTransportVelocityTop(1,iCell) = 0.0
        vertTransportVelocityTop(maxLevelCell(iCell)+1,iCell) = 0.0
        do k=maxLevelCell(iCell),2,-1
           vertTransportVelocityTop(k,iCell) = vertTransportVelocityTop(k+1,iCell) - div_hu(k) - h_tend_col(k)
        end do
      end do

      deallocate(div_hu, h_tend_col)

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

   subroutine ocn_fuperp(s, grid)!{{{
      implicit none

      type (state_type), intent(inout) :: s !< Input/Output: State information
      type (mesh_type), intent(in) :: grid !< Input: Grid information

      integer :: iEdge, cell1, cell2, eoe, i, j, k
      integer :: nEdgesSolve
      real (kind=RKIND), dimension(:), pointer :: fEdge
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge, normalVelocity, normalBaroclinicVelocity
      type (dm_info) :: dminfo

      integer, dimension(:), pointer :: maxLevelEdgeTop, nEdgesOnEdge
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnEdge

      call mpas_timer_start("ocn_fuperp")

      normalVelocity           => s % normalVelocity % array
      normalBaroclinicVelocity        => s % normalBaroclinicVelocity % array
      weightsOnEdge     => grid % weightsOnEdge % array
      fEdge             => grid % fEdge % array
      maxLevelEdgeTop      => grid % maxLevelEdgeTop % array
      cellsOnEdge       => grid % cellsOnEdge % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array

      fEdge       => grid % fEdge % array

      nEdgesSolve = grid % nEdgesSolve

      !
      ! Put f*normalBaroclinicVelocity^{perp} in u as a work variable
      !
      do iEdge=1,nEdgesSolve
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k=1,maxLevelEdgeTop(iEdge)

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
   subroutine ocn_filter_btr_mode_vel(s, grid)!{{{
      implicit none

      type (state_type), intent(inout) :: s
      type (mesh_type), intent(in) :: grid

      integer :: iEdge, k, nEdges
      real (kind=RKIND) :: vertSum, normalThicknessFluxSum, thicknessSum
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge, normalVelocity
      integer, dimension(:), pointer :: maxLevelEdgeTop

      call mpas_timer_start("ocn_filter_btr_mode_vel")

      normalVelocity => s % normalVelocity % array
      layerThicknessEdge => s % layerThicknessEdge % array
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      nEdges      = grid % nEdges

      do iEdge=1,nEdges

        ! thicknessSum is initialized outside the loop because on land boundaries 
        ! maxLevelEdgeTop=0, but I want to initialize thicknessSum with a 
        ! nonzero value to avoid a NaN.
        normalThicknessFluxSum = layerThicknessEdge(1,iEdge) * normalVelocity(1,iEdge)
        thicknessSum  = layerThicknessEdge(1,iEdge)

        do k=2,maxLevelEdgeTop(iEdge)
          normalThicknessFluxSum = normalThicknessFluxSum + layerThicknessEdge(k,iEdge) * normalVelocity(k,iEdge)
          thicknessSum  =  thicknessSum + layerThicknessEdge(k,iEdge)
        enddo

        vertSum = normalThicknessFluxSum/thicknessSum
        do k=1,maxLevelEdgeTop(iEdge)
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
   subroutine ocn_filter_btr_mode_tend_vel(tend, s, grid)!{{{
      implicit none

      type (tend_type), intent(inout) :: tend
      type (state_type), intent(in) :: s
      type (mesh_type), intent(in) :: grid

      integer :: iEdge, k, nEdges
      real (kind=RKIND) :: vertSum, normalThicknessFluxSum, thicknessSum
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge, tend_normalVelocity

      integer, dimension(:), pointer :: maxLevelEdgeTop

      call mpas_timer_start("ocn_filter_btr_mode_tend_vel")

      tend_normalVelocity      => tend % normalVelocity % array
      layerThicknessEdge      => s % layerThicknessEdge % array
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      nEdges      = grid % nEdges

      do iEdge=1,nEdges

        ! thicknessSum is initialized outside the loop because on land boundaries 
        ! maxLevelEdgeTop=0, but I want to initialize thicknessSum with a 
        ! nonzero value to avoid a NaN.
        normalThicknessFluxSum = layerThicknessEdge(1,iEdge) * tend_normalVelocity(1,iEdge)
        thicknessSum  = layerThicknessEdge(1,iEdge)

        do k=2,maxLevelEdgeTop(iEdge)
          normalThicknessFluxSum = normalThicknessFluxSum + layerThicknessEdge(k,iEdge) * tend_normalVelocity(k,iEdge)
          thicknessSum  =  thicknessSum + layerThicknessEdge(k,iEdge)
        enddo

        vertSum = normalThicknessFluxSum/thicknessSum
        do k=1,maxLevelEdgeTop(iEdge)
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
          .or.trim(config_time_integrator) == 'unsplit_explicit') then
            ! For split explicit, PV is eta/h because the Coriolis term 
            ! is added separately to the momentum tendencies.
            fCoef = 0
        end if

    end subroutine ocn_diagnostics_init!}}}

!***********************************************************************

end module ocn_diagnostics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
