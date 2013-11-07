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
   use ocn_thick_ale

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

   subroutine ocn_diagnostic_solve(dt, state, forcing, mesh, diagnostics, scratch)!{{{

      real (kind=RKIND), intent(in) :: dt !< Input: Time step
      type (state_type), intent(inout) :: state !< Input/Output: State information
      type (forcing_type), intent(in) :: forcing !< Input: Forcing information
      type (mesh_type), intent(in) :: mesh !< Input: mesh information
      type (diagnostics_type), intent(inout) :: diagnostics  !< Input/Output: diagnostic fields derived from State
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
        circulation, kineticEnergyCell, montgomeryPotential, vertTransportVelocityTop, zMid, zTop, divergence, &
        relativeVorticity, relativeVorticityCell, &
        normalizedPlanetaryVorticityEdge, normalizedPlanetaryVorticityVertex, &
        normalizedRelativeVorticityEdge, normalizedRelativeVorticityVertex, normalizedRelativeVorticityCell, &
        density, displacedDensity, potentialDensity, temperature, salinity, kineticEnergyVertex, kineticEnergyVertexOnCells, uBolusGM, uTransport, &
        vertVelocityTop, BruntVaisalaFreqTop, &
        vorticityGradientNormalComponent, vorticityGradientTangentialComponent
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers, derivTwo
      character :: c1*6

      real (kind=RKIND), dimension(:,:), pointer :: tracersSurfaceValue ! => diagnostics % tracersSurfaceValue % array

      layerThickness => state % layerThickness % array
      normalVelocity => state % normalVelocity % array
      tracers        => state % tracers % array
      ssh            => state % ssh % array

      zMid                             => diagnostics % zMid % array
      zTop                             => diagnostics % zTop % array
      divergence                       => diagnostics % divergence % array
      circulation                      => diagnostics % circulation % array
      relativeVorticity                => diagnostics % relativeVorticity % array
      relativeVorticityCell            => diagnostics % relativeVorticityCell % array
      normalizedPlanetaryVorticityEdge => diagnostics % normalizedPlanetaryVorticityEdge % array
      normalizedRelativeVorticityEdge  => diagnostics % normalizedRelativeVorticityEdge % array
      normalizedRelativeVorticityCell  => diagnostics % normalizedRelativeVorticityCell % array
      density                          => diagnostics % density % array
      displacedDensity                 => diagnostics % displacedDensity % array
      potentialDensity                 => diagnostics % potentialDensity % array
      montgomeryPotential              => diagnostics % montgomeryPotential % array
      pressure                         => diagnostics % pressure % array
      BruntVaisalaFreqTop              => diagnostics % BruntVaisalaFreqTop % array
      tangentialVelocity               => diagnostics % tangentialVelocity % array
      layerThicknessEdge               => diagnostics % layerThicknessEdge % array
      kineticEnergyCell                => diagnostics % kineticEnergyCell % array
      vertVelocityTop                  => diagnostics % vertVelocityTop % array
      uBolusGM                         => diagnostics % uBolusGM % array
      uTransport                       => diagnostics % uTransport % array

      weightsOnEdge     => mesh % weightsOnEdge % array
      kiteAreasOnVertex => mesh % kiteAreasOnVertex % array
      cellsOnEdge       => mesh % cellsOnEdge % array
      cellsOnVertex     => mesh % cellsOnVertex % array
      verticesOnEdge    => mesh % verticesOnEdge % array
      nEdgesOnCell      => mesh % nEdgesOnCell % array
      nEdgesOnEdge      => mesh % nEdgesOnEdge % array
      edgesOnCell       => mesh % edgesOnCell % array
      edgesOnEdge       => mesh % edgesOnEdge % array
      edgesOnVertex     => mesh % edgesOnVertex % array
      dcEdge            => mesh % dcEdge % array
      dvEdge            => mesh % dvEdge % array
      areaCell          => mesh % areaCell % array
      areaTriangle      => mesh % areaTriangle % array
      bottomDepth       => mesh % bottomDepth % array
      fVertex           => mesh % fVertex % array
      derivTwo          => mesh % derivTwo % array
      maxLevelCell      => mesh % maxLevelCell % array
      maxLevelEdgeTop   => mesh % maxLevelEdgeTop % array
      maxLevelEdgeBot   => mesh % maxLevelEdgeBot % array
      maxLevelVertexBot => mesh % maxLevelVertexBot % array
      kiteIndexOnCell   => mesh % kiteIndexOnCell % array
      verticesOnCell    => mesh % verticesOnCell % array

      seaSurfacePressure => forcing % seaSurfacePressure % array
                  
      nCells      = mesh % nCells
      nEdges      = mesh % nEdges
      nVertices   = mesh % nVertices
      nVertLevels = mesh % nVertLevels
      vertexDegree = mesh % vertexDegree

      boundaryCell => mesh % boundaryCell % array

      edgeSignOnVertex => mesh % edgeSignOnVertex % array
      edgeSignOnCell => mesh % edgeSignOnCell % array

      tracersSurfaceValue  => diagnostics % tracersSurfaceValue % array(:,:)

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
      tracers(state % index_temperature,:,nCells+1) = -1e34
      tracers(state % index_salinity,:,nCells+1) = -1e34

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
      call mpas_timer_start("equation of state", .false., diagEOSTimer)

      ! compute in-place density
      call ocn_equation_of_state_density(state, diagnostics, mesh, 0, 'relative', density, err)

      ! compute potentialDensity, the density displaced adiabatically to the mid-depth of top layer.
      call ocn_equation_of_state_density(state, diagnostics, mesh, 1, 'absolute', potentialDensity, err)

      ! compute displacedDensity, density displaced adiabatically to the mid-depth one layer deeper.  
      ! That is, layer k has been displaced to the depth of layer k+1.
      call ocn_equation_of_state_density(state, diagnostics, mesh, 1, 'relative', displacedDensity, err)

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
        do iCell=1,nCells

           ! assume atmospheric pressure at the surface is zero for now.
           pTop(1) = 0.0
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

      elseif (config_pressure_gradient_type.eq.'pressure_and_zmid') then

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
           ! Compute zTop, the z-coordinate of the top of the layer.
           ! Note the negative sign, since bottomDepth is positive
           ! and z-coordinates are negative below the surface.
           k = maxLevelCell(iCell)
           zMid(k:nVertLevels,iCell) = -bottomDepth(iCell) + 0.5*layerThickness(k,iCell)
           zTop(k:nVertLevels,iCell) = -bottomDepth(iCell) +     layerThickness(k,iCell)

           do k=maxLevelCell(iCell)-1, 1, -1
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
      coef = -gravity/config_density0
      do iCell=1,nCells
         BruntVaisalaFreqTop(1,iCell) = 0.0
         do k=2,maxLevelCell(iCell)
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
      call computeKPPInputFields(state, forcing, mesh, diagnostics, scratch)

      !
      ! Apply the GM closure as a bolus velocity
      !
      if (config_h_kappa .GE. epsilon(0D0)) then
         call ocn_gm_compute_uBolus(state, diagnostics, mesh)
      else
         uBolusGM = 0.0
      end if

#ifdef MPAS_CESM
      do iEdge = 1, mesh % nEdgesSolve
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)

         diagnostics % gradSSH % array(1, iEdge) = (ssh(cell2) - ssh(cell1)) / dcEdge(iEdge)
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
   subroutine ocn_vert_transport_velocity_top(mesh, verticalMesh, oldLayerThickness, layerThicknessEdge, &
     normalVelocity, oldSSH, newHighFreqThickness, dt, vertTransportVelocityTop, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh           !< Input: horizonal mesh information

      type (verticalMesh_type), intent(in) :: &
         verticalMesh   !< Input: vertical mesh information

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
         vertTransportVelocityTop     !< Output: vertical transport at top of cell

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, iCell, k, i, nCells, nVertLevels
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

      err = 0

      nEdgesOnCell      => mesh % nEdgesOnCell % array
      areaCell          => mesh % areaCell % array
      edgesOnCell       => mesh % edgesOnCell % array
      edgeSignOnCell    => mesh % edgeSignOnCell % array
      maxLevelCell      => mesh % maxLevelCell % array
      maxLevelEdgeBot   => mesh % maxLevelEdgeBot % array
      dvEdge            => mesh % dvEdge % array

      nCells      = mesh % nCells
      nVertLevels = mesh % nVertLevels

      if (config_vert_coord_movement.eq.'impermeable_interfaces') then
        vertTransportVelocityTop=0.0
        return
      end if

      allocate(div_hu(nVertLevels,nCells), div_hu_btr(nCells), ALE_Thickness(nVertLevels,nCells))

      !
      ! thickness-weighted divergence and barotropic divergence
      !
      ! See Ringler et al. (2010) jcp paper, eqn 19, 21, and fig. 3.
      do iCell=1,nCells
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
      call ocn_ALE_thickness(mesh, verticalMesh, oldSSH, div_hu_btr, newHighFreqThickness, dt, ALE_thickness, err)

      !
      ! Vertical transport through layer interfaces
      !
      ! Vertical transport through layer interface at top and bottom is zero.
      ! Here we are using solving the continuity equation for vertTransportVelocityTop ($w^t$),
      ! and using ALE_Thickness for thickness at the new time.

      do iCell=1,nCells
         vertTransportVelocityTop(1,iCell) = 0.0
         vertTransportVelocityTop(maxLevelCell(iCell)+1,iCell) = 0.0
         do k=maxLevelCell(iCell),2,-1
            vertTransportVelocityTop(k,iCell) = vertTransportVelocityTop(k+1,iCell) - div_hu(k,iCell) &
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

   subroutine ocn_fuperp(state, mesh)!{{{

      type (state_type), intent(inout) :: state !< Input/Output: State information
      type (mesh_type), intent(in) :: mesh !< Input: mesh information

      integer :: iEdge, cell1, cell2, eoe, i, j, k
      integer :: nEdgesSolve
      real (kind=RKIND), dimension(:), pointer :: fEdge
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge, normalVelocity, normalBaroclinicVelocity
      type (dm_info) :: dminfo

      integer, dimension(:), pointer :: maxLevelEdgeTop, nEdgesOnEdge
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnEdge

      call mpas_timer_start("ocn_fuperp")

      normalVelocity           => state % normalVelocity % array
      normalBaroclinicVelocity        => state % normalBaroclinicVelocity % array
      weightsOnEdge     => mesh % weightsOnEdge % array
      fEdge             => mesh % fEdge % array
      maxLevelEdgeTop      => mesh % maxLevelEdgeTop % array
      cellsOnEdge       => mesh % cellsOnEdge % array
      nEdgesOnEdge      => mesh % nEdgesOnEdge % array
      edgesOnEdge       => mesh % edgesOnEdge % array

      fEdge       => mesh % fEdge % array

      nEdgesSolve = mesh % nEdgesSolve

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
   subroutine ocn_filter_btr_mode_vel(state, diagnostics, mesh)!{{{

      type (state_type), intent(inout) :: state !< Input/Output: State information
      type (diagnostics_type), intent(in) :: diagnostics !< Input: Diagnostics information
      type (mesh_type), intent(in) :: mesh !< Input: Mesh information

      integer :: iEdge, k, nEdges
      real (kind=RKIND) :: vertSum, normalThicknessFluxSum, thicknessSum
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge, normalVelocity
      integer, dimension(:), pointer :: maxLevelEdgeTop

      call mpas_timer_start("ocn_filter_btr_mode_vel")

      normalVelocity => state % normalVelocity % array
      layerThicknessEdge => diagnostics % layerThicknessEdge % array
      maxLevelEdgeTop => mesh % maxLevelEdgeTop % array
      nEdges      = mesh % nEdges

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
   subroutine ocn_filter_btr_mode_tend_vel(tend, state, diagnostics, mesh)!{{{

      type (tend_type), intent(inout) :: tend !< Input/Output: Tendency information
      type (state_type), intent(in) :: state !< Input: State information
      type (diagnostics_type), intent(in) :: diagnostics !< Input: Diagnostics information
      type (mesh_type), intent(in) :: mesh !< Input: Mesh information

      integer :: iEdge, k, nEdges
      real (kind=RKIND) :: vertSum, normalThicknessFluxSum, thicknessSum
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge, tend_normalVelocity

      integer, dimension(:), pointer :: maxLevelEdgeTop

      call mpas_timer_start("ocn_filter_btr_mode_tend_vel")

      tend_normalVelocity => tend % normalVelocity % array
      layerThicknessEdge => diagnostics % layerThicknessEdge % array
      maxLevelEdgeTop => mesh % maxLevelEdgeTop % array
      nEdges      = mesh % nEdges

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

    subroutine computeKPPInputFields(state, forcing, mesh, diagnostics, scratch)!{{{

      type (state_type), intent(inout) :: state !< Input/Output: State information
      type (forcing_type), intent(in) :: forcing !< Input: Forcing information
      type (mesh_type), intent(in) :: mesh !< Input: Mesh information
      type (diagnostics_type), intent(inout) :: diagnostics !< Diagnostics information derived from State
      type (scratch_type), intent(inout) :: scratch !< Input: scratch variables

      ! scalars
      integer :: nCells, nVertLevels

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
      integer :: iCell, iEdge, i, k, err, indexTempFlux, indexSaltFlux
      real (kind=RKIND) :: numerator, denominator, factor, deltaVelocitySquared, turbulentVelocitySquared, delU2, invAreaCell

      ! set the parameter turbulentVelocitySquared
      turbulentVelocitySquared = 0.001

      ! set scalar values
      nCells      = mesh % nCells
      nVertLevels = mesh % nVertLevels
      indexTempFlux = forcing % index_surfaceTemperatureFlux
      indexSaltFlux = forcing % index_surfaceSalinityFlux

      ! set pointers into state, mesh, diagnostics and scratch
      normalVelocity => state % normalVelocity % array
      layerThickness => state % layerThickness % array

      maxLevelCell => mesh % maxLevelCell % array 
      nEdgesOnCell => mesh % nEdgesOnCell % array
      edgesOnCell  => mesh % edgesOnCell % array
      areaCell     => mesh % areaCell % array
      dcEdge       => mesh % dcEdge % array
      dvEdge       => mesh % dvEdge % array

      zMid                          => diagnostics % zMid % array
      zTop                          => diagnostics % zTop % array
      density                       => diagnostics % density % array
      bulkRichardsonNumber          => diagnostics % bulkRichardsonNumber % array
      tracersSurfaceValues          => diagnostics % tracersSurfaceValue  % array
      boundaryLayerDepth            => diagnostics % boundaryLayerDepth % array
      surfaceFrictionVelocity       => diagnostics % surfaceFrictionVelocity % array
      penetrativeTemperatureFluxOBL => diagnostics % penetrativeTemperatureFluxOBL % array
      buoyancyForcingOBL            => diagnostics % buoyancyForcingOBL % array

      normalVelocity => state % normalVelocity % array

      surfaceMassFlux            => forcing % surfaceMassFlux % array
      surfaceTracerFlux          => forcing % surfaceTracerFlux % array
      penetrativeTemperatureFlux => forcing % penetrativeTemperatureFlux % array 
      surfaceWindStressMagnitude => forcing % surfaceWindStressMagnitude % array      

      ! allocate scratch space displaced density computation
      call mpas_allocate_scratch_field(scratch % densitySurfaceDisplaced, .true.)
      call mpas_allocate_scratch_field(scratch % thermalExpansionCoeff, .true.)
      call mpas_allocate_scratch_field(scratch % salineContractionCoeff, .true.)
      densitySurfaceDisplaced => scratch % densitySurfaceDisplaced % array
      thermalExpansionCoeff => scratch % thermalExpansionCoeff % array
      salineContractionCoeff => scratch % salineContractionCoeff % array

      ! compute EOS by displacing SST/SSS to every vertical layer in column
      call ocn_equation_of_state_density(state, diagnostics, mesh, 0, 'surfaceDisplaced', densitySurfaceDisplaced, err, &
              thermalExpansionCoeff, salineContractionCoeff)

      ! set value to out-of-bounds
      bulkRichardsonNumber(:,:) = -1.0e34

      do iCell=1,nCells
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
       do k=1,maxLevelCell(iCell)

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
      call mpas_deallocate_scratch_field(scratch % densitySurfaceDisplaced, .true.)
      call mpas_deallocate_scratch_field(scratch % thermalExpansionCoeff, .true.)
      call mpas_deallocate_scratch_field(scratch % salineContractionCoeff, .true.)

    end subroutine computeKPPInputFields!}}}


end module ocn_diagnostics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
