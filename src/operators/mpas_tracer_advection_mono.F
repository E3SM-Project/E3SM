! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_tracer_advection_mono
!
!> \brief MPAS monotonic tracer advection with FCT
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This module contains routines for monotonic advection of tracers using a FCT
!
!-----------------------------------------------------------------------
module mpas_tracer_advection_mono

   use mpas_kind_types
   use mpas_grid_types
   use mpas_io_units

   use mpas_tracer_advection_helpers

   implicit none
   private
   save 

   real (kind=RKIND) :: coef_3rd_order
   integer :: horizOrder
   logical :: vert2ndOrder, vert3rdOrder, vert4thOrder
   logical :: positiveDzDk, monotonicityCheck

   public :: mpas_tracer_advection_mono_tend, &
             mpas_tracer_advection_mono_init

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_tracer_advection_mono_tend
!
!> \brief MPAS monotonic tracer advection tendency with FCT
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine computes the monotonic tracer advection tendencity using a FCT.
!>  Both horizontal and vertical.
!
!-----------------------------------------------------------------------
   subroutine mpas_tracer_advection_mono_tend(tracers, adv_coefs, adv_coefs_3rd, nAdvCellsForEdge, advCellsForEdge, normalThicknessFlux, &!{{{
                                              w, layerThickness, verticalCellSize, dt, mesh, tend_layerThickness, tend, maxLevelCell_in, & 
                                              maxLevelEdgeTop_in, highOrderAdvectionMask_in, verticalDivergenceFactor_in, edgeSignOnCell_in)

      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers !< Input: current tracer values
      real (kind=RKIND), dimension(:,:), intent(in) :: adv_coefs !< Input: Advection coefficients for 2nd order advection
      real (kind=RKIND), dimension(:,:), intent(in) :: adv_coefs_3rd !< Input: Advection coefficients for blending in 3rd or 4th order advection
      integer, dimension(:), intent(in) :: nAdvCellsForEdge !< Input: Number of advection cells for each edge
      integer, dimension(:,:), intent(in) :: advCellsForEdge !< Input: List of advection cells for each edge
      real (kind=RKIND), dimension(:,:), intent(in) :: normalThicknessFlux !< Input: Thichness weighted velocitiy
      real (kind=RKIND), dimension(:,:), intent(in) :: w !< Input: Vertical velocitiy
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness !< Input: Thickness
      real (kind=RKIND), dimension(:,:), intent(in) :: verticalCellSize !< Input: Distance between vertical interfaces of a cell
      real (kind=RKIND), dimension(:,:), intent(in) :: tend_layerThickness !< Input: Tendency for thickness field
      real (kind=RKIND), intent(in) :: dt !< Input: Timestep
      type (mesh_type), intent(in) :: mesh !< Input: Mesh information
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tend !< Input/Output: Tracer tendency
      integer, dimension(:), pointer, optional :: maxLevelCell_in !< Input - optional: Index to max level at cell center
      integer, dimension(:), pointer, optional :: maxLevelEdgeTop_in !< Input - optional: Index to max level at edge with non-land cells on both sides
      integer, dimension(:,:), pointer, optional :: highOrderAdvectionMask_in !< Input - optional: Mask for high order advection
      real (kind=RKIND), dimension(:), pointer, optional :: verticalDivergenceFactor_in !< Input - optional: Denominator for vertical divergence. Given as an inverse which can be multiplied.
      integer, dimension(:, :), pointer, optional :: edgeSignOnCell_in !< Input - optional: Sign for flux from edge on each cell. Used for bit-reproducibility

      integer :: i, iCell, iEdge, k, iTracer, cell1, cell2
      integer :: nVertLevels, num_tracers, nCells, nEdges, nCellsSolve, maxEdges
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnCell, highOrderAdvectionMask, edgeSignOnCell, edgesOnCell

      real (kind=RKIND) :: flux_upwind, tracer_min_new, tracer_max_new, tracer_upwind_new, scale_factor
      real (kind=RKIND) :: flux, tracer_weight, invAreaCell1, invAreaCell2
      real (kind=RKIND) :: verticalWeightK, verticalWeightKm1
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell, verticalDivergenceFactor
      real (kind=RKIND), dimension(:,:), pointer :: tracer_cur, tracer_new, upwind_tendency, inv_h_new, tracer_max, tracer_min
      real (kind=RKIND), dimension(:,:), pointer :: flux_incoming, flux_outgoing, high_order_horiz_flux, high_order_vert_flux

      real (kind=RKIND), parameter :: eps = 1.e-10_RKIND

      ! Get dimensions
      nCells = mesh % nCells
      nCellsSolve = mesh % nCellsSolve
      nEdges = mesh % nEdges
      maxEdges = mesh % maxEdges
      nVertLevels = mesh % nVertLevels
      num_tracers = size(tracers,dim=1)

      ! Initialize pointers
      dvEdge      => mesh % dvEdge % array
      cellsOnEdge => mesh % cellsOnEdge % array
      edgesOnCell => mesh % edgesOnCell % array
      cellsOnCell => mesh % cellsOnCell % array
      areaCell    => mesh % areaCell % array
      nEdgesOnCell => mesh % nEdgesOnCell % array
      
      ! Setup arrays from optional arguments
      if(present(maxLevelCell_in)) then
        maxLevelCell => maxLevelCell_in
      else
        allocate(maxLevelCell(nCells+1))
        maxLevelCell(:) = nVertLevels
      end if

      if(present(maxLevelEdgeTop_in)) then
        maxLevelEdgeTop => maxLevelEdgeTop_in
      else
        allocate(maxLevelEdgeTop(nEdges+1))
        maxLevelEdgeTop(:) = nVertLevels
      end if

      if(present(highOrderAdvectionMask_in)) then
        highOrderAdvectionMask => highOrderAdvectionMask_in
      else
        allocate(highOrderAdvectionMask(nVertLevels, nEdges+1))
        if(horizOrder == 2) then
          highOrderAdvectionMask = 0
        else
          highOrderAdvectionMask = 1
        end if
      end if

      if(present(verticalDivergenceFactor_in)) then
        verticalDivergenceFactor => verticalDivergenceFactor_in
      else
        allocate(verticalDivergenceFactor(nVertLevels))
        verticalDivergenceFactor = 1.0_RKIND
      end if

      if(present(edgeSignOnCell_in)) then
        edgeSignOnCell => edgeSignOnCell_in
      else
        allocate(edgeSignOnCell(maxEdges, nCells+1))
        do iCell = 1, nCells
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            if(iCell == cellsOnEdge(1, iEdge)) then
              edgeSignOnCell(i, iCell) = -1
            else
              edgeSignOnCell(i, iCell) =  1
            end if
          end do
        end do
      end if

      ! Setup high order horizontal flux field
      allocate(high_order_horiz_flux(nVertLevels, nEdges+1))

      ! allocate nCells arrays
      allocate(tracer_new(nVertLevels, nCells+1))
      allocate(tracer_cur(nVertLevels, nCells+1))
      allocate(upwind_tendency(nVertLevels, nCells+1))
      allocate(inv_h_new(nVertLevels, nCells+1))
      allocate(tracer_max(nVertLevels, nCells+1))
      allocate(tracer_min(nVertLevels, nCells+1))
      allocate(flux_incoming(nVertLevels, nCells+1))
      allocate(flux_outgoing(nVertLevels, nCells+1))

      ! allocate nVertLevels+1 and nCells arrays
      allocate(high_order_vert_flux(nVertLevels+1, nCells))

      do iCell = 1, nCells
        do k=1, maxLevelCell(iCell)
          inv_h_new(k, iCell) = 1.0 / (layerThickness(k, iCell) + dt * tend_layerThickness(k, iCell))
        end do
      end do

      ! Loop over tracers. One tracer is advected at a time. It is copied into a temporary array in order to improve locality
      do iTracer = 1, num_tracers
        ! Initialize variables for use in this iTracer iteration
        do iCell = 1, nCells
          do k=1, maxLevelCell(iCell)
            tracer_cur(k,iCell) = tracers(iTracer,k,iCell)
            upwind_tendency(k, iCell) = 0.0_RKIND

            !tracer_new is supposed to be the "new" tracer state. This allows bounds checks.
            if (monotonicityCheck) then
              tracer_new(k,iCell) = 0.0_RKIND
            end if
          end do ! k loop
        end do ! iCell loop

        high_order_vert_flux = 0.0_RKIND
        high_order_horiz_flux = 0.0_RKIND

        !  Compute the high order vertical flux. Also determine bounds on tracer_cur. 
        do iCell = 1, nCells
          k = 1
          tracer_max(k,iCell) = max(tracer_cur(k,iCell),tracer_cur(k+1,iCell))
          tracer_min(k,iCell) = min(tracer_cur(k,iCell),tracer_cur(k+1,iCell))

          k = max(1, min(maxLevelCell(iCell), 2))
          verticalWeightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
          verticalWeightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
          high_order_vert_flux(k,iCell) = w(k,iCell)*(verticalWeightK*tracer_cur(k,iCell)+verticalWeightKm1*tracer_cur(k-1,iCell))
          tracer_max(k,iCell) = max(tracer_cur(k-1,iCell),tracer_cur(k,iCell),tracer_cur(k+1,iCell))
          tracer_min(k,iCell) = min(tracer_cur(k-1,iCell),tracer_cur(k,iCell),tracer_cur(k+1,iCell))
             
          do k=3,maxLevelCell(iCell)-1
             if(vert4thOrder) then
               high_order_vert_flux(k, iCell) = mpas_tracer_advection_vflux4( tracer_cur(k-2,iCell),tracer_cur(k-1,iCell),  &
                                      tracer_cur(k  ,iCell),tracer_cur(k+1,iCell), w(k,iCell))
             else if(vert3rdOrder) then
               high_order_vert_flux(k, iCell) = mpas_tracer_advection_vflux3( tracer_cur(k-2,iCell),tracer_cur(k-1,iCell),  &
                                      tracer_cur(k  ,iCell),tracer_cur(k+1,iCell), w(k,iCell), coef_3rd_order )
             else if (vert2ndOrder) then
               verticalWeightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
               verticalWeightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
               high_order_vert_flux(k,iCell) = w(k,iCell)*(verticalWeightK*tracer_cur(k,iCell)+verticalWeightKm1*tracer_cur(k-1,iCell))
             end if
             tracer_max(k,iCell) = max(tracer_cur(k-1,iCell),tracer_cur(k,iCell),tracer_cur(k+1,iCell))
             tracer_min(k,iCell) = min(tracer_cur(k-1,iCell),tracer_cur(k,iCell),tracer_cur(k+1,iCell))
          end do 
 
          k = max(1, maxLevelCell(iCell))
          verticalWeightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
          verticalWeightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
          high_order_vert_flux(k,iCell) = w(k,iCell)*(verticalWeightK*tracer_cur(k,iCell)+verticalWeightKm1*tracer_cur(k-1,iCell))
          tracer_max(k,iCell) = max(tracer_cur(k,iCell),tracer_cur(k-1,iCell))
          tracer_min(k,iCell) = min(tracer_cur(k,iCell),tracer_cur(k-1,iCell))

          ! pull tracer_min and tracer_max from the (horizontal) surrounding cells
          do i = 1, nEdgesOnCell(iCell)
            do k=1, min(maxLevelCell(iCell), maxLevelCell(cellsOnCell(i, iCell)))
              tracer_max(k,iCell) = max(tracer_max(k,iCell),tracer_cur(k, cellsOnCell(i,iCell)))
              tracer_min(k,iCell) = min(tracer_min(k,iCell),tracer_cur(k, cellsOnCell(i,iCell)))
            end do ! k loop
          end do ! i loop over nEdgesOnCell
        end do ! iCell Loop

        !  Compute the high order horizontal flux
        do iEdge = 1, nEdges
          cell1 = cellsOnEdge(1, iEdge)
          cell2 = cellsOnEdge(2, iEdge)

          ! Compute 2nd order fluxes where needed.
          do k = 1, maxLevelEdgeTop(iEdge)
            tracer_weight = iand(highOrderAdvectionMask(k, iEdge)+1, 1) * (dvEdge(iEdge) * 0.5_RKIND) * normalThicknessFlux(k, iEdge)

            high_order_horiz_flux(k, iEdge) = high_order_horiz_flux(k, iedge) + tracer_weight * (tracer_cur(k, cell1) + tracer_cur(k, cell2))
          end do ! k loop

          ! Compute 3rd or 4th fluxes where requested.
          do i = 1, nAdvCellsForEdge(iEdge)
            iCell = advCellsForEdge(i,iEdge)
            do k = 1, maxLevelCell(iCell)
              tracer_weight = highOrderAdvectionMask(k, iEdge) * (adv_coefs(i,iEdge) + coef_3rd_order*sign(1.0_RKIND,normalThicknessFlux(k,iEdge))*adv_coefs_3rd(i,iEdge))

              tracer_weight = normalThicknessFlux(k,iEdge)*tracer_weight
              high_order_horiz_flux(k,iEdge) = high_order_horiz_flux(k,iEdge) + tracer_weight * tracer_cur(k,iCell)
            end do ! k loop
          end do ! i loop over nAdvCellsForEdge
        end do ! iEdge loop

        !  low order upwind vertical flux (monotonic and diffused)
        !  Remove low order flux from the high order flux.
        !  Store left over high order flux in high_order_vert_flux array.
        !  Upwind fluxes are accumulated in upwind_tendency
        do iCell = 1, nCells
          do k = 2, maxLevelCell(iCell)
            ! dwj 02/03/12 and Atmosphere are different in vertical
            if(positiveDzDk) then
              flux_upwind = max(0.0_RKIND,w(k,iCell))*tracer_cur(k-1,iCell) + min(0.0_RKIND,w(k,iCell))*tracer_cur(k,iCell)
            else
              flux_upwind = min(0.0_RKIND,w(k,iCell))*tracer_cur(k-1,iCell) + max(0.0_RKIND,w(k,iCell))*tracer_cur(k,iCell)
            end if
            upwind_tendency(k-1,iCell) = upwind_tendency(k-1,iCell) + flux_upwind
            upwind_tendency(k  ,iCell) = upwind_tendency(k  ,iCell) - flux_upwind
            high_order_vert_flux(k,iCell) = high_order_vert_flux(k,iCell) - flux_upwind
          end do ! k loop

          ! flux_incoming contains the total remaining high order flux into iCell
          !          it is positive.
          ! flux_outgoing contains the total remaining high order flux out of iCell
          !           it is negative
          do k = 1, maxLevelCell(iCell)
            ! dwj 02/03/12 and Atmosphere are different in vertical
            if(positiveDzDk) then
              flux_incoming (k,iCell) = -(min(0.0_RKIND,high_order_vert_flux(k+1,iCell))-max(0.0_RKIND,high_order_vert_flux(k,iCell)))
              flux_outgoing(k,iCell) = -(max(0.0_RKIND,high_order_vert_flux(k+1,iCell))-min(0.0_RKIND,high_order_vert_flux(k,iCell)))
            else
              flux_incoming (k, iCell) = max(0.0_RKIND, high_order_vert_flux(k+1, iCell)) - min(0.0_RKIND, high_order_vert_flux(k, iCell))
              flux_outgoing(k, iCell) = min(0.0_RKIND, high_order_vert_flux(k+1, iCell)) - max(0.0_RKIND, high_order_vert_flux(k, iCell))
            end if
          end do ! k Loop
        end do ! iCell Loop

        !  low order upwind horizontal flux (monotinc and diffused)
        !  Remove low order flux from the high order flux
        !  Store left over high order flux in high_order_horiz_flux array
        !  Upwind fluxes are accumulated in upwind_tendency
        do iEdge = 1, nEdges
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)

          invAreaCell1 = 1.0_RKIND / areaCell(cell1)
          invAreaCell2 = 1.0_RKIND / areaCell(cell2)

          do k = 1, maxLevelEdgeTop(iEdge)
            flux_upwind = dvEdge(iEdge) * (max(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracer_cur(k,cell1) + min(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracer_cur(k,cell2))
            high_order_horiz_flux(k,iEdge) = high_order_horiz_flux(k,iEdge) - flux_upwind
          end do ! k loop
        end do ! iEdge loop

        do iCell = 1, nCells
          invAreaCell1 = 1.0_RKIND / areaCell(iCell)
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k = 1, maxLevelEdgeTop(iEdge)
              flux_upwind = dvEdge(iEdge) * (max(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracer_cur(k,cell1) + min(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracer_cur(k,cell2))

              upwind_tendency(k,iCell) = upwind_tendency(k,iCell) + edgeSignOncell(i, iCell) * flux_upwind * invAreaCell1

              ! Accumulate remaining high order fluxes
              flux_outgoing(k,iCell) = flux_outgoing(k,iCell) + min(0.0_RKIND, edgeSignOnCell(i, iCell) * high_order_horiz_flux(k, iEdge)) * invAreaCell1
              flux_incoming(k,iCell) = flux_incoming(k,iCell) + max(0.0_RKIND, edgeSignOnCell(i, iCell) * high_order_horiz_flux(k, iEdge)) * invAreaCell1
            end do
          end do
        end do

        ! Build the factors for the FCT
        ! Computed using the bounds that were computed previously, and the bounds on the newly updated value
        ! Factors are placed in the flux_incoming and flux_outgoing arrays
        do iCell = 1, nCells
          do k = 1, maxLevelCell(iCell)
            tracer_min_new = (tracer_cur(k,iCell)*layerThickness(k,iCell) + dt*(upwind_tendency(k,iCell)+flux_outgoing(k,iCell))) * inv_h_new(k,iCell)
            tracer_max_new = (tracer_cur(k,iCell)*layerThickness(k,iCell) + dt*(upwind_tendency(k,iCell)+flux_incoming(k,iCell))) * inv_h_new(k,iCell)
            tracer_upwind_new = (tracer_cur(k,iCell)*layerThickness(k,iCell) + dt*upwind_tendency(k,iCell)) * inv_h_new(k,iCell)
           
            scale_factor = (tracer_max(k,iCell)-tracer_upwind_new)/(tracer_max_new-tracer_upwind_new+eps)
            flux_incoming(k,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )

            scale_factor = (tracer_upwind_new-tracer_min(k,iCell))/(tracer_upwind_new-tracer_min_new+eps)
            flux_outgoing(k,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )
          end do ! k loop
        end do ! iCell loop

        !  rescale the high order horizontal fluxes
        do iEdge = 1, nEdges
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)
          do k = 1, maxLevelEdgeTop(iEdge)
            flux = high_order_horiz_flux(k,iEdge)
            flux = max(0.0_RKIND,flux) * min(flux_outgoing(k,cell1), flux_incoming(k,cell2)) &
                 + min(0.0_RKIND,flux) * min(flux_incoming(k,cell1), flux_outgoing(k,cell2))
            high_order_horiz_flux(k,iEdge) = flux
          end do ! k loop
        end do ! iEdge loop

        ! rescale the high order vertical flux
        do iCell = 1, nCellsSolve
          do k = 2, maxLevelCell(iCell)
            flux =  high_order_vert_flux(k,iCell)
            ! dwj 02/03/12 and Atmosphere are different in vertical.
            if(positiveDzDk) then
              flux = max(0.0_RKIND,flux) * min(flux_outgoing(k-1,iCell), flux_incoming(k  ,iCell)) &
                   + min(0.0_RKIND,flux) * min(flux_outgoing(k  ,iCell), flux_incoming(k-1,iCell))
            else
              flux = max(0.0_RKIND,flux) * min(flux_outgoing(k  ,iCell), flux_incoming(k-1,iCell)) &
                   + min(0.0_RKIND,flux) * min(flux_outgoing(k-1,iCell), flux_incoming(k  ,iCell))
            end if
            high_order_vert_flux(k,iCell) = flux
          end do ! k loop
        end do ! iCell loop

        ! Accumulate the scaled high order horizontal tendencies
        do iCell = 1, nCells
          invAreaCell1 = 1.0 / areaCell(iCell)
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelEdgeTop(iEdge)
              tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + edgeSignOnCell(i, iCell) * high_order_horiz_flux(k, iEdge) * invAreaCell1

              if(monotonicityCheck) then
                tracer_new(k, iCell) = tracer_new(k, iCell) + edgeSignOnCell(i, iCell) * high_order_horiz_flux(k, iEdge) * invAreaCell1
              end if
            end do
          end do
        end do

        ! Accumulate the scaled high order vertical tendencies, and the upwind tendencies
        do iCell = 1, nCellsSolve
          do k = 1,maxLevelCell(iCell)
            tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + verticalDivergenceFactor(k) * (high_order_vert_flux(k+1, iCell) - high_order_vert_flux(k, iCell)) + upwind_tendency(k,iCell)

            if (monotonicityCheck) then
              !tracer_new holds a tendency for now. Only for a check on monotonicity
              tracer_new(k, iCell) = tracer_new(k, iCell) + verticalDivergenceFactor(k) * (high_order_vert_flux(k+1, iCell) - high_order_vert_flux(k, iCell)) + upwind_tendency(k,iCell)

              !tracer_new is now the new state of the tracer. Only for a check on monotonicity
              tracer_new(k, iCell) = (tracer_cur(k, iCell)*layerThickness(k, iCell) + dt * tracer_new(k, iCell)) * inv_h_new(k, iCell)
            end if
          end do ! k loop
        end do ! iCell loop

        if (monotonicityCheck) then
          !build min and max bounds on old and new tracer for check on monotonicity.
          do iCell = 1, nCellsSolve
            do k = 1, maxLevelCell(iCell)
              if(tracer_new(k,iCell) < tracer_min(k, iCell)-eps) then
                write(stderrUnit,*) 'Minimum out of bounds on tracer ', iTracer, tracer_min(k, iCell), tracer_new(k,iCell)
              end if

              if(tracer_new(k,iCell) > tracer_max(k,iCell)+eps) then
                write(stderrUnit,*) 'Maximum out of bounds on tracer ', iTracer, tracer_max(k, iCell), tracer_new(k,iCell)
              end if
            end do
          end do
        end if
      end do ! iTracer loop

      deallocate(tracer_new)
      deallocate(tracer_cur)
      deallocate(upwind_tendency)
      deallocate(inv_h_new)
      deallocate(tracer_max)
      deallocate(tracer_min)
      deallocate(flux_incoming)
      deallocate(flux_outgoing)
      deallocate(high_order_horiz_flux)
      deallocate(high_order_vert_flux)

      if(.not.present(maxLevelCell_in)) then
        deallocate(maxLevelCell)
      end if

      if(.not.present(maxLevelEdgeTop_in)) then
        deallocate(maxLevelEdgeTop)
      end if

      if(.not.present(highOrderAdvectionMask_in)) then
        deallocate(highOrderAdvectionMask)
      end if

      if(.not.present(verticalDivergenceFactor_in)) then
        deallocate(verticalDivergenceFactor)
      end if
   end subroutine mpas_tracer_advection_mono_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_tracer_advection_mono_init
!
!> \brief MPAS initialize monotonic tracer advection tendency with FCT
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine initializes the monotonic tracer advection tendencity using a FCT.
!
!-----------------------------------------------------------------------
   subroutine mpas_tracer_advection_mono_init(horiz_adv_order, vert_adv_order, coef_3rd_order_in, dzdk_positive, check_monotonicity, err)!{{{

      use mpas_dmpar
      integer, intent(in) :: horiz_adv_order !< Input: Order for horizontal advection
      integer, intent(in) :: vert_adv_order !< Input: Order for vertical advection
      real (kind=RKIND), intent(in) :: coef_3rd_order_in !< Input: coefficient for blending advection orders.
      logical, intent(in) :: dzdk_positive !< Input: Logical flag determining if dzdk is positive or negative.
      logical, intent(in) :: check_monotonicity !< Input: Logical flag determining check on monotonicity of tracers
      integer, intent(inout) :: err !< Input/Output: Error Flag

      err = 0

      vert2ndOrder = .false.
      vert3rdOrder = .false.
      vert4thOrder = .false.

      if ( horiz_adv_order == 3) then
          coef_3rd_order = coef_3rd_order_in
      else if(horiz_adv_order == 2 .or. horiz_adv_order == 4) then
          coef_3rd_order = 0.0_RKIND
      end if

      horizOrder = horiz_adv_order

      if (vert_adv_order == 3) then
          vert3rdOrder = .true.
      else if (vert_adv_order == 4) then
          vert4thOrder = .true.
      else
          vert2ndOrder = .true.
          if(vert_adv_order /= 2) then
            write(stderrUnit,*) 'Invalid value for vert_adv_order, defaulting to 2nd order'
          end if
      end if

      if (config_num_halos < 3) then
        call mpas_dmpar_global_abort('ERROR: Monotonic advection cannot be used with less than 3 halos.')
      end if

      positiveDzDk = dzdk_positive
      monotonicityCheck = check_monotonicity

   end subroutine mpas_tracer_advection_mono_init!}}}

end module mpas_tracer_advection_mono
