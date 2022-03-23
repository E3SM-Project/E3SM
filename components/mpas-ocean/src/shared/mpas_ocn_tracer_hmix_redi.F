! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_hmix_Redi
!
!> \brief MPAS ocean horizontal tracer mixing driver
!> \details
!>  This module contains the main driver routine for computing
!>  horizontal mixing tendencies.
!>
!>  It provides an init and a tend function. Each are described below.
!
!-----------------------------------------------------------------------

module ocn_tracer_hmix_Redi

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants

   use ocn_config
   use ocn_constants

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

   public :: ocn_tracer_hmix_Redi_tend, &
             ocn_tracer_hmix_Redi_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real(kind=RKIND) :: term1TaperFactor
   real(kind=RKIND), dimension(:), allocatable :: rediKappaCell
   real(kind=RKIND), dimension(:, :, :), allocatable :: redi_term2_edge, redi_term3_topOfCell
   real(kind=RKIND), dimension(:, :, :), allocatable :: redi_term1, redi_term2, redi_term3

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_hmix_Redi_tend
!
!> \brief   Computes Laplacian tendency term for horizontal tracer mixing
!> \details
!>  This routine computes the horizontal mixing tendency for tracers
!>  based on current state using a Laplacian parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hmix_Redi_tend(meshPool, layerThickness, layerThicknessEdge, zMid, tracers, &
                                        RediKappa, dt, isActiveTracer, slopeTriadUp, slopeTriadDown, &
                                        RediKappaSfcTaper, RediKappaScaling, rediLimiterCount, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type(mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information

      real(kind=RKIND), dimension(:), intent(in) :: &
         RediKappa

      real(kind=RKIND), dimension(:, :), intent(in) :: &
         layerThicknessEdge, &!< Input: thickness at edge
         zMid, &!< Input: Z coordinate at the center of a cell
         RediKappaSfcTaper, &!< Input: vertical structure of GM/Redi limiter
         RediKappaScaling, &
         layerThickness

      real(kind=RKIND), dimension(:, :, :), intent(in) :: &
         slopeTriadUp, slopeTriadDown, &
         tracers !< Input: tracer quantities

      real(kind=RKIND), intent(in) :: dt

      integer, dimension(:, :), intent(inout) :: rediLimiterCount

      logical, intent(in) :: isActiveTracer

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(:, :, :), intent(inout) :: &
         tend

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

      integer :: iCell, iEdge, cell1, cell2, iCellSelf
      integer :: i, k, iTr, nTracers, nCells, nEdges, nCellsP1
      logical :: use_Redi_diag_terms
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      integer, dimension(:), pointer :: minLevelEdgeBot, maxLevelEdgeTop, nEdgesOnCell, minLevelCell, maxLevelCell
      integer, dimension(:, :), pointer :: cellsOnEdge, edgesOnCell, edgeSignOnCell, cellsOnCell

      real(kind=RKIND) :: tempTracer, invAreaCell
      real(kind=RKIND) :: flux, flux_term2, flux_term3, dTracerDx, coef
      real(kind=RKIND) :: r_tmp, tracer_turb_flux, kappaRediEdgeVal
      real(kind=RKIND), dimension(:), pointer :: areaCell, dvEdge, dcEdge
      real(kind=RKIND), dimension(:), allocatable :: minimumVal
      real(kind=RKIND), dimension(:, :), allocatable :: fluxRediZTop
      err = 0

      call mpas_timer_start("tracer redi")

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      nTracers = size(tracers, dim=1)

      call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

      nCells = nCellsArray(size(nCellsArray))
      nEdges = nEdgesArray(size(nEdgesArray))

      allocate (minimumVal(nTracers))
      allocate (fluxRediZTop(nTracers, nVertLevels + 1))
      allocate (rediKappaCell(nCells))
      allocate (redi_term1(nTracers, nVertLevels, nCells))
      allocate (redi_term2(nTracers, nVertLevels, nCells))
      allocate (redi_term3(nTracers, nVertLevels, nCells))
      allocate (redi_term2_edge(nTracers, nVertLevels, nEdges))
      allocate (redi_term3_topOfCell(nTracers, nVertLevels + 1, nCells))

      nCells = nCellsArray(2)
      !$omp parallel
      !$omp do schedule(runtime) private(i, iEdge)
      do iCell = 1,nCells
        rediKappaCell(iCell) = 0.0_RKIND
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          rediKappaCell(iCell) = rediKappaCell(iCell) + 0.25_RKIND*RediKappa(iEdge) * dvEdge(iEdge) * dcEdge(iEdge)
        enddo
        rediKappaCell(iCell) = rediKappaCell(iCell) / areaCell(iCell)
      enddo
      !$omp end do
      !$omp end parallel

      nCells = nCellsArray(2)
      nCellsP1 = nCellsArray(size(nCellsArray)) + 1

      ! Term 1: this is the "standard" horizontal del2 term, but with RediKappa coefficient.
      ! \kappa_2 \nabla \phi on edge
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(fluxRediZTop, invAreaCell, i, iEdge, cell1, cell2, iCellSelf, r_tmp, coef, k, &
      !$omp         kappaRediEdgeVal, iTr, tracer_turb_flux, flux, flux_term2, flux_term3, &
      !$omp         dTracerDx, use_Redi_diag_terms)
      do iCell = 1, nCells
         if (maxLevelCell(iCell) .ge. config_Redi_min_layers_diag_terms) then
            use_Redi_diag_terms = .true.
         else
            use_Redi_diag_terms = .false.
         endif

         redi_term1(:, :, iCell) = 0.0_RKIND
         redi_term2(:, :, iCell) = 0.0_RKIND
         redi_term3(:, :, iCell) = 0.0_RKIND
         fluxRediZTop(:, :) = 0.0_RKIND
        !rediKappaCell(iCell) = 0.0_RKIND
         invAreaCell = 1.0_RKIND/areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            ! Check if neighboring cell exists
            if (cellsOnCell(i, iCell) .eq. nCellsP1) cycle
            iEdge = edgesOnCell(i, iCell)
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            if (cell1 == iCell) then
               iCellSelf = 1
            else  ! cell2 == iCell
               iCellSelf = 2
            endif

            r_tmp = dvEdge(iEdge)/dcEdge(iEdge)
            coef = dvEdge(iEdge)

            k = minLevelEdgeBot(iEdge)
            kappaRediEdgeVal = 0.25_RKIND*(RediKappaScaling(k, cell1) + RediKappaScaling(k, cell2) + &
                                           RediKappaScaling(k + 1, cell1) + RediKappaScaling(k + 1, cell2))* &
                               0.25_RKIND*(RediKappaSfcTaper(k, cell1) + RediKappaSfcTaper(k, cell2) + &
                                           RediKappaSfcTaper(k + 1, cell1) + RediKappaSfcTaper(k + 1, cell2))
            k = 1
            kappaRediEdgeVal = 0.25_RKIND*(RediKappaScaling(k, cell1) + RediKappaScaling(k, cell2) + &
                                           RediKappaScaling(k + 1, cell1) + RediKappaScaling(k + 1, cell2))* &
                               0.25_RKIND*(RediKappaSfcTaper(k, cell1) + RediKappaSfcTaper(k, cell2) + &
                                           RediKappaSfcTaper(k + 1, cell1) + RediKappaSfcTaper(k + 1, cell2))

            do iTr = 1, nTracers
               tracer_turb_flux = tracers(iTr, k, cell2) - tracers(iTr, k, cell1)
               flux = layerThicknessEdge(k, iEdge)*tracer_turb_flux*r_tmp*RediKappa(iEdge)
               redi_term1(iTr, k, iCell) = redi_term1(iTr, k, iCell) - (1.0_RKIND + term1TaperFactor* &
                                           (kappaRediEdgeVal - 1.0_RKIND))*edgeSignOnCell(i, iCell)*flux*invAreaCell

               if (.not.use_Redi_diag_terms) cycle

               ! Term 2: div( h S dphi/dz)
               flux_term2 = coef*kappaRediEdgeVal*RediKappa(iEdge)*layerThicknessEdge(k, iEdge)* &
                            0.25_RKIND* &
                            (slopeTriadDown(k, 1, iEdge)*(tracers(iTr, k, cell1) - tracers(iTr, k + 1, cell1)) &
                             /(zMid(k, cell1) - zMid(k + 1, cell1)) &
                             + slopeTriadDown(k, 2, iEdge)*(tracers(iTr, k, cell2) - tracers(iTr, k + 1, cell2)) &
                             /(zMid(k, cell2) - zMid(k + 1, cell2)))
               if (minLevelCell(cell1) < minLevelEdgeBot(iEdge)) then
                  flux_term2 = flux_term2 + coef*RediKappa(iEdge)* &
                               layerThicknessEdge(minLevelEdgeBot(iEdge), iEdge)*0.25_RKIND*kappaRediEdgeVal*slopeTriadUp(minLevelEdgeBot(iEdge), 1, iEdge)* &
                               (tracers(iTr, minLevelEdgeBot(iEdge) - 1, cell1) - tracers(iTr, minLevelEdgeBot(iEdge), cell1))/(zMid(minLevelEdgeBot(iEdge) - 1, cell1) - zMid(minLevelEdgeBot(iEdge), cell1))
               endif
               if (minLevelCell(cell2) < minLevelEdgeBot(iEdge)) then
                  flux_term2 = flux_term2 + coef*RediKappa(iEdge)* &
                               layerThicknessEdge(minLevelEdgeBot(iEdge), iEdge)*0.25_RKIND*kappaRediEdgeVal*slopeTriadUp(minLevelEdgeBot(iEdge), 2, iEdge)* &
                               (tracers(iTr, minLevelEdgeBot(iEdge) - 1, cell2) - tracers(iTr, minLevelEdgeBot(iEdge), cell2))/(zMid(minLevelEdgeBot(iEdge) - 1, cell2) - zMid(minLevelEdgeBot(iEdge), cell2))
               endif

               redi_term2(iTr, k, iCell) = redi_term2(iTr, k, iCell) - edgeSignOnCell(i, iCell)*flux_term2*invAreaCell
               redi_term2_edge(iTr, k, iEdge) = -flux_term2
            end do

            do k = minLevelEdgeBot(iEdge) + 1, maxLevelEdgeTop(iEdge) - 1

               kappaRediEdgeVal = 0.25_RKIND*(RediKappaScaling(k, cell1) + RediKappaScaling(k, cell2) + &
                                              RediKappaScaling(k + 1, cell1) + RediKappaScaling(k + 1, cell2))*0.25_RKIND* &
                                  (RediKappaSfcTaper(k, cell1) + RediKappaSfcTaper(k, cell2) + &
                                   RediKappaSfcTaper(k + 1, cell1) + RediKappaSfcTaper(k + 1, cell2))

               do iTr = 1, nTracers
                  ! \kappa_2 \nabla \phi on edge
                  tracer_turb_flux = tracers(iTr, k, cell2) - tracers(iTr, k, cell1)

                  ! div(h \kappa_2 \nabla \phi) at cell center
                  flux = layerThicknessEdge(k, iEdge)*tracer_turb_flux*r_tmp*RediKappa(iEdge)

                  redi_term1(iTr, k, iCell) = redi_term1(iTr, k, iCell) - (1.0_RKIND + term1TaperFactor* &
                                              (kappaRediEdgeVal - 1.0_RKIND))*edgeSignOnCell(i, iCell)*flux*invAreaCell

                  if (.not.use_Redi_diag_terms) cycle

                  flux_term2 = coef*RediKappa(iEdge)*kappaRediEdgeVal*layerThicknessEdge(k, iEdge)* &
                               0.25_RKIND* &
                               (slopeTriadUp(k, 1, iEdge)*(tracers(iTr, k - 1, cell1) - tracers(iTr, k, cell1)) &
                                /(zMid(k - 1, cell1) - zMid(k, cell1)) &
                                + slopeTriadUp(k, 2, iEdge)*(tracers(iTr, k - 1, cell2) - tracers(iTr, k, cell2)) &
                                /(zMid(k - 1, cell2) - zMid(k, cell2)) &
                                + slopeTriadDown(k, 1, iEdge)*(tracers(iTr, k, cell1) - tracers(iTr, k + 1, cell1)) &
                                /(zMid(k, cell1) - zMid(k + 1, cell1)) &
                                + slopeTriadDown(k, 2, iEdge)*(tracers(iTr, k, cell2) - tracers(iTr, k + 1, cell2)) &
                                /(zMid(k, cell2) - zMid(k + 1, cell2)))
                  redi_term2(iTr, k, iCell) = redi_term2(iTr, k, iCell) - edgeSignOnCell(i, iCell)*flux_term2* &
                                invAreaCell
                  redi_term2_edge(iTr, k, iEdge) = -flux_term2

                  dTracerDx = (tracers(iTr, k, cell2) - tracers(iTr, k, cell1))*dvEdge(iEdge)*0.25_RKIND
                  fluxRediZTop(iTr, k) = fluxRediZTop(iTr, k) + &
                                         0.5_RKIND*(slopeTriadUp(k, iCellSelf, iEdge) + slopeTriadDown(k - 1, iCellSelf, iEdge)) &
                                         *dTracerDx
               end do
            end do

            k = maxLevelEdgeTop(iEdge)
            kappaRediEdgeVal = 0.25_RKIND*(RediKappaScaling(k, cell1) + RediKappaScaling(k, cell2) + &
                                           RediKappaScaling(k + 1, cell1) + RediKappaScaling(k + 1, cell2))* &
                               0.25_RKIND*(RediKappaSfcTaper(k, cell1) + RediKappaSfcTaper(k, cell2) + &
                                           RediKappaSfcTaper(k + 1, cell1) + RediKappaSfcTaper(k + 1, cell2))

            do iTr = 1, nTracers
               tracer_turb_flux = tracers(iTr, k, cell2) - tracers(iTr, k, cell1)
               flux = layerThicknessEdge(k, iEdge)*tracer_turb_flux*r_tmp*RediKappa(iEdge)
               redi_term1(iTr, k, iCell) = redi_term1(iTr, k, iCell) - (1.0_RKIND + term1TaperFactor* &
                                           (kappaRediEdgeVal - 1.0_RKIND))*edgeSignOnCell(i, iCell)*flux*invAreaCell

               if (.not.use_Redi_diag_terms) cycle

               ! For bottom layer, only use triads pointing up:
               flux_term2 = coef*kappaRediEdgeVal*RediKappa(iEdge)*layerThicknessEdge(k, iEdge)* &
                            0.25_RKIND* &
                            (slopeTriadUp(k, 1, iEdge)*(tracers(iTr, k - 1, cell1) - tracers(iTr, k, cell1)) &
                             /(zMid(k - 1, cell1) - zMid(k, cell1)) &
                             + slopeTriadUp(k, 2, iEdge)*(tracers(iTr, k - 1, cell2) - tracers(iTr, k, cell2)) &
                             /(zMid(k - 1, cell2) - zMid(k, cell2)))

               if (maxLevelCell(cell1) > maxLevelEdgeTop(iEdge)) then
                  flux_term2 = flux_term2 + coef*RediKappa(iEdge)* &
                               layerThicknessEdge(k, iEdge)*0.25_RKIND*kappaRediEdgeVal*slopeTriadDown(k, 1, iEdge)* &
                               (tracers(iTr, k, cell1) - tracers(iTr, k + 1, cell1))/(zMid(k, cell1) - zMid(k + 1, cell1))
               endif
               if (maxLevelCell(cell2) > maxLevelEdgeTop(iEdge)) then
                  flux_term2 = flux_term2 + coef*RediKappa(iEdge)* &
                               layerThicknessEdge(k, iEdge)*0.25_RKIND*kappaRediEdgeVal*slopeTriadDown(k, 2, iEdge)* &
                               (tracers(iTr, k, cell2) - tracers(iTr, k + 1, cell2))/(zMid(k, cell2) - zMid(k + 1, cell2))
               endif
               redi_term2(iTr, k, iCell) = redi_term2(iTr, k, iCell) - edgeSignOnCell(i, iCell)*flux_term2 &
                                *invAreaCell
               redi_term2_edge(iTr, k, iEdge) = -flux_term2

               dTracerDx = (tracers(iTr, k, cell2) - tracers(iTr, k, cell1))*dvEdge(iEdge)*0.25_RKIND
               fluxRediZTop(iTr, k) = fluxRediZTop(iTr, k) + &
                                      0.5_RKIND*(slopeTriadUp(k, iCellSelf, iEdge) + slopeTriadDown(k - 1, iCellSelf, iEdge)) &
                                      *dTracerDx
            end do
         end do ! nEdgesOnCell(iCell)

         if ( use_Redi_diag_terms) then
            ! impose no-flux boundary conditions at top and bottom of column
            fluxRediZTop(:, 1:minLevelCell(iCell)) = 0.0_RKIND
            fluxRediZTop(:, maxLevelCell(iCell) + 1) = 0.0_RKIND
            redi_term3_topOfCell(:, 1:minLevelCell(iCell), iCell) = 0.0_RKIND
            do k = minLevelCell(iCell) + 1, maxLevelCell(iCell)
               redi_term3_topOfCell(:, k, iCell) = fluxRediZTop(:, k) * invAreaCell
            end do
            redi_term3_topOfCell(:, maxLevelCell(iCell) + 1, iCell) = 0.0_RKIND

            do k = minLevelCell(iCell), maxLevelCell(iCell)
               do iTr = 1, nTracers
                  ! Add tendency for Term 3: d/dz ( h S grad phi) = ( S grad phi) fluxes
                  ! 2.0 in next line is because a dot product on a C-grid
                  ! requires a factor of 1/2 to average to the cell center.
                  flux_term3 = rediKappaCell(iCell)*2.0_RKIND* &
                               (RediKappaSfcTaper(k, iCell)*RediKappaScaling(k, iCell)*fluxRediZTop(iTr, k)  &
                               * invAreaCell - RediKappaSfcTaper(k + 1, iCell)*RediKappaScaling(k + 1, iCell) &
                               *fluxRediZTop(iTr, k + 1) * invAreaCell)

                  redi_term3(iTr, k, iCell) = redi_term3(iTr, k, iCell) + flux_term3
               end do
            end do
         end if

      end do ! iCell
      !$omp end do
      !$omp end parallel

      !loop over cells and check for out of bounds
      if (isActiveTracer) then
         minimumVal(1) = -2.0_RKIND
         minimumVal(2:nTracers) = 0.0_RKIND
         !$omp parallel
         !$omp do schedule(runtime)
         do iCell = 1, nCells
            rediLimiterCount(:, iCell) = 0
         end do
         !$omp end do
         !$omp end parallel
      else
         minimumVal(:) = 0.0_RKIND
      endif

      !$omp parallel
      !$omp do schedule(runtime) private(k, iTr, tempTracer, iEdge)
      do iCell = 1, nCells
         do k = minLevelCell(iCell), maxLevelCell(iCell)
            do iTr = 1, ntracers
               tempTracer = tracers(iTr, k, iCell) + dt*(redi_term1(iTr, k, iCell) + redi_term2(iTr, k, iCell) + &
                                                         redi_term3(iTr, k, iCell))/layerThickness(k, iCell)

               if (tempTracer < minimumVal(iTr)) then

                  do i = 1, nEdgesOnCell(iCell)
                     iEdge = edgesOnCell(i, iCell)
                     redi_term2_edge(iTr, k, iEdge) = 0.0_RKIND
                  end do

                  redi_term3_topOfCell(iTr, k, iCell) = 0.0_RKIND
                  redi_term3_topOfCell(iTr, k + 1, iCell) = 0.0_RKIND
                  if (isActiveTracer) then
                     rediLimiterCount(k, iCell) = rediLimiterCount(k, iCell) + 1
                  endif
               endif
            end do
         end do
      end do ! iCell
      !$omp end do
      !$omp end parallel

      !now go back and reapply all tendencies
      !$omp parallel
      !$omp do schedule(runtime) private(invAreaCell, k, iTr, i, iEdge)
      do iCell = 1, nCells
         invAreaCell = 1.0_RKIND/areaCell(iCell)
         do k = minLevelCell(iCell), maxLevelCell(iCell)
            do iTr = 1, ntracers
               tend(iTr, k, iCell) = tend(iTr, k, iCell) + redi_term1(iTr, k, iCell)
            end do
         end do

         if (maxLevelCell(iCell) .ge. config_Redi_min_layers_diag_terms) then
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               do iTr = 1, ntracers
                  tend(iTr, k, iCell) = tend(iTr, k, iCell) + rediKappaCell(iCell)*2.0_RKIND* &
                                        (RediKappaSfcTaper(k, iCell)*RediKappaScaling(k, iCell)* &
                                         redi_term3_topOfCell(iTr, k, iCell) - RediKappaSfcTaper(k + 1, iCell) &
                                         *RediKappaScaling(k + 1, iCell)*redi_term3_topOfCell(iTr, k + 1, iCell))
               end do
            end do

            do i = 1, nEdgesOnCell(iCell)
               if (cellsOnCell(i, iCell) .eq. nCellsP1) cycle
               iEdge = edgesOnCell(i, iCell)
               do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                  do iTr = 1, ntracers
                     tend(iTr, k, iCell) = tend(iTr, k, iCell) + edgeSignOnCell(i, iCell)*invAreaCell*redi_term2_edge(iTr, k, iEdge)
                  end do
               end do
            end do
         end if
      end do ! iCell
      !$omp end do
      !$omp end parallel

      deallocate (fluxRediZTop)
      deallocate (minimumVal)
      deallocate (rediKappaCell)
      deallocate (redi_term1)
      deallocate (redi_term2)
      deallocate (redi_term3)
      deallocate (redi_term2_edge)
      deallocate (redi_term3_topOfCell)

      call mpas_timer_stop("tracer redi")

      !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_Redi_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_hmix_Redi_init
!
!> \brief   Initializes ocean tracer horizontal mixing quantities
!> \details
!>  This routine initializes a variety of quantities related to
!>  Laplacian horizontal velocity mixing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hmix_Redi_init(domain, err)!{{{

      type(domain_type), intent(in) :: domain

      integer, intent(out) :: err !< Output: error flag

      type(block_type), pointer :: block
      type(mpas_pool_type), pointer :: meshPool
      type(mpas_pool_type), pointer :: diagnosticsPool
      type(mpas_pool_type), pointer :: forcingPool

      real(kind=RKIND), dimension(:), pointer :: RediKappa, RediKappaData
      real(kind=RKIND), dimension(:), pointer :: dcEdge
      integer, dimension(:,:), pointer :: cellsOnEdge

      integer :: k, iEdge
      integer, pointer :: nEdges

      err = 0

      if (.not. config_use_Redi) return

      if (config_Redi_N2_based_taper_limit_term1) then
         term1TaperFactor = 1.0_RKIND
      else
         term1TaperFactor = 0.0_RKIND
      endif

      block => domain%blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block%structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block%structs, 'forcing', forcingPool)
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
         call mpas_pool_get_array(diagnosticsPool, 'RediKappa', RediKappa)
         call mpas_pool_get_array(forcingPool, 'RediKappaData', RediKappaData)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
         ! initialize Redi kappa array
         if (config_Redi_closure == 'constant') then
            !$omp parallel
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               RediKappa(iEdge) = config_Redi_constant_kappa
            end do
            !$omp end do
            !$omp end parallel
          else if (config_Redi_closure == 'data') then
            ! read RediKappa in from input
            call mpas_log_write("config_Redi_closure = 'data'. "// &
                                "Make sure that the variable RediKappa is read in from an input file.")
         else
            call mpas_log_write('Invalid choice of config_Redi_closure.', MPAS_LOG_CRIT)
            err = 1
            call mpas_dmpar_finalize(domain%dminfo)
         end if

         ! Add resolution taper
         if (config_eddying_resolution_taper == 'none') then
            ! Nothing to do, as we just keep the same assignment as above.
         else if (config_eddying_resolution_taper == 'ramp') then
            !$omp parallel
            !$omp do schedule(runtime) 
            do iEdge=1,nEdges
               if (dcEdge(iEdge) <= config_eddying_resolution_ramp_min) then
                  RediKappa(iEdge) = 0.0_RKIND
               else if (dcEdge(iEdge) >= config_eddying_resolution_ramp_max) then
                  ! nothing to do, i.e. RediKappa(iCell) = RediKappa(iCell)
               else
                  RediKappa(iEdge) = RediKappa(iEdge) &
                                     *(dcEdge(iEdge) - config_eddying_resolution_ramp_min) &
                                     /(config_eddying_resolution_ramp_max - config_eddying_resolution_ramp_min)
               end if
            end do ! iCell
            !$omp end do
            !$omp end parallel
         else
            call mpas_log_write('Invalid choice of config_eddying_resolution_taper.', MPAS_LOG_CRIT)
            err = 1
            call mpas_dmpar_finalize(domain%dminfo)
         end if

         block => block%next
      end do

      !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_Redi_init!}}}

!***********************************************************************

end module ocn_tracer_hmix_Redi

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
