! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_advection_std
!
!> \brief MPAS standard tracer advection
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This module contains routines for standard advection of tracers
!
!-----------------------------------------------------------------------
module ocn_tracer_advection_std

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_io_units
   use mpas_threading

   use mpas_tracer_advection_helpers

   implicit none
   private
   save

   real (kind=RKIND) :: coef_3rd_order
   integer :: horizOrder
   logical :: vert2ndOrder, vert3rdOrder, vert4thOrder
   logical :: positiveDzDk, monotonicityCheck

   public :: ocn_tracer_advection_std_tend, &
             ocn_tracer_advection_std_init

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_tracer_advection_std_tend
!
!> \brief MPAS standard tracer advection tendency
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine computes the standard tracer advection tendencity.
!>  Both horizontal and vertical.
!
!-----------------------------------------------------------------------
   subroutine ocn_tracer_advection_std_tend(tracers, adv_coefs, adv_coefs_3rd, nAdvCellsForEdge, advCellsForEdge, &!{{{
                                             normalThicknessFlux, w, layerThickness, verticalCellSize, dt, meshPool, &
                                             scratchPool, tend, maxLevelCell, maxLevelEdgeTop, &
                                             highOrderAdvectionMask, edgeSignOnCell)

      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers !< Input: current tracer values
      real (kind=RKIND), dimension(:,:), intent(in) :: adv_coefs !< Input: Advection coefficients for 2nd order advection
      real (kind=RKIND), dimension(:,:), intent(in) :: adv_coefs_3rd !< Input: Advection coeffs for blending in 3rd/4th order
      integer, dimension(:), intent(in) :: nAdvCellsForEdge !< Input: Number of advection cells for each edge
      integer, dimension(:,:), intent(in) :: advCellsForEdge !< Input: List of advection cells for each edge
      real (kind=RKIND), dimension(:,:), intent(in) :: normalThicknessFlux !< Input: Thichness weighted velocitiy
      real (kind=RKIND), dimension(:,:), intent(in) :: w !< Input: Vertical velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness !< Input: Thickness
      real (kind=RKIND), dimension(:,:), intent(in) :: verticalCellSize !< Input: Distance between vertical interfaces of a cell
      real (kind=RKIND), intent(in) :: dt !< Input: Timestep
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      type (mpas_pool_type), intent(in) :: scratchPool !< Input: Scratch fields
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tend !< Input/Output: Tracer tendency
      integer, dimension(:), pointer :: maxLevelCell !< Input: Index to max level at cell center
      integer, dimension(:), pointer :: maxLevelEdgeTop !< Input: Index to max level at edge with non-land cells on both sides
      integer, dimension(:,:), pointer :: highOrderAdvectionMask !< Input: Mask for high order advection
      integer, dimension(:, :), pointer :: edgeSignOnCell !< Input: Sign for flux from edge on each cell.

      integer :: i, iCell, iEdge, k, iTracer, cell1, cell2
      integer :: nVertLevels, num_tracers
      integer, pointer :: nCells, nEdges, nCellsSolve, maxEdges
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnCell, edgesOnCell

      real (kind=RKIND) :: tracer_weight, invAreaCell1
      real (kind=RKIND) :: verticalWeightK, verticalWeightKm1
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell, verticalDivergenceFactor
      real (kind=RKIND), dimension(:,:), pointer :: tracer_cur, high_order_horiz_flux, high_order_vert_flux

      type (field2DReal), pointer :: highOrderHorizFluxField, tracerCurField, highOrderVertFluxField

      real (kind=RKIND), parameter :: eps = 1.e-10_RKIND

      ! Get dimensions
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdges)
      nVertLevels = size(tracers,dim=2)
      num_tracers = size(tracers,dim=1)

      ! Initialize pointers
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)

      allocate(verticalDivergenceFactor(nVertLevels))
      verticalDivergenceFactor = 1.0_RKIND

      call mpas_pool_get_field(scratchPool, 'highOrderHorizFlux', highOrderHorizFluxField)
      call mpas_pool_get_field(scratchPool, 'tracerValue', tracerCurField, 1)
      call mpas_pool_get_field(scratchPool, 'highOrderVertFlux', highOrderVertFluxField)

      call mpas_allocate_scratch_field(highOrderHorizFluxField, .true.)
      call mpas_allocate_scratch_field(tracerCurField, .true.)
      call mpas_allocate_scratch_field(highOrderVertFluxField, .true.)
      call mpas_threading_barrier()

      high_order_horiz_flux => highOrderHorizFluxField % array
      tracer_cur => tracerCurField % array
      high_order_vert_flux => highOrderVertFluxField % array

      ! Loop over tracers. One tracer is advected at a time. It is copied into a temporary array in order to improve locality
      do iTracer = 1, num_tracers
        ! Initialize variables for use in this iTracer iteration
        !$omp do schedule(runtime)
        do iCell = 1, nCells
           tracer_cur(:, iCell) = tracers(iTracer, :, iCell)

           high_order_vert_flux(:, iCell) = 0.0_RKIND
        end do
        !$omp end do

        !$omp do schedule(runtime)
        do iEdge = 1, nEdges
           high_order_horiz_flux(:, iEdge) = 0.0_RKIND
        end do
        !$omp end do

        !  Compute the high order vertical flux. Also determine bounds on tracer_cur.
        !$omp do schedule(runtime) private(k, verticalWeightK, verticalWeightKm1)
        do iCell = 1, nCells
          k = max(1, min(maxLevelCell(iCell), 2))
          verticalWeightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
          verticalWeightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
          high_order_vert_flux(k,iCell) = w(k,iCell)*(verticalWeightK*tracer_cur(k,iCell)+verticalWeightKm1*tracer_cur(k-1,iCell))

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
               high_order_vert_flux(k,iCell) = w(k, iCell) * (verticalWeightK * tracer_cur(k, iCell) &
                                             + verticalWeightKm1 * tracer_cur(k-1, iCell))
             end if
          end do

          k = max(1, maxLevelCell(iCell))
          verticalWeightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
          verticalWeightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
          high_order_vert_flux(k,iCell) = w(k,iCell)*(verticalWeightK*tracer_cur(k,iCell)+verticalWeightKm1*tracer_cur(k-1,iCell))
        end do ! iCell Loop
        !$omp end do

        !  Compute the high order horizontal flux
        !$omp do schedule(runtime) private(cell1, cell2, k, tracer_weight, i, iCell)
        do iEdge = 1, nEdges
          cell1 = cellsOnEdge(1, iEdge)
          cell2 = cellsOnEdge(2, iEdge)

          ! Compute 2nd order fluxes where needed.
          do k = 1, maxLevelEdgeTop(iEdge)
            tracer_weight = iand(highOrderAdvectionMask(k, iEdge)+1, 1) * (dvEdge(iEdge) * 0.5_RKIND) &
                           * normalThicknessFlux(k, iEdge)

            high_order_horiz_flux(k, iEdge) = high_order_horiz_flux(k, iedge) + tracer_weight &
                                            * (tracer_cur(k, cell1) + tracer_cur(k, cell2))
          end do ! k loop

          ! Compute 3rd or 4th fluxes where requested.
          do i = 1, nAdvCellsForEdge(iEdge)
            iCell = advCellsForEdge(i,iEdge)
            do k = 1, maxLevelCell(iCell)
              tracer_weight = highOrderAdvectionMask(k, iEdge) * (adv_coefs(i,iEdge) + coef_3rd_order &
                            * sign(1.0_RKIND,normalThicknessFlux(k,iEdge))*adv_coefs_3rd(i,iEdge))

              tracer_weight = normalThicknessFlux(k,iEdge)*tracer_weight
              high_order_horiz_flux(k,iEdge) = high_order_horiz_flux(k,iEdge) + tracer_weight * tracer_cur(k,iCell)
            end do ! k loop
          end do ! i loop over nAdvCellsForEdge
        end do ! iEdge loop
        !$omp end do

        ! Accumulate the scaled high order horizontal tendencies
        !$omp do schedule(runtime) private(invAreaCell1, i, iEdge, k)
        do iCell = 1, nCells
          invAreaCell1 = 1.0_RKIND / areaCell(iCell)
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelEdgeTop(iEdge)
              tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + edgeSignOnCell(i, iCell) * high_order_horiz_flux(k, iEdge) &
                                      * invAreaCell1
            end do
          end do
        end do
        !$omp end do

        ! Accumulate the scaled high order vertical tendencies.
        !$omp do schedule(runtime) private(k)
        do iCell = 1, nCellsSolve
          do k = 1,maxLevelCell(iCell)
            tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + verticalDivergenceFactor(k) * (high_order_vert_flux(k+1, iCell) &
                                    - high_order_vert_flux(k, iCell))
          end do ! k loop
        end do ! iCell loop
        !$omp end do
      end do ! iTracer loop

      call mpas_threading_barrier()
      call mpas_deallocate_scratch_field(highOrderHorizFluxField, .true.)
      call mpas_deallocate_scratch_field(tracerCurField, .true.)
      call mpas_deallocate_scratch_field(highOrderVertFluxField, .true.)

      deallocate(verticalDivergenceFactor)

   end subroutine ocn_tracer_advection_std_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_tracer_advection_std_init
!
!> \brief MPAS initialize standard tracer advection tendency.
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine initializes the standard tracer advection tendencity.
!
!-----------------------------------------------------------------------
   subroutine ocn_tracer_advection_std_init(horiz_adv_order, vert_adv_order, coef_3rd_order_in, dzdk_positive, & !{{{
                                            check_monotonicity, err)
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
             call mpas_log_write( &
                'Invalid value for vert_adv_order, defaulting to 2nd order', &
                MPAS_LOG_WARN)
          end if
      end if

      positiveDzDk = dzdk_positive
      monotonicityCheck = check_monotonicity

   end subroutine ocn_tracer_advection_std_init!}}}

end module ocn_tracer_advection_std

