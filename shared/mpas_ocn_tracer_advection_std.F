! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_advection_std
!
!> \brief MPAS standard tracer advection
!> \author Doug Jacobsen, Phil Jones
!> \date   03/09/12, updated May 2019
!> \details
!>  This module contains routines for advection of tracers using a standard
!>  FV algorithm in MPAS discretization.
!
!-------------------------------------------------------------------------------

module ocn_tracer_advection_std

   ! module includes
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_io_units
   use mpas_threading

   use mpas_tracer_advection_helpers

   implicit none
   private
   save

   ! private module variables
   real (kind=RKIND) :: &
      coef3rdOrder       !< coefficient for blending high-order terms

   integer :: vertOrder  !< choice of order for vertical advection
   integer, parameter :: &! enumerator for supported vertical adv order
      vertOrder2=2,      &!< 2nd order
      vertOrder3=3,      &!< 3rd order
      vertOrder4=4        !< 4th order

   logical :: monotonicityCheck !< flag to check monotonicity

   ! public method interfaces
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
                                             tend, maxLevelCell, maxLevelEdgeTop, &
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
      real (kind=RKIND), dimension(:,:), allocatable :: tracer_cur, high_order_horiz_flux, high_order_vert_flux

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

      allocate(high_order_horiz_flux(nVertLevels,nEdges))
      allocate(tracer_cur(nVertLevels,nCells+1))
      allocate(high_order_vert_flux(nVertLevels+1,nCells))

      call mpas_threading_barrier()

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
             select case (vertOrder)
             case (vertOrder4)
               high_order_vert_flux(k, iCell) = mpas_tracer_advection_vflux4( tracer_cur(k-2,iCell),tracer_cur(k-1,iCell),  &
                                      tracer_cur(k  ,iCell),tracer_cur(k+1,iCell), w(k,iCell))
             case (vertOrder3)
               high_order_vert_flux(k, iCell) = mpas_tracer_advection_vflux3( tracer_cur(k-2,iCell),tracer_cur(k-1,iCell),  &
                                      tracer_cur(k  ,iCell),tracer_cur(k+1,iCell), w(k,iCell), coef3rdOrder )
             case (vertOrder2)
               verticalWeightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
               verticalWeightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
               high_order_vert_flux(k,iCell) = w(k, iCell) * (verticalWeightK * tracer_cur(k, iCell) &
                                             + verticalWeightKm1 * tracer_cur(k-1, iCell))
             end select ! vertOrder
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
              tracer_weight = highOrderAdvectionMask(k, iEdge) * (adv_coefs(i,iEdge) + coef3rdOrder &
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

      deallocate(tracer_cur, high_order_horiz_flux, high_order_vert_flux)
      deallocate(verticalDivergenceFactor)

   end subroutine ocn_tracer_advection_std_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_tracer_advection_std_init
!
!> \brief MPAS initialize standard tracer advection tendency.
!> \author Doug Jacobsen, Phil Jones
!> \date   03/09/12, updated May 2019
!> \details
!>  This routine initializes constants and choices for the standard tracer
!>  advection tendency.
!
!-------------------------------------------------------------------------------

   subroutine ocn_tracer_advection_std_init(horzAdvOrder, vertAdvOrder, &
                                            inCoef3rdOrder,             &
                                            checkMonotonicity, err) !{{{

      !*** input parameters

      integer, intent(in) :: &
         horzAdvOrder         !< [in] Order for horizontal advection
      integer, intent(in) :: &
         vertAdvOrder         !< [in] Order for vertical advection
      real (kind=RKIND), intent(in) :: &
         inCoef3rdOrder       !< [in] Coefficient for blending advection orders
      logical, intent(in) :: &
         checkMonotonicity    !< [in] Flag to check on monotonicity of tracers

      !*** output parameters

      integer, intent(out) :: err !< [out] Error Flag

      ! end of preamble
      !----------------
      ! begin code

      err = 0 ! set error code to success

      ! set monotonicity flag
      monotonicityCheck = checkMonotonicity

      ! set 3rd order coefficient based on horizontal order choice
      select case (horzAdvOrder)
      case (2)
         coef3rdOrder = 0.0_RKIND
      case (3)
         coef3rdOrder = inCoef3rdOrder
      case (4)
         coef3rdOrder = 0.0_RKIND
      case default
         call mpas_log_write( &
            'Invalid value for horz advection order, defaulting to 2nd order', &
            MPAS_LOG_WARN)
      end select ! horzAdvOrder

      ! set choice of vertical advection order
      select case (vertAdvOrder)
      case (2)
         vertOrder = vertOrder2
      case (3)
         vertOrder = vertOrder3
      case (4)
         vertOrder = vertOrder4
      case default
         vertOrder = vertOrder2
         call mpas_log_write( &
         'Invalid value for vertical advection order, defaulting to 2nd order',&
         MPAS_LOG_WARN)
      end select ! vertAdvOrder

   end subroutine ocn_tracer_advection_std_init!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

end module ocn_tracer_advection_std

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
