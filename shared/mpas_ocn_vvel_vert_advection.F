! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vvel_vert_advection_std
!
!> \brief MPAS standard vertical velocity advection
!> \author Luke Van Roekel
!> \date   January 2021
!> \details
!>  This module contains routines for advection of vertical velocity using a standard
!>  FV algorithm in MPAS discretization.  This routine has been adapted from the
!>  tracer advection code without FCT
!
!-------------------------------------------------------------------------------

module ocn_vvel_vert_advection

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

   integer :: vertOrder  !< choice of order for vertical advection
   integer, parameter :: &! enumerator for supported vertical adv order
      vertOrder2=2,      &!< 2nd order
      vertOrder3=3,      &!< 3rd order
      vertOrder4=4        !< 4th order

   ! public method interfaces
   public :: ocn_vvel_vert_advection_tend, &
             ocn_vvel_vert_advection_init

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_vvel_vert_advection_tend
!
!> \brief MPAS standard vertical velocity vertical advection tendency
!> \author Luke Van Roekel
!> \date   12/22/2020
!> \details
!>  This routine computes the standard vertical velocity vertical advection tendency.
!>
!
!-----------------------------------------------------------------------
   subroutine ocn_vvel_vert_advection_tend(verticalVelocity, &!{{{
                                             layerThickness, dt, meshPool, &
                                             tend, maxLevelCell, maxLevelEdgeTop)

      real (kind=RKIND), dimension(:,:), intent(in) :: verticalVelocity !< Input: current vertVel values
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness !< Input: Thickness
      real (kind=RKIND), intent(in) :: dt !< Input: Timestep
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend !< Input/Output: VertVel tendency
      integer, dimension(:), pointer :: maxLevelCell !< Input: Index to max level at cell center
      integer, dimension(:), pointer :: maxLevelEdgeTop !< Input: Index to max level at edge with non-land cells on both sides

      integer :: i, iCell, iEdge, k, cell1, cell2
      integer :: nVertLevels
      integer, pointer :: nCells, nEdges, maxEdges
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnCell, edgesOnCell

      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real (kind=RKIND), parameter :: eps = 1.e-10_RKIND

      ! Get dimensions
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdges)

      ! Initialize pointers
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)

        !  Compute the high order vertical flux. now to be consistent with Horizontal
! momentum

        !omp do schedule(runtime)
        do iCell = 1, nCells
          do k=2,maxLevelCell(iCell)
            tend(k,iCell) = tend(k,iCell) - 0.5_RKIND*(verticalVelocity(k-1,iCell)**2.0 - &
                              verticalVelocity(k+1,iCell)**2.0) / (layerThickness(k-1,iCell) + &
                              layerThickness(k,iCell))
          end do
        enddo
        !omp end do

!        !$omp do schedule(runtime) private(k, verticalWeightK, verticalWeightKm1)
!        do iCell = 1, nCells
!          k = max(1, min(maxLevelCell(iCell), 2))
!          verticalWeightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
!          verticalWeightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
!          high_order_vert_flux(k,iCell) = 0.5*(verticalVelocity(k,iCell)**2

!          do k=3,maxLevelCell(iCell)-1
!             select case (vertOrder)
!             case (vertOrder4) !need to inline these
!               high_order_vert_flux(k, iCell) = mpas_tracer_advection_vflux4( tracer_cur(k-2,iCell),tracer_cur(k-1,iCell),  &
!                                      tracer_cur(k  ,iCell),tracer_cur(k+1,iCell), w(k,iCell))
!             case (vertOrder3)
!               high_order_vert_flux(k, iCell) = mpas_tracer_advection_vflux3( tracer_cur(k-2,iCell),tracer_cur(k-1,iCell),  &
!                                      tracer_cur(k  ,iCell),tracer_cur(k+1,iCell), w(k,iCell), coef3rdOrder )
!             case (vertOrder2)
!               verticalWeightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
!               verticalWeightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
!               high_order_vert_flux(k,iCell) = w(k, iCell) * (verticalWeightK * tracer_cur(k, iCell) &
!                                             + verticalWeightKm1 * tracer_cur(k-1, iCell))
!             end select ! vertOrder
!          end do

!          k = max(1, maxLevelCell(iCell))
!          verticalWeightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
!          verticalWeightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
!          high_order_vert_flux(k,iCell) = w(k,iCell)*(verticalWeightK*tracer_cur(k,iCell)+verticalWeightKm1*tracer_cur(k-1,iCell))
!        end do ! iCell Loop
!        !$omp end do


end subroutine ocn_vvel_vert_advection_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_vvel_vert_advection_std_init
!
!> \brief MPAS initialize vertical velocity vertical advection tendency.
!> \author Luke Van Roekel
!> \date   January 2021
!> \details
!>  This routine initializes constants and choices for the standard vertical
!>  velocity vertical advection tendency.
!
!-------------------------------------------------------------------------------

   subroutine ocn_vvel_vert_advection_init(vertAdvOrder, &
                                            inCoef3rdOrder, err) !{{{

      !*** input parameters
!FIXME -- NOTE this routine may not be necessary.
      integer, intent(in) :: &
         vertAdvOrder         !< [in] Order for vertical advection
      real (kind=RKIND), intent(in) :: &
         inCoef3rdOrder       !< [in] Coefficient for blending advection orders

      !*** output parameters

      integer, intent(out) :: err !< [out] Error Flag

      ! end of preamble
      !----------------
      ! begin code

      err = 0 ! set error code to success

      ! set choice of vertical advection order NOT FUNCTIONAL YET
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
         'Invalid value for vertical velocity vertical advection order, defaulting to 2nd order',&
         MPAS_LOG_WARN)
      end select ! vertAdvOrder

   end subroutine ocn_vvel_vert_advection_init!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

end module ocn_vvel_vert_advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
