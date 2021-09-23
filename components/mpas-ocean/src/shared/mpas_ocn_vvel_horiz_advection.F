! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vvel_advection_std
!
!> \brief MPAS standard vertical velocity horizontal advection
!> \author Luke Van Roekel
!> \date   January 2021
!> \details
!>  This module contains routines for horizontal advection of vertical velocity
!>  using a standard FV algorithm in MPAS discretization.  Based on the standard
!>  tracer advection scheme
!
!-------------------------------------------------------------------------------

module ocn_vvel_horiz_advection_vector_invariant

   ! module includes
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_io_units
   use mpas_threading
   use ocn_config
   use mpas_tracer_advection_helpers
   use ocn_diagnostics_variables

   implicit none
   private
   save

   ! private module variables
   real (kind=RKIND) :: &
      coef3rdOrder       !< coefficient for blending high-order terms

   real (kind=RKIND) :: advCoeffOrderVal

   ! public method interfaces
   public :: ocn_vvel_horiz_vector_invariant_advection_tend, &
             ocn_vvel_horiz_vector_invariant_advection_init

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_vvel_horiz_advection_std_tend
!
!> \brief MPAS standard vertical velocity advection tendency
!> \author Luke Van Roekel
!> \date   12/21/20
!> \details
!>  This routine computes the standard vertical velocity advection tendencity.
!>  NOTE this is only active for non hydrostatic mode
!
!-----------------------------------------------------------------------
   subroutine ocn_vvel_horiz_vector_invariant_advection_tend(verticalVelocity, adv_coefs, adv_coefs_3rd, nAdvCellsForEdge, &
                                            advCellsForEdge, &!{{{
                                            normalTransportVelocity, layerThicknessEdge, dt, meshPool, &
                                             tend, maxLevelCell, maxLevelEdgeTop, &
                                             highOrderAdvectionMask, edgeSignOnCell)

      real (kind=RKIND), dimension(:,:), intent(in) :: adv_coefs !< Input: Advection coefficients for 2nd order advection
      real (kind=RKIND), dimension(:,:), intent(in) :: adv_coefs_3rd !< Input: Advection coeffs for blending in 3rd/4th order
      integer, dimension(:), intent(in) :: nAdvCellsForEdge !< Input: Number of advection cells for each edge
      integer, dimension(:,:), intent(in) :: advCellsForEdge !< Input: List of advection cells for each edge
      real (kind=RKIND), dimension(:,:), intent(in) :: normalTransportVelocity !< Input: Thichness weighted velocitiy
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThicknessEdge !< Input: Thickness
      real (kind=RKIND), dimension(:,:), intent(in) :: verticalVelocity
      real (kind=RKIND), intent(in) :: dt !< Input: Timestep
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend !< Input/Output: Vertical Velocity tendency
      integer, dimension(:), pointer :: maxLevelCell !< Input: Index to max level at cell center
      integer, dimension(:), pointer :: maxLevelEdgeTop !< Input: Index to max level at edge with non-land cells on both sides
      integer, dimension(:,:), pointer :: highOrderAdvectionMask !< Input: Mask for high order advection
      integer, dimension(:, :), pointer :: edgeSignOnCell !< Input: Sign for flux from edge on each cell.

      integer :: i, iCell, iEdge, k, cell1, cell2
      integer :: nVertLevels, nCells, nEdges
      integer, pointer :: nCellsSolve, maxEdges
      integer, dimension(:), pointer :: nEdgesOnCell, nCellsArray, nEdgesArray
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnCell, edgesOnCell

      real (kind=RKIND) :: vvel_weight, invAreaCell1
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell, dcEdge

      real (kind=RKIND), dimension(:), allocatable :: fluxTmp

      real (kind=RKIND), dimension(:,:), allocatable :: high_order_horiz_flux, tendtmp

      real (kind=RKIND), parameter :: eps = 1.e-10_RKIND

      ! Get dimensions
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdges)
      nVertLevels = size(verticalVelocity,dim=1)

      ! Initialize pointers
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)

      nCells = nCellsArray( 2 )
      nEdges = nEdgesArray( 3 )
      allocate(high_order_horiz_flux(nVertLevels, nEdges))
      allocate(tendtmp(nVertLevels,nCells))
      allocate(fluxTmp(nVertLevels))

      do iEdge = 1, nEdges
         high_order_horiz_flux(:, iEdge) = 0.0_RKIND
      end do

      !Switch to the vector invariant form, may be buggy.  Likely in the sign.
      !$omp do schedule(runtime) private(cell1, cell2, k)
      do iEdge = 1,nEdges
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)

         do k = 2, maxLevelEdgeTop(iEdge)
            high_order_horiz_flux(k, iEdge) = 0.5_RKIND*(normalTransportVelocity(k-1,iEdge) + &
                  normalTransportVelocity(k,iEdge))*((verticalVelocity(k,cell2) - &
                  verticalVelocity(k,cell1))/dcEdge(iEdge) + (normalTransportVelocity(k-1, iEdge) - &
                  normalTransportVelocity(k,iEdge)) / (0.5_RKIND*(layerThicknessEdge(k-1,iEdge) + &
                  layerThicknessEdge(k,iEdge))))
         end do
      end do
      !$omp end do

      !$omp do schedule(runtime) private(vvel_weight, iEdge, i, k)
      do iCell = 1,nCells
         invAreaCell1 = 1.0_RKIND / areaCell(iCell)
         fluxTmp(:) = 0.0_RKIND
         do k = 2, maxLevelEdgeTop(iEdge)
            do i = 1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
  
               fluxTmp(k) = fluxTmp(k) + high_order_horiz_flux(k,iEdge)*edgeAreaFractionOfCell(i,iCell)* &
                 edgeSignOnCell(i, iCell)
            end do
            tend(k,iCell) = tend(k,iCell) + 2.0_RKIND*fluxTmp(k)
         end do
      end do
      !$omp end do

      deallocate(fluxTmp)
      deallocate(high_order_horiz_flux)
      deallocate(tendtmp)

   end subroutine ocn_vvel_horiz_vector_invariant_advection_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_vvel_horiz_advection_init
!
!> \brief MPAS initialize vertical velocity advection tendency.
!> \author Luke Van Roekel
!> \date   12/21/2020
!> \details
!>  This routine initializes constants and choices for the vertical velocity
!>  advection tendency
!
!-------------------------------------------------------------------------------

   subroutine ocn_vvel_horiz_vector_invariant_advection_init(horzAdvOrder, &
                                            inCoef3rdOrder,             &
                                            err) !{{{

      !*** input parameters

      integer, intent(in) :: &
         horzAdvOrder         !< [in] Order for horizontal advection
      real (kind=RKIND), intent(in) :: &
         inCoef3rdOrder       !< [in] Coefficient for blending advection orders

      !*** output parameters

      integer, intent(out) :: err !< [out] Error Flag

      ! end of preamble
      !----------------
      ! begin code

      err = 0 ! set error code to success

      advCoeffOrderVal = 0.0_RKIND

      ! set 3rd order coefficient based on horizontal order choice
      select case (horzAdvOrder)
      case (2)
         coef3rdOrder = 0.0_RKIND
         advCoeffOrderVal = 1.0_RKIND
      case (3)
         coef3rdOrder = inCoef3rdOrder
      case (4)
         coef3rdOrder = 0.0_RKIND
      case default
         call mpas_log_write( &
            'Invalid value for horz advection order, defaulting to 2nd order', &
            MPAS_LOG_WARN)
      end select ! horzAdvOrder

   end subroutine ocn_vvel_horiz_vector_invariant_advection_init!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

end module ocn_vvel_horiz_advection_vector_invariant

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
