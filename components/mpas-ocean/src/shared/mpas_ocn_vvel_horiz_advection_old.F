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

module ocn_vvel_horiz_advection

   ! module includes
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_io_units
   use mpas_threading
   use ocn_mesh
   use ocn_config
   use mpas_tracer_advection_helpers
   use ocn_tracer_advection_shared

   implicit none
   private
   save

   ! private module variables
   real (kind=RKIND) :: &
      coef3rdOrder       !< coefficient for blending high-order terms

   real (kind=RKIND) :: advCoeffOrderVal

   ! public method interfaces
   public :: ocn_vvel_horiz_advection_tend, &
             ocn_vvel_horiz_advection_init

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
   subroutine ocn_vvel_horiz_advection_tend(verticalVelocity, normalTransportVelocity, &!{{{
                                             layerThickness, dt, tend)

      real (kind=RKIND), dimension(:,:), intent(in) :: normalTransportVelocity !< Input: Thichness weighted velocitiy
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness !< Input: Thickness
      real (kind=RKIND), dimension(:,:), intent(in) :: verticalVelocity
      real (kind=RKIND), intent(in) :: dt !< Input: Timestep
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend !< Input/Output: Vertical Velocity tendency

      integer :: i, iCell, iEdge, k, cell1, cell2
      integer :: nVertLevels, nCells, nEdges

      real (kind=RKIND) :: vvel_weight, invAreaCell1

      real (kind=RKIND), dimension(:,:), allocatable :: high_order_horiz_flux, tendtmp

      real (kind=RKIND), parameter :: eps = 1.e-10_RKIND

      ! Get dimensions
     nVertLevels = size(verticalVelocity,dim=1)

      ! Initialize pointers
      nCells = nCellsHalo( 1 )
      nEdges = nEdgesHalo( 2 )
      allocate(high_order_horiz_flux(nVertLevels, nEdges))
      allocate(tendtmp(nVertLevels,nCells))

       !$omp do schedule(runtime)
        do iEdge = 1, nEdges
           high_order_horiz_flux(:, iEdge) = 0.0_RKIND
        end do
        !$omp end do

        !$omp do schedule(runtime)
        do iCell = 1,nCells
           tendtmp(:, iCell) = tend(:, iCell)
        end do
        !$omp end do

        !FIXME -- will need to update for boundaries, maybe add edgeMask
        !  Compute the high order horizontal flux
        !$omp do schedule(runtime) private(cell1, cell2, k, vvel_weight, i, iCell)
        do iEdge = 1, nEdges
          cell1 = cellsOnEdge(1, iEdge)
          cell2 = cellsOnEdge(2, iEdge)

          ! Compute 2nd order fluxes where needed.
          do k = 2, maxLevelEdgeTop(iEdge)
            vvel_weight = advCoeffOrderVal * (dvEdge(iEdge) * 0.5_RKIND) &
                           * 0.5_RKIND*(normalTransportVelocity(k-1, iEdge) + normalTransportVelocity(k, iEdge))

           vvel_weight = advCoeffOrderVal * 0.5_RKIND*dvEdge(iEdge)*0.5_RKIND*(normalTransportVelocity(k-1, iEdge) + normalTransportVelocity(k, iEdge))
            high_order_horiz_flux(k, iEdge) = high_order_horiz_flux(k, iedge) + vvel_weight &
                                            * (verticalVelocity(k, cell1) + verticalVelocity(k, cell2))
          end do ! k loop

          if(advCoeffOrderVal == 0.0_RKIND) then
          ! Compute 3rd or 4th fluxes where requested.
          do i = 1, nAdvCellsForEdge(iEdge)
            iCell = advCellsForEdge(i,iEdge)
            do k = 2, maxLevelCell(iCell)
              vvel_weight = (1.0_RKIND - advCoeffOrderVal) * (advCoefs(i,iEdge) + coef3rdOrder &
                            * sign(1.0_RKIND,0.5_RKIND*(normalTransportVelocity(k-1,iEdge) +          &
                            normalTransportVelocity(k,iEdge)))*advCoefs3rd(i,iEdge))

              vvel_weight = 0.5*(normalTransportVelocity(k-1,iEdge) + normalTransportVelocity(k,iEdge))*vvel_weight
              high_order_horiz_flux(k,iEdge) = high_order_horiz_flux(k,iEdge) + vvel_weight * verticalVelocity(k,iCell)
            end do ! k loop
          end do ! i loop over nAdvCellsForEdge
        endif

        end do ! iEdge loop
        !$omp end do

        ! Accumulate the scaled high order horizontal tendencies
        !$omp do schedule(runtime) private(invAreaCell1, i, iEdge, k)
        do iCell = 1, nCells
          invAreaCell1 = 1.0_RKIND / areaCell(iCell)
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 2, maxLevelEdgeTop(iEdge)-1
              tend(k, iCell) = tend(k, iCell) + edgeSignOnCell(i, iCell) * high_order_horiz_flux(k, iEdge) &
                                      * invAreaCell1
              tendtmp(k,iCell) = tendtmp(k,iCell) - tend(k,iCell)
            end do
          end do
        end do
        !$omp end do

      deallocate(high_order_horiz_flux)
      deallocate(tendtmp)

   end subroutine ocn_vvel_horiz_advection_tend!}}}

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

   subroutine ocn_vvel_horiz_advection_init(horzAdvOrder, &
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

   end subroutine ocn_vvel_horiz_advection_init!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

end module ocn_vvel_horiz_advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
