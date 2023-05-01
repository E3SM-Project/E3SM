! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the
! LICENSE file
! distributed with this code, or at
! http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vvel_hmix_del4
!
!> \brief MPAS ocean horizontal mixing driver for vertical momentum
!> \author Darren Engwirda, Luke VanRoekel, Sara Calandrini
!> \date   April 2023
!> \details
!>  This module contains the main driver routine for computing
!>  horizontal mixing tendencies in the vertical momentum.
!>
!>  It provides an init and a tend function. Each are described below.
!
!-----------------------------------------------------------------------

module ocn_vvel_hmix_del4

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_threading

   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_diagnostics_variables

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

   public :: ocn_vvel_hmix_del4_tend, &
             ocn_vvel_hmix_del4_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: del4On
   real (kind=RKIND) :: eddyVisc4


!***********************************************************************

contains

!***********************************************************************
!
!  routine vvel_hmix_del4_tend
!
!> \brief   Computes biharmonic tendency term for horizontal vertical
!velocity mixing
!> \author  Darren Engwirda, Sara Calandrini
!> \date    April 2023
!> \details
!>  This routine computes the horizontal mixing tendency for vertical
!velocity
!>  based on current state using a biharmonic parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_vvel_hmix_del4_tend(layerThicknessEdge, vertVelocityNonhydro, tend, err)!{{{

      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdge !< Input: thickness at edges

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertVelocityNonhydro !< Input: tracer quantities

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

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

      integer :: iCell, iEdge, cell1, cell2
      integer :: i, k, nCells, nEdges

      real (kind=RKIND) :: invAreaCell
      real (kind=RKIND) :: vvel_turb_flux, flux, r_tmp
      real (kind=RKIND), dimension(:, :), allocatable :: del2

      err = 0

      if (.not.del4On) return

      call mpas_timer_start("vertical velocity del4")

      nCells = nCellsAll
      nEdges = nEdgesHalo( 1 )

      allocate(del2(nVertLevels+1, nCells))
      del2(:,:) = 0.0_RKIND

      !
      ! compute del^2(vvel)
      !
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(invAreaCell, i, iEdge, cell1, cell2, r_tmp, k, &
      !$omp         vvel_turb_flux, flux)
      do iCell = 1, nCells
        invAreaCell = 1.0_RKIND / areaCell(iCell)
        !del2(:, iCell) = 0.0_RKIND
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)

          r_tmp = dvEdge(iEdge) / dcEdge(iEdge)

          do k = 2, maxLevelEdgeTop(iEdge)
              ! \kappa_2 \nabla \phi on edge
              vvel_turb_flux = vertVelocityNonhydro(k, cell2) - vertVelocityNonhydro(k, cell1)

              ! div(\kappa_2 \nabla \phi) at cell center
              flux = vvel_turb_flux * r_tmp

              del2(k, iCell) = del2(k, iCell) - edgeSignOnCell(i, iCell) * flux * invAreaCell
          end do
        end do
      end do
      !$omp end do

      !
      ! compute del^2(del2)
      !
      !$omp do schedule(runtime) &
      !$omp private(invAreaCell, i, iEdge, cell1, cell2, r_tmp, k, &
      !$omp         vvel_turb_flux, flux)
      do iCell = 1, nCells
        invAreaCell = 1.0_RKIND / areaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)

          r_tmp = meshScalingDel4(iEdge) * eddyVisc4 * dvEdge(iEdge) / dcEdge(iEdge)

          do k = 2, maxLevelEdgeTop(iEdge)
              ! \kappa_2 \nabla \phi on edge
              vvel_turb_flux = del2(k, cell2) - del2(k, cell1)

              ! div(\kappa_2 \nabla \phi) at cell center
              flux = vvel_turb_flux * r_tmp

              ! NB. +ve for del^4
              tend(k, iCell) = tend(k, iCell) + edgeSignOnCell(i, iCell) * flux * invAreaCell
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel

      deallocate(del2)

      call mpas_timer_stop("vertical velocity del4")

   !--------------------------------------------------------------------

   end subroutine ocn_vvel_hmix_del4_tend!}}}

!***********************************************************************
!
!  routine ocn_vvel_hmix_del4_init
!
!> \brief   Initializes ocean vertical velocity horizontal mixing quantities
!> \author  Darren Engwirda, Sara Calandrini
!> \date    April 2023
!> \details
!>  This routine initializes a variety of quantities related to
!>  biharmonic vertical velocity mixing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vvel_hmix_del4_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      del4On = .false.

      if ( config_use_vertMom_del4 ) then
         if ( config_vertMom_del4 > 0.0_RKIND ) then
            del4On = .true.
            eddyVisc4 = config_vertMom_del4
         endif
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_vvel_hmix_del4_init!}}}

!***********************************************************************

end module ocn_vvel_hmix_del4

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
