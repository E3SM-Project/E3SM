! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_hmix_del2
!
!> \brief MPAS ocean horizontal tracer mixing driver
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the main driver routine for computing
!>  horizontal mixing tendencies.
!>
!>  It provides an init and a tend function. Each are described below.
!
!-----------------------------------------------------------------------

module ocn_vvel_hmix_del2

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

   public :: ocn_vvel_hmix_del2_tend, &
             ocn_vvel_hmix_del2_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: del2On
   real (kind=RKIND) :: eddyVisc2
   real (kind=RKIND), dimension(:,:), allocatable :: eddyViscArr


!***********************************************************************

contains

!***********************************************************************
!
!  routine vvel_hmix_del2_tend
!
!> \brief   Computes Laplacian tendency term for horizontal vertical velocity mixing
!> \author  Luke Van Roekel -- Copied from tracer_hmix_del2
!> \date    December 2020
!> \details
!>  This routine computes the horizontal mixing tendency for vertical velocity
!>  based on current state using a Laplacian parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_vvel_hmix_del2_tend(layerThicknessEdge, verticalVelocity, tend, err)!{{{

      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdge !< Input: thickness at edges

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        verticalVelocity !< Input: tracer quantities

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

      err = 0

      if (.not.del2On) return

      call mpas_timer_start("vertical velocity del2")

      nCells = nCellsOwned
      nEdges = nEdgesHalo( 1 )

      allocate(eddyViscArr(nVertLevels+1, nEdges))

      !
      ! compute a boundary mask to enforce insulating boundary conditions in the horizontal
      !
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(invAreaCell, i, iEdge, cell1, cell2, r_tmp, k, &
      !$omp         vvel_turb_flux, flux)
      do iCell = 1, nCells
        invAreaCell = 1.0_RKIND / areaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)

          r_tmp = meshScalingDel2(iEdge) * eddyVisc2 * dvEdge(iEdge) / dcEdge(iEdge)

          do k = 2, maxLevelEdgeTop(iEdge)
              ! \kappa_2 \nabla \phi on edge
              vvel_turb_flux = verticalVelocity(k, cell2) - verticalVelocity(k, cell1)

              ! div(h \kappa_2 \nabla \phi) at cell center
              flux = 0.5_RKIND*(layerThicknessEdge(k, iEdge) + layerThicknessEdge(k-1,iEdge)) &
                        * vvel_turb_flux * r_tmp

              tend(k, iCell) = tend(k, iCell) - edgeSignOnCell(i, iCell) * flux * invAreaCell
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel

      deallocate(eddyViscArr)
      call mpas_timer_stop("vertical velocity del2")

   !--------------------------------------------------------------------

   end subroutine ocn_vvel_hmix_del2_tend!}}}

!***********************************************************************
!
!  routine ocn_vvel_hmix_del2_init
!
!> \brief   Initializes ocean vertical velocity horizontal mixing quantities
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine initializes a variety of quantities related to
!>  Laplacian horizontal velocity mixing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vvel_hmix_del2_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      del2on = .false.

      !if ( config_use_vertMom_del2 .or. config_use_two_equation_turbulence_model ) then
      if ( config_use_vertMom_del2 ) then
         if ( config_vertMom_del2 > 0.0_RKIND ) then
            del2On = .true.
            eddyVisc2 = config_vertMom_del2
         endif
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_vvel_hmix_del2_init!}}}

!***********************************************************************

end module ocn_vvel_hmix_del2

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
