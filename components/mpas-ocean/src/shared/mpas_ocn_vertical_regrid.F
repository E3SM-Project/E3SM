! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vertical_regrid
!
!> \brief MPAS ocean vertical regridding
!> \author Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis
!> \date   July 2021
!> \details
!>  This module contains the vertical regridding routine, used for
!>  vertical Lagrangian remapping.
!
!-----------------------------------------------------------------------

module ocn_vertical_regrid

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_threading
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer
   use mpas_log

   use ocn_constants
   use ocn_config
   use ocn_mesh

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

   public :: ocn_vert_regrid
   public :: ocn_vert_regrid_init

   contains

!***********************************************************************
!
!  routine ocn_vert_regrid
!
!> \brief   Regridding for vertical Lagrangian remapping
!> \author  Carolyn Begeman
!> \date    2021
!> \details
!>  This routine determines the layerThickness to remap to subject to
!>  constraints
!
!-----------------------------------------------------------------------

   subroutine ocn_vert_regrid(restingThickness, layerThicknessLag, &
                              layerThicknessTarget, err)

      real (kind=RKIND), dimension(:, :), intent(in) :: &
         layerThicknessLag, & ! layerThickness after the lagrangian step
         restingThickness
      real (kind=RKIND), dimension(:, :), intent(out) :: &
         layerThicknessTarget ! adjusted target layerThickness for remapping
      integer, intent(out) :: err !< Output: Error flag

      integer :: iCell, k, kmin, kmax

      real (kind=RKIND) :: totalThickness

      err = 0

      ! Calculate z-star layer locations
      layerThicknessTarget = 0.0_RKIND

      !$omp parallel
      !$omp do schedule(runtime) private(k,totalThickness)
      do iCell = 1, nCellsAll
         totalThickness = 0.0_RKIND
         do k = maxLevelCell(iCell), minLevelCell(iCell), -1
            totalThickness = totalThickness + layerThicknessLag(k,iCell)
         end do
         do k = maxLevelCell(iCell), minLevelCell(iCell), -1
            layerThicknessTarget(k, iCell) = restingThickness(k, iCell) * &
                                             totalThickness / bottomDepth(iCell)
         end do
      end do
      !$omp end do
      !$omp end parallel


   end subroutine ocn_vert_regrid


!***********************************************************************
!
!  routine ocn_vert_regrid_init
!
!> \brief   Initializes ocean vertical regridding
!> \author  Carolyn Begeman
!> \date    July 2021
!> \details
!>  This routine initializes parameters required for vertical Lagrangian
!>  regridding
!
!-----------------------------------------------------------------------

   subroutine ocn_vert_regrid_init(err)

      integer, intent(out) :: err !< Output: Error flag

      err = 0

   end subroutine ocn_vert_regrid_init

end module ocn_vertical_regrid
! vim: foldmethod=marker
