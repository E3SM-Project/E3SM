! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vvel_horiz_advection
!
!> \brief MPAS horizontal advection for (non-hydrostatic) vertical velocity
!> \author Luke Van Roekel, Darren Engwirda, Sara Calandrini
!> \date   January 2021, updated Novermber 2021
!> \details
!>  This module contains routines for horizontal advection of vertical velocity
!>  via a simple TRiSK-type scheme at vertical levels (ie. half-layers).
!
!-------------------------------------------------------------------------------

module ocn_vvel_horiz_advection

   ! module includes
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_io_units
   use mpas_threading
   use ocn_config
   use ocn_mesh
   use ocn_tracer_advection_shared
   use ocn_diagnostics_variables

   implicit none
   private
   save 

   ! public method interfaces
   public :: ocn_vvel_horiz_advection_tend, &
             ocn_vvel_horiz_advection_init

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_vvel_horiz_advection_tend
!
!> \brief MPAS horizontal advection tendency for (non-hydrostatic) vertical velocity
!> \author Luke Van Roekel, Darren Engwirda, Sara Calandrini
!> \date   12/21/20, updated 11/19/2021
!> \details
!>  NOTE this is only active for non-hydrostatic mode
!
!-----------------------------------------------------------------------

   subroutine ocn_vvel_horiz_advection_tend(vertVelocityNonhydro, &
      normalThicknessFlux, layerThickness, dt, meshPool, tend)!{{{ 

      implicit none

      !*** input parameters

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertVelocityNonhydro                   !< [in] non-hydrostatic vertical velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalThicknessFlux                !< [in] normal thickness flux uh at edges
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness                     !< [in] layer thickness at cells
      real (kind=RKIND), intent(in) :: dt   !< [in] time-step
      type (mpas_pool_type), intent(in) :: meshPool !< [in] mesh information

      !*** output parameters
      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend                                       !< [out] vvel tendencies

      integer :: iCell, iEdge, jEdge, k, cell1, cell2
      integer :: nVertLevels, nCells, nEdges
      real (kind=RKIND) :: vertVelocityFlux
 
      real (kind=RKIND), dimension(:), allocatable :: invVolCell
      real (kind=RKIND), dimension(:,:), allocatable :: horizFlux

      ! Get dimensions
      nVertLevels = size(vertVelocityNonhydro,dim=1)

      ! Initialize pointers
      nCells = nCellsHalo( 1 ) !1
      nEdges = nEdgesHalo( 2 ) !2
      allocate(invVolCell(nVertLevels))
      allocate(horizFlux(nVertLevels, nEdges)) 

      !$omp parallel
      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
         horizFlux(:, iEdge) = 0.0_RKIND
      end do
      !$omp end do
      !$omp end parallel

      ! compute edge fluxes for div(u*h*w)
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(k, cell1, cell2, vertVelocityFlux)
      do iEdge = 1, nEdges 
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = 1, maxLevelEdgeTop(iEdge)
            if (k .gt. 1) then
               vertVelocityFlux = &
                  0.5_RKIND * (normalThicknessFlux(k-1,iEdge) + &
                               normalThicknessFlux(k-0,iEdge)) * &
                  0.5_RKIND * (vertVelocityNonhydro(k,cell1) + vertVelocityNonhydro(k,cell2))
            else
               vertVelocityFlux = & ! surface
                  1.0_RKIND * (normalThicknessFlux(k-0,iEdge)) * &
                  0.5_RKIND * (vertVelocityNonhydro(k,cell1) + vertVelocityNonhydro(k,cell2))
            end if

            horizFlux(k,iEdge) = &
               edgeMask(k,iEdge) * dvEdge(iEdge) * vertVelocityFlux
         end do
      end do
      !$omp end do
      !$omp end parallel
     
      ! eval. tendencies as div(u*h*w) / (A*h)
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(k, iEdge, jEdge, invVolCell)
      do iCell = 1, nCells
         invVolCell(:) =-1.0e34_RKIND
         invVolCell(1) = 1.0_RKIND / (areaCell(iCell) * layerThickness(1,iCell))
         do k = 2, maxLevelCell(iCell)
            invVolCell(k) = 2.0_RKIND / ( &
               areaCell(iCell) * (layerThickness(k-1,iCell) + layerThickness(k-0,iCell)))
         end do

         do jEdge = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(jEdge,iCell)
            do k = 1, maxLevelEdgeTop(iEdge)
               tend(k,iCell) = tend(k,iCell) + &
                  edgeSignOnCell(jEdge,iCell) * invVolCell(k) * horizFlux(k,iEdge)  
            end do
         end do
      end do
      !$omp end do
      !$omp end parallel

      deallocate(horizFlux)
      deallocate(invVolCell) 

   end subroutine ocn_vvel_horiz_advection_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_vvel_horiz_advection_init
!
!> \brief MPAS initialize vertical velocity advection tendency.
!> \author Luke Van Roekel, Darren Engwirda, Sara Calandrini
!> \date   12/21/2020, updated 11/19/2021
!> \details
!>  This routine initializes constants and choices for the (non-hydrostatic) 
!>  vertical velocity advection tendency
!
!-------------------------------------------------------------------------------

   subroutine ocn_vvel_horiz_advection_init(err) !{{{ 

      implicit none

      !*** input parameters 

      !*** output parameters

      integer, intent(out) :: err !< [out] Error Flag

      err = 0 ! set error code to success
     
      ! don't have actual choices for vvel advection here, unlike tracers       

   end subroutine ocn_vvel_horiz_advection_init!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

end module ocn_vvel_horiz_advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
