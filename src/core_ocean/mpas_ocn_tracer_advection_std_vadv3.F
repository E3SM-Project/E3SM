! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_advection_std_vadv3
!
!> \brief MPAS ocean 3rd order vertical tracer advection driver (non-monotonic/fct)
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This module contains routines for 3rd order vertical tracer advection tendencies.
!
!-----------------------------------------------------------------------
module mpas_ocn_tracer_advection_std_vadv3

   use mpas_kind_types
   use mpas_grid_types
   use mpas_configure
   use mpas_dmpar

   use mpas_ocn_tracer_advection_helpers

   implicit none
   private
   save

   public :: mpas_ocn_tracer_advection_std_vadv3_tend

   real (kind=RKIND) :: coef_3rd_order

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_ocn_tracer_advection_std_vadv3_tend
!
!> \brief MPAS ocean 3rd order standard vertical tracer advection tendency
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine computes the 3rd order vertical tracer advection tendencies.
!
!-----------------------------------------------------------------------
   subroutine mpas_ocn_tracer_advection_std_vadv3_tend(tracers, w, verticalCellSize, grid, tend)!{{{
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tend !< Input/Output: Tracer Tendency
      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers !< Input: Tracer Values
      real (kind=RKIND), dimension(:,:), intent(in) :: w !< Input: Vertical Velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: verticalCellSize !< Input: Distance between vertical interfaces of cell
      type (mesh_type), intent(in) :: grid !< Input: Grid information

      integer :: i, iCell, iEdge, k, iTracer, cell1, cell2

      real (kind=RKIND), dimension(:,:), allocatable :: vert_flux
      real (kind=RKIND) :: weightK, weightKm1
      integer :: nVertLevels, num_tracers
      integer, dimension(:), pointer :: maxLevelCell

      coef_3rd_order = config_coef_3rd_order

      nVertLevels = grid % nVertLevels
      num_tracers = size(tracers, dim=1)
      maxLevelCell => grid % maxLevelCell % array

      allocate(vert_flux(num_tracers, nVertLevels+1))

      vert_flux(:,1) = 0.

      do iCell=1,grid % nCellsSolve

        k = 2
        do iTracer=1,num_tracers
          weightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k-1, iCell) + verticalCellSize(k, iCell))
          weightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k-1, iCell) + verticalCellSize(k, iCell))
          vert_flux(iTracer,k) = w(k,iCell)*(weightK*tracers(iTracer,k,iCell)+weightKm1*tracers(iTracer,k-1,iCell))
        enddo
        
        do k=3,maxLevelCell(iCell)-1
          do iTracer=1,num_tracers
            vert_flux(iTracer,k) = mpas_ocn_tracer_advection_vflux3( tracers(iTracer,k-2,iCell),tracers(iTracer,k-1,iCell),  &
                                     tracers(iTracer,k  ,iCell),tracers(iTracer,k+1,iCell),  &
                                     w(k,iCell), coef_3rd_order )
          end do
        end do

        k = maxLevelCell(iCell)

        do iTracer=1,num_tracers
          weightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k-1, iCell) + verticalCellSize(k, iCell))
          weightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k-1, iCell) + verticalCellSize(k, iCell))
          vert_flux(iTracer,k) = w(k,iCell)*(weightK*tracers(iTracer,k,iCell)+weightKm1*tracers(iTracer,k-1,iCell))
        enddo

        vert_Flux(:, maxLevelCell(iCell)+1) = 0.0

        do k=1,maxLevelCell(iCell)
           do iTracer=1,num_tracers
             tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + (vert_flux(iTracer, k+1) - vert_flux(iTracer, k))
           end do
        end do
      end do

      deallocate(vert_flux)

   end subroutine mpas_ocn_tracer_advection_std_vadv3_tend!}}}

end module mpas_ocn_tracer_advection_std_vadv3
