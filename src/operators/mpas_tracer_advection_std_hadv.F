! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_tracer_advection_std_hadv
!
!> \brief MPAS standard horizontal tracer advection (non-monotonic/fct)
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This module contains routines for horizontal tracer advection tendencies
!>  and initialization of the horizontal advection routines.
!
!-----------------------------------------------------------------------
module mpas_tracer_advection_std_hadv

   use mpas_kind_types
   use mpas_grid_types
   use mpas_configure
   use mpas_dmpar

   use mpas_tracer_advection_helpers

   implicit none
   private
   save

   public :: mpas_tracer_advection_std_hadv_tend, &
             mpas_tracer_advection_std_hadv_init

   real (kind=RKIND) :: coef_3rd_order

   logical :: hadvOn

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_tracer_advection_std_hadv_tend
!
!> \brief MPAS standard horizontal tracer advection tendency
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine computes the tendency for 3rd order horizontal advection of tracers.
!
!-----------------------------------------------------------------------
   subroutine mpas_tracer_advection_std_hadv_tend(tracers, normalThicknessFlux, grid, tend)!{{{
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tend !< Input/output: Tracer tendency
      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers !< Input: Tracer values
      real (kind=RKIND), dimension(:,:), intent(in) :: normalThicknessFlux !< Input: Thickness weighted horizontal velocity
      type (mesh_type), intent(in) :: grid !< Input: Grid information

      integer :: i, iCell, iEdge, k, iTracer, cell1, cell2
      real (kind=RKIND) :: flux, tracer_weight

      real (kind=RKIND), dimension(:), pointer :: areaCell
      integer, dimension(:,:), pointer :: cellsOnEdge

      integer, dimension(:,:), pointer :: advCellsForEdge, highOrderAdvectionMask, lowOrderAdvectionMask
      integer, dimension(:), pointer :: nAdvCellsForEdge
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_2nd, adv_coefs_3rd
      real (kind=RKIND), dimension(:,:), allocatable :: flux_arr
      integer :: nVertLevels, num_tracers

      if (.not. hadvOn) return

      cellsOnEdge => grid % cellsOnEdge % array
      areaCell    => grid % areaCell % array

      nAdvCellsForEdge => grid % nAdvCellsForEdge % array
      advCellsForEdge => grid % advCellsForEdge % array
      adv_coefs => grid % adv_coefs % array
      adv_coefs_2nd => grid % adv_coefs_2nd % array
      adv_coefs_3rd => grid % adv_coefs_3rd % array
      highOrderAdvectionMask => grid % highOrderAdvectionMask % array
      lowOrderAdvectionMask => grid % lowOrderAdvectionMask % array

      nVertLevels = grid % nVertLevels
      num_tracers = size(tracers, dim=1)

      allocate(flux_arr(num_tracers, grid % nVertLevels))

      !  horizontal flux divergence, accumulate in tracer_tend
      do iEdge=1,grid%nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         if (cell1 <= grid%nCellsSolve .or. cell2 <= grid%nCellsSolve) then  ! only for owned cells
            flux_arr(:,:) = 0.
            do i=1,nAdvCellsForEdge(iEdge)
              iCell = advCellsForEdge(i,iEdge)
              do k=1,grid % nVertLevels
              tracer_weight = lowOrderAdvectionMask(k, iEdge) * adv_coefs_2nd(i,iEdge) & 
                            + highOrderAdvectionMask(k, iEdge) * (adv_coefs(i,iEdge) + coef_3rd_order*sign(1.0_RKIND,normalThicknessFlux(k,iEdge))*adv_coefs_3rd(i,iEdge))
                do iTracer=1,num_tracers
                  flux_arr(iTracer,k) = flux_arr(iTracer,k) + tracer_weight* tracers(iTracer,k,iCell)
                end do
              end do
            end do

            do k=1,grid % nVertLevels
               do iTracer=1,num_tracers
                  tend(iTracer,k,cell1) = tend(iTracer,k,cell1) - normalThicknessFlux(k,iEdge)*flux_arr(iTracer,k)/areaCell(cell1)
                  tend(iTracer,k,cell2) = tend(iTracer,k,cell2) + normalThicknessFlux(k,iEdge)*flux_arr(iTracer,k)/areaCell(cell2)
               end do
            end do
         end if
       end do

       deallocate(flux_arr)

   end subroutine mpas_tracer_advection_std_hadv_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_tracer_advection_std_hadv_init
!
!> \brief MPAS standard horizontal tracer advection initialization
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine initializes the 3rd order standard horizontal advection of tracers
!
!-----------------------------------------------------------------------
   subroutine mpas_tracer_advection_std_hadv_init(err)!{{{
      integer, intent(inout) :: err !< Input/Output: Error flag

      err = 0

      hadvOn =.true.

      if ( config_horiz_tracer_adv_order == 3) then
        coef_3rd_order = config_coef_3rd_order
      else if ( config_horiz_tracer_adv_order == 2 .or. config_horiz_tracer_adv_order == 4) then
        coef_3rd_order = 0.0
      end if
   end subroutine mpas_tracer_advection_std_hadv_init!}}}

end module mpas_tracer_advection_std_hadv
