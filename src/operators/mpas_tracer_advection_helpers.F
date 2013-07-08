! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_tracer_advection_helpers
!
!> \brief MPAS tracer advection helper functions
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This module contains helper functions tracer advection.
!
!-----------------------------------------------------------------------
module mpas_tracer_advection_helpers

   use mpas_kind_types
   use mpas_grid_types
   use mpas_configure

   implicit none
   save

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  function mpas_tracer_advection_vflux4
!
!> \brief MPAS 4th order vertical tracer advection stencil
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This function provides the stencil for 4th order vertical advection of tracers.
!
!-----------------------------------------------------------------------
   real (kind=RKIND) function mpas_tracer_advection_vflux4(q_im2, q_im1, q_i, q_ip1, w)!{{{
        real (kind=RKIND), intent(in) :: q_im2 !< Input: Tracer value at index i-2
        real (kind=RKIND), intent(in) :: q_im1 !< Input: Tracer value at index i-1
        real (kind=RKIND), intent(in) :: q_i !< Input: Tracer value at index i
        real (kind=RKIND), intent(in) :: q_ip1 !< Input: Tracer value at index i+1
        real (kind=RKIND), intent(in) :: w !< Input: vertical veloicity
        mpas_tracer_advection_vflux4 = w*( 7.0_RKIND*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0_RKIND
   end function!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  function mpas_tracer_advection_vflux3
!
!> \brief MPAS 3rd order vertical tracer advection stencil
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This function provides the stencil for 3rd order vertical advection of tracers.
!
!-----------------------------------------------------------------------
   real (kind=RKIND) function mpas_tracer_advection_vflux3( q_im2, q_im1, q_i, q_ip1, w, coef)!{{{
        real (kind=RKIND), intent(in) :: q_im2 !< Input: Tracer value at index i-2
        real (kind=RKIND), intent(in) :: q_im1 !< Input: Tracer value at index i-1
        real (kind=RKIND), intent(in) :: q_i !< Input: Tracer value at index i
        real (kind=RKIND), intent(in) :: q_ip1 !< Input: Tracer value at index i+1
        real (kind=RKIND), intent(in) :: w !< Input: vertical veloicity
        real (kind=RKIND), intent(in) :: coef !< Input: Advection coefficient

        mpas_tracer_advection_vflux3 = (w * (7.0_RKIND * (q_i + q_im1) - (q_ip1 + q_im2)) - coef * abs(w) * ((q_ip1 - q_im2) - 3.0_RKIND*(q_i-q_im1)))/12.0_RKIND
   end function!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_tracer_advection_coefficients
!
!> \brief MPAS tracer advection coefficients
!> \author Doug Jacobsen, Bill Skamarock
!> \date   03/09/12
!> \details
!>  This routine precomputes the advection coefficients for horizontal
!>  advection of tracers.
!
!-----------------------------------------------------------------------
   subroutine mpas_tracer_advection_coefficients( grid, err, maxLevelCell_in, highOrderAdvectionMask_in, boundaryCell_in )!{{{

      use mpas_hash

      implicit none
      type (mesh_type) :: grid !< Input: Grid information
      integer, intent(out) :: err !< Input/Output: Error flag
      integer, dimension(:), pointer, optional :: maxLevelCell_in !< Input - optional: Index to last real cell
      integer, dimension(:,:), pointer, optional :: highOrderAdvectionMask_in !< Input - optional: Mask for high order advection
      integer, dimension(:,:), pointer, optional :: boundaryCell_in !< Input- optional: Mask for boundary cells

      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_2nd, adv_coefs_3rd
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, advCellsForEdge, highOrderAdvectionMask, lowOrderAdvectionMask, boundaryCell
      integer, dimension(:), pointer :: nEdgesOnCell, nAdvCellsForEdge, maxLevelCell

      integer, dimension(:), pointer :: cell_indices 
      integer, dimension(:,:), pointer :: sorted_cell_indices
      integer :: cell1, cell2, iEdge, n, i, j, j_in, iCell, k, nVertLevels
      logical :: addcell, highOrderAdvection

      type (hashtable) :: cell_hash

      deriv_two => grid % deriv_two % array
      adv_coefs => grid % adv_coefs % array
      adv_coefs_3rd => grid % adv_coefs_3rd % array
      cellsOnCell => grid % cellsOnCell % array
      cellsOnEdge => grid % cellsOnEdge % array
      advCellsForEdge => grid % advCellsForEdge % array
      nEdgesOnCell => grid % nEdgesOnCell % array
      nAdvCellsForEdge => grid % nAdvCellsForEdge % array

      nVertLevels = grid % nVertLevels

      allocate(cell_indices(grid % maxEdges2 + 2))
      allocate(sorted_cell_indices(2, grid % maxEdges2 + 2))

      err = 0

      if(present(maxLevelCell_in)) then
        maxLevelCell => maxLevelCell_in
      else
        allocate(maxLevelCell(grid % nCells+1))
        maxLevelCell(:) = nVertLevels
      end if

      if(present(highOrderAdvectionMask_in)) then
        highOrderAdvectionMask => highOrderAdvectionMask_in
        highOrderAdvectionMask = 0
      end if

      if(present(boundaryCell_in)) then
        boundaryCell => boundaryCell_in
      else
        allocate(boundaryCell(nVertLevels, grid % nCells+1))
        boundaryCell(:,:) = 0
      end if

      do iEdge = 1, grid % nEdges
        nAdvCellsForEdge(iEdge) = 0
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)

        if(present(highOrderAdvectionMask_in)) then
          do k = 1, nVertLevels
            if (boundaryCell(k, cell1) == 1 .or. boundaryCell(k, cell2) == 1) then
              highOrderAdvectionMask(k, iEdge) = 0
            else
              highOrderAdvectionMask(k, iEdge) = 1
            end if
          end do
        end if

        !
        ! do only if this edge flux is needed to update owned cells
        !
        if (cell1 <= grid % nCells .and. cell2 <= grid % nCells) then
           ! Insert cellsOnEdge to list of advection cells
           call mpas_hash_init(cell_hash)
           call mpas_hash_insert(cell_hash, cell1)
           call mpas_hash_insert(cell_hash, cell2)
           cell_indices(1) = cell1
           cell_indices(2) = cell2
           sorted_cell_indices(1, 1) = grid % indexToCellID % array(cell1)
           sorted_cell_indices(2, 1) = cell1
           sorted_cell_indices(1, 2) = grid % indexToCellID % array(cell2)
           sorted_cell_indices(2, 2) = cell2
           n = 2

           ! Build unique list of cells used for advection on edge
           do i = 1, nEdgesOnCell(cell1)
             if(.not. mpas_hash_search(cell_hash, cellsOnCell(i, cell1))) then
               n = n + 1
               cell_indices(n) = cellsOnCell(i, cell1)
               sorted_cell_indices(1, n) = grid % indexToCellID % array(cellsOnCell(i, cell1))
               sorted_cell_indices(2, n) = cellsOnCell(i, cell1)
               call mpas_hash_insert(cell_hash, cellsOnCell(i, cell1))
             end if
           end do ! loop over i

           do i = 1, nEdgesOnCell(cell2)
             if(.not. mpas_hash_search(cell_hash, cellsOnCell(i, cell2))) then
               n = n + 1
               cell_indices(n) = cellsOnCell(i, cell2)
               sorted_cell_indices(1, n) = grid % indexToCellID % array(cellsOnCell(i, cell2))
               sorted_cell_indices(2, n) = cellsOnCell(i, cell2)
               call mpas_hash_insert(cell_hash, cellsOnCell(i, cell2))
             end if
           end do ! loop over i

           call mpas_hash_destroy(cell_hash)

           call mpas_quicksort(n, sorted_cell_indices)

           nAdvCellsForEdge(iEdge) = n
           do iCell = 1, nAdvCellsForEdge(iEdge)
             advCellsForEdge(iCell, iEdge) = sorted_cell_indices(2, iCell)
           end do ! loop over iCell

           adv_coefs(:,iEdge) = 0.
           adv_coefs_3rd(:,iEdge) = 0.

           k = mpas_binary_search(sorted_cell_indices, 2, 1, nAdvCellsForEdge(iEdge), grid % indexToCellID % array(cell1))
           if(k <= nAdvCellsForEdge(iEdge)) then
             adv_coefs(k, iEdge) = adv_coefs(k, iEdge) + deriv_two(1,1,iEdge)
             adv_coefs_3rd(k, iEdge) = adv_coefs_3rd(k, iEdge) + deriv_two(1,1,iEdge)
           end if

           do iCell = 1, nEdgesOnCell(cell1)
             k = mpas_binary_search(sorted_cell_indices, 2, 1, nAdvCellsForEdge(iEdge), grid % indexToCellID % array(cellsOnCell(iCell,cell1)))
             if(k <= nAdvCellsForEdge(iEdge)) then
               adv_coefs(k, iEdge) = adv_coefs(k, iEdge) + deriv_two(iCell+1, 1, iEdge)
               adv_coefs_3rd(k, iEdge) = adv_coefs_3rd(k, iEdge) + deriv_two(iCell+1, 1, iEdge)
             end if
           end do ! loop over iCell

           k = mpas_binary_search(sorted_cell_indices, 2, 1, nAdvCellsForEdge(iEdge), grid % indexToCellID % array(cell2))
           if(k <= nAdvCellsForEdge(iEdge)) then
             adv_coefs(k, iEdge) = adv_coefs(k, iEdge) + deriv_two(1,2,iEdge)
             adv_coefs_3rd(k, iEdge) = adv_coefs_3rd(k, iEdge) + deriv_two(1,2,iEdge)
           end if

           do iCell = 1, nEdgesOnCell(cell2)
             k = mpas_binary_search(sorted_cell_indices, 2, 1, nAdvCellsForEdge(iEdge), grid % indexToCellID % array(cellsOnCell(iCell,cell2)))
             if(k <= nAdvCellsForEdge(iEdge)) then
               adv_coefs(k, iEdge) = adv_coefs(k, iEdge) + deriv_two(iCell+1, 2, iEdge)
               adv_coefs_3rd(k, iEdge) = adv_coefs_3rd(k, iEdge) + deriv_two(iCell+1, 2, iEdge)
             end if
           end do ! loop over iCell

           do iCell = 1,nAdvCellsForEdge(iEdge)
             adv_coefs    (iCell,iEdge) = - (grid % dcEdge % array (iEdge) **2) * adv_coefs    (iCell,iEdge) / 12.
             adv_coefs_3rd(iCell,iEdge) = - (grid % dcEdge % array (iEdge) **2) * adv_coefs_3rd(iCell,iEdge) / 12.
           end do ! loop over iCell

           k = mpas_binary_search(sorted_cell_indices, 2, 1, nAdvCellsForEdge(iEdge), grid % indexToCellID % array(cell1))
           if(k <= nAdvCellsForEdge(iEdge)) then
             adv_coefs(k, iEdge) = adv_coefs(k, iEdge) + 0.5
             adv_coefs_2nd(k, iEdge) = adv_coefs_2nd(k, iEdge) + 0.5
           end if

           k = mpas_binary_search(sorted_cell_indices, 2, 1, nAdvCellsForEdge(iEdge), grid % indexToCellID % array(cell2))
           if(k <= nAdvCellsForEdge(iEdge)) then
             adv_coefs(k, iEdge) = adv_coefs(k, iEdge) + 0.5
             adv_coefs_2nd(k, iEdge) = adv_coefs_2nd(k, iEdge) + 0.5
           end if

           do iCell=1,nAdvCellsForEdge(iEdge)
             adv_coefs    (iCell,iEdge) = grid % dvEdge % array(iEdge) * adv_coefs    (iCell,iEdge)
             adv_coefs_2nd(iCell,iEdge) = grid % dvEdge % array(iEdge) * adv_coefs_2nd(iCell,iEdge)
             adv_coefs_3rd(iCell,iEdge) = grid % dvEdge % array(iEdge) * adv_coefs_3rd(iCell,iEdge)
           end do ! loop over iCell
        end if
      end do ! end loop over edges

      deallocate(cell_indices)
      deallocate(sorted_cell_indices)

      ! If 2nd order advection, set masks appropriately.
      if(config_horiz_tracer_adv_order == 2 .and. present(highOrderAdvectionMask_in)) then
        highOrderAdvectionMask = 0
      end if

      if(.not.present(maxLevelCell_in)) then
        deallocate(maxLevelCell)
      end if

      if(.not.present(boundaryCell_in)) then
        deallocate(boundaryCell)
      end if

   end subroutine mpas_tracer_advection_coefficients!}}}

end module mpas_tracer_advection_helpers
