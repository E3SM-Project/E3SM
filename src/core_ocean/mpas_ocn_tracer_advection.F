! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_advection
!
!> \brief MPAS ocean tracer advection driver
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This module contains driver routine for tracer advection tendencys
!>  as well as the routines for setting up advection coefficients and 
!>  initialization of the advection routines.
!
!-----------------------------------------------------------------------

module mpas_ocn_tracer_advection

   use mpas_kind_types
   use mpas_grid_types
   use mpas_configure
   use mpas_sort
   use mpas_hash

   use mpas_ocn_tracer_advection_std
   use mpas_ocn_tracer_advection_mono
     
   implicit none
   private
   save

   public :: mpas_ocn_tracer_advection_init,         &
             mpas_ocn_tracer_advection_coefficients, &
             mpas_ocn_tracer_advection_tend

   logical :: tracerAdvOn
   logical :: monotonicOn

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_ocn_tracer_advection_coefficients
!
!> \brief MPAS ocean tracer advection coefficients
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine precomputes the advection coefficients for horizontal
!>  advection of tracers.
!
!-----------------------------------------------------------------------
   subroutine mpas_ocn_tracer_advection_coefficients( grid, err )!{{{

      implicit none
      type (mesh_type) :: grid !< Input: Grid information
      integer, intent(out) :: err

      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_2nd, adv_coefs_3rd
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, advCellsForEdge, highOrderAdvectionMask, lowOrderAdvectionMask, boundaryCell
      integer, dimension(:), pointer :: nEdgesOnCell, nAdvCellsForEdge, maxLevelCell

      integer, dimension(:), pointer :: cell_indices
      integer, dimension(:,:), pointer :: sorted_cell_indices
      integer :: cell1, cell2, iEdge, n, i, j, j_in, iCell, k, nVertLevels, nCells
      logical :: addcell, highOrderAdvection

      type (hashtable) :: cell_hash

      deriv_two => grid % deriv_two % array
      adv_coefs => grid % adv_coefs % array
      adv_coefs_2nd => grid % adv_coefs_2nd % array
      adv_coefs_3rd => grid % adv_coefs_3rd % array
      cellsOnCell => grid % cellsOnCell % array
      cellsOnEdge => grid % cellsOnEdge % array
      advCellsForEdge => grid % advCellsForEdge % array
      boundaryCell => grid % boundaryCell % array
      highOrderAdvectionMask => grid % highOrderAdvectionMask % array
      lowOrderAdvectionMask => grid % lowOrderAdvectionMask % array
      nEdgesOnCell => grid % nEdgesOnCell % array
      maxLevelCell => grid % maxLevelCell % array
      nAdvCellsForEdge => grid % nAdvCellsForEdge % array

      nCells = grid % nCells
      nVertLevels = grid % nVertLevels

      allocate(cell_indices(grid % maxEdges2 + 2))
      allocate(sorted_cell_indices(2, grid % maxEdges2 + 2))

      err = 0

      highOrderAdvectionMask = 0
      lowOrderAdvectionMask = 0

      do iEdge = 1, grid % nEdges
        nAdvCellsForEdge(iEdge) = 0
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)

        do k = 1, nVertLevels
          if (boundaryCell(k, cell1) == 1 .or. boundaryCell(k, cell2) == 1) then
            highOrderAdvectionMask(k, iEdge) = 0
            lowOrderAdvectionMask(k, iEdge) = 1
          else
            highOrderAdvectionMask(k, iEdge) = 1
            lowOrderAdvectionMask(k, iEdge) = 0
          end if
        end do

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
           adv_coefs_2nd(:,iEdge) = 0.
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
      if(config_horiz_tracer_adv_order == 2) then
        lowOrderAdvectionMask = 1
        highOrderAdvectionMask = 0
      end if

      if (maxval(highOrderAdvectionMask+lowOrderAdvectionMask) > 1) then
        write(*,*) "Masks don't sum to 1."
        err = 1
      endif

   end subroutine mpas_ocn_tracer_advection_coefficients!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_ocn_tracer_advection_tend
!
!> \brief MPAS ocean tracer advection tendency
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine is the driver routine for computing the tendency for
!>  advection of tracers.
!
!-----------------------------------------------------------------------
   subroutine mpas_ocn_tracer_advection_tend(tracers, normalThicknessFlux, w, layerThickness, verticalCellSize, dt, grid, tend_layerThickness, tend)!{{{

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tend !< Input/Output: tracer tendency
      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers !< Input/Output: tracer values
      real (kind=RKIND), dimension(:,:), intent(in) :: normalThicknessFlux !< Input: Thickness weighted horizontal velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: w  !< Input: Vertical velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness  !< Input: Thickness field
      real (kind=RKIND), dimension(:,:), intent(in) :: verticalCellSize !< Input: Distance between vertical interfaces of a cell
      real (kind=RKIND), intent(in) :: dt !< Input: Time step
      type (mesh_type), intent(in) :: grid !< Input: grid information
      real (kind=RKIND), dimension(:,:), intent(in) :: tend_layerThickness !< Input: Thickness tendency information

      if(.not. tracerAdvOn) return

      if(monotonicOn) then
         call mpas_ocn_tracer_advection_mono_tend(tracers, normalThicknessFlux, w, layerThickness, verticalCellSize, dt, grid, tend_layerThickness, tend)
      else
         call mpas_ocn_tracer_advection_std_tend(tracers, normalThicknessFlux, w, verticalCellSize, grid, tend)
      endif
   end subroutine mpas_ocn_tracer_advection_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_ocn_tracer_advection_init
!
!> \brief MPAS ocean tracer advection tendency
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine is the driver routine for initialization of 
!>  the tracer advection routines.
!
!-----------------------------------------------------------------------
   subroutine mpas_ocn_tracer_advection_init(err)!{{{

      integer, intent(inout) :: err !< Input/Output: Error flag

      integer :: err_tmp

      err = 0

      tracerAdvOn = .true.

      if(config_disable_tr_adv) tracerAdvOn = .false.

      call mpas_ocn_tracer_advection_std_init(err_tmp)
      call mpas_ocn_tracer_advection_mono_init(err_tmp)

      err = ior(err, err_tmp)

      monotonicOn = .false.

      if(config_monotonic) then
         monotonicOn = .true.
      endif

   end subroutine mpas_ocn_tracer_advection_init!}}}

end module mpas_ocn_tracer_advection
