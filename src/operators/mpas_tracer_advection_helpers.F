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
   use mpas_sort
   use mpas_geometry_utils

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
   subroutine mpas_tracer_advection_coefficients( grid, horiz_adv_order, deriv_two, adv_coefs, adv_coefs_3rd, nAdvCellsForEdge, advCellsForEdge, err, maxLevelCell_in, highOrderAdvectionMask_in, boundaryCell_in )!{{{

      use mpas_hash

      implicit none
      type (mesh_type) :: grid !< Input: Grid information
      integer, intent(in) :: horiz_adv_order !< Input: Order of horizontal advection
      real (kind=RKIND), dimension(:,:,:), intent(in) :: deriv_two !< Input: 2nd derivative values for polynomial fit to tracers
      real (kind=RKIND), dimension(:,:), intent(inout) :: adv_coefs !< Input/Output: Advection coefficients for 2nd order advection
      real (kind=RKIND), dimension(:,:), intent(inout) :: adv_coefs_3rd !< Input/Output: Advection coefficients for blending in 3rd or 4th order advection
      integer, dimension(:), intent(inout) :: nAdvCellsForEdge !< Input/Output: Number of advection cells for each edge
      integer, dimension(:,:), intent(inout) :: advCellsForEdge !< Input/Output: List of advection cells for each edge
      integer, intent(out) :: err !< Input/Output: Error flag
      integer, dimension(:), pointer, optional :: maxLevelCell_in !< Input - optional: Index to last real cell
      integer, dimension(:,:), pointer, optional :: highOrderAdvectionMask_in !< Input - optional: Mask for high order advection
      integer, dimension(:,:), pointer, optional :: boundaryCell_in !< Input- optional: Mask for boundary cells

      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, highOrderAdvectionMask, boundaryCell
      integer, dimension(:), pointer :: nEdgesOnCell, maxLevelCell

      integer, dimension(:), pointer :: cell_indices 
      integer, dimension(:,:), pointer :: sorted_cell_indices
      integer :: cell1, cell2, iEdge, n, i, j, j_in, iCell, k, nVertLevels
      logical :: addcell, highOrderAdvection

      type (hashtable) :: cell_hash

      cellsOnCell => grid % cellsOnCell % array
      cellsOnEdge => grid % cellsOnEdge % array
      nEdgesOnCell => grid % nEdgesOnCell % array

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
           end if

           k = mpas_binary_search(sorted_cell_indices, 2, 1, nAdvCellsForEdge(iEdge), grid % indexToCellID % array(cell2))
           if(k <= nAdvCellsForEdge(iEdge)) then
             adv_coefs(k, iEdge) = adv_coefs(k, iEdge) + 0.5
           end if

           do iCell=1,nAdvCellsForEdge(iEdge)
             adv_coefs    (iCell,iEdge) = grid % dvEdge % array(iEdge) * adv_coefs    (iCell,iEdge)
             adv_coefs_3rd(iCell,iEdge) = grid % dvEdge % array(iEdge) * adv_coefs_3rd(iCell,iEdge)
           end do ! loop over iCell
        end if
      end do ! end loop over edges

      deallocate(cell_indices)
      deallocate(sorted_cell_indices)

      ! If 2nd order advection, set masks appropriately.
      if(horiz_adv_order == 2 .and. present(highOrderAdvectionMask_in)) then
        highOrderAdvectionMask = 0
      end if

      if(.not.present(maxLevelCell_in)) then
        deallocate(maxLevelCell)
      end if

      if(.not.present(boundaryCell_in)) then
        deallocate(boundaryCell)
      end if

   end subroutine mpas_tracer_advection_coefficients!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_initialize_deriv_two
!
!> \brief MPAS deriv two computation
!> \author Doug Jacobsen, Bill Skamarock
!> \date   03/09/12
!> \details
!>  This routine precomputes the second derivative values for tracer advection.
!
!-----------------------------------------------------------------------
   subroutine mpas_initialize_deriv_two( grid, deriv_two, err )!{{{
                                      
!
! compute the cell coefficients for the polynomial fit.
! this is performed during setup for model integration.
! WCS, 31 August 2009
!
! Described in:
! Skamarock, W. C., & Gassmann, A. (2011). 
! Conservative Transport Schemes for Spherical Geodesic Grids: High-Order Flux Operators for ODE-Based Time Integration. 
! Monthly Weather Review, 139(9), 2962-2975. doi:10.1175/MWR-D-10-05056.1
!
      implicit none

      type (mesh_type), intent(in) :: grid
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: deriv_two !< Input/Output: 2nd derivative values of polynomial for tracer fit.
      integer, intent(out) :: err

      integer, dimension(:), pointer :: advCells

!  local variables

      real (kind=RKIND), dimension(2, grid % nEdges) :: thetae
      real (kind=RKIND), dimension(grid % nEdges) :: xe, ye
      real (kind=RKIND), dimension(grid % nCells) :: theta_abs

      real (kind=RKIND), dimension(25) :: xc, yc, zc ! cell center coordinates
      real (kind=RKIND), dimension(25) :: thetav, thetat, dl_sphere
      real (kind=RKIND) :: xm, ym, zm, dl, xec, yec, zec
      real (kind=RKIND) :: thetae_tmp, xe_tmp, ye_tmp
      real (kind=RKIND) :: xv1, xv2, yv1, yv2, zv1, zv2
      integer :: i, j, k, ip1, ip2, m, n, ip1a, ii
      integer :: iCell, iEdge
      real (kind=RKIND) :: pii
      real (kind=RKIND) :: x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5
      real (kind=RKIND) :: pdx1, pdx2, pdx3, pdy1, pdy2, pdy3, dx1, dx2, dy1, dy2
      real (kind=RKIND) :: angv1, angv2, dl1, dl2
      real (kind=RKIND), dimension(25) :: dxe, dye, x2v, y2v, xp, yp
      
      real (kind=RKIND) :: amatrix(25,25), bmatrix(25,25), wmatrix(25,25)
      real (kind=RKIND) :: length_scale
      integer :: ma,na, cell_add, mw, nn
      integer, dimension(25) :: cell_list


      integer :: cell1, cell2
      integer, parameter :: polynomial_order = 2
      logical, parameter :: debug = .false.
      logical, parameter :: least_squares = .true.
      logical :: add_the_cell, do_the_cell

      logical, parameter :: reset_poly = .true.

      real (kind=RKIND) :: rcell, cos2t, costsint, sin2t
      real (kind=RKIND), dimension(grid%maxEdges) :: angle_2d

!---    
      err = 0

      if(polynomial_order > 2) then
        write (*,*) 'Polynomial for second derivitave can only be 2'
        err = 1
        return
      end if

      pii = 2.*asin(1.0)

!     advCells => grid % advCells % array
      allocate(advCells(grid % maxEdges2))
      deriv_two(:,:,:) = 0.

      do iCell = 1, grid % nCells !  is this correct? - we need first halo cell also...

         cell_list(1) = iCell
         do i=2, grid % nEdgesOnCell % array(iCell)+1
            cell_list(i) = grid % CellsOnCell % array(i-1,iCell)
         end do
         n = grid % nEdgesOnCell % array(iCell) + 1

         if ( polynomial_order > 2 ) then
            do i=2,grid % nEdgesOnCell % array(iCell) + 1
               do j=1,grid % nEdgesOnCell % array ( cell_list(i) )
                  cell_add = grid % CellsOnCell % array (j,cell_list(i))
                  add_the_cell = .true.
                  do k=1,n
                     if ( cell_add == cell_list(k) ) add_the_cell = .false.
                  end do
                  if (add_the_cell) then
                     n = n+1
                     cell_list(n) = cell_add
                  end if
               end do
            end do
         end if
 
         advCells(1) = n

!  check to see if we are reaching outside the halo

         do_the_cell = .true.
         do i=1,n
            if (cell_list(i) > grid % nCells) do_the_cell = .false.
         end do


         if ( .not. do_the_cell ) cycle


!  compute poynomial fit for this cell if all needed neighbors exist
         if ( grid % on_a_sphere ) then

            do i=1,n
               advCells(i+1) = cell_list(i)
               xc(i) = grid % xCell % array(advCells(i+1))/grid % sphere_radius
               yc(i) = grid % yCell % array(advCells(i+1))/grid % sphere_radius
               zc(i) = grid % zCell % array(advCells(i+1))/grid % sphere_radius
            end do

            if ( zc(1) == 1.0_RKIND) then
               theta_abs(iCell) = pii/2.0_RKIND
            else
               theta_abs(iCell) =  pii/2.0_RKIND - mpas_sphere_angle( xc(1), yc(1), zc(1),  &
                                                                 xc(2), yc(2), zc(2),  &
                                                                 0.0_RKIND, 0.0_RKIND, 1.0_RKIND ) 
            end if

! angles from cell center to neighbor centers (thetav)

            do i=1,n-1
   
               ip2 = i+2
               if (ip2 > n) ip2 = 2
    
               thetav(i) = mpas_sphere_angle( xc(1),   yc(1),   zc(1),    &
                                         xc(i+1), yc(i+1), zc(i+1),  &
                                         xc(ip2), yc(ip2), zc(ip2)   )

               dl_sphere(i) = grid % sphere_radius*mpas_arc_length( xc(1),   yc(1),   zc(1),  &
                                            xc(i+1), yc(i+1), zc(i+1) )
            end do

            length_scale = 1.
            do i=1,n-1
               dl_sphere(i) = dl_sphere(i)/length_scale
            end do

!            thetat(1) = 0.  !  this defines the x direction, cell center 1 -> 
            thetat(1) = theta_abs(iCell)  !  this defines the x direction, longitude line
            do i=2,n-1
               thetat(i) = thetat(i-1) + thetav(i-1)
            end do
   
            do i=1,n-1
               xp(i) = cos(thetat(i)) * dl_sphere(i)
               yp(i) = sin(thetat(i)) * dl_sphere(i)
            end do

         else     ! On an x-y plane

            do i=1,n-1

               angle_2d(i) = grid%angleEdge%array(grid % EdgesOnCell % array(i,iCell))
               iEdge = grid % EdgesOnCell % array(i,iCell)
               if ( iCell .ne. grid % CellsOnEdge % array(1,iEdge)) &
                  angle_2d(i) = angle_2d(i) - pii

!              xp(i) = grid % xCell % array(cell_list(i)) - grid % xCell % array(iCell)
!              yp(i) = grid % yCell % array(cell_list(i)) - grid % yCell % array(iCell)

               xp(i) = grid % dcEdge % array(grid % EdgesOnCell % array(i,iCell)) * cos(angle_2d(i))
               yp(i) = grid % dcEdge % array(grid % EdgesOnCell % array(i,iCell)) * sin(angle_2d(i))

            end do

         end if


         ma = n-1
         mw = grid % nEdgesOnCell % array (iCell)

         bmatrix = 0.
         amatrix = 0.
         wmatrix = 0.

         if (polynomial_order == 2) then
            na = 6
            ma = ma+1
  
            amatrix(1,1) = 1.
            wmatrix(1,1) = 1.
            do i=2,ma
               amatrix(i,1) = 1.
               amatrix(i,2) = xp(i-1)
               amatrix(i,3) = yp(i-1)
               amatrix(i,4) = xp(i-1)**2
               amatrix(i,5) = xp(i-1) * yp(i-1)
               amatrix(i,6) = yp(i-1)**2
   
               wmatrix(i,i) = 1.
            end do
 
         else if (polynomial_order == 3) then
            na = 10
            ma = ma+1
  
            amatrix(1,1) = 1.
            wmatrix(1,1) = 1.
            do i=2,ma
               amatrix(i,1) = 1.
               amatrix(i,2) = xp(i-1)
               amatrix(i,3) = yp(i-1)
   
               amatrix(i,4) = xp(i-1)**2
               amatrix(i,5) = xp(i-1) * yp(i-1)
               amatrix(i,6) = yp(i-1)**2
   
               amatrix(i,7) = xp(i-1)**3
               amatrix(i,8) = yp(i-1) * (xp(i-1)**2)
               amatrix(i,9) = xp(i-1) * (yp(i-1)**2)
               amatrix(i,10) = yp(i-1)**3
   
               wmatrix(i,i) = 1.
 
            end do

         else
            na = 15
            ma = ma+1
  
            amatrix(1,1) = 1.
            wmatrix(1,1) = 1.
            do i=2,ma
               amatrix(i,1) = 1.
               amatrix(i,2) = xp(i-1)
               amatrix(i,3) = yp(i-1)
   
               amatrix(i,4) = xp(i-1)**2
               amatrix(i,5) = xp(i-1) * yp(i-1)
               amatrix(i,6) = yp(i-1)**2
   
               amatrix(i,7) = xp(i-1)**3
               amatrix(i,8) = yp(i-1) * (xp(i-1)**2)
               amatrix(i,9) = xp(i-1) * (yp(i-1)**2)
               amatrix(i,10) = yp(i-1)**3
   
               amatrix(i,11) = xp(i-1)**4
               amatrix(i,12) = yp(i-1) * (xp(i-1)**3)
               amatrix(i,13) = (xp(i-1)**2)*(yp(i-1)**2)
               amatrix(i,14) = xp(i-1) * (yp(i-1)**3)
               amatrix(i,15) = yp(i-1)**4
   
               wmatrix(i,i) = 1.
  
            end do
 
            do i=1,mw
               wmatrix(i,i) = 1.
            end do
 
         end if
 
         call mpas_poly_fit_2( amatrix, bmatrix, wmatrix, ma, na, 25 )

         do i=1,grid % nEdgesOnCell % array (iCell)
            ip1 = i+1
            if (ip1 > n-1) ip1 = 1
  
            iEdge = grid % EdgesOnCell % array (i,iCell)

            if(grid % on_a_sphere) then
              xv1 = grid % xVertex % array(grid % verticesOnEdge % array (1,iedge))/grid % sphere_radius
              yv1 = grid % yVertex % array(grid % verticesOnEdge % array (1,iedge))/grid % sphere_radius
              zv1 = grid % zVertex % array(grid % verticesOnEdge % array (1,iedge))/grid % sphere_radius
              xv2 = grid % xVertex % array(grid % verticesOnEdge % array (2,iedge))/grid % sphere_radius
              yv2 = grid % yVertex % array(grid % verticesOnEdge % array (2,iedge))/grid % sphere_radius
              zv2 = grid % zVertex % array(grid % verticesOnEdge % array (2,iedge))/grid % sphere_radius
            else
              xv1 = grid % xVertex % array(grid % verticesOnEdge % array (1,iedge))
              yv1 = grid % yVertex % array(grid % verticesOnEdge % array (1,iedge))
              zv1 = grid % zVertex % array(grid % verticesOnEdge % array (1,iedge))
              xv2 = grid % xVertex % array(grid % verticesOnEdge % array (2,iedge))
              yv2 = grid % yVertex % array(grid % verticesOnEdge % array (2,iedge))
              zv2 = grid % zVertex % array(grid % verticesOnEdge % array (2,iedge))
            end if
  
            if ( grid % on_a_sphere ) then
               call mpas_arc_bisect( xv1, yv1, zv1,  &
                                xv2, yv2, zv2,  &
                                xec, yec, zec   )
  
               thetae_tmp = mpas_sphere_angle( xc(1),   yc(1),   zc(1),    &
                                          xc(i+1), yc(i+1), zc(i+1),  &
                                          xec,     yec,     zec       )
               thetae_tmp = thetae_tmp + thetat(i)
               if (iCell == grid % cellsOnEdge % array(1,iEdge)) then
                  thetae(1,grid % EdgesOnCell % array (i,iCell)) = thetae_tmp
               else
                  thetae(2,grid % EdgesOnCell % array (i,iCell)) = thetae_tmp
               end if
!            else
!
!               xe(grid % EdgesOnCell % array (i,iCell)) = 0.5 * (xv1 + xv2)
!               ye(grid % EdgesOnCell % array (i,iCell)) = 0.5 * (yv1 + yv2)

            end if
  
         end do

!  fill second derivative stencil for rk advection 

         do i=1, grid % nEdgesOnCell % array (iCell)
            iEdge = grid % EdgesOnCell % array (i,iCell)
  
  
            if ( grid % on_a_sphere ) then
               if (iCell == grid % cellsOnEdge % array(1,iEdge)) then
  
                  cos2t = cos(thetae(1,grid % EdgesOnCell % array (i,iCell)))
                  sin2t = sin(thetae(1,grid % EdgesOnCell % array (i,iCell)))
                  costsint = cos2t*sin2t
                  cos2t = cos2t**2
                  sin2t = sin2t**2
   
                  do j=1,n
                     deriv_two(j,1,iEdge) =   2.*cos2t*bmatrix(4,j)  &
                                            + 2.*costsint*bmatrix(5,j)  &
                                            + 2.*sin2t*bmatrix(6,j)
                  end do
               else
     
                  cos2t = cos(thetae(2,grid % EdgesOnCell % array (i,iCell)))
                  sin2t = sin(thetae(2,grid % EdgesOnCell % array (i,iCell)))
                  costsint = cos2t*sin2t
                  cos2t = cos2t**2
                  sin2t = sin2t**2
      
                  do j=1,n
                     deriv_two(j,2,iEdge) =   2.*cos2t*bmatrix(4,j)  &
                                            + 2.*costsint*bmatrix(5,j)  &
                                            + 2.*sin2t*bmatrix(6,j)
                  end do
               end if

            else

               cos2t = cos(angle_2d(i))
               sin2t = sin(angle_2d(i))
               costsint = cos2t*sin2t
               cos2t = cos2t**2
               sin2t = sin2t**2

!               do j=1,n
!
!                  deriv_two(j,1,iEdge) =   2.*xe(iEdge)*xe(iEdge)*bmatrix(4,j)  &
!                                         + 2.*xe(iEdge)*ye(iEdge)*bmatrix(5,j)  &
!                                         + 2.*ye(iEdge)*ye(iEdge)*bmatrix(6,j)
!               end do

               if (iCell == grid % cellsOnEdge % array(1,iEdge)) then
                  do j=1,n
                     deriv_two(j,1,iEdge) =   2.*cos2t*bmatrix(4,j)  &
                                            + 2.*costsint*bmatrix(5,j)  &
                                            + 2.*sin2t*bmatrix(6,j)
                  end do
               else
                  do j=1,n
                     deriv_two(j,2,iEdge) =   2.*cos2t*bmatrix(4,j)  &
                                            + 2.*costsint*bmatrix(5,j)  &
                                            + 2.*sin2t*bmatrix(6,j)
                  end do
               end if

            end if
         end do
 
      end do ! end of loop over cells

      if (debug) stop


!      write(0,*) ' check for deriv2 coefficients, iEdge 4 '
!
!      iEdge = 4
!      j = 1
!      iCell = grid % cellsOnEdge % array(1,iEdge)
!      write(0,*) ' j, icell, coef ',j,iCell,deriv_two(j,1,iEdge)
!      do j=2,7
!         write(0,*) ' j, icell, coef ',j,grid % CellsOnCell % array(j-1,iCell),deriv_two(j,1,iEdge)
!      end do
!
!      j = 1
!      iCell = grid % cellsOnEdge % array(2,iEdge)
!      write(0,*) ' j, icell, coef ',j,iCell,deriv_two(j,2,iEdge)
!      do j=2,7
!         write(0,*) ' j, icell, coef ',j,grid % CellsOnCell % array(j-1,iCell),deriv_two(j,2,iEdge)
!      end do
!      stop

   end subroutine mpas_initialize_deriv_two!}}}

end module mpas_tracer_advection_helpers
