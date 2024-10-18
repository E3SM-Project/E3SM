! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_ocn_tracer_advection_shared
!
!> \brief MPAS ocean tracer advection common coefficients and variables
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This module contains the routines and arrays for some common 
!>  tracer advection quantities.
!
!-----------------------------------------------------------------------

module ocn_tracer_advection_shared

   use mpas_kind_types
   use mpas_derived_types
   use mpas_hash
   use mpas_sort
   use mpas_constants
   use mpas_geometry_utils
   use mpas_log

   use ocn_config
   use ocn_constants
   use ocn_mesh

   implicit none
   private
   save

   !--------------------------------------------------------------------
   ! Public variables (should only be used by other advection modules)
   !--------------------------------------------------------------------

   integer, public :: &
      nAdvCellsMax     ! largest number of advection cells for any edge

   integer, public, dimension(:), allocatable :: &
      nAdvCellsForEdge ! number of cells contrib to advection at edge

   integer, public, dimension(:,:), allocatable :: &
      advCellsForEdge ! index of cells contributing to advection at edge

   real (kind=RKIND), public, dimension(:,:), allocatable :: &
      advCoefs,       &! common advection coefficients
      advCoefs3rd,    &! common advection coeffs for high order
      advMaskHighOrder ! mask for high order advection terms

   !--------------------------------------------------------------------
   ! Public member functions 
   !-------------------------------------------------------------------- 

   public :: ocn_tracer_advection_shared_init

!***********************************************************************

   contains

!***********************************************************************
!
!  routine ocn_tracer_advection_shared_init
!
!> \brief MPAS tracer advection coefficients
!> \author Doug Jacobsen, Bill Skamarock
!> \date   03/09/12
!> \details
!>  This routine precomputes advection coefficients for horizontal
!>  advection of tracers.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_advection_shared_init(err)!{{{

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] Error flag

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer :: &
         n, i, k,      &! loop indices for neighbor, vertical loops
         iEdge, iCell, &! loop indices for edge, cell loops
         cell1, cell2   ! neighbor cell indices across edge

      integer, dimension(:), allocatable :: &
         cellIndx       ! cell indices gathered from neighbors 

      integer, dimension(:,:), allocatable :: &
         cellIndxSorted ! nbr cell indices sorted

      real (kind=RKIND), dimension(:,:,:), allocatable :: &
         derivTwo !< 2nd derivative values for polynomial fit to tracers

      type (hashtable) :: cell_hash

      ! End preamble
      !-------------
      ! Begin code

      err = 0 ! initialize error code

      ! define max number of advection cells as max number of edges*2
      nAdvCellsMax = maxEdges2

      ! Allocate common variables
      allocate(nAdvCellsForEdge (             nEdgesAll),   &
                advCellsForEdge (nAdvCellsMax,nEdgesAll), &
                advCoefs        (nAdvCellsMax,nEdgesAll), &
                advCoefs3rd     (nAdvCellsMax,nEdgesAll), &
                advMaskHighOrder(nVertLevels ,nEdgesAll), &
                derivTwo      (nAdvCellsMax,2,nEdgesAll))

      ! Compute derivTwo array
      call computeDerivTwo(derivTwo, err)
      if (err /= 0) then
         call mpas_log_write( &
            'Error computing derivTwo in ocn advect shared init ', &
             MPAS_LOG_CRIT)
         deallocate (derivTwo)
         return
      endif

      allocate(cellIndx      (   maxEdges2 + 2), &
               cellIndxSorted(2, maxEdges2 + 2))

      ! Compute masks and coefficients
      do iEdge = 1, nEdgesAll
         nAdvCellsForEdge(iEdge) = 0
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         ! at boundaries, must stay at low order
         do k = 1, nVertLevels
            if (boundaryCell(k,cell1) == 1 .or. &
                boundaryCell(k,cell2) == 1) then
               advMaskHighOrder(k,iEdge) = 0.0_RKIND
            else
               advMaskHighOrder(k,iEdge) = 1.0_RKIND
            end if
         end do

         ! do only if this edge flux is needed to update owned cells
         if (cell1 <= nCellsAll .and. cell2 <= nCellsAll) then
            ! Insert cellsOnEdge to list of advection cells
            call mpas_hash_init(cell_hash)
            call mpas_hash_insert(cell_hash, cell1)
            call mpas_hash_insert(cell_hash, cell2)
            cellIndx(1) = cell1
            cellIndx(2) = cell2
            cellIndxSorted(1,1) = indexToCellID(cell1)
            cellIndxSorted(2,1) = cell1
            cellIndxSorted(1,2) = indexToCellID(cell2)
            cellIndxSorted(2,2) = cell2
            n = 2

            ! Build unique list of cells used for advection on edge
            ! by expanding to the extended neighbor cells
            do i = 1, nEdgesOnCell(cell1)
               if (.not. mpas_hash_search(cell_hash, &
                                          cellsOnCell(i,cell1))) then
                  n = n + 1
                  cellIndx(n) = cellsOnCell(i, cell1)
                  cellIndxSorted(1,n) = indexToCellID(cellsOnCell(i,cell1))
                  cellIndxSorted(2,n) = cellsOnCell(i,cell1)
                  call mpas_hash_insert(cell_hash,cellsOnCell(i,cell1))
               end if
            end do

            do i = 1, nEdgesOnCell(cell2)
               if(.not. mpas_hash_search(cell_hash, &
                                         cellsOnCell(i, cell2))) then
                  n = n + 1
                  cellIndx(n) = cellsOnCell(i, cell2)
                  cellIndxSorted(1,n) = indexToCellID(cellsOnCell(i,cell2))
                  cellIndxSorted(2,n) = cellsOnCell(i,cell2)
                  call mpas_hash_insert(cell_hash, cellsOnCell(i,cell2))
               end if
            end do

            call mpas_hash_destroy(cell_hash)

            ! sort the cell indices by cellID
            call mpas_quicksort(n, cellIndxSorted)

            ! store local cell indices for high-order calculations
            nAdvCellsForEdge(iEdge) = n
            do iCell = 1, nAdvCellsForEdge(iEdge)
               advCellsForEdge(iCell,iEdge) = cellIndxSorted(2,iCell)
            end do

            ! equation 7 in Skamarock, W. C., & Gassmann, A. (2011):
            ! F(u,psi)_{i+1/2} = u_{i+1/2} *
            !  [   1/2 (psi_{i+1} + psi_i)                       term 1
            !    - 1/12(dx^2psi_{i+1} + dx^2psi_i)               term 2
            !    + sign(u) beta/12 (dx^2psi_{i+1} - dx^2psi_i)]  term 3
            !                                         (note minus sign)
            !
            ! advCoefs accounts for terms 1 and 2 in SG11 equation 7.
            ! Term 1 is the 2nd-order flux-function term. advCoefs
            ! accounts for this with the "+ 0.5" lines below. In the
            ! advection routines that use these coefficients, the
            ! 2nd-order flux loop is then skipped. Term 2 is the
            ! 4th-order flux-function term. advCoefs accounts for
            ! term 3, the beta term. beta > 0 corresponds to the
            ! third-order flux function. The - sign in the derivTwo
            ! accumulation is for the i+1 part of term 3, while
            ! the + sign is for the i part.

            do i=1,nAdvCellsMax
               advCoefs   (i,iEdge) = 0.0_RKIND
               advCoefs3rd(i,iEdge) = 0.0_RKIND
            end do

            ! pull together third and fourth order contributions to the
            ! flux first from cell1
            i = mpas_binary_search(cellIndxSorted, 2, 1, &
                    nAdvCellsForEdge(iEdge), indexToCellID(cell1))
            if (i <= nAdvCellsForEdge(iEdge)) then
               advCoefs   (i,iEdge) = advCoefs   (i,iEdge) &
                                    +  derivTwo(1,1,iEdge)
               advCoefs3rd(i,iEdge) = advCoefs3rd(i,iEdge) &
                                    +  derivTwo(1,1,iEdge)
            end if

            do iCell = 1, nEdgesOnCell(cell1)
               i = mpas_binary_search(cellIndxSorted, 2, 1, &
                              nAdvCellsForEdge(iEdge), &
                              indexToCellID(cellsOnCell(iCell,cell1)))
               if (i <= nAdvCellsForEdge(iEdge)) then
                  advCoefs   (i,iEdge) = advCoefs   (i,iEdge) &
                                       + derivTwo(iCell+1,1,iEdge)
                  advCoefs3rd(i,iEdge) = advCoefs3rd(i,iEdge) &
                                       + derivTwo(iCell+1,1,iEdge)
               end if
            end do

            ! pull together third and fourth order contributions to the
            ! flux now from cell2
            i = mpas_binary_search(cellIndxSorted, 2, 1, &
                       nAdvCellsForEdge(iEdge), indexToCellID(cell2))
            if (i <= nAdvCellsForEdge(iEdge)) then
               advCoefs   (i,iEdge) = advCoefs   (i,iEdge) &
                                    +  derivTwo(1,2,iEdge)
               advCoefs3rd(i,iEdge) = advCoefs3rd(i,iEdge) &
                                    -  derivTwo(1,2,iEdge)
            end if

            do iCell = 1, nEdgesOnCell(cell2)
               i = mpas_binary_search(cellIndxSorted, 2, 1, &
                           nAdvCellsForEdge(iEdge), &
                           indexToCellID(cellsOnCell(iCell,cell2)))
               if (i <= nAdvCellsForEdge(iEdge)) then
                  advCoefs   (i,iEdge) = advCoefs   (i,iEdge) &
                                       + derivTwo(iCell+1,2,iEdge)
                  advCoefs3rd(i,iEdge) = advCoefs3rd(i,iEdge) &
                                       - derivTwo(iCell+1,2,iEdge)
               end if
            end do

            do iCell = 1,nAdvCellsForEdge(iEdge)
               advCoefs   (iCell,iEdge) = - (dcEdge(iEdge)**2)* &
                                     advCoefs(iCell,iEdge) / 12.
               advCoefs3rd(iCell,iEdge) = - (dcEdge(iEdge)**2)* &
                                  advCoefs3rd(iCell,iEdge) / 12.
            end do

            ! 2nd order centered contribution
            ! place this in the main flux weights
            i = mpas_binary_search(cellIndxSorted, 2, 1, &
                                   nAdvCellsForEdge(iEdge), &
                                   indexToCellID(cell1))
            if (i <= nAdvCellsForEdge(iEdge)) then
               advCoefs(i,iEdge) = advCoefs(i, iEdge) + 0.5
            end if

            i = mpas_binary_search(cellIndxSorted, 2, 1, &
                                   nAdvCellsForEdge(iEdge), &
                                   indexToCellID(cell2))
            if (i <= nAdvCellsForEdge(iEdge)) then
               advCoefs(i,iEdge) = advCoefs(i, iEdge) + 0.5
            end if

            ! multiply by edge length - thus the flux is just dt*ru
            ! times the results of the vector-vector multiply
            do iCell=1,nAdvCellsForEdge(iEdge)
               advCoefs   (iCell,iEdge) = dvEdge(iEdge)* &
                                          advCoefs   (iCell,iEdge)
               advCoefs3rd(iCell,iEdge) = dvEdge(iEdge)* &
                                          advCoefs3rd(iCell,iEdge)
            end do
         end if  ! only do for edges of owned-cells
      end do ! end loop over edges

      deallocate(cellIndx, &
                 cellIndxSorted, &
                 derivTwo)

      ! If 2nd order advection, disable high-order terms by
      ! setting mask to zero.
      if (config_horiz_tracer_adv_order == 2) &
         advMaskHighOrder(:,:) = 0.0_RKIND

      ! Copy module variables to device
      !$acc enter data copyin(nAdvCellsForEdge, advCellsForEdge, &
      !$acc                   advCoefs, advCoefs3rd, advMaskHighOrder)

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_advection_shared_init !}}}

!***********************************************************************
!
!  routine computeDerivTwo
!
!> \brief MPAS deriv two computation
!> \author Doug Jacobsen, Bill Skamarock
!> \date   03/09/12
!> \details
!>  This routine precomputes the second derivative values for tracer
!>  advection. It computes cell coefficients for the polynomial fit
!>  as described in:
!>  Skamarock, W. C., & Gassmann, A. (2011). 
!>    Conservative Transport Schemes for Spherical Geodesic Meshs:
!>    High-Order Flux Operators for ODE-Based Time Integration. 
!>    Monthly Weather Review, 139(9), 2962-2975.
!>    doi:10.1175/MWR-D-10-05056.1
!>  This is performed during model initialization.
!
!-----------------------------------------------------------------------

   subroutine computeDerivTwo(derivTwo, err)!{{{
                                      
      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(out) :: &
         derivTwo !< [out] 2nd deriv values of polynomial for tracer fit

      integer, intent(out) :: err !< [out] returned error flag

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      logical :: &
         addCell,      &! flag for adding cell add to list 
         doCell         ! flag for whether to process cell

      integer, parameter :: &
         polynomial_order = 2  ! option for polynomial order, forced to 2

      integer :: &
         i, j, n,   &! loop counters
         ip1, ip2,     &! temps for i+1, i+2
         iCell, iEdge, &! loop indices for cell, edge loops
         ma, na, cell_add, mw

      real (kind=RKIND) :: &
         xec, yec, zec,    &! arc bisection coords
         thetae_tmp,       &! angle
         xv1, xv2, yv1, yv2, zv1, zv2, &! vertex cart coords
         length_scale,     &! length scale
         cos2t, costsint, sin2t  ! trig function temps

      real (kind=RKIND), dimension(:,:), allocatable :: &
         thetae          ! sphere angle between vectors 

      real (kind=RKIND), dimension(:), allocatable :: &
         theta_abs,     &! sphere angle
         angle_2d        ! 2d angle

      integer, dimension(25) :: cell_list

      real (kind=RKIND), dimension(25) :: &
         xc, yc, zc,    &! cell center coordinates
         xp, yp,        &!
         thetav, thetat, dl_sphere

      real (kind=RKIND) :: &
         amatrix(25,25), bmatrix(25,25), wmatrix(25,25)

      ! End preamble
      !-------------
      ! Begin code

      ! set return code and check proper poly order
      err = 0

      if (polynomial_order > 2) then
        call mpas_log_write( &
           'Polynomial for second derivitave can only be 2',  &
           MPAS_LOG_ERR)
        err = 1
        return
      end if

      ! allocate local arrays
      allocate(thetae(2, nEdgesAll))
      allocate(theta_abs(nCellsAll))
      allocate(angle_2d(maxEdges))

      ! Initialize derivTwo and pi

      derivTwo(:,:,:) = 0.0_RKIND

      do iCell = 1, nCellsAll

         cell_list(1) = iCell
         do i=2, nEdgesOnCell(iCell)+1
            cell_list(i) = cellsOnCell(i-1,iCell)
         end do
         n = nEdgesOnCell(iCell) + 1

         !if ( polynomial_order > 2 ) then
         !   do i=2, nEdgesOnCell(iCell) + 1
         !      do j=1, nEdgesOnCell( cell_list(i) )
         !         cell_add = cellsOnCell(j,cell_list(i))
         !         addCell = .true.
         !         do k=1,n
         !            if ( cell_add == cell_list(k) ) addCell = .false.
         !         end do
         !         if (addCell) then
         !            n = n+1
         !            cell_list(n) = cell_add
         !         end if
         !      end do
         !   end do
         !end if
 
         ! check to see if we are reaching outside the halo

         doCell = .true.
         do i=1,n
            if (cell_list(i) > nCellsAll) doCell = .false.
         end do

         if ( .not. doCell ) cycle

         ! compute poynomial fit for this cell if all needed
         ! neighbors exist
         if ( onSphere ) then

            do i=1,n
               j = cell_list(i)
               xc(i) = xCell(j) / sphereRadius
               yc(i) = yCell(j) / sphereRadius
               zc(i) = zCell(j) / sphereRadius
            end do

            if ( zc(1) == 1.0_RKIND) then
               theta_abs(iCell) = pi/2.0_RKIND
            else
               theta_abs(iCell) = pi/2.0_RKIND &
                        - mpas_sphere_angle( xc(1), yc(1), zc(1), &
                                             xc(2), yc(2), zc(2), &
                                   0.0_RKIND, 0.0_RKIND, 1.0_RKIND) 
            end if

            ! angles from cell center to neighbor centers (thetav)

            do i=1,n-1
   
               ip2 = i+2
               if (ip2 > n) ip2 = 2
    
               thetav(i) = mpas_sphere_angle(xc(  1), yc(  1), zc(  1),&
                                             xc(i+1), yc(i+1), zc(i+1),&
                                             xc(ip2), yc(ip2), zc(ip2))

               dl_sphere(i) = sphereRadius* &
                              mpas_arc_length(xc(  1),yc(  1),zc(  1), &
                                              xc(i+1),yc(i+1),zc(i+1))
            end do

            length_scale = 1.0_RKIND
            do i=1,n-1
               dl_sphere(i) = dl_sphere(i)/length_scale
            end do

            ! thetat(1) = 0.  !  this defines the x direction, cell center 1 -> 
            ! this defines the x direction, longitude line
            thetat(1) = theta_abs(iCell)
            do i=2,n-1
               thetat(i) = thetat(i-1) + thetav(i-1)
            end do
   
            do i=1,n-1
               xp(i) = cos(thetat(i)) * dl_sphere(i)
               yp(i) = sin(thetat(i)) * dl_sphere(i)
            end do

         else     ! On an x-y plane

            do i=1,n-1

               angle_2d(i) = angleEdge(edgesOnCell(i,iCell))
               iEdge = edgesOnCell(i,iCell)
               if ( iCell .ne. cellsOnEdge(1,iEdge)) &
                  angle_2d(i) = angle_2d(i) - pi

               xp(i) = dcEdge(edgesOnCell(i,iCell)) * cos(angle_2d(i))
               yp(i) = dcEdge(edgesOnCell(i,iCell)) * sin(angle_2d(i))

            end do

         end if

         ma = n-1
         mw = nEdgesOnCell(iCell)

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

         do i=1, nEdgesOnCell(iCell)
            ip1 = i+1
            if (ip1 > n-1) ip1 = 1
  
            iEdge = edgesOnCell(i,iCell)

            if ( onSphere ) then
              xv1 = xVertex(verticesOnEdge(1,iedge)) / sphereRadius
              yv1 = yVertex(verticesOnEdge(1,iedge)) / sphereRadius
              zv1 = zVertex(verticesOnEdge(1,iedge)) / sphereRadius
              xv2 = xVertex(verticesOnEdge(2,iedge)) / sphereRadius
              yv2 = yVertex(verticesOnEdge(2,iedge)) / sphereRadius
              zv2 = zVertex(verticesOnEdge(2,iedge)) / sphereRadius
            else
              xv1 = xVertex(verticesOnEdge(1,iedge))
              yv1 = yVertex(verticesOnEdge(1,iedge))
              zv1 = zVertex(verticesOnEdge(1,iedge))
              xv2 = xVertex(verticesOnEdge(2,iedge))
              yv2 = yVertex(verticesOnEdge(2,iedge))
              zv2 = zVertex(verticesOnEdge(2,iedge))
            end if
  
            if ( onSphere ) then
               call mpas_arc_bisect( xv1, yv1, zv1,  &
                                     xv2, yv2, zv2,  &
                                     xec, yec, zec   )
  
               thetae_tmp = mpas_sphere_angle( &
                                          xc(  1), yc(  1), zc(  1), &
                                          xc(i+1), yc(i+1), zc(i+1), &
                                          xec,     yec,     zec      )
               thetae_tmp = thetae_tmp + thetat(i)
               if (iCell == cellsOnEdge(1,iEdge)) then
                  thetae(1, edgesOnCell(i,iCell)) = thetae_tmp
               else
                  thetae(2, edgesOnCell(i,iCell)) = thetae_tmp
               end if
            end if
  
         end do

         ! fill second derivative stencil for rk advection 

         do i=1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
  
            if ( onSphere ) then
               if (iCell == cellsOnEdge(1,iEdge)) then
  
                  cos2t = cos(thetae(1, edgesOnCell(i,iCell)))
                  sin2t = sin(thetae(1, edgesOnCell(i,iCell)))
                  costsint = cos2t*sin2t
                  cos2t = cos2t**2
                  sin2t = sin2t**2
   
                  do j=1,n
                     derivTwo(j,1,iEdge) = 2.*cos2t   *bmatrix(4,j) &
                                         + 2.*costsint*bmatrix(5,j) &
                                         + 2.*sin2t   *bmatrix(6,j)
                  end do

               else
     
                  cos2t = cos(thetae(2, edgesOnCell(i,iCell)))
                  sin2t = sin(thetae(2, edgesOnCell(i,iCell)))
                  costsint = cos2t*sin2t
                  cos2t = cos2t**2
                  sin2t = sin2t**2
      
                  do j=1,n
                     derivTwo(j,2,iEdge) = 2.*cos2t   *bmatrix(4,j) &
                                         + 2.*costsint*bmatrix(5,j) &
                                         + 2.*sin2t   *bmatrix(6,j)
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
!                  derivTwo(j,1,iEdge) =   2.*xe(iEdge)*xe(iEdge)*bmatrix(4,j)  &
!                                         + 2.*xe(iEdge)*ye(iEdge)*bmatrix(5,j)  &
!                                         + 2.*ye(iEdge)*ye(iEdge)*bmatrix(6,j)
!               end do

               if (iCell == cellsOnEdge(1,iEdge)) then
                  do j=1,n
                     derivTwo(j,1,iEdge) = 2.*cos2t   *bmatrix(4,j) &
                                         + 2.*costsint*bmatrix(5,j) &
                                         + 2.*sin2t   *bmatrix(6,j)
                  end do
               else
                  do j=1,n
                     derivTwo(j,2,iEdge) = 2.*cos2t   *bmatrix(4,j) &
                                         + 2.*costsint*bmatrix(5,j) &
                                         + 2.*sin2t   *bmatrix(6,j)
                  end do
               end if

            end if
         end do
 
      end do ! end of loop over cells

      deallocate(thetae)
      deallocate(theta_abs)
      deallocate(angle_2d)

   !--------------------------------------------------------------------

   end subroutine computeDerivTwo !}}}

!***********************************************************************

end module ocn_tracer_advection_shared

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
