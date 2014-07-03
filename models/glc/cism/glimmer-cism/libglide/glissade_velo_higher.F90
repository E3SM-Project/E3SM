!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_velo_higher.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!
! This module contains routines for computing the ice velocity
!  using a variational finite-element approach.
!
! See these papers for details:
!
! J.K. Dukowicz, S.F. Price and W.H. Lipscomb, 2010: Consistent
!    approximations and boundary conditions for ice-sheet dynamics
!    using a principle of least action.  J. Glaciology, 56 (197),
!    480-495.
!
! M. Perego, M. Gunzburger, and J. Burkardt, 2012: Parallel
!    finite-element implementation for higher-order ice-sheet models.
!    J. Glaciology, 58 (207), 76-88.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

  module glissade_velo_higher

    use glimmer_global, only: dp
    use glimmer_physcon, only: gn, rhoi, rhoo, grav, scyr
    use glimmer_paramets, only: thk0, len0, tim0, tau0, vel0   !TODO - remove these later
    use glimmer_paramets, only: vel_scale, len_scale   ! used for whichefvs = HO_EFVS_FLOWFACT
    use glimmer_log, only: write_log
    use glimmer_sparse_type
    use glimmer_sparse
     
    use glide_types  ! for HO_EFVS and other options

    use cism_sparse_pcg, only: pcg_solver_structured

    use parallel

    implicit none

    private
    public :: glissade_velo_higher_init, glissade_velo_higher_solve

    !----------------------------------------------------------------
    ! Here are some definitions:
    !
    ! The horizontal mesh is composed of cells and vertices.
    ! All cells are assumed to be quadrilaterals, but the code can be
    !  generalized later to triangles (e.g., for MPAS mesh).
    ! Each cell can be extruded to form a column with a specified number of layers.
    ! 
    ! An element is a layer of a cell, and a node is a corner of an element.
    ! Elements and nodes live in 3D space, whereas cells and vertices live in
    !  the horizontal plane.
    !
    ! Locally owned cells and vertices have indices (nhalo+1:nx-nhalo, nhalo+1,ny-nhalo).
    ! Active cells are cells that (1) contain ice and (2) border locally owned vertices.
    ! Active vertices are all vertices of active cells.
    ! Active nodes are all nodes in the columns associated with active vertices.
    !
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Finite element properties
    !
    ! Assume 3D hexahedral elements.
    !----------------------------------------------------------------

    integer, parameter ::        &
       nNodesPerElement = 8,     & ! 8 nodes for hexahedral elements
       nQuadPoints = 8,          & ! number of quadrature points per element
                                   ! These live at +- 1/sqrt(3) for reference hexahedron
       nNodesPerElement_2d = 4,  & ! 4 nodes for faces of hexahedral elements
       nQuadPoints_2d = 4          ! number of quadrature points per element face
                                   ! These live at +- 1/sqrt(3) for reference square

    real(dp), parameter ::     &
       rsqrt3 = 1.d0/sqrt(3.d0)    ! for quadrature points
         
    !----------------------------------------------------------------
    ! Arrays used for finite-element calculations
    !
    ! Most integals are done over 3D hexahedral elements.
    ! Surface integrals are done over 2D faces of these elements. 
    !----------------------------------------------------------------

    real(dp), dimension(nNodesPerElement, nQuadPoints) ::   & 
       phi,            &    ! trilinear basis function, evaluated at quad pts
       dphi_dxr,       &    ! dphi/dx for reference element, evaluated at quad pts
       dphi_dyr,       &    ! dphi/dy for reference element, evaluated at quad pts
       dphi_dzr             ! dphi/dy for reference element, evaluated at quad pts

    real(dp), dimension(nQuadPoints) ::  &
       xqp, yqp, zqp,  &    ! quad pt coordinates in reference element
       wqp                  ! quad pt weights

    real(dp), dimension(nNodesPerElement_2d, nQuadPoints_2d) ::   & 
       phi_2d,         &    ! bilinear basis function, evaluated at quad pts
       dphi_dxr_2d,    &    ! dphi/dx for reference square, evaluated at quad pts
       dphi_dyr_2d          ! dphi/dy for reference square, evaluated at quad pts

    real(dp), dimension(nQuadPoints_2d) ::  &
       xqp_2d, yqp_2d, &    ! quad pt coordinates in reference square
       wqp_2d               ! quad pt weights

    integer, dimension(nNodesPerElement, nNodesPerElement) ::  &
       ishift, jshift, kshift   ! matrices describing relative indices of nodes in an element

    integer, dimension(-1:1,-1:1,-1:1) :: &
       indxA                 ! maps relative (x,y,z) coordinates to an index between 1 and 27
                             ! index order is (i,j,k)

    real(dp), dimension(3,3) ::  &
       identity3                ! 3 x 3 identity matrix

    real(dp), parameter ::   &
       eps12 = 1.d-12           ! small number

!WHL - TODO - Keep volume scale?
    real(dp) :: vol0    ! volume scale = dx * dy * (1000 m)

    logical, parameter ::  &
       check_symmetry = .true.   ! if true, then check matrix symmetry

!WHL - debug
    logical :: verbose = .false.        ! for debug print statements
!    logical :: verbose = .true.  
    logical :: verbose_init = .false.   
    logical :: verbose_Jac = .false.

    integer, parameter :: &
       itest = 26, jtest = 19, ktest = 1, ptest = 1  ! for dome global (i,j) = (24,17), 1 proc
!       itest = 26, jtest = 19, ktest = 2, ptest = 1
!       itest = 26, jtest = 19, ktest = 10, ptest = 1
!       itest = 17, jtest = 10, ktest = 1, ptest = 1
!       itest = 10, jtest = 17, ktest = 1, ptest = 1

!        itest = 3, jtest = 3, ktest = 1, ptest = 1    ! for ishom.a.80km global (i,j) = (1,1)
!        itest = 22, jtest = 7, ktest = 1, ptest = 1    ! for ishom.a.80km symmetry check

    integer, parameter :: ntest = 2371  ! nodeID for dome global (24,17,1)    
!    integer, parameter :: ntest = 2372  ! nodeID for dome global (24,17,2)
!    integer, parameter :: ntest = 2380  ! nodeID for dome global (24,17,10)
!    integer, parameter :: ntest = 411  ! nodeID for dome global (17,10,1)
!    integer, parameter :: ntest = 1771  ! nodeID for dome global (10,17,1)
    
!    integer, parameter :: ntest = 1    ! for ishom.a.80km global (i,j) = (1,1)
!    integer, parameter :: ntest = 882    ! for ishom.a.80km symmetry check

    integer, parameter :: rtest = 0    ! rank for any single-process run
!    integer, parameter :: rtest = 0    ! rank for (17,10) with 4 procs

!WHL = debug
    integer, parameter :: rowtest = -999
    integer, parameter :: coltest = -999

    integer, parameter :: &
      iAtest=0, jAtest=0, kAtest=0

    contains

!****************************************************************************

  subroutine glissade_velo_higher_init

    !----------------------------------------------------------------
    ! Initial calculations for glissade higher-order solver.
    !----------------------------------------------------------------

    integer :: i, j, k, m, n, p

!WHL - debug
    real(dp) :: sumx, sumy, sumz

    !----------------------------------------------------------------
    ! Initialize some time-independent finite element arrays
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Trilinear basis set for reference hexahedron, x=(-1,1), y=(-1,1), z=(-1,1)             
    ! Indexing is counter-clockwise from SW corner, with 1-4 on lower surface
    !  and 5-8 on upper surface
    ! In the code we use "phi" to denote these basis functions. 
    !
    ! N1 = (1-x)*(1-y)*(1-z)/8             N4----N3
    ! N2 = (1+x)*(1-y)*(1-z)/8             |     |    Lower layer        
    ! N3 = (1+x)*(1+y)*(1-z)/8             |     |
    ! N4 = (1-x)*(1+y)*(1-z)/8             N1----N2

    ! N5 = (1-x)*(1-y)*(1+z)/8             N8----N7
    ! N6 = (1+x)*(1-y)*(1+z)/8             |     |    Upper layer
    ! N7 = (1+x)*(1+y)*(1+z)/8             |     |
    ! N8 = (1-x)*(1+y)*(1+z)/8             N5----N6
    !----------------------------------------------------------------
   
    ! Set coordinates and weights of quadrature points for reference hexahedral element.
    ! Numbering is counter-clockwise from southwest, lower face (1-4) followed by
    !  upper face (5-8).

    xqp(1) = -rsqrt3; yqp(1) = -rsqrt3; zqp(1) = -rsqrt3
    wqp(1) =  1.d0

    xqp(2) =  rsqrt3; yqp(2) = -rsqrt3; zqp(2) = -rsqrt3
    wqp(2) =  1.d0

    xqp(3) =  rsqrt3; yqp(3) =  rsqrt3; zqp(3) = -rsqrt3
    wqp(3) =  1.d0

    xqp(4) = -rsqrt3; yqp(4) =  rsqrt3; zqp(4) = -rsqrt3
    wqp(4) =  1.d0

    xqp(5) = -rsqrt3; yqp(5) = -rsqrt3; zqp(5) =  rsqrt3
    wqp(5) =  1.d0

    xqp(6) =  rsqrt3; yqp(6) = -rsqrt3; zqp(6) =  rsqrt3
    wqp(6) =  1.d0

    xqp(7) =  rsqrt3; yqp(7) =  rsqrt3; zqp(7) =  rsqrt3
    wqp(7) =  1.d0

    xqp(8) = -rsqrt3; yqp(8) =  rsqrt3; zqp(8) =  rsqrt3
    wqp(8) =  1.d0

    if (verbose_init) then
       print*, ' '
       sumx = 0.d0; sumy = 0.d0; sumz = 0.d0
       do p = 1, nQuadPoints
          print*, xqp(p), yqp(p), zqp(p)
          sumx = sumx + xqp(p); sumy = sumy + yqp(p); sumz = sumz + zqp(p)
       enddo
       print*, ' '
       print*, 'sums:', sumx, sumy, sumz
    endif

    ! Evaluate trilinear basis functions and their derivatives at each quad pt

    do p = 1, nQuadPoints

       phi(1,p) = (1.d0 - xqp(p)) * (1.d0 - yqp(p)) * (1.d0 - zqp(p)) / 8.d0
       phi(2,p) = (1.d0 + xqp(p)) * (1.d0 - yqp(p)) * (1.d0 - zqp(p)) / 8.d0
       phi(3,p) = (1.d0 + xqp(p)) * (1.d0 + yqp(p)) * (1.d0 - zqp(p)) / 8.d0
       phi(4,p) = (1.d0 - xqp(p)) * (1.d0 + yqp(p)) * (1.d0 - zqp(p)) / 8.d0
       phi(5,p) = (1.d0 - xqp(p)) * (1.d0 - yqp(p)) * (1.d0 + zqp(p)) / 8.d0
       phi(6,p) = (1.d0 + xqp(p)) * (1.d0 - yqp(p)) * (1.d0 + zqp(p)) / 8.d0
       phi(7,p) = (1.d0 + xqp(p)) * (1.d0 + yqp(p)) * (1.d0 + zqp(p)) / 8.d0
       phi(8,p) = (1.d0 - xqp(p)) * (1.d0 + yqp(p)) * (1.d0 + zqp(p)) / 8.d0

       dphi_dxr(1,p) = -(1.d0 - yqp(p)) * (1.d0 - zqp(p)) / 8.d0 
       dphi_dxr(2,p) =  (1.d0 - yqp(p)) * (1.d0 - zqp(p)) / 8.d0 
       dphi_dxr(3,p) =  (1.d0 + yqp(p)) * (1.d0 - zqp(p)) / 8.d0 
       dphi_dxr(4,p) = -(1.d0 + yqp(p)) * (1.d0 - zqp(p)) / 8.d0
       dphi_dxr(5,p) = -(1.d0 - yqp(p)) * (1.d0 + zqp(p)) / 8.d0 
       dphi_dxr(6,p) =  (1.d0 - yqp(p)) * (1.d0 + zqp(p)) / 8.d0 
       dphi_dxr(7,p) =  (1.d0 + yqp(p)) * (1.d0 + zqp(p)) / 8.d0 
       dphi_dxr(8,p) = -(1.d0 + yqp(p)) * (1.d0 + zqp(p)) / 8.d0

       dphi_dyr(1,p) = -(1.d0 - xqp(p)) * (1.d0 - zqp(p)) / 8.d0 
       dphi_dyr(2,p) = -(1.d0 + xqp(p)) * (1.d0 - zqp(p)) / 8.d0 
       dphi_dyr(3,p) =  (1.d0 + xqp(p)) * (1.d0 - zqp(p)) / 8.d0 
       dphi_dyr(4,p) =  (1.d0 - xqp(p)) * (1.d0 - zqp(p)) / 8.d0 
       dphi_dyr(5,p) = -(1.d0 - xqp(p)) * (1.d0 + zqp(p)) / 8.d0 
       dphi_dyr(6,p) = -(1.d0 + xqp(p)) * (1.d0 + zqp(p)) / 8.d0 
       dphi_dyr(7,p) =  (1.d0 + xqp(p)) * (1.d0 + zqp(p)) / 8.d0 
       dphi_dyr(8,p) =  (1.d0 - xqp(p)) * (1.d0 + zqp(p)) / 8.d0 

       dphi_dzr(1,p) = -(1.d0 - xqp(p)) * (1.d0 - yqp(p)) / 8.d0 
       dphi_dzr(2,p) = -(1.d0 + xqp(p)) * (1.d0 - yqp(p)) / 8.d0 
       dphi_dzr(3,p) = -(1.d0 + xqp(p)) * (1.d0 + yqp(p)) / 8.d0 
       dphi_dzr(4,p) = -(1.d0 - xqp(p)) * (1.d0 + yqp(p)) / 8.d0 
       dphi_dzr(5,p) =  (1.d0 - xqp(p)) * (1.d0 - yqp(p)) / 8.d0 
       dphi_dzr(6,p) =  (1.d0 + xqp(p)) * (1.d0 - yqp(p)) / 8.d0 
       dphi_dzr(7,p) =  (1.d0 + xqp(p)) * (1.d0 + yqp(p)) / 8.d0 
       dphi_dzr(8,p) =  (1.d0 - xqp(p)) * (1.d0 + yqp(p)) / 8.d0 

       if (verbose_init) then
          print*, ' '
          print*, 'Quad point, p =', p
          do n = 1, 8
             print*, n, phi(n,p)
!             print*, n, dphi_dxr(n,p)
!             print*, n, dphi_dyr(n,p)
!             print*, n, dphi_dzr(n,p)
          enddo
          print*, 'sum(phi)', sum(phi(:,p))  ! verified that sum = 1
          print*, 'sum(dphi/dx)', sum(dphi_dxr(:,p))  ! verified that sum = 0 (within roundoff)
          print*, 'sum(dphi/dy)', sum(dphi_dyr(:,p))  ! verified that sum = 0 (within roundoff)
          print*, 'sum(dphi/dz)', sum(dphi_dzr(:,p))  ! verified that sum = 0 (within roundoff)
       endif

    enddo   ! nQuadPoints

    ! Identity matrix
    identity3(1,:) = (/ 1.d0, 0.d0, 0.d0 /)
    identity3(2,:) = (/ 0.d0, 1.d0, 0.d0 /)
    identity3(3,:) = (/ 0.d0, 0.d0, 1.d0 /)

    ! Initialize some matrices that describe how the i, j and k indices of each node
    ! in each element are related to one another.

    ! The ishift matrix describes how the i indices of the 8 nodes are related to one another.
    ! E.g, if ishift (1,2) = 1, this means that node 2 has an i index
    ! one greater than the i index of node 1.

    ishift(1,:) = (/ 0,  1,  1,  0,  0,  1,  1,  0/)   
    ishift(2,:) = (/-1,  0,  0, -1, -1,  0,  0, -1/)   
    ishift(3,:) = ishift(2,:)
    ishift(4,:) = ishift(1,:)
    ishift(5,:) = ishift(1,:)
    ishift(6,:) = ishift(2,:)
    ishift(7,:) = ishift(2,:)
    ishift(8,:) = ishift(1,:)

    ! The jshift matrix describes how the j indices of the 8 nodes are related to one another.
    ! E.g, if jshift (1,4) = 1, this means that node 4 has a j index
    ! one greater than the j index of node 1.

    jshift(1,:) = (/ 0,  0,  1,  1,  0,  0,  1,  1/)   
    jshift(2,:) = jshift(1,:)
    jshift(3,:) = (/-1, -1,  0,  0, -1, -1,  0,  0/)   
    jshift(4,:) = jshift(3,:)
    jshift(5,:) = jshift(1,:)
    jshift(6,:) = jshift(1,:)
    jshift(7,:) = jshift(3,:)
    jshift(8,:) = jshift(3,:)

    ! The kshift matrix describes how the k indices of the 8 nodes are related to one another.
    ! E.g, if kshift (1,5) = -1, this means that node 5 has a k index
    ! one less than the k index of node 1.  (Assume that k increases downward.)

    kshift(1,:) = (/ 0,  0,  0,  0, -1, -1, -1, -1/)   
    kshift(2,:) = kshift(1,:)
    kshift(3,:) = kshift(1,:)
    kshift(4,:) = kshift(1,:)
    kshift(5,:) = (/ 1,  1,  1,  1,  0,  0,  0,  0/)
    kshift(6,:) = kshift(5,:)
    kshift(7,:) = kshift(5,:)
    kshift(8,:) = kshift(5,:)

    if (verbose_init) then
       print*, ' '
       print*, 'ishift:'
       do n = 1, 8
          write (6,'(8i4)') ishift(n,:)
       enddo
       print*, ' '
       print*, 'jshift:'
       do n = 1, 8
          write (6,'(8i4)') jshift(n,:)
       enddo
       print*, ' '
       print*, 'kshift:'
       do n = 1, 8
          write (6,'(8i4)') kshift(n,:)
       enddo
    endif

    !----------------------------------------------------------------
    ! Bilinear basis set for reference square, x=(-1,1), y=(-1,1)             
    ! Indexing is counter-clockwise from SW corner
    ! In the code we use "phi_2d" to denote these basis functions. 
    !
    ! N1 = (1-x)*(1-y)/8             N4----N3
    ! N2 = (1+x)*(1-y)/8             |     |
    ! N3 = (1+x)*(1+y)/8             |     |
    ! N4 = (1-x)*(1+y)/8             N1----N2
    !----------------------------------------------------------------

    ! Set coordinates and weights of quadrature points for reference square.
    ! Numbering is counter-clockwise from southwest

    xqp_2d(1) = -rsqrt3; yqp(1) = -rsqrt3
    wqp_2d(1) =  1.d0

    xqp_2d(2) =  rsqrt3; yqp(2) = -rsqrt3
    wqp_2d(2) =  1.d0

    xqp_2d(3) =  rsqrt3; yqp(3) =  rsqrt3
    wqp_2d(3) =  1.d0

    xqp_2d(4) = -rsqrt3; yqp(4) =  rsqrt3
    wqp_2d(4) =  1.d0

    if (verbose_init) then
       print*, ' '
       sumx = 0.d0; sumy = 0.d0; sumz = 0.d0
       do p = 1, nQuadPoints_2d
          print*, xqp_2d(p), yqp_2d(p)
          sumx = sumx + xqp_2d(p); sumy = sumy + yqp_2d(p)
       enddo
       print*, ' '
       print*, 'sums:', sumx, sumy
    endif

    ! Evaluate bilinear basis functions and their derivatives at each quad pt

    do p = 1, nQuadPoints_2d

       phi_2d(1,p) = (1.d0 - xqp_2d(p)) * (1.d0 - yqp_2d(p)) / 4.d0 
       phi_2d(2,p) = (1.d0 + xqp_2d(p)) * (1.d0 - yqp_2d(p)) / 4.d0
       phi_2d(3,p) = (1.d0 + xqp_2d(p)) * (1.d0 + yqp_2d(p)) / 4.d0 
       phi_2d(4,p) = (1.d0 - xqp_2d(p)) * (1.d0 + yqp_2d(p)) / 4.d0

       dphi_dxr_2d(1,p) = -(1.d0 - yqp_2d(p)) / 4.d0 
       dphi_dxr_2d(2,p) =  (1.d0 - yqp_2d(p)) / 4.d0 
       dphi_dxr_2d(3,p) =  (1.d0 + yqp_2d(p)) / 4.d0 
       dphi_dxr_2d(4,p) = -(1.d0 + yqp_2d(p)) / 4.d0

       dphi_dyr_2d(1,p) = -(1.d0 - xqp_2d(p)) / 4.d0 
       dphi_dyr_2d(2,p) = -(1.d0 + xqp_2d(p)) / 4.d0 
       dphi_dyr_2d(3,p) =  (1.d0 + xqp_2d(p)) / 4.d0 
       dphi_dyr_2d(4,p) =  (1.d0 - xqp_2d(p)) / 4.d0 

       if (verbose_init) then
          print*, ' '
          print*, 'Quad point, p =', p
          do n = 1, 4
             print*, n, phi_2d(n,p)
             print*, n, dphi_dxr_2d(n,p)
             print*, n, dphi_dyr_2d(n,p)
          enddo
          print*, 'sum(phi_2d)', sum(phi_2d(:,p))        ! verified that sum = 1
          print*, 'sum(dphi/dx_2d)', sum(dphi_dxr_2d(:,p))  ! verified that sum = 0 (within roundoff)
          print*, 'sum(dphi/dy_2d)', sum(dphi_dyr_2d(:,p))  ! verified that sum = 0 (within roundoff)
       endif

    enddo   ! nQuadPoints_2d

    ! Compute indxA; maps displacements i,j,k = (-1,0,1) onto an index from 1 to 27
    ! Numbering starts in SW corner of layers k-1, finishes in NE corner of layer k+1
    ! Diagonal term has index 14

    ! Layer k-1:           Layer k:            Layer k+1:
    !
    !   7    8    9          16   17   18        25   26   27 
    !   4    5    6          13   14   15        22   23   24
    !   1    2    3          10   11   12        19   20   21                                                                                               
    m = 0
    do k = -1,1
       do j = -1,1
          do i = -1,1
             m = m + 1
             indxA(i,j,k) = m
          enddo
       enddo
    enddo

  end subroutine glissade_velo_higher_init

!****************************************************************************

  subroutine glissade_velo_higher_solve(nx,         ny,           &
                                        nz,         sigma,        &
                                        nhalo,                    &
                                        dx,         dy,           &
                                        thck,       usrf,         &
                                        topg,       eus,          &
!!                                        stagthck, stagusrf,     &
                                        thklim,                   &
                                        flwa,                     &
                                        uvel,     vvel,           &
!                                        beta,                    &
!                                        whichbabc,               &
                                        whichefvs,                &
                                        whichresid,               &
                                        whichnonlinear,           &
                                        whichsparse,              &
                                        whichapprox)

!TODO - Something like this may be needed if building with Trilinos.
!!    use glam_strs2, only: linearSolveTime, totalLinearSolveTime

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------
    !
    ! Note that the glissade solver uses SI units.
    ! Thus we have grid cell dimensions and ice thickness in meters,
    !  velocity in m/s, and the rate factor in Pa^(-n) s(-1).
    !----------------------------------------------------------------

!TODO - Compute nx, ny, and nlyr based on size(thck) and size(flwa)?

    !----------------------------------------------------------------
    ! Note: nx and ny are the horizontal dimensions of scalar arrays (e.g., thck and temp)
    !       The velocity arrays have horizontal dimensions (nx-1, ny-1)
    !       nz is the number of levels at which uvel and vvel are computed
    !       The scalar variables generally live at layer midpoints and have
    !         vertical dimension nz-1.
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! number of grid cells in each direction
       nz,                   &  ! number of vertical levels where velocity is computed
                                ! (same as model%general%upn)
       nhalo                    ! number of rows/columns of halo cells
 
    real(dp), dimension(:), intent(in) :: &
       sigma

    real(dp), intent(in) ::  &
       dx,  dy                  ! grid cell length and width (m)
                                ! assumed to have the same value for each grid cell

    real(dp), dimension(:,:), intent(in) ::  &
       thck,                 &  ! ice thickness (m)
       usrf,                 &  ! upper surface elevation (m)
       topg                     ! elevation of topography (m)

    real(dp), intent(in) ::   &
       eus                      ! eustatic sea level (m)
                                ! = 0. by default

    real(dp), intent(in) ::   & 
       thklim                 ! minimum ice thickness for active cells (m)

    real(dp), dimension(:,:,:), intent(in) ::  &
       flwa                   ! flow factor in units of Pa^(-n) s^(-1)

    real(dp), dimension(:,:,:), intent(inout) ::  &
       uvel, vvel             ! velocity components (m/s)

    !TODO - Uncomment and pass in
!!    real(dp), dimension(:,:), intent(in) ::  &
!!       beta                   ! basal traction parameter

    integer, intent(in) :: whichefvs      ! option for effective viscosity calculation 
                                          ! (calculate it or make it uniform)
    integer, intent(in) :: whichresid     ! option for method to use when calculating residual
    integer, intent(in) :: whichnonlinear ! options for which nonlinear method (Picard or JFNK)
    integer, intent(in) :: whichsparse    ! options for which method for doing elliptic solve
                                          ! (BiCG, GMRES, standalone Trilinos, etc.)
    integer, intent(in), optional :: whichapprox    ! options for which Stokes approximation to use
                                                    ! 0 = SIA, 1 = SSA, 2 = Blatter-Pattyn HO
                                                    ! default = 2

    !--------------------------------------------------------
    ! Local parameters
    !--------------------------------------------------------

    integer, parameter :: cmax = 300          ! max number of outer iterations

    !WHL - What I call resid_target is called minresid in glam_strs2

    real(dp), parameter :: resid_target = 1.0d-04   ! assume velocity fields have converged below this resid 
    real(dp), parameter :: NL_tol   = 1.0d-06   ! to have same criterion as JFNK

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    real(dp), dimension(nx-1,ny-1) :: &
       xVertex, yVertex       ! x and y coordinates of each vertex

    real(dp), dimension(nx-1,ny-1) ::   &
       stagusrf,            & ! upper surface averaged to vertices
       stagthck               ! ice thickness averaged to vertices

    real(dp), dimension(nx,ny) ::    &
       rmask                  ! = 1. where ice is present, else = 0.

    logical, dimension(nx,ny) ::     &
       active_cell,          &! true for active cells (thck > thklim and border locally owned vertices)
       floating_cell,        &! true for cells where ice is present and is floating
       ocean_cell             ! true for cells where topography is below sea level and ice is absent

    logical, dimension(nx-1,ny-1) :: &
       active_vertex          ! true for vertices of active cells

    real(dp), dimension(nz-1,nx,ny) ::  &
       flwafact           ! temperature-based flow factor, 0.5 * A^(-1/n), 
                          ! used to compute effective viscosity
                          ! units: Pa s^(1/n)

    real(dp), dimension(nz,nx-1,ny-1) ::   &
       usav, vsav,                 &! previous guess for velocity solution
       bu, bv,                     &! assembled load vector, divided into 2 parts
       resid_vec_u, resid_vec_v     ! residual vector, divided into 2 parts

    logical, dimension(nz,nx-1,ny-1) ::    &
       umask_dirichlet    ! Dirichlet mask for velocity (if true, u = v = 0)

    real(dp), dimension(27,nz,nx-1,ny-1) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv           ! 1st dimension = 27 (node and its nearest neighbors in x, y and z direction) 
                          ! other dimensions = (k,i,j)

    real(dp) :: &
       resid_velo,          & ! quantity related to velocity convergence
       L2_norm,             & ! L2 norm of residual
       L2_target,           & ! nonlinear convergence target   !TODO - Make sure this is set correctly
       err,                 & ! solution error from sparse_easy_solve
       outer_it_criterion,  & ! current value of outer (nonlinear) loop converence criterion
       outer_it_target        ! target value for outer-loop convergence

    logical, save ::    &
       converged_soln = .false.    ! true if we get a converged solution for velocity

    integer ::  & 
       counter,         & ! outer (nonlinear) iteration counter
       niters             ! linear iteration count from sparse_easy_solve

    real(dp) ::  &
       sia_factor,      & ! = 1. if SIA terms are included, else = 0.
       ssa_factor         ! = 1. if SSA terms are included, else = 0.

    !TODO - Pass in as an argument instead
    real(dp), dimension(nx-1,ny-1) ::  &
       beta                   ! basal traction parameter

    ! The following are used only for the single-processor SLAP solver

    integer ::            &
       nNodesSolve            ! number of nodes where we solve for velocity

    integer, dimension(nz,nx-1,ny-1) ::  &
       NodeID                 ! ID for each node where we solve for velocity
                              ! For periodic BCs, halo node IDs will be copied
                              !  from the other side of the grid

    integer, dimension((nx-1)*(ny-1)*nz) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of nodes

    type(sparse_matrix_type) ::  &
       matrix             ! sparse matrix for SLAP solver, defined in glimmer_sparse_types
                          ! includes nonzeroes, order, col, row, val 

    real(dp), dimension(:), allocatable ::   &   ! for SLAP solver
       rhs,             & ! right-hand-side (b) in Ax = b
       answer,          & ! answer (x) in Ax = b
       resid_vec          ! residual vector Ax - b

    integer ::    &
       matrix_order,    & ! order of matrix = number of rows
       nNonzero           ! upper bound for number of nonzero entries in sparse matrix

!WHL - debug
    integer :: i, j, k, m, n, r
    integer :: iA, jA, kA, colA
    real(dp) :: ds_dx, ds_dy
    real(dp) :: maxthck, maxusrf
    real(dp) :: sumuvel, sumvvel
    logical, parameter :: test_matrix = .false.
    integer, parameter :: test_order = 20
    integer :: rowi
    logical, parameter :: sia_test = .false.


    integer :: nNonzeros    ! number of nonzero entries in structured matrices

    ! Set volume scale
    ! This is not necessary, but dividing by this scale gives matrix coefficients 
    !  that are not quite so large.

    vol0 = dx * dy * 1000.d0    ! typical volume of ice column (m^3)

    if (present(whichapprox)) then
!!       if (whichapprox == HO_APPROX_SIA) then   ! SIA
       if (whichapprox == 0) then   ! SIA
          if (verbose .and. main_task) print*, 'Solving shallow-ice approximation'
          sia_factor = 1.d0
          ssa_factor = 0.d0
!!       elseif (whichapprox == HO_APPROX_SSA) then  ! SSA
       elseif (whichapprox == 1) then  ! SSA
          if (verbose .and. main_task) print*, 'Solving shallow-shelf approximation'
          sia_factor = 0.d0
          ssa_factor = 1.d0
       else   ! Blatter-Pattyn higher-order 
          if (verbose .and. main_task) print*, 'Solving Blatter-Pattyn higher-order approximation'
          sia_factor = 1.d0
          ssa_factor = 1.d0
       endif
    else   ! Blatter-Pattyn higher-order by default
       if (verbose .and. main_task) print*, 'Solving Blatter-Pattyn higher-order approximation'
       sia_factor = 1.d0
       ssa_factor = 1.d0
    endif

    if (test_matrix) then
       call solve_test_matrix(test_order, whichsparse)
    endif

    ! Make sure that the geometry and flow factor are correct in halo cells.
    ! These calls are commented out, since the halo updates are done in 
    !  module glissade.F90, before calling glissade_velo_higher_solve.

!    call parallel_halo(thck)
!    call parallel_halo(topg)
!    call parallel_halo(usrf)
!    call parallel_halo(flwa)

    !------------------------------------------------------------------------------
    ! Setup for higher-order solver: Compute nodal geometry, allocate storage, etc.
    ! These are quantities that do not change during the outer nonlinear loop. 
    !------------------------------------------------------------------------------

    maxthck = maxval(thck(:,:))
    maxthck = parallel_reduce_max(maxthck)
    maxusrf = maxval(usrf(:,:))
    maxusrf = parallel_reduce_max(maxusrf)

    if (verbose .and. main_task) then
       print*, ' '
       print*, 'In glissade_velo_higher_solve'
       print*, 'nx, ny, nz:', nx, ny, nz
       print*, 'vol0:', vol0
       print*, ' '
       print*, 'thklim:', thklim
       print*, 'max thck:', maxthck
       print*, 'max usrf:', maxusrf
    endif

!WHL - debug
    if (verbose .and. this_rank==rtest) then
       print*, ' '
       print*, 'Thickness field, rank =', rtest
       do j = ny, 1, -1
          write(6,'(30f6.0)') thck(3:32,j)
       enddo
       print*, ' '
       print*, 'Topography field, rank =', rtest
       do j = ny, 1, -1
          write(6,'(30f6.0)') topg(3:32,j)
       enddo
       print*, ' '
       print*, 'Upper surface field, rank =', rtest
       do j = ny, 1, -1
          write(6,'(30f6.0)') usrf(3:32,j)
       enddo
       print*, ' '
       print*, 'flwa, k = 1, rank =', rtest
       do j = ny, 1, -1
          write(6,'(30e12.5)') flwa(1,3:32,j)
       enddo
    endif
 
!WHLt2 - If iterating on the velocity during the remapping transport calculation,
!        we do not want these masks to change.
!TODO - Pull these out into a separate subroutine?
    !------------------------------------------------------------------------------
    ! Compute masks: 
    ! (1) real mask = 1 where ice is present, 0 elsewhere
    ! (2) floating mask = .true. where ice is present and is floating
    ! (3) ocean mask = .true. where topography is below sea level and ice is absent
    !------------------------------------------------------------------------------

    do j = 1, ny
       do i = 1, nx

          if (thck(i,j) > thklim) then
	     rmask(i,j) = 1.d0
          else
             rmask(i,j) = 0.d0
          endif

          if (topg(i,j) - eus < (-rhoi/rhoo)*thck(i,j)) then
             if (thck(i,j) > thklim) then
                floating_cell(i,j) = .true.
                ocean_cell(i,j) = .false.
             else
                floating_cell(i,j) = .false.
                ocean_cell(i,j) = .true.
             endif
          else
             floating_cell(i,j) = .false.
             ocean_cell(i,j) = .false.
          endif

       enddo
    enddo

    !------------------------------------------------------------------------------
    ! Compute ice thickness and upper surface on staggered grid
    ! (requires that thck and usrf are up to date in halo cells)
    !------------------------------------------------------------------------------

    call staggered_scalar(nx,           ny,         &
                          nhalo,        rmask,      &
                          thck,         stagthck)

    call staggered_scalar(nx,           ny,         &
                          nhalo,        rmask,      &
                          usrf,         stagusrf)

!WHL - debug
    if (verbose .and. this_rank==rtest) then
       print*, ' '
       print*, 'stagthck, rank =', rtest
       do j = ny-1, 1, -1
          write(6,'(33f6.0)') stagthck(:,j)
       enddo
       print*, ' '
       print*, 'stagusrf, rank =', rtest
       do j = ny-1, 1, -1
          write(6,'(33f6.0)') stagusrf(:,j)
       enddo
    endif

    !------------------------------------------------------------------------------
    ! Compute vertices of each element.
    ! Identify the active cells (i.e., cells with thck > thklim,
    !  bordering a locally owned vertex) and active vertices (all vertices
    !  of active cells).
    !------------------------------------------------------------------------------

    ! For the SLAP solver, count and assign a unique ID to each active node.
    !TODO - Move SLAP-only calculation to another subroutine?

    call get_nodal_geometry(nx,          ny,         nz,  &   
                            nhalo,       sigma,           &
                            dx,          dy,              &
                            thck,        thklim,          &
                            stagusrf,    stagthck,        &
                            xVertex,     yVertex,         &
                            active_cell, active_vertex,   &
                            nNodesSolve, NodeID,          &
                            iNodeIndex,  jNodeIndex,  kNodeIndex)

!WHL - debug
    if (verbose .and. this_rank==rtest) then
       print*, ' '
       print*, 'NodeID before halo update, k = 1:'
       do j = ny-1, 1, -1
          write(6,'(33i6)') NodeID(1,:,j)
       enddo
    endif

    ! Assign the appropriate ID to nodes in the halo.
    ! Note: This will work for single-processor runs with periodic BCs
    !       (e.g., ISMIP-HOM), but not for multiple processors.

    call staggered_parallel_halo(NodeID)

!WHL - debug
    if (verbose .and. this_rank==rtest) then
       print*, ' '
       print*, 'NodeID after halo update, k = 1:'
       do j = ny-1, 1, -1
          write(6,'(33i6)') NodeID(1,:,j)
       enddo
    endif

    !------------------------------------------------------------------------------
    ! Compute the factor A^(-1/n) appearing in the expression for effective viscosity.
    ! Note: The rate factor (flwa = A) is assumed to have units of Pa^(-n) s^(-1).
    !       Thus flwafact = 0.5 * A^(-1/n) has units Pa s^(1/n).
    !------------------------------------------------------------------------------

    flwafact(:,:,:) = 0.d0

    ! Loop over all cells that border locally owned vertices
    !TODO - Simply compute for all cells?  We should have flwa for all cells.

    do j = 1+nhalo, ny-nhalo+1
       do i = 1+nhalo, nx-nhalo+1
          if (active_cell(i,j)) then
             ! gn = exponent in Glen's flow law (= 3 by default)
             flwafact(:,i,j) = 0.5d0 * flwa(:,i,j)**(-1.d0/real(gn,dp))  
          endif
       enddo
    enddo

    if (verbose .and. this_rank==rtest) then
       n = ntest
       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)
       print*, ' '
       print*, 'ntest, i, j, k:', n, i, j, k
    endif
 
    if (verbose .and. this_rank==rtest) then
       print*, ' '
       i = itest; j = jtest; k = ktest
       print*, 'itest, jtest, ktest:', i, j, k
       print*, 'NodeID =', NodeID(k,i,j)
!       print*, ' '
!       print*, 'thck(j+1):', thck(i:i+1,j+1)
!       print*, 'thck(j)  :', thck(i:i+1,j)
!       print*, 'stagthck :', stagthck(i,j)
       print*, ' '
       print*, 'usrf(j+1):', usrf(i:i+1,j+1)
       print*, 'usrf(j)  :', usrf(i:i+1,j)
       print*, 'stagusrf :', stagusrf(i,j)
       print*, ' '
       ds_dx = 0.5d0 * (stagusrf(i,j) + stagusrf(i,j-1) - stagusrf(i-1,j) - stagusrf(i-1,j-1)) / &
               (xVertex(i+1,j) - xVertex(i,j))
       ds_dy = 0.5d0 * (stagusrf(i,j) - stagusrf(i,j-1) + stagusrf(i-1,j) - stagusrf(i-1,j-1)) / &
               (yVertex(i,j+1) - yVertex(i,j))
       print*, 'Finite-difference ds/dx, ds_dy:', ds_dx, ds_dy
       print*, ' '                
       print*, 'flwa, flwafact:', flwa(k,i,j), flwafact(k,i,j)
    endif

    !------------------------------------------------------------------------------
    ! If using SLAP solver, then allocate space for the sparse matrix (A), rhs (b), 
    !  answer (x), and residual vector (Ax-b).
    !------------------------------------------------------------------------------

    if (whichsparse /= STANDALONE_PCG_STRUC .and.    &
        whichsparse /= STANDALONE_TRILINOS_SOLVER) then   

       matrix_order = 2*nNodesSolve    ! Is this exactly enough?
       nNonzero = matrix_order*54      ! 27 = node plus 26 nearest neighbors in hexahedral lattice   
                                       ! 54 = 2 * 27 (since solving for both u and v)

       allocate(matrix%row(nNonzero), matrix%col(nNonzero), matrix%val(nNonzero))
       allocate(rhs(matrix_order), answer(matrix_order), resid_vec(matrix_order))

       if (verbose) then
          print*, 'matrix_order =', matrix_order
          print*, 'nNonzero = ', nNonzero
       endif

    endif   ! SLAP solver
 
    if (verbose .and. this_rank==rtest) then
       print*, 'size(uvel, vvel) =', size(uvel), size(vvel)
       print*, 'sum (uvel, vvel) =', sum(uvel), sum(vvel)
       print*, 'size(usav, vsav) =', size(usav), size(vsav)
    endif

    !---------------------------------------------------------------
    ! Form and solve Ax = b
    !---------------------------------------------------------------

    if (main_task) then
       ! print some info to the screen to update on iteration progress
       print *, ' '
       print *, 'Running Glissade higher-order dynamics solver'
       print *, ' '
       if (whichresid == HO_RESID_L2NORM) then  ! use L2 norm of residual
          print *, 'iter #     resid (L2 norm)       target resid'
       else                                     ! residual based on velocity
          print *, 'iter #     velo resid            target resid'
       end if
       print *, ' '
    endif

    !TODO - Any ghost preprocessing needed, as in glam_strs2?

    !------------------------------------------------------------------------------
    ! set initial values 
    !------------------------------------------------------------------------------

    counter = 0
    resid_velo = 1.d0
    L2_norm   = 1.0d20
    L2_target = 1.0d-4      !TODO - Is this the best value?  
                            ! Should it stay fixed during the outer loop, or should it evolve?

    outer_it_criterion = 1.0d10   ! guarantees at least one loop
    outer_it_target    = 1.0d-12 

    !------------------------------------------------------------------------------
    ! Assemble the load vector b
    ! This goes before the outer loop because the load vector
    !  does not change from one nonlinear iteration to the next.
    !------------------------------------------------------------------------------

    bu(:,:,:) = 0.d0
    bv(:,:,:) = 0.d0

    !------------------------------------------------------------------------------
    ! gravitational forcing
    !------------------------------------------------------------------------------

    call load_vector_gravity(nx,               ny,              &
                             nz,               nhalo,           &
                             nNodesPerElement,                  &
                             sigma,            active_cell,     &
                             xVertex,          yVertex,         &
                             stagusrf,         stagthck,        &
                             bu,               bv)

    !------------------------------------------------------------------------------
    ! lateral water pressure at edge of ice shelves
    !------------------------------------------------------------------------------

!TODO - Test this subroutine

!    call load_vector_shelf_bc(nx,               ny,              &
!                              nz,               nhalo,           &
!                              nNodesPerElement_2d,               &
!                              sigma,                             &
!                              floating_cell,    ocean_cell,      &
!                              xVertex,          yVertex,         &
!                              stagusrf,         stagthck,        &
!                              bu,               bv)

    !------------------------------------------------------------------------------
    ! main outer loop: iteration to solve the nonlinear problem
    !------------------------------------------------------------------------------

    do while (outer_it_criterion >= outer_it_target .and. counter < cmax)

       ! Advance the iteration counter

       counter = counter + 1

       if (verbose .and. main_task) then
          print*, ' '
          print*, 'Outer counter =', counter
          print*, 'whichresid =', whichresid
          print*, 'L2_norm, L2_target =', L2_norm, L2_target
          print*, 'resid_velo, resid_target =', resid_velo, resid_target
       endif

       ! save current velocity
       usav(:,:,:) = uvel(:,:,:)
       vsav(:,:,:) = vvel(:,:,:)

       !---------------------------------------------------------------------------
       ! Assemble the stiffness matrix A
       !---------------------------------------------------------------------------

       if (verbose .and. this_rank==rtest) print*, 'call assemble_stiffness_matrix'

       call assemble_stiffness_matrix(nx,               ny,              &
                                      nz,               nhalo,           &
                                      sigma,            active_cell,     &
                                      xVertex,          yVertex,         &
                                      uvel,             vvel,            &
                                      stagusrf,         stagthck,        &
                                      flwafact,         whichefvs,       &
                                      Auu,              Auv,             &
                                      Avu,              Avv,             &
                                      sia_factor,       ssa_factor)
       
       if (verbose .and. this_rank==rtest) print*, 'Assembled the stiffness matrix'

!WHL - debug
       if (verbose .and. this_rank==rtest) then

          i = itest
          j = jtest
          k = ktest
          n = ntest
          r = rtest
          print*, ' '
          print*, 'Auu after assembly: rank, i, j, k, n =', r, i, j, k, n
!          print*, 'Auv after assembly: rank, i, j, k, n =', r, i, j, k, n
          do kA = -1, 1
          do jA = -1, 1
          do iA = -1, 1
             m = indxA(iA,jA,kA)
             print*, 'iA, jA, kA, Auu:', iA, jA, kA, Auu(m, k, i, j)
!             if (iA==1 .and. jA==0 .and. kA==1) &
!                print*, 'iA, jA, kA, Auv:', iA, jA, kA, Auvm, k, i, j)
          enddo
          enddo
          enddo

       endif

       !---------------------------------------------------------------------------
       ! Incorporate basal sliding boundary conditions
       ! Assume a linear sliding law for now

       ! Note: We could call this subroutine before the main outer loop if beta
       !       is assumed to be independent of velocity.  Putting the call here,
       !       however, allows for more general sliding laws.
       !---------------------------------------------------------------------------

       ! Set beta to zero to deactivate this subroutine for now
       !TODO = Pass in as argument
       beta(:,:) = 0.d0

!TODO - Test this subroutine

!       call basal_sliding_bc(nx,               ny,              &
!                             nz,               nhalo,           &
!                             nNodesPerElement_2d,               &
!                             active_cell,      beta,            &
!                             xVertex,          yVertex,         &
!                             Auu,              Avv)

          i = itest
          j = jtest
          k = ktest
          n = ntest
!          print*, ' '
!          print*, 'Auu after basal sliding BC: rank, i, j, k, n =', r, i, j, k, n
          do kA = -1, 1
          do jA = -1, 1
          do iA = -1, 1
              m = indxA(iA,jA,kA)
!             print*, 'iA, jA, kA, Auu:', iA, jA, kA, Auu(m, k, i, j)
          enddo
          enddo
          enddo

       !---------------------------------------------------------------------------
       ! Incorporate Dirichlet (u = 0) basal boundary conditions
       !---------------------------------------------------------------------------

       ! For now, simply impose zero sliding everywhere at the bed.
       ! TODO: Modify for confined shelf.  Read into subroutine?                                                                                            
       umask_dirichlet(1:nz-1,:,:) = .false.
       umask_dirichlet(  nz,  :,:) = .true.   ! u = v = 0 at bed

       if (verbose .and. main_task) then
          print*, ' '
          print*, 'Call dirichlet_bc'
       endif

       call dirichlet_boundary_conditions(nx,              ny,                &
                                          nz,              nhalo,             &
                                          active_vertex,   umask_dirichlet,   &
                                          Auu,             Auv,               &
                                          Avu,             Avv,               &
                                          bu,              bv)

       !---------------------------------------------------------------------------
       ! Halo updates for matrices
       !
       ! These updates are not strictly necessary unless we're concerned about
       !  roundoff errors.
       ! But suppose we are comparing two entries that are supposed to be equal
       !  (e.g., to preserve symmetry), where entry 1 is owned by processor A and 
       !  entry 2 is owned by processor B.  
       ! Processor A might compute a local version of entry 2 in its halo, with 
       !  entry 2 = entry 1 locally.  But processor B's entry 2 might be different
       !  because of roundoff.  We need to make sure that processor B's value of 
       !  is communicated to processor A.  If these values are slightly different, 
       !  they will be reconciled by the subroutine check_symmetry_assembled_matrix.
       !---------------------------------------------------------------------------
     

!WHL - debug
    if (verbose .and. this_rank==rtest) then
!       print*, ' '
!       m = indxA(0,0,0)
!       print*, 'Before halo update, Auu(0,0,0,1,:,:), rank =', rtest
!       do j = ny-1, 1, -1
!          write(6,'(23e10.3)') Auu(m,1,:,j)
!       enddo
    endif

        call staggered_parallel_halo(Auu(:,:,:,:))
        call staggered_parallel_halo(Auv(:,:,:,:))
        call staggered_parallel_halo(Avu(:,:,:,:))
        call staggered_parallel_halo(Avv(:,:,:,:))

!WHL - debug
    if (verbose .and. this_rank==rtest) then
!       print*, ' '
!       print*, 'After halo update, Auu_diag(1,:,:), rank =', rtest
!       m = indxA(0,0,0)
!       do j = ny-1, 1, -1
!          write(6,'(23e10.3)') Auu(m,1,:,j)
!       enddo
    endif

       !---------------------------------------------------------------------------
       ! Check symmetry of assembled matrix
       ! 
       ! There may be small differences from perfect symmetry due to roundoff errors.  
       ! If sufficiently small, these differences are fixed by averaging the two values 
       !  that should be symmetric.  Otherwise the code aborts.
       !
       ! Note: It may be OK to skip this check for production code.  However,
       !       small violations of symmetry are not tolerated well by some solvers.
       !       For example, the SLAP PCG solver with incomplete Cholesky preconditioning
       !       can crash if symmetry is not perfect. 
       !---------------------------------------------------------------------------

       if (verbose .and. main_task) print*, 'Check matrix symmetry'

       if (check_symmetry) then

          call check_symmetry_assembled_matrix(nx,          ny,      &
                                               nz,          nhalo,   &
                                               active_vertex,        &
                                               Auu,         Auv,     &
                                               Avu,         Avv)

       endif

!WHL - debug - Try stripping out columns from the matrix and see if it can still solve an SIA problem.

          if (sia_test) then

             if (verbose .and. main_task) then
                print*, ' '
                print*, 'SIA test: Removing Auv, Avu, and iA/jA columns'
             endif

             ! Set Auv = Avu = 0

             Auv(:,:,:,:) = 0.d0 
             Avu(:,:,:,:) = 0.d0 

             ! Remove terms from iA and jA columns

             do j = 1, ny-1
             do i = 1, nx-1
             do k = 1, nz
                do kA = -1,1
                do jA = -1,1
                do iA = -1,1
                   if (iA /= 0 .or. jA /= 0) then
                      m = indxA(iA,jA,kA)
                      Auu(m,k,i,j) = 0.d0
                      Avv(m,k,i,j) = 0.d0
                   endif
                enddo
                enddo
                enddo
             enddo
             enddo
             enddo

          endif   ! sia_test             

       !---------------------------------------------------------------------------
       ! Count nonzero elements in structured matrices (strictly diagnostic)
       !---------------------------------------------------------------------------

       nNonzeros = 0
       do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (active_vertex(i,j)) then
             do k = 1, nz
                do kA = -1, 1
                do jA = -1, 1
                do iA = -1, 1
                   m = indxA(iA,jA,kA)
                   if (Auu(m,k,i,j) /= 0.d0) nNonzeros = nNonzeros + 1
                   if (Auv(m,k,i,j) /= 0.d0) nNonzeros = nNonzeros + 1
                   if (Avu(m,k,i,j) /= 0.d0) nNonzeros = nNonzeros + 1
                   if (Avv(m,k,i,j) /= 0.d0) nNonzeros = nNonzeros + 1
                enddo 
                enddo
                enddo
             enddo  ! k
          endif     ! active_vertex
       enddo        ! i
       enddo        ! j

       if (verbose .and. this_rank==rtest) then
          print*, ' '
          print*, 'nNonzeros in structured matrices =', nNonzeros
       endif

!WHL - debug - print out some matrix values for test point

       if (verbose .and. this_rank==rtest) then
          i = itest
          j = jtest
          k = ktest
          n = ntest

          print*, ' '
      	  print*, 'i,j,k,n =', i, j, k, n
          print*, 'Auu sum =', sum(Auu(:,k,i,j))
          print*, 'Auv sum =', sum(Auv(:,k,i,j))
          print*, 'Avu sum =', sum(Avu(:,k,i,j))
          print*, 'Avv sum =', sum(Avv(:,k,i,j))

          print*, ' '
          print*, 'iA, jA, kA, Auu:'
          do kA = -1, 1
          do jA = -1, 1
          do iA = -1, 1
             m = indxA(iA,jA,kA)
             print*, iA, jA, kA, Auu(m,k,i,j) 
          enddo
          enddo
          enddo
     
          do kA = -1, 1
          do jA = -1, 1
          do iA = -1, 1
             if ( (k+kA >= 1 .and. k+kA <= nz)       &
                             .and.                   &
                  (i+iA >= 1 .and. i+iA <= nx-1)     &
                             .and.                   &
                  (j+jA >= 1 .and. j+jA <= ny-1) ) then
                colA = NodeID(k+kA, i+iA, j+jA)   ! ID for neighboring node
                m = indxA(iA,jA,kA)
                print*, ' '
                print*, 'iA, jA, kA, colA =', iA, jA, kA, colA
                print*, 'i, j, k for colA:', iNodeIndex(colA), jNodeIndex(colA), kNodeIndex(colA)
                print*, 'Auu =', Auu(m,k,i,j)
                print*, 'Auv =', Auv(m,k,i,j)
                print*, 'Avu =', Avu(m,k,i,j)
                print*, 'Avv =', Avv(m,k,i,j)
             endif
          enddo
          enddo
          enddo

       endif   ! verbose, this_rank==rtest

       if (whichsparse == STANDALONE_PCG_STRUC) then   ! standalone PCG for structured grid
                                                       ! works for both serial and parallel runs

          !------------------------------------------------------------------------
          ! Compute the residual vector and its L2 norm
          !------------------------------------------------------------------------

          if (verbose .and. this_rank==rtest) then
             print*, ' '
             print*, 'Compute residual vector'
          endif

          call compute_residual_vector(nx,          ny,            &
                                       nz,          nhalo,         &
                                       active_vertex,              &
                                       Auu,         Auv,           &
                                       Avu,         Avv,           &
                                       bu,          bv,            &
                                       uvel,        vvel,          &
                                       resid_vec_u, resid_vec_v,   &
                                       L2_norm)

          if (verbose .and. this_rank==rtest) then
             print*, ' '
             print*, 'counter, L2_norm =', counter, L2_norm
          endif

          ! preprocessing is not needed because the matrix, rhs, and solution already
          ! have the required data structure

          if (verbose .and. this_rank==rtest) then
             print*, ' '
             print*, 'Calling structured PCG solver'
          endif

          !------------------------------------------------------------------------
          ! Call linear PCG solver, compute uvel and vvel on local processor
          !------------------------------------------------------------------------

          call pcg_solver_structured(nx,        ny,            &
                                     nz,        nhalo,         &
                                     indxA,     active_vertex, &
                                     Auu,       Auv,           &
                                     Avu,       Avv,           &
                                     bu,        bv,            &
                                     uvel,      vvel,          &
                                     err,       niters)

          ! Halo updates for uvel and vvel
          !TODO - Put these after the 'endif'?

          call staggered_parallel_halo(uvel)
          call staggered_parallel_halo(vvel)
!          call horiz_bcs_stag_vector_ew(uvel)
!          call horiz_bcs_stag_vector_ns(vvel)

       elseif (whichsparse /= STANDALONE_TRILINOS_SOLVER) then   ! one-processor SLAP solve   
          
          !------------------------------------------------------------------------
          ! Given the stiffness matrices (Auu, etc.) and load vector (bu, bv) in
          !  structured format, form the global matrix and rhs in SLAP format.
          !------------------------------------------------------------------------

          if (verbose) print*, 'Form global matrix in sparse format'
 
          matrix%order = matrix_order
          matrix%nonzeros = nNonzero
          matrix%symmetric = .false.

          call slap_preprocess(nx,           ny,          &   
                               nz,           nNodesSolve, &
                               NodeID,                    &
                               iNodeIndex,   jNodeIndex,  &
                               kNodeIndex,                &
                               Auu,          Auv,         &
                               Avu,          Avv,         &
                               bu,           bv,          &
                               uvel,         vvel,        &
                               matrix_order, nNonzero,    &
                               matrix,       rhs,         &
                               answer)

          !------------------------------------------------------------------------
          ! Compute the residual vector and its L2_norm
          !------------------------------------------------------------------------

          call slap_compute_residual_vector(matrix,    answer,   rhs,  &
                                            resid_vec, L2_norm)

!WHL - bug check
          if (verbose) then
             print*, ' '
             print*, 'Before linear solve: n, row, col, val:'
             print*, ' '
             do n = 1, 10              
                print*, n, matrix%row(n), matrix%col(n), matrix%val(n)
             enddo

             print*, 'L2_norm of residual =', L2_norm
             print*, 'Call sparse_easy_solve, counter =', counter
          endif

          !------------------------------------------------------------------------
          ! Solve the linear matrix problem
          !------------------------------------------------------------------------

          call sparse_easy_solve(matrix, rhs,    answer,  &
                                 err,    niters, whichsparse)

          if (verbose) then
             print*, 'Solved the linear system, niters, err =', niters, err
             print*, ' '
             print*, 'n, u, v (m/yr):', ntest, scyr*answer(2*ntest-1), scyr*answer(2*ntest)
!!           do n = 1, matrix%order 
!!              print*, n, answer(n)
!!           enddo
          endif

          !------------------------------------------------------------------------
          ! Put the velocity solution back into 3D arrays
          !------------------------------------------------------------------------

          call slap_postprocess(nNodesSolve,  answer,                   &
                                iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                                uvel,         vvel)

          !------------------------------------------------------------------------
          ! Halo updates for uvel and vvel
          !------------------------------------------------------------------------

          call staggered_parallel_halo(uvel)
          call staggered_parallel_halo(vvel)
!          call horiz_bcs_stag_vector_ew(uvel)
!          call horiz_bcs_stag_vector_ns(vvel)

#ifdef TRILINOS
       else    ! solve with Trilinos

!WHL - Commented out for now until testing later with Trilinos
!!           call solvewithtrilinos(rhs, answer, linearSolveTime)
!!           totalLinearSolveTime = totalLinearSolveTime + linearSolveTime
          !  write(*,*) 'Total linear solve time so far', totalLinearSolveTime                                           
#endif

       endif   ! whichsparse (STANDALONE_PCG_STRUC or not) 

       ! Any ghost postprocessing needed, as in glam_strs2?

!WHL - bug check
       if (verbose .and. this_rank==rtest) then
          i = itest
          j = jtest
          print*, ' '
          print*, 'After postprocess: rank, i, j =', this_rank, i, j
          print*, 'k, uvel, vvel (m/yr):'
          do k = 1, nz
             print*, k, scyr*uvel(k,i,j), scyr*vvel(k,i,j)               
          enddo
       endif

       !---------------------------------------------------------------------------
       ! Compute residual quantities based on the velocity solution
       !---------------------------------------------------------------------------

       call compute_residual_velocity(nhalo,  whichresid,   &
                                      uvel,   vvel,        &
                                      usav,   vsav,        &
                                      resid_velo)
  
       if (verbose .and. this_rank==rtest) then
          print*, ' '
          print*, 'whichresid, resid_velo =', whichresid, resid_velo
       endif

       !---------------------------------------------------------------------------
       ! Write diagnostics (iteration number, max residual, and location of max residual
       ! (send output to the screen or to the log file, per whichever line is commented out) 
       !---------------------------------------------------------------------------

       !TODO - Find out if this note from glam_strs2 still applies:
       ! "Can't use main_task flag because main_task is true for all processors in case of parallel_single"

       if (main_task) then
          if (whichresid == HO_RESID_L2NORM) then
             if (verbose) then
                print*, ' '
                print*, 'Using L2 norm convergence criterion'
                print*, 'iter#, L2 norm, L2 target:', counter, L2_norm, L2_target
             else   ! standard short message
                print '(i4,2g20.6)', counter, L2_norm, L2_target
             endif
             !write(message,'(i4,3g20.6)') counter, L2_norm, L2_target
             !call write_log (message)
          else
             if (verbose) then
                print*, ' '
                print*, 'Using velocity residual convergence criterion'
                print*, 'iter#, resid velo, resid target:', counter, resid_velo, resid_target
             else
                print '(i4,2g20.6)', counter, resid_velo, resid_target
             endif
             !write(message,'(i4,2g20.6)') counter, resid_velo, resid_target
             !call write_log (message)
          end if
       endif

       !---------------------------------------------------------------------------
       ! update the outer loop stopping criterion
       !---------------------------------------------------------------------------

       if (whichresid == HO_RESID_L2NORM) then
          outer_it_criterion = L2_norm
          outer_it_target = L2_target      ! L2_target is currently set to 1.d-4 and held constant
       else
          outer_it_criterion = resid_velo
          outer_it_target = resid_target   ! resid_target is currently a parameter = 1.d-4  
       end if

    enddo  ! while (outer_it_criterion >= outer_it_target .and. counter < cmax)

    if (counter < cmax) converged_soln = .true.

    if (verbose .and. converged_soln) then

       if (main_task) then
          print*, ' '
          print*, 'GLISSADE SOLUTION HAS CONVERGED, outer counter =', counter
       endif

    endif

    !------------------------------------------------------------------------------
    ! Clean up
    !------------------------------------------------------------------------------

    if (whichsparse /= STANDALONE_PCG_STRUC .and.   &
        whichsparse /= STANDALONE_TRILINOS_SOLVER) then
       deallocate(matrix%row, matrix%col, matrix%val)
       deallocate(rhs, answer, resid_vec)
    endif

  end subroutine glissade_velo_higher_solve

!****************************************************************************

  subroutine staggered_scalar(nx,           ny,        &
                              nhalo,        rmask,     &
                              var,          stagvar)

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nhalo                    ! number of rows/columns of halo cells

    real(dp), dimension(nx,ny), intent(in) ::        &
       rmask                    ! = 1. where ice is present, else = 0.

    real(dp), dimension(nx,ny), intent(in) ::    &
       var                      ! scalar field defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       stagvar                  ! staggered field defined at cell corners

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j

    real(dp) :: sumvar, summask

    stagvar(:,:) = 0.d0

    ! Loop over all vertices of locally owned cells
    ! Average the input field over the neighboring cells with ice present (rmask = 1)

    do j = 1, ny-1     ! all vertices
       do i = 1, nx-1
          sumvar = rmask(i,j+1)*var(i,j+1) + rmask(i+1,j+1)*var(i+1,j+1)  &
                 + rmask(i,j)  *var(i,j)   + rmask(i+1,j)  *var(i+1,j)	  
          summask = rmask(i,j+1) + rmask(i+1,j+1) + rmask(i,j) + rmask(i+1,j)
          if (summask > 0.d0) stagvar(i,j) = sumvar / summask
       enddo
    enddo  

  end subroutine staggered_scalar

!****************************************************************************

  subroutine get_nodal_geometry(nx,          ny,          nz,      &   
                                nhalo,       sigma,                &
                                dx,          dy,                   &
                                thck,        thklim,               &
                                stagusrf,    stagthck,             &
                                xVertex,     yVertex,              &
                                active_cell, active_vertex,        &
                                nNodesSolve, NodeID,               & 
                                iNodeIndex,  jNodeIndex,  kNodeIndex)
                            
    !----------------------------------------------------------------
    ! Compute coordinates for each vertex.
    ! Identify and count the active cells and vertices for the finite-element calculations.
    ! Active cells include all cells that contain ice (thck > thklin) and border locally owned vertices.
    ! Active vertices include all vertices of active cells.
    !
    ! Also compute some node indices needed for the SLAP single-processor solver.
    !TODO - Move SLAP part to different subroutine?
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx,  ny,              &    ! number of grid cells in each direction
       nz,                   &    ! number of vertical levels where velocity is computed
       nhalo                      ! number of halo layers

    real(dp), dimension(nz), intent(in) :: &
       sigma                  ! sigma vertical coordinate

    real(dp), intent(in) ::  &
       dx,  dy                ! grid cell length and width (m)
                              ! assumed to have the same value for each grid cell

    real(dp), dimension(nx,ny), intent(in) ::  &
       thck                   ! ice thickness

    real(dp), intent(in) ::   & 
       thklim                 ! minimum ice thickness for active cells

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf,            & ! upper surface averaged to vertices
       stagthck               ! ice thickness averaged to vertices

    real(dp), dimension(nx-1,ny-1), intent(out) :: &
       xVertex, yVertex       ! x and y coordinates of each vertex

    logical, dimension(nx,ny), intent(out) :: &
       active_cell            ! true for active cells 
                              ! (thck > thklim, bordering a locally owned vertex)

    logical, dimension(nx-1,ny-1), intent(out) :: &
       active_vertex          ! true for vertices of active cells

    ! The remaining input/output arguments are for the SLAP solver

    integer, intent(out) :: &
       nNodesSolve            ! number of nodes where we solve for velocity

    integer, dimension(nz,nx-1,ny-1), intent(out) ::  &
       NodeID                 ! ID for each node where we solve for velocity

    integer, dimension((nx-1)*(ny-1)*nz), intent(out) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of nodes

    !---------------------------------------------------------
    ! Local variables
    !---------------------------------------------------------

    integer :: i, j, k

    !----------------------------------------------------------------
    ! Compute the x and y coordinates of each vertex.
    ! By convention, vertex (i,j) lies at the NE corner of cell(i,j).
    !----------------------------------------------------------------

    xVertex(:,:) = 0.d0
    yVertex(:,:) = 0.d0
    do j = 1, ny-1
    do i = 1, nx-1
       xVertex(i,j) = dx * i
       yVertex(i,j) = dy * j
    enddo
    enddo

    if (verbose_init) then
       write(6,*) ' '
       write(6,*) 'Vertex coordinates:'
       do j = 1, ny-1
          write(6,*) ' '    
          do i = 1, nx-1
             write(6,*) 'i, j, xV, yV:', i, j, xVertex(i,j), yVertex(i,j)
          enddo
       enddo
    endif

    ! Identify the active cells.
    ! Include all cells that border locally owned vertices and contain ice.

    active_cell(:,:) = .false.

    do j = 1+nhalo, ny-nhalo+1
    do i = 1+nhalo, nx-nhalo+1
       if (thck(i,j) > thklim) then
          active_cell(i,j) = .true.
       endif
    enddo
    enddo

    ! Identify the active vertices.
    ! Include all vertices of active cells.

    active_vertex(:,:) = .false.

    do j = 1+nhalo, ny-nhalo+1  ! include east and north layer of halo cells
    do i = 1+nhalo, nx-nhalo+1
       if (active_cell(i,j)) then
          active_vertex(i-1:i, j-1:j) = .true.  ! all vertices of this cell
       endif
    enddo
    enddo

    ! Identify and count the nodes where we must solve for the velocity.
    ! This indexing is used for pre- and post-processing of the assembled matrix
    !  when we call the SLAP solver (one processor only).
    ! It is not required by the structured solver.
    !TODO - Move to separate subroutine?

    nNodesSolve = 0
    NodeID(:,:,:) = 0
    iNodeIndex(:) = 0
    jNodeIndex(:) = 0
    kNodeIndex(:) = 0

    do j = nhalo+1, ny-nhalo    ! locally owned vertices only
    do i = nhalo+1, nx-nhalo
       if (active_vertex(i,j)) then   ! all nodes in column are active
          do k = 1, nz               
             nNodesSolve = nNodesSolve + 1   
             NodeID(k,i,j) = nNodesSolve   ! unique index for each node
             iNodeIndex(nNodesSolve) = i
             jNodeIndex(nNodesSolve) = j
             kNodeIndex(nNodesSolve) = k

             if (verbose .and. this_rank==rtest .and. nNodesSolve==ntest) then
                print*, ' '
                print*, 'i, j, k, n:', i, j, k, nNodesSolve
                print*, 'sigma, stagusrf, stagthck:', sigma(k), stagusrf(i,j), stagthck(i,j)
                print*, 'dx, dy, dz:', xVertex(i,j) - xVertex(i-1,j), &
                                       yVertex(i,j) - yVertex(i,j-1), &
                                       (sigma(k+1) - sigma(k)) * stagthck(i,j)
             endif

           enddo   ! k
        endif      ! active vertex
    enddo          ! i
    enddo          ! j

    if (verbose .and. this_rank==rtest) then
       print*, ' '
       print*, 'nNodesSolve =', nNodesSolve
    endif

  end subroutine get_nodal_geometry

!****************************************************************************

  subroutine load_vector_gravity(nx,               ny,              &
                                 nz,               nhalo,           &
                                 nNodesPerElement,                  &
                                 sigma,            active_cell,     &
                                 xVertex,          yVertex,         &
                                 stagusrf,         stagthck,        &
                                 bu,               bv)

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical layers at which velocity is computed
                                     ! Note: the number of elements per column is nz-1
       nhalo                         ! number of halo layers

    integer, intent(in) :: nNodesPerElement

    real(dp), dimension(nz), intent(in) ::    &
       sigma                         ! sigma vertical coordinate

    logical, dimension(nx,ny), intent(in) ::  &
       active_cell                   ! true if cell contains ice and borders a locally owned vertex

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xVertex, yVertex     ! x and y coordinates of vertices

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf,       &  ! upper surface elevation on staggered grid (m)
       stagthck           ! ice thickness on staggered grid (m)

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::  &
       bu, bv             ! load vector, divided into u and v components

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    real(dp), dimension(nNodesPerElement) ::     &
       x, y, z,         & ! Cartesian coordinates of nodes
       s                  ! upper surface elevation at nodes

    real(dp)  ::   &
       detJ               ! determinant of Jacobian for the transformation
                          !  between the reference element and true element

    real(dp), dimension(nNodesPerElement) ::   &
       dphi_dx, dphi_dy, dphi_dz   ! derivatives of basis functions, evaluated at quad pts

    real(dp) ::    &
       ds_dx, ds_dy       ! horizontal gradient of upper surface elevation

    integer :: i, j, k, n, p

    integer :: iNode, jNode, kNode

    ! Sum over elements in active cells 
    ! Loop over all cells that border locally owned vertices

    do j = nhalo+1, ny-nhalo+1
    do i = nhalo+1, nx-nhalo+1
       
       if (active_cell(i,j)) then

          do k = 1, nz-1    ! loop over elements in this column 
                            ! assume k increases from upper surface to bed

             if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
                print*, 'i, j, k:', i, j, k
                print*, ' '
             endif

             ! compute spatial coordinates, velocity, and upper surface elevation for each node

             do n = 1, nNodesPerElement

                ! Determine (k,i,j) for this node
                ! The reason for the '7' is that node 7, in the NE corner of the upper layer, has index (k,i,j).
                ! Indices for other nodes are computed relative to this node.
                iNode = i + ishift(7,n)
                jNode = j + jshift(7,n)
                kNode = k + kshift(7,n)

                x(n) = xVertex(iNode,jNode)
                y(n) = yVertex(iNode,jNode)
                z(n) = stagusrf(iNode,jNode) - sigma(kNode)*stagthck(iNode,jNode)
                s(n) = stagusrf(iNode,jNode)
 
                if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
                   print*, 'i, j, k, n, x, y, z, s:', i, j, k, n, x(n), y(n), z(n), s(n)
                endif

             enddo   ! nodes per element

             ! Loop over quadrature points for this element
   
             do p = 1, nQuadPoints

                ! Evaluate the derivatives of the element basis functions
                ! at this quadrature point.

!WHL - debug - Pass in i, j, k, and p for now

                call get_basis_function_derivatives(nNodesPerElement,                             &
                                                    x(:),          y(:),          z(:),           &
                                                    dphi_dxr(:,p), dphi_dyr(:,p), dphi_dzr(:,p),  &
                                                    dphi_dx(:),    dphi_dy(:),    dphi_dz(:),     &
                                                    detJ , i, j, k, p                      )

                ! Evaluate ds_dx and ds_dy at this quadrature point
 
                ds_dx = 0.d0
                ds_dy = 0.d0

                do n = 1, nNodesPerElement
                   ds_dx = ds_dx + dphi_dx(n) * s(n)
                   ds_dy = ds_dy + dphi_dy(n) * s(n)
                enddo

                if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
                   print*, ' '
                   print*, 'Increment load vector, i, j, k, p =', i, j, k, p
                   print*, 'ds/dx, ds/dy =', ds_dx, ds_dy
                   print*, 'detJ/vol0 =', detJ/vol0
                   print*, 'detJ/vol0* (ds/dx, ds/dy) =', detJ/vol0*ds_dx, detJ/vol0*ds_dy
                endif

                ! Increment the load vector with the gravitational contribution from
                ! this quadrature point.

                do n = 1, nNodesPerElement

                   ! Determine (k,i,j) for this node
                   iNode = i + ishift(7,n)
                   jNode = j + jshift(7,n)
                   kNode = k + kshift(7,n)
         
                   bu(kNode,iNode,jNode) = bu(kNode,iNode,jNode) - rhoi*grav * wqp(p) * detJ/vol0 * ds_dx * phi(n,p)
                   bv(kNode,iNode,jNode) = bv(kNode,iNode,jNode) - rhoi*grav * wqp(p) * detJ/vol0 * ds_dy * phi(n,p)

                   if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest .and. p==ptest) then
!                      print*, ' '
                      print*, 'n, phi(n), delta(bu), delta(bv):', n, phi(n,p), &
                               rhoi*grav*wqp(p)*detJ/vol0 * ds_dx * phi(n,p), &
                               rhoi*grav*wqp(p)*detJ/vol0 * ds_dy * phi(n,p)
                   endif

                enddo   ! nNodesPerElement

             enddo      ! nQuadPoints

          enddo         ! k

       endif            ! active cell

    enddo               ! i
    enddo               ! j

  end subroutine load_vector_gravity

!****************************************************************************

  subroutine load_vector_shelf_bc(nx,               ny,              &
                                  nz,               nhalo,           &
                                  nNodesPerElement_2d,               &
                                  sigma,                             &
                                  floating_cell,    ocean_cell,      &
                                  xVertex,          yVertex,         &
                                  stagusrf,         stagthck,        &
                                  bu,               bv)

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical layers at which velocity is computed
                                     ! Note: the number of elements per column is nz-1
       nhalo                         ! number of halo layers

    integer, intent(in) :: nNodesPerElement_2d

    real(dp), dimension(nz), intent(in) ::    &
       sigma                         ! sigma vertical coordinate

    logical, dimension(nx,ny), intent(in) ::  &
       floating_cell,               &! true if ice is present and is floating
       ocean_cell                    ! true if topography is below sea level and ice is absent

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xVertex, yVertex     ! x and y coordinates of vertices

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf,       &  ! upper surface elevation on staggered grid (m)
       stagthck           ! ice thickness on staggered grid (m)

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::  &
       bu, bv             ! load vector, divided into u and v components

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j

    ! Sum over elements in active cells 
    ! Loop over cells that contain locally owned vertices
    !TODO - Make sure we don't step out of bounds
    !       With more care, we could skip some computations for vertices that are not locally owned.

    do j = nhalo+1, ny-nhalo+1
    do i = nhalo+1, nx-nhalo+1
       
       if (floating_cell(i,j)) then   ! ice is present and is floating

          if (ocean_cell(i-1,j)) then ! compute lateral BC for west face

             call lateral_shelf_bc(nx,                   ny,          &
                                   nz,                   sigma,       &
                                   nNodesPerElement_2d,  ' west',     &
                                   i,                    j,           &
                                   stagusrf,             stagthck,    &
                                   xVertex,              yVertex,     &
                                   bu,                   bv)

          endif

          if (ocean_cell(i+1,j)) then ! compute lateral BC for east face

          endif

          if (ocean_cell(i,j-1)) then ! compute lateral BC for south face

          endif

          if (ocean_cell(i+1,j)) then ! compute lateral BC for north face

          endif

       endif      ! floating_cell

    enddo         ! i
    enddo         ! j



!TODO - Add basal BC here (or in a separate subroutine).  Start with linear sliding law.
!       Would be similar to lateral shelf BC in terms of doing surface integrals.
!       We would integrate over all faces that contain at least one node with a stress BC. (Not Dirichlet or free-slip)
!       Note that Dirichlet BCs are enforced after matrix assembly. 

  end subroutine load_vector_shelf_bc

!****************************************************************************

  subroutine lateral_shelf_bc(nx,                  ny,           &
                              nz,                  sigma,        &
                              nNodesPerElement_2d, face,         &
                              iCell,               jCell,        &
                              stagusrf,            stagthck,     &
                              xVertex,             yVertex,      &
                              bu,                  bv)

!TODO - Test this subroutine.
!       Verify that the sum of contributions to bu or bv over an entire column
!        is equal to pw_av times the surface area of the column.

    !----------------------------------------------------------------------------------
    ! Determine the contribution to the load vector from water pressure at shelf edges.
    ! This subroutine computes the vertically integrated contribution from a single shelf face.
    ! At a given point, this contribution is proportional to the water pressure pw.
    !
    ! It can be shown that the vertically averaged water pressure over a floating column is
    ! 
    !              pw_av = 0.5*rhoi*grav*(H - s)
    ! 
    ! where H is the ice thickness, s is the surface elevation, and rhoi is the density of ice.
    !
    ! Here we make the simplifying assumption that the water pressure everywhere in the column
    ! is equal to the vertically averaged water pressure.
    ! Thus when we sum over quadrature points, the contribution from each quadrature point is
    ! proportional to (H - s) at that point.
    !
    ! A more correct calculation would use the actual water pressure at each quadrature point.
    ! This is given by pw = rhoo*grav*(-z), where z is the elevation, with sea level at z = 0.
    !-----------------------------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical layers at which velocity is computed
                                     ! Note: the number of elements per column is nz-1
       iCell, jCell,            &    ! i and j indices for this cell
       nNodesPerElement_2d           ! nodes per 2D (square) element face

    character(len=5), intent(in) ::  &
       face                          ! 'north', 'south', ' east', or ' west'
 
    real(dp), dimension(nz), intent(in) ::    &
       sigma                         ! sigma vertical coordinate

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xVertex, yVertex   ! x and y coordinates of vertices

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf,       &  ! upper surface elevation on staggered grid (m)
       stagthck           ! ice thickness on staggered grid (m)

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::  &
       bu, bv             ! load vector, divided into u and v components

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    real(dp), dimension(nNodesPerElement_2d) ::     &
       x, y,               & ! local coordinates of nodes
       s,                  & ! upper surface elevation at nodes
       h                     ! ice thickness at nodes

    integer, dimension(nNodesPerElement_2d) ::     &
       iNode, jNode, kNode   ! global indices of each node

    !TODO - These are not currently used except as dummy arguments
    real(dp), dimension(nNodesPerElement_2d) ::   &
       dphi_dx_2d, dphi_dy_2d    ! derivatives of basis functions, evaluated at quad pts

    real(dp)  ::        &
       dz,              & ! depth of ice below the sea surface, H - s 
       detJ               ! determinant of Jacobian for the transformation
                          !  between the reference element and true element

    integer :: k, n, p

    ! Compute nodal geometry in a local xy reference system
    ! Note: The local y direction is really the vertical direction
    !       The local x direction depends on the face (N/S/E/W)
    ! The diagrams below show the node indexing convention, along with the true
    !  directions for each face.  The true directions are mapped to local (x,y).

    if (trim(face)=='west') then

       !     4-----3       z
       !     |     |       ^
       !     |     |       |
       !     1-----2       ---> -y

       x(1) = yvertex(iCell-1, jCell)
       x(2) = yvertex(iCell-1, jCell-1)
       x(3) = x(2)
       x(4) = x(1)

       s(1) = stagusrf(iCell-1, jCell)
       s(2) = stagusrf(iCell-1, jCell-1)
       s(3) = s(2)
       s(4) = s(1)

       h(1) = stagthck(iCell-1, jCell)
       h(2) = stagthck(iCell-1, jCell-1)
       h(3) = h(2)
       h(4) = h(1)

       iNode(1) = iCell-1
       jNode(1) = jCell

       iNode(2) = iCell-1
       jNode(2) = jCell-1

       iNode(3) = iCell-1
       jNode(3) = jCell-1

       iNode(4) = iCell-1
       jNode(4) = jCell

    elseif (trim(face)=='east') then

       !     4-----3       z
       !     |     |       ^
       !     |     |       |
       !     1-----2       ---> y

       x(1) = yvertex(iCell, jCell-1)
       x(2) = yvertex(iCell, jCell)
       x(3) = x(2)
       x(4) = x(1)

       s(1) = stagusrf(iCell, jCell-1)
       s(2) = stagusrf(iCell, jCell)
       s(3) = s(2)
       s(4) = s(1)

       h(1) = stagthck(iCell, jCell-1)
       h(2) = stagthck(iCell, jCell)
       h(3) = h(2)
       h(4) = h(1)

       iNode(1) = iCell
       jNode(1) = jCell-1

       iNode(2) = iCell
       jNode(2) = jCell

       iNode(3) = iCell
       jNode(3) = jCell

       iNode(4) = iCell
       jNode(4) = jCell-1

    elseif (trim(face)=='south') then

       !     4-----3       z
       !     |     |       ^
       !     |     |       |
       !     1-----2       ---> x

       x(1) = xvertex(iCell-1, jCell-1)
       x(2) = xvertex(iCell,   jCell-1)
       x(3) = x(2)
       x(4) = x(1)

       s(1) = stagusrf(iCell-1, jCell-1)
       s(2) = stagusrf(iCell,   jCell-1)
       s(3) = s(2)
       s(4) = s(1)
       
       h(1) = stagthck(iCell-1, jCell-1)
       h(2) = stagthck(iCell,   jCell-1)
       h(3) = h(2)
       h(4) = h(1)
    
       iNode(1) = iCell-1
       jNode(1) = jCell-1

       iNode(2) = iCell
       jNode(2) = jCell-1

       iNode(3) = iCell
       jNode(3) = jCell-1

       iNode(4) = iCell-1
       jNode(4) = jCell-1

    elseif (trim(face)=='north') then

       !     4-----3       z
       !     |     |       ^
       !     |     |       |
       !     1-----2       ---> -x

       x(1) = xvertex(iCell,   jCell)
       x(2) = xvertex(iCell-1, jCell)
       x(3) = x(2)
       x(4) = x(1)

       s(1) = stagusrf(iCell,   jCell)
       s(2) = stagusrf(iCell-1, jCell)
       s(3) = s(2)
       s(4) = s(1)

       h(1) = stagthck(iCell,   jCell)
       h(2) = stagthck(iCell-1, jCell)
       h(3) = h(2)
       h(4) = h(1)

       iNode(1) = iCell
       jNode(1) = jCell

       iNode(2) = iCell-1
       jNode(2) = jCell

       iNode(3) = iCell-1
       jNode(3) = jCell

       iNode(4) = iCell
       jNode(4) = jCell

    endif

    ! loop over element faces in column
    ! assume k increases from upper surface to bottom 

    do k = 1, nz-1                          

       ! Compute the local y coordinate (i.e., the actual z coordinate)
       y(1) = s(1) - sigma(k+1)*h(1)   ! lower left
       y(2) = s(2) - sigma(k+1)*h(2)   ! lower right
       y(3) = s(3) - sigma(k)  *h(3)   ! upper right
       y(4) = s(4) - sigma(k)  *h(4)   ! upper left

       ! Set the k index for each node
       kNode(1) = k+1
       kNode(2) = k+1
       kNode(3) = k
       kNode(4) = k

       ! loop over quadrature points

       do p = 1, nQuadPoints_2d  

          ! Compute basis function derivatives and det(J) for this quadrature point
          ! For now, pass in i, j, k, p for debugging
          !TODO - We don't actually need the derivatives, just detJ.  
          !       Should we modify the subroutine so that derivatives are optional?

          call get_basis_function_derivatives_2d(nNodesPerElement_2d,                  &
                                                 x(:),       y(:),                     & 
                                                 dphi_dxr_2d(:,p), dphi_dyr_2d(:,p),   &   
                                                 dphi_dx_2d(:),    dphi_dy_2d(:),      &
                                                 detJ, iCell, jCell, k, p)
          
          ! Evaluate the quantity dz = H - s at this quadrature point
 
          dz = 0.d0
          do n = 1, nNodesPerElement_2d
             dz = dz + phi_2d(n,p) * (h(n) - s(n))
          enddo

          if (verbose .and. this_rank==rtest .and. iCell==itest .and. jCell==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Increment shelf load vector, i, j, face, k, p =', iCell, jCell, trim(face), k, p
             print*, 'dz =', dz
             print*, 'detJ/vol0 =', detJ/vol0
             print*, 'detJ/vol0 * dz =', detJ/vol0 * dz
          endif

          ! Increment the load vector with the shelf water pressure contribution from 
          !  this quadrature point.
          ! Increment bu for east/west faces and bv for north/south faces.
          !TODO - Could eliminate the if's by using multiplying factors of +1, -1, or 0
          !       depending on the face.

          if (trim(face) == 'west') then

             do n = 1, nNodesPerElement_2d
                bu(kNode(n),iNode(n),jNode(n)) = bu(kNode(n),iNode(n),jNode(n))    &
                                               + 0.5d0*rhoi*grav*dz * wqp(p) * detJ/vol0 * phi(n,p)
             enddo

             if (verbose .and. this_rank==rtest .and. iCell==itest .and. jCell==jtest .and. k==ktest .and. p==ptest) then
                print*, 'n, phi(n), delta(bu):', n, phi(n,p), &
                         0.5*rhoi*grav*wqp(p)*detJ/vol0 * dz * phi(n,p)
             endif

          elseif (trim(face) == 'east') then 

             do n = 1, nNodesPerElement_2d
                bu(kNode(n),iNode(n),jNode(n)) = bu(kNode(n),iNode(n),jNode(n))    &
                                               - 0.5d0*rhoi*grav*dz * wqp(p) * detJ/vol0 * phi(n,p)
             enddo

             if (verbose .and. this_rank==rtest .and. iCell==itest .and. jCell==jtest .and. k==ktest .and. p==ptest) then
                print*, 'n, phi(n), delta(bu):', n, phi(n,p), &
                         -0.5*rhoi*grav*wqp(p)*detJ/vol0 * dz * phi(n,p)
             endif

          elseif (trim(face) == 'south') then 

             do n = 1, nNodesPerElement_2d
                bv(kNode(n),iNode(n),jNode(n)) = bv(kNode(n),iNode(n),jNode(n))    &
                                               + 0.5d0*rhoi*grav*dz * wqp(p) * detJ/vol0 * phi(n,p)
             enddo

             if (verbose .and. this_rank==rtest .and. iCell==itest .and. jCell==jtest .and. k==ktest .and. p==ptest) then
                print*, 'n, phi(n), delta(bv):', n, phi(n,p), &
                         0.5*rhoi*grav*wqp(p)*detJ/vol0 * dz * phi(n,p)
             endif

          elseif (trim(face) == 'north') then 
 
             do n = 1, nNodesPerElement_2d
                bv(kNode(n),iNode(n),jNode(n)) = bv(kNode(n),iNode(n),jNode(n))    &
                                               - 0.5d0*rhoi*grav*dz * wqp(p) * detJ/vol0 * phi(n,p)
             enddo

             if (verbose .and. this_rank==rtest .and. iCell==itest .and. jCell==jtest .and. k==ktest .and. p==ptest) then
                print*, 'n, phi(n), delta(bv):', n, phi(n,p), &
                         -0.5*rhoi*grav*wqp(p)*detJ/vol0 * dz * phi(n,p)
             endif

          endif   ! face = N/S/E/W

       enddo   ! nQuadPoints_2d

    enddo   ! k

  end subroutine lateral_shelf_bc

!****************************************************************************

!  Note: This subroutine has to be called once per nonlinear iteration if
!         we're iterating on the effective viscosity.
!        Currently we recompute some geometric info that does not
!         change between nonlinear iterations.  Should we compute this info
!         once and store it?

  subroutine assemble_stiffness_matrix(nx,               ny,              &
                                       nz,               nhalo,           &
                                       sigma,            active_cell,     &
                                       xVertex,          yVertex,         &
                                       uvel,             vvel,            &
                                       stagusrf,         stagthck,        &
                                       flwafact,         whichefvs,       &
                                       Auu,              Auv,             &
                                       Avu,              Avv,             &
                                       sia_factor,       ssa_factor)

    ! Assemble the stiffness matrix A and load vector b in the linear system Ax = b
 
    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical layers at which velocity is computed
                                     ! Note: the number of elements per column is nz-1
       nhalo                         ! number of halo layers

    real(dp), dimension(nz), intent(in) ::    &
       sigma                         ! sigma vertical coordinate

    logical, dimension(nx,ny), intent(in) ::  &
       active_cell                   ! true if cell contains ice and borders a locally owned vertex

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xVertex, yVertex     ! x and y coordinates of vertices

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       uvel, vvel         ! velocity components (m/s)

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf,       &  ! upper surface elevation on staggered grid (m)
       stagthck           ! ice thickness on staggered grid (m)

    real(dp), dimension(nz-1,nx,ny), intent(in) ::  &
       flwafact           ! temperature-based flow factor, 0.5 * A^(-1/n), 
                          ! used to compute the effective viscosity
                          ! units: Pa s^(1/n)

    integer, intent(in) :: whichefvs      ! option for effective viscosity calculation 
                                          ! (calculate it or make it uniform)

    real(dp), dimension(27,nz,nx-1,ny-1), intent(out) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv                                    

    real(dp), intent(in) ::  &
       sia_factor,      & ! = 1. if SIA terms are included, else = 0.
       ssa_factor         ! = 1. if SSA terms are included, else = 0.

    !---------------------------------------------------------
    ! Local variables
    !---------------------------------------------------------

    real(dp) ::   &
       detJ               ! determinant of J

    real(dp), dimension(nNodesPerElement) ::   &
       dphi_dx, dphi_dy, dphi_dz   ! derivatives of basis function, evaluated at quad pts

    !----------------------------------------------------------------
    ! Note: Kuu, Kuv, Kvu, and Kvv are 8x8 components of the stiffness matrix
    !       for the local element.  (The combined stiffness matrix is 16x16.)
    !
    ! Once these matrices are formed, their coefficients are summed into the assembled
    !  matrices Auu, Auv, Avu, Avv.  The A matrices each have as many rows as there are
    !  active nodes, but only 27 columns, corresponding to the 27 vertices that belong to
    !  the 8 elements sharing a given node.
    !
    ! The homegrown structured PCG solver works with the dense A matrices in the form
    ! computed here.  For a SLAP solver, the terms of the A matrices are put
    ! in a sparse matrix during preprocessing.
    !
    ! For a Trilinos solve, the terms of the element matrices should be sent to
    ! Trilinos a row at a time.  This capability has not been implemented as of Jan. 2013.
    !----------------------------------------------------------------

    real(dp), dimension(nNodesPerElement, nNodesPerElement) ::   &   !
       Kuu,          &    ! element stiffness matrix, divided into 4 parts as shown below
       Kuv,          &    !  
       Kvu,          &    !
       Kvv                !    Kuu  | Kuv
                          !    _____|____          
                          !         |
                          !    Kvu  | Kvv
                          !         
                          ! Kvu may not be needed if matrix is symmetric, but is included for now

    real(dp), dimension(nNodesPerElement) ::     &
       x, y, z,         & ! Cartesian coordinates of nodes
       u, v,            & ! u and v at nodes
       s                  ! upper surface elevation at nodes

    real(dp) ::    &
       efvs               ! effective viscosity

    integer :: i, j, k, n, p
    integer :: iA, jA, kA

    integer :: iNode, jNode, kNode

    if (verbose .and. main_task) then
       print*, ' '
       print*, 'In assemble_stiffness_matrix'
    endif

    ! Initialize global stiffness matrix and load vector
    !TODO - Initialize at higher level and pass as inout?
    !       Compute load vector in different subroutine?

    Auu(:,:,:,:) = 0.d0
    Auv(:,:,:,:) = 0.d0
    Avu(:,:,:,:) = 0.d0
    Avv(:,:,:,:) = 0.d0

    ! Sum over elements in active cells 
    ! Loop over all cells that border locally owned vertices.

    do j = nhalo+1, ny-nhalo+1
    do i = nhalo+1, nx-nhalo+1
       
     if (active_cell(i,j)) then

       do k = 1, nz-1    ! loop over elements in this column 
                         ! assume k increases from upper surface to bed

          if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, 'i, j, k:', i, j, k
             print*, ' '
          endif

          ! Initialize element stiffness matrix
          Kuu(:,:) = 0.d0
          Kuv(:,:) = 0.d0
          Kvu(:,:) = 0.d0
          Kvv(:,:) = 0.d0
    
          ! compute spatial coordinates, velocity, and upper surface elevation for each node

          do n = 1, nNodesPerElement

             ! Determine (k,i,j) for this node
             ! The reason for the '7' is that node 7, in the NE corner of the upper layer, has index (k,i,j).
             ! Indices for other nodes are computed relative to this node.
             iNode = i + ishift(7,n)
             jNode = j + jshift(7,n)
             kNode = k + kshift(7,n)

             x(n) = xVertex(iNode,jNode)
             y(n) = yVertex(iNode,jNode)
             z(n) = stagusrf(iNode,jNode) - sigma(kNode)*stagthck(iNode,jNode)
             u(n) = uvel(kNode,iNode,jNode)
             v(n) = vvel(kNode,iNode,jNode)
             s(n) = stagusrf(iNode,jNode)

             if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
                print*, 'i, j, k, n, x, y, z:', i, j, k, n, x(n), y(n), z(n)
                print*, 's, u, v:', s(n), u(n), v(n)
             endif

          enddo   ! nodes per element

          ! Loop over quadrature points for this element
   
          !TODO - Think about how often to compute things.  Geometric quantities need to be
          !       computed only once per nonlinear solve, whereas efvs must be recomputed
          !       after each linear solve.

          do p = 1, nQuadPoints

             ! Evaluate the derivatives of the element basis functions
             ! at this quadrature point.

!WHL - debug - Pass in i, j, k, and p for now

             call get_basis_function_derivatives(nNodesPerElement, &
                                                 x(:),          y(:),          z(:),           &
                                                 dphi_dxr(:,p), dphi_dyr(:,p), dphi_dzr(:,p),  &
                                                 dphi_dx(:),    dphi_dy(:),    dphi_dz(:),     &
                                                 detJ , i, j, k, p                      )

             if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
!!                print*, ' '
!!                print*, 'Derivatives of basis functions, p =', p
                do n = 1, nNodesPerElement
!!                   print*, 'n, dphi_dx, dphi_dy, dphi_dz:', n, dphi_dx(n), dphi_dy(n), dphi_dz(n)
                enddo
             endif


!WHL - debug - Pass in i, j, k, and p for now

             call compute_effective_viscosity(whichefvs,        nNodesPerElement,             &
                                              dphi_dx(:),       dphi_dy(:),     dphi_dz(:),   &
                                              u(:),             v(:),                         & 
                                              flwafact(k,i,j),  efvs,                         &
                                              sia_factor,       ssa_factor,  i, j, k, p )

             if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
!                print*, ' '
!                print*, 'efvs:', efvs
             endif

             ! Increment the element stiffness matrix with the contribution from
             ! this quadrature point.

!WHL - debug - Pass in i, j, k, and p for now

             call element_matrix_blatter_pattyn(nNodesPerElement,                          &
                                                wqp(p),           detJ,        efvs,       &
                                                dphi_dx(:),       dphi_dy(:),  dphi_dz(:), &
                                                Kuu(:,:),         Kuv(:,:),                &
                                                Kvu(:,:),         Kvv(:,:),                &
                                                sia_factor,       ssa_factor,  i, j, k, p )

          enddo   ! nQuadPoints

          if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Kuu:'
             do n = 1, nNodesPerElement
                write(6, '(8e12.4)') Kuu(n,:)
             enddo
             print*, 'Kuv:'
             do n = 1, nNodesPerElement
                write(6, '(8e12.4)') Kuv(n,:)
             enddo
             print*, 'Kvu:'
             do n = 1, nNodesPerElement
                write(6, '(8e12.4)') Kvu(n,:)
             enddo
             print*, 'Kvv:'
             do n = 1, nNodesPerElement
                write(6, '(8e12.4)') Kvv(n,:)
             enddo
          endif

          if (check_symmetry) then

             call check_symmetry_element_matrix(nNodesPerElement,  &
                                                Kuu, Kuv, Kvu, Kvv)

          endif

         ! If solving in parallel with Trilinos, we have the option at this point to 
         ! call a sum_into_global_matrix routine, passing one row at a time of the
         ! element matrix.  Trilinos will handle the rest.
         !

          if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Insert Kuu into Auu'
          endif

          call element_to_global_matrix(nx,           ny,          nz,          &
                                        nNodesPerElement,                       &
                                        i,            j,           k,           &
                                        Kuu,          Auu)

          if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Insert Kuv into Auv'
          endif

          call element_to_global_matrix(nx,           ny,          nz,          &
                                        nNodesPerElement,                       &
                                        i,            j,           k,           &
                                        Kuv,          Auv)

          if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Insert Kvu into Avu'
          endif

          call element_to_global_matrix(nx,           ny,          nz,          &
                                        nNodesPerElement,                       &
                                        i,            j,           k,           &
                                        Kvu,          Avu)

          if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Insert Kvv into Avv'
          endif

          call element_to_global_matrix(nx,           ny,          nz,          &
                                        nNodesPerElement,                       &
                                        i,            j,           k,           &
                                        Kvv,          Avv)

       enddo   ! nz  (loop over elements in this column)

     endif   ! active cell

    enddo      ! i
    enddo      ! j


  end subroutine assemble_stiffness_matrix

!****************************************************************************

!WHL - debug - pass in i, j, k, p for now
  subroutine get_basis_function_derivatives(nNodesPerElement, &
                                            xNode,       yNode,     zNode,    &
                                            dphi_dxr,    dphi_dyr,  dphi_dzr, &
                                            dphi_dx,     dphi_dy,   dphi_dz,  &
                                            detJ,        i, j, k, p)

!TODO - Change subroutine name?

    !------------------------------------------------------------------
    ! Evaluate the x, y and z derivatives of the element basis functions
    ! at a particular quadrature point.
    !
    ! Also determine the Jacobian of the transformation between the
    ! reference element and the true element.
    ! 
    ! This subroutine should work for any 3D element with any number of nodes.
    !------------------------------------------------------------------

    integer, intent(in) :: nNodesPerElement   ! number of nodes per element
 
    real(dp), dimension(nNodesPerElement), intent(in) :: &
                xNode, yNode, zNode,          &! nodal coordinates
                dphi_dxr, dphi_dyr, dphi_dzr   ! derivatives of basis functions at quad pt
                                               !  wrt x, y and z in reference element

    real(dp), dimension(nNodesPerElement), intent(out) :: &
                dphi_dx , dphi_dy, dphi_dz     ! derivatives of basis functions at quad pt
                                               !  wrt x, y and z in true Cartesian coordinates  

    real(dp), intent(out) :: &
                detJ      ! determinant of Jacobian matrix

    real(dp), dimension(3,3) ::  &
                Jac,      &! Jacobian matrix
                Jinv,     &! inverse Jacobian matrix
                cofactor   ! matrix of cofactors

    integer, intent(in) :: i, j, k, p

    integer :: n, row, col

!WHL - debug
    real(dp), dimension(3,3) :: prod     ! Jac * Jinv (should be identity matrix)

    !------------------------------------------------------------------
    ! Compute the Jacobian for the transformation from the reference
    ! coordinates to the true coordinates:
    !
    !                 |                                                                          |
    !                 | sum_n{dphi_n/dxr * xn}   sum_n{dphi_n/dxr * yn}   sum_n{dphi_n/dxr * zn} |
    !   J(xr,yr,zr) = |                                                                          |
    !                 | sum_n{dphi_n/dyr * xn}   sum_n{dphi_n/dyr * yn}   sum_n{dphi_n/dyr * zn} |
    !                 |                                                                          |
    !                 | sum_n{dphi_n/dzr * xn}   sum_n{dphi_n/dzr * yn}   sum_n{dphi_n/dzr * zn} |
    !                 !                                                                          |
    !
    ! where (xn,yn,zn) are the true Cartesian nodal coordinates,
    !       (xr,yr,zr) are the coordinates of the quad point in the reference element,
    !       and sum_n denotes a sum over nodes.
    !------------------------------------------------------------------

    Jac(:,:) = 0.d0

    if (verbose_Jac .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, 'In get_basis_function_derivatives: i, j, k, p =', i, j, k, p
    endif

    do n = 1, nNodesPerElement
       if (verbose_Jac .and. i==itest .and. j==jtest .and. k==ktest) then
          print*, ' '
          print*, 'n, x, y, z:', n, xNode(n), yNode(n), zNode(n)
          print*, 'dphi_dxyz:', dphi_dxr(n), dphi_dyr(n), dphi_dzr(n) 
       endif
       Jac(1,1) = Jac(1,1) + dphi_dxr(n) * xNode(n)
       Jac(1,2) = Jac(1,2) + dphi_dxr(n) * yNode(n)
       Jac(1,3) = Jac(1,3) + dphi_dxr(n) * zNode(n)
       Jac(2,1) = Jac(2,1) + dphi_dyr(n) * xNode(n)
       Jac(2,2) = Jac(2,2) + dphi_dyr(n) * yNode(n)
       Jac(2,3) = Jac(2,3) + dphi_dyr(n) * zNode(n)
       Jac(3,1) = Jac(3,1) + dphi_dzr(n) * xNode(n)
       Jac(3,2) = Jac(3,2) + dphi_dzr(n) * yNode(n)
       Jac(3,3) = Jac(3,3) + dphi_dzr(n) * zNode(n)
    enddo

    !------------------------------------------------------------------
    ! Compute the determinant and inverse of J
    !------------------------------------------------------------------

    cofactor(1,1) =   Jac(2,2)*Jac(3,3) - Jac(2,3)*Jac(3,2)
    cofactor(1,2) = -(Jac(2,1)*Jac(3,3) - Jac(2,3)*Jac(3,1))
    cofactor(1,3) =   Jac(2,1)*Jac(3,2) - Jac(2,2)*Jac(3,1)
    cofactor(2,1) = -(Jac(1,2)*Jac(3,3) - Jac(1,3)*Jac(3,2))
    cofactor(2,2) =   Jac(1,1)*Jac(3,3) - Jac(1,3)*Jac(3,1)
    cofactor(2,3) = -(Jac(1,1)*Jac(3,2) - Jac(1,2)*Jac(3,1))
    cofactor(3,1) =   Jac(1,2)*Jac(2,3) - Jac(1,3)*Jac(2,2)
    cofactor(3,2) = -(Jac(1,1)*Jac(2,3) - Jac(1,3)*Jac(2,1))
    cofactor(3,3) =   Jac(1,1)*Jac(2,2) - Jac(1,2)*Jac(2,1)

    detJ = Jac(1,1)*cofactor(1,1) + Jac(1,2)*cofactor(1,2) + Jac(1,3)*cofactor(1,3)

    if (verbose_Jac .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, ' '
       print*, 'detJ1:', Jac(1,1)*cofactor(1,1) + Jac(1,2)*cofactor(1,2) + Jac(1,3)*cofactor(1,3)
       print*, 'detJ2:', Jac(2,1)*cofactor(2,1) + Jac(2,2)*cofactor(2,2) + Jac(2,3)*cofactor(2,3)
       print*, 'detJ3:', Jac(3,1)*cofactor(3,1) + Jac(3,2)*cofactor(3,2) + Jac(3,3)*cofactor(3,3)
    endif

    if (abs(detJ) > 0.d0) then
       do col = 1, 3
          do row = 1, 3
             Jinv(row,col) = cofactor(col,row)
          enddo
       enddo
       Jinv(:,:) = Jinv(:,:) / detJ
    else
       !WHL - do a proper abort here
       print*, 'stopping, det J = 0'
       print*, 'i, j, k, p:', i, j, k, p
       print*, 'Jacobian matrix:'
       print*, Jac(1,:)
       print*, Jac(2,:)
       print*, Jac(3,:)     
       stop
    endif

    if (verbose_Jac .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, ' '
       print*, 'Jacobian calc, p =', p
       print*, 'det J =', detJ
       print*, ' '
       print*, 'Jacobian matrix:'
       print*, Jac(1,:)
       print*, Jac(2,:)
       print*, Jac(3,:)
       print*, ' '
       print*, 'cofactor matrix:'
       print*, cofactor(1,:)
       print*, cofactor(2,:)
       print*, cofactor(3,:)
       print*, ' '
       print*, 'Inverse matrix:'
       print*, Jinv(1,:)
       print*, Jinv(2,:)
       print*, Jinv(3,:)
       print*, ' '
       prod = matmul(Jac, Jinv)
       print*, 'Jac*Jinv:'
       print*, prod(1,:)
       print*, prod(2,:)
       print*, prod(3,:)
    endif

    ! bug check - Verify that J * Jinv = I

    prod = matmul(Jac,Jinv)
    do col = 1, 3
       do row = 1, 3
          if (abs(prod(row,col) - identity3(row,col)) > 1.d-12) then
             !TODO - do a proper abort here 
             print*, 'stopping, Jac * Jinv /= identity'
             print*, 'i, j, k, p:', i, j, k, p
             print*, 'Jac*Jinv:'
             print*, prod(1,:)
             print*, prod(2,:)
             print*, prod(3,:)
             stop
          endif
       enddo
    enddo

    !------------------------------------------------------------------
    ! Compute the contribution of this quadrature point to dphi/dx and dphi/dy
    ! for each basis function.
    !
    !   | dphi_n/dx |          | dphi_n/dxr |
    !   |           |          |            | 
    !   | dphi_n/dy | = Jinv * | dphi_n/dyr |
    !   |           |          |            |
    !   | dphi_n/dz |          | dphi_n/dzr |
    !
    !------------------------------------------------------------------

    dphi_dx(:) = 0.d0
    dphi_dy(:) = 0.d0
    dphi_dz(:) = 0.d0

    do n = 1, nNodesPerElement
       dphi_dx(n) = dphi_dx(n) + Jinv(1,1)*dphi_dxr(n)  &
                               + Jinv(1,2)*dphi_dyr(n)  &
                               + Jinv(1,3)*dphi_dzr(n)
       dphi_dy(n) = dphi_dy(n) + Jinv(2,1)*dphi_dxr(n)  &
                               + Jinv(2,2)*dphi_dyr(n)  &
                               + Jinv(2,3)*dphi_dzr(n)
       dphi_dz(n) = dphi_dz(n) + Jinv(3,1)*dphi_dxr(n)  &
                               + Jinv(3,2)*dphi_dyr(n)  &
                               + Jinv(3,3)*dphi_dzr(n)
    enddo

    ! debug - Check that the sum of dphi_dx, etc. is close to zero  
    !TODO - Think about why this should be the case

    if (abs(sum(dphi_dx)) > 1.d-12) then
        print*, 'stopping, sum over basis functions of dphi_dx > 0'
        print*, 'dphi_dx =', dphi_dx(:)
        print*, 'i, j, k, p =', i, j, k, p
        stop
    endif

    if (abs(sum(dphi_dy)) > 1.d-12) then
        print*, 'stopping, sum over basis functions of dphi_dy > 0'
        print*, 'dphi_dy =', dphi_dy(:)
        print*, 'i, j, k, p =', i, j, k, p
        stop
    endif

    if (abs(sum(dphi_dz)) > 1.d-12) then
        print*, 'stopping, sum over basis functions of dphi_dz > 0'
        print*, 'dphi_dz =', dphi_dz(:)
        print*, 'i, j, k, p =', i, j, k, p
        stop
    endif

  end subroutine get_basis_function_derivatives

!****************************************************************************

!WHL - Pass i, j, k, and p for now

  subroutine compute_effective_viscosity (whichefvs,   nNodesPerElement,       &
                                          dphi_dx,     dphi_dy,  dphi_dz,      &
                                          uvel,        vvel,                   &
                                          flwafact,    efvs,                   &
                                          sia_factor,  ssa_factor,  i, j, k, p )

    ! Compute effective viscosity at a quadrature point, based on the latest
    ! guess for the velocity field

    integer, intent(in) :: i, j, k, p

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) :: &
       whichefvs       ! method for computing effective viscosity
                       ! 0 = based on effective strain rate
                       ! 1 = constant value
                       !TODO - Let 0 = constant value, 1 = eff str rate?
                       ! Add an option representing a true constant value, 
                       ! rather than temperature-based?

    integer, intent(in) :: nNodesPerElement   ! number of nodes per element

    real(dp), dimension(nNodesPerElement), intent(in) ::  &
       dphi_dx, dphi_dy, dphi_dz         ! derivatives of basis functions,
                                         ! evaluated at this quadrature point

    real(dp), dimension(nNodesPerElement), intent(in) ::  &
       uvel, vvel      ! current guess for velocity at each node of element (m/s)

    real(dp), intent(in) ::  &
       flwafact        ! temperature-based flow factor for this element, 0.5 * A^{-1/n}
                       ! units: Pa s^{1/n}

    real(dp), intent(out) ::   &
       efvs            ! effective viscosity at this quadrature point (Pa s)
                       ! computed as 0.5 * A^{-1/n) * effstrain^{(1-n)/n)}
                       
    real(dp), intent(in) ::  &
       sia_factor,      & ! = 1. if SIA terms are included, else = 0.
       ssa_factor         ! = 1. if SSA terms are included, else = 0.

    !----------------------------------------------------------------
    ! Local parameters
    !----------------------------------------------------------------

    real(dp), parameter ::   &
       effstrain_min = 1.d-20,          &! minimum value of effective strain rate, s^{-1}
                                         ! GLAM uses 1.d-20 s^{-1} for minimum effective strain rate
       p_effstr = (1.d0 - real(gn,dp)) / real(gn,dp)    ! exponent (1-n)/n in effective viscosity relation
                                                               
    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    real(dp) ::            &
       du_dx, dv_dx,       & ! strain rate components
       du_dy, dv_dy,       &
       du_dz, dv_dz,       &
       effstrain,          & ! effective strain rate, s^{-1}
       effstrainsq           ! square of effective strain rate
        
    integer :: n

!WHL - debug
    if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, ' '
       print*, 'Compute efvs, i, j, k, p =', i, j, k, p
    endif

    select case(whichefvs)

!!    case(2)
    case(HO_EFVS_CONSTANT)

       efvs = 1.d7 * scyr   ! Steve Price recommends 10^6 to 10^7 Pa yr

!WHL - This is the glam-type scaling
!!       efvs = efvs * scyr/tim0 / tau0   ! tau0 = rhoi*grav*thk0

       if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
          print*, 'Set efvs = constant:', efvs
       endif

!!    case(1)
    case(HO_EFVS_FLOWFACT)      ! set the effective viscosity to a multiple of the flow factor, 0.5*A^(-1/n)
                 
                 !SCALING: Set the effective strain rate (s^{-1}) based on typical 
                 !         velocity and length scales, to agree with GLAM.
  
       effstrain = vel_scale/len_scale   ! typical strain rate, s^{-1}  
       efvs = flwafact * effstrain**p_effstr  

       if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
          print*, 'flwafact, effstrain, efvs =', flwafact, effstrain, efvs       
       endif

!!    case(0)
    case(HO_EFVS_NONLINEAR)    ! calculate effective viscosity based on effective strain rate, n = 3

       du_dx = 0.d0
       dv_dx = 0.d0
       du_dy = 0.d0
       dv_dy = 0.d0
       du_dz = 0.d0
       dv_dz = 0.d0

       ! Compute strain rate components at this quadrature point
       ! by summing over basis functions

       do n = 1, nNodesPerElement

          du_dx = du_dx + dphi_dx(n)*uvel(n)
          dv_dx = dv_dx + dphi_dx(n)*vvel(n)

          du_dy = du_dy + dphi_dy(n)*uvel(n)
          dv_dy = dv_dy + dphi_dy(n)*vvel(n)

          du_dz = du_dz + dphi_dz(n)*uvel(n)
          dv_dz = dv_dz + dphi_dz(n)*vvel(n)

       enddo   ! nNodesPerElement

       ! Compute effective strain rate at this quadrature point (PGB 2012, eq. 3 and 9)
       ! Units are s^(-1)

       effstrainsq = effstrain_min**2                                      &
                   + ssa_factor * (du_dx**2 + dv_dy**2 + du_dx*dv_dy       &
                                   + 0.25d0 * (dv_dx + du_dy)**2)          &
                   + sia_factor * 0.25d0 * (du_dz**2 + dv_dz**2)

       effstrain = sqrt(effstrainsq)

       ! Compute effective viscosity (PGB 2012, eq. 4)
       ! Units: flwafact has units Pa s^{1/n}
       !        effstrain has units s^{-1}
       !        p_effstr = (1-n)/n 
       !                 = -2/3 for n=3
       ! Thus efvs has units Pa s
 
       efvs = flwafact * effstrain**p_effstr

       if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
!!       if (verbose .and. i==itest .and. j==jtest+1 .and. k==ktest) then
          print*, 'effstrain_min =', effstrain_min
          print*, 'du/dx, du/dy, du/dz =', du_dx, du_dy, du_dz
          print*, 'dv/dx, dv/dy, dv/dz =', dv_dx, dv_dy, dv_dz
          print*, 'flwafact, effstrain (s-1), efvs (Pa s) =', flwafact, effstrain, efvs
       endif

    end select

  end subroutine compute_effective_viscosity

!****************************************************************************

!WHL - debug - Pass in i, j, k, and p for now

  subroutine element_matrix_blatter_pattyn(nNodesPerElement,                &
                                           wqp,         detJ,     efvs,     &
                                           dphi_dx,     dphi_dy,  dphi_dz,  &
                                           Kuu,         Kuv,                &
                                           Kvu,         Kvv,                &
                                           sia_factor,  ssa_factor,  ii, jj, k, p)

    !------------------------------------------------------------------
    ! Increment the stiffness matrices Kuu, Kuv, Kvu, Kvv with the
    ! contribution from a particular quadrature point, 
    ! based on the Blatter-Pattyn first-order equations.
    !
    ! This subroutine should work for any 3D element with any number of nodes.
    !------------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

!WHL - debug - Pass in i, j, k, and p for now
    integer, intent(in) :: ii, jj, k, p

    integer, intent(in) :: nNodesPerElement  ! number of nodes per element

    real(dp), intent(in) ::    &
             wqp,        &! weight for this quadrature point
             detJ,       &! determinant of Jacobian for the transformation
                          !  between the reference element and true element
             efvs         ! effective viscosity at this quadrature point

    real(dp), dimension(nNodesPerElement), intent(in) ::  &
             dphi_dx, dphi_dy, dphi_dz   ! derivatives of basis functions,
                                         ! evaluated at this quadrature point

    real(dp), dimension(nNodesPerElement,nNodesPerElement), intent(inout) :: &
             Kuu, Kuv, Kvu, Kvv     ! components of element stiffness matrix

    real(dp), intent(in) ::  &
       sia_factor,      & ! = 1. if SIA terms are included, else = 0.
       ssa_factor         ! = 1. if SSA terms are included, else = 0.

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j

!WHL - debug
    if (verbose .and. this_rank==rtest .and. ii==itest .and. jj==jtest .and. k==ktest) then
       print*, ' '
       print*, 'Increment element matrix, p =', p
    endif

    ! Increment the element stiffness matrices for the first-order Blatter-Pattyn equations
    ! The terms in parentheses can be derived from PGB 2012, eq. 13 and 15.
    ! The factor of 2 in front of efvs has been absorbed into the quantities in parentheses.
    
    do j = 1, nNodesPerElement      ! columns of K
       do i = 1, nNodesPerElement   ! rows of K

!       if (verbose .and. this_rank==rtest .and. ii==itest .and. jj==jtest .and. k==ktest .and. p==ptest) then
!          print*, 'efvs, wqp, detJ/vol0 =', efvs, wqp, detJ/vol0
!          print*, 'dphi_dz(1) =', dphi_dz(1)
!          print*, 'dphi_dx(1) =', dphi_dx(1)
!          print*, 'Kuu dphi/dz increment(1,1) =', efvs*wqp*detJ/vol0*dphi_dz(1)*dphi_dz(1)
!          print*, 'Kuu dphi/dx increment(1,1) =', efvs*wqp*detJ/vol0*4.d0*dphi_dx(1)*dphi_dx(1)
!       endif

          !WHL - Note volume scaling such that detJ/vol0 is closer to unity

          Kuu(i,j) = Kuu(i,j) + efvs*wqp*detJ/vol0 *                                               &
                              ( ssa_factor * (4.d0*dphi_dx(j)*dphi_dx(i) + dphi_dy(j)*dphi_dy(i))  &
                              + sia_factor * (dphi_dz(j)*dphi_dz(i)) )

          Kuv(i,j) = Kuv(i,j) + efvs*wqp*detJ/vol0 *                                               &
                                ssa_factor * (2.d0*dphi_dx(j)*dphi_dy(i) + dphi_dy(j)*dphi_dx(i))

          Kvu(i,j) = Kvu(i,j) + efvs*wqp*detJ/vol0 *                                               &
                                ssa_factor * (2.d0*dphi_dy(j)*dphi_dx(i) + dphi_dx(j)*dphi_dy(i))

          Kvv(i,j) = Kvv(i,j) + efvs*wqp*detJ/vol0 *                                               &
                              ( ssa_factor * (4.d0*dphi_dy(j)*dphi_dy(i) + dphi_dx(j)*dphi_dx(i))  &
                              + sia_factor * (dphi_dz(j)*dphi_dz(i)) )

       enddo  ! i (rows)
    enddo     ! j (columns)

  end subroutine element_matrix_blatter_pattyn

!****************************************************************************

  subroutine element_to_global_matrix(nx,           ny,          nz,          &
                                      nNodesPerElement,                       &
                                      iElement,     jElement,    kElement,    &
                                      Kmat,         Amat)

    ! Sum terms of element matrix K into dense assembled matrix A
    ! Here we assume that K is partitioned into Kuu, Kuv, Kvu, and Kvv,
    !  and similarly for A.
    ! So this subroutine is called four times per element.

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz                       ! number of vertical levels where velocity is computed

    integer, intent(in) ::   &
       nNodesPerElement         ! number of nodes per element

    integer, intent(in) ::   &
       iElement, jElement, kElement     ! i, j and k indices for this element

    real(dp), dimension(nNodesPerElement,nNodesPerElement), intent(in) ::  &
       Kmat              ! element matrix

    real(dp), dimension(27,nz,nx-1,ny-1), intent(inout) ::    &
       Amat              ! assembled matrix

    integer :: i, j, k, m
    integer :: iA, jA, kA
    integer :: n, nr, nc

!WHL - debug
    if (verbose .and. this_rank==rtest .and. iElement==itest .and. jElement==jtest .and. kElement==ktest) then
       print*, 'Element i, j, k:', iElement, jElement, kElement 
       print*, 'Rows of K:'
       do n = 1, nNodesPerElement
          write(6, '(8e12.4)') Kmat(n,:)
       enddo
    endif

!TODO - Switch loops or switch order of indices in K(nr,nc)?
!       Inner index nr is currently the outer loop

    do nr = 1, nNodesPerElement       ! rows of K

       ! Determine (k,i,j) for this node
       ! The reason for the '7' is that node 7, in the NE corner of the upper layer, has index (k,i,j).
       ! Indices for other nodes are computed relative to this node.
       i = iElement + ishift(7,nr)
       j = jElement + jshift(7,nr)
       k = kElement + kshift(7,nr)
      
       do nc = 1, nNodesPerElement    ! columns of K

          kA = kshift(nr,nc)           ! k index of A into which K(m,n) is summed
          iA = ishift(nr,nc)           ! similarly for i and j indices 
          jA = jshift(nr,nc)           ! these indices can take values -1, 0 and 1

          m = indxA(iA,jA,kA)
          Amat(m,k,i,j) = Amat(m,k,i,j) + Kmat(nr,nc)

!WHL - debug
!          if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest  &
!                      .and. iA == iAtest .and. jA==jAtest .and. kA==kAtest) then
!             print*, ' '
!             print*, 'i, j, k, iA, jA, kA:', i, j, k, iA, jA, kA
!             print*, 'nr, nc, Kmat, new Amat:', nr, nc, Kmat(nr,nc), Amat(kA,iA,jA,k,i,j)
!          endif

       enddo     ! nc

    enddo        ! nr

  end subroutine element_to_global_matrix

!****************************************************************************

  subroutine basal_sliding_bc(nx,               ny,              &
                              nz,               nhalo,           &
                              nNodesPerElement_2d,               &
                              active_cell,      beta,            &
                              xVertex,          yVertex,         &
                              Auu,              Avv)

    !------------------------------------------------------------------------
    ! Increment the Auu and Avv matrices with basal traction terms.
    !
    ! For now, assume a linear sliding law:
    !   tau_x = -beta*u
    !   tau_y = -beta*v
    ! 
    ! where beta is defined at vertices.
    !------------------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical layers at which velocity is computed
                                     ! Note: the number of elements per column is nz-1
       nhalo                         ! number of halo layers

    integer, intent(in) :: nNodesPerElement_2d

    logical, dimension(nx,ny), intent(in) ::  &
       active_cell                   ! true if cell contains ice and borders a locally owned vertex

    real(dp), dimension(nx-1,ny-1), intent(in) ::    &
       beta                          ! basal traction field

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xVertex, yVertex     ! x and y coordinates of vertices

    real(dp), dimension(27,nz,nx-1,ny-1), intent(inout) ::  &
       Auu, Avv             ! parts of stiffness matrix

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j, k, n, p

    integer :: iNode, jNode, kNode   ! global indices for nodes

    real(dp), dimension(nNodesPerElement_2d) ::   &
       x, y,        & ! Cartesian coordinates of nodes
       b              ! beta at nodes

    !TODO - These are not currently used except as dummy arguments
    real(dp), dimension(nNodesPerElement_2d) ::   &
       dphi_dx_2d, dphi_dy_2d    ! derivatives of basis functions, evaluated at quad pts

    real(dp) ::   &
       beta_qp,     & ! beta evaluated at quadrature point
       detJ           ! determinant of Jacobian for the transformation
                      !  between the reference element and true element

    real(dp), dimension(nNodesPerElement_2d,nNodesPerElement_2d) ::   &
       Kuu, Kvv       ! components of element matrix associated with basal sliding
        
    ! Sum over elements in active cells 
    ! Loop over cells that contain locally owned vertices
    !TODO - Make sure we don't step out of bounds
    !       With more care, we could skip some computations for vertices that are not locally owned.

    do j = nhalo+1, ny-nhalo+1
    do i = nhalo+1, nx-nhalo+1
       
!TODO - Should we exclude cells that have Dirichlet basal BCs for all vertices?

       if (active_cell(i,j)) then   ! ice is present

          ! Set x and y for each node

          !     4-----3       y
          !     |     |       ^
          !     |     |       |
          !     1-----2       ---> x

          x(1) = xVertex(i-1,j-1)
          x(2) = xVertex(i,j-1)
          x(3) = xVertex(i,j)
          x(4) = xVertex(i-1,j)

          y(1) = yVertex(i-1,j-1)
          y(2) = yVertex(i,j-1)
          y(3) = yVertex(i,j)
          y(4) = yVertex(i-1,j)

          b(1) = beta(i-1,j-1)          
          b(2) = beta(i,j-1)
          b(3) = beta(i,j)
          b(4) = beta(i-1,j)

          k = nz     ! basal layer only

          ! loop over quadrature points

          do p = 1, nQuadPoints_2d  

             ! Compute basis function derivatives and det(J) for this quadrature point
             ! For now, pass in i, j, k, p for debugging
             !TODO - We don't actually need the derivatives, just detJ.  
             !       Should we modify the subroutine so that derivatives are optional?

             call get_basis_function_derivatives_2d(nNodesPerElement_2d,                  &
                                                    x(:),       y(:),                     & 
                                                    dphi_dxr_2d(:,p), dphi_dyr_2d(:,p),   &   
                                                    dphi_dx_2d(:),    dphi_dy_2d(:),      &
                                                    detJ, i, j, k, p)
          
             ! Evaluate beta at this quadrature point
 
             beta_qp = 0.d0
             do n = 1, nNodesPerElement_2d
                beta_qp = beta_qp + phi_2d(n,p) * b(n)
             enddo

             if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'Increment basal traction, i, j, p =', i, j, p
                print*, 'beta_qp =', beta_qp
                print*, 'detJ/vol0 =', detJ/vol0
                print*, 'detJ/vol0 * beta_qp =', detJ/vol0 * beta_qp
             endif

!TODO - Test this subroutine

             call element_matrix_basal_sliding(nNodesPerElement_2d,    &
                                               wqp(p),      detJ,      & 
                                               beta_qp,     phi(:,p),  &
                                               Kuu(:,:),               &
                                               i, j, k, p)

             Kvv(:,:) = Kuu(:,:)   ! TODO: Is this true for more general sliding laws?

             ! Insert terms of element matrix into global matrices Auu and Avv

             call element_to_global_matrix_2d(nx,    ny,    nz,             &
                                              nNodesPerElement_2d,          &
                                              i,     j,    k,               &
                                              Kuu,   Auu)

             call element_to_global_matrix_2d(nx,    ny,    nz,             &
                                              nNodesPerElement_2d,          &
                                              i,     j,    k,               &
                                              Kvv,   Avv)

          enddo   ! nQuadPoints_2d

       endif      ! active_cell

    enddo         ! i
    enddo         ! j

  end subroutine basal_sliding_bc

!****************************************************************************
!WHL - Pass in i, j, k, p for now

  subroutine get_basis_function_derivatives_2d(nNodesPerElement_2d,     &
                                               xNode,       yNode,      &
                                               dphi_dxr,    dphi_dyr,   &
                                               dphi_dx,     dphi_dy,    &
                                               detJ, i, j, k, p)

!TODO - Change subroutine name?

    !------------------------------------------------------------------
    ! Evaluate the x and y derivatives of 2D element basis functions
    ! at a particular quadrature point.
    !
    ! Also determine the Jacobian of the transformation between the
    ! reference element and the true element.
    ! 
    ! This subroutine should work for any 2D element with any number of nodes.
    !------------------------------------------------------------------

    integer, intent(in) :: nNodesPerElement_2d   ! number of nodes per element
 
    real(dp), dimension(nNodesPerElement_2d), intent(in) :: &
                xNode, yNode,                   &! nodal coordinates
                dphi_dxr, dphi_dyr               ! derivatives of basis functions at quad pt
                                                 !  wrt x and y in reference element

    real(dp), dimension(nNodesPerElement_2d), intent(out) :: &
                dphi_dx , dphi_dy                ! derivatives of basis functions at quad pt
                                                 !  wrt x and y in true Cartesian coordinates  

    real(dp), intent(out) :: &
                detJ      ! determinant of Jacobian matrix

    real(dp), dimension(2,2) ::  &
                Jac,      &! Jacobian matrix
                Jinv       ! inverse Jacobian matrix

    integer, intent(in) :: i, j, k, p

    integer :: n, row, col

!WHL - debug
    real(dp), dimension(2,2) :: prod     ! Jac * Jinv (should be identity matrix)

    !------------------------------------------------------------------
    ! Compute the Jacobian for the transformation from the reference
    ! coordinates to the true coordinates:
    !
    !              |                                                  |
    !              | sum_n{dphi_n/dxr * xn}   sum_n{dphi_n/dxr * yn}  |
    !   J(xr,yr) = |                                                  |
    !              | sum_n{dphi_n/dyr * xn}   sum_n{dphi_n/dyr * yn}  |
    !              |                                                  |
    !
    ! where (xn,yn) are the true Cartesian nodal coordinates,
    !       (xr,yr) are the coordinates of the quad point in the reference element,
    !       and sum_n denotes a sum over nodes.
    !------------------------------------------------------------------

    Jac(:,:) = 0.d0

    if (verbose_Jac .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, 'In get_basis_function_derivatives_2d: i, j, k, p =', i, j, k, p
    endif

    do n = 1, nNodesPerElement_2d
       if (verbose_Jac .and. i==itest .and. j==jtest .and. k==ktest) then
          print*, ' '
          print*, 'n, x, y:', n, xNode(n), yNode(n)
          print*, 'dphi_dxr, dphi_dyr:', dphi_dxr(n), dphi_dyr(n) 
       endif
       Jac(1,1) = Jac(1,1) + dphi_dxr(n) * xNode(n)
       Jac(1,2) = Jac(1,2) + dphi_dxr(n) * yNode(n)
       Jac(2,1) = Jac(2,1) + dphi_dyr(n) * xNode(n)
       Jac(2,2) = Jac(2,2) + dphi_dyr(n) * yNode(n)
    enddo

    !------------------------------------------------------------------
    ! Compute the determinant and inverse of J
    !------------------------------------------------------------------

    detJ = Jac(1,1) * Jac(2,2)

    if (abs(detJ) > 0.d0) then
       Jinv(1,1) =  Jac(2,2)/detJ
       Jinv(1,2) = -Jac(1,2)/detJ
       Jinv(2,1) = -Jac(2,1)/detJ
       Jinv(2,2) =  Jac(1,1)/detJ
    else
       !WHL - do a proper abort here
       print*, 'stopping, det J = 0'
       print*, 'i, j, k, p:', i, j, k, p
       print*, 'Jacobian matrix:'
       print*, Jac(1,:)
       print*, Jac(2,:)
       stop
    endif

    if (verbose_Jac .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, ' '
       print*, 'Jacobian calc, p =', p
       print*, 'det J =', detJ
       print*, ' '
       print*, 'Jacobian matrix:'
       print*, Jac(1,:)
       print*, Jac(2,:)
       print*, ' '
       print*, 'Inverse matrix:'
       print*, Jinv(1,:)
       print*, Jinv(2,:)
       print*, ' '
       prod = matmul(Jac, Jinv)
       print*, 'Jac*Jinv:'
       print*, prod(1,:)
       print*, prod(2,:)
    endif

    ! bug check - Verify that J * Jinv = I

    prod = matmul(Jac,Jinv)
    do col = 1, 2
       do row = 1, 2
          if (abs(prod(row,col) - identity3(row,col)) > 1.d-12) then
             !TODO - do a proper abort here 
             print*, 'stopping, Jac * Jinv /= identity'
             print*, 'i, j, k, p:', i, j, k, p
             print*, 'Jac*Jinv:'
             print*, prod(1,:)
             print*, prod(2,:)
             stop
          endif
       enddo
    enddo

!TODO - Make this part optional?

    !------------------------------------------------------------------
    ! Compute the contribution of this quadrature point to dphi/dx and dphi/dy
    ! for each basis function.
    !
    !   | dphi_n/dx |          | dphi_n/dxr |
    !   |           | = Jinv * |            | 
    !   | dphi_n/dy |          | dphi_n/dyr |
    !
    !------------------------------------------------------------------

    dphi_dx(:) = 0.d0
    dphi_dy(:) = 0.d0

    do n = 1, nNodesPerElement_2d
       dphi_dx(n) = dphi_dx(n) + Jinv(1,1)*dphi_dxr(n)  &
                               + Jinv(1,2)*dphi_dyr(n)
       dphi_dy(n) = dphi_dy(n) + Jinv(2,1)*dphi_dxr(n)  &
                               + Jinv(2,2)*dphi_dyr(n)
    enddo

    ! debug - Check that the sum of dphi_dx, etc. is close to zero  
    !TODO - Think about why this should be the case

    if (abs(sum(dphi_dx)) > 1.d-12) then
        print*, 'stopping, sum over basis functions of dphi_dx > 0'
        print*, 'dphi_dx =', dphi_dx(:)
        print*, 'i, j, k, p =', i, j, k, p
        stop
    endif

    if (abs(sum(dphi_dy)) > 1.d-12) then
        print*, 'stopping, sum over basis functions of dphi_dy > 0'
        print*, 'dphi_dy =', dphi_dy(:)
        print*, 'i, j, k, p =', i, j, k, p
        stop
    endif

  end subroutine get_basis_function_derivatives_2d

!****************************************************************************

!WHL - debug - Pass in i, j, k, and p for now

  subroutine element_matrix_basal_sliding(nNodesPerElement_2d,   &
                                          wqp,         detJ,     &
                                          beta,        phi,      &
                                          Kuu,                   &
                                          ii, jj, k, p)

    !------------------------------------------------------------------
    ! Increment the matrix K corresponding to a 2D element face with the contribution
    ! from a particular quadrature point, given a linear sliding law.
    !
    ! For this sliding law, Kuu = Kvv, and Kvu = Kuv = 0.
    ! So it suffices here to compute Kuu.
    !
    ! TODO: Modify for more general sliding laws.
    !------------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

!WHL - debug - Pass in i, j, k, and p for now
    integer, intent(in) :: ii, jj, k, p

    integer, intent(in) :: nNodesPerElement_2d  ! number of nodes per 2D element face

    real(dp), intent(in) ::    &
             wqp,        &! weight for this quadrature point
             detJ,       &! determinant of Jacobian for the transformation
                          !  between the reference element and true element
             beta         ! basal traction factor beta at this quadrature point

    real(dp), dimension(nNodesPerElement_2d), intent(in) ::  &
             phi          ! basis functions phi evaluated at this quadrature point

    real(dp), dimension(nNodesPerElement_2d,nNodesPerElement_2d), intent(inout) :: &
             Kuu          ! components of stiffness matrix for this element face

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j

    if (verbose .and. this_rank==rtest .and. ii==itest .and. jj==jtest) then
       print*, ' '
       print*, 'Increment element matrix for basal BC, p =', p
    endif

    ! Increment the stiffness matrix for this quadrature point.
    
    do j = 1, nNodesPerElement_2d      ! columns of K
       do i = 1, nNodesPerElement_2d   ! rows of K

       if (verbose .and. this_rank==rtest .and. ii==itest .and. jj==jtest .and. p==ptest &
                   .and.  i==1     .and.  j==1) then
          print*, 'beta, wqp, detJ/vol0 =', beta, wqp, detJ/vol0
          print*, 'phi(1) =', phi(1)
          print*, 'Kuu phi increment(1,1) =', beta*wqp*detJ/vol0*phi(1)*phi(1)
       endif

!WHL - Note volume scaling such that detJ/vol0 is closer to unity
!TODO - Do we want vol0, or rather area0 = dx*dy?

          Kuu(i,j) = Kuu(i,j) + beta*wqp*detJ/vol0 * phi(i)*phi(j)

       enddo  ! i (rows)
    enddo     ! j (columns)

  end subroutine element_matrix_basal_sliding

!****************************************************************************

  subroutine element_to_global_matrix_2d(nx,           ny,          nz,          &
                                         nNodesPerElement_2d,                    &
                                         iElement,     jElement,    kLayer,      &
                                         Kmat,         Amat)

    ! Sum terms of element matrix K into dense assembled matrix A
    ! Here we assume that K is partitioned into Kuu, Kuv, Kvu, and Kvv,
    !  and similarly for A.
    ! For basal BC, the contribution to Kuv and Kvu is typically zero.
    ! So this subroutine is called twice per 2D element face.

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz                       ! number of vertical levels where velocity is computed

    integer, intent(in) ::   &
       nNodesPerElement_2d      ! number of nodes per 2D element face

    integer, intent(in) ::   &
       iElement, jElement,   &  ! i and j indices for this element
       kLayer                   ! k index for the layer corresponding to this 2D face
                                ! kLayer = nz for basal BC

    real(dp), dimension(nNodesPerElement_2d,nNodesPerElement_2d), intent(in) ::  &
       Kmat              ! element matrix

    real(dp), dimension(27,nz,nx-1,ny-1), intent(inout) ::    &
       Amat              ! assembled matrix

    integer :: i, j, k, m
    integer :: iA, jA, kA
    integer :: nr, nc

    if (verbose .and. this_rank==rtest .and. iElement==itest .and. jElement==jtest) then
       print*, 'Basal BC: First row of K:'
       write(6, '(8e12.4)') Kmat(1,:)
    endif

!TODO - Switch loops or switch order of indices in K(m,n)?
!       Inner index m is currently the outer loop

    k = kLayer   ! all nodes lie in one horizontal layer, typically k = nz for basal BC
    kA = 0       ! k index of A into which K(nr,n) is summed
                 ! = 0 since all nodes lie in the same layer

    do nr = 1, nNodesPerElement       ! rows of K

       ! Determine (i,j) for this node
       ! The reason for the '3' is that node 3, in the NE corner of the layer, has horizontal indices (i,j).
       ! Indices for other nodes are computed relative to this node.

       i = iElement + ishift(3,nr)
       j = jElement + jshift(3,nr)
      
       do nc = 1, nNodesPerElement    ! columns of K

          iA = ishift(nr,nc)          ! iA index of A into which K(nr,nc) is summed
          jA = jshift(nr,nc)          ! similarly for jA
                                      ! these indices can take values -1, 0 and 1

          m = indxA(iA,jA,kA)
          Amat(m,k,i,j) = Amat(m,k,i,j) + Kmat(nr,nc)

       enddo     ! n

    enddo        ! m

  end subroutine element_to_global_matrix_2d

!****************************************************************************

  subroutine dirichlet_boundary_conditions(nx,  ny,  nz,    nhalo,            &
                                           active_vertex,   umask_dirichlet,  &
                                           Auu,             Auv,              &
                                           Avu,             Avv,              &
                                           bu,              bv)

    !----------------------------------------------------------------
    ! Modify the global matrix and RHS for Dirichlet boundary conditions.
    ! This subroutine assumes that u = v = 0 at Dirichlet points.
    ! For each such point, we find the corresponding row and column of the global matrix.
    ! We zero out each row and column, except for setting the diagonal term to 1.
    ! We set the corresponding rhs to 0, so that the solution must be 0.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz,                   &  ! number of vertical levels where velocity is computed
       nhalo                    ! number of halo layers

    logical, dimension(nx-1,ny-1), intent(in) ::  &
       active_vertex       ! true for active vertices (vertices of active cells)

    logical, dimension(nz,nx-1,ny-1), intent(in) ::  &
       umask_dirichlet     ! Dirichlet mask for velocity (if true, u = 0)

    real(dp), dimension(27,nz,nx-1,ny-1), intent(inout) ::   &
       Auu, Auv,    &      ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv                                    

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::   &
       bu, bv              ! assembled load vector, divided into 2 parts

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------
    
    integer :: i, j, k     ! Cartesian indices of nodes
    integer :: iA, jA, kA  ! i, j, and k offsets of neighboring nodes 
    integer :: m

!PARALLEL - Make sure not to step out of bounds here.
!           This requires nhalo >= 2.

!WHL - debug
!       integer :: n, r
!       integer, dimension(nz,nx-1,ny-1) :: umask3  ! integer mask for umask_dirichlet
    
!          umask3(:,:,:) = 0
!          do j = nhalo, ny-nhalo+1
!             do i = nhalo, nx-nhalo+1
!                do k = 1, nz
!                   if (umask_dirichlet(k,i,j)) umask3(k,i,j) = 1
!                enddo
!             enddo
!          enddo

!          print*, ' '
!          print*, 'umask_dirichlet, k = 1:'
!          k = 1
!          do j = ny-nhalo, nhalo+1, -1
!             write(6,'(28i4)') umask3(k,3:30,j)
!          enddo

!          print*, ' '
!          print*, 'umask_dirichlet, k =', nz
!          k = nz
!          do j = ny-nhalo, nhalo+1, -1
!             write(6,'(28i4)') umask3(k,3:30,j)
!          enddo

! Loop over all vertices that border locally owned vertices.
! Locally owned vertices are (nhalo+1:ny-nhalo, nhalo+1:nx-nhalo)
!           
     do j = nhalo, ny-nhalo+1
        do i = nhalo, nx-nhalo+1
          if (active_vertex(i,j)) then
             do k = 1, nz
                if (umask_dirichlet(k,i,j)) then

                   ! loop through potential nonzero matrix values in the row associated with this vertex
                   ! Note: row = NodeID(k,i,j) if we were forming a matrix with one row per node

                   do kA = -1,1
                   do jA = -1,1
                   do iA = -1,1

                      if ( (k+kA >= 1 .and. k+kA <= nz)         &
                                      .and.                     &
                           (i+iA >= 1 .and. i+iA <= nx-1)       &
                                      .and.                     &
                           (j+jA >= 1 .and. j+jA <= ny-1) ) then

                         ! zero out A(row,col) and A(col,row)
                         ! Note: col = NodeID(k+kA, i+iA, j+jA) if we were forming a matrix with one column per node

                         ! Note: If we were setting (u,v) to something other than (0,0),
                         !  we would have to move terms to the RHS instead of just zeroing them.

!                            Auu( kA,  iA,  jA, row) = 0.d0
!                            Auu(-kA, -iA, -jA, col) = 0.d0
!                            Avv( kA,  iA,  jA, row) = 0.d0
!                            Avv(-kA, -iA, -jA, col) = 0.d0
!                            Auv( kA,  iA,  jA, row) = 0.d0
!                            Auv(-kA, -iA, -jA, col) = 0.d0
!                            Avu( kA,  iA,  jA, row) = 0.d0
!                            Avu(-kA, -iA, -jA, col) = 0.d0

                         m = indxA(iA,jA,kA)
                         Auu(m, k, i, j) = 0.d0
                         Auv(m, k, i, j) = 0.d0
                         Avu(m, k, i, j) = 0.d0
                         Avv(m, k, i, j) = 0.d0

                         m = indxA(-iA,-jA,-kA)
                         Auu(m, k+kA, i+iA, j+jA) = 0.d0
                         Auv(m, k+kA, i+iA, j+jA) = 0.d0
                         Avu(m, k+kA, i+iA, j+jA) = 0.d0
                         Avv(m, k+kA, i+iA, j+jA) = 0.d0

                      endif     ! i+iA, j+jA, and k+kA in bounds

                   enddo        ! kA
                   enddo        ! iA
                   enddo        ! jA

                   ! put back 1's on the main diagonal
                   m = indxA(0,0,0)
                   Auu(m,k,i,j) = 1.d0
                   Avv(m,k,i,j) = 1.d0

                   ! force u = v = 0 for this node by zeroing the rhs  
                   bu(k,i,j) = 0.d0            
                   bv(k,i,j) = 0.d0

                endif    ! umask_dirichlet
             enddo       ! k
          endif          ! active_vertex
       enddo             ! i
    enddo                ! j

  end subroutine dirichlet_boundary_conditions

!****************************************************************************

  subroutine compute_residual_vector(nx,          ny,            &
                                     nz,          nhalo,         &
                                     active_vertex,              &
                                     Auu,         Auv,           &
                                     Avu,         Avv,           &
                                     bu,          bv,            &
                                     uvel,        vvel,          &
                                     resid_vec_u, resid_vec_v,   &
                                     L2_norm)

    ! Compute the residual vector Ax - b and its L2 norm.
    ! This subroutine assumes that the matrix is stored in structured (x/y/z) format.

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions (for scalars)
       nz,                   &  ! number of vertical levels where velocity is computed
       nhalo                    ! number of halo layers

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F

    real(dp), dimension(27,nz,nx-1,ny-1), intent(in) ::   &
       Auu, Auv, Avu, Avv     ! four components of assembled matrix
                              ! 1st dimension = 3 (node and its nearest neighbors in x, y and z direction)
                              ! other dimensions = (z,x,y) indices
                              !
                              !    Auu  | Auv
                              !    _____|____
                              !    Avu  | Avv
                              !         |

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       bu, bv              ! assembled load (rhs) vector, divided into 2 parts

   real(dp), dimension(nz,nx-1,ny-1), intent(in) ::   &
       uvel, vvel          ! u and v components of velocity (m/s)

    real(dp), dimension(nz,nx-1,ny-1), intent(out) ::   &
       resid_vec_u,      & ! residual vector, divided into 2 parts
       resid_vec_v

    real(dp), intent(out) ::    &
       L2_norm             ! L2 norm of residual vector

    integer :: i, j, k, iA, jA, kA, m 

    if (verbose .and. this_rank==rtest) then
       print*, ' '
       i = itest
       j = jtest
       k = ktest
       print*, 'i, j, k, u answer, u rhs:', i, j, k, uvel(k,i,j), bu(k,i,j)
       print*, 'i, j, k, v answer, v rhs:', i, j, k, vvel(k,i,j), bv(k,i,j)
    endif

    ! Compute u and v components of A*x

    resid_vec_u(:,:,:) = 0.d0
    resid_vec_v(:,:,:) = 0.d0

!    do j = 1, ny-1
!    do i = 1, nx-1

    ! Loop over locally owned vertices

    do j = nhalo+1, ny-nhalo
    do i = nhalo+1, nx-nhalo
       !TODO - Use indirect addressing to avoid an 'if' here?
       if (active_vertex(i,j)) then

          do k = 1, nz

             do kA = -1,1
             do jA = -1,1
             do iA = -1,1

                if ( (k+kA >= 1 .and. k+kA <= nz)     &
                                .and.                     &
                     (i+iA >= 1 .and. i+iA <= nx-1)         &
                             .and.                     &
                     (j+jA >= 1 .and. j+jA <= ny-1) ) then

                   m = indxA(iA,jA,kA)

                   resid_vec_u(k,i,j) = resid_vec_u(k,i,j)                        & 
                                      + Auu(m,k,i,j)*uvel(k+kA,i+iA,j+jA)  &
                                      + Auv(m,k,i,j)*vvel(k+kA,i+iA,j+jA)

                   resid_vec_v(k,i,j) = resid_vec_v(k,i,j)                        &
                                      + Avu(m,k,i,j)*uvel(k+kA,i+iA,j+jA)  &
                                      + Avv(m,k,i,j)*vvel(k+kA,i+iA,j+jA)

                endif   ! in bounds

             enddo   ! kA
             enddo   ! iA
             enddo   ! jA

          enddo   ! k

       endif   ! active_vertex

    enddo   ! i
    enddo   ! j


    if (verbose .and. this_rank==rtest) then
!       print*, ' '
!       print*, 'Compute residual:'
!       print*, 'k, i, j, Axu, bu:'
    endif

    ! Subtract b to get A*x - b
    ! Sum up squared L2 norm as we go

    L2_norm = 0.d0

    ! Loop over locally owned vertices

    do j = nhalo+1, ny-nhalo
    do i = nhalo+1, nx-nhalo
       if (active_vertex(i,j)) then
          do k = 1, nz
             resid_vec_u(k,i,j) = resid_vec_u(k,i,j) - bu(k,i,j)
             resid_vec_v(k,i,j) = resid_vec_v(k,i,j) - bv(k,i,j)
             L2_norm = L2_norm + resid_vec_u(k,i,j)*resid_vec_u(k,i,j)  &
                               + resid_vec_v(k,i,j)*resid_vec_v(k,i,j)
          enddo  ! k
       endif     ! active vertex
    enddo        ! i
    enddo        ! j

    ! Take global sum, then take square root

    L2_norm = parallel_reduce_sum(L2_norm)
    L2_norm = sqrt(L2_norm)

  end subroutine compute_residual_vector

!****************************************************************************

  subroutine compute_residual_velocity(nhalo,  whichresid, &
                                       uvel,   vvel,        &
                                       usav,   vsav,        &
                                       resid_velo)

    integer, intent(in) ::   &
       nhalo,           & ! number of layers of halo cells
       whichresid         ! option for method to use when calculating residual

!TODO - Specify dimensions?
    real(dp), dimension(:,:,:), intent(in) ::  &
       uvel, vvel,      & ! current guess for velocity
       usav, vsav         ! previous guess for velocity

    real(dp), intent(out) ::    &
       resid_velo         ! quantity related to velocity convergence


    integer ::   &
       imaxdiff, jmaxdiff, kmaxdiff   ! location of maximum speed difference
                                      ! currently computed but not used
 
    integer :: i, j, k, count

    real(dp) ::   &
       speed,      &   ! current guess for ice speed
       oldspeed,   &   ! previous guess for ice speed
       diffspeed       ! abs(speed-oldspeed)


  ! Compute a residual quantity based on convergence of the velocity field.
  !TODO - Are all of these methods needed?

  ! options for residual calculation method, as specified in configuration file
  ! (see additional notes in "higher-order options" section of documentation)
  ! case(0): use max of abs( vel_old - vel ) / vel )
  ! case(1): use max of abs( vel_old - vel ) / vel ) but ignore basal vels
  ! case(2): use mean of abs( vel_old - vel ) / vel )
  ! case(3): use max of abs( vel_old - vel ) / vel ) (in addition to L2 norm)

    resid_velo = 0.d0
    imaxdiff = 0
    jmaxdiff = 0
    kmaxdiff = 0

    select case (whichresid)

    case(HO_RESID_MAXU_NO_UBAS)   ! max speed difference, excluding the bed

       ! Loop over locally owned vertices

       do j = 1+nhalo, size(uvel,3)-nhalo
          do i = 1+nhalo, size(uvel,2)-nhalo
             do k = 1, size(uvel,1) - 1         ! ignore bed velocity
                speed = sqrt(uvel(k,i,j)**2 + vvel(k,i,j)**2)
                if (speed /= 0.d0) then
                   oldspeed = sqrt(usav(k,i,j)**2 + vsav(k,i,j)**2)
                   diffspeed = abs((oldspeed - speed)/speed)
                   if (diffspeed > resid_velo) then
                      resid_velo = diffspeed
                      imaxdiff = i
                      jmaxdiff = j
                      kmaxdiff = k
                   endif
                endif
             enddo
          enddo
       enddo
       
       ! take global max
       resid_velo = parallel_reduce_max(resid_velo)

    case(HO_RESID_MEANU)   ! mean relative speed difference

       count = 0

       ! Loop over locally owned vertices

       do j = 1+nhalo, size(uvel,3)-nhalo
          do i = 1+nhalo, size(uvel,2)-nhalo
             do k = 1, size(uvel,1) - 1         ! ignore bed velocity
                speed = sqrt(uvel(k,i,j)**2 + vvel(k,i,j)**2)
                if (speed /= 0.d0) then
                   count = count+1
                   oldspeed = sqrt(usav(k,i,j)**2 + vsav(k,i,j)**2)
                   diffspeed = abs((oldspeed - speed)/speed)                
                   resid_velo = resid_velo + diffspeed
                endif
             enddo
          enddo
       enddo

       if (count > 0) resid_velo = resid_velo / count

!TODO - How to convert this to a global value?
!       Maybe remove this case
       call not_parallel(__FILE__, __LINE__)


   case default    ! max speed difference, including basal speeds
                   ! (case HO_RESID_MAXU or HO_RESID_L2NORM)

       ! Loop over locally owned vertices

       do j = 1+nhalo, size(uvel,3)-nhalo
          do i = 1+nhalo, size(uvel,2)-nhalo
             do k = 1, size(uvel,1)
                speed = sqrt(uvel(k,i,j)**2 + vvel(k,i,j)**2)
                if (speed /= 0.d0) then
                   oldspeed = sqrt(usav(k,i,j)**2 + vsav(k,i,j)**2)
                   diffspeed = abs((oldspeed - speed)/speed)
                   if (diffspeed > resid_velo) then
                      resid_velo = diffspeed
                      imaxdiff = i
                      jmaxdiff = j
                      kmaxdiff = k
                   endif
                endif
             enddo
          enddo
       enddo

       resid_velo = parallel_reduce_max(resid_velo)
       
  end select

  end subroutine compute_residual_velocity

!****************************************************************************
! The next three subroutines are used for the SLAP solver only.
!****************************************************************************

  subroutine slap_preprocess(nx,           ny,          &
                             nz,           nNodesSolve, &
                             NodeID,                    &
                             iNodeIndex,   jNodeIndex,  &
                             kNodeIndex,                &
                             Auu,          Auv,         &
                             Avu,          Avv,         &  
                             bu,           bv,          &
                             uvel,         vvel,        &
                             matrix_order, nNonzero,    &
                             matrix,       rhs,         &
                             answer)

    !----------------------------------------------------------------
    ! Using the intermediate matrices (Auu, Auv, Avu, Avv), load vectors (bu, bv),
    ! and velocity components (uvel, vvel), form the matrix and the rhs and answer
    ! vectors in the desired sparse matrix format.
    !
    ! The matrix is formed in ascending row order, so it can easily be transformed
    ! to compressed sparse row (CSR) format without further sorting.
    !
    ! Note: This works only for single-processor runs with the SLAP solver.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz,                   &  ! number of vertical levels at which velocity is computed
       nNodesSolve              ! number of nodes where we solve for velocity

    integer, dimension(nz,nx-1,ny-1), intent(in) ::  &
       NodeID             ! ID for each node

    integer, dimension(:), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of active nodes

    real(dp), dimension(27,nz,nx-1,ny-1), intent(in) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv           ! 1st dimension = node and its nearest neighbors in x, y and z direction 
                          ! other dimensions = (k,i,j) indices

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       bu, bv             ! assembled load (rhs) vector, divided into 2 parts
                          
    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::   &
       uvel, vvel         ! u and v components of velocity

    integer, intent(in) ::    &
       matrix_order,  &   ! order of matrix = number of rows
       nNonzero           ! upper bound for number of nonzero entries in sparse matrix

    type(sparse_matrix_type), intent(inout) ::  &    ! TODO: inout or out?
       matrix             ! sparse matrix, defined in glimmer_sparse_types
                          ! includes nonzeroes, order, col, row, val 

    real(dp), dimension(:), intent(out) ::   &
       rhs,             & ! right-hand-side (b) in Ax = b
       answer             ! answer (x) in Ax = b

    !---------------------------------------------------------
    ! Local variables
    !---------------------------------------------------------

    integer :: i, j, k, iA, jA, kA, m, mm, n, ct

    integer :: rowA, colA   ! row and column of A submatrices (order = nNodesSolve)
    integer :: row, col     ! row and column of sparse matrix (order = 2*nNodesSolve) 

    real(dp) :: val         ! value of matrix coefficient
    
    ! Set the nonzero coefficients of the sparse matrix 

    ct = 0

    do rowA = 1, nNodesSolve

       i = iNodeIndex(rowA)
       j = jNodeIndex(rowA)
       k = kNodeIndex(rowA)

       ! Load the nonzero values associated with Auu and Auv
       ! These are assigned a value of matrix%row = 2*rowA - 1

!TODO: Make sure that this procedure captures all the nonzero entries

       do kA = -1, 1
       do jA = -1, 1
       do iA = -1, 1

          if ( (k+kA >= 1 .and. k+kA <= nz)         &
                          .and.                     &
               (i+iA >= 1 .and. i+iA <= nx-1)       &
                          .and.                     &
               (j+jA >= 1 .and. j+jA <= ny-1) ) then

             colA = NodeID(k+kA, i+iA, j+jA)   ! ID for neighboring node
             m = indxA(iA,jA,kA)

             ! Auu
             val = Auu(m,k,i,j)
             if (val /= 0.d0) then
                ct = ct + 1
                matrix%row(ct) = 2*rowA - 1
                matrix%col(ct) = 2*colA - 1                    
                matrix%val(ct) = val
             endif

             ! Auv 
             val = Auv(m,k,i,j)
             if (val /= 0.d0) then
                ct = ct + 1
                matrix%row(ct) = 2*rowA - 1
                matrix%col(ct) = 2*colA
                matrix%val(ct) = val
             endif

             if (2*colA==rowtest .and. 2*rowA-1==coltest) then
               print*, ' '
               print*, 'rowA, colA, i, j, k =', rowA, colA, i, j, k
               print*, 'iA, jA, kA:', iA, jA, kA
               print*, 'Auv(iA,jA,kA,i,j,k) =', Auv(m,k,i,j)
               mm = indxA(-iA,-jA,-kA)
               print*, 'Avu(-iA,-jA,-kA,i+iA,j+jA,k+kA) =', Avu(mm, k+kA, i+iA, j+jA)
               print*, ' '
               print*, 'ct, row, col, val:', ct, matrix%row(ct), matrix%col(ct), matrix%val(ct)
             endif

          endif     ! i+iA, j+jA, and k+kA in bounds

       enddo        ! kA
       enddo        ! iA
       enddo        ! jA

       ! Load the nonzero values associated with Avu and Avv
       ! These are assigned a value of matrix%row = 2*rowA

       do kA = -1, 1
       do jA = -1, 1
       do iA = -1, 1

          if ( (k+kA >= 1 .and. k+kA <= nz)         &
                          .and.                     &
               (i+iA >= 1 .and. i+iA <= nx-1)       &
                          .and.                     &
               (j+jA >= 1 .and. j+jA <= ny-1) ) then

             colA = NodeID(k+kA, i+iA, j+jA)   ! ID for neighboring node
             m  = indxA( iA, jA, kA)
             mm = indxA(-iA,-jA,-kA)

             ! Avu 
             val = Avu(m,k,i,j)
             if (val /= 0.d0) then
                ct = ct + 1
                matrix%row(ct) = 2*rowA
                matrix%col(ct) = 2*colA - 1
                matrix%val(ct) = val
             endif

!WHL - debug
             if (2*rowA==rowtest .and. 2*colA-1==coltest) then
               print*, ' '
               print*, 'rowA, colA, i, j, k =', rowA, colA, i, j, k
               print*, 'iA, jA, kA:', iA, jA, kA
               print*, 'Avu(iA,jA,kA,i,j,k) =', Avu(m,k,i,j)
               print*, 'Auv(-iA,-jA,-kA,i+iA,j+jA,k+kA) =', Auv(mm, k+kA, i+iA, j+jA)
               print*, ' '
               print*, 'ct, row, col, val:', ct, matrix%row(ct), matrix%col(ct), matrix%val(ct)
             endif

             ! Avv 
             val = Avv(m,k,i,j)
             if (val /= 0.d0) then
                ct = ct + 1
                matrix%row(ct) = 2*rowA
                matrix%col(ct) = 2*colA
                matrix%val(ct) = val
             endif

          endif     ! i+iA, j+jA, and k+kA in bounds

       enddo        ! kA
       enddo        ! iA
       enddo        ! jA

    enddo           ! rowA 

    ! Set basic matrix parameters.

    matrix%order = matrix_order
    matrix%nonzeros = ct
    matrix%symmetric = .false.

    if (verbose) then
       print*, ' '
       print*, 'solver preprocess'
       print*, 'order, nonzeros =', matrix%order, matrix%nonzeros
    endif

    ! Initialize the answer vector
    ! For efficiency, put the u and v terms for a given node adjacent in storage.

    do n = 1, nNodesSolve
       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)

       answer(2*n-1) = uvel(k,i,j)
       answer(2*n)   = vvel(k,i,j)

       if (verbose .and. n==ntest) then
          print*, ' '
          print*, 'n, initial uvel =', n, answer(2*n-1)
          print*, 'n, initial vvel =', n, answer(2*n)
       endif

    enddo

    ! Set the rhs vector
    ! For efficiency, put the u and v terms for a given node adjacent in storage.

    do n = 1, nNodesSolve
       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)

       rhs(2*n-1) = bu(k,i,j)
       rhs(2*n)   = bv(k,i,j)

       if (verbose .and. n==ntest) then
          print*, ' '
          print*, 'n, initial bu =', n, rhs(2*n-1)
          print*, 'n, initial bv =', n, rhs(2*n)
       endif

    enddo

  end subroutine slap_preprocess

!****************************************************************************

  subroutine slap_postprocess(nNodesSolve,  answer,                   &
                              iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                              uvel,         vvel)

  ! Extract the velocities from the solution vector.
  ! Note: This works only for single-processor runs with the SLAP solver.
                                            
    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) :: nNodesSolve     ! number of nodes where we solve for velocity

    real(dp), dimension(:), intent(in) ::  &
       answer             ! velocity solution vector

    integer, dimension(:), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex  ! i, j and k indices of active nodes

    real(dp), dimension(:,:,:), intent(inout) ::   &
       uvel, vvel         ! u and v components of velocity

    integer :: i, j, k, n

       if (verbose) then
          print*, ' '
          print*, 'solver postprocess: n, i, j, k, u, v'
       endif

    do n = 1, nNodesSolve

       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)

       uvel(k,i,j) = answer(2*n-1)
       vvel(k,i,j) = answer(2*n)

    enddo

  end subroutine slap_postprocess

!****************************************************************************

  subroutine slap_compute_residual_vector(matrix,    answer,   rhs, &
                                          resid_vec, L2_norm)

    ! Compute the residual vector Ax - b and its L2 norm.
    ! This subroutine assumes that the matrix is stored in triad (row/col/val) format.

    type(sparse_matrix_type), intent(in) ::  &
       matrix             ! sparse matrix, defined in glimmer_sparse_types
                          ! includes nonzeroes, order, col, row, val 

    real(dp), dimension(:), intent(in) ::   &
       rhs,             & ! right-hand-side (b) in Ax = b
       answer             ! answer (x) in Ax = b

    real(dp), dimension(:), intent(out) ::   &
       resid_vec          ! residual vector

    real(dp), intent(out) ::    &
       L2_norm             ! L2 norm of residual vector

    integer :: i, j, n

    if (verbose) then
       print*, ' '
       print*, 'Residual vector: order, nonzeros =', matrix%order, matrix%nonzeros 
       print*, 'ntest, u answer, u rhs:', ntest, answer(2*ntest-1), rhs(2*ntest-1)
       print*, 'ntest, v answer, v rhs:', ntest, answer(2*ntest),   rhs(2*ntest)
    endif

    resid_vec(:) = 0.d0

    do n = 1, matrix%nonzeros
       i = matrix%row(n)
       j = matrix%col(n)
       resid_vec(i) = resid_vec(i) + matrix%val(n)*answer(j)
    enddo

    L2_norm = 0.d0
    do i = 1, matrix%order
       resid_vec(i) = resid_vec(i) - rhs(i)
       L2_norm = L2_norm + resid_vec(i)*resid_vec(i)
    enddo

!PARALLEL - Parallel sum needed here?  Or will we never call this subroutine in parallel?
!    L2_norm = parallel_reduce_sum(L2_norm)

    L2_norm = sqrt(L2_norm)

  end subroutine slap_compute_residual_vector

!****************************************************************************
! The remaining subroutines are used for testing and bug-checking but may 
! not need to be called for production runs.
!****************************************************************************

  subroutine check_symmetry_element_matrix(nNodesPerElement,  &
                                           Kuu, Kuv, Kvu, Kvv)

    !------------------------------------------------------------------
    ! Check that the element stiffness matrix is symmetric.
    ! This is true provided that (1) Kuu = (Kuu)^T
    !                            (2) Kvv = (Kvv)^T
    !                            (3) Kuv = (Kvu)^T
    ! This check should not be needed for production runs with a well-tested code,
    !  but is included for now to help with debugging.
    !------------------------------------------------------------------

    integer, intent(in) :: nNodesPerElement  ! number of nodes per element

    real(dp), dimension(nNodesPerElement, nNodesPerElement), intent(in) ::   &
             Kuu, Kuv, Kvu, Kvv     ! component of element stiffness matrix
                                    !
                                    !    Kuu  | Kuv
                                    !    _____|____          
                                    !    Kvu  | Kvv
                                    !         |

    integer :: i, j

    ! check that Kuu = (Kuu)^T

    do j = 1, nNodesPerElement
       do i = j, nNodesPerElement
          if (abs(Kuu(i,j) - Kuu(j,i)) > eps12) then
             print*, 'Kuu is not symmetric'
             print*, 'i, j, Kuu(i,j), Kuu(j,i):', i, j, Kuu(i,j), Kuu(j,i)
             stop
          endif    
       enddo
    enddo

    ! check that Kvv = (Kvv)^T

    do j = 1, nNodesPerElement
       do i = j, nNodesPerElement
          if (abs(Kvv(i,j) - Kvv(j,i)) > eps12) then
             print*, 'Kvv is not symmetric'
             print*, 'i, j, Kvv(i,j), Kvv(j,i):', i, j, Kvv(i,j), Kvv(j,i)
             stop
          endif    
       enddo
    enddo

    ! Check that Kuv = (Kvu)^T

    do j = 1, nNodesPerElement
       do i = 1, nNodesPerElement
          if (abs(Kuv(i,j) - Kvu(j,i)) > eps12) then
             print*, 'Kuv .ne. (Kvu)^T'
             print*, 'i, j, Kuv(i,j), Kvu(j,i):', i, j, Kuv(i,j), Kvu(j,i)
             stop
          endif    
       enddo
    enddo

  end subroutine check_symmetry_element_matrix

!****************************************************************************

  subroutine check_symmetry_assembled_matrix(nx,  ny,  nz, nhalo,   &
                                             active_vertex,         &
                                             Auu, Auv, Avu, Avv)

    !------------------------------------------------------------------
    ! Check that the assembled stiffness matrix is symmetric.
    ! This is true provided that (1) Auu = (Auu)^T
    !                            (2) Avv = (Avv)^T
    !                            (3) Auv = (Avu)^T
    ! The A matrices are assembled in a dense fashion to save storage
    !  and preserve the i/j/k structure of the grid..
    !
    ! There can be small differences from perfect symmetry due to roundoff error.
    ! These differences are fixed provided they are small enough.
    !------------------------------------------------------------------    

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz,                   &  ! number of vertical levels where velocity is computed
       nhalo                    ! number of halo layers

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex            ! T for columns (i,j) where velocity is computed, else F

    real(dp), dimension(27,nz,nx-1,ny-1), intent(inout) ::   &
       Auu, Auv, Avu, Avv       ! components of assembled stiffness matrix
                                !
                                !    Auu  | Auv
                                !    _____|____          
                                !         |
                                !    Avu  | Avv                                    

    integer :: i, j, k, iA, jA, kA, m, mm

    real(dp) :: val1, val2          ! values of matrix coefficients

    real(dp) :: maxdiff, diag_entry, avg_val

    ! Check matrix for symmetry

    ! Here we correct for small differences from symmetry due to roundoff error.
    ! The maximum departure from symmetry is set to be a small fraction (1.e-12) 
    !  of the diagonal entry for the row.
    ! If the departure from symmetry is larger than this, then the model prints a warning 
    !  and/or aborts.

    maxdiff = 0.d0

! Loop over locally owned vertices.
! Each active vertex is associate with 2*nz matrix rows belonging to this processor.
! Locally owned vertices are (nhalo+1:ny-nhalo, nhalo+1:nx-nhalo)
!           
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (active_vertex(i,j)) then
             do k = 1, nz

                ! Check Auu and Auv for symmetry

                m = indxA(0,0,0)
                diag_entry = Auu(m,k,i,j)

                do kA = -1, 1
                do jA = -1, 1
                do iA = -1, 1

                   m =  indxA( iA, jA, kA)
                   mm = indxA(-iA,-jA,-kA)

                   ! Check that Auu = Auu^T

                   if (k+kA >= 1 .and. k+kA <=nz) then  ! to keep k index in bounds

                      val1 = Auu( m, k,    i,    j   )   ! value of Auu(row,col)
                      val2 = Auu(mm, k+kA, i+iA, j+jA)   ! value of Auu(col,row)

                      if (val2 /= val1) then
                          
                         if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)

                         ! if difference is small, then fix the asymmetry by averaging values
                         ! else print a warning and abort

                         if ( abs(val2-val1) < eps12*diag_entry ) then
                            avg_val = 0.5d0 * (val1 + val2)
                            Auu( m, k,   i,   j   ) = avg_val
                            Auu(mm, k+kA,i+iA,j+jA) = avg_val
!!                            print*, 'Auu, i, j, k, iA, jA, kA, val1, val2:', i, j, k, iA, jA, kA, val1, val2
                         else
                            print*, ' '
                            print*, 'Auu is not symmetric: i, j, k, iA, jA, kA =', i, j, k, iA, jA, kA
                            print*, 'Auu(row,col), Auu(col,row), diff:', val1, val2, val2 - val1
                            stop  !TODO - Put in a proper abort
                         endif

                      endif   ! val2 /= val1
                
                      ! Check that Auv = (Avu)^T

                      val1 = Auv( m, k,    i,    j)      ! value of Auv(row,col)
                      val2 = Avu(mm, k+kA, i+iA, j+jA)   ! value of Avu(col,row)

                      if (val2 /= val1) then

                         if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)

                         ! if difference is small, then fix the asymmetry by averaging values
                         ! else print a warning and abort

                         if ( abs(val2-val1) < eps12*diag_entry ) then
                            avg_val = 0.5d0 * (val1 + val2)
                            Auv( m, k,   i,   j   ) = avg_val
                            Avu(mm, k+kA,i+iA,j+jA) = avg_val
!!                            print*, 'Auv, i, j, k, iA, jA, kA, val1, val2:', i, j, k, iA, jA, kA, val1, val2
                         else
                            print*, ' '
                            print*, 'Auv is not equal to (Avu)^T, i, j, k, iA, jA, kA =', i, j, k, iA, jA, kA
                            print*, 'Auv(row,col), Avu(col,row), diff:', val1, val2, val2 - val1
                            stop  !TODO - Put in a proper abort
                         endif

                      endif  ! val2 /= val1

                   endif     ! k+kA in bounds
            
                enddo        ! kA
                enddo        ! iA
                enddo        ! jA

                ! Now check Avu and Avv

                m = indxA(0,0,0)
                diag_entry = Avv(m,k,i,j)

                ! check that Avv = (Avv)^T

                do kA = -1, 1
                do jA = -1, 1
                do iA = -1, 1

                   if (k+kA >= 1 .and. k+kA <=nz) then  ! to keep k index in bounds

                      m  = indxA( iA, jA, kA)
                      mm = indxA(-iA,-jA,-kA)

                      val1 = Avv( m, k,    i,    j)      ! value of Avv(row,col)
                      val2 = Avv(mm, k+kA, i+iA, j+jA)   ! value of Avv(col,row)

                      if (val2 /= val1) then
                          
                         if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)

                         ! if difference is small, then fix the asymmetry by averaging values
                         ! else print a warning and abort

                         if ( abs(val2-val1) < eps12*diag_entry ) then
                            avg_val = 0.5d0 * (val1 + val2)
                            Avv( m, k,   i,   j   ) = avg_val
                            Avv(mm, k+kA,i+iA,j+jA) = avg_val
!!                            print*, 'Avv, i, j, k, iA, jA, kA, val1, val2:', i, j, k, iA, jA, kA, val1, val2
                         else
                            print*, ' '
                            print*, 'Avv is not symmetric: i, j, k, iA, jA, kA =', i, j, k, iA, jA, kA
                            print*, 'Avv(row,col), Avv(col,row), diff:', val1, val2, val2 - val1
                            stop  !TODO - Put in a proper abort
                         endif

                      endif   ! val2 /= val1

                      ! Check that Avu = (Auv)^T

                      val1 = Avu( m, k,    i,    j)      ! value of Avu(row,col)
                      val2 = Auv(mm, k+kA, i+iA, j+jA)   ! value of Auv(col,row)

                      if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)

                      if (val2 /= val1) then

                         ! if difference is small, then fix the asymmetry by averaging values
                         ! else print a warning and abort

                         if ( abs(val2-val1) < eps12*diag_entry ) then
                            avg_val = 0.5d0 * (val1 + val2)
                            Avu( m, k,   i,   j   ) = avg_val
                            Auv(mm, k+kA,i+iA,j+jA) = avg_val
!!                            print*, 'Avu, i, j, k, iA, jA, kA, val1, val2:', i, j, k, iA, jA, kA, val1, val2
                         else
                            print*, ' '
                            print*, 'Avu is not equal to (Auv)^T, i, j, k, iA, jA, kA =', i, j, k, iA, jA, kA
                            print*, 'Avu(row,col), Auv(col,row), diff:', val1, val2, val2 - val1
                            stop  !TODO - Put in a proper abort
                         endif

                      endif  ! val2 /= val1

                   endif     ! k+kA in bounds

                enddo        ! kA
                enddo        ! iA
                enddo        ! jA

             enddo     ! k
          endif        ! active_vertex
       enddo           ! i
    enddo              ! j

    maxdiff = parallel_reduce_max(maxdiff)
    if (verbose .and. main_task) then
       print*, ' '
       print*, 'Max difference from symmetry =', maxdiff
    endif

  end subroutine check_symmetry_assembled_matrix

!****************************************************************************

  subroutine solve_test_matrix (matrix_order, whichsparse)

    ! solve a small test matrix

    integer, intent(in) :: &
       matrix_order,       & ! matrix order
       whichsparse           ! solution method (0=BiCG, 1=GMRES, 2=PCG_DIAG, 3=PCG_INCH, 5 = STANDALONE_PCG)

    logical :: verbose_test = .true.

    type(sparse_matrix_type) ::  &
       matrix             ! sparse matrix, defined in glimmer_sparse_types

    real(dp), dimension(:), allocatable ::   &
       rhs,             & ! right-hand-side (b) in Ax = b
       answer             ! answer (x) in Ax = b

    real(dp), dimension(:,:), allocatable :: Atest

    real(dp) :: err

    integer :: niters, nNonzero_max

    integer :: i, j, n

    print*, 'Solving test matrix, order =', matrix_order

    nNonzero_max = matrix_order*matrix_order    ! not sure how big this must be

    allocate(Atest(matrix_order,matrix_order))
    Atest(:,:) = 0.d0

    allocate(matrix%row(nNonzero_max), matrix%col(nNonzero_max), matrix%val(nNonzero_max))
    allocate(rhs(matrix_order), answer(matrix_order))

    rhs(:) = 0.d0
    answer(:) = 0.d0
    matrix%row(:) = 0
    matrix%col(:) = 0
    matrix%val(:) = 0.d0

    matrix%order = matrix_order
    matrix%symmetric = .false.

    if (matrix%order == 2) then    ! symmetric 2x2
       Atest(1,1:2) = (/3.d0, 2.d0 /)
       Atest(2,1:2) = (/2.d0, 6.d0 /)
       rhs(1:2) = (/2.d0, -8.d0 /)   ! answer = (2 -2) 

    elseif (matrix%order == 3) then

       ! symmetric
       Atest(1,1:3) = (/ 7.d0, -2.d0,  0.d0 /)
       Atest(2,1:3) = (/-2.d0,  6.d0, -2.d0 /)
       Atest(3,1:3) = (/ 0.d0, -2.d0,  5.d0 /)
       rhs(1:3)   =   (/ 3.d0,  8.d0,  1.d0 /)   ! answer = (1 2 1)

        ! non-symmetric
!       Atest(1,1:3) = (/3.d0,   1.d0,  1.d0 /)
!       Atest(2,1:3) = (/2.d0,   2.d0,  5.d0 /)
!       Atest(3,1:3) = (/1.d0,  -3.d0, -4.d0 /)
!       rhs(1:3)   =  (/ 6.d0,  11.d0, -9.d0 /)   ! answer = (1 2 1)

    else if (matrix%order == 4) then

        ! symmetric

       Atest(1,1:4) = (/ 2.d0, -1.d0,  0.d0,  0.d0 /)
       Atest(2,1:4) = (/-1.d0,  2.d0, -1.d0,  0.d0 /)
       Atest(3,1:4) = (/ 0.d0, -1.d0,  2.d0, -1.d0 /)
       Atest(4,1:4) = (/ 0.d0,  0.d0, -1.d0,  2.d0 /)
       rhs(1:4)    = (/  0.d0,  1.d0, -1.d0,  4.d0 /)   ! answer = (1 2 2 3)

        ! non-symmetric
!       Atest(1,1:4) = (/3.d0,  0.d0,  2.d0, -1.d0 /)
!       Atest(2,1:4) = (/1.d0,  2.d0,  0.d0,  2.d0 /)
!       Atest(3,1:4) = (/4.d0,  0.d0,  6.d0, -3.d0 /)
!       Atest(4,1:4) = (/5.d0,  0.d0,  2.d0,  0.d0 /)
!       rhs(1:4)    = (/ 6.d0,  7.d0, 13.d0,  9.d0 /)   ! answer = (1 2 2 1)

    elseif (matrix%order > 4) then
  
        Atest(:,:) = 0.d0
        do n = 1, matrix%order 
           Atest(n,n) = 2.d0
           if (n > 1) Atest(n,n-1) = -1.d0
           if (n < matrix%order) Atest(n,n+1) = -1.d0
        enddo

        rhs(1) = 1.d0
        rhs(matrix%order) = 1.d0
        rhs(2:matrix%order-1) = 0.d0              ! answer = (1 1 1 ... 1 1 1)
        
    endif

    if (verbose_test) then
       print*, ' '
       print*, 'Atest =', Atest
       print*, 'rhs =', rhs
    endif

    ! Put in SLAP triad format (column ascending order)

    n = 0
    do j = 1, matrix%order
       do i = 1, matrix%order
          if (Atest(i,j) /= 0.d0) then 
             n = n + 1
             matrix%row(n) = i
             matrix%col(n) = j
             matrix%val(n) = Atest(i,j)
          endif
       enddo
    enddo

    ! Set number of nonzero values
    matrix%nonzeros = n

    if (verbose_test) then
       print*, ' '
       print*, 'row,       col,       val:'
       do n = 1, matrix%nonzeros
          print*, matrix%row(n), matrix%col(n), matrix%val(n)
       enddo
       print*, 'Call sparse_easy_solve, whichsparse =', whichsparse
    endif

    ! Solve the linear matrix problem


    call sparse_easy_solve(matrix, rhs,    answer,  &
                           err,    niters, whichsparse)

    if (verbose_test) then
       print*, ' '
       print*, 'answer =', answer
       print*, 'err =', err       
       print*, 'niters =', niters
    endif

    stop

  end subroutine solve_test_matrix

!****************************************************************************

!WHL - This subroutine is currently not called.
!      It might be useful if the matrix contains many entries that are
!       close to but not quite equal to zero (e.g., due to roundoff error).
!      It may need to be modified to preserve symmetry (so we don't zero out 
!       an entry without zeroing out its symmetric partner).

  subroutine remove_small_matrix_terms(nx, ny, nz,   &
                                       Auu, Auv,     &
                                       Avu, Avv) 


    !------------------------------------------------------------------------
    ! Get rid of matrix entries that are very small compared to diagonal entries.
    !   
    ! This handles a potential numerical issue: Some entries should be zero because of 
    !  cancellation of terms, but are not exactly zero because of rounding errors.
    !------------------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz                            ! number of vertical layers at which velocity is computed

    real(dp), dimension(27,nz,nx-1,ny-1), intent(inout) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv                                    

    !---------------------------------------------------------
    ! Local variables
    !---------------------------------------------------------

    integer :: i, j, k, m, iA, jA, kA

    real(dp) ::       &
       diag_entry,   &! mean size of diagonal entries
       min_entry      ! minimum size of entries kept in matrix

    do j = 1, ny-1
    do i = 1, nx-1
    do k = 1, nz

       ! Compute threshold value of matrix entries for Auu/Auv row

       m = indxA(0,0,0)
       diag_entry = Auu(m,k,i,j)
       min_entry = eps12 * diag_entry

       ! Remove very small values

       do kA = -1, 1
       do jA = -1, 1
       do iA = -1, 1

          m = indxA(iA,jA,kA)
          if (abs(Auu(m,k,i,j)) < min_entry) then
!             if (abs(Auu(m,k,i,j)) > 0.d0) then
!                print*, 'Dropping Auu term:, i, j, k, iA, jA, kA, val:', i, j, k, iA, jA, kA, Auu(m,k,i,j)
!             endif
             Auu(m,k,i,j) = 0.d0
          endif

          if (abs(Auv(m,k,i,j)) < min_entry) then
!             if (abs(Auv(m,k,i,j)) > 0.d0) then
!                print*, 'Dropping Auv term:, i, j, k, iA, jA, kA, val:', i, j, k, iA, jA, kA, Auv(m,k,i,j)
!             endif
             Auv(m,k,i,j) = 0.d0
          endif

       enddo
       enddo
       enddo

       ! Compute threshold value of matrix entries for Avu/Avv row

       m = indxA(0,0,0)
       diag_entry = Avv(m,k,i,j)
       min_entry = eps12 * diag_entry

       ! Remove very small values

       do kA = -1, 1
       do jA = -1, 1
       do iA = -1, 1

          m = indxA(iA,jA,kA)

          if (abs(Avu(m,k,i,j)) < min_entry) then
!             if (abs(Avu(m,k,i,j)) > 0.d0) then
!                print*, 'Dropping Avu term:, i, j, k, iA, jA, kA, val:', i, j, k, iA, jA, kA, Avu(m,k,i,j)
!             endif
             Avu(m,k,i,j) = 0.d0
          endif

          if (abs(Avv(m,k,i,j)) < min_entry) then
!             if (abs(Avv(m,k,i,j)) > 0.d0) then
!                print*, 'Dropping Avv term:, i, j, k, iA, jA, kA, val:', i, j, k, iA, jA, kA, Avv(m,k,i,j)
!             endif
             Avv(m,k,i,j) = 0.d0
          endif

       enddo
       enddo
       enddo

  enddo
  enddo
  enddo
     
  end subroutine remove_small_matrix_terms

!****************************************************************************

  end module glissade_velo_higher

!****************************************************************************
