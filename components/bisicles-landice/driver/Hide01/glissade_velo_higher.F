!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_velo_higher.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This module contains routines for computing the ice velocity using a 
! variational finite-element approach.  It solves the higher-order Blatter-Pattyn
! approximation for Stokes flow, as well as several simpler approximations
! (L1L2, shallow-shelf approximation, and shallow-ice approximation).
!
! See these papers for details:
!
! J.K. Dukowicz, S.F. Price and W.H. Lipscomb, 2010: Consistent
!    approximations and boundary conditions for ice-sheet dynamics
!    using a principle of least action.  J. Glaciology, 56 (197),
!    480-495.
!
! F. Pattyn, 2003: A new three-dimensional higher-order thermomechanical 
!    ice sheet model: Basic sensitivity, ice stream development, and
!    ice flow across subglacial lakes.  J. Geophys. Res., 108 (B8),
!    2382, doi:10.1029/2002JB002329.
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
    use glimmer_physcon, only: gn, rhoi, rhoo, grav, scyr, pi
    use glimmer_paramets, only: thk0, len0, tim0, tau0, vel0, vis0, evs0
    use glimmer_paramets, only: vel_scale, len_scale   ! used for whichefvs = HO_EFVS_FLOWFACT
    use glimmer_log
    use glimmer_sparse_type
    use glimmer_sparse
    use glissade_grid_operators     
    use glissade_masks, only: glissade_get_masks, glissade_grounded_fraction

    use glide_types

    use glissade_velo_higher_slap, only:   &
         slap_preprocess_3d,   slap_preprocess_2d,   &
         slap_postprocess_3d,  slap_postprocess_2d,  &
         slap_compute_residual_vector, slap_solve_test_matrix

    use glissade_velo_higher_pcg, only:   &
         pcg_solver_standard_3d,   pcg_solver_standard_2d,  &
         pcg_solver_chrongear_3d,  pcg_solver_chrongear_2d, &
         matvec_multiply_structured_3d

#ifdef TRILINOS
    use glissade_velo_higher_trilinos, only: &
         trilinos_fill_pattern_3d,     trilinos_fill_pattern_2d,     &
         trilinos_global_id_3d,        trilinos_global_id_2d,        &
         trilinos_assemble_3d,         trilinos_assemble_2d,         &
         trilinos_init_velocity_3d,    trilinos_init_velocity_2d,    &
         trilinos_extract_velocity_3d, trilinos_extract_velocity_2d, &
         trilinos_test
#endif

    use parallel

    implicit none

    private
    public :: glissade_velo_higher_init, glissade_velo_higher_solve

    !----------------------------------------------------------------
    ! Here are some definitions:
    !
    ! The horizontal mesh is composed of cells and vertices.
    ! The cells are rectangular with uniform dimensions dx and dy.
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
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Finite element properties
    ! Assume 3D hexahedral elements.
    !----------------------------------------------------------------

    integer, parameter ::        &
       nNodesPerElement_3d = 8,  & ! 8 nodes for hexahedral elements
       nQuadPoints_3d = 8,       & ! number of quadrature points per hexahedral element
                                   ! These live at +- 1/sqrt(3) for reference hexahedron
       nNodeNeighbors_3d = 27      ! number of nearest node neighbors in 3D (including the node itself)

    integer, parameter ::        &
       nNodesPerElement_2d = 4,  & ! 4 nodes for faces of hexahedral elements
       nQuadPoints_2d = 4,       & ! number of quadrature points per element face
                                   ! These live at +- 1/sqrt(3) for reference square
       nNodeNeighbors_2d = 9       ! number of nearest node neighbors in 2D (including the node itself)

    real(dp), parameter ::     &
       rsqrt3 = 1.d0/sqrt(3.d0)    ! for quadrature points
         
    !----------------------------------------------------------------
    ! Arrays used for finite-element calculations
    !
    ! Most integals are done over 3D hexahedral elements.
    ! Surface integrals are done over 2D faces of these elements. 
    !----------------------------------------------------------------

    real(dp), dimension(nNodesPerElement_3d, nQuadPoints_3d) ::   & 
       phi_3d,         &    ! trilinear basis function, evaluated at quad pts
       dphi_dxr_3d,    &    ! dphi/dx for reference hexehedral element, evaluated at quad pts
       dphi_dyr_3d,    &    ! dphi/dy for reference hexahedral element, evaluated at quad pts
       dphi_dzr_3d          ! dphi/dy for reference hexahedral element, evaluated at quad pts

    real(dp), dimension(nNodesPerElement_3d) ::   & 
       phi_3d_ctr,         &! trilinear basis function, evaluated at cell ctr
       dphi_dxr_3d_ctr,    &! dphi/dx for reference hexahedral element, evaluated at cell ctr
       dphi_dyr_3d_ctr,    &! dphi/dy for reference hexahedral element, evaluated at cell ctr
       dphi_dzr_3d_ctr      ! dphi/dz for reference hexahedral element, evaluated at cell ctr

    real(dp), dimension(nQuadPoints_3d) ::  &
       xqp_3d, yqp_3d, zqp_3d,  &! quad pt coordinates in reference element
       wqp_3d                    ! quad pt weights

    real(dp), dimension(nNodesPerElement_2d, nQuadPoints_2d) ::   & 
       phi_2d,         &    ! bilinear basis function, evaluated at quad pts
       dphi_dxr_2d,    &    ! dphi/dx for reference rectangular element, evaluated at quad pts
       dphi_dyr_2d          ! dphi/dy for reference rectangular element, evaluated at quad pts

    real(dp), dimension(nNodesPerElement_2d) ::   & 
       phi_2d_ctr,         &! bilinear basis function, evaluated at cell ctr
       dphi_dxr_2d_ctr,    &! dphi/dx for reference rectangular element, evaluated at cell ctr
       dphi_dyr_2d_ctr      ! dphi/dy for reference rectangular element, evaluated at cell ctr

    real(dp), dimension(nQuadPoints_2d) ::  &
       xqp_2d, yqp_2d, &    ! quad pt coordinates in reference square
       wqp_2d               ! quad pt weights

    integer, dimension(nNodesPerElement_3d, nNodesPerElement_3d) ::  &
       ishift, jshift, kshift   ! matrices describing relative indices of nodes in an element

    integer, dimension(-1:1,-1:1,-1:1) :: &
       indxA_3d              ! maps relative (x,y,z) coordinates to an index between 1 and 27
                             ! index order is (i,j,k)

    integer, dimension(-1:1,-1:1) :: &
       indxA_2d              ! maps relative (x,y) coordinates to an index between 1 and 9
                             ! index order is (i,j)

    real(dp), dimension(3,3) ::  &
       identity3             ! 3 x 3 identity matrix

    real(dp), parameter ::   &
       eps08 = 1.d-08,      &! small number
       eps10 = 1.d-10        ! smaller number

    real(dp) :: vol0    ! volume scale (m^3), used to scale 3D matrix values

    logical, parameter ::  &
       check_symmetry = .true.   ! if true, then check symmetry of assembled matrix

    ! various options for turning diagnostic prints on and off
    logical :: verbose = .false.
!    logical :: verbose = .true.  
    logical :: verbose_init = .false.   
!    logical :: verbose_init = .true.   
    logical :: verbose_Jac = .false.
!    logical :: verbose_Jac = .true.
    logical :: verbose_residual = .false.
!    logical :: verbose_residual = .true.
    logical :: verbose_state = .false.
!    logical :: verbose_state = .true.
    logical :: verbose_velo = .false.
!    logical :: verbose_velo = .true.
    logical :: verbose_id = .false.
!    logical :: verbose_id = .true.
    logical :: verbose_load = .false.
!    logical :: verbose_load = .true.
    logical :: verbose_shelf = .false.
!    logical :: verbose_shelf = .true.
    logical :: verbose_matrix = .false.
!    logical :: verbose_matrix = .true.
    logical :: verbose_basal = .false.
!    logical :: verbose_basal = .true.
    logical :: verbose_bfric = .false.
!    logical :: verbose_bfric = .true.
    logical :: verbose_trilinos = .false.
!    logical :: verbose_trilinos = .true.
    logical :: verbose_beta = .false.
!    logical :: verbose_beta = .true.
    logical :: verbose_efvs = .false.
!    logical :: verbose_efvs = .true.
    logical :: verbose_tau = .false.
!    logical :: verbose_tau = .true.
    logical :: verbose_gridop = .false.
!    logical :: verbose_gridop= .true.
    logical :: verbose_dirichlet = .false.
!    logical :: verbose_dirichlet= .true.
    logical :: verbose_L1L2 = .false.
!    logical :: verbose_L1L2 = .true.
    logical :: verbose_diva = .false.
!    logical :: verbose_diva = .true.
    logical :: verbose_glp = .false.
!    logical :: verbose_glp = .true.
    logical :: verbose_pcg = .false.
!    logical :: verbose_pcg = .true.

    integer :: itest, jtest    ! coordinates of diagnostic point
    integer :: rtest           ! task number for processor containing diagnostic point

    integer, parameter :: ktest = 1     ! vertical level of diagnostic point
    integer, parameter :: ptest = 1     ! diagnostic quadrature point

    ! option for writing matrix entries to text files
    logical, parameter :: write_matrix = .false.
!    logical, parameter :: write_matrix = .true.
    character(*), parameter :: matrix_label = 'label_here'  ! choose an appropriate label

    !WHL - debug for efvs
    real(dp), dimension(nNodesPerElement_3d, nQuadPoints_2d) ::   & 
       phi_3d_vav,         &! vertical avg of phi_3d
       dphi_dxr_3d_vav,    &! vertical avg of dphi_dxr_3d
       dphi_dyr_3d_vav,    &! vertical avg of dphi_dyr_3d
       dphi_dzr_3d_vav      ! vertical avg of dphi_dzr_3d

    contains

!****************************************************************************

  subroutine glissade_velo_higher_init

    !----------------------------------------------------------------
    ! Initial calculations for glissade higher-order solver.
    !----------------------------------------------------------------

    integer :: i, j, k, m, n, p
    integer :: pplus
    real(dp) :: xctr, yctr, zctr
    real(dp) :: sumx, sumy, sumz

    !----------------------------------------------------------------
    ! Initialize some time-independent finite element arrays
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Trilinear basis set for reference hexahedron, x=(-1,1), y=(-1,1), z=(-1,1)             
    ! Indexing is counter-clockwise from SW corner, with 1-4 on lower surface
    !  and 5-8 on upper surface
    ! The code uses "phi_3d" to denote these basis functions. 
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

    xqp_3d(1) = -rsqrt3; yqp_3d(1) = -rsqrt3; zqp_3d(1) = -rsqrt3
    wqp_3d(1) =  1.d0

    xqp_3d(2) =  rsqrt3; yqp_3d(2) = -rsqrt3; zqp_3d(2) = -rsqrt3
    wqp_3d(2) =  1.d0

    xqp_3d(3) =  rsqrt3; yqp_3d(3) =  rsqrt3; zqp_3d(3) = -rsqrt3
    wqp_3d(3) =  1.d0

    xqp_3d(4) = -rsqrt3; yqp_3d(4) =  rsqrt3; zqp_3d(4) = -rsqrt3
    wqp_3d(4) =  1.d0

    xqp_3d(5) = -rsqrt3; yqp_3d(5) = -rsqrt3; zqp_3d(5) =  rsqrt3
    wqp_3d(5) =  1.d0

    xqp_3d(6) =  rsqrt3; yqp_3d(6) = -rsqrt3; zqp_3d(6) =  rsqrt3
    wqp_3d(6) =  1.d0

    xqp_3d(7) =  rsqrt3; yqp_3d(7) =  rsqrt3; zqp_3d(7) =  rsqrt3
    wqp_3d(7) =  1.d0

    xqp_3d(8) = -rsqrt3; yqp_3d(8) =  rsqrt3; zqp_3d(8) =  rsqrt3
    wqp_3d(8) =  1.d0

    if (verbose_init) then
       print*, ' '
       print*, 'Hexahedral elements, quad points, x, y, z:'
       sumx = 0.d0; sumy = 0.d0; sumz = 0.d0
       do p = 1, nQuadPoints_3d
          print*, p, xqp_3d(p), yqp_3d(p), zqp_3d(p)
          sumx = sumx + xqp_3d(p); sumy = sumy + yqp_3d(p); sumz = sumz + zqp_3d(p)
       enddo
       print*, ' '
       print*, 'sums:', sumx, sumy, sumz
    endif

    ! Evaluate trilinear basis functions and their derivatives at each quad pt

    do p = 1, nQuadPoints_3d

       phi_3d(1,p) = (1.d0 - xqp_3d(p)) * (1.d0 - yqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0
       phi_3d(2,p) = (1.d0 + xqp_3d(p)) * (1.d0 - yqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0
       phi_3d(3,p) = (1.d0 + xqp_3d(p)) * (1.d0 + yqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0
       phi_3d(4,p) = (1.d0 - xqp_3d(p)) * (1.d0 + yqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0
       phi_3d(5,p) = (1.d0 - xqp_3d(p)) * (1.d0 - yqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0
       phi_3d(6,p) = (1.d0 + xqp_3d(p)) * (1.d0 - yqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0
       phi_3d(7,p) = (1.d0 + xqp_3d(p)) * (1.d0 + yqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0
       phi_3d(8,p) = (1.d0 - xqp_3d(p)) * (1.d0 + yqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0

       dphi_dxr_3d(1,p) = -(1.d0 - yqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0 
       dphi_dxr_3d(2,p) =  (1.d0 - yqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0 
       dphi_dxr_3d(3,p) =  (1.d0 + yqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0 
       dphi_dxr_3d(4,p) = -(1.d0 + yqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0
       dphi_dxr_3d(5,p) = -(1.d0 - yqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0 
       dphi_dxr_3d(6,p) =  (1.d0 - yqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0 
       dphi_dxr_3d(7,p) =  (1.d0 + yqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0 
       dphi_dxr_3d(8,p) = -(1.d0 + yqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0

       dphi_dyr_3d(1,p) = -(1.d0 - xqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0 
       dphi_dyr_3d(2,p) = -(1.d0 + xqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0 
       dphi_dyr_3d(3,p) =  (1.d0 + xqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0 
       dphi_dyr_3d(4,p) =  (1.d0 - xqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0 
       dphi_dyr_3d(5,p) = -(1.d0 - xqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0 
       dphi_dyr_3d(6,p) = -(1.d0 + xqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0 
       dphi_dyr_3d(7,p) =  (1.d0 + xqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0 
       dphi_dyr_3d(8,p) =  (1.d0 - xqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0 

       dphi_dzr_3d(1,p) = -(1.d0 - xqp_3d(p)) * (1.d0 - yqp_3d(p)) / 8.d0 
       dphi_dzr_3d(2,p) = -(1.d0 + xqp_3d(p)) * (1.d0 - yqp_3d(p)) / 8.d0 
       dphi_dzr_3d(3,p) = -(1.d0 + xqp_3d(p)) * (1.d0 + yqp_3d(p)) / 8.d0 
       dphi_dzr_3d(4,p) = -(1.d0 - xqp_3d(p)) * (1.d0 + yqp_3d(p)) / 8.d0 
       dphi_dzr_3d(5,p) =  (1.d0 - xqp_3d(p)) * (1.d0 - yqp_3d(p)) / 8.d0 
       dphi_dzr_3d(6,p) =  (1.d0 + xqp_3d(p)) * (1.d0 - yqp_3d(p)) / 8.d0 
       dphi_dzr_3d(7,p) =  (1.d0 + xqp_3d(p)) * (1.d0 + yqp_3d(p)) / 8.d0 
       dphi_dzr_3d(8,p) =  (1.d0 - xqp_3d(p)) * (1.d0 + yqp_3d(p)) / 8.d0 

       if (verbose_init) then
          print*, ' '
          print*, 'Quad point, p =', p
          print*, 'n, phi_3d, dphi_dxr_3d, dphi_dyr_3d, dphi_dzr_3d:'
          do n = 1, 8
             print*, n, phi_3d(n,p), dphi_dxr_3d(n,p), dphi_dyr_3d(n,p), dphi_dzr_3d(n,p)
          enddo
          print*, ' '
          print*, 'sum(phi_3d)', sum(phi_3d(:,p))  ! verified that sum = 1
          print*, 'sum(dphi/dx)', sum(dphi_dxr_3d(:,p))  ! verified that sum = 0 (within roundoff)
          print*, 'sum(dphi/dy)', sum(dphi_dyr_3d(:,p))  ! verified that sum = 0 (within roundoff)
          print*, 'sum(dphi/dz)', sum(dphi_dzr_3d(:,p))  ! verified that sum = 0 (within roundoff)
       endif

    enddo   ! nQuadPoints_3d

    ! Evaluate trilinear basis functions and their derivatives at cell center
    ! Full formulas are not really needed at (x,y,z) = (0,0,0), but are included for completeness

    xctr = 0.d0
    yctr = 0.d0
    zctr = 0.d0

    phi_3d_ctr(1) = (1.d0 - xctr) * (1.d0 - yctr) * (1.d0 - zctr) / 8.d0
    phi_3d_ctr(2) = (1.d0 + xctr) * (1.d0 - yctr) * (1.d0 - zctr) / 8.d0
    phi_3d_ctr(3) = (1.d0 + xctr) * (1.d0 + yctr) * (1.d0 - zctr) / 8.d0
    phi_3d_ctr(4) = (1.d0 - xctr) * (1.d0 + yctr) * (1.d0 - zctr) / 8.d0
    phi_3d_ctr(5) = (1.d0 - xctr) * (1.d0 - yctr) * (1.d0 + zctr) / 8.d0
    phi_3d_ctr(6) = (1.d0 + xctr) * (1.d0 - yctr) * (1.d0 + zctr) / 8.d0
    phi_3d_ctr(7) = (1.d0 + xctr) * (1.d0 + yctr) * (1.d0 + zctr) / 8.d0
    phi_3d_ctr(8) = (1.d0 - xctr) * (1.d0 + yctr) * (1.d0 + zctr) / 8.d0
    
    dphi_dxr_3d_ctr(1) = -(1.d0 - yctr) * (1.d0 - zctr) / 8.d0 
    dphi_dxr_3d_ctr(2) =  (1.d0 - yctr) * (1.d0 - zctr) / 8.d0 
    dphi_dxr_3d_ctr(3) =  (1.d0 + yctr) * (1.d0 - zctr) / 8.d0 
    dphi_dxr_3d_ctr(4) = -(1.d0 + yctr) * (1.d0 - zctr) / 8.d0
    dphi_dxr_3d_ctr(5) = -(1.d0 - yctr) * (1.d0 + zctr) / 8.d0 
    dphi_dxr_3d_ctr(6) =  (1.d0 - yctr) * (1.d0 + zctr) / 8.d0 
    dphi_dxr_3d_ctr(7) =  (1.d0 + yctr) * (1.d0 + zctr) / 8.d0 
    dphi_dxr_3d_ctr(8) = -(1.d0 + yctr) * (1.d0 + zctr) / 8.d0
    
    dphi_dyr_3d_ctr(1) = -(1.d0 - xctr) * (1.d0 - zctr) / 8.d0 
    dphi_dyr_3d_ctr(2) = -(1.d0 + xctr) * (1.d0 - zctr) / 8.d0 
    dphi_dyr_3d_ctr(3) =  (1.d0 + xctr) * (1.d0 - zctr) / 8.d0 
    dphi_dyr_3d_ctr(4) =  (1.d0 - xctr) * (1.d0 - zctr) / 8.d0 
    dphi_dyr_3d_ctr(5) = -(1.d0 - xctr) * (1.d0 + zctr) / 8.d0 
    dphi_dyr_3d_ctr(6) = -(1.d0 + xctr) * (1.d0 + zctr) / 8.d0 
    dphi_dyr_3d_ctr(7) =  (1.d0 + xctr) * (1.d0 + zctr) / 8.d0 
    dphi_dyr_3d_ctr(8) =  (1.d0 - xctr) * (1.d0 + zctr) / 8.d0 
    
    dphi_dzr_3d_ctr(1) = -(1.d0 - xctr) * (1.d0 - yctr) / 8.d0 
    dphi_dzr_3d_ctr(2) = -(1.d0 + xctr) * (1.d0 - yctr) / 8.d0 
    dphi_dzr_3d_ctr(3) = -(1.d0 + xctr) * (1.d0 + yctr) / 8.d0 
    dphi_dzr_3d_ctr(4) = -(1.d0 - xctr) * (1.d0 + yctr) / 8.d0 
    dphi_dzr_3d_ctr(5) =  (1.d0 - xctr) * (1.d0 - yctr) / 8.d0 
    dphi_dzr_3d_ctr(6) =  (1.d0 + xctr) * (1.d0 - yctr) / 8.d0 
    dphi_dzr_3d_ctr(7) =  (1.d0 + xctr) * (1.d0 + yctr) / 8.d0 
    dphi_dzr_3d_ctr(8) =  (1.d0 - xctr) * (1.d0 + yctr) / 8.d0 

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
    ! The code uses "phi_2d" to denote these basis functions. 
    !
    ! N1 = (1-x)*(1-y)/4             N4----N3
    ! N2 = (1+x)*(1-y)/4             |     |
    ! N3 = (1+x)*(1+y)/4             |     |
    ! N4 = (1-x)*(1+y)/4             N1----N2
    !----------------------------------------------------------------

    ! Set coordinates and weights of quadrature points for reference square.
    ! Numbering is counter-clockwise from southwest

    xqp_2d(1) = -rsqrt3; yqp_2d(1) = -rsqrt3
    wqp_2d(1) =  1.d0

    xqp_2d(2) =  rsqrt3; yqp_2d(2) = -rsqrt3
    wqp_2d(2) =  1.d0

    xqp_2d(3) =  rsqrt3; yqp_2d(3) =  rsqrt3
    wqp_2d(3) =  1.d0

    xqp_2d(4) = -rsqrt3; yqp_2d(4) =  rsqrt3
    wqp_2d(4) =  1.d0

    if (verbose_init) then
       print*, ' '
       print*, ' '
       print*, 'Quadrilateral elements, quad points, x, y:'
       sumx = 0.d0; sumy = 0.d0; sumz = 0.d0
       do p = 1, nQuadPoints_2d
          print*, p, xqp_2d(p), yqp_2d(p)
          sumx = sumx + xqp_2d(p); sumy = sumy + yqp_2d(p)
       enddo
       print*, ' '
       print*, 'sumx, sumy:', sumx, sumy
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
          print*, 'n, phi_2d, dphi_dxr_2d, dphi_dyr_2d:'
          do n = 1, 4
             print*, n, phi_2d(n,p), dphi_dxr_2d(n,p), dphi_dyr_2d(n,p)
          enddo
          print*, 'sum(phi_2d)', sum(phi_2d(:,p))        ! verified that sum = 1
          print*, 'sum(dphi/dx_2d)', sum(dphi_dxr_2d(:,p))  ! verified that sum = 0 (within roundoff)
          print*, 'sum(dphi/dy_2d)', sum(dphi_dyr_2d(:,p))  ! verified that sum = 0 (within roundoff)
       endif

    enddo   ! nQuadPoints_2d

    ! Evaluate bilinear basis functions and their derivatives at cell center
    ! Full formulas are not really needed at (x,y) = (0,0), but are included for completeness

    xctr = 0.d0
    yctr = 0.d0

    phi_2d_ctr(1) = (1.d0 - xctr) * (1.d0 - yctr) / 4.d0 
    phi_2d_ctr(2) = (1.d0 + xctr) * (1.d0 - yctr) / 4.d0
    phi_2d_ctr(3) = (1.d0 + xctr) * (1.d0 + yctr) / 4.d0 
    phi_2d_ctr(4) = (1.d0 - xctr) * (1.d0 + yctr) / 4.d0
    
    dphi_dxr_2d_ctr(1) = -(1.d0 - yctr) / 4.d0 
    dphi_dxr_2d_ctr(2) =  (1.d0 - yctr) / 4.d0 
    dphi_dxr_2d_ctr(3) =  (1.d0 + yctr) / 4.d0 
    dphi_dxr_2d_ctr(4) = -(1.d0 + yctr) / 4.d0

    dphi_dyr_2d_ctr(1) = -(1.d0 - xctr) / 4.d0 
    dphi_dyr_2d_ctr(2) = -(1.d0 + xctr) / 4.d0 
    dphi_dyr_2d_ctr(3) =  (1.d0 + xctr) / 4.d0 
    dphi_dyr_2d_ctr(4) =  (1.d0 - xctr) / 4.d0 

    !----------------------------------------------------------------
    ! Compute indxA_3d; maps displacements i,j,k = (-1,0,1) onto an index from 1 to 27
    ! Numbering starts in SW corner of layers k-1, finishes in NE corner of layer k+1
    ! Diagonal term has index 14
    !----------------------------------------------------------------

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
             indxA_3d(i,j,k) = m
          enddo
       enddo
    enddo

    !----------------------------------------------------------------
    ! Compute indxA_2d; maps displacements i,j = (-1,0,1) onto an index from 1 to 9
    ! Same as indxA_3d, but for a single layer
    !----------------------------------------------------------------

    m = 0
    do j = -1,1
       do i = -1,1
          m = m + 1
          indxA_2d(i,j) = m
       enddo
    enddo

    !WHL - debug for efvs

    ! Evaluate vertical averages of dphi_dxr_3d, dphi_dyr_3d and dphi_dzr_3d at each 2d quad pts.
    ! Using these instead of the full 3d basis functions can result in similar accuracy with
    !  only half as many QP computations.

    do p = 1, nQuadPoints_2d
       pplus = p + nQuadPoints_3d/2  ! p + 4 for hexahedra
       do n = 1, nNodesPerElement_3d
          phi_3d_vav(n,p) = 0.5d0 * (phi_3d(n,p) + phi_3d(n,pplus))
          dphi_dxr_3d_vav(n,p) = 0.5d0 * (dphi_dxr_3d(n,p) + dphi_dxr_3d(n,pplus))
          dphi_dyr_3d_vav(n,p) = 0.5d0 * (dphi_dyr_3d(n,p) + dphi_dyr_3d(n,pplus))
          dphi_dzr_3d_vav(n,p) = 0.5d0 * (dphi_dzr_3d(n,p) + dphi_dzr_3d(n,pplus))
       enddo
    enddo

  end subroutine glissade_velo_higher_init

!****************************************************************************

  subroutine glissade_velo_higher_solve(model,                &
                                        nx,     ny,     nz)

    !TODO - Remove nx, ny, nz from argument list?
    !       Would then have to allocate many local arrays.

    !----------------------------------------------------------------
    ! Solve the ice sheet flow equations for the horizontal velocity (uvel, vvel)
    !  at each node of each grid cell where ice is present.
    ! The standard solver is based on the Blatter-Pattyn first-order approximation
    !  of Stokes flow (which_ho_approx = HO_APPROX_BP).
    ! There are also options to solve the shallow-ice equations (HO_APPROX_SIA),
    !  shallow-shelf equations (HO_APPROX_SIA), or L1L2 equations (HO_APPROX_L1L2).
    ! Note: The SIA solver does a full matrix solution and is much slower than
    !       the local SIA solver (HO_APPROX_LOCAL_SIA) in glissade_velo_sia.F90.
    !----------------------------------------------------------------

    use glissade_basal_traction, only: calcbeta
    use glissade_therm, only: glissade_pressure_melting_point

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    type(glide_global_type), intent(inout) :: model   ! derived type holding ice-sheet info

    !----------------------------------------------------------------
    ! Note that the glissade solver uses SI units.
    ! Thus we have grid cell dimensions and ice thickness in meters,
    !  velocity in m/s, and the rate factor in Pa^(-n) s(-1).
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Note: nx and ny are the horizontal dimensions of scalar arrays (e.g., thck and temp).
    !       The velocity arrays have horizontal dimensions (nx-1, ny-1).
    !       nz is the number of levels at which uvel and vvel are computed.
    !       The scalar variables generally live at layer midpoints and have
    !         vertical dimension nz-1.
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! number of grid cells in each horizontal direction
       nz                       ! number of vertical levels where velocity is computed
                                ! (same as model%general%upn)
 
    !----------------------------------------------------------------
    ! Local variables and pointers set to components of model derived type 
    !----------------------------------------------------------------

    real(dp) ::  &
       dx,  dy                  ! grid cell length and width (m)
                                ! assumed to have the same value for each grid cell

    real(dp), dimension(:), pointer :: &
       sigma,                 & ! vertical sigma coordinate at layer interfaces, [0,1]
       stagsigma,             & ! staggered vertical sigma coordinate at layer midpoints
       stagwbndsigma            ! stagsigma augmented by sigma = 0 and 1 at upper and lower surfaces

    real(dp)  ::   & 
       thklim,               &  ! minimum ice thickness for active cells (m)
       max_slope,            &  ! maximum slope allowed for surface gradient computations (unitless)
       eus,                  &  ! eustatic sea level (m), = 0. by default
       ho_beta_const,        &  ! constant beta value (Pa/(m/yr)) for whichbabc = HO_BABC_CONSTANT
       beta_grounded_min,    &  ! minimum beta value (Pa/(m/yr)) for grounded ice
       efvs_constant            ! constant efvs value (Pa yr) for whichefvs = HO_EFVS_CONSTANT

    real(dp), dimension(:,:), pointer ::  &
       thck,                 &  ! ice thickness (m)
       usrf,                 &  ! upper surface elevation (m)
       topg,                 &  ! elevation of topography (m)
       bwat,                 &  ! basal water depth (m)
       mintauf,              &  ! till yield stress (Pa)
       beta,                 &  ! basal traction parameter (Pa/(m/yr))
       beta_internal,        &  ! beta field weighted by f_ground (such that beta = 0 beneath floating ice)
       bfricflx,             &  ! basal heat flux from friction (W/m^2) 
       f_flotation,          &  ! flotation function = (rhoi*thck) / (-rhoo*(topg-eus)) by default
                                ! used to be f_pattyn = -rhoo*(topg-eus) / (rhoi*thck)
       f_ground                 ! grounded ice fraction, 0 <= f_ground <= 1

    !TODO - Remove dependence on stagmask?  Currently it is needed for input to calcbeta.
    integer, dimension(:,:), pointer ::   &
       stagmask                 ! mask on staggered grid

    real(dp), dimension(:,:,:), pointer ::  &
       uvel, vvel,  &           ! velocity components (m/yr)
       temp,   &                ! ice temperature (deg C)
       flwa,   &                ! flow factor in units of Pa^(-n) yr^(-1)
       efvs,   &                ! effective viscosity (Pa yr)
       resid_u, resid_v,   &    ! u and v components of residual Ax - b (Pa/m)
       bu, bv                   ! right-hand-side vector b, divided into 2 parts

    real(dp), dimension(:,:), pointer ::  &
       uvel_2d, vvel_2d,       &! 2D velocity field; solution for SSA, L1L2 and DIVA 
       btractx, btracty,       &! components of basal traction (Pa)
       taudx, taudy             ! components of driving stress (Pa)

    real(dp), dimension(:,:,:), pointer ::  &
       tau_xz, tau_yz,         &! vertical components of stress tensor (Pa)
       tau_xx, tau_yy, tau_xy, &! horizontal components of stress tensor (Pa)
       tau_eff                  ! effective stress (Pa)

    integer,  dimension(:,:), pointer ::   &
       kinbcmask,              &! = 1 at vertices where u and v are prescribed from input data (Dirichlet BC), = 0 elsewhere
       umask_no_penetration,   &! = 1 at vertices along east/west global boundary where uvel = 0, = 0 elsewhere
       vmask_no_penetration     ! = 1 at vertices along north/south global boundary where vvel = 0, = 0 elsewhere

    integer ::   &
       whichbabc, &             ! option for basal boundary condition
       whichefvs, &             ! option for effective viscosity calculation 
                                ! (calculate it or make it uniform)
       whichresid, &            ! option for method of calculating residual
       whichsparse, &           ! option for method of doing elliptic solve
                                ! (BiCG, GMRES, standalone Trilinos, etc.)
       whichapprox, &           ! option for which Stokes approximation to use
                                ! 0 = SIA, 1 = SSA, 2 = Blatter-Pattyn HO, 3 = L1L2
                                ! default = 2
       whichprecond, &          ! option for which preconditioner to use with 
                                !  structured PCG solver
                                ! 0 = none, 1 = diag, 2 = SIA-based
       whichgradient, &         ! option for gradient operator when computing grad(s)
                                ! 0 = centered, 1 = upstream
       whichgradient_margin, &  ! option for computing gradient at ice margin
                                ! 0 = include all neighbor cells in gradient calculation
                                ! 1 = include ice-covered and/or land cells
                                ! 2 = include ice-covered cells only
       whichassemble_beta,  &   ! 0 = standard finite element assembly
                                ! 1 = apply local value of beta at each vertex
       whichassemble_taud,  &   ! 0 = standard finite element assembly
                                ! 1 = apply local value of driving stress at each vertex
       whichassemble_bfric, &   ! 0 = standard finite element assembly
                                ! 1 = apply local value of basal friction at each vertex
       whichground,  &          ! option for computing grounded fraction of each cell
       whichflotation_function,&! option for computing flotation function at and near each vertex
       maxiter_nonlinear        ! maximum number of nonlinear iterations

    !--------------------------------------------------------
    ! Local parameters
    !--------------------------------------------------------

    real(dp), parameter :: resid_target = 1.0d-04   ! assume velocity fields have converged below this resid 

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    real(dp), dimension(nx-1,ny-1) :: &
       xVertex, yVertex,    & ! x and y coordinates of each vertex (m)
       stagusrf,            & ! upper surface averaged to vertices (m)
       stagthck,            & ! ice thickness averaged to vertices (m)
       dusrf_dx, dusrf_dy,  & ! gradient of upper surface elevation (m/m)
       ubas, vbas             ! basal ice velocity (m/yr); input to calcbeta 

    integer, dimension(nx,ny) ::     &
       ice_mask,             &! = 1 for cells where ice is present (thk > thklim), else = 0
       floating_mask,        &! = 1 for cells where ice is present and is floating
       ocean_mask,           &! = 1 for cells where topography is below sea level and ice is absent
       land_mask              ! = 1 for cells where topography is above sea level

    real(dp), dimension(nx,ny) ::  &
       bedpmp                 ! basal pressure melting point temperature (deg C)

    real(dp), dimension(nx-1,ny-1) :: &
       stagbedtemp,         & ! basal ice temperature averaged to vertices (deg C)
       stagbedpmp             ! bedpmp averaged to vertices (deg C)    

    integer, dimension(nx-1,ny-1) :: &
       pmp_mask               ! = 1 where bed is at pressure melting point, elsewhere = 0

    logical, dimension(nx,ny) ::     &
       active_cell            ! true for active cells (thck > thklim and border locally owned vertices)

    logical, dimension(nx-1,ny-1) :: &
       active_vertex          ! true for vertices of active cells

    real(dp), dimension(nz-1,nx,ny) ::  &
       flwafact               ! temperature-based flow factor, 0.5 * A^(-1/n), 
                              ! used to compute effective viscosity
                              ! units: Pa yr^(1/n)

    real(dp), dimension(nz,nx-1,ny-1) ::   &
       usav, vsav,                 &! previous guess for velocity solution
       loadu, loadv                 ! assembled load vector, divided into 2 parts
                                    ! Note: loadu and loadv are computed only once per nonlinear solve,
                                    !       whereas bu and bv can be set each nonlinear iteration to account 
                                    !       for inhomogeneous Dirichlet BC
  
    integer, dimension(nz,nx-1,ny-1) ::    &
       umask_dirichlet,     & ! Dirichlet mask for u component of velocity, = 1 for prescribed velo, else = 0
       vmask_dirichlet        ! Dirichlet mask for v component of velocity, = 1 for prescribed velo, else = 0

    real(dp) :: &
       resid_velo,          & ! quantity related to velocity convergence
       L2_norm,             & ! L2 norm of residual, |Ax - b|
       L2_target,           & ! nonlinear convergence target for residual
       L2_norm_relative,    & ! L2 norm of residual relative to rhs, |Ax - b| / |b|
       L2_target_relative,  & ! nonlinear convergence target for relative residual
       err,                 & ! solution error from sparse_easy_solve
       outer_it_criterion,  & ! current value of outer (nonlinear) loop converence criterion
       outer_it_target        ! target value for outer-loop convergence

    logical, save ::    &
       converged_soln = .false.    ! true if we get a converged solution for velocity

    integer ::  & 
       counter,         & ! outer (nonlinear) iteration counter
       niters             ! linear iteration count

    integer :: nNonzeros  ! number of nonzero matrix entries

    ! The following large matrix arrays are allocated for a 3D solve (SIA or BP)

    real(dp), dimension(:,:,:,:), allocatable ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv           ! 1st dimension = 27 (node and its nearest neighbors in x, y and z direction) 
                          ! other dimensions = (k,i,j)

    ! The following are used for the SLAP and Trilinos solvers

    integer ::            &
       nNodesSolve        ! number of nodes where we solve for velocity

    integer, dimension(nz,nx-1,ny-1) ::  &
       nodeID             ! local ID for each node where we solve for velocity
                          ! For periodic BCs (as in ISMIP-HOM), halo node IDs will be copied
                          !  from the other side of the grid

    integer, dimension((nx-1)*(ny-1)*nz) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of nodes

    ! The following are used for the Trilinos solver only

    integer, dimension(nx-1,ny-1) ::  &
       global_vertex_id    ! unique global IDs for vertices on this processor

    integer, dimension(nz,nx-1,ny-1) ::  &
       global_node_id      ! unique global IDs for nodes on this processor

    integer, dimension(:), allocatable ::    &
       active_owned_unknown_map    ! maps owned active unknowns (u and v at each active node) to global IDs

    logical, dimension(:,:,:,:), allocatable ::  &
       Afill               ! true wherever the matrix value is potentially nonzero

    real(dp), dimension(:), allocatable ::   &
       velocityResult     ! velocity solution vector from Trilinos

    ! The following are used for the SLAP solver only

    type(sparse_matrix_type) ::  &
       matrix             ! sparse matrix for SLAP solver, defined in glimmer_sparse_types
                          ! includes nonzeroes, order, col, row, val 

    real(dp), dimension(:), allocatable ::   &   ! for SLAP solver
       rhs,             & ! right-hand-side (b) in Ax = b
       answer,          & ! answer (x) in Ax = b
       resid_vec          ! residual vector Ax - b

    integer ::          &
       matrix_order,    & ! order of matrix = number of rows
       max_nonzeros       ! upper bound for number of nonzero entries in sparse matrix

    ! The following arrays are used for a 2D matrix solve (SSA, L1L2 or DIVA)

    logical ::  &
       solve_2d           ! if true, solve a 2D matrix)
                          ! else solve a 3D matrix (SIA, BP)

    integer ::            &
       nVerticesSolve     ! number of vertices where we solve for velocity

    integer, dimension(nx-1,ny-1) ::  &
       vertexID           ! local ID for each vertex where we solve for velocity (in 2d)
    
    integer, dimension((nx-1)*(ny-1)) ::   &
       iVertexIndex, jVertexIndex   ! i and j indices of vertices

    real(dp), dimension(:,:,:), allocatable ::  &
       Auu_2d, Auv_2d,   &! assembled stiffness matrix, divided into 4 parts
       Avu_2d, Avv_2d     ! 1st dimension = 9 (node and its nearest neighbors in x and y direction) 
                          ! other dimensions = (i,j)

    real(dp), dimension(:,:), allocatable ::  &
       bu_2d, bv_2d,     &! right-hand-side vector b, divided into 2 parts
       loadu_2d, loadv_2d ! assembled load vector, divided into 2 parts

    real(dp), dimension(:,:), allocatable ::  &
       usav_2d, vsav_2d

    real(dp), dimension(:,:), allocatable ::  &
       resid_u_2d, resid_v_2d   ! components of 2D solution residual

    logical, dimension(:,:,:), allocatable ::  &
       Afill_2d           ! true wherever the matrix value is potentially nonzero
                          ! 2D Trilinos only

    ! The following are used for the depth-integrated viscosity solve
    real(dp), dimension(:,:), allocatable :: &
       beta_eff,            &! effective beta, defined by Goldberg (2011) eq. 41
                             ! beta*u_b = beta_eff*u_av
       omega,               &! double integral, defined by Goldberg (2011) eq. 35
                             ! Note: omega here is equal to Goldberg's omega/H
       stag_omega            ! omega interpolated to staggered grid

    real(dp), dimension(:,:,:), allocatable :: &
       omega_k,             &! single integral, defined by Goldberg (2011) eq. 32
       stag_omega_k          ! omega_k interpolated to staggered grid

    real(dp), dimension(:,:,:,:), allocatable :: &
       efvs_qp_3d            ! effective viscosity at each QP of each layer of each cell

    integer, parameter :: &
       diva_level_index = 0  ! level for which the DIVA scheme finds the 2D velocity
                             ! 0 = mean, 1 = upper surface
                             ! Results are not very sensitive to this choice                     
    real(dp) :: dsigma
    real(dp) :: maxbeta, minbeta
    integer :: i, j, k, m, n, p, r
    integer :: iA, jA, kA
    real(dp) :: maxthck, maxusrf
    logical, parameter :: test_matrix = .false.
!    logical, parameter :: test_matrix = .true.
    integer, parameter :: test_order = 4

    ! for trilinos test problem
    logical, parameter :: test_trilinos = .false.
!    logical, parameter :: test_trilinos = .true.

    ! for diagnostic prints
    integer, parameter :: xmax_print = 20

    call t_startf('glissade_vhs_init')
    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    if (verbose .and. this_rank==rtest) then
       print*, 'In glissade_velo_higher_solve'
       print*, 'rank, itest, jtest, ktest =', rtest, itest, jtest, ktest
    endif

#ifdef TRILINOS
    if (test_trilinos) then
       call trilinos_test
       stop
    endif
#endif

    !--------------------------------------------------------
    ! Assign local pointers and variables to derived type components
    !--------------------------------------------------------

!    nx = model%general%ewn   ! currently passed in
!    ny = model%general%nsn
!    nz = model%general%upn

     dx = model%numerics%dew
     dy = model%numerics%dns

     !TODO - Remove (:), (:,:) and (:,:,:) from pointer targets?
     sigma    => model%numerics%sigma(:)
     stagsigma=> model%numerics%stagsigma(:)
     stagwbndsigma=> model%numerics%stagwbndsigma(:)
     thck     => model%geometry%thck(:,:)
     usrf     => model%geometry%usrf(:,:)
     topg     => model%geometry%topg(:,:)
     stagmask => model%geometry%stagmask(:,:)
     f_ground => model%geometry%f_ground(:,:)
     f_flotation => model%geometry%f_flotation(:,:)

     temp     => model%temper%temp
     flwa     => model%temper%flwa(:,:,:)
     efvs     => model%stress%efvs(:,:,:)
     beta     => model%velocity%beta(:,:)
     beta_internal => model%velocity%beta_internal(:,:)
     bfricflx => model%temper%bfricflx(:,:)
     bwat     => model%temper%bwat(:,:)
     mintauf  => model%basalproc%mintauf(:,:)

     uvel     => model%velocity%uvel(:,:,:)
     vvel     => model%velocity%vvel(:,:,:)
     uvel_2d  => model%velocity%uvel_2d(:,:)
     vvel_2d  => model%velocity%vvel_2d(:,:)
     resid_u  => model%velocity%resid_u(:,:,:)
     resid_v  => model%velocity%resid_v(:,:,:)
     bu       => model%velocity%rhs_u(:,:,:)
     bv       => model%velocity%rhs_v(:,:,:)

     btractx  => model%stress%btractx(:,:)
     btracty  => model%stress%btracty(:,:)
     taudx    => model%stress%taudx(:,:)
     taudy    => model%stress%taudy(:,:)
     tau_xz   => model%stress%tau%xz(:,:,:)
     tau_yz   => model%stress%tau%yz(:,:,:)
     tau_xx   => model%stress%tau%xx(:,:,:)
     tau_yy   => model%stress%tau%yy(:,:,:)
     tau_xy   => model%stress%tau%xy(:,:,:)
     tau_eff  => model%stress%tau%scalar(:,:,:)

     kinbcmask => model%velocity%kinbcmask(:,:)
     umask_no_penetration => model%velocity%umask_no_penetration(:,:)
     vmask_no_penetration => model%velocity%vmask_no_penetration(:,:)

     thklim    = model%numerics%thklim
     max_slope = model%paramets%max_slope
     eus       = model%climate%eus
     ho_beta_const = model%velocity%ho_beta_const
     beta_grounded_min = model%velocity%beta_grounded_min
     efvs_constant = model%paramets%efvs_constant
 
     whichbabc            = model%options%which_ho_babc
     whichefvs            = model%options%which_ho_efvs
     whichresid           = model%options%which_ho_resid
     whichsparse          = model%options%which_ho_sparse
     whichapprox          = model%options%which_ho_approx
     whichprecond         = model%options%which_ho_precond
     maxiter_nonlinear    = model%options%glissade_maxiter
     whichgradient        = model%options%which_ho_gradient
     whichgradient_margin = model%options%which_ho_gradient_margin
     whichassemble_beta   = model%options%which_ho_assemble_beta
     whichassemble_taud   = model%options%which_ho_assemble_taud
     whichassemble_bfric  = model%options%which_ho_assemble_bfric
     whichground          = model%options%which_ho_ground
     whichflotation_function = model%options%which_ho_flotation_function

    !--------------------------------------------------------
    ! Convert input variables to appropriate units for this solver.
    ! (Mainly SI, except that time units in flwa, velocities,
    !  and beta are years instead of seconds)
    !--------------------------------------------------------

!pw call t_startf('glissade_velo_higher_scale_input')
    call glissade_velo_higher_scale_input(dx,      dy,            &
                                          thck,    usrf,          &
                                          topg,                   &
                                          eus,     thklim,        &
                                          flwa,    efvs,          &
                                          bwat,    mintauf,       &
                                          ho_beta_const,          &
                                          beta_grounded_min,      &
                                          btractx, btracty,       &
                                          uvel,    vvel,          &
                                          uvel_2d, vvel_2d)
!pw call t_stopf('glissade_velo_higher_scale_input')

    ! Set volume scale
    ! This is not strictly necessary, but dividing by this scale gives matrix coefficients 
    !  that are ~1.

    vol0  = 1.0d9    ! volume scale (m^3)

    if (whichapprox == HO_APPROX_SIA) then   ! SIA
!!       if (verbose .and. main_task) print*, 'Solving shallow-ice approximation'
       if (main_task) print*, 'Solving shallow-ice approximation'

    elseif (whichapprox == HO_APPROX_SSA) then  ! SSA
!!       if (verbose .and. main_task) print*, 'Solving shallow-shelf approximation'
       if (main_task) print*, 'Solving shallow-shelf approximation'

    elseif (whichapprox == HO_APPROX_L1L2) then  ! L1L2
!!       if (verbose .and. main_task) print*, 'Solving depth-integrated L1L2 approximation'
       if (main_task) print*, 'Solving depth-integrated L1L2 approximation'

    elseif (whichapprox == HO_APPROX_DIVA) then  ! DIVA, based on Goldberg (2011)
!!       if (verbose .and. main_task) print*, 'Solving depth-integrated viscosity approximation'
       if (main_task) print*, 'Solving depth-integrated viscosity approximation'

    else   ! Blatter-Pattyn higher-order 
!!       if (verbose .and. main_task) print*, 'Solving Blatter-Pattyn higher-order approximation'
       if (main_task) print*, 'Solving Blatter-Pattyn higher-order approximation'
    endif

    if (whichapprox==HO_APPROX_SSA .or. whichapprox==HO_APPROX_L1L2 .or. whichapprox==HO_APPROX_DIVA) then
       solve_2d = .true.
    else   ! 3D solve
       solve_2d = .false.
    endif

    if (solve_2d) then
       ! allocate arrays needed for a 2D solve
       allocate(Auu_2d(nNodeNeighbors_2d,nx-1,ny-1))
       allocate(Auv_2d(nNodeNeighbors_2d,nx-1,ny-1))
       allocate(Avu_2d(nNodeNeighbors_2d,nx-1,ny-1))
       allocate(Avv_2d(nNodeNeighbors_2d,nx-1,ny-1))
       allocate(bu_2d(nx-1,ny-1))
       allocate(bv_2d(nx-1,ny-1))
       allocate(loadu_2d(nx-1,ny-1))
       allocate(loadv_2d(nx-1,ny-1))
       allocate(usav_2d(nx-1,ny-1))
       allocate(vsav_2d(nx-1,ny-1))
       allocate(resid_u_2d(nx-1,ny-1))
       allocate(resid_v_2d(nx-1,ny-1))
    else
       ! These are big, so do not allocate them for the 2D solve
       allocate(Auu(nNodeNeighbors_3d,nz,nx-1,ny-1))
       allocate(Auv(nNodeNeighbors_3d,nz,nx-1,ny-1))
       allocate(Avu(nNodeNeighbors_3d,nz,nx-1,ny-1))
       allocate(Avv(nNodeNeighbors_3d,nz,nx-1,ny-1))
    endif

    if (whichapprox == HO_APPROX_DIVA) then
       allocate(beta_eff(nx-1,ny-1))
       allocate(omega(nx,ny))
       allocate(omega_k(nz,nx,ny))
       allocate(stag_omega(nx-1,ny-1))
       allocate(stag_omega_k(nz,nx-1,ny-1))
       allocate(efvs_qp_3d(nz-1,nQuadPoints_2d,nx,ny))
       beta_eff(:,:) = 0.d0
       omega(:,:) = 0.d0
       omega_k(:,:,:) = 0.d0
       stag_omega(:,:) = 0.d0
       stag_omega_k(:,:,:) = 0.d0
       ! Note: Initializing efvs_qp as efvs is a reasonable first guess that allows us to
       !       write efvs to the restart file instead of efvs_qp (which is 4x larger).
       do p = 1, nQuadPoints_2d
          efvs_qp_3d(:,p,:,:) = efvs(:,:,:)
       enddo
    endif

    if (whichapprox /= HO_APPROX_DIVA) then
       ! Set the 2D velocity to the velocity at the bed
       ! Note: For L1L2 and SSA, this is the 2D velocity solution from the previous solve.
       !       For DIVA, the velocity solution from the previous solve is typically the
       !        mean velocity, which cannot be extracted exactly from the 3D velocity field
       !        and must be stored in a separate array.
       uvel_2d(:,:) = uvel(nz,:,:)
       vvel_2d(:,:) = vvel(nz,:,:)
    endif

    if (test_matrix) then
       if (whichsparse <= HO_SPARSE_GMRES) then   ! this test works for SLAP solver only
          call slap_solve_test_matrix(test_order, whichsparse)
       else
          print*, 'Invalid value for whichsparse with test_matrix subroutine'
          stop
       endif
    endif

    ! Make sure that the geometry and flow factor are correct in halo cells.
    ! These calls are commented out, since the halo updates are done in 
    !  module glissade.F90, before calling glissade_velo_higher_solve.

!    call parallel_halo(thck)
!    call parallel_halo(topg)
!    call parallel_halo(usrf)
!    call parallel_halo(flwa)

!    if (whichbabc == HO_BABC_YIELD_PICARD) then
!       call staggered_parallel_halo(mintauf)   
!       call staggered_parallel_halo_extrapolate(mintauf)
!    endif

    !------------------------------------------------------------------------------
    ! Setup for higher-order solver: Compute nodal geometry, allocate storage, etc.
    ! These are quantities that do not change during the outer nonlinear loop. 
    !------------------------------------------------------------------------------

    if (verbose_state) then
       maxthck = maxval(thck(:,:))
       maxthck = parallel_reduce_max(maxthck)
       maxusrf = maxval(usrf(:,:))
       maxusrf = parallel_reduce_max(maxusrf)

       if (this_rank==rtest) then

          print*, ' '
          print*, 'nx, ny, nz:', nx, ny, nz
          print*, 'vol0:', vol0
          print*, 'thklim:', thklim
          print*, 'max thck, usrf:', maxthck, maxusrf
          
          print*, 'sigma coordinate:'
          do k = 1, nz
             print*, k, sigma(k)
          enddo
          
          print*, ' '
          print*, 'Thickness field, rank =', rtest
          do j = ny, 1, -1
             do i = 1, nx
                write(6,'(f6.0)',advance='no') thck(i,j)
             enddo
             write(6,*) ' '
          enddo
          
          print*, ' '
          print*, 'Topography field, rank =', rtest
          do j = ny, 1, -1
             do i = 1, nx
                write(6,'(f6.0)',advance='no') topg(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          
          print*, 'Upper surface field, rank =', rtest
          do j = ny, 1, -1
             do i = 1, nx
                write(6,'(f6.0)',advance='no') usrf(i,j)
             enddo
             write(6,*) ' '
          enddo
          
          print*, ' '
          print*, 'flwa (Pa-3 yr-1), k = 1, rank =', rtest
          do j = ny, 1, -1
             do i = 1, nx
                write(6,'(e12.5)',advance='no') flwa(1,i,j)
             enddo
             write(6,*) ' '
          enddo

       endif   ! this_rank
    endif      ! verbose_state
 
    !------------------------------------------------------------------------------
    ! Specify Dirichlet boundary conditions (prescribed uvel and vvel)
    !------------------------------------------------------------------------------

    ! initialize
    umask_dirichlet(:,:,:) = 0 
    vmask_dirichlet(:,:,:) = 0   

    if (whichbabc == HO_BABC_NO_SLIP .and. whichapprox /= HO_APPROX_DIVA) then
       ! Impose zero sliding everywhere at the bed
       ! Note: For the DIVA case, this BC is handled by setting beta_eff = 1/omega
       !TODO - Allow application of no-slip BC at selected basal nodes instead of all nodes?
       umask_dirichlet(nz,:,:) = 1    ! u = v = 0 at bed
       vmask_dirichlet(nz,:,:) = 1
    endif
       
    ! Set mask in columns identified in kinbcmask, typically read from file at initialization.
    ! Note: Assuming there is no vertical shear at these points, the bed velocity is the same
    !       as the velocity throughout the column.  This allows us to use the 3D umask_dirichlet
    !       and vmask_dirichlet with a 2D solver.
    ! TODO: Support Dirichlet condition with vertical shear for L1L2 and DIVA?
    !
    ! For a no-penetration global BC, set umask_dirichlet = 0 and uvel = 0.d0 along east/west global boundaries,
    !  and set vmask_dirichlet = 0 and vvel = 0.d0 along north/south global boundaries (based on umask_no_penetration
    !  and vmask_no_penetration, which are computed at initialization). 
    !
    ! For a 2D solve, initialize uvel_2d and vvel_2d at Dirichlet points to the bed velocity.

    do j = 1, ny-1
       do i = 1, nx-1

          ! if kinbcmask = 1, set Dirichlet masks for both uvel and vvel
          if (kinbcmask(i,j) == 1) then
             umask_dirichlet(:,i,j) = 1
             vmask_dirichlet(:,i,j) = 1
             if (solve_2d) then
                uvel_2d(i,j) = uvel(nz,i,j)
                vvel_2d(i,j) = vvel(nz,i,j)
             endif
          endif

          ! for the no-penetration global BC, prescribe zero outflow velocities
          ! (v = 0 at N/S boundaries, u = 0 at E/W boundaries)
          ! for other global BCs (periodic and outflow), umask_no_penetration = vmask_no_penetration = 0 everywhere

          if (umask_no_penetration(i,j) == 1) then
             umask_dirichlet(:,i,j) = 1
             uvel(:,i,j) = 0.d0
             if (solve_2d) uvel_2d(i,j) = 0.d0
          endif

          if (vmask_no_penetration(i,j) == 1) then
             vmask_dirichlet(:,i,j) = 1
             vvel(:,i,j) = 0.d0
             if (solve_2d) vvel_2d(i,j) = 0.d0
          endif

       enddo
    enddo

    !Note: The following halo updates are not needed here, provided that kinbcmask,
    !      umask_no_penetration and vmask_no_penetration receive halo updates
    !      (as done in glissade_initialise)
!    call staggered_parallel_halo(umask_dirichlet)
!    call staggered_parallel_halo(vmask_dirichlet)

    if (verbose_dirichlet .and. this_rank==rtest) then

       print*, ' '
       print*, 'kinbcmask:'
       write(6,'(a6)',advance='no')'        '
       do i = 1, xmax_print
          write(6,'(i6)',advance='no') i
       enddo
       write(6,*) ' '
       do j = ny-1, 1, -1
          write(6,'(i6)',advance='no') j
          do i = 1, xmax_print
             write(6,'(i6)',advance='no') kinbcmask(i,j)
          enddo
          write(6,*) ' '
       enddo

       print*, ' '
       print*, 'umask_no_penetration:'
       write(6,'(a6)',advance='no')'        '
       do i = 1, xmax_print
          write(6,'(i6)',advance='no') i
       enddo
       write(6,*) ' '
       do j = ny-1, 1, -1
          write(6,'(i6)',advance='no') j
          do i = 1, xmax_print
             write(6,'(i6)',advance='no') umask_no_penetration(i,j)
          enddo
          write(6,*) ' '
       enddo

       print*, ' '
       print*, 'vmask_no_penetration:'
       write(6,'(a6)',advance='no')'        '
       do i = 1, xmax_print
          write(6,'(i6)',advance='no') i
       enddo
       write(6,*) ' '
       do j = ny-1, 1, -1
          write(6,'(i6)',advance='no') j
          do i = 1, xmax_print
             write(6,'(i6)',advance='no') vmask_no_penetration(i,j)
          enddo
          write(6,*) ' '
       enddo

       print*, ' '
       print*, 'umask_dirichlet, k = 1:'
       write(6,'(a6)',advance='no') '        '
       do i = 1, xmax_print
          write(6,'(i6)',advance='no') i
       enddo
       write(6,*) ' '
       do j = ny-1, 1, -1
          write(6,'(i6)',advance='no') j
          do i = 1, xmax_print
             write(6,'(i6)',advance='no') umask_dirichlet(1,i,j)
          enddo
          write(6,*) ' '
       enddo

       print*, ' '
       print*, 'vmask_dirichlet, k = 1:'
       write(6,'(a6)',advance='no') '        '
       do i = 1, xmax_print
          write(6,'(i6)',advance='no') i
       enddo
       write(6,*) ' '
       do j = ny-1, 1, -1
          write(6,'(i6)',advance='no') j
          do i = 1, xmax_print
             write(6,'(i6)',advance='no') vmask_dirichlet(1,i,j)
          enddo
          write(6,*) ' '
       enddo

       print*, ' '
       print*, 'uvel, k = 1:'
       write(6,'(a10)',advance='no') '          '
!!       do i = 1, xmax_print
       do i = itest-3, itest+3
          write(6,'(i10)',advance='no') i
       enddo
       write(6,*) ' '
       do j = ny-1, 1, -1
          write(6,'(i10)',advance='no') j
!!       do i = 1, xmax_print
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') uvel(1,i,j)
          enddo
          write(6,*) ' '
       enddo

       print*, ' '
       print*, 'vvel, k = 1:'
       write(6,'(a10)',advance='no') '          '
!!       do i = 1, xmax_print
       do i = itest-3, itest+3
          write(6,'(i10)',advance='no') i
       enddo
       write(6,*) ' '
       do j = ny-1, 1, -1
          write(6,'(i10)',advance='no') j
!!       do i = 1, xmax_print
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') vvel(1,i,j)
          enddo
          write(6,*) ' '
       enddo

    endif   ! verbose_dirichlet

    !------------------------------------------------------------------------------
    ! Compute masks: 
    ! (1) ice mask = 1 in cells where ice is present (thck > thklim), = 0 elsewhere
    ! (2) floating mask = 1 in cells where ice is present and is floating
    ! (3) ocean mask = = 1 in cells where topography is below sea level and ice is absent
    ! (4) land mask = 1 in cells where topography is at or above sea level
    !------------------------------------------------------------------------------

    call glissade_get_masks(nx,          ny,         &
                            thck,        topg,       &
                            eus,         thklim,     &
                            ice_mask,    floating_mask, &
                            ocean_mask,  land_mask)

    !------------------------------------------------------------------------------
    ! Compute fraction of grounded ice in each cell
    ! (requires that thck and topg are up to date in halo cells).
    ! This is used below to compute the basal stress BC.
    !
    ! Three options for whichground:
    ! (0) HO_GROUND_NO_GLP: f_ground = 0 or 1 based on flotation criterion
    ! (1) HO_GROUND_GLP:    0 <= f_ground <= 1 based on grounding-line parameterization
    ! (2) HO_GROUND_ALL:    f_ground = 1 for all cells with ice
    !
    ! Three options for whichflotation_function (applies to whichground = 0 or 1):
    ! (0) HO_FLOTATION_FUNCTION_PATTYN:         f = (-rhow*b/rhoi*H) = f_pattyn; <=1 for grounded, > 1 for floating
    ! (1) HO_FLOTATION_FUNCTION_INVERSE_PATTYN: f = (rhoi*H)/(-rhow*b) = 1/f_pattyn; >=1 for grounded, < 1 for floating
    ! (2) HO_FLOTATION_FUNCTION_LINEAR:         f = -rhow*b - rhoi*H; <= 0 for grounded, > 0 for floating
    !
    ! f_flotation is not needed in further calculations but is output as a diagnostic.
    !------------------------------------------------------------------------------

    call glissade_grounded_fraction(nx,          ny,                      &
                                    thck,        topg,                    &
                                    eus,         ice_mask,                &
                                    whichground, whichflotation_function, &
                                    f_ground,    f_flotation)
    
    !------------------------------------------------------------------------------
    ! Compute ice thickness and upper surface on staggered grid
    ! (requires that thck and usrf are up to date in halo cells).
    ! For stagger_margin_in = 0, all cells (including ice-free) are included in interpolation.
    ! For stagger_margin_in = 1, only ice-covered cells are included.
    !------------------------------------------------------------------------------

!pw call t_startf('glissade_stagger')
    call glissade_stagger(nx,       ny,         &
                          thck,     stagthck,   &
                          ice_mask, stagger_margin_in = 1)

    call glissade_stagger(nx,       ny,         &
                          usrf,     stagusrf,   &
                          ice_mask, stagger_margin_in = 1)
!pw call t_stopf('glissade_stagger')

    !------------------------------------------------------------------------------
    ! Compute surface gradient on staggered grid
    ! (requires that usrf is up to date in halo cells)
    !
    ! Setting gradient_margin_in = 0 takes the gradient over all neighboring cells,
    !  including ice-free cells.  This is what Glide does, but is not appropriate
    !  if we have ice-covered floating cells next to ice-free ocean cells,
    !  because the gradient will be too big.
    ! Setting gradient_margin_in = 1 uses any available ice-covered cells
    !  and/or land cells to compute the gradient.  Requires passing in the surface
    !  elevation and a land mask.
    !  This is appropriate for both land-based problems and problems
    !  with ice shelves.  It is the default setting.
    ! Setting gradient_margin_in = 2 uses only ice-covered cells to compute
    !  the gradient.  This is appropriate for problems with ice shelves, but is
    !  is less accurate than options 0 or 1 for land-based problems (e.g., Halfar SIA).
    !
    ! Passing in max_slope ensures that the surface elevation gradient on the edge
    !  between two cells does not exceed a prescribed value.
    ! Although slope-limiting is not very physical, it helps prevent CFL violations
    !  in regions of steep coastal topography. Some input Greenland data sets have
    !  slopes of up to ~0.3 between adjacent grid cells, leading to very large velocities
    !  even with a no-slip basal boundary condition. 
    !
    ! Both the centered and upstream gradients are 2nd order accurate in space.
    ! The upstream gradient may be preferable for evolution problems using 
    !  whichapprox = HO_APPROX_BP or HO_APPROX_SIA, because in these cases
    !  the centered gradient fails to cancel checkerboard noise.
    ! The L1L2 solver computes 3D velocities in a way that damps checkerboard noise,
    !  so a centered difference should work well (and for the Halfar problem is more 
    !  accurate than upstream).
    !------------------------------------------------------------------------------

!pw call t_startf('glissade_gradient')

    if (whichgradient == HO_GRADIENT_CENTERED) then     ! 2nd order centered

       call glissade_centered_gradient(nx,               ny,         &
                                       dx,               dy,         &
                                       usrf,                         &
                                       dusrf_dx,         dusrf_dy,   &
                                       ice_mask,                     &
                                       gradient_margin_in = whichgradient_margin, &
                                       usrf = usrf,                  &
                                       land_mask = land_mask,        &
                                       max_slope = max_slope)

    else          ! 2nd order upstream

       call glissade_upstream_gradient(nx,             ny,           &
                                       dx,             dy,           &
                                       usrf,                         &
                                       dusrf_dx,       dusrf_dy,     &
                                       ice_mask,                     &
                                       accuracy_flag_in = 2,         &
                                       gradient_margin_in = whichgradient_margin, &
                                       usrf = usrf,                  &
                                       land_mask = land_mask,        &
                                       max_slope = max_slope)

    endif   ! whichgradient

!pw call t_stopf('glissade_gradient')

    if (verbose_glp .and. this_rank==rtest) then
       print*, 'effecpress_stag, rank =', rtest
       do j = jtest+1, jtest-1, -1
          write(6,'(a5)',advance='no') '    '
          do i = itest-3, itest+3
             write(6,'(f10.0)',advance='no') model%basal_physics%effecpress_stag(i,j)
          enddo
          print*, ' '
       enddo       
       print*, ' '
       print*, 'usrf, rank =', rtest
       do j = jtest+1, jtest-1, -1
          do i = itest-3, itest+3
             write(6,'(f10.2)',advance='no') usrf(i,j)
          enddo
          print*, ' '
       enddo       
       print*, ' '
       print*, 'thck, rank =', rtest
       do j = jtest+1, jtest-1, -1
          do i = itest-3, itest+3
             write(6,'(f10.2)',advance='no') thck(i,j)
          enddo
          print*, ' '
       enddo       
       print*, ' '
       print*, 'f_flotation, rank =', rtest
       do j = jtest+1, jtest-1, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') f_flotation(i,j)
          enddo
          print*, ' '
       enddo       
       print*, ' '
       print*, 'f_ground, rank =', rtest
       do j = jtest+1, jtest-1, -1
          write(6,'(a5)',advance='no') '    '
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') f_ground(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'dusrf_dx, rank =', rtest
       do j = jtest+1, jtest-1, -1
          write(6,'(a5)',advance='no') '    '
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') dusrf_dx(i,j)
          enddo
          print*, ' '
       enddo       
    endif

    if (verbose_gridop .and. this_rank==rtest) then
       print*, ' '
       print*, 'thck:'
       do j = ny, 1, -1
          do i = 1, nx
             write(6,'(f7.0)',advance='no') thck(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'stagthck, rank =',rtest
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(f7.0)',advance='no') stagthck(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'usrf:'
       do j = ny, 1, -1
          do i = 1, nx
             write(6,'(f7.0)',advance='no') usrf(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'dusrf_dx:'
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(f7.3)',advance='no') dusrf_dx(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'dusrf_dy:'
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(f7.3)',advance='no') dusrf_dy(i,j)
          enddo
          print*, ' '
       enddo

    endif  ! verbose_gridop

    !------------------------------------------------------------------------------
    ! Compute the vertices of each element.
    ! Identify the active cells (i.e., cells with thck > thklim,
    !  bordering a locally owned vertex) and active vertices (all vertices
    !  of active cells).
    ! Count the number of owned active nodes on this processor, and assign a 
    !  unique local ID to each such node.
    !TODO - Move Trilinos- and SLAP-specific computations to a different subroutine?
    !------------------------------------------------------------------------------

!pw call t_startf('glissade_get_vertex_geom')
    call get_vertex_geometry(nx,           ny,              &   
                             nz,           nhalo,           &
                             dx,           dy,              &
                             ice_mask,                      &
                             xVertex,      yVertex,         &
                             active_cell,  active_vertex,   &
                             nNodesSolve,  nVerticesSolve,  &
                             nodeID,       vertexID,        &
                             iNodeIndex,   jNodeIndex,  kNodeIndex, &
                             iVertexIndex, jVertexIndex)
!pw call t_stopf('glissade_get_vertex_geom')

    ! Zero out the velocity for inactive vertices
    do j = nhalo+1, ny-nhalo    ! locally owned vertices only
       do i = nhalo+1, nx-nhalo
          if (.not.active_vertex(i,j)) then
             uvel(:,i,j) = 0.d0
             vvel(:,i,j) = 0.d0
             if (solve_2d) then
                uvel_2d(i,j) = 0.d0
                vvel_2d(i,j) = 0.d0
             endif
          endif
       enddo
    enddo

    ! Assign the appropriate local ID to vertices and nodes in the halo.
    ! NOTE: This works for single-processor runs with periodic BCs
    !       (e.g., ISMIP-HOM), but not for multiple processors.

    call t_startf('glissade_halo_nodeID')
    call staggered_parallel_halo(nodeID)
    call staggered_parallel_halo(vertexID)
    call t_stopf('glissade_halo_nodeID')

    if (verbose_id .and. this_rank==rtest) then
       print*, ' '
       print*, 'vertexID before after halo update:'
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(i5)',advance='no') vertexID(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'nodeID after halo update, k = 1:'
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(i5)',advance='no') nodeID(1,i,j)
          enddo
          print*, ' '
       enddo
    endif

    ! Initialization for the Trilinos solver
    ! Allocate arrays, initialize the velocity solution, compute an array 
    !  that maps the local index for owned active nodes to a unique global ID,
    !  and communicate this array to Trilinos

#ifdef TRILINOS
    if (whichsparse == HO_SPARSE_TRILINOS) then   

       if (solve_2d) then

          allocate(active_owned_unknown_map(2*nVerticesSolve))
          allocate(velocityResult(2*nVerticesSolve))
          allocate(Afill_2d(nNodeNeighbors_2d,nx-1,ny-1))

          !----------------------------------------------------------------
          ! Compute global IDs needed to initialize the Trilinos solver
          !----------------------------------------------------------------

          call t_startf('glissade_trilinos_glbid')
          call trilinos_global_id_2d(nx,             ny,           &
                                     nVerticesSolve,               &
                                     iVertexIndex,   jVertexIndex, &
                                     global_vertex_id,             &
                                     active_owned_unknown_map)
          call t_stopf('glissade_trilinos_glbid')

          !----------------------------------------------------------------
          ! Send this information to Trilinos (trilinosGlissadeSolver.cpp)
          !----------------------------------------------------------------

          call t_startf('glissade_init_tgs')
          call initializetgs(2*nVerticesSolve, active_owned_unknown_map, comm)
          call t_stopf('glissade_init_tgs')

          !----------------------------------------------------------------
          ! If this is the first outer iteration, then save the pattern of matrix
          ! values that are potentially nonzero and should be sent to Trilinos.
          ! Trilinos requires that this pattern remains fixed during the outer loop.
          !----------------------------------------------------------------

          call t_startf('glissade_trilinos_fill_pattern')
          call trilinos_fill_pattern_2d(nx,            ny,              &
                                        active_vertex, nVerticesSolve,  &
                                        iVertexIndex,  jVertexIndex,    &
                                        indxA_2d,      Afill_2d)
          call t_stopf('glissade_trilinos_fill_pattern')

          !----------------------------------------------------------------
          ! Initialize the solution vector from uvel/vvel.
          !----------------------------------------------------------------

          call trilinos_init_velocity_2d(nx,           ny,           &
                                         nVerticesSolve,             &
                                         iNodeIndex,   jNodeIndex,   &
                                         uvel_2d,      vvel_2d,      &
                                         velocityResult)

       else   ! 3D solve

          allocate(active_owned_unknown_map(2*nNodesSolve))
          allocate(velocityResult(2*nNodesSolve))
          allocate(Afill(nNodeNeighbors_3d,nz,nx-1,ny-1))

          !----------------------------------------------------------------
          ! Compute global IDs needed to initialize the Trilinos solver
          !----------------------------------------------------------------

          call t_startf('glissade_trilinos_glbid')
          call trilinos_global_id_3d(nx,         ny,         nz,   &
                                     nNodesSolve,                  &
                                     iNodeIndex, jNodeIndex, kNodeIndex,  &
                                     global_node_id,               &
                                     active_owned_unknown_map)
          call t_stopf('glissade_trilinos_glbid')

          !----------------------------------------------------------------
          ! Send this information to Trilinos (trilinosGlissadeSolver.cpp)
          !----------------------------------------------------------------

          call t_startf('glissade_init_tgs')
          call initializetgs(2*nNodesSolve, active_owned_unknown_map, comm)
          call t_stopf('glissade_init_tgs')

          !----------------------------------------------------------------
          ! If this is the first outer iteration, then save the pattern of matrix
          ! values that are potentially nonzero and should be sent to Trilinos.
          ! Trilinos requires that this pattern remains fixed during the outer loop.
          !----------------------------------------------------------------

          call t_startf('glissade_trilinos_fill_pattern')
          call trilinos_fill_pattern_3d(nx,            ny,           nz,   &
                                        active_vertex, nNodesSolve,        &
                                        iNodeIndex,    jNodeIndex,   kNodeIndex,  &
                                        indxA_3d,      Afill)
                                     
          call t_stopf('glissade_trilinos_fill_pattern')

          !----------------------------------------------------------------
          ! Initialize the solution vector from uvel/vvel.
          !----------------------------------------------------------------

          call trilinos_init_velocity_3d(nx,           ny,                       &
                                         nz,           nNodesSolve,              &
                                         iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                                         uvel,         vvel,                     &
                                         velocityResult)

       endif   ! whichapprox
    endif      ! whichsparse
#endif

    !------------------------------------------------------------------------------
    ! Initialize the basal traction parameter, beta_internal.
    ! Note: If beta is read from an external file, the external value should not be changed.
    !        This value is saved in model%velocity%beta.
    !       The glissade solver uses a beta field weighted by f_ground.
    !        This field is stored in model%velocity%beta_internal and can change over time.
    !       For a no-slip boundary condition (HO_BABC_NO_SLIP), beta_internal is not computed,
    !        so beta_internal = 0 will be written to output.
    !------------------------------------------------------------------------------

    beta_internal(:,:) = 0.d0

    !------------------------------------------------------------------------------
    ! For the HO_BABC_BETA_TPMP option (lower beta where the bed is at the pressure melting point),
    ! compute the bed temperature and the bed pressure-melting-point temperature at vertices.
    !------------------------------------------------------------------------------
  
    if (whichbabc == HO_BABC_BETA_TPMP) then

       ! interpolate bed temperature to vertices
       ! For stagger_margin_in = 1, only ice-covered cells are included in the interpolation

       call glissade_stagger(nx,           ny,           &
                             temp(nz,:,:), stagbedtemp,  &
                             ice_mask,     stagger_margin_in = 1)
       
       ! compute pressure melting point temperature at bed
       do j = 1, ny
          do i = 1, nx
             if (ice_mask(i,j) == 1) then
                call glissade_pressure_melting_point(thck(i,j), bedpmp(i,j))
             else
                bedpmp(i,j) = 0.d0
             endif
          enddo
       enddo

       ! interpolate bed pmp temperature to vertices
       call glissade_stagger(nx,           ny,           &
                             bedpmp(:,:),  stagbedpmp(:,:), &
                             ice_mask,     stagger_margin_in = 1)

       ! compute a pmp mask at vertices; this mask is passed to calcbeta below
       where (abs(stagbedpmp - stagbedtemp) < 1.d-3 .and. active_vertex)
          pmp_mask = 1
       elsewhere
          pmp_mask = 0
       endwhere

    else    ! not HO_BABC_BETA_TPMP

       pmp_mask(:,:) = 0

    endif

    !------------------------------------------------------------------------------
    ! Compute the factor A^(-1/n) appearing in the expression for effective viscosity.
    ! This factor is often denoted as B in the literature.
    ! Note: The rate factor (flwa = A) is assumed to have units of Pa^(-n) yr^(-1).
    !       Thus flwafact = 0.5 * A^(-1/n) has units Pa yr^(1/n).
    !------------------------------------------------------------------------------

    flwafact(:,:,:) = 0.d0

    ! Loop over all cells that border locally owned vertices
    !TODO - Simply compute flwafact for all cells?  We should have flwa for all cells.

    do j = 1+nhalo, ny-nhalo+1
       do i = 1+nhalo, nx-nhalo+1
          if (active_cell(i,j)) then
             ! gn = exponent in Glen's flow law (= 3 by default)
             flwafact(:,i,j) = 0.5d0 * flwa(:,i,j)**(-1.d0/real(gn,dp))  
          endif
       enddo
    enddo

    !------------------------------------------------------------------------------
    ! If using SLAP solver, then allocate space for the sparse matrix (A), rhs (b), 
    !  answer (x), and residual vector (Ax-b).
    !------------------------------------------------------------------------------

    if (whichsparse <= HO_SPARSE_GMRES) then  ! using SLAP solver

       if (solve_2d) then
          matrix_order = 2*nVerticesSolve
          max_nonzeros = matrix_order*2*nNodeNeighbors_2d  ! nNodeNeighbors_2d = 9
                                                           ! 18 = 2 * 9 (since solving for both u and v)
       else  ! 3D solve
          matrix_order = 2*nNodesSolve
          max_nonzeros = matrix_order*2*nNodeNeighbors_3d  ! nNodeNeighbors_3d = 27
                                                           ! 54 = 2 * 27 (since solving for both u and v)
       endif

       allocate(matrix%row(max_nonzeros), matrix%col(max_nonzeros), matrix%val(max_nonzeros))
       allocate(rhs(matrix_order), answer(matrix_order), resid_vec(matrix_order))

       answer(:) = 0.d0
       rhs(:) = 0.d0
       resid_vec(:) = 0.d0

       if (verbose_matrix) then
          print*, 'matrix_order =', matrix_order
          print*, 'max_nonzeros = ', max_nonzeros
       endif

    endif   ! SLAP solver
 
    !---------------------------------------------------------------
    ! Print some diagnostic info
    !---------------------------------------------------------------

    if (main_task) then
       print *, ' '
       print *, 'Running Glissade higher-order dynamics solver'
       print *, ' '
       if (whichresid == HO_RESID_L2NORM) then  ! use L2 norm of residual
          print *, 'iter #     resid (L2 norm)       target resid'
       elseif (whichresid == HO_RESID_L2NORM_RELATIVE) then  ! relative residual, |Ax-b|/|b|
          print *, 'iter #     resid, |Ax-b|/|b|     target resid'
       else                                     ! residual based on velocity
          print *, 'iter #     velo resid            target resid'
       end if
       print *, ' '
    endif

    !------------------------------------------------------------------------------
    ! Set initial solver values 
    !------------------------------------------------------------------------------

    counter = 0
    resid_velo = 1.d0

    L2_norm   = 1.0d20      ! arbitrary large value
    L2_target = 1.0d-4

    !WHL: For standard test cases (dome, circular shelf), a relative target of 1.0d-7 is 
    !     roughly as stringent as an absolute target of 1.0d-4.
    !       
    L2_norm_relative = 1.0d20
    L2_target_relative = 1.0d-7

    outer_it_criterion = 1.0d10   ! guarantees at least one loop
    outer_it_target    = 1.0d-12

    !------------------------------------------------------------------------------
    ! Assemble the load vector b
    ! This goes before the outer loop because the load vector
    !  does not change from one nonlinear iteration to the next.
    !------------------------------------------------------------------------------

    loadu(:,:,:) = 0.d0
    loadv(:,:,:) = 0.d0

    !------------------------------------------------------------------------------
    ! Gravitational forcing
    !------------------------------------------------------------------------------

    call t_startf('glissade_load_vector_gravity')

    call load_vector_gravity(nx,               ny,              &
                             nz,               sigma,           &
                             nhalo,            active_cell,     &
                             xVertex,          yVertex,         &
                             stagusrf,         stagthck,        &
                             dusrf_dx,         dusrf_dy,        &
                             whichassemble_taud,                &
                             loadu,            loadv)
       
    call t_stopf('glissade_load_vector_gravity')

    ! Compute components of gravitational driving stress
    taudx(:,:) = 0.d0
    taudy(:,:) = 0.d0
    do j = 1, ny-1
       do i = 1, nx-1
          do k = 1, nz
             taudx(i,j) = taudx(i,j) + loadu(k,i,j)
             taudy(i,j) = taudy(i,j) + loadv(k,i,j)
          enddo
       enddo
    enddo
    taudx(:,:) = taudx(:,:) * vol0/(dx*dy)  ! convert from model units to Pa
    taudy(:,:) = taudy(:,:) * vol0/(dx*dy)

    if (verbose_glp .and. this_rank==rtest) then
       ! Note: The first of these quantities is the load vector on the rhs of the matrix
       !       The second is the value that would go on the rhs by simply taking rho*g*H*ds/dx.
       !       These will not agree exactly because of the way H is handled in FE assembly,
       !        but they should be close if which_ho_assemble_taud = HO_ASSEMBLE_TAUD_LOCAL.
       !       If which_ho_assemble_taud = HO_ASSEMBLE_TAUD_STANDARD, they can differ substantially.

       print*, ' '
       print*, 'vert sum of grav load vector, rank =', rtest
       do j = jtest+1, jtest-1, -1
          write(6,'(a5)',advance='no') '    '
          do i = itest-3, itest+3
             write(6,'(f10.0)',advance='no') taudx(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'rho*g*H*ds/dx, rank =', rtest
       do j = jtest+1, jtest-1, -1
          write(6,'(a5)',advance='no') '    '
          do i = itest-3, itest+3
             write(6,'(f10.0)',advance='no') -rhoi*grav*stagthck(i,j)*dusrf_dx(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'Starting uvel_2d, rank =', rtest
       do j = jtest+1, jtest-1, -1
          write(6,'(a5)',advance='no') '    '
          do i = itest-3, itest+3
             write(6,'(f10.2)',advance='no') uvel_2d(i,j)
          enddo
          print*, ' '
       enddo       
    endif

    !------------------------------------------------------------------------------
    ! Lateral pressure at vertical ice edge
    !------------------------------------------------------------------------------

    call t_startf('glissade_load_vector_lateral_bc')
    call load_vector_lateral_bc(nx,               ny,              &
                                nz,               sigma,           &
                                nhalo,                             &
                                floating_mask,    ocean_mask,      &
                                active_cell,                       &
                                xVertex,          yVertex,         &
                                stagusrf,         stagthck,        &
                                loadu,            loadv)
    call t_stopf('glissade_load_vector_lateral_bc')

    call t_stopf('glissade_vhs_init')

    !------------------------------------------------------------------------------
    ! If solving a 2D problem (e.g., SSA at one level), sum the load vector over columns.
    ! Note: It would be slightly more efficient to compute the load vector at a single level
    !       using custom 2D subroutines. However, this would require extra code and would
    !       save little work, since the load vector is computed only once per timestep.
    !------------------------------------------------------------------------------

    if (solve_2d) then

       loadu_2d(:,:) = 0.d0
       loadv_2d(:,:) = 0.d0

       do j = 1, ny-1
          do i = 1, nx-1
             do k = 1, nz
                loadu_2d(i,j) = loadu_2d(i,j) + loadu(k,i,j)
                loadv_2d(i,j) = loadv_2d(i,j) + loadv(k,i,j)
             enddo
          enddo
       enddo

    endif

    !------------------------------------------------------------------------------
    ! Main outer loop: Iterate to solve the nonlinear problem
    !------------------------------------------------------------------------------

    call t_startf('glissade_vhs_nonlinear_loop')
    do while (outer_it_criterion >= outer_it_target .and. counter < maxiter_nonlinear)

       ! Advance the iteration counter

       counter = counter + 1

       !---------------------------------------------------------------------------
       ! Compute or prescribe the basal traction field 'beta'.
       !
       ! Notes:
       ! (1) We could compute beta before the main outer loop if beta
       !     were assumed to be independent of velocity.  Computing beta here,
       !     however, allows for more general sliding laws.
       ! (2) The units of the input arguments in calcbeta are assumed to be the
       !     same as the Glissade units.
       ! (3) The computed beta (called beta_internal) is weighted by f_ground, 
       !     the grounded fraction at each vertex.  With a GLP, f_ground is 
       !     between 0 and 1 for vertices adjacent to the GL, allowing for a smooth 
       !     change in beta as the GL advances and retreats.
       ! (4) The basal velocity is a required input to calcbeta.  
       !     DIVA does not compute the basal velocity in the 2D matrix solve, 
       !     but computes the 3D velocity after each iteration so that
       !     uvel/vvel(nz,:,:) are available here.
       ! (5) For which_ho_babc = HO_BABC_EXTERNAL_BETA, beta is passed in
       !     with dimensionless Glimmer units. Rather than incur roundoff errors by
       !     repeatedly multiplying and dividing by scaling constants, the conversion
       !     to Pa yr/m is done here in the argument list.
       !-------------------------------------------------------------------

       if (whichapprox == HO_APPROX_SSA .or. whichapprox == HO_APPROX_L1L2) then
          ubas(:,:) = uvel_2d(:,:)
          vbas(:,:) = vvel_2d(:,:)
       else  ! 3D solve or DIVA
          ubas(:,:) = uvel(nz,:,:)
          vbas(:,:) = vvel(nz,:,:)
       endif

!!       if (verbose_beta .and. this_rank==rtest) then
       if (verbose_beta .and. this_rank==rtest .and. mod(counter,5)==0) then

          print*, ' '
          print*, 'Before calcbeta, counter =', counter
          print*, ' '
          print*, 'usrf field, itest, jtest, rank =', itest, jtest, rtest
!!          do j = ny-1, 1, -1
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
!!             do i = 1, nx-1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') usrf(i,j)
             enddo
             write(6,*) ' '
          enddo          

          print*, ' '
          print*, 'thck field, itest, jtest, rank =', itest, jtest, rtest
!!          do j = ny-1, 1, -1
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
!!             do i = 1, nx-1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'topg field, itest, jtest, rank =', itest, jtest, rtest
!!          do j = ny-1, 1, -1
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
!!             do i = 1, nx-1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') topg(i,j)
             enddo
             write(6,*) ' '
          enddo          

          print*, ' '
          print*, 'f_flotation, itest, jtest, rank =', itest, jtest, rtest
!!          do j = ny-1, 1, -1
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
!!             do i = 1, nx-1
             do i = itest-3, itest+3
                write(6,'(e10.3)',advance='no') f_flotation(i,j)
             enddo
             write(6,*) ' '
          enddo          

          print*, ' '
          print*, 'f_ground field, itest, jtest, rank =', itest, jtest, rtest
!!          do j = ny-1, 1, -1
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
!!             do i = 1, nx-1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') f_ground(i,j)
             enddo
             write(6,*) ' '
          enddo          

       endif  ! verbose_beta

       call calcbeta (whichbabc,                        &
                      dx,            dy,                &
                      nx,            ny,                &
                      ubas,          vbas,              &
                      bwat,          ho_beta_const,     &
                      mintauf,                          &
                      model%basal_physics,              &
                      flwa(nz-1,:,:),                   &  ! basal flwa layer
                      thck,                             &
                      stagmask,                         &
                      beta*tau0/(vel0*scyr),            &  ! external beta (intent in)
                      beta_internal,                    &  ! beta weighted by f_ground (intent out)
                      f_ground,                         &
                      pmp_mask,                         &  ! needed for HO_BABC_BETA_TPMP option
                      beta_grounded_min)

       call staggered_parallel_halo(beta_internal)

       if (verbose_beta) then
          maxbeta = maxval(beta_internal(:,:))
          maxbeta = parallel_reduce_max(maxbeta)
          minbeta = minval(beta_internal(:,:))
          minbeta = parallel_reduce_min(minbeta)
       endif

       if (verbose_beta .and. main_task) then
          print*, ' '
          print*, 'max, min beta (Pa/(m/yr)) =', maxbeta, minbeta
       endif

!!       if (verbose_beta .and. this_rank==rtest) then
       if (verbose_beta .and. this_rank==rtest .and. mod(counter,5)==0) then

          print*, ' '
          print*, 'Weighted beta field, itest, jtest, rank =', itest, jtest, rtest
!!          do j = ny-1, 1, -1
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
!!             do i = 1, nx-1
             do i = itest-3, itest+3
                write(6,'(e10.3)',advance='no') beta_internal(i,j)
             enddo
             write(6,*) ' '
          enddo          

          print*, ' '
          print*, 'Basal uvel field, itest, jtest, rank =', itest, jtest, rtest
!!          do j = ny-1, 1, -1
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
!!             do i = 1, nx-1
             do i = itest-3, itest+3
                write(6,'(e10.3)',advance='no') uvel(nz,i,j)
             enddo
             write(6,*) ' '
          enddo          

          print*, ' '
          print*, 'Basal vvel field, itest, jtest, rank =', itest, jtest, rtest
!!          do j = ny-1, 1, -1
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
!!             do i = 1, nx-1
             do i = itest-3, itest+3
                write(6,'(e10.3)',advance='no') vvel(nz,i,j)
             enddo
             write(6,*) ' '
          enddo          

          !TODO - Remove the remaining verbose_beta diagnostics?
          !       They are not specific to beta but are useful for diagnosing CFL issues.
          print*, ' '
          print*, 'Sfc uvel field, itest, jtest, rank =', itest, jtest, rtest
!!          do j = ny-1, 1, -1
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
!!             do i = 1, nx-1
             do i = itest-3, itest+3
                write(6,'(e10.3)',advance='no') uvel(1,i,j)
             enddo
             write(6,*) ' '
          enddo          

          print*, ' '
          print*, 'Sfc vvel field, itest, jtest, rank =', itest, jtest, rtest
!!          do j = ny-1, 1, -1
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
!!             do i = 1, nx-1
             do i = itest-3, itest+3
                write(6,'(e10.3)',advance='no') vvel(1,i,j)
             enddo
             write(6,*) ' '
          enddo          

          print*, ' '
          print*, 'dusrf/dx field, itest, jtest, rank =', itest, jtest, rtest
!!          do j = ny-1, 1, -1
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
!!             do i = 1, nx-1
             do i = itest-3, itest+3
                write(6,'(f10.5)',advance='no') dusrf_dx(i,j)
             enddo
             write(6,*) ' '
          enddo          

          print*, ' '
          print*, 'dusrf_dy field, itest, jtest, rank =', itest, jtest, rtest
!!          do j = ny-1, 1, -1
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
!!             do i = 1, nx-1
             do i = itest-3, itest+3
                write(6,'(f10.5)',advance='no') dusrf_dy(i,j)
             enddo
             write(6,*) ' '
          enddo          

          if (whichbabc == HO_BABC_BETA_TPMP) then

             print*, ' '
             print*, 'stagbedtemp, itest, jtest, rank =', itest, jtest, rtest
!!             do j = ny-1, 1, -1
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
!!                do i = 1, nx-1
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') stagbedtemp(i,j)
                enddo
                write(6,*) ' '
             enddo             

             print*, ' '
             print*, 'stagbedpmp, itest, jtest, rank =', itest, jtest, rtest
!!             do j = ny-1, 1, -1
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
!!                do i = 1, nx-1
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') stagbedpmp(i,j)
                enddo
                write(6,*) ' '
             enddo             

             print*, ' '
             print*, 'pmp_mask, itest, jtest, rank =', itest, jtest, rtest
!!             do j = ny-1, 1, -1
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
!!                do i = 1, nx-1
                do i = itest-3, itest+3
                   write(6,'(i10)',advance='no') pmp_mask(i,j)
                enddo
                write(6,*) ' '
             enddo             

          endif  ! HO_BABC_BETA_TPMP

          if (whichbabc == HO_BABC_YIELD_PICARD) then
             print*, ' '
             print*, 'mintauf field, rank =', rtest
             do j = ny-1, 1, -1
                write(6,'(i6)',advance='no') j
                do i = 1, nx-1
                   write(6,'(e10.3)',advance='no') mintauf(i,j)
                enddo
                write(6,*) ' '
             enddo
          endif

          if (whichbabc == HO_BABC_COULOMB_FRICTION .or. &
              whichbabc == HO_BABC_COULOMB_CONST_BASAL_FLWA) then
             print*, ' '
             print*, 'C_space_factor_stag, itest, rank =', itest, rtest
             do j = ny-1, 1, -1
                write(6,'(i6)',advance='no') j
!!                do i = 1, nx-1
                do i = itest-4, itest+4
                   write(6,'(f10.3)',advance='no') model%basal_physics%C_space_factor_stag(i,j)
                enddo
                write(6,*) ' '
             enddo
          endif

       endif   ! verbose_beta

       !-------------------------------------------------------------------
       ! Assemble the linear system Ax = b
       !
       ! Depending on the value of whichapprox, we can assemble either a 2D system
       ! (to solve for uvel and vvel at one level) or a 3D system (to solve for
       !  uvel and vvel at all levels).
       !-------------------------------------------------------------------
       
       if (solve_2d) then  ! assemble 2D matrix

          call t_startf('glissade_assemble_2d')

          ! save current velocity
          usav_2d(:,:) = uvel_2d(:,:)
          vsav_2d(:,:) = vvel_2d(:,:)

          if (whichapprox==HO_APPROX_DIVA .and. verbose_diva .and. this_rank==rtest) then
             i = itest
             j = jtest
             print*, ' '
             print*, 'i, j, uvel_2d, vvel_2d, beta_eff, btractx, btracty:',  &
                      i, j, uvel_2d(i,j), vvel_2d(i,j), beta_eff(i,j), btractx(i,j), btracty(i,j)
          endif

          ! Assemble the matrix
          !TODO - Different calls for SSA, L1L2 and DIVA?
          
          call assemble_stiffness_matrix_2d(nx,               ny,              &
                                            nz,                                &
                                            sigma,            stagsigma,       &
                                            nhalo,            active_cell,     &
                                            xVertex,          yVertex,         &
                                            uvel_2d,          vvel_2d,         &
                                            stagusrf,         stagthck,        &
                                            flwa,             flwafact,        &
                                            whichapprox,                       &
                                            whichefvs,        efvs_constant,   &
                                            efvs,                              &
                                            Auu_2d,           Auv_2d,          &
                                            Avu_2d,           Avv_2d,          &
                                            dusrf_dx,         dusrf_dy,        &
                                            thck,                              &          
                                            btractx,          btracty,         &
                                            omega_k,          omega,   &
                                            efvs_qp_3d)

          if (whichapprox == HO_APPROX_DIVA) then

             ! Halo update for omega
             ! This is needed so that beta_eff, computed below, will be correct in halos

             call parallel_halo(omega)

             ! Interpolate the appropriate integral
             if (diva_level_index == 0) then   ! solving for 2D mean velocity field

                ! Interpolate omega to the staggered grid
                call glissade_stagger(nx,           ny,                       &
                                      omega(:,:),   stag_omega(:,:),  &
                                      ice_mask,     stagger_margin_in = 1)

             else  ! solving for the velocity at level k (k = 1 at upper surface)

                k = diva_level_index

                call parallel_halo(omega_k(k,:,:))

                call glissade_stagger(nx,              ny,                           &
                                      omega_k(k,:,:),  stag_omega(:,:),  &
                                      ice_mask,        stagger_margin_in = 1)
                
             endif
                
             !-------------------------------------------------------------------
             ! Compute effective beta based on Goldberg (2011) eq. 40 and 41
             !
             ! If solving for the depth-integrated velocity u_mean:
             !
             !       beta_eff * u_mean = beta * u_b
             !
             ! where beta_eff = beta / (1 + beta*omega)
             !          omega = int_b^z {[(s-z)/H]^2 * 1/efvs * dz}
             !   
             ! If solving for the surface velocity u_sfc:
             !
             !       beta_eff * u_sfc = beta * u_b
             !
             ! where beta_eff = beta / (1 + beta*omega_1)
             !        omega_1 = int_b^s {[(s-z)/H] * 1/efvs * dz}
             !                = omega_k for k = 1
             !
             ! To implement a no-slip basal BC, set beta_eff = 1/omega
             !--------------------------------------------------------------------
    
             beta_eff(:,:) = 0.d0

             if (whichbabc == HO_BABC_NO_SLIP) then
                where (stag_omega > 0.d0) beta_eff = 1.d0 / stag_omega
             else   ! slip allowed at bed
                beta_eff(:,:) = beta_internal(:,:) / (1.d0 + beta_internal(:,:)*stag_omega(:,:))
             endif

             if (verbose_diva .and. this_rank==rtest) then
                i = itest
                j = jtest
                print*, ' '
                print*, 'uvel, F2, beta_eff, btractx:', uvel_2d(i,j), stag_omega(i,j), beta_eff(i,j), btractx(i,j)
                print*, 'vvel, btracty:', vvel_2d(i,j), btracty(i,j)
                print*, ' '
                print*, 'beta_eff:'
!!                do j = ny-1, 1, -1
!!                   do i = 1, nx-1
                do j = jtest-5, jtest+5, -1
                   do i = itest-5, itest+5
                      write(6,'(e10.3)',advance='no') beta_eff(i,j)
                   enddo
                   write(6,*) ' '
                enddo
             endif

             ! Incorporate basal sliding boundary conditions, based on beta_eff

             call basal_sliding_bc(nx,                ny,              &
                                   nNodeNeighbors_2d, nhalo,           &
                                   active_cell,       beta_eff,        &
                                   xVertex,           yVertex,         &
                                   whichassemble_beta,                 &
                                   Auu_2d,            Avv_2d)

          else    ! L1L2, SSA

             ! Incorporate basal sliding boundary conditions, based on beta_internal

             call basal_sliding_bc(nx,                ny,              &
                                   nNodeNeighbors_2d, nhalo,           &
                                   active_cell,       beta_internal,   &
                                   xVertex,           yVertex,         &
                                   whichassemble_beta,                 &
                                   Auu_2d,            Avv_2d)

          endif    ! whichapprox (SSA, L1L2, DIVA)

          call t_stopf('glissade_assemble_2d')

          if (verbose_matrix .and. this_rank==rtest) print*, 'Assembled the 2D stiffness matrix'

          !---------------------------------------------------------------------------
          ! Set rhs to the load vector
          ! The rhs can be adjusted below to account for inhomogeneous Dirichlet BC
          !---------------------------------------------------------------------------

          bu_2d(:,:) = loadu_2d(:,:)
          bv_2d(:,:) = loadv_2d(:,:)

          !---------------------------------------------------------------------------
          ! Incorporate Dirichlet boundary conditions (prescribed uvel and vvel)
          ! Note: With a no-slip BC, umask_dirichlet(nz,:,:) = vmask_dirichlet(nz,:,:) = .true., 
          !        except for the DIVA scheme.
          !       For DIVA, the no-slip BC is enforced by setting beta_eff = 1/omega.
          !---------------------------------------------------------------------------

          if (verbose_dirichlet .and. main_task) then
             print*, 'Call Dirichlet_bc'
          endif

          call t_startf('glissade_dirichlet_2d')
          call dirichlet_boundary_conditions_2d(nx,                       ny,                      &
                                                nhalo,                                             &
                                                active_vertex,                                     &
                                                umask_dirichlet(nz,:,:),  vmask_dirichlet(nz,:,:), &
                                                uvel_2d,                  vvel_2d,                 &
                                                Auu_2d,                   Auv_2d,                  &
                                                Avu_2d,                   Avv_2d,                  &
                                                bu_2d,                    bv_2d)
          call t_stopf('glissade_dirichlet_2d')

          !---------------------------------------------------------------------------
          ! Halo updates for matrices
          !
          ! These updates are not strictly necessary unless we're concerned about
          !  roundoff errors.  See comments below under 3D assembly.
          !---------------------------------------------------------------------------
     
          call t_startf('glissade_halo_Axxs')
          call staggered_parallel_halo(Auu_2d(:,:,:))
          call staggered_parallel_halo(Auv_2d(:,:,:))
          call staggered_parallel_halo(Avu_2d(:,:,:))
          call staggered_parallel_halo(Avv_2d(:,:,:))
          call t_stopf('glissade_halo_Axxs')
          
          !---------------------------------------------------------------------------
          ! Halo updates for rhs vectors
          ! (Not sure if these are necessary, but leaving them for now)
          !---------------------------------------------------------------------------

          call t_startf('glissade_halo_bxxs')
          call staggered_parallel_halo(bu_2d(:,:))
          call staggered_parallel_halo(bv_2d(:,:))
          call t_stopf('glissade_halo_bxxs')

          !---------------------------------------------------------------------------
          ! Check symmetry of assembled matrix
          ! 
          ! There may be small differences from perfect symmetry due to roundoff errors.  
          ! If sufficiently small, these differences are fixed by averaging the two values 
          !  that should be symmetric.  Otherwise the code aborts.
          !---------------------------------------------------------------------------

          if (check_symmetry) then

             call t_startf('glissade_chk_symmetry')
             call check_symmetry_assembled_matrix_2d(nx,          ny,      &
                                                     nhalo,                &
                                                     active_vertex,        &
                                                     Auu_2d,      Auv_2d,  &
                                                     Avu_2d,      Avv_2d)
             call t_stopf('glissade_chk_symmetry')

          endif

          !---------------------------------------------------------------------------
          ! Count the total number of nonzero entries on all processors.
          !---------------------------------------------------------------------------

          call count_nonzeros_2d(nx,      ny,     &
                                 nhalo,           &
                                 Auu_2d,  Auv_2d, &
                                 Avu_2d,  Avv_2d, &
                                 active_vertex,   &
                                 nNonzeros)

          if (write_matrix) then
             if (counter == 1) then    ! first outer iteration only
 
                call t_startf('glissade_wrt_mat')
                call write_matrix_elements_2d(nx,             ny,            &
                                              nVerticesSolve, vertexID,      &
                                              iVertexIndex,   jVertexIndex,  &
                                              Auu_2d,         Auv_2d,        &
                                              Avu_2d,         Avv_2d,        &
                                              bu_2d,          bv_2d)
                call t_stopf('glissade_wrt_mat')

             endif
          endif   ! write_matrix

          if (verbose_matrix .and. this_rank==rtest) then
             i = itest
             j = jtest
             print*, ' '
             print*, 'After assembly and BC, i, j =', i, j
             print*, 'Auu_2d sum =', sum(Auu_2d(:,i,j))
             print*, 'Auv_2d sum =', sum(Auv_2d(:,i,j))
             print*, 'Avu_2d sum =', sum(Avu_2d(:,i,j))
             print*, 'Avv_2d sum =', sum(Avv_2d(:,i,j))

             m = indxA_2d(0,0)  ! diag entry
             print*, ' '
             print*, 'Matrix row properties, j =', j
             print*, ' '
             print*, 'i, diag, max, min, sum:'
!!             do i = 1, nx-1
             do i = itest-3, itest+3
                print*, ' '
                write(6,'(a8, i4, 4f20.8)') 'Auu_2d:', i, Auu_2d(m,i,j), maxval(Auu_2d(:,i,j)), &
                                                   minval(Auu_2d(:,i,j)),   sum(Auu_2d(:,i,j))
                write(6,'(a8, i4, 4f20.8)') 'Auv_2d:', i, Auv_2d(m,i,j), maxval(Auv_2d(:,i,j)), &
                                                   minval(Auv_2d(:,i,j)),   sum(Auv_2d(:,i,j))
             enddo

             i = itest
             j = jtest
             print*, 'i, j =', i, j
             print*, 'iA, jA, Auu_2d, Auv_2d, Avu_2d, Avv_2d:'
             do jA = -1, 1
                do iA = -1, 1
                   m = indxA_2d(iA,jA)
                   print*, iA, jA, Auu_2d(m,i,j), Auv_2d(m,i,j), Avu_2d(m,i,j), Avv_2d(m,i,j) 
                enddo
             enddo
             print*, ' '
             print*, 'bu_2d =', bu_2d(i,j)
             print*, 'bv_2d =', bv_2d(i,j)
             
          endif  ! verbose_matrix

       else  ! assemble 3D matrix

          ! save current velocity
          usav(:,:,:) = uvel(:,:,:)
          vsav(:,:,:) = vvel(:,:,:)

          !---------------------------------------------------------------------------
          ! Assemble the stiffness matrix A
          !---------------------------------------------------------------------------

          call t_startf('glissade_assemble_3d')
          call assemble_stiffness_matrix_3d(nx,               ny,              &
                                            nz,               sigma,           &
                                            nhalo,            active_cell,     &
                                            xVertex,          yVertex,         &
                                            uvel,             vvel,            &
                                            stagusrf,         stagthck,        &
                                            flwafact,         whichapprox,     &
                                            efvs,             whichefvs,       &
                                            efvs_constant,                     &
                                            Auu,              Auv,             &
                                            Avu,              Avv)
          call t_stopf('glissade_assemble_3d')
          
          if (verbose_matrix .and. this_rank==rtest) print*, 'Assembled the 3D stiffness matrix'

          !---------------------------------------------------------------------------
          ! Incorporate basal sliding boundary conditions, based on beta_internal
          !---------------------------------------------------------------------------

          if (whichbabc /= HO_BABC_NO_SLIP) then

             call basal_sliding_bc(nx,                  ny,              &
                                   nNodeNeighbors_3d,   nhalo,           &
                                   active_cell,         beta_internal,   &
                                   xVertex,             yVertex,         &
                                   whichassemble_beta,                   &
                                   Auu(:,nz,:,:),       Avv(:,nz,:,:))

          endif   ! whichbabc

          !---------------------------------------------------------------------------
          ! Set rhs to the load vector
          ! The rhs can be adjusted below to account for inhomogeneous Dirichlet BC
          !---------------------------------------------------------------------------

          bu(:,:,:) = loadu(:,:,:)
          bv(:,:,:) = loadv(:,:,:)

          !---------------------------------------------------------------------------
          ! Incorporate Dirichlet boundary conditions (prescribed uvel and vvel)
          !---------------------------------------------------------------------------

          if (verbose_dirichlet .and. main_task) print*, 'Call Dirichlet_bc'

          call t_startf('glissade_dirichlet_3d')
          call dirichlet_boundary_conditions_3d(nx,              ny,                &
                                                nz,              nhalo,             &
                                                active_vertex,                      &
                                                umask_dirichlet, vmask_dirichlet,   &
                                                uvel,            vvel,              &
                                                Auu,             Auv,               &
                                                Avu,             Avv,               &
                                                bu,              bv)
          call t_startf('glissade_dirichlet_3d')
          
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
          !  because of roundoff.  Then we need to make sure that processor B's value 
          !  is communicated to processor A.  If these values are slightly different, 
          !  they will be reconciled by the subroutine check_symmetry_assembled_matrix.
          !---------------------------------------------------------------------------
          
          call t_startf('glissade_halo_Axxs')
          call staggered_parallel_halo(Auu(:,:,:,:))
          call staggered_parallel_halo(Auv(:,:,:,:))
          call staggered_parallel_halo(Avu(:,:,:,:))
          call staggered_parallel_halo(Avv(:,:,:,:))
          call t_stopf('glissade_halo_Axxs')
          
          !---------------------------------------------------------------------------
          ! Halo updates for rhs vectors
          ! (Not sure if these are necessary, but leaving them for now)
          !---------------------------------------------------------------------------

          call t_startf('glissade_halo_bxxs')
          call staggered_parallel_halo(bu(:,:,:))
          call staggered_parallel_halo(bv(:,:,:))
          call t_stopf('glissade_halo_bxxs')

          !---------------------------------------------------------------------------
          ! Check symmetry of assembled matrix
          ! 
          ! There may be small differences from perfect symmetry due to roundoff errors.  
          ! If sufficiently small, these differences are fixed by averaging the two values 
          !  that should be symmetric.  Otherwise the code aborts.
          !
          ! Note: It might be OK to skip this check for production code.  However,
          !       small violations of symmetry are not tolerated well by some solvers.
          !       For example, the SLAP PCG solver with incomplete Cholesky preconditioning
          !       can crash if symmetry is not perfect. 
          !---------------------------------------------------------------------------

          if (check_symmetry) then

             call t_startf('glissade_chk_symmetry')
             call check_symmetry_assembled_matrix_3d(nx,          ny,      &
                                                     nz,          nhalo,   &
                                                     active_vertex,        &
                                                     Auu,         Auv,     &
                                                     Avu,         Avv)
             call t_stopf('glissade_chk_symmetry')

          endif

          !---------------------------------------------------------------------------
          ! Count the total number of nonzero entries on all processors.
          !---------------------------------------------------------------------------

          call count_nonzeros_3d(nx,      ny,     &
                                 nz,      nhalo,  &
                                 Auu,     Auv,    &
                                 Avu,     Avv,    &
                                 active_vertex,   &
                                 nNonzeros)

          if (write_matrix) then
             if (counter == 1) then    ! first outer iteration only
 
                call t_startf('glissade_wrt_mat')
                call write_matrix_elements_3d(nx,          ny,         nz,         &
                                              nNodesSolve, nodeID,                 &
                                              iNodeIndex,  jNodeIndex, kNodeIndex, &
                                              Auu,         Auv,                    &
                                              Avu,         Avv,                    &
                                              bu,          bv)
                call t_stopf('glissade_wrt_mat')

             endif
          endif   ! write_matrix

          if (verbose_matrix .and. this_rank==rtest) then
             i = itest
             j = jtest
             k = ktest

             print*, ' '
             print*, 'i,j,k =', i, j, k
             print*, 'Auu sum =', sum(Auu(:,k,i,j))
             print*, 'Auv sum =', sum(Auv(:,k,i,j))
             print*, 'Avu sum =', sum(Avu(:,k,i,j))
             print*, 'Avv sum =', sum(Avv(:,k,i,j))

             print*, ' '
             print*, 'iA, jA, kA, Auu, Auv, Avu, Avv:'
             do kA = -1, 1
                do jA = -1, 1
                   do iA = -1, 1
                      m = indxA_3d(iA,jA,kA)
                      print*, iA, jA, kA, Auu(m,k,i,j), Auv(m,k,i,j), Avu(m,k,i,j), Avv(m,k,i,j) 
                   enddo
                enddo
             enddo
             
             print*, 'i, j, k: ', i, j, k
             print*, 'bu =', bu(k,i,j)
             print*, 'bv =', bv(k,i,j)
             
             j = jtest
             k = ktest
             m = indxA_3d(0,0,0)  ! diag entry
             print*, ' '
             print*, 'Matrix row properties, j, k =', j, k
             print*, ' '
             print*, 'i, diag, max, min, sum:'
             do i = 1, nx-1
                print*, ' '
                write(6,'(a4, i4, 4f16.8)') 'Auu:', i, Auu(m,k,i,j), maxval(Auu(:,k,i,j)), minval(Auu(:,k,i,j)), sum(Auu(:,k,i,j))
                write(6,'(a4, i4, 4f16.8)') 'Auv:', i, Auv(m,k,i,j), maxval(Auv(:,k,i,j)), minval(Auv(:,k,i,j)), sum(Auv(:,k,i,j))
             enddo
             
          endif  ! verbose_matrix

       endif  ! assemble 2d or 3d matrix

       !---------------------------------------------------------------------------
       ! If the matrix has no nonzero entries, then set velocities to zero and exit the solver.
       !---------------------------------------------------------------------------

       if (verbose_matrix .and. main_task) print*, 'nNonzeros in matrix =', nNonzeros

       if (nNonzeros == 0) then  ! clean up and return

          resid_u(:,:,:) = 0.d0
          resid_v(:,:,:) = 0.d0
          bu(:,:,:) = 0.d0
          bv(:,:,:) = 0.d0
          uvel(:,:,:) = 0.d0
          vvel(:,:,:) = 0.d0

          call t_startf('glissade_velo_higher_scale_outp')
          call glissade_velo_higher_scale_output(thck,    usrf,          &
                                                 topg,                   &
                                                 flwa,    efvs,          &
                                                 bwat,    mintauf,       &
                                                 beta_internal,          &
                                                 resid_u, resid_v,       &
                                                 bu,      bv,            &
                                                 uvel,    vvel,          &
                                                 uvel_2d, vvel_2d,       &
                                                 btractx, btracty,       &
                                                 taudx,   taudy,         &
                                                 tau_xz,  tau_yz,        &
                                                 tau_xx,  tau_yy,        &
                                                 tau_xy,  tau_eff)
          call t_stopf('glissade_velo_higher_scale_outp')
          
          if (main_task) print*, 'No nonzeros in matrix; exit glissade_velo_higher_solve'
          return

       endif  ! nNonzeros = 0

       !---------------------------------------------------------------------------
       ! Solve the 2D or 3D matrix system.
       !---------------------------------------------------------------------------

       !---------------------------------------------------------------------------
       ! First, handle a possible problem case: Set uvel_2d = vvel_2d = 0 for the case 
       !  of a Dirichlet no-slip basal BC and a 2D L1L2 solve.
       ! It would be pointless to apply the SSA to a no-slip problem, but this case
       !  is included for completeness.
       ! Note: DIVA computes a nonzero 2D velocity with a no-slip BC.
       !---------------------------------------------------------------------------

       if ((whichapprox==HO_APPROX_L1L2 .or. whichapprox==HO_APPROX_SSA) .and. &
              whichbabc==HO_BABC_NO_SLIP) then

          ! zero out velocity and related fields
          uvel_2d(:,:) = 0.d0
          vvel_2d(:,:) = 0.d0
          resid_u_2d(:,:) = 0.d0
          resid_v_2d(:,:) = 0.d0
          L2_norm = 0.d0   ! to force convergence on first step
          L2_norm_relative = 0.d0

       elseif (whichsparse == HO_SPARSE_PCG_STANDARD .or.   &
               whichsparse == HO_SPARSE_PCG_CHRONGEAR) then   ! native PCG solver
                                                              ! works for both serial and parallel runs

          !------------------------------------------------------------------------
          ! Compute the residual vector and its L2 norm
          !------------------------------------------------------------------------

          if (verbose_residual .and. main_task) then
             print*, 'Compute residual vector'
          endif

          if (solve_2d) then

             call t_startf('glissade_resid_vec')
             call compute_residual_vector_2d(nx,          ny,            &
                                             nhalo,                      &
                                             active_vertex,              &
                                             Auu_2d,      Auv_2d,        &
                                             Avu_2d,      Avv_2d,        &
                                             bu_2d,       bv_2d,         &
                                             uvel_2d,     vvel_2d,       &
                                             resid_u_2d,  resid_v_2d,    &
                                             L2_norm,     L2_norm_relative)
             call t_stopf('glissade_resid_vec')

             !------------------------------------------------------------------------
             ! Call linear PCG solver, compute uvel and vvel on local processor
             !------------------------------------------------------------------------

             !WHL - Passing itest, jtest, rtest for debugging

             call t_startf('glissade_pcg_slv_struct')

             if (whichsparse == HO_SPARSE_PCG_CHRONGEAR) then   ! use Chronopoulos-Gear PCG algorithm
                                                                ! (better scaling for large problems)
                call pcg_solver_chrongear_2d(nx,           ny,            &
                                             nhalo,                       &
                                             indxA_2d,     active_vertex, &
                                             Auu_2d,       Auv_2d,        &
                                             Avu_2d,       Avv_2d,        &
                                             bu_2d,        bv_2d,         &
                                             uvel_2d,      vvel_2d,       &
                                             whichprecond, err,           &
                                             niters,                      &
                                             itest, jtest, rtest, verbose_pcg)

             else   ! use standard PCG algorithm
             
                call pcg_solver_standard_2d(nx,           ny,            &
                                            nhalo,                       &
                                            indxA_2d,     active_vertex, &
                                            Auu_2d,       Auv_2d,        &
                                            Avu_2d,       Avv_2d,        &
                                            bu_2d,        bv_2d,         &
                                            uvel_2d,      vvel_2d,       &
                                            whichprecond, err,           &
                                            niters,                      &
                                            itest, jtest, rtest, verbose_pcg)

             endif  ! whichsparse

          else   ! 3D solve

             call t_startf('glissade_resid_vec')
             call compute_residual_vector_3d(nx,          ny,            &
                                             nz,          nhalo,         &
                                             active_vertex,              &
                                             Auu,         Auv,           &
                                             Avu,         Avv,           &
                                             bu,          bv,            &
                                             uvel,        vvel,          &
                                             resid_u,     resid_v,       &
                                             L2_norm,     L2_norm_relative)
             call t_stopf('glissade_resid_vec')

             !------------------------------------------------------------------------
             ! Call linear PCG solver, compute uvel and vvel on local processor
             !------------------------------------------------------------------------

             !WHL - Passing itest, jtest, rtest for debugging

             call t_startf('glissade_pcg_slv_struct')

             if (whichsparse == HO_SPARSE_PCG_CHRONGEAR) then   ! use Chronopoulos-Gear PCG algorithm
                                                                ! (better scaling for large problems)

                call pcg_solver_chrongear_3d(nx,           ny,            &
                                             nz,           nhalo,         &
                                             indxA_3d,     active_vertex, &
                                             Auu,          Auv,           &
                                             Avu,          Avv,           &
                                             bu,           bv,            &
                                             uvel,         vvel,          &
                                             whichprecond, err,           &
                                             niters,                      &
                                             itest, jtest, rtest, verbose_pcg)

             else   ! use standard PCG algorithm
             
                call pcg_solver_standard_3d(nx,           ny,            &
                                            nz,           nhalo,         &
                                            indxA_3d,     active_vertex, &
                                            Auu,          Auv,           &
                                            Avu,          Avv,           &
                                            bu,           bv,            &
                                            uvel,         vvel,          &
                                            whichprecond, err,           &
                                            niters,                      &
                                            itest, jtest, rtest, verbose_pcg)

             endif   ! whichsparse

          endif      ! whichapprox

          call t_stopf('glissade_pcg_slv_struct')

#ifdef TRILINOS
       elseif (whichsparse == HO_SPARSE_TRILINOS) then   ! solve with Trilinos

          !------------------------------------------------------------------------
          ! Compute the residual vector and its L2 norm
          !------------------------------------------------------------------------

          if (solve_2d) then

             if (verbose_residual .and. main_task) print*, 'Compute 2D residual vector'

             call t_startf('glissade_resid_vec')
             call compute_residual_vector_2d(nx,          ny,            &
                                             nhalo,                      &
                                             active_vertex,              &
                                             Auu_2d,      Auv_2d,        &
                                             Avu_2d,      Avv_2d,        &
                                             bu_2d,       bv_2d,         &
                                             uvel_2d,     vvel_2d,       &
                                             resid_u_2d,  resid_v_2d,    &
                                             L2_norm,     L2_norm_relative)
             call t_stopf('glissade_resid_vec')

             !------------------------------------------------------------------------
             ! Given Auu, bu, etc., assemble the matrix and RHS in a form
             ! suitable for Trilinos
             !------------------------------------------------------------------------

             if (verbose_trilinos .and. main_task) then
                print*, 'L2_norm, L2_target =', L2_norm, L2_target
                print*, 'Assemble matrix for Trilinos'
             endif

             call t_startf('glissade_trilinos_assemble')
             call trilinos_assemble_2d(nx,             ny,               &   
                                       nVerticesSolve, global_vertex_id, &
                                       iVertexIndex,   jVertexIndex,     &
                                       indxA_2d,       Afill_2d,         &
                                       Auu_2d,         Auv_2d,           &
                                       Avu_2d,         Avv_2d,           &
                                       bu_2d,          bv_2d)
             call t_stopf('glissade_trilinos_assemble')

             !------------------------------------------------------------------------
             ! Solve the linear matrix problem
             !------------------------------------------------------------------------

             if (verbose_trilinos .and. main_task) print*, 'Solve the matrix using Trilinos'

             call t_startf('glissade_vel_tgs')
             call solvevelocitytgs(velocityResult)
             call t_stopf('glissade_vel_tgs')

             !------------------------------------------------------------------------
             ! Put the velocity solution back into 2D arrays
             !------------------------------------------------------------------------

             call t_startf('glissade_trilinos_post')
             call trilinos_extract_velocity_2d(nx,            ny,           &
                                               nVerticesSolve,              &
                                               iVertexIndex,  jVertexIndex, &
                                               velocityResult,              &
                                               uvel_2d,       vvel_2d)
             call t_stopf('glissade_trilinos_post')

          else   ! 3D solve

             if (verbose_residual .and. main_task) print*, 'Compute 3D residual vector'

             call t_startf('glissade_resid_vec')
             call compute_residual_vector_3d(nx,          ny,            &
                                             nz,          nhalo,         &
                                             active_vertex,              &
                                             Auu,         Auv,           &
                                             Avu,         Avv,           &
                                             bu,          bv,            &
                                             uvel,        vvel,          &
                                             resid_u,     resid_v,       &
                                             L2_norm,     L2_norm_relative)
             call t_stopf('glissade_resid_vec')

             !------------------------------------------------------------------------
             ! Given Auu, bu, etc., assemble the matrix and RHS in a form
             ! suitable for Trilinos
             !------------------------------------------------------------------------

             if (verbose_trilinos .and. main_task) then
                print*, 'L2_norm, L2_target =', L2_norm, L2_target
                print*, 'Assemble matrix for Trilinos'
             endif

             call t_startf('glissade_trilinos_assemble')
             call trilinos_assemble_3d(nx,           ny,            nz,  &   
                                       nNodesSolve,  global_node_id,     &
                                       iNodeIndex,   jNodeIndex,    kNodeIndex,  &
                                       indxA_3d,     Afill,              &
                                       Auu,          Auv,                &
                                       Avu,          Avv,                &
                                       bu,           bv)
             call t_stopf('glissade_trilinos_assemble')

             !------------------------------------------------------------------------
             ! Solve the linear matrix problem
             !------------------------------------------------------------------------

             if (verbose_trilinos .and. main_task) print*, 'Solve the matrix using Trilinos'

             call t_startf('glissade_vel_tgs')
             call solvevelocitytgs(velocityResult)
             call t_stopf('glissade_vel_tgs')

             !------------------------------------------------------------------------
             ! Put the velocity solution back into 3D arrays
             !------------------------------------------------------------------------

             call t_startf('glissade_trilinos_post')
             call trilinos_extract_velocity_3d(nx,          ny,         nz,  &
                                               nNodesSolve,                  &
                                               iNodeIndex,  jNodeIndex, kNodeIndex, &
                                               velocityResult,               &
                                               uvel,        vvel)
             call t_stopf('glissade_trilinos_post')

          endif  ! whichapprox
#endif

       else   ! one-processor SLAP solve   
          
          !------------------------------------------------------------------------
          ! Given the stiffness matrices (Auu, etc.) and rhs vector (bu, bv) in
          !  structured format, form the global matrix and rhs in SLAP format.
          !------------------------------------------------------------------------

          if (verbose) print*, 'Form global matrix in SLAP sparse format'
 
          matrix%order = matrix_order
          matrix%nonzeros = max_nonzeros
          matrix%symmetric = .false.   ! Although the matrix is symmetric, we don't pass it to SLAP in symmetric form

          call t_startf('glissade_slap_preprocess')
          if (solve_2d) then

             call slap_preprocess_2d(nx,             ny,           &   
                                     nVerticesSolve, vertexID,     &
                                     iVertexIndex,   jVertexIndex, &
                                     indxA_2d,                     &
                                     Auu_2d,         Auv_2d,       &
                                     Avu_2d,         Avv_2d,       &
                                     bu_2d,          bv_2d,        &
                                     uvel_2d,        vvel_2d,      &
                                     matrix_order,                 &
                                     matrix,         rhs,          &
                                     answer)

          else   ! 3D solve

             call slap_preprocess_3d(nx,           ny,          nz, &   
                                     nNodesSolve,  nodeID,      &
                                     iNodeIndex,   jNodeIndex,  &
                                     kNodeIndex,   indxA_3d,    &
                                     Auu,          Auv,         &
                                     Avu,          Avv,         &
                                     bu,           bv,          &
                                     uvel,         vvel,        &
                                     matrix_order,              &
                                     matrix,       rhs,         &
                                     answer)

          endif  ! whichapprox
          call t_stopf('glissade_slap_preprocess')

          !------------------------------------------------------------------------
          ! Compute the residual vector and its L2_norm
          !------------------------------------------------------------------------

          call t_startf('glissade_slap_resid_vec')
          call slap_compute_residual_vector(matrix,  answer,    &
                                            rhs,     resid_vec, &
                                            L2_norm, L2_norm_relative)
          call t_stopf('glissade_slap_resid_vec')

          if (verbose_residual .and. main_task) then
             print*, 'L2_norm of residual =', L2_norm
          endif

          !------------------------------------------------------------------------
          ! Solve the linear matrix problem
          !------------------------------------------------------------------------

          call t_startf('glissade_easy_slv')
          call sparse_easy_solve(matrix, rhs,    answer,  &
                                 err,    niters, whichsparse)
          call t_stopf('glissade_easy_slv')

          !------------------------------------------------------------------------
          ! Put the velocity solution back into the uvel and vvel arrays
          !------------------------------------------------------------------------

          call t_startf('glissade_slap_post')

          if (solve_2d) then

             call slap_postprocess_2d(nVerticesSolve,              &
                                      iVertexIndex, jVertexIndex,  &
                                      answer,       resid_vec,     &
                                      uvel_2d,      vvel_2d,       &
                                      resid_u_2d,   resid_v_2d)

          else   ! 3D solve

             call slap_postprocess_3d(nNodesSolve,                            &
                                      iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                                      answer,       resid_vec,                &
                                      uvel,         vvel,                     &
                                      resid_u,      resid_v)

          endif   ! whichapprox

          call t_stopf('glissade_slap_post')

       endif   ! whichsparse 

       if (whichsparse /= HO_SPARSE_TRILINOS) then
          ! niters isn't set when using the trilinos solver
          if (main_task) then
             print*, 'Solved the linear system, niters, err =', niters, err
          endif
       end if

       if (solve_2d) then

          !------------------------------------------------------------------------
          ! Halo updates for uvel and vvel
          !------------------------------------------------------------------------

          call t_startf('glissade_halo_xvel')
          call staggered_parallel_halo(uvel_2d)
          call staggered_parallel_halo(vvel_2d)
          call t_stopf('glissade_halo_xvel')

          if (verbose_velo .and. this_rank==rtest) then
             i = itest
             j = jtest
             print*, ' '
             print*, 'rank, i, j, uvel_2d, vvel_2d (m/yr):', &
                      this_rank, i, j, uvel_2d(i,j), vvel_2d(i,j)               
          endif

          !---------------------------------------------------------------------------
          ! Compute residual quantities based on the velocity solution
          !---------------------------------------------------------------------------

          call t_startf('glissade_resid_vec2')
          call compute_residual_velocity_2d(nhalo,    whichresid,   &
                                            uvel_2d,  vvel_2d,      &
                                            usav_2d,  vsav_2d,      &
                                            resid_velo)
          call t_stopf('glissade_resid_vec2')

       else   ! 3D solve

          !------------------------------------------------------------------------
          ! Halo updates for uvel and vvel
          !------------------------------------------------------------------------

          call t_startf('glissade_halo_xvel')
          call staggered_parallel_halo(uvel)
          call staggered_parallel_halo(vvel)
          call t_stopf('glissade_halo_xvel')
          
          if (verbose_velo .and. this_rank==rtest) then
             i = itest
             j = jtest
             print*, ' '
             print*, 'rank, i, j:', this_rank, i, j
             print*, 'k, uvel, vvel:'
             do k = 1, nz
                print*, k, uvel(k,i,j), vvel(k,i,j)
             enddo
             print*, ' '
          endif

          !---------------------------------------------------------------------------
          ! Compute residual quantities based on the velocity solution
          !---------------------------------------------------------------------------

          call t_startf('glissade_resid_vec2')
          call compute_residual_velocity_3d(nhalo,  whichresid,   &
                                            uvel,   vvel,        &
                                            usav,   vsav,        &
                                            resid_velo)
          call t_stopf('glissade_resid_vec2')

       endif ! 2D or 3D solve

       !---------------------------------------------------------------------------
       ! Some calculations specific to the DIVA scheme
       !---------------------------------------------------------------------------

       if (whichapprox == HO_APPROX_DIVA) then

          ! Compute the components of basal traction, based on Goldberg (2011) eq. 38-39
          ! These are needed to compute the effective viscosity on the next iteration

          btractx(:,:) = beta_eff(:,:) * uvel_2d(:,:)
          btracty(:,:) = beta_eff(:,:) * vvel_2d(:,:)

          ! Interpolate omega_k to the staggered grid

          do k = 1, nz
             call glissade_stagger(nx,              ny,                           &
                                   omega_k(k,:,:),  stag_omega_k(k,:,:),  &
                                   ice_mask,        stagger_margin_in = 1)
          enddo

          ! Compute the new 3D velocity field
          ! NOTE: The full velocity field is not needed to update efvs and solve 
          !       again for uvel_2d and vvel_2D.  However, the basal velocity
          !       may be needed as an input to calcbeta.  It is possible to
          !       compute the basal velocity without computing the full column
          !       velocity, but it is simpler just to compute over the full column.

          call compute_3d_velocity_diva(nx,              ny,                   &
                                        nz,              sigma,                &
                                        active_vertex,   diva_level_index,     &
                                        stag_omega_k,    stag_omega,           &
                                        btractx,         btracty,              &
                                        uvel_2d,         vvel_2d,              &
                                        uvel,            vvel)

          call staggered_parallel_halo(uvel)
          call staggered_parallel_halo(vvel)

       endif   ! DIVA

       !---------------------------------------------------------------------------
       ! Write diagnostics (iteration number, max residual, and residual target
       !---------------------------------------------------------------------------

       if (main_task) then
          if (whichresid == HO_RESID_L2NORM) then
             print '(i4,2g20.6)', counter, L2_norm, L2_target
          elseif (whichresid == HO_RESID_L2NORM_RELATIVE) then
             print '(i4,2g20.6)', counter, L2_norm_relative, L2_target_relative
          else
             print '(i4,2g20.6)', counter, resid_velo, resid_target
          end if
       endif

       !---------------------------------------------------------------------------
       ! Update the outer loop stopping criterion
       !---------------------------------------------------------------------------

       if (whichresid == HO_RESID_L2NORM) then
          outer_it_criterion = L2_norm
          outer_it_target = L2_target           ! L2_target is currently set to 1.d-4 and held constant
       elseif (whichresid == HO_RESID_L2NORM_RELATIVE) then
          outer_it_criterion = L2_norm_relative
          outer_it_target = L2_target_relative  ! L2_target_relative is currently set to 1.d-7 and held constant
       else
          outer_it_criterion = resid_velo
          outer_it_target = resid_target   ! resid_target is currently a parameter = 1.d-4  
       end if

    enddo  ! while (outer_it_criterion >= outer_it_target .and. counter < maxiter_nonlinear)

    call t_stopf('glissade_vhs_nonlinear_loop')

    if (counter < maxiter_nonlinear) then
       converged_soln = .true.
!!       if (verbose .and. main_task) then
       if (main_task) then
          print*, ' '
          print*, 'GLISSADE SOLUTION HAS CONVERGED, outer counter =', counter
       endif
    else
       converged_soln = .false.
!!       if (verbose .and. main_task) then
       if (main_task) then
          print*, ' '
          print*, 'GLISSADE SOLUTION HAS NOT CONVERGED: counter, err =', counter, L2_norm
          !WHL - debug
!!          stop
       endif
    endif

    if (verbose_glp .and. this_rank==rtest) then 
       print*, ' '
       print*, 'beta_internal, rank =', rtest
       do j = jtest+1, jtest-1, -1
          do i = itest-3, itest+3
             write(6,'(f10.2)',advance='no') beta_internal(i,j)
          enddo
          print*, ' '
       enddo       
    endif

    !------------------------------------------------------------------------------
    ! After a 2D solve, fill in the full 3D velocity arrays.
    ! This is a simple copy for SSA, but required vertical integrals for L1L2 and DIVA. 
    ! Note: We store redundant 3D residual info rather than creating a separate 2D residual array.
    !------------------------------------------------------------------------------

    if (whichapprox == HO_APPROX_SSA) then ! fill the 3D velocity and residual arrays with the 2D values

       do k = 1, nz
          uvel(k,:,:) = uvel_2d(:,:)
          vvel(k,:,:) = vvel_2d(:,:)
          resid_u(k,:,:) = resid_u_2d(:,:)
          resid_v(k,:,:) = resid_v_2d(:,:)
       enddo

    elseif (whichapprox == HO_APPROX_L1L2) then

       if (verbose_L1L2 .and. main_task) print*, 'Compute 3D velocity, L1L2'

       uvel(nz,:,:) = uvel_2d(:,:)
       vvel(nz,:,:) = vvel_2d(:,:)
       do k = 1, nz
          resid_u(k,:,:) = resid_u_2d(:,:)
          resid_v(k,:,:) = resid_v_2d(:,:)
       enddo

       call compute_3d_velocity_L1L2(nx,               ny,              &
                                     nz,               sigma,           &
                                     dx,               dy,              &
                                     nhalo,                             &
                                     ice_mask,         land_mask,       &
                                     active_cell,      active_vertex,   &
                                     umask_dirichlet(nz,:,:),           &
                                     vmask_dirichlet(nz,:,:),           &
                                     xVertex,          yVertex,         &
                                     thck,             stagthck,        &
                                     usrf,                              &
                                     dusrf_dx,         dusrf_dy,        &
                                     flwa,             efvs,            &
                                     whichgradient_margin,              &
                                     max_slope,                         &
                                     uvel,             vvel)

       call staggered_parallel_halo(uvel)
       call staggered_parallel_halo(vvel)

    elseif (whichapprox == HO_APPROX_DIVA) then

       do k = 1, nz
          resid_u(k,:,:) = resid_u_2d(:,:)
          resid_v(k,:,:) = resid_v_2d(:,:)
       enddo

       !WHL - Commented out because the 3D velocity is now computed after each iteration.

!       ! Interpolate omega_k to the staggered grid

!       do k = 1, nz
!          call glissade_stagger(nx,              ny,                           &
!                                omega_k(k,:,:),  stag_omega_k(k,:,:),  &
!                                ice_mask,        stagger_margin_in = 1)
!       enddo

!       call compute_3d_velocity_diva(nx,              ny,                   &
!                                     nz,              sigma,                &
!                                     active_vertex,   diva_level_index,     &
!                                     stag_omega_k,    stag_omega,           &
!                                     btractx,         btracty,              &
!                                     uvel_2d,         vvel_2d,              &
!                                     uvel,            vvel)

!       call staggered_parallel_halo(uvel)
!       call staggered_parallel_halo(vvel)

       if (verbose_diva .and. this_rank==rtest) then
          print*, 'Computed 3D velocity, DIVA'
          i = itest
          j = jtest
          print*, ' '
          print*, 'i, j, beta, beta_eff:', i, j, beta_internal(i,j), beta_eff(i,j)
       endif

    endif   ! whichapprox

    !------------------------------------------------------------------------------
    ! Compute the components of the 3D stress tensor.
    ! These are diagnostic, except that tau_eff is used in the temperature calculation.
    !------------------------------------------------------------------------------

    call compute_internal_stress(nx,            ny,            &
                                 nz,            sigma,         &
                                 nhalo,         active_cell,   &
                                 xVertex,       yVertex,       &
                                 stagusrf,      stagthck,      &
                                 flwafact,      efvs,          &
                                 whichefvs,     efvs_constant, &
                                 whichapprox,                  &
                                 uvel,          vvel,          &
                                 tau_xz,        tau_yz,        &
                                 tau_xx,        tau_yy,        &
                                 tau_xy,        tau_eff)

    !------------------------------------------------------------------------------
    ! Compute the heat flux due to basal friction for each grid cell.
    !------------------------------------------------------------------------------

    call compute_basal_friction_heatflx(nx,            ny,            &
                                        nhalo,         active_cell,   &
                                        xVertex,       yVertex,       &
                                        uvel(nz,:,:),  vvel(nz,:,:),  &
                                        beta_internal, whichassemble_bfric,  &
                                        bfricflx)
                         
    !WHL - debug
    if (verbose_bfric .and. this_rank==rtest) then
       print*, ' '
       print*, 'Basal friction, itest, jtest, rank =', itest, jtest, rtest
!!          do j = ny-1, 1, -1
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
!!             do i = 1, nx-1
          do i = itest-3, itest+3
             write(6,'(e10.3)',advance='no') bfricflx(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    !------------------------------------------------------------------------------
    ! Compute the components of basal traction.
    !------------------------------------------------------------------------------

    btractx(:,:) = beta_internal(:,:) * uvel(nz,:,:)
    btracty(:,:) = beta_internal(:,:) * vvel(nz,:,:)

    ! Debug prints
    if (verbose_velo .and. this_rank==rtest) then
       print*, ' '
       print*, 'uvel, k=1 (m/yr):'
       do j = ny-nhalo, nhalo+1, -1
          do i = nhalo+1, nx-nhalo
             write(6,'(f8.2)',advance='no') uvel(1,i,j)
          enddo
          print*, ' '
       enddo

       print*, ' '
       print*, 'vvel, k=1 (m/yr):'
       do j = ny-nhalo, nhalo+1, -1
          do i = nhalo+1, nx-nhalo
             write(6,'(f8.2)',advance='no') vvel(1,i,j)
          enddo
          print*, ' '
       enddo       
       print*, ' '
       print*, 'max(uvel, vvel) =', maxval(uvel), maxval(vvel)
       print*, ' '

       i = itest
       j = jtest
       print*, 'New velocity: rank, i, j =', this_rank, i, j    
       print*, 'k, uvel, vvel:'
       do k = 1, nz
          print*, k, uvel(k,i,j), vvel(k,i,j)
       enddo
       if (solve_2d) print*, '2D velo:', uvel_2d(i,j), vvel_2d(i,j)
    endif  ! verbose_velo

    !------------------------------------------------------------------------------
    ! Clean up
    !------------------------------------------------------------------------------

    call t_startf('glissade_vhs_cleanup')
    if (whichsparse <= HO_SPARSE_GMRES) then  ! using SLAP solver
       deallocate(matrix%row, matrix%col, matrix%val)
       deallocate(rhs, answer, resid_vec)
    endif

#ifdef TRILINOS
    if (whichsparse == HO_SPARSE_TRILINOS) then
       deallocate(active_owned_unknown_map)
       deallocate(velocityResult)
       if (solve_2d) then
          deallocate(Afill_2d)
       else
          deallocate(Afill)
       endif
    endif
#endif

    if (solve_2d) then
       deallocate(Auu_2d, Auv_2d, Avu_2d, Avv_2d)
       deallocate(bu_2d, bv_2d)
       deallocate(loadu_2d, loadv_2d)
       deallocate(usav_2d, vsav_2d)
       deallocate(resid_u_2d, resid_v_2d)
    else
       deallocate(Auu, Auv, Avu, Avv)
    endif

    !------------------------------------------------------------------------------
    ! Convert output variables to appropriate CISM units (generally dimensionless).
    ! Note: bfricflx already has the desired units (W/m^2).
    !------------------------------------------------------------------------------

!pw call t_startf('glissade_velo_higher_scale_output')
    call glissade_velo_higher_scale_output(thck,    usrf,          &
                                           topg,                   &
                                           flwa,    efvs,          &
                                           bwat,    mintauf,       &
                                           beta_internal,          &
                                           resid_u, resid_v,       &
                                           bu,      bv,            &
                                           uvel,    vvel,          &
                                           uvel_2d, vvel_2d,       &
                                           btractx, btracty,       &
                                           taudx,   taudy,         &
                                           tau_xz,  tau_yz,        &
                                           tau_xx,  tau_yy,        &
                                           tau_xy,  tau_eff)
!pw call t_stopf('glissade_velo_higher_scale_output')
    call t_stopf('glissade_vhs_cleanup')

  end subroutine glissade_velo_higher_solve

!****************************************************************************

  subroutine glissade_velo_higher_scale_input(dx,      dy,            &
                                              thck,    usrf,          &
                                              topg,                   &
                                              eus,     thklim,        &
                                              flwa,    efvs,          &
                                              bwat,    mintauf,       &
                                              ho_beta_const,          &
                                              beta_grounded_min,      &
                                              btractx, btracty,       &
                                              uvel,    vvel,          &
                                              uvel_2d, vvel_2d)

    !--------------------------------------------------------
    ! Convert input variables (generally dimensionless)
    ! to appropriate units for the Glissade solver.
    !--------------------------------------------------------

    real(dp), intent(inout) ::   &
       dx, dy                  ! grid cell length and width 

    real(dp), dimension(:,:), intent(inout) ::   &
       thck,                &  ! ice thickness
       usrf,                &  ! upper surface elevation
       topg                    ! elevation of topography

    real(dp), intent(inout) ::   &
       eus,  &                 ! eustatic sea level (= 0 by default)
       thklim,  &              ! minimum ice thickness for active cells
       ho_beta_const, &        ! constant beta value (Pa/(m/yr)) for whichbabc = HO_BABC_CONSTANT
       beta_grounded_min       ! minimum beta value (Pa/(m/yr)) for grounded ice

    real(dp), dimension(:,:,:), intent(inout) ::  &
       flwa,   &               ! flow factor in units of Pa^(-n) yr^(-1)
       efvs                    ! effective viscosity (Pa yr)

    real(dp), dimension(:,:), intent(inout)  ::  &
       bwat,  &                ! basal water depth (m)
       mintauf, &              ! till yield stress (Pa)
       btractx, btracty,  &    ! components of basal traction (Pa)
       uvel_2d, vvel_2d        ! components of 2D velocity (m/yr)

    real(dp), dimension(:,:,:), intent(inout) ::  &
       uvel, vvel              ! components of 3D velocity (m/yr)

    ! grid cell dimensions: rescale from dimensionless to m
    dx = dx * len0
    dy = dy * len0

    ! ice geometry: rescale from dimensionless to m
    thck = thck * thk0
    usrf = usrf * thk0
    topg = topg * thk0
    eus  = eus  * thk0
    thklim = thklim * thk0

    ! rate factor: rescale from dimensionless to Pa^(-n) yr^(-1)
    flwa = flwa * (vis0*scyr)

    ! effective viscosity: rescale from dimensionless to Pa yr
    efvs = efvs * (evs0/scyr)

    ! bwat: rescale from dimensionless to m
    bwat = bwat * thk0

    ! mintauf: rescale from dimensionless to Pa
    mintauf = mintauf * tau0

    ! beta_parameters: rescale from dimensionless to Pa/(m/yr)
    ho_beta_const = ho_beta_const * tau0/(vel0*scyr)
    beta_grounded_min = beta_grounded_min * tau0/(vel0*scyr)

    ! basal traction: rescale from dimensionless to Pa
    btractx = btractx * tau0
    btracty = btracty * tau0

    ! ice velocity: rescale from dimensionless to m/yr
    uvel = uvel * (vel0*scyr)
    vvel = vvel * (vel0*scyr)
    uvel_2d = uvel_2d * (vel0*scyr)
    vvel_2d = vvel_2d * (vel0*scyr)

  end subroutine glissade_velo_higher_scale_input

!****************************************************************************

  subroutine glissade_velo_higher_scale_output(thck,    usrf,           &
                                               topg,                    &
                                               flwa,    efvs,           &                                       
                                               bwat,    mintauf,        &
                                               beta_internal,           &
                                               resid_u, resid_v,        &
                                               bu,      bv,             &
                                               uvel,    vvel,           &
                                               uvel_2d, vvel_2d,        &
                                               btractx, btracty,        &
                                               taudx,   taudy,          &
                                               tau_xz,  tau_yz,         &
                                               tau_xx,  tau_yy,         &
                                               tau_xy,  tau_eff)

    !--------------------------------------------------------
    ! Convert output variables to appropriate CISM units
    ! (generally dimensionless)
    !--------------------------------------------------------

    real(dp), dimension(:,:), intent(inout) ::  &
       thck,                 &  ! ice thickness
       usrf,                 &  ! upper surface elevation
       topg                     ! elevation of topography

    real(dp), dimension(:,:,:), intent(inout) ::  &
       flwa,   &                ! flow factor in units of Pa^(-n) yr^(-1)
       efvs                     ! effective viscosity (Pa yr)

    real(dp), dimension(:,:), intent(inout)  ::  &
       bwat,  &                 ! basal water depth (m)
       mintauf,  &              ! till yield stress (Pa)
       beta_internal            ! basal traction parameter (Pa/(m/yr))

    real(dp), dimension(:,:,:), intent(inout) ::  &
       uvel, vvel,    &         ! components of 3D velocity (m/yr)
       resid_u, resid_v,  &     ! components of residual Ax - b (Pa/m)
       bu, bv                   ! components of b in Ax = b (Pa/m)

    real(dp), dimension(:,:), intent(inout) ::  &
       uvel_2d, vvel_2d,       &! components of 2D velocity (m/yr)
       btractx, btracty,       &! components of basal traction (Pa)
       taudx,   taudy           ! components of driving stress (Pa)

    real(dp), dimension(:,:,:), intent(inout) ::  &
       tau_xz, tau_yz,         &! vertical components of stress tensor (Pa)
       tau_xx, tau_yy, tau_xy, &! horizontal components of stress tensor (Pa)
       tau_eff                  ! effective stress (Pa)

    ! Convert geometry variables from m to dimensionless units
    thck = thck / thk0
    usrf = usrf / thk0
    topg = topg / thk0

    ! Convert flow factor from Pa^(-n) yr^(-1) to dimensionless units
    flwa = flwa / (vis0*scyr)

    ! Convert effective viscosity from Pa yr to dimensionless units
    efvs = efvs / (evs0/scyr)

    ! Convert bwat from m to dimensionless units
    bwat = bwat / thk0

    ! Convert mintauf from Pa to dimensionless units
    mintauf = mintauf / tau0

    ! Convert beta_internal from Pa/(m/yr) to dimensionless units
    beta_internal = beta_internal / (tau0/(vel0*scyr))

    ! Convert velocity from m/yr to dimensionless units
    uvel = uvel / (vel0*scyr)
    vvel = vvel / (vel0*scyr)
    uvel_2d = uvel_2d / (vel0*scyr)
    vvel_2d = vvel_2d / (vel0*scyr)

    ! Convert residual and rhs from Pa/m to dimensionless units
    resid_u = resid_u / (tau0/len0)
    resid_v = resid_v / (tau0/len0)
    bu = bu / (tau0/len0)
    bv = bv / (tau0/len0)

    ! Convert stresses from Pa to dimensionless units
    btractx = btractx/tau0
    btracty = btracty/tau0
    taudx   = taudx/tau0
    taudy   = taudy/tau0
    tau_xz  = tau_xz/tau0
    tau_yz  = tau_yz/tau0
    tau_xx  = tau_xx/tau0
    tau_yy  = tau_yy/tau0
    tau_xy  = tau_xy/tau0
    tau_eff = tau_eff/tau0

  end subroutine glissade_velo_higher_scale_output

!****************************************************************************

  subroutine get_vertex_geometry(nx,           ny,                   &   
                                 nz,           nhalo,                &
                                 dx,           dy,                   &
                                 ice_mask,                           &
                                 xVertex,      yVertex,              &
                                 active_cell,  active_vertex,        &
                                 nNodesSolve,  nVerticesSolve,       &
                                 nodeID,       vertexID,             & 
                                 iNodeIndex,   jNodeIndex,  kNodeIndex, &
                                 iVertexIndex, jVertexIndex)
                            
    !----------------------------------------------------------------
    ! Compute coordinates for each vertex.
    ! Identify and count the active cells and vertices for the finite-element calculations.
    ! Active cells include all cells that contain ice (thck > thklin) and border locally owned vertices.
    ! Active vertices include all vertices of active cells.
    !
    ! Also compute some indices needed for the SLAP and Trilinos solvers.
    !TODO - Move SLAP/Trilinos part to a different subroutine?
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx,  ny,              &    ! number of grid cells in each direction
       nz,                   &    ! number of vertical levels where velocity is computed
       nhalo                      ! number of halo layers

    real(dp), intent(in) ::  &
       dx,  dy                ! grid cell length and width (m)
                              ! assumed to have the same value for each grid cell

    integer, dimension(nx,ny), intent(in) ::  &
       ice_mask               ! = 1 for cells where ice is present (thk > thklim), else = 0

    real(dp), dimension(nx-1,ny-1), intent(out) :: &
       xVertex, yVertex       ! x and y coordinates of each vertex

    logical, dimension(nx,ny), intent(out) :: &
       active_cell            ! true for active cells 
                              ! (thck > thklim, bordering a locally owned vertex)

    logical, dimension(nx-1,ny-1), intent(out) :: &
       active_vertex          ! true for vertices of active cells

    ! The remaining input/output arguments are for the SLAP and Trilinos solvers

    integer, intent(out) :: &
       nNodesSolve,         & ! number of locally owned nodes where we solve for velocity
       nVerticesSolve         ! number of locally owned vertices where we solve for velocity

    integer, dimension(nz,nx-1,ny-1), intent(out) ::  &
       nodeID                 ! local ID for each node where we solve for velocity

    integer, dimension(nx-1,ny-1), intent(out) ::  &
       vertexID               ! local ID for each vertex where we solve for velocity

    integer, dimension((nx-1)*(ny-1)*nz), intent(out) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of nodes

    integer, dimension((nx-1)*(ny-1)), intent(out) ::   &
       iVertexIndex, jVertexIndex   ! i and j indices of vertices

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

    ! Identify the active cells.
    ! Include all cells that border locally owned vertices and contain ice.

    active_cell(:,:) = .false.

    do j = 1+nhalo, ny-nhalo+1
    do i = 1+nhalo, nx-nhalo+1
       if (ice_mask(i,j) == 1) then
          active_cell(i,j) = .true.
       endif
    enddo
    enddo

    ! Identify the active vertices
    ! Include all vertices of active cells

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
    !  when we call the SLAP or Trilinos solver (one processor only).
    ! It is not required by the native PCG solver.

    nVerticesSolve  = 0
    vertexID(:,:)   = 0
    iVertexIndex(:) = 0
    jVertexIndex(:) = 0

    nNodesSolve   = 0
    nodeID(:,:,:) = 0
    iNodeIndex(:) = 0
    jNodeIndex(:) = 0
    kNodeIndex(:) = 0

    do j = nhalo+1, ny-nhalo    ! locally owned vertices only
    do i = nhalo+1, nx-nhalo
       if (active_vertex(i,j)) then   ! all nodes in column are active
          nVerticesSolve = nVerticesSolve + 1
          vertexID(i,j) = nVerticesSolve     ! unique local index for each vertex
          iVertexIndex(nVerticesSolve) = i
          jVertexIndex(nVerticesSolve) = j
          do k = 1, nz               
             nNodesSolve = nNodesSolve + 1   
             nodeID(k,i,j) = nNodesSolve     ! unique local index for each node
             iNodeIndex(nNodesSolve) = i
             jNodeIndex(nNodesSolve) = j
             kNodeIndex(nNodesSolve) = k
           enddo   ! k
        endif      ! active vertex
    enddo          ! i
    enddo          ! j

    if (verbose .and. this_rank==rtest) then
       print*, ' '
       print*, 'nVerticesSolve, nNodesSolve =', nVerticesSolve, nNodesSolve
    endif

  end subroutine get_vertex_geometry

!****************************************************************************

  subroutine load_vector_gravity(nx,               ny,              &
                                 nz,               sigma,           &
                                 nhalo,            active_cell,     &
                                 xVertex,          yVertex,         &
                                 stagusrf,         stagthck,        &
                                 dusrf_dx,         dusrf_dy,        &
                                 whichassemble_taud,                &
                                 loadu,            loadv)

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical levels at which velocity is computed
                                     ! Note: the number of elements per column is nz-1
       nhalo                         ! number of halo layers

    real(dp), dimension(nz), intent(in) ::    &
       sigma                         ! sigma vertical coordinate

    logical, dimension(nx,ny), intent(in) ::  &
       active_cell                   ! true if cell contains ice and borders a locally owned vertex

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xVertex, yVertex     ! x and y coordinates of vertices

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf,       &  ! upper surface elevation on staggered grid (m)
       stagthck           ! ice thickness on staggered grid (m)

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       dusrf_dx,       &  ! upper surface elevation gradient on staggered grid (m/m)
       dusrf_dy

    integer, intent(in) :: &
       whichassemble_taud   ! = 0 for standard finite element computation of driving stress terms
                            ! = 1 for computation that uses only the local value of the driving stress at each node

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::  &
       loadu, loadv       ! load vector, divided into u and v components

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    real(dp), dimension(nNodesPerElement_3d) ::     &
       x, y, z,         & ! Cartesian coordinates of nodes
       dsdx, dsdy         ! upper surface elevation gradient at nodes

    real(dp)  ::   &
       detJ               ! determinant of Jacobian for the transformation
                          !  between the reference element and true element

    !Note - These are not currently used except as dummy arguments
    real(dp), dimension(nNodesPerElement_3d) ::   &
       dphi_dx_3d, dphi_dy_3d, dphi_dz_3d  ! derivatives of basis functions, evaluated at quad pts

    real(dp) ::    &
       dsdx_qp, dsdy_qp       ! upper surface elevation gradient at quad pt

    integer :: i, j, k, n, p

    integer :: iNode, jNode, kNode

    if (verbose_load) then
       print*, ' '
       print*, 'In load_vector_gravity: itest, jtest, ktest, rtest =', itest, jtest, ktest, rtest
    endif
                
    ! Sum over elements in active cells 
    ! Loop over all cells that border locally owned vertices

    do j = nhalo+1, ny-nhalo+1
    do i = nhalo+1, nx-nhalo+1
       
       if (active_cell(i,j)) then

          do k = 1, nz-1    ! loop over elements in this column 
                            ! assume k increases from upper surface to bed

             ! compute spatial coordinates and upper surface elevation gradient for each node

             do n = 1, nNodesPerElement_3d

                ! Determine (k,i,j) for this node
                ! The reason for the '7' is that node 7, in the NE corner of the upper layer, has index (k,i,j).
                ! Indices for other nodes are computed relative to this node.
                iNode = i + ishift(7,n)
                jNode = j + jshift(7,n)
                kNode = k + kshift(7,n)

                x(n) = xVertex(iNode,jNode)
                y(n) = yVertex(iNode,jNode)
                z(n) = stagusrf(iNode,jNode) - sigma(kNode)*stagthck(iNode,jNode)
                dsdx(n) = dusrf_dx(iNode,jNode)
                dsdy(n) = dusrf_dy(iNode,jNode)

                if (verbose_load .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
                   print*, ' '
                   print*, 'i, j, k, n, x, y, z, dsdx, dsdy:', i, j, k, n, x(n), y(n), z(n), dsdx(n), dsdy(n)
                endif

             enddo   ! nodes per element

             ! Loop over quadrature points for this element
   
             do p = 1, nQuadPoints_3d

                ! Evaluate detJ at the quadrature point.
                ! TODO: The derivatives are not used.  Make these optional arguments?
                !WHL - debug - Pass in i, j, k, and p for now

                call get_basis_function_derivatives_3d(x(:),          y(:),          z(:),                    &
                                                       dphi_dxr_3d(:,p), dphi_dyr_3d(:,p), dphi_dzr_3d(:,p),  &
                                                       dphi_dx_3d(:),    dphi_dy_3d(:),    dphi_dz_3d(:),     &
                                                       detJ , i, j, k, p   )

                ! Increment the load vector with the gravitational contribution from this quadrature point
                ! The standard finite-element treatment (HO_ASSEMBLE_TAUD_STANDARD) is to take a 
                !  phi-weighted sum over neighboring vertices.
                ! For local driving stress (HO_ASSEMBLE_TAUD_LOCAL), use the value at the nearest vertex.
                ! (Note that vertex numbering is the same as QP numbering, CCW from 1 to 4 on bottom face and from 5 to 8 on top face.)

                if (whichassemble_taud == HO_ASSEMBLE_TAUD_LOCAL) then

                   ! Determine (k,i,j) for the node nearest to this quadrature point
                   iNode = i + ishift(7,p)
                   jNode = j + jshift(7,p)
                   kNode = k + kshift(7,p)
         
                   ! Add the ds/dx and ds/dy terms to the load vector for this node                   
                   loadu(kNode,iNode,jNode) = loadu(kNode,iNode,jNode) - rhoi*grav * detJ/vol0 * dsdx(p)
                   loadv(kNode,iNode,jNode) = loadv(kNode,iNode,jNode) - rhoi*grav * detJ/vol0 * dsdy(p)

                   if (verbose_load .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest .and. p==ptest) then
                      print*, ' '
                      print*, 'n, delta(loadu), delta(loadv):', n, rhoi*grav*detJ/vol0 * dsdx_qp, &
                                                                   rhoi*grav*detJ/vol0 * dsdy_qp
                   endif

                else   ! standard FE assembly (HO_ASSEMBLE_TAUD_STANDARD)

                   ! Evaluate dsdx and dsdy at this quadrature point
                   dsdx_qp = 0.d0
                   dsdy_qp = 0.d0
                   do n = 1, nNodesPerElement_3d
                      dsdx_qp = dsdx_qp + phi_3d(n,p) * dsdx(n)
                      dsdy_qp = dsdy_qp + phi_3d(n,p) * dsdy(n)
                   enddo

                   if (verbose_load .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
                      print*, ' '
                      print*, 'Increment load vector, i, j, k, p =', i, j, k, p
                      print*, 'ds/dx, ds/dy =', dsdx_qp, dsdy_qp
                      print*, 'detJ/vol0 =', detJ/vol0
                      print*, 'detJ/vol0* (ds/dx, ds/dy) =', detJ/vol0*dsdx_qp, detJ/vol0*dsdy_qp
                   endif

                   ! Loop over the nodes of the element
                   do n = 1, nNodesPerElement_3d

                      ! Determine (k,i,j) for this node
                      iNode = i + ishift(7,n)
                      jNode = j + jshift(7,n)
                      kNode = k + kshift(7,n)
         
                      ! Add the ds/dx and ds/dy terms to the load vector for this node                   
                      loadu(kNode,iNode,jNode) = loadu(kNode,iNode,jNode) - &
                           rhoi*grav * wqp_3d(p) * detJ/vol0 * dsdx_qp * phi_3d(n,p)
                      loadv(kNode,iNode,jNode) = loadv(kNode,iNode,jNode) - &
                           rhoi*grav * wqp_3d(p) * detJ/vol0 * dsdy_qp * phi_3d(n,p)

                      if (verbose_load .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest .and. p==ptest) then
                         print*, ' '
                         print*, 'n, phi_3d(n), delta(loadu), delta(loadv):', n, phi_3d(n,p), &
                                  rhoi*grav*wqp_3d(p)*detJ/vol0 * dsdx_qp * phi_3d(n,p), &
                                  rhoi*grav*wqp_3d(p)*detJ/vol0 * dsdy_qp * phi_3d(n,p)
                      endif

                   enddo   ! nNodesPerElement_3d

                endif   ! whichassemble_taud

             enddo      ! nQuadPoints_3d

          enddo         ! k

       endif            ! active cell

    enddo               ! i
    enddo               ! j

  end subroutine load_vector_gravity

!****************************************************************************

  subroutine load_vector_lateral_bc(nx,               ny,              &
                                    nz,               sigma,           &
                                    nhalo,                             &
                                    floating_mask,    ocean_mask,      &
                                    active_cell,                       &
                                    xVertex,          yVertex,         &
                                    stagusrf,         stagthck,        &
                                    loadu,            loadv)

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical levels at which velocity is computed
                                     ! Note: the number of elements per column is nz-1
       nhalo                         ! number of halo layers

    real(dp), dimension(nz), intent(in) ::    &
       sigma                         ! sigma vertical coordinate

    logical, dimension(nx,ny), intent(in) ::  &
       active_cell                   ! true if cell contains ice and borders a locally owned vertex

    integer, dimension(nx,ny), intent(in) ::  &
       floating_mask,               &! = 1 if ice is present and is floating
       ocean_mask                    ! = 1 if topography is below sea level and ice is absent

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xVertex, yVertex     ! x and y coordinates of vertices

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf,       &  ! upper surface elevation on staggered grid (m)
       stagthck           ! ice thickness on staggered grid (m)

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::  &
       loadu, loadv       ! load vector, divided into u and v components

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j

    ! Sum over elements in active cells 
    ! Loop over cells that contain locally owned vertices
    ! NOTE: Lateral shelf BCs are currently applied only to floating ice.
    !       I tested them for land-terminating ice, including an outward pressure term from the ice
    !       (with no compensating ocean pressure).  This gave excessive margin velocities.
    !
    ! TODO: Generalize to include marine-based ice that borders the ocean but is not floating?

    do j = nhalo+1, ny-nhalo+1
    do i = nhalo+1, nx-nhalo+1
       
       if (verbose_shelf .and. i==itest .and. j==jtest .and. this_rank==rtest) then
          print*, 'i, j =', i, j
          print*, 'active =', active_cell(i,j)
          print*, 'floating_mask =', floating_mask(i,j)
          print*, 'ocean_mask (i-1:i,j)  =', ocean_mask(i-1:i, j) 
          print*, 'ocean_mask (i-1:i,j-1)=', ocean_mask(i-1:i, j-1) 
       endif

       if (floating_mask(i,j) == 1) then   ! ice is present and is floating

          if (ocean_mask(i-1,j) == 1) then ! compute lateral BC for west face
             call lateral_shelf_bc(nx,            ny,          &
                                   nz,            sigma,       &
                                   'west',                     &
                                   i,             j,           &
                                   stagusrf,      stagthck,    &
                                   xVertex,       yVertex,     &
                                   loadu,         loadv)
          endif

          if (ocean_mask(i+1,j) == 1) then ! compute lateral BC for east face
             call lateral_shelf_bc(nx,            ny,          &
                                   nz,            sigma,       &
                                   'east',                     &
                                   i,             j,           &
                                   stagusrf,      stagthck,    &
                                   xVertex,       yVertex,     &
                                   loadu,         loadv)
          endif

          if (ocean_mask(i,j-1) == 1) then ! compute lateral BC for south face
             call lateral_shelf_bc(nx,            ny,          &
                                   nz,            sigma,       &
                                   'south',                    &
                                   i,             j,           &
                                   stagusrf,      stagthck,    &
                                   xVertex,       yVertex,     &
                                   loadu,         loadv)
          endif

          if (ocean_mask(i,j+1) == 1) then ! compute lateral BC for north face
             call lateral_shelf_bc(nx,            ny,          &
                                   nz,            sigma,       &
                                   'north',                    &
                                   i,             j,           &
                                   stagusrf,      stagthck,    &
                                   xVertex,       yVertex,     &
                                   loadu,         loadv)
          endif

       endif      ! floating_mask

    enddo         ! i
    enddo         ! j

  end subroutine load_vector_lateral_bc

!****************************************************************************

  subroutine lateral_shelf_bc(nx,                  ny,           &
                              nz,                  sigma,        &
                              face,                              &
                              iCell,               jCell,        &
                              stagusrf,            stagthck,     &
                              xVertex,             yVertex,      &
                              loadu,               loadv)

    !----------------------------------------------------------------------------------
    ! Determine the contribution to the load vector from ice and water pressure at the
    !  vertical boundary between ice and ocean (or alternatively, from ice pressure alone
    !  at a vertical boundary between ice and air).
    !
    ! This subroutine computes the vertically averaged hydrostatic pressure at a vertical face
    !  associated with the grid cell column (iCell, jCell).
    !
    ! At a given point, this pressure is proportional to the difference between
    ! (1) the vertically averaged pressure exerted outward (toward the ocean) by the ice front
    ! (2) the vertically averaged pressure exerted by the ocean back toward the ice
    ! 
    ! (1) is given by p_out = 0.5*rhoi*grav*H
    ! (2) is given by p_in  = 0.5*rhoi*grav*H*(rhoi/rhoo) for a floating shelf
    !                       = 0.5*rhoo*grav*H*(1 - s/H)^2 for s <= H but ice not necessarily afloat
    !
    ! The second term goes to zero for a land-terminating cliff. 
    ! The two pressure terms are opposite in sign, so the net vertically averaged pressure,
    !  directed toward the ocean (or air), is given by
    ! 
    !                    p_av = 0.5*rhoi*grav*H*(1 - rhoi/rhoo) for a floating shelf
    !                           0.5*rhoi*grav*H - 0.5*rhoo*grav*H * (1 - min((s/H),1)^2 for ice not necessarily afloat
    !
    ! Here we sum over quadrature points for each ocean-bordering face of each element.
    ! The contribution from each quadrature point to node N is proportional to the product
    !
    !                    p_av(s,H) * detJ * phi(n,p)
    !
    ! where s and H are the surface elevation and ice thickness evaluated at that point,
    !  detJ is the determinant of the transformation linking the reference 2D element coordinates
    !  to the true coordinates at that point, and phi(n,p) is the basis function evaluated at that point.
    !
    !-----------------------------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical levels at which velocity is computed
                                     ! Note: the number of elements per column is nz-1
       iCell, jCell                  ! i and j indices for this cell

    character(len=*), intent(in) ::  &
       face                          ! 'north', 'south', 'east', or 'west'
 
    real(dp), dimension(nz), intent(in) ::    &
       sigma                         ! sigma vertical coordinate

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xVertex, yVertex   ! x and y coordinates of vertices

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf,       &  ! upper surface elevation on staggered grid (m)
       stagthck           ! ice thickness on staggered grid (m)

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::  &
       loadu, loadv          ! load vector, divided into u and v components

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    real(dp), dimension(nNodesPerElement_2d) ::     &
       x, y,               & ! local coordinates of nodes
       s,                  & ! upper surface elevation at nodes
       h                     ! ice thickness at nodes

    integer, dimension(nNodesPerElement_2d) ::     &
       iNode, jNode, kNode   ! global indices of each node

    !Note: These are not currently used except as dummy arguments
    real(dp), dimension(nNodesPerElement_2d) ::   &
       dphi_dx_2d, dphi_dy_2d    ! derivatives of basis functions, evaluated at quad pts

    real(dp)  ::        &
       h_qp,            & ! ice thickness at a given quadrature point (m)
       s_qp,            & ! ice surface elevation at a given quadrature point (m)
       p_av,            & ! net outward pressure from ice, p_out - p_in
       detJ               ! determinant of Jacobian for the transformation
                          !  between the reference element and true element

    integer :: k, n, p

    ! Compute nodal geometry in a local xy reference system
    ! Note: The local y direction is really the vertical direction
    !       The local x direction depends on the face (N/S/E/W)
    ! The diagrams below show the node indexing convention, along with the true
    !  directions for each face.  The true directions are mapped to local (x,y).

    iNode(:) = 0
    jNode(:) = 0

    if (face=='west') then

       !     4-----3       z
       !     |     |       ^
       !     |     |       |
       !     1-----2       ---> -y

       iNode(1) = iCell-1
       jNode(1) = jCell

       iNode(2) = iCell-1
       jNode(2) = jCell-1

       x(1) = yvertex(iNode(1), jNode(1))
       x(2) = yvertex(iNode(2), jNode(2))

    elseif (face=='east') then

       !     4-----3       z
       !     |     |       ^
       !     |     |       |
       !     1-----2       ---> y

       iNode(1) = iCell
       jNode(1) = jCell-1

       iNode(2) = iCell
       jNode(2) = jCell

       x(1) = yvertex(iNode(1), jNode(1))
       x(2) = yvertex(iNode(2), jNode(2))

    elseif (face=='south') then

       !     4-----3       z
       !     |     |       ^
       !     |     |       |
       !     1-----2       ---> x

       iNode(1) = iCell-1
       jNode(1) = jCell-1

       iNode(2) = iCell
       jNode(2) = jCell-1

       x(1) = xvertex(iNode(1), jNode(1))
       x(2) = xvertex(iNode(2), jNode(2))

    elseif (face=='north') then

       !     4-----3       z
       !     |     |       ^
       !     |     |       |
       !     1-----2       ---> -x

       iNode(1) = iCell
       jNode(1) = jCell

       iNode(2) = iCell-1
       jNode(2) = jCell

       x(1) = xvertex(iNode(1), jNode(1))
       x(2) = xvertex(iNode(2), jNode(2))

    endif

    iNode(3) = iNode(2)
    jNode(3) = jNode(2)

    iNode(4) = iNode(1)
    jNode(4) = jNode(1)

    x(3) = x(2)
    x(4) = x(1)

    s(1) = stagusrf(iNode(1), jNode(1))
    s(2) = stagusrf(iNode(2), jNode(2))
    s(3) = s(2)
    s(4) = s(1)

    h(1) = stagthck(iNode(1), jNode(1))
    h(2) = stagthck(iNode(2), jNode(2))
    h(3) = h(2)
    h(4) = h(1)

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
          !TODO - Modify this subroutine to return only detJ, and not the derivatives?

          if (verbose_shelf .and. this_rank==rtest .and. iCell==itest .and. jCell==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Get detJ, i, j, k, p =', iCell, jCell, k, p
             print*, 'x =', x(:)
             print*, 'y =', y(:)
             print*, 'dphi_dxr_2d =', dphi_dxr_2d(:,p)
             print*, 'dphi_dyr_2d =', dphi_dyr_2d(:,p)
          endif

          call get_basis_function_derivatives_2d(x(:),              y(:),               & 
                                                 dphi_dxr_2d(:,p),  dphi_dyr_2d(:,p),   &
                                                 dphi_dx_2d(:),     dphi_dy_2d(:),      &
                                                 detJ, iCell, jCell, p)

          ! For some faces, detJ is computed to be a negative number because the face is
          ! oriented opposite the x or y axis.  Fix this by taking the absolute value.

          detJ = abs(detJ)

          ! Evaluate the ice thickness and surface elevation at this quadrature point

          h_qp = 0.d0
          s_qp = 0.d0
          do n = 1, nNodesPerElement_2d
             h_qp = h_qp + phi_2d(n,p) * h(n)
             s_qp = s_qp + phi_2d(n,p) * s(n)
          enddo

          if (verbose_shelf .and. this_rank==rtest .and. iCell==itest .and. jCell==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Increment shelf load vector, i, j, face, k, p =', iCell, jCell, trim(face), k, p
             print*, 'h_qp, s_qp =', h_qp, s_qp
             print*, 'detJ/vol0 =', detJ/vol0
             print*, 'grav =', grav
          endif

          ! Increment the load vector with the shelf water pressure contribution from 
          !  this quadrature point.
          ! Increment loadu for east/west faces and loadv for north/south faces.

          ! This formula works for ice that either is floating or is partially submerged without floating
!!          p_av = 0.5d0*rhoi*grav*h_qp &                                   ! p_out
!!               - 0.5d0*rhoo*grav*h_qp * (1.d0 - min(s_qp/h_qp,1.d0))**2   ! p_in

          ! This formula works for floating ice.
          p_av = 0.5d0*rhoi*grav*h_qp * (1.d0 - rhoi/rhoo)

          if (trim(face) == 'west') then  ! net force in -x direction

             do n = 1, nNodesPerElement_2d
                loadu(kNode(n),iNode(n),jNode(n)) = loadu(kNode(n),iNode(n),jNode(n))    &
                                                  - p_av * wqp_2d(p) * detJ/vol0 * phi_2d(n,p)
             enddo

          elseif (trim(face) == 'east') then  ! net force in x direction

             do n = 1, nNodesPerElement_2d
                loadu(kNode(n),iNode(n),jNode(n)) = loadu(kNode(n),iNode(n),jNode(n))    &
                                                  + p_av * wqp_2d(p) * detJ/vol0 * phi_2d(n,p)
             enddo

          elseif (trim(face) == 'south') then  ! net force in -y direction

             do n = 1, nNodesPerElement_2d
                loadv(kNode(n),iNode(n),jNode(n)) = loadv(kNode(n),iNode(n),jNode(n))    &
                                                  - p_av * wqp_2d(p) * detJ/vol0 * phi_2d(n,p)
             enddo

          elseif (trim(face) == 'north') then  ! net force in y direction
 
             do n = 1, nNodesPerElement_2d
                loadv(kNode(n),iNode(n),jNode(n)) = loadv(kNode(n),iNode(n),jNode(n))    &
                                                  + p_av * wqp_2d(p) * detJ/vol0 * phi_2d(n,p)
             enddo

          endif   ! face = N/S/E/W

       enddo   ! nQuadPoints_2d

    enddo   ! k (element faces in column)

  end subroutine lateral_shelf_bc

!****************************************************************************

  subroutine assemble_stiffness_matrix_3d(nx,               ny,              &
                                          nz,               sigma,           &
                                          nhalo,            active_cell,     &
                                          xVertex,          yVertex,         &
                                          uvel,             vvel,            &
                                          stagusrf,         stagthck,        &
                                          flwafact,         whichapprox,     &
                                          efvs,             whichefvs,       &
                                          efvs_constant,                     &       
                                          Auu,              Auv,             &
                                          Avu,              Avv)

    !----------------------------------------------------------------
    ! Assemble the stiffness matrix A in the linear system Ax = b.
    ! This subroutine is called for each nonlinear iteration if
    !  we are iterating on the effective viscosity.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical levels at which velocity is computed
                                     ! Note: the number of elements per column is nz-1
       nhalo                         ! number of halo layers

    real(dp), dimension(nz), intent(in) ::    &
       sigma                         ! sigma vertical coordinate

    logical, dimension(nx,ny), intent(in) ::  &
       active_cell                   ! true if cell contains ice and borders a locally owned vertex

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xVertex, yVertex     ! x and y coordinates of vertices

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       uvel, vvel         ! velocity components (m/yr)

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf,       &  ! upper surface elevation on staggered grid (m)
       stagthck           ! ice thickness on staggered grid (m)

    real(dp), dimension(nz-1,nx,ny), intent(in) ::  &
       flwafact           ! temperature-based flow factor, 0.5 * A^(-1/n), 
                          ! used to compute the effective viscosity
                          ! units: Pa yr^(1/n)

    integer, intent(in) ::   &
       whichapprox,     & ! option for Stokes approximation (BP, SSA, SIA)
       whichefvs          ! option for effective viscosity calculation 

    real(dp), intent(in) :: &
       efvs_constant      ! constant value of effective viscosity (Pa yr)

    real(dp), dimension(nz-1,nx,ny), intent(out) ::  &
       efvs               ! effective viscosity (Pa yr)

    real(dp), dimension(nNodeNeighbors_3d,nz,nx-1,ny-1), intent(out) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv                                    

    !---------------------------------------------------------
    ! Local variables
    !---------------------------------------------------------

    real(dp), dimension(nQuadPoints_3d) ::   &
       detJ               ! determinant of J

    real(dp), dimension(nNodesPerElement_3d) ::   &
       dphi_dx_3d, dphi_dy_3d, dphi_dz_3d  ! derivatives of basis function, evaluated at quad pt

    !----------------------------------------------------------------
    ! Note: Kuu, Kuv, Kvu, and Kvv are 8x8 components of the stiffness matrix
    !       for the local element.  (The combined stiffness matrix is 16x16.)
    !
    ! Once these matrices are formed, their coefficients are summed into the assembled
    !  matrices Auu, Auv, Avu, Avv.  The A matrices each have as many rows as there are
    !  active nodes, but only 27 columns, corresponding to the 27 vertices that belong to
    !  the 8 elements sharing a given node.
    !
    ! The native structured PCG solver works with the dense A matrices in the form
    ! computed here.  For the SLAP solver, the terms of the A matrices are put
    ! in a sparse matrix during preprocessing.  For the Trilinos solver, the terms
    ! of the A matrices are passed to Trilinos one row at a time. 
    !----------------------------------------------------------------

    real(dp), dimension(nNodesPerElement_3d, nNodesPerElement_3d) ::   &   !
       Kuu,          &    ! element stiffness matrix, divided into 4 parts as shown below
       Kuv,          &    !  
       Kvu,          &    !
       Kvv                !    Kuu  | Kuv
                          !    _____|____          
                          !         |
                          !    Kvu  | Kvv
                          !         
                          ! Kvu may not be needed if matrix is symmetric, but is included for now

    real(dp), dimension(nNodesPerElement_3d) ::     &
       x, y, z,         & ! Cartesian coordinates of nodes
       u, v,            & ! u and v at nodes
       s                  ! upper surface elevation at nodes

    real(dp), dimension(nQuadPoints_3d) ::    &
       efvs_qp            ! effective viscosity at a quad pt

    logical, parameter ::   &
       check_symmetry_element = .true.  ! if true, then check symmetry of element matrix
                                        !Note: Can speed up assembly a bit by setting to false for production

    integer :: i, j, k, n, p
    integer :: iNode, jNode, kNode

    if (verbose_matrix .and. main_task) then
       print*, ' '
       print*, 'In assemble_stiffness_matrix_3d'
       print*, 'itest, jtest, ktest, rtest =', itest, jtest, ktest, rtest
    endif

    ! Initialize effective viscosity
    efvs(:,:,:) = 0.d0

    ! Initialize global stiffness matrix

    Auu(:,:,:,:) = 0.d0
    Auv(:,:,:,:) = 0.d0
    Avu(:,:,:,:) = 0.d0
    Avv(:,:,:,:) = 0.d0

    ! Sum over elements in active cells 
    ! Loop over all cells that border locally owned vertices.

    do j = nhalo+1, ny-nhalo+1
    do i = nhalo+1, nx-nhalo+1

       if (active_cell(i,j)) then

          !WHL - debug
!!          print*, 'i, j:', i, j

          do k = 1, nz-1    ! loop over elements in this column 
                            ! assume k increases from upper surface to bed

             ! Initialize element stiffness matrix
             Kuu(:,:) = 0.d0
             Kuv(:,:) = 0.d0
             Kvu(:,:) = 0.d0
             Kvv(:,:) = 0.d0
  
             ! compute spatial coordinates, velocity, and upper surface elevation for each node

             do n = 1, nNodesPerElement_3d

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

                if (verbose_matrix .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
                   print*, ' '
                   print*, 'i, j, k, n, x, y, z:', i, j, k, n, x(n), y(n), z(n)
                   print*, 's, u, v:', s(n), u(n), v(n)
                endif

             enddo   ! nodes per element

             ! Loop over quadrature points for this element
   
             do p = 1, nQuadPoints_3d

                ! Evaluate the derivatives of the element basis functions at this quadrature point.
                !WHL - Pass in i, j, k, and p to the following subroutines for debugging.

                call get_basis_function_derivatives_3d(x(:),             y(:),             z(:),              &          
                                                       dphi_dxr_3d(:,p), dphi_dyr_3d(:,p), dphi_dzr_3d(:,p),  &
                                                       dphi_dx_3d(:),    dphi_dy_3d(:),    dphi_dz_3d(:),     &
                                                       detJ(p) , i, j, k, p  )

!          call t_startf('glissade_effective_viscosity')
                call compute_effective_viscosity(whichefvs,        whichapprox,                       &
                                                 efvs_constant,    nNodesPerElement_3d,               &
                                                 dphi_dx_3d(:),    dphi_dy_3d(:),    dphi_dz_3d(:),   &
                                                 u(:),             v(:),                              & 
                                                 flwafact(k,i,j),  efvs_qp(p),                        &
                                                 i, j, k, p )
!          call t_stopf('glissade_effective_viscosity')

                if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest .and. p==ptest) then
                   print*, 'i, j, k, p, efvs (Pa yr):', i, j, k, p, efvs_qp(p)
                endif

                ! Increment the element stiffness matrix with the contribution from each quadrature point.

!          call t_startf('glissade_compute_element_matrix')
                call compute_element_matrix(whichapprox,     nNodesPerElement_3d,               & 
                                            wqp_3d(p),       detJ(p),          efvs_qp(p),      &
                                            dphi_dx_3d(:),   dphi_dy_3d(:),    dphi_dz_3d(:),   &
                                            Kuu(:,:),        Kuv(:,:),                          &
                                            Kvu(:,:),        Kvv(:,:),                          &
                                            i, j, k, p )
!          call t_stopf('glissade_compute_element_matrix')

             enddo   ! nQuadPoints_3d

             ! Compute average of effective viscosity over quad pts
             efvs(k,i,j) = 0.d0

             do p = 1, nQuadPoints_3d
                efvs(k,i,j) = efvs(k,i,j) + efvs_qp(p)
             enddo
             efvs(k,i,j) = efvs(k,i,j) / nQuadPoints_3d
             
             if (check_symmetry_element) then
                call check_symmetry_element_matrix(nNodesPerElement_3d,  &
                                                   Kuu, Kuv, Kvu, Kvv)
             endif

             ! Sum terms of element matrix K into dense assembled matrix A

             call element_to_global_matrix_3d(nx,           ny,          nz, &
                                              i,            j,           k,  &
                                              Kuu,          Kuv,             &
                                              Kvu,          Kvv,             &
                                              Auu,          Auv,             &
                                              Avu,          Avv)

          enddo   ! nz  (loop over elements in this column)

          if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest) then
             print*, ' '
             print*, 'Assembled 3D matrix, i, j =', i, j
             print*, 'k, efvs:'
             do k = 1, nz-1
                print*, k, efvs(k,i,j)
             enddo
          endif

       endif   ! active cell

    enddo      ! i
    enddo      ! j

  end subroutine assemble_stiffness_matrix_3d

!****************************************************************************

  subroutine assemble_stiffness_matrix_2d(nx,               ny,              &
                                          nz,                                &
                                          sigma,            stagsigma,       &
                                          nhalo,            active_cell,     &
                                          xVertex,          yVertex,         &
                                          uvel_2d,          vvel_2d,         &
                                          stagusrf,         stagthck,        &
                                          flwa,             flwafact,        &
                                          whichapprox,                       &
                                          whichefvs,        efvs_constant,   &
                                          efvs,                              &
                                          Auu,              Auv,             &
                                          Avu,              Avv,             &
                                          dusrf_dx,         dusrf_dy,        &
                                          thck,                              &
                                          btractx,          btracty,         &
                                          omega_k,          omega,   &
                                          efvs_qp_3d)
  
    !----------------------------------------------------------------
    ! Assemble the stiffness matrix A in the linear system Ax = b.
    ! This subroutine is called for each nonlinear iteration if
    !  we are iterating on the effective viscosity.
    ! The matrix A can be based on the shallow-shelf approximation or 
    !  the depth-integrated L1L2 approximation (Schoof and Hindmarsh, 2010).
    !----------------------------------------------------------------
 
    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical levels at which velocity is computed
                                     ! (used for flwafact)
       nhalo                         ! number of halo layers

    real(dp), dimension(nz), intent(in) ::    &
       sigma              ! sigma vertical coordinate

    real(dp), dimension(nz-1), intent(in) ::    &
       stagsigma          ! staggered sigma vertical coordinate

    logical, dimension(nx,ny), intent(in) ::  &
       active_cell        ! true if cell contains ice and borders a locally owned vertex

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xVertex, yVertex   ! x and y coordinates of vertices

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       uvel_2d, vvel_2d   ! 2D velocity components (m/yr)

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf,       &  ! upper surface elevation on staggered grid (m)
       stagthck           ! ice thickness on staggered grid (m)

    !TODO - Pass in flwa and compute flwafact here?
    real(dp), dimension(nz-1,nx,ny), intent(in) ::  &
       flwa,             &! temperature-based flow factor A, Pa^{-n} yr^{-1}
       flwafact           ! temperature-based flow factor, 0.5 * A^(-1/n), Pa yr^(1/n)
                          ! used to compute the effective viscosity

    integer, intent(in) ::   &
       whichapprox,     & ! option for Stokes approximation (BP, L1L2, SSA, SIA)
       whichefvs          ! option for effective viscosity calculation 

    real(dp), intent(in) :: &
       efvs_constant      ! constant value of effective viscosity (Pa yr)

    real(dp), dimension(nz-1,nx,ny), intent(out) ::  &
       efvs               ! effective viscosity (Pa yr)

    real(dp), dimension(nNodeNeighbors_2d,nx-1,ny-1), intent(out) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv                                    

    ! The following optional arguments are used for the L1L2 approximation only

    real(dp), dimension(nx-1,ny-1), intent(in), optional ::  &
       dusrf_dx,       &  ! upper surface elevation gradient on staggered grid (m/m)
       dusrf_dy           ! needed for L1L2 assembly only

    ! The following optional arguments are used for DIVA only

    real(dp), dimension(nx,ny), intent(in), optional ::   &
       thck               ! ice thickness (m)

    real(dp), dimension(nx-1,ny-1), intent(in), optional ::   &
       btractx, btracty         ! components of basal traction (Pa)

    real(dp), dimension(nz,nx,ny), intent(out), optional :: &
       omega_k            ! single integral, defined by Goldberg (2011) eq. 32

    real(dp), dimension(nx,ny), intent(out), optional :: &
       omega              ! double integral, defined by Goldberg (2011) eq. 35
                          ! Note: omega here = Goldberg's omega/H

    real(dp), dimension(nz-1,nQuadPoints_2d,nx,ny), intent(inout), optional ::  &
       efvs_qp_3d         ! effective viscosity (Pa yr)

    !---------------------------------------------------------
    ! Local variables
    !---------------------------------------------------------

    real(dp), dimension(nQuadPoints_2d) ::   &
       detJ               ! determinant of J

    real(dp), dimension(nNodesPerElement_2d) ::   &
       dphi_dx_2d, dphi_dy_2d, dphi_dz_2d  ! derivatives of basis function, evaluated at quad pts
                                           ! set dphi_dz = 0 for 2D problem

    !----------------------------------------------------------------
    ! Note: Kuu, Kuv, Kvu, and Kvv are 4x4 components of the stiffness matrix
    !       for the local element.  (The combined stiffness matrix is 8x8.)
    !
    ! Once these matrices are formed, their coefficients are summed into the global
    !  matrices Auu_2d, Auv_2d, Avu_2d, Avv_2d.  The global matrices each have as 
    !  many rows as there are active vertices, but only 9 columns, corresponding to 
    !  the 9 vertices of the 4 elements sharing a given node.
    !
    ! The native structured PCG solver works with the dense A matrices in the form
    ! computed here.  For the SLAP solver, the terms of the A matrices are put
    ! in a sparse matrix format during preprocessing.  For the Trilinos solver, 
    ! the terms of the A matrices are passed to Trilinos one row at a time. 
    !----------------------------------------------------------------

    real(dp), dimension(nNodesPerElement_2d, nNodesPerElement_2d) ::   &   !
       Kuu,          &    ! element stiffness matrix, divided into 4 parts as shown below
       Kuv,          &    !  
       Kvu,          &    !
       Kvv                !    Kuu  | Kuv
                          !    _____|____          
                          !         |
                          !    Kvu  | Kvv
                          !         
                          ! Kvu may not be needed if matrix is symmetric, but is included for now

    real(dp), dimension(nNodesPerElement_2d) ::     &
       x, y,            & ! Cartesian coordinates of vertices
       u, v,            & ! depth-integrated mean velocity at vertices (m/yr)
       h,               & ! thickness at vertices (m)
       s,               & ! upper surface elevation at vertices (m)
       bx, by,          & ! basal traction at vertices (Pa) (DIVA only)
       dsdx, dsdy         ! upper surface elevation gradient at vertices (m/m) (L1L2 only)

    real(dp), dimension(nQuadPoints_2d) ::    &
       efvs_qp_vertavg    ! vertically averaged effective viscosity at a quad pt

    real(dp) ::         &
       h_qp               ! thickness at a quad pt

    real(dp), dimension(nz-1,nQuadPoints_2d) ::    &
       efvs_qp            ! effective viscosity at each layer in a cell column
                          ! corresponding to a quad pt

    logical, parameter ::   &
       check_symmetry_element = .true.  ! if true, then check symmetry of element matrix

    real(dp), dimension(nx,ny) ::  &
       flwafact_2d        ! vertically averaged flow factor

    integer :: i, j, k, n, p
    integer :: iVertex, jVertex

    if (verbose_matrix .and. main_task) then
       print*, ' '
       print*, 'In assemble_stiffness_matrix_2d'
    endif

    ! Initialize effective viscosity
    efvs(:,:,:) = 0.d0

    ! Initialize global stiffness matrix

    Auu(:,:,:) = 0.d0
    Auv(:,:,:) = 0.d0
    Avu(:,:,:) = 0.d0
    Avv(:,:,:) = 0.d0

    ! Compute vertical average of flow factor (SSA only)
    if (whichapprox == HO_APPROX_SSA) then
       call glissade_vertical_average(nx,       ny,      &
                                      nz,       sigma,   &
                                      active_cell,       &
                                      flwafact, flwafact_2d)
    endif

    ! Sum over elements in active cells 
    ! Loop over all cells that border locally owned vertices.

    do j = nhalo+1, ny-nhalo+1
    do i = nhalo+1, nx-nhalo+1
       
       if (active_cell(i,j)) then

          ! Initialize element stiffness matrix
          Kuu(:,:) = 0.d0
          Kuv(:,:) = 0.d0
          Kvu(:,:) = 0.d0
          Kvv(:,:) = 0.d0
  
          ! Compute spatial coordinates, velocity, thickness and surface elevation for each vertex
          ! Also compute surface elevation gradient (for L1L2) and basal traction (for DIVA)
          do n = 1, nNodesPerElement_2d

             ! Determine (i,j) for this vertex
             ! The reason for the '3' is that node 3, in the NE corner of the grid cell, has index (i,j).
             ! Indices for other nodes are computed relative to this vertex.
             iVertex = i + ishift(3,n)
             jVertex = j + jshift(3,n)

             x(n) = xVertex(iVertex,jVertex)
             y(n) = yVertex(iVertex,jVertex)
             u(n) = uvel_2d(iVertex,jVertex)
             v(n) = vvel_2d(iVertex,jVertex)
             s(n) = stagusrf(iVertex,jVertex)
             h(n) = stagthck(iVertex,jVertex)
             if (present(dusrf_dx) .and. present(dusrf_dy)) then  ! L1L2
                dsdx(n) = dusrf_dx(iVertex,jVertex)
                dsdy(n) = dusrf_dy(iVertex,jVertex)
             endif
             if (present(btractx) .and. present(btracty)) then    ! DIVA
                bx(n) = btractx(iVertex,jVertex)
                by(n) = btracty(iVertex,jVertex)
             endif

             if (verbose_matrix .and. this_rank==rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'i, j, n, x, y:', i, j, n, x(n), y(n)
                print*, 's, h, u, v:', s(n), h(n), u(n), v(n)
                if (present(btractx) .and. present(btracty)) print*, 'bx, by:', bx(n), by(n)
             endif

          enddo   ! vertices per element

          ! Loop over quadrature points for this element
   
          do p = 1, nQuadPoints_2d

             ! Evaluate the derivatives of the element basis functions at this quadrature point.

             !WHL - Pass in i, j and p to the following subroutines for debugging

             call get_basis_function_derivatives_2d(x(:),             y(:),          &
                                                    dphi_dxr_2d(:,p), dphi_dyr_2d(:,p), &
                                                    dphi_dx_2d(:),    dphi_dy_2d(:),    &
                                                    detJ(p) , i, j, p)
             dphi_dz_2d(:) = 0.d0

             if (whichapprox == HO_APPROX_L1L2) then

                ! Compute effective viscosity for each layer at this quadrature point
                !TODO - sigma -> stagsigma for L1L2 viscosity?
                call compute_effective_viscosity_L1L2(whichefvs,            efvs_constant,     &
                                                      nz,                   sigma,             &
                                                      nNodesPerElement_2d,  phi_2d(:,p),       &
                                                      dphi_dx_2d(:),        dphi_dy_2d(:),     &
                                                      u(:),                 v(:),              & 
                                                      h(:),                                    &
                                                      dsdx(:),              dsdy(:),           &
                                                      flwa(:,i,j),          flwafact(:,i,j),   &
                                                      efvs_qp(:,p),                            &
                                                      i, j, p)

                ! Compute vertical average of effective viscosity
                efvs_qp_vertavg(p) = 0.d0
                do k = 1, nz-1
                   efvs_qp_vertavg(p) = efvs_qp_vertavg(p) + efvs_qp(k,p) * (sigma(k+1) - sigma(k))
                enddo

             elseif (whichapprox == HO_APPROX_DIVA) then

                ! Copy efvs_qp from global array to local column array
                efvs_qp(:,:) = efvs_qp_3d(:,:,i,j)

                ! Compute effective viscosity for each layer at this quadrature point
                ! Note: efvs_qp_3d is intent(inout); old value is used to compute new value
                call compute_effective_viscosity_diva(whichefvs,            efvs_constant,     &
                                                      nz,                   stagsigma,         &
                                                      nNodesPerElement_2d,  phi_2d(:,p),       &
                                                      dphi_dx_2d(:),        dphi_dy_2d(:),     &
                                                      u(:),                 v(:),              & 
                                                      bx(:),                by(:),             &
                                                      h(:),                                    &
                                                      flwa(:,i,j),          flwafact(:,i,j),   &
                                                      efvs_qp(:,p),                            &
                                                      i, j, p)

                !WHL - Copy local efvs_qp to the global array
                efvs_qp_3d(:,:,i,j) = efvs_qp(:,:)

                ! Compute vertical average of effective viscosity
                efvs_qp_vertavg(p) = 0.d0
                do k = 1, nz-1
                   efvs_qp_vertavg(p) = efvs_qp_vertavg(p) + efvs_qp(k,p)*(sigma(k+1) - sigma(k))
                enddo

             else     ! SSA

                ! Compute vertically averaged effective viscosity at this quadrature point
                call compute_effective_viscosity(whichefvs,        whichapprox,                       &
                                                 efvs_constant,    nNodesPerElement_2d,               &
                                                 dphi_dx_2d(:),    dphi_dy_2d(:),    dphi_dz_2d(:),   &
                                                 u(:),             v(:),                              & 
                                                 flwafact_2d(i,j), efvs_qp_vertavg(p),                &
                                                 i, j, 1, p)

                ! Copy vertically averaged value to all levels
                efvs_qp(:,p) = efvs_qp_vertavg(p)

             endif    ! whichapprox

             ! Compute ice thickness at this quadrature point

             h_qp = 0.d0
             do n = 1, nNodesPerElement_2d
                h_qp = h_qp + phi_2d(n,p) * h(n)
             enddo

             ! Increment the element stiffness matrix with the contribution from each quadrature point.
             ! Note: The effective viscosity is multiplied by thickness since the equation to be solved
             !       is vertically integrated.

             call compute_element_matrix(whichapprox,     nNodesPerElement_2d,               & 
                                         wqp_2d(p),       detJ(p),                           &
                                         h_qp*efvs_qp_vertavg(p),                            &
                                         dphi_dx_2d(:),   dphi_dy_2d(:),    dphi_dz_2d(:),   &
                                         Kuu(:,:),        Kuv(:,:),                          &
                                         Kvu(:,:),        Kvv(:,:),                          &
                                         i, j, 1, p )

          enddo   ! nQuadPoints_2d

          if (whichapprox == HO_APPROX_DIVA) then

             ! Compute vertical integrals needed for the 2D solve and 3D velocity reconstruction
             call compute_integrals_diva(nz,               sigma,                &
                                         thck(i,j),        efvs_qp(:,:),  &
                                         omega_k(:,i,j),   omega(i,j),           &
                                         i, j)

          endif

          ! Compute average of effective viscosity over quad points
          ! For L1L2 and DIVA there is a different efvs in each layer.
          ! For SSA, simply write the vertical average value to each layer.

          efvs(:,i,j) = 0.d0
          do p = 1, nQuadPoints_2d
             do k = 1, nz-1
                efvs(k,i,j) = efvs(k,i,j) + efvs_qp(k,p)
             enddo
          enddo
          efvs(:,i,j) = efvs(:,i,j) / nQuadPoints_2d

          if (check_symmetry_element) then
             call check_symmetry_element_matrix(nNodesPerElement_2d,   &
                                                Kuu, Kuv, Kvu, Kvv)
          endif

          ! Sum the terms of element matrix K into the dense assembled matrix A

          call element_to_global_matrix_2d(nx,           ny,        &
                                           i,            j,         &
                                           Kuu,          Kuv,       &
                                           Kvu,          Kvv,       &
                                           Auu,          Auv,       &
                                           Avu,          Avv)

          if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest) then
             print*, ' '
             print*, 'Assembled 2D matrix, i, j =', i, j
             print*, 'k, efvs:'
             do k = 1, nz-1
                print*, k, efvs(k,i,j)
             enddo
          endif

       endif   ! active cell

    enddo      ! i
    enddo      ! j

  end subroutine assemble_stiffness_matrix_2d

!****************************************************************************

! For now, passing in i and j for debugging

  subroutine compute_integrals_diva(nz,        sigma,         & 
                                    thck,      efvs_qp,       &
                                    omega_k,   omega,   i, j)

    !----------------------------------------------------------------
    ! Compute some integrals used by the DIVA solver to relate velocities
    ! in different parts of the column:
    !
    !    F1(z) = int_b^z {[(s-z)/H] * 1/efvs * dz}
    !    F2    = int_b^s {[(s-z)/H]^2 * 1/efvs * dz}
    !          = int_b^s {F1(z)/H * dz}
    !
    ! Because efvs is highly nonlinear and appears in the denominator,
    ! it should be more accurate to compute the integral at each quadrature
    ! point and then average to the cell center, rather than average efvs 
    ! to the cell center and then integrate.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::    &
       nz                 ! number of vertical levels at which velocity is computed

    real(dp), dimension(nz), intent(in) ::    &
       sigma              ! sigma vertical coordinate

    real(dp), intent(in) ::  &
       thck               ! ice thickness (m)

    real(dp), dimension(nz-1,nQuadPoints_2d), intent(in) ::  &
       efvs_qp            ! effective viscosity (Pa yr) at each quad point in each layer

    real(dp), dimension(nz), intent(out) :: &
       omega_k            ! single integral, defined by Goldberg (2011) eq. 32

    real(dp), intent(out) :: &
       omega              ! double integral, defined by Goldberg (2011) eq. 35
                          ! Note: omega here = Goldberg's omega/H

    integer, intent(in) :: i, j   ! temporary, for debugging

    !---------------------------------------------------------
    ! Local variables
    !---------------------------------------------------------

    integer :: k, p

    real(dp), dimension(nz,nQuadPoints_2d) :: &
       omega_kp   ! omega_k in a column associated with a quad point

    real(dp) :: &
       layer_avg, dz, depth

    !WHL - debug
    real(dp), dimension(nz) :: fact_k

    omega_k(:) = 0.d0
    omega = 0.d0

    ! Compute omega_k in the vertical column at each quad point
    do p = 1, nQuadPoints_2d
       omega_kp(nz,p) = 0.d0
       do k = nz-1, 1, -1
!!          depth = 0.5d0*(sigma(k)+sigma(k+1))/thck
          depth = 0.5d0*(sigma(k)+sigma(k+1))   ! depth/thck
          dz = (sigma(k+1)-sigma(k)) * thck
          omega_kp(k,p) = omega_kp(k+1,p) + depth/efvs_qp(k,p) * dz
       enddo
    enddo

    ! Average from quad points to the cell center
    do k = 1, nz
       omega_k(k) = sum(omega_kp(k,:)) / nQuadPoints_2d
    enddo

    ! Integrate omega_k in the vertical to obtain omega
    omega = 0.d0
    do k = 1, nz-1
       layer_avg = 0.5d0*(omega_k(k) + omega_k(k+1))
!!       dz = (sigma(k+1)-sigma(k)) * thck 
       dz = sigma(k+1)-sigma(k)  ! dz/thck
       omega = omega + layer_avg * dz
    enddo
             
    if (verbose_diva .and. this_rank==rtest .and. i==itest .and. j==jtest) then
       print*, ' '
       print*, 'DIVA integrals, i, j =', i, j
       print*, 'k, integral_k:'
       do k = 1, nz
          print*, k, omega_k(k)
       enddo
       print*, 'omega =', omega
    endif

    !TODO - Test results further with this integral
    !Note - The following code computes the integral Arthern-style.
    !       The resulting omega can vary by ~50%, but code answers change little.

    do p = 1, nQuadPoints_2d
       omega_kp(nz,p) = 0.d0
       do k = 1, nz-1
          depth = 0.5d0*(sigma(k)+sigma(k+1))   ! depth/thck
          dz = (sigma(k+1)-sigma(k)) * thck
          omega_kp(k,p) = omega_kp(k+1,p) + depth**2/efvs_qp(k,p) * dz
       enddo
    enddo

    ! Average from quad points to the cell center
    do k = 1, nz
       fact_k(k) = sum(omega_kp(k,:)) / nQuadPoints_2d
    enddo
!!    omega = fact_k(1)  ! Uncomment to use Arthern value of omega
    
!    if (verbose_diva .and. this_rank==rtest .and. i==itest .and. j==jtest) then
!       print*, ' '
!       print*, 'Arthern integrals, i, j =', i, j
!       print*, 'k, fact_k:'
!       do k = 1, nz
!          print*, k, fact_k(k)
!       enddo
!       print*, 'omega =', omega
!    endif

  end subroutine compute_integrals_diva

!****************************************************************************

  subroutine compute_3d_velocity_diva(nx,               ny,                 &
                                      nz,               sigma,              &
                                      active_vertex,    diva_level_index,   &  
                                      stag_omega_k,     stag_omega,         &
                                      btractx,          btracty,            &
                                      uvel_2d,          vvel_2d,            &
                                      uvel,             vvel)
    
    !----------------------------------------------------------------
    ! Compute the 3D velocity field for the DIVA scheme,
    ! given the 2D velocity solution and the 3D effective viscosity.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz                            ! number of vertical levels at which velocity is computed

    real(dp), dimension(nz), intent(in) ::    &
       sigma              ! sigma vertical coordinate

    logical, dimension(nx-1,ny-1), intent(in) ::  &
       active_vertex      ! true for vertices of active cells

    integer, intent(in) ::   &
       diva_level_index   ! level for which the DIVA scheme finds the 2D velocity
                          ! 0 = mean, 1 = upper surface

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       stag_omega_k       ! single integral, defined by Goldberg eq. 32 (m^2/(Pa yr))
                          ! interpolated to staggered grid

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stag_omega,       &! double integral, defined by Goldberg eq. 35 (m^2/(Pa yr))
                          ! interpolated to staggered grid
                          ! Note: omega here = Goldberg's omega/H
       btractx, btracty, &! components of basal traction (Pa)
       uvel_2d, vvel_2d   ! depth-integrated mean velocity; solution of 2D velocity solve (m/yr)
                            
    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::  &
       uvel, vvel         ! 3D velocity components (m/yr)

    ! Local variables

    integer :: i, j, k

    real(dp), dimension(nx-1,ny-1) ::  &
         stag_integral         ! integral that relates bed velocity to uvel_2d/vvel_2d
                               ! = stag_omega for diva_level_index = 0
                               ! = stag_omega_k(k,:,:) for other values of diva_level_index

    ! Identify the appropriate integral for relating uvel_2d/vvel_2d to the bed velocity

    if (diva_level_index == 0) then  ! solved for mean velocity
       stag_integral(:,:) = stag_omega(:,:)
    else
       k = diva_level_index
       stag_integral(:,:) = stag_omega_k(k,:,:)
    endif

    !----------------------------------------------------------------
    ! Compute the 3D velocity field
    !----------------------------------------------------------------

    do j = 1, ny-1
       do i = 1, nx-1
          if (active_vertex(i,j)) then

             ! basal velocity (Goldberg eq. 34)
             uvel(nz,i,j) = uvel_2d(i,j) - btractx(i,j)*stag_integral(i,j) 
             vvel(nz,i,j) = vvel_2d(i,j) - btracty(i,j)*stag_integral(i,j) 
         
             ! vertical velocity profile (Goldberg eq. 32)
             do k = 1, nz-1
                uvel(k,i,j) = uvel(nz,i,j) + btractx(i,j)*stag_omega_k(k,i,j)
                vvel(k,i,j) = vvel(nz,i,j) + btracty(i,j)*stag_omega_k(k,i,j)
             enddo

          endif   ! active_vertex
       enddo      ! i
    enddo         ! j

    if (verbose_diva .and. this_rank==rtest) then
       print*, ' '
       i = itest
       j = jtest
       print*, 'Computed 3D velocities, i, j =', i, j
       print*, 'k, uvel, vvel:'
       do k = 1, nz
          print*, k, uvel(k,i,j), vvel(k,i,j)
       enddo
       print*, ' '
    endif
       
  end subroutine compute_3d_velocity_diva

!****************************************************************************

  subroutine compute_3d_velocity_L1L2(nx,               ny,              &
                                      nz,               sigma,           &
                                      dx,               dy,              &
                                      nhalo,                             &
                                      ice_mask,         land_mask,       &
                                      active_cell,      active_vertex,   &
                                      umask_dirichlet,  vmask_dirichlet, &
                                      xVertex,          yVertex,         &
                                      thck,             stagthck,        &
                                      usrf,                              &
                                      dusrf_dx,         dusrf_dy,        &
                                      flwa,             efvs,            &
                                      whichgradient_margin,              &
                                      max_slope,                         &
                                      uvel,             vvel)

    !----------------------------------------------------------------
    ! Given the basal velocity and the 3D profile of effective viscosity
    !  and horizontal-plane stresses, construct the 3D stress and velocity
    !  profiles for the L1L2 approximation.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical levels at which velocity is computed
       nhalo                         ! number of halo layers

    real(dp), intent(in) ::     &
       dx, dy             ! grid cell length and width 

    real(dp), dimension(nz), intent(in) ::    &
       sigma              ! sigma vertical coordinate

    integer, dimension(nx,ny), intent(in) ::  &
       ice_mask,        & ! = 1 for cells where ice is present (thk > thklim), else = 0
       land_mask          ! = 1 for cells where topography is above sea level

    logical, dimension(nx,ny), intent(in) ::  &
       active_cell        ! true if cell contains ice and borders a locally owned vertex

    logical, dimension(nx-1,ny-1), intent(in) ::  &
       active_vertex      ! true for vertices of active cells

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xVertex, yVertex   ! x and y coordinates of vertices

    integer, dimension(nx-1,ny-1), intent(in) ::  &
       umask_dirichlet,  &! Dirichlet mask for u velocity, = 1 for prescribed velo, else = 0
       vmask_dirichlet    ! Dirichlet mask for v velocity, = 1 for prescribed velo, else = 0

    real(dp), dimension(nx,ny), intent(in) ::  &
       thck,             &! ice thickness at cell centers (m)
       usrf               ! upper surface elevation at cell centers (m)

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagthck,       &  ! ice thickness at vertices (m)
       dusrf_dx,       &  ! upper surface elevation gradient at cell vertices (m/m)
       dusrf_dy

    real(dp), dimension(nz-1,nx,ny), intent(in) ::  &
       flwa,           &  ! temperature-based flow factor A, Pa^{-n} yr^{-1}
       efvs               ! effective viscosity, Pa yr

    integer, intent(in) ::  &
       whichgradient_margin     ! option for computing gradient at ice margin
                                ! 0 = include all neighbor cells in gradient calculation
                                ! 1 = include ice-covered and/or land cells
                                ! 2 = include ice-covered cells only

    real(dp), intent(in) ::  &
       max_slope          ! maximum slope allowed for surface gradient computations (unitless)

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::  &
       uvel, vvel         ! velocity components (m/yr)
                          ! on input, only the basal component (index nz) is known
                          ! on output, the full 3D velocity field is known

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: iVertex, jVertex  ! indices of element vertices

    real(dp), dimension(nNodesPerElement_2d) ::   &
       x, y,                    &! x and y coordinates of element vertices 
       u, v,                    &! basal velocity components at element vertices
       dphi_dx_2d, dphi_dy_2d    ! derivatives of basis functions, evaluated at cell center

    real(dp) ::   &
       detJ                      ! determinant of J (never used in calculation)

    real(dp), dimension(nx,ny) ::  &
       du_dx, du_dy,            &! basal strain rate components, evaluated at cell centers
       dv_dx, dv_dy,            &!
       work1, work2, work3       ! work arrays for computing tau_xz and tau_yz; located at cell centers

    real(dp), dimension(nz-1,nx,ny) ::   &
       tau_parallel,            &! tau_parallel, evaluated at cell centers
       efvs_integral_z_to_s      ! integral of effective viscosity from base of layer k
                                 ! to the upper surface (Pa yr m)

    ! Note: These L1L2 stresses are located at nodes.
    !       The diagnostic stresses (model%stress%tau%xz, etc.) are located at cell centers.
    real(dp), dimension(nz-1,nx-1,ny-1) ::   &
       tau_xz, tau_yz            ! vertical shear stress components at layer midpoints for each vertex
       
    real(dp), dimension(nx-1,ny-1) ::   &
       dwork1_dx, dwork1_dy,    &! derivatives of work arrays; located at vertices
       dwork2_dx, dwork2_dy,    &!
       dwork3_dx, dwork3_dy,    &!
       stagtau_parallel_sq,     &! tau_parallel^2, interpolated to staggered grid
       stagflwa                  ! flwa, interpolated to staggered grid

    real(dp) ::   &
       depth,                   &! distance from upper surface to midpoint of a given layer
       eps_parallel,            &! parallel effective strain rate, evaluated at cell centers
       tau_eff_sq,              &! square of effective stress (Pa^2)
                                 ! = tau_parallel^2 + tau_perp^2 for L1L2
       fact                      ! factor in velocity integral

    real(dp), dimension(nx-1,ny) ::  &
       dusrf_dx_edge             ! x gradient of upper surface elevation at cell edges (m/m)

    real(dp), dimension(nx,ny-1) ::  &
       dusrf_dy_edge             ! y gradient of upper surface elevation at cell edges (m/m)

    integer :: i, j, k, n

    !-----------------------------------------------------------------------------------------------
    !WHL: I tried two ways to compute the 3D velocity, given tau_perp, tau_xz and tau_yz in each layer:
    ! (1) Compute velocity at vertices using     
    !          u(z) = u_b + 2 * integral_b_to_z [A*tau_eff^(n-1)*tau_xz dz]
    !          v(z) = v_b + 2 * integral_b_to_z [A*tau_eff^(n-1)*tau_yz dz]
    ! (2) Compute velocity at edges using 
    !          uedge(z) =  (vintfact(i,j) + vintfact(i,j-1))/2.d0 * dsdx_edge 
    !          vedge(z) =  (vintfact(i,j) + vintfact(i-1,j))/2.d0 * dsdy_edge 
    !     where vintfact = 2*A*tau_eff^(n-1)*(rho*g*|grad(s)|
    !     Average uedge and vedge to vertices and add to u_b to get 3D uvel and vvel.
    !
    ! Method 2 resembles the methods used by Glide and by the Glissade local SIA solver.
    ! For the no-slip case, method 2 gives the same answers (within roundoff) as the local SIA solver.
    ! However, method 2 does not include the gradient of membrane stresses in the tau_xz and tau_yz terms
    !  (Perego et al. Eq. 27).  It does include tau_parallel in tau_eff.
    ! For the Halfar test, method 1 is slightly more accurate but can give rise to checkerboard noise.
    !   Checkerboard noise can be damped by using an upstream gradient for grad(s), but this
    !   reduces the accuracy for the Halfar test. (Method 2 with centered gradients is more
    !   accurate than method 1 with upstream gradients.)
    !-----------------------------------------------------------------------------------------------

    logical, parameter :: edge_velocity = .false.  ! if false, use method 1 as discussed above 
                                                   ! if true, use method 2

    real(dp), dimension(nx,ny) ::   &
       uedge, vedge        ! velocity components at edges of a layer, relative to bed (m/yr)
                           ! u on E edge, v on N edge (C grid)

    real(dp), dimension(nz,nx-1,ny-1) ::   &
       vintfact            ! vertical integration factor at vertices

    ! Initialize
    efvs_integral_z_to_s(:,:,:) = 0.d0
    tau_parallel(:,:,:) = 0.d0
    du_dx(:,:) = 0.d0
    du_dy(:,:) = 0.d0
    dv_dx(:,:) = 0.d0
    dv_dy(:,:) = 0.d0

    ! Compute viscosity integral and strain rates in elements.
    ! Loop over all cells that border locally owned vertices.

    do j = 1+nhalo, ny-nhalo+1
       do i = 1+nhalo, nx-nhalo+1
       
          if (active_cell(i,j)) then

             ! Load x and y coordinates and basal velocity at cell vertices

             do n = 1, nNodesPerElement_2d

                ! Determine (i,j) for this vertex
                ! The reason for the '3' is that node 3, in the NE corner of the grid cell, has index (i,j).
                ! Indices for other nodes are computed relative to this vertex.
                iVertex = i + ishift(3,n)
                jVertex = j + jshift(3,n)

                x(n) = xVertex(iVertex,jVertex)
                y(n) = yVertex(iVertex,jVertex)

                u(n) = uvel(nz,iVertex,jVertex)   ! basal velocity
                v(n) = vvel(nz,iVertex,jVertex)

             enddo

             ! Compute dphi_dx and dphi_dy at cell center

             call get_basis_function_derivatives_2d(x(:),               y(:),               &
                                                    dphi_dxr_2d_ctr(:), dphi_dyr_2d_ctr(:), &
                                                    dphi_dx_2d(:),      dphi_dy_2d(:),      &
                                                    detJ, i, j, 1)

             ! Compute basal strain rate components at cell center
             
             do n = 1, nNodesPerElement_2d
                du_dx(i,j) = du_dx(i,j) + dphi_dx_2d(n)*u(n)
                du_dy(i,j) = du_dy(i,j) + dphi_dy_2d(n)*u(n)
                
                dv_dx(i,j) = dv_dx(i,j) + dphi_dx_2d(n)*v(n)
                dv_dy(i,j) = dv_dy(i,j) + dphi_dy_2d(n)*v(n)
             enddo

             ! Compute effective strain rate (squared) at cell centers
             ! See Perego et al. eq. 17: 
             !     eps_parallel^2 = eps_xx^2 + eps_yy^2 + eps_xx*eps_yy + eps_xy^2

             eps_parallel = sqrt(du_dx(i,j)**2 + dv_dy(i,j)**2 + du_dx(i,j)*dv_dy(i,j)  &
                                 + 0.25d0*(dv_dx(i,j) + du_dy(i,j))**2)

             ! For each layer k, compute tau_parallel at cell centers
             do k = 1, nz-1
                tau_parallel(k,i,j) = 2.d0 * efvs(k,i,j) * eps_parallel
             enddo

             ! For each layer k, compute the integral of the effective viscosity from
             ! the base of layer k to the upper surface.
 
             efvs_integral_z_to_s(1,i,j) = efvs(1,i,j) * (sigma(2) - sigma(1))*thck(i,j)

             do k = 2, nz-1
                efvs_integral_z_to_s(k,i,j) = efvs_integral_z_to_s(k-1,i,j)  &
                                            + efvs(k,i,j) * (sigma(k+1) - sigma(k))*thck(i,j)
             enddo   ! k

          endif   ! active_cell

       enddo      ! i
    enddo         ! j

    !--------------------------------------------------------------------------------
    ! For each active vertex, compute the vertical shear stresses tau_xz and tau_yz
    ! in each layer of the column.
    !
    ! These stresses are given by (PGB eq. 27)
    !
    !   tau_xz(z) = -rhoi*grav*ds_dx*(s-z) + 2*d/dx[efvs_int(z) * (2*du_dx + dv_dy)]
    !                                      + 2*d/dy[efvs_int(z) *   (du_dy + dv_dx)] 
    !
    !   tau_yz(z) = -rhoi*grav*ds_dy*(s-z) + 2*d/dx[efvs_int(z) *   (du_dy + dv_dx)] 
    !                                      + 2*d/dy[efvs_int(z) * (2*dv_dy + du_dx)]
    !
    ! where efvs_int is the integral of efvs from z to s computed above;
    ! the strain rate components of basal velocity are also as computed above.
    !
    ! There is not a clean way to compute these stresses using finite-element techniques,
    ! because strain rates are discontinuous at cell edges and vertices.  Instead, we use
    ! a standard centered finite difference method to evaluate d/dx and d/dy of the
    ! bracketed terms.
    !--------------------------------------------------------------------------------

    tau_xz(:,:,:) = 0.d0
    tau_yz(:,:,:) = 0.d0

    do k = 1, nz-1   ! loop over layers

       ! Evaluate centered finite differences of bracketed terms above.
       ! We need dwork1_dx, dwork2_dx, dwork2_dy and dwork3_dx.
       ! The calls to glissade_centered_gradient compute a couple of extraneous derivatives,
       !  but these calls are simpler than inlining the gradient code.
       ! Setting gradient_margin_in = HO_GRADIENT_MARGIN_ICE_ONLY uses only ice-covered cells to
       !  compute the gradient.  This is the appropriate flag for these
       !  calls, because efvs and strain rates have no meaning in ice-free cells.

       work1(:,:) = efvs_integral_z_to_s(k,:,:) * (2.d0*du_dx(:,:) + dv_dy(:,:)) 
       work2(:,:) = efvs_integral_z_to_s(k,:,:) *      (du_dy(:,:) + dv_dx(:,:))
       work3(:,:) = efvs_integral_z_to_s(k,:,:) * (2.d0*dv_dy(:,:) + du_dx(:,:)) 

       call glissade_centered_gradient(nx,               ny,         &
                                       dx,               dy,         &
                                       work1,                        &
                                       dwork1_dx,        dwork1_dy,  &
                                       ice_mask,                     &
                                       gradient_margin_in = HO_GRADIENT_MARGIN_ICE_ONLY)

       call glissade_centered_gradient(nx,               ny,         &
                                       dx,               dy,         &
                                       work2,                        &
                                       dwork2_dx,        dwork2_dy,  &
                                       ice_mask,                     &
                                       gradient_margin_in = HO_GRADIENT_MARGIN_ICE_ONLY)
 
       call glissade_centered_gradient(nx,               ny,         &
                                       dx,               dy,         &
                                       work3,                        &
                                       dwork3_dx,        dwork3_dy,  &
                                       ice_mask,                     &
                                       gradient_margin_in = HO_GRADIENT_MARGIN_ICE_ONLY)

       ! Loop over locally owned active vertices, evaluating tau_xz and tau_yz for this layer
       do j = 1+nhalo, ny-nhalo
          do i = 1+nhalo, nx-nhalo
             if (active_vertex(i,j)) then
                depth = 0.5d0*(sigma(k) + sigma(k+1)) * stagthck(i,j)   ! depth at layer midpoint
                tau_xz(k,i,j) = -rhoi*grav*depth*dusrf_dx(i,j)   &
                               + 2.d0*dwork1_dx(i,j) + dwork2_dy(i,j)
                tau_yz(k,i,j) = -rhoi*grav*depth*dusrf_dy(i,j)   &
                               + dwork2_dx(i,j) + 2.d0*dwork3_dy(i,j)
             endif
          enddo   ! i
       enddo      ! j

    enddo         ! k
      
    if ((verbose_L1L2 .or. verbose_tau) .and. this_rank==rtest) then 
       i = itest
       j = jtest
       print*, ' '
       print*, 'L1L2: k, -rho*g*(s-z)*ds/dx, -rho*g*(s-z)*ds/dy:'
       do k = 1, nz-1
          depth = 0.5d0*(sigma(k) + sigma(k+1)) * stagthck(i,j)
          print*, k, -rhoi*grav*depth*dusrf_dx(i,j), -rhoi*grav*depth*dusrf_dy(i,j)
       enddo
       print*, ' '
       print*, 'L1L2: k, tau_xz, tau_yz, tau_parallel:'
       do k = 1, nz-1
          print*, k, tau_xz(k,i,j), tau_yz(k,i,j), tau_parallel(k,i,j)
       enddo
    endif

    !--------------------------------------------------------------------------------
    ! Given the vertical shear stresses tau_xz and tau_yz for each layer k,
    !  compute the velocity components at each level.
    !
    ! These are given by (PGB eq. 30)
    ! 
    !    u(z) = u_b + 2 * integral_b_to_z [A*tau_eff^(n-1)*tau_xz dz]
    !    v(z) = v_b + 2 * integral_b_to_z [A*tau_eff^(n-1)*tau_yz dz]
    ! 
    ! where tau_eff^2 = tau_parallel^2 + tau_perp^2
    !
    !    tau_parallel^2 = (2 * efvs * eps_parallel)^2
    !    tau_perp ^2 = tau_xz^2 + tau_yz^2
    !
    ! See comments above about method 2, with edge_velocity = .true. 
    !--------------------------------------------------------------------------------

    ! initialize uvel = vvel = 0 except at bed
       
    uvel(1:nz-1,:,:) = 0.d0
    vvel(1:nz-1,:,:) = 0.d0
    vintfact(:,:,:) = 0.d0

    ! Compute surface elevation gradient on cell edges.
    ! Setting gradient_margin_in = 0 takes the gradient over both neighboring cells,
    !  including ice-free cells.
    ! Setting gradient_margin_in = 1 computes a gradient if both neighbor cells are
    !  either ice-covered cells or land cells; else gradient = 0.
    ! Setting gradient_margin_in = 2 computes a gradient only if both neighbor cells
    !  are ice-covered.
    ! At a land margin, either 0 or 1 is appropriate, but 2 is inaccurate.
    ! At a shelf margin, either 1 or 2 is appropriate, but 0 is inaccurate.
    ! So HO_GRADIENT_MARGIN_ICE_LAND = 1 is the safest value.

    if (edge_velocity) then

       uedge(:,:) = 0.d0
       vedge(:,:) = 0.d0

       call glissade_gradient_at_edges(nx,               ny,                &
                                       dx,               dy,                &
                                       usrf,                                &
                                       dusrf_dx_edge,    dusrf_dy_edge,     &
                                       gradient_margin_in = whichgradient_margin, &
                                       ice_mask = ice_mask,                 &
                                       land_mask = land_mask,               &
                                       max_slope = max_slope)
    endif

    if (verbose_L1L2 .and. this_rank==rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'i, j =', itest, jtest
       print*, 'k, uvel, vvel:'
    endif

    do k = nz-1, 1, -1   ! loop over velocity levels above the bed
       
       ! Average tau_parallel and flwa to vertices
       ! With stagger_margin_in = 1, only cells with ice are included in the average.

       call glissade_stagger(nx,                   ny,                         &
                             tau_parallel(k,:,:),  stagtau_parallel_sq(:,:),   &
                             ice_mask,             stagger_margin_in = 1)
       stagtau_parallel_sq(:,:) = stagtau_parallel_sq(:,:)**2

       call glissade_stagger(nx,          ny,              &
                             flwa(k,:,:), stagflwa(:,:),   &
                             ice_mask,    stagger_margin_in = 1)
       
       if (edge_velocity) then  ! compute velocity at edges and interpolate to vertices
                                ! (method 2)

          ! Compute vertical integration factor at each active vertex
          ! This is int_b_to_z{-2 * A * tau^2 * rho*g*(s-z) * dz},
          !  similar to the factor computed in Glide and glissade_velo_sia..
          ! Note: tau_xz ~ rho*g*(s-z)*ds_dx; ds_dx term is computed on edges below

          do j = 1, ny-1
          do i = 1, nx-1
             if (active_vertex(i,j)) then

                tau_eff_sq = stagtau_parallel_sq(i,j)   &
                           + tau_xz(k,i,j)**2 + tau_yz(k,i,j)**2

                depth = 0.5d0*(sigma(k) + sigma(k+1)) * stagthck(i,j)

                vintfact(k,i,j) = vintfact(k+1,i,j)     &
                     - 2.d0 * stagflwa(i,j) * tau_eff_sq * rhoi*grav*depth  &
                                  * (sigma(k+1) - sigma(k))*stagthck(i,j)

             endif
          enddo
          enddo

          ! Need to have vintfact at halo nodes to compute uvel/vvel at locally owned nodes  
          call staggered_parallel_halo(vintfact(k,:,:))

          ! loop over cells, skipping outer halo rows

          ! u at east edges
          do j = 2, ny-1
          do i = 1, nx-1
             if (active_vertex(i,j) .and. active_vertex(i,j-1)) then
                uedge(i,j) = (vintfact(k,i,j) + vintfact(k,i,j-1))/2.d0 * dusrf_dx_edge(i,j)
             endif
          enddo
          enddo

          ! v at north edges
          do j = 1, ny-1
          do i = 2, nx-1
             if (active_vertex(i,j) .and. active_vertex(i-1,j)) then
                vedge(i,j) = (vintfact(k,i,j) + vintfact(k,i-1,j))/2.d0 * dusrf_dy_edge(i,j)
             endif
          enddo
          enddo

          ! Average edge velocities to vertices and add to ubas                                                                                                   
          ! Do this for locally owned vertices only
          ! (Halo update is done at a higher level after returning)
          ! Note: Currently do not support Dirichlet BC with depth-varying velocity
          
          do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo

             if (umask_dirichlet(i,j) == 1) then
                uvel(k,i,j) = uvel(nz,i,j)
             else
                uvel(k,i,j) = uvel(nz,i,j) + (uedge(i,j) + uedge(i,j+1)) / 2.d0
             endif

             if (vmask_dirichlet(i,j) == 1) then
                vvel(k,i,j) = vvel(nz,i,j)
             else
                vvel(k,i,j) = vvel(nz,i,j) + (vedge(i,j) + vedge(i+1,j)) / 2.d0
             endif

             if (verbose_L1L2 .and. this_rank==rtest .and. i==itest .and. j==jtest) then
                print*, k, uvel(k,i,j), vvel(k,i,j)
             endif

          enddo
          enddo

       else   ! compute velocity at vertices (method 1)

          ! loop over locally owned active vertices
          do j = 1+nhalo, ny-nhalo  
          do i = 1+nhalo, nx-nhalo

             if (active_vertex(i,j)) then

                tau_eff_sq = stagtau_parallel_sq(i,j)   &
                           + tau_xz(k,i,j)**2 + tau_yz(k,i,j)**2

                ! Note: This formula is correct for any value of Glen's n, but currently efvs is computed
                !       only for gn = 3 (in which case (n-1)/2 = 1).
                fact = 2.d0 * stagflwa(i,j) * tau_eff_sq**((gn-1.d0)/2.d0) * (sigma(k+1) - sigma(k))*stagthck(i,j)

                ! reset velocity to prescribed basal value if Dirichlet condition applies
                ! else compute velocity at this level 
                if (umask_dirichlet(i,j) == 1) then
                   uvel(k,i,j) = uvel(nz,i,j)
                else
                   uvel(k,i,j) = uvel(k+1,i,j) + fact * tau_xz(k,i,j)
                endif

                if (vmask_dirichlet(i,j) == 1) then
                   vvel(k,i,j) = vvel(nz,i,j)
                else
                   vvel(k,i,j) = vvel(k+1,i,j) + fact * tau_yz(k,i,j)
                endif

                if (verbose_L1L2 .and. this_rank==rtest .and. i==itest .and. j==jtest) then
                   print*, k, uvel(k,i,j), vvel(k,i,j)
                endif

             endif

          enddo   ! i
          enddo   ! j

       endif      ! edge_velocity

    enddo         ! k

  end subroutine compute_3d_velocity_L1L2

!****************************************************************************

  subroutine get_basis_function_derivatives_3d(xNode,       yNode,       zNode,       &
                                               dphi_dxr_3d, dphi_dyr_3d, dphi_dzr_3d, &
                                               dphi_dx_3d,  dphi_dy_3d,  dphi_dz_3d,  &
                                               detJ,        i, j, k, p)

    !------------------------------------------------------------------
    ! Evaluate the x, y and z derivatives of the element basis functions
    ! at a particular quadrature point.
    !
    ! Also determine the Jacobian of the transformation between the
    ! reference element and the true element.
    ! 
    ! This subroutine should work for any 3D element with any number of nodes.
    !------------------------------------------------------------------
 
    real(dp), dimension(nNodesPerElement_3d), intent(in) :: &
       xNode, yNode, zNode,          &! nodal coordinates
       dphi_dxr_3d, dphi_dyr_3d, dphi_dzr_3d   ! derivatives of basis functions at quad pt
                                               !  wrt x, y and z in reference element

    real(dp), dimension(nNodesPerElement_3d), intent(out) :: &
       dphi_dx_3d, dphi_dy_3d, dphi_dz_3d      ! derivatives of basis functions at quad pt
                                               !  wrt x, y and z in true Cartesian coordinates  

    real(dp), intent(out) :: &
         detJ      ! determinant of Jacobian matrix

    real(dp), dimension(3,3) ::  &
         Jac,      &! Jacobian matrix
         Jinv,     &! inverse Jacobian matrix
         cofactor   ! matrix of cofactors

    integer, intent(in) :: i, j, k, p   ! indices passed in for debugging

    integer :: n, row, col

    logical, parameter :: Jac_bug_check = .false.   ! set to true for debugging
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

    if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, ' '
       print*, 'In get_basis_function_derivatives_3d: i, j, k, p =', i, j, k, p
    endif

    Jac(:,:) = 0.d0

    do n = 1, nNodesPerElement_3d
       Jac(1,1) = Jac(1,1) + dphi_dxr_3d(n) * xNode(n)
       Jac(1,2) = Jac(1,2) + dphi_dxr_3d(n) * yNode(n)
       Jac(1,3) = Jac(1,3) + dphi_dxr_3d(n) * zNode(n)
       Jac(2,1) = Jac(2,1) + dphi_dyr_3d(n) * xNode(n)
       Jac(2,2) = Jac(2,2) + dphi_dyr_3d(n) * yNode(n)
       Jac(2,3) = Jac(2,3) + dphi_dyr_3d(n) * zNode(n)
       Jac(3,1) = Jac(3,1) + dphi_dzr_3d(n) * xNode(n)
       Jac(3,2) = Jac(3,2) + dphi_dzr_3d(n) * yNode(n)
       Jac(3,3) = Jac(3,3) + dphi_dzr_3d(n) * zNode(n)
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

    if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
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
       print*, 'stopping, det J = 0'
       print*, 'i, j, k, p:', i, j, k, p
       print*, 'Jacobian matrix:'
       print*, Jac(1,:)
       print*, Jac(2,:)
       print*, Jac(3,:) 
       call write_log('Jacobian matrix is singular', GM_FATAL)
    endif

    if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
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

    ! Optional bug check: Verify that J * Jinv = I

    if (Jac_bug_check) then
       prod = matmul(Jac,Jinv)
       do col = 1, 3
          do row = 1, 3
             if (abs(prod(row,col) - identity3(row,col)) > 1.d-11) then
                print*, 'stopping, Jac * Jinv /= identity'
                print*, 'i, j, k, p:', i, j, k, p
                print*, 'Jac*Jinv:'
                print*, prod(1,:)
                print*, prod(2,:)
                print*, prod(3,:)
                call write_log('Jacobian matrix was not correctly inverted', GM_FATAL)
             endif
          enddo
       enddo
    endif  ! Jac_bug_check

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

    dphi_dx_3d(:) = 0.d0
    dphi_dy_3d(:) = 0.d0
    dphi_dz_3d(:) = 0.d0

    do n = 1, nNodesPerElement_3d
       dphi_dx_3d(n) = Jinv(1,1)*dphi_dxr_3d(n)  &
                     + Jinv(1,2)*dphi_dyr_3d(n)  &
                     + Jinv(1,3)*dphi_dzr_3d(n)
       dphi_dy_3d(n) = Jinv(2,1)*dphi_dxr_3d(n)  &
                     + Jinv(2,2)*dphi_dyr_3d(n)  &
                     + Jinv(2,3)*dphi_dzr_3d(n)
       dphi_dz_3d(n) = Jinv(3,1)*dphi_dxr_3d(n)  &
                     + Jinv(3,2)*dphi_dyr_3d(n)  &
                     + Jinv(3,3)*dphi_dzr_3d(n)
    enddo

    if (Jac_bug_check) then

       ! Check that the sum of dphi_dx, etc. is close to zero  

       if (abs( sum(dphi_dx_3d)/maxval(dphi_dx_3d) ) > 1.d-11) then
          print*, 'stopping, sum over basis functions of dphi_dx > 0'
          print*, 'dphi_dx_3d =', dphi_dx_3d(:)
          print*, 'sum =', sum(dphi_dx_3d)
          print*, 'i, j, k, p =', i, j, k, p
          call write_log('Sum over basis functions of dphi_dx /= 0', GM_FATAL)
       endif

       if (abs( sum(dphi_dy_3d)/maxval(dphi_dy_3d) ) > 1.d-11) then
          print*, 'stopping, sum over basis functions of dphi_dy > 0'
          print*, 'dphi_dy_3d =', dphi_dy_3d(:)
          print*, 'sum =', sum(dphi_dy_3d)
          print*, 'i, j, k, p =', i, j, k, p
          call write_log('Sum over basis functions of dphi_dy /= 0', GM_FATAL)
       endif

       if (abs( sum(dphi_dz_3d)/maxval(dphi_dz_3d) ) > 1.d-11) then
          print*, 'stopping, sum over basis functions of dphi_dz > 0'
          print*, 'dphi_dz_3d =', dphi_dz_3d(:)
          print*, 'sum =', sum(dphi_dz_3d)
          print*, 'i, j, k, p =', i, j, k, p
          call write_log('Sum over basis functions of dphi_dz /= 0', GM_FATAL)
       endif

    endif  ! Jac_bug_check

  end subroutine get_basis_function_derivatives_3d

!****************************************************************************

  subroutine get_basis_function_derivatives_2d(xNode,       yNode,         &
                                               dphi_dxr_2d, dphi_dyr_2d,   &
                                               dphi_dx_2d,  dphi_dy_2d,    &
                                               detJ, i, j, p)

    !------------------------------------------------------------------
    ! Evaluate the x and y derivatives of 2D element basis functions
    ! at a particular quadrature point.
    !
    ! Also determine the Jacobian of the transformation between the
    ! reference element and the true element.
    ! 
    ! This subroutine should work for any 2D element with any number of nodes.
    !------------------------------------------------------------------

    real(dp), dimension(nNodesPerElement_2d), intent(in) :: &
       xNode, yNode,                   &! nodal coordinates
       dphi_dxr_2d, dphi_dyr_2d         ! derivatives of basis functions at quad pt
                                        !  wrt x and y in reference element

    real(dp), dimension(nNodesPerElement_2d), intent(out) :: &
       dphi_dx_2d, dphi_dy_2d           ! derivatives of basis functions at quad pt
                                        !  wrt x and y in true Cartesian coordinates  

    real(dp), intent(out) :: &
                detJ      ! determinant of Jacobian matrix

    real(dp), dimension(2,2) ::  &
                Jac,      &! Jacobian matrix
                Jinv       ! inverse Jacobian matrix

    integer, intent(in) :: i, j, p

    integer :: n, row, col

    logical, parameter :: Jac_bug_check = .false.   ! set to true for debugging
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

    if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest) then
       print*, ' '
       print*, 'In get_basis_function_derivatives_2d: i, j, p =', i, j, p
    endif

    do n = 1, nNodesPerElement_2d
       if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest) then
          print*, ' '
          print*, 'n, x, y:', n, xNode(n), yNode(n)
          print*, 'dphi_dxr_2d, dphi_dyr_2d:', dphi_dxr_2d(n), dphi_dyr_2d(n)
       endif
       Jac(1,1) = Jac(1,1) + dphi_dxr_2d(n) * xNode(n)
       Jac(1,2) = Jac(1,2) + dphi_dxr_2d(n) * yNode(n)
       Jac(2,1) = Jac(2,1) + dphi_dyr_2d(n) * xNode(n)
       Jac(2,2) = Jac(2,2) + dphi_dyr_2d(n) * yNode(n)
    enddo

    !------------------------------------------------------------------
    ! Compute the determinant and inverse of J
    !------------------------------------------------------------------

    detJ = Jac(1,1)*Jac(2,2) - Jac(1,2)*Jac(2,1)

    if (abs(detJ) > 0.d0) then
       Jinv(1,1) =  Jac(2,2)/detJ
       Jinv(1,2) = -Jac(1,2)/detJ
       Jinv(2,1) = -Jac(2,1)/detJ
       Jinv(2,2) =  Jac(1,1)/detJ
    else
       print*, 'stopping, det J = 0'
       print*, 'i, j, p:', i, j, p
       print*, 'Jacobian matrix:'
       print*, Jac(1,:)
       print*, Jac(2,:)
       call write_log('Jacobian matrix is singular', GM_FATAL)
    endif

    if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest) then
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

    ! Optional bug check - Verify that J * Jinv = I

    if (Jac_bug_check) then
       prod = matmul(Jac,Jinv)
       do col = 1, 2
          do row = 1, 2
             if (abs(prod(row,col) - identity3(row,col)) > 1.d-12) then
                print*, 'stopping, Jac * Jinv /= identity'
                print*, 'i, j, p:', i, j, p
                print*, 'Jac*Jinv:'
                print*, prod(1,:)
                print*, prod(2,:)
                call write_log('Jacobian matrix was not correctly inverted', GM_FATAL)
             endif
          enddo
       enddo
    endif

    !------------------------------------------------------------------
    ! Compute the contribution of this quadrature point to dphi/dx and dphi/dy
    ! for each basis function.
    !
    !   | dphi_n/dx |          | dphi_n/dxr |
    !   |           | = Jinv * |            |
    !   | dphi_n/dy |          | dphi_n/dyr |
    !
    !------------------------------------------------------------------

    dphi_dx_2d(:) = 0.d0
    dphi_dy_2d(:) = 0.d0

    do n = 1, nNodesPerElement_2d
       dphi_dx_2d(n) = dphi_dx_2d(n) + Jinv(1,1)*dphi_dxr_2d(n)  &
                                     + Jinv(1,2)*dphi_dyr_2d(n)
       dphi_dy_2d(n) = dphi_dy_2d(n) + Jinv(2,1)*dphi_dxr_2d(n)  &
                                     + Jinv(2,2)*dphi_dyr_2d(n)
    enddo

    if (Jac_bug_check) then

       ! Check that the sum of dphi_dx, etc. is close to zero  
       if (abs( sum(dphi_dx_2d)/maxval(dphi_dx_2d) ) > 1.d-11) then
          print*, 'stopping, sum over basis functions of dphi_dx > 0'
          print*, 'dphi_dx_2d =', dphi_dx_2d(:)
          print*, 'i, j, p =', i, j, p
          call write_log('Sum over basis functions of dphi_dx /= 0', GM_FATAL)
       endif

       if (abs( sum(dphi_dy_2d)/maxval(dphi_dy_2d) ) > 1.d-11) then
          print*, 'stopping, sum over basis functions of dphi_dy > 0'
          print*, 'dphi_dy =', dphi_dy_2d(:)
          print*, 'i, j, p =', i, j, p
          call write_log('Sum over basis functions of dphi_dy /= 0', GM_FATAL)
       endif

    endif

  end subroutine get_basis_function_derivatives_2d

!****************************************************************************

  subroutine compute_basal_friction_heatflx(nx,            ny,            &
                                            nhalo,         active_cell,   &
                                            xVertex,       yVertex,       &
                                            uvel,          vvel,          &
                                            beta,          whichassemble_bfric,  &
                                            bfricflx)

    !----------------------------------------------------------------
    ! Compute the heat flux due to basal friction, given the 2D basal
    !  velocity and beta fields.
    !
    ! Assume a sliding law of the form:
    !   tau_x = -beta*u
    !   tau_y = -beta*v
    ! where beta and (u,v) are defined at vertices.
    ! 
    ! The frictional heat flux (W/m^2) is given by q_b = tau_b * u_b,
    ! where tau_b and u_b are the magnitudes of the basal stress
    ! and velocity (e.g., Cuffey & Paterson, p. 418).
    !
    ! Note: There is a choice of two methods for this calculation:
    !       (0) a finite-element method, summing over beta*(u^2 + v^2) at quadrature points
    !       (1) a simple method, computing beta*(u^2 + v^2) at vertices
    !       Method (0) should formally be more accurate, at least where the flow is smooth.
    !       However, it can lead to inaccurate and hugely excessive frictional fluxes where
    !        the flow transitions steeply from high beta/low velo to low beta/high velo
    !        (e.g., at the edge of fjords). In this case there are QPs with relatively
    !        high velocity combined with large beta. 
    !       To choose method (1), set which_ho_assemble_bfric = 1 in the config file.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nhalo                         ! number of halo layers

    logical, dimension(nx,ny), intent(in) ::  &
       active_cell            ! true if cell contains ice and borders a locally owned vertex

    real(dp), dimension(nx-1,ny-1), intent(in) :: &
       xVertex, yVertex       ! x and y coordinates of each vertex (m)

    real(dp), dimension(nx-1,ny-1), intent(in) :: &
       uvel, vvel,          & ! basal velocity components at each vertex (m/yr)
       beta                   ! basal traction parameter (Pa/(m/yr))
                              ! typically = beta_internal (beta weighted by f_ground)

    integer, intent(in) ::  &
       whichassemble_bfric    ! = 0 for standard finite element computation of basal friction
                              ! = 1 for computation that uses only the local value of the basal friction at each vertex

    real(dp), dimension(nx,ny), intent(out) :: &
       bfricflx               ! basal heat flux from friction (W/m^2)

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j, n, p
    integer :: iVertex, jVertex

    real(dp), dimension(nNodesPerElement_2d) ::   &
       x, y,            & ! spatial coordinates of nodes
       u, v,            & ! velocity components at nodes
       b                  ! beta at nodes

    real(dp) ::         &
       u_qp, v_qp,      & ! u and v at quadrature points
       beta_qp,         & ! beta at quadrature points
       sum_wqp            ! sum of weighting factors

    logical, parameter :: bfricflx_finite_element = .false.  ! if true, do a finite-element summation
                                                             ! if false, take beta*(u^2 + v^2) at active vertices
                                                             ! (see comments above)
    ! initialize
    bfricflx(:,:) = 0.d0

    if (whichassemble_bfric == HO_ASSEMBLE_BFRIC_STANDARD) then

       ! do finite-element calculation (can be inaccurate at sharp transitions in beta and velocity)

       ! Loop over local cells
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
       
             if (active_cell(i,j)) then   ! ice is present

                ! Load x and y coordinates, basal velocity, and beta at cell vertices

                do n = 1, nNodesPerElement_2d

                   ! Determine (i,j) for this vertex
                   ! The reason for the '3' is that node 3, in the NE corner of the grid cell, has index (i,j).
                   ! Indices for other nodes are computed relative to this vertex.
                   iVertex = i + ishift(3,n)
                   jVertex = j + jshift(3,n)
                   
                   x(n) = xVertex(iVertex,jVertex)
                   y(n) = yVertex(iVertex,jVertex)
                   u(n) = uvel(iVertex,jVertex)
                   v(n) = vvel(iVertex,jVertex)
                   b(n) = beta(iVertex,jVertex)
                   
                enddo

                sum_wqp = 0.d0

                ! loop over quadrature points
                do p = 1, nQuadPoints_2d
                   
                   ! Evaluate u, v and beta at this quadrature point
                   
                   u_qp = 0.d0
                   v_qp = 0.d0
                   beta_qp = 0.d0
                   do n = 1, nNodesPerElement_2d
                      u_qp = u_qp + phi_2d(n,p) * u(n)
                      v_qp = v_qp + phi_2d(n,p) * v(n)
                      beta_qp = beta_qp + phi_2d(n,p) * b(n)
                   enddo
                   
                   ! Increment basal frictional heating
                   
                   bfricflx(i,j) = bfricflx(i,j) + wqp_2d(p) * beta_qp * (u_qp**2 + v_qp**2)
                   sum_wqp = sum_wqp + wqp_2d(p)
                   
                   if (verbose_bfric .and. this_rank==rtest .and. i==itest .and. j==jtest) then
                      print*, ' '
                      print*, 'Increment basal friction heating, i, j, p =', i, j, p
                      print*, 'u, v, beta_qp =', u_qp, v_qp, beta_qp
                      print*, 'local increment =', beta_qp * (u_qp**2 + v_qp**2) / scyr
                   endif
                   
                enddo   ! nQuadPoints_2d
                
                ! Scale the result:
                ! Divide by sum_wqp to get average of beta*(u^2 + v^2) over cell
                ! Divide by scyr to convert Pa m/yr to Pa m/s = W/m^2
                
                bfricflx(i,j) = bfricflx(i,j) / (sum_wqp * scyr)
                
                if (verbose_bfric .and. this_rank==rtest .and. i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'i, j, bfricflx:', i, j, bfricflx(i,j)
                   print*, 'beta, uvel, vvel:', beta(i,j), uvel(i,j), vvel(i,j)
                endif
                
             endif      ! active_cell
             
          enddo         ! i
       enddo            ! j
       
    else   ! whichassemble_bfric = HO_ASSEMBLE_BFRIC_LOCAL; local calculation at active vertices

       ! Loop over local vertices
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
       
             if (active_cell(i,j)) then   ! ice is present

                bfricflx(i,j) = beta(i,j) * (uvel(i,j)**2 + vvel(i,j)**2)
                bfricflx(i,j) = bfricflx(i,j) / scyr   ! convert Pa m/yr to Pa m/s = W/m^2

                if (verbose_bfric .and. this_rank==rtest .and. i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'i, j, bfricflx:', i, j, bfricflx(i,j)
                   print*, 'beta, uvel, vvel:', beta(i,j), uvel(i,j), vvel(i,j)
                endif
                
             endif      ! active_cell
             
          enddo         ! i
       enddo            ! j

    endif  ! whichassemble_bfric

    ! halo update
    call parallel_halo(bfricflx)

  end subroutine compute_basal_friction_heatflx

!****************************************************************************

  subroutine compute_internal_stress (nx,            ny,            &
                                      nz,            sigma,         &
                                      nhalo,         active_cell,   &
                                      xVertex,       yVertex,       &
                                      stagusrf,      stagthck,      &
                                      flwafact,      efvs,          &
                                      whichefvs,     efvs_constant, &
                                      whichapprox,                  &
                                      uvel,          vvel,          &
                                      tau_xz,        tau_yz,        &
                                      tau_xx,        tau_yy,        &
                                      tau_xy,        tau_eff)

    !----------------------------------------------------------------
    ! Compute internal ice stresses at the center of each element,
    !  given the 3D velocity field and flow factor.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical levels at which velocity is computed
       nhalo                         ! number of halo layers

    real(dp), dimension(nz), intent(in) ::    &
       sigma              ! sigma vertical coordinate

    logical, dimension(nx,ny), intent(in) ::  &
       active_cell        ! true if cell contains ice and borders a locally owned vertex

    real(dp), dimension(nx-1,ny-1), intent(in) :: &
       xVertex, yVertex       ! x and y coordinates of each vertex (m)

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf,       &  ! upper surface elevation on staggered grid (m)
       stagthck           ! ice thickness on staggered grid (m)

    integer, intent(in) ::   &
       whichapprox,     & ! option for Stokes approximation (BP, L1L2, SSA, SIA)
       whichefvs          ! option for effective viscosity calculation 

    real(dp), intent(in) :: &
       efvs_constant      ! constant value of effective viscosity (Pa yr)

    real(dp), dimension(nz-1,nx,ny), intent(in) ::  &
       efvs,           &  ! precomputed effective viscosity
                          ! used for L1L2 only; efvs is recomputed at QPs for other approximations
       flwafact           ! temperature-based flow factor, 0.5 * A^(-1/n), Pa yr^(1/n)
                          ! used to compute the effective viscosity

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       uvel, vvel         ! velocity components at each node (m/yr)

    ! stress tensor components, co-located with efvs at the center of each element
    real(dp), dimension(nz-1,nx,ny), intent(out) ::   &
       tau_xz, tau_yz,         &! vertical components of stress tensor (Pa)
       tau_xx, tau_yy, tau_xy, &! horizontal components of stress tensor (Pa)
       tau_eff                  ! effective stress (Pa)

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    real(dp), dimension(nNodesPerElement_3d) ::  &
       dphi_dx_3d, dphi_dy_3d, dphi_dz_3d   ! derivatives of 3D nodal basis functions at a quadrature point

    real(dp) ::               &
       detJ,                  & ! determinant of Jacobian at a quad pt
                                ! not used but part of interface to get_basis_function_derivatives
       du_dx, du_dy, du_dz,   & ! strain rate components
       dv_dx, dv_dy, dv_dz,   & 
       efvs_qp                  ! effective viscosity at a quad pt (Pa yr)

    real(dp), dimension(nNodesPerElement_3d) ::   &
       x, y, z,         & ! spatial coordinates of nodes
       u, v               ! velocity components at nodes

    integer :: i, j, k, n, p
    integer :: iNode, jNode, kNode
   
    ! initialize stresses
    tau_xz (:,:,:) = 0.d0
    tau_yz (:,:,:) = 0.d0
    tau_xx (:,:,:) = 0.d0
    tau_yy (:,:,:) = 0.d0
    tau_xy (:,:,:) = 0.d0
    tau_eff(:,:,:) = 0.d0

    ! Loop over cells that border locally owned vertices

    do j = 1+nhalo, ny-nhalo+1
       do i = 1+nhalo, nx-nhalo+1
       
          if (active_cell(i,j)) then

             ! Loop over layers
             do k = 1, nz-1

                ! compute spatial coordinates and velocity for each node of this element
                do n = 1, nNodesPerElement_3d

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
                   
                enddo   ! nodes per element

                ! Loop over quadrature points
                do p = 1, nQuadPoints_3d

                   ! Compute derivative of basis functions at this quad pt
                   call get_basis_function_derivatives_3d(x(:),             y(:),             z(:),              &          
                                                          dphi_dxr_3d(:,p), dphi_dyr_3d(:,p), dphi_dzr_3d(:,p),  &
                                                          dphi_dx_3d(:),    dphi_dy_3d(:),    dphi_dz_3d(:),     &
                                                          detJ, i, j, k, p  )

                   ! Compute strain rates at this quadrature point, looping over nodes of element
                   du_dx = 0.d0
                   du_dy = 0.d0
                   du_dz = 0.d0
                   dv_dx = 0.d0
                   dv_dy = 0.d0
                   dv_dz = 0.d0

                   if (whichapprox == HO_APPROX_SIA) then

                      do n = 1, nNodesPerElement_3d
                         du_dz = du_dz + dphi_dz_3d(n)*u(n)
                         dv_dz = dv_dz + dphi_dz_3d(n)*v(n)
                      enddo

                   elseif (whichapprox == HO_APPROX_SSA) then

                      do n = 1, nNodesPerElement_3d
                         du_dx = du_dx + dphi_dx_3d(n)*u(n)
                         du_dy = du_dy + dphi_dy_3d(n)*u(n)
                         dv_dx = dv_dx + dphi_dx_3d(n)*v(n)
                         dv_dy = dv_dy + dphi_dy_3d(n)*v(n)
                      enddo

                   else    !  3D higher-order (BP or L1L2)
 
                      do n = 1, nNodesPerElement_3d
                         du_dx = du_dx + dphi_dx_3d(n)*u(n)
                         du_dy = du_dy + dphi_dy_3d(n)*u(n)
                         du_dz = du_dz + dphi_dz_3d(n)*u(n)
                         dv_dx = dv_dx + dphi_dx_3d(n)*v(n)
                         dv_dy = dv_dy + dphi_dy_3d(n)*v(n)
                         dv_dz = dv_dz + dphi_dz_3d(n)*v(n)
                      enddo

                   endif  ! whichapprox

                   if (whichapprox == HO_APPROX_L1L2) then

                      ! efvs is computed in a complicated way for L1L2.
                      ! Instead of recomputing it here for each QP, simply assume that the value at each QP
                      !  is equal to the average efvs in the element. This will give a small averaging error.

                      efvs_qp = efvs(k,i,j)

                   else  ! other approximations (SIA, SSA, BP)

                      ! Compute the effective viscosity at this quadrature point.

                      call compute_effective_viscosity(whichefvs,        whichapprox,                       &
                                                       efvs_constant,    nNodesPerElement_3d,               &
                                                       dphi_dx_3d(:),    dphi_dy_3d(:),    dphi_dz_3d(:),   &
                                                       u(:),             v(:),                              & 
                                                       flwafact(k,i,j),  efvs_qp,                           &
                                                       i, j, k, p)

                   endif

                   ! Increment stresses, adding the value at this quadrature point

                   tau_xz(k,i,j) = tau_xz(k,i,j) + efvs_qp * du_dz            ! 2 * efvs * eps_xz
                   tau_yz(k,i,j) = tau_yz(k,i,j) + efvs_qp * dv_dz            ! 2 * efvs * eps_yz
                   tau_xx(k,i,j) = tau_xx(k,i,j) + 2.d0 * efvs_qp * du_dx     ! 2 * efvs * eps_xx
                   tau_yy(k,i,j) = tau_yy(k,i,j) + 2.d0 * efvs_qp * dv_dy     ! 2 * efvs * eps_yy
                   tau_xy(k,i,j) = tau_xy(k,i,j) + efvs_qp * (dv_dx + du_dy)  ! 2 * efvs * eps_xy

                enddo     ! p

                ! Final stress tensor components, averaged over quad pts
                tau_xz(k,i,j) = tau_xz(k,i,j) / nQuadPoints_3d
                tau_yz(k,i,j) = tau_yz(k,i,j) / nQuadPoints_3d
                tau_xx(k,i,j) = tau_xx(k,i,j) / nQuadPoints_3d
                tau_yy(k,i,j) = tau_yy(k,i,j) / nQuadPoints_3d
                tau_xy(k,i,j) = tau_xy(k,i,j) / nQuadPoints_3d
                
                ! Effective stress
                tau_eff(k,i,j) = sqrt(tau_xx(k,i,j)**2 + tau_yy(k,i,j)**2             &
                                    + tau_xx(k,i,j)*tau_yy(k,i,j) + tau_xy(k,i,j)**2  &
                                    + tau_xz(k,i,j)**2 + tau_yz(k,i,j)**2)

             enddo  ! k

             if (verbose_tau .and. this_rank==rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'i, j =', i, j
                print*, 'k, tau_xz, tau_yz, tau_xx, tau_yy, tau_xy, tau_eff:'
                do k = 1, nz-1
                   print*, k, tau_xz(k,i,j), tau_yz(k,i,j), tau_xx(k,i,j), &
                              tau_yy(k,i,j), tau_xy(k,i,j), tau_eff(k,i,j)
                enddo
             endif   ! verbose_tau

          endif     ! active cell
       enddo        ! i
    enddo           ! j

  end subroutine compute_internal_stress

!****************************************************************************

  subroutine compute_effective_viscosity (whichefvs,     whichapprox,            &
                                          efvs_constant, nNodesPerElement,       &
                                          dphi_dx,       dphi_dy,    dphi_dz,    &
                                          uvel,          vvel,                   &
                                          flwafact,      efvs,                   &
                                          i, j, k, p )

    ! Compute effective viscosity at a quadrature point, based on the latest
    !  guess for the velocity field
    ! Note: Elements can be either 2D or 3D

    integer, intent(in) :: i, j, k, p

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) :: &
       whichefvs       ! method for computing effective viscosity
                       ! 0 = constant value
                       ! 1 = proportional to flow factor
                       ! 2 = nonlinear function of effective strain rate 

    integer, intent(in) :: &
       whichapprox     ! option for Stokes approximation (BP, SSA, SIA)

    real(dp), intent(in) :: &
       efvs_constant   ! constant value of effective viscosity (Pa yr)

    integer, intent(in) :: nNodesPerElement   ! number of nodes per element
                                              ! = 4 for 2D, = 8 for 3D

    real(dp), dimension(nNodesPerElement), intent(in) ::  &
       dphi_dx, dphi_dy, dphi_dz   ! derivatives of basis functions at this quadrature point
                                   ! dphi_dz = 0 for 2D SSA

    real(dp), dimension(nNodesPerElement), intent(in) ::  &
       uvel, vvel      ! current guess for velocity at each node of element (m/yr)

    real(dp), intent(in) ::  &
       flwafact        ! temperature-based flow factor for this element, 0.5 * A^{-1/n}
                       ! units: Pa yr^{1/n}

    real(dp), intent(out) ::   &
       efvs            ! effective viscosity at this quadrature point (Pa yr)
                       ! computed as 0.5 * A^{-1/n) * effstrain^{(1-n)/n)}
                       
    !----------------------------------------------------------------
    ! Local parameters
    !----------------------------------------------------------------

    !TODO - Test sensitivity of model convergence to effstrain_min
    real(dp), parameter ::   &
!!       effstrain_min = 1.d-20*scyr,     &! minimum value of effective strain rate, yr^{-1}
                                           ! GLAM uses 1.d-20 s^{-1} for minimum effective strain rate
       effstrain_min = 1.d-8,     &! minimum value of effective strain rate, yr^{-1}
                                   ! Mauro Perego suggests 1.d-8 yr^{-1}
       p_effstr  = (1.d0 - real(gn,dp))/real(gn,dp),  &! exponent (1-n)/n in effective viscosity relation
       p2_effstr = p_effstr/2                          ! exponent (1-n)/(2n) in effective viscosity relation

                                                               
    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    real(dp) ::               &
       du_dx, du_dy, du_dz,   & ! strain rate components
       dv_dx, dv_dy, dv_dz,   &
       effstrain,             & ! effective strain rate, yr^{-1}
       effstrainsq              ! square of effective strain rate
        
    integer :: n

    real(dp), parameter :: p2 = -1.d0/3.d0
  
    select case(whichefvs)

    case(HO_EFVS_CONSTANT)

       ! Steve Price recommends 10^6 to 10^7 Pa yr
       ! ISMIP-HOM Test F requires 2336041.42829 Pa yr; this is the default value set in glide_types.F90
       efvs = efvs_constant

       if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
          print*, 'Set efvs = constant (Pa yr):', efvs
       endif

    case(HO_EFVS_FLOWFACT)      ! set the effective viscosity to a multiple of the flow factor, 0.5*A^(-1/n)
                 
       ! Units: flwafact has units Pa yr^{1/n}
       !        effstrain has units yr^{-1}
       !        p_effstr = (1-n)/n 
       !                 = -2/3 for n=3
       ! Thus efvs has units Pa yr
 
       !TODO - Test HO_EFVS_FLOWFACT option and make sure the units and scales are OK

       effstrain = vel_scale/len_scale * scyr  ! typical strain rate, yr^{-1}
       efvs = flwafact * effstrain**p_effstr  

       if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
          print*, 'flwafact, effstrain (yr-1), efvs (Pa yr)=', flwafact, effstrain, efvs
       endif

    case(HO_EFVS_NONLINEAR)    ! compute effective viscosity based on effective strain rate

       ! initialize strain rates
       du_dx = 0.d0
       du_dy = 0.d0
       du_dz = 0.d0
       dv_dx = 0.d0
       dv_dy = 0.d0
       dv_dz = 0.d0
       
       ! Compute effective strain rate (squared) at this quadrature point (PGB 2012, eq. 3 and 9)
       ! Units are yr^(-1)

       if (whichapprox == HO_APPROX_SIA) then

          do n = 1, nNodesPerElement
             du_dz = du_dz + dphi_dz(n)*uvel(n)
             dv_dz = dv_dz + dphi_dz(n)*vvel(n)
          enddo

          effstrainsq = effstrain_min**2          &
                      + 0.25d0 * (du_dz**2 + dv_dz**2)

       elseif (whichapprox == HO_APPROX_SSA) then

          do n = 1, nNodesPerElement

             du_dx = du_dx + dphi_dx(n)*uvel(n)
             du_dy = du_dy + dphi_dy(n)*uvel(n)

             dv_dx = dv_dx + dphi_dx(n)*vvel(n)
             dv_dy = dv_dy + dphi_dy(n)*vvel(n)

          enddo

          effstrainsq = effstrain_min**2          &
                      + (du_dx**2 + dv_dy**2 + du_dx*dv_dy + 0.25d0*(dv_dx + du_dy)**2)

       else   ! 3D higher-order

          do n = 1, nNodesPerElement

             du_dx = du_dx + dphi_dx(n)*uvel(n)
             du_dy = du_dy + dphi_dy(n)*uvel(n)
             du_dz = du_dz + dphi_dz(n)*uvel(n)

             dv_dx = dv_dx + dphi_dx(n)*vvel(n)
             dv_dy = dv_dy + dphi_dy(n)*vvel(n)
             dv_dz = dv_dz + dphi_dz(n)*vvel(n)

          enddo

          effstrainsq = effstrain_min**2                                      &
                      + (du_dx**2 + dv_dy**2 + du_dx*dv_dy + 0.25d0*(dv_dx + du_dy)**2)  &
                      + 0.25d0*(du_dz**2 + dv_dz**2)

       endif  ! whichapprox

       ! Compute effective viscosity (PGB 2012, eq. 4)
       ! Units: flwafact has units Pa yr^{1/n}
       !        effstrain has units yr^{-1}
       !        p2_effstr = (1-n)/(2n) 
       !                  = -1/3 for n=3
       ! Thus efvs has units Pa yr
 
       efvs = flwafact * effstrainsq**p2_effstr

       if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest .and. p==ptest) then
          print*, ' '
          print*, 'i, j, k, p =', i, j, k, p
          print*, 'flwafact, effstrain (yr-1), efvs(Pa yr) =', flwafact, effstrain, efvs
       endif
 
   end select

  end subroutine compute_effective_viscosity

!****************************************************************************

  subroutine compute_effective_viscosity_L1L2(whichefvs,        efvs_constant,      &
                                              nz,               sigma,              &
                                              nNodesPerElement, phi,                &
                                              dphi_dx,          dphi_dy,            &
                                              uvel,             vvel,               &
                                              stagthck,                             &
                                              dsdx,             dsdy,               &
                                              flwa,             flwafact,           &
                                              efvs,                                 &
                                              i, j, p )

    ! Compute the effective viscosity at each layer of an ice column corresponding
    !  to a particular quadrature point, based on the L1L2 formulation.
    ! See PGB(2012), section 2.3

    integer, intent(in) :: i, j, p

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) :: &
       whichefvs          ! method for computing effective viscosity
                          ! 0 = constant value
                          ! 1 = proportional to flow factor
                          ! 2 = nonlinear function of effective strain rate 

    real(dp), intent(in) :: &
       efvs_constant      ! constant value of effective viscosity (Pa yr)
                          ! (used for option HO_EFVS_CONSTANT)

    integer, intent(in) ::  &
       nz,               &! number of vertical levels at which velocity is computed
                          ! Note: The number of layers (or elements in a column) is nz-1
       nNodesPerElement   ! number of nodes per element, = 4 for 2D rectangular faces

    real(dp), dimension(nz), intent(in) ::    &
       sigma              ! sigma vertical coordinate

    real(dp), dimension(nNodesPerElement), intent(in) ::  &
       phi,           &   ! basic functions at this quadrature point
       dphi_dx, dphi_dy   ! derivatives of basis functions at this quadrature point

    real(dp), dimension(nNodesPerElement), intent(in) ::  &
       uvel, vvel,       &! current guess for basal velocity at cell vertices (m/yr)
       dsdx, dsdy,       &! upper surface elevation gradient at vertices (m/m)
       stagthck           ! ice thickness at vertices

    real(dp), dimension(nz-1), intent(in) ::  &
       flwa,             &! temperature-based flow factor A at each layer of this cell column
                          ! units: Pa^{-n} yr^{-1}
       flwafact           ! temperature-based flow factor for this element, 0.5 * A^{-1/n}
                          ! units: Pa yr^{1/n}  (used for option HO_EFVS_FLOWFACT)

    real(dp), dimension(nz-1), intent(out) ::   &
       efvs               ! effective viscosity of each layer corresponding to this quadrature point (Pa yr)
                          ! computed as 1 / (2*A*tau_eff^{(n-1)/2})
                          !           = 1 / (2*A*tau_eff^2) given n = 3
                          ! where tau_eff^2 = tau_parallel^2 + tau_perp^2
 
    !----------------------------------------------------------------
    ! Local parameters
    !----------------------------------------------------------------

    real(dp), parameter ::   &
!!       effstrain_min = 1.d-20*scyr,     &! minimum value of effective strain rate, yr^{-1}
                                           ! GLAM uses 1.d-20 s^{-1} for minimum effective strain rate
       effstrain_min = 1.d-8,     &! minimum value of effective strain rate, yr^{-1}
                                   ! Mauro Perego suggests 1.d-8 yr^{-1}
       p_effstr = (1.d0 - real(gn,dp)) / real(gn,dp)    ! exponent (1-n)/n in effective viscosity relation
                                                               
    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    real(dp) ::            &
       du_dx, du_dy,       & ! horizontal strain rate components at this quadrature point, yr^{-1}
       dv_dx, dv_dy,       &
       ds_dx, ds_dy,       & ! gradient of upper surface elevation at this QP (m/m)
       thck,               & ! ice thickness (m) at this QP
       effstrain,          & ! effective strain rate at QP, yr^{-1}
       effstrainsq,        & ! square of effective strain rate
       tau_parallel,       & ! norm of tau_parallel at each layer of this cell column, 
                             !  where |tau_parallel|^2 = tau_xx^2 + tau_yy^2 + tau_xx*tau_yy + tau_xy^2
                             !  See PGB(2012), eq. 17 and 20
       tau_perp,           & ! norm of tau_perp at a given layer of a cell column,
                             !  where |tau_perp|^2 = [rhoi*grav*(s-z)*|grad(s)|]^2
       grads,              & ! norm of sfc elevation gradient at this QP, sqrt(ds_dx^2 + ds_dy^2)
       depth                 ! distance (m) from surface to level k at this QP 

    real(dp) :: a, b, c, rootA, rootB   ! terms in cubic equation

    integer :: n, k

    select case(whichefvs)

    case(HO_EFVS_CONSTANT)

       ! Steve Price recommends 10^6 to 10^7 Pa yr
       ! ISMIP-HOM Test F requires 2336041.42829 Pa yr; this is the default value set in glide_types.F90
       efvs(:) = efvs_constant

       if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest) then
          print*, 'Set efvs = constant (Pa yr):', efvs
       endif

    case(HO_EFVS_FLOWFACT)      ! set the effective viscosity to a multiple of the flow factor, 0.5*A^(-1/n)
                
       ! Set the effective strain rate (s^{-1}) based on typical velocity and length scales
       !
       ! Units: flwafact has units Pa yr^{1/n}
       !        effstrain has units yr^{-1}
       !        p_effstr = (1-n)/n 
       !                 = -2/3 for n=3
       ! Thus efvs has units Pa yr
   
       effstrain = vel_scale/len_scale * scyr  ! typical strain rate, yr^{-1}
       efvs(:) = flwafact(:) * effstrain**p_effstr  

       if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest) then
          print*, 'flwafact, effstrain (yr-1), efvs (Pa yr)=', flwafact, effstrain, efvs
       endif

    case(HO_EFVS_NONLINEAR)    ! compute effective viscosity based on effective strain rate

       du_dx = 0.d0
       du_dy = 0.d0
       dv_dx = 0.d0
       dv_dy = 0.d0
       ds_dx = 0.d0
       ds_dy = 0.d0
       thck  = 0.d0

       do n = 1, nNodesPerElement

          du_dx = du_dx + dphi_dx(n)*uvel(n)
          du_dy = du_dy + dphi_dy(n)*uvel(n)

          dv_dx = dv_dx + dphi_dx(n)*vvel(n)
          dv_dy = dv_dy + dphi_dy(n)*vvel(n)

          ds_dx = ds_dx + phi(n)*dsdx(n)
          ds_dy = ds_dy + phi(n)*dsdy(n)

          thck = thck + phi(n)*stagthck(n)

       enddo

       ! Compute effective strain rate at this quadrature point (PGB 2012, eq. 17)

       effstrainsq = effstrain_min**2          &
                   + du_dx**2 + dv_dy**2 + du_dx*dv_dy + 0.25d0*(dv_dx + du_dy)**2
       effstrain = sqrt(effstrainsq)

       if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest .and. p==ptest) then
          print*, ' '
          print*, 'i, j, p, effstrain (yr-1):', i, j, p, effstrain
          print*, 'du_dx, du_dy =', du_dx, du_dy
          print*, 'dv_dx, dv_dy =', dv_dx, dv_dy
          print*, 'ds_dx, ds_dy =', ds_dx, ds_dy
!          print*, 'n, phi, dphi_dx, dphi_dy:'
!          do n = 1, nNodesPerElement_2d
!             print*, n, phi(n), dphi_dx(n), dphi_dy(n)
!          enddo
       endif

       !---------------------------------------------------------------------------
       ! Solve for tau_parallel in the relation (PGB 2012, eq. 22)
       !
       !     effstrain = A * (tau_parallel^2 + tau_perp^2)^{(n-1)/2} * tau_parallel
       !
       !     where tau_perp^2 = [(pg)*(s-z)*|grad(s)|]^2 = SIA stress
       !              grad(s) = sqrt(ds_dx^2 + ds_dy^2)
       !                    n = 3, so we have a cubic equation
       !
       ! This relation can be written as a cubic equation of the form
       !
       !            x^3 + a*x + b = 0,
       !
       ! where for this problem, x = tau_parallel > 0,
       !                         a = tau_perp^2 >= 0,
       !                         b = -effstrain/A < 0.
       !
       ! If (b^2)/4 + (a^3)/27 > 0, then there is one real root A + B, where
       ! 
       !     A = [-b/2 + sqrt((b^2)/4 + (a^3)/27)]^(1/3)
       !     B = -[b/2 + sqrt((b^2)/4 + (a^3)/27)]^(1/3)
       !  
       ! There is also a pair of complex conjugate roots that are not of interest here.
       !
       ! Note: If a^3/27 << b^2/4 (as can happen if |grad(s)| is small), then the
       !       bracketed term in B is given to a good approximation by 
       !
       !       b/2 + (|b|/2)*(1 + 2a^3/(27b^2)) = a^3 / (27|b|).
       !
       ! Hence B = -a / (3 * |b|^(1/3)).
       !
       ! We use the alternate expression for B when a^3/27 < 1.d-6 * b^2/4,
       !  so as to avoid roundoff error from subtracting two large numbers of nearly
       !  the same size. 
       !---------------------------------------------------------------------------
       !TODO - Code an iterative solution for tau_parallel, for n /= 3.
       !TODO - Replace sigma with stagsigma?  Not sure if depth should be at layer midpt or base

       do k = 1, nz-1   ! loop over layers
          depth = thck * sigma(k+1)
          grads = sqrt(ds_dx**2 + ds_dy**2)
          tau_perp = rhoi*grav*depth*grads
          a = tau_perp**2
          b = -effstrain / flwa(k)
          c = sqrt(b**2/4.d0 + a**3/27.d0)
          rootA = (-b/2.d0 + c)**(1.d0/3.d0)
          if (a**3/(27.d0) > 1.d-6 * (b**2/4.d0)) then
             rootB = -(b/2.d0 + c)**(1.d0/3.d0)
          else    ! b/2 + c is small; compute solution to first order without subtracting two large, nearly equal numbers
             rootB = -a / (3.d0*(abs(b))**(1.d0/3.d0))
          endif
          tau_parallel = rootA + rootB
          efvs(k) = 1.d0 / (2.d0 * flwa(k) * (tau_parallel**2 + tau_perp**2))  ! given n = 3

          !WHL - debug
          if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest .and. p==ptest) then
             print*, 'i, j, k, p =', i, j, k, p
!             print*, 'a, b, c:', a, b, c
!             print*, '-b/2 + c, -b/2 - c:', -b/2 + c, -b/2 - c
!             print*, 'roots A, B:', rootA, rootB
!             print*, 'tau_perp, tau_parallel:', tau_perp, tau_parallel
!             print*, 'flwa:', flwa(k)
             print*, 'flwafact, effstrain, efvs_BP, efvs:', 0.5d0*flwa(k)**(-1.d0/3.d0), effstrain,  &
                                                            0.5d0*flwa(k)**(-1.d0/3.d0) * effstrain**(-2.d0/3.d0), efvs(k)
          endif

       enddo   ! k

    end select

  end subroutine compute_effective_viscosity_L1L2

!****************************************************************************

  subroutine compute_effective_viscosity_diva(whichefvs,        efvs_constant,      &
                                              nz,               stagsigma,          &
                                              nNodesPerElement, phi,                &
                                              dphi_dx,          dphi_dy,            &
                                              uvel,             vvel,               &
                                              btractx,          btracty,            &
                                              stagthck,                             &
                                              flwa,             flwafact,           &
                                              efvs,                                 &
                                              i, j, p )
    
    ! Compute the effective viscosity at each layer of an ice column corresponding
    !  to a particular quadrature point, based on the depth-integrated formulation.
    ! See Goldberg(2011) for details.

    integer, intent(in) :: i, j, p

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) :: &
       whichefvs          ! method for computing effective viscosity
                          ! 0 = constant value
                          ! 1 = proportional to flow factor
                          ! 2 = nonlinear function of effective strain rate 

    real(dp), intent(in) :: &
       efvs_constant      ! constant value of effective viscosity (Pa yr)
                          ! (used for option HO_EFVS_CONSTANT)

    integer, intent(in) ::  &
       nz,               &! number of vertical levels at which velocity is computed
                          ! Note: The number of layers (or elements in a column) is nz-1
       nNodesPerElement   ! number of nodes per element, = 4 for 2D rectangular faces

    real(dp), dimension(nz-1), intent(in) ::    &
       stagsigma          ! staggered sigma vertical coordinate

    real(dp), dimension(nNodesPerElement), intent(in) ::  &
       phi,           &   ! basic functions at this quadrature point
       dphi_dx, dphi_dy   ! derivatives of basis functions at this quadrature point

    real(dp), dimension(nNodesPerElement), intent(in) ::  &
       uvel, vvel,       &! current guess for depth_integrated mean velocity at cell vertices (m/yr)
       btractx, btracty, &! components of basal traction (Pa)
       stagthck           ! ice thickness at vertices

    real(dp), dimension(nz-1), intent(in) ::  &
       flwa,             &! temperature-based flow factor A at each layer of this cell column
                          ! units: Pa^{-n} yr^{-1}
       flwafact           ! temperature-based flow factor for this element, 0.5 * A^{-1/n}
                          ! units: Pa yr^{1/n}  (used for option HO_EFVS_FLOWFACT)

    !WHL - intent(out) if solving cubic, but (inout) if using old efvs in calculation
    real(dp), dimension(nz-1), intent(inout) ::   &
       efvs               ! effective viscosity of each layer corresponding to this quadrature point (Pa yr)

    !----------------------------------------------------------------
    ! Local parameters
    !----------------------------------------------------------------

    real(dp), parameter ::   &
!!       effstrain_min = 1.d-20*scyr,     &! minimum value of effective strain rate (yr^{-1})
                                           ! GLAM uses 1.d-20 s^{-1} for minimum effective strain rate
       effstrain_min = 1.d-8,     &! minimum value of effective strain rate (yr^{-1})
                                   ! Mauro Perego suggests 1.d-8 yr^{-1}
       p_effstr  = (1.d0 - real(gn,dp))/real(gn,dp), &! exponent (1-n)/n in effective viscosity relation
       p2_effstr = p_effstr/2                         ! exponent (1-n)/(2n) in effective viscosity relation
                                                               
    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    real(dp) ::            &
       du_dx, du_dy,       & ! horizontal strain rate components at this quadrature point (yr^{-1})
       dv_dx, dv_dy,       &
       taux,  tauy,        & ! basal shear stress components at this QP (Pa)
       thck,               & ! ice thickness (m) at this QP
       effstrain,          & ! effective strain rate at QP  (yr^{-1})
       effstrainsq,        & ! square of effective strain rate
       depth                 ! distance (m) from surface to layer k at this QP 

    real(dp) :: facta, factb, a, b, c, rootA, rootB   ! terms in cubic equation

    integer :: n, k
    real(dp) :: du_dz, dv_dz

    !WHL - For ISMIP-HOM, the cubic solve is not robust.  It leads to oscillations
    !      in successive iterations between uvel_2d/vvel_2d and btractx/btracty
    !TODO - Remove the cubic solve for efvs, unless we find a way to make it robust?
    logical, parameter :: cubic = .false.

    select case(whichefvs)

    case(HO_EFVS_CONSTANT)

       ! Steve Price recommends 10^6 to 10^7 Pa yr
       ! ISMIP-HOM Test F requires 2336041.42829 Pa yr; this is the default value set in glide_types.F90
       efvs(:) = efvs_constant

       if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest) then
          print*, 'Set efvs = constant (Pa yr):', efvs
       endif

    case(HO_EFVS_FLOWFACT)      ! set the effective viscosity to a multiple of the flow factor, 0.5*A^(-1/n)
                
       ! Set the effective strain rate (s^{-1}) based on typical velocity and length scales
       !
       ! Units: flwafact has units Pa yr^{1/n}
       !        effstrain has units yr^{-1}
       !        p_effstr = (1-n)/n 
       !                 = -2/3 for n=3
       ! Thus efvs has units Pa yr
   
       effstrain = vel_scale/len_scale * scyr  ! typical strain rate, yr^{-1}
       efvs(:) = flwafact(:) * effstrain**p_effstr  

       if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest) then
          print*, 'flwafact, effstrain (yr-1), efvs (Pa yr)=', flwafact, effstrain, efvs
       endif

    case(HO_EFVS_NONLINEAR)    ! compute effective viscosity based on effective strain rate

       du_dx = 0.d0
       du_dy = 0.d0
       dv_dx = 0.d0
       dv_dy = 0.d0
       thck  = 0.d0
       taux  = 0.d0
       tauy  = 0.d0

       do n = 1, nNodesPerElement

          du_dx = du_dx + dphi_dx(n)*uvel(n)
          du_dy = du_dy + dphi_dy(n)*uvel(n)

          dv_dx = dv_dx + dphi_dx(n)*vvel(n)
          dv_dy = dv_dy + dphi_dy(n)*vvel(n)

          taux = taux + phi(n)*btractx(n)
          tauy = tauy + phi(n)*btracty(n)

          thck = thck + phi(n)*stagthck(n)

       enddo

    if (cubic) then

       ! Compute effective strain rate (squared) at this quadrature point

       effstrainsq = effstrain_min**2          &
                   + du_dx**2 + dv_dy**2 + du_dx*dv_dy + 0.25d0*(dv_dx + du_dy)**2

       if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest .and. p==ptest) then
          print*, ' '
          print*, 'i, j, p, effstrain (yr-1):', i, j, p, sqrt(effstrainsq)
          print*, 'du_dx, du_dy =', du_dx, du_dy
          print*, 'dv_dx, dv_dy =', dv_dx, dv_dy
          print*, 'btractx, btracty =',  btractx, btracty
          print*, 'taux, tauy =', taux, tauy
       endif

       !---------------------------------------------------------------------------
       ! Solve for efvs in the relation
       !
       ! efvs = 1/2 * A^(-1/n) * [effstrainsq + (1/4)*(u_z^2 + v_z^2)]^[(1-n)/(2n)]
       !
       ! where effstrainsq = du_dx**2 + dv_dy**2 + du_dx*dv_dy + (1/4)*(dv_dx + du_dy)**2
       !                     + small regularization term
       !       u_z = tau_x*(s-z) / (H*efvs)
       !       v_z = tau_y*(s-z) / (H*efvs)
       !
       !       tau_x = beta*u_b = beta_eff*u
       !       tau_y = beta*v_b = beta_eff*v
       !
       !       (u,v) is the depth-averaged mean velocity
       !
       ! For n = 3, this relation can be written as a cubic equation of the form
       !
       !       x^3 + a*x + b = 0,
       !
       ! where x = efvs
       !       a = [(tau_x^2 + tau_y^2)*(s-z)^2 / (4*H^2*effstrainsq) >= 0
       !       b = -1/(8*A*effstrainsq) < 0
       !
       ! See comments in compute_effective_viscosity_L1L2 for more details on the cubic solve.
       !
       ! NOTE: This scheme does not reliably converge.
       !
       !       The problem is that taux and tauy are proportional to beta_eff, which is
       !        a function of the old viscosity.  Mixing the old and new viscosity in the
       !        expression for vertical shear can lead to oscillations.
       !---------------------------------------------------------------------------

       facta = (taux**2 + tauy**2) / (4.d0 * thck**2 * effstrainsq)
       factb = -1.d0 / (8.d0 * effstrainsq)
       do k = 1, nz-1   ! loop over layers
          depth = thck * stagsigma(k)
          a = facta * depth**2
          b = factb / flwa(k)
          c = sqrt(b**2/4.d0 + a**3/27.d0)
          rootA = (-b/2.d0 + c)**(1.d0/3.d0)
          if (a**3/(27.d0) > 1.d-6 * (b**2/4.d0)) then
             rootB = -(b/2.d0 + c)**(1.d0/3.d0)
          else    ! b/2 + c is small; compute solution to first order without subtracting two large, nearly equal numbers
             rootB = -a / (3.d0*(abs(b))**(1.d0/3.d0))
          endif
          efvs(k) = rootA + rootB

          if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest .and. p==ptest) then
             print*, ' '
             print*, 'i, j, k, p, depth =', i, j, k, p, depth
             print*, 'a, b, c:', a, b, c
             print*, '-b/2 + c, -b/2 - c:', -b/2 + c, -b/2 - c
             print*, 'roots A, B:', rootA, rootB
             print*, 'flwa:', flwa(k)
             effstrain = sqrt(effstrainsq)
             print*, 'flwafact, effstrain, efvs_SSA, efvs:', flwafact(k), effstrain,  &
                                                             flwafact(k)*effstrain**(-2.d0/3.d0), efvs(k)
          endif

       enddo   ! k

    else  ! solve for efvs, using the old value of efvs to estimate the vertical strain rates

       do k = 1, nz-1   ! loop over layers
          if (efvs(k)==0.d0) then
             efvs(k) = flwafact(k) * effstrain_min**p_effstr  ! efvs associated with minimum strain rate
          endif
          du_dz = taux * stagsigma(k) / efvs(k)   ! old value of efvs on RHS
          dv_dz = tauy * stagsigma(k) / efvs(k)
          effstrainsq = effstrain_min**2          &
                      + du_dx**2 + dv_dy**2 + du_dx*dv_dy + 0.25d0*(dv_dx + du_dy)**2  &
                      + 0.25d0 * (du_dz**2 + dv_dz**2)
          efvs(k) = flwafact(k) * effstrainsq**p2_effstr
       enddo

    endif   ! cubic

    end select

  end subroutine compute_effective_viscosity_diva

!****************************************************************************

  subroutine compute_element_matrix(whichapprox, nNodesPerElement,     &
                                    wqp,         detJ,                 &
                                    efvs,                              &
                                    dphi_dx,     dphi_dy,    dphi_dz,  &
                                    Kuu,         Kuv,                  &
                                    Kvu,         Kvv,                  &
                                    i, j, k, p)

    !------------------------------------------------------------------
    ! Increment the stiffness matrices Kuu, Kuv, Kvu, Kvv with the
    ! contribution from a particular quadrature point, 
    ! based on the Blatter-Pattyn first-order equations.
    !
    ! Note: Elements can be either 2D or 3D
    !------------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) :: i, j, k, p

    integer, intent(in) :: &
         whichapprox     ! which Stokes approximation to use (BP, SIA, SSA)

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

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    real(dp) :: efvs_factor
    integer :: nr, nc

    if (verbose_matrix .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, ' '
       print*, 'Increment element matrix, i, j, k, p =', i, j, k, p
    endif

    ! Increment the element stiffness matrices for the appropriate approximation.

    !Note: Scaling by volume such that detJ/vol0 is close to unity
    efvs_factor = efvs * wqp * detJ/vol0
    
    if (verbose_matrix .and. this_rank==rtest .and. i==itest .and. j==jtest .and. &
         k==ktest .and. p==ptest) then
       print*, ' '
       print*, 'i, j, k, p:', i, j, k, p
       print*, 'efvs, wqp, detJ/vol0 =', efvs, wqp, detJ/vol0
       print*, 'dphi_dz(1) =', dphi_dz(1)
       print*, 'dphi_dx(1) =', dphi_dx(1)
       print*, 'Kuu dphi/dz increment(1,1) =', efvs_factor*dphi_dz(1)*dphi_dz(1)
       print*, 'Kuu dphi/dx increment(1,1) =', efvs_factor*4.d0*dphi_dx(1)*dphi_dx(1)
    endif

    if (whichapprox == HO_APPROX_SIA) then

       do nc = 1, nNodesPerElement      ! columns of K
          do nr = 1, nNodesPerElement   ! rows of K

             Kuu(nr,nc) = Kuu(nr,nc) + efvs_factor * (dphi_dz(nr)*dphi_dz(nc))
             Kvv(nr,nc) = Kvv(nr,nc) + efvs_factor * (dphi_dz(nr)*dphi_dz(nc))             

          enddo  ! row
       enddo     ! column

    elseif (whichapprox == HO_APPROX_SSA) then

       do nc = 1, nNodesPerElement      ! columns of K
          do nr = 1, nNodesPerElement   ! rows of K

             Kuu(nr,nc) = Kuu(nr,nc) + efvs_factor * (4.d0*dphi_dx(nr)*dphi_dx(nc) + dphi_dy(nr)*dphi_dy(nc))
             Kuv(nr,nc) = Kuv(nr,nc) + efvs_factor * (2.d0*dphi_dx(nr)*dphi_dy(nc) + dphi_dy(nr)*dphi_dx(nc))
             Kvu(nr,nc) = Kvu(nr,nc) + efvs_factor * (2.d0*dphi_dy(nr)*dphi_dx(nc) + dphi_dx(nr)*dphi_dy(nc))
             Kvv(nr,nc) = Kvv(nr,nc) + efvs_factor * (4.d0*dphi_dy(nr)*dphi_dy(nc) + dphi_dx(nr)*dphi_dx(nc))

          enddo
       enddo

    else   ! Blatter-Pattyn higher-order
           ! The terms in parentheses can be derived from PGB 2012, eq. 13 and 15.
           ! The factor of 2 in front of efvs has been absorbed into the quantities in parentheses.

       do nc = 1, nNodesPerElement      ! columns of K
          do nr = 1, nNodesPerElement   ! rows of K

             Kuu(nr,nc) = Kuu(nr,nc) + efvs_factor *                                             &
                                    ( 4.d0*dphi_dx(nr)*dphi_dx(nc) + dphi_dy(nr)*dphi_dy(nc)     &
                                    + dphi_dz(nr)*dphi_dz(nc) )

             Kuv(nr,nc) = Kuv(nr,nc) + efvs_factor *                                             &
                                     (2.d0*dphi_dx(nr)*dphi_dy(nc) + dphi_dy(nr)*dphi_dx(nc))

             Kvu(nr,nc) = Kvu(nr,nc) + efvs_factor *                                             &
                                     (2.d0*dphi_dy(nr)*dphi_dx(nc) + dphi_dx(nr)*dphi_dy(nc))

             Kvv(nr,nc) = Kvv(nr,nc) + efvs_factor *                                             &
                                    ( 4.d0*dphi_dy(nr)*dphi_dy(nc) + dphi_dx(nr)*dphi_dx(nc)        &
                                    + dphi_dz(nr)*dphi_dz(nc) )

          enddo  ! nr (rows)
       enddo     ! nc (columns)

    endif  ! whichapprox

  end subroutine compute_element_matrix

!****************************************************************************

  subroutine element_to_global_matrix_3d(nx,           ny,          nz,          &
                                         iElement,     jElement,    kElement,    &
                                         Kuu,          Kuv,                      &
                                         Kvu,          Kvv,                      &
                                         Auu,          Auv,                      &
                                         Avu,          Avv)
             
    ! Sum terms of element matrix K into dense assembled matrix A
    ! K is partitioned into Kuu, Kuv, Kvu, and Kvv, and similarly for A.

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz                       ! number of vertical levels where velocity is computed

    integer, intent(in) ::   &
       iElement, jElement, kElement     ! i, j and k indices for this element

    real(dp), dimension(nNodesPerElement_3d,nNodesPerElement_3d), intent(in) ::  &
       Kuu, Kuv, Kvu, Kvv       ! element matrix

    real(dp), dimension(nNodeNeighbors_3d,nz,nx-1,ny-1), intent(inout) ::    &
       Auu, Auv, Avu, Avv       ! assembled matrix

    integer :: i, j, k, m
    integer :: iA, jA, kA
    integer :: n, nr, nc

    if (verbose_matrix .and. this_rank==rtest .and. iElement==itest .and. jElement==jtest .and. kElement==ktest) then
       print*, 'Element i, j, k:', iElement, jElement, kElement 
       print*, 'Rows of Kuu:'
       do n = 1, nNodesPerElement_3d
          write(6, '(8e12.4)') Kuu(n,:)
       enddo
    endif

    !WHL - On a Mac I tried switching the loops to put nc on the outside, but 
    !      the one with nr on the outside is faster.
    do nr = 1, nNodesPerElement_3d       ! rows of K

       ! Determine row of A to be incremented by finding (k,i,j) for this node
       ! The reason for the '7' is that node 7, in the NE corner of the upper layer, has index (k,i,j).
       ! Indices for other nodes are computed relative to this node.
       i = iElement + ishift(7,nr)
       j = jElement + jshift(7,nr)
       k = kElement + kshift(7,nr)
      
       do nc = 1, nNodesPerElement_3d    ! columns of K

          ! Determine column of A to be incremented
          kA = kshift(nr,nc)           ! k index of A into which K(m,n) is summed
          iA = ishift(nr,nc)           ! similarly for i and j indices 
          jA = jshift(nr,nc)           ! these indices can take values -1, 0 and 1
          m = indxA_3d(iA,jA,kA)

          ! Increment A
          Auu(m,k,i,j) = Auu(m,k,i,j) + Kuu(nr,nc)
          Auv(m,k,i,j) = Auv(m,k,i,j) + Kuv(nr,nc)
          Avu(m,k,i,j) = Avu(m,k,i,j) + Kvu(nr,nc)
          Avv(m,k,i,j) = Avv(m,k,i,j) + Kvv(nr,nc)

       enddo     ! nc

    enddo        ! nr

  end subroutine element_to_global_matrix_3d

!****************************************************************************

  subroutine element_to_global_matrix_2d(nx,           ny,        &
                                         iElement,     jElement,  &
                                         Kuu,          Kuv,       &
                                         Kvu,          Kvv,       &
                                         Auu,          Auv,       &
                                         Avu,          Avv)

    ! Sum terms of element matrix K into dense assembled matrix A
    ! K is partitioned into Kuu, Kuv, Kvu, and Kvv, and similarly for A.

    integer, intent(in) ::   &
       nx, ny                   ! horizontal grid dimensions

    integer, intent(in) ::   &
       iElement, jElement       ! i and j indices for this element

    real(dp), dimension(nNodesPerElement_2d,nNodesPerElement_2d), intent(in) ::  &
       Kuu, Kuv, Kvu, Kvv       ! element matrix

    real(dp), dimension(nNodeNeighbors_2d,nx-1,ny-1), intent(inout) ::    &
       Auu, Auv, Avu, Avv       ! assembled matrix

    integer :: i, j, m
    integer :: iA, jA
    integer :: n, nr, nc

    if (verbose_matrix .and. this_rank==rtest .and. iElement==itest .and. jElement==jtest) then
       print*, 'Element i, j:', iElement, jElement 
       print*, 'Rows of Kuu:'
       do n = 1, nNodesPerElement_2d
          write(6, '(8e12.4)') Kuu(n,:)
       enddo
    endif

    do nr = 1, nNodesPerElement_2d       ! rows of K

       ! Determine row of A to be incremented by finding (i,j) for this node
       ! The reason for the '3' is that node 3, in the NE corner of this gridcell, has index (i,j).
       ! Indices for other nodes are computed relative to this node.
       i = iElement + ishift(3,nr)
       j = jElement + jshift(3,nr)
      
       do nc = 1, nNodesPerElement_2d    ! columns of K

          ! Determine column of A to be incremented
          iA = ishift(nr,nc)           ! similarly for i and j indices 
          jA = jshift(nr,nc)           ! these indices can take values -1, 0 and 1
          m = indxA_2d(iA,jA)

          ! Increment A
          Auu(m,i,j) = Auu(m,i,j) + Kuu(nr,nc)
          Auv(m,i,j) = Auv(m,i,j) + Kuv(nr,nc)
          Avu(m,i,j) = Avu(m,i,j) + Kvu(nr,nc)
          Avv(m,i,j) = Avv(m,i,j) + Kvv(nr,nc)

       enddo     ! nc
    enddo        ! nr

  end subroutine element_to_global_matrix_2d

!****************************************************************************

  subroutine basal_sliding_bc(nx,               ny,              &
                              nNeighbors,       nhalo,           &
                              active_cell,      beta,            &
                              xVertex,          yVertex,         &
                              whichassemble_beta,                &
                              Auu,              Avv)

    !------------------------------------------------------------------------
    ! Increment the Auu and Avv matrices with basal traction terms.
    ! Do a surface integral over all basal faces that contain at least one node with a stress BC. 
    ! (Not Dirichlet or free-slip)
    ! Note: Basal Dirichlet BCs are enforced after matrix assembly. 
    !
    ! Assume a sliding law of the form:
    !   tau_x = -beta*u
    !   tau_y = -beta*v
    ! where beta is defined at vertices (and beta may depend
    ! on the velocity from a previous iteration).
    !
    ! Note: The input beta field should already have been weighted by f_ground. We should have
    !       beta = 0 for floating ice (f_ground = 0). If using a GLP, then beta will
    !       have less than its full value for partially floating ice (0 < f_ground < 1). 
    !------------------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nNeighbors,              &    ! number of neighbors of each node (used for first dimension of Auu/Avv)
                                     ! = 27 for 3D solve, = 9 for 2D solve
       nhalo                         ! number of halo layers

    logical, dimension(nx,ny), intent(in) ::  &
       active_cell                   ! true if cell contains ice and borders a locally owned vertex

    real(dp), dimension(nx-1,ny-1), intent(in) ::    &
       beta                          ! basal traction field (Pa/(m/yr)) at cell vertices
                                     ! typically = beta_internal (beta weighted by f_ground)
                                     ! = beta_eff for DIVA

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xVertex, yVertex     ! x and y coordinates of vertices

    integer, intent(in) :: &
       whichassemble_beta   ! = 0 for standard finite element computation of basal forcing terms
                            ! = 1 for computation that uses only the local value of beta at each node

    real(dp), dimension(nNeighbors,nx-1,ny-1), intent(inout) ::  &
       Auu, Avv             ! parts of stiffness matrix (basal layer only)

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j, n, p, nr, nc, iA, jA, m, ii, jj

    real(dp), dimension(nNodesPerElement_2d) ::   &
       x, y,        & ! Cartesian coordinates of basal nodes
       b              ! beta at basal nodes

    !TODO - These are not currently used except as dummy arguments
    real(dp), dimension(nNodesPerElement_2d) ::   &
       dphi_dx_2d, dphi_dy_2d    ! derivatives of basis functions, evaluated at quad pts

    real(dp) ::   &
       beta_qp,     & ! beta evaluated at quadrature point
       detJ           ! determinant of Jacobian for the transformation
                      !  between the reference element and true element

    real(dp), dimension(nNodesPerElement_2d, nNodesPerElement_2d) ::   &
       Kuu, Kvv       ! components of element matrix associated with basal sliding

    if (verbose_basal .and. this_rank==rtest) then
       print*, 'In basal_sliding_bc: itest, jtest, rank =', itest, jtest, rtest
    endif

    ! Sum over elements in active cells 
    ! Loop over all cells that contain locally owned vertices
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

          ! loop over quadrature points

          do p = 1, nQuadPoints_2d

             ! Compute basis function derivatives and det(J) for this quadrature point
             ! For now, pass in i, j, k, p for debugging
             !TODO - Modify this subroutine so that the output derivatives are optional?

             call get_basis_function_derivatives_2d(x(:),             y(:),               & 
                                                    dphi_dxr_2d(:,p), dphi_dyr_2d(:,p),   &   
                                                    dphi_dx_2d(:),    dphi_dy_2d(:),      &
                                                    detJ, i, j, p)
          
             ! Evaluate beta at this quadrature point
             ! Standard finite-element treatment is to take a phi-weighted sum over neighboring vertices.
             ! For local beta, use the value at the nearest vertex.
             !  (Note that vertex numbering is the same as QP numbering, CCW from 1 to 4 starting at SW corner.)
 
             if (whichassemble_beta == HO_ASSEMBLE_BETA_LOCAL) then
                beta_qp = b(p)
             else
                beta_qp = 0.d0
                do n = 1, nNodesPerElement_2d
                   beta_qp = beta_qp + phi_2d(n,p) * b(n)
                enddo
             endif

             if (verbose_basal .and. this_rank==rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'Increment basal traction, i, j, p =', i, j, p
                print*, 'beta_qp =', beta_qp
                print*, 'detJ/vol0 =', detJ/vol0
             endif

             ! Compute the element matrix for this quadrature point
             ! (Note volume scaling)

             Kuu(:,:) = 0.d0

             if (whichassemble_beta == HO_ASSEMBLE_BETA_LOCAL) then  ! Use the value at the nearest vertex
                                      ! Then Kuu is diagonal, so the traction parameter at a vertex depends only on beta at that vertex
                Kuu(p,p) = beta_qp * (detJ/vol0)   

             else

                do nc = 1, nNodesPerElement_2d      ! columns of K
                   do nr = 1, nNodesPerElement_2d   ! rows of K
                      Kuu(nr,nc) = Kuu(nr,nc) + beta_qp * wqp_2d(p) * detJ/vol0 * phi_2d(nr,p)*phi_2d(nc,p)
                   enddo  ! m (rows)
                enddo     ! n (columns)

             endif        ! local beta

             !Note: Is this true for all sliding laws?
             Kvv(:,:) = Kuu(:,:)

             ! Insert terms of basal element matrices into global matrices Auu and Avv

             do nr = 1, nNodesPerElement_2d     ! rows of K

                ! Determine (i,j) for this node
                ! The reason for the '3' is that node 3, in the NE corner of the cell, has horizontal indices (i,j).
                ! Indices for other nodes are computed relative to this node.

                ii = i + ishift(3,nr)
                jj = j + jshift(3,nr)
      
                do nc = 1, nNodesPerElement_2d ! columns of K

                   iA = ishift(nr,nc)          ! iA index of A into which K(nr,nc) is summed
                   jA = jshift(nr,nc)          ! similarly for jA

                   if (nNeighbors == nNodeNeighbors_3d) then  ! 3D problem
                      m = indxA_3d(iA,jA,0)
                   else  ! 2D problem
                      m = indxA_2d(iA,jA)
                   endif

                   Auu(m,ii,jj) = Auu(m,ii,jj) + Kuu(nr,nc)
                   Avv(m,ii,jj) = Avv(m,ii,jj) + Kvv(nr,nc)

                enddo     ! nc
             enddo        ! nr

             if (verbose_basal .and. this_rank==rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'i, j =', i, j
                print*, 'Kuu:'
                do nr = 1, nNodesPerElement_2d
                   print*, nr, Kuu(nr,:)
                enddo
                print*, ' '
                print*, 'rowsum(Kuu):'
                do nr = 1, nNodesPerElement_2d
                   print*, nr, sum(Kuu(nr,:))
                enddo
                print*, ' '
                print*, 'sum(Kuu):', sum(Kuu(:,:))
             endif

          enddo   ! nQuadPoints_2d

       endif      ! active_cell

    enddo         ! i
    enddo         ! j

    if (verbose_basal .and. this_rank==rtest) then
       i = itest
       j = jtest
       if (nNeighbors == nNodeNeighbors_3d) then  ! 3D problem
          m = indxA_3d(0,0,0)
          print*, 'Diagonal index =', m
       else
          m = indxA_2d(0,0)
          print*, 'Diagonal index =', m
       endif
       print*, ' '
       print*, 'New Auu diagonal:', Auu(m,i,j)
       print*, 'New Avv diagonal:', Avv(m,i,j)
    endif

  end subroutine basal_sliding_bc

!****************************************************************************

  subroutine dirichlet_boundary_conditions_3d(nx,              ny,               &
                                              nz,              nhalo,            &
                                              active_vertex,                     &
                                              umask_dirichlet, vmask_dirichlet,  &
                                              uvel,            vvel,             &
                                              Auu,             Auv,              &
                                              Avu,             Avv,              &
                                              bu,              bv)

    !----------------------------------------------------------------
    ! Modify the global matrix and RHS for Dirichlet boundary conditions,
    !  where uvel and vvel are prescribed at certain nodes.
    ! For each such node, we zero out the row, except for setting the diagonal term to 1.
    ! We also zero out the column, moving terms containing uvel/vvel to the rhs.
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

      integer, dimension(nz,nx-1,ny-1), intent(in) ::  &
       umask_dirichlet,   &! Dirichlet mask for u velocity (if true, u is prescribed)
       vmask_dirichlet     ! Dirichlet mask for v velocity (if true, v is prescribed)

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       uvel, vvel          ! velocity components

    real(dp), dimension(nNodeNeighbors_3d,nz,nx-1,ny-1), intent(inout) ::   &
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

    ! Loop over all vertices that border locally owned vertices.
    ! Locally owned vertices are (nhalo+1:nx-nhalo, nhalo+1:ny-nhalo)
    !Note: Need nhalo >= 2 so as not to step out of bounds.

     do j = nhalo, ny-nhalo+1
        do i = nhalo, nx-nhalo+1
          if (active_vertex(i,j)) then
             do k = 1, nz

                if (umask_dirichlet(k,i,j) == 1) then

                   ! set the rhs to the prescribed velocity
                   bu(k,i,j) = uvel(k,i,j)

                   ! loop through matrix values in the rows associated with this node
                   ! (Auu contains one row, Avu contains a second row)
                   do kA = -1,1
                   do jA = -1,1
                   do iA = -1,1

                      if ( (k+kA >= 1 .and. k+kA <= nz)         &
                                      .and.                     &
                           (i+iA >= 1 .and. i+iA <= nx-1)       &
                                      .and.                     &
                           (j+jA >= 1 .and. j+jA <= ny-1) ) then

                         if (iA==0 .and. jA==0 .and. kA==0) then  ! main diagonal

                            ! Set Auu = 1 on the main diagonal
                            ! Set Auv term = 0; this term is off-diagonal for the fully assembled matrix
                            ! Set Avu term = 0 to preserve matrix symmetry (given that Auv term = 0)
                            m = indxA_3d(0,0,0)
                            Auu(m,k,i,j) = 1.d0
                            Auv(m,k,i,j) = 0.d0
                            Avu(m,k,i,j) = 0.d0

                            !TODO - Set bu above, outside iA/jA loop
                            ! Set the rhs to the prescribed velocity, forcing u = prescribed uvel for this vertex
!!                            bu(k,i,j) = uvel(k,i,j)
                            
                         else     ! not on the diagonal

                            ! Zero out non-diagonal matrix terms in the rows associated with this node
                            m = indxA_3d(iA,jA,kA)
                            Auu(m, k, i, j) = 0.d0
                            Auv(m, k, i, j) = 0.d0

                            ! Shift terms associated with this velocity to the rhs.
                            ! Note: The remaining operations do not change the answer, but do restore symmetry to the matrix.
                            m = indxA_3d(-iA,-jA,-kA)

                            if (umask_dirichlet(k+kA, i+iA, j+jA) /= 1) then
                               ! Move (Auu term) * uvel to rhs
                               bu(k+kA, i+iA, j+jA) = bu(k+kA, i+iA, j+jA) - Auu(m, k+kA, i+iA, j+jA) * uvel(k,i,j) 
                               Auu(m, k+kA, i+iA, j+jA) = 0.d0
                            endif

                            if (vmask_dirichlet(k+kA, i+iA, j+jA) /= 1) then
                               ! Move (Avu term) * uvel to rhs
                               bv(k+kA, i+iA, j+jA) = bv(k+kA, i+iA, j+jA) - Avu(m, k+kA, i+iA, j+jA) * uvel(k,i,j)
                               Avu(m, k+kA, i+iA, j+jA) = 0.d0
                            endif

                         endif  ! on the diagonal

                     endif     ! i+iA, j+jA, and k+kA in bounds

                  enddo        ! kA
                  enddo        ! iA
                  enddo        ! jA

                endif    ! umask_dirichlet

                if (vmask_dirichlet(k,i,j) == 1) then

                   ! set the rhs to the prescribed velocity
                   bv(k,i,j) = vvel(k,i,j)

                   ! loop through matrix values in the rows associated with this node
                   ! (Auu contains one row, Avu contains a second row)
                   do kA = -1,1
                   do jA = -1,1
                   do iA = -1,1

                      if ( (k+kA >= 1 .and. k+kA <= nz)         &
                                      .and.                     &
                           (i+iA >= 1 .and. i+iA <= nx-1)       &
                                      .and.                     &
                           (j+jA >= 1 .and. j+jA <= ny-1) ) then

                         if (iA==0 .and. jA==0 .and. kA==0) then  ! main diagonal

                            ! Set Avv = 1 on the main diagonal
                            ! Set Avu term = 0; this term is off-diagonal for the fully assembled matrix
                            ! Set Auv term = 0 to preserve matrix symmetry (given that Avu term = 0)
                            m = indxA_3d(0,0,0)

                            Auv(m,k,i,j) = 0.d0
                            Avu(m,k,i,j) = 0.d0
                            Avv(m,k,i,j) = 1.d0

                            !TODO - Set bv above, outside iA/jA loop
                            ! Set the rhs to the prescribed velocity, forcing v = prescribed vvel for this node
!!                            bv(k,i,j) = vvel(k,i,j)
                            
                         else     ! not on the diagonal

                            ! Zero out non-diagonal matrix terms in the rows associated with this node
                            m = indxA_3d(iA,jA,kA)
                            Avu(m, k, i, j) = 0.d0
                            Avv(m, k, i, j) = 0.d0

                            ! Shift terms associated with this velocity to the rhs.
                            ! Note: The remaining operations do not change the answer, but do restore symmetry to the matrix.
                            m = indxA_3d(-iA,-jA,-kA)

                            if (umask_dirichlet(k+kA, i+iA, j+jA) /= 1) then
                               ! Move (Auv term) * vvel to rhs
                               bu(k+kA, i+iA, j+jA) = bu(k+kA, i+iA, j+jA) - Auv(m, k+kA, i+iA, j+jA) * vvel(k,i,j)
                               Auv(m, k+kA, i+iA, j+jA) = 0.d0
                            endif

                            if (vmask_dirichlet(k+kA, i+iA, j+jA) /= 1) then
                               ! Move (Avv term) * vvel to rhs
                               bv(k+kA, i+iA, j+jA) = bv(k+kA, i+iA, j+jA) - Avv(m, k+kA, i+iA, j+jA) * vvel(k,i,j)
                               Avv(m, k+kA, i+iA, j+jA) = 0.d0
                            endif

                         endif  ! on the diagonal

                     endif     ! i+iA, j+jA, and k+kA in bounds

                  enddo        ! kA
                  enddo        ! iA
                  enddo        ! jA

                endif    ! vmask_dirichlet

             enddo       ! k
          endif          ! active_vertex
       enddo             ! i
    enddo                ! j

  end subroutine dirichlet_boundary_conditions_3d

!****************************************************************************

  subroutine dirichlet_boundary_conditions_2d(nx,              ny,               &
                                              nhalo,                             &
                                              active_vertex,                     &
                                              umask_dirichlet, vmask_dirichlet,  &
                                              uvel,            vvel,             &
                                              Auu,             Auv,              &
                                              Avu,             Avv,              &
                                              bu,              bv)

    !----------------------------------------------------------------
    ! Modify the global matrix and RHS for Dirichlet boundary conditions,
    !  where uvel and vvel are prescribed at certain nodes.
    ! For each such node, we zero out the row, except for setting the diagonal term to 1.
    ! We also zero out the column, moving terms containing uvel/vvel to the rhs.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nhalo                    ! number of halo layers

    logical, dimension(nx-1,ny-1), intent(in) ::  &
       active_vertex       ! true for active vertices (vertices of active cells)

    integer, dimension(nx-1,ny-1), intent(in) ::  &
       umask_dirichlet,   &! Dirichlet mask for velocity (if true, u is prescribed)
       vmask_dirichlet     ! Dirichlet mask for velocity (if true, v is prescribed)

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       uvel, vvel          ! velocity components

    real(dp), dimension(nNodeNeighbors_2d,nx-1,ny-1), intent(inout) ::   &
       Auu, Auv,    &      ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv                                    

    real(dp), dimension(nx-1,ny-1), intent(inout) ::   &
       bu, bv              ! assembled load vector, divided into 2 parts

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------
    
    integer :: i, j     ! Cartesian indices of nodes
    integer :: iA, jA   ! i and j offsets of neighboring nodes 
    integer :: m, m2

    ! Loop over all vertices that border locally owned vertices.
    ! Locally owned vertices are (nhalo+1:nx-nhalo, nhalo+1:ny-nhalo)
    !Note: Need nhalo >= 2 so as not to step out of bounds.

     do j = nhalo, ny-nhalo+1
        do i = nhalo, nx-nhalo+1
          if (active_vertex(i,j)) then

             if (umask_dirichlet(i,j) == 1) then

                ! set the rhs to the prescribed velocity
                bu(i,j) = uvel(i,j)

                ! loop through matrix values in the rows associated with this vertex
                ! (Auu contains one row, Avu contains a second row)
                do jA = -1,1
                do iA = -1,1

                   if ( (i+iA >= 1 .and. i+iA <= nx-1)       &
                                   .and.                     &
                        (j+jA >= 1 .and. j+jA <= ny-1) ) then

                      if (iA==0 .and. jA==0) then  ! main diagonal

                         ! Set Auu = 1 on the main diagonal
                         ! Set Auv term = 0; this term is off-diagonal for the fully assembled matrix
                         ! Set Avu term = 0 to preserve matrix symmetry (given that Auv term = 0)
                         m = indxA_2d(0,0)
                         Auu(m,i,j) = 1.d0
                         Auv(m,i,j) = 0.d0
                         Avu(m,i,j) = 0.d0

                         !TODO - Set bu above, outside iA/jA loop
                         ! Set the rhs to the prescribed velocity, forcing u = prescribed uvel for this vertex
!!                         bu(i,j) = uvel(i,j)
                            
                      else     ! not on the diagonal

                         ! Zero out non-diagonal matrix terms in the row associated with this vertex
                         m = indxA_2d(iA,jA)
                         Auu(m, i, j) = 0.d0
                         Auv(m, i, j) = 0.d0

                         ! Shift terms associated with this velocity to the rhs.
                         ! Note: The remaining operations do not change the answer, but do restore symmetry to the matrix.
                         m = indxA_2d(-iA,-jA)

                         if (umask_dirichlet(i+iA, j+jA) /= 1) then
                            ! Move (Auu term) * uvel to rhs
                            bu(i+iA, j+jA) = bu(i+iA, j+jA) - Auu(m, i+iA, j+jA) * uvel(i,j)
                            Auu(m, i+iA, j+jA) = 0.d0
                         endif

                         if (vmask_dirichlet(i+iA, j+jA) /= 1) then
                            ! Move (Avu term) * uvel to rhs
                            bv(i+iA, j+jA) = bv(i+iA, j+jA) - Avu(m, i+iA, j+jA) * uvel(i,j)
                            Avu(m, i+iA, j+jA) = 0.d0
                         endif

                      endif  ! on the diagonal

                   endif     ! i+iA and j+jA in bounds

                enddo    ! iA
                enddo    ! jA

             endif       ! umask_dirichlet

             if (vmask_dirichlet(i,j) == 1) then

                ! set the rhs to the prescribed velocity
                bv(i,j) = vvel(i,j)

                ! loop through matrix values in the rows associated with this vertex
                ! (Auv contains one row, Avv contains a second row)
                do jA = -1,1
                do iA = -1,1

                   if ( (i+iA >= 1 .and. i+iA <= nx-1)       &
                                   .and.                     &
                        (j+jA >= 1 .and. j+jA <= ny-1) ) then

                      if (iA==0 .and. jA==0) then  ! main diagonal

                         ! Set Avv = 1 on the main diagonal
                         ! Set Avu term = 0; this term is off-diagonal for the fully assembled matrix
                         ! Set Auv term = 0 to preserve matrix symmetry (given that Avu term = 0)
                         m = indxA_2d(0,0)
                         Auv(m,i,j) = 0.d0
                         Avu(m,i,j) = 0.d0
                         Avv(m,i,j) = 1.d0

                         !TODO - Set bv above, outside iA/jA loop
                         ! Set the rhs to the prescribed velocity, forcing v = prescribed vvel for this vertex
!!                         bv(i,j) = vvel(i,j)
                            
                      else     ! not on the diagonal

                         ! Zero out non-diagonal matrix terms in the rows associated with this vertex
                         m = indxA_2d(iA,jA)
                         Avu(m, i, j) = 0.d0
                         Avv(m, i, j) = 0.d0

                         ! Shift terms associated with this velocity to the rhs.
                         ! Note: The remaining operations do not change the answer, but do restore symmetry to the matrix.
                         m = indxA_2d(-iA,-jA)

                         if (umask_dirichlet(i+iA, j+jA) /= 1) then
                            ! Move (Auv term) * vvel to rhs
                            bu(i+iA, j+jA) = bu(i+iA, j+jA) - Auv(m, i+iA, j+jA) * vvel(i,j)
                            Auv(m, i+iA, j+jA) = 0.d0
                         endif

                         if (vmask_dirichlet(i+iA, j+jA) /= 1) then
                            ! Move (Avv term) * vvel to rhs
                            bv(i+iA, j+jA) = bv(i+iA, j+jA) - Avv(m, i+iA, j+jA) * vvel(i,j)
                            Avv(m, i+iA, j+jA) = 0.d0                                           
                         endif

                      endif  ! on the diagonal

                   endif     ! i+iA and j+jA in bounds

                enddo    ! iA
                enddo    ! jA

             endif       ! vmask_dirichlet

          endif          ! active_vertex
       enddo             ! i
    enddo                ! j

  end subroutine dirichlet_boundary_conditions_2d

!****************************************************************************

  subroutine compute_residual_vector_3d(nx,          ny,            &
                                        nz,          nhalo,         &
                                        active_vertex,              &
                                        Auu,         Auv,           &
                                        Avu,         Avv,           &
                                        bu,          bv,            &
                                        uvel,        vvel,          &
                                        resid_u,     resid_v,       &
                                        L2_norm,     L2_norm_relative)

    ! Compute the residual vector Ax - b and its L2 norm.
    ! This subroutine assumes that the matrix is stored in structured (x/y/z) format.

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions (for scalars)
       nz,                   &  ! number of vertical levels where velocity is computed
       nhalo                    ! number of halo layers

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F

    real(dp), dimension(nNodeNeighbors_3d,nz,nx-1,ny-1), intent(in) ::   &
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
       uvel, vvel          ! u and v components of velocity (m/yr)

    real(dp), dimension(nz,nx-1,ny-1), intent(out) ::   &
       resid_u,      & ! residual vector, divided into 2 parts
       resid_v

    real(dp), intent(out) ::    &
       L2_norm             ! L2 norm of residual vector, |Ax - b|

    real(dp), intent(out), optional ::    &
       L2_norm_relative    ! L2 norm of residual vector relative to rhs, |Ax - b| / |b|

    integer :: i, j, k, iA, jA, kA, m

    real(dp) :: L2_norm_rhs   ! L2 norm of rhs vector, |b|

    ! Compute u and v components of A*x

    resid_u(:,:,:) = 0.d0
    resid_v(:,:,:) = 0.d0

    !TODO - Replace the following by a call to matvec_multiply_structured_3d
    ! Loop over locally owned vertices

    do j = nhalo+1, ny-nhalo
    do i = nhalo+1, nx-nhalo

       if (active_vertex(i,j)) then

          do k = 1, nz

             do kA = -1,1
             do jA = -1,1
             do iA = -1,1

                if ( (k+kA >= 1 .and. k+kA <= nz)      &
                                .and.                  &
                     (i+iA >= 1 .and. i+iA <= nx-1)    &
                             .and.                     &
                     (j+jA >= 1 .and. j+jA <= ny-1) ) then

                   m = indxA_3d(iA,jA,kA)

                   resid_u(k,i,j) = resid_u(k,i,j)                     & 
                                  + Auu(m,k,i,j)*uvel(k+kA,i+iA,j+jA)  &
                                  + Auv(m,k,i,j)*vvel(k+kA,i+iA,j+jA)

                   resid_v(k,i,j) = resid_v(k,i,j)                     &
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

    ! Subtract b to get A*x - b
    ! Sum up squared L2 norm as we go

    L2_norm = 0.d0

    ! Loop over locally owned vertices

    do j = nhalo+1, ny-nhalo
    do i = nhalo+1, nx-nhalo
       if (active_vertex(i,j)) then
          do k = 1, nz
             resid_u(k,i,j) = resid_u(k,i,j) - bu(k,i,j)
             resid_v(k,i,j) = resid_v(k,i,j) - bv(k,i,j)
             L2_norm = L2_norm + resid_u(k,i,j)*resid_u(k,i,j)  &
                               + resid_v(k,i,j)*resid_v(k,i,j)
          enddo  ! k
       endif     ! active vertex
    enddo        ! i
    enddo        ! j

    ! Take global sum, then take square root
    L2_norm = parallel_reduce_sum(L2_norm)
    L2_norm = sqrt(L2_norm)

    if (verbose_residual .and. this_rank==rtest) then
       i = itest
       j = jtest
       k = ktest
       print*, 'In compute_residual_vector_3d: i, j, k =', i, j, k
       print*, 'u,  v :', uvel(k,i,j), vvel(k,i,j)
       print*, 'bu, bv:', bu(k,i,j), bv(k,i,j)
       print*, 'resid_u, resid_v:', resid_u(k,i,j), resid_v(k,i,j)
    endif

    if (present(L2_norm_relative)) then   ! compute L2_norm relative to rhs

       L2_norm_rhs = 0.d0

       do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (active_vertex(i,j)) then
             do k = 1, nz
                L2_norm_rhs = L2_norm_rhs + bu(k,i,j)*bu(k,i,j) + bv(k,i,j)*bv(k,i,j)
             enddo  ! k
          endif     ! active vertex
       enddo        ! i
       enddo        ! j

       ! Take global sum, then take square root
       L2_norm_rhs = parallel_reduce_sum(L2_norm_rhs)
       L2_norm_rhs = sqrt(L2_norm_rhs)

       if (L2_norm_rhs > 0.d0) then
          L2_norm_relative = L2_norm / L2_norm_rhs
       else
          L2_norm_relative = 0.d0
       endif

    endif

  end subroutine compute_residual_vector_3d

!****************************************************************************

  subroutine compute_residual_vector_2d(nx,          ny,            &
                                        nhalo,                      &
                                        active_vertex,              &
                                        Auu,         Auv,           &
                                        Avu,         Avv,           &
                                        bu,          bv,            &
                                        uvel,        vvel,          &
                                        resid_u,     resid_v,       &
                                        L2_norm,     L2_norm_relative)

    ! Compute the residual vector Ax - b and its L2 norm.
    ! This subroutine assumes that the matrix is stored in structured (x/y/z) format.

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions (for scalars)
       nhalo                    ! number of halo layers

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F

    real(dp), dimension(nNodeNeighbors_2d,nx-1,ny-1), intent(in) ::   &
       Auu, Auv, Avu, Avv     ! four components of assembled matrix
                              ! 1st dimension = 3 (node and its nearest neighbors in x, y and z direction)
                              ! other dimensions = (z,x,y) indices
                              !
                              !    Auu  | Auv
                              !    _____|____
                              !    Avu  | Avv
                              !         |

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       bu, bv              ! assembled load (rhs) vector, divided into 2 parts

   real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       uvel, vvel          ! u and v components of velocity (m/yr)

    real(dp), dimension(nx-1,ny-1), intent(out) ::   &
       resid_u,      & ! residual vector, divided into 2 parts
       resid_v

    real(dp), intent(out) ::    &
       L2_norm             ! L2 norm of residual vector, |Ax - b|

    real(dp), intent(out), optional ::    &
       L2_norm_relative    ! L2 norm of residual vector relative to rhs, |Ax - b| / |b|

    integer :: i, j, iA, jA, m 

    real(dp) :: L2_norm_rhs   ! L2 norm of rhs vector, |b|

    ! Compute u and v components of A*x

    resid_u(:,:) = 0.d0
    resid_v(:,:) = 0.d0

    ! Loop over locally owned vertices

    do j = nhalo+1, ny-nhalo
    do i = nhalo+1, nx-nhalo

       if (active_vertex(i,j)) then

          do jA = -1,1
             do iA = -1,1

                if ( (i+iA >= 1 .and. i+iA <= nx-1)    &
                             .and.                     &
                     (j+jA >= 1 .and. j+jA <= ny-1) ) then

                   m = indxA_2d(iA,jA)

                   resid_u(i,j) = resid_u(i,j)                     & 
                                + Auu(m,i,j)*uvel(i+iA,j+jA)  &
                                + Auv(m,i,j)*vvel(i+iA,j+jA)

                   resid_v(i,j) = resid_v(i,j)                     &
                                + Avu(m,i,j)*uvel(i+iA,j+jA)  &
                                + Avv(m,i,j)*vvel(i+iA,j+jA)

                endif   ! in bounds

             enddo   ! iA
          enddo      ! jA

       endif   ! active_vertex

    enddo   ! i
    enddo   ! j

    ! Subtract b to get A*x - b
    ! Sum up squared L2 norm as we go

    L2_norm = 0.d0

    ! Loop over locally owned vertices

    do j = nhalo+1, ny-nhalo
    do i = nhalo+1, nx-nhalo
       if (active_vertex(i,j)) then
          resid_u(i,j) = resid_u(i,j) - bu(i,j)
          resid_v(i,j) = resid_v(i,j) - bv(i,j)
          L2_norm = L2_norm + resid_u(i,j)*resid_u(i,j)  &
                            + resid_v(i,j)*resid_v(i,j)
       endif     ! active vertex
    enddo        ! i
    enddo        ! j

    ! Take global sum, then take square root

    L2_norm = parallel_reduce_sum(L2_norm)
    L2_norm = sqrt(L2_norm)

    if (verbose_residual .and. this_rank==rtest) then
       i = itest
       j = jtest
       print*, 'In compute_residual_vector_2d: i, j =', i, j
       print*, 'u,  v :', uvel(i,j), vvel(i,j)
       print*, 'bu, bv:', bu(i,j), bv(i,j)
       print*, 'resid_u, resid_v:', resid_u(i,j), resid_v(i,j)
       print*, ' '
       print*, 'maxval/minval(resid_u) =', maxval(resid_u), minval(resid_u)
       print*, 'maxval/minval(resid_v) =', maxval(resid_v), minval(resid_v)
    endif

    if (present(L2_norm_relative)) then   ! compute L2_norm relative to rhs

       L2_norm_rhs = 0.d0

       do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (active_vertex(i,j)) then
             L2_norm_rhs = L2_norm_rhs + bu(i,j)*bu(i,j) + bv(i,j)*bv(i,j)
          endif     ! active vertex
       enddo        ! i
       enddo        ! j

       ! Take global sum, then take square root
       L2_norm_rhs = parallel_reduce_sum(L2_norm_rhs)
       L2_norm_rhs = sqrt(L2_norm_rhs)

       if (L2_norm_rhs > 0.d0) then
          L2_norm_relative = L2_norm / L2_norm_rhs
       else
          L2_norm_relative = 0.d0
       endif

    endif

  end subroutine compute_residual_vector_2d

!****************************************************************************

  subroutine compute_residual_velocity_3d(nhalo,  whichresid, &
                                          uvel,   vvel,        &
                                          usav,   vsav,        &
                                          resid_velo)

    integer, intent(in) ::   &
       nhalo,           & ! number of layers of halo cells
       whichresid         ! option for method to use when calculating residual

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
    !TODO - Remove some of these velocity residual methods?  They are rarely if ever used.

    ! options for residual calculation method, as specified in configuration file
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

       !TODO - Need to convert the mean residual to a global value.
       !       (Or simply remove this case, which is rarely if ever used)
       call not_parallel(__FILE__, __LINE__)

   case default    ! max speed difference, including basal speeds
                   ! (case HO_RESID_MAXU or HO_RESID_L2NORM or HO_RESID_L2NORM_RELATIVE)

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

  end subroutine compute_residual_velocity_3d

!****************************************************************************

  subroutine compute_residual_velocity_2d(nhalo,  whichresid, &
                                          uvel,   vvel,        &
                                          usav,   vsav,        &
                                          resid_velo)

    integer, intent(in) ::   &
       nhalo,           & ! number of layers of halo cells
       whichresid         ! option for method to use when calculating residual

    real(dp), dimension(:,:), intent(in) ::  &
       uvel, vvel,      & ! current guess for velocity
       usav, vsav         ! previous guess for velocity

    real(dp), intent(out) ::    &
       resid_velo         ! quantity related to velocity convergence


    integer ::   &
       imaxdiff, jmaxdiff   ! location of maximum speed difference
                            ! currently computed but not used
 
    integer :: i, j, count

    real(dp) ::   &
       speed,      &   ! current guess for ice speed
       oldspeed,   &   ! previous guess for ice speed
       diffspeed       ! abs(speed-oldspeed)


    ! Compute a residual quantity based on convergence of the velocity field.

    ! options for residual calculation method, as specified in configuration file
    ! case(0): use max of abs( vel_old - vel ) / vel )
    ! case(1): use max of abs( vel_old - vel ) / vel ) but ignore basal vels
    ! case(2): use mean of abs( vel_old - vel ) / vel )
    ! case(3): use max of abs( vel_old - vel ) / vel ) (in addition to L2 norm)
    
    resid_velo = 0.d0
    imaxdiff = 0
    jmaxdiff = 0

    select case (whichresid)

    case(HO_RESID_MAXU_NO_UBAS)   ! max speed difference, excluding the bed

       ! Loop over locally owned vertices

       do j = 1+nhalo, size(uvel,2)-nhalo
          do i = 1+nhalo, size(uvel,1)-nhalo
             speed = sqrt(uvel(i,j)**2 + vvel(i,j)**2)
             if (speed /= 0.d0) then
                oldspeed = sqrt(usav(i,j)**2 + vsav(i,j)**2)
                diffspeed = abs((oldspeed - speed)/speed)
                if (diffspeed > resid_velo) then
                   resid_velo = diffspeed
                   imaxdiff = i
                   jmaxdiff = j
                endif
             endif
          enddo
       enddo
       
       ! take global max
       resid_velo = parallel_reduce_max(resid_velo)

    case(HO_RESID_MEANU)   ! mean relative speed difference

       count = 0

       ! Loop over locally owned vertices

       do j = 1+nhalo, size(uvel,2)-nhalo
          do i = 1+nhalo, size(uvel,1)-nhalo
             speed = sqrt(uvel(i,j)**2 + vvel(i,j)**2)
             if (speed /= 0.d0) then
                count = count+1
                oldspeed = sqrt(usav(i,j)**2 + vsav(i,j)**2)
                diffspeed = abs((oldspeed - speed)/speed)                
                resid_velo = resid_velo + diffspeed
             endif
          enddo
       enddo

       if (count > 0) resid_velo = resid_velo / count

       !TODO - Need to convert the mean residual to a global value.
       !       (Or simply remove this case, which is rarely if ever used)
       call not_parallel(__FILE__, __LINE__)

   case default    ! max speed difference, including basal speeds
                   ! (case HO_RESID_MAXU or HO_RESID_L2NORM)

       ! Loop over locally owned vertices

       do j = 1+nhalo, size(uvel,2)-nhalo
          do i = 1+nhalo, size(uvel,1)-nhalo
             speed = sqrt(uvel(i,j)**2 + vvel(i,j)**2)
             if (speed /= 0.d0) then
                oldspeed = sqrt(usav(i,j)**2 + vsav(i,j)**2)
                diffspeed = abs((oldspeed - speed)/speed)
                if (diffspeed > resid_velo) then
                   resid_velo = diffspeed
                   imaxdiff = i
                   jmaxdiff = j
                endif
             endif
          enddo
       enddo

       resid_velo = parallel_reduce_max(resid_velo)
       
    end select

  end subroutine compute_residual_velocity_2d

!****************************************************************************

  subroutine count_nonzeros_3d(nx,      ny,     &
                               nz,      nhalo,  &
                               Auu,     Auv,    &
                               Avu,     Avv,    &
                               active_vertex,   & 
                               nNonzeros)

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx,  ny,              &    ! number of grid cells in each direction
       nz,                   &    ! number of vertical levels where velocity is computed
       nhalo                      ! number of halo layers

    real(dp), dimension(nNodeNeighbors_3d,nz,nx-1,ny-1), intent(in) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv                                    

    logical, dimension(nx-1,ny-1), intent(in) :: &
       active_vertex      ! true for vertices of active cells

    integer, intent(out) ::   &
       nNonzeros          ! number of nonzero matrix elements

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j, k, iA, jA, kA, m

    nNonzeros = 0
    do j = nhalo+1, ny-nhalo  ! loop over locally owned vertices
       do i = nhalo+1, nx-nhalo
          if (active_vertex(i,j)) then
             do k = 1, nz
                do kA = -1, 1
                do jA = -1, 1
                do iA = -1, 1
                   m = indxA_3d(iA,jA,kA)
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
    enddo           ! j

    nNonzeros = parallel_reduce_sum(nNonzeros)

  end subroutine count_nonzeros_3d

!****************************************************************************

  subroutine count_nonzeros_2d(nx,      ny,     &
                               nhalo,           &
                               Auu,     Auv,    &
                               Avu,     Avv,    &
                               active_vertex,   & 
                               nNonzeros)

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx,  ny,              &    ! number of grid cells in each direction
       nhalo                      ! number of halo layers

    real(dp), dimension(nNodeNeighbors_2d,nx-1,ny-1), intent(in) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv                                    

    logical, dimension(nx-1,ny-1), intent(in) :: &
       active_vertex      ! true for vertices of active cells

    integer, intent(out) ::   &
       nNonzeros          ! number of nonzero matrix elements

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j, iA, jA, m

    nNonzeros = 0
    do j = nhalo+1, ny-nhalo  ! loop over locally owned vertices
       do i = nhalo+1, nx-nhalo
          if (active_vertex(i,j)) then
             do jA = -1, 1
                do iA = -1, 1
                   m = indxA_2d(iA,jA)
                   if (Auu(m,i,j) /= 0.d0) nNonzeros = nNonzeros + 1
                   if (Auv(m,i,j) /= 0.d0) nNonzeros = nNonzeros + 1
                   if (Avu(m,i,j) /= 0.d0) nNonzeros = nNonzeros + 1
                   if (Avv(m,i,j) /= 0.d0) nNonzeros = nNonzeros + 1
                enddo 
             enddo
          endif     ! active_vertex
       enddo        ! i
    enddo           ! j

    nNonzeros = parallel_reduce_sum(nNonzeros)

  end subroutine count_nonzeros_2d

!****************************************************************************

  subroutine check_symmetry_element_matrix(nNodesPerElement,  &
                                           Kuu, Kuv, Kvu, Kvv)

    !------------------------------------------------------------------
    ! Check that the element stiffness matrix is symmetric.
    ! This is true provided that (1) Kuu = (Kuu)^T
    !                            (2) Kvv = (Kvv)^T
    !                            (3) Kuv = (Kvu)^T
    ! This subroutine works for either 2D or 3D elements.
    ! A symmetry check should not be needed for production runs with a well-tested code,
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

    ! make sure Kuu = (Kuu)^T

    do j = 1, nNodesPerElement
       do i = j, nNodesPerElement
          if (abs(Kuu(i,j) - Kuu(j,i)) > eps10) then
             print*, 'Kuu is not symmetric'
             print*, 'i, j, Kuu(i,j), Kuu(j,i):', i, j, Kuu(i,j), Kuu(j,i)
             stop
          endif    
       enddo
    enddo

    ! check that Kvv = (Kvv)^T

    do j = 1, nNodesPerElement
       do i = j, nNodesPerElement
          if (abs(Kvv(i,j) - Kvv(j,i)) > eps10) then
             print*, 'Kvv is not symmetric'
             print*, 'i, j, Kvv(i,j), Kvv(j,i):', i, j, Kvv(i,j), Kvv(j,i)
             stop
          endif    
       enddo
    enddo

    ! Check that Kuv = (Kvu)^T

    do j = 1, nNodesPerElement
       do i = 1, nNodesPerElement
          if (abs(Kuv(i,j) - Kvu(j,i)) > eps10) then
             print*, 'Kuv /= (Kvu)^T'
             print*, 'i, j, Kuv(i,j), Kvu(j,i):', i, j, Kuv(i,j), Kvu(j,i)
             stop
          endif    
       enddo
    enddo

  end subroutine check_symmetry_element_matrix

!****************************************************************************

  subroutine check_symmetry_assembled_matrix_3d(nx,  ny,  nz, nhalo,   &
                                                active_vertex,         &
                                                Auu, Auv, Avu, Avv)

    !------------------------------------------------------------------
    ! Check that the assembled stiffness matrix is symmetric.
    ! This is true provided that (1) Auu = (Auu)^T
    !                            (2) Avv = (Avv)^T
    !                            (3) Auv = (Avu)^T
    ! The A matrices are assembled in a dense fashion to save storage
    !  and preserve the i/j/k structure of the grid.
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

    real(dp), dimension(nNodeNeighbors_3d,nz,nx-1,ny-1), intent(inout) ::   &
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
    ! The maximum departure from symmetry is set to be a small fraction 
    !  of the diagonal entry for the row.
    ! If the departure from symmetry is larger than this, then the model prints a warning 
    !  and/or aborts.

    maxdiff = 0.d0

    ! Loop over locally owned vertices.
    ! Each active vertex is associate with 2*nz matrix rows belonging to this processor.
    ! Locally owned vertices are (nhalo+1:ny-nhalo, nhalo+1:nx-nhalo)

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (active_vertex(i,j)) then
             do k = 1, nz

                ! Check Auu and Auv for symmetry

                m = indxA_3d(0,0,0)
                diag_entry = Auu(m,k,i,j)

                do jA = -1, 1
                do iA = -1, 1
                do kA = -1, 1

                   if (k+kA >= 1 .and. k+kA <=nz) then  ! to keep k index in bounds

                      m =  indxA_3d( iA, jA, kA)
                      mm = indxA_3d(-iA,-jA,-kA)

                      ! Check that Auu = Auu^T

                      val1 = Auu( m, k,    i,    j   )   ! value of Auu(row,col)
                      val2 = Auu(mm, k+kA, i+iA, j+jA)   ! value of Auu(col,row)

                      if (val2 /= val1) then
                          
                         if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)

                         ! if difference is small, then fix the asymmetry by averaging values
                         ! else print a warning and abort

                         if ( abs(val2-val1) < eps08*abs(diag_entry) ) then
                            avg_val = 0.5d0 * (val1 + val2)
                            Auu( m, k,   i,   j   ) = avg_val
                            Auu(mm, k+kA,i+iA,j+jA) = avg_val
                         else
                            print*, 'WARNING: Auu is not symmetric: i, j, k, iA, jA, kA =', i, j, k, iA, jA, kA
                            print*, 'Auu(row,col), Auu(col,row), diff/diag:', val1, val2, (val2 - val1)/diag_entry
!!                            stop
                         endif

                      endif   ! val2 /= val1
                
                      ! Check that Auv = (Avu)^T

                      val1 = Auv( m, k,    i,    j)      ! value of Auv(row,col)
                      val2 = Avu(mm, k+kA, i+iA, j+jA)   ! value of Avu(col,row)

                      if (val2 /= val1) then

                         if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)

                         ! if difference is small, then fix the asymmetry by averaging values
                         ! else print a warning and abort

                         if ( abs(val2-val1) < eps08*abs(diag_entry) ) then
                            avg_val = 0.5d0 * (val1 + val2)
                            Auv( m, k,   i,   j   ) = avg_val
                            Avu(mm, k+kA,i+iA,j+jA) = avg_val
                         else
                            print*, 'WARNING: Auv is not equal to (Avu)^T, i, j, k, iA, jA, kA =', i, j, k, iA, jA, kA
                            print*, 'Auv(row,col), Avu(col,row), diff/diag:', val1, val2, (val2 - val1)/diag_entry
!!                            stop
                         endif

                      endif  ! val2 /= val1

                   endif     ! k+kA in bounds
            
                enddo        ! kA
                enddo        ! iA
                enddo        ! jA

                ! Now check Avu and Avv

                m = indxA_3d(0,0,0)
                diag_entry = Avv(m,k,i,j)

                ! check that Avv = (Avv)^T

                do jA = -1, 1
                do iA = -1, 1
                do kA = -1, 1

                   if (k+kA >= 1 .and. k+kA <=nz) then  ! to keep k index in bounds

                      m  = indxA_3d( iA, jA, kA)
                      mm = indxA_3d(-iA,-jA,-kA)

                      val1 = Avv( m, k,    i,    j)      ! value of Avv(row,col)
                      val2 = Avv(mm, k+kA, i+iA, j+jA)   ! value of Avv(col,row)

                      if (val2 /= val1) then
                          
                         if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)

                         ! if difference is small, then fix the asymmetry by averaging values
                         ! else print a warning and abort

                         if ( abs(val2-val1) < eps08*abs(diag_entry) ) then
                            avg_val = 0.5d0 * (val1 + val2)
                            Avv( m, k,   i,   j   ) = avg_val
                            Avv(mm, k+kA,i+iA,j+jA) = avg_val
                         else
                            print*, 'WARNING: Avv is not symmetric: i, j, k, iA, jA, kA =', i, j, k, iA, jA, kA
                            print*, 'Avv(row,col), Avv(col,row), diff/diag:', val1, val2, (val2 - val1)/diag_entry
!!                            stop
                         endif

                      endif   ! val2 /= val1

                      ! Check that Avu = (Auv)^T

                      val1 = Avu( m, k,    i,    j)      ! value of Avu(row,col)
                      val2 = Auv(mm, k+kA, i+iA, j+jA)   ! value of Auv(col,row)

                      if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)

                      if (val2 /= val1) then

                         ! if difference is small, then fix the asymmetry by averaging values
                         ! else print a warning and abort

                         if ( abs(val2-val1) < eps08*abs(diag_entry) ) then
                            avg_val = 0.5d0 * (val1 + val2)
                            Avu( m, k,   i,   j   ) = avg_val
                            Auv(mm, k+kA,i+iA,j+jA) = avg_val
                         else
                            print*, 'WARNING: Avu is not equal to (Auv)^T, i, j, k, iA, jA, kA =', i, j, k, iA, jA, kA
                            print*, 'Avu(row,col), Auv(col,row), diff/diag:', val1, val2, (val2 - val1)/diag_entry
!!                            stop
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

    if (verbose_matrix) maxdiff = parallel_reduce_max(maxdiff)

    if (verbose_matrix .and. main_task) then
       print*, ' '
       print*, 'Max difference from symmetry =', maxdiff
    endif

  end subroutine check_symmetry_assembled_matrix_3d

!****************************************************************************

  subroutine check_symmetry_assembled_matrix_2d(nx,  ny, nhalo,     &
                                                active_vertex,      &
                                                Auu, Auv, Avu, Avv)

    !------------------------------------------------------------------
    ! Check that the assembled stiffness matrix is symmetric.
    ! This is true provided that (1) Auu = (Auu)^T
    !                            (2) Avv = (Avv)^T
    !                            (3) Auv = (Avu)^T
    ! The A matrices are assembled in a dense fashion to save storage
    !  and preserve the i/j/k structure of the grid.
    !
    ! There can be small differences from perfect symmetry due to roundoff error.
    ! These differences are fixed provided they are small enough.
    !------------------------------------------------------------------    

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nhalo                    ! number of halo layers

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex            ! T for columns (i,j) where velocity is computed, else F

    real(dp), dimension(nNodeNeighbors_2d,nx-1,ny-1), intent(inout) ::   &
       Auu, Auv, Avu, Avv       ! components of assembled stiffness matrix
                                !
                                !    Auu  | Auv
                                !    _____|____          
                                !         |
                                !    Avu  | Avv                                    

    integer :: i, j, iA, jA, m, mm

    real(dp) :: val1, val2          ! values of matrix coefficients

    real(dp) :: maxdiff, diag_entry, avg_val

    ! Check matrix for symmetry

    ! Here we correct for small differences from symmetry due to roundoff error.
    ! The maximum departure from symmetry is set to be a small fraction
    !  of the diagonal entry for the row.
    ! If the departure from symmetry is larger than this, then the model prints a warning 
    !  and/or aborts.

    maxdiff = 0.d0

    ! Loop over locally owned vertices.
    ! Each active vertex is associate with 2*nz matrix rows belonging to this processor.
    ! Locally owned vertices are (nhalo+1:ny-nhalo, nhalo+1:nx-nhalo)

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (active_vertex(i,j)) then

             ! Check Auu and Auv for symmetry

             m = indxA_2d(0,0)
             diag_entry = Auu(m,i,j)

             do jA = -1, 1
             do iA = -1, 1

                m =  indxA_2d( iA, jA)
                mm = indxA_2d(-iA,-jA)

                ! Check that Auu = Auu^T

                val1 = Auu( m, i,    j   )   ! value of Auu(row,col)
                val2 = Auu(mm, i+iA, j+jA)   ! value of Auu(col,row)

                if (val2 /= val1) then
                          
                   if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)

                   ! if difference is small, then fix the asymmetry by averaging values
                   ! else print a warning and abort

                   if ( abs(val2-val1) < eps08*abs(diag_entry) ) then
                      avg_val = 0.5d0 * (val1 + val2)
                      Auu( m, i,   j   ) = avg_val
                      Auu(mm, i+iA,j+jA) = avg_val
                   else
                      print*, 'WARNING: Auu is not symmetric: this_rank, i, j, iA, jA =', this_rank, i, j, iA, jA
                      print*, 'Auu(row,col), Auu(col,row), diff/diag:', val1, val2, (val2 - val1)/diag_entry
!!                      stop
                   endif

                endif   ! val2 /= val1
                
                ! Check that Auv = (Avu)^T

                val1 = Auv( m,    i,    j)   ! value of Auv(row,col)
                val2 = Avu(mm, i+iA, j+jA)   ! value of Avu(col,row)

                if (val2 /= val1) then

                   if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)

                   ! if difference is small, then fix the asymmetry by averaging values
                   ! else print a warning and abort

                   if ( abs(val2-val1) < eps08*abs(diag_entry) ) then
                      avg_val = 0.5d0 * (val1 + val2)
                      Auv( m,   i,   j) = avg_val
                      Avu(mm,i+iA,j+jA) = avg_val
                   else
                      print*, 'WARNING: Auv is not equal to (Avu)^T, i, j, iA, jA =', i, j, iA, jA
                      print*, 'Auv(row,col), Avu(col,row), diff/diag:', val1, val2, (val2 - val1)/diag_entry
!!                      stop
                   endif

                endif  ! val2 /= val1

             enddo        ! iA
             enddo        ! jA

             ! Now check Avu and Avv

             m = indxA_2d(0,0)
             diag_entry = Avv(m,i,j)

             ! check that Avv = (Avv)^T

             do jA = -1, 1
             do iA = -1, 1

                m  = indxA_2d( iA, jA)
                mm = indxA_2d(-iA,-jA)

                val1 = Avv( m,    i,    j)   ! value of Avv(row,col)
                val2 = Avv(mm, i+iA, j+jA)   ! value of Avv(col,row)

                if (val2 /= val1) then
                          
                   if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)

                   ! if difference is small, then fix the asymmetry by averaging values
                   ! else print a warning and abort

                   if ( abs(val2-val1) < eps08*abs(diag_entry) ) then
                      avg_val = 0.5d0 * (val1 + val2)
                      Avv( m,   i,   j) = avg_val
                      Avv(mm,i+iA,j+jA) = avg_val
                   else
                      print*, 'WARNING: Avv is not symmetric: i, j, iA, jA =', i, j, iA, jA
                      print*, 'Avv(row,col), Avv(col,row), diff/diag:', val1, val2, (val2 - val1)/diag_entry
!!                      stop
                   endif

                endif   ! val2 /= val1

                ! Check that Avu = (Auv)^T

                val1 = Avu( m,    i,    j)   ! value of Avu(row,col)
                val2 = Auv(mm, i+iA, j+jA)   ! value of Auv(col,row)

                if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)

                if (val2 /= val1) then

                   ! if difference is small, then fix the asymmetry by averaging values
                   ! else print a warning and abort

                   if ( abs(val2-val1) < eps08*abs(diag_entry) ) then
                      avg_val = 0.5d0 * (val1 + val2)
                      Avu( m,   i,   j) = avg_val
                      Auv(mm,i+iA,j+jA) = avg_val
                   else
                      print*, 'WARNING: Avu is not equal to (Auv)^T, i, j, iA, jA =', i, j, iA, jA
                      print*, 'Avu(row,col), Auv(col,row), diff/diag:', val1, val2, (val2 - val1)/diag_entry
!!                      stop
                   endif

                endif  ! val2 /= val1

             enddo     ! iA
             enddo     ! jA

          endif        ! active_vertex
       enddo           ! i
    enddo              ! j

    if (verbose_matrix) maxdiff = parallel_reduce_max(maxdiff)

    if (verbose_matrix .and. main_task) then
       print*, ' '
       print*, 'Max difference from symmetry =', maxdiff
    endif

  end subroutine check_symmetry_assembled_matrix_2d

!****************************************************************************

  subroutine write_matrix_elements_3d(nx,    ny,   nz,     &
                                      nNodesSolve, nodeID, &
                                      iNodeIndex,  jNodeIndex,  &
                                      kNodeIndex,          &
                                      Auu,         Auv,    &
                                      Avu,         Avv,    &
                                      bu,          bv)

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz,                   &  ! number of vertical levels at which velocity is computed
       nNodesSolve              ! number of nodes where we solve for velocity

    integer, dimension(nz,nx-1,ny-1), intent(in) ::  &
       nodeID             ! ID for each node

    integer, dimension(:), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of active nodes

    real(dp), dimension(nNodeNeighbors_3d,nz,nx-1,ny-1), intent(in) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv           ! 1st dimension = node and its nearest neighbors in x, y and z direction 
                          ! other dimensions = (k,i,j) indices

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       bu, bv             ! assembled load (rhs) vector, divided into 2 parts

    ! Local variables

    integer :: rowA, colA
    integer :: i, j, k, m, iA, jA, kA

    real(dp), dimension(nNodesSolve, nNodesSolve) ::   &
       Auu_val, Auv_val, Avu_val, Avv_val   ! dense matrices

    real(dp), dimension(nNodesSolve) :: nonzeros

    Auu_val(:,:) = 0.d0
    Auv_val(:,:) = 0.d0
    Avu_val(:,:) = 0.d0
    Avv_val(:,:) = 0.d0

    do rowA = 1, nNodesSolve

       i = iNodeIndex(rowA)
       j = jNodeIndex(rowA)
       k = kNodeIndex(rowA)

       do kA = -1, 1
       do jA = -1, 1
       do iA = -1, 1

          if ( (k+kA >= 1 .and. k+kA <= nz)         &
                          .and.                     &
               (i+iA >= 1 .and. i+iA <= nx-1)       &
                          .and.                     &
               (j+jA >= 1 .and. j+jA <= ny-1) ) then

             colA = nodeID(k+kA, i+iA, j+jA)   ! ID for neighboring node
             m = indxA_3d(iA,jA,kA)

             if (colA > 0) then 
                Auu_val(rowA, colA) = Auu(m,k,i,j)
                Auv_val(rowA, colA) = Auv(m,k,i,j)
                Avu_val(rowA, colA) = Avu(m,k,i,j)
                Avv_val(rowA, colA) = Avv(m,k,i,j)
             endif

          endif     ! i+iA, j+jA, and k+kA in bounds

       enddo        ! kA
       enddo        ! iA
       enddo        ! jA

    enddo           ! rowA 

    !WHL - bug check
    print*, ' '
    print*, 'nonzeros per row:'
    do rowA = 1, nNodesSolve
       nonzeros(rowA) = 0
       do colA = 1, nNodesSolve
          if (abs(Auu_val(rowA,colA)) > 1.d-11) then
             nonzeros(rowA) = nonzeros(rowA) + 1
          endif
       enddo
!       print*, rowA, nonzeros(rowA)
    enddo

    print*, 'Write matrix elements to file, label =', matrix_label

    ! Write matrices to file (one line of file corresponding to each row of matrix)

    open(unit=10, file='Auu.'//matrix_label, status='unknown')
    open(unit=11, file='Auv.'//matrix_label, status='unknown')
    open(unit=12, file='Avu.'//matrix_label, status='unknown')
    open(unit=13, file='Avv.'//matrix_label, status='unknown')

    do rowA = 1, nNodesSolve
       write(10,'(i6)',advance='no') rowA
       write(11,'(i6)',advance='no') rowA
       write(12,'(i6)',advance='no') rowA
       write(13,'(i6)',advance='no') rowA
       do colA = 1, nNodesSolve
          write(10,'(e16.8)',advance='no') Auu_val(rowA,colA)
          write(11,'(e16.8)',advance='no') Auv_val(rowA,colA)
          write(12,'(e16.8)',advance='no') Avu_val(rowA,colA)
          write(13,'(e16.8)',advance='no') Avv_val(rowA,colA)
       enddo
       write(10,*) ' '
       write(11,*) ' '
       write(12,*) ' '
       write(13,*) ' '
    enddo

    close(10)
    close(11)
    close(12)
    close(13)

    print*, 'Done writing matrix elements'

    ! write load vectors to file
    open(unit=14, file='bu.'//matrix_label, status='unknown')
    open(unit=15, file='bv.'//matrix_label, status='unknown')
    do rowA = 1, nNodesSolve
       i = iNodeIndex(rowA)
       j = jNodeIndex(rowA)
       k = kNodeIndex(rowA)
       write(14,'(i6, e16.8)') rowA, bu(k,i,j)
       write(15,'(i6, e16.8)') rowA, bv(k,i,j)
    enddo
    close(14)
    close(15)

  end subroutine write_matrix_elements_3d
  
!****************************************************************************

  subroutine write_matrix_elements_2d(nx,             ny,            &
                                      nVerticesSolve, vertexID,      &
                                      iVertexIndex,   jVertexIndex,  &
                                      Auu,            Auv,           &
                                      Avu,            Avv,           &
                                      bu,             bv)

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nVerticesSolve           ! number of vertices where we solve for velocity

    integer, dimension(nx-1,ny-1), intent(in) ::  &
       vertexID             ! ID for each vertex

    integer, dimension(:), intent(in) ::   &
       iVertexIndex, jVertexIndex   ! i and j indices of active vertices

    real(dp), dimension(nNodeNeighbors_2d,nx-1,ny-1), intent(in) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv           ! 1st dimension = vertex and its nearest neighbors in x and y direction 
                          ! other dimensions = (i,j) indices

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       bu, bv             ! assembled load (rhs) vector, divided into 2 parts

    ! Local variables

    integer :: rowA, colA
    integer :: i, j, m, iA, jA

    real(dp), dimension(nVerticesSolve, nVerticesSolve) ::   &
       Auu_val, Auv_val, Avu_val, Avv_val   ! dense matrices

    real(dp), dimension(nVerticesSolve) :: nonzeros

    Auu_val(:,:) = 0.d0
    Auv_val(:,:) = 0.d0
    Avu_val(:,:) = 0.d0
    Avv_val(:,:) = 0.d0

    do rowA = 1, nVerticesSolve

       i = iVertexIndex(rowA)
       j = jVertexIndex(rowA)
       do jA = -1, 1
       do iA = -1, 1

          if ( (i+iA >= 1 .and. i+iA <= nx-1)       &
                          .and.                     &
               (j+jA >= 1 .and. j+jA <= ny-1) ) then

             colA = vertexID(i+iA, j+jA)   ! ID for neighboring vertex
             m = indxA_2d(iA,jA)

             if (colA > 0) then 
                Auu_val(rowA, colA) = Auu(m,i,j)
                Auv_val(rowA, colA) = Auv(m,i,j)
                Avu_val(rowA, colA) = Avu(m,i,j)
                Avv_val(rowA, colA) = Avv(m,i,j)
             endif

          endif     ! i+iA and j+jA in bounds

       enddo        ! iA
       enddo        ! jA

    enddo           ! rowA 

    !WHL - bug check
    print*, ' '
    print*, 'nonzeros per row:'
    do rowA = 1, nVerticesSolve
       nonzeros(rowA) = 0
       do colA = 1, nVerticesSolve
          if (abs(Auu_val(rowA,colA)) > 1.d-11) then
             nonzeros(rowA) = nonzeros(rowA) + 1
          endif
       enddo
!       print*, rowA, nonzeros(rowA)
    enddo

    print*, 'Write matrix elements to file, label =', matrix_label

    ! Write matrices to file (one line of file corresponding to each row of matrix)

    open(unit=10, file='Auu.'//matrix_label, status='unknown')
    open(unit=11, file='Auv.'//matrix_label, status='unknown')
    open(unit=12, file='Avu.'//matrix_label, status='unknown')
    open(unit=13, file='Avv.'//matrix_label, status='unknown')

    do rowA = 1, nVerticesSolve
       write(10,'(i6)',advance='no') rowA
       write(11,'(i6)',advance='no') rowA
       write(12,'(i6)',advance='no') rowA
       write(13,'(i6)',advance='no') rowA
       do colA = 1, nVerticesSolve
          write(10,'(e16.8)',advance='no') Auu_val(rowA,colA)
          write(11,'(e16.8)',advance='no') Auv_val(rowA,colA)
          write(12,'(e16.8)',advance='no') Avu_val(rowA,colA)
          write(13,'(e16.8)',advance='no') Avv_val(rowA,colA)
       enddo
       write(10,*) ' '
       write(11,*) ' '
       write(12,*) ' '
       write(13,*) ' '
    enddo

    close(10)
    close(11)
    close(12)
    close(13)

    print*, 'Done writing matrix elements'

    ! write load vectors to file
    open(unit=14, file='bu.'//matrix_label, status='unknown')
    open(unit=15, file='bv.'//matrix_label, status='unknown')
    do rowA = 1, nVerticesSolve
       i = iVertexIndex(rowA)
       j = jVertexIndex(rowA)
       write(14,'(i6, e16.8)') rowA, bu(i,j)
       write(15,'(i6, e16.8)') rowA, bv(i,j)
    enddo
    close(14)
    close(15)

  end subroutine write_matrix_elements_2d

!****************************************************************************

  subroutine compress_3d_to_2d(nx,        ny,      nz,  &
                               Auu,       Auv,          &
                               Avu,       Avv,          &
                               bu,        bv,           &
                               Auu_2d,    Auv_2d,       &
                               Avu_2d,    Avv_2d,       &
                               bu_2d,     bv_2d)

    !----------------------------------------------------------------
    ! Form the 2D matrix and rhs by combining terms from the 3D matrix and rhs.
    ! This combination is based on the assumption of no vertical shear;
    !  i.e., uvel and vvel have the same value at each level in a given column.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz                       ! number of vertical levels where velocity is computed

    real(dp), dimension(nNodeNeighbors_3d,nz,nx-1,ny-1), intent(in) ::  &
       Auu, Auv,    &     ! assembled 3D stiffness matrix, divided into 4 parts
       Avu, Avv           
                          
    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       bu, bv             ! assembled 3D rhs vector, divided into 2 parts
                          
    real(dp), dimension(nNodeNeighbors_2d,nx-1,ny-1), intent(out) ::  &
       Auu_2d, Auv_2d,   &! assembled 2D (SSA) stiffness matrix, divided into 4 parts
       Avu_2d, Avv_2d           
                          
    real(dp), dimension(nx-1,ny-1), intent(out) ::  &
       bu_2d, bv_2d       ! assembled 2D (SSA) rhs vector, divided into 2 parts
                          
    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j, k, iA, jA, kA, m, m2

    ! Initialize 2D matrix and rhs

    Auu_2d(:,:,:) = 0.d0
    Auv_2d(:,:,:) = 0.d0
    Avu_2d(:,:,:) = 0.d0
    Avv_2d(:,:,:) = 0.d0
    bu_2d(:,:) = 0.d0
    bv_2d(:,:) = 0.d0

    ! Form 2D matrix and rhs

    do j = 1, ny-1
       do i = 1, nx-1
          do k = 1, nz

             ! matrix
             do kA = -1,1
                do jA = -1,1
                   do iA = -1,1
                      m = indxA_3d(iA,jA,kA)
                      m2 = indxA_2d(iA,jA)
                      Auu_2d(m2,i,j) = Auu_2d(m2,i,j) + Auu(m,k,i,j)
                      Auv_2d(m2,i,j) = Auv_2d(m2,i,j) + Auv(m,k,i,j)
                      Avu_2d(m2,i,j) = Avu_2d(m2,i,j) + Avu(m,k,i,j)
                      Avv_2d(m2,i,j) = Avv_2d(m2,i,j) + Avv(m,k,i,j)
                   enddo   ! iA
                enddo      ! jA
             enddo         ! kA

             ! rhs
             bu_2d(i,j) = bu_2d(i,j) + bu(k,i,j)
             bv_2d(i,j) = bv_2d(i,j) + bv(k,i,j)

          enddo            ! k
       enddo               ! i
    enddo                  ! j

  end subroutine compress_3d_to_2d

!****************************************************************************

  end module glissade_velo_higher

!****************************************************************************
