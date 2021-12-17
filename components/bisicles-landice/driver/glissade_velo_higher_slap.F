!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_velo_higher_slap.F90 - part of the Community Ice Sheet Model (CISM)  
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
! This module contains subroutines called from glissade_velo_higher.F90
! and used to process data before and after linking to SLAP solver routines.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  module glissade_velo_higher_slap

    use glimmer_global, only: dp
    use glimmer_sparse_type
    use glimmer_sparse, only: sparse_easy_solve

    implicit none

    private
    public :: slap_preprocess_3d,  slap_preprocess_2d,   &
              slap_postprocess_3d, slap_postprocess_2d,  &
              slap_compute_residual_vector, slap_solve_test_matrix

  contains

!****************************************************************************

  subroutine slap_preprocess_3d(nx,           ny,         nz,   &
                                nNodesSolve,  NodeID,      &
                                iNodeIndex,   jNodeIndex,  &
                                kNodeIndex,   indxA,       &
                                Auu,          Auv,         &
                                Avu,          Avv,         &  
                                bu,           bv,          &
                                uvel,         vvel,        &
                                matrix_order,              &
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
       NodeID             ! local ID for each node

    integer, dimension(:), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of active nodes

    integer, dimension(-1:1,-1:1,-1:1), intent(in) :: &
       indxA                 ! maps relative (x,y,z) coordinates to an index between 1 and 27                                                 
                             ! index order is (i,j,k)

    real(dp), dimension(27,nz,nx-1,ny-1), intent(in) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv           ! 1st dimension = node and its nearest neighbors in x, y and z direction 
                          ! other dimensions = (k,i,j) indices

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       bu, bv             ! assembled load (rhs) vector, divided into 2 parts
                          
    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::   &
       uvel, vvel         ! u and v components of velocity

    integer, intent(in) ::    &
       matrix_order       ! order of matrix = number of rows

    type(sparse_matrix_type), intent(inout) ::  &
       matrix             ! sparse matrix, defined in glimmer_sparse_type
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

       do kA = -1, 1
       do jA = -1, 1
       do iA = -1, 1

          if ( (k+kA >= 1 .and. k+kA <= nz)         &
                          .and.                     &
               (i+iA >= 1 .and. i+iA <= nx-1)       &
                          .and.                     &
               (j+jA >= 1 .and. j+jA <= ny-1) ) then

             colA = NodeID(k+kA, i+iA, j+jA)   ! local ID for neighboring node
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
             m = indxA( iA, jA, kA)

             ! Avu 
             val = Avu(m,k,i,j)
             if (val /= 0.d0) then
                ct = ct + 1
                matrix%row(ct) = 2*rowA
                matrix%col(ct) = 2*colA - 1
                matrix%val(ct) = val
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

    ! Initialize the answer vector
    ! For efficiency, put the u and v terms for a given node adjacent in storage.

    do n = 1, nNodesSolve
       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)

       answer(2*n-1) = uvel(k,i,j)
       answer(2*n)   = vvel(k,i,j)

    enddo

    ! Set the rhs vector
    ! For efficiency, put the u and v terms for a given node adjacent in storage.

    do n = 1, nNodesSolve
       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)

       rhs(2*n-1) = bu(k,i,j)
       rhs(2*n)   = bv(k,i,j)

    enddo

  end subroutine slap_preprocess_3d

!****************************************************************************

  subroutine slap_preprocess_2d(nx,              ny,            &
                                nVerticesSolve,  vertexID,      &
                                iVertexIndex,    jVertexIndex,  &
                                indxA,                          &
                                Auu,             Auv,           &
                                Avu,             Avv,           &  
                                bu,              bv,            &
                                uvel,            vvel,          &
                                matrix_order,                   &
                                matrix,          rhs,           &
                                answer)

    !----------------------------------------------------------------
    ! This subroutine is analogous to slap_preprocess above, but for a 2D SSA solve 
    !
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
       nVerticesSolve           ! number of vertices where we solve for velocity

    integer, dimension(nx-1,ny-1), intent(in) ::  &
       vertexID                 ! local ID for each vertex

    integer, dimension(:), intent(in) ::   &
       iVertexIndex, jVertexIndex   ! i and j indices of active vertices

    integer, dimension(-1:1,-1:1), intent(in) :: &
       indxA                 ! maps relative (x,y,z) coordinates to an index between 1 and 9
                             ! index order is (i,j)

    real(dp), dimension(9,nx-1,ny-1), intent(in) ::  &
       Auu, Auv,       &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv              ! 1st dimension = vertex and its nearest neighbors in x and y direction 
                             ! other dimensions = (i,j) indices

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       bu, bv                ! assembled load (rhs) vector, divided into 2 parts
                          
    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       uvel, vvel            ! u and v components of velocity for 2D solve

    integer, intent(in) ::    &
       matrix_order       ! order of matrix = number of rows

    type(sparse_matrix_type), intent(inout) ::  &
       matrix             ! sparse matrix, defined in glimmer_sparse_types
                          ! includes nonzeros, order, col, row, val 

    real(dp), dimension(:), intent(out) ::   &
       rhs,             & ! right-hand-side (b) in Ax = b
       answer             ! answer (x) in Ax = b

    !---------------------------------------------------------
    ! Local variables
    !---------------------------------------------------------

    integer :: i, j, iA, jA, m, mm, n, ct

    integer :: rowA, colA   ! row and column of A submatrices (order = nVerticesSolve)

    real(dp) :: val         ! value of matrix coefficient
    
    ! Set the nonzero coefficients of the sparse matrix 

    ct = 0

    do rowA = 1, nVerticesSolve

       i = iVertexIndex(rowA)
       j = jVertexIndex(rowA)

       ! Load the nonzero values associated with Auu and Auv
       ! These are assigned a value of matrix%row = 2*rowA - 1

       do jA = -1, 1
       do iA = -1, 1

          if ( (i+iA >= 1 .and. i+iA <= nx-1)       &
                          .and.                     &
               (j+jA >= 1 .and. j+jA <= ny-1) ) then

             colA = vertexID(i+iA, j+jA)   ! local ID for neighboring vertex
             m = indxA(iA,jA)

             ! Auu
             val = Auu(m,i,j)
             if (val /= 0.d0) then
                ct = ct + 1
                matrix%row(ct) = 2*rowA - 1
                matrix%col(ct) = 2*colA - 1                    
                matrix%val(ct) = val
             endif

             ! Auv 
             val = Auv(m,i,j)
             if (val /= 0.d0) then
                ct = ct + 1
                matrix%row(ct) = 2*rowA - 1
                matrix%col(ct) = 2*colA
                matrix%val(ct) = val
             endif

          endif     ! i+iA and j+jA in bounds

       enddo        ! iA
       enddo        ! jA

       ! Load the nonzero values associated with Avu and Avv
       ! These are assigned a value of matrix%row = 2*rowA

       do jA = -1, 1
       do iA = -1, 1

          if ( (i+iA >= 1 .and. i+iA <= nx-1)       &
                          .and.                     &
               (j+jA >= 1 .and. j+jA <= ny-1) ) then

             colA = vertexID(i+iA, j+jA)   ! ID for neighboring vertex
             m = indxA(iA, jA)

             ! Avu 
             val = Avu(m,i,j)
             if (val /= 0.d0) then
                ct = ct + 1
                matrix%row(ct) = 2*rowA
                matrix%col(ct) = 2*colA - 1
                matrix%val(ct) = val
             endif

             ! Avv 
             val = Avv(m,i,j)
             if (val /= 0.d0) then
                ct = ct + 1
                matrix%row(ct) = 2*rowA
                matrix%col(ct) = 2*colA
                matrix%val(ct) = val
             endif

          endif     ! i+iA and j+jA in bounds

       enddo        ! iA
       enddo        ! jA

    enddo           ! rowA 

    ! Set basic matrix parameters.

    matrix%order = matrix_order
    matrix%nonzeros = ct
    matrix%symmetric = .false.

    ! Set the answer vector
    ! For efficiency, put the u and v terms for a given node adjacent in storage.

    do n = 1, nVerticesSolve
       i = iVertexIndex(n)
       j = jVertexIndex(n)

       answer(2*n-1) = uvel(i,j)
       answer(2*n)   = vvel(i,j)

    enddo

    ! Set the rhs vector
    ! For efficiency, put the u and v terms for a given node adjacent in storage.

    do n = 1, nVerticesSolve
       i = iVertexIndex(n)
       j = jVertexIndex(n)

       rhs(2*n-1) = bu(i,j)
       rhs(2*n)   = bv(i,j)

    enddo

  end subroutine slap_preprocess_2d

!****************************************************************************

  subroutine slap_postprocess_3d(nNodesSolve,                            &
                                 iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                                 answer,       resid_vec,                &
                                 uvel,         vvel,                     &
                                 resid_u,      resid_v)

  ! Extract the velocities from the SLAP solution vector.
                                            
    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) :: nNodesSolve     ! number of nodes where we solve for velocity

    real(dp), dimension(:), intent(in) ::  &
       answer,           &! velocity solution vector
       resid_vec          ! residual vector

    integer, dimension(:), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex  ! i, j and k indices of active nodes

    real(dp), dimension(:,:,:), intent(inout) ::   &
       uvel, vvel,       &! u and v components of velocity
       resid_u, resid_v   ! u and v components of residual

    integer :: i, j, k, n

    do n = 1, nNodesSolve

       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)

       uvel(k,i,j) = answer(2*n-1)
       vvel(k,i,j) = answer(2*n)

       resid_u(k,i,j) = resid_vec(2*n-1)
       resid_v(k,i,j) = resid_vec(2*n)

    enddo

  end subroutine slap_postprocess_3d

!****************************************************************************

  subroutine slap_postprocess_2d(nVerticesSolve,               &
                                 iVertexIndex, jVertexIndex,   &
                                 answer,       resid_vec,      &
                                 uvel,         vvel,           &
                                 resid_u,      resid_v)

    ! Extract the velocities from the SLAP solution vector.
                                            
    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) :: nVerticesSolve   ! number of vertices where we solve for velocity

    real(dp), dimension(:), intent(in) ::  &
       answer,           &! velocity solution vector
       resid_vec          ! residual vector

    integer, dimension(:), intent(in) ::   &
       iVertexIndex, jVertexIndex  ! i and j indices of active vertices

    real(dp), dimension(:,:), intent(inout) ::   &
       uvel, vvel,       &! u and v components of velocity
       resid_u, resid_v   ! u and v components of residual

    integer :: i, j, n

    do n = 1, nVerticesSolve

       i = iVertexIndex(n)
       j = jVertexIndex(n)

       uvel(i,j) = answer(2*n-1)
       vvel(i,j) = answer(2*n)

       resid_u(i,j) = resid_vec(2*n-1)
       resid_v(i,j) = resid_vec(2*n)

    enddo

  end subroutine slap_postprocess_2d

!****************************************************************************

  subroutine slap_compute_residual_vector(matrix,  answer,    &
                                          rhs,     resid_vec, &
                                          L2_norm, L2_norm_relative)

    ! Compute the residual vector Ax - b and its L2 norm.
    ! This subroutine assumes that the matrix is stored in triad (row/col/val) format.

    type(sparse_matrix_type), intent(in) ::  &
       matrix             ! sparse matrix, defined in glimmer_sparse_types
                          ! includes nonzeros, order, col, row, val 

    real(dp), dimension(:), intent(in) ::   &
       rhs,             & ! right-hand-side (b) in Ax = b
       answer             ! answer (x) in Ax = b

    real(dp), dimension(:), intent(out) ::   &
       resid_vec          ! residual vector

    real(dp), intent(out) ::    &
       L2_norm            ! L2 norm of residual vector, |Ax - b|
  
    real(dp), intent(out), optional ::    &
       L2_norm_relative   ! L2 norm of residual vector relative to rhs, |Ax - b| / |b|

    integer :: n, r, c
    
    real(dp) :: L2_norm_rhs  ! L2 norm of rhs vector, |b|

    resid_vec(:) = 0.d0

    do n = 1, matrix%nonzeros
       r = matrix%row(n)
       c = matrix%col(n)
       resid_vec(r) = resid_vec(r) + matrix%val(n)*answer(c)
    enddo

    L2_norm = 0.d0
    do r = 1, matrix%order
       resid_vec(r) = resid_vec(r) - rhs(r)
       L2_norm = L2_norm + resid_vec(r)*resid_vec(r)
    enddo

    L2_norm = sqrt(L2_norm)

    if (present(L2_norm_relative)) then

       L2_norm_rhs = 0.d0

       do r = 1, matrix%order
          L2_norm_rhs = L2_norm_rhs + rhs(r)*rhs(r)
       enddo

       L2_norm_rhs = sqrt(L2_norm_rhs)

       if (L2_norm_rhs > 0.d0) then
          L2_norm_relative = L2_norm / L2_norm_rhs
       else
          L2_norm_relative = 0.d0
       endif

    endif  ! present(L2_norm_relative)

  end subroutine slap_compute_residual_vector

!****************************************************************************

  subroutine slap_solve_test_matrix (matrix_order, whichsparse)

    ! solve a small test matrix

    integer, intent(in) :: &
       matrix_order,       & ! matrix order
       whichsparse           ! solution method (0=BiCG, 1=GMRES, 2==PCG_INCH)

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

  end subroutine slap_solve_test_matrix

!****************************************************************************

end module glissade_velo_higher_slap

!****************************************************************************
