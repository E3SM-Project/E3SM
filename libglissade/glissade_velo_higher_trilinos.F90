!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_velo_higher_trilinos.F90 - part of the Community Ice Sheet Model (CISM)  
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
! and used to process data before and after linking to Trilinos solver routines.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  module glissade_velo_higher_trilinos

    use glimmer_global, only: dp
!    use glimmer_log, only: write_log
    use parallel

    implicit none
    private

#ifdef TRILINOS

    public :: trilinos_global_id_3d,        trilinos_global_id_2d,         &
              trilinos_fill_pattern_3d,     trilinos_fill_pattern_2d,      &
              trilinos_assemble_3d,         trilinos_assemble_2d,          &
              trilinos_init_velocity_3d,    trilinos_init_velocity_2d,     &
              trilinos_extract_velocity_3d, trilinos_extract_velocity_2d,  &
              trilinos_test

  contains

!****************************************************************************

  subroutine trilinos_global_id_3d(nx,         ny,         nz,   &
                                   nNodesSolve,                  &
                                   iNodeIndex, jNodeIndex, kNodeIndex,  &
                                   global_node_id,               &
                                   active_owned_unknown_map)

    !----------------------------------------------------------------
    ! Compute global IDs needed to initialize the Trilinos solver
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,             &  ! number of grid cells in each direction
       nz                     ! number of vertical levels where velocity is computed

    integer, intent(in) ::             &
       nNodesSolve            ! number of nodes where we solve for velocity

    integer, dimension((nx-1)*(ny-1)*nz), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of nodes

    integer, dimension(nz,nx-1,ny-1), intent(out) ::  &
       global_node_id      ! unique global ID for nodes on this processor

    integer, dimension(2*nNodesSolve), intent(out) ::   &
       active_owned_unknown_map

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer, dimension(nx-1,ny-1) ::   &
       global_vertex_id    ! unique global ID for vertices on this processor

    integer :: gnx, gny

    integer :: i, j, k, n

    !----------------------------------------------------------------
    ! Compute unique global IDs for nodes.
    !----------------------------------------------------------------

    global_vertex_id(:,:) = 0

    do j = nhalo+1, ny-nhalo         ! locally owned vertices only
       do i = nhalo+1, nx-nhalo
          gnx = ewlb + i - 1   ! global x index
          gny = nslb + j - 1   ! global y index
          global_vertex_id(i,j) = (gny-1)*global_ewn + gnx
       enddo
    enddo

    call staggered_parallel_halo(global_vertex_id)

    do j = 1, ny-1         ! loop over all vertices, including halo
       do i = 1, nx-1
          do k = 1, nz
             global_node_id(k,i,j) = (global_vertex_id(i,j)-1)*nz + k
          enddo
       enddo
    enddo

    !----------------------------------------------------------------
    ! Associate a unique global index with each unknown on the active nodes
    ! owned by this processor.
    !----------------------------------------------------------------

    do n = 1, nNodesSolve
       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)
       active_owned_unknown_map(2*n-1) = 2*global_node_id(k,i,j) - 1  ! u unknowns
       active_owned_unknown_map(2*n)   = 2*global_node_id(k,i,j)      ! v unknowns
    enddo

  end subroutine trilinos_global_id_3d

!****************************************************************************

  subroutine trilinos_global_id_2d(nx,             ny,           &
                                   nVerticesSolve,               &
                                   iVertexIndex,   jVertexIndex, &
                                   global_vertex_id,             &
                                   active_owned_unknown_map)

    !----------------------------------------------------------------
    ! Compute global IDs needed to initialize the Trilinos solver
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                  ! number of grid cells in each direction

    integer, intent(in) ::             &
       nVerticesSolve          ! number of nodes where we solve for velocity

    integer, dimension((nx-1)*(ny-1)), intent(in) ::   &
       iVertexIndex, jVertexIndex   ! i and j indices of vertices

    integer, dimension(nx-1,ny-1), intent(out) ::  &
       global_vertex_id        ! unique global ID for nodes on this processor

    integer, dimension(2*nVerticesSolve), intent(out) ::   &
       active_owned_unknown_map

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: gnx, gny

    integer :: i, j, n

    !----------------------------------------------------------------
    ! Compute unique global IDs for vertices.
    !----------------------------------------------------------------

    global_vertex_id(:,:) = 0

    do j = nhalo+1, ny-nhalo         ! locally owned vertices only
       do i = nhalo+1, nx-nhalo
          gnx = ewlb + i - 1   ! global x index
          gny = nslb + j - 1   ! global y index
          global_vertex_id(i,j) = (gny-1)*global_ewn + gnx
       enddo
    enddo

    call staggered_parallel_halo(global_vertex_id)

    !----------------------------------------------------------------
    ! Associate a unique global index with each unknown on the active vertices
    ! owned by this processor.
    !----------------------------------------------------------------

    do n = 1, nVerticesSolve
       i = iVertexIndex(n)
       j = jVertexIndex(n)
       active_owned_unknown_map(2*n-1) = 2*global_vertex_id(i,j) - 1  ! u unknowns
       active_owned_unknown_map(2*n)   = 2*global_vertex_id(i,j)      ! v unknowns
    enddo

  end subroutine trilinos_global_id_2d

!****************************************************************************

  subroutine trilinos_fill_pattern_3d(nx,            ny,          nz,         &
                                      active_vertex, nNodesSolve,             &
                                      iNodeIndex,    jNodeIndex,  kNodeIndex, &
                                      indxA,         Afill)

    !------------------------------------------------------------------------
    ! Construct logical arrays identifying which matrix elements are nonzero.
    ! For the Trilinos solver, the number of matrix entries must be held fixed
    !  from one iteration to the next.  The logical arrays are used to
    !  satisfy this requirement.
    ! For now we simply set A**_fill = .true. everywhere in the row corresponding
    !  to each active node.
    ! Later, we could use boundary logic to set A**_fill = .false for some
    !  columns, to avoid including matrix values that are always zero.
    !------------------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,             &  ! number of grid cells in each direction
       nz                     ! number of vertical levels where velocity is computed

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F

    integer, intent(in) :: nNodesSolve     ! number of nodes where we solve for velocity

    integer, dimension(:), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex  ! i, j and k indices of active nodes

    integer, dimension(-1:1,-1:1,-1:1), intent(in) :: &
       indxA                 ! maps relative (x,y,z) coordinates to an index between 1 and 27                                                 
                             ! index order is (i,j,k)
   
    logical, dimension(27,nz,nx-1,ny-1), intent(out) ::  &
       Afill        ! true wherever the matrix value is potentially nonzero
                    ! and should be sent to Trilinos

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j, k, m, n, iA, jA, kA

    Afill(:,:,:,:) = .false.

    ! Loop over active nodes

    do n = 1, nNodesSolve

       i = iNodeIndex(n)
       j = jNodeIndex(n)

       if (active_vertex(i,j)) then

          k = kNodeIndex(n)

          do kA = -1,1
          do jA = -1,1
          do iA = -1,1

             if ( (k+kA >= 1 .and. k+kA <= nz)         &
                             .and.                     &
                  (i+iA >= 1 .and. i+iA <= nx-1)       &
                             .and.                     &
                  (j+jA >= 1 .and. j+jA <= ny-1) ) then

                if (active_vertex(i+iA,j+jA)) then

                   m = indxA(iA,jA,kA)
                   Afill(m,k,i,j) = .true.
                   
                endif  ! active_vertex(i+iA,j+jA)

             endif     ! neighbor node is in bounds
               
          enddo        ! iA
          enddo        ! jA
          enddo        ! kA

       endif           ! active_vertex(i,j)

    enddo              ! n

  end subroutine trilinos_fill_pattern_3d

!****************************************************************************

  subroutine trilinos_fill_pattern_2d(nx,            ny,             &
                                      active_vertex, nVerticesSolve, &
                                      iVertexIndex,  jVertexIndex,   &
                                      indxA,         Afill)

    !------------------------------------------------------------------------
    ! Construct logical arrays identifying which matrix elements are nonzero.
    !
    ! This subroutine is similar to trilinos_fill_pattern_3d, but modified                                                                       
    !  to solve for x and y at a single horizontal level, as in the                                                                           
    !  shallow-shelf approximation.
    !------------------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                 ! number of grid cells in each direction

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F

    integer, intent(in) :: nVerticesSolve     ! number of vertices where we solve for velocity

    integer, dimension(:), intent(in) ::   &
       iVertexIndex, jVertexIndex    ! i and j indices of active vertices

    integer, dimension(-1:1,-1:1), intent(in) :: &
       indxA             ! maps relative (x,y,z) coordinates to an index between 1 and 9  
                         ! index order is (i,j)
   
    logical, dimension(9,nx-1,ny-1), intent(out) ::  &
       Afill        ! true wherever the matrix value is potentially nonzero
                    ! and should be sent to Trilinos

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j, m, n, iA, jA

    Afill(:,:,:) = .false.

    ! Loop over active vertices

    do n = 1, nVerticesSolve

       i = iVertexIndex(n)
       j = jVertexIndex(n)

       if (active_vertex(i,j)) then

          do jA = -1,1
          do iA = -1,1

             if ( (i+iA >= 1 .and. i+iA <= nx-1)       &
                             .and.                     &
                  (j+jA >= 1 .and. j+jA <= ny-1) ) then

                if (active_vertex(i+iA,j+jA)) then

                   m = indxA(iA,jA)
                   Afill(m,i,j) = .true.
                   
                endif  ! active_vertex(i+iA,j+jA)

             endif     ! neighbor node is in bounds
               
          enddo        ! iA
          enddo        ! jA

       endif           ! active_vertex(i,j)

    enddo              ! n

  end subroutine trilinos_fill_pattern_2d
               
!****************************************************************************

  subroutine trilinos_assemble_3d(nx,           ny,            nz,          &   
                                  nNodesSolve,  global_node_id,             &
                                  iNodeIndex,   jNodeIndex,    kNodeIndex,  &
                                  indxA,        Afill,          &
                                  Auu,          Auv,            &
                                  Avu,          Avv,            &
                                  bu,           bv)

    !------------------------------------------------------------------------
    ! Given Auu, bu, etc., assemble the matrix and RHS in a form
    ! suitable for Trilinos.
    !
    ! Note: Trilinos requires that the matrix fill pattern is unchanged from
    !       one outer iteration to the next. This requirement is currently enforced
    !       by sending all 54 columns to Trilinos for each row (since A**_fill
    !       is true everywhere), even though some columns may always equal zero. 
    !       With some more work, we should be able to remove some of these columns 
    !       for nodes at the boundary.
    !------------------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,             &  ! number of grid cells in each direction
       nz                     ! number of vertical levels where velocity is computed

    integer, intent(in) ::             &
       nNodesSolve            ! number of nodes where we solve for velocity

    integer, dimension(nz,nx-1,ny-1), intent(in) ::  &
       global_node_id      ! unique global ID for nodes on this processor

    integer, dimension((nx-1)*(ny-1)*nz), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of nodes

    integer, dimension(-1:1,-1:1,-1:1), intent(in) :: &
       indxA                 ! maps relative (x,y,z) coordinates to an index between 1 and 27                                                 
                             ! index order is (i,j,k)

    logical, dimension(27,nz,nx-1,ny-1), intent(in) ::  &
       Afill              ! true for matrix values to be sent to Trilinos

    real(dp), dimension(27,nz,nx-1,ny-1), intent(in) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv           ! 1st dimension = node and its nearest neighbors in x, y and z direction 
                          ! other dimensions = (k,i,j) indices

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       bu, bv             ! assembled load (rhs) vector, divided into 2 parts
                          
    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: global_row   ! global ID for this matrix row

    integer :: ncol         ! number of columns with nonzero entries in this row

    integer, dimension(54) ::   &
       global_column        ! global ID for this column
                            ! 54 is max number of columns with nonzero entries

    real(dp), dimension(54) ::  &
       matrix_value         ! matrix value for this column

    real(dp) :: rhs_value   ! right-hand side value (bu or bv)

    integer :: i, j, k, m, n, iA, jA, kA

!WHL - debug
    integer :: nc

    do n = 1, nNodesSolve

       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)

       ! uvel equation for this node

       global_row = 2*global_node_id(k,i,j) - 1

!WHL - debug
!       print*, ' '
!       print*, 'n, i, j, k', n, i, j, k
!       print*, 'global_node_id, global_row:', global_node_id(k,i,j), global_row
      
       ncol = 0
       global_column(:) = 0
       matrix_value(:) = 0.d0

       do kA = -1, 1
       do jA = -1, 1
       do iA = -1, 1

          if ( (k+kA >= 1 .and. k+kA <= nz)         &
                          .and.                     &
               (i+iA >= 1 .and. i+iA <= nx-1)       &
                          .and.                     &
               (j+jA >= 1 .and. j+jA <= ny-1) ) then

             m = indxA(iA,jA,kA)

             if (Afill(m,k,i,j)) then

                ncol = ncol + 1
                global_column(ncol) = 2*global_node_id(k+kA,i+iA,j+jA) - 1
                matrix_value(ncol) = Auu(m,k,i,j)

                ncol = ncol + 1
                global_column(ncol) = 2*global_node_id(k+kA,i+iA,j+jA)
                matrix_value(ncol) = Auv(m,k,i,j)

             endif

          endif   ! i+iA, j+jA, k+kA in bounds
       enddo      ! iA
       enddo      ! jA
       enddo      ! kA

       rhs_value = bu(k,i,j)

       call insertrowtgs(global_row, ncol, global_column, matrix_value, rhs_value)

       ! vvel equation for this node

       global_row = 2*global_node_id(k,i,j)
       
       ncol = 0
       global_column(:) = 0
       matrix_value(:) = 0.d0

       do kA = -1, 1
       do jA = -1, 1
       do iA = -1, 1

          if ( (k+kA >= 1 .and. k+kA <= nz)         &
                          .and.                     &
               (i+iA >= 1 .and. i+iA <= nx-1)       &
                          .and.                     &
               (j+jA >= 1 .and. j+jA <= ny-1) ) then

             m = indxA(iA,jA,kA)

             if (Afill(m,k,i,j)) then

                ncol = ncol + 1
                global_column(ncol) = 2*global_node_id(k+kA,i+iA,j+jA) - 1
                matrix_value(ncol) = Avu(m,k,i,j)

                ncol = ncol + 1
                global_column(ncol) = 2*global_node_id(k+kA,i+iA,j+jA)
                matrix_value(ncol) = Avv(m,k,i,j)

             endif

          endif   ! i+iA, j+jA, k+kA in bounds
       enddo    ! iA
       enddo    ! jA
       enddo    ! kA

       rhs_value = bv(k,i,j)

       call insertrowtgs(global_row, ncol, global_column, matrix_value, rhs_value)

    enddo   ! nNodesSolve

  end subroutine trilinos_assemble_3d

!****************************************************************************

  subroutine trilinos_assemble_2d(nx,             ny,                &   
                                  nVerticesSolve, global_vertex_id,  &
                                  iVertexIndex,   jVertexIndex,      &
                                  indxA,          Afill,             &
                                  Auu,            Auv,               &
                                  Avu,            Avv,               &
                                  bu,             bv)

    !------------------------------------------------------------------------
    ! Given Auu, bu, etc., assemble the matrix and RHS in a form
    ! suitable for Trilinos.
    !
    ! Note: Trilinos requires that the matrix fill pattern is unchanged from
    !       one outer iteration to the next. This requirement is currently enforced
    !       by sending all 18 columns to Trilinos for each row, even though some 
    !       column entries may always equal zero. 
    !       With some more work, we could remove some of these columns 
    !       for nodes at the boundary.
    !------------------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                 ! number of grid cells in each direction

    integer, intent(in) ::             &
       nVerticesSolve         ! number of vertices where we solve for velocity

    integer, dimension(nx-1,ny-1), intent(in) ::  &
       global_vertex_id      ! unique global ID for vertices on this processor

    integer, dimension((nx-1)*(ny-1)), intent(in) ::   &
       iVertexIndex, jVertexIndex   ! i and j indices of vertices

    integer, dimension(-1:1,-1:1), intent(in) :: &
       indxA                 ! maps relative (x,y) coordinates to an index between 1 and 9
                             ! index order is (i,j)

    logical, dimension(9,nx-1,ny-1), intent(in) ::  &
       Afill              ! true for matrix values to be sent to Trilinos

    real(dp), dimension(9,nx-1,ny-1), intent(in) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv           ! 1st dimension = node and its nearest neighbors in x, y and z direction 
                          ! other dimensions = (i,j) indices

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       bu, bv             ! assembled load (rhs) vector, divided into 2 parts
                          
    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: global_row   ! global ID for this matrix row

    integer :: ncol         ! number of columns with nonzero entries in this row

    integer, dimension(18) ::   &
       global_column        ! global ID for this column
                            ! 18 is max number of columns with nonzero entries

    real(dp), dimension(18) ::  &
       matrix_value         ! matrix value for this column

    real(dp) :: rhs_value   ! right-hand side value (bu or bv)

    integer :: i, j, m, n, iA, jA

    do n = 1, nVerticesSolve

       i = iVertexIndex(n)
       j = jVertexIndex(n)

       ! uvel equation for this node

       global_row = 2*global_vertex_id(i,j) - 1

       ncol = 0
       global_column(:) = 0
       matrix_value(:) = 0.d0

       do jA = -1, 1
       do iA = -1, 1

          if ( (i+iA >= 1 .and. i+iA <= nx-1)       &
                          .and.                     &
               (j+jA >= 1 .and. j+jA <= ny-1) ) then

             m = indxA(iA,jA)

             if (Afill(m,i,j)) then

                ncol = ncol + 1
                global_column(ncol) = 2*global_vertex_id(i+iA,j+jA) - 1
                matrix_value(ncol) = Auu(m,i,j)

                ncol = ncol + 1
                global_column(ncol) = 2*global_vertex_id(i+iA,j+jA)
                matrix_value(ncol) = Auv(m,i,j)

             endif

          endif   ! i+iA, j+jA in bounds
       enddo    ! iA
       enddo    ! jA

       rhs_value = bu(i,j)

       call insertrowtgs(global_row, ncol, global_column, matrix_value, rhs_value)

       ! vvel equation for this node

       global_row = 2*global_vertex_id(i,j)
       
       ncol = 0
       global_column(:) = 0
       matrix_value(:) = 0.d0

       do jA = -1, 1
       do iA = -1, 1

          if ( (i+iA >= 1 .and. i+iA <= nx-1)       &
                          .and.                     &
               (j+jA >= 1 .and. j+jA <= ny-1) ) then

             m = indxA(iA,jA)

             if (Afill(m,i,j)) then

                ncol = ncol + 1
                global_column(ncol) = 2*global_vertex_id(i+iA,j+jA) - 1
                matrix_value(ncol) = Avu(m,i,j)

                ncol = ncol + 1
                global_column(ncol) = 2*global_vertex_id(i+iA,j+jA)
                matrix_value(ncol) = Avv(m,i,j)

             endif

          endif   ! i+iA, j+jA in bounds
       enddo      ! iA
       enddo      ! jA

       rhs_value = bv(i,j)

       call insertrowtgs(global_row, ncol, global_column, matrix_value, rhs_value)

    enddo   ! nVerticesSolve

  end subroutine trilinos_assemble_2d

!****************************************************************************

  subroutine trilinos_init_velocity_3d(nx,           ny,                       &
                                       nz,           nNodesSolve,              &
                                       iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                                       uvel,         vvel,                     &
                                       velocityResult)

    ! Copy the initial velocities into the Trilinos solution vector.

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,             &  ! number of grid cells in each direction
       nz                     ! number of vertical levels where velocity is computed

    integer, intent(in) ::             &
       nNodesSolve            ! number of nodes where we solve for velocity

    integer, dimension((nx-1)*(ny-1)*nz), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of nodes

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::   &
       uvel, vvel             ! u and v components of velocity

    real(dp), dimension(2*nNodesSolve), intent(out) :: &
       velocityResult         ! initial velocity solution vector for Trilinos

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j, k, n

    velocityResult(:) = 0.d0

    do n = 1, nNodesSolve
       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)
       velocityResult(2*n-1) = uvel(k,i,j) 
       velocityResult(2*n)   = vvel(k,i,j)
    enddo

  end subroutine trilinos_init_velocity_3d

!****************************************************************************

  subroutine trilinos_init_velocity_2d(nx,             ny,           &
                                       nVerticesSolve,               &
                                       iVertexIndex,   jVertexIndex, &
                                       uvel,           vvel,         &
                                       velocityResult)

    ! Copy the initial velocities into the Trilinos solution vector.

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                 ! number of grid cells in each direction

    integer, intent(in) ::             &
       nVerticesSolve         ! number of vertices where we solve for velocity

    integer, dimension((nx-1)*(ny-1)), intent(in) ::   &
       iVertexIndex, jVertexIndex   ! i and j indices of vertices

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       uvel, vvel             ! u and v components of velocity

    real(dp), dimension(2*nVerticesSolve), intent(out) :: &
       velocityResult         ! initial velocity solution vector for Trilinos

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j, n

    velocityResult(:) = 0.d0

    do n = 1, nVerticesSolve
       i = iVertexIndex(n)
       j = jVertexIndex(n)
       velocityResult(2*n-1) = uvel(i,j) 
       velocityResult(2*n)   = vvel(i,j)
    enddo

  end subroutine trilinos_init_velocity_2d

!****************************************************************************

  subroutine trilinos_extract_velocity_3d(nx,           ny,                       &
                                          nz,           nNodesSolve,              &
                                          iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                                          velocityResult,                         &
                                          uvel,         vvel)

    ! Extract the velocities from the Trilinos solution vector.

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,             &  ! number of grid cells in each direction
       nz                     ! number of vertical levels where velocity is computed

    integer, intent(in) ::             &
       nNodesSolve            ! number of nodes where we solve for velocity

    integer, dimension((nx-1)*(ny-1)*nz), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of nodes

    real(dp), dimension(2*nNodesSolve), intent(in) :: &
       velocityResult         ! velocity solution vector from Trilinos

    real(dp), dimension(nz,nx-1,ny-1), intent(out) ::   &
       uvel, vvel             ! u and v components of velocity

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j, k, n

    uvel(:,:,:) = 0.d0
    vvel(:,:,:) = 0.d0

    do n = 1, nNodesSolve
       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)
       uvel(k,i,j) = velocityResult(2*n-1)
       vvel(k,i,j) = velocityResult(2*n)
    enddo

  end subroutine trilinos_extract_velocity_3d

!****************************************************************************

  subroutine trilinos_extract_velocity_2d(nx,             ny,            &
                                          nVerticesSolve,                &
                                          iVertexIndex,   jVertexIndex,  &
                                          velocityResult,                &
                                          uvel,           vvel)

    ! Extract the velocities from the Trilinos solution vector.

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                 ! number of grid cells in each direction

    integer, intent(in) ::   &
       nVerticesSolve         ! number of vertices where we solve for velocity

    integer, dimension((nx-1)*(ny-1)), intent(in) ::   &
       iVertexIndex, jVertexIndex  ! i and j indices of vertices

    real(dp), dimension(2*nVerticesSolve), intent(in) :: &
       velocityResult         ! velocity solution vector from Trilinos

    real(dp), dimension(nx-1,ny-1), intent(out) ::   &
       uvel, vvel             ! u and v components of velocity

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j, n

    uvel(:,:) = 0.d0
    vvel(:,:) = 0.d0

    do n = 1, nVerticesSolve
       i = iVertexIndex(n)
       j = jVertexIndex(n)
       uvel(i,j) = velocityResult(2*n-1)
       vvel(i,j) = velocityResult(2*n)
    enddo

  end subroutine trilinos_extract_velocity_2d

!****************************************************************************

  subroutine trilinos_test

    !--------------------------------------------------------
    ! Small test matrices for Trilinos solver
    !--------------------------------------------------------
    
    use parallel
      
    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------
 
    integer :: nNodesSolve
    integer, dimension(:), allocatable :: &
       active_owned_unknown_map  ! map of global IDs
    integer :: global_row   ! global ID for this matrix row
    integer :: ncol         ! number of columns with nonzero entries in this row
    integer, dimension(:), allocatable ::   &
       global_column        ! global ID for this column
    real(dp), dimension(:), allocatable ::  &
       matrix_value         ! matrix value for this column
    real(dp) :: rhs_value   ! right-hand side value (bu or bv)
    real(dp), dimension(:), allocatable ::   &
       velocityResult     ! velocity solution vector from Trilinos

    if (main_task) then
       print*, ' '
       print*, 'Solve trilinos test matrix, tasks =', tasks
    endif

    if (tasks == 1) then

       ! Set up 2x2 matrix problem on 1 processor

       ! Here is the problem:
       ! |  1   2  | | 1 |    | 3 |
       ! |  3   4  | | 1 | =  | 7 |

       nNodesSolve = 1
       allocate(active_owned_unknown_map(2*nNodesSolve))
       active_owned_unknown_map(:) = (/ 1,2 /)
       print*, 'initializetgs, rank =', this_rank
       call initializetgs(2*nNodesSolve, active_owned_unknown_map, comm)

       ! insert rows

       allocate(global_column(2))
       allocate(matrix_value(2))

       ! row 1 (global ID = 1)
       global_row = 1   
       ncol = 2
       global_column(:) = (/ 1,2 /)
       matrix_value(:)  = (/ 1,2 /)
       rhs_value = 3
       print*, 'insertrowtgs, rank, row =', this_rank, global_row
       call insertrowtgs(global_row, ncol, global_column, matrix_value, rhs_value)

       ! row 2 (global ID = 2)
       global_row = 2
       ncol = 2
       global_column(:) = (/ 1,2 /)
       matrix_value(:)  = (/ 3,4 /)
       rhs_value = 7
       print*, 'insertrowtgs, rank, row =', this_rank, global_row
       call insertrowtgs(global_row, ncol, global_column, matrix_value, rhs_value)

       ! solve
       allocate(velocityResult(2*nNodesSolve))
       print*, 'solvevelocitytgs, rank =', this_rank
       call solvevelocitytgs(velocityResult)

       ! print solution
       print*, 'rank, solution:', this_rank, velocityResult(:)

    elseif (tasks == 2) then

       ! Set up a 4x4 matrix problem on 2 processors
       ! This one has 16 active unknowns

       ! Here is the problem:
       ! |  1   2   0   3  | | 1 |    | 7  |
       ! |  4   5   6   0  | | 0 |    | 10 |
       ! |  7   8   9  10  | | 1 | =  | 36 |
       ! |  0  11  12  13  | | 2 |    | 38 |

       ! initialize

       allocate(global_column(4))
       allocate(matrix_value(4))

       nNodesSolve = 1
       allocate(active_owned_unknown_map(2*nNodesSolve))
       if (this_rank==0) then
          active_owned_unknown_map(:) = (/ 1,3 /)
       elseif (this_rank==1) then
          active_owned_unknown_map(:) = (/ 4,6 /)
       endif
       print*, 'initializetgs, rank =', this_rank
       call initializetgs(2*nNodesSolve, active_owned_unknown_map, comm)

       ! insert rows

       if (this_rank==0) then

          ! row 1 (global ID = 1)
          global_row = 1   
          ncol = 3
          global_column(:) = (/ 1,3,6,0 /)
          matrix_value(:)  = (/ 1,2,3,0 /)
          rhs_value = 7
          print*, 'insertrowtgs, rank, row =', this_rank, global_row
          call insertrowtgs(global_row, ncol, global_column, matrix_value, rhs_value)
       
          ! row 2 (global ID = 3)
          global_row = 3   
          ncol = 3
          global_column(:) = (/ 1,3,4,0 /)
          matrix_value(:)  = (/ 4,5,6,0 /)
          rhs_value = 10
          print*, 'insertrowtgs, rank, row =', this_rank, global_row
          call insertrowtgs(global_row, ncol, global_column, matrix_value, rhs_value)

       elseif (this_rank==1) then

          ! row 1 (global ID = 4)
          global_row = 4   
          ncol = 4
          global_column(:) = (/ 1,3,4,6 /)
          matrix_value(:)  = (/ 7,8,9,10 /)
          rhs_value = 36
          print*, 'insertrowtgs, rank, row =', this_rank, global_row
          call insertrowtgs(global_row, ncol, global_column, matrix_value, rhs_value)

          ! row 2 (global ID = 6)
          global_row = 6   
          ncol = 3
          global_column(:) = (/ 3,4,6,0 /)
          matrix_value(:)  = (/ 11,12,13,0 /)
          rhs_value = 38
          print*, 'insertrowtgs, rank, row =', this_rank, global_row
          call insertrowtgs(global_row, ncol, global_column, matrix_value, rhs_value)

       endif

       ! solve
       allocate(velocityResult(2*nNodesSolve))
       print*, 'solvevelocitytgs, rank =', this_rank
       call solvevelocitytgs(velocityResult)

       ! print solution
       print*, 'rank, solution:', this_rank, velocityResult(:)

       deallocate(active_owned_unknown_map)
       deallocate(velocityResult)

       ! Set up 4x4 matrix problem on 2 processors
       ! This one has 12 active unknowns

       ! Here is the problem:
       ! |  1   2   3   0  | | 1 |    | 4  |
       ! |  0   4   5   6  | | 0 |    | 17 |
       ! |  7   8   9   0  | | 1 | =  | 16 |
       ! |  0   10  11  12 | | 2 |    | 35 |

       ! initialize

       nNodesSolve = 1
       allocate(active_owned_unknown_map(2*nNodesSolve))
       if (this_rank==0) then
          active_owned_unknown_map(:) = (/ 1,3 /)
       elseif (this_rank==1) then
          active_owned_unknown_map(:) = (/ 4,6 /)
       endif
       print*, 'initializetgs, rank =', this_rank
       call initializetgs(2*nNodesSolve, active_owned_unknown_map, comm)

       ! insert rows

       if (this_rank==0) then

          ! row 1 (global ID = 1)
          global_row = 1   
          ncol = 3
          global_column(:) = (/ 1,3,4,0 /)
          matrix_value(:)  = (/ 1,2,3,0 /)
          rhs_value = 4
          print*, 'insertrowtgs, rank, row =', this_rank, global_row
          call insertrowtgs(global_row, ncol, global_column, matrix_value, rhs_value)
       
          ! row 2 (global ID = 3)
          global_row = 3   
          ncol = 3
          global_column(:) = (/ 3,4,6,0 /)
          matrix_value(:)  = (/ 4,5,6,0 /)
          rhs_value = 17
          print*, 'insertrowtgs, rank, row =', this_rank, global_row
          call insertrowtgs(global_row, ncol, global_column, matrix_value, rhs_value)

       elseif (this_rank==1) then

          ! row 1 (global ID = 4)
          global_row = 4   
          ncol = 3
          global_column(:) = (/ 1,3,4,0 /)
          matrix_value(:)  = (/ 7,8,9,0 /)
          rhs_value = 16
          print*, 'insertrowtgs, rank, row =', this_rank, global_row
          call insertrowtgs(global_row, ncol, global_column, matrix_value, rhs_value)

          ! row 2 (global ID = 6)
          global_row = 6   
          ncol = 3
          global_column(:) = (/ 3,4,6,0 /)
          matrix_value(:)  = (/ 10,11,12,0 /)
          rhs_value = 35
          print*, 'insertrowtgs, rank, row =', this_rank, global_row
          call insertrowtgs(global_row, ncol, global_column, matrix_value, rhs_value)

       endif

       ! solve
       allocate(velocityResult(2*nNodesSolve))
       print*, 'solvevelocitytgs, rank =', this_rank
       call solvevelocitytgs(velocityResult)

       ! print solution
       print*, 'rank, solution:', this_rank, velocityResult(:)

       deallocate(active_owned_unknown_map)
       deallocate(velocityResult)

    else
       print*, 'Error: Trilinos test requires 1 or 2 processors'
       stop
    endif

  end subroutine trilinos_test

#endif

!****************************************************************************

  end module glissade_velo_higher_trilinos
  
!****************************************************************************
