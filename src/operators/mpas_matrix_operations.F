! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_matrix_operations
!
!> \brief MPAS matrix operations
!> \author Mark Petersen
!> \date    April 2013
!> \details
!>  This module contains the routines for matrix operations
!
!-----------------------------------------------------------------------

module mpas_matrix_operations

   use mpas_grid_types
   use mpas_constants

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: mpas_rotate_2D_matrix_2x2, &
             mpas_rotate_2D_matrix_sym3index, &
             mpas_matrix_sym6index_to_3x3, &
             mpas_matrix_3x3_to_sym6index, &
             mpas_matrix_cell_to_edge, &
             mpas_outer_product

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine mpas_rotate_2D_matrix_2x2
!
!> \brief   Rotate a 2D matrix in 2x2 index format
!> \author  Mark Petersen
!> \date    April 2013
!> \details 
!>  Rotate matrix A about the angle theta using rotation matrix R, where
!>  B = RAR'  and R' is the tranpose.
!
!-----------------------------------------------------------------------

   subroutine mpas_rotate_2D_matrix_2x2(A,theta,B)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(2,2), intent(in) :: &
         A   !< Input: 2x2 matrix

      real (kind=RKIND), intent(in) :: &
         theta   !< Input: rotation angle

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(2,2), intent(out) :: &
         B   !< Output: 2x2 matrix

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(2,2) :: &
         R   ! 2x2 rotation matrix

      integer :: i,j,k,l

      real (kind=RKIND) :: sin1, cos1

      cos1 = cos(theta)
      sin1 = sin(theta)

      ! create rotation matrix
      R(1,1) = cos1
      R(1,2) = -sin1
      R(2,1) = sin1
      R(2,2) = cos1

      B = 0.0
      do i=1,2
        do j=1,2
          do k=1,2
            do l=1,2
              ! B = R A R' in matrix format
              ! B(i,j) =  R(i,k) A(k,l) R(j,l) in index format
              ! note indices switched on second R due to transpose.
              B(i,j) = B(i,j) + R(i,k)*A(k,l)*R(j,l)
            enddo
          enddo
        enddo
      enddo

   end subroutine mpas_rotate_2D_matrix_2x2!}}}

!***********************************************************************
!
!  routine mpas_rotate_2D_matrix_sym3index
!
!> \brief   Rotate a 2D matrix in 3-index format
!> \author  Mark Petersen
!> \date    April 2013
!> \details 
!>  Rotate matrix A about the angle theta using rotation matrix R, where
!>  B = RAR'  and R' is the tranpose.  A is a symmetric matrix in 3-index
!>  format, where A(1)=A_{11}, A(2)=A_{22}, and A(3)=A_{12}
!
!-----------------------------------------------------------------------

   subroutine mpas_rotate_2D_matrix_sym3index(A,theta,B)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(3), intent(in) :: &
         A   !< Input: 2x2 symmetric matrix in 3-index format

      real (kind=RKIND), intent(in) :: &
         theta   !< Input: rotation angle

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(3), intent(out) :: &
         B   !< Input: 2x2 symmetric matrix in 3-index format

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND) :: sin1, cos1, sin2, cos2, cossin

      cos1 = cos(theta)
      sin1 = sin(theta)
      cos2 = cos1**2
      sin2 = sin1**2
      cossin = cos1*sin1

      ! Tensor rotation through angle theta:
      ! R = | cos(theta)    -sin(theta)  |
      !     | sin(theta)     cos(theta)  |
      ! is the 2D rotation matrix.  Tensor rotation is
      ! B = R A R'  where R' is the transpose.

      ! A and B are symmetric, in 3-index format:
      ! A(1)=A_{11}, A(2)=A_{22}, and A(3)=A_{12}

      B(1) = A(1)*cos2 + A(2)*sin2 - A(3)*2*cossin
      B(2) = A(1)*sin2 + A(2)*cos2 + A(3)*2*cossin 
      B(3) = (A(1) - A(2))*cossin + A(3)*(cos2-sin2)

   end subroutine mpas_rotate_2D_matrix_sym3index!}}}


!***********************************************************************
!
!  routine mpas_matrix_sym6index_to_3x3
!
!> \brief   Convert a symetric 6-index matrix to 3x3 format
!> \author  Mark Petersen
!> \date    April 2013
!> \details 
!>  Given a symetric 6-index matrix A, produce a 3x3 matrix B, where
!>  A(1)=B_{11}, A(2)=B_{22}, A(3)=A_{33},
!>  A(4)=B_{12}, A(5)=B_{23}, A(6)=A_{13},
!
!-----------------------------------------------------------------------

   subroutine mpas_matrix_sym6index_to_3x3(A,B)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(6), intent(in) :: &
         A   !< Input: 3x3 symmetric matrix in 6-index format

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(3,3), intent(out) :: &
         B   !< Output: 3x3 symmetric matrix in 3x3 format


      B(1,1) = A(1)
      B(2,2) = A(2)
      B(3,3) = A(3)
      B(1,2) = A(4)
      B(2,1) = A(4)
      B(2,3) = A(5)
      B(3,2) = A(5)
      B(1,3) = A(6)
      B(3,1) = A(6)

   end subroutine mpas_matrix_sym6index_to_3x3!}}}


!***********************************************************************
!
!  routine mpas_matrix_3x3_to_sym6index
!
!> \brief   Convert a 3x3 format matrix to a symetric 6-index matrix
!> \author  Mark Petersen
!> \date    April 2013
!> \details 
!>  Given a 3x3 matrix B, produce a symetric 6-index matrix A,
!>  B(1)=A(1,1), B(2)=A(2,2), B(3)=A(3,3)
!>  B(4)=(A(1,2)+A(2,1))/2
!>  B(5)=(A(2,3)+A(3,2))/2
!>  B(6)=(A(1,3)+A(3,1))/2
!
!-----------------------------------------------------------------------

   subroutine mpas_matrix_3x3_to_sym6index(A,B)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(3,3), intent(in) :: &
         A   !< Input: 3x3 symmetric matrix in 3x3 format

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(6), intent(out) :: &
         B   !< Output: 3x3 symmetric matrix in 6-index format


      B(1) = A(1,1)
      B(2) = A(2,2)
      B(3) = A(3,3)
      B(4) = 0.5*(A(1,2)+A(2,1))
      B(5) = 0.5*(A(2,3)+A(3,2))
      B(6) = 0.5*(A(1,3)+A(3,1))

   end subroutine mpas_matrix_3x3_to_sym6index!}}}


!***********************************************************************
!
!  routine mpas_matrix_cell_to_edge
!
!> \brief   Interpolate a matrix from cell to edge
!> \author  Mark Petersen
!> \date    5 April 2013
!> \details 
!>  This routine interpolates a matrix from cell to edge locations.
!
!-----------------------------------------------------------------------

   subroutine mpas_matrix_cell_to_edge(grid, matrixCell, matrixEdge)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         matrixCell   !< Input: matrix located at Cell

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(out) :: &
         matrixEdge   !< Output: matrix located at Edge

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, nEdges, cell1, cell2, k, nVertLevels
      integer, dimension(:,:), pointer :: cellsOnEdge

      nEdges = grid % nEdges
      nVertLevels = grid % nVertLevels

      cellsOnEdge => grid % cellsOnEdge % array

      ! error check that index 1 of matrixEdge and matrixCell are same length?

      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,nVertLevels

            matrixEdge(:,k,iEdge) = 0.5*(matrixCell(:,k,cell1) + matrixCell(:,k,cell2))

         enddo
      enddo

   end subroutine mpas_matrix_cell_to_edge!}}}


!***********************************************************************
!
!  routine mpas_outer_product
!
!> \brief   Compute the outer product of vectors u and v
!> \author  Mark Petersen
!> \date    April 2013
!> \details 
!>  Given two n-length vector u and m-length vector v, compute the 
!>  outer product (tensor product) to compute the nxm matrix A,
!>  A(i,j) = u(i)*v(j)
!
!-----------------------------------------------------------------------

   subroutine mpas_outer_product(u,v,A)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
         u, &!< Input: n-length vector
         v   !< Input: m-length vector

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         A   !< Output: nxm matrix

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: n,m,i,j

      n = size(u)
      m = size(v)

      ! mrp: find out how to do correct error check.
      if (size(A,1).ne.n.or.size(A,2).ne.m) then
        write(*,*) 'error, mpas_outer_product: size of A must match u and v'
        stop
      endif

      do i=1,n
        do j=1,m
          A(i,j) = u(i)*v(j)
        enddo
      enddo

   end subroutine mpas_outer_product!}}}


end module mpas_matrix_operations

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker

