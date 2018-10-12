module Block_Solve_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  private
  
  public :: bl3dfac, &
            bl3dsolf, &
            bl3dsolb
      
contains

! ************************************************************************** !

subroutine bl3dfac(n, k, E, D, F, pivot)

!	This version of bl3dfac dated March 15, 2000.

! Modified from bl3dfac.f, bl3dsol.f - Fortran codes of
! P. Keast, http://www.mscs.dal.ca/ keast/research/pubs.html, 
! Mathematics and Statistics Department, Dalhousie University, Canada.

!************************************************************************

!	bl3dfac performs an LU factorization of the block tridiagonal
!	matrix given by:

!	D(1) E(1)  O   O    O ............  O
!	F(1) D(2) E(2) O    O ............  O
!	 O   F(2) D(3) E(3) O ............  O
!	 ....................................
!
!        O    O   ........F(n-1) D(n-1) E(n-1)
!        O    O   .............. F(n-1)   D(n)

!	where each block, E(j), D(j), F(j) is k by k.
!
!	Pivoting is done only within the diagonal blocks.

!	The L matrix is NOT unit lower triangular, but the U matrix
!	is unit UPPER triangular, with the identity matrix on the
!	diagonal. 

!************************************************************************

!	Uses: Lapack routines dgetrf and dgetrs for Gauss LU factorization
!	and solution, and the BLAS routine dgemm for matrix-matrix 
!       multiplication.

!************************************************************************

!	Input variables:

  PetscInt :: n, k
  PetscInt :: pivot(k,n)
  PetscReal :: E(k,k,n-1), D(k,k,n), F(k,k,n-1)

! Input:

! k           Integer, the size of the blocks.
! n           Integer, the number of blocks.
! E(k,k,n-1)  Double precision, the sub-diagonal blocks.
! D(k,k,n)    Double precision, the diagonal blocks.
! F(k,k,n-1)  Double precision, the super-diagonal blocks.

! Output:

! E(k,k,n-1)  The lower triangular blocks of L in the LU factorization.
! F(k,k,n-1)  The upper triangular blocks of U in the LU factorization.
! D(k,k,n)    The diagonal blocks of L in the LU factorization.
!                   The diagonal blocks of U are unit upper triangular.
! pivot(k,n)  Integer pivot vector used in factoring the diagonal blocks.
! 	    pivot(1:k,p) records the pivoting done in diagonal block p.

! Local variables:

  PetscInt :: j, i, l, ll, info, info1, iwork(k)
  character(len=1) :: trans, norm = '1'
  PetscReal, parameter :: one = 1.d0
  PetscReal anorm, sum, rcond, work(4*k)

!************************************************************************

! Executable statements.

!************************************************************************

! First, in main blocks 1 to n-1:

  trans = 'N'
        
  do j = 1, n-1

!  compute 1-norm: only needed for computing condition nr.
#ifdef CONDNR
    anorm = 0.d0
    do l = 1, k
      sum = 0.d0
      do i = 1, k
        sum = sum + abs(d(i,l,j))
      enddo
      anorm = max(anorm,sum)
    enddo
#endif

!  First, factor D(j).
    call dgetrf(k, k, D(1,1,j), k, pivot(1,j), info)

    if ( info .ne. 0 ) then
!        write(*,'("node ",3i3,1p15e12.4)') j,n,k,(d(i,ll,j),ll=1,k)

      print *,'At block ',j,',',' 1-norm = ',anorm
      print *,'problem, info = ', info   ! make this better later.
      return
    endif

!  Estimate condition number
#ifdef CONDNR
    CALL dgecon(norm,k,D(1,1,j),k,anorm,rcond,work,iwork,info1)
    rcond = 1.e0/rcond
!   if (rcond > 1.e10)
        WRITE (*,999) 'Estimate of condition number =', &
        rcond,info1
    999 FORMAT (1X,A,1P,e11.4,i3)
#endif

!  Now, compute E(j) from D(j) * E(j) = E(j).
    call dgetrs(trans, k, k, D(1,1,j), k, pivot(1,j), &
                       E(1,1,j), k, info)
    if ( info .lt. 0 ) then
      print *,'Illegal parameter number ',-info
      return
    endif

!  Finally, compute D(j+1) = D(j+1) - F(j) * E(j).
    call dgemm(trans, trans, k, k, k, -one, F(1,1,j), k, &
                      E(1,1,j), k, one, D(1,1,j+1), k)
  enddo

! Finally, obtain the LU factorization of D(n).
  call dgetrf(k, k, D(1,1,n), k, pivot(1,n), info)
  if ( info .ne. 0 ) then
     print *,'At block ',j,','
     print *,'problem, info = ', info   ! make this better later.
    return
  endif

  return
  
end subroutine bl3dfac

! ************************************************************************** !

subroutine bl3dsol(n, k, E, D, F, pivot, nrhs, rhs)

!       This version of bl3dsol dated April 20, 2000.

!************************************************************************

!       bl3dsol solves the system of linear equations whose coefficient
!       matrix is the block tridiagonal matrix given by:

!       D(1) E(1)  O   O    O ............  O
!       F(1) D(2) E(2) O    O ............  O
!        O   F(2) D(3) E(3) O ............  O
!        ....................................
!
!        O    O   ........F(n-1) D(n-1) E(n-1)
!        O    O   .............. F(n-1)   D(n)

!       where each block, E(j), D(j), F(j) is k by k, and which has already
!       been factored using bl3dfac.

!       The right hand sides are stored by blocks (see below).

!       The L matrix is NOT unit lower triangular, but the U matrix
!       is unit UPPER triangular.

!************************************************************************

!       Uses: 

!************************************************************************

!       Input variables:

        character(len=1) :: trans
        PetscInt :: n, k
        PetscReal :: E(k,k,n-1), D(k,k,n), F(k,k,n-1)
        PetscInt :: pivot(k,n), nrhs
        PetscReal :: rhs(k,nrhs,n)

!       Input:

!       trans         character*1, 'n' or 'N' is A.x = b is to be solved,
!                     't' or 'T' is A^T.x = b is to be solved.
!       n             Integer, the number of blocks.
!       k             Integer, the size of the blocks.
!       E(k,k,n-1)    Double precision, the sub-diagonal blocks.
!       D(k,k,n)      Double precision, the diagonal blocks, after factoring 
!                     by bl3dfac.
!       F(k,k,n-1)    Double precision, the super-diagonal blocks, after 
!                     modification in bl3dfac.
!       pivot(k,n)    Integer pivot vector returned by bl3dfac, used in 
!                     factoring the diagonal blocks.
!                     pivot(1:k,p) records the pivoting done in diagonal 
!                     block p.
!       nrhs          integer, the number of right hand sides.
!       rhs(k,nrhs,n) Double precision, right hand sides. These are stored
!                     by blocks - that is, the nrhs right hand sides 
!                     associated with each k by k block are stored 
!                     consecutively. For example, for n = 3, k = 2 and
!                     nrhs = 3, the elements of the right hand sides
!                     are stored in the following order:

!                     [1 ] [3 ] [5 ]
!                     [2 ] [4 ] [6 ]  First block.

!                     [7 ] [9 ] [11]
!                     [8 ] [10] [12]  Second block.

!                     [13] [15] [17]  Third block.
!                     [14] [16] [18]

!       Output:

!       rhs(k,nrhs,n) double precision, the solutions.
!       
!       Local variables:

        PetscInt :: j, info
        PetscReal, parameter :: one = 1.d0

!************************************************************************

!       Executable statements.

!************************************************************************
        
!       Modification of the right hand sides, that is, solve the
!       block bi-diagonal lower triangular system L.x = b, overwriting
!       b with the solution.

        trans = 'N'


        if ( trans .eq. 'n' .or. trans .eq. 'N' ) then 

!          First, solve D(1).x(1..k) = b(1..k).
           call dgetrs(trans, k, nrhs, D(1,1,1), k, pivot(1,1), &
                       rhs(1,1,1), k, info)

           if ( info .lt. 0 ) then
              print *,'Illegal parameter number ',-info
              return
           endif

!          Then, for j = 2 to n solve 
!          D(j)x((j-1)*k+1..j*k) = b((j-1)*k+1..j*k) 
!                                - F(j-1)x((j-2)*k+1..(j-1)*k)

           do j = 2, n
           
              call dgemm(trans, trans, k, nrhs, k, -one, F(1,1,j-1), k, &
                         rhs(1,1,j-1), k, one, rhs(1,1,j), k)

              call dgetrs(trans, k, nrhs, D(1,1,j), k, pivot(1,j), &
                          rhs(1,1,j), k, info)
              if ( info .lt. 0 ) then
                 print *,'Illegal parameter number ',-info
                 return
              endif

           enddo

!        Now, the back substitution.

           do j = n-1,1,-1
   
!             Form rhs(:,:,j) = rhs(:,:,j) - E(j)*rhs(:,:,j+1)

              call dgemm(trans, trans, k, nrhs, k, -one, E(1,1,j), k, &
                         rhs(1,1,j+1), k, one, rhs(1,1,j), k)

           enddo

     else

!    Solve the transposed system.

!    First, the right hand side modification.

         do j = 2,n

!        Form rhs(:,:,j) = rhs(:,:,j) - E(j-1)^T*rhs(:,:,j-1)
            call dgemm(trans, 'N', k, nrhs, k, -one, E(1,1,j-1), k, &
                           rhs(1,1,j-1), k, one, rhs(1,1,j), k)
         enddo

!      Now, the back substitution.

!      First, solve D(n)^T.rhs(:,:,n) = rhs(:,:,n)

         call dgetrs(trans, k, nrhs, D(1,1,n), k, pivot(1,n), &
                     rhs(1,1,n), k, info)
         if ( info .lt. 0 ) then
            print *,'Illegal parameter number ',-info
            return
         endif

!          For j = n-1 down to 1, solve 
!          D(j)^T.rhs(:,:,j) = rhs(:,:,j) - F(j)^T.rhs(:,:,j+1)

         do j = n-1,1,-1

!        Form rhs(:,:,j) = rhs(:,:,j) - F(j)^T*rhs(:,:,j+1)
            call dgemm(trans, 'N', k, nrhs, k, -one, F(1,1,j), k, &
                           rhs(1,1,j+1), k, one, rhs(1,1,j), k)

            call dgetrs(trans, k, nrhs, D(1,1,j), k, pivot(1,j), &
                        rhs(1,1,j), k, info)
            if ( info .lt. 0 ) then
               print *,'Illegal parameter number ',-info
               return
            endif
         enddo
      endif

      return
      
end subroutine bl3dsol

! ************************************************************************** !

subroutine bl3dsolf(n, k, E, D, F, pivot, nrhs, rhs)

!       This version of bl3dsol dated April 20, 2000.

!************************************************************************

!       bl3dsol solves the system of linear equations whose coefficient
!       matrix is the block tridiagonal matrix given by:

!       D(1) E(1)  O   O    O ............  O
!       F(1) D(2) E(2) O    O ............  O
!        O   F(2) D(3) E(3) O ............  O
!        ....................................
!
!        O    O   ........F(n-1) D(n-1) E(n-1)
!        O    O   .............. F(n-1)   D(n)

!       where each block, E(j), D(j), F(j) is k by k, and which has already
!       been factored using bl3dfac.

!       The right hand sides are stored by blocks (see below).

!       The L matrix is NOT unit lower triangular, but the U matrix
!       is unit UPPER triangular.

!************************************************************************

!       Uses: 

!************************************************************************

!       Input variables:

        character(len=1) :: trans
        PetscInt :: n, k
        PetscReal :: E(k,k,n-1), D(k,k,n), F(k,k,n-1)
        PetscInt :: pivot(k,n), nrhs
        PetscReal :: rhs(k,nrhs,n)

!       Input:

!       trans         character*1, 'n' or 'N' is A.x = b is to be solved,
!                     't' or 'T' is A^T.x = b is to be solved.
!       n             Integer, the number of blocks.
!       k             Integer, the size of the blocks.
!       E(k,k,n-1)    Double precision, the sub-diagonal blocks.
!       D(k,k,n)      Double precision, the diagonal blocks, after factoring 
!                     by bl3dfac.
!       F(k,k,n-1)    Double precision, the super-diagonal blocks, after 
!                     modification in bl3dfac.
!       pivot(k,n)    Integer pivot vector returned by bl3dfac, used in 
!                     factoring the diagonal blocks.
!                     pivot(1:k,p) records the pivoting done in diagonal 
!                     block p.
!       nrhs          integer, the number of right hand sides.
!       rhs(k,nrhs,n) Double precision, right hand sides. These are stored
!                     by blocks - that is, the nrhs right hand sides 
!                     associated with each k by k block are stored 
!                     consecutively. For example, for n = 3, k = 2 and
!                     nrhs = 3, the elements of the right hand sides
!                     are stored in the following order:

!                     [1 ] [3 ] [5 ]
!                     [2 ] [4 ] [6 ]  First block.

!                     [7 ] [9 ] [11]
!                     [8 ] [10] [12]  Second block.

!                     [13] [15] [17]  Third block.
!                     [14] [16] [18]

!       Output:

!       rhs(k,nrhs,n) double precision, the solutions.
!       
!       Local variables:

        PetscInt :: j, info
        PetscReal, parameter :: one = 1.d0

!************************************************************************

!       Executable statements.

!************************************************************************
        
!       Modification of the right hand sides, that is, solve the
!       block bi-diagonal lower triangular system L.x = b, overwriting
!       b with the solution.

          trans = 'N'


!       if ( trans .eq. 'n' .or. trans .eq. 'N' ) then 

!          First, solve D(1).x(1..k) = b(1..k).
           call dgetrs(trans, k, nrhs, D(1,1,1), k, pivot(1,1), &
                       rhs(1,1,1), k, info)

           if ( info .lt. 0 ) then
              print *,'Illegal parameter number ',-info
              return
           endif

!          Then, for j = 2 to n solve 
!          D(j)x((j-1)*k+1..j*k) = b((j-1)*k+1..j*k) 
!                                - F(j-1)x((j-2)*k+1..(j-1)*k)

           do j = 2, n
           
              call dgemm(trans, trans, k, nrhs, k, -one, F(1,1,j-1), k, &
                         rhs(1,1,j-1), k, one, rhs(1,1,j), k)

              call dgetrs(trans, k, nrhs, D(1,1,j), k, pivot(1,j), &
                          rhs(1,1,j), k, info)
              if ( info .lt. 0 ) then
                 print *,'Illegal parameter number ',-info
                 return
              endif

           enddo

      return

end subroutine bl3dsolf

! ************************************************************************** !

subroutine bl3dsolb(n, k, E, D, F, pivot, nrhs, rhs)

!       This version of bl3dsol dated April 20, 2000.

!************************************************************************

!       bl3dsol solves the system of linear equations whose coefficient
!       matrix is the block tridiagonal matrix given by:

!       D(1) E(1)  O   O    O ............  O
!       F(1) D(2) E(2) O    O ............  O
!        O   F(2) D(3) E(3) O ............  O
!        ....................................
!
!        O    O   ........F(n-1) D(n-1) E(n-1)
!        O    O   .............. F(n-1)   D(n)

!       where each block, E(j), D(j), F(j) is k by k, and which has already
!       been factored using bl3dfac.

!       The right hand sides are stored by blocks (see below).

!       The L matrix is NOT unit lower triangular, but the U matrix
!       is unit UPPER triangular.

!************************************************************************

!       Uses: 

!************************************************************************

!       Input variables:

        character(len=1) :: trans
        PetscInt :: n, k
        PetscReal :: E(k,k,n-1), D(k,k,n), F(k,k,n-1)
        PetscInt :: pivot(k,n), nrhs
        PetscReal :: rhs(k,nrhs,n)

!       Input:

!       trans         character*1, 'n' or 'N' is A.x = b is to be solved,
!                     't' or 'T' is A^T.x = b is to be solved.
!       n             Integer, the number of blocks.
!       k             Integer, the size of the blocks.
!       E(k,k,n-1)    Double precision, the sub-diagonal blocks.
!       D(k,k,n)      Double precision, the diagonal blocks, after factoring 
!                     by bl3dfac.
!       F(k,k,n-1)    Double precision, the super-diagonal blocks, after 
!                     modification in bl3dfac.
!       pivot(k,n)    Integer pivot vector returned by bl3dfac, used in 
!                     factoring the diagonal blocks.
!                     pivot(1:k,p) records the pivoting done in diagonal 
!                     block p.
!       nrhs          integer, the number of right hand sides.
!       rhs(k,nrhs,n) Double precision, right hand sides. These are stored
!                     by blocks - that is, the nrhs right hand sides 
!                     associated with each k by k block are stored 
!                     consecutively. For example, for n = 3, k = 2 and
!                     nrhs = 3, the elements of the right hand sides
!                     are stored in the following order:

!                     [1 ] [3 ] [5 ]
!                     [2 ] [4 ] [6 ]  First block.

!                     [7 ] [9 ] [11]
!                     [8 ] [10] [12]  Second block.

!                     [13] [15] [17]  Third block.
!                     [14] [16] [18]

!       Output:

!       rhs(k,nrhs,n) double precision, the solutions.
!       
!       Local variables:

        PetscInt :: j !, info
        PetscReal, parameter :: one = 1.d0

        trans = 'N'

!        Now, the back substitution.

           do j = n-1,1,-1
   
!             Form rhs(:,:,j) = rhs(:,:,j) - E(j)*rhs(:,:,j+1)

              call dgemm(trans, trans, k, nrhs, k, -one, E(1,1,j), k, &
                         rhs(1,1,j+1), k, one, rhs(1,1,j), k)

           enddo


      return

end subroutine bl3dsolb

!************************ END OF bl3dsol ********************************
end module Block_Solve_module
  
