module Block_Tridiag_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Utility_module, only : Equal
  private

  public :: decbt, &
            solbtf, &
            solbtb
     
contains

! ************************************************************************** !

subroutine decbt (m, n, ndim, a, b, c, ip, ier)
      
  use PFLOTRAN_Constants_module

  implicit none

  PetscInt :: m,n,ndim,ip(ndim,n),ier
  PetscInt :: nm1,nm2,km1,i,j,k,l
  PetscReal :: dp
  PetscReal :: a(ndim,ndim,n), b(ndim,ndim,n), c(ndim,ndim,n)
     
!--------------------------------------
! block-tridiagonal matrix decomposition routine.
! written by a. c. hindmarsh.
! latest revision january 26, 1977  (ag)
! reference]  ucid-30150
!             solution of block-tridiagonal systems of linear
!             algebraic equations
!             a.c. hindmarsh
!             february 1977
! the input matrix contains three blocks of elements in each block-row,
! including blocks in the (1,3) and (n,n-2) block positions.
! decbt uses block gauss elimination and subroutines dec and sol
! for solution of blocks.  partial pivoting is done within
! block-rows only.
! input..
!     m = order of each block.
!     n = number of blocks in each direction of the matrix.
!         n must be 4 or more.  the complete matrix has order m*n.
!     a = m by m by n array containing diagonal blocks.
!         a(i,j,k) contains the (i,j) element of the k-th block.
!     b = m by m by n array containing the super-diagonal blocks
!         (in b(*,*,k) for k = 1,...,n-1) and the block in the (n,n-2)
!         block position (in b(*,*,n)).
!     c = m by m by n array containing the subdiagonal blocks
!         (in c(*,*,k) for k = 2,3,...,n) and the block in the
!         (1,3) block position (in c(*,*,1)).
!    ip = integer array of length m*n for working storage.
! output..
! a,b,c = m by m by n arrays containing the block lu decomposition
!         of the input matrix.
!    ip = m by n array of pivot information.  ip(*,k) contains
!         information for the k-th digonal block.
!   ier = 0  if no trouble occurred, or
!       = -1 if the input value of m or n was illegal, or
!       = k  if a singular matrix was found in the k-th diagonal block.
! use solbt to solve the associated linear system.
! decbt calls subroutines  dec(m,m0,a,ip,ier)  and  sol(m,m0,a,y,ip)
! for solution of m by m linear systems.
!------------------------------------------------------------------
      
  if (m .lt. 1 .or. n .lt. 4) goto 210
    nm1 = n - 1
    nm2 = n - 2

! process the first block-row. ------------------------------
    call dec (m, ndim, a(:,:,1), ip, ier)
    k = 1
    if (ier .ne. 0) goto 200
    do j = 1,m
      call sol (m, ndim, a(:,:,1), b(1,j,1), ip)
      call sol (m, ndim, a(:,:,1), c(1,j,1), ip)
    enddo
! adjust b(*,*,2). -----------------------------------------------
    do j = 1,m
      do i = 1,m
        dp = 0.d0
        do l = 1,m
          dp = dp + c(i,l,2)*c(l,j,1)
        enddo
        b(i,j,2) = b(i,j,2) - dp
      enddo
    enddo
! main loop.  process block-rows 2 to n-1. --------------
    do k = 2,nm1
      km1 = k - 1
      do j = 1,m
        do i = 1,m
          dp = 0.d0
          do l = 1,m
            dp = dp + c(i,l,k)*b(l,j,km1)
          enddo
          a(i,j,k) = a(i,j,k) - dp
        enddo
      enddo

      call dec (m, ndim, a(1,1,k), ip(1,k), ier)

      if (ier .ne. 0) goto 200
      do j = 1,m
        call sol (m, ndim, a(1,1,k), b(1,j,k), ip(1,k))
      enddo
    enddo

! process last block-row and return. ------------------
    do j = 1,m
      do i = 1,m
        dp = 0.d0
        do l = 1,m
          dp = dp + b(i,l,n)*b(l,j,nm2)
        enddo
        c(i,j,n) = c(i,j,n) - dp
      enddo
    enddo
    do j = 1,m
      do i = 1,m
        dp = 0.d0
        do l = 1,m
          dp = dp + c(i,l,n)*b(l,j,nm1)
        enddo
        a(i,j,n) = a(i,j,n) - dp
      enddo
    enddo

    call dec (m, ndim, a(1,1,n), ip(1,n), ier)
    k = n
    if (ier .ne. 0) goto 200
    return

! error returns. ------------------------------------
200  ier = k
     return

210  ier = -1
     return

end subroutine decbt

! ************************************************************************** !

subroutine solbt (m, n, ndim, a, b, c, ip, y)

  implicit none

  PetscInt :: m, n, ndim, ip(ndim,n)
  PetscInt :: nm1, nm2, km1, i, j, k, kb, kp1
  PetscReal :: a(ndim,ndim,n),b(ndim,ndim,n),c(ndim,ndim,n)
  PetscReal :: y(ndim,n),dp

!---------------------
! solution of block-tridiagonal linear system.
! coefficient matrix must have been previously processed by decbt.
! m, n, a, b, c, and ip  must not have been changed since call to decbt
! written by a. c. hindmarsh.
! input..
!     m = order of each block.
!     n = number of blocks in each direction of matrix.
! a,b,c = m by m by n arrays containing block lu decomposition
!         of coefficient matrix from decbt.
!    ip = m by n integer array of pivot information from decbt.
!     y = array of length m*n containg the right-hand side vector
!         (treated as an m by n array here).
! output..
!     y = solution vector, of length m*n.
! solbt makes calls to subroutine sol(m,m0,a,y,ip)
! for solution of m by m linear systems.
!----------------------------------------------------------------

  nm1 = n - 1
  nm2 = n - 2

! forward solution sweep. ---------------------------------
  call sol (m, ndim, a(:,:,1), y, ip)
  do k = 2,nm1
    km1 = k - 1
    do i = 1,m
      dp = 0.d0
      do j = 1,m
        dp = dp + c(i,j,k)*y(j,km1)
      enddo
      y(i,k) = y(i,k) - dp
    enddo
    call sol (m, ndim, a(1,1,k), y(1,k), ip(1,k))
  enddo

  do i = 1,m
    dp = 0.d0
    do j = 1,m
      dp = dp + c(i,j,n)*y(j,nm1) + b(i,j,n)*y(j,nm2)
    enddo
    y(i,n) = y(i,n) - dp
  enddo

  call sol (m, ndim, a(1,1,n), y(1,n), ip(1,n))

! backward solution sweep. ----------------------------------

!     do j=1,n
!       print *,'solbt: ',m,n,ndim,j,y(1,j),ip(1,j),a(1,1,j),b(1,1,j),c(1,1,j)
!     enddo
      
  do kb = 1,nm1
    k = n - kb
    kp1 = k + 1
    do i = 1,m
      dp = 0.d0
      do j = 1,m
        dp = dp + b(i,j,k)*y(j,kp1)
      enddo
      y(i,k) = y(i,k) - dp
    enddo
  enddo

  do i = 1,m
    dp = 0.d0
    do j = 1,m
      dp = dp + c(i,j,1)*y(j,3)
    enddo
    y(i,1) = y(i,1) - dp
  enddo

  return
end subroutine solbt

! ************************************************************************** !

subroutine solbtf (m, n, ndim, a, b, c, ip, y)

  implicit none
  PetscInt :: m, n, ndim, ip(ndim,n)
  PetscInt :: nm1, nm2, km1, i, j, k
  PetscReal :: a(ndim,ndim,n),b(ndim,ndim,n),c(ndim,ndim,n)
  PetscReal :: y(ndim,n),dp

!---------------------
! solution of block-tridiagonal linear system.
! coefficient matrix must have been previously processed by decbt.
! m, n, a, b, c, and ip  must not have been changed since call to decbt
! written by a. c. hindmarsh.
! input..
!     m = order of each block.
!     n = number of blocks in each direction of matrix.
! a,b,c = m by m by n arrays containing block lu decomposition
!         of coefficient matrix from decbt.
!    ip = m by n integer array of pivot information from decbt.
!     y = array of length m*n containg the right-hand side vector
!         (treated as an m by n array here). 
! output..
!     y = solution vector, of length m*n.
! solbt makes calls to subroutine sol(m,m0,a,y,ip)
! for solution of m by m linear systems.
!----------------------------------------------------------------

  nm1 = n - 1
  nm2 = n - 2

! forward solution sweep. ---------------------------------
  call sol (m, ndim, a, y, ip)
  do k = 2,nm1
    km1 = k - 1
    do i = 1,m
      dp = 0.d0
      do j = 1,m
        dp = dp + c(i,j,k)*y(j,km1)
      enddo
      y(i,k) = y(i,k) - dp
    enddo
    call sol (m, ndim, a(1,1,k), y(1,k), ip(1,k))
  enddo

  do i = 1,m
    dp = 0.d0
    do j = 1,m
      dp = dp + c(i,j,n)*y(j,nm1) + b(i,j,n)*y(j,nm2)
    enddo
    y(i,n) = y(i,n) - dp
  enddo

  call sol (m, ndim, a(1,1,n), y(1,n), ip(1,n))

  return
end subroutine solbtf

! ************************************************************************** !

subroutine solbtb (m, n, ndim, a, b, c, ip, y)

  implicit none

  PetscInt :: m, n, ndim, ip(ndim,n)
  PetscInt :: nm1, i, j, k, kb, kp1
  PetscReal :: a(ndim,ndim,n),b(ndim,ndim,n),c(ndim,ndim,n)
  PetscReal :: y(ndim,n),dp

!---------------------
! solution of block-tridiagonal linear system.
! coefficient matrix must have been previously processed by decbt.
! m, n, a, b, c, and ip  must not have been changed since call to decbt
! written by a. c. hindmarsh.
! input..
!     m = order of each block.
!     n = number of blocks in each direction of matrix.
! a,b,c = m by m by n arrays containing block lu decomposition
!         of coefficient matrix from decbt.
!    ip = m by n integer array of pivot information from decbt.
!     y = array of length m*n containg the right-hand side vector
!         (treated as an m by n array here).
! output..
!     y = solution vector, of length m*n.
! solbt makes calls to subroutine sol(m,m0,a,y,ip)
! for solution of m by m linear systems.
!----------------------------------------------------------------

!     call sol (m, ndim, a(1,1,n), y(1,n), ip(1,n))

  nm1 = n - 1

!     do j=1,n
!       print *,'solbtb: ',m,n,ndim,j,y(1,j),ip(1,j),a(1,1,j),b(1,1,j),c(1,1,j)
!     enddo

! backward solution sweep. ----------------------------------
  do kb = 1,nm1
    k = n - kb
    kp1 = k + 1
    do i = 1,m
      dp = 0.d0
      do j = 1,m
        dp = dp + b(i,j,k)*y(j,kp1)
      enddo
      y(i,k) = y(i,k) - dp
    enddo
  enddo

  do i = 1,m
    dp = 0.d0
    do j = 1,m
      dp = dp + c(i,j,1)*y(j,3)
    enddo
    y(i,1) = y(i,1) - dp
  enddo

  return
end subroutine solbtb

! ************************************************************************** !

subroutine dec (n, ndim, a, ip, ier)

  implicit none

  PetscInt :: m, n, ndim, ip(n), ier
  PetscInt :: nm1,i,j,k,kp1
  PetscReal :: a(ndim,n), t

!---------------
!  matrix triangularization by gauss elimination with partial pivoting.
!  input..
!     n = order of matrix.
!     ndim = declared first dimension of array  a.
!     a = matrix to be triangularized.
!  output..
!     a(i,j), i.le.j = upper triangular factor, u .
!     a(i,j), i.gt.j = multipliers = lower triangular factor, i - l.
!     ip(k), k.lt.n = index of k-th pivot row.
!     ier = 0 if matrix a is nonsingular, or k if found to be
!           singular at stage k.
!  row interchanges are finished in u, only partly in l.
!  use  sol  to obtain solution of linear system.
!  if ier .ne. 0, a is singular, sol will divide by 0.d0.
!------------------------------------------------------------
!
!  reference:  a. c. hindmarsh, l. j. sloan, k. w. fong, and
!              g. h. rodrigue,
!              dec/sol:  solution of dense systems of linear
!              algebraic equations,
!              lawrence livermore laboratory report ucid-30137,
!              june 1976.
!----------------------------------------------------------------

  ier = 0
  if (n .eq. 1) goto 70
  nm1 = n - 1
  do k = 1,nm1
    kp1 = k + 1

!  find the pivot in column k.  search rows k to n. -------
    m = k
    do i = kp1,n
      if (abs(a(i,k)) .gt. abs(a(m,k))) m = i
    enddo
    ip(k) = m

!  interchange elements in rows k and m. -------------
    t = a(m,k)
    if (m .eq. k) goto 20
    a(m,k) = a(k,k)
    a(k,k) = t
 20     if (Equal(t,0.d0)) goto 80

!  store multipliers in a(i,k), i = k+1,...,n. -
    t = 1.d0/t
    do i = kp1,n
      a(i,k) = -a(i,k)*t
    enddo

!  apply multipliers to other columns of a.
    do j = kp1,n
      t = a(m,j)
      a(m,j) = a(k,j)
      a(k,j) = t
      if (.not. Equal(t,0.d0)) then
        do i = kp1,n
          a(i,j) = a(i,j) + a(i,k)*t
        enddo
      endif
    enddo
  enddo

70   k = n
     if (Equal(a(n,n),0.d0)) goto 80
     return

80   ier = k
     return

end subroutine dec

! ************************************************************************** !

subroutine sol (n, ndim, a, b, ip)

  implicit none

  PetscInt :: m, n, ndim, ip(n)
  PetscInt :: nm1,i,k,kb,km1,kp1
  PetscReal :: a(ndim,n), b(n), t

!--------------------------------------------------
!  solution of linear system a*x = b using output of dec.
!  input..
!     n = order of matrix.
!     ndim = declared first dimension of array  a.
!     a = triangularized matrix obtained from dec.
!     b = right hand side vector.
!     ip = pivot information vector obtained from dec.
!  do not use if dec has set ier .ne. 0.
!  output..
!     b = solution vector, x .
!--------------------------------------------------------------
!
!  reference:  a. c. hindmarsh, l. j. sloan, k. w. fong, and
!              g. h. rodrigue,
!              dec/sol:  solution of dense systems of linear
!              algebraic equations,
!              lawrence livermore laboratory report ucid-30137,
!              june 1976.
!----------------------------------------------------------------

  if (n .eq. 1) goto 50
  nm1 = n - 1

!  apply row permutations and multipliers to b. --------------
  do k = 1,nm1
    kp1 = k + 1
    m = ip(k)
    t = b(m)
    b(m) = b(k)
    b(k) = t
    do i = kp1,n
      b(i) = b(i) + a(i,k)*t
    enddo
   enddo

!  back solve. -------------------------------------
   do kb = 1,nm1
    km1 = n - kb
    k = km1 + 1
    b(k) = b(k)/a(k,k)
    t = -b(k)
    do i = 1,km1
      b(i) = b(i) + a(i,k)*t
    enddo
   enddo

50   b(1) = b(1)/a(1,1)

  return
end subroutine sol

end module Block_Tridiag_module
