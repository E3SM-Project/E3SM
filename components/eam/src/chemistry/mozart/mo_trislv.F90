
      module mo_trislv

      use shr_kind_mod, only : r8 => shr_kind_r8

      contains

      subroutine tridec( syscnt, order, lower, main, upper )
!--------------------------------------------------------------------
! dimension of           lower(syscnt,order), main(syscnt,order), upper(syscnt,order)
! arguments
!
! latest revision        june 1995
!
! purpose                trdi computes the solution of the tridiagonal
!                        linear system,
!                            b(1)*x(1)+c(1)*x(2)               = y(1)
!                            a(i)*x(i-1)+b(i)*x(i)+c(i)*x(i+1) = y(i)
!                                i=2,3,...,n-1
!                            a(n)*x(n-1)+b(n)*x(n)             = y(n)
!
! usage                  call trdi (n, a, b, c, x )
!
! arguments
!
! on input               n
!                          the number of unknowns.  this subroutine
!                          requires that  n  be greater than  2.
!
!                        a
!                          the subdiagonal of the matrix is stored in
!                          locations a(2) through a(n).
!
!                        b
!                          the diagonal of the matrix is stored in
!                          locations b(1) through b(n).
!
!                        c
!                          the super-diagonal of the matrix is stored in
!                          locations c(1) through c(n-1).
!
!                        x
!                          the right-hand side of the equations is
!                          stored in y(1) through y(n).
!
! on output              x
!                          an array which contains the solution to the
!                          system of equations.
!
! special conditions     this subroutine executes satisfactorily
!                        if the input matrix is diagonally dominant
!                        and non-singular.  the diagonal elements
!                        are used to pivot, and no tests are made to
!                        determine singularity.  if a singular, or
!                        nearly singular, matrix is used as input,
!                        a divide by zero or floating point overflow
!                        may result.
!
! precision              compiler dependent
!
! language               fortran
!
! history                originally written by nancy werner at ncar
!                        in october, 1973. modified by s. walters at
!                        ncar in june 1989 to functionally replace
!                        the cal routine tridla.
!
! portability            fortran 90
!
! algorithm              an lu-decomposition is obtained using the
!                        algorithm described in the reference below.
!
!                        to avoid unnecessary divisions, the alpha
!                        values used in the routine are the
!                        reciprocals of the alpha values described
!                        in the reference below.
!
! accuracy               every component of the residual of the linear
!                        system (i.e. the difference between  y  and
!                        the matrix applied to x) should be less in
!                        magnitude than ten times the machine precision
!                        times the matrix order times the maximum
!                        absolute component of the solution vector
!                        times the largest absolute row sum of the
!                        input matrix.
!
! timing                 the timing is roughly proportional to the
!                        order n of the linear system.
!
! references             analysis of numerical methods by
!                        e. isaacson and h. keller
!                        (john wiley and sons, 1966) pp. 55-58.
!--------------------------------------------------------------------

      implicit none

!--------------------------------------------------------------------
!	... dummy args
!--------------------------------------------------------------------
      integer, intent(in) ::  syscnt, order
      real(r8), intent(in), dimension(syscnt,order)    ::  lower
      real(r8), intent(inout), dimension(syscnt,order) ::  main, upper
      
!--------------------------------------------------------------------
!	... local variables
!--------------------------------------------------------------------
      integer :: i
      
!----------------------------------------------------------------------
!     	... lu-decomposition
!----------------------------------------------------------------------
      main(:,1)  = 1._r8 / main(:,1)
      upper(:,1) = upper(:,1)*main(:,1)
      do i = 2,order-1
         main(:,i)  = 1._r8 / (main(:,i) - lower(:,i)*upper(:,i-1))
         upper(:,i) = upper(:,i)*main(:,i)
      end do

      end subroutine tridec

      subroutine trislv( syscnt, order, lower, main, upper, x )

      implicit none

!----------------------------------------------------------------------
!     	... dummy args
!----------------------------------------------------------------------
      integer, intent(in) ::  syscnt, order
      real(r8), intent(in), dimension(syscnt,order)    ::  lower, &
                                                       main, &
                                                       upper
      real(r8), intent(inout), dimension(syscnt,order) ::  x
      
!----------------------------------------------------------------------
!     	... local variables
!----------------------------------------------------------------------
      integer :: i, im1, j, n, nm1

      nm1 = order - 1
      n   = order
!----------------------------------------------------------------------
!     	... solve the system
!----------------------------------------------------------------------
      x(:,1) = x(:,1)*main(:,1)
      do i = 2,nm1
         x(:,i) = (x(:,i) - lower(:,i)*x(:,i-1))*main(:,i)
      end do

      x(:,n) = (x(:,n) - lower(:,n)*x(:,nm1))/(main(:,n) - lower(:,n)*upper(:,nm1))
      do i = nm1,1,-1
         x(:,i) = x(:,i) - upper(:,i)*x(:,i+1)
      end do

      end subroutine trislv

      end module mo_trislv
