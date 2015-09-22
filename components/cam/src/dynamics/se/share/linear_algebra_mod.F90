#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module linear_algebra_mod
  implicit none
  private
  public :: matvec
  public :: reigen_solver

contains

  ! ======================================
  ! matvec : y = Ax
  ! Hand unrolled matrix vector multiply
  ! =====================================
  subroutine matvec(x,y,A,npts)
    use kinds, only: r8 => real_kind
    integer, intent(in) :: npts
    real (r8), intent(in) :: x(npts)
    real (r8), intent(out) :: y(npts)
    real (r8), intent(in) :: A(npts,npts)

    ! Local variables

    integer i,j,npts4
    real (r8) :: tmp0,tmp1,tmp2,tmp3

    npts4 = npts - MODULO(npts,4)

    do j=1,npts,4
       tmp0=0.0D0
       tmp1=0.0D0
       tmp2=0.0D0
       tmp3=0.0D0
       do i=1,npts
          tmp0=tmp0+A(i,j  )*x(i)
          tmp1=tmp1+A(i,j+1)*x(i)
          tmp2=tmp2+A(i,j+2)*x(i)
          tmp3=tmp3+A(i,j+3)*x(i)
       end do
       y(j  ) = tmp0
       y(j+1) = tmp1
       y(j+2) = tmp2
       y(j+3) = tmp3
    end do

    ! ============================
    ! remainder code
    ! ============================

    do j=npts4+1,npts
       tmp0=0.0D0
       do i=1,npts
          tmp0=tmp0+A(i,j  )*x(i)
       end do
       y(j  ) = tmp0
    end do

  end subroutine matvec

  ! =====================================
  ! Solve for AE = lambda*E
  ! F90 interface to LAPACK DGEEV
  ! =====================================

  function reigen_solver(N,Er,El,lamr) result(info)
    use kinds, only: r8 => real_kind
    use parallel_mod, only : abortmp
    integer,  intent(in)     :: N          ! system size
    real(r8), intent(inout)    :: Er(N,N), El(N,N)
    real(r8), intent(out)    :: lamr(N)
    integer info

    ! Local variables
    integer n2
    real(r8) :: work(4*N)
    integer :: lwork

    real(r8) :: wr(N), wi(N), vr(N,N)
#ifdef CAM
    call abortmp('not supported in cam at this time')
#else
    lwork = 4*N
    n2 = n
!    print *,'calling dgeev'
!    call DSYEV( 'V', 'U', N, Er, N2, lamr, work, lwork, info )
    call DGEEV( 'V', 'V', N, Er, N, lamr, wi, El, N, vr, N, work, lwork, info )
!    print *,'done calling dgeev'
    Er = vr
#endif
  end function reigen_solver

end module linear_algebra_mod
