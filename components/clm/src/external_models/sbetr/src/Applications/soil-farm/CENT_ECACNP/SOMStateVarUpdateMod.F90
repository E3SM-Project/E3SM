module SOMStateVarUpdateMod
  !
  ! DESCRIPTION:
  ! subroutines to update state variables of any
  ! reaction based bgc module
  ! !USES:
  use bshr_kind_mod       , only: r8 => shr_kind_r8
  implicit none

  private
  character(len=*), private, parameter :: filename = &
       __FILE__

  public :: calc_dtrend_som_bgc

contains

  !-----------------------------------------------------------------------
  subroutine calc_dtrend_som_bgc(nx, ny, cascade_matrix, reaction_rates, dxdt)
    !
    ! !DESCRIPTION:
    ! return the temporal trend of the state variables
    implicit none
    ! !ARGUMENTS:
    integer,  intent(in) :: nx, ny
    real(r8), intent(in) :: cascade_matrix(1:nx,1:ny)
    real(r8), intent(in) :: reaction_rates(1:ny)
    real(r8), intent(out):: dxdt(1:nx)

    !intel mkl f90 interface
    !call gemv(cascade_matrix, reaction_rates, dxdt, alpha=1._r8, beta=0._r8)
    ! BLAS INTERFACE
    ! y:= alpha * A * x + beta * y
    ! DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    call dgemv('N', nx, ny, 1._r8, cascade_matrix, nx, reaction_rates, 1, 0._r8, dxdt, 1)

  end subroutine calc_dtrend_som_bgc

end module SOMStateVarUpdateMod
