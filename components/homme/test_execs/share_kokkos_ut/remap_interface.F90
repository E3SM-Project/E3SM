module remap_interface

  use iso_c_binding,  only: c_int
  use dimensions_mod, only: nlev
  use kinds,          only: real_kind 
  use parallel_mod,   only: abortmp
  use control_mod,    only: vert_remap_q_alg
  implicit none

contains

  subroutine remap_Q_ppm_c_callable(Qdp,nx,qsize,dp1,dp2,alg) bind(c)
    use vertremap_base, only: remap_Q_ppm
    !
    ! Inputs
    !
    integer(c_int),        intent(in)    :: nx,qsize,alg
    real (kind=real_kind), intent(inout) :: Qdp(nx,nx,nlev,qsize)
    real (kind=real_kind), intent(in)    :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)

    !aim for alg=1 or alg=2 only
    if(alg == 1 .or. alg == 2 .or. alg == 3) then
      ! Need to set alg in control_mod, b/c fortran reads it from there
      vert_remap_q_alg = alg
      call remap_Q_ppm(Qdp,nx,qsize,dp1,dp2)
    else
      call abortmp('compute_ppm_grids_c_callable: bad alg (not 1 or 2) .')
    endif

  end subroutine remap_Q_ppm_c_callable

  subroutine compute_ppm_grids_c_callable(dx,rslt,alg) bind(c)
    use vertremap_base, only: compute_ppm_grids 
    !
    ! Inputs
    !
    integer(c_int),       intent(in)  :: alg
    real(kind=real_kind), intent(in)  :: dx(-1:nlev+2)
    real(kind=real_kind), intent(out) :: rslt(10,0:nlev+1)  !grid spacings

    !aim for alg=1 or alg=2 only
    if(alg == 1 .or. alg == 2 .or. alg == 3) then
      ! Need to set alg in control_mod, b/c fortran reads it from there
      vert_remap_q_alg = alg
      rslt = compute_ppm_grids(dx)
      rslt(4:10, nlev + 1) = 0.0
    else
      call abortmp('compute_ppm_grids_c_callable: bad alg (not 1, 2, or 3) .')
    endif
  end subroutine compute_ppm_grids_c_callable

  subroutine compute_ppm_c_callable(a,dx,coefs,alg) bind(c)
    use vertremap_base, only: compute_ppm
    !
    ! Inputs
    !
    integer(c_int), intent(in) :: alg
    real(kind=real_kind), intent(in) :: a    (    -1:nlev+2)  !Cell-mean values
    real(kind=real_kind), intent(in) :: dx   (10,  0:nlev+1)  !grid spacings
    real(kind=real_kind), intent(out):: coefs(0:2,   nlev  )  !PPM coefficients (for parabola)

    !aim for alg=1 or alg=2 only
    if(alg == 1 .or. alg == 2 .or. alg == 3) then
      ! Need to set alg in control_mod, b/c fortran reads it from there
      vert_remap_q_alg = alg
      coefs = compute_ppm(a,dx)
    else
      call abortmp('compute_ppm_c_callable: bad alg (not 1, 2, or 3) .')
    endif
  end subroutine compute_ppm_c_callable

end module remap_interface
