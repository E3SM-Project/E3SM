module limiters_interface_mod
  use dimensions_mod,       only : np, nlev
  use kinds,                only : real_kind
  use derivative_mod_base,  only : limiter_optim_iter_full, limiter_clip_and_sum

  implicit none
  private

  public  :: limiter_optim_iter_full_c_callable, limiter_clip_and_sum_c_callable

contains

  subroutine limiter_optim_iter_full_c_callable(ptens,sphweights,minp,maxp,dpmass) bind(c)
    real (kind=real_kind), dimension(np,np,nlev), intent(inout) :: ptens
    real (kind=real_kind), dimension(np,np),      intent(in)    :: sphweights
    real (kind=real_kind), dimension(nlev),       intent(inout) :: minp, maxp
    real (kind=real_kind), dimension(np,np,nlev), intent(in)    :: dpmass
    call limiter_optim_iter_full(ptens,sphweights,minp,maxp,dpmass)
  end subroutine limiter_optim_iter_full_c_callable

  subroutine limiter_clip_and_sum_c_callable(ptens,sphweights,minp,maxp,dpmass) bind(c)
    real (kind=real_kind), dimension(np,np,nlev), intent(inout) :: ptens
    real (kind=real_kind), dimension(np,np),      intent(in)    :: sphweights
    real (kind=real_kind), dimension(nlev),       intent(inout) :: minp, maxp
    real (kind=real_kind), dimension(np,np,nlev), intent(in)    :: dpmass
    call limiter_clip_and_sum(ptens,sphweights,minp,maxp,dpmass)
  end subroutine limiter_clip_and_sum_c_callable

end module limiters_interface_mod
