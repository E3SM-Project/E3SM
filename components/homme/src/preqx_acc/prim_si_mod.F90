
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module prim_si_mod
  use prim_si_mod_base, only: preq_vertadv, preq_vertadv_v, preq_vertadv_upwind, preq_omega_ps, preq_hydrostatic, preq_hydrostatic_v2, &
  geopotential_t, prim_set_mass
#if USE_OPENACC
  use openacc, only: acc_async_sync
  implicit none


contains



  subroutine preq_hydrostatic_openacc(phi,phis,T_v,p,dp,nets,nete,ntl,tl)
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev
    use physical_constants, only : rgas
    implicit none
    real(kind=real_kind), intent(  out) :: phi (np,np,nlev,nets:nete)
    real(kind=real_kind), intent(in   ) :: phis(np,np,nets:nete)
    real(kind=real_kind), intent(in   ) :: T_v (np,np,nlev,nets:nete)
    real(kind=real_kind), intent(in   ) :: p   (np,np,nlev,nets:nete)
    real(kind=real_kind), intent(in   ) :: dp  (np,np,nlev,ntl,nets:nete)
    integer             , intent(in   ) :: nets,nete,ntl,tl
    integer i,j,k,ie                      ! longitude, level indices
    real(kind=real_kind) Hkk,Hkl          ! diagonal term of energy conversion matrix
    real(kind=real_kind), dimension(nlev) :: phii       ! Geopotential at interfaces

    !$acc parallel loop gang vector collapse(3) private(hkk,hkl,phii)
    do ie = nets , nete
      do j=1,np
        do i=1,np
          hkk = dp(i,j,nlev,tl,ie)*0.5d0/p(i,j,nlev,ie)
          hkl = 2*hkk
          phii(nlev)  = Rgas*T_v(i,j,nlev,ie)*hkl
          phi (i,j,nlev,ie) = phis(i,j,ie) + Rgas*T_v(i,j,nlev,ie)*hkk
          do k=nlev-1,2,-1
            hkk = dp(i,j,k,tl,ie)*0.5d0/p(i,j,k,ie)
            hkl = 2*hkk
            phii(k) = phii(k+1) + Rgas*T_v(i,j,k,ie)*hkl
            phi(i,j,k,ie) = phis(i,j,ie) + phii(k+1) + Rgas*T_v(i,j,k,ie)*hkk
          enddo
          hkk = 0.5d0*dp(i,j,1,tl,ie)/p(i,j,1,ie)
          phi(i,j,1,ie) = phis(i,j,ie) + phii(2) + Rgas*T_v(i,j,1,ie)*hkk
        enddo
      enddo
    enddo
  end subroutine preq_hydrostatic_openacc



  subroutine preq_omega_ps_openacc(omega_p,hvcoord,p,vgrad_p,divdp,nets,nete)
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev
    use hybvcoord_mod, only : hvcoord_t
    implicit none
    real(kind=real_kind), intent(in   ) :: divdp  (np,np,nlev,nets:nete)      ! divergence
    real(kind=real_kind), intent(in   ) :: vgrad_p(np,np,nlev,nets:nete) ! v.grad(p)
    real(kind=real_kind), intent(in   ) :: p      (np,np,nlev,nets:nete)     ! layer thicknesses (pressure)
    type (hvcoord_t),     intent(in   ) :: hvcoord
    real(kind=real_kind), intent(  out) :: omega_p(np,np,nlev,nets:nete)   ! vertical pressure velocity
    integer             , intent(in   ) :: nets , nete
    integer i,j,k,ie                      ! longitude, level indices
    real(kind=real_kind) term             ! one half of basic term in omega/p summation
    real(kind=real_kind) Ckk,Ckl          ! diagonal term of energy conversion matrix
    real(kind=real_kind) suml             ! partial sum over l = (1, k-1)
    !$acc parallel loop gang vector collapse(3) private(ckk,term,suml,ckl)
    do ie = nets,nete
      do j=1,np   !   Loop inversion (AAM)
        do i=1,np
          ckk = 0.5d0/p(i,j,1,ie)
          term = divdp(i,j,1,ie)
          omega_p(i,j,1,ie) = vgrad_p(i,j,1,ie)/p(i,j,1,ie) - ckk*term
          suml = term
          do k=2,nlev-1
            ckk = 0.5d0/p(i,j,k,ie)
            ckl = 2*ckk
            term = divdp(i,j,k,ie)
            omega_p(i,j,k,ie) = vgrad_p(i,j,k,ie)/p(i,j,k,ie) - ckl*suml - ckk*term
            suml = suml + term
          enddo
          ckk = 0.5d0/p(i,j,nlev,ie)
          ckl = 2*ckk
          term = divdp(i,j,nlev,ie)
          omega_p(i,j,nlev,ie) = vgrad_p(i,j,nlev,ie)/p(i,j,nlev,ie) - ckl*suml - ckk*term
        enddo
      enddo
    enddo
  end subroutine preq_omega_ps_openacc



  subroutine preq_vertadv_openacc(T, v, eta_dot_dp_deta, rpdel, T_vadv, v_vadv)
    use kinds,              only : real_kind
    use dimensions_mod,     only : nlev, np, nlevp
    implicit none
    real (kind=real_kind), intent(in) :: T(np,np,nlev)
    real (kind=real_kind), intent(in) :: v(np,np,2,nlev)
    real (kind=real_kind), intent(in) :: eta_dot_dp_deta(np,np,nlevp)
    real (kind=real_kind), intent(in) :: rpdel(np,np,nlev)
    real (kind=real_kind), intent(out) :: T_vadv(np,np,nlev)
    real (kind=real_kind), intent(out) :: v_vadv(np,np,2,nlev)
    integer :: i,j,k
    real (kind=real_kind) :: facp, facm
    do j=1,np   !   Loop inversion (AAM)
      ! Compute vertical advection of T and v from eq. (3.b.1)
      ! k = 1 case:
      k=1
      do i=1,np
        facp            = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k+1)
        T_vadv(i,j,k)   = facp*(T(i,j,k+1)- T(i,j,k))
        v_vadv(i,j,1,k) = facp*(v(i,j,1,k+1)- v(i,j,1,k))
        v_vadv(i,j,2,k) = facp*(v(i,j,2,k+1)- v(i,j,2,k))
      enddo
      ! vertical advection
      ! 1 < k < nlev case:
      do k=2,nlev-1
        do i=1,np
          facp            = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k+1)
          facm            = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k)
          T_vadv(i,j,k)   = facp*(T(i,j,k+1)- T(i,j,k)) + facm*(T(i,j,k)- T(i,j,k-1))
          v_vadv(i,j,1,k) = facp*(v(i,j,1,k+1)- v(i,j,1,k)) + facm*(v(i,j,1,k)- v(i,j,1,k-1))
          v_vadv(i,j,2,k) = facp*(v(i,j,2,k+1)- v(i,j,2,k)) + facm*(v(i,j,2,k)- v(i,j,2,k-1))
        enddo
      enddo
      ! vertical advection
      ! k = nlev case:
      k=nlev
      do i=1,np
        facm            = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k)
        T_vadv(i,j,k)   = facm*(T(i,j,k)- T(i,j,k-1))
        v_vadv(i,j,1,k) = facm*(v(i,j,1,k)- v(i,j,1,k-1))
        v_vadv(i,j,2,k) = facm*(v(i,j,2,k)- v(i,j,2,k-1))
      enddo
    enddo
  end subroutine preq_vertadv_openacc



#endif
end module prim_si_mod
