
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



  subroutine preq_hydrostatic_openacc(phi,phis,T_v,p,dp)
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev
    use physical_constants, only : rgas
    implicit none
    real(kind=real_kind), intent(  out) :: phi (np,np,nlev)
    real(kind=real_kind), intent(in   ) :: phis(np,np)
    real(kind=real_kind), intent(in   ) :: T_v (np,np,nlev)
    real(kind=real_kind), intent(in   ) :: p   (np,np,nlev)
    real(kind=real_kind), intent(in   ) :: dp  (np,np,nlev)
    integer i,j,k                         ! longitude, level indices
    real(kind=real_kind) Hkk,Hkl          ! diagonal term of energy conversion matrix
    real(kind=real_kind), dimension(nlev) :: phii       ! Geopotential at interfaces

      do j=1,np
        do i=1,np
          hkk = dp(i,j,nlev)*0.5d0/p(i,j,nlev)
          hkl = 2*hkk
          phii(nlev)  = Rgas*T_v(i,j,nlev)*hkl
          phi (i,j,nlev) = phis(i,j) + Rgas*T_v(i,j,nlev)*hkk
          do k=nlev-1,2,-1
            hkk = dp(i,j,k)*0.5d0/p(i,j,k)
            hkl = 2*hkk
            phii(k) = phii(k+1) + Rgas*T_v(i,j,k)*hkl
            phi(i,j,k) = phis(i,j) + phii(k+1) + Rgas*T_v(i,j,k)*hkk
          end do
          hkk = 0.5d0*dp(i,j,1)/p(i,j,1)
          phi(i,j,1) = phis(i,j) + phii(2) + Rgas*T_v(i,j,1)*hkk
        end do
      end do
  end subroutine preq_hydrostatic_openacc



#endif
end module prim_si_mod
