#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_si_mod
  implicit none
contains
	
! ==========================================================
! Implicit system for semi-implicit primitive equations.
! ==========================================================


  subroutine preq_vertadv(T, v, eta_dot_dp_deta, rpdel, &
       T_vadv, v_vadv)
    use kinds,              only : real_kind
    use dimensions_mod,     only : nlev, np, nlevp
    implicit none
    
    real (kind=real_kind), intent(in) :: T(np,np,nlev)
    real (kind=real_kind), intent(in) :: v(np,np,2,nlev)
    real (kind=real_kind), intent(in) :: eta_dot_dp_deta(np,np,nlevp)
    real (kind=real_kind), intent(in) :: rpdel(np,np,nlev)

    real (kind=real_kind), intent(out) :: T_vadv(np,np,nlev)
    real (kind=real_kind), intent(out) :: v_vadv(np,np,2,nlev)

    ! ========================
    ! Local Variables
    ! ========================

    integer :: i,j,k
    real (kind=real_kind) :: facp, facm

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,j,i,facp,facm)
#endif
    do j=1,np   !   Loop inversion (AAM)

    ! ===========================================================
    ! Compute vertical advection of T and v from eq. (3.b.1)
    !
    ! k = 1 case:
    ! ===========================================================

       k=1
       do i=1,np 
          facp            = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k+1)
          T_vadv(i,j,k)   = facp*(T(i,j,k+1)- T(i,j,k))
          v_vadv(i,j,1,k) = facp*(v(i,j,1,k+1)- v(i,j,1,k))
          v_vadv(i,j,2,k) = facp*(v(i,j,2,k+1)- v(i,j,2,k))
       end do

    ! ===========================================================
    ! vertical advection
    !
    ! 1 < k < nlev case:
    ! ===========================================================

       do k=2,nlev-1
          do i=1,np
             facp            = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k+1)
             facm            = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k)
             T_vadv(i,j,k)   = facp*(T(i,j,k+1)- T(i,j,k)) + &
                  facm*(T(i,j,k)- T(i,j,k-1))
             v_vadv(i,j,1,k) = facp*(v(i,j,1,k+1)- v(i,j,1,k)) + &
                  facm*(v(i,j,1,k)- v(i,j,1,k-1))
             v_vadv(i,j,2,k) = facp*(v(i,j,2,k+1)- v(i,j,2,k)) + &
                  facm*(v(i,j,2,k)- v(i,j,2,k-1))
          end do
       end do

    ! ===========================================================
    ! vertical advection
    !
    ! k = nlev case:
    ! ===========================================================

       k=nlev
       do i=1,np
          facm            = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k)
          T_vadv(i,j,k)   = facm*(T(i,j,k)- T(i,j,k-1))
          v_vadv(i,j,1,k) = facm*(v(i,j,1,k)- v(i,j,1,k-1))
          v_vadv(i,j,2,k) = facm*(v(i,j,2,k)- v(i,j,2,k-1))
       end do

    end do

  end subroutine preq_vertadv







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!  compute omega/p using ps, modeled after CCM3 formulas 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine preq_omega_ps(omega_p,hvcoord,p,vgrad_p,divdp)
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev
    use hybvcoord_mod, only : hvcoord_t
    implicit none


    !------------------------------Arguments---------------------------------------------------------------
    real(kind=real_kind), intent(in) :: divdp(np,np,nlev)      ! divergence
    real(kind=real_kind), intent(in) :: vgrad_p(np,np,nlev) ! v.grad(p)
    real(kind=real_kind), intent(in) :: p(np,np,nlev)     ! layer thicknesses (pressure)
    type (hvcoord_t),     intent(in) :: hvcoord
    real(kind=real_kind), intent(out):: omega_p(np,np,nlev)   ! vertical pressure velocity
    !------------------------------------------------------------------------------------------------------

    !---------------------------Local workspace-----------------------------
    integer i,j,k                         ! longitude, level indices
    real(kind=real_kind) term             ! one half of basic term in omega/p summation 
    real(kind=real_kind) Ckk,Ckl          ! diagonal term of energy conversion matrix
    real(kind=real_kind) suml(np,np)      ! partial sum over l = (1, k-1)
    !-----------------------------------------------------------------------

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,j,i,ckk,term,ckl)
#endif
       do j=1,np   !   Loop inversion (AAM)

          do i=1,np
             ckk = 0.5d0/p(i,j,1)
             term = divdp(i,j,1)
!             omega_p(i,j,1) = hvcoord%hybm(1)*vgrad_ps(i,j,1)/p(i,j,1)
             omega_p(i,j,1) = vgrad_p(i,j,1)/p(i,j,1)
             omega_p(i,j,1) = omega_p(i,j,1) - ckk*term
             suml(i,j) = term
          end do

          do k=2,nlev-1
             do i=1,np
                ckk = 0.5d0/p(i,j,k)
                ckl = 2*ckk
                term = divdp(i,j,k)
!                omega_p(i,j,k) = hvcoord%hybm(k)*vgrad_ps(i,j,k)/p(i,j,k)
                omega_p(i,j,k) = vgrad_p(i,j,k)/p(i,j,k)
                omega_p(i,j,k) = omega_p(i,j,k) - ckl*suml(i,j) - ckk*term
                suml(i,j) = suml(i,j) + term

             end do
          end do

          do i=1,np
             ckk = 0.5d0/p(i,j,nlev)
             ckl = 2*ckk
             term = divdp(i,j,nlev)
!             omega_p(i,j,nlev) = hvcoord%hybm(nlev)*vgrad_ps(i,j,nlev)/p(i,j,nlev)
             omega_p(i,j,nlev) = vgrad_p(i,j,nlev)/p(i,j,nlev)
             omega_p(i,j,nlev) = omega_p(i,j,nlev) - ckl*suml(i,j) - ckk*term
          end do

       end do

  end subroutine preq_omega_ps







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!  CCM3 hydrostatic integral
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine preq_hydrostatic(phi,phis,T_v,p,dp)
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev
    use physical_constants, only : rgas
!    use hybvcoord_mod, only : hvcoord_t
    implicit none


    !------------------------------Arguments---------------------------------------------------------------
    real(kind=real_kind), intent(out) :: phi(np,np,nlev)     
    real(kind=real_kind), intent(in) :: phis(np,np)
    real(kind=real_kind), intent(in) :: T_v(np,np,nlev)
    real(kind=real_kind), intent(in) :: p(np,np,nlev)   
    real(kind=real_kind), intent(in) :: dp(np,np,nlev)  
 !   type (hvcoord_t),     intent(in) :: hvcoord
    !------------------------------------------------------------------------------------------------------

    !---------------------------Local workspace-----------------------------
    integer i,j,k                         ! longitude, level indices
    real(kind=real_kind) Hkk,Hkl          ! diagonal term of energy conversion matrix
    real(kind=real_kind), dimension(np,np,nlev) :: phii       ! Geopotential at interfaces
    !-----------------------------------------------------------------------

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,j,i,hkk,hkl)
#endif
       do j=1,np   !   Loop inversion (AAM)

          do i=1,np
             hkk = dp(i,j,nlev)*0.5d0/p(i,j,nlev)
             hkl = 2*hkk
             phii(i,j,nlev)  = Rgas*T_v(i,j,nlev)*hkl
             phi(i,j,nlev) = phis(i,j) + Rgas*T_v(i,j,nlev)*hkk 
          end do

          do k=nlev-1,2,-1
             do i=1,np
                ! hkk = dp*ckk
                hkk = dp(i,j,k)*0.5d0/p(i,j,k)
                hkl = 2*hkk
                phii(i,j,k) = phii(i,j,k+1) + Rgas*T_v(i,j,k)*hkl
                phi(i,j,k) = phis(i,j) + phii(i,j,k+1) + Rgas*T_v(i,j,k)*hkk
             end do
          end do

          do i=1,np
             ! hkk = dp*ckk
             hkk = 0.5d0*dp(i,j,1)/p(i,j,1)
             phi(i,j,1) = phis(i,j) + phii(i,j,2) + Rgas*T_v(i,j,1)*hkk
          end do

       end do


end subroutine preq_hydrostatic



!
!  The hydrostatic routine from CAM physics.
!  (FV stuff removed)
!  t,q input changed to take t_v
!  removed gravit, so this routine returns PHI, not zm
subroutine geopotential_t(                                 &
       pmid   , pdel   ,  tv      , rair   ,  zm)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the geopotential height (above the surface) at the midpoints and 
! interfaces using the input temperatures and pressures.
!
!-----------------------------------------------------------------------
    use dimensions_mod,     only : nlev, nlevp, np
    use kinds, only : real_kind
    implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments



    real(real_kind), intent(in) :: pmid (np*np,nlev)    ! Midpoint pressures
    real(real_kind), intent(in) :: pdel (np*np,nlev)    ! layer thickness
    real(real_kind), intent(in) :: tv    (np*np,nlev)    ! temperature
    real(real_kind), intent(in) :: rair                 ! Gas constant for dry air
    ! real(real_kind), intent(in) :: gravit               ! Acceleration of gravity
    ! real(real_kind), intent(in) :: zvir                 ! rh2o/rair - 1

! Output arguments

    real(real_kind), intent(out) :: zm(np*np,nlev)      ! Geopotential height at mid level
!
!---------------------------Local variables-----------------------------
    integer :: ncol=np*np             ! Number of longitudes

    integer  :: i,k                ! Lon, level indices
    real(real_kind) :: hkk(np*np)         ! diagonal element of hydrostatic matrix
    real(real_kind) :: hkl(np*np)         ! off-diagonal element
    real(real_kind) :: rog                ! Rair / gravit
    real(real_kind) :: zi(np*np,nlevp)     ! Height above surface at interfaces
!
!-----------------------------------------------------------------------
!
!    rog = rair/gravit
    rog = rair

! The surface height is zero by definition.
    do i = 1,ncol
       zi(i,nlevp) = 0.0_real_kind
    end do

! Compute zi, zm from bottom up. 
! Note, zi(i,k) is the interface above zm(i,k)
    do k = nlev, 1, -1
! First set hydrostatic elements consistent with dynamics
       do i = 1,ncol
          hkl(i) = pdel(i,k) / pmid(i,k)
          hkk(i) = 0.5_real_kind * hkl(i)
       end do

! Now compute tv, zm, zi
       do i = 1,ncol
          ! tvfac   = 1._r8 + zvir * q(i,k)
          ! tv      = t(i,k) * tvfac
          zm(i,k) = zi(i,k+1) + rog * tv(i,k) * hkk(i)
          zi(i,k) = zi(i,k+1) + rog * tv(i,k) * hkl(i)
       end do
    end do

    return
  end subroutine geopotential_t



  subroutine prim_set_mass(elem, tl,hybrid,hvcoord,nets,nete)
  use kinds, only : real_kind
  use control_mod, only : initial_total_mass
  use physical_constants, only : g
  use element_mod, only : element_t
  use time_mod, only : timelevel_t 
  use hybvcoord_mod, only : hvcoord_t 
  use hybrid_mod, only : hybrid_t
  use dimensions_mod, only : np
  use global_norms_mod, only : global_integral 

  type (element_t), intent(inout) :: elem(:)
  type (TimeLevel_t), target, intent(in) :: tl
  type (hybrid_t),intent(in)     :: hybrid
  type (hvcoord_t), intent(in)   :: hvcoord
  integer,intent(in)             :: nets,nete
  
  ! local 
  real (kind=real_kind)  :: tmp(np,np,nets:nete)
  real (kind=real_kind)  :: scale,mass0
  integer :: n0,nm1,np1,ie

  if (initial_total_mass == 0) return;
  
  n0=tl%n0
  nm1=tl%nm1
  np1=tl%np1
  
  scale=1/g                                  ! assume code is using Pa
  if (hvcoord%ps0 <  2000 ) scale=100*scale  ! code is using mb
  ! after scaling, Energy is in J/m**2,  Mass kg/m**2
  
  do ie=nets,nete
     tmp(:,:,ie)=elem(ie)%state%ps_v(:,:,n0)
  enddo
  mass0 = global_integral(elem, tmp(:,:,nets:nete),hybrid,np,nets,nete)
  mass0 = mass0*scale;  
  
  do ie=nets,nete
     elem(ie)%state%ps_v(:,:,n0)=elem(ie)%state%ps_v(:,:,n0)*(initial_total_mass/mass0)
     elem(ie)%state%ps_v(:,:,np1)=elem(ie)%state%ps_v(:,:,n0)
     elem(ie)%state%ps_v(:,:,nm1)=elem(ie)%state%ps_v(:,:,n0)
  enddo
  if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
     write (*,'(a,e24.15)') "Initializing Total Mass (kg/m^2) = ",initial_total_mass
  endif
  end subroutine prim_set_mass





end module prim_si_mod
