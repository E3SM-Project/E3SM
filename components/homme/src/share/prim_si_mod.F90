#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_si_mod
  implicit none
  private
  public :: preq_omegap
  public :: preq_omega_ps
  public :: preq_omega_lnps
  public :: preq_hydrostatic, geopotential_t
  public :: preq_pressure
  public :: preq_vertadv
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




!----------------------------------------------------------------------- 
! preq_omegap:

! Purpose: 
! Calculate (omega/p) needed for the Thermodynamics Equation
! 
! Method: 
! Simplified version in CAM2 for clarity
! 
! Author: Modified by Rich Loft for use in HOMME. 
! 
!-----------------------------------------------------------------------

  subroutine preq_omegap(div     ,vgrad_ps,pdel    ,rpmid, &
       hybm    ,hybd    ,omegap   )
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev
    implicit none


    !------------------------------Arguments---------------------------------------------------------------
    real(kind=real_kind), intent(in) :: div(np,np,nlev)      ! divergence
    real(kind=real_kind), intent(in) :: vgrad_ps(np,np,nlev) ! v.grad(ps)
    real(kind=real_kind), intent(in) :: pdel(np,np,nlev)     ! layer thicknesses (pressure)
    real(kind=real_kind), intent(in) :: rpmid(np,np,nlev)    ! 1./pmid
    real(kind=real_kind), intent(in) :: hybm(nlev)           ! Hybrid B coefficient on mid levels
    real(kind=real_kind), intent(in) :: hybd(nlev)           ! Hybrid dB coefficient on mid levels
    real(kind=real_kind), intent(out):: omegap(np,np,nlev)   ! vertical pressure velocity
    !------------------------------------------------------------------------------------------------------

    !---------------------------Local workspace-----------------------------
    integer i,j,k                         ! longitude, level indices
    real(kind=real_kind) term             ! one half of basic term in omega/p summation 
    real(kind=real_kind) Ckk              ! diagonal term of energy conversion matrix
    real(kind=real_kind) suml(np,np)      ! partial sum over l = (1, k-1)
    !-----------------------------------------------------------------------

    ! =========================
    ! Zero partial sum
    ! =========================

    do j=1,np
       do i=1,np
          suml(i,j)=0.0_real_kind
       end do
    end do

    ! =============================
    ! Compute omegap 
    ! =============================

    do k=1,nlev
       do j=1,np
          do i=1,np
             Ckk       = 0.5_real_kind
             term      = Ckk*(div(i,j,k)*pdel(i,j,k) + vgrad_ps(i,j,k)*hybd(k))
             suml(i,j) = suml(i,j) + term
             omegap(i,j,k) = rpmid(i,j,k)*(hybm(k)*vgrad_ps(i,j,k) - suml(i,j))
             suml(i,j) = suml(i,j) + term
          end do
       end do
    end do

  end subroutine preq_omegap



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
!  compute omega/p using lnps 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine preq_omega_lnps(omega_p,hvcoord,ps,p,dp,vgrad_lnps,div)
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev
    use hybvcoord_mod, only : hvcoord_t
    implicit none


    !------------------------------Arguments---------------------------------------------------------------
    real(kind=real_kind), intent(in) :: div(np,np,nlev)      ! divergence
    real(kind=real_kind), intent(in) :: vgrad_lnps(np,np,nlev) ! v.grad(ps)
    real(kind=real_kind), intent(in) :: p(np,np,nlev)      ! pressure
    real(kind=real_kind), intent(in) :: dp(np,np,nlev)     ! dp/dn
    real(kind=real_kind), intent(in) :: ps(np,np)
    type (hvcoord_t),     intent(in) :: hvcoord
    real(kind=real_kind), intent(out):: omega_p(np,np,nlev)   ! vertical pressure velocity
    !------------------------------------------------------------------------------------------------------

    !---------------------------Local workspace-----------------------------
    integer i,j,k                         ! longitude, level indices
    real(kind=real_kind) term             ! one half of basic term in omega/p summation 
    real(kind=real_kind) Ckk,Ckl          ! diagonal term of energy conversion matrix
    real(kind=real_kind) suml(np,np)      ! partial sum over l = (1, k-1)
    !-----------------------------------------------------------------------
       do j=1,np	
          do i=1,np
             ckk = 0.5d0/p(i,j,1)
             term = div(i,j,1)*dp(i,j,1) + vgrad_lnps(i,j,1)*ps(i,j)*hvcoord%hybd(1)
             omega_p(i,j,1) = hvcoord%hybm(1)*(ps(i,j)/p(i,j,1))*vgrad_lnps(i,j,1)
             omega_p(i,j,1) = omega_p(i,j,1) - ckk*term
             suml(i,j) = term
          end do
       end do

       do k=2,nlev-1

          do j=1,np
             do i=1,np
                ckk = 0.5d0/p(i,j,k)
                ckl = 2*ckk
                term = div(i,j,k)*dp(i,j,k) + vgrad_lnps(i,j,k)*ps(i,j)*hvcoord%hybd(k)
                omega_p(i,j,k) = hvcoord%hybm(k)*(ps(i,j)/p(i,j,k))*vgrad_lnps(i,j,k)
                omega_p(i,j,k) = omega_p(i,j,k) - ckl*suml(i,j) - ckk*term
                suml(i,j) = suml(i,j) + term

             end do
          end do

       end do

       do j=1,np
          do i=1,np
             ckk = 0.5d0/p(i,j,nlev)
             ckl = 2*ckk
             term = div(i,j,nlev)*dp(i,j,nlev) + vgrad_lnps(i,j,nlev)*ps(i,j)*hvcoord%hybd(nlev)
             omega_p(i,j,nlev) = hvcoord%hybm(nlev)*(ps(i,j)/p(i,j,nlev))*vgrad_lnps(i,j,nlev)
             omega_p(i,j,nlev) = omega_p(i,j,nlev) - ckl*suml(i,j) - ckk*term
          end do
       end do

  end subroutine preq_omega_lnps



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





!----------------------------------------------------------------------- 
! preq_pressure:
!
! Purpose: 
! Define the pressures of the interfaces and midpoints from the
! coordinate definitions and the surface pressure. Originally plevs0!
! 
! Method: 
! 
! Author: B. Boville/ Adapted for HOMME by Rich Loft
! 
!-----------------------------------------------------------------------
!
! $Id: prim_si_mod.F90,v 2.10 2005/10/14 20:17:22 jedwards Exp $
! $Author: jedwards $
!
!-----------------------------------------------------------------------

  subroutine preq_pressure (ps0,  ps,               &
       hyai, hybi, hyam, hybm, &
       pint, pmid, pdel)
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nlevp
    implicit none

    !-----------------------------------------------------------------------

    real(kind=real_kind), intent(in)  :: ps0                ! Hybrid coordinate reference pressure (pascals)
    real(kind=real_kind), intent(in)  :: ps(np,np)          ! Surface pressure (pascals)
    real(kind=real_kind), intent(in)  :: hyai(nlevp)        ! Hybrid interface A coefficients
    real(kind=real_kind), intent(in)  :: hybi(nlevp)        ! Hybrid interface B coefficients
    real(kind=real_kind), intent(in)  :: hyam(nlev)         ! Hybrid midpoint  A coefficients
    real(kind=real_kind), intent(in)  :: hybm(nlev)         ! Hybrid midpoint  B coefficients
    real(kind=real_kind), intent(out) :: pint(np,np,nlevp)  ! Pressure at model interfaces
    real(kind=real_kind), intent(out) :: pmid(np,np,nlev)   ! Pressure at model levels
    real(kind=real_kind), intent(out) :: pdel(np,np,nlev)   ! Layer thickness (pint(k+1) - pint(k))
    !-----------------------------------------------------------------------

    !---------------------------Local workspace-----------------------------
    integer i,j,k             ! Horizontal, level indices
    !-----------------------------------------------------------------------
    !
    ! Set interface pressures
    !
    do k=1,nlevp
       do j=1,np
          do i=1,np
             pint(i,j,k) = hyai(k)*ps0 + hybi(k)*ps(i,j)
          end do
       end do
    end do
    !
    ! Set midpoint pressures and layer thicknesses
    !
    do k=1,nlev
       do j=1,np
          do i=1,np
             pmid(i,j,k) = hyam(k)*ps0 + hybm(k)*ps(i,j)
             pdel(i,j,k) = pint(i,j,k+1) - pint(i,j,k)
          end do
       end do
    end do

  end subroutine preq_pressure




end module prim_si_mod
