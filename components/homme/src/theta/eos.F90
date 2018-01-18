#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
!
!  NonHydrostatic Equation of State and inverse EOS routines
!  Note: these routines should all be discrete inverses of each other
!
!  get_pnh_and_exner()       Compute nonhydrostatic pressure as a function of 
!                            potential temperature and geopotential
!
!  get_dry_phinh()           Compute geopotential as a function of potential temperature
!                            and pressure (neglicting water vapor)
!
!  get_wet_phinh()           Compute geopotential as a function of potential temperature
!                            and pressure, including water vapor
!  
!  Original version: Mark Taylor 2017/1
!
module eos

  use dimensions_mod, only: np, nlev, nlevp, nelemd
  use element_mod,    only: element_t
  use element_state,  only: timelevels, elem_state_t
  use hybvcoord_mod,  only: hvcoord_t
  use kinds,          only: real_kind
  use parallel_mod,   only: abortmp
  use physical_constants, only : p0, kappa, g
  use control_mod,    only: theta_hydrostatic_mode
  use prim_si_mod,    only: preq_hydrostatic_v2, preq_omega_ps
  implicit none


contains


  subroutine get_pnh_and_exner(hvcoord,theta_dp_cp,dp3d,phi,phis,kappa_star,pnh,dpnh,exner,exner_i_out,pnh_i_out)
  implicit none
!
! Use Equation of State to compute exner pressure, nh presure
! hydrostatic EOS:
!          compute p, exner  (only input needed id dp3d and kappa_star)
!
! nonhydrostatic EOS:   
!
!       theta_dp_cp * ep = -dphi/ds    OR    rho_R_theta * dphi/ds = -theta_dp_cp*kappa
!
! with  ep = d(exner)/dp = kappa*exner/p  ( used in tech report)
! rho_R_theta = p/exner  (used here)
!
! input:  dp3d, Qdp (if use_moisture), phi, phis, theta
! output:  pnh, dphn, exner, exner_i, pnh_i
!
! NOTE: Exner pressure is defined in terms of p0=1000mb.  Be sure to use global constant p0,
! instead of hvcoord%ps0, which is set by CAM to ~1021mb
!  
  type (hvcoord_t),     intent(in)  :: hvcoord             ! hybrid vertical coordinate struct
  real (kind=real_kind), intent(in) :: theta_dp_cp(np,np,nlev)   
  real (kind=real_kind), intent(in) :: dp3d(np,np,nlev)   
  real (kind=real_kind), intent(in) :: phi(np,np,nlev)
  real (kind=real_kind), intent(in) :: phis(np,np)
  real (kind=real_kind), intent(in) :: kappa_star(np,np,nlev)   
  real (kind=real_kind), intent(out) :: pnh(np,np,nlev)   ! nh nonhyrdo pressure
  real (kind=real_kind), intent(out) :: dpnh(np,np,nlev) ! nh nonhyrdo pressure interfaces
  real (kind=real_kind), intent(out) :: exner(np,np,nlev)  ! exner nh pressure
  real (kind=real_kind), intent(out), optional :: exner_i_out(np,np,nlevp)  ! exner nh pressure interfaces
  real (kind=real_kind), intent(out), optional :: pnh_i_out(np,np,nlevp)  ! pnh on interfaces

  !   local
  real (kind=real_kind) :: kappa_star_i(np,np,nlev)
  real (kind=real_kind) :: pnh_i(np,np,nlevp) 
  real (kind=real_kind) :: exner_i(np,np,nlevp) 
  real (kind=real_kind) :: rho_R_theta(np,np,nlevp) 
  integer :: k

  if (theta_hydrostatic_mode) then
     ! hydrostatic pressure
     pnh_i(:,:,1)=hvcoord%hyai(1)*hvcoord%ps0
     do k=1,nlev
        pnh_i(:,:,k+1)=pnh_i(:,:,k) + dp3d(:,:,k)
     enddo
     do k=1,nlev
        pnh(:,:,k)=pnh_i(:,:,k) + dp3d(:,:,k)/2
     enddo
     exner  = (pnh/p0)**kappa_star
     dpnh = dp3d

     if (present(exner_i_out)) then
        do k=2,nlev
           kappa_star_i(:,:,k) = (kappa_star(:,:,k)+kappa_star(:,:,k-1))/2
           exner_i_out(:,:,k) = (pnh_i(:,:,k)/p0)**kappa_star(:,:,k)
        enddo
        ! how to approximate kappa_star at top and bottom of model?  
        exner_i_out(:,:,1) = (pnh_i(:,:,1)/p0)**kappa_star(:,:,1)
        exner_i_out(:,:,nlev+1) = (pnh_i(:,:,nlev+1)/p0)**kappa_star(:,:,nlev)
     endif
     if (present(pnh_i_out)) pnh_i_out=pnh_i
  else
!==============================================================
!  non-hydrostatic EOS
!==============================================================
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=2,nlev
     kappa_star_i(:,:,k) = (kappa_star(:,:,k)+kappa_star(:,:,k-1))/2

     rho_R_theta(:,:,k) =0.5* (theta_dp_cp(:,:,k)*kappa_star(:,:,k)+ &
                          theta_dp_cp(:,:,k-1)*kappa_star(:,:,k-1))/&
                          (phi(:,:,k-1)-phi(:,:,k))


     if (minval(rho_R_theta(:,:,k))<0) then
        print *,k,minval( (dp3d(:,:,k)+dp3d(:,:,k-1))/2)
        print *,'phi(k-1)-phi(k)',minval(phi(:,:,k-1)-phi(:,:,k))
        call abortmp('error: rho<0')
     endif
    
     ! theta = T/e
     ! exner = (p/p0)**kappa         p = p0*exner**(1/kappa)
     ! p/exner = rho* Rstar * theta 
     ! use this formula which only has 1 exponential:
     exner_i(:,:,k) = (rho_R_theta(:,:,k)/p0)**&
          ( kappa_star_i(:,:,k)/ ( 1-kappa_star_i(:,:,k)))
     pnh_i(:,:,k) = rho_R_theta(:,:,k)*exner_i(:,:,k)
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! boundary terms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   pnh_i(:,:,1) = hvcoord%hyai(1)*hvcoord%ps0   ! hydrostatic ptop
   exner_i(:,:,1) = (hvcoord%hyai(1)*hvcoord%ps0/p0)**kappa_star(:,:,1) 

   ! d_eta at surface = d_eta(nlev), so d_eta scales out below:
   rho_R_theta(:,:,nlev) = theta_dp_cp(:,:,nlev)*kappa_star(:,:,nlev) &
        /  (  (phi(:,:,nlev) - phis(:,:))*2  )
   exner_i(:,:,nlev+1) = (rho_R_theta(:,:,nlev)/p0)**&
        ( kappa_star(:,:,nlev)/ ( 1-kappa_star(:,:,nlev)))
   pnh_i(:,:,nlev+1) = rho_R_theta(:,:,nlev)*exner_i(:,:,nlev+1)


#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     dpnh(:,:,k)=pnh_i(:,:,k+1)-pnh_i(:,:,k)
     exner(:,:,k) = (exner_i(:,:,k)+exner_i(:,:,k+1))/2
     pnh(:,:,k) = p0*exner(:,:,k)**(1/kappa_star(:,:,k))
  enddo

  if (present(exner_i_out))  exner_i_out = exner_i 
  if (present(pnh_i_out)) pnh_i_out=pnh_i


  endif ! hydrostatic/nonhydrostatic version
  end subroutine get_pnh_and_exner



  !_____________________________________________________________________
  subroutine get_dry_phinh(hvcoord,phis,theta_dp_cp,dp,phi)
!
! Use Equation of State to compute geopotential
!
! with  ep = d(exner)/dp = kappa*exner/p            
!
! input:  dp, phis, theta_dp_cp  (assumes constant kappa)
! output:  phi
!
! used to initialize phi for dry test cases
! used to compute background phi for reference state
! in both the above uses, we can assume dry and we dont need to compute kappa_star
!
! NOTE1: dp is pressure layer thickness.  If pnh is used to compute thickness, this
! routine should be the discrete inverse of get_pnh_and_exner().
! This routine is usually called with hydrostatic layer thickness (dp3d), 
! in which case it returns a hydrostatic PHI
!
! NOTE2: Exner pressure is defined in terms of p0=1000mb.  Be sure to use global constant p0,
! instead of hvcoord%ps0, which is set by CAM to ~1021mb
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  
  type (hvcoord_t),      intent(in)  :: hvcoord                      ! hybrid vertical coordinate struct
  real (kind=real_kind), intent(in) :: theta_dp_cp(np,np,nlev)
  real (kind=real_kind), intent(in) :: dp(np,np,nlev)
  real (kind=real_kind), intent(in) :: phis(np,np)
  real (kind=real_kind), intent(out) :: phi(np,np,nlev)
  
  !   local
!  real (kind=real_kind) :: p(np,np,nlev)
  real (kind=real_kind) :: p_i(np,np,nlev+1)
  real (kind=real_kind) :: rho_R_theta(np,np,nlev)
  integer :: k

  p_i(:,:,1) =  hvcoord%hyai(1)*hvcoord%ps0
  do k=1,nlev
     p_i(:,:,k+1) = p_i(:,:,k) + dp(:,:,k)
  enddo

!  integrand(:,:) = dp(:,:,nlev)*Rgas*temperature(:,:,nlev)/p(:,:,nlev)
  phi(:,:,nlev) = phis(:,:) + (&
    kappa*theta_dp_cp(:,:,nlev)*p_i(:,:,nlev+1)**(kappa-1)*p0**(-kappa) )/2

  do k=nlev,2,-1
     rho_R_theta(:,:,k) = &
          (theta_dp_cp(:,:,k)*kappa + theta_dp_cp(:,:,k-1)*kappa)/2
     phi(:,:,k-1) = phi(:,:,k) +&
          rho_R_theta(:,:,k) * (p_i(:,:,k)/p0)**kappa / p_i(:,:,k)
  enddo

  end subroutine

  !_____________________________________________________________________
  subroutine get_moist_phinh(hvcoord,phis,theta_dp_cp,dp,kappa_star,phi)
!
! Use Equation of State to compute geopotential
!
! with  ep = d(exner)/dp = kappa*exner/p            
!
! input:  dp, kappa_star, phis, theta_dp_cp
! output:  phi
!
! used to initialize phi for wet  test cases
!
! NOTE1: dp is pressure layer thickness.  If pnh is used to compute thickness, this
! routine should be the discrete inverse of get_pnh_and_exner().
! This routine is usually called with hydrostatic layer thickness (dp3d), 
! in which case it returns a hydrostatic PHI
!
! NOTE2: Exner pressure is defined in terms of p0=1000mb.  Be sure to use global constant p0,
! instead of hvcoord%ps0, which is set by CAM to ~1021mb
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  
  type (hvcoord_t),      intent(in)   :: hvcoord                      ! hybrid vertical coordinate struct
  real (kind=real_kind), intent(in)   :: theta_dp_cp(np,np,nlev)
  real (kind=real_kind), intent(in)   :: dp(np,np,nlev)
  real (kind=real_kind), intent(in)   :: phis(np,np)
  real (kind=real_kind), intent(in)   :: kappa_star(np,np,nlev)
  real (kind=real_kind), intent(out)  :: phi(np,np,nlev)
  
  real (kind=real_kind) :: p_i(np,np,nlev+1)
  real (kind=real_kind) :: rho_R_theta(np,np,nlev)
  real (kind=real_kind) :: kappa_star_i(np,np,nlevp)

  integer :: k

  p_i(:,:,1) =  hvcoord%hyai(1)*hvcoord%ps0
  do k=1,nlev
     p_i(:,:,k+1) = p_i(:,:,k) + dp(:,:,k)
  enddo

  ! integrand(:,:) = dp(:,:,nlev)*Rgas*temperature(:,:,nlev)/p(:,:,nlev)
  phi(:,:,nlev) = phis(:,:) + ( kappa_star(:,:,nlev)*theta_dp_cp(:,:,nlev)*p_i(:,:,nlev+1)**(kappa_star(:,:,nlev)-1)*p0**(-kappa_star(:,:,nlev)) )/2

  do k=nlev,2,-1
     rho_R_theta(:,:,k) = &
          (theta_dp_cp(:,:,k)*kappa_star(:,:,k) + theta_dp_cp(:,:,k-1)*kappa_star(:,:,k-1))/2

     kappa_star_i(:,:,k) = (kappa_star(:,:,k)+kappa_star(:,:,k-1))/2

     phi(:,:,k-1) = phi(:,:,k) +&
          rho_R_theta(:,:,k) * (p_i(:,:,k)/p0)**kappa_star_i(:,:,k) / p_i(:,:,k)
  enddo

  end subroutine


  subroutine get_dirk_jacobian(JacL,JacD,JacU,dt2,dp3d,phi,phis,kappa_star_i,pnh_i,exact,epsie,hvcoord,dpnh_dp, &
    theta_dp_cp,kappa_star,pnh,exner)
  !================================================================================
  ! This subroutine forms the tridiagonal analytic Jacobian (we actually form the diagonal, sub-, and super-diagonal)
  ! J for use in a LApack tridiagonal LU factorization and solver to solve  J * x = -f either exactly or
  ! approximately
  !
  ! epsie == 1 means exact Jacobian, epsie ~= 1 means finite difference approximate jacobian
  ! exact,epsie,hvcoord,dpnh_dp,theta_dp_cp,kappa_star,pnh,exner,exner_i are only needed as inputs
  ! if epsie ~=1
  !  
  ! The rule-of-thumb optimal epsie  is epsie = norm(elem)*sqrt(macheps)
  !===================================================================================
    real (kind=real_kind), intent(inout) :: JacD(nlev,np,np), pnh_i(np,np,nlevp)
    real (kind=real_kind), intent(inout) :: JacL(nlev-1,np,np),JacU(nlev-1,np,np), phi(np,np,nlev)
    real (kind=real_kind), intent(in)    :: dp3d(np,np,nlev), phis(np,np)
    real (kind=real_kind), intent(in) :: kappa_star_i(np,np,nlevp)
    real (kind=real_kind), intent(in)    :: dt2

    real (kind=real_kind), intent(in), optional :: epsie ! epsie is the differencing size in the approx. Jacobian
    real (kind=real_kind), intent(inout),  optional :: dpnh_dp(np,np,nlev), exner(np,np,nlev)
    real (kind=real_kind), intent(inout),  optional :: kappa_star(np,np,nlev),theta_dp_cp(np,np,nlev), pnh(np,np,nlev)
    type (hvcoord_t)     , intent(in), optional :: hvcoord


    integer, intent(in) :: exact

    ! local
    real (kind=real_kind) :: alpha1(np,np),alpha2(np,np)
    real (kind=real_kind) :: e(np,np,nlev),phitemp(np,np,nlev)
    real (kind=real_kind) :: dpnh2(np,np,nlev),dpnh_dpepsie(np,np,nlev)
    !
    integer :: k,l

    if (exact.eq.1) then ! use exact Jacobian
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,alpha1,alpha2)
#endif
      do k=1,nlev
        ! this code will need to change when the equation of state is changed.
        ! precompute the kappa_star, and add special cases for k==1 and k==nlev+1
        if (k==1) then
          alpha2(:,:)    = 1.d0 + kappa_star_i(:,:,k+1)/(1.d0-kappa_star_i(:,:,k+1))
          JacD(k,:,:)    = (dt2*g)**2.d0 *alpha2(:,:)*pnh_i(:,:,k+1)/((phi(:,:,k)-phi(:,:,k+1))* &
            dp3d(:,:,k))+1.d0
          JacU(k,:,:)    = (dt2*g)**2.d0 *alpha2(:,:)*pnh_i(:,:,k+1)/((phi(:,:,k+1)-phi(:,:,k)) &
            *dp3d(:,:,k))
	    elseif (k==nlev) then
          alpha1(:,:)    = 1.d0 + kappa_star_i(:,:,k)/(1.d0-kappa_star_i(:,:,k))
          JacL(k-1,:,:)  = (dt2*g)**2.d0 *(alpha1(:,:)*pnh_i(:,:,k)/((phi(:,:,k)-phi(:,:,k-1))*dp3d(:,:,k)))
          JacD(k,:,:)    = (dt2*g)**2.d0 *(  alpha1(:,:)*pnh_i(:,:,k+1)/(phi(:,:,k)-phis(:,:) ) +  &
            alpha1(:,:)*pnh_i(:,:,k)/(phi(:,:,k-1)-phi(:,:,k)) )/dp3d(:,:,k) + 1.d0
        else
          alpha1(:,:)   = 1.d0 + kappa_star_i(:,:,k)/(1.d0-kappa_star_i(:,:,k))
          alpha2(:,:)   = 1.d0 + kappa_star_i(:,:,k+1)/(1.d0-kappa_star_i(:,:,k+1))
          JacL(k-1,:,:) = (dt2*g)**2.d0 *alpha1(:,:)*pnh_i(:,:,k)/((phi(:,:,k)-phi(:,:,k-1))*dp3d(:,:,k))
          JacD(k,:,:)   = (dt2*g)**2.d0 *(alpha2(:,:)*pnh_i(:,:,k+1)/(phi(:,:,k)-phi(:,:,k+1)) + &
            alpha1(:,:)*pnh_i(:,:,k)/(phi(:,:,k-1)-phi(:,:,k)))/dp3d(:,:,k)+1.d0
          JacU(k,:,:)   = (dt2*g)**2.d0 *(alpha2(:,:)*pnh_i(:,:,k+1)/((phi(:,:,k+1)-phi(:,:,k))*dp3d(:,:,k)))
        end if
      end do
    else ! use finite difference approximation to Jacobian with differencing size espie
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,phitemp,dpnh_dpepsie)
#endif
      ! compute Jacobian of F(phi) = phi +const + (dt*g)^2 *(1-dp/dpi) column wise
      ! we only form the tridagonal entries and this code can easily be modified to
      ! accomodate sparse non-tridigonal and dense Jacobians, however, testing only
      ! the tridiagonal of a Jacobian is probably sufficient for testing purpose
      do k=1,nlev
        e=0.d0
        e(:,:,k)=1.d0
        phitemp(:,:,:)=phi(:,:,:)
        phitemp(:,:,k) = phi(:,:,k)+epsie*e(:,:,k)
        if (theta_hydrostatic_mode) then
          dpnh_dpepsie(:,:,:)=1.d0
        else
          call get_pnh_and_exner(hvcoord,theta_dp_cp,dp3d,phitemp,phis,&
            kappa_star,pnh,dpnh_dpepsie,exner,pnh_i_out=pnh_i)
          dpnh_dpepsie(:,:,:)=dpnh_dpepsie(:,:,:)/dp3d(:,:,:)
        end if
        if (k.eq.1) then
          JacL(k,:,:) = (g*dt2)**2.d0 * (dpnh_dp(:,:,k+1)-dpnh_dpepsie(:,:,k+1))/epsie
          JacD(k,:,:) = 1.d0 + (g*dt2)**2.d0 * (dpnh_dp(:,:,k)-dpnh_dpepsie(:,:,k))/epsie
        elseif (k.eq.nlev) then
          JacD(k,:,:)   = 1.d0 + (g*dt2)**2.d0 * (dpnh_dp(:,:,k)-dpnh_dpepsie(:,:,k))/epsie
          JacU(k-1,:,:) = (g*dt2)**2.d0 * (dpnh_dp(:,:,k-1)-dpnh_dpepsie(:,:,k-1))/epsie
        else
          JacL(k,:,:)   = (g*dt2)**2.d0 * (dpnh_dp(:,:,k+1)-dpnh_dpepsie(:,:,k+1))/epsie
          JacD(k,:,:)   = 1.d0 + (g*dt2)**2.d0 * (dpnh_dp(:,:,k)-dpnh_dpepsie(:,:,k))/epsie
          JacU(k-1,:,:) = (g*dt2)**2.d0 * (dpnh_dp(:,:,k-1)-dpnh_dpepsie(:,:,k-1))/epsie
        end if
      end do
    end if


  end subroutine get_dirk_jacobian



end module

