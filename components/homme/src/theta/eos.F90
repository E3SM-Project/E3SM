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
!
module eos

  use dimensions_mod, only: np, nlev, nlevp, nelemd
  use element_mod,    only: element_t
  use element_state,  only: timelevels, elem_state_t
  use hybvcoord_mod,  only: hvcoord_t
  use kinds,          only: real_kind
  use parallel_mod,   only: abortmp
  use physical_constants, only : p0, kappa
  use control_mod,    only: theta_hydrostatic_mode
  use prim_si_mod,    only: preq_hydrostatic_v2, preq_omega_ps
  implicit none


contains


  subroutine get_pnh_and_exner(hvcoord,theta_dp_cp,dp3d,phi,phis,kappa_star,pnh,dpnh,exner,exner_i_out,pnh_i_out)
  implicit none
!
! Use Equation of State to compute exner pressure, nh presure
! EOS:   
!
!       theta_dp_cp * ep = -dphi/ds    OR    rho_R_theta * dphi/ds = -theta_dp_cp*kappa
!
! with  ep = d(exner)/dp = kappa*exner/p            
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


end module

