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

  use dimensions_mod, only: np, nlev, nlevp
  use element_mod,    only: element_t
  use element_state,  only: timelevels, elem_state_t
  use hybvcoord_mod,  only: hvcoord_t
  use kinds,          only: real_kind
  use parallel_mod,   only: abortmp
  use physical_constants, only : p0, kappa, g, Rgas
  use control_mod,    only: theta_hydrostatic_mode
  implicit none


contains

subroutine get_pnh_and_exner(hvcoord,vtheta_dp,dp3d,phi_i,pnh,exner,&
     dpnh_dp_i,pnh_i_out)
implicit none
!
! Use Equation of State to compute exner pressure, nh presure
! hydrostatic EOS:
!          compute p, exner  
!
! nonhydrostatic EOS:
!      p_over_exner   =  -R  vtheta_dp / (dphi/ds)
!
! input:  dp3d, phi, phis, vtheta_dp
! output:  pnh, dphn, exner, exner_i, pnh_i
!
! NOTE: Exner pressure is defined in terms of p0=1000mb.  Be sure to use global constant p0,
! instead of hvcoord%ps0, which is set by CAM to ~1021mb
!  
  type (hvcoord_t),     intent(in)  :: hvcoord             ! hybrid vertical coordinate struct
  real (kind=real_kind), intent(in) :: vtheta_dp(np,np,nlev)   
  real (kind=real_kind), intent(in) :: dp3d(np,np,nlev)   
  real (kind=real_kind), intent(in) :: phi_i(np,np,nlevp)
  real (kind=real_kind), intent(out) :: pnh(np,np,nlev)        ! nh nonhyrdo pressure
  real (kind=real_kind), intent(out) :: dpnh_dp_i(np,np,nlevp) ! d(pnh) / d(pi)
  real (kind=real_kind), intent(out) :: exner(np,np,nlev)      ! exner nh pressure
  real (kind=real_kind), intent(out), optional :: pnh_i_out(np,np,nlevp)  ! pnh on interfaces

  !   local
  real (kind=real_kind) :: p_over_exner(np,np,nlev)
  real (kind=real_kind) :: pi(np,np,nlev)
  real (kind=real_kind) :: exner_i(np,np,nlevp) 
  real (kind=real_kind) :: pnh_i(np,np,nlevp)  
  real (kind=real_kind) :: dp3d_i(np,np,nlevp)
  real (kind=real_kind) :: pi_i(np,np,nlevp) 
  integer :: i,j,k,k2

  ! hydrostatic pressure
  pi_i(:,:,1)=hvcoord%hyai(1)*hvcoord%ps0
  do k=1,nlev
     pi_i(:,:,k+1)=pi_i(:,:,k) + dp3d(:,:,k)
  enddo
  do k=1,nlev
     pi(:,:,k)=pi_i(:,:,k) + dp3d(:,:,k)/2
  enddo


  if (theta_hydrostatic_mode) then
     ! hydrostatic pressure
     exner  = (pi/p0)**kappa
     pnh = pi ! copy hydrostatic pressure into output variable
     dpnh_dp_i = 1
     if (present(pnh_i_out)) then  
       pnh_i_out=pi_i 
     endif
  else

!==============================================================
!  non-hydrostatic EOS
!==============================================================
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     p_over_exner(:,:,k) = Rgas*vtheta_dp(:,:,k)/(phi_i(:,:,k)-phi_i(:,:,k+1)) 

     if (minval(p_over_exner(:,:,k))<0) then
        do i=1,np
           do j=1,np
              if ( p_over_exner(i,j,k)<0 ) then
                 print *,'i,j,k, p/exner = ',i,j,k,p_over_exner(i,j,k)
                 print *,'k,phi_i(k),phi_i(k+1):',k,(phi_i(i,j,k)),(phi_i(i,j,k+1))
              endif
           enddo
        enddo
        do i=1,np
           do j=1,np
              if ( p_over_exner(i,j,k)<0 ) then
                 print *,'vertical column:'
                 do k2=1,nlev
                    write(*,'(i3,4f14.4)') k2,phi_i(i,j,k2),dp3d(i,j,k2),vtheta_dp(i,j,k2)
                 enddo
                 call abortmp('error: rho<0')
              endif
           enddo
        enddo
     endif
    
     pnh(:,:,k) = p0 * (p_over_exner(:,:,k)/p0)**(1/(1-kappa))
     exner(:,:,k) =  pnh(:,:,k)/ p_over_exner(:,:,k)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! boundary terms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   pnh_i(:,:,1) = pi_i(:,:,1)   ! hydrostatic ptop    
   ! surface boundary condition pnh_i determined by w equation to enforce
   ! w b.c.  This is computed in the RHS calculation.  Here, we use
   ! an approximation (hydrostatic) so that dpnh/dpi = 1
   pnh_i(:,:,nlevp) = pnh(:,:,nlev) + dp3d(:,:,nlev)/2
   ! extrapolote NH perturbation:
   !pnh_i(:,:,nlevp) = pi_i(:,:,nlevp) + (3*(pnh(:,:,nlev)-pi(:,:,nlev)) - (pnh(:,:,nlev-1)-pi(:,:,nlev-1)) )/2


   ! compute d(pnh)/d(pi) at interfaces
   ! use one-sided differences at boundaries
   dp3d_i(:,:,1) = dp3d(:,:,1)
   dp3d_i(:,:,nlevp) = dp3d(:,:,nlev)
   do k=2,nlev
      dp3d_i(:,:,k)=(dp3d(:,:,k)+dp3d(:,:,k-1))/2
   end do

   dpnh_dp_i(:,:,1)  = 2*(pnh(:,:,1)-pnh_i(:,:,1))/dp3d_i(:,:,1)
   dpnh_dp_i(:,:,nlevp)  = 2*(pnh_i(:,:,nlevp)-pnh(:,:,nlev))/dp3d_i(:,:,nlevp)
   do k=2,nlev
      dpnh_dp_i(:,:,k) = (pnh(:,:,k)-pnh(:,:,k-1))/dp3d_i(:,:,k)        
   end do
   

   if (present(pnh_i_out)) then
      ! boundary values already computed.  interior only:
      do k=2,nlev
         pnh_i(:,:,k)=(hvcoord%d_etam(k)*pnh(:,:,k)+hvcoord%d_etam(k-1)*pnh(:,:,k-1))/&
              (hvcoord%d_etai(k)*2)
      enddo
      pnh_i_out=pnh_i    
   endif
   
  endif ! hydrostatic/nonhydrostatic version
  end subroutine get_pnh_and_exner



subroutine get_theta_from_T(hvcoord,Rstar,temperature,dp3d,phi_i,vtheta_dp)
implicit none
!
! Use Equation of State to compute vtheta_dp from Temperature
!
!      p_over_exner   =  -Rstar T dp / (dphi/ds)
!
! input:  dp3d, phi, phis, T, Rstar
! output:  vtheta_dp
!
  type (hvcoord_t),     intent(in)  :: hvcoord             ! hybrid vertical coordinate struct
  real (kind=real_kind), intent(in) :: Rstar(np,np,nlev)   
  real (kind=real_kind), intent(in) :: temperature(np,np,nlev)   
  real (kind=real_kind), intent(in) :: dp3d(np,np,nlev)   
  real (kind=real_kind), intent(in) :: phi_i(np,np,nlevp)
  real (kind=real_kind), intent(out) :: vtheta_dp(np,np,nlev)   

  !   local
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: pi(np,np,nlev)
  real (kind=real_kind) :: pi_i(np,np,nlevp)
  integer :: k

  if (theta_hydrostatic_mode) then
     ! hydrostatic pressure
     pi_i(:,:,1)=hvcoord%hyai(1)*hvcoord%ps0
     do k=1,nlev
        pi_i(:,:,k+1)=pi_i(:,:,k) + dp3d(:,:,k)
     enddo
     do k=1,nlev
        pi(:,:,k)=pi_i(:,:,k) + dp3d(:,:,k)/2
     enddo
     
     ! hydrostatic pressure
     exner  = (pi/p0)**kappa
  else

     do k=1,nlev
        pi(:,:,k) = Rstar(:,:,k)*temperature(:,:,k)*dp3d(:,:,k)/ &
          (phi_i(:,:,k)-phi_i(:,:,k+1))
     enddo
     exner = (pi/p0)**kappa

  endif ! hydrostatic/nonhydrostatic version
   
  vtheta_dp = &
       (Rstar(:,:,:)/Rgas)*temperature(:,:,:)*dp3d(:,:,:)/exner(:,:,:)

  end subroutine get_theta_from_T



  !_____________________________________________________________________
  subroutine get_phinh(hvcoord,phis,vtheta_dp,dp,phi_i)
!
! Use Equation of State to compute geopotential
!
! input:  dp, phis, vtheta_dp  
! output:  phi
!
! used to initialize phi for dry and wet test cases
! used to compute background phi for reference state
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
  real (kind=real_kind), intent(in) :: vtheta_dp(np,np,nlev)
  real (kind=real_kind), intent(in) :: dp(np,np,nlev)
  real (kind=real_kind), intent(in) :: phis(np,np)
  real (kind=real_kind), intent(out) :: phi_i(np,np,nlevp)
 
  !   local
  real (kind=real_kind) :: p(np,np,nlev) ! pressure at cell centers 
  real (kind=real_kind) :: p_i(np,np,nlevp)  ! pressure on interfaces

  integer :: k

  ! compute pressure on interfaces                                                                                   
  p_i(:,:,1)=hvcoord%hyai(1)*hvcoord%ps0
  do k=1,nlev
     p_i(:,:,k+1)=p_i(:,:,k) + dp(:,:,k)
  enddo
  do k=1,nlev
     p(:,:,k) = (p_i(:,:,k+1)+p_i(:,:,k))/2
  enddo
 
  phi_i(:,:,nlevp) = phis(:,:)
  do k=nlev,1,-1
     phi_i(:,:,k) = phi_i(:,:,k+1)+(Rgas*vtheta_dp(:,:,k)*(p(:,:,k)/p0)**(kappa-1))/p0
  enddo
  end subroutine

  

  subroutine get_dirk_jacobian(JacL,JacD,JacU,dt2,dp3d,phi_i,pnh,exact,&
       epsie,hvcoord,dpnh_dp_i,vtheta_dp)
  !================================================================================
  ! compute Jacobian of F(phi) = phi +const + (dt*g)^2 *(1-dp/dpi) column wise
  ! with respect to phi_i 
  !
  ! This subroutine forms the tridiagonal analytic Jacobian (we actually form the diagonal, sub-, and super-diagonal)
  ! J for use in a LApack tridiagonal LU factorization and solver to solve  J * x = -f either exactly or
  ! approximately
  !
  !  input:  
  !  exact jacobian:        phi_i, dp3d, pnh
  !  matrix free jacobian:  phi_i, dp3d, vtheta_dp, hvcoord, dpnh_dp_i
  !
  ! epsie == 1 means exact Jacobian, epsie ~= 1 means finite difference approximate jacobian
  ! exact,epsie,hvcoord,dpnh_dp,vtheta_dp,pnh,exner,exner_i are only needed as inputs
  ! if epsie ~=1
  !  
  ! The rule-of-thumb optimal epsie  is epsie = norm(elem)*sqrt(macheps)
  !===================================================================================
    real (kind=real_kind), intent(out) :: JacD(nlev,np,np)
    real (kind=real_kind), intent(out) :: JacL(nlev-1,np,np),JacU(nlev-1,np,np)
    real (kind=real_kind), intent(in)    :: dp3d(np,np,nlev), phi_i(np,np,nlevp)
    real (kind=real_kind), intent(inout) :: pnh(np,np,nlev)
    real (kind=real_kind), intent(in)    :: dt2

    real (kind=real_kind), intent(in), optional :: epsie ! epsie is the differencing size in the approx. Jacobian
    real (kind=real_kind), intent(in), optional :: dpnh_dp_i(np,np,nlevp)
    real (kind=real_kind), intent(in), optional :: vtheta_dp(np,np,nlev)
    type (hvcoord_t)     , intent(in),  optional    :: hvcoord

    integer, intent(in) :: exact

    ! local
    real (kind=real_kind) :: alpha1(np,np),alpha2(np,np)
    real (kind=real_kind) :: e(np,np,nlev),phi_i_temp(np,np,nlevp),exner(np,np,nlev)
    real (kind=real_kind) :: dpnh2(np,np,nlev),dpnh_dp_i_epsie(np,np,nlevp)
    real (kind=real_kind) :: dp3d_i(np,np,nlevp)
    !
    integer :: k,l
    if (exact.eq.1) then ! use exact Jacobian
 
       dp3d_i(:,:,1) = dp3d(:,:,1)
       dp3d_i(:,:,nlevp) = dp3d(:,:,nlev)
       do k=2,nlev
          dp3d_i(:,:,k)=(dp3d(:,:,k)+dp3d(:,:,k-1))/2
       end do
  
      do k=1,nlev
        ! this code will need to change when the equation of state is changed.
        ! add special cases for k==1 and k==nlev+1
        if (k==1) then

           JacL(k,:,:) = -(dt2*g)**2*pnh(:,:,k)/&
             ((phi_i(:,:,k)-phi_i(:,:,k+1))*(1-kappa)*dp3d_i(:,:,k+1))

           JacU(k,:,:) = -2*(dt2*g)**2 * pnh(:,:,k)/&
             ((phi_i(:,:,k)-phi_i(:,:,k+1))*(1-kappa)*dp3d_i(:,:,k))

           JacD(k,:,:) = 1+2*(dt2*g)**2 *pnh(:,:,k)/&
             ((phi_i(:,:,k)-phi_i(:,:,k+1))*(1-kappa)*dp3d_i(:,:,k))
          
        else if (k.eq.nlev) then 

           JacD(k,:,:) = 1+(dt2*g)**2 *(pnh(:,:,k)/((phi_i(:,:,k)-phi_i(:,:,k+1))*(1-kappa)) &
             +pnh(:,:,k-1)/( (phi_i(:,:,k-1)-phi_i(:,:,k))*(1-kappa)))/dp3d_i(:,:,k)

        else ! k =2,...,nlev-1

           JacL(k,:,:) = -(dt2*g)**2*pnh(:,:,k)/&
             ((phi_i(:,:,k)-phi_i(:,:,k+1))*(1-kappa)*dp3d_i(:,:,k+1))
         
           JacU(k,:,:) = -(dt2*g)**2 * pnh(:,:,k)/&
             ((phi_i(:,:,k)-phi_i(:,:,k+1))*(1-kappa)*dp3d_i(:,:,k))

           JacD(k,:,:) = 1+(dt2*g)**2 *(pnh(:,:,k)/((phi_i(:,:,k)-phi_i(:,:,k+1))*(1-kappa)) &
             +pnh(:,:,k-1)/( (phi_i(:,:,k-1)-phi_i(:,:,k))*(1-kappa)))/dp3d_i(:,:,k)

        end if
      end do
    else ! use finite difference approximation to Jacobian with differencing size espie
      ! compute Jacobian of F(phi) = phi +const + (dt*g)^2 *(1-dp/dpi) column wise
      ! we only form the tridagonal entries and this code can easily be modified to
      ! accomodate sparse non-tridigonal and dense Jacobians, however, testing only
      ! the tridiagonal of a Jacobian is probably sufficient for testing purpose
      do k=1,nlev
        e=0
        e(:,:,k)=1
        phi_i_temp(:,:,:) = phi_i(:,:,:)
        phi_i_temp(:,:,k) = phi_i(:,:,k) + epsie*e(:,:,k)
        if (theta_hydrostatic_mode) then
          dpnh_dp_i_epsie(:,:,:)=1.d0
        else
         call get_pnh_and_exner(hvcoord,vtheta_dp,dp3d,phi_i_temp,pnh,exner,dpnh_dp_i_epsie)
        end if
        if (k.eq.1) then
          JacL(k,:,:) = (g*dt2)**2 * (dpnh_dp_i(:,:,k+1)-dpnh_dp_i_epsie(:,:,k+1))/epsie
          JacD(k,:,:) = 1 + (g*dt2)**2 * (dpnh_dp_i(:,:,k)-dpnh_dp_i_epsie(:,:,k))/epsie
        elseif (k.eq.nlev) then
          JacD(k,:,:)   = 1 + (g*dt2)**2 * (dpnh_dp_i(:,:,k)-dpnh_dp_i_epsie(:,:,k))/epsie
          JacU(k-1,:,:) = (g*dt2)**2 * (dpnh_dp_i(:,:,k-1)-dpnh_dp_i_epsie(:,:,k-1))/epsie
        else
          JacL(k,:,:)   = (g*dt2)**2 * (dpnh_dp_i(:,:,k+1)-dpnh_dp_i_epsie(:,:,k+1))/epsie
          JacD(k,:,:)   = 1 + (g*dt2)**2 * (dpnh_dp_i(:,:,k)-dpnh_dp_i_epsie(:,:,k))/epsie
          JacU(k-1,:,:) = (g*dt2)**2 * (dpnh_dp_i(:,:,k-1)-dpnh_dp_i_epsie(:,:,k-1))/epsie
        end if
      end do
    end if

  end subroutine get_dirk_jacobian



end module

