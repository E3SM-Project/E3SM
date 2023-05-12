#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
!
!  NonHydrostatic Equation of State and inverse EOS routines
!  Note: these routines should all be discrete inverses of each other
!
!  pnh_and_exner_from_eos()  Compute nonhydrostatic pressure as a function of 
!                            potential temperature and geopotential
!
!  phi_from_eos()            Compute geopotential as a function of potential temperature
!                            and pressure. use virtual potential temperature for wet phi
!
!  Original version: Mark Taylor 2017/1
!  
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
#ifdef HOMMEXX_BFB_TESTING
  use bfb_mod,        only: bfb_pow
#endif
  implicit none


contains

subroutine pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_i,pnh,exner,&
     dpnh_dp_i,pnh_i_out,caller)
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
  character(len=*),      intent(in), optional  :: caller       ! name for error

  !   local
  real (kind=real_kind) :: dphi(np,np,nlev)
  integer :: k

  do k=1,nlev
     dphi(:,:,k)=phi_i(:,:,k+1)-phi_i(:,:,k)
  enddo
  if (present(caller)) then
     call pnh_and_exner_from_eos2(hvcoord,vtheta_dp,dp3d,dphi,pnh,exner,&
          dpnh_dp_i,caller,pnh_i_out)
  else
     call pnh_and_exner_from_eos2(hvcoord,vtheta_dp,dp3d,dphi,pnh,exner,&
          dpnh_dp_i,'not specified',pnh_i_out)
  endif
  end subroutine pnh_and_exner_from_eos



subroutine pnh_and_exner_from_eos2(hvcoord,vtheta_dp,dp3d,dphi,pnh,exner,&
     dpnh_dp_i,caller,pnh_i_out)
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
  real (kind=real_kind), intent(in) :: dphi(np,np,nlev)
  real (kind=real_kind), intent(out) :: pnh(np,np,nlev)        ! nh nonhyrdo pressure
  real (kind=real_kind), intent(out) :: dpnh_dp_i(np,np,nlevp) ! d(pnh) / d(pi)
  real (kind=real_kind), intent(out) :: exner(np,np,nlev)      ! exner nh pressure
  character(len=*),      intent(in)  :: caller       ! name for error
  real (kind=real_kind), intent(out), optional :: pnh_i_out(np,np,nlevp)  ! pnh on interfaces


  !   local
  real (kind=real_kind) :: p_over_exner(np,np,nlev)
  real (kind=real_kind) :: pi(np,np,nlev)
  real (kind=real_kind) :: exner_i(np,np,nlevp) 
  real (kind=real_kind) :: pnh_i(np,np,nlevp)  
  real (kind=real_kind) :: dp3d_i(np,np,nlevp)
  real (kind=real_kind) :: pi_i(np,np,nlevp) 
  integer :: i,j,k,k2
  logical :: ierr


  ! check for bad state that will crash exponential function below
  if (theta_hydrostatic_mode) then
    ierr= any(dp3d(:,:,:) < 0 )
  else
    ierr= any(vtheta_dp(:,:,:) < 0 )  .or. &
          any(dp3d(:,:,:) < 0 ) .or. &
          any(dphi(:,:,:) > 0 )
  endif

  if (ierr) then
     print *,'bad state in EOS, called from: ',caller
     do j=1,np
     do i=1,np
     do k=1,nlev
        if ( (vtheta_dp(i,j,k) < 0) .or. (dp3d(i,j,k)<0)  .or. &
             (dphi(i,j,k)>0)  ) then
           print *,'bad i,j,k=',i,j,k
           print *,'vertical column: dphi,dp3d,vtheta_dp'
           do k2=1,nlev
              write(*,'(i3,4f14.4)') k2,dphi(i,j,k2),dp3d(i,j,k2),vtheta_dp(i,j,k2)
           enddo
           call abortmp('EOS bad state: d(phi), dp3d or vtheta_dp < 0')
        endif
     enddo
     enddo
     enddo
  endif

  if (theta_hydrostatic_mode) then
     ! hydrostatic pressure
     pi_i(:,:,1)=hvcoord%hyai(1)*hvcoord%ps0
     do k=1,nlev
        pi_i(:,:,k+1)=pi_i(:,:,k) + dp3d(:,:,k)
     enddo
#ifdef HOMMEXX_BFB_TESTING
     do k=1,nlev
        pi(:,:,k) = (pi_i(:,:,k+1)+pi_i(:,:,k))/2
     enddo
     exner  = bfb_pow(pi/p0,kappa)
#else
     do k=1,nlev
        pi(:,:,k)=pi_i(:,:,k) + dp3d(:,:,k)/2
     enddo
     exner  = (pi/p0)**kappa
#endif

     pnh = pi ! copy hydrostatic pressure into output variable
     dpnh_dp_i = 1
     if (present(pnh_i_out)) then  
       pnh_i_out=pi_i 
     endif
  else

!==============================================================
!  non-hydrostatic EOS
!==============================================================
  do k=1,nlev
     p_over_exner(:,:,k) = Rgas*vtheta_dp(:,:,k)/(-dphi(:,:,k))
#ifndef HOMMEXX_BFB_TESTING
     pnh(:,:,k) = p0 * (p_over_exner(:,:,k)/p0)**(1/(1-kappa))
#else
     pnh(:,:,k) = p0 * bfb_pow(p_over_exner(:,:,k)/p0,1/(1-kappa))
#endif
     exner(:,:,k) =  pnh(:,:,k)/ p_over_exner(:,:,k)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! boundary terms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   pnh_i(:,:,1) = hvcoord%hyai(1)*hvcoord%ps0  ! hydrostatic ptop    
   ! surface boundary condition pnh_i determined by w equation to enforce
   ! w b.c.  This is computed in the RHS calculation.  Here, we use
   ! an approximation (hydrostatic) so that dpnh/dpi = 1
   ! DO NOT CHANGE this approximation.  it is required by 
   ! compute_andor_apply_rhs()
   pnh_i(:,:,nlevp) = pnh(:,:,nlev) + dp3d(:,:,nlev)/2


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
      ! boundary values already computed. interpolate interior
      ! use linear interpolation in hydrostatic pressure coordinate
      ! if pnh=pi, then pnh_i will recover pi_i
      do k=2,nlev
         pnh_i(:,:,k)=(dp3d(:,:,k-1)*pnh(:,:,k)+dp3d(:,:,k)*pnh(:,:,k-1))/&
              (dp3d(:,:,k-1)+dp3d(:,:,k))
      enddo
      pnh_i_out=pnh_i    
   endif
   
  endif ! hydrostatic/nonhydrostatic version
  end subroutine 



  !_____________________________________________________________________
  subroutine phi_from_eos(hvcoord,phis,vtheta_dp,dp,phi_i)
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
! routine should be the discrete inverse of pnh_and_exner_from_eos().
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

#ifndef NDEBUG
  logical :: ierr
  integer :: i,j,k2

  ierr= any(vtheta_dp(:,:,:) < 0 )  .or. &
          any(dp(:,:,:) < 0 )

  if (ierr) then
     print *,'bad state in phi_from_eos:'
     do j=1,np
     do i=1,np
     do k=1,nlev
        if ( (vtheta_dp(i,j,k) < 0) .or. (dp(i,j,k)<0) ) then
           print *,'bad i,j,k=',i,j,k
           print *,'vertical column: dp,vtheta_dp'
           do k2=1,nlev
              write(*,'(i3,4f14.4)') k2,dp(i,j,k2),vtheta_dp(i,j,k2)
           enddo
           call abortmp('EOS bad state: dp or vtheta_dp < 0')
        endif
     enddo
     enddo
     enddo
  endif
#endif
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
#ifdef HOMMEXX_BFB_TESTING
     phi_i(:,:,k) = phi_i(:,:,k+1)+ (Rgas*vtheta_dp(:,:,k)*bfb_pow(p(:,:,k)/p0,(kappa-1)))/p0
#else
     phi_i(:,:,k) = phi_i(:,:,k+1)+(Rgas*vtheta_dp(:,:,k)*(p(:,:,k)/p0)**(kappa-1))/p0
#endif
  enddo
  end subroutine

end module

