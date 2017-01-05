#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
!
!  getter and setter functions that must be provided by each model
!  
!  IMPORTANT NOTE:  For vertically lagrangian models, these
!  routines should ONLY be used outside the timestepping loop
!  on reference levels.  To compute these fields on floating levels
!  the model should do that directly (or we need to modify the interface)
!
!  see src/preqx/element_ops.F90 for documentation
!     
!  helper fuctions in here not required to be provided by all models:
!     get_p_[non]hydrostatic
!     get_temperature
!     get_pottemp
!    
!
module element_ops

  use dimensions_mod, only: np, nlev, nlevp, nelemd, qsize, max_corner_elem
  use element_mod,    only: element_t
  use element_state,  only: timelevels
  use hybrid_mod,     only: hybrid_t
  use hybvcoord_mod,  only: hvcoord_t
  use kinds,          only: real_kind, iulog
  use perf_mod,       only: t_startf, t_stopf, t_barrierf, t_adj_detailf ! _EXTERNAL
  use parallel_mod,   only: abortmp
  use physical_constants, only : p0, Cp, Rgas, Rwater_vapor, Cpwater_vapor, kappa
  use control_mod, only:    use_moisture, use_cpstar, theta_hydrostatic_mode
  use prim_si_mod, only: preq_hydrostatic_v2
  implicit none


contains

  subroutine get_field(elem,name,field,hvcoord,nt,ntQ)
  implicit none
  type (element_t), intent(in) :: elem
  character(len=*), intent(in) :: name
  real (kind=real_kind), intent(out)  :: field(np,np,nlev)
  type (hvcoord_t),     intent(in)    :: hvcoord          
  integer, intent(in) :: nt
  integer, intent(in) :: ntQ

  select case(name)
  case ('temperature')
     call get_temperature(elem,field,hvcoord,nt,ntQ)
  case ('pottemp')
     call get_pottemp(elem,field,hvcoord,nt,ntQ)
  case ('phi')
     field = elem%state%phi(:,:,:,nt)
     !field = elem%derived%phi(:,:,:)
  case ('dpnh_dp')
     call get_dpnh_dp(elem,field,hvcoord,nt,ntQ)
  case default
     print *,'name = ',trim(name)
     call abortmp('ERROR: get_field name not supported in this model')
  end select

  end subroutine



  subroutine get_pottemp(elem,pottemp,hvcoord,nt,ntQ)
  implicit none
    
  type (element_t), intent(in)        :: elem
  real (kind=real_kind), intent(out)  :: pottemp(np,np,nlev)
  type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
  integer, intent(in) :: nt
  integer, intent(in) :: ntQ
  
  !   local
  real (kind=real_kind) :: p(np,np,nlev)
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: kappa_star(np,np,nlev)
  real (kind=real_kind) :: Qt(np,np,nlev)
  integer :: k
  
  pottemp(:,:,:) = elem%state%theta(:,:,:,nt)
  
  end subroutine get_pottemp
  


  subroutine get_temperature(elem,temperature,hvcoord,nt,ntQ)
  implicit none
  
  type (element_t), intent(in)        :: elem
  real (kind=real_kind), intent(out)  :: temperature(np,np,nlev)
  type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
  integer, intent(in) :: nt
  integer, intent(in) :: ntQ
  
  !   local
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev)
  real (kind=real_kind) :: dpnh(np,np,nlev)
  real (kind=real_kind) :: pnh_i(np,np,nlevp)
  integer :: k
  
  
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem%state%ps_v(:,:,nt)
  enddo

  if (theta_hydrostatic_mode) then
     call get_p_hydrostatic(pnh,pnh_i,exner,hvcoord,dp,elem%state%Qdp(:,:,:,1,ntQ))
  else
     call get_p_nonhydrostatic(pnh,dpnh,exner,hvcoord,elem%state%theta(:,:,:,nt),&
          dp,elem%state%phi(:,:,:,nt),elem%state%phis,elem%state%Qdp(:,:,:,1,ntQ))
  endif


#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     temperature(:,:,k)= elem%state%theta(:,:,k,nt)*exner(:,:,k)
  enddo

  end subroutine get_temperature





  subroutine get_dpnh_dp(elem,dpnh_dp,hvcoord,nt,ntQ)
  implicit none
  
  type (element_t), intent(in)        :: elem
  real (kind=real_kind), intent(out)  :: dpnh_dp(np,np,nlev)
  type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
  integer, intent(in) :: nt
  integer, intent(in) :: ntQ
  
  !   local
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev)
  real (kind=real_kind) :: dpnh(np,np,nlev)
  real (kind=real_kind) :: pnh_i(np,np,nlevp)
  integer :: k
  
  
  do k=1,nlev
     dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem%state%ps_v(:,:,nt)
  enddo

  if (theta_hydrostatic_mode) then
     dpnh_dp=1
     !call get_p_hydrostatic(pnh,pnh_i,exner,hvcoord,dp,elem%state%Qdp(:,:,:,1,ntQ))
  else
     call get_p_nonhydrostatic(pnh,dpnh,exner,hvcoord,elem%state%theta(:,:,:,nt),&
          dp,elem%state%phi(:,:,:,nt),elem%state%phis,elem%state%Qdp(:,:,:,1,ntQ))
     dpnh_dp = dpnh/dp
  endif
  end subroutine 








  subroutine get_p_hydrostatic(pnh,pnh_i,exner,hvcoord,dp3d,Qdp)
  implicit none
!
! compute exner pressure and nh presure
!
! hydrostatic case:
!    input:  dp3d, Qdp (if use_moisture)
!    output:  pnh = compute from dp3d
!             exner = pnh**kappa_star
!
!
!  
  real (kind=real_kind), intent(out) :: exner(np,np,nlev)  ! exner nh pressure
  real (kind=real_kind), intent(out) :: pnh(np,np,nlev)   ! nh nonhyrdo pressure
  real (kind=real_kind), intent(out) :: pnh_i(np,np,nlevp) ! nh nonhyrdo pressure interfaces
  type (hvcoord_t),     intent(in)  :: hvcoord             ! hybrid vertical coordinate struct
  real (kind=real_kind), intent(in) :: dp3d(np,np,nlev)   
  real (kind=real_kind), intent(in) :: Qdp(np,np,nlev)   

  !   local
  real (kind=real_kind) :: kappa_star(np,np,nlev)
  real (kind=real_kind) :: R_star(np,np,nlev)
  real (kind=real_kind) :: Qt(np,np,nlev)
  real (kind=real_kind) :: ptop
  integer :: k


  ptop = hvcoord%hyai(1)*hvcoord%ps0

#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     if (use_moisture) then
        Qt(:,:,k) = Qdp(:,:,k)/dp3d(:,:,k)
        R_star(:,:,k)=(Rgas + (Rwater_vapor - Rgas)*Qt(:,:,k))
        
        if (use_cpstar==1) then
           kappa_star(:,:,k) = (Rgas + (Rwater_vapor - Rgas)*Qt(:,:,k)) / &
                (Cp + (Cpwater_vapor-Cp)*Qt(:,:,k) )
        else
           kappa_star(:,:,k) = (Rgas + (Rwater_vapor - Rgas)*Qt(:,:,k)) / Cp
        endif
     else
        R_star(:,:,k)=Rgas
        kappa_star(:,:,k)=Rgas/Cp
     endif
  enddo

! hydrostatic model:
  pnh_i(:,:,1)=ptop
  do k=1,nlev
     pnh_i(:,:,k+1)=pnh_i(:,:,k) + dp3d(:,:,k)
  enddo
  
  do k=1,nlev
     pnh(:,:,k)=pnh_i(:,:,k) + dp3d(:,:,k)/2
  enddo
  exner  = (pnh/p0)**kappa_star

  end subroutine 





  subroutine get_p_nonhydrostatic(pnh,dpnh,exner,hvcoord,theta,dp3d,phi,phis,Qdp)
  implicit none
!
! compute exner pressure, nh presure
!
! nonhydro formula used in hydrostatic model (used to debug this routine)
!         call this routine with phi computed from hydrostatic exner
!         
! nonhydro case:
! input:  dp3d, Qdp (if use_moisture), phi, phis, theta
! output:  
!       for k=2..nlev, use the equation of state:  pnh/e = rho*Rstar*theta
!                                                    rho   = -dp3d/dphi 
! for k=1, we cant compute rho because we dont know phi at the model top
! (our boundary condition specifies ptop) so we take p(k=1) = 2/3 ptop + 1/3 p(k=2)
!
!  
  real (kind=real_kind), intent(out) :: exner(np,np,nlev)  ! exner nh pressure
  real (kind=real_kind), intent(out) :: pnh(np,np,nlev)   ! nh nonhyrdo pressure
  real (kind=real_kind), intent(out) :: dpnh(np,np,nlev) ! nh nonhyrdo pressure interfaces
  type (hvcoord_t),     intent(in)  :: hvcoord             ! hybrid vertical coordinate struct
  real (kind=real_kind), intent(in) :: theta(np,np,nlev)   
  real (kind=real_kind), intent(in) :: dp3d(np,np,nlev)   
  real (kind=real_kind), intent(in) :: phis(np,np)
  real (kind=real_kind), intent(in) :: phi(np,np,nlev)
  real (kind=real_kind), intent(in) :: Qdp(np,np,nlev)   

  !   local
  real (kind=real_kind) :: ptop
  real (kind=real_kind) :: kappa_star(np,np,nlev)
  real (kind=real_kind) :: kappa_star_i(np,np,nlev)
  real (kind=real_kind) :: R_star(np,np,nlev)
  real (kind=real_kind) :: Qt(np,np,nlev)
  real (kind=real_kind) :: pnh_i(np,np,nlevp) 
  real (kind=real_kind) :: rho_R_theta(np,np,nlevp) 
  real (kind=real_kind) :: rho_i(np,np,nlevp) 
  real (kind=real_kind) :: theta_i(np,np,nlevp) 
  real (kind=real_kind) :: exner_i(np,np,nlevp) 
  integer :: k


  ptop = hvcoord%hyai(1)*hvcoord%ps0

#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     if (use_moisture) then
        Qt(:,:,k) = Qdp(:,:,k)/dp3d(:,:,k)
        R_star(:,:,k)=(Rgas + (Rwater_vapor - Rgas)*Qt(:,:,k))
        
        if (use_cpstar==1) then
           kappa_star(:,:,k) = (Rgas + (Rwater_vapor - Rgas)*Qt(:,:,k)) / &
                (Cp + (Cpwater_vapor-Cp)*Qt(:,:,k) )
        else
           kappa_star(:,:,k) = (Rgas + (Rwater_vapor - Rgas)*Qt(:,:,k)) / Cp
        endif
     else
        R_star(:,:,k)=Rgas
        kappa_star(:,:,k)=Rgas/Cp
     endif
  enddo

  pnh_i(:,:,1) = ptop
  pnh_i(:,:,nlev+1) = ptop + sum(dp3d(:,:,:),3)


#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=2,nlev

     kappa_star_i(:,:,k) = (kappa_star(:,:,k)+kappa_star(:,:,k-1))/2

!     rho_i(:,:,k) = - ((dp3d(:,:,k)+dp3d(:,:,k-1))/2) / &
!          (phi(:,:,k)-phi(:,:,k-1))
!     rho_R_theta(:,:,k) = rho_i(:,:,k)*(theta(:,:,k)+theta(:,:,k-1))*&
!          (R_star(:,:,k)+R_star(:,:,k-1)) / 4
!     rho_R_theta(:,:,k) = rho_i(:,:,k)*&
!          (theta(:,:,k)*R_star(:,:,k) +theta(:,:,k-1)*R_star(:,:,k))/2
     rho_R_theta(:,:,k) = &
          (dp3d(:,:,k)  *theta(:,:,k)  *R_star(:,:,k) +&
            dp3d(:,:,k-1)*theta(:,:,k-1)*R_star(:,:,k-1))/2 / &
            (phi(:,:,k-1)-phi(:,:,k))



     if (minval(rho_R_theta(:,:,k))<0) then
        print *,k,minval( (dp3d(:,:,k)+dp3d(:,:,k-1))/2),&
             minval(phi(:,:,k-1)-phi(:,:,k))
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


#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     dpnh(:,:,k)=pnh_i(:,:,k+1)-pnh_i(:,:,k)
     pnh(:,:,k)=(pnh_i(:,:,k)+pnh_i(:,:,k+1))/2
     exner(:,:,k) = (pnh(:,:,k)/p0)**kappa_star(:,:,k)
  enddo



  end subroutine 






  subroutine copy_state(elem,nin,nout)
  implicit none
  
  type (element_t), intent(inout)   :: elem
  integer :: nin,nout

  elem%state%v(:,:,:,:,nout)=elem%state%v(:,:,:,:,nin)
  elem%state%w(:,:,:,nout)   =elem%state%w(:,:,:,nin)
  elem%state%theta(:,:,:,nout)   =elem%state%theta(:,:,:,nin)
  elem%state%phi(:,:,:,nout)   =elem%state%phi(:,:,:,nin)
  elem%state%dp3d(:,:,:,nout)=elem%state%dp3d(:,:,:,nin)
  elem%state%ps_v(:,:,nout)  =elem%state%ps_v(:,:,nin)
  end subroutine copy_state



  subroutine set_thermostate(elem,temperature,hvcoord,nt,ntQ)
!
! Assuming a hydrostatic intital state and given surface pressure,
! and no moisture, compute theta and phi 
!
! input:  ps_v, temperature
! ouput:  state variables:   theta, phi
!
  implicit none
  
  type (element_t), intent(inout)   :: elem
  real (kind=real_kind), intent(in) :: temperature(np,np,nlev)
  type (hvcoord_t),     intent(in)  :: hvcoord                      ! hybrid vertical coordinate struct
  
  !   local
  real (kind=real_kind) :: p(np,np,nlev)
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: dexner(np,np,nlev)
  real (kind=real_kind) :: integrand(np,np,nlev)
  real (kind=real_kind) :: p_i(np,np,nlev+1)
  real (kind=real_kind) :: phi_i(np,np,nlev+1)

  real (kind=real_kind) :: rho(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev)
  real (kind=real_kind) :: dphi(np,np,nlev)
  real (kind=real_kind) :: theta(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: dp_theta_R(np,np,nlev)
  integer :: k,nt,ntQ

  p_i(:,:,1) = hvcoord%hyai(1)*hvcoord%ps0
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     p(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem%state%ps_v(:,:,nt)
     dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem%state%ps_v(:,:,nt)
     p_i(:,:,k+1) = p_i(:,:,k) + dp(:,:,k)
  enddo

#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     elem%state%theta(:,:,k,nt)=temperature(:,:,k)*(p(:,:,k)/p0)**(-kappa)
  enddo


! use dry formula for exner to initialize model:
  integrand(:,:,:) = dp(:,:,:)*Rgas*temperature(:,:,:)/p(:,:,:)
  phi_i(:,:,nlev+1) = elem%state%phis(:,:)
  do k=nlev,1,-1
     phi_i(:,:,k)=phi_i(:,:,k+1) + integrand(:,:,k)
  enddo

#if 0
!  call preq_hydrostatic_v2(elem%state%phi(:,:,:,nt),elem%state%phis,integrand)
! another version (should be the same as calling preq_hydrostatic)
  do k=1,nlev
     elem%state%phi(:,:,k,nt) = (phi_i(:,:,k+1)+phi_i(:,:,k))/2
  enddo
! if we use this version, we need to define 
! dphi = phi_i(k+1)-phi_i(k)
! phi_i(k) =  2*phi(k) - phi_i(k+1)
! which only needs phi_i at the surface, not the model top
#endif

#if 0
! another version, inverse of phi_i(k) = (phi(k-1)+phi(k)/2
! with b.c. at phi_s
! use this version if we define dphi by averaging phi to interfaces and
! then differencing
  elem%state%phi(:,:,nlev,nt) = phi_i(:,:,nlev+1) + integrand(:,:,nlev)/2
  do k=nlev-1,1,-1
     elem%state%phi(:,:,k,nt) = 2*phi_i(:,:,k+1) - elem%state%phi(:,:,k+1,nt)
  enddo
#endif


#if 1
  ! inverse of get_p_nonhydrostatic
  elem%state%phi(:,:,nlev,nt) = elem%state%phis(:,:) + integrand(:,:,nlev)/2
  do k=nlev,2,-1
     !  invert this equation at interfaces:
     !  p/exner = dp_theta_R / d(phi)    
     ! d(phi) = dp_theta_R*exer/p
     dp_theta_R(:,:,k) = &
          (dp(:,:,k)  *elem%state%theta(:,:,k,nt)  *Rgas +&
           dp(:,:,k-1)*elem%state%theta(:,:,k-1,nt)*Rgas)/2
     ! (phi(:,:,k-1)-phi(:,:,k)) = dp_theta_R * exner/p
     elem%state%phi(:,:,k-1,nt) = elem%state%phi(:,:,k,nt) +&
          dp_theta_R(:,:,k) * (p_i(:,:,k)/p0)**kappa / p_i(:,:,k)
  enddo
#endif

  end subroutine set_thermostate




  subroutine set_state(u,v,T,ps,phis,p,dp,zm, g,  i,j,k,elem,n0,n1)
!
! set state variables at node(i,j,k) at layer midpoints
!
  real(real_kind),  intent(in)    :: u,v,T,ps,phis,p,dp,zm,g
  integer,          intent(in)    :: i,j,k,n0,n1
  type(element_t),  intent(inout) :: elem

  ! set prognostic state variables at level midpoints
  elem%state%v   (i,j,1,k,n0:n1) = u
  elem%state%v   (i,j,2,k,n0:n1) = v
  elem%state%ps_v(i,j,n0:n1)     = ps
  elem%state%phi(i,j,k,n0:n1)      = g*zm
  elem%state%phis(i,j)           = phis

  if (use_moisture) then
     call abortmp('ERROR: thetah set_state not yet coded for moisture')
  else
     elem%state%theta(i,j,k,n0:n1)=T/((p/p0)**kappa)
  endif
  end subroutine set_state




  
end module

