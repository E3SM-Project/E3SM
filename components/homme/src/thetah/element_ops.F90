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
  case ('omega')
!     do k=1,nlev
!        field(:,:,k)=elem%derived%omega_p(:,:,k)*&
!             (hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem%state%ps_v(:,:,nt))
!     end do
     field(:,:,:)=elem%state%w(:,:,:,nt)
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

  call get_pnh_and_exner(hvcoord,elem%state%theta(:,:,:,nt),&
          dp,elem%state%phi(:,:,:,nt),elem%state%phis(:,:),elem%state%Qdp(:,:,:,1,ntQ),&
          pnh,dpnh,exner)


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


  call get_pnh_and_exner(hvcoord,elem%state%theta(:,:,:,nt),&
       dp,elem%state%phi(:,:,:,nt),elem%state%phis(:,:),elem%state%Qdp(:,:,:,1,ntQ),&
       pnh,dpnh,exner)
  dpnh_dp = dpnh/dp
  end subroutine 






  subroutine get_pnh_and_exner(hvcoord,theta,dp3d,phi,phis,Qdp,pnh,dpnh,exner,exner_i_out)
  implicit none
!
! compute exner pressure, nh presure
!
! input:  dp3d, Qdp (if use_moisture), phi, phis, theta
! output:  pnh, dphn, exner, exner_i
!       for k=2..nlev, use the equation of state:  pnh/e = rho*Rstar*theta
!                                                    rho   = -dp3d/dphi 
!
!  
  type (hvcoord_t),     intent(in)  :: hvcoord             ! hybrid vertical coordinate struct
  real (kind=real_kind), intent(in) :: theta(np,np,nlev)   
  real (kind=real_kind), intent(in) :: dp3d(np,np,nlev)   
  real (kind=real_kind), intent(in) :: phi(np,np,nlev)
  real (kind=real_kind), intent(in) :: phis(np,np)
  real (kind=real_kind), intent(in) :: Qdp(np,np,nlev)   
  real (kind=real_kind), intent(out) :: pnh(np,np,nlev)   ! nh nonhyrdo pressure
  real (kind=real_kind), intent(out) :: dpnh(np,np,nlev) ! nh nonhyrdo pressure interfaces
  real (kind=real_kind), intent(out) :: exner(np,np,nlev)  ! exner nh pressure
  real (kind=real_kind), intent(out), optional :: exner_i_out(np,np,nlevp)  ! exner nh pressure interfaces

  !   local
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
     
  else
!==============================================================
!  non-hydrostatic formulas
!==============================================================
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! boundary terms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  compute pnh at midpoint, then extrapolate to surface:
!   rho_R_theta(:,:,nlev) = dp3d(:,:,nlev)*theta(:,:,nlev)*R_star(:,:,nlev) / &
!        (  (phi(:,:,nlev)+phi(:,:,nlev-1))/2 - phis(:,:)  )
   rho_R_theta(:,:,nlev) = dp3d(:,:,nlev)*theta(:,:,nlev)*R_star(:,:,nlev) / &
        (  (phi(:,:,nlev) - phis(:,:))*2  )

   exner(:,:,nlev) = (rho_R_theta(:,:,nlev)/p0)**&
        ( kappa_star(:,:,nlev)/ ( 1-kappa_star(:,:,nlev)))
   pnh(:,:,nlev) = rho_R_theta(:,:,nlev)*exner(:,:,nlev)
!  invert  pnh(:,:,k)=(pnh_i(:,:,k)+pnh_i(:,:,k+1))/2
   pnh_i(:,:,nlev+1) = 2*pnh(:,:,nlev) - pnh_i(:,:,nlev)

   pnh_i(:,:,1) = hvcoord%hyai(1)*hvcoord%ps0   ! hydrostatic ptop
!   pnh_i(:,:,1)      = pnh_i(:,:,2)    - dp3d(:,:,1) 






#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     dpnh(:,:,k)=pnh_i(:,:,k+1)-pnh_i(:,:,k)
     pnh(:,:,k)=(pnh_i(:,:,k)+pnh_i(:,:,k+1))/2
     exner(:,:,k) = (pnh(:,:,k)/p0)**kappa_star(:,:,k)
  enddo

  if (present(exner_i_out)) then
     exner_i_out = exner_i  ! computed above for k=2,nlev
     ! how to approximate kappa_star at top and bottom of model?  
     exner_i_out(:,:,1) = (pnh_i(:,:,1)/p0)**kappa_star(:,:,1)
     exner_i_out(:,:,nlev+1) = (pnh_i(:,:,nlev+1)/p0)**kappa_star(:,:,nlev)
  endif


  endif ! hydrostatic/nonhydrostatic version
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
  real (kind=real_kind) :: pnh(np,np,nlev)
  real (kind=real_kind) :: dpnh(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)
  integer :: k,nt,ntQ

  do k=1,nlev
     p(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem%state%ps_v(:,:,nt)
     dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem%state%ps_v(:,:,nt)
  enddo

  do k=1,nlev
     elem%state%theta(:,:,k,nt)=temperature(:,:,k)*(p(:,:,k)/p0)**(-kappa)
  enddo

  call set_hydrostatic_phi(hvcoord,elem%state%phis,elem%state%theta(:,:,:,nt),dp,&
       elem%state%phi(:,:,:,nt))


  ! debug
  call get_pnh_and_exner(hvcoord,elem%state%theta(:,:,:,nt),dp,&
       elem%state%phi(:,:,:,nt),&
       elem%state%phis(:,:),elem%state%Qdp(:,:,:,1,nt),pnh,dpnh,exner)
  do k=1,nlev
     if (maxval(abs(1-dpnh(:,:,k)/dp(:,:,k))) > 1e-10) then
        print *,'WARNING: hydrostatic inverse FAILED!'
        print *,k,minval(dpnh(:,:,k)/dp(:,:,k)),maxval(dpnh(:,:,k)/dp(:,:,k))
     endif
  enddo

  end subroutine set_thermostate






  subroutine set_hydrostatic_phi(hvcoord,phis,theta,dp,phi)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! inverse of get_pnh_and_exner():
  ! for dry hydrostatic case
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  
  type (hvcoord_t),     intent(in)  :: hvcoord                      ! hybrid vertical coordinate struct
  real (kind=real_kind), intent(in) :: theta(np,np,nlev)
  real (kind=real_kind), intent(in) :: dp(np,np,nlev)
  real (kind=real_kind), intent(in) :: phis(np,np)
  real (kind=real_kind), intent(out) :: phi(np,np,nlev)
  
  !   local
  real (kind=real_kind) :: p(np,np,nlev)
  real (kind=real_kind) :: p_i(np,np,nlev+1)

  real (kind=real_kind) :: dp_theta_R(np,np,nlev)
  integer :: k

  p_i(:,:,1) =  hvcoord%hyai(1)*hvcoord%ps0   
  do k=1,nlev
     p_i(:,:,k+1) = p_i(:,:,k) + dp(:,:,k)
  enddo
  do k=1,nlev
     p(:,:,k) = (p_i(:,:,k) + p_i(:,:,k+1))/2
  enddo

!  integrand(:,:) = dp(:,:,nlev)*Rgas*temperature(:,:,nlev)/p(:,:,nlev)
  phi(:,:,nlev) = phis(:,:) + (&
    dp(:,:,nlev)*Rgas*theta(:,:,nlev)*p(:,:,nlev)**(kappa-1)*p0**(-kappa) )/2

  do k=nlev,2,-1
     !  invert this equation at interfaces:
     !  p/exner = dp_theta_R / d(phi)    
     ! d(phi) = dp_theta_R*exer/p
     dp_theta_R(:,:,k) = &
          (dp(:,:,k)  *theta(:,:,k)  *Rgas +&
           dp(:,:,k-1)*theta(:,:,k-1)*Rgas)/2
     ! (phi(:,:,k-1)-phi(:,:,k)) = dp_theta_R * exner/p
     phi(:,:,k-1) = phi(:,:,k) +&
          dp_theta_R(:,:,k) * (p_i(:,:,k)/p0)**kappa / p_i(:,:,k)
  enddo


  end subroutine



 
  subroutine set_state(u,v,w,T,ps,phis,p,dp,zm, g,  i,j,k,elem,n0,n1)
!
! set state variables at node(i,j,k) at layer midpoints
!
  real(real_kind),  intent(in)    :: u,v,w,T,ps,phis,p,dp,zm,g
  integer,          intent(in)    :: i,j,k,n0,n1
  type(element_t),  intent(inout) :: elem

  ! set prognostic state variables at level midpoints
  elem%state%v   (i,j,1,k,n0:n1) = u
  elem%state%v   (i,j,2,k,n0:n1) = v
  elem%state%w   (i,j,k,n0:n1)   = w
  elem%state%ps_v(i,j,n0:n1)     = ps
  elem%state%phi(i,j,k,n0:n1)      = g*zm
  elem%state%phis(i,j)           = phis
  elem%state%theta(i,j,k,n0:n1)=T/((p/p0)**kappa)

  end subroutine set_state



  subroutine set_forcing_rayleigh_friction(elem, f_d, u0,v0, n)
!
! test cases which use rayleigh friciton will call this with the relaxation coefficient
! f_d, and the reference state u0,v0.  Currently assume w0 = 0
!
  implicit none

  type(element_t),  intent(inout)  :: elem
  real(real_kind):: u0(np,np,nlev)
  real(real_kind):: v0(np,np,nlev)
  real(real_kind):: f_d(nlev)
  integer :: n,k

  do k=1,nlev
     elem%derived%FM(:,:,1,k) = f_d(k) * ( elem%state%v(:,:,1,k,n) - u0(:,:,k) )
     elem%derived%FM(:,:,2,k) = f_d(k) * ( elem%state%v(:,:,2,k,n) - v0(:,:,k) )
     elem%derived%FM(:,:,3,k) = f_d(k) * ( elem%state%w(:,:,k,n)  )
  enddo
  end subroutine 





  subroutine tests_finalize(elem, hvcoord,ns,ne)
!
! Now that all variables have been initialized, set phi to be in hydrostatic balance
!
  implicit none

  type(hvcoord_t),     intent(in)  :: hvcoord
  type(element_t),  intent(inout)  :: elem
  integer,             intent(in)  :: ns,ne

  integer :: k,ie,tl
  real(real_kind):: dp(np,np,nlev)

!  set dp
  do k=1,nlev
    dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k))*elem%state%ps_v(:,:,ns)
  enddo
  call set_hydrostatic_phi(hvcoord,elem%state%phis,elem%state%theta(:,:,:,ns),dp,elem%state%phi(:,:,:,ns))

  do tl = ns+1,ne
    call copy_state(elem,ns,tl)
  enddo

  end subroutine tests_finalize

end module

