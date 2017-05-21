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
! ROUTINES REQUIRED FOR ALL MODELS:
!  get_field() 
!     returns temperature, potential temperature, phi, etc..
!  copy_state()
!     copy state variables from one timelevel to another timelevel 
!
! Test case initializaiton:
!  set_thermostate()    
!     initial condition interface used by DCMIP 2008 tests
!  set_state()
!     initial condition interface used by DCMIP 2012 tests
!  set_forcing_rayleigh_friction()
!     used by dcmip2012 test cases
!  tests_finalize()
!     used by dcmip2012 test cases
!
!
!
! OTHER ROUTINES USED BY THETA MODEL
!     perhaps these should be moved to theta_ops.F90?
!
!     get_pnh_and_exner
!     get_temperature
!     get_pottemp
!    
!
module element_ops

  use dimensions_mod, only: np, nlev, nlevp, nelemd, qsize, max_corner_elem
  use element_mod,    only: element_t
  use element_state,  only: timelevels, elem_state_t
  use hybrid_mod,     only: hybrid_t
  use hybvcoord_mod,  only: hvcoord_t
  use kinds,          only: real_kind, iulog
  use perf_mod,       only: t_startf, t_stopf, t_barrierf, t_adj_detailf ! _EXTERNAL
  use parallel_mod,   only: abortmp
  use physical_constants, only : p0, Cp, Rgas, Rwater_vapor, Cpwater_vapor, kappa, g, dd_pi
  use control_mod,    only: use_moisture, use_cpstar, theta_hydrostatic_mode
  use prim_si_mod,    only: preq_hydrostatic_v2
  implicit none

  type(elem_state_t), dimension(:), allocatable :: state0 ! storage for save_initial_state routine

contains

  recursive subroutine get_field(elem,name,field,hvcoord,nt,ntQ)
  implicit none
  type (element_t),       intent(in) :: elem
  character(len=*),       intent(in) :: name
  real (kind=real_kind),  intent(out):: field(np,np,nlev)
  type (hvcoord_t),       intent(in) :: hvcoord
  integer,                intent(in) :: nt
  integer,                intent(in) :: ntQ

  integer :: k
  real(kind=real_kind), dimension(np,np,nlev) :: tmp,p

  select case(name)
    case ('temperature'); call get_temperature(elem,field,hvcoord,nt,ntQ)
    case ('pottemp');     call get_pottemp(elem,field,hvcoord,nt,ntQ)
    case ('phi');         field = elem%state%phi(:,:,:,nt)
    case ('dpnh_dp');     call get_dpnh_dp(elem,field,hvcoord,nt,ntQ)
    case ('pnh');         call get_nonhydro_pressure(elem,field,tmp  ,hvcoord,nt,ntQ)
    case ('exner');       call get_nonhydro_pressure(elem,tmp  ,field,hvcoord,nt,ntQ)
    case ('omega');
        call get_field(elem,'p',p,hvcoord,nt,ntQ)
        field = elem%derived%omega_p * p

    case ('p');
      do k=1,nlev
        field(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem%state%ps_v(:,:,nt)
      enddo

    case ('w'); field(:,:,:)=elem%state%w(:,:,:,nt)

    case default
       print *,'name = ',trim(name)
       call abortmp('ERROR: get_field name not supported in this model')

  end select

  end subroutine


  !_____________________________________________________________________
  subroutine get_pottemp(elem,pottemp,hvcoord,nt,ntQ)
  !
  ! Should only be called outside timestep loop, state variables on reference levels
  !
  implicit none
    
  type (element_t), intent(in)        :: elem
  real (kind=real_kind), intent(out)  :: pottemp(np,np,nlev)
  type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
  integer, intent(in) :: nt
  integer, intent(in) :: ntQ
  
  !   local
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: cp_star(np,np,nlev)
  integer :: k

  do k=1,nlev
     dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem%state%ps_v(:,:,nt)
  enddo
  call get_cp_star(cp_star,elem%state%Qdp(:,:,:,1,ntQ),dp(:,:,:))
  
  pottemp(:,:,:) = elem%state%theta_dp_cp(:,:,:,nt)/(Cp_star(:,:,:)*dp(:,:,:))
  
  end subroutine get_pottemp
  

  !_____________________________________________________________________
  subroutine get_temperature(elem,temperature,hvcoord,nt,ntQ)
  !
  ! Should only be called outside timestep loop, state variables on reference levels
  !
  implicit none
  
  type (element_t), intent(in)        :: elem
  real (kind=real_kind), intent(out)  :: temperature(np,np,nlev)
  type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
  integer, intent(in) :: nt
  integer, intent(in) :: ntQ
  
  !   local
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: cp_star(np,np,nlev)
  real (kind=real_kind) :: kappa_star(np,np,nlev)
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
  call get_cp_star(cp_star,elem%state%Qdp(:,:,:,1,ntQ),dp)
  call get_kappa_star(kappa_star,elem%state%Qdp(:,:,:,1,ntQ),dp)


  call get_pnh_and_exner(hvcoord,elem%state%theta_dp_cp(:,:,:,nt),&
          dp,elem%state%phi(:,:,:,nt),elem%state%phis(:,:),kappa_star,&
          pnh,dpnh,exner)

  
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     temperature(:,:,k)= elem%state%theta_dp_cp(:,:,k,nt)*exner(:,:,k)/(Cp_star(:,:,k)*dp(:,:,k))
  enddo

  end subroutine get_temperature


  !_____________________________________________________________________
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
  real (kind=real_kind) :: kappa_star(np,np,nlev)
  integer :: k
  
  
  do k=1,nlev
     dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem%state%ps_v(:,:,nt)
  enddo
  call get_kappa_star(kappa_star,elem%state%Qdp(:,:,:,1,ntQ),dp)

  call get_pnh_and_exner(hvcoord,elem%state%theta_dp_cp(:,:,:,nt),&
       dp,elem%state%phi(:,:,:,nt),elem%state%phis(:,:),kappa_star,&
       pnh,dpnh,exner)
  dpnh_dp = dpnh/dp
  end subroutine 

  !_____________________________________________________________________
  subroutine get_nonhydro_pressure(elem,pnh,exner,hvcoord,nt,ntQ)
    implicit none
    
    type (element_t),       intent(in)  :: elem
    real (kind=real_kind),  intent(out) :: pnh(np,np,nlev)
    real (kind=real_kind),  intent(out) :: exner(np,np,nlev)
    type (hvcoord_t),       intent(in)  :: hvcoord
    integer,                intent(in)  :: nt
    integer,                intent(in)  :: ntQ
    
    real (kind=real_kind), dimension(np,np,nlev) :: dp,dpnh,kappa_star
    integer :: k

    do k=1,nlev
      dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
      (hvcoord%hybi(k+1)-hvcoord%hybi(k))*elem%state%ps_v(:,:,nt)
    enddo

    call get_kappa_star(kappa_star,elem%state%Qdp(:,:,:,1,ntQ),dp)

    call get_pnh_and_exner(hvcoord,elem%state%theta_dp_cp(:,:,:,nt),&
         dp,elem%state%phi(:,:,:,nt),elem%state%phis(:,:),kappa_star,&
         pnh,dpnh,exner)

  end subroutine


  !_____________________________________________________________________
  subroutine get_kappa_star(kappa_star,Qdp,dp)
  !
  ! note: interface written in this way so that it can be called outside timelevel loop,
  ! where dp is computed from reverence levels, or inside timelevel loop where dp = prognostic dp3d
  !
  implicit none
  real (kind=real_kind), intent(out)  :: kappa_star(np,np,nlev)
  real (kind=real_kind), intent(in)   :: Qdp(np,np,nlev)
  real (kind=real_kind), intent(in)   :: dp(np,np,nlev)
  !   local
  integer :: k

  if (use_moisture .and. use_cpstar==1) then
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
     do k=1,nlev
        kappa_star(:,:,k) = (Rgas + (Rwater_vapor - Rgas)*Qdp(:,:,k)/dp(:,:,k)) / &
             (Cp + (Cpwater_vapor-Cp)*Qdp(:,:,k)/dp(:,:,k) )
     enddo
  else if (use_moisture .and. use_cpstar==0) then
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
     do k=1,nlev
        kappa_star(:,:,k) = (Rgas + (Rwater_vapor - Rgas)*Qdp(:,:,k)/dp(:,:,k)) / Cp
     enddo
  else
     kappa_star(:,:,:)=Rgas/Cp
  endif
  end subroutine 

  !_____________________________________________________________________
  subroutine get_cp_star(cp_star,Qdp,dp)
  !
  ! note: interface written in this way so that it can be called outside timelevel loop,
  ! where dp is computed from reverence levels, or inside timelevel loop where dp = prognostic dp3d
  !
  implicit none
  real (kind=real_kind), intent(out):: cp_star(np,np,nlev)
  real (kind=real_kind), intent(in) :: Qdp(np,np,nlev)
  real (kind=real_kind), intent(in) :: dp(np,np,nlev)

  integer :: k
  if (use_moisture .and. use_cpstar==1) then
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
     do k=1,nlev
        cp_star(:,:,k) = (Cp + (Cpwater_vapor-Cp)*Qdp(:,:,k)/dp(:,:,k) )
     enddo
  else
     cp_star(:,:,:)=Cp
  endif
  end subroutine






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

! new EOS
  real (kind=real_kind) :: exner_over_p(np,np,nlev)
  real (kind=real_kind) :: exner_over_p_i(np,np,nlevp)
  real (kind=real_kind) :: dphi(np,np,nlev)
  real (kind=real_kind) :: dphi_i(np,np,nlevp) 

! old OES
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
  subroutine copy_state(elem,nin,nout)
  implicit none
  
  type (element_t), intent(inout)   :: elem
  integer :: nin,nout

  elem%state%v(:,:,:,:,nout)  =elem%state%v(:,:,:,:,nin)
  elem%state%w(:,:,:,nout)    =elem%state%w(:,:,:,nin)
  elem%state%theta_dp_cp(:,:,:,nout) =elem%state%theta_dp_cp(:,:,:,nin)
  elem%state%phi(:,:,:,nout)  =elem%state%phi(:,:,:,nin)
  elem%state%dp3d(:,:,:,nout) =elem%state%dp3d(:,:,:,nin)
  elem%state%ps_v(:,:,nout)   =elem%state%ps_v(:,:,nin)
  end subroutine copy_state


  !_____________________________________________________________________
  subroutine set_thermostate(elem,temperature,hvcoord,nt,ntQ)
  !
  ! Assuming a hydrostatic intital state and given surface pressure,
  ! and no moisture, compute theta and phi 
  !
  ! input:  ps_v, temperature
  ! ouput:  state variables:   theta_dp_cp, phi
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
  real (kind=real_kind) :: kappa_star(np,np,nlev)
  integer :: k,nt,ntQ

  do k=1,nlev
     p(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem%state%ps_v(:,:,nt)
     dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem%state%ps_v(:,:,nt)
  enddo

  do k=1,nlev
     elem%state%theta_dp_cp(:,:,k,nt)=temperature(:,:,k)*(p(:,:,k)/p0)**(-kappa)*      &
      (Cp*dp(:,:,k))
  enddo

  call set_hydrostatic_phi(hvcoord,elem%state%phis,elem%state%theta_dp_cp(:,:,:,nt),dp,&
       elem%state%phi(:,:,:,nt))


  ! debug
  kappa_star=kappa   
  call get_pnh_and_exner(hvcoord,elem%state%theta_dp_cp(:,:,:,nt),dp,&
       elem%state%phi(:,:,:,nt),&
       elem%state%phis(:,:),kappa_star,pnh,dpnh,exner)
  do k=1,nlev
     if (maxval(abs(1-dpnh(:,:,k)/dp(:,:,k))) > 1e-10) then
        print *,'WARNING: hydrostatic inverse FAILED!'
        print *, minval(dpnh(:,:,k)), minval(dp(:,:,k))
        print *,k,minval(dpnh(:,:,k)/dp(:,:,k)),maxval(dpnh(:,:,k)/dp(:,:,k))
     endif
  enddo

  end subroutine set_thermostate

  !_____________________________________________________________________
  subroutine set_hydrostatic_phi(hvcoord,phis,theta_dp_cp,dp,phi)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! inverse of get_pnh_and_exner() for dry hydrostatic case
  ! used to initialize phi for dry test cases
  ! used to compute background phi for reference state
  ! in both the above uses, we can assume dry and we dont need to compute kappa_star
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
  subroutine set_moist_hydrostatic_phi(hvcoord,phis,theta_dp_cp,dp,kappa_star,phi)

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

  !_____________________________________________________________________
  subroutine set_theta_ref(hvcoord,dp,theta_ref)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! inverse of get_pnh_and_exner() for dry hydrostatic case
  ! used to initialize phi for dry test cases
  ! used to compute background phi for reference state
  ! in both the above uses, we can assume dry and we dont need to compute kappa_star
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none
  
  type (hvcoord_t),     intent(in)  :: hvcoord                      ! hybrid vertical coordinate struct
  real (kind=real_kind), intent(in) :: dp(np,np,nlev)
  real (kind=real_kind), intent(out) :: theta_ref(np,np,nlev)
  
  !   local
  real (kind=real_kind) :: p_i(np,np,nlevp)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: T0,T1
  integer :: k

  ! reference T = 288K.  reference lapse rate = 6.5K/km   = .0065 K/m
  ! Tref = T0+T1*exner
  ! Thetaref = T0/exner + T1
  T1 = .0065*288d0*Cp/g ! = 191
  T0 = 288d0-T1         ! = 97

  p_i(:,:,1) =  hvcoord%hyai(1)*hvcoord%ps0   
  do k=1,nlev
     p_i(:,:,k+1) = p_i(:,:,k) + dp(:,:,k)
  enddo
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     exner(:,:,k) = ( (p_i(:,:,k) + p_i(:,:,k+1))/(2*p0)) **kappa
     !theta_ref(:,:,k,ie) = (T0/exner(:,:,k) + T1)*Cp*dp_ref(:,:,k,ie)
     theta_ref(:,:,k) = (T0/exner(:,:,k) + T1)
  enddo

  end subroutine


  !_____________________________________________________________________
  subroutine set_state(u,v,w,T,ps,phis,p,dp,zm,g,i,j,k,elem,n0,n1)
  !
  ! set state variables at node(i,j,k) at layer midpoints
  ! used by idealized tests for dry initial conditions
  ! so we use constants cp, kappa  
  !
  real(real_kind),  intent(in)    :: u,v,w,T,ps,phis,p,dp,zm,g
  integer,          intent(in)    :: i,j,k,n0,n1
  type(element_t),  intent(inout) :: elem

  ! set prognostic state variables at level midpoints
  elem%state%v   (i,j,1,k,n0:n1) = u
  elem%state%v   (i,j,2,k,n0:n1) = v
  elem%state%w   (i,j,k,  n0:n1) = w
  elem%state%dp3d(i,j,k,  n0:n1) = dp
  elem%state%ps_v(i,j,    n0:n1) = ps
  elem%state%phi (i,j,k,  n0:n1) = g*zm
  elem%state%phis(i,j)           = phis
  elem%state%theta_dp_cp(i,j,k,n0:n1)=T*Cp*dp*((p/p0)**(-kappa))

  end subroutine set_state

  !_____________________________________________________________________
  subroutine set_elem_state(u,v,w,T,ps,phis,p,dp,zm,g,elem,n0,n1,ntQ)

  ! set element state variables
  ! works for both dry and moist initial conditions

  real(real_kind), dimension(np,np,nlev), intent(in):: u,v,w,T,p,dp,zm
  real(real_kind), dimension(np,np),      intent(in):: ps,phis
  real(real_kind),  intent(in)    :: g
  integer,          intent(in)    :: n0,n1,ntQ
  type(element_t),  intent(inout) :: elem
  integer :: n
  real(real_kind), dimension(np,np,nlev) :: cp_star,kappa_star

  ! get cp and kappa for dry or moist cases
  call get_cp_star(cp_star,elem%state%Qdp(:,:,:,1,ntQ),dp)
  call get_kappa_star(kappa_star,elem%state%Qdp(:,:,:,1,ntQ),dp)

  do n=n0,n1
    ! set prognostic state variables at level midpoints
    elem%state%v   (:,:,1,:,n) = u
    elem%state%v   (:,:,2,:,n) = v
    elem%state%w   (:,:,:,  n) = w
    elem%state%dp3d(:,:,:,  n) = dp
    elem%state%ps_v(:,:,    n) = ps
    elem%state%phi (:,:,:,  n) = g*zm
    elem%state%phis(:,:)       = phis
    elem%state%theta_dp_cp(:,:,:,n)=T*cp_star*dp*((p/p0)**(-kappa_star))
  end do

  end subroutine set_elem_state

  !_____________________________________________________________________
  subroutine get_state(u,v,w,T,theta,exner,pnh,dp,Cp_star,rho,zm,g,i,j,elem,hvcoord,nt,ntQ)

    ! get state variables at layer midpoints
    ! used by idealized tests to compute idealized physics forcing terms

    real(real_kind), dimension(np,np,nlev), intent(inout) :: u,v,w,T,theta,exner,pnh,dp,cp_star,zm,rho
    real(real_kind), intent(in)    :: g
    integer,         intent(in)    :: i,j,nt,ntQ
    type(element_t), intent(inout) :: elem
    type (hvcoord_t),intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct

    real(real_kind) , dimension(np,np,nlev) :: phi,dpnh,kappa_star,Rstar
    real(real_kind) , dimension(np,np) :: phis

    integer :: k

    ! set prognostic state variables at level midpoints
    u   = elem%state%v   (:,:,1,:,nt)
    v   = elem%state%v   (:,:,2,:,nt)
    w   = elem%state%w   (:,:,  :,nt)
    phi = elem%state%phi (:,:,  :,nt)
    phis= elem%state%phis(:,:)
    zm  = phi/g

   ! dp  = elem%state%dp3d(:,:,  :,nt)

    do k=1,nlev
       dp(:,:,k)=(hvcoord%hyai(k+1)-hvcoord%hyai(k))*hvcoord%ps0 +(hvcoord%hybi(k+1)-hvcoord%hybi(k))*elem%state%ps_v(:,:,nt)
    enddo

    call get_cp_star(cp_star,elem%state%Qdp(:,:,:,1,ntQ),dp)
    call get_kappa_star(kappa_star,elem%state%Qdp(:,:,:,1,ntQ),dp)
    call get_pnh_and_exner(hvcoord,elem%state%theta_dp_cp(:,:,:,nt),dp,phi,phis,kappa_star,pnh,dpnh,exner)

    Rstar = kappa_star*cp_star
    theta = elem%state%theta_dp_cp(:,:,:,nt)/(Cp_star*dp)
    T     = theta*exner
    rho   = pnh/(Rstar*T)

  end subroutine get_state

  !_____________________________________________________________________
  subroutine save_initial_state(state,ie)

    ! save the state at time=0 for use with some dcmip tests

    type(elem_state_t), intent(inout):: state
    integer,            intent(in)   :: ie     ! element index

    if(.not. allocated(state0)) allocate( state0(nelemd) )
    state0(ie) = state

  end subroutine


  !_____________________________________________________________________
  subroutine set_forcing_rayleigh_friction(elem, zm, ztop, zc, tau, u0,v0, n)
  !
  ! test cases which use rayleigh friciton will call this with the relaxation coefficient
  ! f_d, and the reference state u0,v0.  Currently assume w0 = 0
  !
  implicit none

  type(element_t), intent(inout):: elem
  real(real_kind), intent(in)   :: zm(np,np,nlev) ! height at layer midpoints
  real(real_kind), intent(in)   :: ztop           ! top of atm height
  real(real_kind), intent(in)   :: zc             ! cutoff height
  real(real_kind), intent(in)   :: tau            ! damping timescale
  real(real_kind), intent(in)   :: u0(np,np,nlev) ! reference u
  real(real_kind), intent(in)   :: v0(np,np,nlev) ! reference v
  integer,         intent(in)   :: n              ! timestep

  real(real_kind):: f_d(np,np,nlev)
  integer :: k

  ! Compute damping as a function of layer-midpoint height
  f_d=0.0d0
  where(zm .ge. zc);   f_d = sin(dd_pi/2 *(zm - zc)/(ztop - zc))**2; end where
  where(zm .ge. ztop); f_d = 1.0d0; end where
  f_d = -f_d/tau

  elem%derived%FM(:,:,1,:) = f_d * ( elem%state%v(:,:,1,:,n) - u0 )
  elem%derived%FM(:,:,2,:) = f_d * ( elem%state%v(:,:,2,:,n) - v0)
  elem%derived%FM(:,:,3,:) = f_d * ( elem%state%w(:,:,:,n)  )

  end subroutine 

  !_____________________________________________________________________
  subroutine tests_finalize(elem, hvcoord,ns,ne,ie)

  ! Now that all variables have been initialized, set phi to be in hydrostatic balance

  implicit none

  type(hvcoord_t),     intent(in)   :: hvcoord
  type(element_t),     intent(inout):: elem
  integer,             intent(in)   :: ns,ne
  integer, optional,   intent(in)   :: ie ! optional element index, to save initial state

  integer :: k,tl, ntQ
  real(real_kind), dimension(np,np,nlev) :: dp, kappa_star

  do k=1,nlev
    dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k))*elem%state%ps_v(:,:,ns)
  enddo

  ntQ=1
  call get_kappa_star(kappa_star,elem%state%Qdp(:,:,:,1,ntQ),dp)
  call set_moist_hydrostatic_phi(hvcoord,elem%state%phis,elem%state%theta_dp_cp(:,:,:,ns),dp,kappa_star,elem%state%phi(:,:,:,ns))

  do tl = ns+1,ne
    call copy_state(elem,ns,tl)
  enddo

  if(present(ie)) call save_initial_state(elem%state,ie)

  end subroutine tests_finalize

end module

