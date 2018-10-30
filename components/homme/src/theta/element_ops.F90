#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
!
!  getter and setter functions that must be provided by each model
!  
!  !!IMPORTANT NOTE!!:  These routines assume we are on REFERENCE levels
!  For vertically lagrangian models, they should only be used outside the dynamics
!  timestep after the vertical remap.  

! ROUTINES REQUIRED FOR ALL MODELS:
!  get_field() 
!     returns temperature, potential temperature, phi, etc..
!  copy_state()
!     copy state variables from one timelevel to another timelevel 
!  set_thermostate()    
!     initial condition interface used by DCMIP 2008 tests, old HOMME tests
!  set_state()
!     initial condition interface used by DCMIP 2012 tests
!  set_elem_state()
!     initial condition interface used by DCMIP 2016 tests
!  get_state()
!     return state variables used by some DCMIP forcing functions
!  save_initial_state()
!     save t=0 in "state0", used by some DCMIP forcing functions       
!  set_forcing_rayleigh_friction()
!     used by dcmip2012 test cases
!  tests_finalize()
!     initialize geopotential to be in hydrostatic balance
!     used by DCMIP2012, 2016 tests
!
!  Accessory routines used by the above:
!  get_pottemp()
!  get_temperature()
!  get_dpnh_dp()
!  get_nonhydro_pressure()
!
!
!  Additional routines that work for both reference levels and Lagrangian levels:
!  (because they accept dp as an input argument)
!  get_kappa_star()
!  get_cp_star()
!  set_theta_ref()
!  
!
module element_ops

  use dimensions_mod, only: np, nlev, nlevp, nelemd
  use element_mod,    only: element_t
  use element_state,  only: elem_state_t
  use hybvcoord_mod,  only: hvcoord_t
  use kinds,          only: real_kind, iulog
  use perf_mod,       only: t_startf, t_stopf, t_barrierf, t_adj_detailf ! _EXTERNAL
  use parallel_mod,   only: abortmp
  use physical_constants, only : p0, Cp, Rgas, Rwater_vapor, Cpwater_vapor, kappa, g, dd_pi
  use control_mod,    only: use_moisture, use_cpstar, theta_hydrostatic_mode
  use eos,            only: get_pnh_and_exner, get_moist_phinh, get_dry_phinh
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
  real(kind=real_kind), dimension(np,np,nlev) :: tmp, p, pnh, dp, omega, rho, T, cp_star, Rstar, kappa_star

  select case(name)
    case ('temperature','T'); call get_temperature(elem,field,hvcoord,nt)
    case ('pottemp','Th');    call get_pottemp(elem,field,hvcoord,nt,ntQ)
    case ('phi','geo');       call get_phi(elem,field,hvcoord,nt,ntQ)
    case ('dpnh_dp');         call get_dpnh_dp(elem,field,hvcoord,nt,ntQ)
    case ('pnh');             call get_nonhydro_pressure(elem,field,tmp  ,hvcoord,nt,ntQ)
    case ('exner');           call get_nonhydro_pressure(elem,tmp  ,field,hvcoord,nt,ntQ)

    case ('p');
      do k=1,nlev
        field(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem%state%ps_v(:,:,nt)
      enddo

    case ('dp');
      do k=1,nlev
        field(:,:,k)=(hvcoord%hyai(k+1)-hvcoord%hyai(k))*hvcoord%ps0 +(hvcoord%hybi(k+1)-hvcoord%hybi(k))*elem%state%ps_v(:,:,nt)
      enddo

    case ('omega');
      call get_field(elem,'p',p,hvcoord,nt,ntQ)
      field = elem%derived%omega_p*p

    case('rho')

      call get_field(elem,'pnh',pnh,hvcoord,nt,ntQ)
      call get_field(elem,'dp',dp,hvcoord,nt,ntQ)
      call get_cp_star(cp_star,elem%state%Q(:,:,:,1))
      call get_kappa_star(kappa_star,elem%state%Q(:,:,:,1))
      call get_temperature(elem,T,hvcoord,nt)

      Rstar = cp_star*kappa_star
      field = pnh/(Rstar*T)

    case ('w');

      if(theta_hydrostatic_mode) then
        call get_field(elem,'omega',omega,hvcoord,nt,ntQ)
        call get_field(elem,'rho'  ,rho  ,hvcoord,nt,ntQ)
        field = -omega/(rho *g)

      else
        field =elem%state%w(:,:,:,nt)
      endif

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
  call get_cp_star(cp_star,elem%state%Q(:,:,:,1))
  
  pottemp(:,:,:) = elem%state%theta_dp_cp(:,:,:,nt)/(Cp_star(:,:,:)*dp(:,:,:))
  
  end subroutine get_pottemp
  

  !_____________________________________________________________________
  subroutine get_temperature(elem,temperature,hvcoord,nt)
  !
  ! Should only be called outside timestep loop, state variables on reference levels
  !
  implicit none
  
  type (element_t), intent(in)        :: elem
  real (kind=real_kind), intent(out)  :: temperature(np,np,nlev)
  type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
  integer, intent(in) :: nt
  
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
  call get_cp_star(cp_star,elem%state%Q(:,:,:,1))
  call get_kappa_star(kappa_star,elem%state%Q(:,:,:,1))


  call get_pnh_and_exner(hvcoord,elem%state%theta_dp_cp(:,:,:,nt),&
          dp,elem%state%phinh(:,:,:,nt),elem%state%phis(:,:),kappa_star,&
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
  call get_kappa_star(kappa_star,elem%state%Q(:,:,:,1))

  call get_pnh_and_exner(hvcoord,elem%state%theta_dp_cp(:,:,:,nt),&
       dp,elem%state%phinh(:,:,:,nt),elem%state%phis(:,:),kappa_star,&
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

    call get_kappa_star(kappa_star,elem%state%Q(:,:,:,1))

    call get_pnh_and_exner(hvcoord,elem%state%theta_dp_cp(:,:,:,nt),&
         dp,elem%state%phinh(:,:,:,nt),elem%state%phis(:,:),kappa_star,&
         pnh,dpnh,exner)

  end subroutine



  subroutine get_phi(elem,phi,hvcoord,nt,ntQ)
    implicit none
    
    type (element_t),       intent(in)  :: elem
    type (hvcoord_t),       intent(in)  :: hvcoord
    real (kind=real_kind),  intent(out) :: phi(np,np,nlev)
    integer,                intent(in)  :: nt
    integer,                intent(in)  :: ntQ
    
    real (kind=real_kind), dimension(np,np,nlev) :: dp,dpnh,kappa_star
    real (kind=real_kind) :: pnh(np,np,nlev)
    real (kind=real_kind) :: exner(np,np,nlev)
    real (kind=real_kind) :: temp(np,np,nlev)
    integer :: k

    if(theta_hydrostatic_mode) then
       do k=1,nlev
          dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               (hvcoord%hybi(k+1)-hvcoord%hybi(k))*elem%state%ps_v(:,:,nt)
       enddo
       
       call get_kappa_star(kappa_star,elem%state%Q(:,:,:,1))
       
       call get_pnh_and_exner(hvcoord,elem%state%theta_dp_cp(:,:,:,nt),&
            dp,elem%state%phinh(:,:,:,nt),elem%state%phis(:,:),kappa_star,&
            pnh,dpnh,exner)

       do k=1,nlev
          !temp(:,:,k) = theta_dp_cp(:,:,k)*(exner_i(:,:,k+1)-exner_i(:,:,k))/dp3d(:,:,k)           
          temp(:,:,k) = kappa_star(:,:,k)*elem%state%theta_dp_cp(:,:,k,nt)*exner(:,:,k)/pnh(:,:,k)
       enddo
       call preq_hydrostatic_v2(phi,elem%state%phis,temp)
       
    else
       phi = elem%state%phinh(:,:,:,nt)
    endif
    
  end subroutine

       

  subroutine get_phi_i(elem,phi_i,hvcoord,nt,ntQ)
    implicit none
    
    type (element_t),       intent(in)  :: elem
    type (hvcoord_t),       intent(in)  :: hvcoord
    real (kind=real_kind),  intent(out) :: phi_i(np,np,nlevp)
    integer,                intent(in)  :: nt
    integer,                intent(in)  :: ntQ
    
    real (kind=real_kind), dimension(np,np,nlev) :: dp,dpnh,kappa_star
    real (kind=real_kind) :: phi(np,np,nlev)
    real (kind=real_kind) :: pnh(np,np,nlev)
    real (kind=real_kind) :: exner(np,np,nlev)
    real (kind=real_kind) :: temp(np,np,nlev)
    integer :: k

    ! compute hydrostatic version first:

    do k=1,nlev
       dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
            (hvcoord%hybi(k+1)-hvcoord%hybi(k))*elem%state%ps_v(:,:,nt)
    enddo
    
    call get_kappa_star(kappa_star,elem%state%Q(:,:,:,1))
    
    call get_pnh_and_exner(hvcoord,elem%state%theta_dp_cp(:,:,:,nt),&
         dp,elem%state%phinh(:,:,:,nt),elem%state%phis(:,:),kappa_star,&
         pnh,dpnh,exner)
    
    do k=1,nlev
       !temp(:,:,k) = theta_dp_cp(:,:,k)*(exner_i(:,:,k+1)-exner_i(:,:,k))/dp3d(:,:,k)           
       temp(:,:,k) = kappa_star(:,:,k)*elem%state%theta_dp_cp(:,:,k,nt)*exner(:,:,k)/pnh(:,:,k)
    enddo
    
    phi_i(:,:,nlevp) = elem%state%phis(:,:)
    ! traditional Hydrostatic integral
    do k=nlev,1,-1
       phi_i(:,:,k)=phi_i(:,:,k+1)+temp(:,:,k)
    enddo
    
    if (.not. theta_hydrostatic_mode) then
       phi = elem%state%phinh(:,:,:,nt)
       ! average phinh 
       phi_i(:,:,2:nlev)= (phi(:,:,2:nlev)+phi(:,:,1:nlev-1))/2

       ! boundaries:
       phi_i(:,:,nlevp) = elem%state%phis(:,:)
       phi_i(:,:,1)     = phi(:,:,1)+temp(:,:,1)/2
    endif
    
  end subroutine

       





  !_____________________________________________________________________
  subroutine copy_state(elem,nin,nout)
  implicit none
  
  type (element_t), intent(inout)   :: elem
  integer :: nin,nout

  elem%state%v(:,:,:,:,nout)  =elem%state%v(:,:,:,:,nin)
  elem%state%w(:,:,:,nout)    =elem%state%w(:,:,:,nin)
  elem%state%theta_dp_cp(:,:,:,nout) =elem%state%theta_dp_cp(:,:,:,nin)
  elem%state%phinh(:,:,:,nout)  =elem%state%phinh(:,:,:,nin)
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

  call get_dry_phinh(hvcoord,elem%state%phis,elem%state%theta_dp_cp(:,:,:,nt),dp,&
       elem%state%phinh(:,:,:,nt))


  ! debug
  kappa_star=kappa   
  call get_pnh_and_exner(hvcoord,elem%state%theta_dp_cp(:,:,:,nt),dp,&
       elem%state%phinh(:,:,:,nt),&
       elem%state%phis(:,:),kappa_star,pnh,dpnh,exner)
  do k=1,nlev
     if (maxval(abs(1-dpnh(:,:,k)/dp(:,:,k))) > 1e-10) then
        write(iulog,*) 'WARNING: hydrostatic inverse FAILED!'
        write(iulog,*) minval(dpnh(:,:,k)), minval(dp(:,:,k))
        write(iulog,*) k,minval(dpnh(:,:,k)/dp(:,:,k)),maxval(dpnh(:,:,k)/dp(:,:,k))
     endif
  enddo

  end subroutine set_thermostate



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
  elem%state%phinh(i,j,k, n0:n1) = g*zm
  elem%state%phis(i,j)           = phis
  elem%state%theta_dp_cp(i,j,k,n0:n1)=T*Cp*dp*((p/p0)**(-kappa))

  end subroutine set_state



  subroutine set_state_i(u,v,w,T,ps,phis,p,dp,zm,g,i,j,k,elem,n0,n1)
  !
  ! set state variables at node(i,j,k) at layer interfaces
  ! preqx model has no such variables, so do nothing
  !
  real(real_kind),  intent(in)    :: u,v,w,T,ps,phis,p,dp,zm,g
  integer,          intent(in)    :: i,j,k,n0,n1
  type(element_t),  intent(inout) :: elem

  end subroutine set_state_i




  !_____________________________________________________________________
  subroutine set_elem_state(u,v,w,w_i,T,ps,phis,p,dp,zm,zi,g,elem,n0,n1,ntQ)

  ! set element state variables
  ! works for both dry and moist initial conditions

  real(real_kind), dimension(np,np,nlev), intent(in):: u,v,w,T,p,dp,zm
  real(real_kind), dimension(np,np,nlevp), intent(in):: w_i,zi
  real(real_kind), dimension(np,np),      intent(in):: ps,phis
  real(real_kind),  intent(in)    :: g
  integer,          intent(in)    :: n0,n1,ntQ
  type(element_t),  intent(inout) :: elem
  integer :: n
  real(real_kind), dimension(np,np,nlev) :: cp_star,kappa_star

  ! get cp and kappa for dry or moist cases
  call get_cp_star(cp_star,elem%state%Q(:,:,:,1))
  call get_kappa_star(kappa_star,elem%state%Q(:,:,:,1))

  do n=n0,n1
    ! set prognostic state variables at level midpoints
    elem%state%v   (:,:,1,:,n) = u
    elem%state%v   (:,:,2,:,n) = v
    elem%state%w   (:,:,:,  n) = w
    elem%state%dp3d(:,:,:,  n) = dp
    elem%state%ps_v(:,:,    n) = ps
    elem%state%phinh(:,:,:, n) = g*zm
    elem%state%phis(:,:)       = phis
    elem%state%theta_dp_cp(:,:,:,n)=T*cp_star*dp*((p/p0)**(-kappa_star))
  end do

  end subroutine set_elem_state

  !_____________________________________________________________________
  subroutine get_state(u,v,w,T,pnh,dp,ps,rho,zm,zi,g,elem,hvcoord,nt,ntQ)

    ! get state variables at layer midpoints
    ! used by idealized tests to compute idealized physics forcing terms

    real(real_kind), dimension(np,np,nlev), intent(inout) :: u,v,w,T,pnh,dp,zm,rho
    real(real_kind), dimension(np,np,nlevp), intent(inout) :: zi
    real(real_kind), dimension(np,np),      intent(inout) :: ps
    real(real_kind), intent(in)    :: g
    integer,         intent(in)    :: nt,ntQ
    type(element_t), intent(inout) :: elem
    type (hvcoord_t),intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct

    real(real_kind) , dimension(np,np,nlev) :: phi,dpnh,kappa_star,exner,cp_star, temp
    real(real_kind) , dimension(np,np) :: phis

    integer :: k

    ! set prognostic state variables at level midpoints
    u   = elem%state%v   (:,:,1,:,nt)
    v   = elem%state%v   (:,:,2,:,nt)
    ps  = elem%state%ps_v(:,:,    nt)
    phis= elem%state%phis(:,:)
    w   = elem%state%w   (:,:,  :,nt)
    phi = elem%state%phinh(:,:, :,nt)

    do k=1,nlev
       dp(:,:,k)=(hvcoord%hyai(k+1)-hvcoord%hyai(k))*hvcoord%ps0 +(hvcoord%hybi(k+1)-hvcoord%hybi(k))*ps(:,:)
    enddo
    call get_cp_star(cp_star,elem%state%Q(:,:,:,1))
    call get_kappa_star(kappa_star,elem%state%Q(:,:,:,1))
    call get_pnh_and_exner(hvcoord,elem%state%theta_dp_cp(:,:,:,nt),dp,phi,phis,kappa_star,pnh,dpnh,exner)

    T     = elem%state%theta_dp_cp(:,:,:,nt)/(Cp_star*dp)*exner
    rho   = pnh/(kappa_star*cp_star*T)

    if(theta_hydrostatic_mode) then
       w = -(elem%derived%omega_p*pnh)/(rho*g)
       
    endif
    call get_phi(elem,zm,hvcoord,nt,ntQ); zm=zm/g
    call get_phi_i(elem,zi,hvcoord,nt,ntQ); zi=zi/g

    

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
  subroutine set_forcing_rayleigh_friction(elem, zm, zi, ztop, zc, tau, u0,v0, n)
  !
  ! test cases which use rayleigh friciton will call this with the relaxation coefficient
  ! f_d, and the reference state u0,v0.  Currently assume w0 = 0
  !
  implicit none

  type(element_t), intent(inout):: elem
  real(real_kind), intent(in)   :: zm(np,np,nlev) ! height at layer midpoints
  real(real_kind), intent(in)   :: zi(nlevp)      ! height at interfaces
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

real(real_kind), dimension(np,np,nlev) :: pnh,dpnh,exner

  do k=1,nlev
    dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k))*elem%state%ps_v(:,:,ns)
  enddo

  ntQ=1
  call get_kappa_star(kappa_star,elem%state%Q(:,:,:,1))
  call get_moist_phinh(hvcoord,elem%state%phis,elem%state%theta_dp_cp(:,:,:,ns),dp,kappa_star,elem%state%phinh(:,:,:,ns))

  ! verify discrete hydrostatic balance
  call get_pnh_and_exner(hvcoord,elem%state%theta_dp_cp(:,:,:,ns),dp,&
       elem%state%phinh(:,:,:,ns),&
       elem%state%phis(:,:),kappa_star,pnh,dpnh,exner)
  
  do k=1,nlev
     if (maxval(abs(1-dpnh(:,:,k)/dp(:,:,k))) > 1e-10) then
        write(iulog,*)'WARNING: hydrostatic inverse FAILED!'
        write(iulog,*) minval(dpnh(:,:,k)), minval(dp(:,:,k))
        write(iulog,*)k,minval(dpnh(:,:,k)/dp(:,:,k)),maxval(dpnh(:,:,k)/dp(:,:,k))
     endif
  enddo
  

  do tl = ns+1,ne
    call copy_state(elem,ns,tl)
  enddo

  if(present(ie)) call save_initial_state(elem%state,ie)


  end subroutine tests_finalize







  !_____________________________________________________________________
  subroutine get_kappa_star(kappa_star,Q)
  !
  ! note: interface written in this way so that it can be called outside
  ! timelevel loop,
  ! where dp is computed from reverence levels, or inside timelevel loop where
  ! dp = prognostic dp3d
  !
  implicit none
  real (kind=real_kind), intent(out)  :: kappa_star(np,np,nlev)
  real (kind=real_kind), intent(in)   :: Q(np,np,nlev)
  !   local
  integer :: k

  if (use_moisture .and. use_cpstar==1) then
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
     do k=1,nlev
        kappa_star(:,:,k) = (Rgas + (Rwater_vapor - Rgas)*Q(:,:,k))/ &
             (Cp + (Cpwater_vapor-Cp)*Q(:,:,k) )
     enddo
  else if (use_moisture .and. use_cpstar==0) then
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
     do k=1,nlev
        kappa_star(:,:,k) = (Rgas + (Rwater_vapor - Rgas)*Q(:,:,k))/ Cp
     enddo
  else
     kappa_star(:,:,:)=Rgas/Cp
  endif
  end subroutine


  !_____________________________________________________________________
  subroutine get_cp_star(cp_star,Q)
  !
  ! note: interface written in this way so that it can be called outside
  ! timelevel loop,
  ! where dp is computed from reverence levels, or inside timelevel loop where
  ! dp = prognostic dp3d
  !
  implicit none
  real (kind=real_kind), intent(out):: cp_star(np,np,nlev)
  real (kind=real_kind), intent(in) :: Q(np,np,nlev)

  integer :: k
  if (use_moisture .and. use_cpstar==1) then
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
     do k=1,nlev
        cp_star(:,:,k) = (Cp + (Cpwater_vapor-Cp)*Q(:,:,k) )
     enddo
  else
     cp_star(:,:,:)=Cp
  endif
  end subroutine






  !_____________________________________________________________________
  subroutine set_theta_ref(hvcoord,dp,theta_ref)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! a reference profile for theta = theta(exner)
  !
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




end module

