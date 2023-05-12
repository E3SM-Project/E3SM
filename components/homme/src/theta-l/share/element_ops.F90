#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
!
!  getter and setter functions that must be provided by each model
!
!  Note: all routines require dp3d() to be valid, with the exception of
!  the initial condition routines which assume reference levels and will initialize dp3d based
!  on their 'ps' input argument
!
!
!
! ROUTINES REQUIRED FOR ALL MODELS:
! Initial condition routines
!  set_thermostate()
!     initial condition interface used by DCMIP 2008 tests, old HOMME tests
!  set_state(), set_state_i()
!     initial condition interface used by DCMIP 2012 tests
!  set_elem_state()
!     initial condition interface used by DCMIP 2016 tests
!  tests_finalize()
!     initialize geopotential to be in hydrostatic balance
!     used by DCMIP2012, 2016 tests
! Other routines:
!  get_field()
!     returns temperature, potential temperature, phi, etc..
!  get_field_i()
!     returns a few quantities on interfaces
!  copy_state()
!     copy state variables from one timelevel to another timelevel
!  get_state()
!     return state variables used by some DCMIP forcing functions
!  save_initial_state()
!     save t=0 in "state0", used by some DCMIP forcing functions
!  set_forcing_rayleigh_friction()
!     used by dcmip2012 test cases
!  set_forcing_rayleigh_friction()
!     apply rayleigh friction to prognostic variables, used by some DCMIP2016 tests
!
!  get_temperature()   used in CAM dp_coupling layer
!
! UTILITY ROUTINES USED BY THETA-L MODEL
!  get_pottemp()
!  get_dpnh_dp()
!  get_hydro_pressure()
!  get_nonhydro_pressure()
!  get_phi()
!  get_cp_star()
!  get_R_star()
!  set_theta_ref()
!
!
module element_ops

  use dimensions_mod, only: np, nlev, nlevp, nelemd
  use element_mod,    only: element_t
  use element_state,  only: elem_state_t, timelevels
  use hybvcoord_mod,  only: hvcoord_t
  use kinds,          only: real_kind, iulog
  use perf_mod,       only: t_startf, t_stopf, t_barrierf, t_adj_detailf ! _EXTERNAL
  use parallel_mod,   only: abortmp
  use physical_constants, only : p0, Cp, Rgas, Rwater_vapor, Cpwater_vapor, kappa, g, dd_pi, TREF
  use control_mod,    only: use_moisture, theta_hydrostatic_mode, hv_ref_profiles
  use eos,            only: pnh_and_exner_from_eos, phi_from_eos
  implicit none
  private

  type(elem_state_t), dimension(:), allocatable :: state0 ! storage for save_initial_state routine

  public get_field, get_field_i, get_state, get_pottemp
  public get_temperature, get_phi, get_R_star, get_hydro_pressure
  public set_thermostate, set_state, set_state_i, set_elem_state
  public set_forcing_rayleigh_friction
  public initialize_reference_states, set_theta_ref
  public copy_state, tests_finalize
  public state0

  ! promote this to _real_kind after V2 code freeze
  real (kind=real_kind), public :: tref_lapse_rate=0.0065D0
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
  real(kind=real_kind), dimension(np,np,nlev) :: tmp, p, pnh, dp, omega, rho, T, cp_star, Rstar
  real(kind=real_kind), dimension(np,np,nlevp) :: phi_i,pi_i


  select case(name)
    case ('temperature','T'); call get_temperature(elem,field,hvcoord,nt)
    case ('pottemp','Th');    call get_pottemp(elem,field,hvcoord,nt,ntQ)
    case ('phi','geo');       call get_phi(elem,field,phi_i,hvcoord,nt)
    case ('dpnh_dp');         call get_dpnh_dp(elem,field,hvcoord,nt)
    case ('pnh');             call get_nonhydro_pressure(elem,field,tmp ,hvcoord,nt)
    case ('exner');           call get_nonhydro_pressure(elem,tmp ,field,hvcoord,nt)

    case ('p');
      call get_hydro_pressure(field,elem%state%dp3d(:,:,:,nt),hvcoord)

    case ('dp');
       field(:,:,:)=elem%state%dp3d(:,:,:,nt)


    case ('omega');
      field = elem%derived%omega_p

    case('rho')

      call get_nonhydro_pressure(elem,pnh,tmp,hvcoord,nt)
      call get_R_star(Rstar,elem%state%Q(:,:,:,1))
      call get_temperature(elem,T,hvcoord,nt)
      field = pnh/(Rstar*T)

    case ('w');

      if(theta_hydrostatic_mode) then
        call get_field(elem,'omega',omega,hvcoord,nt,ntQ)
        call get_field(elem,'rho'  ,rho  ,hvcoord,nt,ntQ)
        field = -omega/(rho *g)

      else
        field =elem%state%w_i(:,:,1:nlev,nt)
      endif

    case default
       print *,'name = ',trim(name)
       call abortmp('ERROR: get_field name not supported in this model')

  end select

  end subroutine


  subroutine get_field_i(elem,name,field,hvcoord,nt)
  implicit none
  type (element_t),       intent(in) :: elem
  character(len=*),       intent(in) :: name
  real (kind=real_kind),  intent(out):: field(np,np,nlevp)
  integer,                intent(in) :: nt
  type (hvcoord_t),       intent(in) :: hvcoord

  select case(name)
    case('w_i');
      if(theta_hydrostatic_mode) then
         call abortmp('ERROR: get_field_i is not supported for w in theta HY')
      else
         field = elem%state%w_i(:,:,1:nlevp,nt)
      endif
    case('geo_i');
      if(theta_hydrostatic_mode) then
         call phi_from_eos(hvcoord,elem%state%phis,elem%state%vtheta_dp(:,:,:,nt),elem%state%dp3d(:,:,:,nt),field)
      else
          field = elem%state%phinh_i(:,:,1:nlevp,nt)
      endif
    case('mu_i');
      call get_dpnh_dp_i(elem,field,hvcoord,nt)
    case('pnh_i');
      call get_nonhydro_pressure_i(elem,field,hvcoord,nt)
    case default
      print *,'name = ',trim(name)
      call abortmp('ERROR: get_field_i name not supported in this model')
  end select

  end subroutine get_field_i


  !_____________________________________________________________________
  subroutine get_pottemp(elem,pottemp,hvcoord,nt,ntQ)
  !
  implicit none

  type (element_t), intent(in)        :: elem
  real (kind=real_kind), intent(out)  :: pottemp(np,np,nlev)
  type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
  integer, intent(in) :: nt
  integer, intent(in) :: ntQ

  !   local
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: Rstar(np,np,nlev)
  integer :: k

  call get_R_star(Rstar,elem%state%Q(:,:,:,1))

  pottemp(:,:,:) = Rgas*elem%state%vtheta_dp(:,:,:,nt)/(Rstar(:,:,:)*elem%state%dp3d(:,:,:,nt))

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
  real (kind=real_kind) :: Rstar(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev)
  real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)
  integer :: k


  dp=elem%state%dp3d(:,:,:,nt)
  call get_R_star(Rstar,elem%state%Q(:,:,:,1))

  call pnh_and_exner_from_eos(hvcoord,elem%state%vtheta_dp(:,:,:,nt),&
          dp,elem%state%phinh_i(:,:,:,nt),pnh,exner,dpnh_dp_i)

  do k=1,nlev
     temperature(:,:,k)= Rgas*elem%state%vtheta_dp(:,:,k,nt)*exner(:,:,k)&
          /(Rstar(:,:,k)*dp(:,:,k))
  enddo

  end subroutine get_temperature

!this routine averages mu to midlevels
  !_____________________________________________________________________
  subroutine get_dpnh_dp(elem,dpnh_dp,hvcoord,nt)
  implicit none

  type (element_t), intent(in)        :: elem
  real (kind=real_kind), intent(out)  :: dpnh_dp(np,np,nlev)
  type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
  integer, intent(in) :: nt

  !   local
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev)
  real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)
  integer :: k


  dp=elem%state%dp3d(:,:,:,nt)
  call pnh_and_exner_from_eos(hvcoord,elem%state%vtheta_dp(:,:,:,nt),&
       dp,elem%state%phinh_i(:,:,:,nt),pnh,exner,dpnh_dp_i)

  do k=1,nlev
     dpnh_dp(:,:,k)=(dpnh_dp_i(:,:,k)+dpnh_dp_i(:,:,k+1))/2
  enddo
  end subroutine

!this routine returns mu at interfaces, mu_surface = 1, will be fixed later
  !_____________________________________________________________________
  subroutine get_dpnh_dp_i(elem,dpnh_dp_i,hvcoord,nt)
  implicit none

  type (element_t), intent(in)        :: elem
  real (kind=real_kind), intent(out)  :: dpnh_dp_i(np,np,nlevp)
  type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
  integer, intent(in) :: nt

  !   local
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev)
  integer :: k

  dp=elem%state%dp3d(:,:,:,nt)
  call pnh_and_exner_from_eos(hvcoord,elem%state%vtheta_dp(:,:,:,nt),&
       dp,elem%state%phinh_i(:,:,:,nt),pnh,exner,dpnh_dp_i)

  end subroutine get_dpnh_dp_i

  !_____________________________________________________________________
  subroutine get_hydro_pressure(p,dp,hvcoord)
  !
  implicit none

  real (kind=real_kind), intent(out)  :: p(np,np,nlev)
  real (kind=real_kind), intent(in)   :: dp(np,np,nlev)
  type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct

  integer :: k
  real(kind=real_kind), dimension(np,np,nlevp) :: p_i

  p_i(:,:,1)=hvcoord%hyai(1)*hvcoord%ps0
  do k=1,nlev  ! SCAN
     p_i(:,:,k+1)=p_i(:,:,k) + dp(:,:,k)
  enddo
#ifdef HOMMEXX_BFB_TESTING
  do k=1,nlev
     p(:,:,k) = (p_i(:,:,k+1)+p_i(:,:,k))/2
  enddo
#else
  do k=1,nlev
     p(:,:,k)=p_i(:,:,k) + dp(:,:,k)/2
  enddo
#endif


  end subroutine get_hydro_pressure


  !_____________________________________________________________________
  subroutine get_nonhydro_pressure(elem,pnh,exner,hvcoord,nt)
    implicit none

    type (element_t),       intent(in)  :: elem
    real (kind=real_kind),  intent(out) :: pnh(np,np,nlev)
    real (kind=real_kind),  intent(out) :: exner(np,np,nlev)
    type (hvcoord_t),       intent(in)  :: hvcoord
    integer,                intent(in)  :: nt

    real (kind=real_kind), dimension(np,np,nlevp) :: dpnh_dp_i

    call pnh_and_exner_from_eos(hvcoord,elem%state%vtheta_dp(:,:,:,nt),&
         elem%state%dp3d(:,:,:,nt),elem%state%phinh_i(:,:,:,nt),&
         pnh,exner,dpnh_dp_i,caller="get_nonhydro_pressire")

  end subroutine


  subroutine get_nonhydro_pressure_i(elem,pnh_i,hvcoord,nt)
    implicit none

    type (element_t),       intent(in)  :: elem
    real (kind=real_kind),  intent(out) :: pnh_i(np,np,nlevp)
    type (hvcoord_t),       intent(in)  :: hvcoord
    integer,                intent(in)  :: nt

    real (kind=real_kind), dimension(np,np,nlevp) :: dpnh_dp_i
    real (kind=real_kind), dimension(np,np,nlev) :: pnh,exner

    call pnh_and_exner_from_eos(hvcoord,elem%state%vtheta_dp(:,:,:,nt),&
         elem%state%dp3d(:,:,:,nt),elem%state%phinh_i(:,:,:,nt),&
         pnh,exner,dpnh_dp_i,pnh_i,"get_nonhydro_pressure_i")

  end subroutine


  subroutine get_phi(elem,phi,phi_i,hvcoord,nt)
    implicit none

    type (element_t),       intent(in)  :: elem
    type (hvcoord_t),       intent(in)  :: hvcoord
    real (kind=real_kind),  intent(out) :: phi(np,np,nlev)
    real (kind=real_kind),  intent(out) :: phi_i(np,np,nlevp)
    integer,                intent(in)  :: nt

    real (kind=real_kind), dimension(np,np,nlev) :: dp
    real (kind=real_kind) :: pnh(np,np,nlev)
    real (kind=real_kind) :: exner(np,np,nlev)
    real (kind=real_kind) :: temp(np,np,nlev)
    real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)
    integer :: k


    if(theta_hydrostatic_mode) then
       dp=elem%state%dp3d(:,:,:,nt)
       call phi_from_eos(hvcoord,elem%state%phis,elem%state%vtheta_dp(:,:,:,nt),dp,phi_i)
    else
       phi_i = elem%state%phinh_i(:,:,:,nt)
    endif

    do k=1,nlev
       phi(:,:,k) = (phi_i(:,:,k)+phi_i(:,:,k+1))/2
    end do

  end subroutine


  !_____________________________________________________________________
  subroutine get_cp_star(cp_star,Q)
  !
  !
  implicit none
  real (kind=real_kind), intent(out):: cp_star(np,np,nlev)
  real (kind=real_kind), intent(in) :: Q(np,np,nlev)

  integer :: k
  if (use_moisture) then
     do k=1,nlev
        cp_star(:,:,k) = (Cp + (Cpwater_vapor-Cp)*Q(:,:,k) )
     enddo
  else
     cp_star(:,:,:)=Cp
  endif
  end subroutine


  !_____________________________________________________________________
  subroutine get_R_star(R_star,Q)
  !
  implicit none
  real (kind=real_kind), intent(out):: R_star(np,np,nlev)
  real (kind=real_kind), intent(in) :: Q(np,np,nlev)

  integer :: k
  if (use_moisture) then
     do k=1,nlev
        R_star(:,:,k) =(Rgas + (Rwater_vapor - Rgas)*Q(:,:,k))
     enddo
  else
     R_star(:,:,:)=Rgas
  endif
  end subroutine


  !_____________________________________________________________________
  subroutine copy_state(elem,nin,nout)
  implicit none

  type (element_t), intent(inout)   :: elem
  integer :: nin,nout

  elem%state%v(:,:,:,:,nout)         =elem%state%v(:,:,:,:,nin)
  elem%state%w_i(:,:,:,nout)         =elem%state%w_i(:,:,:,nin)
  elem%state%vtheta_dp(:,:,:,nout)   =elem%state%vtheta_dp(:,:,:,nin)
  elem%state%phinh_i(:,:,:,nout)     =elem%state%phinh_i(:,:,:,nin)
  elem%state%dp3d(:,:,:,nout)        =elem%state%dp3d(:,:,:,nin)
  elem%state%ps_v(:,:,nout)          =elem%state%ps_v(:,:,nin)
  end subroutine copy_state


  !_____________________________________________________________________
  subroutine set_thermostate(elem,ps,temperature,hvcoord,qv)
  !
  ! Assuming a hydrostatic intital state and given surface pressure,
  ! and no moisture, compute theta and phi
  !
  ! input:  ps_v, temperature
  ! ouput:  state variables:   vtheta_dp, phi
  !
  implicit none

  type (element_t), intent(inout)   :: elem
  real (kind=real_kind), intent(in) :: temperature(np,np,nlev)
  type (hvcoord_t),     intent(in)  :: hvcoord                      ! hybrid vertical coordinate struct
  real (kind=real_kind), intent(in) :: ps(np,np)
  real (kind=real_kind), intent(in), optional :: qv(np,np,nlev)  ! water vapor mixing ratio

  !   local
  real (kind=real_kind) :: p(np,np,nlev)
  real (kind=real_kind) :: pi_i(np,np,nlevp)
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: Rstar(np,np,nlev)
  integer :: k, nt

  nt = 1
  do k=1,nlev
     p(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps(:,:)
     dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*ps(:,:)
  enddo

! recompute pressure using same algrebra as EOS:
! improves precision for initial T/theta conversion
  pi_i(:,:,1)=hvcoord%hyai(1)*hvcoord%ps0
  do k=1,nlev
     pi_i(:,:,k+1)=pi_i(:,:,k) + dp(:,:,k)
  enddo
  do k=1,nlev
     p(:,:,k)=pi_i(:,:,k) + dp(:,:,k)/2
  enddo


!set vtheta
  do k=1,nlev
     elem%state%vtheta_dp(:,:,k,nt)=dp(:,:,k)*temperature(:,:,k)* &
          (p(:,:,k)/p0)**(-kappa)
     elem%state%dp3d(:,:,k,nt)=dp(:,:,k)
  enddo
  elem%state%ps_v(:,:,nt)=ps

  if (present(qv)) then
     call get_R_star(Rstar,qv)
     elem%state%vtheta_dp(:,:,:,nt)=elem%state%vtheta_dp(:,:,:,nt)*Rstar(:,:,:)/Rgas
  endif


  !set phi, copy from 1st timelevel to all
  call tests_finalize(elem,hvcoord)


  ! verify T
  call get_temperature(elem,p,hvcoord,nt)

  dp=p-temperature
  do k=1,nlev
     if (maxval(abs(dp(:,:,k))) > 1e-9) then
        write(iulog,*)'WARNING: T/theta initialization inconsistent!'
        write(iulog,*)k,minval(dp(:,:,k)),maxval(dp(:,:,k))
        write(iulog,*) 'T,Tnew',temperature(1,1,k),p(1,1,k)
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
  elem%state%v   (i,j,1,k,n0:n1)   = u
  elem%state%v   (i,j,2,k,n0:n1)   = v
  elem%state%dp3d(i,j,k,  n0:n1)   = dp
  elem%state%ps_v(i,j,    n0:n1)   = ps
  elem%state%phis(i,j)             = phis
  elem%state%vtheta_dp(i,j,k,n0:n1)=T*dp*((p/p0)**(-kappa))

  end subroutine set_state


  subroutine set_state_i(u,v,w,T,ps,phis,p,zm,g,i,j,k,elem,n0,n1)
  !
  ! set state variables at node(i,j,k) at layer interfaces
  ! used by idealized tests for dry initial conditions
  ! so we use constants cp, kappa
  !
  real(real_kind),  intent(in)    :: u,v,w,T,ps,phis,p,zm,g
  integer,          intent(in)    :: i,j,k,n0,n1
  type(element_t),  intent(inout) :: elem

  ! set prognostic state variables at level midpoints
  elem%state%w_i  (i,j,k,  n0:n1)  = w
  elem%state%phinh_i(i,j,k, n0:n1) = g*zm

  end subroutine set_state_i

  !_____________________________________________________________________
  subroutine set_elem_state(u,v,w,w_i,T,ps,phis,p,dp,zm,zi,g,elem,n0,n1,ntQ)
  !
  ! set element state variables
  ! works for both dry and moist initial conditions
  !
  !
  real(real_kind), dimension(np,np,nlev), intent(in):: u,v,w,T,p,dp,zm
  real(real_kind), dimension(np,np,nlevp), intent(in):: w_i,zi
  real(real_kind), dimension(np,np),      intent(in):: ps,phis
  real(real_kind),  intent(in)    :: g
  integer,          intent(in)    :: n0,n1,ntQ
  type(element_t),  intent(inout) :: elem
  integer :: n
  real(real_kind), dimension(np,np,nlev) :: Rstar

  ! get cp and kappa for dry or moist cases
  call get_R_star(Rstar,elem%state%Q(:,:,:,1))

  do n=n0,n1
    ! set prognostic state variables at level midpoints
    elem%state%v   (:,:,1,:,n)        = u
    elem%state%v   (:,:,2,:,n)        = v
    elem%state%dp3d(:,:,:,  n)        = dp
    elem%state%ps_v(:,:,    n)        = ps
    elem%state%phis(:,:)              = phis
    elem%state%vtheta_dp(:,:,:,n)   = (Rstar/Rgas)*T*dp*((p/p0)**(-kappa))

    elem%state%w_i (:,:,:,  n)   = w_i
    elem%state%phinh_i(:,:,:, n) = g*zi
  end do

  end subroutine set_elem_state

  !_____________________________________________________________________
  subroutine get_state(u,v,w,T,pnh,dp,ps,rho,zm,zi,g,elem,hvcoord,nt,ntQ)
    ! get state variables at layer midpoints
    ! used by idealized tests to compute idealized physics forcing terms
    ! currently all forcing is done on u,v and T/theta - no forcing
    ! is computed for interface variables.   This routine will have to be updated
    ! if we add a test case that computes forcing for interface variables

    real(real_kind), dimension(np,np,nlev), intent(inout) :: u,v,w,T,pnh,dp,zm,rho
    real(real_kind), dimension(np,np,nlevp), intent(inout) :: zi
    real(real_kind), dimension(np,np),      intent(inout) :: ps
    real(real_kind), intent(in)    :: g
    integer,         intent(in)    :: nt,ntQ
    type(element_t), intent(inout) :: elem
    type (hvcoord_t),intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct

    real(real_kind) , dimension(np,np,nlev) :: Rstar,exner,temp
    real(real_kind) , dimension(np,np,nlevp):: dpnh_dp_i,phi_i

    integer :: k

    u   = elem%state%v   (:,:,1,:,nt)
    v   = elem%state%v   (:,:,2,:,nt)
    ps  = elem%state%ps_v(:,:,  nt)
    dp  = elem%state%dp3d(:,:,:,nt)
    phi_i = elem%state%phinh_i(:,:,:,nt)

    do k=1,nlev
       w(:,:,k) = (elem%state%w_i(:,:,k,nt) + elem%state%w_i(:,:,k+1,nt))/2
    end do

    call get_R_star(Rstar,elem%state%Q(:,:,:,1))
    call pnh_and_exner_from_eos(hvcoord,elem%state%vtheta_dp(:,:,:,nt),dp,phi_i,&
         pnh,exner,dpnh_dp_i)

    ! first compute virtual temperature Tv needed for rho
    T     = exner*elem%state%vtheta_dp(:,:,:,nt)/(dp)
    rho   = pnh/(Rgas*T)

    ! convert Tv->T   Tv = R* T / R
    T = Rgas*T/Rstar

    if(theta_hydrostatic_mode) then
       ! overwrite w and phi_i computed above
       w = -(elem%derived%omega_p)/(rho*g)

       do k=nlev,1,-1
          temp(:,:,k) = Rgas*elem%state%vtheta_dp(:,:,k,nt)*exner(:,:,k)/pnh(:,:,k)
          phi_i(:,:,k)=phi_i(:,:,k+1)+temp(:,:,k)
       enddo
    endif

    do k=1,nlev
       zm(:,:,k) = (phi_i(:,:,k)+phi_i(:,:,k+1))/(2*g)
    end do
    do k=1,nlevp
       zi(:,:,k) = phi_i(:,:,k)/g
    end do

  end subroutine get_state

  !_____________________________________________________________________
  subroutine save_initial_state(state,ie)

    ! save the state at time=0 for use with some dcmip tests

    type(elem_state_t), intent(inout):: state
    integer,            intent(in)   :: ie     ! element index

!$OMP BARRIER
!$OMP MASTER
    if(.not. allocated(state0)) allocate( state0(nelemd) )
!$OMP END MASTER
!$OMP BARRIER
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
  real(real_kind), intent(in)   :: zi(np,np,nlevp) ! height at interfaces
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

  ! compute forcing for 1:nlev. w is determined by boundary condition at nlevp:
  f_d=0.0d0
  where(zi(:,:,1:nlev) .ge. zc);   f_d = sin(dd_pi/2 *(zm - zc)/(ztop - zc))**2; end where
  where(zi(:,:,1:nlev) .ge. ztop); f_d = 1.0d0; end where
  f_d = -f_d/tau
  elem%derived%FM(:,:,3,:) = f_d * ( elem%state%w_i(:,:,1:nlev,n)  )
  end subroutine

  !_____________________________________________________________________
  subroutine tests_finalize(elem,hvcoord,ie)

  ! Now that all variables have been initialized, set phi to be in hydrostatic balance

  implicit none

  type(hvcoord_t),     intent(in)   :: hvcoord
  type(element_t),     intent(inout):: elem
  integer, optional,   intent(in)   :: ie ! optional element index, to save initial state

  integer :: k,tl
  real(real_kind), dimension(np,np,nlev) :: pi

  real(real_kind), dimension(np,np,nlev) :: pnh,exner
  real(real_kind), dimension(np,np,nlevp) :: dpnh_dp_i,phi_i

  tl=1

  call phi_from_eos(hvcoord,elem%state%phis,elem%state%vtheta_dp(:,:,:,tl),&
       elem%state%dp3d(:,:,:,tl),elem%state%phinh_i(:,:,:,tl))

  ! Disable the following check in CUDA bfb builds,
  ! since the calls to pow are inexact
#if !(defined(HOMMEXX_BFB_TESTING) && defined(HOMMEXX_ENABLE_GPU))
  ! verify discrete hydrostatic balance
  call pnh_and_exner_from_eos(hvcoord,elem%state%vtheta_dp(:,:,:,tl),&
       elem%state%dp3d(:,:,:,tl),elem%state%phinh_i(:,:,:,tl),pnh,exner,dpnh_dp_i)
  do k=1,nlev
     pi(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem%state%ps_v(:,:,tl)
     if (maxval(abs(1-dpnh_dp_i(:,:,k))) > 1e-10) then
        write(iulog,*)'WARNING: hydrostatic inverse FAILED!'
        write(iulog,*)k,minval(dpnh_dp_i(:,:,k)),maxval(dpnh_dp_i(:,:,k))
        write(iulog,*) 'pnh',pi(1,1,k),pnh(1,1,k)
     endif
  enddo
#endif

  do tl = 2,timelevels
    call copy_state(elem,1,tl)
  enddo

  if(present(ie)) call save_initial_state(elem%state,ie)


  end subroutine tests_finalize

  !_____________________________________________________________________
  subroutine initialize_reference_states(hvcoord, phis, &
                                         dp_ref, theta_ref, phi_ref)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Compute and sets reference profiles
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   implicit none

   ! input variables
   type (hvcoord_t),      intent(in) :: hvcoord     ! hybrid vertical coordinate struct
   real (kind=real_kind), intent(in) :: phis(np,np) ! surface geopotential

   ! output variables
   real (kind=real_kind), intent(out) :: dp_ref(np,np,nlev)
   real (kind=real_kind), intent(out) :: theta_ref(np,np,nlev)
   real (kind=real_kind), intent(out) :: phi_ref(np,np,nlevp)

   ! local variables
   integer :: k
   real (kind=real_kind) :: ps_ref(np,np), temp(np,np,nlev)

   ! compute dp_ref
   ps_ref(:,:) = hvcoord%ps0*exp(-phis(:,:)/(Rgas*TREF))
   do k=1,nlev
     dp_ref(:,:,k) = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + &
                     (hvcoord%hybi(k+1) - hvcoord%hybi(k))*ps_ref(:,:)
   enddo

   ! compute theta_ref
   call set_theta_ref(hvcoord, dp_ref, theta_ref,0)

   ! compute phi_ref
   temp = theta_ref*dp_ref
   call phi_from_eos(hvcoord, phis, temp, dp_ref, phi_ref)

   ! keep profiles, based on the value of hv_ref_profiles
   if (hv_ref_profiles == 0) then
     ! keep phi profile, but dont use theta and dp:
     theta_ref = 0
     dp_ref = 0
   endif
   if (hv_ref_profiles == 1) then
     ! keep all profiles
   endif
   if (hv_ref_profiles == 6) then
      ! keep all profiles, use linearized theta_ref profile
      call set_theta_ref(hvcoord, dp_ref, theta_ref,1)
   endif

   end subroutine initialize_reference_states


  !_____________________________________________________________________
  subroutine set_theta_ref(hvcoord,dp,theta_ref,linear_profile)
#ifdef HOMMEXX_BFB_TESTING
    use bfb_mod, only: bfb_pow
#endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! a reference profile for theta = theta(exner)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none

  type (hvcoord_t),     intent(in)  :: hvcoord                      ! hybrid vertical coordinate struct
  real (kind=real_kind), intent(in) :: dp(np,np,nlev)
  real (kind=real_kind), intent(out) :: theta_ref(np,np,nlev)

  integer :: linear_profile
  !   local
  real (kind=real_kind) :: p_i(np,np,nlevp)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: T0,T1
  integer :: k

  ! reference T = 288K.  reference lapse rate = 6.5K/km   = .0065 K/m
  ! Tref = T0+T1*exner
  ! Thetaref = T0/exner + T1
  T1 = tref_lapse_rate*TREF*Cp/g ! = 191
  T0 = TREF-T1           ! = 97

  p_i(:,:,1) =  hvcoord%hyai(1)*hvcoord%ps0
  do k=1,nlev
     p_i(:,:,k+1) = p_i(:,:,k) + dp(:,:,k)
  enddo
  do k=1,nlev
#ifdef HOMMEXX_BFB_TESTING
     exner(:,:,k) = bfb_pow(((p_i(:,:,k) + p_i(:,:,k+1))/2)/p0,kappa)
#else
     exner(:,:,k) = ( (p_i(:,:,k) + p_i(:,:,k+1))/(2*p0)) **kappa
#endif
     !theta_ref(:,:,k,ie) = (T0/exner(:,:,k) + T1)*Cp*dp_ref(:,:,k,ie)
     if (linear_profile==1) then
        ! linearize around (1-exner)
        theta_ref(:,:,k) = T0 + T0*(1-exner(:,:,k)) + T1   
     else
        theta_ref(:,:,k) = (T0/exner(:,:,k) + T1)
     endif
  enddo

  end subroutine




end module

