#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!
!  getter & setter functions that must be provided by each model
!  
!  IMPORTANT NOTE:  For vertically lagrangian models, these
!  routines should ONLY be used outside the timestepping loop
!  on reference levels.  To compute these fields on floating levels
!  the model should do that directly (or we need to modify the interface)
!
!  get_field() 
!     returns temperature, potential temperature, phi, etc..
!
! These should be unified to a single interface:
!  set_thermostate()    
!     initial condition interface used by DCMIP 2008 tests
!     
!  set_state()
!     initial condition interface used by DCMIP 2012 tests
!
!
module element_ops

  use control_mod,    only: use_moisture, use_cpstar
  use dimensions_mod, only: np, nlev, nlevp, nelemd, qsize, max_corner_elem
  use element_mod,    only: element_t
  use element_state,  only: timelevels, elem_state_t
  use hybrid_mod,     only: hybrid_t
  use hybvcoord_mod,  only: hvcoord_t
  use kinds,          only: real_kind, iulog
  use perf_mod,       only: t_startf, t_stopf, t_barrierf, t_adj_detailf ! _EXTERNAL
  use parallel_mod,   only: abortmp
  use physical_constants, only : kappa, p0, Rgas, cp, g, dd_pi, Rwater_vapor, Cpwater_vapor
  use physics_mod,    only : virtual_temperature
  use prim_si_mod,    only : preq_hydrostatic
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

  real(real_kind), dimension(np,np,nlev) :: p, T, omega, rho
  integer :: k

  select case(name)

    case ('T','temperature'); call get_temperature(elem,field,hvcoord,nt)
    case ('Th','pottemp');    call get_pottemp(elem,field,hvcoord,nt,ntQ)
    case ('geo','phi');       call get_phi(elem,field,hvcoord,nt,ntQ)

    case ('omega')
        call get_field(elem,'p',p,hvcoord,nt,ntQ)
        field =elem%derived%omega_p*p

    case ('p','pnh');
      forall(k=1:nlev) field(:,:,k)=hvcoord%hyam(k)*hvcoord%ps0+hvcoord%hybm(k)*elem%state%ps_v(:,:,nt)

    case('rho'); ! get rho from dry air equation of state
      call get_field(elem,'p',p,hvcoord,nt,ntQ)
      call get_field(elem,'T',T,hvcoord,nt,ntQ)
      field = p/(Rgas*T)

    case ('w'); ! get w from omega using hydrostatic balance condition
      call get_field(elem,'rho',rho,hvcoord,nt,ntQ)
      call get_field(elem,'omega',omega,hvcoord,nt,ntQ)
      field = -omega/(rho*g)

    case('zeta'); ! todo

  case default
     print *,'name = ',trim(name)
     call abortmp('ERROR: get_field name not supported in this model')
  end select

  end subroutine

  !_____________________________________________________________________
  subroutine get_phi(elem,phi,hvcoord,nt,ntQ)
  implicit none
    
  type (element_t), intent(in)        :: elem
  real (kind=real_kind), intent(out)  :: phi(np,np,nlev)
  type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
  integer, intent(in) :: nt
  integer, intent(in) :: ntQ
  
  !   local
  real (kind=real_kind) :: pfull(np,np,nlev)
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: T_v(np,np,nlev)
  real (kind=real_kind) :: Qt
  integer :: k,i,j


  do k=1,nlev
     pfull(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0  &
          + hvcoord%hybm(k)*elem%state%ps_v(:,:,nt)
     dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem%state%ps_v(:,:,nt)


     if (.not. use_moisture ) then
        T_v(:,:,k) = elem%state%T(:,:,k,nt)
     else
        do j=1,np
           do i=1,np
              Qt = elem%state%Qdp(i,j,k,1,ntQ)/dp(i,j,k)
              T_v(i,j,k) = Virtual_Temperature(elem%state%T(i,j,k,nt),Qt)
           end do
        end do
     endif
  enddo
  call preq_hydrostatic(phi,elem%state%phis,T_v,pfull,dp)  

  end subroutine get_phi


  !_____________________________________________________________________
  subroutine get_pottemp(elem,pottemp,hvcoord,nt,ntQ)
  implicit none
    
  type (element_t), intent(in)        :: elem
  real (kind=real_kind), intent(out)  :: pottemp(np,np,nlev)
  type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
  integer, intent(in) :: nt
  integer, intent(in) :: ntQ
  
  !   local
  real (kind=real_kind) :: pfull(np,np,nlev)
  integer :: k


#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     pfull(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0  &
          + hvcoord%hybm(k)*elem%state%ps_v(:,:,nt)
     pottemp(:,:,k)=elem%state%T(:,:,k,nt)* &
          (pfull(:,:,k)/p0)**(-kappa)
  enddo
  
  end subroutine get_pottemp


  !_____________________________________________________________________
  subroutine get_temperature(elem,temperature,hvcoord,nt)
  implicit none
  
  type (element_t), intent(in)        :: elem
  real (kind=real_kind), intent(out)  :: temperature(np,np,nlev)
  type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
  integer, intent(in) :: nt
  
  temperature = elem%state%T(:,:,:,nt)
  
  end subroutine get_temperature

  !_____________________________________________________________________
  subroutine copy_state(elem,nin,nout)
  implicit none
  
  type (element_t), intent(inout)   :: elem
  integer :: nin,nout

  elem%state%v(:,:,:,:,nout) =elem%state%v(:,:,:,:,nin)
  elem%state%T(:,:,:,nout)   =elem%state%T(:,:,:,nin)
  elem%state%dp3d(:,:,:,nout)=elem%state%dp3d(:,:,:,nin)
  elem%state%ps_v(:,:,nout)  =elem%state%ps_v(:,:,nin)
  end subroutine copy_state


  !_____________________________________________________________________
  subroutine set_thermostate(elem,temperature,hvcoord,n0,n0_q)
  implicit none
  
  type (element_t), intent(inout)   :: elem
  real (kind=real_kind), intent(in) :: temperature(np,np,nlev)
  type (hvcoord_t),     intent(in)  :: hvcoord                      ! hybrid vertical coordinate struct
  integer :: n0,n0_q

  elem%state%T(:,:,:,n0)=temperature(:,:,:)

  end subroutine set_thermostate

  !_____________________________________________________________________
  subroutine set_state(u,v,w,T,ps,phis,p,dp,zm, g,  i,j,k,elem,n0,n1)

    ! set state variables at node(i,j,k) at layer midpoints

    real(real_kind),  intent(in)    :: u,v,w,T,ps,phis,p,dp,zm,g
    integer,          intent(in)    :: i,j,k,n0,n1
    type(element_t),  intent(inout) :: elem

    ! set prognostic state variables at level midpoints
    elem%state%v   (i,j,1,k,n0:n1) = u
    elem%state%v   (i,j,2,k,n0:n1) = v
    elem%state%T   (i,j,k,n0:n1)   = T
    elem%state%ps_v(i,j,n0:n1)     = ps
    elem%state%phis(i,j)           = phis

  end subroutine


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

    ! set state variables for entire element

    real(real_kind), dimension(np,np,nlev), intent(in):: u,v,w,T,p,dp,zm
    real(real_kind), dimension(np,np,nlevp), intent(in):: w_i,zi
    real(real_kind), dimension(np,np),      intent(in):: ps,phis
    real(real_kind),  intent(in)    :: g
    integer,          intent(in)    :: n0,n1,ntQ
    type(element_t),  intent(inout) :: elem
    integer :: n
    real(real_kind), dimension(np,np,nlev) :: cp_star,kappa_star

    do n=n0,n1
      ! set prognostic state variables at level midpoints
      elem%state%v    (:,:,1,:,n) = u
      elem%state%v    (:,:,2,:,n) = v
      elem%state%T    (:,:,:,  n) = T
      elem%state%ps_v (:,:,    n) = ps
      elem%state%phis (:,:)       = phis
    end do

  end subroutine set_elem_state

  !_____________________________________________________________________
  subroutine get_state(u,v,w,T,pnh,dp,ps,rho,zm,g,elem,hvcoord,nt,ntQ)

    ! get state variables at layer midpoints
    ! used by tests to compute idealized physics forcing terms

    real(real_kind), dimension(np,np,nlev), intent(inout) :: u,v,w,T,pnh,dp,zm,rho
    real(real_kind), dimension(np,np),      intent(inout) :: ps
    real(real_kind), intent(in)    :: g
    integer,         intent(in)    :: nt,ntQ
    type(element_t), intent(inout) :: elem
    type (hvcoord_t),intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct

    real(real_kind) , dimension(np,np,nlev) :: phi,p,kappa_star, Rstar

    integer :: k
    call  get_phi(elem,phi,hvcoord,nt,ntQ)
    
    ! set prognostic state variables at level midpoints
    u   = elem%state%v   (:,:,1,:,nt)
    v   = elem%state%v   (:,:,2,:,nt)
    w   = 0                             ! todo: w = -omega/(rho g)
    T   = elem%state%T   (:,:,:,nt)
    ps  = elem%state%ps_v(:,:,nt)
    zm  = phi/g

    do k=1,nlev
       p (:,:,k)= hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem%state%ps_v(:,:,nt)
       dp(:,:,k)=(hvcoord%hyai(k+1)-hvcoord%hyai(k))*hvcoord%ps0 +(hvcoord%hybi(k+1)-hvcoord%hybi(k))*elem%state%ps_v(:,:,nt)
    enddo

    Rstar = Rgas
    if (use_moisture) Rstar = (Rgas+(Rwater_vapor -Rgas)*elem%state%Qdp(:,:,:,1,ntQ)/dp)

    pnh = p
    rho = p/(Rstar*T)

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
  real(real_kind), intent(in)   :: zm(nlev)       ! height at layer midpoints
  real(real_kind), intent(in)   :: zi(nlevp)      ! height at interfaces
  real(real_kind), intent(in)   :: ztop           ! top of atm height
  real(real_kind), intent(in)   :: zc             ! cutoff height
  real(real_kind), intent(in)   :: tau            ! damping timescale
  real(real_kind), intent(in)   :: u0(np,np,nlev) ! reference u
  real(real_kind), intent(in)   :: v0(np,np,nlev) ! reference v
  integer,         intent(in)   :: n              ! timestep

  real(real_kind):: f_d(nlev)
  integer :: k

  ! Compute damping as a function of layer-midpoint height
  f_d=0.0d0
  where(zm .ge. zc);   f_d = sin(dd_pi/2 *(zm - zc)/(ztop - zc))**2; end where
  where(zm .ge. ztop); f_d = 1.0d0; end where
  f_d = -f_d/tau

  do k=1,nlev
     elem%derived%FM(:,:,1,k) = f_d(k) * ( elem%state%v(:,:,1,k,n) - u0(:,:,k) )
     elem%derived%FM(:,:,2,k) = f_d(k) * ( elem%state%v(:,:,2,k,n) - v0(:,:,k) )
  enddo
  end subroutine 

  !____________________________________________________________________
  subroutine tests_finalize(elem,hvcoord,ns,ne,ie)
  implicit none

  type(hvcoord_t),     intent(in)  :: hvcoord
  type(element_t),  intent(inout)  :: elem
  integer,             intent(in)  :: ns,ne
  integer, optional,   intent(in)   :: ie ! optional element index, to save initial state

  !do nothing
  end subroutine tests_finalize


end module

