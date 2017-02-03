!
! PESE Dynamics: Primtive Equations, Vertical Spectral Elements
!_______________________________________________________________________
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_advance_mod

  use bndry_mod,			 only: bndry_exchangev
  use control_mod,     only: qsplit, prescribed_wind, use_moisture
  use derivative_mod,  only: derivative_t,  divergence_sphere, gradient_sphere
  use dimensions_mod,  only: np, nlev, nlevp, qsize
  use edge_mod,				 only: initEdgeBuffer
  use edgetype_mod,    only: edgebuffer_t
  use element_mod,		 only: element_t
  use element_state,	 only: elem_state_t, derived_state_t
  use element_ops,     only: pack_edge_data, unpack_edge_data, apply_map, apply_vertical_dss, display_max_and_min
  use hybrid_mod,			 only: hybrid_t
  use hybvcoord_mod,	 only: hvcoord_t
	use kinds,					 only: rl => real_kind, real_kind
  use parallel_mod,    only: abortmp, parallel_t
  use physical_constants, only : cp, cpwater_vapor, Rgas, kappa, p0, g, Rwater_vapor
  use time_mod,        only: timeLevel_t, timelevel_update, timelevel_qdp, nsplit
  use vertical_se,     only: npv, solver_args, pack_solver_args, eta_derivative, advection2,&
                             self_advection2, eta_integral_from_n,eta_integral_from_1, make_vertical_mesh
  use viscosity_mod,    only: apply_hyperviscosity

#ifndef CAM
  use test_mod,        only: set_test_prescribed_wind
#endif

  implicit none

  integer,parameter   :: n_rhs    = 1 + 3*nlev  ! ps + T + u + v          ! num levels in rhs and edge-buffer
  logical,parameter   :: verbose  = .false.                               ! verbose output flag
  integer             :: count    = 1
  type (edgebuffer_t) :: edge_buffer

contains

  !_____________________________________________________________________
  subroutine prim_advance_init(par, elem, integration)
    implicit none
    
    type (parallel_t),    intent(in)              :: par
    type (element_t),     intent(inout), target   :: elem(:)
    character(len=*)    , intent(in)              :: integration
    integer :: i, ie

    ! allocate space for edge-data exchange
    call initEdgeBuffer(par,edge_buffer,elem,n_rhs)

  end subroutine prim_advance_init

  !_____________________________________________________________________
	subroutine prim_advance_init2(elem, nets, nete, hybrid, hvcoord)

    type (element_t),			intent(inout), target :: elem(:)							! array of element_t structures
    integer,							intent(in)		:: nets,nete										! start and end element indices
    type (hybrid_t),			intent(in)		:: hybrid												! mpi/omp data struct
    type (hvcoord_t),			intent(inout)	:: hvcoord											! hybrid vertical coord data struct

    if (hybrid%masterthread) print *,"initializing PESE dynamics solver"

    ! initialize vertical operators and coordinates
    call make_vertical_mesh(hybrid, hvcoord)

	end subroutine

  !_____________________________________________________________________
  subroutine prim_advance_exp(elem,deriv,hvcoord,hybrid,dt,tl,nets,nete,compute_diagnostics)

    type (element_t),   intent(inout), target :: elem(:)
    type (derivative_t),intent(in)            :: deriv
    type (hvcoord_t)                          :: hvcoord
    type (hybrid_t),    intent(in)            :: hybrid
    real (rl),          intent(in)            :: dt
    type (TimeLevel_t), intent(in)            :: tl
    integer,            intent(in)            :: nets, nete
    logical,            intent(in)            :: compute_diagnostics

    integer :: ie, t, q,k,i,j,n, qn0
    integer :: nm1,n0,np1,nstep
    real (rl) ::  eta_ave_w

    qn0   = 1

    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1

    nstep = tl%nstep
    eta_ave_w = 1d0/qsplit  ! set q subcycle time averaging weight

#ifndef CAM
    ! if test uses prescribed wind, set dynamic variables analytically and return
    if (prescribed_wind ==1 ) then
      call set_test_prescribed_wind(elem,deriv,hybrid,hvcoord,dt,tl,nets,nete)
      return
    endif
#endif

    ! integrate dynamics in time using RK3

    call compute_and_apply_rhs(np1,n0,n0,qn0,dt/3,elem,hvcoord,hybrid,deriv,nets,nete,compute_diagnostics,0d0)
    ! u2 = u0 + dt/2 RHS(u1)
    call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,deriv,nets,nete,.false.,0d0)
    ! u3 = u0 + dt RHS(u2)
    call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,deriv,nets,nete,.false.,eta_ave_w)

    ! apply viscosity or hyperviscosity to u,v,T,ps
    call apply_hyperviscosity(elem,edge_buffer,hvcoord,hybrid,deriv,np1,nets,nete,dt,eta_ave_w)

  end subroutine

	!_____________________________________________________________________
  subroutine compute_and_apply_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,deriv,nets,nete,compute_diagnostics,eta_ave_w)

		integer,							intent(in)		:: np1,nm1,n0										! time indices
		integer,							intent(in)		:: qn0													! time level used for virtual temperature
    real*8,               intent(in)    :: dt2
		type (element_t),			intent(inout), target :: elem(:)							! array of element_t structures
		type (hvcoord_t),			intent(inout)	:: hvcoord											! hybrid vertical coord data struct
		type (hybrid_t),			intent(in)		:: hybrid												! mpi/omp data struct
		integer,							intent(in)		:: nets,nete										! start and end element indices
		type (derivative_t),	intent(in)		:: deriv												! horizontal derivative data struct
    logical,              intent(in)    :: compute_diagnostics          ! TODO: enable energy diagnostics
    real (rl),            intent(in)    :: eta_ave_w                    ! TODO: enable qsplit

		type(solver_args) :: a																							! solver arguments
    integer  :: i,k,ie

    count=count+1
    a = pack_solver_args(np1,nm1,n0,qn0,dt2,hvcoord,hybrid,deriv,nets,nete,compute_diagnostics,eta_ave_w)

    do ie=nets,nete
      if (verbose) call display_max_and_min(elem(ie),hybrid,ie,n0,count)
      call apply_rhs(elem(ie),ie, a)
      call apply_vertical_dss(elem(ie),np1)                             ! perform vertical direct stiffness summation
      call apply_map(elem(ie)%spheremp,elem(ie),np1)                    ! apply sphere-map
      call pack_edge_data(edge_buffer,elem(ie),np1,ie)                  ! pack edge data in buffer
    enddo

    call bndry_exchangeV(hybrid, edge_buffer)                           ! exchange edge-data, perform dss

    do ie=nets,nete
      call unpack_edge_data(edge_buffer,elem(ie),np1,ie)                ! unpack edge data from buffer
      call apply_map(elem(ie)%rspheremp, elem(ie), np1)                 ! apply inverse sphere-map
    enddo

 	end subroutine

 	!_____________________________________________________________________
	subroutine apply_rhs(e,ie,a)

		type (element_t), intent(inout), target :: e                        ! element to operate upon
    integer,          intent(in)            :: ie                       ! element index
		type(solver_args),intent(inout)         :: a                        ! solver arguments

		real (rl), dimension(np,np,nlev)	:: u,v,T,p,phi,Tv,qv,alpha        ! local copies of 3d fields
		real (rl), dimension(np,np)				:: ps                             ! local copies of 2d fields
		real (rl), dimension(np,np,nlev)	:: dp_deta, dT_deta								! vertical hydrostatic pressure gradient
    real (rl), dimension(np,np,nlev)	:: du_deta, dv_deta               ! vertical velocity gradients
    real (rl), dimension(np,np,2,nlev):: hflux													! horizontal mass flux
    real (rl), dimension(np,np)	      :: zeros = 0                      ! square matrix of zeros
		real (rl), dimension(np,np,2,nlev):: grad_p                         ! gradient of total pressure
		real (rl), dimension(np,np,2,nlev):: grad_phi												! gradient of geopotential height
		real (rl), dimension(np,np,nlev)	:: eta_flux												! vertical mass flux
    real (rl), dimension(np,np,nlev)	:: eta_dot												! vertical mass flux
		real (rl), dimension(np,np,nlev)	:: div_hflux											! flux of mass out of the column
		real (rl), dimension(np,np,nlev)	:: column_flux										! total mass flux into overlying column
		real (rl), dimension(np,np,nlev)	:: omega													! total time deriv of hydrostatic-pressure
		real (rl), dimension(np,np,2,nlev):: f_corilois											! coriolis force
    real (rl), dimension(np,np,nlev)	:: ddt_T,ddt_u,ddt_v              ! total time derivatives
    real (rl), dimension(np,np,nlev)	:: adv_T                          ! advective terms
    real (rl), dimension(np,np,2,nlev):: adv_v                          ! advective terms 2d
		real (rl), dimension(np,np,nlev)	:: rhs_T,rhs_u,rhs_v              ! partial time derivatives
    real (rl), dimension(np,np)       :: rhs_ps
    real (rl), dimension(np,np,nlev)  :: CpStar
    type (elem_state_t), pointer      :: s															! pointer to element state variables

    integer :: i,j,k                                                    ! loop indices

    ! make local copies of variables for readability
    s   => e%state
    u   = s%v   (:,:,1,:, a%n0)
    v   = s%v   (:,:,2,:, a%n0)
    T   = s%T   (:,:,:,   a%n0)
    ps  = s%ps_v(:,:,     a%n0)
    qv  = s%Q   (:,:,:,1)       ! water vapor mixing ratio

    ! get hydrostatic-pressure at eta levels
    forall(k=1:nlev) p(:,:,k)= a%hvcoord%hyam(k)*p0 + a%hvcoord%hybm(k)*ps

    ! get gradient of hydrostatic pressure
    do k=1,nlev; grad_p (:,:,:,k) = gradient_sphere(p (:,:,k), a%deriv, e%Dinv); enddo

    if (.not. use_moisture ) then
      Tv      = T
      CpStar  = Cp
    else
      ! get virtual temperature (TODO)
      Tv      = T *(1.0_rl + (Rwater_vapor/Rgas - 1.0_rl)*qv)
      CpStar  = Cp*(1.0_rl + (Cpwater_vapor/Cp  - 1.0_rl)*qv)
    endif

    ! get moist density from virtual temperature and pressure
    alpha = rgas*Tv/p

    ! get vertical derivatives
    dp_deta = eta_derivative(p)
    du_deta = eta_derivative(u)
    dv_deta = eta_derivative(v)
    dT_deta = eta_derivative(T)

    ! get geopotential-height by integrating hydrostatic balance from bottom
    phi = eta_integral_from_n( -dp_deta*alpha, s%phis )

    ! get horizontal gradient of geopotential-height
    do k=1,nlev; grad_phi(:,:,:,k) = gradient_sphere(phi(:,:,k),a%deriv,e%Dinv); enddo

    ! get horizontal mass-flux
    hflux(:,:,1,:) = dp_deta*u
    hflux(:,:,2,:) = dp_deta*v

    ! get divergence of horizontal mass flux
    do k=1,nlev;
      div_hflux(:,:,k)= divergence_sphere(hflux(:,:,:,k), a%deriv, e);
    enddo

    ! get net mass-flux out of column, integrated from the top
    column_flux = eta_integral_from_1( div_hflux, zeros)

    ! get vertical mass-flux
    forall(k=1:nlev) eta_flux(:,:,k) = a%hvcoord%hybm(k)*column_flux(:,:,nlev)-column_flux(:,:,k)
    eta_flux(:,:,1   ) = 0.0d0  ! enforce no flux at top
    eta_flux(:,:,nlev) = 0.0d0  ! enforce no flux at bottom

    ! get vertical velocity from vertical flux
    eta_dot = eta_flux/dp_deta

    ! get time derivative of hydrostatic-pressure
    omega = u*grad_p(:,:,1,:)+v*grad_p(:,:,2,:)-column_flux(:,:,:)

    ! get coriolis force
    do k=1,nlev
      f_corilois(:,:,1,k)= -e%fcor * v(:,:,k)
      f_corilois(:,:,2,k)= +e%fcor * u(:,:,k)
    enddo

    ! get total time-derivs
    ddt_u  = -grad_p(:,:,1,:)*alpha -grad_phi(:,:,1,:) -f_corilois(:,:,1,:)
    ddt_v  = -grad_p(:,:,2,:)*alpha -grad_phi(:,:,2,:) -f_corilois(:,:,2,:)
    ddt_T  = omega*alpha/CpStar

    ! get advection terms
    adv_T  = advection2(T,dT_deta,u,v,eta_dot,a,e)
    adv_v  = self_advection2(u,v,du_deta,dv_deta,eta_dot,a,e)

    ! get partial-time derivative of each prognostic variable
    rhs_ps = -column_flux(:,:,nlev)
    rhs_T  = -adv_T          + ddt_T
    rhs_u  = -adv_v(:,:,1,:) + ddt_u
    rhs_v  = -adv_v(:,:,2,:) + ddt_v

    ! apply right-hand-side to prognostics
    s%ps_v(:,:,    a%np1) = s%ps_v(:,:    ,a%nm1) + a%dt * rhs_ps
    s%T   (:,:,:,  a%np1) = s%T   (:,:,:  ,a%nm1) + a%dt * rhs_T
    s%v   (:,:,1,:,a%np1) = s%v   (:,:,1,:,a%nm1) + a%dt * rhs_u
    s%v   (:,:,2,:,a%np1) = s%v   (:,:,2,:,a%nm1) + a%dt * rhs_v

    ! compute derived quantities for output
    !s%dp3d              = 1
    !e%derived%omega_p   = omega/p
    e%derived%phi       = phi
    e%derived%omega     = omega

	end subroutine

  !_____________________________________________________________________
  subroutine applyCAMforcing(elem,hvcoord,np1,np1_qdp,dt_q,nets,nete)

    use physical_constants, only: Cp

  !  implicit none
    type (element_t),       intent(inout) :: elem(:)
    real (kind=real_kind),  intent(in)    :: dt_q
    type (hvcoord_t),       intent(in)    :: hvcoord
    integer,                intent(in)    :: np1,nets,nete,np1_qdp

    ! local
    integer  :: ie,q
    real(rl) :: dp(np,np,nlev)

    do ie=nets,nete

      ! apply forcing to dynamics
      elem(ie)%state%T(:,:,:,  np1) = elem(ie)%state%T(:,:,:,np1)   + dt_q*elem(ie)%derived%FT
      elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + dt_q*elem(ie)%derived%FM

! TODO: enable tracer forcing for PESE

      ! apply forcing to tracer mass
      !dp = elem(ie)%state%dp3d(:,:,:,np1)
      !do q=1,qsize
      !  elem(ie)%state%Qdp(:,:,:,q,np1_qdp) = elem(ie)%state%Qdp(:,:,:,q,np1_qdp)+ dt_q * elem(ie)%derived%FQ(:,:,:,q)
      !  elem(ie)%state%Q(:,:,:,q)           = elem(ie)%state%Qdp(:,:,:,q,np1_qdp)/dp
      !enddo
    enddo
    end subroutine applyCAMforcing

    !___________________________________________________________________
    subroutine applyCAMforcing_dynamics(elem,hvcoord,np1,np1_qdp,dt_q,nets,nete)

    use hybvcoord_mod,  only: hvcoord_t

    implicit none
    type (element_t)     ,  intent(inout) :: elem(:)
    real (kind=real_kind),  intent(in)    :: dt_q
    type (hvcoord_t),       intent(in)    :: hvcoord
    integer,                intent(in)    :: np1,np1_qdp,nets,nete

    integer :: ie,q
    real (kind=real_kind) :: v1,dp

    do ie=nets,nete
       elem(ie)%state%T(:,:,:,np1)   = elem(ie)%state%T(:,:,:,np1)   + dt_q*elem(ie)%derived%FT(:,:,:)
       elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + dt_q*elem(ie)%derived%FM(:,:,:,:)
    enddo

    end subroutine applyCAMforcing_dynamics

end module prim_advance_mod

