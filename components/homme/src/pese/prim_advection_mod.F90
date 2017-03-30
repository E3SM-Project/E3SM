#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_advection_mod

  use prim_advection_base, only: prim_advec_init1_rk2, prim_advec_tracers_remap_rk2

  use bndry_mod,			 only: bndry_exchangev
  use derivative_mod,  only: derivative_t, gradient_sphere
	use dimensions_mod,  only: np, nlev, nlevp, qsize,qsize_d, nelemd
  use edge_mod,				 only: initEdgeBuffer, edgevpack, edgevunpack
  use edgetype_mod,    only: edgebuffer_t
  use element_mod,		 only: element_t
  use element_state,	 only: elem_state_t, derived_state_t
	use hybrid_mod,			 only: hybrid_t
	use hybvcoord_mod,	 only: hvcoord_t
  use parallel_mod,    only: parallel_t
  use perf_mod,        only: t_startf, t_stopf
  use physical_constants, only : p0
  use kinds,					 only: rl => real_kind, real_kind, dd => longdouble_kind
  use time_mod,        only: timelevel_t, timelevel_update, timelevel_qdp, nsplit
  use vertical_se,     only: eta_derivative, vertical_dss

	implicit none

  type(EdgeBuffer_t) :: edge_buffer                                     ! buffer for edge-data exchange

  ! q(np,np,nlev,qsize,nstage,nelemd)
  real(rl), allocatable :: q(:,:,:,:,:,:)                               ! tracer sub-stage storage

  ! rhs(np,np,nlev,qsize,nelemd)
  real(rl), allocatable :: rhs(:,:,:,:,:)                               ! right-hand-side storage

  type (derivative_t) :: deriv
  type (hybrid_t)     :: hybrid                                         ! mpi/omp data struct
  type (hvcoord_t)    :: hv                                             ! hybrid vertical coord data struct
  integer             :: nets,nete                                      ! start and end element indices

  CONTAINS

  !_____________________________________________________________________
  subroutine prim_advect_tracers(elem, deriv_, hybrid_, nets_, nete_, hvcoord_, dt, tl)

    type (element_t),			intent(inout), target :: elem(:)							! array of element_t structures
    type (derivative_t),  intent(in)    :: deriv_
    type (hybrid_t),			intent(in)		:: hybrid_                      ! mpi/omp data struct
    integer,							intent(in)		:: nets_,nete_									! start and end element indices
    type (hvcoord_t),			intent(in)    :: hvcoord_											! hybrid vertical coord data struct
    real(rl),             intent(in)    :: dt                           ! tracer timestep
    type (TimeLevel_t),   intent(inout) :: tl                           ! time-level

    real(rl) :: dp(np,np,nlev), dp_dn

    logical :: initialized = .false.
    integer :: ie,iq,k
    integer :: n0_qdp, np1_qdp

    call t_startf("advect_tracers")

    ! store arguments at module level
    deriv = deriv_; hybrid=hybrid_; nets=nets_; nete=nete_; hv=hvcoord_

    ! allocate edge-buff and rhs-vector based on qsize
    if(.not. initialized) then; initialized = .true.
      call initEdgeBuffer(hybrid%par, edge_buffer, elem, qsize*nlev)    ! allocate edge-data buffer
      allocate( rhs(np,np,nlev,qsize,nelemd))                           ! allocate memmory for euler sub-steps
    endif

    ! advance using 3 stage SSP Runge-Kutta from [Spiteri and Ruuth 2002]
    call RK_SSP3(elem,dt,tl)
    !call euler_forward(elem,dt,tl)

    ! compute Qdp (it is used by many routines in preqx)
    do ie=nets,nete

      dp = elem(ie)%state%dp3d(:,:,:,tl%n0)

      do iq=1,qsize
        elem(ie)%state%Qdp(:,:,:,iq,1) = elem(ie)%state%Q(:,:,:,iq)*dp
        elem(ie)%state%Qdp(:,:,:,iq,2) = elem(ie)%state%Q(:,:,:,iq)*dp
      enddo
    enddo

    call t_stopf("advect_tracers")

  end subroutine

  !_____________________________________________________________________
  subroutine euler_forward(elem,dt,tl)

    ! simplest time-stepping solver: take one explicit step forward

    type (element_t),			intent(inout), target :: elem(:)							! array of element_t structures
    real*8,               intent(in)    :: dt                           ! tracer timestep
    type (TimeLevel_t),   intent(in)    :: tl                           ! time-level
    logical :: initialized = .false.
    integer :: ie

    ! allocate space to store sub-stages
    if(.not. initialized) then; initialized=.true.
      allocate( q(np,np,nlev,qsize,1,nelemd) ); q=0.0d0
    endif

    do ie=nets,nete
      q(:,:,:,:,1,ie) = elem(ie)%state%Q(:,:,:,1:qsize)
    enddo

    call euler_step(elem, 1,1,1 ,dt, tl)        ! q1 = q1 + dt RHS[q1]

    do ie=nets,nete                             ! q  = q1
      elem(ie)%state%Q(:,:,:,1:qsize) = 1.0*q(:,:,:,:,1,ie)
    enddo

  end subroutine
  !_____________________________________________________________________
  subroutine RK3(elem,dt,tl)

    type (element_t),			intent(inout), target :: elem(:)							! array of element_t structures
    real*8,               intent(in)    :: dt                           ! tracer timestep
    type (TimeLevel_t),   intent(in)    :: tl                           ! time-level
    logical :: initialized = .false.
    integer :: ie

    ! allocate space to store sub-stages
    if(.not. initialized) then; initialized=.true.
      allocate( q(np,np,nlev,qsize,4,nelemd) ); q=0
    endif

    ! copy tracer data into sub-stage 1
    do ie=nets,nete
      q(:,:,:,:,1,ie) = elem(ie)%state%Q(:,:,:,1:qsize)
    enddo

    ! RK 3 stage
    call euler_step(elem, 1,1,2 ,dt/3.0d0, tl)  ! q2 = q1 + dt/3 RHS[q1]
    call euler_step(elem, 1,2,3 ,dt/2.0d0, tl)  ! q3 = q1 + dt/2 RHS[q2]
    call euler_step(elem, 1,3,4, dt/1.0d0, tl)  ! q4 = q1 + dt   RHS[q3]

    do ie=nets,nete                             ! q  = q4
      elem(ie)%state%Q(:,:,:,1:qsize) = 1.0*q(:,:,:,:,4,ie)
    enddo

  end subroutine

  !_____________________________________________________________________
  subroutine RK_SSP3(elem,dt,tl)

    type (element_t),			intent(inout), target :: elem(:)							! array of element_t structures
    real*8,               intent(in)    :: dt                           ! tracer timestep
    type (TimeLevel_t),   intent(in)    :: tl                           ! time-level
    logical :: initialized = .false.
    integer :: ie

    ! allocate space to store sub-stages
    if(.not. initialized) then; initialized=.true.
      allocate( q(np,np,nlev,qsize,2,nelemd) ); q=0
    endif

    ! copy tracer data into sub-stage 1
    do ie=nets,nete
      q(:,:,:,:,1,ie) = elem(ie)%state%Q(:,:,:,1:qsize)
    enddo

    ! RK-SSP 3 stage from Spiteri and Ruuth 2002
    call euler_step(elem, 1,1,2 ,dt/2.0d0, tl) ! q2 = q1 + dt/2 RHS[q1]
    call euler_step(elem, 2,2,2 ,dt/2.0d0, tl) ! q2 = q2 + dt/2 RHS[q2]
    call euler_step(elem, 2,2,2 ,dt/2.0d0, tl) ! q2 = q2 + dt/2 RHS[q2]

    do ie=nets,nete;                           ! q  = 1/3*q + 2/3 q2
      elem(ie)%state%Q(:,:,:,1:qsize) = ( q(:,:,:,:,1,ie)+2.0d0*q(:,:,:,:,2,ie) )/3.0d0
    enddo

  end subroutine

	!_____________________________________________________________________
  subroutine euler_step(elem,n0,n1,n2,dt,tl)

    ! solve dq/dt = -u dot grad q -eta_dot dq/dn 

    type (element_t),			intent(inout), target :: elem(:)							! element array
		integer,							intent(in)		:: n0,n1,n2                     ! substep indices: q(n2) = q(n0) + dt * RHS(q(n1))
    real*8,               intent(in)    :: dt                           ! euler forward timestep
    type (TimeLevel_t),   intent(in)    :: tl                           ! time-level

    integer :: k,ie,iq                                                  ! loop indices

    real(rl) :: rhs(np,np,nlev,qsize)
    real(rl), dimension(np,np)       :: ps
    real(rl), dimension(np,np,nlev)  :: p, eta_dot, vflux
    real(rl), dimension(np,np,nlev)  :: dp_dn, div_hflux, ddn_vflux
    real(rl), dimension(np,np,2,nlev):: hflux, V
    real(rl), dimension(np,np,2,nlev):: grad_q
    real(rl), dimension(np,np,nlev)  :: dq_dn
    real(rl), dimension(np,np,nlev)  :: h_adv, v_adv

    do ie=nets,nete

      V   = elem(ie)%state%v(:,:,:,:,tl%np1)
      ps  = elem(ie)%state%ps_v(:,:,tl%np1)

      ! get eta_dot from vertical mass flux
      forall(k=1:nlev) p(:,:,k)= hv%hyam(k)*p0 + hv%hybm(k)*ps
      dp_dn   = eta_derivative(p)
      eta_dot = elem(ie)%derived%eta_dot_dpdn(:,:,2:nlevp) / dp_dn

      ! ensure no flux at top and bottom
      eta_dot(:,:,1)    =0.0
      eta_dot(:,:,nlev) =0.0

      do iq=1,qsize

        ! get horizontal advection term
        do k=1,nlev; grad_q(:,:,:,k) = gradient_sphere( q(:,:,k,iq,n1,ie),deriv, elem(ie)%Dinv); enddo
        h_adv = V(:,:,1,:)*grad_q(:,:,1,:)+V(:,:,2,:)*grad_q(:,:,2,:)

        ! get vertical advection term
        dq_dn = eta_derivative( q(:,:,:,iq,n1,ie) )
        v_adv = eta_dot*dq_dn

        rhs(:,:,:,iq) =-h_adv -v_adv
      enddo

      ! compute euler-forward step
      q(:,:,:,:,n2,ie) = q(:,:,:,:,n0,ie) + dt * rhs

      do iq=1,qsize

        ! apply vertical dss
        call vertical_dss(q(:,:,:,iq,n2,ie))

        do k=1,nlev
          ! apply mass matrix
          q(:,:,k,iq,n2,ie) = q(:,:,k,iq,n2,ie)*elem(ie)%spheremp
        enddo
      enddo

      ! pack data into edge buffer
      call edgeVpack(edge_buffer, q(:,:,:,:,n2,ie), qsize*nlev, 0, ie);
    enddo

    ! exchange edge data and perform direct stiffness summation
    call bndry_exchangeV(hybrid, edge_buffer)

    do ie=nets,nete

      ! unpack data from edge buffer
      call edgeVunpack(edge_buffer, q(:,:,:,:,n2,ie), qsize*nlev, 0, ie)

      do iq=1,qsize
        do k=1,nlev
          ! apply inverse mass matrix
          q(:,:,k,iq,n2,ie) = q(:,:,k,iq,n2,ie)*elem(ie)%rspheremp
        enddo
      enddo
    enddo

	end subroutine

  !_____________________________________________________________________
  subroutine prim_advec_init1(par, elem, n_domains)

    type(parallel_t)    :: par
    integer, intent(in) :: n_domains
    type (element_t)    :: elem(:)

    call prim_advec_init1_rk2(par, elem, n_domains)

  end subroutine prim_advec_init1

  !_____________________________________________________________________
  subroutine prim_advec_init2(elem,hvcoord,hybrid)

    use element_mod   , only : element_t
    use hybvcoord_mod , only : hvcoord_t
    implicit none
    type(element_t)   , intent(in) :: elem(:)
    type(hvcoord_t)   , intent(in) :: hvcoord
    type (hybrid_t)   , intent(in) :: hybrid
  end subroutine prim_advec_init2


end module prim_advection_mod
