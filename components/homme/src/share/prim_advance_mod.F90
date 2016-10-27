#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!#define _DBG_ print *,"File:",__FILE__," at ",__LINE__
!#define _DBG_ !DBG
!
!
module prim_advance_mod

  use control_mod,    only: qsplit,rsplit
  use derivative_mod, only: derivative_t, vorticity, divergence, gradient, gradient_wk
  use dimensions_mod, only: np, nlev, nlevp, nvar, nc, nelemd
  use edgetype_mod,   only: EdgeDescriptor_t, EdgeBuffer_t
  use element_mod,    only: element_t
  use hybrid_mod,     only: hybrid_t
  use hybvcoord_mod,  only: hvcoord_t
  use kinds,          only: real_kind, iulog
  use perf_mod,       only: t_startf, t_stopf, t_barrierf, t_adj_detailf ! _EXTERNAL
  use parallel_mod,   only: abortmp, parallel_t, iam
  use time_mod,       only: timelevel_t

  implicit none
  private
  save
  public :: prim_advance_exp, prim_advance_si, prim_advance_init, preq_robert3,&
       applyCAMforcing_dynamics, applyCAMforcing, smooth_phis, overwrite_SEdensity

#ifdef TRILINOS
  public :: distribute_flux_at_corners
#endif
  
  type (EdgeBuffer_t) :: edge1
  type (EdgeBuffer_t) :: edge2
  type (EdgeBuffer_t) :: edge3p1

  real (kind=real_kind) :: initialized_for_dt   = 0

  real (kind=real_kind), allocatable :: ur_weights(:)

contains

  subroutine prim_advance_init(par, elem,integration)
    use edge_mod, only : initEdgeBuffer
    implicit none
    
    type (parallel_t) :: par
    type (element_t), intent(inout), target   :: elem(:)
    character(len=*)    , intent(in) :: integration
    integer :: i
    integer :: ie

    call initEdgeBuffer(par,edge3p1,elem,4*nlev)

    if(integration == 'semi_imp') then
       call initEdgeBuffer(par,edge1,elem,nlev)
       call initEdgeBuffer(par,edge2,elem,2*nlev)
    end if

    ! compute averaging weights for RK+LF (tstep_type=1) timestepping:
    allocate(ur_weights(qsplit))
    ur_weights(:)=0.0d0

    if(mod(qsplit,2).NE.0)then
       ur_weights(1)=1.0d0/qsplit
       do i=3,qsplit,2
         ur_weights(i)=2.0d0/qsplit
       enddo
    else
       do i=2,qsplit,2
         ur_weights(i)=2.0d0/qsplit
       enddo
    endif

  end subroutine prim_advance_init

#ifndef CAM
  !_____________________________________________________________________
  subroutine set_prescribed_wind(elem,deriv,hybrid,hv,dt,tl,nets,nete,eta_ave_w)

    use test_mod,  only: set_test_prescribed_wind
    use asp_tests, only: asp_advection_vertical

    type (element_t),      intent(inout), target  :: elem(:)
    type (derivative_t),   intent(in)             :: deriv
    type (hvcoord_t),      intent(inout)          :: hv
    type (hybrid_t),       intent(in)             :: hybrid
    real (kind=real_kind), intent(in)             :: dt
    type (TimeLevel_t)   , intent(in)             :: tl
    integer              , intent(in)             :: nets
    integer              , intent(in)             :: nete
    real (kind=real_kind), intent(in)             :: eta_ave_w

    real (kind=real_kind) :: dp(np,np)! pressure thickness, vflux
    real(kind=real_kind)  :: time

    integer :: ie,k,n0,np1

    time  = tl%nstep*dt
    n0    = tl%n0
    np1   = tl%np1

    call set_test_prescribed_wind(elem,deriv,hybrid,hv,dt,tl,nets,nete)

    do ie = nets,nete
      ! asp2008 tests:
      ! call asp_advection_vertical(time,hv,elem(ie)%state%ps_v(:,:,n0),eta_dot_dpdn)
      ! accumulate mean fluxes for advection
      !   if (rsplit==0) then
      !      elem(ie)%derived%eta_dot_dpdn(:,:,:) = &
      !           elem(ie)%derived%eta_dot_dpdn(:,:,:) + eta_dot_dpdn(:,:,:)*eta_ave_w
      !   else
      !      ! lagrangian case.  mean vertical velocity = 0. compute dp3d on floating levels
      !      elem(ie)%derived%eta_dot_dpdn(:,:,:) = 0
      !     do k=1,nlev
      !         elem(ie)%state%dp3d(:,:,k,np1) = elem(ie)%state%dp3d(:,:,k,n0) + dt*(eta_dot_dpdn(:,:,k+1) - eta_dot_dpdn(:,:,k))
      !      enddo
      !   end if

      ! get mean horizontal flux (rho*vel) for tracer advection
      do k=1,nlev
         if (rsplit==0) then
            dp(:,:) =(hv%hyai(k+1)-hv%hyai(k))*hv%ps0 + (hv%hybi(k+1)-hv%hybi(k))*elem(ie)%state%ps_v(:,:,n0)
         else
            dp(:,:) = elem(ie)%state%dp3d(:,:,k,n0)
         end if
         elem(ie)%derived%vn0(:,:,1,k)=elem(ie)%derived%vn0(:,:,1,k) + eta_ave_w*elem(ie)%state%v(:,:,1,k,n0)*dp(:,:)
         elem(ie)%derived%vn0(:,:,2,k)=elem(ie)%derived%vn0(:,:,2,k) + eta_ave_w*elem(ie)%state%v(:,:,2,k,n0)*dp(:,:)
      enddo
   end do

  end subroutine
#endif

  !_____________________________________________________________________
  subroutine prim_advance_exp(elem, deriv, hvcoord, hybrid,dt, tl,  nets, nete, compute_diagnostics)

    use bndry_mod,      only: bndry_exchangev
    use control_mod,    only: prescribed_wind, qsplit, tstep_type, rsplit, qsplit, moisture, integration
    use edge_mod,       only: edgevpack, edgevunpack, initEdgeBuffer
    use edgetype_mod,   only: EdgeBuffer_t
    use reduction_mod,  only: reductionbuffer_ordered_1d_t
    use time_mod,       only: timelevel_qdp, tevolve
    use diffusion_mod,  only: prim_diffusion

#ifdef TRILINOS
    use prim_derived_type_mod ,only : derived_type, initialize
    use, intrinsic :: iso_c_binding
#endif

#ifdef CAM
    use control_mod,    only: prescribed_vertwind
#else
    use asp_tests,      only: asp_advection_vertical
#endif

    implicit none

    type (element_t),      intent(inout), target :: elem(:)
    type (derivative_t),   intent(in)            :: deriv
    type (hvcoord_t)                             :: hvcoord
    type (hybrid_t),       intent(in)            :: hybrid
    real (kind=real_kind), intent(in)            :: dt
    type (TimeLevel_t)   , intent(in)            :: tl
    integer              , intent(in)            :: nets
    integer              , intent(in)            :: nete
    logical,               intent(in)            :: compute_diagnostics

    real (kind=real_kind) ::  dt2, time, dt_vis, x, eta_ave_w
    real (kind=real_kind) ::  tempdp3d(np,np)
    real (kind=real_kind) ::  tempmass(nc,nc)
    real (kind=real_kind) ::  tempflux(nc,nc,4)
    !real (kind=real_kind) ::  dn

    integer :: ie,nm1,n0,np1,nstep,method,qsplit_stage,k, qn0
    integer :: n,i,j,lx,lenx

#ifdef TRILINOS
    real (c_double) ,allocatable, dimension(:) :: xstate(:)

    ! state_object is a derived data type passed thru noxinit as a pointer
    type(derived_type) ,target         :: state_object
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    type(derived_type) ,target         :: pre_object
    type(derived_type) ,pointer         :: pptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_pre
    type(derived_type) ,target         :: jac_object
    type(derived_type) ,pointer         :: jptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_jac

    integer(c_int) :: ierr = 0

  interface

   subroutine noxsolve(vectorSize,vector,v_container,p_container,j_container,ierr) &
     bind(C,name='noxsolve')
     use ,intrinsic :: iso_c_binding
       integer(c_int)                :: vectorSize
       real(c_double)  ,dimension(*) :: vector
       type(c_ptr)                   :: v_container
       type(c_ptr)                   :: p_container  !precon ptr
       type(c_ptr)                   :: j_container  !analytic jacobian ptr
       integer(c_int)                :: ierr         !error flag
    end subroutine noxsolve

  end interface
#endif

    call t_startf('prim_advance_exp')
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

    ! get timelevel for accessing tracer mass Qdp() to compute virtual temperature

    qn0 = -1    ! -1 = disabled (assume dry dynamics)
    if ( moisture /= "dry") then
       call TimeLevel_Qdp(tl, qsplit, qn0)  ! compute current Qdp() timelevel
    endif

    ! get time integration method for this timestep

    method    = tstep_type  ! set default integration method
    eta_ave_w = 1d0/qsplit  ! set default vertical flux weight

    if(tstep_type==0)then
      ! 0: use leapfrog, but RK2 on first step
      if (nstep==0) method=1

    else if (tstep_type==1) then
      ! 1: use leapfrog, but RK2 on first qsplit stage
      method = 0
      qsplit_stage = mod(nstep,qsplit)    ! get qsplit stage
      if (qsplit_stage==0) method=1       ! use RK2 on first stage
      eta_ave_w=ur_weights(qsplit_stage+1)! RK2 + LF scheme has tricky weights
    endif

#ifndef CAM
    ! if "prescribed wind" set dynamics explicitly and skip time-integration
    if (prescribed_wind ==1 ) then
       call set_prescribed_wind(elem,deriv,hybrid,hvcoord,dt,tl,nets,nete,eta_ave_w)
       call t_stopf('prim_advance_exp')
       return
    endif
#endif

    ! integration = "explicit"
    !
    !   tstep_type=0  pure leapfrog except for very first timestep   CFL=1
    !                    typically requires qsplit=4 or 5
    !   tstep_type=1  RK2 followed by qsplit-1 leapfrog steps        CFL=close to qsplit
    !                    typically requires qsplit=4 or 5
    !   tstep_type=2  RK2-SSP 3 stage (as used by tracers)           CFL=.58
    !                    optimal in terms of SSP CFL, but not        CFLSSP=2
    !                    optimal in terms of CFL
    !                    typically requires qsplit=3
    !                    but if windspeed > 340m/s, could use this
    !                    with qsplit=1
    !   tstep_type=3  classic RK3                                    CFL=1.73 (sqrt(3))
    !
    !   tstep_type=4  Kinnmark&Gray RK4 4 stage                      CFL=sqrt(8)=2.8
    !                 should we replace by standard RK4 (CFL=sqrt(8))?
    !                 (K&G 1st order method has CFL=3)
    !   tstep_type=5  Kinnmark&Gray RK3 5 stage 3rd order            CFL=3.87  (sqrt(15))
    !                 From Paul Ullrich.  3rd order for nonlinear terms also
    !                 K&G method is only 3rd order for linear
    !                 optimal: for windspeeds ~120m/s,gravity: 340m/2
    !                 run with qsplit=1
    !                 (K&G 2nd order method has CFL=4. tiny CFL improvement not worth 2nd order)
    !
    ! integration = "full_imp"
    !
    !   tstep_type=1  Backward Euler or BDF2 implicit dynamics
    !

    ! ==================================
    ! Take timestep
    ! ==================================

    dt_vis = dt
    if (method==0) then
       ! regular LF step
       dt2 = 2*dt
       call t_startf("LF_timestep")
       call compute_and_apply_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w)
       call t_stopf("LF_timestep")
       dt_vis = dt2  ! dt to use for time-split dissipation
    else if (method==1) then
       ! RK2
       ! forward euler to u(dt/2) = u(0) + (dt/2) RHS(0)  (store in u(np1))
       call t_startf("RK2_timestep")
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0)
       ! leapfrog:  u(dt) = u(0) + dt RHS(dt/2)     (store in u(np1))
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w)
       call t_stopf("RK2_timestep")
    else if (method==2) then
       ! RK2-SSP 3 stage.  matches tracer scheme. optimal SSP CFL, but
       ! not optimal for regular CFL
       ! u1 = u0 + dt/2 RHS(u0)
       call t_startf("RK2-SSP3_timestep")
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w/3)
       ! u2 = u1 + dt/2 RHS(u1)
       call compute_and_apply_rhs(np1,np1,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w/3)
       ! u3 = u2 + dt/2 RHS(u2)
       call compute_and_apply_rhs(np1,np1,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w/3)
       ! unew = u/3 +2*u3/3  = u + 1/3 (RHS(u) + RHS(u1) + RHS(u2))
       do ie=nets,nete
          elem(ie)%state%v(:,:,:,:,np1)= elem(ie)%state%v(:,:,:,:,n0)/3 &
               + 2*elem(ie)%state%v(:,:,:,:,np1)/3
          elem(ie)%state%T(:,:,:,np1)= elem(ie)%state%T(:,:,:,n0)/3 &
               + 2*elem(ie)%state%T(:,:,:,np1)/3
          elem(ie)%state%dp3d(:,:,:,np1)= elem(ie)%state%dp3d(:,:,:,n0)/3 &
               + 2*elem(ie)%state%dp3d(:,:,:,np1)/3
       enddo
       call t_stopf("RK2-SSP3_timestep")
    else if (method==3) then
       ! classic RK3  CFL=sqrt(3)
       ! u1 = u0 + dt/3 RHS(u0)
       call t_startf("RK3_timestep")
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0)
       ! u2 = u0 + dt/2 RHS(u1)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u3 = u0 + dt RHS(u2)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w)
       call t_stopf("RK3_timestep")
    else if (method==4) then
       ! KG 4th order 4 stage:   CFL=sqrt(8)
       ! low storage version of classic RK4
       ! u1 = u0 + dt/4 RHS(u0)
       call t_startf("RK4_timestep")
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0)
       ! u2 = u0 + dt/3 RHS(u1)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u3 = u0 + dt/2 RHS(u2)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u4 = u0 + dt RHS(u3)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w)
       call t_stopf("RK4_timestep")
    else if (method==5) then
#if 0
       ! KG 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
       ! but nonlinearly only 2nd order
       ! u1 = u0 + dt/5 RHS(u0)
       call t_startf("KG3-5stage_timestep")
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0)
       ! u2 = u0 + dt/5 RHS(u1)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u3 = u0 + dt/3 RHS(u2)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u4 = u0 + dt/2 RHS(u3)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u5 = u0 + dt RHS(u4)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w)
       call t_stopf("KG3-5stage_timestep")
#else
       ! Ullrich 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
       ! u1 = u0 + dt/5 RHS(u0)  (save u1 in timelevel nm1)
       call t_startf("U3-5stage_timestep")
       call compute_and_apply_rhs(nm1,n0,n0,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w/4)
       ! u2 = u0 + dt/5 RHS(u1)
       call compute_and_apply_rhs(np1,n0,nm1,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u3 = u0 + dt/3 RHS(u2)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u4 = u0 + 2dt/3 RHS(u3)
       call compute_and_apply_rhs(np1,n0,np1,qn0,2*dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)

       ! compute (5*u1/4 - u0/4) in timelevel nm1:
       do ie=nets,nete
          elem(ie)%state%v(:,:,:,:,nm1)= (5*elem(ie)%state%v(:,:,:,:,nm1) &
               - elem(ie)%state%v(:,:,:,:,n0) ) /4
          elem(ie)%state%T(:,:,:,nm1)= (5*elem(ie)%state%T(:,:,:,nm1) &
               - elem(ie)%state%T(:,:,:,n0) )/4
          elem(ie)%state%dp3d(:,:,:,nm1)= (5*elem(ie)%state%dp3d(:,:,:,nm1) &
                  - elem(ie)%state%dp3d(:,:,:,n0) )/4
       enddo
       ! u5 = (5*u1/4 - u0/4) + 3dt/4 RHS(u4)
       call compute_and_apply_rhs(np1,nm1,np1,qn0,3*dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,3*eta_ave_w/4)
       ! final method is the same as:
       ! u5 = u0 +  dt/4 RHS(u0)) + 3dt/4 RHS(u4)
       call t_stopf("U3-5stage_timestep")
#endif

    else if ((method==11).or.(method==12)) then
       ! Fully implicit JFNK method (vertically langragian not active yet)
       if (rsplit > 0) then
       call abortmp('ERROR: full_imp integration not yet coded for vert lagrangian adv option')
       end if
!      if (hybrid%masterthread) print*, "fully implicit integration is still under development"

#ifdef TRILINOS
      call t_startf("JFNK_imp_timestep")
      lenx=(np*np*nlev*3 + np*np*1)*(nete-nets+1)  ! 3 3d vars plus 1 2d vars
      allocate(xstate(lenx))
      xstate(:) = 0d0

      call initialize(state_object, method, elem, hvcoord, compute_diagnostics, &
        qn0, eta_ave_w, hybrid, deriv, dt, tl, nets, nete)

      call initialize(pre_object, method, elem, hvcoord, compute_diagnostics, &
        qn0, eta_ave_w, hybrid, deriv, dt, tl, nets, nete)

      call initialize(jac_object, method, elem, hvcoord, compute_diagnostics, &
        qn0, eta_ave_w, hybrid, deriv, dt, tl, nets, nete)

!      pc_elem = elem
!      jac_elem = elem

        fptr => state_object
        c_ptr_to_object =  c_loc(fptr)
        pptr => state_object
        c_ptr_to_pre =  c_loc(pptr)
        jptr => state_object
        c_ptr_to_jac =  c_loc(jptr)

! create flat state vector to pass through NOX
! use previous time step as the first guess for the new one (because with LF time level update n0=np1)

       np1 = n0

       lx = 1
	   do ie=nets,nete
		   do k=1,nlev
			   do j=1,np
				   do i=1,np
					   xstate(lx) = elem(ie)%state%v(i,j,1,k,n0)
					   lx = lx+1
				   end do
			   end do
		   end do
	   end do
	   do ie=nets,nete
		   do k=1,nlev
			   do j=1,np
				   do i=1,np
					   xstate(lx) = elem(ie)%state%v(i,j,2,k,n0)
					   lx = lx+1
				   end do
			   end do
		   end do
	   end do
	   do ie=nets,nete
		   do k=1,nlev
			   do j=1,np
				   do i=1,np
					   xstate(lx) = elem(ie)%state%T(i,j,k,n0)
					   lx = lx+1
				   end do
			   end do
		   end do
	   end do
	   do ie=nets,nete
		   do j=1,np
			   do i=1,np
				   xstate(lx) = elem(ie)%state%ps_v(i,j,n0)
				   lx = lx+1
			   end do
		   end do
	   end do
       
! activate these lines to test infrastructure and still solve with explicit code
!       ! RK2
!       ! forward euler to u(dt/2) = u(0) + (dt/2) RHS(0)  (store in u(np1))
!       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&
!            deriv,nets,nete,compute_diagnostics,0d0)
!       ! leapfrog:  u(dt) = u(0) + dt RHS(dt/2)     (store in u(np1))
!       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
!            deriv,nets,nete,.false.,eta_ave_w)

! interface to use nox and loca solver libraries using JFNK, and returns xstate(n+1)
    call noxsolve(size(xstate), xstate, c_ptr_to_object, c_ptr_to_pre, c_ptr_to_jac, ierr)

    if (ierr /= 0) call abortmp('Error in noxsolve: Newton failed to converge')

      call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr
      elem = fptr%base

	  lx = 1
	  do ie=nets,nete
		  do k=1,nlev
			  do j=1,np
				  do i=1,np
					  elem(ie)%state%v(i,j,1,k,np1) = xstate(lx)
					  lx = lx+1
				  end do
			  end do
		  end do
	  end do
	  do ie=nets,nete
		  do k=1,nlev
			  do j=1,np
				  do i=1,np
					  elem(ie)%state%v(i,j,2,k,np1) = xstate(lx) 
					  lx = lx+1
				  end do
			  end do
		  end do
	  end do
	  do ie=nets,nete
		  do k=1,nlev
			  do j=1,np
				  do i=1,np
					  elem(ie)%state%T(i,j,k,np1) = xstate(lx)
					  lx = lx+1
				  end do
			  end do
		  end do
	  end do
	  do ie=nets,nete
		  do j=1,np
			  do i=1,np
				  elem(ie)%state%ps_v(i,j,np1) = xstate(lx)
				  lx = lx+1
			  end do
		  end do
	  end do
      call t_stopf("JFNK_imp_timestep")
#endif

    else
       call abortmp('ERROR: bad choice of tstep_type')
    endif

!    call prim_printstate(elem,tl,hybrid,hvcoord,nets,nete)

    ! ==============================================
    ! Time-split Horizontal diffusion: nu.del^2 or nu.del^4
    ! U(*) = U(t+1)  + dt2 * HYPER_DIFF_TERM(t+1)
    ! ==============================================
#ifdef ENERGY_DIAGNOSTICS
    if (compute_diagnostics) then
       do ie = nets,nete
          elem(ie)%accum%DIFF(:,:,:,:)=elem(ie)%state%v(:,:,:,:,np1)
          elem(ie)%accum%DIFFT(:,:,:)=elem(ie)%state%T(:,:,:,np1)
       enddo
    endif
#endif

    ! note:time step computes u(t+1)= u(t*) + RHS.
    ! for consistency, dt_vis = t-1 - t*, so this is timestep method dependent
    if (tstep_type==0) then
       ! leapfrog special case
       call advance_hypervis_lf(edge3p1,elem,hvcoord,hybrid,deriv,nm1,n0,np1,nets,nete,dt_vis)
    else if (method<=10) then ! not implicit
       ! forward-in-time, hypervis applied to dp3d
       call advance_hypervis_dp(edge3p1,elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)
    endif

#ifdef ENERGY_DIAGNOSTICS
    if (compute_diagnostics) then
       do ie = nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
         do k=1,nlev  !  Loop index added (AAM)
          elem(ie)%accum%DIFF(:,:,:,k)=( elem(ie)%state%v(:,:,:,k,np1) -&
               elem(ie)%accum%DIFF(:,:,:,k) ) / dt_vis
          elem(ie)%accum%DIFFT(:,:,k)=( elem(ie)%state%T(:,:,k,np1) -&
               elem(ie)%accum%DIFFT(:,:,k) ) / dt_vis
         enddo
       enddo
    endif
#endif

    tevolve=tevolve+dt

    call t_stopf('prim_advance_exp')
!pw call t_adj_detailf(-1)
    end subroutine prim_advance_exp


subroutine prim_advance_si(elem, nets, nete, cg, blkjac, red, &
          refstate, hvcoord, deriv, flt, hybrid, tl, dt)
       use bndry_mod, only : bndry_exchangev
       use cg_mod, only : cg_t, cg_create
       use control_mod, only : filter_freq,debug_level, precon_method
       use derivative_mod, only : derivative_t, vorticity, divergence, gradient, gradient_wk
       use edge_mod, only : edgevpack, edgevunpack, initEdgeBuffer
       use edgetype_mod, only : EdgeBuffer_t
       use filter_mod, only : filter_t, preq_filter
       use hybrid_mod, only : hybrid_t
       use prim_si_ref_mod, only : ref_state_t, set_vert_struct_mat
       use reduction_mod, only : reductionbuffer_ordered_1d_t
       use solver_mod, only : pcg_solver, blkjac_t, blkjac_init
       use prim_si_mod, only : preq_vertadv, preq_omegap, preq_pressure
       use diffusion_mod, only :  prim_diffusion
       use physical_constants, only : kappa, rrearth, rgas, cp, rwater_vapor
       use physics_mod, only : virtual_temperature, virtual_specific_heat
       implicit none

       integer, intent(in)               :: nets,nete
       type (element_t), intent(inout), target :: elem(:)
       type (blkjac_t), allocatable      :: blkjac(:)

       type (cg_t)                       :: cg

       type (ReductionBuffer_ordered_1d_t), intent(inout) :: red

       type (ref_state_t), intent(in), target :: refstate
       type (hvcoord_t), intent(in)      :: hvcoord
       type (derivative_t), intent(in)   :: deriv
       type (filter_t), intent(in)       :: flt
       type (hybrid_t), intent(in)       :: hybrid
       type (TimeLevel_t), intent(in)    :: tl
       real(kind=real_kind), intent(in)  :: dt
       real(kind=real_kind)              :: time_adv
#ifndef _CRAYFTN
       ! ==========================
       ! Local variables...
       ! ==========================

       real(kind=real_kind)                           :: ps0
       real(kind=real_kind)                           :: psref

       real(kind=real_kind), dimension(np,np)         :: ps
       real(kind=real_kind), dimension(np,np)         :: rps
       real(kind=real_kind), dimension(np,np,nlev)    :: rpmid
       real(kind=real_kind), dimension(np,np,nlev)    :: omegap
       real(kind=real_kind), dimension(np,np,nlev)    :: rpdel

       real(kind=real_kind) :: pintref(nlevp)
       real(kind=real_kind) :: pdelref(nlev)
       real(kind=real_kind) :: pmidref(nlev)
       real(kind=real_kind) :: rpdelref(nlev)
       real(kind=real_kind) :: rpmidref(nlev)

       real(kind=real_kind) :: pint(np,np,nlevp)
       real(kind=real_kind) :: pdel(np,np,nlev)
       real(kind=real_kind) :: pmid(np,np,nlev)

       real(kind=real_kind), dimension(np,np,nlevp) :: eta_dot_dp_deta
       real(kind=real_kind), dimension(np,np,nlev)  :: vgrad_ps

       real(kind=real_kind), dimension(np,np,nlev)   :: T_vadv
       real(kind=real_kind), dimension(np,np,2,nlev) :: v_vadv

       real(kind=real_kind), dimension(np,np)      :: HT
       real(kind=real_kind), dimension(np,np)      :: HrefT
       real(kind=real_kind), dimension(np,np)      :: HrefTm1

       real(kind=real_kind), dimension(np,np)      :: Gref0
       real(kind=real_kind), dimension(np,np)      :: Grefm1
       real(kind=real_kind), dimension(np,np)      :: E
       real(kind=real_kind), dimension(np,np)      :: Phi
       real(kind=real_kind), dimension(np,np)      :: dGref

       real(kind=real_kind), dimension(np,np,2)    :: vco
       real(kind=real_kind), dimension(np,np,2)    :: gradT
       real(kind=real_kind), dimension(np,np,2)    :: grad_Phi

       real(kind=real_kind), dimension(:,:), pointer  :: Emat
       real(kind=real_kind), dimension(:,:), pointer  :: Emat_inv
       real(kind=real_kind), dimension(:,:), pointer  :: Amat
       real(kind=real_kind), dimension(:,:), pointer  :: Amat_inv
       real(kind=real_kind), dimension(:), pointer    :: Lambda

       real(kind=real_kind), dimension(:), pointer    :: Tref
       real(kind=real_kind), dimension(:), pointer    :: RTref
       real(kind=real_kind), dimension(:), pointer    :: Pvec
       real(kind=real_kind), dimension(:,:), pointer  :: Href
       real(kind=real_kind), dimension(:,:), pointer  :: Tmat

       real(kind=real_kind) :: Vscript(np,np,2,nlev,nets:nete)
       real(kind=real_kind) :: Tscript(np,np,nlev,nets:nete)
       real(kind=real_kind) :: Pscript(np,np,nets:nete)
       real(kind=real_kind) :: Vtemp(np,np,2,nlev,nets:nete)

       real(kind=real_kind), dimension(np,np)      :: HrefTscript
       real(kind=real_kind), dimension(np,np)      :: suml
       real(kind=real_kind), dimension(np,np,2)    :: gVscript
       real(kind=real_kind), dimension(np,np,nlev) :: div_Vscript

       real(kind=real_kind) :: B(np,np,nlev,nets:nete)
       real(kind=real_kind) :: C(np,np,nlev,nets:nete)
       real(kind=real_kind) :: D(np,np,nlev,nets:nete)

       real(kind=real_kind) :: Gamma_ref(np,np,nlev,nets:nete)

       real(kind=real_kind) :: Gref(np,np,nlev,nets:nete)
       real(kind=real_kind) :: grad_dGref(np,np,2,nlev)
       real(kind=real_kind) :: grad_Gref(np,np,2,nlev)

       real(kind=real_kind) :: div(np,np)
       real(kind=real_kind) :: gv(np,np,2)

       real(kind=real_kind) :: dt2
       real(kind=real_kind) :: rpsref
       real(kind=real_kind) :: rdt
       real(kind=real_kind) :: hkk, hkl
       real(kind=real_kind) :: ddiv

       real(kind=real_kind) :: vgradT
       real(kind=real_kind) :: hybfac
       real(kind=real_kind) :: Crkk
       real(kind=real_kind) :: v1,v2
       real(kind=real_kind) :: term

       real(kind=real_kind) :: Vs1,Vs2
       real(kind=real_kind) :: glnps1, glnps2
       real(kind=real_kind) :: gGr1,gGr2

       real (kind=real_kind),allocatable :: solver_wts(:,:)  ! solver weights array for nonstaggered grid

       integer              :: nm1,n0,np1,nfilt
       integer              :: nstep
       integer              :: i,j,k,l,ie,kptr

!JMD       call t_barrierf('sync_prim_advance_si', hybrid%par%comm)
!pw    call t_adj_detailf(+1)
       call t_startf('prim_advance_si')

       nm1   = tl%nm1
       n0    = tl%n0
       np1   = tl%np1
       nstep = tl%nstep


       if ( dt /= initialized_for_dt ) then
          if(hybrid%par%masterproc) print *,'Initializing semi-implicit matricies for dt=',dt

#if (defined HORIZ_OPENMP)
          !$OMP MASTER
#endif
          call set_vert_struct_mat(dt, refstate, hvcoord, hybrid%masterthread)
#if (defined HORIZ_OPENMP)
          !$OMP END MASTER
#endif

          allocate(solver_wts(np*np,nete-nets+1))
          do ie=nets,nete
             kptr=1
             do j=1,np
                do i=1,np

                   ! so this code is BFB  with old code.  should change to simpler formula below
                   solver_wts(kptr,ie-nets+1) = 1d0/nint(1d0/(elem(ie)%mp(i,j)*elem(ie)%rmp(i,j)))
                   !solver_wts(kptr,ie-nets+1) = elem(ie)%mp(i,j)*elem(ie)%rmp(i,j)

                   kptr=kptr+1
                end do
             end do
          end do
          call cg_create(cg, np*np, nlev, nete-nets+1, hybrid, debug_level, solver_wts)
          deallocate(solver_wts)
          if (precon_method == "block_jacobi") then
             if (.not. allocated(blkjac)) then
                allocate(blkjac(nets:nete))
             endif
             call blkjac_init(elem, deriv,refstate%Lambda,nets,nete,blkjac)
          end if
          initialized_for_dt = dt
       endif


       nfilt = tl%nm1     ! time level at which filter is applied (time level n)
       dt2   = 2.0_real_kind*dt
       rdt   = 1.0_real_kind/dt

       ps0      = hvcoord%ps0
       psref    = refstate%psr

       Emat     => refstate%Emat
       Emat_inv => refstate%Emat_inv
       Amat     => refstate%Amat
       Amat_inv => refstate%Amat_inv
       Lambda   => refstate%Lambda

       RTref    => refstate%RTref
       Tref     => refstate%Tref
       Href     => refstate%Href
       Tmat     => refstate%Tmat
       Pvec     => refstate%Pvec

       ! ============================================================
       ! If the time is right, apply a filter to the state variables
       ! ============================================================

       if (nstep > 0 .and. filter_freq > 0 .and. MODULO(nstep,filter_freq) == 0 ) then
          call preq_filter(elem, edge3p1, flt, cg%hybrid, nfilt, nets, nete)
       end if

       ! ================================================
       ! boundary exchange grad_lnps
       ! ================================================

       do ie = nets, nete

          elem(ie)%derived%grad_lnps(:,:,:) = gradient(elem(ie)%state%lnps(:,:,n0),deriv)*rrearth

          do k=1,nlevp
             pintref(k)  = hvcoord%hyai(k)*ps0 + hvcoord%hybi(k)*psref
          end do

          do k=1,nlev
             pmidref(k)  = hvcoord%hyam(k)*ps0 + hvcoord%hybm(k)*psref
             pdelref(k)  = pintref(k+1) - pintref(k)
             rpmidref(k) = 1.0_real_kind/pmidref(k)
             rpdelref(k) = 1.0_real_kind/pdelref(k)
          end do

          rpsref   = 1.0_real_kind/psref

          ps(:,:) = EXP(elem(ie)%state%lnps(:,:,n0))
          rps(:,:) = 1.0_real_kind/ps(:,:)

          call preq_pressure(ps0,ps,hvcoord%hyai,hvcoord%hybi,hvcoord%hyam,hvcoord%hybm,pint,pmid,pdel)

          rpmid = 1.0_real_kind/pmid
          rpdel = 1.0_real_kind/pdel

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,v1,v2)
#endif
          do k=1,nlev
             do j=1,np
                do i=1,np
                   v1     = elem(ie)%state%v(i,j,1,k,n0)
                   v2     = elem(ie)%state%v(i,j,2,k,n0)

                   ! Contravariant velocities

                   vco(i,j,1) = elem(ie)%Dinv(i,j,1,1)*v1 + elem(ie)%Dinv(i,j,1,2)*v2
                   vco(i,j,2) = elem(ie)%Dinv(i,j,2,1)*v1 + elem(ie)%Dinv(i,j,2,2)*v2

                   vgrad_ps(i,j,k) = ps(i,j)*(vco(i,j,1)*elem(ie)%derived%grad_lnps(i,j,1) + &
                        vco(i,j,2)*elem(ie)%derived%grad_lnps(i,j,2))

                end do
             end do
          end do

          call preq_omegap(elem(ie)%derived%div(:,:,:,n0),vgrad_ps,pdel,rpmid, &
               hvcoord%hybm,hvcoord%hybd,elem(ie)%derived%omega_p)

          Pscript(:,:,ie)        = 0.0_real_kind
          eta_dot_dp_deta(:,:,1) = 0.0_real_kind

          do k=1,nlev
             do j=1,np
                do i=1,np
                   eta_dot_dp_deta(i,j,k+1) = eta_dot_dp_deta(i,j,k) + &
                        vgrad_ps(i,j,k)*hvcoord%hybd(k) + elem(ie)%derived%div(i,j,k,n0)*pdel(i,j,k)
                   ddiv = elem(ie)%derived%div(i,j,k,n0) - 0.5_real_kind*elem(ie)%derived%div(i,j,k,nm1)
                   Pscript(i,j,ie) = Pscript(i,j,ie) + ddiv*pdelref(k)
                end do
             end do
          end do

          do j=1,np
             do i=1,np
                Pscript(i,j,ie) = elem(ie)%state%lnps(i,j,nm1) + &
                     dt2*( rpsref*Pscript(i,j,ie) - rps(i,j)*eta_dot_dp_deta(i,j,nlev+1) )
             end do
          end do

          do k=1,nlev-1
             do j=1,np
                do i=1,np
                   eta_dot_dp_deta(i,j,k+1) = hvcoord%hybi(k+1)*eta_dot_dp_deta(i,j,nlev+1) - &
                        eta_dot_dp_deta(i,j,k+1)
                end do
             end do
          end do

          eta_dot_dp_deta(:,:,nlev+1) = 0.0_real_kind

          call preq_vertadv(elem(ie)%state%T(:,:,:,n0),elem(ie)%state%v(:,:,:,:,n0), &
               eta_dot_dp_deta,rpdel,T_vadv,v_vadv)

          suml(:,:) = 0.0_real_kind

          do k=1,nlev

             gradT(:,:,:) = gradient(elem(ie)%state%T(:,:,k,n0),deriv)*rrearth
             Crkk       = 0.5_real_kind

             do j=1,np
                do i=1,np
                   term = Crkk*(elem(ie)%derived%div(i,j,k,n0) - &
                        0.5_real_kind*elem(ie)%derived%div(i,j,k,nm1))*pdelref(k)
                   suml(i,j)  = suml(i,j) + term

                   v1     = elem(ie)%state%v(i,j,1,k,n0)
                   v2     = elem(ie)%state%v(i,j,2,k,n0)

                   ! Contravariant velocities

                   vco(i,j,1) = elem(ie)%Dinv(i,j,1,1)*v1 + elem(ie)%Dinv(i,j,1,2)*v2
                   vco(i,j,2) = elem(ie)%Dinv(i,j,2,1)*v1 + elem(ie)%Dinv(i,j,2,2)*v2

                   vgradT = vco(i,j,1)*gradT(i,j,1) + vco(i,j,2)*gradT(i,j,2)

                   Tscript(i,j,k,ie) = elem(ie)%state%T(i,j,k,nm1) &
                        + dt2*(- vgradT - T_vadv(i,j,k)           &
                        + kappa*(elem(ie)%state%T(i,j,k,n0)*elem(ie)%derived%omega_p(i,j,k) &
                        + Tref(k)*rpmidref(k)*suml(i,j)))
                   suml(i,j)  = suml(i,j) + term
                end do
             end do
          end do

          HrefT(:,:)   = 0.0_real_kind
          HrefTm1(:,:) = 0.0_real_kind
          HT(:,:)      = 0.0_real_kind

          do k=nlev,1,-1

             do j=1,np
                do i=1,np
                   hkl = rpmidref(k)*pdelref(k)
                   hkk = hkl*0.5_real_kind
                   Gref0(i,j)   = HrefT(i,j) + Rgas*hkk*elem(ie)%state%T(i,j,k,n0)
                   HrefT(i,j)   = HrefT(i,j) + Rgas*hkl*elem(ie)%state%T(i,j,k,n0)
                   Grefm1(i,j)  = HrefTm1(i,j) + Rgas*hkk*elem(ie)%state%T(i,j,k,nm1)
                   HrefTm1(i,j) = HrefTm1(i,j) + Rgas*hkl*elem(ie)%state%T(i,j,k,nm1)
                   hkl = rpmid(i,j,k)*pdel(i,j,k)
                   hkk = hkl*0.5_real_kind
                   Phi(i,j) = HT(i,j) + Rgas*hkk*elem(ie)%state%T(i,j,k,n0)
                   HT(i,j)  = HT(i,j) + Rgas*hkl*elem(ie)%state%T(i,j,k,n0)
                end do
             end do

             do j=1,np
                do i=1,np
                   v1     = elem(ie)%state%v(i,j,1,k,n0)
                   v2     = elem(ie)%state%v(i,j,2,k,n0)

                   ! covariant velocity

                   vco(i,j,1) = elem(ie)%D(i,j,1,1)*v1 + elem(ie)%D(i,j,2,1)*v2
                   vco(i,j,2) = elem(ie)%D(i,j,1,2)*v1 + elem(ie)%D(i,j,2,2)*v2

                   E(i,j) = 0.5_real_kind*( v1*v1 + v2*v2 )

                   Gref0(i,j)  =  Gref0(i,j)  + elem(ie)%state%phis(i,j) + RTref(k)*elem(ie)%state%lnps(i,j,n0)
                   Grefm1(i,j) =  Grefm1(i,j) + elem(ie)%state%phis(i,j) + RTref(k)*elem(ie)%state%lnps(i,j,nm1)

                   Phi(i,j)    =  Phi(i,j) + E(i,j) + elem(ie)%state%phis(i,j)
                   dGref(i,j)  =  -(Gref0(i,j)  - 0.5_real_kind*Grefm1(i,j))
                end do
             end do

             elem(ie)%derived%zeta(:,:,k) = vorticity(vco,deriv)*rrearth
             grad_Phi(:,:,:)     = gradient(Phi,deriv)*rrearth
             grad_dGref(:,:,:,k) = gradient_wk(dGref,deriv)*rrearth

             do j=1,np
                do i=1,np

                   elem(ie)%derived%zeta(i,j,k) = elem(ie)%rmetdet(i,j)*elem(ie)%derived%zeta(i,j,k)
                   hybfac =  hvcoord%hybm(k)*(ps(i,j)*rpmid(i,j,k))

                   glnps1 = elem(ie)%Dinv(i,j,1,1)*elem(ie)%derived%grad_lnps(i,j,1) + &
                        elem(ie)%Dinv(i,j,2,1)*elem(ie)%derived%grad_lnps(i,j,2)
                   glnps2 = elem(ie)%Dinv(i,j,1,2)*elem(ie)%derived%grad_lnps(i,j,1) + &
                        elem(ie)%Dinv(i,j,2,2)*elem(ie)%derived%grad_lnps(i,j,2)

                   v1 = elem(ie)%Dinv(i,j,1,1)*grad_Phi(i,j,1) + elem(ie)%Dinv(i,j,2,1)*grad_Phi(i,j,2)
                   v2 = elem(ie)%Dinv(i,j,1,2)*grad_Phi(i,j,1) + elem(ie)%Dinv(i,j,2,2)*grad_Phi(i,j,2)

                   Vscript(i,j,1,k,ie) = - v_vadv(i,j,1,k) &
                        + elem(ie)%state%v(i,j,2,k,n0) * (elem(ie)%fcor(i,j) + elem(ie)%derived%zeta(i,j,k)) &
                        - v1 - Rgas*hybfac*elem(ie)%state%T(i,j,k,n0)*glnps1

                   Vscript(i,j,2,k,ie) = - v_vadv(i,j,2,k) &
                        - elem(ie)%state%v(i,j,1,k,n0) * (elem(ie)%fcor(i,j) + elem(ie)%derived%zeta(i,j,k)) &
                        - v2 - Rgas*hybfac*elem(ie)%state%T(i,j,k,n0)*glnps2

                end do
             end do

          end do

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,Vs1,Vs2)
#endif
          do k=1,nlev
             do j=1,np
                do i=1,np
                   Vs1 = elem(ie)%Dinv(i,j,1,1)*grad_dGref(i,j,1,k) + elem(ie)%Dinv(i,j,2,1)*grad_dGref(i,j,2,k)
                   Vs2 = elem(ie)%Dinv(i,j,1,2)*grad_dGref(i,j,1,k) + elem(ie)%Dinv(i,j,2,2)*grad_dGref(i,j,2,k)

                   Vscript(i,j,1,k,ie) = elem(ie)%mp(i,j)*Vscript(i,j,1,k,ie) + Vs1
                   Vscript(i,j,2,k,ie) = elem(ie)%mp(i,j)*Vscript(i,j,2,k,ie) + Vs2

                   Vscript(i,j,1,k,ie) = elem(ie)%mp(i,j)*elem(ie)%state%v(i,j,1,k,nm1) + dt2*Vscript(i,j,1,k,ie)
                   Vscript(i,j,2,k,ie) = elem(ie)%mp(i,j)*elem(ie)%state%v(i,j,2,k,nm1) + dt2*Vscript(i,j,2,k,ie)
                end do
             end do

          end do

          HrefTscript(:,:) = 0.0_real_kind

          do k=nlev,1,-1

             do j=1,np
                do i=1,np
                   hkl = rpmidref(k)*pdelref(k)
                   hkk = hkl*0.5_real_kind
                   B(i,j,k,ie)      = HrefTscript(i,j) + Rgas*hkk*Tscript(i,j,k,ie)
                   B(i,j,k,ie)      = B(i,j,k,ie) +  elem(ie)%state%phis(i,j) + RTref(k)*Pscript(i,j,ie)
                   HrefTscript(i,j) = HrefTscript(i,j) + Rgas*hkl*Tscript(i,j,k,ie)
                end do
             end do

          end do

          kptr=0
          call edgeVpack(edge2, Vscript(1,1,1,1,ie),2*nlev,kptr,ie)

       end do

       call t_startf('pasi_bexchV1')
       call bndry_exchangeV(cg%hybrid,edge2)
       call t_stopf('pasi_bexchV1')

       do ie = nets, nete

          kptr=0
          call edgeVunpack(edge2, Vscript(1,1,1,1,ie), 2*nlev, kptr, ie)
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
#if (defined COLUMN_OPENMP)
! Not sure about deriv here.
!$omp parallel do private(k,i,j)
#endif
          do k=1,nlev
             do j=1,np
                do i=1,np
                   Vscript(i,j,1,k,ie) = elem(ie)%rmp(i,j)*Vscript(i,j,1,k,ie)
                   Vscript(i,j,2,k,ie) = elem(ie)%rmp(i,j)*Vscript(i,j,2,k,ie)
                end do
             end do

             do j=1,np
                do i=1,np

                   ! Contravariant Vscript

                   gVscript(i,j,1) = elem(ie)%Dinv(i,j,1,1)*Vscript(i,j,1,k,ie) + &
                        elem(ie)%Dinv(i,j,1,2)*Vscript(i,j,2,k,ie)
                   gVscript(i,j,2) = elem(ie)%Dinv(i,j,2,1)*Vscript(i,j,1,k,ie) + &
                        elem(ie)%Dinv(i,j,2,2)*Vscript(i,j,2,k,ie)

                   gVscript(i,j,1) = elem(ie)%metdet(i,j)*gVscript(i,j,1)
                   gVscript(i,j,2) = elem(ie)%metdet(i,j)*gVscript(i,j,2)

                end do
             end do

             div_Vscript(:,:,k) = divergence(gVscript,deriv)*rrearth

          end do

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,l)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   C(i,j,k,ie) = elem(ie)%metdet(i,j)*B(i,j,k,ie)
                end do
             end do

             do l=1,nlev
                do j=1,np
                   do i=1,np
                      C(i,j,k,ie) = C(i,j,k,ie) - dt*Amat(l,k)*div_Vscript(i,j,l)
                   end do
                end do
             end do

          end do

          ! ===============================================================
          !  Weight C (the RHS of the helmholtz problem) by the mass matrix
          ! ===============================================================

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j)
#endif
          do k=1,nlev
             do j=1,np
                do i=1,np
                   C(i,j,k,ie) = elem(ie)%mp(i,j)*C(i,j,k,ie)
                end do
             end do
          end do

          ! ===================================
          ! Pack C into the edge1 buffer
          ! ===================================

          kptr=0
          call edgeVpack(edge1,C(1,1,1,ie),nlev,kptr,ie)

       end do

       ! ==================================
       ! boundary exchange C
       ! ==================================

       call t_startf('pasi_bexchV2')
       call bndry_exchangeV(cg%hybrid,edge1)
       call t_stopf('pasi_bexchV2')

       do ie=nets,nete

          ! ===================================
          ! Unpack C from the edge1 buffer
          ! ===================================

          kptr=0
          call edgeVunpack(edge1, C(1,1,1,ie), nlev, kptr, ie)

          ! ===============================================
          ! Complete global assembly by normalizing by rmp
          ! ===============================================

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   D(i,j,k,ie) = elem(ie)%rmp(i,j)*C(i,j,k,ie)
                end do
             end do

          end do

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,l)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   C(i,j,k,ie) = 0.0_real_kind
                end do
             end do

             do l=1,nlev
                do j=1,np
                   do i=1,np
                      C(i,j,k,ie) = C(i,j,k,ie) + Emat_inv(l,k)*D(i,j,l,ie)
                   end do
                end do
             end do

          end do

       end do
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
       ! ==========================================
       ! solve for Gamma_ref, given C as RHS input
       ! ==========================================

       Gamma_ref = pcg_solver(elem, &
            C,          &
            cg,         &
            red,        &
            edge1,      &
            edge2,      &
            Lambda,     &
            deriv,      &
            nets,       &
            nete,       &
            blkjac)

       ! ================================================================
       ! Backsubstitute Gamma_ref into semi-implicit system of equations
       ! to find prognostic variables at time level n+1
       ! ================================================================

       kptr=0
       do ie = nets, nete

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,l)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   Gref(i,j,k,ie) = 0.0_real_kind
                end do
             end do

             do l=1,nlev
                do j=1,np
                   do i=1,np
                      Gref(i,j,k,ie) = Gref(i,j,k,ie) + Emat(l,k)*Gamma_ref(i,j,l,ie)
                   end do
                end do
             end do

             do j=1,np
                do i=1,np
                   B(i,j,k,ie) = elem(ie)%mp(i,j) * dt * (B(i,j,k,ie) - Gref(i,j,k,ie))
                end do
             end do

          end do

          call edgeVpack(edge1,B(:,:,:,ie),nlev,kptr,ie)

       end do

       call t_startf('pasi_bexchV3')
       call bndry_exchangeV(cg%hybrid,edge1)
       call t_stopf('pasi_bexchV3')

       do ie = nets, nete

          kptr=0
          call edgeVunpack(edge1, B(:,:,:,ie), nlev, kptr, ie)
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   B(i,j,k,ie) = elem(ie)%rmp(i,j)*B(i,j,k,ie)
                end do
             end do

          end do

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,l)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   D(i,j,k,ie) = 0.0_real_kind
                end do
             end do

             do l=1,nlev
                do j=1,np
                   do i=1,np
                      D(i,j,k,ie) = D(i,j,k,ie) + Emat_inv(l,k)*B(i,j,l,ie)
                   end do
                end do
             end do

          end do

#if 1
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,l)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   elem(ie)%derived%div(i,j,k,np1) = 0.0_real_kind
                end do
             end do

             do l=1,nlev
                do j=1,np
                   do i=1,np
                      elem(ie)%derived%div(i,j,k,np1) = elem(ie)%derived%div(i,j,k,np1) + Emat(l,k)*D(i,j,l,ie)/Lambda(l)
                   end do
                end do
             end do

          end do
#endif

          do k=1,nlev

             grad_Gref(:,:,:,k)=gradient_wk(Gref(:,:,k,ie),deriv)*rrearth

             do j=1,np
                do i=1,np
                   gGr1 = grad_Gref(i,j,1,k)
                   gGr2 = grad_Gref(i,j,2,k)
                   elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%Dinv(i,j,1,1)*gGr1 + elem(ie)%Dinv(i,j,2,1)*gGr2
                   elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%Dinv(i,j,1,2)*gGr1 + elem(ie)%Dinv(i,j,2,2)*gGr2
                   Vtemp(i,j,1,k,ie) = elem(ie)%state%v(i,j,1,k,np1)
                   Vtemp(i,j,2,k,ie) = elem(ie)%state%v(i,j,2,k,np1)
                end do
             end do

             do j=1,np
                do i=1,np
                   Pscript(i,j,ie) = Pscript(i,j,ie) - dt*Pvec(k)*elem(ie)%derived%div(i,j,k,np1)
                end do
             end do


             do l=1,nlev
                do j=1,np
                   do i=1,np
                      Tscript(i,j,k,ie) = Tscript(i,j,k,ie) - dt*Tmat(l,k)*elem(ie)%derived%div(i,j,l,np1)
                   end do
                end do
             end do

          end do

          do j=1,np
             do i=1,np
                Pscript(i,j,ie) = elem(ie)%mp(i,j)*Pscript(i,j,ie)
             end do
          end do
          do k=1,nlev
             do j=1,np
                do i=1,np
                   Tscript(i,j,k,ie) = elem(ie)%mp(i,j)*Tscript(i,j,k,ie)
                end do
             end do
          end do

          ! ===============================================
          ! Pack v at time level n+1 into the edge3p1 buffer
          ! ===============================================

          kptr=0
!          call edgeVpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,ie)
          call edgeVpack(edge3p1, Vtemp(:,:,:,:,ie),2*nlev,kptr,ie)

          kptr=2*nlev
          call edgeVpack(edge3p1, Tscript(:,:,:,ie),nlev,kptr,ie)

          kptr=3*nlev
          call edgeVpack(edge3p1, Pscript(:,:,ie),1,kptr,ie)

       end do

       ! ======================================
       ! boundary exchange v at time level n+1
       ! ======================================

       call t_startf('pasi_bexchV4')
       call bndry_exchangeV(cg%hybrid,edge3p1)
       call t_stopf('pasi_bexchV4')

!KGEN START(prim_advance_si_bug1)
       do ie=nets,nete

          ! ===================================
          ! Unpack v from the edge2 buffer
          ! ===================================

          kptr=0
          call edgeVunpack(edge3p1, Vtemp(:,:,:,:,ie), 2*nlev, kptr, ie)
!JMD          call edgeVunpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, ie)
!JMD          elem(ie)%state%v(:,:,:,:,np1) = Vtemp(:,:,:,:,ie)


          kptr=2*nlev
          call edgeVunpack(edge3p1, Tscript(:,:,:,ie), nlev, kptr, ie)

          kptr=3*nlev
          call edgeVunpack(edge3p1, Pscript(:,:,ie), 1, kptr, ie)

          ! ==========================================================
          ! Complete global assembly by normalizing velocity by rmp
          ! Vscript = Vscript - dt*grad(Gref)
          ! ==========================================================
          if(iam==1) then           
!BUG  There appears to be a bug in the Intel 15.0.1 compiler that generates
!BUG  incorrect code for this loop if the following print * statement is removed.
             print *,'IAM: ',iam, ' prim_advance_si: after SUM(v(np1)) ',sum(elem(ie)%state%v(:,:,:,:,np1)) 
          endif

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
!                   elem(ie)%state%v(i,j,1,k,np1) = Vscript(i,j,1,k,ie) + dt*elem(ie)%rmp(i,j)*elem(ie)%state%v(i,j,1,k,np1)
!                   elem(ie)%state%v(i,j,2,k,np1) = Vscript(i,j,2,k,ie) + dt*elem(ie)%rmp(i,j)*elem(ie)%state%v(i,j,2,k,np1)
                   elem(ie)%state%v(i,j,1,k,np1) = Vscript(i,j,1,k,ie) + dt*elem(ie)%rmp(i,j)*Vtemp(i,j,1,k,ie)
                   elem(ie)%state%v(i,j,2,k,np1) = Vscript(i,j,2,k,ie) + dt*elem(ie)%rmp(i,j)*Vtemp(i,j,2,k,ie)
                   elem(ie)%state%T(i,j,k,np1)   = elem(ie)%rmp(i,j)*Tscript(i,j,k,ie)
                end do
             end do

          end do

          do j=1,np
             do i=1,np
                elem(ie)%state%lnps(i,j,np1) = elem(ie)%rmp(i,j)*Pscript(i,j,ie)
             end do
          end do

       end do
!KGEN END(prim_advance_si_bug1)
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif

       call prim_diffusion(elem, nets,nete,np1,deriv,dt2,cg%hybrid)
       call t_stopf('prim_advance_si')
!pw    call t_adj_detailf(-1)
#endif
  end subroutine prim_advance_si


  subroutine preq_robert3(nm1,n0,np1,elem,hvcoord,nets,nete)
  use dimensions_mod, only : np, nlev, qsize
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use time_mod, only: smooth
  use control_mod, only : integration

  implicit none
  integer              , intent(in) :: nm1,n0,np1,nets,nete
  type (hvcoord_t), intent(in)      :: hvcoord
  type (element_t)     , intent(inout) :: elem(:)


  integer :: i,j,k,ie,q
  real (kind=real_kind) :: dp
  logical :: filter_ps = .false.
  if (integration == "explicit") filter_ps = .true.

!pw  call t_adj_detailf(+1)
  call t_startf('preq_robert')
  do ie=nets,nete
     if (filter_ps) then
        elem(ie)%state%ps_v(:,:,n0) = elem(ie)%state%ps_v(:,:,n0) + smooth*(elem(ie)%state%ps_v(:,:,nm1) &
             - 2.0D0*elem(ie)%state%ps_v(:,:,n0)   + elem(ie)%state%ps_v(:,:,np1))
        elem(ie)%state%lnps(:,:,n0) = LOG(elem(ie)%state%ps_v(:,:,n0))
     else
        elem(ie)%state%lnps(:,:,n0) = elem(ie)%state%lnps(:,:,n0) + smooth*(elem(ie)%state%lnps(:,:,nm1) &
             - 2.0D0*elem(ie)%state%lnps(:,:,n0)   + elem(ie)%state%lnps(:,:,np1))
        elem(ie)%state%ps_v(:,:,n0) = EXP(elem(ie)%state%lnps(:,:,n0))
     endif

     elem(ie)%state%T(:,:,:,n0) = elem(ie)%state%T(:,:,:,n0) + smooth*(elem(ie)%state%T(:,:,:,nm1) &
          - 2.0D0*elem(ie)%state%T(:,:,:,n0)   + elem(ie)%state%T(:,:,:,np1))
     elem(ie)%state%v(:,:,:,:,n0) = elem(ie)%state%v(:,:,:,:,n0) + smooth*(elem(ie)%state%v(:,:,:,:,nm1) &
          - 2.0D0*elem(ie)%state%v(:,:,:,:,n0) + elem(ie)%state%v(:,:,:,:,np1))

  end do
  call t_stopf('preq_robert')
!pw  call t_adj_detailf(-1)

  end subroutine preq_robert3




  subroutine applyCAMforcing(elem,fvm,hvcoord,np1,np1_qdp,dt_q,nets,nete)

  use dimensions_mod, only: np, nc, nlev, qsize, ntrac
  use control_mod,    only: moisture, tracer_grid_type
  use control_mod,    only: TRACER_GRIDTYPE_GLL, TRACER_GRIDTYPE_FVM
  use physical_constants, only: Cp
  use fvm_control_volume_mod, only : fvm_struct, n0_fvm

  implicit none
  type (element_t),       intent(inout) :: elem(:)
  type(fvm_struct),       intent(inout) :: fvm(:)
  real (kind=real_kind),  intent(in)    :: dt_q
  type (hvcoord_t),       intent(in)    :: hvcoord
  integer,                intent(in)    :: np1,nets,nete,np1_qdp

  ! local
  integer :: i,j,k,ie,q
  real (kind=real_kind) :: v1,dp
  real (kind=real_kind) :: beta(np,np),E0(np,np),ED(np,np),dp0m1(np,np),dpsum(np,np)
  logical :: wet

  wet = (moisture /= "dry")

  do ie=nets,nete
     ! apply forcing to Qdp
     elem(ie)%derived%FQps(:,:,1)=0
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k,i,j,v1)
#endif
     !
     ! even when running fvm tracers we need to updates forcing on ps and qv on GLL grid
     !
     ! for fvm tracer qsize is usually 1 (qv)
     !
     do q=1,qsize
        do k=1,nlev
           do j=1,np
              do i=1,np
                 v1 = dt_q*elem(ie)%derived%FQ(i,j,k,q,1)
                 !if (elem(ie)%state%Qdp(i,j,k,q,np1) + v1 < 0 .and. v1<0) then
                 if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) + v1 < 0 .and. v1<0) then
                    !if (elem(ie)%state%Qdp(i,j,k,q,np1) < 0 ) then
                    if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) < 0 ) then
                       v1=0  ! Q already negative, dont make it more so
                    else
                       !v1 = -elem(ie)%state%Qdp(i,j,k,q,np1)
                       v1 = -elem(ie)%state%Qdp(i,j,k,q,np1_qdp)
                    endif
                 endif
                 !elem(ie)%state%Qdp(i,j,k,q,np1) = elem(ie)%state%Qdp(i,j,k,q,np1)+v1
                 elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp)+v1
                 if (q==1) then
                    elem(ie)%derived%FQps(i,j,1)=elem(ie)%derived%FQps(i,j,1)+v1/dt_q
                 endif
              enddo
           enddo
        enddo
     enddo
     ! Repeat for the fvm tracers
     do q = 1, ntrac
        do k = 1, nlev
           do j = 1, nc
              do i = 1, nc
                 v1 = fvm(ie)%fc(i,j,k,q)
                 if (fvm(ie)%c(i,j,k,q,n0_fvm) + v1 < 0 .and. v1<0) then
                    if (fvm(ie)%c(i,j,k,q,n0_fvm) < 0 ) then
                       v1 = 0  ! C already negative, dont make it more so
                    else
                       v1 = -fvm(ie)%c(i,j,k,q,n0_fvm)
                    end if
                 end if
                 fvm(ie)%c(i,j,k,q,np1_qdp) = fvm(ie)%c(i,j,k,q,n0_fvm) + v1
                 !                    if (q == 1) then
                 !!XXgoldyXX: Should update the pressure forcing here??!!??
                 !                    elem(ie)%derived%FQps(i,j,1)=elem(ie)%derived%FQps(i,j,1)+v1/dt_q
                 !                  end if
              end do
           end do
        end do
     end do

     if (wet .and. qsize>0) then
        ! to conserve dry mass in the precese of Q1 forcing:
        elem(ie)%state%ps_v(:,:,np1) = elem(ie)%state%ps_v(:,:,np1) + &
             dt_q*elem(ie)%derived%FQps(:,:,1)
     endif

#if 0
     ! disabled - energy fixers will be moving into CAM physics
     ! energy fixer for FQps term
     ! dp1 = dp0 + d(FQps)
     ! dp0-dp1 = -d(FQps)
     ! E0-E1 = sum( dp0*ED) - sum( dp1*ED) = sum( dp0-dp1) * ED )
     ! compute E0-E1
     E0=0
     do k=1,nlev
        ED(:,:) = ( 0.5d0* &
             (elem(ie)%state%v(:,:,1,k,np1)**2 + elem(ie)%state%v(:,:,2,k,np1)**2)&
             + cp*elem(ie)%state%T(:,:,k,np1)  &
             + elem(ie)%state%phis(:,:) )

        dp0m1(:,:) = -dt_q*( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%derived%FQps(:,:,1)

        E0(:,:) = E0(:,:) + dp0m1(:,:)*ED(:,:)
     enddo
     ! energy fixer:
     ! Tnew = T + beta
     ! cp*dp*beta  = E0-E1   beta = (E0-E1)/(cp*sum(dp))

     dpsum(:,:) = ( hvcoord%hyai(nlev+1) - hvcoord%hyai(1) )*hvcoord%ps0 + &
          ( hvcoord%hybi(nlev+1) - hvcoord%hybi(1) )*elem(ie)%state%ps_v(:,:,np1)

     beta(:,:)=E0(:,:)/(dpsum(:,:)*cp)
     do k=1,nlev
        elem(ie)%state%T(:,:,k,np1)=elem(ie)%state%T(:,:,k,np1)+beta(:,:)
     enddo
#endif

     ! Qdp(np1) and ps_v(np1) were updated by forcing - update Q(np1)
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k,i,j,dp)
#endif
     do q=1,qsize
        do k=1,nlev
           do j=1,np
              do i=1,np
                 dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,np1)
                 elem(ie)%state%Q(i,j,k,q) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp)/dp
              enddo
           enddo
        enddo
     enddo

     elem(ie)%state%T(:,:,:,np1)   = elem(ie)%state%T(:,:,:,np1)   + dt_q*elem(ie)%derived%FT(:,:,:,1)
     elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + dt_q*elem(ie)%derived%FM(:,:,:,:,1)

  enddo
  end subroutine applyCAMforcing



  subroutine applyCAMforcing_dynamics(elem,hvcoord,np1,dt_q,nets,nete)

  use dimensions_mod, only: np, nlev, qsize
  use element_mod,    only: element_t
  use hybvcoord_mod,  only: hvcoord_t

  implicit none
  type (element_t)     ,  intent(inout) :: elem(:)
  real (kind=real_kind),  intent(in)    :: dt_q
  type (hvcoord_t),       intent(in)    :: hvcoord
  integer,                intent(in)    :: np1,nets,nete

  integer :: i,j,k,ie,q
  real (kind=real_kind) :: v1,dp
  logical :: wet

  do ie=nets,nete
     elem(ie)%state%T(:,:,:,np1)  = elem(ie)%state%T(:,:,:,np1)    + dt_q*elem(ie)%derived%FT(:,:,:,1)
     elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + dt_q*elem(ie)%derived%FM(:,:,:,:,1)
  enddo
  end subroutine applyCAMforcing_dynamics



  subroutine advance_hypervis(edge3,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
  !
  !  take one timestep of:
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use dimensions_mod, only : np, np, nlev
  use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
  use edge_mod, only : edgevpack, edgevunpack
  use edgetype_mod, only : EdgeBuffer_t
  use bndry_mod, only : bndry_exchangev
  use viscosity_mod, only : biharmonic_wk
  use physical_constants, only: Cp
  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord

  real (kind=real_kind) :: dt2
  integer :: nets,nete

  ! local
  real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
  real (kind=real_kind) :: nu_scale, dpdn,dpdn0, nu_scale_top
  integer :: k,kptr,i,j,ie,ic,nt
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ptens
  real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
  real (kind=real_kind), dimension(np,np,nlev) :: p
  real (kind=real_kind), dimension(np,np) :: dptemp1,dptemp2


! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

  real (kind=real_kind), dimension(np,np) :: lap_p
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ptens_tmp


  if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
!JMD  call t_barrierf('sync_advance_hypervis', hybrid%par%comm)
!pw   call t_adj_detailf(+1)
  call t_startf('advance_hypervis')


  dt=dt2/hypervis_subcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  regular viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 1) then
     if (nu_p>0) call abortmp( 'ERROR: hypervis_order == 1 not coded for nu_p>0')
     do ic=1,hypervis_subcycle
        do ie=nets,nete

#if (defined COLUMN_OPENMP)
! Not sure about deriv here
!$omp parallel do private(k,lap_p,lap_v,i,j)
#endif
           do k=1,nlev
              lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              ! advance in time.  (note: DSS commutes with time stepping, so we
              ! can time advance and then DSS.  this has the advantage of
              ! not letting any discontinuties accumulate in p,v via roundoff
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)*elem(ie)%spheremp(i,j)  +  dt*nu_s*lap_p(i,j)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,1)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,2)
                 enddo
              enddo
           enddo

           kptr=0
           call edgeVpack(edge3, elem(ie)%state%T(:,:,:,nt),nlev,kptr,ie)
           kptr=nlev
           call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,ie)
        enddo

        call t_startf('ah_bexchV1')
        call bndry_exchangeV(hybrid,edge3)
        call t_stopf('ah_bexchV1')

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, elem(ie)%state%T(:,:,:,nt), nlev, kptr, ie)
           kptr=nlev
           call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, ie)

           ! apply inverse mass matrix
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j)
#endif
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,nt)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
                 enddo
              enddo
           enddo
        enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
     enddo  ! subcycle
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nu_p=0:
!   scale T dissipaton by dp  (conserve IE, dissipate T^2)
! nu_p>0
!   dont scale:  T equation IE dissipation matches (to truncation error)
!                IE dissipation from continuity equation
!                (1 deg: to about 0.1 W/m^2)
!
  if (hypervis_order == 2) then
     do ic=1,hypervis_subcycle
        call biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)
        do ie=nets,nete

           ! comptue mean flux
           if (nu_p>0) then
#if 0
              elem(ie)%derived%psdiss_ave(:,:)=&
                   elem(ie)%derived%psdiss_ave(:,:)+eta_ave_w*elem(ie)%state%ps_v(:,:,nt)/hypervis_subcycle
              elem(ie)%derived%psdiss_biharmonic(:,:)=&
                   elem(ie)%derived%psdiss_biharmonic(:,:)+eta_ave_w*pstens(:,:,ie)/hypervis_subcycle
#else
              do k=1,nlev
                 dptemp1(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,nt)
                 elem(ie)%derived%dpdiss_ave(:,:,k)=elem(ie)%derived%dpdiss_ave(:,:,k)+eta_ave_w*dptemp1(:,:)/hypervis_subcycle

                 dptemp2(:,:) = (hvcoord%hybi(k+1)-hvcoord%hybi(k))*pstens(:,:,ie)
                 elem(ie)%derived%dpdiss_biharmonic(:,:,k)=&
                      elem(ie)%derived%dpdiss_biharmonic(:,:,k)+eta_ave_w*dptemp2(:,:)/hypervis_subcycle
              enddo
#endif
           endif
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,lap_p,lap_v,nu_scale_top,dpdn,dpdn0,nu_scale,utens_tmp,vtens_tmp,ptens_tmp)
#endif
           do k=1,nlev
              nu_scale=1
              ! advace in time.
              ! note: DSS commutes with time stepping, so we can time advance and then DSS.
              ! note: weak operators alreayd have mass matrix "included"

              ! add regular diffusion in top 3 layers:
              if (nu_top>0 .and. k<=3) then
                 lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              endif
              nu_scale_top = 1
              if (k==1) nu_scale_top=4
              if (k==2) nu_scale_top=2

              do j=1,np
                 do i=1,np
                    if (nu_p==0) then
                       ! normalize so as to conserve IE
                       ! scale by 1/rho (normalized to be O(1))
                       ! dp/dn = O(ps0)*O(delta_eta) = O(ps0)/O(nlev)
                       dpdn = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,nt)
                       dpdn0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
                       nu_scale = dpdn0/dpdn
                    endif

                    ! biharmonic terms need a negative sign:
                    if (nu_top>0 .and. k<=3) then
                       utens_tmp=(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                       vtens_tmp=(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
                       ptens_tmp=nu_scale*(-nu_s*ptens(i,j,k,ie) + nu_scale_top*nu_top*lap_p(i,j) )
                    else
                       utens_tmp=-nu*vtens(i,j,1,k,ie)
                       vtens_tmp=-nu*vtens(i,j,2,k,ie)
                       ptens_tmp=-nu_scale*nu_s*ptens(i,j,k,ie)
                    endif

                    ptens(i,j,k,ie) = ptens_tmp
                    vtens(i,j,1,k,ie)=utens_tmp
                    vtens(i,j,2,k,ie)=vtens_tmp
                 enddo
              enddo
           enddo

           pstens(:,:,ie)  =  -nu_p*pstens(:,:,ie)
           kptr=0
           call edgeVpack(edge3, ptens(:,:,:,ie),nlev,kptr,ie)
           kptr=nlev
           call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,ie)
           kptr=3*nlev
           call edgeVpack(edge3,pstens(:,:,ie),1,kptr,ie)
        enddo

        call t_startf('ah_bexchV2')
        call bndry_exchangeV(hybrid,edge3)
        call t_stopf('ah_bexchV2')

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, ptens(:,:,:,ie), nlev, kptr, ie)
           kptr=nlev
           call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)


           ! apply inverse mass matrix, accumulate tendencies
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
           do k=1,nlev
              vtens(:,:,1,k,ie)=dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
              vtens(:,:,2,k,ie)=dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
              ptens(:,:,k,ie)=dt*ptens(:,:,k,ie)*elem(ie)%rspheremp(:,:)
           enddo

           ! apply hypervis to u -> u+utens:
           ! E0 = dpdn * .5*u dot u + dpdn * T  + dpdn*PHIS
           ! E1 = dpdn * .5*(u+utens) dot (u+utens) + dpdn * (T-X) + dpdn*PHIS
           ! E1-E0:   dpdn (u dot utens) + dpdn .5 utens dot utens   - dpdn X
           !      X = (u dot utens) + .5 utens dot utens
           !  alt:  (u+utens) dot utens
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    ! update v first (gives better results than updating v after heating)
                    elem(ie)%state%v(i,j,:,k,nt)=elem(ie)%state%v(i,j,:,k,nt) + &
                         vtens(i,j,:,k,ie)

                    v1=elem(ie)%state%v(i,j,1,k,nt)
                    v2=elem(ie)%state%v(i,j,2,k,nt)
                    heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
                         +ptens(i,j,k,ie)-heating/cp

                 enddo
              enddo
           enddo

           if (nu_p>0) then
              kptr=3*nlev
              call edgeVunpack(edge3, pstens(:,:,ie), 1, kptr, ie)
              pstens(:,:,ie)=dt*pstens(:,:,ie)*elem(ie)%rspheremp(:,:)
              elem(ie)%state%ps_v(:,:,nt)=elem(ie)%state%ps_v(:,:,nt) + pstens(:,:,ie)
           endif

        enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
     enddo
  endif

  call t_stopf('advance_hypervis')
!pw  call t_adj_detailf(-1)

  end subroutine advance_hypervis




  subroutine advance_hypervis_dp(edge3,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
  !
  !  take one timestep of:
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use dimensions_mod, only : np, np, nlev, nc, ntrac, max_corner_elem
  use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis, swest
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
  use derivative_mod, only : subcell_Laplace_fluxes, subcell_dss_fluxes
  use edge_mod, only : edgevpack, edgevunpack, edgeDGVunpack
  use edgetype_mod, only : EdgeBuffer_t, EdgeDescriptor_t
  use bndry_mod, only : bndry_exchangev
  use viscosity_mod, only : biharmonic_wk_dp3d
  use physical_constants, only: Cp
  use derivative_mod, only : subcell_Laplace_fluxes
  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord

  real (kind=real_kind) :: dt2
  integer :: nets,nete

  ! local
  real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
  real (kind=real_kind) :: dpdn,dpdn0, nu_scale_top
  integer :: k,kptr,i,j,ie,ic,nt
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ttens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: dptens
  real (kind=real_kind), dimension(0:np+1,0:np+1,nlev)          :: corners
  real (kind=real_kind), dimension(2,2,2)                       :: cflux
  real (kind=real_kind), dimension(nc,nc,4,nlev,nets:nete)      :: dpflux
  real (kind=real_kind), dimension(np,np,nlev) :: p
  real (kind=real_kind), dimension(np,np) :: dptemp1,dptemp2
  type (EdgeDescriptor_t)                                       :: desc


! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

  real (kind=real_kind), dimension(np,np) :: lap_t,lap_dp
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ttens_tmp,dptens_tmp

  real (kind=real_kind)                     :: temp      (np,np,nlev)
  real (kind=real_kind)                     :: laplace_fluxes(nc,nc,4)



  if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
!JMD  call t_barrierf('sync_advance_hypervis', hybrid%par%comm)
!pw   call t_adj_detailf(+1) 
  call t_startf('advance_hypervis_dp')


  dt=dt2/hypervis_subcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  regular viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 1) then
     if (nu_p>0) call abortmp( 'ERROR: hypervis_order == 1 not coded for nu_p>0')
     do ic=1,hypervis_subcycle
        do ie=nets,nete

#if (defined COLUMN_OPENMP)
! Not sure about deriv here
!$omp parallel do private(k,lap_t,lap_v,i,j)
#endif
           do k=1,nlev
              lap_t=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              ! advace in time.  (note: DSS commutes with time stepping, so we
              ! can time advance and then DSS.  this has the advantage of
              ! not letting any discontinuties accumulate in p,v via roundoff
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)*elem(ie)%spheremp(i,j)  +  dt*nu_s*lap_t(i,j)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,1)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,2)
                 enddo
              enddo
           enddo

           kptr=0
           call edgeVpack(edge3, elem(ie)%state%T(:,:,:,nt),nlev,kptr,ie)
           kptr=nlev
           call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,ie)
        enddo

        call t_startf('ahdp_bexchV1')
        call bndry_exchangeV(hybrid,edge3)
        call t_stopf('ahdp_bexchV1')

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, elem(ie)%state%T(:,:,:,nt), nlev, kptr, ie)
           kptr=nlev
           call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, ie)

           ! apply inverse mass matrix
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j)
#endif
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,nt)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
                 enddo
              enddo
           enddo
        enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
     enddo  ! subcycle
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nu_p=0:
!   scale T dissipaton by dp  (conserve IE, dissipate T^2)
! nu_p>0
!   dont scale:  T equation IE dissipation matches (to truncation error)
!                IE dissipation from continuity equation
!                (1 deg: to about 0.1 W/m^2)
!
  if (hypervis_order == 2) then
     do ic=1,hypervis_subcycle
        call biharmonic_wk_dp3d(elem,dptens,dpflux,ttens,vtens,deriv,edge3,hybrid,nt,nets,nete)

        do ie=nets,nete

           ! comptue mean flux
           if (nu_p>0) then
              elem(ie)%derived%dpdiss_ave(:,:,:)=elem(ie)%derived%dpdiss_ave(:,:,:)+&
                   eta_ave_w*elem(ie)%state%dp3d(:,:,:,nt)/hypervis_subcycle
              elem(ie)%derived%dpdiss_biharmonic(:,:,:)=elem(ie)%derived%dpdiss_biharmonic(:,:,:)+&
                   eta_ave_w*dptens(:,:,:,ie)/hypervis_subcycle
           endif
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,lap_t,lap_dp,lap_v,nu_scale_top,utens_tmp,vtens_tmp,ttens_tmp,dptens_tmp,laplace_fluxes)
#endif
           do k=1,nlev
              ! advace in time.
              ! note: DSS commutes with time stepping, so we can time advance and then DSS.
              ! note: weak operators alreayd have mass matrix "included"

              ! add regular diffusion in top 3 layers:
              if (nu_top>0 .and. k<=3) then
                 lap_t=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_dp=laplace_sphere_wk(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              endif
              nu_scale_top = 1
              if (k==1) nu_scale_top=4
              if (k==2) nu_scale_top=2

              do j=1,np
                 do i=1,np
                    ! biharmonic terms need a negative sign:
                    if (nu_top>0 .and. k<=3) then
                       utens_tmp=(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                       vtens_tmp=(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
                       ttens_tmp=(-nu_s*ttens(i,j,k,ie) + nu_scale_top*nu_top*lap_t(i,j) )
                       dptens_tmp=(-nu_p*dptens(i,j,k,ie) + nu_scale_top*nu_top*lap_dp(i,j) )
                    else
                       utens_tmp=-nu*vtens(i,j,1,k,ie)
                       vtens_tmp=-nu*vtens(i,j,2,k,ie)
                       ttens_tmp=-nu_s*ttens(i,j,k,ie)
                       dptens_tmp=-nu_p*dptens(i,j,k,ie)
                    endif
                    ttens(i,j,k,ie) = ttens_tmp
                    dptens(i,j,k,ie) =dptens_tmp
                    vtens(i,j,1,k,ie)=utens_tmp
                    vtens(i,j,2,k,ie)=vtens_tmp
                 enddo
              enddo
              if (0<ntrac) then 
                elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) - &
                                              eta_ave_w*nu_p*dpflux(:,:,:,k,ie)/hypervis_subcycle
                if (nu_top>0 .and. k<=3) then
                  laplace_fluxes=subcell_Laplace_fluxes(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),np,nc)
                  elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + &
                                           eta_ave_w*nu_scale_top*nu_top*laplace_fluxes/hypervis_subcycle
                endif
              endif

              ! NOTE: we will DSS all tendicies, EXCEPT for dp3d, where we DSS the new state
              elem(ie)%state%dp3d(:,:,k,nt) = elem(ie)%state%dp3d(:,:,k,nt)*elem(ie)%spheremp(:,:)&
                   + dt*dptens(:,:,k,ie)
              
           enddo


           kptr=0
           call edgeVpack(edge3, ttens(:,:,:,ie),nlev,kptr,ie)
           kptr=nlev
           call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,ie)
           kptr=3*nlev
           call edgeVpack(edge3,elem(ie)%state%dp3d(:,:,:,nt),nlev,kptr,ie)
        enddo

        call t_startf('ahdp_bexchV2')
        call bndry_exchangeV(hybrid,edge3)
        call t_stopf('ahdp_bexchV2')

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, ttens(:,:,:,ie), nlev, kptr, ie)
           kptr=nlev
           call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
           kptr=3*nlev
           if (0<ntrac) then
             do k=1,nlev
               temp(:,:,k) = elem(ie)%state%dp3d(:,:,k,nt) / elem(ie)%spheremp  ! STATE before DSS 
             enddo
             corners = 0.0d0
             corners(1:np,1:np,:) = elem(ie)%state%dp3d(:,:,:,nt) ! fill in interior data of STATE*mass
           endif
           call edgeVunpack(edge3, elem(ie)%state%dp3d(:,:,:,nt), nlev, kptr, ie)



           if (0<ntrac) then
             kptr=3*nlev
             desc = elem(ie)%desc
             
             call edgeDGVunpack(edge3, corners, nlev, kptr, ie) 
             corners = corners/dt
             
             do k=1,nlev
               temp(:,:,k) =  elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,nt) - temp(:,:,k)
               temp(:,:,k) =  temp(:,:,k)/dt

               call distribute_flux_at_corners(cflux, corners(:,:,k), desc%getmapP)
 
               cflux(1,1,:)   = elem(ie)%rspheremp(1,  1) * cflux(1,1,:)  
               cflux(2,1,:)   = elem(ie)%rspheremp(np, 1) * cflux(2,1,:) 
               cflux(1,2,:)   = elem(ie)%rspheremp(1, np) * cflux(1,2,:) 
               cflux(2,2,:)   = elem(ie)%rspheremp(np,np) * cflux(2,2,:) 

               elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + &
                 eta_ave_w*subcell_dss_fluxes(temp(:,:,k), np, nc, elem(ie)%metdet,cflux)/hypervis_subcycle
             end do
           endif



           ! apply inverse mass matrix, accumulate tendencies
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
           do k=1,nlev
              vtens(:,:,1,k,ie)=dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
              vtens(:,:,2,k,ie)=dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
              ttens(:,:,k,ie)=dt*ttens(:,:,k,ie)*elem(ie)%rspheremp(:,:)
              
              elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt)*elem(ie)%rspheremp(:,:)
           enddo

           ! apply hypervis to u -> u+utens:
           ! E0 = dpdn * .5*u dot u + dpdn * T  + dpdn*PHIS
           ! E1 = dpdn * .5*(u+utens) dot (u+utens) + dpdn * (T-X) + dpdn*PHIS
           ! E1-E0:   dpdn (u dot utens) + dpdn .5 utens dot utens   - dpdn X
           !      X = (u dot utens) + .5 utens dot utens
           !  alt:  (u+utens) dot utens
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,v1,v2,heating)
#endif
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    ! update v first (gives better results than updating v after heating)
                    elem(ie)%state%v(i,j,:,k,nt)=elem(ie)%state%v(i,j,:,k,nt) + &
                         vtens(i,j,:,k,ie)

                    v1=elem(ie)%state%v(i,j,1,k,nt)
                    v2=elem(ie)%state%v(i,j,2,k,nt)
                    heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
                         +ttens(i,j,k,ie)-heating/cp
                    !elem(ie)%state%dp3d(i,j,k,nt)=elem(ie)%state%dp3d(i,j,k,nt) + &
                    !     dptens(i,j,k,ie)
                 enddo
              enddo
           enddo
        enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
     enddo
  endif

  call t_stopf('advance_hypervis_dp')
!pw  call t_adj_detailf(-1)

  end subroutine advance_hypervis_dp






  subroutine advance_hypervis_lf(edge3,elem,hvcoord,hybrid,deriv,nm1,n0,nt,nets,nete,dt2)
  !
  !  take one timestep of:
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use dimensions_mod, only : np, np, nlev
  use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
  use edge_mod, only : edgevpack, edgevunpack
  use edgetype_mod, only : EdgeBuffer_t
  use bndry_mod, only : bndry_exchangev
  use viscosity_mod, only : biharmonic_wk
  use physical_constants, only: Cp
  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord

  real (kind=real_kind) :: dt2
  integer :: nets,nete

  ! local
  real (kind=real_kind) :: nu_scale, dpdn,dpdn0, nu_scale_top
  integer :: k,kptr,i,j,ie,ic,n0,nt,nm1
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ptens
  real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
  real (kind=real_kind), dimension(np,np,nlev) :: p
  real (kind=real_kind), dimension(np,np) :: dXdp


! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

  real (kind=real_kind), dimension(np,np) :: lap_p
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ptens_tmp


  if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
!JMD  call t_barrierf('sync_advance_hypervis_lf', hybrid%par%comm)
!pw   call t_adj_detailf(+1)
  call t_startf('advance_hypervis_lf')

! for non-leapfrog,nt=n0=nmt
!
!  nm1 = tl%nm1   ! heating term uses U,V at average of nt and nm1 levels
!  n0 = tl%n0     ! timelevel used for ps scaling.  use n0 for leapfrog.
!  nt = tl%np1    ! apply viscosity to this timelevel  (np1)


  dt=dt2/hypervis_subcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  regular viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 1) then
     if (nu_p>0) stop 'ERROR: hypervis_order == 1 not coded for nu_p>0'
     do ic=1,hypervis_subcycle
        do ie=nets,nete

#if (defined COLUMN_OPENMP)
! Not sure about deriv here
!$omp parallel do private(k,lap_p,lap_v,i,j)
#endif
           do k=1,nlev
              lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              ! advace in time.  (note: DSS commutes with time stepping, so we
              ! can time advance and then DSS.  this has the advantage of
              ! not letting any discontinuties accumulate in p,v via roundoff
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)*elem(ie)%spheremp(i,j)  +  dt*nu_s*lap_p(i,j)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,1)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,2)
                 enddo
              enddo
           enddo

           kptr=0
           call edgeVpack(edge3, elem(ie)%state%T(:,:,:,nt),nlev,kptr,ie)
           kptr=nlev
           call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,ie)
        enddo

        call t_startf('ahlf_bexchV1')
        call bndry_exchangeV(hybrid,edge3)
        call t_stopf('ahlf_bexchV1')

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, elem(ie)%state%T(:,:,:,nt), nlev, kptr, ie)
           kptr=nlev
           call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, ie)

           ! apply inverse mass matrix
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j)
#endif
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,nt)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
                 enddo
              enddo
           enddo
        enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
     enddo  ! subcycle
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 2) then
     do ic=1,hypervis_subcycle
        call biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)
        do ie=nets,nete

           nu_scale=1
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,lap_p,lap_v,nu_scale_top,dpdn,dpdn0,nu_scale,utens_tmp,vtens_tmp,ptens_tmp)
#endif
           do k=1,nlev
              ! advace in time.
              ! note: DSS commutes with time stepping, so we can time advance and then DSS.
              ! note: weak operators alreayd have mass matrix "included"

              ! add regular diffusion in top 3 layers:
              if (nu_top>0 .and. k<=3) then
                 lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              endif
              nu_scale_top = 1
              if (k==1) nu_scale_top=4
              if (k==2) nu_scale_top=2

              do j=1,np
                 do i=1,np
                    if (psurf_vis==0) then
                       ! normalize so as to conserve IE  (not needed when using p-surface viscosity)
                       ! scale velosity by 1/rho (normalized to be O(1))
                       ! dp/dn = O(ps0)*O(delta_eta) = O(ps0)/O(nlev)
                       dpdn = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,n0)  ! nt ?
                       dpdn0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
                       nu_scale = dpdn0/dpdn
                    endif

                    ! biharmonic terms need a negative sign:
                    if (nu_top>0 .and. k<=3) then
                       utens_tmp=(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                       vtens_tmp=(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
                       ptens_tmp=nu_scale*(-nu_s*ptens(i,j,k,ie) + nu_scale_top*nu_top*lap_p(i,j) )
                    else
                       utens_tmp=-nu*vtens(i,j,1,k,ie)
                       vtens_tmp=-nu*vtens(i,j,2,k,ie)
                       ptens_tmp=-nu_scale*nu_s*ptens(i,j,k,ie)
                    endif

                    ptens(i,j,k,ie) = ptens_tmp
                    vtens(i,j,1,k,ie)=utens_tmp
                    vtens(i,j,2,k,ie)=vtens_tmp
                 enddo
              enddo
           enddo

           pstens(:,:,ie)  =  -nu_p*pstens(:,:,ie)
           kptr=0
           call edgeVpack(edge3, ptens(:,:,:,ie),nlev,kptr,ie)
           kptr=nlev
           call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,ie)
           kptr=3*nlev
           call edgeVpack(edge3,pstens(:,:,ie),1,kptr,ie)
        enddo

        call t_startf('ahlf_bexchV2')
        call bndry_exchangeV(hybrid,edge3)
        call t_stopf('ahlf_bexchV2')

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, ptens(:,:,:,ie), nlev, kptr, ie)
           kptr=nlev
           call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
           kptr=3*nlev
           call edgeVunpack(edge3, pstens(:,:,ie), 1, kptr, ie)

           if (psurf_vis == 1 ) then
              ! apply p-surface correction
              do k=1,nlev
                 p(:,:,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,nt)
              enddo
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,dXdp)
#endif
              do k=1,nlev
                 if (k.eq.1) then
                    ! no correction needed
                 else if (k.eq.nlev) then
                    ! one-sided difference
                    dXdp = (elem(ie)%state%T(:,:,k,nt) - elem(ie)%state%T(:,:,k-1,nt)) / &
                        (p(:,:,k)-p(:,:,k-1))
                    ptens(:,:,k,ie) = ptens(:,:,k,ie) - dXdp(:,:)*hvcoord%hybm(k)*pstens(:,:,ie)
                 else
                    dXdp = (elem(ie)%state%T(:,:,k+1,nt) - elem(ie)%state%T(:,:,k-1,nt)) / &
                         (p(:,:,k+1)-p(:,:,k-1))
                    ptens(:,:,k,ie) = ptens(:,:,k,ie) - dXdp(:,:)*hvcoord%hybm(k)*pstens(:,:,ie)
                 endif
              enddo
           endif


           ! apply inverse mass matrix, accumulate tendencies
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,v1,v2,heating)
#endif
           do k=1,nlev
              do j=1,np
                 do i=1,np

                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt) + &
                         dt*elem(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt) +  &
                         dt*elem(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

                    ! better E conservation if we use v after adding in vtens:
                    v1=.5*(elem(ie)%state%v(i,j,1,k,nt)+elem(ie)%state%v(i,j,1,k,nm1))
                    v2=.5*(elem(ie)%state%v(i,j,2,k,nt)+elem(ie)%state%v(i,j,2,k,nm1))
                    heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )

                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)     + &
                         dt*elem(ie)%rspheremp(i,j)*(cp*ptens(i,j,k,ie) - heating)/cp

                 enddo
              enddo
           enddo
           elem(ie)%state%ps_v(:,:,nt)=elem(ie)%state%ps_v(:,:,nt) + dt*elem(ie)%rspheremp(:,:)*pstens(:,:,ie)
        enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
     enddo
  endif

  call t_stopf('advance_hypervis_lf')
!pw  call t_adj_detailf(-1)

  end subroutine advance_hypervis_lf


  !
  ! phl notes: output is stored in first argument. Advances from 2nd argument using tendencies evaluated at 3rd rgument: 
  ! phl: for offline winds use time at 3rd argument (same as rhs currently)
  !
  subroutine compute_and_apply_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w)
  ! ===================================
  ! compute the RHS, accumulate into u(np1) and apply DSS
  !
  !           u(np1) = u(nm1) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! This subroutine is normally called to compute a leapfrog timestep
  ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
  ! accomodated.  For example, setting nm1=np1=n0 this routine will
  ! take a forward euler step, overwriting the input with the output.
  !
  !    qn0 = timelevel used to access Qdp() in order to compute virtual Temperature
  !          qn0=-1 for the dry case
  !
  ! if  dt2<0, then the DSS'd RHS is returned in timelevel np1
  !
  ! Combining the RHS and DSS pack operation in one routine
  ! allows us to fuse these two loops for more cache reuse
  !
  ! Combining the dt advance and DSS unpack operation in one routine
  ! allows us to fuse these two loops for more cache reuse
  !
  ! note: for prescribed velocity case, velocity will be computed at
  ! "real_time", which should be the time of timelevel n0.
  !
  !
  ! ===================================
  use kinds, only : real_kind
  use dimensions_mod, only : np, nc, nlev, ntrac, max_corner_elem
  use element_mod, only : element_t,PrintElem
  use derivative_mod, only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use derivative_mod, only : subcell_div_fluxes, subcell_dss_fluxes
  use edge_mod, only : edgevpack, edgevunpack, edgeDGVunpack
  use edgetype_mod, only : edgedescriptor_t
  use bndry_mod, only : bndry_exchangev
  use control_mod, only : moisture, qsplit, use_cpstar, rsplit, swest
  use hybvcoord_mod, only : hvcoord_t

  use physical_constants, only : cp, cpwater_vapor, Rgas, kappa
  use physics_mod, only : virtual_specific_heat, virtual_temperature
  use prim_si_mod, only : preq_vertadv, preq_omega_ps, preq_hydrostatic
#if ( defined CAM )
  use control_mod, only: se_met_nudge_u, se_met_nudge_p, se_met_nudge_t, se_met_tevolve
#endif

  use time_mod, only : tevolve

  implicit none
  integer, intent(in) :: np1,nm1,n0,qn0,nets,nete
  real*8, intent(in) :: dt2
  logical, intent(in)  :: compute_diagnostics

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind) :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux

  ! local
  real (kind=real_kind), pointer, dimension(:,:,:)   :: phi
  real (kind=real_kind), pointer, dimension(:,:,:)   :: dp

  real (kind=real_kind), dimension(np,np,nlev)   :: omega_p
  real (kind=real_kind), dimension(np,np,nlev)   :: T_v
  real (kind=real_kind), dimension(np,np,nlev)   :: divdp
  real (kind=real_kind), dimension(np,np,nlev+1)   :: eta_dot_dpdn  ! half level vertical velocity on p-grid
  real (kind=real_kind), dimension(np,np)      :: sdot_sum   ! temporary field
  real (kind=real_kind), dimension(np,np,2)    :: vtemp     ! generic gradient storage
  real (kind=real_kind), dimension(np,np,2,nlev):: vdp       !                            
  real (kind=real_kind), dimension(np,np,2     ):: v         !                            
  real (kind=real_kind), dimension(np,np)      :: vgrad_T    ! v.grad(T)
  real (kind=real_kind), dimension(np,np)      :: Ephi       ! kinetic energy + PHI term
  real (kind=real_kind), dimension(np,np,2,nlev) :: grad_p
  real (kind=real_kind), dimension(np,np,2,nlev) :: grad_p_m_pmet  ! gradient(p - p_met)
  real (kind=real_kind), dimension(np,np,nlev)   :: vort       ! vorticity
  real (kind=real_kind), dimension(np,np,nlev)   :: p          ! pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: rdp        ! inverse of delta pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: T_vadv     ! temperature vertical advection
  real (kind=real_kind), dimension(np,np,nlev)   :: vgrad_p    ! v.grad(p)
  real (kind=real_kind), dimension(np,np,nlev+1) :: ph               ! half level pressures on p-grid
  real (kind=real_kind), dimension(np,np,2,nlev) :: v_vadv   ! velocity vertical advection
  real (kind=real_kind), dimension(0:np+1,0:np+1,nlev)          :: corners
  real (kind=real_kind), dimension(2,2,2)                         :: cflux
  real (kind=real_kind) ::  kappa_star(np,np,nlev)
  real (kind=real_kind) ::  vtens1(np,np,nlev)
  real (kind=real_kind) ::  vtens2(np,np,nlev)
  real (kind=real_kind) ::  ttens(np,np,nlev)
  real (kind=real_kind) ::  stashdp3d (np,np,nlev)
  real (kind=real_kind) ::  tempdp3d  (np,np)
  real (kind=real_kind) ::  tempflux  (nc,nc,4)
  type (EdgeDescriptor_t)                                       :: desc


  real (kind=real_kind) ::  cp2,cp_ratio,E,de,Qt,v1,v2
  real (kind=real_kind) ::  glnps1,glnps2,gpterm
  integer :: i,j,k,kptr,ie
  real (kind=real_kind) ::  u_m_umet, v_m_vmet, t_m_tmet 

!JMD  call t_barrierf('sync_compute_and_apply_rhs', hybrid%par%comm)

!pw call t_adj_detailf(+1)
  call t_startf('compute_and_apply_rhs')
  do ie=nets,nete
     !ps => elem(ie)%state%ps_v(:,:,n0)
     phi => elem(ie)%derived%phi(:,:,:)
     dp  => elem(ie)%state%dp3d(:,:,:,n0)

! dont thread this because of k-1 dependence:
     p(:,:,1)=hvcoord%hyai(1)*hvcoord%ps0 + dp(:,:,1)/2
     do k=2,nlev
        p(:,:,k)=p(:,:,k-1) + dp(:,:,k-1)/2 + dp(:,:,k)/2
     enddo

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,v1,v2,vtemp)
#endif
     do k=1,nlev
        grad_p(:,:,:,k) = gradient_sphere(p(:,:,k),deriv,elem(ie)%Dinv)
        rdp(:,:,k) = 1.0D0/dp(:,:,k)

        ! ============================
        ! compute vgrad_lnps
        ! ============================
        do j=1,np
           do i=1,np
              v1 = elem(ie)%state%v(i,j,1,k,n0)
              v2 = elem(ie)%state%v(i,j,2,k,n0)
              vgrad_p(i,j,k) = (v1*grad_p(i,j,1,k) + v2*grad_p(i,j,2,k))
              vdp(i,j,1,k) = v1*dp(i,j,k)
              vdp(i,j,2,k) = v2*dp(i,j,k)
           end do
        end do

#if ( defined CAM )
        ! ============================
        ! compute grad(P-P_met)
        ! ============================
        if (se_met_nudge_p.gt.0.D0) then
           grad_p_m_pmet(:,:,:,k) = &
                grad_p(:,:,:,k) - &
                hvcoord%hybm(k)* &
                gradient_sphere( elem(ie)%derived%ps_met(:,:)+tevolve*elem(ie)%derived%dpsdt_met(:,:), &
                                 deriv,elem(ie)%Dinv)
        endif
#endif

        ! ================================
        ! Accumulate mean Vel_rho flux in vn0
        ! ================================
        elem(ie)%derived%vn0(:,:,:,k)=elem(ie)%derived%vn0(:,:,:,k)+eta_ave_w*vdp(:,:,:,k)


        ! =========================================
        !
        ! Compute relative vorticity and divergence
        !
        ! =========================================
        divdp(:,:,k)=divergence_sphere(vdp(:,:,:,k),deriv,elem(ie))
        vort(:,:,k)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,n0),deriv,elem(ie))

     enddo

     ! compute T_v for timelevel n0
     !if ( moisture /= "dry") then
     if (qn0 == -1 ) then
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j)
#endif
        do k=1,nlev
           do j=1,np
              do i=1,np
                 T_v(i,j,k) = elem(ie)%state%T(i,j,k,n0)
                 kappa_star(i,j,k) = kappa
              end do
           end do
        end do
     else
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,Qt)
#endif

        do k=1,nlev
           do j=1,np
              do i=1,np
                 ! Qt = elem(ie)%state%Q(i,j,k,1)
                 Qt = elem(ie)%state%Qdp(i,j,k,1,qn0)/dp(i,j,k)
!!XXgoldyXX
!Qt=0._real_kind
!!XXgoldyXX
                 T_v(i,j,k) = Virtual_Temperature(elem(ie)%state%T(i,j,k,n0),Qt)
                 if (use_cpstar==1) then
                    kappa_star(i,j,k) =  Rgas/Virtual_Specific_Heat(Qt)
                 else
                    kappa_star(i,j,k) = kappa
                 endif
              end do
           end do
        end do
     end if


     ! ====================================================
     ! Compute Hydrostatic equation, modeld after CCM-3
     ! ====================================================
     !call geopotential_t(p,dp,T_v,Rgas,phi)
     call preq_hydrostatic(phi,elem(ie)%state%phis,T_v,p,dp)

     ! ====================================================
     ! Compute omega_p according to CCM-3
     ! ====================================================
     call preq_omega_ps(omega_p,hvcoord,p,vgrad_p,divdp)


     ! ==================================================
     ! zero partial sum for accumulating sum
     !    (div(v_k) + v_k.grad(lnps))*dsigma_k = div( v dp )
     ! used by eta_dot_dpdn and lnps tendency
     ! ==================================================
     sdot_sum=0


     ! ==================================================
     ! Compute eta_dot_dpdn
     ! save sdot_sum as this is the -RHS of ps_v equation
     ! ==================================================
     if (rsplit>0) then
        ! VERTICALLY LAGRANGIAN:   no vertical motion
        eta_dot_dpdn=0
        T_vadv=0
        v_vadv=0
     else
        do k=1,nlev
           ! ==================================================
           ! add this term to PS equation so we exactly conserve dry mass
           ! ==================================================
           sdot_sum(:,:) = sdot_sum(:,:) + divdp(:,:,k)
           eta_dot_dpdn(:,:,k+1) = sdot_sum(:,:)
        end do


        ! ===========================================================
        ! at this point, eta_dot_dpdn contains integral_etatop^eta[ divdp ]
        ! compute at interfaces:
        !    eta_dot_dpdn = -dp/dt - integral_etatop^eta[ divdp ]
        ! for reference: at mid layers we have:
        !    omega = v grad p  - integral_etatop^eta[ divdp ]
        ! ===========================================================
        do k=1,nlev-1
           eta_dot_dpdn(:,:,k+1) = hvcoord%hybi(k+1)*sdot_sum(:,:) - eta_dot_dpdn(:,:,k+1)
        end do

        eta_dot_dpdn(:,:,1     ) = 0.0D0
        eta_dot_dpdn(:,:,nlev+1) = 0.0D0

        ! ===========================================================
        ! Compute vertical advection of T and v from eq. CCM2 (3.b.1)
        ! ==============================================
        call preq_vertadv(elem(ie)%state%T(:,:,:,n0),elem(ie)%state%v(:,:,:,:,n0), &
             eta_dot_dpdn,rdp,T_vadv,v_vadv)
     endif


     ! ================================
     ! accumulate mean vertical flux:
     ! ================================
#if (defined COLUMN_OPENMP)
     !$omp parallel do private(k)
#endif
     do k=1,nlev  !  Loop index added (AAM)
        elem(ie)%derived%eta_dot_dpdn(:,:,k) = &
             elem(ie)%derived%eta_dot_dpdn(:,:,k) + eta_ave_w*eta_dot_dpdn(:,:,k)
        elem(ie)%derived%omega_p(:,:,k) = &
             elem(ie)%derived%omega_p(:,:,k) + eta_ave_w*omega_p(:,:,k)
     enddo
     elem(ie)%derived%eta_dot_dpdn(:,:,nlev+1) = &
          elem(ie)%derived%eta_dot_dpdn(:,:,nlev+1) + eta_ave_w*eta_dot_dpdn(:,:,nlev+1)





     ! ==============================================
     ! Compute phi + kinetic energy term: 10*nv*nv Flops
     ! ==============================================
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,v1,v2,E,Ephi,vtemp,vgrad_T,gpterm,glnps1,glnps2,u_m_umet,v_m_vmet,t_m_tmet)
#endif
     vertloop: do k=1,nlev
        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              E = 0.5D0*( v1*v1 + v2*v2 )
              Ephi(i,j)=E+phi(i,j,k)+elem(ie)%derived%pecnd(i,j,k)
           end do
        end do
        ! ================================================
        ! compute gradp term (ps/p)*(dp/dps)*T
        ! ================================================
        vtemp(:,:,:)   = gradient_sphere(elem(ie)%state%T(:,:,k,n0),deriv,elem(ie)%Dinv)
        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              vgrad_T(i,j) =  v1*vtemp(i,j,1) + v2*vtemp(i,j,2)
           end do
        end do


        ! vtemp = grad ( E + PHI )
        vtemp = gradient_sphere(Ephi(:,:),deriv,elem(ie)%Dinv)

        do j=1,np
           do i=1,np
!              gpterm = hvcoord%hybm(k)*T_v(i,j,k)/p(i,j,k)
!              glnps1 = Rgas*gpterm*grad_ps(i,j,1)
!              glnps2 = Rgas*gpterm*grad_ps(i,j,2)
              gpterm = T_v(i,j,k)/p(i,j,k)
              glnps1 = Rgas*gpterm*grad_p(i,j,1,k)
              glnps2 = Rgas*gpterm*grad_p(i,j,2,k)

              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)

              vtens1(i,j,k) =   - v_vadv(i,j,1,k)                           &
                   + v2*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - vtemp(i,j,1) - glnps1
              !
              ! phl: add forcing term to zonal wind u
              !
              vtens2(i,j,k) =   - v_vadv(i,j,2,k)                            &
                   - v1*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - vtemp(i,j,2) - glnps2
              !
              ! phl: add forcing term to meridional wind v
              !
              ttens(i,j,k)  = - T_vadv(i,j,k) - vgrad_T(i,j) + kappa_star(i,j,k)*T_v(i,j,k)*omega_p(i,j,k)
              !
              ! phl: add forcing term to T
              !

           end do
        end do

     end do vertloop

#ifdef ENERGY_DIAGNOSTICS
     ! =========================================================
     !
     ! diagnostics
     ! recomputes some gradients that were not saved above
     ! uses:  sdot_sum(), eta_dot_dpdn(), grad_ps()
     ! grad_phi(), dp(), p(), T_vadv(), v_vadv(), divdp()
     ! =========================================================

     ! =========================================================
     ! (AAM) - This section has accumulations over vertical levels.
     !   Be careful if implementing OpenMP
     ! =========================================================

     if (compute_diagnostics) then
        elem(ie)%accum%KEhorz1=0
        elem(ie)%accum%KEhorz2=0
        elem(ie)%accum%IEhorz1=0
        elem(ie)%accum%IEhorz2=0
        elem(ie)%accum%IEhorz1_wet=0
        elem(ie)%accum%IEhorz2_wet=0
        elem(ie)%accum%KEvert1=0
        elem(ie)%accum%KEvert2=0
        elem(ie)%accum%IEvert1=0
        elem(ie)%accum%IEvert2=0
        elem(ie)%accum%IEvert1_wet=0
        elem(ie)%accum%IEvert2_wet=0
        elem(ie)%accum%T1=0
        elem(ie)%accum%T2=0
        elem(ie)%accum%T2_s=0
        elem(ie)%accum%S1=0
        elem(ie)%accum%S1_wet=0
        elem(ie)%accum%S2=0

        do j=1,np
           do i=1,np
              elem(ie)%accum%S2(i,j) = elem(ie)%accum%S2(i,j) - &
                   sdot_sum(i,j)*elem(ie)%state%phis(i,j)
           enddo
        enddo

        do k=1,nlev
           ! vtemp = grad_E(:,:,k)
           do j=1,np
              do i=1,np
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 Ephi(i,j)=0.5D0*( v1*v1 + v2*v2 )
              enddo
           enddo
           vtemp = gradient_sphere(Ephi,deriv,elem(ie)%Dinv)
           do j=1,np
              do i=1,np
                 ! dp/dn u dot grad(E)
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 elem(ie)%accum%KEhorz2(i,j) = elem(ie)%accum%KEhorz2(i,j) + &
                      (v1*vtemp(i,j,1)  + v2*vtemp(i,j,2))*dp(i,j,k)
                 ! E div( u dp/dn )
                 elem(ie)%accum%KEhorz1(i,j) = elem(ie)%accum%KEhorz1(i,j) + Ephi(i,j)*divdp(i,j,k)

                 ! Cp T div( u dp/dn)   ! dry horizontal advection component
                 elem(ie)%accum%IEhorz1(i,j) = elem(ie)%accum%IEhorz1(i,j) + Cp*elem(ie)%state%T(i,j,k,n0)*divdp(i,j,k)


              enddo
           enddo


           ! vtemp = grad_phi(:,:,k)
           vtemp = gradient_sphere(phi(:,:,k),deriv,elem(ie)%Dinv)
           do j=1,np
              do i=1,np
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 E = 0.5D0*( v1*v1 + v2*v2 )
                 ! NOTE:  Cp_star = Cp + (Cpv-Cp)*q
                 ! advection terms can thus be broken into two components: dry and wet
                 ! dry components cancel exactly
                 ! wet components should cancel exactly
                 !
                 ! some diagnostics
                 ! e = eta_dot_dpdn()
                 de =  eta_dot_dpdn(i,j,k+1)-eta_dot_dpdn(i,j,k)
                 ! Cp T de/dn, integral dn:
                 elem(ie)%accum%IEvert1(i,j)=elem(ie)%accum%IEvert1(i,j) + Cp*elem(ie)%state%T(i,j,k,n0)*de
                 ! E de/dn
                 elem(ie)%accum%KEvert1(i,j)=elem(ie)%accum%KEvert1(i,j) + E*de
                 ! Cp T_vadv dp/dn
                 elem(ie)%accum%IEvert2(i,j)=elem(ie)%accum%IEvert2(i,j) + Cp*T_vadv(i,j,k)*dp(i,j,k)
                 ! dp/dn V dot V_vadv
                 elem(ie)%accum%KEvert2(i,j)=elem(ie)%accum%KEvert2(i,j) + (v1*v_vadv(i,j,1,k) + v2*v_vadv(i,j,2,k)) *dp(i,j,k)

                 ! IEvert1_wet():  (Cpv-Cp) T Qdp_vadv  (Q equation)
                 ! IEvert2_wet():  (Cpv-Cp) Qdp T_vadv   T equation
                 if (use_cpstar==1) then
                 elem(ie)%accum%IEvert2_wet(i,j)=elem(ie)%accum%IEvert2_wet(i,j) +&
                      (Cpwater_vapor-Cp)*elem(ie)%state%Q(i,j,k,1)*T_vadv(i,j,k)*dp(i,j,k)
                 endif

                 gpterm = T_v(i,j,k)/p(i,j,k)
                 elem(ie)%accum%T1(i,j) = elem(ie)%accum%T1(i,j) - &
                      Rgas*gpterm*(grad_p(i,j,1,k)*v1 + grad_p(i,j,2,k)*v2)*dp(i,j,k)

                 elem(ie)%accum%T2(i,j) = elem(ie)%accum%T2(i,j) - &
                      (vtemp(i,j,1)*v1 + vtemp(i,j,2)*v2)*dp(i,j,k)

                 ! S1 = < Cp_star dp/dn , RT omega_p/cp_star >
                 elem(ie)%accum%S1(i,j) = elem(ie)%accum%S1(i,j) + &
                      Rgas*T_v(i,j,k)*omega_p(i,j,k)*dp(i,j,k)

                 ! cp_star = cp + cp2
                 if (use_cpstar==1) then
                 cp2 = (Cpwater_vapor-Cp)*elem(ie)%state%Q(i,j,k,1)
                 cp_ratio = cp2/(cp+cp2)
                 elem(ie)%accum%S1_wet(i,j) = elem(ie)%accum%S1_wet(i,j) + &
                      cp_ratio*(Rgas*T_v(i,j,k)*omega_p(i,j,k)*dp(i,j,k))
                 endif

                 elem(ie)%accum%CONV(i,j,:,k)=-Rgas*gpterm*grad_p(i,j,:,k)-vtemp(i,j,:)
              enddo
           enddo

           vtemp(:,:,:) = gradient_sphere(elem(ie)%state%phis(:,:),deriv,elem(ie)%Dinv)
           do j=1,np
              do i=1,np
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 elem(ie)%accum%T2_s(i,j) = elem(ie)%accum%T2_s(i,j) - &
                      (vtemp(i,j,1)*v1 + vtemp(i,j,2)*v2)*dp(i,j,k)
              enddo
           enddo

           vtemp(:,:,:)   = gradient_sphere(elem(ie)%state%T(:,:,k,n0),deriv,elem(ie)%Dinv)
           do j=1,np
              do i=1,np
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)

                 ! Cp dp/dn u dot gradT
                 elem(ie)%accum%IEhorz2(i,j) = elem(ie)%accum%IEhorz2(i,j) + &
                      Cp*(v1*vtemp(i,j,1) + v2*vtemp(i,j,2))*dp(i,j,k)

                 if (use_cpstar==1) then
                 elem(ie)%accum%IEhorz2_wet(i,j) = elem(ie)%accum%IEhorz2_wet(i,j) + &
                      (Cpwater_vapor-Cp)*elem(ie)%state%Q(i,j,k,1)*&
                      (v1*vtemp(i,j,1) + v2*vtemp(i,j,2))*dp(i,j,k)
                 endif

              enddo
           enddo

        enddo
     endif
#endif
     ! =========================================================
     ! local element timestep, store in np1.
     ! note that we allow np1=n0 or nm1
     ! apply mass matrix
     ! =========================================================
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,tempflux)
#endif
     do k=1,nlev
        elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%spheremp(:,:)*( elem(ie)%state%v(:,:,1,k,nm1) + dt2*vtens1(:,:,k) )
        elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%spheremp(:,:)*( elem(ie)%state%v(:,:,2,k,nm1) + dt2*vtens2(:,:,k) )
        elem(ie)%state%T(:,:,k,np1) = elem(ie)%spheremp(:,:)*(elem(ie)%state%T(:,:,k,nm1) + dt2*ttens(:,:,k))
        elem(ie)%state%dp3d(:,:,k,np1) = &
             elem(ie)%spheremp(:,:) * (elem(ie)%state%dp3d(:,:,k,nm1) - &
             dt2 * (divdp(:,:,k) + eta_dot_dpdn(:,:,k+1)-eta_dot_dpdn(:,:,k)))
        
        if (0<rsplit.and.0<ntrac.and.eta_ave_w.ne.0.) then
           v(:,:,1) =  elem(ie)%Dinv(:,:,1,1)*vdp(:,:,1,k) + elem(ie)%Dinv(:,:,1,2)*vdp(:,:,2,k)
           v(:,:,2) =  elem(ie)%Dinv(:,:,2,1)*vdp(:,:,1,k) + elem(ie)%Dinv(:,:,2,2)*vdp(:,:,2,k)
           tempflux =  eta_ave_w*subcell_div_fluxes(v, np, nc, elem(ie)%metdet)
           elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) - tempflux
        end if
     enddo


     ! =========================================================
     !
     ! Pack ps(np1), T, and v tendencies into comm buffer
     !
     ! =========================================================
     kptr=0
     call edgeVpack(edge3p1, elem(ie)%state%T(:,:,:,np1),nlev,kptr,ie)

     kptr=kptr+nlev
     call edgeVpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,ie)

     kptr=kptr+2*nlev
     call edgeVpack(edge3p1, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr, ie)
  end do

  ! =============================================================
    ! Insert communications here: for shared memory, just a single
  ! sync is required
  ! =============================================================

  call t_startf('caar_bexchV')
  call bndry_exchangeV(hybrid,edge3p1)
  call t_stopf('caar_bexchV')

  do ie=nets,nete
     ! ===========================================================
     ! Unpack the edges for vgrad_T and v tendencies...
     ! ===========================================================
     kptr=0
     call edgeVunpack(edge3p1, elem(ie)%state%T(:,:,:,np1), nlev, kptr, ie)

     kptr=kptr+nlev
     call edgeVunpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, ie)

     if (0<ntrac.and.eta_ave_w.ne.0.) then
        do k=1,nlev
           stashdp3d(:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)/elem(ie)%spheremp(:,:)
        end do
        corners = 0.0d0
        corners(1:np,1:np,:) = elem(ie)%state%dp3d(:,:,:,np1)
     endif
     
     kptr=kptr+2*nlev
     call edgeVunpack(edge3p1, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr,ie)
     
     if  (0<ntrac.and.eta_ave_w.ne.0.) then
        desc = elem(ie)%desc
        call edgeDGVunpack(edge3p1, corners, nlev, kptr, ie)
        corners = corners/dt2
        
        do k=1,nlev
           tempdp3d = elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
           tempdp3d = tempdp3d - stashdp3d(:,:,k)
           tempdp3d = tempdp3d/dt2
           
           call distribute_flux_at_corners(cflux, corners(:,:,k), desc%getmapP)
           
           cflux(1,1,:)   = elem(ie)%rspheremp(1,  1) * cflux(1,1,:)  
           cflux(2,1,:)   = elem(ie)%rspheremp(np, 1) * cflux(2,1,:) 
           cflux(1,2,:)   = elem(ie)%rspheremp(1, np) * cflux(1,2,:) 
           cflux(2,2,:)   = elem(ie)%rspheremp(np,np) * cflux(2,2,:) 
           
           tempflux =  eta_ave_w*subcell_dss_fluxes(tempdp3d, np, nc, elem(ie)%metdet, cflux)
           elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + tempflux
        end do
     end if


     ! ====================================================
     ! Scale tendencies by inverse mass matrix
     ! ====================================================

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        elem(ie)%state%T(:,:,k,np1)   = elem(ie)%rspheremp(:,:)*elem(ie)%state%T(:,:,k,np1)
        elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,1,k,np1)
        elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,2,k,np1)
     end do

     ! vertically lagrangian: complete dp3d timestep:
     do k=1,nlev
        elem(ie)%state%dp3d(:,:,k,np1)= elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
     enddo
  end do

#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
  call t_stopf('compute_and_apply_rhs')
!pw  call t_adj_detailf(-1)

  end subroutine compute_and_apply_rhs


  subroutine distribute_flux_at_corners(cflux, corners, getmapP)
    use kinds,          only : int_kind, real_kind
    use dimensions_mod, only : np, max_corner_elem
    use control_mod,    only : swest
    implicit none

    real   (kind=real_kind), intent(out)  :: cflux(2,2,2)
    real   (kind=real_kind), intent(in)   :: corners(0:np+1,0:np+1)
    integer(kind=int_kind),  intent(in)   :: getmapP(:)

    cflux = 0.0d0
    if (getmapP(swest+0*max_corner_elem) /= -1) then
      cflux(1,1,1) =                (corners(0,1) - corners(1,1))     
      cflux(1,1,1) = cflux(1,1,1) + (corners(0,0) - corners(1,1)) / 2.0d0
      cflux(1,1,1) = cflux(1,1,1) + (corners(0,1) - corners(1,0)) / 2.0d0
 
      cflux(1,1,2) =                (corners(1,0) - corners(1,1))     
      cflux(1,1,2) = cflux(1,1,2) + (corners(0,0) - corners(1,1)) / 2.0d0
      cflux(1,1,2) = cflux(1,1,2) + (corners(1,0) - corners(0,1)) / 2.0d0
    else
      cflux(1,1,1) =                (corners(0,1) - corners(1,1))     
      cflux(1,1,2) =                (corners(1,0) - corners(1,1))     
    endif
 
    if (getmapP(swest+1*max_corner_elem) /= -1) then
      cflux(2,1,1) =                (corners(np+1,1) - corners(np,1))     
      cflux(2,1,1) = cflux(2,1,1) + (corners(np+1,0) - corners(np,1)) / 2.0d0
      cflux(2,1,1) = cflux(2,1,1) + (corners(np+1,1) - corners(np,0)) / 2.0d0
 
      cflux(2,1,2) =                (corners(np  ,0) - corners(np,  1))     
      cflux(2,1,2) = cflux(2,1,2) + (corners(np+1,0) - corners(np,  1)) / 2.0d0
      cflux(2,1,2) = cflux(2,1,2) + (corners(np  ,0) - corners(np+1,1)) / 2.0d0
    else
      cflux(2,1,1) =                (corners(np+1,1) - corners(np,1))     
      cflux(2,1,2) =                (corners(np  ,0) - corners(np,1))     
    endif
 
    if (getmapP(swest+2*max_corner_elem) /= -1) then
      cflux(1,2,1) =                (corners(0,np  ) - corners(1,np  ))     
      cflux(1,2,1) = cflux(1,2,1) + (corners(0,np+1) - corners(1,np  )) / 2.0d0
      cflux(1,2,1) = cflux(1,2,1) + (corners(0,np  ) - corners(1,np+1)) / 2.0d0
 
      cflux(1,2,2) =                (corners(1,np+1) - corners(1,np  ))     
      cflux(1,2,2) = cflux(1,2,2) + (corners(0,np+1) - corners(1,np  )) / 2.0d0
      cflux(1,2,2) = cflux(1,2,2) + (corners(1,np+1) - corners(0,np  )) / 2.0d0
    else
      cflux(1,2,1) =                (corners(0,np  ) - corners(1,np  ))     
      cflux(1,2,2) =                (corners(1,np+1) - corners(1,np  ))     
    endif
 
    if (getmapP(swest+3*max_corner_elem) /= -1) then
      cflux(2,2,1) =                (corners(np+1,np  ) - corners(np,np  ))     
      cflux(2,2,1) = cflux(2,2,1) + (corners(np+1,np+1) - corners(np,np  )) / 2.0d0
      cflux(2,2,1) = cflux(2,2,1) + (corners(np+1,np  ) - corners(np,np+1)) / 2.0d0
 
      cflux(2,2,2) =                (corners(np  ,np+1) - corners(np,np  ))     
      cflux(2,2,2) = cflux(2,2,2) + (corners(np+1,np+1) - corners(np,np  )) / 2.0d0
      cflux(2,2,2) = cflux(2,2,2) + (corners(np  ,np+1) - corners(np+1,np)) / 2.0d0
    else
      cflux(2,2,1) =                (corners(np+1,np  ) - corners(np,np  ))     
      cflux(2,2,2) =                (corners(np  ,np+1) - corners(np,np  ))     
    endif
  end subroutine



  subroutine smooth_phis(phis,elem,hybrid,deriv,nets,nete,minf,numcycle)
  use dimensions_mod, only : np, np, nlev
  use control_mod, only : smooth_phis_nudt, hypervis_scaling
  use hybrid_mod, only : hybrid_t
  use edge_mod, only : edgevpack, edgevunpack, edgevunpackmax, edgevunpackmin
  use edgetype_mod, only : EdgeBuffer_t
  use bndry_mod, only : bndry_exchangev
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t , laplace_sphere_wk
  use time_mod, only : TimeLevel_t
  implicit none

  integer :: nets,nete
  real (kind=real_kind), dimension(np,np,nets:nete), intent(inout)   :: phis
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind), intent(in)   :: minf
  integer,               intent(in) :: numcycle

  ! local
  real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
  real (kind=real_kind), dimension(nets:nete) :: pmin,pmax
  real (kind=real_kind) :: mx,mn
  integer :: nt,ie,ic,i,j,order,order_max, iuse
  logical :: use_var_coef


  ! compute local element neighbor min/max
  do ie=nets,nete
     pstens(:,:,ie)=minval(phis(:,:,ie))
     call edgeVpack(edge3p1,pstens(:,:,ie),1,0,ie)
  enddo

  call t_startf('smooth_phis_bexchV1')
  call bndry_exchangeV(hybrid,edge3p1)
  call t_stopf('smooth_phis_bexchV1')

  do ie=nets,nete
     call edgeVunpackMin(edge3p1, pstens(:,:,ie), 1, 0, ie)
     pmin(ie)=minval(pstens(:,:,ie))
  enddo
  do ie=nets,nete
     pstens(:,:,ie)=maxval(phis(:,:,ie))
     call edgeVpack(edge3p1,pstens(:,:,ie),1,0,ie)
  enddo

  call t_startf('smooth_phis_bexchV2')
  call bndry_exchangeV(hybrid,edge3p1)
  call t_stopf('smooth_phis_bexchV2')

  do ie=nets,nete
     call edgeVunpackMax(edge3p1, pstens(:,:,ie), 1, 0, ie)
     pmax(ie)=maxval(pstens(:,:,ie))
  enddo

  ! order = 1   grad^2 laplacian
  ! order = 2   grad^4 (need to add a negative sign)
  ! order = 3   grad^6
  ! order = 4   grad^8 (need to add a negative sign)
  order_max = 1


  use_var_coef=.true.
  if (hypervis_scaling/=0) then
     ! for tensorHV option, we turn off the tensor except for *last* laplace operator
     use_var_coef=.false.
     if (hypervis_scaling>=3) then
        ! with a 3.2 or 4 scaling, assume hyperviscosity
        order_max = 2
     endif
  endif


  do ic=1,numcycle
     pstens=phis

     do order=1,order_max-1

        do ie=nets,nete
           pstens(:,:,ie)=laplace_sphere_wk(pstens(:,:,ie),deriv,elem(ie),var_coef=use_var_coef)
           call edgeVpack(edge3p1,pstens(:,:,ie),1,0,ie)
        enddo

        call t_startf('smooth_phis_bexchV3')
        call bndry_exchangeV(hybrid,edge3p1)
        call t_stopf('smooth_phis_bexchV3')

        do ie=nets,nete
           call edgeVunpack(edge3p1, pstens(:,:,ie), 1, 0, ie)
           pstens(:,:,ie)=pstens(:,:,ie)*elem(ie)%rspheremp(:,:)
        enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
     enddo
     do ie=nets,nete
        pstens(:,:,ie)=laplace_sphere_wk(pstens(:,:,ie),deriv,elem(ie),var_coef=.true.)
     enddo
     if (mod(order_max,2)==0) pstens=-pstens

     do ie=nets,nete
        !  ps(t+1) = ps(t) + Minv * DSS * M * RHS
        !  ps(t+1) = Minv * DSS * M [ ps(t) +  RHS ]
        ! but output of biharminc_wk is of the form M*RHS.  rewrite as:
        !  ps(t+1) = Minv * DSS * M [ ps(t) +  M*RHS/M ]
        ! so we can apply limiter to ps(t) +  (M*RHS)/M
#if 1
        mn=pmin(ie)
        mx=pmax(ie)
        iuse = numcycle+1  ! always apply min/max limiter
#endif
        phis(:,:,ie)=phis(:,:,ie) + &
           smooth_phis_nudt*pstens(:,:,ie)/elem(ie)%spheremp(:,:)


        ! remove new extrema.  could use conservative reconstruction from advection
        ! but no reason to conserve mean PHI.
        if (ic < iuse) then
        do i=1,np
        do j=1,np
           if (phis(i,j,ie)>mx) phis(i,j,ie)=mx
           if (phis(i,j,ie)<mn) phis(i,j,ie)=mn
        enddo
        enddo
        endif


        ! user specified minimum
        do i=1,np
        do j=1,np
           if (phis(i,j,ie)<minf) phis(i,j,ie)=minf
        enddo
        enddo

        phis(:,:,ie)=phis(:,:,ie)*elem(ie)%spheremp(:,:)
        call edgeVpack(edge3p1,phis(:,:,ie),1,0,ie)

     enddo

     call t_startf('smooth_phis_bexchV4')
     call bndry_exchangeV(hybrid,edge3p1)
     call t_stopf('smooth_phis_bexchV4')

     do ie=nets,nete
        call edgeVunpack(edge3p1, phis(:,:,ie), 1, 0, ie)
        phis(:,:,ie)=phis(:,:,ie)*elem(ie)%rspheremp(:,:)
     enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
  enddo
  end subroutine smooth_phis

  subroutine overwrite_SEdensity(elem, fvm, dt_q, hybrid,nets,nete, np1)

    use fvm_reconstruction_mod, only: reconstruction
    use fvm_filter_mod, only: monotonic_gradient_cart, recons_val_cart
    use dimensions_mod, only : np, nlev, nc,nhe
    use hybrid_mod, only : hybrid_t
    use edge_mod, only : edgevpack, edgevunpack, edgevunpackmax, edgevunpackmin
    use bndry_mod, only : bndry_exchangev
    use element_mod, only : element_t
    use derivative_mod, only : derivative_t , laplace_sphere_wk
    use time_mod, only : TimeLevel_t
    use fvm_control_volume_mod, only : fvm_struct

    type(element_t) , intent(inout) :: elem(:)
    type(fvm_struct), intent(inout) :: fvm(:)
    type(hybrid_t),   intent(in)    :: hybrid ! distributed parallel structure (shared)
    integer,          intent(in)    :: nets   ! starting thread element number (private)
    integer,          intent(in)    :: nete   ! ending thread element number   (private)
    integer,          intent(in)    :: np1
    integer :: ie, k

    real (kind=real_kind)             :: xp,yp, tmpval, dt_q
    integer                           :: i, j,ix, jy, starti,endi,tmpi

    real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons

    if ((nc .ne. 4) .or. (np .ne. 4)) then
      if(hybrid%masterthread) then
        print *,"You are in OVERWRITE SE AIR DENSITY MODE"
        print *,"This only works for nc=4 and np=4"
        print *,"Write a new search algorithm or pay $10000!"
      endif
      stop
    endif
#if defined(_FVM)
    do ie=nets,nete
      call reconstruction(fvm(ie)%psc, fvm(ie),recons)
      call monotonic_gradient_cart(fvm(ie)%psc, fvm(ie),recons, elem(ie)%desc)
      do j=1,np
        do i=1,np
          xp=tan(elem(ie)%cartp(i,j)%x)
          yp=tan(elem(ie)%cartp(i,j)%y)
          ix=i
          jy=j
          ! Search index along "x"  (bisection method)
!           starti = 1
!           endi = nc+1
!           do
!              if  ((endi-starti) <=  1)  exit
!              tmpi = (endi + starti)/2
!              if (xp  >  fvm%acartx(tmpi)) then
!                 starti = tmpi
!              else
!                 endi = tmpi
!              endif
!           enddo
!           ix = starti
!
!         ! Search index along "y"
!           starti = 1
!           endi = nc+1
!           do
!              if  ((endi-starti) <=  1)  exit
!              tmpi = (endi + starti)/2
!              if (yp  >  fvm%acarty(tmpi)) then
!                 starti = tmpi
!              else
!                 endi = tmpi
!              endif
!           enddo
!           jy = starti

          call recons_val_cart(fvm(ie)%psc, xp,yp,fvm(ie)%spherecentroid,recons,ix,jy,tmpval)
          elem(ie)%state%ps_v(i,j,np1)= elem(ie)%state%ps_v(i,j,np1) +&
               dt_q*(tmpval - elem(ie)%state%ps_v(i,j,np1) )/(7*24*60*60)
        end do
      end do
      elem(ie)%state%ps_v(:,:,np1)=elem(ie)%state%ps_v(:,:,np1)*elem(ie)%spheremp(:,:)
     call edgeVpack(edge3p1,elem(ie)%state%ps_v(:,:,np1),1,0,ie)
  enddo

  call t_startf('overwr_SEdens_bexchV')
  call bndry_exchangeV(hybrid,edge3p1)
  call t_stopf('overwr_SEdens_bexchV')

  do ie=nets,nete
     call edgeVunpack(edge3p1, elem(ie)%state%ps_v(:,:,np1), 1, 0, ie)
     elem(ie)%state%ps_v(:,:,np1)=elem(ie)%state%ps_v(:,:,np1)*elem(ie)%rspheremp(:,:)
  enddo
#endif
  end subroutine overwrite_SEdensity


end module prim_advance_mod

