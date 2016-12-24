#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!#define _DBG_ print *,"File:",__FILE__," at ",__LINE__
!#define _DBG_ !DBG
!
!
module prim_advance_mod

  use control_mod,    only: qsplit,rsplit, use_moisture
  use derivative_mod, only: derivative_t
  use dimensions_mod, only: np, nlev, nlevp, nelemd, qsize, max_corner_elem
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
  public :: prim_advance_exp, prim_advance_init, &
       applyCAMforcing_dynamics, applyCAMforcing, smooth_phis

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
    real(kind=real_kind)  :: eta_dot_dpdn(np,np,nlevp)

    integer :: ie,k,n0,np1

    time  = tl%nstep*dt
    n0    = tl%n0
    np1   = tl%np1

    call set_test_prescribed_wind(elem,deriv,hybrid,hv,dt,tl,nets,nete)
    ! accumulate velocities and fluxes over timesteps
    ! test code only dont bother to openmp thread                                                    
    do ie = nets,nete
       eta_dot_dpdn(:,:,:)=elem(ie)%derived%eta_dot_dpdn_prescribed(:,:,:)
       ! accumulate mean fluxes for advection                                                        
       if (rsplit==0) then
          elem(ie)%derived%eta_dot_dpdn(:,:,:) = &
               elem(ie)%derived%eta_dot_dpdn(:,:,:) + eta_dot_dpdn(:,:,:)*eta_ave_w
       else
          ! lagrangian case.  mean vertical velocity = 0                                             
          elem(ie)%derived%eta_dot_dpdn(:,:,:) = 0
          ! update position of floating levels
          do k=1,nlev
             elem(ie)%state%dp3d(:,:,k,np1) = elem(ie)%state%dp3d(:,:,k,n0)  &
                  + dt*(eta_dot_dpdn(:,:,k+1) - eta_dot_dpdn(:,:,k))
          enddo
       end if
       ! accumulate U*dp 
       do k=1,nlev
          elem(ie)%derived%vn0(:,:,1,k)=elem(ie)%derived%vn0(:,:,1,k)+&
               eta_ave_w*elem(ie)%state%v(:,:,1,k,n0)*elem(ie)%state%dp3d(:,:,k,tl%n0)
          elem(ie)%derived%vn0(:,:,2,k)=elem(ie)%derived%vn0(:,:,2,k)+&
               eta_ave_w*elem(ie)%state%v(:,:,2,k,n0)*elem(ie)%state%dp3d(:,:,k,tl%n0)
       enddo

    enddo
  end subroutine
#endif

  !_____________________________________________________________________
  subroutine prim_advance_exp(elem, deriv, hvcoord, hybrid,dt, tl,  nets, nete, compute_diagnostics)

    use bndry_mod,      only: bndry_exchangev
    use control_mod,    only: prescribed_wind, qsplit, tstep_type, rsplit, qsplit, integration
    use edge_mod,       only: edgevpack, edgevunpack, initEdgeBuffer
    use edgetype_mod,   only: EdgeBuffer_t
    use reduction_mod,  only: reductionbuffer_ordered_1d_t
    use time_mod,       only: timelevel_qdp, tevolve

#ifdef TRILINOS
    use prim_derived_type_mod ,only : derived_type, initialize
    use, intrinsic :: iso_c_binding
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
    call TimeLevel_Qdp(tl, qsplit, qn0)  ! compute current Qdp() timelevel

    ! get time integration method for this timestep

! integration = "explicit"
!
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

! default weights for computing mean dynamics fluxes
    eta_ave_w = 1d0/qsplit

    if (tstep_type==1) then
       method=0                           ! LF
       qsplit_stage = mod(nstep,qsplit)
       if (qsplit_stage==0) method=1      ! RK2 on first of qsplit steps
       ! RK2 + LF scheme has tricky weights:
       eta_ave_w=ur_weights(qsplit_stage+1)
    else
       method = tstep_type                ! other RK variants
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
    if (method<=10) then ! not implicit
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




  subroutine applyCAMforcing(elem,hvcoord,np1,np1_qdp,dt_q,nets,nete)

  use physical_constants, only: Cp

  implicit none
  type (element_t),       intent(inout) :: elem(:)
  real (kind=real_kind),  intent(in)    :: dt_q
  type (hvcoord_t),       intent(in)    :: hvcoord
  integer,                intent(in)    :: np1,nets,nete,np1_qdp

  ! local
  integer :: i,j,k,ie,q
  real (kind=real_kind) :: v1,dp
  real (kind=real_kind) :: beta(np,np),E0(np,np),ED(np,np),dp0m1(np,np),dpsum(np,np)

  do ie=nets,nete
     ! apply forcing to Qdp
     elem(ie)%derived%FQps(:,:,1)=0
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k,i,j,v1)
#endif
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

     if (use_moisture) then
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

  use hybvcoord_mod,  only: hvcoord_t

  implicit none
  type (element_t)     ,  intent(inout) :: elem(:)
  real (kind=real_kind),  intent(in)    :: dt_q
  type (hvcoord_t),       intent(in)    :: hvcoord
  integer,                intent(in)    :: np1,nets,nete

  integer :: i,j,k,ie,q
  real (kind=real_kind) :: v1,dp

  do ie=nets,nete
     elem(ie)%state%T(:,:,:,np1)  = elem(ie)%state%T(:,:,:,np1)    + dt_q*elem(ie)%derived%FT(:,:,:,1)
     elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + dt_q*elem(ie)%derived%FM(:,:,:,:,1)
  enddo
  end subroutine applyCAMforcing_dynamics



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
  use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis, swest
  use hybvcoord_mod, only : hvcoord_t
  use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
  use edge_mod, only : edgevpack, edgevunpack, edgeDGVunpack
  use edgetype_mod, only : EdgeBuffer_t, EdgeDescriptor_t
  use bndry_mod, only : bndry_exchangev
  use viscosity_mod, only : biharmonic_wk_dp3d
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
  real (kind=real_kind) :: nu_scale_top,dpdn0
  integer :: k,kptr,i,j,ie,ic,nt
  real (kind=real_kind), dimension(np,np)      :: dpdn
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ttens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: dptens
  real (kind=real_kind), dimension(0:np+1,0:np+1,nlev)          :: corners
  real (kind=real_kind), dimension(2,2,2)                       :: cflux
  real (kind=real_kind), dimension(np,np,nlev) :: p
  type (EdgeDescriptor_t)                                       :: desc


! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

  real (kind=real_kind), dimension(np,np) :: lap_t,lap_dp
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: v1,v2,dt,heating
  real (kind=real_kind) :: temp(np,np,nlev)



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
        call biharmonic_wk_dp3d(elem,dptens,ttens,vtens,deriv,edge3,hybrid,nt,nets,nete)

        do ie=nets,nete

           ! comptue mean flux
           if (nu_p>0) then
              elem(ie)%derived%dpdiss_ave(:,:,:)=elem(ie)%derived%dpdiss_ave(:,:,:)+&
                   eta_ave_w*elem(ie)%state%dp3d(:,:,:,nt)/hypervis_subcycle
              elem(ie)%derived%dpdiss_biharmonic(:,:,:)=elem(ie)%derived%dpdiss_biharmonic(:,:,:)+&
                   eta_ave_w*dptens(:,:,:,ie)/hypervis_subcycle
           endif
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,lap_t,lap_dp,lap_v,nu_scale_top)
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


              ! biharmonic terms need a negative sign:
              if (nu_top>0 .and. k<=3) then
                 vtens(:,:,:,k,ie)=(-nu*vtens(:,:,:,k,ie) + nu_scale_top*nu_top*lap_v(:,:,:))
                 ttens(:,:,k,ie)  =(-nu_s*ttens(:,:,k,ie) + nu_scale_top*nu_top*lap_t(:,:) )
                 dptens(:,:,k,ie) =(-nu_p*dptens(:,:,k,ie) + nu_scale_top*nu_top*lap_dp(:,:) )
              else
                 vtens(:,:,:,k,ie)=-nu*vtens(:,:,:,k,ie) 
                 ttens(:,:,k,ie)  =-nu_s*ttens(:,:,k,ie) 
                 dptens(:,:,k,ie) =-nu_p*dptens(:,:,k,ie) 
              endif

              if (nu_p==0) then
                 ! nu_p==0 is only for certain regression tests, so perfromance is not an issue
                 ! normalize so as to conserve IE
                 ! scale by 1/rho (normalized to be O(1))
                 ! dp/dn = O(ps0)*O(delta_eta) = O(ps0)/O(nlev)
                 dpdn(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,nt)
                 dpdn0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
                 ttens(:,:,k,ie) = ttens(:,:,k,ie) * dpdn0/dpdn(:,:)
                 dptens(:,:,k,ie) = 0
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
           call edgeVunpack(edge3, elem(ie)%state%dp3d(:,:,:,nt), nlev, kptr, ie)



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
     if (.not. use_moisture ) then
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
              Ephi(i,j)=E+phi(i,j,k)
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
!$omp parallel do private(k)
#endif
     do k=1,nlev
        elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%spheremp(:,:)*( elem(ie)%state%v(:,:,1,k,nm1) + dt2*vtens1(:,:,k) )
        elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%spheremp(:,:)*( elem(ie)%state%v(:,:,2,k,nm1) + dt2*vtens2(:,:,k) )
        elem(ie)%state%T(:,:,k,np1) = elem(ie)%spheremp(:,:)*(elem(ie)%state%T(:,:,k,nm1) + dt2*ttens(:,:,k))
        elem(ie)%state%dp3d(:,:,k,np1) = &
             elem(ie)%spheremp(:,:) * (elem(ie)%state%dp3d(:,:,k,nm1) - &
             dt2 * (divdp(:,:,k) + eta_dot_dpdn(:,:,k+1)-eta_dot_dpdn(:,:,k)))
        
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

     kptr=kptr+2*nlev
     call edgeVunpack(edge3p1, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr,ie)
     
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


  subroutine smooth_phis(phis,elem,hybrid,deriv,nets,nete,minf,numcycle)
  use control_mod, only : smooth_phis_nudt, hypervis_scaling
  use hybrid_mod, only : hybrid_t
  use edge_mod, only : edgevpack, edgevunpack, edgevunpackmax, edgevunpackmin
  use edgetype_mod, only : EdgeBuffer_t
  use bndry_mod, only : bndry_exchangev
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



end module prim_advance_mod

