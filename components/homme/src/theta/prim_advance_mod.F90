#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!#define _DBG_ print *,"File:",__FILE__," at ",__LINE__
!#define _DBG_ !DBG
!
!
module prim_advance_mod

  use control_mod,    only: qsplit,rsplit, use_moisture, theta_hydrostatic_mode
  use derivative_mod, only: derivative_t
  use dimensions_mod, only: np, nlev, nlevp, nelemd, qsize, max_corner_elem
  use edgetype_mod,   only: EdgeDescriptor_t, EdgeBuffer_t
  use element_mod,    only: element_t
  use element_ops,    only: get_pnh_and_exner, set_hydrostatic_phi, get_kappa_star,&
       get_cp_star, get_temperature, set_theta_ref, copy_state
  use hybrid_mod,     only: hybrid_t
  use hybvcoord_mod,  only: hvcoord_t
  use kinds,          only: real_kind, iulog
  use perf_mod,       only: t_startf, t_stopf, t_barrierf, t_adj_detailf ! _EXTERNAL
  use parallel_mod,   only: abortmp, parallel_t, iam
  use time_mod,       only: timelevel_t
  use test_mod,       only: set_prescribed_wind

  implicit none
  private
  save
  public :: prim_advance_exp, prim_advance_init1, &
       applyCAMforcing_dynamics, applyCAMforcing, vertical_mesh_init2

!  type (EdgeBuffer_t) :: edge5
  type (EdgeBuffer_t) :: edge6
  real (kind=real_kind), allocatable :: ur_weights(:)

contains

  subroutine prim_advance_init1(par, elem,integration)
    use edge_mod, only : initEdgeBuffer
    implicit none
    
    type (parallel_t) :: par
    type (element_t), intent(inout), target   :: elem(:)
    character(len=*)    , intent(in) :: integration
    integer :: i
    integer :: ie

!    call initEdgeBuffer(par,edge5,elem,5*nlev)
    call initEdgeBuffer(par,edge6,elem,6*nlev)

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

  end subroutine prim_advance_init1


  subroutine vertical_mesh_init2(elem, nets, nete, hybrid, hvcoord)

    ! additional solver specific initializations (called from prim_init2)

    type (element_t),			intent(inout), target :: elem(:)! array of element_t structures
    integer,				intent(in) :: nets,nete		! start and end element indices
    type (hybrid_t),			intent(in) :: hybrid		! mpi/omp data struct
    type (hvcoord_t),			intent(inout)	:: hvcoord	! hybrid vertical coord data struct

  end subroutine vertical_mesh_init2


  !_____________________________________________________________________
  subroutine prim_advance_exp(elem, deriv, hvcoord, hybrid,dt, tl,  nets, nete, compute_diagnostics)

    use bndry_mod,      only: bndry_exchangev
    use control_mod,    only: prescribed_wind, qsplit, tstep_type, rsplit, qsplit, integration
    use edge_mod,       only: edgevpack, edgevunpack, initEdgeBuffer
    use edgetype_mod,   only: EdgeBuffer_t
    use reduction_mod,  only: reductionbuffer_ordered_1d_t, parallelmax
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
    real (kind=real_kind) ::  itertol,statesave(nets:nete,np,np,nlev,6)
    real (kind=real_kind) ::  statesave2(nets:nete,np,np,nlev,6)
    real (kind=real_kind) ::  statesave3(nets:nete,np,np,nlev,6)
    real (kind=real_kind) ::  itererrmax,gamma,delta
    real (kind=real_kind) ::  fimp_stagevalues(nets:nete,np,np,nlev,6,10)
    real (kind=real_kind) ::  fexp_stagevalues(nets:nete,np,np,nlev,6,10)
    real (kind=real_kind) ::  expvalues(nets:nete,np,np,nlev,6)
 
    integer :: ie,nm1,n0,np1,nstep,method,qsplit_stage,k, qn0
    integer :: n,i,j,lx,lenx,maxiter,itercount

    call t_startf('prim_advance_exp')
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

    ! get timelevel for accessing tracer mass Qdp() to compute virtual temperature
    call TimeLevel_Qdp(tl, qsplit, qn0)  ! compute current Qdp() timelevel

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
          elem(ie)%state%theta_dp_cp(:,:,:,np1)= elem(ie)%state%theta_dp_cp(:,:,:,n0)/3 &
               + 2*elem(ie)%state%theta_dp_cp(:,:,:,np1)/3
          elem(ie)%state%dp3d(:,:,:,np1)= elem(ie)%state%dp3d(:,:,:,n0)/3 &
               + 2*elem(ie)%state%dp3d(:,:,:,np1)/3
          elem(ie)%state%w(:,:,:,np1)= elem(ie)%state%w(:,:,:,n0)/3 &
               + 2*elem(ie)%state%w(:,:,:,np1)/3
          elem(ie)%state%phi(:,:,:,np1)= elem(ie)%state%phi(:,:,:,n0)/3 &
               + 2*elem(ie)%state%phi(:,:,:,np1)/3
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
          elem(ie)%state%theta_dp_cp(:,:,:,nm1)= (5*elem(ie)%state%theta_dp_cp(:,:,:,nm1) &
               - elem(ie)%state%theta_dp_cp(:,:,:,n0) )/4
          elem(ie)%state%dp3d(:,:,:,nm1)= (5*elem(ie)%state%dp3d(:,:,:,nm1) &
                  - elem(ie)%state%dp3d(:,:,:,n0) )/4
          elem(ie)%state%w(:,:,:,nm1)= (5*elem(ie)%state%w(:,:,:,nm1) &
                  - elem(ie)%state%w(:,:,:,n0) )/4
          elem(ie)%state%phi(:,:,:,nm1)= (5*elem(ie)%state%phi(:,:,:,nm1) &
                  - elem(ie)%state%phi(:,:,:,n0) )/4
       enddo
       ! u5 = (5*u1/4 - u0/4) + 3dt/4 RHS(u4)
       call compute_and_apply_rhs(np1,nm1,np1,qn0,3*dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,3*eta_ave_w/4)
       ! final method is the same as:
       ! u5 = u0 +  dt/4 RHS(u0)) + 3dt/4 RHS(u4)
       call t_stopf("U3-5stage_timestep")  
    else if (method==6) then
      call t_startf('implicit Euler')
      itercount=1 
      maxiter=10000
      itertol=1e-1
      itererrmax=2.0*itertol
      do ie=nets,nete
        elem(ie)%state%v(:,:,1,:,np1)   = elem(ie)%state%v(:,:,1,:,n0)
        elem(ie)%state%v(:,:,2,:,np1)   = elem(ie)%state%v(:,:,2,:,n0)
        elem(ie)%state%w(:,:,:,np1)     = elem(ie)%state%w(:,:,:,n0)
        elem(ie)%state%phi(:,:,:,np1)   = elem(ie)%state%phi(:,:,:,n0)
        elem(ie)%state%theta_dp_cp(:,:,:,np1) = elem(ie)%state%theta_dp_cp(:,:,:,n0)
        elem(ie)%state%dp3d(:,:,:,np1)  = elem(ie)%state%dp3d(:,:,:,n0)
      end do
      call fp_iteration_impeuler(elem,np1,n0,nm1,qn0,hvcoord,dt,hybrid,&
            deriv,nets,nete,eta_ave_w,itertol,maxiter,itererrmax,itercount,.false.)
      print*, itercount, itererrmax 
      call t_stopf('implicit Euler')
!================================================================================
    else if (method==7) then ! Imex hevi, implicit euler after the full explicit time-step
    ! it seems to run with ne=16, nlev=26, dt=200 for JW Baro up to 18.8 days at least
    ! Ullrich 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
       ! u1 = u0 + dt/5 RHS(u0)  (save u1 in timelevel nm1)
       call t_startf("U3-5stage_timestep")
       call compute_and_apply_rhs_imex_nonstiff(nm1,n0,n0,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w/4,1.d0,0.d0,1.d0)
       ! u2 = u0 + dt/5 RHS(u1)
       call compute_and_apply_rhs_imex_nonstiff(np1,n0,nm1,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,0.d0,1.d0)
       ! u3 = u0 + dt/3 RHS(u2)
       call compute_and_apply_rhs_imex_nonstiff(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,0.d0,1.d0)
       ! u4 = u0 + 2dt/3 RHS(u3)
       call compute_and_apply_rhs_imex_nonstiff(np1,n0,np1,qn0,2*dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,0.d0,1.d0)
       ! compute (5*u1/4 - u0/4) in timelevel nm1:
       do ie=nets,nete
          elem(ie)%state%v(:,:,:,:,nm1)= (5*elem(ie)%state%v(:,:,:,:,nm1) &
               - elem(ie)%state%v(:,:,:,:,n0) ) /4
          elem(ie)%state%theta_dp_cp(:,:,:,nm1)= (5*elem(ie)%state%theta_dp_cp(:,:,:,nm1) &
               - elem(ie)%state%theta_dp_cp(:,:,:,n0) )/4
          elem(ie)%state%dp3d(:,:,:,nm1)= (5*elem(ie)%state%dp3d(:,:,:,nm1) &
                  - elem(ie)%state%dp3d(:,:,:,n0) )/4
          elem(ie)%state%w(:,:,:,nm1)= (5*elem(ie)%state%w(:,:,:,nm1) &
                  - elem(ie)%state%w(:,:,:,n0) )/4
          elem(ie)%state%phi(:,:,:,nm1)= (5*elem(ie)%state%phi(:,:,:,nm1) &
                  - elem(ie)%state%phi(:,:,:,n0) )/4
       enddo
       ! u5 = (5*u1/4 - u0/4) + 3dt/4 RHS(u4)
       call compute_and_apply_rhs_imex_nonstiff(np1,nm1,np1,qn0,3*dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,3*eta_ave_w/4,1.d0,0.d0,1.d0)
 
      maxiter=1000
      itertol=1e-8
 
      call compute_stage_value_dirk_stiff(np1,n0,n0,qn0,dt,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w,maxiter,itertol)
  !    print *, 'maxiter', maxiter
       call t_stopf("U3-5stage_timestep")
!==========================================================================================
    else if (method==8) then  ! ARS232 from (Ascher et al., 1997), hydrostatic debug
      call t_startf("ARS232_timestep")
      delta = -2.d0*sqrt(2.d0)/3.d0
      gamma = 1.d0 - 1.d0/sqrt(2.d0)

      ! compute g2=un0+dt*gamma*n(un0) and save in unp1
      call compute_and_apply_rhs_imex_nonstiff(np1,n0,n0,qn0,gamma*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w/4,1.d0,0.d0,1.d0)      
     
      ! save un0 as statesave
      call state_save(elem,statesave,n0,nets,nete)

      ! form un0 + dt*delta*n(g1) at save at un0     
      call compute_and_apply_rhs_imex_nonstiff(n0,n0,n0,qn0,delta*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w/4,1.d0,0.d0,1.d0)

      ! compute g3=(un0+dt*delta*n(g1))+dt*(1-delta)*n(g2) and save at unp1
      call compute_and_apply_rhs_imex_nonstiff(n0,n0,np1,qn0,(1-delta)*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w/2,1.d0,0.d0,1.d0)
   
      ! form unp1 = dt*(1-gamma)*n(g2)
      call compute_and_apply_rhs_imex_nonstiff(np1,n0,np1,qn0,(1.d0-gamma)*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w/2,1.d0,0.d0,0.d0)    
      ! form unp1 = dt*(1-gamma)*n(g2)+dt*gamma*n(g3)
      call compute_and_apply_rhs_imex_nonstiff(np1,n0,n0,qn0,(1.d0-gamma)*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w/2,1.d0,0.d0,0.d0)

      ! form unp1 = un0 + dt*(1-gamma)*f(g2)+dt*gamma*f(g3)
      call state_add(elem,statesave,np1,nets,nete,1.d0,1.d0)
    
      maxiter=1000
      itertol=1e-8
      call state_read(elem,statesave,np1,nets,nete)
      call compute_stage_value_dirk_stiff(np1,n0,n0,qn0,dt,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w,maxiter,itertol)

      call state_read(elem,statesave,n0,nets,nete)
      call t_stopf("ARS232_timestep")
!============================================================================================
    else if (method==9) then ! ARS232 from (Ascher et al., 1997), nh-imex
      call t_startf("ARS232_timestep")

      gamma = 1.d0-1.d0/sqrt(2.d0)
      delta = -2.d0*sqrt(2.d0)/3.d0
      
      call state_save(elem,statesave,n0,nets,nete)

      call compute_and_apply_rhs_imex_nonstiff(n0,n0,n0,qn0,gamma*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,1.d0,0.d0,1.d0)

      call state_save(elem,statesave2,n0,nets,nete)
      call state_read(elem,statesave2,np1,nets,nete)

      maxiter = 1000
      itertol = 1e-8
      ! solve for g2
      call compute_stage_value_dirk_stiff(np1,n0,n0,qn0,gamma*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,maxiter,itertol)

     call state_read(elem,statesave,n0,nets,nete)
     call state_save(elem,statesave3,np1,nets,nete)

     call compute_and_apply_rhs_imex_nonstiff(n0,n0,n0,qn0,delta*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,1.d0,0.d0,1.d0)

     call compute_and_apply_rhs_imex_nonstiff(n0,n0,np1,qn0,(1.d0-delta)*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,1.d0,0.d0,1.d0)

     call compute_and_apply_rhs_imex_nonstiff(n0,n0,np1,qn0,(1.d0-gamma)*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,0.d0,1.d0,1.d0)

      maxiter = 1000
      itertol = 1e-8
      ! solve for g3
      call compute_stage_value_dirk_stiff(np1,n0,n0,qn0,gamma*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,maxiter,itertol)

     call state_read(elem,statesave,n0,nets,nete)

     call compute_and_apply_rhs_imex_nonstiff(np1,n0,np1,qn0,gamma*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,1.d0,1.d0,1.d0)

     call state_read(elem,statesave3,n0,nets,nete)

     call compute_and_apply_rhs_imex_nonstiff(np1,np1,n0,qn0,(1.d0-gamma)*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,1.d0,1.d0,1.d0)
     
     call state_read(elem,statesave,n0,nets,nete) 
     
      call t_stopf("ARS232_timestep")
 !================================================================================================
    elseif (method==10) then ! ARK232
      call t_startf("ARK232_timestep")
                        


   
      call t_stopf("ARK232_timestep")
!====================================================================================
    elseif (method==11) then
      call t_startf("imexeuler")
      call state_save(elem,statesave,n0,nets,nete)
      call compute_and_apply_rhs_imex_nonstiff(np1,n0,n0,qn0,dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,1.d0,0.d0,1.d0)
      call state_save(elem,statesave2,np1,nets,nete)
      call state_read(elem,statesave2,n0,nets,nete)
      maxiter = 1000
      itertol = 1e-8
      call compute_stage_value_dirk_stiff(np1,n0,n0,qn0,dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,maxiter,itertol)
      call state_read(elem,statesave,n0,nets,nete)
      call t_stopf("imexeuler")
!======================================================================================
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
          elem(ie)%accum%DIFFTHETA(:,:,:)=elem(ie)%state%theta_dp_cp(:,:,:,np1)
       enddo
    endif
#endif

    ! note:time step computes u(t+1)= u(t*) + RHS.
    ! for consistency, dt_vis = t-1 - t*, so this is timestep method dependent
    ! forward-in-time, hypervis applied to dp3d
    call advance_hypervis(edge6,elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)

#ifdef ENERGY_DIAGNOSTICS
    if (compute_diagnostics) then
       do ie = nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
         do k=1,nlev  !  Loop index added (AAM)
          elem(ie)%accum%DIFF(:,:,:,k)=( elem(ie)%state%v(:,:,:,k,np1) -&
               elem(ie)%accum%DIFF(:,:,:,k) ) / dt_vis
          elem(ie)%accum%DIFFTHETA(:,:,k)=( elem(ie)%state%theta_dp_cp(:,:,k,np1) -&
               elem(ie)%accum%DIFFTHETA(:,:,k) ) / dt_vis
         enddo
       enddo
    endif
#endif

    tevolve=tevolve+dt

    call t_stopf('prim_advance_exp')
    end subroutine prim_advance_exp




  subroutine applyCAMforcing(elem,hvcoord,np1,np1_qdp,dt,nets,nete)

  implicit none
  type (element_t),       intent(inout) :: elem(:)
  real (kind=real_kind),  intent(in)    :: dt
  type (hvcoord_t),       intent(in)    :: hvcoord
  integer,                intent(in)    :: np1,nets,nete,np1_qdp

  ! local
  integer :: i,j,k,ie,q
  real (kind=real_kind) :: v1
  real (kind=real_kind) :: temperature(np,np,nlev)
  real (kind=real_kind) :: kappa_star(np,np,nlev)
  real (kind=real_kind) :: cp_star(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev)
  real (kind=real_kind) :: dpnh(np,np,nlev)

  do ie=nets,nete
     ! apply forcing to Qdp
     elem(ie)%derived%FQps(:,:)=0

     ! apply forcing to temperature
     call get_temperature(elem(ie),temperature,hvcoord,np1,np1_qdp)
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        temperature(:,:,k) = temperature(:,:,k) + elem(ie)%derived%FT(:,:,k)
     enddo

     
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k,i,j,v1)
#endif
     do q=1,qsize
        do k=1,nlev
           do j=1,np
              do i=1,np
                 v1 = dt*elem(ie)%derived%FQ(i,j,k,q)
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
                    elem(ie)%derived%FQps(i,j)=elem(ie)%derived%FQps(i,j)+v1/dt
                 endif
              enddo
           enddo
        enddo
     enddo

     if (use_moisture) then
        ! to conserve dry mass in the precese of Q1 forcing:
        elem(ie)%state%ps_v(:,:,np1) = elem(ie)%state%ps_v(:,:,np1) + &
             dt*elem(ie)%derived%FQps(:,:)
     endif


     ! Qdp(np1) and ps_v(np1) were updated by forcing - update Q(np1)
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k)
#endif
     do k=1,nlev
        dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
             ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
     enddo
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k)
#endif
     do q=1,qsize
        do k=1,nlev
           elem(ie)%state%Q(:,:,k,q) = elem(ie)%state%Qdp(:,:,k,q,np1_qdp)/dp(:,:,k)
        enddo
     enddo

     ! now that we have updated Qdp and dp, compute theta_dp_cp from temperature
     call get_kappa_star(kappa_star,elem(ie)%state%Qdp(:,:,:,1,np1_qdp),dp)
     call get_cp_star(cp_star,elem(ie)%state%Qdp(:,:,:,1,np1_qdp),dp)
     call get_pnh_and_exner(hvcoord,elem(ie)%state%theta_dp_cp(:,:,:,np1),dp,&
          elem(ie)%state%phi(:,:,:,np1),elem(ie)%state%phis(:,:),kappa_star,&
          pnh,dpnh,exner)

     elem(ie)%state%theta_dp_cp(:,:,:,np1) = temperature(:,:,:)*cp_star(:,:,:)&
          *dp(:,:,:)/exner(:,:,:)

  enddo
  call applyCAMforcing_dynamics(elem,hvcoord,np1,np1_qdp,dt,nets,nete)
  end subroutine applyCAMforcing



  subroutine applyCAMforcing_dynamics(elem,hvcoord,np1,np1_qdp,dt,nets,nete)
  use physical_constants, only: Cp
  use hybvcoord_mod,  only: hvcoord_t

  implicit none
  type (element_t)     ,  intent(inout) :: elem(:)
  real (kind=real_kind),  intent(in)    :: dt
  type (hvcoord_t),       intent(in)    :: hvcoord
  integer,                intent(in)    :: np1,nets,nete,np1_qdp

  integer :: k,ie

  do ie=nets,nete
     elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + dt*elem(ie)%derived%FM(:,:,1:2,:)
     elem(ie)%state%w(:,:,:,np1) = elem(ie)%state%w(:,:,:,np1) + dt*elem(ie)%derived%FM(:,:,3,:)
  enddo
  end subroutine applyCAMforcing_dynamics



  subroutine advance_hypervis(edgebuf,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
  !
  !  take one timestep of:
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use control_mod, only : nu, nu_div, hypervis_order, hypervis_subcycle, nu_s, nu_p, nu_top, psurf_vis, swest
  use hybvcoord_mod, only : hvcoord_t
  use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
  use edge_mod, only : edgevpack, edgevunpack, edgeDGVunpack
  use edgetype_mod, only : EdgeBuffer_t, EdgeDescriptor_t
  use bndry_mod, only : bndry_exchangev
  use viscosity_theta, only : biharmonic_wk_theta
  use physical_constants, only: Cp,p0,kappa,g
  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edgebuf
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord

  real (kind=real_kind) :: dt2
  integer :: nets,nete

  ! local
  real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
  real (kind=real_kind) :: nu_scale_top
  integer :: k,kptr,i,j,ie,ic,nt
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,4,nets:nete)      :: stens  ! dp3d,theta,w,phi
!  real (kind=real_kind), dimension(np,np,nlev) :: p


! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

  real (kind=real_kind), dimension(np,np,4) :: lap_s  ! dp3d,theta,w,phi
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: v1(np,np),v2(np,np),heating(np,np)
  real (kind=real_kind) :: dt
  real (kind=real_kind) :: ps_ref(np,np)

  real (kind=real_kind) :: theta_ref(np,np,nlev,nets:nete)
  real (kind=real_kind) :: phi_ref(np,np,nlev,nets:nete)
  real (kind=real_kind) :: dp_ref(np,np,nlev,nets:nete)



  if (nu == 0 .and. nu_p==0 ) return;
  call t_startf('advance_hypervis')


  dt=dt2/hypervis_subcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  regular viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 1) then
     call abortmp( 'ERROR: hypervis_order == 1 not coded')
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTE1:  Diffusion works best when applied to theta.
! It creates some TOM noise when applied to theta_dp_cp in DCMIP 2.0 test
! so we convert from theta_dp_cp->theta, and then convert back at the end of diffusion
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute reference states
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie=nets,nete
     ps_ref(:,:) = sum(elem(ie)%state%dp3d(:,:,:,nt),3)
     do k=1,nlev
        dp_ref(:,:,k,ie) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
             (hvcoord%hybi(k+1)-hvcoord%hybi(k))*ps_ref(:,:)
     enddo

     call set_hydrostatic_phi(hvcoord,elem(ie)%state%phis,&
          elem(ie)%state%theta_dp_cp(:,:,:,nt),elem(ie)%state%dp3d(:,:,:,nt),&
          phi_ref(:,:,:,ie))

     do k=1,nlev
        ! convert theta_dp_cp -> theta
        elem(ie)%state%theta_dp_cp(:,:,k,nt)=&
             elem(ie)%state%theta_dp_cp(:,:,k,nt)/(Cp*elem(ie)%state%dp3d(:,:,k,nt))
     enddo

     call set_theta_ref(hvcoord,elem(ie)%state%dp3d(:,:,:,nt),theta_ref(:,:,:,ie))
#if 0
     theta_ref(:,:,:,ie)=0
     phi_ref(:,:,:,ie)=0
     dp_ref(:,:,:,ie)=0
#endif              
  enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 2) then
     do ic=1,hypervis_subcycle
        do ie=nets,nete

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
           do k=1,nlev
              elem(ie)%state%theta_dp_cp(:,:,k,nt)=elem(ie)%state%theta_dp_cp(:,:,k,nt)-&
                   theta_ref(:,:,k,ie)
              elem(ie)%state%phi(:,:,k,nt)=elem(ie)%state%phi(:,:,k,nt)-&
                   phi_ref(:,:,k,ie)
              elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt)-&
                   dp_ref(:,:,k,ie)
           enddo
        enddo

        call biharmonic_wk_theta(elem,stens,vtens,deriv,edge6,hybrid,nt,nets,nete)

        do ie=nets,nete

           ! comptue mean flux
           if (nu_p>0) then
              elem(ie)%derived%dpdiss_ave(:,:,:)=elem(ie)%derived%dpdiss_ave(:,:,:)+&
                   eta_ave_w*elem(ie)%state%dp3d(:,:,:,nt)/hypervis_subcycle
              elem(ie)%derived%dpdiss_biharmonic(:,:,:)=elem(ie)%derived%dpdiss_biharmonic(:,:,:)+&
                   eta_ave_w*stens(:,:,:,1,ie)/hypervis_subcycle
           endif
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,lap_s,lap_v,nu_scale_top)
#endif
           do k=1,nlev
              ! advace in time.
              ! note: DSS commutes with time stepping, so we can time advance and then DSS.
              ! note: weak operators alreayd have mass matrix "included"

              ! add regular diffusion in top 3 layers:
              if (nu_top>0 .and. k<=3) then
                 lap_s(:,:,1)=laplace_sphere_wk(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_s(:,:,2)=laplace_sphere_wk(elem(ie)%state%theta_dp_cp(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_s(:,:,3)=laplace_sphere_wk(elem(ie)%state%w(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_s(:,:,4)=laplace_sphere_wk(elem(ie)%state%phi(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              endif
              nu_scale_top = 1
              if (k==1) nu_scale_top=4
              if (k==2) nu_scale_top=2


              ! biharmonic terms need a negative sign:
              if (nu_top>0 .and. k<=3) then
                 vtens(:,:,:,k,ie)=(  -nu*vtens(:,:,:,k,ie) + nu_scale_top*nu_top*lap_v(:,:,:))
                 stens(:,:,k,1,ie)=(-nu_p*stens(:,:,k,1,ie) + nu_scale_top*nu_top*lap_s(:,:,1)) ! dp3d
                 stens(:,:,k,2,ie)=(  -nu*stens(:,:,k,2,ie) + nu_scale_top*nu_top*lap_s(:,:,2)) ! theta
                 stens(:,:,k,3,ie)=(  -nu*stens(:,:,k,3,ie) + nu_scale_top*nu_top*lap_s(:,:,3)) ! w
                 stens(:,:,k,4,ie)=(-nu_s*stens(:,:,k,4,ie) + nu_scale_top*nu_top*lap_s(:,:,4)) ! w
              else
                 vtens(:,:,:,k,ie)=  -nu*vtens(:,:,:,k,ie)
                 stens(:,:,k,1,ie)=-nu_p*stens(:,:,k,1,ie)
                 stens(:,:,k,2,ie)=  -nu*stens(:,:,k,2,ie)
                 stens(:,:,k,3,ie)=  -nu*stens(:,:,k,3,ie)
                 stens(:,:,k,4,ie)=-nu_s*stens(:,:,k,4,ie)
              endif
           enddo


           kptr=0
           call edgeVpack(edgebuf,vtens(:,:,:,:,ie),2*nlev,kptr,ie)
           kptr=2*nlev
           call edgeVpack(edgebuf,stens(:,:,:,:,ie),4*nlev,kptr,ie)
        enddo

        call t_startf('ahdp_bexchV2')
        call bndry_exchangeV(hybrid,edgebuf)
        call t_stopf('ahdp_bexchV2')

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edgebuf, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
           kptr=2*nlev
           call edgeVunpack(edgebuf, stens(:,:,:,:,ie), 4*nlev, kptr, ie)


           ! apply inverse mass matrix, accumulate tendencies
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
           do k=1,nlev
              vtens(:,:,1,k,ie)=dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
              vtens(:,:,2,k,ie)=dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
              stens(:,:,k,1,ie)=dt*stens(:,:,k,1,ie)*elem(ie)%rspheremp(:,:)  ! dp3d
              stens(:,:,k,2,ie)=dt*stens(:,:,k,2,ie)*elem(ie)%rspheremp(:,:)  ! theta
              stens(:,:,k,3,ie)=dt*stens(:,:,k,3,ie)*elem(ie)%rspheremp(:,:)  ! w
              stens(:,:,k,4,ie)=dt*stens(:,:,k,4,ie)*elem(ie)%rspheremp(:,:)  ! phi

              !add ref state back
              elem(ie)%state%theta_dp_cp(:,:,k,nt)=elem(ie)%state%theta_dp_cp(:,:,k,nt)+&
                   theta_ref(:,:,k,ie)
              elem(ie)%state%phi(:,:,k,nt)=elem(ie)%state%phi(:,:,k,nt)+&
                   phi_ref(:,:,k,ie)
              elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt)+&
                   dp_ref(:,:,k,ie)

           enddo



#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,v1,v2,heating)
#endif
           do k=1,nlev
              ! update v first (gives better results than updating v after heating)
              elem(ie)%state%v(:,:,:,k,nt)=elem(ie)%state%v(:,:,:,k,nt) + &
                   vtens(:,:,:,k,ie)
              elem(ie)%state%w(:,:,k,nt)=elem(ie)%state%w(:,:,k,nt) &
                   +stens(:,:,k,3,ie)
              !v1=elem(ie)%state%v(:,:,1,k,nt)
              !v2=elem(ie)%state%v(:,:,2,k,nt)
              !                   For the 3D non-hydrostatic with hypervisosity in the horizontal components, 
              !                   the heating term to be added to theta is (dpi/ds)*HVterm/ p^kappa as opposed to 
              !                   HVterm/c_p^*
              !                    
              ! commenting out for now, have to figure out how to get exner in this routine
              !                   Form u*nu div^H(u)
              !                    heating = ( (vtens(:,:,1,k,ie)*v1  + vtens(:,:,2,k,ie)*v2 )
              !                    heating = dpnh(:,:,k)*heating/exner(:,:,k)
              heating(:,:) = 0
              
              elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt) &
                   +stens(:,:,k,1,ie)
              
              elem(ie)%state%theta_dp_cp(:,:,k,nt)=elem(ie)%state%theta_dp_cp(:,:,k,nt) &
                   +stens(:,:,k,2,ie)-heating(:,:)
              
              elem(ie)%state%w(:,:,k,nt)=elem(ie)%state%w(:,:,k,nt) &
                   +stens(:,:,k,3,ie)
              
              elem(ie)%state%phi(:,:,k,nt)=elem(ie)%state%phi(:,:,k,nt) &
                   +stens(:,:,k,4,ie)
           enddo
        enddo
     enddo
  endif

  ! convert theta_dp_cp -> theta
  do ie=nets,nete            
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,v1,v2,heating)
#endif
     do k=1,nlev
        elem(ie)%state%theta_dp_cp(:,:,k,nt)=&
             elem(ie)%state%theta_dp_cp(:,:,k,nt)*Cp*elem(ie)%state%dp3d(:,:,k,nt)
     enddo
  enddo


  call t_stopf('advance_hypervis')

  end subroutine advance_hypervis






  subroutine compute_and_apply_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w)
  ! ===================================
  ! compute the RHS, accumulate into u(np1) and apply DSS
  !
  !           u(np1) = u(nm1) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! This subroutine was orgininally called to compute a leapfrog timestep
  ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
  ! accomodated.  For example, setting nm1=np1=n0 this routine will
  ! take a forward euler step, overwriting the input with the output.
  !
  !    qn0 = timelevel used to access Qdp() in order to compute virtual Temperature
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

  use physical_constants, only : cp, cpwater_vapor, Rgas, kappa, Rwater_vapor,p0, g
  use physics_mod, only : virtual_specific_heat, virtual_temperature
  use prim_si_mod, only : preq_vertadv_v, preq_vertadv_upwind, preq_omega_ps, preq_hydrostatic_v2

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
  real (kind=real_kind), pointer, dimension(:,:,:)   :: dp3d
  real (kind=real_kind), pointer, dimension(:,:,:)   :: theta_dp_cp
   
  real (kind=real_kind) :: kappa_star(np,np,nlev)
  real (kind=real_kind) :: theta_cp(np,np,nlev)
  real (kind=real_kind) :: omega_p(np,np,nlev)
  real (kind=real_kind) :: vort(np,np,nlev)      ! vorticity
  real (kind=real_kind) :: divdp(np,np,nlev)     
  real (kind=real_kind) :: pnh(np,np,nlev)     ! nh (nonydro) pressure
  real (kind=real_kind) :: dpnh(np,np,nlev)    ! 
  real (kind=real_kind) :: exner(np,np,nlev)     ! exner nh pressure
  real (kind=real_kind) :: dpnh_dp(np,np,nlev)   ! dpnh / dp3d  
  real (kind=real_kind) :: eta_dot_dpdn(np,np,nlevp)  ! vertical velocity at interfaces
  real (kind=real_kind) :: KE(np,np,nlev)           ! Kinetic energy
  real (kind=real_kind) :: gradexner(np,np,2,nlev)  ! grad(p^kappa)   
  real (kind=real_kind) :: gradphi(np,np,2,nlev)     
  real (kind=real_kind) :: gradKE(np,np,2,nlev)  ! grad(0.5 u^T u )
  
  real (kind=real_kind) :: v_gradw(np,np,nlev)     
  real (kind=real_kind) :: v_gradtheta(np,np,nlev)     
  real (kind=real_kind) :: v_theta(np,np,2,nlev)
  real (kind=real_kind) :: div_v_theta(np,np,nlev)
  real (kind=real_kind) :: v_gradphi(np,np,nlev)
  real (kind=real_kind) :: v_gradKE(np,np,nlev)     
  real (kind=real_kind) :: vdp(np,np,2,nlev)

  real (kind=real_kind) :: vtens1(np,np,nlev)
  real (kind=real_kind) :: vtens2(np,np,nlev)
  real (kind=real_kind) :: v_vadv(np,np,2,nlev)   ! velocity vertical advection
  real (kind=real_kind) :: s_state(np,np,nlev,2)  ! scalars w,theta,phi
  real (kind=real_kind) :: s_vadv(np,np,nlev,3)   ! scalar vertical advection 
  real (kind=real_kind) :: s_theta_dp_cpadv(np,np,nlev)
  real (kind=real_kind) :: stens(np,np,nlev,3)    ! tendencies w,theta,phi

  real (kind=real_kind) ::  temp(np,np,nlev)
  real (kind=real_kind), dimension(np,np)      :: sdot_sum   ! temporary field
  real (kind=real_kind), dimension(np,np,2)    :: vtemp     ! generic gradient storage
  real (kind=real_kind) ::  v1,v2,w,d_eta_dot_dpdn_dn
  integer :: i,j,k,kptr,ie

  call t_startf('compute_and_apply_rhs')
  do ie=nets,nete
     dp3d  => elem(ie)%state%dp3d(:,:,:,n0)
     theta_dp_cp  => elem(ie)%state%theta_dp_cp(:,:,:,n0)
     theta_cp(:,:,:) = theta_dp_cp(:,:,:)/dp3d(:,:,:)
     phi => elem(ie)%state%phi(:,:,:,n0)

     call get_kappa_star(kappa_star,elem(ie)%state%Qdp(:,:,:,1,qn0),dp3d)

     call get_pnh_and_exner(hvcoord,theta_dp_cp,dp3d,phi,elem(ie)%state%phis,&
             kappa_star,pnh,dpnh,exner) ! ,exner_i)
   
     if (theta_hydrostatic_mode) then
        ! Compute Hydrostatic equation
        do k=1,nlev
           temp(:,:,k) = kappa_star(:,:,k)*theta_dp_cp(:,:,k)*exner(:,:,k)/pnh(:,:,k)
        enddo
        !call preq_hydrostatic(phi,elem(ie)%state%phis,T_v,p,dp)
        call preq_hydrostatic_v2(phi,elem(ie)%state%phis,temp)
        dpnh_dp(:,:,:) = 1
     else
        ! d(p-nh) / d(p-hyrdostatic)
        dpnh_dp(:,:,:) = dpnh(:,:,:)/dp3d(:,:,:)
     endif
     
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        vdp(:,:,1,k) = elem(ie)%state%v(:,:,1,k,n0)*dp3d(:,:,k)
        vdp(:,:,2,k) = elem(ie)%state%v(:,:,2,k,n0)*dp3d(:,:,k)

        ! ================================
        ! Accumulate mean Vel_rho flux in vn0
        ! ================================
        elem(ie)%derived%vn0(:,:,:,k)=elem(ie)%derived%vn0(:,:,:,k)+eta_ave_w*vdp(:,:,:,k)

        divdp(:,:,k)=divergence_sphere(vdp(:,:,:,k),deriv,elem(ie))
        vort(:,:,k)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,n0),deriv,elem(ie))

     enddo


     ! Compute omega_p according to CCM-3
     !call preq_omega_ps(omega_p,hvcoord,p,vgrad_p,divdp)
     ! how will we compute this?  omega_p = 1/p Dp/Dt
     omega_p = 0


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
        s_vadv=0
        s_theta_dp_cpadv=0
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
        ! TODO: remove theta from s_state and s_vadv
        s_state(:,:,:,1)=elem(ie)%state%w(:,:,:,n0)
        s_state(:,:,:,2)=elem(ie)%state%phi(:,:,:,n0)
        call preq_vertadv_v(elem(ie)%state%v(:,:,:,:,n0),s_state,2,eta_dot_dpdn,dp3d,v_vadv,s_vadv)
        !call preq_vertadv_v(elem(ie)%state%v(:,:,:,:,n0),s_state,3,eta_dot_dpdn,dp3d,v_vadv,s_vadv)
        !call preq_vertadv_upwind(elem(ie)%state%v(:,:,:,:,n0),s_state,3,eta_dot_dpdn,dp3d,v_vadv,s_vadv)

        !    this loop constructs d( eta-dot * theta_dp_cp)/deta
        !   d( eta_dot_dpdn * theta*cp)
        !  so we need to compute theta_cp form theta_dp_cp and average to interfaces
        k=1
        s_theta_dp_cpadv(:,:,k)= &
             eta_dot_dpdn(:,:,k+1)* (theta_cp(:,:,k+1)+theta_cp(:,:,k))/2  

        k=nlev
        s_theta_dp_cpadv(:,:,k)= - &
             eta_dot_dpdn(:,:,k)  * (theta_cp(:,:,k)+theta_cp(:,:,k-1))/2  

        do k=2,nlev-1
           s_theta_dp_cpadv(:,:,k)= &
              eta_dot_dpdn(:,:,k+1)* (theta_cp(:,:,k+1)+theta_cp(:,:,k))/2  - &
              eta_dot_dpdn(:,:,k)  * (theta_cp(:,:,k)+theta_cp(:,:,k-1))/2  
        end do
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
!$omp parallel do private(k,i,j,v1,v2,KE,vtemp)
#endif

     vertloop: do k=1,nlev

        ! ================================================
        ! w,theta,phi tendencies:
        ! ================================================

        vtemp(:,:,:)   = gradient_sphere(elem(ie)%state%w(:,:,k,n0),deriv,elem(ie)%Dinv)
        v_gradw(:,:,k) = elem(ie)%state%v(:,:,1,k,n0)*vtemp(:,:,1) &
             +elem(ie)%state%v(:,:,2,k,n0)*vtemp(:,:,2) 
        stens(:,:,k,1) = -s_vadv(:,:,k,1) - v_gradw(:,:,k)  -  g*(1-dpnh_dp(:,:,k) )
        v_theta(:,:,1,k) = elem(ie)%state%v(:,:,1,k,n0)*               &
          elem(ie)%state%theta_dp_cp(:,:,k,n0)
        v_theta(:,:,2,k) =                                             &
          elem(ie)%state%v(:,:,2,k,n0)                                 &
          *elem(ie)%state%theta_dp_cp(:,:,k,n0)
        div_v_theta(:,:,k)=divergence_sphere(v_theta(:,:,:,k),deriv,elem(ie))
        stens(:,:,k,2)=-s_theta_dp_cpadv(:,:,k)-div_v_theta(:,:,k)

        gradphi(:,:,:,k) = gradient_sphere(phi(:,:,k),deriv,elem(ie)%Dinv)

        v_gradphi(:,:,k) = elem(ie)%state%v(:,:,1,k,n0)*gradphi(:,:,1,k) &
             +elem(ie)%state%v(:,:,2,k,n0)*gradphi(:,:,2,k) 
        ! use of s_vadv(:,:,k,2) here is correct since this corresponds to etadot d(phi)/deta
        stens(:,:,k,3) =  -s_vadv(:,:,k,2) - v_gradphi(:,:,k) + g*elem(ie)%state%w(:,:,k,n0)

        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              KE(i,j,k)=0.5D0*( v1*v1 + v2*v2 )
           end do
        end do
        gradKE(:,:,:,k) = gradient_sphere(KE(:,:,k),deriv,elem(ie)%Dinv)
  
     gradexner(:,:,:,k) = gradient_sphere(exner(:,:,k),deriv,elem(ie)%Dinv)        

        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              v_gradKE(i,j,k)=v1*gradKE(i,j,1,k)+v2*gradKE(i,j,2,k)

              vtens1(i,j,k) = -v_vadv(i,j,1,k) &
                   + v2*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - gradKE(i,j,1,k) -gradphi(i,j,1,k)*dpnh_dp(i,j,k) &
                   -theta_cp(i,j,k)*gradexner(i,j,1,k)

              vtens2(i,j,k) = -v_vadv(i,j,2,k) &
                   - v1*(elem(ie)%fcor(i,j) + vort(i,j,k)) &
                   - gradKE(i,j,2,k) -gradphi(i,j,2,k)*dpnh_dp(i,j,k) &
                   -theta_cp(i,j,k)*gradexner(i,j,2,k)
           end do
        end do
     end do vertloop


     
#ifdef ENERGY_DIAGNOSTICS
     ! =========================================================
     !
     ! diagnostics
     ! recomputes some gradients that were not saved above
     ! uses:  sdot_sum(), eta_dot_dpdn(), grad_ps()
     ! grad_phi(), dp3d(), p(), theta_vadv(), v_vadv(), divdp()
     ! =========================================================
 
     ! =========================================================
     ! (AAM) - This section has accumulations over vertical levels.
     !   Be careful if implementing OpenMP
     ! =========================================================

     if (compute_diagnostics) then
        elem(ie)%accum%KEhoriz1=0
        elem(ie)%accum%PEhoriz1=0
        elem(ie)%accum%PEhoriz2=0
        elem(ie)%accum%KE1=0
        elem(ie)%accum%KE2=0
        elem(ie)%accum%KEvert1=0
        elem(ie)%accum%KEvert2=0
        elem(ie)%accum%IEvert1=0
        elem(ie)%accum%IEvert2=0
        elem(ie)%accum%PEvert1=0
        elem(ie)%accum%PEvert2=0
        elem(ie)%accum%T01=0
        elem(ie)%accum%T2=0
        elem(ie)%accum%S1=0
        elem(ie)%accum%S2=0
        elem(ie)%accum%P1=0
        elem(ie)%accum%P2=0
        ! See element_state.F90 for an account of what these variables are defined as
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,v1,v2,vtemp,KE,d_eta_dot_dpdn_dn)
#endif
        do k =1,nlev
          do j=1,np
            do i=1,np                
                  d_eta_dot_dpdn_dn=(eta_dot_dpdn(i,j,k+1)-                    &
                  eta_dot_dpdn(i,j,k))
               !  Form KEhoriz1
                  elem(ie)%accum%KEhoriz1(i,j)=elem(ie)%accum%KEhoriz1(i,j)    &
                  -v_gradKE(i,j,k)*dp3d(i,j,k) - KE(i,j,k)*divdp(i,j,k)
                  elem(ie)%accum%KE1(i,j)=elem(ie)%accum%KE1(i,j)              &
                  -v_gradKE(i,j,k)*dp3d(i,j,k) 
                  elem(ie)%accum%KE2(i,j)=elem(ie)%accum%KE2(i,j)              &
                  -KE(i,j,k)*divdp(i,j,k)
               !  Form KEhoriz2
                  elem(ie)%accum%KEhoriz2(i,j)=elem(ie)%accum%KEhoriz2(i,j)-   &
                  dp3d(i,j,k) * elem(ie)%state%w(i,j,k,n0) * v_gradw(i,j,k)    &
                  -0.5*(elem(ie)%state%w(i,j,k,n0))**2 * divdp(i,j,k) 
               !  Form KEvert1
                  elem(ie)%accum%KEvert1(i,j)=elem(ie)%accum%KEvert1(i,j)-     &
                  (elem(ie)%state%v(i,j,1,k,n0) * v_vadv(i,j,1,k) +            &
                  elem(ie)%state%v(i,j,2,k,n0) *v_vadv(i,j,2,k))*dp3d(i,j,k)-  &      
                  0.5*((elem(ie)%state%v(i,j,1,k,n0))**2 +                     &
                       (elem(ie)%state%v(i,j,2,k,n0))**2)*d_eta_dot_dpdn_dn
               !  Form KEvert2
                  elem(ie)%accum%KEvert2(i,j)=elem(ie)%accum%KEvert2(i,j)      &
                  -s_vadv(i,j,k,1)*(elem(ie)%state%w(i,j,k,n0))*dp3d(i,j,k)    &   
                  -d_eta_dot_dpdn_dn*(elem(ie)%state%w(i,j,k,n0))**2
               !  Form IEvert1
                  elem(ie)%accum%IEvert1(i,j)=elem(ie)%accum%IEvert1(i,j)      &
                  -exner(i,j,k)*s_theta_dp_cpadv(i,j,k)                        
               !  Form IEvert2
                  elem(ie)%accum%IEvert2(i,j)=elem(ie)%accum%IEvert2(i,j)      &
                  +dpnh(i,j,k)*s_vadv(i,j,k,2)
               !  Form PEhoriz1
                  elem(ie)%accum%PEhoriz1(i,j)=(elem(ie)%accum%PEhoriz1(i,j))  &
                  -phi(i,j,k)*divdp(i,j,k) 
               !  Form PEhoriz2                                                &
                  elem(ie)%accum%PEhoriz2(i,j)=elem(ie)%accum%PEhoriz2(i,j)    &
                  -dp3d(i,j,k)*v_gradphi(i,j,k)      
               !  Form PEvert1
                  elem(ie)%accum%PEvert1(i,j) = (elem(ie)%accum%PEvert1(i,j))   &
                  -phi(i,j,k)*d_eta_dot_dpdn_dn                                 
               !  Form PEvert2
                  elem(ie)%accum%PEvert2(i,j) = elem(ie)%accum%PEvert2(i,j)     &
                  -dp3d(i,j,k)*s_theta_dp_cpadv(i,j,k)
               !  Form T01
                  elem(ie)%accum%T01(i,j)=elem(ie)%accum%T01(i,j)               &
                  -(elem(ie)%state%theta_dp_cp(i,j,k,n0))                       &
                  *(gradexner(i,j,1,k)*elem(ie)%state%v(i,j,1,k,n0) +           &
                  gradexner(i,j,2,k)*elem(ie)%state%v(i,j,2,k,n0))              
               !  Form S1 
                  elem(ie)%accum%S1(i,j)=elem(ie)%accum%S1(i,j)                 &
                  -exner(i,j,k)*div_v_theta(i,j,k)
               !  Form T2 
                  elem(ie)%accum%T2(i,j)=elem(ie)%accum%T2(i,j)+                & 
                  (g*(elem(ie)%state%w(i,j,k,n0))-                              &
                  v_gradphi(i,j,k))*dpnh(i,j,k)                                 
               !  Form S2
                  elem(ie)%accum%S2(i,j)=elem(ie)%accum%S2(i,j)                 &
                  -g*(elem(ie)%state%w(i,j,k,n0))+v_gradphi(i,j,k)              &
                  *dpnh(i,j,k)
               !  Form P1
                  elem(ie)%accum%P1(i,j)=elem(ie)%accum%P1(i,j)                 &
                  -g*(elem(ie)%state%w(i,j,k,n0)) * dp3d(i,j,k)
               !  Form P2
                  elem(ie)%accum%P2(i,j)=elem(ie)%accum%P2(i,j)                 &
                  + g * (elem(ie)%state%w(i,j,k,n0)) * dp3d(i,j,k)
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

        elem(ie)%state%w(:,:,k,np1)     = elem(ie)%spheremp(:,:)*(elem(ie)%state%w(:,:,k,nm1)     + dt2*stens(:,:,k,1))
        elem(ie)%state%theta_dp_cp(:,:,k,np1) = elem(ie)%spheremp(:,:)*(elem(ie)%state%theta_dp_cp(:,:,k,nm1) + dt2*stens(:,:,k,2))
        elem(ie)%state%phi(:,:,k,np1)   = elem(ie)%spheremp(:,:)*(elem(ie)%state%phi(:,:,k,nm1)   + dt2*stens(:,:,k,3))

        elem(ie)%state%dp3d(:,:,k,np1) = &
             elem(ie)%spheremp(:,:) * (elem(ie)%state%dp3d(:,:,k,nm1) - &
             dt2 * (divdp(:,:,k) + eta_dot_dpdn(:,:,k+1)-eta_dot_dpdn(:,:,k)))
        
     enddo

     kptr=0
     call edgeVpack(edge6, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr, ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%theta_dp_cp(:,:,:,np1),nlev,kptr,ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%w(:,:,:,np1),nlev,kptr,ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%phi(:,:,:,np1),nlev,kptr,ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,ie)
  end do

  call t_startf('caar_bexchV')
  call bndry_exchangeV(hybrid,edge6)
  call t_stopf('caar_bexchV')

  do ie=nets,nete
     kptr=0
     call edgeVunpack(edge6, elem(ie)%state%dp3d(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%theta_dp_cp(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%w(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%phi(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, ie)

     
     ! ====================================================
     ! Scale tendencies by inverse mass matrix
     ! ====================================================
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        elem(ie)%state%dp3d(:,:,k,np1) =elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
        elem(ie)%state%theta_dp_cp(:,:,k,np1)=elem(ie)%rspheremp(:,:)*elem(ie)%state%theta_dp_cp(:,:,k,np1)
        elem(ie)%state%w(:,:,k,np1)    =elem(ie)%rspheremp(:,:)*elem(ie)%state%w(:,:,k,np1)
        elem(ie)%state%phi(:,:,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%phi(:,:,k,np1)
        elem(ie)%state%v(:,:,1,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,1,k,np1)
        elem(ie)%state%v(:,:,2,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,2,k,np1)
     end do
  end do

  call t_stopf('compute_and_apply_rhs')

  end subroutine compute_and_apply_rhs



  subroutine compute_and_apply_rhs_imex(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w,stiff)
  ! ===================================
  ! compute the RHS, accumulate into u(np1) and apply DSS
  !
  !           u(np1) = u(nm1) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! This subroutine was orgininally called to compute a leapfrog timestep
  ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
  ! accomodated.  For example, setting nm1=np1=n0 this routine will
  ! take a forward euler step, overwriting the input with the output.
  !
  !    qn0 = timelevel used to access Qdp() in order to compute virtual Temperature
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

  use physical_constants, only : cp, cpwater_vapor, Rgas, kappa, Rwater_vapor,p0, g
  use physics_mod, only : virtual_specific_heat, virtual_temperature
  use prim_si_mod, only : preq_vertadv_v, preq_vertadv_upwind, preq_omega_ps, preq_hydrostatic_v2

  use time_mod, only : tevolve

  implicit none
  integer, intent(in) :: np1,nm1,n0,qn0,nets,nete
  real*8, intent(in) :: dt2
  logical, intent(in)  :: compute_diagnostics,stiff

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind) :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux

  ! local
  real (kind=real_kind), pointer, dimension(:,:,:)   :: phi
  real (kind=real_kind), pointer, dimension(:,:,:)   :: dp3d
  real (kind=real_kind), pointer, dimension(:,:,:)   :: theta_dp_cp
   
  real (kind=real_kind) :: kappa_star(np,np,nlev)
  real (kind=real_kind) :: theta_cp(np,np,nlev)
  real (kind=real_kind) :: omega_p(np,np,nlev)
  real (kind=real_kind) :: vort(np,np,nlev)      ! vorticity
  real (kind=real_kind) :: divdp(np,np,nlev)     
  real (kind=real_kind) :: pnh(np,np,nlev)     ! nh (nonydro) pressure
  real (kind=real_kind) :: dpnh(np,np,nlev)    ! 
  real (kind=real_kind) :: exner(np,np,nlev)     ! exner nh pressure
  real (kind=real_kind) :: dpnh_dp(np,np,nlev)   ! dpnh / dp3d  
  real (kind=real_kind) :: eta_dot_dpdn(np,np,nlevp)  ! vertical velocity at interfaces
  real (kind=real_kind) :: KE(np,np,nlev)           ! Kinetic energy
  real (kind=real_kind) :: gradexner(np,np,2,nlev)  ! grad(p^kappa)   
  real (kind=real_kind) :: gradphi(np,np,2,nlev)     
  real (kind=real_kind) :: gradKE(np,np,2,nlev)  ! grad(0.5 u^T u )
  
  real (kind=real_kind) :: v_gradw(np,np,nlev)     
  real (kind=real_kind) :: v_gradtheta(np,np,nlev)     
  real (kind=real_kind) :: v_theta(np,np,2,nlev)
  real (kind=real_kind) :: div_v_theta(np,np,nlev)
  real (kind=real_kind) :: v_gradphi(np,np,nlev)
  real (kind=real_kind) :: v_gradKE(np,np,nlev)     
  real (kind=real_kind) :: vdp(np,np,2,nlev)

  real (kind=real_kind) :: vtens1(np,np,nlev)
  real (kind=real_kind) :: vtens2(np,np,nlev)
  real (kind=real_kind) :: v_vadv(np,np,2,nlev)   ! velocity vertical advection
  real (kind=real_kind) :: s_state(np,np,nlev,2)  ! scalars w,theta,phi
  real (kind=real_kind) :: s_vadv(np,np,nlev,3)   ! scalar vertical advection 
  real (kind=real_kind) :: s_theta_dp_cpadv(np,np,nlev)
  real (kind=real_kind) :: stens(np,np,nlev,3)    ! tendencies w,theta,phi

  real (kind=real_kind) ::  temp(np,np,nlev)
  real (kind=real_kind), dimension(np,np)      :: sdot_sum   ! temporary field
  real (kind=real_kind), dimension(np,np,2)    :: vtemp     ! generic gradient storage
  real (kind=real_kind) ::  v1,v2,w,d_eta_dot_dpdn_dn
  integer :: i,j,k,kptr,ie

  call t_startf('compute_and_apply_rhs_imex')
  do ie=nets,nete
     dp3d  => elem(ie)%state%dp3d(:,:,:,n0)
     theta_dp_cp  => elem(ie)%state%theta_dp_cp(:,:,:,n0)
     theta_cp(:,:,:) = theta_dp_cp(:,:,:)/dp3d(:,:,:)
     phi => elem(ie)%state%phi(:,:,:,n0)

     call get_kappa_star(kappa_star,elem(ie)%state%Qdp(:,:,:,1,qn0),dp3d)

     call get_pnh_and_exner(hvcoord,theta_dp_cp,dp3d,phi,elem(ie)%state%phis,&
             kappa_star,pnh,dpnh,exner) ! ,exner_i)
    
     if (theta_hydrostatic_mode) then
        ! Compute Hydrostatic equation
        do k=1,nlev
           temp(:,:,k) = kappa_star(:,:,k)*theta_dp_cp(:,:,k)*exner(:,:,k)/pnh(:,:,k)
        enddo
        !call preq_hydrostatic(phi,elem(ie)%state%phis,T_v,p,dp)
        call preq_hydrostatic_v2(phi,elem(ie)%state%phis,temp)
        dpnh_dp(:,:,:) = 1
     else
        ! d(p-nh) / d(p-hyrdostatic)
        dpnh_dp(:,:,:) = dpnh(:,:,:)/dp3d(:,:,:)
     endif

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        vdp(:,:,1,k) = elem(ie)%state%v(:,:,1,k,n0)*dp3d(:,:,k)
        vdp(:,:,2,k) = elem(ie)%state%v(:,:,2,k,n0)*dp3d(:,:,k)

        ! ================================
        ! Accumulate mean Vel_rho flux in vn0
        ! ================================
        elem(ie)%derived%vn0(:,:,:,k)=elem(ie)%derived%vn0(:,:,:,k)+eta_ave_w*vdp(:,:,:,k)

        divdp(:,:,k)=divergence_sphere(vdp(:,:,:,k),deriv,elem(ie))
        vort(:,:,k)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,n0),deriv,elem(ie))

     enddo


     ! Compute omega_p according to CCM-3
     !call preq_omega_ps(omega_p,hvcoord,p,vgrad_p,divdp)
     ! how will we compute this?  omega_p = 1/p Dp/Dt
     omega_p = 0


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
        s_vadv=0
        s_theta_dp_cpadv=0
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
        ! TODO: remove theta from s_state and s_vadv
        s_state(:,:,:,1)=elem(ie)%state%w(:,:,:,n0)
        s_state(:,:,:,2)=elem(ie)%state%phi(:,:,:,n0)
        call preq_vertadv_v(elem(ie)%state%v(:,:,:,:,n0),s_state,2,eta_dot_dpdn,dp3d,v_vadv,s_vadv)
        !call preq_vertadv_v(elem(ie)%state%v(:,:,:,:,n0),s_state,3,eta_dot_dpdn,dp3d,v_vadv,s_vadv)
        !call preq_vertadv_upwind(elem(ie)%state%v(:,:,:,:,n0),s_state,3,eta_dot_dpdn,dp3d,v_vadv,s_vadv)

        !    this loop constructs d( eta-dot * theta_dp_cp)/deta
        !   d( eta_dot_dpdn * theta*cp)
        !  so we need to compute theta_cp form theta_dp_cp and average to interfaces
        k=1
        s_theta_dp_cpadv(:,:,k)= &
             eta_dot_dpdn(:,:,k+1)* (theta_cp(:,:,k+1)+theta_cp(:,:,k))/2  

        k=nlev
        s_theta_dp_cpadv(:,:,k)= - &
             eta_dot_dpdn(:,:,k)  * (theta_cp(:,:,k)+theta_cp(:,:,k-1))/2  

        do k=2,nlev-1
           s_theta_dp_cpadv(:,:,k)= &
              eta_dot_dpdn(:,:,k+1)* (theta_cp(:,:,k+1)+theta_cp(:,:,k))/2  - &
              eta_dot_dpdn(:,:,k)  * (theta_cp(:,:,k)+theta_cp(:,:,k-1))/2  
        end do
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
!$omp parallel do private(k,i,j,v1,v2,KE,vtemp)
#endif

     vertloop: do k=1,nlev

        ! ================================================
        ! w,theta,phi tendencies:
        ! ================================================
       if (stiff) then ! return stiff component of RHS
         stens(:,:,k,1) = -g*(1-dpnh_dp(:,:,k))
         stens(:,:,k,2) = 0.d0
         stens(:,:,k,3) = g*elem(ie)%state%w(:,:,k,n0)
         vtens1(:,:,:) = 0.d0 
         vtens2(:,:,:) = 0.d0
       else ! return nonstiff component of RHS
        vtemp(:,:,:)   = gradient_sphere(elem(ie)%state%w(:,:,k,n0),deriv,elem(ie)%Dinv)
        v_gradw(:,:,k) = elem(ie)%state%v(:,:,1,k,n0)*vtemp(:,:,1) &
             +elem(ie)%state%v(:,:,2,k,n0)*vtemp(:,:,2) 
        stens(:,:,k,1) = -s_vadv(:,:,k,1) - v_gradw(:,:,k) !-  g*(1-dpnh_dp(:,:,k) )
        v_theta(:,:,1,k) = elem(ie)%state%v(:,:,1,k,n0)*               &
          elem(ie)%state%theta_dp_cp(:,:,k,n0)
        v_theta(:,:,2,k) =                                             &
          elem(ie)%state%v(:,:,2,k,n0)                                 &
          *elem(ie)%state%theta_dp_cp(:,:,k,n0)
        div_v_theta(:,:,k)=divergence_sphere(v_theta(:,:,:,k),deriv,elem(ie))
        stens(:,:,k,2)=-s_theta_dp_cpadv(:,:,k)-div_v_theta(:,:,k)

        gradphi(:,:,:,k) = gradient_sphere(phi(:,:,k),deriv,elem(ie)%Dinv)

        v_gradphi(:,:,k) = elem(ie)%state%v(:,:,1,k,n0)*gradphi(:,:,1,k) &
             +elem(ie)%state%v(:,:,2,k,n0)*gradphi(:,:,2,k) 
        ! use of s_vadv(:,:,k,2) here is correct since this corresponds to etadot d(phi)/deta
        stens(:,:,k,3) =  -s_vadv(:,:,k,2) - v_gradphi(:,:,k) !+ g*elem(ie)%state%w(:,:,k,n0)

        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              KE(i,j,k)=0.5D0*( v1*v1 + v2*v2 )
           end do
        end do
        gradKE(:,:,:,k) = gradient_sphere(KE(:,:,k),deriv,elem(ie)%Dinv)
  
     gradexner(:,:,:,k) = gradient_sphere(exner(:,:,k),deriv,elem(ie)%Dinv)        

        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              v_gradKE(i,j,k)=v1*gradKE(i,j,1,k)+v2*gradKE(i,j,2,k)

              vtens1(i,j,k) = -v_vadv(i,j,1,k) &
                   + v2*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - gradKE(i,j,1,k) -gradphi(i,j,1,k)*dpnh_dp(i,j,k) &
                   -theta_cp(i,j,k)*gradexner(i,j,1,k)

              vtens2(i,j,k) = -v_vadv(i,j,2,k) &
                   - v1*(elem(ie)%fcor(i,j) + vort(i,j,k)) &
                   - gradKE(i,j,2,k) -gradphi(i,j,2,k)*dpnh_dp(i,j,k) &
                   -theta_cp(i,j,k)*gradexner(i,j,2,k)
           end do
        end do
      end if
     end do vertloop


     
#ifdef ENERGY_DIAGNOSTICS
     ! =========================================================
     !
     ! diagnostics
     ! recomputes some gradients that were not saved above
     ! uses:  sdot_sum(), eta_dot_dpdn(), grad_ps()
     ! grad_phi(), dp3d(), p(), theta_vadv(), v_vadv(), divdp()
     ! =========================================================
 
     ! =========================================================
     ! (AAM) - This section has accumulations over vertical levels.
     !   Be careful if implementing OpenMP
     ! =========================================================

     if (compute_diagnostics) then
        elem(ie)%accum%KEhoriz1=0
        elem(ie)%accum%PEhoriz1=0
        elem(ie)%accum%PEhoriz2=0
        elem(ie)%accum%KE1=0
        elem(ie)%accum%KE2=0
        elem(ie)%accum%KEvert1=0
        elem(ie)%accum%KEvert2=0
        elem(ie)%accum%IEvert1=0
        elem(ie)%accum%IEvert2=0
        elem(ie)%accum%PEvert1=0
        elem(ie)%accum%PEvert2=0
        elem(ie)%accum%T01=0
        elem(ie)%accum%T2=0
        elem(ie)%accum%S1=0
        elem(ie)%accum%S2=0
        elem(ie)%accum%P1=0
        elem(ie)%accum%P2=0
        ! See element_state.F90 for an account of what these variables are defined as
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,v1,v2,vtemp,KE,d_eta_dot_dpdn_dn)
#endif
        do k =1,nlev
          do j=1,np
            do i=1,np                
                  d_eta_dot_dpdn_dn=(eta_dot_dpdn(i,j,k+1)-                    &
                  eta_dot_dpdn(i,j,k))
               !  Form KEhoriz1
                  elem(ie)%accum%KEhoriz1(i,j)=elem(ie)%accum%KEhoriz1(i,j)    &
                  -v_gradKE(i,j,k)*dp3d(i,j,k) - KE(i,j,k)*divdp(i,j,k)
                  elem(ie)%accum%KE1(i,j)=elem(ie)%accum%KE1(i,j)              &
                  -v_gradKE(i,j,k)*dp3d(i,j,k) 
                  elem(ie)%accum%KE2(i,j)=elem(ie)%accum%KE2(i,j)              &
                  -KE(i,j,k)*divdp(i,j,k)
               !  Form KEhoriz2
                  elem(ie)%accum%KEhoriz2(i,j)=elem(ie)%accum%KEhoriz2(i,j)-   &
                  dp3d(i,j,k) * elem(ie)%state%w(i,j,k,n0) * v_gradw(i,j,k)    &
                  -0.5*(elem(ie)%state%w(i,j,k,n0))**2 * divdp(i,j,k) 
               !  Form KEvert1
                  elem(ie)%accum%KEvert1(i,j)=elem(ie)%accum%KEvert1(i,j)-     &
                  (elem(ie)%state%v(i,j,1,k,n0) * v_vadv(i,j,1,k) +            &
                  elem(ie)%state%v(i,j,2,k,n0) *v_vadv(i,j,2,k))*dp3d(i,j,k)-  &      
                  0.5*((elem(ie)%state%v(i,j,1,k,n0))**2 +                     &
                       (elem(ie)%state%v(i,j,2,k,n0))**2)*d_eta_dot_dpdn_dn
               !  Form KEvert2
                  elem(ie)%accum%KEvert2(i,j)=elem(ie)%accum%KEvert2(i,j)      &
                  -s_vadv(i,j,k,1)*(elem(ie)%state%w(i,j,k,n0))*dp3d(i,j,k)    &   
                  -d_eta_dot_dpdn_dn*(elem(ie)%state%w(i,j,k,n0))**2
               !  Form IEvert1
                  elem(ie)%accum%IEvert1(i,j)=elem(ie)%accum%IEvert1(i,j)      &
                  -exner(i,j,k)*s_theta_dp_cpadv(i,j,k)                        
               !  Form IEvert2
                  elem(ie)%accum%IEvert2(i,j)=elem(ie)%accum%IEvert2(i,j)      &
                  +dpnh(i,j,k)*s_vadv(i,j,k,2)
               !  Form PEhoriz1
                  elem(ie)%accum%PEhoriz1(i,j)=(elem(ie)%accum%PEhoriz1(i,j))  &
                  -phi(i,j,k)*divdp(i,j,k) 
               !  Form PEhoriz2                                                &
                  elem(ie)%accum%PEhoriz2(i,j)=elem(ie)%accum%PEhoriz2(i,j)    &
                  -dp3d(i,j,k)*v_gradphi(i,j,k)      
               !  Form PEvert1
                  elem(ie)%accum%PEvert1(i,j) = (elem(ie)%accum%PEvert1(i,j))   &
                  -phi(i,j,k)*d_eta_dot_dpdn_dn                                 
               !  Form PEvert2
                  elem(ie)%accum%PEvert2(i,j) = elem(ie)%accum%PEvert2(i,j)     &
                  -dp3d(i,j,k)*s_theta_dp_cpadv(i,j,k)
               !  Form T01
                  elem(ie)%accum%T01(i,j)=elem(ie)%accum%T01(i,j)               &
                  -(elem(ie)%state%theta_dp_cp(i,j,k,n0))                       &
                  *(gradexner(i,j,1,k)*elem(ie)%state%v(i,j,1,k,n0) +           &
                  gradexner(i,j,2,k)*elem(ie)%state%v(i,j,2,k,n0))              
               !  Form S1 
                  elem(ie)%accum%S1(i,j)=elem(ie)%accum%S1(i,j)                 &
                  -exner(i,j,k)*div_v_theta(i,j,k)
               !  Form T2 
                  elem(ie)%accum%T2(i,j)=elem(ie)%accum%T2(i,j)+                & 
                  (g*(elem(ie)%state%w(i,j,k,n0))-                              &
                  v_gradphi(i,j,k))*dpnh(i,j,k)                                 
               !  Form S2
                  elem(ie)%accum%S2(i,j)=elem(ie)%accum%S2(i,j)                 &
                  -g*(elem(ie)%state%w(i,j,k,n0))+v_gradphi(i,j,k)              &
                  *dpnh(i,j,k)
               !  Form P1
                  elem(ie)%accum%P1(i,j)=elem(ie)%accum%P1(i,j)                 &
                  -g*(elem(ie)%state%w(i,j,k,n0)) * dp3d(i,j,k)
               !  Form P2
                  elem(ie)%accum%P2(i,j)=elem(ie)%accum%P2(i,j)                 &
                  + g * (elem(ie)%state%w(i,j,k,n0)) * dp3d(i,j,k)
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
       if (stiff) then
        elem(ie)%state%w(:,:,k,np1) =  elem(ie)%state%w(:,:,k,nm1) + dt2*stens(:,:,k,1)
        elem(ie)%state%phi(:,:,k,np1) = elem(ie)%state%phi(:,:,k,nm1)   + dt2*stens(:,:,k,3)
       else
        elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%spheremp(:,:)*( elem(ie)%state%v(:,:,1,k,nm1) + dt2*vtens1(:,:,k) )
        elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%spheremp(:,:)*( elem(ie)%state%v(:,:,2,k,nm1) + dt2*vtens2(:,:,k) )

        elem(ie)%state%w(:,:,k,np1)     = elem(ie)%spheremp(:,:)*(elem(ie)%state%w(:,:,k,nm1)     + dt2*stens(:,:,k,1))
        elem(ie)%state%theta_dp_cp(:,:,k,np1) = elem(ie)%spheremp(:,:)*(elem(ie)%state%theta_dp_cp(:,:,k,nm1) + dt2*stens(:,:,k,2))
        elem(ie)%state%phi(:,:,k,np1)   = elem(ie)%spheremp(:,:)*(elem(ie)%state%phi(:,:,k,nm1)   + dt2*stens(:,:,k,3))
       
        elem(ie)%state%dp3d(:,:,k,np1) = &
             elem(ie)%spheremp(:,:) * (elem(ie)%state%dp3d(:,:,k,nm1) - &
             dt2 * (divdp(:,:,k) + eta_dot_dpdn(:,:,k+1)-eta_dot_dpdn(:,:,k)))
       end if
     enddo
    if (stiff) then
 
    else
     kptr=0
     call edgeVpack(edge6, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr, ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%theta_dp_cp(:,:,:,np1),nlev,kptr,ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%w(:,:,:,np1),nlev,kptr,ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%phi(:,:,:,np1),nlev,kptr,ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,ie)
    end if
  end do

  if (stiff) then

  else

    call t_startf('caar_bexchV')
  call bndry_exchangeV(hybrid,edge6)
  call t_stopf('caar_bexchV')

  do ie=nets,nete
     kptr=0
     call edgeVunpack(edge6, elem(ie)%state%dp3d(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%theta_dp_cp(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%w(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%phi(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, ie)

     
     ! ====================================================
     ! Scale tendencies by inverse mass matrix
     ! ====================================================
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        elem(ie)%state%dp3d(:,:,k,np1) =elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
        elem(ie)%state%theta_dp_cp(:,:,k,np1)=elem(ie)%rspheremp(:,:)*elem(ie)%state%theta_dp_cp(:,:,k,np1)
        elem(ie)%state%w(:,:,k,np1)    =elem(ie)%rspheremp(:,:)*elem(ie)%state%w(:,:,k,np1)
        elem(ie)%state%phi(:,:,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%phi(:,:,k,np1)
        elem(ie)%state%v(:,:,1,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,1,k,np1)
        elem(ie)%state%v(:,:,2,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,2,k,np1)
     end do
  end do
  end if
  call t_stopf('compute_and_apply_rhs_imex')

  end subroutine compute_and_apply_rhs_imex


!==================================== non-stiff =================================================================
!================================================================================================================
!================================================================================================================
!================================================================================================================
!================================================================================================================
!================================================================================================================
!================================================================================================================
!================================================================================================================
!================================================================================================================
!================================================================================================================
 

 subroutine compute_and_apply_rhs_imex_nonstiff(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w,scale1,scale2,scale3)
  ! ===================================
  ! compute the RHS, accumulate into u(np1) and apply DSS
  !
  !           u(np1) = u(nm1) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! This subroutine was orgininally called to compute a leapfrog timestep
  ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
  ! accomodated.  For example, setting nm1=np1=n0 this routine will
  ! take a forward euler step, overwriting the input with the output.
  !
  !    qn0 = timelevel used to access Qdp() in order to compute virtual Temperature
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

  use physical_constants, only : cp, cpwater_vapor, Rgas, kappa, Rwater_vapor,p0, g
  use physics_mod, only : virtual_specific_heat, virtual_temperature
  use prim_si_mod, only : preq_vertadv_v, preq_vertadv_upwind, preq_omega_ps, preq_hydrostatic_v2

  use time_mod, only : tevolve

  implicit none
  integer, intent(in) :: np1,nm1,n0,qn0,nets,nete
  real*8, intent(in) :: dt2
  logical, intent(in)  :: compute_diagnostics

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind) :: eta_ave_w,scale1,scale2,scale3  ! weighting for eta_dot_dpdn mean flux, scale of unm1

  ! local
  real (kind=real_kind), pointer, dimension(:,:,:)   :: phi
  real (kind=real_kind), pointer, dimension(:,:,:)   :: dp3d
  real (kind=real_kind), pointer, dimension(:,:,:)   :: theta_dp_cp
   
  real (kind=real_kind) :: kappa_star(np,np,nlev)
  real (kind=real_kind) :: theta_cp(np,np,nlev)
  real (kind=real_kind) :: omega_p(np,np,nlev)
  real (kind=real_kind) :: vort(np,np,nlev)      ! vorticity
  real (kind=real_kind) :: divdp(np,np,nlev)     
  real (kind=real_kind) :: pnh(np,np,nlev)     ! nh (nonydro) pressure
  real (kind=real_kind) :: dpnh(np,np,nlev)    ! 
  real (kind=real_kind) :: exner(np,np,nlev)     ! exner nh pressure
  real (kind=real_kind) :: dpnh_dp(np,np,nlev)   ! dpnh / dp3d  
  real (kind=real_kind) :: eta_dot_dpdn(np,np,nlevp)  ! vertical velocity at interfaces
  real (kind=real_kind) :: KE(np,np,nlev)           ! Kinetic energy
  real (kind=real_kind) :: gradexner(np,np,2,nlev)  ! grad(p^kappa)   
  real (kind=real_kind) :: gradphi(np,np,2,nlev)     
  real (kind=real_kind) :: gradKE(np,np,2,nlev)  ! grad(0.5 u^T u )
  
  real (kind=real_kind) :: v_gradw(np,np,nlev)     
  real (kind=real_kind) :: v_gradtheta(np,np,nlev)     
  real (kind=real_kind) :: v_theta(np,np,2,nlev)
  real (kind=real_kind) :: div_v_theta(np,np,nlev)
  real (kind=real_kind) :: v_gradphi(np,np,nlev)
  real (kind=real_kind) :: v_gradKE(np,np,nlev)     
  real (kind=real_kind) :: vdp(np,np,2,nlev)

  real (kind=real_kind) :: vtens1(np,np,nlev)
  real (kind=real_kind) :: vtens2(np,np,nlev)
  real (kind=real_kind) :: v_vadv(np,np,2,nlev)   ! velocity vertical advection
  real (kind=real_kind) :: s_state(np,np,nlev,2)  ! scalars w,theta,phi
  real (kind=real_kind) :: s_vadv(np,np,nlev,3)   ! scalar vertical advection 
  real (kind=real_kind) :: s_theta_dp_cpadv(np,np,nlev)
  real (kind=real_kind) :: stens(np,np,nlev,3)    ! tendencies w,theta,phi

  real (kind=real_kind) ::  temp(np,np,nlev)
  real (kind=real_kind), dimension(np,np)      :: sdot_sum   ! temporary field
  real (kind=real_kind), dimension(np,np,2)    :: vtemp     ! generic gradient storage
  real (kind=real_kind) ::  v1,v2,w,d_eta_dot_dpdn_dn
  integer :: i,j,k,kptr,ie

  call t_startf('compute_and_apply_rhs_imex_nonstiff')
  do ie=nets,nete
     dp3d  => elem(ie)%state%dp3d(:,:,:,n0)
     theta_dp_cp  => elem(ie)%state%theta_dp_cp(:,:,:,n0)
     theta_cp(:,:,:) = theta_dp_cp(:,:,:)/dp3d(:,:,:)
     phi => elem(ie)%state%phi(:,:,:,n0)

     call get_kappa_star(kappa_star,elem(ie)%state%Qdp(:,:,:,1,qn0),dp3d)

     call get_pnh_and_exner(hvcoord,theta_dp_cp,dp3d,phi,elem(ie)%state%phis,&
             kappa_star,pnh,dpnh,exner) ! ,exner_i)
    
     if (theta_hydrostatic_mode) then
        ! Compute Hydrostatic equation
        do k=1,nlev
           temp(:,:,k) = kappa_star(:,:,k)*theta_dp_cp(:,:,k)*exner(:,:,k)/pnh(:,:,k)
        enddo
        !call preq_hydrostatic(phi,elem(ie)%state%phis,T_v,p,dp)
        call preq_hydrostatic_v2(phi,elem(ie)%state%phis,temp)
        dpnh_dp(:,:,:) = 1
     else
        ! d(p-nh) / d(p-hyrdostatic)
        dpnh_dp(:,:,:) = dpnh(:,:,:)/dp3d(:,:,:)
     endif

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        vdp(:,:,1,k) = elem(ie)%state%v(:,:,1,k,n0)*dp3d(:,:,k)
        vdp(:,:,2,k) = elem(ie)%state%v(:,:,2,k,n0)*dp3d(:,:,k)

        ! ================================
        ! Accumulate mean Vel_rho flux in vn0
        ! ================================
        elem(ie)%derived%vn0(:,:,:,k)=elem(ie)%derived%vn0(:,:,:,k)+eta_ave_w*vdp(:,:,:,k)

        divdp(:,:,k)=divergence_sphere(vdp(:,:,:,k),deriv,elem(ie))
        vort(:,:,k)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,n0),deriv,elem(ie))

     enddo


     ! Compute omega_p according to CCM-3
     !call preq_omega_ps(omega_p,hvcoord,p,vgrad_p,divdp)
     ! how will we compute this?  omega_p = 1/p Dp/Dt
     omega_p = 0


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
        s_vadv=0
        s_theta_dp_cpadv=0
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
        ! TODO: remove theta from s_state and s_vadv
        s_state(:,:,:,1)=elem(ie)%state%w(:,:,:,n0)
        s_state(:,:,:,2)=elem(ie)%state%phi(:,:,:,n0)
        call preq_vertadv_v(elem(ie)%state%v(:,:,:,:,n0),s_state,2,eta_dot_dpdn,dp3d,v_vadv,s_vadv)
        !call preq_vertadv_v(elem(ie)%state%v(:,:,:,:,n0),s_state,3,eta_dot_dpdn,dp3d,v_vadv,s_vadv)
        !call preq_vertadv_upwind(elem(ie)%state%v(:,:,:,:,n0),s_state,3,eta_dot_dpdn,dp3d,v_vadv,s_vadv)

        !    this loop constructs d( eta-dot * theta_dp_cp)/deta
        !   d( eta_dot_dpdn * theta*cp)
        !  so we need to compute theta_cp form theta_dp_cp and average to interfaces
        k=1
        s_theta_dp_cpadv(:,:,k)= &
             eta_dot_dpdn(:,:,k+1)* (theta_cp(:,:,k+1)+theta_cp(:,:,k))/2  

        k=nlev
        s_theta_dp_cpadv(:,:,k)= - &
             eta_dot_dpdn(:,:,k)  * (theta_cp(:,:,k)+theta_cp(:,:,k-1))/2  

        do k=2,nlev-1
           s_theta_dp_cpadv(:,:,k)= &
              eta_dot_dpdn(:,:,k+1)* (theta_cp(:,:,k+1)+theta_cp(:,:,k))/2  - &
              eta_dot_dpdn(:,:,k)  * (theta_cp(:,:,k)+theta_cp(:,:,k-1))/2  
        end do
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
!$omp parallel do private(k,i,j,v1,v2,KE,vtemp)
#endif

     vertloop: do k=1,nlev
                   
        ! ================================================
        ! w,theta,phi tendencies:
        ! ================================================
        vtemp(:,:,:)   = gradient_sphere(elem(ie)%state%w(:,:,k,n0),deriv,elem(ie)%Dinv)
        v_gradw(:,:,k) = elem(ie)%state%v(:,:,1,k,n0)*vtemp(:,:,1) &
             +elem(ie)%state%v(:,:,2,k,n0)*vtemp(:,:,2) 
        stens(:,:,k,1) = (-s_vadv(:,:,k,1) - v_gradw(:,:,k))*scale1 - scale2*g*(1-dpnh_dp(:,:,k) )
        v_theta(:,:,1,k) = elem(ie)%state%v(:,:,1,k,n0)*               &
          elem(ie)%state%theta_dp_cp(:,:,k,n0)
        v_theta(:,:,2,k) =                                             &
          elem(ie)%state%v(:,:,2,k,n0)                                 &
          *elem(ie)%state%theta_dp_cp(:,:,k,n0)
        div_v_theta(:,:,k)=divergence_sphere(v_theta(:,:,:,k),deriv,elem(ie))
        stens(:,:,k,2)=(-s_theta_dp_cpadv(:,:,k)-div_v_theta(:,:,k))*scale1

        gradphi(:,:,:,k) = gradient_sphere(phi(:,:,k),deriv,elem(ie)%Dinv)

        v_gradphi(:,:,k) = elem(ie)%state%v(:,:,1,k,n0)*gradphi(:,:,1,k) &
             +elem(ie)%state%v(:,:,2,k,n0)*gradphi(:,:,2,k) 
        ! use of s_vadv(:,:,k,2) here is correct since this corresponds to etadot d(phi)/deta
        stens(:,:,k,3) =  (-s_vadv(:,:,k,2) - v_gradphi(:,:,k))*scale1 + scale2*g*elem(ie)%state%w(:,:,k,n0)

        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              KE(i,j,k)=0.5D0*( v1*v1 + v2*v2 )
           end do
        end do
        gradKE(:,:,:,k) = gradient_sphere(KE(:,:,k),deriv,elem(ie)%Dinv)
  
     gradexner(:,:,:,k) = gradient_sphere(exner(:,:,k),deriv,elem(ie)%Dinv)        

        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              v_gradKE(i,j,k)=v1*gradKE(i,j,1,k)+v2*gradKE(i,j,2,k)

              vtens1(i,j,k) = (-v_vadv(i,j,1,k) &
                   + v2*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - gradKE(i,j,1,k) - gradphi(i,j,1,k)*dpnh_dp(i,j,k) &
                   -theta_cp(i,j,k)*gradexner(i,j,1,k))*scale1

              vtens2(i,j,k) = (-v_vadv(i,j,2,k) &
                   - v1*(elem(ie)%fcor(i,j) + vort(i,j,k)) &
                   - gradKE(i,j,2,k) - gradphi(i,j,2,k)*dpnh_dp(i,j,k) &
                   -theta_cp(i,j,k)*gradexner(i,j,2,k))*scale1
           end do
        end do     
     end do vertloop

     
#ifdef ENERGY_DIAGNOSTICS
     ! =========================================================
     !
     ! diagnostics
     ! recomputes some gradients that were not saved above
     ! uses:  sdot_sum(), eta_dot_dpdn(), grad_ps()
     ! grad_phi(), dp3d(), p(), theta_vadv(), v_vadv(), divdp()
     ! =========================================================
 
     ! =========================================================
     ! (AAM) - This section has accumulations over vertical levels.
     !   Be careful if implementing OpenMP
     ! =========================================================

     if (compute_diagnostics) then
        elem(ie)%accum%KEhoriz1=0
        elem(ie)%accum%PEhoriz1=0
        elem(ie)%accum%PEhoriz2=0
        elem(ie)%accum%KE1=0
        elem(ie)%accum%KE2=0
        elem(ie)%accum%KEvert1=0
        elem(ie)%accum%KEvert2=0
        elem(ie)%accum%IEvert1=0
        elem(ie)%accum%IEvert2=0
        elem(ie)%accum%PEvert1=0
        elem(ie)%accum%PEvert2=0
        elem(ie)%accum%T01=0
        elem(ie)%accum%T2=0
        elem(ie)%accum%S1=0
        elem(ie)%accum%S2=0
        elem(ie)%accum%P1=0
        elem(ie)%accum%P2=0
        ! See element_state.F90 for an account of what these variables are defined as
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,v1,v2,vtemp,KE,d_eta_dot_dpdn_dn)
#endif
        do k =1,nlev
          do j=1,np
            do i=1,np                
                  d_eta_dot_dpdn_dn=(eta_dot_dpdn(i,j,k+1)-                    &
                  eta_dot_dpdn(i,j,k))
               !  Form KEhoriz1
                  elem(ie)%accum%KEhoriz1(i,j)=elem(ie)%accum%KEhoriz1(i,j)    &
                  -v_gradKE(i,j,k)*dp3d(i,j,k) - KE(i,j,k)*divdp(i,j,k)
                  elem(ie)%accum%KE1(i,j)=elem(ie)%accum%KE1(i,j)              &
                  -v_gradKE(i,j,k)*dp3d(i,j,k) 
                  elem(ie)%accum%KE2(i,j)=elem(ie)%accum%KE2(i,j)              &
                  -KE(i,j,k)*divdp(i,j,k)
               !  Form KEhoriz2
                  elem(ie)%accum%KEhoriz2(i,j)=elem(ie)%accum%KEhoriz2(i,j)-   &
                  dp3d(i,j,k) * elem(ie)%state%w(i,j,k,n0) * v_gradw(i,j,k)    &
                  -0.5*(elem(ie)%state%w(i,j,k,n0))**2 * divdp(i,j,k) 
               !  Form KEvert1
                  elem(ie)%accum%KEvert1(i,j)=elem(ie)%accum%KEvert1(i,j)-     &
                  (elem(ie)%state%v(i,j,1,k,n0) * v_vadv(i,j,1,k) +            &
                  elem(ie)%state%v(i,j,2,k,n0) *v_vadv(i,j,2,k))*dp3d(i,j,k)-  &      
                  0.5*((elem(ie)%state%v(i,j,1,k,n0))**2 +                     &
                       (elem(ie)%state%v(i,j,2,k,n0))**2)*d_eta_dot_dpdn_dn
               !  Form KEvert2
                  elem(ie)%accum%KEvert2(i,j)=elem(ie)%accum%KEvert2(i,j)      &
                  -s_vadv(i,j,k,1)*(elem(ie)%state%w(i,j,k,n0))*dp3d(i,j,k)    &   
                  -d_eta_dot_dpdn_dn*(elem(ie)%state%w(i,j,k,n0))**2
               !  Form IEvert1
                  elem(ie)%accum%IEvert1(i,j)=elem(ie)%accum%IEvert1(i,j)      &
                  -exner(i,j,k)*s_theta_dp_cpadv(i,j,k)                        
               !  Form IEvert2
                  elem(ie)%accum%IEvert2(i,j)=elem(ie)%accum%IEvert2(i,j)      &
                  +dpnh(i,j,k)*s_vadv(i,j,k,2)
               !  Form PEhoriz1
                  elem(ie)%accum%PEhoriz1(i,j)=(elem(ie)%accum%PEhoriz1(i,j))  &
                  -phi(i,j,k)*divdp(i,j,k) 
               !  Form PEhoriz2                                                &
                  elem(ie)%accum%PEhoriz2(i,j)=elem(ie)%accum%PEhoriz2(i,j)    &
                  -dp3d(i,j,k)*v_gradphi(i,j,k)      
               !  Form PEvert1
                  elem(ie)%accum%PEvert1(i,j) = (elem(ie)%accum%PEvert1(i,j))   &
                  -phi(i,j,k)*d_eta_dot_dpdn_dn                                 
               !  Form PEvert2
                  elem(ie)%accum%PEvert2(i,j) = elem(ie)%accum%PEvert2(i,j)     &
                  -dp3d(i,j,k)*s_theta_dp_cpadv(i,j,k)
               !  Form T01
                  elem(ie)%accum%T01(i,j)=elem(ie)%accum%T01(i,j)               &
                  -(elem(ie)%state%theta_dp_cp(i,j,k,n0))                       &
                  *(gradexner(i,j,1,k)*elem(ie)%state%v(i,j,1,k,n0) +           &
                  gradexner(i,j,2,k)*elem(ie)%state%v(i,j,2,k,n0))              
               !  Form S1 
                  elem(ie)%accum%S1(i,j)=elem(ie)%accum%S1(i,j)                 &
                  -exner(i,j,k)*div_v_theta(i,j,k)
               !  Form T2 
                  elem(ie)%accum%T2(i,j)=elem(ie)%accum%T2(i,j)+                & 
                  (g*(elem(ie)%state%w(i,j,k,n0))-                              &
                  v_gradphi(i,j,k))*dpnh(i,j,k)                                 
               !  Form S2
                  elem(ie)%accum%S2(i,j)=elem(ie)%accum%S2(i,j)                 &
                  -g*(elem(ie)%state%w(i,j,k,n0))+v_gradphi(i,j,k)              &
                  *dpnh(i,j,k)
               !  Form P1
                  elem(ie)%accum%P1(i,j)=elem(ie)%accum%P1(i,j)                 &
                  -g*(elem(ie)%state%w(i,j,k,n0)) * dp3d(i,j,k)
               !  Form P2
                  elem(ie)%accum%P2(i,j)=elem(ie)%accum%P2(i,j)                 &
                  + g * (elem(ie)%state%w(i,j,k,n0)) * dp3d(i,j,k)
              enddo
            enddo
          enddo 
     endif
#endif
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%v(:,:,1,k,nm1) &
          + dt2*vtens1(:,:,k) )
        elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%v(:,:,2,k,nm1) &
          +  dt2*vtens2(:,:,k) )
        elem(ie)%state%w(:,:,k,np1)    = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%w(:,:,k,nm1)   &
          + dt2*stens(:,:,k,1))
        elem(ie)%state%theta_dp_cp(:,:,k,np1) = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%theta_dp_cp(:,:,k,nm1) &
          + dt2*stens(:,:,k,2))
        elem(ie)%state%phi(:,:,k,np1)   = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%phi(:,:,k,nm1) & 
          + dt2*stens(:,:,k,3))

        elem(ie)%state%dp3d(:,:,k,np1) = &
             elem(ie)%spheremp(:,:) * (scale3 * elem(ie)%state%dp3d(:,:,k,nm1) - &
             scale1*dt2 * (divdp(:,:,k) + eta_dot_dpdn(:,:,k+1)-eta_dot_dpdn(:,:,k)))
        
     enddo

     kptr=0
     call edgeVpack(edge6, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr, ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%theta_dp_cp(:,:,:,np1),nlev,kptr,ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%w(:,:,:,np1),nlev,kptr,ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%phi(:,:,:,np1),nlev,kptr,ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,ie)
   end do ! end do for the ie=nets,nete loop

  call t_startf('caar_bexchV')
  call bndry_exchangeV(hybrid,edge6)
  call t_stopf('caar_bexchV')

  do ie=nets,nete
     kptr=0
     call edgeVunpack(edge6, elem(ie)%state%dp3d(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%theta_dp_cp(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%w(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%phi(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, ie)

      
     ! ====================================================
     ! Scale tendencies by inverse mass matrix
     ! ====================================================
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        elem(ie)%state%dp3d(:,:,k,np1) =elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
        elem(ie)%state%theta_dp_cp(:,:,k,np1)=elem(ie)%rspheremp(:,:)*elem(ie)%state%theta_dp_cp(:,:,k,np1)
        elem(ie)%state%w(:,:,k,np1)    =elem(ie)%rspheremp(:,:)*elem(ie)%state%w(:,:,k,np1)
        elem(ie)%state%phi(:,:,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%phi(:,:,k,np1)
        elem(ie)%state%v(:,:,1,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,1,k,np1)
        elem(ie)%state%v(:,:,2,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,2,k,np1)
     end do
  end do


  call t_stopf('compute_and_apply_rhs_imex_nonstiff')

  end subroutine compute_and_apply_rhs_imex_nonstiff
 




! ====================================== stiff =================================================================
!================================================================================================================
!================================================================================================================
!================================================================================================================
!================================================================================================================
!================================================================================================================
!================================================================================================================
!================================================================================================================
!================================================================================================================
!================================================================================================================


  subroutine compute_and_apply_rhs_imex_stiff(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w,maxiter,itertol,statesave) !Arki,stagevalues)
  ! ===================================
  ! compute the RHS, accumulate into u(np1) and apply DSS
  !
  !           u(np1) = u(nm1) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! This subroutine was orgininally called to compute a leapfrog timestep
  ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
  ! accomodated.  For example, setting nm1=np1=n0 this routine will
  ! take a forward euler step, overwriting the input with the output.
  !
  !    qn0 = timelevel used to access Qdp() in order to compute virtual Temperature
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

  use physical_constants, only : cp, cpwater_vapor, Rgas, kappa, Rwater_vapor,p0, g
  use physics_mod, only : virtual_specific_heat, virtual_temperature
  use prim_si_mod, only : preq_vertadv_v, preq_vertadv_upwind, preq_omega_ps, preq_hydrostatic_v2

  use time_mod, only : tevolve

  implicit none
  integer, intent(in) :: np1,nm1,n0,qn0,nets,nete
  real*8, intent(in) :: dt2
  logical, intent(in)  :: compute_diagnostics
  integer :: maxiter
  real*8 :: itertol

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind) :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux
  real (kind=real_kind) :: statesave(nets:nete,np,np,nlev,6)

  ! local
  real (kind=real_kind), pointer, dimension(:,:,:)   :: phi
  real (kind=real_kind), pointer, dimension(:,:,:)   :: dp3d
  real (kind=real_kind), pointer, dimension(:,:,:)   :: theta_dp_cp
   
  real (kind=real_kind) :: kappa_star(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev)     ! nh (nonydro) pressure
  real (kind=real_kind) :: dpnh(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)     ! exner nh pressure
  real (kind=real_kind) :: dpnh_dp(np,np,nlev),dpnh_dp2(np,np,nlev)    !    ! dpnh / dp3d  
  real (kind=real_kind) :: temp(np,np,nlev)
  real (kind=real_kind) :: Jac(np,np,2*nlev,2*nlev), Q(np,np,2*nlev,2*nlev)
  real (kind=real_kind) :: R(np,np,2*nlev,2*nlev), Qt(2*nlev,2*nlev)
  real (kind=real_kind) :: e(np,np,2*nlev)
  real (kind=real_kind) :: Fn(np,np,2*nlev,1),x(np,np,2*nlev,1),epsie
  real (kind=real_kind) :: dFn(np,np,2*nlev,1)
  real (kind=real_kind) :: QtFn(np,np,2*nlev,1), Fntemp(2*nlev,1)
  real (kind=real_kind) :: res(np,np,2*nlev),resnorm,resnormmax
  real (kind=real_kind) :: statesavetemp(np,np,nlev,6),linsolveerror

  real (kind=real_kind) ::  itererr, itererrmax
  integer :: i,j,k,l,kptr,ie,itercount,itercountmax

  itercountmax=1
  itererrmax=0.d0
  resnormmax=0.d0

  epsie=1e-6
  call t_startf('compute_and_apply_rhs_imex_stiff')
  do ie=nets,nete 
    itercount=1
    itererr = 2.0*itertol       
    do while ((itercount < maxiter).and.((itererr > itertol).or.(resnorm > 1e-5)) )
      
      dp3d  => elem(ie)%state%dp3d(:,:,:,np1)
      theta_dp_cp  => elem(ie)%state%theta_dp_cp(:,:,:,np1)
      phi => elem(ie)%state%phi(:,:,:,np1)
      call get_kappa_star(kappa_star,elem(ie)%state%Qdp(:,:,:,1,qn0),dp3d)

      call get_pnh_and_exner(hvcoord,theta_dp_cp,dp3d,phi,elem(ie)%state%phis,&
      kappa_star,pnh,dpnh,exner) ! ,exner_i)
      dpnh_dp(:,:,:) = dpnh(:,:,:)/dp3d(:,:,:)
      
      Fn(:,:,1:nlev,1) = elem(ie)%state%w(:,:,:,np1)-statesave(ie,:,:,:,3) &
        +dt2*g*(1.0-dpnh_dp(:,:,:))
     
      Fn(:,:,nlev+1:2*nlev,1) = elem(ie)%state%phi(:,:,:,np1)-statesave(ie,:,:,:,4) &
        -dt2*g*elem(ie)%state%w(:,:,:,np1)      

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
      do k=1,2*nlev
           
        e(:,:,:)=0.0   
        e(:,:,k)=1.0  
           
       ! compute the new dpnh_dp at the perturbed values
       ! use the pointers only  
              
        phi(:,:,:)=phi(:,:,:)+epsie*e(:,:,nlev+1:2*nlev)
                 
        call get_kappa_star(kappa_star,elem(ie)%state%Qdp(:,:,:,1,qn0),dp3d)
                  
        call get_pnh_and_exner(hvcoord,theta_dp_cp,dp3d,phi,elem(ie)%state%phis,&     
          kappa_star,pnh,dpnh,exner) ! ,exner_i)
            
        dpnh_dp2(:,:,:) = dpnh(:,:,:)/dp3d(:,:,:)
        phi(:,:,:)=phi(:,:,:)-epsie*e(:,:,nlev+1:2*nlev)
              
       ! Form the approximate Jacobian
        Jac(:,:,1:nlev,k)=e(:,:,1:nlev)+g*dt2*(dpnh_dp(:,:,:)-dpnh_dp2(:,:,:))/epsie 
        Jac(:,:,nlev+1:2*nlev,k)=e(:,:,nlev+1:2*nlev)-g*dt2*e(:,:,1:nlev)
      end do
 
      call mgs(Jac,Q,R)
    
      do i=1,np
        do j=1,np
          Qt(:,:)=Q(i,j,:,:)
          Fntemp(:,1)=Fn(i,j,:,1)
          Qt=transpose(Qt)
          Fntemp=matmul(Qt,Fntemp)
          QtFn(i,j,:,1) = Fntemp(:,1)
        end do
      end do               
                        
      call backsubstitution(R,-QtFn,x,linsolveerror)
                                  
      elem(ie)%state%w(:,:,1:nlev,np1) = elem(ie)%state%w(:,:,1:nlev,np1)     &
        + x(:,:,1:nlev,1)
      elem(ie)%state%phi(:,:,1:nlev,np1) = elem(ie)%state%phi(:,:,1:nlev,np1) &
        + x(:,:,nlev+1:2*nlev,1)
      itererr=norm2(x)
      resnorm=norm2(Fn)
      itercount=itercount+1
     end do                   

      if (itercount > itercountmax) then
        itercountmax=itercount
      endif
      if (itererr > itererrmax) then 
        itererrmax=itererr
      end if 
      if (resnorm > resnormmax) then 
        resnormmax = resnorm
      end if 
  end do
  maxiter=itercountmax
  if (itererrmax > resnormmax) then 
    itertol=itererrmax
  else
    itertol= resnormmax
  end if 

! now compte the apply the boundary exchange 
! ==========================================  
  do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        elem(ie)%state%dp3d(:,:,k,np1) =elem(ie)%spheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
        elem(ie)%state%theta_dp_cp(:,:,k,np1)=elem(ie)%spheremp(:,:)*elem(ie)%state%theta_dp_cp(:,:,k,np1)
        elem(ie)%state%w(:,:,k,np1)    =elem(ie)%spheremp(:,:)*elem(ie)%state%w(:,:,k,np1)
        elem(ie)%state%phi(:,:,k,np1)  =elem(ie)%spheremp(:,:)*elem(ie)%state%phi(:,:,k,np1)
        elem(ie)%state%v(:,:,1,k,np1)  =elem(ie)%spheremp(:,:)*elem(ie)%state%v(:,:,1,k,np1)
        elem(ie)%state%v(:,:,2,k,np1)  =elem(ie)%spheremp(:,:)*elem(ie)%state%v(:,:,2,k,np1)
     end do

     kptr=0
     call edgeVpack(edge6, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr, ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%theta_dp_cp(:,:,:,np1),nlev,kptr,ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%w(:,:,:,np1),nlev,kptr,ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%phi(:,:,:,np1),nlev,kptr,ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,ie)
  end do
  call t_startf('caar_bexchV')
  call bndry_exchangeV(hybrid,edge6)
  call t_stopf('caar_bexchV')
  do ie=nets,nete
     kptr=0
     call edgeVunpack(edge6, elem(ie)%state%dp3d(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%theta_dp_cp(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%w(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%phi(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, ie)

     
     ! ====================================================
     ! Scale tendencies by inverse mass matrix
     ! ====================================================
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        elem(ie)%state%dp3d(:,:,k,np1) =elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
        elem(ie)%state%theta_dp_cp(:,:,k,np1)=elem(ie)%rspheremp(:,:)*elem(ie)%state%theta_dp_cp(:,:,k,np1)
        elem(ie)%state%w(:,:,k,np1)    =elem(ie)%rspheremp(:,:)*elem(ie)%state%w(:,:,k,np1)
        elem(ie)%state%phi(:,:,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%phi(:,:,k,np1)
        elem(ie)%state%v(:,:,1,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,1,k,np1)
        elem(ie)%state%v(:,:,2,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,2,k,np1)
     end do
  end do

  call t_stopf('compute_and_apply_rhs_imex_stiff')

  end subroutine compute_and_apply_rhs_imex_stiff



  subroutine compute_stage_value_dirk_stiff(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w,maxiter,itertol)
  !===================================================================================
  ! this subroutine solves a stage value equation for a DIRK method which takes the form 
  ! gi = un0 + dt* sum(1:i-1)(aij n(gj)+a2ij s(gj)) + dt *a2ii s(gi) := y + dt a2ii s(gi)
  ! It is assumed that un0 has the value of y and the computed value of gi is stored at 
  ! unp1
  !===================================================================================
  use kinds, only : real_kind
  use derivative_mod, only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use derivative_mod, only : subcell_div_fluxes, subcell_dss_fluxes
  use edge_mod, only : edgevpack, edgevunpack, edgeDGVunpack
  use edgetype_mod, only : edgedescriptor_t
  use bndry_mod, only : bndry_exchangev
  use control_mod, only : moisture, qsplit, use_cpstar, rsplit, swest
  use hybvcoord_mod, only : hvcoord_t

  use physical_constants, only : cp, cpwater_vapor, Rgas, kappa, Rwater_vapor,p0, g
  use physics_mod, only : virtual_specific_heat, virtual_temperature
  use prim_si_mod, only : preq_vertadv_v, preq_vertadv_upwind, preq_omega_ps, preq_hydrostatic_v2

  use time_mod, only : tevolve

  implicit none
  integer, intent(in) :: np1,nm1,n0,qn0,nets,nete
  real*8, intent(in) :: dt2
  logical, intent(in)  :: compute_diagnostics
  integer :: maxiter
  real*8 :: itertol

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind) :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux

  ! local
  real (kind=real_kind), pointer, dimension(:,:,:)   :: phi
  real (kind=real_kind), pointer, dimension(:,:,:)   :: dp3d
  real (kind=real_kind), pointer, dimension(:,:,:)   :: theta_dp_cp
   
  real (kind=real_kind) :: kappa_star(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev)     ! nh (nonydro) pressure
  real (kind=real_kind) :: dpnh(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)     ! exner nh pressure
  real (kind=real_kind) :: dpnh_dp(np,np,nlev),dpnh_dp2(np,np,nlev)    !    ! dpnh / dp3d  
  real (kind=real_kind) :: temp(np,np,nlev)
  real (kind=real_kind) :: Jac(np,np,2*nlev,2*nlev), Q(np,np,2*nlev,2*nlev)
  real (kind=real_kind) :: R(np,np,2*nlev,2*nlev), Qt(2*nlev,2*nlev)
  real (kind=real_kind) :: e(np,np,2*nlev)
  real (kind=real_kind) :: Fn(np,np,2*nlev,1),x(np,np,2*nlev,1),epsie
  real (kind=real_kind) :: dFn(np,np,2*nlev,1)
  real (kind=real_kind) :: QtFn(np,np,2*nlev,1), Fntemp(2*nlev,1)
  real (kind=real_kind) :: res(np,np,2*nlev),resnorm,resnormmax
  real (kind=real_kind) :: stagevaluesum(np,np,nlev,6),linsolveerror

  real (kind=real_kind) ::  itererr, itererrmax
  integer :: i,j,k,l,kptr,ie,itercount,itercountmax

  itercountmax=1
  itererrmax=0.d0
  resnormmax=0.d0

  epsie=1e-4
  call t_startf('compute_stage_value_dirk_stiff')

  do ie=nets,nete 

    itercount=1
    itererr = 2.0*itertol      

    do while ((itercount < maxiter).and.((itererr > itertol).or.(resnorm > itertol*1.d2)) )
  
      dp3d  => elem(ie)%state%dp3d(:,:,:,np1)
      theta_dp_cp  => elem(ie)%state%theta_dp_cp(:,:,:,np1)
      phi => elem(ie)%state%phi(:,:,:,np1)
           
      call get_kappa_star(kappa_star,elem(ie)%state%Qdp(:,:,:,1,qn0),dp3d)
          
      call get_pnh_and_exner(hvcoord,theta_dp_cp,dp3d,phi,elem(ie)%state%phis,&
      kappa_star,pnh,dpnh,exner) ! ,exner_i)
           
      dpnh_dp(:,:,:) = dpnh(:,:,:)/dp3d(:,:,:)
                
      Fn(:,:,1:nlev,1) = elem(ie)%state%w(:,:,:,np1)-elem(ie)%state%w(:,:,:,n0) &
        +dt2*g*(1.0-dpnh_dp(:,:,:))
                
      Fn(:,:,nlev+1:2*nlev,1) = elem(ie)%state%phi(:,:,:,np1)-elem(ie)%state%phi(:,:,:,n0) &
        -dt2*g*elem(ie)%state%w(:,:,:,np1)      
      
!      do i=1,np
!        do j=1,np
!          do l=1,2*nlev
!            print *, 'bad NaN'
!            print *, 'Fn(i,j,k,l,1) ', Fn(i,j,l,1)
!          end do
!        end do
!     end do
  
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
      do k=1,2*nlev
           
        e(:,:,:)=0.0   
        e(:,:,k)=1.0  
           
       ! compute the new dpnh_dp at the perturbed values
       ! use the pointers only  
              
        phi(:,:,:)=phi(:,:,:)+epsie*e(:,:,nlev+1:2*nlev)
                 
        call get_kappa_star(kappa_star,elem(ie)%state%Qdp(:,:,:,1,qn0),dp3d)
                  
        call get_pnh_and_exner(hvcoord,theta_dp_cp,dp3d,phi,elem(ie)%state%phis,&     
          kappa_star,pnh,dpnh,exner) ! ,exner_i)
            
        dpnh_dp2(:,:,:) = dpnh(:,:,:)/dp3d(:,:,:)
        phi(:,:,:)=phi(:,:,:)-epsie*e(:,:,nlev+1:2*nlev)
              
       ! Form the approximate Jacobian
        Jac(:,:,1:nlev,k)=e(:,:,1:nlev)+dt2*g*(dpnh_dp(:,:,:)-dpnh_dp2(:,:,:))/epsie 
        Jac(:,:,nlev+1:2*nlev,k)=e(:,:,nlev+1:2*nlev)-dt2*g*e(:,:,1:nlev)
      end do
 
      call mgs(Jac,Q,R)
    
      do i=1,np
        do j=1,np
          Qt(:,:)=Q(i,j,:,:)
          Fntemp(:,1)=Fn(i,j,:,1)
          Qt=transpose(Qt)
          Fntemp=matmul(Qt,Fntemp)
          QtFn(i,j,:,1) = Fntemp(:,1)
        end do
      end do               
                        
      call backsubstitution(R,-QtFn,x,linsolveerror)
      elem(ie)%state%w(:,:,:,np1)   = elem(ie)%state%w(:,:,:,np1) + x(:,:,1:nlev,1)                          
      elem(ie)%state%phi(:,:,:,np1) = elem(ie)%state%phi(:,:,:,np1) + x(:,:,nlev+1:2*nlev,1)                          
      itererr=norm2(x)
      resnorm=norm2(Fn)
      itercount=itercount+1
!      if ( (isnan(itererr)).or.(isnan(resnorm))) then 
!        print *, 'bad NaN'
!        print *, 'itererr ', itererr
!        print *, 'resnorm ', resnorm
!        stop
!      end if
!      print *, itercount, itererr,resnorm
     end do                   

      if (itercount > itercountmax) then
        itercountmax=itercount
      endif
      if (itererr > itererrmax) then 
        itererrmax=itererr
      end if 
      if (resnorm > resnormmax) then 
        resnormmax = resnorm
      end if 
  end do
  maxiter=itercountmax
  if (itererrmax > resnormmax) then 
    itertol=itererrmax
  else
    itertol= resnormmax
  end if 

! now compute the apply the boundary exchange 
! ==========================================  
  do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        elem(ie)%state%dp3d(:,:,k,np1) =elem(ie)%spheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
        elem(ie)%state%theta_dp_cp(:,:,k,np1)=elem(ie)%spheremp(:,:)*elem(ie)%state%theta_dp_cp(:,:,k,np1)
        elem(ie)%state%w(:,:,k,np1)    =elem(ie)%spheremp(:,:)*elem(ie)%state%w(:,:,k,np1)
        elem(ie)%state%phi(:,:,k,np1)  =elem(ie)%spheremp(:,:)*elem(ie)%state%phi(:,:,k,np1)
        elem(ie)%state%v(:,:,1,k,np1)  =elem(ie)%spheremp(:,:)*elem(ie)%state%v(:,:,1,k,np1)
        elem(ie)%state%v(:,:,2,k,np1)  =elem(ie)%spheremp(:,:)*elem(ie)%state%v(:,:,2,k,np1)
     end do

     kptr=0
     call edgeVpack(edge6, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr, ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%theta_dp_cp(:,:,:,np1),nlev,kptr,ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%w(:,:,:,np1),nlev,kptr,ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%phi(:,:,:,np1),nlev,kptr,ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,ie)
  end do
  call t_startf('caar_bexchV')
  call bndry_exchangeV(hybrid,edge6)
  call t_stopf('caar_bexchV')
  do ie=nets,nete
     kptr=0
     call edgeVunpack(edge6, elem(ie)%state%dp3d(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%theta_dp_cp(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%w(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%phi(:,:,:,np1), nlev, kptr, ie)
     kptr=kptr+nlev
     call edgeVunpack(edge6, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, ie)

     
     ! ====================================================
     ! Scale tendencies by inverse mass matrix
     ! ====================================================
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        elem(ie)%state%dp3d(:,:,k,np1) =elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
        elem(ie)%state%theta_dp_cp(:,:,k,np1)=elem(ie)%rspheremp(:,:)*elem(ie)%state%theta_dp_cp(:,:,k,np1)
        elem(ie)%state%w(:,:,k,np1)    =elem(ie)%rspheremp(:,:)*elem(ie)%state%w(:,:,k,np1)
        elem(ie)%state%phi(:,:,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%phi(:,:,k,np1)
        elem(ie)%state%v(:,:,1,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,1,k,np1)
        elem(ie)%state%v(:,:,2,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,2,k,np1)
     end do
  end do

  call t_stopf('compute_stage_value_dirk_stiff')

  end subroutine compute_stage_value_dirk_stiff


 
!==================================================================================
!==================================================================================
!==================================================================================

 subroutine fp_iteration_impeuler(elem,np1,n0,nm1,qn0,hvcoord,dt,hybrid,&
            deriv,nets,nete,eta_ave_w,itertol,maxiter,itererr,itercount,imex)
    use reduction_mod,  only: parallelmax

    implicit none
    integer, intent(in) :: np1,n0,nm1,qn0,nets,nete,maxiter
    integer, intent(inout) :: itercount
    real*8, intent(inout) :: itererr
    real*8, intent(in) :: dt,itertol
    logical, intent(in) :: imex

    type (hvcoord_t)     , intent(in) :: hvcoord
    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    type (derivative_t)  , intent(in) :: deriv
    real (kind=real_kind) :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux

    integer :: ie,i,j,k
    real (kind=real_kind) :: delta(nets:nete,np,np,nlev,6)
    real (kind=real_kind) :: itererrvec(nets:nete,6),itererrnorm(6)
    real (kind=real_kind) :: statesave(nets:nete,np,np,nlev,6),dampfac

    dampfac=1d-2
    do ie=nets,nete
      elem(ie)%state%v(:,:,1,:,n0)         = elem(ie)%state%v(:,:,1,:,n0)         * dampfac
      elem(ie)%state%v(:,:,2,:,n0)         = elem(ie)%state%v(:,:,2,:,n0)         * dampfac
      elem(ie)%state%w(:,:,:,n0)           = elem(ie)%state%w(:,:,:,n0)           * dampfac
      elem(ie)%state%phi(:,:,:,n0)         = elem(ie)%state%phi(:,:,:,n0)         * dampfac
      elem(ie)%state%theta_dp_cp(:,:,:,n0) = elem(ie)%state%theta_dp_cp(:,:,:,n0) * dampfac
      elem(ie)%state%dp3d(:,:,:,n0)        = elem(ie)%state%dp3d(:,:,:,n0)        * dampfac
    end do
    do while (( itercount < maxiter ).and.( itererr > itertol ))
       do ie=nets,nete
         delta(ie,:,:,:,1) = elem(ie)%state%v(:,:,1,:,np1)
         delta(ie,:,:,:,2) = elem(ie)%state%v(:,:,2,:,np1)
         delta(ie,:,:,:,3) = elem(ie)%state%w(:,:,:,np1)
         delta(ie,:,:,:,4) = elem(ie)%state%phi(:,:,:,np1)
         delta(ie,:,:,:,5) = elem(ie)%state%theta_dp_cp(:,:,:,np1)
         delta(ie,:,:,:,6) = elem(ie)%state%dp3d(:,:,:,np1)
       end do
      if (imex) then 
!        call compute_and_apply_rhs_imex_stiff(np1,n0,np1,qn0,dampfac*dt,elem,hvcoord,hybrid,&
!          deriv,nets,nete,.false.,eta_ave_w,maxiter,itertol,1d0)
     call compute_and_apply_rhs_imex(np1,n0,np1,qn0,dampfac*dt,elem,hvcoord,hybrid,&
          deriv,nets,nete,.false.,eta_ave_w,.true.)
      else
       call compute_and_apply_rhs(np1,n0,np1,qn0,dampfac*dt,elem,hvcoord,hybrid,&
          deriv,nets,nete,.false.,eta_ave_w)        
      end if
     
      do ie=nets,nete
         elem(ie)%state%w(:,:,:,np1)           = (1-dampfac) * delta(ie,:,:,:,3) + elem(ie)%state%w(:,:,:,np1)
         elem(ie)%state%phi(:,:,:,np1)         = (1-dampfac) * delta(ie,:,:,:,4) + elem(ie)%state%phi(:,:,:,np1)
         if (imex) then 
        
         else
           elem(ie)%state%v(:,:,1,:,np1)         = (1-dampfac) * delta(ie,:,:,:,1) + elem(ie)%state%v(:,:,1,:,np1)
           elem(ie)%state%v(:,:,2,:,np1)         = (1-dampfac) * delta(ie,:,:,:,2) + elem(ie)%state%v(:,:,2,:,np1)
           elem(ie)%state%theta_dp_cp(:,:,:,np1) = (1-dampfac) * delta(ie,:,:,:,5) + elem(ie)%state%theta_dp_cp(:,:,:,np1)
           elem(ie)%state%dp3d(:,:,:,np1)        = (1-dampfac) * delta(ie,:,:,:,6) + elem(ie)%state%dp3d(:,:,:,np1)
         end if
      end do

      do ie=nets,nete
         delta(ie,:,:,:,1) = delta(ie,:,:,:,1) - elem(ie)%state%v(:,:,1,:,np1)
         delta(ie,:,:,:,2) = delta(ie,:,:,:,2) - elem(ie)%state%v(:,:,2,:,np1)
         delta(ie,:,:,:,3) = delta(ie,:,:,:,3) - elem(ie)%state%w(:,:,:,np1)
         delta(ie,:,:,:,4) = delta(ie,:,:,:,4) - elem(ie)%state%phi(:,:,:,np1)
         delta(ie,:,:,:,5) = delta(ie,:,:,:,5) - elem(ie)%state%theta_dp_cp(:,:,:,np1)
         delta(ie,:,:,:,6) = delta(ie,:,:,:,6) - elem(ie)%state%dp3d(:,:,:,np1)
       end do

       do ie=nets,nete
         itererrvec(ie,1)=MAXVAL(abs(delta(ie,:,:,:,1)))
         itererrvec(ie,2)=MAXVAL(abs(delta(ie,:,:,:,2)))
         itererrvec(ie,3)=MAXVAL(abs(delta(ie,:,:,:,3)))
         itererrvec(ie,4)=MAXVAL(abs(delta(ie,:,:,:,4)))
         itererrvec(ie,5)=MAXVAL(abs(delta(ie,:,:,:,5)))
         itererrvec(ie,6)=MAXVAL(abs(delta(ie,:,:,:,6)))
       end do
       do i=1,6
         itererrnorm(i)=ParallelMax(abs(itererrvec(:,i)),hybrid)
       end do
       itererr=MAXVAL(itererrnorm(:))
       itercount=itercount+1
       print *, 'maxval v1 ', 'maxval v2 ', 'maxval w '
       print *, itererrnorm(1), itererrnorm(2), itererrnorm(3)      
       print *, 'maxval phi ', 'maxval theta ', 'maxval dp3d '
       print *, itererrnorm(4), itererrnorm(5), itererrnorm(6)
       print *, 'nu. of iter. ' , 'appr. err. '
       print *, itercount, itererr

    end do
    
    do ie=nets,nete
      elem(ie)%state%v(:,:,1,:,n0)         = elem(ie)%state%v(:,:,1,:,n0)         / dampfac
      elem(ie)%state%v(:,:,2,:,n0)         = elem(ie)%state%v(:,:,2,:,n0)         / dampfac
      elem(ie)%state%w(:,:,:,n0)           = elem(ie)%state%w(:,:,:,n0)           / dampfac
      elem(ie)%state%phi(:,:,:,n0)         = elem(ie)%state%phi(:,:,:,n0)         / dampfac
      elem(ie)%state%theta_dp_cp(:,:,:,n0) = elem(ie)%state%theta_dp_cp(:,:,:,n0) / dampfac
      elem(ie)%state%dp3d(:,:,:,n0)        = elem(ie)%state%dp3d(:,:,:,n0)        / dampfac
    end do
    
  end subroutine


  subroutine mgs(A,Q,R)
    
    real (kind=real_kind), intent(inout) :: A(np,np,2*nlev,2*nlev)
    real (kind=real_kind), intent(inout) :: Q(np,np,2*nlev,2*nlev)
    real (kind=real_kind), intent(inout) :: R(np,np,2*nlev,2*nlev)

    ! local variables
    integer :: i,j,k,l
    real (kind=real_kind) :: Atemp(2*nlev,2*nlev)
    real (kind=real_kind) :: Qtemp(2*nlev,2*nlev)
    real (kind=real_kind) :: Rtemp(2*nlev,2*nlev)

    R=0.0
    do i=1,np
      do j=1,np
        Atemp(:,:)=A(i,j,:,:)
        Rtemp=0.0
        do k=1,2*nlev
          Rtemp(k,k)=norm2(Atemp(:,k))
          Qtemp(1:2*nlev,k)=Atemp(1:2*nlev,k)/Rtemp(k,k)
          do l=k+1,2*nlev
            Rtemp(k,l)=dot_product(Qtemp(:,k),Atemp(:,l))
            Atemp(:,l)=Atemp(:,l)-Rtemp(k,l)*Qtemp(:,k)
          end do
        end do
        A(i,j,:,:)=Atemp(:,:)
        Q(i,j,:,:)=Qtemp(:,:)
        R(i,j,:,:)=Rtemp(:,:)
      end do
    end do
  end subroutine



  subroutine backsubstitution(R,b,x,err)
    
    real (kind=real_kind), intent(in) :: R(np,np,2*nlev,2*nlev)
    real (kind=real_kind), intent(in) :: b(np,np,2*nlev)
    real (kind=real_kind), intent(inout) :: x(np,np,2*nlev)
    real (kind=real_kind), intent(inout) :: err

    integer :: i,j
    real (kind=real_kind) :: sum(np,np),error,errortemp
    real (kind=real_kind) :: Rtemp(2*nlev,2*nlev)
    real (kind=real_kind) :: btemp(2*nlev,1)
    real (kind=real_kind) :: xtemp(2*nlev,1)


    do i=2*nlev,1,-1
      sum(:,:) = b(:,:,i)
      do j=i+1,2*nlev
        sum(:,:)=sum(:,:)-R(:,:,i,j)*x(:,:,j)
      end do
      x(:,:,i)=sum(:,:)/R(:,:,i,i)
    end do
    err=0.0
    do i=1,np
      do j=1,np
        Rtemp(:,:)=R(i,j,:,:)
        xtemp(:,1)=x(i,j,:)
        btemp(:,1)=b(i,j,:)
        errortemp = norm2(matmul(Rtemp(:,:),xtemp(:,1))-btemp(:,1))
        if (errortemp > err) then
          err=errortemp
        end if
      end do
   end do
  
  end subroutine

  subroutine state_save(elem,state,n,nets,nete)

    type (element_t),			intent(inout), target :: elem(:)!   
    real (kind=real_kind), intent(inout) :: state(nets:nete,np,np,nlev,6)
    integer                             :: n,nets,nete

    integer :: ie
    do ie=nets,nete
      state(ie,:,:,:,1) = elem(ie)%state%v(:,:,1,:,n)
      state(ie,:,:,:,2) = elem(ie)%state%v(:,:,2,:,n)
      state(ie,:,:,:,3) = elem(ie)%state%w(:,:,:,n)
      state(ie,:,:,:,4) = elem(ie)%state%phi(:,:,:,n)
      state(ie,:,:,:,5) = elem(ie)%state%theta_dp_cp(:,:,:,n)
      state(ie,:,:,:,6) = elem(ie)%state%dp3d(:,:,:,n)
    end do

  end subroutine state_save

  subroutine state_read(elem,state,n,nets,nete)
   
    type (element_t),			intent(inout), target :: elem(:)!   
    real (kind=real_kind), intent(inout) :: state(nets:nete,np,np,nlev,6)
    integer                            :: n,nets,nete

    integer :: ie
    do ie=nets,nete
      elem(ie)%state%v(:,:,1,:,n)         = state(ie,:,:,:,1)
      elem(ie)%state%v(:,:,2,:,n)         = state(ie,:,:,:,2)
      elem(ie)%state%w(:,:,:,n)           = state(ie,:,:,:,3)
      elem(ie)%state%phi(:,:,:,n)         = state(ie,:,:,:,4)
      elem(ie)%state%theta_dp_cp(:,:,:,n) = state(ie,:,:,:,5)
      elem(ie)%state%dp3d(:,:,:,n)        = state(ie,:,:,:,6)
    end do

  end subroutine state_read


  subroutine state_add(elem,state,n,nets,nete,scale1,scale2)

    type (element_t),			intent(inout), target :: elem(:)
    real (kind=real_kind), intent(inout) :: state(nets:nete,np,np,nlev,6)
    real (kind=real_kind), intent(in)    ::  scale1,scale2
    integer                             :: n,nets,nete

    integer :: ie

    do ie=nets,nete
      elem(ie)%state%v(:,:,1,:,n)         = state(ie,:,:,:,1)*scale1 &
        + elem(ie)%state%v(:,:,1,:,n)*scale2
      elem(ie)%state%v(:,:,2,:,n)         = state(ie,:,:,:,2)*scale1 & 
        + elem(ie)%state%v(:,:,2,:,n)*scale2
      elem(ie)%state%w(:,:,:,n)           = state(ie,:,:,:,3)*scale1 &
        + elem(ie)%state%w(:,:,:,n)*scale2
      elem(ie)%state%phi(:,:,:,n)         = state(ie,:,:,:,4)*scale1 &
        + elem(ie)%state%phi(:,:,:,n)*scale2
      elem(ie)%state%theta_dp_cp(:,:,:,n) = state(ie,:,:,:,5)*scale1 &
        + elem(ie)%state%theta_dp_cp(:,:,:,n)*scale2
      elem(ie)%state%dp3d(:,:,:,n)        = state(ie,:,:,:,6)*scale1 &
        + elem(ie)%state%dp3d(:,:,:,n)*scale2
    end do

  end subroutine state_add

end module prim_advance_mod

