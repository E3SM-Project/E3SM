#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
!
!
!  Man dynamics routines for "theta" nonhydrostatic model
!  Original version: Mark Taylor 2017/1
!
module prim_advance_mod

  use bndry_mod,          only: bndry_exchangev
  use control_mod,        only: dcmip16_mu, dcmip16_mu_s, hypervis_order, hypervis_subcycle,&
    integration, nu, nu_div, nu_p, nu_s, nu_top, prescribed_wind, qsplit, rsplit, test_case,&
    theta_hydrostatic_mode, tstep_type, use_moisture, use_cpstar
  use derivative_mod,     only: derivative_t, divergence_sphere, gradient_sphere, laplace_sphere_wk,&
    laplace_z, vorticity_sphere, vlaplace_sphere_wk 
  use derivative_mod,     only: subcell_div_fluxes, subcell_dss_fluxes
  use dimensions_mod,     only: max_corner_elem, nelemd, nlev, nlevp, np, qsize
  use edge_mod,           only: edgeDGVunpack, edgevpack, edgevunpack, initEdgeBuffer
  use edgetype_mod,       only: EdgeBuffer_t,  EdgeDescriptor_t, edgedescriptor_t
  use element_mod,        only: element_t
  use element_ops,        only: copy_state, get_cp_star, get_kappa_star, &
    get_temperature, set_theta_ref, state0
  use eos,                only: get_pnh_and_exner,get_dry_phinh,get_dirk_jacobian
  use hevi_mod,           only: backsubstitution, elemstate_add, mgs, state_save,state_read
  use hybrid_mod,         only: hybrid_t
  use hybvcoord_mod,      only: hvcoord_t
  use kinds,              only: iulog, real_kind
  use perf_mod,           only: t_adj_detailf, t_barrierf, t_startf, t_stopf ! _EXTERNAL
  use parallel_mod,       only: abortmp, global_shared_buf, global_shared_sum, iam, parallel_t
  use physical_constants, only: Cp, cp, cpwater_vapor, g, kappa, Rgas, Rwater_vapor, p0 
  use physics_mod,        only: virtual_specific_heat, virtual_temperature
  use prim_si_mod,        only: preq_vertadv_upwind, preq_vertadv_v, preq_hydrostatic_v2, preq_omega_ps
  use reduction_mod,      only: parallelmax, reductionbuffer_ordered_1d_t
  use time_mod,           only: timelevel_qdp, timelevel_t
  use test_mod,           only: set_prescribed_wind
  use viscosity_theta,    only: biharmonic_wk_theta

#ifdef TRILINOS
    use prim_derived_type_mod ,only : derived_type, initialize
    use, intrinsic :: iso_c_binding
#endif
 
  implicit none
  private
  save
  public :: prim_advance_exp, prim_advance_init1, &
       applyCAMforcing_dp3d, applyCAMforcing_ps

!  type (EdgeBuffer_t) :: edge5
  type (EdgeBuffer_t) :: edge6
  real (kind=real_kind), allocatable :: ur_weights(:)

contains





  subroutine prim_advance_init1(par, elem,integration)
        
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





  !_____________________________________________________________________
  subroutine prim_advance_exp(elem, deriv, hvcoord, hybrid,dt, tl,  nets, nete, compute_diagnostics)

    type (element_t),      intent(inout), target :: elem(:)
    type (derivative_t),   intent(in)            :: deriv
    type (hvcoord_t)                             :: hvcoord
    type (hybrid_t),       intent(in)            :: hybrid
    real (kind=real_kind), intent(in)            :: dt
    type (TimeLevel_t)   , intent(in)            :: tl
    integer              , intent(in)            :: nets
    integer              , intent(in)            :: nete
    logical,               intent(in)            :: compute_diagnostics

    real (kind=real_kind) :: dt2, time, dt_vis, x, eta_ave_w
    real (kind=real_kind) :: itertol,statesave(nets:nete,np,np,nlev,6)
    real (kind=real_kind) :: a1,a2,a3,a4,a5,ahat1,ahat2,ahat3,ahat4,ahat5,dhat1,dhat2
    real (kind=real_kind) :: statesave0(nets:nete,np,np,nlev,6)
    real (kind=real_kind) :: statesave2(nets:nete,np,np,nlev,6)
    real (kind=real_kind) :: statesave3(nets:nete,np,np,nlev,6)

    real (kind=real_kind) ::  gamma,delta

    integer :: ie,nm1,n0,np1,nstep,qsplit_stage,k, qn0
    integer :: n,i,j,maxiter
 

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
!
!   tstep_type=5  Kinnmark&Gray RK3 5 stage 3rd order            CFL=3.87  (sqrt(15))
!                 From Paul Ullrich.  3rd order for nonlinear terms also
!                 K&G method is only 3rd order for linear
!                 optimal: for windspeeds ~120m/s,gravity: 340m/2
!                 run with qsplit=1
!                 (K&G 2nd order method has CFL=4. tiny CFL improvement not worth 2nd order)
!   tstep_type=6  KG with BW Euler implicit step, usful as a debug
!   tstep_type=7  ARS232 ARK-IMEX method with 3 explicit stages and 2 implicit stages, 2nd order 
!                 accurate with stage order 1
!   

! default weights for computing mean dynamics fluxes
    eta_ave_w = 1d0/qsplit
 
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
    if (tstep_type==1) then 
       ! RK2                                                                                                              
       ! forward euler to u(dt/2) = u(0) + (dt/2) RHS(0)  (store in u(np1))                                               
       call t_startf("RK2_timestep")                                                                                      
       call compute_andor_apply_rhs(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&                                              
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)                                                      
       ! leapfrog:  u(dt) = u(0) + dt RHS(dt/2)     (store in u(np1))                                                     
       call compute_andor_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&                                               
            deriv,nets,nete,.false.,eta_ave_w,1.d0,1.d0,1.d0)                                                             
       call t_stopf("RK2_timestep")   
    else if (tstep_type==5) then
       ! Ullrich 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
       ! u1 = u0 + dt/5 RHS(u0)  (save u1 in timelevel nm1)
       call t_startf("U3-5stage_timestep")
       call compute_andor_apply_rhs(nm1,n0,n0,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w/4,1.d0,1.d0,1.d0)
       ! u2 = u0 + dt/5 RHS(u1)
       call compute_andor_apply_rhs(np1,n0,nm1,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,1.d0,1.d0)
       ! u3 = u0 + dt/3 RHS(u2)
       call compute_andor_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,1.d0,1.d0)
       ! u4 = u0 + 2dt/3 RHS(u3)
       call compute_andor_apply_rhs(np1,n0,np1,qn0,2*dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,1.d0,1.d0)
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
          elem(ie)%state%phinh(:,:,:,nm1)= (5*elem(ie)%state%phinh(:,:,:,nm1) &
                  - elem(ie)%state%phinh(:,:,:,n0) )/4
       enddo
       ! u5 = (5*u1/4 - u0/4) + 3dt/4 RHS(u4)
       call compute_andor_apply_rhs(np1,nm1,np1,qn0,3*dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,3*eta_ave_w/4,1.d0,1.d0,1.d0)
       ! final method is the same as:
       ! u5 = u0 +  dt/4 RHS(u0)) + 3dt/4 RHS(u4)
       call t_stopf("U3-5stage_timestep")  
 ! ==============================================================================
    else if (tstep_type==6) then ! Imex hevi, implicit euler after the full explicit time-step
    ! it seems to run with ne=16, nlev=26, dt=200 for JW Baro up to 18.8 days at least
    ! Ullrich 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
       ! u1 = u0 + dt/5 RHS(u0)  (save u1 in timelevel nm1)
       call t_startf("U3-5stage_timestep")
       call compute_andor_apply_rhs(nm1,n0,n0,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w/4,1.d0,0.d0,1.d0)
       ! u2 = u0 + dt/5 RHS(u1)
       call compute_andor_apply_rhs(np1,n0,nm1,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,0.d0,1.d0)
       ! u3 = u0 + dt/3 RHS(u2)
       call compute_andor_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,0.d0,1.d0)
       ! u4 = u0 + 2dt/3 RHS(u3)
       call compute_andor_apply_rhs(np1,n0,np1,qn0,2*dt/3,elem,hvcoord,hybrid,&
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
          elem(ie)%state%phinh(:,:,:,nm1)= (5*elem(ie)%state%phinh(:,:,:,nm1) &
                  - elem(ie)%state%phinh(:,:,:,n0) )/4
       enddo
       ! u5 = (5*u1/4 - u0/4) + 3dt/4 RHS(u4)
       call compute_andor_apply_rhs(np1,nm1,np1,qn0,3*dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,3*eta_ave_w/4,1.d0,0.d0,1.d0)
 
      maxiter=4
      itertol=1e-8
      call elemstate_add(elem,statesave,nets,nete,1,np1,n0,n0,1d0,0d0,0d0)
      call compute_stage_value_dirk(np1,qn0,dt,elem,hvcoord,hybrid,&
       deriv,nets,nete,maxiter,itertol)
       call t_stopf("U3-5stage_timestep")
!============================================================================================
    else if (tstep_type==7) then ! ARS232 from (Ascher et al., 1997), nh-imex
      ! ARS232 is 2nd order, stage order 1, DIRK scheme is A-stable and L-stable
      ! 2 implicit solves and 3 stages total
      call t_startf("ARS232_timestep")
      delta = -2.d0*sqrt(2.d0)/3.d0
      gamma = 1.d0 - 1.d0/sqrt(2.d0)

      ! save un0 as statesave
      call state_save(elem,statesave,n0,nets,nete)
                               
      ! compute dt*n(un0)=dt*n(g1) and save at np1
      call compute_andor_apply_rhs(np1,n0,n0,qn0,dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0.d0,1.d0,0.d0,0.d0)      

      ! form un0+dt*gamma*n(g1) and store at n0
      call elemstate_add(elem,statesave,nets,nete,1,n0,np1,n0,gamma,1.d0,0.d0)
                             
      maxiter=10
      itertol=1e-15
      ! solve g2 = un0 + dt*gamma*n(g1)+dt*gamma*s(g2) for g2 and save at nm1
      call elemstate_add(elem,statesave,nets,nete,1,nm1,n0,n0,1d0,0d0,0d0)
      call compute_stage_value_dirk(nm1,qn0,gamma*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
!      print *, 'num iters  ', maxiter
!=== End of Phase 1 ====
! at this point, g2 is at nm1, un0+dt*gamma*n(g1) is at n0, and dt*n(g1) is at np1
                
      ! Form dt*n(g2) and store at np1
      call compute_andor_apply_rhs(np1,nm1,nm1,qn0,dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,.false.,(1-gamma)*eta_ave_w,1.d0,0.d0,0.d0)

      ! solve g2 = un0 + dt*gamma*n(g1) + dt*gamma*s(g2) for dt*s(g2) and 
      ! store at nm1
      call elemstate_add(elem,statesave,nets,nete,1,nm1,nm1,n0,1.d0/gamma,-1.d0/gamma,0.d0)

      ! Form dt*gamma*n(g1) and store at n0
      call elemstate_add(elem,statesave,nets,nete,2,n0,n0,n0,1.d0,1.d0,-1.d0)                  

      ! Form un0+dt*delta*n(g1) and store at n0
      call elemstate_add(elem,statesave,nets,nete,2,n0,n0,n0,delta/gamma,1.d0,1.d0)                  

      ! Form un0+dt*delta*n(g1)+dt*(1-delta)*n(g2)+dt*(1-gamma)*n(g3)
      call elemstate_add(elem,statesave,nets,nete,4,n0,np1,nm1,1.d0-delta,1.d0-gamma,0.d0)
      
      ! form un0+dt*(1-gamma)*(n(g2)+s(g2)) at nm1
      call elemstate_add(elem,statesave,nets,nete,3,nm1,np1,nm1,1.d0-gamma,1.d0-gamma,1.d0)
                       
      maxiter=10
      itertol=1e-15
      !	solve g3 = (un0+dt*delta*n(g1))+dt*(1-delta)*n(g2)+dt*(1-gamma)*s(g2)+dt*gamma*s(g3)
      ! for g3 using (un0+dt*delta*n(g1))+dt*(1-delta)*n(g2)+dt*(1-gamma)*s(g2) as initial guess
      ! and save at np1
      call elemstate_add(elem,statesave,nets,nete,1,np1,n0,n0,1d0,0d0,0d0)
      call compute_stage_value_dirk(np1,qn0,gamma*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
!      print *, 'num iters  ', maxiter
!=== End of Phase 2 ===
! at this point, un0+dt*(1-gamma)*(n(g2)+s(g2)) is at nm1, g3 is at np1, and n0 is free
       
     ! form unp1 = un0+dt*(1-gamma)*(n(g2)+s(g2))+dt*gamma*(n(g3)+s(g3))
      call compute_andor_apply_rhs(np1,nm1,np1,qn0,gamma*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,.false.,gamma*eta_ave_w,1.d0,1.d0,1.d0)     
            
      call state_read(elem,statesave,n0,nets,nete)
      call t_stopf("ARS232_timestep")
!======================================================================================================
   elseif (tstep_type==8) then ! kgs242
      call t_startf("KGS242_timestep")
      ! denote the stages as k1,...,k4 and note that k1 = un0
      a1 = 0.5
      a2 = 0.5
      a3 = 1.0
      dhat2 = 2.25
      dhat1 = (0.5-dhat1)/(1.0-dhat1)
      ahat2 = 0.5-dhat2
      ahat3 = 1.0

      ! compute un0 + dt*a1*n(k1=u(n0))+dt*ahat1*s(k1) and store at np1
      call compute_andor_apply_rhs(np1,n0,n0,qn0,a1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w*a1,1d0,0d0,1d0)

      maxiter=10
      itertol=1e-12
      ! solve k2 = u(n) + dt*a1*n(k1) + dt*dhat1*s(k2) and store solution at np1
      call compute_stage_value_dirk(np1,qn0,dhat1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
!================= end of phase 1 =========================================

     ! compute u(n)+dt*a2*n(k2) and store at np1
      call compute_andor_apply_rhs(np1,n0,np1,qn0,a2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,.false.,eta_ave_w/2,1d0,ahat2/a2,1d0)

      maxiter=10
      itertol=1e-12
      ! solve k3 = u(n) + dt*a2*n(k2) + dt*ahat2*s(k2) + dt*dhat2*s(k3) and store solution at np1
      call compute_stage_value_dirk(np1,qn0,dhat2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
! ================ end of phase 2 ========================================

     ! compute u(n+1) = k4 =  u(n)+dt*(a3*n(k3)+ahat3*s(k3)) and store at np1
      call compute_andor_apply_rhs(np1,n0,np1,qn0,a3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,.false.,eta_ave_w,1d0,ahat3/a3,1d0)

      call t_stopf("KGS242_timestep")
!===========================================================================================
    elseif (tstep_type==9) then ! kgs252
      call t_startf("KGS252_timestep")
      ! denote the stages as k1,...,k4 and note that k1 = un0
      a1 = 0.25
      a2 = 1d0/3d0
      a3 = 0.5
      a4 = 1.0
      dhat2 = 2.25
      dhat1 = (0.5-dhat1)/(1.0-dhat1)
      ahat3= 0.5-dhat2
      ahat4 = 1.0

     ! ============ first stage is pure explicit =======================

     ! compute k2 = u(n)+dt*a1*n(k1=u(n)) and store at np1
     call compute_andor_apply_rhs(np1,n0,n0,qn0,dt*a1,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w*a1,1d0,0d0,1d0)

     ! compute u(n)+dt*a2*n(k2) and store at np1
     call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a2,elem,hvcoord,hybrid,&
        deriv,nets,nete,.false.,eta_ave_w*a2,1d0,0d0,1d0)

      ! solve k3 = u(n)+dt*a2*n(k2)+dt*dhat1*s(k3) and store at np1
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
 !     print *, maxiter
 !     print *, itertol

     !  ============== end of stage 3 =============================
     ! compute u(n)+dt*a3*n(k3)+dt*ahat3*s(k3) and store at np1
     call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a3,elem,hvcoord,hybrid,&
     deriv,nets,nete,.false.,eta_ave_w*a3,1d0,ahat3/a3,1d0)

      ! solve k4 = u(n)+dt*a3*n(k3)+dt*ahat3*s(k3)+dt*dhat2*s(k4) and store at np1
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
 !     print *, maxiter
 !     print *, itertol

      ! ============= end of stage 4 =================================

     ! compute u(n+1) = k5 = u(n)+dt*a4*n(k4)+dt*ahat4*s(k4) and store at np1
     call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a4,elem,hvcoord,hybrid,&
     deriv,nets,nete,.false.,eta_ave_w*a4,1d0,ahat4/a4,1d0)

      call t_stopf("KGS252_timestep")
!===========================================================================================
    elseif (tstep_type==10)  then ! kgs 262
      call t_startf("KGS262_timestep")

     a1 = .25
     a2 = 1.0/6.0
     a3 = 3.0/8.0
     a4 = .5
     a5 = 1.0
     dhat2 = 2.25
     ahat4 = 0.5-dhat2
     dhat1 = (0.5-dhat2)/(1-dhat2)

    ! ======== first two stages are pure explicit  =============
     call compute_andor_apply_rhs(np1,n0,n0,qn0,dt*a1,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w*a1,1d0,0d0,1d0)

     call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a2,elem,hvcoord,hybrid,&
        deriv,nets,nete,.false.,eta_ave_w*a2,1d0,0d0,1d0)

    ! at this stage k2 is at np1, u(n) is at n0

    ! compute u(n)+dt*a3*n(k2) and store at np1
     call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a3,elem,hvcoord,hybrid,&
        deriv,nets,nete,.false.,eta_ave_w*a3,1d0,0d0,1d0)

      ! solve k3 = u(n)+dt*a3*n(k2)+dt*dhat1*s(k3) and store at np1
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
 !     print *, maxiter
 !     print *, itertol
    ! ========== end of stage 3 =================================

    ! compute u(n)+dt*a4*n(k3) and store at np1
     call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a4,elem,hvcoord,hybrid,&
        deriv,nets,nete,.false.,eta_ave_w*a4,1d0,ahat4/a4,1d0)

     ! solve k4 = u(n)+dt*a4*n(k2)+dt*ahat4*s(k3)+dt*dhat2*s(k4) and store at np1
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
!      print *, maxiter
!      print *, itertol

    ! ============ end of stage 4 ==================================

    ! compute u(n+1) = k5 = u(n)+dt*a5*n(k4), final stage is the solution update
     call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a5,elem,hvcoord,hybrid,&
       deriv,nets,nete,.false.,eta_ave_w*a5,1d0,1d0,1d0)

      call t_stopf("KGS262_timestep")

    else
       call abortmp('ERROR: bad choice of tstep_type')
    endif


    ! ==============================================
    ! Time-split Horizontal diffusion: nu.del^2 or nu.del^4
    ! U(*) = U(t+1)  + dt2 * HYPER_DIFF_TERM(t+1)
    ! ==============================================
    ! note:time step computes u(t+1)= u(t*) + RHS.
    ! for consistency, dt_vis = t-1 - t*, so this is timestep method dependent
    ! forward-in-time, hypervis applied to dp3d
    if (hypervis_order == 2 .and. nu>0) &
         call advance_hypervis(edge6,elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)


    ! warning: advance_physical_vis currently requires levels that are equally spaced in z
    if (dcmip16_mu>0) call advance_physical_vis(edge6,elem,hvcoord,hybrid,deriv,np1,nets,nete,dt,dcmip16_mu_s,dcmip16_mu)

    call t_stopf('prim_advance_exp')
  end subroutine prim_advance_exp



!placeholder
  subroutine applyCAMforcing_dp3d(elem,hvcoord,np1,dt,nets,nete)
  implicit none
  type (element_t),       intent(inout) :: elem(:)
  real (kind=real_kind),  intent(in)    :: dt
  type (hvcoord_t),       intent(in)    :: hvcoord
  integer,                intent(in)    :: np1,nets,nete
  end subroutine applyCAMforcing_dp3d

!temp solution for theta+ftype0
  subroutine applyCAMforcing_ps(elem,hvcoord,np1,np1_qdp,dt,nets,nete)

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
     call get_temperature(elem(ie),temperature,hvcoord,np1)
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        temperature(:,:,k) = temperature(:,:,k) + dt*elem(ie)%derived%FT(:,:,k)
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
     call get_kappa_star(kappa_star,elem(ie)%state%Q(:,:,:,1))
     call get_cp_star(cp_star,elem(ie)%state%Q(:,:,:,1))
     call get_pnh_and_exner(hvcoord,elem(ie)%state%theta_dp_cp(:,:,:,np1),dp,&
          elem(ie)%state%phinh(:,:,:,np1),elem(ie)%state%phis(:,:),kappa_star,&
          pnh,dpnh,exner)

     elem(ie)%state%theta_dp_cp(:,:,:,np1) = temperature(:,:,:)*cp_star(:,:,:)&
          *dp(:,:,:)/exner(:,:,:)

    enddo
    call applyCAMforcing_dynamics(elem,hvcoord,np1,np1_qdp,dt,nets,nete)
  end subroutine applyCAMforcing_ps





  subroutine applyCAMforcing_dynamics(elem,hvcoord,np1,np1_qdp,dt,nets,nete)

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


! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

  real (kind=real_kind), dimension(np,np,4) :: lap_s  ! dp3d,theta,w,phi
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: exner0(nlev)
  real (kind=real_kind) :: heating(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: p_i(np,np,nlevp)    ! pressure on interfaces
  real (kind=real_kind) :: dt
  real (kind=real_kind) :: ps_ref(np,np)

  real (kind=real_kind) :: theta_ref(np,np,nlev,nets:nete)
  real (kind=real_kind) :: phi_ref(np,np,nlev,nets:nete)
  real (kind=real_kind) :: dp_ref(np,np,nlev,nets:nete)

  call t_startf('advance_hypervis')


  dt=dt2/hypervis_subcycle
  
  do k=1,nlev
     exner0(k) = (hvcoord%etam(k)*hvcoord%ps0/p0 )**kappa
  enddo


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
     ps_ref(:,:) = hvcoord%hyai(1)*hvcoord%ps0 + sum(elem(ie)%state%dp3d(:,:,:,nt),3)
     do k=1,nlev
        dp_ref(:,:,k,ie) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
             (hvcoord%hybi(k+1)-hvcoord%hybi(k))*ps_ref(:,:)
     enddo

     call get_dry_phinh(hvcoord,elem(ie)%state%phis,&
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
  do ic=1,hypervis_subcycle
     do ie=nets,nete

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
        do k=1,nlev
           elem(ie)%state%theta_dp_cp(:,:,k,nt)=elem(ie)%state%theta_dp_cp(:,:,k,nt)-&
                theta_ref(:,:,k,ie)
           elem(ie)%state%phinh(:,:,k,nt)=elem(ie)%state%phinh(:,:,k,nt)-&
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
                (eta_ave_w*elem(ie)%state%dp3d(:,:,:,nt)+dp_ref(:,:,:,ie))/hypervis_subcycle
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
                 lap_s(:,:,1)=laplace_sphere_wk(elem(ie)%state%dp3d       (:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_s(:,:,2)=laplace_sphere_wk(elem(ie)%state%theta_dp_cp(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_s(:,:,3)=laplace_sphere_wk(elem(ie)%state%w          (:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_s(:,:,4)=laplace_sphere_wk(elem(ie)%state%phinh      (:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              endif

              nu_scale_top = 1
              if (k==1) nu_scale_top=4
              if (k==2) nu_scale_top=2

              ! biharmonic terms need a negative sign:
              if (nu_top>0 .and. k<=3) then
                 vtens(:,:,:,k,ie)=(  -nu*vtens(:,:,:,k,ie) + nu_scale_top*nu_top*lap_v(:,:,:)) ! u and v
                 stens(:,:,k,1,ie)=(-nu_p*stens(:,:,k,1,ie) + nu_scale_top*nu_top*lap_s(:,:,1)) ! dp3d
                 stens(:,:,k,2,ie)=(  -nu*stens(:,:,k,2,ie) + nu_scale_top*nu_top*lap_s(:,:,2)) ! theta
                 stens(:,:,k,3,ie)=(  -nu*stens(:,:,k,3,ie) + nu_scale_top*nu_top*lap_s(:,:,3)) ! w
                 stens(:,:,k,4,ie)=(-nu_s*stens(:,:,k,4,ie) + nu_scale_top*nu_top*lap_s(:,:,4)) ! phi
              else
                 vtens(:,:,:,k,ie)=-nu  *vtens(:,:,:,k,ie) ! u,v
                 stens(:,:,k,1,ie)=-nu_p*stens(:,:,k,1,ie) ! dp3d
                 stens(:,:,k,2,ie)=-nu  *stens(:,:,k,2,ie) ! theta
                 stens(:,:,k,3,ie)=-nu  *stens(:,:,k,3,ie) ! w
                 stens(:,:,k,4,ie)=-nu_s*stens(:,:,k,4,ie) ! phi
              endif

           enddo

           kptr=0;      call edgeVpack(edgebuf,vtens(:,:,:,:,ie),2*nlev,kptr,ie)
           kptr=2*nlev; call edgeVpack(edgebuf,stens(:,:,:,:,ie),4*nlev,kptr,ie)

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
              vtens(:,:,1,k,ie)=dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)  ! u
              vtens(:,:,2,k,ie)=dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)  ! v
              stens(:,:,k,1,ie)=dt*stens(:,:,k,1,ie)*elem(ie)%rspheremp(:,:)  ! dp3d
              stens(:,:,k,2,ie)=dt*stens(:,:,k,2,ie)*elem(ie)%rspheremp(:,:)  ! theta
              stens(:,:,k,3,ie)=dt*stens(:,:,k,3,ie)*elem(ie)%rspheremp(:,:)  ! w
              stens(:,:,k,4,ie)=dt*stens(:,:,k,4,ie)*elem(ie)%rspheremp(:,:)  ! phi

              !add ref state back
              elem(ie)%state%theta_dp_cp(:,:,k,nt)=elem(ie)%state%theta_dp_cp(:,:,k,nt)+&
                   theta_ref(:,:,k,ie)
              elem(ie)%state%phinh(:,:,k,nt)=elem(ie)%state%phinh(:,:,k,nt)+&
                   phi_ref(:,:,k,ie)
              elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt)+&
                   dp_ref(:,:,k,ie)

           enddo



#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
        do k=1,nlev
           elem(ie)%state%v(:,:,:,k,nt)=elem(ie)%state%v(:,:,:,k,nt) + &
                vtens(:,:,:,k,ie)
           elem(ie)%state%w(:,:,k,nt)=elem(ie)%state%w(:,:,k,nt) &
                +stens(:,:,k,3,ie)
           
           elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt) &
                +stens(:,:,k,1,ie)
           
           elem(ie)%state%w(:,:,k,nt)=elem(ie)%state%w(:,:,k,nt) &
                +stens(:,:,k,3,ie)
           
           elem(ie)%state%phinh(:,:,k,nt)=elem(ie)%state%phinh(:,:,k,nt) &
                +stens(:,:,k,4,ie)
        enddo


        ! apply heating after updating sate.  using updated v gives better results in PREQX model
        !
        ! d(IE)/dt =  exner * d(Theta)/dt + phi d(dp3d)/dt   (Theta = dp3d*cp*theta)
        !   Our eqation:  d(theta)/dt = diss(theta) + heating
        !   Assuming no diffusion on dp3d, we can approximate by:
        !   d(IE)/dt = exner*cp*dp3d * diss(theta)  -   exner*cp*dp3d*heating               
        !
        ! KE dissipaiton will be given by:
        !   d(KE)/dt = dp3d*U dot diss(U)
        ! we want exner*cp*dp3d*heating = dp3d*U dot diss(U)
        ! and thus heating =  U dot diss(U) / exner*cp
        ! 
        ! PE dissipation
        ! d(PE)/dt = dp3d diss(phi) 
        !     we want dp3d diss(phi) = exner*cp*dp3d*heating
        !     heating = diss(phi) / exner*cp
        !
        ! use hydrostatic pressure for simplicity
        p_i(:,:,1)=hvcoord%hyai(1)*hvcoord%ps0
        do k=1,nlev
           p_i(:,:,k+1)=p_i(:,:,k) + elem(ie)%state%dp3d(:,:,k,nt)
        enddo
#if (defined COLUMN_OPENMP)
!$omp parallel do default(shared), private(k)
#endif
        do k=1,nlev
           
           ! p(:,:,k) = (p_i(:,:,k) + elem(ie)%state%dp3d(:,:,k,nt)/2)
           exner(:,:,k)  = ( (p_i(:,:,k) + elem(ie)%state%dp3d(:,:,k,nt)/2) /p0)**kappa
           if (theta_hydrostatic_mode) then
              heating(:,:,k)= (elem(ie)%state%v(:,:,1,k,nt)*vtens(:,:,1,k,ie) + &
                   elem(ie)%state%v(:,:,2,k,nt)*vtens(:,:,2,k,ie) ) / &
                   (exner(:,:,k)*Cp)

           else
              heating(:,:,k)= (elem(ie)%state%v(:,:,1,k,nt)*vtens(:,:,1,k,ie) + &
                   elem(ie)%state%v(:,:,2,k,nt)*vtens(:,:,2,k,ie)  +&
                   elem(ie)%state%w(:,:,k,nt)*stens(:,:,k,3,ie)  +&
                   stens(:,:,k,4,ie) ) / &
                   (exner(:,:,k)*Cp)  
           endif
           
           elem(ie)%state%theta_dp_cp(:,:,k,nt)=elem(ie)%state%theta_dp_cp(:,:,k,nt) &
                +stens(:,:,k,2,ie)*hvcoord%dp0(k)*exner0(k)/(exner(:,:,k)*elem(ie)%state%dp3d(:,:,k,nt))&
                -heating(:,:,k)
        enddo
        
     enddo
  enddo

! convert theta_dp_cp -> theta
  do ie=nets,nete            
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        elem(ie)%state%theta_dp_cp(:,:,k,nt)=&
             elem(ie)%state%theta_dp_cp(:,:,k,nt)*Cp*elem(ie)%state%dp3d(:,:,k,nt)
     enddo
  enddo


  call t_stopf('advance_hypervis')

  end subroutine advance_hypervis





  subroutine advance_physical_vis(edgebuf,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt,mu_s,mu)
  !
  !  take one timestep of of physical viscosity (single laplace operator) for
  !  all state variables in both horizontal and vertical
  !  
  !  as of 2017/5, used only for the supercell test case
  !  so for now:
  !     dont bother to optimize
  !     apply only to perturbation from background state (supercell initial condition)
  !     uniform spacing in z with delz = 20km/nlev
  !
  !

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edgebuf
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord

  real (kind=real_kind) :: dt, mu_s, mu
  integer :: nt,nets,nete

  ! local
  integer :: k,kptr,ie
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,4,nets:nete)      :: stens  ! dp3d,theta,w,phi


  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: delz

  real (kind=real_kind) :: theta_ref(np,np,nlev)

  real (kind=real_kind) :: theta_prime(np,np,nlev)
  real (kind=real_kind) :: phi_prime(np,np,nlev)
  real (kind=real_kind) :: dp_prime(np,np,nlev)
  real (kind=real_kind) :: w_prime(np,np,nlev)
  real (kind=real_kind) :: u_prime(np,np,2,nlev)

  !if(test_case .ne. 'dcmip2016_test3') call abortmp("dcmip16_mu is currently limited to dcmip16 test 3")

  call t_startf('advance_physical_vis')
  delz = 20d3/nlev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute reference states
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie=nets,nete

     theta_ref(:,:,:) = state0(ie)%theta_dp_cp(:,:,:,1)/(cp*state0(ie)%dp3d(:,:,:,1))
     elem(ie)%state%theta_dp_cp(:,:,:,nt)=&
             elem(ie)%state%theta_dp_cp(:,:,:,nt)/(Cp*elem(ie)%state%dp3d(:,:,:,nt))

     ! perturbation variables
     u_prime(:,:,:,:)  = elem(ie)%state%v(:,:,:,:,nt)        -state0(ie)%v(:,:,:,:,1)
     w_prime(:,:,:)    = elem(ie)%state%w(:,:,:,nt)          -state0(ie)%w(:,:,:,1)
     dp_prime(:,:,:)   = elem(ie)%state%dp3d(:,:,:,nt)       -state0(ie)%dp3d(:,:,:,1)
     phi_prime(:,:,:)  = elem(ie)%state%phinh(:,:,:,nt)      -state0(ie)%phinh(:,:,:,1)
     theta_prime(:,:,:)= elem(ie)%state%theta_dp_cp(:,:,:,nt)-theta_ref(:,:,:)

     ! vertical viscosity
     call laplace_z(u_prime,    vtens(:,:,:,:,ie),2,nlev,delz)
     call laplace_z(dp_prime,   stens(:,:,:,1,ie),1,nlev,delz)
     call laplace_z(theta_prime,stens(:,:,:,2,ie),1,nlev,delz)
     call laplace_z(w_prime,    stens(:,:,:,3,ie),1,nlev,delz)
     call laplace_z(phi_prime,  stens(:,:,:,4,ie),1,nlev,delz)

     ! add in horizontal viscosity
     ! multiply by mass matrix for DSS
     ! horiz viscosity already has mass matrix built in
     do k=1,nlev
        lap_v = vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)

        vtens(:,:,1,k,ie) = (vtens(:,:,1,k,ie)*elem(ie)%spheremp(:,:) + lap_v(:,:,1))
        vtens(:,:,2,k,ie) = (vtens(:,:,2,k,ie)*elem(ie)%spheremp(:,:) + lap_v(:,:,2))

        stens(:,:,k,1,ie) = (stens(:,:,k,1,ie)*elem(ie)%spheremp(:,:) + &
             laplace_sphere_wk(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),var_coef=.false.)  )

        stens(:,:,k,2,ie) = (stens(:,:,k,2,ie)*elem(ie)%spheremp(:,:) + &
             laplace_sphere_wk(elem(ie)%state%theta_dp_cp(:,:,k,nt),deriv,elem(ie),var_coef=.false.)  )

        stens(:,:,k,3,ie) = (stens(:,:,k,3,ie)*elem(ie)%spheremp(:,:) + &
             laplace_sphere_wk(elem(ie)%state%w(:,:,k,nt),deriv,elem(ie),var_coef=.false.) )

        stens(:,:,k,4,ie) = (stens(:,:,k,4,ie)*elem(ie)%spheremp(:,:) + &
             laplace_sphere_wk(elem(ie)%state%phinh(:,:,k,nt),deriv,elem(ie),var_coef=.false.) ) 

     enddo

     kptr=0
     call edgeVpack(edgebuf,vtens(:,:,:,:,ie),2*nlev,kptr,ie)
     kptr=2*nlev
     call edgeVpack(edgebuf,stens(:,:,:,:,ie),4*nlev,kptr,ie)
     
  enddo

  call bndry_exchangeV(hybrid,edgebuf)
  
  do ie=nets,nete
     
     kptr=0
     call edgeVunpack(edgebuf, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
     kptr=2*nlev
     call edgeVunpack(edgebuf, stens(:,:,:,:,ie), 4*nlev, kptr, ie)
     
     ! apply inverse mass matrix, accumulate tendencies
     do k=1,nlev
        elem(ie)%state%v(:,:,1,k,nt)=elem(ie)%state%v(:,:,1,k,nt) + &
             mu*dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)

        elem(ie)%state%v(:,:,2,k,nt)=elem(ie)%state%v(:,:,2,k,nt) + &
             mu*dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
        
        elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt) &
             +mu_s*dt*stens(:,:,k,1,ie)*elem(ie)%rspheremp(:,:)
        
        elem(ie)%state%theta_dp_cp(:,:,k,nt)=elem(ie)%state%theta_dp_cp(:,:,k,nt) &
             +mu_s*dt*stens(:,:,k,2,ie)*elem(ie)%rspheremp(:,:)

        elem(ie)%state%w(:,:,k,nt)=elem(ie)%state%w(:,:,k,nt) &
             +mu_s*dt*stens(:,:,k,3,ie)*elem(ie)%rspheremp(:,:)
        
        elem(ie)%state%phinh(:,:,k,nt)=elem(ie)%state%phinh(:,:,k,nt) &
             +mu_s*dt*stens(:,:,k,4,ie)*elem(ie)%rspheremp(:,:)
        
     enddo
  enddo


  ! convert theta_dp_cp -> theta
  do ie=nets,nete            
     do k=1,nlev
        elem(ie)%state%theta_dp_cp(:,:,k,nt)=&
             elem(ie)%state%theta_dp_cp(:,:,k,nt)*Cp*elem(ie)%state%dp3d(:,:,k,nt)
     enddo
  enddo


  call t_stopf('advance_physical_vis')

  end subroutine advance_physical_vis





!============================ stiff and or non-stiff ============================================

 subroutine compute_andor_apply_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w,scale1,scale2,scale3)
  ! ===================================
  ! compute the RHS, accumulate into u(np1) and apply DSS
  !
  !   u(np1) = scale3*u(nm1) + dt2*DSS[ nonstiffRHS(u(n0))*scale1 + stiffRHS(un0)*scale2 ]
  !
  ! This subroutine was orgininally called to compute a leapfrog timestep
  ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
  ! accomodated.  For example, setting nm1=np1=n0 this routine will
  ! take a forward euler step, overwriting the input with the output.
  !
  !    qn0 = timelevel used to access Qdp() in order to compute virtual Temperature
  !
  ! ===================================

  integer,              intent(in) :: np1,nm1,n0,qn0,nets,nete
  real*8,               intent(in) :: dt2
  logical,              intent(in) :: compute_diagnostics
  type (hvcoord_t),     intent(in) :: hvcoord
  type (hybrid_t),      intent(in) :: hybrid
  type (element_t),     intent(inout), target :: elem(:)
  type (derivative_t),  intent(in) :: deriv

  real (kind=real_kind) :: eta_ave_w,scale1,scale2,scale3  ! weighting for eta_dot_dpdn mean flux, scale of unm1

  ! local
  real (kind=real_kind), pointer, dimension(:,:,:) :: phi
  real (kind=real_kind), pointer, dimension(:,:,:) :: dp3d
  real (kind=real_kind), pointer, dimension(:,:,:) :: theta_dp_cp
   
  real (kind=real_kind) :: kappa_star(np,np,nlev)
  real (kind=real_kind) :: theta_cp(np,np,nlev)
  real (kind=real_kind) :: theta_bar(np,np,nlevp)
  real (kind=real_kind) :: omega_p(np,np,nlev)
  real (kind=real_kind) :: vort(np,np,nlev)           ! vorticity
  real (kind=real_kind) :: divdp(np,np,nlev)
  real (kind=real_kind) :: pi(np,np,nlev)             ! hydrostatic pressure
  real (kind=real_kind) :: pi_i(np,np,nlevp)          ! hydrostatic pressure, interfaces
  real (kind=real_kind) :: pnh(np,np,nlev)            ! nh (nonydro) pressure
  real (kind=real_kind) :: dpnh(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)          ! exner nh pressure
  real (kind=real_kind) :: dpnh_dp(np,np,nlev)        ! dpnh / dp3d
  real (kind=real_kind) :: eta_dot_dpdn(np,np,nlevp)  ! vertical velocity at interfaces
  real (kind=real_kind) :: KE(np,np,nlev)             ! Kinetic energy
  real (kind=real_kind) :: gradexner(np,np,2,nlev)    ! grad(p^kappa)
  real (kind=real_kind) :: gradphi(np,np,2,nlev)     
  real (kind=real_kind) :: gradKE(np,np,2,nlev)       ! grad(0.5 u^T u )

  real (kind=real_kind) :: grad_kappastar(np,np,2,nlev)

  real (kind=real_kind) :: v_gradw(np,np,nlev)     
  real (kind=real_kind) :: v_gradtheta(np,np,nlev)     
  real (kind=real_kind) :: v_theta(np,np,2,nlev)
  real (kind=real_kind) :: div_v_theta(np,np,nlev)
  real (kind=real_kind) :: v_gradphi(np,np,nlev)
  real (kind=real_kind) :: v_gradKE(np,np,nlev)     
  real (kind=real_kind) :: wvor(np,np,2,nlev)

  real (kind=real_kind) :: vtens1(np,np,nlev)
  real (kind=real_kind) :: vtens2(np,np,nlev)
  real (kind=real_kind) :: v_vadv(np,np,2,nlev)       ! velocity vertical advection
  real (kind=real_kind) :: s_state(np,np,nlev,2)      ! scalars w,phi
  real (kind=real_kind) :: s_vadv(np,np,nlev,3)       ! w,phi, theta  vertical advection term
  real (kind=real_kind) :: stens(np,np,nlev,3)        ! tendencies w,phi,theta

  real (kind=real_kind) ::  temp(np,np,nlev)
  real (kind=real_kind), dimension(np,np) :: sdot_sum ! temporary field
  real (kind=real_kind), dimension(np,np,2) :: vtemp  ! generic gradient storage
  real (kind=real_kind) ::  v1,v2,w,d_eta_dot_dpdn_dn
  integer :: i,j,k,kptr,ie

  real (kind=real_kind), dimension(np,np,nlev) :: p, vgrad_p


  call t_startf('compute_andor_apply_rhs')

  do ie=nets,nete
     dp3d  => elem(ie)%state%dp3d(:,:,:,n0)
     theta_dp_cp  => elem(ie)%state%theta_dp_cp(:,:,:,n0)
     theta_cp(:,:,:) = theta_dp_cp(:,:,:)/dp3d(:,:,:)
     phi => elem(ie)%state%phinh(:,:,:,n0)

     call get_kappa_star(kappa_star,elem(ie)%state%Q(:,:,:,1))

     call get_pnh_and_exner(hvcoord,theta_dp_cp,dp3d,phi,elem(ie)%state%phis,&
             kappa_star,pnh,dpnh,exner)

     if (theta_hydrostatic_mode) then
        ! traditional Hydrostatic integral
        do k=1,nlev
           !temp(:,:,k) = theta_dp_cp(:,:,k)*(exner_i(:,:,k+1)-exner_i(:,:,k))/dp3d(:,:,k)           
           temp(:,:,k) = kappa_star(:,:,k)*theta_dp_cp(:,:,k)*exner(:,:,k)/pnh(:,:,k)
        enddo
        !call preq_hydrostatic(phi,elem(ie)%state%phis,T_v,p,dp)
        call preq_hydrostatic_v2(phi,elem(ie)%state%phis,temp)
        dpnh_dp(:,:,:) = 1
     else
        ! d(p-nh) / d(p-hyrdostatic)
        dpnh_dp(:,:,:) = dpnh(:,:,:)/dp3d(:,:,:)
     endif

     ! Compute omega_p = 1/pi Dpi/Dt
     ! first compute hydrostatic pressure
     pi_i(:,:,1)=hvcoord%hyai(1)*hvcoord%ps0
     do k=1,nlev
        pi_i(:,:,k+1)=pi_i(:,:,k) + dp3d(:,:,k)
     enddo
     do k=1,nlev
        pi(:,:,k)=pi_i(:,:,k) + dp3d(:,:,k)/2
     enddo
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,vtemp)
#endif
     do k=1,nlev
        ! get hydrostatic pressure flux to compute omega=dp/dt
        vtemp(:,:,:) = gradient_sphere( pi(:,:,k), deriv, elem(ie)%Dinv);
        vgrad_p(:,:,k) = elem(ie)%state%v(:,:,1,k,n0)*vtemp(:,:,1)+&
             elem(ie)%state%v(:,:,2,k,n0)*vtemp(:,:,2)
     enddo        
     call preq_omega_ps(omega_p,hvcoord,pi,vgrad_p,divdp)


#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,vtemp)
#endif
     do k=1,nlev
        vtemp(:,:,1) = elem(ie)%state%v(:,:,1,k,n0)*dp3d(:,:,k)
        vtemp(:,:,2) = elem(ie)%state%v(:,:,2,k,n0)*dp3d(:,:,k)

        ! ================================
        ! Accumulate mean Vel_rho flux in vn0
        ! ================================
        elem(ie)%derived%vn0(:,:,:,k)=elem(ie)%derived%vn0(:,:,:,k)+eta_ave_w*vtemp(:,:,:)

        divdp(:,:,k)=divergence_sphere(vtemp(:,:,:),deriv,elem(ie))
        vort(:,:,k)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,n0),deriv,elem(ie))
     enddo




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
        theta_bar(:,:,1     ) = 0.0D0
        theta_bar(:,:,nlev+1) = 0.0D0

        ! ===========================================================
        ! Compute vertical advection of T and v from eq. CCM2 (3.b.1)
        ! ==============================================
        ! TODO: remove theta from s_state and s_vadv
        s_state(:,:,:,1)=elem(ie)%state%w(:,:,:,n0)
        s_state(:,:,:,2)=phi
        call preq_vertadv_v(elem(ie)%state%v(:,:,:,:,n0),s_state,2,eta_dot_dpdn,dp3d,v_vadv,s_vadv)
        !call preq_vertadv_v(elem(ie)%state%v(:,:,:,:,n0),s_state,3,eta_dot_dpdn,dp3d,v_vadv,s_vadv)
        !call preq_vertadv_upwind(elem(ie)%state%v(:,:,:,:,n0),s_state,3,eta_dot_dpdn,dp3d,v_vadv,s_vadv)

        ! this loop constructs d( eta-dot * theta_dp_cp)/deta
        ! d( eta_dot_dpdn * theta*cp)
        ! so we need to compute theta_cp form theta_dp_cp and average to interfaces
#if 1
        do k=2,nlev  ! E conserving averaging, but much more unstable
           theta_bar(:,:,k) = - ( (dpnh_dp(:,:,k) + dpnh_dp(:,:,k-1))/2 ) * &
                (phi(:,:,k)-phi(:,:,k-1))/(exner(:,:,k)-exner(:,:,k-1))
        enddo
#else
        do k=2,nlev  ! simple averaging
           theta_bar(:,:,k) = (theta_cp(:,:,k)+theta_cp(:,:,k-1))/2
        enddo
#endif

        do k=1,nlev
           s_vadv(:,:,k,3)= &
              eta_dot_dpdn(:,:,k+1)* theta_bar(:,:,k+1)  - &
              eta_dot_dpdn(:,:,k)  * theta_bar(:,:,k)
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
!$omp parallel do private(k,i,j,v1,v2,vtemp)
#endif
     vertloop: do k=1,nlev
        ! w-vorticity correction term added to u momentum equation for E conservation
        vtemp(:,:,:)  = gradient_sphere(elem(ie)%state%w(:,:,k,n0),deriv,elem(ie)%Dinv)

        temp(:,:,k) = (elem(ie)%state%w(:,:,k,n0)**2)/2  
        wvor(:,:,:,k) = gradient_sphere(temp(:,:,k),deriv,elem(ie)%Dinv)
        wvor(:,:,1,k) = wvor(:,:,1,k) - elem(ie)%state%w(:,:,k,n0)*vtemp(:,:,1)
        wvor(:,:,2,k) = wvor(:,:,2,k) - elem(ie)%state%w(:,:,k,n0)*vtemp(:,:,2)

        v_gradw(:,:,k) = elem(ie)%state%v(:,:,1,k,n0)*vtemp(:,:,1) &
             +elem(ie)%state%v(:,:,2,k,n0)*vtemp(:,:,2) 

        ! ================================================
        ! w,theta,phi tendencies:
        ! ================================================
        stens(:,:,k,1) = (-s_vadv(:,:,k,1) - v_gradw(:,:,k))*scale1 - scale2*g*(1-dpnh_dp(:,:,k) )
        v_theta(:,:,1,k) = elem(ie)%state%v(:,:,1,k,n0)*               &
          elem(ie)%state%theta_dp_cp(:,:,k,n0)
        v_theta(:,:,2,k) =                                             &
          elem(ie)%state%v(:,:,2,k,n0)                                 &
          *elem(ie)%state%theta_dp_cp(:,:,k,n0)
        div_v_theta(:,:,k)=divergence_sphere(v_theta(:,:,:,k),deriv,elem(ie))
        stens(:,:,k,3)=(-s_vadv(:,:,k,3)-div_v_theta(:,:,k))*scale1

        gradphi(:,:,:,k) = gradient_sphere(phi(:,:,k),deriv,elem(ie)%Dinv)

        v_gradphi(:,:,k) = elem(ie)%state%v(:,:,1,k,n0)*gradphi(:,:,1,k) &
             +elem(ie)%state%v(:,:,2,k,n0)*gradphi(:,:,2,k) 
        ! use of s_vadv(:,:,k,2) here is correct since this corresponds to etadot d(phi)/deta
        stens(:,:,k,2) =  (-s_vadv(:,:,k,2) - v_gradphi(:,:,k))*scale1 + scale2*g*elem(ie)%state%w(:,:,k,n0)

        KE(:,:,k) = ( elem(ie)%state%v(:,:,1,k,n0)**2 + elem(ie)%state%v(:,:,2,k,n0)**2)/2
        gradKE(:,:,:,k) = gradient_sphere(KE(:,:,k),deriv,elem(ie)%Dinv)
        gradexner(:,:,:,k) = gradient_sphere(exner(:,:,k),deriv,elem(ie)%Dinv)

        grad_kappastar(:,:,:,k) = gradient_sphere(kappa_star(:,:,k),deriv,elem(ie)%Dinv)

        
        

        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              v_gradKE(i,j,k)=v1*gradKE(i,j,1,k)+v2*gradKE(i,j,2,k)

              vtens1(i,j,k) = (-v_vadv(i,j,1,k) &
                   + v2*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - gradKE(i,j,1,k) - gradphi(i,j,1,k)*dpnh_dp(i,j,k) &
                  -theta_cp(i,j,k)*gradexner(i,j,1,k)&
                  +theta_cp(i,j,k)*grad_kappastar(i,j,1,k)*exner(i,j,k)*log(pnh(i,j,k)/p0)&
                  -wvor(i,j,1,k) &
                  )*scale1


              vtens2(i,j,k) = (-v_vadv(i,j,2,k) &
                   - v1*(elem(ie)%fcor(i,j) + vort(i,j,k)) &
                   - gradKE(i,j,2,k) - gradphi(i,j,2,k)*dpnh_dp(i,j,k) &
                  -theta_cp(i,j,k)*gradexner(i,j,2,k) &
                  +theta_cp(i,j,k)*grad_kappastar(i,j,2,k)*exner(i,j,k)*log(pnh(i,j,k)/p0) &
                  -wvor(i,j,2,k) &
                  )*scale1

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
        elem(ie)%accum%KEu_horiz1=0
        elem(ie)%accum%KEu_horiz2=0
        elem(ie)%accum%KEu_vert1=0
        elem(ie)%accum%KEu_vert2=0
        elem(ie)%accum%KEw_horiz1=0
        elem(ie)%accum%KEw_horiz2=0
        elem(ie)%accum%KEw_horiz3=0
        elem(ie)%accum%KEw_vert1=0
        elem(ie)%accum%KEw_vert2=0

        elem(ie)%accum%PEhoriz1=0
        elem(ie)%accum%PEhoriz2=0
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
!$omp parallel do private(k,i,j,d_eta_dot_dpdn_dn)
#endif
        do k =1,nlev
          do j=1,np
            do i=1,np                
                  d_eta_dot_dpdn_dn=(eta_dot_dpdn(i,j,k+1)-eta_dot_dpdn(i,j,k))
               !  Form horiz advection of KE-u
                  elem(ie)%accum%KEu_horiz1(i,j)=elem(ie)%accum%KEu_horiz1(i,j)              &
                  -v_gradKE(i,j,k)*dp3d(i,j,k) 
                  elem(ie)%accum%KEu_horiz2(i,j)=elem(ie)%accum%KEu_horiz2(i,j)              &
                  -KE(i,j,k)*divdp(i,j,k)
               !  Form horiz advection of KE-w
                  elem(ie)%accum%KEw_horiz1(i,j)=elem(ie)%accum%KEw_horiz1(i,j)   &
                  -dp3d(i,j,k) * elem(ie)%state%w(i,j,k,n0) * v_gradw(i,j,k)      
                  elem(ie)%accum%KEw_horiz2(i,j)=elem(ie)%accum%KEw_horiz2(i,j)-   &
                       0.5*(elem(ie)%state%w(i,j,k,n0))**2 * divdp(i,j,k)
                  elem(ie)%accum%KEw_horiz3(i,j)=elem(ie)%accum%KEw_horiz3(i,j)   &
                  -dp3d(i,j,k) * (elem(ie)%state%v(i,j,1,k,n0) * wvor(i,j,1,k) +  &
                                  elem(ie)%state%v(i,j,2,k,n0) * wvor(i,j,2,k))
               !  Form vertical advection of KE-u 
                  elem(ie)%accum%KEu_vert1(i,j)=elem(ie)%accum%KEu_vert1(i,j)- &
                  (elem(ie)%state%v(i,j,1,k,n0) * v_vadv(i,j,1,k) +            &
                  elem(ie)%state%v(i,j,2,k,n0) *v_vadv(i,j,2,k))*dp3d(i,j,k)
                  elem(ie)%accum%KEu_vert2(i,j)=elem(ie)%accum%KEu_vert2(i,j)- &
                  0.5*((elem(ie)%state%v(i,j,1,k,n0))**2 +                     &
                       (elem(ie)%state%v(i,j,2,k,n0))**2)*d_eta_dot_dpdn_dn
               !  Form vertical advection of KE-w
                  elem(ie)%accum%KEw_vert1(i,j)=elem(ie)%accum%KEw_vert1(i,j)      &
                  -s_vadv(i,j,k,1)*elem(ie)%state%w(i,j,k,n0)*dp3d(i,j,k)    
                  elem(ie)%accum%KEw_vert2(i,j)=elem(ie)%accum%KEw_vert2(i,j)      &
                  -0.5*d_eta_dot_dpdn_dn*(elem(ie)%state%w(i,j,k,n0)**2)

               !  Form IEvert1
                  elem(ie)%accum%IEvert1(i,j)=elem(ie)%accum%IEvert1(i,j)      &
                  -exner(i,j,k)*s_vadv(i,j,k,3)                        
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
                  -dp3d(i,j,k)*s_vadv(i,j,k,2)
               !  Form T01
                  elem(ie)%accum%T01(i,j)=elem(ie)%accum%T01(i,j)               &
                  -(elem(ie)%state%theta_dp_cp(i,j,k,n0))                       &
                  *(gradexner(i,j,1,k)*elem(ie)%state%v(i,j,1,k,n0) +           &
                  gradexner(i,j,2,k)*elem(ie)%state%v(i,j,2,k,n0))              
               !  Form S1 
                  elem(ie)%accum%S1(i,j)=elem(ie)%accum%S1(i,j)                 &
                  -exner(i,j,k)*div_v_theta(i,j,k)
               !  Form T2  = -S2 (no reason to compute S2?)
                  elem(ie)%accum%T2(i,j)=elem(ie)%accum%T2(i,j)+                & 
                  (g*elem(ie)%state%w(i,j,k,n0)-v_gradphi(i,j,k))               &
                  *dpnh(i,j,k)                                 
               !  Form S2
                  elem(ie)%accum%S2(i,j)=elem(ie)%accum%S2(i,j)                 &
                  +(v_gradphi(i,j,k)-g*elem(ie)%state%w(i,j,k,n0))              &
                  *dpnh(i,j,k)
               !  Form P1  = -P2  (no reason to compute P2?)
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
          + dt2*stens(:,:,k,3))
        elem(ie)%state%phinh(:,:,k,np1)   = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%phinh(:,:,k,nm1) & 
          + dt2*stens(:,:,k,2))

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
     call edgeVpack(edge6, elem(ie)%state%phinh(:,:,:,np1),nlev,kptr,ie)
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
     call edgeVunpack(edge6, elem(ie)%state%phinh(:,:,:,np1), nlev, kptr, ie)
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
        elem(ie)%state%phinh(:,:,k,np1)=elem(ie)%rspheremp(:,:)*elem(ie)%state%phinh(:,:,k,np1)
        elem(ie)%state%v(:,:,1,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,1,k,np1)
        elem(ie)%state%v(:,:,2,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,2,k,np1)
     end do
  end do

  call t_stopf('compute_andor_apply_rhs')

  end subroutine compute_andor_apply_rhs




 
!===========================================================================================================
!===========================================================================================================
!===========================================================================================================
!===========================================================================================================
!===========================================================================================================

  subroutine compute_stage_value_dirk(np1,qn0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,maxiter,itertol)
  !===================================================================================
  ! this subroutine solves a stage value equation for a DIRK method which takes the form
  !
  ! gi = un0 + dt* sum(1:i-1)(aij n(gj)+a2ij s(gj)) + dt *a2ii s(gi) := y + dt a2ii s(gi)
  !
  ! It is assumed that un0 has the value of y and the computed value of gi is stored at
  ! unp1
  !===================================================================================
  integer, intent(in) :: np1,qn0,nets,nete
  real*8, intent(in) :: dt2
  integer :: maxiter
  real*8 :: itertol

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv

  ! local
  real (kind=real_kind), pointer, dimension(:,:,:)   :: phi_np1
  real (kind=real_kind), pointer, dimension(:,:,:)   :: dp3d
  real (kind=real_kind), pointer, dimension(:,:,:)   :: theta_dp_cp
  real (kind=real_kind), pointer, dimension(:,:)   :: phis
  real (kind=real_kind) :: JacD(nlev,np,np)  , JacL(nlev-1,np,np)
  real (kind=real_kind) :: JacU(nlev-1,np,np), JacU2(nlev-2,np,np)
  real (kind=real_kind) :: kappa_star(np,np,nlev),kappa_star_i(np,np,nlevp)
  real (kind=real_kind) :: pnh(np,np,nlev)     ! nh (nonydro) pressure
  real (kind=real_kind) :: dpnh(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)     ! exner nh pressure
  real (kind=real_kind) :: dpnh_dp(np,np,nlev)    !    ! dpnh / dp3d
  real (kind=real_kind) :: Ipiv(nlev,np,np)
  real (kind=real_kind) :: Fn(np,np,nlev),x(nlev,np,np)
  real (kind=real_kind) :: pnh_i(np,np,nlevp)
  real (kind=real_kind) :: itererr,itererrtemp(np,np)
  real (kind=real_kind) :: itererrmat,itercountmax,itererrmax
  real (kind=real_kind) :: norminfr0(np,np),norminfJ0(np,np)
  real (kind=real_kind) :: maxnorminfJ0r0
  real (kind=real_kind) :: alpha1(np,np),alpha2(np,np)

  real (kind=real_kind) :: Jac2D(nlev,np,np)  , Jac2L(nlev-1,np,np)
  real (kind=real_kind) :: Jac2U(nlev-1,np,np)



  integer :: i,j,k,l,ie,itercount,info(np,np)
  itercountmax=0
  itererrmax=0.d0

  call t_startf('compute_stage_value_dirk')
  do ie=nets,nete
  
    itercount=0

    ! approximate the initial error of f(x) \approx 0
    dp3d  => elem(ie)%state%dp3d(:,:,:,np1)
    theta_dp_cp  => elem(ie)%state%theta_dp_cp(:,:,:,np1)
    phi_np1 => elem(ie)%state%phinh(:,:,:,np1)
    phis => elem(ie)%state%phis(:,:)
    call get_kappa_star(kappa_star,elem(ie)%state%Q(:,:,:,1))
    if (theta_hydrostatic_mode) then
      dpnh_dp(:,:,:)=1.d0
    else
      call get_pnh_and_exner(hvcoord,theta_dp_cp,dp3d,phi_np1,phis,&
        kappa_star,pnh,dpnh,exner,pnh_i_out=pnh_i)
        dpnh_dp(:,:,:) = dpnh(:,:,:)/dp3d(:,:,:)
    end if
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
    do k=1,nlev-1
      kappa_star_i(:,:,k+1) = 0.5D0* (kappa_star(:,:,k+1)+kappa_star(:,:,k))
    end do
    kappa_star_i(:,:,1) = kappa_star(:,:,1)
    kappa_star_i(:,:,nlev+1) = kappa_star(:,:,nlev)

   ! we first compute the initial Jacobian J0 and residual r0 and their infinity norms
     Fn(:,:,1:nlev) = phi_np1-elem(ie)%state%phinh(:,:,:,np1) &
       - dt2*g*elem(ie)%state%w(:,:,:,np1) + (dt2*g)**2 * (1.0-dpnh_dp(:,:,:))

     elem(ie)%state%w(:,:,1:nlev,np1) = elem(ie)%state%w(:,:,1:nlev,np1) - g*dt2 * &
        (1.0-dpnh_dp(:,:,1:nlev))


     norminfr0=0.d0
     norminfJ0=0.d0
     call get_dirk_jacobian(JacL,JacD,JacU,dt2,dp3d,phi_np1,phis,kappa_star_i,pnh_i,1)
    ! compute dp3d-weighted infinity norms of the initial Jacobian and residual
#if (defined COLUMN_OPENMP)
!$omp parallel do private(i,j) collapse(2)
#endif
     do i=1,np
     do j=1,np
       itererrtemp(i,j)=0 
       do k=1,nlev
        norminfr0(i,j)=max(norminfr0(i,j),abs(Fn(i,j,k)) *dp3d(i,j,k))
        if (k.eq.1) then
          norminfJ0(i,j) = max(norminfJ0(i,j),(abs(JacD(k,i,j))+abs(JacU(k,i,j)))*dp3d(i,j,k))
        elseif (k.eq.nlev) then
          norminfJ0(i,j) = max(norminfJ0(i,j),(abs(JacL(k,i,j))+abs(JacD(k,i,j)))*dp3d(i,j,k))
        else
          norminfJ0(i,j) = max(norminfJ0(i,j),(abs(JacL(k,i,j))+abs(JacD(k,i,j))+ &
            abs(JacU(k,i,j)))*dp3d(i,j,k))
        end if
        itererrtemp(i,j)=itererrtemp(i,j)+Fn(i,j,k)**2.d0 *dp3d(i,j,k)
      end do
      itererrtemp(i,j)=sqrt(itererrtemp(i,j))
    end do
    end do

    maxnorminfJ0r0=max(maxval(norminfJ0(:,:)),maxval(norminfr0(:,:)))
    itererr=maxval(itererrtemp(:,:))/maxnorminfJ0r0

    do while ((itercount < maxiter).and.(itererr > itertol))

      info(:,:) = 0
      ! compute the analytic Jacobian
      call get_dirk_jacobian(JacL,JacD,JacU,dt2,dp3d,phi_np1,phis,kappa_star_i,pnh_i,1)
#if (defined COLUMN_OPENMP)
!$omp parallel do private(i,j) collapse(2)
#endif
      do i=1,np
      do j=1,np
        x(1:nlev,i,j) = -Fn(i,j,1:nlev)  !+Fn(i,j,nlev+1:2*nlev,1)/(g*dt2))
        call DGTTRF(nlev, JacL(:,i,j), JacD(:,i,j),JacU(:,i,j),JacU2(:,i,j), Ipiv(:,i,j), info(i,j) )
        ! Tridiagonal solve
        call DGTTRS( 'N', nlev,1, JacL(:,i,j), JacD(:,i,j), JacU(:,i,j), JacU2(:,i,j), Ipiv(:,i,j),x(:,i,j), nlev, info(i,j) )
	! de-update Fn
	Fn(i,j,1:nlev) = Fn(i,j,1:nlev) -  phi_np1(i,j,1:nlev) - (dt2*g)**2 * (1.0-dpnh_dp(i,j,1:nlev))
        ! update approximate solution of phi
        phi_np1(i,j,1:nlev) = phi_np1(i,j,1:nlev) + x(1:nlev,i,j)
      end do
      end do

      !	  de-update w
      elem(ie)%state%w(:,:,1:nlev,np1) = elem(ie)%state%w(:,:,1:nlev,np1) + g*dt2 * &
        (1.0-dpnh_dp(:,:,1:nlev))

      ! compute new dpnh
      call get_pnh_and_exner(hvcoord,theta_dp_cp,dp3d,phi_np1,phis,&
        kappa_star,pnh,dpnh,exner,pnh_i_out=pnh_i)
      dpnh_dp(:,:,:) = dpnh(:,:,:)/dp3d(:,:,:)

      ! update approximate solution of w
      elem(ie)%state%w(:,:,1:nlev,np1) = elem(ie)%state%w(:,:,1:nlev,np1) - g*dt2 * &
        (1.0-dpnh_dp(:,:,1:nlev))
      ! update right-hand side of phi
      Fn(:,:,1:nlev) = Fn(:,:,1:nlev) + phi_np1 + (dt2*g)**2 * (1.0-dpnh_dp(:,:,1:nlev))
 
      ! compute relative errors
      itererrtemp=0.d0
#if (defined COLUMN_OPENMP)
!$omp parallel do private(i,j) collapse(2)
#endif
      do i=1,np
      do j=1,np
        do k=1,nlev
          itererrtemp(i,j)=itererrtemp(i,j)+Fn(i,j,k)**2.d0 *dp3d(i,j,k)
        end do
        itererrtemp(i,j)=sqrt(itererrtemp(i,j))
      end do
      end do
      itererr=maxval(itererrtemp(:,:))/maxnorminfJ0r0

      ! update iteration count and error measure
      itercount=itercount+1
    end do ! end do for the do while loop
!  the following two if-statements are for debugging/testing purposes to track the number of iterations and error attained
!  by the Newton iteration
!      if (itercount > itercountmax) then
!        itercountmax=itercount
!      end if
!      if (itererr > itererrmax) then
!        itererrmax = itererr
!      end if
  end do ! end do for the ie=nets,nete loop
!  maxiter=itercountmax
!  itertol=itererrmax
!  print *, 'max itercount', itercountmax, 'maxitererr ', itererrmax
  call t_stopf('compute_stage_value_dirk')

  end subroutine compute_stage_value_dirk


end module prim_advance_mod

