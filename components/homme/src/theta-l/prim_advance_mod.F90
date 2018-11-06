#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
!
!
!  Man dynamics routines for "theta" nonhydrostatic model
!  Original version: Mark Taylor 2017/1
!  
!  2018/8 TOM sponger layer scaling from P. Lauritzen
!
module prim_advance_mod

  use bndry_mod,          only: bndry_exchangev
  use control_mod,        only: dcmip16_mu, dcmip16_mu_s, hypervis_order, hypervis_subcycle,&
    integration, nu, nu_div, nu_p, nu_s, nu_top, prescribed_wind, qsplit, rsplit, test_case,&
    theta_hydrostatic_mode, tstep_type, use_moisture
  use derivative_mod,     only: derivative_t, divergence_sphere, gradient_sphere, laplace_sphere_wk,&
    laplace_z, vorticity_sphere, vlaplace_sphere_wk 
  use derivative_mod,     only: subcell_div_fluxes, subcell_dss_fluxes
  use dimensions_mod,     only: max_corner_elem, nlev, nlevp, np, qsize
  use edge_mod,           only: edge_g, edgevpack_nlyr, edgevunpack_nlyr
  use edgetype_mod,       only: EdgeBuffer_t,  EdgeDescriptor_t, edgedescriptor_t
  use element_mod,        only: element_t
  use element_state,      only: max_itercnt_perstep,avg_itercnt,max_itererr_perstep, nu_scale_top
  use element_ops,        only: get_temperature, set_theta_ref, state0, get_R_star
  use eos,                only: get_pnh_and_exner,get_theta_from_T,get_phinh,get_dirk_jacobian
  use hybrid_mod,         only: hybrid_t
  use hybvcoord_mod,      only: hvcoord_t
  use kinds,              only: iulog, real_kind
  use perf_mod,           only: t_adj_detailf, t_barrierf, t_startf, t_stopf ! _EXTERNAL
  use parallel_mod,       only: abortmp, global_shared_buf, global_shared_sum, iam, parallel_t
  use physical_constants, only: Cp, cp, cpwater_vapor, g, kappa, Rgas, Rwater_vapor, p0 
  use physics_mod,        only: virtual_specific_heat, virtual_temperature
  use prim_si_mod,        only: preq_vertadv_v1
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
       applycamforcing_ps, applycamforcing_dp3d, &
       applyCAMforcing_dynamics



contains





  subroutine prim_advance_init1(par, elem,integration)
        
    type (parallel_t) :: par
    type (element_t), intent(inout), target   :: elem(:)
    character(len=*)    , intent(in) :: integration
    integer :: i
    integer :: ie


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
    real (kind=real_kind) :: itertol,a1,a2,a3,a4,a5,a6,ahat1,ahat2
    real (kind=real_kind) :: ahat3,ahat4,ahat5,ahat6,dhat1,dhat2,dhat3,dhat4
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
!   tstep_type=6  IMEX-KG243 
!   tstep_type=7  IMEX-KG254
!   

! default weights for computing mean dynamics fluxes
    eta_ave_w = 1d0/qsplit

!   this should not be needed, but in case physics update u without updating w b.c.:
    do ie=nets,nete
       elem(ie)%state%w_i(:,:,nlevp,n0) = (elem(ie)%state%v(:,:,1,nlev,n0)*elem(ie)%derived%gradphis(:,:,1) + &
            elem(ie)%state%v(:,:,2,nlev,n0)*elem(ie)%derived%gradphis(:,:,2))/g
    enddo
 
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
       call compute_andor_apply_rhs(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&                                              
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)                                                      
       ! leapfrog:  u(dt) = u(0) + dt RHS(dt/2)     (store in u(np1))                                                     
       call compute_andor_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&                                               
            deriv,nets,nete,.false.,eta_ave_w,1.d0,1.d0,1.d0)                                                             


    else if (tstep_type==4) then ! explicit table from IMEX-KG254  method                                                              
      call compute_andor_apply_rhs(np1,n0,n0,qn0,dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)
      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt/6,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)
      call compute_andor_apply_rhs(np1,n0,np1,qn0,3*dt/8,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)
      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)
      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w*1d0,1.d0,1.d0,1.d0)



    else if (tstep_type==5) then
       ! Ullrich 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
       ! u1 = u0 + dt/5 RHS(u0)  (save u1 in timelevel nm1)
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
          elem(ie)%state%vtheta_dp(:,:,:,nm1)= (5*elem(ie)%state%vtheta_dp(:,:,:,nm1) &
               - elem(ie)%state%vtheta_dp(:,:,:,n0) )/4
          elem(ie)%state%dp3d(:,:,:,nm1)= (5*elem(ie)%state%dp3d(:,:,:,nm1) &
                  - elem(ie)%state%dp3d(:,:,:,n0) )/4
          elem(ie)%state%w_i(:,:,1:nlevp,nm1)= (5*elem(ie)%state%w_i(:,:,1:nlevp,nm1) &
                  - elem(ie)%state%w_i(:,:,1:nlevp,n0) )/4
          elem(ie)%state%phinh_i(:,:,1:nlev,nm1)= (5*elem(ie)%state%phinh_i(:,:,1:nlev,nm1) &
                  - elem(ie)%state%phinh_i(:,:,1:nlev,n0) )/4
       enddo
       ! u5 = (5*u1/4 - u0/4) + 3dt/4 RHS(u4)
       call compute_andor_apply_rhs(np1,nm1,np1,qn0,3*dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,3*eta_ave_w/4,1.d0,1.d0,1.d0)
       ! final method is the same as:
       ! u5 = u0 +  dt/4 RHS(u0)) + 3dt/4 RHS(u4)
!=========================================================================================
    elseif (tstep_type == 6) then  ! IMEX-KG243

      a1 = 1./4.
      a2 = 1./3.
      a3 = 1./2.
      a4 = 1.0

      ahat4 = 1.
      ahat1 = 0.
      max_itercnt_perstep = 0
      max_itererr_perstep = 0.0  
      ! IMEX-KGNO243
      dhat2 = (1.+sqrt(3.)/3.)/2.
      dhat3 = dhat2
      ahat3 = 1./2.-dhat3
      dhat1 = (ahat3-dhat2+dhat2*dhat3)/(1.-dhat2-dhat3)
      ahat2 = (dhat1-dhat1*dhat3-dhat1*dhat2+dhat1*dhat2*dhat3)/(1.-dhat3)

      call compute_andor_apply_rhs(np1,n0,n0,qn0,dt*a1,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0)
  

      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)
 
      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a2,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat2/a2,1d0)

      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol) 
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)


      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a3,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat3/a3,1d0)

      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a4,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,1d0,ahat4/a4,1d0)

      avg_itercnt = ((nstep)*avg_itercnt + max_itercnt_perstep)/(nstep+1)

!==============================================================================================
    elseif (tstep_type == 7) then 

      max_itercnt_perstep = 0
      max_itererr_perstep = 0.0

      a1 = 1/4d0
      a2 = 1/6d0
      a3 = 3/8d0
      a4 = 1/2d0
      a5 = 1
      ahat1 = 0
      ahat5 = 1

#if 1
      ! IMEX-KGO254 most stable coefficients
      dhat2 = 1d0
      dhat3 = 1d0
      dhat4 = 1d0
      ahat4 = 1d0/2d0-dhat4
      dhat1= (ahat4*ahat5 - ahat5*dhat3 - ahat5*dhat2 + dhat3*dhat2+ dhat3*dhat4 + dhat2*dhat4)/&
        (ahat5-dhat3-dhat2-dhat4)
#endif

#if 0
      ! IMEX-KGO254c coefficients
      dhat2 = 5/6d0
      dhat3 = 5/6d0
      dhat4 = 2/3d0
      ahat4 = 1/2d0-dhat4
      dhat1= (ahat4*ahat5 - ahat5*dhat3 - ahat5*dhat2 + dhat3*dhat2+ dhat3*dhat4 + dhat2*dhat4)/&
           (ahat5-dhat3-dhat2-dhat4)
#endif

#if 0
      ! IMEX-KGO254b coefficients NOT GOOD
      dhat2 = 1./6.
      dhat3 = 1./6.
      dhat4 = 1./6.
      ahat4 = 1./2.-dhat4
      dhat1= (ahat4*ahat5 - ahat5*dhat3 - ahat5*dhat2 + dhat3*dhat2+ dhat3*dhat4 + dhat2*dhat4)/&
           (ahat5-dhat3-dhat2-dhat4)
#endif

      ! IMEX-KG254
      ahat3 = (- ahat4*ahat5*dhat1 - ahat4*ahat5*dhat2+ ahat5*dhat1*dhat2 + ahat5*dhat1*dhat3 +&
        ahat5*dhat2*dhat3- dhat1*dhat2*dhat3 - dhat1*dhat2*dhat4 - dhat1*dhat3*dhat4- &
        dhat2*dhat3*dhat4)/(-ahat4*ahat5)
      ahat2 = ( - ahat3*ahat4*ahat5*dhat1 + ahat4*ahat5*dhat1*dhat2 -&
        ahat5*dhat1*dhat2*dhat3 + dhat1*dhat2*dhat3*dhat4)/(-ahat3*ahat4*ahat5)


      call compute_andor_apply_rhs(np1,n0,n0,qn0,a1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat2/a2,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat3/a3,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a4*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat4/a4,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat4*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a5*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,1d0,ahat5/a5,1d0)

      avg_itercnt = ((nstep)*avg_itercnt + max_itercnt_perstep)/(nstep+1)

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
         call advance_hypervis(elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)




    ! warning: advance_physical_vis currently requires levels that are equally spaced in z
    if (dcmip16_mu>0) call advance_physical_vis(elem,hvcoord,hybrid,deriv,np1,nets,nete,dt,dcmip16_mu_s,dcmip16_mu)

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


!renaming it just to build
!
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
  real (kind=real_kind) :: Rstar(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev)
  real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)

  do ie=nets,nete
     ! apply forcing to Qdp
     elem(ie)%derived%FQps(:,:)=0

     ! apply forcing to temperature
     call get_temperature(elem(ie),temperature,hvcoord,np1)
     do k=1,nlev
        temperature(:,:,k) = temperature(:,:,k) + dt*elem(ie)%derived%FT(:,:,k)
     enddo


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
     do k=1,nlev
        dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
             ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
     enddo
     do q=1,qsize
        do k=1,nlev
           elem(ie)%state%Q(:,:,k,q) = elem(ie)%state%Qdp(:,:,k,q,np1_qdp)/dp(:,:,k)
        enddo
     enddo

     ! now that we have updated Qdp and dp, compute vtheta_dp from temperature
     call get_R_star(Rstar,elem(ie)%state%Q(:,:,:,1))
     call get_theta_from_T(hvcoord,Rstar,temperature,dp,&
          elem(ie)%state%phinh_i(:,:,:,np1),elem(ie)%state%vtheta_dp(:,:,:,np1))


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
     elem(ie)%state%w_i(:,:,1:nlev,np1) = elem(ie)%state%w_i(:,:,1:nlev,np1) + dt*elem(ie)%derived%FM(:,:,3,:)

     ! finally update w at the surface: 
     elem(ie)%state%w_i(:,:,nlevp,np1) = (elem(ie)%state%v(:,:,1,nlev,np1)*elem(ie)%derived%gradphis(:,:,1) + &
          elem(ie)%state%v(:,:,2,nlev,np1)*elem(ie)%derived%gradphis(:,:,2))/g
  enddo
  
  end subroutine applyCAMforcing_dynamics





  subroutine advance_hypervis(elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
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
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord

  real (kind=real_kind) :: dt2
  integer :: nets,nete

  ! local
  real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
  integer :: k2,k,kptr,i,j,ie,ic,nt,nlyr_tot,ssize
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
  real (kind=real_kind) :: phi_ref(np,np,nlevp,nets:nete)
  real (kind=real_kind) :: dp_ref(np,np,nlev,nets:nete)

  call t_startf('advance_hypervis')

  dt=dt2/hypervis_subcycle

  if (theta_hydrostatic_mode) then
     nlyr_tot=4*nlev        ! dont bother to dss w_i and phinh_i
     ssize=2*nlev
  else
     nlyr_tot=6*nlev  ! total amount of data for DSS
     ssize=4*nlev
  endif
  
  do k=1,nlev
     exner0(k) = (hvcoord%etam(k)*hvcoord%ps0/p0 )**kappa
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTE1:  Diffusion works best when applied to theta.
! It creates some TOM noise when applied to vtheta_dp in DCMIP 2.0 test
! so we convert from vtheta_dp->theta, and then convert back at the end of diffusion
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


#if 0
     ! phi_ref,theta_ref depend only on ps:
     call set_theta_ref(hvcoord,dp_ref(:,:,:,ie),theta_ref(:,:,:,ie))
     exner(:,:,:)=theta_ref(:,:,:,ie)*dp_ref(:,:,:,ie) ! use as temp array
     call get_phinh(hvcoord,elem(ie)%state%phis,&
          exner(:,:,:),dp_ref(:,:,:,ie),phi_ref(:,:,:,ie))
#endif
#if 1
     ! phi_ref depends only on ps, theta_ref depends on dp3d
     call set_theta_ref(hvcoord,dp_ref(:,:,:,ie),theta_ref(:,:,:,ie))
     exner(:,:,:)=theta_ref(:,:,:,ie)*dp_ref(:,:,:,ie) ! use as temp array
     call get_phinh(hvcoord,elem(ie)%state%phis,&
          exner(:,:,:),dp_ref(:,:,:,ie),phi_ref(:,:,:,ie))

     call set_theta_ref(hvcoord,elem(ie)%state%dp3d(:,:,:,nt),theta_ref(:,:,:,ie))
#endif
#if 0
     ! no reference state, for testing
     theta_ref(:,:,:,ie)=0
     phi_ref(:,:,:,ie)=0
     dp_ref(:,:,:,ie)=0
#endif


     ! convert vtheta_dp -> theta
     do k=1,nlev
        elem(ie)%state%vtheta_dp(:,:,k,nt)=&
             elem(ie)%state%vtheta_dp(:,:,k,nt)/elem(ie)%state%dp3d(:,:,k,nt)
     enddo
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
           elem(ie)%state%vtheta_dp(:,:,k,nt)=elem(ie)%state%vtheta_dp(:,:,k,nt)-&
                theta_ref(:,:,k,ie)
           elem(ie)%state%phinh_i(:,:,k,nt)=elem(ie)%state%phinh_i(:,:,k,nt)-&
                phi_ref(:,:,k,ie)
           elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt)-&
                dp_ref(:,:,k,ie)
        enddo
     enddo
     
     call biharmonic_wk_theta(elem,stens,vtens,deriv,edge_g,hybrid,nt,nets,nete)
     
     do ie=nets,nete
        
        ! comptue mean flux
        if (nu_p>0) then
           elem(ie)%derived%dpdiss_ave(:,:,:)=elem(ie)%derived%dpdiss_ave(:,:,:)+&
                (eta_ave_w*elem(ie)%state%dp3d(:,:,:,nt)+dp_ref(:,:,:,ie))/hypervis_subcycle
           elem(ie)%derived%dpdiss_biharmonic(:,:,:)=elem(ie)%derived%dpdiss_biharmonic(:,:,:)+&
                eta_ave_w*stens(:,:,:,1,ie)/hypervis_subcycle
        endif
           do k=1,nlev
              ! advace in time.
              ! note: DSS commutes with time stepping, so we can time advance and then DSS.
              ! note: weak operators alreayd have mass matrix "included"

              ! biharmonic terms need a negative sign:
              if (nu_top>0 .and. nu_scale_top(k)>1) then
                 ! add regular diffusion near top
                 lap_s(:,:,1)=laplace_sphere_wk(elem(ie)%state%dp3d       (:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_s(:,:,2)=laplace_sphere_wk(elem(ie)%state%vtheta_dp  (:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_s(:,:,3)=laplace_sphere_wk(elem(ie)%state%w_i        (:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_s(:,:,4)=laplace_sphere_wk(elem(ie)%state%phinh_i    (:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)

                 vtens(:,:,:,k,ie)=(  -nu*vtens(:,:,:,k,ie) + nu_scale_top(k)*nu_top*lap_v(:,:,:)) ! u and v
                 stens(:,:,k,1,ie)=(-nu_p*stens(:,:,k,1,ie) + nu_scale_top(k)*nu_top*lap_s(:,:,1)) ! dp3d
                 stens(:,:,k,2,ie)=(  -nu*stens(:,:,k,2,ie) + nu_scale_top(k)*nu_top*lap_s(:,:,2)) ! theta
                 stens(:,:,k,3,ie)=(  -nu*stens(:,:,k,3,ie) + nu_scale_top(k)*nu_top*lap_s(:,:,3)) ! w
                 stens(:,:,k,4,ie)=(-nu_s*stens(:,:,k,4,ie) + nu_scale_top(k)*nu_top*lap_s(:,:,4)) ! phi
              else
                 vtens(:,:,:,k,ie)=-nu  *vtens(:,:,:,k,ie) ! u,v
                 stens(:,:,k,1,ie)=-nu_p*stens(:,:,k,1,ie) ! dp3d
                 stens(:,:,k,2,ie)=-nu  *stens(:,:,k,2,ie) ! theta
                 stens(:,:,k,3,ie)=-nu  *stens(:,:,k,3,ie) ! w
                 stens(:,:,k,4,ie)=-nu_s*stens(:,:,k,4,ie) ! phi
              endif

           enddo

           kptr=0;      call edgeVpack_nlyr(edge_g,elem(ie)%desc,vtens(:,:,:,:,ie),2*nlev,kptr,nlyr_tot)
           kptr=2*nlev; call edgeVpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,:,ie),ssize,kptr,nlyr_tot)

        enddo

        call t_startf('ahdp_bexchV2')
        call bndry_exchangeV(hybrid,edge_g)
        call t_stopf('ahdp_bexchV2')

        do ie=nets,nete

           kptr=0
           call edgeVunpack_nlyr(edge_g,elem(ie)%desc,vtens(:,:,:,:,ie),2*nlev,kptr,nlyr_tot)
           kptr=2*nlev
           call edgeVunpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,:,ie),ssize,kptr,nlyr_tot)


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
              elem(ie)%state%vtheta_dp(:,:,k,nt)=elem(ie)%state%vtheta_dp(:,:,k,nt)+&
                   theta_ref(:,:,k,ie)
              elem(ie)%state%phinh_i(:,:,k,nt)=elem(ie)%state%phinh_i(:,:,k,nt)+&
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
           elem(ie)%state%w_i(:,:,k,nt)=elem(ie)%state%w_i(:,:,k,nt) &
                +stens(:,:,k,3,ie)
           
           elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt) &
                +stens(:,:,k,1,ie)
           
           elem(ie)%state%phinh_i(:,:,k,nt)=elem(ie)%state%phinh_i(:,:,k,nt) &
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
        ! use hydrostatic pressure for simplicity
        p_i(:,:,1)=hvcoord%hyai(1)*hvcoord%ps0
        do k=1,nlev
           p_i(:,:,k+1)=p_i(:,:,k) + elem(ie)%state%dp3d(:,:,k,nt)
        enddo
        !call get_pnh_and_exner(hvcoord,elem(ie)%state%vtheta_dp(:,:,:,nt),&
        !     elem(ie)%state%dp3d(:,:,:,nt),elem(ie)%state%phinh_i(:,:,:,nt),&
        !     pnh,exner,dpnh_dp_i)

        do k=1,nlev
           ! for w averaging, we didn't compute dissipation at surface, so just use one level
           k2=max(k,nlev)
           ! p(:,:,k) = (p_i(:,:,k) + elem(ie)%state%dp3d(:,:,k,nt)/2)
           exner(:,:,k)  = ( (p_i(:,:,k) + elem(ie)%state%dp3d(:,:,k,nt)/2) /p0)**kappa
           if (theta_hydrostatic_mode) then
              heating(:,:,k)= (elem(ie)%state%v(:,:,1,k,nt)*vtens(:,:,1,k,ie) + &
                   elem(ie)%state%v(:,:,2,k,nt)*vtens(:,:,2,k,ie) ) / &
                   (exner(:,:,k)*Cp)

           else
              heating(:,:,k)= (elem(ie)%state%v(:,:,1,k,nt)*vtens(:,:,1,k,ie) + &
                   elem(ie)%state%v(:,:,2,k,nt)*vtens(:,:,2,k,ie)  +&
                   (elem(ie)%state%w_i(:,:,k,nt)*stens(:,:,k,3,ie)  +&
                     elem(ie)%state%w_i(:,:,k2,nt)*stens(:,:,k2,3,ie))/2 ) /  +&
                   (exner(:,:,k)*Cp)  
           endif
           
           elem(ie)%state%vtheta_dp(:,:,k,nt)=elem(ie)%state%vtheta_dp(:,:,k,nt) &
                +stens(:,:,k,2,ie)*hvcoord%dp0(k)*exner0(k)/(exner(:,:,k)*elem(ie)%state%dp3d(:,:,k,nt)&
                ) ! -heating(:,:,k)
        enddo
        
     enddo
  enddo

! convert vtheta_dp -> theta
  do ie=nets,nete            
     elem(ie)%state%vtheta_dp(:,:,:,nt)=&
          elem(ie)%state%vtheta_dp(:,:,:,nt)*elem(ie)%state%dp3d(:,:,:,nt)
     
     ! finally update w at the surface: 
     elem(ie)%state%w_i(:,:,nlevp,nt) = (elem(ie)%state%v(:,:,1,nlev,nt)*elem(ie)%derived%gradphis(:,:,1) + &
          elem(ie)%state%v(:,:,2,nlev,nt)*elem(ie)%derived%gradphis(:,:,2))/g
  enddo	


  call t_stopf('advance_hypervis')

  end subroutine advance_hypervis





  subroutine advance_physical_vis(elem,hvcoord,hybrid,deriv,nt,nets,nete,dt,mu_s,mu)
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
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord

  real (kind=real_kind) :: dt, mu_s, mu
  integer :: nt,nets,nete

  ! local
  integer :: k,kptr,ie
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,4,nets:nete)      :: stens  ! dp3d,theta,w,phi
  real (kind=real_kind), dimension(np,np,nlevp,2,nets:nete)     :: stens_i 


  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: delz,delz_i

  real (kind=real_kind) :: theta_ref(np,np,nlev)

  real (kind=real_kind) :: theta_prime(np,np,nlev)
  real (kind=real_kind) :: phi_prime(np,np,nlevp)
  real (kind=real_kind) :: dp_prime(np,np,nlev)
  real (kind=real_kind) :: w_prime(np,np,nlevp)
  real (kind=real_kind) :: u_prime(np,np,2,nlev)

  !if(test_case .ne. 'dcmip2016_test3') call abortmp("dcmip16_mu is currently limited to dcmip16 test 3")

  call t_startf('advance_physical_vis')
  delz = 20d3/nlev
  delz_i = 20d3/nlevp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute reference states
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie=nets,nete

     theta_ref(:,:,:) = state0(ie)%vtheta_dp(:,:,:,1)/state0(ie)%dp3d(:,:,:,1)
     elem(ie)%state%vtheta_dp(:,:,:,nt)=&
             elem(ie)%state%vtheta_dp(:,:,:,nt)/elem(ie)%state%dp3d(:,:,:,nt)

     ! perturbation variables
     u_prime(:,:,:,:)  = elem(ie)%state%v(:,:,:,:,nt)          -state0(ie)%v(:,:,:,:,1)
     w_prime(:,:,:)    = elem(ie)%state%w_i(:,:,:,nt)     -state0(ie)%w_i(:,:,:,1)
     dp_prime(:,:,:)   = elem(ie)%state%dp3d(:,:,:,nt)         -state0(ie)%dp3d(:,:,:,1)
     phi_prime(:,:,:)  = elem(ie)%state%phinh_i(:,:,:,nt) -state0(ie)%phinh_i(:,:,:,1)
     theta_prime(:,:,:)= elem(ie)%state%vtheta_dp(:,:,:,nt)-theta_ref(:,:,:)

     ! vertical viscosity
     call laplace_z(u_prime,    vtens(:,:,:,:,ie),2,nlev,delz)
     call laplace_z(dp_prime,   stens(:,:,:,1,ie),1,nlev,delz)
     call laplace_z(theta_prime,stens(:,:,:,2,ie),1,nlev,delz)
     call laplace_z(w_prime,    stens_i(:,:,:,1,ie),1,nlevp,delz_i)
     call laplace_z(phi_prime,  stens_i(:,:,:,2,ie),1,nlevp,delz_i)

     ! add in horizontal viscosity
     ! multiply by mass matrix for DSS
     ! horiz viscosity already has mass matrix built in
     ! for interface quantities, only use 1:nlev (dont apply at surface)
     do k=1,nlev
        lap_v = vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)

        vtens(:,:,1,k,ie) = (vtens(:,:,1,k,ie)*elem(ie)%spheremp(:,:) + lap_v(:,:,1))
        vtens(:,:,2,k,ie) = (vtens(:,:,2,k,ie)*elem(ie)%spheremp(:,:) + lap_v(:,:,2))

        stens(:,:,k,1,ie) = (stens(:,:,k,1,ie)*elem(ie)%spheremp(:,:) + &
             laplace_sphere_wk(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),var_coef=.false.)  )

        stens(:,:,k,2,ie) = (stens(:,:,k,2,ie)*elem(ie)%spheremp(:,:) + &
             laplace_sphere_wk(elem(ie)%state%vtheta_dp(:,:,k,nt),deriv,elem(ie),var_coef=.false.)  )

        stens(:,:,k,3,ie) = (stens_i(:,:,k,1,ie)*elem(ie)%spheremp(:,:) + &
             laplace_sphere_wk(elem(ie)%state%w_i(:,:,k,nt),deriv,elem(ie),var_coef=.false.) )

        stens(:,:,k,4,ie) = (stens_i(:,:,k,2,ie)*elem(ie)%spheremp(:,:) + &
             laplace_sphere_wk(elem(ie)%state%phinh_i(:,:,k,nt),deriv,elem(ie),var_coef=.false.) ) 

     enddo

     kptr=0
     call edgeVpack_nlyr(edge_g,elem(ie)%desc,vtens(:,:,:,:,ie),2*nlev,kptr,6*nlev)
     kptr=2*nlev
     call edgeVpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,:,ie),4*nlev,kptr,6*nlev)
     
  enddo

  call bndry_exchangeV(hybrid,edge_g)
  
  do ie=nets,nete
     
     kptr=0
     call edgeVunpack_nlyr(edge_g,elem(ie)%desc,vtens(:,:,:,:,ie),2*nlev,kptr,6*nlev)
     kptr=2*nlev
     call edgeVunpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,:,ie),4*nlev,kptr,6*nlev)
     
     ! apply inverse mass matrix, accumulate tendencies
     do k=1,nlev
        elem(ie)%state%v(:,:,1,k,nt)=elem(ie)%state%v(:,:,1,k,nt) + &
             mu*dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)

        elem(ie)%state%v(:,:,2,k,nt)=elem(ie)%state%v(:,:,2,k,nt) + &
             mu*dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
        
        elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt) &
             +mu_s*dt*stens(:,:,k,1,ie)*elem(ie)%rspheremp(:,:)
        
        elem(ie)%state%vtheta_dp(:,:,k,nt)=elem(ie)%state%vtheta_dp(:,:,k,nt) &
             +mu_s*dt*stens(:,:,k,2,ie)*elem(ie)%rspheremp(:,:)

        elem(ie)%state%w_i(:,:,k,nt)=elem(ie)%state%w_i(:,:,k,nt) &
             +mu_s*dt*stens(:,:,k,3,ie)*elem(ie)%rspheremp(:,:)
        
        elem(ie)%state%phinh_i(:,:,k,nt)=elem(ie)%state%phinh_i(:,:,k,nt) &
             +mu_s*dt*stens(:,:,k,4,ie)*elem(ie)%rspheremp(:,:)
        
     enddo
  enddo


  ! convert vtheta_dp -> theta
  do ie=nets,nete            
     do k=1,nlev
        elem(ie)%state%vtheta_dp(:,:,k,nt)=&
             elem(ie)%state%vtheta_dp(:,:,k,nt)*elem(ie)%state%dp3d(:,:,k,nt)
     enddo

     ! finally update w at the surface: 
     elem(ie)%state%w_i(:,:,nlevp,nt) = (elem(ie)%state%v(:,:,1,nlev,nt)*elem(ie)%derived%gradphis(:,:,1) + &
          elem(ie)%state%v(:,:,2,nlev,nt)*elem(ie)%derived%gradphis(:,:,2))/g
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
  real (kind=real_kind), pointer, dimension(:,:,:) :: phi_i
  real (kind=real_kind), pointer, dimension(:,:,:) :: dp3d
  real (kind=real_kind), pointer, dimension(:,:,:) :: vtheta_dp
   
  real (kind=real_kind) :: vtheta(np,np,nlev)
  real (kind=real_kind) :: vtheta_i(np,np,nlevp)
  real (kind=real_kind) :: omega_i(np,np,nlevp)
  real (kind=real_kind) :: omega(np,np,nlev)
  real (kind=real_kind) :: vort(np,np,nlev)           ! vorticity
  real (kind=real_kind) :: divdp(np,np,nlev)     
  real (kind=real_kind) :: phi(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev)               ! nh (nonydro) pressure
  real (kind=real_kind) :: dp3d_i(np,np,nlevp)
  real (kind=real_kind) :: exner(np,np,nlev)         ! exner nh pressure
  real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)    ! dpnh / dp3d at interfaces
  real (kind=real_kind) :: eta_dot_dpdn(np,np,nlevp)  ! vertical velocity at interfaces
  real (kind=real_kind) :: KE(np,np,nlev)             ! Kinetic energy
  real (kind=real_kind) :: gradexner(np,np,2,nlev)    ! grad(p^kappa)
  real (kind=real_kind) :: gradphinh_i(np,np,2,nlevp) ! gradphi at interfaces
  real (kind=real_kind) :: mgrad(np,np,2,nlev)        ! gradphi metric term at cell centers
  real (kind=real_kind) :: gradKE(np,np,2,nlev)       ! grad(0.5 u^T u )
  real (kind=real_kind) :: wvor(np,np,2,nlev)         ! w vorticity term

  real (kind=real_kind) :: gradw_i(np,np,2,nlevp)
  real (kind=real_kind) :: v_gradw_i(np,np,nlevp)     
  real (kind=real_kind) :: v_gradtheta(np,np,nlev)     
  real (kind=real_kind) :: v_theta(np,np,2,nlev)
  real (kind=real_kind) :: div_v_theta(np,np,nlev)
  real (kind=real_kind) :: v_gradphinh_i(np,np,nlevp) ! v*gradphi at interfaces
  real (kind=real_kind) :: v_i(np,np,2,nlevp)

  real (kind=real_kind) :: v_vadv(np,np,2,nlev)     ! velocity vertical advection
  real (kind=real_kind) :: theta_vadv(np,np,nlev)   ! w,phi, theta  vertical advection term
  real (kind=real_kind) :: w_vadv_i(np,np,nlevp)      ! w,phi, theta  vertical advection term
  real (kind=real_kind) :: phi_vadv_i(np,np,nlevp)    ! w,phi, theta  vertical advection term

  real (kind=real_kind) :: vtens1(np,np,nlev)
  real (kind=real_kind) :: vtens2(np,np,nlev)
  real (kind=real_kind) :: stens(np,np,nlev,3) ! tendencies w,phi,theta
                                               ! w,phi tendencies not computed at nlevp
  real (kind=real_kind) :: w_tens(np,np,nlevp)  ! need to update w at surface as well
  real (kind=real_kind) :: theta_tens(np,np,nlev)
  real (kind=real_kind) :: phi_tens(np,np,nlevp)
                                               

  real (kind=real_kind) :: pi(np,np,nlev)                ! hydrostatic pressure
  real (kind=real_kind) :: pi_i(np,np,nlevp)             ! hydrostatic pressure interfaces
  real (kind=real_kind), dimension(np,np,nlev) :: vgrad_p

  real (kind=real_kind) ::  temp(np,np,nlev)
  real (kind=real_kind) ::  vtemp(np,np,2,nlev)       ! generic gradient storage
  real (kind=real_kind), dimension(np,np) :: sdot_sum ! temporary field
  real (kind=real_kind) ::  v1,v2,w,d_eta_dot_dpdn_dn
  integer :: i,j,k,kptr,ie, nlyr_tot


  call t_startf('compute_andor_apply_rhs')

  if (theta_hydrostatic_mode) then
     nlyr_tot=4*nlev        ! dont bother to dss w_i and phinh_i
  else
     nlyr_tot=5*nlev+nlevp  ! total amount of data for DSS
  endif
     

  do ie=nets,nete
#if 0
     if (.not. theta_hydrostatic_mode) then
        temp(:,:,1) =  (elem(ie)%state%v(:,:,1,nlev,n0)*elem(ie)%derived%gradphis(:,:,1) + &
             elem(ie)%state%v(:,:,2,nlev,n0)*elem(ie)%derived%gradphis(:,:,2))/g
        if ( maxval(abs(temp(:,:,1)-elem(ie)%state%w_i(:,:,nlevp,n0))) >1e-10) then
           write(iulog,*) 'WARNING: w(n0) does not satisfy b.c.'
           write(iulog,*) 'val1 = ',temp(:,:,1)
           write(iulog,*) 'val2 = ',elem(ie)%state%w_i(:,:,nlevp,n0)
           write(iulog,*) 'diff: ',temp(:,:,1)-elem(ie)%state%w_i(:,:,nlevp,n0)
        endif
        ! w boundary condition. just in case:
        elem(ie)%state%w_i(:,:,nlevp,n0) = (elem(ie)%state%v(:,:,1,nlev,n0)*elem(ie)%derived%gradphis(:,:,1) + &
             elem(ie)%state%v(:,:,2,nlev,n0)*elem(ie)%derived%gradphis(:,:,2))/g
     endif
#endif
     dp3d  => elem(ie)%state%dp3d(:,:,:,n0)
     vtheta_dp  => elem(ie)%state%vtheta_dp(:,:,:,n0)
     vtheta(:,:,:) = vtheta_dp(:,:,:)/dp3d(:,:,:)
     phi_i => elem(ie)%state%phinh_i(:,:,:,n0)

     call get_pnh_and_exner(hvcoord,vtheta_dp,dp3d,phi_i,pnh,exner,dpnh_dp_i)

     dp3d_i(:,:,1) = dp3d(:,:,1)
     dp3d_i(:,:,nlevp) = dp3d(:,:,nlev)
     do k=2,nlev
        dp3d_i(:,:,k)=(dp3d(:,:,k)+dp3d(:,:,k-1))/2
     end do

     ! special averaging for velocity for energy conservation
     v_i(:,:,1:2,1) = elem(ie)%state%v(:,:,1:2,1,n0)  
     v_i(:,:,1:2,nlevp) = elem(ie)%state%v(:,:,1:2,nlev,n0)
     do k=2,nlev
        v_i(:,:,1,k) = (dp3d(:,:,k)*elem(ie)%state%v(:,:,1,k,n0) + &
             dp3d(:,:,k-1)*elem(ie)%state%v(:,:,1,k-1,n0) ) / (2*dp3d_i(:,:,k))
        v_i(:,:,2,k) = (dp3d(:,:,k)*elem(ie)%state%v(:,:,2,k,n0) + &
             dp3d(:,:,k-1)*elem(ie)%state%v(:,:,2,k-1,n0) ) / (2*dp3d_i(:,:,k))
     end do
     
     if (theta_hydrostatic_mode) then
        do k=nlev,1,-1          ! traditional Hydrostatic integral
           phi_i(:,:,k)=phi_i(:,:,k+1)+&
                Rgas*vtheta_dp(:,:,k)*exner(:,:,k)/pnh(:,:,k)
        enddo
        ! in H mode, ignore w contibutions to KE term
        ! set to zero so H and NH can share code and reduce if statements
        elem(ie)%state%w_i(:,:,:,n0)=0   
     endif

     do k=1,nlev
        phi(:,:,k) = (phi_i(:,:,k)+phi_i(:,:,k+1))/2  ! for diagnostics

        ! ================================
        ! Accumulate mean Vel_rho flux in vn0
        ! ================================
        vtemp(:,:,1,k) = elem(ie)%state%v(:,:,1,k,n0)*dp3d(:,:,k)
        vtemp(:,:,2,k) = elem(ie)%state%v(:,:,2,k,n0)*dp3d(:,:,k)
        elem(ie)%derived%vn0(:,:,:,k)=elem(ie)%derived%vn0(:,:,:,k)+eta_ave_w*vtemp(:,:,:,k)

        divdp(:,:,k)=divergence_sphere(vtemp(:,:,:,k),deriv,elem(ie))
        vort(:,:,k)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,n0),deriv,elem(ie))
     enddo

     ! Compute omega =  Dpi/Dt   Used only as a DIAGNOSTIC
     ! for historical reasons, we actually compute w/pi
     pi_i(:,:,1)=hvcoord%hyai(1)*hvcoord%ps0
     omega_i(:,:,1)=0
     do k=1,nlev
        pi_i(:,:,k+1)=pi_i(:,:,k) + dp3d(:,:,k)
        omega_i(:,:,k+1)=omega_i(:,:,k)+divdp(:,:,k)
     enddo
     do k=1,nlev
        pi(:,:,k)=pi_i(:,:,k) + dp3d(:,:,k)/2
        vtemp(:,:,:,k) = gradient_sphere( pi(:,:,k), deriv, elem(ie)%Dinv);
        vgrad_p(:,:,k) = elem(ie)%state%v(:,:,1,k,n0)*vtemp(:,:,1,k)+&
             elem(ie)%state%v(:,:,2,k,n0)*vtemp(:,:,2,k)
        omega(:,:,k) = (vgrad_p(:,:,k) - ( omega_i(:,:,k)+omega_i(:,:,k+1))/2) 
     enddo        

     ! ==================================================
     ! Compute eta_dot_dpdn
     ! save sdot_sum as this is the -RHS of ps_v equation
     ! ==================================================
     if (rsplit>0) then
        ! VERTICALLY LAGRANGIAN:   no vertical motion
        eta_dot_dpdn=0
        w_vadv_i=0
        phi_vadv_i=0
        theta_vadv=0
        v_vadv=0
     else
        sdot_sum=0
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

        eta_dot_dpdn(:,:,1     ) = 0
        eta_dot_dpdn(:,:,nlevp)  = 0
        vtheta_i(:,:,1) =       0
        vtheta_i(:,:,nlevp) =   0


        ! ===========================================================
        ! Compute vertical advection of v from eq. CCM2 (3.b.1)
        ! ==============================================
        call preq_vertadv_v1(elem(ie)%state%v(:,:,:,:,n0),eta_dot_dpdn,dp3d,v_vadv)

        ! compute (cp*theta) at interfaces
        ! for energy conservation, use averaging consistent with EOS
        ! dont bother to compute at surface and top since it will be multiplied by eta-dot
#if 0           
           do k=2,nlev  ! simple averaging
              vtheta_i(:,:,k) = (vtheta(:,:,k)+vtheta(:,:,k-1))/2
           enddo
#else
           ! E conserving average, but much more dissipative
           do k=2,nlev
              vtheta_i(:,:,k) = -dpnh_dp_i(:,:,k)*(phi(:,:,k)-phi(:,:,k-1))/&
                   (exner(:,:,k)-exner(:,:,k-1)) / Cp
           enddo
#endif           



        do k=1,nlev
           ! average interface quantity to midpoints:
           temp(:,:,k) = (( eta_dot_dpdn(:,:,k)+eta_dot_dpdn(:,:,k+1))/2)*&
                (elem(ie)%state%w_i(:,:,k+1,n0)-elem(ie)%state%w_i(:,:,k,n0))
           
           ! theta vadv term at midoints
           theta_vadv(:,:,k)= eta_dot_dpdn(:,:,k+1)*vtheta_i(:,:,k+1) - &
                eta_dot_dpdn(:,:,k)*vtheta_i(:,:,k)
        enddo
        ! compute ave( ave(etadot) d/dx )
        do k=2,nlev
           w_vadv_i(:,:,k)  =(temp(:,:,k-1)+temp(:,:,k))/2
           phi_vadv_i(:,:,k)=eta_dot_dpdn(:,:,k)*(phi(:,:,k)-phi(:,:,k-1))
        end do
        w_vadv_i(:,:,1) = temp(:,:,1)
        w_vadv_i(:,:,nlevp) = temp(:,:,nlev)
        phi_vadv_i(:,:,1) = 0
        phi_vadv_i(:,:,nlevp) = 0

        ! final form of SB81 vertical advection operator:
        w_vadv_i=w_vadv_i/dp3d_i
        phi_vadv_i=phi_vadv_i/dp3d_i
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
             elem(ie)%derived%omega_p(:,:,k) + eta_ave_w*omega(:,:,k)
     enddo
     elem(ie)%derived%eta_dot_dpdn(:,:,nlev+1) = &
             elem(ie)%derived%eta_dot_dpdn(:,:,nlev+1) + eta_ave_w*eta_dot_dpdn(:,:,nlev+1)

#if (defined COLUMN_OPENMP)
 !$omp parallel do private(k)
#endif
     ! ================================================
     ! w,phi tendencies including surface
     ! ================================================  
     do k=1,nlev
        ! compute gradphi at interfaces and then average to levels
        gradphinh_i(:,:,:,k)   = gradient_sphere(phi_i(:,:,k),deriv,elem(ie)%Dinv)   
           
        gradw_i(:,:,:,k)   = gradient_sphere(elem(ie)%state%w_i(:,:,k,n0),deriv,elem(ie)%Dinv)
        v_gradw_i(:,:,k) = v_i(:,:,1,k)*gradw_i(:,:,1,k) + v_i(:,:,2,k)*gradw_i(:,:,2,k)
        ! w - tendency on interfaces 
        w_tens(:,:,k) = (-w_vadv_i(:,:,k) - v_gradw_i(:,:,k))*scale1 - scale2*g*(1-dpnh_dp_i(:,:,k) )

        ! phi - tendency on interfaces
        v_gradphinh_i(:,:,k) = v_i(:,:,1,k)*gradphinh_i(:,:,1,k) &
             +v_i(:,:,2,k)*gradphinh_i(:,:,2,k) 
        phi_tens(:,:,k) =  (-phi_vadv_i(:,:,k) - v_gradphinh_i(:,:,k))*scale1 &
          + scale2*g*elem(ie)%state%w_i(:,:,k,n0)
     end do

     ! k =nlevp case, all terms in the imex methods are treated explicitly at the boundary
     k =nlevp 
    ! compute gradphi at interfaces and then average to levels
    gradphinh_i(:,:,:,k)   = gradient_sphere(phi_i(:,:,k),deriv,elem(ie)%Dinv)
    gradw_i(:,:,:,k)   = gradient_sphere(elem(ie)%state%w_i(:,:,k,n0),deriv,elem(ie)%Dinv)
    v_gradw_i(:,:,k) = v_i(:,:,1,k)*gradw_i(:,:,1,k) + v_i(:,:,2,k)*gradw_i(:,:,2,k)
    ! w - tendency on interfaces
    w_tens(:,:,k) = (-w_vadv_i(:,:,k) - v_gradw_i(:,:,k))*scale1 - scale1*g*(1-dpnh_dp_i(:,:,k) )

    ! phi - tendency on interfaces
    v_gradphinh_i(:,:,k) = v_i(:,:,1,k)*gradphinh_i(:,:,1,k) &
     +v_i(:,:,2,k)*gradphinh_i(:,:,2,k)
    phi_tens(:,:,k) =  (-phi_vadv_i(:,:,k) - v_gradphinh_i(:,:,k))*scale1 &
    + scale1*g*elem(ie)%state%w_i(:,:,k,n0)
    




#if (defined COLUMN_OPENMP)
 !$omp parallel do private(k,i,j,v1,v2)                                                                           
#endif
     ! ================================================                                                                 
     ! v1,v2 tendencies:                                                                                          
     ! ================================================           
     do k=1,nlev
        ! theta - tendency on levels
        v_theta(:,:,1,k)=elem(ie)%state%v(:,:,1,k,n0)*vtheta_dp(:,:,k)
        v_theta(:,:,2,k)=elem(ie)%state%v(:,:,2,k,n0)*vtheta_dp(:,:,k)
        div_v_theta(:,:,k)=divergence_sphere(v_theta(:,:,:,k),deriv,elem(ie))
        theta_tens(:,:,k)=(-theta_vadv(:,:,k)-div_v_theta(:,:,k))*scale1

        ! w vorticity correction term
        temp(:,:,k) = (elem(ie)%state%w_i(:,:,k,n0)**2 + &
             elem(ie)%state%w_i(:,:,k+1,n0)**2)/4
        wvor(:,:,:,k) = gradient_sphere(temp(:,:,k),deriv,elem(ie)%Dinv)
        wvor(:,:,1,k) = wvor(:,:,1,k) - (elem(ie)%state%w_i(:,:,k,n0)*gradw_i(:,:,1,k) +&
             elem(ie)%state%w_i(:,:,k+1,n0)*gradw_i(:,:,1,k+1))/2
        wvor(:,:,2,k) = wvor(:,:,2,k) - (elem(ie)%state%w_i(:,:,k,n0)*gradw_i(:,:,2,k) +&
             elem(ie)%state%w_i(:,:,k+1,n0)*gradw_i(:,:,2,k+1))/2

        KE(:,:,k) = ( elem(ie)%state%v(:,:,1,k,n0)**2 + elem(ie)%state%v(:,:,2,k,n0)**2)/2
        gradKE(:,:,:,k) = gradient_sphere(KE(:,:,k),deriv,elem(ie)%Dinv)
        gradexner(:,:,:,k) = gradient_sphere(exner(:,:,k),deriv,elem(ie)%Dinv)

        ! special averaging of dpnh/dpi grad(phi) for E conservation
        mgrad(:,:,1,k) = (dpnh_dp_i(:,:,k)*gradphinh_i(:,:,1,k)+ &
              dpnh_dp_i(:,:,k+1)*gradphinh_i(:,:,1,k+1))/2
        mgrad(:,:,2,k) = (dpnh_dp_i(:,:,k)*gradphinh_i(:,:,2,k)+ &
              dpnh_dp_i(:,:,k+1)*gradphinh_i(:,:,2,k+1))/2


        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)

              vtens1(i,j,k) = (-v_vadv(i,j,1,k) &
                   + v2*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - gradKE(i,j,1,k) - mgrad(i,j,1,k) &
                  -Cp*vtheta(i,j,k)*gradexner(i,j,1,k)&
                  -wvor(i,j,1,k) )*scale1


              vtens2(i,j,k) = (-v_vadv(i,j,2,k) &
                   - v1*(elem(ie)%fcor(i,j) + vort(i,j,k)) &
                   - gradKE(i,j,2,k) - mgrad(i,j,2,k) &
                  -Cp*vtheta(i,j,k)*gradexner(i,j,2,k) &
                  -wvor(i,j,2,k) )*scale1
           end do
        end do     
     end do 



     
#ifdef ENERGY_DIAGNOSTICS
     ! =========================================================
     ! diagnostics. not performance critical, dont thread
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

        do k =1,nlev
          do j=1,np
            do i=1,np                
               d_eta_dot_dpdn_dn=(eta_dot_dpdn(i,j,k+1)-eta_dot_dpdn(i,j,k))
               !  Form horiz advection of KE-u
               elem(ie)%accum%KEu_horiz1(i,j)=elem(ie)%accum%KEu_horiz1(i,j) &
                    -dp3d(i,j,k)*( &
                    elem(ie)%state%v(i,j,1,k,n0)*gradKE(i,j,1,k) + &
                    elem(ie)%state%v(i,j,2,k,n0)*gradKE(i,j,2,k) )
               elem(ie)%accum%KEu_horiz2(i,j)=elem(ie)%accum%KEu_horiz2(i,j)              &
                    -KE(i,j,k)*divdp(i,j,k)
               !  Form horiz advection of KE-w
               elem(ie)%accum%KEw_horiz1(i,j)=elem(ie)%accum%KEw_horiz1(i,j)-   &
                    dp3d(i,j,k) * (&
                    elem(ie)%state%w_i(i,j,k,n0) * v_gradw_i(i,j,k)    + &
                    elem(ie)%state%w_i(i,j,k+1,n0) * v_gradw_i(i,j,k+1) )/2
               elem(ie)%accum%KEw_horiz2(i,j)=elem(ie)%accum%KEw_horiz2(i,j)-   &
                    divdp(i,j,k)*(elem(ie)%state%w_i(i,j,k,n0)**2 + &
                    elem(ie)%state%w_i(i,j,k+1,n0)**2 ) /4
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
               elem(ie)%accum%KEw_vert1(i,j)=elem(ie)%accum%KEw_vert1(i,j) - &
                    dp3d(i,j,k) * &
                    (w_vadv_i(i,j,k)*elem(ie)%state%w_i(i,j,k,n0)+ &
                    w_vadv_i(i,j,k+1)*elem(ie)%state%w_i(i,j,k+1,n0))/2
               
               elem(ie)%accum%KEw_vert2(i,j)=elem(ie)%accum%KEw_vert2(i,j)      &
                    -d_eta_dot_dpdn_dn* &
                    (.5*elem(ie)%state%w_i(i,j,k,n0)**2 +&
                    .5*elem(ie)%state%w_i(i,j,k+1,n0)**2)/2
               
               !  Form IEvert1
               elem(ie)%accum%IEvert1(i,j)=elem(ie)%accum%IEvert1(i,j)      &
                    -exner(i,j,k)*theta_vadv(i,j,k)                        
               ! Form IEvert2 
               ! here use of dpnh_dp_i on boundry (with incorrect data)
               ! is harmess becuase eta_dot_dpdn=0
               elem(ie)%accum%IEvert2(i,j)=elem(ie)%accum%IEvert2(i,j)      &
                    + ( dpnh_dp_i(i,j,k)*eta_dot_dpdn(i,j,k)+ &
                        dpnh_dp_i(i,j,k+1)*eta_dot_dpdn(i,j,k+1)) &
                    *(phi_i(i,j,k+1)-phi_i(i,j,k))/2
               
               !  Form PEhoriz1
               elem(ie)%accum%PEhoriz1(i,j)=(elem(ie)%accum%PEhoriz1(i,j))  &
                    -phi(i,j,k)*divdp(i,j,k) 
               !  Form PEhoriz2
               elem(ie)%accum%PEhoriz2(i,j)=elem(ie)%accum%PEhoriz2(i,j)    &
                    -dp3d(i,j,k)* &
                    (elem(ie)%state%v(i,j,1,k,n0)*                          &
                    (gradphinh_i(i,j,1,k)+gradphinh_i(i,j,1,k+1))/2  +      &
                    elem(ie)%state%v(i,j,2,k,n0)*                           &
                    (gradphinh_i(i,j,2,k)+gradphinh_i(i,j,2,k+1))/2  )
               
               !  Form PEvert1
               elem(ie)%accum%PEvert1(i,j) = elem(ie)%accum%PEvert1(i,j)    &
                    -phi(i,j,k)*d_eta_dot_dpdn_dn                                 
               elem(ie)%accum%PEvert2(i,j) = elem(ie)%accum%PEvert2(i,j)     &
                    -dp3d(i,j,k)*(phi_vadv_i(i,j,k)+phi_vadv_i(i,j,k+1))/2
               
               !  Form T01
               elem(ie)%accum%T01(i,j)=elem(ie)%accum%T01(i,j)               &
                    -(Cp*elem(ie)%state%vtheta_dp(i,j,k,n0))                       &
                    *(gradexner(i,j,1,k)*elem(ie)%state%v(i,j,1,k,n0) +           &
                    gradexner(i,j,2,k)*elem(ie)%state%v(i,j,2,k,n0))              
               !  Form S1 
               elem(ie)%accum%S1(i,j)=elem(ie)%accum%S1(i,j)                 &
                    -Cp*exner(i,j,k)*div_v_theta(i,j,k)

               !  Form P1  = -P2  (no reason to compute P2?)
               elem(ie)%accum%P1(i,j)=elem(ie)%accum%P1(i,j) -g*dp3d(i,j,k)* &
                    ( elem(ie)%state%w_i(i,j,k,n0) + &
                    elem(ie)%state%w_i(i,j,k+1,n0) )/2
               !  Form P2
               elem(ie)%accum%P2(i,j)=elem(ie)%accum%P2(i,j) + g*dp3d(i,j,k)*&
                    ( elem(ie)%state%w_i(i,j,k,n0) + &
                    elem(ie)%state%w_i(i,j,k+1,n0) )/2
            enddo
         enddo
      enddo

      ! these terms are better easier to compute by summing interfaces
      do k=2,nlev
         elem(ie)%accum%T2(:,:)=elem(ie)%accum%T2(:,:)+                &
              (g*elem(ie)%state%w_i(:,:,k,n0)-v_gradphinh_i(:,:,k)) &
               * dpnh_dp_i(:,:,k)*dp3d_i(:,:,k)
      enddo
      ! boundary terms
      do k=1,nlevp,nlev
         elem(ie)%accum%T2(:,:)=elem(ie)%accum%T2(:,:)+                &
           (g*elem(ie)%state%w_i(:,:,k,n0)-v_gradphinh_i(:,:,k)) &
           * dpnh_dp_i(:,:,k)*dp3d_i(:,:,k)/2
      enddo
      ! boundary term is incorrect.  save the term so we can correct it
      ! once we have coorect value of dpnh_dp_i:
      elem(ie)%accum%T2_nlevp_term(:,:)=&
           (g*elem(ie)%state%w_i(:,:,nlevp,n0)-v_gradphinh_i(:,:,nlevp)) &
           * dp3d_i(:,:,nlevp)/2

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
        elem(ie)%state%w_i(:,:,k,np1)    = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%w_i(:,:,k,nm1)   &
          + dt2*w_tens(:,:,k))
        elem(ie)%state%vtheta_dp(:,:,k,np1) = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%vtheta_dp(:,:,k,nm1) &
          + dt2*theta_tens(:,:,k))
        elem(ie)%state%phinh_i(:,:,k,np1)   = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%phinh_i(:,:,k,nm1) & 
          + dt2*phi_tens(:,:,k))

        elem(ie)%state%dp3d(:,:,k,np1) = &
             elem(ie)%spheremp(:,:) * (scale3 * elem(ie)%state%dp3d(:,:,k,nm1) - &
             scale1*dt2 * (divdp(:,:,k) + eta_dot_dpdn(:,:,k+1)-eta_dot_dpdn(:,:,k)))
     enddo
     k=nlevp
     elem(ie)%state%w_i(:,:,k,np1)    = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%w_i(:,:,k,nm1)   &
          + dt2*w_tens(:,:,k))


     kptr=0
     call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,nlyr_tot)
     kptr=kptr+2*nlev
     call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr, nlyr_tot)
     kptr=kptr+nlev
     call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%vtheta_dp(:,:,:,np1),nlev,kptr,nlyr_tot)
     if (.not. theta_hydrostatic_mode) then
        kptr=kptr+nlev
        call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%w_i(:,:,:,np1),nlevp,kptr,nlyr_tot)
        kptr=kptr+nlevp
        call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%phinh_i(:,:,:,np1),nlev,kptr,nlyr_tot)
     endif

   end do ! end do for the ie=nets,nete loop

  call t_startf('caar_bexchV')
  call bndry_exchangeV(hybrid,edge_g)
  call t_stopf('caar_bexchV')

  do ie=nets,nete
     kptr=0
     call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,nlyr_tot)
     kptr=kptr+2*nlev
     call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr,nlyr_tot)
     kptr=kptr+nlev
     call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%vtheta_dp(:,:,:,np1),nlev,kptr,nlyr_tot)
     if (.not. theta_hydrostatic_mode) then
        kptr=kptr+nlev
        call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%w_i(:,:,:,np1),nlevp,kptr,nlyr_tot)
        kptr=kptr+nlevp
        call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%phinh_i(:,:,:,np1),nlev,kptr,nlyr_tot)
     endif
      
     ! ====================================================
     ! Scale tendencies by inverse mass matrix
     ! ====================================================
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        elem(ie)%state%dp3d(:,:,k,np1) =elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
        elem(ie)%state%vtheta_dp(:,:,k,np1)=elem(ie)%rspheremp(:,:)*elem(ie)%state%vtheta_dp(:,:,k,np1)
        elem(ie)%state%w_i(:,:,k,np1)    =elem(ie)%rspheremp(:,:)*elem(ie)%state%w_i(:,:,k,np1)
        elem(ie)%state%phinh_i(:,:,k,np1)=elem(ie)%rspheremp(:,:)*elem(ie)%state%phinh_i(:,:,k,np1)
        elem(ie)%state%v(:,:,1,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,1,k,np1)
        elem(ie)%state%v(:,:,2,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,2,k,np1)
     end do
     k=nlevp
     elem(ie)%state%w_i(:,:,k,np1)    =elem(ie)%rspheremp(:,:)*elem(ie)%state%w_i(:,:,k,np1)


     ! now we can compute the correct dphn_dp_i() at the surface:
     if (.not. theta_hydrostatic_mode) then
        ! solve for (dpnh_dp_i-1)
        dpnh_dp_i(:,:,nlevp) = 1 + &
             ((elem(ie)%state%v(:,:,1,nlev,np1)*elem(ie)%derived%gradphis(:,:,1) + &
             elem(ie)%state%v(:,:,2,nlev,np1)*elem(ie)%derived%gradphis(:,:,2))/g - &
             elem(ie)%state%w_i(:,:,nlevp,np1)) / &
             (g + ( elem(ie)%derived%gradphis(:,:,1)**2 + &
             elem(ie)%derived%gradphis(:,:,2)**2)/(2*g)) 
        
        ! update solution with new dpnh_dp_i value:
        elem(ie)%state%w_i(:,:,nlevp,np1) = elem(ie)%state%w_i(:,:,nlevp,np1) +&
             scale1*g*(dpnh_dp_i(:,:,nlevp)-1)
        elem(ie)%state%v(:,:,1,nlev,np1) =  elem(ie)%state%v(:,:,1,nlev,np1) -&
             scale1*(dpnh_dp_i(:,:,nlevp)-1)*elem(ie)%derived%gradphis(:,:,1)/2
        elem(ie)%state%v(:,:,2,nlev,np1) =  elem(ie)%state%v(:,:,2,nlev,np1) -&
             scale1*(dpnh_dp_i(:,:,nlevp)-1)*elem(ie)%derived%gradphis(:,:,2)/2
        

#ifdef ENERGY_DIAGNOSTICS
        ! add in boundary term to T2 and S2 diagnostics:
        if (compute_diagnostics) then
           elem(ie)%accum%T2(:,:)=elem(ie)%accum%T2(:,:)+                &
                elem(ie)%accum%T2_nlevp_term(:,:)*(dpnh_dp_i(:,:,nlevp)-1)
           elem(ie)%accum%S2(:,:)=-elem(ie)%accum%T2(:,:)      
        endif
#endif

        temp(:,:,1) =  (elem(ie)%state%v(:,:,1,nlev,np1)*elem(ie)%derived%gradphis(:,:,1) + &
             elem(ie)%state%v(:,:,2,nlev,np1)*elem(ie)%derived%gradphis(:,:,2))/g
        if ( maxval(abs(temp(:,:,1)-elem(ie)%state%w_i(:,:,nlevp,np1))) >1e-10) then
           write(iulog,*) 'WARNING: w(np1) surface b.c. violated'
           write(iulog,*) 'val1 = ',temp(:,:,1)
           write(iulog,*) 'val2 = ',elem(ie)%state%w_i(:,:,nlevp,np1)
           write(iulog,*) 'diff: ',temp(:,:,1)-elem(ie)%state%w_i(:,:,nlevp,np1)
        endif
     endif
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
  real (kind=real_kind), pointer, dimension(:,:,:)   :: vtheta_dp
  real (kind=real_kind), pointer, dimension(:,:)   :: phis
  real (kind=real_kind) :: JacD(nlev,np,np)  , JacL(nlev-1,np,np)
  real (kind=real_kind) :: JacU(nlev-1,np,np), JacU2(nlev-2,np,np)
  real (kind=real_kind) :: pnh(np,np,nlev)     ! nh (nonydro) pressure
  real (kind=real_kind) :: dp3d_i(np,np,nlevp)
  real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)
  real (kind=real_kind) :: exner(np,np,nlev)     ! exner nh pressure
  real (kind=real_kind) :: w_n0(np,np,nlevp)    
  real (kind=real_kind) :: phi_n0(np,np,nlevp)    
  real (kind=real_kind) :: Ipiv(nlev,np,np)
  real (kind=real_kind) :: Fn(np,np,nlev),x(nlev,np,np)
  real (kind=real_kind) :: itererr,itererrtemp(np,np)
  real (kind=real_kind) :: itercountmax,itererrmax
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
    w_n0 = elem(ie)%state%w_i(:,:,:,np1)
    phi_n0 = elem(ie)%state%phinh_i(:,:,:,np1)
    itercount=0

    ! approximate the initial error of f(x) \approx 0
    dp3d  => elem(ie)%state%dp3d(:,:,:,np1)
    vtheta_dp  => elem(ie)%state%vtheta_dp(:,:,:,np1)
    phi_np1 => elem(ie)%state%phinh_i(:,:,:,np1)
    phis => elem(ie)%state%phis(:,:)

    call get_pnh_and_exner(hvcoord,vtheta_dp,dp3d,phi_np1,pnh,exner,dpnh_dp_i)

    dp3d_i(:,:,1) = dp3d(:,:,1)
    dp3d_i(:,:,nlevp) = dp3d(:,:,nlev)
    do k=2,nlev
       dp3d_i(:,:,k)=(dp3d(:,:,k)+dp3d(:,:,k-1))/2
    end do

   ! we first compute the initial Jacobian J0 and residual r0 and their infinity norms
     Fn(:,:,1:nlev) = phi_np1(:,:,1:nlev)-phi_n0(:,:,1:nlev) &
       - dt2*g*w_n0(:,:,1:nlev) + (dt2*g)**2 * (1.0-dpnh_dp_i(:,:,1:nlev))

     norminfr0=0.d0
     norminfJ0=0.d0
      ! Here's how to call inexact Jacobian
!     call get_dirk_jacobian(Jac2L,Jac2D,Jac2U,dt2,dp3d,phi_np1,pnh,0,&
!       1d-6,hvcoord,dpnh_dp_i,vtheta_dp)
      ! here's the call to the exact Jacobian
     call get_dirk_jacobian(JacL,JacD,JacU,dt2,dp3d,phi_np1,pnh,1)

    ! compute dp3d-weighted infinity norms of the initial Jacobian and residual
#if (defined COLUMN_OPENMP)
!$omp parallel do private(i,j) collapse(2)
#endif
     do i=1,np
     do j=1,np
       itererrtemp(i,j)=0 
       do k=1,nlev
        norminfr0(i,j)=max(norminfr0(i,j),abs(Fn(i,j,k)) *dp3d_i(i,j,k))
        if (k.eq.1) then
          norminfJ0(i,j) = max(norminfJ0(i,j),(dp3d_i(i,j,k)*abs(JacD(k,i,j))+dp3d_i(i,j,k+1))*abs(JacU(k,i,j)))
        elseif (k.eq.nlev) then
          norminfJ0(i,j) = max(norminfJ0(i,j),(dp3d_i(i,j,k-1)*abs(JacL(k,i,j))+abs(JacD(k,i,j))*dp3d(i,j,k)))
        else
          norminfJ0(i,j) = max(norminfJ0(i,j),(dp3d_i(i,j,k-1)*abs(JacL(k,i,j))+dp3d_i(i,j,k)*abs(JacD(k,i,j))+ &
            dp3d_i(i,j,k+1)*abs(JacU(k,i,j))))
        end if
        itererrtemp(i,j)=itererrtemp(i,j)+Fn(i,j,k)**2.d0 *dp3d_i(i,j,k)
      end do
      itererrtemp(i,j)=sqrt(itererrtemp(i,j))
    end do
    end do

    maxnorminfJ0r0=max(maxval(norminfJ0(:,:)),maxval(norminfr0(:,:)))
    itererr=maxval(itererrtemp(:,:))/maxnorminfJ0r0


    do while ((itercount < maxiter).and.(itererr > itertol))

      info(:,:) = 0
      ! Here's how to call inexact Jacobian
!      call get_dirk_jacobian(JacL,JacD,JacU,dt2,dp3d,phi_np1,pnh,0,&
!       1d-4,hvcoord,dpnh_dp_i,vtheta_dp)
      ! here's the call to the exact Jacobian
       call get_dirk_jacobian(JacL,JacD,JacU,dt2,dp3d,phi_np1,pnh,1)

 
#if (defined COLUMN_OPENMP)
!$omp parallel do private(i,j) collapse(2)
#endif
      do i=1,np
      do j=1,np
        x(1:nlev,i,j) = -Fn(i,j,1:nlev)  !+Fn(i,j,nlev+1:2*nlev,1)/(g*dt2))
        call DGTTRF(nlev, JacL(:,i,j), JacD(:,i,j),JacU(:,i,j),JacU2(:,i,j), Ipiv(:,i,j), info(i,j) )
        ! Tridiagonal solve
        call DGTTRS( 'N', nlev,1, JacL(:,i,j), JacD(:,i,j), JacU(:,i,j), JacU2(:,i,j), Ipiv(:,i,j),x(:,i,j), nlev, info(i,j) )
        ! update approximate solution of phi
        phi_np1(i,j,1:nlev) = phi_np1(i,j,1:nlev) + x(1:nlev,i,j)
      end do
      end do

      call get_pnh_and_exner(hvcoord,vtheta_dp,dp3d,phi_np1,pnh,exner,dpnh_dp_i)

      ! update approximate solution of w
      elem(ie)%state%w_i(:,:,1:nlev,np1) = w_n0(:,:,1:nlev) - g*dt2 * &
        (1.0-dpnh_dp_i(:,:,1:nlev))
      ! update right-hand side of phi
      Fn(:,:,1:nlev) = phi_np1(:,:,1:nlev)-phi_n0(:,:,1:nlev) &
        - dt2*g*w_n0(:,:,1:nlev) + (dt2*g)**2 * (1.0-dpnh_dp_i(:,:,1:nlev))

      ! compute relative errors
      itererrtemp=0.d0
#if (defined COLUMN_OPENMP)
!$omp parallel do private(i,j) collapse(2)
#endif
      do i=1,np
      do j=1,np
        do k=1,nlev
          itererrtemp(i,j)=itererrtemp(i,j)+Fn(i,j,k)**2.d0 *dp3d_i(i,j,k)
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
  maxiter=itercount
  itertol=itererr
!  print *, 'max itercount', itercountmax, 'maxitererr ', itererrmax

  call t_stopf('compute_stage_value_dirk')

  end subroutine compute_stage_value_dirk




end module prim_advance_mod

