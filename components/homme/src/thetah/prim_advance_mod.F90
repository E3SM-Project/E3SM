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
  use element_ops,    only: get_p_hydrostatic, get_p_nonhydrostatic
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
  public :: prim_advance_exp, prim_advance_init, &
       applyCAMforcing_dynamics, applyCAMforcing

!  type (EdgeBuffer_t) :: edge5
  type (EdgeBuffer_t) :: edge6
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

  end subroutine prim_advance_init


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
          elem(ie)%state%theta(:,:,:,np1)= elem(ie)%state%theta(:,:,:,n0)/3 &
               + 2*elem(ie)%state%theta(:,:,:,np1)/3
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
          elem(ie)%state%theta(:,:,:,nm1)= (5*elem(ie)%state%theta(:,:,:,nm1) &
               - elem(ie)%state%theta(:,:,:,n0) )/4
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
          elem(ie)%accum%DIFFTHETA(:,:,:)=elem(ie)%state%theta(:,:,:,np1)
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
          elem(ie)%accum%DIFFTHETA(:,:,k)=( elem(ie)%state%theta(:,:,k,np1) -&
               elem(ie)%accum%DIFFTHETA(:,:,k) ) / dt_vis
         enddo
       enddo
    endif
#endif

    tevolve=tevolve+dt

    call t_stopf('prim_advance_exp')
    end subroutine prim_advance_exp




  subroutine applyCAMforcing(elem,hvcoord,np1,np1_qdp,dt,nets,nete)

  use physical_constants, only: Cp

  implicit none
  type (element_t),       intent(inout) :: elem(:)
  real (kind=real_kind),  intent(in)    :: dt
  type (hvcoord_t),       intent(in)    :: hvcoord
  integer,                intent(in)    :: np1,nets,nete,np1_qdp

  ! local
  integer :: i,j,k,ie,q
  real (kind=real_kind) :: v1,dp
  real (kind=real_kind) :: beta(np,np),E0(np,np),ED(np,np),dp0m1(np,np),dpsum(np,np)

  do ie=nets,nete
     ! apply forcing to Qdp
     elem(ie)%derived%FQps(:,:)=0
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

  enddo
  call applyCAMforcing_dynamics(elem,hvcoord,np1,np1_qdp,dt,nets,nete)

  end subroutine applyCAMforcing



  subroutine applyCAMforcing_dynamics(elem,hvcoord,np1,np1_qdp,dt,nets,nete)

  use hybvcoord_mod,  only: hvcoord_t

  implicit none
  type (element_t)     ,  intent(inout) :: elem(:)
  real (kind=real_kind),  intent(in)    :: dt
  type (hvcoord_t),       intent(in)    :: hvcoord
  integer,                intent(in)    :: np1,nets,nete,np1_qdp

  integer :: i,j,k,ie,q
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev)
  real (kind=real_kind) :: pnh_i(np,np,nlevp)

  do ie=nets,nete
     do k=1,nlev
        dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
             ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
     enddo
     if (theta_hydrostatic_mode) then
        call get_p_hydrostatic(pnh,pnh_i,exner,hvcoord,&
          dp,elem(ie)%state%Qdp(:,:,:,1,np1_qdp))
     else
        call get_p_nonhydrostatic(pnh,pnh_i,exner,hvcoord,&
             elem(ie)%state%theta(:,:,:,np1),&
             dp,elem(ie)%state%phi(:,:,:,np1),elem(ie)%state%phis,&
             elem(ie)%state%Qdp(:,:,:,1,np1_qdp))
     endif

     elem(ie)%state%theta(:,:,:,np1) = elem(ie)%state%theta(:,:,:,np1) + &
          dt*elem(ie)%derived%FT(:,:,:) / exner(:,:,:)
     elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + dt*elem(ie)%derived%FM(:,:,:,:)
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
  use physical_constants, only: Cp
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
  real (kind=real_kind) :: v1,v2,dt,heating
  real (kind=real_kind) :: temp(np,np,nlev)



  if (nu == 0 .and. nu_p==0 ) return;
  call t_startf('advance_hypervis')


  dt=dt2/hypervis_subcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  regular viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 1) then
     call abortmp( 'ERROR: hypervis_order == 1 not coded')
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
                 lap_s(:,:,2)=laplace_sphere_wk(elem(ie)%state%theta(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
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
           enddo


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
!                    heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )
                    heating = 0

                    elem(ie)%state%dp3d(i,j,k,nt)=elem(ie)%state%dp3d(i,j,k,nt) &
                         +stens(i,j,k,1,ie)

                    elem(ie)%state%theta(i,j,k,nt)=elem(ie)%state%theta(i,j,k,nt) &
                         +stens(i,j,k,2,ie)-heating/cp

                    elem(ie)%state%w(i,j,k,nt)=elem(ie)%state%w(i,j,k,nt) &
                         +stens(i,j,k,3,ie)

                    elem(ie)%state%phi(i,j,k,nt)=elem(ie)%state%phi(i,j,k,nt) &
                         +stens(i,j,k,4,ie)

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
  use prim_si_mod, only : preq_vertadv_v, preq_omega_ps, preq_hydrostatic, preq_hydrostatic_v2

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
  real (kind=real_kind), pointer, dimension(:,:,:)   :: theta

  real (kind=real_kind) :: omega_p(np,np,nlev)
  real (kind=real_kind) :: vort(np,np,nlev)      ! vorticity
  real (kind=real_kind) :: divdp(np,np,nlev)     
  real (kind=real_kind) :: pnh(np,np,nlev)     ! nh (nonydro) pressure
  real (kind=real_kind) :: dpnh(np,np,nlev)    ! 
  real (kind=real_kind) :: pnh_i(np,np,nlevp)  ! nh pressre on interfaces
  real (kind=real_kind) :: exner(np,np,nlev)     ! exner nh pressure
  real (kind=real_kind) :: dpnh_dp(np,np,nlev)   ! dpnh / dp3d  
  real (kind=real_kind) :: grad_exner(np,np,2,nlev)     
  real (kind=real_kind) :: eta_dot_dpdn(np,np,nlevp)  ! vertical velocity at interfaces

  real (kind=real_kind) :: vtens1(np,np,nlev)
  real (kind=real_kind) :: vtens2(np,np,nlev)
  real (kind=real_kind) :: v_vadv(np,np,2,nlev)   ! velocity vertical advection
  real (kind=real_kind) :: s_state(np,np,nlev,3)  ! scalars w,theta,phi
  real (kind=real_kind) :: s_vadv(np,np,nlev,3)   ! scalar vertical advection 
  real (kind=real_kind) :: stens(np,np,nlev,3)    ! tendencies w,theta,phi

  real (kind=real_kind) ::  temp(np,np,nlev)
  real (kind=real_kind), dimension(np,np)      :: sdot_sum   ! temporary field
  real (kind=real_kind), dimension(np,np,2)    :: vtemp     ! generic gradient storage
  real (kind=real_kind), dimension(np,np,2)    :: vtemp2    ! secondary generic gradient storage
  real (kind=real_kind), dimension(np,np,2,nlev):: vdp       !                            
  real (kind=real_kind), dimension(np,np)      :: KE
  real (kind=real_kind), dimension(np,np)      :: Eexner     ! energy diagnostic Exner pressure
  real (kind=real_kind), dimension(np,np,2)      :: thetau  !  theta*u in the diagnostics
  real (kind=real_kind), dimension(np,np)        :: divtemp ! temp divergence in the  in the diagnostics


  real (kind=real_kind) ::  v1,v2,glnps1,glnps2,gpterm
  integer :: i,j,k,kptr,ie

  call t_startf('compute_and_apply_rhs')
  do ie=nets,nete
     dp3d  => elem(ie)%state%dp3d(:,:,:,n0)
     theta  => elem(ie)%state%theta(:,:,:,n0)

     if (theta_hydrostatic_mode) then
        phi => elem(ie)%derived%phi(:,:,:)

        call get_p_hydrostatic(pnh,pnh_i,exner,hvcoord,&
             dp3d,elem(ie)%state%Qdp(:,:,:,1,qn0))
        
        ! Compute Hydrostatic equation, modeld after CCM-3
        do k=1,nlev
           !temp(:,:,k) = Cp*theta(:,:,k)*&
           !  ( (pnh_i(:,:,k+1)/p0)**kappa - (pnh_i(:,:,k)/p0)**kappa )
           temp(:,:,k) = dp3d(:,:,k) * ( Rgas*theta(:,:,k)*exner(:,:,k)/pnh(:,:,k))
        enddo
        call preq_hydrostatic_v2(phi,elem(ie)%state%phis,temp)
        dpnh_dp(:,:,:) = 1
        
        !    as a debug step, use nonhydrostatic formulas to compute exner pressure:
        !    needs 2x smaller timestep, results are ok, maybe noisy
        !call get_p_nonhydrostatic(pnh,dpnh,exner,hvcoord,theta,dp3d,&
        !     phi,elem(ie)%state%phis,elem(ie)%state%Qdp(:,:,:,1,qn0))
     else
        phi => elem(ie)%state%phi(:,:,:,n0)
        call get_p_nonhydrostatic(pnh,dpnh,exner,hvcoord,theta,dp3d,&
             phi,elem(ie)%state%phis,elem(ie)%state%Qdp(:,:,:,1,qn0))
        
        ! d(p-nh) / d(p-hyrdostatic)
        dpnh_dp(:,:,:) = dpnh(:,:,:)/dp3d(:,:,:)
#if 0
        if (hybrid%masterthread) then
           if (ie==1) then
              do k=2,nlev,4
                 write(*,"(i3,4f15.5)") k,(elem(ie)%state%phi(1,1,k,n0)),&
                      (elem(ie)%state%phi(:,:,k-1,n0)),&
                      minval(dpnh_dp(:,:,k)),maxval(dpnh_dp(:,:,k))
              enddo
           endif
        endif
#endif
     endif




#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,v1,v2,vtemp)
#endif
     do k=1,nlev
        grad_exner(:,:,:,k) = gradient_sphere(exner(:,:,k),deriv,elem(ie)%Dinv)        
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
        s_state(:,:,:,1)=elem(ie)%state%w(:,:,:,n0)
        s_state(:,:,:,2)=elem(ie)%state%theta(:,:,:,n0)
        s_state(:,:,:,3)=elem(ie)%state%phi(:,:,:,n0)

        call preq_vertadv_v(elem(ie)%state%v(:,:,:,:,n0),s_state,3,eta_dot_dpdn,dp3d,v_vadv,s_vadv)
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
!$omp parallel do private(k,i,j,v1,v2,KE,vtemp,vtemp2,gpterm,glnps1,glnps2)
#endif
     vertloop: do k=1,nlev

        ! ================================================
        ! w,theta,phi tendencies:
        ! ================================================
        vtemp(:,:,:)   = gradient_sphere(elem(ie)%state%w(:,:,k,n0),deriv,elem(ie)%Dinv)
        stens(:,:,k,1) = -s_vadv(:,:,k,1) &
             -elem(ie)%state%v(:,:,1,k,n0)*vtemp(:,:,1) &
             -elem(ie)%state%v(:,:,2,k,n0)*vtemp(:,:,2) &
             - g *(1-dpnh_dp(:,:,k) )

        vtemp(:,:,:)   = gradient_sphere(theta(:,:,k),deriv,elem(ie)%Dinv)
        stens(:,:,k,2) = -s_vadv(:,:,k,2) &
             -elem(ie)%state%v(:,:,1,k,n0)*vtemp(:,:,1) &
             -elem(ie)%state%v(:,:,2,k,n0)*vtemp(:,:,2)

        vtemp(:,:,:)   = gradient_sphere(elem(ie)%state%phi(:,:,k,n0),deriv,elem(ie)%Dinv)
        stens(:,:,k,3) = -s_vadv(:,:,k,3) &
             -elem(ie)%state%v(:,:,1,k,n0)*vtemp(:,:,1) &
             -elem(ie)%state%v(:,:,2,k,n0)*vtemp(:,:,2) &
             + g*elem(ie)%state%w(:,:,k,n0)                              


        ! vtemp = grad ( E + PHI )
        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              KE(i,j)=0.5D0*( v1*v1 + v2*v2 )
           end do
        end do
        vtemp = gradient_sphere(KE(:,:),deriv,elem(ie)%Dinv)
        vtemp2 = gradient_sphere(phi(:,:,k),deriv,elem(ie)%Dinv)
        vtemp2(:,:,1) = vtemp2(:,:,1)*dpnh_dp(:,:,k)
        vtemp2(:,:,2) = vtemp2(:,:,2)*dpnh_dp(:,:,k)

        do j=1,np
           do i=1,np
              glnps1 = cp*theta(i,j,k)*grad_exner(i,j,1,k)
              glnps2 = cp*theta(i,j,k)*grad_exner(i,j,2,k)
              !gpterm = T_v(i,j,k)/p(i,j,k)
              !glnps1 = Rgas*gpterm*grad_p(i,j,1,k)
              !glnps2 = Rgas*gpterm*grad_p(i,j,2,k)


              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)

              vtens1(i,j,k) =   - v_vadv(i,j,1,k)                           &
                   + v2*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - vtemp(i,j,1) -vtemp2(i,j,1) -glnps1

              vtens2(i,j,k) =   - v_vadv(i,j,2,k)                            &
                   - v1*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - vtemp(i,j,2) -vtemp2(i,j,2) -glnps2

           end do
        end do

     end do vertloop

#ifdef ENERGY_DIAGNOSTICSXXX
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
        elem(ie)%accum%KEhorz1=0
        elem(ie)%accum%KEhorz2=0
        elem(ie)%accum%IEhorz1=0
        elem(ie)%accum%IEhorz1_wet=0
        elem(ie)%accum%IEhorz2_wet=0
 	elem(ie)%accum%PEhorz1=0 
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
        ! The way this is set up desire that 
        !  KEhorz1 = 0 , KEvert1 =small(hopefully), KEhorz2+IEhorz1 = 0
        ! PEhorz1=0
        ! Essentially, we are checking that discrete integration by parts 
        ! holds for the sames terms as in the continuum equations
        !
        ! TODO:  ptop and phitop terms?
        !        

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
		 Eexner(i,j)=(p(i,j,k)/p0)**kappa                 
                 thetau(i,j,1)=theta(i,j,k)*v1
                 thetau(i,j,2)=theta(i,j,k)*v2
              enddo
           enddo
           vtemp = gradient_sphere(Ephi,deriv,elem(ie)%Dinv)
           vtemp2 = gradient_sphere(elem(ie)%state%w(:,:,:,n0),deriv,elem(ie)%Dinv)
           do j=1,np
              do i=1,np
	          v1     = elem(ie)%state%v(i,j,1,k,n0)
                  v2     = elem(ie)%state%v(i,j,2,k,n0)
                  w     = elem(ie)%state%w(i,j,k,n0)
               !  KEhorz1=-0.5*u^2 grad^T ( u dp/ds )-0.5*(dp/ds) grad(u^2)^T 
                  elem(ie)%accum%KEhorz1(i,j)= (v1*vtemp(i,j,1)+v2*vtemp(i,j,2))*dp3d(i,j,k)+Ephi(i,j)*divdp(i,j,k)
	       !  KEvert1= -0.5*w^2 grad^T( u dp/ds)-(dp/ds) u w grad(w)
                  KEvert1=0.5*w*w * divdp(i,j,k)+dp3d(i,j,k)*w*(v1*vtemp2(i,j,1)+v2*vtemp2(i,j,2))          
               ! KEvert1 is not guaranteed vanish in the hydrostatic discretization, we hope that it is small
              enddo
           enddo


           ! vtemp = grad_phi(:,:,k)
           vtemp  =gradient_sphere(phi(:,:,k),deriv,elem(ie)%Dinv)
           vtemp2 =gradient_sphere(Eexner(:,:),deriv,elem(ie)%Dinv)
           divtemp=divergence_sphere(thetau(:,:),deriv,elem(ie))
           do j=1,np
              do i=1,np
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 E = 0.5D0*( v1*v1 + v2*v2 )
               ! KEhorz2 = theta*grad(p^kappa)^T u + (dp/ds)*grad(phi)^T u
                 KEhorz2=theta(i,j,k)*(vtemp2(i,j,1)*v1+vtemp2(i,j,1)*v2)+dp3d(i,j,k)*(vtemp(i,j,1)*v1+vtemp(i,j,2)*v2) 
               ! Form second term of IEhorz1=-p^kappa grad^T(theta*u)-(dp/ds) u^T grad(phi)
                 IEhorz1=Eexner(i,j)*divtemp-dp3d(i,j,k)*(vtemp(i,j,1)*v1+vtemp(i,j,2)*v2)
               ! Form PEhorz1 = phi*grad^T ((dp/ds)*u) + (dp/ds)*u^T grad(phi), should equal zero
		 PEhorz1=phi(i,j,k)*divdp(i,j,k)+dp3d(i,j,k)*(vtemp(i,j,1)*v1+vtemp(i,j,2)*v2)	
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
        elem(ie)%state%theta(:,:,k,np1) = elem(ie)%spheremp(:,:)*(elem(ie)%state%theta(:,:,k,nm1) + dt2*stens(:,:,k,2))
        elem(ie)%state%phi(:,:,k,np1)   = elem(ie)%spheremp(:,:)*(elem(ie)%state%phi(:,:,k,nm1)   + dt2*stens(:,:,k,3))

        elem(ie)%state%dp3d(:,:,k,np1) = &
             elem(ie)%spheremp(:,:) * (elem(ie)%state%dp3d(:,:,k,nm1) - &
             dt2 * (divdp(:,:,k) + eta_dot_dpdn(:,:,k+1)-eta_dot_dpdn(:,:,k)))
        
     enddo

     kptr=0
     call edgeVpack(edge6, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr, ie)
     kptr=kptr+nlev
     call edgeVpack(edge6, elem(ie)%state%theta(:,:,:,np1),nlev,kptr,ie)
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
     call edgeVunpack(edge6, elem(ie)%state%theta(:,:,:,np1), nlev, kptr, ie)
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
        elem(ie)%state%theta(:,:,k,np1)=elem(ie)%rspheremp(:,:)*elem(ie)%state%theta(:,:,k,np1)
        elem(ie)%state%w(:,:,k,np1)    =elem(ie)%rspheremp(:,:)*elem(ie)%state%w(:,:,k,np1)
        elem(ie)%state%phi(:,:,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%phi(:,:,k,np1)
        elem(ie)%state%v(:,:,1,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,1,k,np1)
        elem(ie)%state%v(:,:,2,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,2,k,np1)
     end do
  end do

#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
  call t_stopf('compute_and_apply_rhs')

  end subroutine compute_and_apply_rhs



end module prim_advance_mod

