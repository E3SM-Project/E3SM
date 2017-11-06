
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
  use prim_advance_mod_base, only: applyCAMforcing_dynamics, applyCAMforcing, &
                                   vertical_mesh_init2, edge3p1, ur_weights, set_prescribed_wind
  implicit none
  public :: prim_advance_exp, prim_advance_init1, applyCAMforcing_dynamics, applyCAMforcing, vertical_mesh_init2

  real (kind=real_kind), allocatable :: p           (:,:,:,:)
  real (kind=real_kind), allocatable :: grad_p      (:,:,:,:,:)
  real (kind=real_kind), allocatable :: rdp         (:,:,:,:)
  real (kind=real_kind), allocatable :: vgrad_p     (:,:,:,:)
  real (kind=real_kind), allocatable :: vdp         (:,:,:,:,:)
  real (kind=real_kind), allocatable :: divdp       (:,:,:,:)
  real (kind=real_kind), allocatable :: vort        (:,:,:,:)
  real (kind=real_kind), allocatable :: kappa_star  (:,:,:,:)
  real (kind=real_kind), allocatable :: T_v         (:,:,:,:)
  real (kind=real_kind), allocatable :: phi         (:,:,:,:)
  real (kind=real_kind), allocatable :: omega_p     (:,:,:,:)
  real (kind=real_kind), allocatable :: eta_dot_dpdn(:,:,:,:)
  real (kind=real_kind), allocatable :: T_vadv      (:,:,:,:)
  real (kind=real_kind), allocatable :: v_vadv      (:,:,:,:,:)
  real (kind=real_kind), allocatable :: Ephi        (:,:,:,:)
  real (kind=real_kind), allocatable :: vtemp1      (:,:,:,:,:)
  real (kind=real_kind), allocatable :: vtemp2      (:,:,:,:,:)
  real (kind=real_kind), allocatable :: vtens       (:,:,:,:,:)
  real (kind=real_kind), allocatable :: ttens       (:,:,:,:)
  real (kind=real_kind), allocatable :: dptens      (:,:,:,:)
  real (kind=real_kind), allocatable :: lap_t       (:,:,:,:)
  real (kind=real_kind), allocatable :: lap_dp      (:,:,:,:)
  real (kind=real_kind), allocatable :: lap_v       (:,:,:,:,:)

contains



  subroutine prim_advance_init1(par, elem,integration)
    use edge_mod, only : initEdgeBuffer
    implicit none

    type (parallel_t) :: par
    type (element_t), intent(inout), target   :: elem(:)
    character(len=*), intent(in) :: integration
    integer :: i, ie

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

    allocate(p           (np,np  ,nlev  ,nelemd))
    allocate(grad_p      (np,np,2,nlev  ,nelemd))
    allocate(rdp         (np,np  ,nlev  ,nelemd))
    allocate(vgrad_p     (np,np  ,nlev  ,nelemd))
    allocate(vdp         (np,np,2,nlev  ,nelemd))
    allocate(divdp       (np,np  ,nlev  ,nelemd))
    allocate(vort        (np,np  ,nlev  ,nelemd))
    allocate(kappa_star  (np,np  ,nlev  ,nelemd))
    allocate(T_v         (np,np  ,nlev  ,nelemd))
    allocate(phi         (np,np  ,nlev  ,nelemd))
    allocate(omega_p     (np,np  ,nlev  ,nelemd))
    allocate(eta_dot_dpdn(np,np  ,nlev+1,nelemd))
    allocate(T_vadv      (np,np  ,nlev  ,nelemd))
    allocate(v_vadv      (np,np,2,nlev  ,nelemd))
    allocate(Ephi        (np,np  ,nlev  ,nelemd))
    allocate(vtemp1      (np,np,2,nlev  ,nelemd))
    allocate(vtemp2      (np,np,2,nlev  ,nelemd))
    allocate(vtens       (np,np,2,nlev  ,nelemd))
    allocate(ttens       (np,np  ,nlev  ,nelemd))
    allocate(dptens      (np,np  ,nlev  ,nelemd))
    allocate(lap_t       (np,np  ,nlev  ,nelemd))
    allocate(lap_dp      (np,np  ,nlev  ,nelemd))
    allocate(lap_v       (np,np,2,nlev  ,nelemd))

    !$acc enter data pcreate(edge3p1         )
    !$acc enter data pcreate(edge3p1%buf     )
    !$acc enter data pcreate(edge3p1%receive )
    !$acc enter data pcreate(edge3p1%putmap  )
    !$acc enter data pcreate(edge3p1%getmap  )
    !$acc enter data pcreate(edge3p1%reverse )
    !$acc enter data pcreate(edge3p1%tag     )
    !$acc enter data pcreate(edge3p1%srequest)
    !$acc enter data pcreate(edge3p1%rrequest)

  end subroutine prim_advance_init1



  subroutine prim_advance_exp(elem, deriv, hvcoord, hybrid,dt, tl,  nets, nete, compute_diagnostics)
    use bndry_mod,      only: bndry_exchangev
    use control_mod,    only: prescribed_wind, qsplit, tstep_type, rsplit, qsplit, integration
    use edge_mod,       only: edgevpack, edgevunpack, initEdgeBuffer
    use edgetype_mod,   only: EdgeBuffer_t
    use reduction_mod,  only: reductionbuffer_ordered_1d_t
    use time_mod,       only: timelevel_qdp
    use element_state,  only: state_v, state_T, state_dp3d
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
      !$acc parallel loop gang vector collapse(4)
      do ie = 1 , nelemd
        do k = 1 , nlev
          do j = 1 , np
            do i = 1 , np
              state_v   (i,j,:,k,nm1,ie)= (5*state_v   (i,j,:,k,nm1,ie) - state_v   (i,j,:,k,n0,ie) ) /4
              state_T   (i,j  ,k,nm1,ie)= (5*state_T   (i,j  ,k,nm1,ie) - state_T   (i,j  ,k,n0,ie) )/4
              state_dp3d(i,j  ,k,nm1,ie)= (5*state_dp3d(i,j  ,k,nm1,ie) - state_dp3d(i,j  ,k,n0,ie) )/4
            enddo
          enddo
        enddo
      enddo
      ! u5 = (5*u1/4 - u0/4) + 3dt/4 RHS(u4)
      call compute_and_apply_rhs(np1,nm1,np1,qn0,3*dt/4,elem,hvcoord,hybrid,&
      deriv,nets,nete,.false.,3*eta_ave_w/4)
      ! final method is the same as:
      ! u5 = u0 +  dt/4 RHS(u0)) + 3dt/4 RHS(u4)
      call t_stopf("U3-5stage_timestep")
    else if ((method==11).or.(method==12)) then
      ! Fully implicit JFNK method (vertically langragian not active yet)
      if (rsplit > 0) then
        call abortmp('ERROR: full_imp integration not yet coded for vert lagrangian adv option')
      endif
    else
      call abortmp('ERROR: bad choice of tstep_type')
    endif

    ! ==============================================
    ! Time-split Horizontal diffusion: nu.del^2 or nu.del^4
    ! U(*) = U(t+1)  + dt2 * HYPER_DIFF_TERM(t+1)
    ! ==============================================

    ! note:time step computes u(t+1)= u(t*) + RHS.
    ! for consistency, dt_vis = t-1 - t*, so this is timestep method dependent
    if (method<=10) then ! not implicit
      ! forward-in-time, hypervis applied to dp3d
      call advance_hypervis_dp(edge3p1,elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)
    endif

    call t_stopf('prim_advance_exp')
  end subroutine prim_advance_exp





    subroutine advance_hypervis_dp(edge3,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
      !  take one timestep of:
      !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
      !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
      !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
      use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis, swest
      use hybvcoord_mod, only : hvcoord_t
      use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk, laplace_sphere_wk_openacc, vlaplace_sphere_wk_openacc
      use edge_mod, only : edgevpack, edgevunpack, edgeVpack_openacc, edgeVunpack_openacc
      use edgetype_mod, only : EdgeBuffer_t, EdgeDescriptor_t
      use bndry_mod, only : bndry_exchangev
      use viscosity_mod, only : biharmonic_wk_dp3d, biharmonic_wk_dp3d_openacc
      use physical_constants, only: Cp
      use element_state, only: state_T, state_v, state_dp3d, timelevels
      implicit none
      type (hybrid_t)      , intent(in   ) :: hybrid
      type (element_t)     , intent(inout), target :: elem(:)
      type (EdgeBuffer_t)  , intent(inout) :: edge3
      type (derivative_t)  , intent(in   ) :: deriv
      type (hvcoord_t), intent(in)         :: hvcoord
      real (kind=real_kind) :: dt2
      integer :: nets,nete
      ! local
      real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
      real (kind=real_kind) :: nu_scale_top
      integer :: k,kptr,i,j,ie,ic,nt
      real (kind=real_kind) :: dpdn
      real (kind=real_kind) :: v1,v2,dt,heating

      if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;

      call t_startf('advance_hypervis_dp')
      dt=dt2/hypervis_subcycle
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  regular viscosity
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (hypervis_order == 1) then
        !  if (nu_p>0) call abortmp( 'ERROR: hypervis_order == 1 not coded for nu_p>0')
        !  do ic=1,hypervis_subcycle
        !     do ie=nets,nete
        !        do k=1,nlev
        !           lap_t=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
        !           lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
        !           ! advace in time.  (note: DSS commutes with time stepping, so we
        !           ! can time advance and then DSS.  this has the advantage of
        !           ! not letting any discontinuties accumulate in p,v via roundoff
        !           do j=1,np
        !              do i=1,np
        !                 elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)*elem(ie)%spheremp(i,j)  +  dt*nu_s*lap_t(i,j)
        !                 elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,1)
        !                 elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,2)
        !              enddo
        !           enddo
        !        enddo
        !        kptr=0
        !        call edgeVpack(edge3, elem(ie)%state%T(:,:,:,nt),nlev,kptr,ie)
        !        kptr=nlev
        !        call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,ie)
        !     enddo
        !
        !     call t_startf('ahdp_bexchV1')
        !     call bndry_exchangeV(hybrid,edge3)
        !     call t_stopf('ahdp_bexchV1')
        !
        !     do ie=nets,nete
        !        kptr=0
        !        call edgeVunpack(edge3, elem(ie)%state%T(:,:,:,nt), nlev, kptr, ie)
        !        kptr=nlev
        !        call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, ie)
        !        do k=1,nlev
        !           do j=1,np
        !              do i=1,np
        !                 elem(ie)%state%T(i,j,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,nt)
        !                 elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
        !                 elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
        !              enddo
        !           enddo
        !        enddo
        !     enddo
        !  enddo  ! subcycle
        write(*,*) 'Not implemented for OpenACC yet'
        write(*,*) __FILE__, ': ', __LINE__
        stop
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
      if (hypervis_order == 2) then
        do ic=1,hypervis_subcycle
          call biharmonic_wk_dp3d_openacc(elem,dptens,ttens,vtens,deriv,edge3,hybrid,vtemp1,vort,divdp,nt,1,nelemd)
          !$omp barrier
          !$omp master
          if (nu_p>0) then
            ! comptue mean flux
            !$acc parallel loop gang vector collapse(4)
            do ie = 1 , nelemd
              do k = 1 , nlev
                do j = 1 , np
                  do i = 1 , np
                    elem(ie)%derived%dpdiss_ave       (i,j,k)=elem(ie)%derived%dpdiss_ave       (i,j,k)+eta_ave_w*state_dp3d(i,j,k,nt,ie)/hypervis_subcycle
                    elem(ie)%derived%dpdiss_biharmonic(i,j,k)=elem(ie)%derived%dpdiss_biharmonic(i,j,k)+eta_ave_w*dptens    (i,j,k   ,ie)/hypervis_subcycle
                  enddo
                enddo
              enddo
            enddo
          endif
          if (nu_top>0) then
            ! advace in time.
            ! note: DSS commutes with time stepping, so we can time advance and then DSS.
            ! note: weak operators alreayd have mass matrix "included"
            ! add regular diffusion in top 3 layers:
            call laplace_sphere_wk_openacc (state_T   ,vtemp1    ,deriv,elem,.false.,lap_t ,nlev,1,nelemd,timelevels,nt,1,1  ,klim_in=3)
            call laplace_sphere_wk_openacc (state_dp3d,vtemp1    ,deriv,elem,.false.,lap_dp,nlev,1,nelemd,timelevels,nt,1,1  ,klim_in=3)
            call vlaplace_sphere_wk_openacc(state_v   ,vort,divdp,deriv,elem,.false.,nlev,1,nelemd,timelevels,nt,1,1,lap_v,1._real_kind,klim_in=3)
          endif
          !$acc parallel loop gang vector collapse(4)
          do ie=1,nelemd
            do k=1,nlev
              do j = 1 , np
                do i = 1 , np
                  nu_scale_top = 1
                  if (k==1) nu_scale_top=4
                  if (k==2) nu_scale_top=2
                  ! biharmonic terms need a negative sign:
                  if (nu_top>0 .and. k<=3) then
                    vtens(i,j,:,k,ie)=(-nu*vtens(i,j,:,k,ie) + nu_scale_top*nu_top*lap_v(i,j,:,k,ie))
                    ttens(i,j,k,ie)  =(-nu_s*ttens(i,j,k,ie) + nu_scale_top*nu_top*lap_t(i,j,k,ie) )
                    dptens(i,j,k,ie) =(-nu_p*dptens(i,j,k,ie) + nu_scale_top*nu_top*lap_dp(i,j,k,ie) )
                  else
                    vtens(i,j,:,k,ie)=-nu*vtens(i,j,:,k,ie)
                    ttens(i,j,k,ie)  =-nu_s*ttens(i,j,k,ie)
                    dptens(i,j,k,ie) =-nu_p*dptens(i,j,k,ie)
                  endif
                  if (nu_p==0) then
                    ! nu_p==0 is only for certain regression tests, so perfromance is not an issue
                    ! normalize so as to conserve IE
                    ! scale by 1/rho (normalized to be O(1))
                    ! dp/dn = O(ps0)*O(delta_eta) = O(ps0)/O(nlev)
                    dpdn = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,nt)
                    ttens(i,j,k,ie) = ttens(i,j,k,ie) * hvcoord%dp0(k)/dpdn
                    dptens(i,j,k,ie) = 0
                  endif
                  ! NOTE: we will DSS all tendicies, EXCEPT for dp3d, where we DSS the new state
                  elem(ie)%state%dp3d(i,j,k,nt) = elem(ie)%state%dp3d(i,j,k,nt)*elem(ie)%spheremp(i,j) + dt*dptens(i,j,k,ie)
                enddo
              enddo
            enddo
          enddo
          kptr = 0     ; call edgeVpack_openacc(edge3,ttens     ,  nlev,kptr,elem,1,nelemd,1,1)
          kptr = nlev  ; call edgeVpack_openacc(edge3,vtens     ,2*nlev,kptr,elem,1,nelemd,1,1)
          kptr = 3*nlev; call edgeVpack_openacc(edge3,state_dp3d,  nlev,kptr,elem,1,nelemd,timelevels,nt)
          !$omp end master
          !$omp barrier
          ! do ie=nets,nete
          !   kptr=0
          !   call edgeVpack(edge3, ttens(:,:,:,ie),nlev,kptr,ie)
          !   kptr=nlev
          !   call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,ie)
          !   kptr=3*nlev
          !   call edgeVpack(edge3,elem(ie)%state%dp3d(:,:,:,nt),nlev,kptr,ie)
          ! enddo

          call t_startf('ahdp_bexchV2')
          call bndry_exchangeV(hybrid,edge3)
          call t_stopf('ahdp_bexchV2')

          !$omp barrier
          !$omp master
          kptr = 0     ; call edgeVunpack_openacc(edge3,ttens     ,  nlev,kptr,elem,1,nelemd,1,1)
          kptr = nlev  ; call edgeVunpack_openacc(edge3,vtens     ,2*nlev,kptr,elem,1,nelemd,1,1)
          kptr = 3*nlev; call edgeVunpack_openacc(edge3,state_dp3d,  nlev,kptr,elem,1,nelemd,timelevels,nt)
          !$acc parallel loop gang vector collapse(4)
          do ie = 1 , nelemd
            do k=1,nlev
              do j = 1 , np
                do i = 1 , np
                  ! apply inverse mass matrix, accumulate tendencies
                  vtens     (i,j,1,k   ,ie) = dt*vtens  (i,j,1,k   ,ie)*elem(ie)%rspheremp(i,j)
                  vtens     (i,j,2,k   ,ie) = dt*vtens  (i,j,2,k   ,ie)*elem(ie)%rspheremp(i,j)
                  ttens     (i,j  ,k   ,ie) = dt*ttens  (i,j  ,k   ,ie)*elem(ie)%rspheremp(i,j)
                  state_dp3d(i,j  ,k,nt,ie) = state_dp3d(i,j  ,k,nt,ie)*elem(ie)%rspheremp(i,j)
                  ! apply hypervis to u -> u+utens:
                  ! E0 = dpdn * .5*u dot u + dpdn * T  + dpdn*PHIS
                  ! E1 = dpdn * .5*(u+utens) dot (u+utens) + dpdn * (T-X) + dpdn*PHIS
                  ! E1-E0:   dpdn (u dot utens) + dpdn .5 utens dot utens   - dpdn X
                  !      X = (u dot utens) + .5 utens dot utens
                  !  alt:  (u+utens) dot utens
                  ! update v first (gives better results than updating v after heating)
                  elem(ie)%state%v(i,j,:,k,nt)=elem(ie)%state%v(i,j,:,k,nt) + vtens(i,j,:,k,ie)
                  v1=elem(ie)%state%v(i,j,1,k,nt)
                  v2=elem(ie)%state%v(i,j,2,k,nt)
                  heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )
                  elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)+ttens(i,j,k,ie)-heating/cp
                enddo
              enddo
            enddo
          enddo
          !$omp end master
          !$omp barrier
        enddo
      endif

      call t_stopf('advance_hypervis_dp')
    end subroutine advance_hypervis_dp





  !
  ! phl notes: output is stored in first argument. Advances from 2nd argument using tendencies evaluated at 3rd rgument:
  ! phl: for offline winds use time at 3rd argument (same as rhs currently)
  !
  subroutine compute_and_apply_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,deriv,nets,nete,compute_diagnostics,eta_ave_w)
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
    ! ===================================
    use kinds,          only : real_kind
    use derivative_mod, only : derivative_t, divergence_sphere, divergence_sphere_openacc, gradient_sphere, &
                               gradient_sphere_openacc, vorticity_sphere, vorticity_sphere_openacc
    use derivative_mod, only : subcell_div_fluxes, subcell_dss_fluxes
    use edge_mod,       only : edgevpack, edgevunpack, edgeVpack_openacc, edgeVunpack_openacc
    use edgetype_mod,   only : edgedescriptor_t
    !use bndry_mod,      only : bndry_exchangev
    use bndry_mod    , only: bndry_exchangeV => bndry_exchangeV_simple_overlap
    use control_mod,    only : moisture, qsplit, use_cpstar, rsplit, swest
    use hybvcoord_mod,  only : hvcoord_t
    use physical_constants, only : cp, cpwater_vapor, Rgas, kappa
    use physics_mod,    only : virtual_specific_heat, virtual_temperature
    use prim_si_mod,    only : preq_vertadv, preq_omega_ps, preq_hydrostatic, preq_vertadv_openacc, preq_omega_ps_openacc, preq_hydrostatic_openacc
    use viscosity_base, only: smooth_phis
    use element_state, only: state_v, timelevels, state_dp3d, state_phis, state_t
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
    !real (kind=real_kind), dimension(np,np)        :: sdot_sum   ! temporary field
    !real (kind=real_kind), dimension(np,np,2     ) :: v         !
    !real (kind=real_kind), dimension(np,np,2,nlev) :: grad_p_m_pmet  ! gradient(p - p_met)
    !real (kind=real_kind), dimension(np,np,nlev+1) :: ph               ! half level pressures on p-grid
    real (kind=real_kind) :: vgrad_T, vtens1, vtens2, ttens, cp2,cp_ratio,E,de,Qt,v1,v2, glnps1,glnps2,gpterm, fcor_vort
    integer :: i,j,k,kptr,ie
    type (EdgeDescriptor_t) :: desc

    call t_startf('compute_and_apply_rhs')
    !$omp barrier
    !$omp master

    !$acc parallel loop gang vector collapse(3)
    do ie=1,nelemd
      ! dont thread this because of k-1 dependence:
      do j = 1 , np
        do i = 1 , np
          p(i,j,1,ie)=hvcoord%hyai(1)*hvcoord%ps0 + elem(ie)%state%dp3d(i,j,1,n0)/2
          !$acc loop seq
          do k=2,nlev
            p(i,j,k,ie)=p(i,j,k-1,ie) + elem(ie)%state%dp3d(i,j,k-1,n0)/2 + elem(ie)%state%dp3d(i,j,k,n0)/2
          enddo
        enddo
      enddo
    enddo
    call gradient_sphere_openacc(p,deriv,elem,grad_p,nlev,1,nelemd,1,1,1,1)
    ! ============================
    ! compute vgrad_lnps
    ! ============================
    !$acc parallel loop gang vector collapse(4)
    do ie=1,nelemd
      do k=1,nlev
        do j=1,np
          do i=1,np
            rdp(i,j,k,ie) = 1.0D0/elem(ie)%state%dp3d(i,j,k,n0)
            v1 = elem(ie)%state%v(i,j,1,k,n0)
            v2 = elem(ie)%state%v(i,j,2,k,n0)
            vgrad_p(i,j,k,ie) = (v1*grad_p(i,j,1,k,ie) + v2*grad_p(i,j,2,k,ie))
            vdp(i,j,1,k,ie) = v1*elem(ie)%state%dp3d(i,j,k,n0)
            vdp(i,j,2,k,ie) = v2*elem(ie)%state%dp3d(i,j,k,n0)
            ! Accumulate mean Vel_rho flux in vn0
            elem(ie)%derived%vn0(i,j,:,k)=elem(ie)%derived%vn0(i,j,:,k)+eta_ave_w*vdp(i,j,:,k,ie)
          enddo
        enddo
      enddo
    enddo
    ! =========================================
    ! Compute relative vorticity and divergence
    ! =========================================
    call divergence_sphere_openacc(vdp,deriv,elem,divdp,nlev,1,nelemd,1,1,1,1)
    call vorticity_sphere_openacc(state_v,deriv,elem,vort,nlev,nets,nete,timelevels,n0,1,1)

    !$acc parallel loop gang vector collapse(4) private(Qt)
    do ie=1,nelemd
      ! compute T_v for timelevel n0
      do k=1,nlev
        do j=1,np
          do i=1,np
            if (.not. use_moisture ) then
              T_v(i,j,k,ie) = elem(ie)%state%T(i,j,k,n0)
              kappa_star(i,j,k,ie) = kappa
            else
              Qt = elem(ie)%state%Qdp(i,j,k,1,qn0)/elem(ie)%state%dp3d(i,j,k,n0)
              T_v(i,j,k,ie) = Virtual_Temperature(elem(ie)%state%T(i,j,k,n0),Qt)
              if (use_cpstar==1) then
                kappa_star(i,j,k,ie) =  Rgas/Virtual_Specific_Heat(Qt)
              else
                kappa_star(i,j,k,ie) = kappa
              endif
            endif
          enddo
        enddo
      enddo
    enddo
    ! ====================================================
    ! Compute Hydrostatic equation, modeld after CCM-3
    ! ====================================================
    call preq_hydrostatic_openacc(phi,state_phis,T_v,p,state_dp3d,1,nelemd,timelevels,n0)
    ! ====================================================
    ! Compute omega_p according to CCM-3
    ! ====================================================
    call preq_omega_ps_openacc(omega_p,hvcoord,p,vgrad_p,divdp,1,nelemd)
    ! ==================================================
    ! Compute eta_dot_dpdn
    ! save sdot_sum as this is the -RHS of ps_v equation
    ! ==================================================
    if (rsplit>0) then
      ! VERTICALLY LAGRANGIAN:   no vertical motion
      !$acc parallel loop gang vector collapse(4)
      do ie=1,nelemd
        do k = 1 , nlev+1
          do j = 1 , np
            do i = 1 , np
              eta_dot_dpdn(i,j,k  ,ie)=0
              if (k <= nlev) then
                T_vadv    (i,j,k  ,ie)=0
                v_vadv    (i,j,k,:,ie)=0
              endif
            enddo
          enddo
        enddo
      enddo
    else
      write(*,*) 'Not currently supported in OpenACC port'
      write(*,*) __FILE__, ': ', __LINE__
      stop
      ! do ie=1,nelemd
      !   ! ==================================================
      !   ! zero partial sum for accumulating sum
      !   !    (div(v_k) + v_k.grad(lnps))*dsigma_k = div( v dp )
      !   ! used by eta_dot_dpdn and lnps tendency
      !   ! ==================================================
      !   sdot_sum=0
      !   do k=1,nlev
      !     ! ==================================================
      !     ! add this term to PS equation so we exactly conserve dry mass
      !     ! ==================================================
      !     sdot_sum(:,:) = sdot_sum(:,:) + divdp(:,:,k,ie)
      !     eta_dot_dpdn(:,:,k+1,ie) = sdot_sum(:,:)
      !   enddo
      !   ! ===========================================================
      !   ! at this point, eta_dot_dpdn contains integral_etatop^eta[ divdp ]
      !   ! compute at interfaces:
      !   !    eta_dot_dpdn = -dp/dt - integral_etatop^eta[ divdp ]
      !   ! for reference: at mid layers we have:
      !   !    omega = v grad p  - integral_etatop^eta[ divdp ]
      !   ! ===========================================================
      !   do k=1,nlev-1
      !     eta_dot_dpdn(:,:,k+1,ie) = hvcoord%hybi(k+1)*sdot_sum(:,:) - eta_dot_dpdn(:,:,k+1,ie)
      !   enddo
      !   eta_dot_dpdn(:,:,1     ,ie) = 0.0D0
      !   eta_dot_dpdn(:,:,nlev+1,ie) = 0.0D0
      ! enddo
      ! ! ===========================================================
      ! ! Compute vertical advection of T and v from eq. CCM2 (3.b.1)
      ! ! ==============================================
      ! do ie=1,nelemd
      !   call preq_vertadv(elem(ie)%state%T(:,:,:,n0),elem(ie)%state%v(:,:,:,:,n0),eta_dot_dpdn(:,:,:,ie),rdp(:,:,:,ie),T_vadv(:,:,:,ie),v_vadv(:,:,:,:,ie))
      ! enddo
    endif
    !$acc parallel loop gang vector collapse(4) private(v1,v2,E)
    do ie=1,nelemd
      do k=1,nlev  !  Loop index added (AAM)
        do j = 1 , np
          do i = 1 , np
            ! ================================
            ! accumulate mean vertical flux:
            ! ================================
            elem(ie)%derived%eta_dot_dpdn(i,j,k) = elem(ie)%derived%eta_dot_dpdn(i,j,k) + eta_ave_w*eta_dot_dpdn(i,j,k,ie)
            elem(ie)%derived%omega_p(i,j,k) = elem(ie)%derived%omega_p(i,j,k) + eta_ave_w*omega_p(i,j,k,ie)
            if (k == 1) then
              elem(ie)%derived%eta_dot_dpdn(i,j,nlev+1) = elem(ie)%derived%eta_dot_dpdn(i,j,nlev+1) + eta_ave_w*eta_dot_dpdn(i,j,nlev+1,ie)
            endif
            ! ==============================================
            ! Compute phi + kinetic energy term: 10*nv*nv Flops
            ! ==============================================
            v1 = elem(ie)%state%v(i,j,1,k,n0)
            v2 = elem(ie)%state%v(i,j,2,k,n0)
            E = 0.5D0*( v1*v1 + v2*v2 )
            Ephi(i,j,k,ie)=E+phi(i,j,k,ie)
          enddo
        enddo
      enddo
    enddo
    ! ================================================
    ! compute gradp term (ps/p)*(dp/dps)*T
    ! ================================================
    call gradient_sphere_openacc(state_t,deriv,elem,vtemp1,nlev,1,nelemd,timelevels,n0,1,1)
    call gradient_sphere_openacc(Ephi   ,deriv,elem,vtemp2,nlev,1,nelemd,1         ,1 ,1,1)
    !$acc parallel loop gang vector collapse(4)
    do ie=1,nelemd
      do k=1,nlev
        do j=1,np
          do i=1,np
            v1        = elem(ie)%state%v(i,j,1,k,n0)
            v2        = elem(ie)%state%v(i,j,2,k,n0)
            vgrad_T   = v1*vtemp1(i,j,1,k,ie) + v2*vtemp1(i,j,2,k,ie)
            gpterm    = T_v(i,j,k,ie)/p(i,j,k,ie)
            glnps1    = Rgas*gpterm*grad_p(i,j,1,k,ie)
            glnps2    = Rgas*gpterm*grad_p(i,j,2,k,ie)
            fcor_vort = elem(ie)%fcor(i,j) + vort(i,j,k,ie)
            vtens1    = - v_vadv(i,j,1,k,ie) + v2*fcor_vort - vtemp2(i,j,1,k,ie) - glnps1
            vtens2    = - v_vadv(i,j,2,k,ie) - v1*fcor_vort - vtemp2(i,j,2,k,ie) - glnps2
            ttens     = - T_vadv(i,j  ,k,ie) - vgrad_T + kappa_star(i,j,k,ie)*T_v(i,j,k,ie)*omega_p(i,j,k,ie)
            ! =========================================================
            ! local element timestep, store in np1.
            ! note that we allow np1=n0 or nm1
            ! apply mass matrix
            ! =========================================================
            elem(ie)%state%v   (i,j,1,k,np1) = elem(ie)%spheremp(i,j)*( elem(ie)%state%v   (i,j,1,k,nm1) + dt2*vtens1 )
            elem(ie)%state%v   (i,j,2,k,np1) = elem(ie)%spheremp(i,j)*( elem(ie)%state%v   (i,j,2,k,nm1) + dt2*vtens2 )
            elem(ie)%state%T   (i,j  ,k,np1) = elem(ie)%spheremp(i,j)*( elem(ie)%state%T   (i,j  ,k,nm1) + dt2*ttens  )
            elem(ie)%state%dp3d(i,j  ,k,np1) = elem(ie)%spheremp(i,j)*( elem(ie)%state%dp3d(i,j  ,k,nm1) - dt2* &
                                                                        (divdp(i,j,k,ie) + eta_dot_dpdn(i,j,k+1,ie)-eta_dot_dpdn(i,j,k,ie)) &
                                                                                                                      )
          enddo
        enddo
      enddo
    enddo

    ! =========================================================
    ! Pack ps(np1), T, and v tendencies into comm buffer
    ! =========================================================
    kptr = 0            ; call edgeVpack_openacc(edge3p1,state_T   ,nlev  ,kptr,elem,1,nelemd,timelevels,np1)
    kptr = kptr + nlev  ; call edgeVpack_openacc(edge3p1,state_v   ,nlev*2,kptr,elem,1,nelemd,timelevels,np1)
    kptr = kptr + 2*nlev; call edgeVpack_openacc(edge3p1,state_dp3d,nlev  ,kptr,elem,1,nelemd,timelevels,np1)
    !$omp end master
    !$omp barrier

    ! =============================================================
    ! Insert communications here: for shared memory, just a single
    ! sync is required
    ! =============================================================

    call t_startf('caar_bexchV')
    call bndry_exchangeV(hybrid,edge3p1)
    call t_stopf('caar_bexchV')

    !$omp barrier
    !$omp master
    ! ===========================================================
    ! Unpack the edges for vgrad_T and v tendencies...
    ! ===========================================================
    kptr = 0            ; call edgeVunpack_openacc(edge3p1,state_T   ,nlev  ,kptr,elem,1,nelemd,timelevels,np1)
    kptr = kptr + nlev  ; call edgeVunpack_openacc(edge3p1,state_v   ,nlev*2,kptr,elem,1,nelemd,timelevels,np1)
    kptr = kptr + 2*nlev; call edgeVunpack_openacc(edge3p1,state_dp3d,nlev  ,kptr,elem,1,nelemd,timelevels,np1)
    ! ====================================================
    ! Scale tendencies by inverse mass matrix
    ! ====================================================
    !$acc parallel loop gang vector collapse(4)
    do ie=1,nelemd
      do k=1,nlev
        do j = 1 , np
          do i = 1 , np
            elem(ie)%state%T   (i,j  ,k,np1) = elem(ie)%rspheremp(i,j)*elem(ie)%state%T   (i,j  ,k,np1)
            elem(ie)%state%v   (i,j,1,k,np1) = elem(ie)%rspheremp(i,j)*elem(ie)%state%v   (i,j,1,k,np1)
            elem(ie)%state%v   (i,j,2,k,np1) = elem(ie)%rspheremp(i,j)*elem(ie)%state%v   (i,j,2,k,np1)
            elem(ie)%state%dp3d(i,j  ,k,np1) = elem(ie)%rspheremp(i,j)*elem(ie)%state%dp3d(i,j  ,k,np1)
          enddo
        enddo
      enddo
    enddo
    !$omp end master
    !$omp barrier
    call t_stopf('compute_and_apply_rhs')

  end subroutine compute_and_apply_rhs

end module prim_advance_mod
