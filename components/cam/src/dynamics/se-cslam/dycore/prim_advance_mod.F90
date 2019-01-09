module prim_advance_mod
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use edgetype_mod,   only: EdgeBuffer_t
  use perf_mod,       only: t_startf, t_stopf, t_adj_detailf !, t_barrierf _EXTERNAL
  use cam_abortutils, only: endrun
  use parallel_mod,   only: parallel_t, HME_BNDRY_P2P!,HME_BNDRY_A2A
  use control_mod,    only: se_prescribed_wind_2d
  use thread_mod ,    only: horz_num_threads, vert_num_threads

  implicit none
  private
  save
  
  public :: prim_advance_exp, prim_advance_init, applyCAMforcing, calc_tot_energy_dynamics, compute_omega

  type (EdgeBuffer_t) :: edge3p1
  real (kind=r8), allocatable :: ur_weights(:)
  
contains
  
  subroutine prim_advance_init(par, elem)
    use edge_mod,       only: initEdgeBuffer
    use element_mod,    only: element_t
    use dimensions_mod, only: nlev,hypervis_on_plevs
    use control_mod,    only: qsplit,nu_p
    
    type (parallel_t)                       :: par
    type (element_t), target, intent(inout) :: elem(:)
    integer                                 :: i

    if (hypervis_on_plevs.and.nu_p>0) then
      call initEdgeBuffer(par,edge3p1,elem,5*nlev+1,bndry_type=HME_BNDRY_P2P, nthreads=horz_num_threads)
    else
      call initEdgeBuffer(par,edge3p1,elem,4*nlev+1,bndry_type=HME_BNDRY_P2P, nthreads=horz_num_threads)
    end if
    
    allocate(ur_weights(qsplit))
    ur_weights(:)=0.0_r8
    
    if(mod(qsplit,2).NE.0)then
      ur_weights(1)=1.0_r8/qsplit
      do i=3,qsplit,2
        ur_weights(i)=2.0_r8/qsplit
      enddo
    else
      do i=2,qsplit,2
        ur_weights(i)=2.0_r8/qsplit
      enddo
    endif
  end subroutine prim_advance_init
  
  
  subroutine prim_advance_exp(elem, fvm, deriv, hvcoord, hybrid,dt, tl,  nets, nete)   
    use control_mod,       only: prescribed_wind, tstep_type, qsplit
    use derivative_mod,    only: derivative_t
    use dimensions_mod,    only: np, nlev, ntrac
    use element_mod,       only: element_t
    use hybvcoord_mod,     only: hvcoord_t
    use hybrid_mod,        only: hybrid_t
    use time_mod,          only: TimeLevel_t,  timelevel_qdp, tevolve
    use dimensions_mod,    only: qsize_condensate_loading,qsize_condensate_loading_idx_gll
    use dimensions_mod,    only: qsize_condensate_loading_cp, lcp_moist
    use physconst,         only: cpair
    use fvm_control_volume_mod, only: fvm_struct, n0_fvm
    
    implicit none
    
    type (element_t), intent(inout), target   :: elem(:)
    type(fvm_struct)     , intent(in) :: fvm(:)
    type (derivative_t)  , intent(in) :: deriv
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=r8), intent(in) :: dt
    type (TimeLevel_t)   , intent(in) :: tl
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete
    
    ! Local
    real (kind=r8) :: dt_vis, eta_ave_w
    real (kind=r8) :: dp(np,np)
    integer        :: ie,nm1,n0,np1,k,qn0,i,j,m_cnst, nq
    real (kind=r8) :: inv_cp_full(np,np,nlev,nets:nete), sum_cp, sum_water
    real (kind=r8) :: qwater(np,np,nlev,qsize_condensate_loading,nets:nete)

    call t_startf('prim_advance_exp')
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    
    call TimeLevel_Qdp(tl, qsplit, qn0)  ! compute current Qdp() timelevel
    !
    !   tstep_type=1  RK2-SSP 3 stage (as used by tracers)           CFL=.58
    !                    optimal in terms of SSP CFL, but not        CFLSSP=2
    !                    optimal in terms of CFL
    !                    typically requires qsplit=3
    !                    but if windspeed > 340m/s, could use this
    !                    with qsplit=1
    !   tstep_type=2  classic RK3                                    CFL=1.73 (sqrt(3))
    !
    !   tstep_type=3  Kinnmark&Gray RK4 4 stage                      CFL=sqrt(8)=2.8
    !                 should we replace by standard RK4 (CFL=sqrt(8))?
    !                 (K&G 1st order method has CFL=3)
    !   tstep_type=4  Kinnmark&Gray RK3 5 stage 3rd order            CFL=3.87  (sqrt(15))
    !                 From Paul Ullrich.  3rd order for nonlinear terms also
    !                 K&G method is only 3rd order for linear
    !                 optimal: for windspeeds ~120m/s,gravity: 340m/2
    !                 run with qsplit=1
    !                 (K&G 2nd order method has CFL=4. tiny CFL improvement not worth 2nd order)
    !

#ifdef _OPENMP
    call omp_set_nested(.true.)
#endif

    ! default weights for computing mean dynamics fluxes
    eta_ave_w = 1_r8/qsplit

    if (1==prescribed_wind .and. .not.se_prescribed_wind_2d) then
      do ie=nets,nete
        do k=1,nlev
          elem(ie)%state%dp3d(:,:,k,np1) = elem(ie)%state%dp3d(:,:,k,n0)
        enddo
      end do

      
      
      do ie=nets,nete
        ! subcycling code uses a mean flux to advect tracers
 !$omp parallel do num_threads (vert_num_threads) private(dp)
        do k=1,nlev
          dp(:,:) = elem(ie)%state%dp3d(:,:,k,tl%n0)
          
          elem(ie)%derived%vn0(:,:,1,k)=elem(ie)%derived%vn0(:,:,1,k)+&
               eta_ave_w*elem(ie)%state%v(:,:,1,k,n0)*dp(:,:)
          elem(ie)%derived%vn0(:,:,2,k)=elem(ie)%derived%vn0(:,:,2,k)+&
               eta_ave_w*elem(ie)%state%v(:,:,2,k,n0)*dp(:,:)
        enddo
      end do
      call t_stopf('prim_advance_exp')
      return
    endif    
    
    ! ==================================
    ! Take timestep
    ! ==================================

    do ie=nets,nete
      do nq=1,qsize_condensate_loading
        m_cnst = qsize_condensate_loading_idx_gll(nq)
        !
        ! make sure Q is updated
        !
        qwater(:,:,:,nq,ie)      = elem(ie)%state%Qdp(:,:,:,m_cnst,qn0)/elem(ie)%state%dp3d(:,:,:,n0)
      end do
    end do
    !
    ! compute Cp here and not in RK-stages since Q stays constant Cp also stays constant
    !
    if (lcp_moist) then
      do ie=nets,nete
        do k=1,nlev
          do j=1,np
            do i=1,np
              sum_cp    = cpair
              sum_water = 1.0_r8
              do nq=1,qsize_condensate_loading
                sum_cp    = sum_cp+qsize_condensate_loading_cp(nq)*qwater(i,j,k,nq,ie)
                sum_water = sum_water + qwater(i,j,k,nq,ie)
              end do
              inv_cp_full(i,j,k,ie)   = sum_water/sum_cp
            end do
          end do
        end do
      end do
    else
      do ie=nets,nete
        inv_cp_full(:,:,:,nets:nete) = 1.0_r8/cpair
      end do
    end if
    
    dt_vis = dt
    if (tstep_type==1) then
      ! RK2-SSP 3 stage.  matches tracer scheme. optimal SSP CFL, but
      ! not optimal for regular CFL
      ! u1 = u0 + dt/2 RHS(u0)
      call compute_and_apply_rhs(np1,n0,n0,dt/2,elem,hvcoord,hybrid,&
           deriv,nets,nete,eta_ave_w/3,inv_cp_full,qwater)
      ! u2 = u1 + dt/2 RHS(u1)
      call compute_and_apply_rhs(np1,np1,np1,dt/2,elem,hvcoord,hybrid,&
           deriv,nets,nete,eta_ave_w/3,inv_cp_full,qwater)
      ! u3 = u2 + dt/2 RHS(u2)
      call compute_and_apply_rhs(np1,np1,np1,dt/2,elem,hvcoord,hybrid,&
           deriv,nets,nete,eta_ave_w/3,inv_cp_full,qwater)
      
      ! unew = u/3 +2*u3/3  = u + 1/3 (RHS(u) + RHS(u1) + RHS(u2))
      do ie=nets,nete
        elem(ie)%state%v(:,:,:,:,np1)= elem(ie)%state%v(:,:,:,:,n0)/3 &
             + 2*elem(ie)%state%v(:,:,:,:,np1)/3
        elem(ie)%state%T(:,:,:,np1)= elem(ie)%state%T(:,:,:,n0)/3 &
             + 2*elem(ie)%state%T(:,:,:,np1)/3
        elem(ie)%state%dp3d(:,:,:,np1)= elem(ie)%state%dp3d(:,:,:,n0)/3 &
             + 2*elem(ie)%state%dp3d(:,:,:,np1)/3
      enddo
    else if (tstep_type==2) then
      ! classic RK3  CFL=sqrt(3)
      ! u1 = u0 + dt/3 RHS(u0)
      call compute_and_apply_rhs(np1,n0,n0,dt/3,elem,hvcoord,hybrid,&
           deriv,nets,nete,0.0_r8,inv_cp_full,qwater)
      ! u2 = u0 + dt/2 RHS(u1)
      call compute_and_apply_rhs(np1,n0,np1,dt/2,elem,hvcoord,hybrid,&
           deriv,nets,nete,0.0_r8,inv_cp_full,qwater)
      ! u3 = u0 + dt RHS(u2)
      call compute_and_apply_rhs(np1,n0,np1,dt,elem,hvcoord,hybrid,&
           deriv,nets,nete,eta_ave_w,inv_cp_full,qwater)
    else if (tstep_type==3) then
      ! KG 4th order 4 stage:   CFL=sqrt(8)
      ! low storage version of classic RK4
      ! u1 = u0 + dt/4 RHS(u0)
      call compute_and_apply_rhs(np1,n0,n0,dt/4,elem,hvcoord,hybrid,&
           deriv,nets,nete,0.0_r8,inv_cp_full,qwater)
      ! u2 = u0 + dt/3 RHS(u1)
      call compute_and_apply_rhs(np1,n0,np1,dt/3,elem,hvcoord,hybrid,&
           deriv,nets,nete,0.0_r8,inv_cp_full,qwater)
      ! u3 = u0 + dt/2 RHS(u2)
      call compute_and_apply_rhs(np1,n0,np1,dt/2,elem,hvcoord,hybrid,&
           deriv,nets,nete,0.0_r8,inv_cp_full,qwater)
      ! u4 = u0 + dt RHS(u3)
      call compute_and_apply_rhs(np1,n0,np1,dt,elem,hvcoord,hybrid,&
           deriv,nets,nete,eta_ave_w,inv_cp_full,qwater)
    else if (tstep_type==4) then
      !
      ! Ullrich 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
      ! u1 = u0 + dt/5 RHS(u0)  (save u1 in timelevel nm1)
      ! rhs: t=t
      call compute_and_apply_rhs(nm1,n0,n0,dt/5,elem,hvcoord,hybrid,&
           deriv,nets,nete,eta_ave_w/4,inv_cp_full,qwater)
      !
      ! u2 = u0 + dt/5 RHS(u1); rhs: t=t+dt/5
      !
      call compute_and_apply_rhs(np1,n0,nm1,dt/5,elem,hvcoord,hybrid,&
           deriv,nets,nete,0.0_r8,inv_cp_full,qwater)
      !
      ! u3 = u0 + dt/3 RHS(u2); rhs: t=t+2*dt/5
      !
      call compute_and_apply_rhs(np1,n0,np1,dt/3,elem,hvcoord,hybrid,&
           deriv,nets,nete,0.0_r8,inv_cp_full,qwater)
      !
      ! u4 = u0 + 2dt/3 RHS(u3); rhs: t=t+2*dt/5+dt/3
      !
      call compute_and_apply_rhs(np1,n0,np1,2*dt/3,elem,hvcoord,hybrid,&
           deriv,nets,nete,0.0_r8,inv_cp_full,qwater)
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
      !
      ! phl: rhs: t=t+2*dt/5+dt/3+3*dt/4         -wrong RK times ...
      !
      call compute_and_apply_rhs(np1,nm1,np1,3*dt/4,elem,hvcoord,hybrid,&
           deriv,nets,nete,3*eta_ave_w/4,inv_cp_full,qwater)
      ! final method is the same as:
      ! u5 = u0 +  dt/4 RHS(u0)) + 3dt/4 RHS(u4)
    else
      call endrun('ERROR: bad choice of tstep_type')
    endif

    ! ==============================================
    ! Time-split Horizontal diffusion: nu.del^2 or nu.del^4
    ! U(*) = U(t+1)  + dt2 * HYPER_DIFF_TERM(t+1)
    ! ==============================================

    call t_startf('advance_hypervis')

    ! note:time step computes u(t+1)= u(t*) + RHS.
    ! for consistency, dt_vis = t-1 - t*, so this is timestep method dependent

    ! forward-in-time, hypervis applied to dp3d
    call advance_hypervis_dp(edge3p1,elem,fvm,hybrid,deriv,np1,qn0,nets,nete,dt_vis,eta_ave_w,&
         inv_cp_full,hvcoord)

    call t_stopf('advance_hypervis')
    !
    ! update psdry
    !
    do ie=nets,nete
      elem(ie)%state%psdry(:,:) = hvcoord%hyai(1)*hvcoord%ps0
      do k=1,nlev
        elem(ie)%state%psdry(:,:) = elem(ie)%state%psdry(:,:)+elem(ie)%state%dp3d(:,:,k,np1)
      end do
    end do
    tevolve=tevolve+dt

#ifdef _OPENMP
    call omp_set_nested(.false.)
#endif

    call t_stopf('prim_advance_exp')
  end subroutine prim_advance_exp


  subroutine applyCAMforcing(elem,fvm,np1,np1_qdp,dt_q,nets,nete,nsubstep)
    use dimensions_mod,         only: np, nc, nlev, qsize, ntrac, nelemd
    use element_mod,            only: element_t
    use time_mod,               only: nsplit
    use control_mod,            only: ftype
    use fvm_control_volume_mod, only: fvm_struct, n0_fvm
    
    type (element_t)     , intent(inout) :: elem(:)
    type(fvm_struct)     , intent(inout) :: fvm(:)
    real (kind=r8), intent(in) :: dt_q
    integer,  intent(in) :: np1,nets,nete,np1_qdp,nsubstep
    
    ! local
    integer :: i,j,k,ie,q
    real (kind=r8) :: v1,dt_local, dt_local_tracer,tmp
    real (kind=r8) :: dt_local_tracer_fvm
    real (kind=r8) :: ftmp(np,np,nlev,qsize,nets:nete) !diagnostics
    real (kind=r8), allocatable :: ftmp_fvm(:,:,:,:,:) !diagnostics
    
    if (ntrac>0) allocate(ftmp_fvm(nc,nc,nlev,ntrac,nets:nete))

    if (ftype==0) then
      !
      ! "Dribble" tendencies: divide total adjustment with nsplit and
      !                       add adjustments to state after each
      !                       vertical remap
      !
      dt_local            = dt_q
      dt_local_tracer     = dt_q
      dt_local_tracer_fvm = dt_q
    else if (ftype==1) then
      !
      ! CAM-FV-stype forcing, i.e. equivalent to updating state once in the
      ! beginning of dynamics
      !
      dt_local            = nsplit*dt_q
      dt_local_tracer     = nsplit*dt_q
      dt_local_tracer_fvm = nsplit*dt_q
      if (nsubstep.ne.1) then
        !
        ! do nothing
        !
        dt_local            = 0.0_r8
        dt_local_tracer     = 0.0_r8
        dt_local_tracer_fvm = 0.0_r8
      end if
    else if (ftype==2) then
      !
      ! do state-update for tracers and "dribbling" forcing for u,v,T
      !
      dt_local            = dt_q
      if (ntrac>0) then
        dt_local_tracer     = dt_q
        dt_local_tracer_fvm = nsplit*dt_q
        if (nsubstep.ne.1) then
          dt_local_tracer_fvm = 0.0_r8
        end if
      else
        dt_local_tracer     = nsplit*dt_q
        dt_local_tracer_fvm = nsplit*dt_q
        if (nsubstep.ne.1) then
          dt_local_tracer     = 0.0_r8
          dt_local_tracer_fvm = 0.0_r8
        end if
      end if
    end if

    do ie=nets,nete
      elem(ie)%state%T(:,:,:,np1) = elem(ie)%state%T(:,:,:,np1) + &
           dt_local*elem(ie)%derived%FT(:,:,:)
      elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + &
           dt_local*elem(ie)%derived%FM(:,:,:,:)
      
#if (defined COLUMN_OPENMP)
    !$omp parallel do private(q,k,i,j,v1)
#endif
      !
      ! tracers
      !
      if (qsize>0.and.dt_local_tracer>0) then
        do q=1,qsize
          do k=1,nlev
            do j=1,np
              do i=1,np
                !
                ! FQ holds q-tendency: (qnew-qold)/dt_physics
                !
                v1 = dt_local_tracer*elem(ie)%derived%FQ(i,j,k,q)
                if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) + v1 < 0 .and. v1<0) then
                  if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) < 0 ) then
                    v1=0  ! Q already negative, dont make it more so
                  else
                    v1 = -elem(ie)%state%Qdp(i,j,k,q,np1_qdp)
                  endif
                endif
                elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp)+v1
                ftmp(i,j,k,q,ie) = dt_local_tracer*&
                     elem(ie)%derived%FQ(i,j,k,q)-v1 !Only used for diagnostics!
              enddo
            enddo
          enddo
        enddo
      else
        ftmp(:,:,:,:,ie) = 0.0_r8
      end if
      if (ntrac>0.and.dt_local_tracer_fvm>0) then
        !
        ! Repeat for the fvm tracers: fc holds tendency (fc_new-fc_old)/dt_physics
        !
        do q = 1, ntrac
          do k = 1, nlev
            do j = 1, nc
              do i = 1, nc
                tmp = dt_local_tracer_fvm*fvm(ie)%fc(i,j,k,q)/fvm(ie)%dp_fvm(i,j,k,n0_fvm)
                v1 = tmp
                if (fvm(ie)%c(i,j,k,q,n0_fvm) + v1 < 0 .and. v1<0) then
                  if (fvm(ie)%c(i,j,k,q,n0_fvm) < 0 ) then
                    v1 = 0  ! C already negative, dont make it more so
                  else
                    v1 = -fvm(ie)%c(i,j,k,q,n0_fvm)
                  end if
                end if
                fvm(ie)%c(i,j,k,q,n0_fvm) = fvm(ie)%c(i,j,k,q,n0_fvm)+ v1
                ftmp_fvm(i,j,k,q,ie) = tmp-v1 !Only used for diagnostics!
              end do
            end do
          end do
        end do
      else
        if (ntrac>0) ftmp_fvm(:,:,:,:,ie) = 0.0_r8
      end if
    end do
    if (ntrac>0) then
      call output_qdp_var_dynamics(ftmp_fvm(:,:,:,:,:),nc,ntrac,nets,nete,'PDC')
    else
      call output_qdp_var_dynamics(ftmp(:,:,:,:,:),np,qsize,nets,nete,'PDC')
    end if
    call calc_tot_energy_dynamics(elem,fvm,nets,nete,np1,np1_qdp,n0_fvm,'dBD')
    if (ftype==1.and.nsubstep==1) call calc_tot_energy_dynamics(elem,fvm,nets,nete,np1,np1_qdp,n0_fvm,'p2d')
    if (ntrac>0) deallocate(ftmp_fvm)
  end subroutine applyCAMforcing


  subroutine advance_hypervis_dp(edge3,elem,fvm,hybrid,deriv,nt,qn0,nets,nete,dt2,eta_ave_w,inv_cp_full,hvcoord)
    !
    !  take one timestep of:
    !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
    !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
    !
    !
    !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
    !
    !
    use dimensions_mod, only: np, np, nlev, nc, ntrac
    use dimensions_mod, only: hypervis_on_plevs
    use control_mod,    only: nu, nu_s, hypervis_subcycle, nu_p, nu_top
    use hybrid_mod,     only: hybrid_t!, get_loop_ranges
    use element_mod,    only: element_t
    use derivative_mod, only: derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
    use derivative_mod, only: subcell_Laplace_fluxes, subcell_dss_fluxes
    use edge_mod,       only: edgevpack, edgevunpack, edgeDGVunpack
    use edgetype_mod,   only: EdgeBuffer_t, EdgeDescriptor_t
    use bndry_mod,      only: bndry_exchange
    use viscosity_mod,  only: biharmonic_wk_dp3d
    use hybvcoord_mod,  only: hvcoord_t
    use fvm_control_volume_mod, only: fvm_struct, n0_fvm
    
    type (hybrid_t)    , intent(in)   :: hybrid
    type (element_t)   , intent(inout), target :: elem(:)
    type(fvm_struct)   , intent(in)   :: fvm(:)
    type (EdgeBuffer_t), intent(inout):: edge3
    type (derivative_t), intent(in  ) :: deriv
    integer            , intent(in)   :: nets,nete, nt, qn0
    real (kind=r8)     , intent(in)   :: inv_cp_full(np,np,nlev,nets:nete)
    type (hvcoord_t)   , intent(in)   :: hvcoord
    real (kind=r8) :: eta_ave_w  ! weighting for mean flux terms
    real (kind=r8) :: dt2

    ! local
    real (kind=r8) :: nu_scale_top, ptop, press
    integer :: k,kptr,i,j,ie,ic
    integer :: kbeg, kend, kblk
    real (kind=r8), dimension(np,np,2,nlev,nets:nete)      :: vtens
    real (kind=r8), dimension(np,np,nlev,nets:nete)        :: ttens
    real (kind=r8), dimension(np,np,nlev,nets:nete)        :: dptens
    real (kind=r8), dimension(np,np,nlev,nets:nete)        :: dptens_ref,dp3d_ref
    real (kind=r8), dimension(0:np+1,0:np+1,nlev)          :: corners
    real (kind=r8), dimension(2,2,2)                       :: cflux
    real (kind=r8), dimension(nc,nc,4,nlev,nets:nete)      :: dpflux
    real (kind=r8), dimension(np,np,nlev)                  :: nabla4_pk,pk,inv_dpk
    type (EdgeDescriptor_t)                                :: desc    
    
    ! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
    !       data is incorrect (offset by a few numbers actually)
    !       removed for now.
    !       real (kind=r8), dimension(:,:), pointer :: spheremp,rspheremp
    !       real (kind=r8), dimension(:,:,:), pointer   :: ps
    
    real (kind=r8), dimension(np,np)   :: lap_t,lap_dp
    real (kind=r8), dimension(np,np,2) :: lap_v
    real (kind=r8)                     :: v1,v2,dt,heating
    real (kind=r8)                     :: temp      (np,np,nlev)
    real (kind=r8)                     :: laplace_fluxes(nc,nc,4)
    real (kind=r8)                     :: tempflux(nc,nc,4)
    real (kind=r8)                     :: rhypervis_subcycle
    
    if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
    call t_startf('advance_hypervis_dp')

    if (hypervis_on_plevs.and.nu_p>0)&
         call calc_dp3d_reference(elem,edge3p1,hybrid,nets,nete,nt,hvcoord,dp3d_ref)

    ! call get_loop_ranges(hybrid,kbeg=kbeg,kend=kend)
    kbeg=1; kend=nlev
    
    kblk = kend - kbeg + 1
    
    dt=dt2/hypervis_subcycle
    
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !  hyper viscosity
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do ic=1,hypervis_subcycle
      call calc_tot_energy_dynamics(elem,fvm,nets,nete,nt,qn0,n0_fvm,'dBH')
      
      rhypervis_subcycle=1.0_r8/real(hypervis_subcycle,kind=r8)
      if (hypervis_on_plevs.and.nu_p>0) then
        call biharmonic_wk_dp3d(elem,dptens_ref,dpflux,ttens,vtens,deriv,edge3,hybrid,nt,nets,nete,kbeg,kend,&
             dptens,dp3d_ref)
      else
        call biharmonic_wk_dp3d(elem,dptens,dpflux,ttens,vtens,deriv,edge3,hybrid,nt,nets,nete,kbeg,kend)
      end if
      do ie=nets,nete
        if (hypervis_on_plevs) then
          !
          ! compute \nabla^4 p_k
          !
          nabla4_pk(:,:,1) = 0.5_r8*dptens_ref(:,:,1,ie)
          do k=2,nlev
            nabla4_pk(:,:,k) = nabla4_pk(:,:,k-1)+0.5_r8*dptens_ref(:,:,k-1,ie)+0.5_r8*dptens_ref(:,:,k,ie)
          end do
          pk(:,:,1) = 0.5_r8*elem(ie)%state%dp3d(:,:,1,nt)
          do k=2,nlev
            pk(:,:,k) = pk(:,:,k-1)+0.5_r8*elem(ie)%state%dp3d(:,:,k-1,nt)+0.5_r8*elem(ie)%state%dp3d(:,:,k,nt)
          end do
          inv_dpk(:,:,1) = 1.0_r8/(pk(:,:,2)-pk(:,:,1))
          do k=2,nlev-1
            inv_dpk(:,:,k) = 1.0_r8/(pk(:,:,k+1)-pk(:,:,k-1))
          end do
          inv_dpk(:,:,nlev) = 1.0_r8/(pk(:,:,nlev)-pk(:,:,nlev-1))
          do k=1,nlev
            !
            ! viscosity on approximate pressure levels (section 3.3.6 in CAM5 scientific documentation; NCAR Tech. Note TN-486)
            !
            ttens(:,:,k,ie)   = ttens(:,:,k,ie)   -nabla4_pk(:,:,k)*(&
                 inv_dpk(:,:,k)*(elem(ie)%state%T(:,:,MIN(k+1,nlev),nt)-elem(ie)%state%T(:,:,MAX(k-1,1),nt)))
!                           vtens(:,:,1,k,ie) = vtens(:,:,1,k,ie) -nabla4_pk(:,:,k)*(&
!                                inv_dpk(:,:,k)*(elem(ie)%state%v(:,:,1,MIN(k+1,nlev),nt)-elem(ie)%state%v(:,:,1,MAX(k-1,1),nt)))
!                           vtens(:,:,2,k,ie) = vtens(:,:,2,k,ie) -nabla4_pk(:,:,k)*(&
!                                inv_dpk(:,:,k)*(elem(ie)%state%v(:,:,2,MIN(k+1,nlev),nt)-elem(ie)%state%v(:,:,2,MAX(k-1,1),nt)))
          end do
          if (nu_p>0) dptens(:,:,:,ie) = dptens_ref(:,:,:,ie) !pressure damping will only be on difference between smoothed dp3d and dp3d
        endif! done correction term
        
        ! compute mean flux
        if (nu_p>0) then
          do k=kbeg,kend
            !OMP_COLLAPSE_SIMD
            !DIR_VECTOR_ALIGNED
            do j=1,np
              do i=1,np
                elem(ie)%derived%dpdiss_ave(i,j,k)=elem(ie)%derived%dpdiss_ave(i,j,k)+&
                     rhypervis_subcycle*eta_ave_w*elem(ie)%state%dp3d(i,j,k,nt)
                elem(ie)%derived%dpdiss_biharmonic(i,j,k)=elem(ie)%derived%dpdiss_biharmonic(i,j,k)+&
                     rhypervis_subcycle*eta_ave_w*dptens(i,j,k,ie)
              enddo
            enddo
          enddo
        endif
        !$omp parallel do num_threads(vert_num_threads) private(lap_t,lap_dp,lap_v,laplace_fluxes,nu_scale_top,ptop,press)
        do k=kbeg,kend
          ! advace in time.
          ! note: DSS commutes with time stepping, so we can time advance and then DSS.
          ! note: weak operators alreayd have mass matrix "included"
          
          ! add regular diffusion in top 3 layers:
          if (nu_top>0 .and. k<=3) then
            call laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),lap_t,var_coef=.false.)
            call laplace_sphere_wk(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),lap_dp,var_coef=.false.)
            call vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),lap_v, var_coef=.false.)
          endif
          !
          ! compute scaling of sponge layer damping (following cd_core.F90 in CAM-FV)
          !
          press = (hvcoord%hyam(k)+hvcoord%hybm(k))*hvcoord%ps0
          ptop  = hvcoord%hyai(1)*hvcoord%ps0
          nu_scale_top = 8.0_r8*(1.0_r8+ tanh(1.0_r8*log(ptop/press))) ! tau will be maximum 8 at model top
          if (nu_scale_top < 1.0_r8) nu_scale_top = 0.0_r8
          
          ! biharmonic terms need a negative sign:
          if (nu_top>0 .and. nu_scale_top>0.0_r8) then
            !OMP_COLLAPSE_SIMD
            !DIR_VECTOR_ALIGNED
            do j=1,np
              do i=1,np
                ttens(i,j,k,ie)   = (-nu_s*ttens(i,j,k,ie) + nu_scale_top*nu_top*lap_t(i,j) )
                dptens(i,j,k,ie)  = (-nu_p*dptens(i,j,k,ie) + nu_scale_top*nu_top*lap_dp(i,j) )
                vtens(i,j,1,k,ie) = (-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                vtens(i,j,2,k,ie) = (-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
              enddo
            enddo
          else
            !OMP_COLLAPSE_SIMD
            !DIR_VECTOR_ALIGNED
            do j=1,np
              do i=1,np
                ttens(i,j,k,ie)   = -nu_s*ttens(i,j,k,ie)
                dptens(i,j,k,ie)  = -nu_p*dptens(i,j,k,ie)
                vtens(i,j,1,k,ie) = -nu*vtens(i,j,1,k,ie)
                vtens(i,j,2,k,ie) = -nu*vtens(i,j,2,k,ie)
              enddo
            enddo
          endif
          
          if (ntrac>0) then
            !OMP_COLLAPSE_SIMD
            !DIR_VECTOR_ALIGNED
            do j=1,nc
              do i=1,nc
                elem(ie)%sub_elem_mass_flux(i,j,:,k) = elem(ie)%sub_elem_mass_flux(i,j,:,k) - &
                     rhypervis_subcycle*eta_ave_w*nu_p*dpflux(i,j,:,k,ie)
              enddo
            enddo
          if (nu_top>0 .and. nu_scale_top>0.0_r8) then
              call subcell_Laplace_fluxes(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),np,nc,laplace_fluxes)
              elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + &
                   rhypervis_subcycle*eta_ave_w*nu_scale_top*nu_top*laplace_fluxes
            endif
          endif
          
          ! NOTE: we will DSS all tendicies, EXCEPT for dp3d, where we DSS the new state
          !OMP_COLLAPSE_SIMD
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
              elem(ie)%state%dp3d(i,j,k,nt) = elem(ie)%state%dp3d(i,j,k,nt)*elem(ie)%spheremp(i,j)&
                   + dt*dptens(i,j,k,ie)
            enddo
          enddo
          
        enddo
        
        kptr = kbeg - 1
        call edgeVpack(edge3,ttens(:,:,kbeg:kend,ie),kblk,kptr,ie)
        
        kptr = kbeg - 1 + nlev
        call edgeVpack(edge3,vtens(:,:,1,kbeg:kend,ie),kblk,kptr,ie)
        
        kptr = kbeg - 1 + 2*nlev
        call edgeVpack(edge3,vtens(:,:,2,kbeg:kend,ie),kblk,kptr,ie)
        
        kptr = kbeg - 1 + 3*nlev
        call edgeVpack(edge3,elem(ie)%state%dp3d(:,:,kbeg:kend,nt),kblk,kptr,ie)
      enddo
      
      call bndry_exchange(hybrid,edge3,location='advance_hypervis_dp2')
      
      do ie=nets,nete
        
        kptr = kbeg - 1
        call edgeVunpack(edge3,ttens(:,:,kbeg:kend,ie),kblk,kptr,ie)
        
        kptr = kbeg - 1 + nlev
        call edgeVunpack(edge3,vtens(:,:,1,kbeg:kend,ie),kblk,kptr,ie)
        
        kptr = kbeg - 1 + 2*nlev
        call edgeVunpack(edge3,vtens(:,:,2,kbeg:kend,ie),kblk,kptr,ie)
        
        if (ntrac>0) then
          do k=kbeg,kend
            temp(:,:,k) = elem(ie)%state%dp3d(:,:,k,nt) / elem(ie)%spheremp  ! STATE before DSS
            corners(0:np+1,0:np+1,k) = 0.0_r8
            corners(1:np  ,1:np  ,k) = elem(ie)%state%dp3d(1:np,1:np,k,nt) ! fill in interior data of STATE*mass
          enddo
        endif
        kptr = kbeg - 1 + 3*nlev
        call edgeVunpack(edge3,elem(ie)%state%dp3d(:,:,kbeg:kend,nt),kblk,kptr,ie)
        
        if (ntrac>0) then
          desc = elem(ie)%desc
          
          kptr = kbeg - 1 + 3*nlev
          call edgeDGVunpack(edge3,corners(:,:,kbeg:kend),kblk,kptr,ie)
          do k=kbeg,kend
            corners(:,:,k) = corners(:,:,k)/dt !note: array size is 0:np+1
            !OMP_COLLAPSE_SIMD
            !DIR_VECTOR_ALIGNED
            do j=1,np
              do i=1,np
                temp(i,j,k) =  elem(ie)%rspheremp(i,j)*elem(ie)%state%dp3d(i,j,k,nt) - temp(i,j,k)
                temp(i,j,k) =  temp(i,j,k)/dt
              enddo
            enddo
            
            call distribute_flux_at_corners(cflux, corners(:,:,k), desc%getmapP)
            
            cflux(1,1,:)   = elem(ie)%rspheremp(1,  1) * cflux(1,1,:)
            cflux(2,1,:)   = elem(ie)%rspheremp(np, 1) * cflux(2,1,:)
            cflux(1,2,:)   = elem(ie)%rspheremp(1, np) * cflux(1,2,:)
            cflux(2,2,:)   = elem(ie)%rspheremp(np,np) * cflux(2,2,:)
            
            call subcell_dss_fluxes(temp(:,:,k), np, nc, elem(ie)%metdet,cflux,tempflux)
            elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + &
                 rhypervis_subcycle*eta_ave_w*tempflux
          end do
        endif
        
        ! apply inverse mass matrix, accumulate tendencies
        !$omp parallel do num_threads(vert_num_threads)
        do k=kbeg,kend
          !OMP_COLLAPSE_SIMD
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
              vtens(i,j,1,k,ie)=dt*vtens(i,j,1,k,ie)*elem(ie)%rspheremp(i,j)
              vtens(i,j,2,k,ie)=dt*vtens(i,j,2,k,ie)*elem(ie)%rspheremp(i,j)
              ttens(i,j,k,ie)=dt*ttens(i,j,k,ie)*elem(ie)%rspheremp(i,j)
              elem(ie)%state%dp3d(i,j,k,nt)=elem(ie)%state%dp3d(i,j,k,nt)*elem(ie)%rspheremp(i,j)
            enddo
          enddo
        enddo
        
        ! apply hypervis to u -> u+utens:
        ! E0 = dpdn * .5*u dot u + dpdn * T  + dpdn*PHIS
        ! E1 = dpdn * .5*(u+utens) dot (u+utens) + dpdn * (T-X) + dpdn*PHIS
        ! E1-E0:   dpdn (u dot utens) + dpdn .5 utens dot utens   - dpdn X
        !      X = (u dot utens) + .5 utens dot utens
        !  alt:  (u+utens) dot utens
        !$omp parallel do num_threads(vert_num_threads) private(k,i,j)
        do k=kbeg,kend
          !OMP_COLLAPSE_SIMD
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
              ! update v first (gives better results than updating v after heating)
              elem(ie)%state%v(i,j,:,k,nt)=elem(ie)%state%v(i,j,:,k,nt) + &
                   vtens(i,j,:,k,ie)
              elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
                   +ttens(i,j,k,ie)
            enddo
          enddo
        enddo
      end do
      
      call calc_tot_energy_dynamics(elem,fvm,nets,nete,nt,qn0,n0_fvm,'dCH')
      do ie=nets,nete
        !$omp parallel do num_threads(vert_num_threads), private(k,i,j,v1,v2,heating)
        do k=kbeg,kend
          !OMP_COLLAPSE_SIMD
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
              v1=elem(ie)%state%v(i,j,1,k,nt)
              v2=elem(ie)%state%v(i,j,2,k,nt)
              heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )
              elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
                   -heating*inv_cp_full(i,j,k,ie)
            enddo
          enddo
        enddo
      enddo
      
      call calc_tot_energy_dynamics(elem,fvm,nets,nete,nt,qn0,n0_fvm,'dAH')
    enddo

    call t_stopf('advance_hypervis_dp')
   end subroutine advance_hypervis_dp
   
   subroutine compute_and_apply_rhs(np1,nm1,n0,dt2,elem,hvcoord,hybrid,&
        deriv,nets,nete,eta_ave_w,inv_cp_full,qwater)
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
     use dimensions_mod,  only: np, nc, nlev
     use dimensions_mod,  only: qsize_condensate_loading, ntrac
     use hybrid_mod,      only: hybrid_t
     use element_mod,     only: element_t
     use derivative_mod,  only: derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
     use derivative_mod,  only: subcell_div_fluxes, subcell_dss_fluxes
     use edge_mod,        only: edgevpack, edgevunpack, edgeDGVunpack
     use edgetype_mod,    only: edgedescriptor_t
     use bndry_mod,       only: bndry_exchange
     use hybvcoord_mod,   only: hvcoord_t
     use physconst,       only: rair, epsilo
     use prim_si_mod,     only: preq_hydrostatic
     use control_mod,     only: se_met_nudge_u, se_met_nudge_p, se_met_nudge_t, se_met_tevolve
     
     use time_mod, only : tevolve
     
     implicit none
     integer,        intent(in) :: np1,nm1,n0,nets,nete
     real (kind=r8), intent(in) :: dt2
     
     type (hvcoord_t)     , intent(in) :: hvcoord
     type (hybrid_t)      , intent(in) :: hybrid
     type (element_t)     , intent(inout), target :: elem(:)
     type (derivative_t)  , intent(in) :: deriv
     real (kind=r8)       , intent(in) :: inv_cp_full(np,np,nlev,nets:nete)
     real (kind=r8)       , intent(in) :: qwater(np,np,nlev,qsize_condensate_loading,nets:nete)

     real (kind=r8) :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux
     
     ! local
     real (kind=r8), pointer, dimension(:,:,:)                     :: phi
     real (kind=r8), dimension(np,np,nlev)                         :: omega_full
     real (kind=r8), dimension(np,np,nlev)                         :: divdp_dry
     real (kind=r8), dimension(np,np,nlev)                         :: divdp_full
     real (kind=r8), dimension(np,np,2)                            :: vtemp    
     real (kind=r8), dimension(np,np,2,nlev)                       :: vdp_dry  
     real (kind=r8), dimension(np,np,2,nlev)                       :: vdp_full 
     real (kind=r8), dimension(np,np,nlev)                         :: vgrad_p_full
     real (kind=r8), dimension(np,np,2     )                       :: v            !
     real (kind=r8), dimension(np,np)                              :: vgrad_T      ! v.grad(T)
     real (kind=r8), dimension(np,np)                              :: Ephi         ! kinetic energy + PHI term
     real (kind=r8), dimension(np,np,2,nlev)                       :: grad_p_full
     real (kind=r8), dimension(np,np,2,nlev)                       :: grad_p_m_pmet! gradient(p - p_met)
     real (kind=r8), dimension(np,np,nlev)                         :: vort         ! vorticity
     real (kind=r8), dimension(np,np,nlev)                         :: p_dry        ! pressure dry
     real (kind=r8), dimension(np,np,nlev)                         :: dp_dry       ! delta pressure dry
     real (kind=r8), dimension(np,np,nlev)                         :: p_full       ! pressure
     real (kind=r8), dimension(np,np,nlev)                         :: dp_full      ! delta pressure
     real (kind=r8), dimension(0:np+1,0:np+1,nlev)                 :: corners
     real (kind=r8), dimension(2,2,2)                              :: cflux
     real (kind=r8), dimension(np,np)                              :: suml
     real (kind=r8) :: vtens1(np,np,nlev),vtens2(np,np,nlev),ttens(np,np,nlev)
     real (kind=r8) :: stashdp3d (np,np,nlev),tempdp3d(np,np), tempflux(nc,nc,4)
     real (kind=r8) :: inv_epsilon, ckk, term, T_v(np,np,nlev)
     
     type (EdgeDescriptor_t):: desc
     
     real (kind=r8) :: sum_water(np,np,nlev), density_inv
     real (kind=r8) :: E,v1,v2,glnps1,glnps2
     integer        :: i,j,k,kptr,ie,nq
     real (kind=r8) :: u_m_umet, v_m_vmet, t_m_tmet

!JMD  call t_barrierf('sync_compute_and_apply_rhs', hybrid%par%comm)
     inv_epsilon = 1/Epsilo
     call t_adj_detailf(+1)
     call t_startf('compute_and_apply_rhs')
     do ie=nets,nete
       phi => elem(ie)%derived%phi(:,:,:)
       
       ! ==================================================
       ! compute pressure (p) on half levels from ps
       ! using the hybrid coordinates relationship, i.e.
       ! e.g. equation (3.a.92) of the CCM-2 description,
       ! (NCAR/TN-382+STR), June 1993, p. 24.
       ! ==================================================
       
       ! ============================
       ! compute p and delta p
       ! ============================
       
       do k=1,nlev
         ! vertically lagrangian code: we advect dp3d instead of ps
         ! we also need grad(p) at all levels (not just grad(ps))
         !p(k)= hyam(k)*ps0 + hybm(k)*ps
         !    = .5_r8*(hyai(k+1)+hyai(k))*ps0 + .5_r8*(hybi(k+1)+hybi(k))*ps
         !    = .5_r8*(ph(k+1) + ph(k) )  = ph(k) + dp(k)/2
         !
         ! p(k+1)-p(k) = ph(k+1)-ph(k) + (dp(k+1)-dp(k))/2
         !             = dp(k) + (dp(k+1)-dp(k))/2 = (dp(k+1)+dp(k))/2
         dp_dry(:,:,k) = elem(ie)%state%dp3d(:,:,k,n0)
         if (k==1) then
           p_dry(:,:,k)=hvcoord%hyai(k)*hvcoord%ps0 + dp_dry(:,:,k)/2
         else
           p_dry(:,:,k)=p_dry(:,:,k-1) + dp_dry(:,:,k-1)/2 + dp_dry(:,:,k)/2
         endif
         !
         ! compute virtual temperature
         !
         sum_water(:,:,k) = 1.0_r8
         do nq=1,qsize_condensate_loading
           sum_water(:,:,k) = sum_water(:,:,k) + qwater(:,:,k,nq,ie)
         end do
         do j=1,np
           do i=1,np
             t_v(i,j,k)  = elem(ie)%state%T(i,j,k,n0)*(1+inv_epsilon*qwater(i,j,k,1,ie))/sum_water(i,j,k)
           end do
         end do
         !
         ! convert to gas pressure (dry + water vapor pressure)
         ! (assumes T and q are constant in the layer)
         !
         dp_full(:,:,k)=sum_water(:,:,k)*dp_dry(:,:,k)
         if (k==1) then
           p_full(:,:,k) = hvcoord%hyai(k)*hvcoord%ps0 + dp_full(:,:,k)/2
         else
           p_full(:,:,k)=p_full(:,:,k-1) + dp_full(:,:,k-1)/2 + dp_full(:,:,k)/2
         endif
         call gradient_sphere(p_full(:,:,k),deriv,elem(ie)%Dinv,grad_p_full(:,:,:,k))

         
         ! ==============================
         ! compute vgrad_lnps - for omega_full
         ! ==============================
         !OMP_COLLAPSE_SIMD
         !DIR_VECTOR_ALIGNED
         do j=1,np
           do i=1,np
             v1 = elem(ie)%state%v(i,j,1,k,n0)
             v2 = elem(ie)%state%v(i,j,2,k,n0)
             vgrad_p_full(i,j,k) = (v1*grad_p_full(i,j,1,k) + v2*grad_p_full(i,j,2,k))
             vdp_dry(i,j,1,k) = v1*dp_dry(i,j,k)
             vdp_dry(i,j,2,k) = v2*dp_dry(i,j,k)
             vdp_full(i,j,1,k) = v1*dp_full(i,j,k)
             vdp_full(i,j,2,k) = v2*dp_full(i,j,k)
           end do
         end do
         ! ============================
         ! compute grad(P-P_met)
         ! ============================
         if (se_met_nudge_p.gt.0._r8) then
           suml(:,:) = elem(ie)%derived%ps_met(:,:)+tevolve*elem(ie)%derived%dpsdt_met(:,:)
           call gradient_sphere(suml,deriv,elem(ie)%Dinv,vtemp)
           grad_p_m_pmet(:,:,:,k) = grad_p_full(:,:,:,k) - hvcoord%hybm(k)*vtemp
         endif
         ! ================================
         ! Accumulate mean Vel_rho flux in vn0
         ! ================================
         !OMP_COLLAPSE_SIMD 
         !DIR_VECTOR_ALIGNED
         do j=1,np
           do i=1,np
             elem(ie)%derived%vn0(i,j,1,k)=elem(ie)%derived%vn0(i,j,1,k)+eta_ave_w*vdp_dry(i,j,1,k)
             elem(ie)%derived%vn0(i,j,2,k)=elem(ie)%derived%vn0(i,j,2,k)+eta_ave_w*vdp_dry(i,j,2,k)
           enddo
         enddo
         !divdp_dry(:,:,k)
         ! =========================================
         !
         ! Compute relative vorticity and divergence
         !
         ! =========================================
         call divergence_sphere(vdp_dry(:,:,:,k),deriv,elem(ie),divdp_dry(:,:,k))
         call divergence_sphere(vdp_full(:,:,:,k),deriv,elem(ie),divdp_full(:,:,k))
         call vorticity_sphere(elem(ie)%state%v(:,:,:,k,n0),deriv,elem(ie),vort(:,:,k))
       enddo
       
       ! compute T_v for timelevel n0
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j)
#endif
       !
       ! ====================================================
       ! Compute Hydrostatic equation, modelled after CCM-3
       ! ====================================================
       call preq_hydrostatic(phi,elem(ie)%state%phis,T_v,p_full,dp_full)
       ! ====================================================
       ! Compute omega_full
       ! ====================================================
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,j,i,ckk,term)
#endif
       ckk         = 0.5_r8
       suml(:,:  ) = 0
       do k=1,nlev
        !OMP_COLLAPSE_SIMD 
        !DIR_VECTOR_ALIGNED
         do j=1,np   !   Loop inversion (AAM)
           do i=1,np
             term         = -divdp_full(i,j,k)
             
             v1 = elem(ie)%state%v(i,j,1,k,n0)
             v2 = elem(ie)%state%v(i,j,2,k,n0)
             
             omega_full(i,j,k) = suml(i,j) + ckk*term+vgrad_p_full(i,j,k)
             suml(i,j)    = suml(i,j)   + term
           end do
         end do
       end do       
#if (defined COLUMN_OPENMP)
     !$omp parallel do private(k)
#endif
       do k=1,nlev  !  Loop index added (AAM)
         elem(ie)%derived%omega(:,:,k) = &
              elem(ie)%derived%omega(:,:,k) + eta_ave_w*omega_full(:,:,k)
       enddo
       ! ==============================================
       ! Compute phi + kinetic energy term: 10*nv*nv Flops
       ! ==============================================
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,v1,v2,E,Ephi,vtemp,vgrad_T,gpterm,glnps1,glnps2)
#endif
       vertloop: do k=1,nlev
        !OMP_COLLAPSE_SIMD 
        !DIR_VECTOR_ALIGNED
         do j=1,np
           do i=1,np
             v1     = elem(ie)%state%v(i,j,1,k,n0)
             v2     = elem(ie)%state%v(i,j,2,k,n0)
             E = 0.5_r8*( v1*v1 + v2*v2 )
             Ephi(i,j)=E+phi(i,j,k)
           end do
         end do
         ! ================================================
         ! compute gradp term (ps/p)*(dp/dps)*T
         ! ================================================
         call gradient_sphere(elem(ie)%state%T(:,:,k,n0),deriv,elem(ie)%Dinv,vtemp)
        !OMP_COLLAPSE_SIMD 
        !DIR_VECTOR_ALIGNED
         do j=1,np
           do i=1,np
             v1     = elem(ie)%state%v(i,j,1,k,n0)
             v2     = elem(ie)%state%v(i,j,2,k,n0)
             vgrad_T(i,j) =  v1*vtemp(i,j,1) + v2*vtemp(i,j,2)
           end do
         end do
         
         
         ! vtemp = grad ( E + PHI )
         ! vtemp = gradient_sphere(Ephi(:,:),deriv,elem(ie)%Dinv)
         call gradient_sphere(Ephi(:,:),deriv,elem(ie)%Dinv,vtemp)
         
         do j=1,np
           do i=1,np
             density_inv = Rair*T_v(i,j,k)/p_full(i,j,k)
             
             glnps1  = density_inv*grad_p_full(i,j,1,k)
             glnps2  = density_inv*grad_p_full(i,j,2,k)
             
             v1     = elem(ie)%state%v(i,j,1,k,n0)
             v2     = elem(ie)%state%v(i,j,2,k,n0)
             
             vtens1(i,j,k) =   &
                  + v2*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                  - vtemp(i,j,1) - glnps1
             
             vtens2(i,j,k) =   &
                  - v1*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                  - vtemp(i,j,2) - glnps2
             ttens(i,j,k)  =  - vgrad_T(i,j) + &
                  density_inv*omega_full(i,j,k)*inv_cp_full(i,j,k,ie)
             
             if (se_prescribed_wind_2d) then
               vtens1(i,j,k) = 0._r8
               vtens2(i,j,k) = 0._r8
               ttens(i,j,k) = 0._r8
             else
               if(se_met_nudge_u.gt.0._r8)then
                 u_m_umet = v1 - &
                      elem(ie)%derived%u_met(i,j,k) - &
                      se_met_tevolve*tevolve*elem(ie)%derived%dudt_met(i,j,k)
                 v_m_vmet = v2 - &
                      elem(ie)%derived%v_met(i,j,k) - &
                      se_met_tevolve*tevolve*elem(ie)%derived%dvdt_met(i,j,k)
                 
                 vtens1(i,j,k) =   vtens1(i,j,k) - se_met_nudge_u*u_m_umet * elem(ie)%derived%nudge_factor(i,j,k)
                 
                 elem(ie)%derived%Utnd(i+(j-1)*np,k) = elem(ie)%derived%Utnd(i+(j-1)*np,k) &
                      + se_met_nudge_u*u_m_umet * elem(ie)%derived%nudge_factor(i,j,k)
                 
                 vtens2(i,j,k) =   vtens2(i,j,k) - se_met_nudge_u*v_m_vmet * elem(ie)%derived%nudge_factor(i,j,k)
                 
                 elem(ie)%derived%Vtnd(i+(j-1)*np,k) = elem(ie)%derived%Vtnd(i+(j-1)*np,k) &
                      + se_met_nudge_u*v_m_vmet * elem(ie)%derived%nudge_factor(i,j,k)
                 
               endif
               
               if(se_met_nudge_p.gt.0._r8)then
                 vtens1(i,j,k) =   vtens1(i,j,k) - se_met_nudge_p*grad_p_m_pmet(i,j,1,k)  * elem(ie)%derived%nudge_factor(i,j,k)
                 vtens2(i,j,k) =   vtens2(i,j,k) - se_met_nudge_p*grad_p_m_pmet(i,j,2,k)  * elem(ie)%derived%nudge_factor(i,j,k)
               endif
               
               if(se_met_nudge_t.gt.0._r8)then
                 t_m_tmet = elem(ie)%state%T(i,j,k,n0) - &
                      elem(ie)%derived%T_met(i,j,k) - &
                      se_met_tevolve*tevolve*elem(ie)%derived%dTdt_met(i,j,k)
                 ttens(i,j,k)  = ttens(i,j,k) - se_met_nudge_t*t_m_tmet * elem(ie)%derived%nudge_factor(i,j,k)
                 elem(ie)%derived%Ttnd(i+(j-1)*np,k) = elem(ie)%derived%Ttnd(i+(j-1)*np,k) &
                      + se_met_nudge_t*t_m_tmet * elem(ie)%derived%nudge_factor(i,j,k)
               endif
             endif
             
           end do
         end do
         
       end do vertloop
       
       ! =========================================================
       ! local element timestep, store in np1.
       ! note that we allow np1=n0 or nm1
       ! apply mass matrix
       ! =========================================================
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
       do k=1,nlev
         !OMP_COLLAPSE_SIMD 
         !DIR_VECTOR_ALIGNED
         do j=1,np
           do i=1,np
             elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%spheremp(i,j)*( elem(ie)%state%v(i,j,1,k,nm1) + dt2*vtens1(i,j,k) )
             elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%spheremp(i,j)*( elem(ie)%state%v(i,j,2,k,nm1) + dt2*vtens2(i,j,k) )
             elem(ie)%state%T(i,j,k,np1) = elem(ie)%spheremp(i,j)*(elem(ie)%state%T(i,j,k,nm1) + dt2*ttens(i,j,k))
             elem(ie)%state%dp3d(i,j,k,np1) = &
                  elem(ie)%spheremp(i,j) * (elem(ie)%state%dp3d(i,j,k,nm1) - &
                  dt2 * (divdp_dry(i,j,k)))
           enddo
         enddo
         
         
         if (ntrac>0.and.eta_ave_w.ne.0._r8) then
           !OMP_COLLAPSE_SIMD 
           !DIR_VECTOR_ALIGNED
           do j=1,np
             do i=1,np
               v(i,j,1) =  elem(ie)%Dinv(i,j,1,1)*vdp_dry(i,j,1,k) + elem(ie)%Dinv(i,j,1,2)*vdp_dry(i,j,2,k)
               v(i,j,2) =  elem(ie)%Dinv(i,j,2,1)*vdp_dry(i,j,1,k) + elem(ie)%Dinv(i,j,2,2)*vdp_dry(i,j,2,k)
             enddo
           enddo
           call subcell_div_fluxes(v, np, nc, elem(ie)%metdet,tempflux)
           elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) - eta_ave_w*tempflux
         end if
       enddo
       ! =========================================================
       !
       ! Pack 
       !
       ! =========================================================
       kptr=0
       call edgeVpack(edge3p1, elem(ie)%state%T(:,:,:,np1),nlev,kptr,ie)
       
       kptr=nlev
       call edgeVpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,ie)
       
       kptr=kptr+2*nlev
       call edgeVpack(edge3p1, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr, ie)
     end do
     
     ! =============================================================
     ! Insert communications here: for shared memory, just a single
     ! sync is required
     ! =============================================================
     call bndry_exchange(hybrid,edge3p1)
     do ie=nets,nete
       ! ===========================================================
       ! Unpack the edges for vgrad_T and v tendencies...
       ! ===========================================================
       kptr=0
       call edgeVunpack(edge3p1, elem(ie)%state%T(:,:,:,np1), nlev, kptr, ie)
       
       kptr=nlev
       call edgeVunpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, ie)
       
       if (ntrac>0.and.eta_ave_w.ne.0._r8) then
         do k=1,nlev
           stashdp3d(:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)/elem(ie)%spheremp(:,:)
         end do
       endif
         
       corners = 0.0_r8
       corners(1:np,1:np,:) = elem(ie)%state%dp3d(:,:,:,np1)
       kptr=kptr+2*nlev
       call edgeVunpack(edge3p1, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr,ie)
       
       if  (ntrac>0.and.eta_ave_w.ne.0._r8) then
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
           
           call subcell_dss_fluxes(tempdp3d, np, nc, elem(ie)%metdet, cflux,tempflux)
           elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + eta_ave_w*tempflux
         end do
       end if
       
       ! ====================================================
       ! Scale tendencies by inverse mass matrix
       ! ====================================================
       
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
       do k=1,nlev
         !OMP_COLLAPSE_SIMD 
         !DIR_VECTOR_ALIGNED
         do j=1,np
         do i=1,np
             elem(ie)%state%T(i,j,k,np1)   = elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,np1)
             elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,np1)
             elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,np1)
         enddo
         enddo
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
     call t_adj_detailf(-1)
   end subroutine compute_and_apply_rhs
   

   !
   ! corner fluxes for CSLAM
   !
   subroutine distribute_flux_at_corners(cflux, corners, getmapP)
     use dimensions_mod, only : np, max_corner_elem
     use control_mod,    only : swest
     
     real(r8), intent(out) :: cflux(2,2,2)
     real(r8), intent(in)  :: corners(0:np+1,0:np+1)
     integer,  intent(in)  :: getmapP(:)
     
     cflux = 0.0_r8
     if (getmapP(swest+0*max_corner_elem) /= -1) then
       cflux(1,1,1) =                (corners(0,1) - corners(1,1))
       cflux(1,1,1) = cflux(1,1,1) + (corners(0,0) - corners(1,1)) / 2.0_r8
       cflux(1,1,1) = cflux(1,1,1) + (corners(0,1) - corners(1,0)) / 2.0_r8
       
       cflux(1,1,2) =                (corners(1,0) - corners(1,1))
       cflux(1,1,2) = cflux(1,1,2) + (corners(0,0) - corners(1,1)) / 2.0_r8
       cflux(1,1,2) = cflux(1,1,2) + (corners(1,0) - corners(0,1)) / 2.0_r8
     else
       cflux(1,1,1) =                (corners(0,1) - corners(1,1))
       cflux(1,1,2) =                (corners(1,0) - corners(1,1))
     endif
     
     if (getmapP(swest+1*max_corner_elem) /= -1) then
       cflux(2,1,1) =                (corners(np+1,1) - corners(np,1))
       cflux(2,1,1) = cflux(2,1,1) + (corners(np+1,0) - corners(np,1)) / 2.0_r8
       cflux(2,1,1) = cflux(2,1,1) + (corners(np+1,1) - corners(np,0)) / 2.0_r8
       
       cflux(2,1,2) =                (corners(np  ,0) - corners(np,  1))
       cflux(2,1,2) = cflux(2,1,2) + (corners(np+1,0) - corners(np,  1)) / 2.0_r8
       cflux(2,1,2) = cflux(2,1,2) + (corners(np  ,0) - corners(np+1,1)) / 2.0_r8
     else
       cflux(2,1,1) =                (corners(np+1,1) - corners(np,1))
       cflux(2,1,2) =                (corners(np  ,0) - corners(np,1))
     endif
     
     if (getmapP(swest+2*max_corner_elem) /= -1) then
       cflux(1,2,1) =                (corners(0,np  ) - corners(1,np  ))
       cflux(1,2,1) = cflux(1,2,1) + (corners(0,np+1) - corners(1,np  )) / 2.0_r8
       cflux(1,2,1) = cflux(1,2,1) + (corners(0,np  ) - corners(1,np+1)) / 2.0_r8
       
       cflux(1,2,2) =                (corners(1,np+1) - corners(1,np  ))
       cflux(1,2,2) = cflux(1,2,2) + (corners(0,np+1) - corners(1,np  )) / 2.0_r8
       cflux(1,2,2) = cflux(1,2,2) + (corners(1,np+1) - corners(0,np  )) / 2.0_r8
     else
       cflux(1,2,1) =                (corners(0,np  ) - corners(1,np  ))
       cflux(1,2,2) =                (corners(1,np+1) - corners(1,np  ))
     endif
     
     if (getmapP(swest+3*max_corner_elem) /= -1) then
       cflux(2,2,1) =                (corners(np+1,np  ) - corners(np,np  ))
       cflux(2,2,1) = cflux(2,2,1) + (corners(np+1,np+1) - corners(np,np  )) / 2.0_r8
       cflux(2,2,1) = cflux(2,2,1) + (corners(np+1,np  ) - corners(np,np+1)) / 2.0_r8
       
       cflux(2,2,2) =                (corners(np  ,np+1) - corners(np,np  ))
       cflux(2,2,2) = cflux(2,2,2) + (corners(np+1,np+1) - corners(np,np  )) / 2.0_r8
       cflux(2,2,2) = cflux(2,2,2) + (corners(np  ,np+1) - corners(np+1,np)) / 2.0_r8
     else
       cflux(2,2,1) =                (corners(np+1,np  ) - corners(np,np  ))
       cflux(2,2,2) =                (corners(np  ,np+1) - corners(np,np  ))
     endif
   end subroutine distribute_flux_at_corners


  subroutine calc_tot_energy_dynamics(elem,fvm,nets,nete,tl,tl_qdp,n_fvm,outfld_name_suffix)
    use dimensions_mod,         only: npsq,nlev,np,lcp_moist,nc,ntrac
    use dimensions_mod,         only: qsize_condensate_loading, qsize_condensate_loading_cp
    use dimensions_mod,         only: qsize_condensate_loading_idx_gll
    use physconst,              only: gravit, cpair, rearth,omega
    use element_mod,            only: element_t
    use cam_history,            only: outfld, hist_fld_active
    use constituents,           only: cnst_get_ind
    use hycoef,                 only: hyai, ps0
    use fvm_control_volume_mod, only: fvm_struct, n0_fvm
    !------------------------------Arguments--------------------------------
    
    type (element_t) , intent(in) :: elem(:)
    type(fvm_struct) , intent(in) :: fvm(:)
    integer          , intent(in) :: tl, tl_qdp,nets,nete,n_fvm
    character*(*)    , intent(in) :: outfld_name_suffix ! suffix for "outfld" names
    
    !---------------------------Local storage-------------------------------
    
    real(kind=r8) :: se(npsq)                          ! Dry Static energy (J/m2)
    real(kind=r8) :: ke(npsq)                          ! kinetic energy    (J/m2)

    real(kind=r8) :: cdp_fvm(nc,nc,nlev)
    real(kind=r8) :: se_tmp
    real(kind=r8) :: ke_tmp
    real(kind=r8) :: ps(np,np)
    real(kind=r8) :: pdel
    !
    ! global axial angular momentum (AAM) can be separated into one part (mr) associatedwith the relative motion 
    ! of the atmosphere with respect to the planets surface (also known as wind AAM) and another part (mo) 
    ! associated with the angular velocity OMEGA (2*pi/d, where d is the length of the day) of the planet 
    ! (also known as mass AAM)
    !
    real(kind=r8) :: mr(npsq)  ! wind AAM                        
    real(kind=r8) :: mo(npsq)  ! mass AAM
    real(kind=r8) :: mr_cnst, mo_cnst, cos_lat, mr_tmp, mo_tmp

    integer :: ie,i,j,k,nq,m_cnst
    integer :: ixwv,ixcldice, ixcldliq, ixtt ! CLDICE, CLDLIQ and test tracer indices
    character(len=16) :: name_out1,name_out2,name_out3,name_out4,name_out5,name_out6

    !-----------------------------------------------------------------------

    name_out1 = 'SE_'   //trim(outfld_name_suffix)
    name_out2 = 'KE_'   //trim(outfld_name_suffix)
    name_out3 = 'WV_'   //trim(outfld_name_suffix)
    name_out4 = 'WL_'   //trim(outfld_name_suffix)
    name_out5 = 'WI_'   //trim(outfld_name_suffix)
    name_out6 = 'TT_'   //trim(outfld_name_suffix)

    if ( hist_fld_active(name_out1).or.hist_fld_active(name_out2).or.hist_fld_active(name_out3).or.&
         hist_fld_active(name_out4).or.hist_fld_active(name_out5).or.hist_fld_active(name_out6)) then

      if (ntrac>0) then
        ixwv = 1
        call cnst_get_ind('CLDLIQ' , ixcldliq, abort=.false.)
        call cnst_get_ind('CLDICE' , ixcldice, abort=.false.)
      else
        !
        ! when using CSLAM the condensates on the GLL grid may be located in a different index than in physics
        !
        ixwv = -1!to avoid compiletime error in code below
        ixcldliq = qsize_condensate_loading_idx_gll(2)
        ixcldice = qsize_condensate_loading_idx_gll(3)
      end if
      call cnst_get_ind('TT_LW' , ixtt    , abort=.false.)
      !
      ! Compute frozen static energy in 3 parts:  KE, SE, and energy associated with vapor and liquid
      !
      do ie=nets,nete
        se    = 0.0_r8
        ke    = 0.0_r8

        ps(:,:)    = hyai(1)*ps0
        do k = 1, nlev
          do j=1,np
            do i = 1, np
              pdel = elem(ie)%state%dp3d(i,j,k,tl)
              do nq=1,qsize_condensate_loading
                m_cnst = qsize_condensate_loading_idx_gll(nq)
                pdel = pdel + elem(ie)%state%qdp(i,j,k,m_cnst,tl_qdp)
              end do
              ps(i,j)  = ps(i,j)+pdel
              !
              ! kinetic energy
              !
              ke_tmp   = 0.5_r8*(elem(ie)%state%v(i,j,1,k,tl)**2+ elem(ie)%state%v(i,j,2,k,tl)**2)*pdel/gravit
              if (lcp_moist) then
                !
                ! Internal energy formula including all condensates and corresponding heat capacities
                !
                se_tmp = cpair*elem(ie)%state%dp3d(i,j,k,tl)
                do nq=1,qsize_condensate_loading
                  m_cnst = qsize_condensate_loading_idx_gll(nq)
                  se_tmp = se_tmp+qsize_condensate_loading_cp(nq)*elem(ie)%state%qdp(i,j,k,m_cnst,tl_qdp)
                end do
                se_tmp = se_tmp*elem(ie)%state%T(i,j,k,tl)/gravit
              else
                !
                ! using CAM physics definition of internal energy
                !
                se_tmp   = cpair*elem(ie)%state%T(i,j,k,tl)*pdel/gravit
              end if              
              se   (i+(j-1)*np) = se   (i+(j-1)*np) + se_tmp
              ke   (i+(j-1)*np) = ke   (i+(j-1)*np) + ke_tmp
            end do
          end do
        end do
          
        do j=1,np
          do i = 1, np
            se(i+(j-1)*np) = se(i+(j-1)*np) + elem(ie)%state%phis(i,j)*ps(i,j)/gravit
          end do
        end do
        !
        ! Output energy diagnostics on GLL grid
        !
        call outfld(name_out1  ,se       ,npsq,ie)
        call outfld(name_out2  ,ke       ,npsq,ie)
        !
        ! mass variables are output on CSLAM grid if using CSLAM else GLL grid
        !
        if (ntrac>0) then
          if (ixwv>0) then
            cdp_fvm = fvm(ie)%c(1:nc,1:nc,:,ixwv,n_fvm)*fvm(ie)%dp_fvm(1:nc,1:nc,:,n_fvm)
            call util_function(cdp_fvm,nc,nlev,name_out3,ie)
          end if
          if (ixcldliq>0) then
            cdp_fvm = fvm(ie)%c(1:nc,1:nc,:,ixcldliq,n_fvm)*fvm(ie)%dp_fvm(1:nc,1:nc,:,n_fvm)
            call util_function(cdp_fvm,nc,nlev,name_out4,ie)
          end if
          if (ixcldice>0) then
            cdp_fvm = fvm(ie)%c(1:nc,1:nc,:,ixcldice,n_fvm)*fvm(ie)%dp_fvm(1:nc,1:nc,:,n_fvm) 
            call util_function(cdp_fvm,nc,nlev,name_out5,ie)
          end if
          if (ixtt>0) then
            cdp_fvm = fvm(ie)%c(1:nc,1:nc,:,ixtt,n_fvm)*fvm(ie)%dp_fvm(1:nc,1:nc,:,n_fvm)
            call util_function(cdp_fvm,nc,nlev,name_out6,ie)                        
          end if
        else
          call util_function(elem(ie)%state%qdp(:,:,:,1       ,tl_qdp),np,nlev,name_out3,ie)
          if (ixcldliq>0) call util_function(elem(ie)%state%qdp(:,:,:,ixcldliq,tl_qdp),np,nlev,name_out4,ie)
          if (ixcldice>0) call util_function(elem(ie)%state%qdp(:,:,:,ixcldice,tl_qdp),np,nlev,name_out5,ie)
          if (ixtt>0    ) call util_function(elem(ie)%state%qdp(:,:,:,ixtt    ,tl_qdp),np,nlev,name_out6,ie)
        end if
      end do
    end if
    !
    ! Axial angular momentum diagnostics
    !
    ! Code follows 
    !
    ! Lauritzen et al., (2014): Held-Suarez simulations with the Community Atmosphere Model
    ! Spectral Element (CAM-SE) dynamical core: A global axial angularmomentum analysis using Eulerian
    ! and floating Lagrangian vertical coordinates. J. Adv. Model. Earth Syst. 6,129-140, 
    ! doi:10.1002/2013MS000268
    !
    ! MR is equation (6) without \Delta A and sum over areas (areas are in units of radians**2)
    ! MO is equation (7) without \Delta A and sum over areas (areas are in units of radians**2)
    !
    name_out1 = 'MR_'   //trim(outfld_name_suffix)
    name_out2 = 'MO_'   //trim(outfld_name_suffix)

    if ( hist_fld_active(name_out1).or.hist_fld_active(name_out2)) then
      if (qsize_condensate_loading>1) then
        ixcldliq = qsize_condensate_loading_idx_gll(2)
        ixcldice = qsize_condensate_loading_idx_gll(3)
      else
        ixcldliq = -1
        ixcldice = -1
      end if
      mr_cnst = rearth**3/gravit
      mo_cnst = omega*rearth**4/gravit
      do ie=nets,nete
        mr    = 0.0_r8
        mo    = 0.0_r8
        do k = 1, nlev
          do j=1,np
            do i = 1, np
              pdel = elem(ie)%state%dp3d(i,j,k,tl)
              do nq=1,qsize_condensate_loading
                m_cnst = qsize_condensate_loading_idx_gll(nq)
                pdel = pdel + elem(ie)%state%qdp(i,j,k,m_cnst,tl_qdp)
              end do

              cos_lat = cos(elem(ie)%spherep(i,j)%lat)
              mr_tmp   = mr_cnst*elem(ie)%state%v(i,j,1,k,tl)*pdel*cos_lat
              mo_tmp   = mo_cnst*pdel*cos_lat**2
              
              mr   (i+(j-1)*np) = mr   (i+(j-1)*np) + mr_tmp
              mo   (i+(j-1)*np) = mo   (i+(j-1)*np) + mo_tmp
            end do
          end do
        end do
        call outfld(name_out1  ,mr       ,npsq,ie)
        call outfld(name_out2  ,mo       ,npsq,ie)
      end do
    end if


  end subroutine calc_tot_energy_dynamics


  subroutine output_qdp_var_dynamics(qdp,nx,num_trac,nets,nete,outfld_name)
    use dimensions_mod, only: nlev,ntrac,nelemd
    use physconst     , only: gravit
    use cam_history   , only: outfld, hist_fld_active
    use constituents  , only: cnst_get_ind
    use control_mod,    only: TRACERTRANSPORT_SE_GLL, tracer_transport_type
    use dimensions_mod, only: qsize_condensate_loading,qsize_condensate_loading_idx_gll
    use dimensions_mod, only: qsize_condensate_loading_idx
    !------------------------------Arguments--------------------------------

    integer      ,intent(in) :: nx,num_trac,nets,nete
    real(kind=r8) :: qdp(nx,nx,nlev,num_trac,nets:nete)
    character*(*),intent(in) :: outfld_name

    !---------------------------Local storage-------------------------------

    integer :: ie
    integer :: ixcldice, ixcldliq, ixtt
    character(len=16) :: name_out1,name_out2,name_out3,name_out4
    !-----------------------------------------------------------------------

    name_out1 = 'WV_'   //trim(outfld_name)
    name_out2 = 'WI_'   //trim(outfld_name)
    name_out3 = 'WL_'   //trim(outfld_name)
    name_out4 = 'TT_'   //trim(outfld_name)

    if ( hist_fld_active(name_out1).or.hist_fld_active(name_out2).or.hist_fld_active(name_out3).or.&
         hist_fld_active(name_out4)) then

      if (qsize_condensate_loading>1) then
        if (ntrac>0) then
          ixcldliq = qsize_condensate_loading_idx(2)
          ixcldice = qsize_condensate_loading_idx(3)
        else
          ixcldliq = qsize_condensate_loading_idx_gll(2)
          ixcldice = qsize_condensate_loading_idx_gll(3)
        end if
      else
        ixcldliq = -1
        ixcldice = -1
      end if
      call cnst_get_ind('TT_LW' , ixtt    , abort=.false.)

      do ie=nets,nete
        call util_function(qdp(:,:,:,1,ie),nx,nlev,name_out1,ie)
        if (ixcldice>0) call util_function(qdp(:,:,:,ixcldice,ie),nx,nlev,name_out2,ie)
        if (ixcldliq>0) call util_function(qdp(:,:,:,ixcldliq,ie),nx,nlev,name_out3,ie)
        if (ixtt>0    ) call util_function(qdp(:,:,:,ixtt    ,ie),nx,nlev,name_out4,ie)
      end do
    end if
  end subroutine output_qdp_var_dynamics

  !
  ! column integrate mass-variable and outfld
  !
  subroutine util_function(f_in,nx,nz,name_out,ie)
    use physconst,   only: gravit
    use cam_history, only: outfld, hist_fld_active
    integer,           intent(in) :: nx,nz,ie
    real(kind=r8),     intent(in) :: f_in(nx,nx,nz)
    character(len=16), intent(in) :: name_out
    real(kind=r8)       :: f_out(nx*nx)
    integer             :: i,j,k
    real(kind=r8)       :: inv_g
    if (hist_fld_active(name_out)) then
      f_out = 0.0_r8      
      inv_g = 1.0_r8/gravit
      do k = 1, nz
        do j = 1, nx
          do i = 1, nx
            f_out(i+(j-1)*nx) = f_out(i+(j-1)*nx) + f_in(i,j,k)
          end do
        end do
      end do
      f_out = f_out*inv_g
      call outfld(name_out,f_out,nx*nx,ie)
    end if
  end subroutine util_function


   subroutine compute_omega(hybrid,n0,qn0,elem,deriv,nets,nete,dt,hvcoord)
     use control_mod,    only : nu_p, hypervis_subcycle
     use edge_mod,       only : initEdgeBuffer, FreeEdgeBuffer
     use dimensions_mod, only : np, nlev
     use hybrid_mod,     only : hybrid_t
     use element_mod,    only : element_t
     use derivative_mod, only : divergence_sphere, derivative_t,gradient_sphere
     use hybvcoord_mod,  only : hvcoord_t
     use edge_mod,       only : edgevpack, edgevunpack
     use bndry_mod,      only : bndry_exchange
     use viscosity_mod,  only: biharmonic_wk_omega
     use dimensions_mod, only: qsize_condensate_loading,qsize_condensate_loading_idx_gll
     implicit none
     type (hybrid_t)      , intent(in)            :: hybrid
     type (element_t)     , intent(inout), target :: elem(:)
     type (derivative_t)  , intent(in)            :: deriv
     integer              , intent(in)            :: nets,nete,n0,qn0
     real (kind=r8)       , intent(in)            :: dt
     type (hvcoord_t)     , intent(in)            :: hvcoord

     integer        :: i,j,k,ie,kptr,ic,m_cnst,nq
     real (kind=r8) :: ckk, suml(np,np), v1, v2, term
     real (kind=r8) :: dp_full(np,np,nlev)
     real (kind=r8) :: p_full(np,np,nlev),grad_p_full(np,np,2),vgrad_p_full(np,np,nlev)
     real (kind=r8) :: divdp_full(np,np,nlev),vdp_full(np,np,2)
     real(kind=r8)  :: Otens(np,np  ,nlev,nets:nete), dt_hyper, sum_water(np,np,nlev)

     type (EdgeBuffer_t) :: edgeOmega
     logical, parameter  :: del4omega = .true.

     call initEdgeBuffer(hybrid%par,edgeOmega,elem,nlev,bndry_type=HME_BNDRY_P2P, nthreads=horz_num_threads)
     do ie=nets,nete
        sum_water(:,:,:) = 1.0_r8
        do nq=1,qsize_condensate_loading
          m_cnst = qsize_condensate_loading_idx_gll(nq)
          sum_water(:,:,:) = sum_water(:,:,:) + &
               elem(ie)%state%Qdp(:,:,:,m_cnst,qn0)/elem(ie)%state%dp3d(:,:,:,n0)
        end do
        do k=1,nlev
           dp_full(:,:,k)=sum_water(:,:,k)*elem(ie)%state%dp3d(:,:,k,n0)
           if (k==1) then
              p_full(:,:,k) = hvcoord%hyai(k)*hvcoord%ps0 + dp_full(:,:,k)/2
           else
              p_full(:,:,k)=p_full(:,:,k-1) + dp_full(:,:,k-1)/2 + dp_full(:,:,k)/2
           endif
           call gradient_sphere(p_full(:,:,k),deriv,elem(ie)%Dinv,grad_p_full (:,:,:))
         do j=1,np
           do i=1,np
             v1 = elem(ie)%state%v(i,j,1,k,n0)
             v2 = elem(ie)%state%v(i,j,2,k,n0)
             vdp_full(i,j,1) = dp_full(i,j,k)*v1
             vdp_full(i,j,2) = dp_full(i,j,k)*v2
             vgrad_p_full(i,j,k) = (v1*grad_p_full(i,j,1) + v2*grad_p_full(i,j,2))
           end do
         end do
         call divergence_sphere(vdp_full(:,:,:),deriv,elem(ie),divdp_full(:,:,k))
       end do
       ckk       = 0.5_r8
       suml(:,:  ) = 0
       do k=1,nlev
         do j=1,np   !   Loop inversion (AAM)
           do i=1,np
             term         = -divdp_full(i,j,k)
             
             v1 = elem(ie)%state%v(i,j,1,k,n0)
             v2 = elem(ie)%state%v(i,j,2,k,n0)

             elem(ie)%derived%omega(i,j,k) = suml(i,j) + ckk*term+vgrad_p_full(i,j,k)          
             
             suml(i,j)    = suml(i,j)   + term
           end do
         end do
       end do
     end do
     do ie=nets,nete
       do k=1,nlev
         elem(ie)%derived%omega(:,:,k) = elem(ie)%spheremp(:,:)*elem(ie)%derived%omega(:,:,k)
       end do
       kptr=0
       call edgeVpack(edgeOmega, elem(ie)%derived%omega(:,:,:),nlev,kptr, ie)
     end do
     call bndry_exchange(hybrid,edgeOmega)
     do ie=nets,nete
       kptr=0
       call edgeVunpack(edgeOmega, elem(ie)%derived%omega(:,:,:),nlev,kptr, ie)
       do k=1,nlev
         elem(ie)%derived%omega(:,:,k) = elem(ie)%rspheremp(:,:)*elem(ie)%derived%omega(:,:,k)
       end do
     end do

     if  (del4omega) then
       dt_hyper=dt/hypervis_subcycle
       do ic=1,hypervis_subcycle
         do ie=nets,nete
           Otens(:,:,:,ie) = elem(ie)%derived%omega(:,:,:)
         end do
         call biharmonic_wk_omega(elem,Otens,deriv,edgeOmega,hybrid,nets,nete,1,nlev)
         do ie=nets,nete
           do k=1,nlev
             Otens(:,:,k,ie) = -dt_hyper*nu_p*Otens(:,:,k,ie)
           end do
           kptr=0
           call edgeVpack(edgeOmega,Otens(:,:,:,ie) ,nlev,kptr, ie)
         end do
         call bndry_exchange(hybrid,edgeOmega)
         do ie=nets,nete
           kptr=0
           call edgeVunpack(edgeOmega, Otens(:,:,:,ie),nlev,kptr, ie)
         end do
         do ie=nets,nete
           do k=1,nlev
             elem(ie)%derived%omega(:,:,k) =elem(ie)%derived%omega(:,:,k)+&
                  elem(ie)%rspheremp(:,:)*Otens(:,:,k,ie)
           end do
         end do
       end do
     end if
     call FreeEdgeBuffer(edgeOmega)
   end subroutine compute_omega

  subroutine calc_dp3d_reference(elem,edge3,hybrid,nets,nete,nt,hvcoord,dp3d_ref)
    !
    ! calc_dp3d_reference: When the del^4 horizontal damping is applied to dp3d 
    !                      the values are implicitly affected by natural variations 
    !                      due to surface topography. 
    !
    !                    To account for these physicaly correct variations, use 
    !                    the current state values to compute appropriate 
    !                    reference values for the current (lagrangian) ETA-surfaces.
    !                    Damping should then be applied to values relative to 
    !                    this reference.
    !=======================================================================
    use hybvcoord_mod  ,only: hvcoord_t
    use physconst      ,only: rair,cappa
    use element_mod,    only: element_t
    use dimensions_mod, only: np,nlev
    use hybrid_mod,     only: hybrid_t
    use edge_mod,       only: edgevpack, edgevunpack
    use bndry_mod,      only: bndry_exchange
    !                  
    ! Passed variables
    !-------------------
    type(element_t   ),target,intent(inout):: elem(:)
    type(EdgeBuffer_t)       ,intent(inout):: edge3
    type(hybrid_t    )       ,intent(in   ):: hybrid
    integer                  ,intent(in   ):: nets,nete
    integer                  ,intent(in   ):: nt
    type(hvcoord_t   )       ,intent(in   ):: hvcoord
    real(kind=r8)            ,intent(out  ):: dp3d_ref(np,np,nlev,nets:nete)
    !
    ! Local Values
    !--------------
    real(kind=r8):: Phis_avg(np,np,     nets:nete)
    real(kind=r8):: Phi_avg (np,np,nlev,nets:nete)
    real(kind=r8):: RT_avg  (np,np,nlev,nets:nete)
    real(kind=r8):: P_val   (np,np,nlev)
    real(kind=r8):: Ps_val  (np,np)
    real(kind=r8):: Phi_val (np,np,nlev)
    real(kind=r8):: Phi_ival(np,np)
    real(kind=r8):: I_Phi   (np,np,nlev+1)
    real(kind=r8):: Alpha   (np,np,nlev  )
    real(kind=r8):: I_P     (np,np,nlev+1)
    real(kind=r8):: DP_avg  (np,np,nlev)
    real(kind=r8):: P_avg   (np,np,nlev)
    real(kind=r8):: Ps_avg  (np,np)
    real(kind=r8):: Ps_ref  (np,np)
    real(kind=r8):: RT_lapse(np,np)
    real(kind=r8):: dlt_Ps  (np,np)
    real(kind=r8):: dPhi    (np,np,nlev)
    real(kind=r8):: dPhis   (np,np)
    real(kind=r8):: E_Awgt,E_phis,E_phi(nlev),E_T(nlev),Lapse0,Expon0
    integer      :: ie,ii,jj,kk,kptr

    ! Loop over elements
    !--------------------
    do ie=nets,nete

      ! Calculate Pressure values from dp3dp
      !--------------------------------------
      P_val(:,:,1) = hvcoord%hyai(1)*hvcoord%ps0 + elem(ie)%state%dp3d(:,:,1,nt)*0.5_r8
      do kk=2,nlev
        P_val(:,:,kk) =               P_val(:,:,kk-1)           &
                      + elem(ie)%state%dp3d(:,:,kk-1,nt)*0.5_r8 &
                      + elem(ie)%state%dp3d(:,:,kk  ,nt)*0.5_r8
      end do
      Ps_val(:,:) = P_val(:,:,nlev) + elem(ie)%state%dp3d(:,:,nlev,nt)*0.5_r8

      ! Calculate (dry) geopotential values
      !--------------------------------------
      dPhi    (:,:,:)    = 0.5_r8*(rair*elem(ie)%state%T   (:,:,:,nt) &
                                      *elem(ie)%state%dp3d(:,:,:,nt) &
                                                    /P_val(:,:,:)    )
      Phi_val (:,:,nlev) = elem(ie)%state%phis(:,:) + dPhi(:,:,nlev)
      Phi_ival(:,:)      = elem(ie)%state%phis(:,:) + dPhi(:,:,nlev)*2._r8
      do kk=(nlev-1),1,-1
        Phi_val (:,:,kk) = Phi_ival(:,:)    + dPhi(:,:,kk)
        Phi_ival(:,:)    = Phi_val (:,:,kk) + dPhi(:,:,kk)
      end do

      ! Calculate Element averages
      !----------------------------
      E_Awgt   = 0.0_r8
      E_phis   = 0.0_r8
      E_phi(:) = 0._r8
      E_T  (:) = 0._r8
      do jj=1,np
      do ii=1,np
        E_Awgt    = E_Awgt    + elem(ie)%spheremp(ii,jj)
        E_phis    = E_phis    + elem(ie)%spheremp(ii,jj)*elem(ie)%state%phis(ii,jj)
        E_phi (:) = E_phi (:) + elem(ie)%spheremp(ii,jj)*Phi_val(ii,jj,:)
        E_T   (:) = E_T   (:) + elem(ie)%spheremp(ii,jj)*elem(ie)%state%T(ii,jj,:,nt)
      end do
      end do

      Phis_avg(:,:,ie) = E_phis/E_Awgt
      do kk=1,nlev
        Phi_avg(:,:,kk,ie) = E_phi(kk)     /E_Awgt
        RT_avg (:,:,kk,ie) = E_T  (kk)*rair/E_Awgt
      end do
    end do ! ie=nets,nete

    ! Boundary Exchange of average values
    !-------------------------------------
    do ie=nets,nete
      Phis_avg(:,:,ie) = elem(ie)%spheremp(:,:)*Phis_avg(:,:,ie)
      do kk=1,nlev
        Phi_avg(:,:,kk,ie) = elem(ie)%spheremp(:,:)*Phi_avg(:,:,kk,ie)
        RT_avg (:,:,kk,ie) = elem(ie)%spheremp(:,:)*RT_avg (:,:,kk,ie)
      end do
      kptr = 0
      call edgeVpack(edge3,Phi_avg(:,:,:,ie),nlev,kptr,ie)
      kptr = nlev
      call edgeVpack(edge3,RT_avg (:,:,:,ie),nlev,kptr,ie)
      kptr = 2*nlev
      call edgeVpack(edge3,Phis_avg (:,:,ie),1   ,kptr,ie)
    end do ! ie=nets,nete

    call bndry_exchange(hybrid,edge3)

    do ie=nets,nete
      kptr = 0
      call edgeVunpack(edge3,Phi_avg(:,:,:,ie),nlev,kptr,ie)
      kptr = nlev
      call edgeVunpack(edge3,RT_avg (:,:,:,ie),nlev,kptr,ie)
      kptr = 2*nlev
      call edgeVunpack(edge3,Phis_avg (:,:,ie),1   ,kptr,ie)
      Phis_avg(:,:,ie) = elem(ie)%rspheremp(:,:)*Phis_avg(:,:,ie)
      do kk=1,nlev
        Phi_avg(:,:,kk,ie) = elem(ie)%rspheremp(:,:)*Phi_avg(:,:,kk,ie)
        RT_avg (:,:,kk,ie) = elem(ie)%rspheremp(:,:)*RT_avg (:,:,kk,ie)
      end do
    end do ! ie=nets,nete

    ! Loop over elements
    !--------------------
    do ie=nets,nete

      ! Fill elements with uniformly varying average values
      !-----------------------------------------------------
      call fill_element(Phis_avg(1,1,ie))
      do kk=1,nlev
        call fill_element(Phi_avg(1,1,kk,ie))
        call fill_element(RT_avg (1,1,kk,ie))
      end do

      ! Integrate upward to compute Alpha == (dp3d/P)
      !----------------------------------------------
      I_Phi(:,:,nlev+1) = Phis_avg(:,:,ie)
      do kk=nlev,1,-1
        I_Phi(:,:,kk) = 2._r8* Phi_avg(:,:,kk,ie) - I_Phi(:,:,kk+1)
        Alpha(:,:,kk) = 2._r8*(Phi_avg(:,:,kk,ie) - I_Phi(:,:,kk+1))/RT_avg(:,:,kk,ie)
      end do

      ! Integrate downward to compute corresponding average pressure values 
      !---------------------------------------------------------------------
      I_P(:,:,1) = hvcoord%hyai(1)*hvcoord%ps0
      do kk=1,nlev
        DP_avg(:,:,kk  ) = I_P(:,:,kk)*(2._r8 * Alpha(:,:,kk))/(2._r8 - Alpha(:,:,kk))
        P_avg (:,:,kk  ) = I_P(:,:,kk)*(2._r8                )/(2._r8 - Alpha(:,:,kk))
        I_P   (:,:,kk+1) = I_P(:,:,kk)*(2._r8 + Alpha(:,:,kk))/(2._r8 - Alpha(:,:,kk))
      end do
      Ps_avg(:,:) = I_P(:,:,nlev+1)

      ! Determine an appropriate d<T>/d<PHI> lapse rate near the surface
      ! OPTIONALLY: Use dry adiabatic lapse rate or environmental lapse rate.
      !-----------------------------------------------------------------------
      if(.FALSE.) then
        ! DRY ADIABATIC laspe rate
        !------------------------------
        RT_lapse(:,:) = -cappa
      else
        ! ENVIRONMENTAL (empirical) laspe rate
        !--------------------------------------
        RT_lapse(:,:) =  (RT_avg (:,:,nlev-1,ie)-RT_avg (:,:,nlev,ie)) &
                        /(Phi_avg(:,:,nlev-1,ie)-Phi_avg(:,:,nlev,ie))
      endif

      ! Calcualte reference surface pressure
      !--------------------------------------
      dPhis(:,:) = elem(ie)%state%phis(:,:)-Phis_avg(:,:,ie)
      do jj=1,np
      do ii=1,np
        if (abs(RT_lapse(ii,jj)) .gt. 1.e-3_r8) then
          Lapse0 = RT_lapse(ii,jj)/RT_avg(ii,jj,nlev,ie)
          Expon0 = (-1._r8/RT_lapse(ii,jj))
          Ps_ref(ii,jj) = Ps_avg(ii,jj)*((1._r8 + Lapse0*dPhis(ii,jj))**Expon0)
        else
          Ps_ref(ii,jj) = Ps_avg(ii,jj)*exp(-dPhis(ii,jj)/RT_avg(ii,jj,nlev,ie))
        endif
      end do
      end do

      ! Calculate reference dp3d values
      !---------------------------------
      dlt_Ps(:,:) = Ps_ref(:,:) - Ps_avg(:,:)
      do kk=1,nlev
        dp3d_ref(:,:,kk,ie) = DP_avg(:,:,kk) + (hvcoord%hybi(kk+1)            &
                                               -hvcoord%hybi(kk  ))*dlt_Ps(:,:)
      end do

    end do ! ie=nets,nete

    ! End Routine
    !------------
    return
  end subroutine calc_dp3d_reference
  !=============================================================================


  !=============================================================================
  subroutine fill_element(Eval)
    !
    ! fill_element_bilin: Fill in element gridpoints using local bi-linear
    !                     interpolation of nearby average values.
    ! 
    !                     NOTE: This routine is hard coded for NP=4, if a
    !                           different value of NP is used... bad things 
    !                           will happen.
    !=======================================================================
    use dimensions_mod,only: np
    !                  
    ! Passed variables
    !-------------------
    real(kind=r8),intent(inout):: Eval(np,np)
    !
    ! Local Values
    !--------------
    real(kind=r8):: X0
    real(kind=r8):: S1,S2,S3,S4
    real(kind=r8):: C1,C2,C3,C4
    real(kind=r8):: E1,E2,E3,E4,E0

    X0 = sqrt(1._r8/5._r8)

    ! Set the "known" values Eval
    !----------------------------
    S1 = (Eval(1 ,2 )+Eval(1 ,3 ))/2._r8
    S2 = (Eval(2 ,np)+Eval(3 ,np))/2._r8
    S3 = (Eval(np,2 )+Eval(np,3 ))/2._r8
    S4 = (Eval(2 ,1 )+Eval(3 ,1 ))/2._r8
    C1 = Eval(1 ,1 )
    C2 = Eval(1 ,np)
    C3 = Eval(np,np)
    C4 = Eval(np,1 )

    ! E0 OPTION: Element Center value:
    !---------------------------------
    IF(.FALSE.) THEN
      ! Use ELEMENT AVERAGE value contained in (2,2)
      !----------------------------------------------
      E0 = Eval(2,2)
    ELSE
      ! Use AVG OF SIDE VALUES after boundary exchange of E0 (smooting option)
      !-----------------------------------------------------------------------
      E0 = (S1 + S2 + S3 + S4)/4._r8  
    ENDIF

    ! Calc interior values along center axes
    !----------------------------------------
    E1 = E0 + X0*(S1-E0)
    E2 = E0 + X0*(S2-E0)
    E3 = E0 + X0*(S3-E0)
    E4 = E0 + X0*(S4-E0)

    ! Calculate Side Gridpoint Values for Eval
    !------------------------------------------
    Eval(1 ,2 ) = S1 + X0*(C1-S1)
    Eval(1 ,3 ) = S1 + X0*(C2-S1)
    Eval(2 ,np) = S2 + X0*(C2-S2)
    Eval(3 ,np) = S2 + X0*(C3-S2)
    Eval(np,2 ) = S3 + X0*(C4-S3)
    Eval(np,3 ) = S3 + X0*(C3-S3)
    Eval(2 ,1 ) = S4 + X0*(C1-S4)
    Eval(3 ,1 ) = S4 + X0*(C4-S4)

    ! Calculate interior values
    !---------------------------
    Eval(2 ,2 ) = E1 + X0*(Eval(2 ,1 )-E1)
    Eval(2 ,3 ) = E1 + X0*(Eval(2 ,np)-E1)
    Eval(3 ,2 ) = E3 + X0*(Eval(3 ,1 )-E3)
    Eval(3 ,3 ) = E3 + X0*(Eval(3 ,np)-E3)

    ! End Routine
    !------------
    return
  end subroutine fill_element

end module prim_advance_mod
