#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "omp_config.h"
#define NEWEULER_B4B 1
#define OVERLAP 1

!SUBROUTINES:
!   prim_advec_tracers_remap_rk2()
!      SEM 2D RK2 + monotone remap + hyper viscosity
!      SEM 2D RK2 can use sign-preserving or monotone reconstruction
!
!For RK2 advection of Q:  (example of 2 stage RK for tracers):   dtq = qsplit*dt
!For consistency, if Q=1
!  dp1  = dp(t)- dtq div[ U1 dp(t)]
!  dp2  = dp1  - dtq div[ U2 dp1  ]  + 2*dtq D( dpdiss_ave )
!  dp*  = (dp(t) + dp2 )/2
!       =  dp(t) - dtq  div[ U1 dp(t) + U2 dp1 ]/2   + dtq D( dpdiss_ave )
!
!so we require:
!  U1 = Udp_ave / dp(t)
!  U2 = Udp_ave / dp1
!
!For tracer advection:
!  Qdp1  = Qdp(t)- dtq div[ U1 Qdp(t)]
!  Qdp2  = Qdp1  - dtq div[ U2 Qdp1  ]  + 2*dtq D( Q dpdiss_ave )
!  Qdp*  = (Qdp(t) + Qdp2 )/2
!       =  Qdp(t) - dtq  div[ U1 Qdp(t) + U2 Qdp1 ]   + dtq D( Q dpdiss_ave )
!
!Qdp1:  limit Q, with Q = Qdp1-before-DSS/(dp1-before-DSS)      with dp1 as computed above
!Qdp2:  limit Q, with Q = Qdp2-before-DSS/(dp2-before-DSS)      with dp2 as computed above
!
!For dissipation: Q = Qdp1-after-DSS / dp1-after-DSS
!
!
!last step:
!  remap Qdp* to Qdp(t+1)   [ dp_star(t+1) -> dp(t+1) ]


module prim_advection_mod_base
!
! two formulations.  both are conservative
! u grad Q formulation:
!
!    d/dt[ Q] +  U grad Q  +  eta_dot dp/dn dQ/dp  = 0
!                            ( eta_dot dQ/dn )
!
!    d/dt[ dp/dn ] = div( dp/dn U ) + d/dn ( eta_dot dp/dn )
!
! total divergence formulation:
!    d/dt[dp/dn Q] +  div( U dp/dn Q ) + d/dn ( eta_dot dp/dn Q ) = 0
!
! for convience, rewrite this as dp Q:  (since dn does not depend on time or the horizonal):
! equation is now:
!    d/dt[dp Q] +  div( U dp Q ) + d( eta_dot_dpdn Q ) = 0
!
!
  use kinds, only              : real_kind
  use dimensions_mod, only     : nlev, nlevp, np, qsize, ntrac, nc
  use physical_constants, only : rgas, Rwater_vapor, kappa, g, rearth, rrearth, cp
  use derivative_mod, only     : gradient, vorticity, gradient_wk, derivative_t, divergence, &
                                 gradient_sphere, divergence_sphere
  use element_mod, only        : element_t
  use fvm_control_volume_mod, only        : fvm_struct
  use filter_mod, only         : filter_t, filter_P
  use hybvcoord_mod, only      : hvcoord_t
  use time_mod, only           : TimeLevel_t, smooth, TimeLevel_Qdp
  use prim_si_mod, only        : preq_pressure
  use diffusion_mod, only      : scalar_diffusion, diffusion_init
  use control_mod, only        : integration, test_case, filter_freq_advection,  hypervis_order, &
        statefreq, moisture, TRACERADV_TOTAL_DIVERGENCE, TRACERADV_UGRADQ, &
        nu_q, nu_p, limiter_option, hypervis_subcycle_q, rsplit
  use edge_mod, only           : edgevpack, edgerotate, edgevunpack, initedgebuffer, initedgesbuffer, &
        edgevunpackmin, initghostbuffer3D
 
  use edgetype_mod, only       : EdgeDescriptor_t, EdgeBuffer_t, ghostbuffer3D_t
  use hybrid_mod, only         : hybrid_t
  use bndry_mod, only          : bndry_exchangev
  use viscosity_mod, only      : biharmonic_wk_scalar, biharmonic_wk_scalar_minmax, neighbor_minmax, &
                                 neighbor_minmax_start, neighbor_minmax_finish
  use perf_mod, only           : t_startf, t_stopf, t_barrierf ! _EXTERNAL
  use parallel_mod, only   : abortmp

  implicit none

  private
  save

  public :: Prim_Advec_Init1, Prim_Advec_Init2, prim_advec_init_deriv
  public :: Prim_Advec_Tracers_remap, Prim_Advec_Tracers_remap_rk2, Prim_Advec_Tracers_remap_ALE
  public :: prim_advec_tracers_fvm
  public :: vertical_remap

  type (EdgeBuffer_t)      :: edgeAdv, edgeAdvp1, edgeAdvQminmax, edgeAdv1,  edgeveloc
  type (ghostBuffer3D_t)   :: ghostbuf_tr

  integer,parameter :: DSSeta = 1
  integer,parameter :: DSSomega = 2
  integer,parameter :: DSSdiv_vdp_ave = 3
  integer,parameter :: DSSno_var = -1

  real(kind=real_kind), allocatable :: qmin(:,:,:), qmax(:,:,:)

  type (derivative_t), public, allocatable   :: deriv(:) ! derivative struct (nthreads)

contains

  subroutine Prim_Advec_Init1(par, elem, n_domains)
    use dimensions_mod, only : nlev, qsize, nelemd
    use parallel_mod, only : parallel_t
    use control_mod, only : use_semi_lagrange_transport
    use interpolate_mod,        only : interpolate_tracers_init
    type(parallel_t) :: par
    integer, intent(in) :: n_domains
    type (element_t) :: elem(:)
    type (EdgeDescriptor_t), allocatable :: desc(:)


    integer :: ie
    ! Shared buffer pointers.
    ! Using "=> null()" in a subroutine is usually bad, because it makes
    ! the variable have an implicit "save", and therefore shared between
    ! threads. But in this case we want shared pointers.
    real(kind=real_kind), pointer :: buf_ptr(:) => null()
    real(kind=real_kind), pointer :: receive_ptr(:) => null()


    ! this might be called with qsize=0
    ! allocate largest one first
    ! Currently this is never freed. If it was, only this first one should
    ! be freed, as only it knows the true size of the buffer.
    call initEdgeBuffer(par,edgeAdvp1,elem,qsize*nlev + nlev)
    call initEdgeBuffer(par,edgeAdv,elem,qsize*nlev)
    call initEdgeBuffer(par,edgeAdv1,elem,nlev)
    call initEdgeBuffer(par,edgeveloc,elem,2*nlev)

    ! This is a different type of buffer pointer allocation 
    ! used for determine the minimum and maximum value from 
    ! neighboring  elements
    call initEdgeSBuffer(par,edgeAdvQminmax,elem,qsize*nlev*2)

    ! Don't actually want these saved, if this is ever called twice.
    nullify(buf_ptr)
    nullify(receive_ptr)

    allocate(deriv(0:n_domains-1))

    ! this static array is shared by all threads, so dimension for all threads (nelemd), not nets:nete:
    allocate (qmin(nlev,qsize,nelemd))
    allocate (qmax(nlev,qsize,nelemd))

    if  (use_semi_lagrange_transport) then
       call initghostbuffer3D(ghostbuf_tr,nlev*qsize,np)
       call interpolate_tracers_init()
    endif

  end subroutine Prim_Advec_Init1

  subroutine Prim_Advec_Init_deriv(hybrid,fvm_corners, fvm_points)

    use kinds,          only : longdouble_kind
    use dimensions_mod, only : nc
    use derivative_mod, only : derivinit
    implicit none
    type (hybrid_t), intent(in) :: hybrid
    real(kind=longdouble_kind), intent(in) :: fvm_corners(nc+1)
    real(kind=longdouble_kind), intent(in) :: fvm_points(nc)

    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call derivinit(deriv(hybrid%ithr),fvm_corners, fvm_points)
  end subroutine Prim_Advec_Init_deriv

  subroutine Prim_Advec_Init2(elem,hvcoord,hybrid)
    use element_mod   , only : element_t
    use hybvcoord_mod , only : hvcoord_t
    implicit none
    type(element_t)   , intent(in) :: elem(:)
    type(hvcoord_t)   , intent(in) :: hvcoord
    type (hybrid_t)   , intent(in) :: hybrid
    !Nothing to do
  end subroutine Prim_Advec_Init2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fvm driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Prim_Advec_Tracers_fvm(elem, fvm, deriv,hvcoord,hybrid,&
        dt,tl,nets,nete)
    use perf_mod, only : t_startf, t_stopf            ! _EXTERNAL
    use vertremap_mod, only: remap1_nofilter  ! _EXTERNAL (actually INTERNAL)
    use fvm_mod, only : cslam_runairdensity, cslam_runflux, edgeveloc
    use fvm_mod, only: fvm_mcgregor, fvm_mcgregordss, fvm_rkdss
    use fvm_mod, only : fvm_ideal_test, IDEAL_TEST_OFF, IDEAL_TEST_ANALYTICAL_WINDS
    use fvm_mod, only : fvm_test_type, IDEAL_TEST_BOOMERANG, IDEAL_TEST_SOLIDBODY
    use fvm_bsp_mod, only: get_boomerang_velocities_gll, get_solidbody_velocities_gll
    use fvm_control_volume_mod, only : n0_fvm,np1_fvm,fvm_supercycling
    use control_mod, only : tracer_transport_type
    use control_mod, only : TRACERTRANSPORT_LAGRANGIAN_FVM, TRACERTRANSPORT_FLUXFORM_FVM
    use time_mod,    only : time_at

    implicit none
    type (element_t), intent(inout)   :: elem(:)
    type (fvm_struct), intent(inout)   :: fvm(:)
    type (derivative_t), intent(in)   :: deriv
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t),     intent(in):: hybrid
    type (TimeLevel_t)                :: tl

    real(kind=real_kind) , intent(in) :: dt
    integer,intent(in)                :: nets,nete


    real (kind=real_kind), dimension(np,np,nlev)    :: dp_star
    real (kind=real_kind), dimension(np,np,nlev)    :: dp
    real (kind=real_kind)                           :: eta_dot_dpdn(np,np,nlevp) 
    integer :: np1,ie,k

    real (kind=real_kind)  :: vstar(np,np,2)
    real (kind=real_kind)  :: vhat(np,np,2)
    real (kind=real_kind), dimension(np, np) :: v1, v2


    call t_barrierf('sync_prim_advec_tracers_fvm', hybrid%par%comm)
    call t_startf('prim_advec_tracers_fvm')
    np1 = tl%np1

    ! departure algorithm requires two velocities:
    !
    ! fvm%v0:        velocity at beginning of tracer timestep (time n0_qdp)
    !                this was saved before the (possibly many) dynamics steps
    ! elem%derived%vstar:    
    !                velocity at end of tracer timestep (time np1 = np1_qdp)
    !                for lagrangian dynamics, this is on lagrangian levels
    !                for eulerian dynamcis, this is on reference levels
    !                and it should be interpolated.
    !
    do ie=nets,nete
       elem(ie)%derived%vstar(:,:,:,:)=elem(ie)%state%v(:,:,:,:,np1)
    enddo


    if (rsplit==0) then
       ! interpolate t+1 velocity from reference levels to lagrangian levels
       ! For rsplit=0, we need to first compute lagrangian levels based on vertical velocity
       ! which requires we first DSS mean vertical velocity from dynamics
       ! note: we introduce a local eta_dot_dpdn() variable instead of DSSing elem%eta_dot_dpdn
       ! so as to preserve BFB results in some HOMME regression tests
       do ie=nets,nete
          do k=1,nlevp
             eta_dot_dpdn(:,:,k) = elem(ie)%derived%eta_dot_dpdn(:,:,k)*elem(ie)%spheremp(:,:)
          enddo
          ! eta_dot_dpdn at nlevp is zero, so we dont boundary exchange it:
          call edgeVpack(edgeAdv1,eta_dot_dpdn(:,:,1:nlev),nlev,0,ie)
       enddo

       call t_startf('pat_fvm_bexchV')
       call bndry_exchangeV(hybrid,edgeAdv1)
       call t_stopf('pat_fvm_bexchV')

       do ie=nets,nete
          ! restor interior values.  we could avoid this if we created a global array for eta_dot_dpdn
          do k=1,nlevp
             eta_dot_dpdn(:,:,k) = elem(ie)%derived%eta_dot_dpdn(:,:,k)*elem(ie)%spheremp(:,:)
          enddo
          ! unpack DSSed edge data
          call edgeVunpack(edgeAdv1,eta_dot_dpdn(:,:,1:nlev),nlev,0,ie)
          do k=1,nlevp
             eta_dot_dpdn(:,:,k) = eta_dot_dpdn(:,:,k)*elem(ie)%rspheremp(:,:)
          enddo


          ! SET VERTICAL VELOCITY TO ZERO FOR DEBUGGING
          !        elem(ie)%derived%eta_dot_dpdn(:,:,:)=0
          ! elem%state%u(np1)  = velocity at time t+1 on reference levels
          ! elem%derived%vstar = velocity at t+1 on floating levels (computed below)
!           call remap_UV_ref2lagrange(np1,dt,elem,hvcoord,ie)
          do k=1,nlev
             dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                  ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
             dp_star(:,:,k) = dp(:,:,k) + dt*(eta_dot_dpdn(:,:,k+1) -   eta_dot_dpdn(:,:,k))
             if (fvm_ideal_test == IDEAL_TEST_ANALYTICAL_WINDS) then
               if (fvm_test_type == IDEAL_TEST_BOOMERANG) then
                 elem(ie)%derived%vstar(:,:,:,k)=get_boomerang_velocities_gll(elem(ie), time_at(np1))
               else if (fvm_test_type == IDEAL_TEST_SOLIDBODY) then
                 elem(ie)%derived%vstar(:,:,:,k)=get_solidbody_velocities_gll(elem(ie), time_at(np1))
               else
                 call abortmp('Bad fvm_test_type in prim_step')
               end if
             end if
          enddo
          call remap1_nofilter(elem(ie)%derived%vstar,np,1,dp,dp_star)
       end do
    else
       ! do nothing
       ! for rsplit>0:  dynamics is also vertically Lagrangian, so we do not need 
       ! to interpolate v(np1). 
    endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2D advection step
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------------
!     call t_startf('fvm_depalg')

!     call fvm_mcgregordss(elem,fvm,nets,nete, hybrid, deriv, dt*fvm_supercycling, 3)
    call fvm_rkdss(elem,fvm,nets,nete, hybrid, deriv, dt*fvm_supercycling, 3)
    write(*,*) "fvm_rkdss dt ",dt*fvm_supercycling
!     call t_stopf('fvm_depalg')

!------------------------------------------------------------------------------------
    
    ! fvm departure calcluation should use vstar.
    ! from c(n0) compute c(np1):
    if (tracer_transport_type == TRACERTRANSPORT_FLUXFORM_FVM) then
      call cslam_runflux(elem,fvm,hybrid,deriv,dt,tl,nets,nete,hvcoord%hyai(1)*hvcoord%ps0)
    else if (tracer_transport_type == TRACERTRANSPORT_LAGRANGIAN_FVM) then
      call cslam_runairdensity(elem,fvm,hybrid,deriv,dt,tl,nets,nete,hvcoord%hyai(1)*hvcoord%ps0)
    else
      call abortmp('Bad tracer_transport_type in Prim_Advec_Tracers_fvm')
    end if

    call t_stopf('prim_advec_tracers_fvm')
  end subroutine Prim_Advec_Tracers_fvm



!=================================================================================================!

  subroutine Prim_Advec_Tracers_remap( elem , deriv , hvcoord , flt , hybrid , dt , tl , nets , nete )
    use control_mod   , only : use_semi_lagrange_transport
    implicit none
    type (element_t)     , intent(inout) :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (filter_t)      , intent(in   ) :: flt
    type (hybrid_t)      , intent(in   ) :: hybrid
    real(kind=real_kind) , intent(in   ) :: dt
    type (TimeLevel_t)   , intent(inout) :: tl
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete


  if (use_semi_lagrange_transport) then
    call Prim_Advec_Tracers_remap_ALE( elem , deriv ,                 hybrid , dt , tl , nets , nete )

  else
    call Prim_Advec_Tracers_remap_rk2( elem , deriv , hvcoord , flt , hybrid , dt , tl , nets , nete )
  end if
  end subroutine Prim_Advec_Tracers_remap



subroutine  Prim_Advec_Tracers_remap_ALE( elem , deriv , hybrid , dt , tl , nets , nete )
  use coordinate_systems_mod, only : cartesian3D_t, cartesian2D_t
  use dimensions_mod,         only : max_neigh_edges
  use edge_mod,               only : initghostbuffer3D, ghostVpack_unoriented, ghostVunpack_unoriented
  use bndry_mod,              only : ghost_exchangevfull
  use interpolate_mod,        only : interpolate_tracers, minmax_tracers
  use control_mod   ,         only : qsplit
  use global_norms_mod,       only : wrap_repro_sum
  use parallel_mod,           only: global_shared_buf, global_shared_sum, syncmp
  use control_mod,            only : use_semi_lagrange_transport_local_conservation



  implicit none
  type (element_t)     , intent(inout) :: elem(:)
  type (derivative_t)  , intent(in   ) :: deriv
  type (hybrid_t)      , intent(in   ) :: hybrid
  real(kind=real_kind) , intent(in   ) :: dt
  type (TimeLevel_t)   , intent(in   ) :: tl
  integer              , intent(in   ) :: nets
  integer              , intent(in   ) :: nete


  type(cartesian3D_t)                           :: dep_points  (np,np)
  integer                                       :: elem_indexes(np,np)
  type(cartesian2D_t)                           :: para_coords (np,np)
  real(kind=real_kind)                          :: Que         (np,np,nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: Que_t       (np,np,nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: minq        (np,np,nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: maxq        (np,np,nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: f                      (qsize)
  real(kind=real_kind)                          :: g                      (qsize)
  real(kind=real_kind)                          :: mass              (nlev,qsize)
  real(kind=real_kind)                          :: elem_mass         (nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: rho         (np,np,nlev,      nets:nete)

  real(kind=real_kind)                          :: neigh_q     (np,np,qsize,max_neigh_edges+1)
  real(kind=real_kind)                          :: u           (np,np,qsize)

  integer                                       :: i,j,k,l,n,q,ie,kptr, n0_qdp, np1_qdp
  integer                                       :: num_neighbors

  call t_barrierf('Prim_Advec_Tracers_remap_ALE', hybrid%par%comm)
  call t_startf('Prim_Advec_Tracers_remap_ALE')

  call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)

  ! compute displacements for departure grid
  ! store in elem%derived%vstar
  call ALE_RKdss (elem, nets, nete, hybrid, deriv, dt, tl)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  run ghost exchange to get global ID of all neighbors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  call t_startf('Prim_Advec_Tracers_remap_ALE_ghost_exchange')
  do ie=nets,nete
     kptr=0
     do k=1,nlev
     do q=1,qsize
        ! note: pack so that tracers per level are contiguous so we can unpack into
        ! array neigh_q()
        elem(ie)%state%Q(:,:,k,q) = elem(ie)%state%Qdp(:,:,k,q,n0_qdp) / elem(ie)%derived%dp(:,:,k)
        call ghostVpack_unoriented(ghostbuf_tr, elem(ie)%state%Q(:,:,k,q),np,1,kptr,elem(ie)%desc)
        kptr=kptr+1
     enddo
     enddo
  end do

  call t_startf('pat_remap_ale_gexchV')
  call ghost_exchangeVfull(hybrid%par,hybrid%ithr,ghostbuf_tr)
  call t_stopf('pat_remap_ale_gexchV')

  do ie=nets,nete
     num_neighbors = elem(ie)%desc%actual_neigh_edges+1
     do k=1,nlev

        ! find departure points
        call ALE_departure_from_gll     (dep_points, elem(ie)%derived%vstar(:,:,:,k), elem(ie), dt)

        ! find element containing departure point
        call ALE_elems_with_dep_points  (elem_indexes, dep_points, num_neighbors, elem(ie)%desc%neigh_corners)

        ! compute the parametric points
        call ALE_parametric_coords      (para_coords, elem_indexes, dep_points, num_neighbors, elem(ie)%desc%neigh_corners)

        ! for each level k, unpack all tracer neighbor data on that level
        kptr=(k-1)*qsize
        neigh_q=0
        u(:,:,:) = elem(ie)%state%Q(:,:,k,1:qsize)
        call ghostVunpack_unoriented (ghostbuf_tr, neigh_q, np, qsize, kptr, elem(ie)%desc, elem(ie)%GlobalId, u)

        do i=1,np
        do j=1,np
          ! interpolate tracers to deperature grid
          call interpolate_tracers     (para_coords(i,j), neigh_q(:,:,:,elem_indexes(i,j)),f)
          elem(ie)%state%Q(i,j,k,:) = f 
!         call minmax_tracers          (para_coords(i,j), neigh_q(:,:,:,elem_indexes(i,j)),f,g)
          do q=1,qsize
            f(q) = MINVAL(neigh_q(:,:,q,elem_indexes(i,j)))
            g(q) = MAXVAL(neigh_q(:,:,q,elem_indexes(i,j)))
          end do 
          minq(i,j,k,:,ie) = f 
          maxq(i,j,k,:,ie) = g 
        enddo
        enddo

     end do
  end do


  call t_stopf('Prim_Advec_Tracers_remap_ALE_ghost_exchange')
  ! compute original mass, at tl_1%n0
  elem_mass = 0
  do ie=nets,nete
    n=0
    do k=1,nlev
    do q=1,qsize
      n=n+1
      global_shared_buf(ie,n) = 0
      do j=1,np
        global_shared_buf(ie,n) = global_shared_buf(ie,n) + DOT_PRODUCT(elem(ie)%state%Qdp(:,j,k,q,n0_qdp),elem(ie)%spheremp(:,j))
        elem_mass(k,q,ie)       = elem_mass(k,q,ie)       + DOT_PRODUCT(elem(ie)%state%Qdp(:,j,k,q,n0_qdp),elem(ie)%spheremp(:,j))
      end do
    end do
    end do
  end do
  call wrap_repro_sum(nvars=n, comm=hybrid%par%comm)
  n=0
  do k=1,nlev
  do q=1,qsize
    n=n+1
    mass(k,q) = global_shared_sum(n)
  enddo
  enddo


  do ie=nets,nete
  do k=1,nlev
    rho(:,:,k,ie) = elem(ie)%spheremp(:,:)*elem(ie)%state%dp3d(:,:,k,tl%np1)
  end do
  end do

  do ie=nets,nete
    Que_t(:,:,:,1:qsize,ie) = elem(ie)%state%Q(:,:,:,1:qsize)
  end do
  call t_startf('Prim_Advec_Tracers_remap_ALE_Cobra')
  if (use_semi_lagrange_transport_local_conservation) then
    call Cobra_Elem (Que, Que_t, rho, minq, maxq, elem_mass, hybrid, nets, nete)
  else
    call Cobra_SLBQP(Que, Que_t, rho, minq, maxq, mass, hybrid, nets, nete)
  end if
  call t_stopf('Prim_Advec_Tracers_remap_ALE_Cobra')

  do ie=nets,nete
    elem(ie)%state%Q(:,:,:,1:qsize) =  Que(:,:,:,1:qsize,ie)
  end do


  do ie=nets,nete
  do k=1,nlev
  do q=1,qsize
     ! note: pack so that tracers per level are contiguous so we can unpack into
     ! array neigh_q()
     elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%state%Q(:,:,k,q) * elem(ie)%state%dp3d(:,:,k,tl%np1)
  enddo
  enddo
  end do


! do ie=nets,nete
!    do q = 1 , qsize
!       do k = 1 , nlev    !  Potential loop inversion (AAM)
!          elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%spheremp(:,:)* elem(ie)%state%Qdp(:,:,k,q,np1_qdp)
!       enddo
!    enddo
!     call edgeVpack(edgeAdv    , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
! enddo
! call bndry_exchangeV( hybrid , edgeAdv    )
! do ie = nets , nete
!    call edgeVunpack( edgeAdv    , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
!    do q = 1 , qsize
!       do k = 1 , nlev    !  Potential loop inversion (AAM)
!          elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,np1_qdp)
!       enddo
!    enddo
! enddo
  call t_stopf('Prim_Advec_Tracers_remap_ALE')


end subroutine Prim_Advec_Tracers_remap_ALE

subroutine VDOT(rp,Que,rho,mass,hybrid,nets,nete)
  use parallel_mod,        only: global_shared_buf, global_shared_sum
  use global_norms_mod,    only: wrap_repro_sum

  implicit none
  integer             , intent(in)              :: nets
  integer             , intent(in)              :: nete
  real(kind=real_kind), intent(out)             :: rp                (nlev,qsize)
  real(kind=real_kind), intent(in)              :: Que         (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: rho         (np*np,nlev,      nets:nete)
  real(kind=real_kind), intent(in)              :: mass              (nlev,qsize)
  type (hybrid_t)     , intent(in)              :: hybrid

  integer                                       :: k,n,q,ie

  global_shared_buf = 0 
  do ie=nets,nete
    n=0
    do q=1,qsize
    do k=1,nlev
      n=n+1
      global_shared_buf(ie,n) = global_shared_buf(ie,n) + DOT_PRODUCT(Que(:,k,q,ie), rho(:,k,ie))
    end do
    end do
  end do

  call wrap_repro_sum(nvars=n, comm=hybrid%par%comm)

  n=0
  do q=1,qsize
  do k=1,nlev
    n=n+1
    rp(k,q) = global_shared_sum(n) - mass(k,q)
  enddo
  enddo
  
end subroutine VDOT

subroutine Cobra_SLBQP(Que, Que_t, rho, minq, maxq, mass, hybrid, nets, nete) 

  use parallel_mod,        only: global_shared_buf, global_shared_sum
  use global_norms_mod,    only: wrap_repro_sum

  implicit none
  integer             , intent(in)              :: nets
  integer             , intent(in)              :: nete
  real(kind=real_kind), intent(out)             :: Que         (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: Que_t       (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: rho         (np*np,nlev,      nets:nete)
  real(kind=real_kind), intent(in)              :: minq        (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: maxq        (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: mass              (nlev,qsize)
  type (hybrid_t)     , intent(in)              :: hybrid

  integer,                            parameter :: max_clip = 50
  real(kind=real_kind),               parameter :: eta = 1D-08           
  real(kind=real_kind),               parameter :: hfd = 1D-10             
  real(kind=real_kind)                          :: lambda_p          (nlev,qsize)
  real(kind=real_kind)                          :: lambda_c          (nlev,qsize)
  real(kind=real_kind)                          :: rp                (nlev,qsize)
  real(kind=real_kind)                          :: rc                (nlev,qsize)
  real(kind=real_kind)                          :: rd                (nlev,qsize)
  real(kind=real_kind)                          :: alpha             (nlev,qsize)
  integer                                       :: j,k,n,q,ie
  integer                                       :: nclip

  nclip = 0

  Que(:,:,:,:) = Que_t(:,:,:,:)

  Que = MIN(MAX(Que,minq),maxq)

  call VDOT(rp,Que,rho,mass,hybrid,nets,nete)
  nclip = nclip + 1

  if (MAXVAL(ABS(rp)).lt.eta) return

  do ie=nets,nete
  do q=1,qsize
  do k=1,nlev
     Que(:,k,q,ie) = hfd * rho(:,k,ie) + Que_t(:,k,q,ie)
  enddo
  enddo
  enddo

  Que = MIN(MAX(Que,minq),maxq)

  call VDOT(rc,Que,rho,mass,hybrid,nets,nete)

  rd = rc-rp
  if (MAXVAL(ABS(rd)).eq.0) return 
  
  alpha = 0
  WHERE (rd.ne.0) alpha = hfd / rd 

  lambda_p = 0
  lambda_c =  -alpha*rp

  do while (MAXVAL(ABS(rc)).gt.eta .and. nclip.lt.max_clip)

    do ie=nets,nete
    do q=1,qsize
    do k=1,nlev
       Que(:,k,q,ie) = (lambda_c(k,q) + hfd) * rho(:,k,ie) + Que_t(:,k,q,ie)
    enddo
    enddo
    enddo
    Que = MIN(MAX(Que,minq),maxq)

    call VDOT(rc,Que,rho,mass,hybrid,nets,nete)
    nclip = nclip + 1

    rd = rp-rc

    if (MAXVAL(ABS(rd)).eq.0) exit

    alpha = 0
    WHERE (rd.ne.0) alpha = (lambda_p - lambda_c) / rd 

    rp       = rc
    lambda_p = lambda_c

    lambda_c = lambda_c -  alpha * rc

  enddo
end subroutine Cobra_SLBQP


subroutine Cobra_Elem(Que, Que_t, rho, minq, maxq, mass, hybrid, nets, nete) 

  use parallel_mod,        only: global_shared_buf, global_shared_sum
  use global_norms_mod,    only: wrap_repro_sum

  implicit none
  integer             , intent(in)              :: nets
  integer             , intent(in)              :: nete
  real(kind=real_kind), intent(out)             :: Que         (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: Que_t       (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: rho         (np*np,nlev,      nets:nete)
  real(kind=real_kind), intent(in)              :: minq        (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: maxq        (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: mass              (nlev,qsize,nets:nete)
  type (hybrid_t)     , intent(in)              :: hybrid

  integer,                            parameter :: max_clip = 50
  real(kind=real_kind),               parameter :: eta = 1D-10           
  real(kind=real_kind),               parameter :: hfd = 1D-08             
  real(kind=real_kind)                          :: lambda_p          (nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: lambda_c          (nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: rp                (nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: rc                (nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: rd                (nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: alpha             (nlev,qsize,nets:nete)
  integer                                       :: j,k,n,q,ie
  integer                                       :: nclip
  integer                                       :: mloc(3)

  nclip = 1

  Que(:,:,:,:) = Que_t(:,:,:,:)

  Que = MIN(MAX(Que,minq),maxq)

  do ie=nets,nete
  do q=1,qsize
  do k=1,nlev
    rp(k,q,ie) = DOT_PRODUCT(Que(:,k,q,ie), rho(:,k,ie)) - mass(k,q,ie)
  end do
  end do
  end do

  if (MAXVAL(ABS(rp)).lt.eta) return

  do ie=nets,nete
  do q=1,qsize
  do k=1,nlev
     Que(:,k,q,ie) = hfd * rho(:,k,ie) + Que_t(:,k,q,ie)
  enddo
  enddo
  enddo

  Que = MIN(MAX(Que,minq),maxq)

  do ie=nets,nete
  do q=1,qsize
  do k=1,nlev
    rc(k,q,ie) = DOT_PRODUCT(Que(:,k,q,ie), rho(:,k,ie)) - mass(k,q,ie)
  end do
  end do
  end do

  rd = rc-rp
  if (MAXVAL(ABS(rd)).eq.0) return 
  
  alpha = 0
  WHERE (rd.ne.0) alpha = hfd / rd 

  lambda_p = 0
  lambda_c =  -alpha*rp

! if (hybrid%par%masterproc) print *,__FILE__,__LINE__," mass(20,1,4):",mass(20,1,4)
! do k=1,np*np
!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," maxq(k,20,1,4):", &
!     maxq(k,20,1,4) ,minq(k,20,1,4),maxq(k,20,1,4)-minq(k,20,1,4)
! enddo
! do k=1,np*np
!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," Que(k,20,1,4):",Que(k,20,1,4) ,rho(k,20,4)
! enddo
  do while (MAXVAL(ABS(rc)).gt.eta .and. nclip.lt.max_clip)
    nclip = nclip + 1

    do ie=nets,nete
    do q=1,qsize
    do k=1,nlev
!      Que(:,k,q,ie) = (lambda_c(k,q,ie) + hfd) * rho(:,k,ie) + Que_t(:,k,q,ie)
       Que(:,k,q,ie) = lambda_c(k,q,ie) * rho(:,k,ie) + Que_t(:,k,q,ie)
    enddo
    enddo
    enddo

!   do ie=nets,nete
!   do q=1,qsize
!   do k=1,nlev
!     rc(k,q,ie) = DOT_PRODUCT(Que(:,k,q,ie), rho(:,k,ie)) - mass(k,q,ie)
!   end do
!   end do
!   end do
!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," rc(20,1,4):",rc(20,1,4), DOT_PRODUCT(Que(:,20,1,4), rho(:,20,4))

    Que = MIN(MAX(Que,minq),maxq)

    do ie=nets,nete
    do q=1,qsize
    do k=1,nlev
      rc(k,q,ie) = DOT_PRODUCT(Que(:,k,q,ie), rho(:,k,ie)) - mass(k,q,ie)
    end do
    end do
    end do

    
    mloc = MAXLOC(ABS(rc))
!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," MAXVAL(ABS(rc)):",MAXVAL(ABS(rc)), mloc, nclip

    rd = rp-rc

!   if (MAXVAL(ABS(rd)).eq.0) exit

    alpha = 0
    WHERE (rd.ne.0) alpha = (lambda_p - lambda_c) / rd 
!   WHERE (alpha.eq.0.and.MAXVAL(ABS(rc)).gt.eta) alpha=10;

    rp       = rc
    lambda_p = lambda_c

    lambda_c = lambda_c -  alpha * rc

!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," rc(20,1,4):",rc(20,1,4)
!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," rd(20,1,4):",rd(20,1,4)
!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," lambda_p(20,1,4):",lambda_p(20,1,4)
!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," lambda_c(20,1,4):",lambda_c(20,1,4)
!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," alpha(20,1,4):",alpha(20,1,4)
!   if (hybrid%par%masterproc) print *
  enddo
! if (hybrid%par%masterproc) print *,__FILE__,__LINE__," MAXVAL(ABS(rc)):",MAXVAL(ABS(rc)),eta," nclip:",nclip
end subroutine Cobra_Elem 

! ----------------------------------------------------------------------------------!
!SUBROUTINE ALE_RKDSS-----------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, MARK TAYLOR, 06. December 2012
!
! DESCRIPTION: ! create a runge kutta taylor serios mixture to calculate the departure grid
!
! CALLS:
! INPUT:
!
! OUTPUT:
!-----------------------------------------------------------------------------------!

! this will calculate the velocity at time t+1/2  along the trajectory s(t) given the velocities
! at the GLL points at time t and t+1 using a second order time accurate formulation.


subroutine ALE_RKdss(elem, nets, nete, hy, deriv, dt, tl)

  use derivative_mod,  only : derivative_t, ugradv_sphere
  use edge_mod,        only : edgevpack, edgevunpack
  use edgetype_mod,    only : EdgeBuffer_t
  use bndry_mod,       only : bndry_exchangev
  use kinds,           only : real_kind
  use hybrid_mod,      only : hybrid_t
  use element_mod,     only : element_t
  use dimensions_mod,   only : np, nlev

  implicit none

  type (element_t)     , intent(inout)             :: elem(:)
  integer              , intent(in   )             :: nets
  integer              , intent(in   )             :: nete
  type (hybrid_t)      , intent(in)                :: hy ! distributed parallel structure (shared)
  type (derivative_t)  , intent(in)                :: deriv ! derivative struct
  real (kind=real_kind), intent(in)                :: dt ! timestep
  type (TimeLevel_t)   , intent(in)                :: tl

  integer                                          :: ie, k
  real (kind=real_kind), dimension(np,np,2)        :: vtmp
  integer :: np1

  np1 = tl%np1


! RK-SSP 2 stage 2nd order:
!     x*(t+1) = x(t) + U(x(t),t) dt
!     x(t+1) = x(t) +  1/2 ( U(x*(t+1),t+1) + U(x(t),t) ) dt
! apply taylor series:
!  U(x*(t+1),t+1) = U(x(t),t+1) + (x*(t+1)-x(t)) gradU(x(t),t+1)
!
! x(t+1) = x(t) +  1/2 ( U(x(t),t+1) + (x*(t+1)-x(t)) gradU(x(t),t+1) + U(x(t),t) ) dt
! (x(t+1) - x(t)) / dt =  1/2 ( U(x(t),t+1) + (x*(t+1)-x(t)) gradU(x(t),t+1) + U(x(t),t) )
! (x(t+1) - x(t)) / dt =  1/2 ( U(x(t),t+1) + U(x(t),t) + (x*(t+1)-x(t)) gradU(x(t),t+1) )
! (x(t+1) - x(t)) / dt =  1/2 ( U(x(t),t+1) + U(x(t),t) + U(x(t),t) dt  gradU(x(t),t+1) )
!
!
!  (x(t+1)-x(t))/dt =  1/2(U(x(t),t+1) + U(x(t),t) + dt U(x(t),t) gradU(x(t),t+1))
!
! suppose dt = -ts (we go backward)
!  (x(t-ts)-x(t))/-ts =  1/2( U(x(t),t-ts)+U(x(t),t)) - ts 1/2 U(x(t),t) gradU(x(t),t-ts)
!
!  x(t-ts) = x(t)) -ts * [ 1/2( U(x(t),t-ts)+U(x(t),t)) - ts 1/2 U(x(t),t) gradU(x(t),t-ts) ]
!
!    !------------------------------------------------------------------------------------

  do ie=nets,nete
     ! vstarn0 = U(x,t)
     ! vstar   = U(x,t+1)
    do k=1,nlev
      vtmp(:,:,:)=ugradv_sphere(elem(ie)%state%v(:,:,:,k,np1), elem(ie)%derived%vstar(:,:,:,k),deriv,elem(ie))

      elem(ie)%derived%vstar(:,:,:,k) = &
           (elem(ie)%state%v(:,:,:,k,np1) + elem(ie)%derived%vstar(:,:,:,k))/2 - dt*vtmp(:,:,:)/2

      elem(ie)%derived%vstar(:,:,1,k) = elem(ie)%derived%vstar(:,:,1,k)*elem(ie)%spheremp(:,:)
      elem(ie)%derived%vstar(:,:,2,k) = elem(ie)%derived%vstar(:,:,2,k)*elem(ie)%spheremp(:,:)
    enddo
    call edgeVpack(edgeveloc,elem(ie)%derived%vstar,2*nlev,0,ie)
  enddo

  call t_startf('ALE_RKdss_bexchV')
  call bndry_exchangeV(hy,edgeveloc)
  call t_stopf('ALE_RKdss_bexchV')

  do ie=nets,nete
    call edgeVunpack(edgeveloc,elem(ie)%derived%vstar,2*nlev,0,ie)
    do k=1, nlev
      elem(ie)%derived%vstar(:,:,1,k) = elem(ie)%derived%vstar(:,:,1,k)*elem(ie)%rspheremp(:,:)
      elem(ie)%derived%vstar(:,:,2,k) = elem(ie)%derived%vstar(:,:,2,k)*elem(ie)%rspheremp(:,:)
    end do
  end do
end subroutine ALE_RKdss

! ----------------------------------------------------------------------------------!
!SUBROUTINE FVM_DEP_FROM_GLL----------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, MARK TAYLOR 14. December 2011                            !
! DESCRIPTION: calculates the deparute grid for fvm coming from the gll points      !
!                                                                                   !
! CALLS:
! INPUT:
!
! OUTPUT:
!-----------------------------------------------------------------------------------!
subroutine ALE_departure_from_gll(acart, vstar, elem, dt)
  use physical_constants,     only : rearth
  use coordinate_systems_mod, only : spherical_polar_t, cartesian3D_t, change_coordinates
  use time_mod,               only : timelevel_t
  use element_mod,            only : element_t
  use kinds,                  only : real_kind
  use dimensions_mod,         only : np

  implicit none

  type(cartesian3D_t)     ,intent(out)  :: acart(np,np)
  real (kind=real_kind)   ,intent(in)   :: vstar(np,np,2)
  type (element_t)        ,intent(in)   :: elem
  real (kind=real_kind)   ,intent(in)   :: dt

  integer                               :: i,j

  real (kind=real_kind)                 :: uxyz (np,np,3)

   ! convert velocity from lat/lon to cartesian 3D

  do i=1,3
     ! Summing along the third dimension is a sum over components for each point.
     ! (This is just a faster way of doing a dot product for each grid point,
     ! since reindexing the inputs to use the intrinsic effectively would be
     ! just asking for trouble.)
     uxyz(:,:,i)=sum( elem%vec_sphere2cart(:,:,i,:)*vstar(:,:,:) ,3)
  end do
  ! interpolate velocity to fvm nodes
  ! compute departure point
  ! crude, 1st order accurate approximation.  to be improved
  do i=1,np
     do j=1,np
        acart(i,j) = change_coordinates(elem%spherep(i,j)) 
        acart(i,j)%x = acart(i,j)%x - dt*uxyz(i,j,1)/rearth
        acart(i,j)%y = acart(i,j)%y - dt*uxyz(i,j,2)/rearth
        acart(i,j)%z = acart(i,j)%z - dt*uxyz(i,j,3)/rearth
     enddo
  enddo

end subroutine ALE_departure_from_gll




subroutine ALE_elems_with_dep_points (elem_indexes, dep_points, num_neighbors, ngh_corners)

  use element_mod,            only : element_t
  use dimensions_mod,         only : np
  use coordinate_systems_mod, only : cartesian3D_t, change_coordinates
  use interpolate_mod,        only : point_inside_quad

  implicit none

  ! The ngh_corners array is a list of corners of both elem and all of it's
  ! neighor elements all sorted by global id.
  integer              , intent(in)                :: num_neighbors
  type(cartesian3D_t),intent(in)                   :: ngh_corners(4,num_neighbors)
  integer              , intent(out)               :: elem_indexes(np,np)
  type(cartesian3D_t)  , intent(in)                :: dep_points(np,np)

  integer                                          :: i,j,n
  logical                                          :: inside

  elem_indexes = -1
  do i=1,np
    do j=1,np
! Just itererate the neighbors in global id order to get the same result on every processor.
      do n = 1, num_neighbors
! Mark Taylor's handy dandy point_inside_gc check.
       inside = point_inside_quad (ngh_corners(:,n), dep_points(i,j))
       if (inside) then
         elem_indexes(i,j) = n
         exit
       end if
     end do
    end do
  end do

  if (MINVAL(elem_indexes(:,:))==-1) then
    write (*,*) __FILE__,__LINE__,"Aborting because point not found in neighbor list. Info:"
    do i=1,np
      do j=1,np
        if (elem_indexes(i,j)==-1) then
          write (*,*)   " departure point ",dep_points(i,j)
          do n = 1, num_neighbors
            write (*,*) " quad checked    ",ngh_corners(1,n)
            write (*,*) "                 ",ngh_corners(2,n)
            write (*,*) "                 ",ngh_corners(3,n)
            write (*,*) "                 ",ngh_corners(4,n)
            write (*,*)
          end do
          exit
        end if
      end do
    end do
    call abortmp("ERROR elems_with_dep_points: Can't find departure grid. Time step too long?")
  end if
end subroutine ALE_elems_with_dep_points

function  shape_fcn_deriv(pc) result(dNds)
  real (kind=real_kind), intent(in)  ::  pc(2)
  real (kind=real_kind)              :: dNds(4,2)
 
  dNds(1, 1) = - 0.25 * (1.0 - pc(2))
  dNds(1, 2) = - 0.25 * (1.0 - pc(1))

  dNds(2, 1) =   0.25 * (1.0 - pc(2))
  dNds(2, 2) = - 0.25 * (1.0 + pc(1))

  dNds(3, 1) =   0.25 * (1.0 + pc(2))
  dNds(3, 2) =   0.25 * (1.0 + pc(1))

  dNds(4, 1) = - 0.25 * (1.0 + pc(2))
  dNds(4, 2) =   0.25 * (1.0 - pc(1))
end function   

function inv_2x2(A) result(A_inv)
  real (kind=real_kind), intent(in)  :: A    (2,2)
  real (kind=real_kind)              :: A_inv(2,2)
  real (kind=real_kind) :: det, denom

  det = A(1,1) * A(2,2) - A(2,1) * A(1,2)
  denom = 1/det
  ! inverse:
  A_inv(1,1) =  denom * A(2,2)  !  dxidx
  A_inv(2,1) = -denom * A(2,1)  !  detadx
  A_inv(1,2) = -denom * A(1,2)  !  dxidy
  A_inv(2,2) =  denom * A(1,1)  !  detady
end function

function INV(dxds) result(dsdx)

  real (kind=real_kind), intent(in)  :: dxds(3,2)

  real (kind=real_kind)  ::     dsdx(2,3)
  real (kind=real_kind)  ::      ata(2,2)
  real (kind=real_kind)  ::  ata_inv(2,2)


  !     dxds = | dxdxi   dxdeta |
  !            | dydxi   dydeta |
  !            | dzdxi   dzdeta |
  ata  = MATMUL(TRANSPOSE(dxds), dxds)
  ata_inv = inv_2x2(ata)
  dsdx = MATMUL(ata_inv, TRANSPOSE(dxds))
  !     dsdx = |  dxidx   dxidy   dxidz |
  !            | detadx  detady  detadz |

end function

subroutine shape_fcn(N, pc)
  real (kind=real_kind), intent(out) :: N(4)
  real (kind=real_kind), intent(in)  :: pc(2)

  ! shape function for each node evaluated at param_coords
  N(1) = 0.25 * (1.0 - pc(1)) * (1.0 - pc(2)) 
  N(2) = 0.25 * (1.0 + pc(1)) * (1.0 - pc(2)) 
  N(3) = 0.25 * (1.0 + pc(1)) * (1.0 + pc(2)) 
  N(4) = 0.25 * (1.0 - pc(1)) * (1.0 + pc(2)) 
end subroutine


function F(coords, pc) result(x)
  real (kind=real_kind), intent(in) :: pc(2), coords(4,3)

  real (kind=real_kind)            :: N(4), x(3)
  call shape_fcn(N,pc)
  x = MATMUL(TRANSPOSE(coords), N)
  x = x/SQRT(DOT_PRODUCT(x,x))
end function

function  DF(coords, pc) result(dxds)
  real (kind=real_kind), intent(in)  :: coords(4,3)
  real (kind=real_kind), intent(in)  :: pc(2)
 
  real (kind=real_kind)              :: dxds(3,2)
  real (kind=real_kind)              :: dNds(4,2)
  real (kind=real_kind)              ::  dds(3,2)
  real (kind=real_kind)              ::    c(2)
  real (kind=real_kind)              ::    x(3)
  real (kind=real_kind)              ::   xc(3,2)
  real (kind=real_kind)              :: nx, nx2 
  integer                            :: i,j

  dNds = shape_fcn_deriv  (pc)
  dds  = MATMUL(TRANSPOSE(coords), dNds)

  x    = F(coords, pc)
  nx2  = DOT_PRODUCT(x,x)
  nx   = SQRT(nx2)
  c    = MATMUL(TRANSPOSE(dds), x)
  do j=1,2
    do i=1,3
      xc(i,j) = x(i)*c(j)
    end do
  end do
  dxds = nx2*dds - xc
  dxds = dxds/(nx*nx2)
end function   


function cartesian_parametric_coordinates(sphere, corners3D) result (ref)
  use coordinate_systems_mod, only : cartesian2d_t, cartesian3D_t, spherical_polar_t, spherical_to_cart
  implicit none
  type (spherical_polar_t), intent(in) :: sphere
  type (cartesian3D_t)    , intent(in) :: corners3D(4)  !x,y,z coords of element corners

  type (cartesian2D_t)                 :: ref

  integer,               parameter :: MAXIT = 20
  real (kind=real_kind), parameter :: TOL   = 1.0E-13
  integer,               parameter :: n     = 3

  type (cartesian3D_t)             :: cart
  real (kind=real_kind)            :: coords(4,3), dxds(3,2), dsdx(2,3)
  real (kind=real_kind)            :: p(3), pc(2), dx(3), x(3), ds(2)
  real (kind=real_kind)            :: dist, step                          
  
  integer                          :: i,j,k,iter
  do i=1,4                               
    coords(i,1) = corners3D(i)%x 
    coords(i,2) = corners3D(i)%y 
    coords(i,3) = corners3D(i)%z 
  end do

  pc = 0
  p  = 0
  cart = spherical_to_cart(sphere)

  p(1) = cart%x
  p(2) = cart%y
  p(3) = cart%z 

  dx   = 0
  ds   = 0
  dsdx = 0
  dxds = 0

  !*-------------------------------------------------------------------------*!

  ! Initial guess, center of element
  dist = 9999999.                         
  step = 9999999.                         
  iter = 0 

  do while  (TOL*TOL.lt.dist .and. iter.lt.MAXIT .and. TOL*TOL.lt.step)
    iter = iter + 1

    dxds =  DF (coords, pc)
    x    =   F (coords, pc)
    dsdx = INV (dxds)

    dx   = x - p
    dist = DOT_PRODUCT(dx,dx)
    ds   = MATMUL(dsdx, dx)
    pc   = pc - ds
    step = DOT_PRODUCT(ds,ds)
  enddo

  ref%x = pc(1)
  ref%y = pc(2)
end function


subroutine  ALE_parametric_coords (parametric_coord, elem_indexes, dep_points, num_neighbors, ngh_corners)
  use coordinate_systems_mod, only : cartesian2d_t, cartesian3D_t, spherical_polar_t, change_coordinates, distance
  use interpolate_mod,        only : parametric_coordinates
  use dimensions_mod,         only : np

  implicit none

  type(cartesian2D_t)       , intent(out)       :: parametric_coord(np,np)
  type(cartesian3D_t)       , intent(in)        :: dep_points(np,np)
  integer                   , intent(in)        :: elem_indexes(np,np)
  integer                   , intent(in)        :: num_neighbors
  type(cartesian3D_t)       , intent(in)        :: ngh_corners(4,num_neighbors)

  type (spherical_polar_t)                      :: sphere(np,np)
  integer                                       :: i,j,n
  type(cartesian2D_t)                           :: parametric_test(np,np)
  real(kind=real_kind)                          :: d

  do j=1,np
    sphere(:,j) = change_coordinates(dep_points(:,j))
  end do

!!$  call t_startf('Prim_Advec_Tracers_remap_ALE_parametric_coords_cart')
!!$  do i=1,np
!!$    do j=1,np
!!$      n = elem_indexes(i,j)
!!$      ! Mark will fill in  parametric_coordinates for corners.
!!$      parametric_coord(i,j) = cartesian_parametric_coordinates(sphere(i,j),ngh_corners(:,n))
!!$    end do
!!$  end do
!!$  call t_stopf('Prim_Advec_Tracers_remap_ALE_parametric_coords_cart')
!  call t_startf('Prim_Advec_Tracers_remap_ALE_parametric_coords_dmap')
  do i=1,np
    do j=1,np
      n = elem_indexes(i,j)
      !parametric_test(i,j)  = parametric_coordinates(sphere(i,j),ngh_corners(:,n))
      parametric_coord(i,j)= parametric_coordinates(sphere(i,j),ngh_corners(:,n))
    end do
  end do
!  call t_stopf('Prim_Advec_Tracers_remap_ALE_parametric_coords_dmap')
! do i=1,np
!   do j=1,np
!     d = distance(parametric_coord(i,j),parametric_test(i,j))
!     print *,__LINE__,parametric_coord(i,j), parametric_test(i,j), d
!   end do
! end do
end subroutine ALE_parametric_coords


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! forward-in-time 2 level vertically lagrangian step
!  this code takes a lagrangian step in the horizontal
! (complete with DSS), and then applies a vertical remap
!
! This routine may use dynamics fields at timelevel np1
! In addition, other fields are required, which have to be
! explicitly saved by the dynamics:  (in elem(ie)%derived struct)
!
! Fields required from dynamics: (in
!    omega_p   it will be DSS'd here, for later use by CAM physics
!              we DSS omega here because it can be done for "free"
!    Consistent mass/tracer-mass advection (used if subcycling turned on)
!       dp()   dp at timelevel n0
!       vn0()  mean flux  < U dp > going from n0 to np1
!
! 3 stage
!    Euler step from t     -> t+.5
!    Euler step from t+.5  -> t+1.0
!    Euler step from t+1.0 -> t+1.5
!    u(t) = u(t)/3 + u(t+2)*2/3
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  subroutine Prim_Advec_Tracers_remap_rk2( elem , deriv , hvcoord , flt , hybrid , dt , tl , nets , nete )
    use perf_mod      , only : t_startf, t_stopf            ! _EXTERNAL
    use derivative_mod, only : divergence_sphere
    use control_mod   , only : vert_remap_q_alg, qsplit
    implicit none
    type (element_t)     , intent(inout) :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (filter_t)      , intent(in   ) :: flt
    type (hybrid_t)      , intent(in   ) :: hybrid
    real(kind=real_kind) , intent(in   ) :: dt
    type (TimeLevel_t)   , intent(inout) :: tl
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete

    integer :: i,j,k,l,ie,q,nmin
    integer :: nfilt,rkstage,rhs_multiplier
    integer :: n0_qdp, np1_qdp

    call t_barrierf('sync_prim_advec_tracers_remap_k2', hybrid%par%comm)
    call t_startf('prim_advec_tracers_remap_rk2')
!    call extrae_user_function(1)
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp) !time levels for qdp are not the same
    rkstage = 3 !   3 stage RKSSP scheme, with optimal SSP CFL

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! RK2 2D advection step
    ! note: stage 3 we take the oppertunity to DSS omega
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! use these for consistent advection (preserve Q=1)
    ! derived%vdp_ave        =  mean horiz. flux:   U*dp
    ! derived%eta_dot_dpdn    =  mean vertical velocity (used for remap)
    ! derived%omega_p         =  advection code will DSS this for the physics, but otherwise
    !                            it is not needed
    ! Also: save a copy of div(U dp) in derived%div(:,:,:,1), which will be DSS'd
    !       and a DSS'ed version stored in derived%div(:,:,:,2)

    call t_startf('precomput_divdp')
    call precompute_divdp( elem , hybrid , deriv , dt , nets , nete , n0_qdp )   
    call t_stopf('precomput_divdp')

    !rhs_multiplier is for obtaining dp_tracers at each stage:
    !dp_tracers(stage) = dp - rhs_multiplier*dt*divdp_proj

    call t_startf('euler_step_0')
    rhs_multiplier = 0
    call euler_step( np1_qdp , n0_qdp  , dt/2 , elem , hvcoord , hybrid , deriv , nets , nete , DSSdiv_vdp_ave , rhs_multiplier )
    call t_stopf('euler_step_0')

    call t_startf('euler_step_1')
    rhs_multiplier = 1
    call euler_step( np1_qdp , np1_qdp , dt/2 , elem , hvcoord , hybrid , deriv , nets , nete , DSSeta         , rhs_multiplier )
    call t_stopf('euler_step_1')

    call t_startf('euler_step_2')
    rhs_multiplier = 2
    call euler_step( np1_qdp , np1_qdp , dt/2 , elem , hvcoord , hybrid , deriv , nets , nete , DSSomega       , rhs_multiplier )
    call t_stopf('euler_step_2')

    !to finish the 2D advection step, we need to average the t and t+2 results to get a second order estimate for t+1.
    call t_startf('qdp_tavg')
    call qdp_time_avg( elem , rkstage , n0_qdp , np1_qdp , limiter_option , nu_p , nets , nete )
    call t_stopf('qdp_tavg')

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Dissipation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ( limiter_option == 8  ) then
      ! dissipation was applied in RHS.
    else
      call t_startf('ah_scalar')
      call advance_hypervis_scalar(edgeadv,elem,hvcoord,hybrid,deriv,tl%np1,np1_qdp,nets,nete,dt)
      call t_stopf('ah_scalar')
    endif
!    call extrae_user_function(0)

    call t_stopf('prim_advec_tracers_remap_rk2')

  end subroutine prim_advec_tracers_remap_rk2

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine precompute_divdp( elem , hybrid , deriv , dt , nets , nete , n0_qdp )
    implicit none
    type(element_t)      , intent(inout) :: elem(:)
    type (hybrid_t)      , intent(in   ) :: hybrid
    type (derivative_t)  , intent(in   ) :: deriv
    real(kind=real_kind) , intent(in   ) :: dt
    integer              , intent(in   ) :: nets , nete , n0_qdp
    integer :: ie , k

    do ie = nets , nete 
      do k = 1 , nlev   ! div( U dp Q),
        elem(ie)%derived%divdp(:,:,k) = divergence_sphere(elem(ie)%derived%vn0(:,:,:,k),deriv,elem(ie))
        elem(ie)%derived%divdp_proj(:,:,k) = elem(ie)%derived%divdp(:,:,k)
      enddo  
    enddo
  end subroutine precompute_divdp
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine qdp_time_avg( elem , rkstage , n0_qdp , np1_qdp , limiter_option , nu_p , nets , nete )
    implicit none
    type(element_t)     , intent(inout) :: elem(:)
    integer             , intent(in   ) :: rkstage , n0_qdp , np1_qdp , nets , nete , limiter_option
    real(kind=real_kind), intent(in   ) :: nu_p
    integer :: ie,q,k
    real(kind=real_kind) :: rrkstage
#ifdef NEWEULER_B4B
    do ie=nets,nete
      elem(ie)%state%Qdp(:,:,:,1:qsize,np1_qdp) =               &
                   ( elem(ie)%state%Qdp(:,:,:,1:qsize,n0_qdp) + &
                     (rkstage-1)*elem(ie)%state%Qdp(:,:,:,1:qsize,np1_qdp) ) / rkstage
    enddo
#else
    rrkstage=1.0d0/real(rkstage,kind=real_kind)
    do ie=nets,nete
      do q=1,qsize
        do k=1,nlev
           elem(ie)%state%Qdp(:,:,k,q,np1_qdp) =               &
               rrkstage *( elem(ie)%state%Qdp(:,:,k,q,n0_qdp) + &
               (rkstage-1)*elem(ie)%state%Qdp(:,:,k,q,np1_qdp) )
        enddo
      enddo
    enddo
#endif
  end subroutine qdp_time_avg

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
  ! ===================================
  ! This routine is the basic foward
  ! euler component used to construct RK SSP methods
  !
  !           u(np1) = u(n0) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! n0 can be the same as np1.
  !
  ! DSSopt = DSSeta or DSSomega:   also DSS eta_dot_dpdn or omega
  !
  ! ===================================
  use kinds          , only : real_kind
  use dimensions_mod , only : np, nlev
  use hybrid_mod     , only : hybrid_t
  use element_mod    , only : element_t
  use derivative_mod , only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere, limiter_optim_iter_full
  use edge_mod       , only : edgevpack, edgevunpack
  use bndry_mod      , only : bndry_exchangev
  use hybvcoord_mod  , only : hvcoord_t
  use parallel_mod, only : abortmp, iam
  implicit none
  integer              , intent(in   )         :: np1_qdp, n0_qdp
  real (kind=real_kind), intent(in   )         :: dt
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  integer              , intent(in   )         :: DSSopt
  integer              , intent(in   )         :: rhs_multiplier

  ! local
  real(kind=real_kind), dimension(np,np                       ) :: divdp, dpdiss
  real(kind=real_kind), dimension(np,np,nlev) :: dpdissk
  real(kind=real_kind), dimension(np,np,2                     ) :: gradQ
  real(kind=real_kind), dimension(np,np,2,nlev                ) :: Vstar
  real(kind=real_kind), dimension(np,np  ,nlev                ) :: Qtens
  real(kind=real_kind), dimension(np,np  ,nlev                ) :: dp,dp_star
  real(kind=real_kind), dimension(np,np  ,nlev,qsize,nets:nete) :: Qtens_biharmonic
  real(kind=real_kind), pointer, dimension(:,:,:)               :: DSSvar
  real(kind=real_kind) :: dp0(nlev)
  integer :: ie,q,i,j,k, kptr
  integer :: rhs_viss

!  call t_barrierf('sync_euler_step', hybrid%par%comm)
OMP_SIMD
  do k = 1 , nlev
    dp0(k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
  enddo
!pw  call t_startf('euler_step')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   compute Q min/max values for lim8
  !   compute biharmonic mixing term f
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rhs_viss = 0
  if ( limiter_option == 8  ) then
    call t_startf('bihmix_qminmax')
    ! when running lim8, we also need to limit the biharmonic, so that term needs
    ! to be included in each euler step.  three possible algorithms here:
    ! 1) most expensive:
    !     compute biharmonic (which also computes qmin/qmax) during all 3 stages
    !     be sure to set rhs_viss=1
    !     cost:  3 biharmonic steps with 3 DSS
    !
    ! 2) cheapest:
    !     compute biharmonic (which also computes qmin/qmax) only on first stage
    !     be sure to set rhs_viss=3
    !     reuse qmin/qmax for all following stages (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps with 1 DSS
    !     main concern:  viscosity
    !
    ! 3)  compromise:
    !     compute biharmonic (which also computes qmin/qmax) only on last stage
    !     be sure to set rhs_viss=3
    !     compute qmin/qmax directly on first stage
    !     reuse qmin/qmax for 2nd stage stage (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps, 2 DSS
    !
    !  NOTE  when nu_p=0 (no dissipation applied in dynamics to dp equation), we should
    !        apply dissipation to Q (not Qdp) to preserve Q=1
    !        i.e.  laplace(Qdp) ~  dp0 laplace(Q)
    !        for nu_p=nu_q>0, we need to apply dissipation to Q * diffusion_dp
    !
    ! initialize dp, and compute Q from Qdp (and store Q in Qtens_biharmonic)
    do ie = nets , nete
      ! add hyperviscosity to RHS.  apply to Q at timelevel n0, Qdp(n0)/dp
OMP_SIMD
      do k = 1 , nlev    !  Loop index added with implicit inversion (AAM)
        dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - rhs_multiplier*dt*elem(ie)%derived%divdp_proj(:,:,k)
      enddo
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k) collapse(2)
#endif
      do q = 1 , qsize
        do k=1,nlev
          Qtens_biharmonic(:,:,k,q,ie) = elem(ie)%state%Qdp(:,:,k,q,n0_qdp)/dp(:,:,k)
          if ( rhs_multiplier == 1 ) then
             ! for this stage, we skip neighbor_minmax() call, but update
             ! qmin/qmax with any new local extrema:
              qmin(k,q,ie)=min(qmin(k,q,ie),minval(Qtens_biharmonic(:,:,k,q,ie)))
              qmax(k,q,ie)=max(qmax(k,q,ie),maxval(Qtens_biharmonic(:,:,k,q,ie)))
          else
             ! for rhs_multiplier=0,2 we will call neighbor_minmax and compute
             ! the correct min/max values
              qmin(k,q,ie)=minval(Qtens_biharmonic(:,:,k,q,ie))
              qmax(k,q,ie)=maxval(Qtens_biharmonic(:,:,k,q,ie))
          endif
        enddo
      enddo
    enddo

    ! compute element qmin/qmax
    if ( rhs_multiplier == 0 ) then
      ! update qmin/qmax based on neighbor data for lim8
      call t_startf('eus_neighbor_minmax1')
      call neighbor_minmax(hybrid,edgeAdvQminmax,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
      call t_stopf('eus_neighbor_minmax1')
    endif

    ! get niew min/max values, and also compute biharmonic mixing term
    if ( rhs_multiplier == 2 ) then
      rhs_viss = 3
      ! two scalings depending on nu_p:
      ! nu_p=0:    qtens_biharmonic *= dp0                   (apply viscsoity only to q)
      ! nu_p>0):   qtens_biharmonc *= elem()%psdiss_ave      (for consistency, if nu_p=nu_q)
      if ( nu_p > 0 ) then
        do ie = nets , nete
#ifdef NEWEULER_B4B
#if (defined COLUMN_OPENMP)
       !$omp parallel do private(k, q) collapse(2)
#endif
          do k = 1 , nlev
            do q = 1 , qsize
              ! NOTE: divide by dp0 since we multiply by dp0 below
              Qtens_biharmonic(:,:,k,q,ie)=Qtens_biharmonic(:,:,k,q,ie)&
                *elem(ie)%derived%dpdiss_ave(:,:,k)/dp0(k)
            enddo
          enddo
#else
#if (defined COLUMN_OPENMP)
        !$omp parallel do private(k)
#endif
          do k = 1 , nlev
            ! NOTE: divide by dp0 since we multiply by dp0 below
            dpdissk(:,:,k) = elem(ie)%derived%dpdiss_ave(:,:,k)/dp0(k)
          enddo
#if (defined COLUMN_OPENMP)
        !$omp parallel do private(q,k) collapse(2)
#endif
          do q = 1 , qsize
            do k = 1 , nlev
              Qtens_biharmonic(:,:,k,q,ie)=Qtens_biharmonic(:,:,k,q,ie)*dpdissk(:,:,k)
            enddo
          enddo
#endif
        enddo ! ie loop
      endif ! nu_p > 0

!   Previous version of biharmonic_wk_scalar_minmax included a min/max
!   calculation into the boundary exchange.  This was causing cache issues.
!   Split the single operation into two separate calls 
!      call neighbor_minmax()
!      call biharmonic_wk_scalar() 
! 
!      call biharmonic_wk_scalar_minmax( elem , qtens_biharmonic , deriv , edgeAdvQ3 , hybrid , &
!           nets , nete , qmin(:,:,nets:nete) , qmax(:,:,nets:nete) )
#ifdef OVERLAP 
      call neighbor_minmax_start(hybrid,edgeAdvQminmax,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
      call biharmonic_wk_scalar(elem,qtens_biharmonic,deriv,edgeAdv,hybrid,nets,nete) 
      do ie = nets , nete
#if (defined COLUMN_OPENMP_notB4B)
!$omp parallel do private(k, q)
#endif
        do q = 1 , qsize
          do k = 1 , nlev    !  Loop inversion (AAM)
            ! note: biharmonic_wk() output has mass matrix already applied. Un-apply since we apply again below:
            qtens_biharmonic(:,:,k,q,ie) = &
                     -rhs_viss*dt*nu_q*dp0(k)*Qtens_biharmonic(:,:,k,q,ie) / elem(ie)%spheremp(:,:)
          enddo
        enddo
      enddo
      call neighbor_minmax_finish(hybrid,edgeAdvQminmax,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
#else
      call t_startf('eus_neighbor_minmax2')
      call neighbor_minmax(hybrid,edgeAdvQminmax,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
      call t_stopf('eus_neighbor_minmax2')
      call biharmonic_wk_scalar(elem,qtens_biharmonic,deriv,edgeAdv,hybrid,nets,nete) 

      do ie = nets , nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k, q) collapse(2)
#endif
        do q = 1 , qsize
          do k = 1 , nlev    !  Loop inversion (AAM)
            ! note: biharmonic_wk() output has mass matrix already applied. Un-apply since we apply again below:
            qtens_biharmonic(:,:,k,q,ie) = &
                     -rhs_viss*dt*nu_q*dp0(k)*Qtens_biharmonic(:,:,k,q,ie) / elem(ie)%spheremp(:,:)
          enddo
        enddo
      enddo
#endif

    endif
    call t_stopf('bihmix_qminmax')
  endif  ! compute biharmonic mixing term and qmin/qmax
  ! end of limiter_option == 8 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   2D Advection step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call t_startf('eus_2d_advec')
  do ie = nets , nete
    ! note: eta_dot_dpdn is actually dimension nlev+1, but nlev+1 data is
    ! all zero so we only have to DSS 1:nlev
    if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
    if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
    if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)

    ! Compute velocity used to advance Qdp
#if (defined COLUMN_OPENMP)
    !$omp parallel do private(k,q)
#endif
    do k = 1 , nlev    !  Loop index added (AAM)
      ! derived variable divdp_proj() (DSS'd version of divdp) will only be correct on 2nd and 3rd stage
      ! but that's ok because rhs_multiplier=0 on the first stage:
      dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - rhs_multiplier * dt * elem(ie)%derived%divdp_proj(:,:,k)
      Vstar(:,:,1,k) = elem(ie)%derived%vn0(:,:,1,k) / dp(:,:,k)
      Vstar(:,:,2,k) = elem(ie)%derived%vn0(:,:,2,k) / dp(:,:,k)

      if ( limiter_option == 8) then
        ! Note that the term dpdissk is independent of Q
        ! UN-DSS'ed dp at timelevel n0+1:
        dpdissk(:,:,k) = dp(:,:,k) - dt * elem(ie)%derived%divdp(:,:,k)
        if ( nu_p > 0 .and. rhs_viss /= 0 ) then
          ! add contribution from UN-DSS'ed PS dissipation
!          dpdiss(:,:) = ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) *
!          elem(ie)%derived%psdiss_biharmonic(:,:)
          dpdissk(:,:,k) = dpdissk(:,:,k) - rhs_viss * dt * nu_q &
                           * elem(ie)%derived%dpdiss_biharmonic(:,:,k) / elem(ie)%spheremp(:,:)
        endif
        ! IMPOSE ZERO THRESHOLD.  do this here so it can be turned off for
        ! testing
        do q=1,qsize
          qmin(k,q,ie)=max(qmin(k,q,ie),0d0)
        enddo
      endif  ! limiter == 8

      ! also DSS extra field
      DSSvar(:,:,k) = elem(ie)%spheremp(:,:) * DSSvar(:,:,k)
    enddo
    call edgeVpack( edgeAdvp1 , DSSvar(:,:,1:nlev) , nlev , nlev*qsize , ie)

    ! advance Qdp
#if (defined COLUMN_OPENMP)
 !$omp parallel do private(q,k,gradQ,dp_star,qtens)
#endif
    do q = 1 , qsize
      do k = 1 , nlev  !  dp_star used as temporary instead of divdp (AAM)
        ! div( U dp Q),
        gradQ(:,:,1) = Vstar(:,:,1,k) * elem(ie)%state%Qdp(:,:,k,q,n0_qdp)
        gradQ(:,:,2) = Vstar(:,:,2,k) * elem(ie)%state%Qdp(:,:,k,q,n0_qdp)
        dp_star(:,:,k) = divergence_sphere( gradQ , deriv , elem(ie) )
        Qtens(:,:,k) = elem(ie)%state%Qdp(:,:,k,q,n0_qdp) - dt * dp_star(:,:,k)
        ! optionally add in hyperviscosity computed above:
        if ( rhs_viss /= 0 ) Qtens(:,:,k) = Qtens(:,:,k) + Qtens_biharmonic(:,:,k,q,ie)
      enddo

      if ( limiter_option == 8) then
        ! apply limiter to Q = Qtens / dp_star
        call limiter_optim_iter_full( Qtens(:,:,:) , elem(ie)%spheremp(:,:) , qmin(:,q,ie) , &
                                      qmax(:,q,ie) , dpdissk )
      endif

      ! apply mass matrix, overwrite np1 with solution:
      ! dont do this earlier, since we allow np1_qdp == n0_qdp
      ! and we dont want to overwrite n0_qdp until we are done using it
      do k = 1 , nlev
        elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%spheremp(:,:) * Qtens(:,:,k)
      enddo

      if ( limiter_option == 4 ) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        ! sign-preserving limiter, applied after mass matrix
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        call limiter2d_zero( elem(ie)%state%Qdp(:,:,:,q,np1_qdp))
      endif

      call edgeVpack(edgeAdvp1 , elem(ie)%state%Qdp(:,:,:,q,np1_qdp) , nlev , nlev*(q-1) , ie )
    enddo
  enddo ! ie loop

  call t_startf('eus_bexchV')
  call bndry_exchangeV( hybrid , edgeAdvp1 )
  call t_stopf('eus_bexchV')

  do ie = nets , nete
    if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
    if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
    if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)

    call edgeVunpack( edgeAdvp1 , DSSvar(:,:,1:nlev) , nlev , qsize*nlev , ie )
OMP_SIMD
    do k = 1 , nlev
      DSSvar(:,:,k) = DSSvar(:,:,k) * elem(ie)%rspheremp(:,:)
    enddo

#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k)
#endif
    do q = 1 , qsize
      call edgeVunpack( edgeAdvp1 , elem(ie)%state%Qdp(:,:,:,q,np1_qdp) , nlev , nlev*(q-1) , ie )
      do k = 1 , nlev    !  Potential loop inversion (AAM)
        elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,np1_qdp)
      enddo
    enddo
  enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
  call t_stopf('eus_2d_advec')
!pw call t_stopf('euler_step')
  end subroutine euler_step
!-----------------------------------------------------------------------------



  subroutine limiter2d_zero(Q)
  ! mass conserving zero limiter (2D only).  to be called just before DSS
  !
  ! this routine is called inside a DSS loop, and so Q had already
  ! been multiplied by the mass matrix.  Thus dont include the mass
  ! matrix when computing the mass = integral of Q over the element
  !
  ! ps is only used when advecting Q instead of Qdp
  ! so ps should be at one timelevel behind Q
  implicit none
  real (kind=real_kind), intent(inout) :: Q(np,np,nlev)

  ! local
  real (kind=real_kind) :: dp(np,np)
  real (kind=real_kind) :: mass,mass_new,ml
  integer i,j,k

  do k = nlev , 1 , -1
    mass = 0
    do j = 1 , np
      do i = 1 , np
        !ml = Q(i,j,k)*dp(i,j)*spheremp(i,j)  ! see above
        ml = Q(i,j,k)
        mass = mass + ml
      enddo
    enddo

    ! negative mass.  so reduce all postive values to zero
    ! then increase negative values as much as possible
    if ( mass < 0 ) Q(:,:,k) = -Q(:,:,k)
    mass_new = 0
    do j = 1 , np
      do i = 1 , np
        if ( Q(i,j,k) < 0 ) then
          Q(i,j,k) = 0
        else
          ml = Q(i,j,k)
          mass_new = mass_new + ml
        endif
      enddo
    enddo

    ! now scale the all positive values to restore mass
    if ( mass_new > 0 ) Q(:,:,k) = Q(:,:,k) * abs(mass) / mass_new
    if ( mass     < 0 ) Q(:,:,k) = -Q(:,:,k)
  enddo
  end subroutine limiter2d_zero

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine advance_hypervis_scalar( edgeAdv , elem , hvcoord , hybrid , deriv , nt , nt_qdp , nets , nete , dt2 )
  !  hyperviscsoity operator for foward-in-time scheme
  !  take one timestep of:
  !          Q(:,:,:,np) = Q(:,:,:,np) +  dt2*nu*laplacian**order ( Q )
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  use kinds          , only : real_kind
  use dimensions_mod , only : np, nlev
  use hybrid_mod     , only : hybrid_t
  use element_mod    , only : element_t
  use derivative_mod , only : derivative_t
  use edge_mod       , only : edgevpack, edgevunpack
  use edgetype_mod   , only : EdgeBuffer_t
  use bndry_mod      , only : bndry_exchangev
  use perf_mod       , only : t_startf, t_stopf                          ! _EXTERNAL
  implicit none
  type (EdgeBuffer_t)  , intent(inout)         :: edgeAdv
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nt
  integer              , intent(in   )         :: nt_qdp
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  real (kind=real_kind), intent(in   )         :: dt2

  ! local
  real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: Qtens
  real (kind=real_kind), dimension(np,np,nlev                ) :: dp
  real (kind=real_kind) :: dt,dp0
  integer :: k , i , j , ie , ic , q

  if ( nu_q           == 0 ) return
  if ( hypervis_order /= 2 ) return
!   call t_barrierf('sync_advance_hypervis_scalar', hybrid%par%comm)
  call t_startf('advance_hypervis_scalar')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dt = dt2 / hypervis_subcycle_q

  do ic = 1 , hypervis_subcycle_q
    do ie = nets , nete
      ! Qtens = Q/dp   (apply hyperviscsoity to dp0 * Q, not Qdp)
      ! various options:
      !   1)  biharmonic( Qdp )
      !   2)  dp0 * biharmonic( Qdp/dp )
      !   3)  dpave * biharmonic(Q/dp)
      ! For trace mass / mass consistenciy, we use #2 when nu_p=0
      ! and #e when nu_p>0, where dpave is the mean mass flux from the nu_p
      ! contribution from dynamics.

      if (nu_p>0) then
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k) collapse(2)
#endif
        do q = 1 , qsize
          do k = 1 , nlev
            dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - dt2*elem(ie)%derived%divdp_proj(:,:,k)
            Qtens(:,:,k,q,ie) = elem(ie)%derived%dpdiss_ave(:,:,k)*&
                                elem(ie)%state%Qdp(:,:,k,q,nt_qdp) / dp(:,:,k)
          enddo
        enddo

      else
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k,dp0) collapse(2)
#endif
        do q = 1 , qsize
          do k = 1 , nlev
            dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - dt2*elem(ie)%derived%divdp_proj(:,:,k)
            dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) ) * hvcoord%ps0 + &
                  ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) * hvcoord%ps0
            Qtens(:,:,k,q,ie) = dp0*elem(ie)%state%Qdp(:,:,k,q,nt_qdp) / dp(:,:,k)
          enddo
        enddo
      endif
    enddo ! ie loop

    ! compute biharmonic operator. Qtens = input and output
    call biharmonic_wk_scalar( elem , Qtens , deriv , edgeAdv , hybrid , nets , nete )

    do ie = nets , nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k,j,i)
#endif
      do q = 1 , qsize
        do k = 1 , nlev
          do j = 1 , np
            do i = 1 , np
              ! advection Qdp.  For mass advection consistency:
              ! DIFF( Qdp) ~   dp0 DIFF (Q)  =  dp0 DIFF ( Qdp/dp )
              elem(ie)%state%Qdp(i,j,k,q,nt_qdp) = elem(ie)%state%Qdp(i,j,k,q,nt_qdp) * elem(ie)%spheremp(i,j) &
                                                   - dt * nu_q * Qtens(i,j,k,q,ie)
            enddo
          enddo
        enddo
        
        if (limiter_option .ne. 0 ) then
           ! smooth some of the negativities introduced by diffusion:
           call limiter2d_zero( elem(ie)%state%Qdp(:,:,:,q,nt_qdp) )
        endif

      enddo
      call edgeVpack  ( edgeAdv , elem(ie)%state%Qdp(:,:,:,:,nt_qdp) , qsize*nlev , 0 , ie )
    enddo ! ie loop

    call t_startf('ah_scalar_bexchV')
    call bndry_exchangeV( hybrid , edgeAdv )
    call t_stopf('ah_scalar_bexchV')

    do ie = nets , nete
      call edgeVunpack( edgeAdv , elem(ie)%state%Qdp(:,:,:,:,nt_qdp) , qsize*nlev , 0 , ie )
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k) collapse(2)
#endif
      do q = 1 , qsize
        ! apply inverse mass matrix
        do k = 1 , nlev
          elem(ie)%state%Qdp(:,:,k,q,nt_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,nt_qdp)
        enddo
      enddo
    enddo ! ie loop
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
  enddo
  call t_stopf('advance_hypervis_scalar')
  end subroutine advance_hypervis_scalar





  subroutine vertical_remap(hybrid,elem,fvm,hvcoord,dt,np1,np1_qdp,np1_fvm,nets,nete)

  ! This routine is called at the end of the vertically Lagrangian
  ! dynamics step to compute the vertical flux needed to get back
  ! to reference eta levels
  !
  ! input:
  !     derived%dp()  delta p on levels at beginning of timestep
  !     state%dp3d(np1)  delta p on levels at end of timestep
  ! output:
  !     state%ps_v(np1)          surface pressure at time np1
  !     derived%eta_dot_dpdn()   vertical flux from final Lagrangian
  !                              levels to reference eta levels
  !
  use kinds,          only: real_kind
  use hybvcoord_mod,  only: hvcoord_t
  use vertremap_mod,  only: remap1, remap1_nofilter, remap_q_ppm ! _EXTERNAL (actually INTERNAL)
  use control_mod,    only: rsplit, tracer_transport_type
  use parallel_mod,   only: abortmp
  use hybrid_mod,     only: hybrid_t
  use derivative_mod, only: interpolate_gll2fvm_points

  use fvm_control_volume_mod, only : fvm_struct

  type (hybrid_t),  intent(in)    :: hybrid  ! distributed parallel structure (shared)
  type(fvm_struct), intent(inout) :: fvm(:)
  real (kind=real_kind)           :: cdp(1:nc,1:nc,nlev,ntrac)
  real (kind=real_kind)           :: psc(nc,nc), dpc(nc,nc,nlev),dpc_star(nc,nc,nlev)
  type (element_t), intent(inout) :: elem(:)
  type (hvcoord_t)                :: hvcoord
  real (kind=real_kind)           :: dt

  integer :: ie,i,j,k,np1,nets,nete,np1_qdp,np1_fvm
  integer :: q

  real (kind=real_kind), dimension(np,np,nlev)  :: dp,dp_star
  real (kind=real_kind), dimension(np,np,nlev,2)  :: ttmp

  call t_startf('vertical_remap')

  ! reference levels:
  !   dp(k) = (hyai(k+1)-hyai(k))*ps0 + (hybi(k+1)-hybi(k))*ps_v(i,j)
  !   hybi(1)=0          pure pressure at top of atmosphere
  !   hyai(1)=ptop
  !   hyai(nlev+1) = 0   pure sigma at bottom
  !   hybi(nlev+1) = 1
  !
  ! sum over k=1,nlev
  !  sum(dp(k)) = (hyai(nlev+1)-hyai(1))*ps0 + (hybi(nlev+1)-hybi(1))*ps_v
  !             = -ps0 + ps_v
  !  ps_v =  ps0+sum(dp(k))
  !
  ! reference levels:
  !    dp(k) = (hyai(k+1)-hyai(k))*ps0 + (hybi(k+1)-hybi(k))*ps_v
  ! floating levels:
  !    dp_star(k) = dp(k) + dt_q*(eta_dot_dpdn(i,j,k+1) - eta_dot_dpdn(i,j,k) )
  ! hence:
  !    (dp_star(k)-dp(k))/dt_q = (eta_dot_dpdn(i,j,k+1) - eta_dot_dpdn(i,j,k) )
  !

  do ie=nets,nete
!     ! SET VERTICAL VELOCITY TO ZERO FOR DEBUGGING
!     elem(ie)%derived%eta_dot_dpdn(:,:,:)=0
     if (rsplit==0) then
        ! compute dp_star from eta_dot_dpdn():
        do k=1,nlev
           dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
           dp_star(:,:,k) = dp(:,:,k) + dt*(elem(ie)%derived%eta_dot_dpdn(:,:,k+1) -&
                elem(ie)%derived%eta_dot_dpdn(:,:,k))
        enddo
        if (minval(dp_star)<0) call abortmp('negative layer thickness.  timestep or remap time too large')
     else
        !  REMAP u,v,T from levels in dp3d() to REF levels
        !
        ! update final ps_v
        elem(ie)%state%ps_v(:,:,np1) = hvcoord%hyai(1)*hvcoord%ps0 + &
             sum(elem(ie)%state%dp3d(:,:,:,np1),3)
        do k=1,nlev
           dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
           dp_star(:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)
        enddo
        if (minval(dp_star)<0) call abortmp('negative layer thickness.  timestep or remap time too large')

        ! remap the dynamics:
#undef REMAP_TE
#ifdef REMAP_TE
        ! remap u,v and cp*T + .5 u^2
        ttmp(:,:,:,1)=(elem(ie)%state%v(:,:,1,:,np1)**2 + &
             elem(ie)%state%v(:,:,2,:,np1)**2)/2 + &
             elem(ie)%state%t(:,:,:,np1)*cp
#else
        ttmp(:,:,:,1)=elem(ie)%state%t(:,:,:,np1)
#endif
        ttmp(:,:,:,1)=ttmp(:,:,:,1)*dp_star

        call t_startf('vertical_remap1_1')
        call remap1(ttmp,np,1,dp_star,dp)
        call t_stopf('vertical_remap1_1')

        elem(ie)%state%t(:,:,:,np1)=ttmp(:,:,:,1)/dp

        ttmp(:,:,:,1)=elem(ie)%state%v(:,:,1,:,np1)*dp_star
        ttmp(:,:,:,2)=elem(ie)%state%v(:,:,2,:,np1)*dp_star

        call t_startf('vertical_remap1_2')
        call remap1(ttmp,np,2,dp_star,dp)
        call t_stopf('vertical_remap1_2')

        elem(ie)%state%v(:,:,1,:,np1)=ttmp(:,:,:,1)/dp
        elem(ie)%state%v(:,:,2,:,np1)=ttmp(:,:,:,2)/dp

#ifdef REMAP_TE
        ! back out T from TE
        elem(ie)%state%t(:,:,:,np1) = &
             ( elem(ie)%state%t(:,:,:,np1) - ( (elem(ie)%state%v(:,:,1,:,np1)**2 + &
             elem(ie)%state%v(:,:,2,:,np1)**2)/2))/cp
#endif
     endif

     ! remap the gll tracers from lagrangian levels (dp_star)  to REF levels dp
     if (qsize>0) then

       call t_startf('vertical_remap1_3')
       call remap1(elem(ie)%state%Qdp(:,:,:,:,np1_qdp),np,qsize,dp_star,dp)
       call t_stopf('vertical_remap1_3')

     endif

     if (ntrac>0) then

        ! create local variable  cdp(1:nc,1:nc,nlev,ntrac)
        ! cdp(:,:,:,n) = fvm%c(:,:,:,n,np1_fvm)*fvm%dp_fvm(:,:,:,np1_fvm)
        ! dp(:,:,:) = reference level thicknesses

        ! call remap1(cdp,nc,ntrac-1,fvm%c(:,:,:,1,np1_fvm),dp)

        ! convert back to mass:
        ! fvm%dp_fvm(:,:,:,np1_fvm) = dp(:,:,:) ??XXgoldyXX??
        ! fvm%c(:,:,:,n,np1_fvm) = fvm%c(:,:,:,n,np1_fvm)/dp(:,:,:)

        if (ntrac>0) then
           !
           ! Recompute dp_fvm (this will not be necessary when SE fluxes are coded)
           !
           do k = 1, nlev
              fvm(ie)%dp_fvm(1:nc,1:nc,k,np1_fvm)=interpolate_gll2fvm_points(dp(:,:,k),deriv(hybrid%ithr))
           end do
           !
           !
           !
           do i=1,nc
              do j=1,nc
                 !
                 ! compute surface pressure implied by fvm scheme: psC
                 psc(i,j)=sum(fvm(ie)%dp_fvm(i,j,:,np1_fvm)) +  hvcoord%hyai(1)*hvcoord%ps0
                 !
                 ! compute source (cdp) and target (dpc) pressure grids for vertical remapping
                 !
                 do k=1,nlev
                    dpc(i,j,k) = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + &
                         (hvcoord%hybi(k+1) - hvcoord%hybi(k))*psc(i,j)
                    cdp(i,j,k,1:ntrac)=fvm(ie)%c(i,j,k,1:ntrac,np1_fvm)*fvm(ie)%dp_fvm(i,j,k,np1_fvm)
                 end do
              end do
           end do
           dpc_star=fvm(ie)%dp_fvm(1:nc,1:nc,:,np1_fvm)

           call t_startf('vertical_remap1_5')
           call remap1(cdp,nc,ntrac,dpc_star,dpc)
           call t_stopf('vertical_remap1_5')

           do k=1,nlev
              do j=1,nc
                 do i=1,nc
                    fvm(ie)%dp_fvm(i,j,k,np1_fvm)=dpc(i,j,k) !!XXgoldyXX??
                    fvm(ie)%c(i,j,k,1:ntrac,np1_fvm)=cdp(i,j,k,1:ntrac)/dpc(i,j,k)
                 end do
              end do
           end do
        end if
        !         call remap_velocityC(np1,dt,elem,fvm,hvcoord,ie)
     endif

  enddo
  call t_stopf('vertical_remap')
  end subroutine vertical_remap

end module prim_advection_mod_base
