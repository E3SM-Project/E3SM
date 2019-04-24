#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

module sl_advection
  use kinds, only              : real_kind, int_kind
  use dimensions_mod, only     : nlev, nlevp, np, qsize, qsize_d
  use physical_constants, only : rgas, Rwater_vapor, kappa, g, rearth, rrearth, cp
  use derivative_mod, only     : derivative_t, gradient_sphere, divergence_sphere
  use element_mod, only        : element_t
  use hybvcoord_mod, only      : hvcoord_t
  use time_mod, only           : TimeLevel_t, TimeLevel_Qdp
  use control_mod, only        : integration, test_case, hypervis_order, transport_alg, limiter_option
  use edge_mod, only           : edgevpack_nlyr, edgevunpack_nlyr, edge_g, &
       initghostbuffer3D, ghostVpack_unoriented, ghostVunpack_unoriented
  use edgetype_mod, only       : EdgeDescriptor_t, EdgeBuffer_t, ghostbuffer3D_t
  use hybrid_mod, only         : hybrid_t
  use bndry_mod, only          : bndry_exchangev
  use perf_mod, only           : t_startf, t_stopf, t_barrierf ! _EXTERNAL
  use parallel_mod, only       : abortmp, parallel_t
  use coordinate_systems_mod, only : cartesian3D_t
  use compose_mod

  implicit none

  private

  type (ghostBuffer3D_t)   :: ghostbuf_tr
  integer :: sl_mpi
  type (cartesian3D_t), allocatable :: dep_points_all(:,:,:,:) ! (np,np,nlev,nelemd)

  public :: Prim_Advec_Tracers_remap_ALE, sl_init1

  logical, parameter :: barrier = .false.

contains

  !=================================================================================================!

  subroutine sl_parse_transport_alg(transport_alg, slmm, cisl, qos, sl_test)
    integer, intent(in) :: transport_alg
    logical, intent(out) :: slmm, cisl, qos, sl_test

    slmm = transport_alg > 1
    cisl = transport_alg == 2 .or. transport_alg == 3 .or. transport_alg >= 20
    qos  = cisl .and. (transport_alg == 3 .or. transport_alg == 39)  
    sl_test = (transport_alg >= 17 .and. transport_alg <= 19) .or. &
         transport_alg == 29 .or. transport_alg == 39
  end subroutine sl_parse_transport_alg

  subroutine sl_init1(par, elem)
    use interpolate_mod,        only : interpolate_tracers_init
    use control_mod,            only : transport_alg, semi_lagrange_cdr_alg, cubed_sphere_map, &
         nu_q, semi_lagrange_hv_q_all
    use element_state,          only : timelevels
    use coordinate_systems_mod, only : cartesian3D_t, change_coordinates
    use perf_mod, only : t_startf, t_stopf

    type (parallel_t) :: par
    type (element_t) :: elem(:)
    type (cartesian3D_t) :: pinside
    integer :: nslots, ie, num_neighbors, need_conservation, i, j
    logical :: slmm, cisl, qos, sl_test

#ifdef HOMME_ENABLE_COMPOSE
    call t_startf('sl_init1')
    if (transport_alg > 0) then
       call sl_parse_transport_alg(transport_alg, slmm, cisl, qos, sl_test)
       if (par%masterproc .and. nu_q > 0) &
            print *, 'COMPOSE> use HV; nu_q, all:', nu_q, semi_lagrange_hv_q_all
       nslots = nlev*qsize
       sl_mpi = 0
       call slmm_get_mpi_pattern(sl_mpi)
       if (.not. slmm .or. cisl .or. sl_mpi == 0) then
          call abortmp('Only slmm with SL MPI is supported')
       end if
       call interpolate_tracers_init()
       ! Technically a memory leak, but the array persists for the entire
       ! run, so not a big deal for now.
       allocate(dep_points_all(np,np,nlev,size(elem)))
       do ie = 1, size(elem)
          ! Provide a point inside the target element.
          pinside = change_coordinates(elem(ie)%spherep(2,2))
          num_neighbors = elem(ie)%desc%actual_neigh_edges + 1
          call slmm_init_local_mesh(ie, elem(ie)%desc%neigh_corners, num_neighbors, &
               pinside)
          if (sl_test) then
             do j = 1,np
                do i = 1,np
                   pinside = change_coordinates(elem(ie)%spherep(i,j))
                   call slmm_check_ref2sphere(ie, pinside)
                end do
             end do
          end if
       end do
       call slmm_init_finalize()
       if (semi_lagrange_cdr_alg > 1) then
          need_conservation = 1
          call cedr_sl_init(np, nlev, qsize, qsize_d, timelevels, need_conservation)
       end if
    endif
    call t_stopf('sl_init1')
#endif
  end subroutine sl_init1

  subroutine  Prim_Advec_Tracers_remap_ALE( elem , deriv , hvcoord, hybrid , dt , tl , nets , nete )
    use coordinate_systems_mod, only : cartesian3D_t, cartesian2D_t
    use dimensions_mod,         only : max_neigh_edges
    use bndry_mod,              only : ghost_exchangevfull
    use interpolate_mod,        only : interpolate_tracers, minmax_tracers
    use control_mod,            only : qsplit, nu_q, semi_lagrange_hv_q_all, &
         transport_alg, semi_lagrange_cdr_alg, semi_lagrange_cdr_check
    ! For DCMIP16 supercell test case.
    use control_mod,            only : dcmip16_mu_q
    use prim_advection_base,    only : advance_physical_vis

    implicit none
    type (element_t)     , intent(inout) :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (hybrid_t)      , intent(in   ) :: hybrid
    real(kind=real_kind) , intent(in   ) :: dt
    type (TimeLevel_t)   , intent(in   ) :: tl
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete

    type(cartesian3D_t)   :: dep_points  (np,np)
    real(kind=real_kind)  :: minq        (np,np,nlev,qsize,nets:nete)
    real(kind=real_kind)  :: maxq        (np,np,nlev,qsize,nets:nete)

    integer               :: i,j,k,l,n,q,ie,n0_qdp,np1_qdp
    integer               :: num_neighbors, scalar_q_bounds, info
    logical :: slmm, cisl, qos, sl_test

#ifdef HOMME_ENABLE_COMPOSE
    call t_barrierf('Prim_Advec_Tracers_remap_ALE', hybrid%par%comm)
    call t_startf('Prim_Advec_Tracers_remap_ALE')

    call sl_parse_transport_alg(transport_alg, slmm, cisl, qos, sl_test)

    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)

    ! compute displacements for departure grid store in elem%derived%vstar
    call ALE_RKdss (elem, nets, nete, hybrid, deriv, dt, tl)

    if (barrier) call perf_barrier(hybrid)
    call t_startf('SLMM_v2x')
    do ie = nets, nete
       num_neighbors = elem(ie)%desc%actual_neigh_edges + 1
#if (defined COLUMN_OPENMP)
       !$omp parallel do private(k)
#endif
       do k = 1, nlev
          call ALE_departure_from_gll(dep_points_all(:,:,k,ie), &
               elem(ie)%derived%vstar(:,:,:,k), elem(ie), dt, normalize=.true.)
       end do
    end do
    call t_stopf('SLMM_v2x')

    call t_startf('SLMM_csl')
    !todo Here and in the set-pointer loop for CEDR, do just in the first call.
    do ie = nets, nete
       call slmm_csl_set_elem_data(ie, elem(ie)%metdet, &
            elem(ie)%state%Qdp(:,:,:,:,n0_qdp), &
            elem(ie)%derived%dp, elem(ie)%state%Q, &
            elem(ie)%desc%actual_neigh_edges + 1)
    end do
    call slmm_csl(nets, nete, dep_points_all, minq, maxq, info)
    if (info /= 0) then
       call write_velocity_data(elem, nets, nete, hybrid, deriv, dt, tl)
       call abortmp('slmm_csl returned -1; see output above for more information.')
    end if
    if (barrier) call perf_barrier(hybrid)
    call t_stopf('SLMM_csl')

    if (nu_q > 0) then
       if (semi_lagrange_hv_q_all) then
          n = qsize
       else
          n = 1
       end if
       do ie = nets, nete
          do q = 1, n
             do k = 1, nlev
                elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%state%Q(:,:,k,q) * &
                     elem(ie)%state%dp3d(:,:,k,tl%np1)
             enddo
          enddo
       end do
       call advance_hypervis_scalar(elem, hvcoord, hybrid, deriv, tl%np1, np1_qdp, nets, nete, dt, n)
       do ie = nets, nete
          do q = 1, n
             do k = 1, nlev
                elem(ie)%state%Q(:,:,k,q) = elem(ie)%state%Qdp(:,:,k,q,np1_qdp) / &
                     elem(ie)%state%dp3d(:,:,k,tl%np1)
             enddo
          enddo
       end do
    end if

    ! CEDR works with either classical SL or IR.
    if (semi_lagrange_cdr_alg > 1) then
       scalar_q_bounds = 0
       call cedr_sl_set_pointers_begin(nets, nete)
       do ie = nets, nete
          call cedr_sl_set_spheremp(ie, elem(ie)%spheremp)
          call cedr_sl_set_Qdp(ie, elem(ie)%state%Qdp, n0_qdp, np1_qdp)
          call cedr_sl_set_dp3d(ie, elem(ie)%state%dp3d, tl%np1)
          call cedr_sl_set_Q(ie, elem(ie)%state%Q)
       end do
       call cedr_sl_set_pointers_end()
       call t_startf('CEDR')
       call cedr_sl_run(minq, maxq, nets, nete)
       if (barrier) call perf_barrier(hybrid)
       call t_stopf('CEDR')
       call t_startf('CEDR_local')
       call cedr_sl_run_local(minq, maxq, nets, nete, scalar_q_bounds, limiter_option)
       if (barrier) call perf_barrier(hybrid)
       call t_stopf('CEDR_local')
       call t_startf('SL_dss')
       call dss_Qdp(elem, nets, nete, hybrid, np1_qdp)
       if (barrier) call perf_barrier(hybrid)
       call t_stopf('SL_dss')
    else
       do ie = nets, nete
          do k = 1, nlev
             do q = 1, qsize
                elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%state%Q(:,:,k,q) * &
                     elem(ie)%state%dp3d(:,:,k,tl%np1)
             enddo
          enddo
       end do
    end if
    if (semi_lagrange_cdr_check) then
       call t_startf('CEDR_check')
       call cedr_sl_check(minq, maxq, nets, nete)
       if (barrier) call perf_barrier(hybrid)
       call t_stopf('CEDR_check')
    end if
    ! physical viscosity for supercell test case
    if (dcmip16_mu_q > 0) then
       call advance_physical_vis(elem, hvcoord, hybrid, deriv, tl%np1, np1_qdp, nets, nete, dt, dcmip16_mu_q)
    endif
    call t_stopf('Prim_Advec_Tracers_remap_ALE')
#endif
  end subroutine Prim_Advec_Tracers_remap_ALE

  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE ALE_RKDSS
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
       call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%derived%vstar,2*nlev,0,2*nlev)
    enddo

    call t_startf('ALE_RKdss_bexchV')
    call bndry_exchangeV(hy,edge_g)
    call t_stopf('ALE_RKdss_bexchV')

    do ie=nets,nete
       call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%derived%vstar,2*nlev,0,2*nlev)
       do k=1, nlev
          elem(ie)%derived%vstar(:,:,1,k) = elem(ie)%derived%vstar(:,:,1,k)*elem(ie)%rspheremp(:,:)
          elem(ie)%derived%vstar(:,:,2,k) = elem(ie)%derived%vstar(:,:,2,k)*elem(ie)%rspheremp(:,:)
       end do
    end do
  end subroutine ALE_RKdss

  subroutine write_velocity_data(elem, nets, nete, hy, deriv, dt, tl)
    use derivative_mod,  only : derivative_t, ugradv_sphere
    use edgetype_mod,    only : EdgeBuffer_t
    use bndry_mod,       only : bndry_exchangev
    use kinds,           only : real_kind
    use hybrid_mod,      only : hybrid_t
    use element_mod,     only : element_t
    use dimensions_mod,  only : np, nlev
    implicit none

    type (element_t)     , intent(inout) :: elem(:)
    integer              , intent(in)    :: nets
    integer              , intent(in)    :: nete
    type (hybrid_t)      , intent(in)    :: hy
    type (derivative_t)  , intent(in)    :: deriv
    real (kind=real_kind), intent(in)    :: dt
    type (TimeLevel_t)   , intent(in)    :: tl

    integer :: ie, i, j, k, np1
    real (kind=real_kind) :: max_v, max_vstar

    np1 = tl%np1
    max_v = 0.d0
    max_vstar = 0.d0
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                max_v = max(max_v, sqrt(sum(elem(ie)%state%v(i,j,:,k,np1)**2)))
                max_vstar = max(max_vstar, sqrt(sum(elem(ie)%derived%vstar(i,j,:,k)**2)))
             end do
          end do
       end do
    end do
    print *, 'max_v, max_vstar on rank', max_v, max_vstar, hy%par%rank
  end subroutine write_velocity_data

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
  subroutine ALE_departure_from_gll(acart, vstar, elem, dt, normalize)
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
    logical, intent(in) :: normalize

    integer                               :: i,j

    real (kind=real_kind)                 :: uxyz (np,np,3), norm

    ! convert velocity from lat/lon to cartesian 3D
    do i=1,3
       ! Summing along the third dimension is a sum over components for each point.
       ! (This is just a faster way of doing a dot product for each grid point,
       ! since reindexing the inputs to use the intrinsic effectively would be
       ! just asking for trouble.)
       uxyz(:,:,i)=sum( elem%vec_sphere2cart(:,:,i,:)*vstar(:,:,:) ,3)
    end do
    ! compute departure point
    ! crude, 1st order accurate approximation.  to be improved
    do i=1,np
       do j=1,np
          acart(i,j) = change_coordinates(elem%spherep(i,j)) 
          acart(i,j)%x = acart(i,j)%x - dt*uxyz(i,j,1)/rearth
          acart(i,j)%y = acart(i,j)%y - dt*uxyz(i,j,2)/rearth
          acart(i,j)%z = acart(i,j)%z - dt*uxyz(i,j,3)/rearth
          if (normalize) then
             norm = sqrt(acart(i,j)%x*acart(i,j)%x + acart(i,j)%y*acart(i,j)%y + &
                  acart(i,j)%z*acart(i,j)%z)
             acart(i,j)%x = acart(i,j)%x / norm
             acart(i,j)%y = acart(i,j)%y / norm
             acart(i,j)%z = acart(i,j)%z / norm
          end if
       enddo
    enddo

  end subroutine ALE_departure_from_gll

  subroutine dss_Qdp(elem, nets, nete, hybrid, np1_qdp)
    use edgetype_mod,    only : EdgeBuffer_t
    use bndry_mod,       only : bndry_exchangev
    use hybrid_mod,      only : hybrid_t
    use element_mod,     only : element_t
    implicit none

    type (element_t), intent(inout) :: elem(:)
    integer         , intent(in   ) :: nets, nete, np1_qdp
    type (hybrid_t) , intent(in)    :: hybrid
    integer                         :: ie, q, k

    do ie = nets, nete
       do q = 1, qsize
          do k = 1, nlev
             elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%state%Qdp(:,:,k,q,np1_qdp)*elem(ie)%spheremp(:,:)
          end do
       end do
       call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%Qdp(:,:,:,:,np1_qdp), qsize*nlev, 0, qsize*nlev)
    enddo

    call t_startf('SLMM_bexchV')
    call bndry_exchangeV(hybrid, edge_g)
    call t_stopf('SLMM_bexchV')

    do ie = nets, nete
       call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%Qdp(:,:,:,:,np1_qdp), qsize*nlev, 0, qsize*nlev)
       do q = 1, qsize
          do k = 1, nlev
             elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%state%Qdp(:,:,k,q,np1_qdp)*elem(ie)%rspheremp(:,:)
          end do
       end do
    end do
  end subroutine dss_Qdp

  ! Replacement for edge_mod_base::ghostvpack_unoriented, which has a strange
  ! 'threadsafe' module variable that causes a race condition when HORIZ and
  ! COLUMN threading are on at the same time.
  subroutine amb_ghostvpack_unoriented(edge,v,nc,vlyr,kptr,desc)
    use edgetype_mod, only : edgedescriptor_t, ghostbuffer3d_t 
    implicit none
    type (Ghostbuffer3D_t),intent(inout) :: edge
    integer,              intent(in)   :: vlyr
    integer,              intent(in)   :: nc
    real (kind=real_kind),intent(in)   :: v(nc,nc,vlyr)
    integer,              intent(in)   :: kptr
    type (EdgeDescriptor_t),intent(in) :: desc

    integer :: k,l,l_local,is

    do l_local=1,desc%actual_neigh_edges
       l=desc%loc2buf(l_local)
       is = desc%putmapP_ghost(l)
       do k=1,vlyr
          edge%buf(:,:,kptr+k,is) = v(:,:,k)  
       enddo
    end do
  end subroutine amb_ghostvpack_unoriented

  subroutine amb_ghostvunpack_unoriented(edge,v,nc,vlyr,kptr,desc,GlobalId,u)
    use edgetype_mod, only : Ghostbuffer3d_t, EdgeDescriptor_t
    implicit none

    type (Ghostbuffer3D_t),intent(inout)  :: edge
    integer,               intent(in)     :: vlyr
    integer,               intent(in)     :: nc
    real (kind=real_kind), intent(out)    :: v(nc,nc,vlyr,*)
    integer,               intent(in)     :: kptr
    type (EdgeDescriptor_t),intent(in)    :: desc
    integer(kind=int_kind),intent(in)     :: GlobalId
    real (kind=real_kind), intent(in)     :: u(nc,nc,vlyr)

    integer :: k,l,n,is,m,pid,gid

    m=0
    gid = GlobalID
    do n=1,desc%actual_neigh_edges+1
       l = desc%loc2buf(m+1)
       pid = desc%globalID(l)
       if (m==desc%actual_neigh_edges .OR. pid < gid) then
          gid = -1
          v(:,:,:,n) = u(:,:,:)
       else
          m = m+1
          is = desc%getmapP_ghost(l)
          do k=1,vlyr
             v(:,:,k,n) = edge%buf(:,:,kptr+k,is) 
          enddo
       end if
    end do
  end subroutine amb_ghostvunpack_unoriented

  subroutine perf_barrier(hybrid)
    use hybrid_mod, only : hybrid_t
    implicit none
    type (hybrid_t), intent(in) :: hybrid
    integer :: ierr

#ifdef HORIZ_OPENMP
    !$OMP BARRIER
#endif
    if (hybrid%ithr == 0) call mpi_barrier(hybrid%par%comm, ierr)
#ifdef HORIZ_OPENMP
    !$OMP BARRIER
#endif
  end subroutine perf_barrier

  subroutine advance_hypervis_scalar(elem, hvcoord , hybrid , deriv , nt , nt_qdp , nets , nete , dt2, nq)
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
    use bndry_mod      , only : bndry_exchangev
    use perf_mod       , only : t_startf, t_stopf                          ! _EXTERNAL
    use control_mod    , only : nu_q, nu_p, hypervis_subcycle_q
    implicit none
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)     , intent(in   )         :: hvcoord
    type (hybrid_t)      , intent(in   )         :: hybrid
    type (derivative_t)  , intent(in   )         :: deriv
    integer              , intent(in   )         :: nt
    integer              , intent(in   )         :: nt_qdp
    integer              , intent(in   )         :: nets
    integer              , intent(in   )         :: nete
    integer              , intent(in   )         :: nq
    real (kind=real_kind), intent(in   )         :: dt2

    ! local
    real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: Qtens
    real (kind=real_kind), dimension(np,np,nlev                ) :: dp
    real (kind=real_kind) :: dt
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
             do q = 1 , nq
                do k = 1 , nlev
                   dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - dt2*elem(ie)%derived%divdp_proj(:,:,k)
                   Qtens(:,:,k,q,ie) = elem(ie)%derived%dpdiss_ave(:,:,k)*&
                        elem(ie)%state%Qdp(:,:,k,q,nt_qdp) / dp(:,:,k)
                enddo
             enddo

          else
#if (defined COLUMN_OPENMP)
             !$omp parallel do private(q,k) collapse(2)
#endif
             do q = 1 , nq
                do k = 1 , nlev
                   dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - dt2*elem(ie)%derived%divdp_proj(:,:,k)
                   Qtens(:,:,k,q,ie) = hvcoord%dp0(k)*elem(ie)%state%Qdp(:,:,k,q,nt_qdp) / dp(:,:,k)
                enddo
             enddo
          endif
       enddo ! ie loop

       ! compute biharmonic operator. Qtens = input and output
       call biharmonic_wk_scalar( elem , Qtens , deriv , edge_g , hybrid , nets , nete, nq )

       do ie = nets , nete
#if (defined COLUMN_OPENMP)
          !$omp parallel do private(q,k,j,i)
#endif
          do q = 1 , nq
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

          enddo
          call edgeVpack_nlyr(edge_g , elem(ie)%desc, elem(ie)%state%Qdp(:,:,:,:,nt_qdp) , nq*nlev , 0 , nq*nlev )
       enddo ! ie loop

       call t_startf('ah_scalar_bexchV')
       call bndry_exchangeV( hybrid , edge_g )
       call t_stopf('ah_scalar_bexchV')

       do ie = nets , nete
          call edgeVunpack_nlyr(edge_g , elem(ie)%desc, elem(ie)%state%Qdp(:,:,:,:,nt_qdp) , nq*nlev , 0, nq*nlev)
#if (defined COLUMN_OPENMP)
          !$omp parallel do private(q,k) collapse(2)
#endif
          do q = 1 , nq
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

  subroutine biharmonic_wk_scalar(elem,qtens,deriv,edgeq,hybrid,nets,nete,nq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute weak biharmonic operator
    !    input:  qtens = Q
    !    output: qtens = weak biharmonic of Q
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use control_mod, only : hypervis_scaling
    use derivative_mod, only : laplace_sphere_wk

    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    integer :: nets,nete,nq
    real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: qtens
    type (EdgeBuffer_t)  , intent(inout) :: edgeq
    type (derivative_t)  , intent(in) :: deriv

    ! local
    integer :: k,kptr,i,j,ie,ic,q
    real (kind=real_kind), dimension(np,np) :: lap_p
    logical var_coef1

    !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
    !so tensor is only used on second call to laplace_sphere_wk
    var_coef1 = .true.
    if(hypervis_scaling > 0)    var_coef1 = .false.

    do ie=nets,nete
#if (defined COLUMN_OPENMP)
       !$omp parallel do private(k, q, lap_p)
#endif
       do q=1,nq
          do k=1,nlev    !  Potential loop inversion (AAM)
             lap_p(:,:)=qtens(:,:,k,q,ie)
             ! Original use of qtens on left and right hand sides caused OpenMP errors (AAM)
             qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=var_coef1)
          enddo
          call edgeVpack_nlyr(edgeq, elem(ie)%desc, qtens(:,:,:,q,ie),nlev,nlev*(q-1),nq*nlev)
       enddo
    enddo

    call t_startf('biwksc_bexchV')
    call bndry_exchangeV(hybrid,edgeq)
    call t_stopf('biwksc_bexchV')

    do ie=nets,nete

       ! apply inverse mass matrix, then apply laplace again
#if (defined COLUMN_OPENMP)
       !$omp parallel do private(k, q, lap_p)
#endif
       do q=1,nq      
          call edgeVunpack_nlyr(edgeq,elem(ie)%desc,qtens(:,:,:,q,ie),nlev,nlev*(q-1),nq*nlev)
          do k=1,nlev    !  Potential loop inversion (AAM)
             lap_p(:,:)=elem(ie)%rspheremp(:,:)*qtens(:,:,k,q,ie)
             qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=.true.)
          enddo
       enddo
    enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine biharmonic_wk_scalar

end module sl_advection
