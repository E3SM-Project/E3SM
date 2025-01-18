#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module compose_mod

  implicit none

  logical :: compose_h2d = .false., compose_d2h = .false.
  logical, private :: control_kokkos_init_and_fin = .true.
  
  interface
     
#ifdef HOMME_ENABLE_COMPOSE
     subroutine kokkos_init() bind(c)
     end subroutine kokkos_init

     subroutine kokkos_finalize() bind(c)
     end subroutine kokkos_finalize

     subroutine cedr_unittest(comm, nerr) bind(c)
       use iso_c_binding, only: c_int
       integer(kind=c_int), value, intent(in) :: comm
       integer(kind=c_int), intent(out) :: nerr
     end subroutine cedr_unittest

     subroutine cedr_init_impl(comm, cdr_alg, use_sgi, gid_data, rank_data, &
          ncell, nlclcell, nlev, np, qsize, independent_time_steps, hard_zero, &
          gid_data_sz, rank_data_sz) bind(c)
       use iso_c_binding, only: c_int, c_bool
       integer(kind=c_int), value, intent(in) :: comm, cdr_alg, ncell, nlclcell, nlev, np, &
            qsize, gid_data_sz, rank_data_sz
       logical(kind=c_bool), value, intent(in) :: use_sgi, independent_time_steps, hard_zero
       integer(kind=c_int), intent(in) :: gid_data(gid_data_sz), rank_data(rank_data_sz)
     end subroutine cedr_init_impl

     subroutine slmm_init_impl(comm, transport_alg, np, nlev, qsize, qsize_d, &
          nelem, nelemd, cubed_sphere_map, geometry, lid2gid, lid2facenum, &
          nbr_id_rank, nirptr, sl_halo, sl_traj_3d, sl_traj_nsubstep, sl_nearest_point_lev, &
          lid2gid_sz, lid2facenum_sz, nbr_id_rank_sz, nirptr_sz) bind(c)
       use iso_c_binding, only: c_int
       integer(kind=c_int), value, intent(in) :: comm, transport_alg, np, nlev, qsize, &
            qsize_d, nelem, nelemd, cubed_sphere_map, geometry, sl_halo, sl_traj_3d, &
            sl_traj_nsubstep, sl_nearest_point_lev, lid2gid_sz, lid2facenum_sz, &
            nbr_id_rank_sz, nirptr_sz
       integer(kind=c_int), intent(in) :: lid2gid(lid2gid_sz), lid2facenum(lid2facenum_sz), &
            nbr_id_rank(nbr_id_rank_sz), nirptr(nirptr_sz)
     end subroutine slmm_init_impl

     subroutine slmm_init_plane(Sx, Sy, Lx, Ly) bind(c)
       use iso_c_binding, only: c_double
       real(kind=c_double), value, intent(in) :: Sx, Sy, Lx, Ly
     end subroutine slmm_init_plane

     subroutine cedr_query_bufsz(sendsz, recvsz) bind(c)
       use iso_c_binding, only: c_int
       integer(kind=c_int), intent(out) :: sendsz, recvsz
     end subroutine cedr_query_bufsz

     subroutine cedr_set_bufs(sendbuf, recvbuf, sendsz, recvsz) bind(c)
       use iso_c_binding, only: c_double, c_int
       integer(kind=c_int), value, intent(in) :: sendsz, recvsz
       real(kind=c_double), intent(in) :: sendbuf(sendsz), recvbuf(recvsz)
     end subroutine cedr_set_bufs

     subroutine cedr_set_null_bufs() bind(c)
     end subroutine cedr_set_null_bufs

     subroutine cedr_set_ie2gci(ie, gci) bind(c)
       use iso_c_binding, only: c_int
       integer(kind=c_int), value, intent(in) :: ie, gci
     end subroutine cedr_set_ie2gci

     subroutine cedr_sl_init(np, nlev, qsize, qsize_d, timelevels, &
          need_conservation) bind(c)
       use iso_c_binding, only: c_int
       integer(kind=c_int), value, intent(in) :: np, nlev, qsize, qsize_d, timelevels, &
            need_conservation
     end subroutine cedr_sl_init

     subroutine cedr_sl_set_pointers_begin(nets, nete) bind(c)
       use iso_c_binding, only: c_int
       integer(kind=c_int), value, intent(in) :: nets, nete
     end subroutine cedr_sl_set_pointers_begin

     subroutine cedr_sl_set_spheremp(ie, spheremp) bind(c)
       use iso_c_binding, only: c_int, c_double
       use dimensions_mod, only : np
       integer(kind=c_int), value, intent(in) :: ie
       real(kind=c_double), intent(in) :: spheremp(np,np)
     end subroutine cedr_sl_set_spheremp

     subroutine cedr_sl_set_dp0(dp0) bind(c)
       use iso_c_binding, only: c_double
       use dimensions_mod, only : nlev
       real(kind=c_double), intent(in) :: dp0(nlev)
     end subroutine cedr_sl_set_dp0

     subroutine cedr_sl_set_Qdp(ie, Qdp, n0_qdp, np1_qdp) bind(c)
       use iso_c_binding, only: c_int, c_double
       use dimensions_mod, only : nlev, np, qsize_d
       integer(kind=c_int), value, intent(in) :: ie, n0_qdp, np1_qdp
       real(kind=c_double), intent(in) :: Qdp(np,np,nlev,qsize_d,2)
     end subroutine cedr_sl_set_Qdp

     subroutine cedr_sl_set_dp3d(ie, dp3d, tl_np1) bind(c)
       use iso_c_binding, only: c_int, c_double
       use dimensions_mod, only : nlev, np
       use element_state,  only : timelevels
       integer(kind=c_int), value, intent(in) :: ie, tl_np1
       real(kind=c_double), intent(in) :: dp3d(np,np,nlev,timelevels)
     end subroutine cedr_sl_set_dp3d

     subroutine cedr_sl_set_dp(ie, dp) bind(c)
       use iso_c_binding , only : c_int, c_double
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np
       use element_state,  only : timelevels
       integer(kind=c_int), value, intent(in) :: ie
       real(kind=c_double), intent(in) :: dp(np,np,nlev)
     end subroutine cedr_sl_set_dp

     subroutine cedr_sl_set_Q(ie, Q) bind(c)
       use iso_c_binding, only: c_int, c_double
       use dimensions_mod, only : nlev, np, qsize_d
       integer(kind=c_int), value, intent(in) :: ie
       real(kind=c_double), intent(in) :: Q(np,np,nlev,qsize_d)
     end subroutine cedr_sl_set_Q

     subroutine cedr_sl_set_pointers_end(h2d, d2h) bind(c)
       use iso_c_binding, only: c_bool
       logical(kind=c_bool), value, intent(in) :: h2d, d2h
     end subroutine cedr_sl_set_pointers_end

     subroutine cedr_sl_run_global(minq, maxq, nets, nete) bind(c)
       use iso_c_binding, only: c_int, c_double
       use dimensions_mod, only : nlev, np, qsize
       integer(kind=c_int), value, intent(in) :: nets, nete
       real(kind=c_double), intent(in) :: minq(np,np,nlev,qsize,nets:nete)
       real(kind=c_double), intent(in) :: maxq(np,np,nlev,qsize,nets:nete)
     end subroutine cedr_sl_run_global

     subroutine cedr_sl_run_local(minq, maxq, nets, nete, use_ir, limiter_option) bind(c)
       use iso_c_binding, only: c_int, c_double
       use dimensions_mod, only : nlev, np, qsize
       integer(kind=c_int), value, intent(in) :: nets, nete, use_ir, limiter_option
       real(kind=c_double), intent(in) :: minq(np,np,nlev,qsize,nets:nete)
       real(kind=c_double), intent(in) :: maxq(np,np,nlev,qsize,nets:nete)
     end subroutine cedr_sl_run_local

     subroutine cedr_sl_check(minq, maxq, nets, nete) bind(c)
       use iso_c_binding, only: c_int, c_double
       use dimensions_mod, only : nlev, np, qsize
       integer(kind=c_int), value, intent(in) :: nets, nete
       real(kind=c_double), intent(in) :: minq(np,np,nlev,qsize,nets:nete)
       real(kind=c_double), intent(in) :: maxq(np,np,nlev,qsize,nets:nete)
     end subroutine cedr_sl_check

     subroutine slmm_init_local_mesh(ie, neigh_corners, num_neighbors, pinside, &
          patch_size) bind(c)
       use iso_c_binding, only: c_int
       use coordinate_systems_mod, only : cartesian3D_t
       integer(kind=c_int), value, intent(in) :: ie, num_neighbors, patch_size
       type(cartesian3D_t), intent(in) :: neigh_corners(4,patch_size), pinside
     end subroutine slmm_init_local_mesh

     subroutine slmm_query_bufsz(sendsz, recvsz) bind(c)
       use iso_c_binding, only: c_int
       integer(kind=c_int), intent(out) :: sendsz, recvsz
     end subroutine slmm_query_bufsz

     subroutine slmm_set_bufs(sendbuf, recvbuf, sendsz, recvsz) bind(c)
       use iso_c_binding, only: c_double, c_int
       integer(kind=c_int), intent(in) :: sendsz, recvsz
       real(kind=c_double), intent(in) :: sendbuf(sendsz), recvbuf(recvsz)
     end subroutine slmm_set_bufs

     subroutine slmm_set_null_bufs() bind(c)
     end subroutine slmm_set_null_bufs

     subroutine slmm_init_finalize() bind(c)
     end subroutine slmm_init_finalize

     subroutine slmm_check_ref2sphere(ie, sphere_cart_coord) bind(c)
       use iso_c_binding, only: c_int
       use coordinate_systems_mod, only : cartesian3D_t
       integer(kind=c_int), value, intent(in) :: ie
       type(cartesian3D_t), intent(in) :: sphere_cart_coord
     end subroutine slmm_check_ref2sphere

     subroutine slmm_set_hvcoord(etai_beg, etai_end, etam) bind(c)
       use iso_c_binding, only: c_double
       use dimensions_mod, only : nlev
       real(kind=c_double), value, intent(in) :: etai_beg, etai_end
       real(kind=c_double), intent(in) :: etam(nlev)
     end subroutine slmm_set_hvcoord

     subroutine slmm_calc_v_departure(nets, nete, step, dtsub, dep_points, &
          dep_points_ndim, vnode, vdep, info) bind(c)
       use iso_c_binding, only: c_int, c_double
       use dimensions_mod, only : np, nlev, nelemd, qsize
       use coordinate_systems_mod, only : cartesian3D_t
       integer(kind=c_int), value, intent(in) :: nets, nete, step, dep_points_ndim
       real(kind=c_double), value, intent(in) :: dtsub
       real(kind=c_double), intent(inout) :: dep_points(dep_points_ndim,np,np,nlev,nelemd)
       real(kind=c_double), intent(in) :: vnode(dep_points_ndim,np,np,nlev,nelemd)
       real(kind=c_double), intent(out) :: vdep(dep_points_ndim,np,np,nlev,nelemd)
       integer(kind=c_int), intent(out) :: info
     end subroutine slmm_calc_v_departure

     subroutine slmm_csl_set_elem_data(ie, metdet, qdp, n0_qdp, dp, q, nelem_in_patch, &
          h2d, d2h) bind(c)
       use iso_c_binding, only: c_int, c_double, c_bool
       use dimensions_mod, only : nlev, np, qsize
       real(kind=c_double), intent(in) :: metdet(np,np), qdp(np,np,nlev,qsize,2), &
            dp(np,np,nlev), q(np,np,nlev,qsize)
       integer(kind=c_int), value, intent(in) :: ie, n0_qdp, nelem_in_patch
       logical(kind=c_bool), value, intent(in) :: h2d, d2h
     end subroutine slmm_csl_set_elem_data

     subroutine slmm_csl(nets, nete, dep_points, dep_points_ndim, minq, maxq, info) bind(c)
       use iso_c_binding, only: c_int, c_double
       use dimensions_mod, only : np, nlev, nelemd, qsize
       use coordinate_systems_mod, only : cartesian3D_t
       integer(kind=c_int), value, intent(in) :: nets, nete, dep_points_ndim
       ! dep_points is const in principle, but if lev <=
       ! semi_lagrange_nearest_point_lev, a departure point may be altered if
       ! the winds take it outside of the comm halo.
       real(kind=c_double), intent(inout) :: dep_points(dep_points_ndim,np,np,nlev,nelemd)
       real(kind=c_double), intent(in) :: &
            minq(np,np,nlev,qsize,nelemd), maxq(np,np,nlev,qsize,nelemd)
       integer(kind=c_int), intent(out) :: info
     end subroutine slmm_csl

     ! Temporary for work on and characterization of the SL MPI pattern.
     subroutine slmm_get_mpi_pattern(sl_mpi) bind(c)
       use iso_c_binding, only: c_int
       integer(kind=c_int), intent(out) :: sl_mpi
     end subroutine slmm_get_mpi_pattern
#endif

     subroutine cedr_finalize() bind(c)
     end subroutine cedr_finalize

     subroutine slmm_finalize() bind(c)
     end subroutine slmm_finalize

  end interface

contains

  subroutine compose_control_kokkos_init_and_fin(control)
    logical, intent(in) :: control
    control_kokkos_init_and_fin = control
  end subroutine compose_control_kokkos_init_and_fin

  subroutine compose_init(par, elem, GridVertex, init_kokkos)
    use iso_c_binding, only: c_bool
    use parallel_mod, only: parallel_t, abortmp
    use dimensions_mod, only: np, nlev, qsize, qsize_d, nelem, nelemd, ne_x, ne_y
    use element_mod, only: element_t
    use gridgraph_mod, only: GridVertex_t
    use control_mod, only: semi_lagrange_cdr_alg, transport_alg, cubed_sphere_map, &
         semi_lagrange_halo, semi_lagrange_trajectory_nsubstep, &
         semi_lagrange_nearest_point_lev, dt_remap_factor, dt_tracer_factor, geometry
    use physical_constants, only: Sx, Sy, Lx, Ly
    use scalable_grid_init_mod, only: sgi_is_initialized, sgi_get_rank2sfc, &
         sgi_gid2igv
    use perf_mod, only: t_startf, t_stopf

    type (parallel_t), intent(in) :: par
    type (element_t), intent(in) :: elem(:)
    type (GridVertex_t), intent(in), target :: GridVertex(:)
    logical, optional, intent(in) :: init_kokkos

    integer, allocatable :: &
         nbr_id_rank(:), nirptr(:), & ! (GID, rank) in local mesh patch, starting with own
         ! These are for non-scalable grid initialization, still used for RRM.
         sc2gci(:), sc2rank(:)        ! space curve index -> (GID, rank)
    integer :: lid2gid(nelemd), lid2facenum(nelemd)
    integer :: i, j, k, sfc, gid, igv, sc, geometry_type, sl_traj_3d
    ! To map SFC index to IDs and ranks
    logical(kind=c_bool) :: use_sgi, owned, independent_time_steps, hard_zero
    integer, allocatable :: owned_ids(:)
    integer, pointer :: rank2sfc(:) => null()
    integer, target :: null_target(1)
    logical :: call_init_kokkos

#ifdef HOMME_ENABLE_COMPOSE
    call t_startf('compose_init')

    call_init_kokkos = .true.
    if (present(init_kokkos)) call_init_kokkos = init_kokkos
    if (control_kokkos_init_and_fin .and. call_init_kokkos) call kokkos_init()

    use_sgi = sgi_is_initialized()
    hard_zero = .true.

    independent_time_steps = dt_remap_factor < dt_tracer_factor
    
    if (semi_lagrange_halo < 1) then
       ! For test problems, the relationship between dt_tracer_factor and halo
       ! may not be clear. But for real problems, the advective CFL implies that
       ! a parcel can cross a cell in three time steps. Since this is closely
       ! related to the dynamics' tstep, dt_tracer_factor is meaningful,
       ! implying:
       semi_lagrange_halo = (dt_tracer_factor + 2) / 3
       if (semi_lagrange_halo < 1) semi_lagrange_halo = 1
    end if

    geometry_type = 0 ! sphere
    if (trim(geometry) == "plane") then
       geometry_type = 1
       if (min(ne_x, ne_y) < 5) then
          ! If we really want to, we can support ne := min(ne_x, ne_y) >= 3 just
          ! by setting halo = 1 if ne < 5. But for now that's not important,
          ! support ne no less than 5.
          call abortmp('SL transport for planar geometry does not support min(ne_x, ne_y) < 5.')
       end if
    end if

    if ( semi_lagrange_cdr_alg == 2 .or. semi_lagrange_cdr_alg == 20 .or. &
         semi_lagrange_cdr_alg == 21) then
       if (use_sgi) then
          call sgi_get_rank2sfc(rank2sfc)
          allocate(owned_ids(size(GridVertex)))
          j = rank2sfc(par%rank+1)
          do i = 1, size(GridVertex)
             sfc = GridVertex(i)%SpaceCurve
             owned = sfc >= j .and. sfc < rank2sfc(par%rank+2)
             if (owned) then
                owned_ids(sfc - j + 1) = GridVertex(i)%number - 1
             end if
          end do
       else
          allocate(sc2gci(nelem), sc2rank(nelem))
          do i = 1, nelem
             sc = GridVertex(i)%SpaceCurve + 1
             sc2gci(sc) = i - 1
             sc2rank(sc) = GridVertex(i)%processor_number - 1
          end do
       end if
    else
       ! These lines fix ifort -check catches. The data are not used,
       ! but the function call cedr_init_impl makes -check think they
       ! are used.
       allocate(owned_ids(1), sc2gci(1), sc2rank(1))
       rank2sfc => null_target
    end if
    if (use_sgi) then
       if (.not. allocated(owned_ids)) allocate(owned_ids(1))
       call cedr_init_impl(par%comm, semi_lagrange_cdr_alg, &
            use_sgi, owned_ids, rank2sfc, nelem, nelemd, nlev, np, qsize, &
            independent_time_steps, hard_zero, size(owned_ids), size(rank2sfc))
    else
       if (.not. allocated(sc2gci)) allocate(sc2gci(1), sc2rank(1))
       call cedr_init_impl(par%comm, semi_lagrange_cdr_alg, &
            use_sgi, sc2gci, sc2rank, nelem, nelemd, nlev, np, qsize, &
            independent_time_steps, hard_zero, size(sc2gci), size(sc2rank))
    end if
    if (allocated(sc2gci)) deallocate(sc2gci, sc2rank)
    if (allocated(owned_ids)) deallocate(owned_ids)

    if (transport_alg > 1) then
       k = 0
       do i = 1, nelemd
          k = k + 2*(elem(i)%desc%actual_neigh_edges + 2)
       end do
       allocate(nbr_id_rank(k), nirptr(nelemd+1))
       k = 1
       do i = 1, nelemd
          nirptr(i) = k - 1
          gid = elem(i)%globalID
          lid2gid(i) = gid
          lid2facenum(i) = elem(i)%faceNum
          nbr_id_rank(k) = gid
          if (use_sgi) then
             igv = sgi_gid2igv(gid)
          else
             igv = gid
          end if
          nbr_id_rank(k+1) = GridVertex(igv)%processor_number - 1
          k = k + 2
          do j = 1, size(elem(i)%desc%globalID_neigh_corners)
             gid = elem(i)%desc%globalID_neigh_corners(j)
             if (use_sgi) then
                igv = sgi_gid2igv(gid)
             else
                igv = gid
             end if
             nbr_id_rank(k) = gid
             nbr_id_rank(k+1) = GridVertex(igv)%processor_number - 1
             k = k + 2
          end do
       end do
       nirptr(nelemd+1) = k - 1
       sl_traj_3d = 0
       if (independent_time_steps) sl_traj_3d = 1
       call slmm_init_impl(par%comm, transport_alg, np, nlev, qsize, qsize_d, &
            nelem, nelemd, cubed_sphere_map, geometry_type, lid2gid, lid2facenum, &
            nbr_id_rank, nirptr, semi_lagrange_halo, sl_traj_3d, &
            semi_lagrange_trajectory_nsubstep, semi_lagrange_nearest_point_lev, &
            size(lid2gid), size(lid2facenum), size(nbr_id_rank), size(nirptr))
       if (geometry_type == 1) call slmm_init_plane(Sx, Sy, Lx, Ly)
       deallocate(nbr_id_rank, nirptr)
    end if
    call t_stopf('compose_init')
#endif
  end subroutine compose_init

  subroutine compose_finalize(finalize_kokkos)
    logical, optional, intent(in) :: finalize_kokkos

#ifdef HOMME_ENABLE_COMPOSE
    logical :: call_finalize_kokkos

    call cedr_finalize()
    call slmm_finalize()

    call_finalize_kokkos = .true.
    if (present(finalize_kokkos)) call_finalize_kokkos = finalize_kokkos
    if (control_kokkos_init_and_fin .and. call_finalize_kokkos) call kokkos_finalize()
#endif
  end subroutine compose_finalize
  
  subroutine compose_query_bufsz(sendsz, recvsz)
    integer, intent(out) :: sendsz, recvsz

#ifdef HOMME_ENABLE_COMPOSE
    integer :: ssz, rsz

    call slmm_query_bufsz(sendsz, recvsz)
    call cedr_query_bufsz(ssz, rsz)
    sendsz = max(sendsz, ssz)
    recvsz = max(recvsz, rsz)
#endif
  end subroutine compose_query_bufsz

  subroutine compose_set_bufs(sendbuf, recvbuf)
    use kinds, only: real_kind

    real(kind=real_kind), intent(in) :: sendbuf(:), recvbuf(:)

#ifdef HOMME_ENABLE_COMPOSE
    ! CEDR and SLMM can use the same buffers because they operate in sequence
    ! and never leave persistent state in these buffers between top-level calls.
    call slmm_set_bufs(sendbuf, recvbuf, size(sendbuf), size(recvbuf))
    call cedr_set_bufs(sendbuf, recvbuf, size(sendbuf), size(recvbuf))
#endif
  end subroutine compose_set_bufs

  subroutine compose_set_null_bufs()
    use kinds, only: real_kind

#ifdef HOMME_ENABLE_COMPOSE
    call slmm_set_null_bufs()
    call cedr_set_null_bufs()
#endif
  end subroutine compose_set_null_bufs

  subroutine compose_repro_sum(send, recv, nlocal, nfld, comm) bind(c)
    use iso_c_binding, only: c_int, c_double
#ifdef CAM
    use shr_reprosum_mod, only: repro_sum => shr_reprosum_calc
#else
    use repro_sum_mod, only: repro_sum
#endif

    integer(kind=c_int), value, intent(in) :: nlocal, nfld, comm
    real(kind=c_double), intent(in) :: send(nlocal,nfld)
    real(kind=c_double), intent(out) :: recv(nfld)

    call repro_sum(send, recv, nlocal, nlocal, nfld, commid=comm)
  end subroutine compose_repro_sum

end module compose_mod
