module compose_mod

  implicit none

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
          ncell, nlclcell, nlev, independent_time_steps, hard_zero, &
          gid_data_sz, rank_data_sz) bind(c)
       use iso_c_binding, only: c_int, c_bool
       integer(kind=c_int), value, intent(in) :: comm, cdr_alg, ncell, nlclcell, nlev, &
            gid_data_sz, rank_data_sz
       logical(kind=c_bool), value, intent(in) :: use_sgi, independent_time_steps, hard_zero
       integer(kind=c_int), intent(in) :: gid_data(gid_data_sz), rank_data(rank_data_sz)
     end subroutine cedr_init_impl

     subroutine slmm_init_impl(comm, transport_alg, np, nlev, qsize, qsize_d, &
          nelem, nelemd, cubed_sphere_map, lid2gid, lid2facenum, nbr_id_rank, nirptr, &
          sl_nearest_point_lev, lid2gid_sz, lid2facenum_sz, nbr_id_rank_sz, nirptr_sz) bind(c)
       use iso_c_binding, only: c_int
       integer(kind=c_int), value, intent(in) :: comm, transport_alg, np, nlev, qsize, qsize_d, &
            nelem, nelemd, cubed_sphere_map, sl_nearest_point_lev, lid2gid_sz, lid2facenum_sz, &
            nbr_id_rank_sz, nirptr_sz
       integer(kind=c_int), intent(in) :: lid2gid(lid2gid_sz), lid2facenum(lid2facenum_sz), &
            nbr_id_rank(nbr_id_rank_sz), nirptr(nirptr_sz)
     end subroutine slmm_init_impl

     subroutine cedr_query_bufsz(sendsz, recvsz) bind(c)
       use iso_c_binding, only: c_int
       integer(kind=c_int), intent(out) :: sendsz, recvsz
     end subroutine cedr_query_bufsz

     subroutine cedr_set_bufs(sendbuf, recvbuf, sendsz, recvsz) bind(c)
       use iso_c_binding, only: c_double, c_int
       integer(kind=c_int), value, intent(in) :: sendsz, recvsz
       real(kind=c_double), intent(in) :: sendbuf(sendsz), recvbuf(recvsz)
     end subroutine cedr_set_bufs

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
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np
       use element_state,  only : timelevels
       integer, value, intent(in) :: ie
       real(kind=real_kind), intent(in) :: dp(np,np,nlev)
     end subroutine cedr_sl_set_dp

     subroutine cedr_sl_set_Q(ie, Q) bind(c)
       use iso_c_binding, only: c_int, c_double
       use dimensions_mod, only : nlev, np, qsize_d
       integer(kind=c_int), value, intent(in) :: ie
       real(kind=c_double), intent(in) :: Q(np,np,nlev,qsize_d)
     end subroutine cedr_sl_set_Q

     subroutine cedr_sl_set_pointers_end() bind(c)
     end subroutine cedr_sl_set_pointers_end

     subroutine cedr_sl_run(minq, maxq, nets, nete) bind(c)
       use iso_c_binding, only: c_int, c_double
       use dimensions_mod, only : nlev, np, qsize
       real(kind=c_double), intent(in) :: minq(np,np,nlev,qsize,nets:nete)
       real(kind=c_double), intent(in) :: maxq(np,np,nlev,qsize,nets:nete)
       integer(kind=c_int), value, intent(in) :: nets, nete
     end subroutine cedr_sl_run

     subroutine cedr_sl_run_local(minq, maxq, nets, nete, use_ir, limiter_option) bind(c)
       use iso_c_binding, only: c_int, c_double
       use dimensions_mod, only : nlev, np, qsize
       real(kind=c_double), intent(in) :: minq(np,np,nlev,qsize,nets:nete)
       real(kind=c_double), intent(in) :: maxq(np,np,nlev,qsize,nets:nete)
       integer(kind=c_int), value, intent(in) :: nets, nete, use_ir, limiter_option
     end subroutine cedr_sl_run_local

     subroutine cedr_sl_check(minq, maxq, nets, nete) bind(c)
       use iso_c_binding, only: c_int, c_double
       use dimensions_mod, only : nlev, np, qsize
       real(kind=c_double), intent(in) :: minq(np,np,nlev,qsize,nets:nete)
       real(kind=c_double), intent(in) :: maxq(np,np,nlev,qsize,nets:nete)
       integer(kind=c_int), value, intent(in) :: nets, nete
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

     subroutine slmm_init_finalize() bind(c)
     end subroutine slmm_init_finalize

     subroutine slmm_check_ref2sphere(ie, sphere_cart_coord) bind(c)
       use iso_c_binding, only: c_int
       use coordinate_systems_mod, only : cartesian3D_t
       integer(kind=c_int), value, intent(in) :: ie
       type(cartesian3D_t), intent(in) :: sphere_cart_coord
     end subroutine slmm_check_ref2sphere

     subroutine slmm_csl_set_elem_data(ie, metdet, qdp, dp, q, nelem_in_patch) bind(c)
       use iso_c_binding, only: c_int, c_double
       use dimensions_mod, only : nlev, np, qsize
       real(kind=c_double), intent(in) :: metdet(np,np), qdp(np,np,nlev,qsize), &
            dp(np,np,nlev), q(np,np,nlev,qsize)
       integer(kind=c_int), value, intent(in) :: ie, nelem_in_patch
     end subroutine slmm_csl_set_elem_data

     subroutine slmm_csl(nets, nete, dep_points, minq, maxq, info) bind(c)
       use iso_c_binding, only: c_int, c_double
       use dimensions_mod, only : np, nlev, nelemd, qsize
       use coordinate_systems_mod, only : cartesian3D_t
       integer(kind=c_int), value, intent(in) :: nets, nete
       ! dep_points is const in principle, but if lev <=
       ! semi_lagrange_nearest_point_lev, a departure point may be altered if
       ! the winds take it outside of the comm halo.
       type(cartesian3D_t), intent(inout) :: dep_points(np,np,nlev,nelemd)
       real(kind=c_double), intent(in) :: &
            minq(np,np,nlev,qsize,nets:nete), maxq(np,np,nlev,qsize,nets:nete)
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

  subroutine compose_init(par, elem, GridVertex)
    use iso_c_binding, only: c_bool
    use parallel_mod, only: parallel_t, abortmp
    use dimensions_mod, only: np, nlev, qsize, qsize_d, nelem, nelemd
    use element_mod, only: element_t
    use gridgraph_mod, only: GridVertex_t
    use control_mod, only: semi_lagrange_cdr_alg, transport_alg, cubed_sphere_map, &
         semi_lagrange_nearest_point_lev, dt_remap_factor, dt_tracer_factor
    use scalable_grid_init_mod, only: sgi_is_initialized, sgi_get_rank2sfc, &
         sgi_gid2igv
    use perf_mod, only: t_startf, t_stopf

    type (parallel_t), intent(in) :: par
    type (element_t), intent(in) :: elem(:)
    type (GridVertex_t), intent(in), target :: GridVertex(:)
    integer, allocatable :: &
         nbr_id_rank(:), nirptr(:), & ! (GID, rank) in local mesh patch, starting with own
         ! These are for non-scalable grid initialization, still used for RRM.
         sc2gci(:), sc2rank(:)        ! space curve index -> (GID, rank)
    integer :: lid2gid(nelemd), lid2facenum(nelemd)
    integer :: i, j, k, sfc, gid, igv, sc
    ! To map SFC index to IDs and ranks
    logical(kind=c_bool) :: use_sgi, owned, independent_time_steps, hard_zero
    integer, allocatable :: owned_ids(:)
    integer, pointer :: rank2sfc(:) => null()
    integer, target :: null_target(1)

#ifdef HOMME_ENABLE_COMPOSE
    call t_startf('compose_init')
    call kokkos_init()

    use_sgi = sgi_is_initialized()
    hard_zero = .true.

    independent_time_steps = dt_remap_factor < dt_tracer_factor

    if (semi_lagrange_cdr_alg == 2 .or. semi_lagrange_cdr_alg == 20 .or. &
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
       allocate(owned_ids(1))
       rank2sfc => null_target
    end if
    if (use_sgi) then
       call cedr_init_impl(par%comm, semi_lagrange_cdr_alg, &
            use_sgi, owned_ids, rank2sfc, nelem, nelemd, nlev, &
            independent_time_steps, hard_zero, size(owned_ids), size(rank2sfc))
    else
       call cedr_init_impl(par%comm, semi_lagrange_cdr_alg, &
            use_sgi, sc2gci, sc2rank, nelem, nelemd, nlev, &
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
       call slmm_init_impl(par%comm, transport_alg, np, nlev, qsize, qsize_d, &
            nelem, nelemd, cubed_sphere_map, lid2gid, lid2facenum, &
            nbr_id_rank, nirptr, semi_lagrange_nearest_point_lev, &
            size(lid2gid), size(lid2facenum), size(nbr_id_rank), size(nirptr))
       deallocate(nbr_id_rank, nirptr)
    end if
    call t_stopf('compose_init')
#endif
  end subroutine compose_init

  subroutine compose_finalize()
#ifdef HOMME_ENABLE_COMPOSE
    call cedr_finalize()
    call slmm_finalize()
    call kokkos_finalize()
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
#else
    sendsz = 0
    recvsz = 0
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

  subroutine compose_repro_sum(send, recv, nlocal, nfld, comm) bind(c)
    use iso_c_binding, only: c_int, c_double
#ifdef CAM
    use shr_reprosum_mod, only: repro_sum => shr_reprosum_calc
#else
    use repro_sum_mod, only: repro_sum
#endif

    real(kind=c_double), intent(in) :: send(nlocal,nfld)
    real(kind=c_double), intent(out) :: recv(nfld)
    integer(kind=c_int), value, intent(in) :: nlocal, nfld, comm

    call repro_sum(send, recv, nlocal, nlocal, nfld, commid=comm)
  end subroutine compose_repro_sum

end module compose_mod
