module compose_mod

  implicit none

  interface
     subroutine kokkos_init() bind(c)
     end subroutine kokkos_init

     subroutine kokkos_finalize() bind(c)
     end subroutine kokkos_finalize

     subroutine cedr_unittest(comm, nerr) bind(c)
       integer, value, intent(in) :: comm
       integer, intent(out) :: nerr
     end subroutine cedr_unittest

     subroutine cedr_init_impl(comm, cdr_alg, sc2gci, sc2rank, &
          ncell, nlclcell, nlev) bind(c)
       integer, value, intent(in) :: comm, cdr_alg, ncell, nlclcell, nlev
       integer, intent(in) :: sc2gci(:), sc2rank(:)
     end subroutine cedr_init_impl

     subroutine slmm_init_impl(comm, transport_alg, np, nlev, qsize, qsize_d, &
          nelem, nelemd, cubed_sphere_map, lid2gid, lid2facenum, nbr_id_rank, nirptr, &
          sl_nearest_point_lev) bind(c)
       integer, value, intent(in) :: comm, transport_alg, np, nlev, qsize, qsize_d, &
            nelem, nelemd, cubed_sphere_map, sl_nearest_point_lev
       integer, intent(in) :: lid2gid(:), lid2facenum(:), nbr_id_rank(:), nirptr(:)
     end subroutine slmm_init_impl

     subroutine cedr_set_ie2gci(ie, gci) bind(c)
       integer, value, intent(in) :: ie, gci
     end subroutine cedr_set_ie2gci

     subroutine cedr_sl_init(np, nlev, qsize, qsize_d, timelevels, &
          need_conservation) bind(c)
       integer, value, intent(in) :: np, nlev, qsize, qsize_d, timelevels, &
            need_conservation
     end subroutine cedr_sl_init

     subroutine cedr_sl_set_pointers_begin(nets, nete) bind(c)
       integer, value, intent(in) :: nets, nete
     end subroutine cedr_sl_set_pointers_begin

     subroutine cedr_sl_set_spheremp(ie, spheremp) bind(c)
       use kinds         , only : real_kind
       use dimensions_mod, only : np
       integer, value, intent(in) :: ie
       real(kind=real_kind), intent(in) :: spheremp(np,np)
     end subroutine cedr_sl_set_spheremp

     subroutine cedr_sl_set_Qdp(ie, Qdp, n0_qdp, np1_qdp) bind(c)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np, qsize_d
       integer, value, intent(in) :: ie, n0_qdp, np1_qdp
       real(kind=real_kind), intent(in) :: Qdp(np,np,nlev,qsize_d,2)
     end subroutine cedr_sl_set_Qdp

     subroutine cedr_sl_set_dp3d(ie, dp3d, tl_np1) bind(c)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np
       use element_state,  only : timelevels
       integer, value, intent(in) :: ie, tl_np1
       real(kind=real_kind), intent(in) :: dp3d(np,np,nlev,timelevels)
     end subroutine cedr_sl_set_dp3d

     subroutine cedr_sl_set_Q(ie, Q) bind(c)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np, qsize_d
       integer, value, intent(in) :: ie
       real(kind=real_kind), intent(in) :: Q(np,np,nlev,qsize_d)
     end subroutine cedr_sl_set_Q

     subroutine cedr_sl_set_pointers_end() bind(c)
     end subroutine cedr_sl_set_pointers_end

     subroutine cedr_sl_run(minq, maxq, nets, nete) bind(c)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np, qsize
       real(kind=real_kind), intent(in) :: minq(np,np,nlev,qsize,nets:nete)
       real(kind=real_kind), intent(in) :: maxq(np,np,nlev,qsize,nets:nete)
       integer, value, intent(in) :: nets, nete
     end subroutine cedr_sl_run

     subroutine cedr_sl_run_local(minq, maxq, nets, nete, use_ir, limiter_option) bind(c)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np, qsize
       real(kind=real_kind), intent(in) :: minq(np,np,nlev,qsize,nets:nete)
       real(kind=real_kind), intent(in) :: maxq(np,np,nlev,qsize,nets:nete)
       integer, value, intent(in) :: nets, nete, use_ir, limiter_option
     end subroutine cedr_sl_run_local

     subroutine cedr_sl_check(minq, maxq, nets, nete) bind(c)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np, qsize
       real(kind=real_kind), intent(in) :: minq(np,np,nlev,qsize,nets:nete)
       real(kind=real_kind), intent(in) :: maxq(np,np,nlev,qsize,nets:nete)
       integer, value, intent(in) :: nets, nete
     end subroutine cedr_sl_check

     subroutine slmm_init_local_mesh(ie, neigh_corners, num_neighbors, pinside) bind(c)
       use coordinate_systems_mod, only : cartesian3D_t
       integer, value, intent(in) :: ie, num_neighbors
       type(cartesian3D_t), intent(in) :: neigh_corners(:,:), pinside
     end subroutine slmm_init_local_mesh

     subroutine slmm_check_ref2sphere(ie, sphere_cart_coord) bind(c)
       use coordinate_systems_mod, only : cartesian3D_t
       integer, value, intent(in) :: ie
       type(cartesian3D_t), intent(in) :: sphere_cart_coord
     end subroutine slmm_check_ref2sphere

     subroutine slmm_advect(lev, ie, nnc, np, nlev, qsize, nets, nete, &
          dep_points, Qj_src, metdet, dp3d, tl_np1, q, minq, maxq) bind(c)
       use coordinate_systems_mod, only : cartesian3D_t
       use kinds, only : real_kind
       use dimensions_mod, only : qsize_d, max_neigh_edges
       use element_state, only : timelevels
       integer, value, intent(in) :: lev, ie, nnc, np, nlev, qsize, nets, nete, tl_np1
       type(cartesian3D_t), intent(in) :: dep_points(np,np)
       real(kind=real_kind), intent(in) :: Qj_src(np,np,qsize+1,max_neigh_edges+1), &
            metdet(np,np), dp3d(np,np,nlev,timelevels)
       real(kind=real_kind), intent(out) :: q(np,np,nlev,qsize_d), &
            minq(np,np,nlev,qsize,nets:nete), maxq(np,np,nlev,qsize,nets:nete)
     end subroutine slmm_advect

     subroutine slmm_csl_set_elem_data(ie, metdet, qdp, dp, q, nelem_in_patch) bind(c)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np, qsize
       real(kind=real_kind), intent(in) :: metdet(np,np), qdp(np,np,nlev,qsize), &
            dp(np,np,nlev), q(np,np,nlev,qsize)
       integer, value, intent(in) :: ie, nelem_in_patch
     end subroutine slmm_csl_set_elem_data

     subroutine slmm_csl(nets, nete, dep_points, minq, maxq, info) bind(c)
       use kinds         , only : real_kind
       use dimensions_mod, only : np, nlev, nelemd, qsize
       use coordinate_systems_mod, only : cartesian3D_t
       integer, value, intent(in) :: nets, nete
       ! dep_points is const in principle, but if lev <=
       ! semi_lagrange_nearest_point_lev, a departure point may be altered if
       ! the winds take it outside of the comm halo.
       type(cartesian3D_t), intent(inout) :: dep_points(np,np,nlev,nelemd)
       real(kind=real_kind), intent(in) :: &
            minq(np,np,nlev,qsize,nets:nete), maxq(np,np,nlev,qsize,nets:nete)
       integer, intent(out) :: info
     end subroutine slmm_csl

     ! Temporary for work on and characterization of the SL MPI pattern.
     subroutine slmm_get_mpi_pattern(sl_mpi) bind(c)
       integer, intent(out) :: sl_mpi
     end subroutine slmm_get_mpi_pattern
  end interface

contains

  subroutine compose_init(comm, elem, GridVertex)
    use dimensions_mod, only : np, nlev, qsize, qsize_d, nelem, nelemd
    use element_mod, only : element_t
    use gridgraph_mod, only : GridVertex_t
    use control_mod, only : semi_lagrange_cdr_alg, transport_alg, cubed_sphere_map, &
         semi_lagrange_nearest_point_lev

    integer, intent(in) :: comm
    type (element_t), intent(in) :: elem(:)
    type (GridVertex_t), intent(in), target :: GridVertex(:)
    integer, allocatable :: &
         sc2gci(:), sc2rank(:), &    ! space curve index -> (GID, rank)
         nbr_id_rank(:), nirptr(:)   ! (GID, rank) in local mesh patch, starting with own
    integer :: lid2gid(nelemd), lid2facenum(nelemd)
    integer :: i, j, k, sc, gid

    allocate(sc2gci(nelem), sc2rank(nelem))
    do i = 1, nelem
       sc = GridVertex(i)%SpaceCurve + 1
       sc2gci(sc) = i - 1
       sc2rank(sc) = GridVertex(i)%processor_number - 1
    end do
    call cedr_init_impl(comm, semi_lagrange_cdr_alg, sc2gci, sc2rank, &
         nelem, nelemd, nlev)
    deallocate(sc2gci, sc2rank)

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
          nbr_id_rank(k+1) = GridVertex(gid)%processor_number - 1
          k = k + 2
          do j = 1, size(elem(i)%desc%globalID_neigh_corners)
             gid = elem(i)%desc%globalID_neigh_corners(j)
             nbr_id_rank(k) = gid
             nbr_id_rank(k+1) = GridVertex(gid)%processor_number - 1
             k = k + 2
          end do
       end do
       nirptr(nelemd+1) = k - 1
       call slmm_init_impl(comm, transport_alg, np, nlev, qsize, qsize_d, &
            nelem, nelemd, cubed_sphere_map, lid2gid, lid2facenum, &
            nbr_id_rank, nirptr, semi_lagrange_nearest_point_lev)
       deallocate(nbr_id_rank, nirptr)
    end if
  end subroutine compose_init

  subroutine compose_repro_sum(send, recv, nlocal, nfld, comm) bind(c)
    use kinds, only: real_kind
    use repro_sum_mod, only: repro_sum

    real(kind=real_kind), intent(in) :: send(nlocal,nfld)
    real(kind=real_kind), intent(out) :: recv(nfld)
    integer, value, intent(in) :: nlocal, nfld, comm

    call repro_sum(send, recv, nlocal, nlocal, nfld, commid=comm)
  end subroutine compose_repro_sum

end module compose_mod
