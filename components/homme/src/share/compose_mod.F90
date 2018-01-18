module compose_mod

  interface
     ! For Kokkos.
     subroutine kokkos_init()
     end subroutine kokkos_init

     subroutine kokkos_finalize()
     end subroutine kokkos_finalize

     ! For QLT.
     subroutine cedr_unittest(comm, nerr)
       integer, intent(in) :: comm
       integer, intent(out) :: nerr
     end subroutine cedr_unittest

     subroutine cedr_init_impl(comm, cdr_alg, sc2gci, sc2rank, &
          ncell, nlclcell, nlev)
       integer, intent(in) :: comm, cdr_alg, sc2gci(:), sc2rank(:), &
            ncell, nlclcell, nlev
     end subroutine cedr_init_impl

     subroutine cedr_set_ie2gci(ie, gci)
       integer, intent(in) :: ie, gci
     end subroutine cedr_set_ie2gci

     ! For sl_advection's use of QLT.
     subroutine cedr_sl_init(np, nlev, qsize, qsize_d, timelevels, &
          need_conservation)
       integer, intent(in) :: np, nlev, qsize, qsize_d, timelevels, &
            need_conservation
     end subroutine cedr_sl_init

     subroutine cedr_sl_set_pointers_begin(nets, nete)
       integer, intent(in) :: nets, nete
     end subroutine cedr_sl_set_pointers_begin

     subroutine cedr_sl_set_spheremp(ie, spheremp)
       use kinds         , only : real_kind
       use dimensions_mod, only : np
       integer, intent(in) :: ie
       real(kind=real_kind), intent(in) :: spheremp(np,np)
     end subroutine cedr_sl_set_spheremp

     subroutine cedr_sl_set_Qdp(ie, Qdp, n0_qdp, np1_qdp)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np, qsize_d
       integer, intent(in) :: ie, n0_qdp, np1_qdp
       real(kind=real_kind), intent(in) :: Qdp(np,np,nlev,qsize_d,2)
     end subroutine cedr_sl_set_Qdp

     subroutine cedr_sl_set_dp3d(ie, dp3d, tl_np1)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np
       use element_state,  only : timelevels
       integer, intent(in) :: ie, tl_np1
       real(kind=real_kind), intent(in) :: dp3d(np,np,nlev,timelevels)
     end subroutine cedr_sl_set_dp3d

     subroutine cedr_sl_set_Q(ie, Q)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np, qsize_d
       integer, intent(in) :: ie
       real(kind=real_kind), intent(in) :: Q(np,np,nlev,qsize_d)
     end subroutine cedr_sl_set_Q

     subroutine cedr_sl_set_pointers_end()
     end subroutine cedr_sl_set_pointers_end

     subroutine cedr_sl_run(minq, maxq, nets, nete)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np, qsize
       real(kind=real_kind), intent(in) :: minq(np*np,nlev,qsize,nets:nete)
       real(kind=real_kind), intent(in) :: maxq(np*np,nlev,qsize,nets:nete)
       integer             , intent(in) :: nets, nete
     end subroutine cedr_sl_run

     subroutine cedr_sl_run_local(minq, maxq, nets, nete, use_ir, limiter_option)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np, qsize
       real(kind=real_kind), intent(in) :: minq(np*np,nlev,qsize,nets:nete)
       real(kind=real_kind), intent(in) :: maxq(np*np,nlev,qsize,nets:nete)
       integer             , intent(in) :: nets, nete, use_ir, limiter_option
     end subroutine cedr_sl_run_local

     subroutine cedr_sl_check(minq, maxq, nets, nete)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np, qsize
       real(kind=real_kind), intent(in) :: minq(np*np,nlev,qsize,nets:nete)
       real(kind=real_kind), intent(in) :: maxq(np*np,nlev,qsize,nets:nete)
       integer             , intent(in) :: nets, nete
     end subroutine cedr_sl_check

     ! For sl_advection's use of incremental remap (IR).
     subroutine slmm_init(np, nelem)
       integer, intent(in) :: np, nelem
     end subroutine slmm_init

     subroutine slmm_init_local_mesh(ie, neigh_corners, num_neighbors)
       use coordinate_systems_mod, only : cartesian3D_t
       integer, intent(in) :: ie, num_neighbors
       type(cartesian3D_t), intent(in) :: neigh_corners(:,:)
     end subroutine slmm_init_local_mesh

     subroutine slmm_study(globalID, corners, globalIDs, neigh_corners, n)
       use coordinate_systems_mod, only : cartesian3D_t
       integer, intent(in) :: globalID, globalIDs(:), n
       type(cartesian3D_t), intent(in) :: corners(4), neigh_corners(:,:)
     end subroutine slmm_study

     subroutine slmm_project(lev, ie, nnc, np, nlev, qsize, nets, nete, &
          dep_points, Qj_src, metdet, dp3d, tl_np1, q, minq, maxq)
       use coordinate_systems_mod, only : cartesian3D_t
       use kinds, only : real_kind
       use dimensions_mod, only : qsize_d, max_neigh_edges
       use element_state, only : timelevels
       integer, intent(in) :: lev, ie, nnc, np, nlev, qsize, nets, nete, tl_np1
       type(cartesian3D_t), intent(in) :: dep_points(np,np)
       real(kind=real_kind), intent(in) :: Qj_src(np,np,qsize+1,max_neigh_edges+1), &
            metdet(np,np), dp3d(np,np,nlev,timelevels)
       real(kind=real_kind), intent(out) :: q(np,np,nlev,qsize_d), &
            minq(np,np,nlev,qsize,nets:nete), maxq(np,np,nlev,qsize,nets:nete)
     end subroutine slmm_project
  end interface

contains

  subroutine cedr_init(comm, nlclelem, GridVertex)
    use dimensions_mod, only : nlev
    use gridgraph_mod, only : GridVertex_t
    use control_mod, only : semi_lagrange_cdr_alg
    integer, intent(in) :: comm, nlclelem
    type (GridVertex_t), intent(in), target :: GridVertex(:)
    integer, allocatable :: sc2gci(:), sc2rank(:)
    integer :: nelem, i, sc

    nelem = size(GridVertex)
    allocate(sc2gci(nelem), sc2rank(nelem))
    do i = 1, nelem
       sc = GridVertex(i)%SpaceCurve + 1
       sc2gci(sc) = i - 1
       sc2rank(sc) = GridVertex(i)%processor_number - 1
    end do
    call cedr_init_impl(comm, semi_lagrange_cdr_alg, sc2gci, sc2rank, &
         nelem, nlclelem, nlev)
    deallocate(sc2gci, sc2rank)
  end subroutine cedr_init

end module compose_mod
