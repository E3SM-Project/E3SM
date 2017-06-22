module compose_mod

  interface
     ! For Kokkos.
     subroutine kokkos_init()
     end subroutine kokkos_init

     subroutine kokkos_finalize()
     end subroutine kokkos_finalize

     ! For QLT.
     subroutine qlt_unittest(comm, nerr)
       integer, intent(in) :: comm
       integer, intent(out) :: nerr
     end subroutine qlt_unittest

     subroutine qlt_init_impl(comm, sc2gci, sc2rank, ncell)
       integer, intent(in) :: comm, sc2gci(:), sc2rank(:), ncell
     end subroutine qlt_init_impl

     subroutine qlt_set_ie2gci(ie, gci)
       integer, intent(in) :: ie, gci
     end subroutine qlt_set_ie2gci

     ! For sl_advection's use of QLT.
     subroutine qlt_sl_set_pointers_begin(nets, nete, np, nlev, qsize, qsize_d, &
          timelevels, need_conservation)
       integer, intent(in) :: nets, nete, np, nlev, qsize, qsize_d, timelevels, &
            need_conservation
     end subroutine qlt_sl_set_pointers_begin

     subroutine qlt_sl_set_spheremp(ie, spheremp)
       use kinds         , only : real_kind
       use dimensions_mod, only : np
       integer, intent(in) :: ie
       real(kind=real_kind), intent(in) :: spheremp(np,np)
     end subroutine qlt_sl_set_spheremp

     subroutine qlt_sl_set_Qdp(ie, Qdp, n0_qdp, np1_qdp)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np, qsize_d
       integer, intent(in) :: ie, n0_qdp, np1_qdp
       real(kind=real_kind), intent(in) :: Qdp(np,np,nlev,qsize_d,2)
     end subroutine qlt_sl_set_Qdp

     subroutine qlt_sl_set_dp3d(ie, dp3d, tl_np1)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np
       use element_state,  only : timelevels
       integer, intent(in) :: ie, tl_np1
       real(kind=real_kind), intent(in) :: dp3d(np,np,nlev,timelevels)
     end subroutine qlt_sl_set_dp3d

     subroutine qlt_sl_set_Q(ie, Q)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np, qsize_d
       integer, intent(in) :: ie
       real(kind=real_kind), intent(in) :: Q(np,np,nlev,qsize_d)
     end subroutine qlt_sl_set_Q

     subroutine qlt_sl_set_pointers_end()
     end subroutine qlt_sl_set_pointers_end

     subroutine qlt_sl_run(minq, maxq, nets, nete)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np, qsize
       real(kind=real_kind), intent(in) :: minq(np*np,nlev,qsize,nets:nete)
       real(kind=real_kind), intent(in) :: maxq(np*np,nlev,qsize,nets:nete)
       integer             , intent(in) :: nets, nete
     end subroutine qlt_sl_run

     subroutine qlt_sl_run_local(minq, maxq, nets, nete, use_ir)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np, qsize
       real(kind=real_kind), intent(in) :: minq(np*np,nlev,qsize,nets:nete)
       real(kind=real_kind), intent(in) :: maxq(np*np,nlev,qsize,nets:nete)
       integer             , intent(in) :: nets, nete, use_ir
     end subroutine qlt_sl_run_local

     subroutine qlt_sl_check(minq, maxq, nets, nete)
       use kinds         , only : real_kind
       use dimensions_mod, only : nlev, np, qsize
       real(kind=real_kind), intent(in) :: minq(np*np,nlev,qsize,nets:nete)
       real(kind=real_kind), intent(in) :: maxq(np*np,nlev,qsize,nets:nete)
       integer             , intent(in) :: nets, nete
     end subroutine qlt_sl_check

     ! For sl_advection's use of incremental remap (IR).
     subroutine ir_init(np, nelem)
       integer, intent(in) :: np, nelem
     end subroutine ir_init

     subroutine ir_study(globalID, corners, globalIDs, neigh_corners, n)
       use coordinate_systems_mod, only : cartesian3D_t
       integer, intent(in) :: globalID, globalIDs(:), n
       type(cartesian3D_t), intent(in) :: corners(4), neigh_corners(:,:)
     end subroutine ir_study

     subroutine ir_project(lev, ie, nnc, np, nlev, qsize, nets, nete, &
          dep_points, neigh_corners, &
          Qj_src, &
          metdet, dp3d, tl_np1, &
          q, minq, maxq)
       use coordinate_systems_mod, only : cartesian3D_t
       use kinds, only : real_kind
       use dimensions_mod, only : qsize_d, max_neigh_edges
       use element_state, only : timelevels
       integer, intent(in) :: lev, ie, nnc, np, nlev, qsize, nets, nete, tl_np1
       type(cartesian3D_t), intent(in) :: dep_points(np,np), neigh_corners(:,:)
       real(kind=real_kind), intent(in) :: Qj_src(np,np,qsize+1,max_neigh_edges+1), &
            metdet(np,np), dp3d(np,np,nlev,timelevels)
       real(kind=real_kind), intent(out) :: q(np,np,nlev,qsize_d), &
            minq(np,np,nlev,qsize,nets:nete), maxq(np,np,nlev,qsize,nets:nete)
     end subroutine ir_project
  end interface

contains

  subroutine qlt_init(comm, GridVertex)
    use gridgraph_mod, only : GridVertex_t
    integer, intent(in) :: comm
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
    call qlt_init_impl(comm, sc2gci, sc2rank, nelem)
    deallocate(sc2gci, sc2rank)
  end subroutine qlt_init

end module compose_mod
